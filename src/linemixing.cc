#include <numeric>

#include <Faddeeva/Faddeeva.hh>

#include "lin_alg.h"
#include "linefunctions.h"
#include "linemixing.h"


namespace Absorption::LineMixing {
EquivalentLines::EquivalentLines(const ComplexMatrix& W,
                                 const Vector& pop,
                                 const Vector& dip) noexcept :
  val(pop.nelem(), 0), str(pop.nelem(), 0)
{
  const Index n = pop.nelem();
  
  // Compute values
  ComplexMatrix V(n, n);
  
  // Main computations
  diagonalize(V, val, W);
  
  // Do the matrix forward multiplication
  for (Index i=0; i<n; i++) {
    for (Index j=0; j<n; j++) {
      str[i] += dip[j] * V(j, i);
    }
  }
  
  // Do the matrix backward multiplication
  auto& invV=V.inv();
  for (Index i=0; i<n; i++) {
    Complex z(0, 0);
    for (Index j=0; j<n; j++) {
      z += pop[j] * dip[j] * invV(i, j);
    }
    str[i] *= z;
  }
}

ArrayOfIndex PopulationAndDipole::sort(const AbsorptionLines& band) noexcept {
  const Index N = pop.nelem();
  
  // List that starts as [0, 1, ... N-2, N-1]
  ArrayOfIndex out(N);
  std::iota(out.begin(), out.end(), 0);
  
  // Strength
  Vector s(N);
  for (Index i=0; i<N; i++) {
    s[i] = band.F0(i) * pop[i] * Constant::pow2(dip[i]);
  }
  
  for (Index i=0; i<N; i++) {
    for (Index j=i+1; j<N; j++) {
      if (s[j] > s[i]) {
        std::swap(out[i], out[j]);
        std::swap(dip[i], dip[j]);
        std::swap(pop[i], pop[j]);
        std::swap(s[i], s[j]);
      }
    }
  }
  
  return out;
}

ComplexMatrix relaxation_matrix_makarov2020(const Numeric T,
                                            const Numeric P,
                                            const Vector& vmrs,
                                            const Vector& pop,  // Already sorted
                                            const AbsorptionLines& band,
                                            const ArrayOfIndex& sorting) {
  const Index N = band.NumLines();
  
  // Create output
  ComplexMatrix W(N, N);
  
  // Diagonal and upper triangular values
  for (Index I=0; I<N; I++) {
    const Index i = sorting[I];
    auto Ji = band.UpperQuantumNumber(i, QuantumNumberType::J);
    auto Jf = band.LowerQuantumNumber(i, QuantumNumberType::J);
    auto Ni = band.UpperQuantumNumber(i, QuantumNumberType::N);
    auto Nf = band.LowerQuantumNumber(i, QuantumNumberType::N);
    for (Index J=0; J<N; J++) {
      const Index j = sorting[J];
      auto Ji_p = band.UpperQuantumNumber(j, QuantumNumberType::J);
      auto Jf_p = band.LowerQuantumNumber(j, QuantumNumberType::J);
      auto Ni_p = band.UpperQuantumNumber(j, QuantumNumberType::N);
      auto Nf_p = band.LowerQuantumNumber(j, QuantumNumberType::N);
      if (I not_eq J) {
        W(I, J) = Complex(0, o2_ecs_wigner_symbol_tran(Ji, Jf, Ni, Nf, 1_rat, 1_rat, Ji_p, Jf_p, Ni_p, Nf_p, 1_rat, T));
        W(J, I) = W(I, J) * pop[J] / pop[I];
      } else {
        W(I, J) = Complex(0, band.ShapeParameters(i, T, P, vmrs).G0);
      }
    }
  }
  
  // Weird undocumented minus sign
  for (Index I=0; I<N; I++) {
    for (Index J=0; J<N; J++) {
      if (J not_eq I) {
        W(I, J) = - std::abs(W(I, J));
      }
    }
  }
  
  // Reduced dipole
  Vector dip0(N);
  for (Index I=0; I<N; I++) {
    const Index i = sorting[I];
    dip0[I] = o2_makarov2013_reduced_dipole(band.UpperQuantumNumber(i, QuantumNumberType::J),
                                            band.LowerQuantumNumber(i, QuantumNumberType::J),
                                            band.UpperQuantumNumber(i, QuantumNumberType::N));
  }
  
  // Sum rule correction
  for (Index I=0; I<N; I++) {
    Numeric sumlw = 0.0;
    Numeric sumup = 0.0;
    
    for (Index J=0; J<N; J++) {
      if (J > I) {
        sumlw += std::abs(dip0[J]) * W(J, I).imag();
      } else {
        sumup += std::abs(dip0[J]) * W(J, I).imag();
      }
    }
    
    for (Index J=I+1; J<N; J++) {
      if (sumlw == 0) {
        W(J, I) = 0.0;
        W(I, J) = 0.0;
      } else {
        W(J, I) *= - sumup / sumlw;
        W(I, J) = W(J, I) * pop[I] / pop[J];
      }
    }
  }
  
  // Add the frequency on the real diagonal
  for (Index I=0; I<N; I++) {
    W(I, I) += band.F0(sorting[I]);
  }
  
  return W;
}


ComplexVector linemixing_ecs_makarov2020(const Numeric T,
                                         const Numeric P,
                                         const Vector& vmrs,
                                         const Vector& f_grid,
                                         const AbsorptionLines& band,
                                         const SpeciesAuxData::AuxType& partition_type,
                                         const ArrayOfGriddedField1& partition_data)
{
  // Constants for the band at this temperature
  const Numeric frenorm = band.F_mean();
  const Numeric GD_div_F0 = Linefunctions::DopplerConstant(T, band.SpeciesMass());
  
  // Sorted population
  PopulationAndDipole tp = PopulationAndDipole(T, band, partition_type, partition_data);
  const ArrayOfIndex sorting = tp.sort(band);
  
  // Relaxation matrix
  ComplexMatrix W = relaxation_matrix_makarov2020(T, P, vmrs, tp.pop, band, sorting);
  W.diagonal() -= frenorm;
  
  // Equivalent lines computations
  EquivalentLines eqv(W, tp.pop, tp.dip);
  eqv.val += frenorm;
  
  // Absorption of this band
  ComplexVector a(f_grid.nelem(), 0);
  for (Index i=0; i<band.NumLines(); i++) {
    const Numeric gamd = GD_div_F0 * eqv.val[i].real();
    const Numeric cte = Constant::sqrt_ln_2 / gamd;
    for (Index iv=0; iv<f_grid.nelem(); iv++) {
      const Complex z = (eqv.val[i]-f_grid[iv]) * cte;
      const Complex w = Faddeeva::w(z);
      a[iv] += (eqv.str[i] * w) / gamd;
    }
  }
  return a;
}
}
