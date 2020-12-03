#include <numeric>

#include <Faddeeva/Faddeeva.hh>

#include "lin_alg.h"
#include "linefunctions.h"
#include "linemixing.h"
#include "physics_funcs.h"


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


void relaxation_matrix_makarov2020_offdiagonal(MatrixView W,
                                               const Numeric T,
                                               const AbsorptionLines& band,
                                               const ArrayOfIndex& sorting,
                                               const Numeric& mass1,
                                               const Numeric& mass2)
{
  const Index n = band.NumLines();
  
  // Energy calculations
  const auto Erot = [](auto J) -> Numeric {
    using Constant::h;
    using Constant::pow2;
    using Constant::pow3;
    
    // Constant for O2-66 energy
    constexpr auto B = 43100.44276e6;
    constexpr auto D = 145.1271e3;
    constexpr auto H = 49e-3;
    constexpr auto lB = 59501.3438e6;
    constexpr auto lD = 58.3680e3;
    constexpr auto lH = 290.8e-3;
    constexpr auto gB = -252.58634e6;
    constexpr auto gD = -243.42;
    constexpr auto gH = -1.46e-3;
    
    const auto lJ = Numeric(J * (J + 1));
    return h * (B*lJ - D*pow2(lJ) + H*pow3(lJ) -
               (gB + gD*lJ + gH*pow2(lJ)) +
     2.0/3.0 * (lB + lD*lJ + lH*pow2(lJ)));
  };
  
  // Adiabatic factor
  const auto Omega = [Erot](auto N, auto temp, auto m1, auto m2) -> Numeric {
    using Constant::h;
    using Constant::k;
    using Constant::pi;
    using Constant::h_bar;
    using Constant::pow2;
    
    // Constant from Makarov etal 2020
    constexpr Numeric dc = Conversion::angstrom2meter(0.61);
                         
    const Numeric en = Erot(N);
    const Numeric enm2 = Erot(N-2);
    const Numeric wnnm2 = (en - enm2) / h_bar;
    
    const Numeric mu = 1 / m1 + 1 / m2;
    const Numeric v_bar_pow2 = 8*k*temp*mu/pi;
    const Numeric tauc_pow2 = pow2(dc) / v_bar_pow2;
    
    return 1.0 / pow2(1 + 1.0/24.0 * pow2(wnnm2) * tauc_pow2);
  };
  
  // Basis rate
  const auto Q = [Erot](auto L, auto temp) -> Numeric {
    using Constant::k;
    
    // Constants from Makarov etal 2020
    constexpr Numeric lambda =  0.39;
    constexpr Numeric beta = 0.567;
    
    auto el = Erot(L);
    
    return Numeric(2*L + 1) / pow(L * (L+1), lambda) * std::exp(-beta * el / (k*temp));
  };
  
  const auto zero_dipole = [](auto Jup, auto Jlo, auto N) {
    return (even(Jlo + N) ? 1 : -1) * sqrt(6 * (2*Jlo + 1) * (2*Jup + 1)) * wigner6j(1, 1, 1, Jup, Jlo, N);
  };
  
  Vector dip0(n);
  for (Index i=0; i<n; i++) {
    dip0[i] = std::abs(zero_dipole(band.UpperQuantumNumber(sorting[i], QuantumNumberType::J),
                                   band.LowerQuantumNumber(sorting[i], QuantumNumberType::J),
                                   band.UpperQuantumNumber(sorting[i], QuantumNumberType::N)));
    
    for (Index j=0; j<n; j++) {
      if (i == j) continue;
      
      // Select upper quantum number
      const bool ihigh = band.E0(sorting[i]) > band.E0(sorting[j]);
      const Index k = ihigh ? i : j;
      const Index l = ihigh ? j : i;
      
      // Quantum numbers
      const Rational Jk = band.UpperQuantumNumber(sorting[k], QuantumNumberType::J);
      const Rational Jl = band.UpperQuantumNumber(sorting[l], QuantumNumberType::J);
      const Rational Nk = band.UpperQuantumNumber(sorting[k], QuantumNumberType::N);
      const Rational Nl = band.UpperQuantumNumber(sorting[l], QuantumNumberType::N);
      const Rational Jk_p = band.LowerQuantumNumber(sorting[k], QuantumNumberType::J);
      const Rational Jl_p = band.LowerQuantumNumber(sorting[l], QuantumNumberType::J);
      
      // Makarov 2013 symbol with modifications:
      //    1) Squared scale
      //    2) Squared 3J-symbol
      //    3) Using the updated 2020 constants
      // These are modified after reading the original code
      Numeric sum = 0;
      const Numeric scale = Numeric((2*Nk + 1) * (2*Nl + 1)) * sqrt((2*Jk + 1) * (2*Jl + 1) * (2*Jk_p + 1) * (2* Jl_p + 1));
      const auto [L0, L1] = wigner_limits(wigner3j_limits<3>(Nl, Nk), {Rational(2), 100000});
      for (Rational L=L0; L<L1; L+=2) {
        const Numeric sgn = even(Jk + Jl + L + 1) ? 1 : -1;
        const Numeric a = Constant::pow2(wigner3j(Nl, Nk, L, 0, 0, 0));
        const Numeric b = wigner6j(L, Jk, Jl, 1, Nl, Nk);
        const Numeric c = wigner6j(L, Jk_p, Jl_p, 1, Nl, Nk);
        const Numeric d = wigner6j(L, Jk, Jl, 1, Jl_p, Jk_p);
        sum += sgn * a * b * c * d * Q(L, T) / Omega(L, T, mass1, mass2);
      }
      sum *= scale * Omega(Nk, T, mass1, mass2);
      
      // Add to W and rescale to upwards element by the populations
      W(l, k) = sum;
      W(k, l) = sum * std::exp((Erot(Jl) - Erot(Jk)) / (Constant::k*T));
    }
  }
  
  // Transpose?  Why is this transpose required?
  if constexpr (1) {
    for (Index i=0; i<n; i++) {
      for (Index j=0; j<i; j++) {
        swap(W(i, j), W(j, i));
      }
    }
  }
  
  // Sum rule correction
  for (Index i=0; i<n; i++) {
    Numeric sumlw = 0.0;
    Numeric sumup = 0.0;
    
    for (Index j=0; j<n; j++) {
      if (j > i) {
        sumlw += std::abs(dip0[j]) * W(j, i);
      } else {
        sumup += std::abs(dip0[j]) * W(j, i);
      }
    }
    
    const Rational Ji = band.UpperQuantumNumber(sorting[i], QuantumNumberType::J);
    for (Index j=i+1; j<n; j++) {
      const Rational Jj = band.UpperQuantumNumber(sorting[j], QuantumNumberType::J);
      if (sumlw == 0) {
        W(j, i) = 0.0;
        W(i, j) = 0.0;
      } else {
        W(j, i) *= - sumup / sumlw;
        W(i, j) = W(j, i) * std::exp((Erot(Ji) - Erot(Jj)) / (Constant::k*T));
      }
    }
  }
}


ComplexMatrix relaxation_matrix(const Numeric T,
                                const Numeric P,
                                const Vector& vmrs,
                                const Vector& mass,
                                const AbsorptionLines& band,
                                const ArrayOfIndex& sorting,
                                const Numeric frenorm) {
  const Index N = band.NumLines();
  const Index M = vmrs.nelem();
  
  // Create output
  ComplexMatrix W(N, N, 0);
  
  // Loop over all the broadeners
  for (Index J; J<M; J++) {
    // Create temporary
    ComplexMatrix Wtmp(N, N, 0);
    
    // Fill diagonal keeping track of the reshuffle (sorting)
    for (Index I=0; I<N; I++) {
      const Index i = sorting[I];
      auto shape = band.ShapeParameters(i, T, P, vmrs);
      Wtmp(I, I) = Complex(shape.D0, shape.G0);
    }

    // Set the off-diagonal part of the matrix for this broadener
    relaxation_matrix_makarov2020_offdiagonal(Wtmp.imag(), T, band, sorting,
                                              Constant::m_u * band.SpeciesMass(),
                                              Constant::m_u * mass[J]);
    
    // Sum up all atmospheric components
    for (Index i=0; i<N; i++) {
      for (Index j=0; j<N; j++) {
        W(i, j) += vmrs[J] * Wtmp(i, j);
      }
    }
  }
  
  // Deal with line frequency and its re-normalization
  for (Index I=0; I<N; I++) {
    W(I, I) += band.F0(sorting[I]) - frenorm;
  }
  
  return W;
}


std::pair<ArrayOfIndex, PopulationAndDipole> sorted_population_and_dipole(const Numeric T,
                                                                          const AbsorptionLines& band,
                                                                          const SpeciesAuxData::AuxType& partition_type,
                                                                          const ArrayOfGriddedField1& partition_data) {
  PopulationAndDipole tp(T, band, partition_type, partition_data);
  return {tp.sort(band), tp};
}


ComplexVector linemixing_ecs_absorption(const Numeric T,
                                        const Numeric P,
                                        const Numeric this_vmr,
                                        const Vector& vmrs,
                                        const Vector& mass,
                                        const Vector& f_grid,
                                        const AbsorptionLines& band,
                                        const SpeciesAuxData::AuxType& partition_type,
                                        const ArrayOfGriddedField1& partition_data)
{
  // Weighted center of the band
  const Numeric frenorm = band.F_mean();
  
  // Band Doppler broadening constant
  const Numeric GD_div_F0 = Linefunctions::DopplerConstant(T, band.SpeciesMass());
  
  // Sorted population
  const auto [sorting, tp] = sorted_population_and_dipole(T, band, partition_type, partition_data);
  
  // Relaxation matrix
  const ComplexMatrix W = relaxation_matrix(T, P, vmrs, mass, band, sorting, frenorm);
  
  // Equivalent lines computations
  const EquivalentLines eqv(W, tp.pop, tp.dip);
  
  // Absorption of this band
  ComplexVector absorption(f_grid.nelem(), 0);
  for (Index i=0; i<band.NumLines(); i++) {
    const Numeric gamd = GD_div_F0 * (eqv.val[i].real() + frenorm);
    const Numeric cte = Constant::sqrt_ln_2 / gamd;
    for (Index iv=0; iv<f_grid.nelem(); iv++) {
      const Complex z = (eqv.val[i] + frenorm - f_grid[iv]) * cte;
      const Complex w = Faddeeva::w(z);
      absorption[iv] += eqv.str[i] * w / gamd;
    }
  }
  
  // Adjust by frequency and number density
  const Numeric numdens = this_vmr * number_density(P, T);
  constexpr Numeric sq_ln2pi = Constant::sqrt_ln_2 / Constant::sqrt_pi;
  for (Index iv=0; iv<f_grid.nelem(); iv++) {
    const Numeric f = f_grid[iv];
    const Numeric fact = f * (1 - stimulated_emission(T, f));
    absorption[iv] *= fact * numdens * sq_ln2pi;
  }
  
  return absorption;
}
}
