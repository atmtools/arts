#include <numeric>

#include <Faddeeva/Faddeeva.hh>

#include "lin_alg.h"
#include "linefunctions.h"
#include "linemixing.h"
#include "minimize.h"
#include "physics_funcs.h"

namespace Absorption::LineMixing {
EquivalentLines::EquivalentLines(const ComplexMatrix& W,
                                 const Vector& pop,
                                 const Vector& dip) noexcept :
  val(pop.nelem(), 0), str(pop.nelem(), 0) {
  /* FIXME:  (Added 2021-01-19; Richard Larsson)
    * 
    * This function cannot easily be used with partial derivatives
    * 
    * Doing so would allow the entire ECS approach to be analytical
    * in its derivatives.
    * 
    * The problem in short:
    *    There is an Eigenvalue decomposition happening (W = V diag(e) V^-1)
    * 
    *    We use V and e in a strange way.  We do not know how to get the
    *    partial derivatives of these two variables
    * 
    * The question:
    *    How do we compute dV and de?  We can get dW somewhat easily.
    */
    
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


void EquivalentLines::sort_by_frequency(Vector& f, const ArrayOfIndex& sorting) {
  const Index N = val.nelem();
  
  // First sort by frequency and assume this sorting
  // is a shifted version of the line centers
  for (Index i=0; i<N; i++) {
    for (Index j=i+1; j<N; j++) {
      if (val.real()[j] < val.real()[i]) {
        std::swap(val[i], val[j]);
        std::swap(str[i], str[j]);
      }
    }
  }
  
  const Vector fc = f;
  const ComplexVector strc = str;
  const ComplexVector valc = val;
  for (Index i=0; i<N; i++) {
    f[i] = fc[sorting[i]];
    str[i] = strc[sorting[i]];
    val[i] = valc[sorting[i]];
  }
}



namespace Makarov2020etal {
/*! Returns the reduced dipole
 * 
 * @param[in] Ju Main rotational number with spin of the upper level
 * @param[in] Jl Main rotational number with spin of the lower level
 * @param[in] N Main rotational number of both levels
 * @return The reduced dipole
 */
Numeric zero_dipole(const Rational Ju, const Rational Jl, const Rational N) {
  return (iseven(Jl + N) ? 1 : -1) * sqrt(6 * (2*Jl + 1) * (2*Ju + 1)) * wigner6j(1, 1, 1, Jl, Ju, N);
};
}


PopulationAndDipole::PopulationAndDipole(const Numeric T,
                                         const AbsorptionLines& band) noexcept :
  pop(band.NumLines()), dip(band.NumLines()) {
  const Index N = band.NumLines();
  
  const Numeric QT = single_partition_function(T, band.Isotopologue());
  const Numeric QT0 = single_partition_function(band.T0(), band.Isotopologue());
  const Numeric ratiopart = QT0 / QT;
  
  for (Index i=0; i<N; i++) {
    const Numeric pop0 = (band.g_upp(i) / QT0) * boltzman_factor(band.T0(), band.E0(i));
    pop[i] = pop0 * ratiopart * boltzman_ratio(T, band.T0(), band.E0(i));
    dip[i] = std::sqrt(- band.I0(i)/(pop0 * band.F0(i) * std::expm1(- (Constant::h * band.F0(i)) / (Constant::k * band.T0()))));
    
    // Adjust the sign depending on type
    if (band.Population() == Absorption::PopulationType::ByMakarovFullRelmat) {
        dip[i] *= std::signbit(Makarov2020etal::zero_dipole(band.UpperQuantumNumber(i, QuantumNumberType::J),
                                                            band.LowerQuantumNumber(i, QuantumNumberType::J),
                                                            band.UpperQuantumNumber(i, QuantumNumberType::N))) ? - 1 : 1;
    }
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
  
  // Sorter
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


void PopulationAndDipole::sort(const ArrayOfIndex& presorting) noexcept {
  const Index N = presorting.size();
  Vector dipcopy(dip), popcopy(pop);
  for (Index i=0; i<N; i++) {
    dip[i] = dipcopy[presorting[i]];
    pop[i] = popcopy[presorting[i]];
  }
}


//! Special parameter modifications
ENUMCLASS(SpecialParam, char,
          None,
          ModifyParameter1,
          ModifyParameter2,
          ModifyParameter3
)


namespace Makarov2020etal {
/*! Compute the rotational energy of ground-state O2 at N and J
 * 
 * If the template argument evaluates true, the erot<false>(1, 0)
 * energy is removed from the output of erot<false>(N, J).
 * 
 * @param[in] N Main rotational number
 * @param[in] j Main rotational number plus spin (if j < 0 then J=N)
 * @return Rotational energy in Joule
 */
template<bool rescale_pure_rotational=true> constexpr
Numeric erot(const Rational N, const Rational j=-1) {
  const Rational J = j<0 ? N : j;
  
  if constexpr (rescale_pure_rotational) {
    return erot<false>(N, J) - erot<false>(1, 0);
  } else {
    using Conversion::mhz2joule;
    using Constant::pow2;
    using Constant::pow3;
    
    constexpr Numeric B0=43100.4425e0;
    constexpr Numeric D0=.145123e0;
    constexpr Numeric H0=3.8e-08;
    constexpr Numeric xl0=59501.3435e0;
    constexpr Numeric xg0=-252.58633e0;
    constexpr Numeric xl1=0.058369e0;
    constexpr Numeric xl2=2.899e-07;
    constexpr Numeric xg1=-2.4344e-04;
    constexpr Numeric xg2=-1.45e-09;
    
    const Numeric XN = Numeric(N);
    const Numeric XX = XN * (XN + 1);
    const Numeric xlambda=xl0+xl1*XX+xl2*pow2(XX);
    const Numeric xgama=xg0+xg1*XX+xg2*pow2(XX);
    const Numeric C1=B0*XX-D0*pow2(XX)+H0*pow3(XX);
    
    if (J < N) {
      if (N == 1) {  // erot<false>(1, 0)
        return mhz2joule(C1 - (xlambda+B0*(2.*XN-1.)+xgama*XN));
      } else {
        return mhz2joule(C1 - (xlambda+B0*(2.*XN-1.)+xgama*XN) + std::sqrt(pow2(B0*(2.*XN-1.))+pow2(xlambda)-2.*B0*xlambda));
      }
    } else if (J > N) {
      return mhz2joule(C1 - (xlambda-B0*(2.*XN+3.)-xgama*(XN+1.)) - std::sqrt(pow2(B0*(2.*XN+3.))+pow2(xlambda)-2.*B0*xlambda));
    } else {
      return mhz2joule(C1);
    }
  }
}


/*! Returns the adiabatic factor
 * 
 * @param[in] N Main rotational number
 * @param[in] T The temperature
 * @param[in] m1 Mass of one species in atomic mass units
 * @param[in] m2 Mass of the other species in atomic mass units
 * @return The adiabatic factor
 */
template <SpecialParam param>
constexpr Numeric Omega(const Rational N,
                        const Numeric T,
                        const Numeric m1,
                        const Numeric m2) {
  using Constant::h;
  using Constant::k;
  using Constant::pi;
  using Constant::h_bar;
  using Constant::pow2;
  
  // Constant from Makarov etal 2020 (modified by 1e-6 in special Jacobian case)
  constexpr Numeric dc = param == SpecialParam::ModifyParameter1 ?
    Conversion::angstrom2meter(0.61 + 1e-6) :
    Conversion::angstrom2meter(0.61);
  
  // Constants for the expression
  constexpr Numeric fac = 8 * k / (Constant::m_u * pi);
  
  // nb. Only N=J considered???
  const Numeric en = erot(N);
  const Numeric enm2 = erot(N-2);
  const Numeric wnnm2 = (en - enm2) / h_bar;
  
  const Numeric mu = 1 / m1 + 1 / m2;
  const Numeric v_bar_pow2 = fac*T*mu;
  const Numeric tauc_pow2 = pow2(dc) / v_bar_pow2;
  
  return 1.0 / pow2(1 + 1.0/24.0 * pow2(wnnm2) * tauc_pow2);
};


/*! Returns the basis rate
 * 
 * @param[in] N Main rotational number
 * @param[in] T The temperature
 * @return The basis rate
 */
template <SpecialParam param>
Numeric Q(const Rational N, const Numeric T) {
  using Conversion::kelvin2joule;
  
  // Constants from Makarov etal 2020 (modified by 1e-6 in special Jacobian case)
  constexpr Numeric lambda = param == SpecialParam::ModifyParameter2 ? 0.39 + 1e-6 : 0.39;
  constexpr Numeric beta = param == SpecialParam::ModifyParameter3 ? 0.567 + 1e-6 : 0.567;
  
  // nb. Only N=J considered???
  return Numeric(2*N + 1) / pow(N * (N+1), lambda) * std::exp(-beta * erot(N) / kelvin2joule(T));
};

/*! Computes the off-diagonal elements of the relaxation matrix
 * following Makarov etal 2020.
 * 
 * @param[in,out] W The imaginary part of the relaxation matrix
 * @param[in] T The temperature
 * @param[in] band The absorption band
 * @param[in] sorting The sorting of the band
 * @param[in] mass_self The mass of the self-species
 * @param[in] mass_other The mass of the colliding species
 */
template <SpecialParam param>
void relaxation_matrix_offdiagonal(MatrixView W,
                                   const Numeric T,
                                   const AbsorptionLines& band,
                                   const ArrayOfIndex& sorting,
                                   const Numeric& mass_self,
                                   const Numeric& mass_other)
{
  using Conversion::kelvin2joule;
  
  const Index n = band.NumLines();
  
  Vector dip0(n);
  for (Index i=0; i<n; i++) {
    dip0[i] = zero_dipole(band.UpperQuantumNumber(sorting[i], QuantumNumberType::J),
                          band.LowerQuantumNumber(sorting[i], QuantumNumberType::J),
                          band.UpperQuantumNumber(sorting[i], QuantumNumberType::N));
    
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
      //    4) Use wigner6j(1, 1, 1, Jf, Ji, N) instead of wigner6j(1, 1, 1, Ji, Jf, N), which is in the paper
      // These are modified after reading the original code
      Numeric sum = 0;
      const Numeric scale = Numeric((2*Nk + 1) * (2*Nl + 1)) * sqrt((2*Jk + 1) * (2*Jl + 1) * (2*Jk_p + 1) * (2* Jl_p + 1));
      const auto [L0, L1] = wigner_limits(wigner3j_limits<3>(Nl, Nk), {Rational(2), 100000});
      for (Rational L=L0; L<=L1; L+=2) {
        const Numeric sgn = iseven(Jk + Jl + L + 1) ? 1 : -1;
        const Numeric a = Constant::pow2(wigner3j(Nl, Nk, L, 0, 0, 0));
        const Numeric b = wigner6j(L, Jk, Jl, 1, Nl, Nk);
        const Numeric c = wigner6j(L, Jk_p, Jl_p, 1, Nl, Nk);
        const Numeric d = wigner6j(L, Jk, Jl, 1, Jl_p, Jk_p);
        sum += sgn * a * b * c * d * Q<param>(L, T) / Omega<param>(L, T, mass_self, mass_other);
      }
      sum *= scale * Omega<param>(Nk, T, mass_self, mass_other);
      
      // Add to W and rescale to upwards element by the populations
      W(l, k) = sum;
      W(k, l) = sum * std::exp((erot(Nl, Jl) - erot(Nk, Jk)) / Conversion::kelvin2joule(T));
    }
  }
  
  // Transpose?  Why is this transpose required?
  for (Index i=0; i<n; i++) {
    for (Index j=0; j<i; j++) {
      swap(W(i, j), W(j, i));
    }
  }
  
  // Sum rule correction
  for (Index i=0; i<n; i++) {
    Numeric sumlw = 0.0;
    Numeric sumup = 0.0;
    
    for (Index j=0; j<n; j++) {
      if (j > i) {
        sumlw += dip0[j] * W(j, i);
      } else {
        sumup += dip0[j] * W(j, i);
      }
    }
    
    const Rational Ji = band.UpperQuantumNumber(sorting[i], QuantumNumberType::J);
    const Rational Ni = band.UpperQuantumNumber(sorting[i], QuantumNumberType::N);
    for (Index j=i+1; j<n; j++) {
      const Rational Jj = band.UpperQuantumNumber(sorting[j], QuantumNumberType::J);
      const Rational Nj = band.UpperQuantumNumber(sorting[j], QuantumNumberType::N);
      if (sumlw == 0) {
        W(j, i) = 0.0;
        W(i, j) = 0.0;
      } else {
        W(j, i) *= - sumup / sumlw;
        W(i, j) = W(j, i) * std::exp((erot(Ni, Ji) - erot(Nj, Jj)) / Conversion::kelvin2joule(T));  // This gives LTE
      }
    }
  }
}
}  // Makarov2020etal


/*! Computes the Error Corrected Sudden relaxation matrix for a single species
 * 
 * The band population type is used to select type of off-diagonal element computations
 * 
 * @param[in] band The absorption band
 * @param[in] sorting The sorting of the band
 * @param[in] T The temperature
 * @param[in] P The pressure
 * @param[in] species_mass The mass of the colliding species
 * @param[in] species_pos The index position of the colliding species among the band broadeners
 * @return The single species relaxation matrix
 */
template <SpecialParam param>
ComplexMatrix single_species_ecs_relaxation_matrix(const AbsorptionLines& band,
                                                   const ArrayOfIndex& sorting,
                                                   const Numeric T,
                                                   const Numeric P,
                                                   const Numeric species_mass,
                                                   const Index species_pos) {
  const Index N = band.NumLines();
  
  // Allocate the matrix
  ComplexMatrix W(N, N, 0);
  
  // Fill diagonal keeping track of the reshuffle (sorting)
  for (Index I=0; I<N; I++) {
    const Index i = sorting[I];
    auto shape = band.ShapeParameters(i, T, P, species_pos);
    W(I, I) = Complex(shape.D0, shape.G0);  // nb. This must be fixed if SDVP ???
  }
  
  // Set the off-diagonal part of the matrix for this broadener
  switch (band.Population()) {
    case PopulationType::ByMakarovFullRelmat:
      Makarov2020etal::relaxation_matrix_offdiagonal<param>(W.imag(), T, band, sorting, band.SpeciesMass(), species_mass);
      break;
    default:
      ARTS_ASSERT ("Bad type, we don't support band population type: ", band.Population(),
                   "\nin this code.  It must either be added or computations aborted earlier");
      break;
  }
  
  return W;
}


/*! Computes the Error Corrected Sudden relaxation matrix
 * 
 * @param[in] T The temperature
 * @param[in] P The pressure
 * @param[in] vmrs The VMRs of all broadeners of the absorption band
 * @param[in] mass The mass of all broadeners of the absorption band
 * @param[in] band The absorption band
 * @param[in] sorting The sorting of the band
 * @param[in] frenorm The renormalization of frequency
 * @return The relaxation matrix
 */
template <SpecialParam param>
ComplexMatrix ecs_relaxation_matrix(const Numeric T,
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
  for (Index k=0; k<M; k++) {
    // Create temporary
    const ComplexMatrix Wtmp = single_species_ecs_relaxation_matrix<param>(band, sorting, T, P, mass[k], k);
    
    // Sum up all atmospheric components
    for (Index i=0; i<N; i++) {
      for (Index j=0; j<N; j++) {
        W(i, j) += vmrs[k] * Wtmp(i, j);
      }
    }
  }
  
  // Deal with line frequency and its re-normalization
  for (Index i=0; i<N; i++) {
    W(i, i) += band.F0(sorting[i]) - frenorm;
  }
  
  return W;
}


/*! Returns sorted population distribtions and dipoles and the original sorting
 * 
 * @param[in] T The temperature
 * @param[in] band The absorption band
 * @return {sorting, sorted population distribtions and dipoles}
 */
std::pair<ArrayOfIndex, PopulationAndDipole> sorted_population_and_dipole(const Numeric T,
                                                                          const AbsorptionLines& band) {
  PopulationAndDipole tp(T, band);
  return {tp.sort(band), tp};
}


/*! Returns pre-sorted population distribtions and dipoles
 * 
 * @param[in] T The temperature
 * @param[in] presorting Changes positions from [0 ... N-1] to presorting positions
 * @param[in] band The absorption band
 * @return Pre-sorted  population distribtions and dipoles
 */
PopulationAndDipole presorted_population_and_dipole(const Numeric T,
                                                    const ArrayOfIndex& presorting,
                                                    const AbsorptionLines& band) {
  PopulationAndDipole tp(T, band);
  tp.sort(presorting);
  return tp;
}


template <SpecialParam param>
ComplexVector ecs_absorption_impl(const Numeric T,
                                  const Numeric P,
                                  const Numeric this_vmr,
                                  const Vector& vmrs,
                                  const Vector& mass,
                                  const Vector& f_grid,
                                  const AbsorptionLines& band) {
  constexpr Numeric sq_ln2pi = Constant::sqrt_ln_2 / Constant::sqrt_pi;
  
  // Weighted center of the band
  const Numeric frenorm = band.F_mean();
  
  // Band Doppler broadening constant
  const Numeric GD_div_F0 = Linefunctions::DopplerConstant(T, band.SpeciesMass());
  
  // Sorted population
  auto [sorting, tp] = sorted_population_and_dipole(T, band);
  
  // Relaxation matrix
  const ComplexMatrix W = ecs_relaxation_matrix<param>(T, P, vmrs, mass, band, sorting, frenorm);
  
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
  for (Index iv=0; iv<f_grid.nelem(); iv++) {
    const Numeric f = f_grid[iv];
    const Numeric fact = - f * std::expm1(- (Constant::h * f) / (Constant::k * T));
    absorption[iv] *= fact * numdens * sq_ln2pi;
  }
  
  return absorption;
}


std::pair<ComplexVector, ArrayOfComplexVector> ecs_absorption(const Numeric T,
                                                              const Numeric P,
                                                              const Numeric this_vmr,
                                                              const Vector& vmrs,
                                                              const Vector& mass,
                                                              const Vector& f_grid,
                                                              const AbsorptionLines& band,
                                                              const ArrayOfRetrievalQuantity& jacobian_quantities) {
  const ComplexVector absorption = ecs_absorption_impl<SpecialParam::None>(T, P, this_vmr, vmrs, mass, f_grid, band);
  
  // Start as original, so remove new and divide with the negative to get forward derivative
  ArrayOfComplexVector jacobian(jacobian_quantities.nelem(), absorption);
  
  for (Index i=0; i<jacobian_quantities.nelem(); i++) {
    auto& vec = jacobian[i];
    auto& target = jacobian_quantities[i].Target();
    
    if (target == Jacobian::Atm::Temperature) {
      const Numeric dT = target.Perturbation();
      vec -= ecs_absorption_impl<SpecialParam::None>(T+dT, P, this_vmr, vmrs, mass, f_grid, band);
      vec /= -dT;
    } else if (target.isWind()) {
      const Numeric df = target.Perturbation();
      Vector f_grid_copy = f_grid;
      f_grid_copy += df;
      vec -= ecs_absorption_impl<SpecialParam::None>(T, P, this_vmr, vmrs, mass, f_grid_copy, band);
    } else if (target == Jacobian::Line::VMR) {
      if (Absorption::QuantumIdentifierLineTarget(target.QuantumIdentity(), band) == Absorption::QuantumIdentifierLineTargetType::Isotopologue) {
        Vector vmrs_copy = vmrs;
        Numeric this_vmr_copy = this_vmr;
        const Numeric dvmr = target.Perturbation();
        
        // Alter the VMRs for self 
        if (band.Species() == target.QuantumIdentity().Species() and
          band.Isotopologue() == target.QuantumIdentity().Isotopologue()) {
          this_vmr_copy += dvmr;
          if (band.Self()) vmrs_copy[0] += dvmr;  // First value is self if band has self broadener
        } else {
          for (Index j=band.Self(); j<band.BroadeningSpecies().nelem()-band.Bath(); j++) {
            if (band.BroadeningSpecies()[j] == target.QuantumIdentity().Species()) {
              vmrs_copy[j] += dvmr;
            }
          }
        }
        
        // Computations
        vec -= ecs_absorption_impl<SpecialParam::None>(T, P, this_vmr_copy, vmrs_copy, mass, f_grid, band);
        vec /= -dvmr;
      }
    } else if (target.needQuantumIdentity()) {
      Numeric d=1e-6;
      
      for (Index iline=0; iline<band.NumLines(); iline++) {
        const Absorption::QuantumIdentifierLineTarget qlt(target.QuantumIdentity(), band, iline);
        if (qlt == Absorption::QuantumIdentifierLineTargetType::Line) {
          AbsorptionLines band_copy = band;
          
          const Index pos = band.BroadeningSpeciesPosition(target.QuantumIdentity().Species());
          switch (target.LineType()) {
            case Jacobian::Line::ShapeG0X0:
              d *= band.AllLines()[iline].LineShape()[pos].G0().X0;
              band_copy.AllLines()[iline].LineShape()[pos].G0().X0 += d;
              break;
            case Jacobian::Line::ShapeG0X1:
              d *= band.AllLines()[iline].LineShape()[pos].G0().X1;
              band_copy.AllLines()[iline].LineShape()[pos].G0().X1 += d;
              break;
            case Jacobian::Line::ShapeG0X2:
              d *= band.AllLines()[iline].LineShape()[pos].G0().X2;
              band_copy.AllLines()[iline].LineShape()[pos].G0().X2 += d;
              break;
            case Jacobian::Line::ShapeG0X3:
              d *= band.AllLines()[iline].LineShape()[pos].G0().X3;
              band_copy.AllLines()[iline].LineShape()[pos].G0().X3 += d;
              break;
            case Jacobian::Line::ShapeD0X0:
              d *= band.AllLines()[iline].LineShape()[pos].D0().X0;
              band_copy.AllLines()[iline].LineShape()[pos].D0().X0 += d;
              break;
            case Jacobian::Line::ShapeD0X1:
              d *= band.AllLines()[iline].LineShape()[pos].D0().X1;
              band_copy.AllLines()[iline].LineShape()[pos].D0().X1 += d;
              break;
            case Jacobian::Line::ShapeD0X2:
              d *= band.AllLines()[iline].LineShape()[pos].D0().X2;
              band_copy.AllLines()[iline].LineShape()[pos].D0().X2 += d;
              break;
            case Jacobian::Line::ShapeD0X3:
              d *= band.AllLines()[iline].LineShape()[pos].D0().X3;
              band_copy.AllLines()[iline].LineShape()[pos].D0().X3 += d;
              break;
            case Jacobian::Line::ShapeG2X0:
              d *= band.AllLines()[iline].LineShape()[pos].G2().X0;
              band_copy.AllLines()[iline].LineShape()[pos].G2().X0 += d;
              break;
            case Jacobian::Line::ShapeG2X1:
              d *= band.AllLines()[iline].LineShape()[pos].G2().X1;
              band_copy.AllLines()[iline].LineShape()[pos].G2().X1 += d;
              break;
            case Jacobian::Line::ShapeG2X2:
              d *= band.AllLines()[iline].LineShape()[pos].G2().X2;
              band_copy.AllLines()[iline].LineShape()[pos].G2().X2 += d;
              break;
            case Jacobian::Line::ShapeG2X3:
              d *= band.AllLines()[iline].LineShape()[pos].G2().X3;
              band_copy.AllLines()[iline].LineShape()[pos].G2().X3 += d;
              break;
            case Jacobian::Line::ShapeD2X0:
              d *= band.AllLines()[iline].LineShape()[pos].D2().X0;
              band_copy.AllLines()[iline].LineShape()[pos].D2().X0 += d;
              break;
            case Jacobian::Line::ShapeD2X1:
              d *= band.AllLines()[iline].LineShape()[pos].D2().X1;
              band_copy.AllLines()[iline].LineShape()[pos].D2().X1 += d;
              break;
            case Jacobian::Line::ShapeD2X2:
              d *= band.AllLines()[iline].LineShape()[pos].D2().X2;
              band_copy.AllLines()[iline].LineShape()[pos].D2().X2 += d;
              break;
            case Jacobian::Line::ShapeD2X3:
              d *= band.AllLines()[iline].LineShape()[pos].D2().X3;
              band_copy.AllLines()[iline].LineShape()[pos].D2().X3 += d;
              break;
            case Jacobian::Line::ShapeFVCX0:
              d *= band.AllLines()[iline].LineShape()[pos].FVC().X0;
              band_copy.AllLines()[iline].LineShape()[pos].FVC().X0 += d;
              break;
            case Jacobian::Line::ShapeFVCX1:
              d *= band.AllLines()[iline].LineShape()[pos].FVC().X1;
              band_copy.AllLines()[iline].LineShape()[pos].FVC().X1 += d;
              break;
            case Jacobian::Line::ShapeFVCX2:
              d *= band.AllLines()[iline].LineShape()[pos].FVC().X2;
              band_copy.AllLines()[iline].LineShape()[pos].FVC().X2 += d;
              break;
            case Jacobian::Line::ShapeFVCX3:
              d *= band.AllLines()[iline].LineShape()[pos].FVC().X3;
              band_copy.AllLines()[iline].LineShape()[pos].FVC().X3 += d;
              break;
            case Jacobian::Line::ShapeETAX0:
              d *= band.AllLines()[iline].LineShape()[pos].ETA().X0;
              band_copy.AllLines()[iline].LineShape()[pos].ETA().X0 += d;
              break;
            case Jacobian::Line::ShapeETAX1:
              d *= band.AllLines()[iline].LineShape()[pos].ETA().X1;
              band_copy.AllLines()[iline].LineShape()[pos].ETA().X1 += d;
              break;
            case Jacobian::Line::ShapeETAX2:
              d *= band.AllLines()[iline].LineShape()[pos].ETA().X2;
              band_copy.AllLines()[iline].LineShape()[pos].ETA().X2 += d;
              break;
            case Jacobian::Line::ShapeETAX3:
              d *= band.AllLines()[iline].LineShape()[pos].ETA().X3;
              band_copy.AllLines()[iline].LineShape()[pos].ETA().X3 += d;
              break;
            case Jacobian::Line::Center:
              d *= band.AllLines()[iline].F0();
              band_copy.AllLines()[iline].F0() += d;
              break;
            case Jacobian::Line::Strength:
              d *= band.AllLines()[iline].I0();
              band_copy.AllLines()[iline].I0() += d;
              break;
            case Jacobian::Line::ShapeYX0:
            case Jacobian::Line::ShapeYX1:
            case Jacobian::Line::ShapeYX2:
            case Jacobian::Line::ShapeYX3:
            case Jacobian::Line::ShapeGX0:
            case Jacobian::Line::ShapeGX1:
            case Jacobian::Line::ShapeGX2:
            case Jacobian::Line::ShapeGX3:
            case Jacobian::Line::ShapeDVX0:
            case Jacobian::Line::ShapeDVX1:
            case Jacobian::Line::ShapeDVX2:
            case Jacobian::Line::ShapeDVX3:
            case Jacobian::Line::NLTE:
            case Jacobian::Line::VMR:
            case Jacobian::Line::SpecialParameter1:
            case Jacobian::Line::SpecialParameter2:
            case Jacobian::Line::SpecialParameter3:
            case Jacobian::Line::FINAL: {
              /* do nothing */
            }
          }
          
          // Perform calculations and estimate derivative
          vec -= ecs_absorption_impl<SpecialParam::None>(T, P, this_vmr, vmrs, mass, f_grid, band_copy);
          vec /= -d;
        } else if (qlt == Absorption::QuantumIdentifierLineTargetType::Band) {
          if (target == Jacobian::Line::SpecialParameter1) {
            if (band.Population() == Absorption::PopulationType::ByMakarovFullRelmat) {
              d = Conversion::angstrom2meter(1e-6);
            }
            vec -= ecs_absorption_impl<SpecialParam::ModifyParameter1>(T, P, this_vmr, vmrs, mass, f_grid, band);
          } else if (target == Jacobian::Line::SpecialParameter2) {
            if (band.Population() == Absorption::PopulationType::ByMakarovFullRelmat) {
              d = 1e-6;
            }
            vec -= ecs_absorption_impl<SpecialParam::ModifyParameter2>(T, P, this_vmr, vmrs, mass, f_grid, band);
          } else if (target == Jacobian::Line::SpecialParameter3) {
            if (band.Population() == Absorption::PopulationType::ByMakarovFullRelmat) {
              d = 1e-6;
            }
            vec -= ecs_absorption_impl<SpecialParam::ModifyParameter3>(T, P, this_vmr, vmrs, mass, f_grid, band);
          }
          vec /= -d;
        } else {
          ARTS_ASSERT (false, "Missing Line Derivative");
        }
      }
    } else {
      vec *= 0;  // No derivative, so don't mess around and remove everything
    }
  }
  
  return {absorption, jacobian};
}


template <SpecialParam param>
ComplexVector ecs_absorption_zeeman_impl(const Numeric T,
                                         const Numeric H,
                                         const Numeric P,
                                         const Numeric this_vmr,
                                         const Vector& vmrs,
                                         const Vector& mass,
                                         const Vector& f_grid,
                                         const Zeeman::Polarization zeeman_polarization,
                                         const AbsorptionLines& band) {
  constexpr Numeric sq_ln2pi = Constant::sqrt_ln_2 / Constant::sqrt_pi;
  
  // Weighted center of the band
  const Numeric frenorm = band.F_mean();
  
  // Band Doppler broadening constant
  const Numeric GD_div_F0 = Linefunctions::DopplerConstant(T, band.SpeciesMass());
  
  // Sorted population
  const auto [sorting, tp] = sorted_population_and_dipole(T, band);
  
  // Relaxation matrix
  const ComplexMatrix W = ecs_relaxation_matrix<param>(T, P, vmrs, mass, band, sorting, frenorm);
  
  // Equivalent lines computations
  const EquivalentLines eqv(W, tp.pop, tp.dip);
  
  // Absorption of this band
  ComplexVector absorption(f_grid.nelem(), 0);
  for (Index i=0; i<band.NumLines(); i++) {
    // Zeeman lines if necessary
    const Index nz = band.ZeemanCount(i, zeeman_polarization);
    for (Index j=0; j<nz; j++) {
      const Numeric Sz = band.ZeemanStrength(i, zeeman_polarization, j);
      const Numeric dzeeman = H * band.ZeemanSplitting(i, zeeman_polarization, j);
      
      const Numeric gamd = GD_div_F0 * (eqv.val[i].real() + frenorm + dzeeman);
      const Numeric cte = Constant::sqrt_ln_2 / gamd;
      for (Index iv=0; iv<f_grid.nelem(); iv++) {
        const Complex z = (eqv.val[i] + frenorm + dzeeman - f_grid[iv]) * cte;
        const Complex w = Faddeeva::w(z);
        absorption[iv] += Sz * eqv.str[i] * w / gamd;
      }
    }
  }
  
  // Adjust by frequency and number density
  const Numeric numdens = this_vmr * number_density(P, T);
  for (Index iv=0; iv<f_grid.nelem(); iv++) {
    const Numeric f = f_grid[iv];
    const Numeric fact = - f * std::expm1(- (Constant::h * f) / (Constant::k * T));
    absorption[iv] *= fact * numdens * sq_ln2pi;
  }
  
  return absorption;
}


std::pair<ComplexVector, ArrayOfComplexVector> ecs_absorption_zeeman(const Numeric T,
                                                                     const Numeric H,
                                                                     const Numeric P,
                                                                     const Numeric this_vmr,
                                                                     const Vector& vmrs,
                                                                     const Vector& mass,
                                                                     const Vector& f_grid,
                                                                     const Zeeman::Polarization zeeman_polarization,
                                                                     const AbsorptionLines& band,
                                                                     const ArrayOfRetrievalQuantity& jacobian_quantities) {
  const ComplexVector absorption = ecs_absorption_zeeman_impl<SpecialParam::None>(T, H, P, this_vmr, vmrs, mass, f_grid, zeeman_polarization, band);
  
  // Start as original, so remove new and divide with the negative to get forward derivative
  ArrayOfComplexVector jacobian(jacobian_quantities.nelem(), absorption);
  
  for (Index i=0; i<jacobian_quantities.nelem(); i++) {
    auto& vec = jacobian[i];
    auto& target = jacobian_quantities[i].Target();
    
    if (target == Jacobian::Atm::Temperature) {
      const Numeric dT = target.Perturbation();
      vec -= ecs_absorption_zeeman_impl<SpecialParam::None>(T+dT, H, P, this_vmr, vmrs, mass, f_grid, zeeman_polarization, band);
      vec /= -dT;
    } else if (target.isMagnetic()) {
      const Numeric dH = target.Perturbation();
      vec -= ecs_absorption_zeeman_impl<SpecialParam::None>(T, H+dH, P, this_vmr, vmrs, mass, f_grid, zeeman_polarization, band);
      vec /= -dH;
    } else if (target.isWind()) {
      const Numeric df = target.Perturbation();
      Vector f_grid_copy = f_grid;
      f_grid_copy += df;
      vec -= ecs_absorption_zeeman_impl<SpecialParam::None>(T, H, P, this_vmr, vmrs, mass, f_grid_copy, zeeman_polarization, band);
    } else if (target == Jacobian::Line::VMR) {
      if (Absorption::QuantumIdentifierLineTarget(target.QuantumIdentity(), band) == Absorption::QuantumIdentifierLineTargetType::Isotopologue) {
        Vector vmrs_copy = vmrs;
        Numeric this_vmr_copy = this_vmr;
        const Numeric dvmr = target.Perturbation();
        
        // Alter the VMRs for self 
        if (band.Isotopologue() == target.QuantumIdentity().Isotopologue()) {
          this_vmr_copy += dvmr;
          if (band.Self()) vmrs_copy[0] += dvmr;  // First value is self if band has self broadener
        } else {
          for (Index j=band.Self(); j<band.BroadeningSpecies().nelem()-band.Bath(); j++) {
            if (band.BroadeningSpecies()[j] == target.QuantumIdentity().Species()) {
              vmrs_copy[j] += dvmr;
            }
          }
        }
        
        // Computations
        vec -= ecs_absorption_zeeman_impl<SpecialParam::None>(T, H, P, this_vmr_copy, vmrs_copy, mass, f_grid, zeeman_polarization, band);
        vec /= -dvmr;
      }
    } else if (target.needQuantumIdentity()) {
      Numeric d=1e-6;
      
      for (Index iline=0; iline<band.NumLines(); iline++) {
        const Absorption::QuantumIdentifierLineTarget qlt(target.QuantumIdentity(), band, iline);
        if (qlt == Absorption::QuantumIdentifierLineTargetType::Line) {
          AbsorptionLines band_copy = band;
          
          const Index pos = band.BroadeningSpeciesPosition(target.QuantumIdentity().Species());
          switch (target.LineType()) {
            case Jacobian::Line::ShapeG0X0:
              d *= band.AllLines()[iline].LineShape()[pos].G0().X0;
              band_copy.AllLines()[iline].LineShape()[pos].G0().X0 += d;
              break;
            case Jacobian::Line::ShapeG0X1:
              d *= band.AllLines()[iline].LineShape()[pos].G0().X1;
              band_copy.AllLines()[iline].LineShape()[pos].G0().X1 += d;
              break;
            case Jacobian::Line::ShapeG0X2:
              d *= band.AllLines()[iline].LineShape()[pos].G0().X2;
              band_copy.AllLines()[iline].LineShape()[pos].G0().X2 += d;
              break;
            case Jacobian::Line::ShapeG0X3:
              d *= band.AllLines()[iline].LineShape()[pos].G0().X3;
              band_copy.AllLines()[iline].LineShape()[pos].G0().X3 += d;
              break;
            case Jacobian::Line::ShapeD0X0:
              d *= band.AllLines()[iline].LineShape()[pos].D0().X0;
              band_copy.AllLines()[iline].LineShape()[pos].D0().X0 += d;
              break;
            case Jacobian::Line::ShapeD0X1:
              d *= band.AllLines()[iline].LineShape()[pos].D0().X1;
              band_copy.AllLines()[iline].LineShape()[pos].D0().X1 += d;
              break;
            case Jacobian::Line::ShapeD0X2:
              d *= band.AllLines()[iline].LineShape()[pos].D0().X2;
              band_copy.AllLines()[iline].LineShape()[pos].D0().X2 += d;
              break;
            case Jacobian::Line::ShapeD0X3:
              d *= band.AllLines()[iline].LineShape()[pos].D0().X3;
              band_copy.AllLines()[iline].LineShape()[pos].D0().X3 += d;
              break;
            case Jacobian::Line::ShapeG2X0:
              d *= band.AllLines()[iline].LineShape()[pos].G2().X0;
              band_copy.AllLines()[iline].LineShape()[pos].G2().X0 += d;
              break;
            case Jacobian::Line::ShapeG2X1:
              d *= band.AllLines()[iline].LineShape()[pos].G2().X1;
              band_copy.AllLines()[iline].LineShape()[pos].G2().X1 += d;
              break;
            case Jacobian::Line::ShapeG2X2:
              d *= band.AllLines()[iline].LineShape()[pos].G2().X2;
              band_copy.AllLines()[iline].LineShape()[pos].G2().X2 += d;
              break;
            case Jacobian::Line::ShapeG2X3:
              d *= band.AllLines()[iline].LineShape()[pos].G2().X3;
              band_copy.AllLines()[iline].LineShape()[pos].G2().X3 += d;
              break;
            case Jacobian::Line::ShapeD2X0:
              d *= band.AllLines()[iline].LineShape()[pos].D2().X0;
              band_copy.AllLines()[iline].LineShape()[pos].D2().X0 += d;
              break;
            case Jacobian::Line::ShapeD2X1:
              d *= band.AllLines()[iline].LineShape()[pos].D2().X1;
              band_copy.AllLines()[iline].LineShape()[pos].D2().X1 += d;
              break;
            case Jacobian::Line::ShapeD2X2:
              d *= band.AllLines()[iline].LineShape()[pos].D2().X2;
              band_copy.AllLines()[iline].LineShape()[pos].D2().X2 += d;
              break;
            case Jacobian::Line::ShapeD2X3:
              d *= band.AllLines()[iline].LineShape()[pos].D2().X3;
              band_copy.AllLines()[iline].LineShape()[pos].D2().X3 += d;
              break;
            case Jacobian::Line::ShapeFVCX0:
              d *= band.AllLines()[iline].LineShape()[pos].FVC().X0;
              band_copy.AllLines()[iline].LineShape()[pos].FVC().X0 += d;
              break;
            case Jacobian::Line::ShapeFVCX1:
              d *= band.AllLines()[iline].LineShape()[pos].FVC().X1;
              band_copy.AllLines()[iline].LineShape()[pos].FVC().X1 += d;
              break;
            case Jacobian::Line::ShapeFVCX2:
              d *= band.AllLines()[iline].LineShape()[pos].FVC().X2;
              band_copy.AllLines()[iline].LineShape()[pos].FVC().X2 += d;
              break;
            case Jacobian::Line::ShapeFVCX3:
              d *= band.AllLines()[iline].LineShape()[pos].FVC().X3;
              band_copy.AllLines()[iline].LineShape()[pos].FVC().X3 += d;
              break;
            case Jacobian::Line::ShapeETAX0:
              d *= band.AllLines()[iline].LineShape()[pos].ETA().X0;
              band_copy.AllLines()[iline].LineShape()[pos].ETA().X0 += d;
              break;
            case Jacobian::Line::ShapeETAX1:
              d *= band.AllLines()[iline].LineShape()[pos].ETA().X1;
              band_copy.AllLines()[iline].LineShape()[pos].ETA().X1 += d;
              break;
            case Jacobian::Line::ShapeETAX2:
              d *= band.AllLines()[iline].LineShape()[pos].ETA().X2;
              band_copy.AllLines()[iline].LineShape()[pos].ETA().X2 += d;
              break;
            case Jacobian::Line::ShapeETAX3:
              d *= band.AllLines()[iline].LineShape()[pos].ETA().X3;
              band_copy.AllLines()[iline].LineShape()[pos].ETA().X3 += d;
              break;
            case Jacobian::Line::Center:
              d *= band.AllLines()[iline].F0();
              band_copy.AllLines()[iline].F0() += d;
              break;
            case Jacobian::Line::Strength:
              d *= band.AllLines()[iline].I0();
              band_copy.AllLines()[iline].I0() += d;
              break;
            case Jacobian::Line::ShapeYX0:
            case Jacobian::Line::ShapeYX1:
            case Jacobian::Line::ShapeYX2:
            case Jacobian::Line::ShapeYX3:
            case Jacobian::Line::ShapeGX0:
            case Jacobian::Line::ShapeGX1:
            case Jacobian::Line::ShapeGX2:
            case Jacobian::Line::ShapeGX3:
            case Jacobian::Line::ShapeDVX0:
            case Jacobian::Line::ShapeDVX1:
            case Jacobian::Line::ShapeDVX2:
            case Jacobian::Line::ShapeDVX3:
            case Jacobian::Line::NLTE:
            case Jacobian::Line::VMR:
            case Jacobian::Line::SpecialParameter1:
            case Jacobian::Line::SpecialParameter2:
            case Jacobian::Line::SpecialParameter3:
            case Jacobian::Line::FINAL: {
              /* do nothing */
            }
          }
          
          // Perform calculations and estimate derivative
          vec -= ecs_absorption_zeeman_impl<SpecialParam::None>(T, H, P, this_vmr, vmrs, mass, f_grid, zeeman_polarization, band_copy);
          vec /= -d;
        } else if (qlt == Absorption::QuantumIdentifierLineTargetType::Band) {
          if (target == Jacobian::Line::SpecialParameter1) {
            if (band.Population() == Absorption::PopulationType::ByMakarovFullRelmat) {
              d = Conversion::angstrom2meter(1e-6);
            }
            vec -= ecs_absorption_zeeman_impl<SpecialParam::ModifyParameter1>(T, H, P, this_vmr, vmrs, mass, f_grid, zeeman_polarization, band);
          } else if (target == Jacobian::Line::SpecialParameter2) {
            if (band.Population() == Absorption::PopulationType::ByMakarovFullRelmat) {
              d = 1e-6;
            }
            vec -= ecs_absorption_zeeman_impl<SpecialParam::ModifyParameter2>(T, H, P, this_vmr, vmrs, mass, f_grid, zeeman_polarization, band);
          } else if (target == Jacobian::Line::SpecialParameter3) {
            if (band.Population() == Absorption::PopulationType::ByMakarovFullRelmat) {
              d = 1e-6;
            }
            vec -= ecs_absorption_zeeman_impl<SpecialParam::ModifyParameter3>(T, H, P, this_vmr, vmrs, mass, f_grid, zeeman_polarization, band);
          }
          vec /= -d;
        } else {
          ARTS_ASSERT (false, "Missing Line Derivative");
        }
      }
    } else {
      vec *= 0;  // No derivative, so don't mess around and remove everything
    }
  }
  
  return {absorption, jacobian};
}

Index band_eigenvalue_adaptation(
  AbsorptionLines& band,
  const Tensor4& tempdata,
  const Vector& temperatures,
  const Numeric P0,
  const Index ord)
{
  // Sizes
  const Index S = band.NumBroadeners();
  const Index N = band.NumLines();
  
  for (Index i=0; i<S; i++) {
    // Allocate Vector for fitting the data
    Vector targ(N);
    
    // Rosenkranz 1st order parameter
    for (Index j=0; j<N; j++) {
      targ = tempdata(1, j, i, joker);  // Get values
      
      // Assume linear pressure dependency
      targ /= P0;
      
      // Fit to a polynomial
      auto [fit, c] = Minimize::curve_fit<Minimize::Polynom>(temperatures, targ, LineShape::ModelParameters::N - 1);
      if (not fit) return EXIT_FAILURE;
      band.Line(j).LineShape()[i].Y() = LineShape::ModelParameters(LineShape::TemperatureModel::POLY, c);
    }
    
    // Rosenkranz 2nd order frequency parameter
    if (ord > 1) {
      for (Index j=0; j<N; j++) {
        targ = tempdata(2, j, i, joker);  // Get values
        
        // Assume squared pressure dependency
        targ /= P0 * P0;
        
        // Fit to a polynomial
        auto [fit, c] = Minimize::curve_fit<Minimize::Polynom>(temperatures, targ, LineShape::ModelParameters::N - 1);
        if (not fit) return EXIT_FAILURE;
        band.Line(j).LineShape()[i].DV() = LineShape::ModelParameters(LineShape::TemperatureModel::POLY, c);
      }
    }
    
    // Rosenkranz 2nd order strength parameter
    if (ord > 1) {
      for (Index j=0; j<N; j++) {
        targ = tempdata(0, j, i, joker);  // Get values
        
        // Assume squared pressure dependency
        targ /= P0 * P0;
        
        // Fit to a polynomial
        auto [fit, c] = Minimize::curve_fit<Minimize::Polynom>(temperatures, targ, LineShape::ModelParameters::N - 1);
        if (not fit) return EXIT_FAILURE;
        band.Line(j).LineShape()[i].G() = LineShape::ModelParameters(LineShape::TemperatureModel::POLY, c);
      }
    }
    
    // 3rd order broadening parameter  [[FIXME: UNTESTED]]
    if (ord > 2) {
      for (Index j=0; j<N; j++) {
        targ = tempdata(3, j, i, joker);
        
        // Assume cubic pressure dependency
        targ /= P0 * P0 * P0;
        
        // Fit to a polynomial
        auto [fit, c] = Minimize::curve_fit<Minimize::Polynom>(temperatures, targ, LineShape::ModelParameters::N - 1);
        if (not fit) return EXIT_FAILURE;
        ARTS_USER_ERROR("Not yet implemented: 3rd order line mixing")
        //         band.Line(j).LineShape()[i].DG() = LineShape::ModelParameters(LineShape::TemperatureModel::POLY, c);
      }
    }
  }
  
  // If we reach here, we have to set the band population type to LTE
  band.Population(Absorption::PopulationType::LTE);
  
  return EXIT_SUCCESS;
}


/*! Computes the Rosenkranz first order perturbation
 * 
 * @param[in] dip Reduced dipoles
 * @param[in] W The relaxation matrix
 * @param[in] band The absorption band
 * @return First order coefficients
 */
Vector RosenkranzY(const Vector& dip,
                   const ConstMatrixView& W,
                   const AbsorptionLines& band) {
  const Index N = dip.nelem();
  
  // Output
  Vector Y(N, 0);
  
  // Rosenkranz coefficients
  for (Index k=0; k<N; k++) {
    for (Index j=0; j<N; j++) {
      if (k == j) continue;
      else Y[k] += 2 * std::abs(dip[j] / dip[k]) * W(j, k) / (band.F0(k) - band.F0(j));
    }
  }
  
  return Y;
}


/*! Computes the Rosenkranz second order real perturbation
 * 
 * @param[in] dip Reduced dipoles
 * @param[in] W The relaxation matrix
 * @param[in] band The absorption band
 * @return Second order real coefficients
 */
Vector RosenkranzG(const Vector& dip,
                   const ConstMatrixView& W,
                   const AbsorptionLines& band) {
  const Index N = dip.nelem();
  
  // Output
  Vector G(N, 0);
  
  // Rosenkranz coefficients
  for (Index k=0; k<N; k++) {
    for (Index j=0; j<N; j++) {
      if (k == j) continue;
      else {
        G[k] += W(k, j) * W(j, k) / Constant::pow2(band.F0(j) - band.F0(k));
        G[k] += Constant::pow2(std::abs(dip[j] / dip[k]) * W(j, k) / (band.F0(j) - band.F0(k)));
        G[k] += 2 * std::abs(dip[j] / dip[k]) * W(j, k) * W(k, k) / Constant::pow2(band.F0(j) - band.F0(k));
        for (Index l=0; l<N; l++) {
          if (l == k or l == j) continue;
          else {
            G[k] -= 2 * std::abs(dip[j] / dip[k]) * W(j, l) * W(l, k) / ((band.F0(j) - band.F0(k)) * (band.F0(l) - band.F0(k)));
          }
        }
      }
    }
  }
  
  return G;
}


/*! Computes the Rosenkranz second order imaginary perturbation
 * 
 * @param[in] dip Reduced dipoles
 * @param[in] W The relaxation matrix
 * @param[in] band The absorption band
 * @return Second order imaginary coefficients
 */
Vector RosenkranzDV(const Vector& dip,
                    const ConstMatrixView& W,
                    const AbsorptionLines& band) {
  const Index N = dip.nelem();
  
  // Output
  Vector DV(N, 0);
  
  // Rosenkranz coefficients
  for (Index k=0; k<N; k++) {
    for (Index j=0; j<N; j++) {
      if (k == j) continue;
      else DV[k] += W(k, j) * W(j, k) / (band.F0(j) - band.F0(k));
    }
  }
  
  return DV;
}


EquivalentLines eigenvalue_adaptation_of_relmat(
  const ComplexMatrix& W,
  const Vector& pop,
  const Vector& dip,
  const AbsorptionLines& band,
  const Numeric frenorm,
  const Numeric T,
  const Numeric P,
  const Numeric QT,
  const Numeric QT0,
  const Index broadener) {
  // Compute and adapt values for the Voigt adaptation
  EquivalentLines eig(W, pop, dip);
  
  // Compute all the shifted line centers in band-order
  Vector f0(band.NumLines());
  for (Index i=0; i<pop.size(); i++) {
    f0[i] = band.F0(i) - frenorm + P * band.Line(i).LineShape()[broadener].compute(T, band.T0(), LineShape::Variable::D0);
  }
  
  // Find how the band-orders have shifted line order at the current temp/pres
  ArrayOfIndex sorting(band.NumLines());
  std::iota(sorting.begin(), sorting.end(), 0);
  for (Index i=0; i<band.NumLines(); i++) {
    for (Index j=i+1; j<band.NumLines(); j++) {
      if (f0[j] < f0[i]) {
        std::swap(sorting[i], sorting[j]);
        std::swap(f0[i], f0[j]);
      }
    }
  }
  
  // Sort the eigen-values and strengths in the same order as the band
  eig.sort_by_frequency(f0, sorting);
  
  // Eigenvalue should now be (F0 + P * D0(T) + i * P * G0(T)) + (P * P * DV(T) + i * P * P * P * DG(T)),
  // or in form (1) + (2).  We only want to keep (2), since it is new from the line mixing
  for (Index i=0; i<pop.size(); i++) {
    eig.val[i] -= Complex(f0[i], P * band.Line(i).LineShape()[broadener].compute(T, band.T0(), LineShape::Variable::G0));
  }
  
  // The strength term should now be (1 + y(T) + i * g(T)) * d**2 * rho(T)
  // d**2 * rho(T) * F0 * (1 - exp(-hF0/kT)) should be I0(T)
  // We only want to keep (y(T) + i * g(T)).
  // So we divide by d**2 * rho(T) through the use of I0(T) and remove 1 from the real component
  for (Index i=0; i<pop.size(); i++) {
    const Numeric i0 = Linefunctions::lte_linestrength(band.I0(i), band.E0(i), band.F0(i), QT0, band.T0(), QT, T);
    eig.str[i] *= - band.F0(i) * std::expm1(- (Constant::h * band.F0(i)) / (Constant::k * T)) / i0;
  }
  eig.str -= 1;
  
  return eig;
}


/*! Computes the Eigenvalue adaptation values
 * 
 * The output is sorted based on frequency.  If pressure shifts
 * moves one line's frequency past another, this introduces errors.
 * This is uncommon but happens more easily at higher pressures.
 * 
 * If this happens, the output of this function is not usable, but no
 * safe-guards are in place to guard against this.
 * 
 * @param[in] band The absorption band [N lines]
 * @param[in] temperatures The temperature grid for fitting parameters upon [K temperatures]
 * @param[in] mass The mass of all broadeners of the absorption band [M broadeners]
 * @return Eigenvalue line mixing parameters
 */
Tensor4 ecs_eigenvalue_approximation(const AbsorptionLines& band,
                                     const Vector& temperatures,
                                     const Vector& mass,
                                     const Numeric P) {
  const Index N = band.NumLines();
  const Index M = band.NumBroadeners();
  const Index K = temperatures.nelem();
  
  // Weighted center of the band
  const Numeric frenorm = band.F_mean();
  
  // Need sorting to put weak lines last, but we need the sorting constant or the output jumps
  const ArrayOfIndex sorting = sorted_population_and_dipole(band.T0(), band).first;
  const Numeric QT0 = single_partition_function(band.T0(), band.Isotopologue());
  
  // Output
  Tensor4 out(4, N, M, K);
  
  #pragma omp parallel for collapse(2) if (!arts_omp_in_parallel())
  for (Index m=0; m<M; m++) {
    for (Index k=0; k<K; k++) {
      const Numeric T = temperatures[k];
      const Numeric QT = single_partition_function(T, band.Isotopologue());
      
      // Relaxation matrix of T0 sorting at T
      ComplexMatrix W = single_species_ecs_relaxation_matrix<SpecialParam::None>(band, sorting, T, P, mass[m], m);
      for (Index n=0; n<N; n++) {
        W(n, n) += band.F0(sorting[n]) - frenorm;
      }
      
      // Populations and dipoles of T0 sorting at T
      const auto [pop, dip] = presorted_population_and_dipole(T, sorting, band);
      
      const auto eig = eigenvalue_adaptation_of_relmat(W, pop, dip, band, frenorm, T, P, QT, QT0, m);
      
      out(0, joker, m, k) = eig.str.real();
      out(1, joker, m, k) = eig.str.imag();
      out(2, joker, m, k) = eig.val.real();
      out(3, joker, m, k) = eig.val.imag();
    }
  }
  
  return out;
}

void ecs_eigenvalue_adaptation(AbsorptionLines& band,
                               const Vector& temperatures,
                               const Vector& mass,
                               const Numeric P0,
                               const Index ord
                               ) {
  ARTS_USER_ERROR_IF (P0 <= 0, P0, " Pa is not possible")
  
  ARTS_USER_ERROR_IF(
    not is_sorted(temperatures),
  "The temperature list [", temperatures,  "] K\n"
  "must be fully sorted from low to high"
  )
  
  ARTS_USER_ERROR_IF (ord < 1 or ord > 3, "Order not in list [1, 2, 3], is: ", ord)
  
  if (band_eigenvalue_adaptation(band,
    ecs_eigenvalue_approximation(band, temperatures, mass, P0),
    temperatures, P0, ord)) {
    ARTS_USER_ERROR("Bad eigenvalue adaptation")
  }
}


Tensor5 ecs_eigenvalue_adaptation_test(const AbsorptionLines& band,
                                       const Vector& temperatures,
                                       const Vector& mass,
                                       const Vector& pressures) {
  const Index N = band.NumLines();
  const Index M = band.NumBroadeners();
  const Index K = temperatures.nelem();
  const Index L = pressures.size();
  
  Tensor5 out(4, N, M, K, L);
  for (Index l=0; l<L; l++) {
    out(joker, joker, joker, joker, l) = ecs_eigenvalue_approximation(band, temperatures, mass, pressures[l]);
  }
  return out;
}
}  // Absorption::LineMixing

