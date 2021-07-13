#include <numeric>

#include <Faddeeva/Faddeeva.hh>

#include "lin_alg.h"
#include "linemixing.h"
#include "lineshape.h"
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

std::ostream& operator<<(std::ostream& os, const EquivalentLines& eqv) {
  const Index n = eqv.str.nelem();
  for (Index i=0; i<n; i++) {
    if (i) os  << '\n';
    os << eqv.str[i] << ' ' << eqv.val[i];
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const SpeciesErrorCorrectedSuddenData& srbd) {
  os << Species::toShortName(srbd.spec);
  os << ' ' << srbd.scaling;
  os << ' ' << srbd.beta;
  os << ' ' << srbd.lambda;
  os << ' ' << srbd.collisional_distance;
  os << ' ' << srbd.mass;
  return os;
}

std::istream& operator>>(std::istream& is, SpeciesErrorCorrectedSuddenData& srbd) {
  std::string spec_name;
  is >> spec_name;
  is >> srbd.scaling;
  is >> srbd.beta;
  is >> srbd.lambda;
  is >> srbd.collisional_distance;
  is >> double_imanip() >> srbd.mass;
  srbd.spec = Species::fromShortName(spec_name);
  ARTS_USER_ERROR_IF(not good_enum(srbd.spec), "Cannot recognize species: ", spec_name)
  return is;
}

namespace Makarov2020etal {
/*! Returns the reduced dipole
 * 
 * @param[in] Ju Main rotational number with spin of the upper level
 * @param[in] Jl Main rotational number with spin of the lower level
 * @param[in] N Main rotational number of both levels
 * @return The reduced dipole
 */
Numeric reduced_dipole(const Rational Ju, const Rational Jl, const Rational N) {
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
      dip[i] *= std::signbit(Makarov2020etal::reduced_dipole(band.UpperQuantumNumber(i, QuantumNumberType::J),
                                                             band.LowerQuantumNumber(i, QuantumNumberType::J),
                                                             band.UpperQuantumNumber(i, QuantumNumberType::N))) ? - 1 : 1;
    } else if (band.Population() == Absorption::PopulationType::ByRovibLinearDipoleLineMixing) {
      const Rational& Ji = band.UpperQuantumNumber(i, QuantumNumberType::J);
      const Rational& Jf = band.LowerQuantumNumber(i, QuantumNumberType::J);
      const Rational& li = band.UpperQuantumNumber(i, QuantumNumberType::l2);
      const Rational& lf = band.LowerQuantumNumber(i, QuantumNumberType::l2);
      dip[i] *= std::signbit((iseven(lf + Jf) ? 1 : -1) * sqrt(2* Jf + 1) * wigner3j(Ji, 1, Jf, li, lf-li, lf)) ? -1 : 1;
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

Numeric SpeciesErrorCorrectedSuddenData::Q(const Rational J,
                                           const Numeric T,
                                           const Numeric T0,
                                           const Numeric energy) const noexcept
{
  return std::exp(- beta.at(T, T0) * energy / (Constant::k * T)) * scaling.at(T, T0) / pow(J * (J+1), lambda.at(T, T0));
}

Numeric SpeciesErrorCorrectedSuddenData::Omega(const Numeric T,
                                               const Numeric T0,
                                               const Numeric other_mass,
                                               const Numeric energy_x,
                                               const Numeric energy_xm2) const noexcept
{
  using Constant::h;
  using Constant::k;
  using Constant::pi;
  using Constant::m_u;
  using Constant::h_bar;
  using Constant::pow2;
  
  // Constants for the expression
  constexpr Numeric fac = 8 * k / (m_u * pi);
  
  const Numeric wnnm2 = (energy_x - energy_xm2) / h_bar;
  
  const Numeric inv_eff_mass = 1 / mass + 1 / other_mass;
  const Numeric v_bar_pow2 = fac*T*inv_eff_mass;
  const Numeric tauc_pow2 = pow2(collisional_distance.at(T, T0)) / v_bar_pow2;
  
  return 1.0 / pow2(1 + 1.0/24.0 * pow2(wnnm2) * tauc_pow2);
}

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
void relaxation_matrix_offdiagonal(MatrixView W,
                                   const AbsorptionLines& band,
                                   const ArrayOfIndex& sorting,
                                   const SpeciesErrorCorrectedSuddenData& rovib_data,
                                   const Numeric T)
{
  using Conversion::kelvin2joule;
  
  const auto bk = [](const Rational& r) -> Numeric {return sqrt(2*r + 1);};
  
  const Index n = band.NumLines();
  
  Vector dipr(n);
  for (Index i=0; i<n; i++) {
    dipr[i] = reduced_dipole(band.UpperQuantumNumber(sorting[i], QuantumNumberType::J),
                             band.LowerQuantumNumber(sorting[i], QuantumNumberType::J),
                             band.UpperQuantumNumber(sorting[i], QuantumNumberType::N));
    
    for (Index j=i+1; j<n; j++) {
      // Select upper quantum number
      const bool ihigh = band.E0(sorting[i]) > band.E0(sorting[j]);
      const Index k = ihigh ? i : j;
      const Index l = ihigh ? j : i;
      
      // Quantum numbers
      const Rational Ji = band.UpperQuantumNumber(sorting[k], QuantumNumberType::J);
      const Rational Jf = band.LowerQuantumNumber(sorting[k], QuantumNumberType::J);
      const Rational Ni = band.UpperQuantumNumber(sorting[k], QuantumNumberType::N);
      const Rational Nf = band.LowerQuantumNumber(sorting[k], QuantumNumberType::N);
      const Rational Si = band.UpperQuantumNumber(sorting[k], QuantumNumberType::S);
      const Rational Sf = band.LowerQuantumNumber(sorting[k], QuantumNumberType::S);
      const Rational Ji_p = band.UpperQuantumNumber(sorting[l], QuantumNumberType::J);
      const Rational Jf_p = band.LowerQuantumNumber(sorting[l], QuantumNumberType::J);
      const Rational Ni_p = band.UpperQuantumNumber(sorting[l], QuantumNumberType::N);
      const Rational Nf_p = band.LowerQuantumNumber(sorting[l], QuantumNumberType::N);
      
      // Tran etal 2006 symbol with modifications:
      //    1) [Ji] * [Ji_p] instead of [Ji_p] ^ 2 in partial accordance with Makarov etal 2013
      Numeric sum=0;
      const Numeric scl = (iseven(Ji_p + Ji + 1) ? 1 : -1) * bk(Ni) * bk(Nf) * bk(Nf_p) * bk(Ni_p) * bk(Jf) * bk(Jf_p) * bk(Ji) * bk(Ji_p);
      const auto [L0, L1] = wigner_limits(wigner3j_limits<3>(Ni_p, Ni), {Rational(2), std::numeric_limits<Index>::max()});
      for (Rational L=L0; L<=L1; L+=2) {
        const Numeric a = wigner3j(Ni_p, Ni, L, 0, 0, 0);
        const Numeric b = wigner3j(Nf_p, Nf, L, 0, 0, 0);
        const Numeric c = wigner6j(L, Ji, Ji_p, Si, Ni_p, Ni);
        const Numeric d = wigner6j(L, Jf, Jf_p, Sf, Nf_p, Nf);
        const Numeric e = wigner6j(L, Ji, Ji_p, 1, Jf_p, Jf);
        sum += a * b * c * d * e * Numeric(2*L + 1) * rovib_data.Q(L, T, band.T0(), erot(L)) / rovib_data.Omega(T, band.T0(), band.SpeciesMass(), erot(L), erot(L-2));
      }
      sum *= scl * rovib_data.Omega(T, band.T0(), band.SpeciesMass(), erot(Ni), erot(Ni - 2));
      
      // Add to W and rescale to upwards element by the populations
      W(k, l) = sum;
      W(l, k) = sum * std::exp((erot(Nf_p, Jf_p) - erot(Nf, Jf)) / Conversion::kelvin2joule(T));
    }
  }
  
  // Sum rule correction
  for (Index i=0; i<n; i++) {
    Numeric sumlw = 0.0;
    Numeric sumup = 0.0;
    
    for (Index j=0; j<n; j++) {
      if (j > i) {
        sumlw += dipr[j] * W(j, i);
      } else {
        sumup += dipr[j] * W(j, i);
      }
    }
    
    const Rational Ji = band.LowerQuantumNumber(sorting[i], QuantumNumberType::J);
    const Rational Ni = band.LowerQuantumNumber(sorting[i], QuantumNumberType::N);
    for (Index j=i+1; j<n; j++) {
      const Rational Jj = band.LowerQuantumNumber(sorting[j], QuantumNumberType::J);
      const Rational Nj = band.LowerQuantumNumber(sorting[j], QuantumNumberType::N);
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


/** Generic interface for ECS calculations
 * 
 * Based on Collisional Effects On Molecular Spectra by
 * J.-M. Hartmann, C. Boulet, and D. Robert,
 * equation (IV.109; 1st edition)
 * 
 * It requires inputs for each of the internal components
 */
namespace LinearRovibErrorCorrectedSudden {
/** The entire "sum over L" full expression
 * 
 * The main transition is Jf   ← Ji
 * The virt transition is Jf_p ← Ji_p
 * 
 * @param[in] Ji Main upper state angular momentum quantum number
 * @param[in] Ji_p Virtual upper state angular momentum quantum number
 * @param[in] Jf Main lower state angular momentum quantum number
 * @param[in] Jf_p Virtual lower state angular momentum quantum number
 * @param[in] k 1 means IR absorption, 0 means isotropic Raman, and 2 means anisotropic Raman
 * @param[in] Q The Q function, must accept Q(L) and return a Numeric
 * @param[in] Omega The adiabatic function, must accept Omega(L) and return a Numeric
 * @return full sum result
 */
Numeric upper_offdiagonal_element(const Rational Ji,
                                  const Rational Ji_p,
                                  const Rational Jf,
                                  const Rational Jf_p,
                                  const Rational li,
                                  const Rational lf,
                                  const SpeciesErrorCorrectedSuddenData& rovib_data,
                                  const Rational k,
                                  const Numeric T,
                                  const Numeric T0,
                                  const Numeric self_mass,
                                  const EnergyFunction& erot) ARTS_NOEXCEPT
{
  ARTS_ASSERT(Ji <= Ji_p, "Ji is selected as the upper state in our formalism")
  
  Numeric sum = 0;
  const auto [Ls, Lf] = wigner_limits(wigner3j_limits<2>(Ji_p, Ji, li, 0, -li),
                                      {2, std::numeric_limits<Index>::max()});
  
  for (Rational L=Ls; L<=Lf; L+=2) {
    const Numeric a = wigner3j(Ji_p, L, Ji, li, 0, -li);
    const Numeric b = wigner3j(Jf_p, L, Jf, -lf, 0, lf);
    const Numeric c = wigner6j(Ji, Jf, k, Jf_p, Ji_p, L);
    sum += a * b * c * Numeric(2 * L + 1) * rovib_data.Q(L, T, T0, erot(L)) / rovib_data.Omega(T, T0, self_mass, erot(L), erot(L-2));
  }
  
  return (iseven(li + lf + k) ? -1 : 1) * rovib_data.Omega(T, T0, self_mass, erot(Ji), erot(Ji-2)) *
          Numeric(2 * Ji_p + 1) * sqrt((2 * Jf + 1) * (2 * Jf_p + 1)) * sum;
}

/** Compute off-diagonal elements for a pair of quantum numbers
 * 
 * One transition is Jf1 ← Ji1, and
 * one transition is Jf2 ← Ji2
 * 
 * Selects the largest Ji* (or Ji1 if both are equal) to be the 
 *'reference' off-diagonal element.  The other is going to be a
 * scale of the population ratios of the reference element
 * 
 * @param[in] Ji1 First transition's upper state angular momentum quantum number
 * @param[in] Ji2 Second transition's upper state angular momentum quantum number
 * @param[in] Jf1 First transition's lower state angular momentum quantum number
 * @param[in] Jf2 Second transition's lower state angular momentum quantum number
 * @param[in] k 1 means IR absorption, 0 means isotropic Raman, and 2 means anisotropic Raman
 * @param[in] pop1 First transition's population ratio
 * @param[in] pop2 Second transition's population ratio
 * @param[in] Q The Q function, must accept Q(L) and return a Numeric
 * @param[in] Omega The adiabatic function, must accept Omega(L) and return a Numeric
 * @return The pair of off-diagonal elements W12 and W21 in that order
 */
std::pair<Numeric, Numeric> offdiagonal_elements(const Rational Ji1,
                                                 const Rational Ji2,
                                                 const Rational Jf1,
                                                 const Rational Jf2,
                                                 const Rational li,
                                                 const Rational lf,
                                                 const SpeciesErrorCorrectedSuddenData& rovib_data,
                                                 const Rational k,
                                                 const Numeric pop12,
                                                 const Numeric T,
                                                 const Numeric T0,
                                                 const Numeric self_mass,
                                                 const EnergyFunction& erot) ARTS_NOEXCEPT
{
  if (Ji1 <= Ji2) {
    const Numeric W = upper_offdiagonal_element(Ji1, Ji2, Jf1, Jf2, li, lf, rovib_data,
                                                k, T, T0, self_mass, erot);
    return {W, W / pop12};
  } else {
    const Numeric W = upper_offdiagonal_element(Ji2, Ji1, Jf2, Jf1, li, lf, rovib_data,
                                                k, T, T0, self_mass, erot);
    return {W * pop12, W};
  }
}

/** Returns the off-diagonal relaxation matrix if full
 * 
 * @param[in] band The absorption band
 * @param[in] pop The population level distribtion
 * @param[in] sorting pop[i] corresponds to absorption line band.Lines(sorting[i])
 * @param[in] Q The Q function, must accept Q(L) and return a Numeric
 * @param[in] Omega The adiabatic function, must accept Omega(L) and return a Numeric
 * @return W but imaginary parts and diagonal elements are zero
 */
void real_offdiagonal_relaxation_matrix(MatrixView W,
                                        const AbsorptionLines& band,
                                        const ArrayOfIndex& sorting,
                                        const SpeciesErrorCorrectedSuddenData& rovib_data,
                                        const Rational k,
                                        const Numeric T,
                                        const EnergyFunction& erot) ARTS_NOEXCEPT
{
  const Index N = band.NumLines();
  ARTS_ASSERT(sorting.nelem() == N, "Inconsistent sorting and lines count")
  
  const Rational lf = band.LowerQuantumNumber(0, QuantumNumberType::l2);
  const Rational li = band.UpperQuantumNumber(0, QuantumNumberType::l2);
  
  for (Index i=0; i<N; i++) {
    for (Index j=0; j<i; j++) {
      const Rational Jf1 = band.LowerQuantumNumber(sorting[i], QuantumNumberType::J);
      const Rational Jf2 = band.LowerQuantumNumber(sorting[j], QuantumNumberType::J);
      const Rational Ji1 = band.UpperQuantumNumber(sorting[i], QuantumNumberType::J);
      const Rational Ji2 = band.UpperQuantumNumber(sorting[j], QuantumNumberType::J);
      
      const auto [W12, W21] = offdiagonal_elements(Ji1, Ji2, Jf1, Jf2, li, lf, rovib_data, k,
                                                   std::exp((erot(Jf2) - erot(Jf1)) / Conversion::kelvin2joule(T)),
                                                   T, band.T0(), band.SpeciesMass(), erot);
      W(j, i) = W12;
      W(i, j) = W21;
    }
  }
}

/** Readjust the real part of the relaxation matrix
 * to fulfill the sum rule
 * 
 * @param[in,out] W The real part of the relaxation matrix
 * @param[in] popr The relative population ratios
 * @param[in] dipr The reduced dipole moments
 */
void verify_sum_rule(MatrixView W,
                     const AbsorptionLines& band,
                     const Vector& dipr,
                     const ArrayOfIndex& sorting,
                     const Numeric T,
                     const EnergyFunction& erot) ARTS_NOEXCEPT
{
  const Index N = dipr.nelem();
  ARTS_ASSERT(W.nrows() == N and W.nrows() == W.ncols(), "Bad lines count and matrix size")
  
  // Sum rule correction
  for (Index i=0; i<N; i++) {
    Numeric sumlw = 0.0;
    Numeric sumup = 0.0;
    
    for (Index j=0; j<N; j++) {
      if (j > i) {
        sumlw += dipr[j] * W(j, i);
      } else {
        sumup += dipr[j] * W(j, i);
      }
    }
    
    for (Index j=i+1; j<N; j++) {
      if (sumlw == 0) {
        W(j, i) = 0.0;
        W(i, j) = 0.0;
      } else {
        W(j, i) *= - sumup / sumlw;
        W(i, j) = W(j, i) * std::exp((erot(band.LowerQuantumNumber(sorting[i], QuantumNumberType::J)) - erot(band.LowerQuantumNumber(sorting[j], QuantumNumberType::J))) / Conversion::kelvin2joule(T));  // This gives LTE
      }
    }
  }
}

EnergyFunction erot_selection(const SpeciesIsotopeRecord& isot)
{
  if (isot.spec == Species::Species::CarbonDioxide and isot.isotname == "626") {
    return [](const Rational J) -> Numeric {return Conversion::kaycm2joule(0.39021) * Numeric(J * (J + 1)) - Conversion::kaycm2joule(0.39021);};
  }
  
  ARTS_USER_ERROR(isot.FullName(), " has no rotational energies in ARTS")
  return [](const Rational J) -> Numeric {return Numeric(J) * std::numeric_limits<Numeric>::signaling_NaN();};
}

ENUMCLASS(TypeOfTransition, unsigned char, RamanIsotropic, Dipole, RamanAnisotropic)
constexpr Rational getTypeOfTransition(TypeOfTransition t) noexcept {
  switch(t) {
    case TypeOfTransition::RamanIsotropic: return 0;
    case TypeOfTransition::Dipole: return 1;
    case TypeOfTransition::RamanAnisotropic: return 2;
    case TypeOfTransition::FINAL: {/* leave last */}
  }
  return -1;
}

void relaxation_matrix_offdiagonal(MatrixView W,
                                   const AbsorptionLines& band,
                                   const ArrayOfIndex& sorting,
                                   const SpeciesErrorCorrectedSuddenData& rovib_data,
                                   const Numeric T) ARTS_NOEXCEPT
{
  constexpr Rational k = getTypeOfTransition(TypeOfTransition::Dipole);
  const EnergyFunction erot = erot_selection(band.Isotopologue());
  
  const Index N = band.NumLines();
  Vector dip(N);
  for (Index i=0; i<N; i++) {
    const Rational Ji = band.UpperQuantumNumber(sorting[i], QuantumNumberType::J);
    const Rational li = band.UpperQuantumNumber(sorting[i], QuantumNumberType::l2);
    const Rational Jf = band.LowerQuantumNumber(sorting[i], QuantumNumberType::J);
    const Rational lf = band.LowerQuantumNumber(sorting[i], QuantumNumberType::l2);
    dip[i] = (iseven(lf + Jf) ? 1 : -1) * sqrt(2* Jf + 1) * wigner3j(Ji, k, Jf, li, lf-li, lf);
  }
  
  real_offdiagonal_relaxation_matrix(W, band, sorting, rovib_data, k,  T, erot);
  
  // FIXME: Use local pop and dip instead?
  verify_sum_rule(W, band, dip, sorting, T, erot);  // FIXME: make this use local ratio?
}
}

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
ComplexMatrix single_species_ecs_relaxation_matrix(const AbsorptionLines& band,
                                                   const ArrayOfIndex& sorting,
                                                   const Numeric T,
                                                   const Numeric P,
                                                   const SpeciesErrorCorrectedSuddenData& species_ecs_data,
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
      Makarov2020etal::relaxation_matrix_offdiagonal(W.imag(), band, sorting, species_ecs_data, T);
      break;
    case PopulationType::ByRovibLinearDipoleLineMixing: {
      LinearRovibErrorCorrectedSudden::relaxation_matrix_offdiagonal(W.imag(), band, sorting, species_ecs_data, T);
    } break;
    default:
      ARTS_ASSERT(false, "Bad type, we don't support band population type: ", band.Population(),
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
ComplexMatrix ecs_relaxation_matrix(const Numeric T,
                                    const Numeric P,
                                    const Vector& vmrs,
                                    const ErrorCorrectedSuddenData& ecs_data,
                                    const AbsorptionLines& band,
                                    const ArrayOfIndex& sorting,
                                    const Numeric frenorm) {
  const Index N = band.NumLines();
  const Index M = vmrs.nelem();
  
  // Create output
  ComplexMatrix W(N, N, 0);
  
  // Loop over all the broadeners
  for (Index k=0; k<M; k++) {
    if (vmrs[k] == 0) continue;
    
    // Create temporary
    const ComplexMatrix Wtmp = single_species_ecs_relaxation_matrix(band, sorting, T, P, ecs_data[band.BroadeningSpecies()[k]], k);
    
    // Sum up all atmospheric components
    MapToEigen(W).noalias() += vmrs[k] * MapToEigen(Wtmp);
  }
  
  // Deal with line frequency and its re-normalization
  for (Index i=0; i<N; i++) {
    W(i, i) += band.F0(sorting[i]) - frenorm;
  }
  
  return W;
}


ComplexVector ecs_absorption_impl(const Numeric T,
                                  const Numeric P,
                                  const Numeric this_vmr,
                                  const Vector& vmrs,
                                  const ErrorCorrectedSuddenData& ecs_data,
                                  const Vector& f_grid,
                                  const AbsorptionLines& band) {
  constexpr Numeric sq_ln2pi = Constant::sqrt_ln_2 / Constant::sqrt_pi;
  
  // Weighted center of the band
  const Numeric frenorm = band.F_mean();
  
  // Band Doppler broadening constant
  const Numeric GD_div_F0 = band.DopplerConstant(T);
  
  // Sorted population
  auto [sorting, tp] = sorted_population_and_dipole(T, band);
  
  // Relaxation matrix
  const ComplexMatrix W = ecs_relaxation_matrix(T, P, vmrs, ecs_data, band, sorting, frenorm);
  
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
    
    ARTS_USER_ERROR_IF(isnan(absorption[iv]) or absorption[iv].real() < 0,
                       "There's a bad value in the absorption profile.  The offending band\n"
                       "must not be fulfilling some of the conditions associated with\n"
                       "on-the-fly line mixing, the value is: ", absorption[iv], "\n\n"
                       "The full band metadata:\n", band.MetaData(), '\n',
                       "The full band data:\n", band)
  }
  
  return absorption;
}


std::pair<ComplexVector, ArrayOfComplexVector> ecs_absorption(const Numeric T,
                                                              const Numeric P,
                                                              const Numeric this_vmr,
                                                              const Vector& vmrs,
                                                              const ErrorCorrectedSuddenData& ecs_data,
                                                              const Vector& f_grid,
                                                              const AbsorptionLines& band,
                                                              const ArrayOfRetrievalQuantity& jacobian_quantities) {
  const ComplexVector absorption = ecs_absorption_impl(T, P, this_vmr, vmrs, ecs_data, f_grid, band);
  
  // Start as original, so remove new and divide with the negative to get forward derivative
  ArrayOfComplexVector jacobian(jacobian_quantities.nelem(), absorption);
  
  for (Index i=0; i<jacobian_quantities.nelem(); i++) {
    auto& vec = jacobian[i];
    auto& target = jacobian_quantities[i].Target();
    
    if (target == Jacobian::Atm::Temperature) {
      const Numeric dT = target.Perturbation();
      vec -= ecs_absorption_impl(T+dT, P, this_vmr, vmrs, ecs_data, f_grid, band);
      vec /= -dT;
    } else if (target.isWind()) {
      const Numeric df = target.Perturbation();
      Vector f_grid_copy = f_grid;
      f_grid_copy += df;
      vec -= ecs_absorption_impl(T, P, this_vmr, vmrs, ecs_data, f_grid_copy, band);
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
        vec -= ecs_absorption_impl(T, P, this_vmr_copy, vmrs_copy, ecs_data, f_grid, band);
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
            case Jacobian::Line::ECS_A:
            case Jacobian::Line::ECS_B:
            case Jacobian::Line::ECS_GAMMA:
            case Jacobian::Line::ECS_DC:
            case Jacobian::Line::FINAL: {
              /* do nothing */
            }
          }
          
          // Perform calculations and estimate derivative
          vec -= ecs_absorption_impl(T, P, this_vmr, vmrs, ecs_data, f_grid, band_copy);
          vec /= -d;
        } else if (qlt == Absorption::QuantumIdentifierLineTargetType::Band) {
          /*
          d = target.Perturbation();
          ErrorCorrectedSuddenData ecs_data_copy = ecs_data;
          
          if (target == Jacobian::Line::ECS_A) {
            ecs_data_copy[target.LineSpecies()].a += d;
            vec -= ecs_absorption_impl(T, P, this_vmr, vmrs, ecs_data_copy, f_grid, band);
          } else if (target == Jacobian::Line::ECS_B) {
            ecs_data_copy[target.LineSpecies()].b += d;
            vec -= ecs_absorption_impl(T, P, this_vmr, vmrs, ecs_data, f_grid, band);
          } else if (target == Jacobian::Line::ECS_GAMMA) {
            ecs_data_copy[target.LineSpecies()].gamma += d;
            vec -= ecs_absorption_impl(T, P, this_vmr, vmrs, ecs_data, f_grid, band);
          } else if (target == Jacobian::Line::ECS_DC) {
            ecs_data_copy[target.LineSpecies()].dc += d;
            vec -= ecs_absorption_impl(T, P, this_vmr, vmrs, ecs_data, f_grid, band);
          }
          vec /= -d;
          */
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

ComplexVector ecs_absorption_zeeman_impl(const Numeric T,
                                         const Numeric H,
                                         const Numeric P,
                                         const Numeric this_vmr,
                                         const Vector& vmrs,
                                         const ErrorCorrectedSuddenData& ecs_data,
                                         const Vector& f_grid,
                                         const Zeeman::Polarization zeeman_polarization,
                                         const AbsorptionLines& band) {
  constexpr Numeric sq_ln2pi = Constant::sqrt_ln_2 / Constant::sqrt_pi;
  
  // Weighted center of the band
  const Numeric frenorm = band.F_mean();
  
  // Band Doppler broadening constant
  const Numeric GD_div_F0 = band.DopplerConstant(T);
  
  // Sorted population
  const auto [sorting, tp] = sorted_population_and_dipole(T, band);
  
  // Relaxation matrix
  const ComplexMatrix W = ecs_relaxation_matrix(T, P, vmrs, ecs_data, band, sorting, frenorm);
  
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
    
    ARTS_USER_ERROR_IF(isnan(absorption[iv]) or absorption[iv].real() < 0,
                       "There's a bad value in the absorption profile.  The offending band\n"
                       "must not be fulfilling some of the conditions associated with\n"
                       "on-the-fly line mixing, the value is: ", absorption[iv], "\n\n"
                       "The full band metadata:\n", band.MetaData(), '\n',
                       "The full band data:\n", band)
  }
  
  return absorption;
}


std::pair<ComplexVector, ArrayOfComplexVector> ecs_absorption_zeeman(const Numeric T,
                                                                     const Numeric H,
                                                                     const Numeric P,
                                                                     const Numeric this_vmr,
                                                                     const Vector& vmrs,
                                                                     const ErrorCorrectedSuddenData& ecs_data,
                                                                     const Vector& f_grid,
                                                                     const Zeeman::Polarization zeeman_polarization,
                                                                     const AbsorptionLines& band,
                                                                     const ArrayOfRetrievalQuantity& jacobian_quantities) {
  const ComplexVector absorption = ecs_absorption_zeeman_impl(T, H, P, this_vmr, vmrs, ecs_data, f_grid, zeeman_polarization, band);
  
  // Start as original, so remove new and divide with the negative to get forward derivative
  ArrayOfComplexVector jacobian(jacobian_quantities.nelem(), absorption);
  
  for (Index i=0; i<jacobian_quantities.nelem(); i++) {
    auto& vec = jacobian[i];
    auto& target = jacobian_quantities[i].Target();
    
    if (target == Jacobian::Atm::Temperature) {
      const Numeric dT = target.Perturbation();
      vec -= ecs_absorption_zeeman_impl(T+dT, H, P, this_vmr, vmrs, ecs_data, f_grid, zeeman_polarization, band);
      vec /= -dT;
    } else if (target.isMagnetic()) {
      const Numeric dH = target.Perturbation();
      vec -= ecs_absorption_zeeman_impl(T, H+dH, P, this_vmr, vmrs, ecs_data, f_grid, zeeman_polarization, band);
      vec /= -dH;
    } else if (target.isWind()) {
      const Numeric df = target.Perturbation();
      Vector f_grid_copy = f_grid;
      f_grid_copy += df;
      vec -= ecs_absorption_zeeman_impl(T, H, P, this_vmr, vmrs, ecs_data, f_grid_copy, zeeman_polarization, band);
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
        vec -= ecs_absorption_zeeman_impl(T, H, P, this_vmr_copy, vmrs_copy, ecs_data, f_grid, zeeman_polarization, band);
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
            case Jacobian::Line::ECS_A:
            case Jacobian::Line::ECS_B:
            case Jacobian::Line::ECS_GAMMA:
            case Jacobian::Line::ECS_DC:
            case Jacobian::Line::FINAL: {
              /* do nothing */
            }
          }
          
          // Perform calculations and estimate derivative
          vec -= ecs_absorption_zeeman_impl(T, H, P, this_vmr, vmrs, ecs_data, f_grid, zeeman_polarization, band_copy);
          vec /= -d;
        } else if (qlt == Absorption::QuantumIdentifierLineTargetType::Band) {
          /*
          d = target.Perturbation();
          ErrorCorrectedSuddenData ecs_data_copy = ecs_data;
          
          if (target == Jacobian::Line::ECS_A) {
            ecs_data_copy[target.LineSpecies()].a += d;
            vec -= ecs_absorption_zeeman_impl(T, H, P, this_vmr, vmrs, ecs_data, f_grid, zeeman_polarization, band);
          } else if (target == Jacobian::Line::ECS_B) {
            ecs_data_copy[target.LineSpecies()].b += d;
            vec -= ecs_absorption_zeeman_impl(T, H, P, this_vmr, vmrs, ecs_data, f_grid, zeeman_polarization, band);
          } else if (target == Jacobian::Line::ECS_GAMMA) {
            ecs_data_copy[target.LineSpecies()].gamma += d;
            vec -= ecs_absorption_zeeman_impl(T, H, P, this_vmr, vmrs, ecs_data, f_grid, zeeman_polarization, band);
          } else if (target == Jacobian::Line::ECS_DC) {
            ecs_data_copy[target.LineSpecies()].dc += d;
            vec -= ecs_absorption_zeeman_impl(T, H, P, this_vmr, vmrs, ecs_data, f_grid, zeeman_polarization, band);
          }
          
          vec /= -d;
          */
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
  
  // If we reach here, we have to set the band population type
  // to LTE and renormalize the strength to physical frequency
  band.Normalization(Absorption::NormalizationType::SFS);
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
    f0[i] = band.F0(i) - frenorm + P * band.Line(i).LineShape()[broadener].D0().at(T, band.T0());
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
    eig.val[i] -= Complex(f0[i], P * band.Line(i).LineShape()[broadener].G0().at(T, band.T0()));
  }
  
  // The strength term should now be (1 + i y(T) + g(T)) * d**2 * rho(T)
  // d**2 * rho(T) * F0 * (1 - exp(-hF0/kT)) should be I0(T)
  // We only want to keep (i y(T) + g(T)).
  // So we divide by d**2 * rho(T) through the use of I0(T) and remove 1 from the real component
  for (Index i=0; i<pop.size(); i++) {
    const Numeric i0 = LineShape::LocalThermodynamicEquilibrium(band.I0(i), band.T0(), T, band.F0(i), band.E0(i), QT, QT0, 0, 1, 0, 0).S;
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
                                     const ErrorCorrectedSuddenData& ecs_data,
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
      ComplexMatrix W = single_species_ecs_relaxation_matrix(band, sorting, T, P, ecs_data[band.BroadeningSpecies()[m]], m);
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
                               const ErrorCorrectedSuddenData& ecs_data,
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
    ecs_eigenvalue_approximation(band, temperatures, ecs_data, P0),
    temperatures, P0, ord)) {
    ARTS_USER_ERROR("Bad eigenvalue adaptation")
  }
}


Tensor5 ecs_eigenvalue_adaptation_test(const AbsorptionLines& band,
                                       const Vector& temperatures,
                                       const ErrorCorrectedSuddenData& ecs_data,
                                       const Vector& pressures) {
  const Index N = band.NumLines();
  const Index M = band.NumBroadeners();
  const Index K = temperatures.nelem();
  const Index L = pressures.size();
  
  Tensor5 out(4, N, M, K, L);
  for (Index l=0; l<L; l++) {
    out(joker, joker, joker, joker, l) = ecs_eigenvalue_approximation(band, temperatures, ecs_data, pressures[l]);
  }
  return out;
}
}  // Absorption::LineMixing

