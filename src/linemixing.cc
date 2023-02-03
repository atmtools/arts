#include <iomanip>
#include <iostream>
#include <numeric>

#include <Faddeeva/Faddeeva.hh>
#include <vector>

#include "arts_conversions.h"
#include "debug.h"
#include "lin_alg.h"
#include "linemixing.h"
#include "lineshape.h"
#include "matpack_complex.h"
#include "matpack_data.h"
#include "matpack_eigen.h"
#include "messages.h"
#include "minimize.h"
#include "physics_funcs.h"
#include "quantum_numbers.h"
#include "rational.h"
#include "species.h"
#include "wigner_functions.h"

#if DO_FAST_WIGNER
#define WIGNER3 fw3jja6
#define WIGNER6 fw6jja
#else
#define WIGNER3 wig3jj
#define WIGNER6 wig6jj
#endif

Numeric wig3(const Rational& a,
             const Rational& b,
             const Rational& c,
             const Rational& d,
             const Rational& e,
             const Rational& f) noexcept {
  return WIGNER3(
      a.toInt(2), b.toInt(2), c.toInt(2), d.toInt(2), e.toInt(2), f.toInt(2));
}

Numeric wig6(const Rational& a,
             const Rational& b,
             const Rational& c,
             const Rational& d,
             const Rational& e,
             const Rational& f) noexcept {
  return WIGNER6(
      a.toInt(2), b.toInt(2), c.toInt(2), d.toInt(2), e.toInt(2), f.toInt(2));
}

#undef WIGNER3
#undef WIGNER6

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
  inv(V, V);
  for (Index i=0; i<n; i++) {
    Complex z(0, 0);
    for (Index j=0; j<n; j++) {
      z += pop[j] * dip[j] * V(i, j);
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
} // namespace Makarov2020etal


PopulationAndDipole::PopulationAndDipole(const Numeric T,
                                         const AbsorptionLines& band) :
  pop(band.NumLines()), dip(band.NumLines()) {
  const Index n = band.NumLines();
  
  const Numeric QT = single_partition_function(T, band.Isotopologue());
  const Numeric QT0 = single_partition_function(band.T0, band.Isotopologue());
  const Numeric ratiopart = QT0 / QT;
  
  for (Index i=0; i<n; i++) {
    const Numeric pop0 = (band.lines[i].gupp / QT0) * boltzman_factor(band.T0, band.lines[i].E0);
    pop[i] = pop0 * ratiopart * boltzman_ratio(T, band.T0, band.lines[i].E0);
    dip[i] = std::sqrt(- band.lines[i].I0/(pop0 * band.lines[i].F0 * std::expm1(- (Constant::h * band.lines[i].F0)) / (Constant::k * band.T0)));
    
    // Adjust the sign depending on type
    if (band.population == Absorption::PopulationType::ByMakarovFullRelmat) {
      auto& J = band.lines[i].localquanta.val[QuantumNumberType::J];
      auto& N = band.lines[i].localquanta.val[QuantumNumberType::N];
      dip[i] *= std::signbit(Makarov2020etal::reduced_dipole(J.upp(), J.low(), N.upp())) ? - 1 : 1;
    } else if (band.population == Absorption::PopulationType::ByRovibLinearDipoleLineMixing) {
      // pass
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
    s[i] = band.lines[i].F0 * pop[i] * Math::pow2(dip[i]);
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
                                           const Numeric energy) const {
  return std::exp(- beta.at(T, T0) * energy / (Constant::k * T)) * scaling.at(T, T0) / pow(J * (J+1), lambda.at(T, T0));
}

Numeric SpeciesErrorCorrectedSuddenData::Omega(const Numeric T,
                                               const Numeric T0,
                                               const Numeric other_mass,
                                               const Numeric energy_x,
                                               const Numeric energy_xm2) const {
  using Constant::h;
  using Constant::k;
  using Constant::pi;
  using Constant::m_u;
  using Constant::h_bar;
  using Math::pow2;
  
  // Constants for the expression
  constexpr Numeric fac = 8 * k / (m_u * pi);
  
  const Numeric wnnm2 = (energy_x - energy_xm2) / h_bar;
  
  const Numeric inv_eff_mass = 1 / mass + 1 / other_mass;
  const Numeric v_bar_pow2 = fac*T*inv_eff_mass;
  const Numeric tauc_pow2 = pow2(collisional_distance.at(T, T0)) / v_bar_pow2;
  
  return 1.0 / pow2(1 + pow2(wnnm2) * tauc_pow2 / 24.0);
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
    using Math::pow2;
    using Math::pow3;
    
    constexpr Numeric B0=43100.4425e0;
    constexpr Numeric D0=.145123e0;
    constexpr Numeric H0=3.8e-08;
    constexpr Numeric xl0=59501.3435e0;
    constexpr Numeric xg0=-252.58633e0;
    constexpr Numeric xl1=0.058369e0;
    constexpr Numeric xl2=2.899e-07;
    constexpr Numeric xg1=-2.4344e-04;
    constexpr Numeric xg2=-1.45e-09;
    
    const auto XN = Numeric(N);
    const Numeric XX = XN * (XN + 1);
    const Numeric xlambda=xl0+xl1*XX+xl2*pow2(XX);
    const Numeric xgama=xg0+xg1*XX+xg2*pow2(XX);
    const Numeric C1=B0*XX-D0*pow2(XX)+H0*pow3(XX);
    
    if (J < N) {
      if (N == 1)  // erot<false>(1, 0)
        return mhz2joule(C1 - (xlambda+B0*(2.*XN-1.)+xgama*XN));
      return mhz2joule(C1 - (xlambda+B0*(2.*XN-1.)+xgama*XN) + std::sqrt(pow2(B0*(2.*XN-1.))+pow2(xlambda)-2.*B0*xlambda));
    }
    if (J > N)
      return mhz2joule(C1 - (xlambda-B0*(2.*XN+3.)-xgama*(XN+1.)) - std::sqrt(pow2(B0*(2.*XN+3.))+pow2(xlambda)-2.*B0*xlambda));
    return mhz2joule(C1);
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
void relaxation_matrix_offdiagonal(
    MatrixView W,
    const AbsorptionLines& band,
    const ArrayOfIndex& sorting,
    const SpeciesErrorCorrectedSuddenData& rovib_data,
    const Numeric T) {
  using Conversion::kelvin2joule;

  const auto bk = [](const Rational& r) -> Numeric { return sqrt(2 * r + 1); };

  const Index n = band.NumLines();
  if (n == 0) return;

  auto& S = band.quantumidentity.val[QuantumNumberType::S];
  const Rational Si = S.upp();
  const Rational Sf = S.low();

  Vector dipr(n);
  for (Index i = 0; i < n; i++) {
    auto& J = band.lines[sorting[i]].localquanta.val[QuantumNumberType::J];
    auto& N = band.lines[sorting[i]].localquanta.val[QuantumNumberType::N];

    dipr[i] = reduced_dipole(J.upp(), J.low(), N.upp());
  }

  const auto maxL = temp_init_size(band.max(QuantumNumberType::J), band.max(QuantumNumberType::N));

  const auto Om = [&]() {
    std::vector<Numeric> out(maxL);
    for (Index i = 0; i < maxL; i++)
      out[i] = rovib_data.Omega(
          T, band.T0, band.SpeciesMass(), erot(i), erot(i - 2));
    return out;
  }();

  const auto Q = [&]() {
    std::vector<Numeric> out(maxL);
    for (Index i = 0; i < maxL; i++)
      out[i] = rovib_data.Q(i, T, band.T0, erot(i));
    return out;
  }();

  wig_thread_temp_init(maxL);
  for (Index i = 0; i < n; i++) {
    auto& J = band.lines[sorting[i]].localquanta.val[QuantumNumberType::J];
    auto& N = band.lines[sorting[i]].localquanta.val[QuantumNumberType::N];

    const Rational Ji = J.upp();
    const Rational Jf = J.low();
    const Rational Ni = N.upp();
    const Rational Nf = N.low();

    for (Index j = 0; j < n; j++) {
      if (i == j) continue;

      auto& J_p = band.lines[sorting[j]].localquanta.val[QuantumNumberType::J];
      auto& N_p = band.lines[sorting[j]].localquanta.val[QuantumNumberType::N];

      const Rational Ji_p = J_p.upp();
      const Rational Jf_p = J_p.low();
      const Rational Ni_p = N_p.upp();
      const Rational Nf_p = N_p.low();

      if (Jf_p > Jf) continue;

      // Tran etal 2006 symbol with modifications:
      //    1) [Ji] * [Ji_p] instead of [Ji_p] ^ 2 in partial accordance with Makarov etal 2013
      Numeric sum = 0;
      const Numeric scl = (iseven(Ji_p + Ji + 1) ? 1 : -1) * bk(Ni) * bk(Nf) *
                          bk(Nf_p) * bk(Ni_p) * bk(Jf) * bk(Jf_p) * bk(Ji) *
                          bk(Ji_p);
      const auto [L0, L1] =
          wigner_limits(wigner3j_limits<3>(Ni_p, Ni),
                        {Rational(2), std::numeric_limits<Index>::max()});
      for (Rational L = L0; L <= L1; L += 2) {
        const Numeric a = wig3(Ni_p, Ni, L, 0, 0, 0);
        const Numeric b = wig3(Nf_p, Nf, L, 0, 0, 0);
        const Numeric c = wig6(L, Ji, Ji_p, Si, Ni_p, Ni);
        const Numeric d = wig6(L, Jf, Jf_p, Sf, Nf_p, Nf);
        const Numeric e = wig6(L, Ji, Ji_p, 1, Jf_p, Jf);
        sum += a * b * c * d * e * Numeric(2 * L + 1) * Q[L.toIndex()] / Om[L.toIndex()];
      }
      sum *= scl * Om[Ni.toIndex()];

      // Add to W and rescale to upwards element by the populations
      W(i, j) = sum;
      W(j, i) = sum * std::exp((erot(Nf_p) - erot(Nf)) / kelvin2joule(T));
    }
  }
  wig_temp_free();

  ARTS_USER_ERROR_IF(errno == EDOM, "Cannot compute the wigner symbols")

  // Sum rule correction
  for (Index i = 0; i < n; i++) {
    Numeric sumlw = 0.0;
    Numeric sumup = 0.0;

    for (Index j = 0; j < n; j++) {
      if (j > i) {
        sumlw += dipr[j] * W(j, i);
      } else {
        sumup += dipr[j] * W(j, i);
      }
    }

    const Rational Ni =
        band.lines[sorting[i]].localquanta.val[QuantumNumberType::N].low();
    for (Index j = i + 1; j < n; j++) {
      const Rational Nj =
          band.lines[sorting[j]].localquanta.val[QuantumNumberType::N].low();
      if (sumlw == 0) {
        W(j, i) = 0.0;
        W(i, j) = 0.0;
      } else {
        W(j, i) *= -sumup / sumlw;
        W(i, j) = W(j, i) * std::exp((erot(Ni) - erot(Nj)) /
                                     kelvin2joule(T));  // This gives LTE
      }
    }
  }
}
} // namespace Makarov2020etal


/** Generic interface for ECS calculations
 * 
 * Based on Collisional Effects On Molecular Spectra by
 * J.-M. Hartmann, C. Boulet, and D. Robert,
 * equation (IV.109; 1st edition)
 * 
 * It requires inputs for each of the internal components
 */
namespace LinearRovibErrorCorrectedSudden {

EnergyFunction erot_selection(const SpeciesIsotopeRecord& isot) {
  if (isot.spec == Species::Species::CarbonDioxide and isot.isotname == "626") {
    return [](const Rational J) -> Numeric {return Conversion::kaycm2joule(0.39021) * Numeric(J * (J + 1));};
  }
  
  ARTS_USER_ERROR(isot.FullName(), " has no rotational energies in ARTS")
  return [](const Rational J) -> Numeric {return Numeric(J) * std::numeric_limits<Numeric>::signaling_NaN();};
}

void relaxation_matrix_offdiagonal(MatrixView W,
                                   const AbsorptionLines& band,
                                   const ArrayOfIndex& sorting,
                                   const SpeciesErrorCorrectedSuddenData& rovib_data,
                                   const Numeric T) {
  using Conversion::kelvin2joule;
  
  const Index n = band.NumLines();
  if (not n) return;

  // These are constant for a band
  auto& l2 = band.quantumidentity.val[QuantumNumberType::l2];
  Rational li = l2.upp();
  Rational lf = l2.low();
  
  const bool swap_order = li > lf;
  if (swap_order) swap(li, lf);
  const int sgn = iseven(li + lf + 1) ? -1 : 1;
  if (abs(li - lf) > 1) return;

  const EnergyFunction erot = erot_selection(band.Isotopologue());
  
  Vector dipr(n);
  for (Index i=0; i<n; i++) {
    auto& J = band.lines[i].localquanta.val[QuantumNumberType::J];
    dipr[i] = Absorption::reduced_rovibrational_dipole(J.upp(), J.low(), lf, li);
  }
  
  for (Index i=0; i<n; i++) {
    auto& J = band.lines[sorting[i]].localquanta.val[QuantumNumberType::J];
    Rational Ji = J.upp();
    Rational Jf = J.low();
    if (swap_order) swap(Ji, Jf);

    for (Index j=0; j<n; j++) {
      if(i == j) continue;
      auto& J_p = band.lines[sorting[j]].localquanta.val[QuantumNumberType::J];
      Rational Ji_p = J_p.upp();
      Rational Jf_p = J_p.low();
      if (swap_order) swap(Ji_p, Jf_p);

      // Select upper quantum number
      if (Jf_p > Jf) continue;

      Index L = std::max(std::abs((Ji-Ji_p).toIndex()), std::abs((Jf-Jf_p).toIndex()));
      L += L % 2;
      const Index Lf = std::min((Ji+Ji_p).toIndex(), (Jf+Jf_p).toIndex());
      
      Numeric sum=0;
      for (; L<=Lf; L+=2) {
        const Numeric a = wigner3j(Ji_p, L, Ji, li, 0, -li);
        const Numeric b = wigner3j(Jf_p, L, Jf, lf, 0, -lf);
        const Numeric c = wigner6j(Ji, Jf, 1, Jf_p, Ji_p, L);
        const Numeric QL = rovib_data.Q(L, T, band.T0, erot(L));
        const Numeric ECS = rovib_data.Omega(T, band.T0, band.SpeciesMass(), erot(L), erot(L-2));
        sum += a * b * c * Numeric(2 * L + 1) * QL / ECS;
      }
      const Numeric ECS = rovib_data.Omega(T, band.T0, band.SpeciesMass(), erot(Ji), erot(Ji-2));
      const Numeric scl = sgn * ECS * Numeric(2 * Ji_p + 1) * sqrt((2 * Jf + 1) * (2 * Jf_p + 1));
      sum *= scl;
      
      // Add to W and rescale to upwards element by the populations
      W(j, i) = sum;
      W(i, j) = sum * std::exp((erot(Jf_p) - erot(Jf)) / kelvin2joule(T));
    }
  }

  // Undocumented negative absolute sign 
  for (Index i=0; i<n; i++)
    for (Index j=0; j<n; j++)
      if (j not_eq i and W(i, j) > 0)
        W(i, j) *= -1;
  
  // Sum rule correction
  for (Index i=0; i<n; i++) {
    Numeric sumlw = 0.0;
    Numeric sumup = 0.0;
    
    for (Index j=0; j<n; j++) {
      if (j > i) {
        sumlw += std::abs(dipr[j]) * W(j, i);  // Undocumented abs-sign
      } else {
        sumup += std::abs(dipr[j]) * W(j, i);  // Undocumented abs-sign
      }
    }
    
    const Rational Ji =  band.lines[sorting[i]].localquanta.val[QuantumNumberType::J].low();
    for (Index j=i+1; j<n; j++) {
      const Rational Jj =  band.lines[sorting[j]].localquanta.val[QuantumNumberType::J].low();
      if (sumlw == 0) {
        W(j, i) = 0.0;
        W(i, j) = 0.0;
      } else {
        W(j, i) *= - sumup / sumlw;
        W(i, j) = W(j, i) * std::exp((erot(Ji) - erot(Jj)) / kelvin2joule(T));  // This gives LTE
      }
    }
  }
}
} // namespace LinearRovibErrorCorrectedSudden

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
  switch (band.population) {
    case PopulationType::ByMakarovFullRelmat:
      Makarov2020etal::relaxation_matrix_offdiagonal(W.imag(), band, sorting, species_ecs_data, T);
      break;
    case PopulationType::ByRovibLinearDipoleLineMixing: {
      LinearRovibErrorCorrectedSudden::relaxation_matrix_offdiagonal(W.imag(), band, sorting, species_ecs_data, T);
    } break;
    default:
      ARTS_ASSERT(false, "Bad type, we don't support band population type: ", band.population,
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
    
    // Sum up all atmospheric components
    matpack::eigen::mat(W).noalias() += vmrs[k] * matpack::eigen::mat(
      single_species_ecs_relaxation_matrix(band, sorting, T, P, ecs_data[band.broadeningspecies[k]], k));
  }

  // Deal with line frequency and its re-normalization
  for (Index i=0; i<N; i++) {
   real_val(W(i, i)) += band.lines[sorting[i]].F0 - frenorm;
  }

  return W;
}


std::pair<ComplexVector, bool> ecs_absorption_impl(const Numeric T,
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
  
  // Return the absorption and if this works
  std::pair<ComplexVector, bool> retval{ComplexVector(f_grid.nelem(), 0), false};
  auto& [absorption,works] = retval;

  //  Add the lines
  for (Index i=0; i<band.NumLines(); i++) {
    // Zeeman lines if necessary
    const Index nz = band.ZeemanCount(i, zeeman_polarization);
    for (Index j=0; j<nz; j++) {
      const Numeric Sz = band.ZeemanStrength(i, zeeman_polarization, j);
      const Numeric dzeeman = H * band.ZeemanSplitting(i, zeeman_polarization, j);
      
      if (band.lineshapetype == LineShape::Type::LP) {
        for (Index iv=0; iv<f_grid.nelem(); iv++) {
          absorption[iv] -= 1i * Sz * ((eqv.str[i] / (f_grid[iv] - frenorm - dzeeman - eqv.val[i]))) / (Constant::sqrt_ln_2 * Constant::sqrt_pi);
        }
      } else /*if (band.LineShapeType() == LineShape::Type::VP)*/ {
        const Numeric gamd = GD_div_F0 * (eqv.val[i].real() + frenorm + dzeeman);
        const Numeric cte = Constant::sqrt_ln_2 / gamd;
        for (Index iv=0; iv<f_grid.nelem(); iv++) {
          const Complex z = (eqv.val[i] + frenorm + dzeeman - f_grid[iv]) * cte;
          const Complex w = Faddeeva::w(z);
          absorption[iv] += Sz * eqv.str[i] * w / gamd;
        }
      }
    }
  }
  
  // Adjust by frequency and number density
  const Numeric numdens = this_vmr * number_density(P, T);
  for (Index iv=0; iv<f_grid.nelem(); iv++) {
    const Numeric f = f_grid[iv];
    const Numeric fact = - f * std::expm1(- (Constant::h * f) / (Constant::k * T));
    absorption[iv] *= fact * numdens * sq_ln2pi;

    // correct for too low value
    if (auto& a = real_val(absorption[iv]); not std::isnormal(a) or a < 0) {
      absorption[iv] = 0;
      works = false;
    }
  }

  return retval;
}


EcsReturn ecs_absorption(const Numeric T,
                         const Numeric H,
                         const Numeric P,
                         const Numeric this_vmr,
                         const Vector& vmrs,
                         const ErrorCorrectedSuddenData& ecs_data,
                         const Vector& f_grid,
                         const Zeeman::Polarization zeeman_polarization,
                         const AbsorptionLines& band,
                         const ArrayOfRetrievalQuantity& jacobian_quantities) {
  auto [absorption, work] = ecs_absorption_impl(T, H, P, this_vmr, vmrs, ecs_data, f_grid, zeeman_polarization, band);
  
  // Start as original, so remove new and divide with the negative to get forward derivative
  ArrayOfComplexVector jacobian(jacobian_quantities.nelem(), absorption);
  
  for (Index i=0; i<jacobian_quantities.nelem(); i++) {
    auto& vec = jacobian[i];
    auto& target = jacobian_quantities[i].Target();
    
    if (target == Jacobian::Atm::Temperature) {
      const Numeric dT = target.perturbation;
      const auto [dabs, dwork] = ecs_absorption_impl(T+dT, H, P, this_vmr, vmrs, ecs_data, f_grid, zeeman_polarization, band);
      vec -= dabs;
      vec /= -dT;
      work &= dwork;
    } else if (target.isMagnetic()) {
      const Numeric dH = target.perturbation;
      const auto [dabs, dwork] = ecs_absorption_impl(T, H+dH, P, this_vmr, vmrs, ecs_data, f_grid, zeeman_polarization, band);
      vec -= dabs;
      vec /= -dH;
      work &= dwork;
    } else if (target.isWind()) {
      const Numeric df = target.perturbation;
      Vector f_grid_copy = f_grid;
      f_grid_copy += df;
      const auto [dabs, dwork] = ecs_absorption_impl(T, H, P, this_vmr, vmrs, ecs_data, f_grid_copy, zeeman_polarization, band);
      vec -= dabs;
      vec /= df;
      work &= dwork;
    } else if (target == Jacobian::Line::VMR) {
      if (band.DoVmrDerivative(target.qid)) {
        Vector vmrs_copy = vmrs;
        Numeric this_vmr_copy = this_vmr;
        const Numeric dvmr = target.perturbation;
        
        // Alter the VMRs for self 
        if (band.Isotopologue() == target.qid.Isotopologue()) {
          this_vmr_copy += dvmr;
          if (band.selfbroadening) vmrs_copy[0] += dvmr;  // First value is self if band has self broadener
        } else {
          for (Index j=band.selfbroadening; j<band.broadeningspecies.nelem()-band.bathbroadening; j++) {
            if (band.broadeningspecies[j] == target.qid.Species()) {
              vmrs_copy[j] += dvmr;
            }
          }
        }
        
        // Computations
        const auto [dabs, dwork] = ecs_absorption_impl(T, H, P, this_vmr_copy, vmrs_copy, ecs_data, f_grid, zeeman_polarization, band);
        vec -= dabs;
        vec /= -dvmr;
        work &= dwork;
      }
    } else if (target.needQuantumIdentity()) {
      Numeric d=1e-6;
      
      for (Index iline=0; iline<band.NumLines(); iline++) {
        if (Quantum::Number::StateMatch(target.qid, band.lines[iline].localquanta, band.quantumidentity) == Quantum::Number::StateMatchType::Full) {
          AbsorptionLines band_copy = band;
          
          const Index pos = band.BroadeningSpeciesPosition(target.qid.Species());
          switch (target.line) {
            case Jacobian::Line::ShapeG0X0:
              d *= band.lines[iline].lineshape[pos].G0().X0;
              band_copy.lines[iline].lineshape[pos].G0().X0 += d;
              break;
            case Jacobian::Line::ShapeG0X1:
              d *= band.lines[iline].lineshape[pos].G0().X1;
              band_copy.lines[iline].lineshape[pos].G0().X1 += d;
              break;
            case Jacobian::Line::ShapeG0X2:
              d *= band.lines[iline].lineshape[pos].G0().X2;
              band_copy.lines[iline].lineshape[pos].G0().X2 += d;
              break;
            case Jacobian::Line::ShapeG0X3:
              d *= band.lines[iline].lineshape[pos].G0().X3;
              band_copy.lines[iline].lineshape[pos].G0().X3 += d;
              break;
            case Jacobian::Line::ShapeD0X0:
              d *= band.lines[iline].lineshape[pos].D0().X0;
              band_copy.lines[iline].lineshape[pos].D0().X0 += d;
              break;
            case Jacobian::Line::ShapeD0X1:
              d *= band.lines[iline].lineshape[pos].D0().X1;
              band_copy.lines[iline].lineshape[pos].D0().X1 += d;
              break;
            case Jacobian::Line::ShapeD0X2:
              d *= band.lines[iline].lineshape[pos].D0().X2;
              band_copy.lines[iline].lineshape[pos].D0().X2 += d;
              break;
            case Jacobian::Line::ShapeD0X3:
              d *= band.lines[iline].lineshape[pos].D0().X3;
              band_copy.lines[iline].lineshape[pos].D0().X3 += d;
              break;
            case Jacobian::Line::ShapeG2X0:
              d *= band.lines[iline].lineshape[pos].G2().X0;
              band_copy.lines[iline].lineshape[pos].G2().X0 += d;
              break;
            case Jacobian::Line::ShapeG2X1:
              d *= band.lines[iline].lineshape[pos].G2().X1;
              band_copy.lines[iline].lineshape[pos].G2().X1 += d;
              break;
            case Jacobian::Line::ShapeG2X2:
              d *= band.lines[iline].lineshape[pos].G2().X2;
              band_copy.lines[iline].lineshape[pos].G2().X2 += d;
              break;
            case Jacobian::Line::ShapeG2X3:
              d *= band.lines[iline].lineshape[pos].G2().X3;
              band_copy.lines[iline].lineshape[pos].G2().X3 += d;
              break;
            case Jacobian::Line::ShapeD2X0:
              d *= band.lines[iline].lineshape[pos].D2().X0;
              band_copy.lines[iline].lineshape[pos].D2().X0 += d;
              break;
            case Jacobian::Line::ShapeD2X1:
              d *= band.lines[iline].lineshape[pos].D2().X1;
              band_copy.lines[iline].lineshape[pos].D2().X1 += d;
              break;
            case Jacobian::Line::ShapeD2X2:
              d *= band.lines[iline].lineshape[pos].D2().X2;
              band_copy.lines[iline].lineshape[pos].D2().X2 += d;
              break;
            case Jacobian::Line::ShapeD2X3:
              d *= band.lines[iline].lineshape[pos].D2().X3;
              band_copy.lines[iline].lineshape[pos].D2().X3 += d;
              break;
            case Jacobian::Line::ShapeFVCX0:
              d *= band.lines[iline].lineshape[pos].FVC().X0;
              band_copy.lines[iline].lineshape[pos].FVC().X0 += d;
              break;
            case Jacobian::Line::ShapeFVCX1:
              d *= band.lines[iline].lineshape[pos].FVC().X1;
              band_copy.lines[iline].lineshape[pos].FVC().X1 += d;
              break;
            case Jacobian::Line::ShapeFVCX2:
              d *= band.lines[iline].lineshape[pos].FVC().X2;
              band_copy.lines[iline].lineshape[pos].FVC().X2 += d;
              break;
            case Jacobian::Line::ShapeFVCX3:
              d *= band.lines[iline].lineshape[pos].FVC().X3;
              band_copy.lines[iline].lineshape[pos].FVC().X3 += d;
              break;
            case Jacobian::Line::ShapeETAX0:
              d *= band.lines[iline].lineshape[pos].ETA().X0;
              band_copy.lines[iline].lineshape[pos].ETA().X0 += d;
              break;
            case Jacobian::Line::ShapeETAX1:
              d *= band.lines[iline].lineshape[pos].ETA().X1;
              band_copy.lines[iline].lineshape[pos].ETA().X1 += d;
              break;
            case Jacobian::Line::ShapeETAX2:
              d *= band.lines[iline].lineshape[pos].ETA().X2;
              band_copy.lines[iline].lineshape[pos].ETA().X2 += d;
              break;
            case Jacobian::Line::ShapeETAX3:
              d *= band.lines[iline].lineshape[pos].ETA().X3;
              band_copy.lines[iline].lineshape[pos].ETA().X3 += d;
              break;
            case Jacobian::Line::Center:
              d *= band.lines[iline].F0;
              band_copy.lines[iline].F0 += d;
              break;
            case Jacobian::Line::Strength:
              d *= band.lines[iline].I0;
              band_copy.lines[iline].I0 += d;
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
            case Jacobian::Line::ECS_SCALINGX0:
            case Jacobian::Line::ECS_SCALINGX1:
            case Jacobian::Line::ECS_SCALINGX2:
            case Jacobian::Line::ECS_SCALINGX3:
            case Jacobian::Line::ECS_BETAX0:
            case Jacobian::Line::ECS_BETAX1:
            case Jacobian::Line::ECS_BETAX2:
            case Jacobian::Line::ECS_BETAX3:
            case Jacobian::Line::ECS_LAMBDAX0:
            case Jacobian::Line::ECS_LAMBDAX1:
            case Jacobian::Line::ECS_LAMBDAX2:
            case Jacobian::Line::ECS_LAMBDAX3:
            case Jacobian::Line::ECS_DCX0:
            case Jacobian::Line::ECS_DCX1:
            case Jacobian::Line::ECS_DCX2:
            case Jacobian::Line::ECS_DCX3:
            case Jacobian::Line::FINAL: {
              /* do nothing */
            }
          }
          
          // Perform calculations and estimate derivative
          const auto [dabs, dwork] = ecs_absorption_impl(T, H, P, this_vmr, vmrs, ecs_data, f_grid, zeeman_polarization, band_copy);
          vec -= dabs;
          vec /= -d;
          work &= dwork;
        } else if (Quantum::Number::StateMatch(target.qid, band.quantumidentity) == Quantum::Number::StateMatchType::Full) {
          ErrorCorrectedSuddenData ecs_data_copy = ecs_data;
          
          const auto spec = target.qid.Species();
          ARTS_USER_ERROR_IF(const Index pos = ecs_data.pos(spec); pos == ecs_data.size(), "No data for species ", spec, " in ecs_data:\n", ecs_data)
          switch (target.line) {
            case Jacobian::Line::ECS_SCALINGX0:
              d *= ecs_data_copy[spec].scaling.X0;
              ecs_data_copy[spec].scaling.X0 += d;
              break;
            case Jacobian::Line::ECS_SCALINGX1:
              d *= ecs_data_copy[spec].scaling.X1;
              ecs_data_copy[spec].scaling.X1 += d;
              break;
            case Jacobian::Line::ECS_SCALINGX2:
              d *= ecs_data_copy[spec].scaling.X2;
              ecs_data_copy[spec].scaling.X2 += d;
              break;
            case Jacobian::Line::ECS_SCALINGX3:
              d *= ecs_data_copy[spec].scaling.X3;
              ecs_data_copy[spec].scaling.X3 += d;
              break;
            case Jacobian::Line::ECS_BETAX0:
              d *= ecs_data_copy[spec].beta.X0;
              ecs_data_copy[spec].beta.X0 += d;
              break;
            case Jacobian::Line::ECS_BETAX1:
              d *= ecs_data_copy[spec].beta.X1;
              ecs_data_copy[spec].beta.X1 += d;
              break;
            case Jacobian::Line::ECS_BETAX2:
              d *= ecs_data_copy[spec].beta.X2;
              ecs_data_copy[spec].beta.X2 += d;
              break;
            case Jacobian::Line::ECS_BETAX3:
              d *= ecs_data_copy[spec].beta.X3;
              ecs_data_copy[spec].beta.X3 += d;
              break;
            case Jacobian::Line::ECS_LAMBDAX0:
              d *= ecs_data_copy[spec].lambda.X0;
              ecs_data_copy[spec].lambda.X0 += d;
              break;
            case Jacobian::Line::ECS_LAMBDAX1:
              d *= ecs_data_copy[spec].lambda.X1;
              ecs_data_copy[spec].lambda.X1 += d;
              break;
            case Jacobian::Line::ECS_LAMBDAX2:
              d *= ecs_data_copy[spec].lambda.X1;
              ecs_data_copy[spec].lambda.X1 += d;
              break;
            case Jacobian::Line::ECS_LAMBDAX3:
              d *= ecs_data_copy[spec].lambda.X0;
              ecs_data_copy[spec].lambda.X0 += d;
              break;
            case Jacobian::Line::ECS_DCX0:
              d *= ecs_data_copy[spec].collisional_distance.X0;
              ecs_data_copy[spec].collisional_distance.X0 += d;
              break;
            case Jacobian::Line::ECS_DCX1:
              d *= ecs_data_copy[spec].collisional_distance.X1;
              ecs_data_copy[spec].collisional_distance.X1 += d;
              break;
            case Jacobian::Line::ECS_DCX2:
              d *= ecs_data_copy[spec].collisional_distance.X2;
              ecs_data_copy[spec].collisional_distance.X2 += d;
              break;
            case Jacobian::Line::ECS_DCX3:
              d *= ecs_data_copy[spec].collisional_distance.X3;
              ecs_data_copy[spec].collisional_distance.X3 += d;
              break;
            case Jacobian::Line::ShapeG0X0:
            case Jacobian::Line::ShapeG0X1:
            case Jacobian::Line::ShapeG0X2:
            case Jacobian::Line::ShapeG0X3:
            case Jacobian::Line::ShapeD0X0:
            case Jacobian::Line::ShapeD0X1:
            case Jacobian::Line::ShapeD0X2:
            case Jacobian::Line::ShapeD0X3:
            case Jacobian::Line::ShapeG2X0:
            case Jacobian::Line::ShapeG2X1:
            case Jacobian::Line::ShapeG2X2:
            case Jacobian::Line::ShapeG2X3:
            case Jacobian::Line::ShapeD2X0:
            case Jacobian::Line::ShapeD2X1:
            case Jacobian::Line::ShapeD2X2:
            case Jacobian::Line::ShapeD2X3:
            case Jacobian::Line::ShapeFVCX0:
            case Jacobian::Line::ShapeFVCX1:
            case Jacobian::Line::ShapeFVCX2:
            case Jacobian::Line::ShapeFVCX3:
            case Jacobian::Line::ShapeETAX0:
            case Jacobian::Line::ShapeETAX1:
            case Jacobian::Line::ShapeETAX2:
            case Jacobian::Line::ShapeETAX3:
            case Jacobian::Line::Center:
            case Jacobian::Line::Strength:
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
            case Jacobian::Line::FINAL: {
              /* do nothing */
            }
          }
          
          // Perform calculations and estimate derivative
          const auto [dabs, dwork] = ecs_absorption_impl(T, H, P, this_vmr, vmrs, ecs_data_copy, f_grid, zeeman_polarization, band);
          vec -= dabs;
          vec /= -d;
          work &= dwork;
          break;  // Leave early because it is a band derivative
        }
      }
    } else {
      vec *= 0;  // No derivative, so don't mess around and remove everything
    }
  }
  
  return {std::move(absorption), std::move(jacobian), not work};
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
      band.lines[j].lineshape[i].Y() = LineShape::ModelParameters(LineShape::TemperatureModel::POLY, c);
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
        band.lines[j].lineshape[i].DV() = LineShape::ModelParameters(LineShape::TemperatureModel::POLY, c);
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
        band.lines[j].lineshape[i].G() = LineShape::ModelParameters(LineShape::TemperatureModel::POLY, c);
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
        //         band.lines[j].lineshape[i].DG() = LineShape::ModelParameters(LineShape::TemperatureModel::POLY, c);
      }
    }
  }
  
  // If we reach here, we have to set the band population type
  // to LTE and renormalize the strength to physical frequency
  band.normalization = Absorption::NormalizationType::SFS;
  band.population = Absorption::PopulationType::LTE;
  
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
      Y[k] += 2 * std::abs(dip[j] / dip[k]) * W(j, k) / (band.lines[k].F0 - band.lines[j].F0);
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
      G[k] += W(k, j) * W(j, k) / Math::pow2(band.lines[j].F0 - band.lines[k].F0);
      G[k] += Math::pow2(std::abs(dip[j] / dip[k]) * W(j, k) / (band.lines[j].F0 - band.lines[k].F0));
      G[k] += 2 * std::abs(dip[j] / dip[k]) * W(j, k) * W(k, k) / Math::pow2(band.lines[j].F0 - band.lines[k].F0);
      for (Index l=0; l<N; l++) {
        if (l == k or l == j) continue;
        G[k] -= 2 * std::abs(dip[j] / dip[k]) * W(j, l) * W(l, k) / ((band.lines[j].F0 - band.lines[k].F0) * (band.lines[l].F0 - band.lines[k].F0));
       
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
      DV[k] += W(k, j) * W(j, k) / (band.lines[j].F0 - band.lines[k].F0);
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
    f0[i] = band.lines[i].F0 - frenorm + P * band.lines[i].lineshape[broadener].D0().at(T, band.T0);
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
    eig.val[i] -= Complex(f0[i], P * band.lines[i].lineshape[broadener].G0().at(T, band.T0));
  }
  
  // The strength term should now be (1 + i y(T) + g(T)) * d**2 * rho(T)
  // d**2 * rho(T) * F0 * (1 - exp(-hF0/kT)) should be I0(T)
  // We only want to keep (i y(T) + g(T)).
  // So we divide by d**2 * rho(T) through the use of I0(T) and remove 1 from the real component
  for (Index i=0; i<pop.size(); i++) {
    const Numeric i0 = LineShape::LocalThermodynamicEquilibrium(band.lines[i].I0, band.T0, T, band.lines[i].F0, band.lines[i].E0, QT, QT0, 0, 1, 0, 0).S;
    eig.str[i] *= - band.lines[i].F0 * std::expm1(- (Constant::h * band.lines[i].F0) / (Constant::k * T)) / i0;
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
  const ArrayOfIndex sorting = sorted_population_and_dipole(band.T0, band).first;
  const Numeric QT0 = single_partition_function(band.T0, band.Isotopologue());
  
  // Output
  Tensor4 out(4, N, M, K);
  
  #pragma omp parallel for collapse(2) if (!arts_omp_in_parallel())
  for (Index m=0; m<M; m++) {
    for (Index k=0; k<K; k++) {
      const Numeric T = temperatures[k];
      const Numeric QT = single_partition_function(T, band.Isotopologue());
      
      // Relaxation matrix of T0 sorting at T
      ComplexMatrix W = single_species_ecs_relaxation_matrix(band, sorting, T, P, ecs_data[band.broadeningspecies[m]], m);
      for (Index n=0; n<N; n++) {
        W(n, n) += band.lines[sorting[n]].F0 - frenorm;
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
Tensor4 rosenkranz_approximation(const AbsorptionLines& band,
                                 const Vector& temperatures,
                                 const ErrorCorrectedSuddenData& ecs_data,
                                 const Numeric P) {
  const Index N = band.NumLines();
  const Index M = band.NumBroadeners();
  const Index K = temperatures.nelem();
  
  // Weighted center of the band
  const Numeric frenorm = band.F_mean();
  
  // Need sorting to put weak lines last, but we need the sorting constant or the output jumps
  const ArrayOfIndex sorting = sorted_population_and_dipole(band.T0, band).first;
  
  // Output
  Tensor4 out(4, N, M, K);
  
  #pragma omp parallel for collapse(2) if (!arts_omp_in_parallel())
  for (Index m=0; m<M; m++) {
    for (Index k=0; k<K; k++) {
      const Numeric T = temperatures[k];
      
      // Relaxation matrix of T0 sorting at T
      ComplexMatrix W = single_species_ecs_relaxation_matrix(band, sorting, T, P, ecs_data[band.broadeningspecies[m]], m);
      for (Index n=0; n<N; n++) {
        W(n, n) += band.lines[n].F0 - frenorm;
      }
      
      // Populations and dipoles of T0 sorting at T
      const auto [pop, dip] = presorted_population_and_dipole(T, sorting, band);
      
      out(0, joker, m, k) = RosenkranzG(dip, W.imag(), band);
      out(1, joker, m, k) = RosenkranzY(dip, W.imag(), band);
      out(2, joker, m, k) = RosenkranzDV(dip, W.imag(), band);
      out(3, joker, m, k) = 0;
    }
  }
  
  return out;
}

void ecs_eigenvalue_adaptation(AbsorptionLines& band,
                               const Vector& temperatures,
                               const ErrorCorrectedSuddenData& ecs_data,
                               const Numeric P0,
                               const Index ord,
                               const bool robust,
                               const bool rosenkranz_adaptation,
                               const Verbosity& verbosity) {
  CREATE_OUT3;
  ARTS_USER_ERROR_IF(P0 <= 0, P0, " Pa is not possible")

  ARTS_USER_ERROR_IF(not is_sorted(temperatures),
                     "The temperature list [",
                     temperatures,
                     "] K\n"
                     "must be fully sorted from low to high")

  ARTS_USER_ERROR_IF(
      ord < 1 or ord > 3, "Order not in list [1, 2, 3], is: ", ord)

  if (band_eigenvalue_adaptation(
          band,
          rosenkranz_adaptation not_eq 0
              ? rosenkranz_approximation(band, temperatures, ecs_data, P0)
              : ecs_eigenvalue_approximation(band, temperatures, ecs_data, P0),
          temperatures,
          P0,
          ord)) {
    ARTS_USER_ERROR_IF(not robust,
                       "Bad eigenvalue adaptation for band: ",
                       band.quantumidentity,
                       '\n')
    out3 << "Bad eigenvalue adaptation for band: " << band.quantumidentity
         << '\n';

    band.normalization = Absorption::NormalizationType::SFS;
    band.population = Absorption::PopulationType::LTE;
    for (auto& line : band.lines) {
      for (auto& sm : line.lineshape.Data()) {
        sm.Y() = LineShape::ModelParameters(LineShape::TemperatureModel::None);
        sm.G() = LineShape::ModelParameters(LineShape::TemperatureModel::None);
        sm.DV() = LineShape::ModelParameters(LineShape::TemperatureModel::None);
      }
    }
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

std::ostream& operator<<(std::ostream& os,
                         const ErrorCorrectedSuddenData& rbd) {
  for (Index i = 0; i < rbd.data.nelem(); i++) {
    if (i) os << '\n';
    os << rbd.data[i];
  }
  return os;
}

std::istream& operator>>(std::istream& is, ErrorCorrectedSuddenData& rbd) {
  for (auto& x : rbd.data) is >> x;
  return is;
}

Index ErrorCorrectedSuddenData::pos(Species::Species spec) const {
  return std::distance(data.begin(),
                       std::find(data.cbegin(), data.cend(), spec));
}

const SpeciesErrorCorrectedSuddenData& ErrorCorrectedSuddenData::operator[](
    Species::Species spec) const {
  if (auto ptr = std::find(data.cbegin(), data.cend(), spec);
      ptr not_eq data.cend())
    return *ptr;
  return *std::find(data.cbegin(), data.cend(), Species::Species::Bath);
}

SpeciesErrorCorrectedSuddenData& ErrorCorrectedSuddenData::operator[](
    Species::Species spec) {
  if (auto ptr = std::find(data.begin(), data.end(), spec);
      ptr not_eq data.end())
    return *ptr;
  return data.emplace_back(spec);
}

ErrorCorrectedSuddenData& MapOfErrorCorrectedSuddenData::operator[](
    const QuantumIdentifier& id) {
  if (auto ptr = std::find(begin(), end(), id); ptr not_eq end()) return *ptr;
  return emplace_back(id);
}

const ErrorCorrectedSuddenData& MapOfErrorCorrectedSuddenData::operator[](
    const QuantumIdentifier& id) const {
  if (auto ptr = std::find(cbegin(), cend(), id); ptr not_eq cend()) {
    return *ptr;
  }
  ARTS_USER_ERROR("Cannot find data for QuantumIdentifier:\n",
                  id,
                  '\n',
                  "Available data:\n",
                  *this);
  return front();  // To get rid of potential warnings...
}

std::ostream& operator<<(std::ostream& os,
                         const MapOfErrorCorrectedSuddenData& m) {
  std::for_each(m.cbegin(), m.cend(), [&](auto& x) { os << x << '\n'; });
  return os;
}
}  // namespace Absorption::LineMixing
