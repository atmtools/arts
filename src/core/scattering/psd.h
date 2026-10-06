#ifndef ARTS_CORE_SCATTERING_PSD_H_
#define ARTS_CORE_SCATTERING_PSD_H_

#include <matpack.h>

#include <optional>
#include <unordered_map>

#include "atm.h"
#include "enums.h"
#include "properties.h"
#include "utils.h"

namespace scattering {

struct PSDData {
  Vector                                                values;
  std::unordered_map<ScatteringSpeciesProperty, Vector> derivatives;
};

/** A single-particle population scaled by an atmospheric number-density field. */
struct MonodispersePSD {
  ScatteringSpeciesProperty number_density;
  Numeric                   t_min = 0.0;
  Numeric                   t_max = 350.0;

  MonodispersePSD() = default;
  MonodispersePSD(ScatteringSpeciesProperty number_density_, Numeric t_min_ = 0.0, Numeric t_max_ = 350.0);

  static constexpr SizeParameter get_size_parameter() { return SizeParameter::DVeq; }
  //! evaluate gives the number density [m^-3] of the one particle
  static constexpr bool size_density = false;
  Vector                         evaluate(const AtmPoint&, const Vector&, Numeric, Numeric) const;
  PSDData                        evaluate_with_derivatives(const AtmPoint&, const Vector&, Numeric, Numeric) const;
};

/*** Single-moment modified gamma distribution
 *
 * Implements a modified gamma distribution with a single free moment.
 * Currently this moment is expected to be the mass density.
 */
struct MGDSingleMoment {
  ScatteringSpeciesProperty moment;
  Numeric                   n_alpha;
  Numeric                   n_b;
  Numeric                   mu;
  Numeric                   gamma;
  Numeric                   t_min;
  Numeric                   t_max;
  bool                      picky;

  MGDSingleMoment() = default;

  MGDSingleMoment(ScatteringSpeciesProperty moment_,
                  Numeric                   n_alpha_,
                  Numeric                   n_b_,
                  Numeric                   mu_,
                  Numeric                   gamma_,
                  Numeric                   t_min_,
                  Numeric                   t_max_,
                  bool                      picky_);

  MGDSingleMoment(ScatteringSpeciesProperty moment_, std::string name, Numeric t_min_, Numeric t_max_, bool picky_);

  static constexpr SizeParameter get_size_parameter() { return SizeParameter::DVeq; }
  //! evaluate gives the density per unit size [m^-3 / unit of the size parameter] at each size
  static constexpr bool size_density = true;

  /** Evaluate PSD at given atmospheric point.
   *
   * @param point The atmospheric point.
   * @param particle_size A vector containing the particles sizes at which to
   * evaluate the PSD.
   * @param scat_species_a The a parameter of the mass-size relationship of
   * the particle data.
   * @param scat_species_b The b parameter of the mass-size relationship of
   * the particle data.
   *
   */
  Vector evaluate(const AtmPoint& point,
                  const Vector&   particle_sizes,
                  const Numeric&  scat_species_a,
                  const Numeric&  scat_species_b) const;

  PSDData evaluate_with_derivatives(const AtmPoint&, const Vector&, Numeric, Numeric) const;
};

/** Modified-gamma distribution constrained by particle mass density.
 *
 * Exactly one of n0 and lambda must be NaN.  The missing parameter is derived
 * from the mass moment and the supplied mass-size relationship m=a*x^b.
 */
struct MGDMass {
  ScatteringSpeciesProperty mass;
  Numeric                   n0;
  Numeric                   mu;
  Numeric                   lambda;
  Numeric                   gamma;
  Numeric                   t_min;
  Numeric                   t_max;
  bool                      picky;

  MGDMass() = default;
  MGDMass(ScatteringSpeciesProperty, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, bool);

  static constexpr SizeParameter get_size_parameter() { return SizeParameter::DVeq; }
  //! evaluate gives the density per unit size [m^-3 / unit of the size parameter] at each size
  static constexpr bool size_density = true;
  Vector                         evaluate(const AtmPoint&, const Vector&, Numeric, Numeric) const;
  PSDData                        evaluate_with_derivatives(const AtmPoint&, const Vector&, Numeric, Numeric) const;
};

enum class MGDTwoMomentType : char { MassMeanSize, MassNumberDensity };

/** Modified-gamma distribution constrained by mass and a second moment. */
struct MGDTwoMoment {
  ScatteringSpeciesProperty mass;
  ScatteringSpeciesProperty second_moment;
  MGDTwoMomentType          type;
  Numeric                   mu;
  Numeric                   gamma;
  Numeric                   t_min;
  Numeric                   t_max;
  bool                      picky;

  MGDTwoMoment() = default;
  MGDTwoMoment(
      ScatteringSpeciesProperty, ScatteringSpeciesProperty, MGDTwoMomentType, Numeric, Numeric, Numeric, Numeric, bool);

  static constexpr SizeParameter get_size_parameter() { return SizeParameter::DVeq; }
  //! evaluate gives the density per unit size [m^-3 / unit of the size parameter] at each size
  static constexpr bool size_density = true;
  Vector                         evaluate(const AtmPoint&, const Vector&, Numeric, Numeric) const;
  PSDData                        evaluate_with_derivatives(const AtmPoint&, const Vector&, Numeric, Numeric) const;
};

/** Delanoe et al. (2014) normalized ice PSD. */
struct DelanoeEtAl14 {
  ScatteringSpeciesProperty mass;
  ScatteringSpeciesProperty intercept_parameter;
  ScatteringSpeciesProperty mean_size;
  Numeric                   rho;
  Numeric                   alpha;
  Numeric                   beta;
  Numeric                   t_min;
  Numeric                   t_max;
  Numeric                   dm_min;
  bool                      picky;

  DelanoeEtAl14() = default;
  DelanoeEtAl14(ScatteringSpeciesProperty, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, bool);
  DelanoeEtAl14(ScatteringSpeciesProperty,
                ScatteringSpeciesProperty,
                ScatteringSpeciesProperty,
                Numeric,
                Numeric,
                Numeric,
                Numeric,
                Numeric,
                Numeric,
                bool);

  static constexpr SizeParameter get_size_parameter() { return SizeParameter::DVeq; }
  //! evaluate gives the density per unit size [m^-3 / unit of the size parameter] at each size
  static constexpr bool size_density = true;
  Vector                         evaluate(const AtmPoint&, const Vector&, Numeric, Numeric) const;
  PSDData                        evaluate_with_derivatives(const AtmPoint&, const Vector&, Numeric, Numeric) const;
};

/** Field et al. (2007) tropical or midlatitude ice PSD. */
struct FieldEtAl07 {
  ScatteringSpeciesProperty mass;
  bool                      tropical;
  Numeric                   t_min;
  Numeric                   t_max;
  Numeric                   t_min_psd;
  Numeric                   t_max_psd;
  bool                      picky;

  FieldEtAl07() = default;
  FieldEtAl07(ScatteringSpeciesProperty, std::string, Numeric, Numeric, Numeric, Numeric, bool);

  static constexpr SizeParameter get_size_parameter() { return SizeParameter::DVeq; }
  //! evaluate gives the density per unit size [m^-3 / unit of the size parameter] at each size
  static constexpr bool size_density = true;
  Vector                         evaluate(const AtmPoint&, const Vector&, Numeric, Numeric) const;
  PSDData                        evaluate_with_derivatives(const AtmPoint&, const Vector&, Numeric, Numeric) const;
};

/** McFarquhar and Heymsfield (1997) cloud-ice PSD. */
struct McFarquharHeymsfield97 {
  ScatteringSpeciesProperty mass;
  Numeric                   t_min;
  Numeric                   t_max;
  Numeric                   t_min_psd;
  Numeric                   t_max_psd;
  bool                      picky;

  McFarquharHeymsfield97() = default;
  McFarquharHeymsfield97(ScatteringSpeciesProperty, Numeric, Numeric, Numeric, Numeric, bool);

  static constexpr SizeParameter get_size_parameter() { return SizeParameter::DVeq; }
  //! evaluate gives the density per unit size [m^-3 / unit of the size parameter] at each size
  static constexpr bool size_density = true;
  Vector                         evaluate(const AtmPoint&, const Vector&, Numeric, Numeric) const;
  PSDData                        evaluate_with_derivatives(const AtmPoint&, const Vector&, Numeric, Numeric) const;
};

/*** Binned PSD
 *
 * The BinnedPSD class represents a particle size distribution using a fixed number of particle counts
 * over a sequence of size bins. Particles with sizes outside the size and temperature range are set
 * to zero.
 */
struct BinnedPSD {
  SizeParameter size_parameter = SizeParameter::DVeq;
  Vector        bins;
  Vector        counts;
  Numeric       t_min = 0.0;
  Numeric       t_max = 350.0;

  BinnedPSD() = default;

  BinnedPSD(SizeParameter size_parameter_, Vector bins_, Vector counts_, Numeric t_min_ = 0.0, Numeric t_max_ = 350.0);

  SizeParameter get_size_parameter() const { return size_parameter; }
  //! evaluate gives the count of the bin each size falls in, a number density [m^-3] per particle
  static constexpr bool size_density = false;

  Vector evaluate(const AtmPoint& point,
                  const Vector&   particle_sizes,
                  const Numeric& /*scat_species_a*/,
                  const Numeric& /*scat_species_b*/) const;

  PSDData evaluate_with_derivatives(const AtmPoint&, const Vector&, Numeric, Numeric) const;
};

/** The widths that turn a density per unit size into the number density each size stands for
 *
 * Half the distance between a size's two neighbours, and the distance to the
 * one neighbour at either end.  The sizes must ascend strictly, and there
 * must be at least two.
 */
Vector size_bin_widths(const Vector& sizes);

/** The number density [m^-3] of the particles that each size stands for
 *
 * The PSD's values times size_bin_widths if it gives a density per unit
 * size (PSD::size_density), else its values.
 */
template <typename PSD>
Vector number_densities(const PSD& psd, const AtmPoint& point, const Vector& sizes, Numeric a, Numeric b) {
  Vector out = psd.evaluate(point, sizes, a, b);
  if constexpr (PSD::size_density) {
    const Vector w = size_bin_widths(sizes);
    for (Size i = 0; i < out.size(); i++) out[i] *= w[i];
  }
  return out;
}

//! number_densities and their derivatives
template <typename PSD>
PSDData number_densities_with_derivatives(
    const PSD& psd, const AtmPoint& point, const Vector& sizes, Numeric a, Numeric b) {
  PSDData out = psd.evaluate_with_derivatives(point, sizes, a, b);
  if constexpr (PSD::size_density) {
    const Vector w = size_bin_widths(sizes);
    for (Size i = 0; i < out.values.size(); i++) out.values[i] *= w[i];
    for (auto& [key, d] : out.derivatives)
      for (Size i = 0; i < d.size(); i++) d[i] *= w[i];
  }
  return out;
}
}  // namespace scattering
#endif  // ARTS_CORE_PSD_H_
