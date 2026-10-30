#pragma once

#include <atm.h>
#include <enums.h>
#include <jacobian.h>
#include <path_point.h>
#include <rtepack.h>
#include <scattering_species.h>

#include <functional>
#include <span>
#include <vector>

using ArrayOfRadarAuxiliary = Array<RadarAuxiliary>;

namespace radar {
/** The liquid-water dielectric factor |K|^2 of the Liebe-93 permittivity model.
 *
 * @param frequency The frequency [Hz]
 * @param temperature The temperature [K]
 * @return |K|^2
 */
[[nodiscard]] Numeric liebe93_k2(Numeric frequency, Numeric temperature);

/** The factor that converts backscatter [1/(m sr)] to Ze [mm^6/m^3].
 *
 * @param frequency The frequency [Hz]
 * @param ze_tref The temperature of the Liebe-93 dielectric factor [K]
 * @param k2 The dielectric factor |K|^2, negative selects Liebe-93
 * @return The conversion factor
 */
[[nodiscard]] Numeric ze_factor(Numeric frequency, Numeric ze_tref, Numeric k2);

/** The factor that converts backscatter to the linear part of the unit.
 *
 * dBZe uses the Ze factor, the logarithm is left to the caller.
 */
[[nodiscard]] Numeric unit_factor(RadarUnit unit, Numeric frequency, Numeric ze_tref, Numeric k2);

/** Resolve the Legacy mode from the largest range limit, other modes are returned as is */
[[nodiscard]] RadarRangeMode resolve(RadarRangeMode mode, Numeric largest_limit);

//! The ranges of the points of a radar path
struct PathRange {
  //! The geometric thickness of every layer [np - 1]
  Vector layer;

  //! The one-way distance from the radar [np]
  Vector distance;

  //! The round-trip time from the radar, including the path group refractive indices [np]
  Vector time;
};

/** The ranges of a radar path.
 *
 * The first point of the path need not be at the radar, for instance when
 * the radar is outside the atmosphere.  The range to it is then the geometric
 * distance in vacuum.
 *
 * @param path The propagation path, starting closest to the radar
 * @param radar_pos The position of the radar
 * @param ellipsoid The reference ellipsoid
 */
[[nodiscard]] PathRange path_range(const ArrayOfPropagationPathPoint& path,
                                   const Vector3&                     radar_pos,
                                   const Vector2&                     ellipsoid);

/** The coordinate of the range-gate limits at every path point.
 *
 * @param range The ranges of the path
 * @param path The propagation path
 * @param mode The coordinate, Legacy must already be resolved
 */
[[nodiscard]] Vector range_coordinate(const PathRange&                   range,
                                      const ArrayOfPropagationPathPoint& path,
                                      RadarRangeMode                     mode);

/** Append the quadrature of a range gate.
 *
 * Every layer that overlaps the gate gets one sample at the midpoint of the
 * overlap.  The weights are normalized by the full gate width, so the parts
 * of a gate without any path contribute nothing.
 *
 * @param samples The samples to extend
 * @param coordinate The range coordinate at every path point
 * @param lo The lower gate limit
 * @param hi The upper gate limit
 * @param gate The output index of the gate
 */
void add_gate_samples(
    std::vector<rtepack::radar_sample>& samples, const ConstVectorView& coordinate, Numeric lo, Numeric hi, Size gate);

/** Sort the samples by layer.
 *
 * @return The number of levels the samples need
 */
[[nodiscard]] Size sort_samples(std::vector<rtepack::radar_sample>& samples);

/** The sum of the sample weights of every gate, which is the covered fraction of the gate */
[[nodiscard]] Vector gate_coverage(const std::span<const rtepack::radar_sample>& samples, Size ngate);

/** The weighted sum of level values, linearly interpolated to the samples, of every gate */
[[nodiscard]] Vector gate_sum(const std::span<const rtepack::radar_sample>& samples,
                              const ConstVectorView&                        level,
                              Size                                          ngate);

/** The gas propagation matrices on the frequency grid, and their Jacobian [nq, nf], at a path point */
using GasModel = std::function<void(
    PropmatVector& propmat, PropmatMatrix& propmat_jac, const PropagationPathPoint& point, const AtmPoint& atm_point)>;

//! The optical properties of the levels of a radar path
struct Optics {
  //! Propagation matrices away from and towards the radar [np][nf]
  std::vector<PropmatVector> K_out, K_ret;

  //! Their derivatives [np][nq, nf]
  std::vector<PropmatMatrix> dK_out, dK_ret;

  //! The backscatter matrix [np][nf]
  std::vector<MuelmatVector> Z;

  //! Its derivative [np][nq, nf]
  std::vector<MuelmatMatrix> dZ;

  //! The extinction away from the radar by absorbing species and by (scaled) particles [nf, np]
  Matrix gas_extinction, particle_extinction;

  [[nodiscard]] Size size() const { return K_out.size(); }
};

/** The optical properties of the first np points of a radar path.
 *
 * The scattering properties are evaluated in the monostatic geometry of
 * every path point.  Derivatives are with regards to the atmospheric targets
 * of jac_targets, the gas model provides all others.
 *
 * @param path The propagation path, starting closest to the radar
 * @param np The number of levels to evaluate
 * @param atm_field The atmosphere
 * @param freq_grid The frequency grid
 * @param scattering_species The scattering species
 * @param jac_targets The Jacobian targets
 * @param pext_scaling The scaling of particle extinction
 * @param gas The gas propagation matrices
 */
[[nodiscard]] Optics optics(const ArrayOfPropagationPathPoint& path,
                            Size                               np,
                            const AtmField&                    atm_field,
                            const AscendingGrid&               freq_grid,
                            const ArrayOfScatteringSpecies&    scattering_species,
                            const JacobianTargets&             jac_targets,
                            Numeric                            pext_scaling,
                            const GasModel&                    gas);

/** The single-scattering return of every gate, and its derivatives.
 *
 * Wraps rtepack::radar_single_scattering for the levels of optics.
 *
 * @param[out] I The return of every gate [ngate, nf]
 * @param[out] dI The derivatives [ngate, np, nq, nf]
 * @param transmitted The transmitted Stokes vector
 * @param optics The optical properties of the levels
 * @param layer The layer thicknesses of the path, at least optics.size() - 1
 * @param samples The gate samples, sorted by layer
 * @param ngate The number of gates
 * @param option The transmittance option
 */
void single_scattering(StokvecMatrix&                                I,
                       StokvecTensor4&                               dI,
                       const Stokvec&                                transmitted,
                       const Optics&                                 optics,
                       const ConstVectorView&                        layer,
                       const std::span<const rtepack::radar_sample>& samples,
                       Size                                          ngate,
                       TransmittanceOption                           option);
}  // namespace radar
