#pragma once

#include <sun.h>

class Agenda;
class Workspace;

/*===========================================================================
  === Functions in sun.h
  ===========================================================================*/

/** Calculates the radiance spectrum of sun which is scattered by the atmospheric gases.
 *
 * @param[in,out] ws ARTS workspace.
 * @param[out] scattered_sunlight StokvecVector scattered monochromatic radiance
 *             spectrum of sun.
 * @param[in] f_grid Vector frequency grid.
 * @param[in] p Numeric pressure at location of scattering.
 * @param[in] T Numeric temperature at location of scattering.
 * @param[in] vmr Vector volume mixing ratios of absorption species at location
 *            of scattering.
 * @param[in] transmitted_sunlight Matrix transmitted monochromatic irradiance
 *             spectrum of sun at location of scattering.
 * @param[in] gas_scattering_los_in Vector incoming direction of the transmitted sun irradiance
 *            spectrum.
 * @param[in] gas_scattering_los_out outgoing direction of the transmitted sun irradiance
 *            spectrum.
 * @param[in] gas_scattering_agenda Agenda agenda calculating the gas scattering
 *            cross sectionand matrix.
 */
void get_scattered_sunsource(const Workspace& ws,
                             StokvecVector& scattered_sunlight,
                             const Vector& f_grid,
                             const AtmPoint& atm_point,
                             const Matrix& transmitted_sunlight,
                             const Vector& gas_scattering_los_in,
                             const Vector& gas_scattering_los_out,
                             const Agenda& gas_scattering_agenda);

/** Finds a path from the observer to the sun.
 *
 * Computes the angular offset between the observer and the sun, and returns the
 * the path in output parameter.  The algorithm first checks the path to the sun
 * as if it was geometric.  It then proceeds to look up, down, left, and right,
 * using a multiple of the angular offset from the sun based on the space-facing
 * point in the ray path.
 *
 * This multiple starts a 1x the angular offset and is decreased by a factor of
 * 0.5 per level of refinement.  So a refinement of 2 would look at 1x, 0.5x, and
 * then return. Of 3 would look at 1x, 0.5x, 0.25x, and then return.
 *
 * Two other speed-up parameters are provided.  The first is the angle_cut, which
 * stops the calculations if the angular offset to the sun is larger than 
 * 
 * @param[in] ws ARTS workspace
 * @param[out] sun_path A path to the sun.
 * @param[in] sun A sun object.
 * @param[in] ray_path_observer_agenda As WSV.
 * @param[in] surface_field As WSV.
 * @param[in] observer_pos Position of the observer.
 * @param[in] angle_cut Angular cutoff to return the path, see above.
 * @param[in] refinements Refinements of the resolution, see above.
 * @param[in] just_hit If true, exits the moment a sun is hit.
 */
void find_sun_path(const Workspace& ws,
                   ArrayOfPropagationPathPoint& sun_path,
                   const Sun& sun,
                   const Agenda& ray_path_observer_agenda,
                   const SurfaceField& surface_field,
                   const Vector3 observer_pos,
                   const Numeric angle_cut,
                   const Index refinements,
                   const bool just_hit);

std::pair<Numeric, bool> beta_angle(const Workspace& ws,
                                    ArrayOfPropagationPathPoint& sun_path,
                                    const Sun& sun,
                                    const Vector3& observer_pos,
                                    const Vector2& observer_los,
                                    const Agenda& ray_path_observer_agenda,
                                    const SurfaceField& surface_field,
                                    const Numeric& angle_cut);

Vector2 geometric_los(const Vector3 from, const Vector3 to, const Vector2 ell);
