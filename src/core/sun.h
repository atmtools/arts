/*!
  \file   sun.h
  \author Jon Petersen <jon.petersen@studium.uni-hamburg.de>
          Manfred Brath  <manfred.brath@uni-hamburg.de>
  \date   2021-02-22

  \brief  Declaration of functions in star.cc.
*/

#ifndef star_h
#define star_h

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <matpack.h>
#include <path_point.h>
#include <rtepack.h>
#include <surf.h>

/*===========================================================================
  === structs/classes  in sun.h
  ===========================================================================*/

/** The structure to describe a propagation path and releated quantities.
 *
 *  The fields of the structure are described more in detail inside the ARTS
 *  user guide (AUG).
 */
struct Sun {
  /** Sun description */
  String description;
  /** Sun spectrum, monochrmatic radiance spectrum at the surface of the sun*/
  Matrix spectrum;
  /** Sun radius */
  Numeric radius;
  /** distance from center of planet to center of sun*/
  Numeric distance;
  /** latitude of the sun in the sky of the planet */
  Numeric latitude;
  /** longitude of the sun in the sky of the planet */
  Numeric longitude;

  [[nodiscard]] Numeric sin_alpha_squared(Vector3 pos, Vector2 ell) const;

  friend std::ostream& operator<<(std::ostream& os, const Sun& sun);
};

/** An array of sun. */
using ArrayOfSun = Array<Sun>;

std::ostream& operator<<(std::ostream& os, const ArrayOfSun& a);

/** regrid_sun_spectrum
 *
 * Regrids a given spectrum from a griddedfield2 to the f_grid.
 * if the f_grid covers a larger range as the given one, one
 * can choose between two padding options:
 * zeros: Intensities outside the given spectrum are set to zero 
 * planck: Intensities outside the given spectrum are initilizied 
 *        with the black body value at that frequency.
 *
 * @param[in]  sun_spectrum_raw  gf2 of the given spectrum.
 * @param[in]  f_grid  f_grid for the calculation.
 * @param[in]  temperature  Temperature for the planck padding.
 *
 * @return     interpolated spectrum
 *
 * @author Jon Petersen
 * @date   2022-01-19
 */
Matrix regrid_sun_spectrum(const GriddedField2& sun_spectrum_raw,
                           const Vector& f_grid,
                           const Numeric& temperature);

Vector3 sph2cart(const Vector3 sph);

/** Checks and sets sun radiance if sun is in line of sight.
 *
 * @param[out] spectral_radiance Spectral radiance of sun if set.
 * @param[in] sun Sun-structure.
 * @param[in] propagation_path_point A path poistion and line of sight.
 * @param[in] surface_field The surface for the ellipsoid.
 * @return True if sun is in line of sight and spectral_radiance is set.
  */
bool set_spectral_radiance_if_sun_intersection(
    StokvecVector& spectral_radiance,
    const Sun& sun,
    const PropagationPathPoint& propagation_path_point,
    const SurfaceField& surface_field);

/** Checks if the sun is within the line of sight.
 * 
 * Returns the angle between the sun center and the line of sight
 * as first output, and whether the sun is hit as second output.
 *
 * @param sun 
 * @param pos [h, lat, lon]
 * @param los [za, aa]
 * @param ell [a, b]
 * @return Angle and whether or not the sun is hit
 */
std::pair<Numeric, bool> hit_sun(const Sun& sun,
                                 const Vector3 pos,
                                 const Vector2 los,
                                 const Vector2 ell);

template <>
struct std::formatter<Sun> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Sun& v, FmtContext& ctx) const {
    const std::string_view sep   = tags.sep();
    const std::string_view quote = tags.quote();

    tags.add_if_bracket(ctx, '[');
    tags.format(ctx,
                quote,
                v.description,
                quote,
                sep,
                v.spectrum,
                sep,
                v.radius,
                sep,
                v.distance,
                sep,
                v.latitude,
                sep,
                v.longitude);
    tags.add_if_bracket(ctx, ']');
    return ctx.out();
  }
};

#endif /* star_h */
