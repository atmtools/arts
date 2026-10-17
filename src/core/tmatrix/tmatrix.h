#pragma once

#include <enumsEquivalentRadius.h>
#include <enumsParticleShape.h>
#include <enumsTMatrixSizeDistribution.h>
#include <rtepack.h>

#include <span>
#include <vector>

namespace tmatrix {
//! Whether the optional T-matrix solver is built (ENABLE_TMATRIX).
bool available();

/** The radius is that of the equivalent sphere of the same volume or the same surface area.
 *
 * Shapes: a spheroid with aspect_ratio its horizontal-to-rotational axis ratio, a cylinder
 * with aspect_ratio its diameter-to-length ratio, or a Chebyshev particle
 * r_0 (1 + aspect_ratio cos(degree theta)) with |aspect_ratio| < 1 and degree >= 1.  The
 * degree is 0 for the other shapes.
 */

//! A fixed-orientation geometry, in degrees: the zenith and azimuth angles of the incident and the scattered
//! direction in the laboratory frame, and the Euler angles alpha and beta of the particle's symmetry axis.
struct Geometry {
  Numeric theta_incident{}, theta_scattered{}, phi_incident{}, phi_scattered{}, alpha{}, beta{};
};

//! The largest distance between two points of the surface, for the volume-equivalent radius.
Numeric maximum_dimension(
    Numeric radius, Numeric aspect_ratio, ParticleShape shape, Index degree = 0, const Vector& coefficients = {});

//! Owning results.
struct FixedResult {
  Index   order{};
  Numeric scattering{}, extinction{};
  //! The Jones matrix [[S_vv, S_vh], [S_hv, S_hh]] between ARTS's meridional bases of the two rays
  Jones amplitude{Complex{}};
  //! Its Mueller matrix, rtepack::to_muelmat(amplitude)
  Muelmat phase{0.0};
};

//! Angles in degrees. Radius and wavelength must use the same length unit.
//! Amplitude has that length unit; phase and cross sections its square.
FixedResult fixed(Numeric          radius,
                  Numeric          wavelength,
                  Numeric          aspect_ratio,
                  Numeric          refractive_real,
                  Numeric          refractive_imag,
                  Numeric          theta_incident,
                  Numeric          theta_scattered,
                  Numeric          phi_incident,
                  Numeric          phi_scattered,
                  Numeric          alpha,
                  Numeric          beta,
                  Numeric          accuracy           = 0.001,
                  EquivalentRadius equivalent_radius  = EquivalentRadius::volume,
                  ParticleShape    shape              = ParticleShape::spheroid,
                  Index            degree             = 0,
                  const Vector&    coefficients       = {},
                  Index            surface_quadrature = 2);

//! Compute once, then evaluate each geometry.
//! Caller supplies one output per geometry; results own their data.
void fixed_batch(std::span<FixedResult>    results,
                 Numeric                   radius,
                 Numeric                   wavelength,
                 Numeric                   aspect_ratio,
                 Numeric                   refractive_real,
                 Numeric                   refractive_imag,
                 std::span<const Geometry> geometries,
                 Numeric                   accuracy           = 0.001,
                 EquivalentRadius          equivalent_radius  = EquivalentRadius::volume,
                 ParticleShape             shape              = ParticleShape::spheroid,
                 Index                     degree             = 0,
                 const Vector&             coefficients       = {},
                 Index                     surface_quadrature = 2);

//! Allocating convenience overload; delegates to the caller-owned span overload.
std::vector<FixedResult> fixed_batch(Numeric                   radius,
                                     Numeric                   wavelength,
                                     Numeric                   aspect_ratio,
                                     Numeric                   refractive_real,
                                     Numeric                   refractive_imag,
                                     std::span<const Geometry> geometries,
                                     Numeric                   accuracy           = 0.001,
                                     EquivalentRadius          equivalent_radius  = EquivalentRadius::volume,
                                     ParticleShape             shape              = ParticleShape::spheroid,
                                     Index                     degree             = 0,
                                     const Vector&             coefficients       = {},
                                     Index                     surface_quadrature = 2);

struct RandomResult {
  Numeric effective_radius{}, effective_variance{}, extinction{}, scattering{}, albedo{}, asymmetry{};
  //! The scattering matrix at the equally spaced scattering angles 0..180 degrees,
  //! normalized as in Mishchenko (integral F11 = 4*pi), in ARTS's Stokes
  //! convention: F34 is the negative of Mishchenko's.
  CompactPlanarMuelmatVector phase;
  //! The same scattering matrix as its exact Legendre series, the coefficients
  //! on Y_l0 of scattering/tro_legendre.h for l = 0, ..., l_max.
  CompactPlanarMuelmatVector legendre;
};

//! One size distribution (NPNAX=1). Distribution parameters follow TMD.
//! Default is an effectively monodisperse, equal-volume-radius particle.
//! angles is the number of equally spaced scattering angles of phase, 0 for none.
RandomResult random(Numeric                 radius,
                    Numeric                 wavelength,
                    Numeric                 aspect_ratio,
                    Numeric                 refractive_real,
                    Numeric                 refractive_imag,
                    Index                   angles             = 0,
                    Numeric                 accuracy           = 0.001,
                    EquivalentRadius        equivalent_radius  = EquivalentRadius::volume,
                    ParticleShape           shape              = ParticleShape::spheroid,
                    Index                   degree             = 0,
                    const Vector&           coefficients       = {},
                    TMatrixSizeDistribution distribution       = TMatrixSizeDistribution::gamma,
                    Numeric                 b                  = 0.1,
                    Numeric                 gamma              = 1,
                    Index                   size_quadrature    = -1,
                    Index                   surface_quadrature = 2,
                    Numeric                 lower_radius_ratio = 0.9999999,
                    Numeric                 upper_radius_ratio = 1.0000001);
//! Owning results of azimuthally random orientation, on the grids of the call.
struct AzimuthallyRandomResult {
  Index order{};
  //! The phase matrix [za_inc, delta_aa, za_scat] between ARTS's meridional bases, in length squared
  MuelmatTensor3 phase;
  //! K11, K12 and K34 of the extinction matrix [za_inc], in length squared
  CompactPlanarPropmatVector extinction;
  //! I and Q of the absorption vector [za_inc], in length squared
  CompactPlanarStokvecVector absorption;
};

/** Particles of a tilt distribution with random azimuth (ARO), on ARTS's grids.
 *
 * za_inc and za_scat are the zenith angles of the propagation directions and delta_aa the scattered minus the
 * incident azimuth, in degrees, as ARTS's scattering data hold them.  tilt_angles are the angles beta between the
 * symmetry axis and the vertical [deg], in [0, 180], with tilt_weights, which need not be normalized.  The other
 * arguments are those of fixed.
 */
AzimuthallyRandomResult azimuthally_random(Numeric          radius,
                                           Numeric          wavelength,
                                           Numeric          aspect_ratio,
                                           Numeric          refractive_real,
                                           Numeric          refractive_imag,
                                           const Vector&    za_inc,
                                           const Vector&    delta_aa,
                                           const Vector&    za_scat,
                                           const Vector&    tilt_angles,
                                           const Vector&    tilt_weights,
                                           Numeric          accuracy           = 0.001,
                                           EquivalentRadius equivalent_radius  = EquivalentRadius::volume,
                                           ParticleShape    shape              = ParticleShape::spheroid,
                                           Index            degree             = 0,
                                           const Vector&    coefficients       = {},
                                           Index            surface_quadrature = 2);
}  // namespace tmatrix
