#pragma once

#include <atm.h>
#include <subsurf.h>

#include <concepts>
#include <variant>

/** The media that scatterers are evaluated in
 *
 * A scatterer has explicit methods that take the physical quantities it
 * depends on, and template wrappers that take a point of a medium and its
 * spectral point.  The wrappers get those quantities from the points with the
 * functions below, one overload per medium, and call the explicit methods.  A
 * medium that lacks a quantity throws when it is asked for it.  A user's
 * callback instead takes the points as a MediumPoint and a SpectralMediumPoint
 * and deals with the media itself.
 */
namespace scattering {
//! The point of a medium that a scatterer can be evaluated at
template <typename T>
concept medium_point = std::same_as<T, AtmPoint> or std::same_as<T, SubsurfacePoint>;

//! The spectral point of the medium of Point: AtmSpectralPoint for AtmPoint, SubsurfaceSpectralPoint for SubsurfacePoint
template <typename T, typename Point>
concept spectral_medium_point =
    (std::same_as<Point, AtmPoint> and std::same_as<T, AtmSpectralPoint>) or
    (std::same_as<Point, SubsurfacePoint> and std::same_as<T, SubsurfaceSpectralPoint>);

//! The key of a derivative target in the medium of Point
template <medium_point Point>
using medium_key = std::conditional_t<std::same_as<Point, AtmPoint>, AtmKeyVal, SubsurfaceKeyVal>;

//! The scattering species properties of a point, by property
using ScatteringSpeciesMap = AtmPoint::ScatteringSpeciesMap;

//! The point of any medium, as a callback takes it: whichever medium it is called in
using MediumPoint = std::variant<AtmPoint, SubsurfacePoint>;

//! The spectral point of any medium, as a callback takes it
using SpectralMediumPoint = std::variant<AtmSpectralPoint, SubsurfaceSpectralPoint>;

//! The values of a spectral property at every frequency of an atmospheric spectral point, which has nf of them
[[nodiscard]] const Vector& spectral_property(const AtmSpectralPoint& point, const String& name, Size nf);

//! The values of a spectral property at every frequency of a subsurface spectral point, which has nf of them
[[nodiscard]] const Vector& spectral_property(const SubsurfaceSpectralPoint& point, const String& name, Size nf);

//! The scattering species properties of an atmospheric point
[[nodiscard]] const ScatteringSpeciesMap& scattering_species_properties(const AtmPoint& point);

//! The subsurface holds no scattering species properties, so this throws
[[nodiscard]] const ScatteringSpeciesMap& scattering_species_properties(const SubsurfacePoint& point);

/** The derivatives of the scattering species properties of an atmospheric point by a target
 *
 * @return The property that is the target with derivative 1, or nothing if the target is no property
 */
[[nodiscard]] ScatteringSpeciesMap scattering_species_property_derivatives(const AtmKeyVal& target);

//! The subsurface holds no scattering species properties, so this throws
[[nodiscard]] ScatteringSpeciesMap scattering_species_property_derivatives(const SubsurfaceKeyVal& target);

//! The derivative of the temperature of a medium by a target: 1 if the target is the temperature, else 0
[[nodiscard]] Numeric temperature_derivative(const AtmKeyVal& target);

//! The derivative of the temperature of a medium by a target: 1 if the target is the temperature, else 0
[[nodiscard]] Numeric temperature_derivative(const SubsurfaceKeyVal& target);

//! The number density of the gas at an atmospheric point, from its pressure and temperature [m^-3]
[[nodiscard]] Numeric gas_number_density(const AtmPoint& point);

//! The subsurface holds no gas, so this throws
[[nodiscard]] Numeric gas_number_density(const SubsurfacePoint& point);

//! The derivative of the logarithm of the gas number density at an atmospheric point by a target
[[nodiscard]] Numeric gas_number_density_relative_derivative(const AtmPoint& point, const AtmKeyVal& target);

//! The subsurface holds no gas, so this throws
[[nodiscard]] Numeric gas_number_density_relative_derivative(const SubsurfacePoint&  point,
                                                             const SubsurfaceKeyVal& target);
}  // namespace scattering
