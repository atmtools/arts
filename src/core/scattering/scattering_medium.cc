#include "scattering_medium.h"

#include <debug.h>
#include <physics_funcs.h>

namespace scattering {
namespace {
template <typename SpectralPoint, typename Tag>
const Vector& spectral_property(const SpectralPoint& point, const Tag& tag, Size nf) {
  const Vector& values = point[tag];
  if (values.size() != nf)
    arts::user_error("The spectral property {} has {} values, but there are {} frequencies", tag, values.size(), nf);
  return values;
}
}  // namespace

const Vector& spectral_property(const AtmSpectralPoint& point, const String& name, Size nf) {
  return spectral_property(point, AtmSpectralPropertyTag{name}, nf);
}

const Vector& spectral_property(const SubsurfaceSpectralPoint& point, const String& name, Size nf) {
  return spectral_property(point, SubsurfaceSpectralPropertyTag{name}, nf);
}

const ScatteringSpeciesMap& scattering_species_properties(const AtmPoint& point) { return point.ssprops; }

const ScatteringSpeciesMap& scattering_species_properties(const SubsurfacePoint&) {
  arts::user_error("Scattering species properties are needed, and the subsurface holds none");
}

ScatteringSpeciesMap scattering_species_property_derivatives(const AtmKeyVal& target) {
  if (const auto* property = std::get_if<ScatteringSpeciesProperty>(&target)) return {{*property, 1.0}};
  return {};
}

ScatteringSpeciesMap scattering_species_property_derivatives(const SubsurfaceKeyVal&) {
  arts::user_error("Scattering species properties are needed, and the subsurface holds none");
}

Numeric temperature_derivative(const AtmKeyVal& target) { return target == AtmKeyVal{AtmKey::t} ? 1.0 : 0.0; }

Numeric temperature_derivative(const SubsurfaceKeyVal& target) {
  return target == SubsurfaceKeyVal{SubsurfaceKey::t} ? 1.0 : 0.0;
}

Numeric gas_number_density(const AtmPoint& point) { return number_density(point.pressure, point.temperature); }

Numeric gas_number_density(const SubsurfacePoint&) {
  arts::user_error("Gas scattering needs a gas number density, and the subsurface holds no gas");
}

Numeric gas_number_density_relative_derivative(const AtmPoint& point, const AtmKeyVal& target) {
  if (target == AtmKeyVal{AtmKey::p}) return 1.0 / point.pressure;
  if (target == AtmKeyVal{AtmKey::t}) return -1.0 / point.temperature;
  return 0.0;
}

Numeric gas_number_density_relative_derivative(const SubsurfacePoint&, const SubsurfaceKeyVal&) {
  arts::user_error("Gas scattering needs a gas number density, and the subsurface holds no gas");
}
}  // namespace scattering
