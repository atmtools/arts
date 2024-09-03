#include "properties.h"

Scattering::ParticulateProperty toParticulatePropertyOrThrow(std::string name) {
  if (name == "number_density")
    return Scattering::ParticulateProperty::NumberDensity;
  if (name == "mass_density")
    return Scattering::ParticulateProperty::MassDensity;
  if (name == "d_max") return Scattering::ParticulateProperty::DMax;
  if (name == "d_veq") return Scattering::ParticulateProperty::DVeq;
  if (name == "intercept")
    return Scattering::ParticulateProperty::ShapeParameter;
  if (name == "shape")
    return Scattering::ParticulateProperty::InterceptParameter;
  throw std::runtime_error("The give bulk property name is not valid.");
}
