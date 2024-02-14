#include "properties.h"

ParticulateProperty toParticulatePropertyOrThrow(std::string name) {
  if (name == "number_density") return ParticulateProperty::NumberDensity;
  if (name == "mass_density") return ParticulateProperty::MassDensity;
  if (name == "d_max") return ParticulateProperty::DMax;
  if (name == "d_veq") return ParticulateProperty::DVeq;
  if (name == "intercept") return ParticulateProperty::ShapeParameter;
  if (name == "shape") return ParticulateProperty::InterceptParameter;
  throw std::runtime_error("The give bulk property name is not valid.");
}
