#include <workspace.h>

void atmospheric_fieldCheck(const AtmField& field) {
  for (auto&& key : field.keys()) {
    const auto& data = field[key];
    ARTS_USER_ERROR_IF(not data.ok(),
                       "Invalid atmospheric field data for key {}, data:\n{}",
                       key,
                       data);
  }
}

void surface_fieldCheck(const SurfaceField& field) {
  for (auto&& key : field.keys()) {
    const auto& data = field[key];
    ARTS_USER_ERROR_IF(not data.ok(),
                       "Invalid atmospheric field data for key {}, data:\n{}",
                       key,
                       data);
  }
}

void subsurface_fieldCheck(const SubsurfaceField& field) {
  for (auto&& key : field.keys()) {
    const auto& data = field[key];
    ARTS_USER_ERROR_IF(not data.ok(),
                       "Invalid atmospheric field data for key {}, data:\n{}",
                       key,
                       data);
  }
}

void atmospheric_fieldFixCyclicity(AtmField& field) {
  for (auto&& key : field.keys()) field[key].fix_cyclicity();
}

void surface_fieldFixCyclicity(SurfaceField& field) {
  for (auto&& key : field.keys()) field[key].fix_cyclicity();
}

void subsurface_fieldFixCyclicity(SubsurfaceField& field) {
  for (auto&& key : field.keys()) field[key].fix_cyclicity();
}
