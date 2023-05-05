#include "gridded_fields.h"
#include "messages.h"
#include "surf.h"

void surface_fieldInit(SurfaceField &surface_field, const Verbosity &) {
  surface_field = SurfaceField{};
}

void surface_fieldSet(SurfaceField &surface_field, const Numeric &value,
                      const String &key, const Verbosity &) {
  surface_field[Surf::toKeyOrThrow(key)] = value;
}

void surface_fieldSet(SurfaceField &surface_field, const GriddedField2 &value,
                      const String &key, const Verbosity &) {
  surface_field[Surf::toKeyOrThrow(key)] = value;
}

void surface_fieldSetProp(SurfaceField &surface_field, const Numeric &value,
                          const String &key, const Verbosity &) {
  surface_field[SurfacePropertyTag{key}] = value;
}

void surface_fieldSetProp(SurfaceField &surface_field,
                          const GriddedField2 &value, const String &key,
                          const Verbosity &) {
  surface_field[SurfacePropertyTag{key}] = value;
}

void surface_fieldSetType(SurfaceField &surface_field, const Numeric &value,
                          const String &key, const Verbosity &) {
  surface_field[SurfaceTypeTag{key}] = value;
}

void surface_fieldSetType(SurfaceField &surface_field,
                          const GriddedField2 &value, const String &key,
                          const Verbosity &) {
  surface_field[SurfaceTypeTag{key}] = value;
}
