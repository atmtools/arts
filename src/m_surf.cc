#include "surf.h"

void surface_fieldSet(SurfaceField &surface_field, const Numeric &value,
                      const String &key) {
  surface_field[to<SurfaceKey>(key)] = value;
}

void surface_fieldSet(SurfaceField &surface_field, const GriddedField2 &value,
                      const String &key) {
  surface_field[to<SurfaceKey>(key)] = value;
}

void surface_fieldSetProp(SurfaceField &surface_field, const Numeric &value,
                          const String &key) {
  surface_field[SurfacePropertyTag{key}] = value;
}

void surface_fieldSetProp(SurfaceField &surface_field,
                          const GriddedField2 &value, const String &key) {
  surface_field[SurfacePropertyTag{key}] = value;
}
