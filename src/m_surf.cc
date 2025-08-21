#include <workspace.h>

void surface_fieldSet(SurfaceField &surface_field,
                      const Numeric &value,
                      const String &key) {
  ARTS_TIME_REPORT

  surface_field[to<SurfaceKey>(key)] = value;
}

void surface_fieldSet(SurfaceField &surface_field,
                      const GeodeticField2 &value,
                      const String &key) {
  ARTS_TIME_REPORT

  surface_field[to<SurfaceKey>(key)] = value;
}

void surface_fieldSetProp(SurfaceField &surface_field,
                          const Numeric &value,
                          const String &key) {
  ARTS_TIME_REPORT

  surface_field[SurfacePropertyTag{key}] = value;
}

void surface_fieldSetProp(SurfaceField &surface_field,
                          const GeodeticField2 &value,
                          const String &key) {
  ARTS_TIME_REPORT

  surface_field[SurfacePropertyTag{key}] = value;
}
