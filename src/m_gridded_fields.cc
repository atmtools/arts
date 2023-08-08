#include "gridded_fields.h"
#include "mystring.h"

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void GriddedFieldGetNamer(  // WS Generic Output:
    String& name,
    // WS Generic Input:
    const T& gf) {
  // Return the name of the given GriddedField.
  name = gf.get_name();
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void ArrayOfGriddedFieldGetNamers(  // WS Generic Output:
    ArrayOfString& names,
    // WS Generic Input:
    const Array<T>& aogf) {
  // Return the name of the given GriddedField.
  names.resize(aogf.nelem());
  for (Index i = 0; i < aogf.nelem(); i++) {
    names[i] = aogf[i].get_name();
  }
}

#define gridded_field_get_namer(type) void GriddedFieldGetName(String& name, const type& gf) { \
  GriddedFieldGetNamer(name, gf); \
}

#define array_gridded_field_get_namer(type) void ArrayOfGriddedFieldGetNames(ArrayOfString& name, const ArrayOf##type& gf) { \
  ArrayOfGriddedFieldGetNamers(name, gf); \
}

gridded_field_get_namer(GriddedField1)
gridded_field_get_namer(GriddedField2)
gridded_field_get_namer(GriddedField3)
gridded_field_get_namer(GriddedField4)
gridded_field_get_namer(GriddedField5)
gridded_field_get_namer(GriddedField6)

array_gridded_field_get_namer(GriddedField1)
array_gridded_field_get_namer(GriddedField2)
array_gridded_field_get_namer(GriddedField3)
array_gridded_field_get_namer(GriddedField4)

