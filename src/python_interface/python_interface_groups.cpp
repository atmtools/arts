#include "python_interface_groups.h"

namespace Python {
auto py_manual_staticArrayOfSpeciesTag(py::module& m)
    -> artsclass<ArrayOfSpeciesTag, Array<SpeciesTag>>& {
  static auto g1 = artsarray<Array<SpeciesTag>>(m, "_ArrayOfSpeciesTag");

  bool doc = false;
  if (not doc) {
    g1.doc() = "Internal array type - do not use manually ";
  }
  doc = true;

  static auto g =
      artsclass<ArrayOfSpeciesTag, Array<SpeciesTag>>(m, "ArrayOfSpeciesTag");
  return g;
}

void py_manual_groupsEarlyInit(py::module& m) { py_manual_staticArrayOfSpeciesTag(m); }
}  // namespace Python
