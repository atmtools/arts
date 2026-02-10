#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/bind_map.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>
#include <partfun.h>
#include <python_interface.h>

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "debug.h"
#include "hpy_arts.h"
#include "hpy_vector.h"
#include "isotopologues.h"
#include "nanobind/operators.h"
#include "species.h"
#include "species_tags.h"

std::string docs_isotopes() {
  std::ostringstream os;

  const auto isotrat = Species::isotopologue_ratiosInitFromBuiltin();

  std::print(os,
             R"(The valid isotopologues are:

.. list-table::
  :header-rows: 1

  * - Name
    - Mass (g/mol)
    - Degeneracy
    - Default Ratio
    - Predefined or Joker or Normal
)");

  for (auto& x : Species::Isotopologues) {
    std::print(
        os,
        R"(  * - {}
    - {}
    - {}
    - {}
    - {}
)",
        x.FullName(),
        x.mass,
        x.gi,
        isotrat[x],
        x.is_joker() ? "Joker" : (x.is_predefined() ? "Predefined" : "Normal"));
  }

  os << '\n';
  return os.str();
}

namespace Python {
void py_species(py::module_& m) try {
  py::class_<SpeciesIsotopologueRatios> sirs(m, "SpeciesIsotopologueRatios");
  generic_interface(sirs);
  sirs.doc() = "Isotopologue ratios for a species";
  sirs.def_static("builtin",
                  &Species::isotopologue_ratiosInitFromBuiltin,
                  "Builtin values")
      .def_ro_static("maxsize",
                     &SpeciesIsotopologueRatios::maxsize,
                     "The max size of the data\n\n.. :class:`int`")
      .def_rw("data",
              &SpeciesIsotopologueRatios::data,
              "The isotopologue ratios\n\n.. :class:`list[Numeric]`")
      .def("valueless_isotopes",
           &SpeciesIsotopologueRatios::valueless_isotopes,
           "Get a list of isotopologues without defined ratio");

  auto aose =
      py::bind_vector<ArrayOfSpeciesEnum, py::rv_policy::reference_internal>(
          m, "ArrayOfSpeciesEnum");
  generic_interface(aose);
  vector_interface(aose);
  aose.def("__init__",
           [](ArrayOfSpeciesEnum* self, const std::vector<std::string>& x) {
             new (self) ArrayOfSpeciesEnum();
             self->reserve(x.size());
             std::transform(
                 x.begin(),
                 x.end(),
                 std::back_inserter(*self),
                 [](const std::string& s) { return to<SpeciesEnum>(s); });
           });
  py::implicitly_convertible<std::vector<std::string>, ArrayOfSpeciesEnum>();

  py::class_<SpeciesIsotope> siso(m, "SpeciesIsotope");
  generic_interface(siso);
  siso.doc() = std::format("{}\n\n{}\n\n",
                           PythonWorkspaceGroupInfo<SpeciesIsotope>::desc(),
                           docs_isotopes());
  siso.def(
          "__init__",
          [](SpeciesIsotope* self, Index i) {
            try {
              auto x = Species::Isotopologues.at(i);
              new (self) SpeciesIsotope(x);
            } catch (std::exception& e) {
              throw std::runtime_error(
                  std::format("No species index {}. Error:\n{}", i, e.what()));
            }
          },
          "isot"_a = 0,
          "From position")
      .def(py::init_implicit<const std::string_view>())
      .def(
          "Q",
          [](const SpeciesIsotope& self, Numeric T) {
            return PartitionFunctions::Q(T, self);
          },
          "T"_a,
          "Partition function")
      .def_ro("spec",
              &SpeciesIsotope::spec,
              "The species\n\n.. :class:`~pyarts3.arts.SpeciesEnum`")
      .def_ro(
          "isotname",
          &SpeciesIsotope::isotname,
          "A custom name that is unique for this Species type\n\n.. :class:`str`")
      .def_ro(
          "mass",
          &SpeciesIsotope::mass,
          "The mass of the isotope in units of grams per mol. It is Nan if not defined\n\n.. :class:`float`")
      .def_ro(
          "gi",
          &SpeciesIsotope::gi,
          "The degeneracy of states of the molecule. It is -1 if not defined.\n\n.. :class:`float`")
      .def_prop_ro("name",
                   &SpeciesIsotope::FullName,
                   "The full name\n\n.. :class:`String`")
      .def_prop_ro(
          "predef",
          &SpeciesIsotope::is_predefined,
          "Check if this represents a predefined model\n\n.. :class:`bool`");
  siso.def(py::self == py::self);
  siso.def(py::self != py::self);
  siso.def(py::self <= py::self);
  siso.def(py::self >= py::self);
  siso.def(py::self < py::self);
  siso.def(py::self > py::self);
  siso.def("__hash__", [](const SpeciesIsotope& self) {
    return std::hash<SpeciesIsotope>{}(self);
  });

  auto a1 =
      py::bind_vector<ArrayOfSpeciesIsotope, py::rv_policy::reference_internal>(
          m, "ArrayOfSpeciesIsotope");
  generic_interface(a1);
  vector_interface(a1);

  py::class_<SpeciesTag> stag(m, "SpeciesTag");
  generic_interface(stag);
  stag.def_rw(
          "spec_ind", &SpeciesTag::spec_ind, "Species index\n\n.. :class:`int`")
      .def_prop_ro(
          "spec", &SpeciesTag::Spec, "Species\n\n.. :class:`SpeciesEnum`")
      .def_rw("type",
              &SpeciesTag::type,
              "Type of tag\n\n.. :class:`~pyarts3.arts.SpeciesTagType`")
      .def_rw("cia_2nd_species",
              &SpeciesTag::cia_2nd_species,
              "CIA species\n\n.. :class:`~pyarts3.arts.SpeciesEnum`")
      .def(
          "partfun",
          [](const SpeciesTag& self, Numeric T) {
            return PartitionFunctions::Q(T, self.Isotopologue());
          },
          R"--(Compute the partition function at a given temperature

Parameters
----------
  T : Numeric
    Temperature [K]

Returns
-------
  Q : Numeric
    Partition function [-]
)--",
          "T"_a)
      .def_prop_ro("full_name",
                   &SpeciesTag::FullName,
                   "The full name\n\n.. :class:`~pyarts3.arts.String`")
      .def(py::self == py::self)

      .def(py::init_implicit<std::string>());

  //////////////////////////////////////////////////////////////////////

  auto tmp1_ =
      py::bind_vector<ArrayOfSpeciesTag, py::rv_policy::reference_internal>(
          m, "ArrayOfSpeciesTag");
  vector_interface(tmp1_);
  generic_interface(tmp1_);

  //////////////////////////////////////////////////////////////////////

  auto sev = py::bind_map<SpeciesEnumVectors>(m, "SpeciesEnumVectors");
  generic_interface(sev);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize species\n{}", e.what()));
}
}  // namespace Python
