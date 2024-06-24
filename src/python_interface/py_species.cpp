#include <nanobind/stl/bind_vector.h>
#include <partfun.h>
#include <python_interface.h>

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "debug.h"
#include "hpy_arts.h"
#include "isotopologues.h"
#include "nanobind/nanobind.h"
#include "py_macros.h"
#include "species.h"
#include "species_tags.h"

std::string docs_isotopes() {
  std::ostringstream os;

  os << "The valid isotopologues are:\n\n";
  for (auto& x : Species::Isotopologues) {
    os << "- ``\"" << x.FullName() << "\"``\n";
  }
  os << '\n';
  return os.str();
}

namespace Python {
void py_species(py::module_& m) try {
  py::class_<SpeciesIsotopologueRatios>(m, "SpeciesIsotopologueRatios")
      .def(
          "__init__",
          [](SpeciesIsotopologueRatios* sir) {
            new (sir) SpeciesIsotopologueRatios(
                Species::isotopologue_ratiosInitFromBuiltin());
          },
          "Builtin values")
      .PythonInterfaceCopyValue(SpeciesIsotopologueRatios)
      // .PythonInterfaceFileIO(SpeciesIsotopologueRatios)
      .PythonInterfaceBasicRepresentation(SpeciesIsotopologueRatios)
      .def_ro_static("maxsize",
                     &SpeciesIsotopologueRatios::maxsize,
                     ":class:`int` The max size of the data")
      .def_rw("data",
              &SpeciesIsotopologueRatios::data,
              ":class:`list` The max size of the data")
      .def("__getstate__",
           [](const SpeciesIsotopologueRatios& self) {
             return std::make_tuple(self.data);
           })
      .def("__setstate__",
           [](SpeciesIsotopologueRatios* self,
              const std::tuple<
                  std::array<Numeric, SpeciesIsotopologueRatios::maxsize>>&
                  state) {
             new (self) SpeciesIsotopologueRatios{};
             self->data = std::get<0>(state);
           });

  py::bind_vector<ArrayOfSpeciesEnum>(m, "ArrayOfSpeciesEnum")
      .def("__init__",
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

  py::class_<SpeciesIsotope>(m, "SpeciesIsotope")
      .def(
          "__init__",
          [](SpeciesIsotope* self, Index i) {
            new (self) SpeciesIsotope(Species::Isotopologues.at(i));
          },
          "isot"_a = 0,
          "From position")
      .def(
          "__init__",
          [](SpeciesIsotope* self, const std::string& c) {
            new (self) SpeciesIsotope(
                Species::Isotopologues.at(Species::find_species_index(c)));
          },
          "From :class:`str`")
      .def(
          "Q",
          [](const SpeciesIsotope& self, Numeric T) {
            return PartitionFunctions::Q(T, self);
          },
          py::arg("T"),
          "Partition function")
      .def_ro("spec",
              &SpeciesIsotope::spec,
              ":class:`~pyarts.arts.Species` The species")
      .def_ro("isotname",
              &SpeciesIsotope::isotname,
              ":class:`str` A custom name that is unique for this Species type")
      .def_ro(
          "mass",
          &SpeciesIsotope::mass,
          ":class:`float` The mass of the isotope in units of grams per mol. It is Nan if not defined")
      .def_ro(
          "gi",
          &SpeciesIsotope::gi,
          ":class:`float` The degeneracy of states of the molecule. It is -1 if not defined.")
      .def_prop_ro("name",
                   &SpeciesIsotope::FullName,
                   ":class:`~pyarts.arts.String` The full name")
      .def_prop_ro("predef",
                   &Species::is_predefined_model,
                   ":class:`bool` Check if this represents a predefined model")
      .def("__getstate__",
           [](const SpeciesIsotope& self) {
             return std::make_tuple(
                 self.spec, self.isotname, self.mass, self.gi);
           })
      .def("__setstate__",
           [](SpeciesIsotope* self,
              const std::tuple<SpeciesEnum, std::string, Numeric, Index>&
                  state) {
             new (self) SpeciesIsotope(std::get<0>(state),
                                       std::get<1>(state),
                                       std::get<2>(state),
                                       std::get<3>(state));
           })
      .def(py::init_implicit<std::string>());

  auto a1 = py::bind_vector<ArrayOfSpeciesIsotope>(m, "ArrayOfSpeciesIsotope");
  workspace_group_interface(a1);

  py::class_<SpeciesTag>(m, "SpeciesTag")
      .def(
          "__init__",
          [](SpeciesTag* self, const std::string& s) {
            new (self) SpeciesTag(s);
          },
          "From :class:`str`")
      .def_rw("spec_ind", &SpeciesTag::spec_ind, ":class:`int` Species index")
      .def_rw("type",
              &SpeciesTag::type,
              ":class:`~pyarts.arts.options.SpeciesTagType` Type of tag")
      .def_rw("cia_2nd_species",
              &SpeciesTag::cia_2nd_species,
              ":class:`~pyarts.arts.Species` CIA species")
      .def("partfun",
           &SpeciesTag::Q,
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
                   ":class:`~pyarts.arts.String` The full name")
      .def(py::self == py::self)
      .def("__getstate__",
           [](const SpeciesTag& self) {
             return std::make_tuple(
                 self.spec_ind, self.type, self.cia_2nd_species);
           })
      .def("__setstate__",
           [](SpeciesTag* self,
              const std::tuple<Index, SpeciesTagType, SpeciesEnum>& state) {
             new (self) SpeciesTag{};
             self->spec_ind        = std::get<0>(state);
             self->type            = std::get<1>(state);
             self->cia_2nd_species = std::get<2>(state);
           })
      .def(py::init_implicit<std::string>());

  py::bind_vector<ArrayOfSpeciesTag>(m, "ArrayOfSpeciesTag")
      .def(
          "__init__",
          [](ArrayOfSpeciesTag* self, const std::string& x) {
            new (self) ArrayOfSpeciesTag(x);
          },
          "From :class:`str`")
      .def("__init__",
           [](ArrayOfSpeciesTag* self, const std::vector<std::string>& x) {
             new (self) ArrayOfSpeciesTag(x.size());
             std::transform(
                 x.begin(), x.end(), self->begin(), [](const std::string& s) {
                   return SpeciesTag(s);
                 });
           })
      // .def(py::init([](const std::vector<std::string>& x) {
      //   auto out = std::make_shared<ArrayOfSpeciesTag>(x.size());
      //   std::transform(
      //       x.begin(), x.end(), out->begin(), [](const std::string& s) {
      //         return SpeciesTag(s);
      //       });
      //   return out;
      // }))
      // .PythonInterfaceFileIO(ArrayOfSpeciesTag)
      .PythonInterfaceCopyValue(ArrayOfSpeciesTag)
      .PythonInterfaceBasicRepresentation(ArrayOfSpeciesTag)
      .PythonInterfaceIndexItemAccess(ArrayOfSpeciesTag)
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def("__hash__",
           [](const ArrayOfSpeciesTag& x) {
             return std::hash<ArrayOfSpeciesTag>{}(x);
           })
      .def(
          "__init__",
          [](ArrayOfSpeciesTag* self) { new (self) ArrayOfSpeciesTag(); },
          "Empty list")
      .def("__init__",
           [](ArrayOfSpeciesTag* self, Index a, const SpeciesTag& b) {
             new (self) ArrayOfSpeciesTag(a, b);
           })
      .def(
          "__init__",
          [](ArrayOfSpeciesTag* self, const std::vector<SpeciesTag>& v) {
            new (self) ArrayOfSpeciesTag(v);
          },
          "From :class:`list`")
      .def(
          "append",
          [](ArrayOfSpeciesTag& x, SpeciesTag y) { x.emplace_back(y); },
          "Appends a SpeciesTag at the end of the Array")
      .def(
          "pop",
          [](ArrayOfSpeciesTag& x) {
            SpeciesTag y = x.back();
            x.pop_back();
            return y;
          },
          "Pops a SpeciesTag from the end of the Array")
      .def(
          "__getstate",
          [](const ArrayOfSpeciesTag& x) {
            return std::make_tuple(std::vector<SpeciesTag>(x.begin(), x.end()));
          })
      .def("__setstate__",
           [](ArrayOfSpeciesTag* x,
              const std::tuple<std::vector<SpeciesTag>>& v) {
             new (x) ArrayOfSpeciesTag(std::get<0>(v));
           })
      .def("__init__",
           [](ArrayOfSpeciesTag* self, const std::vector<SpeciesTag>& x) {
             new (self) ArrayOfSpeciesTag(x);
           })
      .def(py::init_implicit<std::string>())
      .def(py::init_implicit<Array<SpeciesTag>>())
      .doc() = "List of :class:`~pyarts.arts.SpeciesTag`";

  auto b1 =
      py::bind_vector<ArrayOfArrayOfSpeciesTag>(m, "ArrayOfArrayOfSpeciesTag");
  workspace_group_interface(b1);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize species\n", e.what()));
}
}  // namespace Python
