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
#include "py_macros.h"
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
  sirs.doc() = "Isotopologue ratios for a species";
  sirs.def(
          "__init__",
          [](SpeciesIsotopologueRatios* sir) {
            new (sir) SpeciesIsotopologueRatios(
                Species::isotopologue_ratiosInitFromBuiltin());
          },
          "Builtin values")
      .PythonInterfaceCopyValue(SpeciesIsotopologueRatios)
      // .PythonInterfaceBasicRepresentation(SpeciesIsotopologueRatios)
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

  auto aose =
      py::bind_vector<ArrayOfSpeciesEnum, py::rv_policy::reference_internal>(
          m, "ArrayOfSpeciesEnum");
  workspace_group_interface(aose);
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
  workspace_group_interface(siso);
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
      .def(
          "__init__",
          [](SpeciesIsotope* self, const std::string& c) {
            try {
              auto i = Species::find_species_index(c);
              auto x = Species::Isotopologues.at(i);
              new (self) SpeciesIsotope(x);
            } catch (std::exception& e) {
              throw std::runtime_error(
                  std::format("No species {}. Error:\n{}", c, e.what()));
            }
          },
          "From :class:`str`")
      .def(
          "Q",
          [](const SpeciesIsotope& self, Numeric T) {
            return PartitionFunctions::Q(T, self);
          },
          "T"_a,
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
                   &SpeciesIsotope::is_predefined,
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

  auto a1 =
      py::bind_vector<ArrayOfSpeciesIsotope, py::rv_policy::reference_internal>(
          m, "ArrayOfSpeciesIsotope");
  workspace_group_interface(a1);
  vector_interface(a1);

  py::class_<SpeciesTag> stag(m, "SpeciesTag");
  workspace_group_interface(stag);
  stag.def_rw("spec_ind", &SpeciesTag::spec_ind, ":class:`int` Species index")
      .def_rw("type",
              &SpeciesTag::type,
              ":class:`~pyarts.arts.SpeciesTagType` Type of tag")
      .def_rw("cia_2nd_species",
              &SpeciesTag::cia_2nd_species,
              ":class:`~pyarts.arts.Species` CIA species")
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

  //////////////////////////////////////////////////////////////////////

  auto tmp1_ =
      py::bind_vector<Array<SpeciesTag>, py::rv_policy::reference_internal>(
          m, "_ArrayOfSpeciesTag");
  vector_interface(tmp1_);

  //////////////////////////////////////////////////////////////////////

  py::class_<ArrayOfSpeciesTag> astag(m, "ArrayOfSpeciesTag");
  workspace_group_interface(astag);
  astag.def(py::init_implicit<std::string>())
      .def(py::init_implicit<Array<SpeciesTag>>())
      .def(
          "_inner",
          [](ArrayOfSpeciesTag& s) -> Array<SpeciesTag>& { return s; },
          py::rv_policy::reference_internal)
      .def("__len__",
           [](py::object& s) { return s.attr("_inner")().attr("__len__")(); })
      .def("__getitem__",
           [](py::object& s, py::object i) {
             return s.attr("_inner")().attr("__getitem__")(i);
           })
      .def("__setitem__",
           [](py::object& s, py::object i, py::object v) {
             s.attr("_inner")().attr("__setitem__")(i, v);
           })
      .def(
          "append",
          [](py::object& s, py::object v) {
            s.attr("_inner")().attr("append")(v);
          },
          "value"_a,
          "Appends a value to the end of the list")
      .def(
          "pop",
          [](py::object& s, py::object i) {
            return s.attr("_inner")().attr("pop")(i);
          },
          "i"_a = -1,
          "Pops an element from the list")
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(
          "__hash__",
          [](const ArrayOfSpeciesTag& x) {
            return std::hash<ArrayOfSpeciesTag>{}(x);
          },
          "Allows hashing")
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
      .def("__init__", [](ArrayOfSpeciesTag* s, py::list l) {
        auto x = py::cast<Array<SpeciesTag>>(l);
        new (s) ArrayOfSpeciesTag();
        for (auto& e : x) {
          s->push_back(e);
        }
      });
  py::implicitly_convertible<py::list, ArrayOfSpeciesTag>();

  //////////////////////////////////////////////////////////////////////

  auto tmp2_ = py::bind_vector<Array<Array<SpeciesTag>>,
                               py::rv_policy::reference_internal>(
      m, "_ArrayOfArrayOfSpeciesTag");
  vector_interface(tmp2_);

  //////////////////////////////////////////////////////////////////////

  auto b1 = py::bind_vector<ArrayOfArrayOfSpeciesTag,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfSpeciesTag");
  workspace_group_interface(b1);
  vector_interface(b1);
  b1.def("__init__",
         [](ArrayOfArrayOfSpeciesTag* s, const Array<Array<SpeciesTag>>& x) {
           new (s) ArrayOfArrayOfSpeciesTag(x.size());
           for (auto& v : x) {
             s->emplace_back();
             for (auto& e : v) {
               s->back().push_back(e);
             }
           }
         });
  b1.def_static(
      "all_isotopologues",
      []() {
        ArrayOfArrayOfSpeciesTag out;
        for (auto& x : Species::Isotopologues) {
          if (not x.is_normal()) continue;
          out.emplace_back(Array<SpeciesTag>{SpeciesTag(x)});
        }
        return out;
      },
      "Return a list of all species tags");
  py::implicitly_convertible<Array<Array<SpeciesTag>>,
                             ArrayOfArrayOfSpeciesTag>();
  b1.def("__init__", [](ArrayOfArrayOfSpeciesTag* s, py::list l) {
    auto x = py::cast<Array<Array<SpeciesTag>>>(l);
    new (s) ArrayOfArrayOfSpeciesTag(x.size());
    for (auto& v : x) {
      s->emplace_back();
      for (auto& e : v) {
        s->back().push_back(e);
      }
    }
  });
  py::implicitly_convertible<py::list, ArrayOfArrayOfSpeciesTag>();

  auto sev = py::bind_map<SpeciesEnumVectors>(m, "SpeciesEnumVectors");
  workspace_group_interface(sev);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize species\n{}", e.what()));
}
}  // namespace Python
