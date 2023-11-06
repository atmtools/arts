#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <python_interface.h>

#include <algorithm>
#include <memory>
#include <vector>

#include "debug.h"
#include "isotopologues.h"
#include "py_macros.h"
#include "species.h"
#include "species_tags.h"

#include <partfun.h>

namespace Python {
void py_species(py::module_& m) try {
  //! This class should behave as an `options` but we need also the "fromShortName" policy
  artsclass<Species::Species>(m, "Species")
      .def(py::init([]() { return std::make_shared<Species::Species>(); }), "Default value")
      .def(py::init([](const std::string& c) {
        if (auto out = Species::fromShortName(c); good_enum(out)) return out;
        return Species::toSpeciesOrThrow(c);
      }), "From :class:`str`")
      .PythonInterfaceCopyValue(Species::Species)
      .PythonInterfaceBasicRepresentation(Species::Species)
      .def(py::pickle(
          [](const Species::Species& t) {
            return py::make_tuple(std::string(Species::toString(t)));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_shared<Species::Species>(
                Species::toSpecies(t[0].cast<std::string>()));
          }));
  py::implicitly_convertible<std::string, Species::Species>();

  artsclass<SpeciesIsotopologueRatios>(m, "SpeciesIsotopologueRatios")
      .def(py::init(&Species::isotopologue_ratiosInitFromBuiltin),
           py::doc("Builtin values"))
      .PythonInterfaceCopyValue(SpeciesIsotopologueRatios)
      .PythonInterfaceFileIO(SpeciesIsotopologueRatios)
      .PythonInterfaceBasicRepresentation(SpeciesIsotopologueRatios)
      .def_readonly_static("maxsize", &SpeciesIsotopologueRatios::maxsize, ":class:`int` The max size of the data")
      .def_readwrite("data", &SpeciesIsotopologueRatios::data, ":class:`list` The max size of the data")
      .def(py::pickle(
          [](const SpeciesIsotopologueRatios& t) {
            return py::make_tuple(t.maxsize, t.data);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 2, "Invalid state!")
            ARTS_USER_ERROR_IF(
                t[0].cast<Index>() != SpeciesIsotopologueRatios::maxsize,
                "Bad version")
            auto v = t[1].cast<
                std::array<Numeric, SpeciesIsotopologueRatios::maxsize>>();
            auto out = std::make_shared<SpeciesIsotopologueRatios>();
            out->data = v;
            return out;
          }));

  artsarray<ArrayOfSpecies>(m, "ArrayOfSpecies")
      .def(py::init([](const std::vector<std::string>& x) {
        ArrayOfSpecies out;
        std::transform(
            x.begin(),
            x.end(),
            std::back_inserter(out),
            [](const std::string& s) { return Species::toSpeciesOrThrow(s); });
        return out;
      }))
      .doc() = "List of :class:`~pyarts.arts.Species`";
  py::implicitly_convertible<std::vector<std::string>, ArrayOfSpecies>();

  artsclass<SpeciesIsotopeRecord>(m, "SpeciesIsotopeRecord")
      .def(py::init([](Index i) {
    return Species::Isotopologues.at(i); }),
           py::arg("isot") = 0, "From position")
      .def(py::init([](const std::string& c) {
        return Species::Isotopologues.at(Species::find_species_index(c));
      }), "From :class:`str`")
      .def("Q", [](const SpeciesIsotopeRecord& self, Numeric T){
        return PartitionFunctions::Q(T, self);
      }, py::arg("T"), "Partition function")
      .PythonInterfaceCopyValue(SpeciesIsotopeRecord)
      .PythonInterfaceBasicRepresentation(SpeciesIsotopeRecord)
      .def_readwrite("spec", &SpeciesIsotopeRecord::spec, ":class:`~pyarts.arts.Species` The species")
      .def_readwrite("isotname", &SpeciesIsotopeRecord::isotname, ":class:`str` A custom name that is unique for this Species type")
      .def_readwrite("mass", &SpeciesIsotopeRecord::mass, ":class:`float` The mass of the isotope in units of grams per mol. It is Nan if not defined")
      .def_readwrite("gi", &SpeciesIsotopeRecord::gi, ":class:`float` The degeneracy of states of the molecule. It is -1 if not defined.")
      .def_property_readonly("name", &SpeciesIsotopeRecord::FullName, ":class:`~pyarts.arts.String` The full name")
      .def_property_readonly("predef", &Species::is_predefined_model, ":class:`bool` Check if this represents a predefined model")
      .def(py::pickle(
          [](const SpeciesIsotopeRecord& t) {
            return py::make_tuple(t.spec, t.isotname, t.mass, t.gi);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 4, "Invalid state!")
            return std::make_shared<SpeciesIsotopeRecord>(t[0].cast<Species::Species>(),
                                            t[1].cast<std::string>(),
                                            t[2].cast<Numeric>(),
                                            t[3].cast<Index>());
          })).doc() = "An isotopologue record entry";
  py::implicitly_convertible<std::string, SpeciesIsotopeRecord>();

  artsarray<ArrayOfIsotopeRecord>(m, "ArrayOfIsotopeRecord")
      .doc() =
      R"(A list of :class:`~pyarts.arts.IsotopeRecord`

Initialize with ``ArrayOfIsotopeRecord(True)`` to get all
available Arts isotopologues
)";

  artsclass<SpeciesTag>(m, "SpeciesTag")
      .def(py::init([]() { return std::make_shared<SpeciesTag>(); }), "Empty tag")
      .def(py::init([](const std::string& s) { return std::make_shared<SpeciesTag>(s); }), "From :class:`str`")
      .PythonInterfaceCopyValue(SpeciesTag)
      .def_readwrite("spec_ind", &SpeciesTag::spec_ind, ":class:`int` Species index")
      .def_readwrite("lower_freq", &SpeciesTag::lower_freq, ":class:`float` Lower cutoff frequency")
      .def_readwrite("upper_freq", &SpeciesTag::upper_freq, ":class:`float` Upper cutoff frequency")
      .def_readwrite("type", &SpeciesTag::type, ":class:`~pyarts.arts.options.SpeciesTagType` Type of tag")
      .def_readwrite("cia_2nd_species", &SpeciesTag::cia_2nd_species, ":class:`~pyarts.arts.Species` CIA species")
      .def("partfun",
           py::vectorize(&SpeciesTag::Q),
           py::doc(R"--(Compute the partition function at a given temperature

Parameters
----------
  T : Numeric
    Temperature [K]

Returns
-------
  Q : Numeric
    Partition function [-]
)--"),
           py::arg("T"))
      .def_property_readonly("full_name", &SpeciesTag::FullName, ":class:`~pyarts.arts.String` The full name")
      .PythonInterfaceBasicRepresentation(SpeciesTag)
      .def(py::self == py::self)
      .def(py::pickle(
          [](const SpeciesTag& t) {
            return py::make_tuple(t.spec_ind,
                                  t.lower_freq,
                                  t.upper_freq,
                                  t.type,
                                  t.cia_2nd_species);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 5, "Invalid state!")
            auto out = std::make_shared<SpeciesTag>();;
            out->spec_ind = t[0].cast<Index>();
            out->lower_freq = t[1].cast<Numeric>();
            out->upper_freq = t[2].cast<Numeric>();
            out->type = t[3].cast<Species::TagType>();
            out->cia_2nd_species = t[4].cast<Species::Species>();
            return out;
          })).doc() = "The tag of a single absorption species";
  py::implicitly_convertible<std::string, SpeciesTag>();

  artsarray<Array<SpeciesTag>>(m, "_ArrayOfSpeciesTag").doc() =
      "Internal array type - do not use manually ";

  artsclass<ArrayOfSpeciesTag, Array<SpeciesTag>>(m, "ArrayOfSpeciesTag")
      .def(py::init([](const std::string& x) {
        return std::make_shared<ArrayOfSpeciesTag>(x);
      }))
      .def(py::init([](const std::vector<std::string>& x) {
        auto out = std::make_shared<ArrayOfSpeciesTag>(x.size());
        std::transform(
            x.begin(), x.end(), out->begin(), [](const std::string& s) {
              return SpeciesTag(s);
            });
        return out;
      }))
      .PythonInterfaceFileIO(ArrayOfSpeciesTag)
      .PythonInterfaceCopyValue(ArrayOfSpeciesTag)
      .PythonInterfaceWorkspaceVariableConversion(ArrayOfSpeciesTag)
      .PythonInterfaceBasicRepresentation(ArrayOfSpeciesTag)
      .PythonInterfaceIndexItemAccess(ArrayOfSpeciesTag)
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def("__hash__",
           [](const ArrayOfSpeciesTag& x) {
             return std::hash<ArrayOfSpeciesTag>{}(x);
           })
      .def(py::init([]() { return std::make_shared<ArrayOfSpeciesTag>(); }),
           "Empty list")
      .def(py::init([](const std::string& s) {
             return std::make_shared<ArrayOfSpeciesTag>(s);
           }),
           "From :class:`str`")
      .def(py::init([](Index a, SpeciesTag b) {
        return std::make_shared<ArrayOfSpeciesTag>(a, b);
      }))
      .def(py::init([](const std::vector<SpeciesTag>& v) {
             return std::make_shared<ArrayOfSpeciesTag>(v);
           }),
           "From :class:`list`")
      .def(
          "append",
          [](ArrayOfSpeciesTag& x, SpeciesTag y) { x.emplace_back(y); },
          py::doc("Appends a SpeciesTag at the end of the Array"))
      .def(
          "pop",
          [](ArrayOfSpeciesTag& x) {
            SpeciesTag y = x.back();
            x.pop_back();
            return y;
          },
          py::doc("Pops a SpeciesTag from the end of the Array"))
      .def(py::pickle(
          [](const ArrayOfSpeciesTag& v) {
            auto n = v.size();
            std::vector<SpeciesTag> out(n);
            std::copy(v.begin(), v.end(), out.begin());
            return py::make_tuple(std::move(out));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_shared<ArrayOfSpeciesTag>(t[0].cast<std::vector<SpeciesTag>>());
          }))
      .def(py::init([](const std::vector<SpeciesTag>& x) {
        return std::make_shared<ArrayOfSpeciesTag>(x);
      }))
      .def(py::init([](const std::vector<py::object>& x) {
        auto out = std::make_shared<ArrayOfSpeciesTag>(x.size());
        std::transform(x.begin(), x.end(), out->begin(), [](const auto& s) {
          return py::cast<SpeciesTag>(s);
        });
        return out;
      }))
      .def(py::init([](const py::list& x) {
        auto out = std::make_shared<ArrayOfSpeciesTag>(x.size());
        std::transform(x.begin(), x.end(), out->begin(), [](const auto& s) {
          return py::cast<SpeciesTag>(s);
        });
        return out;
      }))
      .PythonInterfaceWorkspaceDocumentation(ArrayOfSpeciesTag);
  py::implicitly_convertible<std::string, ArrayOfSpeciesTag>();
  py::implicitly_convertible<Array<SpeciesTag>, ArrayOfSpeciesTag>();
  py::implicitly_convertible<Array<py::object>, ArrayOfSpeciesTag>();
  py::implicitly_convertible<py::list, ArrayOfSpeciesTag>();

  artsarray<ArrayOfArrayOfSpeciesTag>(m, "ArrayOfArrayOfSpeciesTag")
      .PythonInterfaceFileIO(ArrayOfArrayOfSpeciesTag)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfArrayOfSpeciesTag);
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize species\n", e.what()));
}
}  // namespace Python
