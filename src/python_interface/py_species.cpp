#include <python_interface.h>

#include <memory>
#include <vector>

#include "debug.h"
#include "isotopologues.h"
#include "py_macros.h"
#include "species.h"
#include "species_tags.h"

namespace Python {
void py_species(py::module_& m) try {
  //! This class should behave as an `options` but we need also the "fromShortName" policy
  py::class_<Species::Species>(m, "Species")
      .def(py::init([]() { return std::make_unique<Species::Species>(); }), "Default value")
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
            return std::make_unique<Species::Species>(
                Species::toSpecies(t[0].cast<std::string>()));
          }));
  py::implicitly_convertible<std::string, Species::Species>();

  py::class_<SpeciesIsotopologueRatios>(m, "SpeciesIsotopologueRatios")
      .def(py::init(&Species::isotopologue_ratiosInitFromBuiltin),
           py::doc("Builtin values"))
      .PythonInterfaceCopyValue(SpeciesIsotopologueRatios)
      .PythonInterfaceWorkspaceVariableConversion(SpeciesIsotopologueRatios)
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
            auto out = std::make_unique<SpeciesIsotopologueRatios>();
            out->data = v;
            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(SpeciesIsotopologueRatios);

  py::class_<ArrayOfSpecies>(m, "ArrayOfSpecies")
      .PythonInterfaceBasicRepresentation(ArrayOfSpecies)
      .PythonInterfaceArrayDefault(Species::Species).doc() = "List of :class:`~pyarts.arts.Species`";
  py::implicitly_convertible<std::vector<Species::Species>, ArrayOfSpecies>();
  py::implicitly_convertible<std::vector<std::string>, ArrayOfSpecies>();

  py::class_<SpeciesIsotopeRecord>(m, "SpeciesIsotopeRecord")
      .def(py::init([](Index i) { return Species::Isotopologues.at(i); }),
           py::arg("isot") = 0, "From position")
      .def(py::init([](const std::string& c) {
        return Species::Isotopologues.at(Species::find_species_index(c));
      }), "From :class:`str`")
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
            return std::make_unique<SpeciesIsotopeRecord>(t[0].cast<Species::Species>(),
                                            t[1].cast<std::string>(),
                                            t[2].cast<Numeric>(),
                                            t[3].cast<Index>());
          })).doc() = "An isotopologue record entry";
  py::implicitly_convertible<std::string, SpeciesIsotopeRecord>();

  py::class_<ArrayOfIsotopeRecord>(m, "ArrayOfIsotopeRecord")
      .def(py::init([](bool full_list) -> ArrayOfIsotopeRecord {
             if (full_list) return ArrayOfIsotopeRecord{Species::Isotopologues};
             return ArrayOfIsotopeRecord{};
           }),
           py::arg("full_list") = false, "Empty list")
      .PythonInterfaceBasicRepresentation(ArrayOfIsotopeRecord)
      .PythonInterfaceIndexItemAccess(ArrayOfIsotopeRecord)
      .def(py::init([](const std::vector<SpeciesIsotopeRecord>& a) {
        return std::make_unique<ArrayOfIsotopeRecord>(a);
      }), "From :class:`list`")
      .def(
          "append",
          [](ArrayOfIsotopeRecord& x, SpeciesIsotopeRecord y) {
            x.emplace_back(y);
          },
          py::doc("Appends a :class:`pyarts.arts.SpeciesIsotopeRecord` at the end of the array"))
      .def(py::pickle(
          [](const ArrayOfIsotopeRecord& v) {
            auto n = v.size();
            std::vector<SpeciesIsotopeRecord> out(n);
            std::copy(v.begin(), v.end(), out.begin());
            return py::make_tuple(std::move(out));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_unique<ArrayOfIsotopeRecord>(
                t[0].cast<std::vector<SpeciesIsotopeRecord>>());
          }))
      .doc() =
      R"(A list of :class:`~pyarts.arts.IsotopeRecord`

Initialize with ``ArrayOfIsotopeRecord(True)`` to get all
available Arts isotopologues
)";

  py::class_<SpeciesTag>(m, "SpeciesTag")
      .def(py::init([]() { return std::make_unique<SpeciesTag>(); }), "Empty tag")
      .def(py::init([](const std::string& s) { return std::make_unique<SpeciesTag>(s); }), "From :class:`str`")
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
            auto out = std::make_unique<SpeciesTag>();;
            out->spec_ind = t[0].cast<Index>();
            out->lower_freq = t[1].cast<Numeric>();
            out->upper_freq = t[2].cast<Numeric>();
            out->type = t[3].cast<Species::TagType>();
            out->cia_2nd_species = t[4].cast<Species::Species>();
            return out;
          })).doc() = "The tag of a single absorption species";
  py::implicitly_convertible<std::string, SpeciesTag>();

  py::class_<Array<SpeciesTag>>(m, "_ArrayOfSpeciesTag")
      .PythonInterfaceBasicRepresentation(Array<SpeciesTag>)
      .PythonInterfaceArrayDefault(Species::Tag).doc() = "Internal array type - do not use manually ";

  py::class_<ArrayOfSpeciesTag, Array<SpeciesTag>>(m, "ArrayOfSpeciesTag")
      .PythonInterfaceFileIO(ArrayOfSpeciesTag)
      .PythonInterfaceCopyValue(ArrayOfSpeciesTag)
      .PythonInterfaceWorkspaceVariableConversion(ArrayOfSpeciesTag)
      .PythonInterfaceBasicRepresentation(ArrayOfSpeciesTag)
      .PythonInterfaceIndexItemAccess(ArrayOfSpeciesTag)
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def("__hash__",
           [](const ArrayOfSpeciesTag &x) {
             return std::hash<ArrayOfSpeciesTag>{}(x);
           })
      .def(py::init([]() { return std::make_unique<ArrayOfSpeciesTag>(); }), "Empty list")
      .def(py::init(
          [](const std::string& s) { return std::make_unique<ArrayOfSpeciesTag>(s); }), "From :class:`str`")
      .def(py::init([](Index a, SpeciesTag b) {
        return std::make_unique<ArrayOfSpeciesTag>(a, b);
      }))
      .def(py::init([](const std::vector<SpeciesTag>& v) {
        return std::make_unique<ArrayOfSpeciesTag>(v);
      }), "From :class:`list`")
      .def(
          "append",
          [](ArrayOfSpeciesTag &x, SpeciesTag y) { x.emplace_back(y); },
          py::doc("Appends a SpeciesTag at the end of the Array"))
      .def(
          "pop",
          [](ArrayOfSpeciesTag &x) {
            SpeciesTag y = x.back();
            x.pop_back();
            return y;
          },
          py::doc("Pops a SpeciesTag from the end of the Array"))
      .def(py::pickle(
          [](const ArrayOfSpeciesTag &v) {
            auto n = v.size();
            std::vector<SpeciesTag> out(n);
            std::copy(v.begin(), v.end(), out.begin());
            return py::make_tuple(std::move(out));
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_unique<ArrayOfSpeciesTag>(t[0].cast<std::vector<SpeciesTag>>());
          }))
      .PythonInterfaceWorkspaceDocumentation(ArrayOfSpeciesTag);
  py::implicitly_convertible<std::vector<SpeciesTag>, ArrayOfSpeciesTag>();
  py::implicitly_convertible<std::vector<std::string>, ArrayOfSpeciesTag>();
  py::implicitly_convertible<std::string, ArrayOfSpeciesTag>();

  PythonInterfaceWorkspaceArray(ArrayOfSpeciesTag).def(py::self == py::self);
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize species\n", e.what()));
}
}  // namespace Python
