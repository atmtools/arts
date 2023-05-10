#include <py_auto_interface.h>
#include <pybind11/attr.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include <memory>
#include <vector>

#include "debug.h"
#include "isotopologues.h"
#include "py_macros.h"
#include "species.h"
#include "species_tags.h"

namespace Python {
void internal_species(py::module_& m) {
  m.doc() = "For internal functionality dealing with species";

  py::class_<SpeciesTagTypeStatus>(m, "SpeciesTagTypeStatus")
      .def(py::init([](const ArrayOfArrayOfSpeciesTag& s) {
        return std::make_unique<SpeciesTagTypeStatus>(s);
      }))
      .def_readwrite("Plain", &SpeciesTagTypeStatus::Plain)
      .def_readwrite("Zeeman", &SpeciesTagTypeStatus::Zeeman)
      .def_readwrite("Predefined",
                     &SpeciesTagTypeStatus::Predefined)
      .def_readwrite("Cia", &SpeciesTagTypeStatus::Cia)
      .def_readwrite("FreeElectrons", &SpeciesTagTypeStatus::FreeElectrons)
      .def_readwrite("Particles", &SpeciesTagTypeStatus::Particles)
      .def_readwrite("XsecFit", &SpeciesTagTypeStatus::XsecFit)
      .PythonInterfaceBasicRepresentation(SpeciesTagTypeStatus);
}

void py_species(py::module_& m) {
  py::class_<SpeciesIsotopologueRatios>(m, "SpeciesIsotopologueRatios")
      .def(py::init(&Species::isotopologue_ratiosInitFromBuiltin),
           py::doc("Get the builtin values"))
      .PythonInterfaceCopyValue(SpeciesIsotopologueRatios)
      .PythonInterfaceWorkspaceVariableConversion(SpeciesIsotopologueRatios)
      .PythonInterfaceFileIO(SpeciesIsotopologueRatios)
      .PythonInterfaceBasicRepresentation(SpeciesIsotopologueRatios)
      .def_readonly_static("maxsize", &SpeciesIsotopologueRatios::maxsize)
      .def_readwrite("data", &SpeciesIsotopologueRatios::data)
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

  py::class_<Species::Species>(m, "Species")
      .def(py::init([]() { return std::make_unique<Species::Species>(); }))
      .def(py::init([](const std::string& c) {
        if (auto out = Species::fromShortName(c); good_enum(out)) return out;
        return Species::toSpeciesOrThrow(c);
      }))
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

  py::class_<ArrayOfSpecies>(m, "ArrayOfSpecies")
      .PythonInterfaceBasicRepresentation(ArrayOfSpecies)
      .PythonInterfaceArrayDefault(Species::Species);
  py::implicitly_convertible<std::vector<Species::Species>, ArrayOfSpecies>();
  py::implicitly_convertible<std::vector<std::string>, ArrayOfSpecies>();

  py::class_<SpeciesIsotopeRecord>(m, "SpeciesIsotopeRecord")
      .def(py::init([](Index i) { return Species::Isotopologues.at(i); }),
           py::arg("isot") = 0)
      .def(py::init([](const std::string& c) {
        return Species::Isotopologues.at(Species::find_species_index(c));
      }))
      .PythonInterfaceCopyValue(SpeciesIsotopeRecord)
      .PythonInterfaceBasicRepresentation(SpeciesIsotopeRecord)
      .def_readwrite("spec", &SpeciesIsotopeRecord::spec)
      .def_readwrite("isotname", &SpeciesIsotopeRecord::isotname)
      .def_readwrite("mass", &SpeciesIsotopeRecord::mass)
      .def_readwrite("gi", &SpeciesIsotopeRecord::gi)
      .def_property_readonly("name", &SpeciesIsotopeRecord::FullName)
      .def_property_readonly("predef", &Species::is_predefined_model)
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
          }));
  py::implicitly_convertible<std::string, SpeciesIsotopeRecord>();

  py::class_<ArrayOfIsotopeRecord>(m, "ArrayOfIsotopeRecord")
      .def(py::init([](bool full_list) -> ArrayOfIsotopeRecord {
             if (full_list) return ArrayOfIsotopeRecord{Species::Isotopologues};
             return ArrayOfIsotopeRecord{};
           }),
           py::arg("full_list") = false)
      .PythonInterfaceBasicRepresentation(ArrayOfIsotopeRecord)
      .PythonInterfaceIndexItemAccess(ArrayOfIsotopeRecord)
      .def(py::init([](Index a, SpeciesIsotopeRecord b) {
        return std::make_unique<ArrayOfIsotopeRecord>(a, b);
      }))
      .def(py::init([](const std::vector<SpeciesIsotopeRecord>& a) {
        return std::make_unique<ArrayOfIsotopeRecord>(a);
      }))
      .def(
          "append",
          [](ArrayOfIsotopeRecord& x, SpeciesIsotopeRecord y) {
            x.emplace_back(y);
          },
          py::doc("Appends a SpeciesIsotopeRecord at the end of the Array"))
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
      "A list of isotope records\n"
      "\n"
      "Initialize with ArrayOfIsotopeRecord(True) to get all "
      "available Arts isotopologues\n";

  py::class_<Species::TagType>(m, "SpeciesTagType")
      .def(py::init([]() { return std::make_unique<Species::TagType>(); }))
      .def(py::init([](const std::string& c) { return Species::toTagType(c); }))
      .PythonInterfaceCopyValue(Species::TagType)
      .PythonInterfaceBasicRepresentation(Species::TagType)
      .def(py::pickle(
          [](const Species::TagType& t) {
            return py::make_tuple(std::string(Species::toString(t)));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_unique<Species::TagType>(
                Species::toTagType(t[0].cast<std::string>()));
          }));
  py::implicitly_convertible<std::string, Species::TagType>();

  py::class_<SpeciesTag>(m, "SpeciesTag")
      .def(py::init([]() { return std::make_unique<SpeciesTag>(); }))
      .def(py::init([](const std::string& s) { return std::make_unique<SpeciesTag>(s); }))
      .PythonInterfaceCopyValue(SpeciesTag)
      .def_readwrite("spec_ind", &SpeciesTag::spec_ind)
      .def_readwrite("lower_freq", &SpeciesTag::lower_freq)
      .def_readwrite("upper_freq", &SpeciesTag::upper_freq)
      .def_readwrite("type", &SpeciesTag::type)
      .def_readwrite("cia_2nd_species", &SpeciesTag::cia_2nd_species)
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
      .def_property_readonly("full_name", &SpeciesTag::FullName)
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
          }));
  py::implicitly_convertible<std::string, SpeciesTag>();

  py::class_<Array<SpeciesTag>>(m, "ArrayOfSpeciesTagInternal")
      .PythonInterfaceBasicRepresentation(Array<SpeciesTag>)
      .PythonInterfaceArrayDefault(Species::Tag);

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
      .def(py::init([]() { return std::make_unique<ArrayOfSpeciesTag>(); }))
      .def(py::init([](const std::string &s) {
        return std::make_unique<ArrayOfSpeciesTag>(s);
      }))
      .def(py::init([](Index a, SpeciesTag b) {
        return std::make_unique<ArrayOfSpeciesTag>(a, b);
      }))
      .def(py::init([](const std::vector<SpeciesTag> &v) {
        return std::make_unique<ArrayOfSpeciesTag>(v);
      }))
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

  auto species = m.def_submodule("species");
  internal_species(species);
}
}  // namespace Python
