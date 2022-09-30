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
        return new SpeciesTagTypeStatus{s};
      }))
      .def_readwrite("Plain", &SpeciesTagTypeStatus::Plain)
      .def_readwrite("Zeeman", &SpeciesTagTypeStatus::Zeeman)
      .def_readwrite("Predefined",
                     &SpeciesTagTypeStatus::Predefined)
      .def_readwrite("Cia", &SpeciesTagTypeStatus::Cia)
      .def_readwrite("FreeElectrons", &SpeciesTagTypeStatus::FreeElectrons)
      .def_readwrite("Particles", &SpeciesTagTypeStatus::Particles)
      .def_readwrite("XsecFit", &SpeciesTagTypeStatus::XsecFit)
      .def_readwrite("NoLines", &SpeciesTagTypeStatus::NoLines)
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
            auto* out = new SpeciesIsotopologueRatios{};
            out->data = v;
            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(SpeciesIsotopologueRatios);

  py::class_<Species::Species>(m, "Species")
      .def(py::init([]() { return new Species::Species{}; }))
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
            return new Species::Species{
                Species::toSpecies(t[0].cast<std::string>())};
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
      .def(py::pickle(
          [](const SpeciesIsotopeRecord& t) {
            return py::make_tuple(t.spec, t.isotname, t.mass, t.gi);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 4, "Invalid state!")
            return new SpeciesIsotopeRecord{t[0].cast<Species::Species>(),
                                            t[1].cast<std::string>(),
                                            t[2].cast<Numeric>(),
                                            t[3].cast<Index>()};
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
        return new ArrayOfIsotopeRecord{a, b};
      }))
      .def(py::init([](const std::vector<SpeciesIsotopeRecord>& a) {
        return new ArrayOfIsotopeRecord{a};
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
            return new ArrayOfIsotopeRecord{
                t[0].cast<std::vector<SpeciesIsotopeRecord>>()};
          }))
      .doc() =
      "A list of isotope records\n"
      "\n"
      "Initialize with ArrayOfIsotopeRecord(True) to get all "
      "available Arts isotopologues\n";

  py::class_<Species::TagType>(m, "SpeciesTagType")
      .def(py::init([]() { return new Species::TagType{}; }))
      .def(py::init([](const std::string& c) { return Species::toTagType(c); }))
      .PythonInterfaceCopyValue(Species::TagType)
      .PythonInterfaceBasicRepresentation(Species::TagType)
      .def(py::pickle(
          [](const Species::TagType& t) {
            return py::make_tuple(std::string(Species::toString(t)));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return new Species::TagType{
                Species::toTagType(t[0].cast<std::string>())};
          }));
  py::implicitly_convertible<std::string, Species::TagType>();

  py::class_<SpeciesTag>(m, "SpeciesTag")
      .def(py::init([]() { return new SpeciesTag{}; }))
      .def(py::init([](const std::string& s) { return new SpeciesTag{s}; }))
      .PythonInterfaceCopyValue(SpeciesTag)
      .def_readwrite("spec_ind", &SpeciesTag::spec_ind)
      .def_readwrite("lower_freq", &SpeciesTag::lower_freq)
      .def_readwrite("upper_freq", &SpeciesTag::upper_freq)
      .def_readwrite("type", &SpeciesTag::type)
      .def_readwrite("cia_2nd_species", &SpeciesTag::cia_2nd_species)
      .def_readwrite("cia_dataset_index", &SpeciesTag::cia_dataset_index)
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
                                  t.cia_2nd_species,
                                  t.cia_dataset_index);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 6, "Invalid state!")
            auto* out = new SpeciesTag{};
            out->spec_ind = t[0].cast<Index>();
            out->lower_freq = t[1].cast<Numeric>();
            out->upper_freq = t[2].cast<Numeric>();
            out->type = t[3].cast<Species::TagType>();
            out->cia_2nd_species = t[4].cast<Species::Species>();
            out->cia_dataset_index = t[5].cast<Index>();
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
      .def(py::init([]() { return new ArrayOfSpeciesTag{}; }))
      .def(py::init(
          [](const std::string& s) { return new ArrayOfSpeciesTag(s); }))
      .def(py::init(
          [](Index a, SpeciesTag b) { return new ArrayOfSpeciesTag(a, b); }))
      .def(py::init([](const std::vector<SpeciesTag>& v) {
        return new ArrayOfSpeciesTag{v};
      }))
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
            return new ArrayOfSpeciesTag{t[0].cast<std::vector<SpeciesTag>>()};
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