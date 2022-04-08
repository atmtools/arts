#include <py_auto_interface.h>

#include <pybind11/attr.h>
#include <pybind11/operators.h>

#include "py_macros.h"

namespace Python {
void py_species(py::module_& m) {
  py::class_<SpeciesIsotopologueRatios>(m, "SpeciesIsotopologueRatios")
      .def(py::init(&Species::isotopologue_ratiosInitFromBuiltin), py::doc("Get the builtin values"))
      .PythonInterfaceCopyValue(SpeciesIsotopologueRatios)
      .PythonInterfaceWorkspaceVariableConversion(SpeciesIsotopologueRatios)
      .PythonInterfaceFileIO(SpeciesIsotopologueRatios)
      .PythonInterfaceBasicRepresentation(SpeciesIsotopologueRatios)
      .def_readonly_static("maxsize", &SpeciesIsotopologueRatios::maxsize)
      .def_readwrite("data", &SpeciesIsotopologueRatios::data);

  py::class_<Species::Species>(m, "Species")
      .def(py::init([]() { return new Species::Species{}; }))
      .def(py::init([](const std::string& c) {
        if (auto out = Species::fromShortName(c); good_enum(out)) return out;
        return Species::toSpeciesOrThrow(c);
      }))
      .PythonInterfaceCopyValue(Species::Species)
      .PythonInterfaceBasicRepresentation(Species::Species);
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
      .def_readwrite("gi", &SpeciesIsotopeRecord::gi);
  py::implicitly_convertible<std::string, SpeciesIsotopeRecord>();

  py::class_<ArrayOfIsotopeRecord>(m, "ArrayOfIsotopeRecord")
      .def(py::init([](bool full_list) -> ArrayOfIsotopeRecord {
             if (full_list) return ArrayOfIsotopeRecord{Species::Isotopologues};
             return ArrayOfIsotopeRecord{};
           }),
           py::arg("full_list") = false)
      .PythonInterfaceBasicRepresentation(ArrayOfIsotopeRecord)
      .PythonInterfaceIndexItemAccess(ArrayOfIsotopeRecord)
      .def(py::init([](Index a, SpeciesIsotopeRecord b) { return new ArrayOfIsotopeRecord{a, b}; }))
      .def(py::init([](const std::vector<SpeciesIsotopeRecord>& a) { return new ArrayOfIsotopeRecord{a}; }))
      .def(
          "append",
          [](ArrayOfIsotopeRecord& x, SpeciesIsotopeRecord y) {
            x.emplace_back(std::move(y));
          },
          py::doc("Appends a SpeciesIsotopeRecord at the end of the Array"))
      .doc() =
      "A list of isotope records\n"
      "\n"
      "Initialize with ArrayOfIsotopeRecord(True) to get all "
      "available Arts isotopologues\n";

  py::class_<Species::TagType>(m, "SpeciesTagType")
      .def(py::init([]() { return new Species::TagType{}; }))
      .def(py::init([](const std::string& c) { return Species::toTagType(c); }))
      .PythonInterfaceCopyValue(Species::TagType)
      .PythonInterfaceBasicRepresentation(Species::TagType);
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
      .PythonInterfaceBasicRepresentation(SpeciesTag)
      .def(py::self == py::self);
  py::implicitly_convertible<std::string, SpeciesTag>();

  py::class_<Array<SpeciesTag>>(m, "ArrayOfSpeciesTagInternal");
  py::class_<ArrayOfSpeciesTag, Array<SpeciesTag>>(m, "ArrayOfSpeciesTag")
      .PythonInterfaceFileIO(ArrayOfSpeciesTag)
      .PythonInterfaceCopyValue(ArrayOfSpeciesTag)
      .PythonInterfaceWorkspaceVariableConversion(ArrayOfSpeciesTag)
      .PythonInterfaceBasicRepresentation(ArrayOfSpeciesTag)
      .PythonInterfaceIndexItemAccess(ArrayOfSpeciesTag)
      .def(py::self == py::self)
      .def(py::init([]() { return new ArrayOfSpeciesTag{}; }))
      .def(py::init([](const std::string& s) { return new ArrayOfSpeciesTag(s); }))
      .def(py::init([](Index a, SpeciesTag b) { return new ArrayOfSpeciesTag(a, b); }))
      .def(py::init([](const std::vector<SpeciesTag>& v) { return new ArrayOfSpeciesTag{v}; }))
      .def(
          "append",
          [](ArrayOfSpeciesTag& x, SpeciesTag y) {
            x.emplace_back(std::move(y));
          },
          py::doc("Appends a SpeciesTag at the end of the Array"))
      .def(
          "pop",
          [](ArrayOfSpeciesTag& x) {
            SpeciesTag y = x.back();
            x.pop_back();
            return y;
          },
          py::doc("Pops a SpeciesTag from the end of the Array"))
      .doc() = "The Arts ArrayOfArrayOfSpeciesTag class";
  py::implicitly_convertible<std::vector<SpeciesTag>, ArrayOfSpeciesTag>();
  py::implicitly_convertible<std::vector<py::str>, ArrayOfSpeciesTag>();
  py::implicitly_convertible<std::string, ArrayOfSpeciesTag>();

  PythonInterfaceWorkspaceArray(ArrayOfSpeciesTag)
      .def(py::self == py::self);
}
}  // namespace Python