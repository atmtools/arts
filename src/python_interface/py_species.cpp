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
      .def(py::init<>())
      .def(py::init([](const char* c) {
        if (auto out = Species::fromShortName(c); good_enum(out)) return out;
        return Species::toSpeciesOrThrow(c);
      }), py::arg("str").none(false))
      .PythonInterfaceCopyValue(Species::Species)
      .PythonInterfaceBasicRepresentation(Species::Species);
  py::implicitly_convertible<py::str, Species::Species>();

  py::class_<ArrayOfSpecies>(m, "ArrayOfSpecies")
      .PythonInterfaceBasicRepresentation(ArrayOfSpecies)
      .PythonInterfaceArrayDefault(Species::Species);
  py::implicitly_convertible<std::vector<Species::Species>, ArrayOfSpecies>();
  py::implicitly_convertible<std::vector<py::str>, ArrayOfSpecies>();

  py::class_<SpeciesIsotopeRecord>(m, "SpeciesIsotopeRecord")
      .def(py::init([](Index i) { return Species::Isotopologues.at(i); }),
           py::arg("isot") = 0)
      .def(py::init([](const char* c) {
        return Species::Isotopologues.at(Species::find_species_index(c));
      }), py::arg("str").none(false))
      .PythonInterfaceCopyValue(SpeciesIsotopeRecord)
      .PythonInterfaceBasicRepresentation(SpeciesIsotopeRecord)
      .def_readwrite("spec", &SpeciesIsotopeRecord::spec)
      .def_readwrite("isotname", &SpeciesIsotopeRecord::isotname)
      .def_readwrite("mass", &SpeciesIsotopeRecord::mass)
      .def_readwrite("gi", &SpeciesIsotopeRecord::gi);
  py::implicitly_convertible<py::str, SpeciesIsotopeRecord>();

  py::class_<ArrayOfIsotopeRecord>(m, "ArrayOfIsotopeRecord")
      .def(py::init([](bool full_list) -> ArrayOfIsotopeRecord {
             if (full_list) return ArrayOfIsotopeRecord{Species::Isotopologues};
             return ArrayOfIsotopeRecord{};
           }),
           py::arg("full_list") = false)
      .PythonInterfaceBasicRepresentation(ArrayOfIsotopeRecord)
      .PythonInterfaceIndexItemAccess(ArrayOfIsotopeRecord)
      .def(py::init<Index, SpeciesIsotopeRecord>())
      .def(py::init<const std::vector<SpeciesIsotopeRecord>&>())
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
      .def(py::init<>())
      .def(py::init([](const char* c) { return Species::toTagType(c); }), py::arg("str").none(false))
      .PythonInterfaceCopyValue(Species::TagType)
      .PythonInterfaceBasicRepresentation(Species::TagType);
  py::implicitly_convertible<py::str, Species::TagType>();

  py::class_<SpeciesTag>(m, "SpeciesTag")
      .def(py::init<>())
      .def(py::init<const char*>(), py::arg("str").none(false))
      .PythonInterfaceCopyValue(SpeciesTag)
      .def_readwrite("spec_ind", &SpeciesTag::spec_ind)
      .def_readwrite("lower_freq", &SpeciesTag::lower_freq)
      .def_readwrite("upper_freq", &SpeciesTag::upper_freq)
      .def_readwrite("type", &SpeciesTag::type)
      .def_readwrite("cia_2nd_species", &SpeciesTag::cia_2nd_species)
      .def_readwrite("cia_dataset_index", &SpeciesTag::cia_dataset_index)
      .PythonInterfaceBasicRepresentation(SpeciesTag)
      .def(py::self == py::self);
  py::implicitly_convertible<py::str, SpeciesTag>();

  py::class_<ArrayOfSpeciesTag>(m, "ArrayOfSpeciesTag")
      .PythonInterfaceFileIO(ArrayOfSpeciesTag)
      .PythonInterfaceCopyValue(ArrayOfSpeciesTag)
      .PythonInterfaceWorkspaceVariableConversion(ArrayOfSpeciesTag)
      .PythonInterfaceBasicRepresentation(ArrayOfSpeciesTag)
      .PythonInterfaceIndexItemAccess(ArrayOfSpeciesTag)
      .def(py::self == py::self)
      .def(py::init<>())
      .def(py::init<Index>())
      .def(py::init<const char*>(), py::arg("str").none(false))
      .def(py::init<Index, SpeciesTag>())
      .def(py::init<const std::vector<SpeciesTag>&>())
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
  py::implicitly_convertible<py::str, ArrayOfSpeciesTag>();

  PythonInterfaceWorkspaceArray(ArrayOfSpeciesTag)
      .def(py::self == py::self);
}
}  // namespace Python