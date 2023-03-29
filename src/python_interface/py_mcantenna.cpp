#include <py_auto_interface.h>

#include "mc_antenna.h"
#include "py_macros.h"

namespace Python {
void py_mcantenna(py::module_& m) {
  py::enum_<AntennaType>(m, "AntennaType")
      .value("ANTENNA_TYPE_PENCIL_BEAM", AntennaType::ANTENNA_TYPE_PENCIL_BEAM)
      .value("ANTENNA_TYPE_GAUSSIAN", AntennaType::ANTENNA_TYPE_GAUSSIAN)
      .value("ANTENNA_TYPE_LOOKUP", AntennaType::ANTENNA_TYPE_LOOKUP)
      .PythonInterfaceCopyValue(AntennaType)
      .def(py::pickle(
          [](const AntennaType& self) {
            return py::make_tuple(static_cast<Index>(self));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            
            return static_cast<AntennaType>(t[0].cast<Index>());
          }));

  py::class_<MCAntenna>(m, "MCAntenna")
      .def(py::init([]() { return std::make_unique<MCAntenna>(); }))
      .PythonInterfaceCopyValue(MCAntenna)
      .PythonInterfaceWorkspaceVariableConversion(MCAntenna)
      .PythonInterfaceFileIO(MCAntenna)
      .PythonInterfaceBasicRepresentation(MCAntenna)
      .def_readwrite("atype", &MCAntenna::atype)
      .def_readwrite("sigma_aa", &MCAntenna::sigma_aa)
      .def_readwrite("sigma_za", &MCAntenna::sigma_za)
      .def_readwrite("aa_grid", &MCAntenna::aa_grid)
      .def_readwrite("za_grid", &MCAntenna::za_grid)
      .def_readwrite("G_lookup", &MCAntenna::G_lookup)
      .def(py::pickle(
          [](const MCAntenna& self) {
            return py::make_tuple(self.atype,
                                  self.sigma_aa,
                                  self.sigma_za,
                                  self.aa_grid,
                                  self.za_grid,
                                  self.G_lookup);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 6, "Invalid state!")

            auto out = std::make_unique<MCAntenna>();
            out->atype = t[0].cast<AntennaType>();
            out->sigma_aa = t[1].cast<Numeric>();
            out->sigma_za = t[2].cast<Numeric>();
            out->aa_grid = t[3].cast<Vector>();
            out->za_grid = t[4].cast<Vector>();
            out->G_lookup = t[5].cast<Matrix>();
            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(MCAntenna);
}
}  // namespace Python
