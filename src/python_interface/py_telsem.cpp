#include <py_auto_interface.h>

#include "py_macros.h"

namespace Python {
void py_telsem(py::module_& m) {
  py::class_<TelsemAtlas>(m, "TelsemAtlas")
      .def(py::init([]() { return new TelsemAtlas{}; }))
      .PythonInterfaceCopyValue(TelsemAtlas)
      .PythonInterfaceWorkspaceVariableConversion(TelsemAtlas)
      .PythonInterfaceFileIO(TelsemAtlas)
      .PythonInterfaceBasicRepresentation(TelsemAtlas)
      .def(py::pickle(
          [](TelsemAtlas& self) {
            return py::make_tuple(self.DataCount(),
                                  self.ChannelCount(),
                                  self.Name(),
                                  self.Month(),
                                  self.Lat(),
                                  self.Cells(),
                                  self.FirstCells(),
                                  self.Emis(),
                                  self.Emis_err(),
                                  self.Correlations(),
                                  self.Classes1(),
                                  self.Classes2(),
                                  self.Cellnumber(),
                                  self.Correspondance());
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 14, "Invalid state!")
            
            auto* out = new TelsemAtlas{};
            out->DataCount() = t[0].cast<Index>();
            out->ChannelCount() = t[1].cast<Index>();
            out->Name() = t[2].cast<String>();
            out->Month() = t[3].cast<Index>();
            out->Lat() = t[4].cast<Numeric>();
            out->Cells() = t[5].cast<ArrayOfIndex>();
            out->FirstCells() = t[6].cast<ArrayOfIndex>();
            out->Emis() = t[7].cast<Matrix>();
            out->Emis_err() = t[8].cast<Matrix>();
            out->Correlations() = t[9].cast<Tensor3>();
            out->Classes1() = t[10].cast<ArrayOfIndex>();
            out->Classes2() = t[11].cast<ArrayOfIndex>();
            out->Cellnumber() = t[12].cast<ArrayOfIndex>();
            out->Correspondance() = t[13].cast<ArrayOfIndex>();

            return out;
          }));

  PythonInterfaceWorkspaceArray(TelsemAtlas);
}
}  // namespace Python