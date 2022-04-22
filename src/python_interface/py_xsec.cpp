#include <py_auto_interface.h>

#include "hitran_xsec.h"
#include "py_macros.h"
#include "species.h"

namespace Python {
void py_xsec(py::module_& m) {
  py::class_<XsecRecord>(m, "XsecRecord")
      .def(py::init([]() { return new XsecRecord{}; }))
      .PythonInterfaceBasicRepresentation(XsecRecord)
      .def(py::pickle(
          [](const XsecRecord& self) {
            return py::make_tuple(self.Version(),
                                  self.Species(),
                                  self.FitMinPressures(),
                                  self.FitMaxPressures(),
                                  self.FitMinTemperatures(),
                                  self.FitMaxTemperatures(),
                                  self.FitCoeffs());
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 7, "Invalid state!")

            auto* out = new XsecRecord{};
            out->SetVersion(t[0].cast<Index>());
            out->SetSpecies(t[1].cast<Species::Species>());
            out->FitMinPressures() = t[2].cast<Vector>();
            out->FitMaxPressures() = t[3].cast<Vector>();
            out->FitMinTemperatures() = t[4].cast<Vector>();
            out->FitMaxTemperatures() = t[5].cast<Vector>();
            out->FitCoeffs() = t[6].cast<ArrayOfGriddedField2>();

            return out;
          }));

  PythonInterfaceWorkspaceArray(XsecRecord);
}
}  // namespace Python