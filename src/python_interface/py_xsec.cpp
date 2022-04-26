#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include "python_interface.h"

#include <functional>
#include <stdexcept>

#include "hitran_xsec.h"
#include "py_macros.h"
#include "species.h"

namespace Python {
namespace details {
  struct XsecRecord {
    static std::function<py::object(py::object&)> to_xarray;
  };

  static std::function<py::object(py::object&)> to_xarray_static{[](py::object&){throw std::logic_error("Not implemented"); return py::none();}};
  std::function<py::object(py::object&)> XsecRecord::to_xarray = to_xarray_static;
} // namespace details

void py_xsec(py::module_& m) {
  m.add_object("_cleanupXsecRecord", py::capsule([](){
    details::XsecRecord::to_xarray = details::to_xarray_static;
  }));

  py::class_<details::XsecRecord>(m , "XsecRecord::details")
    .def_readwrite_static("to_xarray", &details::XsecRecord::to_xarray);

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
          }))
          .def("to_xarray", [](py::object& xr){return details::XsecRecord::to_xarray(xr);});

  PythonInterfaceWorkspaceArray(XsecRecord);
}
}  // namespace Python