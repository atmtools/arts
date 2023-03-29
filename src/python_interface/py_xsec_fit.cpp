#include <memory>
#include <pybind11/attr.h>
#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include "python_interface.h"

#include <functional>
#include <stdexcept>

#include "details.h"
#include "xsec_fit.h"
#include "py_macros.h"
#include "species.h"

namespace Python {
namespace details {
struct XsecRecord {
  inline static auto to_xarray{one_arg};
  inline static auto from_xarray{two_args};
  inline static auto to_netcdf{two_args};
  inline static auto from_netcdf{two_args};
  inline static auto __eq__{two_args};
};
}  // namespace details

void py_xsec(py::module_& m) {
  m.add_object("_cleanupXsecRecord", py::capsule([](){
    details::XsecRecord::to_xarray = details::one_arg;
    details::XsecRecord::from_xarray = details::two_args;
    details::XsecRecord::to_netcdf = details::two_args;
    details::XsecRecord::from_netcdf = details::two_args;
    details::XsecRecord::__eq__ = details::two_args;
  }));

  py::class_<details::XsecRecord>(m , "XsecRecord::details")
    .def_readwrite_static("to_xarray", &details::XsecRecord::to_xarray)
    .def_readwrite_static("from_xarray", &details::XsecRecord::from_xarray)
    .def_readwrite_static("to_netcdf", &details::XsecRecord::to_netcdf)
    .def_readwrite_static("from_netcdf", &details::XsecRecord::from_netcdf)
    .def_readwrite_static("__eq__", &details::XsecRecord::__eq__);

  py::class_<XsecRecord>(m, "XsecRecord")
      .def(py::init([]() { return std::make_unique<XsecRecord>(); }))
      .PythonInterfaceCopyValue(XsecRecord)
      // .PythonInterfaceWorkspaceVariableConversion(XsecRecord)
      .PythonInterfaceFileIO(XsecRecord)
      .PythonInterfaceBasicRepresentation(XsecRecord)
      .def_property("version", &XsecRecord::Version, &XsecRecord::SetVersion)
      .def_property("species", &XsecRecord::Species, &XsecRecord::SetSpecies)
      .PythonInterfaceBasicReferenceProperty(
          XsecRecord, fitcoeffs, FitCoeffs, FitCoeffs)
      .PythonInterfaceBasicReferenceProperty(
          XsecRecord, fitminpressures, FitMinPressures, FitMinPressures)
      .PythonInterfaceBasicReferenceProperty(
          XsecRecord, fitmaxpressures, FitMaxPressures, FitMaxPressures)
      .PythonInterfaceBasicReferenceProperty(XsecRecord,
                                             fitmintemperatures,
                                             FitMinTemperatures,
                                             FitMinTemperatures)
      .PythonInterfaceBasicReferenceProperty(XsecRecord,
                                             fitmaxtemperatures,
                                             FitMaxTemperatures,
                                             FitMaxTemperatures)
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

            auto out = std::make_unique<XsecRecord>();
            out->SetVersion(t[0].cast<Index>());
            out->SetSpecies(t[1].cast<Species::Species>());
            out->FitMinPressures() = t[2].cast<Vector>();
            out->FitMaxPressures() = t[3].cast<Vector>();
            out->FitMinTemperatures() = t[4].cast<Vector>();
            out->FitMaxTemperatures() = t[5].cast<Vector>();
            out->FitCoeffs() = t[6].cast<ArrayOfGriddedField2>();

            return out;
          }))
      .def(
          "to_xarray",
          [](py::object& xr) { return details::XsecRecord::to_xarray(xr); },
          py::doc(R"--(Convert XsecRecord to :func:`xarray.DataArray`.)--"))
      .def_static(
          "from_xarray",
          [](py::object& v) {
            auto t = py::type::of<XsecRecord>();
            return details::XsecRecord::from_xarray(t, v);
          },
          py::doc(R"--(Create XsecRecord from :func:`xarray.DataArray`.)--"))
      .def(
          "to_netcdf",
          [](py::object& xr, py::object& f) {
            return details::XsecRecord::to_netcdf(xr, f);
          },
          py::doc(R"--(Save XsecRecord to NetCDF file.)--"))
      .def_static(
          "from_netcdf",
          [](py::object& v) {
            auto t = py::type::of<XsecRecord>();
            return details::XsecRecord::from_netcdf(t, v);
          },
          py::doc(R"--(Create an XsecRecord from a NetCDF file.)--"))
      .def(
          "__eq__",
          [](py::object& xr, py::object& other) {
            return details::XsecRecord::__eq__(xr, other);
          },
          py::is_operator());

  PythonInterfaceWorkspaceArray(XsecRecord);
}
}  // namespace Python