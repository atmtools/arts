#include <py_auto_wsg_init.h>
#include <pybind11/attr.h>
#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

#include <algorithm>
#include <functional>
#include <memory>
#include <stdexcept>

#include "configtypes.h"
#include "details.h"
#include "matpack_math.h"
#include "physics_funcs.h"
#include "py_macros.h"
#include "python_interface.h"
#include "species.h"
#include "xsec_fit.h"

namespace Python {
void py_xsec(py::module_& m) try {
  py_staticXsecRecord(m)
      .def_property("version",
                    &XsecRecord::Version,
                    &XsecRecord::SetVersion,
                    ":class:`int` The version")
      .def_property("species",
                    &XsecRecord::Species,
                    &XsecRecord::SetSpecies,
                    ":class:`~pyarts.arts.Species` The species")
      .PythonInterfaceBasicReferenceProperty(
          XsecRecord,
          fitcoeffs,
          FitCoeffs,
          FitCoeffs,
          ":class:`~pyarts.arts.ArrayOfGriddedField2` Fit coefficients")
      .PythonInterfaceBasicReferenceProperty(
          XsecRecord,
          fitminpressures,
          FitMinPressures,
          FitMinPressures,
          ":class:`~pyarts.arts.ArrayOfGriddedField2` Fit coefficients")
      .PythonInterfaceBasicReferenceProperty(
          XsecRecord,
          fitmaxpressures,
          FitMaxPressures,
          FitMaxPressures,
          ":class:`~pyarts.arts.ArrayOfGriddedField2` Fit coefficients")
      .PythonInterfaceBasicReferenceProperty(
          XsecRecord,
          fitmintemperatures,
          FitMinTemperatures,
          FitMinTemperatures,
          ":class:`~pyarts.arts.ArrayOfGriddedField2` Fit coefficients")
      .PythonInterfaceBasicReferenceProperty(
          XsecRecord,
          fitmaxtemperatures,
          FitMaxTemperatures,
          FitMaxTemperatures,
          ":class:`~pyarts.arts.ArrayOfGriddedField2` Fit coefficients")
      .def(
          "compute_abs",
          [](XsecRecord& xsec,
             Numeric T,
             Numeric P,
             Numeric VMR,
             const Vector& f) {
            Vector out(f.nelem(), 0);

            xsec.Extract(out, f, P, T);

            out *= VMR * number_density(P, T);
            return out;
          },
          py::arg("T"),
          py::arg("P"),
          py::arg("VMR"),
          py::arg("f"),
          py::doc(
              R"--(Computes the Hitran cross-section absorption in 1/m

Parameters
----------
T : Numeric
    Temperature [K]
P : Numeric
    Pressure [Pa]
VMR : Numeric
    VMR of species [-]
f : Vector
    Frequency grid [Hz]

Returns
-------
  abs : Vector
    Absorption profile [1/m]

)--"))
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

            auto out = std::make_shared<XsecRecord>();
            out->SetVersion(t[0].cast<Index>());
            out->SetSpecies(t[1].cast<SpeciesEnum>());
            out->FitMinPressures() = t[2].cast<Vector>();
            out->FitMaxPressures() = t[3].cast<Vector>();
            out->FitMinTemperatures() = t[4].cast<Vector>();
            out->FitMaxTemperatures() = t[5].cast<Vector>();
            out->FitCoeffs() = t[6].cast<ArrayOfGriddedField1Named>();

            return out;
          }))
      .def(
          "to_dict",
          [](const XsecRecord& self) {
            py::dict out;

            py::dict attrs;
            attrs["creation_data"] = var_string(Time{});
            attrs["version"] = self.Version();
            attrs["species"] = toString(self.Species());

            std::vector<Size> r(self.FitCoeffs().size());
            std::iota(r.begin(), r.end(), 0);
            auto coords = py::dict{};

            coords["coeffs"] = py::dict{};
            coords["coeffs"]["dims"] = "coeffs";
            coords["coeffs"]["data"] =
                std::vector<std::string>{"p00", "p10", "p01", "p20"};

            coords["bands"] = py::dict{};
            coords["bands"]["dims"] = "bands";
            coords["bands"]["data"] = std::move(r);

            py::dict data_vars{};
            data_vars["fitminpressures"] = py::dict{};
            data_vars["fitminpressures"]["dims"] = "bands";
            data_vars["fitminpressures"]["data"] = self.FitMinPressures();

            data_vars["fitmaxpressures"] = py::dict{};
            data_vars["fitmaxpressures"]["dims"] = "bands";
            data_vars["fitmaxpressures"]["data"] = self.FitMaxPressures();

            data_vars["fitmintemperatures"] = py::dict{};
            data_vars["fitmintemperatures"]["dims"] = "bands";
            data_vars["fitmintemperatures"]["data"] = self.FitMinTemperatures();

            data_vars["fitmaxtemperatures"] = py::dict{};
            data_vars["fitmaxtemperatures"]["dims"] = "bands";
            data_vars["fitmaxtemperatures"]["data"] = self.FitMaxTemperatures();

            for (Size i = 0; i < self.FitCoeffs().size(); ++i) {
              const String band_fgrid = var_string("band", i, "_fgrid");
              const String band_coeffs = var_string("band", i, "_coeffs");
              const auto band_fgrid_key = py::str(band_fgrid);
              const auto band_coeffs_key = py::str(band_coeffs);

              coords[band_fgrid_key] = py::dict{};
              coords[band_fgrid_key]["dims"] = band_fgrid;
              coords[band_fgrid_key]["data"] = self.FitCoeffs()[i].grid<0>();
              coords[band_fgrid_key]["attrs"] = py::dict{};
              coords[band_fgrid_key]["attrs"]["name"] =
                  self.FitCoeffs()[i].gridname<0>();

              data_vars[band_coeffs_key] = py::dict{};
              data_vars[band_coeffs_key]["dims"] =
                  std::vector<std::string>{band_fgrid, "coeffs"};
              data_vars[band_coeffs_key]["data"] = self.FitCoeffs()[i].data;
              data_vars[band_coeffs_key]["attrs"] = py::dict{};
              data_vars[band_coeffs_key]["attrs"]["name"] =
                  self.FitCoeffs()[i].data_name;

              if (i == 0) {
                coords["coeffs"]["attrs"] = py::dict{};
                coords["coeffs"]["attrs"]["name"] =
                    self.FitCoeffs()[i].gridname<1>();
              }
            }

            out["coords"] = coords;
            out["attrs"] = attrs;
            out["data_vars"] = data_vars;
            return out;
          })
      .def(
          "to_xarray",
          [](py::object& xr) {
            py::module_ xarray = py::module_::import("xarray");
            return xarray.attr("Dataset").attr("from_dict")(
                xr.attr("to_dict")());
          },
          py::doc(R"--(Convert XsecRecord to :func:`xarray.DataArray`.)--"))
      .def(
          "to_netcdf",
          [](py::object& xr, py::object& f) {
            return xr.attr("to_xarray")().attr("to_netcdf")(f);
          },
          py::doc(R"--(Save XsecRecord to NetCDF file.)--"))
      .def_static(
          "from_dict",
          [](const py::dict& d) {
            auto attrs = d["attrs"];
            auto coords = d["coords"];
            auto data_vars = d["data_vars"];

            auto out = std::make_shared<XsecRecord>();
            out->SetVersion(attrs["version"].cast<Index>());
            out->SetSpecies(attrs["species"].cast<SpeciesEnum>());
            out->FitMinPressures() =
                data_vars["fitminpressures"]["data"].cast<Vector>();
            out->FitMaxPressures() =
                data_vars["fitmaxpressures"]["data"].cast<Vector>();
            out->FitMinTemperatures() =
                data_vars["fitmintemperatures"]["data"].cast<Vector>();
            out->FitMaxTemperatures() =
                data_vars["fitmaxtemperatures"]["data"].cast<Vector>();

            out->FitCoeffs().reserve(out->FitMinPressures().size());
            for (Index i = 0; i < out->FitMinPressures().size(); ++i) {
              const String band_fgrid = var_string("band", i, "_fgrid");
              const String band_coeffs = var_string("band", i, "_coeffs");
              const auto band_fgrid_key = py::str(band_fgrid);
              const auto band_coeffs_key = py::str(band_coeffs);

              auto& band = out->FitCoeffs().emplace_back();
              band.grid<0>() = coords[band_fgrid_key]["data"].cast<Vector>();
              band.grid<1>() = coords["coeffs"]["data"].cast<ArrayOfString>();
              band.data = data_vars[band_coeffs_key]["data"].cast<Matrix>();

              if (coords[band_fgrid_key].contains("attrs") and
                  coords[band_fgrid_key]["attrs"].contains("name")) {
                band.gridname<0>() =
                    coords[band_fgrid_key]["attrs"]["name"].cast<String>();
              }

              if (coords["coeffs"].contains("attrs") and
                  coords["coeffs"]["attrs"].contains("name")) {
                band.gridname<1>() =
                    coords["coeffs"]["attrs"]["name"].cast<String>();
              }

              if (data_vars[band_coeffs_key].contains("attrs") and
                  data_vars[band_coeffs_key]["attrs"].contains("name")) {
                band.data_name =
                    data_vars[band_coeffs_key]["attrs"]["name"].cast<String>();
              }
            }
            return out;
          })
      .def_static(
          "from_xarray",
          [](py::object& v) {
            return py::type::of<XsecRecord>().attr("from_dict")(
                v.attr("to_dict")());
          },
          py::doc(R"--(Create XsecRecord from :func:`xarray.DataArray`.)--"))
      .def_static(
          "from_netcdf",
          [](py::object& v) {
            py::module_ xarray = py::module_::import("xarray");
            return py::type::of<XsecRecord>().attr("from_dict")(
                xarray.attr("open_dataset")(v).attr("to_dict")());
          },
          py::doc(R"--(Create an XsecRecord from a NetCDF file.)--"))
      .def("__eq__", [](const XsecRecord& self, const XsecRecord& other) {
        return self.Version() == other.Version() and
               self.Species() == other.Species() and
               self.FitMinPressures() == other.FitMinPressures() and
               self.FitMaxPressures() == other.FitMaxPressures() and
               self.FitMinTemperatures() == other.FitMinTemperatures() and
               self.FitMaxTemperatures() == other.FitMaxTemperatures() and
               std::equal(
                   self.FitCoeffs().begin(),
                   self.FitCoeffs().end(),
                   other.FitCoeffs().begin(),
                   other.FitCoeffs().end(),
                   [](const auto& a, const auto& b) {
                     return a.template grid<0>() == b.template grid<0>() and
                            a.template grid<1>() == b.template grid<1>() and
                            a.data == b.data;
                   });
      });

  py_staticArrayOfXsecRecord(m);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize xsec fit\n", e.what()));
}
}  // namespace Python
