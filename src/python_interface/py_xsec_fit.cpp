#include <artstime.h>
#include <nanobind/stl/bind_vector.h>
#include <python_interface.h>
#include <xsec_fit.h>

#include "hpy_arts.h"

namespace Python {
void py_xsec(py::module_& m) try {
  py::class_<XsecRecord> xsec(m, "XsecRecord");
  workspace_group_interface(xsec);
  xsec.def_ro_static(
          "version", &XsecRecord::mversion, ":class:`int` The version")
      .def_rw("species",
              &XsecRecord::mspecies,
              ":class:`~pyarts.arts.Species` The species")
      .def_rw("fitcoeffs",
              &XsecRecord::mfitcoeffs,
              ":class:`~pyarts.arts.ArrayOfGriddedField2` Fit coefficients")
      .def_rw("fitminpressures",
              &XsecRecord::mfitminpressures,
              ":class:`~pyarts.arts.ArrayOfGriddedField2` Fit coefficients")
      .def_rw("fitmaxpressures",
              &XsecRecord::mfitmaxpressures,
              ":class:`~pyarts.arts.ArrayOfGriddedField2` Fit coefficients")
      .def_rw("fitmintemperatures",
              &XsecRecord::mfitmintemperatures,
              ":class:`~pyarts.arts.ArrayOfGriddedField2` Fit coefficients")
      .def_rw("fitmaxtemperatures",
              &XsecRecord::mfitmaxtemperatures,
              ":class:`~pyarts.arts.ArrayOfGriddedField2` Fit coefficients")
      .def(
          "compute_abs",
          [](XsecRecord& self,
             Numeric T,
             Numeric P,
             Numeric VMR,
             const Vector& f) {
            Vector out(f.nelem(), 0);

            self.Extract(out, f, P, T);

            out *= VMR * number_density(P, T);
            return out;
          },
          "T"_a,
          "P"_a,
          "VMR"_a,
          "f"_a,
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

      )--")
      .def("__getstate__",
           [](const XsecRecord& self) {
             return std::make_tuple(self.Version(),
                                    self.Species(),
                                    self.FitMinPressures(),
                                    self.FitMaxPressures(),
                                    self.FitMinTemperatures(),
                                    self.FitMaxTemperatures(),
                                    self.FitCoeffs());
           })
      .def("__setstate__",
           [](XsecRecord* self,
              const std::tuple<Index,
                               SpeciesEnum,
                               Vector,
                               Vector,
                               Vector,
                               Vector,
                               ArrayOfGriddedField1Named>& state) {
             new (self) XsecRecord();
             self->SetVersion(std::get<0>(state));
             self->SetSpecies(std::get<1>(state));
             self->FitMinPressures()    = std::get<2>(state);
             self->FitMaxPressures()    = std::get<3>(state);
             self->FitMinTemperatures() = std::get<4>(state);
             self->FitMaxTemperatures() = std::get<5>(state);
             self->FitCoeffs()          = std::get<6>(state);
           })
      .def(
          "to_dict",
          [](const XsecRecord& self) {
            py::dict out;

            py::dict attrs;
            attrs["creation_data"] = var_string(Time{});
            attrs["version"]       = self.Version();
            attrs["species"]       = toString(self.Species());

            std::vector<Size> r(self.FitCoeffs().size());
            std::iota(r.begin(), r.end(), 0);
            auto coords = py::dict{};

            coords["coeffs"]         = py::dict{};
            coords["coeffs"]["dims"] = "coeffs";
            coords["coeffs"]["data"] =
                std::vector<std::string>{"p00", "p10", "p01", "p20"};

            coords["bands"]         = py::dict{};
            coords["bands"]["dims"] = "bands";
            coords["bands"]["data"] = std::move(r);

            py::dict data_vars{};
            data_vars["fitminpressures"]         = py::dict{};
            data_vars["fitminpressures"]["dims"] = "bands";
            data_vars["fitminpressures"]["data"] = self.FitMinPressures();

            data_vars["fitmaxpressures"]         = py::dict{};
            data_vars["fitmaxpressures"]["dims"] = "bands";
            data_vars["fitmaxpressures"]["data"] = self.FitMaxPressures();

            data_vars["fitmintemperatures"]         = py::dict{};
            data_vars["fitmintemperatures"]["dims"] = "bands";
            data_vars["fitmintemperatures"]["data"] = self.FitMinTemperatures();

            data_vars["fitmaxtemperatures"]         = py::dict{};
            data_vars["fitmaxtemperatures"]["dims"] = "bands";
            data_vars["fitmaxtemperatures"]["data"] = self.FitMaxTemperatures();

            for (Size i = 0; i < self.FitCoeffs().size(); ++i) {
              const String band_fgrid    = var_string("band", i, "_fgrid");
              const String band_coeffs   = var_string("band", i, "_coeffs");
              const auto band_fgrid_key  = py::str(band_fgrid.c_str());
              const auto band_coeffs_key = py::str(band_coeffs.c_str());

              coords[band_fgrid_key]          = py::dict{};
              coords[band_fgrid_key]["attrs"] = py::dict{};
              coords[band_fgrid_key]["dims"]  = band_fgrid;
              coords[band_fgrid_key]["data"]  = self.FitCoeffs()[i].grid<0>();
              coords[band_fgrid_key]["attrs"]["name"] =
                  self.FitCoeffs()[i].gridname<0>();

              data_vars[band_coeffs_key] = py::dict{};
              data_vars[band_coeffs_key]["dims"] =
                  std::vector<std::string>{band_fgrid, "coeffs"};
              data_vars[band_coeffs_key]["data"]  = self.FitCoeffs()[i].data;
              data_vars[band_coeffs_key]["attrs"] = py::dict{};
              data_vars[band_coeffs_key]["attrs"]["name"] =
                  self.FitCoeffs()[i].data_name;

              if (i == 0) {
                coords["coeffs"]["attrs"] = py::dict{};
                coords["coeffs"]["attrs"]["name"] =
                    self.FitCoeffs()[i].gridname<1>();
              }
            }

            out["coords"]    = coords;
            out["attrs"]     = attrs;
            out["data_vars"] = data_vars;
            return out;
          })
      .def(
          "to_xarray",
          [](py::object& xr) {
            py::module_ xarray = py::module_::import_("xarray");
            return xarray.attr("Dataset").attr("from_dict")(
                xr.attr("to_dict")());
          },
          R"--(Convert XsecRecord to :func:`xarray.DataArray`.)--")
      .def(
          "to_netcdf",
          [](py::object& xr, py::object& f) {
            return xr.attr("to_xarray")().attr("to_netcdf")(f);
          },
          R"--(Save XsecRecord to NetCDF file.)--")
      .def_static(
          "from_dict",
          [](const py::dict& d) {
            const py::dict attrs  = d["attrs"];
            const py::dict coords = d["coords"];
            auto data_vars        = d["data_vars"];

            auto out = std::make_shared<XsecRecord>();
            out->SetVersion(py::cast<Index>(attrs["version"]));
            out->SetSpecies(py::cast<SpeciesEnum>(attrs["species"]));
            out->FitMinPressures() =
                py::cast<Vector>(data_vars["fitminpressures"]["data"]);
            out->FitMaxPressures() =
                py::cast<Vector>(data_vars["fitmaxpressures"]["data"]);
            out->FitMinTemperatures() =
                py::cast<Vector>(data_vars["fitmintemperatures"]["data"]);
            out->FitMaxTemperatures() =
                py::cast<Vector>(data_vars["fitmaxtemperatures"]["data"]);

            out->FitCoeffs().reserve(out->FitMinPressures().size());
            for (Index i = 0; i < out->FitMinPressures().size(); ++i) {
              const String band_fgrid    = var_string("band", i, "_fgrid");
              const String band_coeffs   = var_string("band", i, "_coeffs");
              const auto band_fgrid_key  = py::str(band_fgrid.c_str());
              const auto band_coeffs_key = py::str(band_coeffs.c_str());

              auto& band     = out->FitCoeffs().emplace_back();
              band.grid<0>() = py::cast<Vector>(coords[band_fgrid_key]["data"]);
              band.grid<1>() =
                  py::cast<ArrayOfString>(coords["coeffs"]["data"]);
              band.data = py::cast<Matrix>(data_vars[band_coeffs_key]["data"]);

              const py::dict coords_band_fgrid = coords[band_fgrid_key];
              const py::dict coords_band_fgrid_attrs =
                  coords[band_fgrid_key]["attrs"];
              if (coords_band_fgrid.contains("attrs") and
                  coords_band_fgrid_attrs.contains("name")) {
                band.gridname<0>() =
                    py::cast<String>(coords[band_fgrid_key]["attrs"]["name"]);
              }

              const py::dict coords_coeffs       = coords["coeffs"];
              const py::dict coords_coeffs_attrs = coords["coeffs"]["attrs"];
              if (coords_coeffs.contains("attrs") and
                  coords_coeffs_attrs.contains("name")) {
                band.gridname<1>() =
                    py::cast<String>(coords["coeffs"]["attrs"]["name"]);
              }

              const py::dict data_vars_band_coeffs = coords["coeffs"];
              const py::dict data_vars_band_coeffs_attrs =
                  coords["coeffs"]["attrs"];
              if (data_vars_band_coeffs.contains("attrs") and
                  data_vars_band_coeffs_attrs.contains("name")) {
                band.data_name = py::cast<String>(
                    data_vars[band_coeffs_key]["attrs"]["name"]);
              }
            }
            return out;
          })
      .def_static(
          "from_xarray",
          [](py::object& v) {
            return py::type<XsecRecord>().attr("from_dict")(
                v.attr("to_dict")());
          },
          R"--(Create XsecRecord from :func:`xarray.DataArray`.)--")
      .def_static(
          "from_netcdf",
          [](py::object& v) {
            py::module_ xarray = py::module_::import_("xarray");
            return py::type<XsecRecord>().attr("from_dict")(
                xarray.attr("open_dataset")(v).attr("to_dict")());
          },
          R"--(Create an XsecRecord from a NetCDF file.)--")
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

  auto a1 =
      py::bind_vector<ArrayOfXsecRecord, py::rv_policy::reference_internal>(
          m, "ArrayOfXsecRecord");
  workspace_group_interface(a1);
  vector_interface(a1);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize xsec fit\n", e.what()));
}
}  // namespace Python
