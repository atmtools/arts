#include <atm.h>
#include <atm_path.h>
#include <debug.h>
#include <hpy_arts.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <python_interface.h>
#include <quantum.h>
#include <species_tags.h>

#include <stdexcept>
#include <unordered_map>

namespace Python {
void py_atm(py::module_ &m) try {
  py::class_<Atm::Data> atmdata(m, "AtmData");
  generic_interface(atmdata);
  atmdata.def(py::init_implicit<GeodeticField3>())
      .def(py::init_implicit<Numeric>())
      .def(py::init_implicit<Atm::FunctionalData>())
      .def(py::init_implicit<GeodeticField3>())
      .def(
          "__init__",
          [](Atm::Data *a, const GriddedField3 &v) {
            new (a) Atm::Data(GeodeticField3(v));
          },
          "v"_a,
          "Initialize with a gridded field")
      .def(
          "__init__",
          [](Atm::Data *a, const SortedGriddedField3 &v) {
            new (a) Atm::Data(GeodeticField3{
                .data_name  = v.data_name,
                .data       = v.data,
                .grid_names = v.grid_names,
                .grids = {v.grid<0>(), v.grid<1>().vec(), v.grid<2>().vec()}});
          },
          "v"_a,
          "Initialize with a gridded field")
      .def_rw(
          "data",
          &Atm::Data::data,
          "The data.\n\n.. :class:`~pyarts3.arts.GeodeticField3`\n\n.. :class:`~pyarts3.arts.Numeric`\n\n.. :class:`~pyarts3.arts.NumericTernaryOperator`")
      .def_rw(
          "alt_upp",
          &Atm::Data::alt_upp,
          "Upper altitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw(
          "alt_low",
          &Atm::Data::alt_low,
          "Lower altitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw(
          "lat_upp",
          &Atm::Data::lat_upp,
          "Upper latitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw(
          "lat_low",
          &Atm::Data::lat_low,
          "Lower latitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw(
          "lon_upp",
          &Atm::Data::lon_upp,
          "Upper longitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw(
          "lon_low",
          &Atm::Data::lon_low,
          "Lower longitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def(
          "set_extrapolation",
          [](Atm::Data &self, InterpolationExtrapolation x) {
            self.alt_upp = x;
            self.alt_low = x;
            self.lat_upp = x;
            self.lat_low = x;
            self.lon_upp = x;
            self.lon_low = x;
          },
          "extrapolation"_a,
          "Set the extrapolation for all dimensions")
      .def(
          "__call__",
          [](const Atm::Data &d, Numeric alt, Numeric lat, Numeric lon) {
            return d.at(alt, lat, lon);
          },
          "alt"_a,
          "lat"_a,
          "lon"_a,
          "Get a point of data at the position")
      .def(
          "ws",
          [](const Atm::Data &d, Numeric alt, Numeric lat, Numeric lon) {
            return d.flat_weight(alt, lat, lon);
          },
          "alt"_a,
          "lat"_a,
          "lon"_a,
          "Get the weights of neighbors at a position")
      .def_prop_ro("data_type",
                   &Atm::Data::data_type,
                   "The data type\n\n.. :class:`~pyarts3.arts.String`")
      .def(
          "regrid",
          [](AtmData &data,
             const AscendingGrid &alt,
             const LatGrid &lat,
             const LonGrid &lon,
             const InterpolationExtrapolation &extrapolation) {
            data.regrid(alt, lat, lon);
            data.alt_upp = extrapolation;
            data.alt_low = extrapolation;
            data.lat_upp = extrapolation;
            data.lat_low = extrapolation;
            data.lon_upp = extrapolation;
            data.lon_low = extrapolation;
          },
          "alt"_a,
          "lat"_a,
          "lon"_a,
          "extrapolation"_a = InterpolationExtrapolation::Nearest,
          R"(Regrid the data to a new grid.

This will convert the data to a *GeodeticField3* if it is not already.
It will not respect the existing extrapolation settings, so these must be set
manually after calling this method.

Parameters
----------
alt : AscendingGrid
    The new altitude grid.
lat : LatGrid
    The new latitude grid.
lon : LonGrid
    The new longitude grid.
extrapolation : InterpolationExtrapolation
    The extrapolation method to use.  Defaults to Nearest
)");
  py::implicitly_convertible<Atm::FunctionalData::func_t, Atm::Data>();
  py::implicitly_convertible<GriddedField3, Atm::Data>();
  py::implicitly_convertible<SortedGriddedField3, Atm::Data>();

  auto aad =
      py::bind_vector<std::vector<Atm::Data>,
                      py::rv_policy::reference_internal>(m, "ArrayOfAtmData");
  generic_interface(aad);
  aad.doc() = "A list of :class:`~pyarts3.arts.AtmData`";
  vector_interface(aad);

  auto pnt = py::class_<AtmPoint>(m, "AtmPoint");
  generic_interface(pnt);

  auto fld = py::class_<AtmField>(m, "AtmField");
  generic_interface(fld);

  fld.def(
      "regrid",
      [](AtmField &atm,
         const AscendingGrid &alt,
         const LatGrid &lat,
         const LonGrid &lon,
         const InterpolationExtrapolation extrapolation,
         bool do_atmkeys,
         bool do_species,
         bool do_isotopologues,
         bool do_nlte,
         bool do_scattering) {
        if (do_atmkeys) {
          for (auto &val : atm.other | stdv::values) {
            val.regrid(alt, lat, lon);
            val.alt_upp = extrapolation;
            val.alt_low = extrapolation;
            val.lat_upp = extrapolation;
            val.lat_low = extrapolation;
            val.lon_upp = extrapolation;
            val.lon_low = extrapolation;
          }
        }

        if (do_species) {
          for (auto &val : atm.specs | stdv::values) {
            val.regrid(alt, lat, lon);
            val.alt_upp = extrapolation;
            val.alt_low = extrapolation;
            val.lat_upp = extrapolation;
            val.lat_low = extrapolation;
            val.lon_upp = extrapolation;
            val.lon_low = extrapolation;
          }
        }

        if (do_isotopologues) {
          for (auto &val : atm.isots | stdv::values) {
            val.regrid(alt, lat, lon);
            val.alt_upp = extrapolation;
            val.alt_low = extrapolation;
            val.lat_upp = extrapolation;
            val.lat_low = extrapolation;
            val.lon_upp = extrapolation;
            val.lon_low = extrapolation;
          }
        }

        if (do_nlte) {
          for (auto &val : atm.nlte | stdv::values) {
            val.regrid(alt, lat, lon);
            val.alt_upp = extrapolation;
            val.alt_low = extrapolation;
            val.lat_upp = extrapolation;
            val.lat_low = extrapolation;
            val.lon_upp = extrapolation;
            val.lon_low = extrapolation;
          }
        }

        if (do_scattering) {
          for (auto &val : atm.ssprops | stdv::values) {
            val.regrid(alt, lat, lon);
            val.alt_upp = extrapolation;
            val.alt_low = extrapolation;
            val.lat_upp = extrapolation;
            val.lat_low = extrapolation;
            val.lon_upp = extrapolation;
            val.lon_low = extrapolation;
          }
        }
      },
      "alt"_a,
      "lat"_a,
      "lon"_a,
      "extrapolation"_a    = InterpolationExtrapolation::Nearest,
      "do_atmkeys"_a       = true,
      "do_species"_a       = true,
      "do_isotopologues"_a = false,
      "do_nlte"_a          = true,
      "do_scattering"_a    = true,
      R"(Regrid all fields to a new grid.

This will convert the data to a *GeodeticField3* if it is not already.
It will not respect the existing extrapolation settings, so these must be set
manually after calling this method.

Parameters
----------
alt : AscendingGrid
    The new altitude grid.
lat : LatGrid
    The new latitude grid.
lon : LonGrid
    The new longitude grid.
extrapolation : InterpolationExtrapolation
    The extrapolation method to use.  The default is Nearest.
do_atmkeys : bool
    If true, regrid the basic atmospheric keys (temperature, pressure, magnetic field, wind).  Default is true.
do_species : bool
    If true, regrid the species VMRs.  Default is true.
do_isotopologues : bool
    If true, regrid the isotopologue ratios.  Default is false.
do_nlte : bool
    If true, regrid the NLTE values.  Default is true.
do_scattering : bool
    If true, regrid the scattering species properties.  Default is true.
)");

  pnt.def_rw("temperature",
             &AtmPoint::temperature,
             "Temperature [K]\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("pressure",
              &AtmPoint::pressure,
              "Pressure [Pa]\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("mag",
              &AtmPoint::mag,
              "Magnetic field [T]\n\n.. :class:`~pyarts3.arts.Vector3`")
      .def_rw("wind",
              &AtmPoint::wind,
              "Wind field [m/s]\n\n.. :class:`~pyarts3.arts.Vector3`")
      .def_rw(
          "nlte",
          &AtmPoint::nlte,
          "NLTE data\n\n.. :class:`dict[~pyarts3.arts.QuantumIdentifier, ~pyarts3.arts.Numeric]`")
      .def_rw(
          "specs",
          &AtmPoint::specs,
          "Species data\n\n.. :class:`dict[~pyarts3.arts.SpeciesEnum, ~pyarts3.arts.Numeric]`")
      .def_rw(
          "isots",
          &AtmPoint::isots,
          "Isotopologue ratio data\n\n.. :class:`dict[~pyarts3.arts.SpeciesIsotope, ~pyarts3.arts.Numeric]`")
      .def_rw(
          "ssprops",
          &AtmPoint::ssprops,
          "Scattering species properties data\n\n.. :class:`dict[~pyarts3.arts.ScatteringSpeciesProperty, ~pyarts3.arts.Numeric]`")
      .def(
          "number_density",
          [](AtmPoint &atm, std::variant<SpeciesEnum, SpeciesIsotope> key) {
            return std::visit([&atm]<typename T>(
                                  const T &s) { return atm.number_density(s); },
                              key);
          },
          "spec"_a,
          R"(Get the number density [1 / m :sup:`3`] of a species or of a species isotopologue.

The number density is computed as:

.. math::
  n_d = \frac{xp}{kT}

where :math:`x` is the VMR of the species if ``spec`` is a :class:`SpeciesEnum` (in code: ``self[spec]``) or the VMR of the
species isotopolgue if ``spec`` is a :class:`SpeciesIsotope` (in code: ``self[spec] * self[spec.spec]``).  :math:`p` and :math:`T` are the
pressure [Pa] and temperature [K] of the atmospheric point, respectively.  :math:`k` is the
Boltzmann constant.

Parameters
----------
spec : SpeciesEnum | SpeciesIsotope
    The species or isotopologue that is considered.

Returns
-------
nd : float
    The number density [1 / :math:`\textrm{m}^3`] of the species or isotopologue.
)")
      .def("__getitem__",
           [](AtmPoint &self, const AtmKeyVal &key) {
             if (self.contains(key)) return self[key];
             throw py::key_error(std::format("{}", key).c_str());
           })
      .def("__setitem__",
           [](AtmPoint &self, const AtmKeyVal &key, Numeric x) {
             self[key] = x;
           })
      .def("__contains__",
           [](const AtmPoint &self, const AtmKeyVal &x) {
             return self.contains(x);
           })
      .def("keys",
           &AtmPoint::keys,
           "Available keys",
           "keep_basic"_a   = true,
           "keep_specs"_a   = true,
           "keep_isots"_a   = false,
           "keep_nlte"_a    = false,
           "keep_ssprops"_a = true)
      .def(
          "no_isotopologues",
          [](AtmPoint in) {
            in.isots.clear();
            return in;
          },
          "Returns an atmospheric point without isotopologue ratios.")
      .def(
          "flat_values",
          [](const AtmPoint &in) {
            Vector out;
            for (auto &&key : in.keys()) {
              out.push_back(in[key]);
            }
            return out;
          },
          "Returns a flat list of values.");

  fld.def(
         "__init__",
         [](AtmField *x, Numeric toa, IsoRatioOption iso) {
           new (x) AtmField(iso);
           x->top_of_atmosphere = toa;
         },
         "toa"_a = 100e3,
         "iso"_a = IsoRatioOption::Builtin)
      .def(
          "species_keys",
          [](const AtmField &atm) {
            return atm.specs | stdv::keys |
                   stdr::to<std::vector<SpeciesEnum>>();
          },
          "Species keys")
      .def(
          "__call__",
          [](const AtmField &atm,
             const Numeric &h,
             const Numeric &lat,
             const Numeric &lon) { return atm.at(h, lat, lon); },
          "h"_a,
          "lat"_a,
          "lon"_a,
          "Get the data at a point")
      .def(
          "__call__",
          [](const AtmField &atm,
             const Vector &hv,
             const Vector &latv,
             const Vector &lonv) {
            const Size N = hv.size();
            if (latv.size() != lonv.size() or N != latv.size())
              throw std::logic_error(std::format(R"(Not same size:
  h:   {:B,} (size: {})
  lat: {:B,} (size: {})
  lon: {:B,} (size: {})
)",
                                                 hv,
                                                 hv.size(),
                                                 latv,
                                                 latv.size(),
                                                 lonv,
                                                 lonv.size()));
            ArrayOfAtmPoint out;
            out.reserve(N);
            for (Size i = 0; i < N; i++)
              out.emplace_back(atm.at(hv[i], latv[i], lonv[i]));
            return out;
          },
          "h"_a,
          "lat"_a,
          "lon"_a,
          "Get the data as a list")
      .def_rw("top_of_atmosphere",
              &AtmField::top_of_atmosphere,
              "Top of the atmosphere [m]\n\n.. :class:`~pyarts3.arts.Numeric`");

  m.def(
      "frequency_shift",
      [](AscendingGrid freq_grid,
         const AtmPoint &atm_point,
         const PropagationPathPoint &ray_point) {
        Vector3 x;
        freq_gridWindShift(freq_grid, x, atm_point, ray_point);
        return freq_grid;
      },
      "freq_grid"_a,
      "atm_point"_a,
      "ray_point"_a,
      R"(Get the frequency-shifted frequency grid at a point in the atmosphere.

Parameters
----------
  freq_grid : AscendingGrid
    The frequency grid to shift.
  atm_point : AtmPoint
    The point in the atmosphere.
  ray_point : PropagationPathPoint
    The point along the ray path.

Return
------
  freq_grid : AscendingGrid
      The frequency-shifted frequency grid.
)");

  fld.def(
      "to_dict",
      [](const AtmField &atm,
         bool core,
         bool specs,
         bool isots,
         bool nlte,
         bool ssprops) {
        std::unordered_map<Atm::KeyVal, Atm::FieldData> out;

        const auto keys = atm.keys();

        for (const Atm::KeyVal &key : keys) {
          if (std::holds_alternative<SpeciesEnum>(key)) {
            if (specs) out[key] = atm[key].data;
          } else if (std::holds_alternative<SpeciesIsotope>(key)) {
            if (isots) out[key] = atm[key].data;
          } else if (std::holds_alternative<QuantumLevelIdentifier>(key)) {
            if (nlte) out[key] = atm[key].data;
          } else if (std::holds_alternative<ScatteringSpeciesProperty>(key)) {
            if (ssprops) out[key] = atm[key].data;
          } else if (std::holds_alternative<AtmKey>(key)) {
            if (core) out[key] = atm[key].data;
          } else {
            throw std::runtime_error("Unknown key type.");
          }
        }

        return out;
      },
      "core"_a    = true,
      "specs"_a   = true,
      "isots"_a   = true,
      "nlte"_a    = true,
      "ssprops"_a = true,
      R"(Convert an atmospheric field to a dictionary.

.. tip::
  The function :func:`stringify_keys` may be used to convert the dictionary keys to strings.
  This allows for easier manipulation of the data of the dictionary.

Parameters
----------
  core : bool, optional
    If True, the core atmospheric data will be included (i.e., temperature, pressure, etc).  Default is True.
  specs : bool, optional
    If True, the species VMR data will be included.  Default is True.
  isots : bool, optional
    If True, the isotopologue ratio data will be included.  Default is True.
  nlte : bool, optional
    If True, the NLTE data will be included.  Default is True.
  ssprops : bool, optional
    If True, the scattering species properties data will be included.  Default is True.

Return
------
  dict : dict
    The atmospheric field data.
)");

  fld.def(
      "keys",
      [](const AtmField &atm,
         bool core,
         bool specs,
         bool isots,
         bool nlte,
         bool ssprops) {
        using std::back_inserter;

        std::vector<Atm::KeyVal> out;
        if (core) stdr::copy(atm.other | stdv::keys, back_inserter(out));
        if (specs) stdr::copy(atm.specs | stdv::keys, back_inserter(out));
        if (isots) stdr::copy(atm.isots | stdv::keys, back_inserter(out));
        if (nlte) stdr::copy(atm.nlte | stdv::keys, back_inserter(out));
        if (ssprops) stdr::copy(atm.ssprops | stdv::keys, back_inserter(out));

        return out;
      },
      "core"_a    = true,
      "specs"_a   = true,
      "isots"_a   = true,
      "nlte"_a    = true,
      "ssprops"_a = true,
      R"(Get a list of the keys from an atmospheric field.

>>> from pyarts3.arts import AtmField
>>> field = AtmField()
>>> field["t"] = 273
>>> k = AtmField.keys(isots=False)  # Get list of keys ignoring isotopologue ratios
["t"]

Parameters
----------
  core : bool, optional
    If True, the core atmospheric keys will be included (i.e., temperature, pressure, etc).  Default is True.
  specs : bool, optional
    If True, the species VMR keys will be included.  Default is True.
  isots : bool, optional
    If True, the isotopologue ratio keys will be included.  Default is True.
  nlte : bool, optional
    If True, the NLTE keys will be included.  Default is True.
  ssprops : bool, optional
    If True, the scattering species property keys will be included.  Default is True.

Return
------
  keys : list
    A list of the keys in the atmospheric field.
)");

  const auto update_field =
      [](AtmField &atm,
         const std::unordered_map<Atm::KeyVal, Atm::FieldData> &dict,
         const InterpolationExtrapolation &extrap) {
        for (auto &[key, data] : dict) {
          const bool existed = atm.contains(key);
          atm[key].data      = data;
          if (not existed) {
            atm[key].alt_low = extrap;
            atm[key].alt_upp = extrap;
            atm[key].lat_low = extrap;
            atm[key].lat_upp = extrap;
            atm[key].lon_low = extrap;
            atm[key].lon_upp = extrap;
          }
        }
      };

  fld.def_static(
      "from_dict",
      [&update_field](
          const std::unordered_map<Atm::KeyVal, Atm::FieldData> &dict,
          const Numeric &toa,
          const IsoRatioOption &iso,
          const InterpolationExtrapolation &extrap) {
        AtmField out(iso);
        out.top_of_atmosphere = toa;
        update_field(out, dict, extrap);
        return out;
      },
      "data"_a,
      "toa"_a    = 100e3,
      "iso"_a    = IsoRatioOption::Builtin,
      "extrap"_a = InterpolationExtrapolation::Nearest,
      R"(Create an atmospheric field from a dictionary.

.. tip::
  Each key-type of the dictionary is constructible from a :class:`str`.  So using a dictionary with string keys is possible.

Parameters
----------
  data : dict
    A dictionary of atmospheric keys and corresponding data.
  toa : Numeric, optional
    The top of the atmosphere.  Default is 100e3.
  iso : IsoRatioOption, optional
    The isotopologue ratio option to use.  Default is "Builtin". Use "None" to create field without isotopologue ratios.
  extrap : InterpolationExtrapolation, optional
    The extrapolation method to use for the new keys.  Default is "Nearest".

Return
------
  atm_field : AtmField
    The atmospheric field created from the dictionary.
)");

  fld.def("update",
          update_field,
          "data"_a,
          "extrap"_a = InterpolationExtrapolation::Nearest,
          R"(Update the atmospheric field with dictionary values.

.. tip::
  Each key-type of the dictionary is constructible from a :class:`str`.  So using a dictionary with string keys is possible.

Parameters
----------
  data : dict
    A dictionary of atmospheric keys and corresponding data.
  extrap : InterpolationExtrapolation, optional
    The extrapolation method to use for the new keys.  Default is "Nearest".  Ignored by existing keys.
)");

  fld.def(
      "as_gridded",
      &AtmField::gridded,
      "alt"_a,
      "lat"_a,
      "lon"_a,
      R"(Convert all fields of the input atmospheric field to gridded fields.

Parameters
----------
alt : AscendingGrid
  The altitude grid.
lat : AscendingGrid
  The latitude grid.
lon : AscendingGrid
  The longitude grid.

Returns
-------
gridded_atm : AtmField
  An atmospheric field with all fields gridded to the input altitude, latitude, and longitude grids.
)");

  fld.def(
      "to_xarray",
      [](const AtmField &atm, const std::vector<Atm::Field::KeyVal> &keys) {
        auto xa = py::module_::import_("xarray");
        auto np = py::module_::import_("numpy");

        const auto xarray = Atm::Xarr(atm, keys);

        py::dict coords;
        coords["alt"] = np.attr("array")(xarray.altitudes);
        coords["lat"] = np.attr("array")(xarray.latitudes);
        coords["lon"] = np.attr("array")(xarray.longitudes);

        py::list dims;
        dims.append("alt");
        dims.append("lat");
        dims.append("lon");

        py::dict data;
        for (Size i = 0; i < xarray.keys.size(); i++) {
          const std::string key = std::format("{}", xarray.keys[i]);
          const Tensor3 d{xarray.data[i]};
          ARTS_USER_ERROR_IF(data.contains(key), "Duplicate key: {}", key)
          data[key.c_str()] = xa.attr("Variable")(dims, np.attr("array")(d));
        }

        py::dict attrs;
        attrs["top_of_atmosphere"] = xarray.toa;

        return xa.attr("Dataset")(data, coords, attrs);
      },
      "keys"_a = std::vector<Atm::Field::KeyVal>{},
      R"(Convert the atmospheric field to an xarray dataset.

The atmospheric field must be gridded using, e.g., the :func:`~pyarts3.arts.AtmField.as_gridded` method.

Parameters
----------
keys : list, optional
  The keys to include in the xarray dataset.  Default is empty, which includes all keys.

Returns
-------
dataset : xarray.Dataset
  The dataset with the atmospheric field data.  The coordinates are "alt", "lat", and "lon".
)");

  m.def(
      "stringify_keys",
      [](const std::unordered_map<Atm::KeyVal, Atm::FieldData> &dict,
         bool unique) {
        std::unordered_map<std::string, Atm::FieldData> out;

        for (auto &[key, vec] : dict) {
          const std::string strkey = std::format("{}", key);
          if (unique and out.contains(strkey)) {
            throw std::runtime_error(
                std::format("Key \"{}\" us not unique", strkey));
          }
          out[strkey] = vec;
        }

        return out;
      },
      "data"_a,
      "unique"_a = true,
      R"(Convert a dictionary of atmospheric field keys to their string representation.

.. warning::
  Setting `unique` to `False` will overwrite keys with the same string representation.

Parameters
----------
  data : dict
    A dictionary of atmospheric keys and corresponding data.
  unique : bool, optional
    If True, non-unique keys will throw an error.  Default is True.

Return
------
  py_atm_field : dict
    A dictionary of atmospheric keys and corresponding data.  The keys are the string representations of the atmospheric field keys.
)");

  fld.def(
         "__init__",
         [](AtmField *x, Numeric toa, IsoRatioOption iso) {
           new (x) AtmField(iso);
           x->top_of_atmosphere = toa;
         },
         "toa"_a = 100e3,
         "iso"_a = IsoRatioOption::Builtin)
      .def(
          "__getitem__",
          [](AtmField &self, const AtmKeyVal &key) -> Atm::Data & {
            if (self.contains(key)) return self[key];
            throw py::key_error(std::format("{}", key).c_str());
          },
          py::rv_policy::reference_internal)
      .def("__setitem__",
           [](AtmField &self, const AtmKeyVal &key, const Atm::Data &x) {
             self[key] = x;
           })
      .def("__contains__", [](const AtmField &self, const AtmKeyVal &x) {
        return self.contains(x);
      });

  fld.def_rw(
      "nlte",
      &AtmField::nlte,
      "NLTE data\n\n.. :class:`dict[~pyarts3.arts.QuantumIdentifier, ~pyarts3.arts.AtmData]`");

  fld.def_rw(
      "specs",
      &AtmField::specs,
      "Species data\n\n.. :class:`dict[~pyarts3.arts.SpeciesEnum, ~pyarts3.arts.AtmData]`");

  fld.def_rw(
      "isots",
      &AtmField::isots,
      "Isotopologue ratio data\n\n.. :class:`dict[~pyarts3.arts.SpeciesIsotope, ~pyarts3.arts.AtmData]`");

  fld.def_rw(
      "other",
      &AtmField::other,
      "Basic atmospheric data\n\n.. :class:`dict[~pyarts3.arts.AtmKey, ~pyarts3.arts.AtmData]`");

  fld.def_rw(
      "ssprops",
      &AtmField::ssprops,
      "Scattering species properties data\n\n.. :class:`dict[~pyarts3.arts.ScatteringSpeciesProperty, ~pyarts3.arts.AtmData]`");

  pnt.def(
      "to_dict",
      [](const AtmPoint &atm,
         bool core,
         bool specs,
         bool isots,
         bool nlte,
         bool ssprops) {
        std::unordered_map<Atm::KeyVal, Numeric> out;

        const auto keys = atm.keys();

        for (const Atm::KeyVal &key : keys) {
          if (std::holds_alternative<SpeciesEnum>(key)) {
            if (specs) out[key] = atm[key];
          } else if (std::holds_alternative<SpeciesIsotope>(key)) {
            if (isots) out[key] = atm[key];
          } else if (std::holds_alternative<QuantumLevelIdentifier>(key)) {
            if (nlte) out[key] = atm[key];
          } else if (std::holds_alternative<ScatteringSpeciesProperty>(key)) {
            if (ssprops) out[key] = atm[key];
          } else if (std::holds_alternative<AtmKey>(key)) {
            if (core) out[key] = atm[key];
          } else {
            throw std::runtime_error("Unknown key type.");
          }
        }

        return out;
      },
      "core"_a    = true,
      "specs"_a   = true,
      "isots"_a   = true,
      "nlte"_a    = true,
      "ssprops"_a = true,
      R"(Convert an atmospheric point to a dictionary.

.. tip::
  The function :func:`stringify_keys` may be used to convert the dictionary keys to strings.
  This allows for easier manipulation of the data of the dictionary.

Parameters
----------
  core : bool, optional
    If True, the core atmospheric data will be included (i.e., temperature, pressure, etc).  Default is True.
  specs : bool, optional
    If True, the species VMR data will be included.  Default is True.
  isots : bool, optional
    If True, the isotopologue ratio data will be included.  Default is True.
  nlte : bool, optional
    If True, the NLTE data will be included.  Default is True.
  ssprops : bool, optional
    If True, the scattering species properties data will be included.  Default is True.

Return
------
  dict : dict
    A dictionary of atmospheric keys and corresponding data.
)");

  const auto update_point =
      [](AtmPoint &atm, const std::unordered_map<Atm::KeyVal, Numeric> &dict) {
        for (auto &[key, data] : dict) {
          atm[key] = data;
        }
      };

  pnt.def_static(
      "from_dict",
      [&update_point](const std::unordered_map<Atm::KeyVal, Numeric> &dict,
                      const IsoRatioOption &iso) {
        AtmPoint out(iso);
        update_point(out, dict);
        return out;
      },
      "data"_a,
      "iso"_a = IsoRatioOption::Builtin,
      R"(Create an atmospheric point from a dictionary.

.. tip::
  Each key-type of the dictionary is constructible from a :class:`str`.  So using a dictionary with string keys is possible.

Parameters
----------
  data : dict
    A dictionary of atmospheric keys and corresponding data.
  iso : IsoRatioOption, optional
    The isotopologue ratio option to use.  Default is "Builtin". Use "None" to create point without isotopologue ratios.

Return
------
  atm_point : AtmPoint
    The atmospheric point created from the dictionary.
)");

  pnt.def("update",
          update_point,
          "data"_a,
          R"(Update the atmospheric point with dictionary values.

.. tip::
  Each key-type of the dictionary is constructible from a :class:`str`.  So using a dictionary with string keys is possible.

Parameters
----------
  data : dict
    A dictionary of atmospheric keys and corresponding data.
)");

  m.def(
      "stringify_keys",
      [](const std::unordered_map<Atm::KeyVal, Numeric> &dict, bool unique) {
        std::unordered_map<std::string, Numeric> out;

        for (auto &[key, vec] : dict) {
          const std::string strkey = std::format("{}", key);
          if (unique and out.contains(strkey)) {
            throw std::runtime_error(
                std::format("Key \"{}\" is not unique", strkey));
          }
          out[strkey] = vec;
        }

        return out;
      },
      "data"_a,
      "unique"_a = true,
      R"(Convert a dictionary of atmospheric point keys to their string representation.

.. warning::
  Setting `unique` to `False` will overwrite keys with the same string representation.

Parameters
----------
  data : dict
    A dictionary of atmospheric keys and corresponding data.
  unique : bool, optional
    If True, non-unique keys will throw an error.  Default is True.

Return
------
  dict : dict
    A dictionary of atmospheric keys and corresponding data.
)");

  auto aap =
      py::bind_vector<ArrayOfAtmPoint, py::rv_policy::reference_internal>(
          m, "ArrayOfAtmPoint");
  generic_interface(aap);
  aap.def(
      "field1D",
      [](const ArrayOfAtmPoint &atm,
         const AscendingGrid &altitudes,
         const InterpolationExtrapolation &extrap) {
        return Atm::atm_from_profile(atm, altitudes, extrap);
      },
      "altitudes"_a,
      "extrap"_a = InterpolationExtrapolation::Nearest,
      R"(Create an atmospheric field from a profile of atmospheric points.

The field is a 1D profile

Parameters
----------
  atm : ArrayOfAtmPoint
    The atmospheric profile.
  altitudes : AscendingGrid
    The altitude grid
  extrap : InterpolationExtrapolation, optional
    The extrapolation method to use for the new keys.  Default is "Nearest".

Return
------
  atm_field : AtmField
    The atmospheric field created from the profile.
)");

  aap.def("extend_in_pressure",
          &Atm::extend_in_pressure,
          R"(Extend the atmosphere to a new pressure point.

The logarithm of the pressure profile will be used as a pseudo-altitude coordinate
in calling :func:`~pyarts3.arts.ArrayOfAtmPoint.field1D`.  The extrapolation type will be used for the new points,
above and below the original profile (if necessary).

Parameters
----------
  extrapolation_type : InterpolationExtrapolation
    The extrapolation type to use for the new points.  Default is "Nearest".
  pressure : Numeric, optional
    The new pressure point.
  logarithmic : bool, optional
    The grid should be extended as in logarithmic pressure.  Default is True.
)",
          "pressure"_a,
          "extrapolation_type"_a = InterpolationExtrapolation::Nearest,
          "logarithmic"_a        = true);

  aap.def(
      "to_dict",
      [](const ArrayOfAtmPoint &atm,
         bool core,
         bool specs,
         bool isots,
         bool nlte,
         bool ssprops) {
        std::unordered_map<Atm::KeyVal, Vector> out;

        if (not atm.empty()) {
          const auto keys = atm.front().keys();

          for (const AtmPoint &point : atm) {
            if (point.keys() != keys) {
              throw std::runtime_error(std::format(
                  R"(Error for item of index {}.

All atmospheric points must have the same keys.  Indexed item differs from first item in keys.

Front item keys:      {:B,}
Indexed item keys:    {:B,}
)",
                  std::distance(&atm.front(), &point),
                  keys,
                  point.keys()));
            }

            for (const Atm::KeyVal &key : keys) {
              if (std::holds_alternative<SpeciesEnum>(key)) {
                if (specs) out[key].push_back(point[key]);
              } else if (std::holds_alternative<SpeciesIsotope>(key)) {
                if (isots) out[key].push_back(point[key]);
              } else if (std::holds_alternative<QuantumLevelIdentifier>(key)) {
                if (nlte) out[key].push_back(point[key]);
              } else if (std::holds_alternative<ScatteringSpeciesProperty>(
                             key)) {
                if (ssprops) out[key].push_back(point[key]);
              } else if (std::holds_alternative<AtmKey>(key)) {
                if (core) out[key].push_back(point[key]);
              } else {
                throw std::runtime_error("Unknown key type.");
              }
            }
          }
        }

        return out;
      },
      "core"_a    = true,
      "specs"_a   = true,
      "isots"_a   = true,
      "nlte"_a    = true,
      "ssprops"_a = true,
      R"(Convert an array of atmospheric points to a dictionary.

.. tip::
  The function :func:`stringify_keys` may be used to convert the dictionary keys to strings.
  This allows for easier manipulation of the data of the dictionary.

Parameters
----------
  core : bool, optional
    If True, the core atmospheric data will be included (i.e., temperature, pressure, etc).  Default is True.
  specs : bool, optional
    If True, the species VMR data will be included.  Default is True.
  isots : bool, optional
    If True, the isotopologue ratio data will be included.  Default is True.
  nlte : bool, optional
    If True, the NLTE data will be included.  Default is True.
  ssprops : bool, optional
    If True, the scattering species properties data will be included.  Default is True.

Return
------
  dict : dict
    A dictionary of atmospheric keys and corresponding data.  The keys are the atmospheric point keys and the values are the corresponding data.
)");

  const auto update_pointlist = [](ArrayOfAtmPoint &atm,
                                   const std::unordered_map<Atm::KeyVal, Vector>
                                       &dict) {
    if (dict.empty()) return;

    const Size n = atm.size();

    if (stdr::any_of(dict, [n](auto &v) { return v.second.size() != n; })) {
      throw std::runtime_error(
          "All values in the dictionary must have the same length to match the AtmPoint-list length.");
    }

    for (auto &[key, vec] : dict) {
      for (Size i = 0; i < n; i++) {
        atm[i][key] = vec[i];
      }
    }
  };

  aap.def_static(
      "from_dict",
      [&update_pointlist](const std::unordered_map<Atm::KeyVal, Vector> &dict,
                          const IsoRatioOption &iso) {
        ArrayOfAtmPoint out;
        if (dict.empty()) return out;
        out.resize(dict.begin()->second.size());
        for (auto &point : out) point = AtmPoint(iso);
        update_pointlist(out, dict);
        return out;
      },
      "data"_a,
      "iso"_a = IsoRatioOption::Builtin,
      R"(Create an array of atmospheric points from a dictionary.

.. tip::
  Each key-type of the dictionary is constructible from a :class:`str`.  So using a dictionary with string keys is possible.

Parameters
----------
  data : dict
    A dictionary with keys as the atmospheric point keys and values.
  iso : IsoRatioOption, optional
    The isotopologue ratio option to use.  Default is "Builtin". Use "None" to create points without isotopologue ratios.

Return
------
  atm_profile : AtmProfile
    The atmospheric profile created from the dictionary.
)");

  aap.def("update",
          update_pointlist,
          "data"_a,
          R"(Update the array of atmospheric points.

.. tip::
  Each key-type of the dictionary is constructible from a :class:`str`.  So using a dictionary with string keys is possible.

Parameters
----------
  data : dict
    A dictionary with keys as the atmospheric point keys and values.
)");

  m.def(
      "stringify_keys",
      [](const std::unordered_map<Atm::KeyVal, Vector> &dict, bool unique) {
        std::unordered_map<std::string, Vector> out;

        for (auto &[key, vec] : dict) {
          const std::string strkey = std::format("{}", key);
          if (unique and out.contains(strkey)) {
            throw std::runtime_error(
                std::format("Key \"{}\" is not unique", strkey));
          }
          out[strkey] = vec;
        }

        return out;
      },
      "data"_a,
      "unique"_a = true,
      R"(Convert a dictionary of atmospheric point keys to their string representation.

.. warning::
  Setting `unique` to `False` will overwrite keys with the same string representation.

Parameters
----------
  data : dict
    A dictionary with keys as the atmospheric point keys and values.
  unique : bool, optional
    If True, non-unique keys will throw an error.  Default is True.

Return
------
  dict : dict
    A dictionary of atmospheric keys and corresponding data.  The keys are the atmospheric point keys and the values are the corresponding data.
)");
} catch (std::exception &e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize atm\n{}", e.what()));
}
}  // namespace Python
