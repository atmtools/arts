
#include <atm.h>
#include <atm_path.h>
#include <debug.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>
#include <quantum_numbers.h>
#include <species_tags.h>

#include <stdexcept>
#include <unordered_map>

#include "enumsInterpolationExtrapolation.h"
#include "enumsIsoRatioOption.h"
#include "enumsSpeciesEnum.h"
#include "hpy_arts.h"
#include "isotopologues.h"
#include "physics_funcs.h"

namespace Python {
void py_atm(py::module_ &m) try {
  py::class_<Atm::Data> atmdata(m, "AtmData");
  atmdata.def(py::init<Atm::Data>())
      .def(py::init<GriddedField3>())
      .def(py::init<Numeric>())
      .def(py::init<Atm::FunctionalData>())
      .def_rw("data", &Atm::Data::data, "The data")
      .def_rw("alt_upp", &Atm::Data::alt_upp, "Upper altitude limit")
      .def_rw("alt_low", &Atm::Data::alt_low, "Lower altitude limit")
      .def_rw("lat_upp", &Atm::Data::lat_upp, "Upper latitude limit")
      .def_rw("lat_low", &Atm::Data::lat_low, "Lower latitude limit")
      .def_rw("lon_upp", &Atm::Data::lon_upp, "Upper longitude limit")
      .def_rw("lon_low", &Atm::Data::lon_low, "Lower longitude limit")
      .def(
          "at",
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
      .def_prop_ro("data_type", &Atm::Data::data_type, "The data type")
      .def(
          "__call__",
          [](const Atm::Data &d, Numeric alt, Numeric lat, Numeric lon) {
            return d.at(alt, lat, lon);
          },
          "alt"_a,
          "lat"_a,
          "lon"_a,
          "Get a point of data at the position")
      .def("__getstate__",
           [](const Atm::Data &t) {
             return std::make_tuple(t.data,
                                    t.alt_low,
                                    t.alt_upp,
                                    t.lat_low,
                                    t.lat_upp,
                                    t.lon_low,
                                    t.lon_upp);
           })
      .def("__setstate__",
           [](Atm::Data *a,
              const std::tuple<Atm::FieldData,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation> &state) {
             new (a) Atm::Data();
             a->data    = std::get<0>(state);
             a->alt_low = std::get<1>(state);
             a->alt_upp = std::get<2>(state);
             a->lat_low = std::get<3>(state);
             a->lat_upp = std::get<4>(state);
             a->lon_low = std::get<5>(state);
             a->lon_upp = std::get<6>(state);
           })
      .doc() = "Atmospheric data";
  str_interface(atmdata);
  py::implicitly_convertible<GriddedField3, Atm::Data>();
  py::implicitly_convertible<Numeric, Atm::Data>();
  py::implicitly_convertible<Index, Atm::Data>();
  py::implicitly_convertible<Atm::FunctionalData::func_t, Atm::Data>();

  auto pnt = py::class_<AtmPoint>(m, "AtmPoint");
  workspace_group_interface(pnt);

  auto fld = py::class_<AtmField>(m, "AtmField");
  workspace_group_interface(fld);

  pnt.def_rw("temperature", &AtmPoint::temperature, "Temperature [K]")
      .def_rw("pressure", &AtmPoint::pressure, "Pressure [Pa]")
      .def_rw("mag", &AtmPoint::mag, "Magnetic field [T]")
      .def_rw("wind", &AtmPoint::wind, "Wind field [m/s]")
      .def(
          "number_density",
          [](AtmPoint &atm, std::variant<SpeciesEnum, SpeciesIsotope> key) {
            return std::visit(
                [&atm]<typename T>(const T &s) {
                  if (not atm.has(s))
                    throw py::key_error(var_string(s).c_str());

                  Numeric nd =
                      atm[s] * number_density(atm.pressure, atm.temperature);
                  if constexpr (requires { s.spec; }) {
                    nd *= atm[s.spec];
                  }
                  return nd;
                },
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
)")
      .def(
          "species_vmr",
          [](AtmPoint &atm, SpeciesEnum s) {
            if (not atm.has(s)) throw py::key_error(var_string(s).c_str());
            return atm[s];
          },
          "spec"_a,
          "Get the VMR of the species")
      .def(
          "set_species_vmr",
          [](AtmPoint &atm, SpeciesEnum s, Numeric x) { atm[s] = x; },
          "spec"_a,
          "x"_a,
          "Set the VMR of the species")
      .def(
          "isotopologue_ratio",
          [](AtmPoint &atm, SpeciesIsotope s) {
            if (not atm.has(s)) throw py::key_error(var_string(s).c_str());
            return atm[s];
          },
          "isot"_a,
          "Get the isotopologue ratio")
      .def(
          "set_isotopologue_ratio",
          [](AtmPoint &atm, SpeciesIsotope s, Numeric x) { atm[s] = x; },
          "isot"_a,
          "x"_a,
          "Set the isotopologue ratio")
      .def(
          "nlte_value",
          [](AtmPoint &atm, const QuantumIdentifier &s) {
            if (not atm.has(s)) throw py::key_error(var_string(s).c_str());
            return atm[s];
          },
          "qid"_a,
          "Get the NLTE value")
      .def(
          "set_nlte_value",
          [](AtmPoint &atm, const QuantumIdentifier &s, Numeric x) {
            atm[s] = x;
          },
          "qid"_a,
          "x"_a,
          "Set the NLTE value")
      .def(
          "__getitem__",
          [](AtmPoint &atm, const Atm::KeyVal &x) {
            return std::visit(
                [&atm](auto &key) {
                  if (not atm.has(key))
                    throw py::key_error(var_string(key).c_str());
                  return atm[key];
                },
                x);
          },
          "key"_a,
          "Get the value of a key")
      .def(
          "__getitem__",
          [](AtmPoint &atm, const ArrayOfSpeciesTag &x) {
            if (x.empty()) throw py::key_error("empty species tag list");
            const SpeciesEnum key = x.Species();
            if (not atm.has(key)) throw py::key_error(var_string(key).c_str());
            return atm[key];
          },
          "key"_a,
          "Set the value of a key")
      .def(
          "__setitem__",
          [](AtmPoint &atm, const Atm::KeyVal &x, Numeric data) {
            atm[x] = data;
          },
          "key"_a,
          "val"_a,
          "Set the value of a key")
      .def(
          "__setitem__",
          [](AtmPoint &atm, const ArrayOfSpeciesTag &x, Numeric data) {
            if (x.empty()) throw py::key_error("empty species tag list");
            atm[x.Species()] = data;
          },
          "key"_a,
          "val"_a,
          "Set the value of a key")
      .def("keys", &AtmPoint::keys, "Available keys")
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
          "Returns a flat list of values.")
      .def("__getstate__",
           [](const AtmPoint &t) {
             auto k = t.keys();
             std::vector<Numeric> v;
             v.reserve(k.size());
             for (auto &kn : k)
               v.emplace_back(
                   std::visit([&](auto &&key) { return t[key]; }, kn));
             return std::make_tuple(k, v);
           })
      .def("__setstate__",
           [](AtmPoint *a,
              const std::tuple<std::vector<AtmKeyVal>, std::vector<Numeric>>
                  &state) {
             auto k = std::get<0>(state);
             auto v = std::get<1>(state);
             ARTS_USER_ERROR_IF(k.size() != v.size(), "Invalid state!")

             new (a) AtmPoint();
             for (std::size_t i = 0; i < k.size(); i++)
               std::visit(
                   [&](auto &&key) -> Numeric & { return a->operator[](key); },
                   k[i]) = v[i];
           });

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
          [](AtmField &atm, const Atm::KeyVal &x) -> Atm::Data & {
            std::visit(
                [&atm](auto &key) {
                  if (not atm.has(key))
                    throw py::key_error(var_string(key).c_str());
                },
                x);

            return atm[x];
          },
          py::rv_policy::reference_internal,
          "key"_a,
          "Get the value of a key")
      .def(
          "__getitem__",
          [](AtmField &atm, const ArrayOfSpeciesTag &x) -> Atm::Data & {
            if (x.empty()) throw py::key_error("empty species tag list");
            const SpeciesEnum key = x.Species();
            if (not atm.has(key)) throw py::key_error(var_string(key).c_str());
            return atm[key];
          },
          py::rv_policy::reference_internal,
          "key"_a,
          "Get the value of a key")
      .def(
          "__setitem__",
          [](AtmField &atm, const Atm::KeyVal &x, const Atm::FieldData &data) {
            atm[x].data = data;
          },
          "key"_a,
          "val"_a,
          "Set the value of a key")
      .def(
          "__setitem__",
          [](AtmField &atm,
             const ArrayOfSpeciesTag &x,
             const Atm::FieldData &data) {
            if (x.empty()) throw py::key_error("empty species tag list");
            atm[x.Species()].data = data;
          },
          "key"_a,
          "val"_a,
          "Set the value of a key")
      .def(
          "at",
          [](const AtmField &atm,
             const Numeric &h,
             const Numeric &lat,
             const Numeric &lon) { return atm.at(h, lat, lon); },
          "h"_a,
          "lat"_a,
          "lon"_a,
          "Get the data at a point")
      .def_rw("top_of_atmosphere",
              &AtmField::top_of_atmosphere,
              "Top of the atmosphere [m]")
      .def("__getstate__",
           [](const AtmField &t) {
             auto k = t.keys();
             std::vector<Atm::Data> v;
             v.reserve(k.size());
             for (auto &kn : k)
               v.emplace_back(
                   std::visit([&](auto &&key) { return t[key]; }, kn));
             return std::make_tuple(k, v, t.top_of_atmosphere);
           })
      .def(
          "__setstate__",
          [](AtmField *a,
             const std::tuple<std::vector<AtmKeyVal>,
                              std::vector<Atm::Data>,
                              Numeric> &state) {
            auto k = std::get<0>(state);
            auto v = std::get<1>(state);
            ARTS_USER_ERROR_IF(k.size() != v.size(), "Invalid state!")

            new (a) AtmField();
            a->top_of_atmosphere = std::get<2>(state);

            for (std::size_t i = 0; i < k.size(); i++)
              std::visit(
                  [&](auto &&key) -> Atm::Data & { return a->operator[](key); },
                  k[i]) = v[i];
          });

  m.def(
      "frequency_shift",
      [](AscendingGrid frequency_grid,
         const AtmPoint &atmospheric_point,
         const PropagationPathPoint &ray_path_point) {
        Vector3 x;
        frequency_gridWindShift(
            frequency_grid, x, atmospheric_point, ray_path_point);
        return frequency_grid;
      },
      "frequency_grid"_a,
      "atmospheric_point"_a,
      "ray_path_point"_a,
      R"(Get the frequency-shifted frequency grid at a point in the atmosphere.

Parameters
----------
  frequency_grid : AscendingGrid
    The frequency grid to shift.
  atmospheric_point : AtmPoint
    The point in the atmosphere.
  ray_path_point : PropagationPathPoint
    The point along the ray path.
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
          } else if (std::holds_alternative<QuantumIdentifier>(key)) {
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
)");

  const auto update_field =
      [](AtmField &atm,
         const std::unordered_map<Atm::KeyVal, Atm::FieldData> &dict,
         const InterpolationExtrapolation &extrap) {
        for (auto &[key, data] : dict) {
          const bool existed =
              std::visit([&atm](auto &x) { return atm.has(x); }, key);
          atm[key].data = data;
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

  m.def(
      "stringify_keys",
      [](const std::unordered_map<Atm::KeyVal, Atm::FieldData> &dict,
         bool unique) {
        std::unordered_map<std::string, Atm::FieldData> out;

        for (auto &[key, vec] : dict) {
          const std::string strkey = var_string(key);
          if (unique and out.contains(strkey)) {
            throw std::runtime_error(
                var_string("Key \"", strkey, "\" is not unique."));
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
          [](AtmField &atm, const Atm::KeyVal &x) -> Atm::Data & {
            std::visit(
                [&atm](auto &key) {
                  if (not atm.has(key))
                    throw py::key_error(var_string(key).c_str());
                },
                x);

            return atm[x];
          },
          py::rv_policy::reference_internal,
          "Get the value of a key")
      .def(
          "__getitem__",
          [](AtmField &atm, const ArrayOfSpeciesTag &x) -> Atm::Data & {
            if (x.empty()) throw py::key_error("empty species tag list");
            const SpeciesEnum key = x.Species();
            if (not atm.has(key)) throw py::key_error(var_string(key).c_str());
            return atm[key];
          },
          py::rv_policy::reference_internal,
          "Get the value of a key")
      .def(
          "__setitem__",
          [](AtmField &atm, const Atm::KeyVal &x, const Atm::FieldData &data) {
            atm[x].data = data;
          },
          "Set the value of a key")
      .def(
          "__setitem__",
          [](AtmField &atm,
             const ArrayOfSpeciesTag &x,
             const Atm::FieldData &data) {
            if (x.empty()) throw py::key_error("empty species tag list");
            atm[x.Species()].data = data;
          },
          "Set the value of a key")
      .def(
          "at",
          [](const AtmField &atm,
             const Numeric &h,
             const Numeric &lat,
             const Numeric &lon) { return atm.at(h, lat, lon); },
          "h"_a,
          "lat"_a,
          "lon"_a,
          "Get the data at a point")
      .def_rw("top_of_atmosphere",
              &AtmField::top_of_atmosphere,
              "Top of the atmosphere [m]")
      .def("__getstate__",
           [](const AtmField &t) {
             auto k = t.keys();
             std::vector<Atm::Data> v;
             v.reserve(k.size());
             for (auto &kn : k)
               v.emplace_back(
                   std::visit([&](auto &&key) { return t[key]; }, kn));
             return std::make_tuple(k, v, t.top_of_atmosphere);
           })
      .def(
          "__setstate__",
          [](AtmField *a,
             const std::tuple<std::vector<AtmKeyVal>,
                              std::vector<Atm::Data>,
                              Numeric> &state) {
            auto k = std::get<0>(state);
            auto v = std::get<1>(state);
            ARTS_USER_ERROR_IF(k.size() != v.size(), "Invalid state!")

            new (a) AtmField();
            a->top_of_atmosphere = std::get<2>(state);

            for (std::size_t i = 0; i < k.size(); i++)
              std::visit(
                  [&](auto &&key) -> Atm::Data & { return a->operator[](key); },
                  k[i]) = v[i];
          });

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
          } else if (std::holds_alternative<QuantumIdentifier>(key)) {
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
          const std::string strkey = var_string(key);
          if (unique and out.contains(strkey)) {
            throw std::runtime_error(
                var_string("Key \"", strkey, "\" is not unique."));
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
)");

  auto aap =
      py::bind_vector<ArrayOfAtmPoint, py::rv_policy::reference_internal>(
          m, "ArrayOfAtmPoint");
  workspace_group_interface(aap);

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
              } else if (std::holds_alternative<QuantumIdentifier>(key)) {
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
)");

  const auto update_pointlist = [](ArrayOfAtmPoint &atm,
                                   const std::unordered_map<Atm::KeyVal, Vector>
                                       &dict) {
    if (dict.empty()) return;

    const Size n = atm.size();

    if (std::ranges::any_of(dict, [N = static_cast<Index>(n)](auto &v) {
          return v.second.size() != N;
        })) {
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
          const std::string strkey = var_string(key);
          if (unique and out.contains(strkey)) {
            throw std::runtime_error(
                var_string("Key \"", strkey, "\" is not unique."));
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
)");
} catch (std::exception &e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize atm\n", e.what()));
}
}  // namespace Python
