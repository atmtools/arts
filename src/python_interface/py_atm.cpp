
#include <atm.h>
#include <atm_path.h>
#include <debug.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>
#include <quantum_numbers.h>
#include <species_tags.h>

#include "enums.h"
#include "hpy_arts.h"
#include "isotopologues.h"
#include "physics_funcs.h"

namespace Python {

void py_atm(py::module_ &m) try {
  py::class_<Atm::Data>(m, "AtmData")
      .def(py::init<Atm::Data>())
      .def(py::init<GriddedField3>())
      .def(py::init<Numeric>())
      .def(py::init<Atm::FunctionalData>())
      .def_rw("data", &Atm::Data::data)
      .def_rw("alt_upp", &Atm::Data::alt_upp)
      .def_rw("alt_low", &Atm::Data::alt_low)
      .def_rw("lat_upp", &Atm::Data::lat_upp)
      .def_rw("lat_low", &Atm::Data::lat_low)
      .def_rw("lon_upp", &Atm::Data::lon_upp)
      .def_rw("lon_low", &Atm::Data::lon_low)
      .def_prop_ro("data_type", &Atm::Data::data_type)
      .def("__call__",
           [](const Atm::Data &d, Numeric alt, Numeric lat, Numeric lon) {
             return d.at(alt, lat, lon);
           })
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
           });
  py::implicitly_convertible<GriddedField3, Atm::Data>();
  py::implicitly_convertible<Numeric, Atm::Data>();
  py::implicitly_convertible<Index, Atm::Data>();
  py::implicitly_convertible<Atm::FunctionalData, Atm::Data>();

  auto pnt = py::class_<AtmPoint>(m, "AtmPoint");
  workspace_group_interface(pnt);

  auto fld = py::class_<AtmField>(m, "AtmField");
  workspace_group_interface(fld);

  pnt.def_rw("temperature", &AtmPoint::temperature)
      .def_rw("pressure", &AtmPoint::pressure)
      .def_rw("mag", &AtmPoint::mag)
      .def_rw("wind", &AtmPoint::wind)
      .def("number_density",
           [](AtmPoint &atm, SpeciesIsotope s) {
             if (not atm.has(s)) throw py::key_error(var_string(s).c_str());
             if (not atm.has(s.spec))
               throw py::key_error(var_string(s).c_str());
             return atm[s] * atm[s.spec] *
                    number_density(atm.pressure, atm.temperature);
           })
      .def("species_vmr",
           [](AtmPoint &atm, SpeciesEnum s) {
             if (not atm.has(s)) throw py::key_error(var_string(s).c_str());
             return atm[s];
           })
      .def("set_species_vmr",
           [](AtmPoint &atm, SpeciesEnum s, Numeric x) { atm[s] = x; })
      .def("isotopologue_ratio",
           [](AtmPoint &atm, SpeciesIsotope s) {
             if (not atm.has(s)) throw py::key_error(var_string(s).c_str());
             return atm[s];
           })
      .def("set_isotopologue_ratio",
           [](AtmPoint &atm, SpeciesIsotope s, Numeric x) { atm[s] = x; })
      .def("nlte_value",
           [](AtmPoint &atm, const QuantumIdentifier &s) {
             if (not atm.has(s)) throw py::key_error(var_string(s).c_str());
             return atm[s];
           })
      .def("set_nlte_value",
           [](AtmPoint &atm, const QuantumIdentifier &s, Numeric x) {
             atm[s] = x;
           })
      .def(
          "__getitem__",
          [](AtmPoint &atm, AtmKey x) {
            if (not atm.has(x)) throw py::key_error(var_string(x).c_str());
            return atm[x];
          },
          py::rv_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmPoint &atm, const QuantumIdentifier &x) {
            if (not atm.has(x)) throw py::key_error(var_string(x).c_str());
            return atm[x];
          },
          py::rv_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmPoint &atm, const SpeciesEnum &x) {
            if (not atm.has(x)) throw py::key_error(var_string(x).c_str());
            return atm[x];
          },
          py::rv_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmPoint &atm, const SpeciesIsotope &x) {
            if (not atm.has(x)) throw py::key_error(var_string(x).c_str());
            return atm[x];
          },
          py::rv_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmPoint &atm, const ArrayOfSpeciesTag &x) {
            if (not atm.has(x.Species()))
              throw py::key_error(var_string(x).c_str());
            return atm[x.Species()];
          },
          py::rv_policy::reference_internal)
      .def("__setitem__",
           [](AtmPoint &atm, AtmKey x, Numeric data) { atm[x] = data; })
      .def("__setitem__",
           [](AtmPoint &atm, const QuantumIdentifier &x, Numeric data) {
             atm[x] = data;
           })
      .def("__setitem__",
           [](AtmPoint &atm, const SpeciesEnum &x, Numeric data) {
             atm[x] = data;
           })
      .def("__setitem__",
           [](AtmPoint &atm, const SpeciesIsotope &x, Numeric data) {
             atm[x] = data;
           })
      .def("__setitem__",
           [](AtmPoint &atm, const ArrayOfSpeciesTag &x, Numeric data) {
             atm[x.Species()] = data;
           })
      .def("keys", &AtmPoint::keys)
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
         "__getitem__",
         [](AtmField &atm, AtmKey x) -> Atm::Data & {
           if (not atm.has(x)) throw py::key_error(var_string(x).c_str());
           return atm[x];
         },
         py::rv_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmField &atm, const QuantumIdentifier &x) -> Atm::Data & {
            if (not atm.has(x)) throw py::key_error(var_string(x).c_str());
            return atm[x];
          },
          py::rv_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmField &atm, const SpeciesEnum &x) -> Atm::Data & {
            if (not atm.has(x)) throw py::key_error(var_string(x).c_str());
            return atm[x];
          },
          py::rv_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmField &atm, const SpeciesIsotope &x) -> Atm::Data & {
            if (not atm.has(x)) throw py::key_error(var_string(x).c_str());
            return atm[x];
          },
          py::rv_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmField &atm, const ArrayOfSpeciesTag &x) -> Atm::Data & {
            if (not atm.has(x.Species()))
              throw py::key_error(var_string(x.Species()).c_str());
            return atm[x.Species()];
          },
          py::rv_policy::reference_internal)
      .def(
          "__setitem__",
          [](AtmField &atm, AtmKey x, const Atm::Data &data) { atm[x] = data; })
      .def(
          "__setitem__",
          [](AtmField &atm, const QuantumIdentifier &x, const Atm::Data &data) {
            atm[x] = data;
          })
      .def("__setitem__",
           [](AtmField &atm, const SpeciesEnum &x, const Atm::Data &data) {
             atm[x] = data;
           })
      .def(
          "__setitem__",
          [](AtmField &atm, const ArrayOfSpeciesTag &x, const Atm::Data &data) {
            atm[x.Species()] = data;
          })
      .def("__setitem__",
           [](AtmField &atm, const SpeciesIsotope &x, const Atm::Data &data) {
             atm[x] = data;
           })
      .def("at",
           [](const AtmField &atm,
              const Vector &h,
              const Vector &lat,
              const Vector &lon) { return atm.at(h, lat, lon); })
      .def("at",
           [](const AtmField &atm,
              const Numeric &h,
              const Numeric &lat,
              const Numeric &lon) { return atm.at(h, lat, lon); })
      .def_rw("top_of_atmosphere", &AtmField::top_of_atmosphere)
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

  m.def("frequency_shift",
        [](const AscendingGrid &frequency_grid,
           const PropagationPathPoint &ray_path_point,
           const AtmPoint &atmospheric_point,
           const Numeric &rte_alonglos_v) {
          AscendingGrid out = frequency_grid;
          forward_path_freq(out,
                            frequency_grid,
                            ray_path_point,
                            atmospheric_point,
                            rte_alonglos_v);
          return out;
        });

  auto aap = py::bind_vector<ArrayOfAtmPoint>(m, "ArrayOfAtmPoint");
  workspace_group_interface(aap);
} catch (std::exception &e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize atm\n", e.what()));
}
}  // namespace Python
