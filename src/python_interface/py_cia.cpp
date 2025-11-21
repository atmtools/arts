#include <nanobind/stl/array.h>
#include <nanobind/stl/bind_map.h>
#include <python_interface.h>

#include <algorithm>
#include <memory>
#include <string>
#include <utility>

#include "cia.h"
#include "debug.h"
#include "hpy_arts.h"
#include "isotopologues.h"
#include "matpack_mdspan_helpers_grid_t.h"
#include "matpack_mdspan_helpers_gridded_data_t.h"
#include "physics_funcs.h"
#include "species_tags.h"

namespace Python {
void py_cia(py::module_& m) try {
  py::class_<SpeciesEnumPair> sep(m, "SpeciesEnumPair");
  generic_interface(sep);
  sep.def(py::init<SpeciesEnum, SpeciesEnum>(),
          "spec1"_a,
          "spec2"_a,
          "Create a species pair from two species")
      .def(
          "__init__",
          [](SpeciesEnumPair* self, const std::string& s) {
            auto dash_pos = s.find('-');
            if (dash_pos == std::string::npos) {
              throw std::runtime_error(
                  std::format("Invalid SpeciesEnumPair string format: '{}'. "
                              "Expected 'SpeciesA-SpeciesB'",
                              s));
            }
            std::string spec1_str = s.substr(0, dash_pos);
            std::string spec2_str = s.substr(dash_pos + 1);
            SpeciesEnum spec1     = to<SpeciesEnum>(spec1_str);
            SpeciesEnum spec2     = to<SpeciesEnum>(spec2_str);
            new (self) SpeciesEnumPair{.spec1 = spec1, .spec2 = spec2};
          },
          "s"_a,
          "Create a species pair from string 'SpeciesA-SpeciesB'")
      .def_rw("spec1", &SpeciesEnumPair::spec1, "First species")
      .def_rw("spec2", &SpeciesEnumPair::spec2, "Second species")
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def("__hash__",
           [](const SpeciesEnumPair& self) {
             return std::hash<SpeciesEnumPair>{}(self);
           })
      .def("__getstate__",
           [](const SpeciesEnumPair& self) {
             return std::make_tuple(self.spec1, self.spec2);
           })
      .def("__setstate__",
           [](SpeciesEnumPair* self,
              const std::tuple<SpeciesEnum, SpeciesEnum>& state) {
             new (self) SpeciesEnumPair(std::get<0>(state), std::get<1>(state));
           });
  generic_interface(sep);

  // Make SpeciesEnumPair implicitly convertible from string
  py::implicitly_convertible<std::string, SpeciesEnumPair>();

  py::class_<CIARecord> cia(m, "CIARecord");
  generic_interface(cia);
  cia.def(py::init<ArrayOfGriddedField2>())
      .def_prop_rw(
          "data",
          [](const CIARecord& c) { return c.Data(); },
          [](CIARecord& c, const ArrayOfGriddedField2& data) {
            c.Data() = data;
          },
          "Data by bands\n\n.. :class:`~pyarts3.arts.ArrayOfGriddedField2`")
      .def(
          "spectral_propmat",
          [](const CIARecord& self,
             const AscendingGrid& f,
             const AtmPoint& atm,
             const SpeciesEnum spec1,
             const SpeciesEnum spec2,
             const Numeric T_extrapolfac,
             const Index robust) {
            PropmatVector out(f.size());

            const Numeric scl =
                atm.number_density(spec1) * atm.number_density(spec2);
            for (auto& cia_data : self.Data()) {
              Vector result(f.size(), 0);
              cia_interpolation(
                  result, f, atm.temperature, cia_data, T_extrapolfac, robust);

              for (Size i = 0; i < f.size(); i++) {
                out[i].A() += scl * result[i];
              }
            }

            return out;
          },
          "f"_a,
          "atm"_a,
          "spec1"_a         = "AIR"_spec,
          "spec2"_a         = "AIR"_spec,
          "T_extrapolfac"_a = 0.0,
          "robust"_a        = 1,
          R"--(Computes the collision-induced absorption in 1/m

Parameters
----------
f : AscendingGrid
    Frequency grid [Hz]
atm : AtmPoint
    Atmospheric point
spec1 : SpeciesEnum
    First species
spec2 : SpeciesEnum
    Second species
T_extrapolfac : Numeric, optional
    Extrapolation in temperature.  The default is 0
robust : Index, optional
    Returns NaN instead of throwing if it evaluates true.  The default is 1

Returns
-------
  abs : PropmatVector
    Absorption profile [1/m]

)--")
      .def(
          "compute_abs",
          [](CIARecord& cia_,
             Numeric T,
             Numeric P,
             Numeric X0,
             Numeric X1,
             const Vector& f,
             Numeric T_extrapolfac,
             Index robust) {
            Vector out(f.size(), 0);

            for (auto& cia_data : cia_.Data()) {
              Vector result(f.size(), 0);
              cia_interpolation(result, f, T, cia_data, T_extrapolfac, robust);
              out += result;
            }

            out *= Math::pow2(number_density(P, T)) * X0 * X1;
            return out;
          },
          "T"_a,
          "P"_a,
          "X0"_a,
          "X1"_a,
          "f"_a,
          "T_extrapolfac"_a = 0.0,
          "robust"_a        = 1,
          R"--(Computes the collision-induced absorption in 1/m

Parameters
----------
T : Numeric
    Temperature [K]
P : Numeric
    Pressure [Pa]
X0 : Numeric
    VMR of species 1 [-]
X1 : Numeric
    VMR of species 2 [-]
f : Vector
    Frequency grid [Hz]
T_extrapolfac : Numeric, optional
    Extrapolation in temperature.  The default is 0
robust : Index, optional
    Returns NaN instead of throwing if it evaluates true.  The default is 1

Returns
-------
  abs : Vector
    Absorption profile [1/m]

)--")
      .def("__getstate__", [](const CIARecord& t) { return t.Data(); })
      .def("__setstate__", [](CIARecord* c, const ArrayOfGriddedField2& state) {
        new (c) CIARecord();
        c->Data() = state;
      });

  // Bind CIARecords as map
  auto acr = py::bind_map<CIARecords, py::rv_policy::reference_internal>(
      m, "CIARecords");
  generic_interface(acr);

  acr.def(
      "spectral_propmat",
      [](const CIARecords& self,
         const AscendingGrid& f,
         const AtmPoint& atm,
         const SpeciesEnum& spec,
         const Numeric T_extrapolfac,
         const Index ignore_errors,
         const py::kwargs&) {
        PropmatVector spectral_propmat(f.size());
        PropmatMatrix spectral_propmat_jac(0, f.size());
        JacobianTargets jac_targets{};

        spectral_propmatAddCIA(spectral_propmat,
                               spectral_propmat_jac,
                               spec,
                               jac_targets,
                               f,
                               atm,
                               self,
                               T_extrapolfac,
                               ignore_errors);

        return spectral_propmat;
      },
      "f"_a,
      "atm"_a,
      "spec"_a          = SpeciesEnum::Bath,
      "T_extrapolfac"_a = Numeric{0.0},
      "ignore_errors"_a = Index{1},
      "kwargs"_a        = py::kwargs{},
      R"--(Computes the collision-induced absorption in 1/m

Parameters
----------
f : AscendingGrid
    Frequency grid [Hz]
atm : AtmPoint
    Atmospheric point
spec : SpeciesEnum, optional
    Species to use.  Defaults to all.
T_extrapolfac : Numeric, optional
    Extrapolation in temperature.  The default is 0
ignore_errors : Index, optional
    Returns NaN instead of throwing if it evaluates true.  The default is 1

Returns
-------
spectral_propmat : PropmatVector
    Propagation matrix by frequency [1/m]

)--");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize cia\n{}", e.what()));
}
}  // namespace Python
