#include <nanobind/stl/array.h>
#include <nanobind/stl/bind_vector.h>
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
  py::class_<CIARecord> cia(m, "CIARecord");
  generic_interface(cia);
  cia.def(py::init<ArrayOfGriddedField2, SpeciesEnum, SpeciesEnum>())
      .def_prop_ro(
          "specs",
          [](const CIARecord& c) { return c.TwoSpecies(); },
          "The two species\n\n.. :class:`tuple[pyarts3.arts.SpeciesEnum, pyarts3.arts.SpeciesEnum]`")
      .def_prop_rw(
          "data",
          [](const CIARecord& c) { return c.Data(); },
          [](CIARecord& c, const ArrayOfGriddedField2& data) {
            c.Data() = data;
          },
          "Data by bands\n\n.. :class:`~pyarts3.arts.ArrayOfGriddedField2`")
      .def(
          "propagation_matrix",
          [](const CIARecord& self,
             const AscendingGrid& f,
             const AtmPoint& atm,
             const Numeric T_extrapolfac,
             const Index robust) {
            PropmatVector out(f.size());

            const Numeric scl = atm.number_density(self.Species(0)) *
                                atm.number_density(self.Species(1));
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
          "T_extrapolfac"_a = 0.0,
          "robust"_a        = 1,
          R"--(Computes the collision-induced absorption in 1/m

Parameters
----------
f : AscendingGrid
    Frequency grid [Hz]
atm : AtmPoint
    Atmospheric point
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
      .def("__getstate__",
           [](const CIARecord& t) {
             return std::make_tuple(t.Data(), t.TwoSpecies());
           })
      .def(
          "__setstate__",
          [](CIARecord* c,
             const std::tuple<ArrayOfGriddedField2, std::array<SpeciesEnum, 2>>&
                 state) {
            new (c) CIARecord();
            c->Data()       = std::get<0>(state);
            c->TwoSpecies() = std::get<1>(state);
          });

  auto acr =
      py::bind_vector<ArrayOfCIARecord, py::rv_policy::reference_internal>(
          m, "ArrayOfCIARecord");
  generic_interface(acr);
  vector_interface(acr);

  acr.def(
      "propagation_matrix",
      [](const ArrayOfCIARecord& self,
         const AscendingGrid& f,
         const AtmPoint& atm,
         const SpeciesEnum& spec,
         const Numeric T_extrapolfac,
         const Index ignore_errors,
         const py::kwargs&) {
        PropmatVector propagation_matrix(f.size());
        PropmatMatrix propagation_matrix_jacobian(0, f.size());
        JacobianTargets jac_targets{};

        propagation_matrixAddCIA(propagation_matrix,
                                 propagation_matrix_jacobian,
                                 spec,
                                 jac_targets,
                                 f,
                                 atm,
                                 self,
                                 T_extrapolfac,
                                 ignore_errors);

        return propagation_matrix;
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
propagation_matrix : PropmatVector
    Propagation matrix by frequency [1/m]

)--");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize cia\n{}", e.what()));
}
}  // namespace Python
