#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <obsel.h>
#include <python_interface.h>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"
#include "sorted_grid.h"

namespace Python {
void py_sensor(py::module_& m) try {
  py::class_<SensorPosLos> splos(m, "SensorPosLos");
  workspace_group_interface(splos);
  common_ndarray(splos);
  splos
      .def(
          "__array__",
          [](SensorPosLos& x, py::object dtype, py::object copy) {
            std::array<size_t, 1> shape = {5};
            auto np                     = py::module_::import_("numpy");
            auto w =
                py::ndarray<py::numpy, Numeric, py::shape<5>, py::c_contig>(
                    &x, 1, shape.data(), py::cast(x));
            return np.attr("asarray")(w, "dtype"_a = dtype, "copy"_a = copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none(),
          "Returns a :class:`~numpy.ndarray` of the object.")
      .def_prop_rw(
          "value",
          [](py::object& x) { return x.attr("__array__")("copy"_a = false); },
          [](SensorPosLos& a, const SensorPosLos& b) { a = b; },
          "A :class:`~numpy.ndarray` of the object.")
      .def(py::init<Vector3, Vector2>(), "From pos and los")
      .def_rw("pos", &SensorPosLos::pos, "Position")
      .def_rw("los", &SensorPosLos::los, "Line of sight");

  py::class_<SensorPosLosVector> vsplos(m, "SensorPosLosVector");
  workspace_group_interface(vsplos);
  common_ndarray(vsplos);
  vsplos
      .def(
          "__array__",
          [](SensorPosLosVector& x, py::object dtype, py::object copy) {
            std::array<size_t, 2> shape = {static_cast<size_t>(x.size()), 5};
            auto np                     = py::module_::import_("numpy");
            auto w =
                py::ndarray<py::numpy, Numeric, py::shape<-1, 5>, py::c_contig>(
                    x.data_handle(), 2, shape.data(), py::cast(x));
            return np.attr("asarray")(w, "dtype"_a = dtype, "copy"_a = copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none(),
          "Returns a :class:`~numpy.ndarray` of the object.")
      .def_prop_rw(
          "value",
          [](py::object& x) { return x.attr("__array__")("copy"_a = false); },
          [](SensorPosLosVector& x, Matrix& y) {
            if (y.ncols() != 5) {
              throw std::runtime_error("Bad shape");
            }
            x.resize(y.nrows());
            std::transform(y.begin(), y.end(), x.begin(), [](const auto& row) {
              return SensorPosLos{.pos = {row[0], row[1], row[2]},
                                  .los = {row[3], row[4]}};
            });
          },
          "A :class:`~numpy.ndarray` of the object.")
      .def("__getstate__",
           [](const py::object& self) {
             return py::make_tuple(self.attr("value"));
           })
      .def("__setstate__", [](py::object& self, const py::tuple& state) {
        self.attr("value") = state[0];
      });

  py::class_<SensorObsel> so(m, "SensorObsel");
  workspace_group_interface(so);
  so.def(py::init<const AscendingGrid&,
                  const SensorPosLosVector&,
                  StokvecMatrix>())
      .def_prop_ro("f_grid", &SensorObsel::f_grid, "Frequency grid")
      .def_prop_ro(
          "weight_matrix", &SensorObsel::weight_matrix, "Weights matrix")
      .def_prop_ro("poslos",
                   &SensorObsel::poslos_grid,
                   "Position and line of sight grid");

  auto a0 =
      py::bind_vector<Array<SensorPosLosVector>, py::rv_policy::reference_internal>(
          m, "ArrayOfSensorPosLosVector");

  auto a1 =
      py::bind_vector<ArrayOfSensorObsel, py::rv_policy::reference_internal>(
          m, "ArrayOfSensorObsel");
  workspace_group_interface(a1);
  vector_interface(a1);

  a1.def(
      "unique_frequency_grids",
      [](const ArrayOfSensorObsel& x) {
        const SensorSimulations simuls = collect_simulations(x);
        std::vector<AscendingGrid> out;
        out.reserve(simuls.size());
        for (auto& k : simuls | std::views::keys) out.push_back(*k);
        return out;
      },
      "List of the unique frequency grids");

  a1.def(
      "unique_poslos_grids",
      [](const ArrayOfSensorObsel& x) {
        const SensorSimulations simuls = collect_simulations(x);
        std::unordered_set<std::shared_ptr<const SensorPosLosVector>> vecs;
        for (auto& k : simuls | std::views::values) {
          vecs.insert(k.begin(), k.end());
        }

        std::vector<SensorPosLosVector> out;
        out.reserve(vecs.size());
        for (auto& v : vecs) out.push_back(*v);
        return out;
      },
      "List of the unique poslos grids");

  a1.def(
      "collect_frequency_grids",
      [](ArrayOfSensorObsel& x) {
        for (Size i = 0; i < x.size(); i++) {
          auto& f  = x[i].f_grid();
          auto& fp = x[i].f_grid_ptr();

          for (Size j = i + 1; j < x.size(); j++) {
            if (x[j].f_grid_ptr() == fp) continue;
            if (std::ranges::equal(x[j].f_grid(), f)) x[j].set_f_grid_ptr(fp);
          }
        }
      },
      R"(Attempts to collect all frequency grids that overlap 1-to-1 between elements.

Will leave non-unique grids alone.
)");

  a1.def(
      "collect_poslos_grids",
      [](ArrayOfSensorObsel& x) {
        for (Size i = 0; i < x.size(); i++) {
          auto& p  = x[i].poslos_grid();
          auto& pp = x[i].poslos_grid_ptr();

          for (Size j = i + 1; j < x.size(); j++) {
            if (x[j].poslos_grid_ptr() == pp) continue;
            if (std::ranges::equal(x[j].poslos_grid(), p)) x[j].set_poslos_grid_ptr(pp);
          }
        }
      },
      R"(Attempts to collect all poslos grids that overlap 1-to-1 between elements.

Will leave non-unique grids alone.
)");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize sensors\n{}", e.what()));
}
}  // namespace Python
