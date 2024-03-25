#include <obsel.h>
#include <py_auto_wsg_init.h>
#include <python_interface.h>

#include "py_macros.h"

namespace Python {
void py_sensor(py::module_& m) try {
  py_staticSensorPosLos(m)
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](SensorPosLos& x) -> py::buffer_info {
        return py::buffer_info(&x,
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               1,
                               {static_cast<ssize_t>(5)},
                               {static_cast<ssize_t>(sizeof(Numeric))});
      })
      .def(py::init<Vector3, Vector2>(), "From pos and los")
      .def_readwrite("pos", &SensorPosLos::pos, "Position")
      .def_readwrite("los", &SensorPosLos::los, "Line of sight");

  py_staticSensorPosLosVector(m)
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](SensorPosLosVector& x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            2,
            {static_cast<ssize_t>(x.size()), static_cast<ssize_t>(5)},
            {static_cast<ssize_t>(5 * sizeof(Numeric)),
             static_cast<ssize_t>(sizeof(Numeric))});
      })
      .def_property(
          "value",
          py::cpp_function(
              [](SensorPosLosVector& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
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
          py::doc(":class:`~pyarts.arts.Matrix`"))
      .def(py::pickle(
          [](const py::object& self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<SensorPosLosVector>()(t[0])
                .cast<SensorPosLosVector>();
          }));

  py_staticSensorObsel(m)
      .def(py::init<Vector,
                    AscendingGrid,
                    MuelmatVector,
                    SensorPosLosVector,
                    Stokvec>(),
           "From pos and los")
      .def_readwrite("f_grid_w", &SensorObsel::f_grid_w, "Frequency weights")
      .def_readwrite("f_grid", &SensorObsel::f_grid, "Frequency grid")
      .def_readwrite("poslos_w",
                     &SensorObsel::poslos_grid_w,
                     "Position and line of sight weights")
      .def_readwrite("poslos",
                     &SensorObsel::poslos_grid,
                     "Position and line of sight grid")
      .def_readwrite(
          "polarization", &SensorObsel::polarization, "Polarization sampling")
      .def("ok", &SensorObsel::ok, "Check if the obsel is valid");

  py_staticArrayOfSensorObsel(m);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize rtepack\n", e.what()));
}
}  // namespace Python
