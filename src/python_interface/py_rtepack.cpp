#include <pybind11/pybind11.h>
#include <python_interface.h>
#include <rtepack.h>

#include "enums.h"
#include "py_macros.h"

namespace Python {
void py_rtepack(py::module_ &m) try {
  py_staticStokvec(m)
      .def(py::init<Numeric>())
      .def(py::init(
          [](const PolarizationChoice p) { return rtepack::to_stokvec(p); }))
      .def(py::init([](const String &p) {
        return rtepack::to_stokvec(to<PolarizationChoice>(p));
      }))
      .def_static(
          "linpol",
          [](const Numeric angle) {
            return Stokvec{1.0,
                           Conversion::cosd(2.0 * angle),
                           Conversion::sind(2.0 * angle),
                           0.0};
          },
          "Returns [1.0, cos(2*angle), sin(2*angle), 0.0], the linear polarization vector for a given angle",
          py::arg("angle"))
      .def_static(
          "cirpol",
          [](const Numeric angle) {
            return Stokvec{1.0, 0.0, 0.0, Conversion::sind(angle)};
          },
          "Returns [1.0, 0.0, 0.0, sin(angle)], the circular polarization vector for a given phase delay angle",
          py::arg("angle"))
      .def(py::init([](std::array<Numeric, 4> a) {
        return Stokvec{a[0], a[1], a[2], a[3]};
      }))
      .def_buffer([](Stokvec &x) -> py::buffer_info {
        return py::buffer_info(x.data.data(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               1,
                               {4},
                               {sizeof(Numeric)});
      })
      .def_property("value",
                    py::cpp_function(
                        [](Stokvec &x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Stokvec &x, Stokvec &y) { x = y; })
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties;
  py::implicitly_convertible<Numeric, Stokvec>();
  py::implicitly_convertible<std::array<Numeric, 4>, Stokvec>();
  py::implicitly_convertible<PolarizationChoice, Stokvec>();
  py::implicitly_convertible<String, Stokvec>();

  py_staticStokvecVector(m)
      .def(py::init<std::vector<Numeric>>())
      .def(py::init<std::vector<Stokvec>>())
      .def_buffer([](StokvecVector &x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            2,
            {static_cast<ssize_t>(x.nelem()), static_cast<ssize_t>(4)},
            {static_cast<ssize_t>(4 * sizeof(Numeric)),
             static_cast<ssize_t>(sizeof(Numeric))});
      })
      .def_property("value",
                    py::cpp_function(
                        [](StokvecVector &x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](StokvecVector &x, StokvecVector &y) { x = y; })
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties;
  py::implicitly_convertible<std::vector<Numeric>, StokvecVector>();
  py::implicitly_convertible<std::vector<Stokvec>, StokvecVector>();

  py_staticStokvecMatrix(m)
      .def_buffer([](StokvecMatrix &x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            3,
            {static_cast<ssize_t>(x.nrows()),
             static_cast<ssize_t>(x.ncols()),
             static_cast<ssize_t>(4)},
            {static_cast<ssize_t>(x.ncols() * 4 * sizeof(Numeric)),
             static_cast<ssize_t>(4 * sizeof(Numeric)),
             static_cast<ssize_t>(sizeof(Numeric))});
      })
      .def_property("value",
                    py::cpp_function(
                        [](StokvecMatrix &x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](StokvecMatrix &x, StokvecMatrix &y) { x = y; })
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties;

  py_staticStokvecTensor3(m)
      .def_buffer([](StokvecTensor3 &x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            4,
            {static_cast<ssize_t>(x.npages()),
             static_cast<ssize_t>(x.nrows()),
             static_cast<ssize_t>(x.ncols()),
             static_cast<ssize_t>(4)},
            {static_cast<ssize_t>(x.nrows() * x.ncols() * 4 * sizeof(Numeric)),
             static_cast<ssize_t>(x.ncols() * 4 * sizeof(Numeric)),
             static_cast<ssize_t>(4 * sizeof(Numeric)),
             static_cast<ssize_t>(sizeof(Numeric))});
      })
      .def_property("value",
                    py::cpp_function(
                        [](StokvecTensor3 &x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](StokvecTensor3 &x, StokvecTensor3 &y) { x = y; })
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties;

  py_staticStokvecTensor4(m)
      .def_buffer([](StokvecTensor4 &x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            5,
            {static_cast<ssize_t>(x.nbooks()),
             static_cast<ssize_t>(x.npages()),
             static_cast<ssize_t>(x.nrows()),
             static_cast<ssize_t>(x.ncols()),
             static_cast<ssize_t>(4)},
            {static_cast<ssize_t>(x.npages() * x.nrows() * x.ncols() * 4 *
                                  sizeof(Numeric)),
             static_cast<ssize_t>(x.nrows() * x.ncols() * 4 * sizeof(Numeric)),
             static_cast<ssize_t>(x.ncols() * 4 * sizeof(Numeric)),
             static_cast<ssize_t>(4 * sizeof(Numeric)),
             static_cast<ssize_t>(sizeof(Numeric))});
      })
      .def_property("value",
                    py::cpp_function(
                        [](StokvecTensor4 &x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](StokvecTensor4 &x, StokvecTensor4 &y) { x = y; })
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties;

  py_staticStokvecTensor5(m)
      .def_buffer([](StokvecTensor5 &x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            6,
            {static_cast<ssize_t>(x.nshelves()),
             static_cast<ssize_t>(x.nbooks()),
             static_cast<ssize_t>(x.npages()),
             static_cast<ssize_t>(x.nrows()),
             static_cast<ssize_t>(x.ncols()),
             static_cast<ssize_t>(4)},
            {static_cast<ssize_t>(x.nbooks() * x.npages() * x.nrows() *
                                  x.ncols() * 4 * sizeof(Numeric)),
             static_cast<ssize_t>(x.npages() * x.nrows() * x.ncols() * 4 *
                                  sizeof(Numeric)),
             static_cast<ssize_t>(x.nrows() * x.ncols() * 4 * sizeof(Numeric)),
             static_cast<ssize_t>(x.ncols() * 4 * sizeof(Numeric)),
             static_cast<ssize_t>(4 * sizeof(Numeric)),
             static_cast<ssize_t>(sizeof(Numeric))});
      })
      .def_property("value",
                    py::cpp_function(
                        [](StokvecTensor5 &x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](StokvecTensor5 &x, StokvecTensor5 &y) { x = y; })
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties;

  py_staticStokvecTensor6(m)
      .def_buffer([](StokvecTensor6 &x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            7,
            {static_cast<ssize_t>(x.nvitrines()),
             static_cast<ssize_t>(x.nshelves()),
             static_cast<ssize_t>(x.nbooks()),
             static_cast<ssize_t>(x.npages()),
             static_cast<ssize_t>(x.nrows()),
             static_cast<ssize_t>(x.ncols()),
             static_cast<ssize_t>(4)},
            {static_cast<ssize_t>(x.nshelves() * x.nbooks() * x.npages() *
                                  x.nrows() * x.ncols() * 4 * sizeof(Numeric)),
             static_cast<ssize_t>(x.nbooks() * x.npages() * x.nrows() *
                                  x.ncols() * 4 * sizeof(Numeric)),
             static_cast<ssize_t>(x.npages() * x.nrows() * x.ncols() * 4 *
                                  sizeof(Numeric)),
             static_cast<ssize_t>(x.nrows() * x.ncols() * 4 * sizeof(Numeric)),
             static_cast<ssize_t>(x.ncols() * 4 * sizeof(Numeric)),
             static_cast<ssize_t>(4 * sizeof(Numeric)),
             static_cast<ssize_t>(sizeof(Numeric))});
      })
      .def_property("value",
                    py::cpp_function(
                        [](StokvecTensor6 &x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](StokvecTensor6 &x, StokvecTensor6 &y) { x = y; })
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties;

  py_staticArrayOfStokvecVector(m);
  py_staticArrayOfStokvecMatrix(m);
  py_staticArrayOfArrayOfStokvecVector(m);
  py_staticArrayOfArrayOfStokvecMatrix(m);

  py_staticPropmat(m)
      .def(py::init<Numeric>())
      .def(py::init<Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric>())
      .def_buffer([](Propmat &x) -> py::buffer_info {
        return py::buffer_info(x.data.data(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               1,
                               {7},
                               {sizeof(Numeric)});
      })
      .def_property("value",
                    py::cpp_function(
                        [](Propmat &x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Propmat &x, Propmat &y) { x = y; })
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties;
  py::implicitly_convertible<Numeric, Propmat>();

  py_staticPropmatVector(m)
      .def(py::init<std::vector<Numeric>>())
      .def(py::init<std::vector<Propmat>>())
      .def_buffer([](PropmatVector &x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            2,
            {static_cast<ssize_t>(x.nelem()), static_cast<ssize_t>(7)},
            {static_cast<ssize_t>(7 * sizeof(Numeric)),
             static_cast<ssize_t>(sizeof(Numeric))});
      })
      .def_property("value",
                    py::cpp_function(
                        [](PropmatVector &x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](PropmatVector &x, PropmatVector &y) { x = y; })
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties;
  py::implicitly_convertible<std::vector<Numeric>, PropmatVector>();
  py::implicitly_convertible<std::vector<Propmat>, PropmatVector>();

  py_staticPropmatMatrix(m)
      .def_buffer([](PropmatMatrix &x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            3,
            {static_cast<ssize_t>(x.nrows()),
             static_cast<ssize_t>(x.ncols()),
             static_cast<ssize_t>(7)},
            {static_cast<ssize_t>(7 * x.ncols() * sizeof(Numeric)),
             static_cast<ssize_t>(7 * sizeof(Numeric)),
             static_cast<ssize_t>(sizeof(Numeric))});
      })
      .def_property("value",
                    py::cpp_function(
                        [](PropmatMatrix &x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](PropmatMatrix &x, PropmatMatrix &y) { x = y; })
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties;

  py_staticArrayOfPropmatVector(m);
  py_staticArrayOfPropmatMatrix(m);
  py_staticArrayOfArrayOfPropmatVector(m);
  py_staticArrayOfArrayOfPropmatMatrix(m);

  py_staticMuelmat(m)
      .def(py::init<Numeric>())
      .def(py::init<Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric>())
      .def_buffer([](Muelmat &x) -> py::buffer_info {
        return py::buffer_info(x.data.data(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               2,
                               {4, 4},
                               {4 * sizeof(Numeric), sizeof(Numeric)});
      })
      .def_property("value",
                    py::cpp_function(
                        [](Muelmat &x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Muelmat &x, Muelmat &y) { x = y; })
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties;
  py::implicitly_convertible<Numeric, Muelmat>();

  py_staticMuelmatVector(m)
      .def(py::init<std::vector<Numeric>>())
      .def(py::init<std::vector<Muelmat>>())
      .def_buffer([](MuelmatVector &x) -> py::buffer_info {
        return py::buffer_info(x.data_handle(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               3,
                               {static_cast<ssize_t>(x.nelem()),
                                static_cast<ssize_t>(4),
                                static_cast<ssize_t>(4)},
                               {static_cast<ssize_t>(16 * sizeof(Numeric)),
                                static_cast<ssize_t>(4 * sizeof(Numeric)),
                                static_cast<ssize_t>(sizeof(Numeric))});
      })
      .def_property("value",
                    py::cpp_function(
                        [](MuelmatVector &x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](MuelmatVector &x, MuelmatVector &y) { x = y; })
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties;
  py::implicitly_convertible<std::vector<Numeric>, MuelmatVector>();
  py::implicitly_convertible<std::vector<Muelmat>, MuelmatVector>();

  py_staticMuelmatMatrix(m)
      .def_buffer([](MuelmatMatrix &x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            4,
            {static_cast<ssize_t>(x.nrows()),
             static_cast<ssize_t>(x.ncols()),
             static_cast<ssize_t>(4),
             static_cast<ssize_t>(4)},
            {static_cast<ssize_t>(x.ncols() * 16 * sizeof(Numeric)),
             static_cast<ssize_t>(16 * sizeof(Numeric)),
             static_cast<ssize_t>(4 * sizeof(Numeric)),
             static_cast<ssize_t>(sizeof(Numeric))});
      })
      .def_property("value",
                    py::cpp_function(
                        [](MuelmatMatrix &x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](MuelmatMatrix &x, MuelmatMatrix &y) { x = y; })
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties;

  py_staticArrayOfMuelmatVector(m);
  py_staticArrayOfMuelmatMatrix(m);
  py_staticArrayOfArrayOfMuelmatVector(m);
  py_staticArrayOfArrayOfMuelmatMatrix(m);
} catch (std::exception &e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize rtepack\n", e.what()));
}
}  // namespace Python
