#include <python_interface.h>

#include <type_traits>
#include <variant>

#include "py_macros.h"

#include <rtepack.h>

namespace Python {
void py_rtepack(py::module_ &m) try {
  artsclass<Stokvec>(m, "Stokvec", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Numeric>())
      .def(py::init<std::array<Numeric, 4>>())
      .def_buffer([](Stokvec &x) -> py::buffer_info {
        return py::buffer_info(x.data.data(), sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(), 1, {4},
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
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .PythonInterfaceWorkspaceVariableConversion(Stokvec)
      .PythonInterfaceCopyValue(Stokvec)
      .PythonInterfaceBasicRepresentation(Stokvec)
      .PythonInterfaceFileIO(Stokvec)
      .def(py::pickle(
          [](const py::object &self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<Stokvec>()(t[0]).cast<Stokvec>();
          }))
      .PythonInterfaceWorkspaceDocumentation(Stokvec);
  py::implicitly_convertible<Numeric, Stokvec>();
  py::implicitly_convertible<std::array<Numeric, 4>, Stokvec>();

  artsclass<StokvecVector>(m, "StokvecVector", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<std::vector<Numeric>>())
      .def(py::init<std::vector<Stokvec>>())
      .def_buffer([](StokvecVector &x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(), sizeof(Numeric),
            py::format_descriptor<Numeric>::format(), 2,
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
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .PythonInterfaceWorkspaceVariableConversion(StokvecVector)
      .PythonInterfaceCopyValue(StokvecVector)
      .PythonInterfaceBasicRepresentation(StokvecVector)
      .PythonInterfaceFileIO(StokvecVector)
      .def(py::pickle(
          [](const py::object &self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<StokvecVector>()(t[0]).cast<StokvecVector>();
          }))
      .PythonInterfaceWorkspaceDocumentation(StokvecVector);
  py::implicitly_convertible<std::vector<Numeric>, StokvecVector>();
  py::implicitly_convertible<std::vector<Stokvec>, StokvecVector>();

  artsclass<StokvecMatrix>(m, "StokvecMatrix", py::buffer_protocol())
      .def(py::init<>())
      .def_buffer([](StokvecMatrix &x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(), sizeof(Numeric),
            py::format_descriptor<Numeric>::format(), 3,
            {static_cast<ssize_t>(x.nrows()), static_cast<ssize_t>(x.ncols()),
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
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .PythonInterfaceWorkspaceVariableConversion(StokvecMatrix)
      .PythonInterfaceCopyValue(StokvecMatrix)
      .PythonInterfaceBasicRepresentation(StokvecMatrix)
      .PythonInterfaceFileIO(StokvecMatrix)
      .def(py::pickle(
          [](const py::object &self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<StokvecMatrix>()(t[0]).cast<StokvecMatrix>();
          }))
      .PythonInterfaceWorkspaceDocumentation(StokvecMatrix);

  PythonInterfaceWorkspaceArray(StokvecVector);
  PythonInterfaceWorkspaceArray(StokvecMatrix);
  PythonInterfaceWorkspaceArray(ArrayOfStokvecVector);
  PythonInterfaceWorkspaceArray(ArrayOfStokvecMatrix);

  artsclass<Propmat>(m, "Propmat", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Numeric>())
      .def(py::init<Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric>())
      .def_buffer([](Propmat &x) -> py::buffer_info {
        return py::buffer_info(x.data.data(), sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(), 1, {7},
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
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .PythonInterfaceWorkspaceVariableConversion(Propmat)
      .PythonInterfaceCopyValue(Propmat)
      .PythonInterfaceBasicRepresentation(Propmat)
      .PythonInterfaceFileIO(Propmat)
      .def(py::pickle(
          [](const py::object &self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<Propmat>()(t[0]).cast<Propmat>();
          }))
      .PythonInterfaceWorkspaceDocumentation(Propmat);
  py::implicitly_convertible<Numeric, Propmat>();

  artsclass<PropmatVector>(m, "PropmatVector", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<std::vector<Numeric>>())
      .def(py::init<std::vector<Propmat>>())
      .def_buffer([](PropmatVector &x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(), sizeof(Numeric),
            py::format_descriptor<Numeric>::format(), 2,
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
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .PythonInterfaceWorkspaceVariableConversion(PropmatVector)
      .PythonInterfaceCopyValue(PropmatVector)
      .PythonInterfaceBasicRepresentation(PropmatVector)
      .PythonInterfaceFileIO(PropmatVector)
      .def(py::pickle(
          [](const py::object &self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<PropmatVector>()(t[0]).cast<PropmatVector>();
          }))
      .PythonInterfaceWorkspaceDocumentation(PropmatVector);
  py::implicitly_convertible<std::vector<Numeric>, PropmatVector>();
  py::implicitly_convertible<std::vector<Propmat>, PropmatVector>();

  artsclass<PropmatMatrix>(m, "PropmatMatrix", py::buffer_protocol())
      .def(py::init<>())
      .def_buffer([](PropmatMatrix &x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(), sizeof(Numeric),
            py::format_descriptor<Numeric>::format(), 3,
            {static_cast<ssize_t>(x.nrows()), static_cast<ssize_t>(x.ncols()),
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
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .PythonInterfaceWorkspaceVariableConversion(PropmatMatrix)
      .PythonInterfaceCopyValue(PropmatMatrix)
      .PythonInterfaceBasicRepresentation(PropmatMatrix)
      .PythonInterfaceFileIO(PropmatMatrix)
      .def(py::pickle(
          [](const py::object &self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<PropmatMatrix>()(t[0]).cast<PropmatMatrix>();
          }))
      .PythonInterfaceWorkspaceDocumentation(PropmatMatrix);

  PythonInterfaceWorkspaceArray(PropmatVector);
  PythonInterfaceWorkspaceArray(PropmatMatrix);
  PythonInterfaceWorkspaceArray(ArrayOfPropmatVector);
  PythonInterfaceWorkspaceArray(ArrayOfPropmatMatrix);

  artsclass<Muelmat>(m, "Muelmat", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Numeric>())
      .def(py::init<Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric>())
      .def_buffer([](Muelmat &x) -> py::buffer_info {
        return py::buffer_info(x.data.data(), sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(), 2,
                               {4, 4}, {4 * sizeof(Numeric), sizeof(Numeric)});
      })
      .def_property("value",
                    py::cpp_function(
                        [](Muelmat &x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Muelmat &x, Muelmat &y) { x = y; })
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .PythonInterfaceWorkspaceVariableConversion(Muelmat)
      .PythonInterfaceCopyValue(Muelmat)
      .PythonInterfaceBasicRepresentation(Muelmat)
      .PythonInterfaceFileIO(Muelmat)
      .def(py::pickle(
          [](const py::object &self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<Muelmat>()(t[0]).cast<Muelmat>();
          }))
      .PythonInterfaceWorkspaceDocumentation(Muelmat);
  py::implicitly_convertible<Numeric, Muelmat>();

  artsclass<MuelmatVector>(m, "MuelmatVector", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<std::vector<Numeric>>())
      .def(py::init<std::vector<Muelmat>>())
      .def_buffer([](MuelmatVector &x) -> py::buffer_info {
        return py::buffer_info(x.data_handle(), sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(), 3,
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
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .PythonInterfaceWorkspaceVariableConversion(MuelmatVector)
      .PythonInterfaceCopyValue(MuelmatVector)
      .PythonInterfaceBasicRepresentation(MuelmatVector)
      .PythonInterfaceFileIO(MuelmatVector)
      .def(py::pickle(
          [](const py::object &self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<MuelmatVector>()(t[0]).cast<MuelmatVector>();
          }))
      .PythonInterfaceWorkspaceDocumentation(MuelmatVector);
  py::implicitly_convertible<std::vector<Numeric>, MuelmatVector>();
  py::implicitly_convertible<std::vector<Muelmat>, MuelmatVector>();

  artsclass<MuelmatMatrix>(m, "MuelmatMatrix", py::buffer_protocol())
      .def(py::init<>())
      .def_buffer([](MuelmatMatrix &x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(), sizeof(Numeric),
            py::format_descriptor<Numeric>::format(), 4,
            {static_cast<ssize_t>(x.nrows()), static_cast<ssize_t>(x.ncols()),
             static_cast<ssize_t>(4), static_cast<ssize_t>(4)},
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
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .PythonInterfaceWorkspaceVariableConversion(MuelmatMatrix)
      .PythonInterfaceCopyValue(MuelmatMatrix)
      .PythonInterfaceBasicRepresentation(MuelmatMatrix)
      .PythonInterfaceFileIO(MuelmatMatrix)
      .def(py::pickle(
          [](const py::object &self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<MuelmatMatrix>()(t[0]).cast<MuelmatMatrix>();
          }))
      .PythonInterfaceWorkspaceDocumentation(MuelmatMatrix);

  PythonInterfaceWorkspaceArray(MuelmatVector);
  PythonInterfaceWorkspaceArray(MuelmatMatrix);
  PythonInterfaceWorkspaceArray(ArrayOfMuelmatVector);
  PythonInterfaceWorkspaceArray(ArrayOfMuelmatMatrix);
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize rtepack\n", e.what()));
}
} // namespace Python
