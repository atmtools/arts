#include <python_interface.h>

#include "py_macros.h"

namespace Python {
void py_time(py::module_& m) try {
  artsclass<Time>(m, "Time")
      .def(py::init([]() { return std::make_shared<Time>(); }), "Current time")
      .def(py::init([](std::chrono::system_clock::time_point nt) {
        Time t;
        t.time = nt;
        return t;
      }), "From :class:`datetime.datetime`")
      .def(py::init([](Numeric x) {
        Time t;
        t.Seconds(x);
        return t;
      }), "From :class:`float` seconds from Unix time start")
      .PythonInterfaceCopyValue(Time)
      .PythonInterfaceWorkspaceVariableConversion(Time)
      .def(py::init([](const std::string& s) { return std::make_shared<Time>(s); }), "From :class:`str` of form \"YYYY-MM-DD hh:mm:ss\"")
      .PythonInterfaceFileIO(Time)
      .PythonInterfaceBasicRepresentation(Time)
      .def_readwrite("time", &Time::time, ":class:`datetime.datetime` The time")
      .def_property(
          "sec",
          [](const Time& t) { return t.Seconds(); },
          [](Time& t, Numeric n) { return t.Seconds(n); }, ":class:`float` Time from Unix start")
      .def(
          "__add__",
          [](Time& t, Numeric n) {
            Time t2;
            t2.Seconds(t.Seconds() + n);
            return t2;
          },
          py::is_operator())
      .def(
          "__radd__",
          [](Time& t, Numeric n) {
            Time t2;
            t2.Seconds(t.Seconds() + n);
            return t2;
          },
          py::is_operator())
      .def(
          "__sub__",
          [](Time& t, Numeric n) {
            Time t2;
            t2.Seconds(t.Seconds() - n);
            return t2;
          },
          py::is_operator())
      .def(
          "__rsub__",
          [](Time& t, Numeric n) {
            Time t2;
            t2.Seconds(n - t.Seconds());
            return t2;
          },
          py::is_operator())
      .def(py::pickle(
          [](const py::object& self) {
            return py::make_tuple(self.attr("time"));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<Time>()(t[0]).cast<Time>();
          }))
      .PythonInterfaceWorkspaceDocumentation(Time);
  py::implicitly_convertible<std::chrono::system_clock::time_point, Time>();
  py::implicitly_convertible<std::string, Time>();
  py::implicitly_convertible<Numeric, Time>();

  PythonInterfaceWorkspaceArray(Time).def_property_readonly(
      "as_datetime",
      [](const ArrayOfTime& in)
          -> std::vector<std::chrono::system_clock::time_point> {
        const Index n = in.nelem();
        std::vector<std::chrono::system_clock::time_point> out(n);
        for (Index i = 0; i < n; i++) out[i] = in[i].time;
        return out;
      },
      py::doc("A :class:`list` of :class:`datetime.datetime`"));

  PythonInterfaceWorkspaceArray(ArrayOfTime)
      .def(py::init([](const std::vector<std::vector<Time>>& x) {
        ArrayOfArrayOfTime y(x.size());
        std::copy(x.begin(), x.end(), y.begin());
        return y;
      }), "From nested lists");
  py::implicitly_convertible<std::vector<std::vector<Time>>,
                             ArrayOfArrayOfTime>();
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize time\n", e.what()));
}
}  // namespace Python
