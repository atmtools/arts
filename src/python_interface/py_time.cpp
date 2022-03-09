#include <py_auto_interface.h>
#include <pybind11/chrono.h>

#include <vector>

#include "py_macros.h"

namespace Python {
void py_time(py::module_& m) {
  py::class_<Timer>(m, "Timer")
      .def(py::init<>())
      .PythonInterfaceWorkspaceVariableConversion(Timer)
      .PythonInterfaceFileIO(Timer)
      .def("__repr__", [](Timer&) { return "Timer"; })
      .def_readwrite("running", &Timer::running)
      .def_readwrite("finished", &Timer::finished)
#ifdef TIME_SUPPORT
      .def_readwrite("tms cputime_start", &Timer::cputime_start)
      .def_readwrite("realtime_start", &Timer::realtime_start)
      .def_readwrite("tms cputime_end", &Timer::cputime_end)
      .def_readwrite("realtime_end", &Timer::realtime_end)
#endif
      ;

  py::class_<Time>(m, "Time")
      .def(py::init<>())
      .def(py::init([](std::chrono::system_clock::time_point nt) {
        Time t;
        t.time = nt;
        return t;
      }))
      .def(py::init([](Numeric x) {
        Time t;
        t.Seconds(x);
        return t;
      }))
      .PythonInterfaceWorkspaceVariableConversion(Time)
      .def(py::init<const String&>())
      .PythonInterfaceFileIO(Time)
      .PythonInterfaceBasicRepresentation(Time)
      .def_readwrite("time", &Time::time)
      .def_property(
          "sec",
          [](const Time& t) { return t.Seconds(); },
          [](Time& t, Numeric n) { return t.Seconds(n); })
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
          py::is_operator());
  py::implicitly_convertible<std::chrono::system_clock::time_point, Time>();
  py::implicitly_convertible<py::str, Time>();
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
      py::doc("Converts to a list of datetimes (e.g. for plotting)"));

  PythonInterfaceWorkspaceArray(ArrayOfTime)
      .def(py::init([](const std::vector<std::vector<Time>>& x) {
        ArrayOfArrayOfTime y(x.size());
        std::copy(x.begin(), x.end(), y.begin());
        return y;
      }));
  py::implicitly_convertible<std::vector<std::vector<Time>>,
                             ArrayOfArrayOfTime>();
}
}  // namespace Python
