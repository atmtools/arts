#include <py_auto_interface.h>

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
      .def(py::init([](Numeric x){Time t; t.Seconds(x); return t;}))
      .PythonInterfaceWorkspaceVariableConversion(Time)
      .def(py::init<const String&>())
      .PythonInterfaceFileIO(Time)
      .PythonInterfaceBasicRepresentation(Time)
      .def_property("sec", [](const Time& t) { return t.Seconds(); }, [](Time& t, Numeric n) { return t.Seconds(n); })
      .def("__add__", [](Time& t, Numeric n){Time t2; t2.Seconds(t.Seconds() + n); return t2;}, py::is_operator())
      .def("__radd__", [](Time& t, Numeric n){Time t2; t2.Seconds(t.Seconds() + n); return t2;}, py::is_operator())
      .def("__sub__", [](Time& t, Numeric n){Time t2; t2.Seconds(t.Seconds() - n); return t2;}, py::is_operator())
      .def("__rsub__", [](Time& t, Numeric n){Time t2; t2.Seconds(n - t.Seconds()); return t2;}, py::is_operator());
  py::implicitly_convertible<py::str, Time>();
  py::implicitly_convertible<Numeric, Time>();

  PythonInterfaceWorkspaceArray(Time);
  PythonInterfaceWorkspaceArray(ArrayOfTime);
}
}  // namespace Python
