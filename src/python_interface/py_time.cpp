#include <py_auto_interface.h>

#include "py_macros.h"

namespace Python {
void py_time(py::module_& m) {
  py::class_<Timer>(m, "Timer")
      .def(py::init<>())
      .PythonInterfaceWorkspaceVariableConversion(Timer)
      .PythonInterfaceFileIO(Timer)
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
      .PythonInterfaceWorkspaceVariableConversion(Time)
      .def(py::init<const String&>())
      .PythonInterfaceFileIO(Time)
      .PythonInterfaceBasicRepresentation(Time)
      .def("sec", [](const Time& t) { return t.Seconds(); });

  PythonInterfaceWorkspaceArray(Time);
  PythonInterfaceWorkspaceArray(ArrayOfTime);
}
}  // namespace Python