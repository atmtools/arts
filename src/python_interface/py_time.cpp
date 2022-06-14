#include <py_auto_interface.h>
#include <pybind11/chrono.h>
#include <pybind11/pybind11.h>

#include "py_macros.h"

namespace Python {
void py_time(py::module_& m) {
  py::class_<std::clock_t>(m, "clock_t")
      .def(py::init([]() { return std::clock(); }))
      .def(py::pickle([](const clock_t& self) { return py::make_tuple(self); },
                      [](const py::tuple& t) {
                        ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
                        return new clock_t{t[0].cast<clock_t>()};
                      }));

  py::class_<Timer>(m, "Timer")
      .def(py::init([]() { return new Timer{}; }))
      .PythonInterfaceCopyValue(Timer)
      .PythonInterfaceWorkspaceVariableConversion(Timer)
      .PythonInterfaceFileIO(Timer)
      .def("__repr__", [](Timer&) { return "Timer"; })
      .def_readwrite("running", &Timer::running)
      .def_readwrite("finished", &Timer::finished)
      .def_readwrite("cputime_start", &Timer::cputime_start)
      .def_readwrite("realtime_start", &Timer::realtime_start)
      .def_readwrite("cputime_end", &Timer::cputime_end)
      .def_readwrite("realtime_end", &Timer::realtime_end)
      .def(py::pickle(
          [](const Timer& self) {
            return py::make_tuple(self.running,
                                  self.finished,
                                  self.cputime_start,
                                  self.realtime_start,
                                  self.cputime_end,
                                  self.realtime_end
            );
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 6, "Invalid state!")

            auto* out = new Timer{};
            out->running = t[0].cast<bool>();
            out->finished = t[1].cast<bool>();
            out->cputime_start = t[2].cast<std::clock_t>();
            out->realtime_start = t[3].cast<std::chrono::time_point<std::chrono::high_resolution_clock>>();
            out->cputime_end = t[4].cast<std::clock_t>();
            out->realtime_end = t[5].cast<std::chrono::time_point<std::chrono::high_resolution_clock>>();

            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(Timer);

  py::class_<Time>(m, "Time")
      .def(py::init([]() { return new Time{}; }))
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
      .PythonInterfaceCopyValue(Time)
      .PythonInterfaceWorkspaceVariableConversion(Time)
      .def(py::init([](const std::string& s) { return new Time{s}; }))
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
