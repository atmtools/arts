#include <python_interface.h>

#include "py_macros.h"

namespace Python {
void py_time(py::module_& m) try {
  py::class_<std::clock_t>(m, "clock_t")
      .def(py::init([]() { return std::clock(); }), "Default clock")
      .def(py::pickle([](const clock_t& self) { return py::make_tuple(self); },
                      [](const py::tuple& t) {
                        ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
                        return clock_t{t[0].cast<clock_t>()};
                      })).doc() = "The C++ standard clock type";

  py::class_<Timer>(m, "Timer")
      .def(py::init([]() { return std::make_unique<Timer>(); }), "Empty timer")
      .PythonInterfaceCopyValue(Timer)
      .PythonInterfaceWorkspaceVariableConversion(Timer)
      .PythonInterfaceFileIO(Timer)
      .def("__repr__", [](Timer&) { return "Timer"; })
      .def_readwrite("running", &Timer::running, ":class:`bool`")
      .def_readwrite("finished", &Timer::finished, ":class:`bool`")
      .def_readwrite("cputime_start", &Timer::cputime_start, ":class:`~pyarts.arts.clock_t`")
      .def_readwrite("realtime_start", &Timer::realtime_start, ":class:`~pyarts.arts.clock_t`")
      .def_readwrite("cputime_end", &Timer::cputime_end, ":class:`~pyarts.arts.clock_t`")
      .def_readwrite("realtime_end", &Timer::realtime_end, ":class:`~pyarts.arts.clock_t`")
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

            auto out = std::make_unique<Timer>();
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
      .def(py::init([]() { return std::make_unique<Time>(); }), "Current time")
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
      .def(py::init([](const std::string& s) { return std::make_unique<Time>(s); }), "From :class:`str` of form \"YYYY-MM-DD hh:mm:ss\"")
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
