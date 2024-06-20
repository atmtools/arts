#include <python_interface.h>

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include "artstime.h"
#include "nanobind/stl/bind_vector.h"

namespace Python {
void py_time(py::module_& m) try {
  py::class_<Time>(m, "Time")
      .def("__init__", [](Time* t,std::chrono::system_clock::time_point nt) {
        new (t) Time{};
        t -> time = nt;
      }, "From :class:`datetime.datetime`")
      .def("__init__", [](Time* t,Numeric x) {
       new (t) Time{};
        t -> Seconds(x);
      }, "From :class:`float` seconds from Unix time start")
      .def(py::init<std::string>(), "From :class:`str` of form \"YYYY-MM-DD hh:mm:ss\"")
      .def_rw("time", &Time::time, ":class:`datetime.datetime` The time")
      .def_prop_rw(
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
      .def("__getstate__", [](const Time&t ){
        return std::tuple<std::string>{var_string(t)};
      })
      .def("__setstate__", [](Time*t,const std::tuple<std::string>& state){
        new (t) Time{std::get<0>(state)};
      });
  py::implicitly_convertible<std::chrono::system_clock::time_point, Time>();
  py::implicitly_convertible<std::string, Time>();
  py::implicitly_convertible<Numeric, Time>();

  py::bind_vector<ArrayOfTime>(m, "ArrayOfTime")
      .def_prop_ro(
          "as_datetime",
          [](const ArrayOfTime& in)
              -> std::vector<std::chrono::system_clock::time_point> {
            const Index n = in.size();
            std::vector<std::chrono::system_clock::time_point> out(n);
            for (Index i = 0; i < n; i++) out[i] = in[i].time;
            return out;
          },"A :class:`list` of :class:`datetime.datetime`");
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize time\n", e.what()));
}
}  // namespace Python
