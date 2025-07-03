#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/chrono.h>
#include <nanobind/stl/string.h>
#include <python_interface.h>

#include "artstime.h"
#include "hpy_arts.h"
#include "hpy_vector.h"
#include "nanobind/stl/bind_vector.h"

namespace Python {
void py_time(py::module_& m) try {
  py::class_<Time> time(m, "Time");

  generic_interface(time);

  time.def(
          "__init__",
          [](Time* t, std::chrono::system_clock::time_point nt) {
            new (t) Time{};
            t->time = nt;
          },
          "From :class:`datetime.datetime`")
      .def(
          "__init__",
          [](Time* t, Numeric x) {
            new (t) Time{};
            t->Seconds(x);
          },
          "From :class:`float` seconds from Unix time start")
      .def(py::init<std::string>(),
           "From :class:`str` of form \"YYYY-MM-DD hh:mm:ss\"")
      .def_rw("time", &Time::time, ":class:`datetime.datetime` The time")
      .def_prop_rw(
          "sec",
          [](const Time& t) { return t.Seconds(); },
          [](Time& t, Numeric n) { return t.Seconds(n); },
          ":class:`float` Time from Unix start")
      .def(
          "__add__",
          [](Time& t, Numeric n) {
            Time t2;
            t2.Seconds(t.Seconds() + n);
            return t2;
          },
          py::is_operator(),
          "n"_a,
          "Allows `self + n`")
      .def(
          "__radd__",
          [](Time& t, Numeric n) {
            Time t2;
            t2.Seconds(t.Seconds() + n);
            return t2;
          },
          py::is_operator(),
          "n"_a,
          "Allows `n + self`")
      .def(
          "__sub__",
          [](Time& t, Numeric n) {
            Time t2;
            t2.Seconds(t.Seconds() - n);
            return t2;
          },
          py::is_operator(),
          "n"_a,
          "Allows `self - n`")
      .def(
          "__rsub__",
          [](Time& t, Numeric n) {
            Time t2;
            t2.Seconds(n - t.Seconds());
            return t2;
          },
          py::is_operator(),
          "n"_a,
          "Allows `n - self`")
      .def("__getstate__",
           [](const Time& t) {
             return std::tuple<std::string>{std::format("{}", t)};
           })
      .def("__setstate__", [](Time* t, const std::tuple<std::string>& state) {
        new (t) Time{std::get<0>(state)};
      });

  py::implicitly_convertible<std::chrono::system_clock::time_point, Time>();
  py::implicitly_convertible<std::string, Time>();
  py::implicitly_convertible<Numeric, Time>();

  auto a1 = py::bind_vector<ArrayOfTime, py::rv_policy::reference_internal>(
                m, "ArrayOfTime")
                .def_prop_ro(
                    "as_datetime",
                    [](const ArrayOfTime& in)
                        -> std::vector<std::chrono::system_clock::time_point> {
                      const Index n = in.size();
                      std::vector<std::chrono::system_clock::time_point> out(n);
                      for (Index i = 0; i < n; i++) out[i] = in[i].time;
                      return out;
                    },
                    "A :class:`list` of :class:`datetime.datetime`");

  generic_interface(a1);

  vector_interface(a1);

  auto a2 =
      py::bind_vector<ArrayOfArrayOfTime, py::rv_policy::reference_internal>(
          m, "ArrayOfArrayOfTime");
  generic_interface(a2);
  vector_interface(a2);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize time\n{}", e.what()));
}
}  // namespace Python
