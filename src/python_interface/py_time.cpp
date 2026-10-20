#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/chrono.h>
#include <nanobind/stl/string.h>
#include <python_interface.h>

#include <chrono>
#include <cmath>

#include "artstime.h"
#include "hpy_arts.h"
#include "hpy_vector.h"

namespace Python {
//! A Python datetime.datetime, where naive ones are in UTC and aware ones keep their offset from UTC
struct DateTime {
  Time time;
};
}  // namespace Python

namespace nanobind::detail {
template <> struct type_caster<Python::DateTime> {
  NB_TYPE_CASTER(Python::DateTime, const_name("datetime.datetime"))

  bool from_python(handle src, uint32_t, cleanup_list* cleanup) noexcept try {
    using namespace std::chrono;

    if (not isinstance(src, module_::import_("datetime").attr("datetime"))) return false;

    int yy, mon, dd, hh, mm, ss, us;
    if (unpack_datetime(src.ptr(), &yy, &mon, &dd, &hh, &mm, &ss, &us, cleanup) <= 0) {
      PyErr_Clear();
      return false;
    }

    minutes      offset{0};
    const object utc_offset = src.attr("utcoffset")();
    if (not utc_offset.is_none()) {
      const double sec = cast<double>(utc_offset.attr("total_seconds")());
      if (std::fmod(sec, 60.0) != 0.0) return false;
      offset = minutes{static_cast<int>(sec / 60.0)};
    }

    value = Python::DateTime{
        Time{year{yy} / mon / dd,
             duration_cast<Time::InternalTimeStep>(hours{hh} + minutes{mm} + seconds{ss} + microseconds{us}),
             offset}};
    return true;
  } catch (...) {
    PyErr_Clear();
    return false;
  }

  static handle from_cpp(const Python::DateTime& x, rv_policy, cleanup_list*) noexcept try {
    using namespace std::chrono;

    const year_month_day                   date = x.time.date();
    const hh_mm_ss<Time::InternalTimeStep> tod  = x.time.time_of_day();

    const module_ datetime = module_::import_("datetime");
    const object tz = datetime.attr("timezone")(datetime.attr("timedelta")(arg("minutes") = x.time.utc_offset.count()));

    return datetime
        .attr("datetime")(static_cast<int>(date.year()),
                          static_cast<unsigned>(date.month()),
                          static_cast<unsigned>(date.day()),
                          tod.hours().count(),
                          tod.minutes().count(),
                          tod.seconds().count(),
                          duration_cast<microseconds>(tod.subseconds()).count(),
                          arg("tzinfo") = tz)
        .release();
  } catch (python_error& e) {
    e.restore();
    return handle();
  } catch (const std::exception& e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    return handle();
  }
};
}  // namespace nanobind::detail

namespace Python {
void py_time(py::module_& m) try {
  py::class_<Time> time(m, "Time");

  generic_interface(time);

  time.def(
          "__init__",
          [](Time* t, const DateTime& dt) { new (t) Time{dt.time}; },
          "From :class:`datetime.datetime`, in UTC if it has no timezone")
      .def(
          "__init__",
          [](Time* t, Numeric x) {
            new (t) Time{};
            t->Seconds(x);
          },
          "From :class:`float` seconds from Unix time start")
      .def(
          py::init<std::string>(),
          "From :class:`str` of form \"YYYY-MM-DD hh:mm:ss\", optionally with fractional seconds and an offset from UTC, "
          "e.g., \"2000-01-01 12:00:00.5+09:00\".  Times without an offset are UTC")
      .def_prop_rw(
          "time",
          [](const Time& t) { return DateTime{t}; },
          [](Time& t, const DateTime& dt) { t = dt.time; },
          "The time with its timezone\n\n.. :class:`~datetime.datetime`")
      .def_rw("utc_offset", &Time::utc_offset, "The timezone as an offset from UTC\n\n.. :class:`~datetime.timedelta`")
      .def_prop_rw(
          "sec",
          [](const Time& t) { return t.Seconds(); },
          [](Time& t, Numeric n) { return t.Seconds(n); },
          "Time from Unix start\n\n.. :class:`float`");

  py::implicitly_convertible<DateTime, Time>();
  py::implicitly_convertible<std::string, Time>();
  py::implicitly_convertible<Numeric, Time>();

  auto a1 = py::bind_vector<ArrayOfTime, py::rv_policy::reference_internal>(m, "ArrayOfTime")
                .def_prop_ro(
                    "as_datetime",
                    [](const ArrayOfTime& in) {
                      std::vector<DateTime> out;
                      out.reserve(in.size());
                      for (auto& t : in) out.push_back(DateTime{t});
                      return out;
                    },
                    "Convert to time lists\n\n.. :class:`list[~datetime.datetime]`");

  generic_interface(a1);

  vector_interface(a1);

  auto a2 = py::bind_vector<ArrayOfArrayOfTime, py::rv_policy::reference_internal>(m, "ArrayOfArrayOfTime");
  generic_interface(a2);
  vector_interface(a2);
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize time\n{}", e.what()));
}
}  // namespace Python
