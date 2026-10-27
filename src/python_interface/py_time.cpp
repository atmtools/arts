#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/chrono.h>
#include <nanobind/stl/string.h>
#include <python_interface.h>

#include <chrono>
#include <cmath>
#include <optional>

#include "artstime.h"
#include "hpy_arts.h"
#include "hpy_vector.h"

namespace Python {
namespace {
//! The time as a datetime.datetime in its timezone
py::object to_datetime(const Time& t) {
  using namespace std::chrono;

  const year_month_day                   date = t.date();
  const hh_mm_ss<Time::InternalTimeStep> tod  = t.time_of_day();

  const py::module_ datetime = py::module_::import_("datetime");
  const py::object  tz = datetime.attr("timezone")(datetime.attr("timedelta")("minutes"_a = t.utc_offset.count()));

  return datetime.attr("datetime")(static_cast<int>(date.year()),
                                   static_cast<unsigned>(date.month()),
                                   static_cast<unsigned>(date.day()),
                                   tod.hours().count(),
                                   tod.minutes().count(),
                                   tod.seconds().count(),
                                   duration_cast<microseconds>(tod.subseconds()).count(),
                                   "tzinfo"_a = tz);
}

//! The time of a datetime.datetime, in UTC if it has no timezone, or nothing if it is not a datetime.datetime
std::optional<Time> from_datetime(const py::handle& dt) {
  using namespace std::chrono;

  if (not py::isinstance(dt, py::module_::import_("datetime").attr("datetime"))) return std::nullopt;

  const auto field = [&dt](const char* name) { return py::cast<int>(dt.attr(name)); };

  minutes offset{0};
  if (const py::object utc_offset = dt.attr("utcoffset")(); not utc_offset.is_none()) {
    const Numeric sec = py::cast<Numeric>(utc_offset.attr("total_seconds")());
    if (std::fmod(sec, 60.0) != 0.0)
      arts::user_error("The offset from UTC of the datetime must be whole minutes, it is {} s", sec);
    offset = minutes{static_cast<int>(sec / 60.0)};
  }

  return Time{year{field("year")} / field("month") / field("day"),
              duration_cast<Time::InternalTimeStep>(hours{field("hour")} + minutes{field("minute")} +
                                                    seconds{field("second")} + microseconds{field("microsecond")}),
              offset};
}
}  // namespace

void py_time(py::module_& m) try {
  py::class_<Time> time(m, "Time");

  generic_interface(time);

  time.def(
          "__init__",
          [](Time* t, Numeric x) {
            new (t) Time{};
            t->Seconds(x);
          },
          "From :class:`float` seconds from Unix time start")
      .def(py::init<std::string>(),
           "From :class:`str` of form \"YYYY-MM-DD hh:mm:ss\" or \"YYYY-MM-DDThh:mm:ss\", optionally with fractional "
           "seconds and an offset from UTC, \"Z\", \"+hh:mm\" or \"-hh:mm\", e.g., \"2000-01-01T12:00:00.5+09:00\".  "
           "Times without an offset are UTC")
      .def(
          "__init__",
          [](Time* t, const py::object& dt) {
            const std::optional<Time> x = from_datetime(dt);
            if (not x) throw py::next_overload();
            new (t) Time{*x};
          },
          "dt"_a,
          "From :class:`datetime.datetime`, in UTC if it has no timezone")
      .def_prop_rw(
          "time",
          [](const Time& t) { return to_datetime(t); },
          [](Time& t, const py::object& dt) {
            const std::optional<Time> x = from_datetime(dt);
            if (not x) throw py::type_error("The time must be a datetime.datetime");
            t = *x;
          },
          "The time with its timezone\n\n.. :class:`~datetime.datetime`")
      .def_rw("utc_offset", &Time::utc_offset, "The timezone as an offset from UTC\n\n.. :class:`~datetime.timedelta`")
      .def_prop_rw(
          "sec",
          [](const Time& t) { return t.Seconds(); },
          [](Time& t, Numeric n) { return t.Seconds(n); },
          "Time from Unix start\n\n.. :class:`float`");

  py::implicitly_convertible<std::string, Time>();
  py::implicitly_convertible<Numeric, Time>();

  auto a1 = py::bind_vector<ArrayOfTime, py::rv_policy::reference_internal>(m, "ArrayOfTime")
                .def_prop_ro(
                    "as_datetime",
                    [](const ArrayOfTime& in) {
                      py::list out;
                      for (auto& t : in) out.append(to_datetime(t));
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
