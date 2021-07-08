/* Copyright (C) 2020
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/**
 * @file   artstime.cc
 * @author Richard Larsson
 * @date   2020-04-13
 * 
 * @brief  Stuff related to time
*/

#include <ctime>
#include <stdlib.h>

#include "artstime.h"
#include "constants.h"

Time::Time(const String& t) {
  auto s = std::istringstream(t);
  s >> *this;
}

TimeStep time_stepper_selection(const String& time_step)
{
  std::istringstream x(time_step);
  Index length;
  String type;
  x >> length >> type;
  type.tolower();
  
  const Options::TimeStep t = Options::toTimeStep(type);
  check_enum_error(t, "bad time step: ", time_step);
  
  switch (t) {
    case Options::TimeStep::hour:
    case Options::TimeStep::hours:
    case Options::TimeStep::h:
      return TimeStep(std::chrono::hours(length ));
    case Options::TimeStep::minute:
    case Options::TimeStep::minutes:
    case Options::TimeStep::min:
      return TimeStep(std::chrono::minutes(length));
    case Options::TimeStep::second:
    case Options::TimeStep::seconds:
    case Options::TimeStep::s:
      return TimeStep(std::chrono::seconds(length));
    case Options::TimeStep::FINAL: {
      /* Leave empty and last */
    }
  }
  return {};
}

Time next_even(const Time& t, const TimeStep& dt)
{
  auto dt_internal = std::chrono::duration_cast<Time::InternalTimeStep>(dt);
  
  if (dt_internal.count())
    return t + dt_internal - t.EpochTime() % dt_internal;
  else 
    return t;
}

ArrayOfIndex time_steps(const ArrayOfTime& times, const TimeStep& DT)
{
  Index N = times.nelem();
  ARTS_USER_ERROR_IF (N < 2,
    "Can only find time steps for 2-long or longer time grids");
  
  // algorithm only works with absolute times
  const TimeStep dt = std::chrono::abs(DT);
  
  ArrayOfIndex time_steps{0};
  
  if (N > 2) {
    Time tupp = next_even(times[0], dt);
    time_steps.emplace_back(1);
    
    for (; time_steps.back() < N; time_steps.back()++) {
      if (not (times[time_steps.back()] < tupp)) {
        tupp = next_even(tupp, dt);
        time_steps.emplace_back(time_steps.back());
      }
    }
  } else {
    time_steps.emplace_back(N);
  }
  return time_steps;
}

std::ostream& operator<<(std::ostream& os, const Time& t)
{ 
  // FIXME: C++20 has much better calendar handling
  std::tm x = t.toStruct();
  
  // Deal with seconds
  char sec[2+1+9+100];
  Numeric seconds = Numeric(x.tm_sec) + t.PartOfSecond();
  sprintf(sec, "%.9lf", seconds);
  
  // Print based on std::tm specs
  return os << 1900 + x.tm_year << '-' << std::setfill('0') << std::setw(2)
            << 1 + x.tm_mon << '-' << std::setfill('0') << std::setw(2)
            << x.tm_mday << ' ' << std::setfill('0') << std::setw(2)
            << x.tm_hour << ':' << std::setfill('0') << std::setw(2)
            << x.tm_min <<':' << std::setfill('0') << std::setw(12) << sec;
}

std::istream& operator>>(std::istream& is, Time& t)
{
  String ymd, hms;
  is >> ymd >> hms;
  
  ArrayOfString YMD, HMS;
  ymd.split(YMD, "-");
  hms.split(HMS, ":");
  
  ARTS_USER_ERROR_IF (YMD.nelem() not_eq HMS.nelem() and YMD.nelem() not_eq 3,
    "Time stream must look like \"year-month-day hour:min:seconds\"\n"
    "\"year-month-day\"   looks like: ", ymd, '\n',
    "\"hour:min:seconds\" looks like: ", hms);
  
  // FIXME: C++20 has much better calendar handling
  std::tm x;
  x.tm_year = std::stoi(YMD[0]) - 1900;
  x.tm_mon = std::stoi(YMD[1]) - 1;
  x.tm_mday = std::stoi(YMD[2]);
  x.tm_hour = std::stoi(HMS[0]);
  x.tm_min = std::stoi(HMS[1]);
  x.tm_sec = 0;
  x.tm_isdst = -1;
  
  t = Time(x) + TimeStep(std::stod(HMS[2]));
  
  return is;
}

Time mean_time(const ArrayOfTime& ts, Index s, Index E)
{
  Index e=0;
  if (e == -1) e = ts.nelem();
  ARTS_USER_ERROR_IF (e < 0 or e > ts.nelem(),
    "Bad last index, valid options are [-1, ts.nelem()], got: ", E);
  
  ARTS_USER_ERROR_IF (s < 0 or s > ts.nelem(),
    "Bad first index, valid options are [0, ts.nelem()], got: ", s);

  Time::InternalTimeStep dt(0);    
  for (Index i=s+1; i<e; i++)
    dt += (ts[i] - ts[s])  / (e - s);
  return ts[s] + dt;
}

Vector time_vector(const ArrayOfTime& times) {
  Vector t(times.nelem());
  for (Index i=0; i<times.nelem(); i++) t[i] = Numeric(times[i]);
  return t;
}

ArrayOfTime time_vector(const Vector& times) {
  ArrayOfTime t(times.nelem());
  for (Index i=0; i<times.nelem(); i++) t[i].Seconds(times[i]);
  return t;
}

TimeStep median(ArrayOfTimeStep dt) {
  const auto n = dt.size();
  if (n) {
    std::sort(dt.begin(), dt.end());
    if (n % 2) {
      return dt[n / 2];
    } else {
      return (dt[(n-1)/2] + dt[n/2]) / 2;
    }
  } else {
    return TimeStep(0);
  }
}

TimeStep mean(const ArrayOfTimeStep& dt) {
  TimeStep t(0);
  const auto n = dt.size();
  for (std::size_t i=0; i<n; i++) {
    t += dt[i] / n;
  }
  return t;
}
