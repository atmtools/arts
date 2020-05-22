/* Copyright (C) 2020
 * Richard Larsson <larsson@mps.mpg.de>
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

TimeStep time_stepper_selection(const String& time_step)
{
  std::istringstream x(time_step);
  Index length;
  String type;
  x >> length >> type;
  type.tolower();
  
  if (type == "hour" or type == "hours" or type == "h")
    return TimeStep(std::chrono::hours(length ));
  else if (type == "minute" or type == "minutes" or type == "min")
    return TimeStep(std::chrono::minutes(length));
  else if (type == "second" or type == "seconds" or type == "s")
    return TimeStep(std::chrono::seconds(length));
  else 
    throw std::runtime_error("Bad time step definition");
}

Time next_even(const Time& t, const TimeStep& dt)
{
  auto dt_internal = std::chrono::duration_cast<Time::InternalTimeStep>(dt);
  
  if (dt_internal.count())
    return t + dt_internal - t.EpochTime() % dt_internal;
  else 
    return t;
}

ArrayOfIndex time_steps(const ArrayOfTime& times, const String& step)
{
  Index N = times.nelem();
  if (N < 2) {
    throw std::runtime_error("Can only find time steps for 2-long or longer time grids");
  }
  
  auto dt = time_stepper_selection(step);
  if (dt < decltype(dt)(0)) {
    throw std::runtime_error("Must have positive time steps (or 0 for all times)");
  }
  
  ArrayOfIndex time_steps{0};
  
  if (N > 2) {
    Time tupp = next_even(times[0], dt);
    time_steps.push_back(1);
    
    for (; time_steps.back() < N; time_steps.back()++) {
      if (not (times[time_steps.back()] < tupp)) {
        tupp = next_even(tupp, dt);
        time_steps.push_back(time_steps.back());
      }
    }
  } else {
    time_steps.push_back(N);
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
  
  if (YMD.nelem() not_eq HMS.nelem() and YMD.nelem() not_eq 3)
    throw std::runtime_error("Time stream must look like \"year-month-day hour:min:seconds\"");
  
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

Time mean_time(const ArrayOfTime& ts, Index s, Index e)
{
  if (e == -1)
    e = ts.nelem();
  else if (e < 0 or e > ts.nelem())
    throw std::runtime_error("Bad last index, valid options are [-1, ts.nelem()]");
  
  if (s < 0 or s > ts.nelem())
    throw std::runtime_error("Bad first index, valid options are [0, ts.nelem()]");

  Time::InternalTimeStep dt(0);    
  for (Index i=s+1; i<e; i++)
    dt += (ts[i] - ts[s])  / (e - s);
  return ts[s] + dt;
}

