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
 * @file   time.h
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
  
  if (type == "hour" or type == "hours" or type == "h") {
    if (length > 24)
      throw std::runtime_error("Max 24 hours allowed");
    return TimeStep(std::chrono::hours(length ));
  }
  else if (type == "minute" or type == "minutes" or type == "min")
    return TimeStep(std::chrono::minutes(length));
  else if (type == "second" or type == "seconds" or type == "s")
    return TimeStep(std::chrono::seconds(length));
  else 
    return TimeStep(0);
}

Index next_larger_from_time(const ArrayOfTime& times, const TimeStep& step, const Index curpos, const Index start_pos)
{
  if (TimeStep(0) == step) {
    return curpos+1;
  } else {
    Numeric intt0, intt1;
    std::modf((times[curpos] - times[start_pos]) / step, &intt0);
    for (Index i=curpos+1; i<times.nelem(); i++) {
      std::modf((times[i] - times[start_pos]) / step, &intt1);
      if (intt1 > intt0)
        return i;
    }
    return times.nelem();
  }
}

Index start_index_from_time(const ArrayOfTime& times, const TimeStep& step, const bool next_even)
{
  if (not next_even) {
    return 0;
  } else if(TimeStep(0) == step) {
    return 0;
  } else {
    if(times.nelem() == 0)
      return 0;
    else if(times.nelem() > 1 and times[1] - times[0] > step)
      return 0;
    else {
      for (Index i=1; i<times.nelem(); i++)
        if (std::fmod(times[i].seconds_into_day().count(), step.count()) < std::fmod(times[i-1].seconds_into_day().count(), step.count()))
          return i;
      return times.nelem();
    }
  }
}

ArrayOfIndex time_steps(const ArrayOfTime& times, const String& step, const bool start_even)
{
  auto dt = time_stepper_selection(step);
  if (dt < decltype(dt)(0)) {
    throw std::runtime_error("Must have positive time steps (or 0 for all times)");
  }
  
  ArrayOfIndex time_steps{start_index_from_time(times, dt, start_even)};
  
  while (time_steps.back() not_eq times.nelem())
    time_steps.push_back(next_larger_from_time(times, dt, time_steps.back(), time_steps.front()));
    
  return time_steps;
}

std::ostream& operator<<(std::ostream& os, const Time& t)
{ 
  // FIXME: C++20 has much better calendar handling
  std::tm x = t.toStruct();
  
  // Deal with subsections for partial seconds
  
  // Deal with seconds;
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
