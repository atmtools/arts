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
 * @file   m_artstime.cc
 * @author Richard Larsson
 * @date   2020-04-13
 * 
 * @brief  Stuff related to manipulating time
 */

#include <stdexcept>
#include <thread>

#include "debug.h"
#include "matpack_data.h"
#include "artstime.h"
#include "sorting.h"


void timeNow(Time& time)
{
  time = Time();
}


void Duration(Numeric& duration, const Time& t0, const Time& t1)
{
  duration = std::chrono::duration_cast<TimeStep>(t1 - t0).count();
}


void Sleep(const Numeric& duration)
{
  std::this_thread::sleep_for(TimeStep(duration));
}


void timeSleep(const Time& time)
{
  std::this_thread::sleep_until(time.time);
}


void LocalTimeOffset(Numeric& dt)
{
  const Time time;
  const Time localtime = Time(time.toStruct());
  const Time gmtime = Time(time.toGMTStruct());
  dt = std::chrono::duration_cast<TimeStep>(localtime - gmtime).count();
}


void timeOffset(Time& time, const Numeric& offset)
{
  time += TimeStep(offset);
}


void time_gridOffset(ArrayOfTime& time_grid, const Numeric& offset)
{
  for (Time& time: time_grid)
    timeOffset(time, offset);
}


void timeSet(Time& time, const String& time_str) try {
  time = Time{time_str};
} catch(std::runtime_error& e) {
  // We perform some checks on the string's validity
  throw e;
} catch (...) {
  // Still, sring conversions may throw std::invalid_argument or std::out_of_range,
  // we don't care which but just want to warn users about their input
  ARTS_USER_ERROR("Cannot convert ", std::quoted(time_str), " to valid time")
}
