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

#include <thread>

#include "matpackI.h"
#include "artstime.h"
#include "messages.h"
#include "sorting.h"


void timeNow(Time& time, const Verbosity&)
{
  time = Time();
}


void Duration(Numeric& duration, const Time& t0, const Time& t1, const Verbosity&)
{
  duration = std::chrono::duration_cast<TimeStep>(t1 - t0).count();
}


void Sleep(const Numeric& duration, const Verbosity&)
{
  std::this_thread::sleep_for(TimeStep(duration));
}


void timeSleep(const Time& time, const Verbosity&)
{
  std::this_thread::sleep_until(time.Data());
}


void LocalTimeOffset(Numeric& dt, const Verbosity&)
{
  const Time time;
  const Time localtime = Time(time.toStruct());
  const Time gmtime = Time(time.toGMTStruct());
  dt = std::chrono::duration_cast<TimeStep>(localtime - gmtime).count();
}


void timeOffset(Time& time, const Numeric& offset, const Verbosity&)
{
  time += TimeStep(offset);
}


void time_gridOffset(ArrayOfTime& time_grid, const Numeric& offset, const Verbosity& verbosity)
{
  for (Time& time: time_grid)
    timeOffset(time, offset, verbosity);
}
