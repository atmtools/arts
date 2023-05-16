/**
 * @file   m_artstime.cc
 * @author Richard Larsson
 * @date   2020-04-13
 * 
 * @brief  Stuff related to manipulating time
 */

#include <thread>

#include "matpack_data.h"
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
  std::this_thread::sleep_until(time.time);
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
