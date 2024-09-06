/**
 * @file   m_artstime.cc
 * @author Richard Larsson
 * @date   2020-04-13
 * 
 * @brief  Stuff related to manipulating time
 */

#include <stdexcept>
#include <thread>

#include "artstime.h"
#include "debug.h"

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
  ARTS_USER_ERROR("Cannot convert \"{}\" to valid time", time_str)
}
