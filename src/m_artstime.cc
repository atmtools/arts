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


void timeSet(Time& time, const String& time_str, const Verbosity&) try {
  time = Time{time_str};
} catch(std::runtime_error& e) {
  // We perform some checks on the string's validity
  throw e;
} catch (...) {
  // Still, sring conversions may throw std::invalid_argument or std::out_of_range,
  // we don't care which but just want to warn users about their input
  ARTS_USER_ERROR("Cannot convert ", std::quoted(time_str), " to valid time")
}
