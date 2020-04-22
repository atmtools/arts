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
 * @brief  Stuff related to time in ARTS
*/

#ifndef ARTSTIME_H
#define ARTSTIME_H

#include <chrono>
#include <cmath>

#include "array.h"
#include "mystring.h"

/** A duration of time, 1 full tick should be 1 second */
using TimeStep = std::chrono::duration<Numeric>;

/** Class to handle time in ARTS */
class Time {
private:
  std::chrono::system_clock::time_point mtime;

public:
  using InternalTimeStep = decltype(mtime)::duration;
  
  // Version will be updated when C++20 datetime is available... Version 1 will still assume local time at that time
  Index Version() const noexcept {return 1;}
  
  // Construction
  Time() : mtime(std::chrono::system_clock::now()) {}
  Time(std::time_t t) : mtime(std::chrono::system_clock::from_time_t(t)) {}
  Time(std::tm t) : Time(std::mktime(&t)) {}
  
  // Conversions
  std::time_t toTimeT() const {return std::chrono::system_clock::to_time_t(mtime);}
  std::tm toStruct() const {std::time_t x=toTimeT();  std::tm* y = std::localtime(&x); return *y;}
  std::tm toGMTStruct() const {std::time_t x=toTimeT();  std::tm* y = std::gmtime(&x); return *y;}
  TimeStep seconds_into_day() const {std::tm x = toStruct(); return TimeStep(x.tm_hour*3600 + x.tm_min*60 + x.tm_sec + PartOfSecond());}
  InternalTimeStep EpochTime() const {return mtime.time_since_epoch();}
  
  // Operations
  InternalTimeStep operator-(const Time& t) const {return mtime - t.mtime;}
  bool operator<(const Time& t) const {return mtime < t.mtime;}
  template <typename T, typename R> Time& operator+=(const std::chrono::duration<T, R>& dt) {mtime += std::chrono::duration_cast<InternalTimeStep>(dt); return *this;}
  template <typename T, typename R> Time& operator-=(const std::chrono::duration<T, R>& dt) {mtime -= std::chrono::duration_cast<InternalTimeStep>(dt); return *this;}
  template <typename T, typename R> Time operator+(const std::chrono::duration<T, R>& dt) const {return (Time(*this) += dt);}
  template <typename T, typename R> Time operator-(const std::chrono::duration<T, R>& dt) const {return (Time(*this) -= dt);}
  
  // helpers
  Numeric Seconds() const {return std::chrono::duration_cast<TimeStep>(mtime.time_since_epoch()).count();}
  void Seconds(Numeric x) {operator+=(TimeStep(x - Seconds()));}
  Numeric PartOfSecond() const {return std::fmod(Seconds(), 1.0);}
}; // Time

/** List of times */
using ArrayOfTime = Array<Time>;

/** Output for Time */
std::ostream& operator<<(std::ostream& os, const Time& t);

/** Input for Time */
std::istream& operator>>(std::istream& is, Time& t);

/** Debug output for duration */
inline std::ostream& operator<<(std::ostream& os, const TimeStep& dt) {return os << dt.count() << " seconds";}

/** Finds the index matching demands in a list of times
 *
 * The output makes index(i+1)-index(i) mark the range of a time step.  If start_even
 * is true, the time step begins at the beginning of an even step daily (e.g., 2 hour time
 * step means a range starts at 00:00, 02:00, 04:00, and so on)
 * 
 * The last index is times.nelem().  If output has 1 element, no range was found matching the
 * criteria.
 * 
 * @param[in] times Times sorted in ascending order
 * @param[in] step A duration of time
 * @param[in] start_even Flag for the first time
 * @return Starting index of the time-series
 */
ArrayOfIndex time_steps(const ArrayOfTime& time, const String& step, const bool start_even);

#endif  // ARTSTIME_H
