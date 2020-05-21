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
 * @file   artstime.h
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
  explicit Time(std::time_t t) : mtime(std::chrono::system_clock::from_time_t(t)) {}
  explicit Time(std::tm t) : Time(std::mktime(&t)) {}
  
  // Data
  const std::chrono::system_clock::time_point& Data() const {return mtime;}
  
  // Conversions
  std::time_t toTimeT() const {return std::chrono::system_clock::to_time_t(mtime);}
  std::tm toStruct() const {std::time_t x=toTimeT();  std::tm* y = std::localtime(&x); return *y;}
  std::tm toGMTStruct() const {std::time_t x=toTimeT();  std::tm* y = std::gmtime(&x); return *y;}
  TimeStep seconds_into_day() const {std::tm x = toStruct(); return TimeStep(x.tm_hour*3600 + x.tm_min*60 + x.tm_sec + PartOfSecond());}
  InternalTimeStep EpochTime() const {return mtime.time_since_epoch();}
  
  // Operations
  InternalTimeStep operator-(const Time& t) const {return mtime - t.mtime;}
  bool operator<(const Time& t) const {return mtime < t.mtime;}
  bool operator==(const Time& t) const {return mtime == t.mtime;}
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

/** List of times */
using ArrayOfArrayOfTime = Array<ArrayOfTime>;

/** Output for Time */
std::ostream& operator<<(std::ostream& os, const Time& t);

/** Input for Time */
std::istream& operator>>(std::istream& is, Time& t);

/** Debug output for duration */
inline std::ostream& operator<<(std::ostream& os, const TimeStep& dt) {return os << dt.count() << " seconds";}

/** Returns the next time after t with an even time-step
 * 
 * @param[in] t A time
 * @param[in] dt A duration of time
 * @return Next even time, e.g., 14:14:00 with dt as 10 minutes gives 14:20:00
 */
Time next_even(const Time& t, const TimeStep& dt);

/** Finds the index matching demands in a list of times
 * 
 * The first index is 0 and the second index is the start of the first even period of
 * the given stepsize
 * 
 * The last index is times.nelem().  If output has 1 element, no range was found matching the
 * criteria.
 * 
 * @param[in] times Times sorted in ascending order
 * @param[in] step A duration of time
 * @return Starting index of the time-series
 */
ArrayOfIndex time_steps(const ArrayOfTime& time, const String& step);

/** Computes the average time in a list
 * 
 * @param[in] ts A list of time
 * @param[in] s A starting index; valid range [0, ts.nelem())
 * @param[in] e The end+1 index; valid range [-1, ts.nelem()]; -1 is treated as ts.nelem()
 */
Time mean_time(const ArrayOfTime& ts, Index s=0, Index e=-1);

#endif  // ARTSTIME_H
