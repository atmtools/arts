/* Copyright (C) 2020 Richard Larsson <ric.larsson@gmail.com>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*!
  \file   test_times.cc
  \author Richard Larsson
  \date   2020-04-27
  
  \brief  Test Time Functions
*/

#include "artstime.h"

void test01() {
  auto x=ArrayOfTime(40000, Time(4));
  
  // Set the steps 0.5 seconds appard
  for (Index i=0; i<x.nelem(); i++)
    x[i] += TimeStep(0.001*Numeric(i));
  
  // Find every start of 5 seconds
  auto time_str = "5 seconds";
  auto limits = time_steps(x, time_stepper_selection(time_str));
  
  std::cout << "The " << time_str << " intervals are:\n";
  for (Index i=0; i<limits.nelem()-2; i++)
    std::cout << 1+i <<": " << x[limits[i]] << " to "<< x[limits[i+1]] << "\n";
  std::cout << limits.nelem()-1 <<": " << x[limits[limits.nelem()-2]] << " to "<< x.back() << "\n";
}

int main() {
  test01();
  return 0;
}
