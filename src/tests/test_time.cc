/*!
  \file   test_times.cc
  \author Richard Larsson
  \date   2020-04-27
  
  \brief  Test Time Functions
*/

#include "artstime.h"
#include <iostream>

void test01() {
  auto x=ArrayOfTime(40000, Time(4));
  
  // Set the steps 0.5 seconds appard
  for (Size i=0; i<x.size(); i++)
    x[i] += TimeStep(0.001*Numeric(i));
  
  // Find every start of 5 seconds
  auto time_str = "5 seconds";
  auto limits = time_steps(x, time_stepper_selection(time_str));
  
  std::cout << "The " << time_str << " intervals are:\n";
  for (Size i=0; i<limits.size()-2; i++)
    std::cout << 1+i <<": " << x[limits[i]] << " to "<< x[limits[i+1]] << "\n";
  std::cout << limits.size()-1 <<": " << x[limits[limits.size()-2]] << " to "<< x.back() << "\n";
}

int main() {
  test01();
  return 0;
}
