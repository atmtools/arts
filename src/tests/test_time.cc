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
