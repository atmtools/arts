/*!
  \file   arts_omp.h
  \author Stefan Buehler <sbuehler(at)ltu.se>
  \date   Thu Jan 31 10:04:57 2008
  
  \brief  Header file for helper functions for OpenMP
  
  This file contains headers for the wrapper functions for standard
  OMP functions, that work with and without OMP support. This saves
  the use of \#ifdef statements around omp functions in the ARTS main
  code.

  All functions start with arts_omp. Otherwise, names are the same as
  the standard OMP function names.
*/

#ifndef arts_omp_h
#define arts_omp_h

#ifdef _OPENMP
#include <omp.h>
#endif

int arts_omp_get_max_threads();

bool arts_omp_in_parallel();

int arts_omp_get_thread_num();

int arts_omp_get_nested();

void arts_omp_set_nested(int i);

void arts_omp_set_dynamic(int i);

#endif  // arts_omp_h
