/* Copyright (C) 2008 Stefan Buehler <sbuehler(at)ltu.se>
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

/*!
  \file   arts_omp.cc
  \author Stefan Buehler <sbuehler(at)ltu.se>
  \date   Thu Jan 31 10:04:57 2008
  
  \brief  Helper functions for OpenMP
  
  This file contains wrapper functions for standard OMP functions,
  that work with and without OMP support. This saves the use of \#ifdef
  statements around omp functions in the ARTS main code.

  All functions start with arts_omp. Otherwise, names are the same as
  the standard OMP function names.
*/

#include <iostream>
using namespace std;
#include "arts_omp.h"


//! Wrapper for omp_get_max_threads.
/*! 
  This wrapper works with and without OMP support.

  \return Maximum number of OMP threads, or 1 without OMP.
*/
int arts_omp_get_max_threads()
{

#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
#else
  int max_threads = 1;
#endif

  return max_threads;
}


//! Wrapper for omp_in_parallel.
/*! 
  This wrapper works with and without OMP support.

  \return Returns true if the current region is running parallelized.
*/
bool arts_omp_in_parallel()
{

#ifdef _OPENMP
#ifdef __ICC
  return (omp_get_max_threads() > 1);
#else
  return omp_in_parallel();
#endif
#else
  return false;
#endif
}


//! Wrapper for omp_get_thread_num.
/*! 
  This wrapper works with and without OMP support.

  \return ID number of the current thread, or 0 without OMP.
*/
int arts_omp_get_thread_num()
{

#ifdef _OPENMP
  int thread_num = omp_get_thread_num();
#else
  int thread_num = 0;
#endif
  
  return thread_num;
}

