/* Copyright (C) 2012 Stefan Buehler <sbuehler(at)ltu.se>
  
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
