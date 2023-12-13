/* Copyright (C) 2003-2012 Oliver Lemke  <olemke@core-dump.info>

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
   USA.
   */

////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
  \file   sorting.h

  Contains sorting routines.

  \author Oliver Lemke
  \date 2003-08-20
  */

#ifndef sorting_h
#define sorting_h

#include <algorithm>
#include <functional>

#include "array.h"
#include "matpack.h"


/** get_sorted_indexes
 *
 * The output array contains the sorted indexes of the input data.
 *
 * The following member functions must be defined for the input data type:
 * T.begin, T.end, T.operator[], T.operator<
 *
 * Furthermore, it must provide the type T::const_iterator.
 *
 * \param sorted  Output array with sorted indexes
 * \param data    Data to sort
 *
 * \author Oliver Lemke <olemke@core-dump.info>
 * \date   2003-08-20
 */
template <typename T>
void get_sorted_indexes(ArrayOfIndex& sorted, const T& data) {
  sorted.resize(0);

  Index i = 0;
  for (typename T::const_iterator it = data.begin(); it != data.end(); ++it) {
    sorted.push_back(i);
    i++;
  }

  sort(sorted.begin(), sorted.end(), [&data](const Index a, const Index b) {
    return data[a] < data[b];
  });
}

#endif /* sorting_h */
