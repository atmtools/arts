/* Copyright (C) 2003 Oliver Lemke <olemke@uni-bremen.de>

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


////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   bofstream.cc
  \author Oliver Lemke <olemke@uni-bremen.de>
  \date   2003-01-21

  \brief This file contains the class implementation of bofstream.

*/

#include <fstream>
#include "bofstream.h"
using namespace std;


template<typename T>
bofstream &operator<< (bofstream &bof, T n)
{
  bof.write ((char *)&n, sizeof (n));

  return (bof);
}

/* Explicit instantiation of output operators */
template
bofstream &operator<< <double> (bofstream &bof, double n);

template
bofstream &operator<< <float> (bofstream &bof, float n);

template
bofstream &operator<< <int> (bofstream &bof, int n);

template
bofstream &operator<< <long> (bofstream &bof, long n);

