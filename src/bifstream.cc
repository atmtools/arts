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
  \file   bifstream.cc
  \author Oliver Lemke <olemke@uni-bremen.de>
  \date   2003-01-23

  \brief This file contains the class implementation of bifstream.

*/

#include <fstream>
#include "bifstream.h"
using namespace std;


template<typename T>
bifstream &operator>> (bifstream &bif, T &n)
{
  bif.read ((char *)&n, sizeof (n));

  return (bif);
}

/* Explicit instantiation of output operators */
template
bifstream &operator>> <double> (bifstream &bif, double &n);

template
bifstream &operator>> <float> (bifstream &bif, float &n);

template
bifstream &operator>> <int> (bifstream &bif, int &n);

template
bifstream &operator>> <long> (bifstream &bif, long &n);

