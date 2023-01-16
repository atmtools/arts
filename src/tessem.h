/* Copyright (C) 2016
   Oliver Lemke <olemke@core-dump.info>

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

/*!
  \file   tessem.h

  \brief  This file contains functions that are adapted from TESSEM
  code which is used to calculate surface emissivity.
*/

#ifndef tessem_h
#define tessem_h

#include <fstream>
#include "matpack_data.h"

struct TessemNN {
  Index nb_inputs;
  Index nb_outputs;
  Index nb_cache;
  Vector b1;
  Vector b2;
  Matrix w1;
  Matrix w2;
  Vector x_min;
  Vector x_max;
  Vector y_min;
  Vector y_max;

  friend std::ostream& operator<<(std::ostream& os, const TessemNN&) {return os;}
};

void tessem_read_ascii(std::ifstream& is, TessemNN& net);

void tessem_prop_nn(VectorView ny, const TessemNN& net, ConstVectorView nx);

#endif /* tessem_h */
