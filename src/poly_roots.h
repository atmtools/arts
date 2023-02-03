/* Copyright (C) 2002-2012 Oliver Lemke  <olemke@core-dump.info>

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
/**
   \file   poly_roots.h

   Contains the code to determine roots of polynomials.

   \author Oliver Lemke
   \date 2000-03-06
*/

#ifndef poly_roots_h
#define poly_roots_h

#include "arts.h"
#include "matpack_data.h"

int poly_root_solve(Matrix& roots, Vector& coeffs);

#endif /* poly_roots_h */
