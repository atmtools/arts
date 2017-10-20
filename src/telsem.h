/* Copyright (C) 2017
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
  \file   telsem.h

  \brief  This file contains the definition of the TELSEM atlas format.
*/

#ifndef telsem_h
#define telsem_h

#include "matpackIII.h"
#include "mystring.h"


typedef struct
{
    // Stored in file
    Index ndat;
    Index nchan;
    String name;
    Index month;
    Numeric dlat;
    Matrix emis;
    Matrix emis_err;
    Tensor3 correl;
    Vector class1;  // ArrayOfIndex, but stored as Vector for efficiency
    Vector class2;  // ArrayOfIndex, but stored as Vector for efficiency
    Vector cellnum;  // ArrayOfIndex, but stored as Vector for efficiency

    // Derived from file data
    Vector correspondence;  // ArrayOfIndex, but stored as Vector for efficiency
} TelsemAtlas;

typedef Array<TelsemAtlas> ArrayOfTelsemAtlas;

void telsem_calc_correspondence(TelsemAtlas& ta);

Index calc_cellnum(const Numeric lat, const Numeric lon, const TelsemAtlas& ta);

void telsem_read_ascii(std::istream& is, TelsemAtlas& ta);

std::ostream& operator<<(std::ostream& os, const TelsemAtlas& ta);

#endif /* telsem_h */

