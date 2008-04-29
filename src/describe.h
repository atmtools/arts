/* Copyright (C) 2003-2008 Stefan Buehler <sbuehler@ltu.se>

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

/*!
  \file   describe.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Feb 25 15:35:56 2003
  
  \brief  Header file for describe.cc

  For documentation see file describe.cc
*/

#ifndef describe_h
#define describe_h

#include "matpackVII.h"
#include "mystring.h"

string describe( ConstTensor7View x );
string describe( ConstTensor6View x );
string describe( ConstTensor5View x );
string describe( ConstTensor4View x );
string describe( ConstTensor3View x );
string describe( ConstMatrixView x );
string describe( ConstVectorView x );
string describe( const Numeric& x );

#endif // describe_h
