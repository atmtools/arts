/* Copyright (C) 2002-2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>

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



/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
   \file   physics_funcs.h
   \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   \date   2002-05-08

   This file contains declerations of functions of physical character.
*/


#ifndef physics_h
#define physics_h

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts.h"
#include "complex.h"          
#include "matpackI.h"



/*===========================================================================
  === Functions in physics_funcs.h
  ===========================================================================*/

Numeric dinvplanckdI(
        const Numeric&  i,
        const Numeric&  f );

void fresnel(
             Complex&   Rv,
             Complex&   Rh,
       const Complex&   n1,
       const Complex&   n2,
       const Numeric&   theta );

Numeric invplanck(
        const Numeric&  i,
        const Numeric&  f );

Numeric invrayjean(
        const Numeric&  i,
        const Numeric&  f );

Numeric number_density(  
        const Numeric&   p,
        const Numeric&   t );

Numeric planck( 
        const Numeric&   f, 
        const Numeric&   t );

Numeric rayjean(
        const Numeric&  f,
        const Numeric&  t );

#endif // physics_h
