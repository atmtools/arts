/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                      Stefan Buehler   <sbuehler@uni-bremen.de>

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



/*****************************************************************************
 ***  File description 
 *****************************************************************************/

/*!
   \file   physics_funcs.h
   \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
   \date   2002-05-08

   This file contains declerations of functions of physical character.
*/


#ifndef physics_h
#define physics_h

#include "arts.h"
#include "matpackI.h"



void invplanck (
              VectorView   y,
         ConstVectorView   f,
	 ConstVectorView   za );

void invrayjean (
              VectorView   y,
         ConstVectorView   f,
         ConstVectorView   za );

Numeric number_density (
       const Numeric   p,
       const Numeric   t );

Vector number_density (
       ConstVectorView    p,
       ConstVectorView    t );

void planck (
             MatrixView   B, 
        ConstVectorView   f,
        ConstVectorView   t );

void planck (
              VectorView   B, 
         ConstVectorView   f,
         const Numeric     t );


#endif // physics_h
