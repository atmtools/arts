/* Copyright (C) 2003 Mattias Ekström <ekstrom@rss.chalmers.se>

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
  ===  File description
  ===========================================================================*/

/*!
  \file   sensor.h
  \author Mattias Ekström <ekstrom@rss.chalmers.se>
  \date   2003-02-28

  \brief  Sensor modelling functions.

   This file contains the definition of the functions in sensor.cc 
   that are of interest elsewhere.
*/


#ifndef sensor_h
#define sensor_h


#include "agenda_class.h"
#include "array.h"
#include "arts.h"
#include "interpolation.h"
#include "matpackI.h"
#include "mystring.h"


/*===========================================================================
  === Functions from sensor.cc
  ===========================================================================*/

void sensor_integration_vector(
           VectorView   h,
      ConstVectorView   f,
      ConstVectorView   x_ftot,
      ConstVectorView   x_g );

void antenna_transfer_matrix(
           MatrixView   Hb,
      ConstVectorView   m_za,
      ConstMatrixView   a,
      ConstVectorView   x_a,
      ConstVectorView   x_f );
      
void antenna_diagram_gaussian(
           VectorView   a,
      ConstVectorView   theta_grid,
       const Numeric&   theta );


#endif  // sensor_h
