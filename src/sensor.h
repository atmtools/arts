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

#include "arts.h"
#include "matpackI.h"
#include "matpackII.h"
#include "interpolation.h"
//#include "agenda_class.h"
//#include "array.h"
#include "math_funcs.h"
//#include "mystring.h"


/*===========================================================================
  === Functions from sensor.cc
  ===========================================================================*/

void antenna_transfer_matrix(
           Sparse&      H,
      ConstVectorView   m_za,
      ConstMatrixView   srm,
      ConstVectorView   x_f );

void mixer_transfer_matrix(
              Sparse&   H,
              Vector&   f_mixer,
      ConstVectorView   f_grid,
           const bool   is_upper,
        const Numeric   lo,
      ConstMatrixView   sfrm );

void scale_antenna_diagram(
           VectorView   sc,
      ConstMatrixView   srm,
       const Numeric&   f_ref,
       const Numeric&   f_new );

void sensor_integration_vector(
           VectorView   h,
      ConstVectorView   f,
      ConstVectorView   x_ftot,
      ConstVectorView   x_g );

void sensor_summation_vector(
           VectorView   h,
        const Numeric   f,
      ConstVectorView   f_grid,
	    const Numeric   lo,
      ConstMatrixView   sfrm );

void spectrometer_transfer_matrix(
           Sparse&      H,
      ConstMatrixView   srm,
      ConstVectorView   x_s,
      ConstVectorView   x_f );

#endif  // sensor_h
