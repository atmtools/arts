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
  \file   m_sensor.cc
  \author Mattias Ekström <ekstrom@rss.chalmers.se>
  \date   2003-02-13

  \brief  Workspace functions related to sensor modelling variables.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "math_funcs.h"
#include "messages.h"
#include "ppath.h"
#include "special_interp.h"
#include "xml_io.h"
#include "sensor.h"

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! SensorIntegrationVector
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date   2003-02-13
*/
/*
void SensorIntegrationVector(
        // WS Generic Output:
              Vector&   h_out,
        // WS Generic Output Names:
        const String&   h_out_name,
        // WS Generic Input:
        const Vector&   f_values,
        const Vector&   f_grid,
        const Vector&   g_grid,
        // WS Generic Input Names:
        const String&   f_values_name,
        const String&   f_grid_name,
        const String&   g_grid_name )
{
  //Call sensor_integration_vector in sensor.cc
  h_out.resize(g_grid.nelem());
  sensor_integration_vector( h_out, f_values, f_grid, g_grid);

}
*/

//! AntennaTransferMatrix
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date   2003-03-06
*/
void AntennaTransferMatrix(// WS Generic Output:
                                 Matrix&    Hb,
                           // WS Generic Output Names:
                           const String&    Hb_name,
                           // WS Generic Input:
                           const Vector&    mblock_za_grid,
                           const Matrix&    a,
                           const Vector&    a_grid,
                           const Vector&    f_grid,
                           // WS Generic Input Names:
                           const String&    mblock_za_grid_name,
                           const String&    a_name,
                           const String&    a_grid_name,
                           const String&    f_grid_name)
{
  //Call antenna_transfer_matrix in sensor.cc
  Hb.resize( f_grid.nelem(), mblock_za_grid.nelem() * f_grid.nelem() );
  antenna_transfer_matrix( Hb, mblock_za_grid, a, a_grid, f_grid );

}
