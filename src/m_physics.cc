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
  ===  File description
  ===========================================================================*/

/*!
  \file   m_physics.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2002-08-20 

  \brief  Workspace methods of physical character.

  This file includes workspace methods for operations that have some
  connection to basic physics. Example of methods are:  <br>
  1. Setting WSV to hold blackbody radiation. <br>
  2. Conversion to brightness temperature.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "physics_funcs.h"

extern const Numeric COSMIC_BG_TEMP;


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixCBR(// WS Output:
               Matrix&   m,
               // WS Input:
               const Index&    stokes_dim,
               // WS Generic Input:
               const Vector&   f,
               const Verbosity&)
{
  const Index n = f.nelem();

  if( n == 0 )
    throw runtime_error( "The given frequency vector is empty." );

  m.resize(n,stokes_dim);
  m = 0;

  planck( m(joker,0), f, COSMIC_BG_TEMP );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixPlanck(// WS Output:
                  Matrix&   m,
                  // WS Input:
                  const Index&    stokes_dim,
                  // WS Generic Input:
                  const Vector&   f,
                  const Numeric&  t,
                  const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  const Index n = f.nelem();

  if( n == 0 )
    throw runtime_error( "The given frequency vector is empty." );

  out2 << "  Setting blackbody radiation for a temperature of " << t << " K.\n";

  m.resize(n,stokes_dim);
  m = 0;

  planck( m(joker,0), f, t );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixUnitIntensity(// WS Output:
                         Matrix&   m,
                         // WS Input:
                         const Index&    stokes_dim,
                         // WS Generic Input:
                         const Vector&   f,
                         const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  const Index n = f.nelem();

  if( n == 0 )
    throw runtime_error( "The given frequency vector is empty." );

  out2 << "  Setting unpolarised radiation with an intensity of 1.\n";

  m.resize(n,stokes_dim);
  m = 0;

  for( Index i=0; i<n; i++ )
    { m(i,0) = 1.0; }
}
