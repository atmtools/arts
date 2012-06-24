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
extern const Numeric TEMP_0_C;


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void blackbody_radiationPlanck(
           Vector&   blackbody_radiation,
     const Vector&   f,
     const Numeric&  t,
     const Verbosity&)
{
  const Index   n = f.nelem();

  blackbody_radiation.resize(n);

  for( Index i=0; i<n; i++ )
    { blackbody_radiation[i] = planck( f[i], t ); }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void complex_nWaterLiebe93(Matrix&         complex_n,
                           const Vector&   f_grid,
                           const Numeric&  t,
                           const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  chk_if_in_range( "t", t, TEMP_0_C, TEMP_0_C+100 );
  chk_if_in_range( "min of f_grid", min(f_grid), 10e9, 1000e9 );
  chk_if_in_range( "max of f_grid", max(f_grid), 10e9, 1000e9 );

  out2 << "  Sets *complex_n* to model properties of liquid water,\n"
       << "  according to Liebe 1993\n";
  out3 << "     temperature      : " << t << " K.\n";

  const Index   nf = f_grid.nelem();

  complex_n.resize( nf, 2 );

  // Values from epswater93.m (by C. Mätzler), part of Atmlab.
  // The constant e2 is here set to 3.52, which according to Mätzler 
  // corresponds to Liebe 1993.
  const Numeric   theta = 1 - 300 / t;
  const Numeric   e0    = 77.66 - 103.3 * theta;
  const Numeric   e1    = 0.0671 * e0;
  const Numeric   f1    = 20.2 + 146.4 * theta + 316 * theta * theta;
  const Numeric   e2    = 3.52;  
  const Numeric   f2    = 39.8 * f1;

  for( Index iv=0; iv<nf; iv++ )
    { 
      const Complex  ifGHz( 0.0, f_grid[iv]/1e9 );
          
      Complex eps = e2 + (e1-e2) / (Numeric(1.0)-ifGHz/f2) + 
                     (e0-e1) / (Numeric(1.0)-ifGHz/f1);
    
      complex_n(iv,0) = eps.real();
      complex_n(iv,1) = eps.imag();
    }
}


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

  for( Index i=0; i<n; i++ )
    { m(i,0) = planck( f[i], COSMIC_BG_TEMP ); }
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

  for( Index i=0; i<n; i++ )
    { m(i,0) = planck( f[i], t ); }
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
