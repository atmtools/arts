/* Copyright (C) 2000, 2001, 2002 Axel von Engeln <engeln@uni-bremen.de>
                                  Stefan Buehler  <sbuehler@uni-bremen.de>

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
  \file   m_lineshapes.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Thu May 30 10:19:45 2002
  
  \brief  Lineshape Workspace methods.
  
*/

#include <cmath>
#include "auto_md.h"
#include "make_vector.h"
#include "messages.h"


//! The Lorentz line shape.
/*! 
  See online documentation for details.
*/
void elsLorentz(// WS Output:
               Vector&        els,
               // WS Input:
               const Numeric& ls_gamma,
               const Vector&  els_f_grid)
{
  // PI:
  extern const Numeric PI;

  Index nf = els_f_grid.nelem();

  // Resize will do nothing if the size is already correct.
  els.resize(nf);

  Numeric gamma2 = ls_gamma * ls_gamma;
  Numeric fac = ls_gamma/PI;

  for ( Index i=0; i<nf; ++i )
    {
      Numeric deltaf = els_f_grid[i];
      els[i] =  fac / ( deltaf*deltaf + gamma2 );
    }
}

//! The Doppler line shape
/*!
 See online documentation for details.
*/

void elsDoppler(// WS Output:
		Vector&       els,
                // WS Output:
                const Numeric& ls_sigma,
                const Vector& els_f_grid)
{
  //PI:
  extern const Numeric PI;

  Index nf = els_f_grid.nelem();

  // Resize will do nothing if the size is already correct.
  els.resize(nf);

  Numeric sigma2 = ls_sigma * ls_sigma;
  Numeric fac = 1.0 / (ls_sigma * sqrt(PI));

  for ( Index i=0; i<nf; ++i )
    {
      Numeric deltaf = els_f_grid[i];
      els[i] =  fac * exp( - deltaf * deltaf / sigma2 );
    }
}

void lsWithCutoffAdd(// WS Output:
                     Vector&         ls,
                     Vector&         els,
                     Vector&         els_f_grid,
                     // WS Input:
                     const Agenda&   els_agenda,
                     const Numeric&  ls_cutoff,
                     const Numeric&  ls_f0,
                     const Numeric&  ls_gamma,
                     const Numeric&  ls_sigma,
                     const Vector&   f_grid)
{
  // The frequency grid f_grid must be sorted, otherwise this method
  // will not work correctly. We have decided to make this a general
  // demand. Absorption routines should check this on a higher level.

  // Checks of input quantities:
  // Move these to a higher level and replace by asserts?
  if ( ls_cutoff<0 )
    {
      throw runtime_error("ls_cutoff must be >= 0.");
    }

  if ( ls_gamma<0 )
    {
      throw runtime_error("ls_gamma must be >= 0.");
    }

  if ( ls_sigma<0 )
    {
      throw runtime_error("ls_sigma must be >= 0.");
    }

  if ( ls.nelem() != f_grid.nelem() )
    {
      throw runtime_error("Dimension of ls must be same as f_grid.");
    }

  // The upper and lower cutoff frequencies:
  Numeric l_cut = ls_f0 - ls_cutoff;
  Numeric u_cut = ls_f0 + ls_cutoff;

  // The upper and lower end of the calculation frequencies:
  Numeric l_calc = f_grid[0];
  Numeric u_calc = f_grid[f_grid.nelem()-1];

  // Check whether we are completely outside:
  if ( l_calc > u_cut ) return;
  if ( u_calc < l_cut ) return;

  // Find out, where l_cut and u_cut are, relative to f_grid:
  Index l_i = 0;
  while ( f_grid[l_i] < l_cut ) ++l_i;
  Index u_i = f_grid.nelem()-1;
  while ( f_grid[u_i] > u_cut ) --u_i;

  // Make els_f_grid the right size. We need one additional element for
  // the frequency at the cutoff.
  Index n = u_i - l_i + 1;
  els_f_grid.resize(n+1);
  
  // Copy that range from f_grid (last element of els_f_grid stays
  // free): 
  els_f_grid[Range(0,n)] = f_grid[Range(l_i,n)];

  // Subtract center frequency (This subtracts also from the last
  // element, but this doesn't matter, since we'll overwrite it
  // anyway): 
  els_f_grid -= ls_f0;

  // Put cutoff frequency into the last element:
  els_f_grid[n] = ls_cutoff;

  //----------------------------------------------------------------------
  // From here on we deal with els_agenda
  //----------------------------------------------------------------------

  // Check that the agenda takes the right kind of input:

  if ( !els_agenda.is_input(els_f_grid_) )
    {
      throw runtime_error("The agenda els_agenda must use els_f_grid as an input.");
    }

  if ( !els_agenda.is_input(ls_gamma_) )
    {
      static bool have_reported_gamma = false;
      if ( !have_reported_gamma )
	{
	  out0 << "  Warning: The agenda els_agenda is not using ls_f_gamma.\n"
	       << "  This may be ok, for example if you use a pure Doppler shape.\n";
	  have_reported_gamma = true;
	}
    }
  
  if ( !els_agenda.is_input(ls_sigma_) )
    {
      static bool have_reported_sigma = false;
      if ( !have_reported_sigma )
	{
	  out0 << "  Warning: The agenda els_agenda is not using ls_f_sigma.\n"
	       << "  This may be ok, for example if you use a pure Lorentz shape.\n";
	  have_reported_sigma = true;
	}
    }
  
  // Check that the agenda gives the right kind of output:

  if ( !els_agenda.is_output(els_) )
    {
      throw runtime_error("The agenda els_agenda must generate els as output.");
    }

  // Everything seems ok. Lets call the agenda:
  els_agenda.execute();

  // After this, the new line shape should be in the WSV els!

  // Subtract the value at the cutoff (stored in the last element):
  els -= els[n];

  // This means we also subtract the value at the cutoff from itself,
  // so the last element should be zero afterwards. This is ok, sice
  // we don't need it anymore.

  // Add the new shape to the right place in ls:
  ls[Range(l_i,n)] += els[Range(0,n)];

  // Done.
}
