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

#include "auto_md.h"
#include "make_vector.h"
#include "messages.h"

/*! The Lorentz line shape. This is a quick and dirty implementation.

    \param ls         Output: The shape function.
    \param ls_f0      Line center frequency.
    \param ls_gamma   The pressure broadening parameter.
    \param ls_f_grid  The frequency grid.

    \author Stefan Buehler 
    \date 2000-06-16 */
void lsLorentz(// WS Output:
               Vector&        ls,
               // WS Input:
               const Numeric& ls_f0,
               const Numeric& ls_gamma,
               const Vector&  ls_f_grid)
{
  // PI:
  extern const Numeric PI;

  Index nf = ls_f_grid.nelem();

//   if (!is_size(ls,nf))
//     {
//       throw runtime_error("The lineshape vector ls must have the same size as\n"
// 			  "the frequency grid ls_f_grid.");
//     }

  // Resize will do nothing if the size is already correct.
  ls.resize(nf);

  Numeric gamma2 = ls_gamma * ls_gamma;
  Numeric fac = ls_gamma/PI;

  for ( Index i=0; i<nf; ++i )
    {
      Numeric deltaf = ls_f_grid[i]-ls_f0;
      ls[i] =  fac / ( deltaf*deltaf + gamma2 );
    }
}

void lsWithCutoffAdd(// WS Output:
                     Vector&         ls,
                     Vector&         ls_f_grid,
                     // WS Input:
                     const Agenda&   elem_ls_agenda,
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

  // Find out, where these adjusted l_cut and u_cut are, relative to
  // f_grid:  
  Index l_i = 0;
  while ( f_grid[l_i] < l_cut ) ++l_i;
  Index u_i = f_grid.nelem()-1;
  while ( f_grid[u_i] > u_cut ) --u_i;

  // Make ls_f_grid the right size. We need one additional element for
  // the frequency at the cutoff.
  Index n = u_i - l_i + 1;
  ls_f_grid.resize(n+1);
  
  // Copy that range from f_grid (last element of ls_f_grid stays
  // free): 
  ls_f_grid[Range(0,n)] = f_grid[Range(l_i,n)];

  // Put frequency at the cutoff in the last element:
  ls_f_grid[n] = u_cut;

  //----------------------------------------------------------------------
  // From here on we deal with elem_ls_agenda
  //----------------------------------------------------------------------

  // Check that the agenda takes the right kind of input:

  if ( !elem_ls_agenda.is_input(ls_f0_) )
    {
      throw runtime_error("The agenda elem_ls_agenda must use ls_f0 as an input.");
    }

  if ( !elem_ls_agenda.is_input(ls_f_grid_) )
    {
      throw runtime_error("The agenda elem_ls_agenda must use ls_f_grid as an input.");
    }

  if ( !elem_ls_agenda.is_input(ls_gamma_) )
    {
      static bool have_reported_gamma = false;
      if ( !have_reported_gamma )
	{
	  out0 << "  Warning: The agenda elem_ls_agenda is not using ls_f_gamma.\n"
	       << "  This may be ok, for example if you use a pure Doppler shape.\n";
	  have_reported_gamma = true;
	}
    }
  
  if ( !elem_ls_agenda.is_input(ls_sigma_) )
    {
      static bool have_reported_sigma = false;
      if ( !have_reported_sigma )
	{
	  out0 << "  Warning: The agenda elem_ls_agenda is not using ls_f_sigma.\n"
	       << "  This may be ok, for example if you use a pure Lorentz shape.\n";
	  have_reported_sigma = true;
	}
    }
  
  // Check that the agenda gives the right kind of output:

  if ( !elem_ls_agenda.is_output(ls_) )
    {
      throw runtime_error("The agenda elem_ls_agenda must generate ls as output.");
    }

  // Everything seems ok. Lets call the agenda:
  elem_ls_agenda.execute();

  // After this, the new line shape should be in 
  

}
