/* Copyright (C) 2001,2002,2003 
   Thomas Kuhn    <tkuhn@uni-bremen.de>
   Stefan Buehler <sbuehler@uni-bremen.de>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

// FIXME: Port doxygen file header and documentation from arts-1-0 version to here. 

#include "matpackI.h"
#include "mystring.h"
#include "exceptions.h"
#include "messages.h"

//! Calculate water vapor absorption.
/*! 
  See arts -d absMPM02_H2O for detailed documentation.

  Allowed options for abs_model:<br>
  "MPM02"          - Calculate lines and continuum.<br>
  "MPM02Lines"     - Calculate only lines.<br>
  "MPM02Continuum" - Calculate only continuum.<br>
  "user"           - Use parameters given by abs_user_parameters,
		     instead of the predefined settings.

  Meaning of abs_user_parameters:<br>
  Only used if abs_model=="user". In that case, abs_user_parameters must have 3 elements:<br>
  1. Continuum scaling factor<br>
  2. Line strength scaling factor<br>
  3. Line broadening scaling factor<br>
  Setting all scaling factors to 1 gives the same behavior as abs_model=="MPM02".

  \retval abs Absorption coefficients [1/m], dimension: [ f_grid, abs_p (=abs_t) ].

  \param f_grid Frequency grid [Hz].
  \param abs_p  List of pressures [Pa].
  \param abs_t  List of temperatures [K]. Must have same length as abs_p!
  \param abs_vmr List of volume mixing ratios [absolute number]. Must have same length as abs_p!
  \param abs_model String specifying the model to use. 
  \param abs_user_parameters Parameters for abs_model="user".

  \author Thomas Kuhn, Stefan Buehler (new interface)
  \date   2002-05-06,  2003-11-16
*/
void absMPM02_H2O(// WS Output:
                  Matrix& abs,
                  // WS Input:
                  const Vector& f_grid,
                  const Vector& abs_p,
                  const Vector& abs_t,
                  const Vector& abs_vmr,
                  const String& abs_model,
                  const Vector& abs_user_parameters)
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 model (J. Liebe and G. A. Hufford and M. G. Cotton,
  // "Propagation modeling of moist air and suspended water/ice
  // particles at frequencies below 1000 GHz",
  // AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21)
  const Numeric CC_MPM02 = 1.00000;
  const Numeric CL_MPM02 = 1.00000;
  const Numeric CW_MPM02 = 1.00000;
  // ---------------------------------------------------------------------------------------

  // Vector abs_user_parameters must generally be empty. But if abs_model is
  // "user", then it must have exactly 3 elements.
  if ( abs_model == "user" )
    {
      if ( abs_user_parameters.nelem() != 3 )
        throw runtime_error("Vector abs_user_parameters must have 3 elements.");
    }
  else
    {
      if ( abs_user_parameters.nelem() != 0 )
        throw runtime_error("Vector abs_user_parameters must be empty.");
    }   

  // select the parameter set (!!model dominates values!!):
  Numeric CC, CL, CW;
  // number of lines of Liebe line catalog (0-33 lines, 34 cont. pseudo line)
  Index i_first = 0;
  Index i_last  = 34;
  if ( abs_model == "MPM02" )
    {
      CC      = CC_MPM02;
      CL      = CL_MPM02;
      CW      = CW_MPM02;
      i_first = 0;
      i_last  = 34;
    }
  else if ( abs_model == "MPM02Lines" )
    {
      CC      = 0.000;
      CL      = CL_MPM02;
      CW      = CW_MPM02;
      i_first = 0;
      i_last  = 33;
    }
  else if ( abs_model == "MPM02Continuum" )
    {
      CC      = CC_MPM02;
      CL      = 0.000;
      CW      = 0.000;
      i_first = 34;
      i_last  = 34;
    }
  else if ( abs_model == "user" )
    {
      CC      = abs_user_parameters[1];
      CL      = abs_user_parameters[2];
      CW      = abs_user_parameters[3];
      i_first = 0;
      i_last  = 34;
    }
  else
    {
      ostringstream os;
      os << "Wrong abs_model string given.\n"
	 << "Valid abs_models are: 'MPM02', 'MPM02Lines', 'MPM02Continuum', and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out3  << "  H2O-MPM02: (abs_model=" << abs_model << ") parameter values in use:\n" 
	<< "  CC = " << CC << "\n"
	<< "  CL = " << CL << "\n"
	<< "  CW = " << CW << "\n";
  
  
  const Index n_p = abs_p.nelem();	// Number of pressure levels
  const Index n_f = f_grid.nelem();	// Number of frequencies

  // Check that dimensions of abs_p, abs_t, and abs_vmr agree:
  if ( n_p != abs_t.nelem() ||
       n_p != abs_vmr.nelem() )
    {
      ostringstream os;
      os << "The vectors abs_p, abs_t, and abs_vmr must all have the same length!\n"
	 << "Actual lengths: "
         << n_p << ", "
         << abs_t.nelem() << ", "
         << abs_vmr.nelem() << ".";
      throw runtime_error(os.str());
    }

  // Make abs the right dimension and inititalize to zero:
  abs.resize(n_f,n_p);

  // FIXME: Now should come the real work...
  
  // Assign a dummy value to abs for testing:
  abs = 99;

  return;
}
