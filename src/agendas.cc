/* Copyright (C) 2002 Stefan Buehler <sbuehler@uni-bremen.de>

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
  \file   agendas.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Thu Mar 14 08:49:33 2002
  
  \brief  Initialize lookup data for agendas.

  The lookup data mainly contains information on the required output
  and input.  
*/

#include "agenda_record.h"
#include "auto_wsv.h"

// Some #defines and typedefs to make the records better readable:
#define NAME(x) x 
#define DESCRIPTION(x) x
#define OUTPUT   MakeArray<Index>
#define INPUT    MakeArray<Index>

/*! The lookup information for the agendas. */
Array<AgRecord> agenda_data;

void define_agenda_data()
{
  // Initialize to zero, just in case:
  agenda_data.resize(0);

  /*----------------------------------------------------------------------
    Agendas must be put in in alphabetical order. 
    No distinction is made between uppercase and lowercase letters. 
    The sign "_" comes after all letters.
    ----------------------------------------------------------------------*/


 agenda_data.push_back
    (AgRecord
     ( NAME( "abs_vec_agenda" ),
       DESCRIPTION
       (
	"Calculate the absorption vector at each grid point\n"
	"required for solving the raditive transfer equation.\n"
	"\n"
	"This agenda is a set of the following methods:\n"
	"\n"
	"*abs_vec_partCalc* : this method calculates the absorption \n"
	"		     vector for particles\n"
	"*abs_vec_gasCalc*  : this method calculates the absorption\n"
	"		     vector for gaseous species\n"
	"*abs_vecCalc*      : this method sums up the absorption vector\n"
	"                     for particle and gas to give the total \n"
	"		     absorption vector\n"
	" \n"
	"At present the method *abs_vec_gasCalc* is not implemented. For\n"
	"testing the code, we just set the dimensions of *abs_vec_gas*\n"
	"to be same that of *abs_vec_part* with all elements equal to\n"
	"zero.\n"
	"\n"
	"Output :\n"
	"   abs_vec      : The total absorption vector.\n"
	"\n"
	"Input:\n"
	"   abs_vec_spt : Absorption vector for single particle type\n"
	"   pnd_field   : particle number density field\n"
	),
       OUTPUT( abs_vec_ ),
       INPUT(  abs_vec_spt_,
               pnd_field_ )));


 agenda_data.push_back
    (AgRecord
     ( NAME( "convergence_test_agenda" ),
       DESCRIPTION
       (
	"Compute the convergence test.\n"
	"\n"
	"The method *i_fieldIterate* solves the RTE including all Stokes \n"
        "components as well as scattering iterativly. This method requires \n"
        "a convergence test. The user can choose different convergence tests\n"
        "which have to be defined in this agenda.\n"
        "\n"
        "Possible workspace methods are:\n"
        "*convergence_flagAbs*: Calculates the absolute differences, for \n"
        "                       each Stokes component separately.\n"
	),
       OUTPUT( convergence_flag_ ),
       INPUT(  i_field_,
               i_field_old_,
               cloudbox_limits_,
               scat_za_grid_,
	       scat_aa_grid_,
               stokes_dim_)));


 
  agenda_data.push_back
    (AgRecord
     ( NAME( "els_agenda" ),
       DESCRIPTION
       (
	"Compute an elementary lineshape.\n"
	"\n"
	"The elementary lineshape is a simple and symmetric lineshape, for\n"
	"example a Lorentz or Voigt shape. It does not include a cutoff. It\n"
	"also does not include a fore-factor.\n"
	"\n"
	"The method lsWithCutoffAdd uses this agenda to produce a lineshape\n"
	"with cutoff.\n"
	"\n"
	"Not all lineshapes use ls_sigma. (The Lorentz lineshape uses only\n"
	"ls_gamma). \n"
	"\n"
	"Output:    \n"
	"   els        : The lineshape function [1/Hz]  \n"
	"\n"
	"Input:    \n"
	"   ls_gamma   : Pressure broadened line width [Hz].    \n"
	"   ls_sigma   : Doppler broadened line width [Hz]. (Optional)    \n"
	"   els_f_grid : Frequency grid [Hz]."
	),
       OUTPUT( els_ ),
       INPUT(  ls_gamma_,
	       ls_sigma_,
	       els_f_grid_ )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "ext_mat_agenda" ),
       DESCRIPTION
       (
	"Calculate the extinction coefficient matrix at each grid point\n"
	"required for solving the raditive transfer equation.\n"
	"\n"
	"This agenda is a set of the following methods\n"
	"\n"
	"*ext_mat_partCalc* : this method calculates the extinction matrix\n"
	"		      for particles\n"
	"*ext_mat_gasCalc*  : this method calculates the extinction matrix\n"
	"		      for gaseous species\n"
	"*ext_matCalc*      : this method sums up the extinction matrix for\n"
	"		      particle and gas to give the total extinction\n"
	"		      matrix.\n"
	"\n"
	"At present the method *ext_mat_gasCalc* is not implemented. For \n"
	"testing the code, we just set the dimensions of *abs_vec_gas*\n"
	"to be same that of *abs_vec_part* with all elements equal to zero\n"
	"\n"
	"Output :\n"
	"   ext_mat	: The total extinction coefficient matrix.\n"
	"\n"
	"Input:\n"
	"   ext_mat_spt : Extinction coefficient for single particle type\n"
	"   pnd_field   : particle number density field\n"
	),
       OUTPUT(  ext_mat_ ),
       INPUT(  ext_mat_spt_,
	       pnd_field_  )));


  agenda_data.push_back
    (AgRecord
     ( NAME( "main_agenda" ),
       DESCRIPTION
       (
	"The agenda corresponding to the entire controlfile. This is executed\n"
	"when ARTS is run."
	),
       OUTPUT(),
       INPUT()));

 
  agenda_data.push_back
    (AgRecord
     ( NAME( "ppath_step_agenda" ),
       DESCRIPTION
       (
	"Calculation of a propagation path step.\n"
	"\n"
	"A propagation path step is defined as the path between some point to\n"
	"a crossing with either the pressure, latitude or longitude grid, and\n"
	"this agenda performs the calculations to determine such a partial\n"
	"propagation path. The starting point is normally a grid crossing\n"
	"point, but can also be an arbitrary point inside the atmosphere, such\n"
	"as the sensor position. Only points inside the model atmosphere are\n"
	"handled.\n"
	"\n"
	"The communication between this agenda and the calling method is\n"
	"handled by *ppath_step*. That variable is used both as input and\n"
	"output to *ppath_step_agenda*. The agenda gets back *ppath_step*\n" 
	"as returned to the calling method and the last path point hold by\n"
	"the structure is accordingly the starting point for the new \n"
	"calculations. If a total propagation path shall be determined, this\n"
	"agenda is called repeatedly until the starting point of the\n"
	"propagation path is found and *ppath_step* will hold all path\n"
	"steps that together make up *ppath*. The starting point is included\n"
	"in the returned structure.\n"
	"\n"
	"The path is determined by starting at the end point and moving\n"
	"backwards to the starting point. The calculations are initiated by\n"
	"filling *ppath_step* with the practical end point of the path.\n"
	"This is either the position of the sensor (true or hypothetical), or\n"
	"some point at the top of the atmosphere (determined by geometrical\n"
	"calculations starting at the sensor). This initialisation is not\n"
	"handled by *ppath_step_agenda*. All fields of *ppath_step* are set\n"
	"by *ppath_step_agenda*. If the sensor is above the model atmosphere\n"
	"the field *constant* can be initiated by the calling method.\n"
	"Otherwise the field shall be set to negative and it is set to the\n"
	"correct value by *ppath_step* at the first call. This procedure is\n"
	"needed as the path constant changes if refraction is considered, or\n"
	"not, when the sensor is placed inside the atmosphere .\n"
	"\n"
	"The agenda performs only calculations to next crossing of a grid, all\n"
	"other tasks must be performed by the calling method, with one\n"
	"exception. If there is an intersection of a blackbody ground, the\n"
	"calculations stop at this point. This is flagged by setting the\n" 
	"background field of *ppath_step*. Beside this, the calling method\n" 
	"must check if the starting point of the calculations is inside the \n"
	"scattering box or below the ground level, and check if the last point\n"
	"of the path has been reached. The starting point (the end furthest\n"
	"away from the sensor) of a full propagation path can be the top of the\n"
	"atmosphere, a blackbody ground (if *blackbody_ground* = 1) and the\n"
	"cloud box.\n"
	"\n"
	"The *ppath_step_agenda* put in points along the propagation path\n"
	"at all crossings with the grids, tangent points and points of ground\n"
	"reflection. It is also allowed to make agendas that put in additional\n"
	"points to fulfil some criterion, such as a maximum distance along\n"
	"the path between the points. Accordingly, the number of new points of\n"
	"each step can exceed one.\n"
	"\n"
	"For more information read the chapter on propagation paths in the\n"
	"ARTS user guide.\n"
	"\n"
	"Usage:             Called from *ppathCalc*."
        "\n"
        "Output: \n"    
        "   ppath_step    : Variable which includes structure for the \n"
        "                   propagation path step.\n"
        "\n"
        "Input: \n"
        "   ppath_step    : Structure has to be initialized. \n"
        "   atmosphere_dim: Atmospheric dimension (1-3).\n"
        "   p_grid        : Pressure grid.\n"
        "   lat_grid      : Latitude grid. \n"
        "   lon_grid      : Longitude grid. \n"
        "   z_field       : Geometrical altitudes.\n"
        "   r_geoid       : Geoid radius. \n"
        "   z_ground      : Altitude of the ground. \n"
        "   blackbody_ground: Flag (1 if earth is treated as blackbody). \n"
	),
       OUTPUT(ppath_step_),
       INPUT(ppath_step_,
             atmosphere_dim_,
             p_grid_,
             lat_grid_,
             lon_grid_,
             z_field_,
             r_geoid_,
             z_ground_,
             blackbody_ground_)));

  agenda_data.push_back
    (AgRecord
     ( NAME( "rte_agenda" ),
       DESCRIPTION
       (
       "Calculation of monochromatic pencil beam spectra and absorption WFs."
       "\n"
       "Text will be written (PE).\n"
       ""
	),
       OUTPUT(),
       INPUT()));

 agenda_data.push_back
    (AgRecord
     ( NAME( "scat_rte_agenda" ),
       DESCRIPTION
       (
       "Calculate the radiative transfer equation (RTE) with fixed scattering \n "
       "integral.\n"
       "\n"
       "There are different possibilities to compute the RTE with a fixed value \n"
       "for the scattering integral. These are implemented in the workspace \n"
       "methods:"
       "sto_vecGeneral: This method uses the Pade approximation to compute the \n"
       "         matrix exponential function and a LU decomposition method \n"
       "         to compute the matrix inverse. This method should not be \n"
       "         used for the scalar RTE (if *stokes_dim* equals 1), as it\n"
       "         is numerically not efficient in this case.\n" 
       "sto_vecScalar: This method can only be used for the scalar RTE. It \n"
       "         uses the standard exponential function and there is no need \n"
       "         for a LU decomposition. That means that it is much more \n"
       "         efficient for the computation of the scalar RTE. \n"
       "\n"
       "Output and Input:\n"
       "   stokes_vec  : The Stokes vector. \n"
       "\n"
       "Input:\n"
       "   ext_mat     : Extinction coefficient matrix.\n"
       "   abs_vec     : Absorption coefficient vector. \n"
       "   sca_vec     : Scattered field vector. \n"
       "   l_step      : Pathlength through a grid cell/ layer.\n"
       "   planck      : Planck function. \n"
       "   stokes_dim  : Stokes dimension. \n"
       ""
	),
       OUTPUT(stokes_vec_),
       INPUT(ext_mat_,
             abs_vec_,
             sca_vec_,
             l_step_,
             planck_function_,
             stokes_dim_)));

}
