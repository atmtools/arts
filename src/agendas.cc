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
	"This agenda, for example, can a set of the following methods:\n"
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
	"This agenda, for example, can be a set of the following methods\n"
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
     ( NAME( "ground_refl_agenda" ),
       DESCRIPTION
       (
	"Describes the properties of the ground to consider when there is a\n"
	"ground reflection.\n"
	"\n"
	"The ground properties are described by the WSVs *ground_emission*,\n"
	"*ground_los* and *ground_refl_coeffs*.\n"
	"\n"
	"The upwelling radiation from the ground is calculated as the sum of\n"
	"*ground_emission* and the spectra calculated for the directions\n"
	"given by *ground_los*, multiplicated with the weights in\n"
	"*ground_refl_coeffs*. Or (for frequency i): \n"
	"   i_up = i_emission + sum_over_los( W*i_down ) \n"
	"where i_up is the upwelling radiation, i_emission is row i of\n"
	"*ground_emission*, W is the reflection matrix in \n"
	"*ground_refl_coeffs* for the frequency and LOS of concern and \n"
	"i_down is the downwelling radiation for the LOS given in\n"
	"*ground_los*. \n"
	"\n"
	"With other words, the scattering properties of the ground are \n"
	"described by the variables *ground_los* and *ground_refl_coeffs*."
	),
       OUTPUT( ground_emission_, ground_los_, ground_refl_coeffs_  ),
       INPUT(  f_grid_, stokes_dim_, a_gp_p_, a_gp_lat_, a_gp_lon_, a_los_,
	       r_geoid_, z_ground_, t_field_ )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "i_space_agenda" ),
       DESCRIPTION
       (
	"Sets the workspace variable *i_space* to match the assumptions \n"
	"regarding the radiation entering the atmosphere at the start of a \n"
	"propagation path. \n"
	"\n"
	"The main usage of this agenda is to be called from *RteCalc*. \n"
	"\n"
	"If only cosmic background radiation is considered, *i_space* can be\n"
	"set be before *RteCalc* and the agenda only needs to include ????. \n"
	"\n"
	"A function calling this agenda shall set *a_pos* and *a_los* to \n"
	"the position and line-of-sight for which the entering radiation \n"
	"shall be determined. The position and line-of-sight must be known, \n"
	"for example, when radiation from the sun is considered. \n"
	"\n"
	"Usage:   Called from *RteCalc*."
	),
       OUTPUT(  i_space_ ),
       INPUT(  f_grid_, stokes_dim_, a_pos_, a_los_ )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "main_agenda" ),
       DESCRIPTION
       (
	"The agenda corresponding to the entire controlfile. This is\n" 
	"executed when ARTS is run."
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
	"A propagation path step is defined as the path between some point \n"
	"to a crossing with either the pressure, latitude or longitude grid,\n"
	"and this agenda performs the calculations to determine such a \n"
	"partial propagation path. The starting point is normally a grid \n" 
	"crossing point, but can also be an arbitrary point inside the \n" 
	"atmosphere, such as the sensor position. Only points inside the \n"
	"model atmosphere are handled.\n"
	"\n"
	"The communication between this agenda and the calling method is \n"
	"handled by *ppath_step*. That variable is used both as input and \n"
	"output to *ppath_step_agenda*. The agenda gets back *ppath_step* \n" 
	"as returned to the calling method and the last path point hold by \n"
	"the structure is accordingly the starting point for the new \n"
	"calculations. If a total propagation path shall be determined, this\n"
	"agenda is called repeatedly until the starting point of the \n"
	"propagation path is found and *ppath_step* will hold all path \n"
	"steps that together make up *ppath*. The starting point is included\n"
	"in the returned structure. \n"
	"\n"
	"The path is determined by starting at the end point and moving \n"
	"backwards to the starting point. The calculations are initiated by \n"
	"filling *ppath_step* with the practical end point of the path. \n"
	"This is either the position of the sensor (true or hypothetical), \n"
	"or some point at the top of the atmosphere (determined by\n" 
	"geometrical calculations starting at the sensor). This \n" 
	"initialisation is not handled by *ppath_step_agenda*. All fields of\n"
	"*ppath_step* are set by *ppath_step_agenda*. If the sensor is above\n"
	"the model atmosphere the field *constant* can be initiated by the \n" 
	"calling method. Otherwise the field shall be set to negative and it\n"
	"is set to the correct value by *ppath_step* at the first call. This\n"
	"procedure is needed as the path constant changes if refraction is \n"
	"considered, or not, when the sensor is placed inside the\n" 
	"atmosphere.\n"
	"\n"
	"The agenda performs only calculations to next crossing of a grid, \n"
	"all other tasks must be performed by the calling method, with one \n"
	"exception. If there is an intersection of a blackbody ground, the \n"
	"calculations stop at this point. This is flagged by setting the \n" 
	"background field of *ppath_step*. Beside this, the calling method \n" 
	"must check if the starting point of the calculations is inside the \n"
	"scattering box or below the ground level, and check if the last \n"
	"point of the path has been reached. The starting point (the end \n"
	"furthest away from the sensor) of a full propagation path can be \n"
	"the top of the atmosphere, a blackbody ground (if \n"
	"*blackbody_ground* = 1) and the cloud box.\n"
	"\n"
	"The *ppath_step_agenda* put in points along the propagation path \n"
	"at all crossings with the grids, tangent points and points of \n"
	"ground reflection. It is also allowed to make agendas that put in \n"
	"additional points to fulfil some criterion, such as a maximum \n"
	"distance along the path between the points. Accordingly, the \n"
	"number of new points of each step can exceed one.\n"
	"\n"
	"For more information read the chapter on propagation paths in the\n"
	"ARTS user guide.\n"
	"\n"
	"Usage: Called from *PpathCalc* and from functions doing scattering\n"
	"       calculations."
	),
       OUTPUT(ppath_step_),
       INPUT(ppath_step_,
             atmosphere_dim_,
             p_grid_,
             lat_grid_,
             lon_grid_,
             z_field_,
             r_geoid_,
             z_ground_ )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "rte_agenda" ),
       DESCRIPTION
       (
	"Performs monochromatic pencil beam calculations for a single\n"
	"propagation path."
       "\n"
       "More text will be written (PE).\n"
       ""
	),
       OUTPUT( i_rte_ ),
       INPUT(  ppath_ )));

agenda_data.push_back
    (AgRecord
     ( NAME( "scat_mono_agenda" ),
       DESCRIPTION
       (
       "Performs the monochromatic scattering calculation."
       "\n"
       "Normally this agenda consists of three methods: \n"
       "   1. i_fieldSetClearsky or i_fieldSetConst \n"
       "   2. i_fieldIterate (Other solution methods might be implemented\n"
       "      later.) \n"
       "   3. scat_i_Put \n"
       "\n"
       "Output and Input:\n"
       "   scat_i_p, scat_i_lat, scat_i_lon:  Intensity field on the \n"
       "                                      cloudbox boundary. \n"
       "\n"
       "Input:\n"
       "   f_grid:       Frequency grid. \n"
       "   scat_f_index: Frequency index for the ongoing scattering \n"
       "                 calculation. \n"
       ""
	),
       OUTPUT(scat_i_p_,
              scat_i_lat_,
              scat_i_lon_),
       INPUT(scat_i_p_,
             scat_i_lat_,
             scat_i_lon_,
             f_grid_,
             scat_f_index_)));

 agenda_data.push_back
    (AgRecord
     ( NAME( "scat_rte_agenda" ),
       DESCRIPTION
       (
       "Calculate the radiative transfer equation (RTE) with fixed scattering\n "
       "integral.\n"
       "\n"
       "There are different possibilities to compute the RTE with a fixed value \n"
       "for the scattering integral. These are implemented in the workspace \n"
       "methods: \n"
       "sto_vecGeneral: This method uses the Pade approximation to compute the \n"
       "         matrix exponential function and a LU decomposition method \n"
       "         to compute the matrix inverse. This method should not be \n"
       "         used for the scalar RTE (if *stokes_dim* equals 1), as it\n"
       "         is numerically not efficient in this case.\n" 
       "sto_vecScalar: This method can only be used for the scalar RTE. It \n"
       "         uses the standard exponential function and there is no need\n"
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
       "   a_planck_value: Planck function. \n"
       "   stokes_dim  : Stokes dimension. \n"
       ""
	),
       OUTPUT(stokes_vec_),
       INPUT(ext_mat_,
             abs_vec_,
             sca_vec_,
             l_step_,
             a_planck_value_,
             stokes_dim_)));

 agenda_data.push_back
    (AgRecord
     ( NAME( "spt_calc_agenda" ),
       DESCRIPTION
       (
	"Calculates single particle properties from the amplitude matrix.\n"
        "\n"
        "This agenda sets up the methods, which should be used to calculate \n"
        "the particle properties, i.e. the extinction matrix, the absorbtion\n"
        "vector and the phase matrix from the amplitude matrix for each \n"
        "particle type specified in the control file. \n"
        "\n"
        "Normally you have tu use:\n"
        "1. *pha_mat_sptCalc* \n"
        "2. *ext_mat_sptCalc* \n"
        "3. *abs_vec_sptCalc* \n"
        "Note: the order of calling these methods is important.\n"
        "\n"
        "It can be useful to compute the extinction matrix without \n"
        "particle absorption, for examle to do a convergence test. \n"
        "Then the method *ext_mat_sptScat* has to be used. \n"
        "\n"
	),
       OUTPUT( ext_mat_spt_, abs_vec_spt_, pha_mat_spt_),
       INPUT(  pha_mat_spt_, abs_vec_spt_, ext_mat_spt_, amp_mat_,
               scat_za_index_, scat_aa_index_,
               scat_za_grid_, scat_aa_grid_ )));

}
