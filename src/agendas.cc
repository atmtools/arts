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
               i_field_old_)));


 
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
        "described by the variables *ground_los* and *ground_refl_coeffs*.\n"
        "\n"
        "A function calling this agenda shall set *rte_los* and *rte_gp_XXX*\n"
        "to match the line-of-sight and position at the ground reflection\n"
        "point.\n"
        "\n"
        "See further the user guide.\n"
        "\n"
        "Usage:   Called from *RteCalc*."
        ),
       OUTPUT( ground_emission_, ground_los_, ground_refl_coeffs_  ),
       INPUT(  f_grid_, stokes_dim_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_, 
               rte_los_, r_geoid_, z_ground_, t_field_ )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "i_space_agenda" ),
       DESCRIPTION
       (
        "Sets the workspace variable *i_space* to match the assumptions \n"
        "regarding the radiation entering the atmosphere at the start of the\n"
        "propagation path. \n"
        "\n"
        "The main usage of this agenda is to be called from *RteCalc*. \n"
        "\n"
        "A function calling this agenda shall set *rte_pos* and *rte_los* to\n"
        "the position and line-of-sight for which the entering radiation \n"
        "shall be determined. The position and line-of-sight must be known, \n"
        "for example, when radiation from the sun is considered. \n"
        "\n"
        "Usage:   Called from *RteCalc*."
        ),
       OUTPUT(  i_space_ ),
       INPUT(  f_grid_, stokes_dim_, rte_pos_, rte_los_ )));

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
     ( NAME( "met_profile_calc_agenda" ),
       DESCRIPTION
       (
        "This agenda is used for metoffice profile calculations.\n"
	"\n"
	"This agenda is called inside the method *ybatchMetProfiles* which is\n"
	"used to make a batch calculation for the metoffice profiles.   \n"
	"See the documentation of *ybatchMetProfiles* for more information.\n"
	"\n"
	"This agenda can be, for example, set up like this:\n"
	"\n"
	"*AtmFieldsCalc*\n"
	"*gas_abs_lookupAdapt*\n"
	"*ScatteringInit	*\n"
	"*CloudboxGetIncoming*\n"
	"*ScatteringMain*\n"
	"*RteCalc*\n"
	"*yNoPolarisation*\n"
	"\n"
	"For example, if you want the output in brightness temperature unit,\n"
	"then add the method *VectorToTbByPlanck*.  \n"
	"\n"
        ),
       OUTPUT( y_ ),
       INPUT(t_field_raw_, vmr_field_raw_, z_field_raw_)));

  agenda_data.push_back
    (AgRecord
     ( NAME( "opt_prop_gas_agenda" ),
       DESCRIPTION
       (
        "Calculate the optical properties (absorption vector and extinction.\n"
        "matrix) of gaseous species at a given grid point.\n"
        "\n"
        "This agenda, for example, can be defined in the following manner:\n"
        "\n"
        "*ext_matAddGas* : This method calculates the extinction \n"
        "                  matrix for the gaseous species and adds it to \n"
        "                  the workspace variable *ext_mat*.\n"
        "*abs_vecAddGas* : This method calculates the absorption\n"
        "                  vector for the gaseous species and adds it to\n"
        "                  the workspace variables abs_vec.\n"     
        "If the Zeeman effect should be included the following methods have \n"
        "to be added: \n"
        "*ext_matAddZeeman* \n"
        "*abs_vecAddZeeman* \n"
        " \n"
        "Note that the initialization of *abs_vec* is not done inside the\n"
        "agenda, so *abs_vec* has to be initialize before executing the \n"
        "agenda.\n"
        "\n"
        "Output :\n"
        "   ext_mat     : Extinction matrix.\n"
        "   abs_vec     : Absorption vector.\n"
        "\n"
        "Input:\n"
        "   ext_mat     : Extinction matrix.\n"
        "   abs_vec     : Absorption vector. \n"
        "   abs_scalar_gas : Scalar gas absorption. \n"
        ),
       OUTPUT( ext_mat_, abs_vec_ ),
       INPUT( ext_mat_, abs_vec_, 
              abs_scalar_gas_)));


 agenda_data.push_back
    (AgRecord
     ( NAME( "opt_prop_part_agenda" ),
       DESCRIPTION
       (
        "Calculate the optical properties (absorption vector and extinction.\n"
        "matrix) for particles at a given atmospheric grid point.\n"
        "\n"
        "This agenda, for example, can be defined in the following manner:\n"
        "\n"
        "*ext_matAddPart* : This method calculates the extinction \n"
        "                   matrix for particles and adds it to the \n"
        "                   workspace variable *ext_mat*.\n"
        "*abs_vecAddPart* : This method calculates the absorption\n"
        "                   vector for particles and adds it to the\n"
        "                   workspace variables abs_vec.\n"     
        " \n"
        "Note that the initialization of *ext_mat* is not done inside the\n"
        "agenda, so *ext_mat* has to be initialize before executing the \n"
        "agenda.\n"
        "\n"
        "Output :\n"
        "   ext_mat     : Extinction matrix.\n"
        "   abs_vec     : Absorption vector. \n"
        "\n"
        "Input:\n"
        "   ext_mat     : Extinction matrix. \n"
        "   ext_mat_spt : Extinction matrix for single particle type. \n"
        "   abs_vec     : Absorption vector. \n"
        "   abs_vec_spt : Absorption vector for single particle type. \n"
        "   pnd_field   : Particle number density field. \n"
        "   atmosphere_dim: Atmospheric dimension. \n"
        "   scat_p_index : Position. \n:"
        "   scat_lat_index : Position. \n"
        "   scat_lon_index : Position. \n"
        ),
       OUTPUT( ext_mat_, abs_vec_ ),
       INPUT( ext_mat_, abs_vec_, 
              ext_mat_spt_, abs_vec_spt_,
              pnd_field_, atmosphere_dim_, scat_p_index_, scat_lat_index_,
              scat_lon_index_)));

 agenda_data.push_back
    (AgRecord
     ( NAME( "pha_mat_spt_agenda" ),
       DESCRIPTION
       (
        "Calculates the phase matrix for a single particle type.\n"
        "\n"
        "Different options are possible for the usage of this agenda: \n"
        "*pha_mat_sptFromData* or *pha_mat_sptDOITOpt*. \n"
        "\n"
        ),
       OUTPUT( pha_mat_spt_),
       INPUT( pha_mat_spt_, scat_theta_, scat_za_grid_, scat_aa_grid_)));
       

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
     ( NAME( "refr_index_agenda" ),
       DESCRIPTION
       (
        "Calculates the refractive index.\n"
        "\n"
        "This agenda should calculate the summed refractive index for all\n"
	"relevant constituients. The result is returned in *refr_index*, the\n"
        "atmospheric state is speciefied by *rte_pressure*, \n"
        "*rte_temperature* and *rte_vmr_list*."
        ),
       OUTPUT( refr_index_ ),
       INPUT(  rte_pressure_, rte_temperature_, rte_vmr_list_ )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "rte_agenda" ),
       DESCRIPTION
       (
        "Performs monochromatic pencil beam calculations for a single\n"
        "propagation path.\n"
        "\n"
        "When calling the agenda, *i_rte* shall be set to the radiances, or\n"
        "optical thicknesses, at the start of the propagation path described\n"
        "by *ppath*. The agenda then solves the radiative transfer equation\n"
        "along the propagation path and returns the result in *i_rte*."
        ),
       OUTPUT( i_rte_ ),
       INPUT(  ppath_, atmosphere_dim_, stokes_dim_, f_grid_, p_grid_,
               lat_grid_, lon_grid_, t_field_, vmr_field_,
               scalar_gas_absorption_agenda_, opt_prop_gas_agenda_ )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "scalar_gas_absorption_agenda" ),
       DESCRIPTION
       (
        "Calculate scalar gas absorption.\n"
        "\n"
        "This agenda should calculate absorption coefficients for all gas\n"
        "species as a function of the given atmospheric state for one point in\n"
        "the atmosphere. The result is returned in *abs_scalar_gas*, the\n"
        "atmospheric state has to be specified by *rte_pressure*,\n"
        "*rte_temperature*, and *rte_vmr_list*\n"
        "\n"
        "A mandatory input parameter is f_index, which is used as follows:\n"
        "\n"
        "1. f_index < 0 : Return absorption for all frequencies (in f_grid).\n"
        "\n"
        "2. f_index >= 0 : Return absorption for the frequency indicated by\n"
        "   f_index. \n"
        "\n"
        "The methods inside this agenda may require a lot of additional input\n"
        "variables, such as *f_grid*, *species*, etc.."
        ),
       OUTPUT( abs_scalar_gas_ ),
       INPUT(  f_index_,
               rte_pressure_, rte_temperature_, rte_vmr_list_ )));
  
  agenda_data.push_back
    (AgRecord
     ( NAME( "scat_field_agenda" ),
       DESCRIPTION
       (
        "Calculation of the scattering integral. \n"
        "\n"
        "Different methods for calculating the scattering integral have been\n"
        "implemented. Here you can define the method which should be taken \n"
        "in your calculation. \n"
        "\n"
        ),
        OUTPUT( scat_field_ ),
        INPUT(  i_field_, pnd_field_, scat_za_grid_, scat_aa_grid_, p_grid_, 
                lat_grid_, lon_grid_, cloudbox_limits_)));


  agenda_data.push_back
    (AgRecord
     ( NAME( "scat_mono_agenda" ),
       DESCRIPTION
       (
       "Performs the monochromatic scattering calculation."
       "\n"
       "Normally this agenda consists of four methods: \n"
       "   1. amp_matCalc\n"
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
       "   f_index: Frequency index for the ongoing scattering \n"
       "                 calculation. \n"
       "   scat_za_grid: Zenith angle grid. \n"
       "   scat_aa_grid: Azimuthal angle grid. \n"
       "   amp_mat_raw: Amplitude matrix raw data. \n"
       ""
        ),
       OUTPUT(scat_i_p_,
              scat_i_lat_,
              scat_i_lon_),
       INPUT(scat_i_p_,
             scat_i_lat_,
             scat_i_lon_,
             f_grid_,
             f_index_,
             scat_za_grid_,
             scat_aa_grid_)));

 agenda_data.push_back
    (AgRecord
     ( NAME( "scat_rte_agenda" ),
       DESCRIPTION
       (
        "Radiative transfer calculation in cloudbox. \n"
        "\n"
        "Differnt methods for the radiative transfer calculation are:\n"
        "*i_fieldUpdate{1,3}D*: 'Normal' function for successive order of\n"
        "                      scattering method. \n"
        "*i_fieldUpdateSeq{1,3}D*: Seqential update of the radiation field. \n"
        "                      This method is much faster than the 'normal'\n"
        "                      method. \n"
        "*i_fieldUpdate{1,3}DPlaneParallel*: The same method as \n"
        "                      *i_fieldUpdateSeq*** but with plane- \n"
        "                      oarallel geometry. \n"
        "\n"
        ),
        OUTPUT( i_field_),
        INPUT(  scalar_gas_absorption_agenda_, spt_calc_agenda_,
                opt_prop_part_agenda_, opt_prop_gas_agenda_,
                ppath_step_agenda_)));

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
       OUTPUT( ext_mat_spt_, abs_vec_spt_),
       INPUT(  abs_vec_spt_, ext_mat_spt_,
               scat_za_index_, scat_aa_index_,
               scat_za_grid_, scat_aa_grid_ )));

}
