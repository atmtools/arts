/* Copyright (C) 2002-2008 Stefan Buehler <sbuehler@ltu.se>

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
  \author Stefan Buehler <sbuehler@ltu.se>
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

//   agenda_data.push_back
//     (AgRecord
//      ( NAME( "abs_coef_agenda" ),
//        DESCRIPTION
//        (
//         "FIXME: This is just a placeholder for now."
//         ),
//        OUTPUT( abs_coef_ ),
//        INPUT(  abs_species_,
//                abs_lines_per_species_,
//                abs_lineshape_,
//                abs_vmrs_,
//                abs_t_ )));
  
  agenda_data.push_back
    (AgRecord
     ( NAME( "abs_scalar_gas_agenda" ),
       DESCRIPTION
       (
        "Calculate scalar gas absorption.\n"
        "\n"
        "This agenda should calculate absorption coefficients for all gas\n"
        "species as a function of the given atmospheric state for one point\n"
        "in the atmosphere. The result is returned in *abs_scalar_gas*, the\n"
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
        "The methods inside this agenda may require a lot of additional\n"
        "input variables, such as *f_grid*, *species*, etc.."
        ),
       OUTPUT( abs_scalar_gas_ ),
       INPUT(  f_index_,
               rte_pressure_, rte_temperature_, rte_vmr_list_ )));
  
 agenda_data.push_back
    (AgRecord
     ( NAME( "ybatch_calc_agenda" ),
       DESCRIPTION
       (
        "Calculations to perform for each batch case.\n"
        "\n"
        "Must produce a new spectrum vector (y). See further *ybatchCalc*."
        ),
       OUTPUT( y_ ),
       INPUT( ybatch_index_ )));

 agenda_data.push_back
    (AgRecord
     ( NAME( "doit_conv_test_agenda" ),
       DESCRIPTION
       (
        "Compute the convergence test.\n"
        "\n"
        "The method *scat_i_fieldIterate* solves the VRTE iteratively."
        "This method requires \n"
        "a convergence test. The user can choose different convergence tests\n"
        "which are to be defined in this agenda.\n"
        "\n"
        "Possible workspace methods are:\n"
        "*doit_conv_flagAbs*: Calculates the absolute differences \n"
        "  for each Stokes component separately.\n"
        "*doit_conv_flagAbs_BT*: Same as above, but the convergence limit\n"
        "  can be specified in Kelvin BT (Rayleigh Jeans).\n"
        "*doit_conv_flagLsq*: Least square convergence test. Not recommended\n"
        "  because result can be inaccurate.\n"
        "\n"
        ),
       OUTPUT( doit_conv_flag_, doit_iteration_counter_ ),
       INPUT(  doit_conv_flag_, doit_iteration_counter_,
               doit_i_field_, doit_i_field_old_)));

 agenda_data.push_back
    (AgRecord
     ( NAME( "doit_scat_field_agenda" ),
       DESCRIPTION
       (
        "Calculation of the scattering integral field (DOIT). \n"
        "\n"
        "This agenda is called repeatedly in each DOIT iteration.\n"
        "The following methods can be used for calculating the \n"
        "scattering integral field: \n"
        "\n"
        "*doit_scat_fieldCalc*: This method calculates the scattering \n"
        "  integral field by using the angular grids *scat_za_grid* \n"
        "  and *scat_aa_grid*, which are also used in the update of the \n"
        "  radiation field (*doit_rte_agenda*).\n"
        "\n"
        "*doit_scat_fieldCalcLimb*: This method calculates the scattering \n"
        "  integral field.  The difference to the previous method is that \n"
        "  the data is interpolated on equidistant angular grids. \n"
        "  Especially for limb, where a very fine zenith angle grid \n"
        "  resolution is required for the RT transfer part, this method \n"
        "  is much faster than *doit_scat_fieldCalc*. \n"
        "\n"
        ),
        OUTPUT( doit_scat_field_ ),
        INPUT(  doit_scat_field_, doit_i_field_)));

  agenda_data.push_back
    (AgRecord
     ( NAME( "doit_mono_agenda" ),
       DESCRIPTION
       (
       "Performs monochromatic DOIT calculation."
       "\n"
       "This agenda includes for example the following methods:\n"
       "   1. *DoitScatteringDataPrepare* \n"
       "   2. *doit_i_fieldSetClearsky \n"
       "   3. *doit_i_fieldIterate*\n"
       "   4. *DoitCloudboxFieldPut*\n"
       "\n"
       "The result of the agenda is the radiation field inside the \n"
       "cloudbox and on the cloudbox boundary, which can be used \n"
       "as radiative background for a clearsky radiative transfer \n"
       "calculation. \n"
       "\n"
       "See the ArtsWiki page *UsingArtsDoit* and the online documentation\n"
       "for more information about\n"
       "the methods.\n"
       "\n"
        ),
       OUTPUT(doit_i_field_,scat_i_p_,scat_i_lat_, scat_i_lon_, 
              doit_i_field1D_spectrum_),
       INPUT(f_index_,scat_i_p_,scat_i_lat_, scat_i_lon_)));
            
 agenda_data.push_back
    (AgRecord
     ( NAME( "doit_rte_agenda" ),
       DESCRIPTION
       (
        "Radiative transfer calculations in cloudbox.\n"
        "\n"
        "Agenda for radiative transfer step calculations with \n"
        "fixed scattering integral term shoul be specified here.\n"
        "Output is the updated radiation field in the cloudbox. \n"
        "This agenda is called repeatedly in each DOIT iteration.\n"
        "\n"
        "Normally one should use \n"
        "*doit_i_fieldUpdateSeq{1,3}D*: Seqential update of the"
        "radiation field.\n"
        "   This method is the fastest and most accurate method.\n"
        "\n"
        "A very similar method in plane parallel approximation is \n"
        "*doit_i_fieldUpdate{1,3}DPlaneParallel*: This method also \n"
        "   incluldes the sequential update, and it is slightly faster than\n"
        "   *doit_i_fieldUpdateSeq{1,3}D*. The drawback is, that it is less\n"
        "   accurate, especially for limb geometries and for larger \n"
        "   off-nadir viewing angles. \n"
        "\n"
        "The following methods were used before the sequential update\n"
        "was invented. They are very slow and should therefore only \n"
        "be used for test cases.\n"
        "*doit_i_fieldUpdate{1,3}D*: Old function.\n"
        "\n"
        ),
        OUTPUT( doit_i_field_),
       INPUT( doit_i_field_, doit_scat_field_)));
 
  agenda_data.push_back
    (AgRecord
     ( NAME( "emission_agenda" ),
       DESCRIPTION
       (
        "Calculates the thermal emission source term.\n"
        "\n"
        "This agenda shall return the emission at one position along\n"
        "the propagation path. This source term is used for calculations\n"
        "inside *rte_agenda*. Inside scattering methods, such as DOIT,\n"
        "the calculation of the source term can be hard coded.\n"
        "\n"
        "By setting *emission* to zeros and *iy_space* to ones, the\n"
        "obtained spectrum will be the transmission through the atmosphere."
        ),
       OUTPUT( emission_ ),
       INPUT( rte_temperature_ )));
 
  agenda_data.push_back
    (AgRecord
     ( NAME( "geomag_los_calc_agenda" ),
       DESCRIPTION
       (
        "Calculates the magnetic field along a given propagation path.\n"
        "\n"
        "The agenda relates the vector of the geomagnetic field to \n"
        "a specified propagation path. As a result the magnitude of \n"
        "this vector is calculated in each point of the propagation \n"
        "path, alongside with the corresponding angle between the \n"
        "geomagnetic field vector and the propagation direction.  \n"
        "The output is the WSV *geomag_los*, containing the two \n"
        "quantities discussed above. \n"
        "\n"
        "Output:    \n"
        "   geomag_los : Magnetic field along LOS plus angle  \n"
        "\n"
        "Input: ppath_   \n"
        "       geomag_intensitities.xml \n"
        "\n"
          ),
       OUTPUT( geomag_los_ ),
       INPUT(  )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "iy_cloudbox_agenda" ),
       DESCRIPTION
       (
        "Sets *iy* to scattered radiation for given position and LOS.\n"
        "\n"
        "The calculations inside the agenda differ depending on scattering\n"
        "solution method. If DOIT is used, an interpolate of the intensity\n"
        "field shall be performed. Another option is to start backward Monte\n"
        "Carlos calculations from this point.\n"
        "\n"
        "A function calling this agenda shall set *rte_pos* and *rte_los* to\n"
        "the position and line-of-sight for which the scattered radiation \n"
        "shall be determined. \n"
        "\n"
        "Usage:   Called from *RteCalc*."
        ),
       OUTPUT( iy_, ppath_, rte_pos_, rte_los_ ),
       INPUT( ppath_, rte_pos_, rte_los_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_ )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "iy_space_agenda" ),
       DESCRIPTION
       (
        "Sets the workspace variable *iy* to match the assumptions \n"
        "regarding the radiation entering the atmosphere at the start of the\n"
        "propagation path. \n"
        "\n"
        "A function calling this agenda shall set *rte_pos* and *rte_los* to\n"
        "the position and line-of-sight for which the entering radiation \n"
        "shall be determined. The position and line-of-sight must be known, \n"
        "for example, when radiation from the sun is considered. \n"
        "\n"
        "Usage:   Called from *RteCalc*."
        ),
       OUTPUT( iy_ ),
       INPUT( rte_pos_, rte_los_ )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "jacobian_agenda" ),
       DESCRIPTION
       (
        "The agenda controlling the calculation of the Jacobian matrix.\n"
        "This agenda is not supposed to be set by the user, it should be\n"
        "automatically be constructed when defining the jacobian quantities.\n"
        "\n"
        "Usage:   Called from *jacobianCalc*."
       ),
       OUTPUT( jacobian_ ),
       INPUT())),
       
  agenda_data.push_back
    (AgRecord
     ( NAME( "jacobian_particle_update_agenda" ),
       DESCRIPTION
       (
        "The agenda controlling the update of the scattered field due to\n"
        "changes in *pnd_field* when calculating the particle Jacobian.\n"
        "The agenda has to be specified by the user, and should contain\n"
        "the calculations needed before a call to RteCalc.\n"
        "\n"
        "Usage:   Called from *jacobianCalcParticle*."
       ),
       OUTPUT( ),
       INPUT( pnd_field_ ))),
       
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
        "*abs_lookupAdapt*\n"
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
       INPUT(t_field_raw_, vmr_field_raw_, z_field_raw_, pnd_field_raw_,
             p_grid_, sensor_los_, cloudbox_on_, cloudbox_limits_,
             z_surface_)));

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
       INPUT( ext_mat_, abs_vec_, abs_scalar_gas_)));


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
              scat_p_index_, scat_lat_index_,
              scat_lon_index_ )));

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
       INPUT( pha_mat_spt_, scat_za_index_, scat_lat_index_, scat_lon_index_,
              scat_p_index_, scat_aa_index_, rte_temperature_)));
       

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
        "exception. If there is an intersection of a blackbody surface, the \n"
        "calculations stop at this point. This is flagged by setting the \n" 
        "background field of *ppath_step*. Beside this, the calling method \n" 
        "must check if the starting point of the calculations is inside the \n"
        "scattering box or below the surface level, and check if the last \n"
        "point of the path has been reached. The starting point (the end \n"
        "furthest away from the sensor) of a full propagation path can be \n"
        "the top of the atmosphere, the surface and the cloud box.\n"
        "\n"
        "The *ppath_step_agenda* put in points along the propagation path \n"
        "at all crossings with the grids, tangent points and points of \n"
        "surface reflection. It is also allowed to make agendas that put in \n"
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
       INPUT(ppath_step_, atmosphere_dim_, p_grid_, lat_grid_,
             lon_grid_, z_field_, r_geoid_, z_surface_)));
  
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
        "When calling the agenda, *iy* shall be set to the radiances, or\n"
        "optical thicknesses, at the start of the propagation path described\n"
        "by *ppath*. The agenda then solves the radiative transfer equation\n"
        "along the propagation path and returns the result in *iy*."
        ),
       OUTPUT( iy_, diy_dvmr_, diy_dt_ ),
       INPUT( iy_, diy_dvmr_, diy_dt_, ppath_, ppath_array_, ppath_array_index_,
              rte_do_vmr_jacs_, rte_do_t_jacs_, stokes_dim_, f_grid_ )));

 agenda_data.push_back
    (AgRecord
     ( NAME( "spt_calc_agenda" ),
       DESCRIPTION
       (
        "Calculates single particle properties from the amplitude matrix.\n"
        "\n"
        "This agenda sets up the methods, which should be used to calculate \n"
        "the particle properties, i.e. the extinction matrix and the \n"
        "absorbtion vector.\n "
        "\n"
        "Normally you  use:\n"
        " opt_prop_sptFromMonoData{} \n"
        "\n"
        ),
       OUTPUT( ext_mat_spt_, abs_vec_spt_),
       INPUT(  ext_mat_spt_, abs_vec_spt_,
               scat_p_index_, scat_lat_index_, scat_lon_index_,
               rte_temperature_, scat_za_index_, scat_aa_index_
               )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "surface_prop_agenda" ),
       DESCRIPTION
       (
        "Returns RT properties of the surface. \n"
        "\n"
        "See the user guide for closer definitions of the variables \n"
        "that describe the surface properties. These variables are:\n"
        "   *surface_emission*, *surface_los* and *surface_rmatrix* \n"
        "\n"
        //        "A function calling this agenda shall set *rte_gp_p/lat/lon* to\n"
        //"the position of intersection with the surface."
        ),
       OUTPUT( surface_emission_, surface_los_, surface_rmatrix_ ),
       INPUT( rte_gp_p_, rte_gp_lat_, rte_gp_lon_, rte_los_ )));


//  agenda_data.push_back
//     (AgRecord
//      ( NAME( "zeeman_prop_agenda" ),
//        DESCRIPTION
//        (
//         "Calculates extinction matrix and absorption vector due to the \n"
//         "Zeeman effect. \n"
//         "\n"
//         "The agenda calculates the total extinction matrix and absorption \n"
// 	"vector of O2 only (currently) due to the Zeeman effect induced by the \n"
// 	"geomagnetic field. The polarization pattern, resulting from the effect, \n"
// 	"produces additional contribution to these quantities in the unpolarized \n"
// 	"case. This means that the user should take care not to execute both the \n"
// 	"polarized (Zeeman) and unpolarized  calculation for the O2 lines \n"
// 	"affected by the Zeeman effect, otherwise the unpolarized part gets \n"
// 	"wrongly doubled. \n"
//         "\n"
//         "Output:    \n"
//         "   ext_mat_zeeman: Zeeman extinction matrix  \n"
//         "   abs_vec_zeeman: Zeeman absorption vector  \n"
//         "\n"
//         "Input:    \n"
//         "  geomag_los: Magnetic field along LOS plus angle  \n"
//         "\n"
//           ),
//        OUTPUT(  ),
//        INPUT(  )));


}
