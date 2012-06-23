/* Copyright (C) 2002-2012 Stefan Buehler <sbuehler@ltu.se>

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


// Some #defines and typedefs to make the records better readable:
#define NAME(x) x 
#define DESCRIPTION(x) x
#define OUTPUT   MakeArray<String>
#define INPUT    MakeArray<String>

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
     ( NAME( "abs_scalar_gas_agenda" ),
       DESCRIPTION
       (
        "Calculation of scalar gas absorption.\n"
        "\n"
        "This agenda calculates absorption coefficients for all gas species\n"
        "as a function of the given atmospheric state for one point in the\n"
        "atmosphere. The result is returned in *abs_scalar_gas*, the\n"
        "atmospheric state has to be specified by *rte_pressure*,\n"
        "*rte_temperature*, and *rte_vmr_list*.\n"
        "\n"
        "A mandatory input parameter is f_index, which is used as follows:\n"
        "\n"
        "1. f_index < 0 : Return absorption for all frequencies (in f_grid).\n"
        "\n"
        "2. f_index >= 0 : Return absorption for the frequency indicated by\n"
        "   f_index. \n"
        "\n"
        "The methods inside this agenda may require a lot of additional\n"
        "input variables, such as *f_grid*, *abs_species*, etc.\n"
        ),
       OUTPUT( "abs_scalar_gas" ),
       INPUT(  "f_index", "rte_doppler", "rte_pressure", "rte_temperature", 
               "rte_vmr_list" )));
  
 agenda_data.push_back
    (AgRecord
     ( NAME( "doit_conv_test_agenda" ),
       DESCRIPTION
       (
        "Compute the convergence test.\n"
        "\n"
        "The method *doit_i_fieldIterate* solves the VRTE iteratively."
        "This method requires \n"
        "a convergence test. The user can choose different convergence tests\n"
        "which are to be defined in this agenda.\n"
        "\n"
        "Possible workspace methods are:\n"
        "*doit_conv_flagAbs*: Calculates the absolute differences \n"
        "  for each Stokes component separately.\n"
        "*doit_conv_flagAbsBT*: Same as above, but the convergence limit\n"
        "  can be specified in Kelvin BT (Rayleigh Jeans).\n"
        "*doit_conv_flagLsq*: Least square convergence test. Not recommended\n"
        "  because result can be inaccurate.\n"
        ),
       OUTPUT( "doit_conv_flag", "doit_iteration_counter" ),
       INPUT(  "doit_conv_flag", "doit_iteration_counter",
               "doit_i_field", "doit_i_field_old" )));

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
        ),
        OUTPUT( "doit_scat_field" ),
        INPUT(  "doit_scat_field", "doit_i_field")));

  agenda_data.push_back
    (AgRecord
     ( NAME( "doit_mono_agenda" ),
       DESCRIPTION
       (
       "Performs monochromatic DOIT calculation."
       "\n"
       "This agenda includes for example the following methods:\n"
       "   1. *DoitScatteringDataPrepare* \n"
       "   2. *doit_i_fieldSetClearsky* \n"
       "   3. *doit_i_fieldIterate*\n"
       "   4. *DoitCloudboxFieldPut*\n"
       "\n"
       "The result of the agenda is the radiation field inside the \n"
       "cloudbox and on the cloudbox boundary, which can be used \n"
       "as radiative background for a clearsky radiative transfer \n"
       "calculation. \n"
       "\n"
       "See the ArtsWiki page UsingArtsDoit and the online documentation\n"
       "for more information about the methods.\n"
        ),
       OUTPUT( "doit_i_field", "scat_i_p", "scat_i_lat", "scat_i_lon", 
               "doit_i_field1D_spectrum"),
       INPUT("f_index", "scat_i_p", "scat_i_lat", "scat_i_lon")));
            
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
        ),
        OUTPUT( "doit_i_field" ),
       INPUT(   "doit_i_field", "doit_scat_field" )));
 
  agenda_data.push_back
    (AgRecord
     ( NAME( "emission_agenda" ),
       DESCRIPTION
       (
        "Thermal emission source term.\n"
        "\n"
        "This agenda shall return the emission at one position along\n"
        "the propagation path. The source term equals the Planck function as\n"
        "long as thermodynamic equilibrium (LTE) can be assumed, while for\n"
        "non-LTE conditions much more complex calculations are required.\n"
        ),
       OUTPUT( "emission" ),
       INPUT( "rte_temperature", "f_grid" )));
 
  agenda_data.push_back
    (AgRecord
     ( NAME( "forloop_agenda" ),
       DESCRIPTION
       (
        "The body for a for loop.\n"
        "\n"
        "This agenda contains the body of the for loop to be execute by the\n"
        "method *ForLoop*.\n"
        ),
       OUTPUT(),
       INPUT( "forloop_index" )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "fos_y_agenda" ),
       DESCRIPTION
       (
        "Calculation of incoming radiation field for FOS method.\n"
        "\n"
        "The direct task of the agenda is to determine the incoming radiation\n"
        "field, to evaluate of the scattering integral, for each angle in\n"
        "*fos_angle*. The data are packed into *fos_y*.\n"
        "\n"
        "The underlying purpose of this agenda is to allow different\n"
        "shortcuts for estimating the incoming radiation field. For example,\n"
        "calculations can be performed for a few directions and then an\n"
        "interpolation is performed to obtain the intensity for other\n"
        "directions. The data could also be taken from a pre-calculated\n"
        "database.\n"
        ),
       OUTPUT( "fos_y" ),
       INPUT( "rte_pos", "fos_angles", "fos_n", "fos_i" )));
 
//   agenda_data.push_back
//     (AgRecord
//      ( NAME( "geomag_los_calc_agenda" ),
//        DESCRIPTION
//        (
//         "Calculates the magnetic field along a given propagation path.\n"
//         "\n"
//         "The agenda relates the vector of the geomagnetic field to \n"
//         "a specified propagation path. As a result the magnitude of \n"
//         "this vector is calculated in each point of the propagation \n"
//         "path, alongside with the corresponding angle between the \n"
//         "geomagnetic field vector and the propagation direction.  \n"
//         "The output is the WSV *geomag_los*, containing the two \n"
//         "quantities discussed above. \n"
//         "\n"
//         "Output:    \n"
//         "   geomag_los : Magnetic field along LOS plus angle  \n"
//         "\n"
//         "Input: ppath_   \n"
//         "       geomag_intensitities.xml \n"
//         "\n"
//           ),
//        OUTPUT( "geomag_los" ),
//        INPUT(  )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "g0_agenda" ),
       DESCRIPTION
       (
        "Calculation of the gravity at zero altitude.\n"
        "\n"
        "Returns *g0* for given geographical position.\n"
        ),
       OUTPUT( "g0" ),
       INPUT( "lat", "lon" )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "iy_clearsky_agenda" ),
       DESCRIPTION
       (
        "Calculation of a single monochromatic pencil beam spectrum.\n"
        "\n"
        "The task of the agenda is to calculate the monochromatic pencil beam\n"
        "spectrum for the position specified by *rte_pos* and the viewing\n"
        "direction specified by *rte_los*. This includes cases when the\n"
        "propagation path intersects with the surface or the cloudbox.\n"
        ),
       OUTPUT( "iy", "iy_error", "iy_error_type", "iy_aux", "diy_dx" ),
       INPUT( "iy_error", "iy_error_type", "iy_aux", "diy_dx",
              "iy_agenda_call1", "iy_transmission", "rte_pos", "rte_los", 
              "cloudbox_on", "jacobian_do", "t_field", "z_field", "vmr_field", 
              "mblock_index" )));  

  agenda_data.push_back
    (AgRecord
     ( NAME( "iy_clearsky_basic_agenda" ),
       DESCRIPTION
       (
        "As *iy_clearsky_agenda*, but lacks all support for jacobian and\n"
        "auxiliary variables.\n"
        "\n"
        "This agenda is used by scattering methods without support for the\n"
        "jacobian and auxilary parts, in order to calculate unscattered\n"
        "radiation for the speciefied position and line-of-sight. As this\n"
        "agenda has much fewer input variables, this simplifies the interface\n"
        "and maintaince of those scattering methods.\n"
        ),
       OUTPUT( "iy" ),
       INPUT( "rte_pos", "rte_los", "cloudbox_on" )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "iy_cloudbox_agenda" ),
       DESCRIPTION
       (
        "Intensity at cloudbox boundary or interior.\n"
        "\n"
        "The task of the agenda is to determine the intensity at some point\n"
        "at the boundary of inside the cloudbox.  The actual calculations\n"
        "inside the agenda differ depending on scattering solution method.\n"
        "If DOIT is used, an interpolating of the intensity field should be\n"
        "performed. Another option is to start backward Monte Carlo \n"
        "calculations from this point.\n"
        "\n"
        "A function calling this agenda shall set *rte_pos* and *rte_los* to\n"
        "the position and line-of-sight for which the scattered radiation\n"
        "shall be determined.\n"
        ),
       OUTPUT( "iy", "iy_error", "iy_error_type", "iy_aux", "diy_dx" ),
       INPUT( "iy_error", "iy_error_type", "iy_aux", "diy_dx", 
              "iy_transmission", "rte_pos", "rte_los" )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "iy_space_agenda" ),
       DESCRIPTION
       (
        "Downwelling radiation at the top of the atmosphere.\n"
        "\n"
        "Possible terms to include in this agenda include cosmic background\n"
        "radiation and solar radiation.\n"
        "\n"
        "A function calling this agenda shall set *rte_pos* and *rte_los* to\n"
        "the position and line-of-sight for which the entering radiation \n"
        "shall be determined. The position and line-of-sight must be known, \n"
        "for example, when radiation from the sun is considered.\n"
        ),
       OUTPUT( "iy" ),
       INPUT( "rte_pos", "rte_los", "f_grid" )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "iy_surface_agenda" ),
       DESCRIPTION
       (
        "Upwelling radiation from the surface.\n"
        "\n"
        "The task of the agenda is to determine the upwelling intensity from\n"
        "the surface, for given point and direction.\n"
        "\n"
        "The standard choice should be to make use of *surface_rtprop_agenda*\n"
        "throught the WSM *iyFromSurfaceRtpropAgenda*.\n"
        "\n"
        "A function calling this agenda shall set *rte_pos* and *rte_los* to\n"
        "the position and line-of-sight for which the upwelling radiation\n"
        "shall be determined.\n"
        ),
       OUTPUT( "iy", "iy_error", "iy_error_type", "iy_aux", "diy_dx" ),
       INPUT( "iy_error", "iy_error_type", "iy_aux", "diy_dx", 
              "iy_transmission", "rte_pos", "rte_los" )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "jacobian_agenda" ),
       DESCRIPTION
       (
        "Pure numerical Jacobian calculations.\n"
        "\n"
        "Parts of the Jacobian matrix can be determined by (semi-)analytical\n"
        "expressions, while other parts are calculated in apure numerical\n"
        "manner (by perturbations). This agenda describes the calculations to\n"
        "be performed in the later case.\n"
        "\n"
        "This agenda is normally not set directly by the user, but is created\n"
        "by calling the the jacobianAdd set of methods.\n"
       ),
       OUTPUT( "jacobian" ),
       INPUT( "jacobian", "mblock_index", "iyb", "yb" )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "jacobian_y_agenda" ),
       DESCRIPTION
       (
        "Agenda providing recalculated *y* after some perturbation.\n"
        "\n"
        "The purpose of this agenda is to determine some jacobians through\n"
        "perturbation calculations. Accordingly, the agenda shall return *y*\n"
        "for the perturbed input (without doing any unnecessary operations,\n"
        "for efficiency reasons). If unperturbed spectra (and analytical\n"
        "jacobians) are calculated with *RteCalc*, the standard choice for\n"
        "this agenda should be *RteCalcNoJacobians*.\n"
       ),
       OUTPUT( "y" ),
       INPUT( "f_grid", "vmr_field", "t_field", "sensor_los" )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "main_agenda" ),
       DESCRIPTION
       (
        "The agenda corresponding to the entire controlfile. This is\n" 
        "executed when ARTS is run.\n"
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
        "then add the method *VectorToTbByPlanck*.\n"
       ),
       OUTPUT( "y" ),
       INPUT("t_field_raw", "vmr_field_raw", "z_field_raw", "pnd_field_raw",
             "p_grid", "sensor_los", "cloudbox_on", "cloudbox_limits",
             "z_surface")));

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
       OUTPUT( "ext_mat", "abs_vec" ),
       INPUT( "ext_mat", "abs_vec", "f_index", "abs_scalar_gas" )));


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
        "   scat_p_index : Position. \n"
        "   scat_lat_index : Position. \n"
        "   scat_lon_index : Position. \n"
        ),
       OUTPUT( "ext_mat", "abs_vec" ),
       INPUT( "ext_mat", "abs_vec", 
              "ext_mat_spt", "abs_vec_spt",
              "scat_p_index", "scat_lat_index",
              "scat_lon_index" )));

 agenda_data.push_back
    (AgRecord
     ( NAME( "pha_mat_spt_agenda" ),
       DESCRIPTION
       (
        "Calculates the phase matrix for a single particle type.\n"
        "\n"
        "Different options are possible for the usage of this agenda: \n"
        "*pha_mat_sptFromData* or *pha_mat_sptDOITOpt*. \n"
        ),
       OUTPUT( "pha_mat_spt"),
       INPUT( "pha_mat_spt", "scat_za_index", "scat_lat_index", "scat_lon_index",
              "scat_p_index", "scat_aa_index", "rte_temperature")));
       
  agenda_data.push_back
    (AgRecord
     ( NAME( "ppath_agenda" ),
       DESCRIPTION
       (
        "Calculation of complete propagation paths.\n"
        "\n"
        "In contrast to *ppath_step_agenda* that controls the ray tracing\n"
        "inside each grid box, this agenda determines how complete paths are\n"
        "determined. The standard choice is to do this in a step-by-step\n"
        "manner using *ppath_step_agenda*, with this agenda set to call\n" 
        "*ppathStepByStep*.\n" 
        "\n"
        "The WSV *rte_los* is both input and output as in some cases it is\n"
        "determined as part of the propagation path calculations (such as for"
        "for radio link calculations).\n"
        ),
       OUTPUT( "ppath" ),
       INPUT( "rte_pos", "rte_los", "cloudbox_on", "ppath_inside_cloudbox_do", 
              "mblock_index", "t_field", "z_field", "vmr_field",
              "edensity_field", "f_index" )));

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
        "initialisation is not handled by *ppath_step_agenda* (but by \n"
        "the internal function ppath_start_stepping). \n"
        "\n"
        "The *ppath_step_agenda* put in points along the propagation path \n"
        "at all crossings with the grids, tangent points and points of \n"
        "surface reflection. It is also allowed to make agendas that put in \n"
        "additional points to fulfil some criterion, such as a maximum \n"
        "distance along the path between the points. Accordingly, the \n"
        "number of new points of each step can exceed one.\n"
        ),
       OUTPUT( "ppath_step" ),
       INPUT( "ppath_step", "t_field", "z_field", "vmr_field", 
              "edensity_field", "f_index" )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "refr_index_agenda" ),
       DESCRIPTION
       (
        "Calculation of the refractive index of air.\n"
        "\n"
        "This agenda should calculate the summed refractive index for all\n"
        "relevant atmospheric constituents, with respect to both phase and\n"
        "group velocity.\n"
        ),
       OUTPUT( "refr_index", "refr_index_group" ),
       INPUT(  "f_index", "rte_pressure", "rte_temperature", "rte_vmr_list", 
               "rte_edensity" )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "sensor_response_agenda" ),
       DESCRIPTION
       (
        "The sensor response data for present measurement block.\n"
        "\n"
        "This agenda shall provide *sensor_response* and associated variables\n"
        "for the present measurement block (*mblock_index*).\n"
        ),
       OUTPUT( "sensor_response", "sensor_response_f", "sensor_response_pol",
               "sensor_response_za", "sensor_response_aa" ),
       INPUT(  "mblock_index" )));

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
        ),
       OUTPUT( "ext_mat_spt", "abs_vec_spt"),
       INPUT(  "ext_mat_spt", "abs_vec_spt",
               "scat_p_index", "scat_lat_index", "scat_lon_index",
               "rte_temperature", "scat_za_index", "scat_aa_index"
               )));

  agenda_data.push_back
    (AgRecord
     ( NAME( "surface_rtprop_agenda" ),
       DESCRIPTION
       (
        "Surface radiative properties. \n"
        "\n"
        "See the user guide for closer definitions of the variables \n"
        "that describe the surface radiative properties. These variables are:\n"
        "   *surface_emission*, *surface_los* and *surface_rmatrix* \n"
        ),
       OUTPUT( "surface_emission", "surface_los", "surface_rmatrix" ),
       INPUT( "rte_pos", "rte_los" )));

 agenda_data.push_back
    (AgRecord
     ( NAME( "ybatch_calc_agenda" ),
       DESCRIPTION
       (
        "Calculations to perform for each batch case.\n"
        "\n"
        "Must produce a new spectrum vector (*y*) and Jacobi matrix (*jacobian*).\n"
        "See further *ybatchCalc*.\n"
        ),
       OUTPUT( "y", "jacobian" ),
       INPUT( "ybatch_index" )));

  agenda_data.push_back
  (AgRecord
   ( NAME( "test_agenda" ),
     DESCRIPTION
     (
      "Dummy agenda for testing purposes.\n"
     ),
     OUTPUT(),
     INPUT()));
  

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
