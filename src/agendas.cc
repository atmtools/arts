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
#define OUTPUT(...) \
  { __VA_ARGS__ }
#define INPUT(...) \
  { __VA_ARGS__ }

/*! The lookup information for the agendas. */
namespace global_data {
Array<AgRecord> agenda_data;
} // namespace global_data

void define_agenda_data() {
  using global_data::agenda_data;

  // Initialize to zero, just in case:
  agenda_data.resize(0);

  /*----------------------------------------------------------------------
    Agendas must be put in in alphabetical order. 
    No distinction is made between uppercase and lowercase letters. 
    The sign "_" comes after all letters.
    ----------------------------------------------------------------------*/

  agenda_data.push_back(AgRecord(
      NAME("propmat_clearsky_agenda"),
      DESCRIPTION(
          "Calculate the absorption coefficient matrix.\n"
          "\n"
          "This agenda calculates the absorption coefficient matrix for all\n"
          "absorption species as a function of the given atmospheric state for\n"
          "one point in the atmosphere. The result is returned in\n"
          "*propmat_clearsky*. The atmospheric state has to be specified by\n"
          "*rtp_pressure*, *rtp_temperature*, *rtp_mag*, and *rtp_vmr*.\n"
          "\n"
          "The methods inside this agenda may require a lot of additional\n"
          "input variables, such as *abs_species*, etc.\n"
          "\n"
          "The include file 'agendas.arts' predefines some possible agendas\n"
          "that can be used here.\n"),
      OUTPUT("propmat_clearsky",
             "nlte_source",
             "dpropmat_clearsky_dx",
             "dnlte_source_dx"),
      INPUT("jacobian_quantities",
            "select_abs_species",
            "f_grid",
            "rtp_mag",
            "rtp_los",
            "rtp_pressure",
            "rtp_temperature",
            "rtp_nlte",
            "rtp_vmr")));

  agenda_data.push_back(
      AgRecord(NAME("dobatch_calc_agenda"),
               DESCRIPTION("Calculations to perform for each batch case.\n"
                           "\n"
                           "See further *DOBatchCalc*.\n"),
               OUTPUT("spectral_radiance_field",
                      "radiance_field",
                      "irradiance_field",
                      "spectral_irradiance_field"),
               INPUT("ybatch_index")));

  agenda_data.push_back(AgRecord(
      NAME("doit_conv_test_agenda"),
      DESCRIPTION(
          "Compute the convergence test.\n"
          "\n"
          "The method *cloudbox_field_monoIterate* solves the VRTE iteratively."
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
          "  because result can be inaccurate.\n"),
      OUTPUT("doit_conv_flag", "doit_iteration_counter"),
      INPUT("doit_conv_flag",
            "doit_iteration_counter",
            "cloudbox_field_mono",
            "cloudbox_field_mono_old")));

  agenda_data.push_back(AgRecord(
      NAME("doit_mono_agenda"),
      DESCRIPTION(
          "Performs monochromatic DOIT calculation."
          "\n"
          "This agenda includes for example the following methods:\n"
          "   1. *DoitScatteringDataPrepare* \n"
          "   2. *cloudbox_field_monoIterate*\n"
          "\n"
          "The result of the agenda is the radiation field inside the \n"
          "cloudbox and on the cloudbox boundary, which can be used \n"
          "as radiative background for a clearsky radiative transfer \n"
          "calculation. \n"
          "\n"
          "See the Arts online documentation\n"
          "for more information about the methods.\n"),
      OUTPUT("cloudbox_field_mono"),
      INPUT("cloudbox_field_mono", "f_grid", "f_index")));

  agenda_data.push_back(AgRecord(
      NAME("doit_scat_field_agenda"),
      DESCRIPTION(
          "Calculation of the scattering integral field (DOIT). \n"
          "\n"
          "This agenda is called repeatedly in each DOIT iteration.\n"
          "The following methods can be used for calculating the \n"
          "scattering integral field: \n"
          "\n"
          "*doit_scat_fieldCalc*: This method calculates the scattering \n"
          "  integral field by using the angular grids *za_grid* \n"
          "  and *aa_grid*, which are also used in the update of the \n"
          "  radiation field (*doit_rte_agenda*).\n"
          "\n"
          "*doit_scat_fieldCalcLimb*: This method calculates the scattering \n"
          "  integral field.  The difference to the previous method is that \n"
          "  the data is interpolated on equidistant angular grids. \n"
          "  Especially for limb, where a very fine zenith angle grid \n"
          "  resolution is required for the RT transfer part, this method \n"
          "  is much faster than *doit_scat_fieldCalc*. \n"),
      OUTPUT("doit_scat_field"),
      INPUT("doit_scat_field", "cloudbox_field_mono")));

  agenda_data.push_back(AgRecord(
      NAME("doit_rte_agenda"),
      DESCRIPTION(
          "Radiative transfer calculations in cloudbox.\n"
          "\n"
          "Agenda for radiative transfer step calculations with \n"
          "fixed scattering integral term shoul be specified here.\n"
          "Output is the updated radiation field in the cloudbox. \n"
          "This agenda is called repeatedly in each DOIT iteration.\n"
          "\n"
          "Normally one should use\n"
          "*cloudbox_fieldUpdateSeq1D* or *cloudbox_fieldUpdateSeq3D*:\n"
          "Seqential update of the radiation field.\n"
          "   This method is the fastest and most accurate method.\n"
          "\n"
          "A very similar method in plane parallel approximation is\n"
          "*cloudbox_fieldUpdateSeq1DPP*:\n"
          "   This method also includes the sequential update and is slightly\n"
          "   faster than the above one. The drawback is that it is less\n"
          "   accurate, especially for limb geometries and large off-nadir\n"
          "   viewing angles.\n"
          "\n"
          "The following method was used before the sequential update\n"
          "was invented. It is very slow and should therefore only \n"
          "be used for test cases.\n"
          "*cloudbox_fieldUpdate1D*: Old method.\n"),
      OUTPUT("cloudbox_field_mono"),
      INPUT("cloudbox_field_mono", "doit_scat_field")));

  agenda_data.push_back(AgRecord(
      NAME("forloop_agenda"),
      DESCRIPTION(
          "The body for a for loop.\n"
          "\n"
          "This agenda contains the body of the for loop to be execute by the\n"
          "method *ForLoop*.\n"),
      OUTPUT(),
      INPUT("forloop_index")));

  /*
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
  */

  agenda_data.push_back(
      AgRecord(NAME("gas_scattering_agenda"),
               DESCRIPTION("Calculation of the gas scattering extinction and phase matrix.\n"
                           "\n"
                           "This agenda calculates the gas scattering cross\n"
                           "section and the normalized phase matrix for a specific\n"
                           "incoming ( *gas_scattering_los_in* ) and outgoing (*gas_scattering_los_out*) direction.\n"
                           "The scattering cross section is calculated along a\n"
                           "propagtion path given by the propagation path variables\n"
                           "*rtp_pressure*, *rtp_temperature*, and *rtp_vmr*."
                           "If *gas_scattering_los_in* and *gas_scattering_los_out* are empty vectors, then\n"
                           "*gas_scattering_mat* is set empty. If *gas_scattering_los_in* and *gas_scattering_los_out*\n"
                           "are not empty, then the phase matrix is calculated\n"
                           "for the define incoming and outgoing direction.\n"),
               OUTPUT("gas_scattering_coef","gas_scattering_mat","gas_scattering_fct_legendre"),
               INPUT("f_grid",
                     "rtp_pressure",
                     "rtp_temperature",
                     "rtp_vmr",
                     "gas_scattering_los_in",
                     "gas_scattering_los_out",
                     "gas_scattering_output_type")));

  agenda_data.push_back(
      AgRecord(NAME("g0_agenda"),
               DESCRIPTION("Calculation of the gravity at zero altitude.\n"
                           "\n"
                           "Returns *g0* for given geographical position.\n"),
               OUTPUT("g0"),
               INPUT("lat", "lon")));

  agenda_data.push_back(AgRecord(
      NAME("inversion_iterate_agenda"),
      DESCRIPTION(
          "Work in progress ...\n"
          "\n"
          "The WSV *jacobian* is both in- and output. As input variable, *jacobian*\n"
          "is assumed to be valid for the previous iteration. For the first iteration\n"
          "the input *jacobian* shall be set to have size zero, to flag that there\n"
          "is not yet any calculated Jacobian.\n"),
      OUTPUT("yf", "jacobian"),
      INPUT("x", "jacobian_do", "inversion_iteration_counter")));

  agenda_data.push_back(AgRecord(
      NAME("iy_cloudbox_agenda"),
      DESCRIPTION(
          "Intensity at boundary or interior of the cloudbox.\n"
          "\n"
          "The task of the agenda is to determine the intensity at some point\n"
          "at the boundary of or inside the cloudbox. The actual calculations\n"
          "inside the agenda differ depending on scattering solution method.\n"
          "If DOIT is used, an interpolating of the intensity field should be\n"
          "performed. Another option is to start backward Monte Carlo \n"
          "calculations from this point.\n"
          "\n"
          "A function calling this agenda shall set *rte_pos* and *rte_los* to\n"
          "the position and line-of-sight for which the scattered radiation\n"
          "shall be determined.\n"
          "\n"
          "The include-file 'agendas.arts' pre-defines some agendas that can\n"
          "either be used directly, or serve as examples.\n"),
      OUTPUT("iy"),
      INPUT("f_grid", "rtp_pos", "rtp_los")));

  agenda_data.push_back(AgRecord(
      NAME("iy_independent_beam_approx_agenda"),
      DESCRIPTION(
          "Agenda dedicated to *iyIndependentBeamApproximation*.\n"
          "\n"
          "If *iyIndependentBeamApproximation* is used, this agenda basically\n"
          "replaces *iy_main_agenda*. Accordingly, this agenda has exactly the\n"
          "same output as *iy_main_agenda*.\n"),
      OUTPUT("iy", "iy_aux", "ppath", "diy_dx"),
      INPUT("diy_dx",
            "iy_agenda_call1",
            "iy_unit",
            "iy_transmittance",
            "iy_aux_vars",
            "iy_id",
            "atmosphere_dim",
            "p_grid",
            "lat_grid",
            "lon_grid",
            "lat_true",
            "lon_true",
            "t_field",
            "z_field",
            "vmr_field",
            "z_surface",
            "ppath_lmax",
            "ppath_lraytrace",
            "cloudbox_on",
            "cloudbox_limits",
            "pnd_field",
            "jacobian_do",
            "f_grid",
            "rte_pos",
            "rte_los",
            "rte_pos2")));

  agenda_data.push_back(AgRecord(
      NAME("iy_loop_freqs_agenda"),
      DESCRIPTION(
          "Agenda dedicated to *iyLoopFrequencies*.\n"
          "\n"
          "If *iyLoopFrequencies* is used, this agenda basically replaces\n"
          "*iy_main_agenda*.Accordingly, this agenda has exactly the same\n"
          "output as *iy_main_agenda*.\n"),
      OUTPUT("iy", "iy_aux", "ppath", "diy_dx"),
      INPUT("diy_dx",
            "iy_agenda_call1",
            "iy_transmittance",
            "iy_aux_vars",
            "iy_id",
            "f_grid",
            "rte_pos",
            "rte_los",
            "rte_pos2")));

  agenda_data.push_back(AgRecord(
      NAME("iy_main_agenda"),
      DESCRIPTION(
          "Calculation of a single monochromatic pencil beam spectrum.\n"
          "\n"
          "The task of the agenda is to calculate the monochromatic pencil beam\n"
          "spectrum for the position specified by *rte_pos* and the viewing\n"
          "direction specified by *rte_los*.\n"
          "\n"
          "Methods for this agenda can either handle the complete calculation,\n"
          "make use of e.g. *iy_cloudbox_agenda* or be restricted to special\n"
          "cases. See the documentation for the different methods.\n"
          "\n"
          "The include-file 'agendas.arts' predefines some typical alternatives\n"
          "that can be used directly, or adapted for specific applications.\n"),
      OUTPUT("iy", "iy_aux", "ppath", "diy_dx", "geo_pos"),
      INPUT("diy_dx",
            "iy_agenda_call1",
            "iy_transmittance",
            "iy_aux_vars",
            "iy_id",
            "iy_unit",
            "cloudbox_on",
            "jacobian_do",
            "f_grid",
            "nlte_field",
            "rte_pos",
            "rte_los",
            "rte_pos2")));

  agenda_data.push_back(AgRecord(
      NAME("iy_radar_agenda"),
      DESCRIPTION(
          "Calculation of pointwise backscattering.\n"
          "\n"
          "This agenda has a similar role for *yRadar* as *iy_main_agenda*.\n"
          "for *yCalc*.\n"),
      OUTPUT("iy", "iy_aux", "ppath", "diy_dx", "geo_pos"),
      INPUT("iy_aux_vars",
            "iy_id",
            "cloudbox_on",
            "jacobian_do",
            "rte_pos",
            "rte_los")));

  agenda_data.push_back(AgRecord(
      NAME("iy_space_agenda"),
      DESCRIPTION(
          "Downwelling radiation at the top of the atmosphere.\n"
          "\n"
          "Possible terms to include in this agenda include cosmic background\n"
          "radiation and solar radiation.\n"
          "\n"
          "A function calling this agenda shall set *rtp_pos* and *rtp_los* to\n"
          "the position and line-of-sight for which the entering radiation \n"
          "shall be determined. The position and line-of-sight must be known, \n"
          "for example, when radiation from the sun is considered.\n"
          "\n"
          "The include-file 'agendas.arts' predefines an agenda that can be\n"
          "applied directly for most users.\n"),
      OUTPUT("iy"),
      INPUT("f_grid", "rtp_pos", "rtp_los")));

  agenda_data.push_back(AgRecord(
      NAME("iy_surface_agenda"),
      DESCRIPTION(
          "Upwelling radiation from the surface.\n"
          "\n"
          "The task of the agenda is to determine the upwelling intensity from\n"
          "the surface, for given point and direction.\n"
          "\n"
          "The standard choice should be to make use of *surface_rtprop_agenda*\n"
          "through the WSM *iySurfaceRtpropAgenda*.\n"
          "\n"
          "A function calling this agenda shall set *rtp_pos* and *rtp_los* to\n"
          "the position and line-of-sight for which the upwelling radiation\n"
          "shall be determined.\n"
          "\n"
          "See also the include-file 'agendas.arts' for a predefined agenda\n"
          "suitable to be used in most applications.\n"),
      OUTPUT("iy", "diy_dx"),
      INPUT("diy_dx",
            "dsurface_rmatrix_dx",
            "dsurface_emission_dx",
            "iy_unit",
            "iy_transmittance",
            "iy_id",
            "cloudbox_on",
            "jacobian_do",
            "iy_main_agenda",
            "f_grid",
            "nlte_field",
            "rtp_pos",
            "rtp_los",
            "rte_pos2",
            "surface_props_data",
            "dsurface_names")));

  agenda_data.push_back(AgRecord(
      NAME("jacobian_agenda"),
      DESCRIPTION(
          "Pure numerical Jacobian calculations.\n"
          "\n"
          "Parts of the Jacobian matrix can be determined by (semi-)analytical\n"
          "expressions, while other parts are calculated in a pure numerical\n"
          "manner (by perturbations). This agenda describes the calculations to\n"
          "be performed in the later case.\n"
          "\n"
          "This agenda is normally not set directly by the user, but is created\n"
          "by calling the the jacobianAdd set of methods.\n"),
      OUTPUT("jacobian"),
      INPUT("jacobian", "mblock_index", "iyb", "yb")));

  agenda_data.push_back(AgRecord(
      NAME("main_agenda"),
      DESCRIPTION(
          "The agenda corresponding to the entire controlfile. This is\n"
          "executed when ARTS is run.\n"),
      OUTPUT(),
      INPUT()));

  agenda_data.push_back(AgRecord(
      NAME("met_profile_calc_agenda"),
      DESCRIPTION(
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
          "*DoitInit*\n"
          "*DoitGetIncoming*\n"
          "*cloudbox_fieldSetClearsky*\n"
          "*DoitCalc*\n"
          "*yCalc*\n"),
      OUTPUT("y"),
      INPUT("t_field_raw",
            "vmr_field_raw",
            "z_field_raw",
            "pnd_field_raw",
            "p_grid",
            "sensor_los",
            "cloudbox_on",
            "cloudbox_limits",
            "z_surface")));

  agenda_data.push_back(AgRecord(
      NAME("pha_mat_spt_agenda"),
      DESCRIPTION(
          "Calculates the phase matrix for individual scattering elements.\n"
          "\n"
          "Different options are possible for the usage of this agenda: \n"
          "*pha_mat_sptFromData* or *pha_mat_sptFromDataDOITOpt*. \n"),
      OUTPUT("pha_mat_spt"),
      INPUT("pha_mat_spt",
            "za_index",
            "scat_lat_index",
            "scat_lon_index",
            "scat_p_index",
            "aa_index",
            "rtp_temperature")));

  agenda_data.push_back(AgRecord(
      NAME("pnd_agenda_array"),
      DESCRIPTION(
          "Returns particle number density data for each scattering species.\n"
          "\n"
          "This variable is used when mapping data in *particle_bulkprop_field*\n"
          "to *pnd_field*. The variable is also necessary when calculating\n"
          "scattering species weighting functions.\n"
          "\n"
          "Note that content of this agenda array, *scat_species* and\n"
          "*pnd_agenda_array_input_names* must be consistent.\n"),
      OUTPUT("pnd_data", "dpnd_data_dx"),
      INPUT("agenda_array_index",
            "pnd_agenda_input_t",
            "pnd_agenda_input",
            "pnd_agenda_input_names",
            "dpnd_data_dx_names")));

  agenda_data.push_back(AgRecord(
      NAME("ppath_agenda"),
      DESCRIPTION(
          "Calculation of complete propagation paths.\n"
          "\n"
          "In contrast to *ppath_step_agenda* that controls the ray tracing\n"
          "inside each grid box, this agenda determines how complete paths are\n"
          "determined. The standard choice is to do this in a step-by-step\n"
          "manner using *ppath_step_agenda*, with this agenda set to call\n"
          "*ppathStepByStep*.\n"
          "\n"
          "The WSV *rte_los* is both input and output as in some cases it is\n"
          "determined as part of the propagation path calculations (such as\n"
          "radio link calculations).\n"),
      OUTPUT("ppath"),
      INPUT("ppath_lmax",
            "ppath_lraytrace",
            "rte_pos",
            "rte_los",
            "rte_pos2",
            "cloudbox_on",
            "ppath_inside_cloudbox_do",
            "f_grid")));

  agenda_data.push_back(AgRecord(
      NAME("ppath_step_agenda"),
      DESCRIPTION(
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
          "\n"
          "The include file 'agendas.arts' defines some agendas that can be\n"
          "used here."),
      OUTPUT("ppath_step"),
      INPUT("ppath_step", "ppath_lmax", "ppath_lraytrace", "f_grid")));

  agenda_data.push_back(AgRecord(
      NAME("refr_index_air_agenda"),
      DESCRIPTION(
          "Calculation of the refractive index of air.\n"
          "\n"
          "This agenda should calculate the summed refractive index for all\n"
          "relevant atmospheric constituents, with respect to both phase and\n"
          "group velocity.\n"),
      OUTPUT("refr_index_air", "refr_index_air_group"),
      INPUT("rtp_pressure", "rtp_temperature", "rtp_vmr", "f_grid")));
  
  agenda_data.push_back(AgRecord(
      NAME("refr_index_air_ZZZ_agenda"),
      DESCRIPTION(
          "Calculation of the refractive index of air.\n"
          "\n"
          "This agenda shall the refractive index for given atmospheric position.\n"),
      OUTPUT("refr_index_air", "refr_index_air_group"),
      INPUT("rtp_pos")));

  agenda_data.push_back(AgRecord(
      NAME("sensor_response_agenda"),
      DESCRIPTION(
          "This agenda shall provide *sensor_response* and associated variables.\n"
          "\n"
          "So far only required when doing inversions involving some sensor variables.\n"),
      OUTPUT("sensor_response",
             "sensor_response_f",
             "sensor_response_f_grid",
             "sensor_response_pol",
             "sensor_response_pol_grid",
             "sensor_response_dlos",
             "sensor_response_dlos_grid",
             "mblock_dlos"),
      INPUT("f_backend")));

  agenda_data.push_back(AgRecord(
      NAME("spt_calc_agenda"),
      DESCRIPTION(
          "Calculates single scattering properties for individual scattering\n"
          "elements from the amplitude matrix.\n"
          "\n"
          "This agenda sets up the methods, which should be used to calculate \n"
          "the single scattering properties, i.e. the extinction matrix and the \n"
          "absorbtion vector.\n "
          "\n"
          "Normally you  use:\n"
          "*opt_prop_sptFromMonoData*\n"),
      OUTPUT("ext_mat_spt", "abs_vec_spt"),
      INPUT("ext_mat_spt",
            "abs_vec_spt",
            "scat_p_index",
            "scat_lat_index",
            "scat_lon_index",
            "rtp_temperature",
            "za_index",
            "aa_index")));

  agenda_data.push_back(AgRecord(
      NAME("surface_rtprop_agenda"),
      DESCRIPTION(
          "Provides radiative properties of the surface. \n"
          "\n"
          "Provides surface emission and surface reflection coefficient matrix\n"
          "(see user guide for closer definitions of the respective variables\n"
          "*surface_emission*, *surface_los*, and *surface_rmatrix*) according\n"
          "to the characteristics of the surface specified by the methods called\n"
          "within the agenda. Typical methods include *surfaceBlackbody*,\n"
          "*surfaceFlatScalarReflectivity*, *surfaceFlatReflectivity*,\n"
          "*surfaceFlatRefractiveIndex*, and *surfaceLambertianSimple*.\n"),
      OUTPUT("surface_skin_t",
             "surface_emission",
             "surface_los",
             "surface_rmatrix"),
      INPUT("f_grid", "rtp_pos", "rtp_los")));

  agenda_data.push_back(AgRecord(
      NAME("surface_rtprop_agenda_array"),
      DESCRIPTION(
          "Description of surface radiative properties, for each surface type.\n"
          "\n"
          "Each of these agendas shall treat the radiative properties of a\n"
          "surface type. The task of these agendas is equivalent to that of\n"
          "*surface_rtprop_agenda*.\n"
          "\n"
          "The order of the agendas shall match the coding used in\n"
          "*surface_type_mask*.\n"),
      OUTPUT("surface_skin_t",
             "surface_emission",
             "surface_los",
             "surface_rmatrix"),
      INPUT("agenda_array_index",
            "f_grid",
            "rtp_pos",
            "rtp_los")));

  agenda_data.push_back(
      AgRecord(NAME("test_agenda"),
               DESCRIPTION("Dummy agenda for testing purposes.\n"),
               OUTPUT(),
               INPUT()));

  agenda_data.push_back(
      AgRecord(NAME("test_agenda_array"),
               DESCRIPTION("Agenda array for TestArrayOfAgenda test case.\n"),
               OUTPUT(),
               INPUT("agenda_array_index", "iy_unit")));

  agenda_data.push_back(AgRecord(
      NAME("water_p_eq_agenda"),
      DESCRIPTION("Calculation of the saturation pressure of water.\n"),
      OUTPUT("water_p_eq_field"),
      INPUT("t_field")));

  agenda_data.push_back(AgRecord(
      NAME("ybatch_calc_agenda"),
      DESCRIPTION(
          "Calculations to perform for each batch case.\n"
          "\n"
          "Must produce a new spectrum vector (*y*) and Jacobi matrix (*jacobian*).\n"
          "See further *ybatchCalc*.\n"),
      OUTPUT("y", "y_aux", "jacobian"),
      INPUT("ybatch_index")));

  std::sort(agenda_data.begin(), agenda_data.end(), [](auto& a, auto& b) {
    return a.Name() < b.Name();
  });
}
