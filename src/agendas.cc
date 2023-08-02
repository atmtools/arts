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
      DESCRIPTION(R"--(Calculate the absorption coefficient matrix.

This agenda calculates the absorption coefficient matrix for all
absorption species as a function of the given atmospheric state for
one point in the atmosphere. The result is returned in
*propmat_clearsky*. The atmospheric state has to be specified by
*rtp_pressure*, *rtp_temperature*, ``rtp_mag``, and *rtp_vmr*.

The methods inside this agenda may require a lot of additional
input variables, such as *abs_species*, etc.

The include file 'agendas.arts' predefines some possible agendas
that can be used here.
)--"),
      OUTPUT("propmat_clearsky",
             "nlte_source",
             "dpropmat_clearsky_dx",
             "dnlte_source_dx"),
      INPUT("jacobian_quantities",
            "select_abs_species",
            "f_grid",
            "rtp_los",
            "atm_point")));

  agenda_data.push_back(
      AgRecord(NAME("dobatch_calc_agenda"),
               DESCRIPTION(R"--(Calculations to perform for each batch case.

See also: *DOBatchCalc*
)--"),
               OUTPUT("spectral_radiance_field",
                      "radiance_field",
                      "irradiance_field",
                      "spectral_irradiance_field"),
               INPUT("ybatch_index")));

  agenda_data.push_back(AgRecord(
      NAME("doit_conv_test_agenda"),
      DESCRIPTION(R"--(Compute the convergence test.

The method ``cloudbox_field_monoIterate`` solves the VRTE iteratively.
This method requires 
a convergence test. The user can choose different convergence tests
which are to be defined in this agenda.

Possible workspace methods are:
 - *doit_conv_flagAbs*: Calculates the absolute differences 
   for each Stokes component separately.
 - *doit_conv_flagAbsBT*: Same as above, but the convergence limit
   can be specified in Kelvin BT (Rayleigh Jeans).
 - *doit_conv_flagLsq*: Least square convergence test. Not recommended
   because result can be inaccurate.
)--"),
      OUTPUT("doit_conv_flag", "doit_iteration_counter"),
      INPUT("doit_conv_flag",
            "doit_iteration_counter",
            "cloudbox_field_mono",
            "cloudbox_field_mono_old")));

  agenda_data.push_back(AgRecord(
      NAME("doit_mono_agenda"),
      DESCRIPTION(R"--(Performs monochromatic DOIT calculation.

This agenda includes for example the following methods:

1. *DoitScatteringDataPrepare* 
2. ``cloudbox_field_monoIterate``

The result of the agenda is the radiation field inside the 
cloudbox and on the cloudbox boundary, which can be used 
as radiative background for a clearsky radiative transfer 
calculation.

See the Arts online documentation
for more information about the methods.
)--"),
      OUTPUT("cloudbox_field_mono"),
      INPUT("cloudbox_field_mono", "f_grid", "f_index")));

  agenda_data.push_back(AgRecord(
      NAME("doit_scat_field_agenda"),
      DESCRIPTION(R"--(Calculation of the scattering integral field (DOIT). 

This agenda is called repeatedly in each DOIT iteration.
The following methods can be used for calculating the 
scattering integral field: 

-  *doit_scat_fieldCalc*: This method calculates the scattering 
   integral field by using the angular grids *za_grid* 
   and *aa_grid*, which are also used in the update of the 
   radiation field (*doit_rte_agenda*).

-  *doit_scat_fieldCalcLimb*: This method calculates the scattering 
   integral field.  The difference to the previous method is that 
   the data is interpolated on equidistant angular grids. 
   Especially for limb, where a very fine zenith angle grid 
   resolution is required for the RT transfer part, this method 
   is much faster than *doit_scat_fieldCalc*.
)--"),
      OUTPUT("doit_scat_field"),
      INPUT("doit_scat_field", "cloudbox_field_mono")));

  agenda_data.push_back(AgRecord(
      NAME("doit_rte_agenda"),
      DESCRIPTION(R"--(Radiative transfer calculations in cloudbox.

Agenda for radiative transfer step calculations with 
fixed scattering integral term shoul be specified here.
Output is the updated radiation field in the cloudbox. 
This agenda is called repeatedly in each DOIT iteration.

Normally one should use

 - ``cloudbox_fieldUpdateSeq1D`` or ``cloudbox_fieldUpdateSeq3D``:
    Seqential update of the radiation field.
    This method is the fastest and most accurate method.

A very similar method in plane parallel approximation is
 - ``cloudbox_fieldUpdateSeq1DPP``:
    This method also includes the sequential update and is slightly
    faster than the above one. The drawback is that it is less
    accurate, especially for limb geometries and large off-nadir
    viewing angles.

The following method was used before the sequential update
was invented. It is very slow and should therefore only 
be used for test cases:

 - ``cloudbox_fieldUpdate1D``: Old method.
)--"),
      OUTPUT("cloudbox_field_mono"),
      INPUT("cloudbox_field_mono", "doit_scat_field")));

  agenda_data.push_back(AgRecord(
      NAME("forloop_agenda"),
      DESCRIPTION(R"--(The body for a for loop.

This agenda contains the body of the for loop to be execute by the
method *ForLoop*.
)--"),
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
               DESCRIPTION(R"--(Calculation of the gas scattering extinction and phase matrix.

This agenda calculates the gas scattering cross
section and the normalized phase matrix for a specific
incoming ( *gas_scattering_los_in* ) and outgoing (*gas_scattering_los_out*) direction.
The scattering cross section is calculated along a
propagtion path given by the propagation path variables
*rtp_pressure*, *rtp_temperature*, and *rtp_vmr*.
If *gas_scattering_los_in* and *gas_scattering_los_out* are empty vectors, then
*gas_scattering_mat* is set empty. If *gas_scattering_los_in* and *gas_scattering_los_out*
are not empty, then the phase matrix is calculated
for the define incoming and outgoing direction.
)--"),
               OUTPUT("gas_scattering_coef","gas_scattering_mat","gas_scattering_fct_legendre"),
               INPUT("f_grid",
                     "atm_point",
                     "gas_scattering_los_in",
                     "gas_scattering_los_out",
                     "gas_scattering_output_type")));

  agenda_data.push_back(
      AgRecord(NAME("g0_agenda"),
               DESCRIPTION(R"--(Calculation of the gravity at zero altitude.

Returns *g0* for given geographical position.
)--"),
               OUTPUT("g0"),
               INPUT("lat", "lon")));

  agenda_data.push_back(AgRecord(
      NAME("inversion_iterate_agenda"),
      DESCRIPTION(R"--(Work in progress ...

The WSV *jacobian* is both in- and output. As input variable, *jacobian*
is assumed to be valid for the previous iteration. For the first iteration
the input *jacobian* shall be set to have size zero, to flag that there
is not yet any calculated Jacobian.
)--"),
      OUTPUT("yf", "jacobian"),
      INPUT("x", "jacobian_do", "inversion_iteration_counter")));

  agenda_data.push_back(AgRecord(
      NAME("iy_cloudbox_agenda"),
      DESCRIPTION(R"--(Intensity at boundary or interior of the cloudbox.

The task of the agenda is to determine the intensity at some point
at the boundary of or inside the cloudbox. The actual calculations
inside the agenda differ depending on scattering solution method.
If DOIT is used, an interpolating of the intensity field should be
performed. Another option is to start backward Monte Carlo 
calculations from this point.

A function calling this agenda shall set *rte_pos* and *rte_los* to
the position and line-of-sight for which the scattered radiation
shall be determined.

The include-file 'agendas.arts' pre-defines some agendas that can
either be used directly, or serve as examples.
)--"),
      OUTPUT("iy"),
      INPUT("f_grid", "rtp_pos", "rtp_los")));

  agenda_data.push_back(AgRecord(
      NAME("iy_independent_beam_approx_agenda"),
      DESCRIPTION(R"--(Agenda dedicated to *iyIndependentBeamApproximation*.

If *iyIndependentBeamApproximation* is used, this agenda basically
replaces *iy_main_agenda*. Accordingly, this agenda has exactly the
same output as *iy_main_agenda*.
)--"),
      OUTPUT("iy", "iy_aux", "ppath", "diy_dx"),
      INPUT("diy_dx",
            "iy_agenda_call1",
            "iy_unit",
            "iy_transmittance",
            "iy_aux_vars",
            "iy_id",
            "lat_true",
            "lon_true",
            "atm_field",
            "surface_field",
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
      DESCRIPTION(R"--(Agenda dedicated to *iyLoopFrequencies*.

If *iyLoopFrequencies* is used, this agenda basically replaces
*iy_main_agenda*.Accordingly, this agenda has exactly the same
output as *iy_main_agenda*.
)--"),
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
      DESCRIPTION(R"--(Calculation of a single monochromatic pencil beam spectrum.

The task of the agenda is to calculate the monochromatic pencil beam
spectrum for the position specified by *rte_pos* and the viewing
direction specified by *rte_los*.

Methods for this agenda can either handle the complete calculation,
make use of e.g. *iy_cloudbox_agenda* or be restricted to special
cases. See the documentation for the different methods.

The include-file 'agendas.arts' predefines some typical alternatives
that can be used directly, or adapted for specific applications.
)--"),
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
            "atm_field",
            "rte_pos",
            "rte_los",
            "rte_pos2")));

  agenda_data.push_back(AgRecord(
      NAME("iy_radar_agenda"),
      DESCRIPTION(R"--(Calculation of pointwise backscattering.

This agenda has a similar role for *yRadar* as *iy_main_agenda*.
for *yCalc*.
)--"),
      OUTPUT("iy", "iy_aux", "ppath", "diy_dx", "geo_pos"),
      INPUT("iy_aux_vars",
            "iy_id",
            "cloudbox_on",
            "jacobian_do",
            "rte_pos",
            "rte_los")));

  agenda_data.push_back(AgRecord(
      NAME("iy_space_agenda"),
      DESCRIPTION(R"--(Downwelling radiation at the top of the atmosphere.

Possible terms to include in this agenda include cosmic background
radiation and solar radiation.

A function calling this agenda shall set *rtp_pos* and *rtp_los* to
the position and line-of-sight for which the entering radiation 
shall be determined. The position and line-of-sight must be known, 
for example, when radiation from the sun is considered.

The include-file 'agendas.arts' predefines an agenda that can be
applied directly for most users.
)--"),
      OUTPUT("iy"),
      INPUT("f_grid", "rtp_pos", "rtp_los")));

  agenda_data.push_back(AgRecord(
      NAME("iy_surface_agenda"),
      DESCRIPTION(R"--(Upwelling radiation from the surface.

The task of the agenda is to determine the upwelling intensity from
the surface, for given point and direction.

The standard choice should be to make use of *surface_rtprop_agenda*
through the WSM *iySurfaceRtpropAgenda*.

A function calling this agenda shall set *rtp_pos* and *rtp_los* to
the position and line-of-sight for which the upwelling radiation
shall be determined.

See also the include-file 'agendas.arts' for a predefined agenda
suitable to be used in most applications.
)--"),
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
            "atm_field",
            "rtp_pos",
            "rtp_los",
            "rte_pos2",
            "surface_field",
            "dsurface_names")));

  agenda_data.push_back(AgRecord(
      NAME("jacobian_agenda"),
      DESCRIPTION(R"--(Pure numerical Jacobian calculations.

Parts of the Jacobian matrix can be determined by (semi-)analytical
expressions, while other parts are calculated in a pure numerical
manner (by perturbations). This agenda describes the calculations to
be performed in the later case.

This agenda is normally not set directly by the user, but is created
by calling the the jacobianAdd set of methods.
)--"),
      OUTPUT("jacobian"),
      INPUT("jacobian", "mblock_index", "iyb", "yb")));

  agenda_data.push_back(AgRecord(
      NAME("main_agenda"),
      DESCRIPTION(R"--(The agenda corresponding to the entire controlfile. This is
executed when ARTS is run.
)--"),
      OUTPUT(),
      INPUT()));

  agenda_data.push_back(AgRecord(
      NAME("met_profile_calc_agenda"),
      DESCRIPTION(R"--(This agenda is used for metoffice profile calculations.

This agenda is called inside the method ``ybatchMetProfiles`` which is
used to make a batch calculation for the metoffice profiles.   
See the documentation of ``ybatchMetProfiles`` for more information.

This agenda can be, for example, set up like this:

1. ``AtmFieldsCalc``
2. *abs_lookupAdapt*
3. *DoitInit*
4. *DoitGetIncoming*
5. ``cloudbox_fieldSetClearsky``
6. *DoitCalc*
7. *yCalc*
)--"),
      OUTPUT("y"),
      INPUT("atm_field",
            "pnd_field_raw",
            "sensor_los",
            "cloudbox_on",
            "cloudbox_limits",
            "surface_field")));

  agenda_data.push_back(AgRecord(
      NAME("pha_mat_spt_agenda"),
      DESCRIPTION(R"--(Calculates the phase matrix for individual scattering elements.

Different options are possible for the usage of this agenda: 
*pha_mat_sptFromData* or *pha_mat_sptFromDataDOITOpt*.
)--"),
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
      DESCRIPTION(R"--(Returns particle number density data for each scattering species.

This variable is used when mapping data in ``particle_bulkprop_field``
to *pnd_field*. The variable is also necessary when calculating
scattering species weighting functions.

Note that content of this agenda array, *scat_species* and
*pnd_agenda_array_input_names* must be consistent.
)--"),
      OUTPUT("pnd_data", "dpnd_data_dx"),
      INPUT("agenda_array_index",
            "pnd_agenda_input_t",
            "pnd_agenda_input",
            "pnd_agenda_input_names",
            "dpnd_data_dx_names")));

  agenda_data.push_back(AgRecord(
      NAME("ppath_agenda"),
      DESCRIPTION(R"--(Calculation of complete propagation paths.

In contrast to *ppath_step_agenda* that controls the ray tracing
inside each grid box, this agenda determines how complete paths are
determined. The standard choice is to do this in a step-by-step
manner using *ppath_step_agenda*, with this agenda set to call
``ppathStepByStep``.

The WSV *rte_los* is both input and output as in some cases it is
determined as part of the propagation path calculations (such as
radio link calculations).
)--"),
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
      DESCRIPTION(R"--(Calculation of a propagation path step.

A propagation path step is defined as the path between some point 
to a crossing with either the pressure, latitude or longitude grid,
and this agenda performs the calculations to determine such a 
partial propagation path. The starting point is normally a grid 
crossing point, but can also be an arbitrary point inside the 
atmosphere, such as the sensor position. Only points inside the 
model atmosphere are handled.

The communication between this agenda and the calling method is 
handled by *ppath_step*. That variable is used both as input and 
output to *ppath_step_agenda*. The agenda gets back *ppath_step* 
as returned to the calling method and the last path point hold by 
the structure is accordingly the starting point for the new 
calculations. If a total propagation path shall be determined, this
agenda is called repeatedly until the starting point of the 
propagation path is found and *ppath_step* will hold all path 
steps that together make up *ppath*. The starting point is included
in the returned structure. 

The path is determined by starting at the end point and moving 
backwards to the starting point. The calculations are initiated by 
filling *ppath_step* with the practical end point of the path. 
This is either the position of the sensor (true or hypothetical), 
or some point at the top of the atmosphere (determined by
geometrical calculations starting at the sensor). This 
initialisation is not handled by *ppath_step_agenda* (but by 
the internal function ppath_start_stepping). 

The *ppath_step_agenda* put in points along the propagation path 
at all crossings with the grids, tangent points and points of 
surface reflection. It is also allowed to make agendas that put in 
additional points to fulfil some criterion, such as a maximum 
distance along the path between the points. Accordingly, the 
number of new points of each step can exceed one.

The include file 'agendas.arts' defines some agendas that can be
used here.
)--"),
      OUTPUT("ppath_step"),
      INPUT("ppath_step", "ppath_lmax", "ppath_lraytrace", "f_grid")));

  agenda_data.push_back(
      AgRecord(NAME("ppvar_rtprop_agenda"),
               DESCRIPTION(R"--(Setup propagation path variables for RTE.
)--"),
               OUTPUT("ppvar_propmat", "ppvar_dpropmat", "ppvar_src",
                      "ppvar_dsrc", "ppvar_tramat", "ppvar_dtramat",
                      "ppvar_distance", "ppvar_ddistance", "ppvar_cumtramat"),
               INPUT("ppath", "ppvar_atm", "ppvar_f", "jacobian_do")));

  agenda_data.push_back(
      AgRecord(NAME("rte_background_agenda"),
               DESCRIPTION(R"--(Compute the radiative transfer equation through 
the propagation path.
)--"),
               OUTPUT("background_rad", "diy_dx"),
               INPUT("ppath", "atm_field", "f_grid", "iy_transmittance",
                     "background_transmittance", "jacobian_do")));

  agenda_data.push_back(AgRecord(
      NAME("refr_index_air_agenda"),
      DESCRIPTION(R"--(Calculation of the refractive index of air.

This agenda should calculate the summed refractive index for all
relevant atmospheric constituents, with respect to both phase and
group velocity.
)--"),
      OUTPUT("refr_index_air", "refr_index_air_group"),
      INPUT("rtp_pressure", "rtp_temperature", "rtp_vmr", "f_grid")));
  
  agenda_data.push_back(AgRecord(
      NAME("refr_index_air_ZZZ_agenda"),
      DESCRIPTION(R"--(Calculation of the refractive index of air.

This agenda shall the refractive index for given atmospheric position.
)--"),
      OUTPUT("refr_index_air", "refr_index_air_group"),
      INPUT("rtp_pos")));

  agenda_data.push_back(AgRecord(
      NAME("sensor_response_agenda"),
      DESCRIPTION(R"--(This agenda shall provide *sensor_response* and associated variables.

So far only required when doing inversions involving some sensor variables.
)--"),
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
      DESCRIPTION(R"--(Calculates single scattering properties for individual scattering
elements from the amplitude matrix.

This agenda sets up the methods, which should be used to calculate 
the single scattering properties, i.e. the extinction matrix and the 
absorbtion vector. 

Normally you  use:
*opt_prop_sptFromMonoData*
)--"),
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
      DESCRIPTION(R"--(Provides radiative properties of the surface. 

Provides surface emission and surface reflection coefficient matrix
(see user guide for closer definitions of the respective variables
*surface_emission*, *surface_los*, and *surface_rmatrix*) according
to the characteristics of the surface specified by the methods called
within the agenda. Typical methods include *surfaceBlackbody*,
*surfaceFlatScalarReflectivity*, *surfaceFlatReflectivity*,
*surfaceFlatRefractiveIndex*, and *surfaceLambertianSimple*.
)--"),
      OUTPUT("surface_point",
             "surface_emission",
             "surface_los",
             "surface_rmatrix"),
      INPUT("f_grid", "rtp_pos", "rtp_los")));

  agenda_data.push_back(AgRecord(
      NAME("surface_rtprop_agenda_array"),
      DESCRIPTION(R"--(Description of surface radiative properties, for each surface type.

Each of these agendas shall treat the radiative properties of a
surface type. The task of these agendas is equivalent to that of
*surface_rtprop_agenda*.

The order of the agendas shall match the coding used in
*surface_type_mask*.
)--"),
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
               DESCRIPTION(R"--(Dummy agenda for testing purposes.
)--"),
               OUTPUT(),
               INPUT()));

  agenda_data.push_back(
      AgRecord(NAME("test_agenda_array"),
               DESCRIPTION(R"--(Agenda array for TestArrayOfAgenda test case.
)--"),
               OUTPUT(),
               INPUT("agenda_array_index", "iy_unit")));

  agenda_data.push_back(AgRecord(
      NAME("water_p_eq_agenda"),
      DESCRIPTION(R"--(Calculation of the saturation pressure of water.
)--"),
      OUTPUT("water_p_eq_field"),
      INPUT("atm_field")));

  agenda_data.push_back(AgRecord(
      NAME("ybatch_calc_agenda"),
      DESCRIPTION(R"--(Calculations to perform for each batch case.

Must produce a new spectrum vector (*y*) and Jacobi matrix (*jacobian*).
See further *ybatchCalc*.
)--"),
      OUTPUT("y", "y_aux", "jacobian"),
      INPUT("ybatch_index")));

  std::sort(agenda_data.begin(), agenda_data.end(), [](auto& a, auto& b) {
    return a.Name() < b.Name();
  });
}
