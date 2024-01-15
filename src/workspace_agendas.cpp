#include "workspace_agendas.h"

std::unordered_map<std::string, WorkspaceAgendaInternalRecord>
internal_workspace_agendas() {
  std::unordered_map<std::string, WorkspaceAgendaInternalRecord> wsa_data;

  wsa_data["propmat_clearsky_agenda"] = {
      .desc = R"--(Calculate the absorption coefficient matrix.

This agenda calculates the absorption coefficient matrix for all
absorption species as a function of the given atmospheric state for
one point in the atmosphere. The result is returned in
*propmat_clearsky*. The atmospheric state has to be specified by
*rtp_pressure*, *rtp_temperature*, ``rtp_mag``, and *rtp_vmr*.

The methods inside this agenda may require a lot of additional
input variables, such as *abs_species*, etc.

The include file 'agendas.arts' predefines some possible agendas
that can be used here.
)--",
      .output = {"propmat_clearsky",
                 "nlte_source",
                 "dpropmat_clearsky_dx",
                 "dnlte_source_dx"},
      .input = {"jacobian_targets",
                "select_abs_species",
                "f_grid",
                "rtp_los",
                "atm_point"}};

  wsa_data["doit_conv_test_agenda"] = {
      .desc = R"--(Compute the convergence test.

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
)--",
      .output = {"doit_conv_flag", "doit_iteration_counter"},
      .input = {"doit_conv_flag",
                "doit_iteration_counter",
                "cloudbox_field_mono",
                "cloudbox_field_mono_old"}};

  wsa_data["doit_mono_agenda"] = {
      .desc = R"--(Performs monochromatic DOIT calculation.

This agenda includes for example the following methods:

1. ``DoitScatteringDataPrepare`` 
2. ``cloudbox_field_monoIterate``

The result of the agenda is the radiation field inside the 
cloudbox and on the cloudbox boundary, which can be used 
as radiative background for a clearsky radiative transfer 
calculation.

See the Arts online documentation
for more information about the methods.
)--",
      .output = {"cloudbox_field_mono"},
      .input = {"cloudbox_field_mono", "f_grid", "f_index"}};

  wsa_data["doit_scat_field_agenda"] = {
      .desc = R"--(Calculation of the scattering integral field (DOIT). 

This agenda is called repeatedly in each DOIT iteration.
The following methods can be used for calculating the 
scattering integral field: 

-  ``doit_scat_fieldCalc``: This method calculates the scattering 
   integral field by using the angular grids *za_grid* 
   and *aa_grid*, which are also used in the update of the 
   radiation field (*doit_rte_agenda*).

-  ``doit_scat_fieldCalcLimb``: This method calculates the scattering 
   integral field.  The difference to the previous method is that 
   the data is interpolated on equidistant angular grids. 
   Especially for limb, where a very fine zenith angle grid 
   resolution is required for the RT transfer part, this method 
   is much faster than ``doit_scat_fieldCalc``.
)--",
      .output = {"doit_scat_field"},
      .input = {"doit_scat_field", "cloudbox_field_mono"}};

  wsa_data["doit_rte_agenda"] = {
      .desc = R"--(Radiative transfer calculations in cloudbox.

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
)--",
      .output = {"cloudbox_field_mono"},
      .input = {"cloudbox_field_mono", "doit_scat_field"}};

  wsa_data["gas_scattering_agenda"] = {
      .desc =
          R"--(Calculation of the gas scattering extinction and phase matrix.

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
)--",
      .output = {"gas_scattering_coef",
                 "gas_scattering_mat",
                 "gas_scattering_fct_legendre"},
      .input = {"f_grid",
                "atm_point",
                "gas_scattering_los_in",
                "gas_scattering_los_out",
                "gas_scattering_output_type"}};

  wsa_data["g0_agenda"] = {
      .desc = R"--(Calculation of the gravity at zero altitude.

Returns *g0* for given geographical position.
)--",
      .output = {"g0"},
      .input = {"lat", "lon"}};

  wsa_data["pha_mat_spt_agenda"] = {
      .desc =
          R"--(Calculates the phase matrix for individual scattering elements.

Different options are possible for the usage of this agenda: 
*pha_mat_sptFromData* or *pha_mat_sptFromDataDOITOpt*.
)--",
      .output = {"pha_mat_spt"},
      .input = {"pha_mat_spt",
                "za_index",
                "scat_lat_index",
                "scat_lon_index",
                "scat_p_index",
                "aa_index",
                "rtp_temperature"}};

  wsa_data["pnd_agenda_array"] = {
      .desc =
          R"--(Returns particle number density data for each scattering species.

This variable is used when mapping data in ``particle_bulkprop_field``
to *pnd_field*. The variable is also necessary when calculating
scattering species weighting functions.

Note that content of this agenda array, *scat_species* and
*pnd_agenda_array_input_names* must be consistent.
)--",
      .output = {"pnd_data", "dpnd_data_dx"},
      .input = {"agenda_array_index",
                "pnd_agenda_input_t",
                "pnd_agenda_input",
                "pnd_agenda_input_names",
                "dpnd_data_dx_names"},
      .array = true};

  wsa_data["refr_index_air_agenda"] = {
      .desc = R"--(Calculation of the refractive index of air.

This agenda should calculate the summed refractive index for all
relevant atmospheric constituents, with respect to both phase and
group velocity.
)--",
      .output = {"refr_index_air", "refr_index_air_group"},
      .input = {"rtp_pressure", "rtp_temperature", "rtp_vmr", "f_grid"}};

  wsa_data["refr_index_air_ZZZ_agenda"] = {
      .desc = R"--(Calculation of the refractive index of air.

This agenda shall the refractive index for given atmospheric position.
)--",
      .output = {"refr_index_air", "refr_index_air_group"},
      .input = {"rtp_pos"}};

  wsa_data["sensor_response_agenda"] = {
      .desc =
          R"--(This agenda shall provide *sensor_response* and associated variables.

So far only required when doing inversions involving some sensor variables.
)--",
      .output = {"sensor_response",
                 "sensor_response_f",
                 "sensor_response_f_grid",
                 "sensor_response_pol",
                 "sensor_response_pol_grid",
                 "sensor_response_dlos",
                 "sensor_response_dlos_grid",
                 "mblock_dlos"},
      .input = {"f_backend"}};

  wsa_data["spt_calc_agenda"] = {
      .desc =
          R"--(Calculates single scattering properties for individual scattering
elements from the amplitude matrix.

This agenda sets up the methods, which should be used to calculate 
the single scattering properties, i.e. the extinction matrix and the 
absorbtion vector. 

Normally you  use:
*opt_prop_sptFromMonoData*
)--",
      .output = {"ext_mat_spt", "abs_vec_spt"},
      .input = {"ext_mat_spt",
                "abs_vec_spt",
                "scat_p_index",
                "scat_lat_index",
                "scat_lon_index",
                "rtp_temperature",
                "za_index",
                "aa_index"}};

  wsa_data["surface_rtprop_agenda"] = {
      .desc = R"--(Provides radiative properties of the surface. 

Provides surface emission and surface reflection coefficient matrix
(see user guide for closer definitions of the respective variables
*surface_emission*, *surface_los*, and *surface_rmatrix*) according
to the characteristics of the surface specified by the methods called
within the agenda. Typical methods include *surfaceBlackbody*,
*surfaceFlatScalarReflectivity*, *surfaceFlatReflectivity*,
*surfaceFlatRefractiveIndex*, and *surfaceLambertianSimple*.
)--",
      .output = {"surface_point",
                 "surface_emission",
                 "surface_los",
                 "surface_rmatrix"},
      .input = {"f_grid", "rtp_pos", "rtp_los"}};

  wsa_data["surface_rtprop_agenda_array"] = {
      .desc =
          R"--(Description of surface radiative properties, for each surface type.

Each of these agendas shall treat the radiative properties of a
surface type. The task of these agendas is equivalent to that of
*surface_rtprop_agenda*.

The order of the agendas shall match the coding used in
*surface_type_mask*.
)--",
      .output = {"surface_skin_t",
                 "surface_emission",
                 "surface_los",
                 "surface_rmatrix"},
      .input = {"agenda_array_index", "f_grid", "rtp_pos", "rtp_los"},
      .array = true};

  wsa_data["water_p_eq_agenda"] = {
      .desc = R"--(Calculation of the saturation pressure of water.
)--",
      .output = {"water_p_eq_field"},
      .input = {"atm_field"}};

  wsa_data["space_radiation_agenda"] = {
      .desc = R"--(Radiation as seen of space.

This agenda calculates the radiation as seen of space. The
intent is to provide a background radiation from space that
is input to the atmospheric radiative transfer calculations.

The input position and line-of-sight are as if you were looking
at space.

The output must be sized as:

- *spectral_radiance_background* : (*f_grid*)
- *spectral_radiance_background_jacobian* : (*jacobian_targets*, *f_grid*)
)--",
      .output = {"spectral_radiance_background",
                 "spectral_radiance_background_jacobian"},
      .input = {"f_grid", "jacobian_targets", "rtp_pos", "rtp_los"}};

  wsa_data["spectral_radiance_background_space_agenda"] = {
      .desc = R"--(Spectral radiance as seen of space.

This agenda calculates the spectral radiance as seen of space. The
intent is to provide a background spectral radiance from space that
is input to the atmospheric radiative transfer calculations.

The input path point should be as if it is looking at space.

The output must be sized as:

- *spectral_radiance_background* : (*f_grid*)
- *spectral_radiance_background_jacobian* : (*jacobian_targets*, *f_grid*)
)--",
      .output = {"spectral_radiance_background",
                 "spectral_radiance_background_jacobian"},
      .input = {"f_grid", "jacobian_targets", "path_point"}};

  wsa_data["surface_radiation_agenda"] = {
      .desc = R"--(Radiation as seen of the surface.

This agenda calculates the radiation as seen of the surface.
The intent is to provide a background radiation from the
surface that is input to the atmospheric radiative transfer
calculations.

The input path point should be as if it is looking
at the surface.

The output must be sized as:

- *spectral_radiance_background* : (*f_grid*)
- *spectral_radiance_background_jacobian* : (*jacobian_targets*, *f_grid*)
)--",
      .output = {"spectral_radiance_background",
                 "spectral_radiance_background_jacobian"},
      .input = {"f_grid", "jacobian_targets", "rtp_pos", "rtp_los"}};

  wsa_data["spectral_radiance_background_surface_agenda"] = {
      .desc = R"--(Spectral radiance as seen of the surface.

This agenda calculates the spectral radiance as seen of the surface.
The intent is to provide a background spectral radiance from the
surface that is input to the atmospheric radiative transfer
calculations.

The input path point should be as if it is looking at the surface.

The output must be sized as:

- *spectral_radiance_background* : (*f_grid*)
- *spectral_radiance_background_jacobian* : (*jacobian_targets*, *f_grid*)
)--",
      .output = {"spectral_radiance_background",
                 "spectral_radiance_background_jacobian"},
      .input = {"f_grid", "jacobian_targets", "path_point"}};

  wsa_data["stop_distance_radiation_agenda"] = {
      .desc = R"--(Radiation as seen from stopping the ppath calculations.

The output must be sized as:

- *spectral_radiance_background* : (*f_grid*)
- *spectral_radiance_background_jacobian* : (*jacobian_targets*, *f_grid*)
)--",
      .output = {"spectral_radiance_background",
                 "spectral_radiance_background_jacobian"},
      .input = {"f_grid", "jacobian_targets", "rtp_pos", "rtp_los"}};

  return wsa_data;
}
