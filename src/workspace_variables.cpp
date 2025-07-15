#include "workspace_variables.h"

#include "workspace_agendas.h"

void agendas(std::unordered_map<std::string, WorkspaceVariableInternalRecord>&
                 wsv_data) {
  for (auto& [name, ag] : internal_workspace_agendas()) {
    if (ag.enum_default.empty()) {
      wsv_data[name] = {
          .desc = ag.desc,
          .type = "Agenda",
      };
    } else {
      wsv_data[name] = {
          .desc          = ag.desc,
          .type          = "Agenda",
          .default_value = "get_" + name + "(\"" + ag.enum_default + "\"sv)",
      };
    }
  }
}

std::unordered_map<std::string, WorkspaceVariableInternalRecord>
internal_workspace_variables_creator() {
  std::unordered_map<std::string, WorkspaceVariableInternalRecord> wsv_data;

  wsv_data["absorption_bands"] = {
      .desc = R"--(Bands of absorption lines for line-by-line (LBL) calculations.

See methods that consume this variable for more details on its content.
)--",
      .type = "AbsorptionBands",
  };

  wsv_data["absorption_lookup_table"] = {
      .desc =
          R"--(Absorption lookup table for scalar gas absorption coefficients.

Precomputing this table replaces the need for the calculation of scalar gas line-by-line
absorption.
)--",
      .type = "AbsorptionLookupTables",
  };

  wsv_data["absorption_cia_data"] = {
      .desc = R"--(HITRAN Collision-Induced Absorption (CIA) Data.

This variable holds HITRAN CIA data (binary absorption
cross-sections). The data itself is described in: Richard, C. et al.
(2012), New section of the HITRAN database: Collision-induced
absorption (CIA), J. Quant. Spectrosc. Radiat. Transfer, 113,
1276-1285, doi:10.1016/j.jqsrt.2011.11.004.

The binary absorption cross-sections have to be multiplied with the
densities of both molecules to get a scalar absorption coefficient.

Dimensions:

The length of this array should be equal to the number of pairs
of molecules that have CIA data available.
Some methods that split the data might not work as intended otherwise.
)--",
      .type = "ArrayOfCIARecord",
  };

  wsv_data["absorption_species"] = {
      .desc = R"--(Tag groups for gas absorption.

This allows the user to set up groups of species tags that are used
to load the correct data.

It is only used to let data-reading methods know which
species they should read from the available input files.
)--",
      .type = "ArrayOfArrayOfSpeciesTag",
  };

  wsv_data["altitude"] = {
      .desc =
          R"--(A single altitude in the atmosphere.

Unit: m
)--",
      .type          = "Numeric",
      .default_value = "0.0",
  };

  wsv_data["altitude_grid"] = {
      .desc =
          R"--(An ascending list of *altitude*.  Often related to a field or a profile.

Unit: m

.. note::
    There is no global grid system in ARTS, so beware of the local
    nature of all grids.
)--",
      .type = "AscendingGrid",
  };

  wsv_data["atmospheric_field"] = {
      .desc =
          R"--(An atmospheric field in ARTS.

The atmospheric field defines the altitude of the top-of-the-atmosphere,
as well as the variables that are required for the radiative transfer
calculations along a path through the atmosphere.  The field can be 
accessed at any altitude, latitude, longitude path that is within the
atmosphere to access the relevant atmospheric point data (*atmospheric_point*).

Note that constraints on the various field parameters may be imposed by
extrapolation limitations on the field parameter itself, causing some or
large swaths of the atmosphere to be inaccessible.

The atmospheric field may, but does not have to, consist of the following:

-  Temperature                     - Atmospheric temperatures in Kelvin
-  Pressure                        - Atmospheric pressure in Pascal
-  Wind                            - Atmospheric wind field in meters per second
-  Magnetic Field                  - Magnetic field in Tesla
-  Species content                 - Usually the volume-mixing ratio of various species, with some exceptions. See *SpeciesEnum* for more details.
-  Isotopologue ratios             - The isotopologue ratios of various species.  See *SpeciesIsotope* for more details.
-  Non-local thermodynamics ratios - Unitless [pure-style] OR Kelvin [vibrational-style] ratios replacing the Boltzman distribution used in the LTE calculations.
-  Scattering species content      - See user guide for more information.  This is custom data to aid scattering calculations.
)--",
      .type = "AtmField",
  };

  wsv_data["atmospheric_point"] = {
      .desc =
          R"--(An atmospheric point in ARTS.

The atmospheric point consists of all the relevant atmospheric field data
at a discrete point in the atmosphere.  It is often extracted from an *AtmField*
at a single altitude-latitude-longitude but may of course be generated manually.

See *atmospheric_field* for the data that may be available in the atmospheric point.
)--",
      .type = "AtmPoint",
  };

  wsv_data["atmospheric_profile"] = {
      .desc =
          R"--(An atmospheric profile in ARTS.

The atmospheric profile consists of all the relevant atmospheric field data
at a discrete profile in the atmosphere.  It is often extracted from an *AtmField*
at a single latitude-longitude coordinate but may of course be generated manually.

See *atmospheric_field* for the data that may be available in the atmospheric point.

The size of the profile is the same as *altitude_grid*.
)--",
      .type = "ArrayOfAtmPoint",
  };

  wsv_data["propagation_matrix_jacobian"] = {
      .desc =
          R"--(Partial derivative of the *propagation_matrix* with regards to *jacobian_targets*.

The units depend on what is set in *jacobian_targets* [1 / m / jacobian target's unit].
)--",
      .type = "PropmatMatrix",
  };

  wsv_data["propagation_matrix_source_vector_nonlte_jacobian"] = {
      .desc =
          R"--(Partial derivative of the *propagation_matrix_source_vector_nonlte* with regards to *jacobian_targets*.

The units are *spectral_radiance_jacobian* per meter.
)--",
      .type = "StokvecMatrix",
  };

  wsv_data["ecs_data"] = {
      .desc          = R"--(Error corrected sudden data

Dimensions: [num Isotopologues] [num Species]

Used in line-by-line calculations requiring ECS data.
)--",
      .type          = "LinemixingEcsData",
      .default_value = " ",
  };

  wsv_data["frequency_grid_wind_shift_jacobian"] = {
      .desc =
          R"--(The frequency grid wind shift Jacobian.

Used because all methods inside *propagation_matrix_agenda* work on
the frequency grid, not on the actual wind speed for the sake of
wind shift Jacobian calculations.

The order is ``[df_du, df_dv, df_fw]``
)--",
      .type          = "Vector3",
      .default_value = "0.0, 0.0, 0.0",
  };

  wsv_data["ray_path_frequency_grid_wind_shift_jacobian"] = {
      .desc =
          R"--(A list of *frequency_grid_wind_shift_jacobian* for a ray path.
)--",
      .type = "ArrayOfVector3",
  };

  wsv_data["absorption_xsec_fit_data"] = {
      .desc = R"--(Fitting model coefficients for cross section species.

Dimensions: [ n_species ]

XsecRecord:

- species: Name of species
- version: Fit model version
- fitcoeffs:

  -  Fit model coefficients as an *ArrayOfGriddedField2*
  -  Dimensions: [ n_bands ]

    -   GriddedField2: [ n_band_frequencies, n_coeffs ]

      -    The fit model:

        -     z = p00 + p10*x + p01*y + p20*x^2
        -     z = Xsec [m^2]
        -     x = T / T0
        -     y = P / P0
        -     T0 = 1 [K]
        -     P0 = 1 [Pa]
        -     fitcoeffs(:, 0)           p00  [m^2]
        -     fitcoeffs(:, 1)           p10  [m^2]
        -     fitcoeffs(:, 2)           p01  [m^2]
        -     fitcoeffs(:, 3)           p20  [m^2]

- fitminpressures:

  -  Minimum pressure available in source xsec data to generate the fit coefficients.
  -  Dimensions: [ n_bands ]

- fitmaxpressures:

  -  Maximum pressure available in source xsec data to generate the fit coefficients.
  -  Dimensions: [ n_bands ]

- fitmintemperatures:

  -  Minimum temperature available in source xsec data to generate the fit coefficients.
  -  Dimensions: [ n_bands ]

- fitmintemperatures:

  -  Maximum temperature available in source xsec data to generate the fit coefficients.
  -  Dimensions: [ n_bands ]

fitminpressures, fitmaxpressures, fitmintemperatures and fitmaxtemperatures
are not used to apply the model and solely serve for informational purposes.
)--",
      .type = "ArrayOfXsecRecord",
  };

  wsv_data["propagation_matrix_source_vector_nonlte"] = {
      .desc =
          R"--(The part of the source vector that is due to non-LTE.

This is closely related to *propagation_matrix*.

Gven the level source term:

.. math:: \vec{J} = \mathbf{K}^{-1} \left(\vec{\alpha}B + \vec{J}_n + \cdots\right),

this variable holds :math:`\vec{J}_n`.  Here, :math:`\vec{\alpha}` is the first
column of :math:`\mathbf{K}`, which is from the *propagation_matrix* variable.
:math:`B` is the Planck function.  The ellipsis denotes other terms that can
come from more sources, such as scattering and/or transmitting equipment.

The unit is in *spectral_radiance* per meter.
)--",
      .type = "StokvecVector",
  };

  wsv_data["ray_path_atmospheric_point"] = {
      .desc = R"--(Atmospheric points along the propagation path.

See *atmospheric_point* for information about atmospheric points

Dimension: [ ppath.np ]

Usage: Output of radiative transfer methods.
)--",
      .type = "ArrayOfAtmPoint",
  };

  wsv_data["ray_path_frequency_grid"] = {
      .desc = R"--(All *frequency_grid* along the propagation path.
)--",
      .type = "ArrayOfAscendingGrid",
  };

  wsv_data["absorption_predefined_model_data"] = {
      .desc =
          R"--(This contains predefined model data.

Can currently only contain data for new MT CKD models of water.
)--",
      .type = "PredefinedModelData",
  };

  wsv_data["propagation_matrix"] = {
      .desc =
          R"--(This contains the fully polarized propagation matrix for the current path point.

The propagation matrix can be used to computed the transmission matrix as:

.. math:: \mathbf{T} = \exp\left(-\mathbf{K} r\right),

where :math:`\mathbf{K}` is the propagation matrix, and :math:`r` is some distance
over which it is considered constant.

The unit is [1 / m].

Dimension: *frequency_grid*.
)--",
      .type = "PropmatVector",
  };

  wsv_data["propagation_matrix_scattering"] = {
      .desc =
          R"--(This contains the propagation matrix for scattering for the current path point.

This needs to be used when scattering into the line of sight is considered. And it needs then to
also be added to the *propagation_matrix*, which you should see for more information.

The unit is [1 / m].

Dimension: *frequency_grid*.
)--",
      .type = "PropmatVector",
  };

  wsv_data["select_species"] = {
      .desc          = R"--(Species selection.
      
When Bath is selected, all species are used.  Otherwise, this variable should control so that only the selected species is used.
)--",
      .type          = "SpeciesEnum",
      .default_value = "SpeciesEnum::Bath",
  };

  wsv_data["select_species_list"] = {.desc = R"--(Species selection.
)--",
                                     .type = "ArrayOfSpeciesEnum"};

  wsv_data["spectral_radiance_background_jacobian"] = {
      .desc = R"--(Spectral radiance derivative from the background

Shape: NJAC x NFREQ
)--",
      .type = "StokvecMatrix",
  };

  wsv_data["spectral_radiance_background"] = {
      .desc = R"--(Spectral radiance from the background

Shape: NFREQ
)--",
      .type = "StokvecVector",
  };

  wsv_data["transmission_matrix_background"] = {
      .desc = R"--(Transmittance from the background
)--",
      .type = "MuelmatVector",
  };

  wsv_data["propagation_matrix_scattering"] = {
      .desc =
          R"--(The propgation matrix of totally random orientation particles at a single point along a path using spectral representation
)--",
      .type = "PropmatVector",
  };

  wsv_data["ray_path_propagation_matrix_scattering"] = {
      .desc =
          R"--(The propgation matrix of totally random orientation particles along the propagation path using spectral representation
)--",
      .type = "ArrayOfPropmatVector",
  };

  wsv_data["absorption_vector_scattering"] = {
      .desc =
          R"--(The absorption vector of totally random orientation particles at a single point along a path using spectral representation
)--",
      .type = "StokvecVector",
  };

  wsv_data["ray_path_absorption_vector_scattering"] = {
      .desc =
          R"--(The absorption vector of totally random orientation particles along the propagation path using spectral representation
)--",
      .type = "ArrayOfStokvecVector",
  };

  wsv_data["phase_matrix_scattering_spectral"] = {
      .desc =
          R"--(The spectral phase matrix of totally random orientation particles at a single point along a path using spectral representation
)--",
      .type = "SpecmatMatrix",
  };

  wsv_data["ray_path_phase_matrix_scattering_spectral"] = {
      .desc =
          R"--(The spectral phase matrix of totally random orientation particles along the propagation path using spectral representation
)--",
      .type = "ArrayOfSpecmatMatrix",
  };

  wsv_data["scattering_species"] = {
      .desc = R"--(The scattering species
)--",
      .type = "ArrayOfScatteringSpecies",
  };

  wsv_data["ray_path_spectral_radiance_scattering"] = {
      .desc = R"--(Spectral radiance scattered into the propagation path
)--",
      .type = "ArrayOfStokvecVector",
  };

  wsv_data["ray_path_spectral_radiance_jacobian"] = {
      .desc = R"--(Spectral radiance derivative along the propagation path
)--",
      .type = "ArrayOfStokvecMatrix",
  };

  wsv_data["ray_path_propagation_matrix"] = {
      .desc = R"--(Propagation matrices along the propagation path
)--",
      .type = "ArrayOfPropmatVector",
  };

  wsv_data["ray_path_propagation_matrix_scattering"] = {
      .desc =
          R"--(Propagation matrices along the propagation path for scattering
)--",
      .type = "ArrayOfPropmatVector",
  };

  wsv_data["ray_path_propagation_matrix_jacobian"] = {
      .desc = R"--(Propagation derivative matrices along the propagation path
)--",
      .type = "ArrayOfPropmatMatrix",
  };

  wsv_data["ray_path_propagation_matrix_source_vector_nonlte"] = {
      .desc = R"--(Additional non-LTE along the propagation path
)--",
      .type = "ArrayOfStokvecVector",
  };

  wsv_data["ray_path_propagation_matrix_source_vector_nonlte_jacobian"] = {
      .desc = R"--(Additional non-LTE derivative along the propagation path
)--",
      .type = "ArrayOfStokvecMatrix",
  };

  wsv_data["ray_path_spectral_radiance_source"] = {
      .desc = R"--(Source vectors along the propagation path
)--",
      .type = "ArrayOfStokvecVector",
  };

  wsv_data["ray_path_spectral_radiance_source_jacobian"] = {
      .desc = R"--(Source derivative vectors along the propagation path
)--",
      .type = "ArrayOfStokvecMatrix",
  };

  wsv_data["ray_path_transmission_matrix"] = {
      .desc = R"--(Transmission matrices along the propagation path.

The outer dimension is the number of layers.

The inner dimension is the number of frequency points.

The order of the elements is such that index zero is closest to the obeserver.
)--",
      .type = "ArrayOfMuelmatVector",
  };

  wsv_data["ray_path_transmission_matrix_cumulative"] = {
      .desc = R"--(Cumulative transmission matrices along the propagation path
)--",
      .type = "ArrayOfMuelmatVector",
  };

  wsv_data["ray_path_transmission_matrix_jacobian"] = {
      .desc = R"--(Transmission derivative matrices along the propagation path.

The outer dimension is the number of layers.

The inner dimensions are the number of level derivatives,
the number of jacbian targets, and the number of frequency points.
The required number of level derivatives is determined by the appropriate
method (a common value is 2, for the 2 levels surrounding a layer).

The order of the elements is such that index zero is closest to the obeserver.
)--",
      .type = "ArrayOfMuelmatTensor3",
  };

  wsv_data["surface_field"] = {
      .desc          = R"--(The surface field describes the surface properties.

This describes the global surface values, such as elevation and 
temperature but also entirerly abstract properties and types.
)--",
      .type          = "SurfaceField",
      .default_value = "SurfaceField()",
  };

  wsv_data["subsurface_field"] = {
      .desc = R"--(The sub0surface field describes the sub-surface properties.
)--",
      .type = "SubsurfaceField",
      .default_value = "SubsurfaceField()",
  };

  wsv_data["gravity_operator"] = {
      .desc = R"--(The gravity operator.

Usage: gravity = gravity_operator(altitude, latitude, longitude).

Parameters
----------
altitude : Numeric
    Altitude in meters.
latitude : Numeric
    Latitude in degrees.
longitude : Numeric
    Longitude in degrees.

Returns
-------
gravity : Numeric
    The gravity in m/s :math:`^2`.
)--",
      .type = "NumericTernaryOperator",
  };

  wsv_data["water_equivalent_pressure_operator"] = {
      .desc = R"--(The water equivalent pressure operator.

Usage: psat = water_equivalent_pressure_operator(temperature).

Parameters
----------
temperature : Numeric
    Temperature in Kelvin.

Returns
-------
psat : Numeric
    The water equivalent pressure in Pascal.
)--",
      .type = "NumericUnaryOperator",
  };

  wsv_data["spectral_radiance_field"] = {
      .desc = R"(The spectral radiance field.

*spectral_radiance* but for a field.

Dimensions are *altitude_grid* times *latitude_grid* times *longitude_grid* times *zenith_grid* times ``azimuth_grid`` times *frequency_grid*.
)",
      .type = "StokvecSortedGriddedField6",
  };

  wsv_data["spectral_radiance_transform_operator"] = {
      .desc = R"(The spectral radiance transform operator

This is responsible for things like converting the spectral radiance
into a different unit, e.g., from [W / m :math:`^2` sr Hz] to Kelvin.
)",
      .type = "SpectralRadianceTransformOperator",
      .default_value =
          "SpectralRadianceTransformOperator(SpectralRadianceUnitType::unit)"};

  wsv_data["spectral_radiance"] = {
      .desc = R"--(A spectral radiance vector.

This is the representation of the spectral radiances at discrete frequencies for
a discrete viewing direction.

The unit of spectral radiance is [W / m :math:`^2` sr Hz].

Note that there are conversion routines that changes this unit,
e.g., *spectral_radianceApplyUnit*.  After conversion,
the use of *spectral_radiance* in any method no marked as safe for different units,
will lead to undefined behavior with possibly bad values being computed.

The size of this variable should be the size of the local *frequency_grid*.
)--",
      .type = "StokvecVector",
  };

  wsv_data["spectral_radiance_jacobian"] = {
      .desc =
          R"--(Jacobian of *spectral_radiance* with respect to *jacobian_targets*.

The size of this variable should be the local *jacobian_targets* as rows times the
size of the local *spectral_radiance* as columns.
)--",
      .type = "StokvecMatrix",
  };

  wsv_data["jacobian_targets"] = {
      .desc = R"--(A list of targets for the Jacobian Matrix calculations.
)--",
      .type = "JacobianTargets",
      .default_value = " ",
  };

  wsv_data["ray_path_point"] = {
      .desc = R"--(A single path point.
)--",
      .type = "PropagationPathPoint",
  };

  wsv_data["frequency_grid"] = {
      .desc = R"--(A single frequency grid.

Units: Hz

.. note::
    There is no global grid system in ARTS, so beware of the local
    nature of all grids.
)--",
      .type = "AscendingGrid",
  };

  wsv_data["zenith_grid"] = {
      .desc = R"--(A single zenith angle grid.

Units: degrees

.. note::
    There is no global grid system in ARTS, so beware of the local
    nature of all grids.
)--",
      .type = "AscendingGrid",
  };

  wsv_data["ray_path"] = {
      .desc = R"--(A list path points making up a propagation path.
)--",
      .type = "ArrayOfPropagationPathPoint",
  };

  wsv_data["ray_path_field"] = {
      .desc =
          R"--(A list of *ray_path* intended to build up a field of observations.
)--",
      .type = "ArrayOfArrayOfPropagationPathPoint",
  };

  wsv_data["ray_path_observers"] = {
      .desc =
          R"--(A list path points making up the observers of a propagation path.

These can be used directly for *spectral_radiance_observer_position* and *spectral_radiance_observer_line_of_sight*
)--",
      .type = "ArrayOfPropagationPathPoint",
  };

  wsv_data["spectral_flux_profile"] = {
      .desc = R"--(An altitude profile of spectral flux.
)--",
      .type = "Matrix",
  };

  wsv_data["nlte_line_flux_profile"] = {
      .desc = R"--(A per-line flux profile.
)--",
      .type = "QuantumIdentifierVectorMap",
  };

  wsv_data["spectral_radiance_observer_position"] = {
      .desc = R"--(The position of an observer of spectral radiance.

Most likely only makes sense in combination with *spectral_radiance_observer_line_of_sight*.
)--",
      .type = "Vector3",
  };

  wsv_data["spectral_radiance_observer_line_of_sight"] = {
      .desc = R"--(The line-of-sight of the observer of spectral radiance.

Most likely only makes sense in combination with *spectral_radiance_observer_position*.
)--",
      .type = "Vector2",
  };

  wsv_data["spectral_radiance_operator"] = {
      .desc = R"--(The spectral radiance operator.

This is a class that can compute the spectral radiance
along a path for a single viewing direction and frequency.

It provides several methods to get the path of the spectral
radiance.
)--",
      .type = "SpectralRadianceOperator",
  };

  wsv_data["model_state_vector"] = {
      .desc          = R"(A state vector of the model.

This represents the :emphasis:`chosen` state of the model.
In the notation of *measurement_vector* and *OEM*,
:math:`\vec{x}` is the *model_state_vector*.

To choose the state of the model, you must setup *jacobian_targets* to
include the state parameters you want to be able to change.
)",
      .type          = "Vector",
      .default_value = " ",
  };

  wsv_data["model_state_vector_apriori"] = {
      .desc = R"(An apriori state vector of the model.

See *model_state_vector* for more details.
This is the state vector that is assumed to be the a priori state of the model.
In normal circumstances, this is the state vector that is used to
start the inversion process.  In *OEM*, this is :math:`\vec{x}_a`.
)",
      .type = "Vector",
  };

  wsv_data["measurement_vector"] = {
      .desc = R"(The measurment vector for, e.g., a sensor.

This must often be the same size as *measurement_sensor*.

The notation in ARTS, for the purpose of *OEM*, is that

.. math::
    \vec{y} = \mathbf{F}\left(\vec{x}\right) + \vec{y}_\epsilon\left(\vec{x}\right) + \epsilon

where
:math:`\mathbf{F}` is the forward model function of the physics of the simulation space,
:math:`\vec{x}` is the *model_state_vector*,
:math:`\vec{y}_\epsilon` is the *measurement_vector_error*, and
:math:`\epsilon` are any additional errors, such as random noise.

Throughout ARTS, *measurement_vector* have different contextual meanings.
These are:

1. :math:`\vec{y}` - i.e., measured data.
2. :math:`\vec{y} - \epsilon` - e.g., the best fit to measured data, *measurement_vector_fitted*. 
3. :math:`\mathbf{F}\left(\vec{x}\right)` - i.e., the physical model of the measurement.
)",
      .type = "Vector",
  };

  wsv_data["measurement_vector_error"] = {
      .desc = R"(The model measurment vector error for, e.g., a sensor.

This must often be the same size as *measurement_sensor*.

See *measurement_vector* for more details.
In that notation, this is :math:`\vec{y}_\epsilon`.
)",
      .type = "Vector",
  };

  wsv_data["measurement_vector_fitted"] = {
      .desc          = R"(As *measurement_vector*, but fitted to the model.

This must often be the same size as *measurement_sensor*.

See *measurement_vector* for more details.
In that notation, and in the notation of *OEM*,
:math:`\vec{y}_f \approx \vec{y} - \epsilon`.
Or at least this should be the case depending on how good of a fit of :math:`\vec{x}`
has been produced and if the measurement can be understood properly.

.. tip::
    It is often useful to present :math:`\vec{y} -  \vec{y}_\epsilon`
    and :math:`\vec{y}_f - \vec{y}_\epsilon` instead of
    :math:`\vec{y}_f` and :math:`\vec{y}` directly.  This removes the
    known measurement error from both the data and the fit,
    showing the physical signal from the target rather than
    known sensor noise.
)",
      .type          = "Vector",
      .default_value = " ",
  };

  wsv_data["measurement_gain_matrix"] = {
      .desc =R"(Contribution function (or gain) matrix.

This matrix is the partial derivative of the retrieved state vector with respect to the *measurement_vector*.

Usage: Used and set by inversion methods.
)",
      .type = "Matrix",
  };

  wsv_data["measurement_averaging_kernel"] = {
      .desc =
          R"(Averaging kernel matrix.

This matrix is the partial derivative of the retrieved state vector with respect to the *measurement_vector*.

Usage: Used and set by inversion methods.
)",
      .type = "Matrix",
  };

  wsv_data["inversion_iterate_agenda_counter"] = {
      .desc          = R"(A counter for the inversion iterate agenda.
)",
      .type          = "Index",
      .default_value = "0",
  };

  wsv_data["do_jacobian"] = {
      .desc = R"(A boolean calculations related to the *measurement_jacobian* should be ignored.

This variable is limited to very few methods related to the inversion process for *OEM*.
Note that deep code of ARTS will ignore this variable, so it is not a global switch.
Instead, it is used as a switch to clear the *jacobian_targets* variable, which is used
to determine the size of the *measurement_jacobian*.  It is important to be careful
with this, as it will mess with the size of the *measurement_jacobian* and could
thus lead to runtime errors being thrown in places where unexpected sizes are encountered.
)",
      .type = "Index",
      .default_value = "1",
  };

  wsv_data["measurement_vector_error_covariance_matrix"] = {
      .desc = R"(Covariance matrix for observation uncertainties.
)",
      .type = "CovarianceMatrix",
  };

  wsv_data["model_state_covariance_matrix"] = {
      .desc = R"(Covariance matrix of a priori distribution.
)",
      .type = "CovarianceMatrix",
  };

  wsv_data["measurement_jacobian"] = {
      .desc = R"(The first order partial derivatives of the *measurement_vector*.

This variable represents the matrix

.. math::
    \mathbf{J} = \frac{\partial \vec{y}} {\partial \vec{x}},

where :math:`\vec{y}` is the *measurement_vector* and :math:`\vec{x}` is the *model_state_vector*.
The size of this variable should thus be the size of *measurement_vector* times the size of *model_state_vector*.
Please refer to those variables for more information.
)",
      .type = "Matrix",
  };

  wsv_data["measurement_jacobian_error"] = {
      .desc = R"(The partial derivatives of the *measurement_vector_error*.

This is otherwise the same as *measurement_jacobian*.  See it for more details.
)",
      .type = "Matrix",
  };

  wsv_data["measurement_sensor"] = {
      .desc = R"(A list of sensor elements.

Size is number of elements of the sensor.
)",
      .type = "ArrayOfSensorObsel",
  };

  wsv_data["sun"] = {
      .desc = R"(A sun.
)",
      .type = "Sun",
  };

  wsv_data["suns"] = {
      .desc = R"(A list of *Sun*.

Size is number of suns.
)",
      .type = "ArrayOfSun",
  };

  wsv_data["sun_path"] = {
      .desc = R"(A path to a sun if it is visible.

A related variable is *ray_path*

Size is number of path points for the sun.
)",
      .type = "ArrayOfPropagationPathPoint",
  };

  wsv_data["ray_path_suns_path"] = {
      .desc = R"(A list of paths to the suns from the ray path.

Dimensions: *ray_path* x *suns* x *sun_path*
)",
      .type = "ArrayOfArrayOfArrayOfPropagationPathPoint",
  };

  wsv_data["disort_settings"] = {
      .desc = R"(Contains the full settings of spectral Disort calculations.
)",
      .type = "DisortSettings",
  };

  wsv_data["disort_quadrature_angles"] = {
      .desc = R"(The quadrature angles for Disort.

Unit is in degrees.

Size is *disort_quadrature_dimension*.
)",
      .type = "Vector",
  };

  wsv_data["disort_quadrature_weights"] = {
      .desc = R"(The quadrature weights for Disort.

These weights are symmetric for uplooking and downlooking.

In essence, the matching *disort_quadrature_angles* for the weights can
be found as [*disort_quadrature_weights*, *disort_quadrature_weights*].

Size is *disort_quadrature_dimension* / 2
)",
      .type = "Vector",
  };

  wsv_data["disort_quadrature_dimension"] = {
      .desc = R"(The quadrature size for Disort.)",
      .type = "Index",
  };

  wsv_data["disort_fourier_mode_dimension"] = {
      .desc = R"(The number of Fourier modes for Disort.
)",
      .type = "Index",
  };

  wsv_data["disort_legendre_polynomial_dimension"] = {
      .desc = R"(The number of input Legendre polynimials for Disort.
)",
      .type = "Index",
  };

  wsv_data["latitude"] = {
      .desc =
          R"--(A single latitude.

Units: degrees
)--",
      .type          = "Numeric",
      .default_value = "0.0",
  };

  wsv_data["latitude_grid"] = {
      .desc =
          R"--(An ascending list of *latitude*.  Often related to a field or a profile.

Units: degrees

.. note::
    There is no global grid system in ARTS, so beware of the local
    nature of all grids.
)--",
      .type = "AscendingGrid",
  };

  wsv_data["longitude"] = {
      .desc =
          R"--(A single longitude.

Units: degrees
)--",
      .type          = "Numeric",
      .default_value = "0.0",
  };

  wsv_data["longitude_grid"] = {
      .desc =
          R"--(An ascending list of *longitude*.  Often related to a field or a profile.

Units: degrees

.. note::
    There is no global grid system in ARTS, so beware of the local
    nature of all grids.
)--",
      .type = "AscendingGrid",
  };

  wsv_data["legendre_degree"] = {
      .desc = R"(The degree of a Legendre polynimial.
)",
      .type = "Index",
  };

  wsv_data["disort_spectral_radiance_field"] = {
      .desc = R"(The spectral radiance field from Disort.

Size is *frequency_grid* times *ray_path* - 1 times azimuthal angles
times *disort_quadrature_dimension*.
)",
      .type = "Tensor4",
  };

  wsv_data["disort_spectral_flux_field"] = {
      .desc = R"(The spectral flux field from Disort.

Size is *frequency_grid* times 3 times *ray_path* - 1.

The inner "3" is in order: upwelling, diffuse downwelling, and direct downwelling.
)",
      .type = "Tensor3",
  };

  wsv_data["covariance_matrix_diagonal_blocks"] = {
      .desc = R"(A helper map for setting the covariance matrix.
)",
      .type = "JacobianTargetsDiagonalCovarianceMatrixMap",
  };

  agendas(wsv_data);

  return wsv_data;
}

std::string_view any(const std::string& type) {
  if (type == "Any") {
    return "T";
  }
  return type;
}

const std::unordered_map<std::string, WorkspaceVariableInternalRecord>&
internal_workspace_variables() {
  static const auto wsv_data = internal_workspace_variables_creator();
  return wsv_data;
}
