#include "workspace_variables.h"

std::unordered_map<std::string, WorkspaceVariableInternalRecord>
internal_workspace_variables() {
  std::unordered_map<std::string, WorkspaceVariableInternalRecord> wsv_data;

  wsv_data["absorption_bands"] = {
      .desc = R"--(Bands of absorption lines for LBL calculations.
)--",
      .type = "ArrayOfAbsorptionBand",
  };

  wsv_data["absorption_cia_data"] = {
      .desc = R"--(HITRAN Collision Induced Absorption (CIA) Data.

This variable holds HITRAN CIA data (binary absorption
cross-sections). The data itself is described in: Richard, C. et al.
(2012), New section of the HITRAN database: Collision-induced
absorption (CIA), J. Quant. Spectrosc. Radiat. Transfer, 113,
1276-1285, doi:10.1016/j.jqsrt.2011.11.004.

The binary absorption cross-sections have to be multiplied with the
densities of both molecules to get absorption coefficients.

Dimensions:

The outer array dimension in the ArrayOfArrayOfCIARecord is the same
as that of *absorption_species*. There will be CIA data only for those
species that contain a CIA tag, for all other species it will be
empty. The inner array dimension corresponds to the number of CIA tags
for this species (there could be for example both N2-N2 and N2-H2) in
the same species.

The CIA *absorption_species* tags are described in *absorption_speciesSet*.

Each individual CIARecord holds the complete information from one
HITRAN CIA file. For the given pair of molecules A HITRAN CIA data
file can hold several datasets (data for different temperatures but
fixed frequency range).

Units:

 - Frequencies: Hz
 - Binary absorption cross-sections: m^5*molecule^-2
)--",
      .type = "ArrayOfCIARecord",
  };

  wsv_data["absorption_lookup_table_data"] = {
      .desc = R"--(An absorption lookup table.

It holds an absorption lookup table, as well as all information that
is necessary to use the table to extract absorption. Extraction
routines are implemented as member functions. 

It has quite a complicated structure. For details see the Arts User
Guide section \"The gas absorption lookup table\" or the source code
documentation in gas_abs_lookup.h.
)--",
      .type = "GasAbsLookup",
  };

  wsv_data["absorption_species"] = {
      .desc = R"--(Tag groups for gas absorption.

This is an array of arrays of SpeciesTag tag definitions. It defines the
available tag groups for the calculation of scalar gas absorption
coefficients.  See online documentation of method *absorption_speciesSet* for
more detailed information how tag groups work and some examples.
)--",
      .type = "ArrayOfArrayOfSpeciesTag",
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
    Temperature                             - Kelvin
    Pressure                                - Pascal
    Wind                                    - Meters per second
    Magnetic Field                          - Tesla
    Species content                         - See user guide for relevant species
    Isotopologue ratios                     - Unitless
    Non-local thermodynamics ratios         - Unitless [pure-style] OR Kelvin [vibrational-style]
    Scattering species content              - See user guide for relevant species
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

  wsv_data["propagation_matrix_jacobian"] = {
      .desc = R"--(
Partial derivative of the *propagation_matrix* with regards to *jacobian_targets*.

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
      .default_value = LinemixingEcsData{},
  };

  wsv_data["frequency_grid"] = {
      .desc =
          R"--(The discrete frequency grid.

Unit:  Hz
)--",
      .type = "AscendingGrid",
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
      .desc = R"--(Atmospheric frequency grids along the propagation path.

See *frequency_grid* for information about the frequency grid

Dimension: [ ppath.np ]

Usage: Output of radiative transfer methods.
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
          R"--(This contains the propagation matrix for the current path point.

The propagation matrix can be used to computed the transmission through a layer
as:

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

  wsv_data["propagation_matrix_select_species"] = {
      .desc = R"--(A select species tag group from *absorption_species*

If set to empty, this selection is void.  It must otherwise match perfectly a tag inside
*absorption_species* for that to be the selection.
)--",
      .type = "SpeciesEnum",
      .default_value = SpeciesEnum::Bath,
  };

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
      .desc = R"--(The surface field describes the surface properties.

This describes the global surface values, such as elevation and 
temperature but also entirerly abstract properties and types.
)--",
      .type = "SurfaceField",
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

  wsv_data["spectral_radiance_unit"] = {
      .desc          = R"--(The spectral radiance unit after conversion.

If a unit conversion is desired, the user has to set this variable to one
of the valid options in *SpectralRadianceUnitType*.

Internally, it is always assumed that this is set to "1" and that no
unit conversion are taking place. 

Please be aware of limitations of follow-up method calls when using this variable manually.
Unless a method or variable explicitly mention that a unit conversion is supported before
it is called, the use of *spectral_radiance_unit* with a different unit than "1" may lead to
undesired results.
)--",
      .type          = "String",
      .default_value = String{"1"},
  };

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
      .default_value = JacobianTargets{},
  };

  wsv_data["ray_path_point"] = {
      .desc = R"--(A single path point.
)--",
      .type = "PropagationPathPoint",
  };

  wsv_data["ray_path"] = {
      .desc = R"--(A list path points making up a propagation path.
)--",
      .type = "ArrayOfPropagationPathPoint",
  };

  wsv_data["spectral_radiance_observer_position"] = {
      .desc = R"--(The position of an observer of spectral radiance.

Most likely only makes sense in combination with *spectral_radiance_observer_line_of_sight*.
)--",
      .type = "Vector3",
  };

  wsv_data["spectral_radiance_observer_line_of_sight"] = {
      .desc = R"--(The position of the observer of spectral radiance.

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
      .desc = R"(A state vector of the model.

In classical ``F(x) = y``-notation, this is the ``x``.
)",
      .type = "Vector",
  };

  wsv_data["model_state_vector_apriori"] = {
      .desc = R"(An apriori state vector of the model.

In classical ``F(x) = y``-notation, this is the ``x``.
)",
      .type = "Vector",
  };

  wsv_data["measurement_vector"] = {
      .desc = R"(The measurment vector for, e.g., a sensor.

In classical ``F(x) = y``-notation, this is the ``y``.
)",
      .type = "Vector",
  };

  wsv_data["measurement_vector_fitted"] = {
      .desc = R"(As *measurement_vector*, but fitted to the model.
)",
      .type = "Vector",
  };

  wsv_data["inversion_iterate_agenda_counter"] = {
      .desc          = R"(A counter for the inversion iterate agenda.
)",
      .type          = "Index",
      .default_value = Index{0},
  };

  wsv_data["inversion_iterate_agenda_do_jacobian"] = {
      .desc = R"(A boolean for if Jacobian calculations should be done.
)",
      .type = "Index",
      .default_value = Index{1},
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
      .desc = R"(The partial derivatives of the *measurement_vector*.

The size of this variable should be the size *measurement_vector* times the size of *model_state_vector*.
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

  wsv_data["disort_optical_thicknesses"] = {
      .desc = R"(A list of single scattering albedos for DISORT.

Note that every row must be convertible to an *AscendingGrid*.

Size is *frequency_grid* times (*ray_path* - 1).

*ray_path* - 1 is the number of layers.
)",
      .type = "Matrix",
  };

  wsv_data["disort_single_scattering_albedo"] = {
      .desc = R"(A list of single scattering albedos for DISORT.

Note that "turned-off" means all zeroes.

Size is *frequency_grid* times (*ray_path* - 1).

*ray_path* - 1 is the number of layers.
)",
      .type = "Matrix",
  };

  wsv_data["disort_fractional_scattering"] = {
      .desc = R"(A list of single scattering albedos for DISORT.

Note that "turned-off" means all zeroes.

Size is *frequency_grid* times (*ray_path* - 1).

*ray_path* - 1 is the number of layers.
)",
      .type = "Matrix",
  };

  wsv_data["disort_legendre_coefficients"] = {
      .desc = R"(A list of single scattering albedos for DISORT.

Note that "turned-off" means all zeroes except in the first polynomial,
which should be all ones.

Size is *frequency_grid* times (*ray_path* - 1) times
all Legendre polynomials.

*ray_path* - 1 is the number of layers.
)",
      .type = "Tensor3",
  };

  wsv_data["disort_positive_boundary_condition"] = {
      .desc = R"(A list of single scattering albedos for DISORT.

Note that "turned-off" means all zeroes.

Size is *frequency_grid* times number of Fourier modes times directional-quadratures.

.. note:: 

    Physical unit of *disort_source_polynomial*, *disort_solar_source*,
    *disort_negative_boundary_condition*, and *disort_positive_boundary_condition*
    must match.
)",
      .type = "Tensor3",
  };

  wsv_data["disort_negative_boundary_condition"] = {
      .desc = R"(A list of single scattering albedos for DISORT.

Note that "turned-off" means all zeroes.

Size is *frequency_grid* times number of Fourier modes times directional-quadratures.

.. note::

    Physical unit of *disort_source_polynomial*, *disort_solar_source*,
    *disort_negative_boundary_condition*, and *disort_positive_boundary_condition*
    must match.
)",
      .type = "Tensor3",
  };

  wsv_data["disort_source_polynomial"] = {
      .desc = R"(A list of single scattering albedos for DISORT.

Note that "turned-off" means that the last dimension size is zero.

Size is *frequency_grid* times (*ray_path* - 1) times
source function polynomial order.

A first order polynomial gives the source function as:

.. math::

    B = s_{1} + s_{0} \tau 

.. note:: 

    Physical unit of *disort_source_polynomial*, *disort_solar_source*,
    *disort_negative_boundary_condition*, and *disort_positive_boundary_condition*
    must match.

*ray_path* - 1 is the number of layers.
)",
      .type = "Tensor3",
  };

  wsv_data["disort_solar_zenith_angle"] = {
      .desc = R"(A list of single scattering albedos for DISORT.

Unit is in degrees.  Note that there's an error if the solar zenith
angle overlaps with one of the fixed DISORT quadrature angles.

Size is *frequency_grid*.
)",
      .type = "Vector",
  };

  wsv_data["disort_solar_azimuth_angle"] = {
      .desc = R"(A list of single scattering albedos for DISORT.

Unit is in degrees.

Size is *frequency_grid*.
)",
      .type = "Vector",
  };

  wsv_data["disort_solar_source"] = {
      .desc = R"(A list of single scattering albedos for DISORT.

Note that "turned-off" means all zeroes.

Size is *frequency_grid*.

.. note:: 

    Physical unit of *disort_source_polynomial*, *disort_solar_source*,
    *disort_negative_boundary_condition*, and *disort_positive_boundary_condition*
    must match.
)",
      .type = "Vector",
  };

  wsv_data["disort_bidirectional_reflectance_distribution_functions"] = {
      .desc = R"(A list of single scattering albedos for DISORT.

Note that "turned-off" means the column dimension size is zero.

Size is *frequency_grid* time number of Fourier modes for the BDRF:s.
)",
      .type = "MatrixOfDisortBDRF",
  };

  wsv_data["disort_spectral_radiance_field"] = {
      .desc = "The spectral radiance field from DISORT.\n",
      .type = "Tensor3",
  };

  wsv_data["disort_quadrature_angles"] = {
      .desc = "The quadrature angles for DISORT.\n",
      .type = "Vector",
  };

  wsv_data["disort_quadrature_weights"] = {
      .desc = "The quadrature weights for DISORT.\n",
      .type = "Vector",
  };

  return wsv_data;
}

std::string_view any(const std::string& type) {
  if (type == "Any") {
    return "T";
  }
  return type;
}
