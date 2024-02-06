#include "workspace_variables.h"

std::unordered_map<std::string, WorkspaceVariableInternalRecord>
internal_workspace_variables() {
  std::unordered_map<std::string, WorkspaceVariableInternalRecord> wsv_data;

  wsv_data["absorption_bands"] = {
      .desc = R"--(Bands of absorption lines for LBL calculations.
)--",
      .type = "AbsorptionBands"};

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
      .type = "ArrayOfCIARecord"};

  wsv_data["propagation_matrix_absorption_lookup"] = {
      .desc = R"--(An absorption lookup table.

It holds an absorption lookup table, as well as all information that
is necessary to use the table to extract absorption. Extraction
routines are implemented as member functions. 

It has quite a complicated structure. For details see the Arts User
Guide section \"The gas absorption lookup table\" or the source code
documentation in gas_abs_lookup.h.
)--",
      .type = "GasAbsLookup"};

  wsv_data["absorption_species"] = {.desc = R"--(Tag groups for gas absorption.

This is an array of arrays of SpeciesTag tag definitions. It defines the
available tag groups for the calculation of scalar gas absorption
coefficients.  See online documentation of method *absorption_speciesSet* for
more detailed information how tag groups work and some examples.
)--",
                                    .type = "ArrayOfArrayOfSpeciesTag"};

  wsv_data["atmospheric_field"] = {.desc =
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
                                   .type = "AtmField"};

  wsv_data["atmospheric_point"] = {.desc =
                                       R"--(An atmospheric point in ARTS.

The atmospheric point consists of all the relevant atmospheric field data
at a discrete point in the atmosphere.  It is often extracted from an *AtmField*
at a single altitude-latitude-longitude but may of course be generated manually.

See *atmospheric_field* for the data that may be available in the atmospheric point.
)--",
                                   .type = "AtmPoint"};

  wsv_data["propagation_matrix_jacobian"] = {.desc = R"--(
Partial derivative of the *propagation_matrix* with regards to *jacobian_targets*.

The units depend on what is set in *jacobian_targets* [1 / m / jacobian target's unit].
)--",
                                             .type = "PropmatMatrix"};

  wsv_data["propagation_matrix_source_vector_nonlte_jacobian"] = {
      .desc =
          R"--(Partial derivative of the *propagation_matrix_source_vector_nonlte* with regards to *jacobian_targets*.

The units are *spectral_radiance_jacobian* per meter.
)--",
      .type = "StokvecMatrix"};

  wsv_data["ecs_data"] = {.desc = R"--(Error corrected sudden data

Dimensions: [num Isotopologues] [num Species]

Used in line-by-line calculations requiring ECS data.
)--",
                          .type = "LinemixingEcsData",
                          .default_value = LinemixingEcsData{}};

  wsv_data["frequency_grid"] = {
      .desc =
          R"--(The frequency grid for monochromatic pencil beam calculations.

Usage: Set by the user. 

Unit:  Hz
)--",
      .type = "AscendingGrid"};

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
      .type = "ArrayOfXsecRecord"};

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
      .type = "StokvecVector"};

  wsv_data["propagation_path_atmospheric_point"] = {
      .desc = R"--(Atmospheric points along the propagation path.

See *atmospheric_point* for information about atmospheric points

Dimension: [ ppath.np ]

Usage: Output of radiative transfer methods.
)--",
      .type = "ArrayOfAtmPoint"};

  wsv_data["propagation_path_frequency_grid"] = {
      .desc = R"--(Atmospheric frequency grids along the propagation path.

See *frequency_grid* for information about the frequency grid

Dimension: [ ppath.np ]

Usage: Output of radiative transfer methods.
)--",
      .type = "ArrayOfAscendingGrid"};

  wsv_data["absorption_predefined_model_data"] = {
      .desc =
          R"--(This contains predefined model data.

Can currently only contain data for new MT CKD models of water.
)--",
      .type = "PredefinedModelData",
      .default_value = PredefinedModelData{}};

  wsv_data["propagation_matrix"] = {
      .desc =
          R"--(This contains the propagation matrix for the current path point.

The propagation matrix can be used to computed the transmission through a layer
as:

.. math:: \mathbf{T} = \exp\left(-\mathbf{K} r\right),

where :math:`\mathbf{K}` is the propagation matrix, and :math:`r` is some distance
over which it is considered constant.

The unit is [1 / m]
)--",
      .type = "PropmatVector"};

  wsv_data["propagation_matrix_select_species"] = {
      .desc = R"--(A select species tag group from *absorption_species*

If set to empty, this selection is void.  It must otherwise match perfectly a tag inside
*absorption_species* for that to be the selection.
)--",
      .type = "SpeciesEnum",
      .default_value = SpeciesEnum::Bath};

  wsv_data["spectral_radiance_background_jacobian"] = {
      .desc = R"--(Spectral radiance derivative from the background
)--",
      .type = "StokvecMatrix"};

  wsv_data["spectral_radiance_background"] = {
      .desc = R"--(Spectral radiance from the background
)--",
      .type = "StokvecVector"};

  wsv_data["background_transmittance"] = {
      .desc = R"--(Transmittance from the background
)--",
      .type = "MuelmatVector"};

  wsv_data["propagation_path_spectral_radiance"] = {
      .desc = R"--(Spectral radiance along the propagation path
)--",
      .type = "ArrayOfStokvecVector"};

  wsv_data["propagation_path_spectral_radiance_jacobian"] = {
      .desc = R"--(Spectral radiance derivative along the propagation path
)--",
      .type = "ArrayOfStokvecMatrix"};

  wsv_data["propagation_path_propagation_matrix"] = {
      .desc = R"--(Propagation matrices along the propagation path
)--",
      .type = "ArrayOfPropmatVector"};

  wsv_data["propagation_path_propagation_matrix_jacobian"] = {
      .desc = R"--(Propagation derivative matrices along the propagation path
)--",
      .type = "ArrayOfPropmatMatrix"};

  wsv_data["propagation_path_propagation_matrix_source_vector_nonlte"] = {
      .desc = R"--(Additional non-LTE along the propagation path
)--",
      .type = "ArrayOfStokvecVector"};

  wsv_data
      ["propagation_path_propagation_matrix_source_vector_nonlte_jacobian"] = {
          .desc = R"--(Additional non-LTE derivative along the propagation path
)--",
          .type = "ArrayOfStokvecMatrix"};

  wsv_data["propagation_path_spectral_radiance_source"] = {
      .desc = R"--(Source vectors along the propagation path
)--",
      .type = "ArrayOfStokvecVector"};

  wsv_data["propagation_path_spectral_radiance_source_jacobian"] = {
      .desc = R"--(Source derivative vectors along the propagation path
)--",
      .type = "ArrayOfStokvecMatrix"};

  wsv_data["propagation_path_transmission_matrix"] = {
      .desc = R"--(Transmission matrices along the propagation path
)--",
      .type = "ArrayOfMuelmatVector"};

  wsv_data["propagation_path_transmission_matrix_cumulative"] = {
      .desc = R"--(Cumulative transmission matrices along the propagation path
)--",
      .type = "ArrayOfMuelmatVector"};

  wsv_data["propagation_path_transmission_matrix_jacobian"] = {
      .desc = R"--(Transmission derivative matrices along the propagation path
)--",
      .type = "ArrayOfArrayOfMuelmatMatrix"};

  wsv_data["surface_field"] = {
      .desc = R"--(The surface field describes the surface properties.

This describes the global surface values, such as elevation and 
temperature but also entirerly abstract properties and types.
)--",
      .type = "SurfaceField"};

  wsv_data["gravity_operator"] = {.desc = R"--(The gravity operator.

Returns gravity in m/s^2 for a given altitude [m], latitude [deg] and longitude [deg].
)--",
                                  .type = "NumericTernaryOperator"};

  wsv_data["spectral_radiance"] = {.desc = R"--(A spectral radiance vector.
)--",
                                   .type = "StokvecVector"};

  wsv_data["spectral_radiance_jacobian"] = {
      .desc = R"--(A spectral radiance derivative matrix.
)--",
      .type = "StokvecMatrix"};

  wsv_data["jacobian_targets"] = {
      .desc = R"--(A list of targets for the Jacobian Matrix calculations.
)--",
      .type = "JacobianTargets",
      .default_value = JacobianTargets{}};

  wsv_data["propagation_path_point"] = {.desc = R"--(A single path point.
)--",
                                        .type = "PropagationPathPoint"};

  wsv_data["propagation_path"] = {
      .desc = R"--(A list path points making up a propagation path.
)--",
      .type = "ArrayOfPropagationPathPoint"};

  wsv_data["spectral_radiance_observer_position"] = {
      .desc = R"--(The position of an observer of spectral radiance.

Most likely only makes sense in combination with *spectral_radiance_observer_line_of_sight*.
)--",
      .type = "Vector3"};

  wsv_data["spectral_radiance_observer_line_of_sight"] = {
      .desc = R"--(The position of the observer of spectral radiance.

Most likely only makes sense in combination with *spectral_radiance_observer_position*.
)--",
      .type = "Vector2"};

  wsv_data["spectral_radiance_operator"] = {
      .desc = R"--(The spectral radiance operator.

This is a class that can compute the spectral radiance
along a path for a single viewing direction and frequency.

It provides several methods to get the path of the spectral
radiance.
)--",
      .type = "SpectralRadianceOperator"};

  return wsv_data;
}
