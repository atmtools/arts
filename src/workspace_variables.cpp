#include "workspace_variables.h"

std::unordered_map<std::string, WorkspaceVariableInternalRecord>
internal_workspace_variables() {
  std::unordered_map<std::string, WorkspaceVariableInternalRecord> wsv_data;

  wsv_data["absorption_bands"] = {
      .desc = R"--(Bands of absorption lines for LBL calculations.
)--",
      .type = "AbsorptionBands"};

  wsv_data["aa_grid"] = {.desc = R"--(Azimuthal angle grid.

The grid must be sorted in increasing order, with no repetitions.

Usage:      Set by the user.

Unit:       degrees
)--",
                         .type = "Vector"};

  wsv_data["abs_cia_data"] = {
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
as that of *abs_species*. There will be CIA data only for those
species that contain a CIA tag, for all other species it will be
empty. The inner array dimension corresponds to the number of CIA tags
for this species (there could be for example both N2-N2 and N2-H2) in
the same species.

The CIA *abs_species* tags are described in *abs_speciesSet*.

Each individual CIARecord holds the complete information from one
HITRAN CIA file. For the given pair of molecules A HITRAN CIA data
file can hold several datasets (data for different temperatures but
fixed frequency range).

Units:

 - Frequencies: Hz
 - Binary absorption cross-sections: m^5*molecule^-2
)--",
      .type = "ArrayOfCIARecord"};

  wsv_data["abs_hitran_relmat_data"] = {
      .desc = R"--(HITRAN line mixing data to compute the relaxation matrix.

This variable holds HITRAN line mixing data
as per J. Lamouroux, L. Realia, X. Thomas, et al., J.Q.S.R.T. 151 (2015), 88-96

It is used for absorption bands with these population tags:
 - ByHITRANFullRelmat
 - ByHITRANRosenkranzRelmat
)--",
      .type = "HitranRelaxationMatrixData"};

  wsv_data["abs_lines"] = {.desc = R"--(A list of spectral line data.
)--",
                           .type = "ArrayOfAbsorptionLines"};

  wsv_data["abs_lines_per_species"] = {
      .desc = R"--(A list of spectral line data for each tag.

Dimensions: [*abs_species*.nelem()][Depends on how many bands there are in *abs_lines*]
)--",
      .type = "ArrayOfArrayOfAbsorptionLines"};

  wsv_data["abs_lookup"] = {.desc = R"--(An absorption lookup table.

It holds an absorption lookup table, as well as all information that
is necessary to use the table to extract absorption. Extraction
routines are implemented as member functions. 

It has quite a complicated structure. For details see the Arts User
Guide section \"The gas absorption lookup table\" or the source code
documentation in gas_abs_lookup.h.
)--",
                            .type = "GasAbsLookup"};

  wsv_data["abs_lookup_is_adapted"] = {
      .desc = R"--(Flag to indicate whether *abs_lookupAdapt* has already been
called.

Values:
 - 0=false - 1=true.
)--",
      .type = "Index"};

  wsv_data["abs_species"] = {.desc = R"--(Tag groups for gas absorption.

This is an array of arrays of SpeciesTag tag definitions. It defines the
available tag groups for the calculation of scalar gas absorption
coefficients.  See online documentation of method *abs_speciesSet* for
more detailed information how tag groups work and some examples.
)--",
                             .type = "ArrayOfArrayOfSpeciesTag"};

  wsv_data["atm_field"] = {.desc =
                               R"--(An atmospheric field in ARTS.

The atmospheric field defines the altitude of the top-of-the-atmosphere,
as well as the variables that are required for the radiative transfer
calculations along a path through the atmosphere.

This atmospheric field may exist on a regular grid but can also be
free-form. Some methods require the atmospheric field to be on a regular
grid, so we provide INSERT-METHOD-NAME-HERE to regularize the field grid.

The atmospheric field may, but does not have to, consist of the following:
    Temperature         - Kelvin
    Pressure            - Pascal
    Wind                - Meters per second
    Magnetic Field      - Tesla
    Species content     - See user guide for relevant species
    NLTE ratios         - Unitless [pure-style] OR Kelvin [vibrational-style]
)--",
                           .type = "AtmField"};

  wsv_data["atm_point"] = {.desc =
                               R"--(An atmospheric point in ARTS.

The atmospheric point consists of all the relevant atmospheric field data
at a discrete point in the atmosphere.  It is often extracted from an *AtmField*
at a single altitude-latitude-longitude but may of course be generated manually.

See *atm_field* for the data that may be available in the atmospheric point.
)--",
                           .type = "AtmPoint"};

  wsv_data["avk"] = {.desc = R"--(Averaging kernel matrix.

This matrix is the partial derivative of the retrieved state vector
with respect to the measurement vector (``y``).

Usage: Used and set by inversion methods.
)--",
                     .type = "Matrix"};

  wsv_data["backend_channel_response"] = {
      .desc = R"--(The response of each backend channel.

The response is given as an *ArrayOfGriddedField1*. The grid consists of
relative frequencies. These relative frequencies are added to 
*f_backend* to obtain the absolute frequency for each response value.
The actual data are the response at each frequency grid point.

There are here two options. If the array has length 1, the same
response is applied for all channels. Accordingly, this assumes that
all channels have the same response function. The second option is to
specify the response for each channel seperately. This signifies that
the *backend_channel_response* array has either 1 or n elements, where
n is the length of *f_backend*

Usage: Set by the user.

Size:

- Array[N_ch]

  - GriddedField1

    - [N_f]
    - [N_f]
)--",
      .type = "ArrayOfGriddedField1"};

  wsv_data["complex_refr_index"] = {
      .desc = R"--(Complex refractive index (n) data.

The variable works as a lookup-table of complex refractive index.
The matter type (water, ice ...) is unspecified, it is up to the
user to fill the variable with data for the expected matter.
This variable type can be used to describe n of both the surface and
atmospheric particles.

The column dimension has always size 2, where the first and second
column holds the real and imaginary part of n, respectively. The row
dimension matches temperature, and the page dimension is frequency.
Both the temperature and frequency dimensions grids are allowed to
have length 1, which is interpreted as n being constant in that
dimension.

When mapping these data to the required frequencies and temperatures
a bi-linear interpolation is applied.

Unit:       -

Dimensions: 
 - Vector f_grid[N_f]
 - Vector T_grid[N_T]
 - Tensor3 data[N_f][N_T][2]
)--",
      .type = "ComplexGriddedField2"};

  wsv_data["dpropmat_clearsky_dx"] = {.desc = R"--( FIXMEDOC
Partial derivative of absorption coefficients.

This contains the partial derivative of absorption coefficients for
one point in the atmosphere (one set of pressure, temperature, zn
magnetic field, and VMR values) with respect to one of the input
parameters.

Dimension: [ n_quantities ] [naa, nza, nf, f(stokes_dim)]

*jacobian_targets* should be used to set the input variable for
partial derivation

Unit: 1/m/jacobian_targets
)--",
                                      .type = "PropmatMatrix"};

  wsv_data["dnlte_source_dx"] = {
      .desc =
          R"--(NLTE partial derivatives output is two parts:  S * dB/dx + dS/dx * B.

Dimensions: [ quantities ] [nza, naa, nf, stokes_dim] or [0]

Unit: 1/m/jacobian_targets
)--",
      .type = "StokvecMatrix"};

  wsv_data["ecs_data"] = {.desc = R"--(Error corrected sudden data

Dimensions: [num IDs] [num Species]

It is used for absorption bands with these population tags:
 - ByMakarovFullRelmat 
 - ByRovibLinearDipoleLineMixing
)--",
                          .type = "MapOfErrorCorrectedSuddenData"};

  wsv_data["ecs_data2"] = {.desc = R"--(Error corrected sudden data

Dimensions: [num Isotopologues] [num Species]

Used in line-by-line calculations requiring ECS data.
)--",
                           .type = "LinemixingEcsData",
                           .default_value = LinemixingEcsData{}};

  wsv_data["f_backend"] = {
      .desc =
          R"--(The frequency position of each backend (spectrometer) channel.

Usage: Set by the user. 

Unit:  Hz
)--",
      .type = "Vector"};

  wsv_data["f_grid"] = {
      .desc =
          R"--(The frequency grid for monochromatic pencil beam calculations.

Usage: Set by the user. 

Unit:  Hz
)--",
      .type = "Vector"};

  wsv_data["g0"] = {.desc = R"--(Gravity at zero altitude.

This variable is \"little g\" at the reference ellipsiod. That is,
for Earth this is a value around 9.81 m/s2
)--",
                    .type = "Numeric"};

  wsv_data["xsec_fit_data"] = {
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

  wsv_data["irradiance_field"] = {
      .desc = R"--(Irradiance field also known as flux density.

Radiant flux received by a surface per unit area for each hemisphere.
The last dimension denotes the hemispheres. The first component is
the downward irradiance and the second component is the upward irradiance

Units: W m^-2

Size: [ p_grid,  lat_grid,  lon_grid,  2 ]
)--",
      .type = "Tensor4"};

  wsv_data["nlte_do"] = {.desc = R"--(Flag to perform Non-LTE calculations.
)--",
                         .type = "Index",
                         .default_value = Index{0}};

  wsv_data["nlte_source"] = {
      .desc =
          R"--(Variable to contain the additional source function due to NLTE effects.

Dimensions: [nza, naa, nf, stokes_dim]
)--",
      .type = "StokvecVector"};

  wsv_data["ppvar_atm"] = {
      .desc = R"--(Atmospheric points along the propagation path.

ppvar stands for propagation path variable. The variables named in is
way describe the atmosphere and its properties at each point of the
propagation path

See *atm_point* for information about atmospheric points

Dimension: [ ppath.np ]

Usage: Output of radiative transfer methods.
)--",
      .type = "ArrayOfAtmPoint"};

  wsv_data["ppvar_f"] = {
      .desc = R"--(Atmospheric frequency grids along the propagation path.

ppvar stands for propagation path variable. The variables named in is
way describe the atmosphere and its properties at each point of the
propagation path

See *f_grid* for information about the frequency grid

Dimension: [ ppath.np ]

Usage: Output of radiative transfer methods.
)--",
      .type = "ArrayOfVector"};

  wsv_data["predefined_model_data"] = {
      .desc =
          R"--(This contains predefined model data.

Can currently only contain data for new MT CKD models of water.
)--",
      .type = "PredefinedModelData",
      .default_value = PredefinedModelData{}};

  wsv_data["propmat_clearsky"] = {
      .desc =
          R"--(This contains the absorption coefficients for one point in the atmosphere.

Dimensions: [naa, nza, nf, f(stokes_dim)]

Unit: 1/m
)--",
      .type = "PropmatVector"};

  wsv_data["propmat_clearsky_agenda_checked"] = {
      .desc = R"--(OK-flag for *propmat_clearsky_agenda*.

Set by *propmat_clearsky_agenda_checkedCalc*.
)--",
      .type = "Index",
      .default_value = Index{0}};

  wsv_data["select_abs_species"] = {
      .desc = R"--(A select species tag group from *abs_species*

If set to empty, this selection is void.  It must otherwise match perfectly a tag inside
*abs_species* for that to be the selection.
)--",
      .type = "ArrayOfSpeciesTag",
      .default_value = ArrayOfSpeciesTag{}};

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

  wsv_data["spectral_radiance_path"] = {
      .desc = R"--(Radiation along the propagation path
)--",
      .type = "ArrayOfStokvecVector"};

  wsv_data["spectral_radiance_path_jacobian"] = {
      .desc = R"--(Radiation derivative along the propagation path
)--",
      .type = "ArrayOfStokvecMatrix"};

  wsv_data["ppvar_propmat"] = {
      .desc = R"--(Propagation matrices along the propagation path
)--",
      .type = "ArrayOfPropmatVector"};

  wsv_data["ppvar_dpropmat"] = {
      .desc = R"--(Propagation derivative matrices along the propagation path
)--",
      .type = "ArrayOfPropmatMatrix"};

  wsv_data["ppvar_nlte"] = {
      .desc = R"--(Additional NLTE along the propagation path
)--",
      .type = "ArrayOfStokvecVector"};

  wsv_data["ppvar_dnlte"] = {
      .desc = R"--(Additional NLTE derivative along the propagation path
)--",
      .type = "ArrayOfStokvecMatrix"};

  wsv_data["spectral_radiance_path_source"] = {
      .desc = R"--(Source vectors along the propagation path
)--",
      .type = "ArrayOfStokvecVector"};

  wsv_data["spectral_radiance_path_source_jacobian"] = {
      .desc = R"--(Source derivative vectors along the propagation path
)--",
      .type = "ArrayOfStokvecMatrix"};

  wsv_data["ppvar_tramat"] = {
      .desc = R"--(Transmission matrices along the propagation path
)--",
      .type = "ArrayOfMuelmatVector"};

  wsv_data["ppvar_cumtramat"] = {
      .desc = R"--(Cumulative transmission matrices along the propagation path
)--",
      .type = "ArrayOfMuelmatVector"};

  wsv_data["ppvar_dtramat"] = {
      .desc = R"--(Transmission derivative matrices along the propagation path
)--",
      .type = "ArrayOfArrayOfMuelmatMatrix"};

  wsv_data["sensor_response_f"] = {
      .desc =
          R"--(The frequencies associated with the output of ``sensor_response``.

This vector gives the frequency for each element of the measurement
vector produced inside one measurement block. The frequencies of
the total measurement vector, ``y``, are obtained by repeating these
frequencies n times, where n is the number of measurement blocks.

The variable shall not be set manually, it will be set together with
``sensor_response`` by sensor response WSMs.

Usage: Set by sensor response methods.

Unit:  [ Hz ]
)--",
      .type = "Vector"};

  wsv_data["sensor_response_f_grid"] = {
      .desc = R"--(The frequency grid associated with ``sensor_response``.

A variable for communication between sensor response WSMs. Matches
initially *f_grid*, but is later adjusted according to the sensor
specifications. Only defined when a common grid exists. Values are
here not repeated as in *sensor_response_f*

Usage: Set by sensor response methods.

Unit:  [ Hz ]
)--",
      .type = "Vector"};

  wsv_data["spectral_radiance_field"] = {.desc = R"--(Spectral radiance field.

This variable holds a calculation of the radiance field through
the atmosphere, for the directions matching *za_grid* and *aa_grid*.

Units: W / (m^2 Hz sr)

 Size: [f_grid, p_grid,  lat_grid,  lon_grid,  za_grid, aa_grid, stokes_dim ]

Note:
 For 1D, the size of the latitude, longitude and azimuth
 dimension (N_aa) are all 1.
)--",
                                         .type = "Tensor7"};

  wsv_data["suns"] = {.desc = R"--(Array of Sun.

This variable describes a list of suns.
Each sun is described by a struct with its spectrum, radius,
distance from center of planet to center of sun,
temperature (if possible), latitude in the sky of the planet,
longitude in the sky of the planet and the type
)--",
                      .type = "ArrayOfSun",
                      .default_value = ArrayOfSun{}};

  wsv_data["surface_emission"] = {.desc = R"--(The emission from the surface.

See specific methods generating *surface_emission* and the user
guide for more information.

Dimensions: [ f_grid, stokes_dim ]
)--",
                                  .type = "Matrix"};

  wsv_data["surface_field"] = {
      .desc = R"--(The surface field describes the surface properties.

This describes the global surface values, such as elevation and 
temperature but also entirerly abstract properties and types.
)--",
      .type = "SurfaceField"};

  wsv_data["surface_rmatrix"] = {
      .desc = R"--(The reflection coefficients for the surface
      
The directions are given by surface line-of-sight to the direction of interest.

The rows and columns of this tensor holds the reflection
coefficient matrix for one frequency and one LOS. The reflection
coefficients shall take into accound the angular weighting of the
downwelling radiation.

See specific methods generating *surface_rmatrix* and the user guide
for more information.

Units:      -

Dimensions: [ surface_los, f_grid, stokes_dim, stokes_dim ]
)--",
      .type = "Tensor4"};

  wsv_data["time"] = {.desc = R"--(A UTC time point.
)--",
                      .type = "Time"};

  wsv_data["time_grid"] = {.desc = R"--(A grid of times.  Should be increasing
)--",
                           .type = "ArrayOfTime"};

  wsv_data["wigner_initialized"] = {
      .desc = R"--(Indicates if the wigner tables are initialized.
If they are not, computations will be aborted.

Will hold the value of provided maximum factorial value

The developer should always test this variable in functions
that might require computing wigner symbols because the error
handling is otherwise offloaded to third party software...
)--",
      .type = "Index"};

  wsv_data["wmrf_weights"] = {
      .desc = R"--(The weights for a WMRF fast calculation.

Weights are stored in a sparse matrix. This can be used as a
sensor_response matrix.

The dimension of the matrix is (nchan, nfreq), where nchan
is the number of instrument channels and nfreq is the number
of monochromatic frequencies.
)--",
      .type = "Sparse"};

  wsv_data["za_grid"] = {.desc = R"--(Zenith angle grid.

The grid must be sorted in increasing order, with no repetitions.

Usage:      Set by the user.

Unit:       degrees
)--",
                         .type = "Vector"};

  wsv_data["za_grid_weights"] = {.desc = R"--(Zenith angle integration weights.

The integration weight are needed for calculation of radiation fluxes

Unit:  unitless
)--",
                                 .type = "Vector"};

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

  wsv_data["path_point"] = {.desc = R"--(A single path point.
)--",
                            .type = "PropagationPathPoint"};

  wsv_data["propagation_path"] = {
      .desc = R"--(A list path points making up a propagation path.
)--",
      .type = "ArrayOfPropagationPathPoint"};

  return wsv_data;
}
