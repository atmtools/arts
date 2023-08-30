#include "workspace_methods.h"

#include <limits>
#include <optional>

#pragma clang optimize off

std::unordered_map<std::string, WorkspaceMethodInternalRecord>
internal_workspace_methods() {
  std::unordered_map<std::string, WorkspaceMethodInternalRecord> wsm_data;

  wsm_data["AltLatLonFieldSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Fills an altitude-latitude-longitude field with given input.

Grids and data must match in size.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"gfield3"},
      .gout_type = {"GriddedField3"},
      .gout_desc = {R"--(Field to set.)--"},

      .gin =
          {"altitude_grid", "latitude_grid", "longitude_grid", "data", "name"},
      .gin_type = {"Vector", "Vector", "Vector", "Tensor3", "String"},
      .gin_value =
          {std::nullopt, std::nullopt, std::nullopt, std::nullopt, String("")},
      .gin_desc = {R"--(The altitude grid of ``data``.)--",
                   R"--(The latitude grid of ``data``.)--",
                   R"--(The longitude grid of ``data``.)--",
                   R"--(The data of the field (will become gfield2.data).)--",
                   R"--(The name of the field (will become gfield2.name).)--"},

  };

  wsm_data["AltLatLonFieldSetConstant"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets an altitude-latitude-longitude field to have a constant data value.

All three grids grids are set to have length one, with the value 0.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"gfield3"},
      .gout_type = {"GriddedField3"},
      .gout_desc = {R"--(Field to set.)--"},

      .gin = {"value", "name"},
      .gin_type = {"Numeric", "String"},
      .gin_value = {std::nullopt, String("")},
      .gin_desc = {R"--(The value (to place in gfield3.data).)--",
                   R"--(The name of the field (will become gfield3.name).)--"},

  };

  wsm_data["AngularGridsSetFluxCalc"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets the angular grids for the calculation of radiation fluxes. 

This method sets the angular grids for the radiation fluxes type
calculations and calculates the integration weights *za_grid_weights*
for the zenith angle integration. For down- und up-looking geometries
it suffices to use the default values of N_za_grid and N_aa_grid.
From N_aa_grid an equally spaced grid is created and stored in the
WSV *aa_grid*.

Depending on the desired za_grid_type *za_grid* will be
equally spaced ('linear') or unequally ('linear_mu','double_gauss')
Important, N_za_grid must be an even number because for the 
integration over each hemisphere N_za_grid / 2 zenith angles are needed.

Possible zenith angle grid types are:

- ``double_gauss``:
  The zenith grid and the integration weights are set according
  to a gauss-legendre integration for each hemispheres.
- ``linear``: Equally space grid between 0 deg and 180 deg including the poles
- ``linear_mu``:
  Similar to 'linear' but equally spaced for cos(180 deg) to cos(0 deg),
  which results a unequally spaced angular grid
)--",
      .author = {"Manfred Brath"},
      .out = {"za_grid", "aa_grid", "za_grid_weights"},

      .gin = {"N_za_grid", "N_aa_grid", "za_grid_type"},
      .gin_type = {"Index", "Index", "String"},
      .gin_value = {Index{2}, Index{1}, String("linear_mu")},
      .gin_desc = {R"--(Number of zenith angles)--",
                   R"--(Number of azimuth angles)--",
                   R"--(Zenith angle grid type)--"},

  };

  wsm_data["AntennaMultiBeamsToPencilBeams"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Maps a multi-beam case to a matching pencil beam case.

Cases with overlapping beams are most efficiently handled by
letting *antenna_dlos* have several rows. That is, there are
multiple beams for each measurement block. The drawback is that
many variables must be adjusted if the corresponding pencil beam
spectra shall be calculated. This method makes this adjustment.
That is, if you have a control file for a multiple beam case and
for some reason want to avoid the antenna weighting, you add this
method before *sensor_responseInit*, and remove the call of
*sensor_responseAntenna* and you will get the matching pencil beam
spectra.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_pos",
              "sensor_los",
              "antenna_dlos",
              "antenna_dim",
              "mblock_dlos"},

      .in = {"sensor_pos",
             "sensor_los",
             "antenna_dlos",
             "antenna_dim",
             "mblock_dlos"},

  };

  wsm_data["AntennaOff"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets some antenna related variables

Use this method to set *antenna_dim* and *mblock_dlos* to
suitable values (1 and [0], respectively) for cases when a
sensor is included, but the antenna pattern is neglected.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"antenna_dim", "mblock_dlos"},

  };

  wsm_data["Append"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Append one workspace variable to another.

This method can append an array to another array of the same type,
e.g. ArrayOfIndex to ArrayOfIndex. Or a single element to an array
such as a Tensor3 to an ArrayOfTensor3.

Appending two vectors or a numeric to a vector works as for array
variables.

Both another matrix or a vector can be appended to a matrix. In
addition, for matrices, the 'append dimension' can be selected.
The third argument, ``dimension``, indicates how to append, where
"leading" means to append row-wise, and "trailing" means
column-wise.

Other types (TensorX) are currently only implemented for
appending to the leading dimension.

This method is not implemented for all types, just for those that
were thought or found to be useful. (See variable list below.).
)--",
      .author = {"Stefan Buehler, Oliver Lemke"},

      .gout = {"output"},
      .gout_type =
          {"Vector, Vector, Matrix, Matrix, Tensor3, Tensor3, Tensor4, Tensor4, String, ArrayOfAbsorptionLines, ArrayOfAgenda, ArrayOfArrayOfAbsorptionLines, ArrayOfArrayOfGriddedField1, ArrayOfArrayOfGriddedField2, ArrayOfArrayOfGriddedField3, ArrayOfArrayOfIndex, ArrayOfArrayOfMatrix, ArrayOfArrayOfMuelmatMatrix, ArrayOfArrayOfMuelmatVector, ArrayOfArrayOfPropmatMatrix, ArrayOfArrayOfPropmatVector, ArrayOfArrayOfScatteringMetaData, ArrayOfArrayOfSingleScatteringData, ArrayOfArrayOfSpeciesTag, ArrayOfArrayOfStokvecMatrix, ArrayOfArrayOfStokvecVector, ArrayOfArrayOfString, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6, ArrayOfArrayOfTime, ArrayOfArrayOfVector, ArrayOfAtmPoint, ArrayOfCIARecord, ArrayOfGriddedField1, ArrayOfGriddedField2, ArrayOfGriddedField3, ArrayOfGriddedField4, ArrayOfIndex, ArrayOfJacobianTarget, ArrayOfMatrix, ArrayOfMuelmatMatrix, ArrayOfMuelmatVector, ArrayOfPpath, ArrayOfPropmatMatrix, ArrayOfPropmatVector, ArrayOfQuantumIdentifier, ArrayOfRetrievalQuantity, ArrayOfScatteringMetaData, ArrayOfSingleScatteringData, ArrayOfSparse, ArrayOfSpeciesTag, ArrayOfStokvecMatrix, ArrayOfStokvecVector, ArrayOfString, ArrayOfSun, ArrayOfTelsemAtlas, ArrayOfTensor3, ArrayOfTensor4, ArrayOfTensor5, ArrayOfTensor6, ArrayOfTensor7, ArrayOfTime, ArrayOfVector, ArrayOfXsecRecord, ArrayOfAbsorptionLines, ArrayOfAgenda, ArrayOfArrayOfAbsorptionLines, ArrayOfArrayOfGriddedField1, ArrayOfArrayOfGriddedField2, ArrayOfArrayOfGriddedField3, ArrayOfArrayOfIndex, ArrayOfArrayOfMatrix, ArrayOfArrayOfMuelmatMatrix, ArrayOfArrayOfMuelmatVector, ArrayOfArrayOfPropmatMatrix, ArrayOfArrayOfPropmatVector, ArrayOfArrayOfScatteringMetaData, ArrayOfArrayOfSingleScatteringData, ArrayOfArrayOfSpeciesTag, ArrayOfArrayOfStokvecMatrix, ArrayOfArrayOfStokvecVector, ArrayOfArrayOfString, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6, ArrayOfArrayOfTime, ArrayOfArrayOfVector, ArrayOfAtmPoint, ArrayOfCIARecord, ArrayOfGriddedField1, ArrayOfGriddedField2, ArrayOfGriddedField3, ArrayOfGriddedField4, ArrayOfIndex, ArrayOfJacobianTarget, ArrayOfMatrix, ArrayOfMuelmatMatrix, ArrayOfMuelmatVector, ArrayOfPpath, ArrayOfPropmatMatrix, ArrayOfPropmatVector, ArrayOfQuantumIdentifier, ArrayOfScatteringMetaData, ArrayOfSingleScatteringData, ArrayOfSparse, ArrayOfStokvecMatrix, ArrayOfStokvecVector, ArrayOfString, ArrayOfTelsemAtlas, ArrayOfTensor3, ArrayOfTensor4, ArrayOfTensor5, ArrayOfTensor6, ArrayOfTensor7, ArrayOfTime, ArrayOfVector"},
      .gout_desc = {R"--(The variable to append to.)--"},

      .gin = {"input", "dimension"},
      .gin_type = {"Numeric, Vector, Matrix, Vector, Matrix, Tensor3, Tensor3, Tensor4, String, ArrayOfAbsorptionLines, ArrayOfAgenda, ArrayOfArrayOfAbsorptionLines, ArrayOfArrayOfGriddedField1, ArrayOfArrayOfGriddedField2, ArrayOfArrayOfGriddedField3, ArrayOfArrayOfIndex, ArrayOfArrayOfMatrix, ArrayOfArrayOfMuelmatMatrix, ArrayOfArrayOfMuelmatVector, ArrayOfArrayOfPropmatMatrix, ArrayOfArrayOfPropmatVector, ArrayOfArrayOfScatteringMetaData, ArrayOfArrayOfSingleScatteringData, ArrayOfArrayOfSpeciesTag, ArrayOfArrayOfStokvecMatrix, ArrayOfArrayOfStokvecVector, ArrayOfArrayOfString, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6, ArrayOfArrayOfTime, ArrayOfArrayOfVector, ArrayOfAtmPoint, ArrayOfCIARecord, ArrayOfGriddedField1, ArrayOfGriddedField2, ArrayOfGriddedField3, ArrayOfGriddedField4, ArrayOfIndex, ArrayOfJacobianTarget, ArrayOfMatrix, ArrayOfMuelmatMatrix, ArrayOfMuelmatVector, ArrayOfPpath, ArrayOfPropmatMatrix, ArrayOfPropmatVector, ArrayOfQuantumIdentifier, ArrayOfRetrievalQuantity, ArrayOfScatteringMetaData, ArrayOfSingleScatteringData, ArrayOfSparse, ArrayOfSpeciesTag, ArrayOfStokvecMatrix, ArrayOfStokvecVector, ArrayOfString, ArrayOfSun, ArrayOfTelsemAtlas, ArrayOfTensor3, ArrayOfTensor4, ArrayOfTensor5, ArrayOfTensor6, ArrayOfTensor7, ArrayOfTime, ArrayOfVector, ArrayOfXsecRecord, AbsorptionLines, Agenda, ArrayOfAbsorptionLines, ArrayOfGriddedField1, ArrayOfGriddedField2, ArrayOfGriddedField3, ArrayOfIndex, ArrayOfMatrix, ArrayOfMuelmatMatrix, ArrayOfMuelmatVector, ArrayOfPropmatMatrix, ArrayOfPropmatVector, ArrayOfScatteringMetaData, ArrayOfSingleScatteringData, ArrayOfSpeciesTag, ArrayOfStokvecMatrix, ArrayOfStokvecVector, ArrayOfString, ArrayOfTensor3, ArrayOfTensor6, ArrayOfTime, ArrayOfVector, AtmPoint, CIARecord, GriddedField1, GriddedField2, GriddedField3, GriddedField4, Index, JacobianTarget, Matrix, MuelmatMatrix, MuelmatVector, Ppath, PropmatMatrix, PropmatVector, QuantumIdentifier, ScatteringMetaData, SingleScatteringData, Sparse, StokvecMatrix, StokvecVector, String, TelsemAtlas, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7, Time, Vector", "String"},
      .gin_value = {std::nullopt, String("leading")},
      .gin_desc =
          {R"--(The variable to append.)--",
           R"--(Where to append. Could be either the "leading" or "trailing" dimension.)--"},

      .pass_names = true};

  wsm_data["ArrayOfGriddedFieldGetNames"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Get the names of all GriddedFields stored in an Array.

See *GriddedFieldGetName*.
)--",
      .author = {"Lukas Kluft"},

      .gout = {"names"},
      .gout_type = {"ArrayOfString"},
      .gout_desc =
          {R"--(Names of the GriddedFields in the ArrayOfGriddedField.)--"},

      .gin = {"griddedfields"},
      .gin_type =
          {"ArrayOfGriddedField1, ArrayOfGriddedField2, ArrayOfGriddedField3, ArrayOfGriddedField4"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Array of GriddedFields.)--"},

  };

  wsm_data["ArrayOfIndexLinSpace"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initializes an ArrayOfIndex with linear spacing.

The first element equals always the start value, and the spacing
equals always the step value, but the last value can deviate from
the stop value. ``step`` can be both positive and negative.

The created array is [start, start+step, start+2*step, ...]
)--",
      .author = {"Oliver Lemke"},

      .gout = {"output"},
      .gout_type = {"ArrayOfIndex"},
      .gout_desc = {R"--(Output array.)--"},

      .gin = {"start", "stop", "step"},
      .gin_type = {"Index", "Index", "Index"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Start value.)--",
                   R"--(Maximum/minimum value of the end value)--",
                   R"--(Spacing of the array.)--"},

  };

  wsm_data["ArrayOfIndexSetConstant"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Creates an ArrayOfIndex of length *nelem*, with all values
identical.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"ArrayOfIndex"},
      .gout_desc = {R"--(Variable to initialize.)--"},
      .in = {"nelem"},
      .gin = {"value"},
      .gin_type = {"Index"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Array value..)--"},

  };

  wsm_data["ArrayOfQuantumIdentifierFromLines"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets an ArrayOfQuantumIdentifier to all levels in *abs_lines_per_species*
with defined quantum numbers

Lines without defined quantum numbers are ignored
)--",
      .author = {"Richard Larsson"},

      .gout = {"output"},
      .gout_type = {"ArrayOfQuantumIdentifier"},
      .gout_desc =
          {R"--(Identifiers to all levels in *abs_lines_per_species*)--"},
      .in = {"abs_lines_per_species"},
      .gin = {"global"},
      .gin_type = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc = {R"--(Only look at global quantum numbers)--"},

  };

  wsm_data["ArrayOfTimeNLinSpace"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Creates a time array with length *nelem*, equally spaced between the
given end values.

The length (*nelem*) must be larger than 1.
)--",
      .author = {"Richard Larsson"},

      .gout = {"output"},
      .gout_type = {"ArrayOfTime"},
      .gout_desc = {R"--(Variable to initialize.)--"},
      .in = {"nelem"},
      .gin = {"start", "stop"},
      .gin_type = {"String", "String"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Start value.)--", R"--(End value.)--"},

  };

  wsm_data["ArrayOfTimeSetConstant"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Creates an ArrayOfTime and sets all elements to the specified value.

The vector length is determined by *nelem*.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"ArrayOfTime"},
      .gout_desc = {R"--(Variable to initialize.)--"},
      .in = {"nelem"},
      .gin = {"value"},
      .gin_type = {"Time"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Time value.)--"},

  };

  wsm_data["AtmFieldPRegrid"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Interpolates the input field along the pressure dimension from
``p_grid_old`` to to ``p_grid_new``.

Extrapolation is allowed within the common 0.5grid-step margin.
in and out fields can be the same variable.
)--",
      .author = {"Jana Mendrok"},

      .gout = {"output"},
      .gout_type = {"Tensor3, Tensor4"},
      .gout_desc = {R"--(Regridded atmospheric field.)--"},

      .gin = {"input", "p_grid_new", "p_grid_old", "interp_order"},
      .gin_type = {"Tensor3, Tensor4", "Vector", "Vector", "Index"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, Index{1}},
      .gin_desc = {R"--(Input atmospheric field.)--",
                   R"--(Pressure grid to regrid to)--",
                   R"--(Pressure grid of input field)--",
                   R"--(Interpolation order.)--"},

  };

  wsm_data["AtmFieldsAndParticleBulkPropFieldFromCompact"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Extract pressure grid and atmospheric fields from
*atm_fields_compact*.

An atmospheric scenario includes the following data for each
position (pressure, latitude, longitude) in the atmosphere:

1. temperature field
2. the corresponding altitude field
3. vmr fields for the gaseous species
4. scattering species fields

This method splits up the data found in *atm_fields_compact* to
p_grid, lat_grid, lon_grid, vmr_field, particle_bulkprop_field,
and particle_bulkprop_names.
See documentation of *atm_fields_compact* for a definition of the
data.

Compact states are characterized by having all atmospheric fields
already given on identical grids. That is, no interpolation needs
to be and is performed. Keyword ``p_min`` allows to remove atmospheric
levels with pressures lower than the given value (default: no
removal). This reduces computational burden and is useful when
upper atmospheric contributions are negligible.

Possible future extensions: Add a keyword parameter to refine the
pressure grid if it is too coarse. Or a version that interpolates
onto given grids, instead of using and returning the original grids.
)--",
          .author = {"Jana Mendrok, Manfred Brath"},
          .out = {"atm_field", "particle_bulkprop_names"},

          .in = {"abs_species", "atm_fields_compact"},
          .gin = {"delim", "check_gridnames"},
          .gin_type = {"String", "Index"},
          .gin_value = {String("-"), Index{0}},
          .gin_desc =
              {R"--(Delimiter string of *scat_species* elements.)--",
               R"--(A flag with value 1 or 0. If set to one, the gridnames of the *atm_fields_compact* are checked.)--"},

      };

  wsm_data["CIAInfo"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Display information about the given CIA tags.
The CIA tags shown are in the same format as needed by *abs_speciesSet*.
)--",
      .author = {"Oliver Lemke"},

      .gin = {"catalogpath", "cia_tags"},
      .gin_type = {"String", "ArrayOfString"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(Path to the CIA catalog directory.)--",
           R"--(Array of CIA tags to view, e.g. [ "N2-N2", "H2-H2" ])--"},

  };

  wsm_data["CIARecordReadFromFile"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads CIARecord from Hitran-style file.
)--",
      .author = {"Richard Larsson"},

      .gout = {"cia_record"},
      .gout_type = {"CIARecord"},
      .gout_desc = {R"--(CIARecord type variable for input and output.)--"},

      .gin = {"species_tag", "filename"},
      .gin_type = {"String", "String"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(SpeciesTag string to associate with this CIARecord. See *abs_speciesSet* for correct format.)--",
           R"--(Filename of HITRAN CIA data file.)--"},

  };

  wsm_data["CallbackFunctionExecute"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Execute any code in Arts
)--",
      .author = {"Richard Larsson"},

      .gin = {"function"},
      .gin_type = {"CallbackFunction"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(This will execute as "function(current workspace);")--"},
      .pass_workspace = true,

  };

  wsm_data["CheckUnique"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Checks that *abs_lines* contains only unique absorption lines
)--",
      .author = {"Richard Larsson"},

      .in = {"abs_lines"},

  };

  wsm_data["Compare"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Checks the consistency between two variables.

The two variables are checked to not deviate outside the specified
value (``maxabsdiff``). An error is issued if this is not fulfilled.

The main application of this method is to be part of the test
control files, and then used to check that a calculated value
is consistent with an old, reference, value.
)--",
      .author = {"Oliver Lemke"},

      .gin = {"var1", "var2", "maxabsdiff", "error_message"},
      .gin_type =
          {"Numeric, Vector, Matrix, Tensor3, Tensor4, Tensor5, Tensor7, ArrayOfVector, ArrayOfMatrix, ArrayOfTensor7, GriddedField3, Sparse, SingleScatteringData",
           "Numeric, Vector, Matrix, Tensor3, Tensor4, Tensor5, Tensor7, ArrayOfVector, ArrayOfMatrix, ArrayOfTensor7, GriddedField3, Sparse, SingleScatteringData",
           "Numeric",
           "String"},
      .gin_value = {std::nullopt, std::nullopt, Numeric{}, String("")},
      .gin_desc = {R"--(A first variable)--",
                   R"--(A second variable)--",
                   R"--(Threshold for maximum absolute difference.)--",
                   R"--(Additional error message.)--"},

      .pass_names = true};

  wsm_data["CompareRelative"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Checks the consistency between two variables by their relative values.

The two variables are checked to not deviate outside the specified
relative value (``maxabsreldiff``). An error is issued if this is not
fulfilled.

The main application of this method is to be part of the test
control files, and then used to check that a calculated value
is consistent with an old, reference, value.

If either value is 0.0, the relative error is considered as 0
for easier use.  This really means infinite differences, though
allowing zero-crossings is useful for plenty of tests. So Be Aware!

If both ``var1`` and ``var2`` are non-zero, the difference is evaluated
as: abs(var1/var2-1)
That is, ``var2`` is taken as the reference value.
)--",
      .author = {"Oliver Lemke", "Richard Larsson"},

      .gin = {"var1", "var2", "maxabsreldiff", "error_message"},
      .gin_type =
          {"Numeric, Vector, Matrix, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7, ArrayOfVector, ArrayOfMatrix, ArrayOfTensor3, ArrayOfTensor4, ArrayOfTensor6, ArrayOfTensor7, ArrayOfArrayOfVector, ArrayOfArrayOfMatrix, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6",
           "Numeric, Vector, Matrix, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7, ArrayOfVector, ArrayOfMatrix, ArrayOfTensor3, ArrayOfTensor4, ArrayOfTensor6, ArrayOfTensor7, ArrayOfArrayOfVector, ArrayOfArrayOfMatrix, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6",
           "Numeric",
           "String"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, String("")},
      .gin_desc = {R"--(A first variable)--",
                   R"--(A second variable)--",
                   R"--(Threshold for maximum relative difference.)--",
                   R"--(Additional error message.)--"},

      .pass_names = true};

  wsm_data["Copy"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Copy a workspace variable.

This method can copy any workspace variable
to another workspace variable of the same group. (E.g., a Matrix to
another Matrix.)

As always, output comes first in the argument list!

Usage example:

Copy(f_grid, p_grid)

Will copy the content of ``p_grid`` to *f_grid*. The size of *f_grid*
is adjusted automatically (the normal behaviour for workspace
methods).
)--",
      .author = {"Stefan Buehler"},

      .gout = {"output"},
      .gout_type = {"Any"},
      .gout_desc = {R"--(Destination variable.)--"},

      .gin = {"input"},
      .gin_type = {"Any"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Source variable.)--"},

      .pass_names = true};

  wsm_data["DOAngularGridsSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the angular grids for Discrete Ordinate type scattering
calculations.

This method sets the angular grids for the Discrete Ordinate type
scattering calculations (DOIT, DISORT). For down- und up-looking
geometries it suffices to define ``N_za_grid`` (both solvers) and
``N_aa_grid`` (DOIT). From these numbers equally spaced grids are
created and stored in the WSVs *za_grid* and *aa_grid*.

For limb simulations it is important to use an optimized zenith
angle grid with a very fine resolution around the horizon
(za=90 degrees). Such a grid can be generated using
*doit_za_grid_optCalc*. To be applied, the name of the file holding
the optimized angle grid has to be given (``za_grid_opt_file``).

When an optimized grid is present, the equidistant grid is used for
the calculation of the scattering integrals, while the optimized
grid is applied for the integration of the radiative transfer
equation. Otherwise the equidistant grid is used throughout. For
down-looking cases using the equidistant grid typically suffices
and speeds up the calculations.
)--",
      .author = {"Claudia Emde"},
      .out = {"doit_za_grid_size", "aa_grid", "za_grid"},

      .gin = {"N_za_grid", "N_aa_grid", "za_grid_opt_file"},
      .gin_type = {"Index", "Index", "String"},
      .gin_value = {std::nullopt, Index{1}, String("")},
      .gin_desc =
          {R"--(Number of grid points in zenith angle grid. Recommended value is 19.)--",
           R"--(Number of grid points in azimuth angle grid. Recommended value is 37.)--",
           R"--(Name of special grid for RT part.)--"},

  };

  wsm_data["DOBatchCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Performs batch calculations for radiation fields.

We perform *ybatch_n* jobs, starting at index *ybatch_start*. (Zero
based indexing, as usual.) The output arrays will have
ybatch_n elements. Indices in the output array start
with zero, independent of *ybatch_start*.

WARNING, MEMORY INTENSIVE!!!: Since the outputs of this methods can
be very large, make sure you only pass back output you need.
Estimate the size of your output by looking at the dimensions
beforehand. If you only want to pass back some fields, make sure to
empty the others at the end of your *dobatch_calc_agenda*. E.g.:
Tensor7SetConstant(cloudbox_field, 0, 0, 0, 0, 0, 0, 0, 0.)

The method performs the following:
  1. Sets *ybatch_index* = *ybatch_start*.
  2. Performs a-d until *ybatch_index* = *ybatch_start* + *ybatch_n*.

    a. Executes *dobatch_calc_agenda*.
    b. If *ybatch_index* = *ybatch_start*, resizes the output
       arrays based on *ybatch_n*.
    c. Copies calculated fields to *ybatch_index* - *ybatch_start*
       of output arrays.
    d. Adds 1 to *ybatch_index*.

Beside the *dobatch_calc_agenda*, the WSVs *ybatch_start*
and *ybatch_n* must be set before calling this method.

The input variable *ybatch_start* is set to a default of zero.
)--",
      .author = {"Oliver Lemke"},
      .out = {"dobatch_cloudbox_field",
              "dobatch_radiance_field",
              "dobatch_irradiance_field",
              "dobatch_spectral_irradiance_field"},

      .in = {"ybatch_start", "ybatch_n", "dobatch_calc_agenda"},
      .gin = {"robust"},
      .gin_type = {"Index"},
      .gin_value = {Index{0}},
      .gin_desc =
          {R"--(A flag with value 1 or 0. If set to one, the batch calculation will continue, even if individual jobs fail. In that case, a warning message is written to screen and file (out1 output stream), and the output array entry for the failed job in the output fields is left empty.)--"},
      .pass_workspace = true,

  };

  wsm_data["Delete"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Deletes a workspace variable.

The variable is not deleted from the workspace, but it is
reset to its default value. This is useful if you want to
free memory for heavy variables
)--",
      .author = {"Oliver Lemke"},

      .gout = {"v"},
      .gout_type = {"Any"},
      .gout_desc = {R"--(Variable to be deleted.)--"}
  };

  wsm_data["DiagonalMatrix"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Create a diagonal matrix from a vector.

This creates a dense or sparse diagonal matrix with the elements of the given vector
on the diagonal.
)--",
      .author = {"Simon Pfreundschuh"},

      .gout = {"output"},
      .gout_type = {"Matrix, Sparse"},
      .gout_desc = {R"--(The diagonal matrix)--"},

      .gin = {"v"},
      .gin_type = {"Vector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(The vector containing the diagonal elements.)--"},

  };

  wsm_data["DoitCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Main DOIT method.

This method executes *doit_mono_agenda* for each frequency
in *f_grid*. The output is the radiation field inside the cloudbox
(*cloudbox_field*).
)--",
      .author = {"Claudia Emde"},
      .out = {"cloudbox_field"},

      .in = {"cloudbox_field",
             "atmfields_checked",
             "atmgeom_checked",
             "cloudbox_checked",
             "scat_data_checked",
             "cloudbox_on",
             "f_grid",
             "doit_mono_agenda",
             "doit_is_initialized"},

      .pass_workspace = true,

  };

  wsm_data["DoitInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initialises variables for DOIT scattering calculations.

Note that multi-dimensional output variables (Tensors, specifically)
are NaN-initialized. That is, this methods needs to be called
BEFORE other WSMs that provide input to *DoitCalc*, e.g. before
``DoitGetIncoming``.
)--",
      .author = {"Claudia Emde"},
      .out = {"doit_scat_field", "cloudbox_field", "doit_is_initialized"},

      .in = {"f_grid",
             "za_grid",
             "aa_grid",
             "doit_za_grid_size",
             "cloudbox_on",
             "cloudbox_limits"},

  };

  wsm_data["DoitScatteringDataPrepare"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Prepares single scattering data for a DOIT scattering calculation.

First the scattering data is interpolated in frequency using
*scat_data_monoCalc*. Then the phase matrix data is
transformed or interpolated from the raw data to the laboratory frame
for all possible combinations of the angles contained in the angular
grids which are set in *DOAngularGridsSet*. The resulting phase
matrices are stored in *pha_mat_sptDOITOpt*.
)--",
      .author = {"Claudia Emde"},
      .out = {"pha_mat_sptDOITOpt",
              "scat_data_mono",
              "pha_mat_doit",
              "aa_grid"},

      .in = {"doit_za_grid_size",
             "aa_grid",
             "scat_data",
             "scat_data_checked",
             "f_index",
             "atm_field",
             "cloudbox_limits",
             "pnd_field",
             "pha_mat_spt_agenda"},

      .pass_workspace = true,

  };

  wsm_data["DoitWriteIterationFields"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Writes DOIT iteration fields.

This method writes intermediate iteration fields to xml-files. The
method can be used as a part of *doit_conv_test_agenda*.

The iterations to be stored are specified by ``iterations``, e.g. ::

  iterations = [3, 6, 9]

In this case the 3rd, 6th and 9th iterations are stored.
If a number is larger than the total number of iterations, this
number is ignored. If all iterations should be stored set::

  iterations = [-1]

The frequencies to be stored are specified by ``frequencies`` in the
same way as the iterations. The frequency index corresponds to the
order of frequencies in *f_grid*.

The output files are named doit_iteration_fX_iY.xml with X being the
frequency index and iY the iteration counter.
)--",
      .author = {"Claudia Emde"},

      .in = {"doit_iteration_counter", "cloudbox_field_mono", "f_index"},
      .gin = {"iterations", "frequencies"},
      .gin_type = {"ArrayOfIndex", "ArrayOfIndex"},
      .gin_value = {ArrayOfIndex{-1}, ArrayOfIndex{-1}},
      .gin_desc = {R"--(Selection of iterations to store.)--",
                   R"--(Selection of frequencies to store.)--"},

  };

  wsm_data["Duration"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the seconds between two times.
)--",
      .author = {"Richard Larsson"},

      .gout = {"duration"},
      .gout_type = {"Numeric"},
      .gout_desc = {R"--(Time in seconds between ``start`` and ``end``)--"},

      .gin = {"start", "end"},
      .gin_type = {"Time", "Time"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Start time)--", R"--(End time)--"},

  };

  wsm_data["Error"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Issues an error and exits ARTS.

This method can be placed in agendas that must be specified, but
are expected not to be used for the particular case. An inclusion
in *surface_rtprop_agenda* could look like::

  Error{"Surface interceptions of propagation path not expected."}

Ignore and other dummy method calls must still be included.
)--",
      .author = {"Patrick Eriksson"},

      .gin = {"msg"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(String describing the error.)--"},

  };

  wsm_data["Exit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Stops the execution and exits ARTS.

This method is handy if you want to debug one of your control
files. You can insert it anywhere in the control file. When
it is reached, it will terminate the program.
)--",
      .author = {"Patrick Eriksson"},

  };

  wsm_data["Extract"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Extracts an element from an array.

Copies the element with the given Index from the input
variable to the output variable.

For a Tensor3 as an input, it copies the page with the given
Index from the input Tensor3 variable to the output Matrix.

In other words, the selection is always done on the first dimension.
)--",
      .author = {"Oliver Lemke"},

      .gout = {"needle"},
      .gout_type =
          {"Index, ArrayOfIndex, Numeric, Vector, Matrix, Matrix, Tensor3, Tensor4, Tensor4, GriddedField2, GriddedField3, ArrayOfGriddedField3, GriddedField4, String, SingleScatteringData, ArrayOfSingleScatteringData, TelsemAtlas, QuantumIdentifier"},
      .gout_desc = {R"--(Extracted element.)--"},

      .gin = {"haystack", "index"},
      .gin_type =
          {"ArrayOfIndex, ArrayOfArrayOfIndex, Vector, ArrayOfVector, ArrayOfMatrix, Tensor3, Tensor4, ArrayOfTensor4, Tensor5, ArrayOfGriddedField2, ArrayOfGriddedField3, ArrayOfArrayOfGriddedField3, ArrayOfGriddedField4, ArrayOfString, ArrayOfSingleScatteringData, ArrayOfArrayOfSingleScatteringData, ArrayOfTelsemAtlas, ArrayOfQuantumIdentifier",
           "Index"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Variable to extract from.)--",
                   R"--(Position of the element which should be extracted.)--"},

  };

  wsm_data["ExtractFromMetaSingleScatSpecies"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Extract (numeric) parameters from scat_meta of a single scattering
species.

...
)--",
      .author = {"Jana Mendrok"},

      .gout = {"meta_param"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(The extracted meta parameter values.)--"},
      .in = {"scat_meta"},
      .gin = {"meta_name", "scat_species_index"},
      .gin_type = {"String", "Index"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(Name of the meta parameter to extract.)--",
           R"--(Array index of scattering species from which to extract.)--"},

  };

  wsm_data["FastemStandAlone"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Stand-alone usage of FASTEM.

FASTEM is a parameterisation of the emissivity of water surfaces
including the impact of waves, salinity and non-specular effects.
This is more or less direct interface to FASTEM, but slightly
adopted to fit with ARTS. The unit of frequency and salinity
differ, and this version is "vectorised" in frequency.

The output is four emissivity and reflectivity values for each
frequency. These values are defined in Eq. 13 of  "An Improved
Fast Microwave Water Emissivity Model" by Liu, Weng and English,
I3TRGS, 2011. Note that emissivity and reflectivity do not add up
to 1, which is the way FASTEM compensates for non-specular effects.

There is an error if any frequency is above 250 GHz, or if the skin
temperature is below 260 K. If the skin temperature is below 270 K,
it is adjusted to 270 K.

FASTEM returns unphysical values for propagation close to the
horizon, here emissivity and reflectivity can be outside [0,1].
If either emissivity or reflectivity is below/above 0/1, it is
set to 0/1, and the other value is set to 1/0. That is, e+r=1
is enforced. These problems start about 15 degrees from the horizon.
)--",
      .author = {"Oliver Lemke, Patrick Eriksson"},

      .gout = {"emissivity", "reflectivity"},
      .gout_type = {"Matrix", "Matrix"},
      .gout_desc =
          {R"--(Emission values. One row for each frequency. See above.)--",
           R"--(Reflectivity values. One row for each frequency. See above.)--"},
      .in = {"f_grid", "surface_skin_t"},
      .gin = {"za",
              "salinity",
              "wind_speed",
              "rel_aa",
              "transmittance",
              "fastem_version"},
      .gin_type =
          {"Numeric", "Numeric", "Numeric", "Numeric", "Vector", "Index"},
      .gin_value = {std::nullopt,
                    Numeric{0.035},
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    Index{6}},
      .gin_desc =
          {R"--(Zenith angle of line-of-sigh, 90 to 180 deg.)--",
           R"--(Salinity, 0-1. That is, 3% is given as 0.03.)--",
           R"--(Wind speed.)--",
           R"--(Azimuth angle between wind direction and line-of-sight. This angle is measured clockwise from north, i.e. E=90deg.)--",
           R"--(The transmittance of the atmosphere, along the propagation path of the downwelling radiation. One value per frequency.)--",
           R"--(The version of FASTEM to use.)--"},

  };

  wsm_data["FlagOff"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets an index variable that acts as an on/off flag to 0.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"flag"},
      .gout_type = {"Index"},
      .gout_desc = {R"--(Variable to set to 0.)--"},

  };

  wsm_data["FlagOn"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets an index variable that acts as an on/off flag to 1.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"flag"},
      .gout_type = {"Index"},
      .gout_desc = {R"--(Variable to set to 1.)--"},

  };

  wsm_data["Flatten"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Flattens an ArrayOfArray<T> to Array<T> or an Array
of matpack-types to a larger dimension matpack (if dimensions agree)

The intended transformation for arrays is (sub-arrays can have different sizes):
    {{a, b, c}, {d, e}} -> {a, b, c, d, e}

The intended transformation for arrays to matpack types is (sub-types must have same size):
    {{a, b, c}, {d, e, f}} -> {a, b, c, d, e, f}
)--",
      .author = {"Richard Larsson"},

      .gout = {"output"},
      .gout_type =
          {"ArrayOfTime, ArrayOfVector, Matrix, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7"},
      .gout_desc = {R"--(Flatter array/matpack-type)--"},

      .gin = {"input"},
      .gin_type =
          {"ArrayOfArrayOfTime, ArrayOfArrayOfVector, ArrayOfVector, ArrayOfMatrix, ArrayOfTensor3, ArrayOfTensor4, ArrayOfTensor5, ArrayOfTensor6"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(An array)--"},

  };

  wsm_data["ForLoop"] = WorkspaceMethodInternalRecord{
      .desc = R"--(A simple for-loop.

This method is handy when you quickly want to test out a calculation
with a set of different settings.

It does a for-loop from start to stop in steps of step (who would
have guessed that). For each iteration, the agenda *forloop_agenda* is
executed. Inside the agenda, the variable *forloop_index* is available
as index counter.

There are no other inputs to *forloop_agenda*, and also no outputs. That
means, if you want to get any results out of this loop, you have to
save it to files (for example with *WriteXMLIndexed*), since
variables used inside the agenda will only be local.

Note that this kind of for loop is not parallel.

The method is intended for simple testing, not as a replacement of
*ybatchCalc*. However, it is compatible with *ybatchCalc*, in the sense
that *ybatchCalc* may occur inside *forloop_agenda*.
)--",
      .author = {"Stefan Buehler"},

      .in = {"forloop_agenda"},
      .gin = {"start", "stop", "step"},
      .gin_type = {"Index", "Index", "Index"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Start value.)--",
                   R"--(End value.)--",
                   R"--(Step size.)--"},
      .pass_workspace = true,

  };

  wsm_data["FrequencyFromCGSAngularWavenumber"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Convert from angular wavenumber [cm^-1] to frequency [Hz].

This converts angular wavenumber (2*PI/wavelength) into frequency.
)--",
      .author = {"Richard Larsson"},

      .gout = {"frequency"},
      .gout_type = {"Numeric, Vector"},
      .gout_desc = {R"--(frequency [Hz])--"},

      .gin = {"angular_wavenumber"},
      .gin_type = {"Numeric, Vector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(angular wavenumber [cm^-1])--"},

  };

  wsm_data["FrequencyFromCGSKayserWavenumber"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Convert from Kayser wavenumber [cm^-1] to frequency [Hz].

This converts Kayser wavenumber (1/wavelength) into frequency.
)--",
      .author = {"Richard Larsson"},

      .gout = {"frequency"},
      .gout_type = {"Numeric, Vector"},
      .gout_desc = {R"--(frequency [Hz])--"},

      .gin = {"kayser_wavenumber"},
      .gin_type = {"Numeric, Vector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Kayser wavenumber [cm^-1])--"},

  };

  wsm_data["FrequencyFromWavelength"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Convert from wavelength [m] to frequency [Hz].

This is a generic method. It can take a single wavelength value or a wavelength vector as input.
)--",
      .author = {"Claudia Emde"},

      .gout = {"frequency"},
      .gout_type = {"Numeric, Vector"},
      .gout_desc = {R"--(frequency [Hz])--"},

      .gin = {"wavelength"},
      .gin_type = {"Numeric, Vector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(wavelength [m])--"},

  };

  wsm_data["GetEnvironmentVariable"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Copy the contents of an environment variable to an ARTS String or Index.
)--",
      .author = {"Oliver Lemke"},

      .gout = {"output"},
      .gout_type = {"String, Index"},
      .gout_desc = {R"--(Contents of environment variable.)--"},

      .gin = {"input"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Name of environment variable.)--"},

  };

  wsm_data["GetNumberOfThreads"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Returns the number of threads used by ARTS.
)--",
      .author = {"Oliver Lemke"},

      .gout = {"nthreads"},
      .gout_type = {"Index"},
      .gout_desc = {R"--(Number of threads.)--"},

  };

  wsm_data["GriddedFieldGetName"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Get the name of a GriddedField.

See *ArrayOfGriddedFieldGetNames*.
)--",
      .author = {"Lukas Kluft"},

      .gout = {"name"},
      .gout_type = {"String"},
      .gout_desc = {R"--(Name of the GriddedField.)--"},

      .gin = {"griddedfield"},
      .gin_type =
          {"GriddedField1, GriddedField2, GriddedField3, GriddedField4, GriddedField5, GriddedField6"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(GriddedField.)--"},

  };

  wsm_data["GriddedFieldLatLonExpand"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Expands the latitude and longitude grid of the GriddedField to
[-90, 90] and [0,360], respectively.

Expansion is only done in
the dimension(s), where the grid size is 1.
The values from the input data will be duplicated to accomodate
for the larger size of the output field.
output and input can be the same variable.
)--",
      .author = {"Oliver Lemke"},

      .gout = {"output"},
      .gout_type =
          {"GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3"},
      .gout_desc = {R"--(Expanded gridded field.)--"},

      .gin = {"input"},
      .gin_type =
          {"GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Raw input gridded field.)--"},

  };

  wsm_data["GriddedFieldLatLonRegrid"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Interpolates the input field along the latitude and longitude dimensions
to *lat_true* and *lon_true*.

If the input longitude grid is outside of *lon_true* it will be shifted
left or right by 360. If it covers 360 degrees, a cyclic interpolation
will be performed.
input and output fields can be the same variable.
)--",
      .author = {"Oliver Lemke"},

      .gout = {"output"},
      .gout_type =
          {"GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3"},
      .gout_desc = {R"--(Regridded gridded field.)--"},
      .in = {"lat_true", "lon_true"},
      .gin = {"input", "interp_order"},
      .gin_type =
          {"GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3",
           "Index"},
      .gin_value = {std::nullopt, Index{1}},
      .gin_desc = {R"--(Raw input gridded field.)--",
                   R"--(Interpolation order.)--"},

  };

  wsm_data["HydrotableCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Creates a look-up table of scattering properties.

The table produced largely follows the format used in RTTOV-SCATT for
its "hydrotables". The table is returned as a GriddedField4, with
dimensions (in order):

1. Scattering property
2. Frequency (equals WSV f_grid)
3. Temperature (equals GIN T_grid)
4. Particle content [kg/m3]  (equals GIN wc_grid)

Four scattering properties are calculated. They are (in order)

1. Extinction [m-1]
2. Single scattering albedo [-]
3. Asymmetry parameter [-]
4. Radar reflectivity [m2]
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"hydrotable"},
      .gout_type = {"GriddedField4"},
      .gout_desc = {R"--(Generated hydrotable with format described above.)--"},
      .in = {"pnd_agenda_array",
             "pnd_agenda_array_input_names",
             "scat_data",
             "scat_data_checked",
             "f_grid"},
      .gin = {"iss", "T_grid", "wc_grid"},
      .gin_type = {"Index", "Vector", "Vector"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Index of scattering species.)--",
                   R"--(Temperature grid of table.)--",
                   R"--(Water content grid of table.)--"},
      .pass_workspace = true,

  };

  wsm_data["INCLUDE"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes the contents of another controlfile.

The INCLUDE statement inserts the contents of the controlfile
with the given name into the current controlfile.
If the filename is given without path information, ARTS will
first search for the file in all directories specified with the
-I (see arts -h) commandline option and then in directories given
in the environment variable ARTS_INCLUDE_PATH. In the environment
variable multiple paths have to be separated by colons.

Note that INCLUDE is not a workspace method and thus the
syntax is different::

  Arts {
    INCLUDE "agendas.arts"
  }

Includes can also be nested. In the example above agendas.arts
can contain further includes which will then be treated
the same way.

The idea behind this mechanism is that you can write common settings
for a bunch of calculations into one file. Then, you can create
several controlfiles which include the basic settings and tweak them
for different cases. When you decide to make changes to your setup
that should apply to all calculations, you only have to make a
single change in the include file instead of modifying all your
controlfiles.
)--",
      .author = {"Oliver Lemke"},

  };

  wsm_data["Ignore"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Ignore a workspace variable.

This method is handy for use in agendas in order to suppress warnings
about unused input workspace variables. What it does is: Nothing!
In other words, it just ignores the variable it is called on.

This method can ignore any workspace variable you want.
)--",
      .author = {"Stefan Buehler"},

      .gin = {"input"},
      .gin_type = {"Any"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Variable to be ignored.)--"},

  };

  wsm_data["IndexAdd"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds a Index and a value (output = input + value).

The result can either be stored in the same or another Index.
)--",
      .author = {"Patrick Eriksson, Oliver Lemke"},

      .gout = {"output"},
      .gout_type = {"Index"},
      .gout_desc = {R"--(Output Index.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Index", "Index"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Index.)--", R"--(Value to add.)--"},

  };

  wsm_data["IndexDivide"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Integer division of a Index and a value (output = input / value).

Please note that integer divison is applied, and e.g. 5/3=1.

The result can either be stored in the same or another Index.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Index"},
      .gout_desc = {R"--(Output Index.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Index", "Index"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Index (numerator).)--", R"--(Denominator.)--"},

  };

  wsm_data["IndexMultiply"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Multiplies a Index and a value (output = input * value).

The result can either be stored in the same or another Index.
)--",
      .author = {"Patrick Eriksson, Oliver Lemke"},

      .gout = {"output"},
      .gout_type = {"Index"},
      .gout_desc = {R"--(Output index.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Index", "Index"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Index.)--", R"--(Multiplier.)--"},

  };

  wsm_data["IndexSetToLast"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Set an Index to point towards last position of array-type variables.

This method works as nelemGet, but gives the index number of the last
element (which equals nelem-1).
)--",
      .author = {"Patrick Eriksson", "Oliver Lemke"},
      .out = {"nelem"},

      .gin = {"v"},
      .gin_type = {"ArrayOfAbsorptionLines, ArrayOfAgenda, ArrayOfArrayOfAbsorptionLines, ArrayOfArrayOfGriddedField1, ArrayOfArrayOfGriddedField2, ArrayOfArrayOfGriddedField3, ArrayOfArrayOfIndex, ArrayOfArrayOfMatrix, ArrayOfArrayOfMuelmatMatrix, ArrayOfArrayOfMuelmatVector, ArrayOfArrayOfPropmatMatrix, ArrayOfArrayOfPropmatVector, ArrayOfArrayOfScatteringMetaData, ArrayOfArrayOfSingleScatteringData, ArrayOfArrayOfSpeciesTag, ArrayOfArrayOfStokvecMatrix, ArrayOfArrayOfStokvecVector, ArrayOfArrayOfString, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6, ArrayOfArrayOfTime, ArrayOfArrayOfVector, ArrayOfAtmPoint, ArrayOfCIARecord, ArrayOfGriddedField1, ArrayOfGriddedField2, ArrayOfGriddedField3, ArrayOfGriddedField4, ArrayOfIndex, ArrayOfJacobianTarget, ArrayOfMatrix, ArrayOfMuelmatMatrix, ArrayOfMuelmatVector, ArrayOfPpath, ArrayOfPropmatMatrix, ArrayOfPropmatVector, ArrayOfQuantumIdentifier, ArrayOfRetrievalQuantity, ArrayOfScatteringMetaData, ArrayOfSingleScatteringData, ArrayOfSparse, ArrayOfSpeciesTag, ArrayOfStokvecMatrix, ArrayOfStokvecVector, ArrayOfString, ArrayOfSun, ArrayOfTelsemAtlas, ArrayOfTensor3, ArrayOfTensor4, ArrayOfTensor5, ArrayOfTensor6, ArrayOfTensor7, ArrayOfTime, ArrayOfVector, ArrayOfXsecRecord, Vector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(The method is defined for these groups.)--"},

  };

  wsm_data["IndexStepDown"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Performas: output = input - 1

Input and output can be same variable.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Index"},
      .gout_desc = {R"--(Output index variable.)--"},

      .gin = {"input"},
      .gin_type = {"Index"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Input index variable.)--"},

  };

  wsm_data["IndexStepUp"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Performas: output = input + 1

Input and output can be same variable.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Index"},
      .gout_desc = {R"--(Output index variable.)--"},

      .gin = {"input"},
      .gin_type = {"Index"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Input index variable.)--"},

  };

  wsm_data["IndexSubtract"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Subtracts a Index value (output = input - value).

The result can either be stored in the same or another Index.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Index"},
      .gout_desc = {R"--(Output Index.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Index", "Index"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Index.)--", R"--(Subtrahend.)--"},

  };

  wsm_data["InterpAtmFieldToPosition"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Point interpolation of atmospheric fields.

The default way to specify the position is by *rtp_pos*.

Linear interpolation is applied.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"out"},
      .gout_type = {"AtmPoint"},
      .gout_desc = {R"--(Value obtained by the interpolation.)--"},
      .in = {"atm_field", "rtp_pos"},

  };

  wsm_data["InterpGriddedField2ToPosition"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Latitude and longitude interpolation of a GriddedField2.

The default way to specify the position is by *rtp_pos*.

The interpolation is done for the latitude and longitude in
*rtp_pos*. The altitude in *rtp_pos* is completely ignored.
Linear interpolation is applied.

The input field (``gfield2``) is expected to have latitude and
longitude as first and second dimension.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Numeric"},
      .gout_desc = {R"--(Value obtained by interpolation.)--"},
      .in = {"lat_true", "lon_true", "rtp_pos"},
      .gin = {"gfield2"},
      .gin_type = {"GriddedField2"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Gridded field to interpolate.)--"},

  };

  wsm_data["InterpSurfaceFieldToPosition"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Point interpolation of surface fields.

The default way to specify the position is by *rtp_pos*.

Linear interpolation is applied.

The interpolation is done for the latitude and longitude in
*rtp_pos*, while the altitude in *rtp_pos* is not part of the
calculations. However, it is checked that the altitude of *rtp_pos*
is inside the range covered by ``z_surface`` with a 1 m margin, to
give a warning when the specified position is not consistent with
the surface altitudes.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_point"},

      .in = {"rtp_pos", "surface_field", "surface_search_accuracy"},

  };

  wsm_data["IntersectionGeometricAltitude"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates the geometrical intersection with an altitude.

For each observation geometry specified by the combination of
*sensor_pos* and *sensor_los*, the geometrical intersection with
an altitude is determined. The intersections are described by the
GOUT ``pos`` and ``los``.

For cases with no intersection, ``pos`` and ``los`` are filled with NaN.

The GOUT ``pos`` and ``los`` can NOT be *sensor_pos* and *sensor_los*.
If you want to store the intersections in *sensor_pos* and *sensor_los*
use *sensor_pos_losForwardToAltitude*. For *rte_pos* and *rte_los*
you have *rte_pos_losForwardToAltitude*.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"pos", "los"},
      .gout_type = {"Matrix", "Matrix"},
      .gout_desc = {R"--(Position of intersections.)--",
                    R"--(Line-of-sight at intersections.)--"},
      .in = {"sensor_pos", "sensor_los", "surface_field"},
      .gin = {"altitude"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{0}},
      .gin_desc = {R"--(Target altitude.)--"},

  };

  wsm_data["IntersectionGeometricLatitude"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates the geometrical intersection with a latitude.

For each observation geometry specified by the combination of
*sensor_pos* and *sensor_los*, the geometrical intersection with
a latitude is determined. The intersections are described by the
GOUT ``pos`` and ``los``.

For cases with no intersection, ``pos`` and ``los`` are filled with NaN.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"pos", "los"},
      .gout_type = {"Matrix", "Matrix"},
      .gout_desc = {R"--(Position of intersections.)--",
                    R"--(Line-of-sight at intersections.)--"},
      .in = {"sensor_pos", "sensor_los", "surface_field"},
      .gin = {"latitude"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Target latitude.)--"},

  };

  wsm_data["IntersectionGeometricLongitude"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates the geometrical intersection with a longitude.

For each observation geometry specified by the combination of
*sensor_pos* and *sensor_los*, the geometrical intersection with
a longitude is determined. The intersections are described by the
GOUT ``pos`` and *los.

For cases with no intersection, ``pos`` and ``los`` are filled with NaN.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"pos", "los"},
      .gout_type = {"Matrix", "Matrix"},
      .gout_desc = {R"--(Position of intersections.)--",
                    R"--(Line-of-sight at intersections.)--"},
      .in = {"sensor_pos", "sensor_los", "surface_field"},
      .gin = {"longitude"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Target longitude.)--"},

  };

  wsm_data["IntersectionGeometricSurface"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates the geometrical intersection with the surface.

For each observation geometry specified by the combination of
*sensor_pos* and *sensor_los*, the geometrical intersection with
the surface is determined. The intersections are described by the
GOUT ``pos`` and ``los``. For cases with no intersection, ``pos`` and ``los``
are filled with NaN.

If the surface elevation is constant, the intersections are found
analytically. Otherwise a search in terms of distance from the sensor
is applied. The default is to use a bisection algorithm. This option
should suffice in general, but it can fail if the elevation varies
strongly and/or the incidence angle is high. The path can then cross
the surface at several positions and the bisection search does not
guarantee that the correct intersection is found. To avoid this, set
*surface_search_safe* to 1 and a safe, but much more slow option, is
used. In this case the path is simple sampled in a step-by-step
fashion.

To be clear, the faster bisection algorith can fail if the path goes
through a mountain top. For an upward observation inside a valley, the
bisection can also miss if the path touches the side of the valley.

For both algorithms *surface_search_accuracy* governs the accuracy. In
the first case, the bisection is stopped when the set accuracy has
been reached, while in the safe option *surface_search_accuracy* is
the step length used.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"pos", "los"},
      .gout_type = {"Matrix", "Matrix"},
      .gout_desc = {R"--(Position of intersections.)--",
                    R"--(Line-of-sight at intersections.)--"},
      .in = {"sensor_pos",
             "sensor_los",
             "surface_field",
             "surface_search_accuracy",
             "surface_search_safe"},

  };

  wsm_data["LatLonFieldSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Fills a latitude-longitude field with given input.

Grids and data must match in size.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"gfield2"},
      .gout_type = {"GriddedField2"},
      .gout_desc = {R"--(Field to set.)--"},

      .gin = {"latitude_grid", "longitude_grid", "data", "name"},
      .gin_type = {"Vector", "Vector", "Matrix", "String"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, String("")},
      .gin_desc = {R"--(The latitude grid of ``data``.)--",
                   R"--(The longitude grid of ``data``.)--",
                   R"--(The data of the field (will become gfield2.data).)--",
                   R"--(The name of the field (will become gfield2.name).)--"},

  };

  wsm_data["LatLonFieldSetConstant"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets a latitude-longitude field to have a constant data value.

Both latitude and longitude grids are set to have length one,
with the value 0.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"gfield2"},
      .gout_type = {"GriddedField2"},
      .gout_desc = {R"--(Field to set.)--"},

      .gin = {"value", "name"},
      .gin_type = {"Numeric", "String"},
      .gin_value = {std::nullopt, String("")},
      .gin_desc = {R"--(The value (to place in gfield2.data).)--",
                   R"--(The name of the field (will become gfield2.name).)--"},

  };

  wsm_data["LocalTimeOffset"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets the seconds between localtime and gmtime representation of now().
)--",
      .author = {"Richard Larsson"},

      .gout = {"dt"},
      .gout_type = {"Numeric"},
      .gout_desc = {R"--(Time in seconds between local and gmt)--"},

  };

  wsm_data["MCGeneral"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(A generalised 3D reversed Monte Carlo radiative algorithm, that
allows for 2D antenna patterns, surface reflection and arbitrary
sensor positions.

The main output variables *y* and *mc_error* represent the
Stokes vector integrated over the antenna function, and the
estimated error in this vector, respectively.

The WSV *mc_max_iter* describes the maximum number of 'photons'
used in the simulation (more photons means smaller *mc_error*).
*mc_std_err* is the desired value of mc_error. *mc_max_time* is
the maximum allowed number of seconds for MCGeneral. The method
will terminate once any of the max_iter, std_err, max_time
criteria are met. If negative values are given for these
parameters then it is ignored.

The WSV *mc_min_iter* sets the minimum number of photons to apply
before the condition set by *mc_std_err* is considered. Values
of *mc_min_iter* below 100 are not accepted.

Only "1" and "RJBT" are allowed for *iy_unit*. The value of
*mc_error* follows the selection for *iy_unit* (both for in- and
output.
)--",
      .author = {"Cory Davis"},
      .out = {"y",
              "mc_iteration_count",
              "mc_error",
              "mc_points",
              "mc_source_domain",
              "mc_scat_order"},

      .in = {"mc_antenna",
             "f_grid",
             "f_index",
             "sensor_pos",
             "sensor_los",
             "ppath_step_agenda",
             "ppath_lmax",
             "ppath_lraytrace",
             "iy_space_agenda",
             "surface_rtprop_agenda",
             "propmat_clearsky_agenda",
             "surface_field",
             "atm_field",
             "cloudbox_on",
             "cloudbox_limits",
             "pnd_field",
             "scat_data",
             "atmfields_checked",
             "atmgeom_checked",
             "scat_data_checked",
             "cloudbox_checked",
             "iy_unit",
             "mc_seed",
             "mc_std_err",
             "mc_max_time",
             "mc_max_iter",
             "mc_min_iter",
             "mc_taustep_limit"},
      .gin = {"l_mc_scat_order", "t_interp_order"},
      .gin_type = {"Index", "Index"},
      .gin_value = {Index{11}, Index{1}},
      .gin_desc =
          {R"--(The length to be given to *mc_scat_order*. Note that scattering orders equal and above this value will not be counted.)--",
           R"--(Interpolation order of temperature for scattering data (so far only applied in phase matrix, not in extinction and absorption.)--"},
      .pass_workspace = true,

  };

  wsm_data["MCRadar"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(A radar 3D foward Monte Carlo radiative algorithm, that allows 
for 2D antenna patterns and arbitrary sensor positions.
Surface reflections are currently ignored.

The main output variable *y* and *mc_error* represent the
radar reflectivity integrated over the antenna function, and the
estimated error in this vector, respectively.

Unlike with yRadar, the range bins gives the boundaries of 
the range bins as either round-trip time or distance from radar.

The WSV *mc_y_tx* gives the polarization state of the 
transmitter.

The WSV *mc_max_scatorder* prescribes the maximum scattering 
order to consider, after which 'photon'-tracing will be
terminated. A value of one calculates only single scattering.

The WSV *mc_max_iter* describes the maximum number of 'photons'
used in the simulation (more photons means smaller *mc_error* ).
The method will terminate once the max_iter criterium is met.
If negative values are given for these parameters then it is
ignored.

Here "1" and "Ze" are the allowed options for *iy_unit_radar*.
The value of *mc_error* follows the selection for *iy_unit_radar*
(both for in- and output. See *yRadar* for details of the units.
)--",
      .author = {"Ian S. Adams"},
      .out = {"y", "mc_error"},

      .in = {"mc_antenna",
             "f_grid",
             "f_index",
             "sensor_pos",
             "sensor_los",
             "ppath_lmax",
             "ppath_step_agenda",
             "ppath_lraytrace",
             "propmat_clearsky_agenda",
             "surface_field",
             "atm_field",
             "cloudbox_on",
             "cloudbox_limits",
             "pnd_field",
             "scat_data",
             "mc_y_tx",
             "range_bins",
             "atmfields_checked",
             "atmgeom_checked",
             "scat_data_checked",
             "cloudbox_checked",
             "iy_unit_radar",
             "mc_max_scatorder",
             "mc_seed",
             "mc_max_iter"},
      .gin = {"ze_tref", "k2", "t_interp_order"},
      .gin_type = {"Numeric", "Numeric", "Index"},
      .gin_value = {Numeric{273.15}, Numeric{-1}, Index{1}},
      .gin_desc =
          {R"--(Reference temperature for conversion to Ze.)--",
           R"--(Reference dielectric factor.)--",
           R"--(Interpolation order of temperature for scattering data (so far only applied in phase matrix, not in extinction and absorption.)--"},
      .pass_workspace = true,

  };

  wsm_data["MCSetSeedFromTime"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the value of mc_seed from system time
)--",
      .author = {"Cory Davis"},
      .out = {"mc_seed"},

  };

  wsm_data["Matrix1ColFromVector"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Forms a matrix containing one column from a vector.
)--",
      .author = {"Mattias Ekstrom"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Variable to initialize.)--"},

      .gin = {"v"},
      .gin_type = {"Vector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(The vector to be copied.)--"},

  };

  wsm_data["Matrix1RowFromVector"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Forms a matrix containing one row from a vector.
)--",
      .author = {"Mattias Ekstrom"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Variable to initialize.)--"},

      .gin = {"v"},
      .gin_type = {"Vector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(The vector to be copied.)--"},

  };

  wsm_data["Matrix2ColFromVectors"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Forms a matrix containing two columns from two vectors.

The vectors are included as columns in the matrix in the same order
as they are given.
)--",
      .author = {"Mattias Ekstrom"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Variable to initialize.)--"},

      .gin = {"v1", "v2"},
      .gin_type = {"Vector", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(The vector to be copied into the first column.)--",
                   R"--(The vector to be copied into the second column.)--"},

  };

  wsm_data["Matrix2RowFromVectors"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Forms a matrix containing two rows from two vectors.

The vectors are included as rows in the matrix in the same order
as they are given.
)--",
      .author = {"Mattias Ekstrom"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Variable to initialize.)--"},

      .gin = {"v1", "v2"},
      .gin_type = {"Vector", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(The vector to be copied into the first row.)--",
                   R"--(The vector to be copied into the second row.)--"},

  };

  wsm_data["Matrix3ColFromVectors"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Forms a matrix containing three columns from three vectors.

The vectors are included as columns in the matrix in the same order
as they are given.
)--",
      .author = {"Mattias Ekstrom"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Variable to initialize.)--"},

      .gin = {"v1", "v2", "v3"},
      .gin_type = {"Vector", "Vector", "Vector"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(The vector to be copied into the first column.)--",
                   R"--(The vector to be copied into the second column.)--",
                   R"--(The vector to be copied into the third column.)--"},

  };

  wsm_data["Matrix3RowFromVectors"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Forms a matrix containing three rows from three vectors.

The vectors are included as rows in the matrix in the same order
as they are given.
)--",
      .author = {"Mattias Ekstrom"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Variable to initialize.)--"},

      .gin = {"v1", "v2", "v3"},
      .gin_type = {"Vector", "Vector", "Vector"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(The vector to be copied into the first row.)--",
                   R"--(The vector to be copied into the second row.)--",
                   R"--(The vector to be copied into the third row.)--"},

  };

  wsm_data["MatrixAdd"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds a scalar to all elements of a matrix.

The result can either be stored in the same or another matrix.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Output Matrix)--"},

      .gin = {"input", "value"},
      .gin_type = {"Matrix", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Matrix.)--",
                   R"--(The value to be added to the matrix.)--"},

  };

  wsm_data["MatrixCBR"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets a matrix to hold cosmic background radiation (CBR).

The CBR is assumed to be un-polarized and Stokes components 2-4
are zero. Number of Stokes components, that equals the number
of columns in the created matrix, is determined by ``stokes_dim``.
The number of rows in the created matrix equals the length of the
given frequency vector.

The cosmic radiation is modelled as blackbody radiation for the
temperature given by the global constant COSMIC_BG_TEMP, set in
the file constants.cc. The frequencies are taken from the generic
input vector.

The standard definition, in ARTS, of the Planck function is
followed and the unit of the returned data is W/(m3 * Hz * sr).
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Variable to initialize.)--"},

      .gin = {"f"},
      .gin_type = {"Vector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Frequency vector.)--"},

  };

  wsm_data["MatrixCopySparse"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Creates a matrix by copying a variable of type Sparse.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Created (full) matrix.)--"},

      .gin = {"input"},
      .gin_type = {"Sparse"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(The sparse matrix to be copied.)--"},

  };

  wsm_data["MatrixDivide"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Divides all elements of a matrix with the specified value.

The result can either be stored in the same or another
variable.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Output Matrix)--"},

      .gin = {"input", "value"},
      .gin_type = {"Matrix", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Matrix.)--", R"--(Denominator.)--"},

  };

  wsm_data["MatrixExtractFromTensor3"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Extracts a Matrix from a Tensor3.

Copies page or row or column with given Index from input Tensor3
variable to output Matrix.
Higher order equivalent of *VectorExtractFromMatrix*.
)--",
      .author = {"Jana Mendrok"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Extracted matrix.)--"},

      .gin = {"input", "i", "direction"},
      .gin_type = {"Tensor3", "Index", "String"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input matrix.)--",
                   R"--(Index of page or row or column to extract.)--",
                   R"--(Direction. "page" or "row" or "column".)--"},

  };

  wsm_data["MatrixFromCovarianceMatrix"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Turns a covariance matrix into a Matrix.
)--",
      .author = {"Richard Larsson"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Dense Matrix.)--"},

      .gin = {"input"},
      .gin_type = {"CovarianceMatrix"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Input covariance matrix.)--"},

  };

  wsm_data["MatrixGaussian"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Fills a matrix with a Gaussian function.

Works as *VectorGaussian* but grid, mean and si/fwhm must be
specified for each dimension.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"Y"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Output Matrix.)--"},

      .gin = {"x_row",
              "x0_row",
              "si_row",
              "fwhm_row",
              "x_col",
              "x0_col",
              "si_col",
              "fwhm_col"},
      .gin_type = {"Vector",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Vector",
                   "Numeric",
                   "Numeric",
                   "Numeric"},
      .gin_value = {std::nullopt,
                    Numeric{0},
                    Numeric{-1},
                    Numeric{-1},
                    std::nullopt,
                    Numeric{0},
                    Numeric{-1},
                    Numeric{-1}},
      .gin_desc =
          {R"--(Grid of the function for row dimension.)--",
           R"--(Centre/mean point of the function for row dimension.)--",
           R"--(Row standard deviation of the function, ignored if <=0.)--",
           R"--(Row full width at half-max of the function, ignored if <=0.)--",
           R"--(Grid of the function for column dimension.)--",
           R"--(Centre/mean point of the function for column dimension.)--",
           R"--(Column standard deviation of the function, ignored if <=0.)--",
           R"--(Column full width at half-max of the function, ignored if <=0.)--"},

  };

  wsm_data["MatrixIdentity"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Returns the identity matrix.

The size if the matrix created is n x n. Default is to return a
true identity matrix (I), but you can also select another value
along the diagonal by setting ``value``. That is, the output is
value * I.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Output matrix)--"},

      .gin = {"n", "value"},
      .gin_type = {"Index", "Numeric"},
      .gin_value = {std::nullopt, Numeric{1}},
      .gin_desc = {R"--(Size of the matrix)--",
                   R"--(The value along the diagonal.)--"},

  };

  wsm_data["MatrixMatrixMultiply"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Multiply a Matrix with another Matrix and store the result in the result
Matrix.

This just computes the normal Matrix-Matrix product, Y = M * X. It is ok
if Y and X are the same Matrix.
)--",
      .author = {"Stefan Buehler"},

      .gout = {"Y"},
      .gout_type = {"Matrix"},
      .gout_desc =
          {R"--(The result of the multiplication (dimension m x c).)--"},

      .gin = {"M", "X"},
      .gin_type = {"Matrix", "Matrix"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(The Matrix to multiply (dimension m x n).)--",
                   R"--(The original Matrix (dimension n x c).)--"},

  };

  wsm_data["MatrixMultiply"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Multiplies all elements of a matrix with the specified value.

The result can either be stored in the same or another
variable.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Output Matrix)--"},

      .gin = {"input", "value"},
      .gin_type = {"Matrix", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Matrix.)--",
                   R"--(The value to be multiplied with the matrix.)--"},

  };

  wsm_data["MatrixPlanck"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets a matrix to hold blackbody radiation.

The radiation is assumed to be un-polarized and Stokes components
2-4 are zero. Number of Stokes components, that equals the number
of columns in the created matrix, is determined by ``stokes_dim``.
The number of rows in the created matrix equals the length of the
given frequency vector.

The standard definition, in ARTS, of the Planck function is
followed and the unit of the returned data is W/(m3 * Hz * sr).
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Variable to initialize.)--"},

      .gin = {"f", "t"},
      .gin_type = {"Vector", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Frequency vector.)--", R"--(Temperature [K].)--"},

  };

  wsm_data["MatrixReshapeTensor3"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Creates a matrix as reshaped version of a tenor3.

If the size of the tensor is [npages, nrows, ncols], the created
matrix gets size [npages * nrows, ncols]. The matrix is filled with
the tensor's page dimension as the outermost loop.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Matrix to fill.)--"},

      .gin = {"input"},
      .gin_type = {"Tensor3"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Tensor3 to copy.)--"},

  };

  wsm_data["MatrixSetConstant"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Creates a matrix and sets all elements to the specified value.

The size is determined by *ncols* and *nrows*.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Variable to initialize.)--"},
      .in = {"nrows", "ncols"},
      .gin = {"value"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Matrix value.)--"},

  };

  wsm_data["MatrixSubtract"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Subtracts a scalar from all elements of a matrix.

The result can either be stored in the same or another matrix.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Output Matrix)--"},

      .gin = {"input", "value"},
      .gin_type = {"Matrix", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Matrix.)--",
                   R"--(The value to be subtracted from the matrix.)--"},

  };

  wsm_data["MatrixUnitIntensity"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets a matrix to hold unpolarised radiation with unit intensity.

Works as MatrixPlanck where the radiation is set to 1.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Variable to initialize.)--"},

      .gin = {"f"},
      .gin_type = {"Vector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Frequency vector.)--"},

  };

  wsm_data["NumericAdd"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds a Numeric and a value (output = input + value).

The result can either be stored in the same or another Numeric.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Numeric"},
      .gout_desc = {R"--(Output Numeric.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Numeric.)--", R"--(Value to add.)--"},

  };

  wsm_data["NumericClip"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Clipping of a Numeric.

The input value is copied to the output one (that can be same WSV)
but ensures that ``out`` is inside the range [limit_low,limit_high].
When the input value is below ``limit_low``, ``output`` is set to ``limit_low``.
And the same is performed with respect to ``limit_high``.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Numeric"},
      .gout_desc = {R"--(Output Numeric.)--"},

      .gin = {"input", "limit_low", "limit_high"},
      .gin_type = {"Numeric", "Numeric", "Numeric"},
      .gin_value = {std::nullopt, -std::numeric_limits<Numeric>::infinity(), std::numeric_limits<Numeric>::infinity()},
      .gin_desc = {R"--(Input Numeric.)--",
                   R"--(Lower limit for clipping.)--",
                   R"--(Upper limit for clipping.)--"},

  };

  wsm_data["NumericDivide"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Divides a Numeric with a value (output = input / value).

The result can either be stored in the same or another Numeric.
)--",
      .author = {"Jana Mendrok"},

      .gout = {"output"},
      .gout_type = {"Numeric"},
      .gout_desc = {R"--(Output Numeric.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Numeric (numerator).)--", R"--(Denominator.)--"},

  };

  wsm_data["NumericFromVector"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Derivs a Numeric from a vector, following selected operation.

The following operations can be selected:

- ``"first"``: Selects the first element of the vector.
- ``"last"``: Selects the last element of the vector.
- ``"max"``: Selects the maximum element of the vector.
- ``"min"``: Selects the minimum element of the vector.
- ``"mean"``: Calculates the mean of the vector.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Numeric"},
      .gout_desc = {R"--(Output Numeric.)--"},

      .gin = {"input", "op"},
      .gin_type = {"Vector", "String"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input vector.)--", R"--(Selected operation.)--"},

  };

  wsm_data["NumericInterpAltLatLonField"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Interpolates an altitude-latitude-longitiude field.

The gridded field must have "Altitude", "Latitude" and
"Longitude" as dimensions.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"value"},
      .gout_type = {"Numeric"},
      .gout_desc = {R"--(Result of interpolation)--"},

      .gin = {"gfield3", "pos"},
      .gin_type = {"GriddedField3", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Gridded field to be interpolated.)--",
                   R"--(Interpolate to this position [z, lat, lon].)--"},

  };

  wsm_data["NumericInterpLatLonField"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Interpolates a latitude-longitiude field to the selected position.

The gridded field must have "Latitude" and "Longitude" as dimensions.

The position shall be given as a full atmospheric position. The altitude
in ``pos`` is ignored.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"value"},
      .gout_type = {"Numeric"},
      .gout_desc = {R"--(Result of interpolation)--"},

      .gin = {"gfield2", "pos"},
      .gin_type = {"GriddedField2", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Gridded field to be interpolated.)--",
                   R"--(Interpolate to this position [z, lat, lon].)--"},

  };

  wsm_data["NumericInterpVector"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Interpolates a vector to the selected position.

Returns y(xv) given y(x).
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"value"},
      .gout_type = {"Numeric"},
      .gout_desc = {R"--(Result of interpolation)--"},

      .gin = {"gx", "gy", "xv"},
      .gin_type = {"Vector", "Vector", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Positions (grid) where *y* given.)--",
                   R"--(Values of function to interpolate.)--",
                   R"--(Interpolate to this value.)--"},

  };

  wsm_data["NumericMultiply"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Multiplies a Numeric with a value (output = input * value).

The result can either be stored in the same or another Numeric.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Numeric"},
      .gout_desc = {R"--(Output Numeric.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Numeric.)--", R"--(Multiplier.)--"},

  };

  wsm_data["NumericSubtract"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Subtracts a Numeric value (output = input - value).

The result can either be stored in the same or another Numeric.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Numeric"},
      .gout_desc = {R"--(Output Numeric.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Numeric.)--", R"--(Subtrahend.)--"},

  };

  wsm_data["OEM"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Inversion by the so called optimal estimation method (OEM).

Work in progress ...

The cost function to minimise, including a normalisation with lengthof *y*, is::

  cost = cost_y + cost_x

where::

  cost_y = 1/m * [y-yf]' * covmat_se_inv * [y-yf]
  cost_x = 1/m * [x-xa]' * covmat_sx_inv * [x-xa]

The current implementation provides 3 methods for the minimization of
the cost functional: Linear, Gauss-Newton and Levenberg-Marquardt.
The Gauss-Newton minimizer attempts to find a minimum solution by 
fitting a quadratic function to the cost functional. The linear minimizer
is a special case of the Gauss-Newton method, since for a linear forward
model the exact solution of the minimization problem is obtained after
the first step. The Levenberg-Marquardt method adaptively constrains the
search region for the next iteration step by means of the so-called gamma-factor.
This makes the method more suitable for strongly non-linear problems.
If the gamma-factor is 0, Levenberg-Marquardt and Gauss-Newton method
are identical. Each minimization method (li,gn,lm) has an indirect
variant (li_cg,gn_cg,lm_cg), which uses the conjugate gradient solver
for the linear system that has to be solved in each minimzation step.
This of advantage for very large problems, that would otherwise require
the computation of expensive matrix products.

Description of the special input arguments:

- ``method``: One of the following:

  - ``"li"``: A linear problem is assumed and a single iteration is performed.
  - ``"li_cg"``: A linear problem is assumed and solved using the CG solver.
  - ``"gn"``: Non-linear, with Gauss-Newton iteration scheme.
  - ``"gn_cg"``: Non-linear, with Gauss-Newton and conjugate gradient solver.
  - ``"lm"``: Non-linear, with Levenberg-Marquardt (LM) iteration scheme.
  - ``"lm_cg"``: Non-linear, with Levenberg-Marquardt (LM) iteration scheme and conjugate gradient solver.

- ``max_start_cost``:
  No inversion is done if the cost matching the a priori state is above
  this value. If set to a negative value, all values are accepted.
  This argument also controls if the start cost is calculated. If
  set to <= 0, the start cost in *oem_diagnostics* is set to NaN
  when using "li" and "gn".
- ``x_norm``:
  A normalisation vector for *x*. A normalisation of *x* can be needed
  due to limited numerical precision. If this vector is set to be empty
  no normalisation is done (defualt case). Otherwise, this must be a
  vector with same length as *x*, just having values above zero.
  Elementwise division between *x* and ``x_norm`` (x./x_norm) shall give
  a vector where all values are in the order of unity. Maybe the best
  way to set ``x_norm`` is x_norm = sqrt( diag( Sx ) ).
- ``max_iter``:
  Maximum number of iterations to perform. No effect for "li".
- ``stop_dx``:
  Iteration stop criterion. The criterion used is the same as given
  in Rodgers' "Inverse Methods for Atmospheric Sounding"
- ``lm_ga_settings``:
  Settings controlling the gamma factor, part of the "LM" method.
  This is a vector of length 6, having the elements (0-based index):

    0. Start value.
    1. Fractional decrease after succesfull iteration.
    2. Fractional increase after unsuccessful iteration.
    3. Maximum allowed value. If the value is passed, the inversion
       is halted.
    4. Lower treshold. If the threshold is passed, gamma is set to zero.
       If gamma must be increased from zero, gamma is set to this value.
    5. Gamma limit. This is an additional stop criterion. Convergence
       is not considered until there has been one succesful iteration
       having a gamma <= this value.

  The default setting triggers an error if "lm" is selected.
- ``clear matrices``:
   With this flag set to 1, *jacobian* and *dxdy* are returned as empty
   matrices.
- ``display_progress``:
   Controls if there is any screen output. The overall report level
   is ignored by this WSM.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"x",
              "yf",
              "jacobian",
              "dxdy",
              "oem_diagnostics",
              "lm_ga_history",
              "oem_errors"},

      .in = {"xa",
             "x",
             "covmat_sx",
             "yf",
             "y",
             "covmat_se",
             "jacobian",
             "jacobian_quantities",
             "inversion_iterate_agenda"},
      .gin = {"method",
              "max_start_cost",
              "x_norm",
              "max_iter",
              "stop_dx",
              "lm_ga_settings",
              "clear_matrices",
              "display_progress"},
      .gin_type = {"String",
                   "Numeric",
                   "Vector",
                   "Index",
                   "Numeric",
                   "Vector",
                   "Index",
                   "Index"},
      .gin_value = {std::nullopt,
                    std::numeric_limits<Numeric>::infinity(),
                    Vector{},
                    Index{10},
                    Numeric{0.01},
                    Vector{},
                    Index{0},
                    Index{0}},
      .gin_desc =
          {R"--(Iteration method. For this and all options below, see further above.)--",
           R"--(Maximum allowed value of cost function at start.)--",
           R"--(Normalisation of Sx.)--",
           R"--(Maximum number of iterations.)--",
           R"--(Stop criterion for iterative inversions.)--",
           R"--(Settings associated with the ga factor of the LM method.)--",
           R"--(An option to save memory.)--",
           R"--(Flag to control if inversion diagnostics shall be printed on the screen.)--"},
      .pass_workspace = true,

  };

  wsm_data["PlanetSet"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets *g0_agenda*, ``refellipsoid``, *molarmass_dry_air*, and *planet_rotation_period* to default values

*g0_agenda* is set using *g0_agendaSet* with the same option

Note that the default value of *isotopologue_ratios* is set to "Earth" by
default and that we strongly recommend users to update these values if they
are using non-Earth atmospheres.

Options are:

- ``"Earth"``:

    1. Uses ``refellipsoidEarth`` with model="Sphere"
    2. Sets *molarmass_dry_air* to 28.966
    3. Sets *planet_rotation_period* to 86164.1

- ``"Io"``:

    1. Uses ``refellipsoidIo`` with model="Sphere"
    2. Sets *molarmass_dry_air* to 63.110068828000003
    3. Sets *planet_rotation_period* to 152853

- ``"Jupiter"``:

    1. Uses ``refellipsoidJupiter`` with model="Sphere"
    2. Sets *molarmass_dry_air* to 2.22
    3. Sets *planet_rotation_period* to 35730

- ``"Mars"``:

    1. Uses ``refellipsoidMars`` with model="Sphere"
    2. Sets *molarmass_dry_air* to 43.34
    3. Sets *planet_rotation_period* to 88643

- ``"Venus"``:

    1. Uses ``refellipsoidVenus`` with model="Sphere"
    2. Sets *molarmass_dry_air* to 43.45
    3. Sets *planet_rotation_period* to -2.0997e7
)--",
      .author = {"Richard Larsson"},
      .out = {"g0_agenda",
              "surface_field",
              "molarmass_dry_air",
              "planet_rotation_period"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["Print"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Prints a variable on the screen.
)--",
      .author = {"Oliver Lemke"},

      .gin = {"input", "level"},
      .gin_type = {"Any", "Index"},
      .gin_value = {std::nullopt, Index{1}},
      .gin_desc = {R"--(Variable to be printed.)--",
                   R"--(Output level to use.)--"},

  };

  wsm_data["PrintPhysicalConstants"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Prints (most) physical constants used in ARTS.
)--",
      .author = {"Richard Larsson"},

  };

  wsm_data["PrintWorkspace"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Prints a list of the workspace variables.
)--",
      .author = {"Oliver Lemke"},

      .gin = {"only_allocated", "level"},
      .gin_type = {"Index", "Index"},
      .gin_value = {Index{1}, Index{1}},
      .gin_desc =
          {R"--(Flag for printing either all variables (0) or only allocated ones (1).)--",
           R"--(Output level to use.)--"},
      .pass_workspace = true,

  };

  wsm_data["RT4Calc"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Interface to the PolRadTran RT4 scattering solver (by F. Evans).

RT4 provides the radiation field (*cloudbox_field*) from a vector
1D scattering solution assuming a plane-parallel atmosphere (flat
Earth). It calculates up to two Stokes parameters (``stokes_dim`` <= 2),
i.e., all azimuthally randomly oriented particles are allowed (this
also includes macroscopically isotropic particles). Refraction is
not taken into account.

The scattering solution is internally obtained over the full
(plane-parallel) atmosphere, i.e. not confined to the cloudbox.
However, the radiation field output is limited to the cloudbox.
This allows to consider clearsky RT through a non-spherical
atmosphere outside the cloudbox improving the RT solution for
non-plane-parallel media compared to the plain RT4 output.

``nstreams`` is the number of polar angles taken into account
internally in the scattering solution. That is, ``nstreams``
determines the angular resolution, hence the accuracy, of the
scattering solution. The more anisotropic the bulk scattering
matrix, the more streams are required. The computational burden
increases approximately with the third power of ``nstreams``.
The default value (``nstreams`` = 16) was found to be sufficient for
most microwave scattering calculations. It is likely insufficient
for IR calculations involving ice clouds, though.

Here, *za_grid* is NOT an input parameter, but output, and its
size equals ``nstreams`` or ``nstreams`` + 2 (Gauss-Legendre and Double
Gauss quadratures in case ``add_straight_angles`` = 1) (the reason is
that the computational burden is high for additional angles,
regardless whether they are quadrature angles or not; hence the
quadrature angles supplemented with 0 and 180deg are considered to
provide the best radiation field for a given effort).

The ``auto_inc_nstreams`` feature can be used to increase the number
of streams used internally in the scattering solution when found
necessary.
NOTE: this number-of-streams increase is only internally - the
angular dimension of the output *cloudbox_field* is fixed to the
``nstreams`` given as input to this WSM.

Quadrature methods available are: 'L'obatto, 'G'auss-Legendre and
'D'ouble Gauss quadrature.

This WSM applies *surface_rtprop_agenda* to derive reflection
matrix and surface emission vector that are directly feed into
RT4's core solver (instead of their RT4-internal calculation as
used by *RT4CalcWithRT4Surface*).

Known issues of ARTS implementation:

- TOA incoming radiation is so far assumed as blackbody cosmic
  background (temperature taken from the ARTS-internal constant).

The keyword ``pfct_method`` allows to choose how to extract the
scattering matrix, by chosing one specific temperature grid point
from the single scattering data: 'low' choses the lowest T-point,
'high' the highest T-point, and 'median' the median T-point. As
different scattering elements can have different temperature grids,
the actual temperature value used can differ between the scattering
elements. Note that this keyword solely affects the scattering matrix;
extinction matrix and absorption vector are always interpolated to
the actual temperature.
)--",
      .author = {"Jana Mendrok"},
      .out = {"cloudbox_field", "za_grid", "aa_grid"},

      .in = {"atmfields_checked",
             "atmgeom_checked",
             "scat_data_checked",
             "cloudbox_checked",
             "cloudbox_on",
             "cloudbox_limits",
             "propmat_clearsky_agenda",
             "surface_rtprop_agenda",
             "pnd_field",
             "atm_field",
             "scat_data",
             "abs_species",
             "f_grid",
             "surface_field"},
      .gin = {"nstreams",
              "pfct_method",
              "quad_type",
              "add_straight_angles",
              "pfct_aa_grid_size",
              "auto_inc_nstreams",
              "robust",
              "za_interp_order",
              "cos_za_interp",
              "max_delta_tau"},
      .gin_type = {"Index",
                   "String",
                   "String",
                   "Index",
                   "Index",
                   "Index",
                   "Index",
                   "Index",
                   "Index",
                   "Numeric"},
      .gin_value = {Index{16},
                    String("median"),
                    String("D"),
                    Index{1},
                    Index{19},
                    Index{0},
                    Index{0},
                    Index{1},
                    Index{0},
                    Numeric{1e-6}},
      .gin_desc =
          {R"--(Number of polar angle directions (streams) in RT4 solution (must be an even number).)--",
           R"--(Flag which method to apply to derive phase function (for available options see above).)--",
           R"--(Flag which quadrature to apply in RT4 solution (for available options see above).)--",
           R"--(Flag whether to include nadir and zenith as explicit directions (only effective for quad_type G and D).)--",
           R"--(Number of azimuthal angle grid points to consider in Fourier series decomposition of scattering matrix (only applied for randomly oriented scattering elements))--",
           R"--(Flag whether to internally increase nstreams (individually per frequency) if norm of (bulk) scattering matrix is not preserved properly. If 0, no adaptation is done. Else ``auto_inc_nstreams`` gives the maximum number of streams to increase to. Note that the output *cloudbox_field* remains with angular dimension of ``nstreams``, only the internal solution is adapted (and then interpolated to the lower-resolution output angular grid).)--",
           R"--(For ``auto_inc_nstreams``>0, flag whether to not fail even if scattering matrix norm is not preserved when maximum stream number is reached. Internal RT4 calculations is then performed with nstreams=``auto_inc_nstreams``.)--",
           R"--(For ``auto_inc_nstreams``>0, polar angle interpolation order for interpolation from internal increased stream to originally requested nstreams-ifield.)--",
           R"--(For ``auto_inc_nstreams``>0, flag whether to do polar angle interpolation in cosine (='mu') space.)--",
           R"--(Maximum optical depth of infinitesimal layer (where single scattering approximation is assumed to apply).)--"},
      .pass_workspace = true,

  };

  wsm_data["RT4CalcWithRT4Surface"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(As RT4Calc except for using RT4's proprietary surface type handling.

This WSM is only indented for testing purposes.

The following surface type/property methods are available and
require the the following input:

- 'L'ambertian: *surface_scalar_reflectivity*, *surface_skin_t*
- 'F'resnel: *surface_complex_refr_index*, *surface_skin_t*
- 'S'pecular: *surface_reflectivity*, *surface_skin_t*

'L' and 'F' use proprietary RT4 methods, 'S' uses RT4's Fresnel
methods modified to behave similar to ARTS'
*surfaceFlatReflectivity*.
)--",
      .author = {"Jana Mendrok"},
      .out = {"cloudbox_field", "za_grid", "aa_grid"},

      .in = {"atmfields_checked",
             "atmgeom_checked",
             "scat_data_checked",
             "cloudbox_checked",
             "cloudbox_on",
             "cloudbox_limits",
             "propmat_clearsky_agenda",
             "pnd_field",
             "atm_field",
             "scat_data",
             "abs_species",
             "f_grid",
             "surface_field",
             "surface_skin_t",
             "surface_scalar_reflectivity",
             "surface_reflectivity",
             "surface_complex_refr_index"},
      .gin = {"nstreams",
              "pfct_method",
              "ground_type",
              "quad_type",
              "add_straight_angles",
              "pfct_aa_grid_size",
              "auto_inc_nstreams",
              "robust",
              "za_interp_order",
              "cos_za_interp",
              "max_delta_tau"},
      .gin_type = {"Index",
                   "String",
                   "String",
                   "String",
                   "Index",
                   "Index",
                   "Index",
                   "Index",
                   "Index",
                   "Index",
                   "Numeric"},
      .gin_value = {Index{16},
                    String("median"),
                    String("A"),
                    String("D"),
                    Index{1},
                    Index{19},
                    Index{0},
                    Index{0},
                    Index{1},
                    Index{0},
                    Numeric{1e-6}},
      .gin_desc =
          {R"--(Number of polar angle directions (streams) in RT4 solution (must be an even number).)--",
           R"--(Flag which method to apply to derive phase function (for available options see above).)--",
           R"--(Flag which surface type/surface property method to use (for available options see above).)--",
           R"--(Flag which quadrature to apply in RT4 solution (for available options see above).)--",
           R"--(Flag whether to include nadir and zenith as explicit directions (only effective for quad_type G and D).)--",
           R"--(Number of azimuthal angle grid points to consider in Fourier series decomposition of scattering matrix (only applied for randomly oriented scattering elements))--",
           R"--(Flag whether to internally increase nstreams (individually per frequency) if norm of (bulk) scattering matrix is not preserved properly. If 0, no adaptation is done. Else ``auto_inc_nstreams`` gives the maximum number of streams to increase to.)--",
           R"--(For ``auto_inc_nstreams``>0, flag whether to not fail even if scattering matrix norm is not preserved when maximum stream number is reached. Internal RT4 calculations is then performed with nstreams=``auto_inc_nstreams``.)--",
           R"--(For ``auto_inc_nstreams``>0, polar angle interpolation order for interpolation from internal increased stream to originally requested nstreams-ifield.)--",
           R"--(For ``auto_inc_nstreams``>0, flag whether to do polar angle interpolation in cosine (='mu') space.)--",
           R"--(Maximum optical depth of infinitesimal layer (where single scattering approximation is assumed to apply).)--"},
      .pass_workspace = true,

  };

  wsm_data["RT4Test"] = WorkspaceMethodInternalRecord{
      .desc = R"--(RT4 validation test.

Executes test case testc shipped with PolRadTran/RT4 code (but uses
data files converted to arts-xml). Output written to (xml-)file.
)--",
      .author = {"Jana Mendrok"},

      .gout = {"out_rad"},
      .gout_type = {"Tensor4"},
      .gout_desc = {R"--(RT4 testc calculation results.)--"},

      .gin = {"datapath"},
      .gin_type = {"String"},
      .gin_value = {String("artscomponents/polradtran/testdata/")},
      .gin_desc =
          {R"--(Folder containing arts-xml converted test case input data.)--"},

  };

  wsm_data["RadarOnionPeelingTableCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Creates a radar inversion table.

This method is tailored to make inversion tables that fit
*particle_bulkpropRadarOnionPeeling*. See that method for
format of the table.

The method needs to be called twice to form a complete table,
once for liquid and ice hydrometeors. The table can be empty at
the first call.

The input data (*scat_data* etc.) must match two scattering
species and a single frequency (the one of the radar).
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"invtable"},
      .gout_type = {"ArrayOfGriddedField3"},
      .gout_desc = {R"--()--"},
      .in = {"f_grid",
             "scat_species",
             "scat_data",
             "scat_meta",
             "pnd_agenda_array",
             "pnd_agenda_array_input_names"},
      .gin = {"i_species",
              "dbze_grid",
              "t_grid",
              "wc_min",
              "wc_max",
              "ze_tref",
              "k2"},
      .gin_type = {"Index",
                   "Vector",
                   "Vector",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    Numeric{1e-8},
                    Numeric{2e-2},
                    Numeric{273.15},
                    Numeric{-1}},
      .gin_desc =
          {R"--(Index of *scat_species* to do. Can be 0 or 1.)--",
           R"--(Grid of dBZe values to use for the table.)--",
           R"--(Temperature grid to use for the table.)--",
           R"--(A water content value that gives a dBZe smaller than first value of ``dbze_grid``.)--",
           R"--(A water content value that gives a dBZe larger than last value of ``dbze_grid``.)--",
           R"--(Reference temperature for conversion to Ze. See further *yRadar*.)--",
           R"--(Reference dielectric factor. See further *yRadar*.)--"},
      .pass_workspace = true,

  };

  wsm_data["RadiationBackgroundCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Wraps executing *rte_background_agenda*.
)--",
      .author = {"Richard Larsson"},
      .out = {"background_rad", "diy_dx"},

      .in = {"ppath",
             "atm_field",
             "f_grid",
             "iy_transmittance",
             "background_transmittance",
             "jacobian_do",
             "rte_background_agenda"},

      .pass_workspace = true,

  };

  wsm_data["RadiationFieldSpectralIntegrate"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Integrates fields like *spectral_irradiance_field* or *spectral_radiance_field* over frequency.

Important, the first dimension must be the frequency dimension!
If a field  like *spectral_radiance_field* is input, the stokes dimension
is also removed.
)--",
      .author = {"Manfred Brath"},

      .gout = {"radiation_field"},
      .gout_type = {"Tensor4, Tensor5"},
      .gout_desc =
          {R"--(Field similar to irradiance field or spectral irradiance field)--"},
      .in = {"f_grid"},
      .gin = {"spectral_radiation_field"},
      .gin_type = {"Tensor5, Tensor7"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(Field similar to spectral irradiance field, spectral radiance field)--"},

  };

  wsm_data["RadiativePropertiesCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Wraps executing *ppvar_rtprop_agenda*.
)--",
      .author = {"Richard Larsson"},
      .out = {"ppvar_propmat",
              "ppvar_dpropmat",
              "ppvar_src",
              "ppvar_dsrc",
              "ppvar_tramat",
              "ppvar_dtramat",
              "ppvar_distance",
              "ppvar_ddistance",
              "ppvar_cumtramat"},

      .in = {"ppath",
             "ppvar_atm",
             "ppvar_f",
             "jacobian_do",
             "ppvar_rtprop_agenda"},

      .pass_workspace = true,

  };

  wsm_data["RationalAdd"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds a Rational and a value (output = input + value).

The result can either be stored in the same or another Rational.
)--",
      .author = {"Richard Larsson"},

      .gout = {"output"},
      .gout_type = {"Rational"},
      .gout_desc = {R"--(Output Rational.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Rational", "Rational"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Rational.)--", R"--(Value to add.)--"},

  };

  wsm_data["RationalDivide"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Divides a Rational with a value (output = input / value).

The result can either be stored in the same or another Rational.
)--",
      .author = {"Richard Larsson"},

      .gout = {"output"},
      .gout_type = {"Rational"},
      .gout_desc = {R"--(Output Rational.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Rational", "Rational"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Rational.)--", R"--(Denominator.)--"},

  };

  wsm_data["RationalMultiply"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Multiplies a Rational with a value (output = input * value).

The result can either be stored in the same or another Rational.
)--",
      .author = {"Richard Larsson"},

      .gout = {"output"},
      .gout_type = {"Rational"},
      .gout_desc = {R"--(Output Rational.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Rational", "Rational"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Rational.)--", R"--(Multiplier.)--"},

  };

  wsm_data["RationalSubtract"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Subtracts a Rational value (output = input - value).

The result can either be stored in the same or another Rational.
)--",
      .author = {"Richard Larsson"},

      .gout = {"output"},
      .gout_type = {"Rational"},
      .gout_desc = {R"--(Output Rational.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Rational", "Rational"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Rational.)--", R"--(Subtrahend.)--"},

  };

  wsm_data["ReadARTSCAT"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads an old ArrayOfLineRecord ARTSCAT file

Note that the ARTSCAT-5 had quantum numbers and options
stored inside it but that the options will overwrite that
information.  Be careful setting the options!
)--",
      .author = {"Stefan Buehler", "Richard Larsson"},
      .out = {"abs_lines"},

      .gin = {"filename",
              "fmin",
              "fmax",
              "globalquantumnumbers",
              "localquantumnumbers",
              "normalization_option",
              "mirroring_option",
              "population_option",
              "lineshapetype_option",
              "cutoff_option",
              "cutoff_value",
              "linemixinglimit_value"},
      .gin_type = {"String",
                   "Numeric",
                   "Numeric",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "Numeric",
                   "Numeric"},
      .gin_value = {std::nullopt,
                    Numeric{-1e99},
                    Numeric{1e99},
                    String(""),
                    String(""),
                    String("None"),
                    String("None"),
                    String("LTE"),
                    String("VP"),
                    String("None"),
                    Numeric{750e9},
                    Numeric{-1}},
      .gin_desc = {R"--(Name of the ARTSCAT file)--",
                   R"--(Minimum frequency of read lines)--",
                   R"--(Maximum frequency of read lines)--",
                   R"--(Global quantum number list (space-separated))--",
                   R"--(Local quantum number list (space-separated))--",
                   R"--(Normalization option, see *abs_linesNormalization*)--",
                   R"--(Mirroring option, see *abs_linesMirroring*)--",
                   R"--(Population option, see *abs_linesPopulation*)--",
                   R"--(Lineshape option, see *abs_linesLineShapeType*)--",
                   R"--(Cutoff option, see *abs_linesCutoff*)--",
                   R"--(Cutoff value, see *abs_linesCutoff*)--",
                   R"--(Line mixing limit, see *abs_linesLinemixingLimit*)--"},

  };

  wsm_data["ReadArrayOfARTSCAT"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads an old Array<ArrayOfLineRecord> ARTSCAT file.

Note that the ARTSCAT-5 had quantum numbers and options
stored inside it but that the options will overwrite that
information.  Be careful setting the options!
)--",
      .author = {"Stefan Buehler", "Richard Larsson"},
      .out = {"abs_lines"},

      .gin = {"filename",
              "fmin",
              "fmax",
              "globalquantumnumbers",
              "localquantumnumbers",
              "normalization_option",
              "mirroring_option",
              "population_option",
              "lineshapetype_option",
              "cutoff_option",
              "cutoff_value",
              "linemixinglimit_value"},
      .gin_type = {"String",
                   "Numeric",
                   "Numeric",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "Numeric",
                   "Numeric"},
      .gin_value = {std::nullopt,
                    Numeric{-1e99},
                    Numeric{1e99},
                    String(""),
                    String(""),
                    String("None"),
                    String("None"),
                    String("LTE"),
                    String("VP"),
                    String("None"),
                    Numeric{750e9},
                    Numeric{-1}},
      .gin_desc = {R"--(Name of the ARTSCAT file)--",
                   R"--(Minimum frequency of read lines)--",
                   R"--(Maximum frequency of read lines)--",
                   R"--(Global quantum number list (space-separated))--",
                   R"--(Local quantum number list (space-separated))--",
                   R"--(Normalization option, see *abs_linesNormalization*)--",
                   R"--(Mirroring option, see *abs_linesMirroring*)--",
                   R"--(Population option, see *abs_linesPopulation*)--",
                   R"--(Lineshape option, see *abs_linesLineShapeType*)--",
                   R"--(Cutoff option, see *abs_linesCutoff*)--",
                   R"--(Cutoff value, see *abs_linesCutoff*)--",
                   R"--(Line mixing limit, see *abs_linesLinemixingLimit*)--"},

  };

  wsm_data["ReadHITRAN"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads a HITRAN .par file.

The HITRAN type switch can be:

- ``"Pre2004"``: for old format
- ``"Post2004"``: for new format
- ``"Online"``: for the online format with quantum numbers (recommended)

Be careful setting the options!

Note that the isoptopologues in Hitran changes between its versions.
We support only one version of Hitran isoptopologues per commit.
To read an older version of Hitran, you should down-grade to the
correct commit-version.  If you still want to make use of modern features,
you must also store the generated *abs_lines* to file and reload the
catalog in the new version of Arts.  This step needs only be done once per
version of Hitran you are using

The complete flow to downgrade, read Hitran, and update is:

1. 'git checkout <commit hash>' to get the old version of Arts
2. Compile the program
3. Run *ReadHITRAN* to get *abs_lines* of that version of Hitran
4. Run *abs_linesWriteSpeciesSplitCatalog* to store the *abs_lines* to a folder
5. 'git checkout -' to get back to your previous version of Arts
6. Compile the program
7. Use *abs_linesReadSpeciesSplitCatalog* to read what *abs_lines*

The <commit hash> required per version of Hitran are:

- Hitran 2020-XXXX: Your current version is OK.
- Hitran 2004-2016: 60a9664f69f10b3f3eef3d9456282c3638b637fc
- Hitran  pre-2004: d81802cc7fe887446715491ee8a9eab8e370a0c7
)--",
      .author = {"Hermann Berg", "Thomas Kuhn", "Richard Larsson"},
      .out = {"abs_lines"},

      .gin = {"filename",
              "fmin",
              "fmax",
              "globalquantumnumbers",
              "localquantumnumbers",
              "hitran_type",
              "normalization_option",
              "mirroring_option",
              "population_option",
              "lineshapetype_option",
              "cutoff_option",
              "cutoff_value",
              "linemixinglimit_value"},
      .gin_type = {"String",
                   "Numeric",
                   "Numeric",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "Numeric",
                   "Numeric"},
      .gin_value = {std::nullopt,
                    Numeric{-1e99},
                    Numeric{1e99},
                    String("DEFAULT_GLOBAL"),
                    String("DEFAULT_LOCAL"),
                    String("Online"),
                    String("None"),
                    String("None"),
                    String("LTE"),
                    String("VP"),
                    String("None"),
                    Numeric{750e9},
                    Numeric{-1}},
      .gin_desc =
          {R"--(Name of the HITRAN file)--",
           R"--(Minimum frequency of read lines)--",
           R"--(Maximum frequency of read lines)--",
           R"--(Global quantum number list (space-separated, default gives all))--",
           R"--(Local quantum number list (space-separated, default gives all))--",
           R"--(Method to use to read the line data)--",
           R"--(Normalization option, see *abs_linesNormalization*)--",
           R"--(Mirroring option, see *abs_linesMirroring*)--",
           R"--(Population option, see *abs_linesPopulation*)--",
           R"--(Lineshape option, see *abs_linesLineShapeType*)--",
           R"--(Cutoff option, see *abs_linesCutoff*)--",
           R"--(Cutoff value, see *abs_linesCutoff*)--",
           R"--(Line mixing limit, see *abs_linesLinemixingLimit*)--"},

  };

  wsm_data["ReadJPL"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads a JPL file.

Be careful setting the options!
)--",
      .author = {"Thomas Kuhn", "Richard Larsson"},
      .out = {"abs_lines"},

      .gin = {"filename",
              "fmin",
              "fmax",
              "globalquantumnumbers",
              "localquantumnumbers",
              "normalization_option",
              "mirroring_option",
              "population_option",
              "lineshapetype_option",
              "cutoff_option",
              "cutoff_value",
              "linemixinglimit_value"},
      .gin_type = {"String",
                   "Numeric",
                   "Numeric",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "Numeric",
                   "Numeric"},
      .gin_value = {std::nullopt,
                    Numeric{-1e99},
                    Numeric{1e99},
                    String(""),
                    String(""),
                    String("None"),
                    String("None"),
                    String("LTE"),
                    String("VP"),
                    String("None"),
                    Numeric{750e9},
                    Numeric{-1}},
      .gin_desc = {R"--(Name of the JPL file)--",
                   R"--(Minimum frequency of read lines)--",
                   R"--(Maximum frequency of read lines)--",
                   R"--(Global quantum number list (space-separated))--",
                   R"--(Local quantum number list (space-separated))--",
                   R"--(Normalization option, see *abs_linesNormalization*)--",
                   R"--(Mirroring option, see *abs_linesMirroring*)--",
                   R"--(Population option, see *abs_linesPopulation*)--",
                   R"--(Lineshape option, see *abs_linesLineShapeType*)--",
                   R"--(Cutoff option, see *abs_linesCutoff*)--",
                   R"--(Cutoff value, see *abs_linesCutoff*)--",
                   R"--(Line mixing limit, see *abs_linesLinemixingLimit*)--"},

  };

  wsm_data["ReadLBLRTM"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads a LBLRTM file.

Be careful setting the options!
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .gin = {"filename",
              "fmin",
              "fmax",
              "globalquantumnumbers",
              "localquantumnumbers",
              "normalization_option",
              "mirroring_option",
              "population_option",
              "lineshapetype_option",
              "cutoff_option",
              "cutoff_value",
              "linemixinglimit_value"},
      .gin_type = {"String",
                   "Numeric",
                   "Numeric",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "Numeric",
                   "Numeric"},
      .gin_value = {std::nullopt,
                    Numeric{-1e99},
                    Numeric{1e99},
                    String(""),
                    String(""),
                    String("None"),
                    String("None"),
                    String("LTE"),
                    String("VP"),
                    String("None"),
                    Numeric{750e9},
                    Numeric{-1}},
      .gin_desc = {R"--(Name of the LBLRTM file)--",
                   R"--(Minimum frequency of read lines)--",
                   R"--(Maximum frequency of read lines)--",
                   R"--(Global quantum number list (space-separated))--",
                   R"--(Local quantum number list (space-separated))--",
                   R"--(Normalization option, see *abs_linesNormalization*)--",
                   R"--(Mirroring option, see *abs_linesMirroring*)--",
                   R"--(Population option, see *abs_linesPopulation*)--",
                   R"--(Lineshape option, see *abs_linesLineShapeType*)--",
                   R"--(Cutoff option, see *abs_linesCutoff*)--",
                   R"--(Cutoff value, see *abs_linesCutoff*)--",
                   R"--(Line mixing limit, see *abs_linesLinemixingLimit*)--"},

  };

  wsm_data["ReadNetCDF"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads a workspace variable from a NetCDF file.

This method can read variables of any group.

If the filename is omitted, the variable is read
from <basename>.<variable_name>.nc.
)--",
      .author = {"Oliver Lemke"},

      .gout = {"output"},
      .gout_type =
          {"Vector, Matrix, Tensor3, Tensor4, Tensor5, ArrayOfVector, ArrayOfIndex, ArrayOfMatrix, GasAbsLookup"},
      .gout_desc = {R"--(Variable to be read.)--"},

      .gin = {"filename"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Name of the NetCDF file.)--"},

      .pass_names = true};

  wsm_data["ReadSplitARTSCAT"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads several old ArrayOfLineRecord ARTSCAT file

Note that the ARTSCAT-5 had quantum numbers and options
stored inside it but that the options will overwrite that
information.  Be careful setting the options!
)--",
      .author = {"Oliver Lemke", "Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_species"},
      .gin = {"basename",
              "fmin",
              "fmax",
              "globalquantumnumbers",
              "localquantumnumbers",
              "ignore_missing",
              "normalization_option",
              "mirroring_option",
              "population_option",
              "lineshapetype_option",
              "cutoff_option",
              "cutoff_value",
              "linemixinglimit_value"},
      .gin_type = {"String",
                   "Numeric",
                   "Numeric",
                   "String",
                   "String",
                   "Index",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "Numeric",
                   "Numeric"},
      .gin_value = {std::nullopt,
                    Numeric{-1e99},
                    Numeric{1e99},
                    String(""),
                    String(""),
                    Index{0},
                    String("None"),
                    String("None"),
                    String("LTE"),
                    String("VP"),
                    String("None"),
                    Numeric{750e9},
                    Numeric{-1}},
      .gin_desc =
          {R"--(Path to the files)--",
           R"--(Minimum frequency of read lines)--",
           R"--(Maximum frequency of read lines)--",
           R"--(Global quantum number list (space-separated))--",
           R"--(Local quantum number list (space-separated))--",
           R"--(Ignores instead of throws if an *abs_species* is missing)--",
           R"--(Normalization option, see *abs_linesNormalization*)--",
           R"--(Mirroring option, see *abs_linesMirroring*)--",
           R"--(Population option, see *abs_linesPopulation*)--",
           R"--(Lineshape option, see *abs_linesLineShapeType*)--",
           R"--(Cutoff option, see *abs_linesCutoff*)--",
           R"--(Cutoff value, see *abs_linesCutoff*)--",
           R"--(Line mixing limit, see *abs_linesLinemixingLimit*)--"},

  };

  wsm_data["ReadXML"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads a workspace variable from an XML file.

This method can read variables of any group.

If the filename is omitted, the variable is read
from <basename>.<variable_name>.xml.
If the given filename does not exist, this method will
also look for files with an added .xml, .xml.gz and .gz extension
)--",
      .author = {"Oliver Lemke"},

      .gout = {"output"},
      .gout_type = {"Any"},
      .gout_desc = {R"--(Variable to be read.)--"},

      .gin = {"filename"},
      .gin_type = {"String"},
      .gin_value = {String("")},
      .gin_desc = {R"--(Name of the XML file.)--"},

      .pass_names = true};

  wsm_data["ReadXMLIndexed"] = WorkspaceMethodInternalRecord{
      .desc = R"--(As *ReadXML*, but reads indexed file names.

The variable is read from a file with name::

   <filename>.<file_index>.xml.

where <file_index> is the value of ``file_index``.

This means that ``filename`` shall here not include the .xml
extension. Omitting filename works as for *ReadXML*.
)--",
      .author = {"Oliver Lemke"},

      .gout = {"output"},
      .gout_type = {"Any"},
      .gout_desc = {R"--(Workspace variable to be read.)--"},
      .in = {"file_index"},
      .gin = {"filename", "digits"},
      .gin_type = {"String", "Index"},
      .gin_value = {String(""), Index{0}},
      .gin_desc =
          {R"--(File name. See above.)--",
           R"--(Equalize the widths of all numbers by padding with zeros as necessary. 0 means no padding (default).)--"},

      .pass_names = true};

  wsm_data["ReadXsecData"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads HITRAN Crosssection coefficients

Reads coefficient files for HITRAN Xsec species defined
in *abs_species*.
)--",
      .author = {"Oliver Lemke"},
      .out = {"xsec_fit_data"},

      .in = {"abs_species"},
      .gin = {"basename"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Basepath to the files)--"},

  };

  wsm_data["Reduce"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reduces a larger class to a smaller class of same size.

The Reduce command reduces all "1"-dimensions to nil.  Examples:

1) 1 Vector can be reduced to a Numeric
2) 2x1 Matrix can be reduced to 2 Vector
3) 1x3x1 Tensor3 can be reduced to 3 Vector
4) 1x1x1x1 Tensor4 can be reduced to a Numeric
5) 3x1x4x1x5 Tensor5 can only be reduced to 3x4x5 Tensor3
6) 1x1x1x1x2x3 Tensor6 can be reduced to 2x3 Matrix
7) 2x3x4x5x6x7x1 Tensor7 can be reduced to 2x3x4x5x6x7 Tensor6

And so on
)--",
      .author = {"Oliver Lemke", "Richard Larsson"},

      .gout = {"o"},
      .gout_type =
          {"Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Vector, Vector, Vector, Vector, Vector, Vector, Matrix, Matrix, Matrix, Matrix, Matrix, Tensor3, Tensor3, Tensor3, Tensor3, Tensor4, Tensor4, Tensor4, Tensor5, Tensor5, Tensor6"},
      .gout_desc = {R"--(Reduced form of input.)--"},

      .gin = {"i"},
      .gin_type =
          {"Vector, Matrix, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7, Matrix, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7, Tensor4, Tensor5, Tensor6, Tensor7, Tensor5, Tensor6, Tensor7, Tensor6, Tensor7, Tensor7"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Over-dimensioned input)--"},

  };

  wsm_data["ScatElementsPndAndScatAdd"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds single scattering data and particle number density for
individual scattering elements.

The methods reads the specified files and appends the obtained data
to *scat_data* and *pnd_field_raw*. Scattering data is appended to
the current last existing scattering species in *scat_data*.
)--",
      .author = {"Claudia Emde, Jana Mendrok"},
      .out = {"scat_data_raw", "pnd_field_raw"},

      .in = {"scat_data_raw", "pnd_field_raw"},
      .gin = {"scat_data_files", "pnd_field_files"},
      .gin_type = {"ArrayOfString", "ArrayOfString"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(List of names of single scattering data files.)--",
           R"--(List of names of the corresponding pnd_field files.)--"},

  };

  wsm_data["ScatElementsSelect"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Allows to limit considered scattering elements according to size.

Scattering elements of a specified scattering species are removed
from *scat_data_raw* and *scat_meta*, i.e. removed from further
calculations, if their particle size exceeds the specified limits.
Specification of the scattering species is done by name matching the
scattering species name part of *scat_species* tag.
As size parameter, all size parameters reported by the meta data
can be used (see *scat_meta_single* for offered parameters and
their naming).
)--",
      .author = {"Daniel Kreyling, Oliver Lemke, Jana Mendrok"},
      .out = {"scat_data_raw", "scat_meta"},

      .in = {"scat_data_raw", "scat_meta", "scat_species"},
      .gin =
          {"species", "sizeparam", "sizemin", "sizemax", "tolerance", "delim"},
      .gin_type =
          {"String", "String", "Numeric", "Numeric", "Numeric", "String"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    Numeric{0.},
                    Numeric{-1.},
                    Numeric{1e-6},
                    String("-")},
      .gin_desc =
          {R"--(Species on which to apply size selection.)--",
           R"--(Size parameter to apply for size selection.)--",
           R"--(Minimum size [m] of the scattering elements to consider)--",
           R"--(Maximum size [m] of the scattering elements to consider (if negative, no max. limitation is applied).)--",
           R"--(Relative numerical tolerance of size limit values.)--",
           R"--(Delimiter string of *scat_species* elements.)--"},

  };

  wsm_data["ScatElementsToabs_speciesAdd"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Appends scattering elements to *abs_species* including reading
single scattering data and corresponding pnd field.

The methods reads the specified single scattering and pnd_field
data of individual scattering elements and appends the obtained data
to *scat_data* (appending to its last scattering species) and
``vmr_field_raw``. Per scattering element, it also appends one
instance of species 'particles' to *abs_species*.
)--",
      .author = {"Jana Mendrok"},
      .out = {"scat_data_raw",
              "atm_field",
              "abs_species",
              "propmat_clearsky_agenda_checked"},

      .in = {"scat_data_raw",
             "atm_field",
             "abs_species",
             "propmat_clearsky_agenda_checked",
             "f_grid"},
      .gin = {"scat_data_files", "pnd_field_files"},
      .gin_type = {"ArrayOfString", "ArrayOfString"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(List of names of single scattering data files.)--",
           R"--(List of names of the corresponding pnd_field files.)--"},

  };

  wsm_data["ScatSpeciesExtendTemperature"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Extends valid temperature range of single scattering data.

The method allows to extend the temperature range of given single
scattering data by duplicating optical property data at the low
and/or high limits of the associated temperature grid. ``T_low`` and
``T_high`` specify the temperature grid points that are added.
Extension is only performed if ``T_low`` is lower and ``T_high`` is
higher than the original lowest and highest temperatures,
respectively, and if the original data contains more than one
temperature grid point (i.e., when not assumed constant anyways).

The method is thought, e.g., for atmospheric ice falling into
atmospheric layers with temperatures above the melting point of
ice, where ambient and particle temperature deviate (as long as
frozen the ice temperature remains at the melting point
temperature). It is not internally checked, whether the original
data includes the melting point.
The method can be used in a wider sense. However, it remains in the
responsibility of the user to apply the method in a meaningful
sense and on meaningful single scattering data.

The temperature extension is applied on all scattering elements of
a scattering species. If *scat_species* is defined, ``species`` can
be used to select the species on which the extension shall be
applied comparing ``species`` with the scattering species name part
of *scat_species*. If no ``species`` is specified, the method is
applied on the current last existing scattering species in
*scat_data*. Through the latter the method can be applied for cases
when *scat_species* is not defined (e.g. when *pnd_field* data is
created externally instead of from hydrometeor fields 
)--",
      .author = {"Jana Mendrok"},
      .out = {"scat_data_raw"},

      .in = {"scat_data_raw", "scat_species"},
      .gin = {"species", "scat_species_delim", "T_low", "T_high"},
      .gin_type = {"String", "String", "Numeric", "Numeric"},
      .gin_value = {String(""), String("-"), Numeric{-1.}, Numeric{-1.}},
      .gin_desc =
          {R"--(Scattering species to act on (see WSM description for details).)--",
           R"--(Delimiter string of *scat_species* elements.)--",
           R"--(Temperature grid extension point at low temperature limit.)--",
           R"--(Temperature grid extension point at high temperature limit.)--"},

  };

  wsm_data["ScatSpeciesInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initializes the scattering species related data variables.

This method initializes the *scat_species* WSV, the variables that
will hold the raw optical properties and the raw particle number
distributions of the scattering elements (*scat_data_raw* and
*pnd_field_raw*, respectively) as well as the one holding the meta
information about the scattering elements (*scat_meta*).

This method has to be executed before WSM reading/adding to the
said variable, e.g. before *ScatSpeciesPndAndScatAdd*.
)--",
      .author = {"Jana Mendrok"},
      .out = {"scat_species",
              "scat_data_raw",
              "scat_meta",
              "scat_data_checked",
              "pnd_field_raw"},

  };

  wsm_data["ScatSpeciesPndAndScatAdd"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Adds single scattering data and particle number densities for one
scattering species.

The WSV *pnd_field_raw* containing particle number densities for
all scattering species can be generated outside ARTS, for example
by using PyARTS or atmlab. This method reads this data as well as
its corresponding single scattering data, which is added as a new
scattering species to *scat_data*.
This method needs as input an ArrayOfString holding the filenames
of the single scattering data for each scattering element and a
file containing the corresponding *pnd_field_raw*. In contrast to
the scattering data, the pnd-fields are stored in a single XML-file
containing an ArrayofGriddedField3, i.e. holding the pnd-field data
of all scattering elements.

Important note:
The order of the filenames for the scattering data files has to
correspond to the order of the pnd-fields, stored in the variable
*pnd_field_raw*.
)--",
      .author = {"Claudia Emde, Jana Mendrok"},
      .out = {"scat_data_raw", "pnd_field_raw"},

      .in = {"scat_data_raw", "pnd_field_raw"},
      .gin = {"scat_data_files", "pnd_fieldarray_file"},
      .gin_type = {"ArrayOfString", "String"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(Array of names of files containing the single scattering data.)--",
           R"--(Name of file holding the corresponding array of pnd_field data.)--"},

  };

  wsm_data["ScatSpeciesScatAndMetaRead"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads single scattering data and scattering meta data for one
scattering species.

This method takes a string array as input containing the location
(path and filename) of the single scattering data. Location of
corresponding scattering meta data is derived applying a naming
convention: ending '.xml*' is replaced by '.meta.xml' (search for
zipped files is done automatically).

All scattering elements read in one call of the method are assigned
to one and the same scattering species. That is, reading in data for
a bunch of scattering species can be realized by multiple call of
this method. Assignment to scattering species is in the order of the
calls (i.e., first method call reads data for first *scat_species*
entry, second call for second scat_species entry and so on).
Note that no two scattering elements of the same scattering species
are allowed to be equal in size*

Important note:
The order of the filenames for the single scattering data files has to
exactly correspond to the order of the scattering meta data files.
)--",
      .author = {"Daniel Kreyling, Oliver Lemke, Jana Mendrok"},
      .out = {"scat_data_raw", "scat_meta"},

      .in = {"scat_data_raw", "scat_meta"},
      .gin = {"scat_data_files"},
      .gin_type = {"ArrayOfString"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Array of single scattering data file names.)--"},

  };

  wsm_data["ScatSpeciesSizeMassInfo"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Derives size and mass information for a scattering species.
This method assumes that the mass-size relationship can described
by *scat_species_a* and *scat_species_b*. See documentation of 
*scat_species_a* for details.

The quantity to be used as size descriptor is here denoted as x, and
is selected by setting ``x_unit``. The options are:

- ``"dveq"``: The size grid is set to scat_meta.diameter_volume_equ
- ``"dmax"``: The size grid is set to scat_meta.diameter_max
- ``"area"``: The size grid is set to scat_meta.diameter_area_equ_aerodynamical
- ``"mass"``: The size grid is set to scat_meta.mass

This selection determines *scat_species_x*.

The parameters *scat_species_a* and *scat_species_b* are determined by
a numeric fit between *scat_species_x* and corresponding masses in
*scat_meta*. This fit is performed over sizes inside the range
[x_fit_start,x_fit_end]. This range is allowed to be broader than
the coverage of *scat_species_x*. There must be at least two sizes
inside [x_fit_start,x_fit_end].
)--",
      .author = {"Manfred Brath", "Jana Mendrok", "Patrick Eriksson"},
      .out = {"scat_species_x", "scat_species_a", "scat_species_b"},

      .in = {"scat_meta"},
      .gin =
          {"species_index", "x_unit", "x_fit_start", "x_fit_end", "do_only_x"},
      .gin_type = {"Index", "String", "Numeric", "Numeric", "Index"},
      .gin_value =
          {std::nullopt, std::nullopt, Numeric{0}, Numeric{1e9}, Index{0}},
      .gin_desc =
          {R"--(Take data from scattering species of this index (0-based) in *scat_meta*.)--",
           R"--(Unit for size grid, allowed options listed above.)--",
           R"--(Smallest size to consider in fit to determine a and b.)--",
           R"--(Largest size to consider in fit to determine a and b.)--",
           R"--(A flag to deactivate calculation of a and b, to possibly save some time. The a and b parameters are then set to -1.Default is to calculate a and b.)--"},

  };

  wsm_data["Select"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Method to select some elements from one array and copy them to
a new array. (Works also for vectors.)

This works also for higher dimensional objects, where the selection is
always performed in the first dimension.

If needleindexes is set to [-1], all elements are copied.
For example:

Select(y,x,[0,3])

will select the first and fourth row of matrix x and copy them to the
output matrix y.

Note that it is even safe to use this method if needles and haystack
are the same variable.
)--",
      .author = {"Oliver Lemke"},

      .gout = {"needles"},
      .gout_type =
          {"ArrayOfAbsorptionLines, ArrayOfAgenda, ArrayOfArrayOfAbsorptionLines, ArrayOfArrayOfGriddedField1, ArrayOfArrayOfGriddedField2, ArrayOfArrayOfGriddedField3, ArrayOfArrayOfIndex, ArrayOfArrayOfMatrix, ArrayOfArrayOfMuelmatMatrix, ArrayOfArrayOfMuelmatVector, ArrayOfArrayOfPropmatMatrix, ArrayOfArrayOfPropmatVector, ArrayOfArrayOfScatteringMetaData, ArrayOfArrayOfSingleScatteringData, ArrayOfArrayOfSpeciesTag, ArrayOfArrayOfStokvecMatrix, ArrayOfArrayOfStokvecVector, ArrayOfArrayOfString, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6, ArrayOfArrayOfTime, ArrayOfArrayOfVector, ArrayOfAtmPoint, ArrayOfCIARecord, ArrayOfGriddedField1, ArrayOfGriddedField2, ArrayOfGriddedField3, ArrayOfGriddedField4, ArrayOfIndex, ArrayOfJacobianTarget, ArrayOfMatrix, ArrayOfMuelmatMatrix, ArrayOfMuelmatVector, ArrayOfPpath, ArrayOfPropmatMatrix, ArrayOfPropmatVector, ArrayOfQuantumIdentifier, ArrayOfRetrievalQuantity, ArrayOfScatteringMetaData, ArrayOfSingleScatteringData, ArrayOfSparse, ArrayOfSpeciesTag, ArrayOfStokvecMatrix, ArrayOfStokvecVector, ArrayOfString, ArrayOfSun, ArrayOfTelsemAtlas, ArrayOfTensor3, ArrayOfTensor4, ArrayOfTensor5, ArrayOfTensor6, ArrayOfTensor7, ArrayOfTime, ArrayOfVector, ArrayOfXsecRecord, Vector, Matrix, Sparse"},
      .gout_desc =
          {R"--(Selected elements. Must have the same variable type as haystack.)--"},

      .gin = {"haystack", "needleindexes"},
      .gin_type = {"ArrayOfAbsorptionLines, ArrayOfAgenda, ArrayOfArrayOfAbsorptionLines, ArrayOfArrayOfGriddedField1, ArrayOfArrayOfGriddedField2, ArrayOfArrayOfGriddedField3, ArrayOfArrayOfIndex, ArrayOfArrayOfMatrix, ArrayOfArrayOfMuelmatMatrix, ArrayOfArrayOfMuelmatVector, ArrayOfArrayOfPropmatMatrix, ArrayOfArrayOfPropmatVector, ArrayOfArrayOfScatteringMetaData, ArrayOfArrayOfSingleScatteringData, ArrayOfArrayOfSpeciesTag, ArrayOfArrayOfStokvecMatrix, ArrayOfArrayOfStokvecVector, ArrayOfArrayOfString, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6, ArrayOfArrayOfTime, ArrayOfArrayOfVector, ArrayOfAtmPoint, ArrayOfCIARecord, ArrayOfGriddedField1, ArrayOfGriddedField2, ArrayOfGriddedField3, ArrayOfGriddedField4, ArrayOfIndex, ArrayOfJacobianTarget, ArrayOfMatrix, ArrayOfMuelmatMatrix, ArrayOfMuelmatVector, ArrayOfPpath, ArrayOfPropmatMatrix, ArrayOfPropmatVector, ArrayOfQuantumIdentifier, ArrayOfRetrievalQuantity, ArrayOfScatteringMetaData, ArrayOfSingleScatteringData, ArrayOfSparse, ArrayOfSpeciesTag, ArrayOfStokvecMatrix, ArrayOfStokvecVector, ArrayOfString, ArrayOfSun, ArrayOfTelsemAtlas, ArrayOfTensor3, ArrayOfTensor4, ArrayOfTensor5, ArrayOfTensor6, ArrayOfTensor7, ArrayOfTime, ArrayOfVector, ArrayOfXsecRecord, Vector, Matrix, Sparse",
                   "ArrayOfIndex"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(Variable to select from. May be the same variable as needles.)--",
           R"--(The elements to select (zero based indexing, as always.))--"},

  };

  wsm_data["SetNumberOfThreads"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Change the number of threads used by ARTS.
)--",
      .author = {"Oliver Lemke"},

      .gin = {"nthreads"},
      .gin_type = {"Index"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Number of threads.)--"},

  };

  wsm_data["Sleep"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sleeps for a number of seconds
)--",
      .author = {"Richard Larsson"},

      .gin = {"gtime"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Time to sleep for in seconds)--"},

  };

  wsm_data["SparseIdentity"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Returns a sparse dentity matrix.

The size of the matrix created is n x n. Default is to return a
true identity matrix (I), but you can also select another value
along the diagonal be setting ``value``. That is, the output is
value*I.
)--",
      .author = {"Simon Pfreundschuh"},

      .gout = {"output"},
      .gout_type = {"Sparse"},
      .gout_desc = {R"--(Sparse output matrix)--"},

      .gin = {"n", "value"},
      .gin_type = {"Index", "Numeric"},
      .gin_value = {std::nullopt, Numeric{1}},
      .gin_desc = {R"--(Size of the matrix)--",
                   R"--(The value along the diagonal.)--"},

  };

  wsm_data["SparseSparseMultiply"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Multiplies a Sparse with another Sparse, result stored in Sparse.

Makes the calculation: M = M1 * M2
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"M"},
      .gout_type = {"Sparse"},
      .gout_desc =
          {R"--(Product, can be same variable as any of the inputs.)--"},

      .gin = {"M1", "M2"},
      .gin_type = {"Sparse", "Sparse"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Left sparse matrix (dimension m x n).)--",
                   R"--(Right sparse matrix (dimension n x p).)--"},

  };

  wsm_data["StringJoin"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Concatenate two or more strings.

The output string is overwritten, but is allowed to appear
in the input list. Up to 10 strings can be concatenated at once.
)--",
      .author = {"Oliver Lemke"},

      .gout = {"output"},
      .gout_type = {"String"},
      .gout_desc = {R"--(Concatenated string.)--"},

      .gin = {"in1",
              "in2",
              "in3",
              "in4",
              "in5",
              "in6",
              "in7",
              "in8",
              "in9",
              "in10"},
      .gin_type = {"String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String",
                   "String"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    String(""),
                    String(""),
                    String(""),
                    String(""),
                    String(""),
                    String(""),
                    String(""),
                    String("")},
      .gin_desc = {R"--(Input text string.)--",
                   R"--(Input text string.)--",
                   R"--(Input text string.)--",
                   R"--(Input text string.)--",
                   R"--(Input text string.)--",
                   R"--(Input text string.)--",
                   R"--(Input text string.)--",
                   R"--(Input text string.)--",
                   R"--(Input text string.)--",
                   R"--(Input text string.)--"},

  };

  wsm_data["SurfaceDummy"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Dummy method for *iy_surface_agenda*.

If you don't make use of ``surface_props_data`` and associated
variables, include this method *iy_surface_agenda*. The method
just checks that the variables of concern are set to be empty,
and you don't need to include calls of *Ignore* and *Touch* in
the agenda.

If you use a method of SurfaceSomething type, you don't need
this one.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"dsurface_rmatrix_dx", "dsurface_emission_dx"},

      .in = {"dsurface_rmatrix_dx",
             "dsurface_emission_dx",
             "surface_props_names",
             "dsurface_names",
             "jacobian_do"},

  };

  wsm_data["TMatrixTest"] = WorkspaceMethodInternalRecord{
      .desc = R"--(T-Matrix validation test.

Executes the standard test included with the T-Matrix Fortran code.
Should give the same as running the tmatrix_lp executable in
3rdparty/tmatrix/.
)--",
      .author = {"Oliver Lemke"},

  };

  wsm_data["Tensor3Add"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds a scalar value to all elements of a tensor3.

The result can either be stored in the same or another
variable.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Tensor3"},
      .gout_desc = {R"--(Output Tensor.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Tensor3", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Tensor.)--",
                   R"--(The value to be added to the tensor.)--"},

  };

  wsm_data["Tensor3ExtractFromTensor4"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Extracts a Tensor3 from a Tensor4.

Copies book, page, row or column with given Index from input Tensor4
variable to output Tensor3.
Higher order equivalent of *VectorExtractFromMatrix*.
)--",
      .author = {"Oliver Lemke"},

      .gout = {"output"},
      .gout_type = {"Tensor3"},
      .gout_desc = {R"--(Extracted tensor.)--"},

      .gin = {"input", "i", "direction"},
      .gin_type = {"Tensor4", "Index", "String"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Tensor4.)--",
                   R"--(Index of book, page, row or column to extract.)--",
                   R"--(Direction. "book" or "page" or "row" or "column".)--"},

  };

  wsm_data["Tensor3FromVector"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Forms a Tensor3 of size nx1x1 from a vector of length n.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Tensor3"},
      .gout_desc = {R"--(Output tensor.)--"},

      .gin = {"v"},
      .gin_type = {"Vector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Input vector.)--"},

  };

  wsm_data["Tensor3Multiply"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Multiplies all elements of a tensor with the specified value.

The result can either be stored in the same or another
variable.
)--",
      .author = {"Mattias Ekstrom"},

      .gout = {"output"},
      .gout_type = {"Tensor3"},
      .gout_desc = {R"--(Output Tensor.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Tensor3", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Tensor.)--",
                   R"--(The value to be multiplied with the tensor.)--"},

  };

  wsm_data["Tensor3SetConstant"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Creates a tensor and sets all elements to the specified value.

The size is determined by *ncols*, *nrows* etc.
)--",
      .author = {"Claudia Emde"},

      .gout = {"output"},
      .gout_type = {"Tensor3"},
      .gout_desc = {R"--(Variable to initialize.)--"},
      .in = {"npages", "nrows", "ncols"},
      .gin = {"value"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Tensor value.)--"},

  };

  wsm_data["Tensor4Add"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds a scalar value to all elements of a tensor4.

The result can either be stored in the same or another
variable.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Tensor4"},
      .gout_desc = {R"--(Output Tensor.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Tensor4", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Tensor.)--",
                   R"--(The value to be added to the tensor.)--"},

  };

  wsm_data["Tensor4Multiply"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Multiplies all elements of a tensor with the specified value.

The result can either be stored in the same or another
variable.
)--",
      .author = {"Mattias Ekstrom"},

      .gout = {"output"},
      .gout_type = {"Tensor4"},
      .gout_desc = {R"--(Output Tensor.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Tensor4", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Tensor.)--",
                   R"--(The value to be multiplied with the tensor.)--"},

  };

  wsm_data["Tensor4SetConstant"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Creates a tensor and sets all elements to the specified value.

The size is determined by *ncols*, *nrows* etc.
)--",
      .author = {"Claudia Emde"},

      .gout = {"output"},
      .gout_type = {"Tensor4"},
      .gout_desc = {R"--(Variable to initialize.)--"},
      .in = {"nbooks", "npages", "nrows", "ncols"},
      .gin = {"value"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Tensor value.)--"},

  };

  wsm_data["Tensor5Multiply"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Multiplies all elements of a tensor with the specified value.

The result can either be stored in the same or another
variable.
)--",
      .author = {"Mattias Ekstrom"},

      .gout = {"output"},
      .gout_type = {"Tensor5"},
      .gout_desc = {R"--(Output Tensor.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Tensor5", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Tensor.)--",
                   R"--(The value to be multiplied with the tensor.)--"},

  };

  wsm_data["Tensor5SetConstant"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Creates a tensor and sets all elements to the specified value.

The size is determined by *ncols*, *nrows* etc.
)--",
      .author = {"Claudia Emde"},

      .gout = {"output"},
      .gout_type = {"Tensor5"},
      .gout_desc = {R"--(Variable to initialize.)--"},
      .in = {"nshelves", "nbooks", "npages", "nrows", "ncols"},
      .gin = {"value"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Tensor value.)--"},

  };

  wsm_data["Tensor6Multiply"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Multiplies all elements of a tensor with the specified value.

The result can either be stored in the same or another
variable.
)--",
      .author = {"Mattias Ekstrom"},

      .gout = {"output"},
      .gout_type = {"Tensor6"},
      .gout_desc = {R"--(Output Tensor.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Tensor6", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Tensor.)--",
                   R"--(The value to be multiplied with the tensor.)--"},

  };

  wsm_data["Tensor6SetConstant"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Creates a tensor and sets all elements to the specified value.

The size is determined by *ncols*, *nrows* etc.
)--",
      .author = {"Claudia Emde"},

      .gout = {"output"},
      .gout_type = {"Tensor6"},
      .gout_desc = {R"--(Variable to initialize.)--"},
      .in = {"nvitrines", "nshelves", "nbooks", "npages", "nrows", "ncols"},
      .gin = {"value"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Tensor value.)--"},

  };

  wsm_data["Tensor7Multiply"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Multiplies all elements of a tensor with the specified value.

The result can either be stored in the same or another
variable.
)--",
      .author = {"Mattias Ekstrom"},

      .gout = {"output"},
      .gout_type = {"Tensor7"},
      .gout_desc = {R"--(Output Tensor.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Tensor7", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Tensor.)--",
                   R"--(The value to be multiplied with the tensor.)--"},

  };

  wsm_data["Tensor7SetConstant"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Creates a tensor and sets all elements to the specified value.

The size is determined by *ncols*, *nrows* etc.
)--",
      .author = {"Claudia Emde"},

      .gout = {"output"},
      .gout_type = {"Tensor7"},
      .gout_desc = {R"--(Variable to initialize.)--"},
      .in = {"nlibraries",
             "nvitrines",
             "nshelves",
             "nbooks",
             "npages",
             "nrows",
             "ncols"},
      .gin = {"value"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Tensor value.)--"},

  };

  wsm_data["TessemNNReadAscii"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Reads the initialization data for the TESSEM NeuralNet from an ASCII file.
)--",
      .author = {"Oliver Lemke"},

      .gout = {"tessem_nn"},
      .gout_type = {"TessemNN"},
      .gout_desc = {R"--(Tessem NeuralNet configuration.)--"},

      .gin = {"filename"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(NeuralNet parameters file as provided in the TESSEM 2 distribution.)--"},

  };

  wsm_data["TestBasicGeodeticAccuracy"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Tests the basic accuracy of the geodetic calculations.

Basically all geodetic calculations involves conversion to ECEF coordinates
and back. This method tests the accuracy of this conversion. 

A random position and line-of-sights is generated and the conversion are
made. The change of position is calculated as a distance. If the distance
exceeds ``max_allowed_dl`` an error is issued. Otherwise a new test is made.
This is repeated ``ntests`` times. The maximum error is returned as ``max_dl``.
The position the maximum error is returned as *rte_pos*.

Further, the maximum error for altitude, latitude etc. are returned in
GOUTs ``max_dpos`` and ``max_dlos``. This is the max absolute error for each
value separately (i.e. they can come from different tests/positions.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"rte_pos"},
      .gout = {"max_dl", "max_dpos", "max_dlos"},
      .gout_type = {"Numeric", "Vector", "Vector"},
      .gout_desc = {R"--(Maximum error in term of distance.)--",
                    R"--(The maximum error for each position component.)--",
                    R"--(The maximum error for each LOS component.)--"},
      .in = {"surface_field"},
      .gin = {"ntests", "max_allowed_dl"},
      .gin_type = {"Index", "Numeric"},
      .gin_value = {std::nullopt, Numeric{0.1}},
      .gin_desc = {R"--(Number of tests.)--",
                   R"--(Maximum allowed error in term of distance.)--"},

  };

  wsm_data["TestTessem"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Example method for TESSEM2.

When using the default neural network parameter files
from the Tessem 2 distribution, the input Vector should contain
5 elements:

- Frequency (10-700) in GHz.
- Theta (0-90) Incidence angle in degrees.
- Windspeed (0-25) at 10m (m/s)
  Higher wind speed can be used, but without garantee.
- Surface skin temperature (270-310) in K.
- Salinity (0-0.04) in kg/kg
)--",
      .author = {"Oliver Lemke"},

      .gout = {"outvalues"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Tessem output emissivity.)--"},

      .gin = {"net", "invalues"},
      .gin_type = {"TessemNN", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Tessem NeuralNet parameters.)--", R"--(Input data.)--"},

  };

  wsm_data["Touch"] = WorkspaceMethodInternalRecord{
      .desc = R"--(As *Ignore* but for agenda output.

This method is handy for use in agendas in order to suppress
warnings about not-produced output workspace variables.

What it does, in case the variable is initialized already, is:
Nothing!
In case the variable is not yet initialized, it is set to NaN.
)--",
      .author = {"Oliver Lemke"},

      .gout = {"input"},
      .gout_type = {"Any"},
      .gout_desc = {R"--(Variable to do nothing with.)--"},

  };

  wsm_data["Trapz"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Intregrates a vector of over its grid range

The method integrates y(x) by the trapezoidal method.

The vector x is the positions where the integrand, y, is known.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Numeric"},
      .gout_desc = {R"--(Value of integral)--"},

      .gin = {"gx", "gy"},
      .gin_type = {"Vector", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Grid.)--", R"--(Integrand.)--"},

  };

  wsm_data["VectorAdd"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds a scalar to all elements of a vector.

The result can either be stored in the same or another vector.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Output Vector)--"},

      .gin = {"input", "value"},
      .gin_type = {"Vector", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Vector.)--",
                   R"--(The value to be added to the vector.)--"},

  };

  wsm_data["VectorAddElementwise"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Element-wise addition of two vectors.

The method calculates c = a + b.

The variable ``b`` is allowed to have length 1, for any length of
``a``. This single value in ``b`` is then added to every element of ``a``.

The vectors ``a`` and ``c`` can be the same WSV, while ``b`` can not be
the same WSV as any of the the other vector.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"c"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Output vector)--"},

      .gin = {"a", "b"},
      .gin_type = {"Vector", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input vector.)--", R"--(Vector to be added.)--"},

  };

  wsm_data["VectorClip"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Clipping of a vector.

The input vector is copied to the output one (that can be same WSV)
but ensures that all values in ``output`` are inside the range [limit_low,
limit_high]. Where the input vector is below ``limit_low``, ``out`` is set
to ``limit_low``. And the same is performed with respect to ``limit_high``.
That is, the method works as *NumericClip* for each element of the
vector.

The method adopts the length of ``out`` when needed.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Output vector.)--"},

      .gin = {"input", "limit_low", "limit_high"},
      .gin_type = {"Vector", "Numeric", "Numeric"},
      .gin_value = {std::nullopt, -std::numeric_limits<Numeric>::infinity(), std::numeric_limits<Numeric>::infinity()},
      .gin_desc = {R"--(Input vector.)--",
                   R"--(Lower limit for clipping.)--",
                   R"--(Upper limit for clipping.)--"},

  };

  wsm_data["VectorCrop"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Keeps only values of a vector inside the specified range.

All values outside the range [min_value,max-value] are removed.
Note the default values, that basically should act as -+Inf.

The result can either be stored in the same or another vector.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Cropped vector)--"},

      .gin = {"input", "min_value", "max_value"},
      .gin_type = {"Vector", "Numeric", "Numeric"},
      .gin_value = {std::nullopt, Numeric{-99e99}, Numeric{99e99}},
      .gin_desc = {R"--(Original vector)--",
                   R"--(Minimum value to keep)--",
                   R"--(Maximum value to keep)--"},

  };

  wsm_data["VectorDivide"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Divides all elements of a vector with the same value.

The result can either be stored in the same or another vector.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Output Vector.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Vector", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Vector.)--", R"--(Denominator.)--"},

  };

  wsm_data["VectorDivideElementwise"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Element-wise division of two vectors.

The method calculates c = a / b.

The variable ``b`` is allowed to have length 1, for any length of
``a``. This single value in ``b`` is then applied to every element of ``a``.

The vectors ``a`` and ``c`` can be the same WSV, while ``b`` can not be
the same WSV as any of the the other vector.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"c"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Output vector)--"},

      .gin = {"a", "b"},
      .gin_type = {"Vector", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input vector.)--", R"--(Denominator Vector.)--"},

  };

  wsm_data["VectorExtractFromMatrix"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Extracts a Vector from a Matrix.

Copies row or column with given Index from input Matrix variable
to create output Vector.
)--",
      .author = {"Patrick Eriksson, Oliver Lemke, Stefan Buehler"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Extracted vector.)--"},

      .gin = {"input", "i", "direction"},
      .gin_type = {"Matrix", "Index", "String"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input matrix.)--",
                   R"--(Index of row or column.)--",
                   R"--(Direction. "row" or "column".)--"},

  };

  wsm_data["VectorFlip"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Flips a vector.

The output is the input vector in reversed order. The result can
either be stored in the same or another vector.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Output vector.)--"},

      .gin = {"input"},
      .gin_type = {"Vector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Input vector.)--"},

  };

  wsm_data["VectorGaussian"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Fills a vector with a Gaussian function.

The width can be set in two ways, either by standard deviation or
the full width at half maximum. Only one of the corresponding GINs
can be >0 and that value will determine the width.

The vectors *x* and *y* can be the same variable.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"gy"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Output vector.)--"},

      .gin = {"gx", "x0", "si", "fwhm"},
      .gin_type = {"Vector", "Numeric", "Numeric", "Numeric"},
      .gin_value = {std::nullopt, Numeric{0}, Numeric{-1}, Numeric{-1}},
      .gin_desc =
          {R"--(Grid of the function.)--",
           R"--(Centre/mean point of the function.)--",
           R"--(Standard deviation of the function, ignored if <=0.)--",
           R"--(Full width at half-max of the function, ignored if <=0.)--"},

  };

  wsm_data["VectorInsertGridPoints"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Insert some additional points into a grid.

This method can for example be used to add line center frequencies to
a regular frequency grid. If the original grid is [1,2,3], and the
additional points are [2.2,2.4], the result will be [1,2,2.2,2.4,3].

It is assumed that the original grid is sorted, otherwise a runtime
error is thrown. The vector with the points to insert does not have to
be sorted. If some of the input points are already in the grid, these
points are not inserted again. New points outside the original grid are
appended at the appropriate end. Input vector and output vector can be
the same.

Generic output:
 - Vector : The new grid vector.

Generic input:
 - Vector : The original grid vector.
 - Vector : The points to insert.
)--",
      .author = {"Stefan Buehler"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(The new grid vector)--"},

      .gin = {"input", "points"},
      .gin_type = {"Vector", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(The original grid vector)--",
                   R"--(The points to insert)--"},

  };

  wsm_data["VectorLinSpace"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initializes a vector with linear spacing.

The first element equals always the start value, and the spacing
equals always the step value, but the last value can deviate from
the stop value. ``step`` can be both positive and negative.

The created vector is [start, start+step, start+2*step, ...]
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Output vector.)--"},

      .gin = {"start", "stop", "step"},
      .gin_type = {"Numeric", "Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Start value.)--",
                   R"--(Maximum/minimum value of the end value)--",
                   R"--(Spacing of the vector.)--"},

  };

  wsm_data["VectorLogSpace"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initializes a vector with logarithmic spacing.

The first element equals always the start value, and the spacing
equals always the step value, but note that the last value can 
deviate from the stop value. The keyword step can be both positive
and negative.

Note, that although start has to be given in direct coordinates,
step has to be given in log coordinates.

Explicitly, the vector is:
 exp([ln(start), ln(start)+step, ln(start)+2*step, ...])
)--",
      .author = {"Stefan Buehler"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Variable to initialize.)--"},

      .gin = {"start", "stop", "step"},
      .gin_type = {"Numeric", "Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(The start value. (Direct coordinates!))--",
           R"--(The maximum value of the end value. (Direct coordinates!))--",
           R"--(The spacing of the vector. (Log coordinates!))--"},

  };

  wsm_data["VectorMatrixMultiply"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Multiply a Vector with a Matrix and store the result in another
Vector.

This just computes the normal matrix-vector product, y=M*x. It is ok
if input and output Vector are the same.
)--",
      .author = {"Stefan Buehler"},

      .gout = {"gy"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(The result of the multiplication (dimension m).)--"},

      .gin = {"M", "gx"},
      .gin_type = {"Matrix", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(The Matrix to multiply (dimension m x n).)--",
                   R"--(The original Vector (dimension n).)--"},

  };

  wsm_data["VectorMultiply"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Multiplies all elements of a vector with the same value.

The result can either be stored in the same or another vector.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Output Vector.)--"},

      .gin = {"input", "value"},
      .gin_type = {"Vector", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Vector.)--", R"--(Scaling value.)--"},

  };

  wsm_data["VectorMultiplyElementwise"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Element-wise multiplication of two vectors.

The method calculates c = a * b.

The variable ``b`` is allowed to have length 1, for any length of
``a``. This single value in ``b`` is then multiplied with every element
of ``a``.

The vectors ``a`` and ``c`` can be the same WSV, while ``b`` can not be
the same WSV as any of the the other vector.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"c"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Output vector)--"},

      .gin = {"a", "b"},
      .gin_type = {"Vector", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input vector.)--", R"--(Multiplier.)--"},

  };

  wsm_data["VectorNLinSpace"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Creates a vector with length *nelem*, equally spaced between the
given end values.

The length (*nelem*) must be larger than 1.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Variable to initialize.)--"},
      .in = {"nelem"},
      .gin = {"start", "stop"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Start value.)--", R"--(End value.)--"},

  };

  wsm_data["VectorNLinSpaceVector"] = WorkspaceMethodInternalRecord{
      .desc = R"--(As *VectorNLinSpace* but end points taken from a vector.

The method gives a vector with equidistant spacing between
first and last element of the reference vector.

The length (*nelem*) must be larger than 1.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Variable to initialize.)--"},
      .in = {"nelem"},
      .gin = {"gy"},
      .gin_type = {"Vector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Reference vector.)--"},

  };

  wsm_data["VectorNLogSpace"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Creates a vector with length *nelem*, equally logarithmically
spaced between the given end values.

The length (*nelem*) must be larger than 1.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Variable to initialize.)--"},
      .in = {"nelem"},
      .gin = {"start", "stop"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Start value.)--", R"--(End value.)--"},

  };

  wsm_data["VectorPower"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates the power of each element in a vector.

The result can either be stored in the same or another vector.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Output Vector.)--"},

      .gin = {"input", "power"},
      .gin_type = {"Vector", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Vector.)--", R"--(Power (exponent).)--"},

  };

  wsm_data["VectorReshapeMatrix"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Converts a Matrix to a Vector.

The matrix is reshaped into a vector. That is, all elements of the matrix
are kept. The elements can be extracted both in column (default) and row
order. The output vector has the same length for both options.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Created vector.)--"},

      .gin = {"input", "direction"},
      .gin_type = {"Matrix", "String"},
      .gin_value = {std::nullopt, String("column")},
      .gin_desc = {R"--(Input matrix.)--",
                   R"--(Direction. "row" or "column".)--"},

  };

  wsm_data["VectorSetConstant"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Creates a vector and sets all elements to the specified value.

The vector length is determined by *nelem*.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Variable to initialize.)--"},
      .in = {"nelem"},
      .gin = {"value"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Vector value.)--"},

  };

  wsm_data["VectorSparseMultiply"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Multiply a Vector with a Sparse and store the result in another
Vector.

This just computes the normal matrix-vector product, y=M*x, with
m being a Sparse. It is ok if input and output Vector are the same.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"gy"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(The result of the multiplication (dimension m).)--"},

      .gin = {"M", "gx"},
      .gin_type = {"Sparse", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(The Sparse to multiply (dimension m x n).)--",
                   R"--(The original Vector (dimension n).)--"},

  };

  wsm_data["VectorSubtract"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Subtracts a scalar from all elements of a vector.

The result can either be stored in the same or another vector.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"output"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Output Vector)--"},

      .gin = {"input", "value"},
      .gin_type = {"Vector", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input Vector.)--",
                   R"--(The value to be subtracted from the vector.)--"},

  };

  wsm_data["VectorSubtractElementwise"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Element-wise subtraction of two vectors.

The method calculates c = a - b.

The variable ``b`` is allowed to have length 1, for any length of
``a``. This single value in ``b`` is then subtracted to every element of ``a``.

The vectors ``a`` and ``c`` can be the same WSV, while ``b`` can not be
the same WSV as any of the the other vector.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"c"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Output vector)--"},

      .gin = {"a", "b"},
      .gin_type = {"Vector", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Input vector.)--", R"--(Vector to be subtracted.)--"},

  };

  wsm_data["WMRFSelectChannels"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Select some channels for WMRF calculation.

The HIRS fast setup consists of a precalculated frequency grid
covering all HIRS channels, and associated weights for each channel,
stored in a weight matrix. (A *sensor_response* matrix.)

If not all channels are requested for
simulation, then this method can be used to remove the unwanted
channels. It changes a number of variables in consistent fashion:

- Unwanted channels are removed from f_backend. 
- Unwanted channels are removed from wmrf_weights.
- Unnecessary frequencies are removed from f_grid.
- Unnecessary frequencies are removed from wmrf_weights.
)--",
      .author = {"Stefan Buehler"},
      .out = {"f_grid", "wmrf_weights", "f_backend"},

      .in = {"f_grid", "f_backend", "wmrf_weights", "wmrf_channels"},

  };

  wsm_data["Wigner3Init"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initialize the wigner 3 tables

The default values take about 400 Mb memory.
)--",
      .author = {"Richard Larsson"},
      .out = {"wigner_initialized"},

      .gin = {"fast_wigner_stored_symbols", "largest_wigner_symbol_parameter"},
      .gin_type = {"Index", "Index"},
      .gin_value = {Index{20000000}, Index{250}},
      .gin_desc =
          {R"--(Number of stored symbols possible before replacements)--",
           R"--(Largest symbol used for initializing factorials (e.g., largest J or L))--"},

  };

  wsm_data["Wigner3Unload"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Unloads the wigner 3 tables
)--",
      .author = {"Richard Larsson"},
      .out = {"wigner_initialized"},

      .in = {"wigner_initialized"},

  };

  wsm_data["Wigner6Init"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initialize the wigner 3 and 6 tables

The default values take about 1 Gb memory.
)--",
      .author = {"Richard Larsson"},
      .out = {"wigner_initialized"},

      .gin = {"fast_wigner_stored_symbols", "largest_wigner_symbol_parameter"},
      .gin_type = {"Index", "Index"},
      .gin_value = {Index{20000000}, Index{250}},
      .gin_desc =
          {R"--(Number of stored symbols possible before replacements)--",
           R"--(Largest symbol used for initializing factorials (e.g., largest J or L))--"},

  };

  wsm_data["Wigner6Unload"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Unloads the wigner 3 and 6 tables
)--",
      .author = {"Richard Larsson"},
      .out = {"wigner_initialized"},

      .in = {"wigner_initialized"},

  };

  wsm_data["WignerFastInfoPrint"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Prints the fast wigner table information if compiled with this option
)--",
      .author = {"Richard Larsson"},

      .in = {"wigner_initialized"},

  };

  wsm_data["WriteBuiltinPartitionFunctionsXML"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Writes all the builtin partition functions to file.

All available partition functions are written to files in the select format
in the select directory

The temperature will be linearly spaced between [Tlow, Tupp] with N values
)--",
      .author = {"Richard Larsson"},

      .in = {"output_file_format"},
      .gin = {"dir", "Tlow", "Tupp", "N"},
      .gin_type = {"String", "Numeric", "Numeric", "Index"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(The directory to write the data towards)--",
                   R"--(The lowest temperature)--",
                   R"--(The highest temperature)--",
                   R"--(The number of temperature points)--"},

  };

  wsm_data["WriteNetCDF"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Writes a workspace variable to a NetCDF file.

This method can write variables of limited groups.

If the filename is omitted, the variable is written
to <basename>.<variable_name>.nc.
)--",
      .author = {"Oliver Lemke"},

      .gin = {"input", "filename"},
      .gin_type =
          {"Vector, Matrix, Tensor3, Tensor4, Tensor5, ArrayOfVector, ArrayOfIndex, ArrayOfMatrix, GasAbsLookup",
           "String"},
      .gin_value = {std::nullopt, String("")},
      .gin_desc = {R"--(Variable to be saved.)--",
                   R"--(Name of the NetCDF file.)--"},

      .pass_names = true};

  wsm_data["WriteNetCDFIndexed"] = WorkspaceMethodInternalRecord{
      .desc = R"--(As *WriteNetCDF*, but creates indexed file names.

This method can write variables of any group.

If the filename is omitted, the variable is written
to <basename>.<variable_name>.nc.
)--",
      .author = {"Oliver Lemke"},

      .in = {"file_index"},
      .gin = {"input", "filename"},
      .gin_type =
          {"Vector, Matrix, Tensor3, Tensor4, Tensor5, ArrayOfVector, ArrayOfMatrix, GasAbsLookup",
           "String"},
      .gin_value = {std::nullopt, String("")},
      .gin_desc = {R"--(Variable to be saved.)--",
                   R"--(Name of the NetCDF file.)--"},

      .pass_names = true};

  wsm_data["WriteXML"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Writes a workspace variable to an XML file.

This method can write variables of any group.

If the filename is omitted, the variable is written
to <basename>.<variable_name>.xml.
If no_clobber is set to 1, an increasing number will be
appended to the filename if the file already exists.
)--",
      .author = {"Oliver Lemke"},

      .in = {"output_file_format"},
      .gin = {"input", "filename", "no_clobber"},
      .gin_type = {"Any", "String", "Index"},
      .gin_value = {std::nullopt, String(""), Index{0}},
      .gin_desc =
          {R"--(Variable to be saved.)--",
           R"--(Name of the XML file.)--",
           R"--(0: Overwrite existing files, 1: Use unique filenames)--"},

      .pass_names = true};

  wsm_data["WriteXMLIndexed"] = WorkspaceMethodInternalRecord{
      .desc = R"--(As *WriteXML*, but creates indexed file names.

The variable is written to a file with name::

  <filename>.<file_index>.xml.

where <file_index> is the value of ``file_index``.

This means that ``filename`` shall here not include the .xml
extension. Omitting filename works as for *WriteXML*.
)--",
      .author = {"Patrick Eriksson, Oliver Lemke"},

      .in = {"output_file_format", "file_index"},
      .gin = {"input", "filename", "digits"},
      .gin_type = {"Any", "String", "Index"},
      .gin_value = {std::nullopt, String(""), Index{0}},
      .gin_desc =
          {R"--(Workspace variable to be saved.)--",
           R"--(File name. See above.)--",
           R"--(Equalize the widths of all numbers by padding with zeros as necessary. 0 means no padding (default).)--"},

      .pass_names = true};

  wsm_data["abs_cia_dataAddCIARecord"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Takes CIARecord as input and appends the results in the appropriate place.

If CIARecord has same species as species in *abs_cia_data*, then the array
position is used to append all of the CIARecord into the array.  If clobber
evaluates as true, cia_record overwrites the appropriate *abs_cia_data*.  If
species in cia_record are not in *abs_cia_data*, the CIARecord is pushed back.
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_cia_data"},

      .in = {"abs_cia_data"},
      .gin = {"cia_record", "clobber"},
      .gin_type = {"CIARecord", "Index"},
      .gin_value = {std::nullopt, Index{0}},
      .gin_desc = {R"--(CIA record to append to *abs_cia_data*.)--",
                   R"--(If true, the new input clobbers the old cia data.)--"},

  };

  wsm_data["abs_cia_dataReadFromCIA"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Read data from a CIA data file for all CIA molecules defined
in *abs_species*.

The units in the HITRAN file are:
 - Frequency: cm^(-1)
 - Binary absorption cross-section: cm^5 molec^(-2)

Upon reading we convert this to the ARTS internal SI units 
of Hz and m^5 molec^(-2).
)--",
      .author = {"Oliver Lemke"},
      .out = {"abs_cia_data"},

      .in = {"abs_species"},
      .gin = {"catalogpath"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Path to the CIA catalog directory.)--"},

  };

  wsm_data["abs_cia_dataReadFromXML"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Read data from a CIA XML file and check that all CIA tags defined
in *abs_species* are present in the file.

The units of the data are described in *abs_cia_dataReadFromCIA*.
)--",
      .author = {"Oliver Lemke"},
      .out = {"abs_cia_data"},

      .in = {"abs_species"},
      .gin = {"filename"},
      .gin_type = {"String"},
      .gin_value = {String("")},
      .gin_desc = {R"--(Name of the XML file.)--"},

  };

  wsm_data["abs_cia_dataReadSpeciesSplitCatalog"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads a species split CIA dataset.
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_cia_data"},

      .in = {"abs_species"},
      .gin = {"basename", "robust"},
      .gin_type = {"String", "Index"},
      .gin_value = {std::nullopt, Index{0}},
      .gin_desc =
          {R"--(The path to the split catalog files)--",
           R"--(Flag to continue in case nothing is found [0 throws, 1 continues])--"},

  };

  wsm_data["abs_hitran_relmat_dataReadHitranRelmatDataAndLines"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Reads HITRAN line mixing data from a basedir

The basedir must point at line mixing data as provided by HITRAN.
The lines will be changed such that ALL CO2 lines are truncated
before adding the HITRAN line mixing lines.

The available modes are such that "VP*" uses Voigt profiles and
"SDVP*" uses speed-dependent Voigt profiles, where the "_Y"
signifies if Rosenkranz-style line mixing is considered or not, and
the "W" at the end signifies that full calculations are used.  At
the line mixing limit, line mixing is simply turned off.

The "FullW" mode uses Lorentzian calculations with the full relaxation
matrix until the line mixing limit is reached and it switches to Voigt.

The HITRAN LM data is available for download at:
https://hitran.org/supplementary/
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_hitran_relmat_data", "abs_lines_per_species"},

          .in = {"abs_lines_per_species", "abs_species"},
          .gin = {"basedir", "linemixinglimit", "fmin", "fmax", "stot", "mode"},
          .gin_type =
              {"String", "Numeric", "Numeric", "Numeric", "Numeric", "String"},
          .gin_value = {std::nullopt,
                        Numeric{-1},
                        Numeric{-1e99},
                        Numeric{1e99},
                        Numeric{0},
                        String("VP_W")},
          .gin_desc =
              {R"--(Direcory where the linemixing data is to be found)--",
               R"--(Line mixing limit as defined by *AbsorptionLines*)--",
               R"--(Minimum frequency to read from)--",
               R"--(Maximum frequency to read until)--",
               R"--(Minimum integrated band strength to consider)--",
               R"--(Mode of calculations.  The options are: "VP", "VP_Y", "SDVP", "SDVP_Y", "FullW", and "VP_W")--"},

      };

  wsm_data["abs_linesAdaptOnTheFlyLineMixing"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adapts the line-catalog from using *ecs_data* data to.
instead fit ordered parameters to imitate the line mxixing

The order should be 1 or 2.  It will compute at 3 as well, but
there's no support in current ARTS LBL to make use of it so it
will crash at some point
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines", "ecs_data"},
      .gin = {"t_grid", "pressure", "order", "robust", "rosenkranz_adaptation"},
      .gin_type = {"Vector", "Numeric", "Index", "Index", "Index"},
      .gin_value =
          {std::nullopt, std::nullopt, std::nullopt, Index{1}, Index{0}},
      .gin_desc =
          {R"--(The sorted temperature grid)--",
           R"--(The pressure at which the adaptation is made)--",
           R"--(The order of the parameters in adaptation)--",
           R"--(Boolean for failed band adaptation behavior. 0: throw exception. not 0: conversion to line-by-line calculations)--",
           R"--(Apply direct Rosenkranz adaptation instead of computing the Eigenvalues)--"},

  };

  wsm_data["abs_linesBaseParameterMatchingLevel"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Set parameter of all levels in *abs_lines* that match with *QuantumIdentifier*.

Only works for these ``parameter_name``:
 - ``"Statistical Weight"``
 - ``"Zeeman Coefficient"``
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"QI", "parameter_name", "change"},
      .gin_type = {"QuantumIdentifier", "String", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Information to match the level.)--",
                   R"--(Name of parameter to be replaced)--",
                   R"--(Value with which to set matching level's value)--"},

  };

  wsm_data["abs_linesBaseParameterMatchingLevels"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(See *abs_linesBaseParameterMatchingLevel*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines"},

          .in = {"abs_lines"},
          .gin = {"QI", "parameter_name", "change"},
          .gin_type = {"ArrayOfQuantumIdentifier", "String", "Vector"},
          .gin_value = {std::nullopt, std::nullopt, std::nullopt},
          .gin_desc = {R"--(Information to match the level.)--",
                       R"--(Name of parameter to be replaced)--",
                       R"--(Value with which to set matching level's value)--"},

      };

  wsm_data["abs_linesBaseParameterMatchingLines"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Set parameter of all lines in *abs_lines* that match with *QuantumIdentifier*.

Only works for these ``parameter_name``:
 - ``"Central Frequency"``
 - ``"Line Strength"``
 - ``"Lower State Energy"``
 - ``"Einstein Coefficient"``
 - ``"Lower Statistical Weight"``
 - ``"Upper Statistical Weight"``
 - ``"Lower Zeeman Coefficient"``
 - ``"Upper Zeeman Coefficient"``
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"QI", "parameter_name", "change"},
      .gin_type = {"QuantumIdentifier", "String", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Information to match the line/band.)--",
                   R"--(Name of parameter to be replaced)--",
                   R"--(Value with which to change matching line's value)--"},

  };

  wsm_data["abs_linesChangeBaseParameterForMatchingLevel"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Change parameter of all levels in *abs_lines* that match with *QuantumIdentifier*.

Only works for these ``parameter_name``:
 - ``"Statistical Weight"``
 - ``"Zeeman Coefficient"``
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines"},

          .in = {"abs_lines"},
          .gin = {"QI", "parameter_name", "change", "relative"},
          .gin_type = {"QuantumIdentifier", "String", "Numeric", "Index"},
          .gin_value = {std::nullopt, std::nullopt, std::nullopt, Index{0}},
          .gin_desc =
              {R"--(Information to match the level.)--",
               R"--(Name of parameter to be replaced)--",
               R"--(Value with which to change matching level's value)--",
               R"--(Flag for relative change (0 is absolute change))--"},

      };

  wsm_data["abs_linesChangeBaseParameterForMatchingLevels"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(See *abs_linesChangeBaseParameterForMatchingLevel*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines"},

          .in = {"abs_lines"},
          .gin = {"QI", "parameter_name", "change", "relative"},
          .gin_type = {"ArrayOfQuantumIdentifier", "String", "Vector", "Index"},
          .gin_value = {std::nullopt, std::nullopt, std::nullopt, Index{0}},
          .gin_desc =
              {R"--(Information to match the level.)--",
               R"--(Name of parameter to be replaced)--",
               R"--(Value with which to change matching level's value)--",
               R"--(Flag for relative change (0 is absolute change))--"},

      };

  wsm_data["abs_linesChangeBaseParameterForMatchingLines"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Change parameter of all lines in *abs_lines* that match with *QuantumIdentifier*.

Only works for these ``parameter_name``:
 - ``"Central Frequency"``
 - ``"Line Strength"``
 - ``"Lower State Energy"``
 - ``"Einstein Coefficient"``
 - ``"Lower Statistical Weight"``
 - ``"Upper Statistical Weight"``
 - ``"Lower Zeeman Coefficient"``
 - ``"Upper Zeeman Coefficient"``
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines"},

          .in = {"abs_lines"},
          .gin = {"QI", "parameter_name", "change", "relative"},
          .gin_type = {"QuantumIdentifier", "String", "Numeric", "Index"},
          .gin_value = {std::nullopt, std::nullopt, std::nullopt, Index{0}},
          .gin_desc =
              {R"--(Information to match the line/band.)--",
               R"--(Name of parameter to be replaced)--",
               R"--(Value with which to change matching line's value)--",
               R"--(Flag for relative change (0 is absolute change))--"},

      };

  wsm_data["abs_linesCompact"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Removes lines that are unimportant because of their
cutoff frequency range
)--",
      .author = {"Stefan Buehler", "Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines", "f_grid"},

  };

  wsm_data["abs_linesCutoff"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets cutoff type and magnitude for all lines.

The line is cut off when this is active at the given frequency.
The only non-zero range is from this range to its negative equivalent

Available ``option``:

- ``"None"``: No cutoff
- ``"ByLine"``: Cutoff relative to a speed-independent shifted line center, highest frequency: F0+cutoff+D0

For "ByLine", the negative frequency is at F0-cutoff-D0
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"option", "value"},
      .gin_type = {"String", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Method of line shape calculations)--",
                   R"--(Value of cutoff)--"},

  };

  wsm_data["abs_linesCutoffMatch"] = WorkspaceMethodInternalRecord{
      .desc = R"--(As *abs_linesCutoff* but for matching bands
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"option", "value", "ID"},
      .gin_type = {"String", "Numeric", "QuantumIdentifier"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Method of line shape calculations)--",
                   R"--(Value of cutoff)--",
                   R"--(ID of one or more bands)--"},

  };

  wsm_data["abs_linesDeleteBadF0"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Deletes all lines in *abs_lines* that have bad central frequencies

If lower evaluates as true, deletes all lines with a frequency below f0.
Otherwise deletes all lines with a frequency above f0.
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"f0", "lower"},
      .gin_type = {"Numeric", "Index"},
      .gin_value = {std::nullopt, Index{1}},
      .gin_desc = {R"--(Target frequency)--",
                   R"--(Lower or upper flag (eval as boolean))--"},

  };

  wsm_data["abs_linesEmptyBroadeningParameters"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets a broadening parameter to empty if it is effectively empty
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},

  };

  wsm_data["abs_linesFlatten"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Makes *abs_lines* with the same ID share lines
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},

  };

  wsm_data["abs_linesKeepBand"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Keep only ``qid``-match band lines in *abs_lines*

Note that other bands are technically kept but have zero lines
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"qid"},
      .gin_type = {"QuantumIdentifier"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Band ID)--"},

  };

  wsm_data["abs_linesLineShapeModelParametersMatchingLines"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Sets line shape model data parameter in matching lines.

The matching is done so that QI must be in the line identifier

Acceptable ``parameter`` (s) are:
 - ``"G0"``
 - ``"D0"``
 - ``"G2"``
 - ``"D2"``
 - ``"FVC"``
 - ``"ETA"``
 - ``"Y"``
 - ``"G"``
 - ``"DV"``

Acceptable ``temperaturemodel`` (s) are:
 - ``"None"``
 - ``"T0"``
 - ``"T1"``
 - ``"T2"``
 - ``"T3"``
 - ``"T4"``
 - ``"T5"``
 - ``"LM_AER"``
 - ``"DPL"``

Acceptable ``species`` are:
 - ``"AIR"`` (so long as it is the broadening species list)
 - ``"SELF"`` (so long as it is the broadening species list)
 - Any species in the line broadening species

See the user guide for the meanings of all of these keywords
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines"},

          .in = {"abs_lines"},
          .gin =
              {"QI", "parameter", "species", "temperaturemodel", "new_values"},
          .gin_type =
              {"QuantumIdentifier", "String", "String", "String", "Vector"},
          .gin_value = {std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::nullopt},
          .gin_desc = {R"--(Information to match the line.)--",
                       R"--(Name of parameter to be replaced)--",
                       R"--(Species of parameter to be changed)--",
                       R"--(Temperature model for the new values)--",
                       R"--(Sets the values found)--"},

      };

  wsm_data["abs_linesLineShapeType"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets shape calculations type for all lines.

Available ``option``:

- ``"DP"``: Doppler profile
- ``"LP"``: Lorentz profile
- ``"VP"``: Voigt profile
- ``"SDVP"``: Speed-dependent Voigt profile
- ``"HTP"``: Hartman-Tran profile

See the theory guide for more details.
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Method of line shape calculations)--"},

  };

  wsm_data["abs_linesLineShapeTypeMatch"] = WorkspaceMethodInternalRecord{
      .desc = R"--(As *abs_linesLineShapeType* but for matching bands
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"option", "ID"},
      .gin_type = {"String", "QuantumIdentifier"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Method of line shape calculations)--",
                   R"--(ID of one or more bands)--"},

  };

  wsm_data["abs_linesLinemixingLimit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets line mixing limit for all lines.

If value is less than 0, no limit is applied and line mixing is active.
Otherwise, line mixing is inactive if the pressure is below the limit.
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"value"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Value of limit)--"},

  };

  wsm_data["abs_linesLinemixingLimitMatch"] = WorkspaceMethodInternalRecord{
      .desc = R"--(See *abs_linesLinemixingLimit* for values

This function only acts on matches between the bands and input ID
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"value", "ID"},
      .gin_type = {"Numeric", "QuantumIdentifier"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Value of limit)--", R"--(ID of one or more bands)--"},

  };

  wsm_data["abs_linesManualMirroring"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Makes a copy of all lines at negative frequency setting themto manual mirroring mode
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},

  };

  wsm_data["abs_linesMirroring"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets mirroring type for all lines.

Available ``option``:

- ``"None"``: No mirrored line
- ``"SameAsLineShape"``: Mirrored line broadened by line shape
- ``"Manual"``: Manually mirrored line (be careful; allows all frequencies)
- ``"Lorentz"``: Mirrored line broadened by Lorentz

Note that mirroring is never applied for DP line shape

Also note that Lorentz profile is approached by most line shapes at high frequency offset.

Also note that Manual settings are potentially dangerous as other frequency
offsets might not work as hoped.
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Method of line mirroring)--"},

  };

  wsm_data["abs_linesMirroringMatch"] = WorkspaceMethodInternalRecord{
      .desc = R"--(As *abs_linesMirroring* but for matching bands
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"option", "ID"},
      .gin_type = {"String", "QuantumIdentifier"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Method of line mirroring)--",
                   R"--(ID of one or more bands)--"},

  };

  wsm_data["abs_linesNormalization"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets normalization type for all lines

Available ``option``:

- ``"VVH"``: Van Vleck and Huber
- ``"VVW"``: Van Vleck and Weisskopf
- ``"RQ"``: Rosenkranz quadratic
- ``"SFS"``: Simple frequency scaling
- ``"None"``: No extra normalization

See the theory guide for more details.
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Method of line normalizations)--"},

  };

  wsm_data["abs_linesNormalizationMatch"] = WorkspaceMethodInternalRecord{
      .desc = R"--(As *abs_linesNormalization* but for matching bands
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"option", "ID"},
      .gin_type = {"String", "QuantumIdentifier"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Method of line normalizations)--",
                   R"--(ID of one or more bands)--"},

  };

  wsm_data["abs_linesPopulation"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets population type for all lines.

Available ``option``:

- ``"LTE"``: Assume band is in LTE
- ``"NLTE"``: Assume band is in NLTE and the upper-to-lower ratio is known
- ``"VibTemps"``: Assume band is in NLTE described by vibrational temperatures and LTE at other levels
- ``"ByHITRANRosenkranzRelmat"``: Assume band needs to compute relaxation matrix to derive HITRAN Y-coefficients
- ``"ByHITRANFullRelmat"``: Assume band needs to compute and directly use the relaxation matrix according to HITRAN
- ``"ByMakarovFullRelmat"``: Assume band needs to compute and directly use the relaxation matrix according to Makarov et al 2020
- ``"ByRovibLinearDipoleLineMixing"``: Assume band needs to compute and directly use the relaxation matrix according to Hartmann, Boulet, Robert, 2008, 1st edition

You must have set ``nlte_field`` and/or its ilk to use the NLTE methods.

You must have *abs_hitran_relmat_data* for the ByHITRANXX methods.

You must have *ecs_data* for the other two relaxation matrix options
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Method of line population)--"},

  };

  wsm_data["abs_linesPopulationMatch"] = WorkspaceMethodInternalRecord{
      .desc = R"--(As *abs_linesPopulation* but for matching bands
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"option", "ID"},
      .gin_type = {"String", "QuantumIdentifier"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Method of line population)--",
                   R"--(ID of one or more bands)--"},

  };

  wsm_data["abs_linesReadSpeciesSplitCatalog"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads a catalog of absorption lines files in a directory
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .gin = {"basename", "robust"},
      .gin_type = {"String", "Index"},
      .gin_value = {std::nullopt, Index{0}},
      .gin_desc =
          {R"--(The path to the split catalog files)--",
           R"--(Flag to continue in case nothing is found [0 throws, 1 continues])--"},

  };

  wsm_data["abs_linesRemoveBand"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Removes ``qid`` band from *abs_lines*
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"qid"},
      .gin_type = {"QuantumIdentifier"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Band ID)--"},

  };

  wsm_data["abs_linesRemoveEmptyBands"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Removes emtpy bands from *abs_lines*
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},

  };

  wsm_data["abs_linesRemoveLines"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Remove lines *abs_lines* outside of specifications

The specifications are:
 - The lower frequency bound (all lines of this frequency or higher may be kept)
 - The upper frequency bound (all lines of this frequency or lower may be kept)
 - The lower intensity bound (all lines with lower intensity may be removed)

If safe evaluates true, all lines in an absorption band must fail the above
tests to be removed

The frequency filtering can be reversed, from keeping upper_frequency to
lower_frequency, to instead remove lines inside the range by setting
``flip_flims`` to 1.

The method *abs_linesRemoveEmptyBands* is internally applied after the
filtering.
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"lower_frequency",
              "upper_frequency",
              "lower_intensity",
              "safe",
              "flip_flims"},
      .gin_type = {"Numeric", "Numeric", "Numeric", "Index", "Index"},
      .gin_value =
          {Numeric{-1e99}, Numeric{1e99}, Numeric{0}, Index{1}, Index{0}},
      .gin_desc =
          {R"--(The lower frequency bound)--",
           R"--(The upper frequency bound)--",
           R"--(The lower intensity bound)--",
           R"--(Remove only lines from a band if all lines of a band fail)--",
           R"--(Reverse the frequecy filtering, see above)--"},

  };

  wsm_data["abs_linesRemoveLinesFromSpecies"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(As *abs_linesRemoveLines* but only for bands of the given species.

``species`` must be a single entry, and must specify the isotopologue
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"species",
              "lower_frequency",
              "upper_frequency",
              "lower_intensity",
              "safe",
              "flip_flims"},
      .gin_type = {"ArrayOfSpeciesTag",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Index",
                   "Index"},
      .gin_value = {std::nullopt,
                    Numeric{-1e99},
                    Numeric{1e99},
                    Numeric{0},
                    Index{1},
                    Index{0}},
      .gin_desc =
          {R"--(Species to be removed)--",
           R"--(The lower frequency bound)--",
           R"--(The upper frequency bound)--",
           R"--(The lower intensity bound)--",
           R"--(Remove only lines from a band if all lines of a band fail)--",
           R"--(Reverse the frequecy filtering)--"},

  };

  wsm_data["abs_linesReplaceBands"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Replace all bands in *abs_lines* that match with bands in ``replacing_bands``.

Each ``replacing_bands`` must match excatly a single band in *abs_lines*.

The matching requires identical quantum number signatures to work.
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"replacing_bands"},
      .gin_type = {"ArrayOfAbsorptionLines"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Line-array that removes lines from *abs_lines*.)--"},

  };

  wsm_data["abs_linesReplaceLines"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Replace all lines in *abs_lines* that match with lines in replacement_lines.

Each replacement_lines must match excatly a single line in *abs_lines*.

The matching requires identical quantum number signatures to work

Note that lines are identified by their quantum number identifier, and if the broadening or
compute data disagree between two bands, a new band is appended unless we can work around the issue.
This may cause *CheckUnique* to fail after running this method
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"replacing_lines"},
      .gin_type = {"ArrayOfAbsorptionLines"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Line-array that replace lines in *abs_lines*.)--"},

  };

  wsm_data["abs_linesSort"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sorts first the lines then the bands by smallest first
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {String("ByFrequency")},
      .gin_desc = {R"--(Sorting option)--"},

  };

  wsm_data["abs_linesT0"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets reference temperature for all lines.
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"value"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Value of T0)--"},

  };

  wsm_data["abs_linesT0Match"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets reference temperature

This function only acts on matches between the bands and input ID
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"value", "ID"},
      .gin_type = {"Numeric", "QuantumIdentifier"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Value of T0)--", R"--(ID of one or more bands)--"},

  };

  wsm_data["abs_linesTurnOffLineMixing"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets all line mixing parameters to emtpy.
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

  };

  wsm_data["abs_linesWriteSpeciesSplitCatalog"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Writes a split catalog, AbsorptionLines by AbsorptionLines.

There will be one unique file generated per AbsorptionLines in abs_lines.

The names of these files will be::

	basename + "." + AbsorptionLines.SpeciesName() + "." + to_string(N) + ".xml"

where N>=0 and the species name is something line "H2O".
)--",
      .author = {"Richard Larsson"},

      .in = {"output_file_format", "abs_lines"},
      .gin = {"basename"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Path to store the files at)--"},

  };

  wsm_data["abs_linesZeemanCoefficients"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the Zeeman coefficients of the lines by user input

The matching is permissive, all in qid must just match.  If there
are multiple matches, the last match rules
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines"},
      .gin = {"qid", "gs"},
      .gin_type = {"ArrayOfQuantumIdentifier", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(Information to match an energy level of a/many lines.)--",
           R"--(Corresponding value to set as Zeeman coefficient)--"},

  };

  wsm_data["abs_lines_per_speciesAdaptHitranLineMixing"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Adapts the line-catalog from using *abs_hitran_relmat_data* to.
instead fit ordered parameters to imitate the line mxixing

The order should be 1 or 2.  It will compute at 3 as well, but
there's no support in current ARTS LBL to make use of it so it
will crash at some point
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species", "abs_hitran_relmat_data"},
          .gin = {"t_grid", "pressure", "order"},
          .gin_type = {"Vector", "Numeric", "Index"},
          .gin_value = {std::nullopt, std::nullopt, std::nullopt},
          .gin_desc = {R"--(The sorted temperature grid)--",
                       R"--(The pressure at which the adaptation is made)--",
                       R"--(The order of the parameters in adaptation)--"},

      };

  wsm_data["abs_lines_per_speciesAdaptOnTheFlyLineMixing"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calls *abs_linesAdaptOnTheFlyLineMixing* for each internal array
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines_per_species"},

      .in = {"abs_lines_per_species", "ecs_data"},
      .gin = {"t_grid", "pressure", "order", "robust", "rosenkranz_adaptation"},
      .gin_type = {"Vector", "Numeric", "Index", "Index", "Index"},
      .gin_value =
          {std::nullopt, std::nullopt, std::nullopt, Index{1}, Index{0}},
      .gin_desc =
          {R"--(The sorted temperature grid)--",
           R"--(The pressure at which the adaptation is made)--",
           R"--(The order of the parameters in adaptation)--",
           R"--(Boolean for failed band adaptation behavior. 0: throw exception. not 0: conversion to line-by-line calculations)--",
           R"--(Apply direct Rosenkranz adaptation instead of computing the Eigenvalues)--"},

  };

  wsm_data["abs_lines_per_speciesBaseParameterMatchingLevel"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(See *abs_linesBaseParameterMatchingLevel*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"QI", "parameter_name", "change"},
          .gin_type = {"QuantumIdentifier", "String", "Numeric"},
          .gin_value = {std::nullopt, std::nullopt, std::nullopt},
          .gin_desc = {R"--(Information to match the level.)--",
                       R"--(Name of parameter to be replaced)--",
                       R"--(Value with which to set matching level's value)--"},

      };

  wsm_data["abs_lines_per_speciesBaseParameterMatchingLevels"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(See *abs_linesBaseParameterMatchingLevel*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"QI", "parameter_name", "change"},
          .gin_type = {"ArrayOfQuantumIdentifier", "String", "Vector"},
          .gin_value = {std::nullopt, std::nullopt, std::nullopt},
          .gin_desc = {R"--(Information to match the level.)--",
                       R"--(Name of parameter to be replaced)--",
                       R"--(Value with which to set matching level's value)--"},

      };

  wsm_data["abs_lines_per_speciesChangeBaseParameterForMatchingLevel"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(See *abs_linesChangeBaseParameterForMatchingLevel*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"QI", "parameter_name", "change", "relative"},
          .gin_type = {"QuantumIdentifier", "String", "Numeric", "Index"},
          .gin_value = {std::nullopt, std::nullopt, std::nullopt, Index{0}},
          .gin_desc =
              {R"--(Information to match the level.)--",
               R"--(Name of parameter to be replaced)--",
               R"--(Value with which to change matching level's value)--",
               R"--(Flag for relative change (0 is absolute change))--"},

      };

  wsm_data["abs_lines_per_speciesChangeBaseParameterForMatchingLevels"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(See *abs_linesChangeBaseParameterForMatchingLevel*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"QI", "parameter_name", "change", "relative"},
          .gin_type = {"ArrayOfQuantumIdentifier", "String", "Vector", "Index"},
          .gin_value = {std::nullopt, std::nullopt, std::nullopt, Index{0}},
          .gin_desc =
              {R"--(Information to match the level.)--",
               R"--(Name of parameter to be replaced)--",
               R"--(Value with which to change matching level's value)--",
               R"--(Flag for relative change (0 is absolute change))--"},

      };

  wsm_data["abs_lines_per_speciesChangeBaseParameterForMatchingLines"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(See *abs_linesChangeBaseParameterForMatchingLines*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"QI", "parameter_name", "change", "relative"},
          .gin_type = {"QuantumIdentifier", "String", "Numeric", "Index"},
          .gin_value = {std::nullopt, std::nullopt, std::nullopt, Index{0}},
          .gin_desc =
              {R"--(Information to match the line/band.)--",
               R"--(Name of parameter to be replaced)--",
               R"--(Value with which to change matching line's value)--",
               R"--(Flag for relative change (0 is absolute change))--"},

      };

  wsm_data["abs_lines_per_speciesChangeBaseParameterForSpecies"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(See *abs_linesChangeBaseParameterForMatchingLines* but for single species
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species", "abs_species"},
          .gin = {"QI", "parameter_name", "change", "relative", "species_tag"},
          .gin_type =
              {"QuantumIdentifier", "String", "Numeric", "Index", "String"},
          .gin_value = {std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        Index{0},
                        std::nullopt},
          .gin_desc =
              {R"--(Information to match the line/band.)--",
               R"--(Name of parameter to be replaced)--",
               R"--(Value with which to change matching line's value)--",
               R"--(Flag for relative change (0 is absolute change))--",
               R"--(The species tag from *abs_species* to change)--"},

      };

  wsm_data["abs_lines_per_speciesCompact"] = WorkspaceMethodInternalRecord{
      .desc = R"--(See *abs_linesCompact*
)--",
      .author = {"Stefan Buehler", "Richard Larsson"},
      .out = {"abs_lines_per_species"},

      .in = {"abs_lines_per_species", "f_grid"},

  };

  wsm_data["abs_lines_per_speciesCreateFromLines"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Split lines up into the different species.

The order of the splitting will match the outer layer of *abs_species*
There will be no respect for the internal layer of *abs_species*
)--",
          .author = {"Stefan Buehler"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines", "abs_species"},

      };

  wsm_data["abs_lines_per_speciesCutoff"] = WorkspaceMethodInternalRecord{
      .desc = R"--(As *abs_linesCutoff* but for *abs_lines_per_species*
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines_per_species"},

      .in = {"abs_lines_per_species"},
      .gin = {"option", "value"},
      .gin_type = {"String", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Method of line shape calculations)--",
                   R"--(Value of cutoff)--"},

  };

  wsm_data["abs_lines_per_speciesCutoffMatch"] = WorkspaceMethodInternalRecord{
      .desc = R"--(As *abs_lines_per_speciesCutoff* but for matching bands
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines_per_species"},

      .in = {"abs_lines_per_species"},
      .gin = {"option", "value", "ID"},
      .gin_type = {"String", "Numeric", "QuantumIdentifier"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Method of line shape calculations)--",
                   R"--(Value of cutoff)--",
                   R"--(ID of one or more bands)--"},

  };

  wsm_data["abs_lines_per_speciesCutoffSpecies"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(As *abs_lines_per_speciesCutoff* but for matching *abs_species*
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines_per_species"},

      .in = {"abs_lines_per_species", "abs_species"},
      .gin = {"option", "value", "species_tag"},
      .gin_type = {"String", "Numeric", "String"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Method of line shape calculations)--",
                   R"--(Value of cutoff)--",
                   R"--(The species tag from *abs_species* to change)--"},

  };

  wsm_data["abs_lines_per_speciesFlatten"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calls *abs_linesFlatten* per internal set of bands
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines_per_species"},

      .in = {"abs_lines_per_species"},

  };

  wsm_data["abs_lines_per_speciesLineShapeModelParametersMatchingLines"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(See *abs_linesLineShapeModelParametersMatchingLines*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin =
              {"QI", "parameter", "species", "temperaturemodel", "new_values"},
          .gin_type =
              {"QuantumIdentifier", "String", "String", "String", "Vector"},
          .gin_value = {std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::nullopt},
          .gin_desc = {R"--(Information to match the line.)--",
                       R"--(Name of parameter to be replaced)--",
                       R"--(Species of parameter to be changed)--",
                       R"--(Temperature model for the new values)--",
                       R"--(Sets the values found)--"},

      };

  wsm_data["abs_lines_per_speciesLineShapeType"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(As *abs_linesLineShapeType* but for *abs_lines_per_species*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"option"},
          .gin_type = {"String"},
          .gin_value = {std::nullopt},
          .gin_desc = {R"--(Method of line shape calculations)--"},

      };

  wsm_data["abs_lines_per_speciesLineShapeTypeMatch"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(As *abs_lines_per_speciesLineShapeType* but for matching bands
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"option", "ID"},
          .gin_type = {"String", "QuantumIdentifier"},
          .gin_value = {std::nullopt, std::nullopt},
          .gin_desc = {R"--(Method of line shape calculations)--",
                       R"--(ID of one or more bands)--"},

      };

  wsm_data["abs_lines_per_speciesLineShapeTypeSpecies"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(As *abs_lines_per_speciesLineShapeType* but for matching *abs_species*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species", "abs_species"},
          .gin = {"option", "species_tag"},
          .gin_type = {"String", "String"},
          .gin_value = {std::nullopt, std::nullopt},
          .gin_desc = {R"--(Method of line shape calculations)--",
                       R"--(The species tag from *abs_species* to change)--"},

      };

  wsm_data["abs_lines_per_speciesLinemixingLimit"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(See *abs_linesLinemixingLimit*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"value"},
          .gin_type = {"Numeric"},
          .gin_value = {std::nullopt},
          .gin_desc = {R"--(Value of limit)--"},

      };

  wsm_data["abs_lines_per_speciesLinemixingLimitMatch"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(See *abs_linesLinemixingLimit* for values

This function only acts on matches between the bands and input ID
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"value", "ID"},
          .gin_type = {"Numeric", "QuantumIdentifier"},
          .gin_value = {std::nullopt, std::nullopt},
          .gin_desc = {R"--(Value of limit)--",
                       R"--(ID of one or more bands)--"},

      };

  wsm_data["abs_lines_per_speciesLinemixingLimitSpecies"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(See *abs_linesLinemixingLimit* but for single species
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species", "abs_species"},
          .gin = {"value", "species_tag"},
          .gin_type = {"Numeric", "String"},
          .gin_value = {std::nullopt, std::nullopt},
          .gin_desc = {R"--(Value of limit)--",
                       R"--(The species tag from *abs_species* to change)--"},

      };

  wsm_data["abs_lines_per_speciesManualMirroring"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Makes a copy of all lines at negative frequency setting them.
to manual mirroring mode
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},

      };

  wsm_data["abs_lines_per_speciesManualMirroringSpecies"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Calls *abs_linesManualMirroring* for given species in *abs_species*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species", "abs_species"},
          .gin = {"species"},
          .gin_type = {"ArrayOfSpeciesTag"},
          .gin_value = {std::nullopt},
          .gin_desc = {R"--(Species to mirror)--"},

      };

  wsm_data["abs_lines_per_speciesMirroring"] = WorkspaceMethodInternalRecord{
      .desc = R"--(As *abs_linesMirroring* but for *abs_lines_per_species*
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines_per_species"},

      .in = {"abs_lines_per_species"},
      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Method of line mirroring)--"},

  };

  wsm_data["abs_lines_per_speciesMirroringMatch"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(As *abs_lines_per_speciesMirroring* but for matching bands
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"option", "ID"},
          .gin_type = {"String", "QuantumIdentifier"},
          .gin_value = {std::nullopt, std::nullopt},
          .gin_desc = {R"--(Method of line mirroring)--",
                       R"--(ID of one or more bands)--"},

      };

  wsm_data["abs_lines_per_speciesMirroringSpecies"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(As *abs_lines_per_speciesMirroring* but for matching *abs_species*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species", "abs_species"},
          .gin = {"option", "species_tag"},
          .gin_type = {"String", "String"},
          .gin_value = {std::nullopt, std::nullopt},
          .gin_desc = {R"--(Method of line mirroring)--",
                       R"--(The species tag from *abs_species* to change)--"},

      };

  wsm_data["abs_lines_per_speciesNormalization"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(As *abs_linesNormalization* but for *abs_lines_per_species*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"option"},
          .gin_type = {"String"},
          .gin_value = {std::nullopt},
          .gin_desc = {R"--(Method of line normalizations)--"},

      };

  wsm_data["abs_lines_per_speciesNormalizationMatch"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(As *abs_lines_per_speciesNormalization* but for matching bands
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"option", "ID"},
          .gin_type = {"String", "QuantumIdentifier"},
          .gin_value = {std::nullopt, std::nullopt},
          .gin_desc = {R"--(Method of line normalizations)--",
                       R"--(ID of one or more bands)--"},

      };

  wsm_data["abs_lines_per_speciesNormalizationSpecies"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(As *abs_lines_per_speciesNormalization* but for matching *abs_species*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species", "abs_species"},
          .gin = {"option", "species_tag"},
          .gin_type = {"String", "String"},
          .gin_value = {std::nullopt, std::nullopt},
          .gin_desc = {R"--(Method of line normalizations)--",
                       R"--(The species tag from *abs_species* to change)--"},

      };

  wsm_data["abs_lines_per_speciesPopulation"] = WorkspaceMethodInternalRecord{
      .desc = R"--(As *abs_linesPopulation* but for *abs_lines_per_species*
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines_per_species"},

      .in = {"abs_lines_per_species"},
      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Method of line population)--"},

  };

  wsm_data["abs_lines_per_speciesPopulationMatch"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(As *abs_lines_per_speciesPopulation* but for matching bands
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"option", "ID"},
          .gin_type = {"String", "QuantumIdentifier"},
          .gin_value = {std::nullopt, std::nullopt},
          .gin_desc = {R"--(Method of line population)--",
                       R"--(ID of one or more bands)--"},

      };

  wsm_data["abs_lines_per_speciesPopulationNlteField"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Turns on NTLE calculations.

Takes the quantum identifers for NLTE temperatures and matches it to
lines in *abs_lines_per_species*.  *abs_species* must be set and is used
to speed up calculations.  After the function is done,  all affected
lines in *abs_lines_per_species* will have an internal tag to the relevant
quantum identifier, which is a requirement for deeper code.

If vibrational_energies is input it must match ``nlte_level_identifiers``
in length.  The vibrational energies of the affected lines will then be
set by the function.  Otherwise, it is assumed the vibrational energies
are set by another method.  If they are not set, calculations will complain
later on while running deeper code.

For now only vibrational energy states are assumed to be able to be in
non-LTE conditions.  The *QuantumIdentifier* for an energy state in ARTS
can look like::

  "CO2-626 EN v1 0/1 v2 1/1 l2 1/1 v3 0/1 r 1/1"

and the matching will match ALL lines with the above.  Note then that if, e.g.,
the "v1 0/1" term was removed from the above, then ARTS will assume that
"v1" is not part of the level of energy state of interest, so lines
of different "v1" will be matched as the same state.  If a line is matched
to more than one energy state, errors should be thrown, but be careful.

Set type of population to change computations and expected input as:

- ``"LTE"``: Compute population by ratios found from LTE temperatures
- ``"TV"``: Compute population by ratios found from NLTE vibrational temperatures
- ``"ND"``: Compute population by ratios found from NLTE number densities
)--",
          .author = {"Richard Larsson"},
          .out = {"nlte_do", "abs_lines_per_species"},

          .in = {"abs_lines_per_species", "atm_field", "nlte_vib_energies"},

      };

  wsm_data["abs_lines_per_speciesPopulationSpecies"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(As *abs_lines_per_speciesPopulation* but for matching *abs_species*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species", "abs_species"},
          .gin = {"option", "species_tag"},
          .gin_type = {"String", "String"},
          .gin_value = {std::nullopt, std::nullopt},
          .gin_desc = {R"--(Method of line population)--",
                       R"--(The species tag from *abs_species* to change)--"},

      };

  wsm_data["abs_lines_per_speciesReadSpeciesSplitCatalog"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(See *abs_linesReadSpeciesSplitCatalog* but only for *abs_species*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_species"},
          .gin = {"basename", "robust"},
          .gin_type = {"String", "Index"},
          .gin_value = {std::nullopt, Index{0}},
          .gin_desc =
              {R"--(The path to the split catalog files)--",
               R"--(Flag to continue in case nothing is found [0 throws, 1 continues])--"},

      };

  wsm_data["abs_lines_per_speciesRemoveLines"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Repeats *abs_linesRemoveLines* for all inner arrays
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines_per_species"},

      .in = {"abs_lines_per_species"},
      .gin = {"lower_frequency",
              "upper_frequency",
              "lower_intensity",
              "safe",
              "flip_flims"},
      .gin_type = {"Numeric", "Numeric", "Numeric", "Index", "Index"},
      .gin_value =
          {Numeric{-1e99}, Numeric{1e99}, Numeric{0}, Index{1}, Index{0}},
      .gin_desc =
          {R"--(The lower frequency bound)--",
           R"--(The upper frequency bound)--",
           R"--(The lower intensity bound)--",
           R"--(Remove only lines from a band if all lines of a band fail)--",
           R"--(Reverse the frequecy filtering)--"},

  };

  wsm_data["abs_lines_per_speciesRemoveLinesFromSpecies"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Repeats *abs_linesRemoveLinesFromSpecies* for all inner arrays
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"species",
                  "lower_frequency",
                  "upper_frequency",
                  "lower_intensity",
                  "safe",
                  "flip_flims"},
          .gin_type = {"ArrayOfSpeciesTag",
                       "Numeric",
                       "Numeric",
                       "Numeric",
                       "Index",
                       "Index"},
          .gin_value = {std::nullopt,
                        Numeric{-1e99},
                        Numeric{1e99},
                        Numeric{0},
                        Index{1},
                        Index{0}},
          .gin_desc =
              {R"--(Species to be removed)--",
               R"--(The lower frequency bound)--",
               R"--(The upper frequency bound)--",
               R"--(The lower intensity bound)--",
               R"--(Remove only lines from a band if all lines of a band fail)--",
               R"--(Reverse the frequecy filtering)--"},

      };

  wsm_data["abs_lines_per_speciesSetEmpty"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Empties *abs_lines_per_species* at the correct size.
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines_per_species"},

      .in = {"abs_species"},

  };

  wsm_data["abs_lines_per_speciesT0"] = WorkspaceMethodInternalRecord{
      .desc = R"--(See *abs_linesT0*
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines_per_species"},

      .in = {"abs_lines_per_species"},
      .gin = {"value"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Value of T0)--"},

  };

  wsm_data["abs_lines_per_speciesT0Match"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets reference temperature

This function only acts on matches between the bands and input ID
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines_per_species"},

      .in = {"abs_lines_per_species"},
      .gin = {"value", "ID"},
      .gin_type = {"Numeric", "QuantumIdentifier"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Value of T0)--", R"--(ID of one or more bands)--"},

  };

  wsm_data["abs_lines_per_speciesT0Species"] = WorkspaceMethodInternalRecord{
      .desc = R"--(See *abs_linesT0* but for single species
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines_per_species"},

      .in = {"abs_lines_per_species", "abs_species"},
      .gin = {"value", "species_tag"},
      .gin_type = {"Numeric", "String"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Value of T0)--",
                   R"--(The species tag from *abs_species* to change)--"},

  };

  wsm_data["abs_lines_per_speciesTurnOffLineMixing"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Sets all line mixing parameters to emtpy.
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

      };

  wsm_data["abs_lines_per_speciesWriteSpeciesSplitCatalog"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(See *abs_linesWriteSpeciesSplitCatalog*

In addition, the structure of the files generated will not care about
generating identifiers for the order in *abs_species*
)--",
          .author = {"Richard Larsson"},

          .in = {"output_file_format", "abs_lines_per_species"},
          .gin = {"basename"},
          .gin_type = {"String"},
          .gin_value = {std::nullopt},
          .gin_desc = {R"--(Path to store the files at)--"},

      };

  wsm_data["abs_lines_per_speciesZeemanCoefficients"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(See *abs_linesZeemanCoefficients*
)--",
          .author = {"Richard Larsson"},
          .out = {"abs_lines_per_species"},

          .in = {"abs_lines_per_species"},
          .gin = {"qid", "gs"},
          .gin_type = {"ArrayOfQuantumIdentifier", "Vector"},
          .gin_value = {std::nullopt, std::nullopt},
          .gin_desc =
              {R"--(Information to match an energy level of a/many lines.)--",
               R"--(Corresponding value to set as Zeeman coefficient)--"},

      };

  wsm_data["abs_lookupAdapt"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Adapts a gas absorption lookup table to the current calculation.

The lookup table can contain more species and more frequencies than
are needed for the current calculation. This method cuts down the
table in memory, so that it contains just what is needed. Also, the
species in the table are brought in the same order as the species in
the current calculation.

Of course, the method also performs quite a lot of checks on the
table. If something is not ok, a runtime error is thrown.

The method sets a flag *abs_lookup_is_adapted* to indicate that the
table has been checked and that it is ok. Never set this by hand,
always use this method to set it!
)--",
      .author = {"Stefan Buehler"},
      .out = {"abs_lookup", "abs_lookup_is_adapted"},

      .in = {"abs_lookup", "abs_species", "f_grid"},

  };

  wsm_data["abs_lookupInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Creates an empty gas absorption lookup table.

This is mainly there to help developers. For example, you can write
the empty table to an XML file, to see the file format.
)--",
      .author = {"Stefan Buehler"},
      .out = {"abs_lookup"},

  };

  wsm_data["abs_speciesAdd"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds species tag groups to the list of absorption species.

This WSM is similar to *abs_speciesSet*, the only difference is that
this method appends species to an existing list of absorption species instead
of creating the whole list.

See *abs_speciesSet* for details on how tags are defined and examples of
how to input them in the control file.
)--",
      .author = {"Stefan Buehler"},
      .out = {"abs_species", "propmat_clearsky_agenda_checked"},

      .in = {"abs_species"},
      .gin = {"species"},
      .gin_type = {"ArrayOfString"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(Specify one String for each tag group that you want to add. Inside the String, separate the tags by commas (plus optional blanks).)--"},

  };

  wsm_data["abs_speciesAdd2"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Adds a species tag group to the list of absorption species and
jacobian quantities.

The method is basically a combined call of *abs_speciesAdd* and
*jacobianAddAbsSpecies*. In this way it is not needed to specify a
tag group in two different places.

Arguments exactly as for *jacobianAddAbsSpecies*. Note that this
method only handles a single tag group, in contrast to
*abs_speciesAdd*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"abs_species",
              "jacobian_quantities",
              "jacobian_agenda",
              "propmat_clearsky_agenda_checked"},

      .in = {"abs_species"},
      .gin = {"gin1", "gin2", "gin3", "species", "unit"},
      .gin_type = {"Vector", "Vector", "Vector", "String", "String"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    String("vmr")},
      .gin_desc = {R"--(Pressure retrieval grid.)--",
                   R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retreival grid.)--",
                   R"--(The species tag of the retrieval quantity.)--",
                   R"--(Retrieval unit. See above.)--"},
      .pass_workspace = true,

  };

  wsm_data["abs_speciesDefineAll"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *abs_species* [i][0] to all species in ARTS
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_species", "propmat_clearsky_agenda_checked"},

  };

  wsm_data["abs_speciesDefineAllInScenario"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Define one tag group for each species known to ARTS and included in an
atmospheric scenario.

You can use this as an alternative to *abs_speciesSet* if you want to make an
absorption calculation that is as complete as possible. The method
goes through all defined species and tries to open the VMR file. If
this works the tag is included, otherwise it is skipped.
)--",
      .author = {"Stefan Buehler"},
      .out = {"abs_species", "propmat_clearsky_agenda_checked"},

      .gin = {"basename"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(The name and path of a particular atmospheric scenario. For example: /pool/lookup2/arts-data/atmosphere/fascod/tropical)--"},

  };

  wsm_data["abs_speciesInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets  *abs_species* to be empty.
)--",
      .author = {"Stefan Buehler"},
      .out = {"abs_species"},

  };

  wsm_data["abs_speciesSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Set up a list of absorption species tag groups.

Workspace variables like *abs_species* contain several tag
groups. Each tag group contains one or more tags. This method converts
descriptions of tag groups given in the keyword to the ARTS internal
representation (an *ArrayOfArrayOfSpeciesTag*). A tag group selects
spectral features which belong to the same species.

A tag is defined in terms of the name of the species, isotopologue, and a
range of frequencies. Species are named after the standard chemical
names, e.g., ``"O3"``. Isotopologues are given by the last digit of the atomic
weight, i.g., ``"O3-668"`` for the asymmetric ozone molecule including an
oxygen 18 atom. Groups of transitions are specified by giving a lower
and upper limit of a frequency range, e.g., ``"O3-666-500e9-501e9"``.

To turn on Zeeman calculation for a species, ``"-Z"`` may be appended
to its name: ``"O2-Z"`` or ``"O2-Z-66"``

The symbol ``"*"`` acts as a wild card. Furthermore, frequency range or
frequency range and isotopologue may be omitted.

Finally, instead of the isotopologue the special letter ``"nl"`` may be given,
e.g., ``"H2O-nl"``. This means that no absorption at all is associated
with this tag. (It is not quite clear if this feature is useful for
anything right now.)

Example:

>>> species = [ "O3-666-500e9-501e9, O3-686", "O3", "H2O-PWR98" ]

   The first tag group selects all O3-666 lines between 500 and
   501 GHz plus all O3-686 lines. 

   The second tag group selects all remaining O3 transitions.

   The third tag group selects H2O, with one of the complete
   absorption models (Rosenkranz 98). No spectrocopic line catalogue
   data will be used for that third tag group.  For more available full
   absorption models see *propmat_clearskyAddPredefined*

   Note that order of tag groups in the species list matters. In our
   example, changing the order of the first two tag group will give
   different results: as ``"O3"`` already selects all O3 transitions,
   no lines will remain to be selected by the
   ``"O3-666-500e9-501e9, O3-686"`` tag.

For CIA species the tag consists of the two involved species and
a dataset index. CIA species can be defined for multiple regions
The dataset index determines which region to use from the corresponding
CIARecord in *abs_cia_data*.

Example

>>> species = [ "N2-CIA-N2-0, N2-CIA-N2-1" ]

For Hitran cross section species the tag consists of the species and
the tagtype XFIT, e.g. CFC11-XFIT. The data for the species must be
available in the *xsec_fit_data* variable.
*propmat_clearsky_agenda_checked* is set to be false.
)--",
      .author = {"Stefan Buehler"},
      .out = {"abs_species", "propmat_clearsky_agenda_checked"},

      .gin = {"species"},
      .gin_type = {"ArrayOfString"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(Specify one String for each tag group that you want to create. Inside the String, separate the tags by commas (plus optional blanks).)--"},

  };

  wsm_data["abs_vecAddGas"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Add gas absorption to first element of absorption vector.

The task of this method is to sum up the gas absorption of the
different gas species and add the result to the first element of the
absorption vector.
)--",
      .author = {"Stefan Buehler"},
      .out = {"abs_vec"},

      .in = {"abs_vec", "propmat_clearsky"},

  };

  wsm_data["antenna_responseGaussian"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets up a Gaussian antenna response.

This method works as *antenna_responseGaussianConstant* but allows
to inlude a frequency variation of the antenna width. Here the FWHM
is specified at a set of frequencies. These frequencies will also be
the frequency grid of *antenna_response*.

If ``grid_width`` is set to <=0, the grid width will be twice the max
value in ``fwhm``.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"antenna_response"},

      .gin = {"f_points", "fwhm", "grid_width", "grid_npoints", "do_2d"},
      .gin_type = {"Vector", "Vector", "Numeric", "Index", "Index"},
      .gin_value =
          {std::nullopt, std::nullopt, Numeric{-1.0}, Index{21}, Index{0}},
      .gin_desc =
          {R"--(Frequencies at which FWHM is defined.)--",
           R"--(Full width at half-maximum of the Gaussian function.)--",
           R"--(Full width of grid (negative value gives 2*fwhm).)--",
           R"--(Number of points to represent the grid, see above.)--",
           R"--(Set to 1 to create a 2D antenna pattern.)--"},

  };

  wsm_data["antenna_responseGaussianConstant"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets up a Gaussian antenna response, with no frequency variation.

The method assumes that the response is the same for all
frequencies and polarisations, and that it can be modelled
as Gaussian. The width of the Gaussian is specified by its
full width at half maximum (FWHM).

The grid generated has ``grid_npoints`` equidistant values, with
the first one at -grid_width/2 and the last one at grid_width/2.

If ``grid_width`` is set to <= 0, a default of twice the FWMH is
applied. This gives a coverage of about 98% of the response.

The default for ``grid_npoints`` is 21. When the grid width is 2*FWHM,
that default value gives an error < 0.001 of the integrated response
using trapezoidal integration. ``grid_npoints`` must be > 1.

If the 2D option is selected (``do_2d``), a circular antenna is
assumed. The same grid and FWHM is applied in both dimensions.

If the grid has a sufficiently high width the integral of the
response is 1. Otherwise the integral is smaller than 1. That
is, no normalisation is applied.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"antenna_response"},

      .gin = {"fwhm", "grid_width", "grid_npoints", "do_2d"},
      .gin_type = {"Numeric", "Numeric", "Index", "Index"},
      .gin_value = {std::nullopt, Numeric{-1.0}, Index{21}, Index{0}},
      .gin_desc =
          {R"--(Full width at half-maximum of the Gaussian function.)--",
           R"--(Full width of grid (negative value gives 2*fwhm).)--",
           R"--(Number of points to represent the grid, see above.)--",
           R"--(Set to 1 to create a 2D antenna pattern.)--"},

  };

  wsm_data["antenna_responseGaussianEffectiveSize"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Sets up Gaussian antenna responses.

Similar to *antenna_responseGaussianConstant* but allows to set up
responses that varies with frequency. That is, the method assumes
that the response is the same for all polarisations, and that it
can be modelled as a Gaussian function varying with frequency.

The full width at half maximum (FWHM in radians) is calculated as::

  fwhm = lambda / leff

where lambda is the wavelength and ``leff`` is the effective size of
the antenna. Normally, ``leff`` is smaller than the physical antenna
size.

Antenna responses are created for ``nf`` frequencies spanning the
range [``fstart``,``fstop``], with a logarithmic spacing. That is, the
frequency grid of the responses is taken from *VectorNLogSpace*.

The responses have a common angular grid. The parameters to define
the grid are the same as for *antenna_responseGaussianConstant*. If
``grid_width`` is <= 0, it is set to twice the FWHM at the lowest
frequency.
)--",
          .author = {"Patrick Eriksson"},
          .out = {"antenna_response"},

          .gin = {"leff",
                  "grid_width",
                  "grid_npoints",
                  "nf",
                  "fstart",
                  "fstop",
                  "do_2d"},
          .gin_type = {"Numeric",
                       "Numeric",
                       "Index",
                       "Index",
                       "Numeric",
                       "Numeric",
                       "Index"},
          .gin_value = {std::nullopt,
                        Numeric{-1.0},
                        Index{21},
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        Index{0}},
          .gin_desc =
              {R"--(Effective size of the antenna,)--",
               R"--(Full width of grid.)--",
               R"--(Number of points to represent the grid.)--",
               R"--(Number of points in frequency grid (must be >= 2))--",
               R"--(Start point of frequency grid)--",
               R"--(End point of frequency grid)--",
               R"--(Set to 1 to create a 2D antenna pattern.)--"},

      };

  wsm_data["atm_fieldAddCustomDataFile"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Add some custom data from file to the atm_field

The key field is used to determine the type of data that is added by input type.

If the input is a String, the data is added to corresponding atmospheric data,
these strings can be

- "t"      - temperature
- "p"      - pressure
- "wind_u" - wind u component
- "wind_v" - wind v component
- "wind_w" - wind w component
- "mag_u"  - mag u component
- "mag_v"  - mag v component
- "mag_w"  - mag w component

If the input is a QuantumIdentifier, it is assumed this is an energy level
identifier for NLTE calculations.

If the input is an ArrayOfSpeciesTag, it is assumed this is for the species
content (VMR, LWC, etc).

The file can contain any of GriddedField3, Tensor3, or Numeric data.  Note
that the method iterates over these using a slow exception-handling routine,
so it would be much more efficient to use either of
*atm_fieldAddGriddedData*, or *atm_fieldAddNumericData* to set the data.
Nevertheless this method is provided to make it easier to compose atmospheric
reading.
)--",
      .author = {"Richard Larsson"},
      .out = {"atm_field"},

      .in = {"atm_field"},
      .gin = {"key", "filename", "extrapolation_type"},
      .gin_type = {"String, ArrayOfSpeciesTag, QuantumIdentifier",
                   "String",
                   "String"},
      .gin_value = {std::nullopt, std::nullopt, String("Nearest")},
      .gin_desc = {R"--(Atmospheric data key.)--",
                   R"--(Filename)--",
                   R"--(Style of extrapolation)--"},

  };

  wsm_data["atm_fieldAddField"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Add another atm_field from file to the current atm_field

The optional flag set_toa determines if the old (default) or
new (if it evaluates as true) atm_field's top of the atmosphere altitude
is used in the output
)--",
      .author = {"Richard Larsson"},
      .out = {"atm_field"},

      .in = {"atm_field"},
      .gin = {"filename", "set_toa"},
      .gin_type = {"String", "Index"},
      .gin_value = {std::nullopt, Index{0}},
      .gin_desc = {R"--(Filename)--",
                   R"--(Flag for overwriting the top of the atmosphere)--"},

  };

  wsm_data["atm_fieldAddGriddedData"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds data to the atm_field

The field must not be regular
)--",
      .author = {"Richard Larsson"},
      .out = {"atm_field"},

      .in = {"atm_field"},
      .gin = {"key", "data", "extrapolation_type"},
      .gin_type = {"String, ArrayOfSpeciesTag, QuantumIdentifier",
                   "GriddedField3",
                   "String"},
      .gin_value = {std::nullopt, std::nullopt, String("Nearest")},
      .gin_desc = {R"--(See *atm_fieldAddCustomDataFile*)--",
                   R"--(Some data)--",
                   R"--(Style of extrapolation)--"},

  };

  wsm_data["atm_fieldAddNumericData"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds data to the atm_field

The field must not be regular
)--",
      .author = {"Richard Larsson"},
      .out = {"atm_field"},

      .in = {"atm_field"},
      .gin = {"key", "data"},
      .gin_type = {"String, ArrayOfSpeciesTag, QuantumIdentifier", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(See *atm_fieldAddCustomDataFile*)--",
                   R"--(Some data)--"},

  };

  wsm_data["atm_fieldIGRF"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Use IGRF to compute the magnetic field at each point

The flag ``parsafe`` exists if you need the calculations to be safe in parallel
computations.
)--",
      .author = {"Richard Larsson"},
      .out = {"atm_field"},

      .in = {"atm_field", "time"},
      .gin = {"parsafe"},
      .gin_type = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc = {R"--(Flag for parallel safety at 3X slowdown cost)--"},

  };

  wsm_data["atm_fieldInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initialize the atmospheric field with some altitude
)--",
      .author = {"Richard Larsson"},
      .out = {"atm_field"},

      .gin = {"toa"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Top of atmosphere altitude [m].)--"},

  };

  wsm_data["atm_fieldLteExternalPartitionFunction"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Turns on NTLE calculations.

Sets NLTE ratios to those expected for LTE calculations
with a known partition function
)--",
          .author = {"Richard Larsson"},
          .out = {"nlte_do", "atm_field", "abs_lines_per_species"},

          .in = {"atm_field", "abs_lines_per_species"},
          .gin = {"nlte_level_identifiers"},
          .gin_type = {"ArrayOfQuantumIdentifier"},
          .gin_value = {std::nullopt},
          .gin_desc = {R"--(List of levels to compute for)--"},

      };

  wsm_data["atm_fieldLteInternalPartitionFunction"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Turns on NTLE calculations.

Sets NLTE ratios to those expected for LTE calculations
with estimation of the partition function as the sum of all
states of a species
)--",
          .author = {"Richard Larsson"},
          .out = {"nlte_do", "atm_field", "abs_lines_per_species"},

          .in = {"atm_field", "abs_lines_per_species"},
          .gin = {"nlte_level_identifiers"},
          .gin_type = {"ArrayOfQuantumIdentifier"},
          .gin_value = {std::nullopt},
          .gin_desc = {R"--(List of levels to compute for)--"},

      };

  wsm_data["atm_fieldRead"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads a new atm_field from a folder or base

There are several indices to indicate what data should be read by the routine
At the end of the routine, a check is performed, throwing a warning if the
data is not useable.

The basename is used to determine a scenario or a folder depending of the last characther.
If the last character is "/", a folder structure is assumed, otherwise a scenario
structure is assumed.  Note that it is only the last characther and that you can
have longer paths still.  Also note that this method respects the internal ARTS
environmental variables to find files not just relative to the execution path.

For instance, if you have a folder structure, you can give basename="atm/".  Now
all the files are expected to be in that folder, e.g., "atm/t.xml" if read_tp is true.
If instead you have a scenario structure you give this as basename="scen".  Now
all the files are expected to belong to that scenario by appedning the names.  For
example, "scen.t.xml" if read_tp is true.

If the flags evaluates true, they expect some files to exist in the basename

- read_tp - ["t.xml", "p.xml", ]
- read_mag - ["mag_u.xml", "mag_v.xml", "mag_w.xml", ]
- read_wind - ["wind_u.xml", "wind_v.xml", "wind_w.xml", ]
- read_specs - [See below]
- read_nlte - "nlte.xml"

If "read_specs" is true, then all the species of *abs_species* are read and the
basename path is expected to contain a file with short-name version for each
unique species.  Some examples:

- abs_species=["H2O-161", "O2-66"], - ["H2O.xml", "O2.xml"]
- abs_species=["H2O-161", "O2-66", "CO2-626"], - ["H2O.xml", "O2.xml", "CO2.xml"]
- abs_species=["H2O-161", "O2-66", "O2-PWR98"], - ["H2O.xml", "O2.xml"]
)--",
      .author = {"Richard Larsson"},
      .out = {"atm_field"},

      .in = {"abs_species"},
      .gin = {"basename",
              "toa",
              "read_tp",
              "read_mag",
              "read_wind",
              "read_specs",
              "read_nlte"},
      .gin_type =
          {"String", "Numeric", "Index", "Index", "Index", "Index", "Index"},
      .gin_value = {String("./"),
                    std::nullopt,
                    Index{1},
                    Index{0},
                    Index{0},
                    Index{1},
                    Index{0}},
      .gin_desc = {R"--(Base for the name of the data files.)--",
                   R"--(Top of atmosphere altitude [m].)--",
                   R"--(Flag to read pressure and temperature)--",
                   R"--(Flag to read magnetic field)--",
                   R"--(Flag to read wind field)--",
                   R"--(Flag to read species)--",
                   R"--(Flag to read NLTE)--"},

  };

  wsm_data["atm_fieldRescalePopulationLevels"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Rescale NLTE field to expected total distribution amongst levels
)--",
      .author = {"Richard Larsson"},
      .out = {"atm_field"},

      .in = {"atm_field"},
      .gin = {"s"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Scaling (e.g., 0.75 for only orth-water on Earth))--"},

  };

  wsm_data["atm_fieldSave"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Saves an atm_field to a folder or base

The output files are split to fit with what *atm_fieldRead* can read, so see it
for most of the filenames that may be generated, depending on the content of the
atm_field of course

Note that there are some exceptions.  If no_clobber is true, the new files will
not overwrite old files.  Also, if there are more than one species with the same
short-name, e.g., "H2O-161" and "H2O-181" both have short-name "H2O", only one of
these will print the "H2O.xml" file whereas the other will print "H2O.2.xml".
The latter is not read by *atm_fieldRead*.  Even worse, we can give no guarantee
at all for whether it is the values from the "H2O-161" or "H2O-181" tags that
give the "H2O.xml" file because the internal data structure is unordered.
)--",
      .author = {"Richard Larsson"},

      .in = {"atm_field"},
      .gin = {"basename", "filetype", "no_clobber"},
      .gin_type = {"String", "String", "Index"},
      .gin_value = {std::nullopt, String("ascii"), Index{0}},
      .gin_desc = {R"--(Base for the name of the data files.)--",
                   R"--(See *WriteXML*)--",
                   R"--(See *WriteXML*)--"},

  };

  wsm_data["atm_fieldTopOfAtmosphere"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the top of the atmosphere altitude to the field
)--",
      .author = {"Richard Larsson"},
      .out = {"atm_field"},

      .in = {"atm_field"},
      .gin = {"toa"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Top of atmosphere altitude [m].)--"},

  };

  wsm_data["atm_fields_compactAddConstant"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds a constant field to atm_fields_compact.

This is handy, e.g., for nitrogen or oxygen. The constant value can
be appended or prepended as an additional field to the already
existing collection of fields. All dimensions (pressure, latitude,
longitude) are filled up, so this works for 1D, 2D, or 3D
atmospheres.

The passed ``name`` of the field has to be in accordance with the
tagging structure described for *atm_fields_compact*.

A list of condensibles can be optionally specified if the VMR of
the added species is assuming dry air. The VMR of the added species
is then scaled down by the sum of the condensibles' VMR::

  VMR * (1 - VMR_sum_of_condensibles).


For Earth this should be set to ["abs_species-H2O"]
)--",
      .author = {"Stefan Buehler, Oliver Lemke"},
      .out = {"atm_fields_compact"},

      .in = {"atm_fields_compact"},
      .gin = {"name", "value", "prepend", "condensibles"},
      .gin_type = {"String", "Numeric", "Index", "ArrayOfString"},
      .gin_value = {std::nullopt, std::nullopt, Index{0}, ArrayOfString{}},
      .gin_desc =
          {R"--(Name of additional atmospheric field, with constant value.)--",
           R"--(Constant value of additional field.)--",
           R"--(0 = Append to the end, 1 = insert at the beginning.)--",
           R"--(List of condensibles used to scale down the VMR of the added species.)--"},

  };

  wsm_data["atm_fields_compactAddSpecies"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds a field to atm_fields_compact, with interpolation.

This method appends or prepends a *GriddedField3* to *atm_fields_compact*.
The *GriddedField3* is interpolated upon the grid of
*atm_fields_compact*. A typical use case for this method may be to
add a climatology of some gas when this gas is needed for radiative
transfer calculations, but not yet present in *atm_fields_compact*.
One case where this happens is when using the Chevalier91L dataset
for infrared simulations.

The grids in *atm_fields_compact* must fully encompass the grids in
the *GriddedField3* to be added, for interpolation to succeed. If
this is not the case, a RuntimeError is thrown.

The passed ``name`` of the field has to be in accordance with the
tagging structure described for *atm_fields_compact*.
)--",
      .author = {"Gerrit Holl"},
      .out = {"atm_fields_compact"},

      .in = {"atm_fields_compact"},
      .gin = {"name", "value", "prepend"},
      .gin_type = {"String", "GriddedField3", "Index"},
      .gin_value = {std::nullopt, std::nullopt, Index{0}},
      .gin_desc =
          {R"--(Name of additional atmospheric field.)--",
           R"--(Value of additional atmospheric field.)--",
           R"--(0 = Append to the end, 1 = insert at the beginning.)--"},

  };

  wsm_data["atm_fields_compactCleanup"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Removes unrealistically small or erroneous data from
*atm_fields_compact* (or other GriddedField4 data)

This WSM checks if the data in *atm_fields_compact* contains
values smaller than the given ``threshold``. In this case, these
values will be set to zero.

The method should be applied if *atm_fields_compact* contains
unrealistically small or erroneous data (NWP/GCM model data
occassionally contains negative values, which are numerical
artefacts rather than physical values.)
)--",
      .author = {"Jana Mendrok"},
      .out = {"atm_fields_compact"},

      .in = {"atm_fields_compact"},
      .gin = {"threshold"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(Threshold below which *atm_fields_compact* values are set to zero.)--"},

  };

  wsm_data["atm_fields_compactCreateFromField"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initiates *atm_fields_compact* from a field.

*atm_fields_compact* will have the same size and grids as the GriddedField3,
but with one dimension as length 1.
)--",
      .author = {"Richard Larsson"},
      .out = {"atm_fields_compact"},

      .gin = {"name", "field"},
      .gin_type = {"String", "GriddedField3"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Name atmospheric field.)--",
                   R"--(The atmospheric field.)--"},

  };

  wsm_data["atm_fields_compactFromMatrix"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets *atm_fields_compact* from 1D fields given in form of a matrix.

For batch calculations it is handy to store atmospheric
profiles in an array of matrix. We take such a matrix, and create
*atm_fields_compact* from it.

The matrix must contain one row for each pressure level.

Not all fields contained in the matrix must be selected into
*atm_fields_compact*, but the selection must at least contain
fields of pressure, temperature, altitude and one absorption
species.
The matrix can contain some additional fields which are not
directly used by ARTS for calculations but can be required for
further processing, e.g. wind speed and direction. These fields do
not need to be transfered into the *atm_fields_compact* variable.

Selection of fields into *atm_fields_compact* works by providing a
field name tag in ``field_names`` for the selected fields, while
unselected fields are tagged by 'ignore'. Order of tags in
``field_names`` is strictly taken as corresponding to column order in
the matrix.
The pressure fields are by convention the first column of the
matrix, hence must not be tagged. That is, there must be given one
field name tag less than matrix columns.

For detailed tagging conventions see *atm_fields_compact*.

Works only for ``atmosphere_dim`` == 1.
)--",
      .author = {"Stefan Buehler", "Daniel Kreyling", "Jana Mendrok"},
      .out = {"atm_fields_compact"},

      .gin = {"gin1", "field_names"},
      .gin_type = {"Matrix", "ArrayOfString"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(One atmosphere matrix from batch input ArrayOfMatrix.)--",
           R"--(Order/names of atmospheric fields.)--"},

  };

  wsm_data["atmfields_checkedCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Checks consistency of (clear sky) atmospheric fields.

The following WSVs are treated: ``p_grid``, ``lat_grid``, ``lon_grid``,
``t_field``, ``vmr_field``, wind_u/v/w_field and mag_u/v/w_field.

If any of the variables above is changed, then this method shall be
called again (no automatic check that this is fulfilled!).

The tests include that:
 1. Atmospheric grids (p/lat/lon_grid) are OK with respect to
    ``atmosphere_dim`` (and vmr_field also regarding *abs_species*).
 2. Atmospheric fields have sizes consistent with the atmospheric
    grids.
 3. *abs_f_interp_order* is not zero if any wind is nonzero.
 4. All values in ``t_field`` are > 0.

Default is that values in ``vmr_field`` are demanded to be >= 0
(ie. zero allowed, in contrast to ``t_field``), but this
requirement can be removed by the ``negative_vmr_ok`` argument.

If any test fails, there is an error. Otherwise,
*atmfields_checked* is set to 1.

The cloudbox is covered by *cloudbox_checked*, ``z_field`` is
part of the checks done around *atmgeom_checked*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"atmfields_checked"},

      .in = {"abs_species", "atm_field"},

  };

  wsm_data["avkCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculate the averaging kernel matrix.

This is done by describing the sensitivity of the
OEM retrieval with respect to the true state of the system. A prerequisite
for the calculation of the averaging kernel matrix is a successful OEM
calculation in which the *jacobian* and the gain matrix *dxdy* have been calculated.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"avk"},

      .in = {"dxdy", "jacobian"},

  };

  wsm_data["backend_channel_responseFlat"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets up a rectangular channel response.

The method assumes that all channels have the same response.

The response of the backend channels is hee assumed to be constant
inside the resolution width, and zero outside.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"backend_channel_response"},

      .gin = {"resolution"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(The spectrometer resolution.)--"},

  };

  wsm_data["backend_channel_responseGaussian"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets up a Gaussian backend channel response.

The method assumes that all channels have the same response.

This method works as *backend_channel_responseGaussianConstant*
but handles the case where the response of each channel must be
described. Here the FWHM is specified for each *f_backend*.

The GINs ``fwhm`` and ``grid_npoints`` work in the same way as for
*antenna_responseGaussianConstant*. A negative ``grid_width``
gives a grid that is twice the FWHM of each channel.
)--",
      .author = {"Patrick Eriksson, Oliver Lemke"},
      .out = {"backend_channel_response"},

      .in = {"f_backend"},
      .gin = {"fwhm", "grid_width", "grid_npoints"},
      .gin_type = {"Vector", "Numeric", "Index"},
      .gin_value = {std::nullopt, Numeric{-1.0}, Index{21}},
      .gin_desc =
          {R"--(Full width at half-maximum of the Gaussian function.)--",
           R"--(Full width of grid.)--",
           R"--(Number of points to represent the grid.)--"},

  };

  wsm_data["backend_channel_responseGaussianConstant"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Sets up a single Gaussian backend channel response.

The method assumes that all channels have the same response.

The GINs ``fwhm`` and ``grid_npoints`` work in the same way as for
*antenna_responseGaussianConstant*. A negative ``grid_width``
gives a grid that is twice the FWHM.
)--",
          .author = {"Patrick Eriksson, Oliver Lemke"},
          .out = {"backend_channel_response"},

          .gin = {"fwhm", "grid_width", "grid_npoints"},
          .gin_type = {"Numeric", "Numeric", "Index"},
          .gin_value = {std::nullopt, Numeric{-1.0}, Index{21}},
          .gin_desc =
              {R"--(Full width at half-maximum of the Gaussian function.)--",
               R"--(Full width of grid.)--",
               R"--(Number of points to represent the grid.)--"},

      };

  wsm_data["background_radFromMatrix"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *background_rad* to matrix input
)--",
      .author = {"Richard Larsson"},
      .out = {"background_rad"},

      .gin = {"iy_mat"},
      .gin_type = {"Matrix"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Background radiation)--"},

  };

  wsm_data["background_transmittanceFromBack"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *background_transmittance* to back of *ppvar_cumtramat*
)--",
      .author = {"Richard Larsson"},
      .out = {"background_transmittance"},

      .in = {"ppvar_cumtramat"},

  };

  wsm_data["background_transmittanceFromFront"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *background_transmittance* to front of *ppvar_cumtramat*
)--",
      .author = {"Richard Larsson"},
      .out = {"background_transmittance"},

      .in = {"ppvar_cumtramat"},

  };

  wsm_data["batch_atm_fields_compactAddConstant"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds a constant field to batch_atm_fields_compact.

Applies *atm_fields_compactAddConstant* to each batch.
The format is equal to that WSM.
)--",
      .author = {"Gerrit Holl"},
      .out = {"batch_atm_fields_compact"},

      .in = {"batch_atm_fields_compact"},
      .gin = {"name", "value", "prepend", "condensibles"},
      .gin_type = {"String", "Numeric", "Index", "ArrayOfString"},
      .gin_value = {std::nullopt, std::nullopt, Index{0}, ArrayOfString{}},
      .gin_desc =
          {R"--(Name of additional atmospheric field, with constant value.)--",
           R"--(Constant value of additional field.)--",
           R"--(0 = Append to the end, 1 = insert at the beginning.)--",
           R"--(List of condensibles used to scale down the VMR of the added species.)--"},

  };

  wsm_data["batch_atm_fields_compactAddSpecies"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Adds a field to *batch_atm_fields_compact*, with interpolation.

This method appends or prepends a *GriddedField3* to each *atm_fields_compact*.
in *batch_atm_fields_compact*. For details, see *atm_fields_compactAddSpecies*.
)--",
      .author = {"Gerrit Holl"},
      .out = {"batch_atm_fields_compact"},

      .in = {"batch_atm_fields_compact"},
      .gin = {"name", "value", "prepend"},
      .gin_type = {"String", "GriddedField3", "Index"},
      .gin_value = {std::nullopt, std::nullopt, Index{0}},
      .gin_desc =
          {R"--(Name of additional atmospheric field. Use, e.g., vmr_ch4 for methane VMR)--",
           R"--(Value of additional atmospheric field.)--",
           R"--(0 = Append to the end, 1 = insert at the beginning.)--"},

  };

  wsm_data["batch_atm_fields_compactCleanup"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Removes unrealistically small or erroneous data from each data field
of *batch_atm_fields_compact* (or other AerrayOfGriddedField4 data)

This WSM checks if the data in *batch_atm_fields_compact* contains
values smaller than the given ``threshold``. In this case, these
values will be set to zero.

The method should be applied if *batch_atm_fields_compact* contains
unrealistically small or erroneous data (NWP/GCM model data
occassionally contains negative values, which are numerical
artefacts rather than physical values.)
)--",
      .author = {"Jana Mendrok"},
      .out = {"batch_atm_fields_compact"},

      .in = {"batch_atm_fields_compact"},
      .gin = {"threshold"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(Threshold below which *atm_fields_compact* values are set to zero.)--"},

  };

  wsm_data["batch_atm_fields_compactFromArrayOfMatrix"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Expand batch of 1D atmospheric state matrices to batch_atm_fields_compact.

This is used to handle 1D batch cases, e.g. from NWP/GCM model like
the Chevallier91L data set, stored in a matrix (it is preferred,
though, to immediatedly store the model fields as
*ArrayOfGriddedField4* and use *ReadXML* to load them directly into
*batch_atm_fields_compact*).

Works only for ``atmosphere_dim`` == 1.

See *atm_fields_compactFromMatrix* for basic documentation.

See *batch_atm_fields_compactAddConstant* and
*batch_atm_fields_compactAddSpecies* for adding additional fields.
)--",
          .author = {"Stefan Buehler", "Daniel Kreyling", "Jana Mendrok"},
          .out = {"batch_atm_fields_compact"},

          .gin = {"atmospheres_fields", "field_names"},
          .gin_type = {"ArrayOfMatrix", "ArrayOfString"},
          .gin_value = {std::nullopt, std::nullopt},
          .gin_desc =
              {R"--(Batch of atmospheres stored in one array of matrix)--",
               R"--(Order/names of atmospheric fields.)--"},

      };

  wsm_data["cloudbox_fieldDisort"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Interface to the DISORT scattering solver (by Stamnes et al.).

DISCLAIMER: There is a couple of known issues with the current
implementation (see below). Use this WSM with care and only if
these limitations/requirements are fulfilled. Results might be
erroneous otherwise.

DISORT provides the radiation field (*cloudbox_field*) from a scalar
1D scattering solution assuming a plane-parallel atmosphere (flat
Earth). Only totally randomly oriented particles are allowed.
Refraction is not taken into account. Only Lambertian surface
reflection is handled.

``nstreams`` is the number of polar angles taken into account
internally in the scattering solution, *za_grid* is the
polar angle grid on which *cloudbox_field* is provided.
``nstreams`` determines the angular resolution, hence the accuracy,
of the scattering solution. The more anisotropic the bulk scattering
matrix, the more streams are required. The computational burden
increases approximately linearly with ``nstreams``. The default value
(8) is sufficient for most microwave scattering calculations. It is
likely insufficient for IR calculations involving ice clouds,
though.

Further, *za_grid* determines the resolution of the output
radiation field. The size of *za_grid* has no practical
impact on computation time in the case of Disort and higher
resolution generally improves the interpolation results, hence
larger *za_grid* are recommended. To ensure sufficient
interpolation accuracy, we require a (hardcoded) minimum size of 38.

Different sphericity levels are emulated here by embedding DISORT
in different ways and using different output. The available options
(from low to high sphericity level) are:

- Cloudbox extends over whole atmosphere (e.g. by setting cloudbox from ``cloudboxSetFullAtm``).
- Cloudbox extends over a limited part of the atmosphere only
  (e.g. by setting cloudbox from ``cloudboxSetAutomatically`` or ``cloudboxSetManually``).
  Internally, DISORT is run over the whole atmosphere, but only the radiation field within
  the cloudbox is passed on and used further in ARTS (e.g. by *yCalc*).

Some auxiliary quantities can be obtained. Auxiliary
quantities are selected by *disort_aux_vars* and returned by *disort_aux*.
Valid choices for auxiliary data are:

- ``"Layer optical thickness"``: Matrix [f_grid, size of p_grid - 1] layer optical thickness.
- ``"Single scattering albedo"``: Matrix [f_grid, size of p_grid - 1] layer single" scattering albedo.
- ``"Direct beam"``: Matrix [f_grid, p_grid]. Attenuated direct at level. Zero, if no sun is present 
)--",
      .author = {"Claudia Emde, Jana Mendrok", "Manfred Brath"},
      .out = {"cloudbox_field", "disort_aux"},

      .in = {"atmfields_checked",
             "atmgeom_checked",
             "scat_data_checked",
             "cloudbox_checked",
             "cloudbox_on",
             "cloudbox_limits",
             "propmat_clearsky_agenda",
             "gas_scattering_agenda",
             "pnd_field",
             "atm_field",
             "surface_field",
             "lat_true",
             "lon_true",
             "abs_species",
             "scat_data",
             "suns",
             "f_grid",
             "za_grid",
             "aa_grid",
             "surface_skin_t",
             "surface_scalar_reflectivity",
             "gas_scattering_do",
             "suns_do",
             "disort_aux_vars"},
      .gin = {"nstreams",
              "Npfct",
              "only_tro",
              "quiet",
              "emission",
              "intensity_correction"},
      .gin_type = {"Index", "Index", "Index", "Index", "Index", "Index"},
      .gin_value =
          {Index{8}, Index{181}, Index{0}, Index{0}, Index{1}, Index{1}},
      .gin_desc =
          {R"--(Number of polar angle directions (streams) in DISORT solution (must be an even number).)--",
           R"--(Number of angular grid points to calculate bulk phase function on (and derive Legendre polynomials from). If <0, the finest za_grid from scat_data will be used.)--",
           R"--(Set to 1 if the scattering data is just of TRO type. Has effect only if Npfct > 3 or Npfct<0, but then leads to much faster calculations.)--",
           R"--(Silence C Disort warnings.)--",
           R"--(Enables blackbody emission. Set to zero, if no Emission e. g. like in visible regime for earth is needed)--",
           R"--(Enables intensity correction. Importantant for low number of streams. Set to zero, if problems encounter or using a high number of streams (>30))--"},
      .pass_workspace = true,

  };

  wsm_data["cloudbox_fieldDisortWithARTSSurface"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Interface to the DISORT scattering solver (by Stamnes et al.).

As *cloudbox_fieldDisort* but uses *surface_rtprop_agenda*.

The Lambertian surface reflection is set by *surface_rtprop_agenda*.
If the GIN inc_angle is inside of the range [0,90], the reflection is
set according to the result of *surface_rtprop_agenda* for this incidence
angle. Otherwise (default) is to call *surface_rtprop_agenda* for
multiple angles, to estimate the hemispheric mean value.

Some auxiliary quantities can be obtained. Auxiliary
quantities are selected by *disort_aux_vars* and returned by *disort_aux*.
Valid choices for auxiliary data are:

- ``"Layer optical thickness"``: Matrix [f_grid, size of p_grid - 1] layer optical thickness.
- ``"Single scattering albedo"``: Matrix [f_grid, size of p_grid - 1] layer single"scattering albedo.
- ``"Direct beam"``: Matrix [f_grid, p_grid]. Attenuated direct at level.Zero, if no sun is present 
)--",
      .author = {"Claudia Emde, Jana Mendrok", "Manfred Brath"},
      .out = {"cloudbox_field", "disort_aux"},

      .in = {"atmfields_checked",
             "atmgeom_checked",
             "scat_data_checked",
             "cloudbox_checked",
             "cloudbox_on",
             "cloudbox_limits",
             "propmat_clearsky_agenda",
             "surface_rtprop_agenda",
             "gas_scattering_agenda",
             "pnd_field",
             "atm_field",
             "surface_field",
             "lat_true",
             "lon_true",
             "abs_species",
             "scat_data",
             "suns",
             "f_grid",
             "za_grid",
             "aa_grid",
             "gas_scattering_do",
             "suns_do",
             "disort_aux_vars"},
      .gin = {"nstreams",
              "Npfct",
              "only_tro",
              "quiet",
              "emission",
              "intensity_correction",
              "inc_angle"},
      .gin_type =
          {"Index", "Index", "Index", "Index", "Index", "Index", "Numeric"},
      .gin_value = {Index{8},
                    Index{181},
                    Index{0},
                    Index{0},
                    Index{1},
                    Index{1},
                    Numeric{-1}},
      .gin_desc =
          {R"--(Number of polar angle directions (streams) in DISORT  solution (must be an even number).)--",
           R"--(Number of angular grid points to calculate bulk phase function on (and derive Legendre polynomials from). If <0, the finest za_grid from scat_data will be used.)--",
           R"--(Set to 1 if the scattering data is just of TRO type. Has effect only if Npfct > 3 or Npfct<0, but then leads to much faster calculations.)--",
           R"--(Silence C Disort warnings.)--",
           R"--(Enables blackbody emission. Set to zero, if no Emission e. g. like in visible regime for earth is needed)--",
           R"--(Enables intensity correction. Importantant for low number of streams. Set to zero, if problems encounter or using a high number of streams (>30))--",
           R"--(Incidence angle, see above.)--"},
      .pass_workspace = true,

  };

  wsm_data["collision_coefficientsFromSplitFiles"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Reads *collision_coefficients* and *collision_line_identifiers* from files.

The species in in these files must match *abs_species*.  The location
must also contain an *ArrayOfQuantumIdentifier* file ending with ``qid.xml``
)--",
          .author = {"Richard Larsson"},
          .out = {"collision_coefficients", "collision_line_identifiers"},

          .in = {"abs_species"},
          .gin = {"basename"},
          .gin_type = {"String"},
          .gin_value = {String("./")},
          .gin_desc = {R"--(path to files to read)--"},

      };

  wsm_data["complex_refr_indexConstant"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Set complex refractive index to a constant value.

Frequency and temperature grids are set to have length 1 (and
set to the value 0).
)--",
      .author = {"Oliver Lemke"},
      .out = {"complex_refr_index"},

      .gin = {"refr_index_real", "refr_index_imag"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Real part of refractive index)--",
                   R"--(Imag part of refractive index)--"},

  };

  wsm_data["complex_refr_indexIceMatzler06"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Refractive index of ice following Matzler06 parameterization.

Calculates temperature dependent complex refractive index of
hexagonal ice at microwave and sub-mm frequencies (10MHz-3Tz).

This parametrization is also applied by the microwave and
submm-wave part of the Warren08 model.

References:
Matzler, C., 2006: Thermal Microwave Radiation: Application for
Remote Sensing, Microwave dielectric properties of ice, pp. 455-462,
Inst. Eng. Technol., Stevenage, U. K.
Warren, S. G., and R. E. Brandt, 2008: Optical constants of ice
from the ultraviolet to the microwave: A revised compilation,
J. Geophys. Res., 113, D14220, doi:10.1029/2007JD009744.
)--",
      .author = {"Jana Mendrok"},
      .out = {"complex_refr_index"},

      .gin = {"data_f_grid", "data_T_grid"},
      .gin_type = {"Vector", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Frequency grid for refractive index calculation)--",
                   R"--(Temperature grid for refractive index calculation)--"},

  };

  wsm_data["complex_refr_indexTemperatureConstant"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Set frequency dependent complex refractive index.

Temperature grid is set to have length 1 (and
set to the value 0).
)--",
      .author = {"Manfred Brath"},
      .out = {"complex_refr_index"},

      .in = {"f_grid"},
      .gin = {"refr_index_real", "refr_index_imag", "temperature"},
      .gin_type = {"Vector", "Vector", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt, Numeric{273.15}},
      .gin_desc =
          {R"--(Real part of refractive index, Dimension [Number of frequencies])--",
           R"--(Imag part of refractive index, Dimension [Number of frequencies])--",
           R"--(Temperature [K])--"},

  };

  wsm_data["complex_refr_indexWaterLiebe93"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Complex refractive index of liquid water according to Liebe 1993.

The method treats liquid water without salt. Thus, not valid below
10 GHz. Upper frequency limit not known, here set to 1000 GHz.
Model parameters taken from Atmlab function epswater93 (by
C. Maetzler), which refer to Liebe 1993 without closer
specifications.

Temperatures must be between -40 and 100 degrees Celsius. The
accuracy of the parametrization below 0 C is not known by us.
)--",
      .author = {"Patrick Eriksson", "Oliver Lemke"},
      .out = {"complex_refr_index"},

      .gin = {"data_f_grid", "data_T_grid"},
      .gin_type = {"Vector", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Frequency grid for refractive index calculation)--",
                   R"--(Temperature grid for refractive index calculation)--"},

  };

  wsm_data["complex_refr_indexWaterVisibleNIRHarvey98"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Refractive index of water and steam for the optical and near infrared.

Refractive index as function of temparature, frequency and density.
It is limited only to the real part. The imaginary part is 0.

From:
Revised formulation for the Refractive Index of Water and Steam as a Function
of Wavelength, Temperature and Density
Journal of Physical and Chemical Reference Data 27, 761 (1998), 
https://doi.org/10.1063/1.556029 27, 761 

See also: http://www.iapws.org/release.html or https://www.nist.gov
Range of validity:

- 271.15K < temperature < 773.15K
- 0 kg m^-3 < density < 1060 kg m^-3
- 157.785504THz < frequency < 1498.96229THz or  0.2µm < wavelength < 1.9µm

Density can be set as Vector of size 1 or it must have the same size as
as data_t_grid.

IMPORTANT: Though the output is *complex_refr_index*, it only contains
the real part. The imaginry part is zero.
)--",
      .author = {"Manfred Brath"},
      .out = {"complex_refr_index"},

      .in = {"complex_refr_index"},
      .gin = {"data_f_grid",
              "data_t_grid",
              "density_water",
              "only_valid_range"},
      .gin_type = {"Vector", "Vector", "Vector", "Index"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, Index{1}},
      .gin_desc =
          {R"--(Frequency grid for refractive index calculation)--",
           R"--(Temperature grid for refractive index calculation)--",
           R"--(Density of water)--",
           R"--(Flag. If true refractive index is calculated only within range of validity and it will throw an error if outside range of validity. If false no check is made, so use at your own risk.)--"},

  };

  wsm_data["covmat1D"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Create 1D covariance matrix.

Creates a 1D covariance matrix for two retrieval quantities on given 
grids from a given functional form. Elements  of the covariance matrix
are computed as::

  S_{i,j} = sigma_i * sigma_j * f(d_{i,j} / l_{i,j})

where d_{i,j} is the distance between the two grid points and l_{i,j}
the mean of the correlation lengths of the grid points.

If a cutoff value ``co`` is given elements with absolute value less than this 
are set to zero.

The following functional forms are available:

- ``"exp"``: f(x) = exp(-x) 
- ``"lin"``: f(x) = 1.0 - x, for x > 1.0, 0.0 otherwise 
- ``"gauss"``: f(x) = exp(-x^2) 
)--",
      .author = {"Simon Pfreundschuh"},

      .gout = {"output"},
      .gout_type = {"Matrix, Sparse"},
      .gout_desc =
          {R"--(The matrix in which to store the covariance matrix.)--"},

      .gin = {"grid_1",
              "grid_2",
              "sigma_1",
              "sigma_2",
              "cls_1",
              "cls_2",
              "co",
              "fname"},
      .gin_type = {"Vector",
                   "Vector",
                   "Vector",
                   "Vector",
                   "Vector",
                   "Vector",
                   "Numeric",
                   "String"},
      .gin_value = {std::nullopt,
                    Vector{},
                    std::nullopt,
                    Vector{},
                    std::nullopt,
                    Vector{},
                    Numeric{0.0},
                    std::nullopt},
      .gin_desc =
          {R"--(The retrieval grid for the first retrieval quantity.)--",
           R"--(The retrieval grid for the second retrieval quantity. (If empty taken as grid_1))--",
           R"--(The variances of the first retrieval quantity.)--",
           R"--(The variances of the second retrieval quantity.(If empty taken as sigma_1))--",
           R"--(The correlations lengths of the first retrieval quantity.)--",
           R"--(The correlations lengths of the second retrieval quantity.(If empty taken as cls_1))--",
           R"--(The cutoff value for covariance matrix elements.)--",
           R"--(The name of the functional form to use.)--"},

  };

  wsm_data["covmat1DMarkov"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Create Markov Process Covariance Matrix.

Create a markov process covariance matrix for a retrieval quantity on
evenly spaced 1D grid. The correlation between two grid points i,j is 
is computed as::

  cov(i,j) = sigma[i] * sigma[j] * exp(- d(i,j) / lc)

where d(i,j) = abs(grid[i] - grid[j]).

This function also sets covmat_inv_block to the analytically computed inverse
of the covariance matrix of the markov provess, which is tri-diagonal. Note
that this requires the retrieval grid to be evenly spaced.
)--",
      .author = {"Simon Pfreundschuh"},

      .gout = {"output", "out_inverse"},
      .gout_type = {"Matrix, Sparse", "Matrix, Sparse"},
      .gout_desc =
          {R"--(The matrix in which to store the covariance matrix.)--",
           R"--(The matrix in which to store the inverse of the covariance matrix.)--"},

      .gin = {"grid", "sigma", "lc", "co"},
      .gin_type = {"Vector", "Vector", "Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, Numeric{0.0}},
      .gin_desc =
          {R"--(The retrieval grid.)--",
           R"--(The vairance for each grid point.)--",
           R"--(The correlation length of the Markov process.)--",
           R"--(The cutoff value below which elements will be set to 0.0)--"},

  };

  wsm_data["covmatDiagonal"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets the matrix in covmat_block to a diagonal matrix with the variances
provided in ``vars`` as diagonal elements.
Also sets covmat_block_inv to the inverse of the block so that the
computation of the inverse is avoided.
)--",
      .author = {"Simon Pfreundschuh"},

      .gout = {"output", "out_inverse"},
      .gout_type = {"Matrix, Sparse", "Matrix, Sparse"},
      .gout_desc =
          {R"--(The matrix in which to store the covariance matrix.)--",
           R"--(The matrix in which to store the inverse of the covariance matrix.)--"},

      .gin = {"vars"},
      .gin_type = {"Vector"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(Variances to be used as diagonal elements of covmat_block.)--"},

  };

  wsm_data["covmat_seAddBlock"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Add a block to the measurement covariance matrix *covmat_se*

This functions adds a given dense or sparse matrix as block to the covariance
matrix *covmat_sx*. The position of the block can be given by the generic
arguments ``i`` and ``j``. Note that diagonal blocks must be added in order starting from
in the top left corner. If an off-diagonal block is added it must have corresponding
existing blocks on the diagonal and these must be consistent with the dimensions
of the block.  If ``i`` and ``j``  are not provided, the blok will be added
at the first free spot on the diagonal.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_se"},

      .in = {"covmat_se"},
      .gin = {"block", "i", "j"},
      .gin_type = {"Matrix, Sparse", "Index", "Index"},
      .gin_value = {std::nullopt, Index{-1}, Index{-1}},
      .gin_desc =
          {R"--(The block to add to the covariance matrix)--",
           R"--(Index of a retrieval quantity. Must satisfy ``i`` <= ``j``.)--",
           R"--(Index of a retrieval quantity. Must satisfy ``i`` <= ``j``.)--"},

  };

  wsm_data["covmat_seAddInverseBlock"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Add the inverse of a block to covariance matrix *covmat_se*

This functions adds a given matrix as the inverse of a block in the covariance
matrix *covmat_se*. The purpose of this function is to allow the user to
to use a precomputed inverse for this block in the covariance matrix, that may
for example have been obtained analytically.

This function requires the corresponding non-inverse block to already be present
in *covmat_se*

If the 'i' and 'j' input arguments are not given, the inverse block
will be added at the position of the most recently added non-inverse diagonal
block.


Note that for this to work this retrieval quantity must be independent from
other retrieval quantities that do not have an inverse. Otherwise the inverse
will be ignored and recomputed numerically.

For the rest, the same requirements as for *covmat_seAddBlock* apply.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_se"},

      .in = {"covmat_se"},
      .gin = {"block", "i", "j"},
      .gin_type = {"Matrix, Sparse", "Index", "Index"},
      .gin_value = {std::nullopt, Index{-1}, Index{-1}},
      .gin_desc =
          {R"--(The inverse block to add to the covariance matrix)--",
           R"--(Index of a retrieval quantity. Must satisfy ``i`` <= ``j``.)--",
           R"--(Index of a retrieval quantity. Must satisfy ``i`` <= ``j``.)--"},

  };

  wsm_data["covmat_seSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Set covmat_se to a given matrix.

This sets the measurement covariance matrix *covmat_se* to
the matrix given by the generic input ``covmat``. The covariance
matrix can be of type CovarianceMatrix, Matrix or Sparse.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_se"},

      .gin = {"covmat"},
      .gin_type = {"CovarianceMatrix, Matrix, Sparse"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(The matrix to set as the covariance matrix.)--"},

  };

  wsm_data["covmat_soCalc"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculates the covariance matrix describing the error due to uncertainties
in the observation system.

The uncertainties of the observation system are
described by *covmat_se*, which must be set by the user to include the
relevant contributions from the measurement and the forward model.

Prerequisite for the calculation of *covmat_so* is a successful OEM
computation where also the gain matrix has been computed.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_so"},

      .in = {"dxdy", "covmat_se"},

  };

  wsm_data["covmat_ssCalc"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculates the covariance matrix describing the error due to smoothing.

The calculation of *covmat_ss* also requires the averaging kernel matrix *avk*
to be computed after a successful OEM calculation.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_ss"},

      .in = {"avk", "covmat_sx"},

  };

  wsm_data["covmat_sxAddBlock"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Add a block to the a priori covariance matrix *covmat_sx*

This functions adds a given matrix as a block in the covariance
matrix *covmat_sx*. The position of the block can be given by the generic
arguments ``i`` and ``j``, which should give the index of the retrieval quantity in
*jacobian_quantities*, which is given just by the order the quantities have been
added to the retrieval.

If arguments ``i`` and ``j`` are omitted, the block will be added as diagonal block
for the last added retrieval quantity.

If provided, the index ``i`` must be less than or equal to ``j``. Also the provided
block must be consistent with the corresponding retrieval quantities.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx"},

      .in = {"covmat_sx", "jacobian_quantities"},
      .gin = {"block", "i", "j"},
      .gin_type = {"Matrix, Sparse", "Index", "Index"},
      .gin_value = {std::nullopt, Index{-1}, Index{-1}},
      .gin_desc =
          {R"--(The block to add to the covariance matrix)--",
           R"--(Index of a retrieval quantity. Must satisfy ``i`` <= ``j``.)--",
           R"--(Index of a retrieval quantity. Must satisfy ``i`` <= ``j``.)--"},

  };

  wsm_data["covmat_sxAddInverseBlock"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Add the inverse of a block in covariance matrix *covmat_sx*

This functions adds a given matrix as the inverse of a block in the covariance
matrix *covmat_sx*. The purpose of this function is to allow the user to
to use a precomputed inverse for this block in the covariance matrix, the may
for example by obtained analytically.

This function requires the non-inverse block to already be present in *covmat_sx*
If the 'i' and 'j' input arguments are not given, the inverse block
will be added at the position of the most recently added non-inverse diagonal
block.

Note that for this to work this retrieval quantity must be independent from
other retrieval quantities that do not have an inverse. Otherwise the inverse
will be ignored and recomputed numerically.

For the rest, the same requirements as for *covmat_sxAddBlock* apply.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx"},

      .in = {"covmat_sx", "jacobian_quantities"},
      .gin = {"block", "i", "j"},
      .gin_type = {"Matrix, Sparse", "Index", "Index"},
      .gin_value = {std::nullopt, Index{-1}, Index{-1}},
      .gin_desc =
          {R"--(The inverse block to add to the covariance matrix)--",
           R"--(Index of a retrieval quantity. Must satisfy ``i`` <= ``j``.)--",
           R"--(Index of a retrieval quantity. Must satisfy ``i`` <= ``j``.)--"},

  };

  wsm_data["covmat_sxExtractSqrtDiagonal"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Extract the square root of the diagonal of the state space covariance matrix.

This function extracts the diagonal of the state space covariance matrix
*covmat_sx* and computes its square root. The resulting vector can then
be used as ``x_norm`` argument for the OEM method to avoid scaling problems.
)--",
      .author = {"Simon Pfreundschuh"},

      .gout = {"x_norm"},
      .gout_type = {"Vector"},
      .gout_desc =
          {R"--(The vector containing the square root of the diagonal elements of *covmat_sx*)--"},
      .in = {"covmat_sx"},

  };

  wsm_data["covmat_sxSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Set covmat_sx to a given matrix.

This sets the measurement covariance matrix *covmat_sx* to
the matrix given by the generic input ``covmat``. The covariance
matrix can be of type CovarianceMatrix, Matrix or Sparse.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx"},

      .gin = {"covmat"},
      .gin_type = {"CovarianceMatrix, Matrix, Sparse"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(The matrix to set as the covariance matrix.)--"},

  };

  wsm_data["diameter_maxFromDiameter_volume_equ"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates maximum and area equivalent diameters from volume
equivalent diameter.

This is primarily a help function for using the T-matrix method
and only a few particle shapes are handled. 
For shapes handled and further comments on the input arguments, see
*scat_data_singleTmatrix*.

Area equivalent diameter is the equivalent sphere diameter
corresponding to the "maximum axial area". This is the largest
cross-sectional area of the particle, observed either along the
particle's main axis or in the perpendicular direction. That is,
for a cylinder having diameter d and thickness h, this area is
either (pi*d^2)/4 or (h*d).
)--",
      .author = {"Johan Strandgren", "Patrick Eriksson"},

      .gout = {"diameter_max", "diameter_area_equ"},
      .gout_type = {"Numeric", "Numeric"},
      .gout_desc =
          {R"--(Maximum dimension of the particle.)--",
           R"--(Maximum axial area equivalent diameter of the particle, see above.)--"},

      .gin = {"shape", "diameter_volume_equ", "aspect_ratio"},
      .gin_type = {"String", "Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Particle shape.)--",
                   R"--(Particle equivalent volume diameter.)--",
                   R"--(Particle aspect ratio.)--"},

  };

  wsm_data["diameter_volume_equFromDiameter_max"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Converts from maximum to volume equivalent diameter.

This is primarily a help function for using the T-matrix part
and only a few particle shapes are handled. 
For shapes handled and further comments on the input arguments,
see *scat_data_singleTmatrix*.

Also the volume is provided. It is simply sqrt(pi*dveq^3/6).
)--",
          .author = {"Johan Strandgren", "Patrick Eriksson"},

          .gout = {"diameter_volume_equ", "volume"},
          .gout_type = {"Numeric", "Numeric"},
          .gout_desc = {R"--(Particle volume equivalent diameter.)--",
                        R"--(Volume of the particle.)--"},

          .gin = {"shape", "diameter_max", "aspect_ratio"},
          .gin_type = {"String", "Numeric", "Numeric"},
          .gin_value = {std::nullopt, std::nullopt, std::nullopt},
          .gin_desc = {R"--(Particle shape.)--",
                       R"--(Maximum dimension of the particle.)--",
                       R"--(Particle aspect ratio.)--"},

      };

  wsm_data["diy_dxTransform"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Transforms *diy_dpath* and adds it to *diy_dx*.
)--",
      .author = {"Richard Larsson"},
      .out = {"diy_dx", "diy_dpath"},

      .in = {"diy_dx",
             "diy_dpath",
             "ppath",
             "ppvar_atm",
             "abs_species",
             "iy_transmittance",
             "water_p_eq_agenda",
             "jacobian_quantities",
             "jacobian_do",
             "iy_agenda_call1"},

      .pass_workspace = true,

  };

  wsm_data["dlosDiffOfLos"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Derives the difference betwenn zenith and azimuth angles.

Determines the difference between a set of angles (``other_los``)
and a reference direction (``ref_los``). This method reverses the
addition made by *sensor_losAddLosAndDlos*.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"gdlos"},
      .gout_type = {"Matrix"},
      .gout_desc = {R"--(Derived differences in line-of-sight.)--"},

      .gin = {"ref_los", "other_los"},
      .gin_type = {"Vector", "Matrix"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Reference line-of-sight (a single los).)--",
                   R"--(Other line-of-sights (can be multiple los).)--"},

  };

  wsm_data["dlosGauss"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gives a *dlos* suitable for a circular Gaussian response.

The method generates a *dlos* where each direction is meant to have
an equal weight in terms of the product of solid angle and a circular
Gaussian response. That is, the FWHM of the response is equal in zenith
and azimuth directions.

The points have an unequal distribution in radius. The weight in radius
equals radius times the magnitude of the Gaussian response (this product
peaks for a radius around 0.41 * FWHM). The points are distributed in
polar angle simply by adding 208.8 deg from one point to next. There is
no theoretical basis for this step in angle, just found to result in a
relatively uniform distribution over the circle.

The method should mainly be used for ``npoints`` above 10-20. For lower
``npoints``, a rectangular pattern should give a more robust sampling
spatially.

Default is to let *dlos_weight_vector* represent the solid angle of
each dlos direction. With ``include_response_in_weight`` set to 1, all
elements of *dlos_weight_vector* are equal and their sum is 1.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"dlos", "dlos_weight_vector"},

      .gin = {"fwhm", "npoints", "include_response_in_weight"},
      .gin_type = {"Numeric", "Index", "Index"},
      .gin_value = {std::nullopt, std::nullopt, Index{0}},
      .gin_desc =
          {R"--(The full width at half maximum of the Gaussian response.)--",
           R"--(Number of dlos-directions.)--",
           R"--(Set to 1 to include the response values in *dlos_weight_vector*.)--"},

  };

  wsm_data["dlosUniform"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gives *dlos* a rectangular coverage, with uniform spacing.

The directions described by *dlos* are uniform with respect
to relative zenith and azimuth (and thus are NOT uniform in
solid angle). The same angular grid is applied in both angular
dimensions. With width = 1 and npoints = 5, the angular grids
both are [-0.4, -0.2, 0, 0.2, 0.4].

The inner loop in is the zenith direction. That is, first comes
all relative zenith angles for first relative azimuth angle etc.

For default settings, the resulting number of dlos-directions
is npoints * npoints.

If GIN ``crop_circular`` is true, dlos-es at a radius outside of
width/2 are removed. The resulting number of directions then
approaches pi * npoints * npoints / 4, for high values of ``npoints``.
There is no effect of ``crop_circular`` for npoints=2, while for
npoints=3 the corner points are removed (despite being inside
the radius limit) and the number of directions becomes five.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"dlos", "dlos_weight_vector"},

      .gin = {"width", "npoints", "crop_circular"},
      .gin_type = {"Numeric", "Index", "Index"},
      .gin_value = {std::nullopt, std::nullopt, Index{0}},
      .gin_desc =
          {R"--(The full width, in each dimension, in degrees.)--",
           R"--(Number of points over the width, in each dimension (>1).)--",
           R"--(Set to 1, to crop dlos-es to obtain a pseudo-circular pattern.)--"},

  };

  wsm_data["dobatch_calc_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *dobatch_calc_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--",
      .author = {"Richard Larsson"},
      .out = {"dobatch_calc_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["doit_conv_flagAbs"] = WorkspaceMethodInternalRecord{
      .desc = R"--(DOIT convergence test (maximum absolute difference).

The function calculates the absolute differences for two successive
iteration fields. It picks out the maximum values for each Stokes
component separately. The convergence test is fullfilled under the
following conditions:

- ``|I(m+1) - I(m)| < epsilon_1``: Intensity.
- ``|Q(m+1) - Q(m)| < epsilon_2``: The other Stokes components.
- ``|U(m+1) - U(m)| < epsilon_3``: 
- ``|V(m+1) - V(m)| < epsilon_4``: 

These conditions have to be valid for all positions in the
cloudbox and for all directions.
)--",
      .author = {"Claudia Emde"},
      .out = {"doit_conv_flag",
              "doit_iteration_counter",
              "cloudbox_field_mono"},

      .in = {"doit_conv_flag",
             "doit_iteration_counter",
             "cloudbox_field_mono",
             "cloudbox_field_mono_old"},
      .gin = {"epsilon", "max_iterations", "nonconv_return_nan"},
      .gin_type = {"Vector", "Index", "Index"},
      .gin_value = {std::nullopt, Index{100}, Index{0}},
      .gin_desc =
          {R"--(Limits for convergence. A vector with length matching ``stokes_dim`` with unit [W / (m^2 Hz sr)].)--",
           R"--(Maximum number of iterations allowed to reach convergencelimit.)--",
           R"--(Flag whether to accept result at max_iterations (0=default)or whether to return NaNs in case of non-convergence atmax_iterations)--"},

  };

  wsm_data["doit_conv_flagAbsBT"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(DOIT convergence test (maximum absolute difference in Rayleigh Jeans BT)

As *doit_conv_flagAbs* but convergence limits are specified in
Rayleigh-Jeans brighntess temperatures.
)--",
      .author = {"Sreerekha T.R.", "Claudia Emde"},
      .out = {"doit_conv_flag",
              "doit_iteration_counter",
              "cloudbox_field_mono"},

      .in = {"doit_conv_flag",
             "doit_iteration_counter",
             "cloudbox_field_mono",
             "cloudbox_field_mono_old",
             "f_grid",
             "f_index"},
      .gin = {"epsilon", "max_iterations", "nonconv_return_nan"},
      .gin_type = {"Vector", "Index", "Index"},
      .gin_value = {std::nullopt, Index{100}, Index{0}},
      .gin_desc =
          {R"--(Limits for convergence. A vector with length matching ``stokes_dim`` with unit [K].)--",
           R"--(Maximum number of iterations allowed to reach convergencelimit.)--",
           R"--(Flag whether to accept result at max_iterations (0=default)or whether to return NaNs in case of non-convergence atmax_iterations)--"},

  };

  wsm_data["doit_conv_flagLsq"] = WorkspaceMethodInternalRecord{
      .desc = R"--(DOIT convergence test (least squares).

As *doit_conv_flagAbsBT* but applies a least squares convergence
test between two successive iteration fields.

Warning: This method is not recommended because this kind of
convergence test is not sufficiently strict, so that the
DOIT result might be wrong.
)--",
      .author = {"Claudia Emde"},
      .out = {"doit_conv_flag",
              "doit_iteration_counter",
              "cloudbox_field_mono"},

      .in = {"doit_conv_flag",
             "doit_iteration_counter",
             "cloudbox_field_mono",
             "cloudbox_field_mono_old",
             "f_grid",
             "f_index"},
      .gin = {"epsilon", "max_iterations", "nonconv_return_nan"},
      .gin_type = {"Vector", "Index", "Index"},
      .gin_value = {std::nullopt, Index{100}, Index{0}},
      .gin_desc =
          {R"--(Limits for convergence. A vector with length matching ``stokes_dim`` with unit [K].)--",
           R"--(Maximum number of iterations allowed to reach convergencelimit.)--",
           R"--(Flag whether to accept result at max_iterations (0=default)or whether to return NaNs in case of non-convergence atmax_iterations)--"},

  };

  wsm_data["doit_conv_test_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *doit_conv_test_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--",
      .author = {"Richard Larsson"},
      .out = {"doit_conv_test_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["doit_mono_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *doit_mono_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--",
      .author = {"Richard Larsson"},
      .out = {"doit_mono_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["doit_rte_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *doit_rte_agenda* to a default value

Options are:

- There are currently no options, calling this function is an error.
)--",
      .author = {"Richard Larsson"},
      .out = {"doit_rte_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["doit_scat_field_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *doit_scat_field_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--",
      .author = {"Richard Larsson"},
      .out = {"doit_scat_field_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["doit_za_grid_optCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Zenith angle grid optimization for scattering calculation.

This method optimizes the zenith angle grid. As input it requires
a radiation field (*cloudbox_field*) which is calculated on a very
fine zenith angle grid (*za_grid*). Based on this field
zenith angle grid points are selected, such that the maximum
difference between the radiation field represented on the very
fine zenith angle grid and the radiation field represented on the
optimized grid (*doit_za_grid_opt*) is less than the accuracy
(``acc``). Between the grid points the radiation field is interpolated
linearly or polynomially depending on *doit_za_interp*.

Note: The method works only for a 1D atmosphere and for one
frequency.
)--",
      .author = {"Claudia Emde"},
      .out = {"doit_za_grid_opt"},

      .in = {"cloudbox_field_mono", "za_grid", "doit_za_interp"},
      .gin = {"acc"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Accuracy to achieve [%].)--"},

  };

  wsm_data["doit_za_interpSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Define interpolation method for zenith angle dimension.

You can use this method to choose the interpolation method for
interpolations in the zenith angle dimension.
)--",
      .author = {"Claudia Emde"},
      .out = {"doit_za_interp"},

      .gin = {"interp_method"},
      .gin_type = {"String"},
      .gin_value = {String("linear")},
      .gin_desc = {R"--(Interpolation method ("linear" or "polynomial").)--"},

  };

  wsm_data["ecs_dataAddMakarov2020"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the O2-66 microwave band data for ECS.
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data"},

      .in = {"ecs_data", "isotopologue_ratios"},

  };

  wsm_data["ecs_dataAddMeanAir"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets ECS data for air from other data if available.
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data"},

      .in = {"ecs_data"},
      .gin = {"vmrs", "specs"},
      .gin_type = {"Vector", "ArrayOfSpeciesTag"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(VMRs of air species)--", R"--(Air species)--"},

  };

  wsm_data["ecs_dataAddRodrigues1997"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the CO2-626, CO2-636, and CO2-628 IR band data for ECS.

Note that the broadening species has to be N2 and not AIR for the band,
and that N2 VMR must be present
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data"},

      .in = {"ecs_data", "isotopologue_ratios"},

  };

  wsm_data["ecs_dataAddSpeciesData"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets ECS data for one set of species and quantum identifiers.
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data"},

      .in = {"ecs_data", "isotopologue_ratios"},
      .gin = {"qid",
              "species",
              "scaling_type",
              "scaling",
              "beta_type",
              "beta",
              "lambda_type",
              "lambda",
              "collisional_distance_type",
              "collisional_distance"},
      .gin_type = {"QuantumIdentifier",
                   "String",
                   "String",
                   "Vector",
                   "String",
                   "Vector",
                   "String",
                   "Vector",
                   "String",
                   "Vector"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    String("T0"),
                    std::nullopt,
                    String("T0"),
                    std::nullopt,
                    String("T0"),
                    std::nullopt,
                    String("T0"),
                    std::nullopt},
      .gin_desc =
          {R"--(Band identifier)--",
           R"--(Species identifier)--",
           R"--(Temperature model for the main scaling coefficients for Q)--",
           R"--(Main scaling coefficients for Q)--",
           R"--(Temperature model for the energy scaling coefficient for Q)--",
           R"--(Energy scaling coefficient for Q)--",
           R"--(Temperature model for the energy exponent for Q)--",
           R"--(Energy exponent for Q)--",
           R"--(Temperature model for the mean collision interaction distance)--",
           R"--(Mean collision interaction distance)--"},

  };

  wsm_data["ecs_dataAddTran2006"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the O2-66 visible band data for ECS.
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data"},

      .in = {"ecs_data", "isotopologue_ratios"},

  };

  wsm_data["ecs_dataAddTran2011"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the CO2-626, CO2-636, and CO2-628 IR band data for ECS.
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data"},

      .in = {"ecs_data", "isotopologue_ratios"},

  };

  wsm_data["ecs_dataInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Resets/initializes the ECS data.
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data"},

  };

  wsm_data["ext_matAddGas"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Add gas absorption to all diagonal elements of extinction matrix.

The task of this method is to sum up the gas absorption of the
different gas species and add the result to the extinction matrix.
)--",
      .author = {"Stefan Buehler"},
      .out = {"ext_mat"},

      .in = {"ext_mat", "propmat_clearsky"},

  };

  wsm_data["f_gridFromAbsorptionLines"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *f_grid* to a grid relative to *abs_lines_per_species*

Each line will have *abs_lines_per_species* will have a grid
of ``num_freqs`` grid points in [ f0 + ``delta_f_low``, f0 + ``delta_f_upp`` ],
where f0 is the line center.

Before leaving the function, *f_grid* is sorted.

Note that this method could generate significantly large *f_grid*
if used carelessly
)--",
      .author = {"Richard Larsson"},
      .out = {"f_grid"},

      .in = {"abs_lines_per_species"},
      .gin = {"delta_f_low", "delta_f_upp", "num_freqs"},
      .gin_type = {"Numeric", "Numeric", "Index"},
      .gin_value = {Numeric{-5e6}, Numeric{5e6}, std::nullopt},
      .gin_desc = {R"--(Lower range of delta f)--",
                   R"--(Upper range of delta f)--",
                   R"--(Number of frequencies)--"},

  };

  wsm_data["f_gridFromGasAbsLookup"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *f_grid* to the frequency grid of *abs_lookup*.

Must be called between importing/creating raw absorption table and
call of *abs_lookupAdapt*.
)--",
      .author = {"Stefan Buehler"},
      .out = {"f_grid"},

      .in = {"abs_lookup"},

  };

  wsm_data["f_gridFromSensorAMSU"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Automatically calculate f_grid to match the sensor.

This method is handy if you are simulating an AMSU-type instrument,
consisting of a few discrete channels. The case that channels touch,
as for MHS, is handled correctly. But the case that channels overlap
is not (yet) handled and results in an error message.

The method calculates *f_grid* to match the instrument, as given by
the local oscillator frequencies *lo_multi*, the backend
frequencies *f_backend_multi*, and the backend channel
responses *backend_channel_response_multi*.

You have to specify the desired spacing in the keyword ``spacing``,
which has a default value of 100 MHz. (The actual value is 0.1e9,
since our unit is Hz.)

The produced grid will not have exactly the requested spacing, but
will not be coarser than requested. The algorithm starts with the band
edges, then adds additional points until the spacing is at least as
fine as requested.

There is a similar method for HIRS-type instruments,
see *f_gridFromSensorHIRS*.
)--",
      .author = {"Stefan Buehler, Mathias Milz"},
      .out = {"f_grid"},

      .in = {"lo_multi", "f_backend_multi", "backend_channel_response_multi"},
      .gin = {"spacing"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{.1e9}},
      .gin_desc = {R"--(Desired grid spacing in Hz.)--"},

  };

  wsm_data["f_gridFromSensorAMSUgeneric"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Automatcially calculate f_grid to match the sensor. 
This function is based on 'f_gridFromSensorAMSU' 

The method calculates *f_grid* to match the instrument, as given by
the backend frequencies *f_backend*, and the backend channel
responses *backend_channel_response*.

You have to specify the desired spacing in the keyword ``spacing``,
which has a default value of 100 MHz. (The actual value is 0.1e9,
since our unit is Hz.)
The produced grid will not have exactly the requested spacing, but
it will not be coarser than requested. The algorithm starts with the band
edges, then adds additional points until the spacing is at least as
fine as requested.
)--",
      .author = {"Oscar Isoz"},
      .out = {"f_grid"},

      .in = {"f_backend_multi", "backend_channel_response_multi"},
      .gin = {"spacing", "verbosityVect"},
      .gin_type = {"Numeric", "Vector"},
      .gin_value = {Numeric{.1e9}, Vector{}},
      .gin_desc = {R"--(Desired grid spacing in Hz.)--",
                   R"--(Bandwidth adjusted spacing)--"},

  };

  wsm_data["f_gridFromSensorHIRS"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Automatically calculate f_grid to match the sensor.

This method is handy if you are simulating a HIRS-type instrument,
consisting of a few discrete channels.

It calculates f_grid to match the instrument, as given by the nominal
band frequencies *f_backend* and the spectral channel response
functions given by *backend_channel_response*.

You have to specify the desired spacing in the keyword ``spacing``, which
has a default value of 5e8 Hz.

The produced grid will not have exactly the requested spacing, but
will not be coarser than requested. The algorithm starts with the band
edges, then adds additional points until the spacing is at least as
fine as requested.

There is a similar method for AMSU-type instruments, see
*f_gridFromSensorAMSU*.
)--",
      .author = {"Stefan Buehler"},
      .out = {"f_grid"},

      .in = {"f_backend", "backend_channel_response"},
      .gin = {"spacing"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{5e8}},
      .gin_desc = {R"--(Desired grid spacing in Hz.)--"},

  };

  wsm_data["f_gridMetMM"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *f_grid* and associated variables match MetMM settings.

The method calculates *f_grid* to match the specifications of a
*met_mm_backend* table and method arguments.

You have to specify the desired spacing using the keyword ``freq_spacing``.
You can pass a *Vector* with one element to apply the same spacing to all
channels or pass a spacing value for each channel separately.

Optionally, ``freq_number`` can be set to specify the mininum number of
frequencies per passband for each channel. The frequencies are placed
equally spaced in each passband. The minimum spacing resulting from
``freq_number`` and ``freq_spacing`` will be used for the calculation. To
explicitly use ``freq_spacing`` for a channel, ``freq_number`` can be set
to -1 for this channel.

The number of elements in ``freq_number`` can either be the number of
channels or 1. If only one element is given, this number is used for
all channels. If ``freq_number`` is 1 and ``freq_spacing`` is wider than
the bandwidth of the channel, one frequency is placed in the middle of
each passband.

Frequencies that would be closer than ``freq_merge_threshold`` in the
generated *f_grid* are merged together. This value should be left at
the default value. This is only meant to compensate for numerical
inaccuracies in the frequency calculation to merge frequency that are
supposed to be identical.
)--",
      .author = {"Oliver Lemke", "Patrick Eriksson"},
      .out = {"f_grid",
              "f_backend",
              "channel2fgrid_indexes",
              "channel2fgrid_weights"},

      .in = {"met_mm_backend"},
      .gin = {"freq_spacing", "freq_number", "freq_merge_threshold"},
      .gin_type = {"Vector", "ArrayOfIndex", "Numeric"},
      .gin_value = {Vector{.1e9}, ArrayOfIndex{-1}, Numeric{1}},
      .gin_desc =
          {R"--(Desired grid spacing in Hz.)--",
           R"--(Number of frequencies per passband for each channel.)--",
           R"--(Merge frequencies that are closer than this value in Hz.)--"},

  };

  wsm_data["forloop_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *forloop_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--",
      .author = {"Richard Larsson"},
      .out = {"forloop_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["g0Earth"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gravity at zero altitude on Earth.

Sets *g0* for the given latitude using a standard parameterisation.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"g0"},

      .in = {"lat"},

  };

  wsm_data["g0Io"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gravity at zero altitude on Io.

Numeric from Wikipedia.
)--",
      .author = {"Richard Larsson"},
      .out = {"g0"},

  };

  wsm_data["g0Jupiter"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gravity at zero altitude on Jupiter.

Sets *g0*  to mean equatorial gravity on Jupiter. Value provided by
MPS under ESA-planetary study (TN1).
)--",
      .author = {"Jana Mendrok"},
      .out = {"g0"},

  };

  wsm_data["g0Mars"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gravity at zero altitude on Mars.

Sets *g0*  to mean equatorial gravity on Mars. Value provided by
MPS under ESA-planetary study (TN1).
)--",
      .author = {"Jana Mendrok"},
      .out = {"g0"},

  };

  wsm_data["g0Venus"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gravity at zero altitude on Venus.

Sets *g0*  to mean equatorial gravity on Venus. Value from Ahrens
(1995), provided by MPS under ESA-planetary study (TN1).
)--",
      .author = {"Jana Mendrok"},
      .out = {"g0"},

  };

  wsm_data["g0_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *g0_agenda* to a default value

Options are:

- ``"Earth"``: Uses *g0Earth* to set *g0*
- ``"Io"``: Uses *g0Io* to set *g0*
- ``"Jupiter"``: Uses *g0Jupiter* to set *g0*
- ``"Mars"``: Uses *g0Mars* to set *g0*
- ``"Venus"``: Uses *g0Venus* to set *g0*
)--",
      .author = {"Richard Larsson"},
      .out = {"g0_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["gas_scatteringOff"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Deactivates the gas_scattering within radiative transfer calculations.
)--",
      .author = {"Manfred Brath"},
      .out = {"gas_scattering_do", "gas_scattering_agenda"},

      .pass_workspace = true,

  };

  wsm_data["gas_scattering_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *gas_scattering_agenda* to a default value

Options are:

- ``"Dummy"``:

    1. Will *Ignore* all agenda inputs
    2. Uses *Touch* on all agenda outputs
)--",
      .author = {"Richard Larsson"},
      .out = {"gas_scattering_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {String("Dummy")},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["gas_scattering_coefAirSimple"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates of scattering coefficient matrix for air.

This function calculates the scattering coefficient for air using a
fitted version from Stamnes et al., 2017 of the numerical results of
Bates, 1984. Internally it calculates the spectrum of scattering
coefficient matrices from the spectrum of scattering cross section matrices,
atmospheric pressure, temperature for one point in the atmosphere. The
function multiplies the cross sections with the number density of gas
molecules under the assumption of an ideal gas to get the coefficients.
The result is returned in *gas_scattering_coef*. The atmospheric  pressure  and 
temperature  state  has  to  be  specified by  *rtp_pressure*,
*rtp_temperature*. The formula is accurate to 0.3 percent for wavelengths
between 0.205 and 1.05 micrometer.
)--",
      .author = {"Jon Petersen"},
      .out = {"gas_scattering_coef"},

      .in = {"f_grid", "rtp_pressure", "rtp_temperature"},

  };

  wsm_data["gas_scattering_coefXsecConst"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates the spectrum of scattering coefficient matrices.

It calculates the spectrum of scattering coefficient matrices from 
constant spectrum of scattering cross section matrices, atmospheric pressure,
temperature for one point in the atmosphere. Basically, it multiplies
the cross sections with the number density of gas molecules under the
assumption of an ideal gas. The result is returned in *gas_scattering_coef*. The
atmospheric  pressure  and  temperature  state  has  to  be  specified
by  *rtp_pressure*, *rtp_temperature*.
)--",
      .author = {"Manfred Brath"},
      .out = {"gas_scattering_coef"},

      .in = {"f_grid", "rtp_pressure", "rtp_temperature"},
      .gin = {"ConstXsec"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{0.}},
      .gin_desc = {R"--(Constant Xsec value)--"},

  };

  wsm_data["gas_scattering_matIsotropic"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates the spectrum of normalized scattering matrices.
Important, the angular direction are line of sight direction not the
propagation direction.
)--",
      .author = {"Manfred Brath"},
      .out = {"gas_scattering_mat", "gas_scattering_fct_legendre"},

      .in = {"gas_scattering_los_in",
             "gas_scattering_los_out",
             "gas_scattering_output_type"},

  };

  wsm_data["gas_scattering_matRayleigh"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates the normalized Rayleigh scattering matrix.

The phase matrix for anisotropic Rayleigh particles in random orientations.Important, the angular direction are defined as line of sight direction not as
propagation direction.
)--",
      .author = {"Jon Petersen"},
      .out = {"gas_scattering_mat", "gas_scattering_fct_legendre"},

      .in = {"gas_scattering_los_in",
             "gas_scattering_los_out",
             "gas_scattering_output_type"},
      .gin = {"gdepolarization_factor"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{0.03}},
      .gin_desc = {R"--(depolarization factor for air)--"},

  };

  wsm_data["geo_posEndOfPpath"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets geo-position based on *ppath*.

The geo-position is set to the position of the last point in *ppath*.

NaN is returned if *ppath* is totally outside of the atmosphere.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"geo_pos"},

      .in = {"ppath"},

  };

  wsm_data["geo_posLowestAltitudeOfPpath"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets geo-position based on *ppath*.

The geo-position is set to the position of the point in *ppath*
having the lowest altitude.

NaN is returned if *ppath* is totally outside of the atmosphere.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"geo_pos"},

      .in = {"ppath"},

  };

  wsm_data["geo_posWhereAltitudeIsPassed"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets geo-position based on *ppath*.

The geo-position is set to the position where the propagation
path passes the stated altitude. If this altitude is passed
more than once, the passing closest to the sensor is selected.
If the reference altitude is not passed at all, *geo_pos* is
set to NaN.

NaN is also returned if *ppath* is totally outside of the atmosphere.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"geo_pos"},

      .in = {"ppath"},
      .gin = {"altitude"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Altitude defining *geo_pos*.)--"},

  };

  wsm_data["heating_ratesFromIrradiance"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates heating rates from the *irradiance_field*.

The method assumes that the heating rates depend only on the
vertical derivation of the net flux. The net flux is the sum of the
*irradiance_field* in upward direction and the *irradiance_field*
in downward direction
)--",
      .author = {"Manfred Brath"},
      .out = {"heating_rates"},

      .in = {"ppvar_atm", "irradiance_field", "specific_heat_capacity", "g0"},

  };

  wsm_data["inversion_iterate_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *inversion_iterate_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--",
      .author = {"Richard Larsson"},
      .out = {"inversion_iterate_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["irradiance_fieldFromRadiance"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates the irradiance from the *radiance_field*.

The *radiance_field* is integrated over the angular grids according to
the grids set by *AngularGridsSetFluxCalc*. 
See *AngularGridsSetFluxCalc* to set *za_grid*, *aa_grid*, and
*za_grid_weights*
)--",
      .author = {"Manfred Brath"},
      .out = {"irradiance_field"},

      .in = {"radiance_field", "za_grid", "aa_grid", "za_grid_weights"},

  };

  wsm_data["isotopologue_ratiosInitFromBuiltin"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Initialize isotopologue ratios with default values from built-in
species data.  This should be OK for Earth-like atmospheres
)--",
      .author = {"Oliver Lemke"},
      .out = {"isotopologue_ratios"},

  };

  wsm_data["isotopologue_ratiosInitFromHitran"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Initialize isotopologue ratios with default values from built-in
Hitran species data.
)--",
      .author = {"Richard Larsson"},
      .out = {"isotopologue_ratios"},

  };

  wsm_data["iyApplyUnit"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Conversion of *iy* to other spectral units (for passive observations).

The method allows a change of unit, as a post-processing step,
ignoring the n2-law of radiance.

The conversion made inside ``iyEmissionStandard`` is mimiced,
see that method for constraints and selection of output units.
Restricted to that the n2-law can be ignored. This assumption
is valid if the sensor is placed in space, or if the refractive
index only deviates slightly from unity.

It is stressed that there is no automatic check that the method is
applied correctly, it is up to the user to ensure that the input
data are suitable for the conversion.

Beside *iy*, these auxilary quantities are modified:
    "iy", "Error" and "Error (uncorrelated)"

Please note that *diy_dx* is not handled. Also note that this method
considers *iy_unit*, while *iy_unit_radar* is handled directly by
the methods dealing with such simulations.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"iy", "iy_aux"},

      .in = {"iy", "iy_aux", "f_grid", "iy_aux_vars", "iy_unit"},

  };

  wsm_data["iyBackground"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Computes background radiation by wrapping various agendas
)--",
      .author = {"Richard Larsson"},
      .out = {"iy", "diy_dx"},

      .in = {"iy_transmittance",
             "background_transmittance",
             "surface_field",
             "f_grid",
             "rte_pos2",
             "ppath",
             "atm_field",
             "jacobian_quantities",
             "jacobian_do",
             "cloudbox_on",
             "iy_id",
             "iy_agenda_call1",
             "iy_main_agenda",
             "iy_space_agenda",
             "iy_surface_agenda",
             "iy_cloudbox_agenda",
             "iy_unit"},

      .pass_workspace = true,

  };

  wsm_data["iyCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(A single monochromatic pencil beam calculation.

Performs monochromatic radiative transfer calculations for the
specified position (*rte_pos*) and line-of-sight (*rte_pos*).
See *iy* and associated variables for format of output.

Please note that Jacobian type calculations not are supported.
For this use *yCalc*.

No sensor characteristics are applied. These are most easily
incorporated by using *yCalc*
)--",
      .author = {"Patrick Eriksson"},
      .out = {"iy", "iy_aux", "ppath", "geo_pos"},

      .in = {"atmgeom_checked",
             "atmfields_checked",
             "iy_aux_vars",
             "iy_id",
             "cloudbox_on",
             "cloudbox_checked",
             "scat_data_checked",
             "f_grid",
             "atm_field",
             "rte_pos",
             "rte_los",
             "rte_pos2",
             "iy_unit",
             "iy_main_agenda"},

      .pass_workspace = true,

  };

  wsm_data["iyCopyPath"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Copies the radiative transfer properties to their matpack equivalents.
)--",
      .author = {"Richard Larsson"},
      .out = {"iy",
              "ppvar_iy",
              "ppvar_trans_cumulat",
              "ppvar_trans_partial",
              "diy_dpath"},

      .in = {"ppvar_rad",
             "ppvar_drad",
             "ppvar_cumtramat",
             "ppvar_tramat",
             "jacobian_quantities",
             "jacobian_do"},

  };

  wsm_data["iyEmissionHybrid"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Radiative transfer with emission and precalculated radiation field.

This method works largely as ``iyEmissionStandard`` but incorporates
scattering by a precalculated radiation field. It is so far limited
to 1D calculations.

The method integrates the source term along the propagation path. While
``iyEmissionStandard`` only considers local thermal emission, this method
also includes scattering into the line-of-sight in the source term.
The scattering integral is solved with the precalculated field as incoming
radiation. That is, this method extends the integration into the cloudbox,
while ``iyEmissionStandard`` starts at the cloudbox boundary.

The calculate radiance should be as exact as what is produced by the
scattering solver used to calculate the precalculted radiation field,
but the main reason to use this method is to obtain the Jacobian even
in the presence of scattering. The Jacobian with respect to bulk scattering
properties can be obtained, but it is approximate. This is the case as
the incoming radiation field is treated as fixed in the calculation
of the Jacobian. The impact of this approximation increases with the degree
of scattering.
)--",
      .author = {"Patrick Eriksson", "Jana Mendrok", "Richard Larsson"},
      .out = {"iy",
              "iy_aux",
              "diy_dx",
              "ppvar_atm",
              "ppvar_pnd",
              "ppvar_f",
              "ppvar_iy",
              "ppvar_trans_cumulat",
              "ppvar_trans_partial"},

      .in = {"diy_dx",
             "iy_id",
             "f_grid",
             "abs_species",
             "atm_field",
             "cloudbox_on",
             "cloudbox_limits",
             "pnd_field",
             "dpnd_field_dx",
             "scat_species",
             "scat_data",
             "iy_unit",
             "iy_aux_vars",
             "jacobian_do",
             "jacobian_quantities",
             "propmat_clearsky_agenda",
             "water_p_eq_agenda",
             "rt_integration_option",
             "iy_main_agenda",
             "iy_space_agenda",
             "iy_surface_agenda",
             "iy_cloudbox_agenda",
             "iy_agenda_call1",
             "iy_transmittance",
             "ppath",
             "rte_pos2",
             "rte_alonglos_v",
             "surface_field",
             "cloudbox_field",
             "za_grid"},
      .gin = {"Naa_grid", "t_interp_order"},
      .gin_type = {"Index", "Index"},
      .gin_value = {Index{19}, Index{1}},
      .gin_desc =
          {R"--(Number of azimuth angles to consider in scattering source term integral.)--",
           R"--(Interpolation order of temperature for scattering data (so far only applied in phase matrix, not in extinction and absorption.)--"},
      .pass_workspace = true,

  };

  wsm_data["iyIndependentBeamApproximation"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Samples atmosphere along ppath and make 1D-type RT calculation.

The main application of this method should be to apply 1D
scattering solvers on 2D or 3D atmospheres. The 1D calculation
is set up inside *iy_independent_beam_approx_agenda*.

The method calculates the ppath until reaching the surface or the
top-of-the atmosphere. If the sensor is inside the atmosphere the
ppath is extended from the sensor vertically to cover the remaining
part of the atmosphere. All atmospheric fields are interpolated to
the obtain ppath, to form a 1D view of the atmosphere. This 1D
atmosphere forms the input to *iy_independent_beam_approx_agenda*

*lat_true* and *lon_true* for the 1D view is set to the intersection
with the surface of the ppath described above.

The function accepts that the input atmosphere is 1D, as well as
that there is no active cloudbox.

The constructed 1D atmosphere is exported if the GIN ``return_atm1d``
is set to 1. The default then is to include all atmospheric fields,
but ``vmr_field`` and *pnd_field* can be deselected by two of the GIN-s.
The order of the fields is specified by the first grid in the structure
member grids. If *atm_fields_compact* is denoted as A, then
A.grids{0}{i} gives the name of field with index i.
Each book in ``vmr_field`` and *pnd_field* is stored separately. For
example, the first book in *pnd_field* is stored with the name
"Scattering element 0".
)--",
      .author = {"Patrick Eriksson"},
      .out = {"iy", "iy_aux", "ppath", "diy_dx", "atm_fields_compact"},

      .in = {"diy_dx",          "iy_id",
             "f_grid",          "atm_field",
             "cloudbox_on",     "cloudbox_limits",
             "pnd_field",       "particle_masses",
             "ppath_agenda",    "ppath_lmax",
             "ppath_lraytrace", "iy_agenda_call1",
             "iy_unit",         "iy_transmittance",
             "rte_pos",         "rte_los",
             "rte_pos2",        "jacobian_do",
             "iy_aux_vars",     "iy_independent_beam_approx_agenda"},
      .gin = {"return_atm1d", "skip_vmr", "skip_pnd", "return_masses"},
      .gin_type = {"Index", "Index", "Index", "Index"},
      .gin_value = {Index{0}, Index{0}, Index{0}, Index{0}},
      .gin_desc =
          {R"--(Flag to trigger that *atm_fields_compact* is filled. )--",
           R"--(Flag to not include vmr data in *atm_fields_compact*.)--",
           R"--(Flag to not include pnd data in *atm_fields_compact*.)--",
           R"--(Flag to include particle category masses in *atm_fields_compact*.Conversion is done by *particle_masses*.)--"},
      .pass_workspace = true,

  };

  wsm_data["iyLoopFrequencies"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Radiative transfer calculations one frequency at the time.

The method loops the frequencies in *f_grid* and calls
*iy_loop_freqs_agenda* for each individual value. This method is
placed in *iy_main_agenda*, and the actual radiative transfer
method is put in *iy_loop_freqs_agenda*.

A common justification for using the method should be to consider
dispersion. By using this method it is ensured that the propagation
path for each individual frequency is calculated.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"iy", "iy_aux", "ppath", "diy_dx"},

      .in = {"iy_aux_vars",
             "iy_agenda_call1",
             "iy_transmittance",
             "rte_pos",
             "rte_los",
             "rte_pos2",
             "f_grid",
             "iy_loop_freqs_agenda"},

      .pass_workspace = true,

  };

  wsm_data["iyMC"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Interface to Monte Carlo part for *iy_main_agenda*.

Basically an interface to *MCGeneral* for doing monochromatic
pencil beam calculations. This functions allows Monte Carlo (MC)
calculations for sets of frequencies and sensor pos/los in a single
run. Sensor responses can be included in the standard manner
(through *yCalc*).

This function does not apply the MC approach when it comes
to sensor properties. These properties are not considered when
tracking photons, which is done in *MCGeneral* (but then only for
the antenna pattern).

Output unit options  (*iy_unit*) exactly as for *MCGeneral*.

The MC calculation errors are all assumed be uncorrelated and each
have a normal distribution. These properties are of relevance when
weighting the errors with the sensor repsonse matrix. The seed is
reset for each call of *MCGeneral* to obtain uncorrelated errors.

MC control arguments (mc_std_err, mc_max_time, mc_min_iter, mc_max_iter
mc_taustep_limit) as for *MCGeneral*. The arguments are applied
for each monochromatic pencil beam calculation individually.
As for *MCGeneral*, the value of *mc_error* shall be adopted to
*iy_unit*.

The following auxiliary data can be obtained:

- ``"Error (uncorrelated)"``:
    Calculation error. Size: [nf,ns,1,1].
    (The later part of the text string is required. It is used as
    a flag to yCalc for how to apply the sensor data.)

where

- `nf`: Number of frequencies.
- `ns`: Number of Stokes elements.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"iy", "iy_aux", "diy_dx"},

      .in = {"iy_agenda_call1",
             "iy_transmittance",
             "rte_pos",
             "rte_los",
             "iy_aux_vars",
             "jacobian_do",
             "atm_field",
             "surface_field",
             "cloudbox_on",
             "cloudbox_limits",
             "f_grid",
             "scat_data",
             "iy_space_agenda",
             "surface_rtprop_agenda",
             "propmat_clearsky_agenda",
             "ppath_step_agenda",
             "ppath_lmax",
             "ppath_lraytrace",
             "pnd_field",
             "iy_unit",
             "mc_std_err",
             "mc_max_time",
             "mc_max_iter",
             "mc_min_iter",
             "mc_taustep_limit"},
      .gin = {"t_interp_order"},
      .gin_type = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc =
          {R"--(Interpolation order of temperature for scattering data (so far only applied in phase matrix, not in extinction and absorption.)--"},
      .pass_workspace = true,

  };

  wsm_data["iyRadarSingleScat"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Simulation of radar, restricted to single scattering.

The WSM treats e.g. radar measurements of cloud and precipitation,
on the condition that multiple scattering can be ignored. Beside
the direct back-scattering, the two-way attenuation by gases and
particles is considered. Surface scattering/clutter is ignored.

The method could potentially be used for lidars, but multiple
scattering poses here a must stronger constrain for the range of
applications.

The method shall be used with *yRadar*, NOT with *yCalc*.

The *ppath* provided should be calculated including cloudbox interior::

     ppathCalc( cloudbox_on=0 )

The method returns the back-scattering for each point of *ppath*.
Several frequencies can be treated in parallel. The size of *iy*
is [ nf * np, stokes_dim ], where nf is the length of *f_grid* and
np is the number of path points. The data are stored in blocks
of [ np, stokes_dim ]. That is, all the results for the first
frequency occupy the np first rows of *iy* etc.

The polarisation state of the transmitted pulse is taken from
*iy_transmitter*. If the radar transmits several polarisations at
the same frequency, you need to handle this by using two frequencies
in *f_grid*, but these can be almost identical.

This method does not consider *iy_unit_radar*. Unit changes are instead
applied in *yRadar*. The output of this method matches the option "1".

The extinction due to particles can be scaled (by ``pext_scaling``),
which could be of interest when e.g. characterising inversions or
trying to compensate for ignored multiple scattering. The later is
commented further for *particle_bulkpropRadarOnionPeeling*.

For Jacobian calculations the default is to assume that the
transmittance is unaffected by the retrieval quantities. This is
done to save computational time, and should be a valid approximation
for the single-scattering conditions. Set ``trans_in_jacobian`` to 1 to
activate full Jacobian calculations.

Some auxiliary radiative transfer quantities can be obtained. Auxiliary
quantities are selected by *iy_aux_vars* and returned by *iy_aux*.
Valid choices for auxiliary data are:

- ``"Radiative background"``:
    Index value flagging the radiative
    background. The following coding is used: 0=space, 1=surface
    and 2=cloudbox (the last case should not occur!). Only column
    matching first Stokes element filled. Other columns are set to 0.
- ``"Backscattering"``:
    The unattenuated back-scattering. That is, as
    *iy* but with no attenuated applied. Here all columns are filled.
    By combing *iy* and this auxiliary variable, the total two-way
    attenuation can be derived.
- ``"Abs species extinction"``:
    Extinction due to *abs_species* at each
    ppath point, taken as the diagonal of the local extinction matrix.
- ``"Particle extinction"``:
    Extinction due to particles at each
    ppath point, taken as the diagonal of the local extinction matrix.
    The retunred values includes ``pext_scaling``
)--",
      .author = {"Patrick Eriksson"},
      .out = {"iy", "iy_aux", "diy_dx", "ppvar_atm", "ppvar_pnd", "ppvar_f"},

      .in = {"f_grid",
             "abs_species",
             "atm_field",
             "cloudbox_on",
             "cloudbox_limits",
             "pnd_field",
             "dpnd_field_dx",
             "scat_species",
             "scat_data",
             "scat_data_checked",
             "iy_aux_vars",
             "jacobian_do",
             "jacobian_quantities",
             "ppath",
             "iy_transmitter",
             "propmat_clearsky_agenda",
             "water_p_eq_agenda",
             "rte_alonglos_v"},
      .gin = {"trans_in_jacobian", "pext_scaling", "t_interp_order"},
      .gin_type = {"Index", "Numeric", "Index"},
      .gin_value = {Index{0}, Numeric{1}, Index{1}},
      .gin_desc =
          {R"--(Flag determining if change in transmittance is considered in calculation of the Jacobian or not.)--",
           R"--(Particle extinction is scaled with this value. A value inside [0,2]. Set it to 0 if you want to remove particle extinction totally.)--",
           R"--(Interpolation order of temperature for scattering data (so far only applied in phase matrix, not in extinction and absorption.)--"},
      .pass_workspace = true,

  };

  wsm_data["iyReplaceFromAux"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Change of main output variable.

With this method you can replace the content of *iy* with one of
the auxiliary variables. The selected variable (by ``aux_var``) must
be part of *iy_aux_vars*. The corresponding data from *iy_aux* are
copied to form a new *iy* (*iy_aux* is left unchanged). Elements of
*iy* correponding to Stokes elements not covered by the auxiliary
variable are just set to zero.

Jacobian variables are not handled.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"iy"},

      .in = {"iy", "iy_aux", "iy_aux_vars", "jacobian_do"},
      .gin = {"aux_var"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Auxiliary variable to insert as *iy*.)--"},

  };

  wsm_data["iySurfaceFastem"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Usage of FASTEM for emissivity and reflectivity of water surfaces.

This method allows usage of the FASTEM model inside
*iy_surface_agenda*. The aim is to use FASTEM in the exact same
way as done in RTTOV. For example, the transmittance for down-
welling radiation is considered. RTTOV os just 1D. Here 2D and 3D
are handled as the 1D case, the down-welling radiation is just
calculated for the directuon matching specular reflection.

The wind direction is given as the azimuth angle, counted
clockwise from north (i.e. an easterly wind is at 90 deg).
This matches the general definition of azimuth inside ARTS.
For 1D and 2D, the wind direction must be adjusted to match the
fact that the line-of-sight is locked to be at 0 deg (180 for 2D
in the case of a negative zenith angle). For 3D, the true wind
direction shall be used.

FASTEM is called by *FastemStandAlone*. See that WSM for further
comments on variables and limitations.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"iy", "diy_dx"},

      .in = {"diy_dx",
             "iy_transmittance",
             "iy_id",
             "jacobian_do",
             "atm_field",
             "cloudbox_on",
             "f_grid",
             "rtp_pos",
             "rtp_los",
             "rte_pos2",
             "iy_unit",
             "iy_main_agenda",
             "surface_skin_t"},
      .gin = {"salinity", "wind_speed", "wind_direction", "fastem_version"},
      .gin_type = {"Numeric", "Numeric", "Numeric", "Index"},
      .gin_value = {Numeric{0.035}, std::nullopt, Numeric{0}, Index{6}},
      .gin_desc = {R"--(Salinity, 0-1. That is, 3% is given as 0.03.)--",
                   R"--(Wind speed.)--",
                   R"--(Wind direction. See further above.)--",
                   R"--(The version of FASTEM to use.)--"},
      .pass_workspace = true,

  };

  wsm_data["iySurfaceFlatRefractiveIndex"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(This method calculates upwelling radiation for a specular flat surface.

These are due to the reflection of the downgoing diffuse radiation and emission from
the surface using a predefined reflectivity matrix. 

This method is designed to be part of *iy_surface_agenda*

Important this method calculates only the reflection of the diffuse
downward radiation. No direct incoming radiation is considered

Jacobian is supported only for Skin temperature
)--",
      .author = {"Manfred Brath"},
      .out = {"iy", "diy_dx", "dsurface_rmatrix_dx", "dsurface_emission_dx"},

      .in = {"iy",
             "diy_dx",
             "dsurface_rmatrix_dx",
             "dsurface_emission_dx",
             "iy_transmittance",
             "iy_id",
             "jacobian_do",
             "suns_do",
             "atm_field",
             "cloudbox_on",
             "f_grid",
             "surface_field",
             "rtp_pos",
             "rtp_los",
             "rte_pos2",
             "iy_unit",
             "surface_complex_refr_index",
             "surface_props_names",
             "dsurface_names",
             "jacobian_quantities",
             "iy_main_agenda"},

      .pass_workspace = true,

  };

  wsm_data["iySurfaceInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(This method initialize iy.

This method is designed to be part of *iy_surface_agenda*.
Its only prpose is to initialize *iy* properly within the 
*iy_surface_agenda*
)--",
      .author = {"Manfred Brath"},
      .out = {"iy"},

      .in = {"f_grid"},

  };

  wsm_data["iySurfaceLambertian"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(This method calculates upwelling radiation for a lambertian surface.

These are due to the scattering of the downgoing diffuse radiation and emission from
the surface.
This method works only for 1D or 3D atmospheres.
For the integration over the zenith angles a gaussian quadrature with
N_za angles is used.
For 1D atmospheres N_aa is ignored. For 3D atmospheres without clouds
azimuthal dependency can be neglected. N_aa = 1 is sufficient.
For 3D atmospheres with cloudbox on azimuthal dependency needs to be 
accounted. In that case the number of azimuth angles N_aa as a rule ofthumb should be set to 4*N_za.
For the 1D case N_za downwelling streams and 3D case N_za*N_aa downwelling
streams are calculated.

This method is designed to be part of *iy_surface_agenda*

Important this method calculates only the scattering of the diffuse
downward radiation. No direct incoming radiation is considered

Jacobian is supported only for Skin temperature
)--",
      .author = {"Manfred Brath"},
      .out = {"iy", "diy_dx"},

      .in = {"iy",
             "diy_dx",
             "iy_transmittance",
             "iy_id",
             "jacobian_do",
             "suns_do",
             "atm_field",
             "cloudbox_on",
             "f_grid",
             "surface_field",
             "rtp_pos",
             "rtp_los",
             "rte_pos2",
             "iy_unit",
             "surface_scalar_reflectivity",
             "surface_props_names",
             "dsurface_names",
             "jacobian_quantities",
             "iy_main_agenda"},
      .gin = {"N_za", "N_aa"},
      .gin_type = {"Index", "Index"},
      .gin_value = {Index{3}, Index{1}},
      .gin_desc = {R"--(Number of zenith angles.)--",
                   R"--(Number of azimuth angles)--"},
      .pass_workspace = true,

  };

  wsm_data["iySurfaceRtpropAgenda"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Interface to *surface_rtprop_agenda* for *iy_surface_agenda*.

This method is designed to be part of *iy_surface_agenda*. It
determines the radiative properties of the surface by
*surface_rtprop_agenda* and calculates the downwelling radiation
by *iy_main_agenda*, and sums up the terms as described in AUG.
That is, this WSM uses the output from *surface_rtprop_agenda*
in a straightforward fashion.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"iy",
              "diy_dx",
              "surface_point",
              "surface_los",
              "surface_rmatrix",
              "surface_emission"},

      .in = {"diy_dx",
             "iy_transmittance",
             "iy_id",
             "jacobian_do",
             "suns_do",
             "atm_field",
             "cloudbox_on",
             "f_grid",
             "rtp_pos",
             "rtp_los",
             "rte_pos2",
             "iy_unit",
             "iy_main_agenda",
             "surface_rtprop_agenda"},

      .pass_workspace = true,

  };

  wsm_data["iySurfaceRtpropCalc"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Applies *surface_los*, *surface_rmatrix* and *surface_emission*.

This method is designed to be part of *iy_surface_agenda* and
should be mandatory when using methods describing the surface
radiative transfer properties by *surface_los*, *surface_rmatrix*
and *surface_emission*. The task of this method is to apply these
three WSVs to obtain the upwelling radiation from the surface.
This upwelling radiation is the sum of surface emission and
reflected downwelling radiation. The later part is calculated
by calling *iy_main_agenda*. See further AUG.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"iy", "diy_dx"},

      .in = {"diy_dx",
             "surface_los",
             "surface_rmatrix",
             "surface_emission",
             "dsurface_names",
             "dsurface_rmatrix_dx",
             "dsurface_emission_dx",
             "iy_transmittance",
             "iy_id",
             "jacobian_do",
             "suns_do",
             "jacobian_quantities",
             "atm_field",
             "cloudbox_on",
             "f_grid",
             "rtp_pos",
             "rtp_los",
             "rte_pos2",
             "iy_unit",
             "iy_main_agenda"},

      .pass_workspace = true,

  };

  wsm_data["iyUnitConversion"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Perform unit conversions of *iy*, *diy_dx*, and *ppvar_iy*.
)--",
      .author = {"Richard Larsson"},
      .out = {"iy", "diy_dx", "ppvar_iy"},

      .in = {"iy",
             "diy_dx",
             "ppvar_iy",
             "f_grid",
             "ppath",
             "jacobian_quantities",
             "iy_unit",
             "jacobian_do",
             "iy_agenda_call1"},

  };

  wsm_data["iy_auxFromVars"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Set *iy_aux* from list of parameters.
)--",
      .author = {"Richard Larsson"},
      .out = {"iy_aux"},

      .in = {"iy_aux_vars",
             "background_transmittance",
             "ppath",
             "iy_agenda_call1"},

  };

  wsm_data["iy_cloudbox_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *iy_cloudbox_agenda* to a default value

Options are:

- ``"LinInterpField"``: Uses ``iyInterpCloudboxField`` to set *iy*
- ``"QuarticInterpField"``: Uses ``iyInterpCloudboxField`` to set *iy* using ``za_interp_order=4``
)--",
      .author = {"Richard Larsson"},
      .out = {"iy_cloudbox_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["iy_independent_beam_approx_agendaSet"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Sets *iy_independent_beam_approx_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--",
          .author = {"Richard Larsson"},
          .out = {"iy_independent_beam_approx_agenda"},

          .gin = {"option"},
          .gin_type = {"String"},
          .gin_value = {std::nullopt},
          .gin_desc = {R"--(Default agenda option (see description))--"},
      };

  wsm_data["iy_loop_freqs_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *iy_loop_freqs_agenda* to a default value

Options are:

- ``"Emission"``:

    1. Uses ``ppathCalc`` to set *ppath*
    2. Uses ``iyEmissionStandard`` to set *iy*, *iy_aux*, ``ppvar_p``, ``ppvar_t``, ``ppvar_nlte``, ``ppvar_vmr``, ``ppvar_wind``, ``ppvar_mag``, *ppvar_f*, *ppvar_iy*, *ppvar_trans_cumulat*, and *ppvar_trans_partial*, and also  to modify *diy_dx*

- ``"Transmission"``:

    1. Uses ``ppathCalc`` to set *ppath*
    2. Uses ``iyTransmissionStandard`` to set *iy*, *iy_aux*, ``ppvar_p``, ``ppvar_t``, ``ppvar_nlte``, ``ppvar_vmr``, ``ppvar_wind``, ``ppvar_mag``, *ppvar_f*, *ppvar_iy*, *ppvar_trans_cumulat*, and *ppvar_trans_partial*, and also to modify *diy_dx*
)--",
      .author = {"Richard Larsson"},
      .out = {"iy_loop_freqs_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["iy_main_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *iy_main_agenda* to a default value

Options are:

- ``"Emission"``:

    1. Uses ``ppathCalc`` to set *ppath*
    2. Uses ``iyEmissionStandard`` to set *iy*, *iy_aux*, ``ppvar_p``, ``ppvar_t``, ``ppvar_nlte``, ``ppvar_vmr``, ``ppvar_wind``, ``ppvar_mag``, *ppvar_f*, *ppvar_iy*, *ppvar_trans_cumulat*, and *ppvar_trans_partial*, and also  to modify *diy_dx*
    3. Sets *geo_pos* to empty

- ``"EmissionPlaneParallel"``:

    1. Uses ``ppathPlaneParallel`` to set *ppath*
    2. Uses ``iyEmissionStandard`` to set *iy*, *iy_aux*, ``ppvar_p``, ``ppvar_t``, ``ppvar_nlte``, ``ppvar_vmr``, ``ppvar_wind``, ``ppvar_mag``, *ppvar_f*, *ppvar_iy*, *ppvar_trans_cumulat*, and *ppvar_trans_partial*, and also  to modify *diy_dx*
    3. Sets *geo_pos* to empty

- ``"Clearsky"``:

    1. Uses ``ppathCalc`` to set *ppath*
    2. Uses ``iyClearsky`` to set *iy*, *iy_aux*, ``ppvar_p``, ``ppvar_t``, ``ppvar_nlte``, ``ppvar_vmr``, ``ppvar_wind``, ``ppvar_mag``, *ppvar_f*, *ppvar_iy*, *ppvar_trans_cumulat*, and *ppvar_trans_partial*, and also  to modify *diy_dx*
    3. Sets *geo_pos* to empty

- ``"Transmission"``:

    1. Uses ``ppathCalc`` to set *ppath* using *cloudbox_on* = 0
    2. Uses ``iyTransmissionStandard`` to set *iy*, *iy_aux*, ``ppvar_p``, ``ppvar_t``, ``ppvar_nlte``, ``ppvar_vmr``, ``ppvar_wind``, ``ppvar_mag``, *ppvar_f*, *ppvar_iy*, *ppvar_trans_cumulat*, and *ppvar_trans_partial*, and also to modify *diy_dx*
    3. Sets *geo_pos* to empty

- ``"TransmissionUnitUnpolIntensity"``:

    1. Uses *MatrixUnitIntensity* using out = *iy_transmitter*, and f =* f_grid*
    2. Uses ``ppathCalc`` to set *ppath* using *cloudbox_on* = 0
    3. Uses ``iyTransmissionStandard`` to set *iy*, *iy_aux*, ``ppvar_p``, ``ppvar_t``, ``ppvar_nlte``, ``ppvar_vmr``, ``ppvar_wind``, ``ppvar_mag``, *ppvar_f*, *ppvar_iy*, *ppvar_trans_cumulat*, and *ppvar_trans_partial*, and also to modify *diy_dx*
    4. Sets *geo_pos* to empty

- ``"TransmissionUnitPolIntensity"``:

    1. Uses *iy_transmitterSinglePol* to set *iy_transmitter*
    2. Uses ``ppathCalc`` to set *ppath* using *cloudbox_on* = 0
    3. Uses ``iyTransmissionStandard`` to set *iy*, *iy_aux*, ``ppvar_p``, ``ppvar_t``, ``ppvar_nlte``, ``ppvar_vmr``, ``ppvar_wind``, ``ppvar_mag``, *ppvar_f*, *ppvar_iy*, *ppvar_trans_cumulat*, and *ppvar_trans_partial*, and also to modify *diy_dx*
    4. Sets *geo_pos* to empty

- ``"Freqloop"``:

    1. Uses *iyLoopFrequencies* to set *iy*, *iy_aux*, *ppath*, and *diy_dx*
    2. Sets *geo_pos* to empty
    3. Will *Ignore* the *diy_dx* agenda input

- ``"ScattMC"``:

    1. Uses *iyMC* to set *iy*, *iy_aux*, and *diy_dx*
    2. Sets *geo_pos* to empty
    3. Will *Ignore* the *diy_dx* agenda input
)--",
      .author = {"Richard Larsson"},
      .out = {"iy_main_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["iy_radar_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *iy_radar_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--",
      .author = {"Richard Larsson"},
      .out = {"iy_radar_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["iy_space_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *iy_space_agenda* to a default value

Options are:

- ``"CosmicBackground"``:

    1. Uses *MatrixCBR* using out = *iy*, and f = *f_grid*
)--",
      .author = {"Richard Larsson"},
      .out = {"iy_space_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {String("CosmicBackground")},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["iy_surface_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *iy_space_agenda* to a default value

Options are:

- ``"UseSurfaceRtprop"``:

    1. Uses *SurfaceDummy* to modify *dsurface_rmatrix_dx*, and *dsurface_emission_dx*
    2. Uses *iySurfaceRtpropAgenda* to set *iy*, *surface_skin_t*, *surface_los*, *surface_rmatrix*, and *surface_emission*, and also to modify *diy_dx*
)--",
      .author = {"Richard Larsson"},
      .out = {"iy_surface_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {String("UseSurfaceRtprop")},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["iy_transmitterMultiplePol"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Transmitted signal having multiple polarisations.

The method is intended to be used as possible input of ``iyTransmissionStandard``.
It sets *iy_transmitter* to describe the transmitted signal/pulses.
The polarisation state is taken from *instrument_pol*, where
*instrument_pol* must contain an element for each frequency in *f_grid*.
The transmitted signal/pulses are set to be of unit magnitude, such
as [1,1,0,0].
)--",
      .author = {"Patrick Eriksson"},
      .out = {"iy_transmitter"},

      .in = {"f_grid", "instrument_pol"},

  };

  wsm_data["iy_transmitterSinglePol"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Transmitted signal having a single polarisations.

The method is intended to be used as possible input of ``iyTransmissionStandard``.
It sets *iy_transmitter* to describe the transmitted signal/pulses.
The polarisation state is taken from *instrument_pol*, where
*instrument_pol* must contain a single value.
This polarisation state is applied for all frequencies.
The transmitted pulses/signals are set to be of unit
magnitude, such as [1,1,0,0].
)--",
      .author = {"Patrick Eriksson"},
      .out = {"iy_transmitter"},

      .in = {"f_grid", "instrument_pol"},

  };

  wsm_data["jacobianAddAbsSpecies"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes an absorption species in the Jacobian.

For 1D or 2D calculations the latitude and/or longitude grid of
the retrieval field should set to have zero length.

These retrieval units are at hand for all gas species:

- ``"vmr"``: Volume mixing ratio.
- ``"nd"``: Number density.
- ``"rel"``: Relative unit (e.g. 1.1 means 10% more of the gas).

For water vapour, also these units are at hand:

- ``"rh"``: Relative humidity.
- ``"q"``: Specific humidity.

Note that ``for_species_tag`` is used to indicate if species tag VMR,
rather than atmospheric gas VMR is calculated. Set it to 0 and we
calculate the atmospheric gas VMR, but this only works for "analytical".

Note that the Jacobian is set to zero where volume mixing ratio equals zero.

The number of elements added to the state vector (*x*) is::

  n_g1 * n_g2 * n_g3

where n_g1, n_g2 and n_g3 are the length of GIN ``g1``, ``g2`` and ``g3``,
respectively. Here empty vectors should be considered to have a length 1.
The elements are sorted with pressure as innermost loop, followed by
latitude and longitude as outermost loop.
)--",
      .author = {"Mattias Ekstrom", "Patrick Eriksson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda"},
      .gin = {"g1", "g2", "g3", "species", "unit", "for_species_tag"},
      .gin_type = {"Vector", "Vector", "Vector", "String", "String", "Index"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    String("vmr"),
                    Index{1}},
      .gin_desc = {R"--(Pressure retrieval grid.)--",
                   R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retreival grid.)--",
                   R"--(The species tag of the retrieval quantity.)--",
                   R"--(Retrieval unit. See above.)--",
                   R"--(Index-bool for acting on species tags or species.)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddBasicCatalogParameter"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Includes a basic catalog parameter in the Jacobian. These are constant
over all layers and so only a single vector output is returned.

The only basic catalog parameters currently supported are:

* ``"LineStrength"``
* ``"LineCenter"``

The ``catalog_identity`` should be able to identify one or many
lines in the catalog used for calculating the spectral absorption.
Note that partial matching for energy levels are allowed but not
recommended, as it is somewhat nonsensical to add multiple parameters.

Also note *jacobianAddShapeCatalogParameter* as this allows addition
of shape parameters, e.g., pressure broadening coefficients.

Each call to this function adds just a single value to *x*.

Example given the catalog_identity="O2-66 TR UP v1 0 J 1 LO v1 0 J 0",
only the O2 ground-level 119 GHz line can be accessed and only its
catalog_parameter will be accessed.  However, the more lenient
catalog_identity="O2-66 TR UP J 1 LO J 0" may be used, but then the
118 GHz line belonging to v1=1 branch will be added to the same *x*.
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda"},
      .gin = {"catalog_identity", "catalog_parameter"},
      .gin_type = {"QuantumIdentifier", "String"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(The catalog line matching information.)--",
                   R"--(The catalog parameter of the retrieval quantity.)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddBasicCatalogParameters"] = WorkspaceMethodInternalRecord{
      .desc = R"--(See *jacobianAddBasicCatalogParameter*.

This adds a multiple of parameters for first each catalog identity in
``catalog_identities`` and then for each catalog parameter in
``catalog_parameters`` by looping calls to *jacobianAddBasicCatalogParameter*
over these input.
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda"},
      .gin = {"catalog_identities", "catalog_parameters"},
      .gin_type = {"ArrayOfQuantumIdentifier", "ArrayOfString"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(The catalog line matching information.)--",
                   R"--(The catalog parameter of the retrieval quantity.)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddFreqShift"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes a frequency fit of shift type in the Jacobian.

Retrieval of deviations between nominal and actual backend
frequencies can be included by this method. The assumption here is
that the deviation is a constant off-set, a shift, common for all
frequencies (and not varying between measurement blocks).

This method adds one element to the state vector (*x*).
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda", "f_grid"},
      .gin = {"df"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{100e3}},
      .gin_desc = {R"--(Size of perturbation to apply.)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddFreqStretch"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes a frequency fit of stretch type in the Jacobian.

Retrieval of deviations between nominal and actual backend
frequencies can be included by this method. The assumption here is
that the deviation varies linearly over the frequency range
(following ARTS basis function for polynomial order 1).

This method adds one element to the state vector (*x*).
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda", "f_grid"},
      .gin = {"df"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{100e3}},
      .gin_desc = {R"--(Size of perturbation to apply.)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddMagField"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes one magnetic field component in the Jacobian.

The method follows the pattern of other Jacobian methods. The
calculations can only be performed by analytic expressions.

The magnetic field components are retrieved separately, and,
hence, the argument ``component`` can be  "u", "v", "w",
and "strength".

The number of elements added to the state vector (*x*) is::

  n_g1 * n_g2 * n_g3

where n_g1, n_g2 and n_g3 are the length of GIN ``g1``, ``g2`` and ``g3``,
respectively. Here empty vectors should be considered to have a length 1.
The elements are sorted with pressure as innermost loop, followed by
latitude and longitude as outermost loop.

The dB-parameter is only used for Faraday rotation.
)--",
      .author = {"Patrick Eriksson", "Richard Larsson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda"},
      .gin = {"g1", "g2", "g3", "component", "dB"},
      .gin_type = {"Vector", "Vector", "Vector", "String", "Numeric"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    String("v"),
                    Numeric{1.0e-7}},
      .gin_desc = {R"--(Pressure retrieval grid.)--",
                   R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retreival grid.)--",
                   R"--(Magnetic field component to retrieve)--",
                   R"--(Magnetic field perturbation)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddNLTE"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Experimental NLTE Jacobian.

Intention: Adds the nlte_field level distribution per atmospheric grid
to the Jacobian.

The number of elements added to the state vector (*x*) is::

  n_g1 * n_g2 * n_g3

where n_g1, n_g2 and n_g3 are the length of GIN ``g1``, ``g2`` and ``g3``,
respectively. Here empty vectors should be considered to have a length 1.
The elements are sorted with pressure as innermost loop, followed by
latitude and longitude as outermost loop.

The QuantumIdentifier should identify a single energy level, such as:
"H2O-161 EN J 1 Ka 0 Kc 1", for one of the lower levels in the chains
of transitions of water.  Note that using this method directly is not
best practice, as the quantum identifiers of the levels have to be known
at an early stage in NLTE calculations, and will usually populate the
``nlte_level_identifiers`` variable, meaning it is better to use *jacobianAddNLTE*
directly than to individually call this function.
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda"},
      .gin = {"g1", "g2", "g3", "energy_level_identity", "dx"},
      .gin_type =
          {"Vector", "Vector", "Vector", "QuantumIdentifier", "Numeric"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    Numeric{1.0e-3}},
      .gin_desc = {R"--(Pressure retrieval grid.)--",
                   R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retreival grid.)--",
                   R"--(Identifier to the eneregy level)--",
                   R"--(Perturbation of value if required by method)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddNLTEs"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Experimental NLTE Jacobian.  Same as *jacobianAddNLTE* but for
many levels

Adds energy_level_identities.nelem() times as many arguments to *x*
as *jacobianAddNLTE*, ordered as energy_level_identities describes

This method is preferred to *jacobianAddNLTE*, since ``energy_level_identities``
is conveniently almost always the same as ``nlte_level_identifiers``.
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda"},
      .gin = {"g1", "g2", "g3", "energy_level_identities", "dx"},
      .gin_type =
          {"Vector", "Vector", "Vector", "ArrayOfQuantumIdentifier", "Numeric"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    Numeric{1.0e-3}},
      .gin_desc = {R"--(Pressure retrieval grid.)--",
                   R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retreival grid.)--",
                   R"--(Identifiers to the eneregy level)--",
                   R"--(Perturbation of value if required by method)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddPointingZa"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds sensor pointing zenith angle off-set jacobian.

Retrieval of deviations between nominal and actual zenith angle of
the sensor can be included by this method. The weighing functions
can be calculated in several ways:

- ``calcmode = "recalc"``:
    Recalculation of pencil beam spectra,
    shifted with ``dza`` from nominal values. A single-sided
    perturbation is applied (towards higher zenith angles).
- ``calcmode = "interp"``:
    Inter/extrapolation of existing pencil
    beam spectra. For this option, allow some extra margins for
    zenith angle grids, to avoid artifacts when extrapolating
    the data (to shifted zenith angles). The average of a
    negative and a positive shift is taken.

The interp option is recommended. It should in general be both
faster and more accurate (due to the double sided disturbance).
In addition, it is less sensitive to the choice of dza (as long
as a small value is applied).

The pointing off-set can be modelled to be time varying. The time
variation is then described by a polynomial (with standard base
functions). For example, a polynomial order of 0 means that the
off-set is constant in time. If the off-set is totally uncorrelated
between the spectra, set the order to -1.

The number of elements added to the state vector (*x*) is

* if poly_order < 0 : length of *sensor_time*
* otherwise : poly_order+1

In the first case, the order in *x* matches *sensor_time*. In the second
case, the coefficient for polynomial order 0 comes first etc.
)--",
      .author = {"Patrick Eriksson", "Mattias Ekstrom"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities",
             "jacobian_agenda",
             "sensor_pos",
             "sensor_time"},
      .gin = {"poly_order", "calcmode", "dza"},
      .gin_type = {"Index", "String", "Numeric"},
      .gin_value = {Index{0}, String("recalc"), Numeric{0.01}},
      .gin_desc =
          {R"--(Order of polynomial to describe the time variation of pointing off-sets.)--",
           R"--(Calculation method. See above)--",
           R"--(Size of perturbation to apply (when applicable).)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddPolyfit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes polynomial baseline fit in the Jacobian.

This method deals with retrieval of disturbances of the spectra
that can be described by an additive term, a baseline off-set.

The baseline off-set is here modelled as a polynomial. The
polynomial spans the complete frequency range spanned by
*sensor_response_f_grid* and the method should only of interest for
cases with no frequency gap in the spectra. The default assumption
is that the off-set differs between all spectra, but it can also be
assumed that the off-set is common for all e.g. line-of-sights.

If the simulation/retrieval deals with a single spectrum, the number
of elements added to the state vector (*x*) is poly_order+1. The
coefficient for polynomial order 0 comes first etc. The same is true
if ``no_pol_variation``, ``no_los_variation`` and ``no_mblock_variation``
all are set to 1, even if several spectra are involved. Otherwise thenumber of elements added to *x* depends on the number of spectra and
the settings of ``no_pol_variation``, ``no_los_variation`` and 
``no_mblock_variation``. The coefficients of the different polynomial
orders are treated as separate retrieval quantities. That is, the
the elements associated with polynomial order 0 are grouped and form
together a retrieval quantity. The coefficients for higher polynomial
orders are treated in the same way.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities",
             "jacobian_agenda",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid",
             "sensor_pos"},
      .gin = {"poly_order",
              "no_pol_variation",
              "no_los_variation",
              "no_mblock_variation"},
      .gin_type = {"Index", "Index", "Index", "Index"},
      .gin_value = {std::nullopt, Index{0}, Index{0}, Index{0}},
      .gin_desc =
          {R"--(Polynomial order to use for the fit.)--",
           R"--(Set to 1 if the baseline off-set is the same for all Stokes components.)--",
           R"--(Set to 1 if the baseline off-set is the same for all line-of-sights (inside each measurement block).)--",
           R"--(Set to 1 if the baseline off-set is the same for all measurement blocks.)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddScatSpecies"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes a scattering species in the Jacobian.

For 1D or 2D calculations the latitude and/or longitude grid of
the retrieval field should set to have zero length.

The number of elements added to the state vector (*x*) is::

  n_g1 * n_g2 * n_g3

where n_g1, n_g2 and n_g3 are the length of GIN ``g1``, ``g2`` and ``g3``,
respectively. Here empty vectors should be considered to have a length 1.
The elements are sorted with pressure as innermost loop, followed by
latitude and longitude as outermost loop.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda"},
      .gin = {"g1", "g2", "g3", "species", "quantity"},
      .gin_type = {"Vector", "Vector", "Vector", "String", "String"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt},
      .gin_desc =
          {R"--(Pressure retrieval grid.)--",
           R"--(Latitude retrieval grid.)--",
           R"--(Longitude retreival grid.)--",
           R"--(Name of scattering species, must match one element in *scat_species*.)--",
           R"--(Retrieval quantity, e.g. "IWC".)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddShapeCatalogParameter"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Adds a line shape parameter to the Jacobian calculations. These
are constant over all levels so only a single *x*-value is added

Line function parameter assume the derivatives of internal pressure
broadening and line mixing functionality follows a f(T, T0, X0, X1, X2)
format. The shape of the function f() is determined by input
catalog; please see the ARTS documentation for more details.

The input are as follows:

- line_identity:
    Identifier of preferably a single line
- species:
    A SpeciesTag, e.g., "O2" or "H2O" for common species.
    Note that "SELF" and "AIR" tags are used for shape parameters
    affected by self and air-broadening, respectively.
- variable:
    A variable supported by the line, these can be

    - ``"G0"``:  Speed-independent pressure broadening
    - ``"G2"``:  Speed-dependent pressure broadening
    - ``"D0"``:  Speed-independent pressure shift
    - ``"D2"``:  Speed-dependent pressure shift
    - ``"FVC"``: Frequency of velocity changing collisions
    - ``"ETA"``: partial correlation between velocity and rotational state changes due to collisions
    - ``"Y"``:   First order line-mixing parameter
    - ``"G"``:   Second order line-mixing parameter for strength
    - ``"DV"``:  Second order line-mixing parameter for shifting
- coefficient:
    A coefficient in the model to compute the above parameters.

Note that we cannot test if the line in question supports the variable and
coefficient at the level of this function, so many errors will only be reported
at a later stage.

For other spectroscopic parameters, see *jacobianAddBasicCatalogParameter*.
Also see said function for an example of how to set the QuantumIdentifier.
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda"},
      .gin = {"line_identity", "species", "variable", "coefficient"},
      .gin_type = {"QuantumIdentifier", "String", "String", "String"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Line identifier)--",
                   R"--(Species of interest)--",
                   R"--(Variable of interest)--",
                   R"--(Coefficient of interest)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddShapeCatalogParameters"] = WorkspaceMethodInternalRecord{
      .desc = R"--(See *jacobianAddShapeCatalogParameter* for information on
the GIN parameters

This function accepts the same input but for lists of data.
The function loops over each input list
individually and appends the information to *jacobian_quantities*.

Special "ALL" for 1 length ``variables`` and ``coefficients`` are
allowed to compute all variables/coefficients in the order described
in the description of *jacobianAddShapeCatalogParameter*.

For example, if ``line_identities`` have length 5, ``species`` length 4,
``variables`` length 3, and ``coefficients`` length 2, there will be
5*4x3x2 = 120 new additions to *jacobian_quantities* in the order:

- [{line_identities[0], species[0], variables[0] coefficients[0]}]
- [{line_identities[0], species[0], variables[0] coefficients[1]}]
- [{line_identities[0], species[0], variables[1] coefficients[0]}]
- [{line_identities[0], species[0], variables[1] coefficients[1]}]
- [{line_identities[0], species[0], variables[2] coefficients[0]}]
- [{line_identities[0], species[0], variables[2] coefficients[1]}]
- [{line_identities[0], species[1], variables[0] coefficients[0]}]
- ...
- [{line_identities[4], species[3], variables[1] coefficients[1]}]
- [{line_identities[4], species[3], variables[2] coefficients[0]}]
- [{line_identities[4], species[3], variables[2] coefficients[1]}]

or in words: lines first, then species, then variables, then coefficients
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda"},
      .gin = {"line_identities", "species", "variables", "coefficients"},
      .gin_type = {"ArrayOfQuantumIdentifier",
                   "ArrayOfString",
                   "ArrayOfString",
                   "ArrayOfString"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(List of line identifiers)--",
                   R"--(List of species of interest)--",
                   R"--(List of variables of interest)--",
                   R"--(List of coefficients of interest)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddSinefit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes sinusoidal baseline fit in the Jacobian.

Works as *jacobianAddPolyfit*, beside that a series of sine and
cosine terms are used for the baseline fit.

For each value in ``period_lengths`` one sine and one cosine term are
included (in mentioned order). By these two terms the amplitude and
"phase" for each period length can be determined. The sine and
cosine terms have value 0 and 1, respectively, for first frequency.

If the simulation/retrieval deals with a single spectrum, the number
of elements added to the state vector (*x*) is 2 * nperiods, where
nperiods is the length of ``period_lengths``. The same is true
if ``no_pol_variation``, ``no_los_variation`` and ``no_mblock_variation``
all are set to 1, even if several spectra are involved. Otherwise thenumber of elements added to *x* depends on the number of spectra and
the settings of ``no_pol_variation``, ``no_los_variation`` and 
``no_mblock_variation``. The sine and cosine terms for each period
length are treated as a  separate retrieval quantities. That is, the
the elements associated with the first period length are grouped and
form together a retrieval quantity, etc. Inside each retrieval quantity
the pairs of sine and cosine terms are kept together, in given order.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities",
             "jacobian_agenda",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid",
             "sensor_pos"},
      .gin = {"period_lengths",
              "no_pol_variation",
              "no_los_variation",
              "no_mblock_variation"},
      .gin_type = {"Vector", "Index", "Index", "Index"},
      .gin_value = {std::nullopt, Index{0}, Index{0}, Index{0}},
      .gin_desc =
          {R"--(Period lengths of the fit.)--",
           R"--(Set to 1 if the baseline off-set is the same for all Stokes components.)--",
           R"--(Set to 1 if the baseline off-set is the same for all line-of-sights (inside each measurement block).)--",
           R"--(Set to 1 if the baseline off-set is the same for all measurement blocks.)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddSpecialSpecies"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes a special absorption species in the Jacobian.

Similar to *jacobianAddAbsSpecies* but only for number densities.

Species allowed are:

* "electrons"
* "particulates"

Note that the average of all particulates are used to scale its
*jacobian*, so this method works best when only one type of
particulate is being used, i.e., when *scat_data* has only one
scattering species.

The number of elements added to the state vector (*x*) is::

  n_g1 * n_g2 * n_g3

where n_g1, n_g2 and n_g3 are the length of GIN ``g1``, ``g2`` and ``g3``,
respectively. Here empty vectors should be considered to have a length 1.
The elements are sorted with pressure as innermost loop, followed by
latitude and longitude as outermost loop.
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda"},
      .gin = {"g1", "g2", "g3", "species"},
      .gin_type = {"Vector", "Vector", "Vector", "String"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Pressure retrieval grid.)--",
                   R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retreival grid.)--",
                   R"--(The species of the retrieval quantity.)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddSurfaceQuantity"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes a surface quantity in the Jacobian.

The quantity is specified by the GIN-variable ``quantity``. The name
of the quantity must match the name used in *surface_props_names*.

For 1D or 2D calculations the latitude and/or longitude grid of
the retrieval field should set to have zero length.

The number of elements added to the state vector (*x*) is::

  n_g1 * n_g2

where n_g1 and n_g2 are the length of GIN ``g1`` and ``g2``, respectively.
Here empty vectors should be considered to have a length 1.
The elements are sorted with latitude as innermost loop and longitude
as outermost loop.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda"},
      .gin = {"g1", "g2", "quantity"},
      .gin_type = {"Vector", "Vector", "String"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retreival grid.)--",
                   R"--(Retrieval quantity, e.g. "Wind speed".)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddTemperature"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes atmospheric temperatures in the Jacobian.

The calculations are performed by (semi-)analytical expressions.
Hydrostatic equilibrium (HSE) can be included.

The analytical calculation approach neglects so far refraction
totally, but considers the local effect of HSE.
The later should be accaptable for observations around zenith and
nadir. There is no warning if the method is applied incorrectly, 
with respect to these issues. Note that the argument ``hse`` of this
WSM only refers to the Jacobian calculation, if the model and/or
retrieved atmosphere actually fulfils HSE or not is governed in
other manners.

The calculations (both options) assume that gas species are defined
in VMR (a change in temperature then changes the number density). 
This has the consequence that retrieval of temperatures and number
density can not be mixed. Neither any warning here!

The number of elements added to the state vector (*x*) is::

  n_g1 * n_g2 * n_g3

where n_g1, n_g2 and n_g3 are the length of GIN ``g1``, ``g2`` and ``g3``,
respectively. Here empty vectors should be considered to have a length 1.
The elements are sorted with pressure as innermost loop, followed by
latitude and longitude as outermost loop.
)--",
      .author = {"Mattias Ekstrom", "Patrick Eriksson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda"},
      .gin = {"g1", "g2", "g3", "hse"},
      .gin_type = {"Vector", "Vector", "Vector", "String"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, String("on")},
      .gin_desc = {R"--(Pressure retrieval grid.)--",
                   R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retreival grid.)--",
                   R"--(Flag to assume HSE or not ("on" or "off").)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAddWind"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes one atmospheric wind component in the Jacobian.

The method follows the pattern of other Jacobian methods. The
calculations can only be performed by analytic expressions.
Some lower level function depends on frequency perturbations,
however, so therefore a frequency perturbation df is required
and as a consequence *abs_f_interp_order* must be > 0.

The wind field components are retrieved separately, and,
hence, the argument ``component`` can be "u", "v" or "w" 
for vector components, or just "strength" for total wind speed.

The number of elements added to the state vector (*x*) is::

  n_g1 * n_g2 * n_g3

where n_g1, n_g2 and n_g3 are the length of GIN ``g1``, ``g2`` and ``g3``,
respectively. Here empty vectors should be considered to have a length 1.
The elements are sorted with pressure as innermost loop, followed by
latitude and longitude as outermost loop.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities", "jacobian_agenda"},
      .gin = {"g1", "g2", "g3", "component", "dfrequency"},
      .gin_type = {"Vector", "Vector", "Vector", "String", "Numeric"},
      .gin_value =
          {std::nullopt, std::nullopt, std::nullopt, String("v"), Numeric{0.1}},
      .gin_desc = {R"--(Pressure retrieval grid.)--",
                   R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retrieval grid.)--",
                   R"--(Wind component to retrieve)--",
                   R"--(This is the frequency perturbation)--"},
      .pass_workspace = true,

  };

  wsm_data["jacobianAdjustAndTransform"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Applies adjustments and transformations on *jacobian*.

The method handles two tasks:
1. The retrieval transformations set by the user can not be applied
onthe  Jacobian inside *yCalc*. Transformations are instead applied
by calling this method.
2. It applies required adjustments of the Jacoboan. So far there is
only one possible adjustment. If any absorption species uses the "rel"
unit, an adjustment is needed for later iterations of the inversion.

If no tranformations are selected and the "rel" option is not used at
all, there is no need to call this method(, but you can still include it
without causing any error, the calculations will just be a bit slower).
Otherwise, this method should be called, typically as part of
*inversion_iterate_agenda*.

The method accepts if *jacobian* is empty, and then does, nothing.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian"},

      .in = {"jacobian", "jacobian_quantities", "x"},

  };

  wsm_data["jacobianCalcDoNothing"] = WorkspaceMethodInternalRecord{
      .desc = R"--(This function doesn't do anything. It just exists to satisfy
the input and output requirement of the *jacobian_agenda*.

This method is added to *jacobian_agenda* by *jacobianAddAbsSpecies*
and some similar methods, and it should normally not be called by
the user.
)--",
      .author = {"Oliver Lemke"},
      .out = {"jacobian"},

      .in = {"jacobian", "mblock_index", "iyb", "yb"},

  };

  wsm_data["jacobianCalcFreqShift"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates frequency shift jacobians by interpolation
of *iyb*.

This function is added to *jacobian_agenda* by jacobianAddFreqShift
and should normally not be called by the user.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian"},

      .in = {"jacobian",
             "mblock_index",
             "iyb",
             "yb",
             "f_grid",
             "mblock_dlos",
             "sensor_response",
             "jacobian_quantities"},

  };

  wsm_data["jacobianCalcFreqStretch"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates frequency stretch jacobians by interpolation
of *iyb*.

This function is added to *jacobian_agenda* by jacobianAddFreqStretch
and should normally not be called by the user.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian"},

      .in = {"jacobian",
             "mblock_index",
             "iyb",
             "yb",
             "f_grid",
             "mblock_dlos",
             "sensor_response",
             "sensor_response_pol_grid",
             "sensor_response_f_grid",
             "sensor_response_dlos_grid",
             "jacobian_quantities"},

  };

  wsm_data["jacobianCalcPointingZaInterp"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates zenith angle pointing deviation jacobians by
inter-extrapolation of *iyb*.

This function is added to *jacobian_agenda* by
jacobianAddPointingZa and should normally not be
called by the user.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian"},

      .in = {"jacobian",
             "mblock_index",
             "iyb",
             "yb",
             "f_grid",
             "sensor_los",
             "mblock_dlos",
             "sensor_response",
             "sensor_time",
             "jacobian_quantities"},

  };

  wsm_data["jacobianCalcPointingZaRecalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates zenith angle pointing deviation jacobians by
recalulation of *iyb*.

This function is added to *jacobian_agenda* by
jacobianAddPointingZa and should normally not be
called by the user.
)--",
      .author = {"Mattias Ekstrom", "Patrick Eriksson"},
      .out = {"jacobian"},

      .in = {"jacobian",
             "mblock_index",
             "iyb",
             "yb",
             "atm_field",
             "cloudbox_on",
             "f_grid",
             "sensor_pos",
             "sensor_los",
             "transmitter_pos",
             "mblock_dlos",
             "sensor_response",
             "sensor_time",
             "iy_unit",
             "iy_main_agenda",
             "jacobian_quantities"},

      .pass_workspace = true,

  };

  wsm_data["jacobianCalcPolyfit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates jacobians for polynomial baseline fit.

This function is added to *jacobian_agenda* by jacobianAddPolyfit
and should normally not be called by the user.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian"},

      .in = {"jacobian",
             "mblock_index",
             "iyb",
             "yb",
             "sensor_response",
             "sensor_response_pol_grid",
             "sensor_response_f_grid",
             "sensor_response_dlos_grid",
             "jacobian_quantities"},
      .gin = {"poly_coeff"},
      .gin_type = {"Index"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Polynomial coefficient to handle.)--"},

  };

  wsm_data["jacobianCalcSinefit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates jacobians for sinusoidal baseline fit.

This function is added to *jacobian_agenda* by jacobianAddPolyfit
and should normally not be called by the user.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian"},

      .in = {"jacobian",
             "mblock_index",
             "iyb",
             "yb",
             "sensor_response",
             "sensor_response_pol_grid",
             "sensor_response_f_grid",
             "sensor_response_dlos_grid",
             "jacobian_quantities"},
      .gin = {"period_index"},
      .gin_type = {"Index"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(Index among the period length specified for add-method.)--"},

  };

  wsm_data["jacobianClose"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Closes the array of retrieval quantities and prepares for
calculation of the Jacobian matrix.

This function closes the *jacobian_quantities* array and sets
*jacobian_do* to 1.

Retrieval quantities should not be added after a call to this WSM.
No calculations are performed here.
)--",
      .author = {"Mattias Ekstrom"},
      .out = {"jacobian_do", "jacobian_agenda"},

      .in = {"jacobian_agenda", "jacobian_quantities"},

      .pass_workspace = true,

  };

  wsm_data["jacobianFromTwoY"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets *jacobian* based on the difference vetween two measurement vectors.

This function assumes that ``y_pert`` contains a measurement calculated
with some variable perturbed, in comparison to the calculation
behind *y*. The function takes the differences between ``y_pert``
and *y* to form a numerical derived estimate of *jacobian*.
This gives a Jacobian wit a single column.

*jacobian* equals here: (y_pert-y)/pert_size.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian"},

      .in = {"y"},
      .gin = {"y_pert", "pert_size"},
      .gin_type = {"Vector", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Perturbed measurement vector)--",
                   R"--(Size of perturbation behind spectra in *ybatch*.)--"},

  };

  wsm_data["jacobianFromYbatch"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *jacobian* based on perturbation calcuations.

This function assumes that *ybatch* contains spectra calculated
with some variable perturbed, in comparison to the calculation
behind *y*. The function takes the differences between *ybatch*
and *y* to form a numerical derived estimate of *jacobian*.

Column i of *jacobian* equals: (ybatch[i]-y)/pert_size.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian"},

      .in = {"ybatch", "y"},
      .gin = {"pert_size"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Size of perturbation behind spectra in *ybatch*.)--"},

  };

  wsm_data["jacobianInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initialises the variables connected to the Jacobian matrix.

This function initialises the *jacobian_quantities* array so
that retrieval quantities can be added to it. Accordingly, it has
to be called before any calls to jacobianAddTemperature or
similar methods.

The Jacobian quantities are initialised to be empty.
)--",
      .author = {"Mattias Ekstrom"},
      .out = {"jacobian_quantities", "jacobian_agenda"},

      .pass_workspace = true,

  };

  wsm_data["jacobianOff"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Makes mandatory initialisation of some jacobian variables.

Some clear-sky jacobian WSVs must be initialised even if no such
calculations will be performed.  This is handled with this method.
That is, this method must be called when no clear-sky jacobians
will be calculated (even if cloudy-sky jacobians are calculated!).

Sets *jacobian_do* to 0.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"jacobian_do", "jacobian_agenda", "jacobian_quantities"},

      .pass_workspace = true,

  };

  wsm_data["jacobianSetAffineTransformation"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds an affine transformation of the last element of
*jacobian_quantities*.

See *jacobianSetFuncTransformation* for  a general description of how
retrieval transformations are defined. Transformations are not applied by
methods such as *yCalc*. Instead, the method *jacobianAdjustAndTransform*
must be called to activate the transformations.

The affine transformation is specified by a transformation matrix, A,
and an offset vector, b. These two are applied as described in
*jacobianSetFuncTransformation*.

The transformations is applied as::

  x = A * ( z - b )

where z is the retrieval quantity on the standard retrieval grids
and x is the final state vector.

So far, the following must be true for valid A-matrices::

  z = A' * x + b

That is, the reversed transformation is given by A transposed.

This method must only be called if an affine transformation is wanted.
Default is to make no such tranformation at all.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"jacobian_quantities"},

      .in = {"jacobian_quantities"},
      .gin = {"transformation_matrix", "offset_vector"},
      .gin_type = {"Matrix", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(The transformation matrix A)--",
                   R"--(The offset vector b)--"},

  };

  wsm_data["jacobianSetFuncTransformation"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the functional transformation of the last element of
*jacobian_quantities*.

See below for a general description of how retrieval transformations
are defined. Transformations are not applied by methods such as *yCalc*.
Instead, the method *jacobianAdjustAndTransform* must be called to
activate the transformations.

The following transformations can be selected (by ``transformation_func``):

- ``"log"``: The natural logarithm
- ``"log10"``: The base-10 logarithm
- ``"atanh"``: Area hyperbolic tangent 
- ``"none"``: No transformation at all

This method needs only to be called if a functional transformation
is wanted. Default is to make no such tranformation at all (i.e.
the option "none" exists only for reasons of flexibility).

The log-options are applied as log(z-z_min) and log10(z-z_min).
The default for ``z_min`` is zero, but by changing it the lower limit
for z can be changed. Note that ``z_min`` becomes the lower limit for
allowed values of z. The GIN ``z_max`` is here ignored.

For the atanh-option, also ``z_max`` is considered. This transformation
is applied as atanh((2(z-z_min)/(z_max-z_min))-1). As above,``z_min``
is lower limit for allowed values of z. On the other hand, ``z_max``
eines the upper limit for z.

The GIN ``transformation_func`` is so far only used for atanh. The parameter
specifies the maximum allowed value allowed for u. That is, the valid
range for u becomes ]0,tfunc_parameter[. Note that log and log10
demands/ensures that u > 0, but implies no upper limit.

General handling of retrieval units and transformations:
---
Default is that quantities are retrieved as defined in ARTS, but
both some unit conversion and transformations are provided. These
operations are applied as::

  x = A * ( f(u(z)) - b ) 

where

- z is the quantity as defined ARTS
- u represents the change of unit
- f is the transformation function
- A and b define together an affine transformation
- x is the retrieved quantity

For example, this systen allows to retrive a principal component
representation (A and b) of the log (f) of relative humidity (u).

Change of unit is selected by the quantity specific jacobian-add
methods (so far only at hand for gas species). 

Activating a transformation function is done by this method. Note
that the functions are defined as the transformation from z to x.
For more details on affine transformations, see
*jacobianSetAffineTransformation*.
)--",
      .author = {"Patrick Eriksson", "Simon Pfreundschuh"},
      .out = {"jacobian_quantities"},

      .in = {"jacobian_quantities"},
      .gin = {"transformation_func", "z_min", "z_max"},
      .gin_type = {"String", "Numeric", "Numeric"},
      .gin_value = {std::nullopt, Numeric{0}, Numeric{-99e99}},
      .gin_desc = {R"--(The transformation function.)--",
                   R"--(Lower limit of z.)--",
                   R"--(Upper limit of z.)--"},

  };

  wsm_data["jacobian_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *jacobian_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["lbl_checkedCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Checks that the line-by-line parameters are OK.

On failure, will throw.  On success, lbl_checked evals as true

Note that checks may become more stringent as ARTS evolves, especially for
"new" options.  This test might succeed in one version of ARTS but fail
in later versions
)--",
      .author = {"Richard Larsson"},
      .out = {"lbl_checked"},

      .in = {"abs_lines_per_species", "abs_species", "isotopologue_ratios"},

  };

  wsm_data["mblock_dlosFrom1dAntenna"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *mblock_dlos* based on a 1D gaussian antenna response.

The length of *mblock_dlos* is determined by ``npoints``. The end
points of the grid are set to be the same as for the antenna
response. The spacing of the grid follows the magnitude of the
response; the spacing is smaller where the response is high.
More precisely, the grid points are determined by dividing the
cumulative sum of the response in equal steps.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"mblock_dlos"},

      .in = {"antenna_response"},
      .gin = {"npoints"},
      .gin_type = {"Index"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Number of points (>1) to include in *mblock_dlos*.)--"},

  };

  wsm_data["mc_antennaSetGaussian"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Makes mc_antenna (used by MCGeneral) a 2D Gaussian pattern.

The gaussian antenna pattern is determined by ``za_sigma`` and
``aa_sigma``, which represent the standard deviations in the
uncorrelated bivariate normal distribution.
)--",
      .author = {"Cory Davis"},
      .out = {"mc_antenna"},

      .gin = {"za_sigma", "aa_sigma"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(Width in the zenith angle dimension as described above.)--",
           R"--(Width in the azimuth angle dimension as described above.)--"},

  };

  wsm_data["mc_antennaSetGaussianByFWHM"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Makes mc_antenna (used by MCGeneral) a 2D Gaussian pattern.

The gaussian antenna pattern is determined by ``za_fwhm`` and
``aa_fwhm``, which represent the full width half maximum (FWHM)
of the antenna response, in the zenith and azimuthal planes.
)--",
      .author = {"Cory Davis"},
      .out = {"mc_antenna"},

      .gin = {"za_fwhm", "aa_fwhm"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(Width in the zenith angle dimension as described above.)--",
           R"--(Width in the azimuth angle dimension as described above.)--"},

  };

  wsm_data["mc_antennaSetPencilBeam"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Makes mc_antenna (used by MCGeneral) a pencil beam.

This WSM makes the subsequent MCGeneral WSM perform pencil beam
RT calculations.
)--",
      .author = {"Cory Davis"},
      .out = {"mc_antenna"},

  };

  wsm_data["met_profile_calc_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *met_profile_calc_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--",
      .author = {"Richard Larsson"},
      .out = {"met_profile_calc_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["nbooksGet"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Retrieve nbooks from given variable and store the value in the
workspace variable *nbooks*
)--",
      .author = {"Oliver Lemke"},
      .out = {"nbooks"},

      .gin = {"v"},
      .gin_type = {"Tensor4, Tensor5, Tensor6, Tensor7"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Variable to get the number of books from.)--"},

  };

  wsm_data["ncolsGet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Retrieve ncols from given variable and store the value in the
workspace variable *ncols*
)--",
      .author = {"Oliver Lemke"},
      .out = {"ncols"},

      .gin = {"v"},
      .gin_type =
          {"Matrix, Sparse, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Variable to get the number of columns from.)--"},

  };

  wsm_data["nelemGet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Retrieve nelem from given variable and store the value in the
variable *nelem*.
)--",
      .author = {"Oliver Lemke"},
      .out = {"nelem"},

      .gin = {"v"},
      .gin_type =
          {"ArrayOfAbsorptionLines, ArrayOfAgenda, ArrayOfArrayOfAbsorptionLines, ArrayOfArrayOfGriddedField1, ArrayOfArrayOfGriddedField2, ArrayOfArrayOfGriddedField3, ArrayOfArrayOfIndex, ArrayOfArrayOfMatrix, ArrayOfArrayOfMuelmatMatrix, ArrayOfArrayOfMuelmatVector, ArrayOfArrayOfPropmatMatrix, ArrayOfArrayOfPropmatVector, ArrayOfArrayOfScatteringMetaData, ArrayOfArrayOfSingleScatteringData, ArrayOfArrayOfSpeciesTag, ArrayOfArrayOfStokvecMatrix, ArrayOfArrayOfStokvecVector, ArrayOfArrayOfString, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6, ArrayOfArrayOfTime, ArrayOfArrayOfVector, ArrayOfAtmPoint, ArrayOfCIARecord, ArrayOfGriddedField1, ArrayOfGriddedField2, ArrayOfGriddedField3, ArrayOfGriddedField4, ArrayOfIndex, ArrayOfJacobianTarget, ArrayOfMatrix, ArrayOfMuelmatMatrix, ArrayOfMuelmatVector, ArrayOfPpath, ArrayOfPropmatMatrix, ArrayOfPropmatVector, ArrayOfQuantumIdentifier, ArrayOfRetrievalQuantity, ArrayOfScatteringMetaData, ArrayOfSingleScatteringData, ArrayOfSparse, ArrayOfSpeciesTag, ArrayOfStokvecMatrix, ArrayOfStokvecVector, ArrayOfString, ArrayOfSun, ArrayOfTelsemAtlas, ArrayOfTensor3, ArrayOfTensor4, ArrayOfTensor5, ArrayOfTensor6, ArrayOfTensor7, ArrayOfTime, ArrayOfVector, ArrayOfXsecRecord, Vector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Variable to get the number of elements from.)--"},

  };

  wsm_data["nlibrariesGet"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Retrieve nlibraries from given variable and store the value in the
workspace variable *nlibraries*
)--",
      .author = {"Oliver Lemke"},
      .out = {"nlibraries"},

      .gin = {"v"},
      .gin_type = {"Tensor7"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Variable to get the number of libraries from.)--"},

  };

  wsm_data["nlteOff"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Disable Non-LTE calculations.
)--",
      .author = {"Oliver Lemke"},
      .out = {"nlte_do"},

  };

  wsm_data["npagesGet"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Retrieve npages from given variable and store the value in the
workspace variable *npages*
)--",
      .author = {"Oliver Lemke"},
      .out = {"npages"},

      .gin = {"v"},
      .gin_type = {"Tensor3, Tensor4, Tensor5, Tensor6, Tensor7"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Variable to get the number of pages from.)--"},

  };

  wsm_data["nrowsGet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Retrieve nrows from given variable and store the value in the
workspace variable *nrows*
)--",
      .author = {"Oliver Lemke"},
      .out = {"nrows"},

      .gin = {"v"},
      .gin_type =
          {"Matrix, Sparse, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Variable to get the number of rows from.)--"},

  };

  wsm_data["nshelvesGet"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Retrieve nshelves from given variable and store the value in the
workspace variable *nshelves*
)--",
      .author = {"Oliver Lemke"},
      .out = {"nshelves"},

      .gin = {"v"},
      .gin_type = {"Tensor5, Tensor6, Tensor7"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Variable to get the number of shelves from.)--"},

  };

  wsm_data["nvitrinesGet"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Retrieve nvitrines from given variable and store the value in the
workspace variable *nvitrines*
)--",
      .author = {"Oliver Lemke"},
      .out = {"nvitrines"},

      .gin = {"v"},
      .gin_type = {"Tensor6, Tensor7"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Variable to get the number of vitrines from.)--"},

  };

  wsm_data["opt_prop_bulkCalc"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculates bulk absorption extinction at one atmospheric grid point.

This WSM sums up the monochromatic absorption vectors and
extinction matrices of all scattering elements (*abs_vec_spt* and
*ext_mat_spt*, respectively) weighted by their respective
particle number density given by *pnd_field*, for a single location
within the cloudbox, given by *scat_p_index*, *scat_lat_index*, and
*scat_lon_index*.
The resulting  extinction matrix is added to the workspace variable
*ext_mat*.
)--",
      .author = {"Jana Mendrok, Sreerekha T.R."},
      .out = {"ext_mat", "abs_vec"},

      .in = {"ext_mat",
             "abs_vec",
             "ext_mat_spt",
             "abs_vec_spt",
             "pnd_field",
             "scat_p_index",
             "scat_lat_index",
             "scat_lon_index"},

  };

  wsm_data["opt_prop_sptFromData"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculates monochromatic optical properties for all scattering
elements.

In this function the extinction matrix and the absorption vector
are calculated in the laboratory frame. An interpolation of the
data on the actual frequency is the first step in this function.
The next step is a transformation from the database coordinate
system to the laboratory coordinate system.

Output of the function are *ext_mat_spt* and *abs_vec_spt*, which
hold the optical properties for a specified propagation direction
for each scattering element.
)--",
      .author = {"Claudia Emde"},
      .out = {"ext_mat_spt", "abs_vec_spt"},

      .in = {"ext_mat_spt",
             "abs_vec_spt",
             "scat_data",
             "za_grid",
             "aa_grid",
             "za_index",
             "aa_index",
             "f_index",
             "f_grid",
             "rtp_temperature",
             "pnd_field",
             "scat_p_index",
             "scat_lat_index",
             "scat_lon_index"},

  };

  wsm_data["opt_prop_sptFromMonoData"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates optical properties for the scattering elements.

As *opt_prop_sptFromData* but no frequency interpolation is
performed. The single scattering data is here obtained from
*scat_data_mono*, instead of *scat_data*.
)--",
      .author = {"Cory Davis"},
      .out = {"ext_mat_spt", "abs_vec_spt"},

      .in = {"ext_mat_spt",
             "abs_vec_spt",
             "scat_data_mono",
             "za_grid",
             "aa_grid",
             "za_index",
             "aa_index",
             "rtp_temperature",
             "pnd_field",
             "scat_p_index",
             "scat_lat_index",
             "scat_lon_index"},

  };

  wsm_data["opt_prop_sptFromScat_data"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Derives monochromatic optical properties for all scattering
elements.

As *opt_prop_sptFromData*, but using frequency pre-interpolated
data (as produced by *scat_dataCalc*), i.e. in here no frequency
interpolation is done anymore.
)--",
      .author = {"Jana Mendrok, Claudia Emde"},
      .out = {"ext_mat_spt", "abs_vec_spt"},

      .in = {"ext_mat_spt",
             "abs_vec_spt",
             "scat_data",
             "scat_data_checked",
             "za_grid",
             "aa_grid",
             "za_index",
             "aa_index",
             "f_index",
             "rtp_temperature",
             "pnd_field",
             "scat_p_index",
             "scat_lat_index",
             "scat_lon_index"},

  };

  wsm_data["output_file_formatSetAscii"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the output file format to ASCII.
)--",
      .author = {"Oliver Lemke"},
      .out = {"output_file_format"},

  };

  wsm_data["output_file_formatSetBinary"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the output file format to binary.
)--",
      .author = {"Oliver Lemke"},
      .out = {"output_file_format"},

  };

  wsm_data["output_file_formatSetZippedAscii"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the output file format to zipped ASCII.
)--",
      .author = {"Oliver Lemke"},
      .out = {"output_file_format"},

  };

  wsm_data["particle_bulkpropRadarOnionPeeling"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Inverts radar reflectivities by in an onion peeling manner.

The method assumes space-based measurements and invert one altitude
at the time, based on a pre-calculated inversion table (``invtable``)
and starting at the top of the atmosphere. If attenuation is
completely ignored, the table is effectively used as a look-up table
to map dBZe to hydrometeor values. The method considers attenuation
by default, where extinction due to hydrometeors is taken from the
table and the one due to *abs_species* is obtained by
*propmat_clearsky_agenda*.

The inversion table consists of two GriddedField3. The first field
shall match liquid hydrometeors and is applied for temperatures above
``t_phase``. The second field is applied for lower temperatures and
shall thus correspond to ice hydrometeors.

The size of each field is (2,ndb,nt). The two page dimensions match
the hydrometeor property to retrieve and extinction, respectively.
The table shall hold the 10-logarithm of the property, such as
log10(IWC). ndb is the number of dBZe values in the table and nt
the number of temperatures. The table is interpolated in temperature
in a nearest neighbour fashion, while in a linear interpolation is
applied in the dBZe dimension.

The field of radar reflectivities (``dBZe``) shall cover the complete
atmosphere and then match e.g. ``t_field`` in size. The observation
geometry is here specified by giving the incidence angle for each
profile of dBZe values (by ``incangles``). A flat Earth approximation
is applied inside the method.

All values below ``dbze_noise`` are treated as pure noise and
``particle_bulkprop_field`` is set to zero for these positions.
The comparison to ``dbze_noise`` is done with uncorrected values.

Further, all values at altitudes below z_surface + h_clutter are
assumed to be surface clutter and are rejected. If ``fill_clutter``
is set to 1, the retrieval just above the clutter zone is assumed
valid also below and is copied to all altitudes below (also for
altitudes below the surface).

Unfiltered clutter can cause extremely high retrived water contents.
The GIN ``wc_max`` defines an upper limit for reasonable water contents.
Retrievals ending up above this value are set to zero. Values below
``wc_max`` but above ``wc_clip``, are set to ``wc_clip``.

Significant radar echos (>dbze_noise and above clutter zone) are
assumed to match liquid hydrometeors for temperatures >= ``t_phase``
and ice ones for lower temperatures.

Default is to consider attenuation of both hydrometeors and absorption
species. These two sources to attenuation can be ignored by setting
``do_atten_hyd`` and ``do_atten_abs`` to zero, respectively.

Default is to consider hydrometeor attenuation, but there could be
two reasons to ignore it. It can cause a "run away" effect in the
retrievals. Ignoring it can also compensate for impact of multiple
scattering in space-based observations, as shown by: Matrosov and
Battaglia, GRL, 2009. However, ignoring the hydrometeor attenuation
totally gives a too high compensating effect and the GIN
``atten_hyd_scaling`` allows to test intermediate compensations. This
GIN matches the GIN pext_scaling of *iyRadarSingleScat*, but they
have different default values. The default in this method follows the
results for CloudSat in Matrosov and Battaglia. Please note that
``do_atten_hyd`` must be true to apply ``atten_hyd_scaling``.

Even with ``atten_hyd_scaling`` below 1, there could be a run-away in
the estimated attenuation, and ``atten_hyd_max`` stops this by setting
a maximum value to the hydrometeor attenuation.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"atm_field", "particle_bulkprop_names"},

      .in = {"atm_field",
             "surface_field",
             "atmfields_checked",
             "atmgeom_checked",
             "f_grid",
             "propmat_clearsky_agenda",
             "scat_species"},
      .gin = {"invtable",
              "incangles",
              "dBZe",
              "dbze_noise",
              "h_clutter",
              "fill_clutter",
              "t_phase",
              "wc_max",
              "wc_clip",
              "do_atten_abs",
              "do_atten_hyd",
              "atten_hyd_scaling",
              "atten_hyd_max"},
      .gin_type = {"ArrayOfGriddedField3",
                   "Matrix",
                   "Tensor3",
                   "Numeric",
                   "Matrix",
                   "Index",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Index",
                   "Index",
                   "Numeric",
                   "Numeric"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    Numeric{-99},
                    std::nullopt,
                    Index{0},
                    Numeric{273.15},
                    Numeric{10e-3},
                    Numeric{5e-3},
                    Index{1},
                    Index{1},
                    Numeric{0.5},
                    Numeric{3}},
      .gin_desc =
          {R"--(Inversion table, see above.)--",
           R"--(Incidence angles.)--",
           R"--(Field of radar reflectivities, in dBZe.)--",
           R"--(Noise level. See above.)--",
           R"--(Height of clutter zone. Either same size as ``z_surface`` or a single value. In the later case, that value is applied at all positions.)--",
           R"--(Flag to fill clutter zone, by copying retrieval just above it.)--",
           R"--(Phase boundary temperature. See above.)--",
           R"--(Max reasonable water content)--",
           R"--(Clip value for water content retrievals.)--",
           R"--(Flag to consider attenuation due to hydrometeors.)--",
           R"--(Flag to consider attenuation due to absorption species.)--",
           R"--(Hydrometeor attenuation scaling factor.)--",
           R"--(Hydrometeor attenuation not allowed to pass this value [dB].)--"},
      .pass_workspace = true,

  };

  wsm_data["particle_bulkprop_fieldClip"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Clipping of ``particle_bulkprop_field``.

The method allows you to apply hard limits the values of
``particle_bulkprop_field``. All values, of the property selected,
below ``limit_low``, are simply set to ``limit_low``. And the same
is performed with respect to ``limit_high``. That is, the data in x
for the retrieval quantity are forced to be inside the range
[limit_low,limit_high].

Setting species="ALL", is a shortcut for applying the limits on all
properties.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"atm_field"},

      .in = {"atm_field"},
      .gin = {"bulkprop_name", "limit_low", "limit_high"},
      .gin_type = {"String", "Numeric", "Numeric"},
      .gin_value = {std::nullopt, -std::numeric_limits<Numeric>::infinity(), std::numeric_limits<Numeric>::infinity()},
      .gin_desc = {R"--(Name of bulk property to consider, or "ALL".)--",
                   R"--(Lower limit for clipping.)--",
                   R"--(Upper limit for clipping.)--"},

  };

  wsm_data["particle_fieldCleanup"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Removes unrealistically small or erroneous data from particle fields.

This WSM checks if the input particle field (e.g.
``particle_bulkprop_field``) contains values
smaller than the given ``threshold``. In this case, these values will
be set to zero.

The method should be applied if the particle fields contain
unrealistically small or erroneous data (NWP/GCM model data, e.g.
from the Chevallierl_91l sets, often contain very small or even
negative values, which are numerical artefacts rather than physical
values.)
For the scat_species_XXX_fields, it needs to be applied separately
per Tensor4 type field collection. This allows to use different
thresholds for the different types of fields (not for the different
scattering species, though).

*particle_fieldCleanup* shall be called after generation of the
atmopheric fields.
)--",
      .author = {"Daniel Kreyling"},

      .gout = {"particle_field_out"},
      .gout_type = {"Tensor4"},
      .gout_desc =
          {R"--(A particle property field, e.g. ``particle_bulkprop_field``)--"},

      .gin = {"particle_field_in", "threshold"},
      .gin_type = {"Tensor4", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(A particle property field, e.g. ``particle_bulkprop_field``)--",
           R"--(Threshold below which the ``particle_field`` values are set to zero.)--"},

  };

  wsm_data["particle_massesFromMetaData"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Derives *particle_masses* from *scat_meta*.

It extracts the mass information of the scattering elements
from *scat_meta*. Each scattering species is taken as a
separate category of particle_masses, i.e., the resulting
*particle_masses* matrix will contain as many columns as
scattering species are present in *scat_meta*.
)--",
      .author = {"Jana Mendrok"},
      .out = {"particle_masses"},

      .in = {"scat_meta"},

  };

  wsm_data["particle_massesFromMetaDataSingleCategory"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Sets *particle_masses* based on *scat_meta* assuming
all particles are of the same mass category.

This method derives the particle masses from the mass entry
of each scattering element. It is assumed that all scattering
elements represent particles of the same (bulk) matter
(e.g. water or ice). With other words, a single mass category
is assumed (see *particle_masses* for a definition of "mass
category").

If just having clouds, the resulting mass category can be seen as
the total cloud water content, with possible contribution from
both ice and liquid phase.
)--",
          .author = {"Jana Mendrok", "Patrick Eriksson"},
          .out = {"particle_masses"},

          .in = {"scat_meta"},

      };

  wsm_data["pha_matCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates the total phase matrix of all scattering elements.

This function sums up the monochromatic phase matrices of all
scattering elements *pha_mat_spt* weighted with  their respective
particle number density, given by *pnd_field*, for a single location
within the cloudbox, given by *scat_p_index*, *scat_lat_index*, and
*scat_lon_index*.
)--",
      .author = {"Sreerekha T.R."},
      .out = {"pha_mat"},

      .in = {"pha_mat_spt",
             "pnd_field",
             "scat_p_index",
             "scat_lat_index",
             "scat_lon_index"},

  };

  wsm_data["pha_mat_sptFromData"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculation of the phase matrix of the individual scattering elements.

This function can be used in *pha_mat_spt_agenda* as part of
the calculation of the scattering integral.

First, data at the requested frequency (given by *f_grid* and
*f_index*) and temperature (given by *rtp_temperature*) is
extracted. This is followed by a transformation from the database
coordinate system to the laboratory coordinate system.

Frequency extraction is always done by (linear) interpolation.
Temperature is (linearly) interpolated when at least two
temperature grid points are present in the *scat_data* and
*rtp_temperature* is positive. If only a single temperature point
is available, data for this point is used without modification. In
order to speed up calculations, temperature interpolation can be
avoided by passing a *rtp_temperature* < 0. In this case, a specific
temperature grid from the *scat_data* grid is used without
modification. The selection is as follows:

- -10 < *rtp_temperature* <   0   T_grid[0]     lowest temperature
- -20 < *rtp_temperature* < -10   T_grid[nT-1]  highest temperature
- *rtp_temperature* < -20   T_grid[nT/2]  median grid point
)--",
      .author = {"Claudia Emde"},
      .out = {"pha_mat_spt"},

      .in = {"pha_mat_spt",
             "scat_data",
             "za_grid",
             "aa_grid",
             "za_index",
             "aa_index",
             "f_index",
             "f_grid",
             "rtp_temperature",
             "pnd_field",
             "scat_p_index",
             "scat_lat_index",
             "scat_lon_index"},

  };

  wsm_data["pha_mat_sptFromDataDOITOpt"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculation of the phase matrix of the individual scattering elements.

In this function the phase matrix is extracted from
*pha_mat_sptDOITOpt*. It can be used in the agenda
*pha_mat_spt_agenda*. This method must be used in combination with
*DoitScatteringDataPrepare*.

Temperature is considered as described for *pha_mat_sptFromData*
)--",
      .author = {"Claudia Emde"},
      .out = {"pha_mat_spt"},

      .in = {"pha_mat_spt",
             "pha_mat_sptDOITOpt",
             "scat_data_mono",
             "doit_za_grid_size",
             "aa_grid",
             "za_index",
             "aa_index",
             "rtp_temperature",
             "pnd_field",
             "scat_p_index",
             "scat_lat_index",
             "scat_lon_index"},

  };

  wsm_data["pha_mat_sptFromMonoData"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculation of the phase matrix of the individual scattering elements.

This function is the monochromatic version of *pha_mat_sptFromData*.
)--",
      .author = {"Claudia Emde"},
      .out = {"pha_mat_spt"},

      .in = {"pha_mat_spt",
             "scat_data_mono",
             "doit_za_grid_size",
             "aa_grid",
             "za_index",
             "aa_index",
             "rtp_temperature",
             "pnd_field",
             "scat_p_index",
             "scat_lat_index",
             "scat_lon_index"},

  };

  wsm_data["pha_mat_sptFromScat_data"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculation of the phase matrix of the individual scattering elements.

As *pha_mat_sptFromData*, but using frequency pre-interpolated
data (as produced by *scat_dataCalc*), i.e. in here no frequency
interpolation is done anymore.
)--",
      .author = {"Jana Mendrok, Claudia Emde"},
      .out = {"pha_mat_spt"},

      .in = {"pha_mat_spt",
             "scat_data",
             "scat_data_checked",
             "za_grid",
             "aa_grid",
             "za_index",
             "aa_index",
             "f_index",
             "rtp_temperature",
             "pnd_field",
             "scat_p_index",
             "scat_lat_index",
             "scat_lon_index"},

  };

  wsm_data["pha_mat_spt_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *pha_mat_spt_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--",
      .author = {"Richard Larsson"},
      .out = {"pha_mat_spt_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["pndFromPsd"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculates *pnd_data* from given *psd_data* for one scattering species.

Performs integration of the size distribution over the size grid
bin deriving pnd (units #/m3) from psd (units #/m3/m). Some checks
on the sufficiency of the size grid range and coverage are applied.

``quad_order`` can be 0 for rectangular or 1 for trapezoidal
integration. The only difference is the treatment of the start and
end nodes. For trapezoidal their corresponding bins end exactly at
the nodes, while for rectangular they extend further out by the half
distance to the neighbor node (but not beyond 0).

Attempts to check that the size grids and *scat_data* represent the
bulk extinction sufficiently. Specifically, it is tested that

(a) psd*ext is decreasing at the small and large particle size
    ends of the size grid - but only if scattering species bulk
    extinction exceeds 1% of ``threshold_ss_ext``.
(b) removing the smallest and largest particles changes the
    resulting bulk extinction by less then a fraction of
    ``threshold_se_ext`` - but only if scattering species bulk
    extinction exceeds ``threshold_ss_ext`` and number density (pnd)
    of the edge size point at this atmospheric level is larger
    than ``threshold_se_pnd`` times the maximum pnd of this
    scattering element over all atmospheric levels.

Skipping tests in case of low extinction is done in order to
minimize issues arising from very low mass densities,
particularly at single atmospheric levels, and very low bulk
extinctions, i.e. in cases where the effects on the radiance fields
are estimated to be low.
NOTE: The tests are only approximate and do not guarantee the
validity of the resulting bulk properties (and increasing the
thresholds will decrease the reliability of the bulk properties).
)--",
      .author = {"Jana Mendrok, Patrick Eriksson"},
      .out = {"pnd_data", "dpnd_data_dx"},

      .in = {"pnd_size_grid",
             "psd_data",
             "psd_size_grid",
             "dpsd_data_dx",
             "scat_data",
             "f_grid",
             "scat_data_checked"},
      .gin = {"quad_order",
              "scat_index",
              "threshold_se_ext",
              "threshold_ss_ext",
              "threshold_se_pnd"},
      .gin_type = {"Index", "Index", "Numeric", "Numeric", "Numeric"},
      .gin_value =
          {Index{1}, std::nullopt, Numeric{0.02}, Numeric{1e-8}, Numeric{0.02}},
      .gin_desc =
          {R"--(Order of bin quadrature.)--",
           R"--(Take data from scattering species of this index (0-based) in *scat_data*.)--",
           R"--(Maximum allowed extinction fraction in each of the edge size bins.)--",
           R"--(Minimum bulk extinction in the processed scattering species for which to apply size grid representation checks.)--",
           R"--(Minimum ratio of edge point pnd to maximum pnd of this scattering element over all pressure levels.)--"},

  };

  wsm_data["pndFromPsdBasic"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates *pnd_data* from given *psd_data*.

As *pndFromPsdBasic*, but without bulk extinction representation
checks.
)--",
      .author = {"Jana Mendrok, Patrick Eriksson"},
      .out = {"pnd_data", "dpnd_data_dx"},

      .in = {"pnd_size_grid", "psd_data", "psd_size_grid", "dpsd_data_dx"},
      .gin = {"quad_order"},
      .gin_type = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc = {R"--(Order of bin quadrature.)--"},

  };

  wsm_data["pnd_fieldCalcFromParticleBulkProps"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Converts particle bulk property data to *pnd_field*.

In short, the method combines *scat_species*, *pnd_agenda_array*,
``particle_bulkprop_field`` and their associated variables to derive
*pnd_field*.

The method does nothing if cloudbox is inactive.

Otherwise, cloudbox limits must be set before calling the method,
and ``particle_bulkprop_field`` is checked to have non-zero elements
just inside the cloudbox.
)--",
          .author = {"Patrick Eriksson, Jana Mendrok"},
          .out = {"pnd_field", "dpnd_field_dx"},

          .in = {"cloudbox_on",
                 "cloudbox_limits",
                 "scat_species",
                 "scat_data",
                 "scat_meta",
                 "atm_field",
                 "particle_bulkprop_names",
                 "pnd_agenda_array",
                 "pnd_agenda_array_input_names",
                 "jacobian_do",
                 "jacobian_quantities"},

          .pass_workspace = true,

      };

  wsm_data["pnd_fieldExpand1D"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Maps a 1D pnd_field to a (homogeneous) 2D or 3D pnd_field.

This method takes a 1D *pnd_field* and converts it to a 2D or 3D
"cloud". It is assumed that a complete 1D case has been created,
and after this ``atmosphere_dim``, ``lat_grid``, ``lon_grid`` and
*cloudbox_limits* have been changed to a 2D or 3D case (without
changing the vertical extent of the cloudbox.

No modification of *pnd_field* is made for the pressure dimension.
At the latitude and longitude cloudbox edge points *pnd_field* is set to
zero. This corresponds to nzero=1. If you want a larger margin between
the lat and lon cloudbox edges and the "cloud" you increase
``nzero``, where ``nzero`` is the number of grid points for which
*pnd_field* shall be set to 0, counted from each lat and lon edge.

See further ``AtmFieldsExpand1D``.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"pnd_field"},

      .in = {"pnd_field", "cloudbox_on", "cloudbox_limits"},
      .gin = {"nzero"},
      .gin_type = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc = {R"--(Number of zero values inside lat and lon limits.)--"},

  };

  wsm_data["pnd_fieldZero"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *pnd_field* to zero.

Creates an empty *pnd_field* of cloudbox size according to
*cloudbox_limits* and with number of scattering elemements
according to *scat_data*. If *scat_data* is not set yet, it will be
filled with one dummy scattering element.

The method works with both *scat_data* and *scat_data_raw*.
This method primarily exists for testing purposes.
On the one hand, empty *pnd_field* runs can be used to test the
agreement between true clear-sky (``cloudboxOff``) solutions and the
scattering solver solution in factual clear-sky conditions. It is
important to avoid discontinuities when switching from thin-cloud
to clear-sky conditions.
Moreover, scattering calculations using the DOIT method include
interpolation errors. If one is interested in this effect, one
should compare the DOIT result with an empty cloudbox to a clearsky
calculation. That means that the iterative method is performed for
a cloudbox with no particles.
)--",
      .author = {"Claudia Emde, Jana Mendrok"},
      .out = {"pnd_field", "dpnd_field_dx", "scat_data"},

      .in = {"scat_data", "f_grid", "cloudbox_limits", "jacobian_quantities"},

  };

  wsm_data["ppathAddGridCrossings"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds grids crossings to an existing *ppath*.

This method allows including certain altitudes, latitudes and
longitudes in *ppath*, if they are passed. These additional points
to consider should in general constitute one or several merged
atmospheric grids. By including the grid of e.g. a species, it is
ensured that min and max values get represented when interpolating
the atmosphere to the path.

You can set one or several of the grids to be non-empty.

The input *ppath* shall contain enough points that it can be treated
as a straight line between the points. For pure geometrical cases,
the end points of the path are sufficient, as long as altitudes,
latitudes and longitudes are stictly increasing or decreasing along
the path. If the later not is valid, such as when passing a tangent
point or one of the geographical poles, more input path points are
required (max 10 km distance?).
)--",
      .author = {"Patrick Eriksson"},
      .out = {"ppath"},

      .in = {"ppath", "ppath_lstep", "surface_field"},
      .gin = {"z_grid", "lat_grid", "lon_grid"},
      .gin_type = {"Vector", "Vector", "Vector"},
      .gin_value = {Vector{}, Vector{}, Vector{}},
      .gin_desc = {R"--(Grid/set of altitudes to include, if passed.)--",
                   R"--(Grid/set of latitudes to include, if passed.)--",
                   R"--(Grid/set of longitudes to include, if passed.)--"},

  };

  wsm_data["ppathCheckEndPoint"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Checks that a propagation path ends as expected.

Please note that ppaths are stored in observation direction and the "end"
is at the radiative background.

For example, to check the end altitude, set the GIN ``altitude`` to the
expected value and ``daltitude`` to the allowed tolerance. Latitude,
longitude, zenith angle and azimuth angle can be checked in the same way.

A check is done as soon the tolerance value is >= 0. Don't forget to set
the expected value, otherwise 0 will be applied. 

The radiative background and number of points can be checked in the same
way, but here there are no tolarance values and the expected values are
strings. The following coding is used for the radiative background:
  "Undefined"
  "Space"
  "Surface"
  "Cloudbox"
  "Transmitter"
  "StopDistance" - Start point determined by overall length criterion
)--",
      .author = {"Patrick Eriksson"},

      .in = {"ppath"},
      .gin = {"background",
              "np",
              "altitude",
              "daltitude",
              "latitude",
              "dlatitude",
              "longitude",
              "dlongitude",
              "zenith_angle",
              "dzenith_angle",
              "azimuth_angle",
              "dazimuth_angle"},
      .gin_type = {"String",
                   "Index",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric"},
      .gin_value = {String("Undefined"),
                    Index{-1},
                    Numeric{0},
                    Numeric{-1},
                    Numeric{0},
                    Numeric{-1},
                    Numeric{0},
                    Numeric{-1},
                    Numeric{0},
                    Numeric{-1},
                    Numeric{0},
                    Numeric{-1}},
      .gin_desc = {R"--(Expected radiative background. See above.)--",
                   R"--(Expected number of path points.)--",
                   R"--(Expected altitude.)--",
                   R"--(Allowed deviation for altitude.)--",
                   R"--(Expected latitude.)--",
                   R"--(Allowed deviation for latitude.)--",
                   R"--(Expected longitude.)--",
                   R"--(Allowed deviation for longitude.)--",
                   R"--(Expected zenith angle.)--",
                   R"--(Allowed deviation for zenith angle.)--",
                   R"--(Expected azimuth angle.)--",
                   R"--(Allowed deviation for azimuth angle.)--"},

  };

  wsm_data["ppathCheckInsideDomain"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Checks that propagation path is fully inside specified domain.

An error is issued if a point in *ppath* is outside of ranges
[lat_min, lat_max] and [lon_min, lon_max], in latitude and
longitude, respectively. The path is checked starting at the end
furthest away from the sensor and an error is given as soon an
incorrect point is found. This is not guaranteed to be the point
most outside of the domain.
)--",
      .author = {"Patrick Eriksson"},

      .in = {"ppath"},
      .gin = {"lat_min", "lat_max", "lon_min", "lon_max"},
      .gin_type = {"Numeric", "Numeric", "Numeric", "Numeric"},
      .gin_value = {Numeric{-90.0},
                    Numeric{90.0},
                    Numeric{-180.0},
                    Numeric{360.0}},
      .gin_desc = {R"--(Lowest allowed latitude.)--",
                   R"--(Highest allowed latitude.)--",
                   R"--(Lowest allowed longitude.)--",
                   R"--(Highest allowed longitude.)--"},

  };

  wsm_data["ppathCheckInsideGrids"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Checks that propagation path is fully inside grid ranges.

As *ppathCheckInsideDomain* but with the domain specified by a
combination of latitude and longitude grids.
)--",
      .author = {"Patrick Eriksson"},

      .in = {"ppath"},
      .gin = {"latitude_grid", "longitude_grid"},
      .gin_type = {"Vector", "Vector"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Latitude grid to not exceed.)--",
                   R"--(Longitude grid to not exceed.)--"},

  };

  wsm_data["ppathGeometric"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Geometric propagation path with fixed step length.

The propagation path (ppath) from *rte_pos* in the direction of
*rte_los* is determined. Refraction is ignored and the ppath is
denoted as geometrical. With default settings, the ppath ends
either at the surface or in space, depending on the observation
geometry. The points describing the ppath have an equidistant
spacing, that is the highest possible value satisfying *ppath_lstep*.

Possible intersections with the surface are determined following
*IntersectionGeometricSurface*.

It is possible to set a maximum length of the ppath (for the part
inside the atmosphere) by *ppath_ltotal*. When the length of the
ppath is governed by this variable, the end point of the ppath
is then a point inside the atmosphere.

With ``include_specular_ppath`` set to true, the propagation path is
continued at an intersection with the surface. The additional section
is calculated for the specular direction, with surface topography
considered. The surface intersection point will appear twice in
*ppath*. The case of multiple surface intersections is handled.
 
The *atm_field* is only used for its top of the atmosphere altitude
)--",
      .author = {"Patrick Eriksson"},
      .out = {"ppath"},

      .in = {"rte_pos",
             "rte_los",
             "ppath_lstep",
             "ppath_ltotal",
             "surface_field",
             "surface_search_accuracy",
             "surface_search_safe",
             "atm_field"},
      .gin = {"include_specular_ppath"},
      .gin_type = {"Index"},
      .gin_value = {Index{0}},
      .gin_desc = {R"--(Flag to continue path after surface intersection.)--"},

  };

  wsm_data["ppathRefracted"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculates propagation paths, including refraction by a basic approach.

Refraction is taken into account by probably the simplest approach
possible. The ray tracing is made by piece-wise geometric steps.
At the end of each step, the zenith and azimuth propagation angles
are updated following the local gradients of the refractive index.

Otherwise works largely as *ppathGeometric*. However, the spacing
between path points will here in general exactly be *ppath_lstep*.
The exception is the last ppath step, where the length is adjusted
to the remaining distance of the path. To be clear, the ray tracing
step length is governed by *ppath_lraytrace*.

Surface intersections are found in manner matching setting
*surface_search_safe* to 0 (see *IntersectionGeometricSurface*).
This for efficiency reasons, but also as the ray tracing largely
removes the need for a "safe" search.

For more accurate calculations, but slower, consider the two GIN
parameters ``do_horizontal_gradients`` and ``do_twosided_perturb``

Default is to only determine the altitude gradients of the refractive
index, as this is in general the only relevant term. To also calculate
and consider the latitude and longitude gradients, set
``do_horizontal_gradients`` to true.

The gradients of the refractive index are obtained by perturbing the
position of concern with small positive values. With ``do_twosided_perturb``
set to true, there is also a perturbation in the negative direction.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"ppath"},

      .in = {"refr_index_air_ZZZ_agenda",
             "rte_pos",
             "rte_los",
             "ppath_lstep",
             "ppath_ltotal",
             "ppath_lraytrace",
             "surface_field",
             "surface_search_accuracy"},
      .gin = {"z_toa",
              "do_horizontal_gradients",
              "do_twosided_perturb",
              "include_specular_ppath"},
      .gin_type = {"Numeric", "Index", "Index", "Index"},
      .gin_value = {std::nullopt, Index{0}, Index{0}, Index{0}},
      .gin_desc =
          {R"--(Top-of-the-atmosphere altitude.)--",
           R"--(Consider horisontal gradients of refractive index.)--",
           R"--(Perform double-sided perturbations when calculating refractive index gradients.)--",
           R"--(See *ppathGeometric*.)--"},
      .pass_workspace = true,

  };

  wsm_data["ppath_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *ppath_agenda* to a default value

Options are:

- ``"FollowSensorLosPath"``:

    1. Uses ``ppathStepByStep`` to set *ppath*

- ``"PlaneParallel"``:

    1. Uses ``ppathPlaneParallel`` to set *ppath*

- ``"TransmitterReceiverPath"``:

    1. Uses ``rte_losGeometricFromRtePosToRtePos2`` to set *rte_los*
    2. Uses ``ppathFromRtePos2`` to set *ppath*, and also to modify *rte_los*, and *ppath_lraytrace*
)--",
      .author = {"Richard Larsson"},
      .out = {"ppath_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["ppath_step_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *ppath_step_agenda* to a default value

Options are:

- ``"GeometricPath"``:

    1. Uses ``ppath_stepGeometric`` to modify *ppath*
- ``"RefractedPath"``:

    1. Uses ``ppath_stepRefractionBasic`` to modify *ppath*
)--",
      .author = {"Richard Larsson"},
      .out = {"ppath_step_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["ppvar_atmFromPath"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gets the atmospheric points along the path.
)--",
      .author = {"Richard Larsson"},
      .out = {"ppvar_atm"},

      .in = {"ppath", "atm_field"},

  };

  wsm_data["ppvar_cumtramatForward"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *ppvar_cumtramat* by forward iteration of *ppvar_tramat*
)--",
      .author = {"Richard Larsson"},
      .out = {"ppvar_cumtramat"},

      .in = {"ppvar_tramat"},

  };

  wsm_data["ppvar_cumtramatReverse"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *ppvar_cumtramat* by reverse iteration of *ppvar_tramat*
)--",
      .author = {"Richard Larsson"},
      .out = {"ppvar_cumtramat"},

      .in = {"ppvar_tramat"},

  };

  wsm_data["ppvar_fFromPath"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gets the frequency grid along the path.
)--",
      .author = {"Richard Larsson"},
      .out = {"ppvar_f"},

      .in = {"f_grid", "ppath", "ppvar_atm", "rte_alonglos_v"},

  };

  wsm_data["ppvar_optical_depthFromPpvar_trans_cumulat"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Sets *ppvar_optical_depth* according to provided transmittance data.

The values in ppvar_optical_depth are set to
-log( ppvar_trans_cumulat(joker,joker,0,0) ).
)--",
          .author = {"Patrick Eriksson"},
          .out = {"ppvar_optical_depth"},

          .in = {"ppvar_trans_cumulat"},

      };

  wsm_data["ppvar_propmatCalc"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Gets the propagation matrix and NLTE source term along the path.
)--",
      .author = {"Richard Larsson"},
      .out = {"ppvar_propmat", "ppvar_nlte", "ppvar_dpropmat", "ppvar_dnlte"},

      .in = {"propmat_clearsky_agenda",
             "jacobian_quantities",
             "ppvar_f",
             "ppath",
             "ppvar_atm",
             "jacobian_do"},

      .pass_workspace = true,

  };

  wsm_data["ppvar_radCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gets the radiation along the path.
)--",
      .author = {"Richard Larsson"},
      .out = {"ppvar_rad", "ppvar_drad"},

      .in = {"background_rad",
             "ppvar_src",
             "ppvar_dsrc",
             "ppvar_tramat",
             "ppvar_cumtramat",
             "ppvar_dtramat",
             "ppvar_propmat",
             "ppvar_dpropmat",
             "ppvar_distance",
             "ppvar_ddistance",
             "rt_integration_option"},

  };

  wsm_data["ppvar_radCalcEmission"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Gets the radiation along the path by linear emission calculations.
)--",
      .author = {"Richard Larsson"},
      .out = {"ppvar_rad", "ppvar_drad"},

      .in = {"background_rad",
             "ppvar_src",
             "ppvar_dsrc",
             "ppvar_tramat",
             "ppvar_cumtramat",
             "ppvar_dtramat"},

  };

  wsm_data["ppvar_radCalcTransmission"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Gets the radiation along the path by linear emission calculations.
)--",
      .author = {"Richard Larsson"},
      .out = {"ppvar_rad", "ppvar_drad"},

      .in = {"ppvar_tramat", "ppvar_cumtramat", "ppvar_dtramat"},

  };

  wsm_data["ppvar_rtprop_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *ppvar_rtprop_agenda* to a default value

Options are:
    FIXME
)--",
      .author = {"Richard Larsson"},
      .out = {"ppvar_rtprop_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["ppvar_srcFromPropmat"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gets the source term along the path.
)--",
      .author = {"Richard Larsson"},
      .out = {"ppvar_src", "ppvar_dsrc"},

      .in = {"ppvar_propmat",
             "ppvar_nlte",
             "ppvar_dpropmat",
             "ppvar_dnlte",
             "ppvar_f",
             "ppvar_atm",
             "jacobian_quantities",
             "jacobian_do"},

  };

  wsm_data["ppvar_tramatCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gets the transmission matrix in layers along the path.

A layer is defined as made up by the average of 2 levels, thus the outer-most size
of the derivatives out of this function is 2.
)--",
      .author = {"Richard Larsson"},
      .out = {"ppvar_tramat",
              "ppvar_dtramat",
              "ppvar_distance",
              "ppvar_ddistance"},

      .in = {"ppvar_propmat",
             "ppvar_dpropmat",
             "ppath",
             "ppvar_atm",
             "jacobian_quantities",
             "jacobian_do"},

  };

  wsm_data["predefined_model_dataAddWaterMTCKD400"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Sets the data for MT CKD 4.0 Water model

Note that the vectors must have the same length, and that wavenumbers must be growing
at a constant rate.  The minimum length is 4.

Note also that as this is predefined model data, the units of the values of the vectors
must be as described by each vector.
)--",
          .author = {"Richard Larsson"},
          .out = {"predefined_model_data"},

          .in = {"predefined_model_data"},
          .gin = {"ref_temp",
                  "ref_press",
                  "ref_h2o_vmr",
                  "self_absco_ref",
                  "for_absco_ref",
                  "wavenumbers",
                  "self_texp"},
          .gin_type = {"Numeric",
                       "Numeric",
                       "Numeric",
                       "Vector",
                       "Vector",
                       "Vector",
                       "Vector"},
          .gin_value = {std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::nullopt},
          .gin_desc = {R"--(Reference temperature)--",
                       R"--(Reference pressure)--",
                       R"--(Reference volume mixing ratio of water)--",
                       R"--(Self absorption [1/(cm-1 molecules/cm^2])--",
                       R"--(Foreign absorption [1/(cm-1 molecules/cm^2)])--",
                       R"--(Wavenumbers [cm-1])--",
                       R"--(Self temperature exponent [-])--"},

      };

  wsm_data["predefined_model_dataInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initialize the predefined model data
)--",
      .author = {"Richard Larsson"},
      .out = {"predefined_model_data"},

  };

  wsm_data["propmat_clearskyAddCIA"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculate absorption coefficients per tag group for HITRAN CIA continua.

This interpolates the cross sections from *abs_cia_data*.

The robust option is intended only for testing. Do not use for normal
runs, since subsequent functions will not be able to deal with NAN values.
)--",
      .author = {"Stefan Buehler, Oliver Lemke"},
      .out = {"propmat_clearsky", "dpropmat_clearsky_dx"},

      .in = {"propmat_clearsky",
             "dpropmat_clearsky_dx",
             "abs_species",
             "select_abs_species",
             "jacobian_quantities",
             "f_grid",
             "atm_point",
             "abs_cia_data"},
      .gin = {"T_extrapolfac", "ignore_errors"},
      .gin_type = {"Numeric", "Index"},
      .gin_value = {Numeric{0.5}, Index{0}},
      .gin_desc =
          {R"--(Temperature extrapolation factor (relative to grid spacing).)--",
           R"--(Set to 1 to suppress runtime errors (and return NAN values instead).)--"},

  };

  wsm_data["propmat_clearskyAddFaraday"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates absorption matrix describing Faraday rotation.

Faraday rotation is a change of polarization state of an
electromagnetic wave propagating through charged matter by
interaction with a magnetic field. Hence, this method requires
*abs_species* to contain 'free_electrons' and electron content field
(as part of ``vmr_field``) as well as magnetic field (``mag_u_field``,
``mag_v_field``, ``mag_w_field``) to be specified.

Faraday rotation affects Stokes parameters 2 and 3 (but not
intensity!). Therefore, this method requires stokes_dim>2.

Like all 'propmat_clearskyAdd*' methods, the method is additive,
i.e., does not overwrite the propagation matrix *propmat_clearsky*,
but adds further contributions.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"propmat_clearsky", "dpropmat_clearsky_dx"},

      .in = {"propmat_clearsky",
             "dpropmat_clearsky_dx",
             "f_grid",
             "abs_species",
             "select_abs_species",
             "jacobian_quantities",
             "atm_point",
             "rtp_los"},

  };

  wsm_data["propmat_clearskyAddFromLookup"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Extract gas absorption coefficients from lookup table.

This extracts the absorption coefficient for all species from the
lookup table, and adds them to the propagation matrix. Extraction is
for one specific atmospheric condition, i.e., a set of pressure,
temperature, and VMR values.

Some special species are ignored, for example Zeeman species and free
electrons, since their absorption properties are not simple scalars
and cannot be handled by the lookup table.

The interpolation order in T and H2O is given by *abs_t_interp_order*
and *abs_nls_interp_order*, respectively.

Extraction is done for the frequencies in f_grid. Frequency
interpolation is controlled by *abs_f_interp_order*. If this is zero,
then f_grid must either be the same as the internal frequency grid of
the lookup table (for efficiency reasons, only the first and last
element of f_grid are checked), or must have only a single element.
If *abs_f_interp_order* is above zero, then frequency is interpolated
along with the other interpolation dimensions. This is useful for
calculations with Doppler shift.

For Doppler calculations, you should generate the table with a
somewhat larger frequency grid than the calculation itself has, since
the Doppler shift will push the frequency grid out of the table range
on one side.

Some extrapolation is allowed. For pressure and frequency interpolation
the standard extrapolation factor of 0.5 is applied. The factor is the
default for temperature and VMR interpolation, but the extrapolation
limit can here be adjusted by the ``extpolfac`` argument.
)--",
      .author = {"Stefan Buehler, Richard Larsson"},
      .out = {"propmat_clearsky", "dpropmat_clearsky_dx"},

      .in = {"propmat_clearsky",
             "dpropmat_clearsky_dx",
             "abs_lookup",
             "abs_lookup_is_adapted",
             "abs_p_interp_order",
             "abs_t_interp_order",
             "abs_nls_interp_order",
             "abs_f_interp_order",
             "f_grid",
             "atm_point",
             "jacobian_quantities",
             "abs_species",
             "select_abs_species"},
      .gin = {"extpolfac", "no_negatives"},
      .gin_type = {"Numeric", "Index"},
      .gin_value = {Numeric{0.5}, Index{1}},
      .gin_desc =
          {R"--(Extrapolation factor (for temperature and VMR grid edges).)--",
           R"--(Boolean. If it is true negative values due to interpolation are set to zero.)--"},

  };

  wsm_data["propmat_clearskyAddHitranLineMixingLines"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Calculates gas absorption coefficients line-by-line for HITRAN line mixed data.

*Wigner6Init* or *Wigner3Init* must be called before this function.

Note that you need to have *propmat_clearskyAddLines* in addition to this method
to compensate the calculations for the pressure limit

Please ensure you cite the original authors when you use this function:
	J. Lamouroux, L. Realia, X. Thomas, et al., J.Q.S.R.T. 151 (2015), 88-96
)--",
          .author = {"Richard Larsson"},
          .out = {"propmat_clearsky"},

          .in = {"propmat_clearsky",
                 "abs_hitran_relmat_data",
                 "abs_lines_per_species",
                 "isotopologue_ratios",
                 "f_grid",
                 "abs_species",
                 "select_abs_species",
                 "jacobian_quantities",
                 "atm_point"},

      };

  wsm_data["propmat_clearskyAddLines"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Computes the line-by-line unpolarized absorption and adds
it to the diagonal of *propmat_clearsky* and derivates to other variables.
Does the same for NLTE variables if required.

If ``lines_speedup_option`` is not "None", then some speed-up logic is applied.
Valid speed-up logic other than "None" includes:

- ``"LinearIndependent"``:
  Using a sparse-grid, the points are separated as [f0, f0+df[0], f0+df[0], f0+df[1]...]
  until the entire *f_grid* is covered.  All sparse bins are on *f_grid* so df changes.
  A linear interpolation scheme is used between the bins to fill up the dense
  absorption.  The maximum of df[n] is given by ``lines_sparse_df`` and the minimum
  transition between dense-to-sparse grid calculations are given by ``lines_sparse_lim``.

- ``"QuadraticIndependent"``:
  Using a sparse-grid, the points are separated as [f0, f0+0.5*df[0], f0+df[0], f0+df[0], f0+0.5*df[1], f0+df[1]...]
  until the entire *f_grid* is covered.  All sparse bins are on *f_grid* so df changes.
  A quadratic interpolation scheme is used between the bins to fill up the dense
  absorption.  The maximum of df[n] is given by ``lines_sparse_df`` and the minimum
  transition between dense-to-sparse grid calculations are given by ``lines_sparse_lim``.

Please use *sparse_f_gridFromFrequencyGrid* to see the sparse frequency grid

By default we discourage negative values, which are common when using one of the line mixing
approximations.   Change the value of no_negatives to 0 to allow these negative absorptions.
)--",
      .author = {"Richard Larsson"},
      .out = {"propmat_clearsky",
              "nlte_source",
              "dpropmat_clearsky_dx",
              "dnlte_source_dx"},

      .in = {"propmat_clearsky",
             "nlte_source",
             "dpropmat_clearsky_dx",
             "dnlte_source_dx",
             "f_grid",
             "abs_species",
             "select_abs_species",
             "jacobian_quantities",
             "abs_lines_per_species",
             "isotopologue_ratios",
             "atm_point",
             "nlte_vib_energies",
             "nlte_do",
             "lbl_checked"},
      .gin = {"lines_sparse_df",
              "lines_sparse_lim",
              "lines_speedup_option",
              "no_negatives"},
      .gin_type = {"Numeric", "Numeric", "String", "Index"},
      .gin_value = {Numeric{0}, Numeric{0}, String("None"), Index{1}},
      .gin_desc =
          {R"--(The grid sparse separation)--",
           R"--(The dense-to-sparse limit)--",
           R"--(Speedup logic)--",
           R"--(Boolean.  If it is true, line mixed bands each allocate their own compute data to ensure that they cannot produce negative absorption)--"},

  };

  wsm_data["propmat_clearskyAddOnTheFlyLineMixing"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Compute the line mixing of matching lines and add it to the propagation matrix

Each band's Population Type is checked and the calculations are only performed
for bands with matching population types (and a good pressure limits)

Presently only supports one method: ByMakarovFullRelmat, based on Makarov et al 2020

*Wigner6Init* or *Wigner3Init* must be called before this function.

Note that you need to have *propmat_clearskyAddLines* addition to this method
to compensate the calculations for the pressure limit
)--",
          .author = {"Richard Larsson"},
          .out = {"propmat_clearsky", "dpropmat_clearsky_dx"},

          .in = {"propmat_clearsky",
                 "dpropmat_clearsky_dx",
                 "abs_lines_per_species",
                 "ecs_data",
                 "isotopologue_ratios",
                 "f_grid",
                 "abs_species",
                 "select_abs_species",
                 "jacobian_quantities",
                 "atm_point",
                 "lbl_checked"},

      };

  wsm_data["propmat_clearskyAddOnTheFlyLineMixingWithZeeman"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Compute the line mixing of matching lines and add it to the propagation matrix
Also computes Zeeman effect for all the lines in the band

Each band's Population Type is checked and the calculations are only performed
for bands with matching population types (and a good pressure limits)

Presently only supports one method: ByMakarovFullRelmat, based on Makarov et al 2020

*Wigner6Init* or *Wigner3Init* must be called before this function.

Note that you need to have *propmat_clearskyAddLines* in addition to this method
to compensate the calculations for the pressure limit
)--",
          .author = {"Richard Larsson"},
          .out = {"propmat_clearsky", "dpropmat_clearsky_dx"},

          .in = {"propmat_clearsky",
                 "dpropmat_clearsky_dx",
                 "abs_lines_per_species",
                 "ecs_data",
                 "isotopologue_ratios",
                 "f_grid",
                 "abs_species",
                 "select_abs_species",
                 "jacobian_quantities",
                 "atm_point",
                 "rtp_los",
                 "lbl_checked"},

      };

  wsm_data["propmat_clearskyAddParticles"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates absorption coefficients of particles to be used in
clearsky (non-cloudbox) calculations.

This is a method to include particles (neglecting possible
scattering components) in a clearsky calculation, i.e. without
applying the cloudbox and scattering solvers. Particles are handled
as absorbing species with one instance of 'particles' per scattering
element considered added to *abs_species*. Particle absorption cross-
sections at current atmospheric conditions are extracted from the
single scattering data stored in *scat_data*, i.e., one array
element per 'particles' instance in *abs_species* is required. Number
densities are stored in ``vmr_field_raw`` or ``vmr_field`` as for all
*abs_species*, but can be taken from (raw) pnd_field type data.

Note that the absorption coefficient is applied both in the
extinction term (neglecting scattering out of the line of sight)
and the emission term (neglecting the scattering source term, i.e.
scattering into the line of sight).

Optionally, particle extinction (sum of absorption and scattering
coefficient) can be used instead of absorption only. To choose this
case, set the ``use_abs_as_ext`` flag to 0. However, be aware that
this creates some unphysical emission term, hence is only suitable,
where the source term is negligible anyways, e.g. for occultation
simulations.

A line-of-sight direction *rtp_los* is required as particles can
exhibit directional dependent absorption properties, which is taken
into account by this method.
*ScatElementsToabs_speciesAdd* can be used to add all required
settings/data for individual scattering elements at once, i.e. a
'particles' tag to *abs_species*, a set of single scattering data to
*scat_data* and a number density field to ``vmr_field_raw``
(``vmr_field`` is derived applying AtmFieldsCalc once VMRs for all
*abs_species* have been added) is appended for each scattering
element.

Like all 'propmat_clearskyAdd*' methods, the method is additive,
i.e., does not overwrite the propagation matrix *propmat_clearsky*,
but adds further contributions.
)--",
      .author = {"Jana Mendrok"},
      .out = {"propmat_clearsky", "dpropmat_clearsky_dx"},

      .in = {"propmat_clearsky",
             "dpropmat_clearsky_dx",
             "f_grid",
             "abs_species",
             "select_abs_species",
             "jacobian_quantities",
             "rtp_los",
             "atm_point",
             "scat_data",
             "scat_data_checked"},
      .gin = {"use_abs_as_ext"},
      .gin_type = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc =
          {R"--(A flag with value 1 or 0. If set to one, particle absorption is used in extinction and emission parts of the RT equation, and scattering out of LOS as well as into LOS is neglected. Otherwise, particle extinction (absorption+scattering) is applied in both the extinction as well as the emission part of the RT equation. That is, true extinction is applied, but emission also includes a pseudo-emission contribution from the scattering coefficient. )--"},

  };

  wsm_data["propmat_clearskyAddPredefined"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Adds all of the predefined models in *abs_species* to the propmat_clearsky

Only supports temperature and wind speed derivatives

Available models:

- O2-MPM2020:
  60 GHz and 118 GHz lines only (no continua, no higher Hz line centers) 

  Dmitriy S. Makarov, Mikhail Yu. Tretyakov, Philip W. Rosenkranz, JQSRT 243, 2020, Revision of the 
  60-GHz atmospheric oxygen absorption band models for practical use, 
  https://doi.org/10.1016/j.jqsrt.2019.106798

- H2O-ForeignContCKDMT350:
  Foreign continua.  Expects H2O line center cutoff at 25 cm-1

  CKD_MTv3.50 H2O foreign continuum from the FORTRAN77 code written by Atmospheric and Environmental Research Inc. (AER),
  Radiation and Climate Group 131 Hartwell Avenue Lexington, MA 02421, USA
  http://www.rtweb.aer.com/continuum_frame.html

- H2O-SelfContCKDMT350:
  Self continua.  Expects H2O line center cutoff at 25 cm-1

  CKD_MTv3.50 H2O self continuum from the FORTRAN77 code written by Atmospheric and Environmental Research Inc. (AER),
  Radiation and Climate Group 131 Hartwell Avenue Lexington, MA 02421, USA
  http://www.rtweb.aer.com/continuum_frame.html

- H2O-ForeignContCKDMT320:
  Foreign continua.  Expects H2O line center cutoff at 25 cm-1

  CKD_MTv3.20 H2O foreign continuum from the FORTRAN77 code written by Atmospheric and Environmental Research Inc. (AER),
  Radiation and Climate Group 131 Hartwell Avenue Lexington, MA 02421, USA
  http://www.rtweb.aer.com/continuum_frame.html

- H2O-SelfContCKDMT320:
  Self continua.  Expects H2O line center cutoff at 25 cm-1

  CKD_MTv3.20 H2O self continuum from the FORTRAN77 code written by Atmospheric and Environmental Research Inc. (AER),
  Radiation and Climate Group 131 Hartwell Avenue Lexington, MA 02421, USA
  http://www.rtweb.aer.com/continuum_frame.html

- H2O-SelfContCKDMT400:
  Self continuum for water.  General reference: Mlawer et al. (2012), doi:10.1098/rsta.2011.0295

  Our code is reimplemented based on original Fortran90 code that is/was/will-be-made available via hitran.org

  Note that this model comes with the copyright statement [1].

  Note also that this model requires *predefined_model_data* to contain relevant data set either using
  *predefined_model_dataAddWaterMTCKD400* or via some file reading routine.

- H2O-ForeignContCKDMT400:
  Foreign continuum for water.  General reference: Mlawer et al. (2012), doi:10.1098/rsta.2011.0295

  Our code is reimplemented based on original Fortran90 code that is/was/will-be-made available via hitran.org

  Note that this model comes with the copyright statement [1].

  Note also that this model requires *predefined_model_data* to contain relevant data set either using
  *predefined_model_dataAddWaterMTCKD400* or via some file reading routine.

- H2O-ForeignContStandardType:
  Water microwave continua

  P. W. Rosenkranz., Radio Science, 33(4), 919, 1998 and
  Radio Science, Vol. 34(4), 1025, 1999.

- H2O-SelfContStandardType:
  Water microwave continua

  P. W. Rosenkranz., Radio Science, 33(4), 919, 1998 and
  Radio Science, Vol. 34(4), 1025, 1999.

- H2O-MPM89:
  Microwave water absorption model

  H. J. Liebe, Int. J. Infrared and Millimeter Waves, 10(6), 1989, 631.

- H2O-PWR98:
  Microwave water absorption model

  P. W. Rosenkranz., Radio Science, 33(4), 919, 1998 and
  Radio Science, Vol. 34(4), 1025, 1999.

- H2O-PWR2021:
  Microwave water absorption model developed by P.W. Rosenkranz.

  Our code is reimplemented based on the Fortran code available at http://cetemps.aquila.infn.it/mwrnet/lblmrt_ns.html

- H2O-PWR2022:
  Microwave water absorption model developed by P.W. Rosenkranz.

  Our code is reimplemented based on the Fortran code available at http://cetemps.aquila.infn.it/mwrnet/lblmrt_ns.html

- CO2-CKDMT252:
  MT CKD absorption for CO2

  This absorption model is taken from the FORTRAN77 code of
  CKD_MT version 2.50 written by<br>
  Atmospheric and Environmental Research Inc. (AER),<br>
  Radiation and Climate Group<br>
  131 Hartwell Avenue<br>
  Lexington, MA 02421, USA<br>
  http://www.rtweb.aer.com/continuum_frame.html

- O2-CIAfunCKDMT100:
  CIA for oxygen from MT CKD

  F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman,
  J.-M. Hartmann, Ch. Boulet,<br>
  Infrared collision-induced absorption by O2 near 6.4 microns for
  atmospheric applications: measurements and emprirical modeling,<br>
  Appl. Optics, 35, 5911-5917, (1996).

  This absorption model is taken from the FORTRAN77 code of
  CKD_MT version 1.00 written by<br>
  Atmospheric and Environmental Research Inc. (AER),<br>
  Radiation and Climate Group<br>
  131 Hartwell Avenue<br>
  Lexington, MA 02421, USA<br>
  http://www.rtweb.aer.com/continuum_frame.html

- O2-MPM89:
  Oxygen microwave absorption model

  Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
  <i>Propagation modeling of moist air and suspended water/ice
  particles at frequencies below 1000 GHz</i>,<br>
  AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

- O2-PWR98:
  Oxygen microwave absorption model

  P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
  BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993).
  H.J. Liebe et al, JQSRT V.48, PP.629-643 (1992).
  M.J. Schwartz, Ph.D. thesis, M.I.T. (1997).
  SUBMILLIMETER LINE INTENSITIES FROM HITRAN96.

- O2-PWR2021:
  Oxygen microwave absorption model developed by P.W. Rosenkranz.

  Our code is reimplemented based on the Fortran code available at http://cetemps.aquila.infn.it/mwrnet/lblmrt_ns.html

- O2-PWR2022:
  Oxygen microwave absorption model developed by P.W. Rosenkranz.

  Our code is reimplemented based on the Fortran code available at http://cetemps.aquila.infn.it/mwrnet/lblmrt_ns.html

- O2-SelfContStandardType:
  Microwave continua term

  Reference: P. W. Rosenkranz, Chapter 2, in M. A. Janssen, <br>
  <I>Atmospheric Remote Sensing by Microwave Radiometry</i>,<br>
  John Wiley & Sons, Inc., 1993.<br>
  <br>
  Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
  <i>Propagation modeling of moist air and suspended water/ice
  particles at frequencies below 1000 GHz</i>,<br>
  AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

- O2-TRE05:
  Oxygen microwave absorption model

  References: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
  <i>Propagation modeling of moist air and suspended water/ice
  particles at frequencies below 1000 GHz</i>,<br>
  AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

  M.Yu. Tretyakov, M.A. Koshelev, V.V. Dorovskikh,
  D.S. Makarov, P.W. Rosenkranz; 60-GHz oxygen band: precise broadening and central frequencies
  of fine-structure lines, absolute absorption profile
  at atmospheric pressure, and revision of mixing coefficients
  doi:10.1016/j.jms.2004.11.011

- O2-v0v0CKDMT100:
  MT CKD

  CKD_MT 1.00 implementation of oxygen collision induced fundamental model of
  O2 continuum formulated by
  Mate et al. over the spectral region 7550-8486 cm-1:
  B. Mate, C. Lugez, G.T. Fraser, W.J. Lafferty,
  "Absolute Intensities for the O2 1.27 micron
  continuum absorption",
  J. Geophys. Res., 104, 30,585-30,590, 1999.

  Also, refer to the paper "Observed  Atmospheric
  Collision Induced Absorption in Near Infrared Oxygen Bands",
  Mlawer, Clough, Brown, Stephen, Landry, Goldman, & Murcray,
  Journal of Geophysical Research (1997).

  This absorption model is taken from the FORTRAN77 code of
  CKD_MT version 1.00 written by<br>
  Atmospheric and Environmental Research Inc. (AER),<br>
  Radiation and Climate Group<br>
  131 Hartwell Avenue<br>
  Lexington, MA 02421, USA<br>
  http://www.rtweb.aer.com/continuum_frame.html<br>
  <br>

- O2-v1v0CKDMT100:
  MT CKD

  Mlawer, Clough, Brown, Stephen, Landry, Goldman, Murcray,<br>
  Observed  Atmospheric Collision Induced Absorption in Near Infrared Oxygen Bands,<br>
  J. Geophys. Res., 103, D4, 3859-3863, 1998.

  This absorption model is taken from the FORTRAN77 code of
  CKD_MT version 1.00 written by<br>
  Atmospheric and Environmental Research Inc. (AER),<br>
  Radiation and Climate Group<br>
  131 Hartwell Avenue<br>
  Lexington, MA 02421, USA<br>
  http://www.rtweb.aer.com/continuum_frame.html<br>

- O2-visCKDMT252:
  MT CKD

  O2 continuum formulated by Greenblatt et al. over the spectral region
  8797-29870 cm-1:  "Absorption Coefficients of Oxygen Between 
  330 and 1140 nm, G.D. Green blatt, J.J. Orlando, J.B. Burkholder,
  and A.R. Ravishabkara,  J. Geophys. Res., 95, 18577-18582, 1990.

  This absorption model is taken from the FORTRAN77 code of
  CKD_MT version 2.50 written by<br>
  Atmospheric and Environmental Research Inc. (AER),<br>
  Radiation and Climate Group<br>
  131 Hartwell Avenue<br>
  Lexington, MA 02421, USA<br>
  http://www.rtweb.aer.com/continuum_frame.html<br>

- N2-CIAfunCKDMT252:
  MT CKD

  Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and
  J._M. Hartmann,<br>
  Infrared collision-induced absorption by
  N2 near 4.3 microns for atmospheric applications:
  Measurements and emprirical modeling, <br>
  Appl. Optics, 35, 5911-5917, (1996)

  This absorption model is taken from the FORTRAN77 code of
  CKD_MT version 1.00 written by<br>
  Atmospheric and Environmental Research Inc. (AER),<br>
  Radiation and Climate Group<br>
  131 Hartwell Avenue<br>
  Lexington, MA 02421, USA<br>
  http://www.rtweb.aer.com/continuum_frame.html

- N2-CIArotCKDMT252:
  MT CKD

  Borysow, A, and L. Frommhold,<br>
  Collision-induced rototranslational absorption spectra of N2-N2
  pairs for temperatures from 50 to 300 K,<br>
  The Astrophysical Journal, 311, 1043-1057, 1986.

  This absorption model is taken from the FORTRAN77 code of
  CKD_MT version 1.00 written by<br>
  Atmospheric and Environmental Research Inc. (AER),<br>
  Radiation and Climate Group<br>
  131 Hartwell Avenue<br>
  Lexington, MA 02421, USA<br>
  http://www.rtweb.aer.com/continuum_frame.html

- N2-SelfContStandardType:
  Microwave nitrogen absorption continua

  Reference: P. W. Rosenkranz, Chapter 2, in M. A. Janssen, <br>
  <I>Atmospheric Remote Sensing by Microwave Radiometry</i>,<br>
  John Wiley & Sons, Inc., 1993.

- N2-SelfContMPM93:
  Microwave nitrogen absorption continua from MPM93 model

  Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
  <i>Propagation modeling of moist air and suspended water/ice
  particles at frequencies below 1000 GHz</i>,<br>
  AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

- N2-SelfContPWR2021:
  Microwave nitrogen absorption continua developed by P.W. Rosenkranz.

  Note that this also includes O2-N2 and O2-O2 collision-induced absorption and is
  only applicable to Earth

  Our code is reimplemented based on the Fortran code available at http://cetemps.aquila.infn.it/mwrnet/lblmrt_ns.html

- liquidcloud-ELL07:
  Water droplet absorption

  W. J. Ellison, <br>
  <i>Permittivity of Pure Water, at Standard Atmospheric Pressure, over the
  Frequency Range 0-25 THz and Temperature Range 0-100C</i>,<br>
  J. Phys. Chem. Ref. Data, Vol. 36, No. 1, 2007
)--",
      .author = {"Richard Larsson"},
      .out = {"propmat_clearsky", "dpropmat_clearsky_dx"},

      .in = {"propmat_clearsky",
             "dpropmat_clearsky_dx",
             "predefined_model_data",
             "abs_species",
             "select_abs_species",
             "jacobian_quantities",
             "f_grid",
             "atm_point"},

  };

  wsm_data["propmat_clearskyAddScaledSpecies"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Adds a scaled target species absorption to *propmat_clearsky* and *nlte_source*

This recomputes the entire propagation matrix.  There are more efficient ways
to do these calculations but this method exist because of the composability it
offers
)--",
      .author = {"Richard Larsson"},
      .out = {"propmat_clearsky", "nlte_source"},

      .in = {"propmat_clearsky",
             "nlte_source",
             "jacobian_quantities",
             "select_abs_species",
             "f_grid",
             "rtp_los",
             "atm_point",
             "propmat_clearsky_agenda"},
      .gin = {"target", "scale"},
      .gin_type = {"ArrayOfSpeciesTag", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(Target species tags to rescale (must be in *abs_species*)--",
           R"--(Rescaling factor (e.g., 0.1 adds 10% of the species to the absorption))--"},
      .pass_workspace = true,

  };

  wsm_data["propmat_clearskyAddXsecFit"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculate absorption cross sections per tag group for HITRAN xsec species.

This broadens the cross section data from *xsec_fit_data* and
interpolates it onto the current f_grid.

Model data needs to be read in with *ReadXsecData* before calling
this method.
)--",
      .author = {"Oliver Lemke"},
      .out = {"propmat_clearsky", "dpropmat_clearsky_dx"},

      .in = {"propmat_clearsky",
             "dpropmat_clearsky_dx",
             "abs_species",
             "select_abs_species",
             "jacobian_quantities",
             "f_grid",
             "atm_point",
             "xsec_fit_data"},
      .gin = {"force_p", "force_t"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {Numeric{-1}, Numeric{-1}},
      .gin_desc = {R"--(Positive value forces constant pressure [Pa].)--",
                   R"--(Positive value forces constant temperature [K].)--"},

  };

  wsm_data["propmat_clearskyAddZeeman"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculates Zeeman-affected polarized propagation matrix and its
derivatives.

Otherwise as *propmat_clearskyAddFromLookup*
)--",
      .author = {"Richard Larsson"},
      .out = {"propmat_clearsky",
              "nlte_source",
              "dpropmat_clearsky_dx",
              "dnlte_source_dx"},

      .in = {"propmat_clearsky",
             "nlte_source",
             "dpropmat_clearsky_dx",
             "dnlte_source_dx",
             "abs_lines_per_species",
             "f_grid",
             "abs_species",
             "select_abs_species",
             "jacobian_quantities",
             "isotopologue_ratios",
             "atm_point",
             "nlte_vib_energies",
             "rtp_los",
             "nlte_do",
             "lbl_checked"},
      .gin = {"manual_mag_field", "H", "theta", "eta"},
      .gin_type = {"Index", "Numeric", "Numeric", "Numeric"},
      .gin_value = {Index{0}, Numeric{1.0}, Numeric{0.0}, Numeric{0.0}},
      .gin_desc = {R"--(Manual angles tag)--",
                   R"--(Manual Magnetic Field Strength)--",
                   R"--(Manual theta given positive tag)--",
                   R"--(Manual eta given positive tag)--"},

  };

  wsm_data["propmat_clearskyForceNegativeToZero"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Sets *propmat_clearsky* to match zero attenuation
if negative value.  Useful for line mixing in some cases.

Use this method just if you know what you are doing!
)--",
          .author = {"Richard Larsson"},
          .out = {"propmat_clearsky"},

          .in = {"propmat_clearsky"},

      };

  wsm_data["propmat_clearskyInit"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Initialize *propmat_clearsky*, *nlte_source*, and their derivatives to zeroes.

This method must be used inside *propmat_clearsky_agenda* and then be called first.
)--",
      .author = {"Oliver Lemke", "Richard Larsson"},
      .out = {"propmat_clearsky",
              "nlte_source",
              "dpropmat_clearsky_dx",
              "dnlte_source_dx"},

      .in = {"jacobian_quantities",
             "f_grid",
             "propmat_clearsky_agenda_checked"},

  };

  wsm_data["propmat_clearskyZero"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *propmat_clearsky* to match zero attenuation.

Use this method just if you know what you are doing!

If you want to make a calculation with no clear-sky attenuation at
all, fill *propmat_clearsky_agenda* with this method and required
Ignore statements (don't include *propmat_clearskyInit*).
)--",
      .author = {"Patrick Eriksson"},
      .out = {"propmat_clearsky"},

      .in = {"f_grid"},

  };

  wsm_data["propmat_clearsky_agendaAuto"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the *propmat_clearsky_agenda* automatically

This method introspects the input and uses it for generating the
*propmat_clearsky_agenda* automatically.  If ``use_abs_lookup``, all
methods that can be used to generate the absorption lookup table
are ignored and instead the calculations from the absorption
lookup are used.

The following methods are considered for addition:
    1) *propmat_clearskyInit*
    2) *propmat_clearskyAddCIA*
    3) *propmat_clearskyAddLines*
    4) *propmat_clearskyAddZeeman*
    5) *propmat_clearskyAddFaraday*
    6) *propmat_clearskyAddXsecFit*
    7) *propmat_clearskyAddParticles*
    8) *propmat_clearskyAddFromLookup*
    9) *propmat_clearskyAddPredefined*
    10) *propmat_clearskyAddOnTheFlyLineMixing*
    11) *propmat_clearskyAddHitranLineMixingLines*
    12) *propmat_clearskyAddOnTheFlyLineMixingWithZeeman*

To perform absorption lookupo table calculation, call:
    1) *propmat_clearsky_agendaAuto*
    2) ``abs_lookupCalc``  FIXME: HOW TO COMPUTE IT
    3) *propmat_clearsky_agendaAuto* (use_abs_lookup=1)
    4) Perform other calculations
)--",
      .author = {"Richard Larsson"},
      .out = {"propmat_clearsky_agenda", "propmat_clearsky_agenda_checked"},
      .in = {"abs_species", "abs_lines_per_species"},
      .gin = {"H",
              "T_extrapolfac",
              "eta",
              "extpolfac",
              "force_p",
              "force_t",
              "ignore_errors",
              "lines_sparse_df",
              "lines_sparse_lim",
              "lines_speedup_option",
              "manual_mag_field",
              "no_negatives",
              "theta",
              "use_abs_as_ext",
              "use_abs_lookup"},
      .gin_type = {"Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Index",
                   "Numeric",
                   "Numeric",
                   "String",
                   "Index",
                   "Index",
                   "Numeric",
                   "Index",
                   "Index"},
      .gin_value = {Numeric{1.0},
                    Numeric{0.5},
                    Numeric{0.0},
                    Numeric{0.5},
                    Numeric{-1},
                    Numeric{-1},
                    Index{0},
                    Numeric{0},
                    Numeric{0},
                    String("None"),
                    Index{0},
                    Index{1},
                    Numeric{0.0},
                    Index{1},
                    Index{0}},
      .gin_desc =
          {R"--(See *propmat_clearskyAddZeeman*)--",
           R"--(See *propmat_clearskyAddCIA*)--",
           R"--(See *propmat_clearskyAddZeeman*)--",
           R"--(See *propmat_clearskyAddFromLookup*)--",
           R"--(See *propmat_clearskyAddXsecFit*)--",
           R"--(See *propmat_clearskyAddXsecFit*)--",
           R"--(See *propmat_clearskyAddCIA*)--",
           R"--(See *propmat_clearskyAddLines*)--",
           R"--(See *propmat_clearskyAddLines*)--",
           R"--(See *propmat_clearskyAddLines*)--",
           R"--(See *propmat_clearskyAddZeeman*)--",
           R"--(See *propmat_clearskyAddLines*; See *propmat_clearskyAddFromLookup*)--",
           R"--(See *propmat_clearskyAddZeeman*)--",
           R"--(See *propmat_clearskyAddParticles*)--",
           R"--(Uses lookup calculations if true, ignores methods that can be part of the lookup table)--"},
  };

  wsm_data["propmat_clearsky_agendaGUI"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Opens a GUI for running the propagation matrix agenda

Note that this is not thread-safe and should be executed on the main workspace

The values of all non-control flow are automatically loaded from the workspace
if they are defined.  Otherwise some values are just selected
)--",
      .author = {"Richard Larsson"},

      .in = {"propmat_clearsky_agenda", "abs_species"},
      .gin = {"load"},
      .gin_type = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc = {R"--(Load non-logical variables from workspace if true)--"},
      .pass_workspace = true,

  };

  wsm_data["propmat_clearsky_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *propmat_clearsky_agenda* to a default value

Please consider using *propmat_clearsky_agendaAuto* instead of one of these options
as it will ensure you have the best coverage of use cases.  The options below are
available for feature testing

Options are:

- ``"Empty"``:

    1. Uses *propmat_clearskyInit* to set *propmat_clearsky*, *nlte_source*, *dpropmat_clearsky_dx*, and *dnlte_source_dx*
)--",
      .author = {"Richard Larsson"},
      .out = {"propmat_clearsky_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["propmat_clearsky_agenda_checkedCalc"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Checks if the *propmat_clearsky_agenda* contains all necessary
methods to calculate all the species in *abs_species*.

This method should be called if you use a manual *propmat_clearsky_agenda*
)--",
          .author = {"Oliver Lemke"},
          .out = {"propmat_clearsky_agenda_checked"},
          .in = {"abs_species", "propmat_clearsky_agenda"},
          .pass_workspace = true,
      };

  wsm_data["psdAbelBoutle12"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Abel and Boutle [2012] particle size distribution for rain.

Reference: Abel and Boutle, An improved representation of the 
raindrop size distribution for single-moment microphysics schemes,
QJRMS, 2012.

This is a 1-parameter PSD, i.e. *pnd_agenda_input* shall have one
column and *pnd_agenda_input_names* shall contain a single string.
The input data in *pnd_agenda_input* shall be rain mass content in
unit of [kg/m3]. The naming used is *pnd_agenda_input_names* is free
but the same name must be used in *particle_bulkprop_names* and
*dpnd_data_dx_names*.

Particles are assumed to be near-spherical, ie. *psd_size_grid* can
either be in terms of volume (or mass) equivalent diameter or
maximum diameter.

Derivatives are obtained analytically.

The validity range of mass content is not limited. Negative mass
contents will produce negative psd values following a distribution
given by abs(RWC), ie. abs(psd)=f(abs(RWC)).

If temperature is outside [``t_min``,``t_max``] psd=0 and dpsd=0 if
picky=0, or an error is thrown if picky=1.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names",
             "scat_species_a",
             "scat_species_b"},
      .gin = {"t_min", "t_max", "picky"},
      .gin_type = {"Numeric", "Numeric", "Index"},
      .gin_value = {Numeric{273}, Numeric{373}, Index{0}},
      .gin_desc =
          {R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdDelanoeEtAl14"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Normalized PSD as proposed in Delanoë et al. ((2014)),

Title and journal:
'Normalized particle size distribution for remote sensing
application', J. Geophys. Res. Atmos., 119, 4204–422.

The PSD has two independent parameters ``n0Star``, the intercept
parameter, and ``Dm``, the volume-weighted diameter.
This implementation expects as input two out of the following
three quantities: ``iwc``, ``n0Star``, ``Dm``. In this case one of
the input parameters ``iwc``, ``n0Star``, ``Dm`` must be set to -999.
It is also possible to provide only ``iwc``, in which case an a
priori assumption will be used to deduce ``n0Star`` from temperature.
In this case both ``n0Star`` and ``Dm`` must be set to -999.0.

This PSD is not defined for vanishing concentrations of
scatterers as it requires normalization by ``Dm``. It is up
to the user to ensure that the value of ``Dm`` is sufficiently
large. An error is thrown if ``Dm`` is zero or below the value
provided by ``dm_min``.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names"},
      .gin = {"iwc",
              "n0Star",
              "Dm",
              "rho",
              "alpha",
              "beta",
              "t_min",
              "t_max",
              "dm_min",
              "picky"},
      .gin_type = {"Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Index"},
      .gin_value = {std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    Numeric{916.7},
                    Numeric{-0.237},
                    Numeric{1.839},
                    std::nullopt,
                    std::nullopt,
                    Numeric{-1.0},
                    Index{0}},
      .gin_desc =
          {R"--(Ice water content)--",
           R"--(Intercept parameter)--",
           R"--(Volume weighted diameter)--",
           R"--(Density of ice)--",
           R"--(``alpha`` parameter of the shape function)--",
           R"--(``beta`` paramter of the shape function)--",
           R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Lower threshold for ``Dm`` below which an error is thrown.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdFieldEtAl07"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(The Field et al. [2007] particle size distribution for snow and
cloud ice.

This is a 1-parameter PSD, i.e. *pnd_agenda_input* shall have one
column and *pnd_agenda_input_names* shall contain a single string.
The input data in *pnd_agenda_input* shall be ice hydrometeor mass
content in unit of [kg/m3]. The naming used is *pnd_agenda_input_names*
is free but the same name must be used in *particle_bulkprop_names* and
*dpnd_data_dx_names*.

*psd_size_grid* shall contain size in terms of maximum diameter.

Derivatives are obtained by perturbation of 0.1%, but not less than
1e-9 kg/m3.

Both parametrization for tropics and midlatitudes are handled,
governed by setting of ``regime``, where "TR" selectes the tropical
case, and "ML" the midlatitude one.

The validity range of mass content is not limited. Negative mass
contents will produce negative psd values following a distribution
given by abs(IWC), ie. abs(psd)=f(abs(IWC)).

If temperature is outside [``t_min``,``t_max``] psd=0 and dpsd=0 if
picky=0, or an error is thrown if picky=1.

For temperatures below `t_min_psd``, the size distribution is
calculated for T = `t_min_psd``. Likewise, for temperatures above
``t_max_psd``, the distribution is derived for T = ``t_max_psd``.

Defaults of `t_min_psd`` and ``t_max_psd`` were set considering that
the parametrization has been derived from measurements over
temperatures of -60C to 0C.
Checks of the sanity of the mass-dimension relationship are performed
Errors are thrown if:

- Mass-dimension relation exponent *scat_species_b* is outside [``beta_min``, ``beta_max``].
)--",
      .author = {"Jana Mendrok"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names",
             "scat_species_a",
             "scat_species_b"},
      .gin = {"regime",
              "t_min",
              "t_max",
              "t_min_psd",
              "t_max_psd",
              "beta_min",
              "beta_max",
              "picky"},
      .gin_type = {"String",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Index"},
      .gin_value = {std::nullopt,
                    Numeric{0},
                    Numeric{290.},
                    Numeric{200.},
                    Numeric{273.15},
                    Numeric{1.01},
                    Numeric{4},
                    Index{0}},
      .gin_desc =
          {R"--(Parametrization regime ("TR"=tropical or "ML"=midlatitude).)--",
           R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Low temperature limit to use as paramtrization temperature.)--",
           R"--(High temperature limit to use as paramtrization temperature.)--",
           R"--(Low ``b`` limit (only if picky).)--",
           R"--(High ``b`` limit (only if picky).)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdFieldEtAl19"] = WorkspaceMethodInternalRecord{
      .desc = R"--(The Field [2019] particle size distribution for hail.

Reference: Field, Normalized hail particle size distributions from the T-28
storm-penetrating aircraft, JAMC, 2019

This is a 1-parmater PSD i.e. *pnd_agenda_input* shall have one column and
*pnd_agenda_input_names* shall contain a single string.
The input data in *pnd_agenda_input* shall be hail mass content in
unit of [kg/m3]. The naming used is *pnd_agenda_input_names* is free
but the same name must be used in *particle_bulkprop_names* and
*dpnd_data_dx_names*.
The parameters assume a constant effective density, i.e. scat_species_b pprox 3

Derivatives are obtained analytically.

The validity range of mass content is not limited. Negative mass
contents will produce negative psd values following a distribution
given by abs(HWC), ie. abs(psd)=f(abs(HWC)).

If temperature is outside [ ``t_min`` , ``t_max`` ] psd=0 and dpsd=0 if
picky=0, or an error is thrown if picky=1.
)--",
      .author = {"Stuart Fox"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names",
             "scat_species_a",
             "scat_species_b"},
      .gin = {"t_min", "t_max", "picky"},
      .gin_type = {"Numeric", "Numeric", "Index"},
      .gin_value = {std::nullopt, std::nullopt, Index{0}},
      .gin_desc =
          {R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdMcFarquaharHeymsfield97"] = WorkspaceMethodInternalRecord{
      .desc = R"--(McFarquahar and Heymsfield [1997] particle size distribution
for cloud ice.

This is a 1-parameter PSD, i.e. *pnd_agenda_input* shall have one
column and *pnd_agenda_input_names* shall contain a single string.
The input data in *pnd_agenda_input* shall be ice hydrometeor mass
content in unit of [kg/m3]. The naming used is *pnd_agenda_input_names*
is free but the same name must be used in *particle_bulkprop_names* and
*dpnd_data_dx_names*.

*psd_size_grid* shall contain size in terms of volume equivalent diameter.

Derivatives are obtained by perturbation of 0.1%, but not less than
1e-9 kg/m3.

The validity range of mass content is not limited. Negative mass
contents will produce negative psd values following a distribution
given by abs(IWC), ie. abs(psd)=f(abs(IWC)).

If temperature is outside [ ``t_min`` , ``t_max`` ] psd=0 and dpsd=0 if
picky=0, or an error is thrown if picky=1.

For temperatures below `t_min_psd``, the size distribution is
calculated for T = `t_min_psd``. Likewise, for temperatures above
``t_max_psd``, the distribution is derived for T = ``t_max_psd``.

Defaults of ``t_min_psd`` and ``t_max_psd`` were set considering that
the parametrization has been derived from measurements over
temperatures of -70C to -20C.
The noisy option can not be used together with calculation of
derivatives (ie. when *dpnd_data_dx_names* is not empty).
)--",
      .author = {"Patrick Eriksson, Jana Mendrok"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names",
             "scat_species_a",
             "scat_species_b"},
      .gin = {"t_min", "t_max", "t_min_psd", "t_max_psd", "picky", "noisy"},
      .gin_type =
          {"Numeric", "Numeric", "Numeric", "Numeric", "Index", "Index"},
      .gin_value = {Numeric{0},
                    Numeric{280.},
                    Numeric{180},
                    Numeric{273.15},
                    Index{0},
                    Index{0}},
      .gin_desc =
          {R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Low temperature limit to use as paramtrization temperature.)--",
           R"--(High temperature limit to use as paramtrization temperature.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--",
           R"--(Distribution parameter perturbance flag)--"},

  };

  wsm_data["psdMilbrandtYau05"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculates *psd_data* and  *dpsd_data_dx* following Milbrandt and Yau (2005)
two moment particle size distribution for cloud water, cloud ice,
rain, snow, graupel and hail, which is used in the GEM model.

WSM for use in *pnd_agenda_array* for mapping ``particle_bulkprop_field``
to *pnd_field* using *pnd_fieldCalcFromParticleBulkProps*.
Produces the particle size distribution values (dN/dD) and their
derivates with respect to independent variables x by *dpnd_data_dx_names*
over multiple particle sizes and atmospheric levels (or SWC/T
combinations).

*psd_size_grid* is considered to be in terms of maximum diameter.
WC is considered to be in terms of mass content (or mass density),
ie. units of [kg/m3]. N_tot in terms of number density, ie. units of [1/m3].

Derivatives with respect to WC and N_tot are obtained analytically.

Six particle size distributions for the different hydrometeors are handled,
governed by setting of ``hydrometeor_type``, where 

(1) "cloud_water" selects cloud liquid water , 
(2) "cloud_ice" selects cloud ice, 
(3) "snow" selects snow, 
(4) "rain" selects rain, 
(5) "graupel" selects graupel, and 
(6) "hail" selects hail, 

Requirements:

*pnd_agenda_input_names* must include::

    ["X-mass_density", "X-number_density" ]. "X" is an arbitrary name

The entries in  *dpnd_data_dx_names* (ie. the allowed
independent variablea ) can be "X-mass_density" and\or 
"X-number_density".

The validity range of WC is not limited. Negative WC will produce
negative psd values following a distribution given by abs(WC), ie.
abs(psd)=f(abs(WC)).

If temperature is outside [``t_min``,``t_max``] psd=0 and dpsd=0 if
picky=0, or an error is thrown if picky=1.
)--",
      .author = {"Manfred Brath"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names"},
      .gin = {"hydrometeor_type", "t_min", "t_max", "picky"},
      .gin_type = {"String", "Numeric", "Numeric", "Index"},
      .gin_value = {std::nullopt, Numeric{0}, Numeric{999}, Index{0}},
      .gin_desc =
          {R"--(Hydrometeor type (see above description).)--",
           R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdModifiedGamma"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Modified gamma distribution PSD using n0, mu, la and ga as parameters.

The modified gamma distribution is a 4-parameter (n0, mu, la and ga)
distribution [Petty & Huang, JAS, 2011)]::

  n(x) = n0 * x^mu * exp( -la*x^ga )

where x is particle size or mass.

The parameters can be given in two ways, either by *pnd_agenda_input* or
as GIN arguments. The first option allows the parameter to vary, while
in the second case the parameter gets a constant value. If a parameter is
part of *pnd_agenda_input*, the corresponding GIN argument must be set
to NaN (which is default). This means that the number of columns in
*pnd_agenda_input* and the number of non-NaN choices for n0, mu, la and
ga must add up to four.

Data in *pnd_agenda_input* are linked to the MGD parameters in term of
order, the naming in *pnd_agenda_input_names* is free. If all four
parameteras are specified by *pnd_agenda_input*, the data in the first
column are taken as n0, the second column as mu etc. If a parameter is
given as a GIN argument, the columns are just shifted with one position.
For example, if mu and ga are specified as GIN arguments, *pnd_agenda_input*
shall have two columns, with n0-values in the first one and la-values in
the second one.

The GIN route is especially suitable for selecting special cases of MGD.
For example, by setting mu=0 and ga=1, an exponential PSD is obtained::

  n(x) = n0 * exp( -la * x )

With mu=1 and ga=1, the gamma PSD is obtained::

  n(x) = n0 * x^mu * exp( -la * x )

There should be little overhead in using the method for exponential
and gamma PSDs, there is an internal switch to dedicated expressions
for those PSDs.

Derivatives can only be obtained for parameters that are specified by
*pnd_agenda_input*.

If temperature is outside [ ``t_min`` , ``t_max`` ] psd=0 and dpsd=0 if
picky=0, or an error is thrown if picky=1.

These requirements apply to the MGD parameters:

(1) la > 0
(2) ga > 0
)--",
      .author = {"Patrick Eriksson"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names"},
      .gin = {"n0", "mu", "la", "ga", "t_min", "t_max", "picky"},
      .gin_type = {"Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Index"},
      .gin_value = {std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::nullopt,
                    std::nullopt,
                    Index{0}},
      .gin_desc =
          {R"--(n0)--",
           R"--(mu)--",
           R"--(la)--",
           R"--(ga)--",
           R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdModifiedGammaMass"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Modified gamma distribution (MGD) PSD, with mass content as input.

See *psdModifiedGamma* for a defintion of MGD parameters and how
this PSD is handled in ARTS. Only deviations with respect to
*psdModifiedGamma* are described here.

This version of MGD PSD takes mass content as first input argument.
This means that the first column of *pnd_agenda_input* shall hold
mass content data.

The mass content basically replaces one of the standard parameters
(n0, mu, la and ga). This parameter is denoted as the dependent one.
The dependent parameter is selected by setting the corresponding GIN
to -999. So far only n0 and la are allowed to be dependent.

Regarding remaining columns in *pnd_agenda_input* and constant
parameter values (by GIN) follows the same principle as for
*psdModifiedGamma* except that mass is always in column one (as
mentioned) and that there is no position in *pnd_agenda_input*
for the dependent parameter.

These requirements apply to the MGD parameters:
 (1) mu + scat_species_b + 1 > 0
 (2) la > 0
 (3) ga > 0
 (4) If la is the dependent parameter, mass content must be > 0.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names",
             "scat_species_a",
             "scat_species_b"},
      .gin = {"n0", "mu", "la", "ga", "t_min", "t_max", "picky"},
      .gin_type = {"Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Index"},
      .gin_value = {std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::nullopt,
                    std::nullopt,
                    Index{0}},
      .gin_desc =
          {R"--(n0)--",
           R"--(mu)--",
           R"--(la)--",
           R"--(ga)--",
           R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdModifiedGammaMassMeanParticleMass"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Modified gamma distribution PSD, with mass content and mean particle
mass (Mmean) as inputs.

"Mean particle mass" is here defined as the mass content divided with
the total number density.

This version of MGD PSD works as *psdModifiedGammaMass*, but takes
mass content and mean particle mass as first two arguments. This
means that the first and second column of *pnd_agenda_input* shall
hold mass content and Mmean, respectively. Accordingly, the number
of dependent parameters is two.

These requirements apply to the MGD parameters:
 (1) mu + 1 > 0
 (2) la > 0
 (3) ga > 0
 (4) Mmean must be > 0.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names",
             "scat_species_a",
             "scat_species_b"},
      .gin = {"n0", "mu", "la", "ga", "t_min", "t_max", "picky"},
      .gin_type = {"Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Index"},
      .gin_value = {std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::nullopt,
                    std::nullopt,
                    Index{0}},
      .gin_desc =
          {R"--(n0)--",
           R"--(mu)--",
           R"--(la)--",
           R"--(ga)--",
           R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdModifiedGammaMassNtot"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Modified gamma distribution PSD, with mass content and total number
density (Ntot) as inputs.

This version of MGD PSD works as *psdModifiedGammaMass*, but takes
mass content and total number density as first two arguments. This
means that the first and second column of *pnd_agenda_input* shall
hold mass content and Ntot, respectively. Accordingly, the number
of dependent parameters is two.

These requirements apply:
 (1) mu + 1 > 0
 (2) la > 0
 (3) ga > 0
 (4) Ntot must be > 0.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names",
             "scat_species_a",
             "scat_species_b"},
      .gin = {"n0", "mu", "la", "ga", "t_min", "t_max", "picky"},
      .gin_type = {"Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Index"},
      .gin_value = {std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::nullopt,
                    std::nullopt,
                    Index{0}},
      .gin_desc =
          {R"--(n0)--",
           R"--(mu)--",
           R"--(la)--",
           R"--(ga)--",
           R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdModifiedGammaMassSingleMoment"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Modified gamma distribution PSD, with mass content as input.

The intercept parameter N0 is assumed dependent on the slope parameter
lambda, such that N0=N_alpha*lambda^n_b with fixed N_alpha and n_b.
This is a common form for many PSD parametrizations for use with
single-moment mass-based schemes.

This version of MGD PSD takes mass content as first input argument.
This means that the first column of *pnd_agenda_input* shall hold
mass content data. The dependent parameter is assumed to be lambda.
)--",
      .author = {"Stuart Fox"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names",
             "scat_species_a",
             "scat_species_b"},
      .gin = {"n_alpha", "n_b", "mu", "gamma", "t_min", "t_max", "picky"},
      .gin_type = {"Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Index"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    Index{0}},
      .gin_desc =
          {R"--(n_alpha)--",
           R"--(n_b)--",
           R"--(mu)--",
           R"--(gamma)--",
           R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdModifiedGammaMassXmean"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Modified gamma distribution PSD, with mass content and mean size
(Xmean) as inputs.

"Mean size" is here defined as mass weighted size. Remembering that
mass is a*x^b, this mean size can be expressed as M_b+1/M_b where M_b
is b:th moment of the PSD (see e.g. Eq. 17 in Petty&Huang, JAS, 2011).

This version of MGD PSD works as *psdModifiedGammaMass*, but takes
mass content and mass size as first two arguments. This means that
the first and second column of *pnd_agenda_input* shall hold mass
content and Xmean, respectively. Accordingly, the number of dependent
parameters is two.

These requirements apply to the MGD parameters:
 (1) mu + scat_species_b + 1 > 0
 (2) la > 0
 (3) ga > 0
 (4) Xmean must be > 0.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names",
             "scat_species_a",
             "scat_species_b"},
      .gin = {"n0", "mu", "la", "ga", "t_min", "t_max", "picky"},
      .gin_type = {"Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Index"},
      .gin_value = {std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::nullopt,
                    std::nullopt,
                    Index{0}},
      .gin_desc =
          {R"--(n0)--",
           R"--(mu)--",
           R"--(la)--",
           R"--(ga)--",
           R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdModifiedGammaMassXmedian"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Modified gamma distribution PSD, with mass content and median size
(Xmedian) as inputs.


This version of MGD PSD works as *psdModifiedGammaMass*, but takes
mass content and median size as first two arguments. This means that
the first and second column of *pnd_agenda_input* shall hold mass
content and Xmedian, respectively. Accordingly, the number of
dependent parameters is two.

These requirements apply to the MGD parameters:
 (1) mu + scat_species_b + 1 > 0
 (2) la > 0
 (3) ga > 0
 (4) Xmedian must be > 0.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names",
             "scat_species_a",
             "scat_species_b"},
      .gin = {"n0", "mu", "la", "ga", "t_min", "t_max", "picky"},
      .gin_type = {"Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Index"},
      .gin_value = {std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::nullopt,
                    std::nullopt,
                    Index{0}},
      .gin_desc =
          {R"--(n0)--",
           R"--(mu)--",
           R"--(la)--",
           R"--(ga)--",
           R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdMonoDispersive"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Mono-dispersive PSD, with number density given.

This is a 1-parameter PSD, i.e. *pnd_agenda_input* shall have one
column and *pnd_agenda_input_names* shall contain a single string.
The input data in *pnd_agenda_input* shall be number densities, in
unit of [#/m3]. The naming used is *pnd_agenda_input_names* is free
but the same name must be used in *particle_bulkprop_names* and
*dpnd_data_dx_names*.

The method checks that the scattering species indicated (by
``species_index``) has a single element, and just inserts the provided
number density in *psd_data*.

If temperature is outside [ ``t_min`` , ``t_max`` ] psd=0 and dpsd=0 if
picky=0, or an error is thrown if picky=1.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names",
             "scat_meta"},
      .gin = {"species_index", "t_min", "t_max", "picky"},
      .gin_type = {"Index", "Numeric", "Numeric", "Index"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, Index{0}},
      .gin_desc =
          {R"--(The index of the scattering species of concern (0-based).)--",
           R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdMonoMass"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Mono-dispersive PSD, with mass content given.

This is a 1-parameter PSD, i.e. *pnd_agenda_input* shall have one
column and *pnd_agenda_input_names* shall contain a single string.
The input data in *pnd_agenda_input* shall be mass contents, in
unit of [#/m3]. The naming used is *pnd_agenda_input_names* is free
but the same name must be used in *particle_bulkprop_names* and
*dpnd_data_dx_names*.

The method checks that the scattering species indicated (by
``species_index``) has a single element, and sets *psd_data* based
on the mass contents given and the particle mass (derived from
*scat_meta*).

If temperature is outside [ ``t_min`` , ``t_max`` ] psd=0 and dpsd=0 if
picky=0, or an error is thrown if picky=1.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names",
             "scat_meta"},
      .gin = {"species_index", "t_min", "t_max", "picky"},
      .gin_type = {"Index", "Numeric", "Numeric", "Index"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, Index{0}},
      .gin_desc =
          {R"--(The index of the scattering species of concern (0-based).)--",
           R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdSeifertBeheng06"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculates *psd_data* and *dpsd_data_dx* following Seifert and Beheng (2006)
two moment particle size distribution for cloud water, cloud ice,
rain, snow, graupel and hail, which is used in the ICON model.

WSM for use in *pnd_agenda_array* for mapping ``particle_bulkprop_field``
to *pnd_field* using *pnd_fieldCalcFromParticleBulkProps*.
Produces the particle size distribution values (dN/dD) and their
derivates with respect to independent variables x by *dpnd_data_dx_names*
over multiple particle sizes and atmospheric levels (or SWC/T
combinations).

*psd_size_grid* is considered to be in terms of mass.
WC is considered to be in terms of mass content (or mass density),
ie. units of [kg/m3]. N_tot in terms of number density, ie. units of [1/m3] .
Derivatives with respect to WC and N_tot are obtained analytically.

Six particle size distributions for the different hydrometeors are handled,
governed by setting of ``hydrometeor_type``, where 

(1) "cloud_water" selects cloud liquid water , 
(2) "cloud_ice" selects cloud ice, 
(3) "snow" selects snow, 
(4) "rain" selects rain, 
(5) "graupel" selects graupel, and 
(6) "hail" selects hail, 

Requirements:

*pnd_agenda_input_names* must include::

  ["X-mass_density", "X-number_density" ]. "X" is an arbitrary name

The entries in  *dpnd_data_dx_names* (ie. the allowed
independent variablea ) can be "X-mass_density" and\or 
"X-number_density".

The validity range of WC is not limited. Negative WC will produce
negative psd values following a distribution given by abs(WC), ie.
abs(psd)=f(abs(WC)).

If temperature is outside [ ``t_min`` , ``t_max`` ] psd=0 and dpsd=0 if
picky=0, or an error is thrown if picky=1.
)--",
      .author = {"Manfred Brath"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names"},
      .gin = {"hydrometeor_type", "t_min", "t_max", "picky"},
      .gin_type = {"String", "Numeric", "Numeric", "Index"},
      .gin_value = {std::nullopt, Numeric{0}, Numeric{999}, Index{0}},
      .gin_desc =
          {R"--(Hydrometeor type (see above description).)--",
           R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["psdWangEtAl16"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Wang et al. [2016] particle size distribution for rain.

Reference: Wang et al., Investigation of liquid cloud microphysical
properties of deep convective systems: 1. Parameterization raindrop
size distribution and its application ..., 2016.

This is a 1-parameter PSD, i.e. *pnd_agenda_input* shall have one
column and *pnd_agenda_input_names* shall contain a single string.
The input data in *pnd_agenda_input* shall be rain mass content in
unit of [kg/m3]. The naming used is *pnd_agenda_input_names* is free
but the same name must be used in *particle_bulkprop_names* and
*dpnd_data_dx_names*.

Particles are assumed to be near-spherical, ie. *psd_size_grid* can
either be in terms of volume (or mass) equivalent diameter or
maximum diameter.

Derivatives are obtained analytically.

The validity range of mass content is not limited. Negative mass
contents will produce negative psd values following a distribution
given by abs(RWC), ie. abs(psd)=f(abs(RWC)).

If temperature is outside [ ``t_min`` , ``t_max`` ] psd=0 and dpsd=0 if
picky=0, or an error is thrown if picky=1.
)--",
      .author = {"Jana Mendrok, Patrick Eriksson"},
      .out = {"psd_data", "dpsd_data_dx"},

      .in = {"psd_size_grid",
             "pnd_agenda_input_t",
             "pnd_agenda_input",
             "pnd_agenda_input_names",
             "dpnd_data_dx_names",
             "scat_species_a",
             "scat_species_b"},
      .gin = {"t_min", "t_max", "picky"},
      .gin_type = {"Numeric", "Numeric", "Index"},
      .gin_value = {Numeric{273}, Numeric{373}, Index{0}},
      .gin_desc =
          {R"--(Low temperature limit to calculate a psd.)--",
           R"--(High temperature limit to calculate a psd.)--",
           R"--(Flag whether to be strict with parametrization value checks.)--"},

  };

  wsm_data["refr_index_airFreeElectrons"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Microwave refractive index due to free electrons.

The refractive index of free electrons is added to *refr_index_air*.
To obtain the complete value, *refr_index_air* should be set to 1
before calling this WSM. This applies also to *refr_index_air_group*.

The expression applied is n=sqrt(1-wp^2/w^2) where wp is the plasma
frequency, and w is the angular frequency (the function returns
n-1, that here is slightly negative). This expressions is found in
many textbooks, e.g. Rybicki and Lightman (1979). The above refers
to *refr_index_air*. *refr_index_air_group* is sqrt(1+wp^2/w^2).

The expression is dispersive. The frequency applied is the mean of
first and last element of *f_grid* is selected. This frequency must
be at least twice the plasma frequency.

An error is issued if free electrons not are part of *abs_species*
(and there exist a corresponding "vmr"-value). This demand is
removed if ``demand_vmr_value`` is set to 0, but use this option
with care.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"refr_index_air", "refr_index_air_group"},

      .in = {"refr_index_air",
             "refr_index_air_group",
             "f_grid",
             "abs_species",
             "rtp_vmr"},
      .gin = {"demand_vmr_value"},
      .gin_type = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc =
          {R"--(Flag to control if it is demanded that free electrons are in *abs_species*. Default is that this is demanded.)--"},

  };

  wsm_data["refr_index_airInfraredEarth"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates the IR refractive index due to gases in the
Earth's atmosphere.

Only refractivity of dry air is considered. The formula used is
contributed by Michael Hoepfner, Forschungszentrum Karlsruhe.

The refractivity of dry air is added to *refr_index_air*. To obtain
the complete value, *refr_index_air* should be set to 1 before
calling this WSM. This applies also to *refr_index_air_group*.

The expression used is non-dispersive. Hence, *refr_index_air* and
*refr_index_air_group* are identical.
)--",
      .author = {"Mattias Ekstrom"},
      .out = {"refr_index_air", "refr_index_air_group"},

      .in = {"refr_index_air",
             "refr_index_air_group",
             "rtp_pressure",
             "rtp_temperature"},

  };

  wsm_data["refr_index_airMicrowavesEarth"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Microwave refractive index in Earth's atmosphere.

This method just considers pressure, temperature and water
vapour, which should suffice for Earth. For a more general
method, see *refr_index_airMicrowavesGeneral*.

The refractivity of dry air and water vapour is added to
*refr_index_air*. To obtain the complete value, *refr_index_air*
should be set to 1 before calling this WSM. This applies also to
*refr_index_air_group*.

The expression used is non-dispersive. Hence, *refr_index_air*
and *refr_index_air_group* are identical.

The standard expression for Earth and microwaves is used:

  N = k1*(P-e)/T + k2*e/T + k3*e/T^2

where N is refractivity, P is pressure, T is temperature and
e is water vapour partial pressure. The values of k1, k2 and k3
can be modified.

Many different values of k1, k2 and k3 can be found in the
literature. The default values applied here are taken from
Bevis et al., GPS meteorology: Mapping ..., JAM, 1994.
More specifically, these value are found in Table 1, listed
as "Present study". Note that in ARTS Pa is used for pressure
and k1, k2 and k3 must be adjusted accordingly.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"refr_index_air", "refr_index_air_group"},

      .in = {"refr_index_air",
             "refr_index_air_group",
             "rtp_pressure",
             "rtp_temperature",
             "rtp_vmr",
             "abs_species"},
      .gin = {"k1", "k2", "k3"},
      .gin_type = {"Numeric", "Numeric", "Numeric"},
      .gin_value = {Numeric{77.6e-8}, Numeric{70.4e-8}, Numeric{3.739e-3}},
      .gin_desc = {R"--(Coefficient a, see above)--",
                   R"--(Coefficient b, see above)--",
                   R"--(Coefficient c, see above)--"},

  };

  wsm_data["refr_index_airMicrowavesGeneral"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Microwave refractive index due to gases in planetary atmospheres.

The refractivity of a specified gas mixture is calculated and added
to *refr_index_air*. To obtain the complete value, *refr_index_air*
should be set to 1 before calling this WSM. This applies also to
*refr_index_air_group*.

The expression used is non-dispersive. Hence, *refr_index_air* and
*refr_index_air_group* are identical.

Uses the methodology introduced by Newell&Baird (1965) for calculating
refractivity of variable gas mixtures based on refractivity of the
individual gases at reference conditions. Assuming ideal gas law for
converting reference refractivity to actual pressure and temperature
conditions. Reference refractivities are also taken from Newell&Baird (1965)
and are vailable for N2, O2, CO2, H2, and He. Additionally, H2O reference
refractivity has been derived from H2O contribution in Thayer (see
*refr_index_airMicrowavesEarth*) for T0=273.15K. Any mixture of these gases
can be taken into account.
)--",
      .author = {"Jana Mendrok"},
      .out = {"refr_index_air", "refr_index_air_group"},

      .in = {"refr_index_air",
             "refr_index_air_group",
             "rtp_pressure",
             "rtp_temperature",
             "rtp_vmr",
             "abs_species"},

  };

  wsm_data["refr_index_air_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *refr_index_air_agenda* to a default value

Options are:

- ``"NoRefrac"``:

    1. Sets *refr_index_air* to 1
    2. Sets *refr_index_air_group* to 1

- ``"GasMicrowavesEarth"``:

    1. Sets *refr_index_air* to 1
    2. Sets *refr_index_air_group* to 1
    3. Uses *refr_index_airMicrowavesEarth* to modify *refr_index_air*, and *refr_index_air_group*

- ``"GasInfraredEarth"``:

    1. Sets *refr_index_air* to 1
    2. Sets *refr_index_air_group* to 1
    3. Uses *refr_index_airInfraredEarth* to modify *refr_index_air*, and *refr_index_air_group*

- ``"GasMicrowavesGeneral"``:

    1. Sets *refr_index_air* to 1
    2. Sets *refr_index_air_group* to 1
    3. Uses *refr_index_airMicrowavesGeneral* to modify *refr_index_air*, and *refr_index_air_group*

- ``"FreeElectrons"``:

    1. Sets *refr_index_air* to 1
    2. Sets *refr_index_air_group* to 1
    3. Uses *refr_index_airFreeElectrons* to modify *refr_index_air*, and *refr_index_air_group*

- ``"GasMicrowavesGeneralAndElectrons"``:

    1. Sets *refr_index_air* to 1
    2. Sets *refr_index_air_group* to 1
    3. Uses *refr_index_airMicrowavesGeneral* to modify *refr_index_air*, and *refr_index_air_group*
    4. Uses *refr_index_airFreeElectrons* to modify *refr_index_air*, and *refr_index_air_group*

- ``"GasMicrowavesEarthAndElectrons"``:

    1. Sets *refr_index_air* to 1
    2. Sets *refr_index_air_group* to 1
    3. Uses *refr_index_airMicrowavesEarth* to modify *refr_index_air*, and *refr_index_air_group*
    4. Uses *refr_index_airFreeElectrons* to modify *refr_index_air*, and *refr_index_air_group*
)--",
      .author = {"Richard Larsson"},
      .out = {"refr_index_air_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["retrievalAddAbsSpecies"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds an absorption species to the retrieval quantities.

Similar to *jacobianAddAbsSpecies* but also sets the corresponding block in
*covmat_sx* to the matrices provided in *covmat_block* and *covmat_inv_block*.
The dimensions of *covmat_block* are required to agree with the dimensions of the
retrieval grid.

*covmat_inv_block* must be either empty or the same dimension as *covmat_block*.
If provided, this matrix will be used as the inverse for the covariance matrix block
and numerical inversion of this block is thus avoided. Note, however, that this is
only effective if this block is uncorrelated with any other retrieval quantity.

For number and order of elements added to *x*, see *jacobianAddAbsSpecies*.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},

      .in = {"covmat_sx",
             "jacobian_quantities",
             "jacobian_agenda",
             "covmat_block",
             "covmat_inv_block"},
      .gin = {"g1", "g2", "g3", "species", "unit", "for_species_tag"},
      .gin_type = {"Vector", "Vector", "Vector", "String", "String", "Index"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    String("rel"),
                    Index{1}},
      .gin_desc = {R"--(Pressure retrieval grid.)--",
                   R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retreival grid.)--",
                   R"--(The species tag of the retrieval quantity.)--",
                   R"--(Retrieval unit. See above.)--",
                   R"--(Index-bool for acting on species tags or species.)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalAddCatalogParameter"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Similar to *jacobianAddBasicCatalogParameter* but also adds a corresponding
block to *covmat_sx* with the given ``var`` as variance value.

For number and order of elements added to *x*,
see *jacobianAddBasicCatalogParameter*.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},

      .in = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},
      .gin = {"catalog_identity", "catalog_parameter", "var"},
      .gin_type = {"QuantumIdentifier", "String", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(The catalog line matching information.)--",
                   R"--(The catalog parameter of the retrieval quantity.)--",
                   R"--(The variance of the catalog parameter.)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalAddCatalogParameters"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Same as *jacobianAddBasicCatalogParameters* but also adds a new
block to *covmat_sx* using the matrices in *covmat_block* and
*covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*,
see *jacobianAddBasicCatalogParameters*.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},

      .in = {"covmat_sx",
             "jacobian_quantities",
             "jacobian_agenda",
             "covmat_block",
             "covmat_inv_block"},
      .gin = {"catalog_identities", "catalog_parameters"},
      .gin_type = {"ArrayOfQuantumIdentifier", "ArrayOfString"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(The catalog line matching informations.)--",
                   R"--(The catalog parameters of the retrieval quantity.)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalAddFreqShift"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Same as *jacobianAddFreqShift* but also adds the correlation block
contained in *covmat_block* and *covmat_inv_block* to *covmat_sx*.

For number and order of elements added to *x*, see *jacobianAddFreqShift*.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},

      .in = {"covmat_sx",
             "covmat_block",
             "covmat_inv_block",
             "jacobian_quantities",
             "jacobian_agenda",
             "f_grid"},
      .gin = {"df"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{100e3}},
      .gin_desc = {R"--(Size of perturbation to apply.)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalAddFreqStretch"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Same as *jacobianAddFreqShift* but also adds the correlation block
contained in *covmat_block* and *covmat_inv_block* to *covmat_sx*.

For number and order of elements added to *x*, see *jacobianAddFreqStretch*.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},

      .in = {"jacobian_quantities",
             "jacobian_agenda",
             "f_grid",
             "covmat_block",
             "covmat_inv_block"},
      .gin = {"df"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{100e3}},
      .gin_desc = {R"--(Size of perturbation to apply.)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalAddMagField"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Same as *jacobianAddMagField* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddMagField*.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},

      .in = {"covmat_sx",
             "jacobian_quantities",
             "jacobian_agenda",
             "covmat_block",
             "covmat_inv_block"},
      .gin = {"g1", "g2", "g3", "component", "dB"},
      .gin_type = {"Vector", "Vector", "Vector", "String", "Numeric"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    String("v"),
                    Numeric{1.0e-7}},
      .gin_desc = {R"--(Pressure retrieval grid.)--",
                   R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retreival grid.)--",
                   R"--(Magnetic field component to retrieve)--",
                   R"--(Magnetic field perturbation)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalAddPointingZa"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Same as *jacobianAddPointingZa* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddPointingZa*.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},

      .in = {"covmat_sx",
             "jacobian_quantities",
             "jacobian_agenda",
             "covmat_block",
             "covmat_inv_block",
             "sensor_pos",
             "sensor_time"},
      .gin = {"poly_order", "calcmode", "dza"},
      .gin_type = {"Index", "String", "Numeric"},
      .gin_value = {Index{0}, String("recalc"), Numeric{0.01}},
      .gin_desc =
          {R"--(Order of polynomial to describe the time variation of pointing off-sets.)--",
           R"--(Calculation method. See above)--",
           R"--(Size of perturbation to apply (when applicable).)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalAddPolyfit"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Same as *jacobianAddPolyfit* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddPolyfit*.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},

      .in = {"covmat_sx",
             "jacobian_quantities",
             "jacobian_agenda",
             "covmat_block",
             "covmat_inv_block",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid",
             "sensor_pos"},
      .gin = {"poly_order",
              "no_pol_variation",
              "no_los_variation",
              "no_mblock_variation"},
      .gin_type = {"Index", "Index", "Index", "Index"},
      .gin_value = {std::nullopt, Index{0}, Index{0}, Index{0}},
      .gin_desc =
          {R"--(Polynomial order to use for the fit.)--",
           R"--(Set to 1 if the baseline off-set is the same for all Stokes components.)--",
           R"--(Set to 1 if the baseline off-set is the same for all line-of-sights (inside each measurement block).)--",
           R"--(Set to 1 if the baseline off-set is the same for all measurement blocks.)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalAddScatSpecies"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Same as *jacobianAddPolyfit* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddScatSpecies*.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},

      .in = {"covmat_sx",
             "jacobian_quantities",
             "jacobian_agenda",
             "covmat_block",
             "covmat_inv_block"},
      .gin = {"g1", "g2", "g3", "species", "quantity"},
      .gin_type = {"Vector", "Vector", "Vector", "String", "String"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt},
      .gin_desc =
          {R"--(Pressure retrieval grid.)--",
           R"--(Latitude retrieval grid.)--",
           R"--(Longitude retreival grid.)--",
           R"--(Name of scattering species, must match one element in *scat_species*.)--",
           R"--(Retrieval quantity, e.g. "IWC".)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalAddSinefit"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Same as *jacobianAddSinefit* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddSinefit*.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},

      .in = {"covmat_sx",
             "jacobian_quantities",
             "jacobian_agenda",
             "covmat_block",
             "covmat_inv_block",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid",
             "sensor_pos"},
      .gin = {"period_lengths",
              "no_pol_variation",
              "no_los_variation",
              "no_mblock_variation"},
      .gin_type = {"Vector", "Index", "Index", "Index"},
      .gin_value = {std::nullopt, Index{0}, Index{0}, Index{0}},
      .gin_desc =
          {R"--(Period lengths of the fit.)--",
           R"--(Set to 1 if the baseline off-set is the same for all Stokes components.)--",
           R"--(Set to 1 if the baseline off-set is the same for all line-of-sights (inside each measurement block).)--",
           R"--(Set to 1 if the baseline off-set is the same for all measurement blocks.)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalAddSpecialSpecies"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Same as *jacobianAddSpecialSpecies* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddSpecialSpecies*.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},

      .in = {"covmat_sx",
             "jacobian_quantities",
             "jacobian_agenda",
             "covmat_block",
             "covmat_inv_block"},
      .gin = {"g1", "g2", "g3", "species"},
      .gin_type = {"Vector", "Vector", "Vector", "String"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Pressure retrieval grid.)--",
                   R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retreival grid.)--",
                   R"--(The species of the retrieval quantity.)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalAddSurfaceQuantity"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Same as *jacobianAddSurfaceQuantity* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddSurfaceQuantity*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},

      .in = {"covmat_sx",
             "jacobian_quantities",
             "jacobian_agenda",
             "covmat_block",
             "covmat_inv_block"},
      .gin = {"g1", "g2", "quantity"},
      .gin_type = {"Vector", "Vector", "String"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retreival grid.)--",
                   R"--(Retrieval quantity, e.g. "Wind speed".)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalAddTemperature"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Same as *jacobianAddTemperature* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddTemperature*.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},

      .in = {"covmat_sx",
             "jacobian_quantities",
             "jacobian_agenda",
             "covmat_block",
             "covmat_inv_block"},
      .gin = {"g1", "g2", "g3", "hse"},
      .gin_type = {"Vector", "Vector", "Vector", "String"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, String("on")},
      .gin_desc = {R"--(Pressure retrieval grid.)--",
                   R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retreival grid.)--",
                   R"--(Flag to assume HSE or not ("on" or "off").)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalAddWind"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Same as *jacobianAddWind* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddWind*.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_sx", "jacobian_quantities", "jacobian_agenda"},

      .in = {"covmat_sx",
             "jacobian_quantities",
             "jacobian_agenda",
             "covmat_block",
             "covmat_inv_block"},
      .gin = {"g1", "g2", "g3", "component", "dfrequency"},
      .gin_type = {"Vector", "Vector", "Vector", "String", "Numeric"},
      .gin_value =
          {std::nullopt, std::nullopt, std::nullopt, String("v"), Numeric{0.1}},
      .gin_desc = {R"--(Pressure retrieval grid.)--",
                   R"--(Latitude retrieval grid.)--",
                   R"--(Longitude retrieval grid.)--",
                   R"--(Wind component to retrieve)--",
                   R"--(This is the frequency perturbation)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalDefClose"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Closes the definition of retrieval quantities and correlations and
prepares related WSVs for the retrieval.

This function calls jacobianClose and checks that the corvariance matrices
are consistent with the Jacobian.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"jacobian_do", "jacobian_agenda", "retrieval_checked"},

      .in = {"jacobian_agenda", "covmat_sx", "jacobian_quantities"},

      .pass_workspace = true,

  };

  wsm_data["retrievalDefInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Begin retrieval definition section.

This function initialises all variables required for defining
retrieval quantities and corresponding covariance matrices.
By default, Jacobian quantities should be added withing the.
retrieval definition section. If Jacobian quantities are
defined separately ``initialize_jacobian`` must be set to 0,
otherwise the quantities will be discarded.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"covmat_se",
              "covmat_sx",
              "covmat_block",
              "covmat_inv_block",
              "jacobian_quantities",
              "jacobian_agenda"},

      .gin = {"initialize_jacobian"},
      .gin_type = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc =
          {R"--(Flag whether or not to (re)initialize Jacobian-related quantities. Set to 0 if Jacobian is already defined.)--"},
      .pass_workspace = true,

  };

  wsm_data["retrievalErrorsExtract"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Extract retrieval error from covariance matrices.

Extracts the error estimates for the retrieved quantities from the covariance
matrices for the error due to measurement noise *covmat_so* and the error due
to limited resolution of the observation system *covmat_ss* and stores them in
the vectors *retrieval_eo* and *retrieval_ss*, respectively.
To etract these errors, first the convariance matrices of which the errors 
should be extracted have to be computed using the WSMs *covmat_soCalc*
and *covmat_ssCalc* or set to be empty in order to be ignored. Note, however,
that this will also set the corresponding error vector to be empty.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"retrieval_eo", "retrieval_ss"},

      .in = {"covmat_so", "covmat_ss"},

  };

  wsm_data["rte_background_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *rte_background_agenda* to a default value

Options are:
    FIXME
)--",
      .author = {"Richard Larsson"},
      .out = {"rte_background_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["rte_losGeometricToPosition"] = WorkspaceMethodInternalRecord{
      .desc = R"--(The geometric line-of-sight between two points.

The line-of-sight angles from *rte_pos* to ``target_pos`` are calculated
ignoring refraction. This can be done analytically. The angles are set
without any consideration of the surface. The corresponding propagation
path can thus end with a surface intersection.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"rte_los"},

      .in = {"surface_field", "rte_pos"},
      .gin = {"target_pos"},
      .gin_type = {"Vector"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(The atmospheric position that *rte_los* shall match.)--"},

  };

  wsm_data["rte_losRefractedToPosition"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(The refracted line-of-sight and propagation path between two positions.

General remarks: The method assumes geometrical optics. It can fail for
various reasons, especially if there are sharp or multiple gradients in
refractive index. In multi-path situations, the method finds (in best
case) only a single solution. Default is to consider an intersection
with the surface as a failure and give an error about it. To instead
return an empty ppath when the target appears to be behind the horizon,
set GIN ``robust`` to 1.

The line-of-sight connecting the two points cannot be determined
analytically and a search algorithm is needed. The algorithm options are:

"basic":
A simple iteration scheme with no safety measures. In short, the search
is done by trying to establish the geometric target point that gives the
same path length and a *rte_los* that gives a hit when doing a refracted
path. In each iteration, the geometric target is moved according to the
deviation to ``target_pos`` of the latest refracted path. The convergence
can fail for long distances or with strong refractine index gradients.
Accordingly, the algorithm tends to fail for limb sounding below about
6 km, but should be applicable for other observation geometries (but
notice the general remarks).

Yes, so far only one option! Hopefully there will be more as the basic
option does not handle all cases properly.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"rte_los", "ppath"},

      .in = {"refr_index_air_ZZZ_agenda",
             "ppath_lstep",
             "ppath_lraytrace",
             "surface_field",
             "surface_search_accuracy",
             "rte_pos"},
      .gin = {"target_pos",
              "target_dl",
              "algorithm",
              "max_iterations",
              "robust",
              "z_toa",
              "do_horizontal_gradients",
              "do_twosided_perturb"},
      .gin_type = {"Vector",
                   "Numeric",
                   "String",
                   "Index",
                   "Index",
                   "Numeric",
                   "Index",
                   "Index"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    String("basic"),
                    Index{10},
                    Index{0},
                    std::nullopt,
                    Index{0},
                    Index{0}},
      .gin_desc =
          {R"--(The atmospheric position that *ppath* shall reach.)--",
           R"--(The end point of *ppath* shall be inside this distance from ``target_pos`` (deviation can be in any direction).)--",
           R"--(Search algorithm to use.)--",
           R"--(Max number of iterations before giving up.)--",
           R"--(Set to 1 to not give errors, but return empty *ppath* when a path can not be established.)--",
           R"--(Top-of-the-atmosphere altitude.)--",
           R"--(Consider horisontal gradients of refractive index.)--",
           R"--(Perform double-sided perturbations when calculating refractive index gradients.)--"},
      .pass_workspace = true,

  };

  wsm_data["rte_losReverse"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reverses the direction in *rte_los*.

The method updates *rte_los* to have angles of the reversed
direction.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"rte_los"},

      .in = {"rte_los"},

  };

  wsm_data["rte_losSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *rte_los* to the given angles.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"rte_los"},

      .gin = {"za", "aa"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Zenith angle [0, 180])--",
                   R"--(Azimuth angle [-180, 180])--"},

  };

  wsm_data["rte_posSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *rte_pos* to the given coordinates.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"rte_pos"},

      .gin = {"z", "glat", "glon"},
      .gin_type = {"Numeric", "Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(Altitude)--",
                   R"--(Latitude [-90, 90])--",
                   R"--(Longitude [-180, 360])--"},

  };

  wsm_data["rte_pos_losBackwardToAltitude"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Moves *rte_pos* and *rte_los* backwards to the target altitude.

The method gives the *rte_pos* and *rte_los* at the target altitude
to reach the original *rte_pos* and *rte_los* with a geometrical ppath.
That is, the movement is backwards in terms of viewing direction.

If the original *rte_los* is reversed with respect to the line-of-sight
direction, then set the GIN los_reversed to 1. One such case is that
if *rte_los* represents surface incidence angles, i.e. holds the
zenith and nadir angle towards the sensor.

There is also *sensor_pos_losBackwardToAltitude*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"rte_pos", "rte_los"},

      .in = {"rte_pos", "rte_los", "surface_field"},
      .gin = {"altitude", "los_is_reversed"},
      .gin_type = {"Numeric", "Index"},
      .gin_value = {std::nullopt, Index{0}},
      .gin_desc =
          {R"--(Target altitude.)--",
           R"--(Set to 1 if *rte_los* is valid for the reversed direction.)--"},

  };

  wsm_data["rte_pos_losEndOfPpath"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets *rte_pos* and *rte_los* to values for last point in *ppath*.

For example, if the propagation path intersects with the surface,
this method gives you the position and angle of *ppath* at the
surface.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"rte_pos", "rte_los"},

      .in = {"ppath"},

  };

  wsm_data["rte_pos_losForwardToAltitude"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Moves *rte_pos* and *rte_los* forward to the target altitude.

The method gives the *rte_pos* and *rte_los* at the target altitude
when forward-propagating the original *rte_pos* and *rte_los*
geometrically.

There is also *sensor_pos_losForwardToAltitude*. The WSM
*IntersectionGeometricAltitude* performs the same operation
with *sensor_pos* and *sensor_los* as input.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"rte_pos", "rte_los"},

      .in = {"rte_pos", "rte_los", "surface_field"},
      .gin = {"altitude"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Target altitude.)--"},

  };

  wsm_data["scat_dataCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Prepares *scat_data* for the scattering solver.

Derives single scattering data for the frequencies given by
*f_grid* by interpolation from *scat_data_raw*. *f_grid* should be
the actual WSV *f_grid* or a single-element Vector.
)--",
      .author = {"Jana Mendrok"},
      .out = {"scat_data"},

      .in = {"scat_data_raw", "f_grid"},
      .gin = {"interp_order"},
      .gin_type = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc = {R"--(Interpolation order.)--"},

  };

  wsm_data["scat_dataCheck"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Method for checking the validity and consistency of the single
scattering properties in *scat_data*.

It checks that *scat_data* does not contain any invalid values,
that is any NaN elements in K, Z, or a or any negative values in
the 'scalar' properties K11, Z11, and a1.

When ``check_type`` is 'all', it is furthermore checked that the
scattering matrix is properly normalized, that is that the solid
sphere integrated scattering matrix (int_Z11), which is supposed to
be normalized to the scattering cross section, is sufficiently
consistent with the scattering cross section (C_sca) derived from
the difference of extinction (K11) and absorption (a1):
int_z11 ~ C_sca = K11-a1.
Sufficient consistency is defined by the maximum allowed deviation
in single scattering albedo, ``sca_mat_threshold``, testing for::

  ( <int_Z11>/<C_sca>-1. ) * ( <C_sca>/<K11> ) <= sca_mat_threshold.

The check is skipped if ``check_type`` is 'sane'.
)--",
      .author = {"Claudia Emde", "Jana Mendrok"},

      .in = {"scat_data"},
      .gin = {"check_type", "sca_mat_threshold"},
      .gin_type = {"String", "Numeric"},
      .gin_value = {String("all"), Numeric{5e-2}},
      .gin_desc =
          {R"--(The level of checks to apply on scat_data ('sane' or 'all'; see above).)--",
           R"--(Threshold for allowed albedo deviation (see above).)--"},

  };

  wsm_data["scat_dataReduceT"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Reduces temperature dimension of single scattering to a single entry.

FIXME...
Derives single scattering data for the frequencies given by
*f_grid* by interpolation from *scat_data*. *f_grid* should be
the actual WSV *f_grid* or a single-element Vector.
)--",
      .author = {"Jana Mendrok"},
      .out = {"scat_data"},

      .in = {"scat_data"},
      .gin = {"scat_index",
              "temperature",
              "interp_order",
              "phamat_only",
              "sca_mat_threshold"},
      .gin_type = {"Index", "Numeric", "Index", "Index", "Numeric"},
      .gin_value =
          {std::nullopt, std::nullopt, Index{1}, Index{1}, Numeric{5e-2}},
      .gin_desc =
          {R"--(Apply on *scat_data* from scattering species of this index (0-based).)--",
           R"--(Temperature to interpolate *scat_data* to.)--",
           R"--(Interpolation order.)--",
           R"--(Flag whether to apply temperture reduction on phase matrix data only (1) or on all single scattering properties (0).)--",
           R"--(Threshold for allowed albedo deviation.)--"},

  };

  wsm_data["scat_data_checkedCalc"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Checks dimensions, grids and single scattering properties of all
scattering elements in *scat_data*.

Dimension and grid equirements:

- The scattering element's f_grid is either identical to *f_grid* or
  of dimension 1.
- In the latter case, the scattering element's f_grid value must
  not deviate from any of the *f_grid* values by more than a
  fraction of ``dfrel_threshold``.
- The frequency dimension of pha_mat_data, ext_mat_data, and
  abs_vec_data is either equal to the scattering element's f_grid
  or 1.
- The temperature dimension of pha_mat_data, ext_mat_data, and
  abs_vec_data is either equal to the scattering element's T_grid
  or 1.
- The temperature dimension of ext_mat_data, and abs_vec_data is
  identical.

The single scattering property contents are checked using
*scat_dataCheck*. For details, see there. The depth of these checks
and their rigour can adapted (see description of parameters
``check_level`` and ``sca_mat_threshold`` in *scat_dataCheck*) or can
be skipped entirely (setting ``check_level`` to 'none').
NOTE: These test shall only be skipped when one is confident that
the data is correct, e.g. by having run *scat_dataCheck* on the set
of data before, e.g. in a separate ARTS run.
)--",
      .author = {"Jana Mendrok"},
      .out = {"scat_data_checked"},

      .in = {"scat_data", "f_grid"},
      .gin = {"dfrel_threshold", "check_level", "sca_mat_threshold"},
      .gin_type = {"Numeric", "String", "Numeric"},
      .gin_value = {Numeric{0.1}, String("all"), Numeric{5e-2}},
      .gin_desc =
          {R"--(Maximum relative frequency deviation between (single entry) scattering element f_grid values and the RT calculation's *f_grid*.)--",
           R"--(See ``check_level`` in *scat_dataCheck*.)--",
           R"--(See ``sca_mat_threshold`` in *scat_dataCheck*.)--"},

  };

  wsm_data["scat_data_monoCalc"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Interpolates *scat_data* by frequency to give *scat_data_mono*.
)--",
      .author = {"Cory Davis"},
      .out = {"scat_data_mono"},

      .in = {"scat_data", "f_grid", "f_index"},

  };

  wsm_data["scat_data_monoExtract"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Extracts data at *f_index* from *scat_data* to give *scat_data_mono*.
)--",
      .author = {"Jana Mendrok"},
      .out = {"scat_data_mono"},

      .in = {"scat_data", "f_index"},

  };

  wsm_data["scat_data_singleTmatrix"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(A basic interface to Mishchenko's T-matrix code linked to ARTS.

The method performs a T-matrix calculation for a single scattering
element, i.e. a combination of particle shape, size, aspect ratio
and orientation.

Particle shape (``shape``) has two options::

  "spheroidal" and "cylindrical"

Particle size (``diameter_volume_equ``) is given as the equivalent
volume sphere diameter. That is, the diameter obtained if all the
particle's material is rearranged into a (solid) sphere.

Particle aspect ratio ar (``aspect_ratio``) is a numeric value, defined
according to Mishchenko's definition as ratio of horizontal axis a to
vertical (rotational) axis b: ar=a/b. That is, oblates have ar>1,
prolates ar<1.
Perfect spheres (spheroidals with ar=1) can trigger numerical issues.
To avoid these, we internally increase their aspect ratio by 1e-6,
i.e. turning perfect spheres into very light oblates.

Particle type (``ptype``) has two options::

  "totally_random" and "azimuthally_random"

For totally randomly oriented particles, ``data_aa_grid`` is not taken
into account (but a Vector type container needs to be passed).

For further information on how aspect ratio and the different shapes
and orientations are defined, see the documentation of the T-matrix
code found http://www.giss.nasa.gov/staff/mmishchenko/t_matrix.html

Regarding ``ndgs``, we refer to the this comment from the documentation:
   "Parameter controlling the number of division points
   in computing integrals over the particle surface.
   For compact particles, the recommended value is 2.
   For highly aspherical particles larger values (3, 4,...)
   may be necessary to obtain convergence.
   The code does not check convergence over this parameter.
   Therefore, control comparisons of results obtained with
   different NDGS-values are recommended."
)--",
      .author = {"Johan Strandgren", "Patrick Eriksson"},
      .out = {"scat_data_single", "scat_meta_single"},

      .in = {"complex_refr_index"},
      .gin = {"shape",
              "diameter_volume_equ",
              "aspect_ratio",
              "mass",
              "ptype",
              "data_f_grid",
              "data_t_grid",
              "data_za_grid",
              "data_aa_grid",
              "precision",
              "cri_source",
              "ndgs",
              "robust",
              "quiet"},
      .gin_type = {"String",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "String",
                   "Vector",
                   "Vector",
                   "Vector",
                   "Vector",
                   "Numeric",
                   "String",
                   "Index",
                   "Index",
                   "Index"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::numeric_limits<Numeric>::quiet_NaN(),
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    Vector{},
                    Numeric{0.001},
                    String("Set by user, unknown source."),
                    Index{2},
                    Index{0},
                    Index{1}},
      .gin_desc =
          {R"--(Particle shape. Options listed above.)--",
           R"--(Particle volume equivalent diameter [m]. See defintion above.)--",
           R"--(Particle aspect ratio.)--",
           R"--(Particle mass. This information is just included in the meta data, and does not affect the T-matrix calculations.)--",
           R"--(Particle type/orientation. Options listed above.)--",
           R"--(Frequency grid of the scattering data to be calculated.)--",
           R"--(Temperature grid of the scattering data to be calculated.)--",
           R"--(Zenith angle grid of the scattering data to be calculated.)--",
           R"--(Azimuth angle grid of the scattering data to be calculated.)--",
           R"--(Accuracy of the computations.)--",
           R"--(String describing the source of *complex_refr_index*, for inclusion in meta data.)--",
           R"--(See above. So far only applied for random orientation.)--",
           R"--(Continue even if individual T-matrix calculations fail. Respective scattering element data will be NAN.)--",
           R"--(Suppress print output from tmatrix fortran code.)--"},

  };

  wsm_data["sensorOff"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets sensor WSVs to obtain monochromatic pencil beam values.

The variables are set as follows:
 - *mblock_dlos*        : One row with zero(s).
 - *sensor_response*        : As returned by *sensor_responseInit*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_f_grid",
              "sensor_response_pol_grid",
              "sensor_response_dlos_grid",
              "mblock_dlos"},

      .in = {"f_grid"},

  };

  wsm_data["sensor_checkedCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Checks consistency of the sensor variables.

The following WSVs are examined: *f_grid*, *sensor_pos*, *sensor_los*,
*transmitter_pos*, *mblock_dlos*, *antenna_dim*,
*sensor_response*, *sensor_response_f*, *sensor_response_pol*,
and *sensor_response_dlos*.

If any of these variables are changed, then this method shall be
called again (no automatic check that this is fulfilled!).

The main tests are that dimensions of sensor variables agree
with other settings, e.g., the size of f_grid, atmosphere_dim,
stokes_dim, etc.

If any test fails, there is an error. Otherwise, *sensor_checked*
is set to 1.
)--",
      .author = {"Jana Mendrok"},
      .out = {"sensor_checked"},

      .in = {"f_grid",
             "sensor_pos",
             "sensor_los",
             "transmitter_pos",
             "mblock_dlos",
             "sensor_response",
             "sensor_response_f",
             "sensor_response_pol",
             "sensor_response_dlos"},

  };

  wsm_data["sensor_losAddLosAndDlos"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds zenith and azimuth angles.

Adds up a line-of-sights (ref_los), with relative angle off-sets
(dlos).
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_los"},

      .gin = {"ref_los", "gdlos"},
      .gin_type = {"Vector", "Matrix"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Reference line-of-sight (a single los).)--",
                   R"--(Change in line-of-sight (can be multiple dlos).)--"},

  };

  wsm_data["sensor_losGeometricToPosition"] = WorkspaceMethodInternalRecord{
      .desc = R"--(The geometric line-of-sight to a point.

Works as *rte_losGeometricToPosition*, but sets *sensor_los*.

This method handles the case of a single target position. For
multiple target positions, use: *sensor_losGeometricToPositions*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_los"},

      .in = {"surface_field", "sensor_pos"},
      .gin = {"target_pos"},
      .gin_type = {"Vector"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(The atmospheric position that *sensor_los* shall match.)--"},

  };

  wsm_data["sensor_losGeometricToPositions"] = WorkspaceMethodInternalRecord{
      .desc = R"--(The geometric line-of-sight to multiple point.

Works as *rte_losGeometricToPosition*, but sets *sensor_los*. The
number of rows in *sensor_pos* and ``target_pos`` must be equal.

This method handles the case of mutiple target positions. For
a single target positions, use: *sensor_losGeometricToPosition*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_los"},

      .in = {"surface_field", "sensor_pos"},
      .gin = {"target_pos"},
      .gin_type = {"Matrix"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(The atmospheric positions that *sensor_los* shall match.)--"},

  };

  wsm_data["sensor_losRefractedToPosition"] = WorkspaceMethodInternalRecord{
      .desc = R"--(The refracted line-of-sight to a point.

Works as *rte_losRefractedToPosition*, but sets *sensor_los*.

This method handles the case of a single target position. For
multiple target positions, use: *sensor_losRefractedToPositions*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_los"},

      .in = {"refr_index_air_ZZZ_agenda",
             "ppath_lstep",
             "ppath_lraytrace",
             "surface_field",
             "surface_search_accuracy",
             "sensor_pos"},
      .gin = {"target_pos",
              "target_dl",
              "algorithm",
              "max_iterations",
              "robust",
              "z_toa",
              "do_horizontal_gradients",
              "do_twosided_perturb"},
      .gin_type = {"Vector",
                   "Numeric",
                   "String",
                   "Index",
                   "Index",
                   "Numeric",
                   "Index",
                   "Index"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    String("basic"),
                    Index{10},
                    Index{0},
                    std::nullopt,
                    Index{0},
                    Index{0}},
      .gin_desc =
          {R"--(The atmospheric position that *ppath* shall reach.)--",
           R"--(The end point of *ppath* shall be inside this distance from ``target_pos`` (deviation can be in any direction).)--",
           R"--(Search algorithm to use.)--",
           R"--(Max number of iterations before giving up.)--",
           R"--(Set to 1 to not give errors, but return empty *ppath* when a path can not be established.)--",
           R"--(Top-of-the-atmosphere altitude.)--",
           R"--(Consider horisontal gradients of refractive index.)--",
           R"--(Perform double-sided perturbations when calculating refractive index gradients.)--"},
      .pass_workspace = true,

  };

  wsm_data["sensor_losRefractedToPositions"] = WorkspaceMethodInternalRecord{
      .desc = R"--(The refracted line-of-sight to multiple points.

Works as *rte_losRefractedToPosition*, but sets *sensor_los*.

This method handles the case of multiple target positions. For
a single target position, use: *sensor_losRefractedToPosition*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_los"},

      .in = {"refr_index_air_ZZZ_agenda",
             "ppath_lstep",
             "ppath_lraytrace",
             "surface_field",
             "surface_search_accuracy",
             "sensor_pos"},
      .gin = {"target_pos",
              "target_dl",
              "algorithm",
              "max_iterations",
              "robust",
              "z_toa",
              "do_horizontal_gradients",
              "do_twosided_perturb"},
      .gin_type = {"Matrix",
                   "Numeric",
                   "String",
                   "Index",
                   "Index",
                   "Numeric",
                   "Index",
                   "Index"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    String("basic"),
                    Index{10},
                    Index{0},
                    std::nullopt,
                    Index{0},
                    Index{0}},
      .gin_desc =
          {R"--(The atmospheric positions that *ppath* shall reach.)--",
           R"--(The end point of *ppath* shall be inside this distance from ``target_pos`` (deviation can be in any direction).)--",
           R"--(Search algorithm to use.)--",
           R"--(Max number of iterations before giving up.)--",
           R"--(Set to 1 to not give errors, but return empty *ppath* when a path can not be established.)--",
           R"--(Top-of-the-atmosphere altitude.)--",
           R"--(Consider horisontal gradients of refractive index.)--",
           R"--(Perform double-sided perturbations when calculating refractive index gradients.)--"},
      .pass_workspace = true,

  };

  wsm_data["sensor_losReverse"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reverses the directions in *sensor_los*.

The method updates *sensor_los* to have angles of the reversed
direction.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_los"},

      .in = {"sensor_los"},

  };

  wsm_data["sensor_pos_losBackwardToAltitude"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Moves *sensor_pos* and *sensor_los* backwards to the target altitude.

The method gives the *sensor_pos* and *sensor_los* at the target altitude
to reach the original *sensor_pos* and *sensor_los* with a geometrical
ppath. That is, the movement is backwards in terms of viewing direction.

If the original *sensor_los* is reversed with respect to the line-of-sight
direction, then set the GIN los_reversed to 1. One such case is that
if *sensor_los* represents surface incidence angles, i.e. holds the
zenith and nadir angle towards the sensor.

There is also *rte_pos_losBackwardToAltitude*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_pos", "sensor_los"},

      .in = {"sensor_pos", "sensor_los", "surface_field"},
      .gin = {"altitude", "los_is_reversed"},
      .gin_type = {"Numeric", "Index"},
      .gin_value = {std::nullopt, Index{0}},
      .gin_desc =
          {R"--(Target altitude.)--",
           R"--(Set to 1 if *rte_los* is valid for the reversed direction.)--"},

  };

  wsm_data["sensor_pos_losForwardToAltitude"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Moves *sensor_pos* and *sensor_los* forward to the target altitude.

The method gives the *sensor_pos* and *sensor_los* at the target altitude
when forward-propagating the original *sensor_pos* and *sensor_los*
geometrically.

The WSM *IntersectionGeometricAltitude* performs the same operation
but allows to store the new pos and los as other variables. There is
also *rte_pos_losForwardToAltitude*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_pos", "sensor_los"},

      .in = {"sensor_pos", "sensor_los", "surface_field"},
      .gin = {"altitude"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Target altitude.)--"},

  };

  wsm_data["sensor_responseAntenna"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes response of the antenna.

The function returns the sensor response matrix after the antenna
characteristics have been included.

The function handles "multi-beam" cases where the polarisation
coordinate system is the same for all beams.

See *antenna_dim*, *antenna_dlos* and *antenna_response* for
details on how to specify the antenna response.

The text below refers to *mblock_dlos* despite it is not an input
to the method. The method instead uses *sensor_response_dlos_grid*
but the values in this WSV are likely coming from *mblock_dlos*.

One dimensional antenna patterns are handled as other response
functions. That is, both antenna response and radiances are treated
as piece-wise linear functions, and the pencil beam calculations
must cover the full sensor response (i.e. *mblock_dlos* shall be
sufficiently broad).

There exist different options for two dimensional antenna patterns.
(If 2D, the GIN ``option_2d`` must be set, the default setting results
in an error). A normalisation is always applied for 2D antennas.

"interp_response"
Both radiances and the antenna pattern are treated as step-wise
constant functions. The antenna pattern is interpolated to the
*mblock_dlos* directions. At extrapolation, the antenna response
is set to zero. This option considers GIN ``solid_angles``, that
shall be a vector with length matching the rows of *mblock_dlos*.
The values going into *sensor_response* are the interpolated antenna
values times the corresponding solid angle.

"gridded_dlos"
This option is more similar to the 1D case. The radiances are treated
as a bi-linear function, but the antenna response is treated as step-
wise constant function (in contrast to 1D). For this option
*mblock_dlos* must match a combination of zenith and azimuth
grids, and this for a particular order. If the zenith and azimuth
grids have 3 and 2 values, respectively, the order shall be::

  [(za1,aa1); (za2,aa1); (za3,aa1); (za1,aa2); (za2,aa2); (za3,aa2)]

Both these grids must be strictly increasing and as for 1D must cover
the antenna response completely.
)--",
      .author = {"Patrick Eriksson", "Mattias Ekstrom"},
      .out = {"sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_dlos_grid"},

      .in = {"sensor_response",
             "sensor_response_f",
             "sensor_response_pol",
             "sensor_response_dlos",
             "sensor_response_f_grid",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid",
             "antenna_dim",
             "antenna_dlos",
             "antenna_response",
             "sensor_norm"},
      .gin = {"option_2d", "solid_angles"},
      .gin_type = {"String", "Vector"},
      .gin_value = {String("-"), Vector{}},
      .gin_desc =
          {R"--(Calculation option for 2D antenna cases. See above for details.)--",
           R"--(The solid angle of each *mblock_dlos* direction. Only considered for 2D with "interp_response".)--"},

  };

  wsm_data["sensor_responseBackend"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes response of the backend (spectrometer).

The function returns the sensor response matrix after the backend
characteristics have been included.

See *f_backend*, *backend_channel_response* and *sensor_norm* for
details on how to specify the backend response.
)--",
      .author = {"Mattias Ekstrom", "Patrick Eriksson"},
      .out = {"sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_f_grid"},

      .in = {"sensor_response",
             "sensor_response_f",
             "sensor_response_pol",
             "sensor_response_dlos",
             "sensor_response_f_grid",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid",
             "f_backend",
             "backend_channel_response",
             "sensor_norm"},

  };

  wsm_data["sensor_responseBackendFrequencySwitching"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Frequency switching for a pure SSB reciever.

This function can be used for simulation of frequency switching.
That is, when the final spectrum is the difference of two spectra
shifted in frequency. The switching is performed by the LO, but
for a pure singel sideband reciever this is most easily simulated
by instead shifting the backend, as done here.

A strightforward frequency switching is modelled (no folding)
The channel positions for the first measurement cycle are
f_backend+df1, and for the second f_backend+df2. The first
measurement cycle is given the negive weight. That is, the output
is the spectrum for cycle2 minus the spectrum for cycle1.
Output frequency grids are set to *f_backend*.

Use *sensor_responseFrequencySwitching* for double sideband cases.

The method has the same general functionality as, and can replace,
*sensor_responseBackend*.
)--",
          .author = {"Patrick Eriksson"},
          .out = {"sensor_response",
                  "sensor_response_f",
                  "sensor_response_pol",
                  "sensor_response_dlos",
                  "sensor_response_f_grid"},

          .in = {"sensor_response",
                 "sensor_response_f",
                 "sensor_response_pol",
                 "sensor_response_dlos",
                 "sensor_response_f_grid",
                 "sensor_response_pol_grid",
                 "sensor_response_dlos_grid",
                 "f_backend",
                 "backend_channel_response",
                 "sensor_norm"},
          .gin = {"df1", "df2"},
          .gin_type = {"Numeric", "Numeric"},
          .gin_value = {std::nullopt, std::nullopt},
          .gin_desc = {R"--(Frequency throw for cycle1.)--",
                       R"--(Frequency throw for cycle2.)--"},

      };

  wsm_data["sensor_responseBeamSwitching"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Simulation of "beam switching".

The measurement procedure is based on taking the difference between
two spectra measured in different directions, and the calculation
set-up must treat exactly two observation directions.

The returned spectrum is y = w1*y + w2*y2, where y1 and w1 are the
spectrum and weight for the first direction, respectively (y2 and
(w2 defined correspondingly for the second direction).

Zenith and azimuth angles after beam switching are set to the
values of the second direction.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_dlos_grid"},

      .in = {"sensor_response",
             "sensor_response_f",
             "sensor_response_pol",
             "sensor_response_dlos",
             "sensor_response_f_grid",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid"},
      .gin = {"w1", "w2"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {Numeric{-1}, Numeric{1}},
      .gin_desc = {R"--(Weight for values from first viewing direction.)--",
                   R"--(Weight for values from second viewing direction.)--"},

  };

  wsm_data["sensor_responseFillFgrid"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Polynomial frequency interpolation of spectra.

The sensor response methods treat the spectra to be piece-wise linear
functions. This method is a workaround for making methods handling
the spectra in a more elaborate way: it generates spectra on a more
dense grid by polynomial interpolation. The interpolation is not
done explicitly, it is incorporated into *sensor_response*.

This method should in general increase the calculation accuracy for
a given *f_grid*. However, the selection of (original) grid points
becomes more sensitive when using this method. A poor choice of grid
points can result in a decreased accuracy, or generation of negative
radiances. Test calculations indicated that the error easily can
increase with this method close the edge of *f_grid*, and it could
be wise to make *f_grid* a bit wider than actually necessary to avoid
this effect

The method shall be inserted before the antenna stage. That is, this
method shall normally be called directly after *sensor_responseInit*.

Between each neighbouring points of *f_grid*, this method adds
``nfill`` grid points. The polynomial order of the interpolation is
``polyorder``.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_f_grid"},

      .in = {"sensor_response",
             "sensor_response_f",
             "sensor_response_pol",
             "sensor_response_dlos",
             "sensor_response_f_grid",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid"},
      .gin = {"polyorder", "nfill"},
      .gin_type = {"Index", "Index"},
      .gin_value = {Index{3}, Index{2}},
      .gin_desc = {R"--(Polynomial order of interpolation)--",
                   R"--(Number of points to insert in each gap of f_grid)--"},

  };

  wsm_data["sensor_responseFrequencySwitching"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Simulation of "frequency switching".

A general method for frequency switching. The WSM
*sensor_responseBackendFrequencySwitching* gives a description of
this observation technique, and is also a more straightforward
method for pure singel sideband cases.

It is here assume that *sensor_responseMultiMixerBackend* has been
used to calculate the spectrum for two LO positions. This method
calculates the difference between these two spectra, where the
second spectrum gets weight 1 and the first weight -1 (as in
*sensor_responseBackendFrequencySwitching*).

Output frequency grids are taken from the second spectrum.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_f_grid"},

      .in = {"sensor_response",
             "sensor_response_f",
             "sensor_response_pol",
             "sensor_response_dlos",
             "sensor_response_f_grid",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid"},

  };

  wsm_data["sensor_responseGenericAMSU"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Simplified sensor setup for an AMSU-type instrument.

This function is derived from 'sensor_responseSimpleAMSU' 
but is more generalized since the number of passbands in each 
can be in the range from 1 to 4 - in order to correctly simulate
AMSU-A type sensors 

This method allows quick and simple definition of AMSU-type
sensors. Assumptions:

1. Pencil beam antenna.
2. 1-4 Passband/sidebands per channel.
3. Sideband mode "upper"
4. The channel response is rectangular.

Under these assumptions the only inputs needed are the LO positions,
the offsets from the LO, and the IF bandwidths. They are provided
in sensor_description_amsu.
)--",
      .author = {"Oscar Isoz"},
      .out = {"f_grid",
              "antenna_dim",
              "mblock_dlos",
              "sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_f_grid",
              "sensor_response_pol_grid",
              "sensor_response_dlos_grid",
              "sensor_norm"},

      .in = {"sensor_description_amsu"},
      .gin = {"spacing"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{.1e9}},
      .gin_desc = {R"--(Desired grid spacing in Hz.)--"},

  };

  wsm_data["sensor_responseIF2RF"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Converts sensor response variables from IF to RF.

The function converts intermediate frequencies (IF) in
*sensor_response_f* and *sensor_response_f_grid* to radio
frequencies (RF). This conversion is needed if the frequency
translation of a mixer is included and the position of backend
channels are specified in RF.

A direct frequency conversion is performed. Values are not
sorted in any way.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_response_f", "sensor_response_f_grid"},

      .in = {"sensor_response_f",
             "sensor_response_f_grid",
             "lo",
             "sideband_mode"},

  };

  wsm_data["sensor_responseInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initialises the variables summarising the sensor response.

This method sets the variables to match monochromatic pencil beam
calculations, to be further modified by inclusion of sensor
characteristics. Use *sensorOff* if pure monochromatic pencil
beam calculations shall be performed.

The variables are set as follows:

- sensor_response: Identity matrix, with size matching *f_grid*, ``stokes_dim`` and *mblock_dlos*.
- sensor_response_f: Repeated values of *f_grid*.
- sensor_response_pol: Data matching ``stokes_dim``.
- sensor_response_dlos: Repeated values of *mblock_dlos*.
- sensor_response_f_grid: Equal to *f_grid*.
- sensor_response_pol_grid: Set to 1:``stokes_dim``.
- sensor_response_dlos_grid: Equal to *mblock_dlos*.
)--",
      .author = {"Mattias Ekstrom", "Patrick Eriksson"},
      .out = {"sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_f_grid",
              "sensor_response_pol_grid",
              "sensor_response_dlos_grid"},

      .in = {"f_grid", "mblock_dlos", "antenna_dim", "sensor_norm"},

  };

  wsm_data["sensor_responseMetMM"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sensor setup for meteorological millimeter instruments.

This method is handy if you are simulating a passband-type instrument,
consisting of a few discrete channels.

For flexibility, the Met-MM system is seperated in two calculation
steps. To fully use the system, create *f_grid* (and some associated
variables) by *f_gridMetMM* before calling this method. However, it is
possible to use this method with any *f_grid*, as long as matching
*f_backend*, *channel2fgrid_indexes* and *channel2fgrid_weights*
are provided.

Each scan sequence is treated as a measurement block. *sensor_pos* is
set in the standard way. The number of rows in *sensor_pos* determines the
number of scan sequences that will be simulated. On the other hand,
*sensor_los* is handled in a special way. All zenith angles must be set
to 180 deg. For 3D, the given azimuth angles are taken as the direction
of scanning, where the azimuth angle is defined with respect to North
in standard manner. For example, if the scanning happens to move from
SW to NE, the azimuth angle should be set to 45 deg. The angles of the
scanning sequence are taken from *antenna_dlos*. This WSV is here only
allowed to have a single column, holding relative zenith angles. For
3D, the azimuth angles in *antenna_dlos* are hard-coded to zero. As
zenith angles in *sensor_los* are locked to 180 deg, *antenna_dlos*
effectively holds the nadir angles. These angles can be both positive or
negative, where the recommended choice is to operate with negative
to end up with final zenith angles between 0 and 180 deg.

The method does not support 2D atmospheres (across-track scanning is
inconsistent with 2D). For simpler switching between 1D and 3D,
the argument ``mirror_dza`` is at hand. It can only be used for 3D.
If set to true, the zenith angles in *antenna_dlos* are mapped
to also cover the other side of the swath and the simulations will
cover both sides of the swath.
)--",
      .author = {"Oliver Lemke", "Patrick Eriksson"},
      .out = {"antenna_dim",
              "mblock_dlos",
              "sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_f_grid",
              "sensor_response_pol_grid",
              "sensor_response_dlos_grid",
              "sensor_norm"},

      .in = {"f_grid",
             "f_backend",
             "channel2fgrid_indexes",
             "channel2fgrid_weights",
             "iy_unit",
             "antenna_dlos",
             "met_mm_polarisation",
             "met_mm_antenna"},
      .gin = {"use_antenna", "mirror_dza"},
      .gin_type = {"Index", "Index"},
      .gin_value = {Index{0}, Index{0}},
      .gin_desc =
          {R"--(Flag to enable (1) or disable (0) antenna.)--",
           R"--(Flag to include second part of swath (only 3D, see above).)--"},

  };

  wsm_data["sensor_responseMixer"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes response of the mixer of a heterodyne system.

The function returns the sensor response matrix after the mixer
characteristics have been included. Frequency variables are
converted from radio frequency (RF) to intermediate frequency (IF).
The returned frequency grid covers the range [0,max_if], where
max_if is the highest IF covered by the sideband response grid.

See *lo* and *sideband_response* for details on how to specify the
mixer response
)--",
      .author = {"Mattias Ekstrom", "Patrick Eriksson"},
      .out = {"sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_f_grid"},

      .in = {"sensor_response",
             "sensor_response_f",
             "sensor_response_pol",
             "sensor_response_dlos",
             "sensor_response_f_grid",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid",
             "lo",
             "sideband_response",
             "sensor_norm"},

  };

  wsm_data["sensor_responseMixerBackendPrecalcWeights"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Includes pre-calculated response covering mixer and backend.

This method acts similar to *sensor_responseBackend*, but uses
pre-calculated weights. These weights can also include the effect
of mixer and sideband filtering.

As usual, *f_backend* gives the frequency of the channels. This WSM
has no direct influence on the result, but at least representative
values must be set.

The frequency response is defined using *channel2fgrid_indexes* and
*channel2fgrid_weights*.

Both *channel2fgrid_indexes* and *channel2fgrid_weights* are assumed
to be common for all viewing directions.
)--",
          .author = {"Patrick Eriksson"},
          .out = {"sensor_response",
                  "sensor_response_f",
                  "sensor_response_pol",
                  "sensor_response_dlos",
                  "sensor_response_f_grid"},

          .in = {"sensor_response",
                 "sensor_response_f",
                 "sensor_response_pol",
                 "sensor_response_dlos",
                 "sensor_response_f_grid",
                 "sensor_response_pol_grid",
                 "sensor_response_dlos_grid",
                 "f_backend",
                 "channel2fgrid_indexes",
                 "channel2fgrid_weights"},

      };

  wsm_data["sensor_responseMultiMixerBackend"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Handles mixer and backend parts for an instrument having multiple
mixer chains.

The WSMs *sensor_responseMixer*, *sensor_responseIF2RF* and
*sensor_responseBackend* are called for each mixer chain, and a
complete *sensor_response* is assembled. The instrument responses
are described by *lo_multi*, *sideband_response_multi*,
*sideband_mode_multi*, *f_backend_multi* and
*backend_channel_response_multi*. All these WSVs must have same
vector or array length. As *sensor_responseIF2RF* is called,
*f_backend_multi* must hold RF (not IF) and output frequencies
will be in absolute frequency (RF).
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_f_grid"},

      .in = {"sensor_response",
             "sensor_response_f",
             "sensor_response_pol",
             "sensor_response_dlos",
             "sensor_response_f_grid",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid",
             "lo_multi",
             "sideband_response_multi",
             "sideband_mode_multi",
             "f_backend_multi",
             "backend_channel_response_multi",
             "sensor_norm"},

  };

  wsm_data["sensor_responsePolarisation"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Extraction of non-default polarisation components.

The default is to output the Stokes elements I, Q, U and V (up to
``stokes_dim``). This method allows to change the "polarisation" of
the output. Polarisation components to be extracted are selected by
*instrument_pol*. This method can be applied at any step of the sensor
matrix set-up.

The method can only be applied on data for I, Q, U and V. The value
of ``stokes_dim`` must be sufficiently large for the selected
components. For example, I+45 requires that ``stokes_dim`` is at
least 3. 

See *instrument_pol* for coding of polarisation states.

Note that the state of *iy_unit* is considered. This WSV must give
the actual unit of the data. This as, the extraction of components
is slightly different if data are radiances or brightness
temperatures.  In practise this means that *iy_unit* (as to be
applied inside *iy_main_agenda*) must be set before calling this
method.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_pol_grid"},

      .in = {"sensor_response",
             "sensor_response_f",
             "sensor_response_pol",
             "sensor_response_dlos",
             "sensor_response_f_grid",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid",
             "iy_unit",
             "instrument_pol"},

  };

  wsm_data["sensor_responseSimpleAMSU"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Simplified sensor setup for an AMSU-type instrument.

This method allows quick and simple definition of AMSU-type
sensors. Assumptions:

1. Pencil beam antenna.
2. Double sideband receivers.
3. Sideband mode "upper"
4. The channel response is rectangular.

Under these assumptions the only inputs needed are the LO positions,
the offsets from the LO, and the IF bandwidths. They are provieded
in sensor_description_amsu.
)--",
      .author = {"Stefan Buehler"},
      .out = {"f_grid",
              "antenna_dim",
              "mblock_dlos",
              "sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_f_grid",
              "sensor_response_pol_grid",
              "sensor_response_dlos_grid",
              "sensor_norm"},

      .in = {"sensor_description_amsu"},
      .gin = {"spacing"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{.1e9}},
      .gin_desc = {R"--(Desired grid spacing in Hz.)--"},

  };

  wsm_data["sensor_responseStokesRotation"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Includes a rotation of the Stokes H and V directions.

The method applies the rotations implied by *stokes_rotation*.
See the description of that WSV for details.

This method does not change the size of *sensor_response*, and
the auxiliary variables (sensor_response_f etc.) are not changed.

To apply the method, ``stokes_dim`` must be >= 3. The complete effect
of the rotation can not be determibed with lower ``stokes_dim``.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_response"},

      .in = {"sensor_response",
             "sensor_response_f_grid",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid",
             "stokes_rotation"},

  };

  wsm_data["sensor_responseWMRF"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds WMRF weights to sensor response.

This method adds a spectrometer response that has been calculated
with the weighted mean of representative frequencies (WMRF) method. It
consists of a set of selected frequencies, and associated weights.
)--",
      .author =
          {"Stefan Buehler, based on Patrick Erikssons sensor_responseBackend"},
      .out = {"sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_f_grid"},

      .in = {"sensor_response",
             "sensor_response_f",
             "sensor_response_pol",
             "sensor_response_dlos",
             "sensor_response_f_grid",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid",
             "wmrf_weights",
             "f_backend"},

  };

  wsm_data["sensor_response_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *sensor_response_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--",
      .author = {"Richard Larsson"},
      .out = {"sensor_response_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["sparse_f_gridFromFrequencyGrid"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Outputs the sparse frequency grid in *propmat_clearskyAddLines*
)--",
      .author = {"Richard Larsson"},

      .gout = {"sparse_f_grid"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(A sparse frequency grid.)--"},
      .in = {"f_grid"},
      .gin = {"sparse_df", "speedup_option"},
      .gin_type = {"Numeric", "String"},
      .gin_value = {Numeric{0}, String("None")},
      .gin_desc = {R"--(The grid sparse separation)--", R"--(Speedup logic)--"},

  };

  wsm_data["spectral_irradiance_fieldDisort"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Interface to the DISORT scattering solver (by Stamnes et al.).
for running flux (irradiance) calculations

It provides the irradiance field from a scalar
1D scattering solution assuming a plane-parallel atmosphere (flat
Earth). Only totally randomly oriented particles are allowed.
Refraction is not taken into account. Only Lambertian surface
reflection is handled.

``nstreams`` is the number of polar angles taken into account
internally in the scattering solution and for the angular integration.
``nstreams`` determines the angular resolution, hence the accuracy,
of the scattering solution. The more anisotropic the bulk scattering
matrix, the more streams are required. The computational burden
increases approximately linearly with ``nstreams``. The default value
(6) is sufficient for most flux calculations.

Some auxiliary quantities can be obtained. Auxiliary
quantities are selected by *disort_aux_vars* and returned by *disort_aux*.
Valid choices for auxiliary data are:

- ``"Layer optical thickness"``: Matrix [f_grid, size of p_grid - 1] layer optical thickness.
- ``"Single scattering albedo"``: Matrix [f_grid, size of p_grid - 1] layer single scattering albedo.
- ``"Direct downward spectral irradiance"``: Matrix [f_grid, p_grid]. Direct downward spectral irradiance. Zero, if no sun is present. 
- ``"dFdtau"``: Matrix [f_grid, p_grid]. Flux divergence in optical thickness space.
)--",
      .author = {"Manfred Brath"},
      .out = {"spectral_irradiance_field", "disort_aux"},

      .in = {"atmfields_checked",
             "atmgeom_checked",
             "scat_data_checked",
             "propmat_clearsky_agenda",
             "gas_scattering_agenda",
             "pnd_field",
             "atm_field",
             "surface_field",
             "lat_true",
             "lon_true",
             "abs_species",
             "scat_data",
             "suns",
             "f_grid",
             "surface_skin_t",
             "surface_scalar_reflectivity",
             "gas_scattering_do",
             "suns_do",
             "disort_aux_vars"},
      .gin = {"nstreams",
              "Npfct",
              "only_tro",
              "quiet",
              "emission",
              "intensity_correction"},
      .gin_type = {"Index", "Index", "Index", "Index", "Index", "Index"},
      .gin_value =
          {Index{6}, Index{181}, Index{0}, Index{0}, Index{1}, Index{1}},
      .gin_desc =
          {R"--(Number of polar angle directions (streams) in DISORT solution (must be an even number).)--",
           R"--(Number of angular grid points to calculate bulk phase function on (and derive Legendre polynomials from). If <0, the finest za_grid from scat_data will be used.)--",
           R"--(Set to 1 if the scattering data is just of TRO type. Has effect only if Npfct > 3 or Npfct<0, but then leads to much faster calculations.)--",
           R"--(Silence C Disort warnings.)--",
           R"--(Enables blackbody emission. Set to zero, if no  Emission e. g. like in visible regime for earth is needed)--",
           R"--(Enables intensity correction. Importantant for low number of  streams. Set to zero, if problems encounter or using a high number  of streams (>30))--"},
      .pass_workspace = true,

  };

  wsm_data["spectral_irradiance_fieldFromSpectralRadianceField"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Calculates the spectral irradiance from *spectral_radiance_field*.

The *spectral_radiance_field* is integrated over the angular grids
according to the grids set by *AngularGridsSetFluxCalc*.
See *AngularGridsSetFluxCalc* to set *za_grid*, *aa_grid*, and 
*za_grid_weights*.
)--",
          .author = {"Manfred Brath"},
          .out = {"spectral_irradiance_field"},

          .in = {"spectral_radiance_field",
                 "za_grid",
                 "aa_grid",
                 "za_grid_weights"},

      };

  wsm_data["spectral_radiance_fieldClearskyPlaneParallel"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Clear-sky radiance field of a plane parallel atmosphere.

The method assumes a 1D flat planet. Radiances along each direction
given by *za_grid* are calculated using ``ppathPlaneParallel``
and ``iyEmissionStandard``.

Surface properties are defined by *iy_surface_agenda*, i.e. there is no
restriction to e.g. specular surfaces.

Note that the variable *ppath_lmax* is considered, and that it can be
critical for the accuracy for zenith angles close to 90 degrees. That
is, using ppath_lmax=-1 is not recommended for this function.

Information on transmittance is also provided by the GOUT ``trans_field``.
For up-welling radiation (scat_za > 90), this variable holds the
transmittance to space, for considered position and propagation direction.
For down-welling radiation, ``trans_field`` holds instead the transmittance
down to the surface.
)--",
          .author = {"Patrick Eriksson"},
          .out = {"spectral_radiance_field"},
          .gout = {"trans_field"},
          .gout_type = {"Tensor3"},
          .gout_desc =
              {R"--(Dimensions: [f_grid,p_grid,za_grid]. See further above.)--"},
          .in = {"propmat_clearsky_agenda",
                 "water_p_eq_agenda",
                 "iy_space_agenda",
                 "iy_surface_agenda",
                 "iy_cloudbox_agenda",
                 "f_grid",
                 "abs_species",
                 "atm_field",
                 "surface_field",
                 "ppath_lmax",
                 "rte_alonglos_v",
                 "rt_integration_option",
                 "za_grid"},
          .gin = {"use_parallel_za"},
          .gin_type = {"Index"},
          .gin_value = {Index{1}},
          .gin_desc =
              {R"--(Flag to select parallelization over zenith angles.)--"},
          .pass_workspace = true,

      };

  wsm_data["spectral_radiance_fieldDisortClearsky"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Interface to the DISORT scattering solver (by Stamnes et al.).
for running clear-sky cases.

The method runs DISORT with *pnd_field* set to zero.

Note that this version returns *spectral_radiance_field*, i.e.
the solution for the full atmosphere. The standard *cloudbox_fieldDisort*
only returns the field inside the cloudbox.

Some auxiliary quantities can be obtained. Auxiliary
quantities are selected by *disort_aux_vars* and returned by *disort_aux*.
Valid choices for auxiliary data are:

- ``"Layer optical thickness"``: Matrix [f_grid, size of p_grid - 1] layer optical thickness.
- ``"Single scattering albedo"``: Matrix [f_grid, size of p_grid - 1] layer single scattering albedo.
- ``"Direct beam"``: Matrix [f_grid, p_grid]. Level direct spectral radiance. Zero, if no sun is present 
)--",
      .author = {"Patrick Eriksson", "Manfred Brath"},
      .out = {"spectral_radiance_field", "disort_aux"},

      .in = {"atmfields_checked",
             "atmgeom_checked",
             "propmat_clearsky_agenda",
             "gas_scattering_agenda",
             "atm_field",
             "surface_field",
             "lat_true",
             "lon_true",
             "abs_species",
             "suns",
             "f_grid",
             "za_grid",
             "aa_grid",
             "surface_skin_t",
             "surface_scalar_reflectivity",
             "gas_scattering_do",
             "suns_do",
             "disort_aux_vars"},
      .gin = {"nstreams", "quiet", "emission", "intensity_correction"},
      .gin_type = {"Index", "Index", "Index", "Index"},
      .gin_value = {Index{8}, Index{0}, Index{1}, Index{1}},
      .gin_desc =
          {R"--(Number of polar angle directions (streams) in DISORT solution (must be an even number).)--",
           R"--(Silence C Disort warnings.)--",
           R"--(Enables blackbody emission. Set to zero, if no  Emission e. g. like in visible regime for earth is needed)--",
           R"--(Enables intensity correction. Importantant for low number of  streams. Set to zero, if problems encounter or using a high number  of streams (>30))--"},
      .pass_workspace = true,

  };

  wsm_data["spectral_radiance_fieldExpandCloudboxField"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Uses and expands *cloudbox_field* to set *spectral_radiance_field*.

The method demands that *cloudbox_field* starts at the first pressure
level (i.e. cloudbox_limits[0] is 0). The method copies *cloudbox_field*
to fill *spectral_radiance_field* up to the top of the cloudbox.

To fill the remaning part of *spectral_radiance_field*, clear-sky
calculations are performed largely in the same maner as done by
*spectral_radiance_fieldClearskyPlaneParallel*. That is, clear-sky
calculations are done for the upper part of the atmosphere, assuming
a flat planet.

Note that the cloud box constitutes the lower boundary for the later
calculations, and *iy_cloudbox_agenda* must be set to perform an
interpolation of the cloudbox field.
)--",
          .author = {"Patrick Eriksson"},
          .out = {"spectral_radiance_field"},

          .in = {"propmat_clearsky_agenda",
                 "water_p_eq_agenda",
                 "iy_space_agenda",
                 "iy_surface_agenda",
                 "iy_cloudbox_agenda",
                 "f_grid",
                 "abs_species",
                 "atm_field",
                 "surface_field",
                 "cloudbox_on",
                 "cloudbox_limits",
                 "cloudbox_field",
                 "ppath_lmax",
                 "rte_alonglos_v",
                 "rt_integration_option",
                 "za_grid"},
          .gin = {"use_parallel_za"},
          .gin_type = {"Index"},
          .gin_value = {Index{0}},
          .gin_desc =
              {R"--(Flag to select parallelization over zenith angles.)--"},
          .pass_workspace = true,

      };

  wsm_data["spectral_radiance_fieldPlaneParallelSpectralRadianceOperator"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Create a *spectral_radiance_field*

This is an experimental solution.
)--",
          .author = {"Richard Larsson"},
          .out = {"spectral_radiance_field"},

          .in = {"spectral_radiance_profile_operator", "f_grid", "za_grid"},

      };

  wsm_data["specular_losCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates the specular direction of surface reflections.

The default is to consider surface topography when calculating the
specular direction. That is, the variation of ``surface_elevation``
is allowed to affect the angles of *specular_los*. This impact can
be deactivated by setting ``ignore_topography`` to 1. In this case,
the zenith angle of the specular direction is simply 180-rtp_los[0]
and the azimuth angle is the same as the one in *rtp_los*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"specular_los"},

      .in = {"surface_field", "rtp_pos", "rtp_los"},
      .gin = {"ignore_topography"},
      .gin_type = {"Index"},
      .gin_value = {Index{0}},
      .gin_desc =
          {R"--(Flag to control if surface slope is considered or not.)--"},

  };

  wsm_data["specular_losCalcOldNoTopography"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculates the specular direction of surface reflections for horisontal
surfaces.

In contrast to ``specular_losCalcOld``, this method ignores the topography
implied by ``z_surface``. That is, any slope of the surface is ignored.

The typical application of this WSM should be water surfaces (lakes and
oceans).
)--",
      .author = {"Patrick Eriksson"},
      .out = {"specular_los", "surface_normal"},

      .in = {"rtp_los"},

  };

  wsm_data["spt_calc_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *spt_calc_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--",
      .author = {"Richard Larsson"},
      .out = {"spt_calc_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["sunsAddSingleBlackbody"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds a single blackbody to *suns*

Important note:
For a Sol-like sun there are huge differences in the UV-range 
between the actual sun spectrum and the blackbody spectrumwith the effective temperature of the sun. The blackbody sun"strongly overestimates the UV radiation.
)--",
      .author = {"Jon Petersen"},
      .out = {"suns", "suns_do"},

      .in = {"suns", "f_grid"},
      .gin = {"radius", "distance", "temperature", "latitude", "longitude"},
      .gin_type = {"Numeric", "Numeric", "Numeric", "Numeric", "Numeric"},
      .gin_value = {Numeric{6.963242e8},
                    Numeric{1.495978707e11},
                    Numeric{5772},
                    Numeric{0},
                    Numeric{0}},
      .gin_desc =
          {R"--(The radius of the sun in meter. Default is the radius of our sun. )--",
           R"--(The average distance between the sun and the planet in meter. Default value is set to 1 a.u. )--",
           R"--(The effective temperature of the suns photosphere in Kelvin. Default is the temperature of our sun - 5772 Kelvin )--",
           R"--(The latitude or the zenith position of the sun in the sky. )--",
           R"--(The longitude or azimuthal position of the sun in the sky. )--"},

  };

  wsm_data["sunsAddSingleFromGrid"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Extracts a sun spectrum from a field of such data and
adds it to *suns*.

The method allows to obtain the sun spectrum by
interpolation from a field of such data. 
The sun spectrum is expected to be stored as
the irradiance at the suns photosphere.

Unit:

- GriddedField2: [W m-2 Hz-1]

  - Vector *f_grid* [Hz]
  - Vector ``stokes_dim`` [1]

Dimensions: [f_grid, stokes_dim]

This method performs an interpolation onto the f_grid.
The point of *f_grid* that are outside the data frequency grid
are initialized according to planck's law of the temperature variable.
Hence, a temperature of 0 means 0s the edges of the f_grid.
)--",
      .author = {"Jon Petersen"},
      .out = {"suns", "suns_do"},

      .in = {"suns", "f_grid"},
      .gin = {"sun_spectrum_raw",
              "radius",
              "distance",
              "temperature",
              "latitude",
              "longitude",
              "description"},
      .gin_type = {"GriddedField2",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "String"},
      .gin_value = {std::nullopt,
                    Numeric{6.963242e8},
                    Numeric{1.495978707e11},
                    Numeric{-1},
                    Numeric{0},
                    Numeric{0},
                    String("Sun spectrum from Griddedfield.")},
      .gin_desc =
          {R"--(Raw data for monochromatic irradiance spectra. )--",
           R"--(The radius of the sun in meter. Default is the radius of our sun. )--",
           R"--(The average distance between the center of the sun and the  center of the planet in meter. Default value is set to 1 a.u. )--",
           R"--(The temperature of the padding if the f_grid is outside the  sun spectrum data. Choose 0 for 0 at the edges or a effective temperature for a padding using plack's law. )--",
           R"--(The latitude or the zenith position of the sun in the sky. )--",
           R"--(The longitude or azimuthal position of the sun in the sky. )--",
           R"--(The description of the sun. )--"},

  };

  wsm_data["sunsAddSingleFromGridAtLocation"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Extracts a sun spectrum measured at the given location
adds it to *suns*.

The method allows to obtain the sun spectrum by
interpolation from a field of such data. 
The sun spectrum is expected to be stored as
irradiance.
It is coverted to the irradiance at the suns photosphere.

Unit:

- GriddedField2: [W m-2 Hz-1]

  - Vector *f_grid* [Hz]
  - Vector ``stokes_dim`` [1]

Dimensions: [f_grid, stokes_dim]

This method performs an interpolation onto the f_grid.
The point of *f_grid* that are outside the data frequency grid
are initialized according to planck's law of the temperature variable.
Hence, a temperature of 0 means 0s the edges of the f_grid.
)--",
      .author = {"Jon Petersen"},
      .out = {"suns", "suns_do"},

      .in = {"suns", "f_grid", "surface_field"},
      .gin = {"sun_spectrum_raw",
              "radius",
              "distance",
              "temperature",
              "zenith",
              "azimuth",
              "description",
              "location_latitude",
              "location_longitude",
              "location_altitude"},
      .gin_type = {"GriddedField2",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "String",
                   "Numeric",
                   "Numeric",
                   "Numeric"},
      .gin_value = {std::nullopt,
                    Numeric{6.963242e8},
                    Numeric{1.495978707e11},
                    Numeric{-1},
                    Numeric{0},
                    Numeric{0},
                    String("Sun spectrum from Griddedfield."),
                    Numeric{0},
                    Numeric{0},
                    Numeric{1e5}},
      .gin_desc =
          {R"--(Raw data for monochromatic irradiance spectra. )--",
           R"--(The radius of the sun in meter. Default is the radius of our Sun. )--",
           R"--(The distance between the location and the  center of the sun in meter. Default value is set to 1 a.u. )--",
           R"--(The temperature of the padding if the f_grid is outside the  sun spectrum data. Choose 0 for 0 at the edges or a effective temperature for a padding using plack's law. )--",
           R"--(Zenith angle of the sun in the sky. )--",
           R"--(Azimuthal angle of the sun in the sky. )--",
           R"--(The description of the sun. )--",
           R"--(The latitude of the sun spectrum measurement. )--",
           R"--(The longitude of the sun spectrum measurement. )--",
           R"--(The altitude of the sun spectrum measurement. )--"},

  };

  wsm_data["sunsOff"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Turns all calculations with suns off 
)--",
      .author = {"Jon Petersen"},
      .out = {"suns_do", "suns"},

  };

  wsm_data["surfaceBlackbody"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Creates variables to mimic a blackbody surface.

This method sets up *surface_los*, *surface_rmatrix* and
*surface_emission* for *surface_rtprop_agenda*. Here, *surface_los*
and *surface_rmatrix* are set to be empty, and *surface_emission*
to hold blackbody radiation for a temperature of *surface_skin_t*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_los", "surface_rmatrix", "surface_emission"},

      .in = {"f_grid", "rtp_pos", "rtp_los", "surface_point"},

  };

  wsm_data["surfaceFastem"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Usage of FASTEM together with MC and DOIT.

The recommended way to use FASTEM is by *iySurfaceFastem*, but that
is not always possible, such as when using MC and DOIT. This is the
case as those scattering methods use *surface_rtprop_agenda*,
while *iySurfaceFastem* fits with *iy_surface_agenda*. This WSM solves
this by allowing FASTEM to be used inside *surface_rtprop_agenda*.

However, FASTEM is here used in an approximative way. For a correct
usage of FASTEM, the atmospheric transmittance shall be calculated
for the position and direction of concern, but this is not possible
together with DOIT and MC. Instead, the transmittance is an input
to the method, and must either be pre-calculated or set to a
representative value.

See *iySurfaceFastem*, for further details on the special input
arguments.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_los", "surface_rmatrix", "surface_emission"},

      .in = {"f_grid", "rtp_pos", "rtp_los", "surface_skin_t"},
      .gin = {"salinity",
              "wind_speed",
              "wind_direction",
              "transmittance",
              "fastem_version"},
      .gin_type = {"Numeric", "Numeric", "Numeric", "Vector", "Index"},
      .gin_value =
          {Numeric{0.035}, std::nullopt, Numeric{0}, std::nullopt, Index{6}},
      .gin_desc =
          {R"--(Salinity, 0-1. That is, 3% is given as 0.03.)--",
           R"--(Wind speed.)--",
           R"--(Wind direction. See futher above.)--",
           R"--(Transmittance along path of downwelling radiation. A vector with the same length as *f_grid*.)--",
           R"--(The version of FASTEM to use.)--"},

  };

  wsm_data["surfaceFlatReflectivity"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Creates variables to mimic specular reflection by a (flat) surface
where *surface_reflectivity* is specified.

Works basically as *surfaceFlatScalarReflectivity* but is more
general as vector radiative transfer is more properly handled. See
the ARTS theory document (ATD) for details around how
*surface_emission* is determined. In the nomenclature of ATD,
*surface_reflectivity* gives R.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_los", "surface_rmatrix", "surface_emission"},

      .in = {"f_grid",
             "rtp_pos",
             "rtp_los",
             "specular_los",
             "surface_skin_t",
             "surface_reflectivity"},

  };

  wsm_data["surfaceFlatRefractiveIndex"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Creates variables to mimic specular reflection by a (flat) surface
where the complex refractive index is specified.

The dielectric properties of the surface are described by
*surface_complex_refr_index*. The Fresnel equations are used to
calculate amplitude reflection coefficients. The method can thus
result in that the reflection properties differ between frequencies
and polarisations.

Local thermodynamic equilibrium is assumed, which corresponds to
that the reflection and emission coefficients add up to 1.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_los", "surface_rmatrix", "surface_emission"},

      .in = {"f_grid",
             "rtp_pos",
             "rtp_los",
             "specular_los",
             "surface_skin_t",
             "surface_complex_refr_index"},

  };

  wsm_data["surfaceFlatRvRh"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Creates variables to mimic specular reflection by a (flat) surface
where *surface_rv_rh* is specified.

This method assumes that the reflection at vertical and horizontal
polarisation differs. As power reflection coefficients are provided
there is no information at hand on phase shifts between polarisations,
and they are simply assumed to be zero. These assumptions result in
that *surface_emission* is set to zero for positions corresponding to
U and V, and that all diagonal elementsof  *surface_rmatrix* are equal
(the mean of rv and rh). Further, all off-diagonal elements of
*surface_rmatrix* are all zero except for (0,1) and (1,0).
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_los", "surface_rmatrix", "surface_emission"},

      .in = {"f_grid",
             "rtp_pos",
             "rtp_los",
             "specular_los",
             "surface_skin_t",
             "surface_rv_rh"},

  };

  wsm_data["surfaceFlatScalarReflectivity"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Creates variables to mimic specular reflection by a (flat) surface
where *surface_scalar_reflectivity* is specified.

This method assumes that the reflection at vertical and horizontal
polarisation is identical. This assumption includes that there is no
phase shift between polarisations. These assumptions result in that
*surface_emission* is set to zero for positions corresponding to Q,
U and V, and that *surface_rmatrix* becomes a diagonal matrix (with
all elements on the diagonal equal to the specified reflectivity).
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_los", "surface_rmatrix", "surface_emission"},

      .in = {"f_grid",
             "rtp_pos",
             "rtp_los",
             "specular_los",
             "surface_skin_t",
             "surface_scalar_reflectivity"},

  };

  wsm_data["surfaceLambertianSimple"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Creates variables to mimic a Lambertian surface.

A Lambertian surface can be characterised solely by its
reflectivity, here taken from *surface_scalar_reflectivity*.

The down-welling radiation field is estimated by making calculations
for ``lambertian_nza`` directions. The range of zenith angles ([0,90])
is divided in an equidistant manner for 1D. For 2D and 3D see below.
The values for *surface_rmatrix* are assuming a constant radiance
over each zenith angle range. See AUG.

Default is to select the zenith angles for *sensor_los* to be placed
centrally in the grid ranges. For example, if ``lambertian_nza`` is set
to 9, down-welling radiation will be calculated for zenith angles = 
5, 15, ..., 85. The position of these angles can be shifted by
``za_pos``. This variable specifies the fractional distance inside the
ranges. For example, a ``za_pos`` of 0.7 (np still 9) gives the angles
7, 17, ..., 87.

Only upper-left diagonal element of the *surface_rmatrix* is
non-zero. That is, the upwelling radiation is always unpolarised.

Local thermodynamic equilibrium is assumed, which corresponds to
that the reflection and emission coefficients "add up to 1".

For 2D and 3D, the down-welling directions are placed along the
the viewing direction, e.g. for 3D the azimuth angle is kept constant.
In 2D and 3D surface topography can exist, and to avoid getting views
going directly into the surface, angels are not distributed over 90 deg,
but 90-abs(surface_normal[0]).
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_los", "surface_rmatrix", "surface_emission"},

      .in = {"f_grid",
             "rtp_pos",
             "rtp_los",
             "surface_normal",
             "surface_skin_t",
             "surface_scalar_reflectivity"},
      .gin = {"lambertian_nza", "za_pos"},
      .gin_type = {"Index", "Numeric"},
      .gin_value = {Index{9}, Numeric{0.5}},
      .gin_desc =
          {R"--(Number of downwelling streams.)--",
           R"--(Position of angle in *surface_los* inside ranges of zenith angle grid. See above.)--"},

  };

  wsm_data["surfaceMapToLinearPolarisation"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Convert surface RT properties to a linear polarisation.

The properties converted are *surface_emission* and *surface_rmatrix*.

This method allows to set the surface properties to match a specific
linear polarisation for scalar calculation (stokes_dim=1). If you want
this, you have to call the method(s) setting up the surface RT
properties with ``stokes_dim`` set to 2, 3 or 4. This Stokes dimension
is below called local_stokes_dim.

The polarisation to apply is selected by ``pol_angle``. This angle is
defined as in *sensor_pol* (i.e. 0 and 90 equal V and H, respectively).

If local_stokes_dim was set to 2, *surface_rmatrix* is assumed to have
the structure:

.. math:: 
    \begin{array}{cc} (rv+rh)/2 & (rv-rh)/2 \\ (rv-rh)/2 & (rv+rh)/2 \end{array}

while if local_stokes_dim was set to 3 or 4, the mapping involves
several transformation matrices. The later case covers also couplings
between V/H and +-45 deg, and the mapping is described in the ARTS
theory guide, in section "Rotated modified Stokes vector".

In general it should suffice to set local_stokes_dim to 2, that gives
slightly faster calculations. A local_stokes_dim of 3 handles any case
correctly.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_emission", "surface_rmatrix"},

      .in = {"surface_emission", "surface_rmatrix"},
      .gin = {"pol_angle"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Polarisation angle, see above.)--"},

  };

  wsm_data["surfaceTessem"] = WorkspaceMethodInternalRecord{
      .desc = R"--(TESSEM sea surface microwave emissivity parametrization.

This method computes surface emissivity and reflectivity matrices for
ocean surfaces using the TESSEM emissivity model: Prigent, C., et al.
Sea‐surface emissivity parametrization from microwaves to millimetre
waves, QJRMS, 2017, 143.702: 596-605.

The validity range of the parametrization of is 10 to 700 GHz, but for
some extra flexibility frequencies between 5 and 900 GHz are accepted.
The accepted temperaute range for *surface_skin_t* is [260.0 K, 373.0 K]

The model itself is represented by the neural networks in
*tessem_neth* and *tessem_netv*.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"surface_los", "surface_rmatrix", "surface_emission"},

      .in = {"f_grid",
             "rtp_pos",
             "rtp_los",
             "surface_skin_t",
             "tessem_neth",
             "tessem_netv"},
      .gin = {"salinity", "wind_speed"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {Numeric{0.035}, std::nullopt},
      .gin_desc = {R"--(Salinity, 0-1. That is, 3% is given as 0.03.)--",
                   R"--(Wind speed.)--"},

  };

  wsm_data["surface_fieldEarth"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Earth reference ellipsoids.

The reference ellipsoid (``refellipsoid``) is set to model the Earth,
following different models. The options are:

- "Sphere":
    A spherical Earth. The radius is set following
    the value set for the Earth radius in constants.cc.
- "WGS84":
    The reference ellipsoid used by the GPS system.
    Should be the standard choice for a non-spherical Earth.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_field"},

      .gin = {"model"},
      .gin_type = {"String"},
      .gin_value = {String("Sphere")},
      .gin_desc = {R"--(Model ellipsoid to use. Options listed above.)--"},

  };

  wsm_data["surface_fieldEuropa"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Io reference ellipsoids.

The reference ellipsoid (``refellipsoid``) is set to model Io,
folowing different models. The options are:

- "Sphere": A spherical planetesimal. The radius is taken from report of the IAU/IAG Working Group.
)--",
      .author = {"Richard Larsson"},
      .out = {"surface_field"},

      .gin = {"model"},
      .gin_type = {"String"},
      .gin_value = {String("Sphere")},
      .gin_desc = {R"--(Model ellipsoid to use. Options listed above.)--"},

  };

  wsm_data["surface_fieldGanymede"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Ganymede reference ellipsoids.

From Wikipedia
)--",
      .author = {"Takayoshi Yamada"},
      .out = {"surface_field"},

      .gin = {"model"},
      .gin_type = {"String"},
      .gin_value = {String("Sphere")},
      .gin_desc = {R"--(Model ellipsoid to use. Options listed above.)--"},

  };

  wsm_data["surface_fieldInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Manual setting of the reference ellipsoid.

The two values of ``refellipsoid`` can here be set manually. The two
arguments correspond directly to first and second element of
``refellipsoid``.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_field"},

      .gin = {"a", "b"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Average or equatorial radius.)--",
                   R"--(Average or polar radius.)--"},

  };

  wsm_data["surface_fieldIo"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Io reference ellipsoids.

The reference ellipsoid (``refellipsoid``) is set to model Io,
folowing different models. The options are:

- "Sphere": A spherical planetesimal. The radius is taken from report of the IAU/IAG Working Group.
)--",
      .author = {"Richard Larsson"},
      .out = {"surface_field"},

      .gin = {"model"},
      .gin_type = {"String"},
      .gin_value = {String("Sphere")},
      .gin_desc = {R"--(Model ellipsoid to use. Options listed above.)--"},

  };

  wsm_data["surface_fieldJupiter"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Jupiter reference ellipsoids.

The reference ellipsoid (``refellipsoid``) is set to model Jupiter,
folowing different models. The options are:

- "Sphere": A spherical planet. The radius is taken from a report of the IAU/IAG Working Group.
- "Ellipsoid": A reference ellipsoid with parameters taken from a report of the IAU/IAG Working Group.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_field"},

      .gin = {"model"},
      .gin_type = {"String"},
      .gin_value = {String("Sphere")},
      .gin_desc = {R"--(Model ellipsoid to use. Options listed above.)--"},

  };

  wsm_data["surface_fieldMars"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Mars reference ellipsoids.

The reference ellipsoid (``refellipsoid``) is set to model Mars,
folowing different models. The options are:

- "Sphere": A spherical planet. The radius is taken from a report of the IAU/IAG Working Group.
- "Ellipsoid": A reference ellipsoid with parameters taken from a report of the IAU/IAG Working Group.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_field"},

      .gin = {"model"},
      .gin_type = {"String"},
      .gin_value = {String("Sphere")},
      .gin_desc = {R"--(Model ellipsoid to use. Options listed above.)--"},

  };

  wsm_data["surface_fieldMoon"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Moon reference ellipsoids.

The reference ellipsoid (``refellipsoid``) is set to model Moon,
folowing different models. The options are:

- "Sphere":
    A spherical planet. The radius is taken from a
    report of the IAU/IAG Working Group.

- "Ellipsoid":
    A reference ellipsoid with parameters taken from
    Wikepedia (see code for details). The IAU/IAG working group
    defines the Moon ellipsoid to be a sphere.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_field"},

      .gin = {"model"},
      .gin_type = {"String"},
      .gin_value = {String("Sphere")},
      .gin_desc = {R"--(Model ellipsoid to use. Options listed above.)--"},

  };

  wsm_data["surface_fieldSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Make the surface field hold value at the key.
)--",
      .author = {"Richard Larsson"},
      .out = {"surface_field"},

      .in = {"surface_field"},
      .gin = {"value", "key"},
      .gin_type = {"Numeric, GriddedField2", "String"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Value to set)--", R"--(Key to set value at)--"},

  };

  wsm_data["surface_fieldSetProp"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Make the surface field hold value at the key.
)--",
      .author = {"Richard Larsson"},
      .out = {"surface_field"},

      .in = {"surface_field"},
      .gin = {"value", "key"},
      .gin_type = {"Numeric, GriddedField2", "String"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Value to set)--", R"--(Key to set value at)--"},

  };

  wsm_data["surface_fieldSetType"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Make the surface field hold value at the key.
)--",
      .author = {"Richard Larsson"},
      .out = {"surface_field"},

      .in = {"surface_field"},
      .gin = {"value", "key"},
      .gin_type = {"Numeric, GriddedField2", "String"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Value to set)--", R"--(Key to set value at)--"},

  };

  wsm_data["surface_fieldVenus"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Venus reference ellipsoids.

The reference ellipsoid (``refellipsoid``) is set to model Venus,
folowing different models. The options are:

- "Sphere":
      A spherical planet. The radius is taken from a
      report of the IAU/IAG Working Group.

According to the report used above, the Venus ellipsoid lacks
eccentricity and no further models should be required.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_field"},

      .gin = {"model"},
      .gin_type = {"String"},
      .gin_value = {String("Sphere")},
      .gin_desc = {R"--(Model ellipsoid to use. Options listed above.)--"},

  };

  wsm_data["surface_normalCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates the surface's local normal.

The default is to consider surface topography when calculating the
normal direction. That is, the variation of ``surface_elevation``
is allowed to affect the angles of *surface_normal*. This impact can
be deactivated by setting ``ignore_topography`` to 1. In this case,
the zenith angle of *surface_normal* becomes 0.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_normal"},

      .in = {"surface_field", "rtp_pos"},
      .gin = {"ignore_topography"},
      .gin_type = {"Index"},
      .gin_value = {Index{0}},
      .gin_desc =
          {R"--(Flag to control if surface slope is considered or not.)--"},

  };

  wsm_data["surface_rtpropFromTypesManual"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Extracts surface RT properties by manual selection of surface type.

The surface type to apply is selected by the GIN argument.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_type_mix",
              "surface_skin_t",
              "surface_los",
              "surface_rmatrix",
              "surface_emission"},

      .in = {"f_grid", "rtp_pos", "rtp_los", "surface_rtprop_agenda_array"},
      .gin = {"surface_type"},
      .gin_type = {"Index"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Selected surface type)--"},
      .pass_workspace = true,

  };

  wsm_data["surface_rtpropInterpFreq"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Interpolates surface RT properties in frequency.

The WSVs *surface_rmatrix* and *surface_emission* are inter-
polated linearly in frequency. The original frequency is given
by *f_grid*, and there is an interpolation to new frequency grid.
The function resets *f_grid* to the new grid.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"f_grid", "surface_rmatrix", "surface_emission"},

      .in = {"f_grid", "surface_rmatrix", "surface_emission"},
      .gin = {"f_new"},
      .gin_type = {"Vector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(New frequency grid)--"},

  };

  wsm_data["surface_rtprop_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *surface_rtprop_agenda* to a default value

Options are:

- ``"Blackbody_SurfTFromt_surface"``:

    1. Uses *InterpSurfaceFieldToPosition* using out=*surface_skin_t*, and field=``t_surface``
    2. Uses *surfaceBlackbody* to set *surface_los*, *surface_rmatrix*, and *surface_emission*, and also to modify *dsurface_rmatrix_dx*, and *dsurface_emission_dx*

- ``"Blackbody_SurfTFromt_field"``:

    1. Uses *InterpAtmFieldToPosition* using out=*surface_skin_t*, and field=``t_field``
    2. Uses *surfaceBlackbody* to set *surface_los*, *surface_rmatrix*, and *surface_emission*, and also to modify *dsurface_rmatrix_dx*, and *dsurface_emission_dx*

- ``"Specular_NoPol_ReflFix_SurfTFromt_surface"``:

    1. Uses *specular_losCalc* to set *specular_los*, and *surface_normal*
    2. Uses *InterpSurfaceFieldToPosition* using out=*surface_skin_t*, and field=``t_surface``
    3. Uses *surfaceFlatScalarReflectivity* to set *surface_los*, *surface_rmatrix*, and *surface_emission*, and also to modify *dsurface_rmatrix_dx*, and *dsurface_emission_dx*

- ``"Specular_NoPol_ReflFix_SurfTFromt_field"``:

    1. Uses *specular_losCalc* to set *specular_los*, and *surface_normal*
    2. Uses *InterpAtmFieldToPosition* using out=*surface_skin_t*, and field=``t_field``
    3. Uses *surfaceFlatScalarReflectivity* to set *surface_los*, *surface_rmatrix*, and *surface_emission*, and also to modify *dsurface_rmatrix_dx*, and *dsurface_emission_dx*

- ``"Specular_WithPol_ReflFix_SurfTFromt_surface"``:

    1. Uses *specular_losCalc* to set *specular_los*, and *surface_normal*
    2. Uses *InterpSurfaceFieldToPosition* using out=*surface_skin_t*, and field=``t_surface``
    3. Uses *surfaceFlatReflectivity* to set *surface_los*, *surface_rmatrix*, and *surface_emission*

- ``"lambertian_ReflFix_SurfTFromt_surface"``:

    1. Uses *specular_losCalc* to set *specular_los*, and *surface_normal*
    2. Uses *InterpSurfaceFieldToPosition* using out=*surface_skin_t*, and field=``t_surface``
    3. Uses *surfaceLambertianSimple* to set *surface_los*, *surface_rmatrix*, and *surface_emission*

- ``"lambertian_ReflFix_SurfTFromt_field"``:

    1. Uses *specular_losCalc* to set *specular_los*, and *surface_normal*
    2. Uses *InterpAtmFieldToPosition* using out=*surface_skin_t*, and field=``t_field``
    3. Uses *surfaceLambertianSimple* to set *surface_los*, *surface_rmatrix*, and *surface_emission*
)--",
      .author = {"Richard Larsson"},
      .out = {"surface_rtprop_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["surface_scalar_reflectivityFromSurface_rmatrix"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Sets *surface_scalar_reflectivity* based on *surface_rmatrix*.

For each frequency f, *surface_scalar_reflectivity* is set to
the sum of surface_rmatrix(joker,f,0,0).
)--",
          .author = {"Patrick Eriksson"},
          .out = {"surface_scalar_reflectivity"},

          .in = {"surface_rmatrix"},

      };

  wsm_data["telsemAtlasLookup"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Lookup SSMI emissivities from Telsem Atlas.

This returns the emissivities (indices [0,..,6])
for the SSMI channels that are contained in
the Telsem atlas.

If given latitude and longitude are not in the atlas an empty
vector is returned.
)--",
      .author = {"Simon Pfreundschuh"},

      .gout = {"emissivities"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(The SSMI emissivities from the atlas)--"},

      .gin = {"glat", "glon", "atlas"},
      .gin_type = {"Numeric", "Numeric", "TelsemAtlas"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt},
      .gin_desc = {R"--(The latitude for which to compute the emissivities.)--",
                   R"--(The latitude for which to compute the emissivities.)--",
                   R"--(The Telsem atlas to use.)--"},

  };

  wsm_data["telsem_atlasReadAscii"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads single TELSEM atlas.

'directory' needs to contain the original 12 Telsem atlas files
and the correlations file. This WSM reads the atlas for the specified
month and stores the result in the provided output atlas.
)--",
      .author = {"Simon Pfreundschuh"},

      .gout = {"atlas"},
      .gout_type = {"TelsemAtlas"},
      .gout_desc = {R"--(The atlas into which to store the loaded atlas.)--"},

      .gin = {"directory", "month", "filename_pattern"},
      .gin_type = {"String", "Index", "String"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    String("ssmi_mean_emis_climato_@MM@_cov_interpol_M2")},
      .gin_desc =
          {R"--(Directory with TELSEM 2 SSMI atlas files.)--",
           R"--(The month for which the atlas should be read.)--",
           R"--(Filename pattern (@MM@ gets replaced by month number))--"},

  };

  wsm_data["telsem_atlasesReadAscii"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads TELSEM atlas files.

'directory' needs to contain the original 12 Telsem atlas files
and the correlations file.
The whole data is combined into the WSV *telsem_atlases*
)--",
      .author = {"Oliver Lemke"},
      .out = {"telsem_atlases"},

      .gin = {"directory", "filename_pattern"},
      .gin_type = {"String", "String"},
      .gin_value = {std::nullopt,
                    String("ssmi_mean_emis_climato_@MM@_cov_interpol_M2")},
      .gin_desc =
          {R"--(Directory with TELSEM 2 SSMI atlas files.)--",
           R"--(Filename pattern (@MM@ gets replaced by month number))--"},

  };

  wsm_data["timeNow"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets time to system_clock::now().
)--",
      .author = {"Richard Larsson"},
      .out = {"time"},

  };

  wsm_data["timeOffset"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Offsets time for some seconds
)--",
      .author = {"Richard Larsson"},
      .out = {"time"},

      .in = {"time"},
      .gin = {"offset"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Time in seconds)--"},

  };

  wsm_data["timeSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the time.

The time format is similar to C:s strftime format
    %Y-%m-%d %H:%M:%S

The one exception is that our format accepts, but does not need to use,
decimals of the second in addtion.  Note that the native time resolution
for the decimal is the most that can be kept.  For some systems, this is
nano-seconds, and for others that is micro-seconds.  Please see with your
vendor.

A default argument is a close approximation to the formal first commit to
the ARTS codebase.  It is there to give an example of how the format looks.
)--",
      .author = {"Richard Larsson"},
      .out = {"time"},

      .gin = {"time_str"},
      .gin_type = {"String"},
      .gin_value = {String("2000-03-11 14:39:37.0")},
      .gin_desc = {R"--(A time stamp string in the default format)--"},

  };

  wsm_data["timeSleep"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sleeps until time has been reached.
)--",
      .author = {"Richard Larsson"},

      .in = {"time"},

  };

  wsm_data["time_gridOffset"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Offsets a time grid by some seconds.
)--",
      .author = {"Richard Larsson"},
      .out = {"time_grid"},

      .in = {"time_grid"},
      .gin = {"dt"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Time in seconds to add)--"},

  };

  wsm_data["time_stampsSort"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sort ``input`` by *time_stamps* into ``output``.
)--",
      .author = {"Richard Larsson"},

      .gout = {"output"},
      .gout_type = {"ArrayOfTime, ArrayOfVector"},
      .gout_desc = {R"--(Array sorted by time)--"},
      .in = {"time_stamps"},
      .gin = {"input"},
      .gin_type = {"ArrayOfTime, ArrayOfVector"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Array to sort of same size as *time_stamps*)--"},

  };

  wsm_data["timerStart"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initializes the CPU timer.
Use *timerStop* to stop the timer.

Usage example:

 - timerStart
 - ReadXML(f_grid,"frequencies.xml")
 - timerStop
 - Print(timer)
)--",
      .author = {"Oliver Lemke"},
      .out = {"timer"},

  };

  wsm_data["timerStop"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Stops the CPU timer.
See *timerStart* for example usage.
)--",
      .author = {"Oliver Lemke"},
      .out = {"timer"},

      .in = {"timer"},

  };

  wsm_data["transmittanceFromIy_aux"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Creates a vector of transmittance values.

The transmittances are set based on optical depths in *iy_aux*. That is,
one of the quantities in *iy_aux* must be "Optical depth".

The created vector has a length matching *f_grid* and can e.g. be used
as input to some of the FASTEM methods.
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"transmittance"},
      .gout_type = {"Vector"},
      .gout_desc = {R"--(Created vector of transmittance values.)--"},
      .in = {"iy_aux_vars", "iy_aux"},

  };

  wsm_data["water_p_eq_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *water_p_eq_agenda* to a default value

Options are:

- ``"MK05"``:
    1. Uses *water_p_eq_fieldMK05* to set *water_p_eq_field*
)--",
      .author = {"Richard Larsson"},
      .out = {"water_p_eq_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {String("MK05")},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["water_p_eq_fieldMK05"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculates *water_p_eq_field* according to Murphy and Koop, 2005.

Default is setting the saturation pressure to the one with respect
to water at temperatures >= 0C, and to the one with respect to ice
for <0C. The GIN ``only_liquid`` allows you to apply the liquid value
at all temperatures.

The saturation pressure with respect to liquid and ice water is
calculated according to Eq. 10 and 7, respectively, of:
Murphy, D. M., & Koop, T. (2005). Review of the vapour pressures of
ice and supercooled water for atmospheric applications. Quarterly
Journal of the Royal Meteorological Society, 131(608), 1539-1565.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"water_p_eq_field"},

      .in = {"atm_field"},
      .gin = {"only_liquid"},
      .gin_type = {"Index"},
      .gin_value = {Index{0}},
      .gin_desc =
          {R"--(Set to 1 to use liquid saturation pressure at all temperatures.)--"},

  };

  wsm_data["wind_u_fieldIncludePlanetRotation"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Maps the planet's rotation to an imaginary wind.

This method is of relevance if the observation platform is not
following the planet's rotation, and Doppler effects must be
considered. Examples include full disk observations from another
planet or a satellite not in orbit of the observed planet.

The rotation of the planet is not causing any Doppler shift for
1D and 2D simulations, and the method can only be used for 3D.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"atm_field"},

      .in = {"atm_field", "surface_field", "planet_rotation_period"},

  };

  wsm_data["x2artsAtmAndSurf"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Maps *x* to atmospheric and surface variables.

Maps OEM's state vector, *x*, to the matching ARTS variables. This
method handles atmospheric and surface variables. If you retrieve
other variables, make sure that you also call *x2artsSensor* and/or
*x2artsSpectroscopy*.

The following retrieval quantities are handled by this method:

 - Temperature
 - Absorption species
 - Scattering species
 - Winds
 - Surface variables

Should only be used inside *inversion_iterate_agenda*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"atm_field", "surface_field"},

      .in = {"atm_field",
             "surface_field",
             "jacobian_quantities",
             "x",
             "atmfields_checked",
             "atmgeom_checked",
             "abs_species",
             "cloudbox_on",
             "cloudbox_checked",
             "particle_bulkprop_names",
             "surface_props_names",
             "water_p_eq_agenda"},

      .pass_workspace = true,

  };

  wsm_data["x2artsSensor"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Maps *x* to sensor variables.

Maps OEM's state vector, *x*, to the matching ARTS variables. This
method handles variables associated with the sensor. If you retrieve
other variables, make sure that you also call *x2artsAtmAndSurf*
and/or *x2artsSpectroscopy*.

The following retrieval quantities are handled by this method:
 - Pointing
 - Frequency shift and stretch
 - Baseline fits

Should only be used inside *inversion_iterate_agenda*.

Elements in *x* representing pointing corrections are mapped to
*sensor_los*. Elements representing frequency corrections are mapped
to *f_backend*. Baseline variables are mapped to *y_baseline*.

The sensor response is recalculated if there is any non-zero frequency
correction.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"sensor_los",
              "f_backend",
              "y_baseline",
              "sensor_response",
              "sensor_response_f",
              "sensor_response_pol",
              "sensor_response_dlos",
              "sensor_response_f_grid",
              "sensor_response_pol_grid",
              "sensor_response_dlos_grid",
              "mblock_dlos"},

      .in = {"sensor_los",
             "f_backend",
             "sensor_response",
             "sensor_response_f",
             "sensor_response_pol",
             "sensor_response_dlos",
             "sensor_response_f_grid",
             "sensor_response_pol_grid",
             "sensor_response_dlos_grid",
             "mblock_dlos",
             "jacobian_quantities",
             "x",
             "sensor_response_agenda",
             "sensor_checked",
             "sensor_time"},

      .pass_workspace = true,

  };

  wsm_data["x2artsSpectroscopy"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Just defined to indicate a future extensiom.

Don't call the method, it will just generate an error.
)--",
      .author = {"Patrick Eriksson"},

  };

  wsm_data["xClip"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Clipping of the state vector.

The method allows you to apply hard limits the values of a
retrieval quantity. The retrieval quantity is specified by
``ijq``. All values of the quantity below ``limit_low``, are simply
set to ``limit_low``. And the same is performed with respect to
``limit_high``. That is, the data in x for the retrieval quantity
are forced to be inside the range [limit_low,limit_high].

Setting ijq=-1, is a shortcut for applying the limits on all
retrieval quantities.

Notice that limits must be specified in the unit used in *x*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"x"},

      .in = {"x", "jacobian_quantities"},
      .gin = {"ijq", "limit_low", "limit_high"},
      .gin_type = {"Index", "Numeric", "Numeric"},
      .gin_value = {std::nullopt, -std::numeric_limits<Numeric>::infinity(), std::numeric_limits<Numeric>::infinity()},
      .gin_desc = {R"--(Retrieval quantity index (zero-based))--",
                   R"--(Lower limit for clipping.)--",
                   R"--(Upper limit for clipping.)--"},

  };

  wsm_data["xaStandard"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Standard function for creating *xa*.

The method creates *xa* based on *jacobian_quantities* and the various
atmospheric fields. In the case of scattering species, the data are
taken from ``particle_bulkprop_field``. The following retrieval quantities
are handled:

 - Temperature
 - Absorption species
 - Scattering species
 - Pointing
 - Polynomial baseline fit
 - Sinusoidal baseline fit
)--",
      .author = {"Patrick Eriksson"},
      .out = {"xa"},

      .in = {"jacobian_quantities",
             "atmfields_checked",
             "atmgeom_checked",
             "atm_field",
             "abs_species",
             "cloudbox_on",
             "cloudbox_checked",
             "particle_bulkprop_names",
             "surface_field",
             "surface_props_names",
             "water_p_eq_agenda"},

      .pass_workspace = true,

  };

  wsm_data["yApplySensorPol"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Extraction of arbitrary linear polarisation.

This method shall be called after *yCalc* and then applies *sensor_pol*
on the output of *yCalc*. See *sensor_pol* for definition of the
polarisation responses. The *sensor_response* given to *yCalc* can not
contain any polarisation response, it must maintain original Stokes
elements. The value of ``stokes_dim`` must be >= 3.

The values in *sensor_pol* are applied on *y*, and *jacobian* if relevant.
*y_pol* is set following the values in *sensor_pol* but is rounded to
an integer value. Remaining data associated with *y* (e.g. y_pos) are
set to the value matching the first Stokes element.
)--",
      .author = {"Patrick Eriksson"},
      .out =
          {"y", "y_f", "y_pol", "y_pos", "y_los", "y_aux", "y_geo", "jacobian"},

      .in = {"y",
             "y_f",
             "y_pol",
             "y_pos",
             "y_los",
             "y_aux",
             "y_geo",
             "jacobian",
             "jacobian_do",
             "sensor_pos",
             "sensor_pol"},

  };

  wsm_data["yApplyUnit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Conversion of *y* to other spectral units.

Any conversion to brightness temperature is normally made inside
*yCalc*. This method makes it possible to also make this conversion
after *yCalc*, but with restrictions for *jacobian* and with.
respect to the n2-law of radiance.

The conversion made inside ``iyEmissionStandard`` is mimiced
and see that method for constraints and selection of output units.
This with the restriction that the n2-law can be ignored. The later
is the case if the sensor is placed in space, or if the refractive
only devaites slightly from unity.

The method handles *y* and *jacobian* in parallel, where
the last variable is only considered if it is set. The
input data must be in original radiance units. A completely
stringent check of this can not be performed.

The method can not be used with jacobian quantities that are not
obtained through radiative transfer calculations. One example on
quantity that can not be handled is *jacobianAddPolyfit*. There
are no automatic checks warning for incorrect usage!

If you are using this method, *iy_unit* should be set to "1" when
calling *yCalc*, and be changed before calling this method.

Conversion of *y_aux* is not supported.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"y", "jacobian"},

      .in = {"y", "jacobian", "y_f", "y_pol", "iy_unit"},

  };

  wsm_data["yCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculation of complete measurement vectors (y).

The method performs radiative transfer calculations from a sensor
perspective. Radiative transfer calculations are performed for
monochromatic pencil beams, following *iy_main_agenda* and
associated agendas. Obtained radiances are weighted together by
*sensor_response*, to include the characteristics of the sensor.
The measurement vector obtained can contain anything from a single
frequency value to a series of measurement scans (each consisting
of a series of spectra), all depending on the settings. Spectra
and jacobians are calculated in parallel.

The frequency, polarisation etc. for each measurement value is
given by *y_f*, *y_pol*, *y_pos* and *y_los*.

The content of *y_aux* follows *iy_aux_vars*. See the method selected
for *iy_main_agenda* for allowed choices.

The geo-positions (*y_geo*) are set based on *sensor_response*. When
an antenna pattern is considered, there are several pencil beams,
and thus also several goe-positions, associated with each value of *y*.
The geo-position assigned to a value in *y* is the *geo_pos* of the pencil
beam related to the highest value in *sensor_response*. This means that
*mblock_dlos* must contain the bore-sight direction (0,0), if you
want *y_geo* to exactly match the bore-sight direction.

The Jacobian provided (*jacobian*) is adopted to selected retrieval
units, but no transformations are applied. Transformations are
included by calling *jacobianAdjustAndTransform*.
)--",
      .author = {"Patrick Eriksson"},
      .out =
          {"y", "y_f", "y_pol", "y_pos", "y_los", "y_aux", "y_geo", "jacobian"},

      .in = {"atmgeom_checked",
             "atmfields_checked",
             "atm_field",
             "cloudbox_on",
             "cloudbox_checked",
             "scat_data_checked",
             "sensor_checked",
             "f_grid",
             "sensor_pos",
             "sensor_los",
             "transmitter_pos",
             "mblock_dlos",
             "sensor_response",
             "sensor_response_f",
             "sensor_response_pol",
             "sensor_response_dlos",
             "iy_unit",
             "iy_main_agenda",
             "jacobian_agenda",
             "jacobian_do",
             "jacobian_quantities",
             "iy_aux_vars"},

      .pass_workspace = true,

  };

  wsm_data["yCalcAppend"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Replaces *yCalc* if a measurement shall be appended to an
existing one.

The method works basically as *yCalc* but appends the results to
existing data, instead of creating completely new *y* and its
associated variables. This method is required if your measurement
consists of data from two instruments using different observation
techniques (corresponding to different iyCalc-methods). One such
example is if emission and transmittance data are combined into a
joint retrieval. The method can also be used to get around the
constrain that *sensor_response* is required to be the same for
all data.

The new measurement is simply appended to the input *y*, and the
other output variables are treated correspondingly. Data are
appended "blindly" in *y_aux*. That is, data of different type
are appended if *iy_aux_vars* differs between the two measurements,
the data are appended strictly following the order. First variable
of second measurement is appended to first variable of first
measurement, and so on. The number of auxiliary variables can differ
between the measurements. Missing data are set to zero.

The set of retrieval quantities can differ between the two
calculations. If an atmospheric quantity is part of both Jacobians,
the same retrieval grids must be used in both cases.
The treatment of instrument related Jacobians (baseline fits,
pointing ...) follows the ``append_instrument_wfs`` argument.

A difference to *yCalc* is that *jacobian_quantities* is both in-
and output variable. The input version shall match the measurement
to be calculated, while the output version matches the output *y*,
the combined, measurements. A copies of *jacobian_quantities* of the
first measurement must be made and shall be provided to the method
as ``jacobian_quantities_copy``.

As for *yCalc* Jacobian transformations are not handled, and the
the input Jacobian shall not contain transformations. That is
*jacobianAdjustAndTransform* shall be called after this method,
when the complete Jacobian is at hand.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"y",
              "y_f",
              "y_pol",
              "y_pos",
              "y_los",
              "y_aux",
              "y_geo",
              "jacobian",
              "jacobian_quantities"},

      .in = {"y",
             "y_f",
             "y_pol",
             "y_pos",
             "y_los",
             "y_aux",
             "y_geo",
             "jacobian",
             "atmgeom_checked",
             "atmfields_checked",
             "atm_field",
             "cloudbox_on",
             "cloudbox_checked",
             "scat_data_checked",
             "sensor_checked",
             "f_grid",
             "sensor_pos",
             "sensor_los",
             "transmitter_pos",
             "mblock_dlos",
             "sensor_response",
             "sensor_response_f",
             "sensor_response_pol",
             "sensor_response_dlos",
             "iy_unit",
             "iy_main_agenda",
             "jacobian_agenda",
             "jacobian_do",
             "jacobian_quantities",
             "iy_aux_vars"},
      .gin = {"jacobian_quantities_copy", "append_instrument_wfs"},
      .gin_type = {"ArrayOfRetrievalQuantity", "Index"},
      .gin_value = {std::nullopt, Index{0}},
      .gin_desc =
          {R"--(Copy of *jacobian_quantities* of first measurement.)--",
           R"--(Flag controlling if instrumental weighting functions are appended or treated as different retrieval quantities.)--"},
      .pass_workspace = true,

  };

  wsm_data["yColdAtmHot"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Computes *y* from input using standard calibration scheme of
cold-atm-hot observations

If calib evaluates as true:
    y = cold_temp + (hot_temp - cold_temp) * (atm - cold) / (hot - cold)

If calib evaluates as false:
    y = (hot_temp * cold - cold_temp * hot) / (hot - cold)
)--",
      .author = {"Richard Larsson"},
      .out = {"y"},

      .gin = {"cold", "atm", "hot", "cold_temp", "hot_temp", "calib"},
      .gin_type = {"Vector", "Vector", "Vector", "Numeric", "Numeric", "Index"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    Index{1}},
      .gin_desc =
          {R"--(N-elem Vector of cold load linear power)--",
           R"--(N-elem Vector of atmosphere linear power)--",
           R"--(N-elem Vector of hot load linear power)--",
           R"--(Cold load temperature)--",
           R"--(Hot load temperature)--",
           R"--(Flag for calibration scheme, false means system temperature is computed)--"},

  };

  wsm_data["yDoublingMeanFocus"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Focus in on *y* around some *f_grid*, then sets *f_grid* to
same focus

The algorithm will double the frequency spacing DF
steps away from F0, doubling it again DF steps away
from F0+DF, and so on until the end of the range.
The same happens but reversely towards F0-DF, ...

Inside these ranges, the values will be averaged
so that there's one value per original input
between F0-DF and F0+DF, 1 value per 2 original values
between F0+DF and F0+2*DF, and so on ever doubling the
number of original values per output value

F0 and DF are set inside the function to either the values
given by the user, or if they are non-positive as mean(f_grid)
and 10 * (f_grid[1] - f_grid[0]), respectively

Ignores NaNs and infinities in averaging calculations.
)--",
      .author = {"Richard Larsson"},
      .out = {"f_grid", "y"},

      .in = {"f_grid", "y"},
      .gin = {"f0", "df"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {Numeric{-1}, Numeric{-1}},
      .gin_desc = {R"--(User input for F0 [see description for default])--",
                   R"--(User input for DF [see description for default])--"},

  };

  wsm_data["yMaskOutsideMedianRange"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Masks values not within the range as NaN::

  [median(y) - dx, median(y) + dx]

Ignores NaNs in median calculations.
)--",
      .author = {"Richard Larsson"},
      .out = {"y"},

      .in = {"y"},
      .gin = {"dx"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Range plus-minus the median of unmasked values)--"},

  };

  wsm_data["yRadar"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Replaces *yCalc* for radar/lidar calculations.

The output format for *iy* when simulating radars and lidars differs
from the standard one, and *yCalc* can not be used for such simulations.
This method works largely as *yCalc*, but is tailored to handle the
output from *iyRadarSingleScat*. Note that *iy_radar_agenda* replaces
*iy_main_agenda*.

The method requires additional information about the sensor,
regarding its recieving properties. First of all, recieved
polarisation states are taken from *instrument_pol_array*. Note
that this WSV allows to define several measured polarisations
for each transmitted signal. For example, it is possible to
simulate transmittance of V and measuring backsacttered V and H.

Secondly, the range averaging is described by *range_bins*. These
bins can either be specified in altitude or two-way travel time.
In both case, the edges of the range bins shall be specified.
All data (including auxiliary variables) are returned as the
average inside the bins. If a bin is totally outside the model
atmosphere, NaN is returned.

The options for *iy_unit_radar* are:

- ``"1"``:
    Backscatter coefficient. Unit is 1/(m*sr). At zero
    attenuation, this equals the scattering matrix value for
    the backward direction. See further AUG.
- ``"Ze"``: Equivalent reflectivity. Unit is mm^6/m^3. Conversion formula is given below.
- ``"dBZe"``: 10*log10(Ze/Z0), where Z0 is 1 mm^6/m^3.

The conversion from backscatter coefficient to Ze is::

  Ze = 1e18 * lambda^4 / (k2 * pi^5) * sum(sigma)

where sum(sigma) = 4 * pi * b, and b is the backscatter coefficient.

The reference dielectric factor can either specified directly by
the argument ``k2``. For example, to mimic the CloudSat data, ``k2``
shall be set to 0.75 (citaion needed). If ``k2`` is set to be 
negative (which is defualt), k2 is calculated as::

  k2 = abs( (n^2-1)/(n^2+2) )^2

where n is the refractive index of liquid water at temperature
``ze_tref`` and the frequency of the radar, calculated by the MPM93
parameterization.

A lower limit for dBZe is applied (``dbze_min``). The main reason is to
handle the fact that dBZe is not defined for Ze=0, and dBZe is set to
the clip value when Ze < 10^(dbze_min/10).
)--",
      .author = {"Patrick Eriksson"},
      .out =
          {"y", "y_f", "y_pol", "y_pos", "y_los", "y_aux", "y_geo", "jacobian"},

      .in = {"atmgeom_checked",
             "atmfields_checked",
             "iy_unit_radar",
             "iy_aux_vars",
             "f_grid",
             "cloudbox_on",
             "cloudbox_checked",
             "sensor_pos",
             "sensor_los",
             "sensor_checked",
             "jacobian_do",
             "jacobian_quantities",
             "iy_radar_agenda",
             "instrument_pol_array",
             "range_bins"},
      .gin = {"ze_tref", "k2", "dbze_min"},
      .gin_type = {"Numeric", "Numeric", "Numeric"},
      .gin_value = {Numeric{273.15}, Numeric{-1}, Numeric{-99}},
      .gin_desc = {R"--(Reference temperature for conversion to Ze.)--",
                   R"--(Reference dielectric factor.)--",
                   R"--(Clip value for dBZe.)--"},
      .pass_workspace = true,

  };

  wsm_data["ySimpleSpectrometer"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Converts *iy* to *y* assuming a fixed frequency resolution.

This is a short-cut, avoiding *yCalc*, that can be used to convert
monochromatic pencil beam data to spectra with a fixed resolution.

The method mimics a spectrometer with rectangular response
functions, all having the same width (``df``). The position of
the first spectrometer channel is set to f_grid[0] + df / 2.
The centre frequency of channels are returned as *y_f*.

Auxiliary variables and *jacobian* s are not handled.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"y", "y_f"},

      .in = {"iy", "f_grid"},
      .gin = {"df"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Selected frequency resolution.)--"},

  };

  wsm_data["y_geo_seriesFromY_geo"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Fills *y_geo_series* with data from *y_geo*.

The geo-position is taken from the first channel. There is no check
that the other channels have identical data in *y_geo*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"y_geo_series"},

      .in = {"y_geo", "sensor_response_f_grid"},

  };

  wsm_data["y_geo_swathFromY_geo"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Fills *y_geo_series* with data from *y_geo*.

The geo-position is taken from the first channel. There is no check
that the other channels have identical data in *y_geo*.

The method assumes the same order in *y* as *y_swathFromY*.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"y_geo_swath"},

      .in = {"y_geo", "sensor_response_f_grid"},
      .gin = {"npixel"},
      .gin_type = {"Index"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Number of pixels per swath.)--"},

  };

  wsm_data["y_seriesFromY"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Fills *y_series* with data from *y*.

The method basically reshapes *y* to fit *y_series*.

Default is to check that *y_f* does not change between posistions,
i.e. that the channel frequencies do not vary.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"y_series"},

      .in = {"y", "y_f", "sensor_response_f_grid"},
      .gin = {"safe"},
      .gin_type = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc =
          {R"--(Flag for checking that channels do not vary in frequency.)--"},

  };

  wsm_data["y_swathFromY"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Fills *y_swath* with data from *y*.

The method basically reshapes *y* to fit *y_swath*. It is assumed
that swath forms the outermost loop in *y*. That is, first in *y*
are the data for the first swath etc. The number of pixels per swath
must be specified manually by a GIN parameter.

To set *sensor_pos* and *sensor_los* having data organised in swath
format, use *MatrixReshapeTensor3*.

Default is to check that *y_f* does not change between posistions,
i.e. that the channel frequencies do not vary.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"y_swath"},

      .in = {"y", "y_f", "sensor_response_f_grid"},
      .gin = {"npixel", "safe"},
      .gin_type = {"Index", "Index"},
      .gin_value = {std::nullopt, Index{1}},
      .gin_desc =
          {R"--(Number of pixels per swath.)--",
           R"--(Flag for checking that channels do not vary in frequency.)--"},

  };

  wsm_data["ybatchCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Performs batch calculations for the measurement vector y.

We perform *ybatch_n* jobs, starting at index *ybatch_start*. (Zero
based indexing, as usual.) The output array *ybatch* will have
ybatch_n elements. Indices in the output array start
with zero, independent of *ybatch_start*.

The method performs the following:

   (1) Sets *ybatch_index* = *ybatch_start*.
   (2) Performs a-d until *ybatch_index* = *ybatch_start* + *ybatch_n*.
        a. Executes *ybatch_calc_agenda*.
        b. If *ybatch_index* = *ybatch_start*, resizes *ybatch*
           based on *ybatch_n* and length of *y*.
        c. Copies *y* to *ybatch_index* - *ybatch_start*
           of *ybatch*.
        d. Adds 1 to *ybatch_index*.

Beside the *ybatch_calc_agenda*, the WSVs *ybatch_start*
and *ybatch_n* must be set before calling this method.
Further, *ybatch_calc_agenda* is expected to produce a
spectrum and should accordingly include a call of *yCalc*
(or asimilar method).

The input variable *ybatch_start* is set to a default of zero.

Jacobians are also collected, and stored in output variable *ybatch_jacobians*. 
(This will be empty if *yCalc* produces empty Jacobians.)

See the user guide for further practical examples.
)--",
      .author = {"Stefan Buehler"},
      .out = {"ybatch", "ybatch_aux", "ybatch_jacobians"},

      .in = {"ybatch_start", "ybatch_n", "ybatch_calc_agenda"},
      .gin = {"robust"},
      .gin_type = {"Index"},
      .gin_value = {Index{0}},
      .gin_desc =
          {R"--(A flag with value 1 or 0. If set to one, the batch calculation will continue, even if individual jobs fail. In that case, a warning message is written to screen and file (out1 output stream), and the *y* Vector entry for the failed job in *ybatch* is left empty.)--"},
      .pass_workspace = true,

  };

  wsm_data["ybatchColdAtmHotAtmCycle"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Computes *ybatch* from input using standard calibration scheme of
a cycle through cold-atm-hot-atm-cold-... observations

Computes for every full cycle reaching a new hot or cold:
    y = cold_temp + (hot_temp - cold_temp) * (atm - cold) / (hot - cold)

Assumes data is ordered as Cold-Atm-Hot-Atm-Cold-Atm-Hot-Atm-...,
but Cold does not have to be at data[0], instead the first cold
position is set by ``first_c_index``, which defaults to 0 but can be any positive
index so that *level0_data*[``first_c_index``] is a cold-measurements.  Note that if
``first_c_index`` is larger than 1, then the first output data will be around the
observation cycle -HAC-, where H is at ``first_c_index``-2

Also returns the times of the Atm measurements in *sensor_time*
if the measurement's time data is provided
)--",
      .author = {"Richard Larsson"},
      .out = {"ybatch", "sensor_time"},

      .in = {"level0_data", "level0_time"},
      .gin = {"cold_temp", "hot_temp", "first_c_index"},
      .gin_type = {"Vector", "Vector", "Index"},
      .gin_value = {std::nullopt, std::nullopt, Index{0}},
      .gin_desc =
          {R"--(Cold load calibration temperature (must match level0_data length))--",
           R"--(Hot load calibration temperature (must match level0_data length))--",
           R"--(Index offset of the first cold position)--"},

  };

  wsm_data["ybatchDoublingMeanFocus"] = WorkspaceMethodInternalRecord{
      .desc = R"--(See *yDoublingMeanFocus*
)--",
      .author = {"Richard Larsson"},
      .out = {"f_grid", "ybatch"},

      .in = {"f_grid", "ybatch"},
      .gin = {"f0", "df"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {Numeric{-1}, Numeric{-1}},
      .gin_desc = {R"--(User input for F0 [see description for default])--",
                   R"--(User input for DF [see description for default])--"},

  };

  wsm_data["ybatchMaskOutsideMedianRange"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Apply *yMaskOutsideMedianRange* for each *y* in *ybatch*
)--",
      .author = {"Richard Larsson"},
      .out = {"ybatch"},

      .in = {"ybatch"},
      .gin = {"dx"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Range plus-minus the median of unmasked values)--"},

  };

  wsm_data["ybatchMetProfilesClear"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(This method is used for simulating ARTS for metoffice model fields
for clear sky conditions.

This method reads in *met_amsu_data* which contains the
lat-lon of the metoffice profile files as a Matrix. It then
loops over the number of profiles and corresponding to each
longitude create the appropriate profile basename. Then,
Corresponding to each basename we have temperature field, altitude
field, humidity field, and particle number density field. The
temperature field and altitude field are stored in the same dimensions
as ``t_field_raw`` and ``z_field_raw``. The oxygen and nitrogen VMRs are
set to constant values of 0.209 and 0.782, respectively and are used
along with humidity field to generate ``vmr_field_raw``. 

The three fields ``t_field_raw``, ``z_field_raw``, and ``vmr_field_raw`` are
given as input to *met_profile_calc_agenda* which is called in this
method. See documentation of WSM *met_profile_calc_agenda* for more
information on this agenda. 

The method also converts satellite zenith angle to appropriate
*sensor_los*. It also sets the ``p_grid`` and *cloudbox_limits*
from the profiles inside the function
)--",
      .author = {"Seerekha T.R."},
      .out = {"ybatch"},

      .in = {"abs_species",
             "met_profile_calc_agenda",
             "f_grid",
             "met_amsu_data",
             "sensor_pos",
             "surface_field"},
      .gin = {"nelem_p_grid", "met_profile_path"},
      .gin_type = {"Index", "String"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(FIXME DOC)--", R"--(FIXME DOC)--"},
      .pass_workspace = true,

  };

  wsm_data["ybatchTimeAveraging"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Time average of *ybatch* and *sensor_time*
)--",
      .author = {"Richard Larsson"},
      .out = {"ybatch", "sensor_time"},

      .in = {"ybatch", "sensor_time"},
      .gin = {"time_step", "disregard_first", "disregard_last"},
      .gin_type = {"String", "Index", "Index"},
      .gin_value = {std::nullopt, Index{0}, Index{0}},
      .gin_desc =
          {R"--(Time step in the form "INDEX SCALE", where SCALE is "h", "min", or "s" for hours, minutes or seconds)--",
           R"--(Flag to remove first time step (e.g., if it is an incomplete step))--",
           R"--(Flag to remove last time step (e.g., if it is an incomplete step))--"},

  };

  wsm_data["ybatchTroposphericCorrectionNaiveMedianForward"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Performs naive tropospheric corrections on *ybatch*

Sets *ybatch_corr* to be able to perform the inverse of the corrections,
each array-element with 3 entries as [median, part_trans, trop_temp]

Uses the same tropospheric temperature for all values if trop_temp.nelem()==1
)--",
          .author = {"Richard Larsson"},
          .out = {"ybatch_corr", "ybatch"},

          .in = {"ybatch"},
          .gin = {"range", "trop_temp", "targ_temp"},
          .gin_type = {"ArrayOfIndex", "Vector", "Numeric"},
          .gin_value = {ArrayOfIndex{ArrayOfIndex(0)},
                        std::nullopt,
                        Numeric{2.73}},
          .gin_desc =
              {R"--(Positions where the median of the baseline is computed, if empty all is used)--",
               R"--(Radiative temperature of the troposphere [dim: 1 or ybatch.nelem()])--",
               R"--(Temperature target of the baseline)--"},

      };

  wsm_data["ybatchTroposphericCorrectionNaiveMedianInverse"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Performs inverse of naive tropospheric corrections on *ybatch*
)--",
          .author = {"Richard Larsson"},
          .out = {"ybatch"},

          .in = {"ybatch", "ybatch_corr"},

      };

  wsm_data["ybatch_calc_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *ybatch_calc_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--",
      .author = {"Richard Larsson"},
      .out = {"ybatch_calc_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  return wsm_data;
}