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
                   R"--(Additional error message.)--"}};

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
                   R"--(Additional error message.)--"}};

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
      .gin_desc = {R"--(Source variable.)--"}};

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

  wsm_data["Delete"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Deletes a workspace variable.

The variable is not deleted from the workspace, but it is
reset to its default value. This is useful if you want to
free memory for heavy variables
)--",
      .author = {"Oliver Lemke"},

      .gout = {"v"},
      .gout_type = {"Any"},
      .gout_desc = {R"--(Variable to be deleted.)--"}};

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
      .gout_type = {"NamedGriddedField3"},
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

  wsm_data["surface_pointFromAtm"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Extracts a surface point from an atmospheric field

The elevation is from *path_point*, the normal points straight up, and
the atmosphere must contain the temperature.  If the wind is present,
it is also extracted.
)--",
      .author = {"Richard Larsson"},
      .out = {"surface_point"},
      .in = {"atm_field", "path_point"},
  };

  wsm_data["InterpSurfaceFieldToPosition"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Point interpolation of surface fields.

The default way to specify the position is by *path_point*.

Linear interpolation is applied.

The interpolation is done for the latitude and longitude in
*path_point*, while the altitude in *path_point* is not part of the
calculations. However, it is checked that the altitude of *path_point*
is inside the range covered by ``z_surface`` with a 1 m margin, to
give a warning when the specified position is not consistent with
the surface altitudes.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"surface_point"},

      .in = {"path_point", "surface_field", "surface_search_accuracy"},

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
      .gin_value = {std::nullopt,
                    -std::numeric_limits<Numeric>::infinity(),
                    std::numeric_limits<Numeric>::infinity()},
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
      .gin_desc = {R"--(Positions (grid) where ``y`` given.)--",
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

  wsm_data["PlanetSet"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets *g0_agenda*, ``refellipsoid``, *molarmass_dry_air*, and *planet_rotation_period* to default values

*g0_agenda* is set using *g0_agendaSet* with the same option
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
      .desc = R"--(Prints a list of the workspace variables to standard out.

Each will have a short excerpt of the variable's content
)--",
      .author = {"Oliver Lemke"},
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
``particle_bulkpropRadarOnionPeeling``. See that method for
format of the table.

The method needs to be called twice to form a complete table,
once for liquid and ice hydrometeors. The table can be empty at
the first call.

The input data (*scat_data* etc.) must match two scattering
species and a single frequency (the one of the radar).
)--",
      .author = {"Patrick Eriksson"},

      .gout = {"invtable"},
      .gout_type = {"ArrayOfNamedGriddedField2"},
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
           R"--(Reference temperature for conversion to Ze. See further ``yRadar``.)--",
           R"--(Reference dielectric factor. See further ``yRadar``.)--"},
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
      .gin_desc = {R"--(Name of the NetCDF file.)--"}};

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
      .gin_desc = {R"--(Name of the XML file.)--"}};

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
      .gin_desc = {
          R"--(File name. See above.)--",
          R"--(Equalize the widths of all numbers by padding with zeros as necessary. 0 means no padding (default).)--"}};

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
          {"ArrayOfAbsorptionLines, ArrayOfAgenda, ArrayOfArrayOfAbsorptionLines, ArrayOfArrayOfGriddedField1, ArrayOfArrayOfGriddedField2, ArrayOfArrayOfGriddedField3, ArrayOfArrayOfIndex, ArrayOfArrayOfMatrix, ArrayOfArrayOfMuelmatMatrix, ArrayOfArrayOfMuelmatVector, ArrayOfArrayOfPropmatMatrix, ArrayOfArrayOfPropmatVector, ArrayOfArrayOfScatteringMetaData, ArrayOfArrayOfSingleScatteringData, ArrayOfArrayOfSpeciesTag, ArrayOfArrayOfStokvecMatrix, ArrayOfArrayOfStokvecVector, ArrayOfArrayOfString, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6, ArrayOfArrayOfTime, ArrayOfArrayOfVector, ArrayOfAtmPoint, ArrayOfCIARecord, ArrayOfGriddedField1, ArrayOfGriddedField2, ArrayOfGriddedField3, ArrayOfGriddedField4, ArrayOfIndex, ArrayOfMatrix, ArrayOfMuelmatMatrix, ArrayOfMuelmatVector, ArrayOfPropmatMatrix, ArrayOfPropmatVector, ArrayOfQuantumIdentifier, ArrayOfScatteringMetaData, ArrayOfSingleScatteringData, ArrayOfSparse, ArrayOfSpeciesTag, ArrayOfStokvecMatrix, ArrayOfStokvecVector, ArrayOfString, ArrayOfSun, ArrayOfTelsemAtlas, ArrayOfTensor3, ArrayOfTensor4, ArrayOfTensor5, ArrayOfTensor6, ArrayOfTensor7, ArrayOfTime, ArrayOfVector, ArrayOfXsecRecord, Vector, Matrix, Sparse"},
      .gout_desc =
          {R"--(Selected elements. Must have the same variable type as haystack.)--"},

      .gin = {"haystack", "needleindexes"},
      .gin_type = {"ArrayOfAbsorptionLines, ArrayOfAgenda, ArrayOfArrayOfAbsorptionLines, ArrayOfArrayOfGriddedField1, ArrayOfArrayOfGriddedField2, ArrayOfArrayOfGriddedField3, ArrayOfArrayOfIndex, ArrayOfArrayOfMatrix, ArrayOfArrayOfMuelmatMatrix, ArrayOfArrayOfMuelmatVector, ArrayOfArrayOfPropmatMatrix, ArrayOfArrayOfPropmatVector, ArrayOfArrayOfScatteringMetaData, ArrayOfArrayOfSingleScatteringData, ArrayOfArrayOfSpeciesTag, ArrayOfArrayOfStokvecMatrix, ArrayOfArrayOfStokvecVector, ArrayOfArrayOfString, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6, ArrayOfArrayOfTime, ArrayOfArrayOfVector, ArrayOfAtmPoint, ArrayOfCIARecord, ArrayOfGriddedField1, ArrayOfGriddedField2, ArrayOfGriddedField3, ArrayOfGriddedField4, ArrayOfIndex, ArrayOfMatrix, ArrayOfMuelmatMatrix, ArrayOfMuelmatVector, ArrayOfPropmatMatrix, ArrayOfPropmatVector, ArrayOfQuantumIdentifier, ArrayOfScatteringMetaData, ArrayOfSingleScatteringData, ArrayOfSparse, ArrayOfSpeciesTag, ArrayOfStokvecMatrix, ArrayOfStokvecVector, ArrayOfString, ArrayOfSun, ArrayOfTelsemAtlas, ArrayOfTensor3, ArrayOfTensor4, ArrayOfTensor5, ArrayOfTensor6, ArrayOfTensor7, ArrayOfTime, ArrayOfVector, ArrayOfXsecRecord, Vector, Matrix, Sparse",
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
      .gin_value = {std::nullopt,
                    -std::numeric_limits<Numeric>::infinity(),
                    std::numeric_limits<Numeric>::infinity()},
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

The vectors ``x`` and ``y`` can be the same variable.
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
                   R"--(Name of the NetCDF file.)--"}};

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
                   R"--(Name of the NetCDF file.)--"}};

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
      .gin_desc = {
          R"--(Variable to be saved.)--",
          R"--(Name of the XML file.)--",
          R"--(0: Overwrite existing files, 1: Use unique filenames)--"}};

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
      .gin_desc = {
          R"--(Workspace variable to be saved.)--",
          R"--(File name. See above.)--",
          R"--(Equalize the widths of all numbers by padding with zeros as necessary. 0 means no padding (default).)--"}};

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
      .desc =
          R"--(Initialize the atmospheric field with some altitude and isotopologue ratios
)--",
      .author = {"Richard Larsson"},
      .out = {"atm_field"},

      .gin = {"toa", "default_isotopologue"},
      .gin_type = {"Numeric", "String"},
      .gin_value = {std::nullopt, String{"Builtin"}},
      .gin_desc = {R"--(Top of atmosphere altitude [m].)--",
                   "Default option for the isotopologue ratios"},
  };

  wsm_data["atm_pointInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initialize an atmospheric point with some isotopologue ratios
)--",
      .author = {"Richard Larsson"},
      .out = {"atm_point"},

      .gin = {"default_isotopologue"},
      .gin_type = {"String"},
      .gin_value = {String{"Builtin"}},
      .gin_desc = {"Default option for the isotopologue ratios"},
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

  wsm_data["background_transmittanceFromPathPropagationBack"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Sets *background_transmittance* to back of *ppvar_cumtramat*
)--",
          .author = {"Richard Larsson"},
          .out = {"background_transmittance"},

          .in = {"ppvar_cumtramat"},

      };

  wsm_data["background_transmittanceFromPathPropagationFront"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Sets *background_transmittance* to front of *ppvar_cumtramat*
)--",
          .author = {"Richard Larsson"},
          .out = {"background_transmittance"},

          .in = {"ppvar_cumtramat"},

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

  wsm_data["spectral_radiance_background_space_agendaSet"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Sets *spectral_radiance_background_space_agenda*

Options are:

- ``"UniformCosmicBackground:``:

  1. Calls *spectral_radiance_backgroundUniformCosmicBackground*
  2. Calls *spectral_radiance_background_jacobianEmpty*
)--",
          .author = {"Richard Larsson"},
          .out = {"spectral_radiance_background_space_agenda"},

          .gin = {"option"},
          .gin_type = {"String"},
          .gin_value = {String{"UniformCosmicBackground"}},
          .gin_desc = {R"--(Default agenda option (see description))--"},
      };

  wsm_data["spectral_radiance_background_surface_agendaSet"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Sets *spectral_radiance_background_surface_agenda*

Options are:

- ``"Blackbody"``:

  1. Calls *spectral_radiance_backgroundSurfaceBlackbody*
)--",
          .author = {"Richard Larsson"},
          .out = {"spectral_radiance_background_surface_agenda"},

          .gin = {"option"},
          .gin_type = {"String"},
          .gin_value = {String{"Blackbody"}},
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

      .in = {"ecs_data", "atm_field"},

  };

  wsm_data["ecs_dataAddMakarov2020NEWNEW"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the O2-66 microwave band data for ECS.
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data2"},
      .in = {"ecs_data2"},
  };

  wsm_data["ecs_dataAddTran2011NEWNEW"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the CO2-626, CO2-628, and CO2-636 band data for ECS.

Sets CO2 species
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data2"},
      .in = {"ecs_data2"},
  };

  wsm_data["ecs_dataAddRodrigues1997NEWNEW"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the CO2-626, CO2-628, and CO2-636 band data for ECS.

Sets N2 and O2 speces
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data2"},
      .in = {"ecs_data2"},
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

      .in = {"ecs_data", "atm_field"},

  };

  wsm_data["ecs_dataAddSpeciesData"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets ECS data for one set of species and quantum identifiers.
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data"},

      .in = {"ecs_data", "atm_field"},
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

      .in = {"ecs_data", "atm_field"},

  };

  wsm_data["ecs_dataAddTran2011"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the CO2-626, CO2-636, and CO2-628 IR band data for ECS.
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data"},

      .in = {"ecs_data", "atm_field"},

  };

  wsm_data["ecs_dataInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Resets/initializes the ECS data.
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data"},

  };

  wsm_data["ecs_dataInitNEWNEW"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Resets/initializes the ECS data.
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data2"},

  };

  wsm_data["ecs_dataAddMeanAirNEWNEW"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets ECS data for air from other data if available.
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data2"},

      .in = {"ecs_data2"},
      .gin = {"vmrs", "species"},
      .gin_type = {"Vector", "ArrayOfSpecies"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(VMRs of air species)--", R"--(Air species)--"},

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

  wsm_data["lbl_checkedCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Checks that the line-by-line parameters are OK.

On failure, will throw.  On success, lbl_checked evals as true

Note that checks may become more stringent as ARTS evolves, especially for
"new" options.  This test might succeed in one version of ARTS but fail
in later versions
)--",
      .author = {"Richard Larsson"},
      .out = {"lbl_checked"},

      .in = {"abs_lines_per_species", "abs_species", "atm_field"},

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
``DoitScatteringDataPrepare``.

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

      .in = {"scat_data", "f_grid", "cloudbox_limits", "jacobian_targets"},

  };

  wsm_data["ppvar_atmFromPath"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gets the atmospheric points along the path.
)--",
      .author = {"Richard Larsson"},
      .out = {"ppvar_atm"},
      .in = {"propagation_path", "atm_field"},
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
      .in = {"f_grid", "propagation_path", "ppvar_atm", "rte_alonglos_v"},
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
             "jacobian_targets",
             "ppvar_f",
             "propagation_path",
             "ppvar_atm"},
      .pass_workspace = true,
  };

  wsm_data["spectral_radiance_pathCalcEmission"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Gets the radiation along the path by linear emission calculations.
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance_path", "spectral_radiance_path_jacobian"},

      .in = {"spectral_radiance_background",
             "spectral_radiance_path_source",
             "spectral_radiance_path_source_jacobian",
             "ppvar_tramat",
             "ppvar_cumtramat",
             "ppvar_dtramat"},

  };

  wsm_data["spectral_radiance_pathCalcTransmission"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Gets the radiation along the path by linear emission calculations.
)--",
          .author = {"Richard Larsson"},
          .out = {"spectral_radiance_path", "spectral_radiance_path_jacobian"},

          .in = {"ppvar_tramat", "ppvar_cumtramat", "ppvar_dtramat"},

      };

  wsm_data["spectral_radiance_path_sourceFromPropmat"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Gets the source term along the path.
)--",
          .author = {"Richard Larsson"},
          .out = {"spectral_radiance_path_source",
                  "spectral_radiance_path_source_jacobian"},

          .in = {"ppvar_propmat",
                 "ppvar_nlte",
                 "ppvar_dpropmat",
                 "ppvar_dnlte",
                 "ppvar_f",
                 "ppvar_atm",
                 "jacobian_targets"},

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
             "propagation_path",
             "ppvar_atm",
             "surface_field",
             "jacobian_targets"},
      .gin = {"hse_derivative"},
      .gin_type = {"Index"},
      .gin_value = {Index{0}},
      .gin_desc = {"Flag to compute the hypsometric distance derivatives"}};

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
             "jacobian_targets",
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
             "jacobian_targets",
             "atm_point",
             "path_point"},

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
             "jacobian_targets",
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
                 "f_grid",
                 "abs_species",
                 "select_abs_species",
                 "jacobian_targets",
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
             "jacobian_targets",
             "abs_lines_per_species",
             "atm_point",
             "nlte_vib_energies",
             "nlte_do"},
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
                 "f_grid",
                 "abs_species",
                 "select_abs_species",
                 "jacobian_targets",
                 "atm_point"},

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

A line-of-sight direction *path_point* is required as particles can
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
             "jacobian_targets",
             "path_point",
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
             "jacobian_targets",
             "f_grid",
             "atm_point"},
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
             "jacobian_targets",
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
             "jacobian_targets",
             "atm_point",
             "nlte_vib_energies",
             "path_point",
             "nlte_do"},
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

      .in = {"jacobian_targets", "f_grid", "propmat_clearsky_agenda_checked"},

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
to *pnd_field* using ``pnd_fieldCalcFromParticleBulkProps``.
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
to *pnd_field* using ``pnd_fieldCalcFromParticleBulkProps``.
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

  wsm_data["spectral_radiance_fieldPlaneParallelSpectralRadianceOperator"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Create a *spectral_radiance_field*

This is an experimental solution.
)--",
          .author = {"Richard Larsson"},
          .out = {"spectral_radiance_field"},

          .in = {"spectral_radiance_profile_operator", "f_grid", "za_grid"},

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

  wsm_data["sunsOff"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Turns all calculations with suns off 
)--",
      .author = {"Jon Petersen"},
      .out = {"suns_do", "suns"},

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

  wsm_data["water_p_eq_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *water_p_eq_agenda* to a default value

Options are:

- ``"MK05"``:
    1. Uses ``water_p_eq_fieldMK05`` to set *water_p_eq_field*
)--",
      .author = {"Richard Larsson"},
      .out = {"water_p_eq_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {String("MK05")},
      .gin_desc = {R"--(Default agenda option (see description))--"},
  };

  wsm_data["water_equivalent_pressure_operatorMK05"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Calculate equivalent water pressure according to Murphy and Koop, 2005

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
          .author = {"Patrick Eriksson", "Richard Larsson"},
          .out = {"water_equivalent_pressure_operator"},
          .gin = {"only_liquid"},
          .gin_type = {"Index"},
          .gin_value = {Index{0}},
          .gin_desc = {
              "Set to 1 to use liquid saturation pressure at all temperatures"}};

  wsm_data["atm_fieldHydrostaticPressure"] = {
      .desc = R"-x-(Add the hydrostatic pressure to the atmospheric field
    
The field must already be able to compute temperature as a function of
altitude, latitude, and longitude.

If it also contains species data, the species are used to compute the
average mass of the atmospheric molecules to get the specific gas
constant.  Note that this can also be overwritte with a positive
value for the equivalent GIN.

The ``alts`` vector contains the altitude grid values that limits the
extrapolation distance in altitude.  The first altitude in this
list should corresond to the altitude of the ``p0`` grid.  The extrapolation
outside of this range simply uses the hydrostatic equation
$P_1 = P_0 - g * h * \rho$ by means of the specific gas constant omputed
as desribed above and the pressure of the lower or first altitude level.
)-x-",
      .author = {"Richard Larsson"},
      .out = {"atm_field"},
      .in = {"atm_field", "gravity_operator"},
      .gin = {"p0",
              "alts",
              "fixed_specific_gas_constant",
              "fixed_atm_temperature",
              "hydrostatic_option"},
      .gin_type = {"GriddedField2", "Vector", "Numeric", "Numeric", "String"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    Numeric{-1},
                    Numeric{-1},
                    String{"HydrostaticEquation"}},
      .gin_desc = {
          "Lowest altitude pressure field",
          "Altitude vector",
          "Specific gas constant if larger than 0",
          "Constant atmospheric temprature if larger than 0",
          "Computational option for levels, [HydrostaticEquation, HypsometricEquation]"}};

  wsm_data["gravity_operatorFromGM"] = {
      .desc =
          R"-x-(Sets a gravity operator from the gravitational constant and the mass of the planet

Gets the ellispoid from *surface_field*
)-x-",
      .author = {"Richard Larsson"},
      .out = {"gravity_operator"},
      .in = {"surface_field"},
      .gin = {"GM"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {
          "Gravitation constant so that the gravity at radius ``r`` is ``GM / r^2``"}};

  wsm_data["spectral_radiance_backgroundAgendasAtEndOfPath"] = {
      .desc = R"--(Computes the background radiation.
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance_background",
              "spectral_radiance_background_jacobian"},
      .in = {"f_grid",
             "jacobian_targets",
             "propagation_path",
             "spectral_radiance_background_space_agenda",
             "spectral_radiance_background_surface_agenda"},
      .pass_workspace = true};

  wsm_data["spectral_radiance_backgroundUniformCosmicBackground"] = {
      .desc =
          R"--(Background spectral radiance is from a uniform cosmic background temperature.
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance_background"},
      .in = {"f_grid"}};

  wsm_data["spectral_radiance_backgroundSurfaceBlackbody"] = {
      .desc =
          R"--(Set surface spectral radiance from Planck function of the surface temperature
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance_background",
              "spectral_radiance_background_jacobian"},
      .in = {"f_grid", "surface_field", "jacobian_targets", "path_point"}};

  wsm_data["spectral_radiance_background_jacobianEmpty"] = {
      .desc = R"--(Set the cosmic background radiation derivative to empty.

Size : (*jacobian_targets*, *f_grid*)
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance_background_jacobian"},
      .in = {"f_grid", "jacobian_targets"}};

  wsm_data["spectral_radiance_jacobianEmpty"] = {
      .desc = R"--(Set the cosmic background radiation derivative to empty.

Size : (*jacobian_targets*, *f_grid*)
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance_jacobian"},
      .in = {"f_grid", "jacobian_targets"}};

  wsm_data["spectral_radiance_jacobianFromBackground"] = {
      .desc = R"--(Sets *spectral_radiance_jacobian* from the background values
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance_jacobian"},
      .in = {"spectral_radiance_background_jacobian",
             "background_transmittance"}};

  wsm_data["spectral_radiance_jacobianAddPathPropagation"] = {
      .desc = R"--(Adds the propagation variables to *spectral_radiance_jacobian*
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance_jacobian"},
      .in = {"spectral_radiance_jacobian",
             "spectral_radiance_path_jacobian",
             "jacobian_targets",
             "atm_field",
             "propagation_path"}};

  wsm_data["spectral_radianceFromPathPropagation"] = {
      .desc = R"--(Sets *spectral_radiance* from front of *spectral_radiance_path*
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance"},
      .in = {"spectral_radiance_path"}};

  wsm_data["spectral_radianceStandardEmission"] = {
      .desc = R"--(Sets *spectral_radiance* and *spectral_radiance_jacobian* from standard emission calculations
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance", "spectral_radiance_jacobian"},
      .in = {"f_grid",
             "jacobian_targets",
             "atm_field",
             "surface_field",
             "propagation_path",
             "spectral_radiance_background_space_agenda",
             "spectral_radiance_background_surface_agenda",
             "propmat_clearsky_agenda",
             "rte_alonglos_v"},
      .gin = {"hse_derivative"},
      .gin_type = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc =
          {"Whether or not hypsometric balance is assumed in temperature derivatives"},
      .pass_workspace = true};

  wsm_data["absorption_bandsFromAbsorbtionLines"] = {
      .desc = R"--(Gets modern line catalog from old style
)--",
      .author = {"Richard Larsson"},
      .out = {"absorption_bands"},
      .in = {"abs_lines_per_species"}};

  wsm_data["abs_linesFromAbsorptionBands"] = {
      .desc = R"--(Gets old line catalog from modern style
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines"},
      .in = {"absorption_bands"}};

  wsm_data["propmat_clearskyAddLines2"] = {
      .desc = R"--(Modern line-by-line calculations
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
             "jacobian_targets",
             "absorption_bands",
             "ecs_data2",
             "atm_point"},
      .gin = {"no_negative_absorption"},
      .gin_type = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc = {
          "Turn off to allow individual absorbers to have negative absorption"}};

  wsm_data["jacobian_targetsInit"] = {
      .desc = R"--(Initialize or reset the *jacobian_targets*
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_targets"}};

  wsm_data["jacobian_targetsFinalize"] = {
      .desc = R"--(Finalize *jacobian_targets* for use in RT methods
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_targets"},
      .in = {"jacobian_targets", "atm_field"}};

  wsm_data["jacobian_targetsAddTemperature"] = {
      .desc = R"--(Set temperature derivative
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_targets"},
      .in = {"jacobian_targets"},
      .gin = {"d"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{0.1}},
      .gin_desc = {
          "The perturbation used in methods that cannot compute derivatives analytically"}};

  wsm_data["jacobian_targetsAddPressure"] = {
      .desc = R"--(Set pressure derivative
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_targets"},
      .in = {"jacobian_targets"},
      .gin = {"d"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{0.1}},
      .gin_desc = {
          "The perturbation used in methods that cannot compute derivatives analytically"}};

  wsm_data["jacobian_targetsAddMagneticField"] = {
      .desc = R"--(Set magnetic field derivative
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_targets"},
      .in = {"jacobian_targets"},
      .gin = {"component", "d"},
      .gin_type = {"String", "Numeric"},
      .gin_value = {std::nullopt, Numeric{0.1}},
      .gin_desc = {
          "The component to use [u, v, w]",
          "The perturbation used in methods that cannot compute derivatives analytically"}};

  wsm_data["jacobian_targetsAddWindField"] = {
      .desc = R"--(Set wind field derivative

Note that the derivatives from methods that takes the freqeuncy will return
their derivatives as if these were frequency derivatives.
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_targets"},
      .in = {"jacobian_targets"},
      .gin = {"component", "d"},
      .gin_type = {"String", "Numeric"},
      .gin_value = {std::nullopt, Numeric{0.1}},
      .gin_desc = {
          "The component to use [u, v, w]",
          "The perturbation used in methods that cannot compute derivatives analytically"}};

  wsm_data["jacobian_targetsAddSpeciesVMR"] = {
      .desc = R"--(Set volume mixing ratio derivative
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_targets"},
      .in = {"jacobian_targets"},
      .gin = {"species", "d"},
      .gin_type = {"String", "Numeric"},
      .gin_value = {std::nullopt, Numeric{0.1}},
      .gin_desc = {
          "The species of interest (short or long name)",
          "The perturbation used in methods that cannot compute derivatives analytically"}};

  wsm_data["jacobian_targetsAddSpeciesIsotopologueRatio"] = {
      .desc = R"--(Set isotopologue ratio derivative
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_targets"},
      .in = {"jacobian_targets"},
      .gin = {"species", "d"},
      .gin_type = {"String", "Numeric"},
      .gin_value = {std::nullopt, Numeric{0.1}},
      .gin_desc = {
          "The species isotopologue of interest (short name)",
          "The perturbation used in methods that cannot compute derivatives analytically"}};

  wsm_data["jacobian_targetsAddLineParameter"] = {
      .desc = R"--(Set line parameter derivative

Note that empty ``coefficient`` means that the derivative is for a standard
line parameter (e.g., line center), otherwise it is for a line shape parameter.
By default, this variable is set to empty.
)--",
      .author = {"Richard Larsson"},
      .out = {"jacobian_targets"},
      .in = {"jacobian_targets", "absorption_bands"},
      .gin = {"id", "line_index", "parameter", "species", "coefficient", "d"},
      .gin_type = {"QuantumIdentifier",
                   "Index",
                   "String",
                   "String",
                   "String",
                   "Numeric"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    String{""},
                    String{""},
                    Numeric{0.1}},
      .gin_desc = {
          "The quantum identifier of the band",
          "The index of the line in the band",
          "The parameter to compute the derivative for (see options in error message)",
          "The coefficient in question (if non-empty, ``parameter`` refers to be a line shape parameter, otherwise, ``parameter`` referes to a standard absorption line parameter)",
          "The species of interest (long or short name; error message shows only valid long names)",
          "The perturbation used in methods that cannot compute derivatives analytically"}};

  wsm_data["absorption_bandsSelectFrequency"] = {
      .desc =
          R"--(Remove all bands that strictly falls outside a frequency range
)--",
      .author = {"Richard Larsson"},
      .out = {"absorption_bands"},
      .in = {"absorption_bands"},
      .gin = {"fmin", "fmax"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {-std::numeric_limits<Numeric>::infinity(),
                    std::numeric_limits<Numeric>::infinity()},
      .gin_desc = {"Minimum frequency to keep", "Maximum frequency to keep"}};

  wsm_data["absorption_bandsRemoveID"] = {.desc = R"--(Remove first band of ID
)--",
                                          .author = {"Richard Larsson"},
                                          .out = {"absorption_bands"},
                                          .in = {"absorption_bands"},
                                          .gin = {"id"},
                                          .gin_type = {"QuantumIdentifier"},
                                          .gin_value = {std::nullopt},
                                          .gin_desc = {"Identifier to remove"}};

  wsm_data["absorption_bandsKeepID"] = {
      .desc = R"--(Keeps first band of ID

If ``line`` is positive, also keep only the line of this index
)--",
      .author = {"Richard Larsson"},
      .out = {"absorption_bands"},
      .in = {"absorption_bands"},
      .gin = {"id", "line"},
      .gin_type = {"QuantumIdentifier", "Index"},
      .gin_value = {std::nullopt, Index{-1}},
      .gin_desc = {"Band to keep", "Line to keep (if positive)"}};

  wsm_data["SortedQuantumIdentifiersOfBands"] = {
      .desc =
          R"--(Get the sorting of the bands by first quantum identifier then some ``criteria``

The reverse sorting can also be achieved by setting ``reverse``.

Valid ``criteria`` are:
- None: No sorting after the quantum identifier sorting
- IntegratedIntensity: Sum of the intesities of the band at 296 K
- FrontFrequency: By first frequency
)--",
      .author = {"Richard Larsson"},
      .gout = {"sorted_idxs"},
      .gout_type = {"ArrayOfIndex"},
      .gout_desc = {"Sorted band indices (of *absorption_bands*)"},
      .in = {"absorption_bands"},
      .gin = {"criteria", "reverse"},
      .gin_type = {"String", "Index"},
      .gin_value = {String{"None"}, Index{0}},
      .gin_desc = {"Internal sorting criteria",
                   "Sort in reverse order if true"}};

  wsm_data["absorption_bandsReadSplit"] = {
      .desc = R"--(Saves all bands fin *absorption_bands* to a directory

This will create the directory if it does not exist.  It will also create
subdirectories that are the short-form of the isotopologue names.  The bands
will be stored as 0.xml, 1.xml, 2.xml, and so on

The ``dir`` path has to be absolute or relative to the working path, the environment
variables are not considered
)--",
      .author = {"Richard Larsson"},
      .out = {"absorption_bands"},
      .gin = {"dir"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {"Absolute or relative path to the directory"}};

  wsm_data["absorption_bandsSaveSplit"] = {
      .desc = R"--(Saves all bands fin *absorption_bands* to a directory

This will create the directory if it does not exist.  It will also create
subdirectories that are the short-form of the isotopologue names.  The bands
will be stored as 0.xml, 1.xml, 2.xml, and so on

The ``dir`` path has to be absolute or relative to the working path, the environment
variables are not considered
)--",
      .author = {"Richard Larsson"},
      .in = {"absorption_bands"},
      .gin = {"dir"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {"Absolute or relative path to the directory"}};

  wsm_data["propagation_pathGeometric"] = {
      .desc = R"--(Get a geometric radiation path

The path is defined by the origo and the line of sight.

The ``pos`` is either at the end or at the beginning of the path depending 
on the ``as_sensor`` flag.  A value that evaluates to true means that it is
at the end of the path.  If ``as_sensor`` is true, the ``los`` is therefore
looking backwards along the path.  Basically, ``as_sensor`` true means that
``pos`` and ``los`` behaves as sensor pos and los.

The ``max_step`` is the maximum step length in meters.  The path is first
created between the two extremes of either space and/or surface.  Afterwards,
there are additional points added every ``max_step`` meters between these
points until no more fits (the last step is shorter or exactly ``max_step``).

Upon closing the method, the following options are available to modify
the output:

If ``add_limb`` is true, the limb point is added to the path at the end.  It
is computed using bisections to ensure that the zenith angle of the tangent
point is as close to 90 degrees as it can numerically be.

If ``remove_non_atm`` is true, all points that are not in the atmosphere are
removed.  It is recommended to remove these points as multiple methods will
either perform poorly or not at all with these points present.

If ``fix_updown_azimuth`` is true, the azimuthal angle of the path is
fixed to the initial azimuthal angle of the path.  Because calculations
of the azimuth angle makes use of IEEE atan2, some paths may produce
bad angles if this is turned off.
)--",
      .author = {"Richard Larsson"},
      .out = {"propagation_path"},
      .in = {"atm_field", "surface_field"},
      .gin = {"pos",
              "los",
              "max_step",
              "as_sensor",
              "add_limb",
              "remove_non_atm",
              "fix_updown_azimuth"},
      .gin_type =
          {"Vector3", "Vector2", "Numeric", "Index", "Index", "Index", "Index"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    Numeric{1e3},
                    Index{1},
                    Index{0},
                    Index{1},
                    Index{1}},
      .gin_desc = {
          "The origo of the radiation path",
          "The line of sight of the radiation path",
          "The maximum step length",
          "Whether or not the path is as seen by the sensor or by the radiation (see text)",
          "Wheter or not to add the limb point",
          "Wheter or not to keep only atmospheric points",
          "Whether or not to attempt fix a potential issue with the path azimuthal angle"}};

  wsm_data["propagation_pathGeometricTangentAltitude"] = {
      .desc =
          R"--(Get a geometric radiation path that crosses the tangent altitude

The path is defined by an azimuth, a position, and a tangent altitude.
If the path ends up crossing the surface altitude, an error is thrown.

The ``pos`` is either at the end or at the beginning of the path depending 
on the ``as_sensor`` flag.  A value that evaluates to true means that it
is at the end of the path.  If ``as_sensor`` is true, the ``azimuth`` is
therefore looking backwards along the path.  Basically, ``as_sensor`` true
means that ``pos`` and ``azimuth`` behaves as sensor pos and azimuth.

The ``max_step`` is the maximum step length in meters.  The path is first
created between the two extremes of space and space.  Afterwards,
there are additional points added every ``max_step`` meters between these
points until no more fits (the last step is shorter or exactly ``max_step``).

Upon closing the method, the following options are available to modify
the output:

If ``add_limb`` is true, the limb point is added to the path at the end.  It
is computed using bisections to ensure that the zenith angle of the tangent
point is as close to 90 degrees as it can numerically be.

If ``remove_non_atm`` is true, all points that are not in the atmosphere are
removed.  It is recommended to remove these points as multiple methods will
either perform poorly or not at all with these points present.

If ``fix_updown_azimuth`` is true, the azimuthal angle of the path is
fixed to the initial azimuthal angle of the path.  Because calculations
of the azimuth angle makes use of IEEE atan2, some paths may produce
bad angles if this is turned off.
)--",
      .author = {"Richard Larsson"},
      .out = {"propagation_path"},
      .in = {"atm_field", "surface_field"},
      .gin = {"pos",
              "tangent_altitude",
              "azimuth",
              "max_step",
              "as_sensor",
              "add_limb",
              "remove_non_atm",
              "fix_updown_azimuth"},
      .gin_type = {"Vector3",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Index",
                   "Index",
                   "Index",
                   "Index"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    Numeric{1e3},
                    Index{1},
                    Index{1},
                    Index{1},
                    Index{1}},
      .gin_desc = {
          "The origo of the radiation path",
          "The tangent altitude of the radiation path",
          "The azimuth from the origo of the radiation path towards the tangent altitude",
          "The maximum step length",
          "Whether or not the path is as seen by the sensor or by the radiation (see text)",
          "Wheter or not to add the limb point",
          "Wheter or not to keep only atmospheric points",
          "Whether or not to attempt fix a potential issue with the path azimuthal angle"}};

  return wsm_data;
}
