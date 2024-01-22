#include "workspace_methods.h"

#include <limits>
#include <optional>

#include "nlte.h"

#pragma clang optimize off

using std::nullopt;

std::unordered_map<std::string, WorkspaceMethodInternalRecord>
internal_workspace_methods() {
  std::unordered_map<std::string, WorkspaceMethodInternalRecord> wsm_data;

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
      .in = {"f_grid"},
      .gin = {"surface_skin_t",
              "za",
              "salinity",
              "wind_speed",
              "rel_aa",
              "transmittance",
              "fastem_version"},
      .gin_type = {"Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Vector",
                   "Index"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    Numeric{0.035},
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    Index{6}},
      .gin_desc =
          {"Surface skin temperature [K]",
           R"--(Zenith angle of line-of-sigh, 90 to 180 deg.)--",
           R"--(Salinity, 0-1. That is, 3% is given as 0.03.)--",
           R"--(Wind speed.)--",
           R"--(Azimuth angle between wind direction and line-of-sight. This angle is measured clockwise from north, i.e. E=90deg.)--",
           R"--(The transmittance of the atmosphere, along the propagation path of the downwelling radiation. One value per frequency.)--",
           R"--(The version of FASTEM to use.)--"},

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

  wsm_data["surface_pointFromAtm"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Extracts a surface point from an atmospheric field

The elevation is from *path_point*, the normal points straight up, and
the atmosphere must contain the temperature.  If the wind is present,
it is also extracted.
)--",
      .author = {"Richard Larsson"},
      .gout = {"surface_point"},
      .gout_type = {"SurfacePoint"},
      .gout_desc = {R"--(Surface point.)--"},
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
      .gout = {"surface_point"},
      .gout_type = {"SurfacePoint"},
      .gout_desc = {R"--(Surface point.)--"},
      .in = {"path_point", "surface_field"},
      .gin = {"surface_search_accuracy"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{1.0}},
      .gin_desc = {R"--(Accuracy of surface search in meters.)--"},

  };

  wsm_data["surface_fieldSetPlanetEllipsoid"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets the planet base surface field

Options are:

- ``"Earth"``:

    1. Uses *surface_fieldEarth* with model="WGS84"

- ``"Io"``:

    1. Uses ``surface_fieldIo`` with model="Sphere"

- ``"Jupiter"``:

    1. Uses ``surface_fieldJupiter`` with model="Sphere"

- ``"Mars"``:

    1. Uses ``surface_fieldMars`` with model="Sphere"

- ``"Venus"``:

    1. Uses ``surface_fieldVenus`` with model="Sphere"
)--",
      .author = {"Richard Larsson"},
      .out = {"surface_field"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
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

  wsm_data["RadiationFieldSpectralIntegrate"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Integrates fields like ``spectral_irradiance_field`` or *spectral_radiance_field* over frequency.

Important, the first dimension must be the frequency dimension!
If a field like *spectral_radiance_field* is input, the stokes dimension
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
                   R"--(Normalization option)--",
                   R"--(Mirroring option)--",
                   R"--(Population option)--",
                   R"--(Lineshape option)--",
                   R"--(Cutoff option)--",
                   R"--(Cutoff value)--",
                   R"--(Line mixing limit)--"},

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
                   R"--(Normalization option)--",
                   R"--(Mirroring option)--",
                   R"--(Population option)--",
                   R"--(Lineshape option)--",
                   R"--(Cutoff option)--",
                   R"--(Cutoff value)--",
                   R"--(Line mixing limit)--"},

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
           R"--(Normalization option)--",
           R"--(Mirroring option)--",
           R"--(Population option)--",
           R"--(Lineshape option)--",
           R"--(Cutoff option)--",
           R"--(Cutoff value)--",
           R"--(Line mixing limit)--"},

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
                   R"--(Normalization option)--",
                   R"--(Mirroring option)--",
                   R"--(Population option)--",
                   R"--(Lineshape option)--",
                   R"--(Cutoff option)--",
                   R"--(Cutoff value)--",
                   R"--(Line mixing limit)--"},

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
                   R"--(Normalization option)--",
                   R"--(Mirroring option)--",
                   R"--(Population option)--",
                   R"--(Lineshape option)--",
                   R"--(Cutoff option)--",
                   R"--(Cutoff value)--",
                   R"--(Line mixing limit)--"},

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
           R"--(Normalization option)--",
           R"--(Mirroring option)--",
           R"--(Population option)--",
           R"--(Lineshape option)--",
           R"--(Cutoff option)--",
           R"--(Cutoff value)--",
           R"--(Line mixing limit)--"},

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
      .gin = {"file_index", "filename", "digits"},
      .gin_type = {"Index", "String", "Index"},
      .gin_value = {std::nullopt, String(""), Index{0}},
      .gin_desc = {
          R"--(Index of the file to read.)--",
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

  wsm_data["TMatrixTest"] = WorkspaceMethodInternalRecord{
      .desc = R"--(T-Matrix validation test.

Executes the standard test included with the T-Matrix Fortran code.
Should give the same as running the tmatrix_lp executable in
3rdparty/tmatrix/.
)--",
      .author = {"Oliver Lemke"},

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

  wsm_data["WMRFSelectChannels"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Select some channels for WMRF calculation.

The HIRS fast setup consists of a precalculated frequency grid
covering all HIRS channels, and associated weights for each channel,
stored in a weight matrix. (A ``sensor_response`` matrix.)

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

      .in = {"f_grid", "f_backend", "wmrf_weights"},
      .gin = {"wmrf_channels"},
      .gin_type = {"ArrayOfIndex"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(The channels to keep.)--"},
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

  wsm_data["WriteBuiltinPartitionFunctionsXML"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Writes all the builtin partition functions to file.

All available partition functions are written to files in the select format
in the select directory

The temperature will be linearly spaced between [Tlow, Tupp] with N values
)--",
      .author = {"Richard Larsson"},
      .gin = {"output_file_format", "dir", "Tlow", "Tupp", "N"},
      .gin_type = {"String", "String", "Numeric", "Numeric", "Index"},
      .gin_value = {String{"ascii"},
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt},
      .gin_desc = {"The format of the output",
                   R"--(The directory to write the data towards)--",
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

      .in = {},
      .gin = {"file_index", "input", "filename"},
      .gin_type =
          {"Index",
           "Vector, Matrix, Tensor3, Tensor4, Tensor5, ArrayOfVector, ArrayOfMatrix, GasAbsLookup",
           "String"},
      .gin_value = {std::nullopt, std::nullopt, String("")},
      .gin_desc = {R"--(Index number for files.)--",
                   R"--(Variable to be saved.)--",
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
      .gin = {"output_file_format", "input", "filename", "no_clobber"},
      .gin_type = {"String", "Any", "String", "Index"},
      .gin_value = {String("ascii"), std::nullopt, String(""), Index{0}},
      .gin_desc = {
          "The format of the output",
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
      .gin =
          {"output_file_format", "file_index", "input", "filename", "digits"},
      .gin_type = {"String", "Index", "Any", "String", "Index"},
      .gin_value =
          {String("ascii"), std::nullopt, std::nullopt, String(""), Index{0}},
      .gin_desc = {
          "The format of the output",
          R"--(Index number for files.)--",
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

  wsm_data["abs_linesCompact"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Removes lines that are unimportant because of their
cutoff frequency range
)--",
      .author = {"Stefan Buehler", "Richard Larsson"},
      .out = {"abs_lines"},

      .in = {"abs_lines", "f_grid"},

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

      .in = {"abs_lines"},
      .gin = {"output_file_format", "basename"},
      .gin_type = {"String", "String"},
      .gin_value = {String{"ascii"}, std::nullopt},
      .gin_desc = {"Format of the output file.",
                   R"--(Path to store the files at)--"},

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

  wsm_data["abs_lines_per_speciesPopulationNlteField"] =
      WorkspaceMethodInternalRecord{.desc = R"--(Turns on NTLE calculations.

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

                                    .in =
                                        {
                                            "abs_lines_per_species",
                                            "atm_field",
                                        },
                                    .gin = {"nlte_vib_energies"},
                                    .gin_type = {"VibrationalEnergyLevels"},
                                    .gin_value = {VibrationalEnergyLevels{}},
                                    .gin_desc = {"A map of energy levels"}

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

  wsm_data["abs_lines_per_speciesSetEmpty"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Empties *abs_lines_per_species* at the correct size.
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_lines_per_species"},

      .in = {"abs_species"},

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

          .in = {"abs_lines_per_species"},
          .gin = {"output_file_format", "basename"},
          .gin_type = {"String", "String"},
          .gin_value = {String{"ascii"}, std::nullopt},
          .gin_desc = {"Format of the output files.",
                       R"--(Path to store the files at)--"},

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
      .out = {"abs_species"},

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
      .out = {"abs_species"},

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
      .out = {"abs_species"},

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
)--",
      .author = {"Stefan Buehler"},
      .out = {"abs_species"},

      .gin = {"species"},
      .gin_type = {"ArrayOfString"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(Specify one String for each tag group that you want to create. Inside the String, separate the tags by commas (plus optional blanks).)--"},

  };

  wsm_data["antenna_responseGaussian"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets up a Gaussian antenna response.

This method works as *antenna_responseGaussianConstant* but allows
to inlude a frequency variation of the antenna width. Here the FWHM
is specified at a set of frequencies. These frequencies will also be
the frequency grid of ``antenna_response``.

If ``grid_width`` is set to <=0, the grid width will be twice the max
value in ``fwhm``.
)--",
      .author = {"Patrick Eriksson"},
      .gout = {"antenna_response"},
      .gout_type = {"NamedGriddedField3"},
      .gout_desc = {"The antenna pattern/response."},

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
      .gout = {"antenna_response"},
      .gout_type = {"NamedGriddedField3"},
      .gout_desc = {"The antenna pattern/response."},

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
range [``fstart``,``fstop``], with a logarithmic spacing.

The responses have a common angular grid. The parameters to define
the grid are the same as for *antenna_responseGaussianConstant*. If
``grid_width`` is <= 0, it is set to twice the FWHM at the lowest
frequency.
)--",
          .author = {"Patrick Eriksson"},
          .gout = {"antenna_response"},
          .gout_type = {"NamedGriddedField3"},
          .gout_desc = {"The antenna pattern/response."},

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

      .in = {"atm_field"},
      .gin = {"time", "parsafe"},
      .gin_type = {"Time", "Index"},
      .gin_value = {Time{}, Index{1}},
      .gin_desc = {
        "Time of data to use",
        R"--(Flag for parallel safety at 3X slowdown cost)--"},

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

  wsm_data["avkCalc"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculate the averaging kernel matrix.

This is done by describing the sensitivity of the
OEM retrieval with respect to the true state of the system. A prerequisite
for the calculation of the averaging kernel matrix is a successful OEM
calculation in which the ``jacobian`` and the gain matrix ``dxdy`` have been calculated.
)--",
      .author = {"Simon Pfreundschuh"},
      .out = {"avk"},

      .gin = {"dxdy", "jacobian"},
      .gin_type = {"Matrix", "Matrix"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Gain matrix of the OEM calculation.)--",
                   R"--(Jacobian matrix of the OEM calculation.)--"},

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
described by ``covmat_se``, which must be set by the user to include the
relevant contributions from the measurement and the forward model.

The ``covmat_se`` (``Se``) describes the uncertainty of the measurement vector (``y``),
and can be writtenn as::

  Se = Seps + Kb * Sb * Kb'

where Seps describes direct measurement errors (such as thermal noise),
Kb is Jacobian for forward model parameters, and Sb describes the uncertainty
of the forwatrd model parameters.

Prerequisite for the calculation of ``covmat_so`` is a successful OEM
computation where also the gain matrix has been computed. 
That is: So = G * Se * G', where G is the gain matrix (``dxdy``).
)--",
      .author = {"Simon Pfreundschuh"},
      .gout = {"covmat_so"},
      .gout_type = {"Matrix"},
      .gout_desc =
          {R"--(Covariance matrix describing the retrieval error due to uncertainties of the observation system.)--"},
      .gin = {"dxdy", "covmat_se"},
      .gin_type = {"Matrix", "CovarianceMatrix"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(Gain matrix.)--",
                   R"--(Covariance matrix for observation uncertainties.)--"},

  };

  wsm_data["covmat_ssCalc"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculates the covariance matrix describing the error due to smoothing.

The calculation of ``covmat_ss`` also requires the averaging kernel matrix *avk*
to be computed after a successful OEM calculation.

That is: Ss = (A-I) * Sx * (A-I)', where A is the averaging kernel 
matrix (*avk*).

The ``covmat_ss`` covariance matrix describes the Gaussian a priori distribution
for an OEM retrieval. It is represented using a symmetric block matrix.
covmat_sx can be used in two ways: Either with a block for each retrieval
quantity or with a single block containing the full covariance matrix.

Using a single block for each retrieval quantity has is advantageous for
if the retrieval quantities are assumed to be independent. In this case,
the covariance blocks can be added separately for each quantity and will
allow optimizing matrix multiplications and inverses required for the OEM
calculation.

The other case of using a single-block covariance matrix is supported
for convenience as well.
)--",
      .author = {"Simon Pfreundschuh"},
      .gout = {"covmat_ss"},
      .gout_type = {"Matrix"},
      .gout_desc =
          {R"--(Covariance matrix describing the retrieval error due to smoothing.)--"},
      .in = {"avk"},
      .gin = {"covmat_sx"},
      .gin_type = {"CovarianceMatrix"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Covariance matrix for state vector uncertainties.)--"},
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
the local oscillator frequencies ``lo_multi``, the backend
frequencies ``f_backend_multi``, and the backend channel
responses ``backend_channel_response_multi``.

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
      .gin = {"lo_multi",
              "spacing",
              "backend_channel_response_multi",
              "f_backend_multi"},
      .gin_type = {"Vector",
                   "Numeric",
                   "ArrayOfArrayOfGriddedField1",
                   "ArrayOfVector"},
      .gin_value = {std::nullopt, Numeric{.1e9}, std::nullopt, std::nullopt},
      .gin_desc =
          {"The local oscillator frequencies",
           R"--(Desired grid spacing in Hz.)--",
           "As *backend_channel_response* but describes an instrument with muliple mixer/receiver chains.",
           "As *f_backend* but describes an instrument with muliple mixer/receiver chains."},

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
      .gin = {"spacing",
              "verbosityVect",
              "backend_channel_response_multi",
              "f_backend_multi"},
      .gin_type = {"Numeric",
                   "Vector",
                   "ArrayOfArrayOfGriddedField1",
                   "ArrayOfVector"},
      .gin_value = {Numeric{.1e9}, Vector{}, std::nullopt, std::nullopt},
      .gin_desc =
          {R"--(Desired grid spacing in Hz.)--",
           R"--(Bandwidth adjusted spacing)--",
           "As *backend_channel_response* but describes an instrument with muliple mixer/receiver chains.",
           "As *f_backend* but describes an instrument with muliple mixer/receiver chains."},

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
``met_mm_backend`` table and method arguments.

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

``met_mm_backend`` is a compact description of a passband-type sensor, e.g. AMSU-A. The matrix
contains one row for each instrument channel. Each row contains four elements:

 - LO position [Hz]
 - first offset from the LO [Hz]
 - second offset from the LO+offset1 [Hz]
 - channel width [Hz]

Overview::

                            LO
                             |
                 offset1     |    offset1
             ----------------+----------------
             |                               |
             |                               |
    offset2  |  offset2             offset2  |  offset2
    ---------+---------             ---------+---------
    |                 |             |                 |
    |                 |             |                 |
  #####             #####         #####             #####
  width             width         width             width

For a sensor with 1 passband, offset1 and offset2 are zero.
For a sensor with 2 passbands, only offset2 is zero.
)--",
      .author = {"Oliver Lemke", "Patrick Eriksson"},
      .out = {"f_grid", "f_backend"},
      .gout = {"channel2fgrid_indexes", "channel2fgrid_weights"},
      .gout_type = {"ArrayOfArrayOfIndex", "ArrayOfVector"},
      .gout_desc =
          {"Definition of backend frequency response, link to *f_grid*.",
           "Definition of backend frequency response, weighting of *f_grid*."},
      .gin = {"met_mm_backend",
              "freq_spacing",
              "freq_number",
              "freq_merge_threshold"},
      .gin_type = {"Matrix", "Vector", "ArrayOfIndex", "Numeric"},
      .gin_value = {std::nullopt, Vector{.1e9}, ArrayOfIndex{-1}, Numeric{1}},
      .gin_desc =
          {"Backend description for meteorological millimeter sensors with passbands.",
           R"--(Desired grid spacing in Hz.)--",
           R"--(Number of frequencies per passband for each channel.)--",
           R"--(Merge frequencies that are closer than this value in Hz.)--"},

  };

  wsm_data["g0Earth"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gravity at zero altitude on Earth.

Sets *g0* for the given latitude using a standard parameterisation.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"g0"},

      .gin = {"lat"},
      .gin_type = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Latitude [degrees].)--"},

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

  wsm_data["heating_ratesFromIrradiance"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates heating rates from the *irradiance_field*.

The method assumes that the heating rates depend only on the
vertical derivation of the net flux. The net flux is the sum of the
*irradiance_field* in upward direction and the *irradiance_field*
in downward direction
)--",
      .author = {"Manfred Brath"},
      .gout = {"heating_rates"},
      .gout_type = {"Tensor3"},
      .gout_desc = {R"--(Heating rates)--"},
      .in = {"ppvar_atm", "irradiance_field", "g0"},
      .gin = {"specific_heat_capacity"},
      .gin_type = {"Tensor3"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Specific heat capacity)--"},
  };

  wsm_data["irradiance_fieldFromRadiance"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates the irradiance from the ``radiance_field``.

The ``radiance_field`` is integrated over the angular grids according to
the grids set by *AngularGridsSetFluxCalc*. 
See *AngularGridsSetFluxCalc* to set *za_grid*, *aa_grid*, and
*za_grid_weights*
)--",
      .author = {"Manfred Brath"},
      .out = {"irradiance_field"},
      .in = {"za_grid", "aa_grid", "za_grid_weights"},
      .gin = {"radiance_field"},
      .gin_type = {"Tensor5"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Radiance field)--"}};

  wsm_data["nlteOff"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Disable Non-LTE calculations.
)--",
      .author = {"Oliver Lemke"},
      .out = {"nlte_do"},

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
      .in = {"f_grid", "propagation_path", "ppvar_atm"},
      .gin = {"rte_alonglos_v"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{0.0}},
      .gin_desc =
          {R"--(Velocity along the line-of-sight to consider for a RT calculation.)--"},
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
      .out = {"ppvar_tramat", "ppvar_dtramat"},
      .gout = {"ppvar_distance", "ppvar_ddistance"},
      .gout_type = {"Vector", "ArrayOfArrayOfVector"},
      .gout_desc = {R"--(Distance between layers.)--",
                    R"--(Derivative of distance between layers.)--"},
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

The interpolation order in T and H2O is given by ``abs_t_interp_order``
and ``abs_nls_interp_order``, respectively.

Extraction is done for the frequencies in f_grid. Frequency
interpolation is controlled by ``abs_f_interp_order``. If this is zero,
then f_grid must either be the same as the internal frequency grid of
the lookup table (for efficiency reasons, only the first and last
element of f_grid are checked), or must have only a single element.
If ``abs_f_interp_order`` is above zero, then frequency is interpolated
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
             "f_grid",
             "atm_point",
             "jacobian_targets",
             "abs_species",
             "select_abs_species"},
      .gin =
          {
              "extpolfac",
              "no_negatives",
              "abs_p_interp_order",
              "abs_t_interp_order",
              "abs_nls_interp_order",
              "abs_f_interp_order",
          },
      .gin_type = {"Numeric", "Index", "Index", "Index", "Index", "Index"},
      .gin_value =
          {Numeric{0.5}, Index{1}, Index{5}, Index{7}, Index{5}, Index{0}},
      .gin_desc =
          {R"--(Extrapolation factor (for temperature and VMR grid edges).)--",
           R"--(Boolean. If it is true negative values due to interpolation are set to zero.)--",
           "The interpolation order to use when interpolating absorption between pressure levels.",
           "The interpolation order to use when interpolating absorption between the temperature values given by ``abs_t_pert``.",
           "The interpolation order to use when interpolating absorption between the H2O values given by ``abs_nls_pert``.",
           "Frequency interpolation order for absorption lookup table."},

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
             "nlte_do"},
      .gin = {"nlte_vib_energies",
              "lines_sparse_df",
              "lines_sparse_lim",
              "lines_speedup_option",
              "no_negatives"},
      .gin_type =
          {"VibrationalEnergyLevels", "Numeric", "Numeric", "String", "Index"},
      .gin_value = {VibrationalEnergyLevels{},
                    Numeric{0},
                    Numeric{0},
                    String("None"),
                    Index{1}},
      .gin_desc =
          {"A map of energy levels",
           R"--(The grid sparse separation)--",
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
             "path_point",
             "nlte_do"},
      .gin = {"nlte_vib_energies", "manual_mag_field", "H", "theta", "eta"},
      .gin_type =
          {"VibrationalEnergyLevels", "Index", "Numeric", "Numeric", "Numeric"},
      .gin_value = {VibrationalEnergyLevels{},
                    Index{0},
                    Numeric{1.0},
                    Numeric{0.0},
                    Numeric{0.0}},
      .gin_desc = {"A map of energy levels",
                   R"--(Manual angles tag)--",
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

      .in = {"jacobian_targets", "f_grid"},

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
      .out = {"propmat_clearsky_agenda"},
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

      .in = {"sensor_response_f", "sensor_response_f_grid"},
      .gin = {"lo", "sideband_mode"},
      .gin_type = {"Numeric", "String"},
      .gin_value = {nullopt, std::nullopt},
      .gin_desc = {R"--(Local oscillator frequency.)--",
                   R"("lower" or "upper" sideband mode.)"},

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
          .gout = {"spectral_irradiance_field"},
          .gout_type = {"Tensor5"},
          .gout_desc = {R"--(The spectral irradiance field)--"},
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

          .in = {"f_grid", "za_grid"},
          .gin = {"spectral_radiance_profile_operator"},
          .gin_type = {"SpectralRadianceProfileOperator"},
          .gin_value = {std::nullopt},
          .gin_desc = {R"--(The spectral radiance profile operator)--"},
      };

  wsm_data["sunsAddSingleBlackbody"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds a single blackbody to *suns*

Important note:
For a Sol-like sun there are huge differences in the UV-range 
between the actual sun spectrum and the blackbody spectrumwith the effective temperature of the sun. The blackbody sun"strongly overestimates the UV radiation.
)--",
      .author = {"Jon Petersen"},
      .out = {"suns"},

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
      .out = {"suns"},

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
      .out = {"suns"},
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

The two values of the reference ellipsoid are set manually. The two
arguments correspond directly to first and second element of
reference ellipsoid.
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

  wsm_data["SurfaceRadiationPropertyInterpFreq"] = WorkspaceMethodInternalRecord{
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
              R"--(Sets surface scalar reflectivity based on *surface_rmatrix*.

For each frequency f, surface scalar reflectivity is set to
the sum of surface_rmatrix(joker,f,0,0).
)--",
          .author = {"Patrick Eriksson"},
          .gout = {"surface_scalar_reflectivity"},
          .gout_type = {"Vector"},
          .gout_desc = {R"--(The scalar reflectivity)--"},
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

      .gin = {"lat", "lon", "atlas"},
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
          .gout = {"water_equivalent_pressure_operator"},
          .gout_type = {"NumericUnaryOperator"},
          .gout_desc = {"The water equivalent pressure operator."},
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

  wsm_data["gravity_operatorCentralMass"] = {
      .desc =
          R"-x-(Sets a gravity operator from the gravitational constant and the mass of the planet

Gets the ellispoid from *surface_field*
)-x-",
      .author = {"Richard Larsson"},
      .out = {"gravity_operator"},
      .in = {"surface_field"},
      .gin = {"mass"},
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
      .desc =
          R"--(Adds the propagation variables to *spectral_radiance_jacobian*
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance_jacobian"},
      .in = {"spectral_radiance_jacobian",
             "spectral_radiance_path_jacobian",
             "jacobian_targets",
             "atm_field",
             "propagation_path"}};

  wsm_data["spectral_radianceApplyUnit"] = {
      .desc = "Applies a unit to *spectral_radiance*, returning a new field\n",
      .author = {"Richard Larsson"},
      .gout = {"spectral_radiance_with_unit"},
      .gout_type = {"StokvecVector"},
      .gout_desc = {"The spectral radiance with unit"},
      .in = {"spectral_radiance", "f_grid"},
      .gin = {"spectral_radiance_unit"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {"The unit to apply"}};

  wsm_data["spectral_radianceFromPathPropagation"] = {
      .desc =
          R"--(Sets *spectral_radiance* from front of *spectral_radiance_path*
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance"},
      .in = {"spectral_radiance_path"}};

  wsm_data["spectral_radianceStandardEmission"] = {
      .desc =
          R"--(Sets *spectral_radiance* and *spectral_radiance_jacobian* from standard emission calculations
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
             "propmat_clearsky_agenda"},
      .gin = {"rte_alonglos_v", "hse_derivative"},
      .gin_type = {"Numeric", "Index"},
      .gin_value = {Numeric{0.0}, Index{1}},
      .gin_desc =
          {R"--(Velocity along the line-of-sight to consider for a RT calculation.)--",
           "Whether or not hypsometric balance is assumed in temperature derivatives"},
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
