#include "workspace_methods.h"

#include <limits>
#include <optional>

#pragma clang optimize off

using std::nullopt;

std::unordered_map<std::string, WorkspaceMethodInternalRecord>
internal_workspace_methods() {
  std::unordered_map<std::string, WorkspaceMethodInternalRecord> wsm_data;

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
in *absorption_species*.
)--",
      .author = {"Oliver Lemke"},
      .out = {"xsec_fit_data"},

      .in = {"absorption_species"},
      .gin = {"basename"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Basepath to the files)--"},
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

  wsm_data["propagation_matrix_cia_dataAddCIARecord"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Takes CIARecord as input and appends the results in the appropriate place.

If CIARecord has same species as species in *propagation_matrix_cia_data*, then the array
position is used to append all of the CIARecord into the array.  If clobber
evaluates as true, cia_record overwrites the appropriate *propagation_matrix_cia_data*.  If
species in cia_record are not in *propagation_matrix_cia_data*, the CIARecord is pushed back.
)--",
          .author = {"Richard Larsson"},
          .out = {"propagation_matrix_cia_data"},

          .in = {"propagation_matrix_cia_data"},
          .gin = {"cia_record", "clobber"},
          .gin_type = {"CIARecord", "Index"},
          .gin_value = {std::nullopt, Index{0}},
          .gin_desc =
              {R"--(CIA record to append to *propagation_matrix_cia_data*.)--",
               R"--(If true, the new input clobbers the old cia data.)--"},

      };

  wsm_data["propagation_matrix_cia_dataReadFromCIA"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Read data from a CIA data file for all CIA molecules defined
in *absorption_species*.

The units in the HITRAN file are:
 - Frequency: cm^(-1)
 - Binary absorption cross-section: cm^5 molec^(-2)

Upon reading we convert this to the ARTS internal SI units 
of Hz and m^5 molec^(-2).
)--",
          .author = {"Oliver Lemke"},
          .out = {"propagation_matrix_cia_data"},

          .in = {"absorption_species"},
          .gin = {"catalogpath"},
          .gin_type = {"String"},
          .gin_value = {std::nullopt},
          .gin_desc = {R"--(Path to the CIA catalog directory.)--"},

      };

  wsm_data["propagation_matrix_cia_dataReadFromXML"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Read data from a CIA XML file and check that all CIA tags defined
in *absorption_species* are present in the file.

The units of the data are described in *propagation_matrix_cia_dataReadFromCIA*.
)--",
          .author = {"Oliver Lemke"},
          .out = {"propagation_matrix_cia_data"},

          .in = {"absorption_species"},
          .gin = {"filename"},
          .gin_type = {"String"},
          .gin_value = {String("")},
          .gin_desc = {R"--(Name of the XML file.)--"},

      };

  wsm_data["propagation_matrix_cia_dataReadSpeciesSplitCatalog"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Reads a species split CIA dataset.
)--",
          .author = {"Richard Larsson"},
          .out = {"propagation_matrix_cia_data"},

          .in = {"absorption_species"},
          .gin = {"basename", "robust"},
          .gin_type = {"String", "Index"},
          .gin_value = {std::nullopt, Index{0}},
          .gin_desc =
              {R"--(The path to the split catalog files)--",
               R"--(Flag to continue in case nothing is found [0 throws, 1 continues])--"},
      };

  wsm_data["abs_lines_per_speciesReadSpeciesSplitCatalog"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Reads old style catalog but only for *absorption_species*
)--",
          .author = {"Richard Larsson"},
          .gout = {"abs_lines_per_species"},
          .gout_type = {"ArrayOfArrayOfAbsorptionLines"},
          .gout_desc = {R"--(Absorption lines per species)--"},

          .in = {"absorption_species"},

          .gin = {"basename", "robust"},
          .gin_type = {"String", "Index"},
          .gin_value = {std::nullopt, Index{0}},
          .gin_desc =
              {R"--(The path to the split catalog files)--",
               R"--(Flag to continue in case nothing is found [0 throws, 1 continues])--"},
      };

  wsm_data["absorption_bandsFromAbsorbtionLines"] = {
      .desc = R"--(Gets modern line catalog from old style
)--",
      .author = {"Richard Larsson"},
      .out = {"absorption_bands"},
      .in = {"absorption_species"},
      .gin = {"abs_lines_per_species"},
      .gin_type = {"ArrayOfArrayOfAbsorptionLines"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Absorption lines per species)--"},
  };

  wsm_data["propagation_matrix_absorption_lookupAdapt"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Adapts a gas absorption lookup table to the current calculation.

The lookup table can contain more species and more frequencies than
are needed for the current calculation. This method cuts down the
table in memory, so that it contains just what is needed. Also, the
species in the table are brought in the same order as the species in
the current calculation.

Of course, the method also performs quite a lot of checks on the
table. If something is not ok, a runtime error is thrown.
)--",
          .author = {"Stefan Buehler"},
          .out = {"propagation_matrix_absorption_lookup"},

          .in = {"propagation_matrix_absorption_lookup",
                 "absorption_species",
                 "frequency_grid"},

      };

  wsm_data["propagation_matrix_absorption_lookupInit"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Creates an empty gas absorption lookup table.

This is mainly there to help developers. For example, you can write
the empty table to an XML file, to see the file format.
)--",
          .author = {"Stefan Buehler"},
          .out = {"propagation_matrix_absorption_lookup"},

      };

  wsm_data["absorption_speciesAdd"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds species tag groups to the list of absorption species.

This WSM is similar to *absorption_speciesSet*, the only difference is that
this method appends species to an existing list of absorption species instead
of creating the whole list.

See *absorption_speciesSet* for details on how tags are defined and examples of
how to input them in the control file.
)--",
      .author = {"Stefan Buehler"},
      .out = {"absorption_species"},

      .in = {"absorption_species"},
      .gin = {"species"},
      .gin_type = {"ArrayOfString"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(Specify one String for each tag group that you want to add. Inside the String, separate the tags by commas (plus optional blanks).)--"},

  };

  wsm_data["absorption_speciesDefineAll"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *absorption_species* [i][0] to all species in ARTS
)--",
      .author = {"Richard Larsson"},
      .out = {"absorption_species"},

  };

  wsm_data["absorption_speciesDefineAllInScenario"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Define one tag group for each species known to ARTS and included in an
atmospheric scenario.

You can use this as an alternative to *absorption_speciesSet* if you want to make an
absorption calculation that is as complete as possible. The method
goes through all defined species and tries to open the VMR file. If
this works the tag is included, otherwise it is skipped.
)--",
      .author = {"Stefan Buehler"},
      .out = {"absorption_species"},

      .gin = {"basename"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(The name and path of a particular atmospheric scenario. For example: /pool/lookup2/arts-data/atmosphere/fascod/tropical)--"},

  };

  wsm_data["absorption_speciesInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets  *absorption_species* to be empty.
)--",
      .author = {"Stefan Buehler"},
      .out = {"absorption_species"},

  };

  wsm_data["absorption_speciesSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Set up a list of absorption species tag groups.

Workspace variables like *absorption_species* contain several tag
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
   absorption models see *propagation_matrixAddPredefined*

   Note that order of tag groups in the species list matters. In our
   example, changing the order of the first two tag group will give
   different results: as ``"O3"`` already selects all O3 transitions,
   no lines will remain to be selected by the
   ``"O3-666-500e9-501e9, O3-686"`` tag.

For CIA species the tag consists of the two involved species and
a dataset index. CIA species can be defined for multiple regions
The dataset index determines which region to use from the corresponding
CIARecord in *propagation_matrix_cia_data*.

Example

>>> species = [ "N2-CIA-N2-0, N2-CIA-N2-1" ]

For Hitran cross section species the tag consists of the species and
the tagtype XFIT, e.g. CFC11-XFIT. The data for the species must be
available in the *xsec_fit_data* variable.
)--",
      .author = {"Stefan Buehler"},
      .out = {"absorption_species"},

      .gin = {"species"},
      .gin_type = {"ArrayOfString"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(Specify one String for each tag group that you want to create. Inside the String, separate the tags by commas (plus optional blanks).)--"},
  };

  wsm_data["atmospheric_fieldAddCustomDataFile"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Add some custom data from file to the atmospheric_field

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
*atmospheric_fieldAddGriddedData*, or *atmospheric_fieldAddNumericData* to set the data.
Nevertheless this method is provided to make it easier to compose atmospheric
reading.
)--",
          .author = {"Richard Larsson"},
          .out = {"atmospheric_field"},

          .in = {"atmospheric_field"},
          .gin = {"key", "filename", "extrapolation_type"},
          .gin_type = {"String, ArrayOfSpeciesTag, QuantumIdentifier",
                       "String",
                       "String"},
          .gin_value = {std::nullopt, std::nullopt, String("Nearest")},
          .gin_desc = {R"--(Atmospheric data key.)--",
                       R"--(Filename)--",
                       R"--(Style of extrapolation)--"},

      };

  wsm_data["atmospheric_fieldAddField"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Add another atmospheric_field from file to the current atmospheric_field

The optional flag set_toa determines if the old (default) or
new (if it evaluates as true) atmospheric_field's top of the atmosphere altitude
is used in the output
)--",
      .author = {"Richard Larsson"},
      .out = {"atmospheric_field"},

      .in = {"atmospheric_field"},
      .gin = {"filename", "set_toa"},
      .gin_type = {"String", "Index"},
      .gin_value = {std::nullopt, Index{0}},
      .gin_desc = {R"--(Filename)--",
                   R"--(Flag for overwriting the top of the atmosphere)--"},

  };

  wsm_data["atmospheric_fieldAddGriddedData"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds data to the atmospheric_field

The field must not be regular
)--",
      .author = {"Richard Larsson"},
      .out = {"atmospheric_field"},

      .in = {"atmospheric_field"},
      .gin = {"key", "data", "extrapolation_type"},
      .gin_type = {"String, ArrayOfSpeciesTag, QuantumIdentifier",
                   "GriddedField3",
                   "String"},
      .gin_value = {std::nullopt, std::nullopt, String("Nearest")},
      .gin_desc = {R"--(See *atmospheric_fieldAddCustomDataFile*)--",
                   R"--(Some data)--",
                   R"--(Style of extrapolation)--"},

  };

  wsm_data["atmospheric_fieldAddNumericData"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Adds data to the atmospheric_field

The field must not be regular
)--",
      .author = {"Richard Larsson"},
      .out = {"atmospheric_field"},

      .in = {"atmospheric_field"},
      .gin = {"key", "data"},
      .gin_type = {"String, ArrayOfSpeciesTag, QuantumIdentifier", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(See *atmospheric_fieldAddCustomDataFile*)--",
                   R"--(Some data)--"},

  };

  wsm_data["atmospheric_fieldIGRF"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Use IGRF to compute the magnetic field at each point

The flag ``parsafe`` exists if you need the calculations to be safe in parallel
computations.
)--",
      .author = {"Richard Larsson"},
      .out = {"atmospheric_field"},

      .in = {"atmospheric_field"},
      .gin = {"time", "parsafe"},
      .gin_type = {"Time", "Index"},
      .gin_value = {Time{}, Index{1}},
      .gin_desc = {"Time of data to use",
                   R"--(Flag for parallel safety at 3X slowdown cost)--"},

  };

  wsm_data["atmospheric_fieldInit"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Initialize the atmospheric field with some altitude and isotopologue ratios
)--",
      .author = {"Richard Larsson"},
      .out = {"atmospheric_field"},

      .gin = {"toa", "default_isotopologue"},
      .gin_type = {"Numeric", "String"},
      .gin_value = {std::nullopt, String{"Builtin"}},
      .gin_desc = {R"--(Top of atmosphere altitude [m].)--",
                   "Default option for the isotopologue ratios"},
  };

  wsm_data["atmospheric_pointInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Initialize an atmospheric point with some isotopologue ratios
)--",
      .author = {"Richard Larsson"},
      .out = {"atmospheric_point"},

      .gin = {"default_isotopologue"},
      .gin_type = {"String"},
      .gin_value = {String{"Builtin"}},
      .gin_desc = {"Default option for the isotopologue ratios"},
  };

  wsm_data["atmospheric_fieldRead"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Reads a new atmospheric_field from a folder or base

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

If "read_specs" is true, then all the species of *absorption_species* are read and the
basename path is expected to contain a file with short-name version for each
unique species.  Some examples:

- absorption_species=["H2O-161", "O2-66"], - ["H2O.xml", "O2.xml"]
- absorption_species=["H2O-161", "O2-66", "CO2-626"], - ["H2O.xml", "O2.xml", "CO2.xml"]
- absorption_species=["H2O-161", "O2-66", "O2-PWR98"], - ["H2O.xml", "O2.xml"]
)--",
      .author = {"Richard Larsson"},
      .out = {"atmospheric_field"},

      .in = {"absorption_species"},
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

  wsm_data["atmospheric_fieldRescalePopulationLevels"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Rescale NLTE field to expected total distribution amongst levels
)--",
          .author = {"Richard Larsson"},
          .out = {"atmospheric_field"},

          .in = {"atmospheric_field"},
          .gin = {"s"},
          .gin_type = {"Numeric"},
          .gin_value = {std::nullopt},
          .gin_desc =
              {R"--(Scaling (e.g., 0.75 for only orth-water on Earth))--"},

      };

  wsm_data["atmospheric_fieldSave"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Saves an atmospheric_field to a folder or base

The output files are split to fit with what *atmospheric_fieldRead* can read, so see it
for most of the filenames that may be generated, depending on the content of the
atmospheric_field of course

Note that there are some exceptions.  If no_clobber is true, the new files will
not overwrite old files.  Also, if there are more than one species with the same
short-name, e.g., "H2O-161" and "H2O-181" both have short-name "H2O", only one of
these will print the "H2O.xml" file whereas the other will print "H2O.2.xml".
The latter is not read by *atmospheric_fieldRead*.  Even worse, we can give no guarantee
at all for whether it is the values from the "H2O-161" or "H2O-181" tags that
give the "H2O.xml" file because the internal data structure is unordered.
)--",
      .author = {"Richard Larsson"},

      .in = {"atmospheric_field"},
      .gin = {"basename", "filetype", "no_clobber"},
      .gin_type = {"String", "String", "Index"},
      .gin_value = {std::nullopt, String("ascii"), Index{0}},
      .gin_desc = {R"--(Base for the name of the data files.)--",
                   R"--(See *WriteXML*)--",
                   R"--(See *WriteXML*)--"},

  };

  wsm_data["background_transmittanceFromPathPropagationBack"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Sets *background_transmittance* to back of *propagation_path_transmission_matrix_cumulative*
)--",
          .author = {"Richard Larsson"},
          .out = {"background_transmittance"},

          .in = {"propagation_path_transmission_matrix_cumulative"},

      };

  wsm_data["background_transmittanceFromPathPropagationFront"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Sets *background_transmittance* to front of *propagation_path_transmission_matrix_cumulative*
)--",
          .author = {"Richard Larsson"},
          .out = {"background_transmittance"},

          .in = {"propagation_path_transmission_matrix_cumulative"},

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
      .in = {"ecs_data"},
  };

  wsm_data["ecs_dataAddTran2011"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the CO2-626, CO2-628, and CO2-636 band data for ECS.

Sets CO2 species
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data"},
      .in = {"ecs_data"},
  };

  wsm_data["ecs_dataAddRodrigues1997"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the CO2-626, CO2-628, and CO2-636 band data for ECS.

Sets N2 and O2 speces
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data"},
      .in = {"ecs_data"},
  };

  wsm_data["ecs_dataInit"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Resets/initializes the ECS data.
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data"},
  };

  wsm_data["ecs_dataAddMeanAir"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets ECS data for air from other data if available.
)--",
      .author = {"Richard Larsson"},
      .out = {"ecs_data"},

      .in = {"ecs_data"},

      .gin = {"vmrs", "species"},
      .gin_type = {"Vector", "ArrayOfSpecies"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc = {R"--(VMRs of air species)--", R"--(Air species)--"},
  };

  wsm_data["frequency_gridFromGasAbsLookup"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Sets *frequency_grid* to the frequency grid of *propagation_matrix_absorption_lookup*.

Must be called between importing/creating raw absorption table and
call of *propagation_matrix_absorption_lookupAdapt*.
)--",
      .author = {"Stefan Buehler"},
      .out = {"frequency_grid"},

      .in = {"propagation_matrix_absorption_lookup"},

  };

  wsm_data["propagation_path_atmospheric_pointFromPath"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Gets the atmospheric points along the path.
)--",
          .author = {"Richard Larsson"},
          .out = {"propagation_path_atmospheric_point"},
          .in = {"propagation_path", "atmospheric_field"},
      };

  wsm_data["propagation_path_transmission_matrix_cumulativeForward"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Sets *propagation_path_transmission_matrix_cumulative* by forward iteration of *propagation_path_transmission_matrix*
)--",
          .author = {"Richard Larsson"},
          .out = {"propagation_path_transmission_matrix_cumulative"},

          .in = {"propagation_path_transmission_matrix"},

      };

  wsm_data["propagation_path_transmission_matrix_cumulativeReverse"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Sets *propagation_path_transmission_matrix_cumulative* by reverse iteration of *propagation_path_transmission_matrix*
)--",
          .author = {"Richard Larsson"},
          .out = {"propagation_path_transmission_matrix_cumulative"},

          .in = {"propagation_path_transmission_matrix"},

      };

  wsm_data["propagation_path_frequency_gridFromPath"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Gets the frequency grid along the path.
)--",
      .author = {"Richard Larsson"},
      .out = {"propagation_path_frequency_grid"},
      .in = {"frequency_grid",
             "propagation_path",
             "propagation_path_atmospheric_point"},
      .gin = {"rte_alonglos_v"},
      .gin_type = {"Numeric"},
      .gin_value = {Numeric{0.0}},
      .gin_desc =
          {R"--(Velocity along the line-of-sight to consider for a RT calculation.)--"},
  };

  wsm_data["propagation_path_propagation_matrixCalc"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Gets the propagation matrix and NLTE source term along the path.
)--",
      .author = {"Richard Larsson"},
      .out =
          {"propagation_path_propagation_matrix",
           "propagation_path_propagation_matrix_source_vector_nonlte",
           "propagation_path_propagation_matrix_jacobian",
           "propagation_path_propagation_matrix_source_vector_nonlte_jacobian"},
      .in = {"propagation_matrix_agenda",
             "jacobian_targets",
             "propagation_path_frequency_grid",
             "propagation_path",
             "propagation_path_atmospheric_point"},
      .pass_workspace = true,
  };

  wsm_data["propagation_path_spectral_radianceCalcEmission"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Gets the radiation along the path by linear emission calculations.
)--",
          .author = {"Richard Larsson"},
          .out = {"propagation_path_spectral_radiance",
                  "propagation_path_spectral_radiance_jacobian"},

          .in = {"spectral_radiance_background",
                 "propagation_path_spectral_radiance_source",
                 "propagation_path_spectral_radiance_source_jacobian",
                 "propagation_path_transmission_matrix",
                 "propagation_path_transmission_matrix_cumulative",
                 "propagation_path_transmission_matrix_jacobian"},

      };

  wsm_data["propagation_path_spectral_radianceCalcTransmission"] =
      WorkspaceMethodInternalRecord{
          .desc =
              R"--(Gets the radiation along the path by linear emission calculations.
)--",
          .author = {"Richard Larsson"},
          .out = {"propagation_path_spectral_radiance",
                  "propagation_path_spectral_radiance_jacobian"},

          .in = {"propagation_path_transmission_matrix",
                 "propagation_path_transmission_matrix_cumulative",
                 "propagation_path_transmission_matrix_jacobian"},

      };

  wsm_data["propagation_path_spectral_radiance_sourceFromPropmat"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Gets the source term along the path.
)--",
          .author = {"Richard Larsson"},
          .out = {"propagation_path_spectral_radiance_source",
                  "propagation_path_spectral_radiance_source_jacobian"},

          .in =
              {"propagation_path_propagation_matrix",
               "propagation_path_propagation_matrix_source_vector_nonlte",
               "propagation_path_propagation_matrix_jacobian",
               "propagation_path_propagation_matrix_source_vector_nonlte_jacobian",
               "propagation_path_frequency_grid",
               "propagation_path_atmospheric_point",
               "jacobian_targets"},

      };

  wsm_data["propagation_path_transmission_matrixCalc"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Gets the transmission matrix in layers along the path.

A layer is defined as made up by the average of 2 levels, thus the outer-most size
of the derivatives out of this function is 2.
)--",
          .author = {"Richard Larsson"},
          .out = {"propagation_path_transmission_matrix",
                  "propagation_path_transmission_matrix_jacobian"},
          .gout = {"propagation_path_distance",
                   "propagation_path_distance_jacobian"},
          .gout_type = {"Vector", "ArrayOfArrayOfVector"},
          .gout_desc = {R"--(Distance between layers.)--",
                        R"--(Derivative of distance between layers.)--"},
          .in = {"propagation_path_propagation_matrix",
                 "propagation_path_propagation_matrix_jacobian",
                 "propagation_path",
                 "propagation_path_atmospheric_point",
                 "surface_field",
                 "jacobian_targets"},
          .gin = {"hse_derivative"},
          .gin_type = {"Index"},
          .gin_value = {Index{0}},
          .gin_desc = {"Flag to compute the hypsometric distance derivatives"}};

  wsm_data["propagation_matrix_predefined_model_dataAddWaterMTCKD400"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Sets the data for MT CKD 4.0 Water model

Note that the vectors must have the same length, and that wavenumbers must be growing
at a constant rate.  The minimum length is 4.

Note also that as this is predefined model data, the units of the values of the vectors
must be as described by each vector.
)--",
          .author = {"Richard Larsson"},
          .out = {"propagation_matrix_predefined_model_data"},

          .in = {"propagation_matrix_predefined_model_data"},
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

  wsm_data["propagation_matrix_predefined_model_dataInit"] =
      WorkspaceMethodInternalRecord{
          .desc = R"--(Initialize the predefined model data
)--",
          .author = {"Richard Larsson"},
          .out = {"propagation_matrix_predefined_model_data"},

      };

  wsm_data["propagation_matrixAddCIA"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculate absorption coefficients per tag group for HITRAN CIA continua.

This interpolates the cross sections from *propagation_matrix_cia_data*.

The robust option is intended only for testing. Do not use for normal
runs, since subsequent functions will not be able to deal with NAN values.
)--",
      .author = {"Stefan Buehler, Oliver Lemke"},
      .out = {"propagation_matrix", "propagation_matrix_jacobian"},

      .in = {"propagation_matrix",
             "propagation_matrix_jacobian",
             "absorption_species",
             "propagation_matrix_select_species",
             "jacobian_targets",
             "frequency_grid",
             "atmospheric_point",
             "propagation_matrix_cia_data"},
      .gin = {"T_extrapolfac", "ignore_errors"},
      .gin_type = {"Numeric", "Index"},
      .gin_value = {Numeric{0.5}, Index{0}},
      .gin_desc =
          {R"--(Temperature extrapolation factor (relative to grid spacing).)--",
           R"--(Set to 1 to suppress runtime errors (and return NAN values instead).)--"},

  };

  wsm_data["propagation_matrixAddFaraday"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Calculates absorption matrix describing Faraday rotation.

Faraday rotation is a change of polarization state of an
electromagnetic wave propagating through charged matter by
interaction with a magnetic field. Hence, this method requires
*absorption_species* to contain 'free_electrons' and electron content field
(as part of ``vmr_field``) as well as magnetic field (``mag_u_field``,
``mag_v_field``, ``mag_w_field``) to be specified.

Faraday rotation affects Stokes parameters 2 and 3 (but not
intensity!). Therefore, this method requires stokes_dim>2.

Like all 'propagation_matrixAdd*' methods, the method is additive,
i.e., does not overwrite the propagation matrix *propagation_matrix*,
but adds further contributions.
)--",
      .author = {"Patrick Eriksson"},
      .out = {"propagation_matrix", "propagation_matrix_jacobian"},

      .in = {"propagation_matrix",
             "propagation_matrix_jacobian",
             "frequency_grid",
             "absorption_species",
             "propagation_matrix_select_species",
             "jacobian_targets",
             "atmospheric_point",
             "propagation_path_point"},

  };

  wsm_data["propagation_matrixAddFromLookup"] = WorkspaceMethodInternalRecord{
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

Extraction is done for the frequencies in frequency_grid. Frequency
interpolation is controlled by ``abs_f_interp_order``. If this is zero,
then frequency_grid must either be the same as the internal frequency grid of
the lookup table (for efficiency reasons, only the first and last
element of frequency_grid are checked), or must have only a single element.
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
      .out = {"propagation_matrix", "propagation_matrix_jacobian"},

      .in = {"propagation_matrix",
             "propagation_matrix_jacobian",
             "propagation_matrix_absorption_lookup",
             "frequency_grid",
             "atmospheric_point",
             "jacobian_targets",
             "absorption_species",
             "propagation_matrix_select_species"},
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

  wsm_data["propagation_matrixAddPredefined"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Adds all of the predefined models in *absorption_species* to the propagation_matrix

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

  Note also that this model requires *propagation_matrix_predefined_model_data* to contain relevant data set either using
  *propagation_matrix_predefined_model_dataAddWaterMTCKD400* or via some file reading routine.

- H2O-ForeignContCKDMT400:
  Foreign continuum for water.  General reference: Mlawer et al. (2012), doi:10.1098/rsta.2011.0295

  Our code is reimplemented based on original Fortran90 code that is/was/will-be-made available via hitran.org

  Note that this model comes with the copyright statement [1].

  Note also that this model requires *propagation_matrix_predefined_model_data* to contain relevant data set either using
  *propagation_matrix_predefined_model_dataAddWaterMTCKD400* or via some file reading routine.

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
      .out = {"propagation_matrix", "propagation_matrix_jacobian"},

      .in = {"propagation_matrix",
             "propagation_matrix_jacobian",
             "propagation_matrix_predefined_model_data",
             "absorption_species",
             "propagation_matrix_select_species",
             "jacobian_targets",
             "frequency_grid",
             "atmospheric_point"},
  };

  wsm_data["propagation_matrixAddXsecFit"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Calculate absorption cross sections per tag group for HITRAN xsec species.

This broadens the cross section data from *xsec_fit_data* and
interpolates it onto the current frequency_grid.

Model data needs to be read in with *ReadXsecData* before calling
this method.
)--",
      .author = {"Oliver Lemke"},
      .out = {"propagation_matrix", "propagation_matrix_jacobian"},
      .in = {"propagation_matrix",
             "propagation_matrix_jacobian",
             "absorption_species",
             "propagation_matrix_select_species",
             "jacobian_targets",
             "frequency_grid",
             "atmospheric_point",
             "xsec_fit_data"},
      .gin = {"force_p", "force_t"},
      .gin_type = {"Numeric", "Numeric"},
      .gin_value = {Numeric{-1}, Numeric{-1}},
      .gin_desc = {R"--(Positive value forces constant pressure [Pa].)--",
                   R"--(Positive value forces constant temperature [K].)--"},
  };

  wsm_data["propagation_matrixInit"] = WorkspaceMethodInternalRecord{
      .desc =
          R"--(Initialize *propagation_matrix*, *propagation_matrix_source_vector_nonlte*, and their derivatives to zeroes.

This method must be used inside *propagation_matrix_agenda* and then be called first.
)--",
      .author = {"Oliver Lemke", "Richard Larsson"},
      .out = {"propagation_matrix",
              "propagation_matrix_source_vector_nonlte",
              "propagation_matrix_jacobian",
              "propagation_matrix_source_vector_nonlte_jacobian"},

      .in = {"jacobian_targets", "frequency_grid"},
  };

  wsm_data["propagation_matrix_agendaAuto"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets the *propagation_matrix_agenda* automatically

This method introspects the input and uses it for generating the
*propagation_matrix_agenda* automatically.  If ``use_propagation_matrix_absorption_lookup``, all
methods that can be used to generate the absorption lookup table
are ignored and instead the calculations from the absorption
lookup are used.

The following methods are considered for addition:
    1) *propagation_matrixInit*
    2) *propagation_matrixAddCIA*
    3) *propagation_matrixAddLines*
    5) *propagation_matrixAddFaraday*
    6) *propagation_matrixAddXsecFit*
    8) *propagation_matrixAddFromLookup*
    9) *propagation_matrixAddPredefined*

To perform absorption lookupo table calculation, call:
    1) *propagation_matrix_agendaAuto*
    2) ``propagation_matrix_absorption_lookupCalc``  FIXME: HOW TO COMPUTE IT
    3) *propagation_matrix_agendaAuto* (use_propagation_matrix_absorption_lookup=1)
    4) Perform other calculations
)--",
      .author = {"Richard Larsson"},
      .out = {"propagation_matrix_agenda"},
      .in = {"absorption_species", "absorption_bands"},
      .gin = {"T_extrapolfac",
              "extpolfac",
              "force_p",
              "force_t",
              "ignore_errors",
              "no_negatives",
              "use_propagation_matrix_absorption_lookup"},
      .gin_type = {"Numeric",
                   "Numeric",
                   "Numeric",
                   "Numeric",
                   "Index",
                   "Index",
                   "Index"},
      .gin_value = {Numeric{0.5},
                    Numeric{0.5},
                    Numeric{-1},
                    Numeric{-1},
                    Index{0},
                    Index{1},
                    Index{0}},
      .gin_desc =
          {R"--(See *propagation_matrixAddCIA*)--",
           R"--(See *propagation_matrixAddFromLookup*)--",
           R"--(See *propagation_matrixAddXsecFit*)--",
           R"--(See *propagation_matrixAddXsecFit*)--",
           R"--(See *propagation_matrixAddCIA*)--",
           R"--(See *propagation_matrixAddFromLookup*)--",
           R"--(Uses lookup calculations if true, ignores methods that can be part of the lookup table)--"},
  };

  wsm_data["propagation_matrix_agendaGUI"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Opens a GUI for running the propagation matrix agenda

Note that this is not thread-safe and should be executed on the main workspace

The values of all non-control flow are automatically loaded from the workspace
if they are defined.  Otherwise some values are just selected
)--",
      .author = {"Richard Larsson"},

      .in = {"propagation_matrix_agenda", "absorption_species"},
      .gin = {"load"},
      .gin_type = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc = {R"--(Load non-logical variables from workspace if true)--"},
      .pass_workspace = true,

  };

  wsm_data["propagation_matrix_agendaSet"] = WorkspaceMethodInternalRecord{
      .desc = R"--(Sets *propagation_matrix_agenda* to a default value

Please consider using *propagation_matrix_agendaAuto* instead of one of these options
as it will ensure you have the best coverage of use cases.  The options below are
available for feature testing

Options are:

- ``"Empty"``:

    1. Uses *propagation_matrixInit* to set *propagation_matrix*, *propagation_matrix_source_vector_nonlte*, *propagation_matrix_jacobian*, and *propagation_matrix_source_vector_nonlte_jacobian*
)--",
      .author = {"Richard Larsson"},
      .out = {"propagation_matrix_agenda"},

      .gin = {"option"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {R"--(Default agenda option (see description))--"},
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

  wsm_data["atmospheric_fieldHydrostaticPressure"] = {
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
      .out = {"atmospheric_field"},
      .in = {"atmospheric_field", "gravity_operator"},
      .gin = {"p0",
              "alts",
              "fixed_specific_gas_constant",
              "fixed_atmospheric_temperature",
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
      .in = {"frequency_grid",
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
      .in = {"frequency_grid"}};

  wsm_data["spectral_radiance_backgroundSurfaceBlackbody"] = {
      .desc =
          R"--(Set surface spectral radiance from Planck function of the surface temperature
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance_background",
              "spectral_radiance_background_jacobian"},
      .in = {"frequency_grid",
             "surface_field",
             "jacobian_targets",
             "propagation_path_point"}};

  wsm_data["spectral_radiance_background_jacobianEmpty"] = {
      .desc = R"--(Set the cosmic background radiation derivative to empty.

Size : (*jacobian_targets*, *frequency_grid*)
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance_background_jacobian"},
      .in = {"frequency_grid", "jacobian_targets"}};

  wsm_data["spectral_radiance_jacobianEmpty"] = {
      .desc = R"--(Set the cosmic background radiation derivative to empty.

Size : (*jacobian_targets*, *frequency_grid*)
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance_jacobian"},
      .in = {"frequency_grid", "jacobian_targets"}};

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
             "propagation_path_spectral_radiance_jacobian",
             "jacobian_targets",
             "atmospheric_field",
             "propagation_path"}};

  wsm_data["spectral_radianceApplyUnit"] = {
      .desc = "Applies a unit to *spectral_radiance*, returning a new field\n",
      .author = {"Richard Larsson"},
      .gout = {"spectral_radiance_with_unit"},
      .gout_type = {"StokvecVector"},
      .gout_desc = {"The spectral radiance with unit"},
      .in = {"spectral_radiance", "frequency_grid"},
      .gin = {"spectral_radiance_unit"},
      .gin_type = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc = {"The unit to apply"}};

  wsm_data["spectral_radianceFromPathPropagation"] = {
      .desc =
          R"--(Sets *spectral_radiance* from front of *propagation_path_spectral_radiance*
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance"},
      .in = {"propagation_path_spectral_radiance"}};

  wsm_data["spectral_radianceStandardEmission"] = {
      .desc =
          R"--(Sets *spectral_radiance* and *spectral_radiance_jacobian* from standard emission calculations
)--",
      .author = {"Richard Larsson"},
      .out = {"spectral_radiance", "spectral_radiance_jacobian"},
      .in = {"frequency_grid",
             "jacobian_targets",
             "atmospheric_field",
             "surface_field",
             "propagation_path",
             "spectral_radiance_background_space_agenda",
             "spectral_radiance_background_surface_agenda",
             "propagation_matrix_agenda"},
      .gin = {"rte_alonglos_v", "hse_derivative"},
      .gin_type = {"Numeric", "Index"},
      .gin_value = {Numeric{0.0}, Index{1}},
      .gin_desc =
          {R"--(Velocity along the line-of-sight to consider for a RT calculation.)--",
           "Whether or not hypsometric balance is assumed in temperature derivatives"},
      .pass_workspace = true};

  wsm_data["propagation_matrixAddLines"] = {
      .desc = R"--(Modern line-by-line calculations
)--",
      .author = {"Richard Larsson"},
      .out = {"propagation_matrix",
              "propagation_matrix_source_vector_nonlte",
              "propagation_matrix_jacobian",
              "propagation_matrix_source_vector_nonlte_jacobian"},
      .in = {"propagation_matrix",
             "propagation_matrix_source_vector_nonlte",
             "propagation_matrix_jacobian",
             "propagation_matrix_source_vector_nonlte_jacobian",
             "frequency_grid",
             "jacobian_targets",
             "propagation_matrix_select_species",
             "absorption_bands",
             "ecs_data",
             "atmospheric_point"},
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
      .in = {"jacobian_targets", "atmospheric_field"}};

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
      .in = {"atmospheric_field", "surface_field"},
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
      .in = {"atmospheric_field", "surface_field"},
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
