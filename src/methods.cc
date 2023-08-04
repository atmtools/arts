/*!
  \file   methods.cc
  \brief  Definition of method description data.

  This file contains only the definition of the function
  define_md_data, which sets the WSV lookup data. You have to change
  this function each time you add a new method. See methods.h for more
  documentation.

  \author Stefan Buehler
  \date 2000-06-10 */

#include "methods.h"
#include "arts.h"
#include "groups.h"
#include "wsv_aux.h"
#include <algorithm>
#include <array>
#include <string>

namespace global_data {
Array<MdRecord> md_data_raw;
extern const ArrayOfGroupRecord wsv_groups;
}  // namespace global_data

template <typename ... T>
std::array<const char *, sizeof...(T)> string_array(const T*... input)
{
  return {input...};
}

template <size_t NUM_OF_AUTHORS,
          size_t NUM_OF_OUTPUTS,
          size_t NUM_OF_GOUT_ARGS,
          size_t NUM_OF_GOUT_TYPES,
          size_t NUM_OF_GOUT_DESCRIPTIONS,
          size_t NUM_OF_INPUTS,
          size_t NUM_OF_GIN_ARGS,
          size_t NUM_OF_GIN_TYPES,
          size_t NUM_OF_GIN_DEFAULTS,
          size_t NUM_OF_GIN_DESCRIPTIONS,
          typename ... Ts>
MdRecord create_mdrecord(
  const char *name,
  const char *description,
  const std::array<const char *, NUM_OF_AUTHORS>& authors,
  const std::array<const char *, NUM_OF_OUTPUTS>& output,
  const std::array<const char *, NUM_OF_GOUT_ARGS>& gout,
  const std::array<const char *, NUM_OF_GOUT_TYPES>& gouttype,
  const std::array<const char *, NUM_OF_GOUT_DESCRIPTIONS>& goutdesc,
  const std::array<const char *, NUM_OF_INPUTS>& input,
  const std::array<const char *, NUM_OF_GIN_ARGS>& gin,
  const std::array<const char *, NUM_OF_GIN_TYPES>& gintype,
  const std::array<const char *, NUM_OF_GIN_DEFAULTS>& gindefault,
  const std::array<const char *, NUM_OF_GIN_DESCRIPTIONS>& gindesc,
  Ts ... flags)
{
  static_assert(NUM_OF_AUTHORS not_eq 0, "Must have at least one author");
  static_assert(NUM_OF_GOUT_ARGS == NUM_OF_GOUT_TYPES, "GOUT type(s) count does not match number of GOUT");
  static_assert(NUM_OF_GOUT_ARGS == NUM_OF_GOUT_DESCRIPTIONS, "GOUT description(s) count does not match number of GOUT");
  static_assert(NUM_OF_GIN_ARGS == NUM_OF_GIN_TYPES, "GIN type(s) count does not match number of GIN");
  static_assert(NUM_OF_GIN_ARGS == NUM_OF_GIN_DEFAULTS, "GIN default(s) count does not match number of GIN");
  static_assert(NUM_OF_GIN_ARGS == NUM_OF_GIN_DESCRIPTIONS, "GIN description(s) count does not match number of GIN");

  return MdRecord(name,
                  description,
                  ArrayOfString(authors),
                  ArrayOfString(output),
                  ArrayOfString(gout),
                  ArrayOfString(gouttype),
                  ArrayOfString(goutdesc),
                  ArrayOfString(input),
                  ArrayOfString(gin),
                  ArrayOfString(gintype),
                  ArrayOfString(gindefault),
                  ArrayOfString(gindesc),
                  bool(flags) ...);
}

// Some #defines and typedefs to make the records better readable:
#define NAME(x) x
#define DESCRIPTION(x) x
#define AUTHORS(...) \
  string_array( __VA_ARGS__ )
#define OUT(...) \
  string_array( __VA_ARGS__ )
#define GOUT(...) \
  string_array( __VA_ARGS__ )
#define GOUT_TYPE(...) \
  string_array( __VA_ARGS__ )
#define GOUT_DESC(...) \
  string_array( __VA_ARGS__ )
#define IN(...) \
  string_array( __VA_ARGS__ )
#define GIN(...) \
  string_array( __VA_ARGS__ )
#define GIN_TYPE(...) \
  string_array( __VA_ARGS__ )
#define GIN_DEFAULT(...) \
  string_array( __VA_ARGS__ )
#define GIN_DESC(...) \
  string_array( __VA_ARGS__ )
#define SETMETHOD(x) x
#define AGENDAMETHOD(x) x
#define USES_TEMPLATES(x) x
#define PASSWORKSPACE(x) x
#define PASSWSVNAMES(x) x

/* Here's a template record entry:  (PE 2008-09-20)

  md_data_raw.push_back
    ( create_mdrecord
      ( NAME( "MethodName" ),
        DESCRIPTION(R"--(A concise summary of the method.

A more detailed description of the method. Try to describe the
purpose of the method and important considerations. Try to avoid
references to other WSMs as they might change. Refer to the user
guide for more complex information (as long as it exists, or that
you add it to AUG!).

You do not need to describe workspace variables used. That
information is found in workspace.cc. Generic
output and input variables must be described in GIN_DESC and
GOUT_DESC below.
)--"),
        AUTHORS( "Your Name" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "descriptive_name_for_generic_input1" ),
        GIN_TYPE(    "GenericInput1Type" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC(    "Description for Generic Input Variable 1" )
        ));

 For variable descriptions longer than one line, use the following format.
 Don't forget to remove the space in '/ *' and '* /' if you copy this template.
 I had to put it in there because C++ doesn't allow nested comments.

  md_data_raw.push_back
    ( create_mdrecord
      ( NAME( "MethodName" ),
        ...
        ...
        ...
        GIN( gin_var1, gin_var2, gin_var3 )
        GIN_TYPE( "GInput1Type", "GInput2Type", "GInput3Type" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( / * gin_var1 * /
                  "Long description for Generic Input Variable 1 "
                  "which can span multiple lines like for example "
                  "this one. Don't put any \n in the variable descriptions.",
                  / * gin_var2 * /
                  "Long description for Generic Input Variable 2 "
                  "which can span multiple lines like for example "
                  "this one. Don't put any \n in the variable descriptions.",
                  / * gin_var3 * /
                  "Long description for Generic Input Variable 3 "
                  "which can span multiple lines like for example "
                  "this one. Don't put any \n in the variable descriptions."
                  )

*/

void define_md_data_raw() {
  // The variable md_data is defined in file methods_aux.cc.
  using global_data::md_data_raw;

  // Initialise to zero, just in case:
  md_data_raw.resize(0);

  // String with all array groups
  const String ARRAY_GROUPS = get_array_groups_as_string();
  // String with all groups that also exist as an array type
  const String GROUPS_WITH_ARRAY_TYPE = get_array_groups_as_string(true, true);
  // String with all array types whose element type is also available as a group
  const String ARRAY_GROUPS_WITH_BASETYPE =
      get_array_groups_as_string(true, false);

  using global_data::wsv_groups;

  for (const auto & wsv_group : wsv_groups) {
    if (wsv_group != "Any") {
      md_data_raw.push_back(MdRecord(
          NAME(String(wsv_group.name + "Create").c_str()),
          DESCRIPTION(
              String("Creates a variable of group " + wsv_group.name +
                     ".\n"
                     "\n"
                     "After being created, the variable is uninitialized.\n")
                  .c_str()),
          ArrayOfString(AUTHORS("Oliver Lemke")),
          ArrayOfString(OUT()),
          ArrayOfString(GOUT("output")),
          ArrayOfString(GOUT_TYPE(wsv_group.name.c_str())),
          ArrayOfString(GOUT_DESC("Variable to create.")),
          ArrayOfString(IN()),
          ArrayOfString(GIN()),
          ArrayOfString(GIN_TYPE()),
          ArrayOfString(GIN_DEFAULT()),
          ArrayOfString(GIN_DESC()),
          SETMETHOD(false),
          AGENDAMETHOD(false),
          USES_TEMPLATES(false),
          PASSWORKSPACE(wsv_group == "Agenda"),
          PASSWSVNAMES(false)));

      if (wsv_group not_eq "Agenda" and wsv_group not_eq "ArrayOfAgenda") {
        md_data_raw.push_back(MdRecord(
            NAME(String(wsv_group.name + "Set").c_str()),
            DESCRIPTION(R"--(Sets a workspace variable to the given value.
)--"),
            ArrayOfString(AUTHORS("Richard Larsson")),
            ArrayOfString(OUT()),
            ArrayOfString(GOUT("output")),
            ArrayOfString(GOUT_TYPE(wsv_group.name.c_str())),
            ArrayOfString(GOUT_DESC("Variable to initialize.")),
            ArrayOfString(IN()),
            ArrayOfString(GIN("value")),
            ArrayOfString(GIN_TYPE(wsv_group.name.c_str())),
            ArrayOfString(GIN_DEFAULT(NODEF)),
            ArrayOfString(GIN_DESC("The value.")),
            SETMETHOD(true)));
        }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // Let's put in the functions in alphabetical order. This gives a clear rule
  // for where to place a new function and this gives a nicer results when
  // the functions are listed by "arts -m all".
  // No distinction is made between uppercase and lowercase letters.
  // Patrick Eriksson 2002-05-08
  /////////////////////////////////////////////////////////////////////////////

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_cia_dataAddCIARecord"),
      DESCRIPTION(R"--(Takes CIARecord as input and appends the results in the appropriate place.

If CIARecord has same species as species in *abs_cia_data*, then the array
position is used to append all of the CIARecord into the array.  If clobber
evaluates as true, cia_record overwrites the appropriate *abs_cia_data*.  If
species in cia_record are not in *abs_cia_data*, the CIARecord is pushed back.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_cia_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_cia_data"),
      GIN("cia_record", "clobber"),
      GIN_TYPE("CIARecord", "Index"),
      GIN_DEFAULT(NODEF, "0"),
      GIN_DESC("CIA record to append to *abs_cia_data*.",
               "If true, the new input clobbers the old cia data.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_cia_dataReadSpeciesSplitCatalog"),
      DESCRIPTION(R"--(Reads a species split CIA dataset.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_cia_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("basename", "robust"),
      GIN_TYPE("String", "Index"),
      GIN_DEFAULT(NODEF, "0"),
      GIN_DESC("The path to the split catalog files",
               "Flag to continue in case nothing is found [0 throws, 1 continues]")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_cia_dataReadFromCIA"),
      DESCRIPTION(R"--(Read data from a CIA data file for all CIA molecules defined
in *abs_species*.

The units in the HITRAN file are:
 - Frequency: cm^(-1)
 - Binary absorption cross-section: cm^5 molec^(-2)

Upon reading we convert this to the ARTS internal SI units 
of Hz and m^5 molec^(-2).
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("abs_cia_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("catalogpath"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Path to the CIA catalog directory.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_cia_dataReadFromXML"),
      DESCRIPTION(R"--(Read data from a CIA XML file and check that all CIA tags defined
in *abs_species* are present in the file.

The units of the data are described in *abs_cia_dataReadFromCIA*.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("abs_cia_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("filename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(""),
      GIN_DESC("Name of the XML file.")));

  md_data_raw.push_back(create_mdrecord(
    NAME("abs_hitran_relmat_dataReadHitranRelmatDataAndLines"),
      DESCRIPTION(R"--(Reads HITRAN line mixing data from a basedir

The basedir must point at line mixing data as provided by HITRAN.
The lines will be changed such that ALL CO2 lines are truncated
before adding the HITRAN line mixing lines.

The available modes are such that \"VP*\" uses Voigt profiles and
\"SDVP*\" uses speed-dependent Voigt profiles, where the \"_Y\"
signifies if Rosenkranz-style line mixing is considered or not, and
the \"W\" at the end signifies that full calculations are used.  At
the line mixing limit, line mixing is simply turned off.

The \"FullW\" mode uses Lorentzian calculations with the full relaxation
matrix until the line mixing limit is reached and it switches to Voigt.

The HITRAN LM data is available for download at:
https:
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_hitran_relmat_data", "abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species", "abs_species"),
      GIN("basedir", "linemixinglimit", "fmin", "fmax", "stot", "mode"),
      GIN_TYPE("String", "Numeric", "Numeric", "Numeric", "Numeric", "String"),
      GIN_DEFAULT(NODEF, "-1", "-1e99", "1e99", "0", "VP_W"),
      GIN_DESC("Direcory where the linemixing data is to be found",
               "Line mixing limit as defined by *AbsorptionLines*",
               "Minimum frequency to read from",
               "Maximum frequency to read until",
               "Minimum integrated band strength to consider",
               "Mode of calculations.  The options are: \"VP\", \"VP_Y\", \"SDVP\", \"SDVP_Y\", \"FullW\", and \"VP_W\""
              )));

  md_data_raw.push_back(create_mdrecord(
    NAME("abs_linesAdaptOnTheFlyLineMixing"),
      DESCRIPTION(R"--(Adapts the line-catalog from using *ecs_data* data to.
instead fit ordered parameters to imitate the line mxixing

The order should be 1 or 2.  It will compute at 3 as well, but
there's no support in current ARTS LBL to make use of it so it
will crash at some point
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines", "ecs_data"),
      GIN("t_grid", "pressure", "order", "robust", "rosenkranz_adaptation"),
      GIN_TYPE("Vector", "Numeric", "Index", "Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "1", "0"),
      GIN_DESC("The sorted temperature grid",
               "The pressure at which the adaptation is made",
               "The order of the parameters in adaptation",
               "Boolean for failed band adaptation behavior. 0: throw exception. not 0: conversion to line-by-line calculations",
               "Apply direct Rosenkranz adaptation instead of computing the Eigenvalues")));

  md_data_raw.push_back(create_mdrecord(
    NAME("abs_lines_per_speciesAdaptOnTheFlyLineMixing"),
      DESCRIPTION(R"--(Calls *abs_linesAdaptOnTheFlyLineMixing* for each internal array
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species", "ecs_data"),
      GIN("t_grid", "pressure", "order", "robust", "rosenkranz_adaptation"),
      GIN_TYPE("Vector", "Numeric", "Index", "Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "1", "0"),
      GIN_DESC("The sorted temperature grid",
               "The pressure at which the adaptation is made",
               "The order of the parameters in adaptation",
               "Boolean for failed band adaptation behavior. 0: throw exception. not 0: conversion to line-by-line calculations",
               "Apply direct Rosenkranz adaptation instead of computing the Eigenvalues")));

  md_data_raw.push_back(create_mdrecord(
    NAME("abs_lines_per_speciesAdaptHitranLineMixing"),
      DESCRIPTION(R"--(Adapts the line-catalog from using *abs_hitran_relmat_data* to.
instead fit ordered parameters to imitate the line mxixing

The order should be 1 or 2.  It will compute at 3 as well, but
there's no support in current ARTS LBL to make use of it so it
will crash at some point
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species", "abs_hitran_relmat_data"),
      GIN("t_grid", "pressure", "order"),
      GIN_TYPE("Vector", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("The sorted temperature grid",
               "The pressure at which the adaptation is made",
               "The order of the parameters in adaptation")));

  md_data_raw.push_back(create_mdrecord(
    NAME("abs_linesKeepBand"),
      DESCRIPTION(R"--(Keep only ``qid``-match band lines in *abs_lines*

Note that other bands are technically kept but have zero lines
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("qid"),
      GIN_TYPE("QuantumIdentifier"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Band ID")));

  md_data_raw.push_back(create_mdrecord(
    NAME("abs_linesRemoveBand"),
      DESCRIPTION(R"--(Removes ``qid`` band from *abs_lines*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("qid"),
      GIN_TYPE("QuantumIdentifier"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Band ID")));

  md_data_raw.push_back(create_mdrecord(
    NAME("abs_linesRemoveLines"),
      DESCRIPTION(R"--(Remove lines *abs_lines* outside of specifications

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
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("lower_frequency", "upper_frequency", "lower_intensity", "safe", "flip_flims"),
      GIN_TYPE("Numeric", "Numeric", "Numeric", "Index", "Index"),
      GIN_DEFAULT("-1e99", "1e99", "0", "1", "0"),
      GIN_DESC("The lower frequency bound",
        "The upper frequency bound",
        "The lower intensity bound",
        "Remove only lines from a band if all lines of a band fail",
        "Reverse the frequecy filtering, see above")));

  md_data_raw.push_back(create_mdrecord(
    NAME("abs_lines_per_speciesRemoveLines"),
      DESCRIPTION(R"--(Repeats *abs_linesRemoveLines* for all inner arrays
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN("lower_frequency", "upper_frequency", "lower_intensity", "safe", "flip_flims"),
      GIN_TYPE("Numeric", "Numeric", "Numeric", "Index", "Index"),
      GIN_DEFAULT("-1e99", "1e99", "0", "1", "0"),
      GIN_DESC("The lower frequency bound",
        "The upper frequency bound",
        "The lower intensity bound",
        "Remove only lines from a band if all lines of a band fail",
        "Reverse the frequecy filtering")));

  md_data_raw.push_back(create_mdrecord(
    NAME("abs_linesRemoveLinesFromSpecies"),
      DESCRIPTION(R"--(As *abs_linesRemoveLines* but only for bands of the given species.

``species`` must be a single entry, and must specify the isotopologue
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("species", "lower_frequency", "upper_frequency", "lower_intensity",
          "safe", "flip_flims"),
      GIN_TYPE("ArrayOfSpeciesTag", "Numeric", "Numeric", "Numeric", "Index", "Index"),
      GIN_DEFAULT(NODEF, "-1e99", "1e99", "0", "1", "0"),
      GIN_DESC("Species to be removed",
        "The lower frequency bound",
        "The upper frequency bound",
        "The lower intensity bound",
        "Remove only lines from a band if all lines of a band fail",
        "Reverse the frequecy filtering")));

  md_data_raw.push_back(create_mdrecord(
    NAME("abs_lines_per_speciesRemoveLinesFromSpecies"),
      DESCRIPTION(R"--(Repeats *abs_linesRemoveLinesFromSpecies* for all inner arrays
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN("species", "lower_frequency", "upper_frequency", "lower_intensity",
          "safe", "flip_flims"),
      GIN_TYPE("ArrayOfSpeciesTag", "Numeric", "Numeric", "Numeric", "Index", "Index"),
      GIN_DEFAULT(NODEF, "-1e99", "1e99", "0", "1", "0"),
      GIN_DESC("Species to be removed",
        "The lower frequency bound",
        "The upper frequency bound",
        "The lower intensity bound",
        "Remove only lines from a band if all lines of a band fail",
        "Reverse the frequecy filtering")));

  md_data_raw.push_back(create_mdrecord(
    NAME("abs_linesRemoveEmptyBands"),
      DESCRIPTION(R"--(Removes emtpy bands from *abs_lines*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
    NAME("abs_linesFlatten"),
      DESCRIPTION(R"--(Makes *abs_lines* with the same ID share lines
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
    NAME("abs_lines_per_speciesFlatten"),
      DESCRIPTION(R"--(Calls *abs_linesFlatten* per internal set of bands
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("CallbackFunctionExecute"),
      DESCRIPTION(R"--(Execute any code in Arts
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("function"),
      GIN_TYPE("CallbackFunction"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("This will execute as \"function(current workspace);\""),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true),
      PASSWSVNAMES(false)));

  md_data_raw.push_back(create_mdrecord(
      NAME("CheckUnique"),
      DESCRIPTION(R"--(Checks that *abs_lines* contains only unique absorption lines
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_linesReplaceLines"),
      DESCRIPTION(R"--(Replace all lines in *abs_lines* that match with lines in replacement_lines.

Each replacement_lines must match excatly a single line in *abs_lines*.

The matching requires identical quantum number signatures to work

Note that lines are identified by their quantum number identifier, and if the broadening or
compute data disagree between two bands, a new band is appended unless we can work around the issue.
This may cause *CheckUnique* to fail after running this method
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("replacing_lines"),
      GIN_TYPE("ArrayOfAbsorptionLines"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Line-array that replace lines in *abs_lines*.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_linesReplaceBands"),
      DESCRIPTION(R"--(Replace all bands in *abs_lines* that match with bands in ``replacing_bands``.

Each ``replacing_bands`` must match excatly a single band in *abs_lines*.

The matching requires identical quantum number signatures to work.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("replacing_bands"),
      GIN_TYPE("ArrayOfAbsorptionLines"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Line-array that removes lines from *abs_lines*.")));

  md_data_raw.push_back(create_mdrecord(
    NAME("abs_linesDeleteBadF0"),
      DESCRIPTION(R"--(Deletes all lines in *abs_lines* that have bad central frequencies

If lower evaluates as true, deletes all lines with a frequency below f0.
Otherwise deletes all lines with a frequency above f0.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("f0", "lower"),
      GIN_TYPE("Numeric", "Index"),
      GIN_DEFAULT(NODEF, "1"),
      GIN_DESC("Target frequency",
               "Lower or upper flag (eval as boolean)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_linesReadSpeciesSplitCatalog"),
      DESCRIPTION(R"--(Reads a catalog of absorption lines files in a directory
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("basename", "robust"),
      GIN_TYPE("String", "Index"),
      GIN_DEFAULT(NODEF, "0"),
      GIN_DESC("The path to the split catalog files",
               "Flag to continue in case nothing is found [0 throws, 1 continues]")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesReadSpeciesSplitCatalog"),
      DESCRIPTION(R"--(See *abs_linesReadSpeciesSplitCatalog* but only for *abs_species*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("basename", "robust"),
      GIN_TYPE("String", "Index"),
      GIN_DEFAULT(NODEF, "0"),
      GIN_DESC("The path to the split catalog files",
               "Flag to continue in case nothing is found [0 throws, 1 continues]")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesSetEmpty"),
      DESCRIPTION(R"--(Empties *abs_lines_per_species* at the correct size.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesEmptyBroadeningParameters"),
               DESCRIPTION(R"--(Sets a broadening parameter to empty if it is effectively empty
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesNormalization"),
               DESCRIPTION(R"--(Sets normalization type for all lines

Available ``option``:

- ``"VVH"``: Van Vleck and Huber
- ``"VVW"``: Van Vleck and Weisskopf
- ``"RQ"``: Rosenkranz quadratic
- ``"SFS"``: Simple frequency scaling
- ``"None"``: No extra normalization

See the theory guide for more details.
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("option"),
               GIN_TYPE("String"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Method of line normalizations")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesNormalization"),
               DESCRIPTION(R"--(As *abs_linesNormalization* but for *abs_lines_per_species*
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN("option"),
               GIN_TYPE("String"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Method of line normalizations")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesNormalizationMatch"),
               DESCRIPTION(R"--(As *abs_linesNormalization* but for matching bands
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("option", "ID"),
               GIN_TYPE("String", "QuantumIdentifier"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Method of line normalizations",
                        "ID of one or more bands")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesNormalizationMatch"),
               DESCRIPTION(R"--(As *abs_lines_per_speciesNormalization* but for matching bands
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN("option", "ID"),
               GIN_TYPE("String", "QuantumIdentifier"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Method of line normalizations",
                        "ID of one or more bands")));

  md_data_raw.push_back(
    create_mdrecord(NAME("abs_lines_per_speciesNormalizationSpecies"),
             DESCRIPTION(R"--(As *abs_lines_per_speciesNormalization* but for matching *abs_species*
)--"),
             AUTHORS("Richard Larsson"),
             OUT("abs_lines_per_species"),
             GOUT(),
             GOUT_TYPE(),
             GOUT_DESC(),
             IN("abs_lines_per_species", "abs_species"),
             GIN("option", "species_tag"),
             GIN_TYPE("String", "String"),
             GIN_DEFAULT(NODEF, NODEF),
             GIN_DESC("Method of line normalizations",
                      "The species tag from *abs_species* to change")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesMirroring"),
               DESCRIPTION(R"--(Sets mirroring type for all lines.

Available ``option``:

- ``"None"``: No mirrored line
- ``"SameAsLineShape"``: Mirrored line broadened by line shape
- ``"Manual"``: Manually mirrored line (be careful; allows all frequencies)
- ``"Lorentz"``: Mirrored line broadened by Lorentz

Note that mirroring is never applied for DP line shape

Also note that Lorentz profile is approached by most line shapes at high frequency offset.

Also note that Manual settings are potentially dangerous as other frequency
offsets might not work as hoped.
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("option"),
               GIN_TYPE("String"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Method of line mirroring")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesMirroring"),
               DESCRIPTION(R"--(As *abs_linesMirroring* but for *abs_lines_per_species*
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN("option"),
               GIN_TYPE("String"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Method of line mirroring")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesMirroringMatch"),
               DESCRIPTION(R"--(As *abs_linesMirroring* but for matching bands
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("option", "ID"),
               GIN_TYPE("String", "QuantumIdentifier"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Method of line mirroring",
                        "ID of one or more bands")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesMirroringMatch"),
               DESCRIPTION(R"--(As *abs_lines_per_speciesMirroring* but for matching bands
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN("option", "ID"),
               GIN_TYPE("String", "QuantumIdentifier"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Method of line mirroring",
                        "ID of one or more bands")));

  md_data_raw.push_back(
    create_mdrecord(NAME("abs_lines_per_speciesMirroringSpecies"),
             DESCRIPTION(R"--(As *abs_lines_per_speciesMirroring* but for matching *abs_species*
)--"),
             AUTHORS("Richard Larsson"),
             OUT("abs_lines_per_species"),
             GOUT(),
             GOUT_TYPE(),
             GOUT_DESC(),
             IN("abs_lines_per_species", "abs_species"),
             GIN("option", "species_tag"),
             GIN_TYPE("String", "String"),
             GIN_DEFAULT(NODEF, NODEF),
             GIN_DESC("Method of line mirroring",
                      "The species tag from *abs_species* to change")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesSort"),
               DESCRIPTION(R"--(Sorts first the lines then the bands by smallest first
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("option"),
               GIN_TYPE("String"),
               GIN_DEFAULT("ByFrequency"),
               GIN_DESC("Sorting option")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesManualMirroring"),
               DESCRIPTION(R"--(Makes a copy of all lines at negative frequency setting them
to manual mirroring mode
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesManualMirroring"),
               DESCRIPTION(R"--(Makes a copy of all lines at negative frequency setting them.
to manual mirroring mode
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesManualMirroringSpecies"),
               DESCRIPTION(R"--(Calls *abs_linesManualMirroring* for given species in *abs_species*
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species", "abs_species"),
               GIN("species"),
               GIN_TYPE("ArrayOfSpeciesTag"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Species to mirror")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesPopulation"),
               DESCRIPTION(R"--(Sets population type for all lines.

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
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("option"),
               GIN_TYPE("String"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Method of line population")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesPopulation"),
               DESCRIPTION(R"--(As *abs_linesPopulation* but for *abs_lines_per_species*
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN("option"),
               GIN_TYPE("String"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Method of line population")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesPopulationMatch"),
               DESCRIPTION(R"--(As *abs_linesPopulation* but for matching bands
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("option", "ID"),
               GIN_TYPE("String", "QuantumIdentifier"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Method of line population",
                        "ID of one or more bands")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesPopulationMatch"),
               DESCRIPTION(R"--(As *abs_lines_per_speciesPopulation* but for matching bands
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN("option", "ID"),
               GIN_TYPE("String", "QuantumIdentifier"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Method of line population",
                        "ID of one or more bands")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesPopulationSpecies"),
               DESCRIPTION(R"--(As *abs_lines_per_speciesPopulation* but for matching *abs_species*
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species", "abs_species"),
               GIN("option", "species_tag"),
               GIN_TYPE("String", "String"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Method of line population",
                        "The species tag from *abs_species* to change")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesLineShapeType"),
               DESCRIPTION(R"--(Sets shape calculations type for all lines.

Available ``option``:

- ``"DP"``: Doppler profile
- ``"LP"``: Lorentz profile
- ``"VP"``: Voigt profile
- ``"SDVP"``: Speed-dependent Voigt profile
- ``"HTP"``: Hartman-Tran profile

See the theory guide for more details.
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("option"),
               GIN_TYPE("String"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Method of line shape calculations")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesLineShapeType"),
               DESCRIPTION(R"--(As *abs_linesLineShapeType* but for *abs_lines_per_species*
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN("option"),
               GIN_TYPE("String"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Method of line shape calculations")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesLineShapeTypeMatch"),
               DESCRIPTION(R"--(As *abs_linesLineShapeType* but for matching bands
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("option", "ID"),
               GIN_TYPE("String", "QuantumIdentifier"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Method of line shape calculations",
                        "ID of one or more bands")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesLineShapeTypeMatch"),
               DESCRIPTION(R"--(As *abs_lines_per_speciesLineShapeType* but for matching bands
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN("option", "ID"),
               GIN_TYPE("String", "QuantumIdentifier"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Method of line shape calculations",
                        "ID of one or more bands")));

  md_data_raw.push_back(
    create_mdrecord(NAME("abs_lines_per_speciesLineShapeTypeSpecies"),
             DESCRIPTION(R"--(As *abs_lines_per_speciesLineShapeType* but for matching *abs_species*
)--"),
             AUTHORS("Richard Larsson"),
             OUT("abs_lines_per_species"),
             GOUT(),
             GOUT_TYPE(),
             GOUT_DESC(),
             IN("abs_lines_per_species", "abs_species"),
             GIN("option", "species_tag"),
             GIN_TYPE("String", "String"),
             GIN_DEFAULT(NODEF, NODEF),
             GIN_DESC("Method of line shape calculations",
                      "The species tag from *abs_species* to change")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesCutoff"),
               DESCRIPTION(R"--(Sets cutoff type and magnitude for all lines.

The line is cut off when this is active at the given frequency.
The only non-zero range is from this range to its negative equivalent

Available ``option``:

- ``"None"``: No cutoff
- ``"ByLine"``: Cutoff relative to a speed-independent shifted line center, highest frequency: F0+cutoff+D0

For "ByLine", the negative frequency is at F0-cutoff-D0
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("option", "value"),
               GIN_TYPE("String", "Numeric"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Method of line shape calculations",
                        "Value of cutoff")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesCutoff"),
               DESCRIPTION(R"--(As *abs_linesCutoff* but for *abs_lines_per_species*
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN("option", "value"),
               GIN_TYPE("String", "Numeric"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Method of line shape calculations",
                        "Value of cutoff")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesCutoffMatch"),
               DESCRIPTION(R"--(As *abs_linesCutoff* but for matching bands
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("option", "value", "ID"),
               GIN_TYPE("String", "Numeric", "QuantumIdentifier"),
               GIN_DEFAULT(NODEF, NODEF, NODEF),
               GIN_DESC("Method of line shape calculations",
                        "Value of cutoff",
                        "ID of one or more bands")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesCutoffMatch"),
               DESCRIPTION(R"--(As *abs_lines_per_speciesCutoff* but for matching bands
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN("option", "value", "ID"),
               GIN_TYPE("String", "Numeric", "QuantumIdentifier"),
               GIN_DEFAULT(NODEF, NODEF, NODEF),
               GIN_DESC("Method of line shape calculations",
                        "Value of cutoff",
                        "ID of one or more bands")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesCutoffSpecies"),
               DESCRIPTION(R"--(As *abs_lines_per_speciesCutoff* but for matching *abs_species*
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species", "abs_species"),
               GIN("option", "value", "species_tag"),
               GIN_TYPE("String", "Numeric", "String"),
               GIN_DEFAULT(NODEF, NODEF, NODEF),
               GIN_DESC("Method of line shape calculations",
                        "Value of cutoff",
                        "The species tag from *abs_species* to change")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesLinemixingLimit"),
               DESCRIPTION(R"--(Sets line mixing limit for all lines.

If value is less than 0, no limit is applied and line mixing is active.
Otherwise, line mixing is inactive if the pressure is below the limit.
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("value"),
               GIN_TYPE("Numeric"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Value of limit")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesLinemixingLimit"),
               DESCRIPTION(R"--(See *abs_linesLinemixingLimit*
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN("value"),
               GIN_TYPE("Numeric"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Value of limit")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesLinemixingLimitMatch"),
               DESCRIPTION(R"--(See *abs_linesLinemixingLimit* for values

This function only acts on matches between the bands and input ID
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("value","ID"),
               GIN_TYPE("Numeric", "QuantumIdentifier"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Value of limit",
                        "ID of one or more bands")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesLinemixingLimitMatch"),
               DESCRIPTION(R"--(See *abs_linesLinemixingLimit* for values

This function only acts on matches between the bands and input ID
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN("value","ID"),
               GIN_TYPE("Numeric", "QuantumIdentifier"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Value of limit",
                        "ID of one or more bands")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesLinemixingLimitSpecies"),
               DESCRIPTION(R"--(See *abs_linesLinemixingLimit* but for single species
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species", "abs_species"),
               GIN("value", "species_tag"),
               GIN_TYPE("Numeric", "String"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Value of limit",
                        "The species tag from *abs_species* to change")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesT0"),
               DESCRIPTION(R"--(Sets reference temperature for all lines.
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("value"),
               GIN_TYPE("Numeric"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Value of T0")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesT0"),
               DESCRIPTION(R"--(See *abs_linesT0*
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN("value"),
               GIN_TYPE("Numeric"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Value of T0")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_linesT0Match"),
               DESCRIPTION(R"--(Sets reference temperature

This function only acts on matches between the bands and input ID
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("value", "ID"),
               GIN_TYPE("Numeric", "QuantumIdentifier"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Value of T0",
                        "ID of one or more bands")));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_lines_per_speciesT0Match"),
               DESCRIPTION(R"--(Sets reference temperature

This function only acts on matches between the bands and input ID
)--"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines_per_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_species"),
               GIN("value", "ID"),
               GIN_TYPE("Numeric", "QuantumIdentifier"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Value of T0",
                        "ID of one or more bands")));

  md_data_raw.push_back(
    create_mdrecord(NAME("abs_lines_per_speciesT0Species"),
             DESCRIPTION(R"--(See *abs_linesT0* but for single species
)--"),
             AUTHORS("Richard Larsson"),
             OUT("abs_lines_per_species"),
             GOUT(),
             GOUT_TYPE(),
             GOUT_DESC(),
             IN("abs_lines_per_species", "abs_species"),
             GIN("value", "species_tag"),
             GIN_TYPE("Numeric", "String"),
             GIN_DEFAULT(NODEF, NODEF),
             GIN_DESC("Value of T0",
                      "The species tag from *abs_species* to change")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_linesChangeBaseParameterForMatchingLevel"),
      DESCRIPTION(R"--(Change parameter of all levels in *abs_lines* that match with *QuantumIdentifier*.

Only works for these ``parameter_name``:
 - ``\"Statistical Weight\"``
 - ``\"Zeeman Coefficient\"``
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("QI", "parameter_name", "change", "relative"),
      GIN_TYPE("QuantumIdentifier", "String", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "0"),
      GIN_DESC("Information to match the level.",
               "Name of parameter to be replaced",
               "Value with which to change matching level's value",
               "Flag for relative change (0 is absolute change)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_linesChangeBaseParameterForMatchingLevels"),
      DESCRIPTION(R"--(See *abs_linesChangeBaseParameterForMatchingLevel*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("QI", "parameter_name", "change", "relative"),
      GIN_TYPE("ArrayOfQuantumIdentifier", "String", "Vector", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "0"),
      GIN_DESC("Information to match the level.",
               "Name of parameter to be replaced",
               "Value with which to change matching level's value",
               "Flag for relative change (0 is absolute change)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesChangeBaseParameterForMatchingLevel"),
      DESCRIPTION(R"--(See *abs_linesChangeBaseParameterForMatchingLevel*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN("QI", "parameter_name", "change", "relative"),
      GIN_TYPE("QuantumIdentifier", "String", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "0"),
      GIN_DESC("Information to match the level.",
               "Name of parameter to be replaced",
               "Value with which to change matching level's value",
               "Flag for relative change (0 is absolute change)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesChangeBaseParameterForMatchingLevels"),
      DESCRIPTION(R"--(See *abs_linesChangeBaseParameterForMatchingLevel*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN("QI", "parameter_name", "change", "relative"),
      GIN_TYPE("ArrayOfQuantumIdentifier", "String", "Vector", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "0"),
      GIN_DESC("Information to match the level.",
               "Name of parameter to be replaced",
               "Value with which to change matching level's value",
               "Flag for relative change (0 is absolute change)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_linesBaseParameterMatchingLevel"),
      DESCRIPTION(R"--(Set parameter of all levels in *abs_lines* that match with *QuantumIdentifier*.

Only works for these ``parameter_name``:
 - ``\"Statistical Weight\"``
 - ``\"Zeeman Coefficient\"``
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("QI", "parameter_name", "change"),
      GIN_TYPE("QuantumIdentifier", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Information to match the level.",
               "Name of parameter to be replaced",
               "Value with which to set matching level's value")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_linesBaseParameterMatchingLevels"),
      DESCRIPTION(R"--(See *abs_linesBaseParameterMatchingLevel*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("QI", "parameter_name", "change"),
      GIN_TYPE("ArrayOfQuantumIdentifier", "String", "Vector"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Information to match the level.",
               "Name of parameter to be replaced",
               "Value with which to set matching level's value")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesBaseParameterMatchingLevel"),
      DESCRIPTION(R"--(See *abs_linesBaseParameterMatchingLevel*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN("QI", "parameter_name", "change"),
      GIN_TYPE("QuantumIdentifier", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Information to match the level.",
               "Name of parameter to be replaced",
               "Value with which to set matching level's value")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesBaseParameterMatchingLevels"),
      DESCRIPTION(R"--(See *abs_linesBaseParameterMatchingLevel*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN("QI", "parameter_name", "change"),
      GIN_TYPE("ArrayOfQuantumIdentifier", "String", "Vector"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Information to match the level.",
               "Name of parameter to be replaced",
               "Value with which to set matching level's value")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_linesChangeBaseParameterForMatchingLines"),
      DESCRIPTION(R"--(Change parameter of all lines in *abs_lines* that match with *QuantumIdentifier*.

Only works for these ``parameter_name``:
 - ``\"Central Frequency\"``
 - ``\"Line Strength\"``
 - ``\"Lower State Energy\"``
 - ``\"Einstein Coefficient\"``
 - ``\"Lower Statistical Weight\"``
 - ``\"Upper Statistical Weight\"``
 - ``\"Lower Zeeman Coefficient\"``
 - ``\"Upper Zeeman Coefficient\"``
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("QI", "parameter_name", "change", "relative"),
      GIN_TYPE("QuantumIdentifier", "String", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "0"),
      GIN_DESC("Information to match the line/band.",
               "Name of parameter to be replaced",
               "Value with which to change matching line's value",
               "Flag for relative change (0 is absolute change)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesChangeBaseParameterForMatchingLines"),
      DESCRIPTION(R"--(See *abs_linesChangeBaseParameterForMatchingLines*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN("QI", "parameter_name", "change", "relative"),
      GIN_TYPE("QuantumIdentifier", "String", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "0"),
      GIN_DESC("Information to match the line/band.",
               "Name of parameter to be replaced",
               "Value with which to change matching line's value",
               "Flag for relative change (0 is absolute change)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesChangeBaseParameterForSpecies"),
      DESCRIPTION(R"--(See *abs_linesChangeBaseParameterForMatchingLines* but for single species
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species", "abs_species"),
      GIN("QI", "parameter_name", "change", "relative", "species_tag"),
      GIN_TYPE("QuantumIdentifier", "String", "Numeric", "Index", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "0", NODEF),
      GIN_DESC("Information to match the line/band.",
               "Name of parameter to be replaced",
               "Value with which to change matching line's value",
               "Flag for relative change (0 is absolute change)",
               "The species tag from *abs_species* to change")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_linesBaseParameterMatchingLines"),
      DESCRIPTION(R"--(Set parameter of all lines in *abs_lines* that match with *QuantumIdentifier*.

Only works for these ``parameter_name``:
 - ``\"Central Frequency\"``
 - ``\"Line Strength\"``
 - ``\"Lower State Energy\"``
 - ``\"Einstein Coefficient\"``
 - ``\"Lower Statistical Weight\"``
 - ``\"Upper Statistical Weight\"``
 - ``\"Lower Zeeman Coefficient\"``
 - ``\"Upper Zeeman Coefficient\"``
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("QI", "parameter_name", "change"),
      GIN_TYPE("QuantumIdentifier", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Information to match the line/band.",
               "Name of parameter to be replaced",
               "Value with which to change matching line's value")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_linesLineShapeModelParametersMatchingLines"),
      DESCRIPTION(R"--(Sets line shape model data parameter in matching lines.

The matching is done so that QI must be in the line identifier

Acceptable ``parameter`` (s) are:
 - ``\"G0\"``
 - ``\"D0\"``
 - ``\"G2\"``
 - ``\"D2\"``
 - ``\"FVC\"``
 - ``\"ETA\"``
 - ``\"Y\"``
 - ``\"G\"``
 - ``\"DV\"``

Acceptable ``temperaturemodel`` (s) are:
 - ``\"None\"``
 - ``\"T0\"``
 - ``\"T1\"``
 - ``\"T2\"``
 - ``\"T3\"``
 - ``\"T4\"``
 - ``\"T5\"``
 - ``\"LM_AER\"``
 - ``\"DPL\"``

Acceptable ``species`` are:
 - ``\"AIR\"`` (so long as it is the broadening species list)
 - ``\"SELF\"`` (so long as it is the broadening species list)
 - Any species in the line broadening species

See the user guide for the meanings of all of these keywords
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("QI", "parameter", "species", "temperaturemodel", "new_values"),
      GIN_TYPE("QuantumIdentifier",
               "String",
               "String",
               "String",
               "Vector"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, NODEF),
      GIN_DESC("Information to match the line.",
               "Name of parameter to be replaced",
               "Species of parameter to be changed",
               "Temperature model for the new values",
               "Sets the values found")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesLineShapeModelParametersMatchingLines"),
      DESCRIPTION(R"--(See *abs_linesLineShapeModelParametersMatchingLines*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN("QI", "parameter", "species", "temperaturemodel", "new_values"),
      GIN_TYPE("QuantumIdentifier",
               "String",
               "String",
               "String",
               "Vector"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, NODEF),
      GIN_DESC("Information to match the line.",
               "Name of parameter to be replaced",
               "Species of parameter to be changed",
               "Temperature model for the new values",
               "Sets the values found")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_linesZeemanCoefficients"),
      DESCRIPTION(R"--(Sets the Zeeman coefficients of the lines by user input

The matching is permissive, all in qid must just match.  If there
are multiple matches, the last match rules
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("qid", "gs"),
      GIN_TYPE("ArrayOfQuantumIdentifier", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Information to match an energy level of a/many lines.",
               "Corresponding value to set as Zeeman coefficient")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesZeemanCoefficients"),
      DESCRIPTION(R"--(See *abs_linesZeemanCoefficients*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN("qid", "gs"),
      GIN_TYPE("ArrayOfQuantumIdentifier", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Information to match an energy level of a/many lines.",
               "Corresponding value to set as Zeeman coefficient")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_linesCompact"),
      DESCRIPTION(R"--(Removes lines that are unimportant because of their
cutoff frequency range
)--"),
      AUTHORS("Stefan Buehler", "Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines", "f_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesCompact"),
      DESCRIPTION(R"--(See *abs_linesCompact*
)--"),
      AUTHORS("Stefan Buehler", "Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species", "f_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesCreateFromLines"),
      DESCRIPTION(R"--(Split lines up into the different species.

The order of the splitting will match the outer layer of *abs_species*
There will be no respect for the internal layer of *abs_species*
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines", "abs_species"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_linesTurnOffLineMixing"),
      DESCRIPTION(R"--(Sets all line mixing parameters to emtpy.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesTurnOffLineMixing"),
      DESCRIPTION(R"--(Sets all line mixing parameters to emtpy.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lookupAdapt"),
      DESCRIPTION(R"--(Adapts a gas absorption lookup table to the current calculation.

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
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_lookup", "abs_lookup_is_adapted"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lookup", "abs_species", "f_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lookupCalc"),
      DESCRIPTION(R"--(Creates a gas absorption lookup table.

The lookup table stores absorption cross-sections as a function of
pressure. Additionally, absorption can be stored as a function of
temperature for temperature perturbations from a reference
profile.

Additionally, absorption can be stored as a function of water vapor
VMR perturbations from a reference profile. The variable *abs_nls*
specifies, for which species water vapor perturbations should be
generated.

Note, that the absorbing gas can be any gas, but the perturbing gas is
always H2O.
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_lookup", "abs_lookup_is_adapted"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species",
         "abs_nls",
         "f_grid",
         "abs_p",
         "abs_vmrs",
         "abs_t",
         "abs_t_pert",
         "abs_nls_pert",
         "propmat_clearsky_agenda"),
      GIN("lowest_vmr"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("1e-9"),
      GIN_DESC("Lowest possible VMR to compute absorption at")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lookupInit"),
      DESCRIPTION(R"--(Creates an empty gas absorption lookup table.

This is mainly there to help developers. For example, you can write
the empty table to an XML file, to see the file format.
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_lookup"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lookupSetup"),
      DESCRIPTION(R"--(Set up input parameters for abs_lookupCalc.

More information can be found in the documentation for method
*abs_lookupSetupBatch*

Max and min values of H2O and temperature are adjusted to allow for
numerical perturbations in Jacobian calculation.

The input variables *abs_nls_interp_order* and *abs_t_interp_order*
are used to make sure that there are enough points in *abs_nls_pert*
and *abs_t_pert* for the chosen interpolation order.

Note: For homogeneous 1D cases, it can be advantageous to calculate
*abs_lookup* from the 1D atmosphere, and to expand the atmosphere
to 3D only after that. This particularly if nonlinear species
(i.e., H2O) are involved.

See also: *abs_lookupSetupBatch*
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_p",
          "abs_t",
          "abs_t_pert",
          "abs_vmrs",
          "abs_nls",
          "abs_nls_pert"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         "atm_field",
         "atmfields_checked",
         "abs_species",
         "abs_p_interp_order",
         "abs_t_interp_order",
         "abs_nls_interp_order"),
      GIN("p_step", "t_step", "h2o_step"),
      GIN_TYPE("Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT("0.05", "100", "100"),
      GIN_DESC(/* p_step */
               "Maximum step in log10(p[Pa]). If the pressure grid is "
               "coarser than this, additional points are added until each "
               "log step is smaller than this.",
               /* t_step */
               "The temperature variation grid step in Kelvin, "
               "for a 2D or 3D atmosphere. For a 1D atmosphere this "
               "parameter is not used.",
               /* h2o_step */
               "The H2O variation grid step [fractional], if "
               "H2O variations are done (which is determined automatically, "
               "based on abs_species and the atmospheric dimension). For a "
               "1D atmosphere this parameter is not used.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lookupSetupBatch"),
      DESCRIPTION(R"--(Set up input parameters for abs_lookupCalc for batch calculations.

This method performs a similar task as *abs_lookupSetup*, with the
difference that the lookup table setup is not for a single
atmospheric state, but for a whole batch of them, stored in
*batch_atm_fields_compact*.

The method checks *abs_species* to decide which species require
nonlinear treatment in the lookup table.

The method also checks which range of pressures, temperatures, and
VMRs occurs, and sets *abs_p*, *abs_t*, *abs_t_pert*, and *abs_vmrs*
accordingly.

If nonlinear species are present, *abs_nls* and *abs_nls_pert* are also
generated.

Max and min values of H2O and temperature are adjusted to allow for
numerical perturbations in Jacobian calculation.

The input variables *abs_nls_interp_order* and *abs_t_interp_order*
are used to make sure that there are enough points in *abs_nls_pert*
and *abs_t_pert* for the chosen interpolation order.

The method checks each given field using *atmfields_checkedCalc*.
If a field does not pass the check, a run-time error is thrown.
To prevent this, the parameter ``robust`` can be set to one: Invalid 
atmospheres are skipped, but the run continues. This matches the 
robust behaviour of *ybatchCalc*.

See also:
   *abs_lookupSetup*
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_p",
          "abs_t",
          "abs_t_pert",
          "abs_vmrs",
          "abs_nls",
          "abs_nls_pert"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species",
         "batch_atm_fields_compact",
         "abs_p_interp_order",
         "abs_t_interp_order",
         "abs_nls_interp_order"),
      GIN("p_step",
          "t_step",
          "h2o_step",
          "extremes",
          "robust",
          "check_gridnames"),
      GIN_TYPE("Numeric", "Numeric", "Numeric", "Vector", "Index", "Index"),
      GIN_DEFAULT("0.05", "20", "100", "[]", "0", "0"),
      GIN_DESC(/* p_step */
               "Grid step in log10(p[Pa]) (base 10 logarithm).",
               /* t_step */
               "The temperature variation grid step in Kelvin. The true "
               "step can become finer than this, if required by the "
               "interpolation order.",
               /* h2o_step */
               "The H2O variation grid step [fractional], if H2O variations "
               "are done (which is determined automatically, based on "
               "abs_species and the atmospheric dimension). As for T, the true "
               "step can turn out finer if required by the interpolation order.",
               /* extremes */
               "You can give here explicit extreme values to add to "
               "abs_t_pert and abs_nls_pert. The order is [t_pert_min, "
               "t_pert_max, nls_pert_min, nls_pert_max].",
               /* robust */
               "A flag with value 1 or 0. If set to one, the batch"
               " setup will continue, even if individual fields are invalid."
               " This is consistent with the behaviour of *ybatchCalc*.",
               /* check_gridnames */
               "A flag with value 1 or 0. If set to one, the gridnames of"
               " every *atm_fields_compact* are checked.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lookupSetupWide"),
      DESCRIPTION(R"--(Set up input parameters for abs_lookupCalc for a wide range of
atmospheric conditions.

This method can be used to set up parameters for a lookup table that
really covers all reasonable atmospheric conditions.

Reference profiles of T and H2O will be constant, so that the
different dimensions in the lookup table are actually \"orthogonal\",
unlike the traditional case where we have pressure dependent reference
profiles. This makes the table numerically somewhat more robust then
the traditional ones, and it makes it straightforward to calculate the
accuracy for the different interpolations with abs_lookupTestAccuracy.

You can give min an max values for the atmospheric conditions. The
default values are chosen such that they cover the value range over
the complete Chevallier91L data set, and a bit more. The statistics
of the Chevallier91L data are::

  min(p)   / max(p)   [Pa]:  1 / 104960
  min(T)   / max(T)   [K]:   158.21 / 320.39
  min(H2O) / max(H2O) [VMR]: -5.52e-07 / 0.049
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_p",
          "abs_t",
          "abs_t_pert",
          "abs_vmrs",
          "abs_nls",
          "abs_nls_pert"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species",
         "abs_p_interp_order",
         "abs_t_interp_order",
         "abs_nls_interp_order"),
      GIN("p_min", "p_max", "p_step", "t_min", "t_max", "h2o_min", "h2o_max"),
      GIN_TYPE("Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric"),
      GIN_DEFAULT("0.5", "110000", "0.05", "100", "400", "0", "0.05"),
      GIN_DESC("Pressure grid minimum [Pa].",
               "Pressure grid maximum [Pa].",
               "Pressure grid step in log10(p[Pa]) (base 10 logarithm).",
               "Temperature grid minimum [K].",
               "Temperature grid maximum [K].",
               "Humidity grid minimum [fractional].",
               "Humidity grid maximum [fractional].")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_speciesAdd"),
      DESCRIPTION(R"--(Adds species tag groups to the list of absorption species.

This WSM is similar to *abs_speciesSet*, the only difference is that
this method appends species to an existing list of absorption species instead
of creating the whole list.

See *abs_speciesSet* for details on how tags are defined and examples of
how to input them in the control file.
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_species",
          "propmat_clearsky_agenda_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("species"),
      GIN_TYPE("ArrayOfString"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Specify one String for each tag group that you want to"
               " add. Inside the String, separate the tags by commas"
               " (plus optional blanks).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_speciesAdd2"),
      DESCRIPTION(R"--(Adds a species tag group to the list of absorption species and
jacobian quantities.

The method is basically a combined call of *abs_speciesAdd* and
*jacobianAddAbsSpecies*. In this way it is not needed to specify a
tag group in two different places.

Arguments exactly as for *jacobianAddAbsSpecies*. Note that this
method only handles a single tag group, in contrast to
*abs_speciesAdd*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("abs_species",
          "jacobian_quantities",
          "jacobian_agenda",
          "propmat_clearsky_agenda_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("gin1", "gin2", "gin3", "species", "unit"),
      GIN_TYPE("Vector", "Vector", "Vector", "String", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, "vmr"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "The species tag of the retrieval quantity.",
               "Retrieval unit. See above."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_speciesDefineAllInScenario"),
      DESCRIPTION(R"--(Define one tag group for each species known to ARTS and included in an
atmospheric scenario.

You can use this as an alternative to *abs_speciesSet* if you want to make an
absorption calculation that is as complete as possible. The method
goes through all defined species and tries to open the VMR file. If
this works the tag is included, otherwise it is skipped.
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_species",
          "propmat_clearsky_agenda_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("basename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC(
          "The name and path of a particular atmospheric scenario."
          " For example: /pool/lookup2/arts-data/atmosphere/fascod/tropical")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_speciesDefineAll"),
      DESCRIPTION(R"--(Sets *abs_species* [i][0] to all species in ARTS
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_species",
          "propmat_clearsky_agenda_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord(NAME("abs_speciesInit"),
               DESCRIPTION(R"--(Sets  *abs_species* to be empty.
)--"),
               AUTHORS("Stefan Buehler"),
               OUT("abs_species"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_speciesSet"),
      DESCRIPTION(R"--(Set up a list of absorption species tag groups.

Workspace variables like *abs_species* contain several tag
groups. Each tag group contains one or more tags. This method converts
descriptions of tag groups given in the keyword to the ARTS internal
representation (an *ArrayOfArrayOfSpeciesTag*). A tag group selects
spectral features which belong to the same species.

A tag is defined in terms of the name of the species, isotopologue, and a
range of frequencies. Species are named after the standard chemical
names, e.g., ``\"O3\"``. Isotopologues are given by the last digit of the atomic
weight, i.g., ``\"O3-668\"`` for the asymmetric ozone molecule including an
oxygen 18 atom. Groups of transitions are specified by giving a lower
and upper limit of a frequency range, e.g., ``\"O3-666-500e9-501e9\"``.

To turn on Zeeman calculation for a species, ``\"-Z\"`` may be appended
to its name: ``\"O2-Z\"`` or ``\"O2-Z-66\"``

The symbol ``\"*\"`` acts as a wild card. Furthermore, frequency range or
frequency range and isotopologue may be omitted.

Finally, instead of the isotopologue the special letter ``\"nl\"`` may be given,
e.g., ``\"H2O-nl\"``. This means that no absorption at all is associated
with this tag. (It is not quite clear if this feature is useful for
anything right now.)

Example:

>>> species = [ \"O3-666-500e9-501e9, O3-686\", \"O3\", \"H2O-PWR98\" ]

   The first tag group selects all O3-666 lines between 500 and
   501 GHz plus all O3-686 lines. 

   The second tag group selects all remaining O3 transitions.

   The third tag group selects H2O, with one of the complete
   absorption models (Rosenkranz 98). No spectrocopic line catalogue
   data will be used for that third tag group.  For more available full
   absorption models see *propmat_clearskyAddPredefined*

   Note that order of tag groups in the species list matters. In our
   example, changing the order of the first two tag group will give
   different results: as ``\"O3\"`` already selects all O3 transitions,
   no lines will remain to be selected by the
   ``\"O3-666-500e9-501e9, O3-686\"`` tag.

For CIA species the tag consists of the two involved species and
a dataset index. CIA species can be defined for multiple regions
The dataset index determines which region to use from the corresponding
CIARecord in *abs_cia_data*.

Example

>>> species = [ \"N2-CIA-N2-0, N2-CIA-N2-1\" ]

For Hitran cross section species the tag consists of the species and
the tagtype XFIT, e.g. CFC11-XFIT. The data for the species must be
available in the *xsec_fit_data* variable.

*propmat_clearsky_agenda_checked* is set to be false.
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_species",
          "propmat_clearsky_agenda_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("species"),
      GIN_TYPE("ArrayOfString"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Specify one String for each tag group that you want to"
               " create. Inside the String, separate the tags by commas"
               " (plus optional blanks).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_vecAddGas"),
      DESCRIPTION(R"--(Add gas absorption to first element of absorption vector.

The task of this method is to sum up the gas absorption of the
different gas species and add the result to the first element of the
absorption vector.
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_vec"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_vec", "propmat_clearsky"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearsky_agendaGUI"),
      DESCRIPTION(R"--(Opens a GUI for running the propagation matrix agenda

Note that this is not thread-safe and should be executed on the main workspace

The values of all non-control flow are automatically loaded from the workspace
if they are defined.  Otherwise some values are just selected
)--"
      ),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky_agenda", "abs_species"),
      GIN("load"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Load non-logical variables from workspace if true")));

  md_data_raw.push_back(create_mdrecord(
      NAME("predefined_model_dataInit"),
      DESCRIPTION(R"--(Initialize the predefined model data
)--"),
      AUTHORS("Richard Larsson"),
      OUT("predefined_model_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("predefined_model_dataAddWaterMTCKD400"),
      DESCRIPTION(R"--(Sets the data for MT CKD 4.0 Water model

Note that the vectors must have the same length, and that wavenumbers must be growing
at a constant rate.  The minimum length is 4.

Note also that as this is predefined model data, the units of the values of the vectors
must be as described by each vector.
)--"),
      AUTHORS("Richard Larsson"), OUT("predefined_model_data"), GOUT(),
      GOUT_TYPE(), GOUT_DESC(), IN("predefined_model_data"),
      GIN("ref_temp", "ref_press", "ref_h2o_vmr", "self_absco_ref",
          "for_absco_ref", "wavenumbers", "self_texp"),
      GIN_TYPE("Numeric", "Numeric", "Numeric", "Vector", "Vector", "Vector",
               "Vector"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, NODEF, NODEF, NODEF),
      GIN_DESC("Reference temperature", "Reference pressure",
               "Reference volume mixing ratio of water",
               "Self absorption [1/(cm-1 molecules/cm^2]",
               "Foreign absorption [1/(cm-1 molecules/cm^2)]",
               "Wavenumbers [cm-1]", "Self temperature exponent [-]")));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyAddPredefined"),
      DESCRIPTION(R"--(Adds all of the predefined models in *abs_species* to the propmat_clearsky

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
)--"),
      AUTHORS("Richard Larsson"),
      OUT("propmat_clearsky",
          "dpropmat_clearsky_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
         "dpropmat_clearsky_dx",
         "predefined_model_data",
         "abs_species",
         "select_abs_species",
         "jacobian_quantities",
         "f_grid",
         "atm_point"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("AgendaAppend"),
      DESCRIPTION(R"--(Append methods to an agenda.

An agenda is used to store a list of methods that are meant to be
executed sequentially.

This method takes the methods given in the body (in the curly braces)
and appends them to the agenda given by the output argument (in the round
braces).

It also uses the agenda lookup data (defined in file agendas.cc) to
check, whether the given methods use the right input WSVs and produce
the right output WSVs.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Agenda"),
      GOUT_DESC("Target agenda."),
      IN(),
      GIN("input"),
      GIN_TYPE("Agenda"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Source agenda."),
      SETMETHOD(false),
      AGENDAMETHOD(true),
      USES_TEMPLATES(false),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(create_mdrecord(NAME("AgendaExecute"),
                                 DESCRIPTION(R"--(Execute an agenda.
)--"),
                                 AUTHORS("Oliver Lemke"),
                                 OUT(),
                                 GOUT(),
                                 GOUT_TYPE(),
                                 GOUT_DESC(),
                                 IN(),
                                 GIN("a"),
                                 GIN_TYPE("Agenda"),
                                 GIN_DEFAULT(NODEF),
                                 GIN_DESC("Agenda to be executed."),
                                 SETMETHOD(false),
                                 AGENDAMETHOD(false)));

  md_data_raw.push_back(create_mdrecord(
      NAME("AgendaExecuteExclusive"),
      DESCRIPTION(R"--(Execute an agenda exclusively.

Only one call to *AgendaExecuteExclusive* is executed at a time.
Other calls to this function are blocked until the current one
finishes. WARNING: Can cause deadlocks! Use with care.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("a"),
      GIN_TYPE("Agenda"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Agenda to be executed."),
      SETMETHOD(false),
      AGENDAMETHOD(false)));

  md_data_raw.push_back(create_mdrecord(
      NAME("AgendaSet"),
      DESCRIPTION(R"--(Set up an agenda.

An agenda is used to store a list of methods that are meant to be
executed sequentially.

This method takes the methods given in the body (in the curly braces)
and puts them in the agenda given by the output argument (in the round
braces).

It also uses the agenda lookup data (defined in file agendas.cc) to
check, whether the given methods use the right input WSVs and
produce the right output WSVs.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Agenda"),
      GOUT_DESC("The new agenda."),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC(),
      SETMETHOD(false),
      AGENDAMETHOD(true),
      USES_TEMPLATES(false),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("AltLatLonFieldSet"),
      DESCRIPTION(R"--(Fills an altitude-latitude-longitude field with given input.

Grids and data must match in size.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("gfield3"),
      GOUT_TYPE("GriddedField3"),
      GOUT_DESC("Field to set."),
      IN(),
      GIN("altitude_grid", "latitude_grid", "longitude_grid", "data", "name"),
      GIN_TYPE("Vector", "Vector", "Vector", "Tensor3", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, ""),
      GIN_DESC("The altitude grid of ``data``.",
               "The latitude grid of ``data``.",
               "The longitude grid of ``data``.",
               "The data of the field (will become gfield2.data).",
               "The name of the field (will become gfield2.name).")));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("AltLatLonFieldSetConstant"),
      DESCRIPTION(R"--(Sets an altitude-latitude-longitude field to have a constant data value.

All three grids grids are set to have length one, with the value 0.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("gfield3"),
      GOUT_TYPE("GriddedField3"),
      GOUT_DESC("Field to set."),
      IN(),
      GIN("value", "name"),
      GIN_TYPE("Numeric", "String"),
      GIN_DEFAULT(NODEF, ""),
      GIN_DESC("The value (to place in gfield3.data).",
               "The name of the field (will become gfield3.name).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("AngularGridsSetFluxCalc"),
      DESCRIPTION(R"--(Sets the angular grids for the calculation of radiation fluxes. 

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
)--"),
      AUTHORS("Manfred Brath"),
      OUT("za_grid", "aa_grid", "za_grid_weights"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("N_za_grid", "N_aa_grid", "za_grid_type"),
      GIN_TYPE("Index", "Index", "String"),
      GIN_DEFAULT("2", "1", "linear_mu"),
      GIN_DESC("Number of zenith angles",
               "Number of azimuth angles",
               "Zenith angle grid type")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ArrayOfAgendaAppend"),
      DESCRIPTION(R"--(Set up an agenda and append it to the array of agendas.

See *AgendaSet* for details.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("ArrayOfAgenda"),
      GOUT_DESC("The new agenda."),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC(),
      SETMETHOD(false),
      AGENDAMETHOD(true),
      USES_TEMPLATES(false),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("ArrayOfAgendaExecute"),
               DESCRIPTION(R"--(Execute an agenda from an ArrayOfAgenda.
)--"),
               AUTHORS("Oliver Lemke"),
               OUT(),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("agenda_array_index"),
               GIN("agendas"),
               GIN_TYPE("ArrayOfAgenda"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Array of agendas."),
               SETMETHOD(false),
               AGENDAMETHOD(false)));

  md_data_raw.push_back(create_mdrecord(
      NAME("AntennaMultiBeamsToPencilBeams"),
      DESCRIPTION(R"--(Maps a multi-beam case to a matching pencil beam case.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_pos",
          "sensor_los",
          "antenna_dlos",
          "antenna_dim",
          "mblock_dlos"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_pos",
         "sensor_los",
         "antenna_dlos",
         "antenna_dim",
         "mblock_dlos"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("AntennaOff"),
      DESCRIPTION(R"--(Sets some antenna related variables

Use this method to set *antenna_dim* and *mblock_dlos* to
suitable values (1 and [0], respectively) for cases when a
sensor is included, but the antenna pattern is neglected.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("antenna_dim", "mblock_dlos"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("antenna_responseGaussian"),
      DESCRIPTION(R"--(Sets up a Gaussian antenna response.

This method works as *antenna_responseGaussianConstant* but allows
to inlude a frequency variation of the antenna width. Here the FWHM
is specified at a set of frequencies. These frequencies will also be
the frequency grid of *antenna_response*.

If ``grid_width`` is set to <=0, the grid width will be twice the max
value in ``fwhm``.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("antenna_response"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("f_points", "fwhm", "grid_width", "grid_npoints", "do_2d"),
      GIN_TYPE("Vector", "Vector", "Numeric", "Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF, "-1.0", "21", "0"),
      GIN_DESC("Frequencies at which FWHM is defined.",
               "Full width at half-maximum of the Gaussian function.",
               "Full width of grid (negative value gives 2*fwhm).",
               "Number of points to represent the grid, see above.",
               "Set to 1 to create a 2D antenna pattern.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("antenna_responseGaussianConstant"),
      DESCRIPTION(R"--(Sets up a Gaussian antenna response, with no frequency variation.

The method assumes that the response is the same for all
frequencies and polarisations, and that it can be modelled
as Gaussian. The width of the Gaussian is specified by its
full width at half maximum (FWHM).

The grid generated has ``grid_npoints`` equidistant values, with
the first one at -grid_width/2 and the last one at grid_width/2.

If ``grid_width`` is set to <= 0, a default of twice the FWMH is
applied. This gives a coverage of about 98\% of the response.

The default for ``grid_npoints`` is 21. When the grid width is 2*FWHM,
that default value gives an error < 0.001 of the integrated response
using trapezoidal integration. ``grid_npoints`` must be > 1.

If the 2D option is selected (``do_2d``), a circular antenna is
assumed. The same grid and FWHM is applied in both dimensions.

If the grid has a sufficiently high width the integral of the
response is 1. Otherwise the integral is smaller than 1. That
is, no normalisation is applied.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("antenna_response"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("fwhm", "grid_width", "grid_npoints", "do_2d"),
      GIN_TYPE("Numeric", "Numeric", "Index", "Index"),
      GIN_DEFAULT(NODEF, "-1.0", "21", "0"),
      GIN_DESC("Full width at half-maximum of the Gaussian function.",
               "Full width of grid (negative value gives 2*fwhm).",
               "Number of points to represent the grid, see above.",
               "Set to 1 to create a 2D antenna pattern.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("antenna_responseGaussianEffectiveSize"),
      DESCRIPTION(R"--(Sets up Gaussian antenna responses.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("antenna_response"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("leff", "grid_width", "grid_npoints", "nf", "fstart", "fstop", "do_2d"),
      GIN_TYPE("Numeric",
               "Numeric",
               "Index",
               "Index",
               "Numeric",
               "Numeric",
               "Index"),
      GIN_DEFAULT(NODEF, "-1.0", "21", NODEF, NODEF, NODEF, "0"),
      GIN_DESC("Effective size of the antenna,",
               "Full width of grid.",
               "Number of points to represent the grid.",
               "Number of points in frequency grid (must be >= 2)",
               "Start point of frequency grid",
               "End point of frequency grid",
               "Set to 1 to create a 2D antenna pattern.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Append"),
      DESCRIPTION(R"--(Append one workspace variable to another.

This method can append an array to another array of the same type,
e.g. ArrayOfIndex to ArrayOfIndex. Or a single element to an array
such as a Tensor3 to an ArrayOfTensor3.

Appending two vectors or a numeric to a vector works as for array
variables.

Both another matrix or a vector can be appended to a matrix. In
addition, for matrices, the 'append dimension' can be selected.
The third argument, ``dimension``, indicates how to append, where
\"leading\" means to append row-wise, and \"trailing\" means
column-wise.

Other types (TensorX) are currently only implemented for
appending to the leading dimension.

This method is not implemented for all types, just for those that
were thought or found to be useful. (See variable list below.).
)--"),
      AUTHORS("Stefan Buehler, Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE(("Vector, Vector,"
                 "Matrix, Matrix,"
                 "Tensor3, Tensor3,"
                 "Tensor4, Tensor4,"
                 "String, " +
                 ARRAY_GROUPS + ", " + ARRAY_GROUPS_WITH_BASETYPE)
                    .c_str()),
      GOUT_DESC("The variable to append to."),
      IN(),
      GIN("input", "dimension"),
      GIN_TYPE(("Numeric, Vector,"
               "Matrix, Vector,"
               "Matrix, Tensor3,"
               "Tensor3, Tensor4,"
               "String, " +
                   ARRAY_GROUPS + "," + GROUPS_WITH_ARRAY_TYPE).c_str(),
               "String"),
      GIN_DEFAULT(NODEF, "leading"),
      GIN_DESC(
          "The variable to append.",
          "Where to append. Could be either the \"leading\" or \"trailing\" dimension."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("ArrayOfGriddedFieldGetNames"),
      DESCRIPTION(R"--(Get the names of all GriddedFields stored in an Array.

See *GriddedFieldGetName*.
)--"),
      AUTHORS("Lukas Kluft"),
      OUT(),
      GOUT("names"),
      GOUT_TYPE("ArrayOfString"),
      GOUT_DESC("Names of the GriddedFields in the ArrayOfGriddedField."),
      IN(),
      GIN("griddedfields"),
      GIN_TYPE("ArrayOfGriddedField1, ArrayOfGriddedField2,"
               "ArrayOfGriddedField3, ArrayOfGriddedField4"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Array of GriddedFields."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("ArrayOfIndexLinSpace"),
      DESCRIPTION(R"--(Initializes an ArrayOfIndex with linear spacing.

The first element equals always the start value, and the spacing
equals always the step value, but the last value can deviate from
the stop value. ``step`` can be both positive and negative.

The created array is [start, start+step, start+2*step, ...]
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("ArrayOfIndex"),
      GOUT_DESC("Output array."),
      IN(),
      GIN("start", "stop", "step"),
      GIN_TYPE("Index", "Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Start value.",
               "Maximum/minimum value of the end value",
               "Spacing of the array.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ArrayOfIndexSetConstant"),
      DESCRIPTION(R"--(Creates an ArrayOfIndex of length *nelem*, with all values
identical.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("ArrayOfIndex"),
      GOUT_DESC("Variable to initialize."),
      IN("nelem"),
      GIN("value"),
      GIN_TYPE("Index"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Array value.."),
      SETMETHOD(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("ArrayOfTimeNLinSpace"),
      DESCRIPTION(R"--(Creates a time array with length *nelem*, equally spaced between the
given end values.

The length (*nelem*) must be larger than 1.
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("ArrayOfTime"),
      GOUT_DESC("Variable to initialize."),
      IN("nelem"),
      GIN("start", "stop"),
      GIN_TYPE("String", "String"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Start value.", "End value.")));

  md_data_raw.push_back(create_mdrecord(
    NAME("ArrayOfTimeSetConstant"),
      DESCRIPTION(R"--(Creates an ArrayOfTime and sets all elements to the specified value.

The vector length is determined by *nelem*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("ArrayOfTime"),
      GOUT_DESC("Variable to initialize."),
      IN("nelem"),
      GIN("value"),
      GIN_TYPE("Time"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Time value.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Arts"),
      DESCRIPTION(R"--(Runs the agenda that is specified inside the curly braces. ARTS
controlfiles must define this method. It is executed automatically
when ARTS is run on the controlfile and cannot be called by the user.
This methods was used for Arts 1 controlfiles and is now obsolete.
See *Arts2*
)--"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC(),
      SETMETHOD(false),
      AGENDAMETHOD(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("Arts2"),
      DESCRIPTION(R"--(Runs the agenda that is specified inside the curly braces. ARTS
controlfiles must define this method. It is executed automatically
when ARTS is run on the controlfile and cannot be called by the user.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC(),
      SETMETHOD(false),
      AGENDAMETHOD(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("AtmFieldPRegrid"),
      DESCRIPTION(R"--(Interpolates the input field along the pressure dimension from
``p_grid_old`` to to ``p_grid_new``.

Extrapolation is allowed within the common 0.5grid-step margin.
in and out fields can be the same variable.
)--"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor3, Tensor4"),
      GOUT_DESC("Regridded atmospheric field."),
      IN(),
      GIN("input", "p_grid_new", "p_grid_old", "interp_order"),
      GIN_TYPE("Tensor3, Tensor4", "Vector", "Vector", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "1"),
      GIN_DESC("Input atmospheric field.",
               "Pressure grid to regrid to",
               "Pressure grid of input field",
               "Interpolation order.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("AtmFieldsAndParticleBulkPropFieldFromCompact"),
      DESCRIPTION(R"--(Extract pressure grid and atmospheric fields from
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
)--"),
      AUTHORS("Jana Mendrok, Manfred Brath"),
      OUT(
          "atm_field",
          "particle_bulkprop_names"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species", "atm_fields_compact"),
      GIN("delim", "check_gridnames"),
      GIN_TYPE("String", "Index"),
      GIN_DEFAULT("-", "0"),
      GIN_DESC(/* delim */
               "Delimiter string of *scat_species* elements.",
               /* check_gridnames */
               "A flag with value 1 or 0. If set to one, the gridnames of"
               " the *atm_fields_compact* are checked.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atmfields_checkedCalc"),
      DESCRIPTION(R"--(Checks consistency of (clear sky) atmospheric fields.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("atmfields_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         "abs_species",
         "atm_field"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("atmgeom_checkedCalc"),
      DESCRIPTION(R"--(Checks consistency of geometric considerations of the atmosphere.

The following WSVs are checked: ``z_field``, ``refellipsoid``, ``z_surface``,
*lat_true* and *lon_true*. If any of the variables above is changed,
then this method shall be called again (no automatic check that this is
fulfilled!).

The tests include that:
 1. ``refellipsoid`` has correct size, and that eccentricity is
    set to zero if 1D atmosphere.
 2. ``z_field`` and ``z_surface`` have sizes consistent with the
    atmospheric grids.
 3. There is no gap between ``z_surface`` and ``z_field``.
 4. A rough search of maximum gradient of the altitude of the pressure
    level closest to 500 hPa is made. If this value exceeds the GIN
    ``max500hpa_gradient`` an error is issued. Please note that the unit
    of this GIN is m per 100km. For normal conditions on Earth, large
    scale gradients of the 500 hPa level is in the order of 20m/100km.

*lat_true* and *lon_true* are allowed to be empty.

If any test fails, there is an error. Otherwise, *atmgeom_checked*
is set to 1.

See further *atmgeom_checkedCalc*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("atmgeom_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         "atm_field",
         "surface_field",
         "lat_true",
         "lon_true"),
      GIN("max500hpa_gradient"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("500"),
      GIN_DESC("The maximum allowed gradient of 500 hPa pressure level [m/100km].")));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fieldLteInternalPartitionFunction"),
      DESCRIPTION(R"--(Turns on NTLE calculations.

Sets NLTE ratios to those expected for LTE calculations
with estimation of the partition function as the sum of all
states of a species
)--"),
      AUTHORS("Richard Larsson"),
      OUT("nlte_do", "atm_field", "abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field", "abs_lines_per_species"),
      GIN("nlte_level_identifiers"),
      GIN_TYPE("ArrayOfQuantumIdentifier"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("List of levels to compute for")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fieldLteExternalPartitionFunction"),
      DESCRIPTION(R"--(Turns on NTLE calculations.

Sets NLTE ratios to those expected for LTE calculations
with a known partition function
)--"),
      AUTHORS("Richard Larsson"),
      OUT("nlte_do", "atm_field", "abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field", "abs_lines_per_species"),
      GIN("nlte_level_identifiers"),
      GIN_TYPE("ArrayOfQuantumIdentifier"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("List of levels to compute for")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fields_compactAddConstant"),
      DESCRIPTION(R"--(Adds a constant field to atm_fields_compact.

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


For Earth this should be set to [\"abs_species-H2O\"]
)--"),
      AUTHORS("Stefan Buehler, Oliver Lemke"),
      OUT("atm_fields_compact"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_fields_compact"),
      GIN("name", "value", "prepend", "condensibles"),
      GIN_TYPE("String", "Numeric", "Index", "ArrayOfString"),
      GIN_DEFAULT(NODEF, NODEF, "0", "[]"),
      GIN_DESC(
          "Name of additional atmospheric field, with constant value.",
          "Constant value of additional field.",
          "0 = Append to the end, 1 = insert at the beginning.",
          "List of condensibles used to scale down the VMR of the added species.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fields_compactAddSpecies"),
      DESCRIPTION(R"--(Adds a field to atm_fields_compact, with interpolation.

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
)--"),
      AUTHORS("Gerrit Holl"),
      OUT("atm_fields_compact"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_fields_compact"),
      GIN("name", "value", "prepend"),
      GIN_TYPE("String", "GriddedField3", "Index"),
      GIN_DEFAULT(NODEF, NODEF, "0"),
      GIN_DESC("Name of additional atmospheric field.",
               "Value of additional atmospheric field.",
               "0 = Append to the end, 1 = insert at the beginning.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fields_compactCleanup"),
      DESCRIPTION(R"--(Removes unrealistically small or erroneous data from
*atm_fields_compact* (or other GriddedField4 data)

This WSM checks if the data in *atm_fields_compact* contains
values smaller than the given ``threshold``. In this case, these
values will be set to zero.

The method should be applied if *atm_fields_compact* contains
unrealistically small or erroneous data (NWP/GCM model data
occassionally contains negative values, which are numerical
artefacts rather than physical values.)
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("atm_fields_compact"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_fields_compact"),
      GIN("threshold"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC(
          "Threshold below which *atm_fields_compact* values are set to zero.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fields_compactCreateFromField"),
      DESCRIPTION(R"--(Initiates *atm_fields_compact* from a field.

*atm_fields_compact* will have the same size and grids as the GriddedField3,
but with one dimension as length 1.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("atm_fields_compact"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("name", "field"),
      GIN_TYPE("String", "GriddedField3"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Name atmospheric field.", "The atmospheric field.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fields_compactFromMatrix"),
      DESCRIPTION(R"--(Sets *atm_fields_compact* from 1D fields given in form of a matrix.

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
)--"),
      AUTHORS("Stefan Buehler", "Daniel Kreyling", "Jana Mendrok"),
      OUT("atm_fields_compact"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("gin1", "field_names"),
      GIN_TYPE("Matrix", "ArrayOfString"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("One atmosphere matrix from batch input ArrayOfMatrix.",
               "Order/names of atmospheric fields.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fieldTopOfAtmosphere"),
      DESCRIPTION(R"--(Sets the top of the atmosphere altitude to the field
)--"),
      AUTHORS("Richard Larsson"),
      OUT("atm_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field"),
      GIN("toa"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Top of atmosphere altitude [m].")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fieldInit"),
      DESCRIPTION(R"--(Initialize the atmospheric field with some altitude
)--"),
      AUTHORS("Richard Larsson"),
      OUT("atm_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("toa"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Top of atmosphere altitude [m].")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fieldAddCustomDataFile"),
      DESCRIPTION(R"--(Add some custom data from file to the atm_field

The key field is used to determine the type of data that is added by input type.

If the input is a String, the data is added to corresponding atmospheric data,
these strings can be
    "t"      - temperature
    "p"      - pressure
    "wind_u" - wind u component
    "wind_v" - wind v component
    "wind_w" - wind w component
    "mag_u"  - mag u component
    "mag_v"  - mag v component
    "mag_w"  - mag w component

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
)--"),
      AUTHORS("Richard Larsson"),
      OUT("atm_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field"),
      GIN("key", "filename", "extrapolation_type"),
      GIN_TYPE("String,ArrayOfSpeciesTag,QuantumIdentifier", "String", "String"),
      GIN_DEFAULT(NODEF, NODEF, "Nearest"),
      GIN_DESC("Atmospheric data key.", "Filename", "Style of extrapolation")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fieldAddField"),
      DESCRIPTION(R"--(Add another atm_field from file to the current atm_field

The optional flag set_toa determines if the old (default) or
new (if it evaluates as true) atm_field's top of the atmosphere altitude
is used in the output
)--"),
      AUTHORS("Richard Larsson"),
      OUT("atm_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field"),
      GIN("filename", "set_toa"),
      GIN_TYPE("String", "Index"),
      GIN_DEFAULT(NODEF, "0"),
      GIN_DESC("Filename", "Flag for overwriting the top of the atmosphere")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fieldRead"),
      DESCRIPTION(R"--(Reads a new atm_field from a folder or base

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
    read_tp - ["t.xml", "p.xml", ]
    read_mag - ["mag_u.xml", "mag_v.xml", "mag_w.xml", ]
    read_wind - ["wind_u.xml", "wind_v.xml", "wind_w.xml", ]
    read_specs - [See below]
    read_nlte - "nlte.xml"

If "read_specs" is true, then all the species of *abs_species* are read and the
basename path is expected to contain a file with short-name version for each
unique species.  Some examples:
    abs_species=["H2O-161", "O2-66"], - ["H2O.xml", "O2.xml"]
    abs_species=["H2O-161", "O2-66", "CO2-626"], - ["H2O.xml", "O2.xml", "CO2.xml"]
    abs_species=["H2O-161", "O2-66", "O2-PWR98"], - ["H2O.xml", "O2.xml"]
)--"),
      AUTHORS("Richard Larsson"), OUT("atm_field"), GOUT(), GOUT_TYPE(),
      GOUT_DESC(), IN("abs_species"),
      GIN("basename", "toa", "read_tp", "read_mag", "read_wind",
          "read_specs", "read_nlte"),
      GIN_TYPE("String", "Numeric", "Index", "Index", "Index", "Index",
               "Index"),
      GIN_DEFAULT("./", NODEF, "1", "0", "0", "1", "0"),
      GIN_DESC("Base for the name of the data files.",
               "Top of atmosphere altitude [m].",
               "Flag to read pressure and temperature",
               "Flag to read magnetic field", "Flag to read wind field",
               "Flag to read species", "Flag to read NLTE")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fieldSave"),
      DESCRIPTION(R"--(Saves an atm_field to a folder or base

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
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field"),
      GIN("basename", "filetype", "no_clobber"),
      GIN_TYPE("String", "String", "Index"),
      GIN_DEFAULT(NODEF, "ascii", "0"),
      GIN_DESC("Base for the name of the data files.", "See *WriteXML*", "See *WriteXML*")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fieldAddGriddedData"),
      DESCRIPTION(R"--(Adds data to the atm_field

The field must not be regular
)--"),
      AUTHORS("Richard Larsson"),
      OUT("atm_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field"),
      GIN("key", "data", "extrapolation_type"),
      GIN_TYPE("String,ArrayOfSpeciesTag,QuantumIdentifier", "GriddedField3", "String"),
      GIN_DEFAULT(NODEF, NODEF, "Nearest"),
      GIN_DESC("See *atm_fieldAddCustomDataFile*", "Some data", "Style of extrapolation")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fieldAddNumericData"),
      DESCRIPTION(R"--(Adds data to the atm_field

The field must not be regular
)--"),
      AUTHORS("Richard Larsson"),
      OUT("atm_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field"),
      GIN("key", "data"),
      GIN_TYPE("String,ArrayOfSpeciesTag,QuantumIdentifier", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("See *atm_fieldAddCustomDataFile*", "Some data")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fieldIGRF"),
      DESCRIPTION(R"--(Use IGRF to compute the magnetic field at each point

The flag ``parsafe`` exists if you need the calculations to be safe in parallel
computations.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("atm_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field", "time"),
      GIN("parsafe"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Flag for parallel safety at 3X slowdown cost")));

  md_data_raw.push_back(create_mdrecord(
      NAME("backend_channel_responseFlat"),
      DESCRIPTION(R"--(Sets up a rectangular channel response.

The method assumes that all channels have the same response.

The response of the backend channels is hee assumed to be constant
inside the resolution width, and zero outside.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("backend_channel_response"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("resolution"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The spectrometer resolution.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("backend_channel_responseGaussian"),
      DESCRIPTION(R"--(Sets up a Gaussian backend channel response.

The method assumes that all channels have the same response.

This method works as *backend_channel_responseGaussianConstant*
but handles the case where the response of each channel must be
described. Here the FWHM is specified for each *f_backend*.

The GINs ``fwhm`` and ``grid_npoints`` work in the same way as for
*antenna_responseGaussianConstant*. A negative ``grid_width``
gives a grid that is twice the FWHM of each channel.
)--"),
      AUTHORS("Patrick Eriksson, Oliver Lemke"),
      OUT("backend_channel_response"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_backend"),
      GIN("fwhm", "grid_width", "grid_npoints"),
      GIN_TYPE("Vector", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, "-1.0", "21"),
      GIN_DESC("Full width at half-maximum of the Gaussian function.",
               "Full width of grid.",
               "Number of points to represent the grid.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("backend_channel_responseGaussianConstant"),
      DESCRIPTION(R"--(Sets up a single Gaussian backend channel response.

The method assumes that all channels have the same response.

The GINs ``fwhm`` and ``grid_npoints`` work in the same way as for
*antenna_responseGaussianConstant*. A negative ``grid_width``
gives a grid that is twice the FWHM.
)--"),
      AUTHORS("Patrick Eriksson, Oliver Lemke"),
      OUT("backend_channel_response"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("fwhm", "grid_width", "grid_npoints"),
      GIN_TYPE("Numeric", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, "-1.0", "21"),
      GIN_DESC("Full width at half-maximum of the Gaussian function.",
               "Full width of grid.",
               "Number of points to represent the grid.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("batch_atm_fields_compactAddConstant"),
      DESCRIPTION(R"--(Adds a constant field to batch_atm_fields_compact.

Applies *atm_fields_compactAddConstant* to each batch.
The format is equal to that WSM.
)--"),
      AUTHORS("Gerrit Holl"),
      OUT("batch_atm_fields_compact"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("batch_atm_fields_compact"),
      GIN("name", "value", "prepend", "condensibles"),
      GIN_TYPE("String", "Numeric", "Index", "ArrayOfString"),
      GIN_DEFAULT(NODEF, NODEF, "0", "[]"),
      GIN_DESC(
          "Name of additional atmospheric field, with constant value.",
          "Constant value of additional field.",
          "0 = Append to the end, 1 = insert at the beginning.",
          "List of condensibles used to scale down the VMR of the added species.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("batch_atm_fields_compactAddSpecies"),
      DESCRIPTION(R"--(Adds a field to *batch_atm_fields_compact*, with interpolation.

This method appends or prepends a *GriddedField3* to each *atm_fields_compact*.
in *batch_atm_fields_compact*. For details, see *atm_fields_compactAddSpecies*.
)--"),
      AUTHORS("Gerrit Holl"),
      OUT("batch_atm_fields_compact"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("batch_atm_fields_compact"),
      GIN("name", "value", "prepend"),
      GIN_TYPE("String", "GriddedField3", "Index"),
      GIN_DEFAULT(NODEF, NODEF, "0"),
      GIN_DESC(
          "Name of additional atmospheric field. Use, e.g., vmr_ch4 for methane VMR",
          "Value of additional atmospheric field.",
          "0 = Append to the end, 1 = insert at the beginning.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("batch_atm_fields_compactCleanup"),
      DESCRIPTION(R"--(Removes unrealistically small or erroneous data from each data field
of *batch_atm_fields_compact* (or other AerrayOfGriddedField4 data)

This WSM checks if the data in *batch_atm_fields_compact* contains
values smaller than the given ``threshold``. In this case, these
values will be set to zero.

The method should be applied if *batch_atm_fields_compact* contains
unrealistically small or erroneous data (NWP/GCM model data
occassionally contains negative values, which are numerical
artefacts rather than physical values.)
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("batch_atm_fields_compact"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("batch_atm_fields_compact"),
      GIN("threshold"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC(
          "Threshold below which *atm_fields_compact* values are set to zero.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("batch_atm_fields_compactFromArrayOfMatrix"),
      DESCRIPTION(R"--(Expand batch of 1D atmospheric state matrices to batch_atm_fields_compact.

This is used to handle 1D batch cases, e.g. from NWP/GCM model like
the Chevallier91L data set, stored in a matrix (it is preferred,
though, to immediatedly store the model fields as
*ArrayOfGriddedField4* and use *ReadXML* to load them directly into
*batch_atm_fields_compact*).

Works only for ``atmosphere_dim`` == 1.

See *atm_fields_compactFromMatrix* for basic documentation.

See *batch_atm_fields_compactAddConstant* and
*batch_atm_fields_compactAddSpecies* for adding additional fields.
)--"),
      AUTHORS("Stefan Buehler", "Daniel Kreyling", "Jana Mendrok"),
      OUT("batch_atm_fields_compact"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("atmospheres_fields", "field_names"),
      GIN_TYPE("ArrayOfMatrix", "ArrayOfString"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Batch of atmospheres stored in one array of matrix",
               "Order/names of atmospheric fields.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("CIAInfo"),
      DESCRIPTION(R"--(Display information about the given CIA tags.
The CIA tags shown are in the same format as needed by *abs_speciesSet*.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("catalogpath", "cia_tags"),
      GIN_TYPE("String", "ArrayOfString"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Path to the CIA catalog directory.",
               "Array of CIA tags to view, e.g. [ \"N2-N2\", \"H2-H2\" ]")));

  md_data_raw.push_back(create_mdrecord(
      NAME("CIARecordReadFromFile"),
      DESCRIPTION(R"--(Reads CIARecord from Hitran-style file.
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("cia_record"),
      GOUT_TYPE("CIARecord"),
      GOUT_DESC("CIARecord type variable for input and output."),
      IN(),
      GIN("species_tag", "filename"),
      GIN_TYPE("String", "String"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("SpeciesTag string to associate with this CIARecord. See "
               "*abs_speciesSet* for correct format.",
               "Filename of HITRAN CIA data file.")));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("cloudboxOff"),
      DESCRIPTION(R"--(Deactivates the cloud box.

Use this method if no scattering calculations shall be performed.

The function sets *cloudbox_on* to 0, *cloudbox_limits*,
*pnd_field*, *scat_data*, *scat_data_raw*, *iy_cloudbox_agenda*
and *particle_masses* to be empty and sizes *dpnd_field_dx* to be
consitent with *jacobian_quantities*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("cloudbox_on",
          "ppath_inside_cloudbox_do",
          "cloudbox_limits",
          "iy_cloudbox_agenda",
          "pnd_field",
          "dpnd_field_dx",
          "scat_species",
          "scat_data",
          "scat_data_raw",
          "scat_data_checked",
          "particle_masses"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC(),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true),
      PASSWSVNAMES(false)));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudboxSetAutomatically"),
      DESCRIPTION(R"--(Sets the cloud box to encompass the cloud given by the entries
in ``particle_field``.

This WSM handles one *Tensor4* type ``particle_field`` at a time. It can
be used to determine the cloudbox from ``particle_bulkprop_field``

The function must be called before executing any WSM that applies
*cloudbox_limits*.

The function iterates over all 3D fields in ``particle_field`` (which
might correspond to different particle bulk properties as in
``particle_bulkprop_field``). Each field is searched for the first
and last pressure index, where the value is unequal to zero. This
index is then copied to *cloudbox_limits*.
If ``particle_field`` is empty, the cloudbox is switched off
(*cloudbox_on* = 0).

Additionaly the lower cloudbox_limit is altered by ``cloudbox_margin``.
The margin is given as a height difference in meters and transformed
into a pressure (via isothermal barometric height formula). This
alteration is to ensure covering photons that leave the cloud, but
reenter through a limb path.
If ``cloudbox_margin`` is set to -1 (default), the cloudbox will extend
to the surface. Hence, the lower cloudbox_limit is set to 0 (index
of first pressure level).
``cloudbox_margin`` will be applied on each call of the WSM.

Works only for ``atmosphere_dim`` == 1.
)--"),
      AUTHORS("Jana Mendrok, Daniel Kreyling"),
      OUT("cloudbox_on", "cloudbox_limits"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("p_grid", "lat_grid", "lon_grid"),
      GIN("particle_field", "cloudbox_margin"),
      GIN_TYPE("Tensor4", "Numeric"),
      GIN_DEFAULT(NODEF, "-1"),
      GIN_DESC("A collection of particle property fields (e.g."
               " ``particle_bulkprop_field``).",
               "Minimum distance [m] between lowest 'cloudy' level and"
               " cloudbox lower limit. If set to ``-1`` (default), the"
               " cloudbox lower limit is fixed to 0, i.e., corresponds to"
               " the lowest atmospheric level (or the surface).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudboxSetFullAtm"),
      DESCRIPTION(R"--(Sets the cloudbox to cover the full atmosphere.

The cloudbox is always set to fully span the atmosphere vertically.

For the latitide and longitide dimensions, default is to leave room
between the cloudbox an the end of the atmosphere in these dimensions.
This is required for some scattering solvers (MC and DOIT). In other
cases it can be OK to let the cloudbox to fill the atmosphere fully
also in latitude and longitude. This is triggered by setting the GIN
fullfull to 1.
)--"),
      AUTHORS("Claudia Emde, Jana Mendrok, Patrick Eriksson"),
      OUT("cloudbox_on", "cloudbox_limits"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field"),
      GIN("fullfull"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("Flag to let cloudbox reach ends of lat and lon grids.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudboxSetManually"),
      DESCRIPTION(R"--(Sets the cloud box to encompass the given positions.

The function sets *cloudbox_on* to 1 and sets *cloudbox_limits*
following the given pressure, latitude and longitude positions.
The index limits in *cloudbox_limits* are selected to give the
smallest possible cloud box that encompass the given points.

The points must be given in the same order as used in
*cloudbox_limits*. That means that the first keyword argument
shall be a higher pressure than argument two, while the latitude
and longitude points are given in increasing order. Positions
given for dimensions not used by the selected atmospheric
dimensionality are ignored.

The given pressure points can be outside the range of ``p_grid``.
The pressure limit is then set to the end point of ``p_grid``.
The given latitude and longitude points must be inside the range
of the corresponding grid. In addition, the latitude and longitude
points cannot be inside the outermost grid ranges as the latitude
and longitude limits in *cloudbox_limits* are not allowed to be
grid end points.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("cloudbox_on", "cloudbox_limits"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("p_grid", "lat_grid", "lon_grid"),
      GIN("p1", "p2", "lat1", "lat2", "lon1", "lon2"),
      GIN_TYPE(
          "Numeric", "Numeric", "Numeric", "Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, NODEF, NODEF),
      GIN_DESC("Upper pressure point.",
               "Lower pressure point.",
               "Lower latitude point.",
               "Upper latitude point.",
               "Lower longitude point.",
               "Upper longitude point.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudboxSetManuallyAltitude"),
      DESCRIPTION(R"--(Sets the cloud box to encompass the given positions.

As ``cloudboxSetManually`` but uses altitudes instead of pressure.
The given altitude points can be outside the range of ``z_field``.
The altitude limit is then set to the end point of ``p_grid``.
)--"),
      AUTHORS("Claudia Emde"),
      OUT("cloudbox_on", "cloudbox_limits"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field"),
      GIN("z1", "z2", "lat1", "lat2", "lon1", "lon2"),
      GIN_TYPE(
          "Numeric", "Numeric", "Numeric", "Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, NODEF, NODEF),
      GIN_DESC("Lower altitude point.",
               "Upper altitude point.",
               "Lower latitude point.",
               "Upper latitude point.",
               "Lower longitude point.",
               "Upper longitude point.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_checkedCalc"),
      DESCRIPTION(R"--(Checks consistency and validity of the cloudbox governing variables.

The following WSVs are treated: *cloudbox_on*, *cloudbox_limits*,
*pnd_field*, *scat_data*, *scat_species*, *abs_species*, *particle_masses*
``particle_bulkprop_field``, *particle_bulkprop_names* and wind_u/v/w_field.

If any of these variables is changed, then this method shall be
called again (no automatic check that this is fulfilled!).

The main checks are if the cloudbox limits are OK with respect to
the atmospheric dimensionality and the limits of the atmosphere,
and that the scattering element variables *pnd_field* and
*scat_data* match in size.

Default is to demand that there is a margin between the cloudbox
and the ends of latitide and longitude grids. Such margins are
required by MC and DOIT/3D, but are not needed for e.g. IBA.
If the margins not are a demand, set GIN demand_latlon_margin to 0.

Further checks on *scat_data* are performed in *scat_data_checkedCalc*

*scat_species* and *particle_masses* must either be empty or have a
size that matches the other data. If non-empty, some check of these
variables are performed.

If any test fails, there is an error. Otherwise, *cloudbox_checked*
is set to 1.
)--"),
      AUTHORS("Patrick Eriksson, Jana Mendrok"),
      OUT("cloudbox_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmfields_checked",
         "atm_field",
         "surface_field",
         "cloudbox_on",
         "cloudbox_limits",
         "pnd_field",
         "dpnd_field_dx",
         "jacobian_quantities",
         "scat_data",
         "scat_species",
         "particle_masses",
         "abs_species"),
      GIN("demand_latlon_margin", "negative_pnd_ok"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT("1", "0"),
      GIN_DESC("Flag to demand margin or not w.r.t. to ends of lat/lon grids.",
               "Flag whether to accept pnd_field < 0.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_field_monoIterate"),
      DESCRIPTION(R"--(Iterative solution of the VRTE (DOIT method).

A solution for the RTE with scattering is found using the
DOIT method:

1. Calculate scattering integral using *doit_scat_field_agenda*.
2. Calculate RT with fixed scattered field using
   *doit_rte_agenda*.
3. Convergence test using *doit_conv_test_agenda*.

Note:
      The atmospheric dimensionality ``atmosphere_dim`` can be
      either 1 or 3. To these dimensions the method adapts
      automatically. 2D scattering calculations are not
      supported.
)--"),
      AUTHORS("Claudia Emde, Jakob Doerr"),
      OUT("cloudbox_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field_mono",
         "doit_scat_field_agenda",
         "doit_rte_agenda",
         "doit_conv_test_agenda"),
      GIN("accelerated"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC(
          "Index wether to accelerate only the intensity (1) or the whole Stokes Vector (4)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_fieldCrop"),
      DESCRIPTION(R"--(Extracts a part of an existing *cloudbox_field*.

The cropping is defined by defining new cloudbox limits. Note that
``new_limit0`` is an index with respect to ``p_grid``, etc.

The following must be valid:
 * new_limit0 >= cloudbox_limits[0]
 * new_limit1 <= cloudbox_limits[1]
 * new_limit2 >= cloudbox_limits[2]
 * new_limit3 <= cloudbox_limits[3]
 * new_limit4 >= cloudbox_limits[4]
 * new_limit5 <= cloudbox_limits[5]

Indexes for dimensions not used are ignored.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("cloudbox_field", "cloudbox_limits"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_on", "cloudbox_limits", "cloudbox_field"),
      GIN("new_limit0",
          "new_limit1",
          "new_limit2",
          "new_limit3",
          "new_limit4",
          "new_limit5"),
      GIN_TYPE("Index", "Index", "Index", "Index", "Index", "Index"),
      GIN_DEFAULT("0", "0", "0", "0", "0", "0"),
      GIN_DESC("New value for cloudbox_limits[0].",
               "New value for cloudbox_limits[1].",
               "New value for cloudbox_limits[2].",
               "New value for cloudbox_limits[3].",
               "New value for cloudbox_limits[4].",
               "New value for cloudbox_limits[5].")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_fieldInterp2Azimuth"),
      DESCRIPTION(R"--(Reinterpolate a *cloudbox_field* with azimuthal dependency.

Intended use: Call directly after cloudbox_fieldDisort if sun is present and yCalc should be
should be run afterwards.

In ARTS a 1D atmosphere cannot have a azimuth dependency, but if a 
collimated source like a sun is present even a 1D atmosphere has an 
azimuth dependency. To overcome this constraint, the user must set an 
additional local sensor line of sight azimuth angle for the true
geopgraphical location of the atmosphere. For this angle the 
*cloudbox_field* with azimuthal dependency is interpolated to a 
*cloudbox_field* without azimuthal dependency
)--"),
      AUTHORS("Manfred Brath"),
      OUT("cloudbox_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field","cloudbox_on", "aa_grid"),
      GIN("local_los_azimuth_angle","aa_interp_order"),
      GIN_TYPE("Numeric","Index"),
      GIN_DEFAULT(NODEF,"1"),
      GIN_DESC("Local line of sight azimuth angle",
               "Azimuth angle interpolation order.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_fieldSetFromPrecalc"),
      DESCRIPTION(R"--(Sets the initial cloudbox intensity field *cloudbox_field* from a
precalculated field.

This method sets the (monochromatic) first guess radiation field
inside the cloudbox from a precalculated ``cloudbox_field_precalc``,
e.g., from the solution of a similar atmospheric scenario. The
dimensions of ``cloudbox_field_precalc`` have to be consistent with
the DOIT setup in terms of frequencies, pressure levels inside the
cloudbox, polar angles used as well as the stokes dimension.
Incoming field on the cloudbox boundaries is adapted to the actual
clearsky incoming field as, e.g., calculated by *DoitGetIncoming*.
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("cloudbox_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field",
         "za_grid",
         "f_grid",
         
         "cloudbox_limits",
         "doit_is_initialized"),
      GIN("cloudbox_field_precalc"),
      GIN_TYPE("Tensor7"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Precalculated radiation field (of type *cloudbox_field*)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_fieldSetClearsky"),
      DESCRIPTION(R"--(Interpolate clearsky field on all gridpoints in cloudbox.

This method uses a linear 1D/3D interpolation scheme to obtain the
radiation field on all grid points inside the cloud box from the
clear sky field on the cloudbox boundary. This radiation field
is taken as the first guess radiation field in the DOIT module.

Set the ``all_frequencies`` to 1 if the clearsky field shall be used
as initial field for all frequencies. Set it to 0 if the clear sky
field shall be used only for the first frequency in *f_grid*. For
later frequencies, *cloudbox_field* of the previous frequency is then
used.
)--"),
      AUTHORS("Sreerekha T.R. and Claudia Emde"),
      OUT("cloudbox_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "cloudbox_limits",
         "cloudbox_on",
         "doit_is_initialized"),
      GIN("all_frequencies"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("See above.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_field_monoSetConst"),
      DESCRIPTION(R"--(This method sets the initial field inside the cloudbox to a
constant value. The method works only for monochromatic
calculations (number of elements in f_grid=1).

The user can specify a value for each Stokes dimension in the
control file by ``value``.
)--"),
      AUTHORS("Claudia Emde"),
      OUT("cloudbox_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field_mono",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "cloudbox_limits",
         "stokes_dim"),
      GIN("value"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("A vector containing 4 elements with the value of the "
               "initial field for each Stokes dimension.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_fieldSetConst"),
      DESCRIPTION(R"--(This method sets the initial field inside the cloudbox to a
constant value.

The user has to specify a value for each Stokes dimension in the
control file by ``value``.
)--"),
      AUTHORS("Claudia Emde"),
      OUT("cloudbox_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "cloudbox_limits",
         "stokes_dim"),
      GIN("value"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("A vector containing ``stokes_dim`` elements with the value of"
               " the initial field for each Stokes dimension.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_fieldSetConstPerFreq"),
      DESCRIPTION(R"--(This method sets the initial field inside the cloudbox to a
constant value per frequency slice.

The user has specify a value for each frequency and Stokes
dimension in the control file by ``value``.
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("cloudbox_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "cloudbox_limits",
         "stokes_dim"),
      GIN("value"),
      GIN_TYPE("Matrix"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("A matrix containing ``stokes_dim`` elements per frequency"
               " (row) with the value of the initial field for each"
               " frequency and Stokes dimension.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_fieldUpdate1D"),
      DESCRIPTION(R"--(RT calculation in cloudbox with fixed scattering integral (1D).

Updates the radiation field (DOIT method). The method loops
through the cloudbox to update the radiation field for all
positions and directions in the 1D cloudbox.

Note: This method is very inefficient, because the number of
iterations scales with the number of cloudbox pressure levels.
It is recommended to use *cloudbox_fieldUpdateSeq1D*.
)--"),
      AUTHORS("Claudia Emde"),
      OUT("cloudbox_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field_mono",
         "doit_scat_field",
         "cloudbox_limits",
         "propmat_clearsky_agenda",
         "atm_field",
         "surface_field".
         "abs_species",
         "spt_calc_agenda",
         "za_grid",
         "pnd_field",
         "ppath_step_agenda",
         "ppath_lmax",
         "ppath_lraytrace",
         "f_grid",
         "f_index",
         "surface_rtprop_agenda",
         "doit_za_interp"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_fieldUpdateSeq1D"),
      DESCRIPTION(R"--(RT calculation in cloudbox with fixed scattering integral.

Updates radiation field (*cloudbox_field*) in DOIT module.
This method loops through the cloudbox to update the
radiation field for all positions and directions in the 1D
cloudbox. The method applies the sequential update. For more
information refer to AUG.
)--"),
      AUTHORS("Claudia Emde"),
      OUT("cloudbox_field_mono", "doit_scat_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field_mono",
         "doit_scat_field",
         "cloudbox_limits",
         "propmat_clearsky_agenda",
         "atm_field",
         "abs_species",
         "spt_calc_agenda",
         "za_grid",
         "aa_grid",
         "pnd_field",
         "ppath_step_agenda",
         "ppath_lmax",
         "ppath_lraytrace",
         "f_grid",
         "f_index",
         "surface_rtprop_agenda",
         "doit_za_interp"),
      GIN("normalize", "norm_error_threshold", "norm_debug"),
      GIN_TYPE("Index", "Numeric", "Index"),
      GIN_DEFAULT("1", "1.0", "0"),
      GIN_DESC(
          "Apply normalization to scattered field.",
          "Error threshold for scattered field correction factor.",
          "Debugging flag. Set to 1 to output normalization factor to out0.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_fieldUpdateSeq1DPP"),
      DESCRIPTION(R"--(RT calculation in cloudbox with fixed scattering integral.

Update radiation field (*cloudbox_field*) in DOIT module.
This method loops through the cloudbox to update the
radiation field for all
positions and directions in the 1D cloudbox. The method applies
the sequential update and the plane parallel approximation.
This method is only slightly faster than
*cloudbox_fieldUpdateSeq1D* and it is less accurate. It can not
be used for limb simulations.
)--"),
      AUTHORS("Sreerekha T.R."),
      OUT("cloudbox_field_mono", "za_index"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field_mono",
         "doit_scat_field",
         "cloudbox_limits",
         "propmat_clearsky_agenda",
         "atm_field",
         "abs_species",
         "spt_calc_agenda",
         "za_grid",
         "pnd_field",
         "f_grid",
         "f_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_fieldUpdateSeq3D"),
      DESCRIPTION(R"--(RT calculation in cloudbox with fixed scattering integral.

Update radiation field (*cloudbox_field*) in DOIT module.
This method loops through the cloudbox to update the
radiation field for all positions and directions in the 3D
cloudbox. The method applies the sequential update. For more
information please refer to AUG.
Surface reflections are not yet implemented in 3D scattering
calculations.
)--"),
      AUTHORS("Claudia Emde"),
      OUT("cloudbox_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field_mono",
         "doit_scat_field",
         "cloudbox_limits",
         "propmat_clearsky_agenda",
         "atm_field",
         "abs_species",
         "spt_calc_agenda",
         "za_grid",
         "aa_grid",
         "pnd_field",
         "ppath_step_agenda",
         "ppath_lmax",
         "ppath_lraytrace",
         "f_grid",
         "f_index",
         "doit_za_interp"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_field_monoOptimizeReverse"),
      DESCRIPTION(R"--(Interpolate *cloudbox_field_mono* back to the original p_grid.

For detailed description, see *OptimizeDoitPressureGrid*.
)--"),
      AUTHORS("Jakob Doerr"),
      OUT("cloudbox_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field_mono", "p_grid_orig", "p_grid", "cloudbox_limits"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("Compare"),
      DESCRIPTION(R"--(Checks the consistency between two variables.

The two variables are checked to not deviate outside the specified
value (``maxabsdiff``). An error is issued if this is not fulfilled.

The main application of this method is to be part of the test
control files, and then used to check that a calculated value
is consistent with an old, reference, value.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("var1", "var2", "maxabsdiff", "error_message"),
      GIN_TYPE(  // INPUT 1
          "Numeric, Vector, Matrix, Tensor3, Tensor4, Tensor5, Tensor7,"
          "ArrayOfVector, ArrayOfMatrix, ArrayOfTensor7, GriddedField3,"
          "Sparse, SingleScatteringData",
          // INPUT 2
          "Numeric, Vector, Matrix, Tensor3, Tensor4, Tensor5, Tensor7,"
          "ArrayOfVector, ArrayOfMatrix, ArrayOfTensor7, GriddedField3,"
          "Sparse, SingleScatteringData",
          // OTHER INPUT
          "Numeric",
          "String"),
      GIN_DEFAULT(NODEF, NODEF, "", ""),
      GIN_DESC("A first variable",
               "A second variable",
               "Threshold for maximum absolute difference.",
               "Additional error message."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("CompareRelative"),
      DESCRIPTION(R"--(Checks the consistency between two variables by their relative values.

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
)--"),
      AUTHORS("Oliver Lemke", "Richard Larsson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("var1", "var2", "maxabsreldiff", "error_message"),
      GIN_TYPE(
          "Numeric, Vector, Matrix, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7,"
          "ArrayOfVector, ArrayOfMatrix, ArrayOfTensor3, ArrayOfTensor4,"
          "ArrayOfTensor6, ArrayOfTensor7, ArrayOfArrayOfVector,"
          "ArrayOfArrayOfMatrix, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6",
          "Numeric, Vector, Matrix, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7,"
          "ArrayOfVector, ArrayOfMatrix, ArrayOfTensor3, ArrayOfTensor4,"
          "ArrayOfTensor6, ArrayOfTensor7, ArrayOfArrayOfVector,"
          "ArrayOfArrayOfMatrix, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6",
          "Numeric",
          "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, ""),
      GIN_DESC("A first variable",
               "A second variable",
               "Threshold for maximum relative difference.",
               "Additional error message."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("complex_refr_indexConstant"),
      DESCRIPTION(R"--(Set complex refractive index to a constant value.

Frequency and temperature grids are set to have length 1 (and
set to the value 0).
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("complex_refr_index"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("refr_index_real", "refr_index_imag"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Real part of refractive index",
               "Imag part of refractive index")));

  md_data_raw.push_back(create_mdrecord(
      NAME("complex_refr_indexIceMatzler06"),
      DESCRIPTION(R"--(Refractive index of ice following Matzler06 parameterization.

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
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("complex_refr_index"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("data_f_grid", "data_T_grid"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Frequency grid for refractive index calculation",
               "Temperature grid for refractive index calculation")));

  md_data_raw.push_back(create_mdrecord(
      NAME("complex_refr_indexTemperatureConstant"),
      DESCRIPTION(R"--(Set frequency dependent complex refractive index.

Temperature grid is set to have length 1 (and
set to the value 0).
)--"),
      AUTHORS("Manfred Brath"),
      OUT("complex_refr_index"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid"),
      GIN("refr_index_real", "refr_index_imag", "temperature"),
      GIN_TYPE("Vector", "Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF,"273.15"),
      GIN_DESC("Real part of refractive index, Dimension [Number of frequencies]",
               "Imag part of refractive index, Dimension [Number of frequencies]",
               "Temperature [K]")));

  md_data_raw.push_back(create_mdrecord(
      NAME("complex_refr_indexWaterLiebe93"),
      DESCRIPTION(R"--(Complex refractive index of liquid water according to Liebe 1993.

The method treats liquid water without salt. Thus, not valid below
10 GHz. Upper frequency limit not known, here set to 1000 GHz.
Model parameters taken from Atmlab function epswater93 (by
C. Maetzler), which refer to Liebe 1993 without closer
specifications.

Temperatures must be between -40 and 100 degrees Celsius. The
accuracy of the parametrization below 0 C is not known by us.
)--"),
      AUTHORS("Patrick Eriksson", "Oliver Lemke"),
      OUT("complex_refr_index"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("data_f_grid", "data_T_grid"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Frequency grid for refractive index calculation",
               "Temperature grid for refractive index calculation")));

  md_data_raw.push_back(create_mdrecord(
      NAME("complex_refr_indexWaterVisibleNIRHarvey98"),
      DESCRIPTION(R"--(Refractive index of water and steam for the optical and near infrared.

Refractive index as function of temparature, frequency and density.
It is limited only to the real part. The imaginary part is 0.

From:
Revised formulation for the Refractive Index of Water and Steam as a Function
of Wavelength, Temperature and Density
Journal of Physical and Chemical Reference Data 27, 761 (1998), 
https:

See also: http:

Range of validity:

- 271.15K < temperature < 773.15K
- 0 kg m^-3 < density < 1060 kg m^-3
- 157.785504THz < frequency < 1498.96229THz or  0.2µm < wavelength < 1.9µm

Density can be set as Vector of size 1 or it must have the same size as
as data_t_grid.

IMPORTANT: Though the output is *complex_refr_index*, it only contains
the real part. The imaginry part is zero.
)--"),
      AUTHORS("Manfred Brath"),
      OUT("complex_refr_index"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("complex_refr_index"),
      GIN("data_f_grid", "data_t_grid", "density_water", "only_valid_range"),
      GIN_TYPE("Vector", "Vector", "Vector", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "1"),
      GIN_DESC("Frequency grid for refractive index calculation",
               "Temperature grid for refractive index calculation",
               "Density of water",
               "Flag. If true refractive index is calculated only"
               " within range of validity and it will throw an error if outside"
               " range of validity."
               " If false no check is made, so use at your own risk.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Copy"),
      DESCRIPTION(R"--(Copy a workspace variable.

This method can copy any workspace variable
to another workspace variable of the same group. (E.g., a Matrix to
another Matrix.)

As always, output comes first in the argument list!

Usage example:

Copy(f_grid, p_grid)

Will copy the content of ``p_grid`` to *f_grid*. The size of *f_grid*
is adjusted automatically (the normal behaviour for workspace
methods).
)--"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Any"),
      GOUT_DESC("Destination variable."),
      IN(),
      GIN("input"),
      GIN_TYPE("Any"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Source variable."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("covmat1D"),
      DESCRIPTION(R"--(Create 1D covariance matrix.

Creates a 1D covariance matrix for two retrieval quantities on given 
grids from a given functional form. Elements  of the covariance matrix
are computed as::

  S_{i,j} = sigma_i * sigma_j * f(d_{i,j} / l_{i,j})

where d_{i,j} is the distance between the two grid points and l_{i,j}
the mean of the correlation lengths of the grid points.

If a cutoff value ``co`` is given elements with absolute value less than this 
are set to zero.

The following functional forms are available:

- ``\"exp\"``: f(x) = exp(-x) 
- ``\"lin\"``: f(x) = 1.0 - x, for x > 1.0, 0.0 otherwise 
- ``\"gauss\"``: f(x) = exp(-x^2)
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix, Sparse"),
      GOUT_DESC("The matrix in which to store the covariance matrix."),
      IN(),
      GIN("grid_1",
          "grid_2",
          "sigma_1",
          "sigma_2",
          "cls_1",
          "cls_2",
          "co",
          "fname"),
      GIN_TYPE("Vector",
               "Vector",
               "Vector",
               "Vector",
               "Vector",
               "Vector",
               "Numeric",
               "String"),
      GIN_DEFAULT(NODEF, "[]", NODEF, "[]", NODEF, "[]", "0.0", NODEF),
      GIN_DESC("The retrieval grid for the first retrieval quantity.",
               "The retrieval grid for the second retrieval quantity."
               " (If empty taken as grid_1)",
               "The variances of the first retrieval quantity.",
               "The variances of the second retrieval quantity."
               "(If empty taken as sigma_1)",
               "The correlations lengths of the first retrieval quantity.",
               "The correlations lengths of the second retrieval quantity."
               "(If empty taken as cls_1)",
               "The cutoff value for covariance matrix elements.",
               "The name of the functional form to use."),
      PASSWORKSPACE(false),
      SETMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("covmat1DMarkov"),
      DESCRIPTION(R"--(Create Markov Process Covariance Matrix.

Create a markov process covariance matrix for a retrieval quantity on
evenly spaced 1D grid. The correlation between two grid points i,j is 
is computed as::

  cov(i,j) = sigma[i] * sigma[j] * exp(- d(i,j) / lc)

where d(i,j) = abs(grid[i] - grid[j]).

This function also sets covmat_inv_block to the analytically computed inverse
of the covariance matrix of the markov provess, which is tri-diagonal. Note
that this requires the retrieval grid to be evenly spaced.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("output", "out_inverse"),
      GOUT_TYPE("Matrix, Sparse", "Matrix, Sparse"),
      GOUT_DESC(
          "The matrix in which to store the covariance matrix.",
          "The matrix in which to store the inverse of the covariance matrix."),
      IN(),
      GIN("grid", "sigma", "lc", "co"),
      GIN_TYPE("Vector", "Vector", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "0.0"),
      GIN_DESC("The retrieval grid.",
               "The vairance for each grid point.",
               "The correlation length of the Markov process.",
               "The cutoff value below which elements will be set to 0.0"),
      PASSWORKSPACE(false),
      SETMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("covmatDiagonal"),
      DESCRIPTION(R"--(Sets the matrix in covmat_block to a diagonal matrix with the variances
provided in ``vars`` as diagonal elements.

Also sets covmat_block_inv to the inverse of the block so that the
computation of the inverse is avoided.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("output", "out_inverse"),
      GOUT_TYPE("Matrix, Sparse", "Matrix, Sparse"),
      GOUT_DESC(
          "The matrix in which to store the covariance matrix.",
          "The matrix in which to store the inverse of the covariance matrix."),
      IN(),
      GIN("vars"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Variances to be used as diagonal elements of covmat_block."),
      PASSWORKSPACE(false),
      SETMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("covmat_seAddBlock"),
      DESCRIPTION(R"--(Add a block to the measurement covariance matrix *covmat_se*

This functions adds a given dense or sparse matrix as block to the covariance
matrix *covmat_sx*. The position of the block can be given by the generic
arguments ``i`` and ``j``. Note that diagonal blocks must be added in order starting from
in the top left corner. If an off-diagonal block is added it must have corresponding
existing blocks on the diagonal and these must be consistent with the dimensions
of the block.  If ``i`` and ``j``  are not provided, the blok will be added
at the first free spot on the diagonal.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_se"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_se"),
      GIN("block", "i", "j"),
      GIN_TYPE("Matrix, Sparse", "Index", "Index"),
      GIN_DEFAULT(NODEF, "-1", "-1"),
      GIN_DESC("The block to add to the covariance matrix",
               "Index of a retrieval quantity. Must satisfy ``i`` <= ``j``.",
               "Index of a retrieval quantity. Must satisfy ``i`` <= ``j``."),
      PASSWORKSPACE(false),
      SETMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("covmat_seAddInverseBlock"),
      DESCRIPTION(R"--(Add the inverse of a block to covariance matrix *covmat_se*

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
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_se"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_se"),
      GIN("block", "i", "j"),
      GIN_TYPE("Matrix, Sparse", "Index", "Index"),
      GIN_DEFAULT(NODEF, "-1", "-1"),
      GIN_DESC("The inverse block to add to the covariance matrix",
               "Index of a retrieval quantity. Must satisfy ``i`` <= ``j``.",
               "Index of a retrieval quantity. Must satisfy ``i`` <= ``j``."),
      PASSWORKSPACE(false),
      SETMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("covmat_seSet"),
      DESCRIPTION(R"--(Set covmat_se to a given matrix.

This sets the measurement covariance matrix *covmat_se* to
the matrix given by the generic input ``covmat``. The covariance
matrix can be of type CovarianceMatrix, Matrix or Sparse.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_se"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("covmat"),
      GIN_TYPE("CovarianceMatrix, Matrix, Sparse"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The matrix to set as the covariance matrix."),
      PASSWORKSPACE(false),
      SETMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("covmat_sxSet"),
      DESCRIPTION(R"--(Set covmat_sx to a given matrix.

This sets the measurement covariance matrix *covmat_sx* to
the matrix given by the generic input ``covmat``. The covariance
matrix can be of type CovarianceMatrix, Matrix or Sparse.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("covmat"),
      GIN_TYPE("CovarianceMatrix, Matrix, Sparse"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The matrix to set as the covariance matrix."),
      PASSWORKSPACE(false),
      SETMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("covmat_sxAddBlock"),
      DESCRIPTION(R"--(Add a block to the a priori covariance matrix *covmat_sx*

This functions adds a given matrix as a block in the covariance
matrix *covmat_sx*. The position of the block can be given by the generic
arguments ``i`` and ``j``, which should give the index of the retrieval quantity in
*jacobian_quantities*, which is given just by the order the quantities have been
added to the retrieval.

If arguments ``i`` and ``j`` are omitted, the block will be added as diagonal block
for the last added retrieval quantity.

If provided, the index ``i`` must be less than or equal to ``j``. Also the provided
block must be consistent with the corresponding retrieval quantities.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx", "jacobian_quantities"),
      GIN("block", "i", "j"),
      GIN_TYPE("Matrix, Sparse", "Index", "Index"),
      GIN_DEFAULT(NODEF, "-1", "-1"),
      GIN_DESC("The block to add to the covariance matrix",
               "Index of a retrieval quantity. Must satisfy ``i`` <= ``j``.",
               "Index of a retrieval quantity. Must satisfy ``i`` <= ``j``."),
      PASSWORKSPACE(false),
      SETMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("covmat_sxAddInverseBlock"),
      DESCRIPTION(R"--(Add the inverse of a block in covariance matrix *covmat_sx*

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
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx", "jacobian_quantities"),
      GIN("block", "i", "j"),
      GIN_TYPE("Matrix, Sparse", "Index", "Index"),
      GIN_DEFAULT(NODEF, "-1", "-1"),
      GIN_DESC("The inverse block to add to the covariance matrix",
               "Index of a retrieval quantity. Must satisfy ``i`` <= ``j``.",
               "Index of a retrieval quantity. Must satisfy ``i`` <= ``j``."),
      PASSWORKSPACE(false),
      SETMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("covmat_sxExtractSqrtDiagonal"),
      DESCRIPTION(R"--(Extract the square root of the diagonal of the state space covariance matrix.

This function extracts the diagonal of the state space covariance matrix
*covmat_sx* and computes its square root. The resulting vector can then
be used as ``x_norm`` argument for the OEM method to avoid scaling problems.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("x_norm"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("The vector containing the square root of the diagonal elements"
                " of *covmat_sx*"),
      IN("covmat_sx"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC(),
      PASSWORKSPACE(false),
      SETMETHOD(false)));

  md_data_raw.push_back(create_mdrecord(
      NAME("Delete"),
      DESCRIPTION(R"--(Deletes a workspace variable.

The variable is marked as uninitialized and its memory freed.
It is not removed from the workspace though, therefore you
don't need to/can't call Create for this variable again.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("v"),
      GIN_TYPE("Any"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Variable to be deleted."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true),
      PASSWORKSPACE(true),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("diameter_maxFromDiameter_volume_equ"),
      DESCRIPTION(R"--(Calculates maximum and area equivalent diameters from volume
equivalent diameter.

This is primarily a help function for using the T-matrix method
and only a few particle shapes are handled. 

For shapes handled and further comments on the input arguments, see
*scat_data_singleTmatrix*.

Area equivalent diameter is the equivalent sphere diameter
corresponding to the \"maximum axial area\". This is the largest
cross-sectional area of the particle, observed either along the
particle's main axis or in the perpendicular direction. That is,
for a cylinder having diameter d and thickness h, this area is
either (pi*d^2)/4 or (h*d).
)--"),
      AUTHORS("Johan Strandgren", "Patrick Eriksson"),
      OUT(),
      GOUT("diameter_max", "diameter_area_equ"),
      GOUT_TYPE("Numeric", "Numeric"),
      GOUT_DESC(
          "Maximum dimension of the particle.",
          "Maximum axial area equivalent diameter of the particle, see above."),
      IN(),
      GIN("shape", "diameter_volume_equ", "aspect_ratio"),
      GIN_TYPE("String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Particle shape.",
               "Particle equivalent volume diameter.",
               "Particle aspect ratio.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("diameter_volume_equFromDiameter_max"),
      DESCRIPTION(R"--(Converts from maximum to volume equivalent diameter.

This is primarily a help function for using the T-matrix part
and only a few particle shapes are handled. 

For shapes handled and further comments on the input arguments,
see *scat_data_singleTmatrix*.

Also the volume is provided. It is simply sqrt(pi*dveq^3/6).
)--"),
      AUTHORS("Johan Strandgren", "Patrick Eriksson"),
      OUT(),
      GOUT("diameter_volume_equ", "volume"),
      GOUT_TYPE("Numeric", "Numeric"),
      GOUT_DESC("Particle volume equivalent diameter.",
                "Volume of the particle."),
      IN(),
      GIN("shape", "diameter_max", "aspect_ratio"),
      GIN_TYPE("String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Particle shape.",
               "Maximum dimension of the particle.",
               "Particle aspect ratio.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("DiagonalMatrix"),
      DESCRIPTION(R"--(Create a diagonal matrix from a vector.

This creates a dense or sparse diagonal matrix with the elements of the given vector
on the diagonal.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix, Sparse"),
      GOUT_DESC("The diagonal matrix"),
      IN(),
      GIN("v"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The vector containing the diagonal elements.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("dlosDiffOfLos"),
      DESCRIPTION(R"--(Derives the difference betwenn zenith and azimuth angles.

Determines the difference between a set of angles (``other_los``)
and a reference direction (``ref_los``). This method reverses the
addition made by *sensor_losAddLosAndDlos*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("dlos"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Derived differences in line-of-sight."),
      IN(),
      GIN("ref_los", "other_los"),
      GIN_TYPE("Vector", "Matrix"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Reference line-of-sight (a single los).",
               "Other line-of-sights (can be multiple los).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("dlosGauss"),
      DESCRIPTION(R"--(Gives a *dlos* suitable for a circular Gaussian response.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("dlos", "dlos_weight_vector"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("fwhm", "npoints", "include_response_in_weight"),
      GIN_TYPE("Numeric", "Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF, "0"),
      GIN_DESC("The full width at half maximum of the Gaussian response.",
               "Number of dlos-directions.",
               "Set to 1 to include the response values in *dlos_weight_vector*.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("dlosUniform"),
      DESCRIPTION(R"--(Gives *dlos* a rectangular coverage, with uniform spacing.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("dlos", "dlos_weight_vector"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("width", "npoints", "crop_circular"),
      GIN_TYPE("Numeric", "Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF, "0"),
      GIN_DESC("The full width, in each dimension, in degrees.",
               "Number of points over the width, in each dimension (>1).",
               "Set to 1, to crop dlos-es to obtain a pseudo-circular pattern.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_fieldDisort"),
      DESCRIPTION(R"--(Interface to the DISORT scattering solver (by Stamnes et al.).

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

- ``\"Layer optical thickness\"``: Matrix [f_grid, size of p_grid - 1] layer optical thickness.
- ``\"Single scattering albedo\"``: Matrix [f_grid, size of p_grid - 1] layer single\" scattering albedo.
- ``\"Direct beam\"``: Matrix [f_grid, p_grid]. Attenuated direct at level. Zero, if no sun is present
)--"),
      AUTHORS("Claudia Emde, Jana Mendrok", "Manfred Brath"),
      OUT("cloudbox_field","disort_aux"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmfields_checked",
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
         "disort_aux_vars"),
      GIN("nstreams", "Npfct", "only_tro", "quiet", "emission","intensity_correction"),
      GIN_TYPE("Index", "Index", "Index", "Index", "Index", "Index"),
      GIN_DEFAULT("8", "181", "0", "0", "1", "1"),
      GIN_DESC("Number of polar angle directions (streams) in DISORT"
               " solution (must be an even number).",
               "Number of angular grid points to calculate bulk phase"
               " function on (and derive Legendre polynomials from). If <0,"
               " the finest za_grid from scat_data will be used.",
               "Set to 1 if the scattering data is just of TRO type. Has"
               " effect only if Npfct > 3 or Npfct<0, but then leads to"
               " much faster calculations.",
               "Silence C Disort warnings.",
               "Enables blackbody emission. Set to zero, if no"
               " Emission e. g. like in visible regime for earth"
               " is needed",
               "Enables intensity correction. Importantant for low number of"
               " streams. Set to zero, if problems encounter or using a high number"
               " of streams (>30)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("cloudbox_fieldDisortWithARTSSurface"),
      DESCRIPTION(R"--(Interface to the DISORT scattering solver (by Stamnes et al.).

As *cloudbox_fieldDisort* but uses *surface_rtprop_agenda*.

The Lambertian surface reflection is set by *surface_rtprop_agenda*.
If the GIN inc_angle is inside of the range [0,90], the reflection is
set according to the result of *surface_rtprop_agenda* for this incidence
angle. Otherwise (default) is to call *surface_rtprop_agenda* for
multiple angles, to estimate the hemispheric mean value.

Some auxiliary quantities can be obtained. Auxiliary
quantities are selected by *disort_aux_vars* and returned by *disort_aux*.
Valid choices for auxiliary data are:

- ``\"Layer optical thickness\"``: Matrix [f_grid, size of p_grid - 1] layer optical thickness.
- ``\"Single scattering albedo\"``: Matrix [f_grid, size of p_grid - 1] layer single\"scattering albedo.
- ``\"Direct beam\"``: Matrix [f_grid, p_grid]. Attenuated direct at level.Zero, if no sun is present
)--"),
      AUTHORS("Claudia Emde, Jana Mendrok", "Manfred Brath"),
      OUT("cloudbox_field","disort_aux"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmfields_checked",
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
         "disort_aux_vars"),
      GIN("nstreams", "Npfct", "only_tro", "quiet", "emission", "intensity_correction", "inc_angle"),
      GIN_TYPE("Index", "Index", "Index", "Index", "Index", "Index","Numeric"),
      GIN_DEFAULT("8", "181", "0", "0", "1", "1", "-1"),
      GIN_DESC("Number of polar angle directions (streams) in DISORT "
               " solution (must be an even number).",
               "Number of angular grid points to calculate bulk phase"
               " function on (and derive Legendre polynomials from). If <0,"
               " the finest za_grid from scat_data will be used.",
               "Set to 1 if the scattering data is just of TRO type. Has"
               " effect only if Npfct > 3 or Npfct<0, but then leads to"
               " much faster calculations.",
               "Silence C Disort warnings.",
               "Enables blackbody emission. Set to zero, if no"
               " Emission e. g. like in visible regime for earth"
               " is needed",
               "Enables intensity correction. Importantant for low number of"
               " streams. Set to zero, if problems encounter or using a high number"
               " of streams (>30)",
               "Incidence angle, see above.")));


  md_data_raw.push_back(create_mdrecord(
      NAME("spectral_radiance_fieldDisortClearsky"),
      DESCRIPTION(R"--(Interface to the DISORT scattering solver (by Stamnes et al.).
for running clear-sky cases.

The method runs DISORT with *pnd_field* set to zero.

Note that this version returns *spectral_radiance_field*, i.e.
the solution for the full atmosphere. The standard *cloudbox_fieldDisort*
only returns the field inside the cloudbox.

Some auxiliary quantities can be obtained. Auxiliary
quantities are selected by *disort_aux_vars* and returned by *disort_aux*.
Valid choices for auxiliary data are:

- ``\"Layer optical thickness\"``: Matrix [f_grid, size of p_grid - 1] layer optical thickness.
- ``\"Single scattering albedo\"``: Matrix [f_grid, size of p_grid - 1] layer single scattering albedo.
- ``\"Direct beam\"``: Matrix [f_grid, p_grid]. Level direct spectral radiance. Zero, if no sun is present
)--"),
      AUTHORS("Patrick Eriksson", "Manfred Brath"),
      OUT("spectral_radiance_field","disort_aux"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmfields_checked",
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
         "disort_aux_vars"),
      GIN("nstreams", "quiet", "emission", "intensity_correction"),
      GIN_TYPE("Index", "Index", "Index", "Index"),
      GIN_DEFAULT("8", "0", "1", "1"),
      GIN_DESC("Number of polar angle directions (streams) in DISORT"
               " solution (must be an even number).",
               "Silence C Disort warnings.",
               "Enables blackbody emission. Set to zero, if no "
               " Emission e. g. like in visible regime for earth"
               " is needed",
               "Enables intensity correction. Importantant for low number of "
               " streams. Set to zero, if problems encounter or using a high number "
               " of streams (>30)")));


  md_data_raw.push_back(create_mdrecord(
      NAME("spectral_irradiance_fieldDisort"),
      DESCRIPTION(R"--(Interface to the DISORT scattering solver (by Stamnes et al.).
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

- ``\"Layer optical thickness\"``: Matrix [f_grid, size of p_grid - 1] layer optical thickness.
- ``\"Single scattering albedo\"``: Matrix [f_grid, size of p_grid - 1] layer single scattering albedo.
- ``\"Direct downward spectral irradiance\"``: Matrix [f_grid, p_grid]. Direct downward spectral irradiance. Zero, if no sun is present. 
- ``\"dFdtau\"``: Matrix [f_grid, p_grid]. Flux divergence in optical thickness space.
)--"),
      AUTHORS("Manfred Brath"),
      OUT("spectral_irradiance_field","disort_aux"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmfields_checked",
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
         "disort_aux_vars"),
      GIN("nstreams", "Npfct", "only_tro", "quiet", "emission","intensity_correction"),
      GIN_TYPE("Index", "Index", "Index", "Index", "Index", "Index"),
      GIN_DEFAULT("6", "181", "0", "0", "1", "1"),
      GIN_DESC("Number of polar angle directions (streams) in DISORT"
               " solution (must be an even number).",
               "Number of angular grid points to calculate bulk phase"
               " function on (and derive Legendre polynomials from). If <0,"
               " the finest za_grid from scat_data will be used.",
               "Set to 1 if the scattering data is just of TRO type. Has"
               " effect only if Npfct > 3 or Npfct<0, but then leads to"
               " much faster calculations.",
               "Silence C Disort warnings.",
               "Enables blackbody emission. Set to zero, if no "
               " Emission e. g. like in visible regime for earth"
               " is needed",
               "Enables intensity correction. Importantant for low number of "
               " streams. Set to zero, if problems encounter or using a high number "
               " of streams (>30)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("DOBatchCalc"),
      DESCRIPTION(R"--(Performs batch calculations for radiation fields.

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
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("dobatch_cloudbox_field",
          "dobatch_radiance_field",
          "dobatch_irradiance_field",
          "dobatch_spectral_irradiance_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ybatch_start", "ybatch_n", "dobatch_calc_agenda"),
      GIN("robust"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("A flag with value 1 or 0. If set to one, the batch "
               "calculation will continue, even if individual jobs fail. In "
               "that case, a warning message is written to screen and file "
               "(out1 output stream), and the output array entry for the "
               "failed job in the output fields is left empty.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("DOAngularGridsSet"),
      DESCRIPTION(R"--(Sets the angular grids for Discrete Ordinate type scattering
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
)--"),
      AUTHORS("Claudia Emde"),
      OUT("doit_za_grid_size", "aa_grid", "za_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("N_za_grid", "N_aa_grid", "za_grid_opt_file"),
      GIN_TYPE("Index", "Index", "String"),
      GIN_DEFAULT(NODEF, "1", ""),
      GIN_DESC("Number of grid points in zenith angle grid. "
               "Recommended value is 19.",
               "Number of grid points in azimuth angle grid. "
               "Recommended value is 37.",
               "Name of special grid for RT part.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("DoitCalc"),
      DESCRIPTION(R"--(Main DOIT method.

This method executes *doit_mono_agenda* for each frequency
in *f_grid*. The output is the radiation field inside the cloudbox
(*cloudbox_field*).
)--"),
      AUTHORS("Claudia Emde"),
      OUT("cloudbox_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field",
         "atmfields_checked",
         "atmgeom_checked",
         "cloudbox_checked",
         "scat_data_checked",
         "cloudbox_on",
         "f_grid",
         "doit_mono_agenda",
         "doit_is_initialized"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("DoitGetIncoming"),
      DESCRIPTION(R"--(Calculates incoming radiation field of the cloudbox by repeated
radiative transfer calculations.

The method performs monochromatic pencil beam calculations for
all grid positions on the cloudbox boundary, and all directions
given by scattering angle grids (*aa_grid*). Found radiances
are stored in *cloudbox_field* which can be used as boundary
conditions when scattering inside the cloud box is solved by the
*DoitCalc* method.

Note that *cloudbox_field* will always hold intensity in terms of
radiances, regardless of the setting of *iy_unit* (unit conversion
is done within *yCalc* or *iyCalc*, which will provide their output
in terms of the specified *iy_unit*; no explicit unit conversion by
the user necessary.).
)--"),
      AUTHORS("Sreerekha T.R.", "Claudia Emde"),
      OUT("cloudbox_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field",
         "atmfields_checked",
         "atmgeom_checked",
         "cloudbox_checked",
         "doit_is_initialized",
         "iy_main_agenda",
         "atm_field",
         "cloudbox_on",
         "cloudbox_limits",
         "f_grid",
         
         "za_grid",
         "aa_grid"),
      GIN("rigorous", "maxratio"),
      GIN_TYPE("Index", "Numeric"),
      GIN_DEFAULT("1", "100"),
      GIN_DESC(
          "Fail if incoming field is not safely interpolable.",
          "Maximum allowed ratio of two radiances regarded as interpolable.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("DoitGetIncoming1DAtm"),
      DESCRIPTION(R"--(As *DoitGetIncoming* but assumes clear sky part to be 1D.

The incoming field is calculated only for one position and azimuth
angle for each cloud box boundary, and obtained values are used
for all other postions and azimuth angles. This works if a 3D
cloud box is put into an 1D background atmosphere.

This method can only be used for 3D cases.
)--"),
      AUTHORS("Sreerekha T.R.", "Claudia Emde"),
      OUT("cloudbox_field", "cloudbox_on"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field",
         "atmfields_checked",
         "atmgeom_checked",
         "cloudbox_checked",
         "doit_is_initialized",
         "iy_main_agenda",
         "atm_field",
         "cloudbox_on",
         "cloudbox_limits",
         "f_grid",
         
         "za_grid",
         "aa_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("DoitInit"),
      DESCRIPTION(R"--(Initialises variables for DOIT scattering calculations.

Note that multi-dimensional output variables (Tensors, specifically)
are NaN-initialized. That is, this methods needs to be called
BEFORE other WSMs that provide input to *DoitCalc*, e.g. before
*DoitGetIncoming*.
)--"),
      AUTHORS("Claudia Emde"),
      OUT("doit_scat_field", "cloudbox_field", "doit_is_initialized"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         "f_grid",
         "za_grid",
         "aa_grid",
         "doit_za_grid_size",
         "cloudbox_on",
         "cloudbox_limits"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("DoitScatteringDataPrepare"),
      DESCRIPTION(R"--(Prepares single scattering data for a DOIT scattering calculation.

First the scattering data is interpolated in frequency using
*scat_data_monoCalc*. Then the phase matrix data is
transformed or interpolated from the raw data to the laboratory frame
for all possible combinations of the angles contained in the angular
grids which are set in *DOAngularGridsSet*. The resulting phase
matrices are stored in *pha_mat_sptDOITOpt*.
)--"),
      AUTHORS("Claudia Emde"),
      OUT("pha_mat_sptDOITOpt", "scat_data_mono", "pha_mat_doit", "aa_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_za_grid_size",
         "aa_grid",
         "scat_data",
         "scat_data_checked",
         "f_index",
         
         "atm_field",
         "cloudbox_limits",
         "pnd_field",
         "pha_mat_spt_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("DoitWriteIterationFields"),
      DESCRIPTION(R"--(Writes DOIT iteration fields.

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
)--"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_iteration_counter", "cloudbox_field_mono", "f_index"),
      GIN("iterations", "frequencies"),
      GIN_TYPE("ArrayOfIndex", "ArrayOfIndex"),
      GIN_DEFAULT("[-1]", "[-1]"),
      GIN_DESC("Selection of iterations to store.",
               "Selection of frequencies to store.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("doit_conv_flagAbs"),
      DESCRIPTION(R"--(DOIT convergence test (maximum absolute difference).

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
)--"),
      AUTHORS("Claudia Emde"),
      OUT("doit_conv_flag", "doit_iteration_counter", "cloudbox_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_conv_flag",
         "doit_iteration_counter",
         "cloudbox_field_mono",
         "cloudbox_field_mono_old"),
      GIN("epsilon", "max_iterations", "nonconv_return_nan"),
      GIN_TYPE("Vector", "Index", "Index"),
      GIN_DEFAULT(NODEF, "100", "0"),
      GIN_DESC("Limits for convergence. A vector with length matching "
               "``stokes_dim`` with unit [W / (m^2 Hz sr)].",
               "Maximum number of iterations allowed to reach convergence"
               "limit.",
               "Flag whether to accept result at max_iterations (0=default)"
               "or whether to return NaNs in case of non-convergence at"
               "max_iterations")));

  md_data_raw.push_back(create_mdrecord(
      NAME("doit_conv_flagAbsBT"),
      DESCRIPTION(R"--(DOIT convergence test (maximum absolute difference in Rayleigh Jeans 
BT)

As *doit_conv_flagAbs* but convergence limits are specified in
Rayleigh-Jeans brighntess temperatures.
)--"),
      AUTHORS("Sreerekha T.R.", "Claudia Emde"),
      OUT("doit_conv_flag", "doit_iteration_counter", "cloudbox_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_conv_flag",
         "doit_iteration_counter",
         "cloudbox_field_mono",
         "cloudbox_field_mono_old",
         "f_grid",
         "f_index"),
      GIN("epsilon", "max_iterations", "nonconv_return_nan"),
      GIN_TYPE("Vector", "Index", "Index"),
      GIN_DEFAULT(NODEF, "100", "0"),
      GIN_DESC("Limits for convergence. A vector with length matching "
               "``stokes_dim`` with unit [K].",
               "Maximum number of iterations allowed to reach convergence"
               "limit.",
               "Flag whether to accept result at max_iterations (0=default)"
               "or whether to return NaNs in case of non-convergence at"
               "max_iterations")));

  md_data_raw.push_back(create_mdrecord(
      NAME("doit_conv_flagLsq"),
      DESCRIPTION(R"--(DOIT convergence test (least squares).

As *doit_conv_flagAbsBT* but applies a least squares convergence
test between two successive iteration fields.

Warning: This method is not recommended because this kind of
convergence test is not sufficiently strict, so that the
DOIT result might be wrong.
)--"),
      AUTHORS("Claudia Emde"),
      OUT("doit_conv_flag", "doit_iteration_counter", "cloudbox_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_conv_flag",
         "doit_iteration_counter",
         "cloudbox_field_mono",
         "cloudbox_field_mono_old",
         "f_grid",
         "f_index"),
      GIN("epsilon", "max_iterations", "nonconv_return_nan"),
      GIN_TYPE("Vector", "Index", "Index"),
      GIN_DEFAULT(NODEF, "100", "0"),
      GIN_DESC("Limits for convergence. A vector with length matching "
               "``stokes_dim`` with unit [K].",
               "Maximum number of iterations allowed to reach convergence"
               "limit.",
               "Flag whether to accept result at max_iterations (0=default)"
               "or whether to return NaNs in case of non-convergence at"
               "max_iterations")));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("OptimizeDoitPressureGrid"),
      DESCRIPTION(R"--(Optimization of the pressure grid for RT calculation.

The methods consists of three parts:

1. Calculate the single scattering albedo and the scattering optical
   thickness from the scattering and absorption species. 
2. Enhance z_field according to the two thresholds sgl_alb_max and tau_scat_max.
   If the resulting cloudbox size is bigger than cloudbox_size_max, this step is 
   repeated with a higher threshold of tau_scat_max. 
3. Interpolate all variables used in doit_mono_agenda to the new z_field 
   This method should be called inside
   *doit_mono_agenda*, right before ``cloudbox_field_monoIterate``. It can 
   only be used if *ScatSpeciesMerge* has been called and if it is
   called, *cloudbox_field_monoOptimizeReverse* has to be
   called right after ``cloudbox_field_monoIterate`` to interpolate
   *cloudbox_field_mono* back to the original size.

Optimization currently only works with ``stokes_dim`` = 1 .
)--"),
      AUTHORS("Jakob Doerr"),
      OUT("p_grid",
          "pnd_field",
          "t_field",
          "scat_data_mono",
          "z_field",
          "cloudbox_limits",
          "cloudbox_field_mono",
          "pha_mat_doit",
          "atm_field",
          "p_grid_orig"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("p_grid",
         "pnd_field",
         "t_field",
         "scat_data_mono",
         "z_field",
         "cloudbox_limits",
         "cloudbox_field_mono",
         "pha_mat_doit",
         "atm_field",
         "f_grid",
         "f_index",
         "propmat_clearsky_agenda"),
      GIN("tau_scat_max", "sgl_alb_max", "cloudbox_size_max"),
      GIN_TYPE("Numeric", "Numeric", "Index"),
      GIN_DEFAULT("0.1", "0.9", "200"),
      GIN_DESC("Maximum scattering optical thickness",
               "Maximum single scattering albedo",
               "Maximum cloudbox size")));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("doit_scat_fieldCalc"),
      DESCRIPTION(R"--(Calculates the scattering integral field in the DOIT module.

The scattering integral field is generated by integrating
the product of phase matrix and Stokes vector over all incident
angles. For more information please refer to AUG.
)--"),
      AUTHORS("Sreerekha T.R.", "Claudia Emde"),
      OUT("doit_scat_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_scat_field",
         "pha_mat_spt_agenda",
         "cloudbox_field_mono",
         "pnd_field",
         "atm_field",
         "cloudbox_limits",
         "za_grid",
         "aa_grid",
         "doit_za_grid_size",
         "pha_mat_doit"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("doit_scat_fieldCalcLimb"),
      DESCRIPTION(R"--(Calculates the scattering integral field in the DOIT module (limb).

The scattering integral field is the field generated by integrating
the product of phase matrix and the Stokes vector over all incident
angles.

For limb simulations it makes sense to use different
zenith angle grids for the scattering integral part and the RT part,
because the latter part requires a much finer resolution near
90 degrees. Taking an optimized grid for the RT part and an equidistant
grid for the scattering integral part saves very much CPU time.
This method uses the equidistant za_grid defined in
*DOAngularGridsSet* and it should always be used for limb
simulations.

For more information please refer to AUG.
)--"),
      AUTHORS("Claudia Emde"),
      OUT("doit_scat_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_scat_field",
         "pha_mat_spt_agenda",
         "cloudbox_field_mono",
         "pnd_field",
         "atm_field",
         "cloudbox_limits",
         "za_grid",
         "aa_grid",
         "doit_za_grid_size",
         "doit_za_interp",
         "pha_mat_doit"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("doit_za_grid_optCalc"),
      DESCRIPTION(R"--(Zenith angle grid optimization for scattering calculation.

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
)--"),
      AUTHORS("Claudia Emde"),
      OUT("doit_za_grid_opt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field_mono", "za_grid", "doit_za_interp"),
      GIN("acc"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Accuracy to achieve [%].")));

  md_data_raw.push_back(create_mdrecord(
      NAME("doit_za_interpSet"),
      DESCRIPTION(R"--(Define interpolation method for zenith angle dimension.

You can use this method to choose the interpolation method for
interpolations in the zenith angle dimension.
)--"),
      AUTHORS("Claudia Emde"),
      OUT("doit_za_interp"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("interp_method"),
      GIN_TYPE("String"),
      GIN_DEFAULT("linear"),
      GIN_DESC("Interpolation method (\"linear\" or \"polynomial\").")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Duration"),
      DESCRIPTION(R"--(Sets the seconds between two times.
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("duration"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Time in seconds between ``start`` and ``end``"),
      IN(),
      GIN("start", "end"),
      GIN_TYPE("Time", "Time"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Start time", "End time")));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("ecs_dataAddMakarov2020"),
      DESCRIPTION(R"--(Sets the O2-66 microwave band data for ECS.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ecs_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ecs_data", "isotopologue_ratios"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("ecs_dataAddRodrigues1997"),
      DESCRIPTION(R"--(Sets the CO2-626, CO2-636, and CO2-628 IR band data for ECS.

Note that the broadening species has to be N2 and not AIR for the band,
and that N2 VMR must be present
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ecs_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ecs_data", "isotopologue_ratios"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("ecs_dataAddTran2011"),
      DESCRIPTION(R"--(Sets the CO2-626, CO2-636, and CO2-628 IR band data for ECS.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ecs_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ecs_data", "isotopologue_ratios"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("ecs_dataAddTran2006"),
      DESCRIPTION(R"--(Sets the O2-66 visible band data for ECS.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ecs_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ecs_data", "isotopologue_ratios"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("ecs_dataInit"),
      DESCRIPTION(R"--(Resets/initializes the ECS data.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ecs_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("ecs_dataAddMeanAir"),
      DESCRIPTION(R"--(Sets ECS data for air from other data if available.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ecs_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ecs_data"),
      GIN("vmrs", "specs"),
      GIN_TYPE("Vector", "ArrayOfSpeciesTag"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC(
        "VMRs of air species",
        "Air species")));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("ecs_dataAddSpeciesData"),
      DESCRIPTION(R"--(Sets ECS data for one set of species and quantum identifiers.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ecs_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ecs_data", "isotopologue_ratios"),
      GIN("qid", "species", "scaling_type", "scaling", "beta_type", "beta", "lambda_type", "lambda", "collisional_distance_type", "collisional_distance"),
      GIN_TYPE("QuantumIdentifier", "String", "String", "Vector", "String", "Vector", "String", "Vector", "String", "Vector"),
      GIN_DEFAULT(NODEF, NODEF, "T0", NODEF, "T0", NODEF, "T0", NODEF, "T0", NODEF),
      GIN_DESC(
        "Band identifier",
        "Species identifier",
        "Temperature model for the main scaling coefficients for Q",
        "Main scaling coefficients for Q",
        "Temperature model for the energy scaling coefficient for Q",
        "Energy scaling coefficient for Q",
        "Temperature model for the energy exponent for Q",
        "Energy exponent for Q",
        "Temperature model for the mean collision interaction distance",
        "Mean collision interaction distance")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Error"),
      DESCRIPTION(R"--(Issues an error and exits ARTS.

This method can be placed in agendas that must be specified, but
are expected not to be used for the particular case. An inclusion
in *surface_rtprop_agenda* could look like::

  Error{\"Surface interceptions of propagation path not expected.\"}

Ignore and other dummy method calls must still be included.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("msg"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("String describing the error.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Exit"),
      DESCRIPTION(R"--(Stops the execution and exits ARTS.

This method is handy if you want to debug one of your control
files. You can insert it anywhere in the control file. When
it is reached, it will terminate the program.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("Extract"),
      DESCRIPTION(R"--(Extracts an element from an array.

Copies the element with the given Index from the input
variable to the output variable.

For a Tensor3 as an input, it copies the page with the given
Index from the input Tensor3 variable to the output Matrix.

In other words, the selection is always done on the first dimension.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("needle"),
      GOUT_TYPE("Index, ArrayOfIndex, Numeric, Vector,"
                "Matrix, Matrix,"
                "Tensor3, Tensor4, Tensor4,"
		"GriddedField2,"
                "GriddedField3, ArrayOfGriddedField3,"
                "GriddedField4, String,"
                "SingleScatteringData, ArrayOfSingleScatteringData,"
                "TelsemAtlas,"
                "QuantumIdentifier"),
      GOUT_DESC("Extracted element."),
      IN(),
      GIN("haystack", "index"),
      GIN_TYPE("ArrayOfIndex, ArrayOfArrayOfIndex, Vector, ArrayOfVector,"
               "ArrayOfMatrix, Tensor3,"
               "Tensor4, ArrayOfTensor4, Tensor5,"
	       "ArrayOfGriddedField2,"
               "ArrayOfGriddedField3, ArrayOfArrayOfGriddedField3,"
               "ArrayOfGriddedField4, ArrayOfString,"
               "ArrayOfSingleScatteringData,"
               "ArrayOfArrayOfSingleScatteringData,"
               "ArrayOfTelsemAtlas,"
               "ArrayOfQuantumIdentifier",
               "Index"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Variable to extract from.",
               "Position of the element which should be extracted."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("ExtractFromMetaSingleScatSpecies"),
      DESCRIPTION(R"--(Extract (numeric) parameters from scat_meta of a single scattering
species.

...
)--"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("meta_param"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("The extracted meta parameter values."),
      IN("scat_meta"),
      GIN("meta_name", "scat_species_index"),
      GIN_TYPE("String", "Index"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Name of the meta parameter to extract.",
               "Array index of scattering species from which to extract.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ext_matAddGas"),
      DESCRIPTION(R"--(Add gas absorption to all diagonal elements of extinction matrix.

The task of this method is to sum up the gas absorption of the
different gas species and add the result to the extinction matrix.
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("ext_mat"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ext_mat", "propmat_clearsky"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("FastemStandAlone"),
      DESCRIPTION(R"--(Stand-alone usage of FASTEM.

FASTEM is a parameterisation of the emissivity of water surfaces
including the impact of waves, salinity and non-specular effects.
This is more or less direct interface to FASTEM, but slightly
adopted to fit with ARTS. The unit of frequency and salinity
differ, and this version is \"vectorised\" in frequency.

The output is four emissivity and reflectivity values for each
frequency. These values are defined in Eq. 13 of  \"An Improved
Fast Microwave Water Emissivity Model\" by Liu, Weng and English,
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
)--"),
      AUTHORS("Oliver Lemke, Patrick Eriksson"),
      OUT(),
      GOUT("emissivity", "reflectivity"),
      GOUT_TYPE("Matrix", "Matrix"),
      GOUT_DESC("Emission values. One row for each frequency. See above.",
                "Reflectivity values. One row for each frequency. See above."),
      IN("f_grid", "surface_skin_t"),
      GIN("za",
          "salinity",
          "wind_speed",
          "rel_aa",
          "transmittance",
          "fastem_version"),
      GIN_TYPE("Numeric", "Numeric", "Numeric", "Numeric", "Vector", "Index"),
      GIN_DEFAULT(NODEF, "0.035", NODEF, NODEF, NODEF, "6"),
      GIN_DESC("Zenith angle of line-of-sigh, 90 to 180 deg.",
               "Salinity, 0-1. That is, 3% is given as 0.03.",
               "Wind speed.",
               "Azimuth angle between wind direction and line-of-sight. "
               "This angle is measured clockwise from north, i.e. E=90deg.",
               "The transmittance of the atmosphere, along the propagation "
               "path of the downwelling radiation. One value per frequency.",
               "The version of FASTEM to use.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("FlagOff"),
      DESCRIPTION(R"--(Sets an index variable that acts as an on/off flag to 0.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("flag"),
      GOUT_TYPE("Index"),
      GOUT_DESC("Variable to set to 0."),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("FlagOn"),
      DESCRIPTION(R"--(Sets an index variable that acts as an on/off flag to 1.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("flag"),
      GOUT_TYPE("Index"),
      GOUT_DESC("Variable to set to 1."),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("Flatten"),
      DESCRIPTION(R"--(Flattens an ArrayOfArray<T> to Array<T> or an Array
of matpack-types to a larger dimension matpack (if dimensions agree)

The intended transformation for arrays is (sub-arrays can have different sizes):
    {{a, b, c}, {d, e}} -> {a, b, c, d, e}

The intended transformation for arrays to matpack types is (sub-types must have same size):
    {{a, b, c}, {d, e, f}} -> {a, b, c, d, e, f}
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("ArrayOfTime,ArrayOfVector,Matrix,Tensor3,Tensor4,Tensor5,Tensor6,Tensor7"),
      GOUT_DESC("Flatter array/matpack-type"),
      IN(),
      GIN("input"),
      GIN_TYPE("ArrayOfArrayOfTime,ArrayOfArrayOfVector,ArrayOfVector,ArrayOfMatrix,ArrayOfTensor3,ArrayOfTensor4,ArrayOfTensor5,ArrayOfTensor6"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("An array")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ForLoop"),
      DESCRIPTION(R"--(A simple for-loop.

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
)--"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("forloop_agenda"),
      GIN("start", "stop", "step"),
      GIN_TYPE("Index", "Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Start value.", "End value.", "Step size.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("FrequencyFromWavelength"),
      DESCRIPTION(R"--(Convert from wavelength [m] to frequency [Hz].

This is a generic method. It can take a single wavelength value or a wavelength vector as input.
)--"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT("frequency"),
      GOUT_TYPE("Numeric, Vector"),
      GOUT_DESC("frequency [Hz]"),
      IN(),
      GIN("wavelength"),
      GIN_TYPE("Numeric, Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("wavelength [m]")));

  md_data_raw.push_back(create_mdrecord(
      NAME("FrequencyFromCGSAngularWavenumber"),
      DESCRIPTION(R"--(Convert from angular wavenumber [cm^-1] to frequency [Hz].

This converts angular wavenumber (2*PI/wavelength) into frequency.
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("frequency"),
      GOUT_TYPE("Numeric, Vector"),
      GOUT_DESC("frequency [Hz]"),
      IN(),
      GIN("angular_wavenumber"),
      GIN_TYPE("Numeric, Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("angular wavenumber [cm^-1]")));

  md_data_raw.push_back(create_mdrecord(
      NAME("FrequencyFromCGSKayserWavenumber"),
      DESCRIPTION(R"--(Convert from Kayser wavenumber [cm^-1] to frequency [Hz].

This converts Kayser wavenumber (1/wavelength) into frequency.
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("frequency"),
      GOUT_TYPE("Numeric, Vector"),
      GOUT_DESC("frequency [Hz]"),
      IN(),
      GIN("kayser_wavenumber"),
      GIN_TYPE("Numeric, Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Kayser wavenumber [cm^-1]")));

  md_data_raw.push_back(create_mdrecord(
      NAME("f_gridFromAbsorptionLines"),
      DESCRIPTION(R"--(Sets *f_grid* to a grid relative to *abs_lines_per_species*

Each line will have *abs_lines_per_species* will have a grid
of ``num_freqs`` grid points in [ f0 + ``delta_f_low``, f0 + ``delta_f_upp`` ],
where f0 is the line center.

Before leaving the function, *f_grid* is sorted.

Note that this method could generate significantly large *f_grid*
if used carelessly
)--"),
      AUTHORS("Richard Larsson"),
      OUT("f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN("delta_f_low", "delta_f_upp", "num_freqs"),
      GIN_TYPE("Numeric", "Numeric", "Index"),
      GIN_DEFAULT("-5e6", "5e6", NODEF),
      GIN_DESC("Lower range of delta f",
               "Upper range of delta f",
               "Number of frequencies")));

  md_data_raw.push_back(create_mdrecord(
      NAME("f_gridFromGasAbsLookup"),
      DESCRIPTION(R"--(Sets *f_grid* to the frequency grid of *abs_lookup*.

Must be called between importing/creating raw absorption table and
call of *abs_lookupAdapt*.
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lookup"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("f_gridFromSensorAMSU"),
      DESCRIPTION(R"--(Automatically calculate f_grid to match the sensor.

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
)--"),
      AUTHORS("Stefan Buehler, Mathias Milz"),
      OUT("f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("lo_multi", "f_backend_multi", "backend_channel_response_multi"),
      GIN("spacing"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(".1e9"),
      GIN_DESC("Desired grid spacing in Hz.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("f_gridFromSensorAMSUgeneric"),
      DESCRIPTION(R"--(Automatcially calculate f_grid to match the sensor. 
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
)--"),
      AUTHORS("Oscar Isoz"),
      OUT("f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_backend_multi", "backend_channel_response_multi"),
      GIN("spacing", "verbosityVect"),
      GIN_TYPE("Numeric", "Vector"),
      GIN_DEFAULT(".1e9", "[]"),
      GIN_DESC("Desired grid spacing in Hz.", "Bandwidth adjusted spacing")));

  md_data_raw.push_back(create_mdrecord(
      NAME("f_gridFromSensorHIRS"),
      DESCRIPTION(R"--(Automatically calculate f_grid to match the sensor.

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
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_backend", "backend_channel_response"),
      GIN("spacing"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("5e8"),
      GIN_DESC("Desired grid spacing in Hz.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("f_gridMetMM"),
      DESCRIPTION(R"--(Sets *f_grid* and associated variables match MetMM settings.

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
)--"),
      AUTHORS("Oliver Lemke", "Patrick Eriksson"),
      OUT("f_grid",
          "f_backend",
          "channel2fgrid_indexes",
          "channel2fgrid_weights"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("met_mm_backend"),
      GIN("freq_spacing", "freq_number", "freq_merge_threshold"),
      GIN_TYPE("Vector", "ArrayOfIndex", "Numeric"),
      GIN_DEFAULT("[.1e9]", "[-1]", "1"),
      GIN_DESC("Desired grid spacing in Hz.",
               "Number of frequencies per passband for each channel.",
               "Merge frequencies that are closer than this value in Hz.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("yMaskOutsideMedianRange"),
      DESCRIPTION(R"--(Masks values not within the range as NaN::

  [median(y) - dx, median(y) + dx]

Ignores NaNs in median calculations.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("y"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("y"),
      GIN("dx"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Range plus-minus the median of unmasked values")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ybatchMaskOutsideMedianRange"),
      DESCRIPTION(R"--(Apply *yMaskOutsideMedianRange* for each *y* in *ybatch*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ybatch"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ybatch"),
      GIN("dx"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Range plus-minus the median of unmasked values")));

  md_data_raw.push_back(create_mdrecord(
      NAME("yDoublingMeanFocus"),
      DESCRIPTION(R"--(Focus in on *y* around some *f_grid*, then sets *f_grid* to
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
)--"),
      AUTHORS("Richard Larsson"),
      OUT("f_grid", "y"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid", "y"),
      GIN("f0", "df"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT("-1", "-1"),
      GIN_DESC("User input for F0 [see description for default]",
               "User input for DF [see description for default]")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ybatchDoublingMeanFocus"),
      DESCRIPTION(R"--(See *yDoublingMeanFocus*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("f_grid", "ybatch"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid", "ybatch"),
      GIN("f0", "df"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT("-1", "-1"),
      GIN_DESC("User input for F0 [see description for default]",
               "User input for DF [see description for default]")));

  md_data_raw.push_back(create_mdrecord(
      NAME("gas_scatteringOff"),
      DESCRIPTION(R"--(Deactivates the gas_scattering within radiative transfer calculations.
)--"),
      AUTHORS("Manfred Brath"),
      OUT("gas_scattering_do",
          "gas_scattering_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC(),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("gas_scattering_coefAirSimple"),
      DESCRIPTION(R"--(Calculates of scattering coefficient matrix for air.

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
)--"),
      AUTHORS("Jon Petersen"),
      OUT("gas_scattering_coef"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         "rtp_pressure",
         "rtp_temperature"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("gas_scattering_coefXsecConst"),
      DESCRIPTION(R"--(Calculates the spectrum of scattering coefficient matrices.

It calculates the spectrum of scattering coefficient matrices from 
constant spectrum of scattering cross section matrices, atmospheric pressure,
temperature for one point in the atmosphere. Basically, it multiplies
the cross sections with the number density of gas molecules under the
assumption of an ideal gas. The result is returned in *gas_scattering_coef*. The
atmospheric  pressure  and  temperature  state  has  to  be  specified
by  *rtp_pressure*, *rtp_temperature*.
)--"),
      AUTHORS("Manfred Brath"),
      OUT("gas_scattering_coef"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         "rtp_pressure",
         "rtp_temperature"),
      GIN("ConstXsec"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("0."),
      GIN_DESC("Constant Xsec value")));

  md_data_raw.push_back(create_mdrecord(
      NAME("gas_scattering_matIsotropic"),
      DESCRIPTION(R"--(Calculates the spectrum of normalized scattering matrices.
Important, the angular direction are line of sight direction not the
propagation direction.
)--"),
      AUTHORS("Manfred Brath"),
      OUT("gas_scattering_mat",
          "gas_scattering_fct_legendre"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("gas_scattering_los_in",
         "gas_scattering_los_out",
         
         "gas_scattering_output_type"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("gas_scattering_matRayleigh"),
      DESCRIPTION(R"--(Calculates the normalized Rayleigh scattering matrix.

The phase matrix for anisotropic Rayleigh particles in random orientations.
Important, the angular direction are defined as line of sight direction not as
propagation direction.
)--"),
      AUTHORS("Jon Petersen"),
      OUT("gas_scattering_mat", "gas_scattering_fct_legendre"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("gas_scattering_los_in", "gas_scattering_los_out",  "gas_scattering_output_type"),
      GIN("depolarization_factor"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("0.03"),
      GIN_DESC("depolarization factor for air")));

  md_data_raw.push_back(create_mdrecord(
      NAME("g0Earth"),
      DESCRIPTION(R"--(Gravity at zero altitude on Earth.

Sets *g0* for the given latitude using a standard parameterisation.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("g0"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("lat"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(NAME("g0Io"),
                                 DESCRIPTION(R"--(Gravity at zero altitude on Io.

Numeric from Wikipedia.
)--"),
                                 AUTHORS("Richard Larsson"),
                                 OUT("g0"),
                                 GOUT(),
                                 GOUT_TYPE(),
                                 GOUT_DESC(),
                                 IN(),
                                 GIN(),
                                 GIN_TYPE(),
                                 GIN_DEFAULT(),
                                 GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("g0Jupiter"),
      DESCRIPTION(R"--(Gravity at zero altitude on Jupiter.

Sets *g0*  to mean equatorial gravity on Jupiter. Value provided by
MPS under ESA-planetary study (TN1).
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("g0"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("g0Mars"),
      DESCRIPTION(R"--(Gravity at zero altitude on Mars.

Sets *g0*  to mean equatorial gravity on Mars. Value provided by
MPS under ESA-planetary study (TN1).
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("g0"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("g0Venus"),
      DESCRIPTION(R"--(Gravity at zero altitude on Venus.

Sets *g0*  to mean equatorial gravity on Venus. Value from Ahrens
(1995), provided by MPS under ESA-planetary study (TN1).
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("g0"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("geo_posEndOfPpath"),
      DESCRIPTION(R"--(Sets geo-position based on *ppath*.

The geo-position is set to the position of the last point in *ppath*.

NaN is returned if *ppath* is totally outside of the atmosphere.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("geo_pos"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("geo_posLowestAltitudeOfPpath"),
      DESCRIPTION(R"--(Sets geo-position based on *ppath*.

The geo-position is set to the position of the point in *ppath*
having the lowest altitude.

NaN is returned if *ppath* is totally outside of the atmosphere.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("geo_pos"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("geo_posWhereAltitudeIsPassed"),
      DESCRIPTION(R"--(Sets geo-position based on *ppath*.

The geo-position is set to the position where the propagation
path passes the stated altitude. If this altitude is passed
more than once, the passing closest to the sensor is selected.
If the reference altitude is not passed at all, *geo_pos* is
set to NaN.

NaN is also returned if *ppath* is totally outside of the atmosphere.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("geo_pos"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath"),
      GIN("altitude"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Altitude defining *geo_pos*.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("GetEnvironmentVariable"),
      DESCRIPTION(R"--(Copy the contents of an environment variable to an ARTS String or Index.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("String, Index"),
      GOUT_DESC("Contents of environment variable."),
      IN(),
      GIN("input"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Name of environment variable.")));

  md_data_raw.push_back(
      create_mdrecord(NAME("GetNumberOfThreads"),
               DESCRIPTION(R"--(Returns the number of threads used by ARTS.
)--"),
               AUTHORS("Oliver Lemke"),
               OUT(),
               GOUT("nthreads"),
               GOUT_TYPE("Index"),
               GOUT_DESC("Number of threads."),
               IN(),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord(NAME("GriddedFieldGetName"),
               DESCRIPTION(R"--(Get the name of a GriddedField.

See *ArrayOfGriddedFieldGetNames*.
)--"),
               AUTHORS("Lukas Kluft"),
               OUT(),
               GOUT("name"),
               GOUT_TYPE("String"),
               GOUT_DESC("Name of the GriddedField."),
               IN(),
               GIN("griddedfield"),
               GIN_TYPE("GriddedField1, GriddedField2, GriddedField3, "
                        "GriddedField4, GriddedField5, GriddedField6 "),
               GIN_DEFAULT(NODEF),
               GIN_DESC("GriddedField."),
               SETMETHOD(false),
               AGENDAMETHOD(false),
               USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("GriddedFieldLatLonExpand"),
      DESCRIPTION(R"--(Expands the latitude and longitude grid of the GriddedField to
[-90, 90] and [0,360], respectively.
Expansion is only done in
the dimension(s), where the grid size is 1.
The values from the input data will be duplicated to accomodate
for the larger size of the output field.
output and input can be the same variable.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE(
          "GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3"),
      GOUT_DESC("Expanded gridded field."),
      IN(),
      GIN("input"),
      GIN_TYPE(
          "GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Raw input gridded field.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("GriddedFieldLatLonRegrid"),
      DESCRIPTION(R"--(Interpolates the input field along the latitude and longitude dimensions
to *lat_true* and *lon_true*.

If the input longitude grid is outside of *lon_true* it will be shifted
left or right by 360. If it covers 360 degrees, a cyclic interpolation
will be performed.
input and output fields can be the same variable.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE(
          "GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3"),
      GOUT_DESC("Regridded gridded field."),
      IN("lat_true", "lon_true"),
      GIN("input", "interp_order"),
      GIN_TYPE(
          "GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3",
          "Index"),
      GIN_DEFAULT(NODEF, "1"),
      GIN_DESC("Raw input gridded field.", "Interpolation order.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("heating_ratesFromIrradiance"),
      DESCRIPTION(R"--(Calculates heating rates from the *irradiance_field*.

The method assumes that the heating rates depend only on the
vertical derivation of the net flux. The net flux is the sum of the
*irradiance_field* in upward direction and the *irradiance_field*
in downward direction
)--"),
      AUTHORS("Manfred Brath"),
      OUT("heating_rates"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppvar_atm", "irradiance_field", "specific_heat_capacity", "g0"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("HydrotableCalc"),
      DESCRIPTION(R"--(Creates a look-up table of scattering properties.

The table produced largely follows the format used in RTTOV-SCATT for
its \"hydrotables\". The table is returned as a GriddedField4, with
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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("hydrotable"),
      GOUT_TYPE("GriddedField4"),
      GOUT_DESC("Generated hydrotable with format described above."),
      IN("pnd_agenda_array",
         "pnd_agenda_array_input_names",
         "scat_data",
         "scat_data_checked",
         "f_grid"),
      GIN("iss", "T_grid", "wc_grid"),
      GIN_TYPE("Index", "Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Index of scattering species.",
               "Temperature grid of table.",
               "Water content grid of table.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Ignore"),
      DESCRIPTION(R"--(Ignore a workspace variable.

This method is handy for use in agendas in order to suppress warnings
about unused input workspace variables. What it does is: Nothing!
In other words, it just ignores the variable it is called on.

This method can ignore any workspace variable you want.
)--"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("input"),
      GIN_TYPE("Any"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Variable to be ignored."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("INCLUDE"),
      DESCRIPTION(R"--(Includes the contents of another controlfile.

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
    INCLUDE \"agendas.arts\"
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
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("IndexAdd"),
      DESCRIPTION(R"--(Adds a Index and a value (output = input + value).

The result can either be stored in the same or another Index.
)--"),
      AUTHORS("Patrick Eriksson, Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Index"),
      GOUT_DESC("Output Index."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Index.", "Value to add.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("IndexDivide"),
      DESCRIPTION(R"--(Integer division of a Index and a value (output = input / value).

Please note that integer divison is applied, and e.g. 5/3=1.

The result can either be stored in the same or another Index.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Index"),
      GOUT_DESC("Output Index."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Index (numerator).", "Denominator.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("IndexMultiply"),
      DESCRIPTION(R"--(Multiplies a Index and a value (output = input * value).

The result can either be stored in the same or another Index.
)--"),
      AUTHORS("Patrick Eriksson, Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Index"),
      GOUT_DESC("Output index."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Index.", "Multiplier.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("IndexSetToLast"),
      DESCRIPTION(R"--(Set an Index to point towards last position of array-type variables.

This method works as nelemGet, but gives the index number of the last
element (which equals nelem-1).
)--"),
      AUTHORS("Patrick Eriksson", "Oliver Lemke"),
      OUT("nelem"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("v"),
      GIN_TYPE((ARRAY_GROUPS + ", Vector").c_str()),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The method is defined for these groups."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("IndexStepDown"),
      DESCRIPTION(R"--(Performas: output = input - 1

Input and output can be same variable.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Index"),
      GOUT_DESC("Output index variable."),
      IN(),
      GIN("input"),
      GIN_TYPE("Index"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Input index variable.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("IndexStepUp"),
      DESCRIPTION(R"--(Performas: output = input + 1

Input and output can be same variable.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Index"),
      GOUT_DESC("Output index variable."),
      IN(),
      GIN("input"),
      GIN_TYPE("Index"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Input index variable.")));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("IndexSubtract"),
      DESCRIPTION(R"--(Subtracts a Index value (output = input - value).

The result can either be stored in the same or another Index.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Index"),
      GOUT_DESC("Output Index."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Index.", "Subtrahend.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("InterpAtmFieldToPosition"),
      DESCRIPTION(R"--(Point interpolation of atmospheric fields.

The default way to specify the position is by *rtp_pos*.

Linear interpolation is applied.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("AtmPoint"),
      GOUT_DESC("Value obtained by the interpolation."),
      IN("atm_field",
         "rtp_pos"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("InterpGriddedField2ToPosition"),
      DESCRIPTION(R"--(Latitude and longitude interpolation of a GriddedField2.

The default way to specify the position is by *rtp_pos*.

The interpolation is done for the latitude and longitude in
*rtp_pos*. The altitude in *rtp_pos* is completely ignored.
Linear interpolation is applied.

The input field (``gfield2``) is expected to have latitude and
longitude as first and second dimension.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Value obtained by interpolation."),
      IN( "lat_true", "lon_true", "rtp_pos"),
      GIN("gfield2"),
      GIN_TYPE("GriddedField2"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Gridded field to interpolate.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("InterpSurfaceFieldToPosition"),
      DESCRIPTION(R"--(Point interpolation of surface fields.

The default way to specify the position is by *rtp_pos*.

Linear interpolation is applied.

The interpolation is done for the latitude and longitude in
*rtp_pos*, while the altitude in *rtp_pos* is not part of the
calculations. However, it is checked that the altitude of *rtp_pos*
is inside the range covered by ``z_surface`` with a 1 m margin, to
give a warning when the specified position is not consistent with
the surface altitudes.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_point"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN( "rtp_pos", "surface_field", "surface_search_accuracy"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("InterpSurfaceTypeMask"),
      DESCRIPTION(R"--(Interpolation of surface type mask.

The method determines the surface type at the position of concern
(*rtp_pos*) from the provided type mask.

The surface type is taken as the nearest value in *surface_type_mask*.

The altitude in *rtp_pos* is ignored.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("surface_type"),
      GOUT_TYPE("Index"),
      GOUT_DESC("Surface type index"),
      IN(
         "lat_grid",
         "lat_true",
         "lon_true",
         "rtp_pos",
         "surface_type_mask"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("IntersectionGeometricAltitude"),
      DESCRIPTION(R"--(Calculates the geometrical intersection with an altitude.

For each observation geometry specified by the combination of
*sensor_pos* and *sensor_los*, the geometrical intersection with
an altitude is determined. The intersections are described by the
GOUT ``pos`` and ``los``.

For cases with no intersection, ``pos`` and ``los`` are filled with NaN.

The GOUT ``pos`` and ``los`` can NOT be *sensor_pos* and *sensor_los*.
If you want to store the intersections in *sensor_pos* and *sensor_los*
use *sensor_pos_losForwardToAltitude*. For *rte_pos* and *rte_los*
you have *rte_pos_losForwardToAltitude*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("pos", "los"),
      GOUT_TYPE("Matrix", "Matrix"),
      GOUT_DESC("Position of intersections.",
                "Line-of-sight at intersections."),
      IN("sensor_pos", "sensor_los", "surface_field"),
      GIN("altitude"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("0"),
      GIN_DESC("Target altitude.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("IntersectionGeometricLatitude"),
      DESCRIPTION(R"--(Calculates the geometrical intersection with a latitude.

For each observation geometry specified by the combination of
*sensor_pos* and *sensor_los*, the geometrical intersection with
a latitude is determined. The intersections are described by the
GOUT ``pos`` and *los.

For cases with no intersection, ``pos`` and ``los`` are filled with NaN.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("pos", "los"),
      GOUT_TYPE("Matrix", "Matrix"),
      GOUT_DESC("Position of intersections.",
                "Line-of-sight at intersections."),
      IN("sensor_pos", "sensor_los", "surface_field"),
      GIN("latitude"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Target latitude.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("IntersectionGeometricLongitude"),
      DESCRIPTION(R"--(Calculates the geometrical intersection with a longitude.

For each observation geometry specified by the combination of
*sensor_pos* and *sensor_los*, the geometrical intersection with
a longitude is determined. The intersections are described by the
GOUT ``pos`` and *los.

For cases with no intersection, ``pos`` and ``los`` are filled with NaN.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("pos", "los"),
      GOUT_TYPE("Matrix", "Matrix"),
      GOUT_DESC("Position of intersections.",
                "Line-of-sight at intersections."),
      IN("sensor_pos", "sensor_los", "surface_field"),
      GIN("longitude"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Target longitude.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("IntersectionGeometricSurface"),
      DESCRIPTION(R"--(Calculates the geometrical intersection with the surface.

For each observation geometry specified by the combination of
*sensor_pos* and *sensor_los*, the geometrical intersection with
the surface is determined. The intersections are described by the
GOUT ``pos`` and *los. For cases with no intersection, ``pos`` and ``los``
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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("pos", "los"),
      GOUT_TYPE("Matrix", "Matrix"),
      GOUT_DESC("Position of intersections.",
                "Line-of-sight at intersections."),
      IN("sensor_pos",
         "sensor_los",
         "surface_field",
         "surface_search_accuracy",
         "surface_search_safe"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("irradiance_fieldFromRadiance"),
      DESCRIPTION(R"--(Calculates the irradiance from the *radiance_field*.

The *radiance_field* is integrated over the angular grids according to
the grids set by *AngularGridsSetFluxCalc*. 
See *AngularGridsSetFluxCalc* to set *za_grid*, *aa_grid*, and
*za_grid_weights*
)--"),
      AUTHORS("Manfred Brath"),
      OUT("irradiance_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("radiance_field", "za_grid", "aa_grid", "za_grid_weights"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("isotopologue_ratiosInitFromBuiltin"),
      DESCRIPTION(R"--(Initialize isotopologue ratios with default values from built-in
species data.  This should be OK for Earth-like atmospheres
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("isotopologue_ratios"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("isotopologue_ratiosInitFromHitran"),
      DESCRIPTION(R"--(Initialize isotopologue ratios with default values from built-in
Hitran species data.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("isotopologue_ratios"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iyApplyUnit"),
      DESCRIPTION(R"--(Conversion of *iy* to other spectral units (for passive observations).

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
    \"iy\", \"Error\" and \"Error (uncorrelated)\"

Please note that *diy_dx* is not handled. Also note that this method
considers *iy_unit*, while *iy_unit_radar* is handled directly by
the methods dealing with such simulations.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy", "iy_aux"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy", "iy_aux",  "f_grid", "iy_aux_vars", "iy_unit"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord

      (NAME("iyBackground"),
       DESCRIPTION(R"--(Computes background radiation by wrapping various agendas
)--"),
       AUTHORS("Richard Larsson"),
       OUT("iy", "diy_dx"),
       GOUT(),
       GOUT_TYPE(),
       GOUT_DESC(),
       IN("iy_transmittance",
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
          "iy_unit"),
       GIN(),
       GIN_TYPE(),
       GIN_DEFAULT(),
       GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord

      (NAME("iyCalc"),
       DESCRIPTION(R"--(A single monochromatic pencil beam calculation.

Performs monochromatic radiative transfer calculations for the
specified position (*rte_pos*) and line-of-sight (*rte_pos*).
See *iy* and associated variables for format of output.

Please note that Jacobian type calculations not are supported.
For this use *yCalc*.

No sensor characteristics are applied. These are most easily
incorporated by using *yCalc*
)--"),
       AUTHORS("Patrick Eriksson"),
       OUT("iy", "iy_aux", "ppath", "geo_pos"),
       GOUT(),
       GOUT_TYPE(),
       GOUT_DESC(),
       IN("atmgeom_checked",
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
          "iy_main_agenda"),
       GIN(),
       GIN_TYPE(),
       GIN_DEFAULT(),
       GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iyClearsky"),
      DESCRIPTION(R"--(Standard method for radiative transfer calculations with emission
and a direct (solar) source.

Designed to be part of *iy_main_agenda*. That is, only valid
outside the cloudbox (no scattering). For details se the user guide.

The possible choices for *iy_unit* are

- ``\"1\"``: No conversion, i.e. [W/(m^2 Hz sr)] (radiance per frequency unit).
- ``\"RJBT\"``: Conversion to Rayleigh-Jean brightness temperature.
- ``\"PlanckBT\"``: Conversion to Planck brightness temperature.
- ``\"W/(m^2 m sr)\"``: Conversion to [W/(m^2 m sr)] (radiance per wavelength unit).
- ``\"W/(m^2 m-1 sr)\"``: Conversion to [W/(m^2 m-1 sr)] (radiance per wavenumber unit).

Expressions applied and considerations for the unit conversion of
radiances are discussed in Sec. 5.7 of the ARTS-2.0 article.

*iy_unit* is only applied if *iy_agenda_call1* is 1. This means that
no unit ocnversion is applied for internal iterative calls.

Recognised choices for *rt_integration_option* are:

- ``\"first order\"``: A first order integration is applied.
- ``\"second order\"``: A second order integration is applied.
- ``\"default\"``: Another way to select the first order option.

Some auxiliary radiative transfer quantities can be obtained. Auxiliary
quantities are selected by *iy_aux_vars* and returned by *iy_aux*.
Valid choices for auxiliary data are:

- ``\"Radiative background\"``:
    Index value flagging the radiative
    background. The following coding is used: 0=space, 1=surface
    and 2=cloudbox.
- ``\"Optical depth\"``:
    Scalar optical depth between the observation point
    and the end of the present propagation path. Calculated based on
    the (1,1)-element of the transmittance matrix (1-based indexing),
    i.e. only fully valid for scalar RT.
- ``\"Direct radiation\"``:
    Stokes vector of direct radiation. It dimensions
    are number of frequencies and ``stokes_dim``. If no sun is present 
    in the line of sight, it is zero.
- ``\"Radiation Background\"``:
    Stokes vector of the radiation at start of
    the propagation path. It dimensions are number of frequencies and
    ``stokes_dim``.

If nothing else is stated, only the first column of *iy_aux* is filled,
i.e. the column matching Stokes element I, while remaing columns are
are filled with zeros.

IMPORTANT:
    No jacobian calculation is supported when suns or gas 
    scattering is included! This will be implemented in a future version.
)--"),
      AUTHORS("Patrick Eriksson", "Richard Larsson", "Oliver Lemke", "Manfred Brath"),
      OUT("iy",
          "iy_aux",
          "diy_dx",
          "ppvar_atm",
          "ppvar_f",
          "ppvar_iy",
          "ppvar_trans_cumulat",
          "ppvar_trans_partial"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
         "iy_id",
         
         "f_grid",
         "abs_species",
         "atm_field",
         "surface_field",
         "ppath_lmax",
         "ppath_lraytrace",
         "cloudbox_on",
         "gas_scattering_do",
         "suns_do",
         "iy_unit",
         "iy_aux_vars",
         "jacobian_do",
         "jacobian_quantities",
         "ppath",
         "rte_pos2",
         "suns",
         "propmat_clearsky_agenda",
         "water_p_eq_agenda",
         "rt_integration_option",
         "iy_main_agenda",
         "iy_space_agenda",
         "iy_surface_agenda",
         "iy_cloudbox_agenda",
         "gas_scattering_agenda",
         "ppath_step_agenda",
         "iy_agenda_call1",
         "iy_transmittance",
         "rte_alonglos_v"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iyEmissionHybrid"),
      DESCRIPTION(R"--(Radiative transfer with emission and precalculated radiation field.

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
)--"),
      AUTHORS("Patrick Eriksson", "Jana Mendrok", "Richard Larsson"),
      OUT("iy",
          "iy_aux",
          "diy_dx",
          "ppvar_atm",
          "ppvar_pnd",
          "ppvar_f",
          "ppvar_iy",
          "ppvar_trans_cumulat",
          "ppvar_trans_partial"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
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
         "za_grid"),
      GIN("Naa_grid", "t_interp_order"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT("19", "1"),
      GIN_DESC("Number of azimuth angles to consider in scattering source term"
               " integral.",
               "Interpolation order of temperature for scattering data (so"
               " far only applied in phase matrix, not in extinction and"
               " absorption.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("iyIndependentBeamApproximation"),
      DESCRIPTION(R"--(Samples atmosphere along ppath and make 1D-type RT calculation.

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
\"Scattering element 0\".
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy", "iy_aux", "ppath", "diy_dx", "atm_fields_compact"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
         "iy_id",
         "f_grid",
         "atm_field",
         "cloudbox_on",
         "cloudbox_limits",
         "pnd_field",
         "particle_masses",
         "ppath_agenda",
         "ppath_lmax",
         "ppath_lraytrace",
         "iy_agenda_call1",
         "iy_unit",
         "iy_transmittance",
         "rte_pos",
         "rte_los",
         "rte_pos2",
         "jacobian_do",
         "iy_aux_vars",
         "iy_independent_beam_approx_agenda"),
      GIN("return_atm1d", "skip_vmr", "skip_pnd", "return_masses"),
      GIN_TYPE("Index", "Index", "Index", "Index"),
      GIN_DEFAULT("0", "0", "0", "0"),
      GIN_DESC(
          "Flag to trigger that *atm_fields_compact* is filled. ",
          "Flag to not include vmr data in *atm_fields_compact*.",
          "Flag to not include pnd data in *atm_fields_compact*.",
          "Flag to include particle category masses in *atm_fields_compact*."
          "Conversion is done by *particle_masses*.")));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("iyInterpCloudboxField"),
      DESCRIPTION(R"--(Interpolates the intensity field of the cloud box.

Determines the intensity field at the position and direction
specified by *rte_pos* and *rte_los*. The position can be both
inside the cloud box or at its edge.

The interpolation in the spatial dimensions is linear.

For the zenith angle dimensions several options for controlling
the interpolation are at hand. Default is linear interpolation.
Higher order polynomial interpolation is activated by setting
``za_interp_order`` to a value > 1. Default is to perform the
interpolation separately for [0,90[ and ]90,180]. To handle
90 degree or use the full range ([0,180]) as basis for the
interpolation, set ``za_restrict`` to 0. You can select to use
cos(za) as the independent variable (instead of za) by setting
``cos_za_interp`` to 1.

For the azimuth dimension the interpolation order can be
selected, in the same manner as for zenith.
)--"),
      AUTHORS("Claudia Emde", "Patrick Eriksson", "Jana Mendrok"),
      OUT("iy"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("cloudbox_field",
         "rtp_pos",
         "rtp_los",
         "jacobian_do",
         "cloudbox_on",
         "cloudbox_limits",
         "atm_field",
         "surface_field",
         
         "za_grid",
         "aa_grid",
         "f_grid"),
      GIN("za_interp_order",
          "za_restrict",
          "cos_za_interp",
          "za_extpolfac",
          "aa_interp_order"),
      GIN_TYPE("Index", "Index", "Index", "Numeric", "Index"),
      GIN_DEFAULT("1", "1", "0", "0.5", "1"),
      GIN_DESC("Zenith angle interpolation order.",
               "Flag whether to restric zenith angle interpolation to one "
               "hemisphere.",
               "Flag whether to do zenith angle interpolation in cosine space.",
               "Maximum allowed extrapolation range in zenith angle.",
               "Azimuth angle interpolation order.")));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("iyLoopFrequencies"),
      DESCRIPTION(R"--(Radiative transfer calculations one frequency at the time.

The method loops the frequencies in *f_grid* and calls
*iy_loop_freqs_agenda* for each individual value. This method is
placed in *iy_main_agenda*, and the actual radiative transfer
method is put in *iy_loop_freqs_agenda*.

A common justification for using the method should be to consider
dispersion. By using this method it is ensured that the propagation
path for each individual frequency is calculated.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy", "iy_aux", "ppath", "diy_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy_aux_vars",
         "iy_agenda_call1",
         "iy_transmittance",
         "rte_pos",
         "rte_los",
         "rte_pos2",
         
         "f_grid",
         "iy_loop_freqs_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iyMC"),
      DESCRIPTION(R"--(Interface to Monte Carlo part for *iy_main_agenda*.

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

- ``\"Error (uncorrelated)\"``:
    Calculation error. Size: [nf,ns,1,1].
    (The later part of the text string is required. It is used as
    a flag to yCalc for how to apply the sensor data.)

where

- `nf`: Number of frequencies.
- `ns`: Number of Stokes elements.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy", "iy_aux", "diy_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy_agenda_call1",
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
         "mc_taustep_limit"),
      GIN("t_interp_order"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Interpolation order of temperature for scattering data (so"
               " far only applied in phase matrix, not in extinction and"
               " absorption.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("iyRadarSingleScat"),
      DESCRIPTION(R"--(Simulation of radar, restricted to single scattering.

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
applied in *yRadar*. The output of this method matches the option \"1\".

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

- ``\"Radiative background\"``:
    Index value flagging the radiative
    background. The following coding is used: 0=space, 1=surface
    and 2=cloudbox (the last case should not occur!). Only column
    matching first Stokes element filled. Other columns are set to 0.
- ``\"Backscattering\"``:
    The unattenuated back-scattering. That is, as
    *iy* but with no attenuated applied. Here all columns are filled.
    By combing *iy* and this auxiliary variable, the total two-way
    attenuation can be derived.
- ``\"Abs species extinction\"``:
    Extinction due to *abs_species* at each
    ppath point, taken as the diagonal of the local extinction matrix.
- ``\"Particle extinction\"``:
    Extinction due to particles at each
    ppath point, taken as the diagonal of the local extinction matrix.
    The retunred values includes ``pext_scaling``
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy",
          "iy_aux",
          "diy_dx",
          "ppvar_atm",
          "ppvar_pnd",
          "ppvar_f"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         "f_grid",
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
         "rte_alonglos_v"),
      GIN("trans_in_jacobian", "pext_scaling", "t_interp_order"),
      GIN_TYPE("Index", "Numeric", "Index"),
      GIN_DEFAULT("0", "1", "1"),
      GIN_DESC("Flag determining if change in transmittance is considered"
               " in calculation of the Jacobian or not.",
               "Particle extinction is scaled with this value. A value"
               " inside [0,2]. Set it to 0 if you want to remove particle"
               " extinction totally.",
               "Interpolation order of temperature for scattering data (so"
               " far only applied in phase matrix, not in extinction and"
               " absorption.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("iyReplaceFromAux"),
      DESCRIPTION(R"--(Change of main output variable.

With this method you can replace the content of *iy* with one of
the auxiliary variables. The selected variable (by ``aux_var``) must
be part of *iy_aux_vars*. The corresponding data from *iy_aux* are
copied to form a new *iy* (*iy_aux* is left unchanged). Elements of
*iy* correponding to Stokes elements not covered by the auxiliary
variable are just set to zero.

Jacobian variables are not handled.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy", "iy_aux", "iy_aux_vars", "jacobian_do"),
      GIN("aux_var"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Auxiliary variable to insert as *iy*.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("iySurfaceFastem"),
      DESCRIPTION(R"--(Usage of FASTEM for emissivity and reflectivity of water surfaces.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy", "diy_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
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
         "surface_skin_t"),
      GIN("salinity", "wind_speed", "wind_direction", "fastem_version"),
      GIN_TYPE("Numeric", "Numeric", "Numeric", "Index"),
      GIN_DEFAULT("0.035", NODEF, "0", "6"),
      GIN_DESC("Salinity, 0-1. That is, 3% is given as 0.03.",
               "Wind speed.",
               "Wind direction. See further above.",
               "The version of FASTEM to use.")));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("iySurfaceFlatReflectivity"),
      DESCRIPTION(R"--(This method calculates upwelling radiation for a specular flat surface.

These are due to the reflection of the downgoing diffuse radiation and emission from
the surface using a predefined reflectivity matrix. 

This method is designed to be part of *iy_surface_agenda*

Important this method calculates only the reflection of the diffuse
downward radiation. No direct incoming radiation is considered

Jacobian is supported only for Skin temperature
)--"),
      AUTHORS("Manfred Brath"),
      OUT("iy",
          "diy_dx",
          "dsurface_rmatrix_dx",
          "dsurface_emission_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy",
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
         "lat_grid",
         "lon_grid",
         "surface_field",
         "rtp_pos",
         "rtp_los",
         "rte_pos2",
         "iy_unit",
         "surface_reflectivity",
         "surface_field",
         "surface_props_names",
         "dsurface_names",
         "jacobian_quantities",
         "iy_main_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("iySurfaceFlatReflectivityDirect"),
      DESCRIPTION(R"--(This method calculates the specular reflection at a flat 
surface of the direct radiation with a predefined reflectivity matrix.

This method is designed to be part of *iy_surface_agenda*

Important this method calculates only the scattering of the direct
(sun) radiation. No diffuse incoming radiation is considered

This method has no jacobian capability
)--"),
      AUTHORS("Manfred Brath"),
      OUT("iy"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy",
         "rtp_pos",
         "rtp_los",
         
         "f_grid",
         "abs_species",
         "atm_field",
         "surface_field",
         "surface_reflectivity",
         "pnd_field",
         "dpnd_field_dx",
         "scat_species",
         "scat_data",
         "ppath_lmax",
         "ppath_lraytrace",
         "ppath_inside_cloudbox_do",
         "cloudbox_on",
         "cloudbox_limits",
         "suns_do",
         "gas_scattering_do",
         "jacobian_do",
         "jacobian_quantities",
         "suns",
         "rte_alonglos_v",
         "iy_unit",
         "propmat_clearsky_agenda",
         "water_p_eq_agenda",
         "gas_scattering_agenda",
         "ppath_step_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iySurfaceFlatRefractiveIndex"),
      DESCRIPTION(R"--(This method calculates upwelling radiation for a specular flat surface.

These are due to the reflection of the downgoing diffuse radiation and emission from
the surface using a predefined reflectivity matrix. 

This method is designed to be part of *iy_surface_agenda*

Important this method calculates only the reflection of the diffuse
downward radiation. No direct incoming radiation is considered

Jacobian is supported only for Skin temperature
)--"),
      AUTHORS("Manfred Brath"),
      OUT("iy",
          "diy_dx",
          "dsurface_rmatrix_dx",
          "dsurface_emission_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy",
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
         "iy_main_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iySurfaceFlatRefractiveIndexDirect"),
      DESCRIPTION(R"--(This method calculates the specular reflection at a flat 
surface of the direct radiation.

This method is designed to be part of *iy_surface_agenda*

Important this method calculates only the scattering of the direct
(sun) radiation. No diffuse incoming radiation is considered

This method has no jacobian capability
)--"),
      AUTHORS("Manfred Brath"),
      OUT("iy"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy",
         "rtp_pos",
         "rtp_los",
         
         "f_grid",
         "abs_species",
         "atm_field",
         "surface_field",
         "surface_complex_refr_index",
         "pnd_field",
         "dpnd_field_dx",
         "scat_species",
         "scat_data",
         "ppath_lmax",
         "ppath_lraytrace",
         "ppath_inside_cloudbox_do",
         "cloudbox_on",
         "cloudbox_limits",
         "suns_do",
         "gas_scattering_do",
         "jacobian_do",
         "jacobian_quantities",
         "suns",
         "rte_alonglos_v",
         "iy_unit",
         "propmat_clearsky_agenda",
         "water_p_eq_agenda",
         "gas_scattering_agenda",
         "ppath_step_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iySurfaceInit"),
      DESCRIPTION(R"--(This method initialize iy.

This method is designed to be part of *iy_surface_agenda*.
Its only prpose is to initialize *iy* properly within the 
*iy_surface_agenda*
)--"),
      AUTHORS("Manfred Brath"),
      OUT("iy"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iySurfaceLambertian"),
      DESCRIPTION(R"--(This method calculates upwelling radiation for a lambertian surface.

These are due to the scattering of the downgoing diffuse radiation and emission from
the surface.
This method works only for 1D or 3D atmospheres.
For the integration over the zenith angles a gaussian quadrature with
N_za angles is used.
For 1D atmospheres N_aa is ignored. For 3D atmospheres without clouds
azimuthal dependency can be neglected. N_aa = 1 is sufficient.
For 3D atmospheres with cloudbox on azimuthal dependency needs to be 
accounted. In that case the number of azimuth angles N_aa as a rule of
thumb should be set to 4*N_za.
For the 1D case N_za downwelling streams and 3D case N_za*N_aa downwelling
streams are calculated.

This method is designed to be part of *iy_surface_agenda*

Important this method calculates only the scattering of the diffuse
downward radiation. No direct incoming radiation is considered

Jacobian is supported only for Skin temperature
)--"),
      AUTHORS("Manfred Brath"),
      OUT("iy",
          "diy_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy",
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
         "iy_main_agenda"),
      GIN("N_za","N_aa"),
      GIN_TYPE("Index","Index"),
      GIN_DEFAULT("3","1"),
      GIN_DESC("Number of zenith angles.","Number of azimuth angles")));

  md_data_raw.push_back(create_mdrecord(
      NAME("iySurfaceLambertianDirect"),
      DESCRIPTION(R"--(This method calculates the scattering of the direct radiation
for a Lambertian surface.

This method is designed to be part of *iy_surface_agenda*

Important this method calculates only the scattering of the direct
(sun) radiation. No diffuse incoming radiation is considered

This method has no jacobian capability
)--"),
      AUTHORS("Manfred Brath"),
      OUT("iy"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy",
         "rtp_pos",
         
         "f_grid",
         "abs_species",
         "atm_field",
         "surface_field",
         "surface_scalar_reflectivity",
         "pnd_field",
         "dpnd_field_dx",
         "scat_species",
         "scat_data",
         "ppath_lmax",
         "ppath_lraytrace",
         "cloudbox_on",
         "cloudbox_limits",
         "suns_do",
         "gas_scattering_do",
         "jacobian_do",
         "jacobian_quantities",
         "suns",
         "rte_alonglos_v",
         "iy_unit",
         "propmat_clearsky_agenda",
         "water_p_eq_agenda",
         "gas_scattering_agenda",
         "ppath_step_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iySurfaceRtpropAgenda"),
      DESCRIPTION(R"--(Interface to *surface_rtprop_agenda* for *iy_surface_agenda*.

This method is designed to be part of *iy_surface_agenda*. It
determines the radiative properties of the surface by
*surface_rtprop_agenda* and calculates the downwelling radiation
by *iy_main_agenda*, and sums up the terms as described in AUG.
That is, this WSM uses the output from *surface_rtprop_agenda*
in a straightforward fashion.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy",
          "diy_dx",
          "surface_point",
          "surface_los",
          "surface_rmatrix",
          "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
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
         "surface_rtprop_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iySurfaceRtpropCalc"),
      DESCRIPTION(R"--(Applies *surface_los*, *surface_rmatrix* and *surface_emission*.

This method is designed to be part of *iy_surface_agenda* and
should be mandatory when using methods describing the surface
radiative transfer properties by *surface_los*, *surface_rmatrix*
and *surface_emission*. The task of this method is to apply these
three WSVs to obtain the upwelling radiation from the surface.
This upwelling radiation is the sum of surface emission and
reflected downwelling radiation. The later part is calculated
by calling *iy_main_agenda*. See further AUG.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy", "diy_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
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
         "iy_main_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iyTransmissionStandard"),
      DESCRIPTION(R"--(Standard method for handling transmission measurements.

Designed to be part of *iy_main_agenda*. Treatment of the cloudbox
is incorporated (that is, no need to define *iy_cloudbox_agenda*).

The transmitter is assumed to be placed at the end of provided *ppath*.
The transmitted signal is taken from *iy_transmitter*. This
signal is propagated along the path, considering attenuation alone.
That is, the result of the method (*iy*) is the output of
*iy_transmitter* multiplied with the transmittance along the
propagation path.

As mentioned, the given *ppath* determines the position of the
transmitter. For clear-sky and no modification of *ppath*, this
means that the transitter will either be found at the surface or
at the top-of-the-atmosphere. If you want to maintain this even with
an active cloudbox, calculate *ppath* as::

     ppathCalc( cloudbox_on=0 )

Without setting cloudbox_on=0, the transmitter will end up inside or
at the boundary of the cloudbox.

Some auxiliary radiative transfer quantities can be obtained. Auxiliary
quantities are selected by *iy_aux_vars* and returned by *iy_aux*.
Valid choices for auxiliary data are:

- ``\"Radiative background\"``:
    Index value flagging the radiative
    background. The following coding is used: 0=space, 1=surface
    and 2=cloudbox. The value is added to each column.
- ``\"Optical depth\"``:
    Scalar optical depth between the observation point
    and the end of the present propagation path. Calculated based on
    the (1,1)-element of the transmittance matrix (1-based indexing),
    i.e. only fully valid for scalar RT. The value is added to each
    column.

IMPORTANT:
    No jacobian calculation is supported when gas scattering is
    included! This will be implemented in a future version.
)--"),
      AUTHORS("Patrick Eriksson", "Richard Larsson"),
      OUT("iy",
          "iy_aux",
          "diy_dx",
          "ppvar_atm",
          "ppvar_pnd",
          "ppvar_f",
          "ppvar_iy",
          "ppvar_trans_cumulat",
          "ppvar_trans_partial"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
         
         "f_grid",
         "abs_species",
         "atm_field",
         "cloudbox_on",
         "cloudbox_limits",
         "gas_scattering_do",
         "pnd_field",
         "dpnd_field_dx",
         "scat_species",
         "scat_data",
         "iy_aux_vars",
         "jacobian_do",
         "jacobian_quantities",
         "ppath",
         "iy_transmitter",
         "propmat_clearsky_agenda",
         "water_p_eq_agenda",
         "gas_scattering_agenda",
         "iy_agenda_call1",
         "iy_transmittance",
         "rte_alonglos_v"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iy_transmitterMultiplePol"),
      DESCRIPTION(R"--(Transmitted signal having multiple polarisations.

The method is intended to be used as possible input of 
*iyTransmissionStandard*.
It sets *iy_transmitter* to describe the transmitted signal/pulses.
The polarisation state is taken from *instrument_pol*, where
*instrument_pol* must contain an element for each frequency in *f_grid*.
The transmitted signal/pulses are set to be of unit magnitude, such
as [1,1,0,0].
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy_transmitter"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN( "f_grid", "instrument_pol"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iy_transmitterSinglePol"),
      DESCRIPTION(R"--(Transmitted signal having a single polarisations.

The method is intended to be used as possible input of 
*iyTransmissionStandard*.
It sets *iy_transmitter* to describe the transmitted signal/pulses.
The polarisation state is taken from *instrument_pol*, where
*instrument_pol* must contain a single value.
This polarisation state is applied for all frequencies.
The transmitted pulses/signals are set to be of unit
magnitude, such as [1,1,0,0].
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy_transmitter"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN( "f_grid", "instrument_pol"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddAbsSpecies"),
      DESCRIPTION(R"--(Includes an absorption species in the Jacobian.

For 1D or 2D calculations the latitude and/or longitude grid of
the retrieval field should set to have zero length.

These retrieval units are at hand for all gas species:

- ``\"vmr\"``: Volume mixing ratio.
- ``\"nd\"``: Number density.
- ``\"rel\"``: Relative unit (e.g. 1.1 means 10% more of the gas).

For water vapour, also these units are at hand:

- ``\"rh\"``: Relative humidity.
- ``\"q\"``: Specific humidity.

Note that ``for_species_tag`` is used to indicate if species tag VMR,
rather than atmospheric gas VMR is calculated. Set it to 0 and we
calculate the atmospheric gas VMR, but this only works for \"analytical\".

Note that the Jacobian is set to zero where volume mixing ratio equals zero.

The number of elements added to the state vector (*x*) is::

  n_g1 * n_g2 * n_g3

where n_g1, n_g2 and n_g3 are the length of GIN ``g1``, ``g2`` and ``g3``,
respectively. Here empty vectors should be considered to have a length 1.
The elements are sorted with pressure as innermost loop, followed by
latitude and longitude as outermost loop.
)--"),
      AUTHORS("Mattias Ekstrom", "Patrick Eriksson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda"),
      GIN("g1", "g2", "g3", "species", "unit", "for_species_tag"),
      GIN_TYPE("Vector", "Vector", "Vector", "String", "String", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, "vmr", "1"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "The species tag of the retrieval quantity.",
               "Retrieval unit. See above.",
               "Index-bool for acting on species tags or species."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddBasicCatalogParameter"),
      DESCRIPTION(R"--(Includes a basic catalog parameter in the Jacobian. These are constant
over all layers and so only a single vector output is returned.

The only basic catalog parameters currently supported are:

* ``\"LineStrength\"``
* ``\"LineCenter\"``

The ``catalog_identity`` should be able to identify one or many
lines in the catalog used for calculating the spectral absorption.
Note that partial matching for energy levels are allowed but not
recommended, as it is somewhat nonsensical to add multiple parameters.

Also note *jacobianAddShapeCatalogParameter* as this allows addition
of shape parameters, e.g., pressure broadening coefficients.

Each call to this function adds just a single value to *x*.

Example given the catalog_identity=\"O2-66 TR UP v1 0 J 1 LO v1 0 J 0\",
only the O2 ground-level 119 GHz line can be accessed and only its
catalog_parameter will be accessed.  However, the more lenient
catalog_identity=\"O2-66 TR UP J 1 LO J 0\" may be used, but then the
118 GHz line belonging to v1=1 branch will be added to the same *x*.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities", "jacobian_agenda"),
      GIN("catalog_identity", "catalog_parameter"),
      GIN_TYPE("QuantumIdentifier", "String"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The catalog line matching information.",
               "The catalog parameter of the retrieval quantity.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddBasicCatalogParameters"),
      DESCRIPTION(R"--(See *jacobianAddBasicCatalogParameter*.

This adds a multiple of parameters for first each catalog identity in
``catalog_identities`` and then for each catalog parameter in
``catalog_parameters`` by looping calls to *jacobianAddBasicCatalogParameter*
over these input.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities", "jacobian_agenda"),
      GIN("catalog_identities", "catalog_parameters"),
      GIN_TYPE("ArrayOfQuantumIdentifier", "ArrayOfString"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The catalog line matching information.",
               "The catalog parameter of the retrieval quantity.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddFreqShift"),
      DESCRIPTION(R"--(Includes a frequency fit of shift type in the Jacobian.

Retrieval of deviations between nominal and actual backend
frequencies can be included by this method. The assumption here is
that the deviation is a constant off-set, a shift, common for all
frequencies (and not varying between measurement blocks).

This method adds one element to the state vector (*x*).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities", "jacobian_agenda", "f_grid"),
      GIN("df"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("100e3"),
      GIN_DESC("Size of perturbation to apply.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddFreqStretch"),
      DESCRIPTION(R"--(Includes a frequency fit of stretch type in the Jacobian.

Retrieval of deviations between nominal and actual backend
frequencies can be included by this method. The assumption here is
that the deviation varies linearly over the frequency range
(following ARTS basis function for polynomial order 1).

This method adds one element to the state vector (*x*).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities", "jacobian_agenda", "f_grid"),
      GIN("df"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("100e3"),
      GIN_DESC("Size of perturbation to apply.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddShapeCatalogParameter"),
      DESCRIPTION(R"--(Adds a line shape parameter to the Jacobian calculations. These
are constant over all levels so only a single *x*-value is added

Line function parameter assume the derivatives of internal pressure
broadening and line mixing functionality follows a f(T, T0, X0, X1, X2)
format. The shape of the function f() is determined by input
catalog; please see the ARTS documentation for more details.

The input are as follows:

- line_identity:
    Identifier of preferably a single line
- species:
    A SpeciesTag, e.g., \"O2\" or \"H2O\" for common species.
    Note that \"SELF\" and \"AIR\" tags are used for shape parameters
    affected by self and air-broadening, respectively.
- variable:
    A variable supported by the line, these can be

    - ``\"G0\"``:  Speed-independent pressure broadening
    - ``\"G2\"``:  Speed-dependent pressure broadening
    - ``\"D0\"``:  Speed-independent pressure shift
    - ``\"D2\"``:  Speed-dependent pressure shift
    - ``\"FVC\"``: Frequency of velocity changing collisions
    - ``\"ETA\"``: partial correlation between velocity and rotational state changes due to collisions
    - ``\"Y\"``:   First order line-mixing parameter
    - ``\"G\"``:   Second order line-mixing parameter for strength
    - ``\"DV\"``:  Second order line-mixing parameter for shifting
- coefficient:
    A coefficient in the model to compute the above parameters.

Note that we cannot test if the line in question supports the variable and
coefficient at the level of this function, so many errors will only be reported
at a later stage.

For other spectroscopic parameters, see *jacobianAddBasicCatalogParameter*.
Also see said function for an example of how to set the QuantumIdentifier.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities", "jacobian_agenda"),
      GIN("line_identity", "species", "variable", "coefficient"),
      GIN_TYPE("QuantumIdentifier", "String", "String", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF),
      GIN_DESC("Line identifier",
               "Species of interest",
               "Variable of interest",
               "Coefficient of interest")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddShapeCatalogParameters"),
      DESCRIPTION(R"--(See *jacobianAddShapeCatalogParameter* for information on
the GIN parameters

This function accepts the same input but for lists of data.
The function loops over each input list
individually and appends the information to *jacobian_quantities*.

Special \"ALL\" for 1 length ``variables`` and ``coefficients`` are
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
)--"),
      AUTHORS("Richard Larsson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities", "jacobian_agenda"),
      GIN("line_identities", "species", "variables", "coefficients"),
      GIN_TYPE("ArrayOfQuantumIdentifier",
               "ArrayOfString",
               "ArrayOfString",
               "ArrayOfString"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF),
      GIN_DESC("List of line identifiers",
               "List of species of interest",
               "List of variables of interest",
               "List of coefficients of interest")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddMagField"),
      DESCRIPTION(R"--(Includes one magnetic field component in the Jacobian.

The method follows the pattern of other Jacobian methods. The
calculations can only be performed by analytic expressions.

The magnetic field components are retrieved separately, and,
hence, the argument ``component`` can be  \"u\", \"v\", \"w\",
and \"strength\".

The number of elements added to the state vector (*x*) is::

  n_g1 * n_g2 * n_g3

where n_g1, n_g2 and n_g3 are the length of GIN ``g1``, ``g2`` and ``g3``,
respectively. Here empty vectors should be considered to have a length 1.
The elements are sorted with pressure as innermost loop, followed by
latitude and longitude as outermost loop.

The dB-parameter is only used for Faraday rotation.
)--"),
      AUTHORS("Patrick Eriksson", "Richard Larsson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda"),
      GIN("g1", "g2", "g3", "component", "dB"),
      GIN_TYPE("Vector", "Vector", "Vector", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "v", "1.0e-7"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Magnetic field component to retrieve",
               "Magnetic field perturbation")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddNLTE"),
      DESCRIPTION(R"--(Experimental NLTE Jacobian.

Intention: Adds the nlte_field level distribution per atmospheric grid
to the Jacobian.

The number of elements added to the state vector (*x*) is::

  n_g1 * n_g2 * n_g3

where n_g1, n_g2 and n_g3 are the length of GIN ``g1``, ``g2`` and ``g3``,
respectively. Here empty vectors should be considered to have a length 1.
The elements are sorted with pressure as innermost loop, followed by
latitude and longitude as outermost loop.

The QuantumIdentifier should identify a single energy level, such as:
\"H2O-161 EN J 1 Ka 0 Kc 1\", for one of the lower levels in the chains
of transitions of water.  Note that using this method directly is not
best practice, as the quantum identifiers of the levels have to be known
at an early stage in NLTE calculations, and will usually populate the
``nlte_level_identifiers`` variable, meaning it is better to use *jacobianAddNLTE*
directly than to individually call this function.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda"),
      GIN("g1", "g2", "g3", "energy_level_identity", "dx"),
      GIN_TYPE("Vector", "Vector", "Vector", "QuantumIdentifier", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, "1.0e-3"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Identifier to the eneregy level",
               "Perturbation of value if required by method")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddNLTEs"),
      DESCRIPTION(R"--(Experimental NLTE Jacobian.  Same as *jacobianAddNLTE* but for
many levels

Adds energy_level_identities.nelem() times as many arguments to *x*
as *jacobianAddNLTE*, ordered as energy_level_identities describes

This method is preferred to *jacobianAddNLTE*, since ``energy_level_identities``
is conveniently almost always the same as ``nlte_level_identifiers``.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda"),
      GIN("g1", "g2", "g3", "energy_level_identities", "dx"),
      GIN_TYPE(
          "Vector", "Vector", "Vector", "ArrayOfQuantumIdentifier", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, "1.0e-3"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Identifiers to the eneregy level",
               "Perturbation of value if required by method")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddPointingZa"),
      DESCRIPTION(R"--(Adds sensor pointing zenith angle off-set jacobian.

Retrieval of deviations between nominal and actual zenith angle of
the sensor can be included by this method. The weighing functions
can be calculated in several ways:

- ``calcmode = \"recalc\"``:
    Recalculation of pencil beam spectra,
    shifted with ``dza`` from nominal values. A single-sided
    perturbation is applied (towards higher zenith angles).
- ``calcmode = \"interp\"``:
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
)--"),
      AUTHORS("Patrick Eriksson", "Mattias Ekstrom"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities", "jacobian_agenda", "sensor_pos", "sensor_time"),
      GIN("poly_order", "calcmode", "dza"),
      GIN_TYPE("Index", "String", "Numeric"),
      GIN_DEFAULT("0", "recalc", "0.01"),
      GIN_DESC("Order of polynomial to describe the time variation of "
               "pointing off-sets.",
               "Calculation method. See above",
               "Size of perturbation to apply (when applicable).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddPolyfit"),
      DESCRIPTION(R"--(Includes polynomial baseline fit in the Jacobian.

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
all are set to 1, even if several spectra are involved. Otherwise the
number of elements added to *x* depends on the number of spectra and
the settings of ``no_pol_variation``, ``no_los_variation`` and 
``no_mblock_variation``. The coefficients of the different polynomial
orders are treated as separate retrieval quantities. That is, the
the elements associated with polynomial order 0 are grouped and form
together a retrieval quantity. The coefficients for higher polynomial
orders are treated in the same way.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid",
         "sensor_pos"),
      GIN("poly_order",
          "no_pol_variation",
          "no_los_variation",
          "no_mblock_variation"),
      GIN_TYPE("Index", "Index", "Index", "Index"),
      GIN_DEFAULT(NODEF, "0", "0", "0"),
      GIN_DESC("Polynomial order to use for the fit.",
               "Set to 1 if the baseline off-set is the same for all "
               "Stokes components.",
               "Set to 1 if the baseline off-set is the same for all "
               "line-of-sights (inside each measurement block).",
               "Set to 1 if the baseline off-set is the same for all "
               "measurement blocks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddScatSpecies"),
      DESCRIPTION(R"--(Includes a scattering species in the Jacobian.

For 1D or 2D calculations the latitude and/or longitude grid of
the retrieval field should set to have zero length.

The number of elements added to the state vector (*x*) is::

  n_g1 * n_g2 * n_g3

where n_g1, n_g2 and n_g3 are the length of GIN ``g1``, ``g2`` and ``g3``,
respectively. Here empty vectors should be considered to have a length 1.
The elements are sorted with pressure as innermost loop, followed by
latitude and longitude as outermost loop.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda"),
      GIN("g1", "g2", "g3", "species", "quantity"),
      GIN_TYPE("Vector", "Vector", "Vector", "String", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, NODEF),
      GIN_DESC(
          "Pressure retrieval grid.",
          "Latitude retrieval grid.",
          "Longitude retreival grid.",
          "Name of scattering species, must match one element in *scat_species*.",
          "Retrieval quantity, e.g. \"IWC\"."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddSinefit"),
      DESCRIPTION(R"--(Includes sinusoidal baseline fit in the Jacobian.

Works as *jacobianAddPolyfit*, beside that a series of sine and
cosine terms are used for the baseline fit.

For each value in ``period_lengths`` one sine and one cosine term are
included (in mentioned order). By these two terms the amplitude and
\"phase\" for each period length can be determined. The sine and
cosine terms have value 0 and 1, respectively, for first frequency.

If the simulation/retrieval deals with a single spectrum, the number
of elements added to the state vector (*x*) is 2 * nperiods, where
nperiods is the length of ``period_lengths``. The same is true
if ``no_pol_variation``, ``no_los_variation`` and ``no_mblock_variation``
all are set to 1, even if several spectra are involved. Otherwise the
number of elements added to *x* depends on the number of spectra and
the settings of ``no_pol_variation``, ``no_los_variation`` and 
``no_mblock_variation``. The sine and cosine terms for each period
length are treated as a  separate retrieval quantities. That is, the
the elements associated with the first period length are grouped and
form together a retrieval quantity, etc. Inside each retrieval quantity
the pairs of sine and cosine terms are kept together, in given order.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid",
         "sensor_pos"),
      GIN("period_lengths",
          "no_pol_variation",
          "no_los_variation",
          "no_mblock_variation"),
      GIN_TYPE("Vector", "Index", "Index", "Index"),
      GIN_DEFAULT(NODEF, "0", "0", "0"),
      GIN_DESC("Period lengths of the fit.",
               "Set to 1 if the baseline off-set is the same for all "
               "Stokes components.",
               "Set to 1 if the baseline off-set is the same for all "
               "line-of-sights (inside each measurement block).",
               "Set to 1 if the baseline off-set is the same for all "
               "measurement blocks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddSpecialSpecies"),
      DESCRIPTION(R"--(Includes a special absorption species in the Jacobian.

Similar to *jacobianAddAbsSpecies* but only for number densities.

Species allowed are:

* \"electrons\"
* \"particulates\"

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
)--"),
      AUTHORS("Richard Larsson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda"),
      GIN("g1", "g2", "g3", "species"),
      GIN_TYPE("Vector", "Vector", "Vector", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "The species of the retrieval quantity."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddSurfaceQuantity"),
      DESCRIPTION(R"--(Includes a surface quantity in the Jacobian.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda"),
      GIN("g1", "g2", "quantity"),
      GIN_TYPE("Vector", "Vector", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Retrieval quantity, e.g. \"Wind speed\"."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddTemperature"),
      DESCRIPTION(R"--(Includes atmospheric temperatures in the Jacobian.

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
)--"),
      AUTHORS("Mattias Ekstrom", "Patrick Eriksson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda"),
      GIN("g1", "g2", "g3", "hse"),
      GIN_TYPE("Vector", "Vector", "Vector", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "on"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Flag to assume HSE or not (\"on\" or \"off\").")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAddWind"),
      DESCRIPTION(R"--(Includes one atmospheric wind component in the Jacobian.

The method follows the pattern of other Jacobian methods. The
calculations can only be performed by analytic expressions.
Some lower level function depends on frequency perturbations,
however, so therefore a frequency perturbation df is required
and as a consequence *abs_f_interp_order* must be > 0.

The wind field components are retrieved separately, and,
hence, the argument ``component`` can be \"u\", \"v\" or \"w\" 
for vector components, or just \"strength\" for total wind speed.

The number of elements added to the state vector (*x*) is::

  n_g1 * n_g2 * n_g3

where n_g1, n_g2 and n_g3 are the length of GIN ``g1``, ``g2`` and ``g3``,
respectively. Here empty vectors should be considered to have a length 1.
The elements are sorted with pressure as innermost loop, followed by
latitude and longitude as outermost loop.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda"),
      GIN("g1", "g2", "g3", "component", "dfrequency"),
      GIN_TYPE("Vector", "Vector", "Vector", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "v", "0.1"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retrieval grid.",
               "Wind component to retrieve",
               "This is the frequency perturbation")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianAdjustAndTransform"),
      DESCRIPTION(R"--(Applies adjustments and transformations on *jacobian*.

The method handles two tasks:
1. The retrieval transformations set by the user can not be applied
onthe  Jacobian inside *yCalc*. Transformations are instead applied
by calling this method.
2. It applies required adjustments of the Jacoboan. So far there is
only one possible adjustment. If any absorption species uses the \"rel\"
unit, an adjustment is needed for later iterations of the inversion.

If no tranformations are selected and the \"rel\" option is not used at
all, there is no need to call this method(, but you can still include it
without causing any error, the calculations will just be a bit slower).
Otherwise, this method should be called, typically as part of
*inversion_iterate_agenda*.

The method accepts if *jacobian* is empty, and then does, nothing.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian", "jacobian_quantities", "x"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianCalcDoNothing"),
      DESCRIPTION(R"--(This function doesn't do anything. It just exists to satisfy
the input and output requirement of the *jacobian_agenda*.

This method is added to *jacobian_agenda* by *jacobianAddAbsSpecies*
and some similar methods, and it should normally not be called by
the user.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian", "mblock_index", "iyb", "yb"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianCalcFreqShift"),
      DESCRIPTION(R"--(Calculates frequency shift jacobians by interpolation
of *iyb*.

This function is added to *jacobian_agenda* by jacobianAddFreqShift
and should normally not be called by the user.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian",
         "mblock_index",
         "iyb",
         "yb",
         
         "f_grid",
         "mblock_dlos",
         "sensor_response",
         "jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianCalcFreqStretch"),
      DESCRIPTION(R"--(Calculates frequency stretch jacobians by interpolation
of *iyb*.

This function is added to *jacobian_agenda* by jacobianAddFreqStretch
and should normally not be called by the user.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian",
         "mblock_index",
         "iyb",
         "yb",
         
         "f_grid",
         "mblock_dlos",
         "sensor_response",
         "sensor_response_pol_grid",
         "sensor_response_f_grid",
         "sensor_response_dlos_grid",
         "jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianCalcPointingZaInterp"),
      DESCRIPTION(R"--(Calculates zenith angle pointing deviation jacobians by
inter-extrapolation of *iyb*.

This function is added to *jacobian_agenda* by
jacobianAddPointingZa and should normally not be
called by the user.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian",
         "mblock_index",
         "iyb",
         "yb",
         
         "f_grid",
         "sensor_los",
         "mblock_dlos",
         "sensor_response",
         "sensor_time",
         "jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianCalcPointingZaRecalc"),
      DESCRIPTION(R"--(Calculates zenith angle pointing deviation jacobians by
recalulation of *iyb*.

This function is added to *jacobian_agenda* by
jacobianAddPointingZa and should normally not be
called by the user.
)--"),
      AUTHORS("Mattias Ekstrom", "Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian",
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
         "jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianCalcPolyfit"),
      DESCRIPTION(R"--(Calculates jacobians for polynomial baseline fit.

This function is added to *jacobian_agenda* by jacobianAddPolyfit
and should normally not be called by the user.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian",
         "mblock_index",
         "iyb",
         "yb",
         "sensor_response",
         "sensor_response_pol_grid",
         "sensor_response_f_grid",
         "sensor_response_dlos_grid",
         "jacobian_quantities"),
      GIN("poly_coeff"),
      GIN_TYPE("Index"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Polynomial coefficient to handle."),
      SETMETHOD(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianCalcSinefit"),
      DESCRIPTION(R"--(Calculates jacobians for sinusoidal baseline fit.

This function is added to *jacobian_agenda* by jacobianAddPolyfit
and should normally not be called by the user.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian",
         "mblock_index",
         "iyb",
         "yb",
         "sensor_response",
         "sensor_response_pol_grid",
         "sensor_response_f_grid",
         "sensor_response_dlos_grid",
         "jacobian_quantities"),
      GIN("period_index"),
      GIN_TYPE("Index"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Index among the period length specified for add-method."),
      SETMETHOD(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianClose"),
      DESCRIPTION(R"--(Closes the array of retrieval quantities and prepares for
calculation of the Jacobian matrix.

This function closes the *jacobian_quantities* array and sets
*jacobian_do* to 1.

Retrieval quantities should not be added after a call to this WSM.
No calculations are performed here.
)--"),
      AUTHORS("Mattias Ekstrom"),
      OUT("jacobian_do", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_agenda", "jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianFromTwoY"),
      DESCRIPTION(R"--(Sets *jacobian* based on the difference vetween two measurement vectors.

This function assumes that ``y_pert`` contains a measurement calculated
with some variable perturbed, in comparison to the calculation
behind *y*. The function takes the differences between ``y_pert``
and *y* to form a numerical derived estimate of *jacobian*.
This gives a Jacobian wit a single column.

*jacobian* equals here: (y_pert-y)/pert_size.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("y"),
      GIN("y_pert", "pert_size"),
      GIN_TYPE("Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Perturbed measurement vector",
               "Size of perturbation behind spectra in *ybatch*.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianFromYbatch"),
      DESCRIPTION(R"--(Sets *jacobian* based on perturbation calcuations.

This function assumes that *ybatch* contains spectra calculated
with some variable perturbed, in comparison to the calculation
behind *y*. The function takes the differences between *ybatch*
and *y* to form a numerical derived estimate of *jacobian*.

Column i of *jacobian* equals: (ybatch[i]-y)/pert_size.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ybatch", "y"),
      GIN("pert_size"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Size of perturbation behind spectra in *ybatch*.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianInit"),
      DESCRIPTION(R"--(Initialises the variables connected to the Jacobian matrix.

This function initialises the *jacobian_quantities* array so
that retrieval quantities can be added to it. Accordingly, it has
to be called before any calls to jacobianAddTemperature or
similar methods.

The Jacobian quantities are initialised to be empty.
)--"),
      AUTHORS("Mattias Ekstrom"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC(),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true),
      PASSWSVNAMES(false)));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianOff"),
      DESCRIPTION(R"--(Makes mandatory initialisation of some jacobian variables.

Some clear-sky jacobian WSVs must be initialised even if no such
calculations will be performed.  This is handled with this method.
That is, this method must be called when no clear-sky jacobians
will be calculated (even if cloudy-sky jacobians are calculated!).

Sets *jacobian_do* to 0.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian_do", "jacobian_agenda", "jacobian_quantities"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC(),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true),
      PASSWSVNAMES(false)));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianSetAffineTransformation"),
      DESCRIPTION(R"--(Adds an affine transformation of the last element of
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
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("jacobian_quantities"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities"),
      GIN("transformation_matrix", "offset_vector"),
      GIN_TYPE("Matrix", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The transformation matrix A", "The offset vector b")));

  md_data_raw.push_back(create_mdrecord(
      NAME("jacobianSetFuncTransformation"),
      DESCRIPTION(R"--(Sets the functional transformation of the last element of
*jacobian_quantities*.

See below for a general description of how retrieval transformations
are defined. Transformations are not applied by methods such as *yCalc*.
Instead, the method *jacobianAdjustAndTransform* must be called to
activate the transformations.

The following transformations can be selected (by ``transformation_func``):

- ``\"log\"``: The natural logarithm
- ``\"log10\"``: The base-10 logarithm
- ``\"atanh\"``: Area hyperbolic tangent 
- ``\"none\"``: No transformation at all

This method needs only to be called if a functional transformation
is wanted. Default is to make no such tranformation at all (i.e.
the option \"none\" exists only for reasons of flexibility).

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
)--"),
      AUTHORS("Patrick Eriksson", "Simon Pfreundschuh"),
      OUT("jacobian_quantities"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities"),
      GIN("transformation_func", "z_min", "z_max"),
      GIN_TYPE("String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "0", "-99e99"),
      GIN_DESC("The transformation function.",
               "Lower limit of z.",
               "Upper limit of z.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("LatLonFieldSet"),
      DESCRIPTION(R"--(Fills a latitude-longitude field with given input.

Grids and data must match in size.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("gfield2"),
      GOUT_TYPE("GriddedField2"),
      GOUT_DESC("Field to set."),
      IN(),
      GIN("latitude_grid", "longitude_grid", "data", "name"),
      GIN_TYPE("Vector", "Vector", "Matrix", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, ""),
      GIN_DESC("The latitude grid of ``data``.",
               "The longitude grid of ``data``.",
               "The data of the field (will become gfield2.data).",
               "The name of the field (will become gfield2.name).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("LatLonFieldSetConstant"),
      DESCRIPTION(R"--(Sets a latitude-longitude field to have a constant data value.

Both latitude and longitude grids are set to have length one,
with the value 0.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("gfield2"),
      GOUT_TYPE("GriddedField2"),
      GOUT_DESC("Field to set."),
      IN(),
      GIN("value", "name"),
      GIN_TYPE("Numeric", "String"),
      GIN_DEFAULT(NODEF, ""),
      GIN_DESC("The value (to place in gfield2.data).",
               "The name of the field (will become gfield2.name).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("lbl_checkedCalc"),
      DESCRIPTION(R"--(Checks that the line-by-line parameters are OK.

On failure, will throw.  On success, lbl_checked evals as true

Note that checks may become more stringent as ARTS evolves, especially for
\"new\" options.  This test might succeed in one version of ARTS but fail
in later versions
)--"),
      AUTHORS("Richard Larsson"),
      OUT("lbl_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species", "abs_species", "isotopologue_ratios"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("LocalTimeOffset"),
      DESCRIPTION(R"--(Sets the seconds between localtime and gmtime representation of now().
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("dt"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Time in seconds between local and gmt"),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("MatrixAdd"),
      DESCRIPTION(R"--(Adds a scalar to all elements of a matrix.

The result can either be stored in the same or another matrix.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Output Matrix"),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Matrix", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Matrix.", "The value to be added to the matrix.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("MatrixCBR"),
      DESCRIPTION(R"--(Sets a matrix to hold cosmic background radiation (CBR).

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("f"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Frequency vector.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("MatrixCopySparse"),
      DESCRIPTION(R"--(Creates a matrix by copying a variable of type Sparse.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Created (full) matrix."),
      IN(),
      GIN("input"),
      GIN_TYPE("Sparse"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The sparse matrix to be copied.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("MatrixDivide"),
      DESCRIPTION(R"--(Divides all elements of a matrix with the specified value.

The result can either be stored in the same or another
variable.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Output Matrix"),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Matrix", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Matrix.","Denominator.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("MatrixExtractFromTensor3"),
      DESCRIPTION(R"--(Extracts a Matrix from a Tensor3.

Copies page or row or column with given Index from input Tensor3
variable to output Matrix.
Higher order equivalent of *VectorExtractFromMatrix*.
)--"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Extracted matrix."),
      IN(),
      GIN("input", "i", "direction"),
      GIN_TYPE("Tensor3", "Index", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Input matrix.",
               "Index of page or row or column to extract.",
               "Direction. \"page\" or \"row\" or \"column\".")));

  md_data_raw.push_back(
      create_mdrecord(NAME("MatrixFromCovarianceMatrix"),
               DESCRIPTION(R"--(Turns a covariance matrix into a Matrix.
)--"),
               AUTHORS("Richard Larsson"),
               OUT(),
               GOUT("output"),
               GOUT_TYPE("Matrix"),
               GOUT_DESC("Dense Matrix."),
               IN(),
               GIN("input"),
               GIN_TYPE("CovarianceMatrix"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Input covariance matrix.")));


  md_data_raw.push_back(create_mdrecord(
      NAME("MatrixGaussian"),
      DESCRIPTION(R"--(Fills a matrix with a Gaussian function.

Works as *VectorGaussian* but grid, mean and si/fwhm must be
specified for each dimension.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("Y"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Output Matrix."),
      IN(),
      GIN("x_row", "x0_row", "si_row", "fwhm_row",
          "x_col", "x0_col", "si_col", "fwhm_col"),
      GIN_TYPE("Vector", "Numeric", "Numeric", "Numeric",
               "Vector", "Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "0", "-1", "-1",
                  NODEF, "0", "-1", "-1"),
      GIN_DESC("Grid of the function for row dimension.",
               "Centre/mean point of the function for row dimension.",
               "Row standard deviation of the function, ignored if <=0.",
               "Row full width at half-max of the function, ignored if <=0.",
               "Grid of the function for column dimension.",
               "Centre/mean point of the function for column dimension.",
               "Column standard deviation of the function, ignored if <=0.",
               "Column full width at half-max of the function, ignored if <=0.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("MatrixIdentity"),
      DESCRIPTION(R"--(Returns the identity matrix.

The size if the matrix created is n x n. Default is to return a
true identity matrix (I), but you can also select another value
along the diagonal by setting ``value``. That is, the output is
value * I.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Output matrix"),
      IN(),
      GIN("n", "value"),
      GIN_TYPE("Index", "Numeric"),
      GIN_DEFAULT(NODEF, "1"),
      GIN_DESC("Size of the matrix", "The value along the diagonal.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("MatrixMatrixMultiply"),
      DESCRIPTION(R"--(Multiply a Matrix with another Matrix and store the result in the result
Matrix.

This just computes the normal Matrix-Matrix product, Y = M * X. It is ok
if Y and X are the same Matrix.
)--"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT("Y"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("The result of the multiplication (dimension m x c)."),
      IN(),
      GIN("M", "X"),
      GIN_TYPE("Matrix", "Matrix"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The Matrix to multiply (dimension m x n).",
               "The original Matrix (dimension n x c).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("MatrixPlanck"),
      DESCRIPTION(R"--(Sets a matrix to hold blackbody radiation.

The radiation is assumed to be un-polarized and Stokes components
2-4 are zero. Number of Stokes components, that equals the number
of columns in the created matrix, is determined by ``stokes_dim``.
The number of rows in the created matrix equals the length of the
given frequency vector.

The standard definition, in ARTS, of the Planck function is
followed and the unit of the returned data is W/(m3 * Hz * sr).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("f", "t"),
      GIN_TYPE("Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Frequency vector.", "Temperature [K].")));

  md_data_raw.push_back(create_mdrecord(
      NAME("MatrixMultiply"),
      DESCRIPTION(R"--(Multiplies all elements of a matrix with the specified value.

The result can either be stored in the same or another
variable.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Output Matrix"),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Matrix", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Matrix.",
               "The value to be multiplied with the matrix.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("MatrixReshapeTensor3"),
      DESCRIPTION(R"--(Creates a matrix as reshaped version of a tenor3.

If the size of the tensor is [npages, nrows, ncols], the created
matrix gets size [npages * nrows, ncols]. The matrix is filled with
the tensor's page dimension as the outermost loop.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Matrix to fill."),
      IN(),
      GIN("input"),
      GIN_TYPE("Tensor3"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Tensor3 to copy.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("MatrixSetConstant"),
      DESCRIPTION(R"--(Creates a matrix and sets all elements to the specified value.

The size is determined by *ncols* and *nrows*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN("nrows", "ncols"),
      GIN("value"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Matrix value.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("MatrixSubtract"),
      DESCRIPTION(R"--(Subtracts a scalar from all elements of a matrix.

The result can either be stored in the same or another matrix.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Output Matrix"),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Matrix", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Matrix.", "The value to be subtracted from the matrix.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("MatrixUnitIntensity"),
      DESCRIPTION(R"--(Sets a matrix to hold unpolarised radiation with unit intensity.

Works as MatrixPlanck where the radiation is set to 1.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("f"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Frequency vector.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Matrix1ColFromVector"),
      DESCRIPTION(R"--(Forms a matrix containing one column from a vector.
)--"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("v"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The vector to be copied.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Matrix2ColFromVectors"),
      DESCRIPTION(R"--(Forms a matrix containing two columns from two vectors.

The vectors are included as columns in the matrix in the same order
as they are given.
)--"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("v1", "v2"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The vector to be copied into the first column.",
               "The vector to be copied into the second column.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Matrix3ColFromVectors"),
      DESCRIPTION(R"--(Forms a matrix containing three columns from three vectors.

The vectors are included as columns in the matrix in the same order
as they are given.
)--"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("v1", "v2", "v3"),
      GIN_TYPE("Vector", "Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("The vector to be copied into the first column.",
               "The vector to be copied into the second column.",
               "The vector to be copied into the third column.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Matrix1RowFromVector"),
      DESCRIPTION(R"--(Forms a matrix containing one row from a vector.
)--"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("v"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The vector to be copied.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Matrix2RowFromVectors"),
      DESCRIPTION(R"--(Forms a matrix containing two rows from two vectors.

The vectors are included as rows in the matrix in the same order
as they are given.
)--"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("v1", "v2"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The vector to be copied into the first row.",
               "The vector to be copied into the second row.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Matrix3RowFromVectors"),
      DESCRIPTION(R"--(Forms a matrix containing three rows from three vectors.

The vectors are included as rows in the matrix in the same order
as they are given.
)--"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("v1", "v2", "v3"),
      GIN_TYPE("Vector", "Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("The vector to be copied into the first row.",
               "The vector to be copied into the second row.",
               "The vector to be copied into the third row.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("mblock_dlosFrom1dAntenna"),
      DESCRIPTION(R"--(Sets *mblock_dlos* based on a 1D gaussian antenna response.

The length of *mblock_dlos* is determined by ``npoints``. The end
points of the grid are set to be the same as for the antenna
response. The spacing of the grid follows the magnitude of the
response; the spacing is smaller where the response is high.
More precisely, the grid points are determined by dividing the
cumulative sum of the response in equal steps.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("mblock_dlos"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("antenna_response"),
      GIN("npoints"),
      GIN_TYPE("Index"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Number of points (>1) to include in *mblock_dlos*.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("mc_antennaSetGaussian"),
      DESCRIPTION(R"--(Makes mc_antenna (used by MCGeneral) a 2D Gaussian pattern.

The gaussian antenna pattern is determined by ``za_sigma`` and
``aa_sigma``, which represent the standard deviations in the
uncorrelated bivariate normal distribution.
)--"),
      AUTHORS("Cory Davis"),
      OUT("mc_antenna"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("za_sigma", "aa_sigma"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Width in the zenith angle dimension as described above.",
               "Width in the azimuth angle dimension as described above.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("mc_antennaSetGaussianByFWHM"),
      DESCRIPTION(R"--(Makes mc_antenna (used by MCGeneral) a 2D Gaussian pattern.

The gaussian antenna pattern is determined by ``za_fwhm`` and
``aa_fwhm``, which represent the full width half maximum (FWHM)
of the antenna response, in the zenith and azimuthal planes.
)--"),
      AUTHORS("Cory Davis"),
      OUT("mc_antenna"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("za_fwhm", "aa_fwhm"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Width in the zenith angle dimension as described above.",
               "Width in the azimuth angle dimension as described above.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("mc_antennaSetPencilBeam"),
      DESCRIPTION(R"--(Makes mc_antenna (used by MCGeneral) a pencil beam.

This WSM makes the subsequent MCGeneral WSM perform pencil beam
RT calculations.
)--"),
      AUTHORS("Cory Davis"),
      OUT("mc_antenna"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("MCGeneral"),
      DESCRIPTION(R"--(A generalised 3D reversed Monte Carlo radiative algorithm, that
allows for 2D antenna patterns, surface reflection and arbitrary
sensor positions.

The main output variables *y* and *mc_error* represent the
Stokes vector integrated over the antenna function, and the
estimated error in this vector, respectively.

The WSV *mc_max_iter* describes the maximum number of \'photons\'
used in the simulation (more photons means smaller *mc_error*).
*mc_std_err* is the desired value of mc_error. *mc_max_time* is
the maximum allowed number of seconds for MCGeneral. The method
will terminate once any of the max_iter, std_err, max_time
criteria are met. If negative values are given for these
parameters then it is ignored.

The WSV *mc_min_iter* sets the minimum number of photons to apply
before the condition set by *mc_std_err* is considered. Values
of *mc_min_iter* below 100 are not accepted.

Only \"1\" and \"RJBT\" are allowed for *iy_unit*. The value of
*mc_error* follows the selection for *iy_unit* (both for in- and
output.
)--"),
      AUTHORS("Cory Davis"),
      OUT("y",
          "mc_iteration_count",
          "mc_error",
          "mc_points",
          "mc_source_domain",
          "mc_scat_order"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("mc_antenna",
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
         "mc_taustep_limit"),
      GIN("l_mc_scat_order", "t_interp_order"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT("11", "1"),
      GIN_DESC("The length to be given to *mc_scat_order*. Note that"
               " scattering orders equal and above this value will not"
               " be counted.",
               "Interpolation order of temperature for scattering data (so"
               " far only applied in phase matrix, not in extinction and"
               " absorption.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("MCRadar"),
      DESCRIPTION(R"--(A radar 3D foward Monte Carlo radiative algorithm, that allows 
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
order to consider, after which \'photon\'-tracing will be
terminated. A value of one calculates only single scattering.

The WSV *mc_max_iter* describes the maximum number of \'photons\'
used in the simulation (more photons means smaller *mc_error* ).
The method will terminate once the max_iter criterium is met.
If negative values are given for these parameters then it is
ignored.

Here \"1\" and \"Ze\" are the allowed options for *iy_unit_radar*.
The value of *mc_error* follows the selection for *iy_unit_radar*
(both for in- and output. See *yRadar* for details of the units.
)--"),
      AUTHORS("Ian S. Adams"),
      OUT("y", "mc_error"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("mc_antenna",
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
         "mc_max_iter"),
      GIN("ze_tref", "k2", "t_interp_order"),
      GIN_TYPE("Numeric", "Numeric", "Index"),
      GIN_DEFAULT("273.15", "-1", "1"),
      GIN_DESC("Reference temperature for conversion to Ze.",
               "Reference dielectric factor.",
               "Interpolation order of temperature for scattering data (so"
               " far only applied in phase matrix, not in extinction and"
               " absorption.")));

  md_data_raw.push_back(
      create_mdrecord(NAME("MCSetSeedFromTime"),
               DESCRIPTION(R"--(Sets the value of mc_seed from system time
)--"),
               AUTHORS("Cory Davis"),
               OUT("mc_seed"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fieldRescalePopulationLevels"),
      DESCRIPTION(R"--(Rescale NLTE field to expected total distribution amongst levels
)--"),
      AUTHORS("Richard Larsson"),
      OUT("atm_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field"),
      GIN("s"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Scaling (e.g., 0.75 for only orth-water on Earth)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("atm_fieldForSingleSpeciesNonOverlappingLines"),
      DESCRIPTION(R"--(NLTE field for a simple setup.

This will solve for ``nlte_field`` in the input atmosphere.
The solver depends on the lines not overlapping and that there
is only a single species in the atmosphere.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("atm_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field",
         "abs_species",
         "abs_lines_per_species",
         "collision_coefficients",
         "collision_line_identifiers",
         "isotopologue_ratios",
         "iy_main_agenda",
         "ppath_agenda",
         "iy_space_agenda",
         "iy_surface_agenda",
         "iy_cloudbox_agenda",
         "propmat_clearsky_agenda",
         "surface_field",
         "nlte_do"),
      GIN("df", "convergence_limit", "nz", "nf", "dampened", "iteration_limit"),
      GIN_TYPE("Numeric", "Numeric", "Index", "Index", "Index", "Index"),
      GIN_DEFAULT(NODEF, "1e-6", NODEF, NODEF, NODEF, "20"),
      GIN_DESC("relative frequency to line center",
               "max relative change in ratio of level to stop iterations",
               "number of zenith angles",
               "number of frequency grid-points per line",
               "use transmission dampening or not",
               "max number of iterations before defaul break of iterations")));

  md_data_raw.push_back(create_mdrecord(
      NAME("collision_coefficientsFromSplitFiles"),
      DESCRIPTION(R"--(Reads *collision_coefficients* and *collision_line_identifiers* from files.
The species in in these files must match *abs_species*.  The location
must also contain an *ArrayOfQuantumIdentifier* file ending with ``qid.xml``
)--"),
      AUTHORS("Richard Larsson"),
      OUT("collision_coefficients", "collision_line_identifiers"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("basename"),
      GIN_TYPE("String"),
      GIN_DEFAULT("./"),
      GIN_DESC("path to files to read")));

  md_data_raw.push_back(create_mdrecord(
      NAME("NumericAdd"),
      DESCRIPTION(R"--(Adds a Numeric and a value (output = input + value).

The result can either be stored in the same or another Numeric.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Output Numeric."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Numeric.", "Value to add.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("NumericClip"),
      DESCRIPTION(R"--(Clipping of a Numeric.

The input value is copied to the output one (that can be same WSV)
but ensures that ``out`` is inside the range [limit_low,limit_high].
When the input value is below ``limit_low``, ``output`` is set to ``limit_low``.
And the same is performed with respect to ``limit_high``.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Output Numeric."),
      IN(),
      GIN("input", "limit_low", "limit_high"),
      GIN_TYPE("Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "-Inf", "Inf"),
      GIN_DESC("Input Numeric.",
               "Lower limit for clipping.",
               "Upper limit for clipping.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("NumericDivide"),
      DESCRIPTION(R"--(Divides a Numeric with a value (output = input / value).

The result can either be stored in the same or another Numeric.
)--"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Output Numeric."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Numeric (numerator).", "Denominator.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("NumericFromVector"),
      DESCRIPTION(R"--(Derivs a Numeric from a vector, following selected operation.

The following operations can be selected:

- ``\"first\"``: Selects the first element of the vector.
- ``\"last\"``: Selects the last element of the vector.
- ``\"max\"``: Selects the maximum element of the vector.
- ``\"min\"``: Selects the minimum element of the vector.
- ``\"mean\"``: Calculates the mean of the vector.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Output Numeric."),
      IN(),
      GIN("input", "op"),
      GIN_TYPE("Vector", "String"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input vector.", "Selected operation.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("NumericInterpAltLatLonField"),
      DESCRIPTION(R"--(Interpolates an altitude-latitude-longitiude field.

The gridded field must have \"Altitude\", \"Latitude\" and
\"Longitude\" as dimensions.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("value"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Result of interpolation"),
      IN(),
      GIN("gfield3", "pos"),
      GIN_TYPE("GriddedField3", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Gridded field to be interpolated.", 
               "Interpolate to this position [z, lat, lon].")));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("NumericInterpLatLonField"),
      DESCRIPTION(R"--(Interpolates a latitude-longitiude field to the selected position.

The gridded field must have \"Latitude\" and \"Longitude\" as dimensions.

The position shall be given as a full atmospheric position. The altitude
in ``pos`` is ignored.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("value"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Result of interpolation"),
      IN(),
      GIN("gfield2", "pos"),
      GIN_TYPE("GriddedField2", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Gridded field to be interpolated.",
               "Interpolate to this position [z, lat, lon].")));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("NumericInterpVector"),
      DESCRIPTION(R"--(Interpolates a vector to the selected position.

Returns y(xv) given y(x).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("value"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Result of interpolation"),
      IN(),
      GIN("x", "y", "xv"),
      GIN_TYPE("Vector", "Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Positions (grid) where *y* given.",
               "Values of function to interpolate.",
               "Interpolate to this value.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("NumericMultiply"),
      DESCRIPTION(R"--(Multiplies a Numeric with a value (output = input * value).

The result can either be stored in the same or another Numeric.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Output Numeric."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Numeric.", "Multiplier.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("NumericSubtract"),
      DESCRIPTION(R"--(Subtracts a Numeric value (output = input - value).

The result can either be stored in the same or another Numeric.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Output Numeric."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Numeric.", "Subtrahend.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("nelemGet"),
      DESCRIPTION(R"--(Retrieve nelem from given variable and store the value in the
variable *nelem*.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("nelem"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("v"),
      GIN_TYPE((ARRAY_GROUPS + ", Vector").c_str()),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Variable to get the number of elements from."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("ncolsGet"),
      DESCRIPTION(R"--(Retrieve ncols from given variable and store the value in the
workspace variable *ncols*
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("ncols"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("v"),
      GIN_TYPE("Matrix, Sparse, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Variable to get the number of columns from."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("nrowsGet"),
      DESCRIPTION(R"--(Retrieve nrows from given variable and store the value in the
workspace variable *nrows*
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("nrows"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("v"),
      GIN_TYPE("Matrix, Sparse, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Variable to get the number of rows from."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("npagesGet"),
      DESCRIPTION(R"--(Retrieve npages from given variable and store the value in the
workspace variable *npages*
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("npages"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("v"),
      GIN_TYPE("Tensor3, Tensor4, Tensor5, Tensor6, Tensor7"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Variable to get the number of pages from."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("nbooksGet"),
      DESCRIPTION(R"--(Retrieve nbooks from given variable and store the value in the
workspace variable *nbooks*
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("nbooks"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("v"),
      GIN_TYPE("Tensor4, Tensor5, Tensor6, Tensor7"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Variable to get the number of books from."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("nshelvesGet"),
      DESCRIPTION(R"--(Retrieve nshelves from given variable and store the value in the
workspace variable *nshelves*
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("nshelves"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("v"),
      GIN_TYPE("Tensor5, Tensor6, Tensor7"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Variable to get the number of shelves from."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("nvitrinesGet"),
      DESCRIPTION(R"--(Retrieve nvitrines from given variable and store the value in the
workspace variable *nvitrines*
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("nvitrines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("v"),
      GIN_TYPE("Tensor6, Tensor7"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Variable to get the number of vitrines from."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("nlibrariesGet"),
      DESCRIPTION(R"--(Retrieve nlibraries from given variable and store the value in the
workspace variable *nlibraries*
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("nlibraries"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("v"),
      GIN_TYPE("Tensor7"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Variable to get the number of libraries from."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("nlteOff"),
               DESCRIPTION(R"--(Disable Non-LTE calculations.
)--"),
               AUTHORS("Oliver Lemke"),
               OUT("nlte_do"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesPopulationNlteField"),
      DESCRIPTION(R"--(Turns on NTLE calculations.

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

  \"CO2-626 EN v1 0/1 v2 1/1 l2 1/1 v3 0/1 r 1/1\"

and the matching will match ALL lines with the above.  Note then that if, e.g.,
the \"v1 0/1\" term was removed from the above, then ARTS will assume that
\"v1\" is not part of the level of energy state of interest, so lines
of different \"v1\" will be matched as the same state.  If a line is matched
to more than one energy state, errors should be thrown, but be careful.

Set type of population to change computations and expected input as:

- ``\"LTE\"``: Compute population by ratios found from LTE temperatures
- ``\"TV\"``: Compute population by ratios found from NLTE vibrational temperatures
- ``\"ND\"``: Compute population by ratios found from NLTE number densities
)--"),
      AUTHORS("Richard Larsson"),
      OUT("nlte_do", "abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species", "atm_field", "nlte_vib_energies"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("ArrayOfQuantumIdentifierFromLines"),
      DESCRIPTION(R"--(Sets an ArrayOfQuantumIdentifier to all levels in *abs_lines_per_species*
with defined quantum numbers

Lines without defined quantum numbers are ignored
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("ArrayOfQuantumIdentifier"),
      GOUT_DESC("Identifiers to all levels in *abs_lines_per_species*"),
      IN("abs_lines_per_species"),
      GIN("global"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Only look at global quantum numbers")));

  md_data_raw.push_back(create_mdrecord(
      NAME("timeNow"),
      DESCRIPTION(R"--(Sets time to system_clock::now().
)--"),
      AUTHORS("Richard Larsson"),
      OUT("time"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("timeOffset"),
      DESCRIPTION(R"--(Offsets time for some seconds
)--"),
      AUTHORS("Richard Larsson"),
      OUT("time"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("time"),
      GIN("offset"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Time in seconds")));

  md_data_raw.push_back(create_mdrecord(
      NAME("OEM"),
      DESCRIPTION(R"--(Inversion by the so called optimal estimation method (OEM).

Work in progress ...

The cost function to minimise, including a normalisation with length
of *y*, is::

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

  - ``\"li\"``: A linear problem is assumed and a single iteration is performed.
  - ``\"li_cg\"``: A linear problem is assumed and solved using the CG solver.
  - ``\"gn\"``: Non-linear, with Gauss-Newton iteration scheme.
  - ``\"gn_cg\"``: Non-linear, with Gauss-Newton and conjugate gradient solver.
  - ``\"lm\"``: Non-linear, with Levenberg-Marquardt (LM) iteration scheme.
  - ``\"lm_cg\"``: Non-linear, with Levenberg-Marquardt (LM) iteration scheme and conjugate gradient solver.

- ``max_start_cost``:
  No inversion is done if the cost matching the a priori state is above
  this value. If set to a negative value, all values are accepted.
  This argument also controls if the start cost is calculated. If
  set to <= 0, the start cost in *oem_diagnostics* is set to NaN
  when using \"li\" and \"gn\".
- ``x_norm``:
  A normalisation vector for *x*. A normalisation of *x* can be needed
  due to limited numerical precision. If this vector is set to be empty
  no normalisation is done (defualt case). Otherwise, this must be a
  vector with same length as *x*, just having values above zero.
  Elementwise division between *x* and ``x_norm`` (x./x_norm) shall give
  a vector where all values are in the order of unity. Maybe the best
  way to set ``x_norm`` is x_norm = sqrt( diag( Sx ) ).
- ``max_iter``:
  Maximum number of iterations to perform. No effect for \"li\".
- ``stop_dx``:
  Iteration stop criterion. The criterion used is the same as given
  in Rodgers\' \"Inverse Methods for Atmospheric Sounding\"
- ``lm_ga_settings``:
  Settings controlling the gamma factor, part of the \"LM\" method.
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

  The default setting triggers an error if \"lm\" is selected.
- ``clear matrices``:
   With this flag set to 1, *jacobian* and *dxdy* are returned as empty
   matrices.
- ``display_progress``:
   Controls if there is any screen output. The overall report level
   is ignored by this WSM.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("x",
          "yf",
          "jacobian",
          "dxdy",
          "oem_diagnostics",
          "lm_ga_history",
          "oem_errors"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("xa",
         "x",
         "covmat_sx",
         "yf",
         "y",
         "covmat_se",
         "jacobian",
         "jacobian_quantities",
         "inversion_iterate_agenda"),
      GIN("method",
          "max_start_cost",
          "x_norm",
          "max_iter",
          "stop_dx",
          "lm_ga_settings",
          "clear_matrices",
          "display_progress"),
      GIN_TYPE("String",
               "Numeric",
               "Vector",
               "Index",
               "Numeric",
               "Vector",
               "Index",
               "Index"),
      GIN_DEFAULT(NODEF, "Inf", "[]", "10", "0.01", "[]", "0", "0"),
      GIN_DESC("Iteration method. For this and all options below, see "
               "further above.",
               "Maximum allowed value of cost function at start.",
               "Normalisation of Sx.",
               "Maximum number of iterations.",
               "Stop criterion for iterative inversions.",
               "Settings associated with the ga factor of the LM method.",
               "An option to save memory.",
               "Flag to control if inversion diagnostics shall be printed "
               "on the screen.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("avkCalc"),
      DESCRIPTION(R"--(Calculate the averaging kernel matrix.
This is done by describing the sensitivity of the
OEM retrieval with respect to the true state of the system. A prerequisite
for the calculation of the averaging kernel matrix is a successful OEM
calculation in which the *jacobian* and the gain matrix *dxdy* have been calculated.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("avk"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("dxdy", "jacobian"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("covmat_soCalc"),
      DESCRIPTION(R"--(Calculates the covariance matrix describing the error due to uncertainties
in the observation system.
The uncertainties of the observation system are
described by *covmat_se*, which must be set by the user to include the
relevant contributions from the measurement and the forward model.

Prerequisite for the calculation of *covmat_so* is a successful OEM
computation where also the gain matrix has been computed.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_so"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("dxdy", "covmat_se"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("covmat_ssCalc"),
      DESCRIPTION(R"--(Calculates the covariance matrix describing the error due to smoothing.

The calculation of *covmat_ss* also requires the averaging kernel matrix *avk*
to be computed after a successful OEM calculation.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_ss"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("avk", "covmat_sx"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("opt_prop_bulkCalc"),
      DESCRIPTION(R"--(Calculates bulk absorption extinction at one atmospheric grid point.

This WSM sums up the monochromatic absorption vectors and
extinction matrices of all scattering elements (*abs_vec_spt* and
*ext_mat_spt*, respectively) weighted by their respective
particle number density given by *pnd_field*, for a single location
within the cloudbox, given by *scat_p_index*, *scat_lat_index*, and
*scat_lon_index*.
The resulting  extinction matrix is added to the workspace variable
*ext_mat*.
)--"),
      AUTHORS("Jana Mendrok, Sreerekha T.R."),
      OUT("ext_mat", "abs_vec"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ext_mat",
         "abs_vec",
         "ext_mat_spt",
         "abs_vec_spt",
         "pnd_field",
         "scat_p_index",
         "scat_lat_index",
         "scat_lon_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("opt_prop_sptFromData"),
      DESCRIPTION(R"--(Calculates monochromatic optical properties for all scattering
elements.

In this function the extinction matrix and the absorption vector
are calculated in the laboratory frame. An interpolation of the
data on the actual frequency is the first step in this function.
The next step is a transformation from the database coordinate
system to the laboratory coordinate system.

Output of the function are *ext_mat_spt* and *abs_vec_spt*, which
hold the optical properties for a specified propagation direction
for each scattering element.
)--"),
      AUTHORS("Claudia Emde"),
      OUT("ext_mat_spt", "abs_vec_spt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ext_mat_spt",
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
         "scat_lon_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("opt_prop_sptFromScat_data"),
      DESCRIPTION(R"--(Derives monochromatic optical properties for all scattering
elements.

As *opt_prop_sptFromData*, but using frequency pre-interpolated
data (as produced by *scat_dataCalc*), i.e. in here no frequency
interpolation is done anymore.
)--"),
      AUTHORS("Jana Mendrok, Claudia Emde"),
      OUT("ext_mat_spt", "abs_vec_spt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ext_mat_spt",
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
         "scat_lon_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("opt_prop_sptFromMonoData"),
      DESCRIPTION(R"--(Calculates optical properties for the scattering elements.

As *opt_prop_sptFromData* but no frequency interpolation is
performed. The single scattering data is here obtained from
*scat_data_mono*, instead of *scat_data*.
)--"),
      AUTHORS("Cory Davis"),
      OUT("ext_mat_spt", "abs_vec_spt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ext_mat_spt",
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
         "scat_lon_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord(NAME("output_file_formatSetAscii"),
               DESCRIPTION(R"--(Sets the output file format to ASCII.
)--"),
               AUTHORS("Oliver Lemke"),
               OUT("output_file_format"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord(NAME("output_file_formatSetBinary"),
               DESCRIPTION(R"--(Sets the output file format to binary.
)--"),
               AUTHORS("Oliver Lemke"),
               OUT("output_file_format"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord(NAME("output_file_formatSetZippedAscii"),
               DESCRIPTION(R"--(Sets the output file format to zipped ASCII.
)--"),
               AUTHORS("Oliver Lemke"),
               OUT("output_file_format"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("particle_bulkpropRadarOnionPeeling"),
      DESCRIPTION(R"--(Inverts radar reflectivities by in an onion peeling manner.

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
two reasons to ignore it. It can cause a \"run away\" effect in the
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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("atm_field", "particle_bulkprop_names"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         "atm_field",
         "surface_field",
         "atmfields_checked",
         "atmgeom_checked",
         "f_grid",
         "propmat_clearsky_agenda",
         "scat_species"),
      GIN("invtable",
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
          "atten_hyd_max"),
      GIN_TYPE("ArrayOfGriddedField3", "Matrix", "Tensor3", "Numeric",
               "Matrix", "Index", "Numeric", "Numeric", "Numeric",
               "Index", "Index", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "-99", NODEF, "0", "273.15",
                  "10e-3", "5e-3", "1", "1", "0.5", "3"),
      GIN_DESC("Inversion table, see above.",
               "Incidence angles.",
               "Field of radar reflectivities, in dBZe.",
               "Noise level. See above.",
               "Height of clutter zone. Either same size as ``z_surface`` or a single "
               "value. In the later case, that value is applied at all positions.",
               "Flag to fill clutter zone, by copying retrieval just above it.",
               "Phase boundary temperature. See above.",
               "Max reasonable water content",
               "Clip value for water content retrievals.",
               "Flag to consider attenuation due to hydrometeors.",
               "Flag to consider attenuation due to absorption species.",
               "Hydrometeor attenuation scaling factor.",
               "Hydrometeor attenuation not allowed to pass this value [dB].")));

  md_data_raw.push_back(create_mdrecord(
      NAME("particle_bulkprop_fieldClip"),
      DESCRIPTION(R"--(Clipping of ``particle_bulkprop_field``.

The method allows you to apply hard limits the values of
``particle_bulkprop_field``. All values, of the property selected,
below ``limit_low``, are simply set to ``limit_low``. And the same
is performed with respect to ``limit_high``. That is, the data in x
for the retrieval quantity are forced to be inside the range
[limit_low,limit_high].

Setting species=\"ALL\", is a shortcut for applying the limits on all
properties.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("atm_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field"),
      GIN("bulkprop_name", "limit_low", "limit_high"),
      GIN_TYPE("String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "-Inf", "Inf"),
      GIN_DESC("Name of bulk property to consider, or \"ALL\".",
               "Lower limit for clipping.",
               "Upper limit for clipping.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("particle_massesFromMetaDataSingleCategory"),
      DESCRIPTION(R"--(Sets *particle_masses* based on *scat_meta* assuming
all particles are of the same mass category.

This method derives the particle masses from the mass entry
of each scattering element. It is assumed that all scattering
elements represent particles of the same (bulk) matter
(e.g. water or ice). With other words, a single mass category
is assumed (see *particle_masses* for a definition of \"mass
category\").

If just having clouds, the resulting mass category can be seen as
the total cloud water content, with possible contribution from
both ice and liquid phase.
)--"),
      AUTHORS("Jana Mendrok", "Patrick Eriksson"),
      OUT("particle_masses"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_meta"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("particle_massesFromMetaData"),
      DESCRIPTION(R"--(Derives *particle_masses* from *scat_meta*.

It extracts the mass information of the scattering elements
from *scat_meta*. Each scattering species is taken as a
separate category of particle_masses, i.e., the resulting
*particle_masses* matrix will contain as many columns as
scattering species are present in *scat_meta*.
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("particle_masses"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_meta"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("pha_matCalc"),
      DESCRIPTION(R"--(Calculates the total phase matrix of all scattering elements.

This function sums up the monochromatic phase matrices of all
scattering elements *pha_mat_spt* weighted with  their respective
particle number density, given by *pnd_field*, for a single location
within the cloudbox, given by *scat_p_index*, *scat_lat_index*, and
*scat_lon_index*.
)--"),
      AUTHORS("Sreerekha T.R."),
      OUT("pha_mat"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pha_mat_spt",
         "pnd_field",
         "scat_p_index",
         "scat_lat_index",
         "scat_lon_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("pha_mat_sptFromData"),
      DESCRIPTION(R"--(Calculation of the phase matrix of the individual scattering elements.

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
)--"),
      AUTHORS("Claudia Emde"),
      OUT("pha_mat_spt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pha_mat_spt",
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
         "scat_lon_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("pha_mat_sptFromDataDOITOpt"),
      DESCRIPTION(R"--(Calculation of the phase matrix of the individual scattering elements.

In this function the phase matrix is extracted from
*pha_mat_sptDOITOpt*. It can be used in the agenda
*pha_mat_spt_agenda*. This method must be used in combination with
*DoitScatteringDataPrepare*.

Temperature is considered as described for *pha_mat_sptFromData*
)--"),
      AUTHORS("Claudia Emde"),
      OUT("pha_mat_spt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pha_mat_spt",
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
         "scat_lon_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("pha_mat_sptFromMonoData"),
      DESCRIPTION(R"--(Calculation of the phase matrix of the individual scattering elements.

This function is the monochromatic version of *pha_mat_sptFromData*.
)--"),
      AUTHORS("Claudia Emde"),
      OUT("pha_mat_spt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pha_mat_spt",
         "scat_data_mono",
         "doit_za_grid_size",
         "aa_grid",
         "za_index",
         "aa_index",
         "rtp_temperature",
         "pnd_field",
         "scat_p_index",
         "scat_lat_index",
         "scat_lon_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("pha_mat_sptFromScat_data"),
      DESCRIPTION(R"--(Calculation of the phase matrix of the individual scattering elements.

As *pha_mat_sptFromData*, but using frequency pre-interpolated
data (as produced by *scat_dataCalc*), i.e. in here no frequency
interpolation is done anymore.
)--"),
      AUTHORS("Jana Mendrok, Claudia Emde"),
      OUT("pha_mat_spt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pha_mat_spt",
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
         "scat_lon_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("pndFromPsd"),
      DESCRIPTION(R"--(Calculates *pnd_data* from given *psd_data* for one scattering species.

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
)--"),
      AUTHORS("Jana Mendrok, Patrick Eriksson"),
      OUT("pnd_data", "dpnd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pnd_size_grid",
         "psd_data",
         "psd_size_grid",
         "dpsd_data_dx",
         "scat_data",
         "f_grid",
         "scat_data_checked"),
      GIN("quad_order",
          "scat_index",
          "threshold_se_ext",
          "threshold_ss_ext",
          "threshold_se_pnd"),
      GIN_TYPE("Index", "Index", "Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT("1", NODEF, "0.02", "1e-8", "0.02"),
      GIN_DESC("Order of bin quadrature.",
               "Take data from scattering species of this index (0-based) in"
               " *scat_data*.",
               "Maximum allowed extinction fraction in each of the edge size"
               " bins.",
               "Minimum bulk extinction in the processed scattering species"
               " for which to apply size grid representation checks.",
               "Minimum ratio of edge point pnd to maximum pnd of this"
               " scattering element over all pressure levels.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("pndFromPsdBasic"),
      DESCRIPTION(R"--(Calculates *pnd_data* from given *psd_data*.

As *pndFromPsdBasic*, but without bulk extinction representation
checks.
)--"),
      AUTHORS("Jana Mendrok, Patrick Eriksson"),
      OUT("pnd_data", "dpnd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pnd_size_grid", "psd_data", "psd_size_grid", "dpsd_data_dx"),
      GIN("quad_order"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Order of bin quadrature.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("pnd_fieldCalcFromParticleBulkProps"),
      DESCRIPTION(R"--(Converts particle bulk property data to *pnd_field*.

In short, the method combines *scat_species*, *pnd_agenda_array*,
``particle_bulkprop_field`` and their associated variables to derive
*pnd_field*.

The method does nothing if cloudbox is inactive.

Otherwise, cloudbox limits must be set before calling the method,
and ``particle_bulkprop_field`` is checked to have non-zero elements
just inside the cloudbox.
)--"),
      AUTHORS("Patrick Eriksson, Jana Mendrok"),
      OUT("pnd_field", "dpnd_field_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         "cloudbox_on",
         "cloudbox_limits",
         "scat_species",
         "scat_data",
         "scat_meta",
         "atm_field",
         "particle_bulkprop_names",
         "pnd_agenda_array",
         "pnd_agenda_array_input_names",
         "jacobian_do",
         "jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("pnd_fieldCalcFrompnd_field_raw"),
      DESCRIPTION(R"--(Interpolation of particle number density fields to calculation grid
inside cloudbox.

This method interpolates the particle number density field
from the raw data *pnd_field_raw* to obtain *pnd_field*.
For 1D cases, where internally *GriddedFieldPRegrid* and
*GriddedFieldLatLonRegrid* are applied, ``zeropadding`` = 1 sets the
*pnd_field* at pressure levels levels exceeding pnd_field_raw's
pressure grid to 0 (not implemented for 2D and 3D yet). Default:
zeropadding = 0, which throws an error if the calculation pressure grid
``p_grid`` is not completely covered by pnd_field_raw's pressure grid.
)--"),
      AUTHORS("Sreerekha T.R.", "Claudia Emde", "Oliver Lemke"),
      OUT("pnd_field", "dpnd_field_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("p_grid",
         "lat_grid",
         "lon_grid",
         "pnd_field_raw",
         "cloudbox_limits",
         "jacobian_quantities"),
      GIN("zeropadding"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("Allow zeropadding of pnd_field.")));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("pnd_fieldExpand1D"),
      DESCRIPTION(R"--(Maps a 1D pnd_field to a (homogeneous) 2D or 3D pnd_field.

This method takes a 1D *pnd_field* and converts it to a 2D or 3D
\"cloud\". It is assumed that a complete 1D case has been created,
and after this ``atmosphere_dim``, ``lat_grid``, ``lon_grid`` and
*cloudbox_limits* have been changed to a 2D or 3D case (without
changing the vertical extent of the cloudbox.

No modification of *pnd_field* is made for the pressure dimension.
At the latitude and longitude cloudbox edge points *pnd_field* is set to
zero. This corresponds to nzero=1. If you want a larger margin between
the lat and lon cloudbox edges and the \"cloud\" you increase
``nzero``, where ``nzero`` is the number of grid points for which
*pnd_field* shall be set to 0, counted from each lat and lon edge.

See further ``AtmFieldsExpand1D``.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("pnd_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pnd_field", "cloudbox_on", "cloudbox_limits"),
      GIN("nzero"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Number of zero values inside lat and lon limits.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("pnd_fieldZero"),
      DESCRIPTION(R"--(Sets *pnd_field* to zero.

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
)--"),
      AUTHORS("Claudia Emde, Jana Mendrok"),
      OUT("pnd_field", "dpnd_field_dx", "scat_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data",
         "f_grid",
         "cloudbox_limits",
         "jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("ppathFromRtePos2"),
      DESCRIPTION(R"--(Determines the propagation path from *rte_pos2* to *rte_pos*.

The propagation path linking *rte_pos* and *rte_pos2* is calculated
and returned. The method determines the path in a pure numerical
manner, where a simple algorithm is applied. The task is to find
the value of *rte_los* (at *rte_pos*) linking the two positions.

See the user guide for a description of the search algorithm,
including a more detailed definition of ``za_accuracy``, 
``pplrt_factor`` and ``pplrt_lowest``.

The standard application of this method should be to radio link
calculations, where *rte_pos2* corresponds to a transmitter, and
*rte_pos* to the receiver/sensor.

The details of the ray tracing is controlled by *ppath_step_agenda*
as usual.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppath", "rte_los", "ppath_lraytrace"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath_step_agenda",
         "atm_field",
         "f_grid",
         "surface_field",
         "rte_pos",
         "rte_pos2",
         "rte_los",
         "ppath_lmax",
         "ppath_lraytrace"),
      GIN("za_accuracy", "pplrt_factor", "pplrt_lowest"),
      GIN_TYPE("Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT("2e-5", "5", "0.5"),
      GIN_DESC("Required accuracy, in form of the maximum allowed angular "
               "off-set [deg].",
               "The factor with which ppath_lraytrace is decreased if "
               "no solution is found.",
               "Lowest value ppath_lraytrace to consider. The calculations "
               "are halted if this length is passed.")));
*/
/*
  md_data_raw.push_back(create_mdrecord(
      NAME("ppathPlaneParallel"),
      DESCRIPTION(R"--(Propagation path calculations for a plane parallel atmosphere.

This method basically assumes that the planet's radius is infinite,
i.e. the planet surface has no curvature. Some consequences of this
assumption:

- the mathod can only be used for 1D
- zenith angles between 89.9 and 90.1 deg are not allowed
- refraction is always neglected
- radii in ppath are set to Inf

Notice that the method provides full propagation paths. This means
that *ppath_step_agenda* is ignored (and thus also refraction).
On the other hand, the method considers the cloudbox exactly as
the standard path calculations.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppath"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         "atm_field",
         "surface_field",
         "cloudbox_on",
         "cloudbox_limits",
         "ppath_inside_cloudbox_do",
         "rte_pos",
         "rte_los",
         "ppath_lmax"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
*/
/*
  md_data_raw.push_back(create_mdrecord(
      NAME("ppathStepByStep"),
      DESCRIPTION(R"--(Standard method for calculation of propagation paths.

This method calculates complete propagation paths in a stepwise
manner. Each step is denoted as a \"ppath_step\" and is the path
through/inside a single grid box.

The definition of a propgation path cannot be accommodated here.
For more information read the chapter on propagation paths in the
ARTS user guide.

This method should never be called directly. Use ``ppathCalc`` instead
if you want to extract propagation paths.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppath"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath_step_agenda",
         "ppath_inside_cloudbox_do",
         "atm_field",
         "f_grid",
         "surface_field",
         "cloudbox_on",
         "cloudbox_limits",
         "rte_pos",
         "rte_los",
         "ppath_lmax",
         "ppath_lraytrace"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
*/
/*
  md_data_raw.push_back(create_mdrecord(
      NAME("ppath_stepGeometric"),
      DESCRIPTION(R"--(Calculates a geometrical propagation path step.

This function determines a propagation path step by pure
geometrical calculations. That is, refraction is neglected. Path
points are always included for crossings with the grids, tangent
points and intersection points with the surface. The WSV *ppath_lmax*
gives the option to include additional points to ensure that the
distance along the path between the points does not exceed the
selected maximum length. No additional points are included if
*ppath_lmax* is set to <= 0.

For further information, type see the on-line information for
*ppath_step_agenda*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppath_step"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath_step",
         "atm_field",
         "surface_field",
         "ppath_lmax"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
*/

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("ppath_stepRefractionBasic"),
      DESCRIPTION(R"--(Calculates a propagation path step, considering refraction by a
basic approach.

Refraction is taken into account by probably the simplest approach
possible. The path is treated to consist of piece-wise geometric
steps. A geometric path step is calculated from each point by
using the local line-of-sight. Snell's law for spherical symmetry
is used for 1D to determine the zenith angle at the new point.
For 2D and 3D, the zenith angle is calculated using the average
gradient of the refractive index between the two points. For 3D,
the azimuth angle is treated in the same way as the zenith one.

The maximum length of each ray tracing step is given by the WSV
*ppath_lraytrace*. The length will never exceed the given maximum,
but it can be smaller. The ray tracing steps are only used to
determine the path. Points to describe the path are included as
for ``ppath_stepGeometric``, this including the functionality of
*ppath_lmax*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppath_step"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("refr_index_air_agenda",
         "ppath_step",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "atm_field",
         "surface_field",
         "f_grid",
         "ppath_lmax",
         "ppath_lraytrace"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("ppathAddGridCrossings"),
      DESCRIPTION(R"--(Adds grids crossings to an existing *ppath*.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppath"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath", "ppath_lstep", "surface_field"),
      GIN("z_grid",
          "lat_grid",
          "lon_grid"),
      GIN_TYPE("Vector", "Vector", "Vector"),
      GIN_DEFAULT("[]", "[]", "[]"),
      GIN_DESC("Grid/set of altitudes to include, if passed.",
               "Grid/set of latitudes to include, if passed.",
               "Grid/set of longitudes to include, if passed.")));

    md_data_raw.push_back(create_mdrecord(
      NAME("ppathCheckEndPoint"),
      DESCRIPTION(R"--(Checks that a propagation path ends as expected.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath"),
      GIN("background", "np",
          "altitude", "daltitude",
          "latitude", "dlatitude",
          "longitude", "dlongitude",
          "zenith_angle", "dzenith_angle",
          "azimuth_angle", "dazimuth_angle"),
      GIN_TYPE("String", "Index",
               "Numeric", "Numeric",
               "Numeric", "Numeric",
               "Numeric", "Numeric",
               "Numeric", "Numeric",
               "Numeric", "Numeric"),
      GIN_DEFAULT("Undefined", "-1","0","-1","0","-1","0","-1","0","-1","0","-1"),
      GIN_DESC("Expected radiative background. See above.",
               "Expected number of path points.",
               "Expected altitude.",
               "Allowed deviation for altitude.",
               "Expected latitude.",
               "Allowed deviation for latitude.",
               "Expected longitude.",
               "Allowed deviation for longitude.",
               "Expected zenith angle.",
               "Allowed deviation for zenith angle.",
               "Expected azimuth angle.",
               "Allowed deviation for azimuth angle.")));

    md_data_raw.push_back(create_mdrecord(
      NAME("ppathCheckInsideDomain"),
      DESCRIPTION(R"--(Checks that propagation path is fully inside specified domain.

An error is issued if a point in *ppath* is outside of ranges
[lat_min, lat_max] and [lon_min, lon_max], in latitude and
longitude, respectively. The path is checked starting at the end
furthest away from the sensor and an error is given as soon an
incorrect point is found. This is not guaranteed to be the point
most outside of the domain.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath"),
      GIN("lat_min", "lat_max", "lon_min", "lon_max"),
      GIN_TYPE("Numeric", "Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT("-90.0", "90.0", "-180.0", "360.0"),
      GIN_DESC("Lowest allowed latitude.",
               "Highest allowed latitude.",
               "Lowest allowed longitude.",
               "Highest allowed longitude.")));

    md_data_raw.push_back(create_mdrecord(
      NAME("ppathCheckInsideGrids"),
      DESCRIPTION(R"--(Checks that propagation path is fully inside grid ranges.

As *ppathCheckInsideDomain* but with the domain specified by a
combination of latitude and longitude grids.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath"),
      GIN("latitude_grid", "longitude_grid"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Latitude grid to not exceed.",
               "Longitude grid to not exceed.")));

    md_data_raw.push_back(create_mdrecord(
      NAME("ppathGeometric"),
      DESCRIPTION(
         R"(Geometric propagation path with fixed step length.

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
 )"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppath"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("rte_pos",
         "rte_los",
         "ppath_lstep",
         "ppath_ltotal",
         "surface_field",
         "surface_search_accuracy",
         "surface_search_safe",
         "atm_field"),
      GIN("include_specular_ppath"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("Flag to continue path after surface intersection.")));

    md_data_raw.push_back(create_mdrecord(
      NAME("ppathRefracted"),
      DESCRIPTION(R"--(Calculates propagation paths, including refraction by a basic approach.

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
removes the need for a \"safe\" search.

For more accurate calculations, but slower, consider the two GIN
parameters ``do_horizontal_gradients`` and ``do_twosided_perturb``

Default is to only determine the altitude gradients of the refractive
index, as this is in general the only relevant term. To also calculate
and consider the latitude and longitude gradients, set
``do_horizontal_gradients`` to true.

The gradients of the refractive index are obtained by perturbing the
position of concern with small positive values. With ``do_twosided_perturb``
set to true, there is also a perturbation in the negative direction.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppath"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("refr_index_air_ZZZ_agenda",
         "rte_pos",
         "rte_los",
         "ppath_lstep",
         "ppath_ltotal",
         "ppath_lraytrace",
         "surface_field",
         "surface_search_accuracy"),
      GIN("z_toa",
          "do_horizontal_gradients",
          "do_twosided_perturb",
          "include_specular_ppath"),
      GIN_TYPE("Numeric", "Index", "Index", "Index"),
      GIN_DEFAULT(NODEF, "0", "0", "0"),
      GIN_DESC("Top-of-the-atmosphere altitude.",
               "Consider horisontal gradients of refractive index.",
               "Perform double-sided perturbations when calculating "
               "refractive index gradients.",
               "See *ppathGeometric*.")));
               
  md_data_raw.push_back(create_mdrecord(
      NAME("ppvar_optical_depthFromPpvar_trans_cumulat"),
      DESCRIPTION(R"--(Sets *ppvar_optical_depth* according to provided transmittance data.

The values in ppvar_optical_depth are set to
-log( ppvar_trans_cumulat(joker,joker,0,0) ).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppvar_optical_depth"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppvar_trans_cumulat"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("ppvar_atmFromPath"),
      DESCRIPTION(R"--(Gets the atmospheric points along the path.
)--"),
      AUTHORS("Richard Larsson"), OUT("ppvar_atm"), GOUT(), GOUT_TYPE(),
      GOUT_DESC(), IN("ppath", "atm_field"), GIN(), GIN_TYPE(), GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("ppvar_fFromPath"),
      DESCRIPTION(R"--(Gets the frequency grid along the path.
)--"),
      AUTHORS("Richard Larsson"), OUT("ppvar_f"), GOUT(), GOUT_TYPE(),
      GOUT_DESC(), IN("f_grid", "ppath", "ppvar_atm", "rte_alonglos_v"), GIN(),
      GIN_TYPE(), GIN_DEFAULT(), GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("ppvar_propmatCalc"),
      DESCRIPTION(R"--(Gets the propagation matrix and NLTE source term along the path.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ppvar_propmat", "ppvar_nlte", "ppvar_dpropmat", "ppvar_dnlte"),
      GOUT(), GOUT_TYPE(), GOUT_DESC(),
      IN("propmat_clearsky_agenda", "jacobian_quantities", "ppvar_f", "ppath",
         "ppvar_atm", "jacobian_do"),
      GIN(), GIN_TYPE(), GIN_DEFAULT(), GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("ppvar_radCalc"),
      DESCRIPTION(R"--(Gets the radiation along the path.
)--"),
      AUTHORS("Richard Larsson"), OUT("ppvar_rad", "ppvar_drad"), GOUT(),
      GOUT_TYPE(), GOUT_DESC(),
      IN("background_rad", "ppvar_src", "ppvar_dsrc", "ppvar_tramat",
         "ppvar_cumtramat", "ppvar_dtramat", "ppvar_propmat", "ppvar_dpropmat",
         "ppvar_distance", "ppvar_ddistance", "rt_integration_option"),
      GIN(), GIN_TYPE(), GIN_DEFAULT(), GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("ppvar_radCalcEmission"),
      DESCRIPTION(R"--(Gets the radiation along the path by linear emission calculations.
)--"),
      AUTHORS("Richard Larsson"), OUT("ppvar_rad", "ppvar_drad"), GOUT(),
      GOUT_TYPE(), GOUT_DESC(),
      IN("background_rad", "ppvar_src", "ppvar_dsrc", "ppvar_tramat",
         "ppvar_cumtramat", "ppvar_dtramat"),
      GIN(), GIN_TYPE(), GIN_DEFAULT(), GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("ppvar_radCalcTransmission"),
      DESCRIPTION(R"--(Gets the radiation along the path by linear emission calculations.
)--"),
      AUTHORS("Richard Larsson"), OUT("ppvar_rad", "ppvar_drad"), GOUT(),
      GOUT_TYPE(), GOUT_DESC(),
      IN("ppvar_tramat",
         "ppvar_cumtramat", "ppvar_dtramat"),
      GIN(), GIN_TYPE(), GIN_DEFAULT(), GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("ppvar_srcFromPropmat"),
      DESCRIPTION(R"--(Gets the source term along the path.
)--"),
      AUTHORS("Richard Larsson"), OUT("ppvar_src", "ppvar_dsrc"), GOUT(),
      GOUT_TYPE(), GOUT_DESC(),
      IN("ppvar_propmat", "ppvar_nlte", "ppvar_dpropmat", "ppvar_dnlte",
         "ppvar_f", "ppvar_atm", "jacobian_quantities", "jacobian_do"),
      GIN(), GIN_TYPE(), GIN_DEFAULT(), GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("ppvar_tramatCalc"),
      DESCRIPTION(
          R"--(Gets the transmission matrix in layers along the path.

A layer is defined as made up by the average of 2 levels, thus the outer-most size
of the derivatives out of this function is 2.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ppvar_tramat", "ppvar_dtramat", "ppvar_distance", "ppvar_ddistance"),
      GOUT(), GOUT_TYPE(), GOUT_DESC(),
      IN("ppvar_propmat", "ppvar_dpropmat", "ppath", "ppvar_atm",
         "jacobian_quantities", "jacobian_do"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("background_radFromMatrix"),
      DESCRIPTION(
          R"--(Sets *background_rad* to matrix input
)--"),
      AUTHORS("Richard Larsson"),
      OUT("background_rad"),
      GOUT(), GOUT_TYPE(), GOUT_DESC(),
      IN(),
      GIN("iy_mat"),
      GIN_TYPE("Matrix"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Background radiation")));

  md_data_raw.push_back(create_mdrecord(
      NAME("background_transmittanceFromBack"),
      DESCRIPTION(
          R"--(Sets *background_transmittance* to back of *ppvar_cumtramat*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("background_transmittance"),
      GOUT(), GOUT_TYPE(), GOUT_DESC(),
      IN("ppvar_cumtramat"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("background_transmittanceFromFront"),
      DESCRIPTION(
          R"--(Sets *background_transmittance* to front of *ppvar_cumtramat*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("background_transmittance"),
      GOUT(), GOUT_TYPE(), GOUT_DESC(),
      IN("ppvar_cumtramat"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("ppvar_cumtramatForward"),
      DESCRIPTION(
          R"--(Sets *ppvar_cumtramat* by forward iteration of *ppvar_tramat*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ppvar_cumtramat"),
      GOUT(), GOUT_TYPE(), GOUT_DESC(),
      IN("ppvar_tramat"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("ppvar_cumtramatReverse"),
      DESCRIPTION(
          R"--(Sets *ppvar_cumtramat* by reverse iteration of *ppvar_tramat*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ppvar_cumtramat"),
      GOUT(), GOUT_TYPE(), GOUT_DESC(),
      IN("ppvar_tramat"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iyCopyPath"),
      DESCRIPTION(
          R"--(Copies the radiative transfer properties to their matpack equivalents.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("iy", "ppvar_iy", "ppvar_trans_cumulat", "ppvar_trans_partial",
          "diy_dpath"),
      GOUT(), GOUT_TYPE(), GOUT_DESC(),
      IN("ppvar_rad", "ppvar_drad", "ppvar_cumtramat", "ppvar_tramat",
          "jacobian_quantities", "jacobian_do"),
      GIN(), GIN_TYPE(), GIN_DEFAULT(), GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("diy_dxTransform"),
      DESCRIPTION(
          R"--(Transforms *diy_dpath* and adds it to *diy_dx*.
)--"),
      AUTHORS("Richard Larsson"), OUT("diy_dx", "diy_dpath"), GOUT(),
      GOUT_TYPE(), GOUT_DESC(),
      IN("diy_dx", "diy_dpath", "ppath", "ppvar_atm", "abs_species",
         "iy_transmittance", "water_p_eq_agenda", "jacobian_quantities",
         "jacobian_do", "iy_agenda_call1"),
      GIN(), GIN_TYPE(), GIN_DEFAULT(), GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iyUnitConversion"),
      DESCRIPTION(
          R"--(Perform unit conversions of *iy*, *diy_dx*, and *ppvar_iy*.
)--"),
      AUTHORS("Richard Larsson"), OUT("iy", "diy_dx", "ppvar_iy"), GOUT(), GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy", "diy_dx", "ppvar_iy", "f_grid", "ppath", "jacobian_quantities",
         "iy_unit", "jacobian_do", "iy_agenda_call1"),
      GIN(), GIN_TYPE(), GIN_DEFAULT(), GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("iy_auxFromVars"),
      DESCRIPTION(
          R"--(Set *iy_aux* from list of parameters.
)--"),
      AUTHORS("Richard Larsson"), OUT("iy_aux"), GOUT(), GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy_aux_vars", "background_transmittance", "ppath", "iy_agenda_call1"),
      GIN(), GIN_TYPE(), GIN_DEFAULT(), GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord(NAME("Print"),
               DESCRIPTION(R"--(Prints a variable on the screen.
)--"),
               AUTHORS("Oliver Lemke"),
               OUT(),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN("input", "level"),
               GIN_TYPE("Any", "Index"),
               GIN_DEFAULT(NODEF, "1"),
               GIN_DESC("Variable to be printed.", "Output level to use."),
               SETMETHOD(false),
               AGENDAMETHOD(false),
               USES_TEMPLATES(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("PrintPhysicalConstants"),
               DESCRIPTION(R"--(Prints (most) physical constants used in ARTS.
)--"),
               AUTHORS("Richard Larsson"),
               OUT(),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord(NAME("PrintWorkspace"),
               DESCRIPTION(R"--(Prints a list of the workspace variables.
)--"),
               AUTHORS("Oliver Lemke"),
               OUT(),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN("only_allocated", "level"),
               GIN_TYPE("Index", "Index"),
               GIN_DEFAULT("1", "1"),
               GIN_DESC("Flag for printing either all variables (0) or only "
                        "allocated ones (1).",
                        "Output level to use."),
               SETMETHOD(false),
               AGENDAMETHOD(false),
               USES_TEMPLATES(true),
               PASSWORKSPACE(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyAddCIA"),
      DESCRIPTION(R"--(Calculate absorption coefficients per tag group for HITRAN CIA continua.

This interpolates the cross sections from *abs_cia_data*.

The robust option is intended only for testing. Do not use for normal
runs, since subsequent functions will not be able to deal with NAN values.
)--"),
      AUTHORS("Stefan Buehler, Oliver Lemke"),
      OUT("propmat_clearsky", "dpropmat_clearsky_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
         "dpropmat_clearsky_dx",
         "abs_species",
         "select_abs_species",
         "jacobian_quantities",
         "f_grid",
         "atm_point",
         "abs_cia_data"),
      GIN("T_extrapolfac", "ignore_errors"),
      GIN_TYPE("Numeric", "Index"),
      GIN_DEFAULT("0.5", "0"),
      GIN_DESC(
          "Temperature extrapolation factor (relative to grid spacing).",
          "Set to 1 to suppress runtime errors (and return NAN values instead).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyAddFaraday"),
      DESCRIPTION(R"--(Calculates absorption matrix describing Faraday rotation.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("propmat_clearsky", "dpropmat_clearsky_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
         "dpropmat_clearsky_dx",
         
         "f_grid",
         "abs_species",
         "select_abs_species",
         "jacobian_quantities",
         "atm_point",
         "rtp_los"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyAddFromLookup"),
      DESCRIPTION(R"--(Extract gas absorption coefficients from lookup table.

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
)--"),
      AUTHORS("Stefan Buehler, Richard Larsson"),
      OUT("propmat_clearsky", "dpropmat_clearsky_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
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
         "select_abs_species"),
      GIN("extpolfac","no_negatives"),
      GIN_TYPE("Numeric","Index"),
      GIN_DEFAULT("0.5","1"),
      GIN_DESC("Extrapolation factor (for temperature and VMR grid edges).",
               "Boolean. If it is true negative values due to interpolation "
               "are set to zero.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyAddHitranLineMixingLines"),
      DESCRIPTION(R"--(Calculates gas absorption coefficients line-by-line for HITRAN line mixed data.

*Wigner6Init* or *Wigner3Init* must be called before this function.

Note that you need to have *propmat_clearskyAddLines* in addition to this method
to compensate the calculations for the pressure limit

Please ensure you cite the original authors when you use this function:
\tJ. Lamouroux, L. Realia, X. Thomas, et al., J.Q.S.R.T. 151 (2015), 88-96
)--"),
      AUTHORS("Richard Larsson"),
      OUT("propmat_clearsky"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
         "abs_hitran_relmat_data",
         "abs_lines_per_species",
         "isotopologue_ratios",
         "f_grid",
         "abs_species",
         "select_abs_species",
         "jacobian_quantities",
         "atm_point"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyAddScaledSpecies"),
      DESCRIPTION(
          R"--(Adds a scaled target species absorption to *propmat_clearsky* and *nlte_source*

This recomputes the entire propagation matrix.  There are more efficient ways
to do these calculations but this method exist because of the composability it
offers
)--"),
      AUTHORS("Richard Larsson"),
      OUT("propmat_clearsky", "nlte_source"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
         "nlte_source",
         "jacobian_quantities",
         "select_abs_species",
         "f_grid",
         "rtp_los",
         "atm_point",
         "propmat_clearsky_agenda"),
      GIN("target", "scale"),
      GIN_TYPE("ArrayOfSpeciesTag", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC(
          "Target species tags to rescale (must be in *abs_species*",
          R"(Rescaling factor (e.g., 0.1 adds 10% of the species to the absorption))")));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyAddXsecFit"),
      DESCRIPTION(R"--(Calculate absorption cross sections per tag group for HITRAN xsec species.

This broadens the cross section data from *xsec_fit_data* and
interpolates it onto the current f_grid.

Model data needs to be read in with *ReadXsecData* before calling
this method.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("propmat_clearsky", "dpropmat_clearsky_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
         "dpropmat_clearsky_dx",
         "abs_species",
         "select_abs_species",
         "jacobian_quantities",
         "f_grid",
         "atm_point",
         "xsec_fit_data"),
      GIN("force_p",
          "force_t"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT("-1", "-1"),
      GIN_DESC("Positive value forces constant pressure [Pa].",
               "Positive value forces constant temperature [K].")));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyAddLines"),
      DESCRIPTION(
        R"--(Computes the line-by-line unpolarized absorption and adds
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
)--"
      ),
      AUTHORS("Richard Larsson"),
      OUT("propmat_clearsky", "nlte_source", "dpropmat_clearsky_dx", "dnlte_source_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
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
         "lbl_checked"),
      GIN("lines_sparse_df", "lines_sparse_lim", "lines_speedup_option", "no_negatives"),
      GIN_TYPE("Numeric", "Numeric", "String", "Index"),
      GIN_DEFAULT("0", "0", "None", "1"),
      GIN_DESC(
        "The grid sparse separation",
        "The dense-to-sparse limit",
        "Speedup logic",
        "Boolean.  If it is true, line mixed bands each allocate their own compute data to ensure that they cannot produce negative absorption"
      )));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyAddOnTheFlyLineMixing"),
      DESCRIPTION(R"--(Compute the line mixing of matching lines and add it to the propagation matrix

Each band's Population Type is checked and the calculations are only performed
for bands with matching population types (and a good pressure limits)

Presently only supports one method: ByMakarovFullRelmat, based on Makarov et al 2020

*Wigner6Init* or *Wigner3Init* must be called before this function.

Note that you need to have *propmat_clearskyAddLines* addition to this method
to compensate the calculations for the pressure limit
)--"),
      AUTHORS("Richard Larsson"),
      OUT("propmat_clearsky", "dpropmat_clearsky_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
         "dpropmat_clearsky_dx",
         "abs_lines_per_species",
         "ecs_data",
         "isotopologue_ratios",
         "f_grid",
         "abs_species",
         "select_abs_species",
         "jacobian_quantities",
         "atm_point",
         "lbl_checked"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyAddOnTheFlyLineMixingWithZeeman"),
      DESCRIPTION(R"--(Compute the line mixing of matching lines and add it to the propagation matrix
Also computes Zeeman effect for all the lines in the band

Each band's Population Type is checked and the calculations are only performed
for bands with matching population types (and a good pressure limits)

Presently only supports one method: ByMakarovFullRelmat, based on Makarov et al 2020

*Wigner6Init* or *Wigner3Init* must be called before this function.

Note that you need to have *propmat_clearskyAddLines* in addition to this method
to compensate the calculations for the pressure limit
)--"),
      AUTHORS("Richard Larsson"),
      OUT("propmat_clearsky", "dpropmat_clearsky_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
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
         "lbl_checked"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyAddParticles"),
      DESCRIPTION(R"--(Calculates absorption coefficients of particles to be used in
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
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("propmat_clearsky", "dpropmat_clearsky_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
         "dpropmat_clearsky_dx",
         
         "f_grid",
         "abs_species",
         "select_abs_species",
         "jacobian_quantities",
         "rtp_los",
         "atm_point",
         "scat_data",
         "scat_data_checked"),
      GIN("use_abs_as_ext"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("A flag with value 1 or 0. If set to one, particle absorption "
               "is used in extinction and emission parts of the RT equation, "
               "and scattering out of LOS as well as into LOS is neglected. "
               "Otherwise, particle extinction (absorption+scattering) is "
               "applied in both the extinction as well as the emission part "
               "of the RT equation. That is, true extinction is applied, but "
               "emission also includes a pseudo-emission contribution from "
               "the scattering coefficient. ")));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyAddZeeman"),
      DESCRIPTION(R"--(Calculates Zeeman-affected polarized propagation matrix and its
derivatives.

Otherwise as *propmat_clearskyAddFromLookup*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("propmat_clearsky",
          "nlte_source",
          "dpropmat_clearsky_dx",
          "dnlte_source_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
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
         "lbl_checked"),
      GIN("manual_mag_field", "H", "theta", "eta"),
      GIN_TYPE("Index", "Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT("0", "1.0", "0.0", "0.0"),
      GIN_DESC("Manual angles tag",
               "Manual Magnetic Field Strength",
               "Manual theta given positive tag",
               "Manual eta given positive tag")));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyInit"),
      DESCRIPTION(R"--(Initialize *propmat_clearsky*, *nlte_source*, and their derivatives to zeroes.

This method must be used inside *propmat_clearsky_agenda* and then be called first.
)--"),
      AUTHORS("Oliver Lemke", "Richard Larsson"),
      OUT("propmat_clearsky",
          "nlte_source",
          "dpropmat_clearsky_dx",
          "dnlte_source_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "f_grid",
         
         "propmat_clearsky_agenda_checked"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyZero"),
      DESCRIPTION(R"--(Sets *propmat_clearsky* to match zero attenuation.

Use this method just if you know what you are doing!

If you want to make a calculation with no clear-sky attenuation at
all, fill *propmat_clearsky_agenda* with this method and required
Ignore statements (don't include *propmat_clearskyInit*).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("propmat_clearsky"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearskyForceNegativeToZero"),
      DESCRIPTION(R"--(Sets *propmat_clearsky* to match zero attenuation
if negative value.  Useful for line mixing in some cases.

Use this method just if you know what you are doing!
)--"),
      AUTHORS("Richard Larsson"),
      OUT("propmat_clearsky"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearsky_agenda_checkedCalc"),
      DESCRIPTION(R"--(Checks if the *propmat_clearsky_agenda* contains all necessary
methods to calculate all the species in *abs_species*.

This method should be called just before the *propmat_clearsky_agenda*
is used, e.g. *DoitGetIncoming*, *ybatchCalc*, *yCalc*
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("propmat_clearsky_agenda_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species", "propmat_clearsky_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("propmat_clearsky_fieldCalc"),
      DESCRIPTION(R"--(Calculate (vector) gas absorption coefficients for all points in the
atmosphere.

This is useful in two different contexts:

1. For testing and plotting gas absorption. (For RT calculations, gas
   absorption is calculated or extracted locally, therefore there is no
   need to calculate a global field. But this method is handy for easy
   plotting of absorption vs. pressure, for example.)

2. Inside the scattering region, monochromatic absorption is
   pre-calculated for the entire atmospheric field.

The calculation itself is performed by the
*propmat_clearsky_agenda*.
)--"),
      AUTHORS("Stefan Buehler, Richard Larsson"),
      OUT("propmat_clearsky_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmfields_checked",
         "f_grid",
         "atm_field",
         "abs_species",
         "propmat_clearsky_agenda"),
      GIN("doppler", "los"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT("[]", "[]"),
      GIN_DESC("A vector of doppler shift values in Hz. Must either be "
               "empty or have same dimension as p_grid.",
               "Line of sight")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdAbelBoutle12"),
      DESCRIPTION(R"--(Abel and Boutle [2012] particle size distribution for rain.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names",
         "scat_species_a",
         "scat_species_b"),
      GIN("t_min", "t_max", "picky"),
      GIN_TYPE("Numeric", "Numeric", "Index"),
      GIN_DEFAULT("273", "373", "0"),
      GIN_DESC(
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdDelanoeEtAl14"),
      DESCRIPTION(R"--(Normalized PSD as proposed in Delanoë et al. ((2014)),

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
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names"),
      GIN("iwc",
          "n0Star",
          "Dm",
          "rho",
          "alpha",
          "beta",
          "t_min",
          "t_max",
          "dm_min",
          "picky"),
      GIN_TYPE("Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Index"),
      GIN_DEFAULT("NaN",
                  "NaN",
                  "NaN",
                  "916.7",
                  "-0.237",
                  "1.839",
                  NODEF,
                  NODEF,
                  "-1.0",
                  "0"),
      GIN_DESC(
          "Ice water content",
          "Intercept parameter",
          "Volume weighted diameter",
          "Density of ice",
          "``alpha`` parameter of the shape function",
          "``beta`` paramter of the shape function",
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Lower threshold for ``Dm`` below which an error is thrown.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdFieldEtAl07"),
      DESCRIPTION(R"--(The Field et al. [2007] particle size distribution for snow and
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
governed by setting of ``regime``, where \"TR\" selectes the tropical
case, and \"ML\" the midlatitude one.

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
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names",
         "scat_species_a",
         "scat_species_b"),
      GIN("regime",
          "t_min",
          "t_max",
          "t_min_psd",
          "t_max_psd",
          "beta_min",
          "beta_max",
          "picky"),
      GIN_TYPE("String",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Index"),
      GIN_DEFAULT(NODEF, "0", "290.", "200.", "273.15", "1.01", "4", "0"),
      GIN_DESC(
          "Parametrization regime (\"TR\"=tropical or \"ML\"=midlatitude).",
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Low temperature limit to use as paramtrization temperature.",
          "High temperature limit to use as paramtrization temperature.",
          "Low ``b`` limit (only if picky).",
          "High ``b`` limit (only if picky).",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdFieldEtAl19"),
      DESCRIPTION(R"--(The Field [2019] particle size distribution for hail.

Reference: Field, Normalized hail particle size distributions from the T-28
storm-penetrating aircraft, JAMC, 2019

This is a 1-parmater PSD i.e. *pnd_agenda_input* shall have one column and
*pnd_agenda_input_names* shall contain a single string.
The input data in *pnd_agenda_input* shall be hail mass content in
unit of [kg/m3]. The naming used is *pnd_agenda_input_names* is free
but the same name must be used in *particle_bulkprop_names* and
*dpnd_data_dx_names*.
The parameters assume a constant effective density, i.e. scat_species_b \approx 3

Derivatives are obtained analytically.

The validity range of mass content is not limited. Negative mass
contents will produce negative psd values following a distribution
given by abs(HWC), ie. abs(psd)=f(abs(HWC)).

If temperature is outside [ ``t_min`` , ``t_max`` ] psd=0 and dpsd=0 if
picky=0, or an error is thrown if picky=1.
)--"),
      AUTHORS("Stuart Fox"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names",
         "scat_species_a",
         "scat_species_b"),
      GIN("t_min", "t_max", "picky"),
      GIN_TYPE("Numeric", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, "0"),
      GIN_DESC(
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdMcFarquaharHeymsfield97"),
      DESCRIPTION(R"--(McFarquahar and Heymsfield [1997] particle size distribution
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
)--"),
      AUTHORS("Patrick Eriksson, Jana Mendrok"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names",
         "scat_species_a",
         "scat_species_b"),
      GIN("t_min", "t_max", "t_min_psd", "t_max_psd", "picky", "noisy"),
      GIN_TYPE("Numeric", "Numeric", "Numeric", "Numeric", "Index", "Index"),
      GIN_DEFAULT("0", "280.", "180", "273.15", "0", "0"),
      GIN_DESC("Low temperature limit to calculate a psd.",
               "High temperature limit to calculate a psd.",
               "Low temperature limit to use as paramtrization temperature.",
               "High temperature limit to use as paramtrization temperature.",
               "Flag whether to be strict with parametrization value checks.",
               "Distribution parameter perturbance flag")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdMilbrandtYau05"),
      DESCRIPTION(R"--(Calculates *psd_data* and  *dpsd_data_dx* following Milbrandt and Yau (2005)
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

(1) \"cloud_water\" selects cloud liquid water , 
(2) \"cloud_ice\" selects cloud ice, 
(3) \"snow\" selects snow, 
(4) \"rain\" selects rain, 
(5) \"graupel\" selects graupel, and 
(6) \"hail\" selects hail, 

Requirements:

*pnd_agenda_input_names* must include::

    [\"X-mass_density\", \"X-number_density\" ]. \"X\" is an arbitrary name

The entries in  *dpnd_data_dx_names* (ie. the allowed
independent variablea ) can be \"X-mass_density\" and\\or 
\"X-number_density\".

The validity range of WC is not limited. Negative WC will produce
negative psd values following a distribution given by abs(WC), ie.
abs(psd)=f(abs(WC)).

If temperature is outside [``t_min``,``t_max``] psd=0 and dpsd=0 if
picky=0, or an error is thrown if picky=1.
)--"),
      AUTHORS("Manfred Brath"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names"),
      GIN("hydrometeor_type", "t_min", "t_max", "picky"),
      GIN_TYPE("String", "Numeric", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, "0", "999", "0"),
      GIN_DESC(
          "Hydrometeor type (see above description).",
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdModifiedGamma"),
      DESCRIPTION(R"--(Modified gamma distribution PSD using n0, mu, la and ga as parameters.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names"),
      GIN("n0", "mu", "la", "ga", "t_min", "t_max", "picky"),
      GIN_TYPE("Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Index"),
      GIN_DEFAULT("NaN", "NaN", "NaN", "NaN", NODEF, NODEF, "0"),
      GIN_DESC(
          "n0",
          "mu",
          "la",
          "ga",
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdModifiedGammaMass"),
      DESCRIPTION(R"--(Modified gamma distribution (MGD) PSD, with mass content as input.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names",
         "scat_species_a",
         "scat_species_b"),
      GIN("n0", "mu", "la", "ga", "t_min", "t_max", "picky"),
      GIN_TYPE("Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Index"),
      GIN_DEFAULT("NaN", "NaN", "NaN", "NaN", NODEF, NODEF, "0"),
      GIN_DESC(
          "n0",
          "mu",
          "la",
          "ga",
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdModifiedGammaMassNtot"),
      DESCRIPTION(R"--(Modified gamma distribution PSD, with mass content and total number
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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names",
         "scat_species_a",
         "scat_species_b"),
      GIN("n0", "mu", "la", "ga", "t_min", "t_max", "picky"),
      GIN_TYPE("Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Index"),
      GIN_DEFAULT("NaN", "NaN", "NaN", "NaN", NODEF, NODEF, "0"),
      GIN_DESC(
          "n0",
          "mu",
          "la",
          "ga",
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdModifiedGammaMassMeanParticleMass"),
      DESCRIPTION(R"--(Modified gamma distribution PSD, with mass content and mean particle
mass (Mmean) as inputs.

\"Mean particle mass\" is here defined as the mass content divided with
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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names",
         "scat_species_a",
         "scat_species_b"),
      GIN("n0", "mu", "la", "ga", "t_min", "t_max", "picky"),
      GIN_TYPE("Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Index"),
      GIN_DEFAULT("NaN", "NaN", "NaN", "NaN", NODEF, NODEF, "0"),
      GIN_DESC(
          "n0",
          "mu",
          "la",
          "ga",
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdModifiedGammaMassSingleMoment"),
      DESCRIPTION(R"--(Modified gamma distribution PSD, with mass content as input.

The intercept parameter N0 is assumed dependent on the slope parameter
lambda, such that N0=N_alpha*lambda^n_b with fixed N_alpha and n_b.
This is a common form for many PSD parametrizations for use with
single-moment mass-based schemes.

This version of MGD PSD takes mass content as first input argument.
This means that the first column of *pnd_agenda_input* shall hold
mass content data. The dependent parameter is assumed to be lambda.
)--"),
      AUTHORS("Stuart Fox"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names",
         "scat_species_a",
         "scat_species_b"),
      GIN("n_alpha", "n_b", "mu", "gamma", "t_min", "t_max", "picky"),
      GIN_TYPE("Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, NODEF, NODEF, "0"),
      GIN_DESC(
          "n_alpha",
          "n_b",
          "mu",
          "gamma",
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdModifiedGammaMassXmean"),
      DESCRIPTION(R"--(Modified gamma distribution PSD, with mass content and mean size
(Xmean) as inputs.

\"Mean size\" is here defined as mass weighted size. Remembering that
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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names",
         "scat_species_a",
         "scat_species_b"),
      GIN("n0", "mu", "la", "ga", "t_min", "t_max", "picky"),
      GIN_TYPE("Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Index"),
      GIN_DEFAULT("NaN", "NaN", "NaN", "NaN", NODEF, NODEF, "0"),
      GIN_DESC(
          "n0",
          "mu",
          "la",
          "ga",
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdModifiedGammaMassXmedian"),
      DESCRIPTION(R"--(Modified gamma distribution PSD, with mass content and median size
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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names",
         "scat_species_a",
         "scat_species_b"),
      GIN("n0", "mu", "la", "ga", "t_min", "t_max", "picky"),
      GIN_TYPE("Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Index"),
      GIN_DEFAULT("NaN", "NaN", "NaN", "NaN", NODEF, NODEF, "0"),
      GIN_DESC(
          "n0",
          "mu",
          "la",
          "ga",
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdMonoDispersive"),
      DESCRIPTION(R"--(Mono-dispersive PSD, with number density given.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names",
         "scat_meta"),
      GIN("species_index", "t_min", "t_max", "picky"),
      GIN_TYPE("Index", "Numeric", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "0"),
      GIN_DESC(
          "The index of the scattering species of concern (0-based).",
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdMonoMass"),
      DESCRIPTION(R"--(Mono-dispersive PSD, with mass content given.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names",
         "scat_meta"),
      GIN("species_index", "t_min", "t_max", "picky"),
      GIN_TYPE("Index", "Numeric", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "0"),
      GIN_DESC(
          "The index of the scattering species of concern (0-based).",
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdSeifertBeheng06"),
      DESCRIPTION(R"--(Calculates *psd_data* and *dpsd_data_dx* following Seifert and Beheng (2006)
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

(1) \"cloud_water\" selects cloud liquid water , 
(2) \"cloud_ice\" selects cloud ice, 
(3) \"snow\" selects snow, 
(4) \"rain\" selects rain, 
(5) \"graupel\" selects graupel, and 
(6) \"hail\" selects hail, 

Requirements:

*pnd_agenda_input_names* must include::

  [\"X-mass_density\", \"X-number_density\" ]. \"X\" is an arbitrary name

The entries in  *dpnd_data_dx_names* (ie. the allowed
independent variablea ) can be \"X-mass_density\" and\\or 
\"X-number_density\".

The validity range of WC is not limited. Negative WC will produce
negative psd values following a distribution given by abs(WC), ie.
abs(psd)=f(abs(WC)).

If temperature is outside [ ``t_min`` , ``t_max`` ] psd=0 and dpsd=0 if
picky=0, or an error is thrown if picky=1.
)--"),
      AUTHORS("Manfred Brath"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names"),
      GIN("hydrometeor_type", "t_min", "t_max", "picky"),
      GIN_TYPE("String", "Numeric", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, "0", "999", "0"),
      GIN_DESC(
          "Hydrometeor type (see above description).",
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("psdWangEtAl16"),
      DESCRIPTION(R"--(Wang et al. [2016] particle size distribution for rain.

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
)--"),
      AUTHORS("Jana Mendrok, Patrick Eriksson"),
      OUT("psd_data", "dpsd_data_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("psd_size_grid",
         "pnd_agenda_input_t",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "dpnd_data_dx_names",
         "scat_species_a",
         "scat_species_b"),
      GIN("t_min", "t_max", "picky"),
      GIN_TYPE("Numeric", "Numeric", "Index"),
      GIN_DEFAULT("273", "373", "0"),
      GIN_DESC(
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("RadarOnionPeelingTableCalc"),
      DESCRIPTION(R"--(Creates a radar inversion table.

This method is tailored to make inversion tables that fit
*particle_bulkpropRadarOnionPeeling*. See that method for
format of the table.

The method needs to be called twice to form a complete table,
once for liquid and ice hydrometeors. The table can be empty at
the first call.

The input data (*scat_data* etc.) must match two scattering
species and a single frequency (the one of the radar).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("invtable"),
      GOUT_TYPE("ArrayOfGriddedField3"),
      GOUT_DESC(""),
      IN("f_grid",
         "scat_species",
         "scat_data",
         "scat_meta",
         "pnd_agenda_array",
         "pnd_agenda_array_input_names"),
      GIN("i_species",
          "dbze_grid",
          "t_grid",
          "wc_min",
          "wc_max",
          "ze_tref",
          "k2"),
      GIN_TYPE("Index", "Vector", "Vector", "Numeric", "Numeric",
               "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "1e-8", "2e-2", "273.15", "-1" ),
      GIN_DESC("Index of *scat_species* to do. Can be 0 or 1.",
               "Grid of dBZe values to use for the table.",
               "Temperature grid to use for the table.",
               "A water content value that gives a dBZe smaller than first "
               "value of ``dbze_grid``.",
               "A water content value that gives a dBZe larger than last "
               "value of ``dbze_grid``.",
               "Reference temperature for conversion to Ze. See further *yRadar*.",
               "Reference dielectric factor. See further *yRadar*.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("RadiationFieldSpectralIntegrate"),
      DESCRIPTION(R"--(Integrates fields like *spectral_irradiance_field* or *spectral_radiance_field* over frequency.

Important, the first dimension must be the frequency dimension!
If a field  like *spectral_radiance_field* is input, the stokes dimension
is also removed.
)--"),
      AUTHORS("Manfred Brath"),
      OUT(),
      GOUT("radiation_field"),
      GOUT_TYPE("Tensor4, Tensor5"),
      GOUT_DESC("Field similar to irradiance field or spectral irradiance field"),
      IN("f_grid"),
      GIN("spectral_radiation_field"),
      GIN_TYPE("Tensor5, Tensor7"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Field similar to spectral irradiance field, spectral radiance field")));

  md_data_raw.push_back(create_mdrecord(
      NAME("line_irradianceCalcForSingleSpeciesNonOverlappingLinesPseudo2D"),
      DESCRIPTION(R"--(Computes the line irradiance and line transmission

Presently only works for 1D atmospheres
)--"),
      AUTHORS("Richard Larsson"),
      OUT("line_irradiance", "line_transmission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species",
         "abs_lines_per_species",
         "atm_field",
         "surface_field",
         "ppath_agenda",
         "iy_main_agenda",
         "iy_space_agenda",
         "iy_surface_agenda",
         "iy_cloudbox_agenda",
         "propmat_clearsky_agenda"),
      GIN("df", "nz", "nf", "r"),
      GIN_TYPE("Numeric", "Index", "Index", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "1.0"),
      GIN_DESC("relative frequency to line center",
               "number of zeniths",
               "number of frequencies per line",
               "Distance assumed when computing local (1-T)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("RadiativePropertiesCalc"),
      DESCRIPTION(R"(Wraps executing *ppvar_rtprop_agenda*.
)"),
      AUTHORS("Richard Larsson"),
      OUT("ppvar_propmat", "ppvar_dpropmat", "ppvar_src", "ppvar_dsrc",
          "ppvar_tramat", "ppvar_dtramat", "ppvar_distance", "ppvar_ddistance",
          "ppvar_cumtramat"),
      GOUT(), GOUT_TYPE(), GOUT_DESC(),
      IN("ppath", "ppvar_atm", "ppvar_f", "jacobian_do", "ppvar_rtprop_agenda"),
      GIN(), GIN_TYPE(), GIN_DEFAULT(), GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("RadiationBackgroundCalc"),
      DESCRIPTION(R"(Wraps executing *rte_background_agenda*.
)"),
      AUTHORS("Richard Larsson"), OUT("background_rad", "diy_dx"), GOUT(),
      GOUT_TYPE(), GOUT_DESC(),
      IN("ppath", "atm_field", "f_grid", "iy_transmittance",
         "background_transmittance", "jacobian_do", "rte_background_agenda"),
      GIN(), GIN_TYPE(), GIN_DEFAULT(), GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("RationalAdd"),
      DESCRIPTION(R"--(Adds a Rational and a value (output = input + value).

The result can either be stored in the same or another Rational.
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Rational"),
      GOUT_DESC("Output Rational."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Rational", "Rational"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Rational.", "Value to add.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("RationalDivide"),
      DESCRIPTION(R"--(Divides a Rational with a value (output = input / value).

The result can either be stored in the same or another Rational.
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Rational"),
      GOUT_DESC("Output Rational."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Rational", "Rational"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Rational.", "Denominator.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("RationalMultiply"),
      DESCRIPTION(R"--(Multiplies a Rational with a value (output = input * value).

The result can either be stored in the same or another Rational.
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Rational"),
      GOUT_DESC("Output Rational."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Rational", "Rational"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Rational.", "Multiplier.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("RationalSubtract"),
      DESCRIPTION(R"--(Subtracts a Rational value (output = input - value).

The result can either be stored in the same or another Rational.
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Rational"),
      GOUT_DESC("Output Rational."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Rational", "Rational"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Rational.", "Subtrahend.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ReadArrayOfARTSCAT"),
      DESCRIPTION(R"--(Reads an old Array<ArrayOfLineRecord> ARTSCAT file.

Note that the ARTSCAT-5 had quantum numbers and options
stored inside it but that the options will overwrite that
information.  Be careful setting the options!
)--"),
      AUTHORS("Stefan Buehler", "Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("filename", "fmin", "fmax", "globalquantumnumbers",
          "localquantumnumbers", "normalization_option", "mirroring_option",
          "population_option", "lineshapetype_option", "cutoff_option",
          "cutoff_value", "linemixinglimit_value"),
      GIN_TYPE("String", "Numeric", "Numeric", "String", "String", "String",
               "String", "String", "String", "String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "-1e99", "1e99", "", "", "None", "None", "LTE", "VP",
                 "None", "750e9", "-1"),
      GIN_DESC("Name of the ARTSCAT file",
               "Minimum frequency of read lines",
               "Maximum frequency of read lines",
               "Global quantum number list (space-separated)",
               "Local quantum number list (space-separated)",
               "Normalization option, see *abs_linesNormalization*",
               "Mirroring option, see *abs_linesMirroring*",
               "Population option, see *abs_linesPopulation*",
               "Lineshape option, see *abs_linesLineShapeType*",
               "Cutoff option, see *abs_linesCutoff*",
               "Cutoff value, see *abs_linesCutoff*",
               "Line mixing limit, see *abs_linesLinemixingLimit*")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ReadSplitARTSCAT"),
      DESCRIPTION(R"--(Reads several old ArrayOfLineRecord ARTSCAT file

Note that the ARTSCAT-5 had quantum numbers and options
stored inside it but that the options will overwrite that
information.  Be careful setting the options!
)--"),
      AUTHORS("Oliver Lemke", "Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("basename", "fmin", "fmax", "globalquantumnumbers", "localquantumnumbers",
          "ignore_missing", "normalization_option", "mirroring_option",
          "population_option", "lineshapetype_option", "cutoff_option",
          "cutoff_value", "linemixinglimit_value"),
      GIN_TYPE("String", "Numeric", "Numeric", "String", "String", "Index", "String",
               "String", "String", "String", "String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "-1e99", "1e99", "", "", "0", "None", "None", "LTE", "VP",
                 "None", "750e9", "-1"),
      GIN_DESC("Path to the files",
               "Minimum frequency of read lines",
               "Maximum frequency of read lines",
               "Global quantum number list (space-separated)",
               "Local quantum number list (space-separated)",
               "Ignores instead of throws if an *abs_species* is missing",
               "Normalization option, see *abs_linesNormalization*",
               "Mirroring option, see *abs_linesMirroring*",
               "Population option, see *abs_linesPopulation*",
               "Lineshape option, see *abs_linesLineShapeType*",
               "Cutoff option, see *abs_linesCutoff*",
               "Cutoff value, see *abs_linesCutoff*",
               "Line mixing limit, see *abs_linesLinemixingLimit*")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ReadARTSCAT"),
      DESCRIPTION(R"--(Reads an old ArrayOfLineRecord ARTSCAT file

Note that the ARTSCAT-5 had quantum numbers and options
stored inside it but that the options will overwrite that
information.  Be careful setting the options!
)--"),
      AUTHORS("Stefan Buehler", "Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("filename", "fmin", "fmax", "globalquantumnumbers",
          "localquantumnumbers", "normalization_option", "mirroring_option",
          "population_option", "lineshapetype_option", "cutoff_option",
          "cutoff_value", "linemixinglimit_value"),
      GIN_TYPE("String", "Numeric", "Numeric", "String", "String", "String",
               "String", "String", "String", "String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "-1e99", "1e99", "", "", "None", "None", "LTE", "VP",
                 "None", "750e9", "-1"),
      GIN_DESC("Name of the ARTSCAT file",
               "Minimum frequency of read lines",
               "Maximum frequency of read lines",
               "Global quantum number list (space-separated)",
               "Local quantum number list (space-separated)",
               "Normalization option, see *abs_linesNormalization*",
               "Mirroring option, see *abs_linesMirroring*",
               "Population option, see *abs_linesPopulation*",
               "Lineshape option, see *abs_linesLineShapeType*",
               "Cutoff option, see *abs_linesCutoff*",
               "Cutoff value, see *abs_linesCutoff*",
               "Line mixing limit, see *abs_linesLinemixingLimit*")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ReadHITRAN"),
      DESCRIPTION(R"--(Reads a HITRAN .par file.

The HITRAN type switch can be:

- ``\"Pre2004\"``: for old format
- ``\"Post2004\"``: for new format
- ``\"Online\"``: for the online format with quantum numbers (recommended)

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
)--"),
      AUTHORS("Hermann Berg", "Thomas Kuhn", "Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("filename", "fmin", "fmax", "globalquantumnumbers", "localquantumnumbers",
          "hitran_type", "normalization_option", "mirroring_option",
          "population_option", "lineshapetype_option", "cutoff_option",
          "cutoff_value", "linemixinglimit_value"),
      GIN_TYPE("String", "Numeric", "Numeric", "String", "String", "String", "String",
               "String", "String", "String", "String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "-1e99", "1e99", "DEFAULT_GLOBAL", "DEFAULT_LOCAL", "Online", "None", "None", "LTE", "VP",
                 "None", "750e9", "-1"),
      GIN_DESC("Name of the HITRAN file",
               "Minimum frequency of read lines",
               "Maximum frequency of read lines",
               "Global quantum number list (space-separated, default gives all)",
               "Local quantum number list (space-separated, default gives all)",
               "Method to use to read the line data",
               "Normalization option, see *abs_linesNormalization*",
               "Mirroring option, see *abs_linesMirroring*",
               "Population option, see *abs_linesPopulation*",
               "Lineshape option, see *abs_linesLineShapeType*",
               "Cutoff option, see *abs_linesCutoff*",
               "Cutoff value, see *abs_linesCutoff*",
               "Line mixing limit, see *abs_linesLinemixingLimit*")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ReadLBLRTM"),
      DESCRIPTION(R"--(Reads a LBLRTM file.

Be careful setting the options!
)--"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("filename", "fmin", "fmax", "globalquantumnumbers",
          "localquantumnumbers", "normalization_option", "mirroring_option",
          "population_option", "lineshapetype_option", "cutoff_option",
          "cutoff_value", "linemixinglimit_value"),
      GIN_TYPE("String", "Numeric", "Numeric", "String", "String", "String",
               "String", "String", "String", "String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "-1e99", "1e99", "", "", "None", "None", "LTE", "VP",
                 "None", "750e9", "-1"),
      GIN_DESC("Name of the LBLRTM file",
               "Minimum frequency of read lines",
               "Maximum frequency of read lines",
               "Global quantum number list (space-separated)",
               "Local quantum number list (space-separated)",
               "Normalization option, see *abs_linesNormalization*",
               "Mirroring option, see *abs_linesMirroring*",
               "Population option, see *abs_linesPopulation*",
               "Lineshape option, see *abs_linesLineShapeType*",
               "Cutoff option, see *abs_linesCutoff*",
               "Cutoff value, see *abs_linesCutoff*",
               "Line mixing limit, see *abs_linesLinemixingLimit*")));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("ReadJPL"),
      DESCRIPTION(R"--(Reads a JPL file.

Be careful setting the options!
)--"),
      AUTHORS("Thomas Kuhn", "Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("filename", "fmin", "fmax", "globalquantumnumbers",
          "localquantumnumbers", "normalization_option", "mirroring_option",
          "population_option", "lineshapetype_option", "cutoff_option",
          "cutoff_value", "linemixinglimit_value"),
      GIN_TYPE("String", "Numeric", "Numeric", "String", "String", "String",
               "String", "String", "String", "String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "-1e99", "1e99", "", "", "None", "None", "LTE", "VP",
                 "None", "750e9", "-1"),
      GIN_DESC("Name of the JPL file",
               "Minimum frequency of read lines",
               "Maximum frequency of read lines",
               "Global quantum number list (space-separated)",
               "Local quantum number list (space-separated)",
               "Normalization option, see *abs_linesNormalization*",
               "Mirroring option, see *abs_linesMirroring*",
               "Population option, see *abs_linesPopulation*",
               "Lineshape option, see *abs_linesLineShapeType*",
               "Cutoff option, see *abs_linesCutoff*",
               "Cutoff value, see *abs_linesCutoff*",
               "Line mixing limit, see *abs_linesLinemixingLimit*")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_linesWriteSpeciesSplitCatalog"),
      DESCRIPTION(R"--(Writes a split catalog, AbsorptionLines by AbsorptionLines.

There will be one unique file generated per AbsorptionLines in abs_lines.

The names of these files will be::

	basename + "." + AbsorptionLines.SpeciesName() + "." + to_string(N) + ".xml"

where N>=0 and the species name is something line "H2O".
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("output_file_format", "abs_lines"),
      GIN("basename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Path to store the files at")));

  md_data_raw.push_back(create_mdrecord(
      NAME("abs_lines_per_speciesWriteSpeciesSplitCatalog"),
      DESCRIPTION(R"--(See *abs_linesWriteSpeciesSplitCatalog*

In addition, the structure of the files generated will not care about
generating identifiers for the order in *abs_species*
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("output_file_format", "abs_lines_per_species"),
      GIN("basename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Path to store the files at")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ReadXsecData"),
      DESCRIPTION(R"--(Reads HITRAN Crosssection coefficients

Reads coefficient files for HITRAN Xsec species defined
in *abs_species*.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("xsec_fit_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("basename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Basepath to the files")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ReadNetCDF"),
      DESCRIPTION(R"--(Reads a workspace variable from a NetCDF file.

This method can read variables of any group.

If the filename is omitted, the variable is read
from <basename>.<variable_name>.nc.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector, Matrix, Tensor3, Tensor4, Tensor5, ArrayOfVector,"
                "ArrayOfIndex, ArrayOfMatrix, GasAbsLookup"),
      GOUT_DESC("Variable to be read."),
      IN(),
      GIN("filename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Name of the NetCDF file."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("ReadXML"),
      DESCRIPTION(R"--(Reads a workspace variable from an XML file.

This method can read variables of any group.

If the filename is omitted, the variable is read
from <basename>.<variable_name>.xml.
If the given filename does not exist, this method will
also look for files with an added .xml, .xml.gz and .gz extension
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Any"),
      GOUT_DESC("Variable to be read."),
      IN(),
      GIN("filename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(""),
      GIN_DESC("Name of the XML file."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("ReadXMLIndexed"),
      DESCRIPTION(R"--(As *ReadXML*, but reads indexed file names.

The variable is read from a file with name::

   <filename>.<file_index>.xml.

where <file_index> is the value of ``file_index``.

This means that ``filename`` shall here not include the .xml
extension. Omitting filename works as for *ReadXML*.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Any"),
      GOUT_DESC("Workspace variable to be read."),
      IN("file_index"),
      GIN("filename", "digits"),
      GIN_TYPE("String", "Index"),
      GIN_DEFAULT("", "0"),
      GIN_DESC(
          "File name. See above.",
          "Equalize the widths of all numbers by padding with zeros as necessary. "
          "0 means no padding (default)."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("Reduce"),
      DESCRIPTION(R"--(Reduces a larger class to a smaller class of same size.

The Reduce command reduces all \"1\"-dimensions to nil.  Examples:

1) 1 Vector can be reduced to a Numeric
2) 2x1 Matrix can be reduced to 2 Vector
3) 1x3x1 Tensor3 can be reduced to 3 Vector
4) 1x1x1x1 Tensor4 can be reduced to a Numeric
5) 3x1x4x1x5 Tensor5 can only be reduced to 3x4x5 Tensor3
6) 1x1x1x1x2x3 Tensor6 can be reduced to 2x3 Matrix
7) 2x3x4x5x6x7x1 Tensor7 can be reduced to 2x3x4x5x6x7 Tensor6

And so on
)--"),
      AUTHORS("Oliver Lemke", "Richard Larsson"),
      OUT(),
      GOUT("o"),
      GOUT_TYPE("Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric,"
                "Vector, Vector, Vector, Vector, Vector, Vector,"
                "Matrix, Matrix, Matrix, Matrix, Matrix,"
                "Tensor3, Tensor3, Tensor3, Tensor3,"
                "Tensor4, Tensor4, Tensor4,"
                "Tensor5, Tensor5,"
                "Tensor6"),
      GOUT_DESC("Reduced form of input."),
      IN(),
      GIN("i"),
      GIN_TYPE("Vector, Matrix, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7,"
               "Matrix, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7,"
               "Tensor3, Tensor4, Tensor5, Tensor6, Tensor7,"
               "Tensor4, Tensor5, Tensor6, Tensor7,"
               "Tensor5, Tensor6, Tensor7,"
               "Tensor6, Tensor7,"
               "Tensor7"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Over-dimensioned input"),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false)));

  md_data_raw.push_back(create_mdrecord(
      NAME("surface_fieldEarth"),
      DESCRIPTION(R"--(Earth reference ellipsoids.

The reference ellipsoid (``refellipsoid``) is set to model the Earth,
following different models. The options are:

- \"Sphere\":
    A spherical Earth. The radius is set following
    the value set for the Earth radius in constants.cc.
- \"WGS84\":
    The reference ellipsoid used by the GPS system.
    Should be the standard choice for a non-spherical Earth.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("model"),
      GIN_TYPE("String"),
      GIN_DEFAULT("Sphere"),
      GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(
      create_mdrecord(NAME("surface_fieldGanymede"),
               DESCRIPTION(R"--(Ganymede reference ellipsoids.

From Wikipedia
)--"),
               AUTHORS("Takayoshi Yamada"),
               OUT("surface_field"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN("model"),
               GIN_TYPE("String"),
               GIN_DEFAULT("Sphere"),
               GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("surface_fieldEuropa"),
      DESCRIPTION(R"--(Io reference ellipsoids.

The reference ellipsoid (``refellipsoid``) is set to model Io,
folowing different models. The options are:

- \"Sphere\": A spherical planetesimal. The radius is taken from report of the IAU/IAG Working Group.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("surface_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("model"),
      GIN_TYPE("String"),
      GIN_DEFAULT("Sphere"),
      GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("surface_fieldIo"),
      DESCRIPTION(R"--(Io reference ellipsoids.

The reference ellipsoid (``refellipsoid``) is set to model Io,
folowing different models. The options are:

- \"Sphere\": A spherical planetesimal. The radius is taken from report of the IAU/IAG Working Group.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("surface_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("model"),
      GIN_TYPE("String"),
      GIN_DEFAULT("Sphere"),
      GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("surface_fieldJupiter"),
      DESCRIPTION(R"--(Jupiter reference ellipsoids.

The reference ellipsoid (``refellipsoid``) is set to model Jupiter,
folowing different models. The options are:

- \"Sphere\": A spherical planet. The radius is taken from a report of the IAU/IAG Working Group.
- \"Ellipsoid\": A reference ellipsoid with parameters taken from a report of the IAU/IAG Working Group.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("model"),
      GIN_TYPE("String"),
      GIN_DEFAULT("Sphere"),
      GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("surface_fieldMars"),
      DESCRIPTION(R"--(Mars reference ellipsoids.

The reference ellipsoid (``refellipsoid``) is set to model Mars,
folowing different models. The options are:

- \"Sphere\": A spherical planet. The radius is taken from a report of the IAU/IAG Working Group.
- \"Ellipsoid\": A reference ellipsoid with parameters taken from a report of the IAU/IAG Working Group.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("model"),
      GIN_TYPE("String"),
      GIN_DEFAULT("Sphere"),
      GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("surface_fieldMoon"),
      DESCRIPTION(R"--(Moon reference ellipsoids.

The reference ellipsoid (``refellipsoid``) is set to model Moon,
folowing different models. The options are:

- \"Sphere\":
    A spherical planet. The radius is taken from a
    report of the IAU/IAG Working Group.

- \"Ellipsoid\":
    A reference ellipsoid with parameters taken from
    Wikepedia (see code for details). The IAU/IAG working group
    defines the Moon ellipsoid to be a sphere.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("model"),
      GIN_TYPE("String"),
      GIN_DEFAULT("Sphere"),
      GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("surface_fieldInit"),
      DESCRIPTION(R"--(Manual setting of the reference ellipsoid.

The two values of ``refellipsoid`` can here be set manually. The two
arguments correspond directly to first and second element of
``refellipsoid``.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("a", "b"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Average or equatorial radius.", "Average or polar radius.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("surface_fieldVenus"),
      DESCRIPTION(R"--(Venus reference ellipsoids.

The reference ellipsoid (``refellipsoid``) is set to model Venus,
folowing different models. The options are:

- \"Sphere\":
      A spherical planet. The radius is taken from a
      report of the IAU/IAG Working Group.

According to the report used above, the Venus ellipsoid lacks
eccentricity and no further models should be required.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("model"),
      GIN_TYPE("String"),
      GIN_DEFAULT("Sphere"),
      GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("refr_index_airFreeElectrons"),
      DESCRIPTION(R"--(Microwave refractive index due to free electrons.

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
(and there exist a corresponding \"vmr\"-value). This demand is
removed if ``demand_vmr_value`` is set to 0, but use this option
with care.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("refr_index_air", "refr_index_air_group"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("refr_index_air",
         "refr_index_air_group",
         "f_grid",
         "abs_species",
         "rtp_vmr"),
      GIN("demand_vmr_value"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Flag to control if it is demanded that free electrons are "
               "in *abs_species*. Default is that this is demanded.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("refr_index_airInfraredEarth"),
      DESCRIPTION(R"--(Calculates the IR refractive index due to gases in the
Earth's atmosphere.

Only refractivity of dry air is considered. The formula used is
contributed by Michael Hoepfner, Forschungszentrum Karlsruhe.

The refractivity of dry air is added to *refr_index_air*. To obtain
the complete value, *refr_index_air* should be set to 1 before
calling this WSM. This applies also to *refr_index_air_group*.

The expression used is non-dispersive. Hence, *refr_index_air* and
*refr_index_air_group* are identical.
)--"),
      AUTHORS("Mattias Ekstrom"),
      OUT("refr_index_air", "refr_index_air_group"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("refr_index_air",
         "refr_index_air_group",
         "rtp_pressure",
         "rtp_temperature"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("refr_index_airMicrowavesEarth"),
      DESCRIPTION(R"--(Microwave refractive index in Earth's atmosphere.

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
as \"Present study\". Note that in ARTS Pa is used for pressure
and k1, k2 and k3 must be adjusted accordingly.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("refr_index_air", "refr_index_air_group"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("refr_index_air",
         "refr_index_air_group",
         "rtp_pressure",
         "rtp_temperature",
         "rtp_vmr",
         "abs_species"),
      GIN("k1", "k2", "k3"),
      GIN_TYPE("Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT("77.6e-8", "70.4e-8", "3.739e-3"),
      GIN_DESC("Coefficient a, see above",
               "Coefficient b, see above",
               "Coefficient c, see above")));

  md_data_raw.push_back(create_mdrecord(
      NAME("refr_index_airMicrowavesGeneral"),
      DESCRIPTION(R"--(Microwave refractive index due to gases in planetary atmospheres.

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
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("refr_index_air", "refr_index_air_group"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("refr_index_air",
         "refr_index_air_group",
         "rtp_pressure",
         "rtp_temperature",
         "rtp_vmr",
         "abs_species"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalDefClose"),
      DESCRIPTION(R"--(Closes the definition of retrieval quantities and correlations and
prepares related WSVs for the retrieval.

This function calls jacobianClose and checks that the corvariance matrices
are consistent with the Jacobian.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("jacobian_do", "jacobian_agenda", "retrieval_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_agenda", "covmat_sx", "jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalAddAbsSpecies"),
      DESCRIPTION(R"--(Adds an absorption species to the retrieval quantities.

Similar to *jacobianAddAbsSpecies* but also sets the corresponding block in
*covmat_sx* to the matrices provided in *covmat_block* and *covmat_inv_block*.
The dimensions of *covmat_block* are required to agree with the dimensions of the
retrieval grid.

*covmat_inv_block* must be either empty or the same dimension as *covmat_block*.
If provided, this matrix will be used as the inverse for the covariance matrix block
and numerical inversion of this block is thus avoided. Note, however, that this is
only effective if this block is uncorrelated with any other retrieval quantity.

For number and order of elements added to *x*, see *jacobianAddAbsSpecies*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "covmat_block",
         "covmat_inv_block"),
      GIN("g1", "g2", "g3", "species", "unit", "for_species_tag"),
      GIN_TYPE("Vector", "Vector", "Vector", "String", "String", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, "rel", "1"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "The species tag of the retrieval quantity.",
               "Retrieval unit. See above.",
               "Index-bool for acting on species tags or species."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalAddFreqShift"),
      DESCRIPTION(R"--(Same as *jacobianAddFreqShift* but also adds the correlation block
contained in *covmat_block* and *covmat_inv_block* to *covmat_sx*.

For number and order of elements added to *x*, see *jacobianAddFreqShift*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "covmat_block",
         "covmat_inv_block",
         "jacobian_quantities",
         "jacobian_agenda",
         "f_grid"),
      GIN("df"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("100e3"),
      GIN_DESC("Size of perturbation to apply.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalAddFreqStretch"),
      DESCRIPTION(R"--(Same as *jacobianAddFreqShift* but also adds the correlation block
contained in *covmat_block* and *covmat_inv_block* to *covmat_sx*.

For number and order of elements added to *x*, see *jacobianAddFreqStretch*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda",
         "f_grid",
         "covmat_block",
         "covmat_inv_block"),
      GIN("df"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("100e3"),
      GIN_DESC("Size of perturbation to apply.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalDefInit"),
      DESCRIPTION(R"--(Begin retrieval definition section.

This function initialises all variables required for defining
retrieval quantities and corresponding covariance matrices.
By default, Jacobian quantities should be added withing the.
retrieval definition section. If Jacobian quantities are
defined separately ``initialize_jacobian`` must be set to 0,
otherwise the quantities will be discarded.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_se",
          "covmat_sx",
          "covmat_block",
          "covmat_inv_block",
          "jacobian_quantities",
          "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("initialize_jacobian"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Flag whether or not to (re)initialize Jacobian-related "
               "quantities. Set to 0 if Jacobian is already defined."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true),
      PASSWSVNAMES(false)));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalAddCatalogParameter"),
      DESCRIPTION(R"--(Similar to *jacobianAddBasicCatalogParameter* but also adds a corresponding
block to *covmat_sx* with the given ``var`` as variance value.

For number and order of elements added to *x*,
see *jacobianAddBasicCatalogParameter*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GIN("catalog_identity", "catalog_parameter", "var"),
      GIN_TYPE("QuantumIdentifier", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("The catalog line matching information.",
               "The catalog parameter of the retrieval quantity.",
               "The variance of the catalog parameter.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalAddCatalogParameters"),
      DESCRIPTION(R"--(Same as *jacobianAddBasicCatalogParameters* but also adds a new
block to *covmat_sx* using the matrices in *covmat_block* and
*covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*,
see *jacobianAddBasicCatalogParameters*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "covmat_block",
         "covmat_inv_block"),
      GIN("catalog_identities", "catalog_parameters"),
      GIN_TYPE("ArrayOfQuantumIdentifier", "ArrayOfString"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The catalog line matching informations.",
               "The catalog parameters of the retrieval quantity.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalAddMagField"),
      DESCRIPTION(R"--(Same as *jacobianAddMagField* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddMagField*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "covmat_block",
         "covmat_inv_block"),
      GIN("g1", "g2", "g3", "component", "dB"),
      GIN_TYPE("Vector", "Vector", "Vector", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "v", "1.0e-7"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Magnetic field component to retrieve",
               "Magnetic field perturbation")));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalAddPointingZa"),
      DESCRIPTION(R"--(Same as *jacobianAddPointingZa* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddPointingZa*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "covmat_block",
         "covmat_inv_block",
         "sensor_pos",
         "sensor_time"),
      GIN("poly_order", "calcmode", "dza"),
      GIN_TYPE("Index", "String", "Numeric"),
      GIN_DEFAULT("0", "recalc", "0.01"),
      GIN_DESC("Order of polynomial to describe the time variation of "
               "pointing off-sets.",
               "Calculation method. See above",
               "Size of perturbation to apply (when applicable).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalAddPolyfit"),
      DESCRIPTION(R"--(Same as *jacobianAddPolyfit* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddPolyfit*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "covmat_block",
         "covmat_inv_block",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid",
         "sensor_pos"),
      GIN("poly_order",
          "no_pol_variation",
          "no_los_variation",
          "no_mblock_variation"),
      GIN_TYPE("Index", "Index", "Index", "Index"),
      GIN_DEFAULT(NODEF, "0", "0", "0"),
      GIN_DESC("Polynomial order to use for the fit.",
               "Set to 1 if the baseline off-set is the same for all "
               "Stokes components.",
               "Set to 1 if the baseline off-set is the same for all "
               "line-of-sights (inside each measurement block).",
               "Set to 1 if the baseline off-set is the same for all "
               "measurement blocks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalAddScatSpecies"),
      DESCRIPTION(R"--(Same as *jacobianAddPolyfit* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddScatSpecies*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "covmat_block",
         "covmat_inv_block"),
      GIN("g1", "g2", "g3", "species", "quantity"),
      GIN_TYPE("Vector", "Vector", "Vector", "String", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, NODEF),
      GIN_DESC(
          "Pressure retrieval grid.",
          "Latitude retrieval grid.",
          "Longitude retreival grid.",
          "Name of scattering species, must match one element in *scat_species*.",
          "Retrieval quantity, e.g. \"IWC\"."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalAddSinefit"),
      DESCRIPTION(R"--(Same as *jacobianAddSinefit* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddSinefit*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "covmat_block",
         "covmat_inv_block",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid",
         "sensor_pos"),
      GIN("period_lengths",
          "no_pol_variation",
          "no_los_variation",
          "no_mblock_variation"),
      GIN_TYPE("Vector", "Index", "Index", "Index"),
      GIN_DEFAULT(NODEF, "0", "0", "0"),
      GIN_DESC("Period lengths of the fit.",
               "Set to 1 if the baseline off-set is the same for all "
               "Stokes components.",
               "Set to 1 if the baseline off-set is the same for all "
               "line-of-sights (inside each measurement block).",
               "Set to 1 if the baseline off-set is the same for all "
               "measurement blocks.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalAddSpecialSpecies"),
      DESCRIPTION(R"--(Same as *jacobianAddSpecialSpecies* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddSpecialSpecies*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "covmat_block",
         "covmat_inv_block"),
      GIN("g1", "g2", "g3", "species"),
      GIN_TYPE("Vector", "Vector", "Vector", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "The species of the retrieval quantity."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalAddSurfaceQuantity"),
      DESCRIPTION(R"--(Same as *jacobianAddSurfaceQuantity* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddSurfaceQuantity*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "covmat_block",
         "covmat_inv_block"),
      GIN("g1", "g2", "quantity"),
      GIN_TYPE("Vector", "Vector", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Retrieval quantity, e.g. \"Wind speed\".")));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalAddTemperature"),
      DESCRIPTION(R"--(Same as *jacobianAddTemperature* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddTemperature*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "covmat_block",
         "covmat_inv_block"),
      GIN("g1", "g2", "g3", "hse"),
      GIN_TYPE("Vector", "Vector", "Vector", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "on"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Flag to assume HSE or not (\"on\" or \"off\").")));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalAddWind"),
      DESCRIPTION(R"--(Same as *jacobianAddWind* but also adds a new block to *covmat_sx*
using the matrices in *covmat_block* and *covmat_inv_block*.

If *covmat_inv_block* is non-empty, it is used as inverse for the added block
which avoids its numerical computation.

For number and order of elements added to *x*, see *jacobianAddWind*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "covmat_block",
         "covmat_inv_block"),
      GIN("g1", "g2", "g3", "component", "dfrequency"),
      GIN_TYPE("Vector", "Vector", "Vector", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "v", "0.1"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retrieval grid.",
               "Wind component to retrieve",
               "This is the frequency perturbation")));

  md_data_raw.push_back(create_mdrecord(
      NAME("retrievalErrorsExtract"),
      DESCRIPTION(R"--(Extract retrieval error from covariance matrices.

Extracts the error estimates for the retrieved quantities from the covariance
matrices for the error due to measurement noise *covmat_so* and the error due
to limited resolution of the observation system *covmat_ss* and stores them in
the vectors *retrieval_eo* and *retrieval_ss*, respectively.

To etract these errors, first the convariance matrices of which the errors 
should be extracted have to be computed using the WSMs *covmat_soCalc*
and *covmat_ssCalc* or set to be empty in order to be ignored. Note, however,
that this will also set the corresponding error vector to be empty.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("retrieval_eo", "retrieval_ss"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_so", "covmat_ss"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("RT4Calc"),
      DESCRIPTION(R"--(Interface to the PolRadTran RT4 scattering solver (by F. Evans).

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
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("cloudbox_field", "za_grid", "aa_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmfields_checked",
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
         
         "surface_field"),
      GIN("nstreams",
          "pfct_method",
          "quad_type",
          "add_straight_angles",
          "pfct_aa_grid_size",
          "auto_inc_nstreams",
          "robust",
          "za_interp_order",
          "cos_za_interp",
          "max_delta_tau"),
      GIN_TYPE("Index",
               "String",
               "String",
               "Index",
               "Index",
               "Index",
               "Index",
               "Index",
               "Index",
               "Numeric"),
      GIN_DEFAULT("16", "median", "D", "1", "19", "0", "0", "1", "0", "1e-6"),
      GIN_DESC("Number of polar angle directions (streams) in RT4"
               " solution (must be an even number).",
               "Flag which method to apply to derive phase function (for"
               " available options see above).",
               "Flag which quadrature to apply in RT4 solution (for"
               " available options see above).",
               "Flag whether to include nadir and zenith as explicit"
               " directions (only effective for quad_type G and D).",
               "Number of azimuthal angle grid points to consider in"
               " Fourier series decomposition of scattering matrix (only"
               " applied for randomly oriented scattering elements)",
               "Flag whether to internally increase nstreams (individually"
               " per frequency) if norm of (bulk) scattering matrix is not"
               " preserved properly. If 0, no adaptation is done. Else"
               " ``auto_inc_nstreams`` gives the maximum number of streams to"
               " increase to. Note that the output *cloudbox_field* remains"
               " with angular dimension of ``nstreams``, only the internal"
               " solution is adapted (and then interpolated to the"
               " lower-resolution output angular grid).",
               "For ``auto_inc_nstreams``>0, flag whether to not fail even if"
               " scattering matrix norm is not preserved when maximum stream"
               " number is reached. Internal RT4 calculations is then"
               " performed with nstreams=``auto_inc_nstreams``.",
               "For ``auto_inc_nstreams``>0, polar angle interpolation order"
               " for interpolation from internal increased stream to"
               " originally requested nstreams-ifield.",
               "For ``auto_inc_nstreams``>0, flag whether to do polar angle"
               " interpolation in cosine (='mu') space.",
               "Maximum optical depth of infinitesimal layer (where single"
               " scattering approximation is assumed to apply).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("RT4CalcWithRT4Surface"),
      DESCRIPTION(R"--(As RT4Calc except for using RT4's proprietary surface type handling.

This WSM is only indented for testing purposes.

The following surface type/property methods are available and
require the the following input:

- 'L'ambertian: *surface_scalar_reflectivity*, *surface_skin_t*
- 'F'resnel: *surface_complex_refr_index*, *surface_skin_t*
- 'S'pecular: *surface_reflectivity*, *surface_skin_t*

'L' and 'F' use proprietary RT4 methods, 'S' uses RT4's Fresnel
methods modified to behave similar to ARTS'
*surfaceFlatReflectivity*.
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("cloudbox_field", "za_grid", "aa_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmfields_checked",
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
         "surface_complex_refr_index"),
      GIN("nstreams",
          "pfct_method",
          "ground_type",
          "quad_type",
          "add_straight_angles",
          "pfct_aa_grid_size",
          "auto_inc_nstreams",
          "robust",
          "za_interp_order",
          "cos_za_interp",
          "max_delta_tau"),
      GIN_TYPE("Index",
               "String",
               "String",
               "String",
               "Index",
               "Index",
               "Index",
               "Index",
               "Index",
               "Index",
               "Numeric"),
      GIN_DEFAULT(
          "16", "median", "A", "D", "1", "19", "0", "0", "1", "0", "1e-6"),
      GIN_DESC("Number of polar angle directions (streams) in RT4"
               " solution (must be an even number).",
               "Flag which method to apply to derive phase function (for"
               " available options see above).",
               "Flag which surface type/surface property method to use"
               " (for available options see above).",
               "Flag which quadrature to apply in RT4 solution (for"
               " available options see above).",
               "Flag whether to include nadir and zenith as explicit"
               " directions (only effective for quad_type G and D).",
               "Number of azimuthal angle grid points to consider in"
               " Fourier series decomposition of scattering matrix (only"
               " applied for randomly oriented scattering elements)",
               "Flag whether to internally increase nstreams (individually"
               " per frequency) if norm of (bulk) scattering matrix is not"
               " preserved properly. If 0, no adaptation is done. Else"
               " ``auto_inc_nstreams`` gives the maximum number of streams to"
               " increase to.",
               "For ``auto_inc_nstreams``>0, flag whether to not fail even if"
               " scattering matrix norm is not preserved when maximum stream"
               " number is reached. Internal RT4 calculations is then"
               " performed with nstreams=``auto_inc_nstreams``.",
               "For ``auto_inc_nstreams``>0, polar angle interpolation order"
               " for interpolation from internal increased stream to"
               " originally requested nstreams-ifield.",
               "For ``auto_inc_nstreams``>0, flag whether to do polar angle"
               " interpolation in cosine (='mu') space.",
               "Maximum optical depth of infinitesimal layer (where single"
               " scattering approximation is assumed to apply).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("RT4Test"),
      DESCRIPTION(R"--(RT4 validation test.

Executes test case testc shipped with PolRadTran/RT4 code (but uses
data files converted to arts-xml). Output written to (xml-)file.
)--"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("out_rad"),
      GOUT_TYPE("Tensor4"),
      GOUT_DESC("RT4 testc calculation results."),
      IN(),
      GIN("datapath"),
      GIN_TYPE("String"),
      GIN_DEFAULT("artscomponents/polradtran/testdata/"),
      GIN_DESC("Folder containing arts-xml converted test case input data.")));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("rte_losGeometricFromRtePosToRtePos2"),
      DESCRIPTION(R"--(The geometric line-of-sight between two points.

Old! Use *rte_losGeometricToPosition* ZZZ

The method sets *rte_los* to the line-of-sight, at *rte_pos*,
that matches the geometrical propagation path between *rte_pos*
and *rte_pos2*.

The standard case should be that *rte_pos2* corresponds to a
transmitter, and *rte_pos* to the receiver/sensor.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("rte_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         "lat_grid",
         "lon_grid",
         "rte_pos",
         "rte_pos2"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("rte_losGeometricToPosition"),
      DESCRIPTION(R"--(The geometric line-of-sight between two points.

The line-of-sight angles from *rte_pos* to ``target_pos`` are calculated
ignoring refraction. This can be done analytically. The angles are set
without any consideration of the surface. The corresponding propagation
path can thus end with a surface intersection.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("rte_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("surface_field",
         "rte_pos"),
      GIN("target_pos"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The atmospheric position that *rte_los* shall match.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("rte_losRefractedToPosition"),
      DESCRIPTION(R"--(The refracted line-of-sight and propagation path between two positions.

General remarks: The method assumes geometrical optics. It can fail for
various reasons, especially if there are sharp or multiple gradients in
refractive index. In multi-path situations, the method finds (in best
case) only a single solution. Default is to consider an intersection
with the surface as a failure and give an error about it. To instead
return an empty ppath when the target appears to be behind the horizon,
set GIN ``robust`` to 1.

The line-of-sight connecting the two points cannot be determined
analytically and a search algorithm is needed. The algorithm options are:

\"basic\":
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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("rte_los", "ppath"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("refr_index_air_ZZZ_agenda",
         "ppath_lstep",
         "ppath_lraytrace",
         "surface_field",
         "surface_search_accuracy",
         "rte_pos"),
      GIN("target_pos",
          "target_dl",
          "algorithm",
          "max_iterations",
          "robust",
          "z_toa",
          "do_horizontal_gradients",
          "do_twosided_perturb"),
      GIN_TYPE("Vector",
               "Numeric",
               "String",
               "Index",
               "Index",
               "Numeric",
               "Index",
               "Index"),
      GIN_DEFAULT(NODEF, NODEF, "basic", "10", "0", NODEF, "0", "0"),
      GIN_DESC("The atmospheric position that *ppath* shall reach.",
               "The end point of *ppath* shall be inside this distance "
               "from ``target_pos`` (deviation can be in any direction).",
               "Search algorithm to use.",
               "Max number of iterations before giving up.",
               "Set to 1 to not give errors, but return empty *ppath* "
               "when a path can not be established.",
               "Top-of-the-atmosphere altitude.",
               "Consider horisontal gradients of refractive index.",
               "Perform double-sided perturbations when calculating "
               "refractive index gradients.")));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("rte_losReverse"),
      DESCRIPTION(R"--(Reverses the direction in *rte_los*.

The method updates *rte_los* to have angles of the reversed
direction.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("rte_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("rte_los"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord(NAME("rte_losSet"),
               DESCRIPTION(R"--(Sets *rte_los* to the given angles.
)--"),
               AUTHORS("Patrick Eriksson"),
               OUT("rte_los"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN("za", "aa"),
               GIN_TYPE("Numeric", "Numeric"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Zenith angle [0, 180]",
                        "Azimuth angle [-180, 180]")));

  md_data_raw.push_back(create_mdrecord(
      NAME("rte_posSet"),
      DESCRIPTION(R"--(Sets *rte_pos* to the given coordinates.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("rte_pos"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("z", "lat", "lon"),
      GIN_TYPE("Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Altitude",
               "Latitude [-90, 90]",
               "Longitude [-180, 360]")));

  md_data_raw.push_back(create_mdrecord(
      NAME("rte_pos_losBackwardToAltitude"),
      DESCRIPTION(R"--(Moves *rte_pos* and *rte_los* backwards to the target altitude.

The method gives the *rte_pos* and *rte_los* at the target altitude
to reach the original *rte_pos* and *rte_los* with a geometrical ppath.
That is, the movement is backwards in terms of viewing direction.

If the original *rte_los* is reversed with respect to the line-of-sight
direction, then set the GIN los_reversed to 1. One such case is that
if *rte_los* represents surface incidence angles, i.e. holds the
zenith and nadir angle towards the sensor.

There is also *sensor_pos_losBackwardToAltitude*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("rte_pos", "rte_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("rte_pos", "rte_los", "surface_field"),
      GIN("altitude", "los_is_reversed"),
      GIN_TYPE("Numeric", "Index"),
      GIN_DEFAULT(NODEF, "0"),
      GIN_DESC("Target altitude.",
               "Set to 1 if *rte_los* is valid for the reversed direction.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("rte_pos_losEndOfPpath"),
      DESCRIPTION(R"--(Sets *rte_pos* and *rte_los* to values for last point in *ppath*.

For example, if the propagation path intersects with the surface,
this method gives you the position and angle of *ppath* at the
surface.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("rte_pos", "rte_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("rte_pos_losForwardToAltitude"),
      DESCRIPTION(R"--(Moves *rte_pos* and *rte_los* forward to the target altitude.

The method gives the *rte_pos* and *rte_los* at the target altitude
when forward-propagating the original *rte_pos* and *rte_los*
geometrically.

There is also *sensor_pos_losForwardToAltitude*. The WSM
*IntersectionGeometricAltitude* performs the same operation
with *sensor_pos* and *sensor_los* as input.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("rte_pos", "rte_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("rte_pos", "rte_los", "surface_field"),
      GIN("altitude"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Target altitude.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ScatElementsPndAndScatAdd"),
      DESCRIPTION(R"--(Adds single scattering data and particle number density for
individual scattering elements.

The methods reads the specified files and appends the obtained data
to *scat_data* and *pnd_field_raw*. Scattering data is appended to
the current last existing scattering species in *scat_data*.
)--"),
      AUTHORS("Claudia Emde, Jana Mendrok"),
      OUT("scat_data_raw", "pnd_field_raw"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data_raw", "pnd_field_raw"),
      GIN("scat_data_files", "pnd_field_files"),
      GIN_TYPE("ArrayOfString", "ArrayOfString"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("List of names of single scattering data files.",
               "List of names of the corresponding pnd_field files.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ScatElementsSelect"),
      DESCRIPTION(R"--(Allows to limit considered scattering elements according to size.

Scattering elements of a specified scattering species are removed
from *scat_data_raw* and *scat_meta*, i.e. removed from further
calculations, if their particle size exceeds the specified limits.
Specification of the scattering species is done by name matching the
scattering species name part of *scat_species* tag.
As size parameter, all size parameters reported by the meta data
can be used (see *scat_meta_single* for offered parameters and
their naming).
)--"),
      AUTHORS("Daniel Kreyling, Oliver Lemke, Jana Mendrok"),
      OUT("scat_data_raw", "scat_meta"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data_raw", "scat_meta", "scat_species"),
      GIN("species", "sizeparam", "sizemin", "sizemax", "tolerance", "delim"),
      GIN_TYPE("String", "String", "Numeric", "Numeric", "Numeric", "String"),
      GIN_DEFAULT(NODEF, NODEF, "0.", "-1.", "1e-6", "-"),
      GIN_DESC("Species on which to apply size selection.",
               "Size parameter to apply for size selection.",
               "Minimum size [m] of the scattering elements to consider",
               "Maximum size [m] of the scattering elements to consider (if "
               "negative, no max. limitation is applied).",
               "Relative numerical tolerance of size limit values.",
               "Delimiter string of *scat_species* elements.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ScatElementsToabs_speciesAdd"),
      DESCRIPTION(R"--(Appends scattering elements to *abs_species* including reading
single scattering data and corresponding pnd field.

The methods reads the specified single scattering and pnd_field
data of individual scattering elements and appends the obtained data
to *scat_data* (appending to its last scattering species) and
``vmr_field_raw``. Per scattering element, it also appends one
instance of species 'particles' to *abs_species*.
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("scat_data_raw",
          "atm_field",
          "abs_species",
          "propmat_clearsky_agenda_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data_raw",
         "atm_field",
         "abs_species",
         "propmat_clearsky_agenda_checked",
         "f_grid"),
      GIN("scat_data_files", "pnd_field_files"),
      GIN_TYPE("ArrayOfString", "ArrayOfString"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("List of names of single scattering data files.",
               "List of names of the corresponding pnd_field files.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ScatSpeciesExtendTemperature"),
      DESCRIPTION(R"--(Extends valid temperature range of single scattering data.

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
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("scat_data_raw"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data_raw", "scat_species"),
      GIN("species", "scat_species_delim", "T_low", "T_high"),
      GIN_TYPE("String", "String", "Numeric", "Numeric"),
      GIN_DEFAULT("", "-", "-1.", "-1."),
      GIN_DESC(
          "Scattering species to act on (see WSM description for details).",
          "Delimiter string of *scat_species* elements.",
          "Temperature grid extension point at low temperature limit.",
          "Temperature grid extension point at high temperature limit.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ScatSpeciesInit"),
      DESCRIPTION(R"--(Initializes the scattering species related data variables.

This method initializes the *scat_species* WSV, the variables that
will hold the raw optical properties and the raw particle number
distributions of the scattering elements (*scat_data_raw* and
*pnd_field_raw*, respectively) as well as the one holding the meta
information about the scattering elements (*scat_meta*).

This method has to be executed before WSM reading/adding to the
said variable, e.g. before *ScatSpeciesPndAndScatAdd*.
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("scat_species",
          "scat_data_raw",
          "scat_meta",
          "scat_data_checked",
          "pnd_field_raw"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("ScatSpeciesPndAndScatAdd"),
      DESCRIPTION(R"--(Adds single scattering data and particle number densities for one
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
)--"),
      AUTHORS("Claudia Emde, Jana Mendrok"),
      OUT("scat_data_raw", "pnd_field_raw"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data_raw", "pnd_field_raw"),
      GIN("scat_data_files", "pnd_fieldarray_file"),
      GIN_TYPE("ArrayOfString", "String"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC(
          "Array of names of files containing the single scattering data.",
          "Name of file holding the corresponding array of pnd_field data.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ScatSpeciesScatAndMetaRead"),
      DESCRIPTION(R"--(Reads single scattering data and scattering meta data for one
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
)--"),
      AUTHORS("Daniel Kreyling, Oliver Lemke, Jana Mendrok"),
      OUT("scat_data_raw", "scat_meta"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data_raw", "scat_meta"),
      GIN("scat_data_files"),
      GIN_TYPE("ArrayOfString"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Array of single scattering data file names.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("scat_data_singleTmatrix"),
      DESCRIPTION(R"--(A basic interface to Mishchenko's T-matrix code linked to ARTS.

The method performs a T-matrix calculation for a single scattering
element, i.e. a combination of particle shape, size, aspect ratio
and orientation.

Particle shape (``shape``) has two options::

  \"spheroidal\" and \"cylindrical\"

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

  \"totally_random\" and \"azimuthally_random\"

For totally randomly oriented particles, ``data_aa_grid`` is not taken
into account (but a Vector type container needs to be passed).

For further information on how aspect ratio and the different shapes
and orientations are defined, see the documentation of the T-matrix
code found http:

Regarding ``ndgs``, we refer to the this comment from the documentation:
   \"Parameter controlling the number of division points
   in computing integrals over the particle surface.
   For compact particles, the recommended value is 2.
   For highly aspherical particles larger values (3, 4,...)
   may be necessary to obtain convergence.
   The code does not check convergence over this parameter.
   Therefore, control comparisons of results obtained with
   different NDGS-values are recommended.\
)--"),
      AUTHORS("Johan Strandgren", "Patrick Eriksson"),
      OUT("scat_data_single", "scat_meta_single"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("complex_refr_index"),
      GIN("shape",
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
          "quiet"),
      GIN_TYPE("String",
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
               "Index"),
      GIN_DEFAULT(NODEF,
                  NODEF,
                  NODEF,
                  "NaN",
                  NODEF,
                  NODEF,
                  NODEF,
                  NODEF,
                  "[]",
                  "0.001",
                  "Set by user, unknown source.",
                  "2",
                  "0",
                  "1"),
      GIN_DESC("Particle shape. Options listed above.",
               "Particle volume equivalent diameter [m]. See defintion above.",
               "Particle aspect ratio.",
               "Particle mass. This information is just included in the meta"
               " data, and does not affect the T-matrix calculations.",
               "Particle type/orientation. Options listed above.",
               "Frequency grid of the scattering data to be calculated.",
               "Temperature grid of the scattering data to be calculated.",
               "Zenith angle grid of the scattering data to be calculated.",
               "Azimuth angle grid of the scattering data to be calculated.",
               "Accuracy of the computations.",
               "String describing the source of *complex_refr_index*, for"
               " inclusion in meta data.",
               "See above. So far only applied for random orientation.",
               "Continue even if individual T-matrix calculations fail. "
               "Respective scattering element data will be NAN.",
               "Suppress print output from tmatrix fortran code.")));

  /*
   md_data_raw.push_back
    ( create_mdrecord
    ( NAME( "scat_metaAddTmatrix" ),
      DESCRIPTION
      (
       "This method adds scattering element meta data to the workspace variable\n"
       "*scat_meta*.\n"
       "\n"
       "One set of meta data is created and added to the array for each\n"
       "combination of maximum diameter and aspect ratio in the GINs\n"
       "diameter_max_grid and aspect_ratio_grid. The size of *scat_meta*\n"
       "and hence the usage has been extended. For that reason, a short summary\n"
       "below tells which input parameters are required for certain further\n"
       "calculations.\n"
       "\n"
       "String[description]\t\tNot used for any particular calculations\n"
       "String[material]\t\tUsed for PND calculations\n"
       "String[shape]\t\t\tUsed for scattering and PND calculations\n"
       "Numeric[ptype]\t\tUsed for scattering calculations\n"
       "Numeric[density]\t\tUsed for PND calculations\n"
       "Vector[diameter_max_grid]\t\tUsed for both scattering and PND calculations\n"
       "Vector[aspect_ratio_grid]\t\tUsed for scattering calculations and PND calculations\n"
       "Vector[scat_f_grid]\t\tUsed for scattering calculations\n"
       "Vector[scat_T_grid]\t\tUsed for scattering calculations\n"
       "Tensor3[complex_refr_index]\tUsed for scattering calculations\n"
      ),
      AUTHORS( "Johan Strandgren" ),
      OUT("scat_meta"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN( "scat_meta", "complex_refr_index" ),
      GIN( "description", "material", "shape", "ptype", "density",
           "aspect_ratio_grid", "diameter_max_grid", "scat_f_grid", "scat_T_grid" ),
      GIN_TYPE( "String", "String", "String", "String", "Numeric", "Vector",
           "Vector", "Vector", "Vector" ),
      GIN_DEFAULT( "", "undefined", NODEF, NODEF, "-999", NODEF, NODEF,
                   NODEF, NODEF ),
      GIN_DESC( "Particle description", "Water or Ice", "spheroidal or cylinder",
                "Particle Type: "totally_random" (20) or "azimuthally_random" (30)",
                "Particle mass density",
                "Particle aspect ratio vector",
                "Maximum diameter vector (diameter of a sphere that fully"
                "encloses the particle)",
                "Frequency grid vector",
                "Temperature grid vector" )
      ));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("scat_data_checkedCalc"),
      DESCRIPTION(R"--(Checks dimensions, grids and single scattering properties of all
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
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("scat_data_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data", "f_grid"),
      GIN("dfrel_threshold", "check_level", "sca_mat_threshold"),
      GIN_TYPE("Numeric", "String", "Numeric"),
      GIN_DEFAULT("0.1", "all", "5e-2"),
      GIN_DESC("Maximum relative frequency deviation between (single entry)"
               " scattering element f_grid values and the RT calculation's"
               " *f_grid*.",
               "See ``check_level`` in *scat_dataCheck*.",
               "See ``sca_mat_threshold`` in *scat_dataCheck*.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("scat_data_monoCalc"),
      DESCRIPTION(R"--(Interpolates *scat_data* by frequency to give *scat_data_mono*.
)--"),
      AUTHORS("Cory Davis"),
      OUT("scat_data_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data", "f_grid", "f_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("scat_data_monoExtract"),
      DESCRIPTION(R"--(Extracts data at *f_index* from *scat_data* to give *scat_data_mono*.
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("scat_data_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data", "f_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("scat_dataCalc"),
      DESCRIPTION(R"--(Prepares *scat_data* for the scattering solver.

Derives single scattering data for the frequencies given by
*f_grid* by interpolation from *scat_data_raw*. *f_grid* should be
the actual WSV *f_grid* or a single-element Vector.
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("scat_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data_raw", "f_grid"),
      GIN("interp_order"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Interpolation order.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("scat_dataCheck"),
      DESCRIPTION(R"--(Method for checking the validity and consistency of the single
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
)--"),
      AUTHORS("Claudia Emde", "Jana Mendrok"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data"),
      GIN("check_type", "sca_mat_threshold"),
      GIN_TYPE("String", "Numeric"),
      GIN_DEFAULT("all", "5e-2"),
      GIN_DESC("The level of checks to apply on scat_data ('sane' or 'all';"
               " see above).",
               "Threshold for allowed albedo deviation (see above).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("scat_dataReduceT"),
      DESCRIPTION(R"--(Reduces temperature dimension of single scattering to a single entry.

FIXME...
Derives single scattering data for the frequencies given by
*f_grid* by interpolation from *scat_data*. *f_grid* should be
the actual WSV *f_grid* or a single-element Vector.
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("scat_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data"),
      GIN("scat_index",
          "temperature",
          "interp_order",
          "phamat_only",
          "sca_mat_threshold"),
      GIN_TYPE("Index", "Numeric", "Index", "Index", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, "1", "1", "5e-2"),
      GIN_DESC("Apply on *scat_data* from scattering species of this index"
               " (0-based).",
               "Temperature to interpolate *scat_data* to.",
               "Interpolation order.",
               "Flag whether to apply temperture reduction on phase matrix"
               " data only (1) or on all single scattering properties (0).",
               "Threshold for allowed albedo deviation.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ScatSpeciesSizeMassInfo"),
      DESCRIPTION(R"--(Derives size and mass information for a scattering species.

This method assumes that the mass-size relationship can described
by *scat_species_a* and *scat_species_b*. See documentation of 
*scat_species_a* for details.

The quantity to be used as size descriptor is here denoted as x, and
is selected by setting ``x_unit``. The options are:

- ``\"dveq\"``: The size grid is set to scat_meta.diameter_volume_equ
- ``\"dmax\"``: The size grid is set to scat_meta.diameter_max
- ``\"area\"``: The size grid is set to scat_meta.diameter_area_equ_aerodynamical
- ``\"mass\"``: The size grid is set to scat_meta.mass

This selection determines *scat_species_x*.

The parameters *scat_species_a* and *scat_species_b* are determined by
a numeric fit between *scat_species_x* and corresponding masses in
*scat_meta*. This fit is performed over sizes inside the range
[x_fit_start,x_fit_end]. This range is allowed to be broader than
the coverage of *scat_species_x*. There must be at least two sizes
inside [x_fit_start,x_fit_end].
)--"),
      AUTHORS("Manfred Brath", "Jana Mendrok", "Patrick Eriksson"),
      OUT("scat_species_x", "scat_species_a", "scat_species_b"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_meta"),
      GIN("species_index", "x_unit", "x_fit_start", "x_fit_end", "do_only_x"),
      GIN_TYPE("Index", "String", "Numeric", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, "0", "1e9", "0"),
      GIN_DESC("Take data from scattering species of this index (0-based) in"
               " *scat_meta*.",
               "Unit for size grid, allowed options listed above.",
               "Smallest size to consider in fit to determine a and b.",
               "Largest size to consider in fit to determine a and b.",
               "A flag to deactivate calculation of a and b, to possibly "
               "save some time. The a and b parameters are then set to -1."
               "Default is to calculate a and b.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("particle_fieldCleanup"),
      DESCRIPTION(R"--(Removes unrealistically small or erroneous data from particle fields.

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
)--"),
      AUTHORS("Daniel Kreyling"),
      OUT(),
      GOUT("particle_field_out"),
      GOUT_TYPE("Tensor4"),
      GOUT_DESC("A particle property field, e.g. ``particle_bulkprop_field``"),
      IN(),
      GIN("particle_field_in", "threshold"),
      GIN_TYPE("Tensor4", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("A particle property field, e.g. ``particle_bulkprop_field``",
               "Threshold below which the ``particle_field`` values are set to"
               " zero.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Select"),
      DESCRIPTION(R"--(Method to select some elements from one array and copy them to
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
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("needles"),
      GOUT_TYPE((ARRAY_GROUPS + ", Vector, Matrix, Sparse").c_str()),
      GOUT_DESC("Selected elements. Must have the same variable type as "
                "haystack."),
      IN(),
      GIN("haystack", "needleindexes"),
      GIN_TYPE((ARRAY_GROUPS + ", Vector, Matrix, Sparse").c_str(), "ArrayOfIndex"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Variable to select from. May be the same variable as needles.",
               "The elements to select (zero based indexing, as always.)"),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_checkedCalc"),
      DESCRIPTION(R"--(Checks consistency of the sensor variables.

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
)--"),
      AUTHORS("Jana Mendrok"),
      OUT("sensor_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         
         "f_grid",
         "sensor_pos",
         "sensor_los",
         "transmitter_pos",
         "mblock_dlos",
         "sensor_response",
         "sensor_response_f",
         "sensor_response_pol",
         "sensor_response_dlos"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensorOff"),
      DESCRIPTION(R"--(Sets sensor WSVs to obtain monochromatic pencil beam values.

The variables are set as follows:
 - *mblock_dlos*        : One row with zero(s).
 - *sensor_response*        : As returned by *sensor_responseInit*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid",
          "sensor_response_pol_grid",
          "sensor_response_dlos_grid",
          "mblock_dlos"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN( "f_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_losAddLosAndDlos"),
      DESCRIPTION(R"--(Adds zenith and azimuth angles.

Adds up a line-of-sights (ref_los), with relative angle off-sets
(dlos).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("ref_los", "dlos"),
      GIN_TYPE("Vector", "Matrix"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Reference line-of-sight (a single los).",
               "Change in line-of-sight (can be multiple dlos).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_losGeometricToPosition"),
      DESCRIPTION(R"--(The geometric line-of-sight to a point.

Works as *rte_losGeometricToPosition*, but sets *sensor_los*.

This method handles the case of a single target position. For
multiple target positions, use: *sensor_losGeometricToPositions*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("surface_field",
         "sensor_pos"),
      GIN("target_pos"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The atmospheric position that *sensor_los* shall match.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_losGeometricToPositions"),
      DESCRIPTION(R"--(The geometric line-of-sight to multiple point.

Works as *rte_losGeometricToPosition*, but sets *sensor_los*. The
number of rows in *sensor_pos* and ``target_pos`` must be equal.

This method handles the case of mutiple target positions. For
a single target positions, use: *sensor_losGeometricToPosition*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("surface_field",
         "sensor_pos"),
      GIN("target_pos"),
      GIN_TYPE("Matrix"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The atmospheric positions that *sensor_los* shall match.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_losRefractedToPosition"),
      DESCRIPTION(R"--(The refracted line-of-sight to a point.

Works as *rte_losRefractedToPosition*, but sets *sensor_los*.

This method handles the case of a single target position. For
multiple target positions, use: *sensor_losRefractedToPositions*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("refr_index_air_ZZZ_agenda",
         "ppath_lstep",
         "ppath_lraytrace",
         "surface_field",
         "surface_search_accuracy",
         "sensor_pos"),
      GIN("target_pos",
          "target_dl",
          "algorithm",
          "max_iterations",
          "robust",
          "z_toa",
          "do_horizontal_gradients",
          "do_twosided_perturb"),
      GIN_TYPE("Vector",
               "Numeric",
               "String",
               "Index",
               "Index",
               "Numeric",
               "Index",
               "Index"),
      GIN_DEFAULT(NODEF, NODEF, "basic", "10", "0", NODEF, "0", "0"),
      GIN_DESC("The atmospheric position that *ppath* shall reach.",
               "The end point of *ppath* shall be inside this distance "
               "from ``target_pos`` (deviation can be in any direction).",
               "Search algorithm to use.",
               "Max number of iterations before giving up.",
               "Set to 1 to not give errors, but return empty *ppath* "
               "when a path can not be established.",
               "Top-of-the-atmosphere altitude.",
               "Consider horisontal gradients of refractive index.",
               "Perform double-sided perturbations when calculating "
               "refractive index gradients.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_losRefractedToPositions"),
      DESCRIPTION(R"--(The refracted line-of-sight to multiple points.

Works as *rte_losRefractedToPosition*, but sets *sensor_los*.

This method handles the case of multiple target positions. For
a single target position, use: *sensor_losRefractedToPosition*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("refr_index_air_ZZZ_agenda",
         "ppath_lstep",
         "ppath_lraytrace",
         "surface_field",
         "surface_search_accuracy",
         "sensor_pos"),
      GIN("target_pos",
          "target_dl",
          "algorithm",
          "max_iterations",
          "robust",
          "z_toa",
          "do_horizontal_gradients",
          "do_twosided_perturb"),
      GIN_TYPE("Matrix",
               "Numeric",
               "String",
               "Index",
               "Index",
               "Numeric",
               "Index",
               "Index"),
      GIN_DEFAULT(NODEF, NODEF, "basic", "10", "0", NODEF, "0", "0"),
      GIN_DESC("The atmospheric positions that *ppath* shall reach.",
               "The end point of *ppath* shall be inside this distance "
               "from ``target_pos`` (deviation can be in any direction).",
               "Search algorithm to use.",
               "Max number of iterations before giving up.",
               "Set to 1 to not give errors, but return empty *ppath* "
               "when a path can not be established.",
               "Top-of-the-atmosphere altitude.",
               "Consider horisontal gradients of refractive index.",
               "Perform double-sided perturbations when calculating "
               "refractive index gradients.")));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_losReverse"),
      DESCRIPTION(R"--(Reverses the directions in *sensor_los*.

The method updates *sensor_los* to have angles of the reversed
direction.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_los"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_pos_losBackwardToAltitude"),
      DESCRIPTION(R"--(Moves *sensor_pos* and *sensor_los* backwards to the target altitude.

The method gives the *sensor_pos* and *sensor_los* at the target altitude
to reach the original *sensor_pos* and *sensor_los* with a geometrical
ppath. That is, the movement is backwards in terms of viewing direction.

If the original *sensor_los* is reversed with respect to the line-of-sight
direction, then set the GIN los_reversed to 1. One such case is that
if *sensor_los* represents surface incidence angles, i.e. holds the
zenith and nadir angle towards the sensor.

There is also *rte_pos_losBackwardToAltitude*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_pos", "sensor_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_pos", "sensor_los", "surface_field"),
      GIN("altitude", "los_is_reversed"),
      GIN_TYPE("Numeric", "Index"),
      GIN_DEFAULT(NODEF, "0"),
      GIN_DESC("Target altitude.",
               "Set to 1 if *rte_los* is valid for the reversed direction.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_pos_losForwardToAltitude"),
      DESCRIPTION(R"--(Moves *sensor_pos* and *sensor_los* forward to the target altitude.

The method gives the *sensor_pos* and *sensor_los* at the target altitude
when forward-propagating the original *sensor_pos* and *sensor_los*
geometrically.

The WSM *IntersectionGeometricAltitude* performs the same operation
but allows to store the new pos and los as other variables. There is
also *rte_pos_losForwardToAltitude*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_pos", "sensor_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_pos", "sensor_los", "surface_field"),
      GIN("altitude"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Target altitude.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseAntenna"),
      DESCRIPTION(R"--(Includes response of the antenna.

The function returns the sensor response matrix after the antenna
characteristics have been included.

The function handles \"multi-beam\" cases where the polarisation
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

\"interp_response\"
Both radiances and the antenna pattern are treated as step-wise
constant functions. The antenna pattern is interpolated to the
*mblock_dlos* directions. At extrapolation, the antenna response
is set to zero. This option considers GIN ``solid_angles``, that
shall be a vector with length matching the rows of *mblock_dlos*.
The values going into *sensor_response* are the interpolated antenna
values times the corresponding solid angle.

\"gridded_dlos\"
This option is more similar to the 1D case. The radiances are treated
as a bi-linear function, but the antenna response is treated as step-
wise constant function (in contrast to 1D). For this option
*mblock_dlos* must match a combination of zenith and azimuth
grids, and this for a particular order. If the zenith and azimuth
grids have 3 and 2 values, respectively, the order shall be::

  [(za1,aa1); (za2,aa1); (za3,aa1); (za1,aa2); (za2,aa2); (za3,aa2)]

Both these grids must be strictly increasing and as for 1D must cover
the antenna response completely.
)--"),
      AUTHORS("Patrick Eriksson", "Mattias Ekstrom"),
      OUT("sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_dlos_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_response",
         "sensor_response_f",
         "sensor_response_pol",
         "sensor_response_dlos",
         "sensor_response_f_grid",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid",
         "antenna_dim",
         "antenna_dlos",
         "antenna_response",
         "sensor_norm"),
      GIN("option_2d", "solid_angles"),
      GIN_TYPE("String", "Vector"),
      GIN_DEFAULT("-", "[]"),
      GIN_DESC("Calculation option for 2D antenna cases. See above for details.",
               "The solid angle of each *mblock_dlos* direction. Only considered "
               "for 2D with \"interp_response\".")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseBackend"),
      DESCRIPTION(R"--(Includes response of the backend (spectrometer).

The function returns the sensor response matrix after the backend
characteristics have been included.

See *f_backend*, *backend_channel_response* and *sensor_norm* for
details on how to specify the backend response.
)--"),
      AUTHORS("Mattias Ekstrom", "Patrick Eriksson"),
      OUT("sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_response",
         "sensor_response_f",
         "sensor_response_pol",
         "sensor_response_dlos",
         "sensor_response_f_grid",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid",
         "f_backend",
         "backend_channel_response",
         "sensor_norm"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseBackendFrequencySwitching"),
      DESCRIPTION(R"--(Frequency switching for a pure SSB reciever.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_response",
         "sensor_response_f",
         "sensor_response_pol",
         "sensor_response_dlos",
         "sensor_response_f_grid",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid",
         "f_backend",
         "backend_channel_response",
         "sensor_norm"),
      GIN("df1", "df2"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Frequency throw for cycle1.", "Frequency throw for cycle2.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseBeamSwitching"),
      DESCRIPTION(R"--(Simulation of \"beam switching\".

The measurement procedure is based on taking the difference between
two spectra measured in different directions, and the calculation
set-up must treat exactly two observation directions.

The returned spectrum is y = w1*y + w2*y2, where y1 and w1 are the
spectrum and weight for the first direction, respectively (y2 and
(w2 defined correspondingly for the second direction).

Zenith and azimuth angles after beam switching are set to the
values of the second direction.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_dlos_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_response",
         "sensor_response_f",
         "sensor_response_pol",
         "sensor_response_dlos",
         "sensor_response_f_grid",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid"),
      GIN("w1", "w2"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT("-1", "1"),
      GIN_DESC("Weight for values from first viewing direction.",
               "Weight for values from second viewing direction.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseFillFgrid"),
      DESCRIPTION(R"--(Polynomial frequency interpolation of spectra.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_response",
         "sensor_response_f",
         "sensor_response_pol",
         "sensor_response_dlos",
         "sensor_response_f_grid",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid"),
      GIN("polyorder", "nfill"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT("3", "2"),
      GIN_DESC("Polynomial order of interpolation",
               "Number of points to insert in each gap of f_grid")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseFrequencySwitching"),
      DESCRIPTION(R"--(Simulation of \"frequency switching\".

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_response",
         "sensor_response_f",
         "sensor_response_pol",
         "sensor_response_dlos",
         "sensor_response_f_grid",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseIF2RF"),
      DESCRIPTION(R"--(Converts sensor response variables from IF to RF.

The function converts intermediate frequencies (IF) in
*sensor_response_f* and *sensor_response_f_grid* to radio
frequencies (RF). This conversion is needed if the frequency
translation of a mixer is included and the position of backend
channels are specified in RF.

A direct frequency conversion is performed. Values are not
sorted in any way.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_response_f", "sensor_response_f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_response_f", "sensor_response_f_grid", "lo", "sideband_mode"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseInit"),
      DESCRIPTION(R"--(Initialises the variables summarising the sensor response.

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
)--"),
      AUTHORS("Mattias Ekstrom", "Patrick Eriksson"),
      OUT("sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid",
          "sensor_response_pol_grid",
          "sensor_response_dlos_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         "mblock_dlos",
         "antenna_dim",
         
         "sensor_norm"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseMetMM"),
      DESCRIPTION(R"--(Sensor setup for meteorological millimeter instruments.

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
)--"),
      AUTHORS("Oliver Lemke", "Patrick Eriksson"),
      OUT("antenna_dim",
          "mblock_dlos",
          "sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid",
          "sensor_response_pol_grid",
          "sensor_response_dlos_grid",
          "sensor_norm"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         
         "f_grid",
         "f_backend",
         "channel2fgrid_indexes",
         "channel2fgrid_weights",
         "iy_unit",
         "antenna_dlos",
         "met_mm_polarisation",
         "met_mm_antenna"),
      GIN("use_antenna", "mirror_dza"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT("0", "0"),
      GIN_DESC("Flag to enable (1) or disable (0) antenna.",
               "Flag to include second part of swath (only 3D, see above).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseMixer"),
      DESCRIPTION(R"--(Includes response of the mixer of a heterodyne system.

The function returns the sensor response matrix after the mixer
characteristics have been included. Frequency variables are
converted from radio frequency (RF) to intermediate frequency (IF).
The returned frequency grid covers the range [0,max_if], where
max_if is the highest IF covered by the sideband response grid.

See *lo* and *sideband_response* for details on how to specify the
mixer response
)--"),
      AUTHORS("Mattias Ekstrom", "Patrick Eriksson"),
      OUT("sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_response",
         "sensor_response_f",
         "sensor_response_pol",
         "sensor_response_dlos",
         "sensor_response_f_grid",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid",
         "lo",
         "sideband_response",
         "sensor_norm"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseMixerBackendPrecalcWeights"),
      DESCRIPTION(R"--(Includes pre-calculated response covering mixer and backend.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_response",
         "sensor_response_f",
         "sensor_response_pol",
         "sensor_response_dlos",
         "sensor_response_f_grid",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid",
         "f_backend",
         "channel2fgrid_indexes",
         "channel2fgrid_weights"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseMultiMixerBackend"),
      DESCRIPTION(R"--(Handles mixer and backend parts for an instrument having multiple
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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_response",
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
         "sensor_norm"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responsePolarisation"),
      DESCRIPTION(R"--(Extraction of non-default polarisation components.

The default is to output the Stokes elements I, Q, U and V (up to
``stokes_dim``). This method allows to change the \"polarisation\" of
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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_pol_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_response",
         "sensor_response_f",
         "sensor_response_pol",
         "sensor_response_dlos",
         "sensor_response_f_grid",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid",
         
         "iy_unit",
         "instrument_pol"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseStokesRotation"),
      DESCRIPTION(R"--(Includes a rotation of the Stokes H and V directions.

The method applies the rotations implied by *stokes_rotation*.
See the description of that WSV for details.

This method does not change the size of *sensor_response*, and
the auxiliary variables (sensor_response_f etc.) are not changed.

To apply the method, ``stokes_dim`` must be >= 3. The complete effect
of the rotation can not be determibed with lower ``stokes_dim``.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_response"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_response",
         "sensor_response_f_grid",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid",
         
         "stokes_rotation"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseSimpleAMSU"),
      DESCRIPTION(R"--(Simplified sensor setup for an AMSU-type instrument.

This method allows quick and simple definition of AMSU-type
sensors. Assumptions:

1. Pencil beam antenna.
2. Double sideband receivers.
3. Sideband mode \"upper\"
4. The channel response is rectangular.

Under these assumptions the only inputs needed are the LO positions,
the offsets from the LO, and the IF bandwidths. They are provieded
in sensor_description_amsu.
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("f_grid",
          "antenna_dim",
          "mblock_dlos",
          "sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid",
          "sensor_response_pol_grid",
          "sensor_response_dlos_grid",
          "sensor_norm"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(  "sensor_description_amsu"),
      GIN("spacing"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(".1e9"),
      GIN_DESC("Desired grid spacing in Hz.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseGenericAMSU"),
      DESCRIPTION(R"--(Simplified sensor setup for an AMSU-type instrument.

This function is derived from 'sensor_responseSimpleAMSU' 
but is more generalized since the number of passbands in each 
can be in the range from 1 to 4 - in order to correctly simulate
AMSU-A type sensors 

This method allows quick and simple definition of AMSU-type
sensors. Assumptions:

1. Pencil beam antenna.
2. 1-4 Passband/sidebands per channel.
3. Sideband mode \"upper\"
4. The channel response is rectangular.

Under these assumptions the only inputs needed are the LO positions,
the offsets from the LO, and the IF bandwidths. They are provided
in sensor_description_amsu.
)--"),
      AUTHORS("Oscar Isoz"),
      OUT("f_grid",
          "antenna_dim",
          "mblock_dlos",
          "sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid",
          "sensor_response_pol_grid",
          "sensor_response_dlos_grid",
          "sensor_norm"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(  "sensor_description_amsu"),
      GIN("spacing"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(".1e9"),
      GIN_DESC("Desired grid spacing in Hz.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sensor_responseWMRF"),
      DESCRIPTION(R"--(Adds WMRF weights to sensor response.

This method adds a spectrometer response that has been calculated
with the weighted mean of representative frequencies (WMRF) method. It
consists of a set of selected frequencies, and associated weights.
)--"),
      AUTHORS(
          "Stefan Buehler, based on Patrick Erikssons sensor_responseBackend"),
      OUT("sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_response",
         "sensor_response_f",
         "sensor_response_pol",
         "sensor_response_dlos",
         "sensor_response_f_grid",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid",
         "wmrf_weights",
         "f_backend"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord(NAME("SetNumberOfThreads"),
               DESCRIPTION(R"--(Change the number of threads used by ARTS.
)--"),
               AUTHORS("Oliver Lemke"),
               OUT(),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN("nthreads"),
               GIN_TYPE("Index"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Number of threads.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Sleep"),
      DESCRIPTION(R"--(Sleeps for a number of seconds
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("time"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Time to sleep for in seconds")));

  md_data_raw.push_back(create_mdrecord(
      NAME("timeSleep"),
      DESCRIPTION(R"--(Sleeps until time has been reached.
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("time"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("timeSet"),
      DESCRIPTION(R"(Sets the time.

The time format is similar to C:s strftime format
    %Y-%m-%d %H:%M:%S

The one exception is that our format accepts, but does not need to use,
decimals of the second in addtion.  Note that the native time resolution
for the decimal is the most that can be kept.  For some systems, this is
nano-seconds, and for others that is micro-seconds.  Please see with your
vendor.

A default argument is a close approximation to the formal first commit to
the ARTS codebase.  It is there to give an example of how the format looks.
)"),
      AUTHORS("Richard Larsson"),
      OUT("time"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("time_str"),
      GIN_TYPE("String"),
      GIN_DEFAULT("2000-03-11 14:39:37.0"),
      GIN_DESC("A time stamp string in the default format")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sparse_f_gridFromFrequencyGrid"),
      DESCRIPTION(R"--(Outputs the sparse frequency grid in *propmat_clearskyAddLines*
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("sparse_f_grid"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("A sparse frequency grid."),
      IN("f_grid"),
      GIN("sparse_df", "speedup_option"),
      GIN_TYPE("Numeric", "String"),
      GIN_DEFAULT("0", "None"),
      GIN_DESC("The grid sparse separation", "Speedup logic")));

  md_data_raw.push_back(create_mdrecord(
      NAME("SparseSparseMultiply"),
      DESCRIPTION(R"--(Multiplies a Sparse with another Sparse, result stored in Sparse.

Makes the calculation: M = M1 * M2
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("M"),
      GOUT_TYPE("Sparse"),
      GOUT_DESC("Product, can be same variable as any of the inputs."),
      IN(),
      GIN("M1", "M2"),
      GIN_TYPE("Sparse", "Sparse"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Left sparse matrix (dimension m x n).",
               "Right sparse matrix (dimension n x p).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("SparseIdentity"),
      DESCRIPTION(R"--(Returns a sparse dentity matrix.

The size of the matrix created is n x n. Default is to return a
true identity matrix (I), but you can also select another value
along the diagonal be setting ``value``. That is, the output is
value*I.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Sparse"),
      GOUT_DESC("Sparse output matrix"),
      IN(),
      GIN("n", "value"),
      GIN_TYPE("Index", "Numeric"),
      GIN_DEFAULT(NODEF, "1"),
      GIN_DESC("Size of the matrix", "The value along the diagonal.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("spectral_irradiance_fieldFromSpectralRadianceField"),
      DESCRIPTION(R"--(Calculates the spectral irradiance from *spectral_radiance_field*.

The *spectral_radiance_field* is integrated over the angular grids
according to the grids set by *AngularGridsSetFluxCalc*.
See *AngularGridsSetFluxCalc* to set *za_grid*, *aa_grid*, and 
*za_grid_weights*.
)--"),
      AUTHORS("Manfred Brath"),
      OUT("spectral_irradiance_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("spectral_radiance_field", "za_grid", "aa_grid", "za_grid_weights"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("spectral_radiance_fieldClearskyPlaneParallel"),
      DESCRIPTION(R"--(Clear-sky radiance field of a plane parallel atmosphere.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("spectral_radiance_field"),
      GOUT("trans_field"),
      GOUT_TYPE("Tensor3"),
      GOUT_DESC("Dimensions: [f_grid,p_grid,za_grid]. See further above."),
      IN("propmat_clearsky_agenda",
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
         "za_grid"),
      GIN("use_parallel_za"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Flag to select parallelization over zenith angles.")));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("spectral_radiance_fieldCopyCloudboxField"),
      DESCRIPTION(R"--(Set *spectral_radiance_field* to be a copy of *cloudbox_field*.

This method can only be used for 1D atmospheres and if the cloud
box covers the complete atmosphere. For such case, the two fields
cover the same atmospheric volume and a direct copying can be made.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("spectral_radiance_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         "p_grid",
         "cloudbox_on",
         "cloudbox_limits",
         "cloudbox_field"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("spectral_radiance_fieldExpandCloudboxField"),
      DESCRIPTION(R"--(Uses and expands *cloudbox_field* to set *spectral_radiance_field*.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("spectral_radiance_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky_agenda",
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
         "za_grid"),
      GIN("use_parallel_za"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("Flag to select parallelization over zenith angles.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("specular_losCalc"),
      DESCRIPTION(R"--(Calculates the specular direction of surface reflections.

The default is to consider surface topography when calculating the
specular direction. That is, the variation of ``surface_elevation``
is allowed to affect the angles of *specular_los*. This impact can
be deactivated by setting ``ignore_topography`` to 1. In this case,
the zenith angle of the specular direction is simply 180-rtp_los[0]
and the azimuth angle is the same as the one in *rtp_los*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("specular_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("surface_field", "rtp_pos", "rtp_los"),
      GIN("ignore_topography"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("Flag to control if surface slope is considered or not.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("specular_losCalcOldNoTopography"),
      DESCRIPTION(R"--(Calculates the specular direction of surface reflections for horisontal
surfaces.

In contrast to ``specular_losCalcOld``, this method ignores the topography
implied by ``z_surface``. That is, any slope of the surface is ignored.

The typical application of this WSM should be water surfaces (lakes and
oceans).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("specular_los", "surface_normal"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("rtp_los"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("sunsAddSingleBlackbody"),
      DESCRIPTION(R"--(Adds a single blackbody to *suns*

Important note:
For a Sol-like sun there are huge differences in the UV-range 
between the actual sun spectrum and the blackbody spectrum
with the effective temperature of the sun. The blackbody sun\"
strongly overestimates the UV radiation.
)--"),
      AUTHORS("Jon Petersen"),
      OUT("suns",
          "suns_do"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("suns",
         "f_grid"),
      GIN("radius",
          "distance",
          "temperature",
          "latitude",
          "longitude"),
      GIN_TYPE("Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric"),
      GIN_DEFAULT("6.963242e8",
                  "1.495978707e11",
                  "5772",
                  "0",
                  "0"),
      GIN_DESC("The radius of the sun in meter. "
               "Default is the radius of our sun. ",
               "The average distance between the sun and the planet in meter. "
               "Default value is set to 1 a.u. ",
               "The effective temperature of the suns photosphere in Kelvin. "
               "Default is the temperature of our sun - 5772 Kelvin ",
               "The latitude or the zenith position of the sun in the sky. ",
               "The longitude or azimuthal position of the sun in the sky. ")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sunsAddSingleFromGrid"),
      DESCRIPTION(R"--(Extracts a sun spectrum from a field of such data and
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
)--"),
      AUTHORS("Jon Petersen"),
      OUT("suns",
          "suns_do"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("suns",
         "f_grid"),
      GIN("sun_spectrum_raw",
          "radius",
          "distance",
          "temperature",
          "latitude",
          "longitude",
          "description"),
      GIN_TYPE("GriddedField2",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "String"),
      GIN_DEFAULT(NODEF,
                  "6.963242e8",
                  "1.495978707e11",
                  "-1",
                  "0",
                  "0",
                  "Sun spectrum from Griddedfield."),
      GIN_DESC("Raw data for monochromatic irradiance spectra. ",
               "The radius of the sun in meter. "
               "Default is the radius of our sun. ",
               "The average distance between the center of the sun and the  "
               "center of the planet in meter. "
               "Default value is set to 1 a.u. ",
               "The temperature of the padding if the f_grid is outside the  "
               "sun spectrum data. Choose 0 for 0 at the edges or a effective "
               "temperature for a padding using plack's law. ",
               "The latitude or the zenith position of the sun in the sky. ",
               "The longitude or azimuthal position of the sun in the sky. ",
               "The description of the sun. ")));

  md_data_raw.push_back(create_mdrecord(
      NAME("sunsAddSingleFromGridAtLocation"),
      DESCRIPTION(R"--(Extracts a sun spectrum measured at the given location
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
)--"),
      AUTHORS("Jon Petersen"),
      OUT("suns",
          "suns_do"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("suns",
          "f_grid",
          
          "surface_field"),
      GIN("sun_spectrum_raw",
          "radius",
          "distance",
          "temperature",
          "zenith",
          "azimuth",
          "description",
          "location_latitude",
          "location_longitude",
          "location_altitude"),
      GIN_TYPE("GriddedField2",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "Numeric",
               "String",
               "Numeric",
               "Numeric",
               "Numeric"),
      GIN_DEFAULT(NODEF,
                  "6.963242e8",
                  "1.495978707e11",
                  "-1",
                  "0",
                  "0",
                  "Sun spectrum from Griddedfield.",
                  "0",
                  "0",
                  "1e5"),
      GIN_DESC("Raw data for monochromatic irradiance spectra. ",
               "The radius of the sun in meter. "
               "Default is the radius of our Sun. ",
               "The distance between the location and the  "
               "center of the sun in meter. "
               "Default value is set to 1 a.u. ",
               "The temperature of the padding if the f_grid is outside the  "
               "sun spectrum data. Choose 0 for 0 at the edges or a effective "
               "temperature for a padding using plack's law. ",
               "Zenith angle of the sun in the sky. ",
               "Azimuthal angle of the sun in the sky. ",
               "The description of the sun. ",
               "The latitude of the sun spectrum measurement. ",
               "The longitude of the sun spectrum measurement. ",
               "The altitude of the sun spectrum measurement. ")));  
      
  md_data_raw.push_back(create_mdrecord(
      NAME("sunsOff"),
      DESCRIPTION(R"--(Turns all calculations with suns off
)--"),
      AUTHORS("Jon Petersen"),
      OUT("suns_do",
          "suns"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("StringJoin"),
      DESCRIPTION(R"--(Concatenate two or more strings.

The output string is overwritten, but is allowed to appear
in the input list. Up to 10 strings can be concatenated at once.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("String"),
      GOUT_DESC("Concatenated string."),
      IN(),
      GIN("in1",
          "in2",
          "in3",
          "in4",
          "in5",
          "in6",
          "in7",
          "in8",
          "in9",
          "in10"),
      GIN_TYPE("String",
               "String",
               "String",
               "String",
               "String",
               "String",
               "String",
               "String",
               "String",
               "String"),
      GIN_DEFAULT(NODEF, NODEF, "", "", "", "", "", "", "", ""),
      GIN_DESC("Input text string.",
               "Input text string.",
               "Input text string.",
               "Input text string.",
               "Input text string.",
               "Input text string.",
               "Input text string.",
               "Input text string.",
               "Input text string.",
               "Input text string.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("surface_normalCalc"),
      DESCRIPTION(R"--(Calculates the surface's local normal.

The default is to consider surface topography when calculating the
normal direction. That is, the variation of ``surface_elevation``
is allowed to affect the angles of *surface_normal*. This impact can
be deactivated by setting ``ignore_topography`` to 1. In this case,
the zenith angle of *surface_normal* becomes 0.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_normal"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("surface_field", "rtp_pos"),
      GIN("ignore_topography"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("Flag to control if surface slope is considered or not.")));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("z_surfaceFromFileAndGrid"),
      DESCRIPTION(R"--(Sets the surface altitude for a given latitude and longitude grid.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("surface_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("lat_grid", "lon_grid"),
      GIN("filename", "interp_order", "set_lowest_altitude_to_zero"),
      GIN_TYPE("String", "Index", "Index"),
      GIN_DEFAULT(NODEF, "1", "0"),
      GIN_DESC(
          "File of GriddedField2 with surface altitudes gridded",
          "Interpolation order",
          "Index that sets the lowest altitude to 0 to ignore sub-surface pressures/altitudes")));
*/

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("z_surfaceConstantAltitude"),
      DESCRIPTION(R"--(Sets the surface altitude to a constant. Defaults to zero.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("surface_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("lat_grid", "lon_grid"),
      GIN("altitude"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("0"),
      GIN_DESC("The constant altitude.")));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("surface_fieldSet"),
      DESCRIPTION(R"--(Make the surface field hold value at the key.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("surface_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("surface_field"),
      GIN("value", "key"),
      GIN_TYPE("Numeric,GriddedField2", "String"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Value to set", "Key to set value at")));

  md_data_raw.push_back(create_mdrecord(
      NAME("surface_fieldSetProp"),
      DESCRIPTION(R"--(Make the surface field hold value at the key.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("surface_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("surface_field"),
      GIN("value", "key"),
      GIN_TYPE("Numeric,GriddedField2", "String"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Value to set", "Key to set value at")));

  md_data_raw.push_back(create_mdrecord(
      NAME("surface_fieldSetType"),
      DESCRIPTION(R"--(Make the surface field hold value at the key.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("surface_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("surface_field"),
      GIN("value", "key"),
      GIN_TYPE("Numeric,GriddedField2", "String"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Value to set", "Key to set value at")));

  md_data_raw.push_back(create_mdrecord(
      NAME("surfaceBlackbody"),
      DESCRIPTION(R"--(Creates variables to mimic a blackbody surface.

This method sets up *surface_los*, *surface_rmatrix* and
*surface_emission* for *surface_rtprop_agenda*. Here, *surface_los*
and *surface_rmatrix* are set to be empty, and *surface_emission*
to hold blackbody radiation for a temperature of *surface_skin_t*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         "f_grid",
         
         "rtp_pos",
         "rtp_los",
         "surface_point"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("surfaceFastem"),
      DESCRIPTION(R"--(Usage of FASTEM together with MC and DOIT.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         
         "f_grid",
         "rtp_pos",
         "rtp_los",
         "surface_skin_t"),
      GIN("salinity",
          "wind_speed",
          "wind_direction",
          "transmittance",
          "fastem_version"),
      GIN_TYPE("Numeric", "Numeric", "Numeric", "Vector", "Index"),
      GIN_DEFAULT("0.035", NODEF, "0", NODEF, "6"),
      GIN_DESC("Salinity, 0-1. That is, 3% is given as 0.03.",
               "Wind speed.",
               "Wind direction. See futher above.",
               "Transmittance along path of downwelling radiation. A vector "
               "with the same length as *f_grid*.",
               "The version of FASTEM to use.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("surfaceFlatRefractiveIndex"),
      DESCRIPTION(R"--(Creates variables to mimic specular reflection by a (flat) surface
where the complex refractive index is specified.

The dielectric properties of the surface are described by
*surface_complex_refr_index*. The Fresnel equations are used to
calculate amplitude reflection coefficients. The method can thus
result in that the reflection properties differ between frequencies
and polarisations.

Local thermodynamic equilibrium is assumed, which corresponds to
that the reflection and emission coefficients add up to 1.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         
         "rtp_pos",
         "rtp_los",
         "specular_los",
         "surface_skin_t",
         "surface_complex_refr_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("surfaceFlatReflectivity"),
      DESCRIPTION(R"--(Creates variables to mimic specular reflection by a (flat) surface
where *surface_reflectivity* is specified.

Works basically as *surfaceFlatScalarReflectivity* but is more
general as vector radiative transfer is more properly handled. See
the ARTS theory document (ATD) for details around how
*surface_emission* is determined. In the nomenclature of ATD,
*surface_reflectivity* gives R.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         
         "rtp_pos",
         "rtp_los",
         "specular_los",
         "surface_skin_t",
         "surface_reflectivity"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("surfaceFlatRvRh"),
      DESCRIPTION(R"--(Creates variables to mimic specular reflection by a (flat) surface
where *surface_rv_rh* is specified.

This method assumes that the reflection at vertical and horizontal
polarisation differs. As power reflection coefficients are provided
there is no information at hand on phase shifts between polarisations,
and they are simply assumed to be zero. These assumptions result in
that *surface_emission* is set to zero for positions corresponding to
U and V, and that all diagonal elementsof  *surface_rmatrix* are equal
(the mean of rv and rh). Further, all off-diagonal elements of
*surface_rmatrix* are all zero except for (0,1) and (1,0).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         
         "rtp_pos",
         "rtp_los",
         "specular_los",
         "surface_skin_t",
         "surface_rv_rh"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("surfaceFlatScalarReflectivity"),
      DESCRIPTION(R"--(Creates variables to mimic specular reflection by a (flat) surface
where *surface_scalar_reflectivity* is specified.

This method assumes that the reflection at vertical and horizontal
polarisation is identical. This assumption includes that there is no
phase shift between polarisations. These assumptions result in that
*surface_emission* is set to zero for positions corresponding to Q,
U and V, and that *surface_rmatrix* becomes a diagonal matrix (with
all elements on the diagonal equal to the specified reflectivity).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         
         "rtp_pos",
         "rtp_los",
         "specular_los",
         "surface_skin_t",
         "surface_scalar_reflectivity"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("surfaceLambertianSimple"),
      DESCRIPTION(R"--(Creates variables to mimic a Lambertian surface.

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
that the reflection and emission coefficients \"add up to 1\".

For 2D and 3D, the down-welling directions are placed along the
the viewing direction, e.g. for 3D the azimuth angle is kept constant.
In 2D and 3D surface topography can exist, and to avoid getting views
going directly into the surface, angels are not distributed over 90 deg,
but 90-abs(surface_normal[0]).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         
         "rtp_pos",
         "rtp_los",
         "surface_normal",
         "surface_skin_t",
         "surface_scalar_reflectivity"),
      GIN("lambertian_nza", "za_pos"),
      GIN_TYPE("Index", "Numeric"),
      GIN_DEFAULT("9", "0.5"),
      GIN_DESC("Number of downwelling streams.",
               "Position of angle in *surface_los* inside ranges of zenith "
               "angle grid. See above.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("surfaceMapToLinearPolarisation"),
      DESCRIPTION(R"--(Convert surface RT properties to a linear polarisation.

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
theory guide, in section \"Rotated modified Stokes vector\".

In general it should suffice to set local_stokes_dim to 2, that gives
slightly faster calculations. A local_stokes_dim of 3 handles any case
correctly.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_emission", "surface_rmatrix"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("surface_emission", "surface_rmatrix"),
      GIN("pol_angle"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Polarisation angle, see above.")));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("surfaceTelsem"),
      DESCRIPTION(R"--(Compute surface emissivities using the TELSEM 2 model.

This method uses second version of the TELSEM model for calculating
land surface emissivities (F. Aires et al, \"A Tool to Estimate 
Land‐Surface Emissivities at Microwave frequencies (TELSEM) for use
in numerical weather prediction\" Quarterly Journal of the Royal
Meteorological Society, vol. 137, (656), pp. 690-699, 2011.)
This methods computes land surface emissivities for a given pencil beam
using a given TELSEM2 atlas.

The input must satisfy the following conditions, otherwise an error is thrown:
 - The input frequencies (*f_grid*) must be within the range [5 GHz, 900 GHz]
 - The skin temperature (*surface_skin_t*) must be within the range
   [180 K, 360 K]

A TELSEM atlas contains only suface emissivities for locations that are
classified as land. By default this WSM will throw an error if the
pencil beam hits the surface at a position that is not contained in the
given atlas.

The above behavior can be avoided by setting ``d_max`` to a positive value.
This enables nearest neighbor interpolation, which assigns the emissivities
of the nearest found cell in the atlas to the given position. In this case,
an error is only thrown if the distance of the found neighbor is higher
than the provided value of ``d_max``.

You can limit the final reflectivity applied by setting ``r_min`` and ``r_max``.

To extract a land-sea mask from a given telsem atlas see the WSM
*telsemSurfaceTypeLandSea*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         
         "f_grid",
         "lat_grid",
         "lat_true",
         "lon_true",
         "rtp_pos",
         "rtp_los",
         "specular_los",
         "surface_skin_t"),
      GIN("atlas", "r_min", "r_max", "d_max"),
      GIN_TYPE("TelsemAtlas", "Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "0", "1", "-1.0"),
      GIN_DESC("The Telsem atlas to use for the emissivity calculation.",
               "Minimum allowed value for reflectivity to apply.",
               "Maximum allowed value for reflectivity to apply.",
               "Maximum allowed distance in meters for nearest neighbor"
               " interpolation in meters. Set to a negative value or zero "
               " to disable interpolation.")));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("surfaceTessem"),
      DESCRIPTION(R"--(TESSEM sea surface microwave emissivity parametrization.

This method computes surface emissivity and reflectivity matrices for
ocean surfaces using the TESSEM emissivity model: Prigent, C., et al.
Sea‐surface emissivity parametrization from microwaves to millimetre
waves, QJRMS, 2017, 143.702: 596-605.

The validity range of the parametrization of is 10 to 700 GHz, but for
some extra flexibility frequencies between 5 and 900 GHz are accepted.
The accepted temperaute range for *surface_skin_t* is [260.0 K, 373.0 K]

The model itself is represented by the neural networks in
*tessem_neth* and *tessem_netv*.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         
         "f_grid",
         "rtp_pos",
         "rtp_los",
         "surface_skin_t",
         "tessem_neth",
         "tessem_netv"),
      GIN("salinity", "wind_speed"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT("0.035", NODEF),
      GIN_DESC("Salinity, 0-1. That is, 3% is given as 0.03.", "Wind speed.")));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("surface_complex_refr_indexFromGriddedField5"),
      DESCRIPTION(R"--(Extracts complex refractive index from a field of such data.

The method allows to obtain *surface_complex_refr_index* by
interpolation of a geographical field of such data. The position
for which refraction shall be extracted is given by *rtp_pos*.
The refractive index field is expected to be stored as:

- GriddedField5:

  - Vector f_grid[N_f]
  - Vector T_grid[N_T]
  - ArrayOfString Complex[2]
  - Vector \"Latitude\"  [N_lat]
  - Vector \"Longitude\" [N_lon]
  - Tensor5 data[N_f][N_T][2][N_lat][N_lon]

Definition and treatment of the three first dimensions follows
*complex_refr_index*, e.g. the temperature grid is allowed
to have length 1. The grids for latitude and longitude must have
a length of >= 2 (ie. no automatic expansion).

Hence, this method performs an interpolation only in the lat and
lon dimensions, to a single point. The remaining GriddedField3 is
simply returned as *surface_complex_refr_index*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_complex_refr_index"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN( "lat_grid", "lat_true", "lon_true", "rtp_pos"),
      GIN("complex_refr_index_field"),
      GIN_TYPE("GriddedField5"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("A field of complex refractive index.")));
*/
/*
  md_data_raw.push_back(create_mdrecord(
      NAME("surface_reflectivityFromGriddedField6"),
      DESCRIPTION(R"--(Extracts surface reflectivities from a field of such data.

This method allows to specify a field of surface reflectivity for
automatic interpolation to points of interest. The position and
direction for which the reflectivity shall be extracted are given
by *rtp_pos* and *rtp_los*. The reflectivity field is expected to
be stored as:

- GriddedField6:

  - Vector \"Frequency\"       [N_f]
  - Vector \"Stokes element\"  [N_s1]
  - Vector \"Stokes_element\"  [N_s2]
  - Vector \"Incidence angle\" [N_ia]
  - Vector \"Latitude\"        [N_lat]
  - Vector \"Longitude\"       [N_lon]
  - Tensor6 data[N_f][N_s1][N_s2][N_ia][N_lat][N_lon]

Grids for incidence angle, latitude and longitude must have a
length of >= 2 (ie. no automatic expansion). If the frequency grid
has length 1, this is taken as that the reflectivity is constant,
following the definition of *surface_scalar_reflectivity*.
The data can cover higher Stokes dimensionalities than set by
``stokes_dim``. Data for non-used Stokes elements are just cropped.
The order between the two Stokes dimensions is the same as in
*surface_reflectivity* and surface_rmatrix*.

The interpolation is done in steps:
   (1) Linear interpolation for lat and lon (std. extrapolation).
   (2) Interpolation in incidence angle (std. extrapolation).
       If the grid has a length of >= 4, cubic interpolation is
       applied. Otherwise linear interpolation.
   (3) Linear interpolation in frequency (if input data have more
       than one frequency).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_reflectivity"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         "f_grid",
         "lat_grid",
         "lat_true",
         "lon_true",
         "rtp_pos",
         "rtp_los"),
      GIN("r_field"),
      GIN_TYPE("GriddedField6"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("A field of surface reflectivities")));
*/
/*
  md_data_raw.push_back(create_mdrecord(
      NAME("surface_rtpropFromTypesAverage"),
      DESCRIPTION(R"--(Extracts surface RT properties by averaging.

This method allows to let one pencil beam calculation represent
an area when it comes to surface RT properties. The surface is
sampled at a set of positions. The sampling is defined as a set
of angles by the WSV *dlos*.

A weight must be specified for each angle by *dlos_weight_vector*.
These weights should represent the solid angle each *dlos* covers,
but can also include other weighting factors such as antenna pattern.
The sum of *dlos_weight_vector* is internally normalised to 1.

All the output variables are a weighted average between the surface
types inside the area sampled.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_type_mix",
          "surface_skin_t",
          "surface_los",
          "surface_rmatrix",
          "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         
         "lat_grid",
         "lat_true",
         "lon_true",
         "rtp_pos",
         "rtp_los",
         "surface_type_mask",
         "surface_rtprop_agenda_array",
         "z_sensor",
         "dlos",
         "dlos_weight_vector"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
*/
  md_data_raw.push_back(create_mdrecord(
      NAME("surface_rtpropFromTypesManual"),
      DESCRIPTION(R"--(Extracts surface RT properties by manual selection of surface type.

The surface type to apply is selected by the GIN argument.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_type_mix",
          "surface_skin_t",
          "surface_los",
          "surface_rmatrix",
          "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         "rtp_pos",
         "rtp_los",
         "surface_rtprop_agenda_array"),
      GIN("surface_type"),
      GIN_TYPE("Index"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Selected surface type")));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("surface_rtpropFromTypesNearest"),
      DESCRIPTION(R"--(Extracts surface RT properties from nearest surface type.

The surface type is set by nearest interpolation of *surface_type_mask*
and the corresponding agenda in *surface_rtprop_agenda_array* is
called to obtain the local radiative properties of the surface.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_type_mix",
          "surface_skin_t",
          "surface_los",
          "surface_rmatrix",
          "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         "lat_grid",
         "lat_true",
         "lon_true",
         "rtp_pos",
         "rtp_los",
         "surface_type_mask",
         "surface_rtprop_agenda_array"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("surface_rtpropInterpFreq"),
      DESCRIPTION(R"--(Interpolates surface RT properties in frequency.

The WSVs *surface_rmatrix* and *surface_emission* are inter-
polated linearly in frequency. The original frequency is given
by *f_grid*, and there is an interpolation to new frequency grid.
The function resets *f_grid* to the new grid.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("f_grid",
          "surface_rmatrix",
          "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         "surface_rmatrix",
         "surface_emission"),
      GIN("f_new"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("New frequency grid")));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("surface_scalar_reflectivityFromGriddedField4"),
      DESCRIPTION(R"--(Extracts scalar surface reflectivities from a field of such data.

This method allows to specify a field of surface reflectivity for
automatic interpolation to points of interest. The position and
direction for which the reflectivity shall be extracted are given
by *rtp_pos* and *rtp_los*. The reflectivity field is expected to
be stored as:

- GriddedField4:

  - Vector \"Frequency\"       [N_f]
  - Vector \"Incidence angle\" [N_ia]
  - Vector \"Latitude\"        [N_lat]
  - Vector \"Longitude\"       [N_lon]
  - Tensor4 data[N_f][N_ia][N_lat][N_lon]

Grids for incidence angle, latitude and longitude must have a
length of >= 2 (ie. no automatic expansion). If the frequency grid
has length 1, this is taken as the reflectivity is constant,
following the definition of *surface_scalar_reflectivity*.

The interpolation is done in steps:

1. Linear interpolation for lat and lon (std. extrapolation).
2. Interpolation in incidence angle (std. extrapolation).
   If the grid has a length of >= 4, cubic interpolation is
   applied. Otherwise linear interpolation.
3. Linear interpolation if frequency (if input data have more
   than one frequency).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_scalar_reflectivity"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         "f_grid",
         "lat_grid",
         "lat_true",
         "lon_true",
         "rtp_pos",
         "rtp_los"),
      GIN("r_field"),
      GIN_TYPE("GriddedField4"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("A field of scalar surface reflectivities")));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("surface_scalar_reflectivityFromSurface_rmatrix"),
      DESCRIPTION(R"--(Sets *surface_scalar_reflectivity* based on *surface_rmatrix*.

For each frequency f, *surface_scalar_reflectivity* is set to
the sum of surface_rmatrix(joker,f,0,0).
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_scalar_reflectivity"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("surface_rmatrix"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("SurfaceBlackbody"),
      DESCRIPTION(R"--(Blackbody surface, with support for Jacobian calculations.

See *surfaceBlackbody* and *SurfaceFastem* for complementary
information.

For this method, ``surface_props_data`` must contain these data:
  \"Skin temperature\"

*dsurface_emission_dx* is calculated analytically.
*surface_rmatrix* and *dsurface_rmatrix_dx* are set to 0.
)--"),
      AUTHORS("Marc Prange"),
      OUT("surface_los",
          "surface_rmatrix",
          "dsurface_rmatrix_dx",
          "surface_emission",
          "dsurface_emission_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("dsurface_rmatrix_dx",
         "dsurface_emission_dx",
         
         "lat_grid",
         "lon_grid",
         "f_grid",
         "rtp_pos",
         "rtp_los",
         "surface_field",
         "surface_props_names",
         "dsurface_names",
         "jacobian_do"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("SurfaceDummy"),
      DESCRIPTION(R"--(Dummy method for *iy_surface_agenda*.

If you don't make use of ``surface_props_data`` and associated
variables, include this method *iy_surface_agenda*. The method
just checks that the variables of concern are set to be empty,
and you don't need to include calls of *Ignore* and *Touch* in
the agenda.

If you use a method of SurfaceSomething type, you don't need
this one.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("dsurface_rmatrix_dx", "dsurface_emission_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("dsurface_rmatrix_dx",
         "dsurface_emission_dx",
         "surface_props_names",
         "dsurface_names",
         "jacobian_do"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("SurfaceFastem"),
      DESCRIPTION(R"--(FASTEM sea surface microwave emissivity parametrization.

This method allows to use FASTEM in retrievals requiring the
Jacobian. Otherwise as *surfaceFastem*. See *SurfaceBlackbody
for general remarks about methods of SurfaceSomething type.

This is an example on a surface method starting with a capital S,
e.g. SurfaceSomething. These methods differ from the methods
named as surfaceSomething in two ways:
  1. The surface properties to apply are taken directly from
     ``surface_props_data``.
  2. The Jacobian with respect to the surface properties can be
     obtained.
The Jacobian can be obtained for all variables in ``surface_props_data``
that the method of concern is using.

For this method, ``surface_props_data`` must contain these data:
  \"Water skin temperature\"
  \"Wind speed\"
  \"Wind direction\"
  \"Salinity\
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los",
          "surface_rmatrix",
          "dsurface_rmatrix_dx",
          "surface_emission",
          "dsurface_emission_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("dsurface_rmatrix_dx",
         "dsurface_emission_dx",
         
         "lat_grid",
         "lon_grid",
         "f_grid",
         "rtp_pos",
         "rtp_los",
         "surface_field",
         "surface_props_names",
         "dsurface_names",
         "jacobian_do"),
      GIN("transmittance", "fastem_version"),
      GIN_TYPE("Vector", "Index"),
      GIN_DEFAULT(NODEF, "6"),
      GIN_DESC("Transmittance along path of downwelling radiation. A vector "
               "with the same length as *f_grid*.",
               "The version of FASTEM to use.")));
*/
/*
  md_data_raw.push_back(create_mdrecord(
      NAME("SurfaceFlatScalarReflectivity"),
      DESCRIPTION(R"--(Piecewise linear scalar surface reflectivity.

This method is similar to *surfaceFlatScalarReflectivity* but the
reflectivities are specified differently and Jacobian calculations
are supported. See *SurfaceFastem* for general remarks about
methods of SurfaceSomething type.

The method works with scalar reflectivies, i.e. it is assumed that
the reflection at vertical and horizontal polarisation is identical.
The scalar reflectivity is given at a number of frequencies,
specified by the GIN ``f_reflectivities``. The reflectivity at the
first frequency is denoted as \"Scalar reflectivity 0\" etc. Between
the frequencies in ``f_reflectivities``, the reflectivity is treated
to vary linearly. The reflectivity is assumed to be constant outside
of ``f_reflectivities``, and the end points in ``f_reflectivities`` can be
both inside and outside of the range of *f_grid*. Setting
``f_reflectivities`` to have a single value, implies that the reflectivity
is constant over *f_grid*.

For this method, ``surface_props_data`` must contain these data:
  \"Skin temperature\"
  \"Scalar reflectivity 0\"
  \"Scalar reflectivity 1\"
  ...
  \"Scalar reflectivity N\"
where N is the length of ``f_reflectivities``-1.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los",
          "surface_rmatrix",
          "dsurface_rmatrix_dx",
          "surface_emission",
          "dsurface_emission_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("dsurface_rmatrix_dx",
         "dsurface_emission_dx",
         
         "lat_grid",
         "lon_grid",
         "f_grid",
         "rtp_pos",
         "rtp_los",
         "specular_los",
         "surface_field",
         "surface_props_names",
         "dsurface_names",
         "jacobian_do"),
      GIN("f_reflectivities"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Frequency retrieval grid, see above.")));
      */
  /*
  md_data_raw.push_back(create_mdrecord(
      NAME("SurfaceTessem"),
      DESCRIPTION(R"--(TESSEM sea surface microwave emissivity parametrization.

This method allows to use TESSEM in retrievals requiring the
Jacobian. Otherwise as *surfaceTessem*. See *SurfaceFastem*
for general remarks about methods of SurfaceSomething type.

For this method, ``surface_props_data`` must contain these data:
  \"Water skin temperature\"
  \"Wind speed\"
  \"Salinity\
)--"),
      AUTHORS("Simon Pfreundschuh", "Patrick Eriksson"),
      OUT("surface_los",
          "surface_rmatrix",
          "dsurface_rmatrix_dx",
          "surface_emission",
          "dsurface_emission_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("dsurface_rmatrix_dx",
         "dsurface_emission_dx",
         
         "lat_grid",
         "lon_grid",
         "f_grid",
         "rtp_pos",
         "rtp_los",
         "tessem_neth",
         "tessem_netv",
         "surface_field",
         "surface_props_names",
         "dsurface_names",
         "jacobian_do"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
  md_data_raw.push_back(create_mdrecord(
      NAME("telsemStandalone"),
      DESCRIPTION(R"--(Stand-alone evaluation of the Telsem model.

This evaluates the Telsem land surface emissivity
model using the data from the provided atlas.

Since TELSEM atlases do not contain data for all locations
this function allows for nearest neighbor interpolation, which
can be enabled by setting the ``d_max`` GIN to a positive value.

This WSM throws a runtime error if the queried location is not
contained in the atlas or the distance of the neighboring cell
exceeds the given ``d_max`` value.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("emissivities"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("The computed h and v emissivites"),
      IN(),
      GIN("lat", "lon", "theta", "f", "ta", "d_max"),
      GIN_TYPE(
          "Numeric", "Numeric", "Numeric", "Vector", "TelsemAtlas", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, NODEF, "-1"),
      GIN_DESC("The latitude for which to compute the emissivities.",
               "The latitude for which to compute the emissivities.",
               "The incidence angle.",
               "The frequencies for which to compute the emissivities.",
               "The Telsem atlas to use.",
               "The maximum allowed distance for nearest neighbor"
               " interpolation in meters.")));
*/
  md_data_raw.push_back(create_mdrecord(
      NAME("telsemAtlasLookup"),
      DESCRIPTION(R"--(Lookup SSMI emissivities from Telsem Atlas.

This returns the emissivities (indices [0,..,6])
for the SSMI channels that are contained in
the Telsem atlas.

If given latitude and longitude are not in the atlas an empty
vector is returned.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("emissivities"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("The SSMI emissivities from the atlas"),
      IN(),
      GIN("lat", "lon", "atlas"),
      GIN_TYPE("Numeric", "Numeric", "TelsemAtlas"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("The latitude for which to compute the emissivities.",
               "The latitude for which to compute the emissivities.",
               "The Telsem atlas to use.")));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("telsemSurfaceTypeLandSea"),
      DESCRIPTION(R"--(TELSEM based land sea mask.

This method determines whether the position in *rtp_pos* is
of type ocean or land depending on whether a corresponding
cell is contained in the provided TELSEM atlas.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("surface_type"),
      GOUT_TYPE("Index"),
      GOUT_DESC("Surface type flag"),
      IN( "lat_grid", "lat_true", "lon_true", "rtp_pos"),
      GIN("atlas"),
      GIN_TYPE("TelsemAtlas"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The telsem atlas from which to lookup the surface type.")));
*/
  md_data_raw.push_back(create_mdrecord(
      NAME("telsem_atlasReadAscii"),
      DESCRIPTION(R"--(Reads single TELSEM atlas.

'directory' needs to contain the original 12 Telsem atlas files
and the correlations file. This WSM reads the atlas for the specified
month and stores the result in the provided output atlas.
)--"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("atlas"),
      GOUT_TYPE("TelsemAtlas"),
      GOUT_DESC("The atlas into which to store the loaded atlas."),
      IN(),
      GIN("directory", "month", "filename_pattern"),
      GIN_TYPE("String", "Index", "String"),
      GIN_DEFAULT(NODEF, NODEF, "ssmi_mean_emis_climato_@MM@_cov_interpol_M2"),
      GIN_DESC("Directory with TELSEM 2 SSMI atlas files.",
               "The month for which the atlas should be read.",
               "Filename pattern (@MM@ gets replaced by month number)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("telsem_atlasesReadAscii"),
      DESCRIPTION(R"--(Reads TELSEM atlas files.

'directory' needs to contain the original 12 Telsem atlas files
and the correlations file.
The whole data is combined into the WSV *telsem_atlases*
)--"),
      AUTHORS("Oliver Lemke"),
      OUT("telsem_atlases"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("directory", "filename_pattern"),
      GIN_TYPE("String", "String"),
      GIN_DEFAULT(NODEF, "ssmi_mean_emis_climato_@MM@_cov_interpol_M2"),
      GIN_DESC("Directory with TELSEM 2 SSMI atlas files.",
               "Filename pattern (@MM@ gets replaced by month number)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Tensor3Add"),
      DESCRIPTION(R"--(Adds a scalar value to all elements of a tensor3.

The result can either be stored in the same or another
variable.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor3"),
      GOUT_DESC("Output Tensor."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Tensor3", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Tensor.", "The value to be added to the tensor.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Tensor3ExtractFromTensor4"),
      DESCRIPTION(R"--(Extracts a Tensor3 from a Tensor4.

Copies book, page, row or column with given Index from input Tensor4
variable to output Tensor3.
Higher order equivalent of *VectorExtractFromMatrix*.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor3"),
      GOUT_DESC("Extracted tensor."),
      IN(),
      GIN("input", "i", "direction"),
      GIN_TYPE("Tensor4", "Index", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Input Tensor4.",
               "Index of book, page, row or column to extract.",
               "Direction. \"book\" or \"page\" or \"row\" or \"column\".")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Tensor3FromVector"),
      DESCRIPTION(R"--(Forms a Tensor3 of size nx1x1 from a vector of length n.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor3"),
      GOUT_DESC("Output tensor."),
      IN(),
      GIN("v"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Input vector.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Tensor3Multiply"),
      DESCRIPTION(R"--(Multiplies all elements of a tensor with the specified value.

The result can either be stored in the same or another
variable.
)--"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor3"),
      GOUT_DESC("Output Tensor."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Tensor3", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Tensor.",
               "The value to be multiplied with the tensor.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Tensor3SetConstant"),
      DESCRIPTION(R"--(Creates a tensor and sets all elements to the specified value.

The size is determined by *ncols*, *nrows* etc.
)--"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor3"),
      GOUT_DESC("Variable to initialize."),
      IN("npages", "nrows", "ncols"),
      GIN("value"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Tensor value.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Tensor4Add"),
      DESCRIPTION(R"--(Adds a scalar value to all elements of a tensor4.

The result can either be stored in the same or another
variable.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor4"),
      GOUT_DESC("Output Tensor."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Tensor4", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Tensor.", "The value to be added to the tensor.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Tensor4Multiply"),
      DESCRIPTION(R"--(Multiplies all elements of a tensor with the specified value.

The result can either be stored in the same or another
variable.
)--"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor4"),
      GOUT_DESC("Output Tensor."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Tensor4", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Tensor.",
               "The value to be multiplied with the tensor.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Tensor4SetConstant"),
      DESCRIPTION(R"--(Creates a tensor and sets all elements to the specified value.

The size is determined by *ncols*, *nrows* etc.
)--"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor4"),
      GOUT_DESC("Variable to initialize."),
      IN("nbooks", "npages", "nrows", "ncols"),
      GIN("value"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Tensor value.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Tensor5Multiply"),
      DESCRIPTION(R"--(Multiplies all elements of a tensor with the specified value.

The result can either be stored in the same or another
variable.
)--"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor5"),
      GOUT_DESC("Output Tensor."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Tensor5", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Tensor.",
               "The value to be multiplied with the tensor.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Tensor5SetConstant"),
      DESCRIPTION(R"--(Creates a tensor and sets all elements to the specified value.

The size is determined by *ncols*, *nrows* etc.
)--"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor5"),
      GOUT_DESC("Variable to initialize."),
      IN("nshelves", "nbooks", "npages", "nrows", "ncols"),
      GIN("value"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Tensor value.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Tensor6Multiply"),
      DESCRIPTION(R"--(Multiplies all elements of a tensor with the specified value.

The result can either be stored in the same or another
variable.
)--"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor6"),
      GOUT_DESC("Output Tensor."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Tensor6", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Tensor.",
               "The value to be multiplied with the tensor.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Tensor6SetConstant"),
      DESCRIPTION(R"--(Creates a tensor and sets all elements to the specified value.

The size is determined by *ncols*, *nrows* etc.
)--"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor6"),
      GOUT_DESC("Variable to initialize."),
      IN("nvitrines", "nshelves", "nbooks", "npages", "nrows", "ncols"),
      GIN("value"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Tensor value.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Tensor7Multiply"),
      DESCRIPTION(R"--(Multiplies all elements of a tensor with the specified value.

The result can either be stored in the same or another
variable.
)--"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor7"),
      GOUT_DESC("Output Tensor."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Tensor7", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Tensor.",
               "The value to be multiplied with the tensor.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Tensor7SetConstant"),
      DESCRIPTION(R"--(Creates a tensor and sets all elements to the specified value.

The size is determined by *ncols*, *nrows* etc.
)--"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Tensor7"),
      GOUT_DESC("Variable to initialize."),
      IN("nlibraries",
         "nvitrines",
         "nshelves",
         "nbooks",
         "npages",
         "nrows",
         "ncols"),
      GIN("value"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Tensor value.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("TestArrayOfAgenda"),
      DESCRIPTION(R"--(A method that is used for the TestArrayOfAgenda test case.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("test_agenda_array"),
      GIN("index"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("Index of agenda in array to execute.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("TestBasicGeodeticAccuracy"),
      DESCRIPTION(R"--(Tests the basic accuracy of the geodetic calculations.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("rte_pos"),
      GOUT("max_dl", "max_dpos", "max_dlos"),
      GOUT_TYPE("Numeric", "Vector", "Vector"),
      GOUT_DESC("Maximum error in term of distance.",
                "The maximum error for each position component.",
                "The maximum error for each LOS component."),
      IN("surface_field"),
      GIN("ntests","max_allowed_dl"),
      GIN_TYPE("Index","Numeric"),
      GIN_DEFAULT(NODEF,"0.1"),
      GIN_DESC("Number of tests.", "Maximum allowed error in term of distance.")));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("TessemNNReadAscii"),
      DESCRIPTION(R"--(Reads the initialization data for the TESSEM NeuralNet from an ASCII file.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("tessem_nn"),
      GOUT_TYPE("TessemNN"),
      GOUT_DESC("Tessem NeuralNet configuration."),
      IN(),
      GIN("filename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC(
          "NeuralNet parameters file as provided in the TESSEM 2 distribution.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("TestTessem"),
      DESCRIPTION(R"--(Example method for TESSEM2.

When using the default neural network parameter files
from the Tessem 2 distribution, the input Vector should contain
5 elements:

- Frequency (10-700) in GHz.
- Theta (0-90) Incidence angle in degrees.
- Windspeed (0-25) at 10m (m/s)
  Higher wind speed can be used, but without garantee.
- Surface skin temperature (270-310) in K.
- Salinity (0-0.04) in kg/kg
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("outvalues"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Tessem output emissivity."),
      IN(),
      GIN("net", "invalues"),
      GIN_TYPE("TessemNN", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Tessem NeuralNet parameters.", "Input data.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("time_gridOffset"),
      DESCRIPTION(R"--(Offsets a time grid by some seconds.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("time_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("time_grid"),
      GIN("dt"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Time in seconds to add")));

  md_data_raw.push_back(
      create_mdrecord(NAME("timerStart"),
               DESCRIPTION(R"--(Initializes the CPU timer.

Use *timerStop* to stop the timer.

Usage example:

 - timerStart
 - ReadXML(f_grid,\"frequencies.xml\")
 - timerStop
 - Print(timer)
)--"),
               AUTHORS("Oliver Lemke"),
               OUT("timer"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(
      create_mdrecord(NAME("timerStop"),
               DESCRIPTION(R"--(Stops the CPU timer.

See *timerStart* for example usage.
)--"),
               AUTHORS("Oliver Lemke"),
               OUT("timer"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("timer"),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("time_stampsSort"),
      DESCRIPTION(R"--(Sort ``input`` by *time_stamps* into ``output``.
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("ArrayOfTime,ArrayOfVector"),
      GOUT_DESC("Array sorted by time"),
      IN("time_stamps"),
      GIN("input"),
      GIN_TYPE("ArrayOfTime,ArrayOfVector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Array to sort of same size as *time_stamps*")));

  md_data_raw.push_back(create_mdrecord(
      NAME("TMatrixTest"),
      DESCRIPTION(R"--(T-Matrix validation test.

Executes the standard test included with the T-Matrix Fortran code.
Should give the same as running the tmatrix_lp executable in
3rdparty/tmatrix/.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("Touch"),
      DESCRIPTION(R"--(As *Ignore* but for agenda output.

This method is handy for use in agendas in order to suppress
warnings about not-produced output workspace variables.

What it does, in case the variable is initialized already, is:
Nothing!
In case the variable is not yet initialized, it is set to NaN.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("input"),
      GOUT_TYPE("Any"),
      GOUT_DESC("Variable to do nothing with."),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC(),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("transmittanceFromIy_aux"),
      DESCRIPTION(R"--(Creates a vector of transmittance values.

The transmittances are set based on optical depths in *iy_aux*. That is,
one of the quantities in *iy_aux* must be \"Optical depth\".

The created vector has a length matching *f_grid* and can e.g. be used
as input to some of the FASTEM methods.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("transmittance"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Created vector of transmittance values."),
      IN("iy_aux_vars", "iy_aux"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("Trapz"),
      DESCRIPTION(R"--(Intregrates a vector of over its grid range

The method integrates y(x) by the trapezoidal method.

The vector x is the positions where the integrand, y, is known.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Value of integral"),
      IN(),
      GIN("x", "y"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Grid.", "Integrand.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorAdd"),
      DESCRIPTION(R"--(Adds a scalar to all elements of a vector.

The result can either be stored in the same or another vector.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output Vector"),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Vector.", "The value to be added to the vector.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorAddElementwise"),
      DESCRIPTION(R"--(Element-wise addition of two vectors.

The method calculates c = a + b.

The variable ``b`` is allowed to have length 1, for any length of
``a``. This single value in ``b`` is then added to every element of ``a``.

The vectors ``a`` and ``c`` can be the same WSV, while ``b`` can not be
the same WSV as any of the the other vector.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("c"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output vector"),
      IN(),
      GIN("a", "b"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input vector.", "Vector to be added.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorClip"),
      DESCRIPTION(R"--(Clipping of a vector.

The input vector is copied to the output one (that can be same WSV)
but ensures that all values in ``output`` are inside the range [limit_low,
limit_high]. Where the input vector is below ``limit_low``, ``out`` is set
to ``limit_low``. And the same is performed with respect to ``limit_high``.
That is, the method works as *NumericClip* for each element of the
vector.

The method adopts the length of ``out`` when needed.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output vector."),
      IN(),
      GIN("input", "limit_low", "limit_high"),
      GIN_TYPE("Vector", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "-Inf", "Inf"),
      GIN_DESC("Input vector.",
               "Lower limit for clipping.",
               "Upper limit for clipping.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorCrop"),
      DESCRIPTION(R"--(Keeps only values of a vector inside the specified range.

All values outside the range [min_value,max-value] are removed.
Note the default values, that basically should act as -+Inf.

The result can either be stored in the same or another vector.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Cropped vector"),
      IN(),
      GIN("input", "min_value", "max_value"),
      GIN_TYPE("Vector", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "-99e99", "99e99"),
      GIN_DESC("Original vector",
               "Minimum value to keep",
               "Maximum value to keep")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorDivide"),
      DESCRIPTION(R"--(Divides all elements of a vector with the same value.

The result can either be stored in the same or another vector.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output Vector."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Vector.", "Denominator.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorDivideElementwise"),
      DESCRIPTION(R"--(Element-wise division of two vectors.

The method calculates c = a / b.

The variable ``b`` is allowed to have length 1, for any length of
``a``. This single value in ``b`` is then applied to every element of ``a``.

The vectors ``a`` and ``c`` can be the same WSV, while ``b`` can not be
the same WSV as any of the the other vector.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("c"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output vector"),
      IN(),
      GIN("a", "b"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input vector.", "Denominator Vector.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorExtractFromMatrix"),
      DESCRIPTION(R"--(Extracts a Vector from a Matrix.

Copies row or column with given Index from input Matrix variable
to create output Vector.
)--"),
      AUTHORS("Patrick Eriksson, Oliver Lemke, Stefan Buehler"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Extracted vector."),
      IN(),
      GIN("input", "i", "direction"),
      GIN_TYPE("Matrix", "Index", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Input matrix.",
               "Index of row or column.",
               "Direction. \"row\" or \"column\".")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorFlip"),
      DESCRIPTION(R"--(Flips a vector.

The output is the input vector in reversed order. The result can
either be stored in the same or another vector.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output vector."),
      IN(),
      GIN("input"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Input vector.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorGaussian"),
      DESCRIPTION(R"--(Fills a vector with a Gaussian function.

The width can be set in two ways, either by standard deviation or
the full width at half maximum. Only one of the corresponding GINs
can be >0 and that value will determine the width.

The vectors *x* and *y* can be the same variable.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("y"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output vector."),
      IN(),
      GIN("x", "x0", "si", "fwhm"),
      GIN_TYPE("Vector", "Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "0", "-1", "-1"),
      GIN_DESC("Grid of the function.",
               "Centre/mean point of the function.",
               "Standard deviation of the function, ignored if <=0.",
               "Full width at half-max of the function, ignored if <=0.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorInsertGridPoints"),
      DESCRIPTION(R"--(Insert some additional points into a grid.

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
)--"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("The new grid vector"),
      IN(),
      GIN("input", "points"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The original grid vector", "The points to insert")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorLinSpace"),
      DESCRIPTION(R"--(Initializes a vector with linear spacing.

The first element equals always the start value, and the spacing
equals always the step value, but the last value can deviate from
the stop value. ``step`` can be both positive and negative.

The created vector is [start, start+step, start+2*step, ...]
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output vector."),
      IN(),
      GIN("start", "stop", "step"),
      GIN_TYPE("Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Start value.",
               "Maximum/minimum value of the end value",
               "Spacing of the vector.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorLogSpace"),
      DESCRIPTION(R"--(Initializes a vector with logarithmic spacing.

The first element equals always the start value, and the spacing
equals always the step value, but note that the last value can 
deviate from the stop value. The keyword step can be both positive
and negative.

Note, that although start has to be given in direct coordinates,
step has to be given in log coordinates.

Explicitly, the vector is:
 exp([ln(start), ln(start)+step, ln(start)+2*step, ...])
)--"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("start", "stop", "step"),
      GIN_TYPE("Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("The start value. (Direct coordinates!)",
               "The maximum value of the end value. (Direct coordinates!)",
               "The spacing of the vector. (Log coordinates!)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorMatrixMultiply"),
      DESCRIPTION(R"--(Multiply a Vector with a Matrix and store the result in another
Vector.

This just computes the normal matrix-vector product, y=M*x. It is ok
if input and output Vector are the same.
)--"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT("y"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("The result of the multiplication (dimension m)."),
      IN(),
      GIN("M", "x"),
      GIN_TYPE("Matrix", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The Matrix to multiply (dimension m x n).",
               "The original Vector (dimension n).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorMultiply"),
      DESCRIPTION(R"--(Multiplies all elements of a vector with the same value.

The result can either be stored in the same or another vector.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output Vector."),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Vector.", "Scaling value.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorMultiplyElementwise"),
      DESCRIPTION(R"--(Element-wise multiplication of two vectors.

The method calculates c = a * b.

The variable ``b`` is allowed to have length 1, for any length of
``a``. This single value in ``b`` is then multiplied with every element
of ``a``.

The vectors ``a`` and ``c`` can be the same WSV, while ``b`` can not be
the same WSV as any of the the other vector.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("c"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output vector"),
      IN(),
      GIN("a", "b"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input vector.", "Multiplier.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorNLinSpace"),
      DESCRIPTION(R"--(Creates a vector with length *nelem*, equally spaced between the
given end values.

The length (*nelem*) must be larger than 1.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Variable to initialize."),
      IN("nelem"),
      GIN("start", "stop"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Start value.", "End value.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorNLinSpaceVector"),
      DESCRIPTION(R"--(As *VectorNLinSpace* but end points taken from a vector.

The method gives a vector with equidistant spacing between
first and last element of the reference vector.

The length (*nelem*) must be larger than 1.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Variable to initialize."),
      IN("nelem"),
      GIN("y"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Reference vector.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorNLogSpace"),
      DESCRIPTION(R"--(Creates a vector with length *nelem*, equally logarithmically
spaced between the given end values.

The length (*nelem*) must be larger than 1.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Variable to initialize."),
      IN("nelem"),
      GIN("start", "stop"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Start value.", "End value.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorPower"),
      DESCRIPTION(R"--(Calculates the power of each element in a vector.

The result can either be stored in the same or another vector.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output Vector."),
      IN(),
      GIN("input", "power"),
      GIN_TYPE("Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Vector.", "Power (exponent).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorReshapeMatrix"),
      DESCRIPTION(R"--(Converts a Matrix to a Vector.

The matrix is reshaped into a vector. That is, all elements of the matrix
are kept. The elements can be extracted both in column (default) and row
order. The output vector has the same length for both options.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Created vector."),
      IN(),
      GIN("input", "direction"),
      GIN_TYPE("Matrix", "String"),
      GIN_DEFAULT(NODEF, "column"),
      GIN_DESC("Input matrix.", "Direction. \"row\" or \"column\".")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorSetConstant"),
      DESCRIPTION(R"--(Creates a vector and sets all elements to the specified value.

The vector length is determined by *nelem*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Variable to initialize."),
      IN("nelem"),
      GIN("value"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Vector value.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorSparseMultiply"),
      DESCRIPTION(R"--(Multiply a Vector with a Sparse and store the result in another
Vector.

This just computes the normal matrix-vector product, y=M*x, with
m being a Sparse. It is ok if input and output Vector are the same.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("y"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("The result of the multiplication (dimension m)."),
      IN(),
      GIN("M", "x"),
      GIN_TYPE("Sparse", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The Sparse to multiply (dimension m x n).",
               "The original Vector (dimension n).")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorSubtract"),
      DESCRIPTION(R"--(Subtracts a scalar from all elements of a vector.

The result can either be stored in the same or another vector.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("output"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output Vector"),
      IN(),
      GIN("input", "value"),
      GIN_TYPE("Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Vector.", "The value to be subtracted from the vector.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("VectorSubtractElementwise"),
      DESCRIPTION(R"--(Element-wise subtraction of two vectors.

The method calculates c = a - b.

The variable ``b`` is allowed to have length 1, for any length of
``a``. This single value in ``b`` is then subtracted to every element of ``a``.

The vectors ``a`` and ``c`` can be the same WSV, while ``b`` can not be
the same WSV as any of the the other vector.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("c"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output vector"),
      IN(),
      GIN("a", "b"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input vector.", "Vector to be subtracted.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("water_p_eq_fieldMK05"),
      DESCRIPTION(R"--(Calculates *water_p_eq_field* according to Murphy and Koop, 2005.

Default is setting the saturation pressure to the one with respect
to water at temperatures >= 0C, and to the one with respect to ice
for <0C. The GIN ``only_liquid`` allows you to apply the liquid value
at all temperatures.

The saturation pressure with respect to liquid and ice water is
calculated according to Eq. 10 and 7, respectively, of:
Murphy, D. M., & Koop, T. (2005). Review of the vapour pressures of
ice and supercooled water for atmospheric applications. Quarterly
Journal of the Royal Meteorological Society, 131(608), 1539-1565.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("water_p_eq_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field"),
      GIN("only_liquid"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("Set to 1 to use liquid saturation pressure at all temperatures.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Wigner6Init"),
      DESCRIPTION(R"--(Initialize the wigner 3 and 6 tables

The default values take about 1 Gb memory.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("wigner_initialized"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("fast_wigner_stored_symbols", "largest_wigner_symbol_parameter"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT("20000000", "250"),
      GIN_DESC(
          "Number of stored symbols possible before replacements",
          "Largest symbol used for initializing factorials (e.g., largest J or L)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("Wigner3Init"),
      DESCRIPTION(R"--(Initialize the wigner 3 tables

The default values take about 400 Mb memory.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("wigner_initialized"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("fast_wigner_stored_symbols", "largest_wigner_symbol_parameter"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT("20000000", "250"),
      GIN_DESC(
          "Number of stored symbols possible before replacements",
          "Largest symbol used for initializing factorials (e.g., largest J or L)")));

  md_data_raw.push_back(
      create_mdrecord(NAME("Wigner6Unload"),
               DESCRIPTION(R"--(Unloads the wigner 3 and 6 tables
)--"),
               AUTHORS("Richard Larsson"),
               OUT("wigner_initialized"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("wigner_initialized"),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(NAME("Wigner3Unload"),
                                 DESCRIPTION(R"--(Unloads the wigner 3 tables
)--"),
                                 AUTHORS("Richard Larsson"),
                                 OUT("wigner_initialized"),
                                 GOUT(),
                                 GOUT_TYPE(),
                                 GOUT_DESC(),
                                 IN("wigner_initialized"),
                                 GIN(),
                                 GIN_TYPE(),
                                 GIN_DEFAULT(),
                                 GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("WignerFastInfoPrint"),
      DESCRIPTION(R"--(Prints the fast wigner table information if compiled with this option
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("wigner_initialized"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("wind_u_fieldIncludePlanetRotation"),
      DESCRIPTION(R"--(Maps the planet's rotation to an imaginary wind.

This method is of relevance if the observation platform is not
following the planet's rotation, and Doppler effects must be
considered. Examples include full disk observations from another
planet or a satellite not in orbit of the observed planet.

The rotation of the planet is not causing any Doppler shift for
1D and 2D simulations, and the method can only be used for 3D.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("atm_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field", "surface_field",
         "planet_rotation_period"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("WMRFSelectChannels"),
      DESCRIPTION(R"--(Select some channels for WMRF calculation.

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
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("f_grid", "wmrf_weights", "f_backend"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid", "f_backend", "wmrf_weights", "wmrf_channels"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("WriteMolTau"),
      DESCRIPTION(R"--(Writes a 'molecular_tau_file' as required for libRadtran.

The libRadtran (www.libradtran.org) radiative transfer package is a 
comprehensive package for various applications, it can be used to 
compute radiances, irradiances, actinic fluxes, ... for the solar 
and the thermal spectral ranges. Absorption is usually treated using 
k-distributions or other parameterizations. For calculations with high 
spectral resolution it requires absorption coefficients from an external 
line-by-line model. Using this method, arts generates a file that can be 
used by libRadtran (option molecular_tau_file).

)--"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid", "atm_field", "propmat_clearsky_field"),
      GIN("filename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Name of the ``molecular_tau_file``.")));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("WriteNetCDF"),
      DESCRIPTION(R"--(Writes a workspace variable to a NetCDF file.

This method can write variables of limited groups.

If the filename is omitted, the variable is written
to <basename>.<variable_name>.nc.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("input", "filename"),
      GIN_TYPE("Vector, Matrix, Tensor3, Tensor4, Tensor5, ArrayOfVector,"
               "ArrayOfIndex, ArrayOfMatrix, GasAbsLookup",
               "String"),
      GIN_DEFAULT(NODEF, ""),
      GIN_DESC("Variable to be saved.", "Name of the NetCDF file."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("WriteNetCDFIndexed"),
      DESCRIPTION(R"--(As *WriteNetCDF*, but creates indexed file names.

This method can write variables of any group.

If the filename is omitted, the variable is written
to <basename>.<variable_name>.nc.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("file_index"),
      GIN("input", "filename"),
      GIN_TYPE("Vector, Matrix, Tensor3, Tensor4, Tensor5, ArrayOfVector,"
               "ArrayOfMatrix, GasAbsLookup",
               "String"),
      GIN_DEFAULT(NODEF, ""),
      GIN_DESC("Variable to be saved.", "Name of the NetCDF file."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("WriteBuiltinPartitionFunctionsXML"),
      DESCRIPTION(R"--(Writes all the builtin partition functions to file.

All available partition functions are written to files in the select format
in the select directory

The temperature will be linearly spaced between [Tlow, Tupp] with N values
)--"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("output_file_format"),
      GIN("dir", "Tlow", "Tupp", "N"),
      GIN_TYPE("String", "Numeric", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF),
      GIN_DESC("The directory to write the data towards",
               "The lowest temperature",
               "The highest temperature",
               "The number of temperature points")));

  md_data_raw.push_back(create_mdrecord(
      NAME("WriteXML"),
      DESCRIPTION(R"--(Writes a workspace variable to an XML file.

This method can write variables of any group.

If the filename is omitted, the variable is written
to <basename>.<variable_name>.xml.
If no_clobber is set to 1, an increasing number will be
appended to the filename if the file already exists.
)--"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("output_file_format"),
      GIN("input", "filename", "no_clobber"),
      GIN_TYPE("Any", "String", "Index"),
      GIN_DEFAULT(NODEF, "", "0"),
      GIN_DESC("Variable to be saved.",
               "Name of the XML file.",
               "0: Overwrite existing files, 1: Use unique filenames"),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("WriteXMLIndexed"),
      DESCRIPTION(R"--(As *WriteXML*, but creates indexed file names.

The variable is written to a file with name::

  <filename>.<file_index>.xml.

where <file_index> is the value of ``file_index``.

This means that ``filename`` shall here not include the .xml
extension. Omitting filename works as for *WriteXML*.
)--"),
      AUTHORS("Patrick Eriksson, Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("output_file_format", "file_index"),
      GIN("input", "filename", "digits"),
      GIN_TYPE("Any", "String", "Index"),
      GIN_DEFAULT(NODEF, "", "0"),
      GIN_DESC(
          "Workspace variable to be saved.",
          "File name. See above.",
          "Equalize the widths of all numbers by padding with zeros as necessary. "
          "0 means no padding (default)."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(create_mdrecord(
      NAME("xaStandard"),
      DESCRIPTION(R"--(Standard function for creating *xa*.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("xa"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "atmfields_checked",
         "atmgeom_checked",
         "atm_field",
         "abs_species",
         "cloudbox_on",
         "cloudbox_checked",
         "particle_bulkprop_names",
         "surface_field",
         "surface_props_names",
         "water_p_eq_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("xClip"),
      DESCRIPTION(R"--(Clipping of the state vector.

The method allows you to apply hard limits the values of a
retrieval quantity. The retrieval quantity is specified by
``ijq``. All values of the quantity below ``limit_low``, are simply
set to ``limit_low``. And the same is performed with respect to
``limit_high``. That is, the data in x for the retrieval quantity
are forced to be inside the range [limit_low,limit_high].

Setting ijq=-1, is a shortcut for applying the limits on all
retrieval quantities.

Notice that limits must be specified in the unit used in *x*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("x"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("x", "jacobian_quantities"),
      GIN("ijq", "limit_low", "limit_high"),
      GIN_TYPE("Index", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "-Inf", "Inf"),
      GIN_DESC("Retrieval quantity index (zero-based)",
               "Lower limit for clipping.",
               "Upper limit for clipping.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("x2artsAtmAndSurf"),
      DESCRIPTION(R"--(Maps *x* to atmospheric and surface variables.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("atm_field",
          "surface_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atm_field",
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
         "water_p_eq_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("x2artsSensor"),
      DESCRIPTION(R"--(Maps *x* to sensor variables.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_los",
          "f_backend",
          "y_baseline",
          "sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid",
          "sensor_response_pol_grid",
          "sensor_response_dlos_grid",
          "mblock_dlos"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_los",
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
         "sensor_time"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("x2artsSpectroscopy"),
      DESCRIPTION(R"--(Just defined to indicate a future extensiom.

Don't call the method, it will just generate an error.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("yApplySensorPol"),
      DESCRIPTION(R"--(Extraction of arbitrary linear polarisation.

This method shall be called after *yCalc* and then applies *sensor_pol*
on the output of *yCalc*. See *sensor_pol* for definition of the
polarisation responses. The *sensor_response* given to *yCalc* can not
contain any polarisation response, it must maintain original Stokes
elements. The value of ``stokes_dim`` must be >= 3.

The values in *sensor_pol* are applied on *y*, and *jacobian* if relevant.
*y_pol* is set following the values in *sensor_pol* but is rounded to
an integer value. Remaining data associated with *y* (e.g. y_pos) are
set to the value matching the first Stokes element.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("y", "y_f", "y_pol", "y_pos", "y_los", "y_aux", "y_geo", "jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("y",
         "y_f",
         "y_pol",
         "y_pos",
         "y_los",
         "y_aux",
         "y_geo",
         "jacobian",
         
         "jacobian_do",
         "sensor_pos",
         "sensor_pol"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("yApplyUnit"),
      DESCRIPTION(R"--(Conversion of *y* to other spectral units.

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

If you are using this method, *iy_unit* should be set to \"1\" when
calling *yCalc*, and be changed before calling this method.

Conversion of *y_aux* is not supported.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("y", "jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("y", "jacobian", "y_f", "y_pol", "iy_unit"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
  
  md_data_raw.push_back(create_mdrecord(
      NAME("ybatchColdAtmHotAtmCycle"),
      DESCRIPTION(R"--(Computes *ybatch* from input using standard calibration scheme of
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
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ybatch", "sensor_time"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("level0_data", "level0_time"),
      GIN("cold_temp", "hot_temp", "first_c_index"),
      GIN_TYPE("Vector", "Vector", "Index"),
      GIN_DEFAULT(NODEF, NODEF, "0"),
      GIN_DESC(
        "Cold load calibration temperature (must match level0_data length)",
        "Hot load calibration temperature (must match level0_data length)",
        "Index offset of the first cold position")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ybatchCalc"),
      DESCRIPTION(R"--(Performs batch calculations for the measurement vector y.

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
)--"),
      AUTHORS("Stefan Buehler"),
      OUT("ybatch", "ybatch_aux", "ybatch_jacobians"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ybatch_start", "ybatch_n", "ybatch_calc_agenda"),
      GIN("robust"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("A flag with value 1 or 0. If set to one, the batch "
               "calculation will continue, even if individual jobs fail. In "
               "that case, a warning message is written to screen and file "
               "(out1 output stream), and the *y* Vector entry for the "
               "failed job in *ybatch* is left empty.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("yColdAtmHot"),
      DESCRIPTION(R"--(Computes *y* from input using standard calibration scheme of
cold-atm-hot observations

If calib evaluates as true:
    y = cold_temp + (hot_temp - cold_temp) * (atm - cold) / (hot - cold)

If calib evaluates as false:
    y = (hot_temp * cold - cold_temp * hot) / (hot - cold)
)--"),
      AUTHORS("Richard Larsson"),
      OUT("y"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("cold", "atm", "hot", "cold_temp", "hot_temp", "calib"),
      GIN_TYPE("Vector", "Vector", "Vector", "Numeric", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, NODEF, "1"),
      GIN_DESC("N-elem Vector of cold load linear power",
               "N-elem Vector of atmosphere linear power",
               "N-elem Vector of hot load linear power",
               "Cold load temperature",
               "Hot load temperature",
               "Flag for calibration scheme, false means system temperature is computed")));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("ybatchMetProfiles"),
      DESCRIPTION(R"--(This method is used for simulating ARTS for metoffice model fields

This method reads in *met_amsu_data* which contains the
lat-lon of the metoffice profile files as a Matrix. It then
loops over the number of profiles and corresponding to each
longitude create the appropriate profile basename. Then,
corresponding to each basename we have temperature field, altitude
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
)--"),
      AUTHORS("Sreerekha T.R."),
      OUT("ybatch"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species",
         "met_profile_calc_agenda",
         "f_grid",
         "met_amsu_data",
         "sensor_pos",
         "lat_grid",
         "lon_grid",
         "scat_data"),
      GIN("nelem_z_grid", "met_profile_path", "met_profile_pnd_path"),
      GIN_TYPE("Index", "String", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("FIXME DOC", "FIXME DOC", "FIXME DOC")));
*/
  md_data_raw.push_back(create_mdrecord(
      NAME("ybatchMetProfilesClear"),
      DESCRIPTION(R"--(This method is used for simulating ARTS for metoffice model fields
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
)--"),
      AUTHORS("Seerekha T.R."),
      OUT("ybatch"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species",
         "met_profile_calc_agenda",
         "f_grid",
         "met_amsu_data",
         "sensor_pos",
         "surface_field"),
      GIN("nelem_p_grid", "met_profile_path"),
      GIN_TYPE("Index", "String"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("FIXME DOC", "FIXME DOC")));


  md_data_raw.push_back(create_mdrecord(
      NAME("ybatchTimeAveraging"),
      DESCRIPTION(R"--(Time average of *ybatch* and *sensor_time*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ybatch", "sensor_time"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ybatch", "sensor_time"),
      GIN("time_step", "disregard_first", "disregard_last"),
      GIN_TYPE("String", "Index", "Index"),
      GIN_DEFAULT(NODEF, "0", "0"),
      GIN_DESC("Time step in the form \"INDEX SCALE\", where SCALE is \"h\", \"min\", or \"s\" for hours, minutes or seconds",
               "Flag to remove first time step (e.g., if it is an incomplete step)",
               "Flag to remove last time step (e.g., if it is an incomplete step)")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ybatchTroposphericCorrectionNaiveMedianForward"),
      DESCRIPTION(R"--(Performs naive tropospheric corrections on *ybatch*

Sets *ybatch_corr* to be able to perform the inverse of the corrections,
each array-element with 3 entries as [median, part_trans, trop_temp]

Uses the same tropospheric temperature for all values if trop_temp.nelem()==1
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ybatch_corr", "ybatch"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ybatch"),
      GIN("range", "trop_temp", "targ_temp"),
      GIN_TYPE("ArrayOfIndex", "Vector", "Numeric"),
      GIN_DEFAULT("ArrayOfIndex(0)", NODEF, "2.73"),
      GIN_DESC("Positions where the median of the baseline is computed, if empty all is used",
               "Radiative temperature of the troposphere [dim: 1 or ybatch.nelem()]",
               "Temperature target of the baseline")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ybatchTroposphericCorrectionNaiveMedianInverse"),
      DESCRIPTION(R"--(Performs inverse of naive tropospheric corrections on *ybatch*
)--"),
      AUTHORS("Richard Larsson"),
      OUT("ybatch"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ybatch", "ybatch_corr"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("yCalc"),
      DESCRIPTION(R"--(Calculation of complete measurement vectors (y).

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("y", "y_f", "y_pol", "y_pos", "y_los", "y_aux", "y_geo", "jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmgeom_checked",
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
         "iy_aux_vars"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("yCalcAppend"),
      DESCRIPTION(R"--(Replaces *yCalc* if a measurement shall be appended to an
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
appended \"blindly\" in *y_aux*. That is, data of different type
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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("y",
          "y_f",
          "y_pol",
          "y_pos",
          "y_los",
          "y_aux",
          "y_geo",
          "jacobian",
          "jacobian_quantities"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("y",
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
         "iy_aux_vars"),
      GIN("jacobian_quantities_copy", "append_instrument_wfs"),
      GIN_TYPE("ArrayOfRetrievalQuantity", "Index"),
      GIN_DEFAULT(NODEF, "0"),
      GIN_DESC("Copy of *jacobian_quantities* of first measurement.",
               "Flag controlling if instrumental weighting functions are "
               "appended or treated as different retrieval quantities.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("yRadar"),
      DESCRIPTION(R"--(Replaces *yCalc* for radar/lidar calculations.

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

- ``\"1\"``:
    Backscatter coefficient. Unit is 1/(m*sr). At zero
    attenuation, this equals the scattering matrix value for
    the backward direction. See further AUG.
- ``\"Ze\"``: Equivalent reflectivity. Unit is mm^6/m^3. Conversion formula is given below.
- ``\"dBZe\"``: 10*log10(Ze/Z0), where Z0 is 1 mm^6/m^3.

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
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("y", "y_f", "y_pol", "y_pos", "y_los", "y_aux", "y_geo", "jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmgeom_checked",
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
         "range_bins"),
      GIN("ze_tref", "k2", "dbze_min"),
      GIN_TYPE("Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT("273.15", "-1", "-99"),
      GIN_DESC("Reference temperature for conversion to Ze.",
               "Reference dielectric factor.",
               "Clip value for dBZe.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("ySimpleSpectrometer"),
      DESCRIPTION(R"--(Converts *iy* to *y* assuming a fixed frequency resolution.

This is a short-cut, avoiding *yCalc*, that can be used to convert
monochromatic pencil beam data to spectra with a fixed resolution.

The method mimics a spectrometer with rectangular response
functions, all having the same width (``df``). The position of
the first spectrometer channel is set to f_grid[0] + df / 2.
The centre frequency of channels are returned as *y_f*.

Auxiliary variables and *jacobian* s are not handled.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("y", "y_f"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy",  "f_grid"),
      GIN("df"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Selected frequency resolution.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("y_geo_seriesFromY_geo"),
      DESCRIPTION(R"--(Fills *y_geo_series* with data from *y_geo*.

The geo-position is taken from the first channel. There is no check
that the other channels have identical data in *y_geo*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("y_geo_series"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("y_geo", "sensor_response_f_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(create_mdrecord(
      NAME("y_geo_swathFromY_geo"),
      DESCRIPTION(R"--(Fills *y_geo_series* with data from *y_geo*.

The geo-position is taken from the first channel. There is no check
that the other channels have identical data in *y_geo*.

The method assumes the same order in *y* as *y_swathFromY*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("y_geo_swath"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("y_geo", "sensor_response_f_grid"),
      GIN("npixel"),
      GIN_TYPE("Index"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Number of pixels per swath.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("y_seriesFromY"),
      DESCRIPTION(R"--(Fills *y_series* with data from *y*.

The method basically reshapes *y* to fit *y_series*.

Default is to check that *y_f* does not change between posistions,
i.e. that the channel frequencies do not vary.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("y_series"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("y", "y_f", "sensor_response_f_grid"),
      GIN("safe"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Flag for checking that channels do not vary in frequency.")));

  md_data_raw.push_back(create_mdrecord(
      NAME("y_swathFromY"),
      DESCRIPTION(R"--(Fills *y_swath* with data from *y*.

The method basically reshapes *y* to fit *y_swath*. It is assumed
that swath forms the outermost loop in *y*. That is, first in *y*
are the data for the first swath etc. The number of pixels per swath
must be specified manually by a GIN parameter.

To set *sensor_pos* and *sensor_los* having data organised in swath
format, use *MatrixReshapeTensor3*.

Default is to check that *y_f* does not change between posistions,
i.e. that the channel frequencies do not vary.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("y_swath"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("y", "y_f", "sensor_response_f_grid"),
      GIN("npixel", "safe"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT(NODEF, "1"),
      GIN_DESC("Number of pixels per swath.",
               "Flag for checking that channels do not vary in frequency.")));

/*
  md_data_raw.push_back(create_mdrecord(
      NAME("z_fieldFromHSE"),
      DESCRIPTION(R"--(Force altitudes to fulfil hydrostatic equilibrium.

The method applies hydrostatic equilibrium. A mixture of \"dry
air\" and water vapour (if present as *abs_species* tag) is assumed.
That is, the air is assumed to be well mixed and its weight, apart
from the water vapour, is constant (*molarmass_dry_air*). In
addition, the effect of any particles (including liquid and ice
particles) is neglected.

The output is an update of ``z_field``. This variable is expected to
contain approximative altitudes when calling the function. The
altitude matching *p_hse* is kept constant. Other input altitudes can
basically be arbitrary, but good estimates give quicker calculations.

The calculations are repeated until the change in altitude is below
*z_hse_accuracy*. An iterative process is needed as gravity varies
with altitude.

For 1D and 2D, the geographical position is taken from *lat_true*
and *lon_true*.
)--"),
      AUTHORS("Patrick Eriksson"),
      OUT("z_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(
         "p_grid",
         "lat_grid",
         "lon_grid",
         "lat_true",
         "lon_true",
         "abs_species",
         "atm_field",
         "surface_field",
         "atmfields_checked",
         "g0_agenda",
         "molarmass_dry_air",
         "p_hse",
         "z_hse_accuracy"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
*/

  md_data_raw.push_back(create_mdrecord(
      NAME("spectral_radiance_fieldPlaneParallelSpectralRadianceOperator"),
      DESCRIPTION(R"--(Create a *spectral_radiance_field*

This is an experimental solution.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("spectral_radiance_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("spectral_radiance_profile_operator", "f_grid", "za_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  /*
  md_data_raw.push_back(create_mdrecord(
      NAME("spectral_radiance_profile_operatorPlaneParallel"),
      DESCRIPTION(R"--(Create a radiance profile operator

This is an experimental solution.
)--"),
      AUTHORS("Richard Larsson"),
      OUT("spectral_radiance_profile_operator"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("z_field",
         "ppath_lmax",
         "atmosphere_dim",
         "p_grid",
         "t_field",
         "nlte_field",
         "vmr_field",
         "wind_u_field",
         "wind_v_field",
         "wind_w_field",
         "mag_u_field",
         "mag_v_field",
         "mag_w_field",
         "abs_species",
         "predefined_model_data",
         "abs_cia_data",
         "xsec_fit_data",
         "isotopologue_ratios",
         "abs_lines_per_species"),
      GIN("T_extrapolfac", "ignore_errors"),
      GIN_TYPE("Numeric", "Index"),
      GIN_DEFAULT("0.5", "0"),
      GIN_DESC(
          "Temperature extrapolation factor (relative to grid spacing).",
          "Set to 1 to suppress runtime errors (and return NAN values instead)."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(false),
      PASSWORKSPACE(true),
      PASSWSVNAMES(false)));
*/
  md_data_raw.push_back(
      create_mdrecord(NAME("PlanetSet"),
                      DESCRIPTION(R"--(Sets *g0_agenda*, ``refellipsoid``, *molarmass_dry_air*, and *planet_rotation_period* to default values

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
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("g0_agenda", "surface_field", "molarmass_dry_air", "planet_rotation_period"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  //! Add all the agenda-setting methods below here:

  md_data_raw.push_back(
      create_mdrecord(NAME("dobatch_calc_agendaSet"),
                      DESCRIPTION(R"--(Sets *dobatch_calc_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("dobatch_calc_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("doit_conv_test_agendaSet"),
                      DESCRIPTION(R"--(Sets *doit_conv_test_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("doit_conv_test_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("doit_mono_agendaSet"),
                      DESCRIPTION(R"--(Sets *doit_mono_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("doit_mono_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("doit_rte_agendaSet"),
                      DESCRIPTION(R"--(Sets *doit_rte_agenda* to a default value

Options are:

- There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("doit_rte_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("doit_scat_field_agendaSet"),
                      DESCRIPTION(R"--(Sets *doit_scat_field_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("doit_scat_field_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("forloop_agendaSet"),
                      DESCRIPTION(R"--(Sets *forloop_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("forloop_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("g0_agendaSet"),
                      DESCRIPTION(R"--(Sets *g0_agenda* to a default value

Options are:

- ``"Earth"``: Uses *g0Earth* to set *g0*
- ``"Io"``: Uses *g0Io* to set *g0*
- ``"Jupiter"``: Uses *g0Jupiter* to set *g0*
- ``"Mars"``: Uses *g0Mars* to set *g0*
- ``"Venus"``: Uses *g0Venus* to set *g0*
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("g0_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("gas_scattering_agendaSet"),
                      DESCRIPTION(R"--(Sets *gas_scattering_agenda* to a default value

Options are:

- ``"Dummy"``:

    1. Will *Ignore* all agenda inputs
    2. Uses *Touch* on all agenda outputs
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("gas_scattering_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT("Dummy"),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("inversion_iterate_agendaSet"),
                      DESCRIPTION(R"--(Sets *inversion_iterate_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("inversion_iterate_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("iy_cloudbox_agendaSet"),
                      DESCRIPTION(R"--(Sets *iy_cloudbox_agenda* to a default value

Options are:

- ``"LinInterpField"``: Uses ``iyInterpCloudboxField`` to set *iy*
- ``"QuarticInterpField"``: Uses ``iyInterpCloudboxField`` to set *iy* using ``za_interp_order=4``
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("iy_cloudbox_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("iy_independent_beam_approx_agendaSet"),
                      DESCRIPTION(R"--(Sets *iy_independent_beam_approx_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("iy_independent_beam_approx_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("iy_loop_freqs_agendaSet"),
                      DESCRIPTION(R"--(Sets *iy_loop_freqs_agenda* to a default value

Options are:

- ``"Emission"``:

    1. Uses ``ppathCalc`` to set *ppath*
    2. Uses ``iyEmissionStandard`` to set *iy*, *iy_aux*, ``ppvar_p``, ``ppvar_t``, ``ppvar_nlte``, ``ppvar_vmr``, ``ppvar_wind``, ``ppvar_mag``, *ppvar_f*, *ppvar_iy*, *ppvar_trans_cumulat*, and *ppvar_trans_partial*, and also  to modify *diy_dx*

- ``"Transmission"``:

    1. Uses ``ppathCalc`` to set *ppath*
    2. Uses *iyTransmissionStandard* to set *iy*, *iy_aux*, ``ppvar_p``, ``ppvar_t``, ``ppvar_nlte``, ``ppvar_vmr``, ``ppvar_wind``, ``ppvar_mag``, *ppvar_f*, *ppvar_iy*, *ppvar_trans_cumulat*, and *ppvar_trans_partial*, and also to modify *diy_dx*
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("iy_loop_freqs_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("iy_main_agendaSet"),
                      DESCRIPTION(R"--(Sets *iy_main_agenda* to a default value

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
    2. Uses *iyClearsky* to set *iy*, *iy_aux*, ``ppvar_p``, ``ppvar_t``, ``ppvar_nlte``, ``ppvar_vmr``, ``ppvar_wind``, ``ppvar_mag``, *ppvar_f*, *ppvar_iy*, *ppvar_trans_cumulat*, and *ppvar_trans_partial*, and also  to modify *diy_dx*
    3. Sets *geo_pos* to empty

- ``"Transmission"``:

    1. Uses ``ppathCalc`` to set *ppath* using *cloudbox_on* = 0
    2. Uses *iyTransmissionStandard* to set *iy*, *iy_aux*, ``ppvar_p``, ``ppvar_t``, ``ppvar_nlte``, ``ppvar_vmr``, ``ppvar_wind``, ``ppvar_mag``, *ppvar_f*, *ppvar_iy*, *ppvar_trans_cumulat*, and *ppvar_trans_partial*, and also to modify *diy_dx*
    3. Sets *geo_pos* to empty

- ``"TransmissionUnitUnpolIntensity"``:

    1. Uses *MatrixUnitIntensity* using out = *iy_transmitter*, and f =* f_grid*
    2. Uses ``ppathCalc`` to set *ppath* using *cloudbox_on* = 0
    3. Uses *iyTransmissionStandard* to set *iy*, *iy_aux*, ``ppvar_p``, ``ppvar_t``, ``ppvar_nlte``, ``ppvar_vmr``, ``ppvar_wind``, ``ppvar_mag``, *ppvar_f*, *ppvar_iy*, *ppvar_trans_cumulat*, and *ppvar_trans_partial*, and also to modify *diy_dx*
    4. Sets *geo_pos* to empty

- ``"TransmissionUnitPolIntensity"``:

    1. Uses *iy_transmitterSinglePol* to set *iy_transmitter*
    2. Uses ``ppathCalc`` to set *ppath* using *cloudbox_on* = 0
    3. Uses *iyTransmissionStandard* to set *iy*, *iy_aux*, ``ppvar_p``, ``ppvar_t``, ``ppvar_nlte``, ``ppvar_vmr``, ``ppvar_wind``, ``ppvar_mag``, *ppvar_f*, *ppvar_iy*, *ppvar_trans_cumulat*, and *ppvar_trans_partial*, and also to modify *diy_dx*
    4. Sets *geo_pos* to empty

- ``"Freqloop"``:

    1. Uses *iyLoopFrequencies* to set *iy*, *iy_aux*, *ppath*, and *diy_dx*
    2. Sets *geo_pos* to empty
    3. Will *Ignore* the *diy_dx* agenda input

- ``"ScattMC"``:

    1. Uses *iyMC* to set *iy*, *iy_aux*, and *diy_dx*
    2. Sets *geo_pos* to empty
    3. Will *Ignore* the *diy_dx* agenda input
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("iy_main_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("iy_radar_agendaSet"),
                      DESCRIPTION(R"--(Sets *iy_radar_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("iy_radar_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("iy_space_agendaSet"),
                      DESCRIPTION(R"--(Sets *iy_space_agenda* to a default value

Options are:

- ``"CosmicBackground"``:

    1. Uses *MatrixCBR* using out = *iy*, and f = *f_grid*
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("iy_space_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT("CosmicBackground"),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("iy_surface_agendaSet"),
                      DESCRIPTION(R"--(Sets *iy_space_agenda* to a default value

Options are:

- ``"UseSurfaceRtprop"``:

    1. Uses *SurfaceDummy* to modify *dsurface_rmatrix_dx*, and *dsurface_emission_dx*
    2. Uses *iySurfaceRtpropAgenda* to set *iy*, *surface_skin_t*, *surface_los*, *surface_rmatrix*, and *surface_emission*, and also to modify *diy_dx*
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("iy_surface_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT("UseSurfaceRtprop"),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("jacobian_agendaSet"),
                      DESCRIPTION(R"--(Sets *jacobian_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("jacobian_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("main_agendaSet"),
                      DESCRIPTION(R"--(Calling this is - and always will be - an error
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("main_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("met_profile_calc_agendaSet"),
                      DESCRIPTION(R"--(Sets *met_profile_calc_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("met_profile_calc_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("pha_mat_spt_agendaSet"),
                      DESCRIPTION(R"--(Sets *pha_mat_spt_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("pha_mat_spt_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("ppath_agendaSet"),
                      DESCRIPTION(R"--(Sets *ppath_agenda* to a default value

Options are:

- ``"FollowSensorLosPath"``:

    1. Uses ``ppathStepByStep`` to set *ppath*

- ``"PlaneParallel"``:

    1. Uses ``ppathPlaneParallel`` to set *ppath*

- ``"TransmitterReceiverPath"``:

    1. Uses ``rte_losGeometricFromRtePosToRtePos2`` to set *rte_los*
    2. Uses ``ppathFromRtePos2`` to set *ppath*, and also to modify *rte_los*, and *ppath_lraytrace*
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("ppath_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("ppath_step_agendaSet"),
                      DESCRIPTION(R"--(Sets *ppath_step_agenda* to a default value

Options are:

- ``"GeometricPath"``:

    1. Uses ``ppath_stepGeometric`` to modify *ppath*
- ``"RefractedPath"``:

    1. Uses ``ppath_stepRefractionBasic`` to modify *ppath*
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("ppath_step_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("ppvar_rtprop_agendaSet"),
                      DESCRIPTION(R"--(Sets *ppvar_rtprop_agenda* to a default value

Options are:
    FIXME
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("ppvar_rtprop_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("propmat_clearsky_agendaSet"),
                      DESCRIPTION(R"--(Sets *propmat_clearsky_agenda* to a default value

Please consider using *propmat_clearsky_agendaAuto* instead of one of these options
as it will ensure you have the best coverage of use cases.  The options below are
available for feature testing

Options are:

- ``"Empty"``:

    1. Uses *propmat_clearskyInit* to set *propmat_clearsky*, *nlte_source*, *dpropmat_clearsky_dx*, and *dnlte_source_dx*
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("propmat_clearsky_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("refr_index_air_agendaSet"),
                      DESCRIPTION(R"--(Sets *refr_index_air_agenda* to a default value

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
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("refr_index_air_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("rte_background_agendaSet"),
                      DESCRIPTION(R"--(Sets *rte_background_agenda* to a default value

Options are:
    FIXME
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("rte_background_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("sensor_response_agendaSet"),
                      DESCRIPTION(R"--(Sets *sensor_response_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("sensor_response_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("spt_calc_agendaSet"),
                      DESCRIPTION(R"--(Sets *spt_calc_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("spt_calc_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("surface_rtprop_agendaSet"),
                      DESCRIPTION(R"--(Sets *surface_rtprop_agenda* to a default value

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
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("surface_rtprop_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("test_agendaSet"),
                      DESCRIPTION(R"--(Sets *test_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("test_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("water_p_eq_agendaSet"),
                      DESCRIPTION(R"--(Sets *water_p_eq_agenda* to a default value

Options are:

- ``"MK05"``:
    1. Uses *water_p_eq_fieldMK05* to set *water_p_eq_field*
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("water_p_eq_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT("MK05"),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("ybatch_calc_agendaSet"),
                      DESCRIPTION(R"--(Sets *ybatch_calc_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("ybatch_calc_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("option"),
                      GIN_TYPE("String"),
                      GIN_DEFAULT(NODEF),
                      GIN_DESC("Default agenda option (see description)"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  md_data_raw.push_back(
      create_mdrecord(NAME("iy_main_agendaSetByPart"),
                      DESCRIPTION(R"--(Sets *ybatch_calc_agenda* to a default value

Options are:
    There are currently no options, calling this function is an error.
)--"),
                      AUTHORS("Richard Larsson"),
                      OUT("iy_main_agenda"),
                      GOUT(),
                      GOUT_TYPE(),
                      GOUT_DESC(),
                      IN(),
                      GIN("rte_option", "propagation_properties_option", "background_option", "ppath_option"),
                      GIN_TYPE("String","String","String","String"),
                      GIN_DEFAULT(NODEF,NODEF,NODEF,NODEF),
                      GIN_DESC("Choice for RTE calculations", "Choice of propagation properties calculations", "Choice of background radiation calculations", "Choice of propagation path calculations"),
                      SETMETHOD(false),
                      AGENDAMETHOD(false),
                      USES_TEMPLATES(false),
                      PASSWORKSPACE(true)));

  //! Special method that has to look through some of the above methods for changes
  {
    //! Special method that has to look through some of the above methods for changes
    ArrayOfString gin;
    ArrayOfString gintype;
    ArrayOfString gindefault;
    ArrayOfString gindesc;
    const ArrayOfString targets = {
        "propmat_clearskyInit",
        "propmat_clearskyAddCIA",
        "propmat_clearskyAddLines",
        "propmat_clearskyAddZeeman",
        "propmat_clearskyAddFaraday",
        "propmat_clearskyAddXsecFit",
        "propmat_clearskyAddParticles",
        "propmat_clearskyAddFromLookup",
        "propmat_clearskyAddPredefined",
        "propmat_clearskyAddOnTheFlyLineMixing",
        "propmat_clearskyAddHitranLineMixingLines",
        "propmat_clearskyAddOnTheFlyLineMixingWithZeeman",
        };
    Index method_counter = 0;
    for (auto& m : md_data_raw) {
      if (std::find(targets.cbegin(), targets.cend(), m.Name()) not_eq
          targets.cend()) {
        method_counter++;
        for (auto& x : m.GIn()) gin.push_back(x);
        for (auto& x : m.GInType()) {
          gintype.push_back(global_data::wsv_groups.at(x).name);
          gindesc.push_back("See *" + m.Name() + "*");
        }
        for (auto& x : m.GInDefault()) gindefault.push_back(x);
      }
    }
    if (method_counter not_eq targets.nelem()) throw std::logic_error("Lacking functions");
    String doc{R"--(Sets the *propmat_clearsky_agenda* automatically

This method introspects the input and uses it for generating the
*propmat_clearsky_agenda* automatically.  If ``use_abs_lookup``, all
methods that can be used to generate the absorption lookup table
are ignored and instead the calculations from the absorption
lookup are used.

The following methods are considered for addition:
)--"};
    Index count=1;
    for (auto& m: targets) {
        doc += "    " + std::to_string(count++) + ") *" + m + "*\n";
    }

    doc += R"--(
To perform absorption lookupo table calculation, call:
    1) *propmat_clearsky_agendaAuto*
    2) *abs_lookupCalc*
    3) *propmat_clearsky_agendaAuto* (use_abs_lookup=1)
    4) Perform other calculations
)--";

    //Remove duplicate GINs.
    for (Index i=0; i<gin.nelem(); i++) {
      for (Index j=i+1; j<gin.nelem(); j++) {
        if (gin[j] < gin[i]) {
          std::swap(gin[i], gin[j]);
          std::swap(gintype[i], gintype[j]);
          std::swap(gindefault[i], gindefault[j]);
          std::swap(gindesc[i], gindesc[j]);
        }
      }
    }

    auto ptr = std::adjacent_find(gin.begin(), gin.end());
    while (ptr != gin.end()){
      const Index i = std::distance(gin.begin(),ptr);
      ARTS_ASSERT(gintype[i]==gintype[i+1], "Same name variable must be of the same type");
      ARTS_ASSERT(gindefault[i]==gindefault[i+1], "Same name variable must have the same default");
      gindesc[i]+="; "+gindesc[i+1];

      gin.erase(gin.begin()+i+1);
      gintype.erase(gintype.begin()+i+1);
      gindefault.erase(gindefault.begin()+i+1);
      gindesc.erase(gindesc.begin()+i+1);
      ptr = std::adjacent_find(gin.begin(), gin.end());
    }

    gin.push_back("use_abs_lookup");
    gintype.push_back("Index");
    gindefault.push_back("0");
    gindesc.push_back("Uses lookup calculations if true, ignores methods that can be part of the lookup table");

    md_data_raw.push_back(
        MdRecord("propmat_clearsky_agendaAuto",
                 doc.c_str(),
                 {"Richard Larsson"},
                 {"propmat_clearsky_agenda", "propmat_clearsky_agenda_checked"},
                 {},
                 {},
                 {},
                 {"abs_species", "abs_lines_per_species"},
                 gin,
                 gintype,
                 gindefault,
                 gindesc,
                 false,
                 false,
                 false,
                 true,
                 false));
  }

  std::sort(md_data_raw.begin(), md_data_raw.end(), [](auto& a, auto& b) {
    return a.Name() < b.Name();
  });
}
