/* Copyright (C) 2000-2012
   Stefan Buehler <sbuehler@uni-bremen.de>
   Patrick Eriksson <patrick.eriksson@chalmers.se>
   Oliver Lemke <olemke@ltu.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

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
#include "wsv_aux.h"

namespace global_data {
Array<MdRecord> md_data_raw;
extern const ArrayOfString wsv_group_names;
}  // namespace global_data

// Some #defines and typedefs to make the records better readable:
#define NAME(x) x
#define DESCRIPTION(x) x
#define AUTHORS(...) \
  { __VA_ARGS__ }
#define OUT(...) \
  { __VA_ARGS__ }
#define GOUT(...) \
  { __VA_ARGS__ }
#define GOUT_TYPE(...) \
  { __VA_ARGS__ }
#define GOUT_DESC(...) \
  { __VA_ARGS__ }
#define IN(...) \
  { __VA_ARGS__ }
#define GIN(...) \
  { __VA_ARGS__ }
#define GIN_TYPE(...) \
  { __VA_ARGS__ }
#define GIN_DEFAULT(...) \
  { __VA_ARGS__ }
#define GIN_DESC(...) \
  { __VA_ARGS__ }
#define SETMETHOD(x) x
#define AGENDAMETHOD(x) x
#define USES_TEMPLATES(x) x
#define PASSWORKSPACE(x) x
#define PASSWSVNAMES(x) x

/* Here's a template record entry:  (PE 2008-09-20)

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MethodName" ),
        DESCRIPTION
        (
        "A concise summary of the method.\n"
        "\n"
        "A more detailed description of the method. Try to describe the\n"
        "purpose of the method and important considerations. Try to avoid\n"
        "references to other WSMs as they might change. Refer to the user\n"
        "guide for more complex information (as long as it exists, or that\n"
        "you add it to AUG!).\n"
        "\n"
        "You do not need to describe workspace variables used. That\n"
        "information is found in workspace.cc. Generic\n"
        "output and input variables must be described in GIN_DESC and\n"
        "GOUT_DESC below.\n"
        ),
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
    ( MdRecord
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

  using global_data::wsv_group_names;

  for (ArrayOfString::const_iterator it = wsv_group_names.begin();
       it != wsv_group_names.end();
       it++) {
    if (*it != "Any") {
      md_data_raw.push_back(MdRecord(
          NAME(String(*it + "Create").c_str()),
          DESCRIPTION(
              String("Creates a variable of group " + *it +
                     ".\n"
                     "\n"
                     "After being created, the variable is uninitialized.\n")
                  .c_str()),
          AUTHORS("Oliver Lemke"),
          OUT(),
          GOUT("out"),
          GOUT_TYPE((*it).c_str()),
          GOUT_DESC("Variable to create."),
          IN(),
          GIN(),
          GIN_TYPE(),
          GIN_DEFAULT(),
          GIN_DESC()));
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // Let's put in the functions in alphabetical order. This gives a clear rule
  // for where to place a new function and this gives a nicer results when
  // the functions are listed by "arts -m all".
  // No distinction is made between uppercase and lowercase letters.
  // Patrick Eriksson 2002-05-08
  /////////////////////////////////////////////////////////////////////////////

  md_data_raw.push_back(MdRecord(
      NAME("AbsInputFromAtmFields"),
      DESCRIPTION("Initialises the WSVs *abs_p*, *abs_t* and *abs_vmrs* from\n"
                  "*p_grid, *t_field* and *vmr_field*.\n"
                  "\n"
                  "This only works for a 1D atmosphere!\n"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_p", "abs_t", "abs_vmrs"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim", "p_grid", "t_field", "vmr_field"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("AbsInputFromRteScalars"),
      DESCRIPTION(
          "Initialize absorption input WSVs from local atmospheric conditions.\n"
          "\n"
          "The purpose of this method is to allow an explicit line-by-line\n"
          "calculation, e.g., by *abs_coefCalc*, to be put inside the\n"
          "*propmat_clearsky_agenda*. What the method does is to prepare absorption\n"
          "input parameters (pressure, temperature, VMRs), from the input\n"
          "parameters to *propmat_clearsky_agenda*.\n"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_p", "abs_t", "abs_nlte", "abs_vmrs"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("rtp_pressure", "rtp_temperature", "rtp_nlte", "rtp_vmr"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("abs_cia_dataAddCIARecord"),
      DESCRIPTION(
          "Takes CIARecord as input and appends the results in the appropriate place.\n"
          "\n"
          "If CIARecord has same species as species in *abs_cia_data*, then the array\n"
          "position is used to append all of the CIARecord into the array.  If clobber\n"
          "evaluates as true, cia_record overwrites the appropriate *abs_cia_data*.  If\n"
          "species in cia_record are not in *abs_cia_data*, the CIARecord is pushed back.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("abs_cia_dataReadFromCIA"),
      DESCRIPTION(
          "Read data from a CIA data file for all CIA molecules defined\n"
          "in *abs_species*.\n"
          "\n"
          "The units in the HITRAN file are:\n"
          "Frequency: cm^(-1)\n"
          "Binary absorption cross-section: cm^5 molec^(-2)\n"
          "\n"
          "Upon reading we convert this to the ARTS internal SI units \n"
          "of Hz and m^5 molec^(-2).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("abs_cia_dataReadFromXML"),
      DESCRIPTION(
          "Read data from a CIA XML file and check that all CIA tags defined\n"
          "in *abs_species* are present in the file.\n"
          "\n"
          "The units of the data are described in *abs_cia_dataReadFromCIA*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("abs_coefCalcFromXsec"),
      DESCRIPTION("Calculate absorption coefficients from cross sections.\n"
                  "\n"
                  "This calculates both the total absorption and the\n"
                  "absorption per species.\n"
                  "\n"
                  "Cross sections are multiplied by n*VMR.\n"),
      AUTHORS("Stefan Buehler", "Axel von Engeln"),
      OUT("abs_coef",
          "src_coef",
          "dabs_coef_dx",
          "dsrc_coef_dx",
          "abs_coef_per_species",
          "src_coef_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_xsec_per_species",
         "src_xsec_per_species",
         "dabs_xsec_per_species_dx",
         "dsrc_xsec_per_species_dx",
         "abs_species",
         "jacobian_quantities",
         "abs_vmrs",
         "abs_p",
         "abs_t"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("abs_cont_descriptionAppend"),
      DESCRIPTION(
          "Appends the description of a continuum model or a complete absorption\n"
          "model to *abs_cont_names* and *abs_cont_parameters*.\n"
          "\n"
          "See online documentation for *abs_cont_names* for a list of\n"
          "allowed models and for information what parameters they require. See\n"
          "file includes/continua.arts for default parameters for the various models.\n"),
      AUTHORS("Thomas Kuhn", "Stefan Buehler"),
      OUT("abs_cont_names", "abs_cont_models", "abs_cont_parameters"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_cont_names", "abs_cont_models", "abs_cont_parameters"),
      GIN("tagname", "model", "userparam"),
      GIN_TYPE("String", "String", "Vector"),
      GIN_DEFAULT(NODEF, NODEF, "[]"),
      GIN_DESC(
          "The name (species tag) of a continuum model. Must match one\n"
          "of the models implemented in ARTS.\n",
          "A string selecting a particular continuum/full model under this\n"
          "species tag.\n",
          "A Vector containing the required parameters for the selected model.\n"
          "The meaning of the parameters and how many parameters are required\n"
          "depends on the model.\n")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_cont_descriptionInit"),
      DESCRIPTION(
          "Initializes the two workspace variables for the continuum description,\n"
          "*abs_cont_names* and *abs_cont_parameters*.\n"
          "\n"
          "This method does not really do anything, except setting the two\n"
          "variables to empty Arrays. It is just necessary because the method\n"
          "*abs_cont_descriptionAppend* wants to append to the variables.\n"
          "\n"
          "Formally, the continuum description workspace variables are required\n"
          "by the absorption calculation methods (e.g., *abs_coefCalc*). Therefore you\n"
          "always have to call at least *abs_cont_descriptionInit*, even if you do\n"
          "not want to use any continua.\n"),
      AUTHORS("Thomas Kuhn", "Stefan Buehler"),
      OUT("abs_cont_names", "abs_cont_models", "abs_cont_parameters"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lineshapeDefine"),
      DESCRIPTION(
          "Set the lineshape for all calculated lines.\n"
          "\n"
          "Sets the lineshape function. Beside the lineshape function itself, you\n"
          "also have so select a forefactor and a frequency cutoff. The\n"
          "forefactor is later multiplied with the lineshape function.\n"
          "\n"
          "The cutoff frequency is used to make lineshapes finite in frequency,\n"
          "the response outside the cutoff is set to zero, and the lineshape\n"
          "value at the cutoff frequency is subtracted from the overall lineshape\n"
          "as a constant offset. This ensures that the lineshape goes to zero at\n"
          "the cutoff frequency without a discontinuity.\n"
          "\n"
          "We generate only one copy of the lineshape settings. Absorption\n"
          "routines check for this case and use it for all species.\n"
          "\n"
          "The allowed values for the input parameters are:\n"
          "\n"
          "shape:\n"
          "   no_shape:                 no specified shape\n"
          "   Doppler:                  Doppler lineshape\n"
          "   Lorentz:                  Lorentz lineshape\n"
          "   Voigt_Kuntz6:             Kuntz approximation to the Voigt lineshape,\n"
          "                             accuracy > 2x10^(-6)\n"
          "   Voigt_Drayson:            Drayson approximation to the Voigt lineshape\n"
          "   Rosenkranz_Voigt_Drayson: Rosenkrantz oxygen absortion with overlap correction\n"
          "                             on the basis of Drayson routine\n"
          "   Rosenkranz_Voigt_Kuntz6:  Rosenkrantz oxygen absortion with overlap correction\n"
          "                             on the basis of Kuntz routine, accuracy > 2x10^(-6)\n"
          "   CO2_Lorentz:              Lorentz multiplied with Cousin's chi factors\n"
          "   CO2_Drayson:              Drayson multiplied with Cousin's chi factors\n"
          "   Faddeeva_Algorithm_916:   Faddeeva function based on Zaghloul, M.R. and\n"
          "                             A.N. Ali (2011).  Implementation by Steven G. Johnson\n"
          "                             under the MIT License (attainable through\n"
          "                             http://ab-initio.mit.edu/Faddeeva)\n"
          "   Hartmann-Tran:            Line shape is considered as described by the\n"
          "                             Hartmann-Tran profile.\n"
          "\n"
          "forefactor:\n"
          "   no_norm:                  1\n"
          "   Rosenkranz_quadratic:     (f/f0)^2 * h*f/(2k*T)/sinh(h*f/(2k*T))\n"
          "   VVH:                      (f*tanh(h*f/(2k*T))) / (f0*tanh(h*f0/(2k*T)))\n"
          "   VVW:                      (f*f) / (f0*f0)\n"
          "\n"
          "cutoff:\n"
          "    -1:                      no cutoff\n"
          "   <Number>:                 positive cutoff frequency in Hz\n"),
      AUTHORS("Axel von Engeln", "Stefan Buehler"),
      OUT("abs_lineshape"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("shape", "forefactor", "cutoff"),
      GIN_TYPE("String", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Line shape function.",
               "Normalization factor.",
               "Cutoff frequency [Hz].")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lineshape_per_tgDefine"),
      DESCRIPTION(
          "Set the lineshape, separately for each absorption species.\n"
          "\n"
          "This method is similar to *abs_lineshapeDefine*, except that a\n"
          "different lineshape can be set for each absorption species (see\n"
          "*abs_species*). For example, you might want to use different values of\n"
          "the cutoff frequency for different species.\n"
          "\n"
          "For detailed documentation on the available options for the input\n"
          "parameters see documentation of method *abs_lineshapeDefine*.\n"),
      AUTHORS("Axel von Engeln", "Stefan Buehler"),
      OUT("abs_lineshape"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("shape", "normalizationfactor", "cutoff"),
      GIN_TYPE("ArrayOfString", "ArrayOfString", "Vector"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Line shape function for each species.",
               "Normalization factor for each species.",
               "Cutoff frequency [Hz] for each species.")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesFromSplitLines"),
      DESCRIPTION(
          "Combines *abs_lines_per_species* in to *abs_lines* without sorting\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesNewestArtscatFromLegacyCatalog"),
      DESCRIPTION("Sets the version of the catalog to the latest version\n"
                  "Internally, nothing should change.  Writing of catalog\n"
                  "will change\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesReadFromArts"),
      DESCRIPTION("Read all the lines from an Arts catalogue file in the\n"
                  "given frequency range. Otherwise a runtime error will be\n"
                  "thrown\n"
                  "\n"
                  "Please note that all lines must correspond\n"
                  "to legal species / isotopologue combinations\n"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("filename", "fmin", "fmax"),
      GIN_TYPE("String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Name (and path) of the catalogue file.",
               "Minimum frequency for lines to read [Hz].",
               "Maximum frequency for lines to read [Hz].")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesReadFromHitran"),
      DESCRIPTION(
          "Read all the lines from HITRAN 2004 and later catalogue file in\n"
          "the given frequency range. Otherwise a runtime error is thrown.\n"
          "\n"
          "Records of molecules unknown to ARTS are ignored but a\n"
          "warning is issued. In particular this happens for CH3OH\n"
          "(HITRAN molecule number 39) because there is no total internal\n"
          "partition sum available.\n"
          "\n"
          "The database must be sorted by increasing frequency!\n"
          "\n"
          "WWW access of the HITRAN catalogue: http://www.hitran.com/\n"
          "\n"
          "For data in the Hitran 1986-2001 format use the workspace\n"
          "method *abs_linesReadFromHitranPre2004*\n"),
      AUTHORS("Hermann Berg", "Thomas Kuhn"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("filename", "fmin", "fmax"),
      GIN_TYPE("String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Name (and path) of the catalogue file.",
               "Minimum frequency for lines to read [Hz].",
               "Maximum frequency for lines to read [Hz].")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesReadFromHitranPre2004"),
      DESCRIPTION(
          "Read all the lines from a HITRAN 1986-2001 catalogue file in\n"
          "the given frequency range. Otherwise a runtime error will be\n"
          "thrown. For HITRAN 2004 and later line data use the workspace\n"
          "method *abs_linesReadFromHitran*.\n"
          "\n"
          "Please note that all lines must correspond to legal\n"
          "species / isotopologue combinations and that the line data\n"
          "file must be sorted by increasing frequency\n"
          "\n"
          "WWW access of the HITRAN catalogue: http://www.hitran.com/\n"),
      AUTHORS("Thomas Kuhn"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("filename", "fmin", "fmax"),
      GIN_TYPE("String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Name (and path) of the catalogue file.",
               "Minimum frequency for lines to read [Hz].",
               "Maximum frequency for lines to read [Hz].")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesReadFromLBLRTM"),
      DESCRIPTION(
          "Read all the lines from a LBLRTM catalogue file in\n"
          "the given frequency range. Otherwise a runtime error will be\n"
          "thrown.  \n"
          "\n"
          "Note that this method will automatically set line mixing format\n"
          "to LBLRTM for the lines affected, but the abs_species tag must\n"
          "be set to activate line mixing for the species.\n"
          "\n"
          "Access data via: http://rtweb.aer.com/\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("filename", "fmin", "fmax"),
      GIN_TYPE("String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Name (and path) of the catalogue file.",
               "Minimum frequency for lines to read [Hz].",
               "Maximum frequency for lines to read [Hz].")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesReadFromJpl"),
      DESCRIPTION(
          "Read all the lines from a JPL catalogue file in the\n"
          "given frequency range. Otherwise a runtime error will be\n"
          "thrown\n"
          "\n"
          "Please note that all lines must correspond\n"
          "to legal species / isotopologue combinations.\n"
          "\n"
          "WWW access of the JPL catalogue: http://spec.jpl.nasa.gov/\n"),
      AUTHORS("Thomas Kuhn"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("filename", "fmin", "fmax"),
      GIN_TYPE("String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Name (and path) of the catalogue file.",
               "Minimum frequency for lines to read [Hz].",
               "Maximum frequency for lines to read [Hz].")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesReadFromMytran2"),
      DESCRIPTION("Read all the lines from a MYTRAN2 catalogue file in the\n"
                  "given frequency range. Otherwise a runtime error will be\n"
                  "thrown\n"
                  "\n"
                  "Please note that all lines must correspond\n"
                  "to legal species / isotopologue combinations\n"),
      AUTHORS("Axel von Engeln", "Stefan Buehler"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("filename", "fmin", "fmax"),
      GIN_TYPE("String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Name (and path) of the catalogue file.",
               "Minimum frequency for lines to read [Hz].",
               "Maximum frequency for lines to read [Hz].")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesReadFromSplitArtscat"),
      DESCRIPTION(
          "Read all the lines in the given frequency range from a split\n"
          "Arts catalogue file.\n"
          "\n"
          "Please note that all lines must correspond\n"
          "to legal species / isotopologue combinations\n"),
      AUTHORS("Oliver Lemke"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("basename", "fmin", "fmax"),
      GIN_TYPE("String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Basename of the catalogue.",
               "Minimum frequency for lines to read [Hz].",
               "Maximum frequency for lines to read [Hz].")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesRelativeLineStrengthShift"),
      DESCRIPTION(
          "Shift line strength of all lines in *abs_lines*.\n"
          "\n"
          "The new line strengths are the old line strengths times (1 +\n"
          "relative_line_strength_shift).\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("relative_line_strength_shift"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("0"),
      GIN_DESC("Relative change in line strengths.")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesReplaceWithLines"),
      DESCRIPTION(
          "Replace all lines in *abs_lines* that match with lines in replacement_lines.\n"
          "\n"
          "Each replacement_lines must match at most a single line in *abs_lines*.\n"
          "\n"
          "The matching required identical quantum number signatures to work\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("replacement_lines"),
      GIN_TYPE("ArrayOfLineRecord"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Line-array that replace lines in *abs_lines*.")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesReplaceParameterWithLinesParameter"),
      DESCRIPTION(
          "Replace parameter of all lines in *abs_lines* that match with lines in replacement_lines.\n"
          "Only works for:\n"
          "parameter_name = \"Central Frequency\"\n"
          "parameter_name = \"Line Strength\"\n"
          "parameter_name = \"Line Shape Model\"\n"
          "parameter_name = \"Lower State Energy\"\n"
          "\n"
          "Each replacement_lines must match at most a single line in *abs_lines*.\n"
          "\n"
          "The matching required identical quantum number signatures to work\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("replacement_lines", "parameter_name"),
      GIN_TYPE("ArrayOfLineRecord", "String"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Line-array that replace lines in *abs_lines*.",
               "Name of parameter to be replaced")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lines_per_bandSetLineMixingFromRelmat"),
      DESCRIPTION("A dummy method to test line mixing.\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_band"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_band", "relmat_per_band", "partition_functions"),
      GIN("temperatures", "linemixing_type", "do_g", "do_dv"),
      GIN_TYPE("Vector", "String", "Index", "Index"),
      GIN_DEFAULT(NODEF, "LM2", "1", "1", ),
      GIN_DESC(
          "Vector of temperatures to compute the relaxation matrix at",
          "String describing type of line mixing adaptation in linerecord after computations",
          "Index to indicate if g is to be left as zero",
          "Index to indicate if dv is to be left as zero")));

  md_data_raw.push_back(
      MdRecord(NAME("abs_linesSetNormalizationForAll"),
               DESCRIPTION("Sets normalization type for all lines.\n"
                           "\n"
                           "Available options:\n"
                           "   VVH\n"
                           "   VVW\n"
                           "   RosenkranzQuadratic\n"
                           "   None\n"),
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
      MdRecord(NAME("abs_lines_per_speciesSetNormalizationForAll"),
               DESCRIPTION("See *abs_linesSetNormalizationForAll*\n"),
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

  md_data_raw.push_back(MdRecord(NAME("abs_linesSetNlteOffForAll"),
                                 DESCRIPTION("Sets no NLTE for all lines\n"),
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
      MdRecord(NAME("abs_lines_per_speciesSetNlteOffForAll"),
               DESCRIPTION("See *abs_linesSetNlteOffForAll*\n"),
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
      MdRecord(NAME("abs_linesSetMirroringForAll"),
               DESCRIPTION("Sets mirroring type for all lines.\n"
                           "\n"
                           "Available options:\n"
                           "   Lorentz\n"
                           "   SameAsLineShape\n"
                           "   None\n"),
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
      MdRecord(NAME("abs_lines_per_speciesSetMirroringForAll"),
               DESCRIPTION("See *abs_linesSetMirroringForAll*\n"),
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
      MdRecord(NAME("abs_linesCutOffForAll"),
               DESCRIPTION("Sets cutoff frequency for all lines.\n"
                           "\n"
                           "Available options:\n"
                           "   above 0: cutoff frequency\n"
                           "   -1: no cutoff frequency\n"),
               AUTHORS("Richard Larsson"),
               OUT("abs_lines"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines"),
               GIN("option"),
               GIN_TYPE("Numeric"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Method of line mirroring")));

  md_data_raw.push_back(MdRecord(NAME("abs_lines_per_speciesCutOffForAll"),
                                 DESCRIPTION("See *abs_linesCutOffForAll*\n"),
                                 AUTHORS("Richard Larsson"),
                                 OUT("abs_lines_per_species"),
                                 GOUT(),
                                 GOUT_TYPE(),
                                 GOUT_DESC(),
                                 IN("abs_lines_per_species"),
                                 GIN("option"),
                                 GIN_TYPE("Numeric"),
                                 GIN_DEFAULT(NODEF),
                                 GIN_DESC("Method of line mirroring")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesChangeBaseParameterForMatchingLines"),
      DESCRIPTION(
          "Change parameter of all lines in *abs_lines* that match with *QuantumIdentifier*.\n"
          "Only works for these parameters:\n"
          "parameter_name = \"Central Frequency\"\n"
          "parameter_name = \"Line Strength\"\n"
          "parameter_name = \"Lower State Energy\"\n"
          "\n"
          "Note that loose_matching:=0 means only identical quantum identifiers are accepted,\n"
          "otherwise the numbers in QI must just be contained in the line identifier\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("QI", "parameter_name", "change", "relative", "loose_matching"),
      GIN_TYPE("QuantumIdentifier", "String", "Numeric", "Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "0", "0"),
      GIN_DESC("Information to match the line.",
               "Name of parameter to be replaced",
               "Value with which to change matching line{'s,s'}",
               "Flag for relative change (0 is absolute change)",
               "Flag for loose match (0 means only complete matches)")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesSetBaseParameterForMatchingLines"),
      DESCRIPTION(
          "Same as *abs_linesChangeBaseParameterForMatchingLines* but sets the parameter\n"
          "instead of changing it.  See the other function for inputs.\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("QI", "parameter_name", "new_value", "loose_matching"),
      GIN_TYPE("QuantumIdentifier", "String", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "0"),
      GIN_DESC("Information to match the line.",
               "Name of parameter to be replaced",
               "New value of parameter for matching line(s)",
               "Flag for loose match (0 means only complete matches)")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesSetLineShapeModelParameterForMatchingLines"),
      DESCRIPTION(
        "Same as *abs_linesChangeLineShapeModelParameterForMatchingLines* but\n"
        "sets the line paramater to the provided new_value\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("QI", "parameter", "coefficient", "species", "new_value"),
      GIN_TYPE("QuantumIdentifier", "String", "String", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, "0"),
      GIN_DESC("Information to match the line.",
               "Name of parameter to be replaced",
               "Coefficient of the parameter to be changed",
               "Species of parameter to be changed",
               "New value of parameter for matching line(s)")));

  md_data_raw.push_back(MdRecord(
    NAME("abs_linesChangeLineShapeModelParameterForMatchingLines"),
      DESCRIPTION("Change line shape model data parameter in matching lines.\n"
        "\n"
        "The matching is done so that QI must be in the line identifier\n"
        "\n"
        "Acceptable parameter(s) are:\n"
        "\t\"G0\"\n"
        "\t\"D0\"\n"
        "\t\"G2\"\n"
        "\t\"D2\"\n"
        "\t\"FVC\"\n"
        "\t\"ETA\"\n"
        "\t\"Y\"\n"
        "\t\"G\"\n"
        "\t\"DV\"\n"
        "\n"
        "Acceptable coefficient(s) are:\n"
        "\t\"X0\"\n"
        "\t\"X1\"\n"
        "\t\"X2\"\n"
        "\n"
        "Acceptable species are:\n"
        "\tAIR (so long as it is the broadening species list)\n"
        "\tSELF (so long as it is the broadening species list)\n"
        "\tAny species in the line broadening species\n"
        "\n"
        "The line parameter will have its old value plus the change if\n"
        "relative is false, else it will have its old value times\n"
        "(1+change).\n"
        "\n"
        "Throws an error if it cannot find any targets to change\n"
      ),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("QI", "parameter", "coefficient", "species", "change", "relative"),
      GIN_TYPE("QuantumIdentifier",
               "String",
               "String",
               "String",
               "Numeric",
               "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, NODEF, "0"),
      GIN_DESC("Information to match the line.",
               "Name of parameter to be replaced",
               "Coefficient of the parameter to be changed",
               "Species of parameter to be changed",
               "Change in the value found",
               "Flag for relative change (0 is absolute change)")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesSetQuantumNumberForAll"),
      DESCRIPTION(
          "Sets a quantum number of all lines in abs_lines\n"
          "Useful when a catalog misses a number that is required in ARTS\n"
          "\n"
          "Warning: the only test performed inside is if the name of the\n"
          "provided quantum number is valid.  It is up to the user to provide\n"
          "a valid quantum number.\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("where", "quantum_number_name", "quantum_number_value"),
      GIN_TYPE("Index", "String", "Rational"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC(
          "Set to zero for both, negative for only lower, or positive for only upper levels.",
          "Name of the quantum number",
          "Value of the quantum number")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_linesShiftFrequency"),
      DESCRIPTION(
          "Shift line center of all lines in *abs_lines*.\n"
          "\n"
          "The new frequencies are the old frequencies plus the input frequency shift.\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("frequency_shift"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("0"),
      GIN_DESC("Frequency to shift line centers [Hz].")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lines_per_bandFromband_identifiers"),
      DESCRIPTION("Sets abs_lines_per_band and related variables from\n"
                  "*band_identifiers*.\n"
                  "\n"
                  "If *band_identifiers* identifies bands that are not in\n"
                  "the line data, this method remains silent.  This allows\n"
                  "the user to use autmated identification methods such as:\n"
                  "*SetBandIdentifiersAuto* to identify more bands than will\n"
                  "be used in subsequent calculations.\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_band",
          "abs_species_per_band",
          "abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species", "abs_species", "band_identifiers"),
      GIN("change_linemixing_to_bandwise"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC(
          "(Boolean) Set linemixing to need *abs_xsec_per_speciesAddLineMixedBands* for these lines?")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lines_per_speciesAddMirrorLines"),
      DESCRIPTION(
          "Adds mirror lines at negative frequencies to *abs_lines_per_species*.\n"
          "\n"
          "For each line at frequency +f in *abs_lines_per_species* a corresponding\n"
          "entry at frequency -f is added to *abs_lines_per_species*. The mirror\n"
          "lines are appended to the line list after the original lines.\n"),
      AUTHORS("Axel von Engeln", "Stefan Buehler", "Patrick Eriksson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN("max_f"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("-1"),
      GIN_DESC("Limit for mirroring, ie. lines above this frequency do "
               "not generate a mirror line. All lines mirrored if *max_f* "
               "is < 0, that is the default setting.")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lines_per_speciesCompact"),
      DESCRIPTION(
          "Removes all lines outside the defined lineshape cutoff frequencies\n"
          "from *abs_lines_per_species*. This can save computation time.\n"
          "It should be particularly useful to call this method after\n"
          "*abs_lines_per_speciesAddMirrorLines*.\n"),
      AUTHORS("Axel von Engeln", "Stefan Buehler"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species", "abs_lineshape", "f_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lines_per_speciesCreateFromLines"),
      DESCRIPTION(
          "Split lines up into the different species.\n"
          "\n"
          "The species are tested in the order in which they are specified in the\n"
          "controlfile. Lines are assigned to the first species that\n"
          "matches. That means if the list of species is [\"O3-666\",\"O3\"], then\n"
          "the last group O3 gets assigned all the O3 lines that do not fit in\n"
          "the first group (all other isotopologues than the main isotopologue).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("abs_lines_per_speciesReadFromCatalogues"),
      DESCRIPTION(
          "Read spectral line data from different line catalogues.\n"
          "\n"
          "For each absorption species, you can specify which catalogue to\n"
          "use. Because the method creates *abs_lines_per_species* directly, it\n"
          "replaces for example the following two method calls:\n"
          "\n"
          "  - abs_linesReadFromHitran\n"
          "  - abs_lines_per_speciesCreateFromLines\n"
          "\n"
          "This method needs as input WSVs the list of species\n"
          "*abs_species*. Generic input parameters must specify the names of the\n"
          "catalogue files to use and the matching formats.  Names can be\n"
          "anything, formats can currently be HITRAN96 (for HITRAN 1986-2001\n"
          "databases), HITRAN04 (for HITRAN 2004 database), MYTRAN2, JPL, or\n"
          "ARTS.  Furthermore, you have to specify minimum and maximum frequency\n"
          "for each species. To safe typing, if there are less elements in the\n"
          "keyword parameters than there are species, the last parameters are\n"
          "applied to all following species.\n"
          "\n"
          "Example usage:\n"
          "\n"
          "abs_lines_per_speciesReadFromCatalogues(\n"
          "   [ \"../data/cat1.dat\", \"../data/cat2.dat\" ]\n"
          "   [ \"MYTRAN2\",          \"HITRAN96\"         ]\n"
          "   [ 0,                  0                  ]\n"
          "   [ 2000e9,             100e9              ]\n"
          ")\n"
          "\n"
          "In this example, lines for the first species will be taken from cat1,\n"
          "lines for all other species will be taken from cat2. This allows you\n"
          "for example to use a special line file just for water vapor lines.\n"
          "\n"
          "Catalogues are only read once, even if several tag groups have the\n"
          "same catalogue. However, in that case the frequency ranges MUST be the\n"
          "same. (If you want to do fine-tuning of the frequency ranges, you can\n"
          "do this inside the tag definitions, e.g., \"H2O-*-0-2000e9\".)\n"
          "\n"
          "This function uses the various reading routines\n"
          "(*abs_linesReadFromHitran*, etc.), as well as\n"
          "*abs_lines_per_speciesCreateFromLines*.\n"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("filenames", "formats", "fmin", "fmax"),
      GIN_TYPE("ArrayOfString", "ArrayOfString", "Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF),
      GIN_DESC("Name (and path) of the catalogue files.",
               "Format of each file. (Allowed formats are\n"
               "HITRAN96, HITRAN04, MYTRAN2, JPL, ARTS.",
               "Minimum frequency for lines to read [Hz].",
               "Maximum frequency for lines to read [Hz].")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lines_per_speciesRelativeLineStrengthShift"),
      DESCRIPTION(
        "As *abs_linesRelativeLineStrengthShift* but for all lines in *abs_lines_per_species*\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN("relative_line_strength_shift"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("0"),
      GIN_DESC("Relative change in line strengths.")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lines_per_speciesSetEmpty"),
      DESCRIPTION(
          "Sets abs_lines_per_species to empty line lists.\n"
          "\n"
          "You can use this method to set *abs_lines_per_species* if you do not\n"
          "really want to compute line spectra. Formally, abs_coefCalc will still\n"
          "require *abs_lines_per_species* to be set.\n"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lines_per_speciesShiftFrequency"),
      DESCRIPTION(
        "As *abs_linesShiftFrequency* but for all lines in *abs_lines_per_species*\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species"),
      GIN("frequency_shift"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("0"),
      GIN_DESC("Frequency to shift line centers [Hz].")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lines_per_speciesWriteToSplitArtscat"),
      DESCRIPTION("Write each species to a separate catalogue file.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("output_file_format", "abs_lines_per_species"),
      GIN("basename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(""),
      GIN_DESC("Basename of the catalogue.")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lookupAdapt"),
      DESCRIPTION(
          "Adapts a gas absorption lookup table to the current calculation.\n"
          "\n"
          "The lookup table can contain more species and more frequencies than\n"
          "are needed for the current calculation. This method cuts down the\n"
          "table in memory, so that it contains just what is needed. Also, the\n"
          "species in the table are brought in the same order as the species in\n"
          "the current calculation.\n"
          "\n"
          "Of course, the method also performs quite a lot of checks on the\n"
          "table. If something is not ok, a runtime error is thrown.\n"
          "\n"
          "The method sets a flag *abs_lookup_is_adapted* to indicate that the\n"
          "table has been checked and that it is ok. Never set this by hand,\n"
          "always use this method to set it!\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("abs_lookupCalc"),
      DESCRIPTION(
          "Creates a gas absorption lookup table.\n"
          "\n"
          "The lookup table stores absorption cross-sections as a function of\n"
          "pressure. Additionally, absorption can be stored as a function of\n"
          "temperature for temperature perturbations from a reference\n"
          "profile.\n"
          "\n"
          "Additionally, absorption can be stored as a function of water vapor\n"
          "VMR perturbations from a reference profile. The variable *abs_nls*\n"
          "specifies, for which species water vapor perturbations should be\n"
          "generated.\n"
          "\n"
          "Note, that the absorbing gas can be any gas, but the perturbing gas is\n"
          "always H2O.\n"
          "\n"
          "In contrast to other absorption functions, this method does not use\n"
          "the input variable *abs_h2o*. This is because *abs_h2o* has to be set\n"
          "interally to allow perturbations. If there are more than one H2O\n"
          "species, the first is assumed to be the main one.\n"),
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
         "abs_xsec_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lookupInit"),
      DESCRIPTION(
          "Creates an empty gas absorption lookup table.\n"
          "\n"
          "This is mainly there to help developers. For example, you can write\n"
          "the empty table to an XML file, to see the file format.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("abs_lookupSetup"),
      DESCRIPTION(
          "Set up input parameters for abs_lookupCalc.\n"
          "\n"
          "More information can be found in the documentation for method\n"
          "*abs_lookupSetupBatch*\n"
          "\n"
          "Max and min values of H2O and temperature are adjusted to allow for\n"
          "numerical perturbations in Jacobian calculation.\n"
          "\n"
          "The input variables *abs_nls_interp_order* and *abs_t_interp_order*\n"
          "are used to make sure that there are enough points in *abs_nls_pert*\n"
          "and *abs_t_pert* for the chosen interpolation order.\n"
          "\n"
          "Note: For homogeneous 1D cases, it can be advantageous to calculate\n"
          "*abs_lookup* from the 1D atmosphere, and to expand the atmosphere\n"
          "to 3D only after that. This particularly if nonlinear species\n"
          "(i.e., H2O) are involved."
          "\n"
          "See also:\n"
          "   *abs_lookupSetupBatch*\n"),
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
      IN("atmosphere_dim",
         "p_grid",
         //            "lat_grid",
         //            "lon_grid",
         "t_field",
         "vmr_field",
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

  md_data_raw.push_back(MdRecord(
      NAME("abs_lookupSetupBatch"),
      DESCRIPTION(
          "Set up input parameters for abs_lookupCalc for batch calculations.\n"
          "\n"
          "This method performs a similar task as *abs_lookupSetup*, with the\n"
          "difference, that the lookup table setup is not for a single\n"
          "atmospheric state, but for a whole batch of them, stored in\n"
          "*batch_atm_fields_compact*.\n"
          "\n"
          "The method checks *abs_species* to decide, which species depend on\n"
          "*abs_h2o*, and hence require nonlinear treatment in the lookup table.\n"
          "\n"
          "The method also checks which range of pressures, temperatures, and\n"
          "VMRs occurs, and sets *abs_p*, *abs_t*, *abs_t_pert*, and *abs_vmrs*\n"
          "accordingly.\n"
          "\n"
          "If nonlinear species are present, *abs_nls* and *abs_nls_pert* are also\n"
          "generated.\n"
          "\n"
          "Max and min values of H2O and temperature are adjusted to allow for\n"
          "numerical perturbations in Jacobian calculation.\n"
          "\n"
          "The input variables *abs_nls_interp_order* and *abs_t_interp_order*\n"
          "are used to make sure that there are enough points in *abs_nls_pert*\n"
          "and *abs_t_pert* for the chosen interpolation order.\n"
          "\n"
          "The method checks each given field using *atmfields_checkedCalc*.\n"
          "If a field does not pass the check, a run-time error is thrown.\n"
          "To prevent this, the parameter *robust* can be set to one: Invalid \n"
          "atmospheres are skipped, but the run continues. This matches the \n"
          "robust behaviour of *ybatchCalc*.\n"
          "\n"
          "See also:\n"
          "   *abs_lookupSetup*\n"),
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
         "abs_nls_interp_order",
         "atmosphere_dim"),
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
               "A flag with value 1 or 0. If set to one, the batch\n"
               "setup will continue, even if individual fields are invalid.\n"
               "This is consistent with the behaviour of *ybatchCalc*.",
               /* check_gridnames */
               "A flag with value 1 or 0. If set to one, the gridnames of \n"
               " every *atm_fields_compact* are checked.")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lookupSetupWide"),
      DESCRIPTION(
          "Set up input parameters for abs_lookupCalc for a wide range of\n"
          "atmospheric conditions.\n"
          "\n"
          "This method can be used to set up parameters for a lookup table that\n"
          "really covers all reasonable atmospheric conditions.\n"
          "\n"
          "Reference profiles of T and H2O will be constant, so that the\n"
          "different dimensions in the lookup table are actually \"orthogonal\",\n"
          "unlike the traditional case where we have pressure dependent reference\n"
          "profiles. This makes the table numerically somewhat more robust then\n"
          "the traditional ones, and it makes it straightforward to calculate the\n"
          "accuracy for the different interpolations with abs_lookupTestAccuracy.\n"
          "\n"
          "You can give min an max values for the atmospheric conditions. The\n"
          "default values are chosen such that they cover the value range over\n"
          "the complete Chevallier91L data set, and a bit more. The statistics\n"
          "of the Chevallier91L data are:\n"
          "\n"
          "min(p)   / max(p)   [Pa]:  1 / 104960\n"
          "min(T)   / max(T)   [K]:   158.21 / 320.39\n"
          "min(H2O) / max(H2O) [VMR]: -5.52e-07 / 0.049\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("abs_lookupTestAccuracy"),
      DESCRIPTION(
          "Test accuracy of absorption lookup table.\n"
          "\n"
          "Explicitly compare absorption from the lookup table with line-by-line\n"
          "calculations for strategically selected conditions (in-between the\n"
          "lookup table grid points).\n"
          "\n"
          "For error units see *abs_lookupTestAccMC*\n"
          "\n"
          "Produces no workspace output, only output to the output streams.\n"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lookup",
         "abs_lookup_is_adapted",
         "abs_p_interp_order",
         "abs_t_interp_order",
         "abs_nls_interp_order",
         "abs_xsec_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("abs_lookupTestAccMC"),
      DESCRIPTION(
          "Test accuracy of absorption lookup table with Monte Carlo Algorithm.\n"
          "\n"
          "Explicitly compare absorption from the lookup table with line-by-line\n"
          "calculations for random conditions.\n"
          "\n"
          "The quantities returned are the mean value and standard deviation of\n"
          "the absolute value of the relative error in percent.\n"
          "The relative error itself is computed for a large number of cases\n"
          "(pressure, temperature, and H2O VMR combinations). In the frequency\n"
          "dimension the maximum value is taken for each case.\n"
          "\n"
          "Produces no workspace output, only output to the output streams.\n"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lookup",
         "abs_lookup_is_adapted",
         "abs_p_interp_order",
         "abs_t_interp_order",
         "abs_nls_interp_order",
         "mc_seed",
         "abs_xsec_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("abs_speciesAdd"),
      DESCRIPTION(
          "Adds species tag groups to the list of absorption species.\n"
          "\n"
          "This WSM is similar to *abs_speciesSet*, the only difference is that\n"
          "this method appends species to an existing list of absorption species instead\n"
          "of creating the whole list.\n"
          "\n"
          "See *abs_speciesSet* for details on how tags are defined and examples of\n"
          "how to input them in the control file.\n"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_species",
          "propmat_clearsky_agenda_checked",
          "abs_xsec_agenda_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("species"),
      GIN_TYPE("ArrayOfString"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Specify one String for each tag group that you want to\n"
               "add. Inside the String, separate the tags by commas\n"
               "(plus optional blanks).\n")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_speciesAdd2"),
      DESCRIPTION(
          "Adds a species tag group to the list of absorption species and\n"
          "jacobian quantities.\n"
          "\n"
          "The method is basically a combined call of *abs_speciesAdd* and\n"
          "*jacobianAddAbsSpecies*. In this way it is not needed to specify a\n"
          "tag group in two different places.\n"
          "\n"
          "Arguments exactly as for *jacobianAddAbsSpecies*. Note that this\n"
          "method only handles a single tag group, in contrast to\n"
          "*abs_speciesAdd*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("abs_species",
          "jacobian_quantities",
          "jacobian_agenda",
          "propmat_clearsky_agenda_checked",
          "abs_xsec_agenda_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species", "atmosphere_dim", "p_grid", "lat_grid", "lon_grid"),
      GIN("gin1", "gin2", "gin3", "species", "unit"),
      GIN_TYPE("Vector",
               "Vector",
               "Vector",
               "String",
               "String"),
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

  md_data_raw.push_back(MdRecord(
      NAME("abs_speciesDefineAllInScenario"),
      DESCRIPTION(
          "Define one tag group for each species known to ARTS and included in an\n"
          "atmospheric scenario.\n"
          "\n"
          "You can use this as an alternative to *abs_speciesSet* if you want to make an\n"
          "absorption calculation that is as complete as possible. The method\n"
          "goes through all defined species and tries to open the VMR file. If\n"
          "this works the tag is included, otherwise it is skipped.\n"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_species",
          "propmat_clearsky_agenda_checked",
          "abs_xsec_agenda_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("basename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC(
          "The name and path of a particular atmospheric scenario.\n"
          "For example: /pool/lookup2/arts-data/atmosphere/fascod/tropical")));

  md_data_raw.push_back(
      MdRecord(NAME("abs_speciesInit"),
               DESCRIPTION("Sets  *abs_species* to be empty.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("abs_speciesSet"),
      DESCRIPTION(
          "Set up a list of absorption species tag groups.\n"
          "\n"
          "Workspace variables like *abs_species* contain several tag\n"
          "groups. Each tag group contains one or more tags. This method converts\n"
          "descriptions of tag groups given in the keyword to the ARTS internal\n"
          "representation (an *ArrayOfArrayOfSpeciesTag*). A tag group selects\n"
          "spectral features which belong to the same species.\n"
          "\n"
          "A tag is defined in terms of the name of the species, isotopologue, and a\n"
          "range of frequencies. Species are named after the standard chemical\n"
          "names, e.g., \"O3\". Isotopologues are given by the last digit of the atomic\n"
          "weight, i.g., \"O3-668\" for the asymmetric ozone molecule including an\n"
          "oxygen 18 atom. Groups of transitions are specified by giving a lower\n"
          "and upper limit of a frequency range, e.g., \"O3-666-500e9-501e9\".\n"
          "\n"
          "To turn on Zeeman calculation for a species, \"-Z\" may be appended\n"
          "to its name: \"O2-Z\" or \"O2-Z-66\"\n"
          "\n"
          "To turn on line mixing calculation for a species, \"-LM\" may be appended\n"
          "to its name (or after the Zeeman tag): \"O2-LM\" or \"O2-Z-LM-66\"\n"
          "\n"
          "The symbol \"*\" acts as a wild card. Furthermore, frequency range or\n"
          "frequency range and isotopologue may be omitted.\n"
          "\n"
          "Finally, instead of the isotopologue the special letter \"nl\" may be given,\n"
          "e.g., \"H2O-nl\". This means that no absorption at all is associated\n"
          "with this tag. (It is not quite clear if this feature is useful for\n"
          "anything right now.)\n"
          "\n"
          "Example:\n"
          "\n"
          "   species = [ \"O3-666-500e9-501e9, O3-686\",\n"
          "               \"O3\",\n"
          "               \"H2O-PWR98\" ]\n"
          "\n"
          "   The first tag group selects all O3-666 lines between 500 and\n"
          "   501 GHz plus all O3-686 lines. \n"
          "\n"
          "   The second tag group selects all remaining O3 transitions.\n"
          "\n"
          "   The third tag group selects H2O, with one of the complete\n"
          "   absorption models (Rosenkranz 98). No spectrocopic line catalogue\n"
          "   data will be used for that third tag group.\n"
          "\n"
          "   Note that order of tag groups in the species list matters. In our\n"
          "   example, changing the order of the first two tag group will give\n"
          "   different results: as \"O3\" already selects all O3 transitions,\n"
          "   no lines will remain to be selected by the\n"
          "   \"O3-666-500e9-501e9, O3-686\" tag.\n"
          "\n"
          "For CIA species the tag consists of the two involved species and\n"
          "a dataset index. CIA species can be defined for multiple regions\n"
          "The dataset index determines which region to use from the corresponding\n"
          "CIARecord in *abs_cia_data*.\n"
          "\n"
          "Example\n"
          "\n"
          "species = [ \"N2-CIA-N2-0, N2-CIA-N2-1\" ]\n"
          "\n"
          "For Hitran cross section species the tag consists of the species and\n"
          "the tagtype HXSEC, e.g. CFC11-HXSEC. The data for the species must be\n"
          "available in the *hitran_xsec_data* variable."
          "\n"
          "*abs_xsec_agenda_checked* and *propmat_clearsky_agenda_checked*\n"
          "are set to be false.\n"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_species",
          "abs_xsec_agenda_checked",
          "propmat_clearsky_agenda_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("species"),
      GIN_TYPE("ArrayOfString"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Specify one String for each tag group that you want to\n"
               "create. Inside the String, separate the tags by commas\n"
               "(plus optional blanks).\n")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_vecAddGas"),
      DESCRIPTION(
          "Add gas absorption to first element of absorption vector.\n"
          "\n"
          "The task of this method is to sum up the gas absorption of the\n"
          "different gas species and add the result to the first element of the\n"
          "absorption vector.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("abs_xsec_agenda_checkedCalc"),
      DESCRIPTION(
          "Checks if the *abs_xsec_agenda* contains all necessary\n"
          "methods to calculate all the species in *abs_species*.\n"
          "\n"
          "This method should be called just before the *abs_xsec_agenda*\n"
          "is used, e.g. *abs_lookupCalc*, *ybatchCalc*, *yCalc*\n"),
      AUTHORS("Oliver Lemke"),
      OUT("abs_xsec_agenda_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species", "abs_xsec_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("abs_xsec_per_speciesAddCIA"),
      DESCRIPTION(
          "Calculate absorption cross sections per tag group for HITRAN CIA continua.\n"
          "\n"
          "This interpolates the cross sections from *abs_cia_data*.\n"
          "\n"
          "The robust option is intended only for testing. Do not use for normal\n"
          "runs, since subsequent functions will not be able to deal with NAN values.\n"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_xsec_per_species", "dabs_xsec_per_species_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_xsec_per_species",
         "dabs_xsec_per_species_dx",
         "abs_species",
         "jacobian_quantities",
         "abs_species_active",
         "f_grid",
         "abs_p",
         "abs_t",
         "abs_vmrs",
         "abs_cia_data"),
      GIN("T_extrapolfac", "robust"),
      GIN_TYPE("Numeric", "Index"),
      GIN_DEFAULT("0.5", "0"),
      GIN_DESC(
          "Temperature extrapolation factor (relative to grid spacing).",
          "Set to 1 to suppress runtime errors (and return NAN values instead).")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_xsec_per_speciesAddHitranXsec"),
      DESCRIPTION(
          "Calculate absorption cross sections per tag group for HITRAN xsec species.\n"
          "\n"
          "This broadens the cross section data from *hitran_xsec_data* and\n"
          "interpolates it onto the current f_grid.\n"
          "\n"
          "apply_tfit turns of the temperature fit. It is only meant for testing\n"
          "and should alwasy be kept on for real calculations.\n"
          "\n"
          "This method depends on the FFTW-3 library.\n"),
      AUTHORS("Oliver Lemke"),
      OUT("abs_xsec_per_species", "dabs_xsec_per_species_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_xsec_per_species",
         "dabs_xsec_per_species_dx",
         "abs_species",
         "jacobian_quantities",
         "abs_species_active",
         "f_grid",
         "abs_p",
         "abs_t",
         "hitran_xsec_data"),
      GIN("apply_tfit"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Apply temperature fit.")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_xsec_per_speciesAddConts"),
      DESCRIPTION(
          "Calculate absorption cross sections per tag group for continua.\n"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_xsec_per_species", "dabs_xsec_per_species_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_xsec_per_species",
         "dabs_xsec_per_species_dx",
         "abs_species",
         "jacobian_quantities",
         "abs_species_active",
         "f_grid",
         "abs_p",
         "abs_t",
         "abs_vmrs",
         "abs_cont_names",
         "abs_cont_parameters",
         "abs_cont_models"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("abs_xsec_per_speciesAddLines"),
      DESCRIPTION(
          "Calculates the line spectrum for both attenuation and phase\n"
          "for each tag group and adds it to abs_xsec_per_species.\n"),
      AUTHORS("Stefan Buehler", "Axel von Engeln"),
      OUT("abs_xsec_per_species",
          "src_xsec_per_species",
          "dabs_xsec_per_species_dx",
          "dsrc_xsec_per_species_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_xsec_per_species",
         "src_xsec_per_species",
         "dabs_xsec_per_species_dx",
         "dsrc_xsec_per_species_dx",
         "abs_species",
         "jacobian_quantities",
         "abs_species_active",
         "f_grid",
         "abs_p",
         "abs_t",
         "abs_nlte",
         "abs_vmrs",
         "abs_lines_per_species",
         "abs_lineshape",
         "isotopologue_ratios",
         "partition_functions"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("abs_xsec_per_speciesAddLines2"),
      DESCRIPTION(
          "Calculates the line spectrum for both attenuation and phase\n"
          "for each tag group and adds it to abs_xsec_per_species.\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_xsec_per_species",
          "src_xsec_per_species",
          "dabs_xsec_per_species_dx",
          "dsrc_xsec_per_species_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_xsec_per_species",
         "src_xsec_per_species",
         "dabs_xsec_per_species_dx",
         "dsrc_xsec_per_species_dx",
         "abs_species",
         "jacobian_quantities",
         "abs_species_active",
         "f_grid",
         "abs_p",
         "abs_t",
         "abs_nlte",
         "abs_vmrs",
         "abs_lines_per_species",
         "isotopologue_ratios",
         "partition_functions"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("abs_xsec_per_speciesAddLineMixedLines"),
      DESCRIPTION("Calculates the band-wise cross-section TEST FUNCTION\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_xsec_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_xsec_per_species",
         "abs_species_per_band",
         "abs_species",
         "f_grid",
         "abs_p",
         "abs_t",
         "relmat_per_band",
         "abs_lines_per_band",
         "isotopologue_ratios",
         "partition_functions"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("abs_xsec_per_speciesAddLineMixedLinesInAir"),
      DESCRIPTION("Calculates the band-wise cross-section TEST FUNCTION\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_xsec_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_xsec_per_species",
         "abs_species_per_band",
         "abs_species",
         "f_grid",
         "abs_p",
         "abs_t",
         "abs_lines_per_band",
         "isotopologue_ratios",
         "partition_functions",
         "wigner_initialized"),
      GIN("minimum_line_count"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("10"),
      GIN_DESC("If less than this number of lines in a \"band\", "
               "relaxation matrix is set diagonal")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_xsec_per_speciesAddLineMixedBands"),
      DESCRIPTION(
          "Calculate absorption cross sections per band for Relmat line mixing.\n"
          "\n"
          "This requires setting *abs_lines_per_band* and *abs_species_per_band*.\n"
          "\n"
          "Will call the Relmat fortran library for the bulk of the calculations.\n"
          "There are to date two libraries available to compute the relaxation\n"
          "matrix.  Please set your method using *relmat_type_per_band*.\n"
          "Helpful method: *SetRelaxationMatrixCalcType*.\n"
          "\n"
          "Calculations of *dabs_xsec_per_species_dx* is done through perturbations\n"
          "for temperature and will therefore be twice as slow.\n"
          "\n"
          "Control mechanism for output happens through various GIN:"
          "   pressure_rule_limit:   Will be compared to |PWij/(fj-fi)|. If any\n"
          "                          value is larger than this value, Wij becomes\n"
          "                          diagonal.  Useful for ordered approximations\n"
          "                          of line mixing because these can easily produce\n"
          "                          absurd results if this test fails.  Set to a\n"
          "                          large value when using full inversions and a\n"
          "                          small value when using perturbation inversion\n"
          "\n"
          "   write_relmat_per_band: Will use order_of_linemixing input to fill\n"
          "                          *relmat_per_band* with relevant line mixing\n"
          "                          information\n"
          "\n"
          "   debug:                 Writes Relmat debug information to screen and\n"
          "                          throws errors if negative.  Otherwise, errors\n"
          "                          are instead printed by out3-level of verbosity\n"
          "\n"
          "   order_of_linemixing:   Sets order of line mixing.  0 means no line\n"
          "                          mixing at all and should reproduce similar\n"
          "                          results as one of the pure LBL methods. 1 or 2\n"
          "                          means that the first or second order line\n"
          "                          mixing coefficients are derived and then calls\n"
          "                          the Voigt line shape to compute cross-section.\n"
          "                          Negative values means that the full inversion\n"
          "                          is used for all frequencies.  This results in\n"
          "                          slow computations that should, in theory, be\n"
          "                          more accurate, though this remains to be\n"
          "                          fully tested\n"
          "\n"
          "   use_adiabatic_factor:  If true, the adiabatic factor is used. This\n"
          "                          means you have to have the intermolecular\n"
          "                          distances defined in the code to use the line\n"
          "                          mixing computations\n"
          "\n"
          "Note that for pressures below *lm_p_lim* there are no calls to Relmat\n"
          "if order_of_linemixing is positive.  However, if *lm_p_lim* is negative,\n"
          "Relmat is called with abs(order_of_linemixing).  This allows using\n"
          "full inversions of the relaxation matrix at higher pressures but also use\n"
          "reduced inversions at lower pressures to reproduce Doppler broadening\n"
          "effects.  This means negative order_of_linemixing should never be absolutely\n"
          "larger than the highest supported positive order_of_linemixing.\n"
          "\n"
          "If there is an error using Relmat, please turn on the debug GIN for output\n"
          "to screen and rerun the code to identify the problem.  Relmat makes many\n"
          "assumptions on the user input.  It will sometimes return a diagonal matrix\n"
          "when the theoretical limits are breached by the user input. Please see\n"
          "published literature for examples of such theoretical breaches.\n"
          "\n"
          "More information later.\n"),
      AUTHORS("Teresa Mendaza", "Richard Larsson"),
      OUT("abs_xsec_per_species",
          "dabs_xsec_per_species_dx",
          "relmat_per_band"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_xsec_per_species",
         "dabs_xsec_per_species_dx",
         "abs_lines_per_band",
         "abs_species_per_band",
         "band_identifiers",
         "abs_species",
         "isotopologue_ratios",
         "partition_functions",
         "jacobian_quantities",
         "f_grid",
         "abs_p",
         "abs_t",
         "lm_p_lim",
         "relmat_type_per_band",
         "wigner_initialized"),
      GIN("pressure_rule_limit",
          "write_relmat_per_band",
          "debug",
          "order_of_linemixing",
          "use_adiabatic_factor"),
      GIN_TYPE("Numeric", "Index", "Index", "Index", "Index"),
      GIN_DEFAULT("1e100", "0", "0", "-1", "1"),
      GIN_DESC(
          "Limit when perturbation theory is assumed to work",
          "Writes the relaxation operators instead of the cross-section if true",
          "Lets relmat know it is to print debug information if true.",
          "Choice of order of linemixing",
          "Truth-value if we should use the precomputed adiabatic factors")));

  md_data_raw.push_back(MdRecord(
      NAME("abs_xsec_per_speciesInit"),
      DESCRIPTION(
          "Initialize *abs_xsec_per_species*.\n"
          "\n"
          "The initialization is\n"
          "necessary, because methods *abs_xsec_per_speciesAddLines*\n"
          "and *abs_xsec_per_speciesAddConts* just add to *abs_xsec_per_species*.\n"
          "The size is determined from *abs_species*.\n"),
      AUTHORS("Stefan Buehler"),
      OUT("abs_xsec_per_species",
          "src_xsec_per_species",
          "dabs_xsec_per_species_dx",
          "dsrc_xsec_per_species_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species",
         "jacobian_quantities",
         "abs_species_active",
         "f_grid",
         "abs_p",
         "abs_xsec_agenda_checked",
         "nlte_do"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("AddZaAa"),
      DESCRIPTION(
          "Adds zenith and azimuth angles.\n"
          "\n"
          "Adds up line-of-sights (LOS). In short, *dlos* is added to *ref_los*,\n"
          "assuming that a unit changes in zenith and azimuth are equal where\n"
          "dlos=(0,0).\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("new_los"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("End line-of-sights."),
      IN(),
      GIN("ref_los", "dlos"),
      GIN_TYPE("Vector", "Matrix"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Reference line-of-sight (a single LOS).",
               "Change in line-of-sight (can be multiple LOS).")));

  md_data_raw.push_back(MdRecord(
      NAME("AgendaAppend"),
      DESCRIPTION(
          "Append methods to an agenda.\n"
          "\n"
          "An agenda is used to store a list of methods that are meant to be\n"
          "executed sequentially.\n"
          "\n"
          "This method takes the methods given in the body (in the curly braces)\n"
          "and appends them to the agenda given by the output argument (in the round\n"
          "braces).\n"
          "\n"
          "It also uses the agenda lookup data (defined in file agendas.cc) to\n"
          "check, whether the given methods use the right input WSVs and produce\n"
          "the right output WSVs.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Agenda"),
      GOUT_DESC("Target agenda."),
      IN(),
      GIN("in"),
      GIN_TYPE("Agenda"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Source agenda."),
      SETMETHOD(false),
      AGENDAMETHOD(true),
      USES_TEMPLATES(false),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(MdRecord(NAME("AgendaExecute"),
                                 DESCRIPTION("Execute an agenda.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("AgendaExecuteExclusive"),
      DESCRIPTION(
          "Execute an agenda exclusively.\n"
          "\n"
          "Only one call to *AgendaExecuteExclusive* is executed at a time.\n"
          "Other calls to this function are blocked until the current one\n"
          "finishes. WARNING: Can cause deadlocks! Use with care.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("AgendaSet"),
      DESCRIPTION(
          "Set up an agenda.\n"
          "\n"
          "An agenda is used to store a list of methods that are meant to be\n"
          "executed sequentially.\n"
          "\n"
          "This method takes the methods given in the body (in the curly braces)\n"
          "and puts them in the agenda given by the output argument (in the round\n"
          "braces).\n"
          "\n"
          "It also uses the agenda lookup data (defined in file agendas.cc) to\n"
          "check, whether the given methods use the right input WSVs and\n"
          "produce the right output WSVs.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
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

  md_data_raw.push_back(MdRecord(
      NAME("AngularGridsSetFluxCalc"),
      DESCRIPTION(
          "Sets the angular grids for the calculation of radiation fluxes\n"
          "(irradiance) per hemispheres and heating rates\n"
          "\n"
          "This method sets the angular grids for the radiation fluxes type\n"
          "calculations and calculates the integration weights *za_grid_weights*\n"
          "for the zenith angle integration. For down- und up-looking\n"
          "geometries it suffices to define *N_za_grid* and\n"
          "*N_aa_grid*. From *N_aa_grid* an equally spaced grid is\n"
          "created and stored in the WSV*scat_aa_grid*.\n"
          "Depending on the desired *za_grid_type* *scat_za_grid* will be\n"
          "equally spaced ('linear') or unequally ('linear_mu','double_gauss')\n"
          "Important, *N_za_grid* must be an even number because for the \n"
          "integration over each hemisphere *N_za_grid* / 2 zenith angles are needed.\n"
          "\n"
          "Possible zenith angle grid types are:\n"
          "double_gauss:     The zenith grid and the integration weights are set according\n"
          "                  to a gauss-legendre integration for each hemispheres.\n"
          "linear:           Equally space grid between 0 deg and 180 deg including the poles\n"
          "linear_mu:        Similar to 'linear' but equally spaced for cos(180 deg) to cos(0 deg),\n"
          "                  which results a unequally spaced angular grid\n"

          ),
      AUTHORS("Manfred Brath"),
      OUT("scat_za_grid", "scat_aa_grid", "za_grid_weights"),
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

  md_data_raw.push_back(MdRecord(
      NAME("ArrayOfAgendaAppend"),
      DESCRIPTION("Set up an agenda and append it to the array of agendas.\n"
                  "\n"
                  "See *AgendaSet* for details.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
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
      MdRecord(NAME("ArrayOfAgendaExecute"),
               DESCRIPTION("Execute an agenda from an ArrayOfAgenda.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("AntennaConstantGaussian1D"),
      DESCRIPTION(
          "Sets up a 1D gaussian antenna response and a matching\n"
          "*mblock_dlos_grid*.\n"
          "\n"
          "As *antenna_responseGaussian*, but also creates *mblock_dlos_grid*.\n"
          "For returned antenna response, see *antenna_responseGaussian*.\n"
          "\n"
          "The size of *mblock_dlos_grid* is determined by *n_za_grid*.\n"
          "The end points of the grid are set to be the same as for the\n"
          "antenna response. The spacing of the grid follows the magnitude of\n"
          "the response; the spacing is smaller where the response is high.\n"
          "More precisely, the grid points are determined by dividing\n"
          "the cumulative sum of the response in equal steps. This makes sense\n"
          "if the representation error of the radiance (as a function of\n"
          "zenith angle) increases linearly with the grid spacing.\n"
          "\n"
          "The WSV *antenna_dlos* is set to [0].\n"
          "\n"
          "The antenna repsonse is not normalised.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("antenna_dim",
          "mblock_dlos_grid",
          "antenna_response",
          "antenna_dlos"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("n_za_grid", "fwhm", "xwidth_si", "dx_si"),
      GIN_TYPE("Index", "Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, "3", "0.1"),
      GIN_DESC("Number of points to include in *mblock_dlos_grid*.",
               "Full width at half-maximum of antenna beam [deg].",
               "Half-width of response, in terms of std. dev.",
               "Grid spacing, in terms of std. dev.")));

  md_data_raw.push_back(MdRecord(
      NAME("AntennaMultiBeamsToPencilBeams"),
      DESCRIPTION(
          "Maps a multi-beam case to a matching pencil beam case.\n"
          "\n"
          "Cases with overlapping beams are most efficiently handled by\n"
          "letting *antenna_dlos* have several rows. That is, there are\n"
          "multiple beams for each measurement block. The drawback is that\n"
          "many variables must be adjusted if the corresponding pencil beam\n"
          "spectra shall be calculated. This method makes this adjustment.\n"
          "That is, if you have a control file for a multiple beam case and\n"
          "for some reason want to avoid the antenna weighting, you add this\n"
          "method before *sensor_responseInit*, and remove the call of\n"
          "*sensor_responseAntenna* and you will get the matching pencil beam\n"
          "spectra.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_pos",
          "sensor_los",
          "antenna_dlos",
          "antenna_dim",
          "mblock_dlos_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_pos",
         "sensor_los",
         "antenna_dlos",
         "antenna_dim",
         "mblock_dlos_grid",
         "atmosphere_dim"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("AntennaOff"),
      DESCRIPTION(
          "Sets some antenna related variables\n"
          "\n"
          "Use this method to set *antenna_dim* and *mblock_dlos_grid* to\n"
          "suitable values (1 and [0], respectively) for cases when a\n"
          "sensor is included, but the antenna pattern is neglected.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("antenna_dim", "mblock_dlos_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("antenna_responseGaussian"),
      DESCRIPTION(
          "Sets up a gaussian antenna response.\n"
          "\n"
          "The method assumes that the response is the same for all\n"
          "frequencies and polarisations, and that it can be modelled as\n"
          "gaussian.\n"
          "\n"
          "The grid generated is approximately\n"
          "   si * [-xwidth_si:dx_si:xwidth_si]\n"
          "where si is the standard deviation corresponding to the FWHM.\n"
          "That is, width and spacing of the grid is specified in terms of\n"
          "number of standard deviations. If xwidth_si is set to 2, the\n"
          "response will cover about 95% the complete response. For\n"
          "xwidth_si=3, about 99% is covered. If xwidth_si/dx_si is not\n"
          "an integer, the end points of the grid are kept and the spacing\n"
          "of the grid is reduced (ie. spacing is equal or smaller *dx_si*).\n"
          "\n"
          "If the 2D option is selected (*do_2d*), a circular antenna is\n"
          "assumed and the response is any direction follows the 1D case.\n"
          "\n"
          "The antenna repsonse is not normalised.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("antenna_response"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("fwhm", "xwidth_si", "dx_si", "do_2d"),
      GIN_TYPE("Numeric", "Numeric", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, "3", "0.1", "0"),
      GIN_DESC("Full width at half-maximum",
               "Half-width of response, in terms of std. dev.",
               "Grid spacing, in terms of std. dev.",
               "Set to 1 to create a 2D antenna pattern.")));

  md_data_raw.push_back(MdRecord(
      NAME("antenna_responseVaryingGaussian"),
      DESCRIPTION(
          "Sets up gaussian antenna responses.\n"
          "\n"
          "Similar to *antenna_responseGaussian* but allows to set up\n"
          "responses that varies with frequency. That is, the method assumes\n"
          "that the response is the same for all polarisations, and that it\n"
          "can be modelled as a gaussian function varying with frequency.\n"
          "\n"
          "The full width at half maximum (FWHM in radians) is calculated as:\n"
          "    fwhm = lambda / leff\n"
          "where lambda is the wavelength and *leff* is the effective size of\n"
          "the antenna. Normally, *leff* is smaller than the physical antenna\n"
          "size.\n"
          "\n"
          "Antenna responses are created for *nf* frequencies spanning the\n"
          "range [*fstart*,*fstop*], with a logarithmic spacing. That is, the\n"
          "frequency grid of the responses is taken from *VectorNLogSpace*.\n"
          "\n"
          "The responses have a common angular grid. The width, determined by\n"
          "*xwidth_si*, is set for the lowest frequency, while the spacing\n"
          "(*dx_si*) is set for the highest frequency. This ensures that both\n"
          "the width and spacing are equal or better than *xwidth_si* and\n"
          "*dx_si*, respectively, for all frequencies.\n"
          "\n"
          "If the 2D option is selected (*do_2d*), a circular antenna is\n"
          "assumed and the response is any direction follows the 1D case.\n"
          "\n"
          "The antenna repsonse is not normalised.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("antenna_response"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("leff", "xwidth_si", "dx_si", "nf", "fstart", "fstop", "do_2d"),
      GIN_TYPE("Numeric",
               "Numeric",
               "Numeric",
               "Index",
               "Numeric",
               "Numeric",
               "Index"),
      GIN_DEFAULT(NODEF, "3", "0.1", NODEF, NODEF, NODEF, "0"),
      GIN_DESC("Effective size of the antenna",
               "Half-width of response, in terms of std. dev.",
               "Grid spacing, in terms of std. dev.",
               "Number of points in frequency grid (must be >= 2)",
               "Start point of frequency grid",
               "End point of frequency grid",
               "Set to 1 to create a 2D antenna pattern.")));

  md_data_raw.push_back(MdRecord(
      NAME("Append"),
      DESCRIPTION(
          "Append one workspace variable to another.\n"
          "\n"
          "This method can append an array to another array of the same type,\n"
          "e.g. ArrayOfIndex to ArrayOfIndex. Or a single element to an array\n"
          "such as a Tensor3 to an ArrayOfTensor3.\n"
          "\n"
          "Appending two vectors or a numeric to a vector works as for array\n"
          "variables.\n"
          "\n"
          "Both another matrix or a vector can be appended to a matrix. In\n"
          "addition, for matrices, the 'append dimension' can be selected.\n"
          "The third argument, *dimension*, indicates how to append, where\n"
          "\"leading\" means to append row-wise, and \"trailing\" means\n"
          "column-wise.\n"
          "\n"
          "Other types (TensorX) are currently only implemented for\n"
          "appending to the leading dimension.\n"
          "\n"
          "This method is not implemented for all types, just for those that\n"
          "were thought or found to be useful. (See variable list below.).\n"),
      AUTHORS("Stefan Buehler, Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector, Vector,"
                "Matrix, Matrix,"
                "Tensor3, Tensor3,"
                "Tensor4, Tensor4,"
                "String, " +
                ARRAY_GROUPS + ", " + ARRAY_GROUPS_WITH_BASETYPE),
      GOUT_DESC("The variable to append to."),
      IN(),
      GIN("in", "dimension"),
      GIN_TYPE("Numeric, Vector,"
               "Matrix, Vector,"
               "Matrix, Tensor3,"
               "Tensor3, Tensor4,"
               "String, " +
                   ARRAY_GROUPS + "," + GROUPS_WITH_ARRAY_TYPE,
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

  md_data_raw.push_back(MdRecord(
      NAME("ArrayOfGriddedFieldGetNames"),
      DESCRIPTION("Get the names of all GriddedFields stored in an Array.\n"
                  "\n"
                  "See *GriddedFieldGetName*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("ArrayOfIndexLinSpace"),
      DESCRIPTION(
          "Initializes an ArrayOfIndex with linear spacing.\n"
          "\n"
          "The first element equals always the start value, and the spacing\n"
          "equals always the step value, but the last value can deviate from\n"
          "the stop value. *step* can be both positive and negative.\n"
          "\n"
          "The created array is [start, start+step, start+2*step, ...]\n "),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("ArrayOfIndex"),
      GOUT_DESC("Output array."),
      IN(),
      GIN("start", "stop", "step"),
      GIN_TYPE("Index", "Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Start value.",
               "Maximum/minimum value of the end value",
               "Spacing of the array.")));

  md_data_raw.push_back(MdRecord(
      NAME("ArrayOfIndexSet"),
      DESCRIPTION("Creates an ArrayOfIndex from the given list of numbers.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("ArrayOfIndex"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("value"),
      GIN_TYPE("ArrayOfIndex"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Indexes for initializiation."),
      SETMETHOD(true)));

  md_data_raw.push_back(MdRecord(
      NAME("ArrayOfIndexSetConstant"),
      DESCRIPTION("Creates an ArrayOfIndex of length *nelem*, with all values\n"
                  "identical.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("ArrayOfIndex"),
      GOUT_DESC("Variable to initialize."),
      IN("nelem"),
      GIN("value"),
      GIN_TYPE("Index"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Array value.."),
      SETMETHOD(true)));

  md_data_raw.push_back(MdRecord(
      NAME("ArrayOfStringSet"),
      DESCRIPTION("Sets a String array according the given text.\n"
                  "The format is text = [\"String1\",\"String2\",...]\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("ArrayOfString"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("value"),
      GIN_TYPE("ArrayOfString"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Strings for initialization."),
      SETMETHOD(true)));

  md_data_raw.push_back(MdRecord(
      NAME("Arts"),
      DESCRIPTION(
          "Runs the agenda that is specified inside the curly braces. ARTS\n"
          "controlfiles must define this method. It is executed automatically\n"
          "when ARTS is run on the controlfile and cannot be called by the user.\n"
          "This methods was used for Arts 1 controlfiles and is now obsolete.\n"
          "See *Arts2*\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("Arts2"),
      DESCRIPTION(
          "Runs the agenda that is specified inside the curly braces. ARTS\n"
          "controlfiles must define this method. It is executed automatically\n"
          "when ARTS is run on the controlfile and cannot be called by the user.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("AtmFieldPerturb"),
      DESCRIPTION(
          "Adds a perturbation to an atmospheric field.\n"
          "\n"
          "The shape and position of the perturbation follow the retrieval grids.\n"
          "That is, the shape of the perturbation has a traingular shape, \n"
          "with breake points at the retrieval grid points. The position is\n"
          "given as an index. This index matches the column in the Jacobian\n"
          "for the selected grid position.\n"
          "\n"
          "If the retrieval grids fully match the atmospheric grids, you can\n"
          "use *AtmFieldPerturbAtmGrids*, that is faster. The description of\n"
          "that method can help to understand this method.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("perturbed_field"),
      GOUT_TYPE("Tensor3"),
      GOUT_DESC("Perturbed/modified field."),
      IN("atmosphere_dim","p_grid","lat_grid","lon_grid"),
      GIN("original_field","p_ret_grid","lat_ret_grid","lon_ret_grid",
          "pert_index","pert_size","pert_mode"),
      GIN_TYPE("Tensor3","Vector","Vector","Vector",
               "Index","Numeric","String"),
      GIN_DEFAULT(NODEF,NODEF,NODEF,NODEF,NODEF,NODEF,"absolute"),
      GIN_DESC("Original field, e.g. *t_field*.",
               "Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retrieval grid.",
               "Index of position where the perturbation shall be performed.",
               "Size of perturbation.",
               "Type of perturbation, ""ansolute"" or ""relative"".")));

  md_data_raw.push_back(MdRecord(
      NAME("AtmFieldPerturbAtmGrids"),
      DESCRIPTION(
          "As *AtmFieldPerturb*, but perturbation follows the atmospheric grids.\n"
          "\n"
          "The method effectively performs this\n"
          "  perturbed_field = original_field\n"
          "  perturbed_field(p_index,lat_index,lon_index) += pert_size\n"
          "if not *pert_mode* is set to ""relative"" when this is done\n"
          "  perturbed_field = original_field\n"
          "  perturbed_field(p_index,lat_index,lon_index) *= 1*pert_size\n"
          "where p_index etc. are derived from *pert_index*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("perturbed_field"),
      GOUT_TYPE("Tensor3"),
      GOUT_DESC("Perturbed/modified field."),
      IN("atmosphere_dim","p_grid","lat_grid","lon_grid"),
      GIN("original_field","pert_index","pert_size","pert_mode"),
      GIN_TYPE("Tensor3","Index","Numeric","String"),
      GIN_DEFAULT(NODEF,NODEF,NODEF,"absolute"),
      GIN_DESC("Original field, e.g. *t_field*.",
               "Index of position where the perturbation shall be performed.",
               "Size of perturbation.",
               "Type of perturbation, ""ansolute"" or ""relative"".")));

    md_data_raw.push_back(MdRecord(
      NAME("AtmFieldPRegrid"),
      DESCRIPTION(
          "Interpolates the input field along the pressure dimension from\n"
          "*p_grid_old* to to *p_grid_new*.\n"
          "\n"
          "Extrapolation is allowed within the common 0.5grid-step margin.\n"
          "in and out fields can be the same variable.\n"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Tensor3, Tensor4"),
      GOUT_DESC("Regridded atmospheric field."),
      IN(),
      GIN("in", "p_grid_new", "p_grid_old", "interp_order"),
      GIN_TYPE("Tensor3, Tensor4", "Vector", "Vector", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "1"),
      GIN_DESC("Input atmospheric field.",
               "Pressure grid to regrid to",
               "Pressure grid of input field",
               "Interpolation order.")));

  md_data_raw.push_back(MdRecord(
      NAME("AtmFieldsCalc"),
      DESCRIPTION(
          "Interpolation of raw atmospheric T, z, VMR, and NLTE T/r fields to\n"
          "calculation grids.\n"
          "\n"
          "An atmospheric scenario includes the following data for each\n"
          "position (pressure, latitude, longitude) in the atmosphere:\n"
          "   1. temperature field\n"
          "   2. the corresponding altitude field\n"
          "   3. vmr fields for the gaseous species\n"
          "This method interpolates the fields of raw data (*t_field_raw*,\n"
          "*z_field_raw*, *vmr_field_raw*) which can be stored on arbitrary\n"
          "grids to the calculation grids (*p_grid*, *lat_grid*, *lon_grid*).\n"
          "If *nlte_field_raw* is empty, it is assumed to be so because LTE is\n"
          "assumed by the user and *nlte_field* will be empty.\n"
          "\n"
          "Internally, *AtmFieldsCalc* applies *GriddedFieldPRegrid* and\n"
          "*GriddedFieldLatLonRegrid*. Generally, 'half-grid-step' extrapolation\n"
          "is allowed and applied. However, if *vmr_zeropadding*=1 then VMRs at\n"
          "*p_grid* levels exceeding the raw VMRs' pressure grid are set to 0\n"
          "(applying the *vmr_zeropadding* option of *GriddedFieldPRegrid*).\n"
          "\n"
          "Default is to just accept obtained VMRs. If you want to enforce\n"
          "that all VMR created are >= 0, set *vmr_nonegative* to 1. Negative\n"
          "values are then set 0. Beside being present in input data, negative\n"
          "VMR can be generated from the interpolation if *interp_order* is\n"
          "above 1.\n"),
      AUTHORS("Claudia Emde", "Stefan Buehler"),
      OUT("t_field", "z_field", "vmr_field", "nlte_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("p_grid",
         "lat_grid",
         "lon_grid",
         "t_field_raw",
         "z_field_raw",
         "vmr_field_raw",
         "nlte_field_raw",
         "atmosphere_dim"),
      GIN("interp_order",
          "vmr_zeropadding",
          "vmr_nonegative",
          "nlte_when_negative"),
      GIN_TYPE("Index", "Index", "Index", "Index"),
      GIN_DEFAULT("1", "0", "0", "0"),
      GIN_DESC("Interpolation order (1=linear interpolation).",
               "Pad VMRs with zeroes to fit the pressure grid if necessary.",
               "If set to 1, negative VMRs are set to 0.",
               "-1: Skip step. 0: Negative is 0. Else: Negative is t.")));

  md_data_raw.push_back(MdRecord(
      NAME("AtmFieldsCalcExpand1D"),
      DESCRIPTION(
          "Interpolation of 1D raw atmospheric fields to create 2D or 3D\n"
          "homogeneous atmospheric fields.\n"
          "\n"
          "The method works as *AtmFieldsCalc*, but accepts only raw 1D\n"
          "atmospheres. The raw atmosphere is interpolated to *p_grid* and\n"
          "the obtained values are applied for all latitudes, and also\n"
          "longitudes for 3D, to create a homogeneous atmosphere.\n"
          "\n"
          "The method deals only with the atmospheric fields, and to create\n"
          "a true 2D or 3D version of a 1D case, a demand is also that the\n"
          "ellipsoid is set to be a sphere.\n"),
      AUTHORS("Patrick Eriksson", "Claudia Emde", "Stefan Buehler"),
      OUT("t_field", "z_field", "vmr_field", "nlte_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("p_grid",
         "lat_grid",
         "lon_grid",
         "t_field_raw",
         "z_field_raw",
         "vmr_field_raw",
         "nlte_field_raw",
         "atmosphere_dim"),
      GIN("interp_order",
          "vmr_zeropadding",
          "vmr_nonegative",
          "nlte_when_negative"),
      GIN_TYPE("Index", "Index", "Index", "Index"),
      GIN_DEFAULT("1", "0", "0", "0"),
      GIN_DESC("Interpolation order (1=linear interpolation).",
               "Pad VMRs with zeroes to fit the pressure grid if necessary.",
               "If set to 1, negative VMRs are set to 0.",
               "-1: Skip step. 0: Negative is 0. Else: Negative is t.")));

  md_data_raw.push_back(MdRecord(
      NAME("AtmFieldsExpand1D"),
      DESCRIPTION(
          "Maps a 1D case to 2D or 3D homogeneous atmospheric fields.\n"
          "\n"
          "This method takes a 1D atmospheric case and converts it to the\n"
          "corresponding case for 2D or 3D. The atmospheric fields (t_field,\n"
          "z_field and vmr_field) must be 1D and match *p_grid*. The size of\n"
          "the new data is determined by *atmosphere_dim*, *lat_grid* and\n"
          "*lon_grid*. That is, these later variables have been changed since\n"
          "the original fields were created.\n"
          "\n"
          "The method deals only with the atmospheric fields, and to create\n"
          "a true 2D or 3D version of a 1D case, a demand is also that the\n"
          "ellipsoid is set to be a sphere.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("t_field", "z_field", "vmr_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("t_field",
         "z_field",
         "vmr_field",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "atmosphere_dim"),
      GIN("chk_vmr_nan"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC(
          "Flag to determine if a search for NaN shall be performed or not.")));

  md_data_raw.push_back(MdRecord(
      NAME("AtmFieldsExtract1D"),
      DESCRIPTION(
          "Converts 2D or 3D homogeneous atmospheric fields to a 1D case.\n"
          "\n"
          "The method extracts data for given latitude and longitude index\n"
          "to create a 1D atmosphere. *AtmosphereSet1D* is called to set\n"
          "output values of *atmosphere_dim*, *lat_grid* and *lon_grid*.\n"
          "Nothing is done if *atmosphere_dim* alöready is 1.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("atmosphere_dim",
          "lat_grid",
          "lon_grid",
          "t_field",
          "z_field",
          "vmr_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "lat_grid",
         "lon_grid",
         "t_field",
         "z_field",
         "vmr_field"),
      GIN("ilat", "ilon"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT("0", "0"),
      GIN_DESC("Pick data having this latitude index (0-based).",
               "Pick data having this longitude index (0-based).")));

  md_data_raw.push_back(MdRecord(
      NAME("AtmFieldsRefinePgrid"),
      DESCRIPTION(
          "Refines the pressure grid and regrids the clearsky atmospheric\n"
          "fields accordingly.\n"
          "\n"
          "This method is, e.g., used for absorption lookup table testing. It\n"
          "can also be used to refine the *p_grid* and atmospheric fields from\n"
          "compact state atmospheres.\n"
          "\n"
          "It adds additional vertical grid points to the atmospheric fields, by\n"
          "interpolating them in the usual ARTS way (linear in log pressure).\n"
          "\n"
          "How fine the new grid will be is determined by the keyword parameter\n"
          "p_step. The definition of p_step, and the default interpolation\n"
          "behavior, is consistent with *abs_lookupSetup* and\n"
          "*abs_lookupSetupBatch* (new points are added between the original\n"
          "ones, so that the spacing is always below p_step.)\n"
          "\n"
          "Internally, *AtmFieldsRefinePgrid* applies *p_gridRefine* and\n"
          "*AtmFieldPRegrid* to the clearsky atmospheric fields (T, z, vmr).\n"
          "\n"
          "Atmospheric field related check WSV are reset to 0 (unchecked),\n"
          "i.e., the corresponding checkedCalc methods have to be performed\n"
          "(again) before *yCalc* or similar methods can be executed.\n"),
      AUTHORS("Stefan Buehler"),
      OUT("p_grid",
          "t_field",
          "z_field",
          "vmr_field",
          "atmfields_checked",
          "atmgeom_checked",
          "cloudbox_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("p_grid",
         "lat_grid",
         "lon_grid",
         "t_field",
         "z_field",
         "vmr_field",
         "atmosphere_dim"),
      GIN("p_step", "interp_order"),
      GIN_TYPE("Numeric", "Index"),
      GIN_DEFAULT(NODEF, "1"),
      GIN_DESC("Maximum step in log(p[Pa]) (natural logarithm, as always). If\n"
               "the pressure grid is coarser than this, additional points\n"
               "are added until each log step is smaller than this.\n",
               "Interpolation order.")));

  md_data_raw.push_back(MdRecord(
      NAME("AtmFieldsFromCompact"),
      DESCRIPTION(
          "Extract pressure grid and atmospheric fields from\n"
          "*atm_fields_compact*.\n"
          "\n"
          "An atmospheric scenario includes the following data for each\n"
          "position (pressure, latitude, longitude) in the atmosphere:\n"
          "           1. temperature field\n"
          "           2. the corresponding altitude field\n"
          "           3. vmr fields for the gaseous species\n"
          "           4. scattering species fields\n"
          "\n"
          "This method splits up the data found in *atm_fields_compact* to\n"
          "p_grid, lat_grid, lon_grid, and the various fields.\n"
          "See documentation of *atm_fields_compact* for a definition of the\n"
          "data.\n"
          "\n"
          "Compact states are characterized by having all atmospheric fields\n"
          "already given on identical grids. That is, no interpolation needs\n"
          "to be and is performed. Keyword *p_min* allows to remove atmospheric\n"
          "levels with pressures lower than the given value (default: no\n"
          "removal). This reduces computational burden and is useful when\n"
          "upper atmospheric contributions are negligible.\n"
          "\n"
          "Possible future extensions: Add a keyword parameter to refine the\n"
          "pressure grid if it is too coarse. Or a version that interpolates\n"
          "onto given grids, instead of using and returning the original grids.\n"),
      AUTHORS("Jana Mendrok"),
      OUT("p_grid",
          "lat_grid",
          "lon_grid",
          "t_field",
          "z_field",
          "vmr_field",
          "scat_species_mass_density_field",
          "scat_species_mass_flux_field",
          "scat_species_number_density_field",
          "scat_species_mean_mass_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species", "scat_species", "atm_fields_compact", "atmosphere_dim"),
      GIN("delim", "p_min", "check_gridnames"),
      GIN_TYPE("String", "Numeric", "Index"),
      GIN_DEFAULT("-", "0", "0"),
      GIN_DESC(/* delim */
               "Delimiter string of *scat_species* elements.",
               /* p_min */
               "Minimum-pressure level to consider (for TOA).",
               /* check_gridnames */
               "A flag with value 1 or 0. If set to one, the gridnames of \n"
               " the *atm_fields_compact* are checked.")));

  md_data_raw.push_back(MdRecord(
      NAME("AtmFieldsAndParticleBulkPropFieldFromCompact"),
      DESCRIPTION(
          "Extract pressure grid and atmospheric fields from\n"
          "*atm_fields_compact*.\n"
          "\n"
          "An atmospheric scenario includes the following data for each\n"
          "position (pressure, latitude, longitude) in the atmosphere:\n"
          "           1. temperature field\n"
          "           2. the corresponding altitude field\n"
          "           3. vmr fields for the gaseous species\n"
          "           4. scattering species fields\n"
          "\n"
          "This method splits up the data found in *atm_fields_compact* to\n"
          "p_grid, lat_grid, lon_grid, vmr_field, particle_bulkprop_field,\n"
          "and particle_bulkprop_names.\n"
          "See documentation of *atm_fields_compact* for a definition of the\n"
          "data.\n"
          "\n"
          "Compact states are characterized by having all atmospheric fields\n"
          "already given on identical grids. That is, no interpolation needs\n"
          "to be and is performed. Keyword *p_min* allows to remove atmospheric\n"
          "levels with pressures lower than the given value (default: no\n"
          "removal). This reduces computational burden and is useful when\n"
          "upper atmospheric contributions are negligible.\n"
          "\n"
          "Possible future extensions: Add a keyword parameter to refine the\n"
          "pressure grid if it is too coarse. Or a version that interpolates\n"
          "onto given grids, instead of using and returning the original grids.\n"),
      AUTHORS("Jana Mendrok, Manfred Brath"),
      OUT("p_grid",
          "lat_grid",
          "lon_grid",
          "t_field",
          "z_field",
          "vmr_field",
          "particle_bulkprop_field",
          "particle_bulkprop_names"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species", "atm_fields_compact", "atmosphere_dim"),
      GIN("delim", "p_min", "check_gridnames"),
      GIN_TYPE("String", "Numeric", "Index"),
      GIN_DEFAULT("-", "0", "0"),
      GIN_DESC(/* delim */
               "Delimiter string of *scat_species* elements.",
               /* p_min */
               "Minimum-pressure level to consider (for TOA).",
               /* check_gridnames */
               "A flag with value 1 or 0. If set to one, the gridnames of \n"
               " the *atm_fields_compact* are checked.")));

  md_data_raw.push_back(MdRecord(
      NAME("atmfields_checkedCalc"),
      DESCRIPTION(
          "Checks consistency of (clear sky) atmospheric fields.\n"
          "\n"
          "The following WSVs are treated: *p_grid*, *lat_grid*, *lon_grid*,\n"
          "*t_field*, *vmr_field*, wind_u/v/w_field and mag_u/v/w_field.\n"
          "\n"
          "If any of the variables above is changed, then this method shall be\n"
          "called again (no automatic check that this is fulfilled!).\n"
          "\n"
          "The tests include that:\n"
          " 1. Atmospheric grids (p/lat/lon_grid) are OK with respect to\n"
          "    *atmosphere_dim* (and vmr_field also regarding *abs_species*).\n"
          " 2. Atmospheric fields have sizes consistent with the atmospheric\n"
          "    grids.\n"
          " 3. *abs_f_interp_order* is not zero if any wind is nonzero.\n"
          " 4. All values in *t_field* are > 0.\n"
          "\n"
          "Default is that values in *vmr_field* are demanded to be >= 0\n"
          "(ie. zero allowed, in contrast to *t_field*), but this\n"
          "requirement can be removed by the *negative_vmr_ok* argument.\n"
          "\n"
          "If any test fails, there is an error. Otherwise,\n"
          "*atmfields_checked* is set to 1.\n"
          "\n"
          "The cloudbox is covered by *cloudbox_checked*, *z_field* is\n"
          "part of the checks done around *atmgeom_checked*.\n"
          "\n"
          "If you choose to use *bad_partition_functions_ok* please note that\n"
          "this is done on your own risk and that it could introduce hard-to-\n"
          "track errors into your calculations.  Do not use this for anything\n"
          "important.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("atmfields_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "abs_species",
         "t_field",
         "vmr_field",
         "wind_u_field",
         "wind_v_field",
         "wind_w_field",
         "mag_u_field",
         "mag_v_field",
         "mag_w_field",
         "partition_functions",
         "abs_f_interp_order"),
      GIN("negative_vmr_ok", "bad_partition_functions_ok"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT("0", "0"),
      GIN_DESC("Flag whether to accept vmr_field < 0.",
               "Flag whether to accept partition functions not covering"
               " *t_field* range.")));

  md_data_raw.push_back(MdRecord(
      NAME("atmgeom_checkedCalc"),
      DESCRIPTION(
          "Checks consistency of geometric considerations of the atmosphere.\n"
          "\n"
          "The following WSVs are checked: *z_field*, *refellipsoid*, *z_surface*,\n"
          "*lat_true* and *lon_true*. If any of the variables above is changed,\n"
          "then this method shall be called again (no automatic check that this is\n"
          "fulfilled!).\n"
          "\n"
          "The tests include that:\n"
          " 1. *refellipsoid* has correct size, and that eccentricity is\n"
          "    set to zero if 1D atmosphere.\n"
          " 2. *z_field* and *z_surface* have sizes consistent with the\n"
          "    atmospheric grids.\n"
          " 3. There is no gap between *z_surface* and *z_field*.\n"
          "\n"
          "*lat_true* and *lon_true* are allowed to be empty.\n"
          "\n"
          "If any test fails, there is an error. Otherwise, *atmgeom_checked*\n"
          "is set to 1.\n"
          "\n"
          "See further *atmgeom_checkedCalc*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("atmgeom_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "z_field",
         "refellipsoid",
         "z_surface",
         "lat_true",
         "lon_true"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("AtmosphereSet1D"),
      DESCRIPTION(
          "Sets the atmospheric dimension to 1D.\n"
          "\n"
          "Sets *atmosphere_dim* to 1, and the latitude and longitude grids\n"
          "are set to be empty.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("atmosphere_dim", "lat_grid", "lon_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("AtmosphereSet2D"),
      DESCRIPTION(
          "Sets the atmospheric dimension to be 2D.\n"
          "\n"
          "Sets *atmosphere_dim* to 2 and the longitude grid to be empty.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("atmosphere_dim", "lon_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("AtmosphereSet3D"),
      DESCRIPTION(
          "Sets the atmospheric dimension to 3D.\n"
          "\n"
          "Sets *atmosphere_dim* to 3, and *lat_true* and *lon_true* are\n"
          "set to be empty.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("atmosphere_dim", "lat_true", "lon_true"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("AtmRawRead"),
      DESCRIPTION(
          "Reads atmospheric data from a scenario.\n"
          "\n"
          "An atmospheric scenario includes the following data for each\n"
          "position (pressure, latitude, longitude) in the atmosphere:\n"
          "   1. temperature field\n"
          "   2. the corresponding altitude field\n"
          "   3. vmr fields for the absorption species\n"
          "The vmr fields read are governed by the species given in\n"
          "*abs_species*. Beside gaseous species, these can also contain\n"
          "purely absorbing particulate matter. In the latter case the\n"
          "profiles are supposed to provide the mass content (unit kg/m3) for\n"
          "clouds and precipitation rate (unit kg/m2/s) for precipitation\n"
          "instead of the vmr.\n"
          "\n"
          "The data is stored in different files. This methods reads all\n"
          "files and creates the variables *t_field_raw*, *z_field_raw* and\n"
          "*vmr_field_raw*.  *nlte_field_raw* is set to empty.\n"
          "\n"
          "Files in a scenarios should be named matching the pattern of:\n"
          "basename.speciesname.xml\n (for temperature and altitude the\n"
          "expected 'speciesname' are 't' and'z', respectivly)."
          "\n"
          "The files can be anywhere, but they must all be in the same\n"
          "directory, selected by 'basename'. The files are chosen by the\n"
          "species name. If you have more than one tag group for the same\n"
          "species, the same profile will be used.\n"),
      AUTHORS("Claudia Emde"),
      OUT("t_field_raw",
          "z_field_raw",
          "vmr_field_raw",
          "nlte_field_raw",
          "nlte_level_identifiers"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("basename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Name of scenario, probably including the full path. For "
               "example: \"/smiles_local/arts-data/atmosphere/fascod/"
               "tropical\"")));

  md_data_raw.push_back(MdRecord(
      NAME("AtmWithNLTERawRead"),
      DESCRIPTION(
          "Reads atmospheric data from a scenario.\n"
          "\n"
          "An atmospheric scenario includes the following data for each\n"
          "position (pressure, latitude, longitude) in the atmosphere:\n"
          "   1. temperature field\n"
          "   2. the corresponding altitude field\n"
          "   3. vmr fields for the gaseous species\n"
          "   4. Non-LTE temperature fields and matching identifiers\n"
          "The data is stored in different files. This method reads all\n"
          "files and creates the variables *t_field_raw*, *z_field_raw*,\n"
          "*vmr_field_raw*, *nlte_field_raw*, and *nlte_level_identifiers*.\n"
          "\n"
          "Files in a scenarios should be named matching the pattern of:\n"
          "tropical.H2O.xml\n"
          "\n"
          "The files can be anywhere, but they must be all in the same\n"
          "directory, selected by 'basename'. The files are chosen by the\n"
          "species name. If you have more than one tag group for the same\n"
          "species, the same profile will be used.\n"),
      AUTHORS("Claudia Emde", "Richard Larsson"),
      OUT("t_field_raw",
          "z_field_raw",
          "vmr_field_raw",
          "nlte_field_raw",
          "nlte_level_identifiers"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN("basename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Name of scenario, probably including the full path. For "
               "example: \"/smiles_local/arts-data/atmosphere/fascod/"
               "tropical\"")));

  md_data_raw.push_back(MdRecord(
      NAME("atm_fields_compactAddConstant"),
      DESCRIPTION(
          "Adds a constant field to atm_fields_compact.\n"
          "\n"
          "This is handy, e.g., for nitrogen or oxygen. The constant value can\n"
          "be appended or prepended as an additional field to the already\n"
          "existing collection of fields. All dimensions (pressure, latitude,\n"
          "longitude) are filled up, so this works for 1D, 2D, or 3D\n"
          "atmospheres.\n"
          "\n"
          "The passed *name* of the field has to be in accordance with the\n"
          "tagging structure described for *atm_fields_compact*.\n"
          "\n"
          "A list of condensibles can be optionally specified if the VMR of\n"
          "the added species is assuming dry air. The VMR of the added species\n"
          "is then scaled down by the sum of the condensibles' VMR:\n"
          "VMR * (1 - VMR_sum_of_condensibles).\n"
          "For Earth this should be set to [\"abs_species-H2O\"]\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("atm_fields_compactAddSpecies"),
      DESCRIPTION(
          "Adds a field to atm_fields_compact, with interpolation.\n"
          "\n"
          "This method appends or prepends a *GriddedField3* to *atm_fields_compact*.\n"
          "The *GriddedField3* is interpolated upon the grid of\n"
          "*atm_fields_compact*. A typical use case for this method may be to\n"
          "add a climatology of some gas when this gas is needed for radiative\n"
          "transfer calculations, but not yet present in *atm_fields_compact*.\n"
          "One case where this happens is when using the Chevalier91L dataset\n"
          "for infrared simulations.\n"
          "\n"
          "The grids in *atm_fields_compact* must fully encompass the grids in\n"
          "the *GriddedField3* to be added, for interpolation to succeed. If\n"
          "this is not the case, a RuntimeError is thrown.\n"
          "\n"
          "The passed *name* of the field has to be in accordance with the\n"
          "tagging structure described for *atm_fields_compact*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("atm_fields_compactCleanup"),
      DESCRIPTION(
          "Removes unrealistically small or erroneous data from\n"
          "*atm_fields_compact* (or other GriddedField4 data)\n"
          "\n"
          "This WSM checks if the data in *atm_fields_compact* contains\n"
          "values smaller than the given *threshold*. In this case, these\n"
          "values will be set to zero.\n"
          "\n"
          "The method should be applied if *atm_fields_compact* contains\n"
          "unrealistically small or erroneous data (NWP/GCM model data\n"
          "occassionally contains negative values, which are numerical\n"
          "artefacts rather than physical values.)\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("atm_fields_compactCreateFromField"),
      DESCRIPTION(
          "Initiates *atm_fields_compact* from a field.\n"
          "\n"
          "*atm_fields_compact* will have the same size and grids as the GriddedField3,\n"
          "but with one dimension as length 1.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("atm_fields_compactFromMatrix"),
      DESCRIPTION(
          "Sets *atm_fields_compact* from 1D fields given in form of a matrix.\n"
          "\n"
          "For batch calculations it is handy to store atmospheric\n"
          "profiles in an array of matrix. We take such a matrix, and create\n"
          "*atm_fields_compact* from it.\n"
          "\n"
          "The matrix must contain one row for each pressure level.\n"
          "\n"
          "Not all fields contained in the matrix must be selected into\n"
          "*atm_fields_compact*, but the selection must at least contain\n"
          "fields of pressure, temperature, altitude and one absorption\n"
          "species.\n"
          "The matrix can contain some additional fields which are not\n"
          "directly used by ARTS for calculations but can be required for\n"
          "further processing, e.g. wind speed and direction. These fields do\n"
          "not need to be transfered into the *atm_fields_compact* variable.\n"
          "\n"
          "Selection of fields into *atm_fields_compact* works by providing a\n"
          "field name tag in *field_names* for the selected fields, while\n"
          "unselected fields are tagged by 'ignore'. Order of tags in\n"
          "*field_names* is strictly taken as corresponding to column order in\n"
          "the matrix.\n"
          "The pressure fields are by convention the first column of the\n"
          "matrix, hence must not be tagged. That is, there must be given one\n"
          "field name tag less than matrix columns.\n"
          "\n"
          "For detailed tagging conventions see *atm_fields_compact*.\n"
          "\n"
          "Works only for *atmosphere_dim==1.*\n"),
      AUTHORS("Stefan Buehler", "Daniel Kreyling", "Jana Mendrok"),
      OUT("atm_fields_compact"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim"),
      GIN("gin1", "field_names"),
      GIN_TYPE("Matrix", "ArrayOfString"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("One atmosphere matrix from batch input ArrayOfMatrix.",
               "Order/names of atmospheric fields.")));

  md_data_raw.push_back(MdRecord(
      NAME("backend_channel_responseFlat"),
      DESCRIPTION(
          "Sets up a rectangular channel response.\n"
          "\n"
          "The response of the backend channels is hee assumed to be constant\n"
          "inside the resolution width, and zero outside.\n"
          "\n"
          "The method assumes that all channels have the same response.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("backend_channel_responseGaussian"),
      DESCRIPTION(
          "Sets up a gaussian backend channel response.\n"
          "\n"
          "The method assumes that all channels can be modelled as gaussian.\n"
          "\n"
          "If *fwhm* has only one element, all channels are assumed to be equal.\n"
          "If *fwhm* has multiple elements, *xwidth_si* and *dx_si* must have one\n"
          "element or the same number of elements as *fwhm*. If one element is given,\n"
          "this value will be used for all channels.\n"
          "\n"
          "The grid generated can be written as\n"
          "   si * [-xwidth_si:dx_si:xwidth_si]\n"
          "where si is the standard deviation corresponding to the FWHM.\n"
          "That is, width and spacing of the grid is specified in terms of\n"
          "number of standard deviations. If xwidth_si is set to 2, the\n"
          "response will cover about 95% the complete response. For\n"
          "xwidth_si=3, about 99% is covered. If xwidth_si/dx_si is not\n"
          "an integer, the end points of the grid are kept and the spacing\n"
          "if the grid is adjusted in the downward direction (ie. spacing is.\n"
          "is max *dx_si*).\n"),
      AUTHORS("Patrick Eriksson, Oliver Lemke"),
      OUT("backend_channel_response"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("fwhm", "xwidth_si", "dx_si"),
      GIN_TYPE("Vector", "Vector", "Vector"),
      GIN_DEFAULT(NODEF, "[3]", "[0.1]"),
      GIN_DESC("Full width at half-maximum",
               "Half-width of response, in terms of std. dev.",
               "Grid spacing, in terms of std. dev.")));

  md_data_raw.push_back(MdRecord(
      NAME("batch_atm_fields_compactAddConstant"),
      DESCRIPTION("Adds a constant field to batch_atm_fields_compact.\n"
                  "\n"
                  "Applies *atm_fields_compactAddConstant* to each batch.\n"
                  "The format is equal to that WSM.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("batch_atm_fields_compactAddSpecies"),
      DESCRIPTION(
          "Adds a field to *batch_atm_fields_compact*, with interpolation.\n"
          "\n"
          "This method appends or prepends a *GriddedField3* to each *atm_fields_compact*.\n"
          "in *batch_atm_fields_compact*. For details, see *atm_fields_compactAddSpecies*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("batch_atm_fields_compactCleanup"),
      DESCRIPTION(
          "Removes unrealistically small or erroneous data from each data field\n"
          "of *batch_atm_fields_compact* (or other AerrayOfGriddedField4 data)\n"
          "\n"
          "This WSM checks if the data in *batch_atm_fields_compact* contains\n"
          "values smaller than the given *threshold*. In this case, these\n"
          "values will be set to zero.\n"
          "\n"
          "The method should be applied if *batch_atm_fields_compact* contains\n"
          "unrealistically small or erroneous data (NWP/GCM model data\n"
          "occassionally contains negative values, which are numerical\n"
          "artefacts rather than physical values.)\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("batch_atm_fields_compactFromArrayOfMatrix"),
      DESCRIPTION(
          "Expand batch of 1D atmospheric state matrices to batch_atm_fields_compact.\n"
          "\n"
          "This is used to handle 1D batch cases, e.g. from NWP/GCM model like\n"
          "the Chevallier91L data set, stored in a matrix (it is preferred,\n"
          "though, to immediatedly store the model fields as\n"
          "*ArrayOfGriddedField4* and use *ReadXML* to load them directly into\n"
          "*batch_atm_fields_compact*).\n"
          "\n"
          "Works only for *atmosphere_dim==1.*\n"
          "\n"
          "See *atm_fields_compactFromMatrix* for basic documentation.\n"
          "\n"
          "See *batch_atm_fields_compactAddConstant* and\n"
          "batch_atm_fields_compactAddSpecies* for adding additional fields.\n"),
      AUTHORS("Stefan Buehler", "Daniel Kreyling", "Jana Mendrok"),
      OUT("batch_atm_fields_compact"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim"),
      GIN("atmospheres_fields", "field_names"),
      GIN_TYPE("ArrayOfMatrix", "ArrayOfString"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Batch of atmospheres stored in one array of matrix",
               "Order/names of atmospheric fields.")));

  md_data_raw.push_back(MdRecord(
      NAME("CIAInfo"),
      DESCRIPTION(
          "Display information about the given CIA tags.\n"
          "The CIA tags shown are in the same format as needed by *abs_speciesSet*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("CIARecordReadFromFile"),
      DESCRIPTION("Reads CIARecord from Hitran-style file.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("cloudboxOff"),
      DESCRIPTION(
          "Deactivates the cloud box.\n"
          "\n"
          "Use this method if no scattering calculations shall be performed.\n"
          "\n"
          "The function sets *cloudbox_on* to 0, *cloudbox_limits*,\n"
          "*pnd_field*, *scat_data*, *scat_data_raw*, *iy_cloudbox_agenda*\n"
          "and *particle_masses* to be empty and sizes *dpnd_field_dx* to be\n"
          "consitent with *jacobian_quantities*.\n"),
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
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("cloudboxSetAutomatically"),
      DESCRIPTION(
          "Sets the cloud box to encompass the cloud given by the entries\n"
          "in *particle_field*.\n"
          "\n"
          "This WSM handles one *Tensor4* type *particle_field* at a time. It can\n"
          "be used to determine the cloudbox from *particle_bulkprop_field*,\n"
          "but also from the various *scat_species_XXX_field* (or even from a\n"
          "read-in *pnd_field*). For the latter, the WSM needs to be called\n"
          "once per each *scat_species_XXX_field*, with previously determined\n"
          "*cloudbox_limits* provided through *cloudbox_limits_old*.\n"
          "\n"
          "The function must be called before executing any WSM that applies\n"
          "*cloudbox_limits*.\n"
          "\n"
          "The function iterates over all 3D fields in *particle_field* (which\n"
          "might correspond to the different scattering species as in\n"
          "*scat_species_XX_field* or to different particle bulk properties as\n"
          "in *particle_bulkprop_field*). Each field is searched for the first\n"
          "and last pressure index, where the value is unequal to zero. This\n"
          "index is then copied to *cloudbox_limits*.\n"
          "If *particle_field* is empty, the cloudbox is switched off\n"
          "(*cloudbox_on*=0).\n"
          "\n"
          "Additionaly the lower cloudbox_limit is altered by *cloudbox_margin*.\n"
          "The margin is given as a height difference in meters and transformed\n"
          "into a pressure (via isothermal barometric height formula). This\n"
          "alteration is to ensure covering photons that leave the cloud, but\n"
          "reenter through a limb path.\n"
          "If *cloudbox_margin* is set to -1 (default), the cloudbox will extend\n"
          "to the surface. Hence, the lower cloudbox_limit is set to 0 (index\n"
          "of first pressure level).\n"
          "*cloudbox_margin* will be applied on each call of the WSM. Hence,\n"
          "if called successively, e.g. for several *scat_species_XXX_field*,\n"
          "and *cloudbox_margin* is not -1, it is suggested to apply the\n"
          "desired *cloudbox_margin* only for the last WSM call, while for the\n"
          "others set *cloudbox_margin* to 0.\n"
          "\n"
          "Works only for *atmosphere_dim==1.*\n"),
      AUTHORS("Jana Mendrok, Daniel Kreyling"),
      OUT("cloudbox_on", "cloudbox_limits"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim", "p_grid", "lat_grid", "lon_grid"),
      GIN("particle_field", "cloudbox_limits_old", "cloudbox_margin"),
      GIN_TYPE("Tensor4", "ArrayOfIndex", "Numeric"),
      GIN_DEFAULT(NODEF, "[-1]", "-1"),
      GIN_DESC("A collection of particle property fields (e.g."
               " *particle_bulkprop_field*,"
               " *scat_species_mass_density_field*).",
               "Preset cloudbox limits, e.g. resulting from a previous run"
               " of *cloudboxSetAutomatically*.",
               "Minimum distance [m] between lowest 'cloudy' level and"
               " cloudbox lower limit. If set to *-1* (default), the"
               " cloudbox lower limit is fixed to 0, i.e., corresponds to"
               " the lowest atmospheric level (or the surface).")));

  md_data_raw.push_back(
      MdRecord(NAME("cloudboxSetFullAtm"),
               DESCRIPTION("Sets the cloudbox to cover the full atmosphere.\n"),
               AUTHORS("Claudia Emde, Jana Mendrok"),
               OUT("cloudbox_on", "cloudbox_limits"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("atmosphere_dim", "p_grid", "lat_grid", "lon_grid"),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("cloudboxSetManually"),
      DESCRIPTION(
          "Sets the cloud box to encompass the given positions.\n"
          "\n"
          "The function sets *cloudbox_on* to 1 and sets *cloudbox_limits*\n"
          "following the given pressure, latitude and longitude positions.\n"
          "The index limits in *cloudbox_limits* are selected to give the\n"
          "smallest possible cloud box that encompass the given points.\n"
          "\n"
          "The points must be given in the same order as used in\n"
          "*cloudbox_limits*. That means that the first keyword argument\n"
          "shall be a higher pressure than argument two, while the latitude\n"
          "and longitude points are given in increasing order. Positions\n"
          "given for dimensions not used by the selected atmospheric\n"
          "dimensionality are ignored.\n"
          "\n"
          "The given pressure points can be outside the range of *p_grid*.\n"
          "The pressure limit is then set to the end point of *p_grid*.\n"
          "The given latitude and longitude points must be inside the range\n"
          "of the corresponding grid. In addition, the latitude and longitude\n"
          "points cannot be inside the outermost grid ranges as the latitude\n"
          "and longitude limits in *cloudbox_limits* are not allowed to be\n"
          "grid end points.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("cloudbox_on", "cloudbox_limits"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim", "p_grid", "lat_grid", "lon_grid"),
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

  md_data_raw.push_back(MdRecord(
      NAME("cloudboxSetManuallyAltitude"),
      DESCRIPTION(
          "Sets the cloud box to encompass the given positions.\n"
          "\n"
          "As *cloudboxSetManually* but uses altitudes instead of pressure.\n"
          "The given altitude points can be outside the range of *z_field*.\n"
          "The altitude limit is then set to the end point of *p_grid*.\n"),
      AUTHORS("Claudia Emde"),
      OUT("cloudbox_on", "cloudbox_limits"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim", "z_field", "lat_grid", "lon_grid"),
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

  md_data_raw.push_back(MdRecord(
      NAME("cloudbox_checkedCalc"),
      DESCRIPTION(
          "Checks consistency and validity of the cloudbox governing variables.\n"
          "\n"
          "The following WSVs are treated: *cloudbox_on*, *cloudbox_limits*,\n"
          "*pnd_field*, *scat_data*, *scat_species*, *abs_species*, *particle_masses*\n"
          "*particle_bulkprop_field*, *particle_bulkprop_names* and wind_u/v/w_field.\n"
          "\n"
          "If any of these variables is changed, then this method shall be\n"
          "called again (no automatic check that this is fulfilled!).\n"
          "\n"
          "The main checks are if the cloudbox limits are OK with respect to\n"
          "the atmospheric dimensionality and the limits of the atmosphere,\n"
          "and that the scattering element variables *pnd_field* and\n"
          "*scat_data* match in size.\n"
          "\n"
          "Further checks on *scat_data* are performed in *scat_data_checkedCalc*\n"
          "\n"
          "*scat_species* and *particle_masses* must either be empty or have a\n"
          "size that matches the other data. If non-empty, some check of these\n"
          "variables are performed.\n"
          "\n"
          "If any test fails, there is an error. Otherwise, *cloudbox_checked*\n"
          "is set to 1.\n"),
      AUTHORS("Patrick Eriksson, Jana Mendrok"),
      OUT("cloudbox_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmfields_checked",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "z_field",
         "z_surface",
         "wind_u_field",
         "wind_v_field",
         "wind_w_field",
         "cloudbox_on",
         "cloudbox_limits",
         "pnd_field",
         "dpnd_field_dx",
         "jacobian_quantities",
         "scat_data",
         "scat_species",
         "particle_masses",
         "abs_species"),
      GIN("negative_pnd_ok"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("Flag whether to accept pnd_field < 0.")));

  md_data_raw.push_back(MdRecord(
      NAME("Compare"),
      DESCRIPTION(
          "Checks the consistency between two variables.\n"
          "\n"
          "The two variables are checked to not deviate outside the specified\n"
          "value (*maxabsdiff*). An error is issued if this is not fulfilled.\n"
          "\n"
          "The main application of this method is to be part of the test\n"
          "control files, and then used to check that a calculated value\n"
          "is consistent with an old, reference, value.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("CompareRelative"),
      DESCRIPTION(
          "Checks the consistency between two variables by their relative values.\n"
          "\n"
          "The two variables are checked to not deviate outside the specified\n"
          "relative value (*maxabsreldiff*). An error is issued if this is not\n"
          "fulfilled.\n"
          "\n"
          "The main application of this method is to be part of the test\n"
          "control files, and then used to check that a calculated value\n"
          "is consistent with an old, reference, value.\n"
          "\n"
          "If either value is 0.0, the relative error is considered as 0\n"
          "for easier use.  This really means infinite differences, though\n"
          "allowing zero-crossings is useful for plenty of tests. So Be Aware!\n"
          "\n"
          "If both *var1* and *var2* are non-zero, the difference is evaluated\n"
          "as: abs(var1/var2-1)\n"
          "That is, *var2* is taken as the reference value.\n"),
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
          "ArrayOfArrayOfMatrix, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6,"
          "ArrayOfPropagationMatrix, ArrayOfArrayOfPropagationMatrix,"
          "ArrayOfStokesVector, ArrayOfArrayOfStokesVector,",
          "Numeric, Vector, Matrix, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7,"
          "ArrayOfVector, ArrayOfMatrix, ArrayOfTensor3, ArrayOfTensor4,"
          "ArrayOfTensor6, ArrayOfTensor7, ArrayOfArrayOfVector,"
          "ArrayOfArrayOfMatrix, ArrayOfArrayOfTensor3, ArrayOfArrayOfTensor6,"
          "ArrayOfPropagationMatrix, ArrayOfArrayOfPropagationMatrix,"
          "ArrayOfStokesVector, ArrayOfArrayOfStokesVector,",
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

  md_data_raw.push_back(MdRecord(
      NAME("complex_refr_indexConstant"),
      DESCRIPTION(
          "Set complex refractive index to a constant value.\n"
          "\n"
          "Frequency and temperature grids are set to have length 1 (and\n"
          "set to the value 0).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("complex_refr_indexIceMatzler06"),
      DESCRIPTION(
          "Refractive index of ice following Matzler06 parameterization.\n"
          "\n"
          "Calculates temperature dependent complex refractive index of\n"
          "hexagonal ice at microwave and sub-mm frequencies (10MHz-3Tz).\n"
          "\n"
          "This parametrization is also applied by the microwave and\n"
          "submm-wave part of the Warren08 model.\n"
          "\n"
          "References:\n"
          "Matzler, C., 2006: Thermal Microwave Radiation: Application for\n"
          "Remote Sensing, Microwave dielectric properties of ice, pp. 455-462,\n"
          "Inst. Eng. Technol., Stevenage, U. K.\n"
          "Warren, S. G., and R. E. Brandt, 2008: Optical constants of ice\n"
          "from the ultraviolet to the microwave: A revised compilation,\n"
          "J. Geophys. Res., 113, D14220, doi:10.1029/2007JD009744.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("complex_refr_indexIceWarren84"),
      DESCRIPTION(
          "Refractive index of ice following Warren84 parameterization.\n"
          "\n"
          "Calculates complex refractive index of Ice 1H for wavelengths\n"
          "between 45 nm and 8.6 m.\n"
          "For wavelengths above 167 microns, temperature dependence is\n"
          "included for temperatures between 213 and 272K.\n"
          "Mainly intended for applications in Earth ice\n"
          "clouds and snow, not other planets or interstellar space;\n"
          "the temperature dependence or crystalline form of ice may be\n"
          "incorrect for these latter applications.\n"
          "\n"
          "Authors of Fortran function:\n"
          "Stephen Warren, Univ. of Washington (1983)\n"
          "Bo-Cai Gao, JCESS, Univ. of Maryland (1995)\n"
          "Warren Wiscombe, NASA Goddard (1995)\n"
          "\n"
          "References:\n"
          "Warren, S., 1984: Optical Constants of Ice from the Ultraviolet\n"
          "to the Microwave, Appl. Opt. 23, 1206-1225\n"
          "\n"
          "Kou, L., D. Labrie, and P. Chylek, 1994: Refractive indices\n"
          "of water and ice in the 0.65- to 2.5-micron spectral range,\n"
          "Appl. Opt. 32, 3531-3540\n"
          "\n"
          "Perovich, D., and J. Govoni, 1991: Absorption Coefficients\n"
          "of Ice from 250 to 400 nm, Geophys. Res. Lett. 18, 1233-1235\n"),
      AUTHORS("Oliver Lemke"),
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

  md_data_raw.push_back(MdRecord(
      NAME("complex_refr_indexWaterLiebe93"),
      DESCRIPTION(
          "Complex refractive index of liquid water according to Liebe 1993.\n"
          "\n"
          "The method treats liquid water without salt. Thus, not valid below\n"
          "10 GHz. Upper frequency limit not known, here set to 1000 GHz.\n"
          "Model parameters taken from Atmlab function epswater93 (by\n"
          "C. Maetzler), which refer to Liebe 1993 without closer\n"
          "specifications.\n"
          "\n"
          "Temperatures must be between -40 and 100 degrees Celsius. The\n"
          "accuracy of the parametrization below 0 C is not known by us.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("Copy"),
      DESCRIPTION(
          "Copy a workspace variable.\n"
          "\n"
          "This method can copy any workspace variable\n"
          "to another workspace variable of the same group. (E.g., a Matrix to\n"
          "another Matrix.)\n"
          "\n"
          "As always, output comes first in the argument list!\n"
          "\n"
          "Usage example:\n"
          "\n"
          "Copy(f_grid, p_grid)\n"
          "\n"
          "Will copy the content of *p_grid* to *f_grid*. The size of *f_grid*\n"
          "is adjusted automatically (the normal behaviour for workspace\n"
          "methods).\n"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Any"),
      GOUT_DESC("Destination variable."),
      IN(),
      GIN("in"),
      GIN_TYPE("Any"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Source variable."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(MdRecord(
      NAME("covmat1D"),
      DESCRIPTION(
          "Create 1D covariance matrix.\n"
          "\n"
          "Creates a 1D covariance matrix for two retrieval quantities on given \n"
          " grids from a given functional form. Elements  of the covariance matrix\n"
          "are computed as\n"
          " S_{i,j} = sigma_i * sigma_j * f(d_{i,j} / l_{i,j}) \n"
          " where d_{i,j} is the distance between the two grid points and l_{i,j}\n"
          " the mean of the correlation lengths of the grid points.\n"
          "\n"
          " If a cutoff value co is given elements with absolute value less than this \n"
          " are set to zero.\n"
          "\n"
          "The following functional forms are available:\n"
          "  \"exp\": f(x) = exp(-x) \n"
          "  \"lin\": f(x) = 1.0 - x, for x > 1.0, 0.0 otherwise \n"
          "  \"gauss\": f(x) = exp(-x^2) \n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("out"),
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
               "The name of the functional form to use.", ),
      PASSWORKSPACE(false),
      SETMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(MdRecord(
      NAME("covmat1DMarkov"),
      DESCRIPTION(
          "Create Markov Process Covariance Matrix.\n"
          "\n"
          "Create a markov process covariance matrix for a retrieval quantity on \n"
          " evenly spaced 1D grid. The correlation between two grid points i,j is \n"
          " is computed as \n"
          " cov(i,j) = sigma[i] * sigma[j] * exp(- d(i,j) / lc)\n"
          " where d(i,j) = abs(grid[i] - grid[j]).\n"
          "\n"
          "This function also sets covmat_inv_block to the analytically computed inverse\n"
          "of the covariance matrix of the markov provess, which is tri-diagonal. Note\n"
          "that this requires the retrieval grid to be evenly spaced.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("out", "out_inverse"),
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

  md_data_raw.push_back(MdRecord(
      NAME("covmatDiagonal"),
      DESCRIPTION(
          "Sets the matrix in covmat_block to a diagonal matrix with the variances\n"
          "provided in *vars* as diagonal elements."
          "\n"
          "Also sets covmat_block_inv to the inverse of the block so that the\n"
          "computation of the inverse is avoided.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("out", "out_inverse"),
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

  md_data_raw.push_back(MdRecord(
      NAME("covmat_seAddBlock"),
      DESCRIPTION(
          "Add a block to the measurement covariance matrix *covmat_se*\n"
          "\n"
          "This functions adds a given dense or sparse matrix as block to the covariance\n"
          "matrix *covmatrix_sx*. The position of the block can be given by the generic\n"
          "arguments *i* and *j*. Note that diagonal blocks must be added in order starting from\n"
          " in  the top left corner. If an off-diagonal block is added it must have corresponding\n"
          " existing blocks on the diagonal and these must be consistent with the dimensions\n"
          " of the block.  If *i* and *j*  are not provided, the blok will be added\n"
          "at the first free spot on the diagonal.\n"),
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
               "Index of a retrieval quantity. Must satisfy *i* <= *j*.",
               "Index of a retrieval quantity. Must satisfy *i* <= *j*."),
      PASSWORKSPACE(false),
      SETMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(MdRecord(
      NAME("covmat_seAddInverseBlock"),
      DESCRIPTION(
          "Add the inverse of a block to covariance matrix *covmat_se*\n"
          "\n"
          "This functions adds a given matrix as the inverse of a block in the covariance\n"
          "matrix *covmatrix_se*. The purpose of this function is to allow the user to\n"
          "to use a precomputed inverse for this block in the covariance matrix, that may\n"
          "for example have been obtained analytically.\n"
          "\n"
          "This function requires the corresponding non-inverse block to already be present in *covmat_se*"
          "\n"
          "\n Note that for this to work this retrieval quantity must be independent from\n"
          "other retrieval quantities that do not have an inverse. Otherwise the inverse\n"
          "will be ignored and recomputed numerically.\n"
          "\n"
          "For the rest, the same requirements as for *covmat_seAddBlock* apply.\n"),
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
               "Index of a retrieval quantity. Must satisfy *i* <= *j*.",
               "Index of a retrieval quantity. Must satisfy *i* <= *j*."),
      PASSWORKSPACE(false),
      SETMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(MdRecord(
      NAME("covmat_seSet"),
      DESCRIPTION(
          "Set covmat_se to a given matrix.\n"
          "\n"
          "This sets the measurement covariance matrix *covmat_se* to\n"
          "the matrix given by the generic input *covmat*. The covariance\n"
          "matrix can be of type CovarianceMatrix, Matrix or Sparse.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("covmat_sxSet"),
      DESCRIPTION(
          "Set covmat_sx to a given matrix.\n"
          "\n"
          "This sets the measurement covariance matrix *covmat_sx* to\n"
          "the matrix given by the generic input *covmat*. The covariance\n"
          "matrix can be of type CovarianceMatrix, Matrix or Sparse.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("covmat_sxAddBlock"),
      DESCRIPTION(
          "Add a block to the a priori covariance matrix *covmat_sx*\n"
          "\n"
          "This functions adds a given matrix as a block in the covariance\n"
          "matrix *covmatrix_sx*. The position of the block can be given by the generic\n"
          "arguments *i* and *j*, which should give the index of the retrieval quantity in\n"
          "*jacobian_quantities*, which is given just by the order the quantities have been\n"
          "added to the retrieval.\n"
          "\n"
          "If arguments *i* and *j* are omitted, the block will be added as diagonal block\n"
          "for the last added retrieval quantity.\n"
          "\n"
          "If provided, the index *i* must be less than or equal to *j*. Also the provided\n"
          "block must be consistent with the corresponding retrieval quantities.\n"),
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
               "Index of a retrieval quantity. Must satisfy *i* <= *j*.",
               "Index of a retrieval quantity. Must satisfy *i* <= *j*."),
      PASSWORKSPACE(false),
      SETMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(MdRecord(
      NAME("covmat_sxAddInverseBlock"),
      DESCRIPTION(
          "Add the inverse of a block in covariance matrix *covmat_sx*\n"
          "\n"
          "This functions adds a given matrix as the inverse of a block in the covariance\n"
          "matrix *covmatrix_sx*. The purpose of this function is to allow the user to\n"
          "to use a precomputed inverse for this block in the covariance matrix, the may\n"
          "for example by obtained analytically.\n"
          "\n"
          "This function requires the non-inverse block to already be present in *covmat_sx*"
          "\n"
          "\n Note that for this to work this retrieval quantity must be independent from\n"
          "other retrieval quantities that do not have an inverse. Otherwise the inverse\n"
          "will be ignored and recomputed numerically.\n"
          "\n"
          "For the rest, the same requirements as for *covmat_sxAddBlock* apply.\n"),
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
               "Index of a retrieval quantity. Must satisfy *i* <= *j*.",
               "Index of a retrieval quantity. Must satisfy *i* <= *j*."),
      PASSWORKSPACE(false),
      SETMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(MdRecord(
      NAME("covmat_sxExtractSqrtDiagonal"),
      DESCRIPTION(
          "Extract the square root of the diagonal of the state space covariance matrix."
          "\n"
          "This function extracts the diagonal of the state space covariance matrix\n"
          "*covmatrix_sx* and computes its square root. The resulting vector can then\n"
          "be used as *x_norm* argument for the OEM method to avoid scaling problems.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("Delete"),
      DESCRIPTION(
          "Deletes a workspace variable.\n"
          "\n"
          "The variable is marked as uninitialized and its memory freed.\n"
          "It is not removed from the workspace though, therefore you\n"
          "don't need to/can't call Create for this variable again.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("diameter_maxFromDiameter_volume_equ"),
      DESCRIPTION(
          "Calculates maximum and area equivalent diameters from volume\n"
          "equivalent diameter.\n"
          "\n"
          "This is primarily a help function for using the T-matrix method\n"
          "and only a few particle shapes are handled. "
          "\n"
          "For shapes handled and further comments on the input arguments, see\n"
          "*scat_data_singleTmatrix*.\n"
          "\n"
          "Area equivalent diameter is the equivalent sphere diameter\n"
          "corresponding to the \"maximum axial area\". This is the largest\n"
          "cross-sectional area of the particle, observed either along the\n"
          "particle's main axis or in the perpendicular direction. That is,\n"
          "for a cylinder having diameter d and thickness h, this area is\n"
          "either (pi*d^2)/4 or (h*d).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("diameter_volume_equFromDiameter_max"),
      DESCRIPTION(
          "Converts from maximum to volume equivalent diameter.\n"
          "\n"
          "This is primarily a help function for using the T-matrix part\n"
          "and only a few particle shapes are handled. "
          "\n"
          "For shapes handled and further comments on the input arguments,\n"
          "see *scat_data_singleTmatrix*.\n"
          "\n"
          "Also the volume is provided. It is simply sqrt(pi*dveq^3/6).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("DiffZaAa"),
      DESCRIPTION(
          "Derives the difference betwenn zenith and azimuth angles.\n"
          "\n"
          "Determines the difference between a set of angles (*other_los*)\n"
          "and a reference direction (*ref_los*). This method reverses the\n"
          "addition made by *AddZaAa*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("dlos"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Derived differences in line-of-sight."),
      IN(),
      GIN("ref_los", "other_los"),
      GIN_TYPE("Vector", "Matrix"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Reference line-of-sight (a single LOS).",
               "Other line-of-sights (can be multiple LOS).")));

  md_data_raw.push_back(MdRecord(
      NAME("DisortCalc"),
      DESCRIPTION(
          "Interface to the DISORT scattering solver (Stamnes et al.).\n"
          "\n"
          "DISCLAIMER: There is a couple of known issues with the current\n"
          "implementation (see below). Use this WSM with care and only if\n"
          "these limitations/requirements are fulfilled. Results might be\n"
          "erroneous otherwise.\n"
          "\n"
          "DISORT provides the radiation field (*doit_i_field*) from a scalar\n"
          "1D scattering solution assuming a plane-parallel atmosphere (flat\n"
          "Earth). Only totally randomly oriented particles are allowed.\n"
          "Refraction is not taken into account. Only Lambertian surface\n"
          "reflection is handled.\n"
          "\n"
          "*nstreams* is the number of polar angles taken into account\n"
          "internally in the scattering solution, *scat_za_grid* is the\n"
          "polar angle grid on which *doit_i_field* is provided.\n"
          "*nstreams* determines the angular resolution, hence the accuracy,\n"
          "of the scattering solution. The more anisotropic the bulk scattering\n"
          "matrix, the more streams are required. The computational burden\n"
          "increases approximately linearly with *nstreams*. The default value\n"
          "(8) is sufficient for most microwave scattering calculations. It is\n"
          "likely insufficient for IR calculations involving ice clouds,\n"
          "though.\n"
          "*scat_za_grid* determines the resolution of the output radiation\n"
          "field. The size of *scat_za_grid* has no practical impact on\n"
          "computation time in the case of Disort and higher resolution\n"
          "generally improves the interpolation results, hence larger\n"
          "*scat_za_grid* are recommended. To ensure sufficient interpolation\n"
          "accuracy, we require a (hardcoded) minimum size of 38.\n"
          "\n"
          "ARTS-DISORT can be run with different levels of (pseudo-)sphericity,\n"
          "determined by the cloudbox settings.\n"
          "The higher the sphericity level is, the more accurate are the\n"
          "results, but the longer the calculation takes (typically, for\n"
          "downlooking cases - even 50deg off-nadir ones - the differences\n"
          "between all three levels are in the order of sub-K only;\n"
          "differences for limb cases can be significant, though.).\n"
          "DISORT itself always assumes a plane-parallel atmosphere. Different\n"
          "sphericity levels are emulated here by embedding DISORT in\n"
          "different ways and using different output. The available options\n"
          "(from low to high sphericity level) are:\n"
          "- Cloudbox extends over whole atmosphere (e.g. by setting cloudbox\n"
          "  from *cloudboxSetFullAtm*).\n"
          "- Cloudbox extends over a limited part of the atmosphere only (e.g.\n"
          "  by setting cloudbox from *cloudboxSetAutomatically* or\n"
          "  *cloudboxSetManually*). Internally, DISORT is run over the whole\n"
          "  atmosphere, but only the radiation field within the cloudbox is\n"
          "  passed on and used further in ARTS (e.g. by *yCalc*).\n"
          "  This incoming field is internally calculated by ARTS clearsky \n"
          "  methods that take atmospheric sphericity and refractivity fully\n"
          "  into account.\n"
          "\n"
          "Known issues of ARTS implementation:\n"
          "- Surface altitude is not an interface parameter. Surface is\n"
          "  implicitly assumed to be at the lowest atmospheric level.\n"
          "- Scattering angle grids of all scattering elements have to be\n"
          "  identical (except if *pfct_method* is 'interpolate').\n"
          "\n"
          "Keyword *pfct_method* allows to chose the method to extract phase\n"
          "function. 'interpolate' considers temperature dependence. Others\n"
          "neglect it by chosing one specific temperature grid point from the\n"
          "single scattering data: 'low' choses the lowest T-point, 'high' the\n"
          "highest T-point, and 'median' the median T-point. As different\n"
          "scattering elements can have different temperature grids, the actual\n"
          "temperature value used can differ between the scattering elements.\n"
          "Currently, other methods than 'interpolate' require all scattering\n"
          "elements to be given on identical scattering angle grids.\n"
          "Note that this keyword solely affects the phase function;\n"
          "extinction/absorption/scattering cross sections are always\n"
          "interpolated to the actual temperature.\n"),
      AUTHORS("Claudia Emde, Jana Mendrok"),
      OUT("doit_i_field"),
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
         "atmosphere_dim",
         "pnd_field",
         "t_field",
         "z_field",
         "vmr_field",
         "p_grid",
         "scat_data",
         "f_grid",
         "scat_za_grid",
         "stokes_dim",
         "surface_skin_t",
         "surface_scalar_reflectivity"),
      GIN("nstreams", "pfct_method", "Npfct", "quiet"),
      GIN_TYPE("Index", "String", "Index", "Index"),
      GIN_DEFAULT("8", "median", "181", "0"),
      GIN_DESC("Number of polar angle directions (streams) in DISORT "
               "solution (must be an even number).",
               "Flag which method to apply to derive phase function.",
               "Number of angular grid points to calculate bulk phase"
               " function on (and derive Legendre polnomials from). If <0,"
               " the finest za_grid from scat_data will be used.",
               "Silence C Disort warnings.")));

  md_data_raw.push_back(MdRecord(
      NAME("dNdD_F07"),
      DESCRIPTION(
          "Calculation of particle size distribution (dN/dD) following\n"
          "Field et al. (2007) parametrization.\n"
          "\n"
          "A wrapper to internal particle size distribution calculation.\n"
          "Field et al. (2007) provides parametrizations for snow/cloud ice in\n"
          "the tropics and the midlatitudes. The parametrization is in snow\n"
          "(or ice) water mass content (SWC) and ambient atmospheric\n"
          "temperature over particle size in terms of maximum diameter of the\n"
          "snow particles. It furthermore depends on the parameters of the\n"
          "mass-dimension relation m = alpha * Dmax^beta.\n"
          "\n"
          "Negative SWC trigger an error (unless robust=1, where SWC=0 is used\n"
          "internally, hence dNdD=0 is returned).\n"
          "Negative temperatures always trigger an error. Currently, no further\n"
          "temperature validity limits are implemented (although the\n"
          "parametrization itself has been derived from measurements limited\n"
          "to -60C<=T<=0C, i.e. shall strictly only be applied then).\n"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("dNdD"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("size distribution number density [#/m3/m]"),
      IN(),
      GIN("diameter_max", "SWC", "T", "regime", "alpha", "beta", "robust"),
      GIN_TYPE("Vector",
               "Numeric",
               "Numeric",
               "String",
               "Numeric",
               "Numeric",
               "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, "0.0257", "2.0", "0"),
      GIN_DESC(
          "Maximum diameter of the particles [m].",
          "Atmospheric snow water mass content [kg/m3].",
          "Ambient atmospheric temperature [K].",
          "Parametrization regime (\"TR\"=tropical or \"ML\"=midlatitude).",
          "Mass-dimension relationship scaling factor [kg].",
          "Mass-dimension relationship exponent [-].",
          "Flag whether to ignore parametrization value checks.")));

  md_data_raw.push_back(MdRecord(
      NAME("dNdD_H11"),
      DESCRIPTION(
          "Calculation of particle size distribution (dN/dD) following\n"
          "a parametrization by Heymsfield (2011; unpublished).\n"
          "\n"
          "A wrapper to internal particle size distribution calculation.\n"
          "Heymsfield (2011) is an unpublished pre-version of\n"
          "Heymsfield (2013). It is a globally valid parametrization for cloud\n"
          "ice. The parametrization is in ambient atmospheric temperature\n"
          "over particle size in terms of maximum dimension of the\n"
          "scattering elements).\n"
          "\n"
          "Provides only the shape of the number density distribution\n"
          "function. When used by *pnd_fieldCalcFromscat_speciesFields*, it is\n"
          "rescaled to the given ice water content [kg/m3].\n"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("dNdD"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("size distribution number density [#/m3/m]"),
      IN(),
      GIN("Dmax", "t"),
      GIN_TYPE("Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Maximum dimension of the particles [m]",
               "Ambient atmospheric temperature [K]")));

  md_data_raw.push_back(MdRecord(
      NAME("dNdD_H13_Ar"),
      DESCRIPTION(
          "Calculation of particle size and shape distribution (dN/dD, area\n"
          "ratio) following Heymsfield (2013).\n"
          "\n"
          "A wrapper to internal particle size and shape distribution\n"
          "calculation. Heymsfield (2013) is a globally valid parametrization\n"
          "for cloud ice. The parametrization is in ambient atmospheric\n"
          "temperature over particle size in terms of maximum dimension.\n"
          "\n"
          "Provides only the shape of the distribution function of both number\n"
          "density and area ratio. When used by *pnd_fieldCalcFromscat_speciesFields*,\n"
          "number density is rescaled to the given ice water content [kg/m3].\n"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("dNdD", "Ar"),
      GOUT_TYPE("Vector", "Vector"),
      GOUT_DESC("size distribution number density [#/m3/m]",
                "area ratio distribution"),
      IN(),
      GIN("Dmax", "t"),
      GIN_TYPE("Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Maximum dimension of the particles [m]",
               "Ambient atmospheric temperature [K]")));

  md_data_raw.push_back(MdRecord(
      NAME("dNdD_H98"),
      DESCRIPTION(
          "Calculation of particle size shape distribution (dN/dD) following\n"
          "continental stratus case of Hess et al. (1998).\n"
          "\n"
          "A wrapper to internal particle size distribution calculation. The\n"
          "distribution implemented here is for cloud liquid water,\n"
          "specifically for continental stratus. The parametrization is over\n"
          "diameter of spherical droplets.\n"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("dNdD"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("size distribution number density [#/m3/m]"),
      IN(),
      GIN("diameter_volume_equivalent", "LWC"),
      GIN_TYPE("Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Volume equivalent sphere diameter of the particles [m]",
               "Atmospheric liquid water content [kg/m3]")));

  md_data_raw.push_back(MdRecord(
      NAME("dNdD_MGD_IWC"),
      DESCRIPTION(
          "Calculation of particle size distribution (dN/dD) according\n"
          "to the modified gamma distribution for cloud ice\n"
          "inside of Geer and Baordo (2014).\n"
          "\n"
          "A wrapper to internal particle size distribution calculation.\n"
          "MDG_IWC is a parametrization for cloud ice water. It is a\n"
          "modified gamma distribution with the coefficients of\n"
          "Geer and Baordo (2014). It assumes spherical particles of constant\n"
          "density.\n"),
      AUTHORS("Manfred Brath"),
      OUT(),
      GOUT("dNdD"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("size distribution number density [#/m3/m]"),
      IN(),
      GIN("diameter_volume_equ", "rho", "IWC"),
      GIN_TYPE("Vector", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Volume equivalent diameter of the particles [m]",
               "Density of the particles [kg/m^3]",
               "Atmospheric ice water content [kg/m3]")));

  md_data_raw.push_back(MdRecord(
      NAME("dNdD_MGD_LWC"),
      DESCRIPTION(
          "Calculation of particle size distribution (dN/dD) according\n"
          "to the fixed modified gamma distribution for liquid cloud water\n"
          "inside of Geer and Baordo (2014).\n"
          "\n"
          "A wrapper to internal particle size distribution calculation.\n"
          "MDG_LWC is a parametrization for cloud liquid water. It is a\n"
          "modified gamma distribution with the coefficients of\n"
          " Geer and Baordo (2014). It assumes spherical particles of constant\n"
          "density.\n"),
      AUTHORS("Manfred Brath"),
      OUT(),
      GOUT("dNdD"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("size distribution number density [#/m3/m]"),
      IN(),
      GIN("diameter_volume_equ", "rho", "LWC"),
      GIN_TYPE("Vector", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Volume equivalent diameter of the particles [m]",
               "Density of the particles [kg/m^3]",
               "Atmospheric ice water content [kg/m3]")

          ));

  md_data_raw.push_back(MdRecord(
      NAME("dNdD_MH97"),
      DESCRIPTION(
          "Calculation of particle size distribution (dN/dD) following\n"
          "McFarquahar and Heymsfield (1997) parametrization.\n"
          "\n"
          "A wrapper to the internal particle size distribution calculation.\n"
          "McFarquhar and Heymsfield (1997) is a parametrization for cloud\n"
          "ice in the tropics. Parametrization is in ice water content (IWC)\n"
          "and ambient atmospheric temperature over particle size in terms of\n"
          "mass equivalent sphere diameter of the ice particles. McFarquhar\n"
          "and Heymsfield (1997) additionally provide uncertainties of the\n"
          "distribution's parameters, which can be used here to create\n"
          "perturbed distributions (set *noisy* to 1).\n"
          "\n"
          "Negative IWC trigger an error (unless robust=1, where IWC=0 is used\n"
          "internally, hence dNdD=0 is returned).\n"
          "Negative temperatures always trigger an error, temperatures >280K\n"
          "are only accepted if robust=1.\n"
          "For temperatures >273.15K (=0C), the distribution is evaluated\n"
          "assuming T=273.15K.\n"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("dNdD"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("size distribution number density [#/m3/m]"),
      IN(),
      GIN("diameter_mass_equivalent", "IWC", "T", "noisy", "robust"),
      GIN_TYPE("Vector", "Numeric", "Numeric", "Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "0", "0"),
      GIN_DESC("Mass equivalent sphere diameter of the particles [m].",
               "Atmospheric ice water content [kg/m3].",
               "Ambient atmospheric temperature [K].",
               "Distribution parameter perturbance flag.",
               "Flag whether to ignore parametrization value checks.")));

  md_data_raw.push_back(MdRecord(
      NAME("dNdD_MP48"),
      DESCRIPTION(
          "Calculation of particle size distribution (dN/dD) following\n"
          "Marshall and Palmer (1948) parametrization.\n"
          "\n"
          "A wrapper to internal particle size distribution calculation.\n"
          "Marshall and Palmer (1948) is a parametrization for precipitating\n"
          "hydrometeors, e.g., rain and snow. Parametrization is in\n"
          "precipitation rate (PR) over particle size, here taken in terms of\n"
          "melted equivalent sphere diameter.\n"
          "\n"
          "Precipitation rate can be given in different units, indicated by\n"
          "*PRunit*. Allowed are: 'mm/h' and 'kg/m2/s' and 'SI', where the\n"
          "latter two are equivalent. The latter two also require\n"
          "specification of *density* assuming a constant or mean density over\n"
          "all particle sizes.\n"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("dNdD"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("size distribution number density [#/m3/m]"),
      IN(),
      GIN("diameter_melted_equivalent", "PR", "PRunit", "density"),
      GIN_TYPE("Vector", "Numeric", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, "SI", NODEF),
      GIN_DESC("Melted equivalent sphere diameter of the particles [m]",
               "Precipitation rate [mm/h or kg/m2/s]",
               "Precipitation rate unit",
               "Material density [kg/m3] (required if PRunit!='mm/h')")));

  md_data_raw.push_back(MdRecord(
      NAME("dNdD_SB06"),
      DESCRIPTION(
          "Calculation of particle size distribution (dN/dm) following\n"
          "the two moment scheme of Seifert and Beheng, 2006b,a that is\n"
          "used in the ICON model.\n"
          "\n"
          "A wrapper to internal particle size distribution calculation. The Seifert two\n"
          "moment scheme is a parametrization for 6 different hydrometeors in total number\n"
          "density  (zeroth moment) and the mass concentration (first moment) as function\n"
          "of the  particle mass.\n"
          "Important, depending on the hydrometeor type, there is a lower and a upper\n"
          "boundary for valid particle mass. For masses below and above these boundary\n"
          "the number density is set to zero. Within this function the mass is NOT\n"
          "conserved.\n"
          "The 6 hydrometeors are:\n"
          "'SB06_LWC' - cloud liquid water\n"
          "'SB06_IWC' - cloud ice water\n"
          "'SB06_RWC' - rain\n"
          "'SB06_SWC' - snow\n"
          "'SB06_GWC' - graupel\n"
          "'SB06_HWC' - hail\n"),
      AUTHORS("Manfred Brath"),
      OUT(),
      GOUT("dNdD"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("size distribution number density [#/m3/m]"),
      IN(),
      GIN("mass", "N_tot", "M", "psd_type"),
      GIN_TYPE("Vector", "Numeric", "Numeric", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "SB06_LWC"),
      GIN_DESC("Mass of each particle [kg]",
               "Total number density [1/m3]",
               "Mass concentration [kg/m^3]",
               "Type of particle size distribution (hydrometeor type)")));

  md_data_raw.push_back(MdRecord(
      NAME("dNdD_SB06_M"),
      DESCRIPTION(
          "Calculation of particle size distribution (dN/dm) following\n"
          "the two moment scheme of Seifert and Beheng, 2006b,a that is\n"
          "used in the ICON model.\n"
          "A wrapper to internal particle size distribution calculation. The Seifert two\n"
          "moment scheme is originally a parametrization for 6 different hydrometeors in total number\n"
          "density  (zeroth moment) and the mass concentration (first moment) as function\n"
          "of the  particle mass. Here it is changed to a parametrization in \n"
          "mean particle mass and the mass concentration\n"
          "Important, depending on the hydrometeor type, there is a lower and a upper\n"
          "boundary for valid particle mass. For masses below and above these boundary\n"
          "the number density is set to zero. Within this function the mass is NOT\n"
          "conserved.\n"
          "The 6 hydrometeors are:\n"
          "'SB06_LWC' - cloud liquid water\n"
          "'SB06_IWC' - cloud ice water\n"
          "'SB06_RWC' - rain\n"
          "'SB06_SWC' - snow\n"
          "'SB06_GWC' - graupel\n"
          "'SB06_HWC' - hail\n"),
      AUTHORS("Manfred Brath"),
      OUT(),
      GOUT("dNdD"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("size distribution number density [#/m3/m]"),
      IN(),
      GIN("mass", "mean_mass", "M", "psd_type"),
      GIN_TYPE("Vector", "Numeric", "Numeric", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "SB06_LWC"),
      GIN_DESC("Mass of each particle [kg]",
               "Mean particle mass [kg]",
               "Mass concentration [kg/m^3]",
               "Type of particle size distribution (hydrometeor type)")));

  md_data_raw.push_back(MdRecord(
      NAME("dNdD_MY05"),
      DESCRIPTION(
          "Calculation of particle size distribution (dN/dm) following\n"
          "the Milbrandt and Yau two moment scheme, which is used in the GEM model.\n"
          "See also milbrandt and yau, 2005.\n"
          "\n"
          "A wrapper to internal particle size distribution calculation. The Milbrandt-Yau two\n"
          "moment scheme is a parametrization for 6 different hydrometeors in total number\n"
          "density  (zeroth moment) and the mass concentration (first moment) as function\n"
          "of the  mximum diameter.\n"
          "The 6 hydrometeors are:\n"
          "'MY05_LWC' - cloud liquid water\n"
          "'MY05_IWC' - cloud ice water\n"
          "'MY05_RWC' - rain\n"
          "'MY05_SWC' - snow\n"
          "'MY05_GWC' - graupel\n"
          "'MY05_HWC' - hail\n"),
      AUTHORS("Manfred Brath"),
      OUT(),
      GOUT("dNdD"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("size distribution number density [#/m3/m]"),
      IN(),
      GIN("diameter_max", "N_tot", "M", "psd_type"),
      GIN_TYPE("Vector", "Numeric", "Numeric", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "MY05_LWC"),
      GIN_DESC("Maximum diameter of each particle [m]",
               "Total number density [1/m3]",
               "Mass concentration [kg/m^3]",
               "Type of particle size distribution (hydrometeor type)")));

  md_data_raw.push_back(MdRecord(
      NAME("dNdD_MY05_M"),
      DESCRIPTION(
          "Calculation of particle size distribution (dN/dm) following\n"
          "the Milbrandt and Yau two moment scheme, which is used in the GEM model.\n"
          "See also milbrandt and yau, 2005.\n"
          "\n"
          "A wrapper to internal particle size distribution calculation. The Milbrandt-Yau two\n"
          "moment scheme is a parametrization for 6 different hydrometeors in total number\n"
          "density  (zeroth moment) and the mass concentration (first moment) as function\n"
          "of the  mximum diameter.\n"
          "The 6 hydrometeors are:\n"
          "'MY05_LWC' - cloud liquid water\n"
          "'MY05_IWC' - cloud ice water\n"
          "'MY05_RWC' - rain\n"
          "'MY05_SWC' - snow\n"
          "'MY05_GWC' - graupel\n"
          "'MY05_HWC' - hail\n"),
      AUTHORS("Manfred Brath"),
      OUT(),
      GOUT("dNdD"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("size distribution number density [#/m3/m]"),
      IN(),
      GIN("diameter_max", "N_tot", "M", "psd_type"),
      GIN_TYPE("Vector", "Numeric", "Numeric", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "MY05_LWC"),
      GIN_DESC("Maximum diameter of each particle [m]",
               "Total number density [1/m3]",
               "Mass concentration [kg/m^3]",
               "Type of particle size distribution (hydrometeor type)")));

  md_data_raw.push_back(MdRecord(
      NAME("dNdD_W16"),
      DESCRIPTION(
          "Calculation of raindrop size distribution (dN/dD) following\n"
          "Wang et al. (2016) parametrization.\n"
          "\n"
          "A wrapper to the internal particle size distribution calculation.\n"
          "Wang et al. (2016) is a parametrization for stratiform rain.\n"
          "Parametrization is in rain water mass content (RWC) over particle\n"
          "size in terms of mass equivalent sphere diameter of the liquid\n"
          "drops.\n"
          "\n"
          "Negative RWC trigger an error (unless robust=1, where RWC=0 is used\n"
          "internally, hence dNdD=0 is returned).\n"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("dNdD"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("size distribution number density [#/m3/m]"),
      IN(),
      GIN("diameter_mass_equivalent", "RWC", "robust"),
      GIN_TYPE("Vector", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, "0"),
      GIN_DESC("Mass equivalent sphere diameter of the particles [m].",
               "Rain water mass content [kg/m3].",
               "Flag whether to ignore parametrization value checks.")));

  md_data_raw.push_back(MdRecord(
      NAME("DOBatchCalc"),
      DESCRIPTION(
          "Performs batch calculations for radiation fields.\n"
          "\n"
          "We perform *ybatch_n* jobs, starting at index *ybatch_start*. (Zero\n"
          "based indexing, as usual.) The output arrays will have\n"
          "ybatch_n elements. Indices in the output array start\n"
          "with zero, independent of *ybatch_start*.\n"
          "\n"
          "WARNING, MEMORY INTENSIVE!!!: Since the outputs of this methods can\n"
          "be very large, make sure you only pass back output you need.\n"
          "Estimate the size of your output by looking at the dimensions\n"
          "beforehand. If you only want to pass back some fields, make sure to\n"
          "empty the others at the end of your *dobatch_calc_agenda*. E.g.:\n"
          "Tensor7SetConstant(doit_i_field, 0, 0, 0, 0, 0, 0, 0, 0.)\n"
          "\n"
          "The method performs the following:\n"
          "   1. Sets *ybatch_index* = *ybatch_start*.\n"
          "   2. Performs a-d until\n"
          "      *ybatch_index* = *ybatch_start* + *ybatch_n*.\n"
          "        a. Executes *dobatch_calc_agenda*.\n"
          "        b. If *ybatch_index* = *ybatch_start*, resizes the output\n"
          "           arrays based on *ybatch_n*.\n"
          "        c. Copies calculated fields to *ybatch_index* - *ybatch_start*\n"
          "           of output arrays.\n"
          "        d. Adds 1 to *ybatch_index*.\n"
          "\n"
          "Beside the *dobatch_calc_agenda*, the WSVs *ybatch_start*\n"
          "and *ybatch_n* must be set before calling this method.\n"
          "\n"
          "The input variable *ybatch_start* is set to a default of zero in\n"
          "*general.arts*.\n"),
      AUTHORS("Oliver Lemke"),
      OUT("dobatch_doit_i_field",
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
      GIN_DESC("A flag with value 1 or 0. If set to one, the batch\n"
               "calculation will continue, even if individual jobs fail. In\n"
               "that case, a warning message is written to screen and file\n"
               "(out1 output stream), and the output array entry for the\n"
               "failed job in the output fields is left empty.")));

  md_data_raw.push_back(MdRecord(
      NAME("DOAngularGridsSet"),
      DESCRIPTION(
          "Sets the angular grids for Discrete Ordinate type scattering\n"
          "calculations.\n"
          "\n"
          "This method sets the angular grids for the Discrete Ordinate type\n"
          "scattering calculations (DOIT, DISORT). For down- und up-looking\n"
          "geometries it suffices to define *N_za_grid* (both solvers) and\n"
          "*N_aa_grid* (DOIT). From these numbers equally spaced grids are\n"
          "created and stored in the WSVs *scat_za_grid* and *scat_aa_grid*.\n"
          "\n"
          "For limb simulations it is important to use an optimized zenith\n"
          "angle grid with a very fine resolution around the horizon\n"
          "(za=90 degrees). Such a grid can be generated using\n"
          "*doit_za_grid_optCalc*. To be applied, the name of the file holding\n"
          "the optimized angle grid has to be given (*za_grid_opt_file*).\n"
          "\n"
          "When an optimized grid is present, the equidistant grid is used for\n"
          "the calculation of the scattering integrals, while the optimized\n"
          "grid is applied for the integration of the radiative transfer\n"
          "equation. Otherwise the equidistant grid is used throughout. For\n"
          "down-looking cases using the equidistant grid typically suffices\n"
          "and speeds up the calculations.\n"),
      AUTHORS("Claudia Emde"),
      OUT("doit_za_grid_size", "scat_aa_grid", "scat_za_grid"),
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

  md_data_raw.push_back(MdRecord(
      NAME("DoitCalc"),
      DESCRIPTION(
          "Main DOIT method.\n"
          "\n"
          "This method executes *doit_mono_agenda* for each frequency\n"
          "in *f_grid*. The output is the radiation field inside the cloudbox\n"
          "(*doit_i_field*).\n"),
      AUTHORS("Claudia Emde"),
      OUT("doit_i_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field",
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

  md_data_raw.push_back(MdRecord(
      NAME("DoitGetIncoming"),
      DESCRIPTION(
          "Calculates incoming radiation field of the cloudbox by repeated\n"
          "radiative transfer calculations.\n"
          "\n"
          "The method performs monochromatic pencil beam calculations for\n"
          "all grid positions on the cloudbox boundary, and all directions\n"
          "given by scattering angle grids (*scat_za/aa_grid*). Found radiances\n"
          "are stored in *doit_i_field* which can be used as boundary\n"
          "conditions when scattering inside the cloud box is solved by the\n"
          "*DoitCalc* method.\n"
          "\n"
          "Note that *doit_i_field* will always hold intensity in terms of\n"
          "radiances, regardless of the setting of *iy_unit* (unit conversion\n"
          "is done within *yCalc* or *iyCalc*, which will provide their output\n"
          "in terms of the specified *iy_unit*; no explicit unit conversion by\n"
          "the user necessary.).\n"),
      AUTHORS("Sreerekha T.R.", "Claudia Emde"),
      OUT("doit_i_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field",
         "atmfields_checked",
         "atmgeom_checked",
         "cloudbox_checked",
         "doit_is_initialized",
         "iy_main_agenda",
         "atmosphere_dim",
         "lat_grid",
         "lon_grid",
         "z_field",
         "nlte_field",
         "cloudbox_on",
         "cloudbox_limits",
         "f_grid",
         "stokes_dim",
         "scat_za_grid",
         "scat_aa_grid"),
      GIN("rigorous", "maxratio"),
      GIN_TYPE("Index", "Numeric"),
      GIN_DEFAULT("1", "100"),
      GIN_DESC(
          "Fail if incoming field is not safely interpolable.",
          "Maximum allowed ratio of two radiances regarded as interpolable.")));

  md_data_raw.push_back(MdRecord(
      NAME("DoitGetIncoming1DAtm"),
      DESCRIPTION(
          "As *DoitGetIncoming* but assumes clear sky part to be 1D."
          "\n"
          "The incoming field is calculated only for one position and azimuth\n"
          "angle for each cloud box boundary, and obtained values are used\n"
          "for all other postions and azimuth angles. This works if a 3D\n"
          "cloud box is put into an 1D background atmosphere.\n"
          "\n"
          "This method can only be used for 3D cases.\n"),
      AUTHORS("Sreerekha T.R.", "Claudia Emde"),
      OUT("doit_i_field", "cloudbox_on"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field",
         "atmfields_checked",
         "atmgeom_checked",
         "cloudbox_checked",
         "doit_is_initialized",
         "iy_main_agenda",
         "atmosphere_dim",
         "lat_grid",
         "lon_grid",
         "z_field",
         "nlte_field",
         "cloudbox_on",
         "cloudbox_limits",
         "f_grid",
         "stokes_dim",
         "scat_za_grid",
         "scat_aa_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("DoitInit"),
      DESCRIPTION(
          "Initialises variables for DOIT scattering calculations.\n"
          "\n"
          "Note that multi-dimensional output variables (Tensors, specifically)\n"
          "are NaN-initialized. That is, this methods needs to be called\n"
          "BEFORE other WSMs that provide input to *DoitCalc*, e.g. before\n"
          "*DoitGetIncoming*.\n"),
      AUTHORS("Claudia Emde"),
      OUT("doit_scat_field", "doit_i_field", "doit_is_initialized"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("stokes_dim",
         "atmosphere_dim",
         "f_grid",
         "scat_za_grid",
         "scat_aa_grid",
         "doit_za_grid_size",
         "cloudbox_on",
         "cloudbox_limits"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("DoitScatteringDataPrepare"),
      DESCRIPTION(
          "Prepares single scattering data for a DOIT scattering calculation.\n"
          "\n"
          "First the scattering data is interpolated in frequency using\n"
          "*scat_data_monoCalc*. Then the phase matrix data is\n"
          "transformed or interpolated from the raw data to the laboratory frame\n"
          "for all possible combinations of the angles contained in the angular\n"
          "grids which are set in *DOAngularGridsSet*. The resulting phase\n"
          "matrices are stored in *pha_mat_sptDOITOpt*.\n"),
      AUTHORS("Claudia Emde"),
      OUT("pha_mat_sptDOITOpt",
          "scat_data_mono",
          "pha_mat_doit",
          "scat_aa_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_za_grid_size",
         "scat_aa_grid",
         "scat_data",
         "scat_data_checked",
         "f_index",
         "atmosphere_dim",
         "stokes_dim",
         "t_field",
         "cloudbox_limits",
         "pnd_field",
         "pha_mat_spt_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("DoitWriteIterationFields"),
      DESCRIPTION(
          "Writes DOIT iteration fields.\n"
          "\n"
          "This method writes intermediate iteration fields to xml-files. The\n"
          "method can be used as a part of *doit_conv_test_agenda*.\n"
          "\n"
          "The iterations to be stored are specified by *iterations*, e.g.:\n"
          "    iterations = [3, 6, 9]\n"
          "In this case the 3rd, 6th and 9th iterations are stored.\n"
          "If a number is larger than the total number of iterations, this\n"
          "number is ignored. If all iterations should be stored set\n"
          "   iterations = [-1]\n"
          "\n"
          "The frequencies to be stored are specified by *frequencies* in the\n"
          "same way as the iterations. The frequency index corresponds to the\n"
          "order of frequencies in *f_grid*.\n"
          "\n"
          "The output files are named doit_iteration_fX_iY.xml with X being the\n"
          "frequency index and iY the iteration counter.\n"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_iteration_counter", "doit_i_field_mono", "f_index"),
      GIN("iterations", "frequencies"),
      GIN_TYPE("ArrayOfIndex", "ArrayOfIndex"),
      GIN_DEFAULT("[-1]", "[-1]"),
      GIN_DESC("Selection of iterations to store.",
               "Selection of frequencies to store.")));

  md_data_raw.push_back(MdRecord(
      NAME("doit_conv_flagAbs"),
      DESCRIPTION(
          "DOIT convergence test (maximum absolute difference).\n"
          "\n"
          "The function calculates the absolute differences for two successive\n"
          "iteration fields. It picks out the maximum values for each Stokes\n"
          "component separately. The convergence test is fullfilled under the\n"
          "following conditions:\n"
          "   |I(m+1) - I(m)| < epsilon_1     Intensity.\n"
          "   |Q(m+1) - Q(m)| < epsilon_2     The other Stokes components.\n"
          "   |U(m+1) - U(m)| < epsilon_3   \n"
          "   |V(m+1) - V(m)| < epsilon_4   \n"
          "These conditions have to be valid for all positions in the\n"
          "cloudbox and for all directions.\n"),
      AUTHORS("Claudia Emde"),
      OUT("doit_conv_flag", "doit_iteration_counter", "doit_i_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_conv_flag",
         "doit_iteration_counter",
         "doit_i_field_mono",
         "doit_i_field_mono_old"),
      GIN("epsilon", "max_iterations", "nonconv_return_nan"),
      GIN_TYPE("Vector", "Index", "Index"),
      GIN_DEFAULT(NODEF, "100", "0"),
      GIN_DESC("Limits for convergence. A vector with length matching "
               "*stokes_dim* with unit [W / (m^2 Hz sr)].",
               "Maximum number of iterations allowed to reach convergence"
               "limit.",
               "Flag whether to accept result at max_iterations (0=default)"
               "or whether to return NaNs in case of non-convergence at"
               "max_iterations")));

  md_data_raw.push_back(MdRecord(
      NAME("doit_conv_flagAbsBT"),
      DESCRIPTION(
          "DOIT convergence test (maximum absolute difference in Rayleigh Jeans "
          "BT)\n"
          "\n"
          "As *doit_conv_flagAbs* but convergence limits are specified in\n"
          "Rayleigh-Jeans brighntess temperatures.\n"),
      AUTHORS("Sreerekha T.R.", "Claudia Emde"),
      OUT("doit_conv_flag", "doit_iteration_counter", "doit_i_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_conv_flag",
         "doit_iteration_counter",
         "doit_i_field_mono",
         "doit_i_field_mono_old",
         "f_grid",
         "f_index"),
      GIN("epsilon", "max_iterations", "nonconv_return_nan"),
      GIN_TYPE("Vector", "Index", "Index"),
      GIN_DEFAULT(NODEF, "100", "0"),
      GIN_DESC("Limits for convergence. A vector with length matching "
               "*stokes_dim* with unit [K].",
               "Maximum number of iterations allowed to reach convergence"
               "limit.",
               "Flag whether to accept result at max_iterations (0=default)"
               "or whether to return NaNs in case of non-convergence at"
               "max_iterations")));

  md_data_raw.push_back(MdRecord(
      NAME("doit_conv_flagLsq"),
      DESCRIPTION(
          "DOIT convergence test (least squares).\n"
          "\n"
          "As *doit_conv_flagAbsBT* but applies a least squares convergence\n"
          "test between two successive iteration fields.\n"
          "\n"
          "Warning: This method is not recommended because this kind of\n"
          "convergence test is not sufficiently strict, so that the\n"
          "DOIT result might be wrong.\n"),
      AUTHORS("Claudia Emde"),
      OUT("doit_conv_flag", "doit_iteration_counter", "doit_i_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_conv_flag",
         "doit_iteration_counter",
         "doit_i_field_mono",
         "doit_i_field_mono_old",
         "f_grid",
         "f_index"),
      GIN("epsilon", "max_iterations", "nonconv_return_nan"),
      GIN_TYPE("Vector", "Index", "Index"),
      GIN_DEFAULT(NODEF, "100", "0"),
      GIN_DESC("Limits for convergence. A vector with length matching "
               "*stokes_dim* with unit [K].",
               "Maximum number of iterations allowed to reach convergence"
               "limit.",
               "Flag whether to accept result at max_iterations (0=default)"
               "or whether to return NaNs in case of non-convergence at"
               "max_iterations")));

  md_data_raw.push_back(MdRecord(
      NAME("doit_i_field_monoIterate"),
      DESCRIPTION(
          "Iterative solution of the VRTE (DOIT method).\n"
          "\n"
          "A solution for the RTE with scattering is found using the\n"
          "DOIT method:\n"
          " 1. Calculate scattering integral using *doit_scat_field_agenda*.\n"
          " 2. Calculate RT with fixed scattered field using\n"
          "    *doit_rte_agenda*.\n"
          " 3. Convergence test using *doit_conv_test_agenda*.\n"
          "\n"
          "Note: The atmospheric dimensionality *atmosphere_dim* can be\n"
          "      either 1 or 3. To these dimensions the method adapts\n"
          "      automatically. 2D scattering calculations are not\n"
          "      supported.\n"),
      AUTHORS("Claudia Emde, Jakob Doerr"),
      OUT("doit_i_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field_mono",
         "doit_scat_field_agenda",
         "doit_rte_agenda",
         "doit_conv_test_agenda"),
      GIN("accelerated"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC(
          "Index wether to accelerate only the intensity (1) or the whole Stokes Vector (4)")));

  md_data_raw.push_back(MdRecord(
      NAME("doit_i_fieldClearskyPlaneParallel"),
      DESCRIPTION(
          "Clear-sky radiation field of a plane parallel atmosphere.\n"
          "\n"
          "The method assumes a 1D flat planet. Radiances along each direction\n"
          "given by *scat_za_grid* are calculated using *ppathPlaneParallel*\n"
          "and *iyEmissionStandard*.\n"
          "\n"
          "Surface properties are defined by *iy_surface_agenda*, i.e. there is no\n"
          "restriction e.g. specular surfaces. On the other hand, the method demands\n"
          "that the surface is placed exactly at the first pressure level.\n"
          "\n"
          "Note that the variable *ppath_lmax* is considered, and that it can be\n"
          "crucial for the accuracy for zenith angles close to 90 degrees. That\n"
          "is, using ppath_lmax=-1 is not recommended for this function.\n"
          "\n"
          "Information on transmission is also provided by the GOUT *trans_field*.\n"
          "For up-welling radiation (scat_za > 90), this variable holds the\n"
          "transmission to space, for considered position and propagation direction.\n"
          "For down-welling radiation, *trans_field* holds instead the transmission\n"
          "down to the surface.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("doit_i_field"),
      GOUT("trans_field"),
      GOUT_TYPE("Tensor3"),
      GOUT_DESC("Dimensions: [f_grid,p_grid,scat_za_grid]. See further above."),
      IN("propmat_clearsky_agenda",
         "water_p_eq_agenda",
         "iy_space_agenda",
         "iy_surface_agenda",
         "iy_cloudbox_agenda",
         "stokes_dim",
         "f_grid",
         "atmosphere_dim",
         "p_grid",
         "z_field",
         "t_field",
         "nlte_field",
         "vmr_field",
         "abs_species",
         "wind_u_field",
         "wind_v_field",
         "wind_w_field",
         "mag_u_field",
         "mag_v_field",
         "mag_w_field",
         "z_surface",
         "ppath_lmax",
         "rte_alonglos_v",
         "surface_props_data",
         "scat_za_grid"),
      GIN("use_parallel_iy"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("0: Parallelize over zenith angles\n"
               "1: Use more memory intensiv iyEmissionStandardParallel*")));

  md_data_raw.push_back(MdRecord(
      NAME("doit_i_fieldSetFromPrecalc"),
      DESCRIPTION(
          "Sets the initial cloudbox intensity field *doit_i_field* from a\n"
          "precalculated field.\n"
          "\n"
          "This method sets the (monochromatic) first guess radiation field\n"
          "inside the cloudbox from a precalculated *doit_i_field_precalc*,\n"
          "e.g., from the solution of a similar atmospheric scenario. The\n"
          "dimensions of *doit_i_field_precalc* have to be consistent with\n"
          "the DOIT setup in terms of frequencies, pressure levels inside the\n"
          "cloudbox, polar angles used as well as the stokes dimension.\n"
          "Incoming field on the cloudbox boundaries is adapted to the actual\n"
          "clearsky incoming field as, e.g., calculated by *DoitGetIncoming*.\n"),
      AUTHORS("Jana Mendrok"),
      OUT("doit_i_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field",
         "scat_za_grid",
         "f_grid",
         "atmosphere_dim",
         "stokes_dim",
         "cloudbox_limits",
         "doit_is_initialized"),
      GIN("doit_i_field_precalc"),
      GIN_TYPE("Tensor7"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Precalculated radiation field (of type *doit_i_field*)")));

  md_data_raw.push_back(MdRecord(
      NAME("doit_i_fieldSetClearsky"),
      DESCRIPTION(
          "Interpolate clearsky field on all gridpoints in cloudbox.\n"
          "\n"
          "This method uses a linear 1D/3D interpolation scheme to obtain the\n"
          "radiation field on all grid points inside the cloud box from the\n"
          "clear sky field on the cloudbox boundary. This radiation field\n"
          "is taken as the first guess radiation field in the DOIT module.\n"
          "\n"
          "Set the *all_frequencies* to 1 if the clearsky field shall be used\n"
          "as initial field for all frequencies. Set it to 0 if the clear sky\n"
          "field shall be used only for the first frequency in *f_grid*. For\n"
          "later frequencies, *doit_i_field* of the previous frequency is then\n"
          "used.\n"),
      AUTHORS("Sreerekha T.R. and Claudia Emde"),
      OUT("doit_i_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "cloudbox_limits",
         "atmosphere_dim",
         "cloudbox_on",
         "doit_is_initialized"),
      GIN("all_frequencies"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("See above.")));

  md_data_raw.push_back(MdRecord(
      NAME("doit_i_field_monoSetConst"),
      DESCRIPTION(
          "This method sets the initial field inside the cloudbox to a\n"
          "constant value. The method works only for monochromatic\n"
          "calculations (number of elements in f_grid=1).\n"
          "\n"
          "The user can specify a value for each Stokes dimension in the\n"
          "control file by *value*.\n"),
      AUTHORS("Claudia Emde"),
      OUT("doit_i_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field_mono",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "cloudbox_limits",
         "atmosphere_dim",
         "stokes_dim"),
      GIN("value"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("A vector containing 4 elements with the value of the "
               "initial field for each Stokes dimension.")));

  md_data_raw.push_back(MdRecord(
      NAME("doit_i_fieldSetConst"),
      DESCRIPTION(
          "This method sets the initial field inside the cloudbox to a\n"
          "constant value.\n"
          "\n"
          "The user has to specify a value for each Stokes dimension in the\n"
          "control file by *value*.\n"),
      AUTHORS("Claudia Emde"),
      OUT("doit_i_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "cloudbox_limits",
         "atmosphere_dim",
         "stokes_dim"),
      GIN("value"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("A vector containing *stokes_dim* elements with the value of"
               " the initial field for each Stokes dimension.")));

  md_data_raw.push_back(MdRecord(
      NAME("doit_i_fieldSetConstPerFreq"),
      DESCRIPTION(
          "This method sets the initial field inside the cloudbox to a\n"
          "constant value per frequency slice.\n"
          "\n"
          "The user has specify a value for each frequency and Stokes\n"
          "dimension in the control file by *value*.\n"),
      AUTHORS("Jana Mendrok"),
      OUT("doit_i_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "cloudbox_limits",
         "atmosphere_dim",
         "stokes_dim"),
      GIN("value"),
      GIN_TYPE("Matrix"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("A matrix containing *stokes_dim* elements per frequency"
               " (row) with the value of the initial field for each"
               " frequency and Stokes dimension.")));

  md_data_raw.push_back(MdRecord(
      NAME("doit_i_fieldUpdate1D"),
      DESCRIPTION(
          "RT calculation in cloudbox with fixed scattering integral (1D).\n"
          "\n"
          "Updates the radiation field (DOIT method). The method loops\n"
          "through the cloudbox to update the radiation field for all\n"
          "positions and directions in the 1D cloudbox.\n"
          "\n"
          "Note: This method is very inefficient, because the number of\n"
          "iterations scales with the number of cloudbox pressure levels.\n"
          "It is recommended to use *doit_i_fieldUpdateSeq1D*.\n"),
      AUTHORS("Claudia Emde"),
      OUT("doit_i_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field_mono",
         "doit_scat_field",
         "cloudbox_limits",
         "propmat_clearsky_agenda",
         "vmr_field",
         "spt_calc_agenda",
         "scat_za_grid",
         "pnd_field",
         "ppath_step_agenda",
         "ppath_lmax",
         "ppath_lraytrace",
         "p_grid",
         "z_field",
         "refellipsoid",
         "t_field",
         "f_grid",
         "f_index",
         "surface_rtprop_agenda",
         "doit_za_interp"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("doit_i_fieldUpdateSeq1D"),
      DESCRIPTION(
          "RT calculation in cloudbox with fixed scattering integral.\n"
          "\n"
          "Updates radiation field (*doit_i_field*) in DOIT module.\n"
          "This method loops through the cloudbox to update the\n"
          "radiation field for all positions and directions in the 1D\n"
          "cloudbox. The method applies the sequential update. For more\n"
          "information refer to AUG.\n"),
      AUTHORS("Claudia Emde"),
      OUT("doit_i_field_mono", "doit_scat_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field_mono",
         "doit_scat_field",
         "cloudbox_limits",
         "propmat_clearsky_agenda",
         "vmr_field",
         "spt_calc_agenda",
         "scat_za_grid",
         "scat_aa_grid",
         "pnd_field",
         "ppath_step_agenda",
         "ppath_lmax",
         "ppath_lraytrace",
         "p_grid",
         "z_field",
         "refellipsoid",
         "t_field",
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

  md_data_raw.push_back(MdRecord(
      NAME("doit_i_fieldUpdateSeq1DPP"),
      DESCRIPTION(
          "RT calculation in cloudbox with fixed scattering integral.\n"
          "\n "
          "Update radiation field (*doit_i_field*) in DOIT module.\n"
          "This method loops through the cloudbox to update the\n"
          "radiation field for all\n"
          "positions and directions in the 1D cloudbox. The method applies\n"
          "the sequential update and the plane parallel approximation.\n"
          "This method is only slightly faster than\n"
          "*doit_i_fieldUpdateSeq1D* and it is less accurate. It can not\n"
          "be used for limb simulations.\n"),
      AUTHORS("Sreerekha T.R."),
      OUT("doit_i_field_mono", "scat_za_index"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field_mono",
         "doit_scat_field",
         "cloudbox_limits",
         "propmat_clearsky_agenda",
         "vmr_field",
         "spt_calc_agenda",
         "scat_za_grid",
         "pnd_field",
         "p_grid",
         "z_field",
         "t_field",
         "f_grid",
         "f_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("doit_i_fieldUpdateSeq3D"),
      DESCRIPTION(
          "RT calculation in cloudbox with fixed scattering integral.\n"
          "\n"
          "Update radiation field (*doit_i_field*) in DOIT module.\n"
          "This method loops through the cloudbox to update the\n"
          "radiation field for all positions and directions in the 3D\n"
          "cloudbox. The method applies the sequential update. For more\n"
          "information please refer to AUG.\n"
          "Surface reflections are not yet implemented in 3D scattering\n"
          "calculations.\n"),
      AUTHORS("Claudia Emde"),
      OUT("doit_i_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field_mono",
         "doit_scat_field",
         "cloudbox_limits",
         "propmat_clearsky_agenda",
         "vmr_field",
         "spt_calc_agenda",
         "scat_za_grid",
         "scat_aa_grid",
         "pnd_field",
         "ppath_step_agenda",
         "ppath_lmax",
         "ppath_lraytrace",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "z_field",
         "refellipsoid",
         "t_field",
         "f_grid",
         "f_index",
         "doit_za_interp"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("doit_i_field_monoOptimizeReverse"),
      DESCRIPTION(
          "Interpolate *doit_i_field_mono* back to the original p_grid.\n"
          "For detailed description, see *OptimizeDoitPressureGrid*. \n"),
      AUTHORS("Jakob Doerr"),
      OUT("doit_i_field_mono"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field_mono", "p_grid_orig", "p_grid", "cloudbox_limits"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("OptimizeDoitPressureGrid"),
      DESCRIPTION(
          "Optimization of the pressure grid for RT calculation.\n"
          "The methods consists of three parts:\n"
          "1) Calculate the single scattering albedo and the scattering optical"
          "thickness from the scattering and absorption species. \n"
          "2) Enhance z_field according to the two thresholds sgl_alb_max and tau_scat_max."
          "If the resulting cloudbox size is bigger than cloudbox_size_max, this step is \n"
          "repeated with a higher threshold of tau_scat_max. \n"
          "3) Interpolate all variables used in doit_mono_agenda to the new z_field \n"
          "This method should be called inside\n"
          "*doit_mono_agenda*, right before *doit_i_field_monoIterate*. It can \n"
          "only be used if *ScatSpeciesMerge* has been called and if it is\n"
          "called, *doit_i_field_monoOptimizeReverse* has to be\n"
          "called right after *doit_i_field_monoIterate* to interpolate\n"
          "*doit_i_field_mono* back to the original size.\n"
          "Optimization currently only works with *stokes_dim* = 1 .\n"),
      AUTHORS("Jakob Doerr"),
      OUT("p_grid",
          "pnd_field",
          "t_field",
          "scat_data_mono",
          "z_field",
          "cloudbox_limits",
          "doit_i_field_mono",
          "pha_mat_doit",
          "vmr_field",
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
         "doit_i_field_mono",
         "pha_mat_doit",
         "vmr_field",
         "f_grid",
         "f_index",
         "propmat_clearsky_agenda"),
      GIN("tau_scat_max", "sgl_alb_max", "cloudbox_size_max"),
      GIN_TYPE("Numeric", "Numeric", "Index"),
      GIN_DEFAULT("0.1", "0.9", "200"),
      GIN_DESC("Maximum scattering optical thickness",
               "Maximum single scattering albedo",
               "Maximum cloudbox size")));

  md_data_raw.push_back(MdRecord(
      NAME("doit_scat_fieldCalc"),
      DESCRIPTION(
          "Calculates the scattering integral field in the DOIT module.\n"
          "\n"
          "The scattering integral field is generated by integrating\n"
          "the product of phase matrix and Stokes vector over all incident\n"
          "angles. For more information please refer to AUG.\n"),
      AUTHORS("Sreerekha T.R.", "Claudia Emde"),
      OUT("doit_scat_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_scat_field",
         "pha_mat_spt_agenda",
         "doit_i_field_mono",
         "pnd_field",
         "t_field",
         "atmosphere_dim",
         "cloudbox_limits",
         "scat_za_grid",
         "scat_aa_grid",
         "doit_za_grid_size",
         "pha_mat_doit"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("doit_scat_fieldCalcLimb"),
      DESCRIPTION(
          "Calculates the scattering integral field in the DOIT module (limb).\n"
          "\n"
          "The scattering integral field is the field generated by integrating\n"
          "the product of phase matrix and the Stokes vector over all incident\n"
          "angles.\n"
          "\n"
          "For limb simulations it makes sense to use different\n"
          "zenith angle grids for the scattering integral part and the RT part,\n"
          "because the latter part requires a much finer resolution near\n"
          "90 degrees. Taking an optimized grid for the RT part and an equidistant\n"
          "grid for the scattering integral part saves very much CPU time.\n"
          "This method uses the equidistant za_grid defined in\n"
          "*DOAngularGridsSet* and it should always be used for limb\n"
          "simulations.\n"
          "\n"
          "For more information please refer to AUG.\n"),
      AUTHORS("Claudia Emde"),
      OUT("doit_scat_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_scat_field",
         "pha_mat_spt_agenda",
         "doit_i_field_mono",
         "pnd_field",
         "t_field",
         "atmosphere_dim",
         "cloudbox_limits",
         "scat_za_grid",
         "scat_aa_grid",
         "doit_za_grid_size",
         "doit_za_interp",
         "pha_mat_doit"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("doit_za_grid_optCalc"),
      DESCRIPTION(
          "Zenith angle grid optimization for scattering calculation.\n"
          "\n"
          "This method optimizes the zenith angle grid. As input it requires\n"
          "a radiation field (*doit_i_field*) which is calculated on a very\n"
          "fine zenith angle grid (*scat_za_grid*). Based on this field\n"
          "zenith angle grid points are selected, such that the maximum\n"
          "difference between the radiation field represented on the very\n"
          "fine zenith angle grid and the radiation field represented on the\n"
          "optimized grid (*doit_za_grid_opt*) is less than the accuracy\n"
          "(*acc*). Between the grid points the radiation field is interpolated\n"
          "linearly or polynomially depending on *doit_za_interp*.\n"
          "\n"
          "Note: The method works only for a 1D atmosphere and for one\n"
          "frequency.\n"),
      AUTHORS("Claudia Emde"),
      OUT("doit_za_grid_opt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field_mono", "scat_za_grid", "doit_za_interp"),
      GIN("acc"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Accuracy to achieve [%].")));

  md_data_raw.push_back(MdRecord(
      NAME("doit_za_interpSet"),
      DESCRIPTION(
          "Define interpolation method for zenith angle dimension.\n"
          "\n"
          "You can use this method to choose the interpolation method for\n"
          "interpolations in the zenith angle dimension.\n"),
      AUTHORS("Claudia Emde"),
      OUT("doit_za_interp"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim"),
      GIN("interp_method"),
      GIN_TYPE("String"),
      GIN_DEFAULT("linear"),
      GIN_DESC("Interpolation method (\"linear\" or \"polynomial\").")));

  md_data_raw.push_back(MdRecord(
      NAME("Error"),
      DESCRIPTION(
          "Issues an error and exits ARTS.\n"
          "\n"
          "This method can be placed in agendas that must be specified, but\n"
          "are expected not to be used for the particular case. An inclusion\n"
          "in *surface_rtprop_agenda* could look like:\n   "
          "Error{\"Surface interceptions of propagation path not expected.\"}\n"
          "\n"
          "Ignore and other dummy method calls must still be included.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("Exit"),
      DESCRIPTION(
          "Stops the execution and exits ARTS.\n"
          "\n"
          "This method is handy if you want to debug one of your control\n"
          "files. You can insert it anywhere in the control file. When\n"
          "it is reached, it will terminate the program.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("Extract"),
      DESCRIPTION(
          "Extracts an element from an array.\n"
          "\n"
          "Copies the element with the given Index from the input\n"
          "variable to the output variable.\n"
          "\n"
          "For a Tensor3 as an input, it copies the page with the given\n"
          "Index from the input Tensor3 variable to the output Matrix.\n"
          "\n"
          "In other words, the selection is always done on the first dimension.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("needle"),
      GOUT_TYPE("Index, ArrayOfIndex, Numeric, Vector,"
                "Matrix, Matrix,"
                "Tensor3, Tensor4, Tensor4,"
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

  md_data_raw.push_back(MdRecord(
      NAME("ExtractFromMetaSingleScatSpecies"),
      DESCRIPTION(
          "Extract (numeric) parameters from scat_meta of a single scattering\n"
          "species.\n"
          "\n"
          "...\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("ext_matAddGas"),
      DESCRIPTION(
          "Add gas absorption to all diagonal elements of extinction matrix.\n"
          "\n"
          "The task of this method is to sum up the gas absorption of the\n"
          "different gas species and add the result to the extinction matrix.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("FastemStandAlone"),
      DESCRIPTION(
          "Stand-alone usage of FASTEM.\n"
          "\n"
          "FASTEM is a parameterisation of the emissivity of water surfaces\n"
          "including the impact of waves, salinity and non-specular effects.\n"
          "This is more or less direct interface to FASTEM, but slightly\n"
          "adopted to fit with ARTS. The unit of *frequency* and salinity\n"
          "differ, and this version is \"vectorised\" in frequency.\n"
          "\n"
          "The output is four emissivity and reflectivity values for each\n"
          "frequency. These values are defined in Eq. 13 of  \"An Improved\n"
          "Fast Microwave Water Emissivity Model\" by Liu, Weng and English,\n"
          "I3TRGS, 2011. Note that emissivity and reflectivity do not add up\n"
          "to 1, which is the way FASTEM compensates for non-specular effects.\n"
          "\n"
          "There is an error if any frequency is above 250 GHz, or if the skin\n"
          "temperature is below 260 K. If the skin temperature is below 270 K,\n"
          "it is adjusted to 270 K.\n"
          "\n"
          "FASTEM returns unphysical values for propagation close to the\n"
          "horizon, here emissivity and reflectivity can be outside [0,1].\n"
          "If either emissivity or reflectivity is below/above 0/1, it is\n"
          "set to 0/1, and the other value is set to 1/0. That is, e+r=1\n"
          "is enforced. These problems start about 15 degrees from the horizon.\n"),
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
               "The transmission of the atmosphere, along the propagation "
               "path of the downwelling radiation. One value per frequency.",
               "The version of FASTEM to use.")));

  md_data_raw.push_back(MdRecord(
      NAME("FieldFromGriddedField"),
      DESCRIPTION("Extract the data from a GriddedField.\n"
                  "\n"
                  "A check is performed that the grids from the\n"
                  "GriddedField match *p_grid*, *lat_grid* and *lon_grid*.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix, Tensor3, Tensor4, Tensor4"),
      GOUT_DESC("Extracted field."),
      IN("p_grid", "lat_grid", "lon_grid"),
      GIN("in"),
      GIN_TYPE(
          "GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Raw input gridded field.")));

  md_data_raw.push_back(MdRecord(
      NAME("FlagOff"),
      DESCRIPTION("Sets an index variable that acts as an on/off flag to 0.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("FlagOn"),
      DESCRIPTION("Sets an index variable that acts as an on/off flag to 1.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("ForLoop"),
      DESCRIPTION(
          "A simple for-loop.\n"
          "\n"
          "This method is handy when you quickly want to test out a calculation\n"
          "with a set of different settings.\n"
          "\n"
          "It does a for-loop from start to stop in steps of step (who would\n"
          "have guessed that). For each iteration, the agenda *forloop_agenda* is\n"
          "executed. Inside the agenda, the variable *forloop_index* is available\n"
          "as index counter.\n"
          "\n"
          "There are no other inputs to *forloop_agenda*, and also no outputs. That\n"
          "means, if you want to get any results out of this loop, you have to\n"
          "save it to files (for example with *WriteXMLIndexed*), since\n"
          "variables used inside the agenda will only be local.\n"
          "\n"
          "Note that this kind of for loop is not parallel.\n"
          "\n"
          "The method is intended for simple testing, not as a replacement of\n"
          "*ybatchCalc*. However, it is compatible with *ybatchCalc*, in the sense\n"
          "that *ybatchCalc* may occur inside *forloop_agenda*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("FrequencyFromWavelength"),
      DESCRIPTION(
          "Convert from wavelength [m] to frequency [Hz].\n"
          "\n"
          "This is a generic method. It can take a single wavelength value or a wavelength vector as input.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("FrequencyFromCGSAngularWavenumber"),
      DESCRIPTION(
          "Convert from angular wavenumber [cm^-1] to frequency [Hz].\n"
          "\n"
          "This converts angular wavenumber (2*PI/wavelength) into frequency.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("FrequencyFromCGSKayserWavenumber"),
      DESCRIPTION(
          "Convert from Kayser wavenumber [cm^-1] to frequency [Hz].\n"
          "\n"
          "This converts Kayser wavenumber (1/wavelength) into frequency.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("f_gridFromabs_linesSet"),
      DESCRIPTION(
          "Sets *f_grid* to a frequency grid determined by the available lines.\n"
          "Can set either from full catalog or from single line entry.\n"
          "\n"
          "If line_nr is below zero, then all lines are used.  If above, then the line"
          "in order of appearance in *abs_lines* is used.  Warning: there can be a lot\n"
          "of lines in *abs_lines*.\n"
          "\n"
          "The frequency grid is decided from the line's frequency plus/minus half_width,\n"
          "and the number of frequency bins between these limits is determined by nr_f_per_line.\n"
          "\n"
          "Note that frequencies below 0 Hz are not allowed by test, and that overlaps of\n"
          "half widths (when all lines are selected) is also not allowed.\n"),
      AUTHORS("Richard Larsson"),
      OUT("f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("half_width", "nr_f_per_line", "line_nr"),
      GIN_TYPE("Numeric", "Index", "Index"),
      GIN_DEFAULT("1", "0", "-1"),
      GIN_DESC(
          "Half width in Hz",
          "Number of frequency bins from F0-half_width to F0+half_width",
          "Index of line (0 is first; if negative then all lines are used)")));

  md_data_raw.push_back(MdRecord(
      NAME("f_gridFromabs_lines_per_speciesSetFromSpeciesTag"),
      DESCRIPTION(
          "Sets *f_grid* to a frequency grid determined by the lines in a SpeciesTag.\n"
          "\n"
          "The frequency grid is decided from the line's frequency plus/minus half_width,\n"
          "and the number of frequency bins between these limits is determined by nr_f_per_line.\n"
          "\n"
          "Note that frequencies below 0 Hz are not allowed by test, and that overlaps of\n"
          "half widths is also not allowed.\n"),
      AUTHORS("Richard Larsson"),
      OUT("f_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species", "abs_species"),
      GIN("half_width", "nr_f_per_line", "species_tag"),
      GIN_TYPE("Numeric", "Index", "String"),
      GIN_DEFAULT("1", "0", ""),
      GIN_DESC("Half width in Hz",
               "Number of frequency bins from F0-half_width to F0+half_width",
               "Species tag as found in *abs_species*")));

  md_data_raw.push_back(MdRecord(
      NAME("f_gridFromGasAbsLookup"),
      DESCRIPTION(
          "Sets *f_grid* to the frequency grid of *abs_lookup*.\n"
          "\n"
          "Must be called between importing/creating raw absorption table and\n"
          "call of *abs_lookupAdapt*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("f_gridFromSensorAMSU"),
      DESCRIPTION(
          "Automatically calculate f_grid to match the sensor.\n"
          "\n"
          "This method is handy if you are simulating an AMSU-type instrument,\n"
          "consisting of a few discrete channels. The case that channels touch,\n"
          "as for MHS, is handled correctly. But the case that channels overlap\n"
          "is not (yet) handled and results in an error message.\n"
          "\n"
          "The method calculates *f_grid* to match the instrument, as given by\n"
          "the local oscillator frequencies *lo_multi*, the backend\n"
          "frequencies *f_backend_multi*, and the backend channel\n"
          "responses *backend_channel_response_multi*.\n"
          "\n"
          "You have to specify the desired spacing in the keyword *spacing*,\n"
          "which has a default value of 100 MHz. (The actual value is 0.1e9,\n"
          "since our unit is Hz.)\n"
          "\n"
          "The produced grid will not have exactly the requested spacing, but\n"
          "will not be coarser than requested. The algorithm starts with the band\n"
          "edges, then adds additional points until the spacing is at least as\n"
          "fine as requested.\n"
          "\n"
          "There is a similar method for HIRS-type instruments,\n"
          "see *f_gridFromSensorHIRS*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("f_gridFromSensorAMSUgeneric"),
      DESCRIPTION(
          "Automatcially calculate f_grid to match the sensor. \n"
          "This function is based on 'f_gridFromSensorAMSU' \n"
          "\n"
          "The method calculates *f_grid* to match the instrument, as given by\n"
          "the backend frequencies *f_backend*, and the backend channel\n"
          "responses *backend_channel_response*.\n"
          "\n"
          "You have to specify the desired spacing in the keyword *spacing*,\n"
          "which has a default value of 100 MHz. (The actual value is 0.1e9,\n"
          "since our unit is Hz.)"
          "\n"
          "The produced grid will not have exactly the requested spacing, but\n"
          "it will not be coarser than requested. The algorithm starts with the band\n"
          "edges, then adds additional points until the spacing is at least as\n"
          "fine as requested.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("f_gridFromSensorHIRS"),
      DESCRIPTION(
          "Automatically calculate f_grid to match the sensor.\n"
          "\n"
          "This method is handy if you are simulating a HIRS-type instrument,\n"
          "consisting of a few discrete channels.\n"
          "\n"
          "It calculates f_grid to match the instrument, as given by the nominal\n"
          "band frequencies *f_backend* and the spectral channel response\n"
          "functions given by *backend_channel_response*.\n"
          "\n"
          "You have to specify the desired spacing in the keyword *spacing*, which\n"
          "has a default value of 5e8 Hz.\n"
          "\n"
          "The produced grid will not have exactly the requested spacing, but\n"
          "will not be coarser than requested. The algorithm starts with the band\n"
          "edges, then adds additional points until the spacing is at least as\n"
          "fine as requested.\n"
          "\n"
          "There is a similar method for AMSU-type instruments, see\n"
          "*f_gridFromSensorAMSU*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("f_gridMetMM"),
      DESCRIPTION(
          "Sets *f_grid* and associated variables match MetMM settings.\n"
          "\n"
          "The method calculates *f_grid* to match the specifications of a\n"
          "*met_mm_backend* table and method arguments.\n"
          "\n"
          "You have to specify the desired spacing using the keyword *freq_spacing*.\n"
          "You can pass a *Vector* with one element to apply the same spacing to all\n"
          "channels or pass a spacing value for each channel separately.\n"
          "\n"
          "Optionally, *freq_number* can be set to specify the mininum number of\n"
          "frequencies per passband for each channel. The frequencies are placed\n"
          "equally spaced in each passband. The minimum spacing resulting from\n"
          "*freq_number* and *freq_spacing* will be used for the calculation. To\n"
          "explicitly use *freq_spacing* for a channel, *freq_number* can be set\n"
          "to -1 for this channel.\n"
          "\n"
          "The number of elements in *freq_number* can either be the number of\n"
          "channels or 1. If only one element is given, this number is used for\n"
          "all channels. If *freq_number* is 1 and *freq_spacing* is wider than\n"
          "the bandwidth of the channel, one frequency is placed in the middle of\n"
          "each passband.\n"
          "\n"
          "Frequencies that would be closer than *freq_merge_threshold* in the\n"
          "generated *f_grid* are merged together. This value should be left at\n"
          "the default value. This is only meant to compensate for numerical\n"
          "inaccuracies in the frequency calculation to merge frequency that are\n"
          "supposed to be identical.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("g0Earth"),
      DESCRIPTION(
          "Gravity at zero altitude on Earth.\n"
          "\n"
          "Sets *g0* for the given latitude using a standard parameterisation.\n"),
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

  md_data_raw.push_back(MdRecord(NAME("g0Io"),
                                 DESCRIPTION("Gravity at zero altitude on Io.\n"
                                             "\n"
                                             "Numeric from Wikipedia.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("g0Jupiter"),
      DESCRIPTION(
          "Gravity at zero altitude on Jupiter.\n"
          "\n"
          "Sets *g0*  to mean equatorial gravity on Jupiter. Value provided by\n"
          "MPS under ESA-planetary study (TN1).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("g0Mars"),
      DESCRIPTION(
          "Gravity at zero altitude on Mars.\n"
          "\n"
          "Sets *g0*  to mean equatorial gravity on Mars. Value provided by\n"
          "MPS under ESA-planetary study (TN1).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("g0Venus"),
      DESCRIPTION(
          "Gravity at zero altitude on Venus.\n"
          "\n"
          "Sets *g0*  to mean equatorial gravity on Venus. Value from Ahrens\n"
          "(1995), provided by MPS under ESA-planetary study (TN1).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("geo_posEndOfPpath"),
      DESCRIPTION("Sets geo-position based on *ppath*.\n"
                  "\n"
                  "The geo-position is set to the position of the last point\n"
                  "of the present propagation path. This will be the surface,\n"
                  "top-of-the atmosphere or cloudbox position, depending of\n"
                  "observation geometry and if the cloudbox is active.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("geo_posLowestAltitudeOfPpath"),
      DESCRIPTION(
          "Sets geo-position based on *ppath*.\n"
          "\n"
          "The geo-position is set to the position of the last point\n"
          "of the present propagation path having the lowest altitude.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("geo_posWherePpathPassesZref"),
      DESCRIPTION(
          "Sets geo-position based on *ppath*.\n"
          "\n"
          "The geo-position is set to the position where the propagation\n"
          "path passes the reference altitude. If this altitude is passes\n"
          "more than once, the passing closest to the sensor is selected.\n"
          "If the reference altitude is not passed at all, *geo*pos* is\n"
          "set to NaN.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("geo_pos"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath"),
      GIN("z_ref"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Reference altitude.")));

  md_data_raw.push_back(MdRecord(
      NAME("GetEnvironmentVariable"),
      DESCRIPTION(
          "Copy the contents of an environment variable to an ARTS String or Index.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("String, Index"),
      GOUT_DESC("Contents of environment variable."),
      IN(),
      GIN("in"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Name of environment variable.")));

  md_data_raw.push_back(
      MdRecord(NAME("GriddedFieldGetName"),
               DESCRIPTION("Get the name of a GriddedField.\n"
                           "\n"
                           "See *ArrayOfGriddedFieldGetNames*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("GriddedFieldLatLonExpand"),
      DESCRIPTION(
          "Expands the latitude and longitude grid of the GriddedField to\n"
          "[-90, 90] and [0,360], respectively. Expansion is only done in\n"
          "the dimension(s), where the grid size is 1.\n"
          "The values from the input data will be duplicated to accomodate\n"
          "for the larger size of the output field.\n"
          "gfield_raw_out and gfield_raw_in can be the same variable.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE(
          "GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3"),
      GOUT_DESC("Expanded gridded field."),
      IN(),
      GIN("in"),
      GIN_TYPE(
          "GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Raw input gridded field.")));

  md_data_raw.push_back(MdRecord(
      NAME("GriddedFieldLatLonRegrid"),
      DESCRIPTION(
          "Interpolates the input field along the latitude and longitude dimensions\n"
          "to *lat_true* and *lon_true*.\n"
          "\n"
          "If the input longitude grid is outside of *lon_true* it will be shifted\n"
          "left or right by 360. If it covers 360 degrees, a cyclic interpolation\n"
          "will be performed.\n"
          "in and out fields can be the same variable.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE(
          "GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3"),
      GOUT_DESC("Regridded gridded field."),
      IN("lat_true", "lon_true"),
      GIN("in", "interp_order"),
      GIN_TYPE(
          "GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3",
          "Index"),
      GIN_DEFAULT(NODEF, "1"),
      GIN_DESC("Raw input gridded field.", "Interpolation order.")));

  md_data_raw.push_back(MdRecord(
      NAME("GriddedFieldPRegrid"),
      DESCRIPTION(
          "Interpolates the input field along the pressure dimension to *p_grid*.\n"
          "\n"
          "If zero-padding is applied (zeropadding=1), pressures that are\n"
          "outside the *p_grid* are set to 0. This is thought, e.g., for VMR\n"
          "fields that outside the given pressure can safely be assumed to be\n"
          "zero.\n"
          "Note: Using zeropadding for altitude and temperature fields is\n"
          "strongly discouraged (it will work here, though, but likely trigger\n"
          "errors later on).\n"
          "Extrapolation is allowed within the common 0.5grid-step margin,\n"
          "but is overruled by zeropadding.\n"
          "in and out fields can be the same variable.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("GriddedField3, GriddedField4, ArrayOfGriddedField3"),
      GOUT_DESC("Regridded gridded field."),
      IN("p_grid"),
      GIN("in", "interp_order", "zeropadding"),
      GIN_TYPE("GriddedField3, GriddedField4, ArrayOfGriddedField3",
               "Index",
               "Index"),
      GIN_DEFAULT(NODEF, "1", "0"),
      GIN_DESC("Raw input gridded field.",
               "Interpolation order.",
               "Apply zero-padding.")));

  md_data_raw.push_back(MdRecord(
      NAME("GriddedFieldZToPRegrid"),
      DESCRIPTION(
          "Interpolates the input field along the vertical dimension to *p_grid*.\n"
          "\n"
          "This is done from z_field, and thus requires the atmosphere to be set \n"
          "beforehand.\n"
          "\n"
          "The latitude and longitude grid of the input field must match *lat_grid*\n"
          "and *lon_grid* for the method to work.\n"
          "\n"
          "BETA mode.\n"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("GriddedField3"),
      GOUT_DESC("Regridded output; Pressure-gridded field."),
      IN("p_grid", "lat_grid", "lon_grid", "z_field"),
      GIN("in", "interp_order", "zeropadding"),
      GIN_TYPE("GriddedField3", "Index", "Index"),
      GIN_DEFAULT(NODEF, "1", "0"),
      GIN_DESC("Raw input; Altitude-gridded field.",
               "Interpolation order.",
               "Apply zero-padding.")));
  md_data_raw.push_back(MdRecord(
      NAME("heating_ratesFromIrradiance"),
      DESCRIPTION(
          "Calculates heating rates. It assumes that the heating rates\n"
          "depend only on the vertical derivation of the net flux.\n"
          "The net flux is the sum of the irradiance field in upward \n"
          " direction and the irradiance field in downward direction\n"),
      AUTHORS("Manfred Brath"),
      OUT("heating_rates"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("p_grid", "irradiance_field", "specific_heat_capacity", "g0"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("Ignore"),
      DESCRIPTION(
          "Ignore a workspace variable.\n"
          "\n"
          "This method is handy for use in agendas in order to suppress warnings\n"
          "about unused input workspace variables. What it does is: Nothing!\n"
          "In other words, it just ignores the variable it is called on.\n"
          "\n"
          "This method can ignore any workspace variable you want.\n"
          "\n"
          "Usage example:\n"
          "\n"
          "AgendaSet(els_agenda){\n"
          "  Ignore(ls_sigma)\n"
          "  elsLorentz\n"
          "}\n"
          "\n"
          "Without Ignore you would get an error message, because 'els_agenda' is\n"
          "supposed to use the Doppler width 'ls_sigma', but the Lorentz lineshape\n"
          "'elsLorentz' does not need it.\n"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("in"),
      GIN_TYPE("Any"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Variable to be ignored."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(MdRecord(
      NAME("INCLUDE"),
      DESCRIPTION(
          "Includes the contents of another controlfile.\n"
          "\n"
          "The INCLUDE statement inserts the contents of the controlfile\n"
          "with the given name into the current controlfile.\n"
          "If the filename is given without path information, ARTS will\n"
          "first search for the file in all directories specified with the\n"
          "-I (see arts -h) commandline option and then in directories given\n"
          "in the environment variable ARTS_INCLUDE_PATH. In the environment\n"
          "variable multiple paths have to be separated by colons.\n"
          "\n"
          "Note that INCLUDE is not a workspace method and thus the\n"
          "syntax is different:\n"
          "\n"
          "Arts {\n"
          "  INCLUDE \"general.arts\"\n"
          "}\n"
          "\n"
          "Includes can also be nested. In the example above general.arts\n"
          "can contain further includes which will then be treated\n"
          "the same way.\n"
          "\n"
          "The idea behind this mechanism is that you can write common settings\n"
          "for a bunch of calculations into one file. Then, you can create\n"
          "several controlfiles which include the basic settings and tweak them\n"
          "for different cases. When you decide to make changes to your setup\n"
          "that should apply to all calculations, you only have to make a\n"
          "single change in the include file instead of modifying all your\n"
          "controlfiles.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("IndexAdd"),
      DESCRIPTION(
          "Adds a index and a value (out = in+value).\n"
          "\n"
          "The result can either be stored in the same or another index.\n"
          "(in and out can be the same variable, but not out and value)\n"),
      AUTHORS("Patrick Eriksson, Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Index"),
      GOUT_DESC("Output numeric."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Index.", "Value to add.")));

  md_data_raw.push_back(MdRecord(
      NAME("IndexNumberOfAtmosphericPoints"),
      DESCRIPTION(
        "Counts number of points in the atmosphere.\n"
        "\n"
        "For a 3D atmosphere the method sets *n* to:\n"
        "  p_grid.nelem()*lat_grid.nelem()*lon_grid.nelem()\n"
        "For 1D and 2D the same calculation is done, but ignoring dimensions\n"
        "not active.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("n"),
      GOUT_TYPE("Index"),
      GOUT_DESC("Variable to set with number of points."),
      IN("atmosphere_dim", "p_grid", "lat_grid", "lon_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("IndexSet"),
      DESCRIPTION("Sets an index workspace variable to the given value.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Index"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("value"),
      GIN_TYPE("Index"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Value."),
      SETMETHOD(true)));

  md_data_raw.push_back(MdRecord(
      NAME("IndexSetToLast"),
      DESCRIPTION(
          "Set an Index to point towards last position of array-type variables.\n"
          "\n"
          "This method works as nelemGet, but gives the index number of the last\n"
          "element (which equals nelem-1).\n"),
      AUTHORS("Patrick Eriksson", "Oliver Lemke"),
      OUT("nelem"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("v"),
      GIN_TYPE(ARRAY_GROUPS + ", Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The method is defined for these groups."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(
      MdRecord(NAME("IndexStepDown"),
               DESCRIPTION("Performas: out = in - 1\n"
                           "\n"
                           "Input and output can be same variable.\n"),
               AUTHORS("Patrick Eriksson"),
               OUT(),
               GOUT("out"),
               GOUT_TYPE("Index"),
               GOUT_DESC("Output index variable."),
               IN(),
               GIN("in"),
               GIN_TYPE("Index"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Input index variable.")));

  md_data_raw.push_back(
      MdRecord(NAME("IndexStepUp"),
               DESCRIPTION("Performas: out = in + 1\n"
                           "\n"
                           "Input and output can be same variable.\n"),
               AUTHORS("Patrick Eriksson"),
               OUT(),
               GOUT("out"),
               GOUT_TYPE("Index"),
               GOUT_DESC("Output index variable."),
               IN(),
               GIN("in"),
               GIN_TYPE("Index"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Input index variable.")));

  md_data_raw.push_back(MdRecord(
      NAME("InterpAtmFieldToPosition"),
      DESCRIPTION("Point interpolation of atmospheric fields.\n"
                  "\n"
                  "The default way to specify the position is by *rtp_pos*.\n"
                  "\n"
                  "Linear interpolation is applied.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Value obtained by the interpolation."),
      IN("atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "z_field",
         "rtp_pos"),
      GIN("field"),
      GIN_TYPE("Tensor3"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Field to interpolate.")));

  md_data_raw.push_back(MdRecord(
      NAME("InterpGriddedField2ToPosition"),
      DESCRIPTION(
          "Latitude and longitude interpolation of a GriddedField2.\n"
          "\n"
          "The default way to specify the position is by *rtp_pos*.\n"
          "\n"
          "The interpolation is done for the latitude and longitude in\n"
          "*rtp_pos*. The altitude in *rtp_pos* is completely ignored.\n"
          "Linear interpolation is applied.\n"
          "\n"
          "The input field (*gfield2*) is expected to have latitude and\n"
          "longitude as first and second dimension.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Value obtained by interpolation."),
      IN("atmosphere_dim", "lat_grid", "lat_true", "lon_true", "rtp_pos"),
      GIN("gfield2"),
      GIN_TYPE("GriddedField2"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Gridded field to interpolate.")));

  md_data_raw.push_back(MdRecord(
      NAME("InterpSurfaceFieldToPosition"),
      DESCRIPTION(
          "Point interpolation of surface fields.\n"
          "\n"
          "The default way to specify the position is by *rtp_pos*.\n"
          "\n"
          "Linear interpolation is applied.\n"
          "\n"
          "The interpolation is done for the latitude and longitude in\n"
          "*rtp_pos*, while the altitude in *rtp_pos* is not part of the\n"
          "calculations. However, it is checked that the altitude of *rtp_pos*\n"
          "is inside the range covered by *z_surface* with a 1 m margin, to\n"
          "give a warning when the specified position is not consistent with\n"
          "the surface altitudes.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Value obtained by interpolation."),
      IN("atmosphere_dim", "lat_grid", "lon_grid", "rtp_pos", "z_surface"),
      GIN("field"),
      GIN_TYPE("Matrix"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Field to interpolate.")));

  md_data_raw.push_back(MdRecord(
      NAME("irradiance_fieldFromRadiance"),
      DESCRIPTION(
          "Calculate the irradiance also known as flux density from the *radiance_field* .\n"
          "by integrating over the angular grids according to the grids set\n"
          "by *AngularGridsSetForFluxCalc* \n"
          "See *AngularGridsSetForFluxCalc to set \n"
          "*scat_za_grid, scat_aa_grid, and za_grid_weights*\n"),
      AUTHORS("Manfred Brath"),
      OUT("irradiance_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("radiance_field", "scat_za_grid", "scat_aa_grid", "za_grid_weights"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("isotopologue_ratiosInitFromBuiltin"),
      DESCRIPTION(
          "Initialize isotopologue ratios with default values from built-in\n"
          "species data.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("iyActiveSingleScat"),
      DESCRIPTION(
          "Simulation of radar/lidar, restricted to single scattering.\n"
          "\n"
          "The WSM treats e.g. radar measurements of cloud and precipitation,\n"
          "on the condition that multiple scattering can be ignored. Beside\n"
          "the direct backsacttering, the two-way attenuation by gases and\n"
          "particles is considered. Surface scattering is ignored.\n"
          "\n"
          "The method could potentially be used for lidars, but multiple\n"
          "scattering poses here a must stronger constrain for the range of\n"
          "applications.\n"
          "\n"
          "The method can be used with *iyCalc*, but not with *yCalc*. In the\n"
          "later case, use instead *yActive*.\n"
          "\n"
          "The method returns the backscattering for each point of *ppath*.\n"
          "Several frequencies can be treated in parallel. The size of *iy*\n"
          "is [ nf*np, stokes_dim ], where nf is the length of *f_grid* and\n"
          "np is the number of path points. The data are stored in blocks\n"
          "of [ np, stokes_dim ]. That is, all the results for the first\n"
          "frequency occupy the np first rows of *iy* etc.\n"
          "\n"
          "The polarisation state of the transmitted pulse is taken from\n"
          "*iy_transmitter_agenda*. If the radar transmits several polarisations\n"
          "at the same frequency, you need to handle this by using two frequencies\n"
          "in *f_grid*, but these can be almost identical.\n"
          "\n"
          "This method does not consider *iy_unit*. Unit changes are insted applied\n"
          "in *yActive. The output of this method matches the option \"1\".\n"
          "\n"
          "Transmission is handled in a slightly simplified manner for efficiency\n"
          "reasons. First of all, the transmission matrix is assumed to be the same\n"
          "in both directions between the sensor and the point of back-scattering.\n"
          "This should in general be true, but exceptions could exist. The extinction\n"
          "due to particles can also be scaled, which could be of interest when e.g.\n"
          "characterising inversions.\n"
          "\n"
          "Further, for Jacobian calculations the default is to assume that the\n"
          "transmission is unaffected by the retrieval quantities. This is done\n"
          "to save computational time, and should be a valid approximation for the\n"
          "single-scattering conditions. Set *trans_in_jacobian* to 1 to obtain\n"
          "the more accurate Jacobian.\n"
          "\n"
          "Some auxiliary radiative transfer quantities can be obtained. Auxiliary\n"
          "quantities are selected by *iy_aux_vars* and returned by *iy_aux*.\n"
          "Valid choices for auxiliary data are:\n"
          " \"Radiative background\": Index value flagging the radiative\n"
          "    background. The following coding is used: 0=space, 1=surface\n"
          "    and 2=cloudbox.\n"
          " \"Backscattering\": The unattenuated backscattering. That is, as\n"
          "    *iy* but with no attenuated applied. Here all columns are filled.\n"
          " \"Optical depth\": Scalar, total and two-way, optical depth between\n"
          "    sensor and each point of the propagation path. Calculated based on\n"
          "    the (1,1)-element of the transmission matrix (1-based indexing),\n"
          "    i.e. only fully valid for scalar RT.\n"
          " \"Particle extinction\": As \"Optical depth\", but only with particle\n"
          "    attenuation included. That is, gas absorption is ignored.\n"
          "If nothing else is stated, only the first column of *iy_aux* is filled,\n"
          "i.e. the column matching Stokes element I, while remaing columns are\n"
          "are filled with zeros.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy",
          "iy_aux",
          "diy_dx",
          "ppvar_p",
          "ppvar_t",
          "ppvar_nlte",
          "ppvar_vmr",
          "ppvar_wind",
          "ppvar_mag",
          "ppvar_pnd",
          "ppvar_f",
          "ppvar_trans_cumulat"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
         "stokes_dim",
         "f_grid",
         "atmosphere_dim",
         "p_grid",
         "t_field",
         "nlte_field",
         "vmr_field",
         "abs_species",
         "wind_u_field",
         "wind_v_field",
         "wind_w_field",
         "mag_u_field",
         "mag_v_field",
         "mag_w_field",
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
         "propmat_clearsky_agenda",
         "water_p_eq_agenda",
         "iy_transmitter_agenda",
         "iy_agenda_call1",
         "iy_transmission",
         "rte_alonglos_v"),
      GIN("trans_in_jacobian", "pext_scaling", "t_interp_order"),
      GIN_TYPE("Index", "Numeric", "Index"),
      GIN_DEFAULT("0", "1", "1"),
      GIN_DESC("Flag determining if change in transmission is considered"
               " in calculation of the Jacobian or not.",
               "Particle extinction is scaled with this value. A value"
               " inside [0,2]. Set it to 0 if you want to remove particle"
               " extinction totally.",
               "Interpolation order of temperature for scattering data (so"
               " far only applied in phase matrix, not in extinction and"
               " absorption.")));

  md_data_raw.push_back(MdRecord(
      NAME("iyActiveSingleScat2"),
      DESCRIPTION("TESTING\n"),
      AUTHORS("Patrick Eriksson", "Richard Larsson"),
      OUT("iy",
          "iy_aux",
          "diy_dx",
          "ppvar_p",
          "ppvar_t",
          "ppvar_nlte",
          "ppvar_vmr",
          "ppvar_wind",
          "ppvar_mag",
          "ppvar_pnd",
          "ppvar_f",
          "ppvar_trans_cumulat"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
         "stokes_dim",
         "f_grid",
         "atmosphere_dim",
         "p_grid",
         "t_field",
         "nlte_field",
         "vmr_field",
         "abs_species",
         "wind_u_field",
         "wind_v_field",
         "wind_w_field",
         "mag_u_field",
         "mag_v_field",
         "mag_w_field",
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
         "propmat_clearsky_agenda",
         "water_p_eq_agenda",
         "iy_transmitter_agenda",
         "iy_agenda_call1",
         "iy_transmission",
         "rte_alonglos_v"),
      GIN("trans_in_jacobian", "pext_scaling", "t_interp_order"),
      GIN_TYPE("Index", "Numeric", "Index"),
      GIN_DEFAULT("0", "1", "1"),
      GIN_DESC("Flag determining if change in transmission is considered"
               " in calculation of the Jacobian or not.",
               "Particle extinction is scaled with this value. A value"
               " inside [0,2]. Set it to 0 if you want to remove particle"
               " extinction totally.",
               "Interpolation order of temperature for scattering data (so"
               " far only applied in phase matrix, not in extinction and"
               " absorption.")));

  md_data_raw.push_back(MdRecord(
      NAME("iyApplyUnit"),
      DESCRIPTION(
          "Conversion of *iy* to other spectral units.\n"
          "\n"
          "The method allows a change of unit, as a post-processing step,\n"
          "ignoring the n2-law of radiance.\n"
          "\n"
          "The conversion made inside *iyEmissionStandard* is mimiced,\n"
          "see that method for constraints and selection of output units.\n"
          "Restricted to that the n2-law can be ignored. This assumption\n"
          "is valid if the sensor is placed in space, or if the refractive\n"
          "index only deviates slightly from unity.\n"
          "\n"
          "It is stressed that there is no automatic check that the method is\n"
          "applied correctly, it is up to the user to ensure that the input\n"
          "data are suitable for the conversion.\n"
          "\n"
          "Beside *iy*, these auxilary quantities are modified:\n"
          "    \"iy\", \"Error\" and \"Error (uncorrelated)\"\n"
          "\n"
          "Please note that *diy_dx* is not handled.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy", "iy_aux"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy", "iy_aux", "stokes_dim", "f_grid", "iy_aux_vars", "iy_unit"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(
      MdRecord

      (NAME("iyCalc"),
       DESCRIPTION(
           "A single monochromatic pencil beam calculation.\n"
           "\n"
           "Performs monochromatic radiative transfer calculations for the\n"
           "specified position (*rte_pos*) and line-of-sight (*rte_pos*).\n"
           "See *iy* and associated variables for format of output.\n"
           "\n"
           "Please note that Jacobian type calculations not are supported.\n"
           "For this use *yCalc*.\n"
           "\n"
           "No sensor characteristics are applied. These are most easily\n"
           "incorporated by using *yCalc*\n"),
       AUTHORS("Patrick Eriksson"),
       OUT("iy", "iy_aux", "ppath"),
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
          "nlte_field",
          "rte_pos",
          "rte_los",
          "rte_pos2",
          "iy_unit",
          "iy_main_agenda"),
       GIN(),
       GIN_TYPE(),
       GIN_DEFAULT(),
       GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("iyEmissionStandard"),
      DESCRIPTION(
          "Standard method for radiative transfer calculations with emission.\n"
          "\n"
          "Designed to be part of *iy_main_agenda*. That is, only valid\n"
          "outside the cloudbox (no scattering). For details se the user guide.\n"
          "\n"
          "The possible choices for *iy_unit* are\n"
          " \"1\"             : No conversion, i.e. [W/(m^2 Hz sr)] (radiance per\n"
          "                     frequency unit).\n"
          " \"RJBT\"          : Conversion to Rayleigh-Jean brightness\n"
          "                     temperature.\n"
          " \"PlanckBT\"      : Conversion to Planck brightness temperature.\n"
          " \"W/(m^2 m sr)\"  : Conversion to [W/(m^2 m sr)] (radiance per\n"
          "                     wavelength unit).\n"
          " \"W/(m^2 m-1 sr)\": Conversion to [W/(m^2 m-1 sr)] (radiance per\n"
          "                     wavenumber unit).\n"
          "Expressions applied and considerations for the unit conversion of\n"
          "radiances are discussed in Sec. 5.7 of the ARTS-2.0 article.\n"
          "\n"
          "*iy_unit* is only applied if *iy_agenda_call1* is 1. This means that\n"
          "no unit ocnversion is applied for internal iterative calls.\n"
          "\n"
          "Some auxiliary radiative transfer quantities can be obtained. Auxiliary\n"
          "quantities are selected by *iy_aux_vars* and returned by *iy_aux*.\n"
          "Valid choices for auxiliary data are:\n"
          " \"Radiative background\": Index value flagging the radiative\n"
          "    background. The following coding is used: 0=space, 1=surface\n"
          "    and 2=cloudbox.\n"
          " \"Optical depth\": Scalar optical depth between the observation point\n"
          "    and the end of the present propagation path. Calculated based on\n"
          "    the (1,1)-element of the transmission matrix (1-based indexing),\n"
          "    i.e. only fully valid for scalar RT.\n"
          "If nothing else is stated, only the first column of *iy_aux* is filled,\n"
          "i.e. the column matching Stokes element I, while remaing columns are\n"
          "are filled with zeros.\n"),
      AUTHORS("Patrick Eriksson", "Richard Larsson"),
      OUT("iy",
          "iy_aux",
          "diy_dx",
          "ppvar_p",
          "ppvar_t",
          "ppvar_nlte",
          "ppvar_vmr",
          "ppvar_wind",
          "ppvar_mag",
          "ppvar_f",
          "ppvar_iy",
          "ppvar_trans_cumulat",
          "ppvar_trans_partial"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
         "iy_id",
         "stokes_dim",
         "f_grid",
         "atmosphere_dim",
         "p_grid",
         "t_field",
         "nlte_field",
         "vmr_field",
         "abs_species",
         "wind_u_field",
         "wind_v_field",
         "wind_w_field",
         "mag_u_field",
         "mag_v_field",
         "mag_w_field",
         "cloudbox_on",
         "iy_unit",
         "iy_aux_vars",
         "jacobian_do",
         "jacobian_quantities",
         "ppath",
         "rte_pos2",
         "propmat_clearsky_agenda",
         "water_p_eq_agenda",
         "iy_main_agenda",
         "iy_space_agenda",
         "iy_surface_agenda",
         "iy_cloudbox_agenda",
         "iy_agenda_call1",
         "iy_transmission",
         "rte_alonglos_v",
         "surface_props_data"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("iyEmissionStandardParallel"),
      DESCRIPTION(
          "Parallelized version of *iyEmissionStandard* at the expense of\n"
          "higher memory usage.\n"
          "\n"
          "For documentation see *iyEmissionStandard*.\n"),
      AUTHORS("Patrick Eriksson", "Richard Larsson", "Oliver Lemke"),
      OUT("iy",
          "iy_aux",
          "diy_dx",
          "ppvar_p",
          "ppvar_t",
          "ppvar_nlte",
          "ppvar_vmr",
          "ppvar_wind",
          "ppvar_mag",
          "ppvar_f",
          "ppvar_iy",
          "ppvar_trans_cumulat",
          "ppvar_trans_partial"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
         "iy_id",
         "stokes_dim",
         "f_grid",
         "atmosphere_dim",
         "p_grid",
         "t_field",
         "nlte_field",
         "vmr_field",
         "abs_species",
         "wind_u_field",
         "wind_v_field",
         "wind_w_field",
         "mag_u_field",
         "mag_v_field",
         "mag_w_field",
         "cloudbox_on",
         "iy_unit",
         "iy_aux_vars",
         "jacobian_do",
         "jacobian_quantities",
         "ppath",
         "rte_pos2",
         "propmat_clearsky_agenda",
         "water_p_eq_agenda",
         "iy_main_agenda",
         "iy_space_agenda",
         "iy_surface_agenda",
         "iy_cloudbox_agenda",
         "iy_agenda_call1",
         "iy_transmission",
         "rte_alonglos_v",
         "surface_props_data"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));
  /*
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "iyFOS" ),
        DESCRIPTION
        (
         "Method in development. Don't use without contacting Patrick.\n"
         "\n"
         "Regarding radiance unit, works exactly as *iyEmissionStandard*.\n"
         "\n"
         "The *fos_n* argument determines the maximum scattering order that\n"
         "will be considered. For example, 1 corresponds to that only single\n"
         "scattering is considered. The value 0 is accepted and results\n"
         "in calculations of clear-sky type. In the later case, particle\n"
         "absorption/emission is considered if cloudbox is active. If\n"
         "cloudbox is not active,clear-sky results are returned for all\n"
         "values of *fos_n*.\n"
         "\n"
         "The following auxiliary data can be obtained:\n"
         "  \"Pressure\": The pressure along the propagation path.\n"
         "     Size: [1,1,1,np].\n"
         "  \"Temperature\": The temperature along the propagation path.\n"
         "     Size: [1,1,1,np].\n"
         "  \"VMR, species X\": VMR of the species with index X (zero based).\n"
         "     For example, adding the string \"VMR, species 0\" extracts the\n"
         "     VMR of the first species. Size: [1,1,1,np].\n"
         "  \"Absorption, summed\": The total absorption matrix along the\n"
         "     path. Size: [nf,ns,ns,np].\n"
         "  \"Absorption, species X\": The absorption matrix along the path\n"
         "     for an individual species (X works as for VMR).\n"
         "     Size: [nf,ns,ns,np].\n"
         "  \"PND, type X\": The particle number density for scattering element\n"
         "       type X (ie. corresponds to book X in pnd_field).\n"
         "       Size: [1,1,1,np].\n"
         "  \"Mass content, X\": The mass content for scattering element X.\n"
         "       This corresponds to column X in *particle_masses* (zero-\n"
         "       based indexing). Size: [1,1,1,np].\n"
         "* \"Radiative background\": Index value flagging the radiative\n"
         "     background. The following coding is used: 0=space and\n"
         "     and 1=surface. Size: [nf,1,1,1].\n"
         "  \"iy\": The radiance at each point along the path (*iy_unit* is.\n"
         "     considered). Size: [nf,ns,1,np].\n"
         "* \"Optical depth\": The scalar optical depth between the\n"
         "     observation point and the end of the primary propagation path\n"
         "     (ie. the optical depth to the surface or space.). Calculated\n"
         "     in a pure scalar manner, and not dependent on direction.\n"
         "     Size: [nf,1,1,1].\n"
         "where\n"
         "  nf: Number of frequencies.\n"
         "  ns: Number of Stokes elements.\n"
         "  np: Number of propagation path points.\n"
         "\n"
         "The auxiliary data are returned in *iy_aux* with quantities\n"
         "selected by *iy_aux_vars*. Most variables require that the method\n"
         "is called directly or by *iyCalc*. For calculations using *yCalc*,\n"
         "the selection is restricted to the variables marked with *.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "iy", "iy_aux", "ppath", "diy_dx" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "diy_dx", "stokes_dim", "f_grid", "atmosphere_dim",
            "p_grid", "z_field", "t_field", "vmr_field", "abs_species", 
            "wind_u_field", "wind_v_field", "wind_w_field", "mag_u_field",
            "mag_v_field", "mag_w_field", "cloudbox_on", "cloudbox_limits",
            "pnd_field", "scat_data",
            "particle_masses", "iy_unit", "iy_aux_vars", "jacobian_do", 
            "ppath_agenda", 
            "propmat_clearsky_agenda", "iy_main_agenda", "iy_space_agenda", 
            "iy_surface_agenda", "iy_agenda_call1", "iy_transmission", 
            "rte_pos", "rte_los", "rte_pos2", "rte_alonglos_v",
            "ppath_lmax", "ppath_lraytrace",
            "fos_scatint_angles", "fos_iyin_za_angles"
            ),
        GIN( "fos_za_interporder", "fos_n" ),
        GIN_TYPE( "Index", "Index" ),
        GIN_DEFAULT( "1", "1" ),
        GIN_DESC( "Polynomial order for zenith angle interpolation.",
                  "Max scattering order to consider." )
        ));
  */

  md_data_raw.push_back(MdRecord(
      NAME("iyHybrid"),
      DESCRIPTION("So far just for testing.\n"),
      AUTHORS("Patrick Eriksson", "Jana Mendrok", "Richard Larsson"),
      OUT("iy",
          "iy_aux",
          "diy_dx",
          "ppvar_p",
          "ppvar_t",
          "ppvar_nlte",
          "ppvar_vmr",
          "ppvar_wind",
          "ppvar_mag",
          "ppvar_pnd",
          "ppvar_f",
          "ppvar_iy",
          "ppvar_trans_cumulat"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
         "iy_id",
         "stokes_dim",
         "f_grid",
         "atmosphere_dim",
         "p_grid",
         "t_field",
         "nlte_field",
         "vmr_field",
         "abs_species",
         "wind_u_field",
         "wind_v_field",
         "wind_w_field",
         "mag_u_field",
         "mag_v_field",
         "mag_w_field",
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
         "iy_main_agenda",
         "iy_space_agenda",
         "iy_surface_agenda",
         "iy_cloudbox_agenda",
         "iy_agenda_call1",
         "iy_transmission",
         "ppath",
         "rte_pos2",
         "rte_alonglos_v",
         "surface_props_data",
         "doit_i_field",
         "scat_za_grid"),
      GIN("Naa_grid", "t_interp_order"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT("19", "1"),
      GIN_DESC("Number of azimuth angles to consider in scattering source term"
               " integral.",
               "Interpolation order of temperature for scattering data (so"
               " far only applied in phase matrix, not in extinction and"
               " absorption.")));

  md_data_raw.push_back(MdRecord(
      NAME("iyHybrid2"),
      DESCRIPTION("So far just for even more testing.\n"),
      AUTHORS("Patrick Eriksson", "Jana Mendrok", "Richard Larsson"),
      OUT("iy",
          "iy_aux",
          "diy_dx",
          "ppvar_p",
          "ppvar_t",
          "ppvar_nlte",
          "ppvar_vmr",
          "ppvar_wind",
          "ppvar_mag",
          "ppvar_pnd",
          "ppvar_f",
          "ppvar_iy",
          "ppvar_trans_cumulat"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
         "iy_id",
         "stokes_dim",
         "f_grid",
         "atmosphere_dim",
         "p_grid",
         "t_field",
         "nlte_field",
         "vmr_field",
         "abs_species",
         "wind_u_field",
         "wind_v_field",
         "wind_w_field",
         "mag_u_field",
         "mag_v_field",
         "mag_w_field",
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
         "iy_main_agenda",
         "iy_space_agenda",
         "iy_surface_agenda",
         "iy_cloudbox_agenda",
         "iy_agenda_call1",
         "iy_transmission",
         "ppath",
         "rte_pos2",
         "rte_alonglos_v",
         "surface_props_data",
         "doit_i_field",
         "scat_za_grid"),
      GIN("Naa_grid", "t_interp_order"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT("19", "1"),
      GIN_DESC("Number of azimuth angles to consider in scattering source term"
               " integral.",
               "Interpolation order of temperature for scattering data (so"
               " far only applied in phase matrix, not in extinction and"
               " absorption.")));

  md_data_raw.push_back(MdRecord(
      NAME("iyIndependentBeamApproximation"),
      DESCRIPTION(
         "In development ....\n"
         "\n"
         "Describe how *atm_fields_compact* is filled.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy", "iy_aux", "ppath", "diy_dx", "atm_fields_compact"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
         "iy_id",
         "f_grid",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "lat_true",
         "lon_true",
         "t_field",
         "z_field",
         "vmr_field",
         "nlte_field",
         "wind_u_field",
         "wind_v_field",
         "wind_w_field",
         "mag_u_field",
         "mag_v_field",
         "mag_w_field",
         "cloudbox_on",
         "cloudbox_limits",
         "pnd_field",
         "particle_masses",
         "ppath_agenda",
         "ppath_lmax",
         "ppath_lraytrace",
         "iy_agenda_call1",
         "iy_unit",
         "iy_transmission",
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

  md_data_raw.push_back(MdRecord(
      NAME("iyInterpCloudboxField"),
      DESCRIPTION(
          "Interpolates the intensity field of the cloud box.\n"
          "\n"
          "Revised method using polynomials of user-defined order for polar\n"
          "angle interpolation. It can be (and should be) restricted to one\n"
          "hemisphere (default is currently 'not restricted' in order to\n"
          "reflect previous behaviour of the WSM, which is now\n"
          "*iyInterpLinCloudboxField*).\n"
          "\n"
          "Spatial interpolation so far hardcoded as linear.\n"),
      AUTHORS("Jana Mendrok"),
      OUT("iy"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field",
         "rtp_pos",
         "rtp_los",
         "jacobian_do",
         "cloudbox_on",
         "cloudbox_limits",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "z_field",
         "stokes_dim",
         "scat_za_grid",
         "scat_aa_grid",
         "f_grid"),
      GIN("za_interp_order",
          "za_restrict",
          "cos_za_interp",
          "za_extpolfac",
          "aa_interp_order"),
      GIN_TYPE("Index", "Index", "Index", "Numeric", "Index"),
      GIN_DEFAULT("1", "0", "0", "1.0", "1"),
      GIN_DESC("Zenith angle interpolation order.",
               "Flag whether to restric zenith angle interpolation to one"
               " hemisphere.",
               "Flag whether to do zenith angle interpolation in cosine space.",
               "Maximum allowed extrapolation range in zenith angle.",
               "Azimuth angle interpolation order.")));

  md_data_raw.push_back(MdRecord(
      NAME("iyInterpLinCloudboxField"),
      DESCRIPTION(
          "Interpolates the intensity field of the cloud box.\n"
          "\n"
          "This is the standard method to put in *iy_cloudbox_agenda* if the\n"
          "the scattering inside the cloud box is handled by the DOIT method.\n"
          "\n"
          "The intensity field is interpolated to the position (specified by\n"
          "*rtp_pos*) and direction (specified by *rtp_los*) given. Linear\n"
          "interpolation is used for all dimensions.\n"
          "\n"
          "The intensity field on the cloux box boundaries is provided by\n"
          "*doit_i_field* and interpolated if the given position is at any boundary.\n"),
      AUTHORS("Claudia Emde"),
      OUT("iy"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field",
         "rtp_pos",
         "rtp_los",
         "jacobian_do",
         "cloudbox_on",
         "cloudbox_limits",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "z_field",
         "stokes_dim",
         "scat_za_grid",
         "scat_aa_grid",
         "f_grid"),
      GIN("rigorous", "maxratio"),
      GIN_TYPE("Index", "Numeric"),
      GIN_DEFAULT("1", "3"),
      GIN_DESC(
          "Fail if cloudbox field is not safely interpolable.",
          "Maximum allowed ratio of two radiances regarded as interpolable.")));

  md_data_raw.push_back(MdRecord(
      NAME("iyInterpPolyCloudboxField"),
      DESCRIPTION(
          "As *iyInterpCloudboxField* but performs quadratic interpolation.\n"
          "\n"
          "Works so far only for 1D cases, and accordingly a quadratic\n"
          "interpolation along *scat_za_grid* is performed.\n"),
      AUTHORS("Claudia Emde"),
      OUT("iy"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field",
         "rtp_pos",
         "rtp_los",
         "jacobian_do",
         "cloudbox_on",
         "cloudbox_limits",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "z_field",
         "stokes_dim",
         "scat_za_grid",
         "scat_aa_grid",
         "f_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("iyLoopFrequencies"),
      DESCRIPTION(
          "Radiative transfer calculations one frequency at the time.\n"
          "\n"
          "The method loops the frequencies in *f_grid* and calls\n"
          "*iy_loop_freqs_agenda* for each individual value. This method is\n"
          "placed in *iy_main_agenda*, and the actual radiative transfer\n"
          " method is put in *iy_loop_freqs_agenda*.\n"
          "\n"
          "A common justification for using the method should be to consider\n"
          "dispersion. By using this method it is ensured that the propagation\n"
          "path for each individual frequency is calculated.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy", "iy_aux", "ppath", "diy_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy_aux_vars",
         "iy_agenda_call1",
         "iy_transmission",
         "rte_pos",
         "rte_los",
         "rte_pos2",
         "stokes_dim",
         "f_grid",
         "iy_loop_freqs_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("iyMC"),
      DESCRIPTION(
          "Interface to Monte Carlo part for *iy_main_agenda*.\n"
          "\n"
          "Basically an interface to *MCGeneral* for doing monochromatic\n"
          "pencil beam calculations. This functions allows Monte Carlo (MC)\n"
          "calculations for sets of frequencies and sensor pos/los in a single\n"
          "run. Sensor responses can be included in the standard manner\n"
          "(through *yCalc*).\n"
          "\n"
          "This function does not apply the MC approach when it comes\n"
          "to sensor properties. These properties are not considered when\n"
          "tracking photons, which is done in *MCGeneral* (but then only for\n"
          "the antenna pattern).\n"
          "\n"
          "Output unit options  (*iy_unit*) exactly as for *MCGeneral*.\n"
          "\n"
          "The MC calculation errors are all assumed be uncorrelated and each\n"
          "have a normal distribution. These properties are of relevance when\n"
          "weighting the errors with the sensor repsonse matrix. The seed is\n"
          "reset for each call of *MCGeneral* to obtain uncorrelated errors.\n"
          "\n"
          "MC control arguments (mc_std_err, mc_max_time, mc_min_iter, mc_max_iter\n"
          "mc_taustep_limit) as for *MCGeneral*. The arguments are applied\n"
          "for each monochromatic pencil beam calculation individually.\n"
          "As for *MCGeneral*, the value of *mc_error* shall be adopted to\n"
          "*iy_unit*.\n"
          "\n"
          "The following auxiliary data can be obtained:\n"
          "  \"Error (uncorrelated)\": Calculation error. Size: [nf,ns,1,1].\n"
          "    (The later part of the text string is required. It is used as\n"
          "    a flag to yCalc for how to apply the sensor data.)\n"
          "where\n"
          "  nf: Number of frequencies.\n"
          "  ns: Number of Stokes elements.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy", "iy_aux", "diy_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy_agenda_call1",
         "iy_transmission",
         "rte_pos",
         "rte_los",
         "iy_aux_vars",
         "jacobian_do",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "z_field",
         "t_field",
         "vmr_field",
         "refellipsoid",
         "z_surface",
         "cloudbox_on",
         "cloudbox_limits",
         "stokes_dim",
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

  /*
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "iyRadioLink" ),
        DESCRIPTION
        (
         "Radiative transfer for (active) radio links.\n"
         "\n"
         "The method assumes that *ppath_agenda* is set up to return the\n"
         "propagation path between the transmitter and the receiver. The\n" 
         "position of the transmitter is given as *rte_pos*, and the\n"
         "\"sensor\" is taken as the receiver.\n"
         "\n"
         "The primary output (*y*) is the received signal, where the signal\n"
         "transmitted is taken from *iy_transmitter_agenda*. That is, *y*\n"
         "is a Stokes vector for each frequency considered. Several other\n"
         "possible measurements quantities, such as the bending angle, can\n"
         "be obtained as the auxiliary data (see lost below).\n"
         "\n"
         "If it is found that no link can be obtained due to intersection of\n"
         "the ground, all data are set to zero. If no link could be\n"
         "determined for other reasons (due to critical refraction or\n"
         "numerical problems), all data are set to NaN.\n"
         "\n"
         "This method is just intended for approximative calculations for\n"
         "cases corresponding to relatively simple ray tracing. A detailed,\n"
         "and more exact, treatment of several effects require more advanced\n"
         "calculation approaches. Here a simple geometrical optics approach\n"
         "is followed. See the user guide for details.\n"
         "\n"
         "Defocusing is a special consideration for radio links. Two\n"
         "algorithms are at hand for estimating defocusing, simply denoted\n"
         "as method 1 and 2:\n"
         " 1: This algorithm is of general character. Defocusing is estimated\n"
         "    by making two path calculations with slightly shifted zenith\n"
         "    angles.\n"
         " 2: This method is restricted to satellite-to-satellite links, and\n"
         "    using a standard expression for such links, based on the\n"
         "    vertical gradient of the bending angle.\n"
         "Both methods are described more in detail in the user guide.\n"
         "The argument *defocus_shift* is used by both methods.\n"
         "\n"
         "The following auxiliary data can be obtained:\n"
         "  \"Pressure\": The pressure along the propagation path.\n"
         "     Size: [1,1,1,np].\n"
         "  \"Temperature\": The temperature along the propagation path.\n"
         "     Size: [1,1,1,np].\n"
         "  \"VMR, species X\": VMR of the species with index X (zero based).\n"
         "     For example, adding the string \"VMR, species 0\" extracts the\n"
         "     VMR of the first species. Size: [1,1,1,np].\n"
         "  \"Absorption, summed\": The total absorption matrix along the\n"
         "     path. Size: [nf,ns,ns,np].\n"
         "  \"Absorption, species X\": The absorption matrix along the path\n"
         "     for an individual species (X works as for VMR).\n"
         "     Size: [nf,ns,ns,np].\n"
         "  \"Particle extinction, summed\": The total extinction matrix over\n"
         "       all scattering elements along the path. Size: [nf,ns,ns,np].\n"
         "  \"PND, type X\": The particle number density for scattering element\n"
         "       type X (ie. corresponds to book X in pnd_field).\n"
         "       Size: [1,1,1,np].\n"
         "  \"Mass content, X\": The mass content for scattering element X.\n"
         "       This corresponds to column X in *particle_masses* (zero-\n"
         "       based indexing). Size: [1,1,1,np].\n"
         "* \"Impact parameter\": As normally defined for GNRSS radio\n"
         "       occultations (this equals the propagation path constant,\n"
         "       r*n*sin(theta)). Size: [1,1,1,1].\n"
         "* \"Free space loss\": The total loss due to the inverse square\n"
         "       law. Size: [1,1,1,1].\n"
         "  \"Free space attenuation\": The local attenuation due to the\n"
         "       inverse square law. Size: [1,1,1,np].\n"
         "* \"Atmospheric loss\": Total atmospheric attenuation, reported as\n"
         "       the transmission. Size: [nf,1,1,1].\n"
         "* \"Defocusing loss\": The total loss between the transmitter and\n"
         "       receiver due to defocusing. Given as a transmission.\n"
         "       Size: [1,1,1,1].\n"
         "* \"Faraday rotation\": Total rotation [deg] along the path, for\n"
         "     each frequency. Size: [nf,1,1,1].\n"
         "* \"Faraday speed\": The rotation per length unit [deg/m], at each\n"
         "     path point and each frequency. Size: [nf,1,1,np].\n"
         "* \"Extra path delay\": The time delay of the signal [s], compared\n"
         "       to the case of propagation through vacuum. Size: [1,1,1,1].\n"
         "* \"Bending angle\": As normally defined for GNRSS radio\n"
         "       occultations, in [deg]. Size: [1,1,1,1].\n"
         "where\n"
         "  nf: Number of frequencies.\n"
         "  ns: Number of Stokes elements.\n"
         "  np: Number of propagation path points.\n"
         "\n"
         "The auxiliary data are returned in *iy_aux* with quantities\n"
         "selected by *iy_aux_vars*. Most variables require that the method\n"
         "is called directly or by *iyCalc*. For calculations using *yCalc*,\n"
         "the selection is restricted to the variables marked with *.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "iy", "iy_aux", "ppath", "diy_dx" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "stokes_dim", "f_grid", "atmosphere_dim",
            "p_grid", "lat_grid", "lon_grid",
            "z_field", "t_field", "vmr_field", "abs_species",
            "wind_u_field", "wind_v_field", "wind_w_field", "mag_u_field",
            "mag_v_field", "mag_w_field", 
            "refellipsoid", "z_surface", "cloudbox_on", "cloudbox_limits", 
            "pnd_field", "scat_data", 
            "particle_masses", "iy_aux_vars", "jacobian_do", 
            "ppath_agenda", "ppath_step_agenda",
            "propmat_clearsky_agenda", "iy_transmitter_agenda",
            "iy_agenda_call1", "iy_transmission", "rte_pos", "rte_los", 
            "rte_pos2", "rte_alonglos_v", "ppath_lmax", "ppath_lraytrace" ),
        GIN(      "defocus_method", "defocus_shift" ),
        GIN_TYPE( "Index", "Numeric" ),
        GIN_DEFAULT( "1", "3e-3" ),
        GIN_DESC( "Selection of defocusing calculation method. See above.",
                  "Angular shift to apply in defocusing estimates." )
        ));
  */

  md_data_raw.push_back(MdRecord(
      NAME("iyReplaceFromAux"),
      DESCRIPTION(
          "Change of main output variable.\n"
          "\n"
          "With this method you can replace the content of *iy* with one of\n"
          "the auxiliary variables. The selected variable (by *aux_var*) must\n"
          "be part of *iy_aux_vars*. The corresponding data from *iy_aux* are\n"
          "copied to form a new *iy* (*iy_aux* is left unchanged). Elements of\n"
          "*iy* correponding to Stokes elements not covered by the auxiliary\n"
          "variable are just set to zero.\n"
          "\n"
          "Jacobian variables are not handled.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("iySurfaceCallAgendaX"),
      DESCRIPTION(
          "Switch between the elements of *iy_surface_agenda_array*.\n"
          "\n"
          "This method simply calls the agenda matching *surface_type* and\n"
          "returns the results. That is, the agenda in *iy_surface_agenda_array*\n"
          "with index *surface_type* (0-based) is called.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy", "diy_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
         "iy_unit",
         "iy_transmission",
         "iy_id",
         "cloudbox_on",
         "jacobian_do",
         "f_grid",
         "iy_main_agenda",
         "rtp_pos",
         "rtp_los",
         "rte_pos2",
         "iy_surface_agenda_array",
         "surface_type",
         "surface_type_aux"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("iySurfaceFastem"),
      DESCRIPTION(
          "Usage of FASTEM for emissivity and reflectivity of water surfaces.\n"
          "\n"
          "This method allows usage of the FASTEM model inside\n"
          "*iy_surface_agenda*. The aim is to use FASTEM in the exact same\n"
          "way as done in RTTOV. For example, the transmittance for down-\n"
          "welling radiation is considered. RTTOV os just 1D. Here 2D and 3D\n"
          "are handled as the 1D case, the down-welling radiation is just\n"
          "calculated for the directuon matching specular reflection.\n"
          "\n"
          "The wind direction is given as the azimuth angle, counted\n"
          "clockwise from north (i.e. an easterly wind is at 90 deg).\n"
          "This matches the general definition of azimuth inside ARTS.\n"
          "For 1D and 2D, the wind direction must be adjusted to match the\n"
          "fact that the line-of-sight is locked to be at 0 deg (180 for 2D\n"
          "in the case of a negative zenith angle). For 3D, the true wind\n"
          "direction shall be used.\n"
          "\n"
          "FASTEM is called by *FastemStandAlone*. See that WSM for further\n"
          "comments on variables and limitations.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy", "diy_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
         "iy_transmission",
         "iy_id",
         "jacobian_do",
         "atmosphere_dim",
         "nlte_field",
         "cloudbox_on",
         "stokes_dim",
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

  md_data_raw.push_back(MdRecord(
      NAME("iySurfaceRtpropAgenda"),
      DESCRIPTION(
          "Interface to *surface_rtprop_agenda* for *iy_surface_agenda*.\n"
          "\n"
          "This method is designed to be part of *iy_surface_agenda*. It\n"
          "determines the radiative properties of the surface by\n"
          "*surface_rtprop_agenda* and calculates the downwelling radiation\n"
          "by *iy_main_agenda*, and sums up the terms as described in AUG.\n"
          "That is, this WSM uses the output from *surface_rtprop_agenda*\n"
          "in a straightforward fashion.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy", "diy_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
         "iy_transmission",
         "iy_id",
         "jacobian_do",
         "atmosphere_dim",
         "nlte_field",
         "cloudbox_on",
         "stokes_dim",
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

  md_data_raw.push_back(MdRecord(
      NAME("iySurfaceRtpropCalc"),
      DESCRIPTION(
          "Applies *surface_los*, *surface_rmatrix* and *surface_emission*.\n"
          "\n"
          "This method is designed to be part of *iy_surface_agenda* and\n"
          "should be mandatory when using methods describing the surface\n"
          "radiative transfer properties by *surface_los*, *surface_rmatrix*\n"
          "and *surface_emission*. The task of this method is to apply these\n"
          "three WSVs to obtain the upwelling radiation from the surface.\n"
          "This upwelling radiation is the sum of surface emission and\n"
          "reflected downwelling radiation. The later part is calculated\n"
          "by calling *iy_main_agenda*. See further AUG.\n"),
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
         "iy_transmission",
         "iy_id",
         "jacobian_do",
         "jacobian_quantities",
         "atmosphere_dim",
         "nlte_field",
         "cloudbox_on",
         "stokes_dim",
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

  md_data_raw.push_back(MdRecord(
      NAME("iyTransmissionStandard"),
      DESCRIPTION(
          "Standard method for handling transmission measurements.\n"
          "\n"
          "Designed to be part of *iy_main_agenda*. Treatment of the cloudbox\n"
          "is incorporated (that is, no need to define *iy_cloudbox_agenda*).\n"
          "\n"
          "The transmitter is assumed to be placed at the end of provided *ppath*.\n"
          "The transmitted signal is taken from *iy_transmitter_agenda*. This\n"
          "signal is propagated along the path, considering attenuation alone.\n"
          "That is, the result of the method (*iy*) is the output of\n"
          "*iy_transmitter_agenda* multiplied with the transmission along the\n"
          "propagation path.\n"
          "\n"
          "As mentioned, the given *ppath* determines the position of the\n"
          "transmitter. For clear-sky and no modification of *ppath*, this\n"
          "means that the transitter will either be found at the surface or\n"
          "at the top-of-the-atmosphere. If you want to maintain this even with\n"
          "an active cloudbox, calculate *ppath* as\n"
          "     ppathCalc( cloudbox_on=0 )\n"
          "Without setting cloudbox_on=0, the transmitter will end up inside or\n"
          "at the boundary of the cloudbox.\n"
          "\n"
          "Some auxiliary radiative transfer quantities can be obtained. Auxiliary\n"
          "quantities are selected by *iy_aux_vars* and returned by *iy_aux*.\n"
          "Valid choices for auxiliary data are:\n"
          " \"Radiative background\": Index value flagging the radiative\n"
          "    background. The following coding is used: 0=space, 1=surface\n"
          "    and 2=cloudbox. The value is added to each column.\n"
          " \"Optical depth\": Scalar optical depth between the observation point\n"
          "    and the end of the present propagation path. Calculated based on\n"
          "    the (1,1)-element of the transmission matrix (1-based indexing),\n"
          "    i.e. only fully valid for scalar RT. The value is added to each\n"
          "    column.\n"),
      AUTHORS("Patrick Eriksson", "Richard Larsson"),
      OUT("iy",
          "iy_aux",
          "diy_dx",
          "ppvar_p",
          "ppvar_t",
          "ppvar_nlte",
          "ppvar_vmr",
          "ppvar_wind",
          "ppvar_mag",
          "ppvar_pnd",
          "ppvar_f",
          "ppvar_iy",
          "ppvar_trans_cumulat"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("diy_dx",
         "stokes_dim",
         "f_grid",
         "atmosphere_dim",
         "p_grid",
         "t_field",
         "nlte_field",
         "vmr_field",
         "abs_species",
         "wind_u_field",
         "wind_v_field",
         "wind_w_field",
         "mag_u_field",
         "mag_v_field",
         "mag_w_field",
         "cloudbox_on",
         "cloudbox_limits",
         "pnd_field",
         "dpnd_field_dx",
         "scat_species",
         "scat_data",
         "iy_aux_vars",
         "jacobian_do",
         "jacobian_quantities",
         "ppath",
         "propmat_clearsky_agenda",
         "water_p_eq_agenda",
         "iy_transmitter_agenda",
         "iy_agenda_call1",
         "iy_transmission",
         "rte_alonglos_v"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("iy_transmitterMultiplePol"),
      DESCRIPTION(
          "Transmitted signal having multiple polarisations.\n"
          "\n"
          "The method is intended to be part of *iy_transmitter_agenda*. It\n"
          "sets *iy* to describe the transmitted signal/pulses. The polarisation\n"
          "state is taken from *instrument_pol*, where *instrument_pol* must\n"
          "contain an element for each frequency in *f_grid*. The transmitted\n"
          "signal/pulses are set to be of unit magnitude, such as [1,1,0,0].\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("stokes_dim", "f_grid", "instrument_pol"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("iy_transmitterSinglePol"),
      DESCRIPTION(
          "Transmitted signal having a single polarisations.\n"
          "\n"
          "The method is intended to be part of *iy_transmitter_agenda*. It\n"
          "sets *iy* to describe the transmitted pulses/signal. The polarisation\n"
          "state is taken from *instrument_pol*, where *instrument_pol* must contain\n"
          "a single value. This polarisation state is applied for all\n"
          "frequencies. The transmitted pulses/signals are set to be of unit\n"
          "magnitude, such as [1,1,0,0].\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("iy"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("stokes_dim", "f_grid", "instrument_pol"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddAbsSpecies"),
      DESCRIPTION(
          "Includes an absorption species in the Jacobian.\n"
          "\n"
          "For 1D or 2D calculations the latitude and/or longitude grid of\n"
          "the retrieval field should set to have zero length.\n"
          "\n"
          "These retrieval units are at hand for all gas species:\n"
          "   \"vmr\"    : Volume mixing ratio.\n"
          "   \"nd\"     : Number density.\n"
          "   \"rel\"    : Relative unit (e.g. 1.1 means 10% more of the gas).\n"
          "\n"
          "For water vapour, also these units are at hand:\n"
          "   \"rh\"     : Relative humidity.\n"
          "   \"q\"      : Specific humidity.\n"
          "\n"
          "Note that *for_species_tag* is used to indicate if species tag VMR,\n"
          "rather than atmospheric gas VMR is calculated. Set it to 0 and we\n"
          "calculate the atmospheric gas VMR, but this only works for \"analytical\".\n"
          "\n"
          "Note that the Jacobian is set to zero where volume mixing ratio equals zero.\n"
          "\n"
          "The number of elements added to the state vector (*x*) is:\n"
          "   n_g1 * n_g2 * n_g3\n"
          "where n_g1, n_g2 and n_g3 are the length of GIN *g1*, *g2* and *g3*,\n"
          "respectively. Here empty vectors should be considered to have a length 1.\n"
          "The elements are sorted with pressure as innermost loop, followed by\n"
          "latitude and longitude as outermost loop.\n"),
      AUTHORS("Mattias Ekstrom", "Patrick Eriksson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid"),
      GIN("g1",
          "g2",
          "g3",
          "species",
          "unit",
          "for_species_tag"),
      GIN_TYPE("Vector",
               "Vector",
               "Vector",
               "String",
               "String",
               "Index"),
      GIN_DEFAULT(
          NODEF, NODEF, NODEF, NODEF, "vmr", "1"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddBasicCatalogParameter"),
      DESCRIPTION(
          "Includes a basic catalog parameter in the Jacobian. These are constant\n"
          "over all layers and so only a single vector output is returned.\n"
          "\n"
          "The only basic catalog parameters currently supported are:\n"
          "   \"Line Strength\"\n"
          "   \"Line Center\"\n"
          "\n"
          "The *catalog_identity* should be able to identify one or many\n"
          "lines in the catalog used for calculating the spectral absorption.\n"
          "Note that partial matching for energy levels are allowed but not\n"
          "recommended, as it is somewhat nonsensical to add multiple parameters\n"
          "\n"
          "Also note *jacobianAddShapeCatalogParameter* as this allows addition\n"
          "of shape parameters, e.g., pressure broadening coefficients\n"
          "\n"
          "Each call to this function adds just a single value to *x*.\n"
          "\n"
          "Example given the catalog_identity=\"O2-66 TR UP v1 0 J 1 LO v1 0 J 0\",\n"
          "only the O2 ground-level 119 GHz line can be accessed and only its\n"
          "catalog_parameter will be accessed.  However, the more lenient\n"
          "catalog_identity=\"O2-66 TR UP J 1 LO J 0\" may be used, but then the\n"
          "118 GHz line belonging to v1=1 branch will be added to the same *x*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddBasicCatalogParameters"),
      DESCRIPTION(
          "See *jacobianAddBasicCatalogParameter*.\n"
          "\n"
          "This adds a multiple of parameters for first each catalog_identity in\n"
          "catalog_identities and then for each catalog_parameter in catalog_parameters\n"
          "by looping calls to *jacobianAddBasicCatalogParameter* over these input\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddFreqShift"),
      DESCRIPTION(
          "Includes a frequency fit of shift type in the Jacobian.\n"
          "\n"
          "Retrieval of deviations between nominal and actual backend\n"
          "frequencies can be included by this method. The assumption here is\n"
          "that the deviation is a constant off-set, a shift, common for all\n"
          "frequencies (and not varying between measurement blocks).\n"
          "\n"
          "This method adds one element to the state vector (*x*).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddFreqStretch"),
      DESCRIPTION(
          "Includes a frequency fit of stretch type in the Jacobian.\n"
          "\n"
          "Retrieval of deviations between nominal and actual backend\n"
          "frequencies can be included by this method. The assumption here is\n"
          "that the deviation varies linearly over the frequency range\n"
          "(following ARTS basis function for polynomial order 1).\n"
          "\n"
          "This method adds one element to the state vector (*x*).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddShapeCatalogParameter"),
      DESCRIPTION(
          "Adds a line shape parameter to the Jacobian calculations. These\n"
          "are constant over all levels so only a single *x*-value is added\n"
          "\n"
          "Line function parameter assume the derivatives\n"
          "of internal pressure broadening and line mixing\n"
          "functionality follows a simply f(T, T0, X0, X1, X2)\n"
          "format.  The shape of the function f() is determined by\n"
          "input catalog; please see the ARTS documentation for more\n"
          "details\n"
          "\n"
          "The input are as follows:\n"
          "    line_identity: Identifier of preferably a single line\n"
          "    species:       A SpeciesTag, e.g., \"O2\" or \"H2O\" for common species.\n"
          "                   Note that \"SELF\" and \"AIR\" tags are used for shape parameters\n"
          "                   affected by self and air-broadening, respectively.\n"
          "    variable:      A variable supported by the line, these can be\n"
          "                      \"G0\":  Speed-independent pressure broadening\n"
          "                      \"G2\":  Speed-dependent pressure broadening\n"
          "                      \"D0\":  Speed-independent pressure shift\n"
          "                      \"D2\":  Speed-dependent pressure shift\n"
          "                      \"FVC\": Frequency of velocity changing collisions\n"
          "                      \"ETA\": partial correlation between velocity and\n"
          "                               rotational state changes due to collisions\n"
          "                      \"Y\":   First order line-mixing parameter\n"
          "                      \"G\":   Second order line-mixing parameter for strength\n"
          "                      \"DV\":  Second order line-mixing parameter for shifting\n"
          "    coefficient:   A coefficient in the model to compute the above parameters.\n"
          "\n"
          "Note that we cannot test if the line in question supports the variable and\n"
          "coefficient at the level of this function, so many errors will only be reported\n"
          "at a later stage\n"
          "\n"
          "For other spectroscopic parameters, see *jacobianAddBasicCatalogParameter*.\n"
          "Also see said function for an example of how to set the QuantumIdentifier\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddShapeCatalogParameters"),
      DESCRIPTION(
          "See *jacobianAddShapeCatalogParameter* for information on\n"
          "the GIN parameters\n"
          "\n"
          "This function accepts the same input but for lists of data.\n"
          "The function loops over each input list\n"
          "individually and appends the information to *jacobian_quantities*.\n"
          "\n"
          "Special \"ALL\" for 1 length *variables* and *coefficients* are\n"
          "allowed to compute all variables/coefficients in the order described\n"
          "in the description of *jacobianAddShapeCatalogParameter*\n"
          "\n"
          "For example, if *line_identities* have length 5, *species* length 4,\n"
          "*variables* length 3, and *coefficients* length 2, there will be\n"
          "5*4x3x2 = 120 new additions to *jacobian_quantities* in the order:\n"
          "\t[{line_identities[0], species[0], variables[0] coefficients[0]}]\n"
          "\t[{line_identities[0], species[0], variables[0] coefficients[1]}]\n"
          "\t[{line_identities[0], species[0], variables[1] coefficients[0]}]\n"
          "\t[{line_identities[0], species[0], variables[1] coefficients[1]}]\n"
          "\t[{line_identities[0], species[0], variables[2] coefficients[0]}]\n"
          "\t[{line_identities[0], species[0], variables[2] coefficients[1]}]\n"
          "\t[{line_identities[0], species[1], variables[0] coefficients[0]}]\n"
          "\t...\n"
          "\t[{line_identities[4], species[3], variables[1] coefficients[1]}]\n"
          "\t[{line_identities[4], species[3], variables[2] coefficients[0]}]\n"
          "\t[{line_identities[4], species[3], variables[2] coefficients[1]}]\n"
          "or in words: lines first, then species, then variables, then coefficients\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddMagField"),
      DESCRIPTION(
          "Includes one magnetic field component in the Jacobian.\n"
          "\n"
          "The method follows the pattern of other Jacobian methods. The\n"
          "calculations can only be performed by analytic expressions.\n"
          "\n"
          "The magnetic field components are retrieved separately, and,\n"
          "hence, the argument *component* can be  \"u\", \"v\", \"w\",\n"
          "and \"strength\".\n"
          "\n"
          "The number of elements added to the state vector (*x*) is:\n"
          "   n_g1 * n_g2 * n_g3\n"
          "where n_g1, n_g2 and n_g3 are the length of GIN *g1*, *g2* and *g3*,\n"
          "respectively. Here empty vectors should be considered to have a length 1.\n"
          "The elements are sorted with pressure as innermost loop, followed by\n"
          "latitude and longitude as outermost loop.\n"
          "\n"
          "The dB-parameter is only used for Faraday rotation\n"),
      AUTHORS("Patrick Eriksson", "Richard Larsson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid"),
      GIN("g1", "g2", "g3", "component", "dB"),
      GIN_TYPE("Vector", "Vector", "Vector", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "v", "1.0e-7"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Magnetic field component to retrieve",
               "Magnetic field perturbation")));

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddNLTE"),
      DESCRIPTION(
          "Experimental NLTE Jacobian.\n"
          "\n"
          "Intention: Adds the nlte_field level distribution per atmospheric grid\n"
          "to the Jacobian.\n"
          "\n"
          "The number of elements added to the state vector (*x*) is:\n"
          "   n_g1 * n_g2 * n_g3\n"
          "where n_g1, n_g2 and n_g3 are the length of GIN *g1*, *g2* and *g3*,\n"
          "respectively. Here empty vectors should be considered to have a length 1.\n"
          "The elements are sorted with pressure as innermost loop, followed by\n"
          "latitude and longitude as outermost loop.\n"
          "\n"
          "The QuantumIdentifier should identify a single energy level, such as:\n"
          "\"H2O-161 EN J 1 Ka 0 Kc 1\", for one of the lower levels in the chains\n"
          "of transitions of water.  Note that using this method directly is not\n"
          "best practice, as the quantum identifiers of the levels have to be known\n"
          "at an early stage in NLTE calculations, and will usually populate the\n"
          "*nlte_level_identifiers* variable, meaning it is better to use *jacobianAddNLTE*\n"
          "directly than to individually call this function\n"),
      AUTHORS("Richard Larsson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid"),
      GIN("g1", "g2", "g3", "energy_level_identity", "dx"),
      GIN_TYPE("Vector", "Vector", "Vector", "QuantumIdentifier", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, "1.0e-3"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Identifier to the eneregy level",
               "Perturbation of value if required by method")));

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddNLTEs"),
      DESCRIPTION(
          "Experimental NLTE Jacobian.  Same as *jacobianAddNLTE* but for\n"
          "many levels\n"
          "\n"
          "Adds energy_level_identities.nelem() times as many arguments to *x*\n"
          "as *jacobianAddNLTE*, ordered as energy_level_identities describes\n"
          "\n"
          "This method is preferred to *jacobianAddNLTE*, since *energy_level_identities*\n"
          "is conveniently almost always the same as *nlte_level_identifiers*\n"),
      AUTHORS("Richard Larsson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid"),
      GIN("g1", "g2", "g3", "energy_level_identities", "dx"),
      GIN_TYPE(
          "Vector", "Vector", "Vector", "ArrayOfQuantumIdentifier", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF, "1.0e-3"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Identifiers to the eneregy level",
               "Perturbation of value if required by method")));

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddPointingZa"),
      DESCRIPTION(
          "Adds sensor pointing zenith angle off-set jacobian.\n"
          "\n"
          "Retrieval of deviations between nominal and actual zenith angle of\n"
          "the sensor can be included by this method. The weighing functions\n"
          "can be calculated in several ways:\n"
          "   calcmode = \"recalc\": Recalculation of pencil beam spectra,\n"
          "      shifted with *dza* from nominal values. A single-sided\n"
          "      perturbation is applied (towards higher zenith angles).\n"
          "   calcmode = \"interp\": Inter/extrapolation of existing pencil\n"
          "       beam spectra. For this option, allow some extra margins for\n"
          "       zenith angle grids, to avoid artifacts when extrapolating\n"
          "       the data (to shifted zenith angles). The average of a\n"
          "       negative and a positive shift is taken."
          "\n"
          "The interp option is recommended. It should in general be both\n"
          "faster and more accurate (due to the double sided disturbance).\n"
          "In addition, it is less sensitive to the choice of dza (as long\n"
          "as a small value is applied).\n"
          "\n"
          "The pointing off-set can be modelled to be time varying. The time\n"
          "variation is then described by a polynomial (with standard base\n"
          "functions). For example, a polynomial order of 0 means that the\n"
          "off-set is constant in time. If the off-set is totally uncorrelated\n"
          "between the spectra, set the order to -1.\n"
          "\n"
          "The number of elements added to the state vector (*x*) is\n"
          "  if poly_order < 0 : length of *sensor_time*\n"
          "         otherwise : poly_order+1\n"
          "In the first case, the order in *x* matches *sensor_time*. In the second\n"
          "case, the coefficient for polynomial order 0 comes first etc.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddPolyfit"),
      DESCRIPTION(
          "Includes polynomial baseline fit in the Jacobian.\n"
          "\n"
          "This method deals with retrieval of disturbances of the spectra\n"
          "that can be described by an additive term, a baseline off-set.\n"
          "\n"
          "The baseline off-set is here modelled as a polynomial. The\n"
          "polynomial spans the complete frequency range spanned by\n"
          "*sensor_response_f_grid* and the method should only of interest for\n"
          "cases with no frequency gap in the spectra. The default assumption\n"
          "is that the off-set differs between all spectra, but it can also be\n"
          "assumed that the off-set is common for all e.g. line-of-sights.\n"
          "\n"
          "If the simulation/retrieval deals with a single spectrum, the number\n"
          "of elements added to the state vector (*x*) is poly_order+1. The\n"
          "coefficient for polynomial order 0 comes first etc. The same is true\n"
          "if *no_pol_variation*, *no_los_variation* and *no_mblock_variation*\n"
          "all are set to 1, even if several spectra are involved. Otherwise the"
          "number of elements added to *x* depends on the number of spectra and\n"
          "the settings of *no_pol_variation*, *no_los_variation* and \n"
          "*no_mblock_variation*. The coefficients of the different polynomial\n"
          "orders are treated as separate retrieval quantities. That is, the\n"
          "the elements associated with polynomial order 0 are grouped and form\n"
          "together a retrieval quantity. The coefficients for higher polynomial\n"
          "orders are treated in the same way.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddScatSpecies"),
      DESCRIPTION(
          "Includes a scattering species in the Jacobian.\n"
          "\n"
          "For 1D or 2D calculations the latitude and/or longitude grid of\n"
          "the retrieval field should set to have zero length.\n"
          "\n"
          "The number of elements added to the state vector (*x*) is:\n"
          "   n_g1 * n_g2 * n_g3\n"
          "where n_g1, n_g2 and n_g3 are the length of GIN *g1*, *g2* and *g3*,\n"
          "respectively. Here empty vectors should be considered to have a length 1.\n"
          "The elements are sorted with pressure as innermost loop, followed by\n"
          "latitude and longitude as outermost loop.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddSinefit"),
      DESCRIPTION(
          "Includes sinusoidal baseline fit in the Jacobian.\n"
          "\n"
          "Works as *jacobianAddPolyFit*, beside that a series of sine and\n"
          "cosine terms are used for the baseline fit.\n"
          "\n"
          "For each value in *period_lengths one sine and one cosine term are\n"
          "included (in mentioned order). By these two terms the amplitude and\n"
          "\"phase\" for each period length can be determined. The sine and\n"
          "cosine terms have value 0 and 1, respectively, for first frequency.\n"
          "\n"
          "If the simulation/retrieval deals with a single spectrum, the number\n"
          "of elements added to the state vector (*x*) is 2*nperiods, where\n"
          "*nperiods* is the length of *period_lengths*. The same is true\n"
          "if *no_pol_variation*, *no_los_variation* and *no_mblock_variation*\n"
          "all are set to 1, even if several spectra are involved. Otherwise the"
          "number of elements added to *x* depends on the number of spectra and\n"
          "the settings of *no_pol_variation*, *no_los_variation* and \n"
          "*no_mblock_variation*. The sine and cosine terms for each period\n"
          "length are treated as a  separate retrieval quantities. That is, the\n"
          "the elements associated with the first period length are grouped and\n"
          "form together a retrieval quantity, etc. Inside each retrieval quantity\n"
          "the pairs of sine and cosine terms are kept together, in given order.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddSpecialSpecies"),
      DESCRIPTION(
          "Includes a special absorption species in the Jacobian.\n"
          "\n"
          "Similar to *jacobianAddAbsSpecies* but only for number densities.\n"
          "\n"
          "Species allowed are:\n"
          "    \"electrons\"\n"
          "    \"particulates\"\n"
          "\n"
          "Note that the average of all particulates are used to scale its\n"
          "*jacobian*, so this method works best when only one type of\n"
          "particulate is being used, i.e., when *scat_data* has only one\n"
          "scattering species.\n"
          "\n"
          "The number of elements added to the state vector (*x*) is:\n"
          "   n_g1 * n_g2 * n_g3\n"
          "where n_g1, n_g2 and n_g3 are the length of GIN *g1*, *g2* and *g3*,\n"
          "respectively. Here empty vectors should be considered to have a length 1.\n"
          "The elements are sorted with pressure as innermost loop, followed by\n"
          "latitude and longitude as outermost loop.\n"),
      AUTHORS("Richard Larsson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddSurfaceQuantity"),
      DESCRIPTION(
          "Includes a surface quantity in the Jacobian.\n"
          "\n"
          "The quantity is specified by the GIN-variable *quantity*. The name\n"
          "of the quantity must match the name used in *surface_prop_names*.\n"
          "\n"
          "For 1D or 2D calculations the latitude and/or longitude grid of\n"
          "the retrieval field should set to have zero length.\n"
          "\n"
          "The number of elements added to the state vector (*x*) is:\n"
          "   n_g1 * n_g2\n"
          "where n_g1 and n_g2 are the length of GIN *g1* and *g2*, respectively.\n"
          "Here empty vectors should be considered to have a length 1.\n"
          "The elements are sorted with latitude as innermost loop and longitude\n"
          "as outermost loop.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "lat_grid",
         "lon_grid"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddTemperature"),
      DESCRIPTION(
          "Includes atmospheric temperatures in the Jacobian.\n"
          "\n"
          "The calculations are performed by (semi-)analytical expressions.\n"
          "Hydrostatic equilibrium (HSE) can be included.\n"
          "\n"
          "The analytical calculation approach neglects so far refraction\n"
          "totally, but considers the local effect of HSE.\n"
          "The later should be accaptable for observations around zenith and\n"
          "nadir. There is no warning if the method is applied incorrectly, \n"
          "with respect to these issues. Note that the argument *hse* of this\n"
          "WSM only refers to the Jacobian calculation, if the model and/or\n"
          "retrieved atmosphere actually fulfils HSE or not is governed in\n"
          "other manners.\n"
          "\n"
          "The calculations (both options) assume that gas species are defined\n"
          "in VMR (a change in temperature then changes the number density). \n"
          "This has the consequence that retrieval of temperatures and number\n"
          "density can not be mixed. Neither any warning here!\n"
          "\n"
          "The number of elements added to the state vector (*x*) is:\n"
          "   n_g1 * n_g2 * n_g3\n"
          "where n_g1, n_g2 and n_g3 are the length of GIN *g1*, *g2* and *g3*,\n"
          "respectively. Here empty vectors should be considered to have a length 1.\n"
          "The elements are sorted with pressure as innermost loop, followed by\n"
          "latitude and longitude as outermost loop.\n"),
      AUTHORS("Mattias Ekstrom", "Patrick Eriksson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid"),
      GIN("g1", "g2", "g3", "hse"),
      GIN_TYPE("Vector", "Vector", "Vector", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "on"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Flag to assume HSE or not (\"on\" or \"off\").")));

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAddWind"),
      DESCRIPTION(
          "Includes one atmospheric wind component in the Jacobian.\n"
          "\n"
          "The method follows the pattern of other Jacobian methods. The\n"
          "calculations can only be performed by analytic expressions.\n"
          "Some lower level function depends on frequency perturbations,\n"
          "however, so therefore a frequency perturbation *df* is required\n"
          "and as a consequence *abs_f_interp_order* must be > 0.\n"
          "\n"
          "The wind field components are retrieved separately, and,\n"
          "hence, the argument *component* can be \"u\", \"v\" or \"w\" \n"
          "for vector components, or just \"strength\" for total wind speed.\n"
          "\n"
          "The number of elements added to the state vector (*x*) is:\n"
          "   n_g1 * n_g2 * n_g3\n"
          "where n_g1, n_g2 and n_g3 are the length of GIN *g1*, *g2* and *g3*,\n"
          "respectively. Here empty vectors should be considered to have a length 1.\n"
          "The elements are sorted with pressure as innermost loop, followed by\n"
          "latitude and longitude as outermost loop.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid"),
      GIN("g1", "g2", "g3", "component", "dfrequency"),
      GIN_TYPE("Vector", "Vector", "Vector", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "v", "0.1"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retrieval grid.",
               "Wind component to retrieve",
               "This is the frequency perturbation")));

  md_data_raw.push_back(MdRecord(
      NAME("jacobianAdjustAndTransform"),
      DESCRIPTION(
          "Applies adjustments and transformations on *jacobian*.\n"
          "\n"
          "The method handles two tasks:\n"
          "1. The retrieval transformations set by the user can not be applied\n"
          "onthe  Jacobian inside *yCalc*. Transformations are instead applied\n"
          "by calling this method.\n"
          "2. It applies required adjustments of the Jacoboan. So far there is\n"
          "only one possible adjustment. If any absorption species uses the \"rel\"\n"
          "unit, an adjustment is needed for later iterations of the inversion.\n"
          "\n"
          "If no tranformations are selected and the \"rel\" option is not used at\n"
          "all, there is no need to call this method(, but you can still include it\n"
          "without causing any error, the calculations will just be a bit slower).\n"
          "Otherwise, this method should be called, typically as part of\n"
          "*inversion_iteration_agenda*.\n"
          "\n"
          "The method accepts if *jacobian* is empty, and then does, nothing.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianCalcDoNothing"),
      DESCRIPTION(
          "This function doesn't do anything. It just exists to satisfy\n"
          "the input and output requirement of the *jacobian_agenda*.\n"
          "\n"
          "This method is added to *jacobian_agenda* by *jacobianAddAbsSpecies*\n"
          "and some similar methods, and it should normally not be called by\n"
          "the user.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianCalcFreqShift"),
      DESCRIPTION(
          "Calculates frequency shift jacobians by interpolation\n"
          "of *iyb*.\n"
          "\n"
          "This function is added to *jacobian_agenda* by jacobianAddFreqShift\n"
          "and should normally not be called by the user.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian",
         "mblock_index",
         "iyb",
         "yb",
         "stokes_dim",
         "f_grid",
         "mblock_dlos_grid",
         "sensor_response",
         "jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("jacobianCalcFreqStretch"),
      DESCRIPTION(
          "Calculates frequency stretch jacobians by interpolation\n"
          "of *iyb*.\n"
          "\n"
          "This function is added to *jacobian_agenda* by jacobianAddFreqStretch\n"
          "and should normally not be called by the user.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian",
         "mblock_index",
         "iyb",
         "yb",
         "stokes_dim",
         "f_grid",
         "mblock_dlos_grid",
         "sensor_response",
         "sensor_response_pol_grid",
         "sensor_response_f_grid",
         "sensor_response_dlos_grid",
         "jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("jacobianCalcPointingZaInterp"),
      DESCRIPTION("Calculates zenith angle pointing deviation jacobians by\n"
                  "inter-extrapolation of *iyb*.\n"
                  "\n"
                  "This function is added to *jacobian_agenda* by\n"
                  "jacobianAddPointingZa and should normally not be\n"
                  "called by the user.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian",
         "mblock_index",
         "iyb",
         "yb",
         "stokes_dim",
         "f_grid",
         "sensor_los",
         "mblock_dlos_grid",
         "sensor_response",
         "sensor_time",
         "jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("jacobianCalcPointingZaRecalc"),
      DESCRIPTION("Calculates zenith angle pointing deviation jacobians by\n"
                  "recalulation of *iyb*.\n"
                  "\n"
                  "This function is added to *jacobian_agenda* by\n"
                  "jacobianAddPointingZa and should normally not be\n"
                  "called by the user.\n"),
      AUTHORS("Mattias Ekstrom", "Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian",
         "mblock_index",
         "iyb",
         "yb",
         "atmosphere_dim",
         "nlte_field",
         "cloudbox_on",
         "stokes_dim",
         "f_grid",
         "sensor_pos",
         "sensor_los",
         "transmitter_pos",
         "mblock_dlos_grid",
         "sensor_response",
         "sensor_time",
         "iy_unit",
         "iy_main_agenda",
         "geo_pos_agenda",
         "jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("jacobianCalcPolyfit"),
      DESCRIPTION(
          "Calculates jacobians for polynomial baseline fit.\n"
          "\n"
          "This function is added to *jacobian_agenda* by jacobianAddPolyfit\n"
          "and should normally not be called by the user.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianCalcSinefit"),
      DESCRIPTION(
          "Calculates jacobians for sinusoidal baseline fit.\n"
          "\n"
          "This function is added to *jacobian_agenda* by jacobianAddPolyfit\n"
          "and should normally not be called by the user.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianClose"),
      DESCRIPTION(
          "Closes the array of retrieval quantities and prepares for\n"
          "calculation of the Jacobian matrix.\n"
          "\n"
          "This function closes the *jacobian_quantities* array and sets\n"
          "*jacobian_do* to 1.\n"
          "\n"
          "Retrieval quantities should not be added after a call to this WSM.\n"
          "No calculations are performed here.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianFromTwoY"),
      DESCRIPTION(
          "Sets *jacobian* based on the difference vetween two measurement vectors.\n"
          "\n"
          "This function assumes that *y_pert* contains a measurement calculated\n"
          "with some variable perturbed, in comparison to the calculation\n"
          "behind *y*. The function takes the differences between *y_pert*\n"
          "and *y* to form a numerical derived estimate of *jacobian*.\n"
          "This gives a Jacobian wit a single column.\n"
          "\n"
          "*jacobian* equals here: (y_pert-y)/pert_size.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("y"),
      GIN("y_pert","pert_size"),
      GIN_TYPE("Vector","Numeric"),
      GIN_DEFAULT(NODEF,NODEF),
      GIN_DESC("Perturbed measurement vector",
               "Size of perturbation behind spectra in *ybatch*.")));

  md_data_raw.push_back(MdRecord(
      NAME("jacobianFromYbatch"),
      DESCRIPTION(
          "Sets *jacobian* based on perturbation calcuations.\n"
          "\n"
          "This function assumes that *ybatch* contains spectra calculated\n"
          "with some variable perturbed, in comparison to the calculation\n"
          "behind *y*. The function takes the differences between *ybatch*\n"
          "and *y* to form a numerical derived estimate of *jacobian*.\n"
          "\n"
          "Column i of *jacobian* equals: (ybatch[i]-y)/pert_size.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ybatch","y"),
      GIN("pert_size"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Size of perturbation behind spectra in *ybatch*.")));

  md_data_raw.push_back(MdRecord(
      NAME("jacobianInit"),
      DESCRIPTION(
          "Initialises the variables connected to the Jacobian matrix.\n"
          "\n"
          "This function initialises the *jacobian_quantities* array so\n"
          "that retrieval quantities can be added to it. Accordingly, it has\n"
          "to be called before any calls to jacobianAddTemperature or\n"
          "similar methods.\n"
          "\n"
          "The Jacobian quantities are initialised to be empty.\n"),
      AUTHORS("Mattias Ekstrom"),
      OUT("jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("jacobianOff"),
      DESCRIPTION(
          "Makes mandatory initialisation of some jacobian variables.\n"
          "\n"
          "Some clear-sky jacobian WSVs must be initialised even if no such\n"
          "calculations will be performed.  This is handled with this method.\n"
          "That is, this method must be called when no clear-sky jacobians\n"
          "will be calculated (even if cloudy-sky jacobians are calculated!).\n"
          "\n"
          "Sets *jacobian_do* to 0.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("jacobian_do",
          "jacobian_agenda",
          "jacobian_quantities"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("jacobianSetAffineTransformation"),
      DESCRIPTION(
          "Adds an affine transformation of the last element of\n"
          "*jacobian_quantities*.\n"
          "\n"
          "See *jacobianSetFuncTransformation* for  a general description of how\n"
          "retrieval transformations are defined. Transformations are not applied by\n"
          "methods such as*yCalc*. Instead, the method *jacobianAdjustAndTransform*\n"
          "must be called to activate the transformations.\n"
          "\n"
          "The affine transformation is specified by a transformation matrix, A,\n"
          "and an offset vector, b. These two are applied as described in\n"
          "*jacobianSetFuncTransformation*.\n"
          "\n"
          "The transformations is applied as\n"
          "   x = A * ( z - b )\n"
          "where z is the retrieval quantity on the standard retrieval grids\n"
          "and x is the final state vector.\n"
          "\n"
          "So far, the following must be true for valid A-matrices\n"
          "   z = A'*x + b\n"
          "That is, the reversed transformation is given by A transposed.\n"
          "\n"
          "This method must only be called if an affine transformation is wanted.\n"
          "Default is to make no such tranformation at all.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("jacobianSetFuncTransformation"),
      DESCRIPTION(
          "Sets the functional transformation of the last element of\n"
          "*jacobian_quantities*.\n"
          "\n"
          "See below for a general description of how retrieval transformations\n"
          "are defined. Transformations are not applied by methods such as*yCalc*.\n"
          "Instead, the method *jacobianAdjustAndTransform* must be called to\n"
          "activate the transformations.\n"
          "\n"
          "The following transformations can be selected (by *transformation_func*):\n"
          "   log   : The natural logarithm\n"
          "   log10 : The base-10 logarithm\n"
          "   atanh : Area hyperbolic tangent \n"
          "   none  : No transformation at all\n"
          "\n"
          "This method needs only to be called if a functional transformation\n"
          "is wanted. Default is to make no such tranformation at all (i.e.\n"
          "the option \"none\" exists only for reasons of flexibility).\n"
          "\n"
          "The log-options are applied as log(z-z_min) and log10(z-z_min).\n"
          "The default for *z_min* is zero, but by changing it the lower limit\n"
          "for z can be changed. Note that *z_min* becomes the lower limit for\n"
          "allowed values of z. The GIN *z_max* is here ignored.\n"
          "\n"
          "For the atanh-option, also *z_max* is considered. This transformation\n"
          "is applied as atanh((2(z-z_min)/(z_max-z_min))-1). As above,*z_min*\n"
          "is lower limit for allowed values of z. On the other hand, *z_max*\n"
          "eines the upper limit for z.\n"
          "\n"
          "The GIN *tfunc_parameter* is so far only used for atanh. The parameter\n"
          "specifies the maximum allowed value allowed for u. That is, the valid\n"
          "range for u becomes ]0,tfunc_parameter[. Note that log and log10\n"
          "demands/ensures that u > 0, but implies no upper limit.\n"
          "\n"
          "General handling of retrieval units and transformations:\n"
          "---\n"
          "Default is that quantities are retrieved as defined in ARTS, but\n"
          "both some unit conversion and transformations are provided. These\n"
          "operations are applied as:\n"
          "   x = A * ( f(u(z)) - b ) \n"
          "where\n"
          "   z is the quantity as defined ARTS\n"
          "   u represents the change of unit\n"
          "   f is the transformation function\n"
          "   A and b define together an affine transformation\n"
          "   x is the retrieved quantity\n"
          "For example, this systen allows to retrive a principal component\n"
          "representation (A and b) of the log (f) of relative humidity (u).\n"
          "\n"
          "Change of unit is selected by the quantity specific jacobian-add\n"
          "methods (so far only at hand for gas species). \n"
          "\n"
          "Activating a transformation function is done by this method. Note\n"
          "that the functions are defined as the transformation from z to x.\n"
          "For more details on affine transformations, see\n"
          "*jacobiaSetAffineTransformation*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("lat_gridFromRawField"),
      DESCRIPTION(
          "Sets *lat_grid* according to given raw atmospheric field's lat_grid.\n"
          "Similar to *p_gridFromZRaw*, but acting on a generic *GriddedField3*\n"
          "(e.g., a wind or magnetic field component).\n"),
      AUTHORS("Jana Mendrok"),
      OUT("lat_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("field_raw"),
      GIN_TYPE("GriddedField3"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("A raw atmospheric field.")));

  md_data_raw.push_back(MdRecord(
      NAME("lon_gridFromRawField"),
      DESCRIPTION(
          "Sets *lon_grid* according to given raw atmospheric field's lat_grid.\n"
          "Similar to *p_gridFromZRaw*, but acting on a generic *GriddedField3*\n"
          "(e.g., a wind or magnetic field component).\n"),
      AUTHORS("Jana Mendrok"),
      OUT("lon_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("field_raw"),
      GIN_TYPE("GriddedField3"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("A raw atmospheric field.")));

  md_data_raw.push_back(MdRecord(
      NAME("MagFieldsCalc"),
      DESCRIPTION(
          "Interpolation of raw magnetic fields to calculation grids.\n"
          "Heritage from *AtmFieldsCalc*\n"
          "\n"
          "Internally, *MagFieldsCalc* applies *GriddedFieldPRegrid* and\n"
          "*GriddedFieldLatLonRegrid*. Generally, 'half-grid-step' extrapolation\n"
          "is allowed and applied.\n"),
      AUTHORS("Richard Larsson"),
      OUT("mag_u_field", "mag_v_field", "mag_w_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("p_grid",
         "lat_grid",
         "lon_grid",
         "mag_u_field_raw",
         "mag_v_field_raw",
         "mag_w_field_raw",
         "atmosphere_dim"),
      GIN("interp_order"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Interpolation order (1=linear interpolation).")));

  md_data_raw.push_back(MdRecord(
      NAME("MagFieldsCalcExpand1D"),
      DESCRIPTION(
          "Interpolation of 1D raw atmospheric fields to create 2D or 3D\n"
          "homogeneous magnetic fields.  Derived from *AtmFieldsCalcExpand1D*\n"
          "\n"
          "The method works as *MagFieldsCalc*, but accepts only raw 1D\n"
          "magnetic fields. The raw data is interpolated to *p_grid* and\n"
          "the obtained values are applied for all latitudes, and also\n"
          "longitudes for 3D, to create a homogeneous atmosphere.\n"),
      AUTHORS("Richard Larsson"),
      OUT("mag_u_field", "mag_v_field", "mag_w_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("p_grid",
         "lat_grid",
         "lon_grid",
         "mag_u_field_raw",
         "mag_v_field_raw",
         "mag_w_field_raw",
         "atmosphere_dim"),
      GIN("interp_order"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Interpolation order (1=linear interpolation).")));

  md_data_raw.push_back(MdRecord(
      NAME("MagFieldsFromAltitudeRawCalc"),
      DESCRIPTION(
          "Regrids the rawfield by lat-lon and interpolates to z_field.\n"),
      AUTHORS("Richard Larsson"),
      OUT("mag_u_field", "mag_v_field", "mag_w_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("lat_grid",
         "lon_grid",
         "z_field",
         "mag_u_field_raw",
         "mag_v_field_raw",
         "mag_w_field_raw", ),
      GIN("interp_order", "extrapolating"),
      GIN_TYPE("Index", "Numeric"),
      GIN_DEFAULT("1", "1e99"),
      GIN_DESC("Interpolation order (1=linear interpolation).",
               "Extrapolation allowed in interpolation of altitude.")));

  md_data_raw.push_back(MdRecord(
      NAME("MagRawRead"),
      DESCRIPTION(
          "Reads magnetic field data from a scenario.\n"
          "\n"
          "A full set of field components is read (NOTE: fails if scenario\n"
          "only contains selected field components). The files can be\n"
          "anywhere, but must all be in the same directory specified by\n"
          "'basename'. Naming convention for the field component files is\n"
          "basename.mag_u.xml for the u-component, v- and w-components\n"
          "accordingly.\n"),
      AUTHORS("Richard Larsson"),
      OUT("mag_u_field_raw", "mag_v_field_raw", "mag_w_field_raw"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("basename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Name of scenario, probably including the full path. For "
               "example: \"/data/magnetic_field\"")));

  md_data_raw.push_back(MdRecord(
      NAME("MatrixAddScalar"),
      DESCRIPTION(
          "Adds a scalar to all elements of a matrix.\n"
          "\n"
          "The result can either be stored in the same or another matrix.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Output matrix"),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Matrix", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input matrix.", "The value to be added to the matrix.")));

  md_data_raw.push_back(MdRecord(
      NAME("MatrixCBR"),
      DESCRIPTION(
          "Sets a matrix to hold cosmic background radiation (CBR).\n"
          "\n"
          "The CBR is assumed to be un-polarized and Stokes components 2-4\n"
          "are zero. Number of Stokes components, that equals the number\n"
          "of columns in the created matrix, is determined by *stokes_dim*.\n"
          "The number of rows in the created matrix equals the length of the\n"
          "given frequency vector.\n"
          "\n"
          "The cosmic radiation is modelled as blackbody radiation for the\n"
          "temperature given by the global constant COSMIC_BG_TEMP, set in\n"
          "the file constants.cc. The frequencies are taken from the generic\n"
          "input vector.\n"
          "\n"
          "The standard definition, in ARTS, of the Planck function is\n"
          "followed and the unit of the returned data is W/(m3 * Hz * sr).\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN("stokes_dim"),
      GIN("f"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Frequency vector.")));

  md_data_raw.push_back(MdRecord(
      NAME("MatrixCopySparse"),
      DESCRIPTION("Creates a matrix by copying a variable of type Sparse.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Created (full) matrix."),
      IN(),
      GIN("in"),
      GIN_TYPE("Sparse"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The sparse matrix to be copied.")));

  md_data_raw.push_back(MdRecord(
      NAME("MatrixExtractFromTensor3"),
      DESCRIPTION(
          "Extracts a Matrix from a Tensor3.\n"
          "\n"
          "Copies page or row or column with given Index from input Tensor3\n"
          "variable to output Matrix.\n"
          "Higher order equivalent of *VectorExtractFromMatrix*.\n"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Extracted matrix."),
      IN(),
      GIN("in", "i", "direction"),
      GIN_TYPE("Tensor3", "Index", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Input matrix.",
               "Index of page or row or column to extract.",
               "Direction. \"page\" or \"row\" or \"column\".")));

  md_data_raw.push_back(
      MdRecord(NAME("MatrixFromCovarianceMatrix"),
               DESCRIPTION("Turns a covariance matrix into a Matrix.\n"),
               AUTHORS("Richard Larsson"),
               OUT(),
               GOUT("out"),
               GOUT_TYPE("Matrix"),
               GOUT_DESC("Dense Matrix."),
               IN(),
               GIN("in"),
               GIN_TYPE("CovarianceMatrix"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Input covariance matrix.")));

  md_data_raw.push_back(MdRecord(
      NAME("MatrixIdentity"),
      DESCRIPTION(
          "Returns an identity matrix.\n"
          "\n"
          "The size if the matrix created is n x n. Default is to return a\n"
          "true identity matrix (I), but you can also select another value\n"
          "along the diagonal by setting *value*. That is, the output is\n"
          "value*I.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Output matrix"),
      IN(),
      GIN("n", "value"),
      GIN_TYPE("Index", "Numeric"),
      GIN_DEFAULT(NODEF, "1"),
      GIN_DESC("Size of the matrix", "The value along the diagonal.")));

  md_data_raw.push_back(MdRecord(
      NAME("MatrixMatrixMultiply"),
      DESCRIPTION(
          "Multiply a Matrix with another Matrix and store the result in the result\n"
          "Matrix.\n"
          "\n"
          "This just computes the normal Matrix-Matrix product, Y=M*X. It is ok\n"
          "if Y and X are the same Matrix.\n"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("The result of the multiplication (dimension mxc)."),
      IN(),
      GIN("m", "x"),
      GIN_TYPE("Matrix", "Matrix"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The Matrix to multiply (dimension mxn).",
               "The original Matrix (dimension nxc).")));

  md_data_raw.push_back(MdRecord(
      NAME("MatrixPlanck"),
      DESCRIPTION(
          "Sets a matrix to hold blackbody radiation.\n"
          "\n"
          "The radiation is assumed to be un-polarized and Stokes components\n"
          "2-4 are zero. Number of Stokes components, that equals the number\n"
          "of columns in the created matrix, is determined by *stokes_dim*.\n"
          "The number of rows in the created matrix equals the length of the\n"
          "given frequency vector.\n"
          "\n"
          "The standard definition, in ARTS, of the Planck function is\n"
          "followed and the unit of the returned data is W/(m3 * Hz * sr).\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN("stokes_dim"),
      GIN("f", "t"),
      GIN_TYPE("Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Frequency vector.", "Temperature [K].")));

  md_data_raw.push_back(MdRecord(
      NAME("MatrixScale"),
      DESCRIPTION("Scales all elements of a matrix with the specified value.\n"
                  "\n"
                  "The result can either be stored in the same or another\n"
                  "variable.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Output matrix"),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Matrix", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input matrix.",
               "The value to be multiplied with the matrix.")));

  md_data_raw.push_back(MdRecord(
      NAME("MatrixSet"),
      DESCRIPTION("Initialize a Matrix from the given list of numbers.\n"
                  "\n"
                  "Usage:\n"
                  "   MatrixSet(m1, [1, 2, 3; 4, 5, 6])\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("The newly created matrix"),
      IN(),
      GIN("value"),
      GIN_TYPE("Matrix"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The values of the newly created matrix. Elements are separated "
               "by commas, rows by semicolons."),
      SETMETHOD(true)));

  md_data_raw.push_back(MdRecord(
      NAME("MatrixSetConstant"),
      DESCRIPTION(
          "Creates a matrix and sets all elements to the specified value.\n"
          "\n"
          "The size is determined by *ncols* and *nrows*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN("nrows", "ncols"),
      GIN("value"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Matrix value.")));

  md_data_raw.push_back(MdRecord(
      NAME("MatrixUnitIntensity"),
      DESCRIPTION(
          "Sets a matrix to hold unpolarised radiation with unit intensity.\n"
          "\n"
          "Works as MatrixPlanck where the radiation is set to 1.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN("stokes_dim"),
      GIN("f"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Frequency vector.")));

  md_data_raw.push_back(MdRecord(
      NAME("MatrixVectorMultiply"),
      DESCRIPTION(
          "Multiply a Matrix with a Vector\n"
          "\n"
          "Computes the normal Matrix-Vector product, out=m*v. It is ok if out and v\n"
          "are the same Vector.\n"),
      AUTHORS("Stefan Buehler and Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("The result of the multiplication (length m)."),
      IN(),
      GIN("m", "v"),
      GIN_TYPE("Matrix", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The Matrix to multiply (dimension mxn).",
               "The original Vector (length n).")));

  md_data_raw.push_back(MdRecord(
      NAME("Matrix1ColFromVector"),
      DESCRIPTION("Forms a matrix containing one column from a vector.\n"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("v"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The vector to be copied.")));

  md_data_raw.push_back(MdRecord(
      NAME("Matrix2ColFromVectors"),
      DESCRIPTION(
          "Forms a matrix containing two columns from two vectors.\n"
          "\n"
          "The vectors are included as columns in the matrix in the same order\n"
          "as they are given.\n"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("v1", "v2"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The vector to be copied into the first column.",
               "The vector to be copied into the second column.")));

  md_data_raw.push_back(MdRecord(
      NAME("Matrix3ColFromVectors"),
      DESCRIPTION(
          "Forms a matrix containing three columns from three vectors.\n"
          "\n"
          "The vectors are included as columns in the matrix in the same order\n"
          "as they are given.\n"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("v1", "v2", "v3"),
      GIN_TYPE("Vector", "Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("The vector to be copied into the first column.",
               "The vector to be copied into the second column.",
               "The vector to be copied into the third column.")));

  md_data_raw.push_back(MdRecord(
      NAME("Matrix1RowFromVector"),
      DESCRIPTION("Forms a matrix containing one row from a vector.\n"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("v"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The vector to be copied.")));

  md_data_raw.push_back(MdRecord(
      NAME("Matrix2RowFromVectors"),
      DESCRIPTION(
          "Forms a matrix containing two rows from two vectors.\n"
          "\n"
          "The vectors are included as rows in the matrix in the same order\n"
          "as they are given.\n"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("v1", "v2"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The vector to be copied into the first row.",
               "The vector to be copied into the second row.")));

  md_data_raw.push_back(MdRecord(
      NAME("Matrix3RowFromVectors"),
      DESCRIPTION(
          "Forms a matrix containing three rows from three vectors.\n"
          "\n"
          "The vectors are included as rows in the matrix in the same order\n"
          "as they are given.\n"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("v1", "v2", "v3"),
      GIN_TYPE("Vector", "Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("The vector to be copied into the first row.",
               "The vector to be copied into the second row.",
               "The vector to be copied into the third row.")));

  md_data_raw.push_back(MdRecord(
      NAME("mblock_dlos_gridUniformCircular"),
      DESCRIPTION(
          "Gives *mblock_dlos_grid* roughly circular coverage, with uniform spacing.\n"
          "\n"
          "The method considers points on a regular grid with a spacing set by\n"
          "GIN *spacing*. All points inside *width* from (0,0) are included in\n"
          "*mblock_dlos_grid*. The positions in *mblock_dlos_grid* thus covers\n"
          "a roughly circular domain, and cover the same solid beam angle.\n"
          "Note that the method assumes that width is small and the solid beam\n"
          "angle does not change with distance from (0.0).\n"
          "\n"
          "Defualt is to consider grid positions of ..., -spacing/2, spacing/2, ...\n"
          "If you want to have (0,0) as a point in *mblock_dlos_grid*, change\n"
          "*centre* from its default value.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("mblock_dlos_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("spacing", "width", "centre"),
      GIN_TYPE("Numeric", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, "0"),
      GIN_DESC("The angular spacing between points.",
               "The maximum half-width to include.",
               "Set to 1 to place a point at (0,0).")));

  md_data_raw.push_back(MdRecord(
      NAME("mblock_dlos_gridUniformRectangular"),
      DESCRIPTION(
          "Gives *mblock_dlos_grid* rectangular coverage, with uniform spacing.\n"
          "\n"
          "The method creates an equidistant rectangular grid. The width is zenith\n"
          "and azimuth can differ. Note that selected widths are half-widths, and\n"
          "refers to the maximum value allowed. The actual width depends on values\n"
          "selected for *spacing* and *centre*.\n"
          "\n"
          "Defualt is to consider grid positions of ..., -spacing/2, spacing/2, ...\n"
          "If you want to have (0,0) as a point in *mblock_dlos_grid*, change\n"
          "*centre* from its default value.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("mblock_dlos_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("spacing", "za_width", "aa_width", "centre"),
      GIN_TYPE("Numeric", "Numeric", "Numeric", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "0"),
      GIN_DESC("The angular spacing between points.",
               "Max value of half-width in zenith angle direction.",
               "Max value of half-width in azimuth angle direction.",
               "Set to 1 to place a point at (0,0).")));

  md_data_raw.push_back(MdRecord(
      NAME("mc_antennaSetGaussian"),
      DESCRIPTION(
          "Makes mc_antenna (used by MCGeneral) a 2D Gaussian pattern.\n"
          "\n"
          "The gaussian antenna pattern is determined by *za_sigma* and\n"
          "*aa_sigma*, which represent the standard deviations in the\n"
          "uncorrelated bivariate normal distribution.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("mc_antennaSetGaussianByFWHM"),
      DESCRIPTION(
          "Makes mc_antenna (used by MCGeneral) a 2D Gaussian pattern.\n"
          "\n"
          "The gaussian antenna pattern is determined by *za_fwhm* and\n"
          "*aa_fwhm*, which represent the full width half maximum (FWHM)\n"
          "of the antenna response, in the zenith and azimuthal planes.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("mc_antennaSetPencilBeam"),
      DESCRIPTION(
          "Makes mc_antenna (used by MCGeneral) a pencil beam.\n"
          "\n"
          "This WSM makes the subsequent MCGeneral WSM perform pencil beam\n"
          "RT calculations.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("MCGeneral"),
      DESCRIPTION(
          "A generalised 3D reversed Monte Carlo radiative algorithm, that\n"
          "allows for 2D antenna patterns, surface reflection and arbitrary\n"
          "sensor positions.\n"
          "\n"
          "The main output variables *y* and *mc_error* represent the\n"
          "Stokes vector integrated over the antenna function, and the\n"
          "estimated error in this vector, respectively.\n"
          "\n"
          "The WSV *mc_max_iter* describes the maximum number of `photons\'\n"
          "used in the simulation (more photons means smaller *mc_error*).\n"
          "*mc_std_err* is the desired value of mc_error. *mc_max_time* is\n"
          "the maximum allowed number of seconds for MCGeneral. The method\n"
          "will terminate once any of the max_iter, std_err, max_time\n"
          "criteria are met. If negative values are given for these\n"
          "parameters then it is ignored.\n"
          "\n"
          "The WSV *mc_min_iter* sets the minimum number of photons to apply\n"
          "before the condition set by *mc_std_err* is considered. Values\n"
          "of *mc_min_iter* below 100 are not accepted.\n"
          "\n"
          "Only \"1\" and \"RJBT\" are allowed for *iy_unit*. The value of\n"
          "*mc_error* follows the selection for *iy_unit* (both for in- and\n"
          "output.\n"),
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
         "stokes_dim",
         "atmosphere_dim",
         "ppath_step_agenda",
         "ppath_lmax",
         "ppath_lraytrace",
         "iy_space_agenda",
         "surface_rtprop_agenda",
         "propmat_clearsky_agenda",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "z_field",
         "refellipsoid",
         "z_surface",
         "t_field",
         "vmr_field",
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

  md_data_raw.push_back(MdRecord(
      NAME("MCRadar"),
      DESCRIPTION(
          "A radar 3D foward Monte Carlo radiative algorithm, that allows \n"
          "for 2D antenna patterns and arbitrary sensor positions.\n"
          "Surface reflections are currently ignored.\n"
          "\n"
          "The main output variable *y* and *mc_error* represent the\n"
          "radar reflectivity integrated over the antenna function, and the\n"
          "estimated error in this vector, respectively.\n"
          "\n"
          "Unlike with yActive, the range bins gives the boundaries of \n"
          "the range bins as either round-trip time or distance from radar.\n"
          "\n"
          "The WSV *mc_y_tx* gives the polarization state of the \n"
          "transmitter.\n"
          "\n"
          "The WSV *mc_max_scatorder* prescribes the maximum scattering \n"
          "order to consider, after which `photon\'-tracing will be\n"
          "terminated. A value of one calculates only single scattering.\n"
          "\n"
          "The WSV *mc_max_iter* describes the maximum number of `photons\'\n"
          "used in the simulation (more photons means smaller *mc_error*).\n"
          "The method will terminate once the max_iter criterium is met.\n"
          "If negative values are given for these parameters then it is\n"
          "ignored.\n"
          "\n"
          "Only \"1\" and \"Ze\" are allowed for *iy_unit*. The value of\n"
          "*mc_error* follows the selection for *iy_unit* (both for in- and\n"
          "output.\n"),
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
         "stokes_dim",
         "atmosphere_dim",
         "ppath_lmax",
         "ppath_step_agenda",
         "ppath_lraytrace",
         "propmat_clearsky_agenda",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "z_field",
         "refellipsoid",
         "z_surface",
         "t_field",
         "vmr_field",
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
         "iy_unit",
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
      MdRecord(NAME("MCSetSeedFromTime"),
               DESCRIPTION("Sets the value of mc_seed from system time\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("nlte_fieldRescalePopulationLevels"),
      DESCRIPTION(
          "Rescale NLTE field to expected total distribution amongst levels\n"),
      AUTHORS("Richard Larsson"),
      OUT("nlte_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("nlte_field"),
      GIN("s"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Scaling (e.g., 0.75 for only orth-water on Earth)")));

  md_data_raw.push_back(MdRecord(
      NAME("nlte_fieldForSingleSpeciesNonOverlappingLines"),
      DESCRIPTION("NLTE field for a simple setup.\n"
        "\n"
        "This will solve for *nlte_field* in the input atmosphere.\n"
        "The solver depends on the lines not overlapping and that there\n"
        "is only a single species in the atmosphere.\n"
      ),
      AUTHORS("Richard Larsson"),
      OUT("nlte_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("nlte_field",
         "abs_species",
         "abs_lines_per_species",
         "nlte_level_identifiers",
         "collision_coefficients",
         "collision_line_identifiers",
         "isotopologue_ratios",
         "iy_main_agenda",
         "ppath_agenda",
         "iy_space_agenda",
         "iy_surface_agenda",
         "iy_cloudbox_agenda",
         "propmat_clearsky_agenda",
         "water_p_eq_agenda",
         "vmr_field",
         "t_field",
         "z_field",
         "p_grid",
         "atmosphere_dim",
         "refellipsoid",
         "surface_props_data",
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

  md_data_raw.push_back(MdRecord(
      NAME("collision_coefficientsFromSplitFiles"),
      DESCRIPTION(
          "Reads *collision_coefficients* and *collision_line_identifiers* from location on filesystem\n"
          "with many species.  The species in this location must match *abs_species*.  The location\n"
          "must also contain an ArrayOfQuantumIdentifier file ending with qid.xml\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("NumericAdd"),
      DESCRIPTION(
          "Adds a numeric and a value (out = in+value).\n"
          "\n"
          "The result can either be stored in the same or another numeric.\n"
          "(in and out can be the same varible, but not out and value)\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Output numeric."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input numeric.", "Value to add.")));

  md_data_raw.push_back(MdRecord(
      NAME("NumericFromVector"),
      DESCRIPTION(
          "Derivs a numeric from a vector, following selected operation.\n"
          "\n"
          "The following operations can be selected:\n"
          "  first : Selects the first element of the vector.\n"
          "   last : Selects the last element of the vector.\n"
          "    max : Selects the maximum element of the vector.\n"
          "    min : Selects the minimum element of the vector.\n"
          "   mean : Calculates the mean of the vector.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Output numeric."),
      IN(),
      GIN("in", "op"),
      GIN_TYPE("Vector", "String"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input vector.", "Selected operation.")));

  md_data_raw.push_back(MdRecord(
      NAME("NumericInvScale"),
      DESCRIPTION(
          "Inversely scales/divides a numeric with a value (out = in/value).\n"
          "\n"
          "The result can either be stored in the same or another numeric.\n"
          "(in and out can be the same varible, but not out and value)\n"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Output numeric."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input numeric.", "Scaling value.")));

  md_data_raw.push_back(MdRecord(
      NAME("NumericScale"),
      DESCRIPTION(
          "Scales/multiplies a numeric with a value (out = in*value).\n"
          "\n"
          "The result can either be stored in the same or another numeric.\n"
          "(in and out can be the same varible, but not out and value)\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Output numeric."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input numeric.", "Scaling value.")));

  md_data_raw.push_back(MdRecord(
      NAME("NumericSet"),
      DESCRIPTION("Sets a numeric workspace variable to the given value.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Numeric"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("value"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The value."),
      SETMETHOD(true)));

  md_data_raw.push_back(MdRecord(
      NAME("QuantumIdentifierSet"),
      DESCRIPTION(
          "Sets a QuantumIdentifier workspace variable to the given value\n"
          "by converting the input String\n"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("QuantumIdentifier"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("string_initializer"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The string representing the value."),
      SETMETHOD(true)));

  md_data_raw.push_back(MdRecord(
      NAME("ArrayOfQuantumIdentifierSet"),
      DESCRIPTION(
          "Sets an ArrayOfQuantumIdentifier workspace variable to the given value\n"
          "by converting the input ArrayOfString\n"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("ArrayOfQuantumIdentifier"),
      GOUT_DESC("Variables to initialize."),
      IN(),
      GIN("string_initializers"),
      GIN_TYPE("ArrayOfString"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The array of string representing the values of the array."),
      SETMETHOD(true)));

  md_data_raw.push_back(MdRecord(
      NAME("nelemGet"),
      DESCRIPTION(
          "Retrieve nelem from given variable and store the value in the\n"
          "variable *nelem*.\n"),
      AUTHORS("Oliver Lemke"),
      OUT("nelem"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("v"),
      GIN_TYPE(ARRAY_GROUPS + ", Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Variable to get the number of elements from."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(MdRecord(
      NAME("ncolsGet"),
      DESCRIPTION(
          "Retrieve ncols from given variable and store the value in the\n"
          "workspace variable *ncols*\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("nrowsGet"),
      DESCRIPTION(
          "Retrieve nrows from given variable and store the value in the\n"
          "workspace variable *nrows*\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("npagesGet"),
      DESCRIPTION(
          "Retrieve npages from given variable and store the value in the\n"
          "workspace variable *npages*\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("nbooksGet"),
      DESCRIPTION(
          "Retrieve nbooks from given variable and store the value in the\n"
          "workspace variable *nbooks*\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("nshelvesGet"),
      DESCRIPTION(
          "Retrieve nshelves from given variable and store the value in the\n"
          "workspace variable *nshelves*\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("nvitrinesGet"),
      DESCRIPTION(
          "Retrieve nvitrines from given variable and store the value in the\n"
          "workspace variable *nvitrines*\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("nlibrariesGet"),
      DESCRIPTION(
          "Retrieve nlibraries from given variable and store the value in the\n"
          "workspace variable *nlibraries*\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("nlte_sourceFromTemperatureAndSrcCoefPerSpecies"),
      DESCRIPTION(
          "Turn NLTE absorption per species into the source function by multiplying\n"
          "NLTE absorption per species with the LTE Planck source function.\n"),
      AUTHORS("Richard Larsson"),
      OUT("nlte_source", "dnlte_dx_source", "nlte_dsource_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("src_coef_per_species",
         "dsrc_coef_dx",
         "jacobian_quantities",
         "f_grid",
         "rtp_temperature"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(
      MdRecord(NAME("nlteOff"),
               DESCRIPTION("Disable Non-LTE calculations.\n"
                           "\n"
                           "The variables are set as follows:\n"
                           "   nlte_field             : Empty.\n"
                           "   nlte_level_identifiers : Empty.\n"),
               AUTHORS("Oliver Lemke"),
               OUT("nlte_do", "nlte_field", "nlte_level_identifiers"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN(),
               GIN_TYPE(),
               GIN_DEFAULT(),
               GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("nlteSetByQuantumIdentifiers"),
      DESCRIPTION(
          "Turns on NTLE calculations.\n"
          "\n"
          "Takes the quantum identifers for NLTE temperatures and matches it to\n"
          "lines in *abs_lines_per_species*.  *abs_species* must be set and is used\n"
          "to speed up calculations.  After the function is done,  all affected\n"
          "lines in *abs_lines_per_species* will have an internal tag to the relevant\n"
          "quantum identifier, which is a requirement for deeper code.\n"
          "\n"
          "If vibrational_energies is input it must match *nlte_level_identifiers*\n"
          "in length.  The vibrational energies of the affected lines will then be\n"
          "set by the function.  Otherwise, it is assumed the vibrational energies\n"
          "are set by another method.  If they are not set, calculations will complain\n"
          "later on while running deeper code.\n"
          "\n"
          "For now only vibrational energy states are assumed to be able to be in\n"
          "non-LTE conditions.  The *QuantumIdentifier* for an energy state in ARTS\n"
          "can look like:\n"
          "\t\"CO2-626 EN v1 0/1 v2 1/1 l2 1/1 v3 0/1 r 1/1\"\n"
          "and the matching will match ALL lines with the above.  Note then that if, e.g.,\n"
          "the \"v1 0/1\" term was removed from the above, then ARTS will assume that\n"
          "\"v1\" is not part of the level of energy state of interest, so lines\n"
          "of different \"v1\" will be matched as the same state.  If a line is matched\n"
          "to more than one energy state, errors should be thrown, but be careful.\n"
          "\n"
          "Set type of population to change computations and expected input as:\n"
          "\tLTE: Compute population by ratios found from LTE temperatures\n"
          "\tTV:  Compute population by ratios found from NLTE vibrational temperatures\n"
          "\tND:  Compute population by ratios found from NLTE number densities\n"),
      AUTHORS("Richard Larsson"),
      OUT("nlte_do", "abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species", "nlte_level_identifiers", "abs_species"),
      GIN("vibrational_energies", "population_type"),
      GIN_TYPE("Vector", "String"),
      GIN_DEFAULT("[]", "TV"),
      GIN_DESC(
          "Vector of vibrational energies.  If empty, assume known or unimportant vibrational energies.",
          "String for setting the type of population.")));

  md_data_raw.push_back(MdRecord(
      NAME("nlte_fieldSetLteExternalPartitionFunction"),
      DESCRIPTION("Turns on NTLE calculations.\n"
                  "\n"
                  "Sets NLTE ratios to those expected for LTE calculations\n"
                  "with a known partition function\n"),
      AUTHORS("Richard Larsson"),
      OUT("nlte_do", "nlte_field", "abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species",
         "nlte_level_identifiers",
         "partition_functions",
         "t_field"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("ArrayOfQuantumIdentifierFromLines"),
      DESCRIPTION(
          "Sets an ArrayOfQuantumIdentifier to all levels in *abs_lines_per_species*\n"
          "with defined quantum numbers\n"
          "\n"
          "Lines without defined quantum numbers are ignored\n"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("ArrayOfQuantumIdentifier"),
      GOUT_DESC("Identifiers to all levels in *abs_lines_per_species*"),
      IN("abs_lines_per_species"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("nlte_fieldSetLteInternalPartitionFunction"),
      DESCRIPTION(
          "Turns on NTLE calculations.\n"
          "\n"
          "Sets NLTE ratios to those expected for LTE calculations\n"
          "with estimation of the partition function as the sum of all\n"
          "states of a species\n"),
      AUTHORS("Richard Larsson"),
      OUT("nlte_do", "nlte_field", "abs_lines_per_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_species", "nlte_level_identifiers", "t_field"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("OEM"),
      DESCRIPTION(
          "Inversion by the so called optimal estimation method (OEM).\n"
          "\n"
          "Work in progress ...\n"
          "\n"
          "The cost function to minimise, including a normalisation with length"
          "of *y*, is:\n"
          "   cost = cost_y + cost_x\n"
          "where\n"
          "   cost_y = 1/m * [y-yf]' * covmat_se_inv * [y-yf]\n"
          "   cost_x = 1/m * [x-xa]' * covmat_sx_inv * [x-xa]\n"
          "\n"
          " The current implementation provides 3 methods for the minimization of\n"
          "the cost functional: Linear, Gauss-Newton and Levenberg-Marquardt.\n"
          "The Gauss-Newton minimizer attempts to find a minimum solution by \n"
          "fitting a quadratic function to the cost functional. The linear minimizer\n"
          "is a special case of the Gauss-Newton method, since for a linear forward\n"
          "model the exact solution of the minimization problem is obtained after\n"
          "the first step. The Levenberg-Marquardt method adaptively constrains the\n"
          "search region for the next iteration step by means of the so-called gamma-factor.\n"
          "This makes the method more suitable for strongly non-linear problems.\n"
          "If the gamma-factor is 0, Levenberg-Marquardt and Gauss-Newton method\n"
          "are identical. Each minimization method (li,gn,lm) has an indirect\n"
          "variant (li_cg,gn_cg,lm_cg), which uses the conjugate gradient solver\n"
          "for the linear system that has to be solved in each minimzation step.\n"
          "This of advantage for very large problems, that would otherwise require\n"
          "the computation of expensive matrix products.\n"
          "\n"
          "Description of the special input arguments:\n"
          "\n"
          "*method*\n"
          "  \"li\": A linear problem is assumed and a single iteration is performed.\n"
          "  \"li_cg\": A linear problem is assumed and solved using the CG solver.\n"
          "  \"gn\": Non-linear, with Gauss-Newton iteration scheme.\n"
          "  \"gn_cg\": Non-linear, with Gauss-Newton and conjugate gradient solver.\n"
          "  \"lm\": Non-linear, with Levenberg-Marquardt (LM) iteration scheme.\n"
          "  \"lm_cg\": Non-linear, with Levenberg-Marquardt (LM) iteration scheme and conjugate gradient solver.\n"
          "*max_start_cost*\n"
          "  No inversion is done if the cost matching the a priori state is above\n"
          "  this value. If set to a negative value, all values are accepted.\n"
          "  This argument also controls if the start cost is calculated. If\n"
          "  set to <= 0, the start cost in *oem_diagnostics* is set to NaN\n"
          "  when using \"li\" and \"gn\".\n"
          "*x_norm*\n"
          "  A normalisation vector for *x*. A normalisation of *x* can be needed\n"
          "  due to limited numerical precision. If this vector is set to be empty\n"
          "  no normalisation is done (defualt case). Otherwise, this must be a\n"
          "  vector with same length as *x*, just having values above zero.\n"
          "  Elementwise division between *x* and *x_norm* (x./x_norm) shall give\n"
          "  a vector where all values are in the order of unity. Maybe the best\n"
          "  way to set *x_norm* is x_norm = sqrt( diag( Sx ) ).\n"
          "*max_iter*\n"
          "  Maximum number of iterations to perform. No effect for \"li\".\n"
          "*stop_dx*\n"
          "  Iteration stop criterion. The criterion used is the same as given\n"
          "  in Rodgers\' \"Inverse Methods for Atmospheric Sounding\"\n"
          "*lm_ga_settings*\n"
          "  Settings controlling the gamma factor, part of the \"LM\" method.\n"
          "  This is a vector of length 6, having the elements (0-based index):\n"
          "    0: Start value.\n"
          "    1: Fractional decrease after succesfull iteration.\n"
          "    2: Fractional increase after unsuccessful iteration.\n"
          "    3: Maximum allowed value. If the value is passed, the inversion\n"
          "       is halted.\n"
          "    4: Lower treshold. If the threshold is passed, gamma is set to zero.\n"
          "       If gamma must be increased from zero, gamma is set to this value.\n"
          "    5: Gamma limit. This is an additional stop criterion. Convergence\n"
          "       is not considered until there has been one succesful iteration\n"
          "       having a gamma <= this value.\n"
          "  The default setting triggers an error if \"lm\" is selected.\n"
          "*clear matrices*\n"
          "   With this flag set to 1, *jacobian* and *dxdy* are returned as empty\n"
          "   matrices.\n"
          "*display_progress*\n"
          "   Controls if there is any screen output. The overall report level\n"
          "   is ignored by this WSM.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("avkCalc"),
      DESCRIPTION(
          "Calculates the averaging kernel matrix describing the sensitivity of the\n"
          "OEM retrieval with respect to the true state of the system. A prerequisite\n"
          "for the calculation of the averaging kernel matrix is a successful OEM\n"
          "calculation in which the jacobian and the gain matrix dxdy have been calculated.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("covmat_soCalc"),
      DESCRIPTION(
          "Calculates the covariance matrix describing the error due to uncertainties\n"
          "in the observation system. The uncertainties of the observation system are\n"
          "described by *covmat_se*, which must be set by the user to include the\n"
          "relevant contributions from the measurement and the forward model.\n"
          "\n"
          "Prerequisite for the calculation of *covmat_so* is a successful OEM\n"
          "computation where also the gain matrix has been computed.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("covmat_ssCalc"),
      DESCRIPTION(
          "Calculates the covariance matrix describing the error due to smoothing.\n"
          ""
          "The calculation of *covmat_ss* also requires the averaging kernel matrix *avk*\n"
          "to be computed after a successful OEM calculation.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("opt_prop_bulkCalc"),
      DESCRIPTION(
          "Calculates bulk absorption extinction at one atmospheric grid point.\n"
          "\n"
          "This WSM sums up the monochromatic absorption vectors and\n"
          "extinction matrices of all scattering elements (*abs_vec_spt* and\n"
          "*ext_mat_spt*, respectively) weighted by their respective\n"
          "particle number density given by *pnd_field*, for a single location\n"
          "within the cloudbox, given by *scat_p_index*, *scat_lat_index*, and\n"
          "*scat_lon_index*.\n"
          "The resulting  extinction matrix is added to the workspace variable\n"
          "*ext_mat*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("opt_prop_sptFromData"),
      DESCRIPTION(
          "Calculates monochromatic optical properties for all scattering\n"
          "elements.\n"
          "\n"
          "In this function the extinction matrix and the absorption vector\n"
          "are calculated in the laboratory frame. An interpolation of the\n"
          "data on the actual frequency is the first step in this function.\n"
          "The next step is a transformation from the database coordinate\n"
          "system to the laboratory coordinate system.\n"
          "\n"
          "Output of the function are *ext_mat_spt* and *abs_vec_spt*, which\n"
          "hold the optical properties for a specified propagation direction\n"
          "for each scattering element.\n"),
      AUTHORS("Claudia Emde"),
      OUT("ext_mat_spt", "abs_vec_spt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ext_mat_spt",
         "abs_vec_spt",
         "scat_data",
         "scat_za_grid",
         "scat_aa_grid",
         "scat_za_index",
         "scat_aa_index",
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

  md_data_raw.push_back(MdRecord(
      NAME("opt_prop_sptFromScat_data"),
      DESCRIPTION(
          "Derives monochromatic optical properties for all scattering\n"
          "elements.\n"
          "\n"
          "As *opt_prop_sptFromData*, but using frequency pre-interpolated\n"
          "data (as produced by *scat_dataCalc*), i.e. in here no frequency\n"
          "interpolation is done anymore.\n"),
      AUTHORS("Jana Mendrok, Claudia Emde"),
      OUT("ext_mat_spt", "abs_vec_spt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ext_mat_spt",
         "abs_vec_spt",
         "scat_data",
         "scat_data_checked",
         "scat_za_grid",
         "scat_aa_grid",
         "scat_za_index",
         "scat_aa_index",
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

  md_data_raw.push_back(MdRecord(
      NAME("opt_prop_sptFromMonoData"),
      DESCRIPTION(
          "Calculates optical properties for the scattering elements.\n"
          "\n"
          "As *opt_prop_sptFromData* but no frequency interpolation is\n"
          "performed. The single scattering data is here obtained from\n"
          "*scat_data_mono*, instead of *scat_data*.\n"),
      AUTHORS("Cory Davis"),
      OUT("ext_mat_spt", "abs_vec_spt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ext_mat_spt",
         "abs_vec_spt",
         "scat_data_mono",
         "scat_za_grid",
         "scat_aa_grid",
         "scat_za_index",
         "scat_aa_index",
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
      MdRecord(NAME("output_file_formatSetAscii"),
               DESCRIPTION("Sets the output file format to ASCII.\n"),
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
      MdRecord(NAME("output_file_formatSetBinary"),
               DESCRIPTION("Sets the output file format to binary.\n"),
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
      MdRecord(NAME("output_file_formatSetZippedAscii"),
               DESCRIPTION("Sets the output file format to zipped ASCII.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("particle_bulkprop_fieldClip"),
      DESCRIPTION(
          "Clipping of *particle_bulkprop_field*.\n"
          "\n"
          "The method allows you to apply hard limits the values of\n"
          "*particle_bulkprop_field*. All values, of the property selected,\n"
          "below *limit_low*, are simply set to *limit_low*. And the same\n"
          "is performed with respect to *limit_high*. That is, the data in x\n"
          "for the retrieval quantity are forced to be inside the range\n"
          "[limit_low,limit_high].\n"
          "\n"
          "Setting species=\"ALL\", is a shortcut for applying the limits on all\n"
          "properties.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("particle_bulkprop_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("particle_bulkprop_field", "particle_bulkprop_names"),
      GIN("bulkprop_name", "limit_low", "limit_high"),
      GIN_TYPE("String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "-Inf", "Inf"),
      GIN_DESC("Name of bulk property to consider, or \"ALL\".",
               "Lower limit for clipping.",
               "Upper limit for clipping.")));

  md_data_raw.push_back(MdRecord(
      NAME("particle_bulkprop_fieldPerturb"),
      DESCRIPTION(
          "Adds a perturbation to *particle_bulkprop_field*.\n"
          "\n"
          "Works as *AtmFieldPerturb* but acts on *particle_bulkprop_field*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("particle_bulkprop_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("particle_bulkprop_field","atmosphere_dim",
         "p_grid","lat_grid","lon_grid","particle_bulkprop_names"),
      GIN("particle_type","p_ret_grid","lat_ret_grid","lon_ret_grid",
          "pert_index","pert_size","pert_mode"),
      GIN_TYPE("String","Vector","Vector","Vector",
               "Index","Numeric","String"),
      GIN_DEFAULT(NODEF,NODEF,NODEF,NODEF,NODEF,NODEF,"absolute"),
      GIN_DESC("Name of field to perturb, such as ""IWC"".",
               "Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retrieval grid.",
               "Index of position where the perturbation shall be performed.",
               "Size of perturbation.",
               "Type of perturbation, ""ansolute"" or ""relative"".")));

  md_data_raw.push_back(MdRecord(
      NAME("particle_bulkprop_fieldPerturbAtmGrids"),
      DESCRIPTION(
          "Adds a perturbation to *particle_bulkprop_field*.\n"
          "\n"
          "Works as *AtmFieldPerturbAtmGrids* but acts on *particle_bulkprop_field*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("particle_bulkprop_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("particle_bulkprop_field","atmosphere_dim",
         "p_grid","lat_grid","lon_grid","particle_bulkprop_names"),
      GIN("particle_type","pert_index","pert_size","pert_mode"),
      GIN_TYPE("String","Index","Numeric","String"),
      GIN_DEFAULT(NODEF,NODEF,NODEF,"absolute"),
      GIN_DESC("Name of field to perturb, such as ""IWC"".",
               "Index of position where the perturbation shall be performed.",
               "Size of perturbation.",
               "Type of perturbation, ""ansolute"" or ""relative"".")));

  md_data_raw.push_back(MdRecord(
      NAME("particle_massesFromMetaDataSingleCategory"),
      DESCRIPTION(
          "Sets *particle_masses* based on *scat_meta* assuming\n"
          "all particles are of the same mass category.\n"
          "\n"
          "This method derives the particle masses from the mass entry\n"
          "of each scattering element. It is assumed that all scattering\n"
          "elements represent particles of the same (bulk) matter\n"
          "(e.g. water or ice). With other words, a single mass category\n"
          "is assumed (see *particle_masses* for a definition of \"mass\n"
          "category\").\n"
          "\n"
          "To be clear, the above are assumptions of the method, the user\n"
          "is free to work with any scattering element. For Earth and just having\n"
          "cloud and particles, the resulting mass category can be seen as\n"
          "the total cloud water content, with possible contribution from\n"
          "both ice and liquid phase.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("particle_massesFromMetaData"),
      DESCRIPTION(
          "Derives *particle_masses* from *scat_meta*.\n"
          "\n"
          "This method is supposed to be used to derive *particle_masses*\n"
          "when *pnd_field* is internally calculated using\n"
          "*pnd_fieldCalcFromscat_speciesFields* (in contrast to reading it\n"
          "from external sources using *ScatElementsPndAndScatAdd* and\n"
          "*pnd_fieldCalcFrompnd_field_raw*).\n"
          "It extracts the mass information of the scattering elements from\n"
          "*scat_meta*. Each scattering species is taken as a separate\n"
          "category of particle_masses, i.e., the resulting\n"
          "*particle_masses* matrix will contain as many columns as\n"
          "scattering species are present in *scat_meta*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("partition_functionsInitFromBuiltin"),
      DESCRIPTION(
          "Initialize partition functions with default values from built-in\n"
          "species data.\n"),
      AUTHORS("Oliver Lemke"),
      OUT("partition_functions"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("pha_matCalc"),
      DESCRIPTION(
          "Calculates the total phase matrix of all scattering elements.\n"
          "\n"
          "This function sums up the monochromatic phase matrices of all\n"
          "scattering elements *pha_mat_spt* weighted with  their respective\n"
          "particle number density, given by *pnd_field*, for a single location\n"
          "within the cloudbox, given by *scat_p_index*, *scat_lat_index*, and\n"
          "*scat_lon_index*.\n"),
      AUTHORS("Sreerekha T.R."),
      OUT("pha_mat"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pha_mat_spt",
         "pnd_field",
         "atmosphere_dim",
         "scat_p_index",
         "scat_lat_index",
         "scat_lon_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("pha_mat_sptFromData"),
      DESCRIPTION(
          "Calculation of the phase matrix of the individual scattering elements.\n"
          "\n"
          "This function can be used in *pha_mat_spt_agenda* as part of\n"
          "the calculation of the scattering integral.\n"
          "\n"
          "First, data at the requested frequency (given by *f_grid* and\n"
          "*f_index*) and temperature (given by *rtp_temperature*) is\n"
          "extracted. This is followed by a transformation from the database\n"
          "coordinate system to the laboratory coordinate system.\n"
          "\n"
          "Frequency extraction is always done by (linear) interpolation.\n"
          "Temperature is (linearly) interpolated when at least two\n"
          "temperature grid points are present in the *scat_data* and\n"
          "*rtp_temperature* is positive. If only a single temperature point\n"
          "is available, data for this point is used without modification. In\n"
          "order to speed up calculations, temperature interpolation can be\n"
          "avoided by passing a *rtp_temperature*<0. In this case, a specific\n"
          "temperature grid from the *scat_data* grid is used without\n"
          "modification. The selection is as follows:\n"
          "  -10 < *rtp_temperature * <   0   T_grid[0]     lowest temperature\n"
          "  -20 < *rtp_temperature * < -10   T_grid[nT-1]  highest temperature\n"
          "        *rtp_temperature*  < -20   T_grid[nT/2]  median grid point\n"),
      AUTHORS("Claudia Emde"),
      OUT("pha_mat_spt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pha_mat_spt",
         "scat_data",
         "scat_za_grid",
         "scat_aa_grid",
         "scat_za_index",
         "scat_aa_index",
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

  md_data_raw.push_back(MdRecord(
      NAME("pha_mat_sptFromDataDOITOpt"),
      DESCRIPTION(
          "Calculation of the phase matrix of the individual scattering elements.\n"
          "\n"
          "In this function the phase matrix is extracted from\n"
          "*pha_mat_sptDOITOpt*. It can be used in the agenda\n"
          "*pha_mat_spt_agenda*. This method must be used in combination with\n"
          "*DoitScatteringDataPrepare*.\n"
          "\n"
          "Temperature is considered as described for *pha_mat_sptFromData*\n"),
      AUTHORS("Claudia Emde"),
      OUT("pha_mat_spt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pha_mat_spt",
         "pha_mat_sptDOITOpt",
         "scat_data_mono",
         "doit_za_grid_size",
         "scat_aa_grid",
         "scat_za_index",
         "scat_aa_index",
         "rtp_temperature",
         "pnd_field",
         "scat_p_index",
         "scat_lat_index",
         "scat_lon_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("pha_mat_sptFromMonoData"),
      DESCRIPTION(
          "Calculation of the phase matrix of the individual scattering elements.\n"
          "\n"
          "This function is the monochromatic version of *pha_mat_sptFromData*.\n"),
      AUTHORS("Claudia Emde"),
      OUT("pha_mat_spt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pha_mat_spt",
         "scat_data_mono",
         "doit_za_grid_size",
         "scat_aa_grid",
         "scat_za_index",
         "scat_aa_index",
         "rtp_temperature",
         "pnd_field",
         "scat_p_index",
         "scat_lat_index",
         "scat_lon_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("pha_mat_sptFromScat_data"),
      DESCRIPTION(
          "Calculation of the phase matrix of the individual scattering elements.\n"
          "\n"
          "As *pha_mat_sptFromData*, but using frequency pre-interpolated\n"
          "data (as produced by *scat_dataCalc*), i.e. in here no frequency\n"
          "interpolation is done anymore.\n"),
      AUTHORS("Jana Mendrok, Claudia Emde"),
      OUT("pha_mat_spt"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pha_mat_spt",
         "scat_data",
         "scat_data_checked",
         "scat_za_grid",
         "scat_aa_grid",
         "scat_za_index",
         "scat_aa_index",
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

  md_data_raw.push_back(MdRecord(
      NAME("pndAdjustFromScatMeta"),
      DESCRIPTION("Rescales *pnd_data* to match expected material content.\n"
                  "\n"
                  "A temporary method...\n"
                  "Reproduces rescaling (using chk_pndsum) applied in\n"
                  "*pnd_fieldCalcFromscat_speciesFields*\n"),
      AUTHORS("Jana Mendrok"),
      OUT("pnd_data"),  //, "dpnd_data_dx" ),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pnd_data",  // "dpnd_data_dx",
         "pnd_agenda_input",
         "pnd_agenda_input_names",
         "scat_meta"),
      GIN("pnd_agenda_input_tag", "scat_index"),
      GIN_TYPE("String", "Index"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Name of field in *pnd_agenda_input* to which to adjust"
               " pnd_data to.",
               "Take data from scattering species of this index (0-based) in"
               " *scat_meta*.")));

  md_data_raw.push_back(MdRecord(
      NAME("pndFromdNdD"),
      DESCRIPTION(
          "Calculates pnds from given dNdD.\n"
          "\n"
          "The method mimics what happens inside\n"
          "*pnd_fieldCalcFromscat_speciesFields*, but for a single size\n"
          "distribution. It is supposed to be used with the *dNdD* methods.\n"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("pnd"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("The pnd vector (pnd as function of particle size)"),
      IN(),
      GIN("dNdD", "diameter", "total_content", "scatelem_content"),
      GIN_TYPE("Vector", "Vector", "Numeric", "Vector"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, NODEF),
      GIN_DESC("Size distribution number density",
               "Size of the particles (the same as used in deriving dNdD",
               "Material content that should be contained in the"
               "distribution. E.g., Mass density, mass flux, total number"
               "density. If dNdD was derived from a content-dependent PSD,"
               "then this value should correspond to the one used there.",
               "Material content of each particle (scattering element)."
               "Needs to be the same content type as the total content"
               "above.")));

  md_data_raw.push_back(MdRecord(
      NAME("pndFromPsd"),
      DESCRIPTION(
          "Calculates *pnd_data* from given *psd_data* for one scattering species.\n"
          "\n"
          "Performs integration of the size distribution over the size grid\n"
          "bin deriving pnd (units #/m3) from psd (units #/m3/m). Some checks\n"
          "on the sufficiency of the size grid range and coverage are applied.\n"
          "\n"
          "*quad_order* can be 0 for rectangular or 1 for trapezoidal\n"
          "integration. The only difference is the treatment of the start and\n"
          "end nodes. For trapezoidal their corresponding bins end exactly at\n"
          "the nodes, while for rectangular they extend further out by the half\n"
          "distance to the neighbor node (but not beyond 0).\n"
          "\n"
          "Attempts to check that the size grids and *scat_data* represent the\n"
          "bulk extinction sufficiently. Specifically, it is tested that\n"
          " (a) psd*ext is decreasing at the small and large particle size\n"
          "     ends of the size grid - but only if scattering species bulk\n"
          "     extinction exceeds 1% of *threshold_ss_ext*.\n"
          " (b) removing the smallest and largest particles changes the\n"
          "     resulting bulk extinction by less then a fraction of\n"
          "     *threshold_se_ext* - but only if scattering species bulk\n"
          "     extinction exceeds *threshold_ss_ext* and number density (pnd)\n"
          "     of the edge size point at this atmospheric level is larger\n"
          "     than *threshold_se_pnd* times the maximum pnd of this\n"
          "     scattering element over all atmospheric levels.\n"
          "Skipping tests in case of low extinction is done in order to\n"
          "minimize issues arising from very low mass densities,\n"
          "particularly at single atmospheric levels, and very low bulk\n"
          "extinctions, i.e. in cases where the effects on the radiance fields\n"
          "are estimated to be low."
          "\n"
          "NOTE: The tests are only approximate and do not guarantee the\n"
          "validity of the resulting bulk properties (and increasing the\n"
          "thresholds will decrease the reliability of the bulk properties).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("pndFromPsdBasic"),
      DESCRIPTION(
          "Calculates *pnd_data* from given *psd_data*.\n"
          "\n"
          "As *pndFromPsdBasic*, but without bulk extinction representation\n"
          "checks.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("pnd_fieldCalcFromParticleBulkProps"),
      DESCRIPTION(
          "Converts particle bulk property data to *pnd_field*.\n"
          "\n"
          "In short, the method combines *scat_species*, *pnd_agenda_array*,\n"
          "*particle_bulkprop_field* and their associated variables to derive\n"
          "*pnd_field*.\n"
          "\n"
          "The method does nothing if cloudbox is inactive.\n"
          "\n"
          "Otherwise, cloudbox limits must be set before calling the method,\n"
          "and *particle_bulkprop_field* is checked to have non-zero elements\n"
          "just inside the cloudbox.\n"),
      AUTHORS("Patrick Eriksson, Jana Mendrok"),
      OUT("pnd_field", "dpnd_field_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "t_field",
         "cloudbox_on",
         "cloudbox_limits",
         "scat_species",
         "scat_data",
         "scat_meta",
         "particle_bulkprop_field",
         "particle_bulkprop_names",
         "pnd_agenda_array",
         "pnd_agenda_array_input_names",
         "jacobian_do",
         "jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("pnd_fieldCalcFrompnd_field_raw"),
      DESCRIPTION(
          "Interpolation of particle number density fields to calculation grid\n"
          "inside cloudbox.\n"
          "\n"
          "This method interpolates the particle number density field\n"
          "from the raw data *pnd_field_raw* to obtain *pnd_field*.\n"
          "For 1D cases, where internally *GriddedFieldPRegrid* and\n"
          "*GriddedFieldLatLonRegrid* are applied, *zeropadding*=1 sets the\n"
          "*pnd_field* at pressure levels levels exceeding pnd_field_raw's\n"
          "pressure grid to 0 (not implemented for 2D and 3D yet). Default:\n"
          "zeropadding=0, which throws an error if the calculation pressure grid\n"
          "*p_grid* is not completely covered by pnd_field_raw's pressure grid.\n"),
      AUTHORS("Sreerekha T.R.", "Claudia Emde", "Oliver Lemke"),
      OUT("pnd_field", "dpnd_field_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("p_grid",
         "lat_grid",
         "lon_grid",
         "pnd_field_raw",
         "atmosphere_dim",
         "cloudbox_limits",
         "jacobian_quantities"),
      GIN("zeropadding"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("Allow zeropadding of pnd_field.")));

  md_data_raw.push_back(MdRecord(
      NAME("pnd_fieldCalcFromscat_speciesFields"),
      DESCRIPTION(
          "Calculation of *pnd_field* from *scat_meta* and\n"
          "the scattering species fields.\n"
          "\n"
          "The method calculates the number densities (pnd_field values) of\n"
          "all individual scattering elements covered by the considered\n"
          "scattering species for all grid points in the cloudbox. The pnds\n"
          "represent the mass density, mass flux, an/or number density fields\n"
          "of the scattering species according to the particle size\n"
          "distribution (PSD) chosen for the respective scattering species.\n"
          "\n"
          "The following PSDs are available (for further information,\n"
          "including parametrization limits and internal assumptions, see\n"
          "their corresponding PSD WSM):\n"
          "\n"
          "Tag        PSD WSM          fields(s) used           Target         Notes\n"
          "MH97       *dNdD_MH97*      mass density             cloud ice\n"
          "H11        *dNdD_H11*       mass density             cloud ice\n"
          "H13        *dNdD_H13_Ar*    mass density             cloud ice      neglects shape information\n"
          "H13Shape   *dNdD_H13_Ar*    mass density             cloud ice\n"
          "MGD_IWC    *dNdD_MGD_IWC*   mass density             cloud ice      fixed modified gamma psd\n"
          "SB06_IWC   *dNdD_SB06*      number & mass density    cloud ice      two moment scheme psd\n"
          "SB06_IWC_M *dNdD_SB06_M*    mean mass & mass density cloud ice      two moment scheme psd\n"
          "MY05_IWC   *dNdD_MY05*      number & mass density    cloud ice      two moment scheme psd\n"
          "MY05_IWC_M *dNdD_MY05_M*    mean mass & mass density cloud ice      two moment scheme psd\n"
          "F07TR      *dNdD_F07*       mass density             snow           for tropics\n"
          "F07ML      *dNdD_F07*       mass density             snow           for midlatitudes\n"
          "SB06_SWC   *dNdD_SB06*      number & mass density    snow           two moment scheme psd\n"
          "SB06_SWC_M *dNdD_SB06_M*    mean mass & mass density snow           two moment scheme psd\n"
          "MY05_SWC   *dNdD_MY05*      number & mass density    snow           two moment scheme psd\n"
          "MY05_SWC_M *dNdD_MY05_M*    mean mass & mass density snow           two moment scheme psd\n"
          "MP48       *dNdD_MP48*      mass flux                precipitation  rain in particular\n"
          "W16        *dNdD_W16*       mass density             rain\n"
          "SB06_RWC   *dNdD_SB06*      number & mass density    rain           two moment scheme psd\n"
          "SB06_RWC_M *dNdD_SB06_M*    mean mass & mass density rain           two moment scheme psd\n"
          "MY05_RWC   *dNdD_MY05*      number & mass density    rain           two moment scheme psd\n"
          "MY05_RWC_M *dNdD_MY05_M*    mean mass & mass density rain           two moment scheme psd\n"
          "H98_STCO   *dNdD_H98*       mass density             cloud liquid   specifically continental stratus\n"
          "MGD_LWC    *dNdD_MGD_LWC*   mass density             cloud liquid   fixed modified gamma psd\n"
          "SB06_LWC   *dNdD_SB06*      number & mass density    cloud liquid   two moment scheme psd\n"
          "SB06_LWC_M *dNdD_SB06_M*    mean mass & mass density cloud liquid   two moment scheme psd\n"
          "MY05_LWC   *dNdD_MY05*      number & mass density    cloud liquid   two moment scheme psd\n"
          "MY05_LWC_M *dNdD_MY05_M*    mean mass & mass density cloud liquid   two moment scheme psd\n"
          "SB06_GWC   *dNdD_SB06*      number & mass density    graupel        two moment scheme psd\n"
          "SB06_GWC_M *dNdD_SB06_M*    mean mass & mass density graupel        two moment scheme psd\n"
          "MY05_GWC   *dNdD_MY05*      number & mass density    graupel        two moment scheme psd\n"
          "MY05_GWC_M *dNdD_MY05_M*    mean mass & mass density graupel        two moment scheme psd\n"
          "SB06_HWC   *dNdD_SB06*      number & mass density    hail           two moment scheme psd\n"
          "SB06_HWC_M *dNdD_SB06_M*    mean mass & mass density hail           two moment scheme psd\n"
          "MY05_HWC   *dNdD_MY05*      number & mass density    hail           two moment scheme psd\n"
          "MY05_HWC_M *dNdD_MY05_M*    mean mass & mass density hail           two moment scheme psd\n"
          "\n"
          "NOTE: The number and order of the scattering species in the\n"
          "scattering species fields (*scat_species_mass_density_field*,\n"
          "*scat_species_mass_flux_field*, *scat_species_number_density_field*,\n"
          "*scat_species_mean_mass_field*) has to fit number and order of the\n"
          "*scat_species* tags.\n"
          "Moreover, the order of *scat_species* tags has to fit the order of\n"
          "scattering species in the *scat_meta* array, i.e.,\n"
          "*ScatSpeciesScatAndMetaRead* with the respective scattering\n"
          "data and meta data files has to be applied in the right order!\n"),
      AUTHORS("Daniel Kreyling, Jana Mendrok, Manfred Brath"),
      OUT("pnd_field", "dpnd_field_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "cloudbox_on",
         "cloudbox_limits",
         "scat_species_mass_density_field",
         "scat_species_mass_flux_field",
         "scat_species_number_density_field",
         "scat_species_mean_mass_field",
         "t_field",
         "scat_meta",
         "scat_species",
         "jacobian_quantities"),
      GIN("delim"),
      GIN_TYPE("String"),
      GIN_DEFAULT("-"),
      GIN_DESC("Delimiter string of *scat_species* elements.")));

  md_data_raw.push_back(MdRecord(
      NAME("pnd_fieldExpand1D"),
      DESCRIPTION(
          "Maps a 1D pnd_field to a (homogeneous) 2D or 3D pnd_field.\n"
          "\n"
          "This method takes a 1D *pnd_field* and converts it to a 2D or 3D\n"
          "\"cloud\". It is assumed that a complete 1D case has been created,\n"
          "and after this *atmosphere_dim*, *lat_grid*, *lon_grid* and\n"
          "*cloudbox_limits* have been changed to a 2D or 3D case (without\n"
          "changing the vertical extent of the cloudbox.\n"
          "\n"
          "No modification of *pnd_field* is made for the pressure dimension.\n"
          "At the latitude and longitude cloudbox edge points *pnd_field* is set to\n"
          "zero. This corresponds to nzero=1. If you want a larger margin between\n"
          "the lat and lon cloudbox edges and the \"cloud\" you increase\n"
          "*nzero*, where *nzero* is the number of grid points for which\n"
          "*pnd_field* shall be set to 0, counted from each lat and lon edge.\n"
          "\n"
          "See further *AtmFieldsExpand1D*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("pnd_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pnd_field", "atmosphere_dim", "cloudbox_on", "cloudbox_limits"),
      GIN("nzero"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Number of zero values inside lat and lon limits.")));

  md_data_raw.push_back(MdRecord(
      NAME("pnd_fieldZero"),
      DESCRIPTION(
          "Sets *pnd_field* to zero.\n"
          "\n"
          "Creates an empty *pnd_field* of cloudbox size according to\n"
          "*cloudbox_limits* and with number of scattering elemements\n"
          "according to *scat_data*. If *scat_data* is not set yet, it will be\n"
          "filled with one dummy scattering element.\n"
          "\n"
          "The method works with both *scat_data* and *scat_data_raw*."
          "\n"
          "This method primarily exists for testing purposes.\n"
          "On the one hand, empty *pnd_field* runs can be used to test the\n"
          "agreement between true clear-sky (*cloudboxOff*) solutions and the\n"
          "scattering solver solution in factual clear-sky conditions. It is\n"
          "important to avoid discontinuities when switching from thin-cloud\n"
          "to clear-sky conditions.\n"
          "Moreover, scattering calculations using the DOIT method include\n"
          "interpolation errors. If one is interested in this effect, one\n"
          "should compare the DOIT result with an empty cloudbox to a clearsky\n"
          "calculation. That means that the iterative method is performed for\n"
          "a cloudbox with no particles.\n"),
      AUTHORS("Claudia Emde, Jana Mendrok"),
      OUT("pnd_field", "dpnd_field_dx", "scat_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data",
         "atmosphere_dim",
         "f_grid",
         "cloudbox_limits",
         "jacobian_quantities"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("ppath_fieldFromDownUpLimbGeoms"),
      DESCRIPTION(
          "Computes ppath_field from \"standalone\" sensors looking upwards from\n"
          "0 m altitude with zenith angles range [0, 90], downwards from the top\n"
          "of the atmosphere covering the zenith angle range from 180 degrees to\n"
          "the surface tangent minus 1e-4 degrees, and through the limb covering\n"
          "at the same position as the downwards looking sensor covering the zenith\n"
          "angle range from the surface tangent plus 1e-4 degrees to 90 degrees minus\n"
          "1e-4 degrees.\n"
          "\n"
          "The top of the atmosphere is from *z_field*(-1, 0, 0) [python range notation].\n"
          "\n"
          "The field will consist of 3*nz arrays structured as [up, limb, down]\n"
          "\n"
          "The intent of this function is to generate a field so that calculations\n"
          "of *ppvar_iy* of all the fields will cover the zenith angle space\n"
          "of all positions in *z_field*.\n"
          "\n"
          "Only works for *atmosphere_dim* 1, spherical planets, and *ppath_lmax*<0\n"),
      AUTHORS("Richard Larsson"),
      OUT("ppath_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath_agenda",
         "ppath_lmax",
         "ppath_lraytrace",
         "atmgeom_checked",
         "z_field",
         "f_grid",
         "cloudbox_on",
         "cloudbox_checked",
         "ppath_inside_cloudbox_do",
         "rte_pos",
         "rte_los",
         "rte_pos2",
         "refellipsoid",
         "atmosphere_dim"),
      GIN("nz"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("3"),
      GIN_DESC("Number of zenith angles per position")));

  md_data_raw.push_back(MdRecord(
      NAME("ppathCalc"),
      DESCRIPTION(
          "Stand-alone calculation of propagation paths.\n"
          "\n"
          "Beside a few checks of input data, the only operation of this\n"
          "method is to execute *ppath_agenda*.\n"
          "\n"
          "Propagation paths are normally calculated as part of the radiative\n"
          "transfer calculations, and this method is not part of the control\n"
          "file. A reason to call this function directly would be to obtain a\n"
          "propagation path for plotting. Anyhow, use this method instead\n"
          "of calling e.g.*ppathStepByStep directly.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppath"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath_agenda",
         "ppath_lmax",
         "ppath_lraytrace",
         "atmgeom_checked",
         "f_grid",
         "cloudbox_on",
         "cloudbox_checked",
         "ppath_inside_cloudbox_do",
         "rte_pos",
         "rte_los",
         "rte_pos2"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("ppath_fieldCalc"),
      DESCRIPTION(
          "Stand-alone calculation of propagation path field from sensors.\n"
          "\n"
          "Uses *ppathCalc* internally.\n"),
      AUTHORS("Richard Larsson"),
      OUT("ppath_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath_agenda",
         "ppath_lmax",
         "ppath_lraytrace",
         "atmgeom_checked",
         "f_grid",
         "cloudbox_on",
         "cloudbox_checked",
         "ppath_inside_cloudbox_do",
         "sensor_pos",
         "sensor_los",
         "rte_pos2"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("ppathCalcFromAltitude"),
      DESCRIPTION(
          "Moves *rte_pos* forwards to near altitude before calling *ppathCalc*\n"
          "to compute a different *ppath*.  The accuracy-variable gives minimum\n"
          "distance before the input altitude.\n"
          "\n"
          "The forward-moving algorithm calls *ppathCalc* several\n"
          "times at reduced maximum distances.  The intention is to maintain\n"
          "the correct *rte_los* for a given *rte_pos* at all altitudes.  The\n"
          "method is thus relatively slow, and VERY memory intense at low\n"
          "accuracy.\n"
          "\n"
          "Intended to be used with \"tropospheric corrections\" from ground\n"
          "geometry.  Not well-tested\n"
          "\n"
          "Throws error if no altitude is in line of sight.\n"),
      AUTHORS("Richard Larsson"),
      OUT("ppath"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath_agenda",
         "ppath_lmax",
         "ppath_lraytrace",
         "atmgeom_checked",
         "f_grid",
         "cloudbox_on",
         "cloudbox_checked",
         "ppath_inside_cloudbox_do",
         "rte_pos",
         "rte_los",
         "rte_pos2"),
      GIN("altitude", "accuracy"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "0.5"),
      GIN_DESC("Altitude to move forward towards", "Accuracy of altitude")));

  md_data_raw.push_back(MdRecord(
      NAME("ppathFromRtePos2"),
      DESCRIPTION(
          "Determines the propagation path from *rte_pos2* to *rte_pos*.\n"
          "\n"
          "The propagation path linking *rte_pos* and *rte_pos2* is calculated\n"
          "and returned. The method determines the path in a pure numerical\n"
          "manner, where a simple algorithm is applied. The task is to find\n"
          "the value of *rte_los* (at *rte_pos*) linking the two positions.\n"
          "\n"
          "See the user guide for a description of the search algorithm,\n"
          "including a more detailed definition of *za_accuracy*, \n"
          "*pplrt_factor* and *pplrt_lowest*.\n"
          "\n"
          "The standard application of this method should be to radio link\n"
          "calculations, where *rte_pos2* corresponds to a transmitter, and\n"
          "*rte_pos* to the receiver/sensor.\n"
          "\n"
          "The details of the ray tracing is controlled by *ppath_step_agenda*\n"
          "as usual.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppath", "rte_los", "ppath_lraytrace"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath_step_agenda",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "z_field",
         "f_grid",
         "refellipsoid",
         "z_surface",
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

  md_data_raw.push_back(MdRecord(
      NAME("ppathPlaneParallel"),
      DESCRIPTION(
          "Propagation path calculations for a plane parallel atmosphere.\n"
          "\n"
          "This method basically assumes that the planet's radius is infinite,\n"
          "i.e. the planet surface has no curvature. Some consequences of this\n"
          "assumption:\n"
          "   - the mathod can only be used for 1D\n"
          "   - zenith angles between 89.9 and 90.1 deg are not allowed\n"
          "   - refraction is always neglected\n"
          "   - radii in ppath are set to Inf\n"
          "\n"
          "Notice that the method provides full propagation paths. This means\n"
          "that *ppath_step_agenda* is ignored (and thus also refraction).\n"
          "On the other hand, the method considers the cloudbox exactly as\n"
          "the standard path calculations.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppath"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "z_field",
         "z_surface",
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

  md_data_raw.push_back(MdRecord(
      NAME("ppathStepByStep"),
      DESCRIPTION(
          "Standard method for calculation of propagation paths.\n"
          "\n"
          "This method calculates complete propagation paths in a stepwise\n"
          "manner. Each step is denoted as a \"ppath_step\" and is the path\n"
          "through/inside a single grid box.\n"
          "\n"
          "The definition of a propgation path cannot be accommodated here.\n"
          "For more information read the chapter on propagation paths in the\n"
          "ARTS user guide.\n"
          "\n"
          "This method should never be called directly. Use *ppathCalc* instead\n"
          "if you want to extract propagation paths.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppath"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath_step_agenda",
         "ppath_inside_cloudbox_do",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "z_field",
         "f_grid",
         "refellipsoid",
         "z_surface",
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

  md_data_raw.push_back(MdRecord(
      NAME("ppathWriteXMLPartial"),
      DESCRIPTION(
          "WSM to only write a reduced Ppath, omitting grid positions.\n"
          "\n"
          "The following fields are set to be empty: gp_p, gp_lat and gp_lon.\n"
          "This cam drastically decrease the time for reading the structure\n"
          "by some external software.\n"
          "\n"
          "If *file_index is >= 0, the variable is written to a file with name:\n"
          "   <filename>.<file_index>.xml.\n"
          "where <file_index> is the value of *file_index*.\n"
          "\n"
          "This means that *filename* shall here not include the .xml\n"
          "extension. Omitting filename works as for *WriteXML*.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("output_file_format", "ppath"),
      GIN("filename", "file_index"),
      GIN_TYPE("String", "Index"),
      GIN_DEFAULT("", "-1"),
      GIN_DESC("File name. See above.",
               "Optional file index to append to filename.")));

  md_data_raw.push_back(MdRecord(
      NAME("ppath_stepGeometric"),
      DESCRIPTION(
          "Calculates a geometrical propagation path step.\n"
          "\n"
          "This function determines a propagation path step by pure\n"
          "geometrical calculations. That is, refraction is neglected. Path\n"
          "points are always included for crossings with the grids, tangent\n"
          "points and intersection points with the surface. The WSV *ppath_lmax*\n"
          "gives the option to include additional points to ensure that the\n"
          "distance along the path between the points does not exceed the\n"
          "selected maximum length. No additional points are included if\n"
          "*ppath_lmax* is set to <= 0.\n"
          "\n"
          "For further information, type see the on-line information for\n"
          "*ppath_step_agenda*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppath_step"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath_step",
         "atmosphere_dim",
         "lat_grid",
         "lon_grid",
         "z_field",
         "refellipsoid",
         "z_surface",
         "ppath_lmax"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("ppath_stepRefractionBasic"),
      DESCRIPTION(
          "Calculates a propagation path step, considering refraction by a\n"
          "basic approach.\n"
          "\n"
          "Refraction is taken into account by probably the simplest approach\n"
          "possible. The path is treated to consist of piece-wise geometric\n"
          "steps. A geometric path step is calculated from each point by\n"
          "using the local line-of-sight. Snell's law for spherical symmetry\n"
          "is used for 1D to determine the zenith angle at the new point.\n"
          "For 2D and 3D, the zenith angle is calculated using the average\n"
          "gradient of the refractive index between the two points. For 3D,\n"
          "the azimuth angle is treated in the same way as the zenith one.\n"
          "\n"
          "The maximum length of each ray tracing step is given by the WSV\n"
          "*ppath_lraytrace*. The length will never exceed the given maximum,\n"
          "but it can be smaller. The ray tracing steps are only used to\n"
          "determine the path. Points to describe the path are included as\n"
          "for *ppath_stepGeometric*, this including the functionality of\n"
          "*ppath_lmax*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("ppath_step"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("refr_index_air_agenda",
         "ppath_step",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "z_field",
         "t_field",
         "vmr_field",
         "refellipsoid",
         "z_surface",
         "f_grid",
         "ppath_lmax",
         "ppath_lraytrace"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("ppvar_optical_depthFromPpvar_trans_cumulat"),
      DESCRIPTION(
          "Sets *ppvar_optical_depth* according to provided transmission data.\n"
          "\n"
          "The values in ppvar_optical_depth are set to\n"
          "-log( ppvar_trans_cumulat(joker,joker,0,0) ).\n"),
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

  md_data_raw.push_back(
      MdRecord(NAME("Print"),
               DESCRIPTION("Prints a variable on the screen.\n"),
               AUTHORS("Oliver Lemke"),
               OUT(),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN("in", "level"),
               GIN_TYPE("Any", "Index"),
               GIN_DEFAULT(NODEF, "1"),
               GIN_DESC("Variable to be printed.", "Output level to use."),
               SETMETHOD(false),
               AGENDAMETHOD(false),
               USES_TEMPLATES(true)));

  md_data_raw.push_back(
      MdRecord(NAME("PrintPhysicalConstants"),
               DESCRIPTION("Prints (most) physical constants used in ARTS.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("PrintSelfLineMixingStatus"),
      DESCRIPTION("Test function for printing status of linemixing.\n"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_band", "abs_species_per_band"),
      GIN("temperature"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("296"),
      GIN_DESC("Temperature to evaluate at")));

  md_data_raw.push_back(
      MdRecord(NAME("PrintWorkspace"),
               DESCRIPTION("Prints a list of the workspace variables.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("ZFromPSimple"),
      DESCRIPTION(
          "Simple conversion from pressure to altitude.\n"
          "\n"
          "This function converts a vector of pressure values to an approximate vector\n"
          "of corresponding heights. The formula used to convert pressure to height is:\n"
          "z = 16000 * (5.0 - log10(p))"
          "That is, a pressure is  assumed to decrease by a factor of 10 every 16km.\n"
          "\n"
          "Note that all pressure values in the vector must be greater than 0.01.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("z_grid"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Approximate heights of pressure grid points."),
      IN(),
      GIN("p_grid"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Pressure grid."),
      SETMETHOD(false),
      AGENDAMETHOD(false)));

  md_data_raw.push_back(MdRecord(
      NAME("PFromZSimple"),
      DESCRIPTION(
          "Simple conversion from altitude to pressure.\n"
          "\n"
          "This function converts a vector of altitudes to an approximate vector\n"
          "of corresponding pressures. The formula used to convert altitide z to height\n"
          " is:\n"
          "p = 10.0^(5.0 - z / 1600)\n"
          "\n"
          "Note that all altitude values in the vector must be less than 120 km, \n"
          " otherwise an error will be thrown.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("p_grid"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Approximate pressures of corresponding to given altitudes."),
      IN(),
      GIN("z_grid"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Altitude grid."),
      SETMETHOD(false),
      AGENDAMETHOD(false)));

  md_data_raw.push_back(MdRecord(
      NAME("propmat_clearskyAddFaraday"),
      DESCRIPTION(
          "Calculates absorption matrix describing Faraday rotation.\n"
          "\n"
          "Faraday rotation is a change of polarization state of an\n"
          "electromagnetic wave propagating through charged matter by\n"
          "interaction with a magnetic field. Hence, this method requires\n"
          "*abs_species* to contain 'free_electrons' and electron content field\n"
          "(as part of *vmr_field*) as well as magnetic field (*mag_u_field*,\n"
          "*mag_v_field*, *mag_w_field*) to be specified.\n"
          "\n"
          "Faraday rotation affects Stokes parameters 2 and 3 (but not\n"
          "intensity!). Therefore, this method requires stokes_dim>2.\n"
          "\n"
          "Like all 'propmat_clearskyAdd*' methods, the method is additive,\n"
          "i.e., does not overwrite the propagation matrix *propmat_clearsky*,\n"
          "but adds further contributions.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("propmat_clearsky", "dpropmat_clearsky_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
         "dpropmat_clearsky_dx",
         "stokes_dim",
         "atmosphere_dim",
         "f_grid",
         "abs_species",
         "jacobian_quantities",
         "rtp_vmr",
         "rtp_los",
         "rtp_mag"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("propmat_clearskyAddFromAbsCoefPerSpecies"),
      DESCRIPTION(
          "Copy *propmat_clearsky* from *abs_coef_per_species*. This is handy for putting an\n"
          "explicit line-by-line calculation into the\n"
          "*propmat_clearsky_agenda*. This method is also used internally by.\n"
          "*propmat_clearskyAddOnTheFly*.\n"
          "Like all other propmat_clearsky methods, this method does not overwrite\n"
          "prior content of *propmat_clearsky*, but adds to it.\n"),
      AUTHORS("Stefan Buehler"),
      OUT("propmat_clearsky", "dpropmat_clearsky_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky", "abs_coef_per_species", "dabs_coef_dx"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("propmat_clearskyAddFromLookup"),
      DESCRIPTION(
          "Extract gas absorption coefficients from lookup table.\n"
          "\n"
          "This extracts the absorption coefficient for all species from the\n"
          "lookup table, and adds them to the propagation matrix. Extraction is\n"
          "for one specific atmospheric condition, i.e., a set of pressure,\n"
          "temperature, and VMR values.\n"
          "\n"
          "Some special species are ignored, for example Zeeman species and free\n"
          "electrons, since their absorption properties are not simple scalars\n"
          "and cannot be handled by the lookup table.\n"
          "\n"
          "The interpolation order in T and H2O is given by *abs_t_interp_order*\n"
          "and *abs_nls_interp_order*, respectively.\n"
          "\n"
          "Extraction is done for the frequencies in f_grid. Frequency\n"
          "interpolation is controlled by *abs_f_interp_order*. If this is zero,\n"
          "then f_grid must either be the same as the internal frequency grid of\n"
          "the lookup table (for efficiency reasons, only the first and last\n"
          "element of f_grid are checked), or must have only a single element.\n"
          "If *abs_f_interp_order* is above zero, then frequency is interpolated\n"
          "along with the other interpolation dimensions. This is useful for\n"
          "calculations with Doppler shift.\n"
          "\n"
          "For Doppler calculations, you should generate the table with a\n"
          "somewhat larger frequency grid than the calculation itself has, since\n"
          "the Doppler shift will push the frequency grid out of the table range\n"
          "on one side.\n"
          "\n"
          "Some extrapolation is allowed. For pressure and frequency interpolation\n"
          "the standard extrapolation factor of 0.5 is applied. The factor is the\n"
          "default for temperature and VMR interpolation, but the extrapolation\n"
          "limit can here be adjusted by the *extpolfac* argument.\n"
          "\n"
          "See also: *propmat_clearskyAddOnTheFly*.\n"),
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
         "rtp_pressure",
         "rtp_temperature",
         "rtp_vmr",
         "jacobian_quantities"),
      GIN("extpolfac"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("0.5"),
      GIN_DESC("Extrapolation factor (for temperature and VMR grid edges).")));

  md_data_raw.push_back(MdRecord(
      NAME("propmat_clearskyAddOnTheFly"),
      DESCRIPTION(
          "Calculates gas absorption coefficients line-by-line.\n"
          "\n"
          "This method can be used inside *propmat_clearsky_agenda* just like\n"
          "*propmat_clearskyAddFromLookup*. It is a shortcut for putting in some\n"
          "other methods explicitly, namely:\n"
          "\n"
          "  1. *AbsInputFromRteScalars*\n"
          "  2. Execute *abs_xsec_agenda*\n"
          "  3. *abs_coefCalcFromXsec*\n"
          "  4. *propmat_clearskyAddFromAbsCoefPerSpecies*\n"
          "\n"
          "The calculation is for one specific atmospheric condition, i.e., a set\n"
          "of pressure, temperature, and VMR values.\n"),
      AUTHORS("Stefan Buehler, Richard Larsson"),
      OUT("propmat_clearsky",
          "nlte_source",
          "dpropmat_clearsky_dx",
          "dnlte_dx_source",
          "nlte_dsource_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
         "nlte_source",
         "dpropmat_clearsky_dx",
         "dnlte_dx_source",
         "nlte_dsource_dx",
         "f_grid",
         "abs_species",
         "jacobian_quantities",
         "rtp_pressure",
         "rtp_temperature",
         "rtp_nlte",
         "rtp_vmr",
         "abs_xsec_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("propmat_clearskyAddParticles"),
      DESCRIPTION(
              "Calculates absorption coefficients of particles to be used in\n"
              "clearsky (non-cloudbox) calculations.\n"
              "\n"
              "This is a method to include particles (neglecting possible\n"
              "scattering components) in a clearsky calculation, i.e. without\n"
              "applying the cloudbox and scattering solvers. Particles are handled\n"
              "as absorbing species with one instance of 'particles' per scattering\n"
              "element considered added to *abs_species*. Particle absorption cross-\n"
              "sections at current atmospheric conditions are extracted from the\n"
              "single scattering data stored in *scat_data*, i.e., one array\n"
              "element per 'particles' instance in *abs_species* is required. Number\n"
              "densities are stored in *vmr_field_raw* or *vmr_field* as for all\n"
              "*abs_species*, but can be taken from (raw) pnd_field type data.\n"
              "\n"
              "Note that the absorption coefficient is applied both in the\n"
              "extinction term (neglecting scattering out of the line of sight)\n"
              "and the emission term (neglecting the scattering source term, i.e.\n"
              "scattering into the line of sight).\n"
              "\n"
              "Optionally, particle extinction (sum of absorption and scattering\n"
              "coefficient) can be used instead of absorption only. To choose this\n"
              "case, set the *use_abs_as_ext* flag to 0. However, be aware that\n"
              "this creates some unphysical emission term, hence is only suitable,\n"
              "where the source term is negligible anyways, e.g. for occultation\n"
              "simulations.\n"
              "\n"
              "A line-of-sight direction *rtp_los* is required as particles can\n"
              "exhibit directional dependent absorption properties, which is taken\n"
              "into account by this method."
              "\n"
              "*ScatElementsToabs_speciesAdd* can be used to add all required\n"
              "settings/data for individual scattering elements at once, i.e. a\n"
              " 'particles' tag to *abs_species*, a set of single scattering data to\n"
              "*scat_data* and a number density field to *vmr_field_raw*\n"
              "(*vmr_field* is derived applying AtmFieldsCalc once VMRs for all\n"
              "*abs_species* have been added) is appended for each scattering\n"
              "element.\n"
              "\n"
              "Like all 'propmat_clearskyAdd*' methods, the method is additive,\n"
              "i.e., does not overwrite the propagation matrix *propmat_clearsky*,\n"
              "but adds further contributions.\n"),
      AUTHORS("Jana Mendrok"),
      OUT("propmat_clearsky", "dpropmat_clearsky_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
         "dpropmat_clearsky_dx",
         "stokes_dim",
         "atmosphere_dim",
         "f_grid",
         "abs_species",
         "jacobian_quantities",
         "rtp_vmr",
         "rtp_los",
         "rtp_temperature",
         "scat_data",
         "scat_data_checked"),
      GIN("use_abs_as_ext"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("A flag with value 1 or 0. If set to one, particle absorption\n"
               "is used in extinction and emission parts of the RT equation,\n"
               "and scattering out of LOS as well as into LOS is neglected.\n"
               "Otherwise, particle extinction (absorption+scattering) is\n"
               "applied in both the extinction as well as the emission part\n"
               "of the RT equation. That is, true extinction is applied, but\n"
               "emission also includes a pseudo-emission contribution from\n"
               "the scattering coefficient.\n")));

  md_data_raw.push_back(MdRecord(
      NAME("propmat_clearskyAddZeeman"),
      DESCRIPTION(
          "Calculates Zeeman-affected polarized propagation matrix and its\n"
          "derivatives.\n"
          "\n"
          "Otherwise as *propmat_clearskyAddFromLookup* except line cutoff,\n"
          "shape, and, normalization factors have to be set in the LineRecord\n"),
      AUTHORS("Richard Larsson"),
      OUT("propmat_clearsky",
          "nlte_source",
          "dpropmat_clearsky_dx",
          "dnlte_dx_source",
          "nlte_dsource_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("propmat_clearsky",
         "nlte_source",
         "dpropmat_clearsky_dx",
         "dnlte_dx_source",
         "nlte_dsource_dx",
         "zeeman_linerecord_precalc",
         "f_grid",
         "abs_species",
         "jacobian_quantities",
         "isotopologue_ratios",
         "partition_functions",
         "rtp_pressure",
         "rtp_temperature",
         "rtp_nlte",
         "rtp_vmr",
         "rtp_mag",
         "rtp_los",
         "atmosphere_dim"),
      GIN("manual_zeeman_tag",
          "manual_zeeman_magnetic_field_strength",
          "manual_zeeman_theta",
          "manual_zeeman_eta"),
      GIN_TYPE("Index", "Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT("0", "1.0", "0.0", "0.0"),
      GIN_DESC("Manual angles tag",
               "Manual Magnetic Field Strength",
               "Manual theta given positive tag",
               "Manual eta given positive tag")));

  md_data_raw.push_back(MdRecord(
      NAME("propmat_clearskyInit"),
      DESCRIPTION(
          "Initialize *propmat_clearsky* and *nlte_source*.\n"
          "\n"
          "This method must be used inside *propmat_clearsky_agenda* and then\n"
          "be called first.\n"),
      AUTHORS("Oliver Lemke, Richard Larsson"),
      OUT("propmat_clearsky",
          "nlte_source",
          "dpropmat_clearsky_dx",
          "dnlte_dx_source",
          "nlte_dsource_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species",
         "jacobian_quantities",
         "f_grid",
         "stokes_dim",
         "propmat_clearsky_agenda_checked",
         "nlte_do"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("propmat_clearskyZero"),
      DESCRIPTION(
          "Sets *propmat_clearsky* to match zero attenuation.\n"
          "\n"
          "Use this method just if you know what you are doing!\n"
          "\n"
          "If you want to make a calculation with no clear-sky attenuation at\n"
          "all, fill *propmat_clearsky_agenda* with this method and required\n"
          "Ignore statements (don't include *propmat_clearskyInit*).\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("propmat_clearsky"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid", "stokes_dim"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("propmat_clearskyForceNegativeToZero"),
      DESCRIPTION("Sets *propmat_clearsky* to match zero attenuation\n"
                  "if negative value.  Useful for line mixing in some cases.\n"
                  "\n"
                  "Use this method just if you know what you are doing!\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("propmat_clearsky_agenda_checkedCalc"),
      DESCRIPTION(
          "Checks if the *propmat_clearsky_agenda* contains all necessary\n"
          "methods to calculate all the species in *abs_species*.\n"
          "\n"
          "This method should be called just before the *propmat_clearsky_agenda*\n"
          "is used, e.g. *DoitGetIncoming*, *ybatchCalc*, *yCalc*\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("propmat_clearsky_fieldCalc"),
      DESCRIPTION(
          "Calculate (vector) gas absorption coefficients for all points in the\n"
          "atmosphere.\n"
          "\n"
          "This is useful in two different contexts:\n"
          "\n"
          "1. For testing and plotting gas absorption. (For RT calculations, gas\n"
          "absorption is calculated or extracted locally, therefore there is no\n"
          "need to calculate a global field. But this method is handy for easy\n"
          "plotting of absorption vs. pressure, for example.)\n"
          "\n"
          "2. Inside the scattering region, monochromatic absorption is\n"
          "pre-calculated for the entire atmospheric field.\n"
          "\n"
          "The calculation itself is performed by the\n"
          "*propmat_clearsky_agenda*.\n"),
      AUTHORS("Stefan Buehler, Richard Larsson"),
      OUT("propmat_clearsky_field", "nlte_source_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmfields_checked",
         "f_grid",
         "stokes_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "t_field",
         "vmr_field",
         "nlte_field",
         "mag_u_field",
         "mag_v_field",
         "mag_w_field",
         "propmat_clearsky_agenda"),
      GIN("doppler", "los"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT("[]", "[]"),
      GIN_DESC("A vector of doppler shift values in Hz. Must either be "
               "empty or have same dimension as p_grid.",
               "Line of sight")));

  md_data_raw.push_back(MdRecord(
      NAME("psdA12"),
      DESCRIPTION(
          "Abel and Boutle [2012] (A12) particle size distribution for rain.\n"
          "\n"
          "Reference: Abel and Boutle, An improved representation of the \n"
          "raindrop size distribution for single-moment microphysics schemes,\n"
          "QJRMS, 2012.\n"
          "\n"
          "This is a 1-parameter PSD, i.e. *pnd_agenda_input* shall have one\n"
          "column and *pnd_agenda_input_names* shall contain a single string.\n"
          "The input data in *pnd_agenda_input* shall be rain mass content in\n"
          "unit of [kg/m3]. The naming used is *pnd_agenda_input_names* is free\n"
          "but the same name must be used in *particle_bulkprop_names* and\n"
          "*dpnd_data_dx_names*.\n"
          "\n"
          "Particles are assumed to be near-spherical, ie. *psd_size_grid* can\n"
          "either be in terms of volume (or mass) equivalent diameter or\n"
          "maximum diameter.\n"
          "\n"
          "Derivatives are obtained by perturbation of 0.1%, but not less than\n"
          "1e-9 kg/m3.\n"
          "\n"
          "The validity range of mass content is not limited. Negative mass\n"
          "contents wil produce negative psd values following a distribution\n"
          "given by abs(IWC), ie. abs(psd)=f(abs(IWC)).\n"
          "\n"
          "If temperature is outside [*t_min*,*t_max*] psd=0 and dpsd=0 if\n"
          "picky=0, or an error is thrown if picky=1.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("psdF07"),
      DESCRIPTION(
          "The Field et al. [2007] (F07) particle size distribution for snow and\n"
          "cloud ice.\n"
          "\n"
          "This is a 1-parameter PSD, i.e. *pnd_agenda_input* shall have one\n"
          "column and *pnd_agenda_input_names* shall contain a single string.\n"
          "The input data in *pnd_agenda_input* shall be ice hydrometeor mass\n"
          "content in unit of [kg/m3]. The naming used is *pnd_agenda_input_names*\n"
          "is free but the same name must be used in *particle_bulkprop_names* and\n"
          "*dpnd_data_dx_names*.\n"
          "\n"
          "*psd_size_grid* shall contain size in terms of maximum diameter.\n"
          "\n"
          "Derivatives are obtained by perturbation of 0.1%, but not less than\n"
          "1e-9 kg/m3.\n"
          "\n"
          "Both parametrization for tropics and midlatitudes are handled,\n"
          "governed by setting of *regime*, where \"TR\" selectes the tropical\n"
          "case, and \"ML\" the midlatitude one.\n"
          "\n"
          "The validity range of mass content is not limited. Negative mass\n"
          "contents will produce negative psd values following a distribution\n"
          "given by abs(IWC), ie. abs(psd)=f(abs(IWC)).\n"
          "\n"
          "If temperature is outside [*t_min*,*t_max*] psd=0 and dpsd=0 if\n"
          "picky=0, or an error is thrown if picky=1.\n"
          "\n"
          "For temperatures below *t_min_psd*, the size distribution is\n"
          "calculated for T = *t_min_psd*. Likewise, for temperatures above\n"
          "*t_max_psd*, the distribution is derived for T = *t_max_psd*.\n"
          "\n"
          "Defaults of *t_min_psd* and *t_max_psd* were set considering that\n"
          "the parametrization has been derived from measurements over\n"
          "temperatures of -60C to 0C."
          "\n"
          "Checks of the sanity of the mass-dimension relationship are performed\n"
          "Errors are thrown if:\n"
          "- Mass-dimension relation exponent *scat_species_b* is outside\n"
          "  [*b_min*, *b_max*].\n"),
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
          "Low *b* limit (only if picky).",
          "High *b* limit (only if picky).",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(MdRecord(
      NAME("psdMgd"),
      DESCRIPTION(
          "Modified gamma distribution PSD using n0, mu, la and ga as parameters.\n"
          "\n"
          "The modified gamma distribution is a 4-parameter (n0, mu, la and ga)\n"
          "distribution [Petty & Huang, JAS, 2011)]:\n"
          "   n(x) = n0 * x^mu * exp( -la*x^ga )\n"
          "where x is particle size or mass.\n"
          "\n"
          "The parameters can be given in two ways, either by *pnd_agenda_input* or\n"
          "as GIN arguments. The first option allows the parameter to vary, while\n"
          "in the second case the parameter gets a constant value. If a parameter is\n"
          "part of *pnd_agenda_input*, the corresponding GIN argument must be set\n"
          "to NaN (which is default). This means that the number of columns in\n"
          "*pnd_agenda_input* and the number of non-Inf choices for n0, mu, la and\n"
          "ga must add up to four.\n"
          "\n"
          "Data in *pnd_agenda_input* are linked to the MGD parameters in term of\n"
          "order, the naming in *pnd_agenda_input_names* is free. If all four\n"
          "parameteras are specified by *pnd_agenda_input*, the data in the first\n"
          "column are taken as n0, the second column as mu etc. If a parameter\n"
          "is given as a GIN argument, the columns are just shifted with one position.\n"
          "For example, if mu and ga are specified as GIN arguments, *pnd_agenda_input*\n"
          "shall have two columns, with n0-values in the first one and la-values in\n"
          "the second one.\n"
          "\n"
          "The GIN route is especially suitable for selecting special cases of MGD.\n"
          "For example, by setting mu=0 and ga=1, an exponential PSD is obtained:\n"
          "   n(x) = n0 * exp( -la*x )\n"
          "With mu=1 and ga=1, the gamma PSD is obtained:\n"
          "   n(x) = n0 * x^mu *exp( -la*x )\n"
          "There should be little overhead in using the method for exponential\n"
          "and gamma PSD, there is an internal switch to dedicated expressions for\n"
          "those PSDs.\n"
          "\n"
          "Derivatives can only be obtained for parameters that are specified by\n"
          "*pnd_agenda_input*. That is, parameters that are set to a value beside\n"
          "Inf are treated as fixed, and can not be subject to a rterieval.\n"
          "\n"
          "If temperature is outside [*t_min*,*t_max*] psd=0 and dpsd=0 if\n"
          "picky=0, or an error is thrown if picky=1.\n"
          "\n"
          "These requirements apply to the MGD parameters:\n"
          "  la > 0\n"
          "  ga > 0\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("psdMgdMass"),
      DESCRIPTION(
          "Modified gamma distribution PSD, with mass content as input.\n"
          "\n"
          "See *psdMgd* for a defintion of MGD parameters and how this PSD is\n"
          "handled in ARTS. Only deviations with respect to *psdMgd* are\n"
          "described here.\n"
          "\n"
          "This version of MGD PSD takes mass content as first input argument.\n"
          "This means that the first column of *pnd_agenda_input* shall hold\n"
          "mass content data.\n"
          "\n"
          "The mass content basically replaces one of the standard parameters\n"
          "(n0, mu, la and ga). This parameter is denoted as the dependent one.\n"
          "The dependent parameter is selected by setting the corresponding GIN\n"
          "to -999. So far only n0 and la are allowed to be dependent.\n"
          "\n"
          "Regarding remaining columns in *pnd_agenda_input* and constant\n"
          "parameter values (by GIN) follows the same principle as for *psdMgd*,\n"
          "except that mass is always on column (as mentioned) and that there is\n"
          "no position in *pnd_agenda_input* for the dependent parameter.\n"
          "\n"
          "These requirements apply to the MGD parameters:\n"
          "  mu + scat_species_b + 1 > 0\n"
          "  la > 0\n"
          "  ga > 0\n"
          "  If la is the dependent parameter, mass content must be > 0.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("psdMgdMassNtot"),
      DESCRIPTION(
          "Modified gamma distribution PSD, with mass content and total number\n"
          "density (Ntot) as inputs.\n"
          "\n"
          "This version of MGD PSD works as *psdMgdMass*, but takes mass content\n"
          "and total number density as first two arguments. This means that the first\n"
          "and second column of *pnd_agenda_input* shall hold mass content and Ntot,\n"
          "respectively. Accordingly, the number of dependent parameters is two.\n"
          "\n"
          "These requirements apply:\n"
          "  mu + 1 > 0\n"
          "  la > 0\n"
          "  ga > 0\n"
          "  Ntot must be > 0.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("psdMgdMassMeanParticleMass"),
      DESCRIPTION(
          "Modified gamma distribution PSD, with mass content and mean particle\n"
          "mass (Mmean) as inputs.\n"
          "\n"
          "This version of MGD PSD works as *psdMgdMass*, but takes mass content\n"
          "and mean particle mass as first two arguments. This means that the first\n"
          "and second column of *pnd_agenda_input* shall hold mass content and\n"
          "Mmean, respectively. Accordingly, the number of dependent parameters\n"
          "is two.\n"
          "\n"
          "\"Mean particle mass\" is here defined as the mass content divided with\n"
          "the total number density.\n"
          "\n"
          "These requirements apply to the MGD parameters:\n"
          "  mu + 1 > 0\n"
          "  la > 0\n"
          "  ga > 0\n"
          "  Mmean must be > 0.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("psdMgdMassXmean"),
      DESCRIPTION(
          "Modified gamma distribution PSD, with mass content and mean size\n"
          "(Xmean) as inputs.\n"
          "\n"
          "This version of MGD PSD works as *psdMgdMass*, but takes mass content\n"
          "and mass size as first two arguments. This means that the first and\n"
          "second column of *pnd_agenda_input* shall hold mass content and Xmean,\n"
          "respectively. Accordingly, the number of dependent parameters is two.\n"
          "\n"
          "\"Mean size\" is here defined as mass weighted size. Remembering that\n"
          "mass is a*x^b, this mean size can be expressed as M_b+1/M_b where M_b\n"
          "is b:th moment of the PSD (see e.g. Eq. 17 in Petty&Huang, JAS, 2011).\n"
          "\n"
          "These requirements apply to the MGD parameters:\n"
          "  mu + scat_species_b + 1 > 0\n"
          "  la > 0\n"
          "  ga > 0\n"
          "  Xmean must be > 0.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("psdMgdMassXmedian"),
      DESCRIPTION(
          "Modified gamma distribution PSD, with mass content and median size\n"
          "(Xmedian) as inputs.\n"
          "\n"
          "This version of MGD PSD works as *psdMgdMass*, but takes mass content\n"
          "and median size as first two arguments. This means that the first and\n"
          "second column of *pnd_agenda_input* shall hold mass content and Xmedian,\n"
          "respectively. Accordingly, the number of dependent parameters is two.\n"
          "\n"
          "These requirements apply to the MGD parameters:\n"
          "  mu + scat_species_b + 1 > 0\n"
          "  la > 0\n"
          "  ga > 0\n"
          "  Xmedian must be > 0.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("psdD14"),
      DESCRIPTION(
          "Normalized PSD as proposed in Delanoë et al. ((2014)),\n"
          "'Normalized particle size distribution for remote sensing\n"
          "application', J. Geophys. Res. Atmos., 119, 4204–422.\n"
          "\n"
          "The PSD has two independent parameters *N0star*, the intercept\n"
          " parameter, and *Dm*, the volume-weighted diameter.\n"
          "This implementation expects as input two out of the following\n"
          "three quantities: *iwc*, *n0Star*, *Dm*. In this case one of\n"
          "the input parameters *iwc*, *N0start*, *Dm* must be set to -999.\n*"
          "It is also possible to provide only *iwc*, in which case an a\n"
          "priori assumption will be used to deduce *N0star* from *temperature*.\n"
          "In this case both *N0star* and *Dm* must be set to -999.0.\n"
          "\n"
          "This PSD is not defined for vanishing concentrations of\n"
          "scatterers as it requires normalization by *Dm*. It is up\n"
          "to the user to ensure that the value of *Dm* is sufficiently\n"
          "large. An error is thrown if *Dm* is zero or below the value\n"
          "provided by *Dm_min*.\n"),
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
                  "917.6",
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
          "*alpha* parameter of the shape function",
          "*beta* paramter of the shape function",
          "Low temperature limit to calculate a psd.",
          "High temperature limit to calculate a psd.",
          "Lower threshold for *Dm* below which an error is thrown.",
          "Flag whether to be strict with parametrization value checks.")));

  md_data_raw.push_back(MdRecord(
      NAME("psdMH97"),
      DESCRIPTION(
          "McFarquahar and Heymsfield [1997] (MH97) particle size distribution\n"
          "for cloud ice.\n"
          "\n"
          "This is a 1-parameter PSD, i.e. *pnd_agenda_input* shall have one\n"
          "column and *pnd_agenda_input_names* shall contain a single string.\n"
          "The input data in *pnd_agenda_input* shall be ice hydrometeor mass\n"
          "content in unit of [kg/m3]. The naming used is *pnd_agenda_input_names*\n"
          "is free but the same name must be used in *particle_bulkprop_names* and\n"
          "*dpnd_data_dx_names*.\n"
          "\n"
          "*psd_size_grid* shall contain size in terms of volume equivalent diameter.\n"
          "\n"
          "Derivatives are obtained by perturbation of 0.1%, but not less than\n"
          "1e-9 kg/m3.\n"
          "\n"
          "The validity range of mass content is not limited. Negative mass\n"
          "contents wil produce negative psd values following a distribution\n"
          "given by abs(IWC), ie. abs(psd)=f(abs(IWC)).\n"
          "\n"
          "If temperature is outside [*t_min*,*t_max*] psd=0 and dpsd=0 if\n"
          "picky=0, or an error is thrown if picky=1.\n"
          "\n"
          "For temperatures below *t_min_psd*, the size distribution is\n"
          "calculated for T = *t_min_psd*. Likewise, for temperatures above\n"
          "*t_max_psd*, the distribution is derived for T = *t_max_psd*.\n"
          "\n"
          "Defaults of *t_min_psd* and *t_max_psd* were set considering that\n"
          "the parametrization has been derived from measurements over\n"
          "temperatures of -70C to -20C."
          "\n"
          "The noisy option can not be used together with calculation of\n"
          "derivatives (ie. when *dpnd_data_dx_names* is not empty).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("psdMono"),
      DESCRIPTION(
          "Mono-dispersive PSD, with number density given.\n"
          "\n"
          "This is a 1-parameter PSD, i.e. *pnd_agenda_input* shall have one\n"
          "column and *pnd_agenda_input_names* shall contain a single string.\n"
          "The input data in *pnd_agenda_input* shall be number densities, in\n"
          "unit of [#/m3]. The naming used is *pnd_agenda_input_names* is free\n"
          "but the same name must be used in *particle_bulkprop_names* and\n"
          "*dpnd_data_dx_names*.\n"
          "\n"
          "The method checks that the scattering species indicated (by\n"
          "*species_index*) has a single element, and just inserts the provided\n"
          "number density in *psd_data*.\n"
          "\n"
          "If temperature is outside [*t_min*,*t_max*] psd=0 and dpsd=0 if\n"
          "picky=0, or an error is thrown if picky=1.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("psdMonoMass"),
      DESCRIPTION(
          "Mono-dispersive PSD, with mass content given.\n"
          "\n"
          "This is a 1-parameter PSD, i.e. *pnd_agenda_input* shall have one\n"
          "column and *pnd_agenda_input_names* shall contain a single string.\n"
          "The input data in *pnd_agenda_input* shall be mass contents, in\n"
          "unit of [#/m3]. The naming used is *pnd_agenda_input_names* is free\n"
          "but the same name must be used in *particle_bulkprop_names* and\n"
          "*dpnd_data_dx_names*.\n"
          "\n"
          "The method checks that the scattering species indicated (by\n"
          "*species_index*) has a single element, and sets *psd_data* based\n"
          "on the mass contents given and the particle mass (derived from\n"
          "*scat_meta*).\n"
          "\n"
          "If temperature is outside [*t_min*,*t_max*] psd=0 and dpsd=0 if\n"
          "picky=0, or an error is thrown if picky=1.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("psdMY05"),
      DESCRIPTION(
          "Calculates *psd_data* and  *dpsd_data_dx* following Milbrandt and Yau (2005)\n"
          "two moment particle size distribution (MY05) for cloud water, cloud ice,\n"
          "rain, snow, graupel and hail, which is used in the GEM model.\n"
          "\n"
          "WSM for use in *pnd_agenda_array* for mapping *particle_bulkprop_field*\n"
          "to *pnd_field* using *pnd_fieldCalcFromParticleBulkProps*.\n"
          "Produces the particle size distribution values (dN/dD) and their\n"
          "derivates with respect to independent variables x by *dpnd_data_dx_names*\n"
          "over multiple particle sizes and atmospheric levels (or SWC/T\n"
          "combinations).\n"
          "\n"
          "*psd_size_grid* is considered to be in terms of maximum diameter.\n"
          "WC is considered to be in terms of mass content (or mass density),\n"
          "ie. units of [kg/m3]. N_tot in terms of number density, ie. units of [1/m3] ."
          "\n"
          "Derivatives with respect to WC and N_tot are obtained analytically.\n"
          "\n"
          "Six particle size distributions for the different hydrometeors are handled,\n"
          "governed by setting of *hydrometeor_type*, where \n"
          "    \"cloud_water\" selects cloud liquid water , \n"
          "    \"cloud_ice\" selects cloud ice, \n"
          "    \"snow\" selects snow, \n"
          "    \"rain\" selects rain, \n"
          "    \"graupel\" selects graupel, and \n"
          "    \"hail\" selects hail, \n"
          "\n"
          "Requirements:\n"
          "\n"
          "*pnd_agenda_input_names* must include :\n"
          "    [\"X-mass_density\", \"X-number_density\" ]. \"X\" is an arbitrary name\n"
          "The entries in  *dpnd_data_dx_names* (ie. the allowed\n"
          "independent variablea ) can be \"X-mass_density\" and\\or \n"
          "\"X-number_density\".\n"
          "\n"
          "The validity range of WC is not limited. Negative WC will produce\n"
          "negative psd values following a distribution given by abs(WC), ie.\n"
          "abs(psd)=f(abs(WC)).\n"
          "\n"
          "If temperature is outside [*t_min*,*t_max*] psd=0 and dpsd=0 if\n"
          "picky=0, or an error is thrown if picky=1.\n"

          ),
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

  md_data_raw.push_back(MdRecord(
      NAME("psdSB06"),
      DESCRIPTION(
          "Calculates *psd_data* and  *dpsd_data_dx* following Seifert and Beheng (2006)\n"
          "two moment particle size distribution (SB06) for cloud water, cloud ice,\n"
          "rain, snow, graupel and hail, which is used in the ICON model.\n"
          "\n"
          "WSM for use in *pnd_agenda_array* for mapping *particle_bulkprop_field*\n"
          "to *pnd_field* using *pnd_fieldCalcFromParticleBulkProps*.\n"
          "Produces the particle size distribution values (dN/dD) and their\n"
          "derivates with respect to independent variables x by *dpnd_data_dx_names*\n"
          "over multiple particle sizes and atmospheric levels (or SWC/T\n"
          "combinations).\n"
          "\n"
          "*psd_size_grid* is considered to be in terms of mass.\n"
          "WC is considered to be in terms of mass content (or mass density),\n"
          "ie. units of [kg/m3]. N_tot in terms of number density, ie. units of [1/m3] ."
          "\n"
          "Derivatives with respect to WC and N_tot are obtained analytically.\n"
          "\n"
          "Six particle size distributions for the different hydrometeors are handled,\n"
          "governed by setting of *hydrometeor_type*, where \n"
          "    \"cloud_water\" selects cloud liquid water , \n"
          "    \"cloud_ice\" selects cloud ice, \n"
          "    \"snow\" selects snow, \n"
          "    \"rain\" selects rain, \n"
          "    \"graupel\" selects graupel, and \n"
          "    \"hail\" selects hail, \n"
          "\n"
          "Requirements:\n"
          "\n"
          "*pnd_agenda_input_names* must include :\n"
          "    [\"X-mass_density\", \"X-number_density\" ]. \"X\" is an arbitrary name\n"
          "The entries in  *dpnd_data_dx_names* (ie. the allowed\n"
          "independent variablea ) can be \"X-mass_density\" and\\or \n"
          "\"X-number_density\".\n"
          "\n"
          "The validity range of WC is not limited. Negative WC will produce\n"
          "negative psd values following a distribution given by abs(WC), ie.\n"
          "abs(psd)=f(abs(WC)).\n"
          "\n"
          "If temperature is outside [*t_min*,*t_max*] psd=0 and dpsd=0 if\n"
          "picky=0, or an error is thrown if picky=1.\n"

          ),
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

  md_data_raw.push_back(MdRecord(
      NAME("psdW16"),
      DESCRIPTION(
          "Wang et al. [2016] (W16) particle size distribution for rain.\n"
          "\n"
          "Reference: Wang et al., Investigation of liquid cloud microphysical\n"
          "properties of deep convective systems: 1. Parameterization raindrop\n"
          "size distribution and its application ..., 2016.\n"
          "\n"
          "This is a 1-parameter PSD, i.e. *pnd_agenda_input* shall have one\n"
          "column and *pnd_agenda_input_names* shall contain a single string.\n"
          "The input data in *pnd_agenda_input* shall be rain mass content in\n"
          "unit of [kg/m3]. The naming used is *pnd_agenda_input_names* is free\n"
          "but the same name must be used in *particle_bulkprop_names* and\n"
          "*dpnd_data_dx_names*.\n"
          "\n"
          "Particles are assumed to be near-spherical, ie. *psd_size_grid* can\n"
          "either be in terms of volume (or mass) equivalent diameter or\n"
          "maximum diameter.\n"
          "\n"
          "Derivatives are obtained by perturbation of 0.1%, but not less than\n"
          "1e-9 kg/m3.\n"
          "\n"
          "The validity range of mass content is not limited. Negative mass\n"
          "contents wil produce negative psd values following a distribution\n"
          "given by abs(IWC), ie. abs(psd)=f(abs(IWC)).\n"
          "\n"
          "If temperature is outside [*t_min*,*t_max*] psd=0 and dpsd=0 if\n"
          "picky=0, or an error is thrown if picky=1.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("p_gridDensify"),
      DESCRIPTION(
          "A simple way to make *p_grid* more dense.\n"
          "\n"
          "The method includes new values in *p_grid*. For each intermediate\n"
          "pressure range, *nfill* points are added. That is, setting *nfill*\n"
          "to zero returns an unmodified copy of *p_grid_old*. The number of\n"
          "elements of the new *p_grid* is (n0-1)*(1+nfill)+1, where n0 is the\n"
          "length of *p_grid_old*.\n"
          "\n"
          "The new points are distributed equidistant in log(p).\n"
          "\n"
          "For safety, new grid and old grid Vectors are not allowed to be the\n"
          "same variable (both will be needed later on for regridding of the\n"
          "atmospheric fields), and atmospheric field related *checked WSV are\n"
          "reset to 0 (unchecked).\n"),
      AUTHORS("Patrick Eriksson, Jana Mendrok"),
      OUT("p_grid", "atmfields_checked", "atmgeom_checked", "cloudbox_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("p_grid_old", "nfill"),
      GIN_TYPE("Vector", "Index"),
      GIN_DEFAULT(NODEF, "-1"),
      GIN_DESC(/* p_grid_old */
               "A copy of the current (the old) p_grid. Not allowed to be "
               "the same variable as the output *p_grid*.",
               /* nfill */
               "Number of points to add between adjacent pressure points."
               "The default value (-1) results in an error.")));

  md_data_raw.push_back(MdRecord(
      NAME("p_gridFromZRaw"),
      DESCRIPTION(
          "Sets *p_grid* according to input atmosphere's raw z_field, derived\n"
          "e.g. from *AtmRawRead*.\n"
          "Attention: as default only pressure values for altitudes >= 0 are\n"
          "extracted. If negative altitudes shall also be selected, set no_neg=0.\n"),
      AUTHORS("Claudia Emde, Jana Mendrok"),
      OUT("p_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("z_field_raw"),
      GIN("no_negZ"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Exclude negative altitudes.")));

  md_data_raw.push_back(MdRecord(
      NAME("p_gridFromGasAbsLookup"),
      DESCRIPTION("Sets *p_grid* to the pressure grid of *abs_lookup*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("p_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lookup"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("p_gridRefine"),
      DESCRIPTION(
          "Provides refined pressure grid.\n"
          "\n"
          "Created new pressure grid has (log10) spacings below a given\n"
          "threshold.\n"
          "\n"
          "For safety, new grid and old grid Vectors are not allowed to be the\n"
          "same variable (both will be needed later on for regridding of the\n"
          "atmospheric fields), and atmospheric field related *checked WSV are\n"
          "reset to 0 (unchecked).\n"),
      AUTHORS("Stefan Buehler, Jana Mendrok"),
      OUT("p_grid", "atmfields_checked", "atmgeom_checked", "cloudbox_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("p_grid_old", "p_step"),
      GIN_TYPE("Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC(/* p_grid_old */
               "A copy of the current (the old) p_grid. Not allowed to be "
               "the same variable as the output *p_grid*.",
               /* p_step */
               "Maximum step in log10(p[Pa]). If the pressure grid is "
               "coarser than this, additional points are added until each "
               "log step is smaller than this.")));

  md_data_raw.push_back(MdRecord(
      NAME("RadiationFieldSpectralIntegrate"),
      DESCRIPTION(
          "Integrates fields like *spectral_irradiance_field* or *doit_i_field*\n"
          "over frequency.\n"
          "Important, the first dimension must be the frequency dimension!\n"
          "If a field  like *doit_i_field* is input, the stokes dimension\n"
          "is also removed.\n"),
      AUTHORS("Manfred Brath"),
      OUT(),
      GOUT("radiation_field"),
      GOUT_TYPE("Tensor4, Tensor5"),
      GOUT_DESC("TBD"),
      IN("f_grid"),
      GIN("spectral_radiation_field"),
      GIN_TYPE("Tensor5, Tensor7"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("TBD")));

  md_data_raw.push_back(MdRecord(
      NAME("line_irradianceCalcForSingleSpeciesNonOverlappingLines"),
      DESCRIPTION("Computes the line irradiance and line transmission\n"
                  "\n"
                  "Presently only works for 1D atmospheres\n"),
      AUTHORS("Richard Larsson"),
      OUT("line_irradiance", "line_transmission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species",
         "abs_lines_per_species",
         "nlte_field",
         "vmr_field",
         "t_field",
         "z_field",
         "p_grid",
         "atmosphere_dim",
         "surface_props_data",
         "iy_space_agenda",
         "iy_surface_agenda",
         "iy_cloudbox_agenda",
         "propmat_clearsky_agenda",
         "water_p_eq_agenda"),
      GIN("df", "nz", "nf"),
      GIN_TYPE("Numeric", "Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("relative frequency to line center",
               "number of zeniths",
               "number of frequencies per line")));

  md_data_raw.push_back(MdRecord(
      NAME("line_irradianceCalcForSingleSpeciesNonOverlappingLinesPseudo2D"),
      DESCRIPTION("Computes the line irradiance and line transmission\n"
                  "\n"
                  "Presently only works for 1D atmospheres\n"),
      AUTHORS("Richard Larsson"),
      OUT("line_irradiance", "line_transmission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species",
         "abs_lines_per_species",
         "nlte_field",
         "vmr_field",
         "t_field",
         "z_field",
         "p_grid",
         "refellipsoid",
         "surface_props_data",
         "iy_main_agenda",
         "ppath_agenda",
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

  md_data_raw.push_back(MdRecord(
      NAME("RationalAdd"),
      DESCRIPTION(
          "Adds a Rational and a value (out = in+value).\n"
          "\n"
          "The result can either be stored in the same or another Rational.\n"
          "(in and out can be the same varible, but not out and value)\n"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Rational"),
      GOUT_DESC("Output Rational."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Rational", "Rational"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Rational.", "Value to add.")));

  md_data_raw.push_back(MdRecord(
      NAME("RationalInvScale"),
      DESCRIPTION(
          "Inversely scales/divides a Rational with a value (out = in/value).\n"
          "\n"
          "The result can either be stored in the same or another Rational.\n"
          "(in and out can be the same varible, but not out and value)\n"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Rational"),
      GOUT_DESC("Output Rational."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Rational", "Rational"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Rational.", "Scaling Rational.")));

  md_data_raw.push_back(MdRecord(
      NAME("RationalScale"),
      DESCRIPTION(
          "Scales/multiplies a Rational with a value (out = in*value).\n"
          "\n"
          "The result can either be stored in the same or another Rational.\n"
          "(in and out can be the same varible, but not out and value)\n"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Rational"),
      GOUT_DESC("Output Rational."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Rational", "Rational"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input Rational.", "Scaling value.")));

  md_data_raw.push_back(MdRecord(
      NAME("RationalSet"),
      DESCRIPTION("Sets a Rational workspace variable to the given value.\n"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Rational"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("numerator", "denominator"),
      GIN_TYPE("Index", "Index"),
      GIN_DEFAULT(NODEF, "1"),
      GIN_DESC("The numerator.", "The denominator.")));

  md_data_raw.push_back(MdRecord(
      NAME("ReadNetCDF"),
      DESCRIPTION("Reads a workspace variable from a NetCDF file.\n"
                  "\n"
                  "This method can read variables of any group.\n"
                  "\n"
                  "If the filename is omitted, the variable is read\n"
                  "from <basename>.<variable_name>.nc.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector, Matrix, Tensor3, Tensor4, Tensor5, ArrayOfVector,"
                "ArrayOfMatrix, GasAbsLookup"),
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

  md_data_raw.push_back(MdRecord(
      NAME("ReadXML"),
      DESCRIPTION(
          "Reads a workspace variable from an XML file.\n"
          "\n"
          "This method can read variables of any group.\n"
          "\n"
          "If the filename is omitted, the variable is read\n"
          "from <basename>.<variable_name>.xml.\n"
          "If the given filename does not exist, this method will\n"
          "also look for files with an added .xml, .xml.gz and .gz extension\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
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

  md_data_raw.push_back(MdRecord(
      NAME("ReadXMLIndexed"),
      DESCRIPTION("As *ReadXML*, but reads indexed file names.\n"
                  "\n"
                  "The variable is read from a file with name:\n"
                  "   <filename>.<file_index>.xml.\n"
                  "where <file_index> is the value of *file_index*.\n"
                  "\n"
                  "This means that *filename* shall here not include the .xml\n"
                  "extension. Omitting filename works as for *ReadXML*.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Any"),
      GOUT_DESC("Workspace variable to be read."),
      IN("file_index"),
      GIN("filename", "digits"),
      GIN_TYPE("String", "Index"),
      GIN_DEFAULT("", "0"),
      GIN_DESC(
          "File name. See above.",
          "Equalize the widths of all numbers by padding with zeros as necessary.\n"
          "0 means no padding (default)."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(MdRecord(
      NAME("Reduce"),
      DESCRIPTION(
          "Reduces a larger class to a smaller class of same size.\n"
          "\n"
          "The Reduce command reduces all \"1\"-dimensions to nil.  Examples:\n"
          "\t1) 1 Vector can be reduced to a Numeric\n"
          "\t2) 2x1 Matrix can be reduced to 2 Vector\n"
          "\t3) 1x3x1 Tensor3 can be reduced to 3 Vector\n"
          "\t4) 1x1x1x1 Tensor4 can be reduced to a Numeric\n"
          "\t5) 3x1x4x1x5 Tensor5 can only be reduced to 3x4x5 Tensor3\n"
          "\t6) 1x1x1x1x2x3 Tensor6 can be reduced to 2x3 Matrix\n"
          "\t7) 2x3x4x5x6x7x1 Tensor7 can be reduced to 2x3x4x5x6x7 Tensor6\n"
          "And so on\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("refellipsoidEarth"),
      DESCRIPTION(
          "Earth reference ellipsoids.\n"
          "\n"
          "The reference ellipsoid (*refellipsoid*) is set to model the Earth,\n"
          "following different models. The options are:\n"
          "\n"
          "   \"Sphere\" : A spherical Earth. The radius is set following\n"
          "      the value set for the Earth radius in constants.cc.\n"
          "\n"
          "   \"WGS84\" : The reference ellipsoid used by the GPS system.\n"
          "      Should be the standard choice for a non-spherical Earth.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("refellipsoid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("model"),
      GIN_TYPE("String"),
      GIN_DEFAULT("Sphere"),
      GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(
      MdRecord(NAME("refellipsoidGanymede"),
               DESCRIPTION("Ganymede reference ellipsoids.\n"
                           "\n"
                           "From Wikipedia\n"),
               AUTHORS("Takayoshi Yamada"),
               OUT("refellipsoid"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN(),
               GIN("model"),
               GIN_TYPE("String"),
               GIN_DEFAULT("Sphere"),
               GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(MdRecord(
      NAME("refellipsoidForAzimuth"),
      DESCRIPTION(
          "Conversion of 3D ellipsoid to 1D curvature radius.\n"
          "\n"
          "Calculates the curvature radius for the given latitude and azimuth\n"
          "angle, and uses this to set a spherical reference ellipsoid\n"
          "suitable for 1D calculations. The curvature radius is a better\n"
          "local approximation than using the local ellipsoid radius.\n"
          "\n"
          "The used expression assumes a geodetic latitude, but also\n"
          "latitudes should be OK as using this method anyhow signifies\n"
          "an approximation.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("refellipsoid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("refellipsoid"),
      GIN("latitude", "azimuth"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Latitude.", "Azimuth angle.")));

  md_data_raw.push_back(MdRecord(
      NAME("refellipsoidEuropa"),
      DESCRIPTION(
          "Io reference ellipsoids.\n"
          "\n"
          "The reference ellipsoid (*refellipsoid*) is set to model Io,\n"
          "folowing different models. The options are:\n"
          "\n"
          "   \"Sphere\" : A spherical planetesimal. The radius is taken from\n"
          "      report of the IAU/IAG Working Group.\n"),
      AUTHORS("Richard Larsson"),
      OUT("refellipsoid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("model"),
      GIN_TYPE("String"),
      GIN_DEFAULT("Sphere"),
      GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(MdRecord(
      NAME("refellipsoidIo"),
      DESCRIPTION(
          "Io reference ellipsoids.\n"
          "\n"
          "The reference ellipsoid (*refellipsoid*) is set to model Io,\n"
          "folowing different models. The options are:\n"
          "\n"
          "   \"Sphere\" : A spherical planetesimal. The radius is taken from\n"
          "      report of the IAU/IAG Working Group.\n"),
      AUTHORS("Richard Larsson"),
      OUT("refellipsoid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("model"),
      GIN_TYPE("String"),
      GIN_DEFAULT("Sphere"),
      GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(MdRecord(
      NAME("refellipsoidJupiter"),
      DESCRIPTION(
          "Jupiter reference ellipsoids.\n"
          "\n"
          "The reference ellipsoid (*refellipsoid*) is set to model Jupiter,\n"
          "folowing different models. The options are:\n"
          "\n"
          "   \"Sphere\" : A spherical planet. The radius is taken from a\n"
          "      report of the IAU/IAG Working Group.\n"
          "\n"
          "   \"Ellipsoid\" : A reference ellipsoid with parameters taken from\n"
          "      a report of the IAU/IAG Working Group.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("refellipsoid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("model"),
      GIN_TYPE("String"),
      GIN_DEFAULT("Sphere"),
      GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(MdRecord(
      NAME("refellipsoidMars"),
      DESCRIPTION(
          "Mars reference ellipsoids.\n"
          "\n"
          "The reference ellipsoid (*refellipsoid*) is set to model Mars,\n"
          "folowing different models. The options are:\n"
          "\n"
          "   \"Sphere\" : A spherical planet. The radius is taken from a\n"
          "      report of the IAU/IAG Working Group.\n"
          "\n"
          "   \"Ellipsoid\" : A reference ellipsoid with parameters taken from\n"
          "      a report of the IAU/IAG Working Group.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("refellipsoid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("model"),
      GIN_TYPE("String"),
      GIN_DEFAULT("Sphere"),
      GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(MdRecord(
      NAME("refellipsoidMoon"),
      DESCRIPTION(
          "Moon reference ellipsoids.\n"
          "\n"
          "The reference ellipsoid (*refellipsoid*) is set to model Moon,\n"
          "folowing different models. The options are:\n"
          "\n"
          "   \"Sphere\" : A spherical planet. The radius is taken from a\n"
          "      report of the IAU/IAG Working Group.\n"
          "\n"
          "   \"Ellipsoid\" : A reference ellipsoid with parameters taken from\n"
          "      Wikepedia (see code for details). The IAU/IAG working group\n"
          "      defines the Moon ellipsoid to be a sphere.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("refellipsoid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("model"),
      GIN_TYPE("String"),
      GIN_DEFAULT("Sphere"),
      GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(MdRecord(
      NAME("refellipsoidOrbitPlane"),
      DESCRIPTION(
          "Conversion of 3D ellipsoid to 2D orbit track geometry.\n"
          "\n"
          "Determines an approximate reference ellipsoid following an orbit\n"
          "track. The new ellipsoid is determined simply, by determining the\n"
          "radius at the maximum latitude and from this value calculate a new\n"
          "new eccentricity. The orbit is specified by giving the orbit\n"
          "inclination (*orbitinc*), that is normally a value around 100 deg\n"
          "for polar sun-synchronous orbits.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("refellipsoid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("refellipsoid"),
      GIN("orbitinc"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Orbit inclination.")));

  md_data_raw.push_back(MdRecord(
      NAME("refellipsoidSet"),
      DESCRIPTION(
          "Manual setting of the reference ellipsoid.\n"
          "\n"
          "The two values of *refellipsoid* can here be set manually. The two\n"
          "arguments correspond directly to first and second element of\n"
          "*refellipsoid*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("refellipsoid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("re", "e"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "0"),
      GIN_DESC("Average or equatorial radius.", "Eccentricity")));

  md_data_raw.push_back(MdRecord(
      NAME("refellipsoidVenus"),
      DESCRIPTION(
          "Venus reference ellipsoids.\n"
          "\n"
          "The reference ellipsoid (*refellipsoid*) is set to model Venus,\n"
          "folowing different models. The options are:\n"
          "\n"
          "   \"Sphere\" : A spherical planet. The radius is taken from a\n"
          "      report of the IAU/IAG Working Group.\n"
          "\n"
          "According to the report used above, the Venus ellipsoid lacks\n"
          "eccentricity and no further models should be required.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("refellipsoid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("model"),
      GIN_TYPE("String"),
      GIN_DEFAULT("Sphere"),
      GIN_DESC("Model ellipsoid to use. Options listed above.")));

  md_data_raw.push_back(MdRecord(
      NAME("refr_index_airFreeElectrons"),
      DESCRIPTION(
          "Microwave refractive index due to free electrons.\n"
          "\n"
          "The refractive index of free electrons is added to *refr_index_air*.\n"
          "To obtain the complete value, *refr_index_air* should be set to 1\n"
          "before calling this WSM. This applies also to *refr_index_air_group*.\n"
          "\n"
          "The expression applied is n=sqrt(1-wp^2/w^2) where wp is the plasma\n"
          "frequency, and w is the angular frequency (the function returns\n"
          "n-1, that here is slightly negative). This expressions is found in\n"
          "many textbooks, e.g. Rybicki and Lightman (1979). The above refers\n"
          "to *refr_index*. *refr_index_group* is sqrt(1+wp^2/w^2).\n"
          "\n"
          "The expression is dispersive. The frequency applied is the mean of\n"
          "first and last element of *f_grid* is selected. This frequency must\n"
          "be at least twice the plasma frequency.\n"
          "\n"
          "An error is issued if free electrons not are part of *abs_species*\n"
          "(and there exist a corresponding \"vmr\"-value). This demand is\n"
          "removed if *demand_vmr_value* is set to 0, but use this option\n"
          "with care.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("refr_index_airInfraredEarth"),
      DESCRIPTION(
          "Calculates the IR refractive index due to gases in the\n"
          "Earth's atmosphere.\n"
          "\n"
          "Only refractivity of dry air is considered. The formula used is\n"
          "contributed by Michael Hoepfner, Forschungszentrum Karlsruhe.\n"
          "\n"
          "The refractivity of dry air is added to *refr_index_air*. To obtain\n"
          "the complete value, *refr_index_air* should be set to 1 before\n"
          "calling this WSM. This applies also to *refr_index_air_group*.\n"
          "\n"
          "The expression used is non-dispersive. Hence, *refr_index_air* and\n"
          "*refr_index_air_group* are identical.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("refr_index_airMicrowavesEarth"),
      DESCRIPTION(
          "Microwave refractive index in Earth's atmosphere.\n"
          "\n"
          "This method just considers pressure, temperature and water\n"
          "vapour, which should suffice for Earth. For a more general\n"
          "method, see *refr_index_airMicrowavesGeneral*.\n"
          "\n"
          "The refractivity of dry air and water vapour is added to\n"
          "*refr_index_air*. To obtain the complete value, *refr_index_air*\n"
          "should be set to 1 before calling this WSM. This applies also to\n"
          "*refr_index_air_group.\n"
          "\n"
          "The expression used is non-dispersive. Hence, *refr_index_air*\n"
          "and *refr_index_air_group* are identical.\n"
          "\n"
          "The standard expression for Earth and microwaves is used:\n"
          "   N = k1*(P-e)/T + k2*e/T + k3*e/T^2\n"
          "where N is refractivity, P is pressure, T is temperature and\n"
          "e is water vapour partial pressure. The values of k1, k2 and k3\n"
          "can be modified.\n"
          "\n"
          "Many different values of k1, k2 and k3 can be found in the\n"
          "literature. The default values applied here are taken from\n"
          "Bevis et al., GPS meteorology: Mapping ..., JAM, 1994.\n"
          "More specifically, these value are found in Table 1, listed\n"
          "as \"Present study\". Note that in ARTS Pa is used for pressure\n"
          "and k1, k2 and k3 must be adjusted accordingly.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("refr_index_airMicrowavesGeneral"),
      DESCRIPTION(
          "Microwave refractive index due to gases in planetary atmospheres.\n"
          "\n"
          "The refractivity of a specified gas mixture is calculated and added\n"
          "to *refr_index_air*. To obtain the complete value, *refr_index_air*\n"
          "should be set to 1 before calling this WSM. This applies also to\n"
          "*refr_index_air_group.\n"
          "\n"
          "The expression used is non-dispersive. Hence, *refr_index_air* and\n"
          "*refr_index_air_group* are identical.\n"
          "\n"
          "Uses the methodology introduced by Newell&Baird (1965) for calculating\n"
          "refractivity of variable gas mixtures based on refractivity of the\n"
          "individual gases at reference conditions. Assuming ideal gas law for\n"
          "converting reference refractivity to actual pressure and temperature\n"
          "conditions. Reference refractivities are also taken from Newell&Baird (1965)\n"
          "and are vailable for N2, O2, CO2, H2, and He. Additionally, H2O reference\n"
          "refractivity has been derived from H2O contribution in Thayer (see\n"
          "*refr_index_airThayer*) for T0=273.15K. Any mixture of these gases\n"
          "can be taken into account.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("retrievalDefClose"),
      DESCRIPTION(
          "Closes the definition of retrieval quantities and correlations and\n"
          "prepares related WSVs for the retrieval.\n"
          "\n"
          "This function calls jacobianClose and checks that the corvariance matrices\n"
          "are consistent with the Jacobian.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("retrievalAddAbsSpecies"),
      DESCRIPTION(
          "Adds an absorption species to the retrieval quantities.\n"
          "\n"
          "Similar to *jacobianAddAbsSpecies* but also sets the corresponding block in\n"
          "*covmat_sx* to the matrices provided in *covmat_block* and *covmat_inv_block*.\n"
          "The dimensions of *covmat_block* are required to agree with the dimensions of the\n"
          "retrieval grid.\n"
          "\n"
          "*covmat_inv_block* must be either empty or the same dimension as *covmat_block*.\n"
          "If provided, this matrix will be used as the inverse for the covariance matrix block\n"
          "and numerical inversion of this block is thus avoided. Note, however, that this is\n"
          "only effective if this block is uncorrelated with any other retrieval quantity.\n"
          "\n"
          "For number and order of elements added to *x*, see *jacobianAddAbsSpecies*.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "covmat_block",
         "covmat_inv_block",
         "p_grid",
         "lat_grid",
         "lon_grid"),
      GIN("g1",
          "g2",
          "g3",
          "species",
          "unit",
          "for_species_tag"),
      GIN_TYPE("Vector",
               "Vector",
               "Vector",
               "String",
               "String",
               "Index"),
      GIN_DEFAULT(
          NODEF, NODEF, NODEF, NODEF, "rel", "1"),
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

  md_data_raw.push_back(MdRecord(
      NAME("retrievalAddFreqShift"),
      DESCRIPTION(
          "Same as *jacobianAddFreqShift* but also adds the correlation block\n"
          "contained in *covmat_block* and *covmat_inv_block* to *covmat_sx*.\n"
          "\n"
          "For number and order of elements added to *x*, see *jacobianAddFreqShift*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("retrievalAddFreqStretch"),
      DESCRIPTION(
          "Same as *jacobianAddFreqShift* but also adds the correlation block\n"
          "contained in *covmat_block* and *covmat_inv_block* to *covmat_sx*.\n"
          "\n"
          "For number and order of elements added to *x*, see *jacobianAddFreqStretch*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("retrievalDefInit"),
      DESCRIPTION(
          "Begin retrieval definition section.\n"
          "\n"
          "This function initialises all variables required for defining\n"
          "retrieval quantities and corresponding covariance matrices.\n"
          "By default, Jacobian quantities should be added withing the.\n"
          "retrieval definition section. If Jacobian quantities are\n"
          "defined separately *initialize_jacobian* must be set to 0,\n"
          "otherwise the quantities will be discarded.\n"),
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
      GIN_DESC("Flag whether or not to (re)initialize Jacobian-related\n"
               "quantities. Set to 0 if Jacobian is already defined.")));

  md_data_raw.push_back(MdRecord(
      NAME("retrievalAddCatalogParameter"),
      DESCRIPTION(
          "Siminlar to *jacobianAddCatalogParameter* but also adds a corresponding\n"
          "block to *covmat_sx* with the given *var* as variance value.\n"
          "\n"
          "For number and order of elements added to *x*, see *jacobianAddCatalogParameter*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("retrievalAddCatalogParameters"),
      DESCRIPTION(
          "Same as *jacobianAddCatalogParameters* but also adds a new block to *covmat_sx*\n"
          "using the matrices in *covmat_block* and *covmat_inv_block*.\n"
          "\n"
          "If *covmat_inv_block* is non-empty, it is used as inverse for the added block\n"
          "which avoids its numerical computation.\n"
          "\n"
          "For number and order of elements added to *x*, see *jacobianAddCatalogParameters*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("retrievalAddMagField"),
      DESCRIPTION(
          "Same as *jacobianAddMagField* but also adds a new block to *covmat_sx*\n"
          "using the matrices in *covmat_block* and *covmat_inv_block*.\n"
          "\n"
          "If *covmat_inv_block* is non-empty, it is used as inverse for the added block\n"
          "which avoids its numerical computation.\n"
          "\n"
          "For number and order of elements added to *x*, see *jacobianAddAbsMagField*.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "covmat_block",
         "covmat_inv_block",
         "p_grid",
         "lat_grid",
         "lon_grid"),
      GIN("g1", "g2", "g3", "component", "dB"),
      GIN_TYPE("Vector", "Vector", "Vector", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "v", "1.0e-7"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Magnetic field component to retrieve",
               "Magnetic field perturbation")));

  md_data_raw.push_back(MdRecord(
      NAME("retrievalAddPointingZa"),
      DESCRIPTION(
          "Same as *jacobianAddPointingZa* but also adds a new block to *covmat_sx*\n"
          "using the matrices in *covmat_block* and *covmat_inv_block*.\n"
          "\n"
          "If *covmat_inv_block* is non-empty, it is used as inverse for the added block\n"
          "which avoids its numerical computation.\n"
          "\n"
          "For number and order of elements added to *x*, see *jacobianAddPointingZa*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("retrievalAddPolyfit"),
      DESCRIPTION(
          "Same as *jacobianAddPolyfit* but also adds a new block to *covmat_sx*\n"
          "using the matrices in *covmat_block* and *covmat_inv_block*.\n"
          "\n"
          "If *covmat_inv_block* is non-empty, it is used as inverse for the added block\n"
          "which avoids its numerical computation.\n"
          "\n"
          "For number and order of elements added to *x*, see *jacobianAddPolyfit*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("retrievalAddScatSpecies"),
      DESCRIPTION(
          "Same as *jacobianAddPolyfit* but also adds a new block to *covmat_sx*\n"
          "using the matrices in *covmat_block* and *covmat_inv_block*.\n"
          "\n"
          "If *covmat_inv_block* is non-empty, it is used as inverse for the added block\n"
          "which avoids its numerical computation.\n"
          "\n"
          "For number and order of elements added to *x*, see *jacobianAddScatSpecies*.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "covmat_block",
         "covmat_inv_block",
         "p_grid",
         "lat_grid",
         "lon_grid"),
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

  md_data_raw.push_back(MdRecord(
      NAME("retrievalAddSinefit"),
      DESCRIPTION(
          "Same as *jacobianAddSinefit* but also adds a new block to *covmat_sx*\n"
          "using the matrices in *covmat_block* and *covmat_inv_block*.\n"
          "\n"
          "If *covmat_inv_block* is non-empty, it is used as inverse for the added block\n"
          "which avoids its numerical computation.\n"
          "\n"
          "For number and order of elements added to *x*, see *jacobianAddSinefit*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("retrievalAddSpecialSpecies"),
      DESCRIPTION(
          "Same as *jacobianAddSpecialSpecies* but also adds a new block to *covmat_sx*\n"
          "using the matrices in *covmat_block* and *covmat_inv_block*.\n"
          "\n"
          "If *covmat_inv_block* is non-empty, it is used as inverse for the added block\n"
          "which avoids its numerical computation.\n"
          "\n"
          "For number and order of elements added to *x*, see *jacobianAddSpecialSpecies*.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "covmat_block",
         "covmat_inv_block",
         "p_grid",
         "lat_grid",
         "lon_grid"),
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

  md_data_raw.push_back(MdRecord(
      NAME("retrievalAddSurfaceQuantity"),
      DESCRIPTION(
          "Same as *jacobianAddSurfaceQuantity* but also adds a new block to *covmat_sx*\n"
          "using the matrices in *covmat_block* and *covmat_inv_block*.\n"
          "\n"
          "If *covmat_inv_block* is non-empty, it is used as inverse for the added block\n"
          "which avoids its numerical computation.\n"
          "\n"
          "For number and order of elements added to *x*, see *jacobianAddSurfaceQuantity*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "covmat_block",
         "covmat_inv_block",
         "atmosphere_dim",
         "lat_grid",
         "lon_grid"),
      GIN("g1", "g2", "quantity"),
      GIN_TYPE("Vector", "Vector", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Retrieval quantity, e.g. \"Wind speed\".")));

  md_data_raw.push_back(MdRecord(
      NAME("retrievalAddTemperature"),
      DESCRIPTION(
          "Same as *jacobianAddTemperature* but also adds a new block to *covmat_sx*\n"
          "using the matrices in *covmat_block* and *covmat_inv_block*.\n"
          "\n"
          "If *covmat_inv_block* is non-empty, it is used as inverse for the added block\n"
          "which avoids its numerical computation.\n"
          "\n"
          "For number and order of elements added to *x*, see *jacobianAddTemperature*.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "covmat_block",
         "covmat_inv_block",
         "p_grid",
         "lat_grid",
         "lon_grid"),
      GIN("g1", "g2", "g3", "hse"),
      GIN_TYPE("Vector", "Vector", "Vector", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "on"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retreival grid.",
               "Flag to assume HSE or not (\"on\" or \"off\").")));

  md_data_raw.push_back(MdRecord(
      NAME("retrievalAddWind"),
      DESCRIPTION(
          "Same as *jacobianAddWind* but also adds a new block to *covmat_sx*\n"
          "using the matrices in *covmat_block* and *covmat_inv_block*.\n"
          "\n"
          "If *covmat_inv_block* is non-empty, it is used as inverse for the added block\n"
          "which avoids its numerical computation.\n"
          "\n"
          "For number and order of elements added to *x*, see *jacobianAddWind*.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("covmat_sx", "jacobian_quantities", "jacobian_agenda"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("covmat_sx",
         "jacobian_quantities",
         "jacobian_agenda",
         "atmosphere_dim",
         "covmat_block",
         "covmat_inv_block",
         "p_grid",
         "lat_grid",
         "lon_grid"),
      GIN("g1", "g2", "g3", "component", "dfrequency"),
      GIN_TYPE("Vector", "Vector", "Vector", "String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF, "v", "0.1"),
      GIN_DESC("Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retrieval grid.",
               "Wind component to retrieve",
               "This is the frequency perturbation")));

  md_data_raw.push_back(MdRecord(
      NAME("retrievalErrorsExtract"),
      DESCRIPTION(
          "Extract retrieval error from covariance matrices.\n"
          "\n"
          "Extracts the error estimates for the retrieved quantities from the covariance\n"
          "matrices for the error due to measurement noise *covmat_so* and the error due\n"
          "to limited resolution of the observation system *covmat_ss* and stores them in\n"
          "the vectors *retrieval_eo* and *retrieval_ss*, respectively."
          "\n"
          "To etract these errors, first the convariance matrices of which the errors \n"
          "should be extracted have to be computed using the WSMs *covmat_soCalc*\n"
          "and *covmat_ssCalc* or set to be empty in order to be ignored. Note, however,\n"
          "that this will also set the corresponding error vector to be empty.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("RT4Calc"),
      DESCRIPTION(
          "Interface to the PolRadTran RT4 scattering solver (Evans).\n"
          "\n"
          "DISCLAIMER: There is a number of known issues with the current\n"
          "implementation (see below). Use this WSM with care and only if\n"
          "these limitations/requirements are fulfilled. Results might be\n"
          "erroneous otherwise.\n"
          "\n"
          "RT4 provides the radiation field (*doit_i_field*) from a vector\n"
          "1D scattering solution assuming a plane-parallel atmosphere (flat\n"
          "Earth). It calculates up to two Stokes parameters (*stokes_dim*<=2),\n"
          "i.e., all azimuthally randomly oriented particles are allowed (this\n"
          "also includes macroscopically isotropic particles). Refraction is\n"
          "not taken into account.\n"
          "\n"
          "The scattering solution is internally obtained over the full\n"
          "(plane-parallel) atmosphere, i.e. not confined to the cloudbox.\n"
          "However, the radiation field output is limited to the cloudbox.\n"
          "This allows to consider clearsky RT through a non-spherical\n"
          "atmosphere outside the cloudbox improving the RT solution for\n"
          "non-plane-parallel media compared to the plain RT4 output.\n"
          "\n"
          "*nstreams* is the number of polar angles taken into account\n"
          "internally in the scattering solution. That is, *nstreams*\n"
          "determines the angular resolution, hence the accuracy, of the\n"
          "scattering solution. The more anisotropic the bulk scattering\n"
          "matrix, the more streams are required. The computational burden\n"
          "increases approximately with the third power of *nstreams*.\n"
          "The default value (*nstreams*=16) was found to be sufficient for\n"
          "most microwave scattering calculations. It is likely insufficient\n"
          "for IR calculations involving ice clouds, though.\n"
          "\n"
          "Here, *scat_za_grid* is NOT an input parameter, but output, and its\n"
          "size equals *nstreams* or *nstreams*+2 (Gauss-Legendre and Double\n"
          "Gauss quadratures in case *add_straight_angles*=1) (the reason is\n"
          "that the computational burden is high for additional angles,\n"
          "regardless whether they are quadrature angles or not; hence the\n"
          "quadrature angles supplemented with 0 and 180deg are considered to\n"
          "provide the best radiation field for a given effort).\n"
          "\n"
          "The *auto_inc_streams* feature can be used to increase the number\n"
          "of streams used internally in the scattering solution when found\n"
          "necessary.\n"
          "NOTE: this number-of-streams increase is only internally - the\n"
          "angular dimension of the output *doit_i_field* is fixed to the\n"
          "*nstreams* given as input to this WSM.\n"
          "\n"
          "Quadrature methods available are: 'L'obatto, 'G'auss-Legendre and\n"
          "'D'ouble Gauss quadrature.\n"
          "\n"
          "This WSM applies *surface_rtprop_agenda* to derive reflection\n"
          "matrix and surface emission vector that are directly feed into\n"
          "RT4's core solver (instead of their RT4-internal calculation as\n"
          "used by *RT4CalcWithRT4Surface*).\n"
          "\n"
          "Known issues of ARTS implementation:\n"
          "- Surface altitude is not an interface parameter. Surface is\n"
          "  implicitly assumed to be at the lowest atmospheric level.\n"
          "- TOA incoming radiation is so far assumed as blackbody cosmic\n"
          "  background (temperature taken from the ARTS-internal constant).\n"
          "- *pfct_method* 'interpolate' currently not implemented here.\n"
          "\n"
          "The keyword *pfct_method* allows to choose the method to extract the\n"
          "scattering matrix. 'interpolate' considers temperature dependence,\n"
          "others neglect it by chosing one specific temperature grid point\n"
          "from the single scattering data: 'low' choses the lowest T-point,\n"
          "'high' the highest T-point, and 'median' the median T-point. As\n"
          "different scattering elements can have different temperature grids,\n"
          "the actual temperature value used can differ between the scattering\n"
          "elements.\n"
          "Note that this keyword solely affects the scattering matrix;\n"
          "extinction matrix and absorption vector are always interpolated to\n"
          "the actual temperature.\n"),
      AUTHORS("Jana Mendrok"),
      OUT("doit_i_field", "scat_za_grid", "scat_aa_grid"),
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
         "atmosphere_dim",
         "pnd_field",
         "t_field",
         "z_field",
         "vmr_field",
         "p_grid",
         "scat_data",
         "f_grid",
         "stokes_dim"),
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
      GIN_DEFAULT(
          "16", "median", "D", "1", "19", "0", "0", "1", "0", "1e-6"),
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
               " *auto_inc_nstreams* gives the maximum number of streams to"
               " increase to. Note that the output *doit_i_field* remains"
               " with angular dimension of *nstreams*, only the internal"
               " solution is adapted (and then interpolated to the"
               " lower-resolution output angular grid).",
               "For *auto_inc_nstreams*>0, flag whether to not fail even if"
               " scattering matrix norm is not preserved when maximum stream"
               " number is reached. Internal RT4 calculations is then"
               " performed with nstreams=*auto_inc_nstreams*.",
               "For *auto_inc_nstreams*>0, polar angle interpolation order"
               " for interpolation from internal increased stream to"
               " originally requested nstreams-ifield.",
               "For *auto_inc_nstreams*>0, flag whether to do polar angle"
               " interpolation in cosine (='mu') space.",
               "Maximum optical depth of infinitesimal layer (where single"
               " scattering approximation is assumed to apply).")));

  md_data_raw.push_back(MdRecord(
      NAME("RT4CalcWithRT4Surface"),
      DESCRIPTION(
          "As RT4Calc except for using RT4's proprietary surface type handling.\n"
          "\n"
          "This WSM is only indented for testing purposes.\n"
          "\n"
          "The following surface type/property methods are available and\n"
          "require the the following input:\n"
          "- 'L'ambertian: *surface_scalar_reflectivity*, *surface_skin_t*\n"
          "- 'F'resnel: *surface_complex_refr_index*, *surface_skin_t*\n"
          "- 'S'pecular: *surface_reflectivity*, *surface_skin_t*\n"
          "'L' and 'F' use proprietary RT4 methods, 'S' uses RT4's Fresnel\n"
          "methods modified to behave similar to ARTS'\n"
          "*surfaceFlatReflectivity*.\n"),
      AUTHORS("Jana Mendrok"),
      OUT("doit_i_field", "scat_za_grid", "scat_aa_grid"),
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
         "atmosphere_dim",
         "pnd_field",
         "t_field",
         "z_field",
         "vmr_field",
         "p_grid",
         "scat_data",
         "f_grid",
         "stokes_dim",
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
               " *auto_inc_nstreams* gives the maximum number of streams to"
               " increase to.",
               "For *auto_inc_nstreams*>0, flag whether to not fail even if"
               " scattering matrix norm is not preserved when maximum stream"
               " number is reached. Internal RT4 calculations is then"
               " performed with nstreams=*auto_inc_nstreams*.",
               "For *auto_inc_nstreams*>0, polar angle interpolation order"
               " for interpolation from internal increased stream to"
               " originally requested nstreams-ifield.",
               "For *auto_inc_nstreams*>0, flag whether to do polar angle"
               " interpolation in cosine (='mu') space.",
               "Maximum optical depth of infinitesimal layer (where single"
               " scattering approximation is assumed to apply).")));

  md_data_raw.push_back(MdRecord(
      NAME("RT4Test"),
      DESCRIPTION(
          "RT4 validation test.\n"
          "\n"
          "Executes test case testc shipped with PolRadTran/RT4 code (but uses\n"
          "data files converted to arts-xml). Output written to (xml-)file.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("rte_losGeometricFromRtePosToRtePos2"),
      DESCRIPTION(
          "The geometric line-of-sight between two points.\n"
          "\n"
          "The method sets *rte_los* to the line-of-sight, at *rte_pos*,\n"
          "that matches the geometrical propagation path between *rte_pos*\n"
          "and *rte_pos2*.\n"
          "\n"
          "The standard case should be that *rte_pos2* corresponds to a\n"
          "transmitter, and *rte_pos* to the receiver/sensor.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("rte_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "lat_grid",
         "lon_grid",
         "refellipsoid",
         "rte_pos",
         "rte_pos2"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(
      MdRecord(NAME("rte_losSet"),
               DESCRIPTION("Sets *rte_los* to the given angles.\n"
                           "\n"
                           "The azimuth angle is ignored for 1D and 2D.\n"),
               AUTHORS("Patrick Eriksson"),
               OUT("rte_los"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("atmosphere_dim"),
               GIN("za", "aa"),
               GIN_TYPE("Numeric", "Numeric"),
               GIN_DEFAULT(NODEF, NODEF),
               GIN_DESC("Zenith angle of sensor line-of-sight.",
                        "Azimuth angle of sensor line-of-sight.")));

  md_data_raw.push_back(MdRecord(
      NAME("rte_posSet"),
      DESCRIPTION(
          "Sets *rte_pos* to the given co-ordinates.\n"
          "\n"
          "The longitude is ignored for 1D and 2D, and the latitude is also \n"
          "ignored for 1D.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("rte_pos"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim"),
      GIN("z", "lat", "lon"),
      GIN_TYPE("Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Geometrical altitude of sensor position.",
               "Latitude of sensor position.",
               "Longitude of sensor position.")));

  md_data_raw.push_back(MdRecord(
      NAME("rte_pos_losMoveToStartOfPpath"),
      DESCRIPTION(
          "Sets *rte_pos* and *rte_los* to values for last point in *ppath*.\n"
          "\n"
          "For example, if the propagation path intersects with the surface,\n"
          "this method gives you the position and angle of *ppath* at the\n"
          "surface.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("rte_pos", "rte_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim", "ppath"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("ScatElementsPndAndScatAdd"),
      DESCRIPTION(
          "Adds single scattering data and particle number density for\n"
          "individual scattering elements.\n"
          "\n"
          "The methods reads the specified files and appends the obtained data\n"
          "to *scat_data* and *pnd_field_raw*. Scattering data is appended to\n"
          "the current last existing scattering species in *scat_data*.\n"),
      AUTHORS("Claudia Emde, Jana Mendrok"),
      OUT("scat_data_raw", "pnd_field_raw"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data_raw", "pnd_field_raw", "atmosphere_dim"),
      GIN("scat_data_files", "pnd_field_files"),
      GIN_TYPE("ArrayOfString", "ArrayOfString"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("List of names of single scattering data files.",
               "List of names of the corresponding pnd_field files.")));

  md_data_raw.push_back(MdRecord(
      NAME("ScatElementsSelect"),
      DESCRIPTION(
          "Allows to limit considered scattering elements according to size.\n"
          "\n"
          "Scattering elements of a specified scattering species are removed\n"
          "from *scat_data_raw* and *scat_meta*, i.e. removed from further\n"
          "calculations, if their particle size exceeds the specified limits.\n"
          "Specification of the scattering species is done by name matching the\n"
          "scattering species name part of *scat_species* tag.\n"
          "As size parameter, all size parameters reported by the meta data\n"
          "can be used (see *scat_meta_single* for offered parameters and\n"
          "their naming).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("ScatElementsToabs_speciesAdd"),
      DESCRIPTION(
          "Appends scattering elements to *abs_species* including reading\n"
          "single scattering data and corresponding pnd field.\n"
          "\n"
          "The methods reads the specified single scattering and pnd_field\n"
          "data of individual scattering elements and appends the obtained data\n"
          "to *scat_data* (appending to its last scattering species) and\n"
          "*vmr_field_raw*. Per scattering element, it also appends one\n"
          "instance of species 'particles' to *abs_species*.\n"),
      AUTHORS("Jana Mendrok"),
      OUT("scat_data_raw",
          "vmr_field_raw",
          "abs_species",
          "propmat_clearsky_agenda_checked",
          "abs_xsec_agenda_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data_raw",
         "vmr_field_raw",
         "abs_species",
         "propmat_clearsky_agenda_checked",
         "abs_xsec_agenda_checked",
         "atmosphere_dim",
         "f_grid"),
      GIN("scat_data_files", "pnd_field_files"),
      GIN_TYPE("ArrayOfString", "ArrayOfString"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("List of names of single scattering data files.",
               "List of names of the corresponding pnd_field files.")));

  md_data_raw.push_back(MdRecord(
      NAME("ScatSpeciesExtendTemperature"),
      DESCRIPTION(
          "Extends valid temperature range of single scattering data.\n"
          "\n"
          "The method allows to extend the temperature range of given single\n"
          "scattering data by duplicating optical property data at the low\n"
          "and/or high limits of the associated temperature grid. *T_low* and\n"
          "*T_high* specify the temperature grid points that are added.\n"
          "Extension is only performed if *T_low* is lower and *T_high* is\n"
          "higher than the original lowest and highest temperatures,\n"
          "respectively, and if the original data contains more than one\n"
          "temperature grid point (i.e., when not assumed constant anyways).\n"
          "\n"
          "The method is thought, e.g., for atmospheric ice falling into\n"
          "atmospheric layers with temperatures above the melting point of\n"
          "ice, where ambient and particle temperature deviate (as long as\n"
          "frozen the ice temperature remains at the melting point\n"
          "temperature). It is not internally checked, whether the original\n"
          "data includes the melting point.\n"
          "The method can be used in a wider sense. However, it remains in the\n"
          "responsibility of the user to apply the method in a meaningful\n"
          "sense and on meaningful single scattering data.\n"
          "\n"
          "The temperature extension is applied on all scattering elements of\n"
          "a scattering species. If *scat_species* is defined, *species* can\n"
          "be used to select the species on which the extension shall be\n"
          "applied comparing *species* with the scattering species name part\n"
          "of *scat_species*. If no *species* is specified, the method is\n"
          "applied on the current last existing scattering species in\n"
          "*scat_data*. Through the latter the method can be applied for cases\n"
          "when *scat_species* is not defined (e.g. when *pnd_field* data is\n"
          "created externally instead of from hydrometeor fields using\n"
          "*pnd_fieldCalcFromscat_speciesFields*).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("ScatSpeciesInit"),
      DESCRIPTION(
          "Initializes the scattering species related data variables.\n"
          "\n"
          "This method initializes the *scat_species* WSV, the variables that\n"
          "will hold the raw optical properties and the raw particle number\n"
          "distributions of the scattering elements (*scat_data_raw* and\n"
          "*pnd_field_raw*, respectively) as well as the one holding the meta\n"
          "information about the scattering elements (*scat_meta*).\n"
          "\n"
          "This method has to be executed before WSM reading/adding to the\n"
          "said variable, e.g. before *ScatSpeciesPndAndScatAdd*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("ScatSpeciesMerge"),
      DESCRIPTION(
          "Merges single scattering data of all scattering elements into one\n"
          "element of bulk properties.\n"
          "\n"
          "Before entering the scattering solver, this method prepares the\n"
          "effective bulk single scattering properties of all scattering\n"
          "elements. Done by calculating the particle number density weighted\n"
          "sum of the single scattering properties of all scattering elements\n"
          "per pressure level. Accordingly, *pnd_field* is resized to\n"
          "[np, np, 1, 1], where np is the number of pressure levels inside\n"
          "the cloudbox. The diagonal elements of the new *pnd_field* are set\n"
          "to 1, all others to 0. *scat_data* is resized to np. Each new\n"
          "scattering element represents the weighted sum of all particles at\n"
          "one presssure level.\n"
          "\n"
          "The method also adapts *scat_species* and *scat_meta* such that\n"
          "they remain consistent with *pnd_field* and can pass\n"
          "*cloudbox_checkedCalc*.\n"
          "\n"
          "The method is suggested to be called directly after\n"
          "*pnd_fieldCalcFromscat_speciesFields* (but also after\n"
          "*cloudbox_checkedCalc*).\n"
          "Its purpose is to speed up the scattering calculations.\n"
          "\n"
          "This is an experimental method currently only working for limited\n"
          "cases. All scattering elements must be of the same ptype and must\n"
          "share the same *f_grid*, *za_grid*, and *aa_grid*. That is, the\n"
          "scattering matrix, extinction matrix, and absorption vector of all\n"
          "scattering elements must have the same dimensions. No interpolation\n"
          "(apart from temperature) is performed.\n"
          "\n"
          "This method can only be used with a 1D atmosphere.\n"),
      AUTHORS("Oliver Lemke"),
      OUT("pnd_field",
          "scat_data",
          "scat_meta",
          "scat_species",
          "cloudbox_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("pnd_field",
         "scat_data",
         "scat_meta",
         "scat_species",
         "cloudbox_checked",
         "atmosphere_dim",
         "cloudbox_on",
         "cloudbox_limits",
         "t_field",
         "z_field",
         "z_surface"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("ScatSpeciesPndAndScatAdd"),
      DESCRIPTION(
          "Adds single scattering data and particle number densities for one\n"
          "scattering species.\n"
          "\n"
          "The WSV *pnd_field_raw* containing particle number densities for\n"
          "all scattering species can be generated outside ARTS, for example\n"
          "by using PyARTS or atmlab. This method reads this data as well as\n"
          "its corresponding single scattering data, which is added as a new\n"
          "scattering species to *scat_data*.\n"
          "This method needs as input an ArrayOfString holding the filenames\n"
          "of the single scattering data for each scattering element and a\n"
          "file containing the corresponding *pnd_field_raw*. In contrast to\n"
          "the scattering data, the pnd-fields are stored in a single XML-file\n"
          "containing an ArrayofGriddedField3, i.e. holding the pnd-field data\n"
          "of all scattering elements.\n"
          "\n"
          "Important note:\n"
          "The order of the filenames for the scattering data files has to\n"
          "correspond to the order of the pnd-fields, stored in the variable\n"
          "*pnd_field_raw*.\n"),
      AUTHORS("Claudia Emde, Jana Mendrok"),
      OUT("scat_data_raw", "pnd_field_raw"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("scat_data_raw", "pnd_field_raw", "atmosphere_dim"),
      GIN("scat_data_files", "pnd_fieldarray_file"),
      GIN_TYPE("ArrayOfString", "String"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC(
          "Array of names of files containing the single scattering data.",
          "Name of file holding the corresponding array of pnd_field data.")));

  md_data_raw.push_back(MdRecord(
      NAME("ScatSpeciesScatAndMetaRead"),
      DESCRIPTION(
          "Reads single scattering data and scattering meta data for one\n"
          "scattering species.\n"
          "\n"
          "This method takes a string array as input containing the location\n"
          "(path and filename) of the single scattering data. Location of\n"
          "corresponding scattering meta data is derived applying a naming\n"
          "convention: ending '.xml*' is replaced by '.meta.xml' (search for\n"
          "zipped files is done automatically).\n"
          "\n"
          "All scattering elements read in one call of the method are assigned\n"
          "to one and the same scattering species. That is, reading in data for\n"
          "a bunch of scattering species can be realized by multiple call of\n"
          "this method. Assignment to scattering species is in the order of the\n"
          "calls (i.e., first method call reads data for first *scat_species*\n"
          "entry, second call for second scat_species entry and so on).\n"
          "Note that no two scattering elements of the same scattering species\n"
          "are allowed to be equal in size; this will result in an error in\n"
          "*pnd_fieldCalcFromscat_speciesFields*\n"
          "\n"
          "Important note:\n"
          "The order of the filenames for the single scattering data files has to\n"
          "exactly correspond to the order of the scattering meta data files.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("scat_data_singleTmatrix"),
      DESCRIPTION(
          "A basic interface to Mishchenko's T-matrix code linked to ARTS.\n"
          "\n"
          "The method performs a T-matrix calculation for a single scattering\n"
          "element, i.e. a combination of particle shape, size, aspect ratio\n"
          "and orientation.\n"
          "\n"
          "Particle shape (*shape*) has two options:\n"
          "   \"spheroidal\" and \"cylindrical\"\n"
          "\n"
          "Particle size (*diameter_volume_equ*) is given as the equivalent\n"
          "volume sphere diameter. That is, the diameter obtained if all the\n"
          "particle's material is rearranged into a (solid) sphere.\n"
          "\n"
          "Particle aspect ratio ar (*aspect_ratio*) is a numeric value, defined\n"
          "according to Mishchenko's definition as ratio of horizontal axis a to\n"
          "vertical (rotational) axis b: ar=a/b. That is, oblates have ar>1,\n"
          "prolates ar<1.\n"
          "Perfect spheres (spheroidals with ar=1) can trigger numerical issues.\n"
          "To avoid these, we internally increase their aspect ratio by 1e-6,\n"
          "i.e. turning perfect spheres into very light oblates.\n"
          "\n"
          "Particle type (*ptype*) has two options:\n"
          "   \"totally_random\" and \"azimuthally_random\"\n"
          "For totally randomly oriented particles, *data_aa_grid* is not taken\n"
          "into account (but a Vector type container needs to be passed).\n"
          "\n"
          "For further information on how aspect ratio and the different shapes\n"
          "and orientations are defined, see the documentation of the T-matrix\n"
          "code found http://www.giss.nasa.gov/staff/mmishchenko/t_matrix.html\n"
          "\n"
          "Regarding *ndgs*, we refer to the this comment from the documentation:\n"
          "   \"Parameter controlling the number of division points\n"
          "   in computing integrals over the particle surface.\n"
          "   For compact particles, the recommended value is 2.\n"
          "   For highly aspherical particles larger values (3, 4,...)\n"
          "   may be necessary to obtain convergence.\n"
          "   The code does not check convergence over this parameter.\n"
          "   Therefore, control comparisons of results obtained with\n"
          "   different NDGS-values are recommended.\"\n"),
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
    ( MdRecord
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

  md_data_raw.push_back(MdRecord(
      NAME("scat_data_checkedCalc"),
      DESCRIPTION(
          "Checks dimensions, grids and single scattering properties of all\n"
          "scattering elements in *scat_data*.\n"
          "\n"
          "Dimension and grid equirements:\n"
          "- The scattering element's f_grid is either identical to *f_grid* or\n"
          "  of dimension 1.\n"
          "- In the latter case, the scattering element's f_grid value must\n"
          "  not deviate from any of the *f_grid* values by more than a\n"
          "  fraction of *dfrel_threshold*.\n"
          "- The frequency dimension of pha_mat_data, ext_mat_data, and\n"
          "  abs_vec_data is either equal to the scattering element's f_grid\n"
          "  or 1.\n"
          "- The temperature dimension of pha_mat_data, ext_mat_data, and\n"
          "  abs_vec_data is either equal to the scattering element's T_grid\n"
          "  or 1.\n"
          "- The temperature dimension of ext_mat_data, and abs_vec_data is\n"
          "  identical.\n"
          "\n"
          "The single scattering property contents are checked using\n"
          "*scat_dataCheck*. For details, see there. The depth of these checks\n"
          "and their rigour can adapted (see description of parameters\n"
          "*check_level* and *sca_mat_threshold* in *scat_dataCheck*) or can\n"
          "be skipped entirely (setting *check_level* to 'none').\n"
          "NOTE: These test shall only be skipped when one is confident that\n"
          "the data is correct, e.g. by having run *scat_dataCheck* on the set\n"
          "of data before, e.g. in a separate ARTS run.\n"),
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
               "See *check_level* in *scat_dataCheck*.",
               "See *sca_mat_threshold* in *scat_dataCheck*.")));

  md_data_raw.push_back(MdRecord(
      NAME("scat_data_monoCalc"),
      DESCRIPTION(
          "Interpolates *scat_data* by frequency to give *scat_data_mono*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("scat_data_monoExtract"),
      DESCRIPTION(
          "Extracts data at *f_index* from *scat_data* to give *scat_data_mono*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("scat_dataCalc"),
      DESCRIPTION(
          "Prepares *scat_data* for the scattering solver.\n"
          "\n"
          "Derives single scattering data for the frequencies given by\n"
          "*f_grid* by interpolation from *scat_data_raw*. *f_grid* should be\n"
          "the actual WSV *f_grid* or a single-element Vector.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("scat_dataCheck"),
      DESCRIPTION(
          "Method for checking the validity and consistency of the single\n"
          "scattering properties in *scat_data*.\n"
          "\n"
          "It checks that *scat_data* does not contain any invalid values,\n"
          "that is any NaN elements in K, Z, or a or any negative values in\n"
          "the 'scalar' properties K11, Z11, and a1.\n"
          "\n"
          "When *check_type* is 'all', it is furthermore checked that the\n"
          "scattering matrix is properly normalized, that is that the solid\n"
          "sphere integrated scattering matrix (int_Z11), which is supposed to\n"
          "be normalized to the scattering cross section, is sufficiently\n"
          "consistent with the scattering cross section (C_sca) derived from\n"
          "the difference of extinction (K11) and absorption (a1):\n"
          "int_z11 ~ C_sca = K11-a1.\n"
          "Sufficient consistency is defined by the maximum allowed deviation\n"
          "in single scattering albedo, *sca_mat_threshold*, testing for\n"
          "  ( <int_Z11>/<C_sca>-1. ) * ( <C_sca>/<K11> ) <= sca_mat_threshold.\n"
          "The check is skipped if *check_type* is 'sane'.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("scat_dataReduceT"),
      DESCRIPTION(
          "Reduces temperature dimension of single scattering to a single entry.\n"
          "\n"
          "FIXME...\n"
          "Derives single scattering data for the frequencies given by\n"
          "*f_grid* by interpolation from *scat_data*. *f_grid* should be\n"
          "the actual WSV *f_grid* or a single-element Vector.\n"),
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

  /*
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "scat_dataFromMeta" ),
        DESCRIPTION
        (
         "This workspace method calculates scattering data and adds it to\n"
         "*scat_data* using the *scat_meta* data.\n"
         "The scattering data is calculated with the T-matrix method.\n"
         "\n"
         "One set of scattering data is calculated for each scattering element\n"
         "in *scat_meta*.\n"
         ),
        AUTHORS( "Johan Strandgren, Oliver Lemke" ),
        OUT("scat_data"),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN("scat_meta"),
        GIN( "za_grid", "aa_grid", "precision" ),
        GIN_TYPE("Vector", "Vector", "Numeric" ),
        GIN_DEFAULT(NODEF, NODEF, NODEF ),
        GIN_DESC("Zenith angle grid",
                 "Azimuth angle grid",
                 "Precision"
                 )
        ));
*/

  md_data_raw.push_back(MdRecord(
      NAME("scat_speciesSet"),
      DESCRIPTION(
          "Sets the WSV *scat_species*."
          "\n"
          "With this function, the user specifies settings for the particle\n"
          "number density calculations by *pnd_fieldCalcFromscat_speciesFields*.\n"
          "The input is an ArrayOfString that needs to be in a specific format.\n"
          "For details, see WSV *scat_species*.\n"
          "\n"
          "*Example:*  ['IWC-MH97', 'LWC-H98_STCO']\n"),
      AUTHORS("Daniel Kreyling"),
      OUT("scat_species"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("scat_tags", "delim"),
      GIN_TYPE("ArrayOfString", "String"),
      GIN_DEFAULT(NODEF, "-"),
      GIN_DESC("Array of pnd calculation parameters.",
               "Delimiter string of *scat_species* elements.")));

  md_data_raw.push_back(MdRecord(
      NAME("ScatSpeciesSizeMassInfo"),
      DESCRIPTION(
          "Derives size and mass information for a scattering species."
          "\n"
          "This method assumes that the mass-size relationship can described\n"
          "by *scat_species_a* and *scat_species_b*. See documentation of \n"
          "*scat_species_a* for details.\n"
          "\n"
          "The quantity to be used as size descriptor is here denoted as x, and\n"
          "is selected by setting *x_unit*. The options are:\n"
          " \"dveq\" : The size grid is set to scat_meta.diameter_volume_equ\n"
          " \"dmax\" : The size grid is set to scat_meta.diameter_max\n"
          " \"area\" : The size grid is set to scat_meta.diameter_area_equ_aerodynamical\n"
          " \"mass\" : The size grid is set to scat_meta.mass\n"
          "This selection determines *scat_species_x*.\n"
          "\n"
          "The parameters *scat_species_a* and *scat_species_b* are determined by\n"
          "a numeric fit between *scat_species_x* and corresponding masses in\n"
          "*scat_meta*. This fit is performed over sizes inside the range\n"
          "[x_fit_start,x_fit_end]. This range is allowed to be broader than\n"
          "the coverage of *scat_species_x*. There must be at least two sizes\n"
          "inside [x_fit_start,x_fit_end].\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("particle_fieldCleanup"),
      DESCRIPTION(
          "Removes unrealistically small or erroneous data from particle fields.\n"
          "\n"
          "This WSM checks if the input particle field (e.g.\n"
          "*particle_bulkprop_field*, *scat_species_XXX_field) contains values\n"
          "smaller than the given *threshold*. In this case, these values will\n"
          "be set to zero.\n"
          "\n"
          "The method should be applied if the particle fields contain\n"
          "unrealistically small or erroneous data (NWP/GCM model data, e.g.\n"
          "from the Chevallierl_91l sets, often contain very small or even\n"
          "negative values, which are numerical artefacts rather than physical\n"
          "values.)\n"
          "For the scat_species_XXX_fields, it needs to be applied separately\n"
          "per Tensor4 type field collection. This allows to use different\n"
          "thresholds for the different types of fields (not for the different\n"
          "scattering species, though).\n"
          "\n"
          "*particle_fieldCleanup* shall be called after generation of the\n"
          "atmopheric fields.\n"),
      AUTHORS("Daniel Kreyling"),
      OUT(),
      GOUT("particle_field_out"),
      GOUT_TYPE("Tensor4"),
      GOUT_DESC("A particle property field, e.g. *particle_bulkprop_field*"
                " or *scat_species_mass_density_field*"),
      IN(),
      GIN("particle_field_in", "threshold"),
      GIN_TYPE("Tensor4", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("A particle property field, e.g. *particle_bulkprop_field* or"
               " *scat_species_mass_density_field*",
               "Threshold below which the *particle_field* values are set to"
               " zero.")));

  md_data_raw.push_back(MdRecord(
      NAME("Select"),
      DESCRIPTION(
          "Method to select some elements from one array and copy them to\n"
          "a new array. (Works also for vectors.)\n"
          "\n"
          "This works also for higher dimensional objects, where the selection is\n"
          "always performed in the first dimension.\n"
          "\n"
          "If needleindexes is set to [-1], all elements are copied."
          "\n"
          "For example:\n"
          "\n"
          "Select(y,x,[0,3])\n"
          "\n"
          "will select the first and fourth row of matrix x and copy them to the\n"
          "output matrix y.\n"
          "\n"
          "Note that it is even safe to use this method if needles and haystack\n"
          "are the same variable.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("needles"),
      GOUT_TYPE(ARRAY_GROUPS + ", Vector, Matrix, Sparse"),
      GOUT_DESC("Selected elements. Must have the same variable type as "
                "haystack."),
      IN(),
      GIN("haystack", "needleindexes"),
      GIN_TYPE(ARRAY_GROUPS + ", Vector, Matrix, Sparse", "ArrayOfIndex"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Variable to select from. May be the same variable as needles.",
               "The elements to select (zero based indexing, as always.)"),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true)));

  md_data_raw.push_back(MdRecord(
      NAME("sensor_checkedCalc"),
      DESCRIPTION(
          "Checks consistency of the sensor variables.\n"
          "\n"
          "The following WSVs are examined: *f_grid*, *sensor_pos*, *sensor_los*,\n"
          "*transmitter_pos*, *mblock_dlos_grid*, *antenna_dim*,\n"
          "*sensor_response*, *sensor_response_f*, *sensor_response_pol*,\n"
          "and *sensor_response_dlos*.\n"
          "\n"
          "If any of these variables are changed, then this method shall be\n"
          "called again (no automatic check that this is fulfilled!).\n"
          "\n"
          "The main tests are that dimensions of sensor variables agree\n"
          "with other settings, e.g., the size of f_grid, atmosphere_dim,\n"
          "stokes_dim, etc.\n"
          "\n"
          "If any test fails, there is an error. Otherwise, *sensor_checked*\n"
          "is set to 1.\n"),
      AUTHORS("Jana Mendrok"),
      OUT("sensor_checked"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "stokes_dim",
         "f_grid",
         "sensor_pos",
         "sensor_los",
         "transmitter_pos",
         "mblock_dlos_grid",
         "sensor_response",
         "sensor_response_f",
         "sensor_response_pol",
         "sensor_response_dlos"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("sensorOff"),
      DESCRIPTION(
          "Sets sensor WSVs to obtain monochromatic pencil beam values.\n"
          "\n"
          "The variables are set as follows:\n"
          "   mblock_dlos_grid        : One row with zero(s).\n"
          "   sensor_response*        : As returned by *sensor_responseInit*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_response",
          "sensor_response_f",
          "sensor_response_pol",
          "sensor_response_dlos",
          "sensor_response_f_grid",
          "sensor_response_pol_grid",
          "sensor_response_dlos_grid",
          "mblock_dlos_grid"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("stokes_dim", "f_grid"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("sensor_losGeometricFromSensorPosToOtherPositions"),
      DESCRIPTION(
          "The geometric line-of-sight between pair of points.\n"
          "\n"
          "The method sets *sensor_los* to the line-of-sights, that matches the\n"
          "geometrical propagation path from *sensor_pos* to *target_pos*. This\n"
          "is done for pair of positions, i.e. the two matrices shall have the same\n"
          "number of rows. The number of columns in *target_pos* shall be two for\n"
          "1D and 2D and two for 3D, exactly as for *rte_pos2*.\n"
          "\n"
          "See also *rte_losGeometricFromRtePosToRtePos2*. This method calls that\n"
          "method for each pair of positions, where values in *sensor_pos* matches\n"
          "*rte_pos and values in *target_pos* matches *rte_pos2*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_los"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "lat_grid",
         "lon_grid",
         "refellipsoid",
         "sensor_pos"),
      GIN("target_pos"),
      GIN_TYPE("Matrix"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Target position, for each position in *sensor_pos*.")));

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseAntenna"),
      DESCRIPTION(
          "Includes response of the antenna.\n"
          "\n"
          "The function returns the sensor response matrix after the antenna\n"
          "characteristics have been included.\n"
          "\n"
          "The function handles \"multi-beam\" cases where the polarisation\n"
          "coordinate system is the same for all beams.\n"
          "\n"
          "See *antenna_dim*, *antenna_dlos* and *antenna_response* for\n"
          "details on how to specify the antenna response.\n"
          "\n"
          "One dimensional antenna patterns are handled as other response\n"
          "functions. That is. both antenna repsonse and radiances are treated\n"
          "piece-wise linear functions, and the pencil beam calculations must\n"
          "cover the full sensor reponse (i.e. *mblock_dlos_grid* must be\n"
          "sufficiently braod).\n"
          "\n"
          "On the other hand, 2D antenna patterns are so far handled in a\n"
          "simplified manner. For 2D, the antenna pattern is simply sampled at\n"
          "the points specified by *mblock_dlos_grid* and each pencil beam\n"
          "direction is considered to represent the same size in terms of solid\n"
          "beam ansgle. In addition, there is no check at all on how well\n"
          "*mblock_dlos_grid* covers the antenna response.\n"),
      AUTHORS("Mattias Ekstrom", "Patrick Eriksson"),
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
         "atmosphere_dim",
         "antenna_dim",
         "antenna_dlos",
         "antenna_response",
         "sensor_norm"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseBackend"),
      DESCRIPTION(
          "Includes response of the backend (spectrometer).\n"
          "\n"
          "The function returns the sensor response matrix after the backend\n"
          "characteristics have been included.\n"
          "\n"
          "See *f_backend*, *backend_channel_response* and *sensor_norm* for\n"
          "details on how to specify the backend response.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseBackendFrequencySwitching"),
      DESCRIPTION(
          "Frequency switching for a pure SSB reciever.\n"
          "\n"
          "This function can be used for simulation of frequency switching.\n"
          "That is, when the final spectrum is the difference of two spectra\n"
          "shifted in frequency. The switching is performed by the LO, but\n"
          "for a pure singel sideband reciever this is most easily simulated\n"
          "by instead shifting the backend, as done here.\n"
          "\n"
          "A strightforward frequency switching is modelled (no folding)\n"
          "The channel positions for the first measurement cycle are\n"
          "f_backend+df1, and for the second f_backend+df2. The first\n"
          "measurement cycle is given the negive weight. That is, the output\n"
          "is the spectrum for cycle2 minus the spectrum for cycle1.\n"
          "Output frequency grids are set to *f_backend*.\n"
          "\n"
          "Use *sensor_responseFrequencySwitching* for double sideband cases.\n"
          "\n"
          "The method has the same general functionality as, and can replace,\n"
          "*sensor_responseBackend*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseBeamSwitching"),
      DESCRIPTION(
          "Simulation of \"beam switching\".\n"
          "\n"
          "The measurement procedure is based on taking the difference between\n"
          "two spectra measured in different directions, and the calculation\n"
          "set-up must treat exactly two observation directions.\n"
          "\n"
          "The returned spectrum is y = w1*y + w2*y2, where y1 and w1 are the\n"
          "spectrum and weight for the first direction, respectively (y2 and\n"
          "(w2 defined correspondingly for the second direction).\n"
          "\n"
          "Zenith and azimuth angles after beam switching are set to the\n"
          "values of the second direction.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseFillFgrid"),
      DESCRIPTION(
          "Polynomial frequency interpolation of spectra.\n"
          "\n"
          "The sensor response methods treat the spectra to be piece-wise linear\n"
          "functions. This method is a workaround for making methods handling\n"
          "the spectra in a more elaborate way: it generates spectra on a more\n"
          "dense grid by polynomial interpolation. The interpolation is not\n"
          "done explicitly, it is incorporated into *sensor_response*.\n"
          "\n"
          "This method should in general increase the calculation accuracy for\n"
          "a given *f_grid*. However, the selection of (original) grid points\n"
          "becomes more sensitive when using this method. A poor choice of grid\n"
          "points can result in a decreased accuracy, or generation of negative\n"
          "radiances. Test calculations indicated that the error easily can\n"
          "increase with this method close the edge of *f_grid*, and it could\n"
          "be wise to make *f_grid* a bit wider than actually necessary to avoid\n"
          "this effect\n"
          "\n"
          "The method shall be inserted before the antenna stage. That is, this\n"
          "method shall normally be called directly after *sensor_responseInit*.\n"
          "\n"
          "Between each neighbouring points of *f_grid*, this method adds\n"
          "*nfill* grid points. The polynomial order of the interpolation is\n"
          "*polyorder*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseFrequencySwitching"),
      DESCRIPTION(
          "Simulation of \"frequency switching\".\n"
          "\n"
          "A general method for frequency switching. The WSM\n"
          "*sensor_responseBackendFrequencySwitching* gives a description of\n"
          "this observation technique, and is also a more straightforward\n"
          " method for pure singel sideband cases.\n"
          "\n"
          "It is here assume that *sensor_responseMultiMixerBackend* has been\n"
          "used to calculate the spectrum for two LO positions. This method\n"
          "calculates the difference between these two spectra, where the\n"
          "second spectrum gets weight 1 and the first weight -1 (as in\n"
          "*sensor_responseBackendFrequencySwitching*).\n"
          "\n"
          "Output frequency grids are taken from the second spectrum.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseIF2RF"),
      DESCRIPTION(
          "Converts sensor response variables from IF to RF.\n"
          "\n"
          "The function converts intermediate frequencies (IF) in\n"
          "*sensor_response_f* and *sensor_response_f_grid* to radio\n"
          "frequencies (RF). This conversion is needed if the frequency\n"
          "translation of a mixer is included and the position of backend\n"
          "channels are specified in RF.\n"
          "\n"
          "A direct frequency conversion is performed. Values are not\n"
          "sorted in any way.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseInit"),
      DESCRIPTION(
          "Initialises the variables summarising the sensor response.\n"
          "\n"
          "This method sets the variables to match monochromatic pencil beam\n"
          "calculations, to be further modified by inclusion of sensor\n"
          "characteristics. Use *sensorOff* if pure monochromatic pencil\n"
          "beam calculations shall be performed.\n"
          "\n"
          "The variables are set as follows:\n"
          "   sensor_response : Identity matrix, with size matching *f_grid*,\n"
          "                     *stokes_dim* and *mblock_dlos_grid*.\n"
          "   sensor_response_f       : Repeated values of *f_grid*.\n"
          "   sensor_response_pol     : Data matching *stokes_dim*.\n"
          "   sensor_response_dlos    : Repeated values of *mblock_dlos_grid*.\n"
          "   sensor_response_f_grid  : Equal to *f_grid*.\n"
          "   sensor_response_pol_grid: Set to 1:*stokes_dim*.\n"
          "   sensor_response_dlos_grid : Equal to *mblock_dlos_grid*.\n"),
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
         "mblock_dlos_grid",
         "antenna_dim",
         "atmosphere_dim",
         "stokes_dim",
         "sensor_norm"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseMetMM"),
      DESCRIPTION(
          "Sensor setup for meteorological millimeter instruments.\n"
          "\n"
          "This method is handy if you are simulating a passband-type instrument,\n"
          "consisting of a few discrete channels.\n"
          "\n"
          "For flexibility, the Met-MM system is seperated in two calculation\n"
          "steps. To fully use the system, create *f_grid* (and some associated\n"
          "variables) by *f_gridMetMM* before calling this method. However, it is\n"
          "possible to use this method with any *f_grid*, as long as matching\n"
          "*f_backend*, *channel2fgrid_indexes* and *channel2fgrid_weights*\n"
          "are provided.\n"
          "\n"
          "Each scan sequence is treated as a measurement block. *sensor_pos* is\n"
          "set in the standard way. The number of rows in *sensor_pos* determines the\n"
          "number of scan sequences that will be simulated. On the other hand,\n"
          "*sensor_los* is handled in a special way. All zenith angles must be set\n"
          "to 180 deg. For 3D, the given azimuth angles are taken as the direction\n"
          "of scanning, where the azimuth angle is defined with respect to North\n"
          "in standard manner. For example, if the scanning happens to move from\n"
          "SW to NE, the azimuth angle should be set to 45 deg. The angles of the\n"
          "scanning sequence are taken from *antenna_dlos*. This WSV is here only\n"
          "allowed to have a single column, holding relative zenith angles. For\n"
          "3D, the azimuth angles in *antenna_dlos* are hard-coded to zero. As\n"
          "zenith angles in *sensor_los* are locked to 180 deg, *antenna_dlos*\n"
          "effectively holds the nadir angles. These angles can be both positive or\n"
          "negative, where the recommended choice is to operate with negative\n"
          "to end up with final zenith angles between 0 and 180 deg.\n"
          "\n"
          "The method does not support 2D atmospheres (across-track scanning is\n"
          "inconsistent with 2D). For simpler switching between 1D and 3D,\n"
          "the argument *mirror_dza* is at hand. It can only be used for 3D.\n"
          "If set to true, the zenith angles in *antenna_dlos* are mapped\n"
          "to also cover the other side of the swath and the simulations will\n"
          "cover both sides of the swath.\n"),
      AUTHORS("Oliver Lemke", "Patrick Eriksson"),
      OUT("antenna_dim",
          "mblock_dlos_grid",
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
      IN("atmosphere_dim",
         "stokes_dim",
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

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseMixer"),
      DESCRIPTION(
          "Includes response of the mixer of a heterodyne system.\n"
          "\n"
          "The function returns the sensor response matrix after the mixer\n"
          "characteristics have been included. Frequency variables are\n"
          "converted from radio frequency (RF) to intermediate frequency (IF).\n"
          "The returned frequency grid covers the range [0,max_if], where\n"
          "max_if is the highest IF covered by the sideband response grid.\n"
          "\n"
          "See *lo* and *sideband_response* for details on how to specify the\n"
          "mixer response\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseMixerBackendPrecalcWeights"),
      DESCRIPTION(
          "Includes pre-calculated response covering mixer and backend.\n"
          "\n"
          "This method acts similar to *sensor_responseBackend*, but uses\n"
          "pre-calculated weights. These weights can also include the effect\n"
          "of mixer and sideband filtering.\n"
          "\n"
          "As usual, *f_backend* gives the frequency of the channels. This WSM\n"
          "has no direct influence on the result, but at least representative\n"
          "values must be set.\n"
          "\n"
          "The frequency response is defined using *channel2fgrid_indexes* and\n"
          "*channel2fgrid_weights*.\n"
          "\n"
          "Both *frequency_index* and *frequency_weights* are assumed to be common\n"
          "for all viewing directions.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseMultiMixerBackend"),
      DESCRIPTION(
          "Handles mixer and backend parts for an instrument having multiple\n"
          "mixer chains.\n"
          "\n"
          "The WSMs *sensor_responseMixer*, *sensor_responseIF2RF* and\n"
          "*sensor_responseBackend* are called for each mixer chain, and a\n"
          "complete *sensor_response* is assembled. The instrument responses\n"
          "are described by *lo_multi*, *sideband_response_multi*,\n"
          "*sideband_mode_multi*, *f_backend_multi* and\n"
          "*backend_channel_response_multi*. All these WSVs must have same\n"
          "vector or array length. As *sensor_responseIF2RF* is called,\n"
          "*f_backend_multi* must hold RF (not IF) and output frequencies\n"
          "will be in absolute frequency (RF).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responsePolarisation"),
      DESCRIPTION(
          "Extraction of non-default polarisation components.\n"
          "\n"
          "The default is to output the Stokes elements I, Q, U and V (up to\n"
          "*stokes_dim*). This method allows to change the \"polarisation\" of\n"
          "the output. Polarisation components to be extracted are selected by\n"
          "*instrument_pol*. This method can be applied at any step of the sensor\n"
          "matrix set-up.\n"
          "\n"
          "The method can only be applied on data for I, Q, U and V. The value\n"
          "of *stokes_dim* must be sufficiently large for the selected\n"
          "components. For example, I+45 requires that *stokes_dim* is at\n"
          "least 3. \n"
          "\n"
          "See *instrument_pol* for coding of polarisation states.\n"
          "\n"
          "Note that the state of *iy_unit* is considered. This WSV must give\n"
          "the actual unit of the data. This as, the extraction of components\n"
          "is slightly different if data are radiances or brightness\n"
          "temperatures.  In practise this means that *iy_unit* (as to be\n"
          "applied inside *iy_main_agenda*) must be set before calling this\n"
          "method.\n"),
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
         "stokes_dim",
         "iy_unit",
         "instrument_pol"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseStokesRotation"),
      DESCRIPTION(
          "Includes a rotation of the Stokes H and V directions.\n"
          "\n"
          "The method applies the rotations implied by *stokes_rotation*.\n"
          "See the description of that WSV for details.\n"
          "\n"
          "This method does not change the size of *sensor_response*, and\n"
          "the auxiliary variables (sensor_response_f etc.) are not changed.\n"
          "\n"
          "To apply the method, *stokes_dim* must be >= 3. The complete effect\n"
          "of the rotation can not be determibed with lower *stokes_dim*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("sensor_response"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("sensor_response",
         "sensor_response_f_grid",
         "sensor_response_pol_grid",
         "sensor_response_dlos_grid",
         "stokes_dim",
         "stokes_rotation"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseSimpleAMSU"),
      DESCRIPTION(
          "Simplified sensor setup for an AMSU-type instrument.\n"
          "\n"
          "This method allows quick and simple definition of AMSU-type\n"
          "sensors. Assumptions:\n"
          "\n"
          "1. Pencil beam antenna.\n"
          "2. Double sideband receivers.\n"
          "3. Sideband mode \"upper\"\n"
          "4. The channel response is rectangular.\n"
          "\n"
          "Under these assumptions the only inputs needed are the LO positions,\n"
          "the offsets from the LO, and the IF bandwidths. They are provieded\n"
          "in sensor_description_amsu.\n"),
      AUTHORS("Stefan Buehler"),
      OUT("f_grid",
          "antenna_dim",
          "mblock_dlos_grid",
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
      IN("atmosphere_dim", "stokes_dim", "sensor_description_amsu"),
      GIN("spacing"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(".1e9"),
      GIN_DESC("Desired grid spacing in Hz.")));

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseGenericAMSU"),
      DESCRIPTION(
          "Simplified sensor setup for an AMSU-type instrument.\n"
          "\n"
          "This function is derived from 'sensor_responseSimpleAMSU' \n"
          "but is more generalized since the number of passbands in each \n"
          "can be in the range from 1 to 4 - in order to correctly simulate\n"
          "AMSU-A type sensors \n"
          "\n"
          "This method allows quick and simple definition of AMSU-type\n"
          "sensors. Assumptions:\n"
          "\n"
          "1. Pencil beam antenna.\n"
          "2. 1-4 Passband/sidebands per channel.\n"
          "3. Sideband mode \"upper\"\n"
          "4. The channel response is rectangular.\n"
          "\n"
          "Under these assumptions the only inputs needed are the LO positions,\n"
          "the offsets from the LO, and the IF bandwidths. They are provided\n"
          "in sensor_description_amsu.\n"),
      AUTHORS("Oscar Isoz"),
      OUT("f_grid",
          "antenna_dim",
          "mblock_dlos_grid",
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
      IN("atmosphere_dim", "stokes_dim", "sensor_description_amsu"),
      GIN("spacing"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(".1e9"),
      GIN_DESC("Desired grid spacing in Hz.")));

  md_data_raw.push_back(MdRecord(
      NAME("sensor_responseWMRF"),
      DESCRIPTION(
          "Adds WMRF weights to sensor response.\n"
          "\n"
          "This method adds a spectrometer response that has been calculated\n"
          "with the weighted mean of representative frequencies (WMRF) method. It\n"
          "consists of a set of selected frequencies, and associated weights.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("SetLineMixingCoefficinetsFromRelmat"),
      DESCRIPTION("Sets Rosenkranz coefficients for the line by non-linear\n"
                  "regression from values computed by relmat\n"),
      AUTHORS("Richard Larsson"),
      OUT("abs_lines_per_band", "relmat_per_band"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_band",
         "abs_species_per_band",
         "band_identifiers",
         "abs_species",
         "isotopologue_ratios",
         "partition_functions",
         "rtp_pressure",
         "abs_t",
         "relmat_type_per_band",
         "wigner_initialized"),
      GIN("pressure_rule_limit",
          "debug",
          "order_of_linemixing",
          "use_adiabatic_factor"),
      GIN_TYPE("Numeric", "Index", "Index", "Index"),
      GIN_DEFAULT("0.1", "0", "1", "1"),
      GIN_DESC(
          "Limit when perturbation theory is assumed to work",
          "Lets relmat know it is to print debug information if true.",
          "Choice of order of linemixing",
          "Truth-value: if we should use the precomputed adiabatic factors")));

  md_data_raw.push_back(
      MdRecord(NAME("SetRelaxationMatrixCalcType"),
               DESCRIPTION("Sets *relmat_type_per_band* (see for types).\n"
                           "Set input to a 1-long array for the same type\n"
                           "for all relaxation matrices.\n"),
               AUTHORS("Richard Larsson"),
               OUT("relmat_type_per_band"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("abs_lines_per_band"),
               GIN("type"),
               GIN_TYPE("ArrayOfIndex"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Type of relaxation matrix calculations")));

  md_data_raw.push_back(MdRecord(
      NAME("SetBandIdentifiersAuto"),
      DESCRIPTION("Sets *band_identifiers* to multiple O2-66 and\n"
                  "CO2-* bands.  This is not an exhaustive list so\n"
                  "if you need more bands to identify, please add them\n"
                  "to the list or set *band_identifiers* manually.\n"),
      AUTHORS("Richard Larsson"),
      OUT("band_identifiers"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("SetBandIdentifiersFromLines"),
      DESCRIPTION(
          "Sets (not adds to) *band_identifiers* to all uniques in *abs_lines*\n"
          "\n"
          "Will set all Quantum Numbers that are in the line as undefined\n"
          "if they do not exist in *band_quantums* before comparing the resulting\n"
          "QuantumIdentifier to the ones already in *band_identifiers*.  If new, push_back\n"
          "is called adding the ID.  The functions does not check the species or isotopologue,\n"
          "but you might have to define one anyways for the *band_quantums* input to work.\n"
          "\n"
          "Example *band_quantums*=\"O2-66 EN v1 0\" means all quantum numbers that\n"
          "are not v1 will be removed, and *band_identifiers* will consist of all the\n"
          "identifiers of combinations of v1, e.g., v1 0 to v1 1, v1 0 to v1 0, and so on.\n"),
      AUTHORS("Richard Larsson"),
      OUT("band_identifiers"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines"),
      GIN("band_quantums"),
      GIN_TYPE("QuantumIdentifier"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Quantum numbers that are defined for the band")));

  md_data_raw.push_back(MdRecord(
      NAME("SparseSparseMultiply"),
      DESCRIPTION(
          "Multiplies a Sparse with another Sparse, result stored in Sparse.\n"
          "\n"
          "Makes the calculation: out = m1 * m2\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Sparse"),
      GOUT_DESC("Product, can be same variable as any of the inputs."),
      IN(),
      GIN("m1", "m2"),
      GIN_TYPE("Sparse", "Sparse"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Left sparse matrix.", "Right sparse matrix.")));

  md_data_raw.push_back(MdRecord(
      NAME("SparseMatrixIdentity"),
      DESCRIPTION(
          "Returns a sparse dentity matrix.\n"
          "\n"
          "The size of the matrix created is n x n. Default is to return a\n"
          "true identity matrix (I), but you can also select another value\n"
          "along the diagonal be setting *value*. That is, the output is\n"
          "value*I.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Sparse"),
      GOUT_DESC("Sparse output matrix"),
      IN(),
      GIN("n", "value"),
      GIN_TYPE("Index", "Numeric"),
      GIN_DEFAULT(NODEF, "1"),
      GIN_DESC("Size of the matrix", "The value along the diagonal.")));

  md_data_raw.push_back(MdRecord(
      NAME("DiagonalMatrix"),
      DESCRIPTION(
          "Create a diagonal matrix from a vector."
          "\n"
          "This creates a dense or sparse diagonal matrix with the elements of the given vector\n"
          " on the diagonal.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Matrix, Sparse"),
      GOUT_DESC("The diagonal matrix"),
      IN(),
      GIN("v"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The vector containing the diagonal elements.")));

  md_data_raw.push_back(MdRecord(
      NAME("spectral_irradiance_fieldFromiyField"),
      DESCRIPTION(
          "Calculates the spectral irradiance from *doit_i_field* .\n"
          "by integrating over the angular grids according to the grids set\n"
          "by *AngularGridsSetForFluxCalc* \n"
          "See *AngularGridsSetForFluxCalc to set \n"
          "*scat_za_grid, scat_aa_grid, and za_grid_weights*\n"),
      AUTHORS("Manfred Brath"),
      OUT("spectral_irradiance_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("doit_i_field", "scat_za_grid", "scat_aa_grid", "za_grid_weights"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("specular_losCalc"),
      DESCRIPTION(
          "Calculates the specular direction of surface reflections.\n"
          "\n"
          "A help method to set up the surface properties. This method\n"
          "calculates *specular_los*, that is required in several methods\n"
          "to convert zenith angles to incidence angles.\n"
          "\n"
          "The method also returns the line-of-sight matching the surface normal.\n"
          "\n"
          "The default is to consider the surface slope when calculating the\n"
          "specular direction. That is, the variation of *z_surface* (as well as\n"
          "the geoid radius) is considered and the specular direction is calculated\n"
          "including the specified topography. This part can be deactivated by\n"
          "setting *ignore_surface_slope* to 1. In this case, the zenith angle of\n"
          "the specular direction is simply 180-rtp_los[0]. *ignore_surface_slope*\n"
          "has only an effect for 2D and 3D, as 1D implies a constant radius of\n"
          "the surface (i.e. no topography).\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("specular_los", "surface_normal"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("rtp_pos",
         "rtp_los",
         "atmosphere_dim",
         "lat_grid",
         "lon_grid",
         "refellipsoid",
         "z_surface"),
      GIN("ignore_surface_slope"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("Flag to control if surface slope is consdered or not.")));

  md_data_raw.push_back(MdRecord(
      NAME("specular_losCalcNoTopography"),
      DESCRIPTION(
          "Calculates the specular direction of surface reflections for horisontal\n"
          "surfaces.\n"
          "\n"
          "In contrast to *specular_losCalc*, this method ignores the topography\n"
          "implied by *z_surface*. That is, any slope of the surface is ignored.\n"
          "\n"
          "The typical application of this WSM should be water surfaces (lakes and\n"
          "oceans).\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("specular_los", "surface_normal"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("rtp_pos", "rtp_los", "atmosphere_dim"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("StringJoin"),
      DESCRIPTION(
          "Concatenate two or more strings.\n"
          "\n"
          "The output string is overwritten, but is allowed to appear\n"
          "in the input list. Up to 10 strings can be concatenated at once.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
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

  md_data_raw.push_back(
      MdRecord(NAME("StringSet"),
               DESCRIPTION("Sets a String to the given text string.\n"),
               AUTHORS("Patrick Eriksson"),
               OUT(),
               GOUT("out"),
               GOUT_TYPE("String"),
               GOUT_DESC("Variable to initialize."),
               IN(),
               GIN("text"),
               GIN_TYPE("String"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Input text string."),
               SETMETHOD(true)));

  md_data_raw.push_back(MdRecord(
      NAME("z_surfaceFromFileAndGrid"),
      DESCRIPTION(
          "Sets the surface altitude for a given latitude and longitude grid.\n"),
      AUTHORS("Richard Larsson"),
      OUT("z_surface"),
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

  md_data_raw.push_back(MdRecord(
      NAME("z_surfaceConstantAltitude"),
      DESCRIPTION(
          "Sets the surface altitude to a constant. Defaults to zero.\n"),
      AUTHORS("Richard Larsson"),
      OUT("z_surface"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("lat_grid", "lon_grid"),
      GIN("altitude"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT("0"),
      GIN_DESC("The constant altitude.")));

  md_data_raw.push_back(MdRecord(
      NAME("surfaceBlackbody"),
      DESCRIPTION(
          "Creates variables to mimic a blackbody surface.\n"
          "\n"
          "This method sets up *surface_los*, *surface_rmatrix* and\n"
          "*surface_emission* for *surface_rtprop_agenda*. Here, *surface_los*\n"
          "and *surface_rmatrix* are set to be empty, and *surface_emission*\n"
          "to hold blackbody radiation for a temperature of *surface_skin_t*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "f_grid",
         "stokes_dim",
         "rtp_pos",
         "rtp_los",
         "surface_skin_t"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("surfaceFastem"),
      DESCRIPTION(
          "Usage of FASTEM together with MC and DOIT.\n"
          "\n"
          "The recommended way to use FASTEM is by *iySurfaceFastem*, but that\n"
          "is not always possible, such as when using MC and DOIT. This is the\n"
          "case as those scattering methods use *surface_rtprop_agenda*,\n"
          "while *iySurfaceFastem* fits with *iy_surface_agenda*. This WSM solves\n"
          "this by allowing FASTEM to be used inside *surface_rtprop_agenda*.\n"
          "\n"
          "However, FASTEM is here used in an approximative way. For a correct\n"
          "usage of FASTEM, the atmospheric transmittance shall be calculated\n"
          "for the position and direction of concern, but this is not possible\n"
          "together with DOIT and MC. Instead, the transmittance is an input\n"
          "to the method, and must either be pre-calculated or set to a\n"
          "representative value.\n"
          "\n"
          "See *iySurfaceFastem*, for further details on the special input\n"
          "arguments.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "stokes_dim",
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

  md_data_raw.push_back(MdRecord(
      NAME("surfaceTessem"),
      DESCRIPTION(
          "TESSEM sea surface microwave emissivity parametrization.\n"
          "\n"
          "This method computes surface emissivity and reflectivity matrices for\n"
          "ocean surfaces using the TESSEM emissivity model: Prigent, C., et al.\n"
          "Sea‐surface emissivity parametrization from microwaves to millimetre\n"
          "waves, QJRMS, 2017, 143.702: 596-605.\n"
          "\n"
          "The validity range of the parametrization of is 10 to 700 GHz, but for\n"
          "some extra flexibility frequencies between 5 and 900 GHz are accepted.\n"
          "The accepted temperaute range for *surface_skin_t* is [260.0 K, 373.0 K]\n"
          "\n"
          "The model itself is represented by the neural networks in\n"
          "*tessem_neth* and *tessem_netv*.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "stokes_dim",
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

  md_data_raw.push_back(MdRecord(
      NAME("surfaceTelsem"),
      DESCRIPTION(
          "Compute surface emissivities using the TELSEM 2 model.\n"
          "\n"
          "This method uses second version of the TELSEM model for calculating\n"
          "land surface emissivities (F. Aires et al, \"A Tool to Estimate \n"
          " Land‐Surface Emissivities at Microwave frequencies (TELSEM) for use\n"
          " in numerical weather prediction\" Quarterly Journal of the Royal\n"
          "Meteorological Society, vol. 137, (656), pp. 690-699, 2011.)\n"
          "This methods computes land surface emissivities for a given pencil beam\n"
          "using a given TELSEM2 atlas.\n"
          "The input must satisfy the following conditions, otherwise an error is thrown:\n"
          " - The input frequencies (*f_grid*) must be within the range [5 GHz, 900 GHz]\n"
          " - The skin temperature (*surface_skin_t*) must be within the range\n"
          "   [180 K, 360 K]\n"
          "\n"
          "A TELSEM atlas contains only suface emissivities for locations that are\n"
          "classified as land. By default this WSM will throw an error if the\n"
          "pencil beam hits the surface at a position that is not contained in the\n"
          "given atlas.\n"
          "\n"
          "The above behavior can be avoided by setting *d_max* to a positive value.\n"
          "This enables nearest neighbor interpolation, which assigns the emissivities\n"
          "of the nearest found cell in the atlas to the given position. In this case,\n"
          "an error is only thrown if the distance of the found neighbor is higher\n"
          "than the provided value of *d_max.\n"
          "\n"
          "You can limit the final reflectivity applied by setting *r_min* and *r_max*.\n"
          "\n"
          "To extract a land-sea mask from a given telsem atlas see the WSM\n"
          "*telsemSurfaceTypeLandSea*.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "stokes_dim",
         "f_grid",
         "lat_grid",
         "lat_true",
         "lon_true",
         "rtp_pos",
         "rtp_los",
         "surface_skin_t"),
      GIN("atlas", "r_min", "r_max", "d_max"),
      GIN_TYPE("TelsemAtlas", "Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "0", "1", "-1.0"),
      GIN_DESC("The Telsem atlas to use for the emissivity calculation.",
               "Minimum allowed value for reflectivity to apply.",
               "Maximum allowed value for reflectivity to apply.",
               "Maximum allowed distance in meters for nearest neighbor"
               " interpolation. Set to a negative value or zero to disable"
               " interpolation.")));

  md_data_raw.push_back(MdRecord(
      NAME("surfaceFlatRefractiveIndex"),
      DESCRIPTION(
          "Creates variables to mimic specular reflection by a (flat) surface\n"
          "where the complex refractive index is specified.\n"
          "\n"
          "The dielectric properties of the surface are described by\n"
          "*surface_complex_refr_index*. The Fresnel equations are used to\n"
          "calculate amplitude reflection coefficients. The method can thus\n"
          "result in that the reflection properties differ between frequencies\n"
          "and polarisations.\n"
          "\n"
          "Local thermodynamic equilibrium is assumed, which corresponds to\n"
          "that the reflection and emission coefficients add up to 1.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         "stokes_dim",
         "atmosphere_dim",
         "rtp_pos",
         "rtp_los",
         "specular_los",
         "surface_skin_t",
         "surface_complex_refr_index"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("surfaceFlatReflectivity"),
      DESCRIPTION(
          "Creates variables to mimic specular reflection by a (flat) surface\n"
          "where *surface_reflectivity* is specified.\n"
          "\n"
          "Works basically as *surfaceFlatScalarReflectivity* but is more\n"
          "general as vector radiative transfer is more properly handled. See\n"
          "the ARTS theory document (ATD) for details around how\n"
          "*surface_emission* is determined. In the nomenclature of ATD,\n"
          "*surface_reflectivity* gives R.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         "stokes_dim",
         "atmosphere_dim",
         "rtp_pos",
         "rtp_los",
         "specular_los",
         "surface_skin_t",
         "surface_reflectivity"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("surfaceFlatRvRh"),
      DESCRIPTION(
          "Creates variables to mimic specular reflection by a (flat) surface\n"
          "where *surface_rv_rh* is specified.\n"
          "\n"
          "This method assumes that the reflection at vertical and horizontal\n"
          "polarisation differs. As power reflection coefficients are provided\n"
          "there is no information at hand on phase shifts between polarisations,\n"
          "and they are simply assumed to be zero. These assumptions result in\n"
          "that *surface_emission* is set to zero for positions corresponding to\n"
          "U and V, and that all diagonal elementsof  *surface_rmatrix* are equal\n"
          "(the mean of rv and rh). Further, all off-diagonal elements of\n"
          "*surface_rmatrix* are all zero except for (0,1) and (1,0).\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         "stokes_dim",
         "atmosphere_dim",
         "rtp_pos",
         "rtp_los",
         "specular_los",
         "surface_skin_t",
         "surface_rv_rh"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("surfaceFlatScalarReflectivity"),
      DESCRIPTION(
          "Creates variables to mimic specular reflection by a (flat) surface\n"
          "where *surface_scalar_reflectivity* is specified.\n"
          "\n"
          "This method assumes that the reflection at vertical and horizontal\n"
          "polarisation is identical. This assumption includes that there is no\n"
          "phase shift between polarisations. These assumptions result in that\n"
          "*surface_emission* is set to zero for positions corresponding to Q,\n"
          "U and V, and that *surface_rmatrix* becomes a diagonal matrix (with\n"
          "all elements on the diagonal equal to the specified reflectivity).\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         "stokes_dim",
         "atmosphere_dim",
         "rtp_pos",
         "rtp_los",
         "specular_los",
         "surface_skin_t",
         "surface_scalar_reflectivity"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("surfaceLambertianSimple"),
      DESCRIPTION(
          "Creates variables to mimic a Lambertian surface.\n"
          "\n"
          "A Lambertian surface can be characterised solely by its\n"
          "reflectivity, here taken from *surface_scalar_reflectivity*.\n"
          "\n"
          "The down-welling radiation field is estimated by making calculations\n"
          "for *lambertian_nza* directions. The range of zenith angles ([0,90])\n"
          "is divided in an equidistant manner for 1D. For 2D and 3D see below.\n"
          "The values for *surface_rmatrix* are assuming a constant radiance\n"
          "over each zenith angle range. See AUG.\n"
          "\n"
          "Default is to select the zenith angles for *sensor_los* to be placed\n"
          "centrally in the grid ranges. For example, if *lambertian_nza* is set\n"
          "to 9, down-welling radiation will be calculated for zenith angles = \n"
          "5, 15, ..., 85. The position of these angles can be shifted by\n"
          "*za_pos*. This variable specifies the fractional distance inside the\n"
          "ranges. For example, a *za_pos* of 0.7 (np still 9) gives the angles\n"
          "7, 17, ..., 87.\n"
          "\n"
          "Only upper-left diagonal element of the *surface_rmatrix* is\n"
          "non-zero. That is, the upwelling radiation is always unpolarised.\n"
          "\n"
          "Local thermodynamic equilibrium is assumed, which corresponds to\n"
          "that the reflection and emission coefficients \"add up to 1\".\n"
          "\n"
          "For 2D and 3D, the down-welling directions are placed along the\n"
          "the viewing direction, e.g. for 3D the azimuth angle is kept constant.\n"
          "In 2D and 3D surface topography can exist, and to avoid getting views\n"
          "going directly into the surface, angels are not distributed over 90 deg,\n"
          "but 90-abs(surface_normal[0]).\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix", "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         "stokes_dim",
         "atmosphere_dim",
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

  md_data_raw.push_back(MdRecord(
      NAME("surfaceSemiSpecularBy3beams"),
      DESCRIPTION(
          "A simplistic treatment of semi-specular surfaces.\n"
          "\n"
          "This method has no strong physical basis but could be used for simpler\n"
          "testing or as starting point for more advanced methods.\n"
          "\n"
          "This method assumes that the surface can be treated to have three facets,\n"
          "all lacking surface roughness. One facet is assumed to give standard\n"
          "specular reflection, while the two other facets are tilted with +dza\n"
          "and -dza, respectively. The tilt is assumed to only affect the zenith\n"
          "angle of the reflected direction (azimuth same as for specular direction).\n"
          "The area ratio of the non-tilted facet is set by *specular_factor*.\n"
          "That is, the specular beam is given weight w, while the other two beams\n"
          "each get weight (1-w)/2.\n"
          "\n"
          "If a facet tilts away from the viewing direction in such way that\n"
          "the surface is observed from below, the tilt of the facet is decreased\n"
          "in steps of 1 degree until a successful calculation is obtained. If this\n"
          "turns out to require a tilt of zero, this facete is merged with\n"
          "the specular direction.\n"
          "\n"
          "The pure specular properties of the surface shall be described by\n"
          "*surface_rtprop_sub_agenda*. That is, if you have specular surface\n"
          "described and you want to make it semi-specular by this method, you\n"
          "move the content of the existing *surface_rtprop_agenda* to\n"
          "*surface_rtprop_sub_agenda* and instead fill *surface_rtprop_agenda*\n"
          "with this method.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_skin_t",
          "surface_los",
          "surface_rmatrix",
          "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "f_grid",
         "rtp_pos",
         "rtp_los",
         "surface_rtprop_sub_agenda"),
      GIN("specular_factor", "dza"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The weight given to the specular direction. Denoted as w above."
               " A value between 1/3 and 1.",
               "Zenith angle seperation to each secondary direction. A "
               "between 0 and 45 degrees.")));

  md_data_raw.push_back(MdRecord(
      NAME("surfaceSplitSpecularTo3beams"),
      DESCRIPTION(
          "A very simple approximation of a semi-specular surface.\n"
          "\n"
          "This method has no direct physical basis but could be used for simpler\n"
          "testing or as starting point for more advanced methods.\n"
          "\n"
          "The method requires that the surface RT properties (e.g. *surface_los*)\n"
          "have been set up to mimic a specular surface. This method splits the down-\n"
          "welling radiation into three directions. The specular direction is given\n"
          "weight w, while the other two beams each get weight (1-w)/2. The basic\n"
          "polarised reflectivity from the specular calculations is maintained\n"
          "for each beam. The beams are just separated in zenith angle, with a\n"
          "separation of *dza*. The lowermost beam is not allowed to be closer to\n"
          "the surface than 1 degree. If there is no room for the lowermost beam,\n"
          "it is merged with the main beam.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_los", "surface_rmatrix"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("surface_los",
         "surface_rmatrix",
         "atmosphere_dim",
         "rtp_pos",
         "rtp_los"),
      GIN("specular_factor", "dza"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The weight given to the specular direction. Denoted as w above."
               " A value between 1/3 and 1.",
               "Zenith angle seperation to each secondary direction. A "
               "between 0 and 45 degrees.")));

  md_data_raw.push_back(MdRecord(
      NAME("surface_complex_refr_indexFromGriddedField5"),
      DESCRIPTION(
          "Extracts complex refractive index from a field of such data.\n"
          "\n"
          "The method allows to obtain *surface_complex_refr_index* by\n"
          "interpolation of a geographical field of such data. The position\n"
          "for which refraction shall be extracted is given by *rtp_pos*.\n"
          "The refractive index field is expected to be stored as:\n"
          "   GriddedField5:\n"
          "      Vector f_grid[N_f]\n"
          "      Vector T_grid[N_T]\n"
          "      ArrayOfString Complex[2]\n"
          "      Vector \"Latitude\"  [N_lat]\n"
          "      Vector \"Longitude\" [N_lon]\n"
          "      Tensor5 data[N_f][N_T][2][N_lat][N_lon]\n"
          "\n"
          "Definition and treatment of the three first dimensions follows\n"
          "*complex_refr_index*, e.g. the temperature grid is allowed\n"
          "to have length 1. The grids for latitude and longitude must have\n"
          "a length of >= 2 (ie. no automatic expansion).\n"
          "\n"
          "Hence, this method performs an interpolation only in the lat and\n"
          "lon dimensions, to a single point. The remaining GriddedField3 is\n"
          "simply returned as *surface_complex_refr_index*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_complex_refr_index"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim", "lat_grid", "lat_true", "lon_true", "rtp_pos"),
      GIN("complex_refr_index_field"),
      GIN_TYPE("GriddedField5"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("A field of complex refractive index.")));

  md_data_raw.push_back(MdRecord(
      NAME("surface_reflectivityFromGriddedField6"),
      DESCRIPTION(
          "Extracts surface reflectivities from a field of such data.\n"
          "\n"
          "This method allows to specify a field of surface reflectivity for\n"
          "automatic interpolation to points of interest. The position and\n"
          "direction for which the reflectivity shall be extracted are given\n"
          "by *rtp_pos* and *rtp_los*. The reflectivity field is expected to\n"
          "be stored as:\n"
          "   GriddedField6:\n"
          "      Vector \"Frequency\"       [N_f]\n"
          "      Vector \"Stokes element\"  [N_s1]\n"
          "      Vector \"Stokes_element\"  [N_s2]\n"
          "      Vector \"Incidence angle\" [N_ia]\n"
          "      Vector \"Latitude\"        [N_lat]\n"
          "      Vector \"Longitude\"       [N_lon]\n"
          "      Tensor6 data[N_f][N_s1][N_s2][N_ia][N_lat][N_lon]\n"
          "\n"
          "Grids for incidence angle, latitude and longitude must have a\n"
          "length of >= 2 (ie. no automatic expansion). If the frequency grid\n"
          "has length 1, this is taken as that the reflectivity is constant,\n"
          "following the definition of *surface_scalar_reflectivity*.\n"
          "The data can cover higher Stokes dimensionalities than set by\n"
          "*stokes_dim*. Data for non-used Stokes elements are just cropped.\n"
          "The order between the two Stokes dimensions is the same as in\n"
          "*surface_reflectivity* and surface_rmatrix*.\n"
          "\n"
          "The interpolation is done in steps:\n"
          "   1: Linear interpolation for lat and lon (std. extrapolation).\n"
          "   2: Interpolation in incidence angle (std. extrapolation).\n"
          "      If the grid has a length of >= 4, cubic interpolation is\n"
          "      applied. Otherwise linear interpolation.\n"
          "   3. Linear interpolation in frequency (if input data have more\n"
          "      than one frequency).\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_reflectivity"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("stokes_dim",
         "f_grid",
         "atmosphere_dim",
         "lat_grid",
         "lat_true",
         "lon_true",
         "rtp_pos",
         "rtp_los"),
      GIN("r_field"),
      GIN_TYPE("GriddedField6"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("A field of surface reflectivities")));

  md_data_raw.push_back(MdRecord(
      NAME("surface_rtpropCallAgendaX"),
      DESCRIPTION(
          "Switch between the elements of *surfacertprop__agenda_array*.\n"
          "\n"
          "This method simply calls the agenda matching *surface_type* and\n"
          "returns the results. That is, the agenda in *surface_rtprop_agenda_array*\n"
          "with index *surface_type* (0-based) is called.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_skin_t",
          "surface_los",
          "surface_rmatrix",
          "surface_emission"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid",
         "rtp_pos",
         "rtp_los",
         "surface_rtprop_agenda_array",
         "surface_type",
         "surface_type_aux"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("surface_scalar_reflectivityFromGriddedField4"),
      DESCRIPTION(
          "Extracts scalar surface reflectivities from a field of such data.\n"
          "\n"
          "This method allows to specify a field of surface reflectivity for\n"
          "automatic interpolation to points of interest. The position and\n"
          "direction for which the reflectivity shall be extracted are given\n"
          "by *rtp_pos* and *rtp_los*. The reflectivity field is expected to\n"
          "be stored as:\n"
          "   GriddedField4:\n"
          "      Vector \"Frequency\"       [N_f]\n"
          "      Vector \"Incidence angle\" [N_ia]\n"
          "      Vector \"Latitude\"        [N_lat]\n"
          "      Vector \"Longitude\"       [N_lon]\n"
          "      Tensor4 data[N_f][N_ia][N_lat][N_lon]\n"
          "\n"
          "Grids for incidence angle, latitude and longitude must have a\n"
          "length of >= 2 (ie. no automatic expansion). If the frequency grid\n"
          "has length 1, this is taken as the reflectivity is constant,\n"
          "following the definition of *surface_scalar_reflectivity*.\n"
          "\n"
          "The interpolation is done in steps:\n"
          "   1: Linear interpolation for lat and lon (std. extrapolation).\n"
          "   2: Interpolation in incidence angle (std. extrapolation).\n"
          "      If the grid has a length of >= 4, cubic interpolation is\n"
          "      applied. Otherwise linear interpolation.\n"
          "   3. Linear interpolation if frequency (if input data have more\n"
          "      than one frequency).\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_scalar_reflectivity"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("stokes_dim",
         "f_grid",
         "atmosphere_dim",
         "lat_grid",
         "lat_true",
         "lon_true",
         "rtp_pos",
         "rtp_los"),
      GIN("r_field"),
      GIN_TYPE("GriddedField4"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("A field of scalar surface reflectivities")));

  /*
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "surface_reflectivityFromSurface_rmatrix" ),
        DESCRIPTION
        (
         "As *surface_scalar_reflectivityFromSurface_rmatrix*, but for (vectorRT)"
         "*surface_reflectivity*.\n"
         "\n"
         "For each frequency f and stokes parameter combi st_in and st_out,"
         "*surface_reflectivity* is set to the sum of"
         "surface_rmatrix(joker,f,st_in,st_out).\n"
        ),
        AUTHORS( "Jana Mendrok" ),
        OUT( "surface_reflectivity" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "surface_rmatrix" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
*/

  md_data_raw.push_back(MdRecord(
      NAME("surface_scalar_reflectivityFromSurface_rmatrix"),
      DESCRIPTION(
          "Sets *surface_scalar_reflectivity* based on *surface_rmatrix*.\n"
          "\n"
          "For each frequency f, *surface_scalar_reflectivity* is set to\n"
          "the sum of surface_rmatrix(joker,f,0,0).\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("surface_typeInterpTypeMask"),
      DESCRIPTION(
          "Closest neighbour interpolation of surface type mask.\n"
          "\n"
          "The method determines the surface type at the position of concern\n"
          "(*rtp_pos*) from the provided type mask (*surface_type_mask*).\n"
          "The closest point in the mask is selected. The surface type\n"
          "is set to the integer part of the value at the found point, while\n"
          "*surface_type_aux* is set to the reminder. For example, if the\n"
          "mask value at closest point is 2.23, *surface_type* is set to 2\n"
          "*surface_type_aux* becomes 0.23.\n"
          "\n"
          "The altitude in *rtp_pos* is ignored.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("surface_type", "surface_type_aux"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "lat_grid",
         "lat_true",
         "lon_true",
         "rtp_pos",
         "surface_type_mask"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("SurfaceDummy"),
      DESCRIPTION(
          "Dummy method for *iy_surface_agenda*.\n"
          "\n"
          "If you don't make use of *surface_props_data* and associated\n"
          "variables, include this method *iy_surface_agenda*. The method\n"
          "just checks that the variables of concern are set to be empty,\n"
          "and you don't need to include calls of *Ignore* and *Touch* in\n"
          "the agenda.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("dsurface_rmatrix_dx", "dsurface_emission_dx"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("dsurface_rmatrix_dx",
         "dsurface_emission_dx",
         "atmosphere_dim",
         "lat_grid",
         "lon_grid",
         "surface_props_data",
         "surface_props_names",
         "dsurface_names",
         "jacobian_do"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("SurfaceFastem"),
      DESCRIPTION(
          "FASTEM sea surface microwave emissivity parametrization.\n"
          "\n"
          "The variable *surface_props_data* must contain these data:\n"
          "  \"Water skin temperature\"\n"
          "  \"Wind speed\"\n"
          "  \"Wind direction\"\n"
          "  \"Salinity\"\n"
          "\n"
          "For some details and comments see *FastemStandAlone* and *surfaceFastem*.\n"),
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
         "stokes_dim",
         "atmosphere_dim",
         "lat_grid",
         "lon_grid",
         "f_grid",
         "rtp_pos",
         "rtp_los",
         "surface_props_data",
         "surface_props_names",
         "dsurface_names",
         "jacobian_do"),
      GIN("transmittance", "fastem_version"),
      GIN_TYPE("Vector", "Index"),
      GIN_DEFAULT(NODEF, "6"),
      GIN_DESC("Transmittance along path of downwelling radiation. A vector "
               "with the same length as *f_grid*.",
               "The version of FASTEM to use.")));

  md_data_raw.push_back(MdRecord(
      NAME("SurfaceTessem"),
      DESCRIPTION(
          "TESSEM sea surface microwave emissivity parametrization.\n"
          "\n"
          "The variable *surface_props_data* must contain these data:\n"
          "  \"Water skin temperature\"\n"
          "  \"Wind speed\"\n"
          "  \"Salinity\"\n"
          "\n"
          "This method computes surface emissivity and reflectivity matrices for\n"
          "ocean surfaces using the TESSEM emissivity model: Prigent, C., et al.\n"
          "Sea-surface emissivity parametrization from microwaves to millimetre\n"
          "waves, QJRMS, 2017, 143.702: 596-605.\n"
          "\n"
          "The validity range of the parametrization of is 10 to 700 GHz, but for\n"
          "some extra flexibility frequencies between 5 and 900 GHz are accepted.\n"
          "The accepted temperaute range for water skin temperature is\n"
          "[260.0 K, 373.0 K]. Salinity shall be in the range [0,1]. That is, a\n"
          "salinity of 3% is given as 0.03.\n"
          "\n"
          "The model itself is represented by the neural networks in\n"
          "*tessem_neth* and *tessem_netv*.\n"),
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
         "stokes_dim",
         "atmosphere_dim",
         "lat_grid",
         "lon_grid",
         "f_grid",
         "rtp_pos",
         "rtp_los",
         "tessem_neth",
         "tessem_netv",
         "surface_props_data",
         "surface_props_names",
         "dsurface_names",
         "jacobian_do"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("TangentPointExtract"),
      DESCRIPTION(
          "Finds the tangent point of a propagation path.\n"
          "\n"
          "The tangent point is here defined as the point with the lowest\n"
          "altitude (which differes from the definition used in the code\n"
          "where it is the point with the lowest radius, or equally the point\n"
          "with a zenith angle of 90 deg.)\n"
          "\n"
          "The tangent point is returned as a vector, with columns matching\n"
          "e.g. *rte_pos*. If the propagation path has no tangent point, the\n"
          "vector is set to NaN.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("tan_pos"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("The position vector of the tangent point."),
      IN("ppath"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("TangentPointPrint"),
      DESCRIPTION(
          "Prints information about the tangent point of a propagation path.\n"
          "\n"
          "The tangent point is here defined as the point with the lowest\n"
          "altitude (which differes from the definition used in the code\n"
          "where it is the point with the lowest radius, or equally the point\n"
          "with a zenith angle of 90 deg.)\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ppath"),
      GIN("level"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Output level to use.")));

  md_data_raw.push_back(MdRecord(
      NAME("telsemStandalone"),
      DESCRIPTION(
          "Stand-alone evaluation of the Telsem model.\n"
          "\n"
          "This evaluates the Telsem land surface emissivity\n"
          "model using the data from the provided atlas.\n"
          "\n"
          "Since TELSEM atlases do not contain data for all locations\n"
          "this function allows for nearest neighbor interpolation, which\n"
          "can be enabled by setting the *d_max* GIN to a positive value.\n"
          "\n"
          "This WSM throws a runtime error if the queried location is not\n"
          "contained in the atlas or the distance of the neighboring cell\n"
          "exceeds the given *d_max* value.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("telsemAtlasLookup"),
      DESCRIPTION(
          "Lookup SSMI emissivities from Telsem Atlas.\n"
          "\n"
          "This returns the emissivities (indices [0,..,6])\n"
          " for the SSMI channels that are contained in\n"
          "the Telsem atlas.\n"
          "\n"
          "If given latitude and longitude are not in the atlas an empty\n"
          "vector is returned.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("telsemSurfaceTypeLandSea"),
      DESCRIPTION(
          "TELSEM based land sea mask.\n"
          "\n"
          "This method determines whether the position in *rtp_pos* is\n"
          "of type ocean or land depending on whether a corresponding\n"
          "cell is contained in the provided TELSEM atlas.\n"
          "In combination with the WSM *surface_rtpropCallAgendaX* this\n"
          "can be used to used different methods to compute surface radiative\n"
          "properties.\n"),
      AUTHORS("Simon Pfreundschuh"),
      OUT("surface_type"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim", "lat_grid", "lat_true", "lon_true", "rtp_pos"),
      GIN("atlas"),
      GIN_TYPE("TelsemAtlas"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The telsem atlas from which to lookup the surface type.")));

  md_data_raw.push_back(MdRecord(
      NAME("telsem_atlasReadAscii"),
      DESCRIPTION(
          "Reads single TELSEM atlas.\n"
          "\n"
          "'directory' needs to contain the original 12 Telsem atlas files\n"
          "and the correlations file. This WSM reads the atlas for the specified\n"
          "month and stores the result in the provided output atlas.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("telsem_atlasesReadAscii"),
      DESCRIPTION(
          "Reads TELSEM atlas files.\n"
          "\n"
          "'directory' needs to contain the original 12 Telsem atlas files\n"
          "and the correlations file.\n"
          "The whole data is combined into the WSV *telsem_atlases*\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("Tensor3AddScalar"),
      DESCRIPTION("Adds a scalar value to all elements of a tensor3.\n"
                  "\n"
                  "The result can either be stored in the same or another\n"
                  "variable.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Tensor3"),
      GOUT_DESC("Output tensor."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Tensor3", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input tensor.", "The value to be added to the tensor.")));

  md_data_raw.push_back(MdRecord(
      NAME("Tensor3ExtractFromTensor4"),
      DESCRIPTION(
          "Extracts a Tensor3 from a Tensor4.\n"
          "\n"
          "Copies book, page, row or column with given Index from input Tensor4\n"
          "variable to output Tensor3.\n"
          "Higher order equivalent of *VectorExtractFromMatrix*.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Tensor3"),
      GOUT_DESC("Extracted tensor."),
      IN(),
      GIN("in", "i", "direction"),
      GIN_TYPE("Tensor4", "Index", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Input Tensor4.",
               "Index of book, page, row or column to extract.",
               "Direction. \"book\" or \"page\" or \"row\" or \"column\".")));

  md_data_raw.push_back(MdRecord(
      NAME("Tensor3Scale"),
      DESCRIPTION("Scales all elements of a tensor with the specified value.\n"
                  "\n"
                  "The result can either be stored in the same or another\n"
                  "variable.\n"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Tensor3"),
      GOUT_DESC("Output tensor."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Tensor3", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input tensor.",
               "The value to be multiplied with the tensor.")));

  md_data_raw.push_back(MdRecord(
      NAME("Tensor3SetConstant"),
      DESCRIPTION(
          "Creates a tensor and sets all elements to the specified value.\n"
          "\n"
          "The size is determined by *ncols*, *nrows* etc.\n"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Tensor3"),
      GOUT_DESC("Variable to initialize."),
      IN("npages", "nrows", "ncols"),
      GIN("value"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Tensor value.")));

  md_data_raw.push_back(MdRecord(
      NAME("Tensor4AddScalar"),
      DESCRIPTION("Adds a scalar value to all elements of a tensor4.\n"
                  "\n"
                  "The result can either be stored in the same or another\n"
                  "variable.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Tensor4"),
      GOUT_DESC("Output tensor."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Tensor4", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input tensor.", "The value to be added to the tensor.")));
  /*
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor4Clip" ),
        DESCRIPTION
        (
         "Clipping of e.g. *vmr_field* and *particle_bulkprop_field*.\n"
         "\n"
         "The method allows you to apply hard limits the values of a\n"
         "Tensor4. The quantati (book dimension) is specified by *iq*.\n"
         "*All values of the quantity below *limit_low*, are simply\n"
         "set to *limit_low*. And the same is performed with respect to\n"
         "*limit_high*. That is, the data in x for the quantity are\n"
         "forced to be inside the range [limit_low,limit_high].\n"
         "\n"
         "Setting iq=-1, is a shortcut for applying the limits on all\n"
         "quantities.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT( "x" ),
        GOUT_TYPE( "Tensor4" ),
        GOUT_DESC( "A Tensor4 holding data, with quantity as book-dimension,"
                   "such as *vmr_field*." ),
        IN(  ),
        GIN( "x", "iq", "limit_low", "limit_high" ),
        GIN_TYPE( "Tensor4", "Index", "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF, NODEF, "-Inf", "Inf" ),
        GIN_DESC( "See GOUT for a defintion.",
                  "Quantity index (zero-based)",
                  "Lower limit for clipping.",
                  "Upper limit for clipping." )
        ));
  */

  md_data_raw.push_back(MdRecord(
      NAME("Tensor4Scale"),
      DESCRIPTION("Scales all elements of a tensor with the specified value.\n"
                  "\n"
                  "The result can either be stored in the same or another\n"
                  "variable.\n"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Tensor4"),
      GOUT_DESC("Output tensor."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Tensor4", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input tensor.",
               "The value to be multiplied with the tensor.")));

  md_data_raw.push_back(MdRecord(
      NAME("Tensor4SetConstant"),
      DESCRIPTION(
          "Creates a tensor and sets all elements to the specified value.\n"
          "\n"
          "The size is determined by *ncols*, *nrows* etc.\n"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Tensor4"),
      GOUT_DESC("Variable to initialize."),
      IN("nbooks", "npages", "nrows", "ncols"),
      GIN("value"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Tensor value.")));

  md_data_raw.push_back(MdRecord(
      NAME("Tensor5Scale"),
      DESCRIPTION("Scales all elements of a tensor with the specified value.\n"
                  "\n"
                  "The result can either be stored in the same or another\n"
                  "variable.\n"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Tensor5"),
      GOUT_DESC("Output tensor."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Tensor5", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input tensor.",
               "The value to be multiplied with the tensor.")));

  md_data_raw.push_back(MdRecord(
      NAME("Tensor5SetConstant"),
      DESCRIPTION(
          "Creates a tensor and sets all elements to the specified value.\n"
          "\n"
          "The size is determined by *ncols*, *nrows* etc.\n"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Tensor5"),
      GOUT_DESC("Variable to initialize."),
      IN("nshelves", "nbooks", "npages", "nrows", "ncols"),
      GIN("value"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Tensor value.")));

  md_data_raw.push_back(MdRecord(
      NAME("Tensor6Scale"),
      DESCRIPTION("Scales all elements of a tensor with the specified value.\n"
                  "\n"
                  "The result can either be stored in the same or another\n"
                  "variable.\n"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Tensor6"),
      GOUT_DESC("Output tensor."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Tensor6", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input tensor.",
               "The value to be multiplied with the tensor.")));

  md_data_raw.push_back(MdRecord(
      NAME("Tensor6SetConstant"),
      DESCRIPTION(
          "Creates a tensor and sets all elements to the specified value.\n"
          "\n"
          "The size is determined by *ncols*, *nrows* etc.\n"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Tensor6"),
      GOUT_DESC("Variable to initialize."),
      IN("nvitrines", "nshelves", "nbooks", "npages", "nrows", "ncols"),
      GIN("value"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Tensor value.")));

  md_data_raw.push_back(MdRecord(
      NAME("Tensor7Scale"),
      DESCRIPTION("Scales all elements of a tensor with the specified value.\n"
                  "\n"
                  "The result can either be stored in the same or another\n"
                  "variable.\n"),
      AUTHORS("Mattias Ekstrom"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Tensor7"),
      GOUT_DESC("Output tensor."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Tensor7", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input tensor.",
               "The value to be multiplied with the tensor.")));

  md_data_raw.push_back(MdRecord(
      NAME("Tensor7SetConstant"),
      DESCRIPTION(
          "Creates a tensor and sets all elements to the specified value.\n"
          "\n"
          "The size is determined by *ncols*, *nrows* etc.\n"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT("out"),
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

  md_data_raw.push_back(MdRecord(
      NAME("TestArrayOfAgenda"),
      DESCRIPTION(
          "A method that is used for the TestArrayOfAgenda test case.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("TessemNNReadAscii"),
      DESCRIPTION(
          "Reads the initialization data for the TESSEM NeuralNet from an ASCII file.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("TestTessem"),
      DESCRIPTION(
          "Example method for TESSEM2.\n"
          "\n"
          "When using the default neural network parameter files\n"
          "from the Tessem 2 distribution, the input Vector should contain\n"
          "5 elements:\n"
          "   - Frequency (10-700) in GHz.\n"
          "   - Theta (0-90) Incidence angle in degrees.\n"
          "   - Windspeed (0-25) at 10m (m/s)\n"
          "     Higher wind speed can be used, but without garantee.\n"
          "   - Surface skin temperature (270-310) in K.\n"
          "   - Salinity (0-0.04) in kg/kg\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("relmat_per_bandInAir"),
      DESCRIPTION("A dummy method to test line mixing.\n"),
      AUTHORS("Richard Larsson"),
      OUT("relmat_per_band"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_lines_per_band",
         "abs_species_per_band",
         "partition_functions",
         "wigner_initialized"),
      GIN("temperatures"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Vector of temperatures to compute the relaxation matrix at")));

  md_data_raw.push_back(MdRecord(
      NAME("Test"),
      DESCRIPTION(
          "A dummy method that can be used for test purposes.\n"
          "\n"
          "This method can be used by ARTS developers to quickly test stuff.\n"
          "The implementation is in file m_general.cc. This just saves you the\n"
          "trouble of adding a dummy method everytime you want to try\n"
          "something out quickly.\n"),
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
  /*
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Test" ),
        DESCRIPTION
        (
         "A dummy method that can be used for test purposes.\n"
         "\n"
         "This method can be used by ARTS developers to quickly test stuff.\n"
         "The implementation is in file m_general.cc. This just saves you the\n"
         "trouble of adding a dummy method everytime you want to try\n"
         "something out quickly.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN("abs_lines_per_species"),
        GIN("aqi"),
        GIN_TYPE("ArrayOfQuantumIdentifier"),
        GIN_DEFAULT(NODEF),
        GIN_DESC("")
        ));
  */

  md_data_raw.push_back(
      MdRecord(NAME("timerStart"),
               DESCRIPTION("Initializes the CPU timer."
                           "\n"
                           "Use *timerStop* to stop the timer.\n"
                           "\n"
                           "Usage example:\n"
                           "   timerStart\n"
                           "   ReadXML(f_grid,\"frequencies.xml\")\n"
                           "   timerStop\n"
                           "   Print(timer)\n"),
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
      MdRecord(NAME("timerStop"),
               DESCRIPTION("Stops the CPU timer."
                           "\n"
                           "See *timerStart* for example usage.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("TMatrixTest"),
      DESCRIPTION(
          "T-Matrix validation test.\n"
          "\n"
          "Executes the standard test included with the T-Matrix Fortran code.\n"
          "Should give the same as running the tmatrix_lp executable in\n"
          "3rdparty/tmatrix/.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("Touch"),
      DESCRIPTION(
          "As *Ignore* but for agenda output.\n"
          "\n"
          "This method is handy for use in agendas in order to suppress\n"
          "warnings about not-produced output workspace variables.\n"
          "\n"
          "What it does, in case the variable is initialized already, is:\n"
          "Nothing!\n"
          "In case the variable is not yet initialized, it is set to NaN.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT("in"),
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

  md_data_raw.push_back(MdRecord(
      NAME("transmittanceFromIy_aux"),
      DESCRIPTION(
          "Creates a vector of transmittance values.\n"
          "\n"
          "The transmittances are set based on optical depths in *iy_aux*. That is,\n"
          "one of the quantities in *iy_aux* must be \"Optical depth\".\n"
          "\n"
          "The created vector has a length matching *f_grid* and can e.g. be used\n"
          "as input to some of the FASTEM methods.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("VectorAddScalar"),
      DESCRIPTION(
          "Adds a scalar to all elements of a vector.\n"
          "\n"
          "The result can either be stored in the same or another vector.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output vector"),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input vector.", "The value to be added to the vector.")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorAddVector"),
      DESCRIPTION(
          "Element-wise addition of two vectors.\n"
          "\n"
          "The method calculates c = a + b.\n"
          "\n"
          "The variable *b* is allowed to have length 1, for any length of\n"
          "*a*. This single value in *b* is then added to every element of *a*.\n"
          "\n"
          "The vectors *a* and *c* can be the same WSV, while *b* can not be\n"
          "the same WSV as any of the the other vector.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("VectorCrop"),
      DESCRIPTION(
          "Keeps only values of a vector inside the specified range.\n"
          "\n"
          "All values outside the range [min_value,max-value] are removed.\n"
          "Note the default values, that basically should act as -+Inf.\n"
          "\n"
          "The result can either be stored in the same or another vector.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Cropped vector"),
      IN(),
      GIN("in", "min_value", "max_value"),
      GIN_TYPE("Vector", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "-99e99", "99e99"),
      GIN_DESC("Original vector",
               "Minimum value to keep",
               "Maximum value to keep")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorExtractFromMatrix"),
      DESCRIPTION(
          "Extracts a Vector from a Matrix.\n"
          "\n"
          "Copies row or column with given Index from input Matrix variable\n"
          "to create output Vector.\n"),
      AUTHORS("Patrick Eriksson, Oliver Lemke, Stefan Buehler"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Extracted vector."),
      IN(),
      GIN("in", "i", "direction"),
      GIN_TYPE("Matrix", "Index", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Input matrix.",
               "Index of row or column.",
               "Direction. \"row\" or \"column\".")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorFlip"),
      DESCRIPTION(
          "Flips a vector.\n"
          "\n"
          "The output is the input vector in reversed order. The result can\n"
          "either be stored in the same or another vector.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output vector."),
      IN(),
      GIN("in"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Input vector.")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorInsertGridPoints"),
      DESCRIPTION(
          "Insert some additional points into a grid.\n"
          "\n"
          "This method can for example be used to add line center frequencies to\n"
          "a regular frequency grid. If the original grid is [1,2,3], and the\n"
          "additional points are [2.2,2.4], the result will be [1,2,2.2,2.4,3].\n"
          "\n"
          "It is assumed that the original grid is sorted, otherwise a runtime\n"
          "error is thrown. The vector with the points to insert does not have to\n"
          "be sorted. If some of the input points are already in the grid, these\n"
          "points are not inserted again. New points outside the original grid are\n"
          "appended at the appropriate end. Input vector and output vector can be\n"
          "the same.\n"
          "\n"
          "Generic output:\n"
          "  Vector : The new grid vector.\n"
          "\n"
          "Generic input:\n"
          "  Vector : The original grid vector.\n"
          "  Vector : The points to insert.\n"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("The new grid vector"),
      IN(),
      GIN("in", "points"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The original grid vector", "The points to insert")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorLinSpace"),
      DESCRIPTION(
          "Initializes a vector with linear spacing.\n"
          "\n"
          "The first element equals always the start value, and the spacing\n"
          "equals always the step value, but the last value can deviate from\n"
          "the stop value. *step* can be both positive and negative.\n"
          "\n"
          "The created vector is [start, start+step, start+2*step, ...]\n "),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output vector."),
      IN(),
      GIN("start", "stop", "step"),
      GIN_TYPE("Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Start value.",
               "Maximum/minimum value of the end value",
               "Spacing of the vector.")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorLogSpace"),
      DESCRIPTION(
          "Initializes a vector with logarithmic spacing.\n"
          "\n"
          "The first element equals always the start value, and the spacing\n"
          "equals always the step value, but note that the last value can \n"
          "deviate from the stop value. The keyword step can be both positive\n"
          "and negative.\n"
          "\n"
          "Note, that although start has to be given in direct coordinates,\n"
          "step has to be given in log coordinates.\n"
          "\n"
          "Explicitly, the vector is:\n"
          " exp([ln(start), ln(start)+step, ln(start)+2*step, ...])\n"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("start", "stop", "step"),
      GIN_TYPE("Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("The start value. (Direct coordinates!)",
               "The maximum value of the end value. (Direct coordinates!)",
               "The spacing of the vector. (Log coordinates!)")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorMatrixMultiply"),
      DESCRIPTION(
          "Multiply a Vector with a Matrix and store the result in another\n"
          "Vector.\n"
          "\n"
          "This just computes the normal Matrix-Vector product, y=M*x. It is ok\n"
          "if input and output Vector are the same. This function is handy for\n"
          "multiplying the H Matrix to spectra.\n"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("The result of the multiplication (dimension m)."),
      IN(),
      GIN("m", "v"),
      GIN_TYPE("Matrix", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("The Matrix to multiply (dimension mxn).",
               "The original Vector (dimension n).")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorNLinSpace"),
      DESCRIPTION(
          "Creates a vector with length *nelem*, equally spaced between the\n"
          "given end values.\n"
          "\n"
          "The length (*nelem*) must be larger than 1.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Variable to initialize."),
      IN("nelem"),
      GIN("start", "stop"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Start value.", "End value.")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorNLogSpace"),
      DESCRIPTION(
          "Creates a vector with length *nelem*, equally logarithmically\n"
          "spaced between the given end values.\n"
          "\n"
          "The length (*nelem*) must be larger than 1.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Variable to initialize."),
      IN("nelem"),
      GIN("start", "stop"),
      GIN_TYPE("Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Start value.", "End value.")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorReshapeMatrix"),
      DESCRIPTION(
          "Converts a Matrix to a Vector.\n"
          "\n"
          "The matrix is reshaped into a vector. That is, all elements of the matrix\n"
          "are kept. The elements can be extracted both in column (default) and row\n"
          "order. The ouput vector has the same length for both options.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Created vector."),
      IN(),
      GIN("in", "direction"),
      GIN_TYPE("Matrix", "String"),
      GIN_DEFAULT(NODEF, "column"),
      GIN_DESC("Input matrix.", "Direction. \"row\" or \"column\".")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorScale"),
      DESCRIPTION(
          "Scales all elements of a vector with the same value.\n"
          "\n"
          "The result can either be stored in the same or another vector.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Output vector."),
      IN(),
      GIN("in", "value"),
      GIN_TYPE("Vector", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Input vector.", "Scaling value.")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorSetConstant"),
      DESCRIPTION(
          "Creates a vector and sets all elements to the specified value.\n"
          "\n"
          "The vector length is determined by *nelem*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Variable to initialize."),
      IN("nelem"),
      GIN("value"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Vector value.")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorSet"),
      DESCRIPTION(
          "Create a vector from the given list of numbers.\n"
          "\n"
          "   VectorSet(p_grid, [1000, 100, 10] )\n"
          "   Will create a p_grid vector with these three elements.\n"),
      AUTHORS("Stefan Buehler"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Variable to initialize."),
      IN(),
      GIN("value"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("The vector elements."),
      SETMETHOD(true)));

  md_data_raw.push_back(MdRecord(
      NAME("VectorSubtractVector"),
      DESCRIPTION(
          "Element-wise subtraction of two vectors.\n"
          "\n"
          "The method calculates c = a - b.\n"
          "\n"
          "The variable *b* is allowed to have length 1, for any length of\n"
          "*a*. This single value in *b* is then added to every element of *a*.\n"
          "\n"
          "The vectors *a* and *c* can be the same WSV, while *b* can not be\n"
          "the same WSV as any of the the other vector.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("VectorVectorMultiply"),
      DESCRIPTION(
          "Multiply a Vector with another Vector and store result in a third one.\n"
          "\n"
          "This is an element-wise multiplication. It is ok if output Vector is\n"
          "the same as any of the input Vectors.\n"),
      AUTHORS("Jana Mendrok"),
      OUT(),
      GOUT("out"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("The result of the multiplication (dimension m)."),
      IN(),
      GIN("v1", "v2"),
      GIN_TYPE("Vector", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Original Vector #1 (dimension m).",
               "Original Vector #2 (dimension m).")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorZtanToZaRefr1D"),
      DESCRIPTION(
          "Converts a set of true tangent altitudes to zenith angles.\n"
          "\n"
          "The tangent altitudes are given to the function as a vector, which\n"
          "are converted to a generic vector of zenith angles. The position of\n"
          "the sensor is given by the WSV *sensor_pos*. The function works\n"
          "only for 1D. The zenith angles are always set to be positive.\n"),
      AUTHORS("Patrick Eriksson", "Mattias Ekstrom"),
      OUT(),
      GOUT("v_za"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Vector with zenith angles."),
      IN("refr_index_air_agenda",
         "sensor_pos",
         "p_grid",
         "t_field",
         "z_field",
         "vmr_field",
         "refellipsoid",
         "atmosphere_dim",
         "f_grid"),
      GIN("v_ztan"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Vector with tangent altitudes.")));

  md_data_raw.push_back(MdRecord(
      NAME("VectorZtanToZa1D"),
      DESCRIPTION(
          "Converts a set of geometrical tangent altitudes to zenith angles.\n"
          "\n"
          "The tangent altitudes are given to the function as a vector, which\n"
          "are converted to a generic vector of zenith angles. The position of\n"
          "the sensor is given by the WSV *sensor_pos*. The function works\n"
          "only for 1D. The zenith angles are always set to be positive.\n"),
      AUTHORS("Patrick Eriksson", "Mattias Ekstrom"),
      OUT(),
      GOUT("v_za"),
      GOUT_TYPE("Vector"),
      GOUT_DESC("Vector with zenith angles."),
      IN("sensor_pos", "refellipsoid", "atmosphere_dim"),
      GIN("v_ztan"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Vector with tangent altitudes.")));

  md_data_raw.push_back(MdRecord(
      NAME("verbosityInit"),
      DESCRIPTION(
          "Initializes the verbosity levels.\n"
          "\n"
          "Sets verbosity to defaults or the levels specified by -r on the command line.\n"),
      AUTHORS("Oliver Lemke"),
      OUT("verbosity"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("verbosity"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("verbositySet"),
      DESCRIPTION(
          "Sets the verbosity levels.\n"
          "\n"
          "Sets the reporting level for agenda calls, screen and file.\n"
          "All reporting levels can reach from 0 (only error messages)\n"
          "to 3 (everything). The agenda setting applies in addition\n"
          "to both screen and file output.\n"),
      AUTHORS("Oliver Lemke"),
      OUT("verbosity"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("verbosity"),
      GIN("agenda", "screen", "file"),
      GIN_TYPE("Index", "Index", "Index"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("Agenda verbosity level",
               "Screen verbosity level",
               "Report file verbosity level")));

  md_data_raw.push_back(
      MdRecord(NAME("verbositySetAgenda"),
               DESCRIPTION("Sets the verbosity level for agenda output.\n"
                           "\n"
                           "See *verbositySet*\n"),
               AUTHORS("Oliver Lemke"),
               OUT("verbosity"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("verbosity"),
               GIN("level"),
               GIN_TYPE("Index"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Agenda verbosity level")));

  md_data_raw.push_back(
      MdRecord(NAME("verbositySetFile"),
               DESCRIPTION("Sets the verbosity level for report file output.\n"
                           "\n"
                           "See *verbositySet*\n"),
               AUTHORS("Oliver Lemke"),
               OUT("verbosity"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("verbosity"),
               GIN("level"),
               GIN_TYPE("Index"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Report file verbosity level")));

  md_data_raw.push_back(
      MdRecord(NAME("verbositySetScreen"),
               DESCRIPTION("Sets the verbosity level for screen output.\n"
                           "\n"
                           "See *verbositySet*\n"),
               AUTHORS("Oliver Lemke"),
               OUT("verbosity"),
               GOUT(),
               GOUT_TYPE(),
               GOUT_DESC(),
               IN("verbosity"),
               GIN("level"),
               GIN_TYPE("Index"),
               GIN_DEFAULT(NODEF),
               GIN_DESC("Screen verbosity level")));

  md_data_raw.push_back(MdRecord(
      NAME("vmr_fieldClip"),
      DESCRIPTION(
          "Clipping of *vmr_field*.\n"
          "\n"
          "The method allows you to apply hard limits the values of *vmr_field*.\n"
          "All values, of the species selected, below *limit_low*, are simply\n"
          "set to *limit_low*. And the same is performed with respect to\n"
          "*limit_high*. That is, the data in x for the retrieval quantity are\n"
          "forced to be inside the range [limit_low,limit_high].\n"
          "\n"
          "Setting species=\"ALL\", is a shortcut for applying the limits on all\n"
          "species.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("vmr_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("vmr_field", "abs_species"),
      GIN("species", "limit_low", "limit_high"),
      GIN_TYPE("String", "Numeric", "Numeric"),
      GIN_DEFAULT(NODEF, "-Inf", "Inf"),
      GIN_DESC("Name of species to consider, or \"ALL\".",
               "Lower limit for clipping.",
               "Upper limit for clipping.")));

  md_data_raw.push_back(MdRecord(
      NAME("vmr_fieldPerturb"),
      DESCRIPTION(
          "Adds a perturbation to *vmr_field*.\n"
          "\n"
          "Works as *AtmFieldPerturb* but acts on *vmr_field*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("vmr_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("vmr_field","atmosphere_dim","p_grid","lat_grid","lon_grid","abs_species"),
      GIN("species","p_ret_grid","lat_ret_grid","lon_ret_grid",
          "pert_index","pert_size","pert_mode"),
      GIN_TYPE("String","Vector","Vector","Vector",
               "Index","Numeric","String"),
      GIN_DEFAULT(NODEF,NODEF,NODEF,NODEF,NODEF,NODEF,"absolute"),
      GIN_DESC("Name of species to perturb.",
               "Pressure retrieval grid.",
               "Latitude retrieval grid.",
               "Longitude retrieval grid.",
               "Index of position where the perturbation shall be performed.",
               "Size of perturbation.",
               "Type of perturbation, ""ansolute"" or ""relative"".")));

  md_data_raw.push_back(MdRecord(
      NAME("vmr_fieldPerturbAtmGrids"),
      DESCRIPTION(
          "Adds a perturbation to *vmr_field*.\n"
          "\n"
          "Works as *AtmFieldPerturbAtmGrids* but acts on *vmr_field*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("vmr_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("vmr_field","atmosphere_dim","p_grid","lat_grid","lon_grid","abs_species"),
      GIN("species","pert_index","pert_size","pert_mode"),
      GIN_TYPE("String","Index","Numeric","String"),
      GIN_DEFAULT(NODEF,NODEF,NODEF,"absolute"),
      GIN_DESC("Name of species to perturb.",
               "Index of position where the perturbation shall be performed.",
               "Size of perturbation.",
               "Type of perturbation, ""ansolute"" or ""relative"".")));

  md_data_raw.push_back(MdRecord(
      NAME("vmr_fieldSetAllConstant"),
      DESCRIPTION(
          "Sets the VMR of all species to a select constant value.\n"
          "\n"
          "The *vmr_field* WSM must have a correct size before calling this method.\n"
          "The length of vmr_values and of abs_species must match.\n"),
      AUTHORS("Richard Larsson"),
      OUT("vmr_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("vmr_field", "abs_species"),
      GIN("vmr_values"),
      GIN_TYPE("Vector"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("VMR value to apply for each abs_species.")));

  md_data_raw.push_back(MdRecord(
      NAME("water_p_eq_fieldMK05"),
      DESCRIPTION(
          "Calculates *water_p_eq_field* according to Murphy and Koop, 2005.\n"
          "\n"
          "The saturation pressure is set to the one with respect to water at\n"
          "temperatures >= 0C, and to the one with respect to ice for <0C.\n"
          "\n"
          "The saturation pressure with respect to liquid and ice water is\n"
          "calculated according to Eq. 10 and 7, respectively, of:\n"
          "Murphy, D. M., & Koop, T. (2005). Review of the vapour pressures of\n"
          "ice and supercooled water for atmospheric applications. Quarterly\n"
          "Journal of the Royal Meteorological Society, 131(608), 1539-1565.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("water_p_eq_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("t_field"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("vmr_fieldSetConstant"),
      DESCRIPTION(
          "Sets the VMR of a species to a constant value.\n"
          "\n"
          "The *vmr_field* WSM must have a correct size before calling this method.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("vmr_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("vmr_field", "abs_species"),
      GIN("species", "vmr_value"),
      GIN_TYPE("String", "Numeric"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Species to set.",
               "VMR value to apply for the selected species.")));

  md_data_raw.push_back(MdRecord(
      NAME("Wigner6Init"),
      DESCRIPTION("Initialize the wigner 3 and 6 tables\n"
                  "\n"
                  "The default values take about 1 Gb memory.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("Wigner3Init"),
      DESCRIPTION("Initialize the wigner 3 tables\n"
                  "\n"
                  "The default values take about 400 Mb memory.\n"),
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
      MdRecord(NAME("Wigner6Unload"),
               DESCRIPTION("Unloads the wigner 3 and 6 tables\n"),
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

  md_data_raw.push_back(MdRecord(NAME("Wigner3Unload"),
                                 DESCRIPTION("Unloads the wigner 3 tables\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("WignerFastInfoPrint"),
      DESCRIPTION(
          "Prints the fast wigner table information if compiled with this option\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("WindFieldsCalc"),
      DESCRIPTION(
          "Interpolation of raw wind fields to calculation grids.\n"
          "Heritage from *AtmFieldsCalc*\n"
          "\n"
          "Internally, *WindFieldsCalc* applies *GriddedFieldPRegrid* and\n"
          "*GriddedFieldLatLonRegrid*. Generally, 'half-grid-step' extrapolation\n"
          "is allowed and applied.\n"),
      AUTHORS("Richard Larsson"),
      OUT("wind_u_field", "wind_v_field", "wind_w_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("p_grid",
         "lat_grid",
         "lon_grid",
         "wind_u_field_raw",
         "wind_v_field_raw",
         "wind_w_field_raw",
         "atmosphere_dim"),
      GIN("interp_order"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Interpolation order (1=linear interpolation).")));

  md_data_raw.push_back(MdRecord(
      NAME("WindFieldsCalcExpand1D"),
      DESCRIPTION(
          "Interpolation of 1D raw atmospheric fields to create 2D or 3D\n"
          "homogeneous wind fields.  Derived from *AtmFieldsCalcExpand1D*\n"
          "\n"
          "The method works as *WindFieldsCalc*, but accepts only raw 1D\n"
          "wind fields. The raw data is interpolated to *p_grid* and\n"
          "the obtained values are applied for all latitudes, and also\n"
          "longitudes for 3D, to create a homogeneous atmosphere.\n"),
      AUTHORS("Richard Larsson"),
      OUT("wind_u_field", "wind_v_field", "wind_w_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("p_grid",
         "lat_grid",
         "lon_grid",
         "wind_u_field_raw",
         "wind_v_field_raw",
         "wind_w_field_raw",
         "atmosphere_dim"),
      GIN("interp_order"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("1"),
      GIN_DESC("Interpolation order (1=linear interpolation).")));

  md_data_raw.push_back(MdRecord(
      NAME("WindRawRead"),
      DESCRIPTION(
          "Reads wind field data from a scenario.\n"
          "\n"
          "A full set of field components is read (NOTE: fails if scenario\n"
          "only contains selected field components). The files can be\n"
          "anywhere, but must all be in the same directory specified by\n"
          "'basename'. Naming convention for the field component files is\n"
          "basename.wind_u.xml for the u-component, v- and w-components\n"
          "accordingly.\n"),
      AUTHORS("Richard Larsson"),
      OUT("wind_u_field_raw", "wind_v_field_raw", "wind_w_field_raw"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("basename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Name of scenario, probably including the full path. For "
               "example: \"/data/wind_field\"")));

  md_data_raw.push_back(MdRecord(
      NAME("wind_u_fieldIncludePlanetRotation"),
      DESCRIPTION(
          "Maps the planet's rotation to an imaginary wind.\n"
          "\n"
          "This method is of relevance if the observation platform is not\n"
          "following the planet's rotation, and Doppler effects must be\n"
          "considered. Examples include full disk observations from another\n"
          "planet or a satellite not in orbit of the observed planet.\n"
          "\n"
          "The rotation of the planet is not causing any Doppler shift for\n"
          "1D and 2D simulations, and the method can only be used for 3D.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("wind_u_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("wind_u_field",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "refellipsoid",
         "z_field",
         "planet_rotation_period"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("WMRFSelectChannels"),
      DESCRIPTION(
          "Select some channels for WMRF calculation.\n"
          "\n"
          "The HIRS fast setup consists of a precalculated frequency grid\n"
          "covering all HIRS channels, and associated weights for each channel,\n"
          "stored in a weight matrix. (A *sensor_response* matrix.)\n"
          "\n"
          "If not all channels are requested for\n"
          "simulation, then this method can be used to remove the unwanted\n"
          "channels. It changes a number of variables in consistent fashion:\n"
          "\n"
          "- Unwanted channels are removed from f_backend. \n"
          "- Unwanted channels are removed from wmrf_weights.\n"
          "- Unnecessary frequencies are removed from f_grid.\n"
          "- Unnecessary frequencies are removed from wmrf_weights.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("WriteMolTau"),
      DESCRIPTION(
          "Writes a 'molecular_tau_file' as required for libRadtran.\n"
          "\n"
          "The libRadtran (www.libradtran.org) radiative transfer package is a \n"
          "comprehensive package for various applications, it can be used to \n"
          "compute radiances, irradiances, actinic fluxes, ... for the solar \n"
          "and the thermal spectral ranges. Absorption is usually treated using \n"
          "k-distributions or other parameterizations. For calculations with high \n"
          "spectral resolution it requires absorption coefficients from an external \n"
          "line-by-line model. Using this method, arts generates a file that can be \n"
          "used by libRadtran (option molecular_tau_file)."
          "\n"),
      AUTHORS("Claudia Emde"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("f_grid", "z_field", "propmat_clearsky_field", "atmosphere_dim"),
      GIN("filename"),
      GIN_TYPE("String"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Name of the *molecular_tau_file*.")));

  md_data_raw.push_back(MdRecord(
      NAME("WriteNetCDF"),
      DESCRIPTION("Writes a workspace variable to a NetCDF file.\n"
                  "\n"
                  "This method can write variables of limited groups.\n"
                  "\n"
                  "If the filename is omitted, the variable is written\n"
                  "to <basename>.<variable_name>.nc.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN(),
      GIN("in", "filename"),
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

  md_data_raw.push_back(MdRecord(
      NAME("WriteNetCDFIndexed"),
      DESCRIPTION("As *WriteNetCDF*, but creates indexed file names.\n"
                  "\n"
                  "This method can write variables of any group.\n"
                  "\n"
                  "If the filename is omitted, the variable is written\n"
                  "to <basename>.<variable_name>.nc.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("file_index"),
      GIN("in", "filename"),
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

  md_data_raw.push_back(MdRecord(
      NAME("WriteXML"),
      DESCRIPTION("Writes a workspace variable to an XML file.\n"
                  "\n"
                  "This method can write variables of any group.\n"
                  "\n"
                  "If the filename is omitted, the variable is written\n"
                  "to <basename>.<variable_name>.xml.\n"
                  "If no_clobber is set to 1, an increasing number will be\n"
                  "appended to the filename if the file already exists.\n"),
      AUTHORS("Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("output_file_format"),
      GIN("in", "filename", "no_clobber"),
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

  md_data_raw.push_back(MdRecord(
      NAME("WriteXMLIndexed"),
      DESCRIPTION("As *WriteXML*, but creates indexed file names.\n"
                  "\n"
                  "The variable is written to a file with name:\n"
                  "   <filename>.<file_index>.xml.\n"
                  "where <file_index> is the value of *file_index*.\n"
                  "\n"
                  "This means that *filename* shall here not include the .xml\n"
                  "extension. Omitting filename works as for *WriteXML*.\n"),
      AUTHORS("Patrick Eriksson, Oliver Lemke"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("output_file_format", "file_index"),
      GIN("in", "filename", "digits"),
      GIN_TYPE("Any", "String", "Index"),
      GIN_DEFAULT(NODEF, "", "0"),
      GIN_DESC(
          "Workspace variable to be saved.",
          "File name. See above.",
          "Equalize the widths of all numbers by padding with zeros as necessary.\n"
          "0 means no padding (default)."),
      SETMETHOD(false),
      AGENDAMETHOD(false),
      USES_TEMPLATES(true),
      PASSWORKSPACE(false),
      PASSWSVNAMES(true)));

  md_data_raw.push_back(MdRecord(
      NAME("xaStandard"),
      DESCRIPTION(
          "Standard function for creating *xa*.\n"
          "\n"
          "The method creates *xa* based on *jacobian_quantities* and the various\n"
          "atmospheric fields. In the case of scattering species, the data are\n"
          "taken from *particle_bulkprop_field*. The following retrieval quantities\n"
          "are handled:\n"
          "   Temperature\n"
          "   Absorption species\n"
          "   Scattering species\n"
          "   Pointing\n"
          "   Polynomial baseline fit\n"
          "   Sinusoidal baseline fit\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("xa"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("jacobian_quantities",
         "atmfields_checked",
         "atmgeom_checked",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "t_field",
         "vmr_field",
         "abs_species",
         "cloudbox_on",
         "cloudbox_checked",
         "particle_bulkprop_field",
         "particle_bulkprop_names",
         "wind_u_field",
         "wind_v_field",
         "wind_w_field",
         "mag_u_field",
         "mag_v_field",
         "mag_w_field",
         "surface_props_data",
         "surface_props_names",
         "water_p_eq_agenda"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("xClip"),
      DESCRIPTION(
          "Clipping of the state vector.\n"
          "\n"
          "The method allows you to apply hard limits the values of a\n"
          "retrieval quantity. The retrieval quantity is specified by\n"
          "*ijq*. All values of the quantity below *limit_low*, are simply\n"
          "set to *limit_low*. And the same is performed with respect to\n"
          "*limit_high*. That is, the data in x for the retrieval quantity\n"
          "are forced to be inside the range [limit_low,limit_high].\n"
          "\n"
          "Setting ijq=-1, is a shortcut for applying the limits on all\n"
          "retrieval quantities.\n"
          "\n"
          "Notice that limits must be specified in the unit used in *x*.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("x2artsAtmAndSurf"),
      DESCRIPTION(
          "Maps *x* to atmospheric and surface variables.\n"
          "\n"
          "Maps OEM's state vector, *x*, to the matching ARTS variables. This\n"
          "method handles atmospheric and surface variables. If you retrieve\n"
          "other variables, make sure that you also call *x2artsSensor* and/or\n"
          "*x2artsSpectroscopy*.\n"
          "\n"
          "The following retrieval quantities are handled by this method:\n"
          "   Temperature\n"
          "   Absorption species\n"
          "   Scattering species\n"
          "   Winds\n"
          "   Surface variables\n"
          "\n"
          "Should only be used inside *inversion_iterate_agenda*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("vmr_field",
          "t_field",
          "particle_bulkprop_field",
          "wind_u_field",
          "wind_v_field",
          "wind_w_field",
          "mag_u_field",
          "mag_v_field",
          "mag_w_field",
          "surface_props_data"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("vmr_field",
         "t_field",
         "particle_bulkprop_field",
         "wind_u_field",
         "wind_v_field",
         "wind_w_field",
         "mag_u_field",
         "mag_v_field",
         "mag_w_field",
         "surface_props_data",
         "jacobian_quantities",
         "x",
         "atmfields_checked",
         "atmgeom_checked",
         "atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
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

  md_data_raw.push_back(MdRecord(
      NAME("x2artsSensor"),
      DESCRIPTION(
          "Maps *x* to sensor variables.\n"
          "\n"
          "Maps OEM's state vector, *x*, to the matching ARTS variables. This\n"
          "method handles variables associated with the sensor. If you retrieve\n"
          "other variables, make sure that you also call *x2artsAtmAndSurf*\n"
          " and/or *x2artsSpectroscopy*.\n"
          "\n"
          "The following retrieval quantities are handled by this method:\n"
          "   Pointing\n"
          "   Frequency shift and stretch\n"
          "   Baseline fits\n"
          "\n"
          "Should only be used inside *inversion_iterate_agenda*.\n"
          "\n"
          "Elements in *x* representing pointing corrections are mapped to\n"
          "*sensor_los*. Elements representing frequency corrections are mapped\n"
          "to *f_backend*. Baseline variables are mapped to *y_baseline*.\n"
          "\n"
          "The sensor response is recalculated if there is any non-zero frequency\n"
          "correction.\n"),
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
          "mblock_dlos_grid"),
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
         "mblock_dlos_grid",
         "jacobian_quantities",
         "x",
         "sensor_response_agenda",
         "sensor_checked",
         "sensor_time"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("x2artsSpectroscopy"),
      DESCRIPTION("Just defined to indicate a future extensiom.\n"
                  "\n"
                  "Don't call the method, it will just generate an error.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("yApplySensorPol"),
      DESCRIPTION(
          "Extraction of arbitrary linear polarisation.\n"
          "\n"
          "This method shall be called after *yCalc* and then applies *sensor_pol*\n"
          "on the outout of *yCalc*. See *sensor_pol* for definition of the\n"
          "polarisation responses. THe *sensor_response* give to *yCalc* can not\n"
          "contain any polarisation response, it must maintain original Stokes\n"
          "elelemnts. The value of *stokes_dim* muist be >= 3.\n"
          "\n"
          "The values in *sensor_pol* are applied on *y*, and *jacobian* if relevant.\n"
          "*y_pol* is set following the values in *sensor_pol* but is rounded to\n"
          "an integer value. Remaining data associated with *y* (e.g. y_pos) are\n"
          "set to the value matching the first Stokes element.\n"),
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
         "stokes_dim",
         "jacobian_do",
         "sensor_pos",
         "sensor_pol"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("yApplyUnit"),
      DESCRIPTION(
          "Conversion of *y* to other spectral units.\n"
          "\n"
          "Any conversion to brightness temperature is normally made inside\n"
          "*yCalc*. This method makes it possible to also make this conversion\n"
          "after *yCalc*, but with restrictions for *jacobian* and with.\n"
          "respect to the n2-law of radiance.\n"
          "\n"
          "The conversion made inside *iyEmissionStandard* is mimiced\n"
          "and see that method for constraints and selection of output units.\n"
          "This with the restriction that the n2-law can be ignored. The later\n"
          "is the case if the sensor is placed in space, or if the refractive\n"
          "only devaites slightly from unity.\n"
          "\n"
          "The method handles *y* and *jacobian* in parallel, where\n"
          "the last variable is only considered if it is set. The\n"
          "input data must be in original radiance units. A completely\n"
          "stringent check of this can not be performed.\n"
          "\n"
          "The method can not be used with jacobian quantities that are not\n"
          "obtained through radiative transfer calculations. One example on\n"
          "quantity that can not be handled is *jacobianAddPolyfit*. There\n"
          "are no automatic checks warning for incorrect usage!\n"
          "\n"
          "If you are using this method, *iy_unit* should be set to \"1\" when\n"
          "calling *yCalc*, and be changed before calling this method.\n"
          "\n"
          "Conversion of *y_aux* is not supported.\n"),
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

  md_data_raw.push_back(MdRecord(
      NAME("ybatchCalc"),
      DESCRIPTION(
          "Performs batch calculations for the measurement vector y.\n"
          "\n"
          "We perform *ybatch_n* jobs, starting at index *ybatch_start*. (Zero\n"
          "based indexing, as usual.) The output array *ybatch* will have\n"
          "ybatch_n elements. Indices in the output array start\n"
          "with zero, independent of *ybatch_start*.\n"
          "\n"
          "The method performs the following:\n"
          "   1. Sets *ybatch_index* = *ybatch_start*.\n"
          "   2. Performs a-d until\n"
          "      *ybatch_index* = *ybatch_start* + *ybatch_n*.\n"
          "        a. Executes *ybatch_calc_agenda*.\n"
          "        b. If *ybatch_index* = *ybatch_start*, resizes *ybatch*\n"
          "           based on *ybatch_n* and length of *y*.\n"
          "        c. Copies *y* to *ybatch_index* - *ybatch_start*\n"
          "           of *ybatch*.\n"
          "        d. Adds 1 to *ybatch_index*.\n"
          "\n"
          "Beside the *ybatch_calc_agenda*, the WSVs *ybatch_start*\n"
          "and *ybatch_n* must be set before calling this method.\n"
          "Further, *ybatch_calc_agenda* is expected to produce a\n"
          "spectrum and should accordingly include a call of *yCalc*\n"
          "(or asimilar method).\n"
          "\n"
          "The input variable *ybatch_start* is set to a default of zero in\n"
          "*general.arts*.\n"
          "\n"
          "An agenda that calculates spectra for different temperature profiles\n"
          "could look like this:\n"
          "\n"
          "   AgendaSet(ybatch_calc_agenda){\n"
          "      Extract(t_field,tensor4_1,ybatch_index)\n"
          "      yCalc\n"
          "   }\n"
          "\n"
          "Jacobians are also collected, and stored in output variable *ybatch_jacobians*. \n"
          "(This will be empty if yCalc produces empty Jacobians.)\n"
          "\n"
          "See the user guide for further practical examples.\n"),
      AUTHORS("Stefan Buehler"),
      OUT("ybatch", "ybatch_aux", "ybatch_jacobians"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("ybatch_start", "ybatch_n", "ybatch_calc_agenda"),
      GIN("robust"),
      GIN_TYPE("Index"),
      GIN_DEFAULT("0"),
      GIN_DESC("A flag with value 1 or 0. If set to one, the batch\n"
               "calculation will continue, even if individual jobs fail. In\n"
               "that case, a warning message is written to screen and file\n"
               "(out1 output stream), and the *y* Vector entry for the\n"
               "failed job in *ybatch* is left empty.")));

  md_data_raw.push_back(MdRecord(
      NAME("ybatchMetProfiles"),
      DESCRIPTION(
          "This method is used for simulating ARTS for metoffice model fields"
          "\n"
          "This method reads in *met_amsu_data* which contains the\n"
          "lat-lon of the metoffice profile files as a Matrix. It then\n"
          "loops over the number of profiles and corresponding to each\n"
          "longitude create the appropriate profile basename. Then,\n"
          "corresponding to each basename we have temperature field, altitude\n"
          "field, humidity field, and particle number density field. The\n"
          "temperature field and altitude field are stored in the same dimensions\n"
          "as *t_field_raw* and *z_field_raw*. The oxygen and nitrogen VMRs are\n"
          "set to constant values of 0.209 and 0.782, respectively and are used\n"
          "along with humidity field to generate *vmr_field_raw*. \n"
          "\n"
          "The three fields *t_field_raw*, *z_field_raw*, and *vmr_field_raw* are\n"
          "given as input to *met_profile_calc_agenda* which is called in this\n"
          "method. See documentation of WSM *met_profile_calc_agenda* for more\n"
          "information on this agenda. \n"
          "\n"
          "The method also converts satellite zenith angle to appropriate\n"
          "*sensor_los*. It also sets the *p_grid* and *cloudbox_limits*\n"
          "from the profiles inside the function\n"),
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
         "refellipsoid",
         "lat_grid",
         "lon_grid",
         "atmosphere_dim",
         "scat_data"),
      GIN("nelem_p_grid", "met_profile_path", "met_profile_pnd_path"),
      GIN_TYPE("Index", "String", "String"),
      GIN_DEFAULT(NODEF, NODEF, NODEF),
      GIN_DESC("FIXME DOC", "FIXME DOC", "FIXME DOC")));

  md_data_raw.push_back(MdRecord(
      NAME("ybatchMetProfilesClear"),
      DESCRIPTION(
          "This method is used for simulating ARTS for metoffice model fields\n"
          "for clear sky conditions.\n"
          "\n"
          "This method reads in *met_amsu_data* which contains the\n"
          "lat-lon of the metoffice profile files as a Matrix. It then\n"
          "loops over the number of profiles and corresponding to each\n"
          "longitude create the appropriate profile basename. Then,\n"
          "Corresponding to each basename we have temperature field, altitude\n"
          "field, humidity field, and particle number density field. The\n"
          "temperature field and altitude field are stored in the same dimensions\n"
          "as *t_field_raw* and *z_field_raw*. The oxygen and nitrogen VMRs are\n"
          "set to constant values of 0.209 and 0.782, respectively and are used\n"
          "along with humidity field to generate *vmr_field_raw*. \n"
          "\n"
          "The three fields *t_field_raw*, *z_field_raw*, and *vmr_field_raw* are\n"
          "given as input to *met_profile_calc_agenda* which is called in this\n"
          "method. See documentation of WSM *met_profile_calc_agenda* for more\n"
          "information on this agenda. \n"
          "\n"
          "The method also converts satellite zenith angle to appropriate\n"
          "*sensor_los*. It also sets the *p_grid* and *cloudbox_limits*\n"
          "from the profiles inside the function\n"),
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
         "refellipsoid"),
      GIN("nelem_p_grid", "met_profile_path"),
      GIN_TYPE("Index", "String"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("FIXME DOC", "FIXME DOC")));

  md_data_raw.push_back(MdRecord(
      NAME("yCalc"),
      DESCRIPTION(
          "Calculation of complete measurement vectors (y).\n"
          "\n"
          "The method performs radiative transfer calculations from a sensor\n"
          "perspective. Radiative transfer calculations are performed for\n"
          "monochromatic pencil beams, following *iy_main_agenda* and\n"
          "associated agendas. Obtained radiances are weighted together by\n"
          "*sensor_response*, to include the characteristics of the sensor.\n"
          "The measurement vector obtained can contain anything from a single\n"
          "frequency value to a series of measurement scans (each consisting\n"
          "of a series of spectra), all depending on the settings. Spectra\n"
          "and jacobians are calculated in parallel.\n"
          "\n"
          "The frequency, polarisation etc. for each measurement value is\n"
          "given by *y_f*, *y_pol*, *y_pos* and *y_los*.\n"
          "\n"
          "The content of *y_aux* follows *iy_aux_vars. See the method selected\n"
          "for *iy_main_agenda* for allowed choices.\n"
          "\n"
          "The geo-positions (*y_geo*) are set based on *sensor_reponse*. When\n"
          "an antenna pattern is considered, there are several pencil beams,\n"
          "and thus also several goe-positions, associated with each value of *y*.\n"
          "The geo-position assigned to a value in *y* is the *geo_pos* of the pencil\n"
          "beam related to the highest value in *sensor_response*. This means that\n"
          "*mblock_dlos_grid* must contain the bore-sight direction (0,0), if you\n"
          "want *y_geo* to exactly match the bore-sight direction.\n"
          "\n"
          "The Jacobian provided (*jacobian*) is adopted to selected retrieval\n"
          "units, but no transformations are applied. Transformations are\n"
          "included by calling *jacobianAdjustAndTranform*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("y", "y_f", "y_pol", "y_pos", "y_los", "y_aux", "y_geo", "jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmgeom_checked",
         "atmfields_checked",
         "atmosphere_dim",
         "nlte_field",
         "cloudbox_on",
         "cloudbox_checked",
         "scat_data_checked",
         "sensor_checked",
         "stokes_dim",
         "f_grid",
         "sensor_pos",
         "sensor_los",
         "transmitter_pos",
         "mblock_dlos_grid",
         "sensor_response",
         "sensor_response_f",
         "sensor_response_pol",
         "sensor_response_dlos",
         "iy_unit",
         "iy_main_agenda",
         "geo_pos_agenda",
         "jacobian_agenda",
         "jacobian_do",
         "jacobian_quantities",
         "iy_aux_vars"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("yCalcAppend"),
      DESCRIPTION(
          "Replaces *yCalc* if a measurement shall be appended to an\n"
          "existing one.\n"
          "\n"
          "The method works basically as *yCalc* but appends the results to\n"
          "existing data, instead of creating completely new *y* and its\n"
          "associated variables. This method is required if your measurement\n"
          "consists of data from two instruments using different observation\n"
          "techniques (corresponding to different iyCalc-methods). One such\n"
          "example is if emission and transmission data are combined into a\n"
          "joint retrieval. The method can also be used to get around the\n"
          "constrain that *sensor_response* is required to be the same for\n"
          "all data.\n"
          "\n"
          "The new measurement is simply appended to the input *y*, and the\n"
          "other output variables are treated correspondingly. Data are\n"
          "appended \"blindly\" in *y_aux*. That is, data of different type\n"
          "are appended if *iy_aux_vars* differs between the two measurements,\n"
          "the data are appended strictly following the order. First variable\n"
          "of second measurement is appended to first variable of first\n"
          "measurement, and so on. The number of auxiliary variables can differ\n"
          "between the measurements. Missing data are set to zero.\n"
          "\n"
          "The set of retrieval quantities can differ between the two\n"
          "calculations. If an atmospheric quantity is part of both Jacobians,\n"
          "the same retrieval grids must be used in both cases.\n"
          "The treatment of instrument related Jacobians (baseline fits,\n"
          "pointing ...) follows the *append_instrument_wfs* argument.\n"
          "\n"
          "A difference to *yCalc* is that *jacobian_quantities* is both in-\n"
          "and output variable. The input version shall match the measurement\n"
          "to be calculated, while the output version matches the output *y*,\n"
          "the combined, measurements. A copies of *jacobian_quantities* of the\n"
          "first measurement must be made and shall be provided to the method\n"
          "as *jacobian_quantities_copy*.\n"
          "\n"
          "As for *yCalc* Jacobian transformations are not handled, and the\n"
          "the input Jacobian shall not contain transformations. That is\n"
          "*jacobainAdjustAndTranform* shall be called after this method,\n"
          "when the complete Jacobian is at hand.\n"),
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
         "atmosphere_dim",
         "nlte_field",
         "cloudbox_on",
         "cloudbox_checked",
         "scat_data_checked",
         "sensor_checked",
         "stokes_dim",
         "f_grid",
         "sensor_pos",
         "sensor_los",
         "transmitter_pos",
         "mblock_dlos_grid",
         "sensor_response",
         "sensor_response_f",
         "sensor_response_pol",
         "sensor_response_dlos",
         "iy_unit",
         "iy_main_agenda",
         "geo_pos_agenda",
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

  md_data_raw.push_back(MdRecord(
      NAME("yActive"),
      DESCRIPTION(
          "Replaces *yCalc* for radar/lidar calculations.\n"
          "\n"
          "The output format for *iy* when simulating radars and lidars differs\n"
          "from the standard one, and *yCalc* can not be used for such simulations.\n"
          "This method works largely as *yCalc*, but is tailored to handle the\n"
          "output from *iyActiveSingleScattering*.\n"
          "\n"
          "The method requires additional information about the sensor,\n"
          "regarding its recieving properties. First of all, recieved\n"
          "polarisation states are taken from *instrument_pol_array*. Note\n"
          "that this WSV allows to define several measured polarisations\n"
          "for each transmitted signal. For example, it is possible to\n"
          "simulate transmission of V and measuring backsacttered V and H.\n"
          "\n"
          "Secondly, the range averaging is described by *range_bins*. These\n"
          "bins can either be specified in altitude or two-way travel time.\n"
          "In both case, the edges of the range bins shall be specified.\n"
          "All data (including auxiliary variables) are returned as the\n"
          "average inside the bins. If a bin is totally outside the model\n"
          "atmosphere, NaN is returned.\n"
          "\n"
          "The options for *iy_unit* are:\n"
          " \"1\"   : Backscatter coefficient. Unit is 1/(m*sr). At zero\n"
          "           attenuation, this equals the scattering matrix value for\n"
          "           the backward direction. See further AUG.\n"
          " \"Ze\"  : Equivalent reflectivity. Unit is mm^6/m^3. Conversion\n"
          "           formula is given below.\n"
          " \"dBZe\": 10*log10(Ze/Z0), where Z0 is 1 mm^6/m^3.\n"
          "\n"
          "The conversion from backscatter coefficient to Ze is:\n"
          "   Ze = 1e18 * lambda^4 / (k2 * pi^5) * sum(sigma),\n"
          "where sum(sigma) = 4*pi*b, and b is the backscatter coefficient.\n"
          "\n"
          "The reference dielectric factor can either specified directly by\n"
          "the argument *k2*. For example, to mimic the CloudSat data, *k2*\n"
          "shall be set to 0.75 (citaion needed). If *k2* is set to be \n"
          "negative (which is defualt), k2 is calculated as:\n"
          "   k2 = abs( (n^2-1)/(n^2+2) )^2,\n"
          "where n is the refractive index of liquid water at temperature\n"
          "*ze_tref* and the frequency of the radar, calculated by the MPM93\n"
          "parameterization.\n"
          "\n"
          "A lower limit for dBZe is applied (*dbze_min*). The main reason is to\n"
          "handle the fact dBZe is not defined for Ze=0, and dBZe is set to the\n"
          "clip value when Ze < 10^(dbze_min/10).\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("y", "y_f", "y_pol", "y_pos", "y_los", "y_aux", "y_geo", "jacobian"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmgeom_checked",
         "atmfields_checked",
         "iy_unit",
         "iy_aux_vars",
         "stokes_dim",
         "f_grid",
         "atmosphere_dim",
         "nlte_field",
         "cloudbox_on",
         "cloudbox_checked",
         "sensor_pos",
         "sensor_los",
         "sensor_checked",
         "jacobian_do",
         "jacobian_quantities",
         "iy_main_agenda",
         "geo_pos_agenda",
         "instrument_pol_array",
         "range_bins"),
      GIN("ze_tref", "k2", "dbze_min"),
      GIN_TYPE("Numeric", "Numeric", "Numeric"),
      GIN_DEFAULT("273.15", "-1", "-99"),
      GIN_DESC("Reference temperature for conversion to Ze.",
               "Reference dielectric factor.",
               "Clip value for dBZe.")));

  md_data_raw.push_back(MdRecord(
      NAME("ySimpleSpectrometer"),
      DESCRIPTION(
          "Converts *iy* to *y* assuming a fixed frequency resolution.\n"
          "\n"
          "This is a short-cut, avoiding *yCalc*, that can be used to convert\n"
          "monochromatic pencil beam data to spectra with a fixed resolution.\n"
          "\n"
          "The method mimics a spectrometer with rectangular response\n"
          "functions, all having the same width (*df*). The position of\n"
          "the first spectrometer channel is set to f_grid[0]+df/2.\n"
          "The centre frequency of channels are returned as *y_f*.\n"
          "\n"
          "Auxiliary variables and *jacobian*s are not handled.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("y", "y_f"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("iy", "stokes_dim", "f_grid"),
      GIN("df"),
      GIN_TYPE("Numeric"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Selected frequency resolution.")));

  md_data_raw.push_back(MdRecord(
      NAME("z_fieldFromHSE"),
      DESCRIPTION(
          "Force altitudes to fulfil hydrostatic equilibrium.\n"
          "\n"
          "The method applies hydrostatic equilibrium. A mixture of \"dry\n"
          "air\" and water vapour (if present as *abs_species* tag) is assumed.\n"
          "That is, the air is assumed to be well mixed and its weight, apart\n"
          "from the water vapour, is constant (*molarmass_dry_air*). In\n"
          "addition, the effect of any particles (including liquid and ice\n"
          "particles) is neglected.\n"
          "\n"
          "The output is an update of *z_field*. This variable is expected to\n"
          "contain approximative altitudes when calling the function. The\n"
          "altitude matching *p_hse* is kept constant. Other input altitudes can\n"
          "basically be arbitrary, but good estimates give quicker calculations.\n"
          "\n"
          "The calculations are repeated until the change in altitude is below\n"
          "*z_hse_accuracy*. An iterative process is needed as gravity varies\n"
          "with altitude.\n"
          "\n"
          "For 1D and 2D, the geographical position is taken from *lat_true*\n"
          "and *lon_true*.\n"),
      AUTHORS("Patrick Eriksson"),
      OUT("z_field"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("atmosphere_dim",
         "p_grid",
         "lat_grid",
         "lon_grid",
         "lat_true",
         "lon_true",
         "abs_species",
         "t_field",
         "z_field",
         "vmr_field",
         "refellipsoid",
         "z_surface",
         "atmfields_checked",
         "g0_agenda",
         "molarmass_dry_air",
         "p_hse",
         "z_hse_accuracy"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("zeeman_linerecord_precalcCreateFromLines"),
      DESCRIPTION(
          "Creates a Zeeman ArrayOfArrayOfLineRecord prior to atmospheric looping.\n"
          "\n"
          "Will reset the *zeeman_linerecord_precalc* at every call.\n"),
      AUTHORS("Richard Larsson"),
      OUT("zeeman_linerecord_precalc"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species", "abs_lines_per_species", "wigner_initialized"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("zeeman_linerecord_precalcCreateWithZeroSplitting"),
      DESCRIPTION(
          "Creates a Zeeman ArrayOfArrayOfLineRecord prior to atmospheric looping.\n"
          "\n"
          "Will reset the *zeeman_linerecord_precalc* at every call.\n"
          "\n"
          "The splitting coefficients are set to zero.  Use a modify-function to adjust these\n"
          "for proper calculations\n"),
      AUTHORS("Richard Larsson"),
      OUT("zeeman_linerecord_precalc"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("abs_species", "abs_lines_per_species", "wigner_initialized"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()));

  md_data_raw.push_back(MdRecord(
      NAME("zeeman_linerecord_precalcModifyFromData"),
      DESCRIPTION(
          "Modifies *zeeman_linerecord_precalc* by setting its g value directly.\n"
          "\n"
          "This will set all unmatched energy levels to g=0 for all lines that\n"
          "have at least a single level match\n"),
      AUTHORS("Richard Larsson"),
      OUT("zeeman_linerecord_precalc"),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("zeeman_linerecord_precalc"),
      GIN("keys", "data"),
      GIN_TYPE("ArrayOfQuantumIdentifier", "Vector"),
      GIN_DEFAULT(NODEF, NODEF),
      GIN_DESC("Keys for energy levels in the line array",
               "Matching g-values to set for levels")));

  md_data_raw.push_back(MdRecord(
      NAME("zeeman_linerecord_precalcPrintMissing"),
      DESCRIPTION("Support for *zeeman_linerecord_precalcModifyFromData*\n"
                  "\n"
                  "Prints all lines not having any partial match for *keys*\n"
                  "in *zeeman_linerecord_precalcModifyFromData*\n"),
      AUTHORS("Richard Larsson"),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN("zeeman_linerecord_precalc"),
      GIN("keys"),
      GIN_TYPE("ArrayOfQuantumIdentifier"),
      GIN_DEFAULT(NODEF),
      GIN_DESC("Keys for energy levels in the line array")));

}
