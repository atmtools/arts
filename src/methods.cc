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

#include "arts.h"
#include "make_array.h"
#include "methods.h"
#include "wsv_aux.h"

// Some #defines and typedefs to make the records better readable:
#define NAME(x) x
#define DESCRIPTION(x) x
#define AUTHORS     MakeArray<String>
#define OUT         MakeArray<String>
#define GOUT        MakeArray<String>
#define GOUT_TYPE   MakeArray<String>
#define GOUT_DESC   MakeArray<String>
#define IN          MakeArray<String>
#define GIN         MakeArray<String>
#define GIN_TYPE    MakeArray<String>
#define GIN_DEFAULT MakeArray<String>
#define GIN_DESC    MakeArray<String>
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



void define_md_data_raw()
{
  // The variable md_data is defined in file methods_aux.cc.
  extern Array<MdRecord> md_data_raw;

  // Initialise to zero, just in case:
  md_data_raw.resize(0);

  // String with all array groups
  const String ARRAY_GROUPS = get_array_groups_as_string();
  // String with all groups that also exist as an array type
  const String GROUPS_WITH_ARRAY_TYPE = get_array_groups_as_string(true, true);
  // String with all array types whose element type is also available as a group
  const String ARRAY_GROUPS_WITH_BASETYPE = get_array_groups_as_string(true, false);

  extern const ArrayOfString wsv_group_names;
  
  for (ArrayOfString::const_iterator it = wsv_group_names.begin();
       it != wsv_group_names.end(); it++)
  {
    if (*it != "Any")
    {
      md_data_raw.push_back
      (MdRecord
       (NAME( String(*it + "Create").c_str() ),
        DESCRIPTION
        (
         String("Creates a variable of group " + *it + ".\n"
                "\n"
                "If the variable already exists, it'll be reinitialized.\n").c_str()
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(      "var" ),
        GOUT_TYPE( (*it).c_str() ),
        GOUT_DESC( "Variable to create." ),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        )
       );
    }
  }
  
  /////////////////////////////////////////////////////////////////////////////
  // Let's put in the functions in alphabetical order. This gives a clear rule
  // for where to place a new function and this gives a nicer results when
  // the functions are listed by "arts -m all".
  // No distinction is made between uppercase and lowercase letters. The sign
  // "_" comes after all letters.
  // Patrick Eriksson 2002-05-08
  /////////////////////////////////////////////////////////////////////////////

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AbsInputFromAtmFields" ),
        DESCRIPTION
        (
         "Initialises the WSVs *abs_p*, *abs_t* and *abs_vmrs* from\n"
         "*p_grid, *t_field* and *vmr_field*.\n"
         "\n"
         "This only works for a 1D atmosphere!\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_p", "abs_t", "abs_vmrs" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim", "p_grid", "t_field", "vmr_field" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
  
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AbsInputFromRteScalars" ),
        DESCRIPTION
        (
         "Initialize absorption input WSVs from local atmospheric conditions.\n"
         "\n"
         "The purpose of this method is to allow an explicit line-by-line\n"
         "calculation, e.g., by *abs_coefCalc*, to be put inside the\n"
         "*abs_mat_per_species_agenda*. What the method does is to prepare absorption\n"
         "input parameters (pressure, temperature, VMRs), from the input\n"
         "parameters to *abs_mat_per_species_agenda*.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_p", "abs_t", "abs_vmrs" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "rte_pressure", "rte_temperature", "rte_vmr_list" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_coefCalc" ),
        DESCRIPTION
        (
         "Calculate absorption coefficients.\n"
         "\n"
         "This function calculates both the total absorption (*abs_coef*), and\n"
         "the absorption per species (*abs_coef_per_species*).\n"
         "\n"
         "The method calls four other  methods:\n"
         "\n"
         "1. *abs_xsec_per_speciesInit*:\n"
         "   Initialize *abs_xsec_per_species*\n"
         "\n"
         "2. *abs_xsec_per_speciesAddLines*:\n"
         "   Calculate cross sections per tag group for line spectra.\n"
         "\n"
         "3. *abs_xsec_per_speciesAddConts*:\n"
         "   Calculate cross sections per tag group for continua.\n"
         "\n"
         "4. *abs_coefCalcFromXsec*:\n"
         "   Calculate absorption coefficients from the cross sections by\n"
         "   multiplying each cross section by n*VMR.\n"
         "\n"
         "This is done once for each tag group (output *abs_coef_per_species*),\n"
         "and for the sum of all tag groups (output *abs_coef*).\n"
         ),
        AUTHORS( "Axel von Engeln", "Stefan Buehler" ),
        OUT( "abs_coef"  , "abs_coef_per_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species", "f_grid", "abs_p", "abs_t", "abs_n2", "abs_h2o",
            "abs_vmrs", "abs_lines_per_species", "abs_lineshape",
            "abs_cont_names", "abs_cont_models", 
            "abs_cont_parameters" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_coefCalcFromXsec" ),
        DESCRIPTION
        (
         "Calculate absorption coefficients from cross sections.\n"
         "\n"
         "This calculates both the total absorption and the\n"
         "absorption per species.\n"
         "\n"
         "Cross sections are multiplied by n*VMR.\n"
         ),
        AUTHORS( "Stefan Buehler", "Axel von Engeln" ),
        OUT( "abs_coef", "abs_coef_per_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_xsec_per_species", "abs_vmrs", "abs_p", "abs_t" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_coefCalcSaveMemory" ),
        DESCRIPTION
        (
         "Calculate absorption coefficients, trying to conserve memory.\n"
         "\n"
         "This function calculates only the total absorption (*abs_coef*),\n"
         "NOT the absorption per tag group (*abs_coef_per_species*).\n"
         "\n"
         "This means you cannot use it if you want to calculate Jacobians\n"
         "later.\n"
         "\n"
         "The implementation follows abs_coefCalc.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_coef" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species", "f_grid", "abs_p", "abs_t", "abs_n2", "abs_h2o",
            "abs_vmrs", "abs_lines_per_species", "abs_lineshape",
            "abs_cont_names", "abs_cont_models", 
            "abs_cont_parameters" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_cont_descriptionAppend" ),
        DESCRIPTION
        (
         "Appends the description of a continuum model or a complete absorption\n"
         "model to *abs_cont_names* and *abs_cont_parameters*.\n"
         "\n"
         "See online documentation for *abs_cont_names* for a list of\n"
         "allowed models and for information what parameters they require. See\n"
         "file includes/continua.arts for default parameters for the various models.\n"
         ),
        AUTHORS( "Thomas Kuhn", "Stefan Buehler" ),
        OUT( "abs_cont_names", 
             "abs_cont_models",
             "abs_cont_parameters" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_cont_names", 
            "abs_cont_models",
            "abs_cont_parameters" ),
        GIN( "tagname", "model",  "userparameters" ),
        GIN_TYPE(    "String",  "String", "Vector" ),
        GIN_DEFAULT( NODEF,     NODEF,    NODEF),
        GIN_DESC(
         "The name (species tag) of a continuum model. Must match one\n"
         "of the models implemented in ARTS.\n",
         "A string selecting a particular continuum/full model under this\n"
         "species tag.\n",
         "A Vector containing the required parameters for the selected model.\n"
         "The meaning of the parameters and how many parameters are required\n"
         "depends on the model.\n" )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_cont_descriptionInit" ),
        DESCRIPTION
        (
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
         "not want to use any continua.\n"
         ),
        AUTHORS( "Thomas Kuhn", "Stefan Buehler" ),
        OUT( "abs_cont_names", 
             "abs_cont_models",
             "abs_cont_parameters" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_h2oSet" ),
        DESCRIPTION
        (
         "Sets abs_h2o to the profile of the first tag group containing\n"
         "water.\n" 
         "\n"
         "This is necessary, because for example *abs_coefCalc* requires abs_h2o\n"
         "to contain the water vapour profile(the reason for this is the\n"
         "calculation of oxygen line broadening requires water vapour profile).\n"
         "Then this function can be used to copy the profile of the first tag\n"
         "group of water.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_h2o" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species", "abs_vmrs" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_lineshapeDefine" ),
        DESCRIPTION
        (
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
         "   Voigt_Kuntz3:             Kuntz approximation to the Voigt lineshape,\n"
         "                             accuracy > 2x10^(-3)\n"
         "   Voigt_Kuntz4:             Kuntz approximation to the Voigt lineshape,\n"
         "                             accuracy > 2x10^(-4)\n"
         "   Voigt_Kuntz6:             Kuntz approximation to the Voigt lineshape,\n"
         "                             accuracy > 2x10^(-6)\n"
         "   Voigt_Drayson:            Drayson approximation to the Voigt lineshape\n"
         "   Rosenkranz_Voigt_Drayson: Rosenkrantz oxygen absortion with overlap correction\n"
         "                             on the basis of Drayson routine\n"
         "   Rosenkranz_Voigt_Kuntz6 : Rosenkrantz oxygen absortion with overlap correction\n"
         "                             on the basis of Kuntz routine, accuracy > 2x10^(-6)\n"
         "   CO2_Lorentz:              Lorentz multiplied with Cousin's chi factors\n"
         "   CO2_Drayson:              Drayson multiplied with Cousin's chi factors\n"
         "\n"
         "forefactor:\n"
         "   no_norm:                  1\n"
         "   quadratic:                (f/f0)^2\n"
         "   VVH:                      (f*tanh(h*f/(2k*T))) / (f0*tanh(h*f0/(2k*T)))\n"
         "\n"
         "cutoff:\n"
         "    -1:                      no cutoff\n"
         "   <Number>:                 positive cutoff frequency in Hz\n"
         ),
        AUTHORS( "Axel von Engeln", "Stefan Buehler" ),
        OUT( "abs_lineshape" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN( "shape",  "normalizationfactor", "cutoff" ),
        GIN_TYPE(    "String", "String",              "Numeric" ),
        GIN_DEFAULT( NODEF,    NODEF,                 NODEF ),
        GIN_DESC( "Line shape function.",
                  "Normalization factor.",
                  "Cutoff frequency [Hz]." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_lineshape_per_tgDefine" ),
        DESCRIPTION
        (
         "Set the lineshape, separately for each absorption species.\n"
         "\n"
         "This method is similar to *abs_lineshapeDefine*, except that a\n"
         "different lineshape can be set for each absorption species (see\n"
         "*abs_species*). For example, you might want to use different values of\n"
         "the cutoff frequency for different species.\n"
         "\n"
         "For detailed documentation on the available options for the input\n"
         "parameters see documentation of method *abs_lineshapeDefine*.\n"
         ),
        AUTHORS( "Axel von Engeln", "Stefan Buehler" ),
        OUT( "abs_lineshape" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species" ),
        GIN( "shape",        "normalizationfactor", "cutoff" ),
        GIN_TYPE(    "ArrayOfString", "ArrayOfString",        "Vector" ),
        GIN_DEFAULT( NODEF,          NODEF,                 NODEF ),
        GIN_DESC( "Line shape function for each species.",
                  "Normalization factor for each species.",
                  "Cutoff frequency [Hz] for each species." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_linesReadFromArts" ),
        DESCRIPTION
        (
         "Read all the lines from an Arts catalogue file in the\n"
         "given frequency range. Otherwise a runtime error will be\n"
         "thrown\n"
         "\n"
         "Please note that all lines must correspond\n"
         "to legal species / isotope combinations\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_lines" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN( "filename", "fmin",    "fmax" ),
        GIN_TYPE(    "String",   "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF,      NODEF,     NODEF ),
        GIN_DESC( "Name (and path) of the catalogue file.",
                  "Minimum frequency for lines to read [Hz].",
                  "Maximum frequency for lines to read [Hz]." )
        ));

  md_data_raw.push_back
  ( MdRecord
   ( NAME( "abs_linesReadFromSplitArtscat" ),
    DESCRIPTION
    (
     "Read all the lines in the given frequency range from a split\n"
     "Arts catalogue file.\n"
     "\n"
     "Please note that all lines must correspond\n"
     "to legal species / isotope combinations\n"
     ),
    AUTHORS( "Oliver Lemke" ),
    OUT( "abs_lines" ),
    GOUT(),
    GOUT_TYPE(),
    GOUT_DESC(),
    IN( "abs_species" ),
    GIN(         "basename", "fmin",    "fmax" ),
    GIN_TYPE(    "String",   "Numeric", "Numeric" ),
    GIN_DEFAULT( NODEF,      NODEF,     NODEF ),
    GIN_DESC("Basename of the catalogue.",
             "Minimum frequency for lines to read [Hz].",
             "Maximum frequency for lines to read [Hz]." )
    ));
  
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_linesReadFromHitranPre2004" ),
        DESCRIPTION
        (
         "Read all the lines from a HITRAN 1986-2001 catalogue file in\n"
         "the given frequency range. Otherwise a runtime error will be\n"
         "thrown. For HITRAN 2004 and later line data use the workspace\n"
         "method *abs_linesReadFromHitran*.\n"
         "\n"
         "Please note that all lines must correspond to legal\n"
         "species / isotope combinations and that the line data\n"
         "file must be sorted by increasing frequency\n"
         "\n"
         "WWW access of the HITRAN catalogue: http://www.hitran.com/\n"
         ),
        AUTHORS( "Thomas Kuhn" ),
        OUT( "abs_lines" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN( "filename",  "fmin",    "fmax" ),
        GIN_TYPE(    "String",    "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF,       NODEF,     NODEF ),
        GIN_DESC( "Name (and path) of the catalogue file.",
                  "Minimum frequency for lines to read [Hz].",
                  "Maximum frequency for lines to read [Hz]." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_linesReadFromHitran" ),
        DESCRIPTION
        (
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
         "method: abs_linesReadFromHitranPre2004\n"
         ),
        AUTHORS( "Hermann Berg", "Thomas Kuhn" ),
        OUT( "abs_lines" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN( "filename",  "fmin",    "fmax" ),
        GIN_TYPE( "String", "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF,       NODEF,     NODEF ),
        GIN_DESC( "Name (and path) of the catalogue file.",
                  "Minimum frequency for lines to read [Hz].",
                  "Maximum frequency for lines to read [Hz]." )
        ));
  
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_linesReadFromJpl" ),
        DESCRIPTION
        (
         "Read all the lines from a JPL catalogue file in the\n"
         "given frequency range. Otherwise a runtime error will be\n"
         "thrown\n"
         "\n"
         "Please note that all lines must correspond\n"
         "to legal species / isotope combinations.\n"
         "\n"
         "WWW access of the JPL catalogue: http://spec.jpl.nasa.gov/\n"
         ),
        AUTHORS( "Thomas Kuhn" ),
        OUT( "abs_lines" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN( "filename",  "fmin", "fmax" ),
        GIN_TYPE( "String", "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF,       NODEF,     NODEF ),
        GIN_DESC( "Name (and path) of the catalogue file.",
                  "Minimum frequency for lines to read [Hz].",
                  "Maximum frequency for lines to read [Hz]." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_linesReadFromMytran2" ),
        DESCRIPTION
        (
         "Read all the lines from a MYTRAN2 catalogue file in the\n"
         "given frequency range. Otherwise a runtime error will be\n"
         "thrown\n"
         "\n"
         "Please note that all lines must correspond\n"
         "to legal species / isotope combinations\n"
         ),
        AUTHORS( "Axel von Engeln", "Stefan Buehler" ),
        OUT( "abs_lines" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN( "filename", "fmin", "fmax" ),
        GIN_TYPE( "String", "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF,       NODEF,     NODEF ),
        GIN_DESC( "Name (and path) of the catalogue file.",
                  "Minimum frequency for lines to read [Hz].",
                  "Maximum frequency for lines to read [Hz]." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_lines_per_speciesAddMirrorLines" ),
        DESCRIPTION
        (
         "Adds mirror lines at negative frequencies to *abs_lines_per_species*.\n"
         "\n"
         "For each line at frequency +f in *abs_lines_per_species* a corresponding\n"
         "entry at frequency -f is added to *abs_lines_per_species*. The mirror\n"
         "lines are appended to the line list after the original lines.\n" 
         ),
        AUTHORS( "Axel von Engeln", "Stefan Buehler" ),
        OUT( "abs_lines_per_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_lines_per_species" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_lines_per_speciesCompact" ),
        DESCRIPTION
        (
         "Removes all lines outside the defined lineshape cutoff frequencies\n"
         "from *abs_lines_per_species*. This can save computation time.\n"
         "It should be particularly useful to call this method after\n"
         "*abs_lines_per_speciesAddMirrorLines*.\n" 
         ),
        AUTHORS( "Axel von Engeln", "Stefan Buehler" ),
        OUT( "abs_lines_per_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_lines_per_species", "abs_lineshape", "f_grid" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_lines_per_speciesCreateFromLines" ),
        DESCRIPTION
        (
         "Split lines up into the different species.\n"
         "\n"
         "The species are tested in the order in which they are specified in the\n"
         "controlfile. Lines are assigned to the first species that\n"
         "matches. That means if the list of species is [\"O3-666\",\"O3\"], then\n"
         "the last group O3 gets assigned all the O3 lines that do not fit in\n"
         "the first group (all other isotopes than the main isotope).\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_lines_per_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_lines", "abs_species" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_lines_per_speciesReadFromCatalogues" ),
        DESCRIPTION
        (
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
         "*abs_lines_per_speciesCreateFromLines*.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_lines_per_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species" ),
        GIN(         "filenames",     "formats",       "fmin",   "fmax" ),
        GIN_TYPE(    "ArrayOfString", "ArrayOfString", "Vector", "Vector" ),
        GIN_DEFAULT( NODEF,           NODEF,           NODEF,    NODEF ),
        GIN_DESC( "Name (and path) of the catalogue files.",
                  "Format of each file. (Allowed formats are\n"
                  "HITRAN96, HITRAN04, MYTRAN2, JPL, ARTS.",
                  "Minimum frequency for lines to read [Hz].",
                  "Maximum frequency for lines to read [Hz]." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_lines_per_speciesSetEmpty" ),
        DESCRIPTION
        (
         "Sets abs_lines_per_species to empty line lists.\n"
         "\n"
         "You can use this method to set *abs_lines_per_species* if you do not\n"
         "really want to compute line spectra. Formally, abs_coefCalc will still\n"
         "require *abs_lines_per_species* to be set.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_lines_per_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
  ( MdRecord
   ( NAME( "abs_lines_per_speciesWriteToSplitArtscat" ),
    DESCRIPTION
    (
     "Write each species to a separate catalogue file.\n"
     ),
    AUTHORS( "Oliver Lemke" ),
    OUT(),
    GOUT(),
    GOUT_TYPE(),
    GOUT_DESC(),
    IN( "output_file_format", "abs_lines_per_species" ),
    GIN(         "basename" ),
    GIN_TYPE(    "String" ),
    GIN_DEFAULT( "" ),
    GIN_DESC(    "Basename of the catalogue." )
    ));
  
  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "abs_lookupAdapt" ),
        DESCRIPTION
        (
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
         "always use this method to set it!\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_lookup", "abs_lookup_is_adapted" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_lookup", "abs_species", "f_grid" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "abs_lookupCreate" ),
        DESCRIPTION
        (
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
         "species, the first is assumed to be the main one.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_lookup", "abs_lookup_is_adapted" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species", 
            "abs_lines_per_species",
            "abs_lineshape",
            "abs_nls",
            "f_grid",
            "abs_p",
            "abs_vmrs",
            "abs_t", 
            "abs_t_pert", 
            "abs_nls_pert",
            "abs_n2",
            "abs_cont_names",
            "abs_cont_models", 
            "abs_cont_parameters"
            ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "abs_lookupInit" ),
        DESCRIPTION
        (
         "Creates an empty gas absorption lookup table.\n"
         "\n"
         "This is mainly there to help developers. For example, you can write\n"
         "the empty table to an XML file, to see the file format.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_lookup" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "abs_lookupSetup" ),
        DESCRIPTION
        (
         "Set up input parameters for abs_lookupCreate.\n"
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
         "See also:\n"
         "   *abs_lookupSetupBatch*\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_p",
             "abs_t", 
             "abs_t_pert", 
             "abs_vmrs",
             "abs_nls",
             "abs_nls_pert" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim",
            "p_grid",
            "lat_grid",
            "lon_grid",
            "t_field",
            "vmr_field",
            "abs_species",
            "abs_p_interp_order",
            "abs_t_interp_order",
            "abs_nls_interp_order" ),
        GIN( "p_step",  "t_step",  "h2o_step" ),
        GIN_TYPE(    "Numeric", "Numeric", "Numeric" ),
        GIN_DEFAULT( "0.05",    "100",     "100" ),
        GIN_DESC( /* p_step */
                  "Maximum step in log10(p[Pa]) (base 10 logarithm)."
                  "If the pressure grid is coarser than this, additional "
                  "points are added until each log step is smaller than this.",
                  /* t_step */
                  "The temperature variation grid step in Kelvin, "
                  "for a 2D or 3D atmosphere. For a 1D atmosphere this "
                  "parameter is not used.",
                  /* h2o_step */
                  "The H2O variation grid step [fractional], if "
                  "H2O variations are done (which is determined automatically, "
                  "based on abs_species and the atmospheric dimension). For a "
                  "1D atmosphere this parameter is not used."
                  )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "abs_lookupSetupBatch" ),
        DESCRIPTION
        (
         "Set up input parameters for abs_lookupCreate for batch calculations.\n"
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
         "See also:\n"
         "   *abs_lookupSetup*\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_p",
             "abs_t", 
             "abs_t_pert", 
             "abs_vmrs",
             "abs_nls",
             "abs_nls_pert" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species",
            "batch_atm_fields_compact",
            "abs_p_interp_order",
            "abs_t_interp_order",
            "abs_nls_interp_order" ),
        GIN( "p_step",  "t_step",  "h2o_step", "extremes" ),
        GIN_TYPE(    "Numeric", "Numeric", "Numeric",  "Vector" ),
        GIN_DEFAULT( "0.05",    "20",       "100",      "[]" ),
        GIN_DESC( /* p_step */ 
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
                  "t_pert_max, nls_pert_min, nls_pert_max]."
                  )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "abs_lookupSetupWide" ),
        DESCRIPTION
        (
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
         "default values are chosen such that they cover all Chevallier data set\n"
         "cases, and a bit more. The statistics of the Chevallier data are:\n"
         "\n"
         "min(p)   / max(p)   [Pa]:  1 / 104960\n"
         "min(T)   / max(T)   [K]:   158.21 / 320.39\n"
         "min(H2O) / max(H2O) [VMR]: -5.52e-07 / 0.049\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_p",
             "abs_t", 
             "abs_t_pert", 
             "abs_vmrs",
             "abs_nls",
             "abs_nls_pert" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species",
            "abs_p_interp_order",
            "abs_t_interp_order",
            "abs_nls_interp_order" ),
        GIN( "p_min",   "p_max",   "p_step",  "t_min",   "t_max",   "h2o_min", "h2o_max" ),
        GIN_TYPE(    "Numeric", "Numeric", "Numeric", "Numeric", "Numeric", "Numeric", "Numeric" ),
        GIN_DEFAULT( "0.5",  "110000",  "0.05",    "100",     "400",     "0",       "0.05" ),
        GIN_DESC( "Pressure grid minimum [Pa].",
                  "Pressure grid maximum [Pa].",
                  "Pressure grid step in log10(p[Pa]) (base 10 logarithm).",
                  "Temperature grid minimum [K].",
                  "Temperature grid maximum [K].",
                  "Humidity grid minimum [fractional].",
                  "Humidity grid maximum [fractional]." )
        ));
  
  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "abs_lookupTestAccuracy" ),
        DESCRIPTION
        (
         "Test accuracy of absorption lookup table.\n"
         "\n"
         "Explicitly compare absorption from the lookup table with line-by-line\n"
         "calculations for strategically selected conditions (in-between the\n"
         "lookup table grid points).\n"
         "\n"
         "For error units see *abs_lookupTestAccMC*\n"
         "\n"
         "Produces no workspace output, only output to the output streams.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_lookup",
            "abs_lookup_is_adapted",
            "abs_p_interp_order",
            "abs_t_interp_order",
            "abs_nls_interp_order",
            "abs_n2",
            "abs_lines_per_species", 
            "abs_lineshape", 
            "abs_cont_names", 
            "abs_cont_models", 
            "abs_cont_parameters" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back     
    ( MdRecord
     ( NAME( "abs_lookupTestAccMC" ),
      DESCRIPTION
      (
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
       "Produces no workspace output, only output to the output streams.\n"
       ),
      AUTHORS( "Stefan Buehler" ),
      OUT(),
      GOUT(),
      GOUT_TYPE(),
      GOUT_DESC(),
      IN( "abs_lookup",
         "abs_lookup_is_adapted",
         "abs_p_interp_order",
         "abs_t_interp_order",
         "abs_nls_interp_order",
         "abs_n2",
         "abs_lines_per_species", 
         "abs_lineshape", 
         "abs_cont_names", 
         "abs_cont_models", 
         "abs_cont_parameters",
         "mc_seed"),
      GIN(),
      GIN_TYPE(),
      GIN_DEFAULT(),
      GIN_DESC()
      ));
    
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_n2Set" ),
        DESCRIPTION
        (
         "Sets abs_n2 to the profile of the first tag group containing\n"
         "molecular nitrogen. See *abs_h2oSet* for more details.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_n2" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species", "abs_vmrs" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_scalar_gasCalcLBL" ),
        DESCRIPTION
        (
         "Calculates scalar gas absorption coefficients line-by-line.\n"
         "\n"
         "This method can be used inside *abs_scalar_gas_agenda* just like\n"
         "*abs_scalar_gasExtractFromLookup*. It is a shortcut for putting in some\n"
         "other methods explicitly, namely:\n"
         "\n"
         "  1. *f_gridSelectFIndex*\n"
         "  2. *NumericScale*( rte_doppler, rte_doppler, -1 )\n"
         "  3. *VectorAddScalar*( f_grid, f_grid, rte_doppler )\n"
         "  4. *AbsInputFromRteScalars*\n"
         "  5. *abs_h2oSet*\n"
         "  6. *abs_coefCalc*\n"
         "  7. *abs_scalar_gasFromAbsCoef*\n"
         "\n"
         "Sub-methods 2 and 3 are called only if rte_doppler is not zero.\n"
         "The treatment of the Doppler-shift here is exact, since the underlying\n"
         "frequency grid is shifted.\n"
         "\n"
         "The calculation is for one specific atmospheric condition, i.e., a set\n"
         "of pressure, temperature, VMR values, and Doppler shift. It can be\n"
         "either for a single frequency (f_index>=0), or for all frequencies\n"
         "(f_index<0). The dimension of the output abs_scalar_gas is adjusted\n"
         "accordingly.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_scalar_gas" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "f_grid",
            "abs_species",
            "abs_n2",
            "abs_lines_per_species",
            "abs_lineshape",
            "abs_cont_names",
            "abs_cont_models",
            "abs_cont_parameters",
            "f_index",
            "rte_pressure", "rte_temperature", "rte_vmr_list", "rte_doppler" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_mat_per_speciesCalcLBL" ),
        DESCRIPTION
        (
         "Calculates gas absorption coefficients line-by-line.\n"
         "\n"
         "This method can be used inside *abs_mat_per_species_agenda* just like\n"
         "*abs_mat_per_speciesExtractFromLookup*. It is a shortcut for putting in some\n"
         "other methods explicitly, namely:\n"
         "\n"
         "  1. *f_gridSelectFIndex*\n"
         "  2. *NumericScale*( rte_doppler, rte_doppler, -1 )\n"
         "  3. *VectorAddScalar*( f_grid, f_grid, rte_doppler )\n"
         "  4. *AbsInputFromRteScalars*\n"
         "  5. *abs_h2oSet*\n"
         "  6. *abs_coefCalc*\n"
         "  7. *abs_scalar_gasFromAbsCoef*\n"
         "\n"
         "Sub-methods 2 and 3 are called only if rte_doppler is not zero.\n"
         "The treatment of the Doppler-shift here is exact, since the underlying\n"
         "frequency grid is shifted.\n"
         "\n"
         "The calculation is for one specific atmospheric condition, i.e., a set\n"
         "of pressure, temperature, VMR values, and Doppler shift. It can be\n"
         "either for a single frequency (f_index>=0), or for all frequencies\n"
         "(f_index<0). The dimension of the output abs_mat_per_species is adjusted\n"
         "accordingly.\n"
         ),
        AUTHORS( "Stefan Buehler, Richard Larsson" ),
        OUT( "abs_mat_per_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_mat_per_species",
            "f_grid",
            "abs_species",
            "abs_n2",
            "abs_lines_per_species",
            "abs_lineshape",
            "abs_cont_names",
            "abs_cont_models",
            "abs_cont_parameters",
            "f_index",
            "rte_pressure", "rte_temperature", "rte_vmr_list", "rte_doppler"
           ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_mat_per_speciesInit" ),
        DESCRIPTION
        (
         "Initialize stokes dim gas absorption coefficients line-by-line.\n"
         "\n"
         "This method must be used inside *abs_mat_per_species_agenda*\n"
         ),
        AUTHORS( "Oliver Lemke, Richard Larsson" ),
        OUT( "abs_mat_per_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species",
            "f_grid",
            "f_index",
            "stokes_dim"
        ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
    
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_mat_per_speciesExtractFromLookup" ),
        DESCRIPTION
        (
            //FIXME: Richard
         "Extract gas absorption coefficients from lookup table.\n"
         "\n"
         "This extracts the absorption coefficient for all non-Zeeman species in\n"
         "the current calculation from the lookup table. Extraction is for one\n"
         "specific atmospheric condition, i.e., a set of pressure, temperature,\n"
         "VMR values, and Doppler shift.\n"
         "\n"
         "Extraction can be either for a single frequency (f_index>=0), or for\n"
         "all frequencies (f_index<0). The dimension of the output\n"
         "abs_mat_per_species is adjusted accordingly.\n"
         "\n"
         "The interpolation order in T and H2O is given by *abs_t_interp_order*\n"
         "and *abs_nls_interp_order*, respectively.\n"
         "\n"
         "Note that the treatment of the Doppler-shift here is approximate, since\n"
         "there is a linear interpolation of absorption to a shifted frequency grid.\n"
         "Due to this, with Doppler shift there will be an extrapolation on one edge\n"
         "of the grid, where the spectrum is pushed out of the calculated range.\n"
         "Use extpolfac to control how much extrapolation to tolerate before throwing\n"
         "a runtime error. Default is to allow ten times the outermost grid distance.\n"
         "\n"
         "See also: *abs_mat_per_speciesCalcLBL*.\n"
         ),
        AUTHORS( "Stefan Buehler, Richard Larsson" ),
        OUT( "abs_mat_per_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_mat_per_species", "abs_lookup", "abs_lookup_is_adapted",
            "abs_p_interp_order", "abs_t_interp_order", "abs_nls_interp_order",
            "f_index",
            "rte_pressure", "rte_temperature", "rte_vmr_list", "rte_doppler", "stokes_dim" ),
        GIN("extpolfac"),
        GIN_TYPE("Numeric"),
        GIN_DEFAULT("10"),
        GIN_DESC("Extrapolation factor (for grid edge).")
        ));
    
  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "abs_mat_fieldCalc" ),
        DESCRIPTION
        (
         "Calculate gas absorption for all points in the atmosphere.\n"
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
         "Because of the different contexts, the method can calculate absorption\n"
         "either for all frequencies in the frequency grid (f_index<0), or just\n"
         "for the frequency indicated by f_index (f_index>=0).\n"
         "\n"
         "The calculation itself is performed by the\n"
         "*abs_mat_per_species_agenda*.\n"
         ),
        AUTHORS( "Stefan Buehler, Richard Larsson" ),
        OUT( "abs_mat_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_mat_per_species_agenda",
            "f_index",
            "f_grid",
            "atmosphere_dim",
            "p_grid", "lat_grid", "lon_grid",
            "t_field", "vmr_field" ),
        GIN("doppler", "stokes_dim"),
        GIN_TYPE("Vector", "Index"),
        GIN_DEFAULT("[]", "1"),
        GIN_DESC("A vector of doppler shift values in Hz. Must either be\n"
                 "empty or have same dimension as p_grid\n", "stokes_dim since scalar case is treated differently.")
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_speciesAdd" ),
        DESCRIPTION
        (
         "Adds species tag groups to the list of absorption species.\n"
         "\n"
         "This WSM is similar to *SpeciesSet*, the only difference is that\n"
         "this method appends species to an existing list of absorption species instead\n"
         "of creating the whole list.\n"
         "\n"
         "See *SpeciesSet* for details on how tags are defined and examples of\n"
         "how to input them in the control file.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species" ),
        GIN( "species" ),
        GIN_TYPE(    "ArrayOfString"   ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Specify one String for each tag group that you want to\n"
                  "add. Inside the String, separate the tags by commas\n"
                  "(plus optional blanks).\n")
        ));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_speciesAdd2" ),
        DESCRIPTION
        (
         "Adds a species tag group to the list of absorption species and\n"
         "jacobian quantities.\n"
         "\n"
         "The method is basically a combined call of *abs_speciesAdd* and\n"
         "*jacobianAddAbsSpecies*. In this way it is not needed to specify a\n"
         "tag group in two different places.\n"
         "\n"
         "Arguments exactly as for *jacobianAddAbsSpecies*. Note that this\n"
         "method only handles a single tag group, in contrast to\n"
         "*abs_speciesAdd*\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "abs_species", "jacobian_quantities", "jacobian_agenda" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species", "atmosphere_dim", "p_grid", "lat_grid", 
            "lon_grid" ),
        GIN( "gin1"      , "gin2"      , "gin3"      ,
             "species", "method", "unit", "dx" ),
        GIN_TYPE(    "Vector", "Vector", "Vector",
                     "String", "String", "String", "Numeric" ),
        GIN_DEFAULT( NODEF   , NODEF   , NODEF   ,
                     NODEF,     NODEF,    NODEF,  NODEF ),
        GIN_DESC( "Pressure retrieval grid.",
                  "Latitude retrieval grid.",
                  "Longitude retreival grid.",
                  "The species tag of the retrieval quantity.",
                  "Calculation method. See above.",
                  "Retrieval unit. See above.",
                  "Size of perturbation."
                  ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false  ),
        USES_TEMPLATES( false ),
        PASSWORKSPACE(  true )
        ));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_speciesDefineAllInScenario" ),
        DESCRIPTION
        (
         "Define one tag group for each species known to ARTS and included in an\n"
         "atmospheric scenario.\n"
         "\n"
         "You can use this as an alternative to *SpeciesSet* if you want to make an\n"
         "absorption calculation that is as complete as possible. The method\n"
         "goes through all defined species and tries to open the VMR file. If\n"
         "this works the tag is included, otherwise it is skipped.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN( "basename" ),
        GIN_TYPE( "String" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "The name and path of a particular atmospheric scenario.\n"
                  "For example: /pool/lookup2/arts-data/atmosphere/fascod/tropical" )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_speciesInit" ),
        DESCRIPTION
        (
         "Sets  *abs_species* to be empty.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_vecAddGas" ),
        DESCRIPTION
        (
         "Add gas absorption to first element of absorption vector.\n"
         "\n"
         "The task of this method is to sum up the gas absorption of the\n"
         "different gas species and add the result to the first element of the\n"
         "absorption vector.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_vec" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_vec", "abs_mat_per_species" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_vecAddPart" ),
        DESCRIPTION
        (
         "The particle absorption is added to *abs_vec*\n"
         "\n"
         "This function sums up the absorption vectors for all particle\n"
         "types weighted with particle number density.\n"
         "The resluling absorption vector is added to the workspace\n"
         "variable *abs_vec*\n"
         "Output and input of this method is *abs_vec* (stokes_dim).\n"
         "The inputs are the absorption vector for the single particle type\n"
         "*abs_vec_spt* (N_particletypes, stokes_dim) and the local particle\n"
         " number densities for all particle types namely the\n"
         "*pnd_field* (N_particletypes, p_grid, lat_grid, lon_grid, ) for given\n"
         "*p_grid*, *lat_grid*, and *lon_grid*. The particle types required\n"
         "are specified in the control file.\n"
         ),
        AUTHORS( "Sreerekha T.R." ),
        OUT( "abs_vec" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_vec", "abs_vec_spt", "pnd_field", "atmosphere_dim",
            "scat_p_index",  "scat_lat_index", "scat_lon_index" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_vecInit" ),
        DESCRIPTION
        (
         "Initialize absorption vector.\n"
         "\n"
         "This method is necessary, because all other absorption methods just\n"
         "add to the existing absorption vector.\n"
         "\n"
         "So, here we have to make it the right size and fill it with 0.\n"
         "\n"
         "Note, that the vector is not really a vector, because it has a\n"
         "leading frequency dimension.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_vec" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "f_grid", "stokes_dim", "f_index" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_xsec_per_speciesAddConts" ),
        DESCRIPTION
        (
         "Calculate absorption cross sections per tag group for continua.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_xsec_per_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species", "f_grid", "abs_p", "abs_t", "abs_n2", "abs_h2o",
            "abs_vmrs", "abs_cont_names", "abs_cont_parameters",
            "abs_cont_models" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_xsec_per_speciesAddLines" ),
        DESCRIPTION
        (
         "Calculates the line spectrum for each tag group and adds\n"
         "it to abs_xsec_per_species.\n"
         ),
        AUTHORS( "Stefan Buehler", "Axel von Engeln" ),
        OUT( "abs_xsec_per_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species", "f_grid", "abs_p", "abs_t", "abs_h2o",
            "abs_vmrs", "abs_lines_per_species", "abs_lineshape" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_xsec_per_speciesInit" ),
        DESCRIPTION
        (
         "Initialize *abs_xsec_per_species*.\n"
         "\n"
         "The initialization is\n"
         "necessary, because methods *abs_xsec_per_speciesAddLines*\n"
         "and *abs_xsec_per_speciesAddConts* just add to *abs_xsec_per_species*.\n"
         "The size is determined from *tgs*.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "abs_xsec_per_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species", "f_grid", "abs_p" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AgendaExecute" ),
        DESCRIPTION
        ( 
         "Execute an agenda.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "a" ),
        GIN_TYPE(    "Agenda" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC(    "Agenda to be executed." ),
        SETMETHOD(    false ),
        AGENDAMETHOD( false )
        ));
      
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AgendaAppend" ),
        DESCRIPTION
        ( 
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
         "the right output WSVs.\n"
          ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(        "aout" ),
        GOUT_TYPE(   "Agenda" ),
        GOUT_DESC(   "Target agenda." ),
        IN(),
        GIN(         "ain" ),
        GIN_TYPE(    "Agenda" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC(    "Source agenda." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   true  ),
        USES_TEMPLATES( false ),
        PASSWORKSPACE(  false ),
        PASSWSVNAMES(   true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AgendaSet" ),
        DESCRIPTION
        ( 
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
         "produce the right output WSVs.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(      "a"       ),
        GOUT_TYPE( "Agenda" ),
        GOUT_DESC( "The new agenda." ),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC(),
        SETMETHOD(      false ),
        AGENDAMETHOD(   true  ),
        USES_TEMPLATES( false ),
        PASSWORKSPACE(  false ),
        PASSWSVNAMES(   true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AntennaConstantGaussian1D" ),
        DESCRIPTION
        (
         "Sets up a 1D gaussian antenna response and a matching\n"
         "*mblock_za_grid*.\n"
         "\n"
         "As *antenna_responseGaussian*, but alsp creates *mblock_za_grid*.\n"
         "For returned antenna response, see *antenna_responseGaussian*.\n"
         "\n"
         "The length of *mblock_za_grid* is determined by *n_za_grid*.\n"
         "The end points of the grid are set to be the same as for the\n"
         "antenna response. The spacing of the grid follows the magnitude of\n"
         "the response; the spacing is smaller where the response is high.\n"
         "More precisely, the grid points are determined by dividing\n"
         "the cumulative sum of the response in equal steps. This makes sense\n"
         "if the representation error of the radiance (as a function of\n"
         " zenith angle) increases linearly with the grid spacing.\n"
         "\n"
         "The WSV *antenna_los* is set to 0.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "antenna_dim", "mblock_za_grid", "mblock_aa_grid", 
             "antenna_response", "antenna_los" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(  ),
        GIN( "n_za_grid", "fwhm", "xwidth_si", "dx_si" ),
        GIN_TYPE( "Index", "Numeric", "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF, NODEF, "3", "0.1" ),
        GIN_DESC( "Number of poits to include in*mblock_za_grid*.",
                  "Full width at half-maximum", 
                  "Half-width of response, in terms of std. dev.", 
                  "Grid spacing, in terms of std. dev." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AntennaMultiBeamsToPencilBeams" ),
        DESCRIPTION
        (
         "Maps a multi-beam case to a matching pencil beam case.\n"
         "\n"
         "Cases with overlapping beams are most efficiently handled by\n"
         "letting *antenna_los* have several rows. That is, there are\n"
         "multiple beams for each measurement block. The drawback is that\n"
         "many variables must be adjusted if the corresponding pencil beam\n"
         "spectra shall be calculated. This method makes this adjustment.\n"
         "That is, if you have a control file for a multiple beam case and\n"
         "for some reason want to avoid the antenna weighting, you add this\n"
         "method before *sensor_responseInit*, and remove the call of\n"
         "*sensor_responseAntenna* and you will get the matching pencil beam\n"
         "spectra.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "sensor_pos", "sensor_los", "antenna_los", "antenna_dim", 
             "mblock_za_grid", "mblock_aa_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "sensor_pos", "sensor_los", "antenna_los", "antenna_dim", 
            "mblock_za_grid", "mblock_aa_grid", "atmosphere_dim" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AntennaOff" ),
        DESCRIPTION
        (
         "Sets some antenna related variables\n"
         "\n"
         "Use this method to set *antenna_dim*, *mblock_za_grid* and\n"
         "*mblock_aa_grid* to suitable values (1, [0] and [], respectively)\n"
         "for cases when a sensor is included, but the antenna pattern is\n"
         "neglected.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "antenna_dim", "mblock_za_grid", "mblock_aa_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AntennaSet1D" ),
        DESCRIPTION
        (
         "Sets the antenna dimension to 1D.\n"
         "\n"
         "Sets *antenna_dim* to 1 and sets *mblock_aa_grid* to be empty.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "antenna_dim", "mblock_aa_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AntennaSet2D" ),
        DESCRIPTION
        (
         "Sets the antenna dimension to 2D.\n"
         "\n"
         "Sets *antenna_dim* to 2.\n"
         "\n"
         "It is only allowed to set *antenna_dim* to 2 when *atmosphere_dim*\n"
         "equals 3.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "antenna_dim" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "antenna_responseGaussian" ),
        DESCRIPTION
        (
         "Sets up a gaussian antenna response.\n"
         "\n"
         "The method assumes that the response is the same for all\n"
         "frequencies and polarisations, and that it can be modelled as\n"
         "gaussian.\n"
         "\n"
         "The grid generated can be written as\n"
         "   si * [-xwidth_si:dx_si:xwidth_si]\n"
         "where si is the standard deviation corresponding to the FWHM.\n"
         "That is, width and spacing of the grid is specified in terms of\n"
         "number of standard deviations. If xwidth_si is set to 2, the\n"
         "response will cover about 95% the complete response. For\n"
         "xwidth_si=3, about 99% is covered.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "antenna_response" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( ),
        GIN( "fwhm", "xwidth_si", "dx_si" ),
        GIN_TYPE( "Numeric", "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF, "3", "0.1" ),
        GIN_DESC( "Full width at half-maximum", 
                  "Half-width of response, in terms of std. dev.", 
                  "Grid spacing, in terms of std. dev." )
        ));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Append" ),
        DESCRIPTION
        (
         "Append a workspace variable to another workspace variable.\n"
         "\n"
         "This method can append an array to an array of the same type,\n"
         "e.g. ArrayOfIndex to ArrayOfIndex. Or a single element to an array\n"
         "such as a Tensor3 to an ArrayOfTensor3.\n"
         "\n"
         "In addition to that Vector and Matrices are also supported.\n"
         "For Matrices the third argument indicates the dimension to\n"
         "append to. 'leading' means to append to the leftmost dimension\n"
         "(row-wise for Matrix), 'trailing' to the rightmost dimension\n"
         "(right-most for matrices).\n"
         "\n"
         "This method is not implemented for all types, just for those where an\n"
         "append makes sense. (See variable list below.).\n"
         ),
        AUTHORS( "Stefan Buehler, Oliver Lemke" ),
        OUT(),
        GOUT( "out" ),
        GOUT_TYPE( "Vector, Matrix, String, " +
                   ARRAY_GROUPS + ", " + ARRAY_GROUPS_WITH_BASETYPE ),
        GOUT_DESC( "The variable to append to." ),
        IN(),
        GIN( "in",
             "dimension" ),
        GIN_TYPE( "Vector, Matrix, String, " +
                  ARRAY_GROUPS + "," + GROUPS_WITH_ARRAY_TYPE,
                  "String" ),
        GIN_DEFAULT( NODEF,
                     "leading" ),
        GIN_DESC( "The variable to append.",
                  "Where to append. Could be either the \"leading\" or \"trailing\" dimension." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ArrayOfIndexSet" ),
        DESCRIPTION
        (
         "Creates an ArrayOfIndex from the given list of numbers.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(      "ai"       ),
        GOUT_TYPE( "ArrayOfIndex" ),
        GOUT_DESC( "Variable to initialize." ),
        IN(),
        GIN(         "values" ),
        GIN_TYPE(    "ArrayOfIndex" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Indexes for initializiation." ),
        SETMETHOD( true )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ArrayOfIndexSetConstant" ),
        DESCRIPTION
        (
         "Creates an ArrayOfIndex of length *nelem*, with all values\n"
         "identical.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "ai"       ),
        GOUT_TYPE( "ArrayOfIndex" ),
        GOUT_DESC( "Variable to initialize." ),
        IN( "nelem" ),
        GIN(         "value" ),
        GIN_TYPE(    "Index" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Array value.." ),
        SETMETHOD( true )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ArrayOfStringSet" ),
        DESCRIPTION
        (
         "Sets a String array according the given text.\n"
         "The format is text = [\"String1\",\"String2\",...]\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(      "as" ),
        GOUT_TYPE( "ArrayOfString" ),
        GOUT_DESC( "Variable to initialize." ),
        IN(),
        GIN( "text" ),
        GIN_TYPE(    "ArrayOfString" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Strings for initialization." ),
        SETMETHOD( true )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Arts" ),
        DESCRIPTION
        ( 
         "Runs the agenda that is specified inside the curly braces. ARTS\n"
         "controlfiles must define this method. It is executed automatically\n"
         "when ARTS is run on the controlfile and cannot be called by the user.\n"
         "This methods was used for Arts 1 controlfiles and is now obsolete.\n"
         "See *Arts2*\n"
          ),
        AUTHORS( "Stefan Buehler" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC(),
        SETMETHOD(    false ),
        AGENDAMETHOD( true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Arts2" ),
        DESCRIPTION
        ( 
         "Runs the agenda that is specified inside the curly braces. ARTS\n"
         "controlfiles must define this method. It is executed automatically\n"
         "when ARTS is run on the controlfile and cannot be called by the user.\n"
          ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC(),
        SETMETHOD(    false ),
        AGENDAMETHOD( true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AtmFieldsCalc" ),
        DESCRIPTION
        (
         "Interpolation of raw atmospheric fields.\n"
         "\n"
         "An atmospheric scenario includes the following data for each\n"
         "position (pressure, latitude, longitude) in the atmosphere:\n"
         "   1. temperature field\n"
         "   2. the corresponding altitude field\n"
         "   3. vmr fields for the gaseous species\n"
         "This method interpolates the fields of raw data (*t_field_raw*,\n"
         "*z_field_raw*) which can be stored on arbitrary\n"
         "grids to the calculation grids (*p_grid*, *lat_grid*, *lon_grid*).\n"
         "\n"
         "With parameter interp_order you can control the order of \n"
         "interpolation. The default is 1 (linear interpolation).\n"
         ),
        AUTHORS( "Claudia Emde", "Stefan Buehler" ),
        OUT( "t_field", "z_field", "vmr_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "p_grid", "lat_grid", "lon_grid", "t_field_raw", "z_field_raw", 
            "vmr_field_raw", "atmosphere_dim" ),
        GIN( "interp_order" ),
        GIN_TYPE( "Index" ),
        GIN_DEFAULT( "1" ),
        GIN_DESC( "Interpolation order." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AtmFieldsCalcExpand1D" ),
        DESCRIPTION
        (
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
         "ellipsoid is set to be a sphere.\n"
         "\n"
         "With parameter interp_order you can control the order of \n"
         "interpolation. The default is 1 (linear interpolation).\n"
         ),
        AUTHORS( "Patrick Eriksson", "Claudia Emde", "Stefan Buehler" ),
        OUT( "t_field", "z_field", "vmr_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "p_grid", "lat_grid", "lon_grid", "t_field_raw", "z_field_raw", 
            "vmr_field_raw", "atmosphere_dim" ),
        GIN( "interp_order" ),
        GIN_TYPE( "Index" ),
        GIN_DEFAULT( "1" ),
        GIN_DESC( "Interpolation order." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AtmFieldsExpand1D" ),
        DESCRIPTION
        (
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
         "ellipsoid is set to be a sphere.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "t_field", "z_field", "vmr_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "t_field", "z_field", "vmr_field", "p_grid", "lat_grid", 
            "lon_grid", "atmosphere_dim" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AtmFieldsRefinePgrid" ),
        DESCRIPTION
        (
         "Refine the pressure grid in the atmospheric fields.\n"
         "\n"
         "This method is used for absorption lookup table testing. It probably\n"
         "has no other application.\n"
         "\n"
         "It adds additional vertical grid points to the atmospheric fields, by\n"
         "interpolating them in the usual ARTS way (linear in log pressure).\n"
         "\n"
         "How fine the new grid will be is determined by the keyword parameter\n"
         "p_step. The definition of p_step, and the interpolation behavior, is\n"
         "consistent with *abs_lookupSetup* and *abs_lookupSetupBatch*. (New\n"
         "points are added between the original ones, so that the spacing is\n"
         "always below p_step.)\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "p_grid",
             "t_field", "z_field", "vmr_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "p_grid", "lat_grid", "lon_grid",
            "t_field", "z_field", "vmr_field", "atmosphere_dim" ),
        GIN( "p_step" ),
        GIN_TYPE(    "Numeric" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC("Maximum step in log(p[Pa]) (natural logarithm, as always). If\n"
                 "the pressure grid is coarser than this, additional points\n"
                 "are added until each log step is smaller than this.\n")
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "atm_fields_compactAddConstant" ),
        DESCRIPTION
        (
         "Adds a constant field to atm_fields_compact.\n"
         "\n"
         "This is handy for nitrogen or oxygen. The constant value is\n"
         "appended at the end of the fields that are already there. All\n"
         "dimensions (pressure, latitude, longitude) are filled up, so this\n"
         "works for 1D, 2D, or 3D atmospheres.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "atm_fields_compact" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atm_fields_compact" ),
        GIN( "name",   "value" ),
        GIN_TYPE(    "String", "Numeric" ),
        GIN_DEFAULT( NODEF,    NODEF ),
        GIN_DESC( "Name of additional atmospheric field, with constant value.",
                  "Constant value of additional field." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "atm_fields_compactAddSpecies" ),
        DESCRIPTION
        (
         "Adds a field to atm_fields_compact, with interpolation.\n"
         "\n"
         "This method appends a *GriddedField3* to *atm_fields_compact*.\n"
         "The *GriddedField3* is interpolated upon the grid of *atm_fields_compact*.\n"
         "A typical use case for this method may be to add a climatology of some gas\n"
         "when this gas is needed for radiative transfer calculations, but\n"
         "not yet present in *atm_fields_compact*. One case where this happens\n"
         "is when using the Chevalier dataset for infrared simulations.\n"
         "\n"
         "The grids in *atm_fields_compact* must fully encompass the grids in\n"
         "the *GriddedField3* to be added, for interpolation to succeed. If\n"
         "this is not the case, a RuntimeError is thrown.\n"
         ),
        AUTHORS( "Gerrit Holl" ),
        OUT( "atm_fields_compact" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atm_fields_compact" ),
        GIN( "name",   "value" ),
        GIN_TYPE(    "String", "GriddedField3" ),
        GIN_DEFAULT( NODEF,    NODEF ),
        GIN_DESC( "Name of additional atmospheric field.",
                  "Value of additional atmospheric field." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "batch_atm_fields_compactAddConstant" ),
        DESCRIPTION
        (
         "Adds a constant field to batch_atm_fields_compact.\n"
         "\n"
         "Applies *atm_fields_compactAddConstant* to each batch.\n"
         "The format is equal to that WSM.\n"
         ),
        AUTHORS( "Gerrit Holl" ),
        OUT( "batch_atm_fields_compact" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "batch_atm_fields_compact" ),
        GIN( "name",   "value" ),
        GIN_TYPE(    "String", "Numeric" ),
        GIN_DEFAULT( NODEF,    NODEF ),
        GIN_DESC( "Name of additional atmospheric field, with constant value.",
                  "Constant value of additional field." )
        ));


   md_data_raw.push_back
    ( MdRecord
      ( NAME( "batch_atm_fields_compactAddSpecies" ),
        DESCRIPTION
        (
         "Adds a field to *batch_atm_fields_compact*, with interpolation.\n"
         "\n"
         "This method appends a *GriddedField3* to each *atm_fields_compact*.\n"
         "in *batch_atm_fields_compact*. For details, see *atm_fields_compactAddSpecies*.\n"
         ),
        AUTHORS( "Gerrit Holl" ),
        OUT( "batch_atm_fields_compact" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "batch_atm_fields_compact" ),
        GIN( "name",   "value" ),
        GIN_TYPE(    "String", "GriddedField3" ),
        GIN_DEFAULT( NODEF,    NODEF ),
        GIN_DESC( "Name of additional atmospheric field. Use, e.g., vmr_ch4 for methane VMR",
                  "Value of additional atmospheric field." )
        ));


  md_data_raw.push_back
    ( MdRecord
      ( NAME( "atm_fields_compactFromMatrix" ),
        DESCRIPTION
        (
         "Set *atm_fields_compact* from 1D profiles in a matrix.\n"
         "\n"
         "For clear-sky batch calculations it is handy to store atmospheric\n"
         "profiles in an array of matrix. We take such a matrix, and create\n"
         "*atm_fields_compact* from it.\n"
         "\n"
         "The matrix must contain one row for each pressure level.\n"
         "The matrix can contain some additional fields which are not directly used\n"
         "by ARTS for calculations but can be required for further processing,\n"
         "for e.g. wind speed and direction. In this case, additional fields must\n"
         "be put at the end of the matrix and they must be flagged by 'ignore',\n"
         "large or small letters, in the field names.\n"
         "Recommended row format:\n"
         "\n"
         "p[Pa] T[K] z[m] VMR_1[fractional] ... VMR[fractional] IGNORE ... IGNORE\n"
         "\n"
         "Works only for *atmosphere_dim==1.*\n"         
         "\n"
         "Keywords:\n"
         "   field_names : Field names to store in atm_fields_compact.\n"
         "                 This should be, e.g.:\n"
         "                 [\"T[K]\", \"z[m]\", \"vmr_h2o[fractional]\", \"ignore\"]\n"
         "                 There must be one name less than matrix columns,\n"
         "                 because the first column must contain pressure.\n"
         ),
        AUTHORS( "Stefan Buehler", "Daniel Kreyling", "Jana Mendrok" ),
        OUT( "atm_fields_compact" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim" ),
        GIN(      "gin1"      ,
                  "field_names" ),
        GIN_TYPE(    "Matrix",
                     "ArrayOfString" ),
        GIN_DEFAULT( NODEF   ,
                     NODEF ),
        GIN_DESC( "One atmosphere matrix from batch input ArrayOfMatrix.",
                  "Order/names of atmospheric fields." )
        ));
    

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AtmFieldsFromCompact" ),
        DESCRIPTION
        (
         "Extract pressure grid and atmospheric fields from\n"
         "*atm_fields_compact*.\n"
         "\n"
         "An atmospheric scenario includes the following data for each\n"
         "position (pressure, latitude, longitude) in the atmosphere:\n"
         "           1. temperature field\n"
         "           2. the corresponding altitude field\n"
         "           3. vmr fields for the gaseous species\n"
         "\n"
         "This method just splits up the data found in *atm_fields_compact* to\n"
         "p_grid, lat_grid, lon_grid, and the various fields. No interpolation.\n"
         "See documentation of *atm_fields_compact* for a definition of the data.\n"
         "\n"
         "There are some safety checks on the names of the fields: The first\n"
         "field must be called \"T\", the second \"z\"*. Remaining fields must be\n"
         "trace gas species volume mixing ratios, named for example \"H2O\", \"O3\",\n"
         "and so on. The species names must fit the species in *abs_species*.\n"
         "(Same species in same order.) Only the species name must fit, not the\n"
         "full tag.\n"
         "\n"
         "Possible future extensions: Add a keyword parameter to refine the\n"
         "pressure grid if it is too coarse. Or a version that interpolates onto\n"
         "given grids, instead of using and returning the original grids.\n"
         ),
        AUTHORS( "Stefan Buehler", "Daniel Kreyling", "Jana Mendrok" ),
        OUT( "p_grid", "lat_grid", "lon_grid", "t_field", "z_field",
             "vmr_field", "massdensity_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species", "part_species", "atm_fields_compact", "atmosphere_dim" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
    
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AtmosphereSet1D" ),
        DESCRIPTION
        (
         "Sets the atmospheric dimension to 1D.\n"
         "\n"
         "Sets *atmosphere_dim* to 1 and gives some variables dummy values.\n"
         "\n"
         "The latitude and longitude grids are set to be empty.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "atmosphere_dim", "lat_grid", "lon_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AtmosphereSet2D" ),
        DESCRIPTION
        (
         "Sets the atmospheric dimension to be 2D.\n"
         "\n"
         "Sets *atmosphere_dim* to 2 and the longitude grid to be empty.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "atmosphere_dim", "lon_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AtmosphereSet3D" ),
        DESCRIPTION
        (
         "Sets the atmospheric dimension to 3D.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "atmosphere_dim" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AtmRawRead" ),
        DESCRIPTION
        (
         "Reads atmospheric data from a scenario.\n"
         "\n"
         "An atmospheric scenario includes the following data for each\n"
         "position (pressure, latitude, longitude) in the atmosphere:\n"
         "   1. temperature field\n"
         "   2. the corresponding altitude field\n"
         "   3. vmr fields for the gaseous species\n"
         "The data is stored in different files. This methods reads all\n"
         "files and creates the variables *t_field_raw*, *z_field_raw* and\n"
         "*vmr_field_raw*.\n"
         "\n"
         "Files in a scenarios should be named matching the pattern of:\n"
         "tropical.H2O.xml\n"
         "\n"
         "The files can be anywhere, but they must be all in the same\n"
         "directory, selected by 'basename'. The files are chosen by the\n"
         "species name. If you have more than one tag group for the same\n"
         "species, the same profile will be used.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "t_field_raw", "z_field_raw", "vmr_field_raw" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species" ),
        GIN( "basename" ),
        GIN_TYPE( "String" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Name of scenario, probably including the full path. For "
                  "example: \"/smiles_local/arts-data/atmosphere/fascod/"
                  "tropical\"" )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "backend_channel_responseFlat" ),
        DESCRIPTION
        (
         "Sets up a rectangular channel response.\n"
         "\n"
         "The response of the backend channels is hee assumed to be constant\n"
         "inside the resolution width, and zero outside.\n"
         "\n"
         "The method assumes that all channels have the same response.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "backend_channel_response" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( ),
        GIN( "resolution" ),
        GIN_TYPE( "Numeric" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "The spectrometer resolution." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "backend_channel_responseGaussian" ),
        DESCRIPTION
        (
         "Sets up a gaussian backend channel response.\n"
         "\n"
         "The method assumes that all channels have the same response, and\n"
         "that it can be modelled as gaussian.\n"
         "\n"
         "The grid generated can be written as\n"
         "   si * [-xwidth_si:dx_si:xwidth_si]\n"
         "where si is the standard deviation corresponding to the FWHM.\n"
         "That is, width and spacing of the grid is specified in terms of\n"
         "number of standard deviations. If xwidth_si is set to 2, the\n"
         "response will cover about 95% the complete response. For\n"
         "xwidth_si=3, about 99% is covered.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "backend_channel_response" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( ),
        GIN( "fwhm", "xwidth_si", "dx_si" ),
        GIN_TYPE( "Numeric", "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF, "3", "0.1" ),
        GIN_DESC( "Full width at half-maximum", 
                  "Half-width of response, in terms of std. dev.", 
                  "Grid spacing, in terms of std. dev." )
        ));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "basics_checkedCalc" ),
        DESCRIPTION
        (
         "Checks consistency of the (clear sky) atmosphere.\n"
         "\n"
         "The following WSVs are treated: f_grid, stokes_dim, p_grid,\n"
         "lat_grid, lon_grid, t_field, z_field, vmr_field, wind_u/v/w_field,\n"
         "refellipsoid and z_surface.\n"
         "If any of these variables are changed, then this method shall be\n"
         "called again (no automatic check that this is fulfilled!).\n"
         "\n"
         "The tests include:\n"
         " 1. That basic control variables *stokes_dim* and *atmosphere_dim*\n"
         "    are inside defined ranges.\n"
         " 2. That *f_grid* is sorted and increasing.\n"
         " 3. If atmospheric grids (p/lat/lon_grid) are OK with respect to\n"
         "    *atmosphere_dim*.\n"
         " 4. If *refellipsoid* has correct size, and that eccentricity is\n"
         "    set to zero if 1D atmosphere.\n"
         " 5. If atmospheric fields, and *z_surface* have sizes consistent\n"
         "    with the atmospheric grids.\n"
         " 6. There is no gap between *z_surface* and *z_field*.\n"
         "\n"
         "If any test fails, there is an error. Otherwise, *basics_checked*\n"
         "is set to 1.\n"
         "\n"
         "The cloudbox is covered by *cloudbox_checked*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "basics_checked" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim", "p_grid", "lat_grid", "lon_grid", "abs_species",
            "z_field", "t_field", "vmr_field", "wind_u_field", "wind_v_field",
            "wind_w_field", "mag_u_field", "mag_v_field", "mag_w_field",
            "edensity_field", "refellipsoid", "z_surface", 
            "stokes_dim", "f_grid" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "batch_atm_fields_compactFromArrayOfMatrix" ),
        DESCRIPTION
        (
         "Expand batch of 1D atmospheric states to a batch_atm_fields_compact.\n"
         "\n"
         "This is used to handle 1D batch cases, for example from the Chevallier\n"
         "data set, stored in a matrix.\n"
         "\n"
         "The matrix must contain one row for each pressure level.\n"
         "The matrix can contain some additional fiels which are not directly used\n"
         "by ARTS for calculations but can be required for further processing,\n"
         "for e.g. wind speed and direction. In this case, additional fields must\n"
         "be put at the end of the matrix and they must be flagged by 'ignore',\n"
         "large or small letters, in the field names.\n"
         "Row format:\n"
         "\n" 
         "p[Pa] T[K] z[m] VMR_1[fractional] ... VMR[fractional] IGNORE ... IGNORE\n"
         "\n"
         "Keywords:\n"
         "   field_names : Field names to store in atm_fields_compact.\n"
         "                 This should be, e.g.:\n"
         "                 [\"T\", \"z\", \"H2O\", \"O3\", \"ignore\"]\n"
         "                 There must be one name less than matrix columns,\n"
         "                 because the first column must contain pressure.\n"
         "\n"
         "   extra_field_names : You can add additional constant VMR fields,\n"
         "                       which is handy for O2 and N2. Give here the\n"
         "                       field name, e.g., \"O2\". Default: Empty.\n"
         "\n"
         "   extra_field_values : Give here the constant field value. Default:\n"
         "                        Empty. Dimension must match extra_field_names.\n"
         ),
        AUTHORS( "Stefan Buehler", "Daniel Kreyling", "Jana Mendrok" ),
        OUT( "batch_atm_fields_compact" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim" ),
        GIN(      "gin1"             ,
                  "field_names", "extra_field_names", "extra_field_values" ),
        GIN_TYPE(    "ArrayOfMatrix",
                     "ArrayOfString", "ArrayOfString",     "Vector" ),
        GIN_DEFAULT( NODEF          ,
                     NODEF,         "[]",                "[]" ),
        //KW_DEFAULT( NODEF,         NODEF,                NODEF ),
        GIN_DESC( "Batch of atmospheres stored in one array of matrix",
                  "Order/names of atmospheric fields.",
                  "Names of additional atmospheric fields, with constant values.",
                  "Constant values of additional fields.")
        ));
    
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "blackbody_radiationPlanck" ),
        DESCRIPTION
        (
         "The Planck function (frequency version).\n"
         "\n"
         "The standard function for *blackbody_radiation_agenda*.\n"
         "\n"
         "The is considered as the standard version inside ARTS of the Planck\n"
         "function. The unit of the returned data is W/(m^2 Hz sr).\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "blackbody_radiation" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "f_grid", "rte_temperature" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
    
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "CloudboxGetIncoming" ),
        DESCRIPTION
        (
         "Calculates incoming radiation field of the cloudbox by repeated\n"
         "radiative transfer calculations.\n"
         "\n"
         "The method performs monochromatic pencil beam calculations for\n"
         "all grid positions on the cloudbox boundary, and all directions\n"
         "given by scattering angle grids (*scat_za/aa_grid*). Found radiances\n"
         "are stored in *scat_i_p/lat/lon* which can be used as boundary\n"
         "conditions when scattering inside the cloud box is solved by the\n"
         "DOIT method.\n"
         ),
        AUTHORS( "Sreerekha T.R.", "Claudia Emde" ),
        OUT( "scat_i_p", "scat_i_lat", "scat_i_lon" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "iy_main_agenda", "atmosphere_dim", "lat_grid", 
            "lon_grid", "z_field", "t_field", "vmr_field", "cloudbox_on", 
            "cloudbox_limits", "basics_checked", "cloudbox_checked", 
            "f_grid", "stokes_dim", "scat_za_grid", "scat_aa_grid" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "CloudboxGetIncoming1DAtm" ),
        DESCRIPTION
        (
         "As *CloudboxGetIncoming* but assumes clear sky part to be 1D."
         "\n"
         "The incoming field is calculated only for one position and azimuth\n"
         "angle for each cloud box boundary, and obtained values are used\n"
         "for all other postions and azimuth angles. This works if a 3D\n"
         "cloud box is put into an 1D background atmosphere.\n"
         "\n"
         "This method can only be used for 3D cases.\n"
         ),
        AUTHORS( "Sreerekha T.R.", "Claudia Emde" ),
        OUT( "scat_i_p", "scat_i_lat", "scat_i_lon", "cloudbox_on" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "iy_main_agenda", "atmosphere_dim", "lat_grid", 
            "lon_grid", "z_field", "t_field", "vmr_field", "cloudbox_on", 
            "cloudbox_limits", "f_grid", "stokes_dim", 
            "scat_za_grid", "scat_aa_grid" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "cloudboxOff" ),
        DESCRIPTION
        (
         "Deactivates the cloud box.\n"
         "\n"
         "Use this method if no scattering calculations shall be performed.\n"
         "The function sets *cloudbox_on* to 0, *cloudbox_limits* to be an\n"
         "empty vector and *iy_cloudbox_agenda* to an empty agenda.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "cloudbox_on", "cloudbox_limits", "iy_cloudbox_agenda" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
    
    
    md_data_raw.push_back
    ( MdRecord
      ( NAME( "cloudboxSetAutomatically" ),
        DESCRIPTION
        (
         "Sets the cloud box to encompass the cloud given by the entries\n"
         "in *massdensity_field*. \n"
         "\n"
         "The function must be called before any *cloudbox_limits* using\n"
         "WSMs.\n"
         "NOTE: only 1-dim case is handeled in the moment!\n"
         "\n"
         "The function iterates over all *part_species* and performs a \n"
         "check, to see if the corresponding scattering particle profiles do not\n"
         "contain a cloud (all values equal zero). If, after all iterations,\n"
         "all the considrered profiles proove to contain no cloud,\n"
         "the cloudbox is switched off! (see WSM *cloudboxOff*)\n"
         "\n"
         "Each scattering particle profile is searched for the first and last\n"
         "pressure index, where the value is unequal to zero. This index\n"
         "is then copied to *cloudbox_limits*.\n"
         "\n"
         "Additionaly the lower cloudbox_limit is altered by\n" 
         "*cloudbox_margin*.\n"
         "The margin is given as a height difference in meters and\n"
         "trasformed into a pressure.(via isothermal barometric heightformula)\n"
         "This alteration is needed to ensure, that scattered photons\n"
         "do not leave and re-enter the cloudbox, due to its convex\n"
         "shape.\n"
         "If *cloudbox_margin* is set to -1 (default), the cloudbox will extend to\n" 
         "the surface. Hence the lower cloudbox_limit is set to 0 (index\n"
         "of first pressure level).\n"
         ),
        AUTHORS( "Daniel Kreyling" ),
        OUT( "cloudbox_on", "cloudbox_limits"),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim", "p_grid", "lat_grid", "lon_grid", "massdensity_field"),
        GIN( "cloudbox_margin"),
        GIN_TYPE( "Numeric" ),
        GIN_DEFAULT( "-1" ),
        GIN_DESC( "The margin alters the lower vertical\n"
                  "cloudbox limit. Value must be given in [m].\n"
                  "If cloudbox_margin is set to *-1* (default), the lower\n" 
                  "cloudbox limit equals 0, what corresponds to the surface !\n"
        )
        ));
  
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "cloudboxSetDisort" ),
        DESCRIPTION
        (
         "For Disort calculation the cloudbox must be extended to\n"
         "cover the full atmosphere.\n"
         "This method sets *cloudbox_limits* accordingly.\n"
         ), 
        AUTHORS( "Claudia Emde" ),
        OUT( "cloudbox_on", "cloudbox_limits" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "p_grid" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));


  md_data_raw.push_back
    ( MdRecord
      ( NAME( "cloudboxSetManually" ),
        DESCRIPTION
        (
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
         "grid end points.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "cloudbox_on", "cloudbox_limits" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim", "p_grid", "lat_grid", "lon_grid" ),
        GIN( "p1",      "p2",      "lat1",    "lat2",    "lon1",
             "lon2" ),
        GIN_TYPE(    "Numeric", "Numeric", "Numeric", "Numeric", "Numeric", 
                     "Numeric" ),
        GIN_DEFAULT( NODEF,     NODEF,     NODEF,     NODEF,     NODEF,
                     NODEF ),
        GIN_DESC( "Upper pressure point.",
                  "Lower pressure point.",
                  "Lower latitude point.",
                  "Upper latitude point.",
                  "Lower longitude point.",
                  "Upper longitude point." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "cloudboxSetManuallyAltitude" ),
        DESCRIPTION
        (
         "Sets the cloud box to encompass the given positions.\n"
         "\n"
         "As *cloudboxSetManually* but uses altitudes instead of pressure.\n"
         "The given altitude points can be outside the range of *z_field*.\n"
         "The altitude limit is then set to the end point of *p_grid*.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "cloudbox_on", "cloudbox_limits" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim", "z_field", "lat_grid", "lon_grid" ),
        GIN( "z1",      "z2",      "lat1",    "lat2",    "lon1",
             "lon2" ),
        GIN_TYPE(    "Numeric", "Numeric", "Numeric", "Numeric", "Numeric", 
                     "Numeric" ),
        GIN_DEFAULT( NODEF,     NODEF,     NODEF,     NODEF,     NODEF,
                     NODEF ),
        GIN_DESC( "Lower altitude point.",
                  "Upper altitude point.",
                  "Lower latitude point.",
                  "Upper latitude point.",
                  "Lower longitude point.",
                  "Upper longitude point." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "cloudbox_checkedCalc" ),
        DESCRIPTION
        (
         "Checks consistency between the cloudbox and other variables.\n"
         "\n"
         "The following WSVs are treated: cloudbox_on, cloudbox_limits and\n"
         "wind_u/v/w_field.\n"
         "If any of these variables are changed, then this method shall be\n"
         "called again (no automatic check that this is fulfilled!).\n"
         "\n"
         "The main check is if the cloudbox limits are OK with respect to\n"
         "the atmospheric dimensionality and the limits of the atmosphere.\n"
         "\n"
         "If any test fails, there is an error. Otherwise, *cloudbox_checked*\n"
         "is set to 1.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "cloudbox_checked" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "basics_checked", "atmosphere_dim", "p_grid", "lat_grid", 
            "lon_grid", "wind_u_field", "wind_v_field", "wind_w_field", 
            "cloudbox_on", "cloudbox_limits" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "complex_nFromGriddedField4" ),
        DESCRIPTION
        (
         "Extracts complex refractive index from a field of such data.\n"
         "\n"
         "This method allows to specify a field of *complex_n* for\n"
         "automatic interpolation to points of interest. The position and\n"
         "direction for which the reflectivity shall be extracted are given\n"
         "by *rte_pos* and *rte_los*. The reflectivity field is expected to\n"
         "be stored as:\n"
         "   GriddedField4:\n"
         "      Vector f_grid[N_f]\n"
         "      Vector real/imaginary[2]\n"
         "      Vector lat_grid[N_lat]\n"
         "      Vector lon_grid[N_lon]\n"
         "      Tensor4 data[N_f][2][N_lat][N_lon]\n"
         "\n"
         "Grids for latitude and longitude must have a length of >= 2 (ie.\n"
         "no automatic expansion). If the frequency grid has length 1, this\n"
         "is taken as the refractive index is constant, following the\n"
         "definition of *complex_n*. The remaining dimension must have size\n"
         "two, where the first element corresponds to the real part and the\n"
         "second elemt to the the imaginary part.\n"
         "\n"
         "The interpolation is done in steps:\n"
         "   1: Linear interpolation for lat and lon (std. extrapolation).\n"
         "   2. Linear interpolation if frequency (if input data have more\n"
         "      than one frequency).\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "complex_n" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "stokes_dim", "f_grid", "atmosphere_dim", "lat_grid", "lat_true", 
            "lon_true", "rte_pos", "rte_los" ),
        GIN( "n_field" ),
        GIN_TYPE( "GriddedField4" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "A field of complex refractive index." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "complex_nWaterLiebe93" ),
        DESCRIPTION
        (
         "Complex refractive index of liquid water according to Liebe 1993.\n"
         "\n"
         "The method treats liquid water without salt. Thus, not valid below\n"
         "10 GHz. Upper frequency limit not known, here set to 1000 GHz.\n"
         "Model parameters taken from Atmlab function epswater93 (by\n"
         "C. Maetzler), which refer to Liebe 1993 without closer\n"
         "specifications.\n"
         "\n"
         "Temperature must be between 0 and 100 degrees Celsius.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "complex_n" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "f_grid" ),
        GIN( "t" ),
        GIN_TYPE( "Numeric" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Temperature [K]." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Copy" ),
        DESCRIPTION
        (
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
         "methods).\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT(),
        GOUT(      "gout1"    ),
        GOUT_TYPE( "Any" ),
        GOUT_DESC( "Destination variable." ),
        IN(),
        GIN(      "gin1"    ),
        GIN_TYPE(    "Any" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Source variable." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  ),
        PASSWORKSPACE(  false ),
        PASSWSVNAMES(   true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Delete" ),
        DESCRIPTION
        (
         "Deletes a workspace variable.\n"
         "\n"
         "The variable is marked as uninitialized and its memory freed.\n"
         "It is not removed from the workspace though.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "v" ),
        GIN_TYPE(    "Any" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC(    "Variable to be deleted." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  ),
        PASSWORKSPACE(  true  ),
        PASSWSVNAMES(   true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "DoitAngularGridsSet" ),
        DESCRIPTION
        (
         "Sets the angular grids for DOIT calculation."
         "\n"
         "In this method the angular grids for a DOIT calculation are\n"
         "specified. For down-looking geometries it is sufficient to define\n"
         "*N_za_grid* and *N_aa_grid*. From these numbers equally spaced\n"
         "grids are created and stored in the WSVs *scat_za_grid* and\n"
         "*scat_aa_grid*.\n" 
         "\n"
         "For limb simulations it is important to use an optimized zenith \n"
         "angle grid with a very fine resolution about 90 degrees. Such a grid can be\n"
         "generated using *doit_za_grid_optCalc*. The filename of an optimized\n"
         "zenith angle grid can be given as a keyword (*za_grid_opt_file*).\n"
         "\n"
         "If a filename is given, the equidistant grid is used for the\n"
         "calculation of the scattering integrals and the optimized grid is\n"
         "applied for integration of the radiative transfer equation. \n"
         "\n"
         "For down-looking cases no filename should be specified (za_grid_opt_file = \"\" ) \n"
         "Using only the equidistant grid makes sense to speed up the calculation.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "doit_za_grid_size", "scat_aa_grid", "scat_za_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN( "N_za_grid", "N_aa_grid", "za_grid_opt_file" ),
        GIN_TYPE(    "Index",     "Index",     "String" ),
        GIN_DEFAULT( NODEF,       NODEF,       NODEF ),
        GIN_DESC( "Number of grid points in zenith angle grid. "
                  "Recommended value is 19.",
                  "Number of grid points in azimuth angle grid. "
                  "Recommended value is 37.",
                  "Name of special grid for RT part." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "DoitCloudboxFieldPut" ),
        DESCRIPTION
        (
         "Method for the DOIT communication between cloudbox and clearsky.\n"
         "\n"
         "This method puts the scattered radiation field into the interface\n"
         "variables between the cloudbox and the clearsky, which are\n"
         "*scat_i_p*, *scat_i_lat* and *scat_i_lon*.\n"
         "\n"
         "The best way to calculate spectra including the influence of\n" 
         "scattering is to set up the *doit_mono_agenda* where this method\n"
         "can be included.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "scat_i_p", "scat_i_lat", "scat_i_lon",
             "doit_i_field1D_spectrum" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "scat_i_p", "doit_i_field", "f_grid", "f_index",   "p_grid", "lat_grid", 
            "lon_grid", "scat_za_grid", "scat_aa_grid", "stokes_dim",
            "atmosphere_dim", "cloudbox_limits", "sensor_pos", "z_field" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "doit_conv_flagAbs" ),
        DESCRIPTION
        (
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
         "cloudbox and for all directions.\n"  
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "doit_conv_flag", "doit_iteration_counter" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "doit_conv_flag", "doit_iteration_counter",
            "doit_i_field", "doit_i_field_old" ),
        GIN( "epsilon" ),
        GIN_TYPE( "Vector" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Limits for convergence. A vector with length matching "
                  "*stokes_dim* with unit [W / (m^2 Hz sr)]."
                  )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "doit_conv_flagLsq" ),
        DESCRIPTION
        (
         "DOIT convergence test (least squares).\n"
         "\n"
         "As *doit_conv_flagAbsBT* but applies a least squares convergence\n"
         "test between two successive iteration fields.\n"
         "\n"
         "Warning: This method is not recommended because this kind of\n"
         "convergence test is not sufficiently strict, so that the\n"
         "DOIT result might be wrong.\n" 
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "doit_conv_flag", "doit_iteration_counter" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "doit_conv_flag", "doit_iteration_counter", 
            "doit_i_field", "doit_i_field_old", "f_grid", "f_index" ),
        GIN( "epsilon" ),
        GIN_TYPE( "Vector" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Limits for convergence. A vector with length matching "
                  "*stokes_dim* with unit [K]."
                  )
        ));
  
  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "doit_conv_flagAbsBT" ),
        DESCRIPTION
        (
         "DOIT convergence test (maximum absolute difference in Rayleigh Jeans "
         "BT)\n"
         "\n"
         "As *doit_conv_flagAbs* but convergence limits are specified in\n"
         "Rayleigh-Jeans brighntess temperatures.\n"
         ),
        AUTHORS( "Sreerekha T.R.", "Claudia Emde" ),
        OUT( "doit_conv_flag", "doit_iteration_counter" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "doit_conv_flag", "doit_iteration_counter",
            "doit_i_field", "doit_i_field_old", "f_grid", "f_index" ),
        GIN( "epsilon" ),
        GIN_TYPE(    "Vector" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Limits for convergence. A vector with length matching "
                  "*stokes_dim* with unit [K]."
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "DoitInit" ),
        DESCRIPTION
        (
         "Initialises variables for DOIT scattering calculations.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "scat_p_index", "scat_lat_index", "scat_lon_index", 
             "scat_za_index", "scat_aa_index", "doit_scat_field",
             "doit_i_field", "doit_is_initialized" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "stokes_dim", "atmosphere_dim", "scat_za_grid", "scat_aa_grid",
            "doit_za_grid_size", "cloudbox_on", "cloudbox_limits", "scat_data_raw" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_i_fieldIterate" ),
        DESCRIPTION
        (
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
         "      supported.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "doit_i_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "doit_i_field", "doit_scat_field_agenda", "doit_rte_agenda", 
            "doit_conv_test_agenda" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_i_fieldSetClearsky" ),
        DESCRIPTION
        (
         "Interpolate clearsky field on all gridpoints in cloudbox.\n"
         "\n"
         "This method uses a linear 1D/3D interpolation scheme to obtain the\n"
         "radiation field on all grid points inside the cloud box from the\n"
         "clear sky field on the cloud box boundary. This radiation field\n"
         "is taken as the first guess radiation field in the DOIT module.\n"
         "\n"
         "Set the *all_frequencies* to 1 if the clearsky field shall be used\n"
         "as initial field for all frequencies. Set it to 0 if the clear sky\n"
         "field shall be used only for the first frequency in *f_grid*. For\n"
         "later frequencies, *doit_i_field* of the previous frequency is then\n"
         "used.\n"
         ),
        AUTHORS( "Sreerekha T.R. and Claudia Emde" ),
        OUT( "doit_i_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "scat_i_p", "scat_i_lat", "scat_i_lon", "f_grid", 
            "f_index", "p_grid", "lat_grid", "lon_grid", 
            "cloudbox_limits", "atmosphere_dim" ),
        GIN( "all_frequencies" ),
        GIN_TYPE( "Index" ),
        GIN_DEFAULT( "1" ),
        GIN_DESC( "See above." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_i_fieldSetConst" ),
        DESCRIPTION
        (
         "This method sets the initial field inside the cloudbox to a\n"
         "constant value. The method works only for monochromatic\n"
         "calculations (number of elements in f_grid=1).\n"
         "\n"
         "The user can specify a value for each Stokes dimension in the\n"
         "control file by *value*.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "doit_i_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "scat_i_p", "scat_i_lat", "scat_i_lon", "p_grid", "lat_grid", 
            "lon_grid", 
            "cloudbox_limits", "atmosphere_dim", "stokes_dim" ),
        GIN( "value" ),
        GIN_TYPE(    "Vector" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "A vector containing 4 elements with the value of the "
                  "initial field for each Stokes dimension."
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_i_fieldUpdate1D" ),
        DESCRIPTION
        (
         "RT calculation in cloudbox with fixed scattering integral (1D).\n"
         "\n"
         "Updates the radiation field (DOIT method). The method loops\n"
         "through the cloudbox to update the radiation field for all\n"
         "positions and directions in the 1D cloudbox.\n"
         "\n"
         "Note: This method is very inefficient, because the number of\n"
         "iterations scales with the number of cloudbox pressure levels.\n"
         "It is recommended to use *doit_i_fieldUpdateSeq1D*.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "doit_i_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "doit_i_field_old", "doit_scat_field", "cloudbox_limits", 
            "abs_mat_per_species_agenda",
            "vmr_field", "spt_calc_agenda", "scat_za_grid", "pnd_field", 
            "opt_prop_part_agenda", 
            "ppath_step_agenda", "p_grid", "z_field", "refellipsoid", 
            "t_field", "edensity_field", "f_grid", "f_index", 
            "surface_rtprop_agenda", "doit_za_interp" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_i_fieldUpdateSeq1D" ),
        DESCRIPTION
        (
         "RT calculation in cloudbox with fixed scattering integral.\n"
         "\n"
         "Updates radiation field (*doit_i_field*) in DOIT module.\n"
         "This method loops through the cloudbox to update the\n"
         "radiation field for all positions and directions in the 1D\n"
         "cloudbox. The method applies the sequential update. For more\n"
         "information refer to AUG.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "doit_i_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "doit_i_field", "doit_scat_field", "cloudbox_limits", 
            "abs_mat_per_species_agenda",
            "vmr_field", "spt_calc_agenda", "scat_za_grid", "pnd_field",
            "opt_prop_part_agenda", 
            "ppath_step_agenda", "p_grid", "z_field", "refellipsoid", 
            "t_field", "edensity_field", "f_grid", "f_index", 
            "surface_rtprop_agenda", "doit_za_interp" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_i_fieldUpdateSeq1DPP" ),
        DESCRIPTION
        (
         "RT calculation in cloudbox with fixed scattering integral.\n"
         "\n " 
         "Update radiation field (*doit_i_field*) in DOIT module.\n"
         "This method loops through the cloudbox to update the\n"
         "radiation field for all\n"
         "positions and directions in the 1D cloudbox. The method applies\n"
         "the sequential update and the plane parallel approximation.\n"
         "This method is only slightly faster than\n"
         "*doit_i_fieldUpdateSeq1D* and it is less accurate. It can not\n"
         "be used for limb simulations.\n"
         ),
        AUTHORS( "Sreerekha T.R." ),
        OUT( "doit_i_field", "scat_za_index" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "doit_scat_field", "cloudbox_limits", 
            "abs_mat_per_species_agenda",
            "vmr_field", "spt_calc_agenda", "scat_za_grid", "pnd_field", 
            "opt_prop_part_agenda",
            "p_grid", "z_field", "t_field", "f_grid", "f_index" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_i_fieldUpdateSeq3D" ),
        DESCRIPTION
        (
         "RT calculation in cloudbox with fixed scattering integral.\n"
         "\n"
         "Update radiation field (*doit_i_field*) in DOIT module.\n"
         "This method loops through the cloudbox to update the\n"
         "radiation field for all positions and directions in the 3D\n"
         "cloudbox. The method applies the sequential update. For more\n"
         "information please refer to AUG.\n"
         "Surface reflections are not yet implemented in 3D scattering\n"
         "calculations.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "doit_i_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "doit_i_field", "doit_scat_field", "cloudbox_limits", 
            "abs_mat_per_species_agenda",
            "vmr_field", "spt_calc_agenda", "scat_za_grid", "scat_aa_grid",
            "pnd_field",
            "opt_prop_part_agenda", 
            "ppath_step_agenda", "p_grid", "lat_grid", "lon_grid", "z_field",
            "refellipsoid", "t_field", "edensity_field",
            "f_grid", "f_index", "doit_za_interp" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
  
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_scat_fieldCalc" ),
        DESCRIPTION
        (
         "Calculates the scattering integral field in the DOIT module.\n"
         "\n"
         "The scattering integral field is generated by integrating\n"
         "the product of phase matrix and Stokes vector over all incident\n"
         "angles. For more information please refer to AUG.\n" 
         ),
        AUTHORS( "Sreerekha T.R.", "Claudia Emde" ),
        OUT( "doit_scat_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "doit_scat_field", "pha_mat_spt_agenda",
            "doit_i_field", "pnd_field", "t_field", "atmosphere_dim", 
            "cloudbox_limits", "scat_za_grid", "scat_aa_grid",  
            "doit_za_grid_size" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_scat_fieldCalcLimb" ),
        DESCRIPTION
        (
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
         "*DoitAngularGridsSet* and it should always be used for limb\n"
         "simulations.\n"
         "\n"
         "For more information please refer to AUG.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "doit_scat_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "doit_scat_field", "pha_mat_spt_agenda",
            "doit_i_field", "pnd_field", "t_field", "atmosphere_dim", 
            "cloudbox_limits", "scat_za_grid", "scat_aa_grid",  
            "doit_za_grid_size", "doit_za_interp" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "DoitScatteringDataPrepare" ),
        DESCRIPTION
        (
         "Prepares single scattering data for a DOIT scattering calculation.\n"
         "\n"
         "First the scattering data is interpolated in frequency using\n"
         "*scat_data_monoCalc*. Then the phase matrix data is\n"
         "transformed or interpolated from the raw data to the laboratory frame\n"
         "for all possible combinations of the angles contained in the angular\n"
         "grids which are set in *DoitAngularGridsSet*. The resulting phase\n"
         "matrices are stored in *pha_mat_sptDOITOpt*.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "pha_mat_sptDOITOpt", "scat_data_mono" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "doit_za_grid_size", "scat_aa_grid", "scat_data_raw", "f_grid", 
            "f_index", "atmosphere_dim", "stokes_dim" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "DoitWriteIterationFields" ),
        DESCRIPTION
        (
         "Writes DOIT iteration fields.\n"
         "\n"
         "This method writes intermediate iteration fields to xml-files. The\n"
         "method can be used as a part of *doit_conv_test_agenda*.\n"
         "\n"
         "The iterations to be stored are specified by *iterations*, e.g.:\n"
         "    iterations = [3, 6, 9]\n"
         "In this case the 3rd, 6th and 9th iterations are stored in the\n"
         "files 'doit_iteration_3.xml', 'doit_iteration_6.xml' ...\n"
         "If a number is larger than the total number of iterations, this\n" 
         "number is ignored. If all iterations should be stored set\n"
         "   iterations = [0]\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "doit_iteration_counter", "doit_i_field" ),
        GIN( "iterations" ),
        GIN_TYPE( "ArrayOfIndex" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "See above." )
        ));

   md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_za_grid_optCalc" ),
        DESCRIPTION
        (
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
         "frequency.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "doit_za_grid_opt" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "doit_i_field", "scat_za_grid", "doit_za_interp" ),
        GIN( "acc" ),
        GIN_TYPE( "Numeric" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Accuracy to achieve [%]." )
        ));
                                                                               
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_za_interpSet" ),
        DESCRIPTION
        (
         "Define interpolation method for zenith angle dimension.\n"
         "\n"
         "You can use this method to choose the interpolation method for\n"
         "interpolations in the zenith angle dimension.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "doit_za_interp" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim" ),
        GIN( "interp_method" ),
        GIN_TYPE( "String" ),
        GIN_DEFAULT( "linear" ),
        GIN_DESC( "Interpolation method (\"linear\" or \"polynomial\")." )
        ));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Error" ),
        DESCRIPTION
        (
         "Issues an error and exits ARTS.\n"
         "\n"
         "This method can be placed in agendas that must be specified, but\n"
         "are expected not to be used for the particular case. An inclusion\n"
         "in *surface_rtprop_agenda* could look like:\n   "
         "Error{\"Surface interceptions of propagation path not expected.\"}\n"
         "\n"
         "Ignore and other dummy method calls must still be included.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN( "msg" ),
        GIN_TYPE( "String" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "String describing the error." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Exit" ),
        DESCRIPTION
        (
         "Stops the execution and exits ARTS.\n"
         "\n"
         "This method is handy if you want to debug one of your control\n"
         "files. You can insert it anywhere in the control file. When\n"
         "it is reached, it will terminate the program.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Extract" ),
        DESCRIPTION
        (
         "Extract an element from an array.\n"
         "\n"
         "Copies the element with the given Index from the input\n"
         "variable to the output variable.\n"
         "\n"
         "For a Tensor3 as an input, it copies the page with the given\n"
         "Index from the input Tensor3 variable to the output Matrix.\n"
         "\n"
         "In other words, the selection is always done on the first dimension.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT( "needle" ),
        GOUT_TYPE( "ArrayOfIndex, Numeric, Vector, Matrix, Matrix, Tensor3, Tensor4,"
                   "Tensor4, ArrayOfGriddedField3, GriddedField4, String" ),
        GOUT_DESC( "Extracted element." ),
        IN(),
        GIN( "haystack", "index" ),
        GIN_TYPE( "ArrayOfArrayOfIndex, Vector, ArrayOfVector, ArrayOfMatrix, Tensor3,"
                  "Tensor4, ArrayOfTensor4, Tensor5, ArrayOfArrayOfGriddedField3,"
                  "ArrayOfGriddedField4, ArrayOfString",
                  "Index" ),
        GIN_DEFAULT( NODEF, NODEF ),
        GIN_DESC( "Variable to extract from.",
                  "Position of the element which should be extracted." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ext_matAddGas" ),
        DESCRIPTION
        (
         "Add gas absorption to all diagonal elements of extinction matrix.\n"
         "\n"
         "The task of this method is to sum up the gas absorption of the\n"
         "different gas species and add the result to the extinction matrix.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "ext_mat" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "ext_mat", "abs_mat_per_species" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ext_matAddPart" ),
        DESCRIPTION
        (
         "The particle extinction is added to *ext_mat*\n"
         "\n"
         "This function sums up the extinction matrices for all particle\n"
         "types weighted with particle number density.\n"
         "The resulting extinction matrix is added to the workspace\n"
         "variable *ext_mat*\n"
         "The output of this method is *ext_mat* (stokes_dim, stokes_dim).\n"
         "The inputs are the extinction matrix for the single particle type\n"
         "*ext_mat_spt* (N_particletypes, stokes_dim, stokes_dim) and the local\n"
         "particle number densities for all particle types namely the\n"
         "*pnd_field* (N_particletypes, p_grid, lat_grid, lon_grid ) for given\n"
         "*p_grid*, *lat_grid*, and *lon_grid*. The particle types required\n"
         "are specified in the control file.\n"
         ),
        AUTHORS( "Sreerekha T.R." ),
        OUT( "ext_mat" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "ext_mat", "ext_mat_spt", "pnd_field", "atmosphere_dim", 
            "scat_p_index", "scat_lat_index", "scat_lon_index" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ext_matInit" ),
        DESCRIPTION
        (
         "Initialize extinction matrix.\n"
         "\n"
         "This method is necessary, because all other extinction methods just\n"
         "add to the existing extinction matrix.\n"
         "\n"
         "So, here we have to make it the right size and fill it with 0.\n"
         "\n"
         "Note, that the matrix is not really a matrix, because it has a\n"
         "leading frequency dimension.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "ext_mat" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "f_grid", "stokes_dim", "f_index" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "FrequencyFromWavelength" ),
        DESCRIPTION
        (
         "Convert from wavelength [m] to frequency [Hz].\n"
         "\n"
         "This is a generic method. It can take a single wavelength value or a wavelength vector as input.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT(),
        GOUT("frequency"),
        GOUT_TYPE("Numeric, Vector"),
        GOUT_DESC("frequency [Hz]"),
        IN(),
        GIN( "wavelength"),
        GIN_TYPE("Numeric, Vector" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC("wavelength [m]" )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "FlagOff" ),
        DESCRIPTION
        (
         "Sets an index variable that acts as an on/off flag to 0.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT( "i" ),
        GOUT_TYPE( "Index" ),
        GOUT_DESC( "Variable to set to 0." ),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "FlagOn" ),
        DESCRIPTION
        (
         "Sets an index variable that acts as an on/off flag to 1.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT( "i" ),
        GOUT_TYPE( "Index" ),
        GOUT_DESC( "Variable to set to 1." ),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "ForLoop" ),
        DESCRIPTION
        (
         "A simple for loop.\n"
         "\n"
         "This method is handy when you quickly want to test out a calculation\n"
         "with a set of different settings.\n"
         "\n"
         "It does a for loop from start to stop in steps of step. (Who would\n"
         "have guessed that.) For each iteration, the agenda *forloop_agenda* is\n"
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
         "that *ybatchCalc* may occur inside *forloop_agenda*.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "forloop_agenda" ),
        GIN( "start", "stop",  "step" ),
        GIN_TYPE(    "Index", "Index", "Index" ),
        GIN_DEFAULT( NODEF,   NODEF,   NODEF ),
        GIN_DESC( "Start value.",
                  "End value.",
                  "Step size." )
        ));
  /*
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "fos_yStandard" ),
        DESCRIPTION
        (
         "DO NOT USE!\n"
         "\n"
         "This method is used by Patrick for testing a new scattering scheme.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "fos_y", "diy_dx" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "rte_pos", "atmosphere_dim", "p_grid", "lat_grid", "z_field", 
            "t_field", "vmr_field", "edensity_field", "cloudbox_on", 
            "cloudbox_limits", "stokes_dim", "f_grid", "ppath_agenda", 
            "blackbody_radiation_agenda", "abs_scalar_gas_agenda", "iy_main_agenda",
            "iy_transmission", "pnd_field", "scat_data_raw", 
            "fos_y_agenda", "fos_angles",
            "use_mean_scat_data", "fos_n", "fos_i" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
  */

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "f_gridFromGasAbsLookup" ),
        DESCRIPTION
        (
         "Sets *f_grid* to the frequency grid of *abs_lookup*.\n"
         "\n"
         "Must be called between importing/creating raw absorption table and\n"
         "call of *abs_lookupAdapt*.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "f_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_lookup" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "f_gridFromSensorAMSU" ),
        DESCRIPTION
        (
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
         "see *f_gridFromSensorHIRS*.\n"
         ),
        AUTHORS( "Stefan Buehler, Mathias Milz" ),
        OUT( "f_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "lo_multi", "f_backend_multi", "backend_channel_response_multi" ),
        GIN( "spacing" ),
        GIN_TYPE(    "Numeric" ),
        GIN_DEFAULT( ".1e9" ),
        GIN_DESC( "Desired grid spacing in Hz." )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "f_gridFromSensorHIRS" ),
        DESCRIPTION
        (
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
         "*f_gridFromSensorAMSU*.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "f_grid" ),
        GOUT(      ),
        GOUT_TYPE( ),
        GOUT_DESC(),
        IN( "f_backend", "backend_channel_response" ),
        GIN( "spacing" ),
        GIN_TYPE( "Numeric" ),
        GIN_DEFAULT( "5e8" ),
        GIN_DESC( "Desired grid spacing in Hz." )
        ));
  
  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "f_gridSelectFIndex" ),
        DESCRIPTION
        (
         "Reduce f_grid to the frequency given by f_index.\n"
         "\n"
         "This is one of the methods necessary to do line by line absorption\n"
         "calculations inside *abs_mat_per_species_agenda*.\n"
         "\n"
         "It reduces the f_grid to only one frequency, the one given by\n"
         "f_index. If f_index is -1, then all frequencies are kept. This\n"
         "behavior is consistent with *abs_mat_per_speciesExtractFromLookup*.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "f_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "f_grid", "f_index" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "FieldFromGriddedField" ),
        DESCRIPTION
        (
         "FIXME: OLE\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT( "field" ),
        GOUT_TYPE( "Matrix, Tensor3, Tensor4, Tensor4" ),
        GOUT_DESC( "Extracted field." ),
        IN( "p_grid", "lat_grid", "lon_grid" ),
        GIN( "gfield" ),
        GIN_TYPE( "GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Raw input gridded field." )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "g0Earth" ),
        DESCRIPTION
        (
         "Gravity at zero altitude on Earth.\n"
         "\n"
         "Sets *g0* for the given latitude using a standard parameterisation.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "g0" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "lat" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "GriddedFieldLatLonExpand" ),
        DESCRIPTION
        (
         "Expands the latitude and longitude grid of the GriddedField to\n"
         "[-90, 90] and [0,360]. Lat and lon input grids must have size 1.\n"
         "The values from the input data will be duplicated to accomodate\n"
         "for the larger size of the output field.\n"
         "gfield_raw_out and gfield_raw_in can be the same variable.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT( "gfield_raw_out" ),
        GOUT_TYPE( "GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3" ),
        GOUT_DESC( "Expanded gridded field." ),
        IN(),
        GIN( "gfield_raw_in" ),
        GIN_TYPE( "GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Raw input gridded field." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "GriddedFieldPRegrid" ),
        DESCRIPTION
        (
         "Interpolates the input field along the pressure dimension to *p_grid*.\n"
         "gfield_raw_out and gfield_raw_in can be the same variable.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT( "gfield_raw_out" ),
        GOUT_TYPE( "GriddedField3, GriddedField4, ArrayOfGriddedField3" ),
        GOUT_DESC( "Regridded gridded field." ),
        IN( "p_grid" ),
        GIN( "gfield_raw_in", "interp_order", "zeropadding" ),
        GIN_TYPE( "GriddedField3, GriddedField4, ArrayOfGriddedField3",
                  "Index",
                  "Index" ),
        GIN_DEFAULT( NODEF, "1", "0" ),
        GIN_DESC( "Raw input gridded field.",
                  "Interpolation order.",
                  "Apply zero-padding." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "GriddedFieldLatLonRegrid" ),
        DESCRIPTION
        (
         "Interpolates the input field along the latitude and longitude dimensions\n"
         "to *lat_true* and *lon_true*. If the input longitude grid is outside\n"
         "of *lon_true* it will be shifted left or right by 360.\n"
         "If the input longitude grid covers 360 degrees, a cyclic interpolation\n"
         "will be performed.\n"
         "gfield_raw_out and gfield_raw_in can be the same variable.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT( "gfield_raw_out" ),
        GOUT_TYPE( "GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3" ),
        GOUT_DESC( "Regridded gridded field." ),
        IN( "lat_true", "lon_true" ),
        GIN( "gfield_raw_in", "interp_order" ),
        GIN_TYPE( "GriddedField2, GriddedField3, GriddedField4, ArrayOfGriddedField3",
                  "Index" ),
        GIN_DEFAULT( NODEF, "1" ),
        GIN_DESC( "Raw input gridded field.",
                  "Interpolation order." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Massdensity_cleanup" ),
        DESCRIPTION
        (
         "This WSM checks if *massdensity_field* contains values smaller than\n"
         "*massdensity_threshold*. In this case, these values will be set to zero.\n"
         "\n"
         "The Method should be applied if *massdensity_field* contains unrealistic small\n"
         "or erroneous data. (e.g. the chevallierl_91l data sets contain these small values)\n"
         "\n"
         "*Massdensity_cleanup* is called after generation of atmopheric fields.\n"
         "\n"
         "*Default value*:\t1e-15\n"
         ),
        AUTHORS( "Daniel Kreyling" ),
        OUT( "massdensity_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "massdensity_field" ),
        GIN( "massdensity_threshold" ),
        GIN_TYPE( "Numeric" ),
        GIN_DEFAULT( "1e-15" ),
        GIN_DESC( "Values in *massdensity_field* smaller than *massdensity_threshold* will be set to zero." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Ignore" ),
        DESCRIPTION
        (
         "Ignore a workspace variable.\n"
         "\n"
         "This method is handy for use in agendas in order to suppress warnings\n"
         "about unused input workspace variables. What it does is: Nothing!\n"
         "In other words, it just ignores the variable it is called on.\n"
         "\n"
         "This method can ignore any workspace variable\n"
         "you want.\n"
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
         "'elsLorentz' does not need it.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(       "gin1"    ),
        GIN_TYPE(     "Any" ),
        GIN_DEFAULT(  NODEF ),
        GIN_DESC( "Variable to be ignored." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "INCLUDE" ),
        DESCRIPTION
        (
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
         "controlfiles.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "IndexSet" ),
        DESCRIPTION
        (
         "Sets an index workspace variable to the given value.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "i"     ),
        GOUT_TYPE( "Index" ),
        GOUT_DESC( "Variable to initialize." ),
        IN(),
        GIN(       "value" ),
        GIN_TYPE(  "Index" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Value." ),
        SETMETHOD( true )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "IndexStepDown" ),
        DESCRIPTION
        (
         "Performas: i2 = i1 - 1\n"
         "\n"
         "Input and output can be same variable.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "i2"    ),
        GOUT_TYPE( "Index" ),
        GOUT_DESC( "Output index variable." ),
        IN(),
        GIN(       "i1"  ),
        GIN_TYPE(  "Index" ),
        GIN_DEFAULT( NODEF   ),
        GIN_DESC( "Input index variable." )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "IndexStepUp" ),
        DESCRIPTION
        (
         "Performas: i2 = i1 + 1\n"
         "\n"
         "Input and output can be same variable.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "i2"    ),
        GOUT_TYPE( "Index" ),
        GOUT_DESC( "Output index variable." ),
        IN(),
        GIN(       "i1"  ),
        GIN_TYPE(  "Index" ),
        GIN_DEFAULT( NODEF   ),
        GIN_DESC( "Input index variable." )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "InterpAtmFieldToRtePos" ),
        DESCRIPTION
        (
         "Point interpolation of atmospheric fields.\n" 
         "\n"
         "The standard choice for *pos* should be *rte_pos*.\n"
         "\n"
         "Linear interpolation is applied.\n"         
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "x"       ),
        GOUT_TYPE( "Numeric" ),
        GOUT_DESC( "Value obtained by the interpolation." ),
        IN( "atmosphere_dim", "p_grid", "lat_grid", "lon_grid", "z_field" ),
        GIN(      "pos",    "field"   ),
        GIN_TYPE( "Vector", "Tensor3" ),
        GIN_DEFAULT( NODEF, NODEF     ),
        GIN_DESC( "Position vector.", "Field to interpolate." )
        ));
  
  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "InterpSurfaceFieldToRtePos" ),
        DESCRIPTION
        (
         "Point interpolation of surface fields.\n" 
         "\n"
         "The standard choice for *pos* should be *rte_pos*.\n"
         "\n"
         "Linear interpolation is applied.\n"         
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "x"       ),
        GOUT_TYPE( "Numeric" ),
        GOUT_DESC( "Value obtained by interpolation." ),
        IN( "atmosphere_dim", "lat_grid", "lon_grid" ),
        GIN(      "pos",    "field"   ),
        GIN_TYPE( "Vector", "Matrix" ),
        GIN_DEFAULT( NODEF, NODEF ),
        GIN_DESC( "Position vector.", "Field to interpolate." )
        ));
  
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "iyApplyYunit" ),
        DESCRIPTION
        (
         "Conversion of *iy* to other spectral units.\n"
         "\n"
         "No automatic unit conversion is applied on the iy-variables, but\n"
         "this method allows for a seperate conversion of *iy* and *iy_aux*\n"
         "to the unit selected by *y_unit*.\n"
         "\n"
         "Only these auxilary quantities are modified:\n"
         "    \"iy\", \"Error\" and \"Error (uncorrelated)\"\n"
         "\n"
         "Please note that *diy*dx* is not handled and it is demanded that\n"
         "*jacobian_do* is 0.\n"
         "\n"
         "See further *y_unit*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "iy", "iy_aux" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "iy", "iy_aux", "stokes_dim", "f_grid", "jacobian_do", 
            "iy_aux_vars", "y_unit" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord

      ( NAME( "iyCalc" ),
        DESCRIPTION
        (
         "A single monochromatic pencil beam calculation.\n"
         "\n"
         "Performs monochromatic radiative transfer calculations for the\n"
         "specified position (*rte_pos*) and line-of-sight (*rte_pos*).\n"
         "See *iy* and associated variables for format of output.\n"
         "\n"
         "Only analytical Jacobian calculations are performed. Hence,\n"
         "elements of *diy_dx* can be empty.\n"
         "\n"
         "No unit conversion is applied (but can be done as post-processing).\n"
         "\n"
         "No sensor charactersitcs are applied. That is most easily\n"
         "incorporated by using *yCalc*\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "iy", "iy_aux", "ppath", "diy_dx" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "basics_checked", "iy_aux_vars", "t_field", 
            "z_field", "vmr_field", "cloudbox_on", "cloudbox_checked", 
            "rte_pos", "rte_los", "jacobian_do", "mblock_index", 
            "iy_main_agenda" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "iyEmissionStandard" ),
        DESCRIPTION
        (
         "Standard method for radiative transfer calculations with emission.\n"
         "\n"
         "Designed to be part of *iy_main_agenda*. That is, only valid\n"
         "outside the cloudbox (no scattering or polarised absorption).\n"
         "Assumes local thermodynamic equilibrium for emission.\n"
         "\n"
         "The overall strategy is to take the average of the absorption and\n"
         "the emission source function at the end points of each step of\n"
         "the propagation path. See further the user guide.\n" 
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
         "* \"Radiative background\": Index value flagging the radiative\n"
         "     background. The following coding is used: 0=space, 1=surface\n"
         "     and 2=cloudbox. Size: [nf,1,1,1].\n"
         "  \"iy\": The radiance at each point along the path.\n"
         "     Size: [nf,ns,1,np].\n"
         "* \"Optical depth\": The scalar optical depth between the\n"
         "     observation point and the end of the primary propagation path\n"
         "     (ie. the optical depth to the surface, cloudbox or space.)\n"
         "     Size: [nf,1,1,1].\n"
         "where\n"
         "  nf: Number of freqiencies.\n"
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
        IN( "stokes_dim", "f_grid", "atmosphere_dim", "p_grid", "z_field",
            "t_field", "vmr_field", "abs_species", 
            "wind_u_field", "wind_v_field", "wind_w_field", "mag_u_field",
            "mag_v_field", "mag_w_field", "edensity_field",
            "cloudbox_on", "iy_aux_vars", "jacobian_do", 
            "jacobian_quantities", "jacobian_indices", 
            "ppath_agenda", "blackbody_radiation_agenda",
            "abs_mat_per_species_agenda", "iy_main_agenda", 
            "iy_space_agenda", "iy_surface_agenda", "iy_cloudbox_agenda", 
            "iy_agenda_call1", "iy_transmission", "mblock_index", 
            "rte_pos", "rte_los" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  /*
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "iyFOS" ),
        DESCRIPTION
        (
         "DO NOT USE! So far just used for testing by Patrick.\n"
         "\n"
         "A fixed order of scattering (FOS) scheme.\n"
         "\n"
         "The scattering integral is here solved by calculating the incoming\n"
         "radiation and phase matrix for a set of directions, that can be\n"
         "distributed freely over the sphere of integration. These directions\n"
         "are specified by a combination of zenith angd azimuth angles for\n"
         "each integration point. The product of incoming radiation and phase\n"
         "matrix is summed up, with a weight specified for each product.\n"
         "These angles and integration weights are packed into the WSV\n" 
         "*fos_angles*.\n"
         "\n"
         "If *fos_n* equals 1, the incoming radiation is calculated without\n"
         "scattering. This is thus a singel scattering scheme.\n"
         "For *fos_n* = n, the incoming radiation for a point at the\n"
         "unscattered propgation path is calculated as for fos_n = n-1. This\n"
         "gives a scheme where n scattering events are considered.\n"
         "\n"
         "For regions of no scattering there is a marginal difference between\n"
         "running this method and the corresponding pure clear-sky agenda,\n"
         "and there is not critical to set the cloudbox limits in a \"tight\"\n"
         "manner.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "iy", "diy_dx" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "diy_dx", "iy_transmission",
            "rte_pos", "rte_los", "jacobian_do","atmosphere_dim", "p_grid", 
            "z_field", "t_field", "vmr_field", "edensity_field", "cloudbox_on", 
            "cloudbox_limits", "stokes_dim", "f_grid", "ppath_agenda", 
            "blackbody_radiation_agenda", "abs_mat_per_species_agenda", "iy_main_agenda",
            "pnd_field", "scat_data_raw", "opt_prop_gas_agenda", "fos_y_agenda",
            "fos_angles", "use_mean_scat_data", "fos_n", "fos_i" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
  */

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "iyMC" ),
        DESCRIPTION
        (
         "Interface to Monte Carlo part for *iy_main_agenda*.\n"
         "\n"
         "Basically an interface to *MCGeneral* for doing monochromatic\n"
         "pencil beam calculations. This functions allows Monte Carlo (MC)\n"
         "calculations for sets of frequencies and sensor pos/los in a single\n"
         "run. Sensor responses can be included in the standard manner\n" 
         "(through *yCalc*).\n"
         "\n"
         "MC unit is set as for *MCGeneral*. No antenna pattern is included.\n"
         "\n"
         "This function does not apply the MC approach when it comes\n"
         "to sensor properties. These properties are not considered when\n"
         "tracking photons, which is done in *MCGeneral* (but then only for\n"
         "the antenna pattern).\n"
         "\n"
         "The MC calculation errors are all assumed be uncorrelated and each\n"
         "have a normal distribution. These properties are of relevance when\n"
         "weighting the errors with the sensor repsonse matrix. The seed is\n"
         "reset for each call of *MCGeneral* to obtain uncorrelated errors.\n"
         "\n"
         "MC control arguments (mc_std_err, mc_max_time, mc_max_iter and\n"
         "mc_z_field_is_1D) as for *MCGeneral*. The arguments are applied\n"
         "for each monochromatic pencil beam calculation individually.\n"
         "As or *MCGeneral*, the value of *mc_error* shall be adopted to\n"
         "*y_unit*. However, to be consistent with other iy_methods, this\n"
         "method returns radiance data (*y_unit*) is applied inside yCalc).\n"
         "\n"
         "The following auxiliary data can be obtained:\n"
         "  \"Error (uncorrelated)\": Calculation error. Size: [nf,ns,1,1].\n"
         "    (The later part of the text string is required. It is used as\n"
         "    a flag to yCalc for how to apply the sensor data.)\n"
         "where\n"
         "  nf: Number of freqiencies.\n"
         "  ns: Number of Stokes elements.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "iy", "iy_aux", "diy_dx" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "iy_agenda_call1", "iy_transmission", 
            "rte_pos", "rte_los", "iy_aux_vars",  
            "jacobian_do", "atmosphere_dim", "p_grid", "lat_grid",
            "lon_grid", "z_field", "t_field", "vmr_field", "refellipsoid", 
            "z_surface", "cloudbox_on", "cloudbox_limits", "cloudbox_checked",
            "stokes_dim", "f_grid", "scat_data_raw", 
            "iy_space_agenda", "surface_rtprop_agenda", "abs_mat_per_species_agenda",
            "pnd_field", "y_unit",
            "mc_std_err", "mc_max_time", "mc_max_iter"),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "iyInterpCloudboxField" ),
        DESCRIPTION
        (
         "Interpolates the intensity field of the cloud box.\n"
         "\n"
         "This is the standard method to put in *iy_cloudbox_agenda* if the\n"
         "the scattering inside the cloud box is handled by the DOIT method.\n"
         "\n"
         "The intensity field is interpolated to the position (specified by\n"
         "*rte_pos*) and direction (specified by *rte_los*) given. A linear\n"
         "interpolation is used for all dimensions.\n"
         "\n"
         "The intensity field on the cloux box boundaries is provided by\n"
         "*scat_i_p/lat/lon* and these variables are interpolated if the\n"
         "given position is at any boundary.\n"
         "\n"
         "Interpolation of the internal field is not yet possible.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "iy" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "scat_i_p", "scat_i_lat", "scat_i_lon", "doit_i_field1D_spectrum",
            "rte_pos", "rte_los", "jacobian_do","cloudbox_on", 
            "cloudbox_limits", "atmosphere_dim", "p_grid", "lat_grid",
            "lon_grid", "z_field", "stokes_dim", 
            "scat_za_grid", "scat_aa_grid", "f_grid" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "iyInterpPolyCloudboxField" ),
        DESCRIPTION
        (
         "As *iyInterpCloudboxField* but performs cubic interpolation.\n"
         "\n"
         "Works so far only for 1D cases, and accordingly a cubic\n"
         "interpolation along *scat_za_grid* is performed.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "iy" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "scat_i_p", "scat_i_lat", "scat_i_lon", "doit_i_field1D_spectrum",
            "rte_pos", "rte_los", "jacobian_do", "cloudbox_on", 
            "cloudbox_limits", "atmosphere_dim", "p_grid", "lat_grid",
            "lon_grid", "z_field", "stokes_dim", "scat_za_grid", 
            "scat_aa_grid", "f_grid" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "iyRadioLink" ),
        DESCRIPTION
        (
         "Method in development ...\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "iy", "iy_aux", "ppath", "diy_dx" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "iy_agenda_call1", "iy_transmission", "rte_pos", "iy_aux_vars", 
            "jacobian_do", "atmosphere_dim", 
            "p_grid", "lat_grid", "lon_grid", "z_field", "t_field", "vmr_field",
            "wind_u_field", "wind_v_field", "wind_w_field", "mag_u_field",
            "mag_v_field", "mag_w_field", "edensity_field", 
            "refellipsoid", "z_surface", "cloudbox_on", "stokes_dim", "f_grid",
            "dispersion_do", "mblock_index", "ppath_agenda", 
            "ppath_step_agenda", "abs_mat_per_species_agenda", "iy_transmitter_agenda" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "iySurfaceRtpropAgenda" ),
        DESCRIPTION
        (
         "Interface to *surface_rtprop_agenda* for *iy_surface_agenda*.\n"
         "\n"
         "This method is designed to be part of *iy_surface_agenda*. It\n"
         "determines the radiative properties of the surface by\n"
         "*surface_rtprop_agenda* and calculates the downwelling radiation\n"
         "by *iy_main_agenda*, and sums up the terms as described in AUG.\n"
         "That is, this WSM uses the output from *surface_rtprop_agenda*\n"
         "in a straightforward fashion.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "iy", "diy_dx" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "iy_transmission", "jacobian_do", "atmosphere_dim", "t_field", 
            "z_field", "vmr_field", "cloudbox_on", "stokes_dim", "f_grid", 
            "rte_pos", "rte_los", "iy_main_agenda", "surface_rtprop_agenda"
          ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "iyTransmissionStandard" ),
        DESCRIPTION
        (
         "Standard method for handling (direct) transmission measurements.\n"
         "\n"
         "Designed to be part of *iy_main_agenda*. For cases with no\n"
         "cloudbox, this method works basically as *iyEmissionStandard*,\n"
         "beside that no emission is generated.\n"
         "\n"         
         "For propagation paths above the surface, the result is the output\n"
         "of *iy_space_agenda* times the transmission through the atmosphere.\n"
         "That is, a standard Beer-Lambert calculation.\n"
         "\n"
         "If the propagation paths intercepts with the surface, the surface\n"
         "variables are used as for *iyEmissionStandard*, and \n"
         "*surface_emission* should probably be set to 0. Multiple-beams\n"
         "will be followed if *surface_los* has length > 1.\n"
         "\n"
         "The following auxiliary data can be obtained:\n"
         "  \"Temperature\": The temperature along the propagation path.\n"
         "  \"VMR, species X\": VMR of the species with index X (zero based).\n"
         "     For example, adding the string \"VMR, species 0\" extracts the\n"
         "     VMR of the first species.\n"
         "  \"Absorption, summed\": The total absorption along the path.\n"
         "     If stokes_dim>1, the absorption vector is stored.\n"
         "  \"Absorption, species X\": The absorption vector along the path\n"
         "     for an individual species (X works as for VMR).\n"
         "  \"Radiance\": The state of *iy* at each point along the path.\n"
         "\n"
         "The auxiliary data are returned in *iy_aux* with quantities\n"
         "selected by *iy_aux_vars*. All variables require that the method\n"
         "is called directly or by *iyCalc*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "iy", "iy_aux", "ppath", "diy_dx" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "stokes_dim", "f_grid", "atmosphere_dim", "p_grid", "z_field", 
            "t_field", "vmr_field", "abs_species", 
            "wind_u_field", "wind_v_field", "wind_w_field", "mag_u_field",
            "mag_v_field", "mag_w_field", "edensity_field",
            "cloudbox_on", "cloudbox_limits", "iy_aux_vars", 
            "jacobian_do", "jacobian_quantities", "jacobian_indices", 
            "ppath_agenda", "abs_mat_per_species_agenda", "iy_transmitter_agenda",
            "iy_agenda_call1", "iy_transmission", "mblock_index", 
            "rte_pos", "rte_los" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianAddAbsSpecies" ),
        DESCRIPTION
        (
         "Includes an absorption species in the Jacobian.\n"
         "\n"
         "Details are given in the user guide.\n"
         "\n"         
         "For 1D or 2D calculations the latitude and/or longitude grid of\n"
         "the retrieval field should set to have zero length.\n"
         "\n"
         "There are two possible calculation methods:\n"
         "   \"analytical\"   : (semi-)analytical expressions are used\n"
         "   \"perturbation\" : pure numerical perturbations are used\n"
         "\n"
         "The retrieval unit can be:\n"
         "   \"vmr\"    : Volume mixing ratio.\n"
         "   \"nd\"     : Number density.\n"
         "   \"rel\"    : Relative unit (e.g. 1.1 means 10% more of the gas).\n"
         "   \"logrel\" : The retrieval is performed with the logarithm of\n"
         "                the \"rel\" option.\n"
         "\n"
         "For perturbation calculations the size of the perturbation is set\n"
         "by the user. The unit for the perturbation is the same as for the\n"
         "retrieval unit.\n"
         ),
        AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
        OUT( "jacobian_quantities", "jacobian_agenda" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "jacobian_quantities", "jacobian_agenda",
            "atmosphere_dim", "p_grid", "lat_grid", "lon_grid" ),
        GIN( "g1", "g2", "g3", "species", "method", "unit","dx" ),
        GIN_TYPE( "Vector", "Vector", "Vector", "String", "String", "String", 
                  "Numeric" ),
        GIN_DEFAULT( NODEF, NODEF, NODEF, NODEF, "analytical", "rel", "0.001" ),
        GIN_DESC( "Pressure retrieval grid.",
                  "Latitude retrieval grid.",
                  "Longitude retreival grid.",
                  "The species tag of the retrieval quantity.",
                  "Calculation method. See above.",
                  "Retrieval unit. See above.",
                  "Size of perturbation." 
                  ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( false ),
        PASSWORKSPACE(  true  )
        ));
         
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianAddFreqShiftAndStretch" ),
        DESCRIPTION
        (
         "Includes a frequency fit in the Jacobian.\n"
         "\n"
         "Retrieval of deviations between nominal and actual backend\n"
         "frequencies can be included by this method. The calculations can be\n"
         "performed in the following ways:\n"
         "   calcmode = \"interp\": Interpolation of monochromatic spectra,\n"
         "      shifted with *df* from nominal values.\n"
         "\n"
         "The frequencies can be fitted with 1 or 2 variables. The first one\n"
         "is a \"shift\". That is, an off-set common for all backend channels.\n"
         "The second variable is \"stretch\", that is included only if\n"
         "*do_stretch* is set to 1. The stretch is a frequency shift that goes\n"
         "from -x at the first channel and increases linearly to reach x at\n"
         "the last channel, where x is the retrieved value.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "jacobian_quantities", "jacobian_agenda" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "jacobian_quantities", "jacobian_agenda", "f_grid" ),
        GIN( "calcmode", "df", "do_stretch" ),
        GIN_TYPE( "String", "Numeric", "Index" ),
        GIN_DEFAULT( "interp", "100e3", "0" ),
        GIN_DESC( "Calculation method. See above",
                  "Size of perturbation to apply.", 
                  "Flag to also include frequency stretch."
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianAddPointingZa" ),
        DESCRIPTION
        (
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
         ),
        AUTHORS( "Patrick Eriksson", "Mattias Ekstrom" ),
        OUT( "jacobian_quantities", "jacobian_agenda" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "jacobian_quantities", "jacobian_agenda", "sensor_pos", 
            "sensor_time" ),
        GIN( "poly_order", "calcmode", "dza" ),
        GIN_TYPE( "Index", "String", "Numeric" ),
        GIN_DEFAULT( "0", "recalc", "0.01" ),
        GIN_DESC( "Order of polynomial to describe the time variation of "
                  "pointing off-sets.",
                  "Calculation method. See above",
                  "Size of perturbation to apply (when applicable)."
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianAddPolyfit" ),
        DESCRIPTION
        (
         "Includes polynomial baseline fit in the Jacobian.\n"
         "\n"
         "This method deals with retrieval of disturbances of the spectra\n"
         "that can be described by an addidative term, a baseline off-set.\n"
         "\n"
         "The baseline off-set is here modelled as a polynomial. The\n"
         "polynomial spans the complete frequency range spanned by\n"
         "*sensor_response_f_grid* and the method should only of interest for\n"
         "cases with no frequency gap in the spectra. The default assumption\n"
         "is that the off-set differs between all spectra, but it can also be\n"
         "assumed that the off-set is common for all e.g. line-of-sights.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "jacobian_quantities", "jacobian_agenda" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "jacobian_quantities", "jacobian_agenda", 
            "sensor_response_pol_grid", "sensor_response_f_grid",
            "sensor_response_za_grid", "sensor_pos" ),
        GIN( "poly_order", "no_pol_variation", "no_los_variation", 
             "no_mblock_variation" ),
        GIN_TYPE( "Index", "Index", "Index", "Index" ),
        GIN_DEFAULT( NODEF, "0", "0", "0" ),
        GIN_DESC( "Polynomial order to use for the fit.", 
                  "Set to 1 if the baseline off-set is the same for all "
                  "Stokes components.", 
                  "Set to 1 if the baseline off-set is the same for all "
                  "line-of-sights (inside each measurement block).", 
                  "Set to 1 if the baseline off-set is the same for all "
                  "measurement blocks." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianAddTemperature" ),
        DESCRIPTION
        (
         "Includes atmospheric temperatures in the Jacobian.\n"
         "\n"
         "The calculations can be performed by (semi-)analytical expressions\n"
         "or by perturbations. Hydrostatic equilibrium (HSE) can be included.\n"
         "For perturbation calculations, all possible effects are included\n"
         "(but is a costly option). The analytical calculation approach\n"
         "neglects refraction totally, but considers the local effect of HSE.\n"
         "The later should be accaptable for observations around zenith and\n"
         "nadir. There is no warning if the method is applied incorrectly, \n"
         "with respect to these issues.\n"
         "\n"
         "The calculations (both options) assume that gas species are defined\n"
         "in VMR (a change in temperature then changes the number density). \n"
         "This has the consequence that retrieval of temperatures and number\n" 
         "density can not be mixed. Neither any warning here!\n"
         "\n"
         "The choices for *method* are:\n"
         "   \"analytical\"   : (semi-)analytical expressions are used\n"
         "   \"perturbation\" : pure numerical perturbations are used\n"
         ),
        AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
        OUT( "jacobian_quantities", "jacobian_agenda" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "jacobian_quantities", "jacobian_agenda", 
            "atmosphere_dim", "p_grid", "lat_grid", "lon_grid" ),
        GIN( "g1", "g2", "g3", "hse", "method", "dt" ),
        GIN_TYPE( "Vector", "Vector", "Vector", "String", "String", "Numeric" ),
        GIN_DEFAULT( NODEF, NODEF, NODEF, "off", "analytical", "1" ),
        GIN_DESC( "Pressure retrieval grid.",
                  "Latitude retrieval grid.",
                  "Longitude retreival grid.",
                  "Flag to assume HSE or not (\"on\" or \"off\").",
                  "Calculation method. See above.",
                  "Size of perturbation [K]." 
                  )
        ));

  /*
    md_data_raw.push_back
    ( MdRecord
    ( NAME( "jacobianAddParticle" ),
    DESCRIPTION
    (
    "Add particle number density as retrieval quantity to the Jacobian.\n"
    "\n"
    "The Jacobian is done by perturbation calculation by adding elements\n"
    "of *pnd_field_perturb* to *pnd_field*. Only 1D and 3D atmospheres\n"
    "can be handled by this method.\n"
    "\n"
    "The perturbation field and the unit of it are defined outside ARTS.\n"
    "This method only returns the difference between the reference and\n"
    "perturbed spectra. The division by the size of the perturbation\n"
    "also has to be done outside ARTS.\n"
    "The unit of the particle jacobian is the same as for *y*.\n"
    "\n"
    "Generic input:\n"
    "  Vector : The pressure grid of the retrieval field.\n"
    "  Vector : The latitude grid of the retrieval field.\n"
    "  Vector : The longitude grid of the retrieval field.\n"
    ),
    AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
    OUT( "jacobian_quantities", "jacobian_agenda" ),
    GOUT(),
    GOUT_TYPE(),
    GOUT_DESC(),
    IN( "jacobian", "atmosphere_dim", "p_grid", "lat_grid", "lon_grid", 
    "pnd_field", "pnd_field_perturb", "cloudbox_limits" ),
    GIN(      "gin1"      , "gin2"      , "gin3"       ),
    GIN_TYPE(    "Vector", "Vector", "Vector" ),
    GIN_DEFAULT( NODEF   , NODEF   , NODEF    ),
    GIN_DESC( "FIXME DOC",
    "FIXME DOC",
    "FIXME DOC" )
    ));
  */
  
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianCalcAbsSpeciesAnalytical" ),
        DESCRIPTION
        (
         "This function doesn't do anything. It just exists to satisfy\n"
         "the input and output requirement of the *jacobian_agenda*.\n"
         "\n"
         "This function is added to *jacobian_agenda* by\n"
         "jacobianAddAbsSpecies and should normally not be called\n"
         "by the user.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "jacobian" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "jacobian",
            "mblock_index", "iyb", "yb" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianCalcAbsSpeciesPerturbations" ),
        DESCRIPTION
        (
         "Calculates absorption species jacobians by perturbations.\n"
         "\n"
         "This function is added to *jacobian_agenda* by\n"
         "jacobianAddAbsSpecies and should normally not be called\n"
         "by the user.\n"
         ),
        AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
        OUT( "jacobian" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "jacobian",
            "mblock_index", "iyb", "yb", "atmosphere_dim", "p_grid", "lat_grid",
            "lon_grid", "t_field", "z_field", "vmr_field", "abs_species", 
            "cloudbox_on", "stokes_dim", 
            "f_grid", "sensor_pos", "sensor_los", "mblock_za_grid", 
            "mblock_aa_grid", "antenna_dim", "sensor_response",
            "iy_main_agenda", "y_unit", "jacobian_quantities",
            "jacobian_indices" ),
        GIN( "species" ),
        GIN_TYPE(    "String" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Species of interest." ),
        SETMETHOD( true )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianCalcFreqShiftAndStretchInterp" ),
        DESCRIPTION
        (
         "Calculates frequency shift and stretch jacobians by interpolation\n"
         "of *iyb*.\n"
         "\n"
         "This function is added to *jacobian_agenda* by\n"
         "jacobianAddFreqShiftAndStretch and should normally\n"
         "not be called by the user.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "jacobian" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "jacobian",
            "mblock_index", "iyb", "yb", "stokes_dim", "f_grid", 
            "mblock_za_grid", "mblock_aa_grid", "antenna_dim", 
            "sensor_response", "sensor_response_pol_grid", 
            "sensor_response_f_grid", "sensor_response_za_grid", 
            "jacobian_quantities", "jacobian_indices" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianCalcPointingZaInterp" ),
        DESCRIPTION
        (
         "Calculates zenith angle pointing deviation jacobians by\n"
         "inter-extrapolation of *iyb*.\n"
         "\n"
         "This function is added to *jacobian_agenda* by\n"
         "jacobianAddPointingZa and should normally not be\n"
         "called by the user.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "jacobian" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "jacobian", "mblock_index", "iyb", "yb", "stokes_dim", "f_grid", 
            "sensor_los", "mblock_za_grid", "mblock_aa_grid", "antenna_dim", 
            "sensor_response", "sensor_time", 
            "jacobian_quantities", "jacobian_indices" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianCalcPointingZaRecalc" ),
        DESCRIPTION
        (
         "Calculates zenith angle pointing deviation jacobians by\n"
         "recalulation of *iyb*.\n"
         "\n"
         "This function is added to *jacobian_agenda* by\n"
         "jacobianAddPointingZa and should normally not be\n"
         "called by the user.\n"
         ),
        AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
        OUT( "jacobian" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "jacobian",
            "mblock_index", "iyb", "yb", "atmosphere_dim",
            "t_field", "z_field", "vmr_field", "cloudbox_on", 
            "stokes_dim", "f_grid", "sensor_pos", "sensor_los", 
            "mblock_za_grid", "mblock_aa_grid", "antenna_dim", 
            "sensor_response", "sensor_time", 
            "iy_main_agenda", "y_unit", "jacobian_quantities",
            "jacobian_indices" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianCalcPolyfit" ),
        DESCRIPTION
        (
         "Calculates jacobians for polynomial baseline fit.\n"
         "\n"
         "This function is added to *jacobian_agenda* by jacobianAddPolyfit\n"
         "and should normally not be called by the user.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "jacobian" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "jacobian", "mblock_index", "iyb", "yb", "sensor_response",
            "sensor_response_pol_grid", "sensor_response_f_grid", 
            "sensor_response_za_grid", 
            "jacobian_quantities", "jacobian_indices" ),
        GIN( "poly_coeff" ),
        GIN_TYPE( "Index" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Polynomial coefficient to handle." ),
        SETMETHOD( true )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianCalcTemperatureAnalytical" ),
        DESCRIPTION
        (
         "This function doesn't do anything. It just exists to satisfy\n"
         "the input and output requirement of the *jacobian_agenda*.\n"
         "\n"
         "This function is added to *jacobian_agenda* by\n"
         "jacobianAddTemperature and should normally not be called\n"
         "by the user.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "jacobian" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "jacobian",
            "mblock_index", "iyb", "yb" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianCalcTemperaturePerturbations" ),
        DESCRIPTION
        (
         "Calculates atmospheric temperature jacobians by perturbations.\n"
         "\n"
         "This function is added to *jacobian_agenda* by\n"
         "jacobianAddTemperature and should normally not be called\n"
         "by the user.\n"
         ),
        AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
        OUT( "jacobian" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "jacobian",
            "mblock_index", "iyb", "yb", "atmosphere_dim", "p_grid", "lat_grid",
            "lon_grid", "lat_true", "lon_true", "t_field", "z_field", 
            "vmr_field", "abs_species", "refellipsoid", "z_surface", 
            "cloudbox_on", "stokes_dim", "f_grid", "sensor_pos", "sensor_los", 
            "mblock_za_grid", "mblock_aa_grid", "antenna_dim", 
            "sensor_response", "iy_main_agenda", "y_unit", "g0_agenda", 
            "molarmass_dry_air", "p_hse", "z_hse_accuracy", 
            "jacobian_quantities", "jacobian_indices" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));


 /*        
           md_data_raw.push_back
            ( MdRecord
            ( NAME( "jacobianCalcParticle" ),
            DESCRIPTION
            (
            "Calculates particle number densities jacobians by perturbations\n"
            "\n"
            "This function is added to *jacobian_agenda* by jacobianAddParticle\n"
            "and should normally not be called by the user.\n"
            ),
            AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
            OUT( "jacobian" ),
            GOUT(),
            GOUT_TYPE(),
            GOUT_DESC(),
            IN( "y", "jacobian_quantities", "jacobian_indices", "pnd_field_perturb",
            "jacobian_particle_update_agenda",
            "ppath_step_agenda", "rte_agenda", "iy_space_agenda", 
            "surface_rtprop_agenda", "iy_cloudbox_agenda", "atmosphere_dim", 
            "p_grid", "lat_grid", "lon_grid", "z_field", "t_field", "vmr_field",
            "refellipsoid", "z_surface", 
            "cloudbox_on", "cloudbox_limits", "pnd_field",
            "sensor_response", "sensor_pos", "sensor_los", "f_grid", 
            "stokes_dim", "antenna_dim", "mblock_za_grid", "mblock_aa_grid" ),
            GIN(),
            GIN_TYPE(),
            GIN_DEFAULT(),
            GIN_DESC()
            ));
        
 */

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianClose" ),
        DESCRIPTION
        (
         "Closes the array of retrieval quantities and prepares for\n" 
         "calculation of the Jacobian matrix.\n"
         "\n"
         "This function closes the *jacobian_quantities* array, sets the\n"
         "correct size of *jacobian* and sets *jacobian_do* to 1.\n"
         "\n"
         "Retrieval quantities should not be added after a call to this WSM.\n"
         "No calculations are performed here.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUT( "jacobian_do", "jacobian", "jacobian_indices", "jacobian_agenda" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "jacobian_agenda", "jacobian_quantities", "sensor_pos", "sensor_response" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianInit" ),
        DESCRIPTION
        (
         "Initialises the variables connected to the Jacobian matrix.\n"
         "\n"
         "This function initialises the *jacobian_quantities* array so\n"
         "that retrieval quantities can be added to it. Accordingly, it has\n"
         "to be called before any calls to jacobianAddTemperature or\n"
         "similar methods.\n"
         "\n"
         "The Jacobian quantities are initialised to be empty.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUT( "jacobian_quantities", "jacobian_agenda" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianOff" ),
        DESCRIPTION
        (
         "Makes mandatory initialisation of some jacobian variables.\n"
         "\n"
         "Some jacobian WSVs must be initilised even if no such calculations\n"
         "will be performed and this is handled with this method. That is,\n"
         "this method must be called when no jacobians will be calculated.\n"
         "Sets *jacobian_on* to 0.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "jacobian_do", "jacobian_agenda", "jacobian_quantities", 
             "jacobian_indices" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "MatrixCBR" ),
        DESCRIPTION
        (
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
         "followed and the unit of the returned data is W/(m3 * Hz * sr).\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "m"       ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( "Variable to initialize." ),
        IN( "stokes_dim" ),
        GIN(         "f"      ),
        GIN_TYPE(    "Vector" ),
        GIN_DEFAULT( NODEF    ),
        GIN_DESC( "Frequency vector." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MatrixMatrixMultiply" ),
        DESCRIPTION
        (
         "Multiply a Matrix with another Matrix and store the result in the result\n"
         "Matrix.\n"
         "\n"
         "This just computes the normal Matrix-Matrix product, Y=M*X. It is ok\n"
         "if Y and X are the same Matrix. This function is handy for\n"
         "multiplying the H Matrix to batch spectra.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT(),
        GOUT(      "gout1"       ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( "The result of the multiplication (dimension mxc)." ),
        IN(),
        GIN(      "gin1"      , "gin2"       ),
        GIN_TYPE(    "Matrix", "Matrix" ),
        GIN_DEFAULT( NODEF   , NODEF    ),
        GIN_DESC( "The Matrix to multiply (dimension mxn).",
                  "The original Matrix (dimension nxc)." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MatrixPlanck" ),
        DESCRIPTION
        (
         "Sets a matrix to hold blackbody radiation.\n"
         "\n"
         "The radiation is assumed to be un-polarized and Stokes components\n"
         "2-4 are zero. Number of Stokes components, that equals the number\n"
         "of columns in the created matrix, is determined by *stokes_dim*.\n"
         "The number of rows in the created matrix equals the length of the\n"
         "given frequency vector.\n"
         "\n"
         "The standard definition, in ARTS, of the Planck function is\n"
         "followed and the unit of the returned data is W/(m3 * Hz * sr).\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "m"      ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( "Variable to initialize." ),
        IN( "stokes_dim" ),
        GIN(        "gin1"  , "gin2"    ),
        GIN_TYPE(   "Vector", "Numeric" ),
        GIN_DEFAULT( NODEF   , NODEF    ),
        GIN_DESC( "Frequency vector.",
                  "Temperature [K]." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MatrixScale" ),
        DESCRIPTION
        (
         "Scales all elements of a matrix with the specified value.\n"
         "\n"
         "The result can either be stored in the same or another\n"
         "variable.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "mout"   ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( "Output matrix" ),
        IN(),
        GIN(         "min"   , "value"   ),
        GIN_TYPE(    "Matrix", "Numeric" ),
        GIN_DEFAULT( NODEF   , NODEF     ),
        GIN_DESC( "Input matrix.",
                  "The value to be multiplied with the matrix." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MatrixSetConstant" ),
        DESCRIPTION
        (
         "Creates a matrix and sets all elements to the specified value.\n"
         "\n"
         "The size is determined by *ncols* and *nrows*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "m"      ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( "Variable to initialize." ),
        IN( "nrows", "ncols" ),
        GIN(         "value"   ),
        GIN_TYPE(    "Numeric" ),
        GIN_DEFAULT( NODEF     ),
        GIN_DESC( "Matrix value." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MatrixSet" ),
        DESCRIPTION
        (
         "Create a Matrix from the given list of numbers.\n"
         "\n"
         "Usage:\n"
         "   MatrixSet(m1, [1, 2, 3; 4, 5, 6])\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(      "gout1"       ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( "The newly created matrix" ),
        IN(),
        GIN( "values"   ),
        GIN_TYPE(    "Matrix" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "The values of the newly created matrix. Elements are separated "
                  "by commas, rows by semicolons."),
        SETMETHOD( true )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MatrixUnitIntensity" ),
        DESCRIPTION
        (
         "Sets a matrix to hold unpolarised radiation with unit intensity.\n"
         "\n"
         "Works as MatrixPlanck where the radiation is set to 1.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "m"       ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( "Variable to initialize." ),
        IN( "stokes_dim" ),
        GIN(         "f"      ),
        GIN_TYPE(    "Vector" ),
        GIN_DEFAULT( NODEF    ),
        GIN_DESC( "Frequency vector." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Matrix1ColFromVector" ),
        DESCRIPTION
        (
         "Forms a matrix containing one column from a vector.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUT(),
        GOUT(      "m"      ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( "Variable to initialize." ),
        IN(),
        GIN(         "v"       ),
        GIN_TYPE(    "Vector" ),
        GIN_DEFAULT( NODEF    ),
        GIN_DESC( "The vector to be copied." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Matrix2ColFromVectors" ),
        DESCRIPTION
        (
         "Forms a matrix containing two columns from two vectors.\n"
         "\n"
         "The vectors are included as columns in the matrix in the same order\n"
         "as they are given.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUT(),
        GOUT(      "m"      ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( "Variable to initialize." ),
        IN(),
        GIN(         "v1"    , "v2"     ),
        GIN_TYPE(    "Vector", "Vector" ),
        GIN_DEFAULT( NODEF   , NODEF    ),
        GIN_DESC( "The vector to be copied into the first column.",
                  "The vector to be copied into the second column." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Matrix3ColFromVectors" ),
        DESCRIPTION
        (
         "Forms a matrix containing three columns from three vectors.\n"
         "\n"
         "The vectors are included as columns in the matrix in the same order\n"
         "as they are given.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUT(),
        GOUT(      "m"      ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( "Variable to initialize." ),
        IN(),
        GIN(         "v1"    , "v2"    , "v3"     ),
        GIN_TYPE(    "Vector", "Vector", "Vector" ),
        GIN_DEFAULT( NODEF   , NODEF   , NODEF    ),
        GIN_DESC( "The vector to be copied into the first column.",
                  "The vector to be copied into the second column.",
                  "The vector to be copied into the third column." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MatrixCompare" ),
        DESCRIPTION
        (
         "Checks the consistency between two matrices.\n" 
         "\n"
         "The two matrices are checked to not deviate outside the specified\n"
         "value (*maxdev*). An error is issued if this is not fulfilled.\n"
         "\n"
         "The main application of this method is to be part of the test\n"
         "control files, and then used to check that a calculated Jacobian\n"
         "is consistent with an old, reference, calculation.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( ),
        GIN( "matrix1", "matrix2", "maxabsdiff", "error_message" ),
        GIN_TYPE( "Matrix", "Matrix", "Numeric", "String" ),
        GIN_DEFAULT( NODEF, NODEF, NODEF, "" ),
        GIN_DESC( "A first matrix", "A second matrix", 
                  "Threshold for maximum absolute difference.",
                  "Additional error message.")
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Matrix1RowFromVector" ),
        DESCRIPTION
        (
         "Forms a matrix containing one row from a vector.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUT(),
        GOUT(      "m"      ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( "Variable to initialize." ),
        IN(),
        GIN(         "v"       ),
        GIN_TYPE(    "Vector" ),
        GIN_DEFAULT( NODEF    ),
        GIN_DESC( "The vector to be copied." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Matrix2RowFromVectors" ),
        DESCRIPTION
        (
         "Forms a matrix containing two rows from two vectors.\n"
         "\n"
         "The vectors are included as rows in the matrix in the same order\n"
         "as they are given.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUT(),
        GOUT(      "m"      ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( "Variable to initialize." ),
        IN(),
        GIN(         "v1"    , "v2"     ),
        GIN_TYPE(    "Vector", "Vector" ),
        GIN_DEFAULT( NODEF   , NODEF    ),
        GIN_DESC( "The vector to be copied into the first row.",
                  "The vector to be copied into the second row." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Matrix3RowFromVectors" ),
        DESCRIPTION
        (
         "Forms a matrix containing three rows from three vectors.\n"
         "\n"
         "The vectors are included as rows in the matrix in the same order\n"
         "as they are given.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUT(),
        GOUT(      "m"      ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( "Variable to initialize." ),
        IN(),
        GIN(         "v1"    , "v2"    , "v3"     ),
        GIN_TYPE(    "Vector", "Vector", "Vector" ),
        GIN_DEFAULT( NODEF   , NODEF   , NODEF    ),
        GIN_DESC( "The vector to be copied into the first row.",
                  "The vector to be copied into the second row.",
                  "The vector to be copied into the third row." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "mc_antennaSetGaussian" ),
        DESCRIPTION
        (
         "Makes mc_antenna (used by MCGeneral) a 2D Gaussian pattern.\n"
         "\n"
         "The gaussian antenna pattern is determined by *za_sigma* and\n"
         "*aa_sigma*, which represent the standard deviations in the\n"
         "uncorrelated bivariate normal distribution.\n"
         ),
        AUTHORS( "Cory Davis" ),
        OUT( "mc_antenna" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "za_sigma", "aa_sigma" ),
        GIN_TYPE(    "Numeric",  "Numeric" ),
        GIN_DEFAULT( NODEF,      NODEF ),
        GIN_DESC( "Width in the zenith angle dimension as described above.",
                  "Width in the azimuth angle dimension as described above." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "mc_antennaSetGaussianByFWHM" ),
        DESCRIPTION
        (
         "Makes mc_antenna (used by MCGeneral) a 2D Gaussian pattern.\n"
         "\n"
         "The gaussian antenna pattern is determined by *za_fwhm* and\n"
         "*aa_fwhm*, which represent the full width half maximum (FWHM)\n"
         "of the antenna response, in the zenith and azimuthal planes.\n"
         ),
        AUTHORS( "Cory Davis" ),
        OUT( "mc_antenna" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "za_fwhm", "aa_fwhm" ),
        GIN_TYPE(    "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF,     NODEF ),
        GIN_DESC( "Width in the zenith angle dimension as described above.",
                  "Width in the azimuth angle dimension as described above." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "mc_antennaSetPencilBeam" ),
        DESCRIPTION
        (
         "Makes mc_antenna (used by MCGeneral) a pencil beam.\n"
         "\n"
         "This WSM makes the subsequent MCGeneral WSM perform pencil beam\n"
         "RT calculations.\n" 
         ),
        AUTHORS( "Cory Davis" ),
        OUT( "mc_antenna" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "mc_IWP_cloud_opt_pathCalc" ),
        DESCRIPTION
        (
         "Calculates the FOV averaged ice water path and cloud optical path\n"
         "for a given viewing direction\n"
         ),
        AUTHORS( "Cory Davis" ),
        OUT( "mc_IWP", "mc_cloud_opt_path", "mc_IWP_error", 
             "mc_cloud_opt_path_error", "mc_iteration_count" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "mc_antenna", "sensor_pos", "sensor_los", "ppath_step_agenda", 
            "p_grid", "lat_grid", "lon_grid", "refellipsoid", "z_surface", 
            "z_field", "t_field", "vmr_field", "edensity_field", "f_index",
            "cloudbox_limits", "pnd_field", 
            "scat_data_mono", "particle_masses", "mc_seed" ),
        GIN( "max_iter" ),
        GIN_TYPE(    "Index" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Maximum number of iterations." )
        ));
  
  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "MCGeneral" ),
        DESCRIPTION
        ( "A generalised 3D reversed Monte Carlo radiative algorithm, that\n"
          "allows for 2D antenna patterns, surface reflection and arbitrary\n"
          "sensor positions.\n"
          "\n"
          "The main output variables *y* and *mc_error* represent the\n"
          "Stokes vector integrated over the antenna function, and the\n"
          "estimated error in this vector respectively.\n"
          "\n"
          "The WSV *mc_max_iter* describes the number of `photons\'\n"
          "used in the simulation (more photons means smaller *mc_error*).\n"
          "*mc_std_err* is the desired value of mc_error. *mc_max_time* is\n"
          "the maximum allowed number of seconds for MCGeneral. The method\n"
          "will terminate once any of the max_iter, std_err, max_time\n"
          "criteria are met. If negative values are given for these\n"
          "parameters then it isignored.\n"
          "\n"
          "Negative values of *mc_seed* seed the random number generator\n"
          "according to system time, positive *mc_seed* values are taken\n"
          "literally.\n"
          "\n"
          "Only \"1\" and \"RJBT\" are allowed for *y_unit*. The value of\n"
          "*mc_error* follows the selection for *y_unit* (both for in- and\n"
          "output.\n"
          ),
        AUTHORS( "Cory Davis" ),
        OUT( "y", "mc_iteration_count", "mc_error", "mc_points" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "mc_antenna", "f_grid", "f_index", "sensor_pos", "sensor_los", 
            "stokes_dim", "atmosphere_dim", "iy_space_agenda", "surface_rtprop_agenda", 
            "abs_mat_per_species_agenda", "p_grid",
            "lat_grid", "lon_grid", "z_field", "refellipsoid", "z_surface", 
            "t_field", "vmr_field", "cloudbox_on", "cloudbox_limits", 
            "pnd_field", "scat_data_mono", "basics_checked", "cloudbox_checked",
            "mc_seed", "y_unit", 
            "mc_std_err", "mc_max_time", "mc_max_iter"),//, "mc_z_field_is_1D" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "MCIPA" ),
        DESCRIPTION
        ( "A specialised 3D reversed Monte Carlo radiative algorithm, that\n"
          "mimics independent pixel appoximation simulations.\n"
          ),
        AUTHORS( "Cory Davis" ),
        OUT( "y", "mc_iteration_count", "mc_error", "mc_points" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "mc_antenna", "f_grid", "f_index", "sensor_pos", "sensor_los", 
            "stokes_dim", "atmosphere_dim", "iy_space_agenda", 
            "surface_rtprop_agenda",  
            "abs_mat_per_species_agenda", "ppath_step_agenda", "p_grid", "lat_grid",
            "lon_grid", "z_field", "refellipsoid", "z_surface", "t_field", 
            "vmr_field", "edensity_field", "cloudbox_limits", "pnd_field", 
            "scat_data_mono", "mc_seed", "y_unit",
            "mc_std_err", "mc_max_time", "mc_max_iter", "mc_z_field_is_1D" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "MCSetSeedFromTime" ),
        DESCRIPTION
        ( "Sets the value of mc_seed from system time\n" ),
        AUTHORS( "Cory Davis" ),
        OUT( "mc_seed" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "NumericAdd" ),
        DESCRIPTION
        (
         "Adds a numeric and a value (b = a+v).\n"
         "\n"
         "The result can either be stored in the same or another numeric.\n"
         "(a and b can be the same varible, but not b and v)\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "b"       ),
        GOUT_TYPE( "Numeric" ),
        GOUT_DESC( "Output numeric." ),
        IN(),
        GIN(      "a"      ,
                  "v" ),
        GIN_TYPE(    "Numeric",
                     "Numeric" ),
        GIN_DEFAULT( NODEF   ,
                     NODEF ),
        GIN_DESC( "Input numeric.",
                  "Value to add." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "NumericCompare" ),
        DESCRIPTION
        (
         "Checks the consistency between two numerics.\n" 
         "\n"
         "The two numerics are checked to not deviate outside the specified\n"
         "value (*maxabsdiff*). An error is issued if this is not fulfilled.\n"
         "\n"
         "The main application of this method is to be part of the test\n"
         "control files, and then used to check that a calculated value\n"
         "is consistent with an old, reference, value.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN( "n1", "n2", "maxabsdiff", "error_message" ),
        GIN_TYPE( "Numeric", "Numeric", "Numeric", "String" ),
        GIN_DEFAULT( NODEF, NODEF, "", "" ),
        GIN_DESC( "A first vector", "A second vector", 
                  "Threshold for maximum absolute difference.",
                  "Additional error message.")
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "NumericScale" ),
        DESCRIPTION
        (
         "Scales/multiplies a numeric with a value (b = a*v).\n"
         "\n"
         "The result can either be stored in the same or another numeric.\n"
         "(a and b can be the same varible, but not b and v)\n" 
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "b"       ),
        GOUT_TYPE( "Numeric" ),
        GOUT_DESC( "Output numeric." ),
        IN(),
        GIN(      "a"      ,
                  "v" ),
        GIN_TYPE(    "Numeric",
                     "Numeric" ),
        GIN_DEFAULT( NODEF   ,
                     NODEF ),
        GIN_DESC( "Input numeric.",
                  "Scaling value." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "NumericSet" ),
        DESCRIPTION
        (
         "Sets a numeric workspace variable to the given value.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "n"        ),
        GOUT_TYPE( "Numeric" ),
        GOUT_DESC( "Variable to initialize." ),
        IN(),
        GIN(         "value"   ),
        GIN_TYPE(    "Numeric" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "The value." ),
        SETMETHOD( true )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "nelemGet" ),
        DESCRIPTION
        (
         "Retrieve nelem from given variable and store the value in the\n"
         "workspace variable *nelem*\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "nelem" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "v"    ),
        GIN_TYPE(    ARRAY_GROUPS + ", Vector" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC(    "Variable to get the number of elements from." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "ncolsGet" ),
        DESCRIPTION
        (
         "Retrieve ncols from given variable and store the value in the\n"
         "workspace variable *ncols*\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "ncols" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "v" ),
        GIN_TYPE(    "Matrix, Sparse, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC(    "Variable to get the number of columns from." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "nrowsGet" ),
        DESCRIPTION
        (
         "Retrieve nrows from given variable and store the value in the\n"
         "workspace variable *nrows*\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "nrows" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "v" ),
        GIN_TYPE(    "Matrix, Sparse, Tensor3, Tensor4, Tensor5, Tensor6, Tensor7" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC(    "Variable to get the number of rows from." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "npagesGet" ),
        DESCRIPTION
        (
         "Retrieve npages from given variable and store the value in the\n"
         "workspace variable *npages*\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "npages" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "v" ),
        GIN_TYPE(    "Tensor3, Tensor4, Tensor5, Tensor6, Tensor7" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC(    "Variable to get the number of pages from." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "nbooksGet" ),
        DESCRIPTION
        (
         "Retrieve nbooks from given variable and store the value in the\n"
         "workspace variable *nbooks*\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "nbooks" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "v" ),
        GIN_TYPE(    "Tensor4, Tensor5, Tensor6, Tensor7" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC(    "Variable to get the number of books from." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "nshelvesGet" ),
        DESCRIPTION
        (
         "Retrieve nshelves from given variable and store the value in the\n"
         "workspace variable *nshelves*\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "nshelves" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "v" ),
        GIN_TYPE(    "Tensor5, Tensor6, Tensor7" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC(    "Variable to get the number of shelves from." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "nvitrinesGet" ),
        DESCRIPTION
        (
         "Retrieve nvitrines from given variable and store the value in the\n"
         "workspace variable *nvitrines*\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "nvitrines" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "v" ),
        GIN_TYPE(    "Tensor6, Tensor7" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC(    "Variable to get the number of vitrines from." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "nlibrariesGet" ),
        DESCRIPTION
        (
         "Retrieve nlibraries from given variable and store the value in the\n"
         "workspace variable *nlibraries*\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "nlibraries" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "v" ),
        GIN_TYPE(    "Tensor7" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC(    "Variable to get the number of libraries from." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "opt_prop_sptFromData" ),
        DESCRIPTION
        (
         "Calculates opticle properties for the single particle types.\n"
         "\n"
         "In this function the extinction matrix and the absorption vector\n"
         "are calculated in the laboratory frame. An interpolation of the\n"
         "data on the actual frequency is the first step in this function.\n"
         "The next step is a transformation from the database coordinate\n"
         "system to the laboratory coordinate system.\n"
         "\n"
         "Output of the function are *ext_mat_spt* and *abs_vec_spt* which\n"
         "hold the optical properties for a specified propagation direction\n"
         "for each particle type.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "ext_mat_spt", "abs_vec_spt" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "ext_mat_spt", "abs_vec_spt", "scat_data_raw", "scat_za_grid", 
            "scat_aa_grid", "scat_za_index", "scat_aa_index", 
            "f_index", "f_grid", "rte_temperature",
            "pnd_field", "scat_p_index", "scat_lat_index", "scat_lon_index" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "opt_prop_sptFromMonoData" ),
        DESCRIPTION
        (
         "Calculates optical properties for the single particle types.\n"
         "\n"
         "As *opt_prop_sptFromData* but no frequency interpolation is\n"
         "performed. The single scattering data is here obtained from\n"
         "*scat_data_mono*, instead of *scat_data_raw*.\n"
         ),
        AUTHORS( "Cory Davis" ),
        OUT( "ext_mat_spt", "abs_vec_spt" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "ext_mat_spt", "abs_vec_spt", "scat_data_mono", "scat_za_grid", 
            "scat_aa_grid", "scat_za_index", "scat_aa_index", "rte_temperature",
            "pnd_field", "scat_p_index", "scat_lat_index", "scat_lon_index" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "output_file_formatSetAscii" ),
        DESCRIPTION
        (
         "Sets the output file format to ASCII.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "output_file_format" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "output_file_formatSetBinary" ),
        DESCRIPTION
        (
         "Sets the output file format to binary.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "output_file_format" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "output_file_formatSetZippedAscii" ),
        DESCRIPTION
        (
         "Sets the output file format to zipped ASCII.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "output_file_format" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
    
    md_data_raw.push_back
    ( MdRecord
      ( NAME( "ParticleSpeciesInit" ),
        DESCRIPTION
        (
         "Initializes empty *part_species* array.\n"
         ),
        AUTHORS( "Daniel Kreyling" ),
        OUT( "part_species" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
    
    
    md_data_raw.push_back
    ( MdRecord
      ( NAME( "ParticleSpeciesSet" ),
        DESCRIPTION
        (
         "With this function, the user specifies settings for the \n"
         "particle number density calculations using *pnd_fieldSetup*.\n"
         "The input is an ArrayOfString that needs to be in a specific format,\n"
         "for details, see WSV *part_species*.\n"
         "\n"         
         "*Example:* \t ['IWC-MH97-Ice-0.1-200', 'LWC-H98_STCO-Water-0.1-50'] \n"
         "\n"
         "NOTE: The order of the Strings need to match the order of the\n"
         "*atm_fields_compact* field names, their number determines how many fields\n"
         "of *atm_fields_compact* are considered particle profiles.\n"
         ),
        AUTHORS( "Daniel Kreyling" ),
        OUT( "part_species" ),
        GOUT( ),
        GOUT_TYPE( ),
        GOUT_DESC( ),
        IN(),
        GIN( "names" ),
        GIN_TYPE(  "ArrayOfString" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC("Array of pnd calculation parameters." )
        ));
 
    
   md_data_raw.push_back
    ( MdRecord
      ( NAME( "ParticleTypeAddAll" ),
        DESCRIPTION
        (
         "Reads single scattering data and particle number densities.\n"
         "\n"
         "The WSV *pnd_field_raw* containing particle number densities for all\n"
         "scattering particle species can be generated outside ARTS, for example by using\n"
         "PyARTS. This method needs as input an XML-file containing an array of filenames\n"
         "(ArrayOfString) of single scattering data and a file containing the corresponding\n"
         "*pnd_field_raw*. In contrast to the scattering data, all corresponding pnd-fields\n"
         "are stored in a single XML-file containing an ArrayofGriddedField3\n"
         "\n"
         "Very important note:\n"
         "The order of the filenames for the scattering data files has to\n"
         "correspond to the order of the pnd-fields, stored in the variable\n"
         "*pnd_field_raw*.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "scat_data_raw", "pnd_field_raw" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim", "f_grid", "p_grid", "lat_grid", "lon_grid", 
            "cloudbox_limits" ),
        GIN(         "filename_scat_data", "filename_pnd_field" ),
        GIN_TYPE(    "String",             "String"             ),
        GIN_DEFAULT( NODEF,                NODEF                ),
        GIN_DESC( "File containing single scattering data.",
                  "File including *pnd_field_raw*." 
                  )
        ));


  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ParticleTypeAdd" ),
        DESCRIPTION
        (
         "This method reads single scattering data and the corresonding\n"
         "particle number density fields.\n"
         "\n"
         "The methods reads the  specified files and appends the obtained data\n"
         "to *scat_data_raw* and *pnd_field_raw*.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "scat_data_raw", "pnd_field_raw" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim", "f_grid", "p_grid", "lat_grid", "lon_grid", 
            "cloudbox_limits" ),
        GIN(         "filename_scat_data", "filename_pnd_field" ),
        GIN_TYPE(    "String",             "String"             ),
        GIN_DEFAULT( NODEF,                NODEF                ),
        GIN_DESC( "Filename of single scattering data.",
                  "Filename of the corresponding pnd_field." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ParticleTypeInit" ),
        DESCRIPTION
        (
         "Initializes *scat_data_raw* and *pnd_field_raw*.\n"
         "\n"
         "This method initializes variables containing data about the\n"
         "optical properties of particles (*scat_data_raw*) and about the\n"
         "particle number distribution (*pnd_field_raw*)\n"
         "\n"
         "This method has to be executed before executing e.g.\n"
         "*ParticleTypeAdd*.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "scat_data_raw", "pnd_field_raw" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(), 
        GIN_TYPE(), 
        GIN_DEFAULT(),
        GIN_DESC()
        ));
    
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_matCalc" ),
        DESCRIPTION
        (
         "This function sums up the phase matrices for all particle\n"
         "types weighted with particle number density.\n"
         ),
        AUTHORS( "Sreerekha T.R." ),
        OUT( "pha_mat" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "pha_mat_spt", "pnd_field", "atmosphere_dim", "scat_p_index",
            "scat_lat_index", "scat_lon_index" ),
        GIN(),
        GIN_TYPE(), 
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_matExtractManually" ),
        DESCRIPTION
        (
         "A simple function for manually extract a single phase matrix.\n"
         "\n"
         "The function returns the phase matrix for a single particle, for\n"
         "scattering from (za_in,aa_in) to (za_out,aa_out).\n"
         "\n"
         "Only a single particle type is handled and *scat_data_raw* must\n"
         "have length 1. The frequency is selected by *f_grid* and *f_index*.\n"
         "The temperature is set by *rte_temperature*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( ),
        GOUT( "pha_mat_single" ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( 
            "Phase matrix for a single frequency and combination of angles" ),
        IN( "f_grid", "f_index", "stokes_dim", "scat_data_raw", 
            "rte_temperature" ),
        GIN( "za_out", "aa_out", "za_in", "aa_in" ),
        GIN_TYPE( "Numeric", "Numeric", "Numeric", "Numeric" ), 
        GIN_DEFAULT( NODEF, NODEF, NODEF, NODEF ),
        GIN_DESC( "Outgoing zenith angle", "Outgoing azimuth angle",
                  "Incoming zenith angle", "Incoming azimuth angle" )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_mat_sptFromData" ),
        DESCRIPTION
        (
         "Calculation of the phase matrix for the single particle types.\n"
         "\n"
         "This function can be used in *pha_mat_spt_agenda* as part of\n"
         "the calculation of the scattering integral.\n"
         "\n"
         "The interpolation of the data on the actual frequency is the first\n"
         "step in this function. This is followed by a transformation from the\n"
         "database coordinate system to the laboratory coordinate system.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "pha_mat_spt" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "pha_mat_spt", "scat_data_raw", "scat_za_grid", "scat_aa_grid", 
            "scat_za_index", "scat_aa_index", "f_index", "f_grid",
            "rte_temperature", "pnd_field", "scat_p_index", "scat_lat_index",
            "scat_lon_index" ),
        GIN(),
        GIN_TYPE(), 
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_mat_sptFromMonoData" ),
        DESCRIPTION
        (
         "Calculation of the phase matrix for the single particle types.\n"
         "\n"
         "This function is the monochromatic version of *pha_mat_sptFromData*.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "pha_mat_spt" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "pha_mat_spt", "scat_data_mono", "doit_za_grid_size",
            "scat_aa_grid", "scat_za_index", "scat_aa_index", "rte_temperature",
            "pnd_field", "scat_p_index", "scat_lat_index", "scat_lon_index" ),
        GIN(),
        GIN_TYPE(), 
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_mat_sptFromDataDOITOpt" ),
        DESCRIPTION
        (
         "Calculation of the phase matrix for the single particle types.\n"
         "\n"
         "In this function the phase matrix is extracted from\n"
         "*pha_mat_sptDOITOpt*. It can be used in the agenda\n"
         "*pha_mat_spt_agenda*. This method must be used in \n "
         "combination with *DoitScatteringDataPrepare*.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "pha_mat_spt" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "pha_mat_spt", "pha_mat_sptDOITOpt", "scat_data_mono", 
            "doit_za_grid_size",
            "scat_aa_grid", 
            "scat_za_index", "scat_aa_index", "rte_temperature",
            "pnd_field", "scat_p_index", "scat_lat_index", "scat_lon_index" ),
        GIN(),
        GIN_TYPE(), 
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pnd_fieldCalc" ),
        DESCRIPTION
        ( "Interpolates the particle number density fields.\n"
          "\n"
          "This methods interpolates the particle number density field\n"
          "from the raw data *pnd_field_raw* to obtain *pnd_field*.\n"
          ),
        AUTHORS( "Sreerekha T.R.", "Claudia Emde" ),
        OUT( "pnd_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "p_grid", "lat_grid", "lon_grid", "pnd_field_raw", "atmosphere_dim",
            "cloudbox_limits" ),
        GIN(),
        GIN_TYPE(), 
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pnd_fieldExpand1D" ),
        DESCRIPTION
        (
         "Maps a 1D pnd_field to 2D or 3D pnd_field.\n"
         "\n"
         "This method takes a 1D *pnd_field* and converts it to 2D or 3D\n"
         "\"cloud\". It is assumed that a complet 1D case has been created\n"
         "and after this *atmosphere_dim*, *lat_grid*, *lon_grid* and\n"
         "*cloudbox_limits* have been changed to a 2D or 3D case. This\n"
         "without changing the vertical extension of the cloudbox.\n"
         "\n"
         "No modification of *pnd_field* is made for the pressure dimension.\n"
         "At the latitude and longitude cloudbox edges *pnd_field* is set to\n"
         "zero. This corresponds to nzero=1. If you want a larger margin between\n"
         "the lat and lon cloudbox edgess and the \"cloud\" you increase\n"
         "*nzero*, where *nzero* is the number of grid points for which\n"
         "*pnd_field* shall be set to 0, counted from each lat and lon edge.\n"
         "\n"
         "See further *AtmFieldsExpand1D*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "pnd_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "pnd_field", "atmosphere_dim", "cloudbox_checked", "cloudbox_on", 
            "cloudbox_limits" ),
        GIN( "nzero" ),
        GIN_TYPE( "Index"),
        GIN_DEFAULT( "1" ),
        GIN_DESC( "Number of zero values inside lat and lon limits." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pnd_fieldSetup" ),
        DESCRIPTION
        (
         "Calculation of *pnd_field* using ScatteringMetaData and *massdensity_field*.\n"
         "\n"
         "The WSM first checks if cloudbox is empty. If so, the pnd calculations\n"
         "will be skipped.\n"
         "The *cloudbox_limits* are used to determine the p, lat and lon size for\n"
         "the *pnd_field* tensor.\n"
         "Currently there are three particle size distribution (PSD) parameterisations\n"
         "implemented:\n"
         "\t1. 'MH97' for ice particles. Parameterisation in temperature and mass content.\n"
         "\t Using a first-order gamma distribution for particles smaller than \n"
         "\t 100 microns (melted diameter) and a lognormal distribution for\n"
         "\t particles bigger 100 microns. Values from both modes are cumulative.\n"
         "\t See internal function 'IWCtopnd_MH97' for implementation/units/output.\n"
         "\t (src.: McFarquhar G.M., Heymsfield A.J., 1997)"
         "\n"
	 "\t2. 'H11' for cloud ice and precipitating ice (snow). H11 is NOT dependent\n"
	 "\t on mass content of ice/snow, but only on atmospheric temperature.\n"
	 "\t The PSD is scaled to the current IWC/Snow density in an additional step.\n"
	 "\t See internal function 'pnd_H11' and 'scale_H11' for implementation/units/output.\n"
	 "\t (src.: Heymsfield A.J., 2011, not published yet)\n"
         "\t3. 'H98_STCO' for liquid water clouds. Using a gamma distribution with\n"
         "\t parameters from Hess et al., 1998, continental stratus.\n"
         "\t See internal function 'LWCtopnd' for implementation/units/output.\n"
         "\t (src.: Deirmendjian D., 1963 and Hess M., et al 1998)\n"
         "\n"
         "According to the selection criteria in *part_species*, the first specified\n" 
         "psd parametrisation is selected together with all particles of specified phase\n"
         "and size. Then pnd calculations are performed on all levels inside the cloudbox.\n"
         "The *massdensity_field* input weights the pnds by the amount of scattering\n" 
         "particles in each gridbox inside the cloudbox. Where *massdensity_field* is zero,\n"
         "the *pnd_field* will be zero as well.\n"
         "Subsequently the pnd values get written to *pnd_field*.\n"
         "\n"
         "Now the next selection criteria string in *part_species* is used to repeat\n"
         "the process.The new pnd values will be appended to the existing *pnd_field*.\n"
         "And so on...\n"
         "\n"
         "NOTE: the order of scattering particle profiles in *massdensity_field* has to\n"
         "fit the order of part_species tags!\n"
         ),
        AUTHORS( "Daniel Kreyling" ),
        OUT( "pnd_field"),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim","cloudbox_on", "cloudbox_limits", "massdensity_field", "t_field", "scat_data_meta_array", "part_species", "scat_data_nelem" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));  
    
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pnd_fieldZero" ),
        DESCRIPTION
        (
         "Sets *pnd_field* to hold only zeros.\n"
         "\n"
         "Scattering calculations using the DOIT method include\n"
         "interpolation errors. If one is interested in this effect, one\n"
         "should compare the DOIT result with a clearsky calculation using\n"
         "an empty cloudbox. That means that the iterative method is\n"
         "performed for a cloudbox including no particles. This method sets\n"
         "the particle number density field to zero and creates a\n"
         "dummy *scat_data_raw* structure. \n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "pnd_field", "scat_data_raw" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "p_grid", "lat_grid", "lon_grid" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ppathCalc" ),
        DESCRIPTION
        (
         "Stand-alone calculation of propagation paths.\n"
         "\n"
         "Beside a few checks of input data, the only operation of this\n"
         "method is to execute *ppath_agenda*.\n"
         "\n"
         "Propagation paths are normally calculated as part of the radiative\n"
         "transfer calculations, and this method is not part of the control\n"
         "file. A reason to call this function directly would be to obtain a\n"
         "propagation path for plotting. Anyhow, use this method instead\n"
         "of calling e.g.*ppathStepByStep directly.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "ppath" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "ppath_agenda", "basics_checked", "t_field", "z_field", "vmr_field",
            "edensity_field", "f_index", "cloudbox_on", "cloudbox_checked", 
            "ppath_inside_cloudbox_do", "mblock_index", "rte_pos", "rte_los" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ppathFromRtePos2" ),
        DESCRIPTION
        (
         "Determines the propagation path from *rte_pos2* to *rte_pos*.\n"
         "\n"
         "FIXME: Need an update!\n"
         "\n"
         "The propagation path linking *rte_pos* and *rte_pos2* is calculated\n"
         "and returned. The method assumes that refraction is considered and\n"
         "determines the path in a pure numerical manner. A simple algorithm\n"
         "is applied. Repeated propagation path calculations (starting at\n"
         "*rte_pos*) are performed. The closest distance between the path and"
         "*rte_pos2* is converted to a correction for *rte_los* and a new\n"
         "path is calculated. This is repeated until the convergence\n"
         "criterion is fulfilled.\n"
         "\n"
         "The standard application of this method should be to radio link\n"
         "calculations, where *rte_pos2* corresponds to a transmitter, and\n"
         "*rte_pos* to the receiver/sensor.\n"
         "\n"
         "The WSV *rte_los* is both input and output. The input *rte_los*\n"
         "shall contain an useful \"first guess\" for the line-of-sight at\n"
         "*rte_pos*. The output *rte_los* is the line-of-sight for the\n"
         "returned *ppath*. The best choice should be to first set *rte_los*\n"
         "to the geometrical line-of-sight (can be done by e.g. \n"
         "*rte_losGeometricFromRtePosToRtePos2*). If this is done, the\n"
         "difference between in- and output *rte_los* can be used to obtain\n"
         "the so called bending angle.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "ppath", "rte_los" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "ppath_step_agenda", "basics_checked", "atmosphere_dim", "p_grid", 
            "lat_grid", "lon_grid", "t_field", "z_field", "vmr_field", 
            "edensity_field", "f_index", "refellipsoid", "z_surface", 
            "rte_pos", "rte_pos2" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
  
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ppathStepByStep" ),
        DESCRIPTION
        (
         "Standard method for calculation of propagation paths.\n"
         "\n"
         "This method calculates complete propagation paths in a stepwise\n"
         "manner. Each step is denoted as a \"ppath_step\" and is the path\n"
         "through/inside a single grid box.\n"
         "\n"
         "The definition of a propgation path cannot be accomodated here.\n"
         "For more information read the chapter on propagation paths in the\n"
         "ARTS user guide.\n"
         "\n"
         "This method shuld never be called directly. Use instead *ppathCalc*\n"
         "if you want to extract propagation paths.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "ppath" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "ppath_step_agenda", "ppath_inside_cloudbox_do", "atmosphere_dim", 
            "p_grid", "lat_grid", "lon_grid", "t_field", "z_field", "vmr_field",
            "edensity_field", "f_index", "refellipsoid", "z_surface", 
            "cloudbox_on", "cloudbox_limits", "rte_pos", "rte_los" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ppath_stepGeometric" ),
        DESCRIPTION
        (
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
         "*ppath_step_agenda*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "ppath_step" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "ppath_step", "atmosphere_dim", "lat_grid", "lon_grid", 
            "z_field", "refellipsoid", "z_surface", "ppath_lmax" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ppath_stepRefractionBasic" ),
        DESCRIPTION
        (
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
         "*ppath_lmax*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "ppath_step" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "refr_index_agenda", "ppath_step", "atmosphere_dim", "p_grid", 
            "lat_grid", "lon_grid", "z_field", "t_field", "vmr_field", 
            "edensity_field", "refellipsoid", "z_surface", "f_index",
            "ppath_lmax", "ppath_lraytrace" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "Print" ),
        DESCRIPTION
        (
         "Prints a variable on the screen.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "v"   ,
                     "level" ),
        GIN_TYPE(    "Any",
                     "Index" ),
        GIN_DEFAULT( NODEF,
                     "1" ),
        GIN_DESC(    "Variable to be printed.",
                     "Output level to use." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "PrintWorkspace" ),
        DESCRIPTION
        (
         "Prints a list of the workspace variables.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN( "only_allocated", "level" ),
        GIN_TYPE(    "Index",          "Index" ),
        GIN_DEFAULT( "1",              "1" ),
        GIN_DESC( "Flag for printing either all variables (0) or only "
                  "allocated ones (1).",
                  "Output level to use." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  ),
        PASSWORKSPACE( true  )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "p_gridFromZRaw" ),
        DESCRIPTION
        (
         "Sets *p_grid* according to input atmosphere's raw z_field, derived\n"
         "e.g. from *AtmRawRead*.\n"
         "Attention: as default only pressure values for altitudes >= 0 are\n"
         "extracted. If negative altitudes shall also be selected, set no_neg=0.\n"
         ),
        AUTHORS( "Claudia Emde, Jana Mendrok" ),
        OUT( "p_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "z_field_raw" ),
        GIN(         "no_negZ" ),
        GIN_TYPE(    "Index" ),
        GIN_DEFAULT( "1" ),
        GIN_DESC(    "Exclude negative altitudes." )
        ));
 
  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "p_gridFromGasAbsLookup" ),
        DESCRIPTION
        (
         "Sets *p_grid* to the pressure grid of *abs_lookup*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "p_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_lookup" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ReadNetCDF" ),
        DESCRIPTION
        (
         "Reads a workspace variable from a NetCDF file.\n"
         "\n"
         "This method can read variables of any group.\n"
         "\n"
         "If the filename is omitted, the variable is read\n"
         "from <basename>.<variable_name>.nc.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(      "v"    ),
        GOUT_TYPE( "Vector, Matrix, Tensor3, Tensor4, Tensor5, ArrayOfVector,"
                   "ArrayOfMatrix, GasAbsLookup" ),
        GOUT_DESC( "Variable to be read." ),
        IN(),
        GIN(         "filename" ),
        GIN_TYPE(    "String"   ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC(    "Name of the NetCDF file." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  ),
        PASSWORKSPACE(  false ),
        PASSWSVNAMES(   true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ReadXML" ),
        DESCRIPTION
        (
         "Reads a workspace variable from an XML file.\n"
         "\n"
         "This method can read variables of any group.\n"
         "\n"
         "If the filename is omitted, the variable is read\n"
         "from <basename>.<variable_name>.xml.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(      "v"    ),
        GOUT_TYPE( "Any" ),
        GOUT_DESC( "Variable to be read." ),
        IN(),
        GIN(         "filename" ),
        GIN_TYPE(    "String"   ),
        GIN_DEFAULT( "" ),
        GIN_DESC(    "Name of the XML file." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  ),
        PASSWORKSPACE(  false ),
        PASSWSVNAMES(   true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "refellipsoidEarth" ),
        DESCRIPTION
        (
         "Earth reference ellipsoids.\n"
         "\n"
         "The reference ellipsoid (*refellipsoid*) is set to model the Earth,\n"
         "following different models. The options are:\n"
         "\n"
         "   \"Sphere\" : A spherical Earth. The radius is set following\n"
         "      the value set for the Earth radius in constants.cc.\n"
         "\n"
         "   \"WGS84\" : The reference ellipsoid used by the GPS system.\n"
         "      Should be the standard choice for a non-spherical Earth.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "refellipsoid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(  ),
        GIN( "model" ),
        GIN_TYPE(    "String" ),
        GIN_DEFAULT( "Sphere" ),
        GIN_DESC( "Model ellipsoid to use. Options listed above." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "refellipsoidForAzimuth" ),
        DESCRIPTION
        (
         "Conversion of 3D ellipsoid to 1D curvature radius.\n"
         "\n"
         "Calculates the curvature radius for the given latitude and azimuth\n"
         "angle, and uses this to set a spherical reference ellipsoid\n"
         "suitable for 1D calculations. The curvature radius is a better\n"
         "local approximation than using the local ellipsoid radius.\n"
         "\n"
         "The used expression assumes a geodetic latitude, but also\n"
         "latitudes should be OK as using this method anyhow signifies\n"
         "an approximation.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "refellipsoid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "refellipsoid" ),
        GIN( "latitude", "azimuth" ),
        GIN_TYPE( "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF, NODEF ),
        GIN_DESC( "Latitude.", "Azimuth angle." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "refellipsoidJupiter" ),
        DESCRIPTION
        (
         "Jupiter reference ellipsoids.\n"
         "\n"
         "The reference ellipsoid (*refellipsoid*) is set to model Jupiter,\n"
         "folowing different models. The options are:\n"
         "\n"
         "   \"Sphere\" : A spherical planet. The radius is taken from a\n"
         "      report of the IAU/IAG Working Group.\n"
         "\n"
         "   \"Ellipsoid\" : A reference ellipsoid with parameters taken from\n"
         "      a report of the IAU/IAG Working Group.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "refellipsoid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(  ),
        GIN( "model" ),
        GIN_TYPE(    "String" ),
        GIN_DEFAULT( "Sphere" ),
        GIN_DESC( "Model ellipsoid to use. Options listed above." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "refellipsoidMars" ),
        DESCRIPTION
        (
         "Mars reference ellipsoids.\n"
         "\n"
         "The reference ellipsoid (*refellipsoid*) is set to model Mars,\n"
         "folowing different models. The options are:\n"
         "\n"
         "   \"Sphere\" : A spherical planet. The radius is taken from a\n"
         "      report of the IAU/IAG Working Group.\n"
         "\n"
         "   \"Ellipsoid\" : A reference ellipsoid with parameters taken from\n"
         "      a report of the IAU/IAG Working Group.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "refellipsoid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(  ),
        GIN( "model" ),
        GIN_TYPE(    "String" ),
        GIN_DEFAULT( "Sphere" ),
        GIN_DESC( "Model ellipsoid to use. Options listed above." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "refellipsoidMoon" ),
        DESCRIPTION
        (
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
         "      defines the Moon ellipsoid to be a sphere.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "refellipsoid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(  ),
        GIN( "model" ),
        GIN_TYPE(    "String" ),
        GIN_DEFAULT( "Sphere" ),
        GIN_DESC( "Model ellipsoid to use. Options listed above." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "refellipsoidOrbitPlane" ),
        DESCRIPTION
        (
         "Conversion of 3D ellipsoid to 2D orbit track geometry.\n"
         "\n"
         "Determines an approximate reference ellipsoid following an orbit\n"
         "track. The new ellipsoid is determined simply, by determining the\n"
         "radius at the maximum latitude and from this value calculate a new\n"
         "new eccentricity. The orbit is specified by giving the orbit\n"
         "inclination (*orbitinc*), that is normally a value around 100 deg\n"
         "for polar sun-synchronous orbits.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "refellipsoid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "refellipsoid" ),
        GIN( "re" ),
        GIN_TYPE(    "Numeric" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Orbit inclination." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "refellipsoidSet" ),
        DESCRIPTION
        (
         "Manual setting of the reference ellipsoid.\n"
         "\n"
         "The two values of *refellipsoid* can here be set manually. The two\n"
         "arguments correspond directly to first and second element of\n"
         "*refellipsoid*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "refellipsoid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(  ),
        GIN( "re", "e" ),
        GIN_TYPE(    "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF, "0" ),
        GIN_DESC( "Average or equatorial radius.", "Eccentricity" )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "refellipsoidVenus" ),
        DESCRIPTION
        (
         "Venus reference ellipsoids.\n"
         "\n"
         "The reference ellipsoid (*refellipsoid*) is set to model Venus,\n"
         "folowing different models. The options are:\n"
         "\n"
         "   \"Sphere\" : A spherical planet. The radius is taken from a\n"
         "      report of the IAU/IAG Working Group.\n"
         "\n"
         "According to the report used above, the Venus ellipsoid lacks\n"
         "eccentricity and no further models should be required.\n"         
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "refellipsoid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(  ),
        GIN( "model" ),
        GIN_TYPE(    "String" ),
        GIN_DEFAULT( "Sphere" ),
        GIN_DESC( "Model ellipsoid to use. Options listed above." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "refr_indexFreeElectrons" ),
        DESCRIPTION
        (
         "Microwave refractive index due to free electrons.\n"
         "\n"
         "The refractive index of free electrons is added to *refr_index*.\n"
         "To obtain the complete value, *refr_index* should be set to 1\n"
         "before calling this WSM. This applies also to *refr_index_group*.\n"
         "\n"
         "The expression applied is n=sqrt(1-wp^2/w^2) where wp is the plasma\n"
         "frequency, and w is the angular frequency (the function returns\n"
         "n-1, that here is slightly negative). This expressions is found in\n"
         "many textbooks, e.g. Rybicki and Lightman (1979). The above refers\n"
         "to *refr_index*. *refr_index_group* is sqrt(1+wp^2/w^2).\n"
         "\n"
         "The expression is dispersive. The frequency applied is selected as\n"
         "follows. If *f_index* < 0, the mean of first and last element of\n"
         "*f_grid* is selected. Otherwise, *f_index* specifies the element\n"
         "of *f_grid* to extract. The applied frequency must be at least\n"
         "twice the plasma frequency.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "refr_index", "refr_index_group" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "refr_index", "refr_index_group", "f_grid", "f_index", 
            "rte_edensity" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "refr_indexIR" ),
        DESCRIPTION
        (
         "Calculates the IR refractive index due to gases in the\n"
         "Earth's atmosphere.\n"
         "\n"
         "Only refractivity of dry air is considered. The formula used is\n"
         "contributed by Michael Hoefner, Forschungszentrum Karlsruhe.\n"
         "\n"
         "The refractivity of dry air is added to *refr_index*. To obtain\n"
         "the complete value, *refr_index* should be set to 1 before\n"
         "calling this WSM. This applies also to *refr_index_group*.\n"
         "\n"
         "The expression used is non-dispersive. Hence, *refr_index* and\n"
         "*refr_index_group* are identical.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUT( "refr_index", "refr_index_group" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "refr_index", "refr_index_group", "rte_pressure", 
            "rte_temperature" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "refr_indexThayer" ),
        DESCRIPTION
        (
         "Microwave refractive index due to gases in the Earth's atmosphere.\n"
         "\n"
         "The refractivity of dry air and water vapour is added to\n"
         "*refr_index*. To obtain the complete value, *refr_index* should\n"
         "be set to 1 before calling this WSM. This applies also to\n"
         "*refr_index_group.\n"
         "\n"
         "The expression used is non-dispersive. Hence, *refr_index* and\n"
         "*refr_index_group* are identical.\n"
         "\n"
         "The parameterisation of Thayer (Radio Science, 9, 803-807, 1974)\n"
         "is used. See also Eq. 3 and 5 of Solheim et al. (JGR, 104,\n"
         "pp. 9664).\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "refr_index", "refr_index_group" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "refr_index", "refr_index_group", "rte_pressure", 
            "rte_temperature", "rte_vmr_list", "abs_species" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_losGeometricFromRtePosToRtePos2" ),
        DESCRIPTION
        (
         "The geometric line-of-sight between two points.\n"
         "\n"
         "The method sets *rte_los* to the line-of-sight, at *rte_pos*,\n"
         "that matches the geometrical propagation path between *rte_pos*\n"
         "and *rte_pos2*.\n"
         "\n"
         "The standard case should be that *rte_pos2* corresponds to a\n"
         "transmitter, and *rte_pos* to the receiver/sensor.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "rte_los" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim", "lat_grid", "lon_grid", "refellipsoid", 
            "rte_pos", "rte_pos2" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_losSet" ),
        DESCRIPTION
        (
         "Sets *rte_los* to the given angles.\n"
         "\n"
         "The azimuth angle is ignored for 1D and 2D.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "rte_los" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim" ),
        GIN( "za",      "aa"      ),
        GIN_TYPE(    "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF,     NODEF ),
        GIN_DESC( "Zenith angle of sensor line-of-sight.",
                  "Azimuth angle of sensor line-of-sight." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_posSet" ),
        DESCRIPTION
        (
         "Sets *rte_pos* to the given co-ordinates.\n"
         "\n"
         "The longitude is ignored for 1D and 2D, and the latitude is also \n"
         "ignored for 1D.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "rte_pos" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim" ),
        GIN( "z",  "lat",     "lon"     ),
        GIN_TYPE(    "Numeric", "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF,     NODEF,     NODEF ),
        GIN_DESC( "Geometrical altitude of sensor position.",
                  "Latitude of sensor position.",
                  "Longitude of sensor position." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ScatteringDisort" ),
        DESCRIPTION
        (
         "Calls DISORT RT solver from ARTS.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "scat_i_p", "scat_i_lat", "scat_i_lon", 
             "f_index", "scat_data_mono", "doit_i_field1D_spectrum" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "cloudbox_limits", "stokes_dim", "opt_prop_part_agenda", 
            "abs_mat_per_species_agenda", "spt_calc_agenda", "pnd_field", "t_field",
            "z_field", "p_grid", "vmr_field", "scat_data_raw", "f_grid", 
            "scat_za_grid", "surface_emissivity_DISORT" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ScatteringDoit" ),
        DESCRIPTION
        (
         "Main DOIT method.\n"
         "\n"
         "This method executes *doit_mono_agenda* for each frequency\n"
         "in *f_grid*. The output is the radiation field inside the cloudbox\n"
         "(*doit_i_field*) and on the cloudbox boundary (*scat_i_p* (1D),\n"
         "*scat_i_lat* and *scat_i_lon* (3D)).\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT( "doit_i_field", "scat_i_p", "scat_i_lat", "scat_i_lon",
             "doit_i_field1D_spectrum" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "cloudbox_on", "f_grid", "scat_i_p", "scat_i_lat", "scat_i_lon",
            "doit_mono_agenda", "doit_is_initialized" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
    
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ScatteringParticleTypeAndMetaRead" ),
        DESCRIPTION
        (
         "Reads single scattering data and scattering meta data.\n"
         "\n"
         "This method's input needs two XML-files, one containing an array \n"
         "of path/filenames (*ArrayOfString*) of single scattering data and the \n"
         "corresponding path/filenames to scattering meta data.\n"
         "For each single scattering file, there needs to be exactly one\n"
         "scattering meta data file.\n"
         "\n"
         "Currently particles of phase ice and/or water can be added for the same calculation.\n"
         "It is also possible to read *SingleScatteringData* for different shapes of\n"
         "ice particles. But all ice particels will share the same IWC, while performing\n"
         "the *pnd_field* calculations with *pnd_fieldSetup*.\n"
         "Also make sure, that two scattering particles of the same phase are never equal\n"
         "in size. This will break the calculations in *pnd_fieldSetup*\n"
         "\n"
         "Very important note:\n"
         "The order of the filenames for the single scattering data files has to\n"
         "exactly correspond to the order of the scattering meta data files.\n"
         ),
        AUTHORS( "Daniel Kreyling" ),
        OUT( "scat_data_raw", "scat_data_meta_array" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "f_grid" ),
        GIN(         "filename_scat_data", "filename_scat_meta_data" ),
        GIN_TYPE(    "String",             "String"             ),
        GIN_DEFAULT( NODEF,                NODEF                ),
        GIN_DESC( "File containing single scattering data file names.",
                  "File containing scattering meta data file names." 
                  )
        ));

 
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ScatteringParticlesSelect" ),
        DESCRIPTION
        (
         "This method is a selection function for scattering particles.\n"
         "\n"
         "In *part_species* the user defines selection criteria for:\n"
         "\t...which type of scattering particle profile\n"
         "\t...what particle size ditribution parametrisation\n"
         "\t...the minimum and maximum size of the particle (in terms of volume equivalent radius)\n"
         "to use in the scattering calculations.\n"
   "Minimum and maximum size may be omitted or symbol \"*\" be used as a wildcard.\n"
         "\n"
         "The scattering particle arrays, *scat_data_raw* and *scat_data_meta_array*\n"
         "are searched for particles, that fullfill the selection criteria. \n"
         "Only these particles will be used for scattering calculations.\n"
         "\n"
         "Additionaly an *ArrayOfIndex* *scat_data_nelem* is created. This Array\n"
         "stores the number of scattering particles, that have been selected by each\n"
         "selection string in *part_species*\n"
         ),
        AUTHORS( "Daniel Kreyling" ),
        OUT( "scat_data_raw", "scat_data_meta_array", "scat_data_nelem" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "part_species", "scat_data_raw", "scat_data_meta_array" ),
        GIN(    ),
        GIN_TYPE(),
        GIN_DEFAULT(   ),
        GIN_DESC(   )
         ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "scat_data_monoCalc" ),
        DESCRIPTION
        (
         "Interpolates *scat_data_raw* by frequency to give *scat_data_mono*.\n"
         ),
        AUTHORS( "Cory Davis" ),
        OUT( "scat_data_mono" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "scat_data_raw", "f_grid", "f_index" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "scat_data_rawCheck" ),
        DESCRIPTION
        (
         "Method for checking the consistency of the optical properties\n"
         "in the database.\n"
         "\n"
         "This function can be used to check datafiles containing data for\n"
         "randomly oriented scattering media.\n"
         "It is checked whether the data is consistent. The integral over\n"
         "the phase matrix should result the scattering cross section\n"
         "<C_sca>.\n"
         "\n"
         "The check is if:\n"
         "<C_ext> - <C_sca> = <C_abs>\n"
         "\n"
         "The result is printed on the screen.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "scat_data_raw" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Select" ),
        DESCRIPTION
        (
         "Method to select some elements from one array and copy them to\n"
         "a new array. (Works also for vectors.)\n"
         "\n"
         "This works also for higher dimensional objects, where the selection is\n"
         "always performed in the first dimension.\n"
         "\n"
         "For example:\n"
         "\n"
         "Select(y,x,[0,3])\n"
         "\n"
         "will select the first and fourth row of matrix x and copy them to the\n"
         "output matrix y.\n"
         "\n"
         "Note that it is even safe to use this method if needles and haystack\n"
         "are the same variable.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(      "needles" ),
        GOUT_TYPE( ARRAY_GROUPS + ", Vector, Matrix, Sparse" ),
        GOUT_DESC( "Selected elements. Must have the same variable type as "
                   "haystack." ),
        IN(),
        GIN(       "haystack", "needleindexes" ),
        GIN_TYPE(  ARRAY_GROUPS + ", Vector, Matrix, Sparse",
                   "ArrayOfIndex" ),
        GIN_DEFAULT( NODEF, NODEF ),
        GIN_DESC( "Variable to select from. May be the same variable as needles.",
                  "The elements to select (zero based indexing, as always.)" ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  )
        ));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensorOff" ),
        DESCRIPTION
        (
         "Sets sensor WSVs to obtain monochromatic pencil beam values.\n"
         "\n"
         "A 1D antenna pattern is assumed. The variables are set as follows:\n"
         "   antenna_dim             : 1.\n"
         "   mblock_za_grid          : Length 1, value 0.\n"
         "   mblock_aa_grid          : Empty.\n"
         "   sensor_response*        : As returned by *sensor_responseInit*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "sensor_response", "sensor_response_f", 
             "sensor_response_pol", "sensor_response_za",
             "sensor_response_aa", 
             "sensor_response_f_grid", "sensor_response_pol_grid",
             "sensor_response_za_grid", "sensor_response_aa_grid",
             "antenna_dim", "mblock_za_grid", "mblock_aa_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim", "stokes_dim", "f_grid" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_responseAntenna" ),
        DESCRIPTION
        (
         "Includes response of the antenna.\n"
         "\n"
         "The function returns the sensor response matrix after the antenna\n" 
         "characteristics have been included.\n"
         "\n"
         "The function handles \"multi-beam\" cases where the polarisation\n"
         "coordinate system is the same for all beams.\n"
         "\n"         
         "See *antenna_dim*, *antenna_los* and *antenna_response* for\n"
         "details on how to specify the antenna response.\n"
         ),
        AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
        OUT( "sensor_response", "sensor_response_f", "sensor_response_pol",
             "sensor_response_za", "sensor_response_aa", 
             "sensor_response_za_grid", "sensor_response_aa_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "sensor_response", "sensor_response_f", "sensor_response_pol",
            "sensor_response_za", "sensor_response_aa", "sensor_response_f_grid",
            "sensor_response_pol_grid", "sensor_response_za_grid",
            "sensor_response_aa_grid", "atmosphere_dim", "antenna_dim", 
            "antenna_los", "antenna_response", "sensor_norm" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_responseBackend" ),
        DESCRIPTION
        (
         "Includes response of the backend (spectrometer).\n"
         "\n"
         "The function returns the sensor response matrix after the backend\n" 
         "characteristics have been included.\n"
         "\n"
         "See *f_backend*, *backend_channel_response* and *sensor_norm* for\n"
         "details on how to specify the backend response.\n"
         ),
        AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
        OUT( "sensor_response", "sensor_response_f", "sensor_response_pol",
             "sensor_response_za", "sensor_response_aa", 
             "sensor_response_f_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "sensor_response", "sensor_response_f", "sensor_response_pol",
            "sensor_response_za", "sensor_response_aa", 
            "sensor_response_f_grid", "sensor_response_pol_grid", 
            "sensor_response_za_grid", "sensor_response_aa_grid",
            "f_backend", "backend_channel_response", "sensor_norm" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_responseFillFgrid" ),
        DESCRIPTION
        (
         "Polynomial frequency interpolation of spectra.\n"
         "z\n"
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
         "*polyorder*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "sensor_response", "sensor_response_f", "sensor_response_pol",
             "sensor_response_za", "sensor_response_aa", 
             "sensor_response_f_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "sensor_response", "sensor_response_f", "sensor_response_pol",
            "sensor_response_za", "sensor_response_aa", 
            "sensor_response_f_grid", "sensor_response_pol_grid", 
            "sensor_response_za_grid", "sensor_response_aa_grid" ),
        GIN( "polyorder", "nfill" ),
        GIN_TYPE( "Index", "Index" ),
        GIN_DEFAULT( "3", "2" ),
        GIN_DESC( "Polynomial order of interpolation", 
                  "Number of points to insert in each gap of f_grid" )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_responseBackendFrequencySwitching" ),
        DESCRIPTION
        (
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
         "*sensor_responseBackend*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "sensor_response", "sensor_response_f", "sensor_response_pol",
             "sensor_response_za", "sensor_response_aa", 
             "sensor_response_f_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "sensor_response", "sensor_response_f", "sensor_response_pol",
            "sensor_response_za", "sensor_response_aa", 
            "sensor_response_f_grid", "sensor_response_pol_grid", 
            "sensor_response_za_grid", "sensor_response_aa_grid",
            "f_backend", "backend_channel_response", "sensor_norm" ),
        GIN(    "df_1", "df2" ),
        GIN_TYPE(   "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF, NODEF ),
        GIN_DESC( "Frequency throw for cycle1.", "Frequency throw for cycle2.")
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_responseBeamSwitching" ),
        DESCRIPTION
        (
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
         "values of the second direction.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "sensor_response", "sensor_response_f", "sensor_response_pol",
             "sensor_response_za", "sensor_response_aa", 
             "sensor_response_za_grid", "sensor_response_aa_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "sensor_response", "sensor_response_f", "sensor_response_pol",
            "sensor_response_za", "sensor_response_aa",
            "sensor_response_f_grid", "sensor_response_pol_grid", 
            "sensor_response_za_grid", "sensor_response_aa_grid" ),
        GIN( "w1", "w2" ),
        GIN_TYPE( "Numeric", "Numeric" ),
        GIN_DEFAULT( "-1", "1" ),
        GIN_DESC( "Weight for values from first viewing direction.", 
                  "Weight for values from second viewing direction." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_responseFrequencySwitching" ),
        DESCRIPTION
        (
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
         "Output frequency grids are taken from the second spectrum..\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "sensor_response", "sensor_response_f", "sensor_response_pol",
             "sensor_response_za", "sensor_response_aa", 
             "sensor_response_f_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "sensor_response", "sensor_response_f", "sensor_response_pol",
            "sensor_response_za", "sensor_response_aa", 
            "sensor_response_f_grid", "sensor_response_pol_grid", 
            "sensor_response_za_grid", "sensor_response_aa_grid" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_responseFromArrayData" ),
        DESCRIPTION
        (
         "Sets up *sensor_response_array* from an existing *sensor_response*.\n"
         "\n"
         "Fills *sensor_response_array* and associated variables with\n"
         "corresponding non-array data. Hence, the array variables get all a\n"
         "length of 1.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "sensor_response", "sensor_response_f", "sensor_response_pol", 
             "sensor_response_za", "sensor_response_aa" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "mblock_index", "sensor_response_array", "sensor_response_f_array", 
            "sensor_response_pol_array", "sensor_response_za_array", 
            "sensor_response_aa_array", "sensor_response_index" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_responseIF2RF" ),
        DESCRIPTION
        (
         "Converts sensor response variables from IF to RF.\n"
         "\n"
         "The function converts intermediate frequencies (IF) in\n"
         "*sensor_response_f* and *sensor_response_f_grid* to radio\n"
         "frequencies (RF). This conversion is needed if the frequency\n"
         "translation of a mixer is included and the position of backend\n"
         "channels are specified in RF.\n"
         "\n"
         "A direct frequency conversion is performed. Values are not\n"
         "sorted in any way.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "sensor_response_f", "sensor_response_f_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "sensor_response_f", "sensor_response_f_grid", "lo", 
            "sideband_mode" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_responseInit" ),
        DESCRIPTION
        (
         "Initialises the variables summarising the sensor response.\n"
         "\n"
         "This method sets the variables to match monochromatic pencil beam\n"
         "calculations, to be further modified by inclusion of sensor\n"
         "characteristics. Use *sensorOff* if pure monochromatic pencil\n"
         "beam calculations shall be performed.\n"
         "\n"
         "The variables are set as follows:\n"
         "   sensor_response : Identity matrix, with size matching *f_grid*,\n"
         "                     *stokes_dim* *mblock_za_grid* and\n"
         "                     *mblock_aa_grid*.\n"
         "   sensor_response_f       : Repeated values of *f_grid*.\n"
         "   sensor_response_pol     : Data matching *stokes_dim*.\n"
         "   sensor_response_za      : Repeated values of *mblock_za_grid*.\n"
         "   sensor_response_aa      : Repeated values of *mblock_aa_grid*.\n"
         "   sensor_response_f_grid  : Equal to *f_grid*.\n"
         "   sensor_response_pol_grid: Set to 1:*stokes_dim*.\n"
         "   sensor_response_za_grid : Equal to *mblock_za_grid*.\n"
         "   sensor_response_aa_grid : Equal to *mblock_aa_grid*.\n"
         ),
        AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
        OUT( "sensor_response", "sensor_response_f", "sensor_response_pol", 
             "sensor_response_za", "sensor_response_aa", 
             "sensor_response_f_grid", "sensor_response_pol_grid",
             "sensor_response_za_grid", "sensor_response_aa_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "f_grid", "mblock_za_grid", "mblock_aa_grid", "antenna_dim",
            "atmosphere_dim", "stokes_dim", "sensor_norm" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_responseMixer" ),
        DESCRIPTION
        (
         "Includes response of the mixer of a heterodyne system.\n"
         "\n"
         "The function returns the sensor response matrix after the mixer\n" 
         "characteristics have been included. Frequency variables are\n"
         "converted from radio frequency (RF) to intermediate frequency (IF).\n"
         "The returned frequency grid covers the range [0,max_if], where\n"
         "max_if is the highest IF covered by the sideband response grid.\n" 
         "\n"
         "See *lo* and *sideband_response* for details on how to specify the\n"
         "mixer response\n"
         ),
        AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
        OUT( "sensor_response", "sensor_response_f", "sensor_response_pol",
             "sensor_response_za", "sensor_response_aa", 
             "sensor_response_f_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "sensor_response", "sensor_response_f", "sensor_response_pol",
            "sensor_response_za", "sensor_response_aa", "sensor_response_f_grid",
            "sensor_response_pol_grid", "sensor_response_za_grid",
            "sensor_response_aa_grid", "lo", "sideband_response", "sensor_norm" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_responseMultiMixerBackend" ),
        DESCRIPTION
        (
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
         "will be in absolute frequency (RF).\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "sensor_response", "sensor_response_f", "sensor_response_pol",
             "sensor_response_za", "sensor_response_aa", 
             "sensor_response_f_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "sensor_response", "sensor_response_f", "sensor_response_pol",
            "sensor_response_za", "sensor_response_aa", 
            "sensor_response_f_grid", "sensor_response_pol_grid", 
            "sensor_response_za_grid", "sensor_response_aa_grid",
            "lo_multi", "sideband_response_multi", 
            "sideband_mode_multi", "f_backend_multi",
            "backend_channel_response_multi", "sensor_norm" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_responsePolarisation" ),
        DESCRIPTION
        (
         "Extraction of non-default polarisation components.\n"
         "\n"
         "The default is to output the Stokes elements I, Q, U and V (up to\n" 
         "*stokes_dim*). This method allows to change the \"polarisation\" of\n"
         "the output. Polarisation components to be extracted are selected by\n"
         "*sensor_pol*. This method can be applied at any step of the sensor\n"
         "matrix set-up.\n"
         "\n"
         "The method can only be applied on data for I, Q, U and V. The value\n"
         "of *stokes_dim* must be sufficiently large for the selected\n"
         "components. For example, I+45 requires that *stokes_dim* is at\n"
         "least 3. \n"
         "\n"
         "See *sensor_pol* for coding of polarisation states.\n"
         "\n"
         "Note that the state of *y_unit* is considered. This WSV must give\n"
         "the actual unit of the data. This as, the extraction of components\n"
         "is slightly different if data are radiances or brightness\n"
         "temperatures.  In practise this means that *y_unit* (as to be\n"
         "applied inside *yCalc*) must be set before calling this method.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "sensor_response", "sensor_response_f", "sensor_response_pol",
             "sensor_response_za", "sensor_response_aa", 
             "sensor_response_pol_grid" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "sensor_response", "sensor_response_f", "sensor_response_pol",
            "sensor_response_za", "sensor_response_aa", "sensor_response_f_grid",
            "sensor_response_pol_grid", "sensor_response_za_grid",
            "sensor_response_aa_grid", "stokes_dim", "y_unit", "sensor_pol" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
  ( MdRecord
   ( NAME( "sensor_responseSimpleAMSU" ),
        DESCRIPTION
        (
         "Simplified sensor setup for an AMSU-type instrument.\n"
         "\n"
         "This method allows quick and simple definition of AMSU-type\n"
         "sensors. Assumptions:\n"
         "\n"
     "1. Pencil beam antenna.\n"
         "2. Douple sideband receivers.\n"
         "3. Sideband mode \"upper\"\n"
         "4. The channel response is rectangular.\n"
         "\n"
         "Under these assumptions the only inputs needed are the LO positions,\n"
         "the offsets from the LO, and the IF bandwidths. They are provieded\n"
         "in sensor_description_amsu.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "f_grid", 
         "antenna_dim", 
         "mblock_za_grid", 
         "mblock_aa_grid",
         "sensor_response", 
         "sensor_response_f", 
         "sensor_response_pol", 
         "sensor_response_za", 
         "sensor_response_aa", 
         "sensor_response_f_grid", 
         "sensor_response_pol_grid", 
         "sensor_response_za_grid", 
         "sensor_response_aa_grid", 
         "sensor_norm"
        ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim",
        "stokes_dim", 
        "sensor_description_amsu" ),
    GIN( "spacing" ),
    GIN_TYPE(    "Numeric" ),
    GIN_DEFAULT( ".1e9" ),
    GIN_DESC( "Desired grid spacing in Hz." )
        ));

        /* Not yet updated
     md_data_raw.push_back
     ( MdRecord
     ( NAME( "sensor_responsePolarisation" ),
     DESCRIPTION
     (
     "Adds polarisation to the response matrix.\n"
     "\n"
     "The output polarisations are given by matrix *sensor_pol*.\n"
     ),
     AUTHORS( "Mattias Ekstrom" ),
     OUT( "sensor_response", "sensor_response_pol" ),
     GOUT(),
     GOUT_TYPE(),
     GOUT_DESC(),
     IN( "sensor_pol", "sensor_response_za", "sensor_response_aa",
     "sensor_response_f", "stokes_dim" ),
     GIN(),
     GIN_TYPE(),
     GIN_DEFAULT(),
     GIN_DESC()
     ));
  */

  /* Not yet updated
     md_data_raw.push_back
     ( MdRecord
     ( NAME( "sensor_responseRotation" ),
     DESCRIPTION
     (
     "Adds rotation to the response matrix.\n"
     "\n"
     "The rotations are given by *sensor_rot* combined with *antenna_los*.\n"
     "The rotations are performed within each measurement block for the\n"
     "individual antennae.\n"
     "\n"
     "If used this method has to be run after the antenna response\n"
     "function and prior to sensor_responsePolarisation.\n"
     ),
     AUTHORS( "Mattias Ekstrom" ),
     OUT( "sensor_response" ),
     GOUT(),
     GOUT_TYPE(),
     GOUT_DESC(),
     IN( "sensor_rot", "antenna_los", "antenna_dim", "stokes_dim",
     "sensor_response_f", "sensor_response_za" ),
     GIN(),
     GIN_TYPE(),
     GIN_DEFAULT(),
     GIN_DESC()
     ));
  */

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_responseWMRF" ),
        DESCRIPTION
        (
         "Adds WMRF weights to sensor response.\n"
         "\n"
         "This method adds a spectrometer response that has been calculated\n"
         "with the weighted mean of representative frequencies (WMRF) method. It\n"
         "consists of a set of selected frequencies, and associated weights.\n"
         ),
        AUTHORS( "Stefan Buehler, based on Patrick Erikssons sensor_responseBackend" ),
        OUT( "sensor_response",
             "sensor_response_f",
             "sensor_response_pol",
             "sensor_response_za",
             "sensor_response_aa", 
             "sensor_response_f_grid"),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "sensor_response", "sensor_response_f", "sensor_response_pol",
            "sensor_response_za", "sensor_response_aa", 
            "sensor_response_f_grid", "sensor_response_pol_grid", 
            "sensor_response_za_grid", "sensor_response_aa_grid",
            "wmrf_weights",
            "f_backend" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_response_arraySingle" ),
        DESCRIPTION
        (
         "Sets up *sensor_response_array* from an existing *sensor_response*.\n"
         "\n"
         "Fills *sensor_response_array* and associated variables with\n"
         "corresponding non-array data. Hence, the array variables get all a\n"
         "length of 1.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "sensor_response_array", "sensor_response_f_array",
             "sensor_response_pol_array", "sensor_response_za_array",
             "sensor_response_aa_array" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "sensor_response", "sensor_response_f", "sensor_response_pol",
            "sensor_response_za", "sensor_response_aa" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "SparseSparseMultiply" ),
        DESCRIPTION
        (
         "Multiplies a Sparse with another Sparse, result stored in Sparse.\n"
         "\n"
         "Makes the calculation gout: = gin1 * gin2\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "gout"       ),
        GOUT_TYPE( "Sparse" ),
        GOUT_DESC( "Product, can be same variable as any of the inputs." ),
        IN(),
        GIN(      "gin1"      , "gin2"       ),
        GIN_TYPE(    "Sparse", "Sparse" ),
        GIN_DEFAULT( NODEF   , NODEF    ),
        GIN_DESC( "Left sparse matrix.",
                  "Right sparse matrix." )
        ));

  // This is duplicate with the 1-0 method tgsDefine. Merge!
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "SpeciesSet" ),
        DESCRIPTION
        (
         "Set up a list of absorption species tag groups.\n"
         "\n"
         "Workspace variables like *abs_species* contain several tag\n"
         "groups. Each tag group contains one or more tags. This method converts\n"
         "descriptions of tag groups given in the keyword to the ARTS internal\n"
         "representation (an *ArrayOfArrayOfSpeciesTag*). A tag group selects\n"
         "spectral features which belong to the same species.\n"
         "\n"
         "A tag is defined in terms of the name of the species, isotope, and a\n"
         "range of frequencies. Species are named after the standard chemical\n"
         "names, e.g., \"O3\". Isotopes are given by the last digit of the atomic\n"
         "weight, i.g., \"O3-668\" for the asymmetric ozone molecule including an\n"
         "oxygen 18 atom. Groups of transitions are specified by giving a lower\n"
         "and upper limit of a frequency range, e.g., \"O3-666-500e9-501e9\".\n"
         "\n"
         "To turn on Zeeman calculation for a Species, \"-Z\" may be appended\n"
         "to its name: \"O2-Z\" or \"O2-Z-66\"\n"
         "\n"
         "The symbol \"*\" acts as a wild card. Furthermore, frequency range or\n"
         "frequency range and isotope may be omitted.\n"
         "\n"
         "Finally, instead of the isotope the special letter \"nl\" may be given,\n"
         "e.g., \"H2O-nl\". This means that no absorption at all is associated\n"
         "with this tag. (It is not quite clear if this feature is useful for\n"
         "anything right now.)\n"
         "\n"
         "This method used to be a specific method for *abs_species*. Now it is\n"
         "generic, so that it can also be used to set *abs_nls* and *abs_pts*.\n"
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
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT(),
        GOUT(       "gout1"                         ),
        GOUT_TYPE( "ArrayOfArrayOfSpeciesTag" ),
        GOUT_DESC( "Output tag groups" ),
        IN(),
        GIN( "species" ),
        GIN_TYPE(    "ArrayOfString"   ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC("Specify one String for each tag group that you want to\n"
                 "create. Inside the String, separate the tags by commas\n"
                 "(plus optional blanks).\n")
        ));


  md_data_raw.push_back
    ( MdRecord
      ( NAME( "StringSet" ),
        DESCRIPTION
        (
         "Sets a String to the given text string.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "s"      ),
        GOUT_TYPE( "String" ),
        GOUT_DESC( "Variable to initialize." ),
        IN(),
        GIN(         "text"   ),
        GIN_TYPE(    "String" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Input text string." ),
        SETMETHOD( true )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "surfaceBlackbody" ),
        DESCRIPTION
        (
         "Creates variables to mimic a blackbody surface.\n"
         "\n"
         "This method sets up *surface_los*, *surface_rmatrix* and\n"
         "*surface_emission* for *surface_rtprop_agenda*. Here, *surface_los*\n"
         "and *surface_rmatrix* are set to be empty, and *surface_emission*\n"
         "to hold blackbody radiation for a temperature of *surface_skin_t*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "surface_los", "surface_rmatrix", "surface_emission" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "f_grid", "stokes_dim", "surface_skin_t", 
            "blackbody_radiation_agenda" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "surfaceFlatRefractiveIndex" ),
        DESCRIPTION
        (
         "Creates variables to mimic specular reflection by a (flat) surface\n"
         "where the refractive index is specified.\n"
         "\n"
         "The dielectric properties of the surface are described by\n"
         "*complex_n*. The Fresnel equations are used to calculate\n"
         "amplitude reflection coefficients. The method can thus result\n"
         "in that the reflection properties differ between frequencies\n"
         "and polarisations."
         "\n"
         "Local thermodynamic equilibrium is assumed, which corresponds to\n"
         "that the reflection and emission coefficients add up to 1.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "surface_los", "surface_rmatrix", "surface_emission" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "f_grid", "stokes_dim", "atmosphere_dim", "lat_grid", "lon_grid",
            "refellipsoid", "z_surface", "rte_pos", "rte_los", 
            "surface_skin_t", "complex_n", "blackbody_radiation_agenda" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "surfaceFlatReflectivity" ),
        DESCRIPTION
        (
         "Creates variables to mimic specular reflection by a (flat) surface\n"
         "where *surface_reflectivity* is specified.\n"
         "\n"
         "Works basically as *surfaceFlatScalarReflectivity* buit is more\n"
         "general as all also vector radiative transfer is calculated. See\n"
         "the ARTS theory document (ATD) for deatils around how\n"
         "*surface_emission* is determined. In the nomenclature of ATD,\n"
         "*surface_reflectivity* gives R.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "surface_los", "surface_rmatrix", "surface_emission" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "f_grid", "stokes_dim", "atmosphere_dim", "lat_grid", "lon_grid",
            "refellipsoid", "z_surface", "rte_pos", "rte_los", 
            "surface_skin_t", "surface_reflectivity", 
            "blackbody_radiation_agenda" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "surfaceFlatScalarReflectivity" ),
        DESCRIPTION
        (
         "Creates variables to mimic specular reflection by a (flat) surface\n"
         "where *surface_scalar_reflectivity* is specified.\n"
         "\n"
         "The method can only be used for *stokes_dim* equal to 1. Local\n"
         "thermodynamic equilibrium is assumed, which corresponds to that\n"
         "reflectivity and emissivity add up to 1.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "surface_los", "surface_rmatrix", "surface_emission" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "f_grid", "stokes_dim", "atmosphere_dim", "lat_grid", "lon_grid",
            "refellipsoid", "z_surface", "rte_pos", "rte_los", 
            "surface_skin_t", "surface_scalar_reflectivity",
            "blackbody_radiation_agenda" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "surfaceLambertianSimple" ),
        DESCRIPTION
        (
        "Creates variables to mimic a Lambertian surface.\n"
        "\n"
        "The method can only be used for 1D calculations.\n"        
        "\n"
        "A Lambertian surface can be characterised solely by its\n"
        "reflectivity, here taken from *surface_scalar_reflectivity*.\n"
        "\n"
        "The down-welling radiation field is estimated by making calculations\n"
        "for *np* directions. The range of zenith angles ([0,90]) is divided\n"
        "in an equidistant manner. The values for *surface_rmatrix* are\n"
        "assuming a constant radiance over each zenith angle range. See AUG.\n"
        "\n"
        "Default is to select the zenith angles for *sensor_los* to be placed\n"
        "centrally in the grid ranges. For example, if *np* is set to 9,\n"
        "down-welling radiation will be calculated for zenith angles = \n"
        "5, 15, ..., 85. The position of these angles can be shifted by\n"
        "*za_pos*. This variable specifies the fractional distance inside the\n"
        "ranges. For example, a *za_pos* of 0.7 (np still 9) gives the angles\n"
        "7, 17, ..., 87.\n"
        "\n"
        "Local thermodynamic equilibrium is assumed, which corresponds to\n"
        "that the reflection and emission coefficients \"add up to 1\".\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "surface_los", "surface_rmatrix", "surface_emission" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "f_grid", "stokes_dim", "atmosphere_dim", "rte_los", 
            "surface_skin_t", "surface_scalar_reflectivity",
            "blackbody_radiation_agenda" ),
        GIN(         "np",    "za_pos"  ),
        GIN_TYPE(    "Index", "Numeric" ),
        GIN_DEFAULT( NODEF,   "0.5"     ),
        GIN_DESC( "Number of zenith angles for calculation of down-welling " 
                  "radition.",
                  "Position of angle in *surface_los* inside ranges of zenith "
                  "angle grid. See above."
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "surface_reflectivityFromGriddedField6" ),
        DESCRIPTION
        (
         "Extracts surface reflectivities from a field of such data.\n"
         "\n"
         "This method allows to specify a field of surface reflectivity for\n"
         "automatic interpolation to points of interest. The position and\n"
         "direction for which the reflectivity shall be extracted are given\n"
         "by *rte_pos* and *rte_los*. The reflectivity field is expected to\n"
         "be stored as:\n"
         "   GriddedField4:\n"
         "      Vector f_grid[N_f]\n"
         "      Vector stokes_elements[N_s1]\n"
         "      Vector stokes_elements[N_s2]\n"
         "      Vector incidence_angle_grid[N_ia]\n"
         "      Vector lat_grid[N_lat]\n"
         "      Vector lon_grid[N_lon]\n"
         "      Tensor4 data[N_f][N_s1][N_s1][N_ia][N_lat][N_lon]\n"
         "\n"
         "Grids for incidence angle, latitude and longitude must have a\n"
         "length of >= 2 (ie. no automatic expansion). If the frequency grid\n"
         "has length 1, this is taken as the reflectivity is constant,\n"
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
         "   3. Linear interpolation if frequency (if input data have more\n"
         "      than one frequency).\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "surface_reflectivity" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "stokes_dim", "f_grid", "atmosphere_dim", "lat_grid", "lat_true", 
            "lon_true", "rte_pos", "rte_los" ),
        GIN( "r_field" ),
        GIN_TYPE( "GriddedField6" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "A field of surface reflectivities" )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "surface_scalar_reflectivityFromGriddedField4" ),
        DESCRIPTION
        (
         "Extracts scalar surface reflectivities from a field of such data.\n"
         "\n"
         "This method allows to specify a field of surface reflectivity for\n"
         "automatic interpolation to points of interest. The position and\n"
         "direction for which the reflectivity shall be extracted are given\n"
         "by *rte_pos* and *rte_los*. The reflectivity field is expected to\n"
         "be stored as:\n"
         "   GriddedField4:\n"
         "      Vector f_grid[N_f]\n"
         "      Vector incidence_angle_grid[N_ia]\n"
         "      Vector lat_grid[N_lat]\n"
         "      Vector lon_grid[N_lon]\n"
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
         "      than one frequency).\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "surface_scalar_reflectivity" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "stokes_dim", "f_grid", "atmosphere_dim", "lat_grid", "lat_true", 
            "lon_true", "rte_pos", "rte_los" ),
        GIN( "r_field" ),
        GIN_TYPE( "GriddedField4" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "A field of scalar surface reflectivities" )
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "TangentPointExtract" ),
        DESCRIPTION
        (
         "Finds the tangent point of a propagation path.\n"
         "\n"
         "The tangent point is here defined as the point with the lowest\n"
         "altitude (which differes from the definition used in the code\n"
         "where it is the point with the lowest radius, or equally the point\n"
         "with a zenith angle of 90 deg.)\n"
         "\n"
         "The tangent point is returned as a vector, with columns matching\n"
         "e.g. *rte_pos*. If the propagation path has no tangent point, the\n"
         "vector is set to NaN.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT( "tan_pos" ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "The position vector of the tangent point." ),
        IN( "ppath" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "TangentPointPrint" ),
        DESCRIPTION
        (
         "Prints information about the tangent point of a propagation path.\n"
         "\n"
         "The tangent point is here defined as the point with the lowest\n"
         "altitude (which differes from the definition used in the code\n"
         "where it is the point with the lowest radius, or equally the point\n"
         "with a zenith angle of 90 deg.)\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "ppath" ),
        GIN( "level" ),
        GIN_TYPE(    "Index" ),
        GIN_DEFAULT( "1" ),
        GIN_DESC( "Output level to use." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor3AddScalar" ),
        DESCRIPTION
        (
         "Adds a scalar value to all elements of a tensor3.\n"
         "\n"
         "The result can either be stored in the same or another\n"
         "variable.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "tout"    ),
        GOUT_TYPE( "Tensor3" ),
        GOUT_DESC( "Output tensor." ),
        IN(),
        GIN(         "tin",     "value"   ),
        GIN_TYPE(    "Tensor3", "Numeric" ),
        GIN_DEFAULT( NODEF    , NODEF     ),
        GIN_DESC( "Input tensor.",
                  "The value to be added to the tensor." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor3Scale" ),
        DESCRIPTION
        (
         "Scales all elements of a tensor with the specified value.\n"
         "\n"
         "The result can either be stored in the same or another\n"
         "variable.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUT(),
        GOUT(      "tout"    ),
        GOUT_TYPE( "Tensor3" ),
        GOUT_DESC( "Output tensor." ),
        IN(),
        GIN(         "tin",     "value"   ),
        GIN_TYPE(    "Tensor3", "Numeric" ),
        GIN_DEFAULT( NODEF    , NODEF     ),
        GIN_DESC( "Input tensor.",
                  "The value to be multiplied with the tensor." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor3SetConstant" ),
        DESCRIPTION
        (
         "Creates a tensor and sets all elements to the specified value.\n"
         "\n"
         "The size is determined by *ncols*, *nrows* etc.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT(),
        GOUT(      "t"   ),
        GOUT_TYPE( "Tensor3" ),
        GOUT_DESC( "Variable to initialize." ),
        IN( "npages", "nrows", "ncols" ),
        GIN(         "value"   ),
        GIN_TYPE(    "Numeric" ),
        GIN_DEFAULT( NODEF     ),
        GIN_DESC( "Tensor value." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor4Scale" ),
        DESCRIPTION
        (
         "Scales all elements of a tensor with the specified value.\n"
         "\n"
         "The result can either be stored in the same or another\n"
         "variable.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUT(),
        GOUT(      "tout"    ),
        GOUT_TYPE( "Tensor4" ),
        GOUT_DESC( "Output tensor." ),
        IN(),
        GIN(         "tin",     "value"   ),
        GIN_TYPE(    "Tensor4", "Numeric" ),
        GIN_DEFAULT( NODEF    , NODEF     ),
        GIN_DESC( "Input tensor.",
                  "The value to be multiplied with the tensor." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor4SetConstant" ),
        DESCRIPTION
        (
         "Creates a tensor and sets all elements to the specified value.\n"
         "\n"
         "The size is determined by *ncols*, *nrows* etc.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT(),
        GOUT(      "t"       ),
        GOUT_TYPE( "Tensor4" ),
        GOUT_DESC( "Variable to initialize." ),
        IN( "nbooks", "npages", "nrows", "ncols" ),
        GIN(         "value"   ),
        GIN_TYPE(    "Numeric" ),
        GIN_DEFAULT( NODEF     ),
        GIN_DESC( "Tensor value." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor5Scale" ),
        DESCRIPTION
        (
         "Scales all elements of a tensor with the specified value.\n"
         "\n"
         "The result can either be stored in the same or another\n"
         "variable.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUT(),
        GOUT(      "tout"    ),
        GOUT_TYPE( "Tensor5" ),
        GOUT_DESC( "Output tensor." ),
        IN(),
        GIN(         "tin",     "value"   ),
        GIN_TYPE(    "Tensor5", "Numeric" ),
        GIN_DEFAULT( NODEF    , NODEF     ),
        GIN_DESC( "Input tensor.",
                  "The value to be multiplied with the tensor." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor5SetConstant" ),
        DESCRIPTION
        (
         "Creates a tensor and sets all elements to the specified value.\n"
         "\n"
         "The size is determined by *ncols*, *nrows* etc.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT(),
        GOUT(      "t"       ),
        GOUT_TYPE( "Tensor5" ),
        GOUT_DESC( "Variable to initialize." ),
        IN( "nshelves", "nbooks", "npages", "nrows", "ncols" ),
        GIN(         "value"   ),
        GIN_TYPE(    "Numeric" ),
        GIN_DEFAULT( NODEF     ),
        GIN_DESC( "Tensor value." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor6Scale" ),
        DESCRIPTION
        (
         "Scales all elements of a tensor with the specified value.\n"
         "\n"
         "The result can either be stored in the same or another\n"
         "variable.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUT(),
        GOUT(      "tout"    ),
        GOUT_TYPE( "Tensor6" ),
        GOUT_DESC( "Output tensor." ),
        IN(),
        GIN(         "tin",     "value"   ),
        GIN_TYPE(    "Tensor6", "Numeric" ),
        GIN_DEFAULT( NODEF    , NODEF     ),
        GIN_DESC( "Input tensor.",
                  "The value to be multiplied with the tensor." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor6SetConstant" ),
        DESCRIPTION
        (
         "Creates a tensor and sets all elements to the specified value.\n"
         "\n"
         "The size is determined by *ncols*, *nrows* etc.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT(),
        GOUT(      "t"       ),
        GOUT_TYPE( "Tensor6" ),
        GOUT_DESC( "Variable to initialize." ),
        IN( "nvitrines", "nshelves", "nbooks", "npages", "nrows", "ncols" ),
        GIN(         "value"   ),
        GIN_TYPE(    "Numeric" ),
        GIN_DEFAULT( NODEF     ),
        GIN_DESC( "Tensor value." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor7Scale" ),
        DESCRIPTION
        (
         "Scales all elements of a tensor with the specified value.\n"
         "\n"
         "The result can either be stored in the same or another\n"
         "variable.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUT(),
        GOUT(      "tout"    ),
        GOUT_TYPE( "Tensor7" ),
        GOUT_DESC( "Output tensor." ),
        IN(),
        GIN(         "tin",     "value"   ),
        GIN_TYPE(    "Tensor7", "Numeric" ),
        GIN_DEFAULT( NODEF    , NODEF     ),
        GIN_DESC( "Input tensor.",
                  "The value to be multiplied with the tensor." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor7SetConstant" ),
        DESCRIPTION
        (
         "Creates a tensor and sets all elements to the specified value.\n"
         "\n"
         "The size is determined by *ncols*, *nrows* etc.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT(),
        GOUT(      "t"       ),
        GOUT_TYPE( "Tensor7" ),
        GOUT_DESC( "Variable to initialize." ),
        IN( "nlibraries", "nvitrines", "nshelves", "nbooks", "npages", "nrows",
            "ncols" ),
        GIN(         "value"   ),
        GIN_TYPE(    "Numeric" ),
        GIN_DEFAULT( NODEF     ),
        GIN_DESC( "Tensor value." )
        ));

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
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_mat_per_speciesAddZeeman" ),
        DESCRIPTION
        (
        "This function will, for each Zeeman species, make a local\n"
        "ArrayOfLineRecord for the various transition types with Zeeman\n"
        "altered LineRecord(s).  These are then composed into a single\n"
        "ArrayOfArrayOfLineRecord which is processed as per the scalar case.\n"
        "\n"
        "The line broadened absorption coefficients are finally multiplied with\n"
        "the transition type rotation matrix and the new variable is inserted into\n"
        "the out variable. Only -Z- species are treated.\n"
        "\n"
        "Note that between 55 GHz and 65 GHz there is usually ~700 O_2 lines,\n"
        "however, when this Zeeman splitting method is used, the number of\n"
        "lines is increased to about 45,000. This is a time consuming method.\n"
        "\n"
        "If rte_mag is of length 1 and contains -1.0 scalar calculations will follow.\n"
         ),
        AUTHORS( "Richard Larsson" ),
        OUT("abs_mat_per_species"),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN("abs_mat_per_species",
           "f_grid",
           "abs_species",
           "abs_n2",
           "abs_lines_per_species",
           "abs_lineshape",
           "abs_cont_names",
           "abs_cont_models",
           "abs_cont_parameters",
           "f_index",
           "rte_pressure", "rte_temperature", "rte_vmr_list", "rte_doppler",
           "rte_los", "rte_mag"),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "timerStart" ),
        DESCRIPTION
        (
         "Initializes the CPU timer."
         "\n"
         "Use *timerStop* to output the consumed cpu time since *timerStart*.\n"
         "\n"
         "Usage example:\n"
         "   timerStart\n"
         "   ReadXML(f_grid,\"frequencies.xml\")\n"
         "   timerStop\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "timer" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "timerStop" ),
        DESCRIPTION
        (
         "Stops the CPU timer."
         "\n"
         "Use *timerStop* to output the consumed cpu time since *timerStart*.\n"
         "See *timerStart* for example usage.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "timer" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Touch" ),
        DESCRIPTION
        (
         "As *Ignore* but for agenda output.\n"
         "\n"
         "This method is handy for use in agendas in order to suppress\n"
         "warnings about unused output workspace variables. What it does is:\n"
         "Nothing!\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(      "gout1"    ),
        GOUT_TYPE( "Any" ),
        GOUT_DESC( "Variable to do nothing with." ),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC(),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorAddScalar" ),
        DESCRIPTION
        (
         "Adds a scalar to all elements of a vector.\n"
         "\n"
         "The result can either be stored in the same or another vector.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "v1"     ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "Input vector" ),
        IN(),
        GIN(         "v2"    , "value"   ),
        GIN_TYPE(    "Vector", "Numeric" ),
        GIN_DEFAULT( NODEF   , NODEF     ),
        GIN_DESC( "Output vector", "The value to be added to the vector." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorCompare" ),
        DESCRIPTION
        (
         "Checks the consistency between two vectors.\n" 
         "\n"
         "The two vectors are checked to not deviate outside the specified\n"
         "value (*maxabsdiff*). An error is issued if this is not fulfilled.\n"
         "\n"
         "The main application of this method is to be part of the test\n"
         "control files, and then used to check that a calculated spectrum\n"
         "is consistent with an old, reference, calculation.\n"
         "\n"
         "The default value for *maxabsdiff* is adopted for comparing two\n"
         "spectra with brightness temperature as unit.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN( "vector1", "vector2", "maxabsdiff", "error_message" ),
        GIN_TYPE( "Vector", "Vector", "Numeric", "String" ),
        GIN_DEFAULT( NODEF, NODEF, "0.01", "" ),
        GIN_DESC( "A first vector", "A second vector", 
                  "Threshold for maximum absolute difference.",
                  "Additional error message.")
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorExtractFromMatrix" ),
        DESCRIPTION
        (
         "Extract a Vector from a Matrix.\n"
         "\n"
         "Copies row or column with given Index from input Matrix variable\n"
         "to create output Vector.\n"
         ),
        AUTHORS( "Patrick Eriksson, Oliver Lemke, Stefan Buehler" ),
        OUT(),
        GOUT(      "v"      ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "Extracted vector." ),
        IN(),
        GIN(          "m"     , "i"    , "direction" ),
        GIN_TYPE(     "Matrix", "Index", "String"    ),
        GIN_DEFAULT(  NODEF   , NODEF  , NODEF       ),
        GIN_DESC( "Input matrix.",
                  "Index of row or column.",
                  "Direction. \"row\" or \"column\"." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorFlip" ),
        DESCRIPTION
        (
         "Flips a vector.\n"
         "\n"
         "The output is the input vector in reversed order. The result can\n"
         "either be stored in the same or another vector.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "gout1"       ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "Output vector." ),
        IN(),
        GIN(      "gin1"     ),
        GIN_TYPE( "Vector" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Input vector." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorInsertGridPoints" ),
        DESCRIPTION
        (
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
         "  Vector : The points to insert.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT(),
        GOUT(      "gout1"       ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "The new grid vector" ),
        IN(),
        GIN(       "gin1"      , "gin2"       ),
        GIN_TYPE(     "Vector", "Vector" ),
        GIN_DEFAULT(  NODEF   , NODEF    ),
        GIN_DESC( "The original grid vector",
                  "The points to insert" )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorLinSpace" ),
        DESCRIPTION
        (
         "Creates a vector with linear spacing.\n"
         "\n"
         "The first element equals always the start value, and the spacing\n"
         "equals always the step value, but the last value can deviate from\n"
         "the stop value. *step* can be both positive and negative.\n"
         "\n"
         "The created vector is [start, start+step, start+2*step, ...]\n "
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "v"      ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "Output vector." ),
        IN(),
        GIN(         "start",   "stop",    "step"    ),
        GIN_TYPE(    "Numeric", "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF,     NODEF,     NODEF     ),
        GIN_DESC( "Start value.",
                  "Maximum/minimum value of the end value",
                  "Spacing of the vector." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorLogSpace" ),
        DESCRIPTION
        (
         "Creates a vector with logarithmic spacing.\n"
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
         " exp([ln(start), ln(start)+step, ln(start)+2*step, ...])\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT(),
        GOUT(      "gout1"       ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "Variable to initialize." ),
        IN(),
        GIN( "start",   "stop",    "step"    ),
        GIN_TYPE(    "Numeric", "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF,     NODEF,     NODEF ),
        GIN_DESC( "The start value. (Direct coordinates!)",
                  "The maximum value of the end value. (Direct coordinates!)",
                  "The spacing of the vector. (Log coordinates!)" )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorMatrixMultiply" ),
        DESCRIPTION
        (
         "Multiply a Vector with a Matrix and store the result in another\n"
         "Vector.\n"
         "\n"
         "This just computes the normal Matrix-Vector product, y=M*x. It is ok\n"
         "if input and output Vector are the same. This function is handy for\n"
         "multiplying the H Matrix to spectra.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT(),
        GOUT(      "gout1"       ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "The result of the multiplication (dimension m)." ),
        IN(),
        GIN(       "gin1"      , "gin2"       ),
        GIN_TYPE(     "Matrix", "Vector" ),
        GIN_DEFAULT(  NODEF   , NODEF    ),
        GIN_DESC( "The Matrix to multiply (dimension mxn).",
                  "The original Vector (dimension n)." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorNLinSpace" ),
        DESCRIPTION
        (
         "Creates a vector with length *nelem*, equally spaced between the\n"
         "given end values.\n"
         "\n"
         "The length (*nelem*) must be larger than 1.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "v"      ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "Variable to initialize." ),
        IN( "nelem" ),
        GIN(         "start",   "stop"    ),
        GIN_TYPE(    "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF,     NODEF     ),
        GIN_DESC( "Start value.",
                  "End value." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorNLogSpace" ),
        DESCRIPTION
        (
         "Creates a vector with length *nelem*, equally logarithmically\n"
         "spaced between the given end values.\n"
         "\n"
         "The length (*nelem*) must be larger than 1.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "v"      ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "Variable to initialize." ),
        IN( "nelem" ),
        GIN(         "start",   "stop"    ),
        GIN_TYPE(    "Numeric", "Numeric" ),
        GIN_DEFAULT( NODEF,     NODEF     ),
        GIN_DESC( "Start value.",
                  "End value." 
                  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorScale" ),
        DESCRIPTION
        (
         "Scales all elements of a vector with the same value.\n"
         "\n"
         "The result can either be stored in the same or another vector.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "gout1"       ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "Output vector." ),
        IN(),
        GIN(      "gin1"      ,
                  "value" ),
        GIN_TYPE(    "Vector",
                     "Numeric" ),
        GIN_DEFAULT( NODEF   ,
                     NODEF ),
        GIN_DESC( "Input vector.",
                  "Scaling value." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorSetConstant" ),
        DESCRIPTION
        (
         "Creates a vector and sets all elements to the specified value.\n"
         "\n"
         "The vector length is determined by *nelem*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(      "v"      ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "Variable to initialize." ),
        IN( "nelem" ),
        GIN(         "value"   ),
        GIN_TYPE(    "Numeric" ),
        GIN_DEFAULT( NODEF     ),
        GIN_DESC( "Vector value." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorSet" ),
        DESCRIPTION
        (
         "Create a vector from the given list of numbers.\n"
         "\n"
         "   VectorSet(p_grid, [1000, 100, 10] )\n"
         "   Will create a p_grid vector with these three elements.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT(),
        GOUT(      "gout1"       ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "Variable to initialize." ),
        IN(),
        GIN( "values"   ),
        GIN_TYPE(    "Vector" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "The vector elements." ),
        SETMETHOD( true )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorZtanToZaRefr1D" ),
        DESCRIPTION
        (
         "Converts a set of true tangent altitudes to zenith angles.\n"
         "\n"
         "The tangent altitudes are given to the function as a vector, which\n"
         "are converted to a generic vector of zenith angles. The position of\n"
         "the sensor is given by the WSV *sensor_pos*. The function works\n"
         "only for 1D. The zenith angles are always set to be positive.\n"
         ),
        AUTHORS( "Patrick Eriksson", "Mattias Ekstrom" ),
        OUT(),
        GOUT(      "v_za"       ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "Vector with zenith angles." ),
        IN( "refr_index_agenda", "sensor_pos", "p_grid", "t_field", "z_field",
            "vmr_field", "edensity_field", "refellipsoid", "atmosphere_dim", 
            "f_index" ),
        GIN(         "v_ztan" ),
        GIN_TYPE(    "Vector" ),
        GIN_DEFAULT( NODEF    ),
        GIN_DESC( "Vector with tangent altitudes." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorZtanToZa1D" ),
        DESCRIPTION
        (
         "Converts a set of geometrical tangent altitudes to zenith angles.\n"
         "\n"
         "The tangent altitudes are given to the function as a vector, which\n"
         "are converted to a generic vector of zenith angles. The position of\n"
         "the sensor is given by the WSV *sensor_pos*. The function works\n"
         "only for 1D. The zenith angles are always set to be positive.\n"
         ),
        AUTHORS( "Patrick Eriksson", "Mattias Ekstrom" ),
        OUT(),
        GOUT(      "v_za"       ),
        GOUT_TYPE( "Vector" ),
        GOUT_DESC( "Vector with zenith angles." ),
        IN( "sensor_pos", "refellipsoid", "atmosphere_dim" ),
        GIN(         "v_ztan" ),
        GIN_TYPE(    "Vector" ),
        GIN_DEFAULT( NODEF    ),
        GIN_DESC( "Vector with tangent altitudes." )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "verbosityInit" ),
        DESCRIPTION
        (
         "Initializes the verbosity levels.\n"
         "\n"
         "Sets verbosity to defaults or the levels specified by -r on the command line.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "verbosity" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "verbositySet" ),
        DESCRIPTION
        (
         "Sets the verbosity levels.\n"
         "\n"
         "Sets the reporting level for agenda calls, screen and file.\n"
         "All reporting levels can reach from 0 (only error messages)\n"
         "to 3 (everything). The agenda setting applies in addition\n"
         "to both screen and file output.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "verbosity" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(         "agenda", "screen", "file" ),
        GIN_TYPE(    "Index",  "Index",  "Index" ),
        GIN_DEFAULT( NODEF,    NODEF,    NODEF),
        GIN_DESC(    "Agenda verbosity level",
                     "Screen verbosity level",
                     "Report file verbosity level")
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "verbositySetAgenda" ),
        DESCRIPTION
        (
         "Sets the verbosity level for agenda output.\n"
         "\n"
         "See *verbositySet*\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "verbosity" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "verbosity" ),
        GIN(         "level" ),
        GIN_TYPE(    "Index" ),
        GIN_DEFAULT( NODEF),
        GIN_DESC(    "Agenda verbosity level")
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "verbositySetFile" ),
        DESCRIPTION
        (
         "Sets the verbosity level for report file output.\n"
         "\n"
         "See *verbositySet*\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "verbosity" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "verbosity" ),
        GIN(         "level" ),
        GIN_TYPE(    "Index" ),
        GIN_DEFAULT( NODEF),
        GIN_DESC(    "Report file verbosity level")
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "verbositySetScreen" ),
        DESCRIPTION
        (
         "Sets the verbosity level for screen output.\n"
         "\n"
         "See *verbositySet*\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT( "verbosity" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "verbosity" ),
        GIN(         "level" ),
        GIN_TYPE(    "Index" ),
        GIN_DEFAULT( NODEF),
        GIN_DESC(    "Screen verbosity level")
        ));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "WMRFSelectChannels" ),
        DESCRIPTION
        (
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
         "- Unnecessary frequencies are removed from wmrf_weights.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "f_grid", "wmrf_weights",
             "f_backend" ),
        GOUT(      ),
        GOUT_TYPE( ),
        GOUT_DESC(),
        IN( "f_grid", "f_backend", 
            "wmrf_weights", "wmrf_channels"  ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "WriteMolTau" ),
        DESCRIPTION
        (
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
         "\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN("f_grid", "z_field", "abs_mat_field", "atmosphere_dim" ),
        GIN("filename"),
        GIN_TYPE("String"),
        GIN_DEFAULT( NODEF),
        GIN_DESC("Name of the *molecular_tau_file*." )
        ));
  
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "WriteNetCDF" ),
        DESCRIPTION
        (
         "Writes a workspace variable to a NetCDF file.\n"
         "\n"
         "This method can write variables of any group.\n"
         "\n"
         "If the filename is omitted, the variable is written\n"
         "to <basename>.<variable_name>.nc.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN(),
        GIN(          "v",
                      "filename" ),
        GIN_TYPE(     "Vector, Matrix, Tensor3, Tensor4, Tensor5, ArrayOfVector,"
                      "ArrayOfMatrix, GasAbsLookup",
                      "String" ),
        GIN_DEFAULT(  NODEF,
                      "" ),
        GIN_DESC(     "Variable to be saved.",
                      "Name of the NetCDF file." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  ),
        PASSWORKSPACE(  false ),
        PASSWSVNAMES(   true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "WriteXML" ),
        DESCRIPTION
        (
         "Writes a workspace variable to an XML file.\n"
         "\n"
         "This method can write variables of any group.\n"
         "\n"
         "If the filename is omitted, the variable is written\n"
         "to <basename>.<variable_name>.xml.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "output_file_format" ),
        GIN(         "v",
                     "filename" ),
        GIN_TYPE(    "Any",
                     "String" ),
        GIN_DEFAULT( NODEF,
                     "" ),
        GIN_DESC(    "Variable to be saved.",
                     "Name of the XML file." ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  ),
        PASSWORKSPACE(  false ),
        PASSWSVNAMES(   true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "WriteXMLIndexed" ),
        DESCRIPTION
        (
         "As *WriteXML*, but creates indexed file names.\n"
         "\n"
         "The variable is written to a file with name:\n"
         "   <filename>.<file_index>.xml.\n"
         "where <file_index> is the value of *file_index*.\n"
         "\n"
         "This means that *filename* shall here not include the .xml\n"
         "extension. Omitting filename works as for *WriteXML*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT(),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "output_file_format", "file_index" ),
        GIN(          "wsv", "filename" ),
        GIN_TYPE(     "Any", "String"   ),
        GIN_DEFAULT(  NODEF, ""         ),
        GIN_DESC( "Workspace variable to be saved.",
                  "File name. See above." 
                  ),
        SETMETHOD(      false ),
        AGENDAMETHOD(   false ),
        USES_TEMPLATES( true  ),
        PASSWORKSPACE(  false ),
        PASSWSVNAMES(   true  )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "yApplyYunit" ),
        DESCRIPTION
        (
         "Conversion of *y* to other spectral units.\n"
         "\n"
         "Any conversion to brightness temperature is normally made inside\n"
         "*yCalc*. This method makes it possible to also make this conversion\n"
         "after *yCalc*, but with restrictions for *jacobian*.\n"
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
         "If you are using this method, *y_unit* should be set to \"1\" when\n"
         "calling *yCalc*, and be changed before calling this method.\n"
         "\n"         
         "See further *y_unit*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "y", "jacobian" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "y", "jacobian", "y_f", "y_pol", "y_unit" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "yCalc" ),
        DESCRIPTION
        (
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
         "The unit of output radiances and jacobians follow *y_unit*. The\n"
         "conversion is applied on monochromatic pencil beam values. That\n"
         "is, before any sensor responses have been included.\n"
         "The frequency, polarisation etc. for each measurement value is\n" 
         "given by *y_f*, *y_pol*, *y_pos* and *y_los*.\n"
         "\n"
         "See the method selected for *iy_main_agenda* for quantities\n"
         "that can be obtained by *y_aux*. However, in no case data of\n"
         "along-the-path type can be extracted. *y_unit* is applied on the\n"
         "following aux data:\n"
         "    \"iy\", \"Error\" and \"Error (uncorrelated)\"\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "y", "y_f", "y_pol", "y_pos", "y_los", "y_aux", "jacobian" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "basics_checked", "atmosphere_dim",
            "t_field", "z_field", "vmr_field", "cloudbox_on", 
            "cloudbox_checked", "stokes_dim", "f_grid", "sensor_pos", 
            "sensor_los", "mblock_za_grid", "mblock_aa_grid", "antenna_dim", 
            "sensor_response", "sensor_response_f",
            "sensor_response_pol", "sensor_response_za", "sensor_response_aa",
            "iy_main_agenda", "y_unit", 
            "jacobian_agenda", "jacobian_do", "jacobian_quantities",
            "jacobian_indices", "iy_aux_vars" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ybatchCalc" ),
        DESCRIPTION
        (
         "Performs batch calculations for the measurement vector y.\n"
         "\n"
         "We perform *ybatch_n* jobs, starting at index *ybatch_start*. (Zero\n"
         "based indexing, as usual.) The output matrix *ybatch* will have\n"
         "dimension (y.nelem(),ybatch_n). So, indices in the output matrix start\n"
         "with zero, independent of *ybatch_start*.\n"
         "\n"
         "The method performs the following:\n"
         "   1. Sets *ybatch_index* = *ybatch_start*.\n"
         "   2. Performs a-d until\n"
         "      *ybatch_index* = *ybatch_start* + *ybatch_n*.\n"
         "        a. Executes *ybatch_calc_agenda*.\n"
         "        b. If *ybatch_index* = *ybatch_start*, resizes *ybatch*\n"
         "           based on *ybatch_n* and length of *y*.\n"
         "        c. Copies *y* to column *ybatch_index* - *ybatch_start*\n"
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
         "See the user guide for further practical examples.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUT( "ybatch", "ybatch_jacobians" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "ybatch_start", "ybatch_n", "ybatch_calc_agenda" ), 
        GIN( "robust" ),
        GIN_TYPE(    "Index" ),
        GIN_DEFAULT( "0" ),
        GIN_DESC(
                 "A flag with value 1 or 0. If set to one, the batch\n"
                 "calculation will continue, even if individual jobs\n"
                 "fail. In that case, a warning message is written to\n"
                 "screen and file (out1 output stream), and ybatch for the\n"
                 "failed job is set to -1. The robust behavior does only work\n"
                 "properly if your control file is run single threaded.\n"
                 "Set \"--numthreads 1\". See \"arts --help\"."
                 )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ybatchMetProfiles" ),
        DESCRIPTION
        (
         "This method is used for simulating ARTS for metoffice model fields"
         "\n"
         "This method reads in *met_amsu_data* which contains the\n"
         "lat-lon of the metoffice profile files as a Matrix. It then\n"
         "loops over the number of profiles and corresponding to each\n"
         "longitude create the appropriate profile basename. Then,\n"
         "corresponding to each basename we have temperature field, altitude\n"
         "field, humidity field and particle number density field. The\n"
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
         "from the profiles inside the function\n"
         ),
        AUTHORS( "Sreerekha T.R." ),
        OUT( "ybatch" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species", "met_profile_calc_agenda", "f_grid", "met_amsu_data",
            "sensor_pos", "refellipsoid", "lat_grid", "lon_grid", 
            "atmosphere_dim", "scat_data_raw" ),
        GIN( "nelem_p_grid", "met_profile_path", "met_profile_pnd_path" ),
        GIN_TYPE(    "Index",        "String",           "String" ),
        GIN_DEFAULT( NODEF,          NODEF,              NODEF ),
        GIN_DESC( "FIXME DOC",
                  "FIXME DOC",
                  "FIXME DOC" )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ybatchMetProfilesClear" ),
        DESCRIPTION
        (
         "This method is used for simulating ARTS for metoffice model fields\n"
         "for clear sky conditions.\n"
         "\n"
         "This method reads in *met_amsu_data* which contains the\n"
         "lat-lon of the metoffice profile files as a Matrix. It then\n"
         "loops over the number of profiles and corresponding to each\n"
         "longitude create the appropriate profile basename. Then,\n"
         "Corresponding to each basename we have temperature field, altitude\n"
         "field, humidity field and particle number density field. The\n"
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
         "from the profiles inside the function\n"
         ),
        AUTHORS( "Seerekha T.R." ),
        OUT( "ybatch" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "abs_species", "met_profile_calc_agenda", 
            "f_grid", "met_amsu_data", "sensor_pos", "refellipsoid" ),
        GIN( "nelem_p_grid", "met_profile_path" ),
        GIN_TYPE(    "Index",        "String" ),
        GIN_DEFAULT( NODEF,          NODEF ),
        GIN_DESC( "FIXME DOC",
                  "FIXME DOC" )
        ));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "z_fieldFromHSE" ),
        DESCRIPTION
        (
         "Forse altitudes to fulfil hydrostatic equilibrium.\n"
         "\n"
         "The method applies hydrostatic equilibrium. A mixture of \"dry\n"
         "air\" and water vapour is assumed. That is, the air is assumed to\n"
         "be well mixed and its weight, beside water vapour, is constant\n"
         "(*molarmass_dry_air*). In addition, the effect of any particles\n"
         "(including liquid and ice particles) is neglected.\n"
         "\n"
         "The output is an update of *z_field*. This variable is expected to\n"
         "contain approximative altitudes when calling the function. The\n"
         "altitude matching *p_hse* is kept constant. Other altitudes are\n"
         "basically arbitrary, but good estimates give quicker calculations.\n"
         "\n"
         "The calculations are repeated until the change in altitude is below\n"
         "*z_hse_accuracy*. An iterative process is needed as gravity varies\n"
         "with altitude.\n"
         "\n"
         "For 1D and 2D, the geographical position is taken from *lat_true*\n"
         "and *lon_true*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( "z_field" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "atmosphere_dim", "p_grid", "lat_grid", "lon_grid", "lat_true", 
            "lon_true", "abs_species", "t_field", "z_field", "vmr_field", 
            "refellipsoid", "z_surface", "basics_checked", "g0_agenda",
            "molarmass_dry_air", "p_hse", "z_hse_accuracy" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
}

