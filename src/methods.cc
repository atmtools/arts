/* Copyright (C) 2000-2008
   Stefan Buehler <buehler@uni-bremen.de>
   Patrick Eriksson <patrick@rss.chalmers.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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
#include "auto_wsv.h"
#include "methods.h"
#include "auto_wsv_groups.h"

// Some #defines and typedefs to make the records better readable:
#define NAME(x) x
#define DESCRIPTION(x) x
#define AUTHORS  MakeArray<String>
#define OUTPUT   MakeArray<Index>
#define INPUT    MakeArray<Index>
#define GOUTPUT  MakeArray<Index>
#define GINPUT   MakeArray<Index>
#define KEYWORDS MakeArray<String>
#define DEFAULTS MakeArray<String>
#define TYPES    MakeArray<TokValType>
#define AGENDAMETHOD(x) x
#define SUPPRESSHEADER(x) x
#define PASSWORKSPACE(x) x
#define PASSWSVNAMES(x) x


/* Here's a template record entry:  (PE 2001-09-18)

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "FunctionName" ),
        DESCRIPTION
        (
         "A summary of the function in one sentence.\n"
         "\n"
         "A detailed description of the function. Please, try to be as \n"
         "clear and detailed as possible, this will help both you and \n"
         "others in the long run. \n"
         "\n"
         "Paragraphs are seperated with blank lines.\n"
         "\n"
         "The names of workspace variables and other methods are marked by\n"
         "stars, for example *z_plat*.\n"
         "\n"
         "Generic input and output, and keywords shall be described \n"
         "as exemplified below. If there is no variables of a group, \n"
         "(e.g. generic input) remove that part totally. Note that the \n"
         "on-line help just gives the type of generic input/output and the \n"
         "keyword names, and additional information is for sure needed.\n"
         "\n"
         "Leave space and brake lines when listing input and output \n"
         "variabales to make the code easier to read. See example below. \n"
         "\n"
         "Generic input: \n"
         "   Vector : Vector giving some very important input. Don't \n"
         "            be too short. Use the type of indention used here. \n"
         "\n"
         "Generic output: \n"
         "   Vector : Return vector for the zenith angles. The normal \n"
         "            options are ZA_PENCIL and ZA_SENSOR. \n"
         "\n"
         "Keywords:\n"
         "   delta_t   : Time increment between observations.\n"
         "   z_tan_lim : Vector with start and stop tangent altitudes.\n"
        ),
        AUTHORS( "unknown" ),
        OUTPUT(),
        INPUT( z_plat_, abs_p_, z_abs_, l_step_, refr_, refr_lfac_,
               refr_index_, r_geoid_, z_surface_ ),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "delta_t", "z_tan_lim" ),
        DEFAULTS( NODEF,     NODEF ),
        TYPES(    Numeric_t, Vector_t    )));
  */

  /* Here's an empty record entry:  (PE 2001-09-18)

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "" ),
        DESCRIPTION
        (
         "\n"
         "\n"
         "Generic input: \n"
         "   \n"
         "\n"
         "Generic output: \n"
         "   \n"
         "\n"
         "Keywords:\n"
         "   \n"
        ),
        AUTHORS( "unknown" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));
  */




void define_md_data_raw()
{
  // The variable md_data is defined in file methods_aux.cc.
  extern Array<MdRecord> md_data_raw;


  // Initialize to zero, just in case:
  md_data_raw.resize(0);

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
        DESCRIPTION(
                    "Initialize the WSVs *abs_p*, *abs_t* and *abs_vmrs* from\n"
                    "*p_grid, *t_field* and *vmr_field*.\n"
                    "\n"
                    "This only works for a 1D atmosphere!\n"
                   ) ,
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_p_, abs_t_, abs_vmrs_ ),
        INPUT( atmosphere_dim_, p_grid_, t_field_, vmr_field_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "AbsInputFromRteScalars" ),
        DESCRIPTION(
                    "Initialize absorption input WSVs from local atmospheric conditions.\n"
                    "\n"
                    "The purpose of this method is to allow an explicit line-by-line\n"
                    "calculation, e.g., by *abs_coefCalc*, to be put inside the\n"
                    "*abs_scalar_gas_agenda*. What the method does is to prepare absorption\n"
                    "input parameters (pressure, temperature, VMRs, frequency grid), from\n"
                    "the input parameters to *abs_scalar_gas_agenda*.\n"
                    "There is a matching method to turn the output of *abs_coefCalc*\n"
                    "into what the agenda expects (*abs_scalar_gasFromAbsCoef*).\n"
                    "\n"
                    "Note that the original *f_grid* is distroyed. (This is not a problem\n"
                    "if the method is used inside an agenda.)\n"
                   ) ,
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( f_grid_, abs_p_, abs_t_, abs_vmrs_ ),
        INPUT( f_index_, f_grid_, rte_pressure_, rte_temperature_, rte_vmr_list_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_coefCalc" ),
        DESCRIPTION(
                    "Calculate absorption coefficients. \n"
                    "\n"
                    "This function calculates both, the total absorption (*abs_coef*), and\n"
                    "the absorption per species (*abs_coef_per_species*).\n"
                    "\n"
                    "The method calls four other  methods:\n"
                    "\n"
                    "1. *abs_xsec_per_speciesInit*:\n"
                    "   Initialize *abs_xsec_per_species* \n"
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
                   ) ,
        AUTHORS( "Axel von Engeln", "Stefan Buehler" ),
        OUTPUT( abs_coef_  , abs_coef_per_species_ ),
        INPUT( abs_species_, f_grid_, abs_p_, abs_t_, abs_n2_, abs_h2o_, abs_vmrs_, 
               abs_lines_per_species_, abs_lineshape_,
               abs_cont_names_, abs_cont_models_, 
               abs_cont_parameters_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_coefCalcFromXsec"),
        DESCRIPTION(
                    "Calculate absorption coefficients from cross sections.\n"
                    "\n"
                    "This calculates both the total absorption and the\n"
                    "absorption per tag group. \n"
                    "\n"
                    "Cross sections are multiplied by n*VMR.\n"
                   ),
        AUTHORS( "Stefan Buehler", "Axel von Engeln" ),
        OUTPUT( abs_coef_, abs_coef_per_species_ ),
        INPUT( abs_xsec_per_species_, abs_vmrs_, abs_p_, abs_t_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_coefCalcSaveMemory" ),
        DESCRIPTION(
                    "Calculate absorption coefficients, trying to conserve memory. \n"
                    "\n"
                    "This function calculates only the total absorption (*abs_coef*),\n"
                    "NOT the absorption per tag group (*abs_coef_per_species*).\n"
                    "\n"
                    "This means you cannot use it if you want to calculate Jacobians\n"
                    "later.\n"
                    "\n"
                    "The implementation follows abs_coefCalc.\n"
                   ) ,
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_coef_ ),
        INPUT( abs_species_, f_grid_, abs_p_, abs_t_, abs_n2_, abs_h2o_, abs_vmrs_, 
               abs_lines_per_species_, abs_lineshape_,
               abs_cont_names_, abs_cont_models_, 
               abs_cont_parameters_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_cont_descriptionAppend"),
        DESCRIPTION
        (
         "Appends the description of a continuum model or a complete absorption\n"
         "model to *abs_cont_names* and *abs_cont_parameters*.\n"
         "\n"
         "See online documentation for *abs_cont_names* for a list of\n"
         "allowed models and for information what parameters they require. See\n"
         "file cont.arts in the doc/examples directory for usage examples and\n"
         "default parameters for the various models. \n"
         "\n"
         "Keywords:\n"
         "   name       : The name of a continuum model. Must match one of the models\n"
         "                implemented in ARTS. \n"
         "   option     : give here the option of this continuum/full model.\n"
         "   parameters : A Vector containing the required number of parameters\n"
         "                for the model given. The meaning of the parameters and\n"
         "                how many parameters are required depends on the model.\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_cont_names_, 
                abs_cont_models_,
                abs_cont_parameters_ ),
        INPUT(  abs_cont_names_, 
                abs_cont_models_,
                abs_cont_parameters_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "tagname", "model",  "userparameters" ),
        DEFAULTS( NODEF,     NODEF,    NODEF),
        TYPES(    String_t,  String_t, Vector_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_cont_descriptionInit"),
        DESCRIPTION
        (
         "Initializes the two workspace variables for the continuum description,\n"
         "*abs_cont_names* and *abs_cont_parameters*.\n"
         " \n"
         "This method does not really do anything, except setting the two\n"
         "variables to empty Arrays. It is just necessary because the method\n"
         "*abs_cont_descriptionAppend* wants to append to the variables.\n"
         "   Formally, the continuum description workspace variables are required\n"
         "by the absorption calculation methods (e.g., *abs_coefCalc*). Therefore you\n"
         "always have to call at least *abs_cont_descriptionInit*, even if you do\n"
         "not want to use any continua.\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_cont_names_, 
                abs_cont_models_,
                abs_cont_parameters_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_h2oSet"),
        DESCRIPTION(
                    "Sets abs_h2o to the profile of the first tag group containing\n"
                    "water.\n" 
                    "\n"
                    "This is necessary, because for example *abs_coefCalc* requires abs_h2o\n"
                    "to contain the water vapour profile(the reason for this is the\n"
                    "calculation of oxygen line brodening requires water vapour profile).\n"
                    "Then this function can be used to copy the profile of the first tag\n"
                    "group of water.\n"
                   ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( abs_h2o_ ),
        INPUT( abs_species_, abs_vmrs_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_lineshapeDefine"),
        DESCRIPTION(
                    "Sets the lineshape for all calculated lines.\n"
                    "\n"
                    "   A general lineshape profile is specified, according to a given  \n"
                    "approximation. Alongside a normalization factor is to be set - a \n"
                    "multiplicative forefactor through which the profile can be \n"
                    "modified. This factor is just the 0th or 1st, or 2nd power of the\n"
                    "ratio between the frequency of calculation f and the center frequency\n"
                    "for a specific line f0. A cutoff frequency must also be specified in\n"
                    "order to restrict the calculation within a desired frequency region or\n"
                    "not, when there's no such region.\n"
                    "   The general lineshape profile is given by the keyword shape,\n"
                    "while the normalization factor and the cutoff frequency by\n"
                    "normalizationfactor and cutoff respectively.\n"
                    "   We generate only 1 copy of the lineshape settings. Absorption\n"
                    "routines check for this case and use it for all species.\n"
                    "\n"
                    "   The available values for these keywords are given below.\n"
                    "shape - \"no_shape\" : no specified shape\n"
                    "        \"Doppler\" : Doppler lineshape\n"
                    "        \"Lorentz\" : Lorentz lineshape\n"
                    "        \"Voigt_Kuntz3\" : Kuntz approximation to the Voigt profile,\n"
                    "                         accuracy > 2x10^(-3)\n"
                    "        \"Voigt_Kuntz4\" : Kuntz approximation to the Voigt profile,\n"
                    "                         accuracy > 2x10^(-4)\n"
                    "        \"Voigt_Kuntz6\" : Kuntz approximation to the Voigt profile,\n"
                    "                         accuracy > 2x10^(-6)\n"   
                    "        \"Voigt_Drayson\" : Drayson approximation to the Voigt profile \n"
                    "        \"Rosenkranz_Voigt_Drayson\" : Rosenkrantz oxygen absortion with overlap correction\n" 
                    "                                     on the basis of Drayson routine\n"                                    
                    "        \"Rosenkranz_Voigt_Kuntz6\" : Rosenkrantz oxygen absortion with overlap correction\n"
                    "                                    on the basis of Kuntz routine, accuracy > 2x10^(-6)\n"
                    "        \"CO2_Lorentz\" : Lorentz multiplicated with Cousin's chi factors\n"
                    "        \"CO2_Drayson\" : Drayson multiplicated with Cousin's chi factors\n"
                    "\n"
                    "normalizationfactor - \"no_norm\": 1\n"
                    "                      \"linear\": f/f0\n" 
                    "                      \"quadratic\": (f/f0)^2.\n"
                    "                      \"VVH\": (f*tanh(h*f/(2*k*T))) / (f0*tanh(h*f0/(2*k*T))).\n"
                    "\n"
                    "cutoff - \" -1\" : no cutoff\n"
                    "         \"Number\": positive cutoff frequency in Hz.\n"
                    "\n"
                    "Example usage:\n"
                    "shape=[\"Lorentz\"]\n"
                    "normalizationfactor=[\"linear\"]\n"
                    "cutoff= [650e9]"
                    "\n"
                    "Keywords:\n"
                    "   shape               : The general profile according to an approximation.\n"
                    "   normalizationfactor : The multiplicative forefactor for the general profile.\n"
                    "   cutoff              : The frequency at which a cutoff can be made.\n"),
        AUTHORS( "Axel von Engeln", "Stefan Buehler" ),
        OUTPUT( abs_lineshape_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "shape",  "normalizationfactor", "cutoff" ),
        DEFAULTS( NODEF,    NODEF,                 NODEF ),
        TYPES(    String_t, String_t,              Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_lineshape_per_tgDefine"),
        DESCRIPTION(
                    "Sets the lineshape per tag group for all calculated lines.\n\n"
                    "\n" 
                    "   A general lineshape profile is specified, according to a given  \n"
                    "approximation for each tag group. Alongside a normalization factor\n" 
                    "is to be set also for each tag group - a multiplicative forefactor through\n"
                    "which the profile can be modified. This factor is just the 0th or 1st,\n"
                    "or 2nd power of the ratio between the frequency of calculation f and\n"
                    "the center frequency for a specific line f0. A cutoff frequency must also be\n"
                    "specified for each of the tags in  order to restrict the calculation within\n" 
                    "a desired region or not, when there's no such region.\n"
                    "   The general lineshape profile is given by the keyword shape,\n"
                    "while the normalization factor and the cutoff frequency by\n"
                    "normalizationfactor and cutoff respectively.\n"
                    "\n"
                    "   The available values for these keywords are given below.\n"
                    "shape - \"no_shape\" : no specified shape\n"
                    "        \"Doppler\" : Doppler lineshape\n"
                    "        \"Lorentz\" : Lorentz lineshape\n"
                    "        \"Voigt_Kuntz3\" : Kuntz approximation to the Voigt profile,\n"
                    "                        accuracy > 2x10^(-3)\n"
                    "        \"Voigt_Kuntz4\" : Kuntz approximation to the Voigt profile,\n"
                    "                         accuracy > 2x10^(-4)\n"
                    "        \"Voigt_Kuntz6\" : Kuntz approximation to the Voigt profile,\n"
                    "                         accuracy > 2x10^(-6)\n"   
                    "        \"Voigt_Drayson\" : Drayson approximation to the Voigt profile \n"
                    "        \"Rosenkranz_Voigt_Drayson\" : Rosenkrantz oxygen absortion with overlap correction\n" 
                    "                                     on the basis of Drayson routine\n"                                    
                    "        \"Rosenkranz_Voigt_Kuntz6\" : Rosenkrantz oxygen absortion with overlap correction\n"
                    "                                    on the basis of Kuntz routine, accuracy > 2x10^(-6)\n"
                    "normalizationfactor - \"no_norm\": 1\n"
                    "                      \"linear\": f/f0\n" 
                    "                      \"quadratic\": (f/f0)^2.\n"
                    "cutoff - \" -1\" : no cutoff\n"
                    "           \"Number\": positive cutoff frequency in Hz.\n"
                    "\n"
                    "Example usage:\n"
                    "shape = [\"Lorentz\",\"Voigt_Kuntz6\"] \n"
                    "normalizationfactor= [\"linear\", \"quadratic\"] \n"
                    "cutoff = [ 650e9, -1 ]"
                    "\n"
                    "Keywords:\n"
                    "   shape               : The general profile according to an approximation.\n"
                    "   normalizationfactor : The multiplicative forefactor for the general profile.\n"
                    "   cutoff              : The frequency at which a cutoff can be made.\n"),
        AUTHORS( "Axel von Engeln", "Stefan Buehler" ),
        OUTPUT( abs_lineshape_ ),
        INPUT( abs_species_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "shape",        "normalizationfactor", "cutoff" ),
        DEFAULTS( NODEF,          NODEF,                 NODEF ),
        TYPES(    Array_String_t, Array_String_t,        Vector_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_linesReadFromArts"),
        DESCRIPTION(
                    "Read all the lines from an Arts catalogue file in the \n"
                    "given frequency range. Otherwise a runtime error will be\n"
                    "thrown \n"
                    "\n"
                    "Please note that all lines must correspond\n"
                    "to the legal species / isotope combinations\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name (and path) of the catalogue file.\n"
                    "   fmin     : Minimum frequency for lines to read in Hz.\n"
                    "   fmax     : Maximum frequency for lines to read in Hz.\n"),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_lines_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "filename", "fmin",    "fmax" ),
        DEFAULTS( NODEF,      NODEF,     NODEF ),
        TYPES(    String_t,   Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_linesReadFromArtsObsolete"),
        DESCRIPTION(
                    "Read all the lines from an Arts catalogue file in the \n"
                    "given frequency range. Otherwise a runtime error will be\n"
                    "thrown \n"
                    "\n"
                    "Please note that all lines must correspond\n"
                    "to the legal species / isotope combinations\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name (and path) of the catalogue file.\n"
                    "   fmin     : Minimum frequency for lines to read in Hz.\n"
                    "   fmax     : Maximum frequency for lines to read in Hz.\n"),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_lines_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "filename", "fmin",    "fmax" ),
        DEFAULTS( NODEF,      NODEF,     NODEF ),
        TYPES(    String_t,   Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_linesReadFromHitran"),
        DESCRIPTION
        (
         "Read all the lines from a HITRAN 1986-2001 catalogue file in\n"
         "the given frequency range. Otherwise a runtime error will be\n"
         "thrown. For HITRAN 2004 line data use the workspace method \n"
         "abs_linesReadFromHitran. \n"
         "\n"
         "Please note that all lines must correspond to the legal\n"
         "species / isotope combinations and that the line data \n"
         "file must be sorted by increasing frequency\n"
         "\n"
         "WWW access of the HITRAN catalog: http://www.hitran.com/\n"
         "\n"
         "Keywords: \n"
         "   filename : Name (and path) of the catalogue file.\n"
         "   fmin     : Minimum frequency for lines to read in Hz.\n"
         "   fmax     : Maximum frequency for lines to read in Hz.\n"),
        AUTHORS( "Thomas Kuhn" ),
        OUTPUT( abs_lines_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "filename",  "fmin",    "fmax" ),
        DEFAULTS( NODEF,       NODEF,     NODEF ),
        TYPES(    String_t,    Numeric_t, Numeric_t)));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_linesReadFromHitran2004"),
        DESCRIPTION(
                    "Read all the lines from a HITRAN 2004 catalogue file in the \n"
                    "given frequency range. Otherwise a runtime error is thrown. \n"
                    "\n"
                    "Records of molecules unknown to ARTS are ignored but a \n"
                    "warning is issued. In particular this happens to CH3OH \n"
                    "(HITRAN molecule number 39) because there is no total internal \n"
                    "partition sum available. \n"
                    "\n"
                    "The database must be sorted by increasing frequency!\n"
                    "\n"
                    "WWW access of the HITRAN catalog: http://www.hitran.com/\n"
                    "\n"
                    "For data in the Hitran 1986-2001 format use the workspace \n"
                    "method: abs_linesReadFromHitran\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name (and path) of the catalogue file.\n"
                    "   fmin     : Minimum frequency for lines to read in Hz.\n"
                    "   fmax     : Maximum frequency for lines to read in Hz.\n"),
        AUTHORS( "Hermann Berg", "Thomas Kuhn" ),
        OUTPUT( abs_lines_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "filename",  "fmin",    "fmax" ),
        DEFAULTS( NODEF,       NODEF,     NODEF ),
        TYPES( String_t, Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_linesReadFromJpl"),
        DESCRIPTION(
                    "Read all the lines from a JPL catalogue file in the \n"
                    "given frequency range. Otherwise a runtime error will be\n"
                    "thrown\n"
                    "\n"
                    "Please note that all lines must correspond\n"
                    "to the legal species / isotope combinations.\n"
                    "\n"
                    "WWW access of the JPL catalog: http://spec.jpl.nasa.gov/\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name (and path) of the catalogue file.\n"
                    "   fmin     : Minimum frequency for lines to read in Hz.\n"
                    "   fmax     : Maximum frequency for lines to read in Hz.\n"),
        AUTHORS( "Thomas Kuhn" ),
        OUTPUT( abs_lines_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "filename",  "fmin", "fmax" ),
        DEFAULTS( NODEF,       NODEF,     NODEF ),
        TYPES( String_t, Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_linesReadFromMytran2"),
        DESCRIPTION(
                    "Read all the lines from a MYTRAN2 catalogue file in the \n"
                    "given frequency range. Otherwise a runtime error will be\n"
                    "thrown\n"
                    "\n"
                    "Please note that all lines must correspond\n"
                    "to the legal species / isotope combinations\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name (and path) of the catalogue file.\n"
                    "   fmin     : Minimum frequency for lines to read in Hz.\n"
                    "   fmax     : Maximum frequency for lines to read in Hz.\n"),
        AUTHORS( "Axel von Engeln", "Stefan Buehler" ),
        OUTPUT( abs_lines_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "filename", "fmin", "fmax"),
        DEFAULTS( NODEF,       NODEF,     NODEF ),
        TYPES( String_t, Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_lines_per_speciesAddMirrorLines"),
        DESCRIPTION(
                    "Adds mirror lines at negative frequencies to the *abs_lines_per_species*.\n"
                    "\n"
                    "For each line at frequency +f in *abs_lines_per_species* a corresponding\n"
                    "entry at frequency -f is added to *abs_lines_per_species*.The mirror \n"
                    "lines are appended to the line lists after the original lines.\n"),
        AUTHORS( "Axel von Engeln", "Stefan Buehler" ),
        OUTPUT( abs_lines_per_species_ ),
        INPUT( abs_lines_per_species_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_lines_per_speciesCompact"),
        DESCRIPTION(
                    "Removes all lines outside the defined lineshape cutoff frequency\n"
                    "from the *abs_lines_per_species*. This can save computation time.\n"
                    "It should be particularly useful to call this method after\n"
                    "*abs_lines_per_speciesAddMirrorLines*.\n"),
        AUTHORS( "Axel von Engeln", "Stefan Buehler" ),
        OUTPUT( abs_lines_per_species_ ),
        INPUT( abs_lines_per_species_, abs_lineshape_, f_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_lines_per_speciesCreateFromLines"),
        DESCRIPTION(
                    "Split lines up into the different tag groups.\n"
                    "\n"
                    "The tag groups are tested in the order in which they are\n" 
                    "specified in the controlfile. The lines are assigned to \n"
                    "the tag groups in the order as the groups  are specified.\n"
                    "That means if you do [\"O3-666\",\"O3\"],the last group O3 \n"
                    "gets assigned all the O3 lines that do not fit in the first group.\n"),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_lines_per_species_ ),
        INPUT( abs_lines_, abs_species_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_lines_per_speciesReadFromCatalogues"),
        DESCRIPTION(
                    "This method can read lines from different line \n"
                    "catalogues.\n"
                    "\n"
                    "For each tag group, you can specify which catalogue\n"
                    "to use. Because the method creates abs_lines_per_species directly,\n"
                    "it replaces for example thefollowing two method calls:\n"
                    "  - abs_linesReadFromHitran\n"
                    "  - abs_lines_per_speciesCreateFromLines\n"
                    "   This method needs as input WSVs the list of tag \n"
                    "groups. Keyword parameters must specify the names of\n"
                    "the catalogue files to use and the matching formats.\n"
                    "Names can be anything, formats can currently be \n"
                    "HITRAN96 (for HITRAN 1986-2001 databases), HITRAN04 \n"
                    "(for HITRAN 2004 database), MYTRAN2, JPL, or ARTS. \n"
                    "Furthermore, keyword parameters have to specify minimum \n"
                    "and maximum frequency for each tag group. To safe typing, \n"
                    "if there are less elements in the keyword parameters than \n"
                    "there are tag groups, the last parameters are applied to \n"
                    "all following tag groups.\n"
                    "\n"
                    "Example usage:\n"
                    "\n"
                    "abs_lines_per_speciesReadFromCatalogues{\n"
                    "  filenames = [ \"../data/cat1.dat\", \"../data/cat2.dat\" ]\n"
                    "  formats   = [ \"MYTRAN2\",          \"HITRAN96\"         ]\n"
                    "  fmin      = [ 0,                  0                  ]\n"
                    "  fmax      = [ 2000e9,             100e9              ]\n"
                    "}\n"
                    "   In this example, lines for the first tag group will\n"
                    "be taken from cat1, lines for all other tag groups \n"
                    "will be taken from cat2.\n"
                    "   This methods allows you for example to use a \n"
                    "special line file just for water vapor lines. This\n"
                    "could be the  improved water vapor line file \n"
                    "generated by Thomas Kuhn.\n"
                    "   Catalogues are only read once, even if several tag\n"
                    "groups have the same catalogue. However, in that case\n"
                    "the frequency ranges MUST be the same. (If you want \n"
                    "to do fine-tuning of the frequency ranges, you can do \n"
                    "this inside the tag definitions, e.g., \"H2O-*-0-2000e9\".)\n"
                    "   This function uses the various reading routines\n"
                    "(abs_linesReadFromHitran, etc.), as well as\n"
                    "abs_lines_per_speciesCreateFromLines.\n"
                    "\n"
                    "Keywords: \n"
                    "   filenames : Name (and path) of the catalogue files.\n"
                    "   formats   : allowed formats are HITRAN96,MYTRAN2,JPL,ARTS \n"
                    "   fmin      : Minimum frequency for lines to read in Hz.\n"
                    "   fmax      : Maximum frequency for lines to read in Hz.\n"),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_lines_per_species_ ),
        INPUT( abs_species_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "filenames",    "formats",      "fmin",   "fmax" ),
        DEFAULTS( NODEF,          NODEF,          NODEF,    NODEF ),
        TYPES(    Array_String_t, Array_String_t, Vector_t, Vector_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_lines_per_speciesSetEmpty"),
        DESCRIPTION
        (
         "Sets abs_lines_per_species to empty line lists.\n"
         "\n"
         "You can use this method to set lines per tag if you do not reall want\n"
         "to compute line spectra. Formally, abs_coefCalc will still require\n"
         "abs_lines_per_species to be set.\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_lines_per_species_ ),
        INPUT( abs_species_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  // New name: abs_lookupAdapt
  md_data_raw.push_back     
    ( MdRecord
      ( NAME("abs_lookupAdapt"),
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
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_lookup_, abs_lookup_is_adapted_ ),
        INPUT(  abs_lookup_, abs_species_, f_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("abs_lookupCreate"),
        DESCRIPTION
        (
         "Creates a gas absorption lookup table.\n"
         "\n"
         "The lookup table stores absorption cross-sections as a function of\n"
         "pressure. Additionally, absorption can be stored as a function of\n"
         "temperature for temperature perturbations from a reference\n"
         "profile. \n"
         "\n"
         "Additionally, absorption can be stored as a function of water vapor\n"
         "VMR perturbations from a reference profile. The variable *abs_nls*\n"
         "specifies, for which species water vapor perturbations should be\n"
         "generated. \n"
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
        OUTPUT( abs_lookup_, abs_lookup_is_adapted_ ),
        INPUT( abs_species_, 
               abs_lines_per_species_,
               abs_lineshape_,
               abs_nls_,
               f_grid_,
               abs_p_,
               abs_vmrs_,
               abs_t_, 
               abs_t_pert_, 
               abs_nls_pert_,
               abs_n2_,
               abs_cont_names_,
               abs_cont_models_, 
               abs_cont_parameters_
               ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  // New name: abs_lookupInit
  md_data_raw.push_back     
    ( MdRecord
      ( NAME("abs_lookupInit"),
        DESCRIPTION
        (
         "Creates an empty gas absorption lookup table.\n"
         "\n"
         "This is mainly there to help developers. For example, you can write\n"
         "the empty table to an XML file, to see the file format.\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_lookup_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ))) ;

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("abs_lookupSetup"),
        DESCRIPTION
        (
         "Set up input parameters for abs_lookupCreate.\n"
         "\n"
         "More information can be found in the documentation for method\n"
         "*abs_lookupSetupBatch*\n"
         "\n"
         "Keywords:\n"
         "   p_step   : Maximum step in log(p[Pa]) (natural logarithm, as always). If\n"
         "              the pressure grid is coarser than this, additional points\n"
         "              are added until each log step is smaller than this.\n"
         "              Has a default value.\n"
         "   t_step   : The temperature variation grid step in Kelvin, for a 2D\n"
         "              or 3D atmosphere. For a 1D atmosphere this parameter is\n"
         "              not used. Has a default value.\n"
         "   h2o_step : The H2O variation grid step [fractional], if H2O variations are done\n"
         "              (which is determined automatically, based on abs_species\n"
         "              and the atmospheric dimension). For a 1D atmosphere this parameter is\n"
         "              not used. Has a default value.\n"
         "\n"
         "See also: \n"
         "   *abs_lookupSetupBatch*\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT(  abs_p_,
                 abs_t_, 
                 abs_t_pert_, 
                 abs_vmrs_,
                 abs_nls_,
                 abs_nls_pert_ ),
        INPUT(   atmosphere_dim_,
                 p_grid_,
                 lat_grid_,
                 lon_grid_,
                 t_field_,
                 vmr_field_,
                 abs_species_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "p_step",  "t_step",  "h2o_step" ),
        DEFAULTS( "0.05",    "5",       "0.5" ),
        TYPES(    Numeric_t, Numeric_t, Numeric_t )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("abs_lookupSetupBatch"),
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
         "*h2o_abs*, and hence require nonlinear treatment in the lookup table.\n"
         "\n"
         "The method also checks which range of pressures, temperatures, and\n"
         "VMRs occurs, and sets *abs_p*, *abs_t*, *abs_t_pert*, and *abs_vmrs*\n"
         "accordingly.\n"
         "\n"
         "If nonlinear species are present, *abs_nls* and *abs_nls_pert* are also\n"
         "generated. \n"
         "\n"
         "Keywords:\n"
         "   p_step   : Maximum step in log(p[Pa]) (natural logarithm, as always). If\n"
         "              the pressure grid is coarser than this, additional points\n"
         "              are added until each log step is smaller than this.\n"
         "              Has a default value.\n"
         "   t_step   : The temperature variation grid step in Kelvin, for a 2D\n"
         "              or 3D atmosphere. For a 1D atmosphere this parameter is\n"
         "              not used. Has a default value.\n"
         "   h2o_step : The H2O variation grid step [fractional], if H2O variations are done\n"
         "              (which is determined automatically, based on abs_species\n"
         "                and the atmospheric dimension). For a 1D atmosphere this parameter is\n"
         "              not used. Has a default value.\n"
         "   extremes : You can give here explicit extreme values to add to\n"
         "              abs_t_pert and abs_nls_pert. The order is [t_pert_min,\n"
         "              t_pert_max, nls_pert_min, nls_pert_max]. Has a default value of empty.\n"
         "See also: \n"
         "   *abs_lookupSetup*\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT(  abs_p_,
                 abs_t_, 
                 abs_t_pert_, 
                 abs_vmrs_,
                 abs_nls_,
                 abs_nls_pert_ ),
        INPUT(   abs_species_,
                 batch_atm_fields_compact_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "p_step",  "t_step",  "h2o_step", "extremes" ),
        DEFAULTS( "0.05",    "5",       "0.5",      "[]" ),
        TYPES(    Numeric_t, Numeric_t, Numeric_t,  Vector_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_n2Set"),
        DESCRIPTION(
                    "Sets abs_n2 to the profile of the first tag group containing\n"
                    "molecular nitrogen. See *abs_h2oSet* for more details.\n"
                   ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( abs_n2_ ),
        INPUT( abs_species_, abs_vmrs_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("abs_scalar_gasExtractFromLookup"),
        DESCRIPTION
        (
         "Extract scalar gas absorption coefficients from lookup table.\n"
         "\n"
         "This extracts the absorption coefficient for all species in the\n"
         "current calculation from the lookup table. Extraction is for one\n"
         "specific atmospheric condition, i.e., a set of pressure, temperature,\n"
         "and VMR values.\n"
         "\n"
         "Extraction can be either for a single frequency (f_index>=0), or for\n"
         "all frequencies (f_index<0). The dimension of the output\n"
         "abs_scalar_gas is adjusted accordingly.\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_scalar_gas_ ),
        INPUT(  abs_lookup_, abs_lookup_is_adapted_,
                f_index_, 
                rte_pressure_, rte_temperature_, rte_vmr_list_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("abs_scalar_gasFromAbsCoef"),
        DESCRIPTION
        (
         "Copy *abs_scalar_gas* from *abs_coef*. This is handy for putting an\n"
         "explicit line-by-line calculation into the\n"
         "*abs_scalar_gas_agenda*. See also method *AbsInputFromRteScalars*.\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_scalar_gas_ ),
        INPUT(  abs_coef_per_species_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("abs_fieldCalc"),
        DESCRIPTION
        (
         "Calculate scalar gas absorption for all points in the atmosphere.\n"
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
         "*abs_scalar_gas_agenda*.\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_field_ ),
        INPUT(  abs_scalar_gas_agenda_,
                f_index_,
                f_grid_,
                atmosphere_dim_,
                p_grid_, lat_grid_, lon_grid_,
                t_field_, vmr_field_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_speciesAdd"),
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
         "\n"
         "Keywords:\n"
         "   species : Specify one String for each tag group that you want to\n"
         "             add. Inside the String, separate the tags by commas\n"
         "             (plus optional blanks).\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_species_ ),
        INPUT(  abs_species_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "species" ),
        DEFAULTS( NODEF ),
        TYPES(    Array_String_t   )));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_speciesAdd2"),
        DESCRIPTION
        (
         "Adds a species tag group to the list of absorption species and \n"
         "jacobian quantities.\n"
         "\n"
         "The method is basically a combined call of *abs_speciesAdd* and\n"
         "*jacobianAddAbsSpecies*. In this way it is not needed to specify a\n"
         "tag group in two different places. \n"
         "\n"
         "Arguments exactly as for *jacobianAddAbsSpecies*. Note that this\n"
         "method only handles a single tag group, in contrast to \n"
         "*abs_speciesAdd*\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( abs_species_, jacobian_quantities_, jacobian_agenda_ ),
        INPUT( abs_species_, jacobian_, atmosphere_dim_, p_grid_, lat_grid_, 
               lon_grid_ ),
        GOUTPUT(),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS( "species", "method", "unit", "dx" ),
        DEFAULTS( NODEF,     NODEF,    NODEF,  NODEF ),
        TYPES( String_t, String_t, String_t, Numeric_t )));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_speciesDefineAllInScenario"),
        DESCRIPTION
        (
         "Define one tag group for each species known to ARTS and included in an\n"
         "atmospheric scenario.\n"
         "\n"
         "You can use this as an alternative to tgsDefine if you want to make an\n"
         "absorption calculation that is as complete as possible. The method\n"
         "goes through all defined species and tries to open the VMR file. If\n"
         "this works the tag is included, otherwise it is skipped.\n"
         "\n"
         "Keywords:\n"
         "   basename : The name and path of a particular atmospheric scenario.\n"
         "              For example: /pool/lookup2/arts-data/atmosphere/fascod/tropical\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_species_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "basename" ),
        DEFAULTS( NODEF ),
        TYPES( String_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_speciesInit"),
        DESCRIPTION
        (
         "Sets  *abs_speciesSet* to be empty.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_species_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  // This is duplicate with the 1-0 method tgsDefine. Merge!
  md_data_raw.push_back
    ( MdRecord
      ( NAME("SpeciesSet"),
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
         "names, e.g., \"O3\".  Isotopes are given by the last digit of the atomic\n"
         "weight, i.g., \"O3-668\" for the asymmetric ozone molecule including an\n"
         "oxygen 18 atom. Groups of transitions are specified by giving a lower\n"
         "and upper limit of a frequency range, e.g., \"O3-666-500e9-501e9\".\n"
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
         "Generic Output:\n"
         "   ArrayOfArrayOfSpeciesTag : Output tag groups.\n"
         "\n"
         "Keywords:\n"
         "   species : Specify one String for each tag group that you want to\n"
         "             create. Inside the String, separate the tags by commas\n"
         "             (plus optional blanks).\n"
         "\n"
         "Example:\n"
         "\n"
         "   species = [ \"O3-666-500e9-501e9, O3-686\",\n"
         "               \"O3\",\n"
         "               \"H2O-PWR98\" ]\n"
         "\n"
         "   The first tag group selects all O3-666 lines between 500 and\n"
         "   501 GHz plus all O3-686 lines.  \n"
         "\n"
         "   The second tag group selects all remaining O3 transitions.\n"
         "\n"
         "   The third tag group selects H2O, with one of the complete\n"
         "   absorption models (Rosenkranz 98). No spectrocopic line catalogue\n"
         "   data will be used for that third tag group.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT(  ArrayOfArrayOfSpeciesTag_ ),
        GINPUT( ),
        KEYWORDS( "species" ),
        DEFAULTS( NODEF ),
        TYPES(    Array_String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_vecAddGas"),
        DESCRIPTION
        (
         "Add gas absorption to first element of absorption vector.\n"
         "\n"
         "The task of this method is to sum up the gas absorption of the\n"
         "different gas species and add the result to the first element of the\n"
         "absorption vector.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT(abs_vec_),
        INPUT(abs_vec_, abs_scalar_gas_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_vecAddPart"),
        DESCRIPTION
        (
         "The particle absorption is added to *abs_vec* \n"
         "\n"
         "This function sums up the absorption vectors for all particle \n"
         "types weighted with particle number density.\n"
         "The resluling absorption vector is added to the workspace \n"
         "variable *abs_vec* \n"
         "Output and input of this method is *abs_vec* (stokes_dim).\n"
         "The inputs are the absorption vector for the single particle type \n"
         "*abs_vec_spt* (part_types, stokes_dim) and the local particle\n"
         " number densities for all particle types namely the \n"
         "*pnd_field* (part_types, p_grid, lat_grid, lon_grid, ) for given \n"
         "*p_grid*, *lat_grid*, and *lon_grid*. The particle types required \n"
         "are specified in the control file.\n"
         ),
        AUTHORS( "Sreerekha T.R." ),
        OUTPUT(abs_vec_),
        INPUT(abs_vec_, abs_vec_spt_, pnd_field_, atmosphere_dim_,
              scat_p_index_,  scat_lat_index_, scat_lon_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_vecInit"),
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
        OUTPUT(abs_vec_),
        INPUT(f_grid_, stokes_dim_, f_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_xsec_per_speciesAddConts"),
        DESCRIPTION(
                    "Calculate absorption cross sections per tag group for continua.\n"
                   ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_xsec_per_species_ ),
        INPUT( abs_species_, f_grid_, abs_p_, abs_t_, abs_n2_, abs_h2o_, abs_vmrs_,
               abs_cont_names_, abs_cont_parameters_,
               abs_cont_models_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_xsec_per_speciesAddLines"),
        DESCRIPTION(
                    "Calculates the line spectrum for each tag group and adds\n"
                    "it to abs_xsec_per_species.\n"
                   ),
        AUTHORS( "Stefan Buehler", "Axel von Engeln" ),
        OUTPUT( abs_xsec_per_species_ ),
        INPUT( abs_species_, f_grid_, abs_p_, abs_t_, abs_h2o_, abs_vmrs_, 
               abs_lines_per_species_, abs_lineshape_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "abs_xsec_per_speciesInit" ),
        DESCRIPTION(
                    "Initialize *abs_xsec_per_species*.\n"
                    "\n"
                    "The initialization is\n"
                    "necessary, because methods *abs_xsec_per_speciesAddLines*\n"
                    "and *abs_xsec_per_speciesAddConts* just add to *abs_xsec_per_species*.\n"
                    "The size is determined from *tgs*.\n"
                   ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( abs_xsec_per_species_ ),
        INPUT( abs_species_, f_grid_, abs_p_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AgendaExecute"),
        DESCRIPTION
        ( 
         "Execute an agenda.\n"
         "\n"
         "Generic input:\n"
         "   Agenda : The agenda.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Agenda_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES(),
        AGENDAMETHOD( false )));
      
  md_data_raw.push_back
    ( MdRecord
      ( NAME("AgendaSet"),
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
         "check, whether the given methods use the right input WSVs and produce\n"
         "the right output WSVs.\n"
         " \n"
         "Generic output:\n"
         "   Agenda : The new agenda.\n"
         "\n"
         "Keywords:\n"
         "   No keywords, but other methods can appear in the method body.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(  ),
        INPUT(  ),
        GOUTPUT( Agenda_ ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD(   true  ),
        SUPPRESSHEADER( false ),
        PASSWORKSPACE(  false ),
        PASSWSVNAMES(   true  )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AntennaOff"),
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
        OUTPUT( antenna_dim_, mblock_za_grid_, mblock_aa_grid_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AntennaSet1D"),
        DESCRIPTION
        (
         "Sets the antenna dimension to 1D.\n"
         "\n"
         "Sets *antenna_dim* to 1 and sets *mblock_aa_grid* to be empty.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( antenna_dim_, mblock_aa_grid_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AntennaSet2D"),
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
        OUTPUT( antenna_dim_ ),
        INPUT( atmosphere_dim_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("antenna_diagramAppendArray"),
        DESCRIPTION
        (
         "Appends a ArrayOfMatrix to *antenna_diagram*.\n"
         "\n"
         "This method can be used both to initialise and expand\n"
         "the viewing angles of *antenna_diagram*. At least one viewing\n"
         "angle must be given in *antenna_diagram*, and the array that\n"
         "is appended must have at least one element but not more than\n"
         "the number of polarisation given by *sensor_pol*.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( antenna_diagram_ ),
        INPUT( sensor_pol_ ),
        GOUTPUT( ),
        GINPUT( ArrayOfMatrix_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Append"),
        DESCRIPTION
        (
         "Append a workspace variable to another workspace variable.\n"
         "\n"
         "This is a supergeneric method. It can append a workspace variable\n"
         "to another workspace variable of the same group. (E.g., a Matrix to\n"
         "another Matrix.)\n"
         "\n"
         "This method is not implemented for all types, just for those where an\n"
         "append makes sense. A runtime error is thrown if one attempts to use\n"
         "it on types that are not implemented.\n"         
         "\n"
         "As allways, output comes first in the argument list!\n"
         "\n"
         "Usage example:\n"
         "\n"
         "Append(array_of_matrix_1,array_of_matrix_2){}\n"
         "\n"
         "Will append the matrix array 2 to matrix array 1.\n"
         "\n"
         "Supergeneric output:\n"
         "   Any_ : The output variable.\n"
         "\n"
         "Supergeneric input:\n"
         "   Any_ : The input variable.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Any_ ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfGriddedField3ExtractFromArrayOfArrayOfGriddedField3"),
        DESCRIPTION
        (
         "Extract an ArrayOfGriddedField3 from an array of arrays of GriddedField3.\n"
         "\n"
         "Copies *ArrayOfGriddedField3* with given Index from input\n"
         "*ArrayOfArrayOfGriddedField3* variable to create output\n"
         "*ArrayOfGriddedField3*.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ArrayOfGriddedField3_ ),
        GINPUT(  ArrayOfArrayOfGriddedField3_, Index_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexExtractFromArrayOfArrayOfIndex"),
        DESCRIPTION
        (
         "Extract an ArrayOfIndex from an array of arrays of Index.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ArrayOfIndex_ ),
        GINPUT(  ArrayOfArrayOfIndex_, Index_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ArrayOfMatrixCreate" ),
        DESCRIPTION
        (
         "Creates an empty ArrayOfMatrix.\n"
         "\n"
         "If the variable already exists, it will be reset.\n"
         "\n"
         "Generic output: \n"
         "   ArrayOfMatrix: New empty ArrayOfMatrix.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfMatrix_ ),
        GINPUT(),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixInsert"),
        DESCRIPTION
        (
         "Inserts a Matrix in an ArrayOfMatrix.\n"
         "\n"
         "The keyword can be used to chose which element will be set, If a\n"
         "negative number is given, the matrix will be appended to the array.\n"
         "Note that zero-based indexing is used.\n"
         "\n"
         "Generic output:\n"
         "  ArrayOfMatrix : The new array.\n"
         "\n"
         "Generic input:\n"
         "  ArrayOfMatrix : The original array.\n"
         "         Matrix : The matrix to insert.\n"
         "Keywords:\n"
         "        element : The index to be set.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ArrayOfMatrix_ ),
        GINPUT( ArrayOfMatrix_, Matrix_ ),
        KEYWORDS( "element" ),
        DEFAULTS( NODEF ),
        TYPES( Index_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ArrayOfStringCreate" ),
        DESCRIPTION
        (
         "Creates an empty ArrayOfString.\n"
         "\n"
         "If the variable already exists, it'll be reset.\n"
         "\n"
         "Generic output: \n"
         "   ArrayOfString: New empty ArrayOfString.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfString_ ),
        GINPUT(),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfStringSet"),
        DESCRIPTION
        (
         "Sets a String array according the given text.\n"
         "The format is text = [\"String1\",\"String2\",...]\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfString_ ),
        GINPUT(),
        KEYWORDS( "text" ),
        DEFAULTS( NODEF ),
        TYPES(    Array_String_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Arts"),
        DESCRIPTION
        ( 
         "Run the agenda that is specified inside the curly braces. ARTS\n"
         "controlfiles must define this method. It is executed automatically\n"
         "when ARTS is run on the controlfile.\n" 
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD( true )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmFieldsCalc"),
        DESCRIPTION
        (
         "Interpolate the atmospheric fields.\n"
         "\n"
         "An atmospheric scenario includes the following data for each \n"
         "position (pressure, latitude, longitude) in the atmosphere: \n"
         "           1. temperature field \n"
         "           2. the corresponding altitude field \n"
         "           3. vmr fields for the gaseous species \n"
         "This method interpolates the fields from the raw data\n"
         "(*t_field_raw*, *z_field_raw*) which can be stored on \n"
         "arbitrary grids on the grids for the calculation\n"
         "(*p_grid*, *lat_grid*, *lon_grid*).\n"
        ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT(t_field_, z_field_, vmr_field_),
        INPUT(p_grid_, lat_grid_, lon_grid_, t_field_raw_, z_field_raw_, 
              vmr_field_raw_, atmosphere_dim_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmFieldsCalcExpand1D"),
        DESCRIPTION
        (
         "Interpolate 1D raw atmospheric fields to create 2D or 3D \n"
         "homogeneous atmospheric fields.\n"
         "\n"
         "The method works as *AtmFieldsCalc* but accepts only raw 1D\n"
         "atmospheres. The raw atmosphere is interpolated to *p_grid* and \n"
         "the obtained values are applied for all latitudes, and also \n"
         "longitudes for 3D, to create a homogeneous atmosphere. \n"
         "\n"
         "The method deals only with the atmospheric fields, and to create\n"
         "a true 2D or 3D version of a 1D case, a demand is also that the\n"
         "geoid radius is set to be constant for all latitudes/longitudes.\n"
        ),
        AUTHORS( "Patrick Eriksson", "Claudia Emde" ),
        OUTPUT( t_field_, z_field_, vmr_field_ ),
        INPUT( p_grid_, lat_grid_, lon_grid_, t_field_raw_, z_field_raw_, 
               vmr_field_raw_, atmosphere_dim_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmFieldsRefinePgrid"),
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
         "\n"
         "Keywords:\n"
         "   p_step   : Maximum step in log(p[Pa]) (natural logarithm, as always). If\n"
         "              the pressure grid is coarser than this, additional points\n"
         "              are added until each log step is smaller than this.\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT(p_grid_,
               t_field_, z_field_, vmr_field_),
        INPUT( p_grid_, lat_grid_, lon_grid_,
               t_field_, z_field_, vmr_field_, atmosphere_dim_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "p_step" ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("atm_fields_compactAddConstant"),
        DESCRIPTION
        (
         "Adds a constant field to atm_fields_compact. \n"
         "\n"
         "This is handy for nitrogen or oxygen. The constant value is\n"
         "appended at the end of the fields that are already there.  All\n"
         "dimensions (pressure, latitude, longitude) are filled up, so this\n"
         "works for 1D, 2D, or 3D atmospheres.\n"
         "\n"
         "Keywords:\n"
         "   name  : The field name. Use, e.g., vmr_o2 for oxygen VMR.\n"
         "   value : The constant value of this field.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( atm_fields_compact_ ),
        INPUT(  atm_fields_compact_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "name",   "value" ),
        DEFAULTS( NODEF,    NODEF ),
        TYPES(    String_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("atm_fields_compactFromMatrix"),
        DESCRIPTION
        (
         "Set atm_fields_compact from 1D profiles in a matrix.\n"
         "\n"
         "For clear-sky batch calculations it is handy to store atmospheric\n"
         "profiles in an array of matrix. We take such a matrix, and create\n"
         "*atm_fields_compact* from it. \n"
         "\n"
         "The matrix must contain one row for each pressure level. Recommended\n"
         "row format:\n"
         "\n"
         "p[Pa] T[K] z[m] VMR_1[1] ... VMR[2]\n"
         "\n"
         "Works only for *atmosphere_dim==1.*\n"         
         "\n"
         "Keywords:\n"
         "   field_names : Field names to store in atm_fields_compact.\n"
         "                 This should be, e.g.:\n"
         "                 [\"T[K]\", \"z[m]\", \"vmr_h2o[1]\"]\n"
         "                 There must be one name less than matrix columns,\n"
         "                 because the first column must contain pressure.\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( atm_fields_compact_ ),
        INPUT(  atmosphere_dim_ ),
        GOUTPUT( ),
        GINPUT( Matrix_ ),
        KEYWORDS( "field_names" ),
        DEFAULTS( NODEF ),
        TYPES(    Array_String_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("batch_atm_fields_compactFromArrayOfMatrix"),
        DESCRIPTION
        (
         "Expand batch of 1D atmospheric states to a batch_atm_fields_compact.\n"
         "\n"
         "This is used to handle 1D batch cases, for example from the Chevallier\n"
         "data set, stored in a matrix. \n"
         "\n"
         "The matrix must contain one row for each pressure level. Row format:\n"
         "\n"
         "p[Pa] T[K] z[m] VMR_1[1] ... VMR_N[1]\n"
         "\n"
         "Keywords:\n"
         "   field_names : Field names to store in atm_fields_compact.\n"
         "                 This should be, e.g.:\n"
         "                 [\"T\", \"z\", \"H2O\", \"O3\"]\n"
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
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( batch_atm_fields_compact_ ),
        INPUT(  atmosphere_dim_ ),
        GOUTPUT( ),
        GINPUT( ArrayOfMatrix_ ),
        KEYWORDS( "field_names", "extra_field_names", "extra_field_values" ),
        DEFAULTS( NODEF,         "[]",                "[]" ),
        //        DEFAULTS( NODEF,         NODEF,                NODEF ),
        TYPES(    Array_String_t, Array_String_t,     Vector_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmFieldsFromCompact"),
        DESCRIPTION
        (
         "Extract pressure grid and atmospheric fields from\n"
         "*atm_fields_compact*.\n"
         "\n"
         "An atmospheric scenario includes the following data for each \n"
         "position (pressure, latitude, longitude) in the atmosphere: \n"
         "           1. temperature field \n"
         "           2. the corresponding altitude field \n"
         "           3. vmr fields for the gaseous species \n"
         "\n"
         "This method just splits up the data found in *atm_fields_compact* to\n"
         "p_grid, lat_grid, lon_grid, and the various fields. No interpolation.\n"
         "See documentation of *atm_fields_compact* for a definition of the data.\n"
         "\n"
         "There are some safety checks on the names of the fields: The first\n"
         "field must be called *T*, the second *z*.  Remaining fields must be\n"
         "trace gas species volume mixing ratios, named for example \"H2O\", \"O3\",\n"
         "and so on. The species names must fit the species in *abs_species*.\n"
         "(Same species in same order.) Only the species name must fit, not the\n"
         "full tag.\n"
         "\n"
         "Possible future extensions: Add a keyword parameter to refine the\n"
         "pressure grid if it is too coarse. Or a version that interpolates onto\n"
         "given grids, instead of using and returning the original grids.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( p_grid_, lat_grid_, lon_grid_, t_field_, z_field_, vmr_field_ ),
        INPUT(  abs_species_, atm_fields_compact_, atmosphere_dim_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmosphereSet1D"),
        DESCRIPTION
        (
         "Sets the atmospheric dimension to 1D.\n"
         "\n"
         "Sets *atmosphere_dim* to 1 and gives some variables dummy values.\n"
         "\n"
         "The latitude and longitude grids are set to be empty.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( atmosphere_dim_, lat_grid_, lon_grid_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmosphereSet2D"),
        DESCRIPTION
        (
         "Sets the atmospheric dimension to be 2D.\n"
         "\n"
         "Sets *atmosphere_dim* to 2 and gives some variables dummy values.\n"
         "\n"
         "The longitude grid is set to be empty. The variables *lat_1d*\n"
         "and *meridian_angle_1d* are given values that cause an error\n"
         "message if used.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( atmosphere_dim_, lon_grid_, lat_1d_, meridian_angle_1d_),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmosphereSet3D"),
        DESCRIPTION
        (
         "Sets the atmospheric dimension to 3D.\n"
         "\n"
         "Sets *atmosphere_dim* to 3 and gives some variables dummy values.\n"
         "\n"
         "The variables *lat_1d* and *meridian_angle_1d* are given\n"
         "values that cause an error message if used.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( atmosphere_dim_, lat_1d_, meridian_angle_1d_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmRawRead"),
        DESCRIPTION
        (
         "Reads atmospheric data from a scenario.\n"
         "\n"
         "An atmospheric scenario includes the following data for each \n"
         "position (pressure, latitude, longitude) in the atmosphere: \n"
         "           1. temperature field \n"
         "           2. the corresponding altitude field \n"
         "           3. vmr fields for the gaseous species \n"
         "The data is stored in different files. This methods reads all \n"
         "files and creates the variables *t_field_raw*, *z_field_raw* \n"
         "\n"
         "Different atmospheric scenarios are available in arts data:\n"
         "For example tropical and midlatitude-summer. 3D scenarios are \n"
         "not available yet.\n"
         "\n"
         "Files in the scenarios look like this: tropical.H2O.xml \n"
         "\n"
         "The basename must include the path, i.e., the files can be \n"
         "anywhere, but they must be all in the same directory.\n"
         "The profile is chosen by the species name. If you have more than \n"
         "one tag group for the same species, the same profile will be \n"
         "used.\n"
         "\n"
         "Keywords: \n"
         "basename :The name and path of a particular atmospheric scenario.\n"
         "For example:\n"
         "/smiles_local/arts-data/atmosphere/fascod/tropical \n"
        ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT(t_field_raw_, z_field_raw_, vmr_field_raw_),
        INPUT(abs_species_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "basename" ),
        DEFAULTS( NODEF ),
        TYPES(    String_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "CloudboxGetIncoming" ),
        DESCRIPTION
        (
         "Calculates incoming radiation field of cloudbox by repeated\n"
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
        OUTPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, cloudbox_on_),
        INPUT( ppath_step_agenda_, rte_agenda_, iy_space_agenda_,
               surface_prop_agenda_, iy_cloudbox_agenda_,
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_, 
               t_field_, vmr_field_, r_geoid_, z_surface_, 
               cloudbox_limits_, f_grid_, stokes_dim_, 
               scat_za_grid_, scat_aa_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

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
        OUTPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, cloudbox_on_),
        INPUT( ppath_step_agenda_, rte_agenda_, iy_space_agenda_,
               surface_prop_agenda_, iy_cloudbox_agenda_,
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_, 
               t_field_, vmr_field_, r_geoid_, z_surface_, 
               cloudbox_limits_, f_grid_, stokes_dim_, 
               scat_za_grid_, scat_aa_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("cloudboxOff"),
        DESCRIPTION
        (
         "Deactivates the cloud box. \n"
         "\n"
         "The function sets *cloudbox_on* to 0, *cloudbox_limits* to be an\n"
         "empty vector and *iy_cloudbox_agenda* to an empty agenda.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( cloudbox_on_, cloudbox_limits_, iy_cloudbox_agenda_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));
  
   md_data_raw.push_back
    ( MdRecord
      ( NAME( "cloudboxSetDisort" ),
        DESCRIPTION
        (
         "For Disort calculation the cloudbox must be extended to \n"
         "cover the full atmosphere.\n"
         "This method sets *cloudbox_limits* accordingly. \n"
         ), 
        AUTHORS( "Claudia Emde" ),
        OUTPUT(cloudbox_on_, cloudbox_limits_),
        INPUT(p_grid_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));


  md_data_raw.push_back
    ( MdRecord
      ( NAME( "cloudboxSetEmpty" ),
        DESCRIPTION
        (
         "Sets the cloudbox empty for clearsky DOIT calculations. \n"
         "\n"
         "Scattering calculations using the DOIT method include \n"
         "interpolation errors. If one is interested in the cloud effect, \n"
         "should compare the DOIT result with a clearsky calculation using \n"
         "an empty cloudbox. That means that the iterative method is \n"
         "performed for a cloudbox including no particles. This method sets \n"
         "the particle number density field to zero and creates a \n"
         "dummy *scat_data_raw* structure. For a cleasky calculation, \n"
         "the methods *ParticleTypeAdd(All)* and *pnd_fieldCalc* can be \n"
         "replaced by this method. \n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( pnd_field_, scat_data_raw_),
        INPUT( p_grid_, lat_grid_, lon_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

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
         "smallest possible cloud box that encompass the given points. \n"
         "\n"
         "The points must be given in the same order as used in\n"
         "*cloudbox_limits*. That means that the first keyword argument \n"
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
         "\n"
         "Keywords: \n"
         "   p1   : Upper pressure point.\n"
         "   p2   : Lower pressure point.\n"
         "   lat1 : Lower latitude point.\n"
         "   lat2 : Upper latitude point.\n"
         "   lon1 : Lower longitude point.\n"
         "   lon2 : Upper longitude point.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( cloudbox_on_, cloudbox_limits_),
        INPUT( atmosphere_dim_, p_grid_, lat_grid_, lon_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "p1",      "p2",      "lat1",    "lat2",    "lon1",
                  "lon2" ),
        DEFAULTS( NODEF,     NODEF,     NODEF,     NODEF,     NODEF,
                  NODEF ),
        TYPES(    Numeric_t, Numeric_t, Numeric_t, Numeric_t, Numeric_t, 
                  Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "cloudboxSetManuallyAltitude" ),
        DESCRIPTION
        (
         "Sets the cloud box to encompass the given positions.\n"
         "\n"
         "The function sets *cloudbox_on* to 1 and sets *cloudbox_limits*\n"
         "following the given altitude, latitude and longitude positions.\n"
         "The index limits in *cloudbox_limits* are selected to give the\n" 
         "smallest possible cloud box that encompass the given points. \n"
         "\n"
         "The points must be given in the same order as used in\n"
         "*cloudbox_limits*. That means that altitude, latitude\n"
         "and longitude points are given in increasing order. Positions\n"
         "given for dimensions not used by the selected atmospheric\n"
         "dimensionality are ignored.\n"
         "\n"
         "The given altitude points can be outside the range of *z_field*.\n"
         "The altitude limit is then set to the end point of *p_grid*.\n"
         "The given latitude and longitude points must be inside the range\n"
         "of the corresponding grid. In addition, the latitude and longitude\n"
         "points cannot be inside the outermost grid ranges as the latitude\n"
         "and longitude limits in *cloudbox_limits* are not allowed to be\n"
         "grid end points.\n"
         "\n"
         "Keywords: \n"
         "   z1   : Lower altitude point.\n"
         "   z2   : Upper altitude point.\n"
         "   lat1 : Lower latitude point.\n"
         "   lat2 : Upper latitude point.\n"
         "   lon1 : Lower longitude point.\n"
         "   lon2 : Upper longitude point.\n"
        ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( cloudbox_on_, cloudbox_limits_),
        INPUT( atmosphere_dim_, z_field_, lat_grid_, lon_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "z1",      "z2",      "lat1",    "lat2",    "lon1",
                  "lon2" ),
        DEFAULTS( NODEF,     NODEF,     NODEF,     NODEF,     NODEF,
                  NODEF ),
        TYPES(    Numeric_t, Numeric_t, Numeric_t, Numeric_t, Numeric_t, 
                  Numeric_t )));

  

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ConvertIFToRF"),
        DESCRIPTION
        (
         "Convert *sensor_response_f* from IF to RF. The function also\n"
         "unfolds the measurement spectra *y*.\n"
         "\n"
         "Type of receiver (DSB/SSB) and main band are set by\n"
         "*sideband_mode*.\n"
         "\n"
         "This function should be used when the sensor configuration contains\n"
         "a mixer and the spectra should be given in brightness temperature.\n"
         "The reason is that the mixer converts the RF to IF, and to be able\n"
         "to perform the Rayleigh-Jeans conversion from radiance to\n"
         "brightness temperature, the radiance needs to be given in RF.\n"
         "\n"
         "Note that the number of elements in both *sensor_response_f* and\n"
         "*y* will potentially increase, since the IF is mapped to both the\n"
         "lower and upper sidebands.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( sensor_response_f_, y_ ),
        INPUT( sensor_pol_, sensor_response_za_, sensor_response_aa_, lo_,
               atmosphere_dim_, sensor_pos_, sideband_mode_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Copy"),
        DESCRIPTION
        (
         "Copy a workspace variable.\n"
         "\n"
         "This is a supergeneric method. It can copy any workspace variable\n"
         "to another workspace variable of the same group. (E.g., a Matrix to\n"
         "another Matrix.)\n"
         "\n"
         "As allways, output comes first in the argument list!\n"
         "\n"
         "Usage example:\n"
         "\n"
         "Copy(f_grid,p_grid){}\n"
         "\n"
         "Will copy the content of *p_grid* to *f_grid*. The size of *f_grid*\n"
         "is adjusted automatically (the normal behaviour for workspace\n"
         "methods).\n"
         "\n"
         "Supergeneric output:\n"
         "   Any_ : The output variable.\n"
         "\n"
         "Supergeneric input:\n"
         "   Any_ : The input variable.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Any_ ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Delete"),
        DESCRIPTION
        (
         "Deletes a workspace variable.\n"
         "\n"
         "Supergeneric input:\n"
         "   Any_     : The variable to delete.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT(  Any_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  ),
        PASSWORKSPACE(  true  ),
        PASSWSVNAMES(   true  )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ScatteringDisort"),
        DESCRIPTION
        (
         "Calls DISORT RT solver from ARTS.\n"
         "Detailed documentation to be added.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT(scat_i_p_, scat_i_lat_, scat_i_lon_, 
               f_index_, scat_data_mono_, doit_i_field1D_spectrum_),
        INPUT(cloudbox_limits_, stokes_dim_, opt_prop_part_agenda_, 
              abs_scalar_gas_agenda_, spt_calc_agenda_, 
              pnd_field_, t_field_, 
              z_field_, p_grid_, vmr_field_, scat_data_raw_, f_grid_, 
              scat_za_grid_, surface_emissivity_field_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));


   md_data_raw.push_back
    ( MdRecord
      ( NAME("DoitAngularGridsSet"),
        DESCRIPTION
        (
         "Set angular grids for DOIT calculation."
         "\n"
         "In this method the angular grids for a DOIT calculation are "
         "specified.\n"
         "For down-looking geometries it is sufficient to define\n"
         "\n"  
         "N_za_grid: number of grid points in zenith angle grid\n"
         "           recommended value: 19\n"
         "N_aa_grid: number of grid points in zenith angle grid\n"
         "           recommended value: 37\n"
         "\n"
         "From these numbers equally spaced grids are created and stored in "
         "the\n" 
         "WSVs *scat_za_grid* and *scat_aa_grid*.\n" 
         "\n"
         "For limb simulations it is important to use an optimized zenith "
         "angle \n"
         "grid with a very fine resolution about 90° for the "
         "RT calculations.\n"
         "Such a grid can be generated using *doit_za_grid_optCalc*. \n"
         "The filename of the optimized zenith angle grid can be given \n"
         "as a keyword. If a filename is given, the equidistant grid is \n"
         "taken for the calculation of the scattering integrals and the\n"
         "optimized grid is taken for the radiative transfer part.\n"
         "Otherwise, if no filename is specified \n"
         "( za_grid_opt_file = \"\" ) the equidistant grid is \n"
         "taken for the calculation of the scattering integrals and for \n"
         "the RT calculations. This option makes sense for down-looking \n"
         "cases to speed up the calculation. \n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( doit_za_grid_size_, scat_aa_grid_, scat_za_grid_),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "N_za_grid", "N_aa_grid", "za_grid_opt_file"),
        DEFAULTS( NODEF,       NODEF,       NODEF ),
        TYPES(    Index_t,     Index_t,     String_t)));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "DoitCloudboxFieldPut" ),
        DESCRIPTION
        (
         "Method for the communication between cloudbox and clearsky.\n"
         "\n"
         "This method puts the scattered radiation field into the interface\n"
         "variables between the cloudbox and the clearsky, which are \n"
         "*scat_i_p*, *scat_i_lat* and *scat_i_lon*.\n"
         "\n"
         "The best way to calculate spectra including the influence of\n" 
         "scattering is to set up the *scat_mono_agenda* where this method \n"
         "can be included.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, doit_i_field1D_spectrum_),
        INPUT( doit_i_field_, f_grid_, f_index_,   p_grid_, lat_grid_, 
               lon_grid_, scat_za_grid_, scat_aa_grid_, stokes_dim_,
               atmosphere_dim_, cloudbox_limits_, sensor_pos_, z_field_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

md_data_raw.push_back     
    ( MdRecord
      ( NAME("doit_conv_flagAbs"),
        DESCRIPTION
        (
         "Convergence test (maximum absolute difference).\n"
         "\n"
         "The function calculates the absolute differences for two " 
         "successive\n"
         "iteration fields. It picks out the maximum values for each Stokes \n"
         "component separately. The convergence test is fullfilled under the\n"
         "following conditions: \n"
         "|I(m+1) - I(m)| < epsilon_1     Intensity.\n"
         "|Q(m+1) - Q(m)| < epsilon_2     The other Stokes components.\n" 
         "|U(m+1) - U(m)| < epsilon_3    \n"
         "|V(m+1) - V(m)| < epsilon_4    \n" 
         "These conditions have to be valid for all positions in the \n"
         "cloudbox and for all directions.\n"  
         "\n"
         "The limits for convergence are set in the controlfile by \n"
         "setting the vector *epsilon* to appropriate values.\n"
         "The unit of *epsilon* is radiance [W / (m^2 Hz sr)].\n"
         "\n"
         "This method can be used in *doit_convergence_test_agenda*.\n"
        ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT(doit_conv_flag_, doit_iteration_counter_),
        INPUT(doit_conv_flag_, doit_iteration_counter_,
              doit_i_field_, doit_i_field_old_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "epsilon" ),
        DEFAULTS( NODEF ),
        TYPES(    Vector_t )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("doit_conv_flagLsq"),
        DESCRIPTION
        (
         "Convergence test (Least square).\n"
         "\n"
         "This method performs a least square convergence test for two \n"
         "successive iteration fields.\n"
         "\n"
         "The limits for convergence are set in the controlfile by \n"
         "setting the vector *epsilon* to appropriate values.\n"
         "The unit of *epsilon* is Rayleigh Jeans BT [K].\n"
         "\n"
         "Warning: This method is not recommended because this kind of \n"
         "convergence test is not sufficiently strict, so that the \n"
         "DOIT result might be wrong. \n" 
        ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT(doit_conv_flag_, doit_iteration_counter_),
        INPUT(doit_conv_flag_, doit_iteration_counter_, 
              doit_i_field_, doit_i_field_old_, f_grid_, f_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "epsilon" ),
        DEFAULTS( NODEF ),
        TYPES(    Vector_t )));
  
  md_data_raw.push_back     
    ( MdRecord
      ( NAME("doit_conv_flagAbsBT"),
        DESCRIPTION
        (
         "Convergence test (maximum absolute difference in Rayleigh Jeans "
         "BT)\n"
         "\n"
         "The function calculates the absolute differences for two " 
         "successive\n"
         "iteration fields. It picks out the maximum values for each Stokes \n"
         "component separately. The convergence test is fullfilled under the\n"
         "following conditions: \n"
         "|I(m+1) - I(m)| < epsilon_1     Intensity.\n"
         "|Q(m+1) - Q(m)| < epsilon_2     The other Stokes components.\n" 
         "|U(m+1) - U(m)| < epsilon_3    \n"
         "|V(m+1) - V(m)| < epsilon_4    \n" 
         "These conditions have to be valid for all positions in the \n"
         "cloudbox and for all directions.\n"  
         "\n"
         "The limits for convergence are set in the controlfile by \n"
         "setting the vector *epsilon* to appropriate values.\n"
         "The unit of *epsilon* is Rayleigh Jeans BT [K].\n"
         "\n"
         "This method can be used in *doit_convergence_test_agenda*.\n"
         ),
        AUTHORS( "Sreerekha T.R.", "Claudia Emde" ),
        OUTPUT(doit_conv_flag_, doit_iteration_counter_),
        INPUT(doit_conv_flag_, doit_iteration_counter_,
              doit_i_field_, doit_i_field_old_, f_grid_, f_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "epsilon" ),
        DEFAULTS( NODEF ),
        TYPES(    Vector_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "DoitInit" ),
        DESCRIPTION
        (
         "Initialize variables for DOIT scattering calculations. \n"
         "\n"
         "Before using the WSM *ScatteringDOIT*, please use this method \n"
         "to initialize the required WSVs. \n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( scat_p_index_, scat_lat_index_, scat_lon_index_, 
                scat_za_index_, scat_aa_index_, doit_scat_field_,
                doit_i_field_, doit_za_interp_, doit_is_initialized_ ),
        INPUT( stokes_dim_, atmosphere_dim_, scat_za_grid_, scat_aa_grid_,
               doit_za_grid_size_, cloudbox_limits_, scat_data_raw_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_i_fieldIterate" ),
        DESCRIPTION
        (
         "Iterative solution of the VRTE (DOIT method).\n"
         "\n"
         "A solution for the RTE with scattering is found using the\n"
         "DOIT method:\n"
         "\n"
         "1. Calculate scattering integral using *doit_scat_field_agenda*.\n"
         "2. Calculate RT with fixed scattered field using \n"
         "   *doit_rte_agenda*.\n"
         "3. Convergence test using *doit_conv_test_agenda*.\n"
         "\n"
         "Note: The atmospheric dimensionality *atmosphere_dim* can be \n"
         "      either 1 or 3. To these dimensions the method adapts \n"
         "      automatically. 2D scattering calculations are not \n"
         "      supported.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT(doit_i_field_),
        INPUT( doit_i_field_, doit_scat_field_agenda_, doit_rte_agenda_, 
               doit_conv_test_agenda_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_i_fieldSetClearsky" ),
        DESCRIPTION
        (
         "Interpolate clearsky field on all gridpoints in cloudbox. \n"
         "\n"
         "This method uses a linear 1D/3D interpolation scheme to obtain the\n"
         "radiation field on all grid points inside the cloud box from the\n"
         "clear sky field on the cloud box boundary. This radiation field\n"
         "is taken as the first guess radiation field in the DOIT module. \n"
         "\n"
         "The inputs to this method are *scat_i_p*, *scat_i_lat*\n"
         "*scat_i_lon*.  The method picks the \n"
         "monochromatic radiation field out of these variables.  The \n"
         "output of the method is the first guess field stored in the \n"
         "workspace variable *doit_i_field*.\n"
         "\n"
         "Set keyword *all_frequencies* to 1 if for each frequency the \n"
         "clearsky field should be used as initial field. Set it to 0 if \n"
         "only for the first frequency in *f_grid* the clearsky field should\n"
         "be used and for the next frequencies *doit_i_field* of the \n"
         "previous frequency should be used. Default is 1. \n"
         ),
        AUTHORS( "Sreerekha T.R. and Claudia Emde" ),
        OUTPUT(doit_i_field_),
        INPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, f_grid_, 
               f_index_, p_grid_, lat_grid_, lon_grid_, 
               cloudbox_limits_, atmosphere_dim_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("all_frequencies"),
        DEFAULTS( "1" ),
        TYPES(Index_t)));

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
         "control file by the variable doit_i_field_value, which is a\n"
         "vector containing 4 elements, the value of the initial field\n"
         "for each Stokes dimension.\n"
         "\n"
         "Output of the method is the first guess field stored in the \n"
         "workspace variable *doit_i_field*.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT(doit_i_field_),
        INPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, p_grid_, lat_grid_, 
               lon_grid_, 
               cloudbox_limits_, atmosphere_dim_, stokes_dim_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "value" ),
        DEFAULTS( NODEF ),
        TYPES(    Vector_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_i_fieldUpdate1D" ),
        DESCRIPTION
        (
         "RT calculation in cloudbox with fixed scattering integral (1D). \n"
         "\n"
         "Update the radiation field (DOIT method). The method loops\n"
         "through the cloudbox to update the radiation field for all \n"
         "positions and directions in the 1D cloudbox.\n"
         "\n"
         "Note: This method is very inefficient, because the number of \n"
         "iterations scales with the number of cloudbox pressure levels.\n"
         "It is recommended to use *doit_i_fieldUpdateSeq1D*.\n"
        ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT(doit_i_field_),
        INPUT(doit_i_field_old_, doit_scat_field_, cloudbox_limits_, 
              abs_scalar_gas_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, pnd_field_, 
              opt_prop_part_agenda_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, z_field_, r_geoid_, z_surface_,
              t_field_, f_grid_, f_index_, surface_prop_agenda_,
              doit_za_interp_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_i_fieldUpdateSeq1D" ),
        DESCRIPTION
        (
         "RT calculation in cloudbox with fixed scattering integral. \n"
         "\n"
         "Update radiation field (*doit_i_field*) in DOIT module.\n"
         "This method loops through the cloudbox to update the \n"
         "radiation field for all positions and directions in the 1D \n"
         "cloudbox. The method applies the sequential update. For more \n"
         "information refer to AUG.\n"
         "\n"
         "Note: This is the commonly used WSM for the radiation field \n"
         "update (can be used in *doit_rte_agenda*). It is recommended \n"
         "because it is the most efficient and accurate method.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT(doit_i_field_),
        INPUT(doit_i_field_, doit_scat_field_, cloudbox_limits_, 
              abs_scalar_gas_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, pnd_field_,
              opt_prop_part_agenda_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, z_field_, r_geoid_, z_surface_,
              t_field_, f_grid_, f_index_, surface_prop_agenda_,
              doit_za_interp_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_i_fieldUpdateSeq1DPP" ),
        DESCRIPTION
        (
         "RT calculation in cloudbox with fixed scattering integral. \n"
         "\n " 
         "Update radiation field (*doit_i_field*) in DOIT module.\n"
         "This method loops through the cloudbox to update the \n"
         "radiation field for all \n"
         "positions and directions in the 1D cloudbox. The method applies\n"
         "the sequential update and the plane parallel approximation.\n"
         "This method is only slightly faster than \n"
         "*doit_i_fieldUpdateSeq1D* and it is less accurate. It can not \n"
         "be used for limb simulations. \n"
         ),
        AUTHORS( "Sreerekha T.R." ),
        OUTPUT(doit_i_field_, scat_za_index_),
        INPUT(doit_scat_field_, cloudbox_limits_, 
              abs_scalar_gas_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, pnd_field_, 
              opt_prop_part_agenda_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, z_field_, r_geoid_, t_field_,
              f_grid_, f_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_i_fieldUpdateSeq3D" ),
        DESCRIPTION
        (
         "RT calculation in cloudbox with fixed scattering integral. \n"
         "\n"
         "Update radiation field (*doit_i_field*) in DOIT module.\n"
         "This method loops through the cloudbox to update the \n"
         "radiation field for all positions and directions in the 3D \n"
         "cloudbox. The method applies the sequential update. For more \n"
         "information please refer to AUG.\n"
         "Surface reflections are not yet implemented in 3D scattering \n"
         "calculations.\n"
         "\n"
         "Note: DOIT calculations in 3D are computationally expensive. \n"
         "For large 3D cloud fields it is recommended to use the Monte Carlo \n"
         "module or an independent pixel approach applying DOIT-1D. \n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT(doit_i_field_),
        INPUT(doit_i_field_, doit_scat_field_, cloudbox_limits_, 
              abs_scalar_gas_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, scat_aa_grid_,
              pnd_field_,
              opt_prop_part_agenda_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, lat_grid_, lon_grid_, z_field_,
              r_geoid_, z_surface_, t_field_,
              f_grid_, f_index_, doit_za_interp_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));
  
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_scat_fieldCalc" ),
        DESCRIPTION
        (
         "This method calculates the scattering integral field in the DOIT module.\n"
         "\n"
         "The scattering integral field is generated by integrating\n"
         "the product of phase matrix and Stokes vector over all incident\n"
         "angles. For more information please refer to AUG.\n" 
         "\n"
         "The output of this method is *doit_scat_field*\n"
         "which is used in the radiative transfer part (*doit_i_fieldUpdateXXX*).\n"
         ),
        AUTHORS( "Sreerekha T.R.", "Claudia Emde" ),
        OUTPUT( doit_scat_field_ ),
        INPUT( doit_scat_field_, pha_mat_spt_agenda_,
               doit_i_field_, pnd_field_, t_field_, atmosphere_dim_, 
               cloudbox_limits_, scat_za_grid_, scat_aa_grid_,  
               doit_za_grid_size_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_scat_fieldCalcLimb" ),
        DESCRIPTION
        (
         "This method calculates the scattering integral field in the DOIT module (Limb).\n"
         "\n"
         "The scattering integral field is the field generated by integrating\n"
         "the product of phase matrix and the Stokes vector over all incident \n"
         "angles.\n"
         "\n"
         "The output of this method is the scattering integral field "
         "*doit_scat_field*\n"
         "which is used in the radiative transfer part (*doit_i_fieldUpdateXXX*). \n"
         "For limb simulations it "
         "makes sense to use different \n"
         "zenith angle grids for the scattering integral part and the RT part, \n"
         "because the latter part requires a much finer resolution about \n"
         "90°. Taking an optimized grid for the RT part and an equidistant \n"
         "grid for the scattering integral part saves very much CPU time.\n"
         "This method uses the equidistant za_grid defined in \n"
         "*doit_angular_gridsSet* and it should always be used for limb \n"
         "simulations.\n"
         "\n"
         "For more information please refer to AUG.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( doit_scat_field_),
        INPUT( doit_scat_field_, pha_mat_spt_agenda_,
               doit_i_field_, pnd_field_, t_field_, atmosphere_dim_, 
               cloudbox_limits_, scat_za_grid_, scat_aa_grid_,  
               doit_za_grid_size_, doit_za_interp_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("DoitScatteringDataPrepare"),
        DESCRIPTION
        (
         "Prepare single scattering data for a DOIT scattering calculation.\n"
         "\n"
         "This function can be used for scattering calculations using the \n" 
         "DOIT method. \n"
         "\n"
         "First the scattering data is interpolated on the frequency using\n"
         "*scat_data_monoCalc*. Then the phase matrix data is \n"
         "transformed or interpolated from the raw data to the laboratory frame \n"
         "for all possible combinations of the angles contained in the "
         "angular\n"
         "grids which are set in *doit_angulat_gridsSet*."
         "The resultung phase matrices are \n"
         "stored in *pha_mat_sptDOITOpt*, "
         "which is used in the method\n"
         "*pha_mat_sptFromDataDOITOpt*. \n"
          ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( pha_mat_sptDOITOpt_, scat_data_mono_),
        INPUT( doit_za_grid_size_, scat_aa_grid_, scat_data_raw_, f_grid_, 
               f_index_, atmosphere_dim_, stokes_dim_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("DoitWriteIterationFields"),
        DESCRIPTION
        (
         "Write DOIT iteration fields.\n"
         "\n"
         "This method writes intermediate iteration fields \n"
         "to xml-files which is useful to check the DOIT method. \n"
         "It can be interesting to look how the radiation fields \n"
         "(*doit_i_field*) behave. The method can be used as a part of \n"
         "*convergence_test_agenda*.\n"
         "\n"
         "The keyword 'iterations' includes the numbers of the iterations\n"
         "which should be stored, e.g.,\n"
         "    'iterations = [3, 6, 9]'. \n"
         "In this case the 3rd, 6th and 9th iterations are \n"
         "stored in the files 'doit_iteration_3.xml', \n"
         "'doit_iteration_6.xml' ...\n If a number is larger than the \n"
         "total number of iterations, this number is ignored. \n"
         "\n"
         "If all iterations should be stored please set the keyword \n"
         "   'iterations = [0]'.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( ),
        INPUT( doit_iteration_counter_, doit_i_field_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "iterations" ),
        DEFAULTS( NODEF ),
        TYPES(    Array_Index_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_za_grid_optCalc" ),
        DESCRIPTION
        (
         "Zenith angle grid optimization for scattering calculation.\n"
         "\n"
         "This method optimizes the zenith angle grid. As input it requires \n"
         "a radiation field (*doit_i_field*) which is calculated on a very \n"
         "fine zenith angle grid (*scat_za_grid*).  Based on this field \n"
         "zenith angle grid points are selected, such that the maximum\n"
         "difference between the radiation field represented on the very \n"
         "fine zenith angle grid and the radiation field represented on the\n"
         "optimized grid (*doit_za_grid_opt*) is less than the accuracy \n"
         "(*acc*) provided as a keyword. The accuracy must be given in %.\n"
         "Between the grid points theradiation field is interpolated \n"
         "linearly or polynomially depending on *doit_za_interp*.\n"
         "\n"
         "Note: The method works only for a 1D atmosphere and for one \n"
         "frequency.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( doit_za_grid_opt_ ),
        INPUT( doit_i_field_, scat_za_grid_, doit_za_interp_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "acc" ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));
                                                                               
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_za_interpSet" ),
        DESCRIPTION
        (
         "Define interpolation method for zenith angle dimension.\n"
         "\n"
         "You can use this method to choose the interpolation method for \n"
         "interpolations in the zenith angle dimension. \n"
         "By default, linear interpolation is used.\n"
         "\n"
         "Keyword:\n"
         "  interp_method - 'linear' or 'polynomial' \n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( doit_za_interp_ ),
        INPUT( atmosphere_dim_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "interp_method" ),
        DEFAULTS( NODEF ),
        TYPES(    String_t )));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME("DoNothing"),
        DESCRIPTION
        (
         "As *Ignore* but for agenda output.\n"
         "\n"
         "This method is handy for use in agendas in order to suppress\n"
         "warnings about unused output workspace variables. What it does is:\n"
         "Nothing!\n"
         "\n"
         "To ensure that the variable is already set, the variable must be \n"
         "given is both global input and output. An example:\n"
         "   DoNothing(emission,emission){}\n"
         "Input and output MUST be identical.\n"
         "\n"
         "Supergeneric output:\n"
         "   Any_ : The input variable.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Any_ ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("emissionPlanck"),
        DESCRIPTION
        (
         "Emission source term for LTE.\n"
         "\n"
         "Sets *emission* for cases when emission is considered and local\n"
         "thermodynamic equilibrium is valid.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( emission_ ),
        INPUT( f_grid_, rte_temperature_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Error"),
        DESCRIPTION
        (
         "Issues an error and exits ARTS.\n"
         "\n"
         "This method can be placed in agendas that must be specified , but\n"
         "are expected not to be used for the particular case. An inclusion\n"
         "in *surface_prop_agenda* could look like:\n   "
         "Error{\"Surface interceptions of propagation path not expected.\"}\n"
         "(ignore and other dummy method calls must still be included)\n"
         "\n"
         "Keywords: \n"
         "   msg : String describing the error.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "msg" ),
        DEFAULTS( NODEF ),
        TYPES(    String_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Exit"),
        DESCRIPTION
        (
         "Stops the execution and exits ARTS.\n"
         "\n"
         "This method is handy if you want to debug one of your control\n"
         "files. You can insert it anywhere in the control file. When\n"
         "it is reached, it will terminate the program.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ext_matAddGas"),
        DESCRIPTION
        (
         "Add gas absorption to all diagonal elements of extinction matrix.\n"
         " \n"
         "The task of this method is to sum up the gas absorption of the\n"
         "different gas species and add the result to the extinction matrix.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( ext_mat_ ),
        INPUT( ext_mat_, abs_scalar_gas_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ext_matAddPart"),
        DESCRIPTION
        (
         "The particle extinction is added to *ext_mat* \n"
         "\n"
         "This function sums up the extinction matrices for all particle \n"
         "types weighted with particle number density.\n"
         "The resulting extinction matrix is added to the workspace \n"
         "variable *ext_mat* \n"
         "The output of this method is *ext_mat* (stokes_dim, stokes_dim).\n"
         "The inputs are the extinction matrix for the single particle type \n"
         "*ext_mat_spt* (part_types, stokes_dim, stokes_dim) and the local \n"
         "particle number densities for all particle types namely the \n"
         "*pnd_field* (part_types, p_grid, lat_grid, lon_grid ) for given \n"
         "*p_grid*, *lat_grid*, and *lon_grid*. The particle types required \n"
         "are specified in the control file.\n"
         ),
        AUTHORS( "Sreerekha T.R." ),
        OUTPUT( ext_mat_ ),
        INPUT( ext_mat_, ext_mat_spt_, pnd_field_, atmosphere_dim_, 
               scat_p_index_, scat_lat_index_, scat_lon_index_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME("ext_matInit"),
        DESCRIPTION
        (
         "Initialize extinction matrix.\n"
         "\n"
         "This method is necessary, because all other extinction methods just\n"
         "add to the existing extinction matrix. \n"
         "\n"
         "So, here we have to make it the right size and fill it with 0.\n"
         "\n"
         "Note, that the matrix is not really a matrix, because it has a\n"
         "leading frequency dimension.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( ext_mat_ ),
        INPUT( f_grid_, stokes_dim_, f_index_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("f_gridFromSensor"),
        DESCRIPTION
        (
         "Automatically calculate f_grid to match the sensor.\n"
         "\n"
         "This method is handy if you are simulating an AMSU-type instrument,\n"
         "consisting of a few discrete channels.\n"
         "\n"
         "It calculates f_grid to match the instrument, as given by the local\n"
         "oscillator frequencies *lo*, the backend frequencies *f_backend*, and\n"
         "the backend channel responses *backend_channel_response*.\n"
         "\n"
         "You have to specify the desired spacing in the keyword *spacing*, which\n"
         "has a default value of 100 MHz. (The actual value is 0.1e9, since our\n"
         "unit is Hz.)\n"
         "\n"
         "The produced grid will not have exactly the requested spacing, but\n"
         "will not be coarser than requested. The algorithm starts with the band\n"
         "edges, then adds additional points until the spacing is at least as\n"
         "fine as requested.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( f_grid_ ),
        INPUT( lo_, f_backend_, backend_channel_response_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "spacing" ),
        DEFAULTS( ".1e9"),
        TYPES(    Numeric_t )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("FlagOff"),
        DESCRIPTION
        (
         "Sets an index variable that acts as an on/off flag to 0.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Index_ ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("FlagOn"),
        DESCRIPTION
        (
         "Sets an index variable that acts as an on/off flag to 1.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Index_ ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("f_gridFromGasAbsLookup"),
        DESCRIPTION
        (
         "Sets *f_grid* to the frequency grid of *abs_lookup*.\n"
         "\n"
         "Must be called between importing/creating raw absorption table and\n"
         "call of *abs_lookupAdapt*.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( f_grid_  ),
        INPUT(  abs_lookup_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("GaussianResponse"),
        DESCRIPTION
        (
         "Creates a matrix with a relative grid and a gaussian response.\n"
         "\n"
         "The grid is a relative grid so therefore it is centred around\n"
         "zero, the TotWidth keyword describes the difference between the\n"
         "maximum grid point and the minimum. The grid range is then divided\n"
         "into grid points equally spaced with max distance equal or less\n"
         "than MaxSpacing. The results are the stored in columns in the matrix.\n"
         "Grid points in the first and values in the second column.\n"
         "\n"
         "Generic output: \n"
         "   Matrix     : The matrix containing the grid and response values.\n"
         "\n"
         "Keywords:\n"
         "   fwhm       : The Full Width at Half Mean value for the response.\n"
         "   tot_width   : The total width of the relative grid.\n"
         "   max_spacing : The maximum step between grid points.\n"
         ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Matrix_ ),
        GINPUT( ),
        KEYWORDS( "fwhm",    "tot_width", "max_spacing" ),
        DEFAULTS( NODEF,     NODEF,       NODEF ),
        TYPES(    Numeric_t, Numeric_t,   Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("GriddedField4ExtractFromArrayOfGriddedField4"),
        DESCRIPTION
        (
         "Extract a *GriddedField4* from an array of *GriddedField4*.\n"
         "This is useful for example for extracting *atm_fields_compact*\n"
         "from *batch_atm_fields_compact*.\n"
         "\n"
         "Copies element with given *Index* from input array\n"
         "to create output *GriddedField4*.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( GriddedField4_ ),
        GINPUT(  ArrayOfGriddedField4_, Index_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Ignore"),
        DESCRIPTION
        (
         "Ignore a workspace variable.\n"
         "\n"
         "This method is handy for use in agendas in order to suppress warnings\n"
         "about unused input workspace variables. What it does is: Nothing!\n"
         "In other words, it just ignores the variable it is called on.\n"
         "\n"
         "This is a supergeneric method. It can ignore any workspace variable\n"
         "you want.\n"
         "\n"
         "Usage example:\n"
         "\n"
         "AgendaSet(els_agenda){\n"
         "  Ignore(ls_sigma){}\n"
         "  elsLorentz{}\n"
         "}\n"
         "\n"
         "Without Ignore you would get an error message, because els_agenda is\n"
         "supposed to use the Doppler width *ls_sigma*, but the Lorentz lineshape\n"
         "*elsLorentz* does not need it.\n"
         "\n"
         "Supergeneric input:\n"
         "   Any_ : The input variable.\n"
         ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT(  Any_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("INCLUDE"),
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
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "IndexCreate" ),
        DESCRIPTION
        (
         "Creates an Index variable.\n"
         "\n"
         "If the variable already exists, it'll be reset.\n"
         "\n"
         "Generic output: \n"
         "   Index: New Index variable.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Index_ ),
        GINPUT(),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("IndexSet"),
        DESCRIPTION
        (
         "Sets an index workspace variable to the given value. \n"
         "\n"
         "Generic output: \n"
         "   Index : The index variable to be set. \n"
         "\n"
         "Keywords:\n"
         "   value : A positive integer.\n" 
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Index_ ),
        GINPUT( ),
        KEYWORDS( "value" ),
        DEFAULTS( NODEF ),
        TYPES(    Index_t )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("IndexStep"),
        DESCRIPTION
        (
         "Performs GOUTPUT = GINPUT + 1\n"
         "\n"
         "Input and output can be same variable.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Index_ ),
        GINPUT( Index_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("InterpAtmFieldToRteGps"),
        DESCRIPTION
        (
         "Scalar interpolation of atmospheric fields.\n" 
         "\n"
         "The position is specified by the combination of *rte_gp_p*, \n"
         "*rte_gp_lat* and *rte_gp_lon*.\n"
         "\n"
         "Generic output: \n"
         "   Numeric : Value obtained by interpolation. \n"
         "\n"
         "Generic input:\n"
         "   Tensor3 : Field to interpolate.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, 
               rte_gp_p_, rte_gp_lat_, rte_gp_lon_ ),
        GOUTPUT( Numeric_ ),
        GINPUT( Tensor3_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));
  
  md_data_raw.push_back     
    ( MdRecord
      ( NAME("InterpSurfaceFieldToRteGps"),
        DESCRIPTION
        (
         "Scalar interpolation of surface fields.\n" 
         "\n"
         "The position is specified by the combination of *rte_gp_lat* and \n"
         "*rte_gp_lon*.\n"
         "\n"
         "Generic output: \n"
         "   Numeric : Value obtained by interpolation. \n"
         "\n"
         "Generic input:\n"
         "   Matrix : Field to interpolate.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( atmosphere_dim_, lat_grid_, lon_grid_, 
               rte_gp_lat_, rte_gp_lon_ ),
        GOUTPUT( Numeric_ ),
        GINPUT( Matrix_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));
  
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
         "The intensity field is interpolated to the position and direction\n"
         "given (specified by *rte_XXX*). A linear interpolation is used for\n"
         "all dimensions.\n"
         "\n"
         "The intensity field on the cloux box boundaries is provided by\n"
         "*scat_i_p/lat/lon* and these variables are interpolated if the.\n"
         "given position is at any boundary.\n"
         "\n"
         "Interpolation of the internal field is not yet possible.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( iy_ ),
        INPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, doit_i_field1D_spectrum_,
               rte_gp_p_, rte_gp_lat_, rte_gp_lon_, rte_los_,  cloudbox_on_,
               cloudbox_limits_, atmosphere_dim_, stokes_dim_, scat_za_grid_,
               scat_aa_grid_, f_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

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
        OUTPUT( iy_ ),
        INPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, doit_i_field1D_spectrum_,
               rte_gp_p_, rte_gp_lat_,
               rte_gp_lon_, rte_los_,  cloudbox_on_, cloudbox_limits_,
               atmosphere_dim_, stokes_dim_, scat_za_grid_, scat_aa_grid_, 
               f_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianAddAbsSpecies" ),
        DESCRIPTION
        (
         "Adds a absorption species as a retrieval quantity to the Jacobian.\n"
         "\n"
         "For 1D or 2D calculations the latitude and/or longitude grid of\n"
         "the retrieval field should be set to zero length.\n"
         "\n"
         "There are two possible calculation methods:\n"
         "   \"analytical\"   : (semi-)analytical expressions are used\n"
         "   \"perturbation\" : pure numerical perturbations are used\n"
         "\n"
         "The retrieval unit can be:\n"
         "   \"vmr\" : volume mixing ratio \n"
         "   \"nd\"  : number density \n"
         "   \"rel\" : relative unit (e.g. 1.1 means 10% more of the gas) \n"
         "\n"
         "For perturbation calculations the size of the perturbation is set\n"
         "by the user. The unit of the perturbation size is identical to \n"
         "the retrieval unit.\n"
         "\n"
         "Generic input:\n"
         "  Vector : The pressure grid of the retrieval field.\n"
         "  Vector : The latitude grid of the retrieval field.\n"
         "  Vector : The longitude grid of the retrieval field.\n"
         "\n"
         "Keywords:\n"
         "  species : The SpeciesTag of the retrieval quantity.\n"
         "  method  : Calculation method. See above.\n"
         "  unit    : Retrieval unit. See above.\n"
         "  dx      : Size of perturbation.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( jacobian_quantities_, jacobian_agenda_ ),
        INPUT( jacobian_, atmosphere_dim_, p_grid_, lat_grid_, lon_grid_ ),
        GOUTPUT( ),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS( "species", "method", "unit",   "dx" ),
        DEFAULTS( NODEF,     NODEF,    NODEF,    NODEF),
        TYPES(    String_t,  String_t, String_t, Numeric_t )));
         
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
        OUTPUT( jacobian_quantities_, jacobian_agenda_ ),
        INPUT( jacobian_, atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, 
               pnd_field_, pnd_field_perturb_, cloudbox_limits_ ),
        GOUTPUT( ),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));
         
  md_data_raw.push_back
    ( MdRecord
      ( NAME("jacobianAddPointing"),
        DESCRIPTION
        (
         "Add a pointing as a retrieval quantity to the Jacobian.\n"
         "\n"
         "This function adds a pointing offset described over time by a\n"
         "polynomial or a gitter. By setting the polynomial order to -1,\n"
         "each spectra is treated separately. The WSV *sensor_time* is\n"
         "used to assign a timestamp for each sensor position.\n"
         "\n"
         "The WSM *jacobianCalcPointing is automatically added to\n"
         "*jacobian_agenda*.\n"
         "\n"
         "The perturbation is defined as an absolute perturbation, this\n"
         "perturbation is then applied to the sensor line-of-sight angles.\n"
         "\n"
         "The unit of the Jacobian is the unit of *y* per degree.\n"
         "\n"
         "NOTE: So far this function only treats zenith angle offsets.\n"
         "\n"
         "Keywords:\n"
         "  dza                 : The size of the perturbation.\n"
         "  poly_order          : Order of the polynomial.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( jacobian_quantities_, jacobian_agenda_ ),
        INPUT( jacobian_, sensor_pos_, sensor_time_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "dza",     "poly_order" ),
        DEFAULTS( NODEF,     NODEF ),
        TYPES(    Numeric_t, Index_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianAddTemperature" ),
        DESCRIPTION
        (
         "Add the temperature as a retrieval quantity to the Jacobian.\n"
         "\n"
         "For 1D or 2D calculations the latitude and/or longitude grid of\n"
         "the retrieval field should be set to zero length.\n"
         "The WSM *jacobianCalcTemperature* is automatically added to\n"
         "*jacobian_agenda*.\n"
         "\n" 
         "The perturbation can either be given as an absolute or a relative\n"
         "perturbation, this perturbation is then applied to the temperature\n"
         "field at each retrieval grid point.\n"
         "\n"
         "Unit of the Jacobian is the unit of *y* per Kelvin.\n"
         "\n"
         "NOTE: So far only \"perturbation\" method implemented without\n"
         "hydrostatic equilibrium.\n"
         "\n"
         "Generic input:\n"
         "  Vector : The pressure grid of the retrieval field.\n"
         "  Vector : The latitude grid of the retrieval field.\n"
         "  Vector : The longitude grid of the retrieval field.\n"
         "\n"
         "Keywords:\n"
         "  hse     : \"on\" or \"off\".\n"
         "  method  : \"analytic\" or \"perturbation\".\n"
         "  dx      : Size of perturbation.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( jacobian_quantities_, jacobian_agenda_ ),
        INPUT( jacobian_, atmosphere_dim_, p_grid_, lat_grid_, lon_grid_ ),
        GOUTPUT( ),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS( "hse",    "method", "dx" ),
        DEFAULTS( NODEF,    NODEF,    NODEF ),
        TYPES(    String_t, String_t, Numeric_t )));
  
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianCalc" ),
        DESCRIPTION
        (
         "Executes *jacobian_agenda* to calculate (parts of) *jacobian*."
         "\n"
         "It is important that *y* holds the original output of *RteCalc*\n"
         "as the methods called performs perturbations to obtain cahnges in\n"
         "*y*.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( jacobian_ ),
        INPUT( jacobian_agenda_, jacobian_indices_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));
        
  md_data_raw.push_back
    ( MdRecord
      ( NAME("jacobianCalcAbsSpecies"),
        DESCRIPTION
        (
        "Calculates absorption species jacobians by perturbations.\n"
        "\n"
        "This function is added to *jacobian_agenda* by jacobianAddAbsSpecies\n"
        "and should normally not be called by the user.\n"
        ),
        AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
        OUTPUT( jacobian_ ),
        INPUT( y_, jacobian_quantities_, jacobian_indices_, abs_species_, 
               ppath_step_agenda_, 
               rte_agenda_, iy_space_agenda_, surface_prop_agenda_, 
               iy_cloudbox_agenda_, atmosphere_dim_, p_grid_, lat_grid_, 
               lon_grid_, z_field_, t_field_, vmr_field_, 
               r_geoid_, z_surface_, cloudbox_on_, 
               cloudbox_limits_, sensor_response_, sensor_pos_, sensor_los_, 
               f_grid_, stokes_dim_, antenna_dim_, mblock_za_grid_, 
               mblock_aa_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "species" ),
        DEFAULTS( NODEF ),
        TYPES(    String_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("jacobianCalcParticle"),
        DESCRIPTION
        (
        "Calculates particle number densities jacobians by perturbations\n"
        "\n"
        "This function is added to *jacobian_agenda* by jacobianAddParticle\n"
        "and should normally not be called by the user.\n"
        ),
        AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
        OUTPUT( jacobian_ ),
        INPUT( y_, jacobian_quantities_, jacobian_indices_, pnd_field_perturb_,
               jacobian_particle_update_agenda_,
               ppath_step_agenda_, rte_agenda_, iy_space_agenda_, 
               surface_prop_agenda_, iy_cloudbox_agenda_, atmosphere_dim_, 
               p_grid_, lat_grid_, lon_grid_, z_field_, t_field_, vmr_field_,
               r_geoid_, z_surface_, 
               cloudbox_on_, cloudbox_limits_, pnd_field_,
               sensor_response_, sensor_pos_, sensor_los_, f_grid_, 
               stokes_dim_, antenna_dim_, mblock_za_grid_, mblock_aa_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));
        
  md_data_raw.push_back
    ( MdRecord
      ( NAME("jacobianCalcPointing"),
        DESCRIPTION
        (
        "Calculates pointing deviation jacobians by perturnbations.\n"
        "\n"
        "This function is added to *jacobian_agenda* by jacobianAddPointing\n"
        "and should normally not be called by the user.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( jacobian_ ),
        INPUT( y_, jacobian_quantities_, jacobian_indices_, 
               sensor_time_, ppath_step_agenda_, 
               rte_agenda_, iy_space_agenda_, surface_prop_agenda_, 
               iy_cloudbox_agenda_, atmosphere_dim_, p_grid_, lat_grid_, 
               lon_grid_, z_field_, t_field_, vmr_field_,
               r_geoid_, z_surface_, cloudbox_on_, 
               cloudbox_limits_, sensor_response_, sensor_pos_, sensor_los_, 
               f_grid_, stokes_dim_, antenna_dim_, mblock_za_grid_, 
               mblock_aa_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("jacobianCalcTemperature"),
        DESCRIPTION
        (
        "Calculates temperature jacobians by perturbations..\n"
        "\n"
        "This function is added to *jacobian_agenda* by jacobianAddTemperature\n"
        "and should normally not be called by the user.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( jacobian_ ),
        INPUT( y_, jacobian_quantities_, jacobian_indices_, ppath_step_agenda_,
               rte_agenda_, 
               iy_space_agenda_, surface_prop_agenda_, iy_cloudbox_agenda_, 
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_, 
               t_field_, vmr_field_, r_geoid_, z_surface_, 
               cloudbox_on_, cloudbox_limits_, 
               sensor_response_, sensor_pos_, sensor_los_, f_grid_, 
               stokes_dim_, antenna_dim_, mblock_za_grid_, mblock_aa_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("jacobianClose"),
        DESCRIPTION
        (
         "Close the array of retrieval quantities and prepare for calculation\n"
         "of the Jacobian matrix.\n"
         "\n"
         "This function closes the *jacobian_quantities* array and sets the\n"
         "correct size of *jacobian*. Retrieval quantities should not be\n"
         "added after a call to this WSM.\n"
         "\n"
         "To define the final *jacobian* the number of spectra is needed.\n"
         "Therefor the number of measurement blocks, taken from *sensor_pos*\n"
         "and the size of *sensor_response* has to be defined.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( jacobian_, jacobian_indices_ ),
        INPUT( jacobian_quantities_, sensor_pos_, sensor_response_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("jacobianInit"),
        DESCRIPTION
        (
         "Initialise the variables connected to the Jacobian matrix.\n"
         "\n"
         "This function initialises the *jacobian_quantities* array so\n"
         "that retrieval quantities can be added to it, therefor it has\n"
         "to be called before any subsequent calls to jacobianAddGas\n"
         "jacobianAddTemperature or jacobianAddPointing.\n"
         "\n"
         "The Jacobian quantities are initialised to be empty.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( jacobian_, jacobian_quantities_, jacobian_indices_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("jacobianOff"),
        DESCRIPTION
        (
         "Makes mandatory initialisation of some jacobian variables.\n"
         "\n"
         "Some jacobian WSVs must be initilised even if no such calculations\n"
         "will be performed and this is handled with this method. That is,\n"
         "this method must be called when no jacobians will be calculated.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( jacobian_, jacobian_quantities_, jacobian_indices_, 
                jacobian_unit_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianUnit" ),
        DESCRIPTION
        (
         "Conversion of *jacobian* to other spectral units.\n"
         "\n"
         "Works as *yUnit* but operates on *jacobian* and conversion\n "
         "determined by *jacobian_unit*.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT(jacobian_ ),
        INPUT( jacobian_, jacobian_unit_, y_unit_, sensor_pos_, sensor_los_, 
               sensor_response_f_, sensor_response_za_, sensor_response_aa_,
               sensor_response_pol_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("MatrixCBR"),
        DESCRIPTION
        (
         "Sets a matrix to hold cosmic background radiation (CBR).\n"
         "\n"
         "The CBR is assumed to be un-polarized and Stokes components 2-4\n"
         "are zero. Number of Stokes components, that equals the number \n"
         "of columns in the created matrix, is determined by *stokes_dim. \n"
         "The number of rows in the created matrix equals the length of the \n"
         "given frequency vector. \n"
         "\n"
         "The cosmic radiation is modelled as blackbody radiation for the \n"
         "temperature given by the global constant COSMIC_BG_TEMP, set in \n"
         "the file constants.cc. The frequencies are taken from the generic \n"
         "input vector:\n"
         "   MatrixCBR(iy_space,f_grid){} \n"
         "\n"
         "Generic output: \n"
         "   Matrix : Matrix with cosmic background radiation. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A set of frequencies.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( stokes_dim_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MatrixCreate" ),
        DESCRIPTION
        (
         "Creates an empty Matrix.\n"
         "\n"
         "If the variable already exists, it'll be reset.\n"
         "\n"
         "Generic output: \n"
         "   Matrix: New empty Matrix.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT(),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixMatrixMultiply"),
        DESCRIPTION
        (
         "Multiply a Matrix with another Matrix and store the result in the result\n"
         "Matrix.\n"
         "\n"
         "This just computes the normal Matrix-Matrix product, Y=M*X. It is ok\n"
         "if Y and X are the same Matrix. This function is handy for\n"
         "multiplying the H Matrix to batch spectra.\n"
         "\n"
         "Generic output:\n"
         "   Matrix : The result of the multiplication (dimension mxc).\n"
         "\n"
         "Generic input:\n"
         "   Matrix : The Matrix to multiply (dimension mxn).\n"
         "   Matrix : The original Matrix (dimension nxc).\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Matrix_ ),
        GINPUT( Matrix_, Matrix_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Matrix1ColFromVector"),
        DESCRIPTION
        (
         "Forms a matrix containing 1 column from a vector.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Matrix2ColFromVectors"),
        DESCRIPTION
        (
         "Forms a matrix containing 2 columns from two vectors.\n"
         "\n"
         "The vectors are put as columns in the matrix in the same order\n"
         "as they are given.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied into the first column. \n"
         "   Vector : The vector to be copied into the second column.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_, Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Matrix3ColFromVectors"),
        DESCRIPTION
        (
         "Forms a matrix containing 3 columns from three vectors.\n"
         "\n"
         "The vectors are put as columns in the matrix in the same order\n"
         "as they are given.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied into the first column. \n"
         "   Vector : The vector to be copied into the second column. \n"
         "   Vector : The vector to be copied into the third column.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Matrix1RowFromVector"),
        DESCRIPTION
        (
         "Forms a matrix containing 1 row from a vector.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Matrix2RowFromVectors"),
        DESCRIPTION
        (
         "Forms a matrix containing 2 rows from two vectors.\n"
         "\n"
         "The vectors are put as rows in the matrix in the same order\n"
         "as they are given.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied into the first row. \n"
         "   Vector : The vector to be copied into the second row.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_, Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Matrix3RowFromVectors"),
        DESCRIPTION
        (
         "Forms a matrix containing 3 rows from three vectors.\n"
         "\n"
         "The vectors are put as rows in the matrix in the same order\n"
         "as they are given.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied into the first row. \n"
         "   Vector : The vector to be copied into the second row. \n"
         "   Vector : The vector to be copied into the third row.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixExtractFromArrayOfMatrix"),
        DESCRIPTION
        (
         "Extract a Matrix from an array of matrices.\n"
         "\n"
         "Copies *Matrix* with given Index from input *ArrayOfMatrix*\n"
         "variable to create output *Matrix*.\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Matrix_ ),
        GINPUT(  ArrayOfMatrix_, Index_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixExtractFromTensor3"),
        DESCRIPTION
        (
         "Extract a Matrix from a Tensor3.\n"
         "\n"
         "Copies page with given Index from input Tensor3 variable to create \n"
         "output Matrix.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Matrix_ ),
        GINPUT(  Tensor3_, Index_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixPlanck"),
        DESCRIPTION
        (
         "Sets a matrix to hold blackbody radiation.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : Matrix with blackbody radiation. \n"
         "\n"
         "Generic input: \n"
         "   Vector  : A set of frequencies. \n"
         "   Numeric : Blackbody temperature. \n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( stokes_dim_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_, Numeric_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixScale"),
        DESCRIPTION
        (
         "Scales all elements of a matrix with the same value. \n"
         "\n"
         "The result can either be stored in the same or another matrix. \n"
         "\n"
         "Generic output: \n"
         "   Matrix : Return matrix. \n"
         "\n"
         "Generic input: \n"
         "   Matrix : Original matrix. \n"
         "\n"
         "Keywords: \n"
         "   value : The value to be multiplied with the matrix.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Matrix_ ),
        GINPUT( Matrix_ ),
        KEYWORDS( "value" ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t   )));

   md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixSet"),
        DESCRIPTION
        (
         "Creates a matrix and sets all elements to the specified value.\n"
         "The size is determined by the variables *ncols* and *nrows*.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Keywords:\n"
         "   value : The value of the matrix elements.\n" 
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( nrows_, ncols_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( ),
        KEYWORDS( "value"   ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MatrixToPlanckBT" ),
        DESCRIPTION
        (
         "Converts a matrix of radiances to brightness temperatures by \n"
         "inverting the Planck function.\n"
         "\n"
         "This function works as *MatrixToRJBT*. However, this function \n"
         "is not recommended in connection with inversions, but can be used \n"
         "to display calculated spectra in a temperature scale.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : A matrix with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Matrix : A matrix with radiance values.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( sensor_pos_, sensor_los_, sensor_response_f_,
               sensor_response_za_, sensor_response_aa_,
               sensor_response_pol_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Matrix_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MatrixToRJBT" ),
        DESCRIPTION
        (
         "Converts a matrix of radiances to brightness temperatures by \n"
         "the Rayleigh-Jeans approximation of the Planck function.\n"
         "\n"
         "This function works as *VectorToRJBT*, but operates on a matrix.\n"
         "Each column of the matrix is treated to contain a spectral vector,\n"
         "with frequencies repeated as assumed in *VectorToRJBT*. \n"
         "\n"
         "Generic output: \n"
         "   Matrix : A matrix with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Matrix : A matrix with radiance values.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( sensor_pos_, sensor_los_, sensor_response_f_,
               sensor_response_za_, sensor_response_aa_,
               sensor_response_pol_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Matrix_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixUnitIntensity"),
        DESCRIPTION
        (
         "Sets a matrix to hold unpolarised radiation with unit intensity.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : Matrix with unit radiation. \n"
         "\n"
         "Generic input: \n"
         "   Vector  : A set of frequencies.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( stokes_dim_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("mc_antennaSetGaussian"),
        DESCRIPTION
        (
         "Makes mc_antenna (used by MCGeneral) a 2D Gaussian pattern.\n"
         "\n"
         "The gaussian antenna pattern is determined by the keyword parameters\n"
         "za_sigma, and aa_sigma, which represent the standard deviations in the\n"
         "uncorrelated bivariate normal distribution\n"
        ),
        AUTHORS( "Cory Davis" ),
        OUTPUT( mc_antenna_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "za_sigma", "aa_sigma" ),
        DEFAULTS( NODEF,      NODEF ),
        TYPES(    Numeric_t,  Numeric_t)));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("mc_antennaSetGaussianByFWHM"),
        DESCRIPTION
        (
         "Makes mc_antenna (used by MCGeneral) a 2D Gaussian pattern.\n"
         "\n"
         "The gaussian antenna pattern is determined by the keyword parameters\n"
         "za_fwhm, and aa_fwhm, which represent the full width half maximum (FWHM)\n"
         "of the antenna response, in the zenith and azimuthal planes.\n"
        ),
        AUTHORS( "Cory Davis" ),
        OUTPUT( mc_antenna_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "za_fwhm", "aa_fwhm" ),
        DEFAULTS( NODEF,     NODEF ),
        TYPES(    Numeric_t, Numeric_t)));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("mc_antennaSetPencilBeam"),
        DESCRIPTION
        (
         "Makes mc_antenna (used by MCGeneral) a pencil beam.\n"
         "\n"
         "This WSM makes the subsequent MCGeneral WSM perform pencil beam\n"
         "RT calculations.\n" 
        ),
        AUTHORS( "Cory Davis" ),
        OUTPUT( mc_antenna_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("mc_IWP_cloud_opt_pathCalc"),
        DESCRIPTION
        (
         "Calculates the FOV averaged ice water path and cloud optical path\n"
         "for a given viewing direction\n"
        ),
        AUTHORS( "Cory Davis" ),
        OUTPUT( mc_IWP_, mc_cloud_opt_path_, mc_IWP_error_, mc_cloud_opt_path_error_, 
                mc_iteration_count_),
        INPUT( mc_antenna_, sensor_pos_, sensor_los_, ppath_step_agenda_, p_grid_, 
               lat_grid_, lon_grid_, r_geoid_, z_surface_, z_field_, t_field_, vmr_field_, 
               cloudbox_limits_, pnd_field_, scat_data_mono_, particle_masses_, mc_seed_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "max_iter" ),
        DEFAULTS( NODEF ),
        TYPES(    Index_t )));
  
  md_data_raw.push_back     
    ( MdRecord
      ( NAME("MCGeneral"),
        DESCRIPTION
        ("A generalised 3D reversed Monte Carlo radiative algorithm, that \n"
         "allows for 2D antenna patterns, surface reflection and arbitrary\n"
         "sensor positions.\n\n"
         "The main output variables *y* and *mc_error* represent the \n"
         "Stokes vector integrated over the antenna function, and the \n"
         "estimated error in this vector respectively.\n"
         "The keyword parameter `maxiter\' describes the number of `photons\'\n"
         "used in the simulation (more photons means smaller *mc_error*).\n"
         "std_err is the desired value of mc_error, and max_time is the maximum\n"
         "allowed number of seconds for MCGeneral.  MCGeneral\n"
         "will terminate once any of the max_iter, std_err, max_time criteria are\n"
         "met.  If negative values are given for these parameters then it is ignored.\n"
         " Negative values of rng_seed seed the random number generator \n "
         "according to system time, positive rng_seed values are taken literally.\n"),
        AUTHORS( "Cory Davis" ),
        OUTPUT( y_, mc_iteration_count_, mc_error_, mc_points_ ),
        INPUT( mc_antenna_, f_grid_, f_index_, sensor_pos_, sensor_los_, 
               stokes_dim_, iy_space_agenda_, surface_prop_agenda_, 
               opt_prop_gas_agenda_, abs_scalar_gas_agenda_, p_grid_, lat_grid_, 
               lon_grid_, z_field_, r_geoid_, z_surface_, t_field_, vmr_field_, 
               cloudbox_limits_, pnd_field_, scat_data_mono_, 
               mc_seed_, y_unit_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "std_err", "max_time", "max_iter", "z_field_is_1D" ),
        DEFAULTS( NODEF,     NODEF,      NODEF,      NODEF ),
        TYPES(    Numeric_t, Index_t,    Index_t,    Index_t )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("MCIPA"),
        DESCRIPTION
        ("A specialised 3D reversed Monte Carlo radiative algorithm, that \n"
         "mimics independent pixel appoximation simulations .  Probably temporary.\n"),
        AUTHORS( "Cory Davis" ),
        OUTPUT( y_, mc_iteration_count_, mc_error_, mc_points_),
        INPUT( mc_antenna_, f_grid_, f_index_, sensor_pos_, sensor_los_, 
               stokes_dim_, iy_space_agenda_, surface_prop_agenda_, 
               opt_prop_gas_agenda_, abs_scalar_gas_agenda_, ppath_step_agenda_,
               p_grid_, lat_grid_, lon_grid_, z_field_, r_geoid_, z_surface_, 
               t_field_, vmr_field_, cloudbox_limits_, pnd_field_, 
               scat_data_mono_, mc_seed_, y_unit_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "std_err", "max_time", "max_iter", "z_field_is_1D" ),
        DEFAULTS( NODEF,     NODEF,      NODEF,      NODEF ),
        TYPES(    Numeric_t, Index_t,    Index_t,    Index_t )));

  /* Removed as ScatteringMonteCarlo is not working
  md_data_raw.push_back     
    ( MdRecord
      ( NAME("MCSetIncomingEmpty"),
        DESCRIPTION
        ("Sets mc_incoming to be empty.  \n"
         "This is needed when using ScatteringMonteCarlo with incoming_lookup=0\n"
        ),
        AUTHORS( "Cory Davis" ),
        OUTPUT( mc_incoming_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));
  */

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("MCSetSeedFromTime"),
        DESCRIPTION
        ("Sets the value of mc_seed from system time\n"),
        AUTHORS( "Cory Davis" ),
        OUTPUT( mc_seed_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "NumericCreate" ),
        DESCRIPTION
        (
         "Creates a Numeric variable.\n"
         "\n"
         "If the variable already exists, it'll be reset.\n"
         "\n"
         "Generic output: \n"
         "   Numeric: New Numeric variable.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Numeric_ ),
        GINPUT(),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("NumericSet"),
        DESCRIPTION
        (
         "Sets a numeric workspace variable to the given value. \n"
         "\n"
         "Generic output: \n"
         "   Numeric : The numeric variable to be set. \n"
         "\n"
         "Keywords:\n"
         "   value : The value.\n" 
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Numeric_ ),
        GINPUT(),
        KEYWORDS( "value"   ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("nelemGet"),
        DESCRIPTION
        (
         "Retrieve nelem from given variable and store the value in the \n"
         "workspace variable *nelem*\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( nelem_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("ncolsGet"),
        DESCRIPTION
        (
         "Retrieve ncols from given variable and store the value in the \n"
         "workspace variable *ncols*\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( ncols_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("nrowsGet"),
        DESCRIPTION
        (
         "Retrieve nrows from given variable and store the value in the \n"
         "workspace variable *nrows*\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( nrows_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("npagesGet"),
        DESCRIPTION
        (
         "Retrieve npages from given variable and store the value in the \n"
         "workspace variable *npages*\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( npages_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("nbooksGet"),
        DESCRIPTION
        (
         "Retrieve nbooks from given variable and store the value in the \n"
         "workspace variable *nbooks*\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( nbooks_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("nshelvesGet"),
        DESCRIPTION
        (
         "Retrieve nshelves from given variable and store the value in the \n"
         "workspace variable *nshelves*\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( nshelves_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("nvitrinesGet"),
        DESCRIPTION
        (
         "Retrieve nvitrines from given variable and store the value in the \n"
         "workspace variable *nvitrines*\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( nvitrines_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("nlibrariesGet"),
        DESCRIPTION
        (
         "Retrieve nlibraries from given variable and store the value in the \n"
         "workspace variable *nlibraries*\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( nlibraries_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("NumericExtractFromVector"),
        DESCRIPTION
        (
         "Extract a Numeric from a Vector.\n"
         "\n"
         "Copies element with given Index from input Vector variable to \n"
         "create output Numeric.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Numeric_ ),
        GINPUT(  Vector_, Index_ ),
        KEYWORDS(),
        DEFAULTS( ),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("opt_prop_sptFromData"),
        DESCRIPTION
        (
         "Calculates opticle properties for the single particle types.\n"
         "\n"
         "In this function extinction matrix and absorption vector are \n"
         "calculated in the laboratory frame. These properties are required\n"
         "for the RT calculation, inside the the i_fieldUpdateXXX \n"
         "functions.\n" 
         "\n"
         "The interpolation of the data on the actual frequency is the \n"
         "first step in this function. \n"
         "\n"
         "Then the transformation from the database coordinate system to to \n"
         "laboratory coordinate system is done.\n"
         "\n"
         "Output of the function are *ext_mat_spt*, and *abs_vec_spt* which\n"
         "hold the optical properties for a specified propagation direction\n"
         "for each particle type. \n"
        ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( ext_mat_spt_, abs_vec_spt_ ),
        INPUT(  ext_mat_spt_, abs_vec_spt_, scat_data_raw_,
                scat_za_grid_, 
                scat_aa_grid_, scat_za_index_, scat_aa_index_, 
                f_index_, f_grid_, rte_temperature_,
                pnd_field_, scat_p_index_, scat_lat_index_, scat_lon_index_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("opt_prop_sptFromMonoData"),
        DESCRIPTION
        (
         "Calculates optical properties for the single particle types.\n"
         "\n"
         "In this function extinction matrix and absorption vector are \n"
         "calculated in the laboratory frame. "
         "\n"
         "The single scattering data is obtained from scat_data_mono, so\n"
         "frequency interpolation is not required\n"
         "\n"
         "Output of the function are *ext_mat_spt*, and *abs_vec_spt* which\n"
         "hold the optical properties for a specified propagation direction\n"
         "for each particle type. \n"
        ),
        AUTHORS( "Cory Davis" ),
        OUTPUT( ext_mat_spt_, abs_vec_spt_ ),
        INPUT(  ext_mat_spt_, abs_vec_spt_, scat_data_mono_,
                scat_za_grid_, 
                scat_aa_grid_, scat_za_index_, scat_aa_index_, rte_temperature_,
                pnd_field_, scat_p_index_, scat_lat_index_, scat_lon_index_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME("output_file_formatSetAscii"),
        DESCRIPTION
        (
         "Sets the output file format to ASCII.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( output_file_format_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("output_file_formatSetBinary"),
        DESCRIPTION
        (
         "Sets the output file format to binary.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( output_file_format_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("output_file_formatSetZippedAscii"),
        DESCRIPTION
        (
         "Sets the output file format to zipped ASCII.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( output_file_format_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ParticleTypeAddAll" ),
        DESCRIPTION
        (
         "Read single scattering data and particle number densities.\n"
         "\n"
         "The WSV *pnd_field_raw* containing particle number densities for all \n"
         "hydrometeor species can be generated outside ARTS, for example by using\n"
         "PyARTS. This method needs as input a file containing filenames of \n"
         "single scattering data and a file containing the corresponding \n"
         "*pnd_field_raw*.\n"
         "\n"
         "Very important note: \n"
         "The order of the filenames for the scattering data files has to\n"
         "correspond to the order of the particle types in the file \n"
         "including the variable *pnd_field_raw*!\n" 
         "\n"
         "Keywords:\n"
         "   filename_scat_data : File containing an \n" 
         "                        ArrayOfString of filenames of single scattering data files \n"
         "                        corresponding the the particle number densities in \n"
         "                        *pnd_field_raw*.\n"
         "   filename_pnd_field : File including  the WSV *pnd_field_raw*.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( scat_data_raw_, pnd_field_raw_ ),
        INPUT( atmosphere_dim_, f_grid_, p_grid_, lat_grid_, lon_grid_, 
               cloudbox_limits_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "filename_scat_data", "filename_pnd_field" ),
        DEFAULTS( NODEF,                NODEF ),
        TYPES(    String_t,             String_t )));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ParticleTypeAdd" ),
        DESCRIPTION
        (
         "This method reads single scattering data and the corresonding\n"
         "particle number density fields. \n"
         "\n"
         "The methods reads the  specified files \n"
         "and appends the variables *scat_data_raw* and *pnd_field_raw*. \n"
         "\n"
         "Keywords:\n"
         "   filename_scat_data : Filename of single scattering data \n"
         "   filename_pnd_field : Filename of the corresponding pnd_field \n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( scat_data_raw_, pnd_field_raw_ ),
        INPUT( atmosphere_dim_, f_grid_, p_grid_, lat_grid_, lon_grid_, 
               cloudbox_limits_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "filename_scat_data", "filename_pnd_field" ),
        DEFAULTS( NODEF,                NODEF ),
        TYPES(    String_t,             String_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ParticleTypeInit" ),
        DESCRIPTION
        (
         "This method initializes variables containing data about the \n"
         "optical properties of particles (*scat_data_raw*) and about the \n"
         "particle number distribution (*pnd_field_raw*)\n"
         "\n"
         "*ParticleTypeInit* has to be executed before executing \n"
         "*ParticleTypeAdd(All)*.\n"
        ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( scat_data_raw_, pnd_field_raw_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ), 
        DEFAULTS( ),
        TYPES( ))); 

   md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_matCalc" ),
        DESCRIPTION
        (
         "This function sums up the phase matrices for all particle\n"
         "types weighted with particle number density.\n"
         "\n"
         "The output of this method is *pha_mat* (Nza, Naa, stokes_dim,\n"
         "stokes_dim). The inputs are the phase matrix for the single particle\n"
         "type *pha_mat_spt* (part_types, Nza, Naa, stokes_dim, stokes_dim)\n"
         "and the local particle  number densities for all particle types namely \n"
         "the *pnd_field* (part_types, p_grid, lat_grid, lon_grid ) for given\n"
         "*p_grid*, *lat_grid*, and *lon_grid*. The particle types required \n"
         "are specified in the control file.\n"
         ),
        AUTHORS( "Sreerekha T.R." ),
        OUTPUT( pha_mat_ ),
        INPUT( pha_mat_spt_, pnd_field_, atmosphere_dim_, scat_p_index_,
               scat_lat_index_, scat_lon_index_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ))); 

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_mat_sptFromData" ),
        DESCRIPTION
        (
         "Calculation of the phase matrix for the single particle types.\n"
         "\n"
         "This function can be used in *pha_mat_spt_agenda* as part of \n"
         "the calculation\n"
         "of the scattering integral.\n"
         "\n"
         "The interpolation of the data on the actual frequency is the first\n"
         "step in this function. \n"
         "\n"
         "Then the transformation from the database coordinate system to to \n"
         "laboratory coordinate system is done.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( pha_mat_spt_ ),
        INPUT( pha_mat_spt_, scat_data_raw_, scat_za_grid_, scat_aa_grid_, 
               scat_za_index_, scat_aa_index_, f_index_, f_grid_,
               rte_temperature_, pnd_field_, scat_p_index_, scat_lat_index_,
               scat_lon_index_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ))); 

   md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_mat_sptFromMonoData" ),
        DESCRIPTION
        (
         "Calculation of the phase matrix for the single particle types.\n"
         "\n"
         "This function is the monchromatic version of *pha_mat_sptFromData*.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( pha_mat_spt_ ),
        INPUT( pha_mat_spt_, scat_data_mono_, doit_za_grid_size_,
               scat_aa_grid_, scat_za_index_, scat_aa_index_, rte_temperature_,
               pnd_field_, scat_p_index_, scat_lat_index_, scat_lon_index_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ))); 

   md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_mat_sptFromDataDOITOpt" ),
        DESCRIPTION
        (
         "Calculation of the phase matrix for the single particle types.\n"
         "\n"
         "In this function the phase matrix is extracted from \n"
         "*pha_mat_sptDOITOpt*. It can be used in the agenda\n"
         "*pha_mat_spt_agenda*. This method must be used in \n "
         "conbination with *ScatteringDataPrepareDOITOpt*. \n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( pha_mat_spt_ ),
        INPUT( pha_mat_spt_, pha_mat_sptDOITOpt_, scat_data_mono_, 
               doit_za_grid_size_,
               scat_aa_grid_, 
               scat_za_index_, scat_aa_index_, rte_temperature_,
               pnd_field_, scat_p_index_, scat_lat_index_, scat_lon_index_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS(),
        DEFAULTS( ),
        TYPES( ))); 

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pnd_fieldCalc" ),
        DESCRIPTION
        ("Interpolate the particle number density fields.\n"
         "\n"
         "This methods interpolates the particle number density field\n"
         "from the raw data *pnd_field_raw* to pnd_field* which is definded \n"
         "on sub-grids of *p_grid*, *lat_grid*, *lon_grid*, exactly on the \n"
         "part of the atmosphere where the cloudbox is defined. \n"
         "\n"
         "The method takes as input the *pnd_field_raw* \n"
         "which contains the particle number density for each\n"
         "particle type. \n"
         ),
        AUTHORS( "Sreerekha T.R.", "Claudia Emde" ),
        OUTPUT( pnd_field_ ),
        INPUT( p_grid_, lat_grid_, lon_grid_, pnd_field_raw_, atmosphere_dim_,
               cloudbox_limits_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( ))); 

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ppathCalc" ),
        DESCRIPTION
        (
         "Main function for calculation of propagation paths.\n"
         "\n"
         "There exists only one function to calculate total propagation\n"
         "paths and this is that function. The function is normally not\n"
         "visible in the control file, it is called from inside *RteCalc*.\n"
         "A reason to call this function directly would be to plot a\n"
         "propgation path.\n"
         "\n"
         "The definition of a propgation path cannot be accomodated here.\n"
         "For more information read the chapter on propagation paths in the\n"
         "ARTS user guide and read the  on-line information for\n"
         "*ppath_step_agenda* (type \"arts -d ppath_step_agenda\").\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ppath_ ),
        INPUT( ppath_step_agenda_, atmosphere_dim_, p_grid_, lat_grid_, 
               lon_grid_, z_field_, r_geoid_, z_surface_, 
               cloudbox_on_, cloudbox_limits_, rte_pos_, rte_los_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));


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
         "points and points of surface intersections. The WSV *ppath_lmax* \n"
         "gives the option to include additional points to ensure that the\n"
         "distance along the path between the points does not exceed the \n"
         "selected maximum length. No additional points are included if\n"
         "*ppath_lmax* is set to <= 0.\n"
         "\n"
         "As functions of this kind should very seldom be called directly,\n"
         "and that the functions can be called a high number of times, these\n"
         "functions do not perform any checks of the input that give\n" 
         "detailed error messages, but asserts are performed (if turned on).\n"
         "\n"
         "For further information, type see the on-line information for\n"
         "*ppath_step_agenda* (type \"arts -d ppath_step_agenda\").\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ppath_step_ ),
        INPUT( ppath_step_, atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, 
               z_field_, r_geoid_, z_surface_, ppath_lmax_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ppath_stepRefractionEuler" ),
        DESCRIPTION
        (
         "Calculates a propagation path step, considering refraction by a\n"
         "straightforward Euler approach.\n"
         "\n"
         "Refraction is taken into account by probably the simplest approach\n"
         "possible. The path is treated to consist of piece-wise geometric \n"
         "steps. A geometric path step is calculated from each point by \n"
         "using the local line-of-sight. Except for 1D zenith angles, the\n"
         "path quantities are propagated by solving the differential \n"
         "equations by the Euler method. Snell's law for spherical symmetry\n"
         "is used for 1D to update the zenith angles. \n"
         "\n"
         "See further the on-line information for *ppath_stepGeometric*\n"
         "(type \"arts -d ppath_stepGeometric\") and the user guide for more\n"
         "details on the algorithms used.\n"
         "\n"
         "The maximum length of each ray tracing step is given by the WSV\n"
         "*ppath_lraytrace*. The length will never exceed the \n" 
         "given maximum value, but can be smaller. The ray tracing steps are\n"
         "only used to determine the path. Points to describe the path for \n" 
         "*RteCalc* are included as for *ppath_stepGeometric*, this\n"
         "including the functionality of *ppath_lmax*.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ppath_step_, rte_pressure_, rte_temperature_, rte_vmr_list_, 
                refr_index_ ),
        INPUT( refr_index_agenda_, ppath_step_, atmosphere_dim_, p_grid_, 
               lat_grid_, lon_grid_, z_field_, t_field_, vmr_field_, r_geoid_,
               z_surface_, ppath_lmax_, ppath_lraytrace_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("Print"),
        DESCRIPTION
        (
         "Prints a variable on the screen.\n"
         "\n"
         "Keywords:\n"
         "   level : Output level to use. \n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( "level" ),
        DEFAULTS( "1" ),
        TYPES(    Index_t ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("PrintWorkspace"),
        DESCRIPTION
        (
         "Prints a list of initialized workspace variables.\n"
         "\n"
         "Keywords:\n"
         "   level : Output level to use. \n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "only_allocated", "level" ),
        DEFAULTS( "1",              "1" ),
        TYPES(    Index_t,          Index_t ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  ),
        PASSWORKSPACE( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("p_gridFromGasAbsLookup"),
        DESCRIPTION
        (
         "Sets *p_grid* to the frequency grid of *abs_lookup*.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( p_grid_  ),
        INPUT(  abs_lookup_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ReadXML"),
        DESCRIPTION
        (
         "Reads a workspace variable from an XML file.\n"
         "\n"
         "This is a supergeneric method. It can read variables of any group.\n"
         "\n"
         "If the filename is omitted, the variable is read\n"
         "from <basename>.<variable_name>.xml.\n"
         "\n"
         "Usage example:\n"
         "\n"
         "ReadXML(f_grid){\"frequencies.xml\"}\n"
         "Will read the frequency grid *f_grid* from the specified file.\n"
         "\n"
         "Supergeneric output:\n"
         "   Any_     : The variable to read.\n"
         "\n"
         "Keywords:\n"
         "   filename : Name of the input file.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Any_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        DEFAULTS( NODEF ),
        TYPES(    String_t   ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  ),
        PASSWORKSPACE(  false ),
        PASSWSVNAMES(   true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("refr_indexFieldAndGradients"),
        DESCRIPTION
        (
         "Calculates the field and gradients of the refractive index.\n"
         "\n"
         "This function calculates the refractive index and its gradients\n"
         "for a rectangular grid. \n"
         "\n"
         "Calculations are performed for all combinations of the given \n"
         "vectors, where the first vector shall contain pressure values, the\n"
         "second latitude values, and the last longitude values. For \n"
         "dimensions not used, the corresponding position vector is ignored.\n"
         "\n"
         "The calculated values form a Tensor4, with size:\n"
         "   [atmosphere_dim+1, np, nlat, nlon] \n"
         "where np is the number of pressures given etc. The book of the\n"
         "tensor with the following index holds:\n"
         "   0: the refractive index \n"
         "   1: radial gradient of the refractive index \n"
         "   2: latitude gradient of the refractive index \n"
         "   3: longitude gradient of the refractive index \n"
         "\n"
         "To calculate these quantities for the atmsopheric mesh, execute:\n"
         "   RefrIndexFieldAndGradients(tensor4_1,p_grid,lat_grid,lon_grid)\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( refr_index_, rte_pressure_, rte_temperature_, rte_vmr_list_ ),
        INPUT( refr_index_agenda_, atmosphere_dim_, p_grid_, lat_grid_, 
               lon_grid_, r_geoid_, z_field_, t_field_, vmr_field_ ),
        GOUTPUT( Tensor4_ ),
        GINPUT( Vector_, Vector_, Vector_  ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("refr_indexIR"),
        DESCRIPTION
        (
         "Calculates the IR refractive index due to gases in the\n"
         "Earth's atmosphere. \n"
         "\n"
         "Only refractivity of dry air is considered. The formula used is\n"
         "contributed by Michael Hoefner,bForschungszentrum Karlsruhe.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( refr_index_ ),
        INPUT( rte_pressure_, rte_temperature_, rte_vmr_list_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("refr_indexThayer"),
        DESCRIPTION
        (
         "Calculates the microwave refractive index due to gases in the\n"
         "Earth's atmosphere. \n"
         "\n"
         "The refractivity of dry air and water vapour is summed. All\n"
         "other gases are assumed ti have a negligible contribution.  \n"
         "\n"
         "The parameterisation of Thayer (Radio Science, 9, 803-807, 1974)\n"
         "is used. See also Eq. 3 and 5 of Solheim et al. (JGR, 104, \n"
         "pp. 9664). \n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( refr_index_ ),
        INPUT( rte_pressure_, rte_temperature_, rte_vmr_list_, abs_species_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("refr_indexUnit"),
        DESCRIPTION
        (
         "Sets the refractive index to 1.\n"
         "\n"
         "If this method is used, the obtained path should be identical to\n"
         "the geomtrical path.\n"
         "\n"
         "As this function does not need any input, you have to include call\n"
         "of *Ignore* for all variables expected to be used by\n"
         "*refr_index_agenda*.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( refr_index_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "RteCalc" ),
        DESCRIPTION
        (
         "Main function for calculation of spectra.\n"
         "\n"
         "The overall scheme to solve the radiative transfer equation (RTE)\n"
         "is fixed and found in this method. In short, the method calculates\n"
         "monochromatic spectra for all pencil beam directions and applies\n"
         "the sensor response on obtained radiances.\n"
         "\n"
         "The first step is to calculate the propagation path through the\n"
         "atmosphere for the considered viewing direction. The next step is\n"
         "to determine the spectrum at the starting point of the propagation\n"
         "path. The start point of the propagation path can be found at the\n"
         "top of the atmosphere, the surface, or at the boundary or inside\n"
         "the cloud box. To determine the start spectrum can involve a\n"
         "recursive call of RteCalc (for example to calculate the radiation\n"
         "reflected by the surface). After this, the vector radiative\n"
         "transfer equation is solved to the end point of the propagation\n"
         "path. Finally, the response of the sensor is applied.\n"
         "\n"
         "Analytical jacobians for gas species and temperature can be \n"
         "calcultaed along with the spectrum.\n"
         "\n"        
         "See further the user guide.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( y_, jacobian_  ),
        INPUT( ppath_step_agenda_, rte_agenda_, iy_space_agenda_,
               surface_prop_agenda_, iy_cloudbox_agenda_,
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_, 
               t_field_, vmr_field_, abs_species_, r_geoid_, z_surface_, 
               cloudbox_on_, cloudbox_limits_, sensor_response_, sensor_pos_, 
               sensor_los_, f_grid_, stokes_dim_, 
               antenna_dim_, mblock_za_grid_, mblock_aa_grid_, 
               jacobian_, jacobian_quantities_, jacobian_indices_,
               y_unit_, jacobian_unit_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "RteCalcMC" ),
        DESCRIPTION
        (
         "As *RteCalc* but using *MCGeneral* for doing monochromatic pencil\n"
         "beam calculations.\n"
         "\n"
         "This functions allows Monte Carlo (MC) calculations for sets of \n"
         "frequencies and sensor pos/los in a single run. Sensor responses\n"
         "can be included in the standard manner (as in *RteCalc*).\n"
         "\n"
         "MC unit is set as for *MCGeneral*.No antenna pattern is included.\n"
         "\n"
         "This function does not apply the MC approach when it comes\n"
         "to sensor properties. These properties are not considered when\n"
         "tracking photons, which is done in *MCGeneral* (but only for the\n"
         "antenna pattern).\n"
         "\n"
         "The MC calculation errors are all assumed be uncorrelated and each\n"
         "have a normal distribution. These properties are of relevance when\n"
         "weighting the errors with the sensor repsonse matrix. The seed is\n"
         "reset for each call of *MCGeneral* to obtain uncorrelated errors.\n"
         "\n"
         "Keyword arguments as for *MCGeneral*. The arguments are applied\n"
         "for each monochromatic pencil beam calculation individually.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( y_, mc_error_ ),
        INPUT( iy_space_agenda_, surface_prop_agenda_, opt_prop_gas_agenda_,
               abs_scalar_gas_agenda_, atmosphere_dim_,
               p_grid_, lat_grid_, lon_grid_, z_field_, 
               t_field_, vmr_field_, r_geoid_, z_surface_, 
               cloudbox_on_, cloudbox_limits_, pnd_field_, scat_data_raw_,
               sensor_response_, sensor_pos_, 
               sensor_los_, f_grid_, stokes_dim_, 
               antenna_dim_, mblock_za_grid_, mblock_aa_grid_, y_unit_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "std_err", "max_time", "max_iter", "z_field_is_1D" ),
        DEFAULTS( NODEF,     NODEF,      NODEF,      NODEF           ),
        TYPES(    Numeric_t, Index_t,    Index_t,    Index_t         )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "RteCalcNoJacobian" ),
        DESCRIPTION
        (
         "As *RteCalc* but throughout ignores jacobians.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( y_ ),
        INPUT( ppath_step_agenda_, rte_agenda_, iy_space_agenda_,
               surface_prop_agenda_, iy_cloudbox_agenda_,
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_, 
               t_field_, vmr_field_, r_geoid_, z_surface_, 
               cloudbox_on_, cloudbox_limits_, sensor_response_, sensor_pos_, 
               sensor_los_, f_grid_, stokes_dim_, 
               antenna_dim_, mblock_za_grid_, mblock_aa_grid_, y_unit_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "RteStd" ),
        DESCRIPTION
        (
         "Standard RTE function.\n"
         "\n"
         "This function does a clearsky radiative transfer calculation for\n"
         "a given propagation path. Designed to be part of *rte_agenda*.\n"
         "\n"
         "The overall strategy is to average basic atmospheric quantities\n"
         "(such as temperature) between the end points of each step of \n"
         "the propagation path, and to calculate source term and absorption\n"
         "for these averaged values.\n" 
         "\n"
         "See further the user guide.\n"
        ),
        AUTHORS( "Claudia Emde", "Patrick Eriksson" ),
        OUTPUT( iy_, diy_dvmr_, diy_dt_ ),
        INPUT( iy_, diy_dvmr_, diy_dt_, ppath_, ppath_array_,
               ppath_array_index_, f_grid_, stokes_dim_, emission_agenda_,
               abs_scalar_gas_agenda_, rte_do_vmr_jacs_,
               rte_do_t_jacs_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "RteStdWithTransmissions" ),
        DESCRIPTION
        (
         "As *RteStd*, but also returns path transmissions.\n"
         "\n"
         "The transmission to each point of the propagation path is returned\n"
         "in *ppath_transmissions*.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( iy_, ppath_transmissions_, diy_dvmr_, diy_dt_ ),
        INPUT( iy_, diy_dvmr_, diy_dt_, ppath_, ppath_array_,
               ppath_array_index_, f_grid_, stokes_dim_,
               emission_agenda_, abs_scalar_gas_agenda_,
               rte_do_vmr_jacs_, rte_do_t_jacs_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_losSet" ),
        DESCRIPTION
        (
         "Sets *rte_los* to the given angles.\n"
         "\n"
         "The keyword argument *za* is put in as first element of *rte_los*\n"
         "and *aa* as the second element. However, when *atmosphere_dim* is\n"
         "set to 1D or 2D, the length of *rte_los* is set to 1 and only the\n"
         "given zenith angle is considered.\n"
         "\n"
         "Keywords: \n"
         "   za : Zenith angle of sensor line-of-sight.\n"
         "   aa : Azimuth angle of sensor line-of-sight.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( rte_los_ ),
        INPUT( atmosphere_dim_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "za",      "aa"      ),
        DEFAULTS( NODEF,     NODEF ),
        TYPES(    Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_posAddGeoidWGS84" ),
        DESCRIPTION
        (
         "Adds a geoid radius according to WGS-84 to a geometric altitude.\n"
         "\n"
         "This function assumes that the first element of *rte_pos* is set\n"
         "to the geometric altitude for the position of the sensor. \n"
         "The variable *rte_pos* shall contain the radius instead of the\n"
         "altitude and that can be achieved by this function. The function\n"
         "adds a geoid radius to the given altitude. The geoid radius is\n"
         "taken from the WGS-84 reference ellipsoid.\n"
         "\n"
         "For 1D, the geoid radius is set to the radius of curvature of the\n"
         "WGS-84 ellipsoid for the position and observation direction \n"
         "described with *lat_1d* and *meridian_angle_1d*.\n"
         "For 2D and 3D, the geoid radius is set to the radius of the WGS-84\n"
         "ellipsoid for the latitude value in *rte_pos*.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( rte_pos_ ),
        INPUT( rte_pos_, atmosphere_dim_, lat_1d_, meridian_angle_1d_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_posAddRgeoid" ),
        DESCRIPTION
        (
         "Adds a geoid radius by interpolating *r_geoid*.\n"
         "\n"
         "This function assumes that the first element of *rte_pos* is set\n"
         "to the geometric altitude for the position of the sensor. \n"
         "The variable *rte_pos* shall contain the radius instead of the\n"
         "altitude and that can be achieved by this function. The function\n"
         "adds a geoid radius to the given altitude. The geoid radius is\n"
         "obtained by interpolation of *r_geoid*. There is an error if the\n"
         "given position is outside the latitude and longitude grids.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( rte_pos_ ),
        INPUT( rte_pos_, atmosphere_dim_, lat_grid_, lon_grid_, r_geoid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_posSet" ),
        DESCRIPTION
        (
         "Sets *rte_pos* to the given co-ordinates.\n"
         "\n"
         "The keyword argument *r_or_z* is put in as first element of\n"
         "*rte_pos*, *lat* as the second element and *lon* as third element.\n"
         "However, the length of *rte_pos* is set to *atmosphere_dim* and\n"
         "keyword arguments for dimensions not used are ignored.\n"
         "\n"
         "The first keyword argument can either be a radius, or an altitude\n"
         "above the geoid. In the latter case, a function such as\n"
         "*rte_posAddGeoidWGS84* could be called to obtain a radius as\n"
         "first element of *rte_pos*.\n"
         "\n"
         "Keywords: \n"
         "   r_or_z : Radius or geometrical altitude of sensor position.\n"
         "   lat : Latitude of sensor position.\n"
         "   lon : Longitude of sensor position.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( rte_pos_ ),
        INPUT( atmosphere_dim_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "r_or_z",  "lat",     "lon"     ),
        DEFAULTS( NODEF,     NODEF,     NODEF ),
        TYPES(    Numeric_t, Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_posShift" ),
        DESCRIPTION
        (
         "Shifts rte_pos and rte_los, and rte_gp_XXX to the end of ppath.\n"
        ),
        AUTHORS( "Cory Davis" ),
        OUTPUT( rte_pos_, rte_los_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_ ),
        INPUT( ppath_, atmosphere_dim_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS(),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_pos_and_losFromTangentPressure" ),
        DESCRIPTION
        (
         "If you are doing limb calculations it can be useful to specify\n"
         "viewing direction and sensor position by the tangent pressure.\n"
         "This function takes tan_p as a keyword argument and sets rte_los\n"
         "and rte_pos to the apropriate position on the edge of the modelled\n"
         "atmosphere\n\n"
         "This function is a work in progress. Only 1D is currently supported\n"
        ),
        AUTHORS( "Cory Davis" ),
        OUTPUT( rte_pos_, rte_los_, ppath_ ),
        INPUT( atmosphere_dim_, p_grid_, z_field_, lat_grid_, lon_grid_,
               ppath_step_agenda_, r_geoid_, z_surface_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "tan_p" ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "r_geoidSpherical" ),
        DESCRIPTION
        (
         "Sets the geoid to be a perfect sphere.\n"
         "\n"
         "The radius of the sphere is selected by the keyword argument *r*.\n"
         "If the keyword is set to be negative, the radius is set to the\n"
         "global internal variable *EARTH_RADIUS*, defined in constants.cc.\n"
         "\n"
         "Keywords:\n"
         "   r : Radius of geoid sphere. See further above.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( r_geoid_ ),
        INPUT( atmosphere_dim_, lat_grid_, lon_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "r" ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "r_geoidWGS84" ),
        DESCRIPTION
        (
         "Sets the geoid radius to match the WGS-84 reference ellipsoid.\n"
         "\n"
         "For 1D, the geoid radius is set to the radius of curvature of the\n"
         "WGS-84 ellipsoid for the position and observation direction \n"
         "described with *lat_1d* and *meridian_angle_1d*.\n"
         "For 2D and 3D, *r_geoid* is set to the radius of the WGS-84\n"
         "ellipsoid for the crossing points of the latitude and longitude\n"
         "grids.\n"
         "\n"
         "Please note that the latitude grid must contain true latitudes\n"
         "if the function shall give correct result, and not just arbitrary\n"
         "orbit angles which is allowed for 2D cases.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( r_geoid_ ),
        INPUT( atmosphere_dim_, lat_grid_, lon_grid_, lat_1d_,
               meridian_angle_1d_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ScatteringDoit" ),
        DESCRIPTION
        (
         "This method executes *doit_mono_agenda* for each frequency \n"
         "in *f_grid*. The output is the radiation field inside the cloudbox\n"
         "(*doit_i_field*) and on the cloudbox boundary (*scat_i_p* (1D), \n"
         "*scat_i_lat* and *scat_i_lon* (3D)).\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( doit_i_field_, scat_i_p_, scat_i_lat_, scat_i_lon_,
                doit_i_field1D_spectrum_ ),
        INPUT( f_grid_, scat_i_p_, scat_i_lat_, scat_i_lon_,
               doit_mono_agenda_, doit_is_initialized_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));
 
  /* Has been found to not work. Probably caused by changes in agenda 
     functionality.
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ScatteringMonteCarlo" ),
        DESCRIPTION
        (
         "This method performs a single pencil beam monochromatic scattering\n"
         "calculation using a Monte Carlo algorithm \n"
         "\n"
         "The main output variables *iy* and *mc_error* represent the \n"
         "Stokes vector leaving the cloudbox, and the estimated error in the \n"
         "Stokes vector (at the sensor!- not the cloudbox exit) respectively.\n"
         "The keyword parameter `maxiter\' describes the number of `photons\'\n"
         "used in the simulation (more photons means smaller *mc_error*).\n"
         "std_err is the desired value of mc_error, and max_time is the maximum\n"
         "allowed number of seconds for ScatteringMonteCarlo.  ScatteringMonteCarlo\n"
         "will terminate once any of the max_iter, std_err, max_time criteria are\n"
         "met.  If negative values are given for these parameters then it is ignored.\n"
         " Negative values of rng_seed seed the random number generator \n "
         "according to system time, positive rng_seed values are taken literally.\n"
         "The incoming_lookup keyword determines if incoming radiance is obtained from\n"
         "a precalculated grid (mc_incoming) or calculated on the fly\n"
          ),
        AUTHORS( "Cory Davis" ),
        OUTPUT( ppath_, ppath_step_, 
                mc_error_, mc_iteration_count_, rte_pos_, rte_los_, iy_, 
                rte_pressure_, rte_temperature_, 
                rte_vmr_list_, ext_mat_, abs_vec_ ),
        INPUT( ppath_, rte_pos_, rte_los_, ppath_step_agenda_, atmosphere_dim_,
               p_grid_, lat_grid_, lon_grid_, z_field_, r_geoid_, z_surface_,
               cloudbox_limits_, stokes_dim_, rte_agenda_, iy_space_agenda_,
               surface_prop_agenda_, t_field_, f_grid_, opt_prop_gas_agenda_,
               abs_scalar_gas_agenda_, vmr_field_,
               scat_data_mono_, pnd_field_, mc_seed_, f_index_ , mc_incoming_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "std_err", "max_time", "max_iter", "incoming_lookup",
                  "z_field_is_1D"),
        DEFAULTS( NODEF,     NODEF,      NODEF,      NODEF,
                  NODEF ),
        TYPES(    Numeric_t, Index_t,    Index_t,    Index_t,
                  Index_t )));
  */

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "scat_data_monoCalc" ),
        DESCRIPTION
        (
         "Interpolates scat_data_raw by frequency to give scat_data_mono\n"
         ),
        AUTHORS( "Cory Davis" ),
        OUTPUT( scat_data_mono_ ),
        INPUT( scat_data_raw_ ,f_grid_, f_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "scat_data_rawCheck" ),
        DESCRIPTION
        (
         "Method for checking the consistency of the optical properties\n"
         "in the database. \n"
         "\n"
         "This function can be used to check datafiles containing data for \n"
         "randomly oriented scattering media.\n"
         "It is checked whether the data is consistent. The integral over \n"
         "the phase matrix should result the scattering cross section \n"
         "<C_sca>.\n"
         "\n"
         "The check is if:\n"
         "<C_ext> - <C_sca> = <C_abs>\n"
         "\n"
         "The result is printed on the screen.\n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( ),
        INPUT( scat_data_raw_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS(),
        DEFAULTS( ),
        TYPES( )));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensorOff_NEW" ),
        DESCRIPTION
        (
         "Sets sensor WSVs to obtain monochromatic pencil beam values.\n"
         "\n"
         "The variables are set as follows:\n"
         "   antenna_dim             : 1.\n"
         "   mblock_za_grid          : Length 1, value 0.\n"
         "   mblock_aa_grid          : Empty.\n"
         "   sensor_response         : As returned by *sensor_responseInit*.\n"
         "   sensor_response_f       : As returned by *sensor_responseInit*.\n"
         "   sensor_response_pol     : As returned by *sensor_responseInit*.\n"
         "   sensor_response_za      : As returned by *sensor_responseInit*.\n"
         "   sensor_response_aa      : As returned by *sensor_responseInit*.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( sensor_response_, sensor_response_f_NEW_, 
                sensor_response_pol_NEW_, sensor_response_za_NEW_,
                sensor_response_aa_NEW_, 
                antenna_dim_, mblock_za_grid_, mblock_aa_grid_ ),
        INPUT( atmosphere_dim_, stokes_dim_, sensor_pos_, sensor_los_,
               f_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensorOff" ),
        DESCRIPTION
        (
         "Sets sensor WSVs to obtain monochromatic pencil beam values.\n"
         "\n"
         "The variables are set as follows:\n"
         "   antenna_dim        : 1.\n"
         "   mblock_za_grid     : Length 1, value 0.\n"
         "   mblock_aa_grid     : Empty.\n"
         "   sensor_response    : As returned by *sensor_responseInit*.\n"
         "   sensor_response_f  : As returned by *sensor_responseInit*.\n"
         "   sensor_response_za : As returned by *sensor_responseInit*.\n"
         "   sensor_response_aa : As returned by *sensor_responseInit*.\n"
         "   sensor_response_pol: As returned by *sensor_responseInit*.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( sensor_response_, sensor_response_f_, sensor_response_za_,
                sensor_response_aa_, sensor_response_pol_,
                antenna_dim_, mblock_za_grid_, mblock_aa_grid_ ),
        INPUT( atmosphere_dim_, stokes_dim_, sensor_pos_, sensor_los_,
               f_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_posAddGeoidWGS84" ),
        DESCRIPTION
        (
         "Adds a geoid radius according to WGS-84 to a geometric altitude.\n"
         "\n"
         "This function assumes that the first element of *sensor_pos* is\n"
         "set to the geometric altitude for the positions of the sensor. \n"
         "The variable *sensor_pos* shall contain the radius instead of the\n"
         "altitude and that can be achieved by this function. The function\n"
         "adds a geoid radius to the given altitude. The geoid radius is\n"
         "taken from the WGS-84 reference ellipsoid.\n"
         "\n"
         "For 1D, the geoid radius is set to the radius of curvature of the\n"
         "WGS-84 ellipsoid for the position and observation direction \n"
         "described with *lat_1d* and *meridian_angle_1d*.\n"
         "For 2D and 3D, the geoid radius is set to the radius of the WGS-84\n"
         "ellipsoid for the latitude values in *sensor_pos*.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( sensor_pos_ ),
        INPUT( sensor_pos_, atmosphere_dim_, lat_1d_, meridian_angle_1d_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_posAddRgeoid" ),
        DESCRIPTION
        (
         "Adds a geoid radius by interpolating *r_geoid*.\n"
         "\n"
         "This function assumes that the first element of *rte_pos* is set\n"
         "to the geometric altitude for the position of the sensor. \n"
         "The variable *rte_pos* shall contain the radius instead of the\n"
         "altitude and that can be achieved by this function. The function\n"
         "adds a geoid radius to the given altitude. The geoid radius is\n"
         "obtained by interpolation of *r_geoid*. There is an error if the\n"
         "given position is outside the latitude and longitude grids.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( sensor_pos_ ),
        INPUT( sensor_pos_, atmosphere_dim_, lat_grid_, lon_grid_, r_geoid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responseAntenna1D"),
        DESCRIPTION
        (
         "Returns the response block matrix after it has been modified by\n"
         "a 1D antenna response.\n"
         "\n"
         "The antenna diagram patterns are given as the input variable\n"
         "*antenna_diagram*, which is an array of ArrayOfMatrix. The\n"
         "structure of this variable is that the Matrix describes the\n"
         "antenna diagram values by a relative zenith angle grid, the\n"
         "ArrayOfMatrix then contains antenna diagrams for each polarisation\n"
         "given by the rows of *sensor_pol* and at the top level, the\n"
         "*antenna_diagram* contains antenna diagrams for each viewing angle\n"
         "of the antennas/beams described by *antenna_los*.\n"
         "\n"
         "The individual antenna diagrams, described by the matrices,\n"
         "contain at least two columns where the first column describes a\n"
         "relative grid of angles and the following column(s) describe\n"
         "the antenna diagram.\n"
         "\n"
         "For each level in the antenna diagram there exist two cases,\n"
         "either only one element/column of antenna gain values are given,\n"
         "this element/column will then be used for all directions/-\n"
         "polarisations/frequencies. Or else each direction/polarisation/-\n"
         "frequency is given its individual element/column.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( sensor_response_, sensor_response_za_ ),
        INPUT( sensor_response_f_, sensor_response_pol_, mblock_za_grid_,
               antenna_dim_, antenna_diagram_, sensor_norm_, antenna_los_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responseBackend"),
        DESCRIPTION
        (
         "Returns the response block matrix after it has been modified by\n"
         "a spectrometer backend response.\n"
         "\n"
         "The channel response is given as the generic input array of\n"
         "matrices, where each element in the array represent different\n"
         "polarisations given by *sensor_pol*. The individual matrices\n"
         "describe the channel responses as function of frequency, where the\n"
         "first column describes a relative grid of frequencies and the rest\n"
         "of the columns describe the backend response.\n"
         "\n"
         "For each level, the response can be described in two ways. Either\n"
         "one single array element/matrix column is given and then used for\n"
         "each polarisation/frequency. Or a complete set of array\n"
         "elements/matrix columns covering all polarisations/frequencies are\n"
         "given and in each case a individual response will be used.\n"
         "Note that for both cases there must allways be a column in the\n"
         "matrices, the first, of a relative frequency grid.\n"
         "\n"
         "Generic Input: \n"
         "  ArrayOfMatrix : The backend channel response.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( sensor_response_, sensor_response_f_ ),
        INPUT( f_backend_, sensor_response_pol_, sensor_response_za_,
               sensor_norm_ ),
        GOUTPUT( ),
        GINPUT( ArrayOfMatrix_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responseInit_NEW"),
        DESCRIPTION
        (
         "Initialises the variables summarising the sensor response.\n"
         "\n"
         "This method sets the variables to match monochromatic pencil beam\n"
         "calculations, to be further modified by inclusion of sensor\n"
         "characteristics. If pure monochromatic pencil beam calculations\n"
         "shall be performed use *sensorOff*.\n"
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
         "\n"
         "The standard order of WSM calls for creating *sensor_response* is:\n"
         "   sensor_responseInit\n"
         "   sensor_responseAntenna1D\n"
         "   sensor_responseRotation\n"
         "   sensor_responsePolarisation\n"
         "   sensor_responseMixer\n"
         "   sensor_responseBackend\n"
         "It is not necessary to include a method for all sensor responses.\n"
         "There exist several method versions for some responses.\n"
        ),
        AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
        OUTPUT( sensor_response_, sensor_response_f_NEW_, 
                sensor_response_pol_NEW_, sensor_response_za_NEW_,
                sensor_response_aa_NEW_, 
                sensor_response_f_grid_NEW_, sensor_response_pol_grid_NEW_,
                sensor_response_za_grid_NEW_, sensor_response_aa_grid_NEW_ ),
        INPUT( f_grid_, mblock_za_grid_, mblock_aa_grid_, antenna_dim_,
               atmosphere_dim_, stokes_dim_, sensor_pos_, sensor_los_, 
               sensor_norm_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responseInit"),
        DESCRIPTION
        (
         "Initialises some sensor response variables.\n"
         "\n"
         "This method sets some variables to match monochromatic pencil beam\n"
         "calculations, to be further modified by inclusion of sensor\n"
         "characteristics. If pure monochromatic pencil beam calculations\n"
         "shall be performed use *sensorOff*.\n"
         "\n"
         "The variables are set as follows:\n"
         "   sensor_response : Identity matrix, with size matching *f_grid*,\n"
         "                     *mblock_za_grid*, *mblock_aa_grid* and \n"
         "                     *sensor_pol*.\n"
         "   sensor_response_f  : Equal to *f_grid*.\n"
         "   sensor_response_za : Equal to *mblock_za_grid*.\n"
         "   sensor_response_aa : Equal to *mblock_aa_grid*.\n"
         "   sensor_response_pol: Equal to *stokes_dim*.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( sensor_response_, sensor_response_f_, sensor_response_za_,
                sensor_response_aa_, sensor_response_pol_  ),
        INPUT( f_grid_, mblock_za_grid_, mblock_aa_grid_, antenna_dim_,
               atmosphere_dim_, stokes_dim_, sensor_norm_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responseMixer_NEW"),
        DESCRIPTION
        (
         "Returns the response block matrix after it has been modified by\n"
         "the mixer and sideband filter. The returned matrix converts RF to\n"
         "IF.\n"
        ),
        AUTHORS( "Mattias Ekstrom", "Patrick Eriksson" ),
        OUTPUT( sensor_response_, sensor_response_f_NEW_, 
                sensor_response_f_grid_NEW_ ),
        INPUT( sensor_response_, sensor_response_f_NEW_, 
               sensor_response_f_grid_NEW_,
               sensor_response_pol_grid_NEW_, sensor_response_za_grid_NEW_,
               sensor_response_aa_grid_NEW_,
               lo_NEW_, sideband_response_NEW_, sensor_norm_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responseMixer"),
        DESCRIPTION
        (
         "Returns the response block matrix after it has been modified by\n"
         "the mixer and sideband filter. The returned matrix converts RF to\n"
         "IF.\n"
         "\n"
         "The generic input matrix is a two-column matrix where the first\n"
         "column should be equal to *f_grid* and the second column desrcibes\n"
         "the sideband filter function.\n"
         "\n"
         "The local oscillator frequency is set by the keyword *lo*.\n"
         "\n"
         "Generic Input: \n"
         "       Matrix : The sideband filter response matrix.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( sensor_response_, sensor_response_f_, f_mixer_ ),
        INPUT( sensor_response_pol_, sensor_response_za_, lo_, sensor_norm_ ),
        GOUTPUT( ),
        GINPUT( Matrix_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responseMultiMixerBackend"),
        DESCRIPTION
        (
         "Returns the response block matrix after it has been modified by\n"
         "a mixer-sideband filter-spectrometer configuration, where several\n"
         "mixers are allowed.\n"
         "\n"
         "The sideband filter is represented by the generic input matrix,\n"
         "where the first column should hold frequencies and the following\n"
         "columns describe the sideband filter function for each polarisation."
         "\n\n"
         "The local oscillator frequencies is set by the WSV *lo* and the\n"
         "backend channel frequencies by *f_backend* which should both have\n"
         "the same length as the number of rows of *sensor_pol*.\n"
         "\n"
         "The channel response is given by the WSV *backend_channel_response*.\n"
         "\n"
         "Generic Input: \n"
         "         Matrix : The sideband filter response matrix.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( sensor_response_, sensor_response_f_, f_mixer_, 
                sensor_response_pol_ ),
        INPUT( sensor_response_za_, sensor_response_aa_,
               lo_, sensor_norm_, f_backend_, sensor_pol_,
               backend_channel_response_ ),
        GOUTPUT( ),
        GINPUT( Matrix_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responsePolarisation"),
        DESCRIPTION
        (
         "Adds polarisation to the response matrix.\n"
         "\n"
         "The output polarisations are given by matrix *sensor_pol*.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( sensor_response_, sensor_response_pol_ ),
        INPUT( sensor_pol_, sensor_response_za_, sensor_response_aa_,
               sensor_response_f_, stokes_dim_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responseRotation"),
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
        OUTPUT( sensor_response_ ),
        INPUT( sensor_rot_, antenna_los_, antenna_dim_, stokes_dim_,
               sensor_response_f_, sensor_response_za_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "SparseCreate" ),
        DESCRIPTION
        (
         "Creates an empty Sparse matrix.\n"
         "\n"
         "If the variable already exists, it'll be reset.\n"
         "\n"
         "Generic output: \n"
         "   Sparse: New empty Sparse matrix.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Sparse_ ),
        GINPUT(),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "StringCreate" ),
        DESCRIPTION
        (
         "Creates an empty String.\n"
         "\n"
         "If the variable already exists, it'll be reset.\n"
         "\n"
         "Generic output: \n"
         "   String: New empty String.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( String_ ),
        GINPUT(),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("StringSet"),
        DESCRIPTION
        (
         "Sets a String to the given text string.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( String_ ),
        GINPUT( ),
        KEYWORDS( "text"   ),
        DEFAULTS( NODEF ),
        TYPES(    String_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "surfaceBlackbody" ),
        DESCRIPTION
        (
         "Creates variables to mimic a blackbody surface.\n"
         "\n"
         "This method sets up *surface_los*, *surface_rmatrix* and\n"
         "*surface_emission* for *surfaceCalc*. In this case, *surface_los*\n"
         "and *surface_rmatrix* are set to be empty, and *surface_emission*\n"
         "to hold blackbody radiation for a temperature of *surface_skin_t*.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( surface_los_, surface_rmatrix_, surface_emission_ ),
        INPUT( f_grid_, stokes_dim_, surface_skin_t_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME( "surfaceEmissivityInterpolate" ),
//         DESCRIPTION
//         (
//          "Creates variables to mimic specular reflection by a surface with\n"
//          "emissivity interpolated from WSV surface_emissivity_field by lat..\n"
//          "and lon.\n"
//          "A constant emissivity is assumed as a function of frequency and\n"
//          "polarisation (vertical and horisontal reflection coefficients are\n"
//          "equal. The number of directions in *surface_los* is one.\n"
//          ),
//         AUTHORS( "Cory Davis" ),
//         OUTPUT( surface_los_, surface_rmatrix_, surface_emission_ ),
//         INPUT( f_grid_, rte_gp_lat_,rte_gp_lon_,stokes_dim_, atmosphere_dim_, rte_los_, 
//                surface_skin_t_, surface_emissivity_field_ ),
//         GOUTPUT( ),
//         GINPUT( ),
//         KEYWORDS( ),
//         DEFAULTS( ),
//         TYPES( )));


  md_data_raw.push_back
    ( MdRecord
      ( NAME( "surfaceFlat" ),
        DESCRIPTION
        (
         "Creates variables to mimic specular reflection by a surface with\n"
         "dielectric constant following an internal model.\n"
         "\n"
         "The method results in that the reflection properties differ\n"
         "between frequencies and polarizations. The properties of the\n"
         "surface medium are determined by the model for dielectric constant\n"
         "selected. Local thermodynamic equilibrium is assumed, which\n"
         "corresponds to that the reflection and emission coefficients add up\n"
         "to 1.\n"
         "\n"
         "Available dielectric models:\n"
         "\n"
         " \"water-liebe93\"\n"
         "   Treats liquid water without salt. Not valid below 10 GHz.\n"
         "   Upper frequency limit not known. Model parameters taken from\n"
         "   Atmlab function epswater93 (by C. Maetzler), which refer to\n"
         "   Liebe 93 without closer specifications.\n"
         "\n"
         "Keyword: \n"
         "   epsmodel : Name of model for dielectric constant.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( surface_los_, surface_rmatrix_, surface_emission_ ),
        INPUT( f_grid_, stokes_dim_, atmosphere_dim_, rte_los_, 
               surface_skin_t_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "epsmodel" ),
        DEFAULTS( NODEF ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "surfaceSimple" ),
        DESCRIPTION
        (
         "Creates variables to mimic specular reflection by a surface with\n"
         "a single emissivity.\n"
         "\n"
         "A constant emissivity is assumed as a function of frequency and\n"
         "polarization (vertical and horizontal reflection coefficients are\n"
         "equal. The number of directions in *surface_los* is one.\n"
         "\n"
         "Surface properties are specified by *surface_emissivity* and \n"
         "*surface_skin_t*.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( surface_los_, surface_rmatrix_, surface_emission_ ),
        INPUT( f_grid_, stokes_dim_, atmosphere_dim_, rte_los_, 
               surface_emissivity_, surface_skin_t_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor3Create" ),
        DESCRIPTION
        (
         "Creates an empty Tensor3.\n"
         "\n"
         "If the variable already exists, it'll be reset.\n"
         "\n"
         "Generic output: \n"
         "   Tensor3: New empty Tensor3.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor3_ ),
        GINPUT(),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor3ExtractFromTensor4"),
        DESCRIPTION
        (
         "Extract a Tensor3 from a Tensor4.\n"
         "\n"
         "Copies book with given Index from input Tensor4 variable to create \n"
         "output Tensor3.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Tensor3_ ),
        GINPUT(  Tensor4_, Index_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor3FillWithVector"),
        DESCRIPTION
        (
         "Forms a tensor of order 3 by repeating a vector.\n"
         "\n"
         "The direction of the vector inside the tensor is selected by\n"
         "setting the size determined by the vector length to 0. For \n"
         "example, if the keyword *ncols* is set to 0, the vector will be\n"
         "put in as rows on every page. The remaining sizes are taken from \n"
         "the keyword arguments. \n"
         "\n"
         "One, but only one, keyword argument must be 0.\n"
         "\n"
         "Generic output: \n"
         "   Tensor3 : The tensor to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied. \n"
         "Keyword: \n"
         "   npages : Number of pages in the tensor.\n"
         "   nrows  : Number of rows in the tensor.\n"
         "   ncols  : Number of columns in the tensor.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Tensor3_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "npages", "nrows", "ncols"   ),
        DEFAULTS( NODEF,    NODEF,   NODEF ),
        TYPES(    Index_t,  Index_t, Index_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor3Scale"),
        DESCRIPTION
        (
         "Scales a workspace tensor3 with the specified value. \n"
         "\n"
         "The result can either be stored in the input tensor3 or\n"
         "in a new tensor3.\n"
         "\n"
         "Generic output: \n"
         "   Tensor3 : The scaled tensor3. \n"
         "\n"
         "Generic input: \n"
         "   Tensor3 : The tensor3 to be scaled.\n"
         "\n"
         "Keywords:\n"
         "   value  : The scale factor.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Tensor3_ ),
        GINPUT( Tensor3_ ),
        KEYWORDS( "value"   ),
        DEFAULTS( NODEF ),
        TYPES( Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor3Set"),
        DESCRIPTION
        (
         "Creates a workspace tensor3 and sets all elements of the \n"
         "tensor3 to the specified value. The size is determined by \n"
         "the variables *ncols*, *nrows*, and *npages* \n"
         "\n"
         "Generic output: \n"
         "   Tensor3 : The tensor3 to be created. \n"
         "\n"
         "Keywords:\n"
         "   value  : The value of the tensor3 elements.\n"
        ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( ),
        INPUT( npages_, nrows_, ncols_ ),
        GOUTPUT( Tensor3_ ),
        GINPUT( ),
        KEYWORDS( "value"   ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor4Create" ),
        DESCRIPTION
        (
         "Creates an empty Tensor4.\n"
         "\n"
         "If the variable already exists, it'll be reset.\n"
         "\n"
         "Generic output: \n"
         "   Tensor4: New empty Tensor4.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor4_ ),
        GINPUT(),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor4ExtractFromTensor5"),
        DESCRIPTION
        (
         "Extract a Tensor4 from a Tensor5.\n"
         "\n"
         "Copies shelf with given Index from input Tensor5 variable to \n"
         "create output Tensor4.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Tensor4_ ),
        GINPUT(  Tensor5_, Index_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor4ExtractFromArrayOfTensor4"),
        DESCRIPTION
        (
         "Extract a Tensor4 from an ArrayOfTensor4 .\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Tensor4_ ),
        GINPUT(  ArrayOfTensor4_, Index_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor4Scale"),
        DESCRIPTION
        (
         "Scales a workspace tensor4 with the specified value. \n"
         "\n"
         "The result can either be stored in the input tensor4 or\n"
         "in a new tensor4.\n"
         "\n"
         "Generic output: \n"
         "   Tensor4 : The scaled tensor4. \n"
         "\n"
         "Generic input: \n"
         "   Tensor4 : The tensor4 to be scaled.\n"
         "\n"
         "Keywords:\n"
         "   value  : The scale factor.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Tensor4_ ),
        GINPUT( Tensor4_ ),
        KEYWORDS( "value"   ),
        DEFAULTS( NODEF ),
        TYPES( Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor4Set"),
        DESCRIPTION
        (
         "Creates a workspace tensor4 and sets all elements of the \n"
         "tensor4 to the specified value. The size is determined by \n"
         "the variables *ncols*, *nrows*, *npages*, and *nbooks*. \n"
         "\n"
         "Generic output: \n"
         "   Tensor4 : The tensor4 to be created. \n"
         "\n"
         "Keywords:\n"
         "   value  : The value of the tensor4 elements.\n"
        ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( ),
        INPUT( nbooks_, npages_, nrows_, ncols_ ),
        GOUTPUT( Tensor4_ ),
        GINPUT( ),
        KEYWORDS( "value"   ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor5Create" ),
        DESCRIPTION
        (
         "Creates an empty Tensor5.\n"
         "\n"
         "If the variable already exists, it'll be reset.\n"
         "\n"
         "Generic output: \n"
         "   Tensor5: New empty Tensor5.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor5_ ),
        GINPUT(),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor5Scale"),
        DESCRIPTION
        (
         "Scales a workspace tensor5 with the specified value. \n"
         "\n"
         "The result can either be stored in the input tensor5 or\n"
         "in a new tensor5.\n"
         "\n"
         "Generic output: \n"
         "   Tensor5 : The scaled tensor5. \n"
         "\n"
         "Generic input: \n"
         "   Tensor5 : The tensor5 to be scaled.\n"
         "\n"
         "Keywords:\n"
         "   value  : The scale factor.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Tensor5_ ),
        GINPUT( Tensor5_ ),
        KEYWORDS( "value"   ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor5Set"),
        DESCRIPTION
        (
         "Creates a workspace tensor5 and sets all elements of the \n"
         "tensor5 to the specified value. The size is determined by the \n"
         "variables *ncols*, *nrows*, *npages*, *nbooks*, and *nshelves*. \n"
         "\n"
         "Generic output: \n"
         "   Tensor5 : The tensor5 to be created. \n"
         "\n"
         "Keywords:\n"
         "   value   : The value of the tensor5 elements.\n"
        ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( ),
        INPUT( nshelves_, nbooks_, npages_, nrows_, ncols_ ),
        GOUTPUT( Tensor5_ ),
        GINPUT( ),
        KEYWORDS( "value" ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor6Create" ),
        DESCRIPTION
        (
         "Creates an empty Tensor6.\n"
         "\n"
         "If the variable already exists, it'll be reset.\n"
         "\n"
         "Generic output: \n"
         "   Tensor6: New empty Tensor6.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor6_ ),
        GINPUT(),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor6Scale"),
        DESCRIPTION
        (
         "Scales a workspace tensor6 with the specified value. \n"
         "\n"
         "The result can either be stored in the input tensor6 or\n"
         "in a new tensor6.\n"
         "\n"
         "Generic output: \n"
         "   Tensor6 : The scaled tensor6. \n"
         "\n"
         "Generic input: \n"
         "   Tensor6 : The tensor6 to be scaled.\n"
         "\n"
         "Keywords:\n"
         "   value  : The scale factor.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Tensor6_ ),
        GINPUT( Tensor6_ ),
        KEYWORDS( "value"   ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor6Set"),
        DESCRIPTION
        (
         "Creates a workspace tensor6 and sets all elements of the \n"
         "tensor6 to the specified value. The size is determined by the \n"
         "variables *ncols*, *nrows*, *npages*, *nbooks*, *nshelves*, \n"
         "and *nvitrines*. \n"
         "\n"
         "Generic output: \n"
         "   Tensor6 : The tensor6 to be created. \n"
         "\n"
         "Keywords:\n"
         "   value     : The value of the tensor6 elements.\n"
        ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( ),
        INPUT( nvitrines_, nshelves_, nbooks_, npages_, nrows_, ncols_ ),
        GOUTPUT( Tensor6_ ),
        GINPUT( ),
        KEYWORDS( "value" ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor6ToPlanckBT" ),
        DESCRIPTION
        (
         "Converts a Tensor6 of radiances to brightness temperatures by \n"
         "inverting the Planck function. \n"
         "\n"
         "Generic output: \n"
         "   Tensor6 : A Tensor6 with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Tenosr6 : A Tensor6 with radiance values. \n"
         ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( ),
        INPUT(f_index_, f_grid_),
        GOUTPUT( Tensor6_ ),
        GINPUT( Tensor6_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor7Create" ),
        DESCRIPTION
        (
         "Creates an empty Tensor7.\n"
         "\n"
         "If the variable already exists, it'll be reset.\n"
         "\n"
         "Generic output: \n"
         "   Tensor7: New empty Tensor7.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor7_ ),
        GINPUT(),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor7Scale"),
        DESCRIPTION
        (
         "Scales a workspace tensor7 with the specified value. \n"
         "\n"
         "The result can either be stored in the input tensor7 or\n"
         "in a new tensor7.\n"
         "\n"
         "Generic output: \n"
         "   Tensor7 : The scaled tensor7. \n"
         "\n"
         "Generic input: \n"
         "   Tensor7 : The tensor7 to be scaled.\n"
         "\n"
         "Keywords:\n"
         "   value  : The scale factor.\n"
        ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Tensor7_ ),
        GINPUT( Tensor7_ ),
        KEYWORDS( "value" ),
        DEFAULTS( NODEF ),
        TYPES( Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor7Set"),
        DESCRIPTION
        (
         "Creates a workspace tensor7 and sets all elements of the \n"
         "tensor7 to the specified value. The size is determined by the \n"
         "variables *ncols*, *nrows*, *npages*, *nbooks*, *nshelves*, \n"
         "*nvitrines*, and *nlibraries*. \n"
         "\n"
         "Generic output: \n"
         "   Tensor7 : The tensor7 to be created. \n"
         "\n"
         "Keywords:\n"
         "   value      : The value of the tensor7 elements.\n"
        ),
        AUTHORS( "Claudia Emde" ),
        OUTPUT( ),
        INPUT( nlibraries_, nvitrines_, nshelves_, nbooks_, npages_, nrows_,
               ncols_ ),
        GOUTPUT( Tensor7_ ),
        GINPUT( ),
        KEYWORDS( "value" ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Test"),
        DESCRIPTION
        (
         "A dummy method that can be used for test purposes.\n"
         "\n"
         "This method can be used by ARTS developers to quickly test stuff.\n"
         "The implementation is in file m_io.cc. This just saves you the \n"
         "trouble of adding a dummy method everytime you want to try \n"
         "something out quickly.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("timerStart"),
        DESCRIPTION
        (
         "Initializes the CPU timer."
         "\n"
         "Use *timerStop* to output the consumed cpu time\n"
         "since *timerStart*.\n"
         "\n"
         "Usage example:\n"
         "\n"
         "timerStart()\n"
         "ReadXML(f_grid){\"frequencies.xml\"}\n"
         "timerStop()\n"
         "Prints the CPU time spent for reading the XML file\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( timer_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));


  md_data_raw.push_back
    ( MdRecord
      ( NAME("timerStop"),
        DESCRIPTION
        (
         "Stops the CPU timer."
         "\n"
         "Use *timerStop* to output the consumed cpu time\n"
         "since *timerStart*. See *timerStart* for example\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( ),
        INPUT( timer_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorAddScalar"),
        DESCRIPTION
        (
         "Adds a scalar to all elements of a vector. \n"
         "\n"
         "The result can either be stored in the same or another vector. \n"
         "\n"
         "Generic output: \n"
         "   Vector : Return vector. \n"
         "\n"
         "Generic input: \n"
         "   Vector : Original vector. \n"
         "\n"
         "Keywords:\n"
         "   value : The value to be added to the vector.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "value" ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorCreate" ),
        DESCRIPTION
        (
         "Creates an empty Vector.\n"
         "\n"
         "If the variable already exists, it'll be reset.\n"
         "\n"
         "Generic output: \n"
         "   Vector: New empty Vector.\n"
        ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS(),
        DEFAULTS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorExtractFromMatrix"),
        DESCRIPTION
        (
         "Extract a Vector from a Matrix.\n"
         "\n"
         "Copies row or column with given Index from input Matrix variable\n"
         "to create output Vector.\n"
         "\n"
         "Keywords:\n"
         "   direction : Must be either *row* or *column*.\n"
        ),
        AUTHORS( "Patrick Eriksson, Oliver Lemke, Stefan Buehler" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Vector_ ),
        GINPUT(  Matrix_, Index_ ),
        KEYWORDS( "direction" ),
        DEFAULTS( NODEF ),
        TYPES(    String_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorInsertGridPoints"),
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
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Vector_ ),
        GINPUT(  Vector_, Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorLinSpace"),
        DESCRIPTION
        (
         "Creates a vector with linear spacing.\n"
         "\n"
         "The first element equals always the start value, and the spacing\n"
         "equals always the step value, but the last value can deviate from\n"
         "the stop value. The keyword step can be both positive and\n"
         "negative.\n"
         "\n"
         "The vector is [start, start+step, start+2*step, ...]\n "
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Keywords:\n"
         "   start : The start value. \n"
         "    stop : The maximum value of the end value. \n"
         "    step : The spacing of the vector.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Vector_ ),
        GINPUT( ),
        KEYWORDS( "start",   "stop",    "step"    ),
        DEFAULTS( NODEF,     NODEF,     NODEF ),
        TYPES(    Numeric_t, Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorLogSpace"),
        DESCRIPTION
        (
         "Creates a vector with logarithmic spacing.\n"
         "\n"
         "The first element equals always the start value, and the spacing\n"
         "equals always the step value, but note that the last value can  \n"
         "deviate from the stop value. The keyword step can be both positive\n"
         "and negative.\n"
         "\n"
         "Note, that although start has to be given in direct coordinates,\n"
         "step has to be given in log coordinates.\n"
         "\n"
         "Explicitly, the vector is:\n"
         " exp([ln(start), ln(start)+step, ln(start)+2*step, ...])\n "
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Keywords:\n"
         "   start : The start value. (Direct coordinates!)\n"
         "    stop : The maximum value of the end value. (Direct coordinates!)\n"
         "    step : The spacing of the vector. (Log coordinates!)\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Vector_ ),
        GINPUT( ),
        KEYWORDS( "start",   "stop",    "step"    ),
        DEFAULTS( NODEF,     NODEF,     NODEF ),
        TYPES(    Numeric_t, Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorMatrixMultiply"),
        DESCRIPTION
        (
         "Multiply a Vector with a Matrix and store the result in another\n"
         "Vector.\n"
         "\n"
         "This just computes the normal Matrix-Vector product, y=M*x. It is ok\n"
         "if input and output Vector are the same. This function is handy for\n"
         "multiplying the H Matrix to spectra.\n"
         "\n"
         "Generic output:\n"
         "   Vector : The result of the multiplication (dimension m).\n"
         "\n"
         "Generic input:\n"
         "   Matrix : The Matrix to multiply (dimension mxn).\n"
         "   Vector : The original Vector (dimension n).\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Vector_ ),
        GINPUT(  Matrix_, Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorNLinSpace"),
        DESCRIPTION
        (
         "Creates a vector with length *nelem*, equally spaced between the \n"
         "given end values. \n"
         "\n"
         "The length must be larger than 1. \n"
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Keywords:\n"
         "   start : The start value. \n"
         "    stop : The end value. \n"  
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( nelem_ ),
        GOUTPUT(Vector_),
        GINPUT( ),
        KEYWORDS( "start",   "stop"    ),
        DEFAULTS( NODEF,     NODEF     ),
        TYPES(    Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorNLogSpace"),
        DESCRIPTION
        (
         "Creates a vector with length *nelem*, equally logarithmically \n"
         "spaced between the given end values. \n"
         "\n"
         "The length must be larger than 1. \n"
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Keywords:\n"
         "   start : The start value. \n"
         "    stop : The end value. \n"  
         "       n : Number of elements of the vector.\n" 
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( nelem_ ),
        GOUTPUT(Vector_),
        GINPUT( ),
        KEYWORDS( "start",   "stop"    ),
        DEFAULTS( NODEF,     NODEF     ),
        TYPES(    Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorScale"),
        DESCRIPTION
        (
         "Scales all elements of a vector with the same value. \n"
         "\n"
         "The result can either be stored in the same or another vector. \n"
         "\n"
         "Generic output: \n"
         "   Vector : Return vector. \n"
         "\n"
         "Generic input: \n"
         "   Vector : Original vector. \n"
         "\n"
         "Keywords:\n"
         "   value : The value to be multiplicated with the vector.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "value" ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorSet"),
        DESCRIPTION
        (
         "Creates a workspace vector and sets all elements of the \n"
         "vector to the specified value. The length of the vector is \n"
         "determined by the variable *nelem*. \n"
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Keywords:\n"
         "   value  : The value of the vector elements.\n" 
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( nelem_ ),
        GOUTPUT( Vector_ ),
        GINPUT( ),
        KEYWORDS( "value"   ),
        DEFAULTS( NODEF ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorSetExplicitly"),
        DESCRIPTION
        (
         "Create a vector from the given list of numbers.\n"
         "\n"
         "Generic output:\n"
         "   Vector : The vector to be created.\n"
         "\n"
         "Keywords:\n"
         "   values  : The vector elements.\n"
         "\n"
         "Usage:\n"
         "   VectorSetExplicitly(p_grid){[1000, 100, 10]}\n"
         "   Will create a p_grid vector with these three elements.\n"
        ),
        AUTHORS( "Stefan Buehler" ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Vector_ ),
        GINPUT( ),
        KEYWORDS( "values"   ),
        DEFAULTS( NODEF ),
        TYPES(    Vector_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorToPlanckBT" ),
        DESCRIPTION
        (
         "Converts a vector of radiances to brightness temperatures by \n"
         "inverting the Planck function.\n"
         "\n"
         "This function works as *VectorToRJBT*. However, this function \n"
         "is not recommended in connection with inversions, but can be used \n"
         "to display calculated spectra in a temperature scale.\n"
         "\n"
         "Generic output: \n"
         "   Vector : A vector with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A vector with radiance values.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( sensor_pos_, sensor_los_, sensor_response_f_,
               sensor_response_za_, sensor_response_aa_,
               sensor_response_pol_ ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorToRJBT" ),
        DESCRIPTION
        (
         "Converts a vector of radiances to brightness temperatures by \n"
         "the Rayleigh-Jeans approximation of the Planck function.\n"
         "\n"
         "This function performs a linear transformation of spectral \n"
         "radiances to an approximative temperature scale. The advantage \n"
         "of this linear transformation is that the obtained values can be \n"
         "used for retrievals if the weighting functions are handled \n"
         "likewise (by *MatrixToRJBT*). This is not the case if the \n"
         "radiances are converted to temparatures by the Planck function \n"
         "directly. \n"
         "\n"
         "The conversion assumes that the elements of the input vector are\n"
         "stored in standard order and correspond to *sensor_response_f*,\n"
         "*sensor_response_za*, *sensor_response_aa* and *sensor_pol*.\n"
         "The standard option shall accordingly be to perform this \n"
         "conversion directly after *RteCalc*.\n"
         "\n"
         "If *y* shall be converted from radiances to brightness \n"
         "temperatures: \n"
         "   VectorToRJBT(y,y){} \n"
         "\n"
         "Generic output: \n"
         "   Vector : A vector with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A vector with radiance values.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( sensor_pos_, sensor_los_, sensor_response_f_,
               sensor_response_za_, sensor_response_aa_,
               sensor_response_pol_ ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

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
         "The tangent altitudes are given as the altitude above the geoid.\n"
         "\n"
         "Generic output: \n"
         "   Vector : A vector with zenith angles. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A vector with true tangent altitudes\n"
        ),
        AUTHORS( "Patrick Eriksson", "Mattias Ekstrom" ),
        OUTPUT( refr_index_, rte_pressure_, rte_temperature_, rte_vmr_list_ ),
        INPUT( refr_index_agenda_, sensor_pos_, p_grid_, t_field_, z_field_,
                           vmr_field_, r_geoid_, atmosphere_dim_ ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

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
         "only for 1D, where the geoid radius is taken from *r_geoid*. The\n"
         "zenith angles are always set to be positive. The tangent altitudes\n"
         "are given as the altitude above the geoid.\n"
         "\n"
         "Generic output: \n"
         "   Vector : A vector with zenith angles. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A vector with geometric tangent altitudes\n"
        ),
        AUTHORS( "Patrick Eriksson", "Mattias Ekstrom" ),
        OUTPUT( ),
        INPUT( sensor_pos_, r_geoid_, atmosphere_dim_ ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("WriteXML"),
        DESCRIPTION
        (
         "Writes a workspace variable to an XML file.\n"
         "\n"
         "This is a supergeneric method. It can write variables of any group.\n"
         "\n"
         "If the filename is omitted, the variable is written\n"
         "to <basename>.<variable_name>.xml.\n"
         "\n"
         "Usage example:\n"
         "\n"
         "WriteXML(f_grid){\"\"}\n"
         "Will write the frequency grid *f_grid* to the default file.\n"
         "\n"
         "Supergeneric input:\n"
         "   Any_     : The variable to write.\n"
         "\n"
         "Keywords:\n"
         "   filename : Name of the output file.\n"
         ),
        AUTHORS( "Oliver Lemke" ),
        OUTPUT( ),
        INPUT( output_file_format_ ),
        GOUTPUT( ),
        GINPUT(  Any_ ),
        KEYWORDS( "filename" ),
        DEFAULTS( "" ),
        TYPES(    String_t   ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  ),
        PASSWORKSPACE(  false ),
        PASSWSVNAMES(   true  )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("WriteXMLIndexed"),
        DESCRIPTION
        (
         "As *WriteXML*, but creates indexed file names.\n"
         "\n"
         "The variable is written to a file with name:\n"
         "   <filename>.<file_index>.xml.\n"
         "where <file_index> is the value of *file_index*. This:\n"
         "  IndexSet(file_index){0} \n"
         "  IndexStep(file_index){} \n"
         "  WriteXML(ppath){\"ppath\"}\n"
         "will create the file ppath.1.xml.\n"
         "\n"
         "This means that *filename* shall here not include the .xml\n"
         "extension. Omitting filename works as for *WriteXML*.\n"
         "\n"
         "Supergeneric input:\n"
         "   Any_     : The variable to write.\n"
         "\n"
         "Keywords:\n"
         "   filename : Name of the output file.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( ),
        INPUT( output_file_format_, file_index_ ),
        GOUTPUT( ),
        GINPUT(  Any_ ),
        KEYWORDS( "filename" ),
        DEFAULTS( NODEF ),
        TYPES(    String_t   ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  ),
        PASSWORKSPACE(  false ),
        PASSWSVNAMES(   true  )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ybatchCalc" ),
        DESCRIPTION
        (
         "Performs batch calculations.\n"
         "\n"
         "The method performs the following:\n"
         "   1. Sets *ybatch_index* = 0.\n"
         "   2. Performs a-d until *ybatch_index* = *ybatch_n*.\n"
         "    a. Executes *ybatch_calc_agenda*.\n"
         "    b. If *ybatch_index* = 0, resizes *ybatch* based\n"
         "       on *ybatch_n* and length of *y*.\n"
         "    c. Copies *y* to column *ybatch_index* of *ybatch*.\n"
         "    d. Adds 1 to *ybatch_index*.\n"
         "\n"
         "This means that, beside the *ybatch_calc_agenda*, the WSV\n"
         "*ybatch_n* must be set before calling this method. Further,\n"
         "*ybatch_calc_agenda* is expected to produce a spectrum and should\n"
         "accordingly include a call of *RteCalc* (or a similar method). \n"
         "\n"
         "An agenda that calculates spectra for different temperature profiles\n"
         "could look like this:\n"
         "\n"
         "   AgendaSet(ybatch_calc_agenda){\n"
         "      Tensor3ExtractFromTensor4(t_field,tensor4_1,ybatch_index){}\n"
         "      RteCalc{}\n"
         "   }\n"
         "\n"
         "See the user guide for further practical examples.\n"
         "\n"
         "Keywords:\n"
         "   robust : a flag with value 1 or 0. If set to one, the batch\n"
         "            calculation will continue, even if individual jobs\n"
         "            fail. In that case, a warning message is written to \n"
         "            screen and file (out1 output stream), and ybatch for the\n"
         "            failed job is set to -1. The robust behavior does only work\n"
         "            properly if you have compiled the program without OpenMP!\n"
         "            (Use the configure option \"--disable-vectorize\".)\n"
         ),
        AUTHORS( "Patrick Eriksson, Stefan Buehler" ),
        OUTPUT( ybatch_ ),
        INPUT( ybatch_n_, ybatch_calc_agenda_ ), 
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "robust" ),
        DEFAULTS( "0"),
        TYPES(    Index_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ybatchMetProfiles" ),
        DESCRIPTION
        (
         "This method is used for simulating ARTS for metoffice model fields"
         "\n"
         "This method reads in *met_amsu_data* which contains the\n"
         "lat-lon of the metoffice profile files as a Matrix. It then \n"
         "loops over the number of profiles and corresponding to each \n"
         "longitude create the appropriate profile basename. Then, \n"
         "Corresponding to each basename we have temperature field, altitude\n"
         "field, humidity field and particle number density field.  The\n"
         "temperature field and altitude field are stored in the same dimensions\n"
         "as *t_field_raw* and *z_field_raw*.  The oxygen and nitrogen VMRs are\n"
         "set to constant values of 0.209 and 0.782, respectively and are used\n"
         "along with humidity field to generate *vmr_field_raw*.  \n"
         "\n"
         "The three fields *t_field_raw*, *z_field_raw*, and *vmr_field_raw* are\n"
         "given as input to *met_profile_calc_agenda* which is called in this\n"
         "method.  See documentation of WSM *met_profile_calc_agenda* for more\n"
         "information on this agenda.  \n"
         "\n"
         "The method also converts satellite zenith angle to appropriate \n"
         "*sensor_los*.  It also sets the *p_grid* and *cloudbox_limits* \n"
         "from the profiles inside the function\n"
         ),
        AUTHORS( "Sreerekha T.R." ),
        OUTPUT( ybatch_ ),
        INPUT( abs_species_, met_profile_calc_agenda_, f_grid_, met_amsu_data_,
               sensor_pos_, r_geoid_, lat_grid_, lon_grid_, atmosphere_dim_,
               scat_data_raw_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "nelem_p_grid", "met_profile_path", "met_profile_pnd_path" ),
        DEFAULTS( NODEF,          NODEF,              NODEF ),
        TYPES(    Index_t,        String_t,           String_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ybatchMetProfilesClear" ),
        DESCRIPTION
        (
         "This method is used for simulating ARTS for metoffice model fields\n"
         "for clear sky conditions.\n"
         "\n"
         "This method reads in *met_amsu_data* which contains the\n"
         "lat-lon of the metoffice profile files as a Matrix. It then \n"
         "loops over the number of profiles and corresponding to each \n"
         "longitude create the appropriate profile basename. Then, \n"
         "Corresponding to each basename we have temperature field, altitude\n"
         "field, humidity field and particle number density field.  The\n"
         "temperature field and altitude field are stored in the same dimensions\n"
         "as *t_field_raw* and *z_field_raw*.  The oxygen and nitrogen VMRs are\n"
         "set to constant values of 0.209 and 0.782, respectively and are used\n"
         "along with humidity field to generate *vmr_field_raw*.  \n"
         "\n"
         "The three fields *t_field_raw*, *z_field_raw*, and *vmr_field_raw* are\n"
         "given as input to *met_profile_calc_agenda* which is called in this\n"
         "method.  See documentation of WSM *met_profile_calc_agenda* for more\n"
         "information on this agenda.  \n"
         "\n"
         "The method also converts satellite zenith angle to appropriate \n"
         "*sensor_los*.  It also sets the *p_grid* and *cloudbox_limits* \n"
         "from the profiles inside the function\n"
         ),
        AUTHORS( "Seerekha T.R." ),
        OUTPUT( ybatch_ ),
        INPUT( abs_species_, met_profile_calc_agenda_, 
               f_grid_, met_amsu_data_, sensor_pos_, r_geoid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "nelem_p_grid", "met_profile_path" ),
        DEFAULTS( NODEF,          NODEF ),
        TYPES(    Index_t,        String_t)));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "yUnit" ),
        DESCRIPTION
        (
         "Conversion of *y* to other spectral units.\n"
         "\n"
         "The conversion specified by *y_unit* is applied. This function can\n"
         "be used if the standard way of making the conversion inside the\n"
         "radiative transfer function does not work. The WSV *y_unit* should\n"
         "then be set to \"1\" when performing the radiative transfer\n" 
         "calculations, and be changed before calling this method.\n"
        ),
        AUTHORS( "Patrick Eriksson" ),
        OUTPUT( y_ ),
        INPUT( y_, y_unit_, sensor_pos_, sensor_los_, sensor_response_f_,
               sensor_response_za_, sensor_response_aa_,
               sensor_response_pol_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        DEFAULTS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ZaSatOccultation" ),
        DESCRIPTION
        (
         "Calculates zenith angles for satellite occultations.\n"
         "\n"
         "The zenith angles are calculated with an interval of *t_sample\n"
         "with the recieving satellite at height *z_recieve* above the geoid\n"
         "and the transmitting satellite at height *z_send*.\n"
         "The zenith angles are restricted by the two tangent altitudes\n"
         "*z_scan_low* and *z_scan_high*.\n"
        ),
        AUTHORS( "Mattias Ekstrom" ),
        OUTPUT(),
        INPUT( ppath_step_agenda_, atmosphere_dim_, p_grid_, lat_grid_,
               lon_grid_, z_field_, r_geoid_, z_surface_ ),
        GOUTPUT( Vector_ ),
        GINPUT( ),
        KEYWORDS( "z_recieve", "z_send",  "t_sample", 
                  "z_scan_low", "z_scan_high" ),
        DEFAULTS( NODEF,       NODEF,     NODEF,
                  NODEF,        NODEF ),
        TYPES(    Numeric_t,   Numeric_t, Numeric_t,
                  Numeric_t,    Numeric_t )));


  //--------------------------------------------------------------------------------
  // Zeeman Methods:
  // These were all commented out when the absorption part was
  // back-ported from ARTS-1-0 to ARTS-1-1.  

//   md_data_raw.push_back     
//     ( MdRecord
//       ( NAME("absO2ZeemanModel"),
//         DESCRIPTION
//         (
//          "Calculate oxygen absorption in the 1-1000GHz range from  the absorption"
//          " model based on P.W.Rosenkranz and H. J. Liebe (MPM).\n"
//          "Output:\n"
//          "   abs_coef    : absorption coefficients [1/m], \n"
//          "            dimension: [ f_grid, abs_p (=abs_t) ]\n"
//          "\n"
//          "Input:\n"
//          "   geomag_los          : magnetic filed strength and angle with respect to the LOS\n"
//          "   f_grid              : Frequency grid [Hz].\n"
//          "   abs_p               : List of pressures [Pa].\n"
//          "   abs_t               : List of temperatures [K].\n"
//          "                         (Must have same length as abs_p!)\n"
//          "   abs_vmr             : List of H2O volume mixing ratios [absolute number].\n"
//          "                         Must have same length as abs_p!\n"
//          "   abs_model           : String specifying the model to use.\n"
//          "                         Allowed options are:\n"
//          "   abs_user_parameters : Only used if abs_model==\"user\" or \"*Scaling\". \n"
//          "                         In that case, abs_user_parameters must have \n"
//          "                         4 or 1 element(s) respectively.\n"
//          "                         abs_model==\"user\": \n"
//          "                         1. O2 Continuum scaling factor\n"
//          "                         2. O2 line strength scaling factor\n"
//          "                         3. O2 line pressure broadening scaling factor\n"
//          "                         4. O2 line mixing scaling factor\n"
//          "                         abs_model==\"*Scaling\": \n"
//          "                         1. O2 absorption overall scaling factor\n"
//          "                         Note:\n"
//          "                         abs_user_parameters must be empty if one of the\n"
//          "                         pre-defined models is used.\n"
//         ),
//         AUTHORS( "Thomas Kuhn", "Axel von Engeln" ),
//         OUTPUT( ext_mat_zee_, abs_vec_zee_),
//         INPUT(  geomag_los_, f_grid_, 
//                 zeeman_o2_onoff_, zeeman_o2_pressure_limit_, zeeman_o2_line_,
//                 ppath_index_, rte_pressure_, 
//                 rte_temperature_,rte_vmr_list_, species_index_,
//                 abs_model_, abs_user_parameters_, stokes_dim_ ),
//         GOUTPUT( ),
//         GINPUT( ),
//         KEYWORDS( ),
//         DEFAULTS( ),
//         TYPES( )));
  
//  md_data_raw.push_back
//     ( MdRecord
//       ( NAME("abs_vecAddGasZeeman"),
//         DESCRIPTION
//         (
//          "Add zeeman absorption to the elements of absorption vector.\n"
//          "\n"
//          "The task of this method is to sum up the gas absorption of the\n"
//          "different gas species and add the result to the first element\n"
//          "of the absorption vector.\n"
//         ),
//         AUTHORS( "Sreerekha T.R." ),
//         OUTPUT( abs_vec_ ),
//         INPUT( abs_vec_, abs_vec_zee_ ),
//         GOUTPUT( ),
//         GINPUT( ),
//         KEYWORDS( ),
//         DEFAULTS( ),
//         TYPES( )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("ext_matAddGasZeeman"),
//         DESCRIPTION
//         (
//          "Add Zeeman extinction  to the elements of extinction matrix.\n"
//          " \n"
//          ),
//         AUTHORS( "Sreerekha T.R." ),
//         OUTPUT( ext_mat_ ),
//         INPUT( ext_mat_, ext_mat_zee_ ),
//         GOUTPUT( ),
//         GINPUT( ),
//         KEYWORDS( ),
//         DEFAULTS( ),
//         TYPES( )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("ZeemanO2Settings"),
//         DESCRIPTION
//         (
//          "Make the Zeeman specific settings for O2 Zeeman spectral line\n"
//          "splitting for the microwave range (1-1000 GHz)\n"
//          "\n"
//          "Keywords:\n"
//          "   ZeemanO2OnOff         : The start value. \n"
//          "   ZeemanO2PressureLimit : The end value.   \n"  
//          "   ZeemanO2Line          : Line to do with Zeeman. \n"
//         ),
//         AUTHORS( "Thomas Kuhn", "Axel von Engeln" ),
//         OUTPUT( zeeman_o2_onoff_, zeeman_o2_pressure_limit_, zeeman_o2_line_),
//         INPUT( ),
//         GOUTPUT( ),
//         GINPUT( ),
//         KEYWORDS( "ZeemanO2OnOff", "ZeemanO2PressureLimit", "ZeemanO2Line" ),
//         DEFAULTS( NODEF,           NODEF,                   NODEF ),
//         TYPES(    Index_t,         Numeric_t,               Index_t)));


//  md_data_raw.push_back
//     ( MdRecord
//       ( NAME("test_zeeman"),
//         DESCRIPTION(
//                     "\n"
//                     ),
//         AUTHORS( "Thomas Kuhn" ),
//         OUTPUT( ),
//         INPUT( opt_prop_gas_agenda_ ),
//         GOUTPUT( ),
//         GINPUT( ),
//         KEYWORDS( ),
//         DEFAULTS( ),
//         TYPES( )));


}

