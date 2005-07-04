/* Copyright (C) 2000, 2001, 2002, 2003, 2004
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
#define OUTPUT   MakeArray<Index>
#define INPUT    MakeArray<Index>
#define GOUTPUT  MakeArray<Index>
#define GINPUT   MakeArray<Index>
#define KEYWORDS MakeArray<String>
#define TYPES    MakeArray<TokValType>
#define AGENDAMETHOD(x) x
#define SUPPRESSHEADER(x) x


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
         "   z_tan_lim : Vector with start and stop tangent altitudes."
        ),
        OUTPUT(),
        INPUT( z_plat_, p_abs_, z_abs_, l_step_, refr_, refr_lfac_,
               refr_index_, r_geoid_, z_surface_ ),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "delta_t", "z_tan_lim" ),
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
         "   "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
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
  
  // New name: abs_lookupAdapt
  md_data_raw.push_back     
    ( MdRecord
      ( NAME("gas_abs_lookupAdapt"),
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
         "table. If something is not ok, a runtime error is thrown."
        ),
        OUTPUT( gas_abs_lookup_, gas_abs_lookup_is_adapted_ ),
        INPUT(  gas_abs_lookup_, gas_species_, f_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  // New name: abs_lookupInit
  md_data_raw.push_back     
    ( MdRecord
      ( NAME("gas_abs_lookupInit"),
        DESCRIPTION
        (
         "Creates an empty gas absorption lookup table.\n"
         "\n"
         "This is mainly there to help developers. For example, you can write\n"
         "the empty table to an XML file, to see the file format."
        ),
        OUTPUT( gas_abs_lookup_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "abs_scalar_gas is adjusted accordingly."
        ),
        OUTPUT( abs_scalar_gas_ ),
        INPUT(  gas_abs_lookup_, gas_abs_lookup_is_adapted_,
                f_index_, 
                rte_pressure_, rte_temperature_, rte_vmr_list_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("abs_scalar_gas_fieldCalc"),
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
         "*scalar_gas_absorption_agenda*, which needs the input variables\n"
         "*rte_pressure*, *rte_temperature*, and *rte_vmr_list*, and returns the\n"
         "output variable *abs_scalar_gas*."
        ),
        OUTPUT( abs_scalar_gas_field_,
                abs_scalar_gas_,
                rte_pressure_, rte_temperature_, rte_vmr_list_),
        INPUT(  scalar_gas_absorption_agenda_,
                f_index_,
                f_grid_,
                atmosphere_dim_,
                p_grid_, lat_grid_, lon_grid_,
                t_field_, vmr_field_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  // New name: abs_speciesAdd
  md_data_raw.push_back
    ( MdRecord
      ( NAME("gas_speciesAdd"),
        DESCRIPTION
        (
         "Adds species tag groups to the list of gas species.\n"
         "\n"
         "This WSM is similar to *gas_speciesSet*, the only difference is that\n"
         "this method appends species to an existing list of gas species instead\n"
         "of creating the whole list.\n"
         "\n"
         "See *gas_speciesSet* for details on how tags are defined and examples of\n"
         "how to input them in the control file.\n"
         "\n"
         "Keywords:\n"
         "   species : Specify one String for each tag group that you want to\n"
         "             add. Inside the String, separate the tags by commas\n"
         "             (plus optional blanks).\n"
         ),
        OUTPUT( gas_species_ ),
        INPUT(  gas_species_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "species" ),
        TYPES(    Array_String_t   )));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME("gas_speciesInit"),
        DESCRIPTION
        (
         "Sets  *gas_speciesSet* to be empty."
         ),
        OUTPUT( gas_species_ ),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));
 
  // New name: abs_speciesSet
  // This is duplicate with the 1-0 method tgsDefine. Merge!
  md_data_raw.push_back
    ( MdRecord
      ( NAME("gas_speciesSet"),
        DESCRIPTION
        (
         "Set up the list of gas species tag groups.\n"
         "\n"
         "The workspace variable *gas_species* contains several tag groups. Each\n"
         "tag group contain one or more tags. This method converts descriptions\n"
         "of tag groups given in the keyword to the internal representation of\n"
         "*gas_species*. A tag group selects spectral features which belong to\n"
         "the same species.\n"
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
         "e.g., \"H2O-nl\". This means that no lines of this species should be\n"
         "included in the general line-by-line calculation. This feature is\n"
         "useful if you want to define a tag group just for a continuum, or for\n"
         "a complete absorption model.\n"
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
         "               \"H2O-nl\" ]\n"
         "\n"
         "   The first tag group selects all O3-666 lines between 500 and\n"
         "   501 GHz plus all O3-686 lines.  \n"
         "\n"
         "   The second tag group selects all remaining O3 transitions.\n"
         "\n"
         "   The third tag group selects H2O, but will not put any lines in the\n"
         "   line list for this species. Presumably, we are using a complete\n"
         "   absorption model like MPM89 for H2O in this case."
         ),
        OUTPUT( gas_species_ ),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "species" ),
        TYPES(    Array_String_t   )));
 
  // New name: abs_speciesSet
  // This is duplicate with the 1-1 method gas_species set. Merge!
//   md_data.push_back
//     ( MdRecord
//       ( NAME("tgsDefine"),
//   	DESCRIPTION(
// 		    "Set up the list of tag groups.\n"
// 		    "\n"
// 		    "The workspace variable *tgs* contains several tag groups. Each \n"
// 		    "tag group contain one or more tags. This method converts \n"
// 		    "description of tag groups  given in the keyword to the internal \n"
// 		    "representation *tgs*. A tag group selects spectral features which \n"
// 		    "belong to the same species. \n"
// 		    "   A tag group can contain a mixture of general and special \n"
// 		    "tags.  All the continuum tags belong to the special tags and \n"
// 		    "the rest come under the general tags.\n"
// 		    "   A general tag is defined in terms of the name of the species,\n"
// 		    "isotope and a range of frequencies. Species are named after the \n"
// 		    "standard chemical names,e.g., \"O3\".  Isotopes are given by the \n"
// 		    "last digit of the atomic weight, i.e., \"O3-668\" for the \n"
// 		    "asymmetric ozone molecule including an oxygen 18 atom.  Groups\n"
// 		    "of transitions are specified by giving a lower and upper limit \n"
// 		    "of a frequency range,\"O3-666-500e9-501e9\".Moreover the symbol\n"
// 		    "'*' acts as a wild card. Furthermore, frequency range or frequency\n"
// 		    "range and isotope may be omitted.\n"
// 		    "Example for some tag groups containing only general tags:\n"
// 		    "tags = [\"O3-666-500e9-501e9, O3-686\",\"O3\"]\n"
// 		    "The first tag group consist of all O3-666 lines between 500 and\n"
// 		    "501 GHz plus all O3-686 lines.  The second tag group will contain\n"
// 		    "all remaining O3 transitions.\n"
// 		    "\n"
// 		    "Keywords:\n"
// 		    "   tags : Specify one String for each tag group that you want to create.\n"
// 		    "   Inside the String, separate the tags by comma (plus optional blanks).\n"
// 		    "   Example:\n"
// 		    "   tag = [\"O3-686\",\"H2O\"]"),
// 	OUTPUT( tgs_ ),
// 	INPUT(),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS( "tags" ),
// 	TYPES(    Array_String_t   )));
  
 md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_vecAddGas"),
        DESCRIPTION
        (
         "Add gas absorption to first element of absorption vector.\n"
         "\n"
         "The task of this method is to sum up the gas absorption of the\n"
         "different gas species and add the result to the first element of the\n"
         "absorption vector."
         ),
        OUTPUT(abs_vec_),
        INPUT(abs_vec_, abs_scalar_gas_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
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
         "are specified in the control file.  \n"
         ),
        OUTPUT(abs_vec_),
        INPUT(abs_vec_, abs_vec_spt_, pnd_field_, atmosphere_dim_,
              scat_p_index_,  scat_lat_index_, scat_lon_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "leading frequency dimension."
         ),
        OUTPUT(abs_vec_),
        INPUT(f_grid_, stokes_dim_, f_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Agenda_ ),
        KEYWORDS(),
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
         "   No keywords, but other methods can appear in the method body."
        ),
        OUTPUT(  ),
        INPUT(  ),
        GOUTPUT( Agenda_ ),
        GINPUT(),
        KEYWORDS(),
        TYPES(),
        AGENDAMETHOD( true )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AntennaSet1D"),
        DESCRIPTION
        (
         "To set the antenna dimension to be 1D.\n"
         "\n"
         "Sets *antenna_dim* to 1 and sets *mblock_aa_grid* to be empty."
        ),
        OUTPUT( antenna_dim_, mblock_aa_grid_ ),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AntennaSet2D"),
        DESCRIPTION
        (
         "To set the antenna dimension to be 2D.\n"
         "\n"
         "Sets *antenna_dim* to 2.\n"
         "\n"
         "It is only allowed to set *antenna_dim* to 2 when *atmosphere_dim*\n"
         "equals 3."
        ),
        OUTPUT( antenna_dim_ ),
        INPUT( atmosphere_dim_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "the number of polarisation given by *sensor_pol*."
        ),
        OUTPUT( antenna_diagram_ ),
        INPUT( sensor_pol_ ),
        GOUTPUT(),
        GINPUT( ArrayOfMatrix_ ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixSet"),
        DESCRIPTION
        (
         "Sets a element in an ArrayOfMatrix to a Matrix.\n"
         "\n"
         "The keyword can be used to chose which element will be set, If a\n"
         "negative number is given, the matrix will be appended to the array.\n"
         "Note that zero-based indexing is used.\n"
         "\n"
         "Generic output:\n"
         "  ArrayOfMatrix : The array to be expanded.\n"
         "\n"
         "Generic input:\n"
         "         Matrix : The matrix to append.\n"
         "Keywords:\n"
         "        element : The index to be set."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfMatrix_ ),
        GINPUT( Matrix_ ),
        KEYWORDS( "element" ),
        TYPES( Index_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfStringSet"),
        DESCRIPTION
        (
         "Sets a String array according the given text.\n"
         "The format is text = [\"String1\",\"String2\",...]"
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfString_ ),
        GINPUT(),
        KEYWORDS( "text"         ),
        TYPES(    Array_String_t )));

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
         "(*p_grid*, *lat_grid*, *lon_grid*). "
        ),
        OUTPUT(t_field_, z_field_, vmr_field_),
        INPUT(p_grid_, lat_grid_, lon_grid_, t_field_raw_, z_field_raw_, 
              vmr_field_raw_, atmosphere_dim_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmFieldsCalcExpand1D"),
        DESCRIPTION
        (
         "Interpolate 1D raw atmospheric fields to create 2D or 3D \n"
         "homogenous atmospheric fields.\n"
         "\n"
         "The method works as *AtmFieldsCalc* but accepts only raw 1D\n"
         "atmsopheres. The raw atmsophere is interpolated to *p_grid* and \n"
         "the obtained values are applied for all latitudes, and also \n"
         "longitudes for 3D, to create a homogenous atmsophere. \n"
         "\n"
         "Note that the method only deals with the atmospheric fields, and\n"
         "to create a 2D or 3D version of a 1D case, a demand is also that\n"
         "the geoid radius is set to be constant for all latitudes/longitudes."
        ),
        OUTPUT( t_field_, z_field_, vmr_field_ ),
        INPUT( p_grid_, lat_grid_, lon_grid_, t_field_raw_, z_field_raw_, 
               vmr_field_raw_, atmosphere_dim_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmosphereSet1D"),
        DESCRIPTION
        (
         "To set the atmosheric dimension to be 1D.\n"
         "\n"
         "Sets *atmosphere_dim* to 1 and gives some variables dummy values.\n"
         "\n"
         "The latitude and longitude grids are set to be empty."
        ),
        OUTPUT( atmosphere_dim_, lat_grid_, lon_grid_ ),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmosphereSet2D"),
        DESCRIPTION
        (
         "To set the atmosheric dimension to be 2D.\n"
         "\n"
         "Sets *atmosphere_dim* to 2 and gives some variables dummy values.\n"
         "\n"
         "The longitude grid is set to be empty. The variables *lat_1d*\n"
         "and *meridian_angle_1d* are given values that cause an error\n"
         "message if used."
        ),
        OUTPUT( atmosphere_dim_, lon_grid_, lat_1d_, meridian_angle_1d_),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmosphereSet3D"),
        DESCRIPTION
        (
         "To set the atmosheric dimension to be 3D.\n"
         "\n"
         "Sets *atmosphere_dim* to 3 and gives some variables dummy values.\n"
         "\n"
         "The variables *lat_1d* and *meridian_angle_1d* are given\n"
         "values that cause an error message if used."
        ),
        OUTPUT( atmosphere_dim_, lat_1d_, meridian_angle_1d_ ),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "\n"
        ),
        OUTPUT(t_field_raw_, z_field_raw_, vmr_field_raw_),
        INPUT(gas_species_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("basename"),
        TYPES(String_t)));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("BatchUpdateMatrix"),
        DESCRIPTION
        (
         "Updates a WSV of type Matrix for batch calculations.\n"
         "\n"
         "Copies page *ybatch_index* from input Tensor3 variable to create \n"
         "output Matrix."
        ),
        OUTPUT( ),
        INPUT( ybatch_index_ ),
        GOUTPUT( Matrix_ ),
        GINPUT(  Tensor3_ ),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("BatchUpdateNumeric"),
        DESCRIPTION
        (
         "Updates a WSV of type Numeric for batch calculations.\n"
         "\n"
         "Copies column *ybatch_index* from input Vector variable to create \n"
         "output Numeric."
        ),
        OUTPUT( ),
        INPUT( ybatch_index_ ),
        GOUTPUT( Numeric_ ),
        GINPUT(  Vector_ ),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("BatchUpdateTensor3"),
        DESCRIPTION
        (
         "Updates a WSV of type Tensor3 for batch calculations.\n"
         "\n"
         "Copies book *ybatch_index* from input Tensor4 variable to create \n"
         "output Tensor3."
        ),
        OUTPUT( ),
        INPUT( ybatch_index_ ),
        GOUTPUT( Tensor3_ ),
        GINPUT(  Tensor4_ ),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("BatchUpdateTensor4"),
        DESCRIPTION
        (
         "Updates a WSV of type Tensor4 for batch calculations.\n"
         "\n"
         "Copies shelf *ybatch_index* from input Tensor5 variable to create \n"
         "output Tensor4."
        ),
        OUTPUT( ),
        INPUT( ybatch_index_ ),
        GOUTPUT( Tensor4_ ),
        GINPUT(  Tensor5_ ),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("BatchUpdateVector"),
        DESCRIPTION
        (
         "Updates a WSV of type Vector for batch calculations.\n"
         "\n"
         "Copies row *ybatch_index* from input Matrix variable to create \n"
         "output Vector."
        ),
        OUTPUT( ),
        INPUT( ybatch_index_ ),
        GOUTPUT( Vector_ ),
        GINPUT(  Matrix_ ),
        KEYWORDS(),
        TYPES()));

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
         "DOIT method."
         ),
        OUTPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, iy_, ppath_, ppath_step_, 
                rte_pos_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_, rte_los_ ),
        INPUT( ppath_step_agenda_, rte_agenda_, iy_space_agenda_,
               iy_surface_agenda_, iy_cloudbox_agenda_,
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_, 
               t_field_, vmr_field_, r_geoid_, z_surface_, cloudbox_on_, 
               cloudbox_limits_, f_grid_, stokes_dim_, 
               scat_za_grid_, scat_aa_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "This method can only be used for 3D cases."
         ),
        OUTPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, iy_, ppath_, ppath_step_, 
                rte_pos_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_, rte_los_ ),
        INPUT( ppath_step_agenda_, rte_agenda_, iy_space_agenda_,
               iy_surface_agenda_, iy_cloudbox_agenda_,
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_, 
               t_field_, vmr_field_, r_geoid_, z_surface_, cloudbox_on_, 
               cloudbox_limits_, f_grid_, stokes_dim_, 
               scat_za_grid_, scat_aa_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("cloudboxOff"),
        DESCRIPTION
        (
         "Deactivates the cloud box. \n"
         "\n"
         "The function sets *cloudbox_on* to 0, *cloudbox_limits* to be an\n"
         "empty vector and *iy_cloudbox_agenda* to an empty agenda."
        ),
        OUTPUT( cloudbox_on_, cloudbox_limits_, iy_cloudbox_agenda_ ),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "\n"
         ),
        OUTPUT( pnd_field_, scat_data_raw_),
        INPUT( p_grid_, lat_grid_, lon_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "   lon2 : Upper longitude point."
        ),
        OUTPUT( cloudbox_on_, cloudbox_limits_),
        INPUT( atmosphere_dim_, p_grid_, lat_grid_, lon_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "p1", "p2", "lat1", "lat2", "lon1", "lon2" ),
        TYPES( Numeric_t, Numeric_t, Numeric_t, Numeric_t, Numeric_t, 
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
         "   lon2 : Upper longitude point."
        ),
        OUTPUT( cloudbox_on_, cloudbox_limits_),
        INPUT( atmosphere_dim_, z_field_, lat_grid_, lon_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "z1", "z2", "lat1", "lat2", "lon1", "lon2" ),
        TYPES( Numeric_t, Numeric_t, Numeric_t, Numeric_t, Numeric_t, 
               Numeric_t )));

  

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ConvertIFToRF"),
        DESCRIPTION
        (
         "Convert *sensor_response_f* from IF to RF, the function also\n"
         "unfolds the measurement spectra *y*.\n"
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
         "\n"
         "Keyword: \n"
         "   output : Which sideband(s) to output, \"lower\", \"upper\" or\n"
         "            \"double\""
         ),
         OUTPUT( sensor_response_f_, y_ ),
         INPUT( sensor_pol_, sensor_response_za_, sensor_response_aa_, lo_,
                atmosphere_dim_, sensor_pos_ ),
         GOUTPUT( ),
         GINPUT( ),
         KEYWORDS( "output" ),
         TYPES( String_t )));

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
         "   Any_ : The input variable."
         ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Any_ ),
        GINPUT(  Any_ ),
        KEYWORDS(),
        TYPES(),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

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
         "           recommended value: 10\n"
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
         "\n"
         ),
        OUTPUT( doit_za_grid_size_, scat_aa_grid_, scat_za_grid_),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("N_za_grid", "N_aa_grid", "za_grid_opt_file"),
        TYPES(Index_t, Index_t, String_t)));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "DoitCloudboxFieldPut" ),
        DESCRIPTION
        (
         "Method for the communication between cloudbox and clearsky.\n"
         "\n"
         "This method puts the scattered radiation field into the interface\n"
         "variables between the cloudbox and the clearsky, which are \n"
         "*scat_i_p*, *scat_i_lat* and *scat_i_lon*."
         "The best way to calculate spectra including the influence of\n" 
         "scattering is to set up the *scat_mono_agenda* where this method \n"
         "can be included.\n"
         "\n"
         ),
        OUTPUT( scat_i_p_, scat_i_lat_, scat_i_lon_ ),
        INPUT( doit_i_field_, f_grid_, f_index_,   p_grid_, lat_grid_, 
               lon_grid_, scat_za_grid_, scat_aa_grid_, stokes_dim_,
               atmosphere_dim_, cloudbox_limits_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "\n"
        ),
        OUTPUT(doit_conv_flag_, doit_iteration_counter_),
        INPUT(doit_conv_flag_, doit_iteration_counter_,
              doit_i_field_, doit_i_field_old_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS("epsilon"),
        TYPES(Vector_t)));

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
         "\n"
        ),
        OUTPUT(doit_conv_flag_, doit_iteration_counter_),
        INPUT(doit_conv_flag_, doit_iteration_counter_, 
              doit_i_field_, doit_i_field_old_, f_grid_, f_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS("epsilon"),
        TYPES(Vector_t)));
  
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
         "\n"
         ),
        OUTPUT(doit_conv_flag_, doit_iteration_counter_),
        INPUT(doit_conv_flag_, doit_iteration_counter_,
              doit_i_field_, doit_i_field_old_, f_grid_, f_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS("epsilon"),
        TYPES(Vector_t)));

 md_data_raw.push_back
   ( MdRecord
      ( NAME( "DoitInit" ),
        DESCRIPTION
        (
         "Initialize variables for DOIT scattering calculations. \n"
         "\n"
         "Before using the WSM *ScatteringDOIT*, please use this method \n"
         "to initialize the required WSVs. \n"
         "\n"
         ),
        OUTPUT(scat_p_index_, scat_lat_index_, scat_lon_index_, 
               scat_za_index_, scat_aa_index_, doit_scat_field_,
               doit_i_field_, doit_za_interp_),
        INPUT(stokes_dim_, atmosphere_dim_, scat_za_grid_, scat_aa_grid_,
              doit_za_grid_size_, 
              cloudbox_limits_, scat_data_raw_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "\n"
         ),
        OUTPUT(doit_i_field_),
        INPUT( doit_i_field_, doit_scat_field_agenda_, doit_rte_agenda_, 
               doit_conv_test_agenda_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         ),
        OUTPUT(doit_i_field_),
        INPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, f_grid_, 
               f_index_, p_grid_, lat_grid_, lon_grid_, 
               cloudbox_limits_, atmosphere_dim_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_i_fieldSetConst" ),
        DESCRIPTION
        (
         "Set the initial field inside the cloudbox to a constant value.\n"
         "\n"
         "The keyword *value* is specified by the user.\n"
         "\n"
         "Output of the method is the first guess field stored in the \n"
         "workspace variable *doit_i_field*.\n"
         "\n"
         ),
        OUTPUT(doit_i_field_),
        INPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, p_grid_, lat_grid_, 
               lon_grid_, 
               cloudbox_limits_, atmosphere_dim_, stokes_dim_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("value"),
        TYPES(Vector_t)));

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
         "\n" 
        ),
        OUTPUT(doit_i_field_, ppath_step_, rte_los_, rte_pos_, rte_gp_p_),
        INPUT(doit_i_field_old_, doit_scat_field_, cloudbox_limits_, 
              scalar_gas_absorption_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, pnd_field_, 
              opt_prop_part_agenda_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, z_field_, r_geoid_, t_field_,
              f_grid_, f_index_, iy_surface_agenda_, doit_za_interp_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "\n"
         ),
        OUTPUT(doit_i_field_, ppath_step_, rte_los_, rte_pos_, rte_gp_p_),
        INPUT(doit_i_field_, doit_scat_field_, cloudbox_limits_, 
              scalar_gas_absorption_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, pnd_field_,
              opt_prop_part_agenda_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, z_field_, r_geoid_, t_field_,
              f_grid_, f_index_, iy_surface_agenda_, doit_za_interp_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "\n"
         ),
        OUTPUT(doit_i_field_, rte_pressure_, rte_temperature_,
               rte_vmr_list_, scat_za_index_, ext_mat_, abs_vec_,
               scat_p_index_, ppath_step_, rte_los_, rte_pos_, rte_gp_p_),
        INPUT(doit_scat_field_, cloudbox_limits_, 
              scalar_gas_absorption_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, pnd_field_, 
              opt_prop_part_agenda_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, z_field_, r_geoid_, t_field_,
              f_grid_, f_index_, iy_surface_agenda_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 // md_data_raw.push_back
//     ( MdRecord
//       ( NAME( "doit_i_fieldUpdate3D" ),
//         DESCRIPTION
//         (
//          "RT calculation in cloudbox with fixed scattering integral. \n"
//          "\n"
//          "Update the radiation field (DOIT). This method loops\n"
//          "through the cloudbox to update the radiation field for all \n"
//          "positions and directions in the 3D cloudbox.\n"
//          "\n"
//          "Notes: This method is very inefficient, because the number of \n"
//          "iterations scales with the number of cloudbox pressure levels.\n"
//          "It is recommended to take *doit_i_fieldUpdateSeq3D*.\n"
//          "\n"
//          "Surface reflection is not yet implemented in 3D scattering \n"
//          "calculations.\n"
//          "\n " 
//         ),
//         OUTPUT(doit_i_field_, rte_pressure_, rte_temperature_,
//                rte_vmr_list_, scat_za_index_, scat_aa_index_, ext_mat_,
//                abs_vec_,
//                scat_p_index_, scat_lat_index_, scat_lon_index_,  ppath_step_),
//         INPUT(doit_i_field_old_, doit_scat_field_, cloudbox_limits_, 
//               scalar_gas_absorption_agenda_,
//               vmr_field_, spt_calc_agenda_, scat_za_grid_, scat_aa_grid_,
//               opt_prop_part_agenda_, opt_prop_gas_agenda_,
//               ppath_step_agenda_, p_grid_, lat_grid_, lon_grid_, z_field_,
//               r_geoid_, t_field_,
//               f_grid_, f_index_),
//         GOUTPUT(),
//         GINPUT(),
//         KEYWORDS(),
//         TYPES()));

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
         "\n"
         ),
        OUTPUT(doit_i_field_, ppath_step_),
        INPUT(doit_i_field_, doit_scat_field_, cloudbox_limits_, 
              scalar_gas_absorption_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, scat_aa_grid_,
              pnd_field_,
              opt_prop_part_agenda_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, lat_grid_, lon_grid_, z_field_,
              r_geoid_, t_field_,
              f_grid_, f_index_, doit_za_interp_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));
  
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
         "\n"
         ),
        OUTPUT( doit_scat_field_ ),
        INPUT( doit_scat_field_, pha_mat_spt_agenda_,
               doit_i_field_, pnd_field_, t_field_, atmosphere_dim_, 
               cloudbox_limits_, scat_za_grid_, scat_aa_grid_,  
               doit_za_grid_size_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "doit_scat_fieldCalc1D" ),
        DESCRIPTION
        (
         "This method calculates the scattering integral field.\n"
         "\n"
         "By scattering integral field we mean the field generated by integrating\n"
         "the product of phase matrix and Stokes vector over all incident \n"
         "angles. \n"
         "\n"
         "The output of this method is  *doit_scat_field*\n"
         "which is used in the radiative transfer part (*doit_i_fieldUpdateXXX*). \n"
         "\n"
         ),
        OUTPUT( doit_scat_field_, pha_mat_, pha_mat_spt_, scat_za_index_,
                scat_aa_index_),
        INPUT( pha_mat_spt_agenda_, doit_i_field_, pnd_field_, scat_za_grid_, 
               scat_aa_grid_,  
               atmosphere_dim_, cloudbox_limits_, doit_za_grid_size_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "\n"
         ),
        OUTPUT( doit_scat_field_),
        INPUT( doit_scat_field_, pha_mat_spt_agenda_,
               doit_i_field_, pnd_field_, t_field_, atmosphere_dim_, 
               cloudbox_limits_, scat_za_grid_, scat_aa_grid_,  
               doit_za_grid_size_, doit_za_interp_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("DoitScatteringDataPrepare"),
        DESCRIPTION
        (
         "Prepare single scattering data for a DOIT scattering calculation.\n"
         "\n"
         "This function can be used for scattering calcualtions using the \n" 
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
         "\n" 
          ),
        OUTPUT( pha_mat_sptDOITOpt_, scat_data_mono_),
        INPUT( doit_za_grid_size_, scat_aa_grid_, scat_data_raw_, f_grid_, 
               f_index_, atmosphere_dim_, stokes_dim_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "\n"
         ),
        OUTPUT( ),
        INPUT(doit_iteration_counter_, doit_i_field_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS("iterations"),
        TYPES(Array_Index_t )));

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
        OUTPUT(doit_za_grid_opt_),
        INPUT(doit_i_field_, scat_za_grid_, doit_za_interp_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("acc"),
        TYPES(Numeric_t)));
                                                                               
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
         "\n"
         ),
        OUTPUT( doit_za_interp_ ),
        INPUT(atmosphere_dim_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("interp_method"),
        TYPES(String_t)));
 
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
         "   Any_ : The input variable."
         ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Any_ ),
        GINPUT( Any_ ),
        KEYWORDS(),
        TYPES(),
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
         "thermodynamic equilibrium is valid."
        ),
        OUTPUT( emission_ ),
        INPUT( f_grid_, rte_temperature_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Error"),
        DESCRIPTION
        (
         "Issues an error and exits ARTS.\n"
         "\n"
         "This method can be placed in agendas that must be specified , but\n"
         "are expected not to be used for the particular case. An inclusion\n"
         "in *iy_surface_agenda* could look like:\n   "
         "Error{\"Surface interceptions of propagation path not expected.\"}\n"
         "(ignore and other dummy metho calls must still be included)\n"
         "\n"
         "Keywords: \n"
         "   msg : String describing the error."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "msg" ),
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
         "it is reached, it will terminate the program."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ext_matAddGas"),
        DESCRIPTION
        (
         "Add gas absorption to all diagonal elements of extinction matrix.\n"
         " \n"
         "The task of this method is to sum up the gas absorption of the\n"
         "different gas species and add the result to the extinction matrix."
         ),
        OUTPUT(ext_mat_),
        INPUT(ext_mat_, abs_scalar_gas_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "are specified in the control file.  \n"
         ),
        OUTPUT( ext_mat_  ),
        INPUT( ext_mat_, ext_mat_spt_, pnd_field_, atmosphere_dim_, 
               scat_p_index_, scat_lat_index_, scat_lon_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "leading frequency dimension."
         ),
        OUTPUT(ext_mat_ ),
        INPUT(f_grid_, stokes_dim_, f_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("FlagOff"),
        DESCRIPTION
        (
         "Sets an index variable that acts as an on/off flag to 0."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Index_ ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("FlagOn"),
        DESCRIPTION
        (
         "Sets an index variable that acts as an on/off flag to 1."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Index_ ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  // New name: f_gridFromAbsLookup
  md_data_raw.push_back     
    ( MdRecord
      ( NAME("f_gridFromGasAbsLookup"),
        DESCRIPTION
        (
         "Sets *f_grid* to the frequency grid of *gas_abs_lookup*.\n"
         "\n"
         "Must be called between importing/creating raw absorption table and\n"
         "call of *gas_abs_lookupAdapt*."
        ),
        OUTPUT( f_grid_  ),
        INPUT(  gas_abs_lookup_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "   max_spacing : The maximum step between grid points."
         ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT(),
        KEYWORDS( "fwhm", "tot_width", "max_spacing" ),
        TYPES( Numeric_t, Numeric_t, Numeric_t )));

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
         "   Any_ : The input variable."
         ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT(  Any_ ),
        KEYWORDS(),
        TYPES(),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

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
         "   value : A positive integer." 
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Index_ ),
        GINPUT( ),
        KEYWORDS( "value" ),
        TYPES(     Index_t   )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("IndexStep"),
        DESCRIPTION
        (
         "Performs GOUTPUT = GINPUT + 1\n"
         "\n"
         "Input and output can be same variable."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Index_ ),
        GINPUT( Index_ ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("InterpAtmFieldToRteGps"),
        DESCRIPTION
        (
         "Scalar interpolation  of atmospheric fields.\n" 
         "\n"
         "The position is specified by the combination of *rte_gp_p*, \n"
         "*rte_gp_lat* and *rte_gp_lon*.\n"
         "\n"
         "Generic output: \n"
         "   Numeric : Value obtained by interpolation. \n"
         "\n"
         "Generic input:\n"
         "   Tensor3 : Field to interpolate." 
        ),
        OUTPUT( ),
        INPUT( atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, 
               rte_gp_p_, rte_gp_lat_, rte_gp_lon_ ),
        GOUTPUT( Numeric_ ),
        GINPUT( Tensor3_ ),
        KEYWORDS(),
        TYPES()));
  
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
         "given (specified by *rte_XXX*). A linear interpolation is used in\n"
         "all dimensions.\n"
         "\n"
         "The intensity field on the cloux box boundaries is provided by\n"
         "*scat_i_p/lat/lon* and these variables are interpolated if the.\n"
         "given position is at any boundary.\n"
         "\n"
         "Interpolation of the internal field is not yet possible."
         ),
        OUTPUT( iy_ ),
        INPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, rte_gp_p_, rte_gp_lat_,
               rte_gp_lon_, rte_los_,  cloudbox_on_, cloudbox_limits_,
               atmosphere_dim_, stokes_dim_, scat_za_grid_, scat_aa_grid_, 
               f_grid_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "iyInterpPolyCloudboxField" ),
        DESCRIPTION
        (
         "As *iyInterpCloudboxField* but performs cubic interpolation.\n"
         "\n"
         "Works so far only for 1D cases, and accordingly a cubic\n"
         "interpolation along *scat_za_grid* is performed."
         ),
        OUTPUT( iy_ ),
        INPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, rte_gp_p_, rte_gp_lat_,
               rte_gp_lon_, rte_los_,  cloudbox_on_, cloudbox_limits_,
               atmosphere_dim_, stokes_dim_, scat_za_grid_, scat_aa_grid_, 
               f_grid_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianAddGas" ),
        DESCRIPTION
        (
         "Adds a gas species as a retrieval quantity to the Jacobian.\n"
         "\n"
         "This functions also adds the SpeciesTag of the given species to\n"
         "*gas_species*. This way the treated gas species only need to be\n"
         "given at one place in the control file. It also appends\n"
         "*jacobianCalcGas* with the given gas species to *jacobian_agenda*.\n"
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
         "   \"rel\" : relative unit (e.g. 1.1 means 10\% more of the gas) \n"
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
         "  dx      : Size of perturbation."
        ),
        OUTPUT( jacobian_quantities_, jacobian_agenda_, gas_species_ ),
        INPUT( jacobian_, atmosphere_dim_, p_grid_, lat_grid_, lon_grid_,
               gas_species_ ),
        GOUTPUT(),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS( "species", "method", "unit", "dx" ),
        TYPES( String_t, String_t, String_t, Numeric_t )));
         
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianAddGasAnalytical" ),
        DESCRIPTION
        (
         "As *jacobianAddGas* but assumes analytical calculations."
        ),
        OUTPUT( jacobian_quantities_, jacobian_agenda_, gas_species_ ),
        INPUT( jacobian_, atmosphere_dim_, p_grid_, lat_grid_, lon_grid_,
               gas_species_ ),
        GOUTPUT(),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS( "species", "method" ),
        TYPES( String_t, String_t )));
         
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
        OUTPUT( jacobian_quantities_, jacobian_agenda_ ),
        INPUT( jacobian_, atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, 
               pnd_field_, pnd_field_perturb_, cloudbox_limits_ ),
        GOUTPUT(),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS(),
        TYPES()));
         
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
         "  poly_order          : Order of the polynomial."
        ),
        OUTPUT( jacobian_quantities_, jacobian_agenda_ ),
        INPUT( jacobian_, sensor_pos_, sensor_time_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "dza", "poly_order" ),
        TYPES( Numeric_t, Index_t )));

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
         "  dx      : Size of perturbation."
        ),
        OUTPUT( jacobian_quantities_, jacobian_agenda_ ),
        INPUT( jacobian_, atmosphere_dim_, p_grid_, lat_grid_, lon_grid_ ),
        GOUTPUT(),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS( "hse", "method", "dx" ),
        TYPES( String_t, String_t, Numeric_t )));
  
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "jacobianCalc" ),
        DESCRIPTION
        (
         "Executes *jacobian_agenda* to calculate (parts of) *jacobian*.\n"
        ),
        OUTPUT( jacobian_ ),
        INPUT( jacobian_agenda_, jacobian_indices_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));
        
  md_data_raw.push_back
    ( MdRecord
      ( NAME("jacobianCalcGas"),
        DESCRIPTION
        (
        "Calculates the contribution of a gas species to the Jacobian.\n"
        "\n"
        "This function is added to *jacobian_agenda* by jacobianAddGas\n"
        "and should normally not be called by the user.\n"
        ),
        OUTPUT( jacobian_, vmr_field_, y_, ppath_, ppath_step_, 
                iy_, rte_pos_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_, rte_los_ ),
        INPUT( jacobian_quantities_, jacobian_indices_, gas_species_, 
               ppath_step_agenda_, 
               rte_agenda_, iy_space_agenda_, iy_surface_agenda_, 
               iy_cloudbox_agenda_, atmosphere_dim_, p_grid_, lat_grid_, 
               lon_grid_, z_field_, t_field_, 
               r_geoid_, z_surface_, cloudbox_on_, 
               cloudbox_limits_, sensor_response_, sensor_pos_, sensor_los_, 
               f_grid_, stokes_dim_, antenna_dim_, mblock_za_grid_, 
               mblock_aa_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "species" ),
        TYPES( String_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("jacobianCalcParticle"),
        DESCRIPTION
        (
        "Calculates the contribution of particle number densities to the\n"
        "Jacobian.\n"
        "\n"
        "This function is added to *jacobian_agenda* by jacobianAddParticle\n"
        "and should normally not be called by the user.\n"
        ),
        OUTPUT( jacobian_, pnd_field_, y_, ppath_, ppath_step_,
                iy_, rte_pos_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_, rte_los_ ),
        INPUT( jacobian_quantities_, jacobian_indices_, pnd_field_perturb_, 
               jacobian_particle_update_agenda_,
               ppath_step_agenda_, rte_agenda_, iy_space_agenda_, 
               iy_surface_agenda_, iy_cloudbox_agenda_, atmosphere_dim_, 
               p_grid_, lat_grid_, lon_grid_, z_field_, t_field_, vmr_field_,
               r_geoid_, z_surface_, cloudbox_on_, cloudbox_limits_, 
               sensor_response_, sensor_pos_, sensor_los_, f_grid_, 
               stokes_dim_, antenna_dim_, mblock_za_grid_, mblock_aa_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));
        
  md_data_raw.push_back
    ( MdRecord
      ( NAME("jacobianCalcPointing"),
        DESCRIPTION
        (
        "Calculates the contribution of a pointing deviation to the Jacobian.\n"
        "\n"
        "This function is added to *jacobian_agenda* by jacobianAddPointing\n"
        "and should normally not be called by the user.\n"
        ),
        OUTPUT( jacobian_, y_, ppath_, ppath_step_, iy_, rte_pos_, rte_gp_p_,
                rte_gp_lat_, rte_gp_lon_, rte_los_ ),
        INPUT( jacobian_quantities_, jacobian_indices_, 
               sensor_time_, ppath_step_agenda_, 
               rte_agenda_, iy_space_agenda_, iy_surface_agenda_, 
               iy_cloudbox_agenda_, atmosphere_dim_, p_grid_, lat_grid_, 
               lon_grid_, z_field_, t_field_, vmr_field_,
               r_geoid_, z_surface_, cloudbox_on_, 
               cloudbox_limits_, sensor_response_, sensor_pos_, sensor_los_, 
               f_grid_, stokes_dim_, antenna_dim_, mblock_za_grid_, 
               mblock_aa_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("jacobianCalcTemperature"),
        DESCRIPTION
        (
        "Calculates the contribution the temperature to the Jacobian.\n"
        "\n"
        "This function is added to *jacobian_agenda* by jacobianAddTemperature\n"
        "and should normally not be called by the user.\n"
        ),
        OUTPUT( jacobian_, t_field_, y_, ppath_, ppath_step_, iy_, 
                rte_pos_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_, rte_los_ ),
        INPUT( jacobian_quantities_, jacobian_indices_, ppath_step_agenda_, 
               rte_agenda_, 
               iy_space_agenda_, iy_surface_agenda_, iy_cloudbox_agenda_, 
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_, 
               vmr_field_, r_geoid_, z_surface_, 
               cloudbox_on_, cloudbox_limits_, 
               sensor_response_, sensor_pos_, sensor_los_, f_grid_, 
               stokes_dim_, antenna_dim_, mblock_za_grid_, mblock_aa_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         OUTPUT( jacobian_, jacobian_indices_ ),
         INPUT( jacobian_quantities_, sensor_pos_, sensor_response_ ),
         GOUTPUT(),
         GINPUT(),
         KEYWORDS(),
         TYPES()));

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
         "The Jacobian quantities are initialised to be empty."
        ),
         OUTPUT( jacobian_, jacobian_quantities_, jacobian_indices_,
                 ppath_array_do_, ppath_array_, ppath_array_index_,
                 rte_do_vmr_jacs_, rte_do_t_jacs_ ),
         INPUT(),
         GOUTPUT(),
         GINPUT(),
         KEYWORDS(),
         TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("jacobianOff"),
        DESCRIPTION
        (
         "Makes mandatory initialisation of some jacobian variables.\n"
         "\n"
         "Some jacobian WSVs must be initilised even if no such calculations\n"
         "will be performed and this is handled with this method. That is,\n"
         "this method must be called when no jacobians will be calculated"
        ),
         OUTPUT( jacobian_, jacobian_quantities_, jacobian_indices_,
                 ppath_array_do_, ppath_array_, ppath_array_index_,
                 rte_do_vmr_jacs_, rte_do_t_jacs_ ),
         INPUT(),
         GOUTPUT(),
         GINPUT(),
         KEYWORDS(),
         TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Main"),
        DESCRIPTION
        ( 
         "Run the agenda that is specified inside the curly braces. ARTS\n"
         "controlfiles must define this method. It is executed automatically\n"
         "when ARTS is run on the controlfile." 
        ),
        OUTPUT( ),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES(),
        AGENDAMETHOD( true )));

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
         "   Vector : A set of frequencies. "
        ),
        OUTPUT(),
        INPUT( stokes_dim_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
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
         "   Vector : The vector to be copied."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
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
         "   Vector : The vector to be copied into the second column."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_, Vector_ ),
        KEYWORDS( ),
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
         "   Vector : The vector to be copied into the third column."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS( ),
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
         "   Vector : The vector to be copied."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
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
         "   Vector : The vector to be copied into the second row."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_, Vector_ ),
        KEYWORDS( ),
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
         "   Vector : The vector to be copied into the third row."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS( ),
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
         "   Numeric : Blackbody temperature. "
        ),
        OUTPUT(),
        INPUT( stokes_dim_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_, Numeric_ ),
        KEYWORDS(),
        TYPES()));

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
         "   value : The value to be multiplicated with the matrix."  
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Matrix_ ),
        KEYWORDS( "value" ),
        TYPES(    Numeric_t   )));

   md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixSet"),
        DESCRIPTION
        (
         "Creates a workspace matrix and sets all elements of the \n"
         "matrix to the specified value. The size is determined by \n"
         "the variables *ncols* and *nrows*.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Keywords:\n"
         "   value : The value of the matrix elements. " 
        ),
        OUTPUT(),
        INPUT( nrows_, ncols_ ),
        GOUTPUT( Matrix_ ),
        GINPUT(),
        KEYWORDS( "value"   ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MatrixToTbByPlanck" ),
        DESCRIPTION
        (
         "Converts a matrix of radiances to brightness temperatures by \n"
         "inverting the Planck function.\n"
         "\n"
         "This function works as *MatrixToTbByRJ*. However, this function \n"
         "is not recommended in connection with inversions, but can be used \n"
         "to display calculated spectra in a temperature scale.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : A matrix with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Matrix : A matrix with radiance values."
        ),
        OUTPUT(),
        INPUT( sensor_pos_, sensor_los_, sensor_response_f_,
               sensor_response_za_, sensor_response_aa_,
               sensor_response_pol_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Matrix_ ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MatrixToTbByRJ" ),
        DESCRIPTION
        (
         "Converts a matrix of radiances to brightness temperatures by \n"
         "the Rayleigh-Jeans approximation of the Planck function.\n"
         "\n"
         "This function works as *VectorToTbByRJ*, but operates on a matrix.\n"
         "Each column of the matrix is treated to contain a spectral vector,\n"
         "with frequencies repeated as assumed in *VectorToTbByRJ*. \n"
         "\n"
         "Generic output: \n"
         "   Matrix : A matrix with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Matrix : A matrix with radiance values."
        ),
        OUTPUT(),
        INPUT( sensor_pos_, sensor_los_, sensor_response_f_,
               sensor_response_za_, sensor_response_aa_,
               sensor_response_pol_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Matrix_ ),
        KEYWORDS(),
        TYPES()));

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
         "   Vector  : A set of frequencies. "
        ),
        OUTPUT(),
        INPUT( stokes_dim_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_ ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("mc_errorApplySensor"),
        DESCRIPTION
        (
         "Applies the sensor response on *mc_error*.\n"
         "\n"
         "The values in *mc_error* are weighted together assuming that the\n"
         "calculation error is independent between Stokes components, and \n"
         "is normal distributed. That is, if e=mc_error and H is the sensor\n"
         "response, the following calculation is performed (Matlab syntax):\n"
         "   e = sqrt( abs(H)*(e.*e) ) " 
        ),
        OUTPUT( mc_error_ ),
        INPUT( mc_error_, sensor_response_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( ),
        TYPES( )));
  
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
         "according to system time, positive rng_seed values are taken literally."),
        OUTPUT( y_, mc_iteration_count_, mc_error_ ),
        INPUT( f_grid_, rte_pos_, rte_los_, stokes_dim_, iy_space_agenda_,
               surface_prop_agenda_, opt_prop_gas_agenda_, 
               scalar_gas_absorption_agenda_, p_grid_, lat_grid_, lon_grid_, 
               z_field_, r_geoid_, z_surface_, t_field_, 
               vmr_field_, cloudbox_limits_, pnd_field_, scat_data_mono_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "std_err", "max_time", "max_iter", "rng_seed", "z_field_is_1D"),
        TYPES( Numeric_t, Index_t, Index_t, Index_t, Index_t)));

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
         "   value : The value." 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Numeric_ ),
        GINPUT(),
        KEYWORDS( "value"   ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("nelemGet"),
        DESCRIPTION
        (
         "Retrieve nelem from given variable and store the value in the \n"
         "workspace variable *nelem*"
        ),
        OUTPUT( nelem_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("ncolsGet"),
        DESCRIPTION
        (
         "Retrieve ncols from given variable and store the value in the \n"
         "workspace variable *ncols*"
        ),
        OUTPUT( ncols_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("nrowsGet"),
        DESCRIPTION
        (
         "Retrieve nrows from given variable and store the value in the \n"
         "workspace variable *nrows*"
        ),
        OUTPUT( nrows_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("npagesGet"),
        DESCRIPTION
        (
         "Retrieve npages from given variable and store the value in the \n"
         "workspace variable *npages*"
        ),
        OUTPUT( npages_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("nbooksGet"),
        DESCRIPTION
        (
         "Retrieve nbooks from given variable and store the value in the \n"
         "workspace variable *nbooks*"
        ),
        OUTPUT( nbooks_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("nshelvesGet"),
        DESCRIPTION
        (
         "Retrieve nshelves from given variable and store the value in the \n"
         "workspace variable *nshelves*"
        ),
        OUTPUT( nshelves_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("nvitrinesGet"),
        DESCRIPTION
        (
         "Retrieve nvitrines from given variable and store the value in the \n"
         "workspace variable *nvitrines*"
        ),
        OUTPUT( nvitrines_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("nlibrariesGet"),
        DESCRIPTION
        (
         "Retrieve nlibraries from given variable and store the value in the \n"
         "workspace variable *nlibraries*"
        ),
        OUTPUT( nlibraries_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( ),
        TYPES( ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

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
         "\n"
        ),
        OUTPUT( ext_mat_spt_, abs_vec_spt_ ),
        INPUT(  ext_mat_spt_, abs_vec_spt_, scat_data_raw_,
                scat_za_grid_, 
                scat_aa_grid_, scat_za_index_, scat_aa_index_, 
                f_index_, f_grid_, rte_temperature_,
                pnd_field_, scat_p_index_, scat_lat_index_, scat_lon_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "\n"
        ),
        OUTPUT( ext_mat_spt_, abs_vec_spt_ ),
        INPUT(  ext_mat_spt_, abs_vec_spt_, scat_data_mono_,
                scat_za_grid_, 
                scat_aa_grid_, scat_za_index_, scat_aa_index_, rte_temperature_,
                pnd_field_, scat_p_index_, scat_lat_index_, scat_lon_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME("output_file_formatSetAscii"),
        DESCRIPTION
        (
         "Sets the output file format to ASCII."
        ),
        OUTPUT( output_file_format_),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("output_file_formatSetBinary"),
        DESCRIPTION
        (
         "Sets the output file format to binary."
        ),
        OUTPUT( output_file_format_),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "\n"
         ),
        OUTPUT(scat_data_raw_, pnd_field_raw_),
        INPUT(atmosphere_dim_, f_grid_, p_grid_, lat_grid_, lon_grid_, 
              cloudbox_limits_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("filename_scat_data", "filename_pnd_field"),
        TYPES(String_t, String_t)));
 
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
         "\n"
         ),
        OUTPUT(scat_data_raw_, pnd_field_raw_),
        INPUT(atmosphere_dim_, f_grid_, p_grid_, lat_grid_, lon_grid_, 
              cloudbox_limits_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("filename_scat_data", "filename_pnd_field"),
        TYPES(String_t, String_t)));

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
       OUTPUT(scat_data_raw_, pnd_field_raw_),
       INPUT(),
       GOUTPUT(),
       GINPUT(),
       KEYWORDS(), 
       TYPES())); 

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
        OUTPUT(pha_mat_),
        INPUT(pha_mat_spt_, pnd_field_, atmosphere_dim_, scat_p_index_,
              scat_lat_index_, scat_lon_index_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES())); 

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
         "\n"
         ),
        OUTPUT(pha_mat_spt_),
        INPUT(pha_mat_spt_, scat_data_raw_, scat_za_grid_, scat_aa_grid_, 
              scat_za_index_, scat_aa_index_, f_index_, f_grid_, rte_temperature_,
              pnd_field_, scat_p_index_, scat_lat_index_, scat_lon_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES())); 

   md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_mat_sptFromMonoData" ),
        DESCRIPTION
        (
         "Calculation of the phase matrix for the single particle types.\n"
         "\n"
         "This function is the monchromatic version of *pha_mat_sptFromData*.\n"
         "\n"
         ),
        OUTPUT(pha_mat_spt_),
        INPUT(pha_mat_spt_, scat_data_mono_, doit_za_grid_size_, scat_aa_grid_, 
              scat_za_index_, scat_aa_index_, rte_temperature_,
              pnd_field_, scat_p_index_, scat_lat_index_, scat_lon_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES())); 

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
         "\n"
         ),
        OUTPUT(pha_mat_spt_),
        INPUT(pha_mat_spt_, pha_mat_sptDOITOpt_, scat_data_mono_, 
              doit_za_grid_size_,
              scat_aa_grid_, 
              scat_za_index_, scat_aa_index_, rte_temperature_,
              pnd_field_, scat_p_index_, scat_lat_index_, scat_lon_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES())); 

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
         "\n"
         ),
        OUTPUT(pnd_field_),
        INPUT(p_grid_, lat_grid_, lon_grid_, pnd_field_raw_, atmosphere_dim_,
              cloudbox_limits_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES())); 

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
         "*ppath_step_agenda* (type \"arts -d ppath_step_agenda\")."
        ),
        OUTPUT( ppath_, ppath_step_ ),
        INPUT( ppath_step_agenda_, atmosphere_dim_, p_grid_, lat_grid_, 
               lon_grid_, z_field_, r_geoid_, z_surface_, 
               cloudbox_on_, cloudbox_limits_, rte_pos_, rte_los_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));


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
         "points and points of surface intersections. The keyword *lmax* \n"
         "gives the option to include additional points to ensure that the\n"
         "distance along the path between the points does not exceed the \n"
         "selected maximum length. No additional points are included if\n"
         "*lmax* is set to <= 0.\n"
         "\n"
         "As functions of this kind should very seldom be called directly,\n"
         "and that the functions can be called a high number of times, these\n"
         "functions do not perform any checks of the input that give\n" 
         "detailed error messages, but asserts are performed (if turned on).\n"
         "\n"
         "For further information, type see the on-line information for\n"
         "*ppath_step_agenda* (type \"arts -d ppath_step_agenda\").\n"
         "\n"
         "Keywords: \n"
         "   lmax      : Maximum allowed length between path points."
        ),
        OUTPUT( ppath_step_ ),
        INPUT( ppath_step_, atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, 
               z_field_, r_geoid_, z_surface_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "lmax" ),
        TYPES(    Numeric_t )));

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
         "equations by the Euler method. Snell's law for a case with \n"
         "spherical symmetry is used for 1D to update the zenith angles. \n"
         "\n"
         "See further the on-line information for *ppath_stepGeometric*\n"
         "(type \"arts -d ppath_stepGeometric\") and the user guide for more\n"
         "details on the algorithms used.\n"
         "\n"
         "The maximum length of each ray tracing step is given by the \n"
         "keyword argument *lraytrace*. The length will never exceed the \n" 
         "given maximum value, but can be smaller. The ray tracing steps are\n"
         "only used to determine the path, points to describe the path for \n" 
         "*RteCalc* are included as for *ppath_stepGeometric*, this\n"
         "including the functionality for the keyword *lmax*.\n"
         "\n"
         "Keywords: \n"
         "   lraytrace : Maximum length of ray tracing steps.\n"
         "   lmax      : Maximum allowed length between path points."
        ),
        OUTPUT( ppath_step_, rte_pressure_, rte_temperature_, rte_vmr_list_, 
                refr_index_ ),
        INPUT( refr_index_agenda_, ppath_step_, atmosphere_dim_, p_grid_, 
               lat_grid_, lon_grid_, z_field_, t_field_, vmr_field_, r_geoid_,
               z_surface_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "lraytrace", "lmax"    ),
        TYPES(    Numeric_t,   Numeric_t )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("Print"),
        DESCRIPTION
        (
         "Prints a variable on the screen."
         "\n"
         "Keywords:\n"
         "   level : Output level to use. \n"
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( "level" ),
        TYPES( Index_t ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  // New name: p_gridFromAbsLookup
  md_data_raw.push_back     
    ( MdRecord
      ( NAME("p_gridFromGasAbsLookup"),
        DESCRIPTION
        (
         "Sets *p_grid* to the frequency grid of *gas_abs_lookup*."
        ),
        OUTPUT( p_grid_  ),
        INPUT(  gas_abs_lookup_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "   filename : Name of the input file."
         ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Any_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

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
         "The calculated values for a Tensor4, with size:\n"
         "   [atmosphere_dim+1, np, nlat, nlon] \n"
         "where np is the number of pressures given etc. The book of the\n"
         "tensor with the following index holds:\n"
         "   0: the refractive index \n"
         "   1: radial gradient of the refractive index \n"
         "   2: latitude gradient of the refractive index \n"
         "   3: longitude gradient of the refractive index \n"
         "\n"
         "To calculate these quantities for the atmsopheric mesh, execute:\n"
         "   RefrIndexFieldAndGradients(tensor4_1,p_grid,lat_grid,lon_grid)"
        ),
        OUTPUT( refr_index_, rte_pressure_, rte_temperature_, rte_vmr_list_ ),
        INPUT( refr_index_agenda_, atmosphere_dim_, p_grid_, lat_grid_, 
               lon_grid_, r_geoid_, z_field_, t_field_, vmr_field_ ),
        GOUTPUT( Tensor4_ ),
        GINPUT( Vector_, Vector_, Vector_  ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("refr_indexIR"),
        DESCRIPTION
        (
         "Calculates the microwave refractive index due to gases in the\n"
         "Earth's atmosphere. \n"
         "\n"
         "Only refractivity of dry air is considered. All other gases has\n"
                 "a negligible contribution.  \n"
         "\n"
         "The formula used is contributed by Michael Hï¿½fner,\n"
                 "Forschungszentrum Karlsruhe."
        ),
        OUTPUT( refr_index_ ),
        INPUT( rte_pressure_, rte_temperature_, rte_vmr_list_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "other gases has a negligible contribution.  \n"
         "\n"
         "The parameterisation of Thayer (Radio Science, 9, 803-807, 1974)\n"
         "is used. See also Eq. 3 and 5 of Solheim et al. (JGR, 104, \n"
         "pp. 9664). "
        ),
        OUTPUT( refr_index_ ),
        INPUT( rte_pressure_, rte_temperature_, rte_vmr_list_, gas_species_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "*refr_index_agenda*."
        ),
        OUTPUT( refr_index_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "See further the user guide."
        ),
        OUTPUT( y_, ppath_, ppath_step_, 
                iy_, rte_pos_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_, rte_los_,
                jacobian_, rte_do_vmr_jacs_, diy_dvmr_, rte_do_t_jacs_, 
                diy_dt_, ppath_array_do_, ppath_array_, ppath_array_index_  ),
        INPUT( ppath_step_agenda_, rte_agenda_, iy_space_agenda_,
               iy_surface_agenda_, iy_cloudbox_agenda_,
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_, 
               t_field_, vmr_field_, gas_species_, r_geoid_, z_surface_, 
               cloudbox_on_, cloudbox_limits_, sensor_response_, sensor_pos_, 
               sensor_los_, f_grid_, stokes_dim_, 
               antenna_dim_, mblock_za_grid_, mblock_aa_grid_, 
               jacobian_, jacobian_quantities_, jacobian_indices_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "RteCalcNoJacobian" ),
        DESCRIPTION
        (
         "As *RteCalc* but thorughout ignores jacobians."
        ),
        OUTPUT( y_, ppath_, ppath_step_, 
                iy_, rte_pos_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_, rte_los_ ),
        INPUT( ppath_step_agenda_, rte_agenda_, iy_space_agenda_,
               iy_surface_agenda_, iy_cloudbox_agenda_,
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_, 
               t_field_, vmr_field_, r_geoid_, z_surface_, 
               cloudbox_on_, cloudbox_limits_, sensor_response_, sensor_pos_, 
               sensor_los_, f_grid_, stokes_dim_, 
               antenna_dim_, mblock_za_grid_, mblock_aa_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
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
         "See further the user guide."
        ),
        OUTPUT( iy_, emission_, abs_scalar_gas_, 
                rte_pressure_, rte_temperature_, rte_vmr_list_, f_index_, 
                diy_dvmr_, diy_dt_ ),
        INPUT( iy_, ppath_, ppath_array_, ppath_array_index_, f_grid_, 
               stokes_dim_, emission_agenda_, scalar_gas_absorption_agenda_, 
               rte_do_vmr_jacs_, rte_do_t_jacs_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "RteStdWithTransmissions" ),
        DESCRIPTION
        (
         "As *RteStd*, but also returns path transmissions.\n"
         "\n"
         "The transmission to each point of the propagation path is returned\n"
         "in *ppath_transmissions*."
        ),
        OUTPUT( iy_, emission_, abs_scalar_gas_,
                rte_pressure_, rte_temperature_, rte_vmr_list_, f_index_, 
                ppath_transmissions_, diy_dvmr_, diy_dt_ ),
        INPUT( iy_, ppath_, ppath_array_, ppath_array_index_, f_grid_, 
               stokes_dim_, emission_agenda_, scalar_gas_absorption_agenda_, 
               rte_do_vmr_jacs_, rte_do_t_jacs_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "   aa : Azimuth angle of sensor line-of-sight."
        ),
        OUTPUT( rte_los_ ),
        INPUT( atmosphere_dim_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "za",      "aa"      ),
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
         "taken from the WGS-84 reference ellipsiod.\n"
         "\n"
         "For 1D, the geoid radius is set to the radius of curvature of the\n"
         "WGS-84 ellipsiod for the position and observation direction \n"
         "described with *lat_1d* and *meridian_angle_1d*.\n"
         "For 2D and 3D, the geoid radius is set to the radius of the WGS-84\n"
         "ellipsiod for the latitude value in *rte_pos*."
        ),
        OUTPUT( rte_pos_ ),
        INPUT( rte_pos_, atmosphere_dim_, lat_1d_, meridian_angle_1d_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "given position is outside the latitude and longitude grids."
        ),
        OUTPUT( rte_pos_ ),
        INPUT( rte_pos_, atmosphere_dim_, lat_grid_, lon_grid_, r_geoid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "*rte_posAddGeoidWGS84* should also be called to obtain a radius as\n"
         "first element of *rte_pos*.\n"
         "\n"
         "Keywords: \n"
         "   r_or_z : Radius or geometrical altitude of sensor position.\n"
         "   lat : Latitude of sensor position.\n"
         "   lon : Longitude of sensor position."
        ),
        OUTPUT( rte_pos_ ),
        INPUT( atmosphere_dim_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "r_or_z",  "lat",     "lon"     ),
        TYPES(    Numeric_t, Numeric_t, Numeric_t )));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_posShift" ),
        DESCRIPTION
        (
         "shifts rte_pos and rte_los, and rte_gp_XXX to the end of ppath."
        ),
        OUTPUT( rte_pos_, rte_los_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_ ),
        INPUT( ppath_, atmosphere_dim_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "This function is a work in progress. Only 1D is currently supported"
        ),
        OUTPUT( rte_pos_, rte_los_, ppath_, ppath_step_ ),
        INPUT( atmosphere_dim_, p_grid_, z_field_, lat_grid_, lon_grid_,
               ppath_step_agenda_, r_geoid_, z_surface_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("tan_p"),
        TYPES(Numeric_t)));

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
         "   r : Radius of geoid sphere. See further above."
        ),
        OUTPUT( r_geoid_ ),
        INPUT( atmosphere_dim_, lat_grid_, lon_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "r" ),
        TYPES( Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "r_geoidWGS84" ),
        DESCRIPTION
        (
         "Sets the geoid radius to match the WGS-84 reference ellipsiod.\n"
         "\n"
         "For 1D, the geoid radius is set to the radius of curvature of the\n"
         "WGS-84 ellipsiod for the position and observation direction \n"
         "described with *lat_1d* and *meridian_angle_1d*.\n"
         "For 2D and 3D, *r_geoid* is set to the radius of the WGS-84\n"
         "ellipsiod for the crossing points of the latitude and longitude\n"
         "grids.\n"
         "\n"
         "Please note that the latitude grid must contain true latitudes\n"
         "if the function shall give correct result, and not just arbitrary\n"
         "orbit angles which is allowed for 2D cases."
        ),
        OUTPUT( r_geoid_ ),
        INPUT( atmosphere_dim_, lat_grid_, lon_grid_, lat_1d_,
               meridian_angle_1d_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
   ( MdRecord
      ( NAME( "ScatteringDoit" ),
        DESCRIPTION
        (
         "This method executes *doit_mono_agenda* for each frequency \n"
         "in *f_grid*. The output is the radiation field inside the cloudbox\n"
         "(*doit_i_field*) and on the cloudbox boundary (*scat_i_p* (1D), \n"
         "*scat_i_lat* and *scat_i_lon* (3D)).\n"
         "\n"
         ),
        OUTPUT(doit_i_field_, scat_i_p_, scat_i_lat_, scat_i_lon_),
        INPUT(f_grid_, scat_i_p_, scat_i_lat_, scat_i_lon_, doit_mono_agenda_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));
 
 md_data_raw.push_back
    ( MdRecord
      ( NAME( "ScatteringMonteCarlo" ),
        DESCRIPTION
        (
         "This method performs a single pencil beam monochromatic scattering\n"
         "calculation using a Monte Carlo algorithm \n"
         "\n"
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
         "The incoming_lookup keyword determines if incoming radiance is obtained from"
         "a precalculated grid (mc_incoming) or calculated on the fly"
          ),
        OUTPUT(ppath_, ppath_step_, 
               mc_error_, mc_iteration_count_, 
               rte_pos_, rte_los_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_, iy_, 
               rte_pressure_, rte_temperature_, 
               rte_vmr_list_, ext_mat_, abs_vec_, f_index_,mc_incoming_),
        INPUT(ppath_step_agenda_, atmosphere_dim_, p_grid_, lat_grid_,
              lon_grid_, z_field_, r_geoid_, z_surface_, cloudbox_limits_,
              stokes_dim_, rte_agenda_, iy_space_agenda_, iy_surface_agenda_,
              t_field_, f_grid_, opt_prop_gas_agenda_,
              scalar_gas_absorption_agenda_, vmr_field_,
              scat_data_mono_, pnd_field_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("std_err","max_time","max_iter","rng_seed","incoming_lookup","z_field_is_1D"),
        TYPES( Numeric_t, Index_t, Index_t, Index_t, Index_t, Index_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "scat_data_monoCalc" ),
        DESCRIPTION
        (
         "Interpolates scat_data_raw by frequency to give scat_data_mono\n"
         "\n"
         ),
        OUTPUT( scat_data_mono_ ),
        INPUT( scat_data_raw_ ,f_grid_, f_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "\n"
         ),
        OUTPUT( ),
        INPUT( scat_data_raw_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));
 
 md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensorOff" ),
        DESCRIPTION
        (
         "Sets sensor WSVs to obtain monochromatic pencil beam values.\n"
         "\n"
         "The variables are set as follows:\n"
         "   sensor_response : As returned by *sensor_responseInit*.\n"
         "   sensor_pol      : Identity matrix, with size following\n"
         "                     *stokes_dim*.\n"
         "   sensor_rot      : Length matching *sensor_pos/los*. All values\n"
         "                     are set 0.\n"
         "   antenna_dim     : 1.\n"
         "   mblock_za_grid  : Length 1, value 0."
         "   mblock_aa_grid  : Empty."
        ),
        OUTPUT( sensor_response_, sensor_response_f_, sensor_response_za_,
                sensor_response_aa_, sensor_response_pol_, sensor_rot_,
                antenna_dim_, mblock_za_grid_, mblock_aa_grid_, sensor_norm_ ),
        INPUT( atmosphere_dim_, stokes_dim_, sensor_pos_, sensor_los_,
               f_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "taken from the WGS-84 reference ellipsiod.\n"
         "\n"
         "For 1D, the geoid radius is set to the radius of curvature of the\n"
         "WGS-84 ellipsiod for the position and observation direction \n"
         "described with *lat_1d* and *meridian_angle_1d*.\n"
         "For 2D and 3D, the geoid radius is set to the radius of the WGS-84\n"
         "ellipsiod for the latitude values in *sensor_pos*."
        ),
        OUTPUT( sensor_pos_ ),
        INPUT( sensor_pos_, atmosphere_dim_, lat_1d_, meridian_angle_1d_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "given position is outside the latitude and longitude grids."
        ),
        OUTPUT( sensor_pos_ ),
        INPUT( sensor_pos_, atmosphere_dim_, lat_grid_, lon_grid_, r_geoid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
        OUTPUT( sensor_response_, sensor_response_za_ ),
        INPUT( sensor_response_f_, sensor_response_pol_, mblock_za_grid_,
               antenna_dim_, antenna_diagram_, sensor_norm_, antenna_los_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "  ArrayOfMatrix : The backend channel response."
        ),
        OUTPUT( sensor_response_, sensor_response_f_ ),
        INPUT( f_backend_, sensor_response_pol_, sensor_response_za_,
               sensor_norm_ ),
        GOUTPUT( ),
        GINPUT( ArrayOfMatrix_ ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responseInit"),
        DESCRIPTION
        (
         "Initialises the response block matrix to an identity matrix.\n"
         "\n"
         "The initialised matrix is a quadratic matrix with sidelength equal\n"
         "to the product of the length of *f_grid*, *mblock_za_grid*,\n"
         "*mblock_aa_grid* and the columns of *sensor_pol*."
        ),
        OUTPUT( sensor_response_, sensor_response_f_, sensor_response_za_,
                sensor_response_aa_, sensor_response_pol_  ),
        INPUT( f_grid_, mblock_za_grid_, mblock_aa_grid_, antenna_dim_,
               atmosphere_dim_, stokes_dim_, sensor_norm_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "       Matrix : The sideband filter response matrix."
        ),
        OUTPUT( sensor_response_, sensor_response_f_, f_mixer_ ),
        INPUT( sensor_response_pol_, sensor_response_za_, lo_, sensor_norm_ ),
        GOUTPUT( ),
        GINPUT( Matrix_ ),
        KEYWORDS( ),
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
         "The local oscillator frequencies is set by the WSV *lo* which\n"
         "should have the same length as number of polarisations.\n"
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
         "         Matrix : The sideband filter response matrix."
         "  ArrayOfMatrix : The backend channel response."
        ),
        OUTPUT( sensor_response_, sensor_response_f_, f_mixer_ ),
        INPUT( sensor_response_pol_, sensor_response_za_, sensor_response_aa_,
               lo_, sensor_norm_, f_backend_, sensor_pol_ ),
        GOUTPUT( ),
        GINPUT( Matrix_, ArrayOfMatrix_ ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responsePolarisation"),
        DESCRIPTION
        (
         "Adds polarisation to the response matrix.\n"
         "\n"
         "The output polarisations are given by matrix *sensor_pol*."
        ),
        OUTPUT( sensor_response_, sensor_response_pol_ ),
        INPUT( sensor_pol_, sensor_response_za_, sensor_response_aa_,
               sensor_response_f_, stokes_dim_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "function and prior to sensor_responsePolarisation."
        ),
        OUTPUT( sensor_response_ ),
        INPUT( sensor_rot_, antenna_los_, antenna_dim_, stokes_dim_,
               sensor_response_f_, sensor_response_za_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "SparseToTbByPlanck" ),
        DESCRIPTION
        (
         "As *MatrixToTbByPlanck* but operates on a sparse matrix.\n"
         "\n"
         "Generic output: \n"
         "   Sparse : A matrix with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Sparse : A matrix with radiance values."
        ),
        OUTPUT(),
        INPUT( sensor_pos_, sensor_los_, sensor_response_f_,
               sensor_response_za_, sensor_response_aa_,
               sensor_response_pol_ ),
        GOUTPUT( Sparse_ ),
        GINPUT( Sparse_ ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "SparseToTbByRJ" ),
        DESCRIPTION
        (
         "As *MatrixToTbByRJ* but operates on a sparse matrix.\n"
         "\n"
         "Generic output: \n"
         "   Spasr : A sparse matrix with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Sparse : A sparse matrix with radiance values."
        ),
        OUTPUT(),
        INPUT( sensor_pos_, sensor_los_, sensor_response_f_,
               sensor_response_za_, sensor_response_aa_,
               sensor_response_pol_ ),
        GOUTPUT( Sparse_ ),
        GINPUT( Sparse_ ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("StringSet"),
        DESCRIPTION
        (
         "Sets a String to the given text String."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( String_ ),
        GINPUT(),
        KEYWORDS( "text"   ),
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
         "to hold blackbody radiation for a temperature of *surface_skin_t*."
        ),
        OUTPUT( surface_los_, surface_rmatrix_, surface_emission_ ),
        INPUT( f_grid_, stokes_dim_, surface_skin_t_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "surfaceCalc" ),
        DESCRIPTION
        (
         "General method for inclusion in *iy_surface_agenda*.\n"
         "\n"
         "The upwelling radiation, put into *iy*, is calculated as:\n"
         "   iy = *surface_emission* + sum_i[ i_down(i) *r(i) ]\n"
         "where i is LOS index, i_down(i) is downwelling radiation for \n"
         "direction i in *surface_los*, and r(i) is reflection matrix i in\n"
         "*surface_rmatrix*.\n"
         "\n"
         "See further the user guide."
        ),
        OUTPUT( iy_, ppath_, ppath_step_, 
                rte_pos_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_, rte_los_,
                ppath_array_do_, ppath_array_, ppath_array_index_,
                diy_dvmr_, diy_dt_ ),
        INPUT( ppath_step_agenda_, rte_agenda_, iy_space_agenda_,
               iy_surface_agenda_, iy_cloudbox_agenda_,
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_,
               t_field_, vmr_field_, r_geoid_, z_surface_, 
               cloudbox_on_, cloudbox_limits_, f_grid_, stokes_dim_, 
               rte_do_vmr_jacs_, rte_do_t_jacs_,
               surface_los_, surface_rmatrix_, surface_emission_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "surfaceFlat" ),
        DESCRIPTION
        (
         "Creates variables to mimic specular reflection by a surface with\n"
         "dielectric constant following an intrnal model.\n"
         "\n"
         "The method results in that the reflection properties differ\n"
         "between frequencies and polarisations. The properties of the\n"
         "surface medium are dtermined by the model for dielectric constant\n"
         "selected. Thermodynamic equilibrium is assumed, which in the\n"
         "corresponds to the reflection and emission coefficients add up\n"
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
         "   epsmodel : Name of model for dielectric constant."
        ),
        OUTPUT( surface_los_, surface_rmatrix_, surface_emission_ ),
        INPUT( f_grid_, stokes_dim_, atmosphere_dim_, rte_los_, 
               surface_skin_t_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "epsmodel" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "surfaceSingleEmissivity" ),
        DESCRIPTION
        (
         "Creates variables to mimic specular reflection by a surface with\n"
         "a single emissivity.\n"
         "\n"
         "A constant emissivity is assumed as a function of frequency and\n"
         "polarisation (vertical and horisontal reflection coefficients are\n"
         "equal. The number of directions in *surface_los* is one.\n"
         "\n"
         "Generic Input: \n"
         "   Numeric : Surface emissivity (a value between 0 and 1)."
        ),
        OUTPUT( surface_los_, surface_rmatrix_, surface_emission_ ),
        INPUT( f_grid_, stokes_dim_, atmosphere_dim_, rte_los_, 
               surface_skin_t_ ),
        GOUTPUT(),
        GINPUT( Numeric_t ),
        KEYWORDS(),
        TYPES()));

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
         "   ncols  : Number of columns in the tensor. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor3_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "npages", "nrows", "ncols"   ),
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
         "   value  : The scale factor. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor3_ ),
        GINPUT( Tensor3_ ),
        KEYWORDS( "value"   ),
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
         "   value  : The value of the tensor3 elements. "
        ),
        OUTPUT(),
        INPUT( npages_, nrows_, ncols_ ),
        GOUTPUT( Tensor3_ ),
        GINPUT(),
        KEYWORDS( "value"   ),
        TYPES(    Numeric_t )));

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
         "   value  : The scale factor. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor4_ ),
        GINPUT( Tensor4_ ),
        KEYWORDS( "value"   ),
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
         "   value  : The value of the tensor4 elements. "
        ),
        OUTPUT(),
        INPUT( nbooks_, npages_, nrows_, ncols_ ),
        GOUTPUT( Tensor4_ ),
        GINPUT(),
        KEYWORDS( "value"   ),
        TYPES(    Numeric_t )));

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
         "   value  : The scale factor. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor5_ ),
        GINPUT( Tensor5_ ),
        KEYWORDS( "value"   ),
        TYPES( Numeric_t )));

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
         "   value   : The value of the tensor5 elements. "
        ),
        OUTPUT(),
        INPUT( nshelves_, nbooks_, npages_, nrows_, ncols_ ),
        GOUTPUT( Tensor5_ ),
        GINPUT(),
        KEYWORDS( "value" ),
        TYPES(    Numeric_t )));

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
         "   value  : The scale factor. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor6_ ),
        GINPUT( Tensor6_ ),
        KEYWORDS( "value"   ),
        TYPES( Numeric_t )));

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
         "   value     : The value of the tensor6 elements. "
        ),
        OUTPUT(),
        INPUT( nvitrines_, nshelves_, nbooks_, npages_, nrows_, ncols_ ),
        GOUTPUT( Tensor6_ ),
        GINPUT(),
        KEYWORDS( "value" ),
        TYPES(    Numeric_t )));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor6ToTbByPlanck" ),
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
        OUTPUT(),
        INPUT(f_index_, f_grid_),
        GOUTPUT( Tensor6_ ),
        GINPUT( Tensor6_ ),
        KEYWORDS(),
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
         "   value  : The scale factor. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor7_ ),
        GINPUT( Tensor7_ ),
        KEYWORDS( "value"   ),
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
         "   value      : The value of the tensor7 elements. "
        ),
        OUTPUT(),
        INPUT( nlibraries_, nvitrines_, nshelves_, nbooks_, npages_, nrows_,
               ncols_ ),
        GOUTPUT( Tensor7_ ),
        GINPUT(),
        KEYWORDS( "value" ),
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
         "something out quickly."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
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
         "timerStart()"
         "ReadXML(f_grid){\"frequencies.xml\"}\n"
         "timerStop()\n"
         "Prints the CPU time spent for reading the XML file"
         ),
        OUTPUT(timer_),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));


  md_data_raw.push_back
    ( MdRecord
      ( NAME("timerStop"),
        DESCRIPTION
        (
         "Stops the CPU timer."
         "\n"
         "Use *timerStop* to output the consumed cpu time\n"
         "since *timerStart*. See *timerStart* for example"
         ),
        OUTPUT(),
        INPUT(timer_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
         "   value : The value to be added to the vector."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "value" ),
        TYPES( Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorLinSpace"),
        DESCRIPTION
        (
         "Creates a vector with linear spacing.\n"
         "\n"
         "The first element equals always the start value, and the spacing\n"
         "equlas always the step value, but note that the last value can  \n"
         "deviate from the stop value. The keyword step can be both positive\n"
         "and negative. \n"
         "   The vector is [start, start+step, start+2*step, ...]\n "
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Keywords:\n"
         "   start : The start value. \n"
         "    stop : The maximum value of the end value. \n"
         "    step : The spacing of the vector. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "start",   "stop",    "step"    ),
        TYPES(    Numeric_t, Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorNLinSpace"),
        DESCRIPTION
        (
         "Creates a vector with defined length, equally spaced between the \n"
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
         "       n : Number of elements of the vector. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(Vector_),
        GINPUT(),
        KEYWORDS( "start",   "stop",    "n"   ),
        TYPES(    Numeric_t, Numeric_t, Index_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorNLogSpace"),
        DESCRIPTION
        (
         "Creates a vector with defined length, equally logarithmically \n"
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
         "       n : Number of elements of the vector. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(Vector_),
        GINPUT(),
        KEYWORDS( "start",   "stop",    "n"   ),
        TYPES(    Numeric_t, Numeric_t, Index_t )));

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
         "   value : The value to be multiplicated with the vector."  
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "value" ),
        TYPES( Numeric_t )));

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
         "   value  : The value of the vector elements. " 
        ),
        OUTPUT(),
        INPUT( nelem_ ),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "value"   ),
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
         "   VectorSetExplicitly(p_grid){[1000 100 10]}\n"
         "   Will create a p_grid vector with these three elements."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "values"   ),
        TYPES(    Vector_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorToTbByPlanck" ),
        DESCRIPTION
        (
         "Converts a vector of radiances to brightness temperatures by \n"
         "inverting the Planck function.\n"
         "\n"
         "This function works as *VectorToTbByRJ*. However, this function \n"
         "is not recommended in connection with inversions, but can be used \n"
         "to display calculated spectra in a temperature scale.\n"
         "\n"
         "Generic output: \n"
         "   Vector : A vector with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A vector with radiance values. "
        ),
        OUTPUT(),
        INPUT( sensor_pos_, sensor_los_, sensor_response_f_,
               sensor_response_za_, sensor_response_aa_,
               sensor_response_pol_ ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorToTbByRJ" ),
        DESCRIPTION
        (
         "Converts a vector of radiances to brightness temperatures by \n"
         "the Rayleigh-Jeans approximation of the Planck function.\n"
         "\n"
         "This function performs a linear transformation of spectral \n"
         "radiances to an approximative temperature scale. The advantage \n"
         "of this linear transformation is that the obtained values can be \n"
         "used for retrievals if the weighting functions are handled \n"
         "likewise (by *MatrixToTbByRJ*). This is not the case if the \n"
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
         "   VectorToTbByRJ(y,y){} \n"
         "\n"
         "Generic output: \n"
         "   Vector : A vector with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A vector with radiance values."
        ),
        OUTPUT(),
        INPUT( sensor_pos_, sensor_los_, sensor_response_f_,
               sensor_response_za_, sensor_response_aa_,
               sensor_response_pol_ ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorZtanToZaRefr" ),
        DESCRIPTION
        (
         "Converts a set of true tangent altitudes to zenith angles.\n"
         "\n"
         "The tangent altitudes are given to the function as a vector, which\n"
         "are converted to a generic vector of zenith angles. The position of\n"
         "the sensor is given by the WSV *sensor_pos*. The function works only\n"
         "for a spherical geoid. The zenith angles are always set to be\n"
                 "positive.\n"
         "The tangent altitudes are given as the altitude above the geoid.\n"
         "\n"
         "Generic output: \n"
         "   Vector : A vector with zenith angles. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A vector with true tangent altitudes\n"
        ),
        OUTPUT( refr_index_, rte_pressure_, rte_temperature_, rte_vmr_list_ ),
        INPUT( refr_index_agenda_, sensor_pos_, p_grid_, t_field_, z_field_,
                           vmr_field_, r_geoid_, atmosphere_dim_ ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS(),
        TYPES()));

md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorZtanToZa" ),
        DESCRIPTION
        (
         "Converts a set of geometrical tangent altitudes to zenith angles.\n"
         "\n"
         "The tangent altitudes are given to the function as a vector, which\n"
         "are converted to a generic vector of zenith angles. The position of\n"
         "the sensor is given by the WSV *sensor_pos*. The function works only\n"
         "for a spherical geoid, where the geoid radius is given as a keyword\n"
         "argument (*r_geoid*). The zenith angles are always set to be positive.\n"
         "The tangent altitudes are given as the altitude above the geoid.\n"
         "\n"
         "Generic output: \n"
         "   Vector : A vector with zenith angles. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A vector with geometric tangent altitudes\n"
         "\n"
         "Keywords:\n"
         "   r_geoid : The geoid radius for the given tangent altitudes."
        ),
        OUTPUT(),
        INPUT( sensor_pos_ ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "r_geoid" ),
        TYPES( Numeric_t )));

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
         "   filename : Name of the output file."
         ),
        OUTPUT(),
        INPUT( output_file_format_ ),
        GOUTPUT( ),
        GINPUT(  Any_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

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
         "   filename : Name of the output file."
         ),
        OUTPUT(),
        INPUT( output_file_format_, file_index_ ),
        GOUTPUT( ),
        GINPUT(  Any_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ybatchCalc" ),
        DESCRIPTION
        (
         "Performs batch calculations."
         "\n"
         "The method performs the following:"
         "   1. Performs a-d with *ybatch_index* = 0 : (*ybatch_n*-1).\n"
         "    a. Executes *batch_update_agenda*.\n"
         "    b. Executes *batch_calc_agenda*.\n"
         "    c. If *ybatch_index* = 0, allocates memory for *ybatch* based\n"
         "       on *ybatch_n* and length of *y*.\n"
         "    d. Makes copy of *y* in column *ybatch_index* of *ybatch*.\n"
         "   2. Executes *batch_post_agenda*.\n"
         "This means that, beside involved agendas, the WSV *ybatch_n* must\n"
         "be set before calling this method.\n"
         "\n"
         "See the user guide for practical examples.\n"
         "\n"
         "Note that *y* and other variables modified by the involved agendas\n"
         "(e.g. *vmr_field*) are not restored by the method. If original\n"
         "values must be preserved, a copy must be made before calling this"
         "method, and a reversed copy must be made in *batch_post_agenda*\n"
         "or after executing *ybatchCalc*."
         ),
        OUTPUT( ybatch_, ybatch_index_, ybatch_n_, y_  ),
        INPUT( batch_update_agenda_, batch_calc_agenda_, batch_post_agenda_ ), 
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
        OUTPUT( ybatch_, y_, t_field_raw_, z_field_raw_, vmr_field_raw_, 
                pnd_field_raw_,  p_grid_, sensor_los_,cloudbox_on_, 
                cloudbox_limits_, z_surface_),
        INPUT(gas_species_, met_profile_calc_agenda_, f_grid_, met_amsu_data_,
              sensor_pos_, r_geoid_, lat_grid_, lon_grid_, atmosphere_dim_,
              scat_data_raw_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("nelem_p_grid", "met_profile_path", "met_profile_pnd_path"),
        TYPES(Index_t, String_t, String_t)));

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
        OUTPUT( ybatch_, t_field_raw_, z_field_raw_, vmr_field_raw_, 
                 y_, p_grid_, sensor_los_, z_surface_),
        INPUT(gas_species_, met_profile_calc_agenda_, 
              f_grid_, met_amsu_data_, sensor_pos_, r_geoid_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("nelem_p_grid","met_profile_path" ),
        TYPES(Index_t, String_t)));

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
         "*z_scan_low* and *z_scan_high*."
        ),
        OUTPUT( ppath_step_),
        INPUT( ppath_step_agenda_, atmosphere_dim_, p_grid_, lat_grid_,
               lon_grid_, z_field_, r_geoid_, z_surface_ ),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "z_recieve", "z_send", "t_sample", 
                  "z_scan_low", "z_scan_high" ),
        TYPES( Numeric_t, Numeric_t, Numeric_t,
               Numeric_t, Numeric_t )));


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
//          "   abs    : absorption coefficients [1/m], \n"
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
//          "                         pre-defined models is used."
//         ),
//         OUTPUT( ext_mat_zee_, abs_vec_zee_),
//         INPUT(  geomag_los_, f_grid_, 
//                 zeeman_o2_onoff_, zeeman_o2_pressure_limit_, zeeman_o2_line_,
//                 ppath_index_, rte_pressure_, 
//                 rte_temperature_,rte_vmr_list_, species_index_,
//                 abs_model_, abs_user_parameters_, stokes_dim_ ),
//         GOUTPUT( ),
//         GINPUT( ),
//         KEYWORDS( ),
//         TYPES( )));
  
//  md_data_raw.push_back
//     ( MdRecord
//       ( NAME("abs_vecAddGasZeeman"),
//         DESCRIPTION
//         (
//          "Add zeeman absorption to the elements of absorption vector.\n"
//          "\n"
//          ),
//         OUTPUT(abs_vec_),
//         INPUT(abs_vec_, abs_vec_zee_),
//         GOUTPUT(),
//         GINPUT(),
//         KEYWORDS(),
//         TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("ext_matAddGasZeeman"),
//         DESCRIPTION
//         (
//          "Add Zeeman extinction  to the elements of extinction matrix.\n"
//          " \n"
//          ),
//         OUTPUT(ext_mat_),
//         INPUT(ext_mat_, ext_mat_zee_),
//         GOUTPUT( ),
//         GINPUT( ),
//         KEYWORDS( ),
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
//         OUTPUT( zeeman_o2_onoff_, zeeman_o2_pressure_limit_, zeeman_o2_line_),
//         INPUT(),
//         GOUTPUT(),
//         GINPUT(),
//         KEYWORDS( "ZeemanO2OnOff", "ZeemanO2PressureLimit","ZeemanO2Line" ),
//         TYPES(    Index_t,       Numeric_t,               Index_t)));


//  md_data_raw.push_back
//     ( MdRecord
//       ( NAME("test_zeeman"),
//         DESCRIPTION(
//                     "\n"
//                     ),
//         OUTPUT(),
//         INPUT(opt_prop_gas_agenda_),
//         GOUTPUT( ),
//         GINPUT(),
//         KEYWORDS( ),
//         TYPES()));


}

