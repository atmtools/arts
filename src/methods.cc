/* Copyright (C) 2000, 2001, 2002
   Stefan Buehler <sbuehler@uni-bremen.de>
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
         "(e.g. global input) remove that part totally. Note that the \n"
         "on-line help just gives the type of global input/output and the \n"
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
               refr_index_, r_geoid_, z_ground_ ),
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


 //  md_data_raw.push_back
//     ( MdRecord
//       ( NAME("abs_vec_gasExample"),
//      DESCRIPTION
//      (
//       "This is only an example method created to perform test \n"
//       "calculations.\n"
//       ),
//      OUTPUT(abs_vec_),
//      INPUT(p_grid_,  atmosphere_dim_, stokes_dim_, scat_p_index_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));
  
  md_data_raw.push_back
    ( MdRecord
      ( NAME("amp_matCalc"),
        DESCRIPTION
        (
         "This method converts the raw amplitude matrix data namely \n"
         "amp_mat_raw to workspace variable amp_mat which can be directly\n"
         "used for calculating scattering properties of particles.  \n"
         "\n"
         "The data type of amp_mat_raw is an ArrayOfArrayOfTensor6 which \n"
         "contains one gridded field for each particle type the user wants.\n"
         "One data file contains amp_mat_raw for one particle type.  One \n"
         "data file has amp_mat_raw for one particle for a set of \n"
         "frequencies, incoming and outgoing zenith angles and azimuth \n"
         "angles. The frequencies, angles, and the amplitude matrix are \n"
         "each a Tensor 6. The size of amp_mat_raw is amp_mat_raw[Npt][7],\n"
         "where Npt is the number of particle types. amp_mat_raw[Npt] [0]\n"
         "gives the frequency tensor [N_f, 1, 1, 1, 1, 1] where N_f gives\n"
         "the number of frequencies considered in that particular database \n"
         "file. Similarly, amp_mat_raw[ Npt ][ 1 ] gives the outgoing zenith\n"
         "angle tensor [1, Nza, 1, 1, 1, 1], amp_mat_raw[ Npt ][ 2 ] gives \n"
         "the outgoing azimuth angle tensor [1, 1, Naa, 1, 1, 1], \n"
         "amp_mat_raw[ Npt ][ 3 ] gives the incoming zentih angle tensor\n"
         "[1, 1, 1, Nza, 1, 1], amp_mat_raw[ Npt ][ 4 ] gives the incoming\n"
         "azimuth angle tensor [1, 1, 1, 1, Naa, 1], amp_mat_raw[ Npt ][ 5 ]\n"
         "is a dummy tensor6 and amp_mat_raw[ Npt ][ 6 ] gives amplitude\n"
         "matrix which is also a tensor6 of size \n"
         "[N_f, N_za, N_aa, N_za, N_aa, 8]. Here, Nza is the number of \n"
         "zenith angles, Naa is the number of azimuth angles and 8 denotes \n"
         "the amplitude matrix elements.  \n"
         "\n"
         "In this method, we have to interpolate the raw data calculated on \n"
         "specific angular and frequency grids onto a grid which is set by \n"
         "the user. Since we decide that frequency should be the outermost \n"
         "loop for our radiative transfer calculation, the frequency grid \n"
         "contains just one value specified by the index scat_f_index. The\n"
         "angles for which the calculations are to be done are specified by\n"
         "scat_za_grid and scat_aa_grid. The interpolation has to be done \n"
         "for the frequency grid, zenith angle grid and azimuth angle grid. \n"
         "Since this interpolation is from a gridded field to a new field, \n"
         "we have to perform a green interpolation. For more insight into \n"
         "the interpolation schemes refer to Chapter 8-Interpolation of AUG.\n"
         "\n"
         "The output of this method is amp_mat has to be a Tensor6 with the \n"
         "first dimension being that of the particle type, then the angles \n"
         "and finally the amplitude matrix element 8. The size of amp_mat is\n"
         "(Npt, Nza, Naa, Nza, Naa, 8).  Note that the dimension frequency \n"
         "is taken out.\n"
         ),
        OUTPUT(amp_mat_),
        INPUT(amp_mat_raw_, f_index_, f_grid_, scat_za_grid_, 
              scat_aa_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES())); 

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
                a_pressure_, a_temperature_, a_vmr_list_ ),
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
         "*a_pressure*, *a_temperature*, and *a_vmr_list*, and returns the\n"
         "output variable *abs_scalar_gas*."
        ),
        OUTPUT( abs_scalar_gas_field_,
                abs_scalar_gas_,
                a_pressure_, a_temperature_, a_vmr_list_),
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
      ( NAME("abs_vec_sptCalc"),
        DESCRIPTION
        (
         "This method calculates the absorption vector for a single particle\n"
         "type.\n"
         "\n"
         "All elements of absorption vector for a particle can be expressed\n"
         "in terms of extinction matrix and phase matrix. So it is necessary\n"
         "that the methods ext_mat_sptCalc and pha_mat_sptCalc are done \n"
         "before calling this method. \n"
         "\n"
         "The output of the method *abs_vec_sptCalc* is *abs_vec_spt* \n"
         "(Matrix, size: [Npt,stokes_dim]). The input to the method \n"
         "*abs_vec_sptCalc are abs_vec_spt, *pha_mat_spt*(Tensor 5,*\n"
         "size: [Npt,Nza,Naa,stokes_dim,stokes_dim]), *ext_mat_spt* \n"
         "(Tensor 3,size = [Npt,stokes_dim,stokes_dim]), *scat_za_grid*\n"
         "(Vector,size = [Nza]) and *scat_aa_grid* (Vector,size = [Naa]).\n"
         "This method calls the function amp2abs which does the actual\n"
         "physics, that of computing the elements of absorption vector \n"
         "from the elements of extinction matrix and phase matrix.\n"
         ),
        OUTPUT(abs_vec_spt_ ),
        INPUT(abs_vec_spt_, ext_mat_spt_, pha_mat_spt_, scat_za_grid_,
              scat_aa_grid_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
      ( NAME("ArrayOfIndexPrint"),
        DESCRIPTION
        (
         "Prints the value of an ArrayOfIndex variable on the screen."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ArrayOfIndex_ ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("ArrayOfStringPrint"),
        DESCRIPTION
        (
         "Prints the value of an ArrayOfString variable on the screen."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ArrayOfString_ ),
        KEYWORDS( ),
        TYPES( )));
 
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
      ( NAME("AtmosphereSet1D"),
        DESCRIPTION
        (
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
      ( NAME("AtmFieldsCalc"),
        DESCRIPTION
        (
         "Interpolated the atmospheric fields.\n"
         "\n"
         "An atmospheric scenario includes the following data for each \n"
         "position (pressure, latitude, longitude) in the atmosphere: \n"
         "           1. temperature field \n"
         "           2. the corresponding altitude field \n"
         "           3. vmr fields for the gaseous species \n"
         "This methods interpolates the fields from the raw data\n"
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
         "Interpolate 1D raw atmospheric fields to create  2D or 3D \n"
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
      ( NAME("AtmFieldsFromAscii2Xml"),
        DESCRIPTION
        (
         "Converts atmospheric data profiles (ascii-format) to xml format.\n"
         "\n"
         "This function can be used to convert data into the right format.\n"
         "The input files are matrices, holding the pressure grid and the \n"
         "data and output is a 3D gridded field (ArrayOfTensor3). We need \n"
         "this format, because generally we store the data not only on the \n"
         "pressure grid, but also on latitude and longitude grid. \n"
         "\n"
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("path", "basename"),
        TYPES(String_t, String_t)));

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
      ( NAME( "a_losSet" ),
        DESCRIPTION
        (
         "Sets *a_los* to the given angles.\n"
         "\n"
         "The keyword argument *za* is put in as first element of *a_los*\n"
         "and *aa* as the second element. However, when *atmosphere_dim* is\n"
         "set to 1D or 2D, the length of *a_los* is set to 1 and only the\n"
         "given zenith angle is considered.\n"
         "\n"
         "Keywords: \n"
         "   za : Zenith angle of sensor line-of-sight.\n"
         "   aa : Azimuth angle of sensor line-of-sight."
        ),
        OUTPUT( a_los_ ),
        INPUT( atmosphere_dim_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "za",      "aa"      ),
        TYPES(    Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "a_posAddGeoidWGS84" ),
        DESCRIPTION
        (
         "Adds a geoid radius according to WGS-84 to a geometric altitude.\n"
         "\n"
         "This function assumes that the first element of *a_pos* is set\n"
         "to the geometric altitude for the position of the sensor. \n"
         "The variable *a_pos* shall contain the radius instead of the\n"
         "altitude and that can be achieved by this function. The function\n"
         "adds a geoid radius to the given altitude. The geoid radius is\n"
         "taken from the WGS-84 reference ellipsiod.\n"
         "\n"
         "For 1D, the geoid radius is set to the radius of curvature of the\n"
         "WGS-84 ellipsiod for the position and observation direction \n"
         "described with *lat_1d* and *meridian_angle_1d*.\n"
         "For 2D and 3D, the geoid radius is set to the radius of the WGS-84\n"
         "ellipsiod for the latitude value in *a_pos*."
        ),
        OUTPUT( a_pos_ ),
        INPUT( a_pos_, atmosphere_dim_, lat_1d_, meridian_angle_1d_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "a_posAddRgeoid" ),
        DESCRIPTION
        (
         "Adds a geoid radius by interpolating *r_geoid*.\n"
         "\n"
         "This function assumes that the first element of *a_pos* is set\n"
         "to the geometric altitude for the position of the sensor. \n"
         "The variable *a_pos* shall contain the radius instead of the\n"
         "altitude and that can be achieved by this function. The function\n"
         "adds a geoid radius to the given altitude. The geoid radius is\n"
         "obtained by interpolation of *r_geoid*. There is an error if the\n"
         "given position is outside the latitude and longitude grids."
        ),
        OUTPUT( a_pos_ ),
        INPUT( a_pos_, atmosphere_dim_, lat_grid_, lon_grid_, r_geoid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "a_posSet" ),
        DESCRIPTION
        (
         "Sets *a_pos* to the given co-ordinates.\n"
         "\n"
         "The keyword argument *r_or_z* is put in as first element of\n"
         "*a_pos*, *lat* as the second element and *lon* as third element.\n"
         "However, the length of *a_pos* is set to *atmosphere_dim* and\n"
         "keyword arguments for dimensions not used are ignored.\n"
         "\n"
         "The first keyword argument can either be a radius, or an altitude\n"
         "above the geoid. In the latter case, a function such as\n"
         "*a_posAddGeoidWGS84* should also be called to obtain a radius as\n"
         "first element of *a_pos*.\n"
         "\n"
         "Keywords: \n"
         "   r_or_z : Radius or geometrical altitude of sensor position.\n"
         "   lat : Latitude of sensor position.\n"
         "   lon : Longitude of sensor position."
        ),
        OUTPUT( a_pos_ ),
        INPUT( atmosphere_dim_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "r_or_z",  "lat",     "lon"     ),
        TYPES(    Numeric_t, Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "CloudboxGetIncoming" ),
        DESCRIPTION
        (
         "This method gives the intensity field at the boundary of the\n"
         "cloudbox.\n"
         "\n"
         "The method uses *PpathCalc* and *rte_agenda*.  The input to\n"
         "this method is the position of the cloudbox given by the \n"
         "variable *cloudbox_limits*. Then for each propagation direction\n"
         "it calls the function *PpathCalc* and executes the agenda \n"
         "*rte_agenda*.  This gives *i_rte* which holds the Stokes vector\n"
         "for all frequencies btained by the monochromatic pencil beam \n"
         "calculations performed by *rte_agenda*. Then this is copied to\n"
         "the interface variable. \n"
         "\n"
         "This is implemented only for the case of 1D atmosphere at present\n"
         "\n"
         ),
        OUTPUT(scat_i_p_, scat_i_lat_, scat_i_lon_, ppath_, ppath_step_,
               i_rte_, y_rte_, i_space_, ground_emission_, ground_los_, 
               ground_refl_coeffs_,
               mblock_index_, a_los_, a_pos_, a_gp_p_, a_gp_lat_, a_gp_lon_),
        INPUT( cloudbox_limits_, atmosphere_dim_, stokes_dim_, scat_za_grid_,
                scat_aa_grid_, f_grid_, ppath_step_agenda_,  rte_agenda_,
                i_space_agenda_, ground_refl_agenda_, p_grid_, lat_grid_,
                lon_grid_, z_field_, t_field_, r_geoid_, z_ground_),
              
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "CloudboxGetOutgoing" ),
        DESCRIPTION
        (
         "Scattered radiance on the cloudbox boundary.\n"
         "\n"
         "This method returns the radiances for a given direction and \n"
         "position on the boundary of the cloudbox. The variable *y_scat* \n"
         "is a matrix with the dimensions [f_grid, stokes_dim].\n"
         "\n"
         ),
        OUTPUT(),
        INPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, a_gp_p_, a_gp_lat_, 
               a_gp_lon_, a_los_,  cloudbox_on_, cloudbox_limits_, 
               atmosphere_dim_, stokes_dim_, 
               scat_za_grid_, scat_aa_grid_, f_grid_),
        GOUTPUT(Matrix_),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("CloudboxOff"),
        DESCRIPTION
        (
         "Deactivates the cloud box. \n"
         "\n"
         "The function sets *cloudbox_on* to 0, and *cloudbox_limits* to be\n"
         "an empty vector. The variables *scat_i_p*, *scat_i_lat*,  \n"
         "*scat_i_lon*, *scat_za_grid* and *scat_aa_grid* are also set to be\n"
         "empty."
        ),
        OUTPUT( cloudbox_on_, cloudbox_limits_, scat_i_p_, scat_i_lat_,
                scat_i_lon_, scat_za_grid_, scat_aa_grid_ ),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "CloudboxSetManually" ),
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
        OUTPUT( cloudbox_on_, cloudbox_limits_ ),
        INPUT( atmosphere_dim_, p_grid_, lat_grid_, lon_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "p1", "p2", "lat1", "lat2", "lon1", "lon2" ),
        TYPES( Numeric_t, Numeric_t, Numeric_t, Numeric_t, Numeric_t, 
               Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("CloneSize"),
        DESCRIPTION
        (
         "Clone the size of a workspace variable.\n"
         "\n"
         "This is a supergeneric method. It can be applied to any two workspace\n"
         "variables of the same group. It makes the output variable the same\n"
         "size as the input variable and initializes it. (This means setting all\n"
         "elements to zero for numeric variables like vectors or matrices.)\n"
         "\n"
         "As always, output comes first in the argument list!\n"
         "\n"
         "Usage example:\n"
         "\n"
         "CloneSize(f_grid,p_grid){}\n"
         "\n"
         "Will make the Vector *f_grid* the same length as the Vector *p_grid*\n"
         "and set all elements to zero.\n"
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
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Copy"),
        DESCRIPTION
        (
         "Copy a workspace variable.\n"
         "\n"
         "This is a supergeneric method. It can copy any workspace variable to\n"
         "another workspace variable of the same group. (E.g., a Matrix to\n"
         "another Matrix.)\n"
         "\n"
         "As allways, output comes first in the argument list!\n"
         "\n"
         "Usage example:\n"
         "\n"
         "Copy(f_grid,p_grid){}\n"
         "\n"
         "Will copy the content of *p_grid* to *f_grid*. The size of *f_grid* is\n"
         "adjusted automatically (the normal behaviour for workspace methods).\n"
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
      ( NAME("convergence_flagAbs"),
        DESCRIPTION
        (
         "Convergence test (maximum absolute difference).\n"
         "\n"
         "The function calculates the absolute differences for two successive\n"
         "iteration fields. It picks out the maximum values for each Stokes \n"
         "component separately. The convergence test is fullfilled under the\n"
         "following conditions: \n"
         "|I(m+1) - I(m)| < epsilon_1     Intensity.\n"
         "|Q(m+1) - Q(m)| < epsilon_2     The other Stokes components.\n" 
         "|U(m+1) - U(m)| < epsilon_3    \n"
         "|V(m+1) - V(m)| < epsilon_4    \n" 
         "\n"
         "The limits for convergence have to be set in the controlfile by \n"
         "setting the vector *epsilon* to appropriate values.\n"
         "\n"
         "The conditions have to be valid for all positions in the cloudbox \n"
         "and for all directions.\n"  
         "Then *convergence_flag* is set to 1.\n" 
        ),
        OUTPUT(convergence_flag_),
        INPUT(i_field_, i_field_old_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS("epsilon"),
        TYPES(Vector_t)));

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
         "The resluling extinction matrix is added to the workspace \n"
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

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("ext_mat_partScat"),
//      DESCRIPTION
//      (
//       "This function sums up the convergence extinction matrices \n"
//          "for all particle types weighted with particle number density.\n"
//       "\n"
//          "This method is only used for the convergence test, where particle\n"
//          "absorption is set to be 0. \n"
//       "The output of this method is *ext_mat_part* (stokes_dim, stokes_dim).\n"
//       "The inputs are the convergence extinction matrix for the single\n"
//          "particle type *ext_mat_conv_spt* \n"
//          "(part_types, stokes_dim, stokes_dim) and the local \n"
//       "particle number densities for all particle types namely the \n"
//       "*pnd_field* (part_types, p_grid, lat_grid, lon_grid ) for given \n"
//       "*p_grid*, *lat_grid*, and *lon_grid*. The particle types required \n"
//       "are specified in the control file.  \n"
//       ),
//      OUTPUT( ext_mat_  ),
//         INPUT( ext_mat_spt_, pnd_field_, atmosphere_dim_, scat_p_index_, 
//             scat_lat_index_, scat_lon_index_),
//      GOUTPUT( ),
//      GINPUT( ),
//      KEYWORDS( ),
//      TYPES( )));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("ext_mat_sptScat"),
        DESCRIPTION
        ( 
         "This method calculates the convergence extinction matrix.\n"
         "\n"
         "This function calculates extinction due to scattering (without \n"
         "absorption). It is used only for testing, if the iterative method \n"
         "converges towards the right solution.\n"
         "\n"
         ),
        OUTPUT( ext_mat_spt_  ),
        INPUT( pha_mat_spt_, scat_za_grid_, scat_aa_grid_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));
 
 md_data_raw.push_back
    ( MdRecord
      ( NAME("ext_mat_sptCalc"),
        DESCRIPTION
        ( 
         "This method calculates the extinction matrix for a single particle\n"
         "type.\n"
         "\n"
         "Extinction matrix describes the total attenuation of the incident \n"
         "radiation resulting from the combined effect of scattering and\n"
         "absorption by the particle.  It is a 4X4 matrix and all the\n"
         "elements of extinction matrix for a particle can be expressed in\n"
         "terms of the elements of the forward scattering amplitude matrix.\n"
         "\n"
         "The output of this method is \n"
         "*ext_mat_spt*(Tensor 3, size = [Npt,stokes_dim,stokes_dim])and the\n"
         "inputs are *ext_mat_spt*,*amp_mat*(Tensor 6, Size=[Npt,Nza,Naa,Nza,Naa,8]), \n"
         "*scat_za_index*,*scat_aa_index*,*f_index* and *scat_f_grid*. \n"
         "\n"
         "The variables *scat_za_index* and *scat_aa_index picks the right \n"
         "element of the Tensor *amp_mat*. *f_grid* and *f_index* picks \n"
         "the right frequeny for calculation. Frequeny is needed because the\n"
         "computation of extinction matrix from amplitude matrix involves \n"
         "multiplication by wavelength.  Then this method calls the \n"
         "function amp2ext which does the actual physics, that of computing\n"
         "the elements of extinction matrix from the elements of amplitude \n"
         "matrix.\n"
         ),
        OUTPUT( ext_mat_spt_  ),
        INPUT( ext_mat_spt_,amp_mat_, scat_za_index_, scat_aa_index_,
               f_index_, f_grid_),
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
      ( NAME( "i_fieldIterate" ),
        DESCRIPTION
        (
         "Iterative solution of the RTE.\n"
         "\n"
         "A solution for the RTE with scattering is found using an iterative\n"
         "scheme:\n"
         "\n"
         "1. Calculate scattering integral.\n"
         "   To calculate the scattering integral the total scattering \n"
         "   matrix is required. This is obtained using the functions \n"
         "   *pha_mat_sptCalc* and  *pha_matCalc*. The method *.....* \n"
         "   performs the integration.\n"
         "2. Calculate RT with fixed scattered field.\n"
         "   The radiative transfer equation with fixed scattering integral\n"
         "   can be solved analytically if the coefficients are assumed to \n"
         "   be constant.\n"
         "   According to *atmosphere_dim* either *i_fieldUpdate1D* or \n"
         "   *i_fieldUpdate2D* are called to perform the calculation. Inside\n"
         "   these methods the agenda *scat_rte_agenda* is executed. \n"
         "3. Convergence test.\n"
         "   Here the *convergence_test_agenda* is executed.\n"
         "\n"
         "Note: The atmospheric dimensionality *atmosphere_dim* can be \n"
         "      either 1 or 3. To these dimensions the method adapts \n"
         "      automatically. \n"
         "      If *atmosphere_dim* equals 2, it returns an error message,\n"
         "      as 2D scattering calculations can not be performed.\n"
         ),
        OUTPUT(i_field_, i_field_old_, convergence_flag_, scat_field_, 
               pha_mat_, pha_mat_spt_, abs_scalar_gas_, a_pressure_,
               a_temperature_, a_vmr_list_, scat_za_index_, scat_aa_index_,
               ext_mat_, abs_vec_, scat_p_index_, scat_lat_index_, 
               scat_lon_index_, ppath_step_),
        INPUT(cloudbox_limits_, atmosphere_dim_, part_types_, amp_mat_, 
              p_grid_, lat_grid_, lon_grid_, convergence_test_agenda_, 
              scalar_gas_absorption_agenda_, vmr_field_, spt_calc_agenda_,
              scat_za_grid_, scat_aa_grid_, opt_prop_part_agenda_, 
              pnd_field_, opt_prop_gas_agenda_, ppath_step_agenda_, z_field_,
              r_geoid_, t_field_, f_grid_, f_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "i_fieldSetClearsky" ),
        DESCRIPTION
        (
         "This method uses a linear 3D interpolation scheme to obtain the\n"
         "radiation field on all grid points inside the cloud box from the\n"
         "clear sky field on the cloud bod boundary. This radiation field\n"
         "is taken as the first guess radiation field for the iterative \n"
         "solution method of the RTE.\n"
         "\n"
         "The inputs to this method are *scat_i_p*[ N_f, 2, N_lat, N_lon,\n"
         "N_za, N_aa, stokes_dim], *scat_i_lat*[ N_f, N_p, 2, N_lon, \n"
         "N_za, N_aa, stokes_dim] and *scat_i_lon*[ N_f, N_p, N_lat, 2,\n"
         "N_za, N_aa, stokes_dim].  The method has to pick the \n"
         "monochromatic radiation field out of these variables.  The \n"
         "output of the method is the initial field stored in the \n"
         "workspace variable *i_field*[ N_p, N_lat, N_lon, N_za, N_aa,\n"
         "stokes_dim].\n"
         ),
        OUTPUT(i_field_),
        INPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, f_grid_, 
               f_index_, p_grid_, lat_grid_, lon_grid_, 
               cloudbox_limits_, atmosphere_dim_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "i_fieldSetConst" ),
        DESCRIPTION
        (
         "This method can be used to set the initial field inside the \n"
         "cloudbox to a constant value. \n"
         "The value can be given by the user .\n"
         "\n"
         "The inputs to this method are *scat_i_p*[ N_f, 2, N_lat, N_lon,\n"
         "N_za, N_aa, stokes_dim], *scat_i_lat*[ N_f, N_p, 2, N_lon, \n"
         "N_za, N_aa, stokes_dim] and *scat_i_lon*[ N_f, N_p, N_lat, 2,\n"
         "N_za, N_aa, stokes_dim].  The method has to pick the \n"
         "monochromatic radiation field out of these variables.  The \n"
         "output of the method is the initial field stored in the \n"
         "workspace variable *i_field*[ N_p, N_lat, N_lon, N_za, N_aa,\n"
         "stokes_dim].\n"
         ),
        OUTPUT(i_field_),
        INPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, f_grid_, 
               f_index_, p_grid_, lat_grid_, lon_grid_, 
               cloudbox_limits_, atmosphere_dim_, stokes_dim_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("value"),
        TYPES(Vector_t)));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "i_fieldUpdate1D" ),
        DESCRIPTION
        (
         "Updates the i_field during the iteration. It performs the RT \n"
         "calculation using a fixed value for the scattering integral stored \n"
         "in *scat_field*.\n"
         "\n" 
        ),
        OUTPUT(i_field_, abs_scalar_gas_, a_pressure_, a_temperature_,
               a_vmr_list_, scat_za_index_, ext_mat_, abs_vec_,
               scat_p_index_, ppath_step_),
        INPUT(i_field_old_, scat_field_, cloudbox_limits_, 
              scalar_gas_absorption_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, 
              opt_prop_part_agenda_, pnd_field_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, z_field_, r_geoid_, t_field_,
              f_grid_, f_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "i_fieldUpdate3D" ),
        DESCRIPTION
        (
         "Updates the i_field during the iteration. It performs the RT \n"
         "calculation using a fixed value for the scattering integral stored \n"
         "in *scat_field*.\n"
         "\n " 
        ),
        OUTPUT(i_field_, abs_scalar_gas_, a_pressure_, a_temperature_,
               a_vmr_list_, scat_za_index_, scat_aa_index_, ext_mat_, abs_vec_,
               scat_p_index_, scat_lat_index_, scat_lon_index_,  ppath_step_),
        INPUT(i_field_old_, scat_field_, cloudbox_limits_, 
              scalar_gas_absorption_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, scat_aa_grid_,
              opt_prop_part_agenda_, pnd_field_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, lat_grid_, lon_grid_, z_field_,
              r_geoid_, t_field_,
              f_grid_, f_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "i_fieldUpdate1D_PlaneParallel" ),
        DESCRIPTION
        (
         "Updates the i_field during the iteration. It performs the RT \n"
         "calculation using a fixed value for the scattering integral stored \n"
         "in *scat_field*.\n"
         "   " 
        ),
        OUTPUT(i_field_, ppath_step_, stokes_vec_, 
               sca_vec_, a_planck_value_, l_step_,
               abs_vec_spt_, ext_mat_spt_, pha_mat_spt_, ext_mat_, abs_vec_,
               scat_p_index_, scat_za_index_, scat_aa_index_, abs_scalar_gas_),
        INPUT(spt_calc_agenda_, opt_prop_part_agenda_, opt_prop_gas_agenda_,
              scalar_gas_absorption_agenda_, ppath_step_agenda_,
              amp_mat_, scat_field_, cloudbox_limits_,
              scat_za_grid_, scat_aa_grid_, p_grid_, t_field_, z_field_, 
              r_geoid_, f_grid_, f_index_, 
              pnd_field_, stokes_dim_, atmosphere_dim_, part_types_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));




 md_data_raw.push_back     
    ( MdRecord
      ( NAME("IndexPrint"),
        DESCRIPTION
        (
         "Prints the value of an Index variable on the screen."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Index_ ),
        KEYWORDS( ),
        TYPES( )));

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
      ( NAME("iteration_counterIncrease"),
        DESCRIPTION
        (
         "Increase iteration counter. \n"
         "\n"
         "This function can be used for writing the separate iteration \n"
         "fields into differtent files using *Tensor6WriteIteration*.\n"
         ),
        OUTPUT(iteration_counter_),
        INPUT(iteration_counter_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));
  
   md_data_raw.push_back
    ( MdRecord
      ( NAME("elsDoppler"),
        DESCRIPTION
        (
         "The Doppler lineshape.\n"
         "\n"
         "This computes the Doppler lineshape as:\n"
         "\n"
         "els[i] = 1/(sqrt(PI) * ls_sigma) /\n"
         "         exp(-(els_f_grid[i]/ls_sigma)^2)\n"
         "\n"
         "Note that the frequency grid els_f_grid must hold\n"
         "offset frequencies from line center. Hence, the\n"
         "line center frequency is not needed as input.\n"
         "\n"
         "Output:\n"
         "   els        : The lineshape function [1/Hz]\n"
         "\n"
         "Input:\n"
         "   ls_sigma   : Line width [Hz].\n"
         "   els_f_grid : Frequency grid [Hz]."
        ),
        OUTPUT( els_ ),
        INPUT(  ls_sigma_, els_f_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES())); 

  md_data_raw.push_back
    ( MdRecord
      ( NAME("elsLorentz"),
        DESCRIPTION
        (
         "The Lorentz lineshape.\n"
         "\n"
         "This computes the simple Lorentz lineshape as:\n"
         "\n"
         "els[i] = 1/PI * ls_gamma /\n"
         "         ( (els_f_grid[i])^2 + ls_gamma^2 )\n"
         "\n"
         "Note that the frequency grid els_f_grid must hold\n"
         "offset frequencies from line center. Hence, the\n"
         "line center frequency is not needed as input.\n"
         "\n"
         "Output:\n"
         "   els        : The lineshape function [1/Hz]\n"
         "\n"
         "Input:\n"
         "   ls_gamma   : Line width [Hz].\n"
         "   els_f_grid : Frequency grid [Hz]."
        ),
        OUTPUT( els_ ),
        INPUT(  ls_gamma_, els_f_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("elsVoigt_Drayson"),
        DESCRIPTION
        (
         "The Voigt Drayson linshape.\n"
         "\n"
         "This computes the Voigt profile\n"
         "through the Drayson approximation.\n"
         "\n"
         "Note that the frequency grid els_f_grid must hold\n"
         "offset frequencies from line center. Hence, the\n"
         "line center frequency is not needed as input.\n"
         "\n"
         "Output:\n"
         "   els        : The lineshape function [1/Hz]\n"
         "\n"
         "Input:\n"
         "   ls_sigma   : Lorentz width [Hz].\n"
         "   ls_gamma   : Doppler width [Hz].\n"
         "   els_f_grid : Frequency grid [Hz]."
        ),
        OUTPUT( els_ ),
        INPUT(  ls_sigma_, ls_gamma_, els_f_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("elsVoigt_Kuntz6"),
        DESCRIPTION
        (
         "The Voigt-Kuntz linshape.\n"
         "\n"
         "This computes the Voigt profile\n"
         "through the Kuntz algorithm with\n"
         "a relative accuracy better than 2*10-6.\n"
         "\n"
         "Note that the frequency grid els_f_grid must hold\n"
         "offset frequencies from line center. Hence, the\n"
         "line center frequency is not needed as input.\n"
         "\n"
         "Output:\n"
         "   els        : The lineshape function [1/Hz]\n"
         "\n"
         "Input:\n"
         "   ls_sigma   : Lorentz width [Hz].\n"
         "   ls_gamma   : Doppler width [Hz].\n"
         "   els_f_grid : Frequency grid [Hz]."
        ),
        OUTPUT( els_ ),
        INPUT(  ls_sigma_, ls_gamma_, els_f_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("GroundNoScatteringSingleEmissivity"),
        DESCRIPTION
        (
         "Treats the ground to not cause any scattering, and to have a\n"
         "reflection coefficient of 1-e. \n"
         "\n"
         "The size of *ground_emission* is set to [ nf, stokes_dim ] where \n"
         "nf is the length of *f_grid*. Columns 2-4 are set to zero.\n"
         "The temperature of the ground is obtained by interpolating \n"
         "*t_field* to the position of the ground reflection. The obtained \n"
         "temperature and *f_grid* are then used as input to the Planck\n"
         "function. The emission from the ground is then calculated as eB,\n"
         "where B is the Planck function.\n"
         "\n"
         "It is here assumed that the downwelling radiation to consider\n"
         "comes from a single direction and the returned *ground_los*\n"
         "contains only one LOS. The slope of the ground is considered\n"
         "when calculating the LOS for the downwelling radiation. The\n"
         "reflection matrices in *ground_refl_coeffs* are all set to be\n"
         "diagonal matrices, where all diagonal elements are 1-e.\n"
         "\n"
         "Keywords: \n"
         "   e : Ground emissivity. Must be a value in the range [0,1].\n"
         "       All frequencies are assumed to have the same e."
        ),
        OUTPUT( ground_emission_, ground_los_, ground_refl_coeffs_ ),
        INPUT( f_grid_, stokes_dim_, a_gp_p_, a_gp_lat_, a_gp_lon_, a_los_,
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, 
               r_geoid_,z_ground_, t_field_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS(    "e"    ),
        TYPES(    Numeric_t )));
 
  md_data_raw.push_back     
    ( MdRecord
      ( NAME("GroundTreatAsBlackbody"),
        DESCRIPTION
        (
         "Sets the ground variables (see below) to model a blackbdoy ground.\n"
         "\n"
         "The function creates the variables *ground_emission*, *ground_los*\n"
         "and *ground_refl_coeffs*. When the ground is treated to act as a\n"
         "blackbody, no downwelling radiation needs to be calculated and\n"
         "*ground_los* and *ground_refl_coeffs* are set to be empty.\n"
         "\n"
         "The size of *ground_emission* is set to [ nf, stokes_dim ] where \n"
         "nf is the length of *f_grid*. Columns 2-4 are set to zero.\n"
         "\n"
         "The temperature of the ground is obtained by interpolating \n"
         "*t_field* to the position of the ground reflection. The obtained \n"
         "temperature and *f_grid* are then used as input to the Planck\n"
         "function and the calculated blackbody radiation is put into the\n"
         "first column of *ground_emission*.\n"
         "\n"
         "Note that this function does not use *a_los*, *r_geoid* and\n"
         "*z_ground* as input, and if used inside *ground_refl_agenda*,\n"
         "ignore commands for those variables must be added to the agenda."
        ),
        OUTPUT( ground_emission_, ground_los_, ground_refl_coeffs_ ),
        INPUT( f_grid_, stokes_dim_, a_gp_p_, a_gp_lat_, a_gp_lon_, 
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, t_field_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME("lsWithCutoffAdd"),
        DESCRIPTION
        (
         "Adds lineshape with cutoff to ls.\n"
         "\n"
         "The output variable ls has to exist before and has to have the correct\n"
         "size. Use this method only if you want to have a cutoff, otherwise use\n"
         "the simpler method lsAdd. What this method does is the following:\n"
         "\n"
         "0. Check if any part of f_grid is inside the cutoff. If not, return\n"
         "   immediately.\n"
         "\n"
         "1. Create els_f_grid from f_grid. Only frequencies inside the cutoff\n"
         "   are used, but the cutoff frequency itself is also added. Note that\n"
         "   els_f_grid is relative to ls_f0.\n"
         "\n"
         "2. Use els_agenda to compute a lineshape. (This\n"
         "   should be something very simple, e.g., Lorentz or Voigt.)\n"
         "\n"
         "3. Find the generated lineshape in WSV els.\n"
         "\n"
         "4. Subtract the value at the cutoff.\n"
         "\n"
         "5. Add to ls.\n"
         "\n"
         "Output:\n"
         "   ls         : The lineshape added to previous content.\n"
         "   els        : The output of els_agenda.\n"
         "                (Comunication variable with called agenda.)\n"
         "   els_f_grid : The frequency grid for els_agenda.\n"
         "                (Comunication variable with called agenda.)\n"
         "\n"
         "Input:\n"
         "   ls         : The lineshape to add to.\n"
         "   els_agenda : Used to comute lineshape.\n"
         "   ls_cutoff  : Cutoff frequency (must be > 0).\n"
         "   ls_f0      : Line center frequency (may be < 0).\n"
         "   ls_gamma   : The pressure broadening parameter (must be > 0).\n"
         "   ls_sigma   : The Doppler broadening parameter (must be > 0).\n"
         "   f_grid     : Global frequency grid."
        ),
        OUTPUT( ls_,
                els_,
                els_f_grid_ ),
        INPUT(  ls_,
                els_agenda_,
                ls_cutoff_,
                ls_f0_, ls_gamma_, ls_sigma_, f_grid_ ),
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
         "input vector. A standard way to use the function should be (as \n"
         "part of *i_space_agenda*): \n"
         "   MatrixCBR(i_space,f_grid){} \n"
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
         "This function works as MatrixCBR but the temperature for which \n"
         "(unpolarised) balckbody radiation shall be calculated is selected \n"
         "as a keyword argument.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : Matrix with cosmic background radiation. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A set of frequencies. \n"
         "\n"
         "Keyword: \n"
         "   t : Temperature for the balckbody radiation. "
        ),
        OUTPUT(),
        INPUT( stokes_dim_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "t" ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("MatrixPrint"),
        DESCRIPTION
        (
         "Prints a matrix variable on the screen."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Matrix_ ),
        KEYWORDS( ),
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
         "Creates a workspace matrix with the specified size and sets \n"
         "all values of the matrix to the specified value. \n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Keywords:\n"
         "   nrows : The number of rows of the matrix to create. \n"
         "   ncols : The number of columns of the matrix to create. \n"
         "   value : The value of the matrix elements. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT(),
        KEYWORDS( "nrows", "ncols", "value"   ),
        TYPES(    Index_t, Index_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixSetTakingSizeFromMatrix"),
        DESCRIPTION
        (
         "Creates a workspace vector with the same size as another matrix,\n"
         "and sets all values of the new matrix to the specified value. \n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Matrix : The matrix specifying the size. \n"
         "\n"
         "Keywords:\n"
         "   value  : The value of the matrix elements. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Matrix_ ),
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
         "   Matrix : A matrix with radiance values. \n"
         "   Vector : A set of frequencies. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Matrix_, Vector_ ),
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
         "   Matrix : A matrix with radiance values. \n"
         "   Vector : A set of frequencies. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Matrix_, Vector_ ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("NumericPrint"),
        DESCRIPTION
        (
         "Prints the value of a Numeric variable on the screen."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Numeric_ ),
        KEYWORDS( ),
        TYPES( )));

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
      ( NAME( "ParticleTypeAdd" ),
        DESCRIPTION
        (
         "This method reads the amplitute matrix and the particle number\n"
         "density field from a data base. \n"
         "\n"
         "The method allows the user to chose particle types and particle \n"
         "number density fields. The methods reads from the chosen files \n"
         "and appends the variables *amp_mat_raw* and *pnd_field_raw*. \n"
         "There is one database for particle number density fields ( ....),\n"
         "which includes the following particle types:\n"
         "\n"
         "Another database (....) containes the amplitude matrices for \n"
         "those particle types from which all optical properties can be \n"
         "derived.\n"
         ),
        OUTPUT(amp_mat_raw_, pnd_field_raw_),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("filename_amp_mat", "filename_pnd_field"),
        TYPES(String_t, String_t)));

   md_data_raw.push_back
    ( MdRecord
      ( NAME( "ParticleTypeInit" ),
        DESCRIPTION
        (
         "This method initializes variables containing data about the \n"
         "optical properties of particles (*amp_mat_raw*) and about the \n"
         "particle number distribution (*pnd_field_raw*)\n"
         "\n"
         "*ParticleTypeInit* has to be executed before executing \n"
         "*ParticleTypeAdd*.\n"
         ),
        OUTPUT(amp_mat_raw_, pnd_field_raw_),
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
      ( NAME( "pha_mat_sptCalc" ),
        DESCRIPTION
        (
         "This method calculates the phase matrix for a single particle type.\n"
         "\n"
         "Phase matrix explains the tranformation of Stokes parameters of \n"
         "incident plane wave into those of the scattered spherical wave \n"
         "due to scattering of radiation by the particle. All the elements \n"
         "of phase matrix can be expressed in terms of the elements of the\n"
         "amplitude matrix.\n"
         "\n"
         "The output of the method pha_mat_sptCalc is pha_mat_spt(Tensor 5,\n"
         "size: [Npt,Nza,Naa,stokes_dim,stokes_dim]). The input to the method\n"
         "pha_mat_sptCalc are *pha_mat_spt*,*amp_mat*(Tensor 6,\n"
         "Size=[Npt,Nza,Naa,Nza,Naa,8]), *scat_za_index* and scat_aa_index.\n"
         "\n"
         "The variables *scat_za_index* and *scat_aa_index picks the right \n"
         "element of the Tensor *amp_mat*.Then this method calls the \n"
         "function amp2pha which does the actual physics, that of computing\n"
         "the elements of phase matrix from the elements of amplitude \n"
         "matrix."
         ),
        OUTPUT(pha_mat_spt_),
        INPUT(pha_mat_spt_, amp_mat_, scat_za_index_, scat_aa_index_),
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
               lon_grid_, z_field_, r_geoid_, z_ground_, 
               cloudbox_on_, cloudbox_limits_, a_pos_, a_los_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("PpathPrint"),
        DESCRIPTION
        (
         "Prints a variable of type Ppath on the screen."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Ppath_ ),
        KEYWORDS( ),
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
         "points and points of ground intersections. The keyword *lmax* \n"
         "gives the option to include additional points to ensure that the\n"
         "distance along the path between the points does not exceed the \n"
         "selected maximum length. No additional points are included if\n"
         "*lmax* is set to be negative.\n"
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
               z_field_, r_geoid_, z_ground_ ),
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
        OUTPUT( ppath_step_, a_pressure_, a_temperature_, a_vmr_list_, 
                refr_index_ ),
        INPUT( refr_index_agenda_, ppath_step_, atmosphere_dim_, p_grid_, 
               lat_grid_, lon_grid_, z_field_, t_field_, vmr_field_, r_geoid_,
               z_ground_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "lraytrace", "lmax"    ),
        TYPES(    Numeric_t,   Numeric_t )));

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
      ( NAME("RefrIndexFieldAndGradients"),
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
        OUTPUT( refr_index_, a_pressure_, a_temperature_, a_vmr_list_ ),
        INPUT( refr_index_agenda_, atmosphere_dim_, p_grid_, lat_grid_, 
               lon_grid_, r_geoid_, z_field_, t_field_, vmr_field_ ),
        GOUTPUT( Tensor4_ ),
        GINPUT( Vector_, Vector_, Vector_  ),
        KEYWORDS(),
        TYPES()));

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
        INPUT( a_pressure_, a_temperature_, a_vmr_list_, gas_species_ ),
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
         "Main function for calculation of spectra and WFs.\n"
         "\n"
         "More text will be written (PE)."
        ),
        OUTPUT( y_rte_, ppath_, ppath_step_, i_rte_, mblock_index_,
                a_pos_, a_los_, a_gp_p_, a_gp_lat_, a_gp_lon_, i_space_, 
                ground_emission_, ground_los_, ground_refl_coeffs_ ),
        INPUT( ppath_step_agenda_, rte_agenda_, i_space_agenda_,
               ground_refl_agenda_,
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_, 
               t_field_, r_geoid_, z_ground_, cloudbox_on_, cloudbox_limits_, 
               scat_i_p_, scat_i_lat_, scat_i_lon_, 
               scat_za_grid_, scat_aa_grid_,
               sensor_pos_, sensor_los_, f_grid_, stokes_dim_,
               antenna_dim_, mblock_za_grid_, mblock_aa_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "RteEmissionStd" ),
        DESCRIPTION
        (
         "Standard RTE function with emission.\n"
         "\n"
         "This function does a clearsky radiative transfer calculation for\n"
         "a given propagation path. \n"
         "Gaseous emission and absorption is calculated for each propagation\n"
         " path point using the agenda *gas_absorption_agenda*. \n"
         "Absorption vector and extinction matrix are created using \n" 
         "*opt_prop_part_agenda*.\n"
         "The coefficients for the radiative transfer are averaged between\n" 
         "two successive propagation path points. \n"
        ),
        OUTPUT( i_rte_, abs_vec_, ext_mat_, a_pressure_, a_temperature_,
                a_vmr_list_, f_index_ ),
        INPUT( i_rte_, ppath_, f_grid_, stokes_dim_, 
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, t_field_,
               vmr_field_, scalar_gas_absorption_agenda_, 
               opt_prop_gas_agenda_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "r_geoidSpherical" ),
        DESCRIPTION
        (
         "Sets the geoid to be a perfect sphere.\n"
         "\n"
         "The radius of the sphere is selected by the keyword argument *r*.\n"
         "If the keyword is set to be negative, the radius is set to the\n"
         "global internal varaible *EARTH_RADIUS*, defined in constants.cc.\n"
         "\n"
         "Keywords:\n"
         "   r : Radius of geoid sphere. See further above."
        ),
        OUTPUT( r_geoid_ ),
        INPUT( atmosphere_dim_, lat_grid_, lon_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "r" ),
        TYPES(    Numeric_t )));

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
      ( NAME( "scat_fieldCalc" ),
        DESCRIPTION
        (
         "This method calculates the scattering integral.\n"
         "\n"
         "By scattering integral we mean the field generated by integrating\n"
         "the product of intensity field and phase matrix over all incident \n"
         "angles. This term results solely due to the scattering properties\n"
         "of ice particles in cirrus clouds.  \n"
         "\n"
         "The output of this method is the scattered field *scat_field*\n"
         "which is used in the radiative transfer equation to give a new\n"
         "radiation field *i_field*. The dimensions of *scat_field* and \n"
         "*i_field* are the same. This resultant field is again given as \n"
         "input to this method, which calculates a new *scat_field*.  The \n"
         "iteration continues till the field converges.  Another important\n"
         "requirement for this method is the phase matrix.  For this we \n"
         "give as input to this method *pha_mat_spt* and *pnd_field*. From\n"
         "these two workspace variables we calculate *pha_mat* with\n"
         "the method *pha_matCalc*.  \n"
         ),
        OUTPUT( scat_field_, pha_mat_, pha_mat_spt_ ),
        INPUT( amp_mat_, i_field_, pnd_field_, scat_za_grid_, 
               scat_aa_grid_, p_grid_, lat_grid_, lon_grid_, 
               atmosphere_dim_, cloudbox_limits_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "scat_iPut" ),
        DESCRIPTION
        (
         "Method for the communication between cloudbox and clearsky.\n"
         "\n"
         "This method puts the scattered radiation field into the interface\n"
         "variables between the cloudbox and the clearsky, which are \n"
         "*scat_i_p*, *scat_i_lat* and *scat_i_lon*. As i_field is only\n"
         "stored for one frequency given by *f_index* this method has\n" 
         "to be\n"
         "executed after each scattering calculation to store the scattered\n"
         "field on the boundary of the cloudbox.\n"
         "\n"
         "The best way to calculate spectra including the influence of\n" 
         "scattering is to set up the *scat_mono_agenda* where this method \n"
         "can be included.\n"
         "\n"
         ),
        OUTPUT( scat_i_p_, scat_i_lat_, scat_i_lon_ ),
        INPUT( i_field_, f_grid_, f_index_,   p_grid_, lat_grid_, 
               lon_grid_, scat_za_grid_, scat_aa_grid_, stokes_dim_,
               atmosphere_dim_, cloudbox_limits_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
   ( MdRecord
      ( NAME( "ScatteringInit" ),
        DESCRIPTION
        (
         "Initialize variables for a scattering calculation. \n"
         "\n"
         "Variables needed in the scattering calculations are initialzed\n"
         "here. This method has to be executed before using \n"
         "*ScatteringMain*.\n"
         "\n"
         ),
        OUTPUT(scat_p_index_, scat_lat_index_, scat_lon_index_, 
               scat_za_index_, scat_aa_index_, iteration_counter_, pha_mat_,
               pha_mat_spt_, ext_mat_spt_, abs_vec_spt_, scat_field_,
               i_field_),
        INPUT(stokes_dim_, atmosphere_dim_, scat_za_grid_, scat_aa_grid_,
              cloudbox_limits_, part_types_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "ScatteringMain" ),
        DESCRIPTION
        (
         "This method executes *scat_mono_agenda* for each frequency defined\n"
         "in *f_grid* \n"
         "\n"
         "If the number of frequencies is only two, it is assumed that the \n"
         "user is only interested in a monochromatic scattering calculation\n"
         "and executes the agenda only for the first frequency.\n"
         "\n"
         ),
        OUTPUT(f_index_),
        INPUT(f_grid_, scat_mono_agenda_),
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
         "This function assumes that the first element of *a_pos* is set\n"
         "to the geometric altitude for the position of the sensor. \n"
         "The variable *a_pos* shall contain the radius instead of the\n"
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
      ( NAME("SensorIntegrationVector"),
        DESCRIPTION
        (
         "Returns a vector that approximates an integral of products.\n"
         "\n"
         "The output (row) vector multiplied with an unknown (column)\n"
         "vector approximates the integral of the product between the two\n"
         "functions corresponding to the two vectors.\n"
         "\n"
         "See ARTS User Guide (ver. 1.0.64), chapter 7 Sensor modelling,\n"
         "section 7.2 Integration as vector multiplication, for more details.\n"
         "\n"
         "Generic output: \n"
         "   Vector: the (row) vector.\n"
         "Generic input: \n"
         "   Vector: the values of the known function.\n"
         "   Vector: the grid of the known function.\n"
         "   Vector: the grid of the unknown function."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("StringPrint"),
        DESCRIPTION
        (
         "Prints the value of a String variable on the screen."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( String_ ),
        KEYWORDS( ),
        TYPES( )));

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
         "Global output: \n"
         "   Tensor3 : The scaled tensor3. \n"
         "\n"
         "Global input: \n"
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
         "Creates a workspace tensor3 with the specified size and sets \n"
         "all values of the tensor3 to the specified value. \n"
         "\n"
         "Global output: \n"
         "   Tensor3 : The tensor3 to be created. \n"
         "\n"
         "Keywords:\n"
         "   npages : The number of pages of the tensor3 to create. \n"
         "   nrows  : The number of rows of the tensor3 to create. \n"
         "   ncols  : The number of columns of the tensor3 to create. \n"
         "   value  : The value of the tensor3 elements. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor3_ ),
        GINPUT(),
        KEYWORDS( "npages", "nrows", "ncols", "value"   ),
        TYPES( Index_t, Index_t, Index_t, Numeric_t )));

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
         "Global output: \n"
         "   Tensor4 : The scaled tensor4. \n"
         "\n"
         "Global input: \n"
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
         "Creates a workspace tensor4 with the specified size and sets \n"
         "all values of the tensor4 to the specified value. \n"
         "\n"
         "Global output: \n"
         "   Tensor4 : The tensor4 to be created. \n"
         "\n"
         "Keywords:\n"
         "   nbooks : The number of books of the tensor4 to create. \n"
         "   npages : The number of pages of the tensor4 to create. \n"
         "   nrows  : The number of rows of the tensor4 to create. \n"
         "   ncols  : The number of columns of the tensor4 to create. \n"
         "   value  : The value of the tensor4 elements. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor4_ ),
        GINPUT(),
        KEYWORDS( "nbooks", "npages", "nrows", "ncols", "value"   ),
        TYPES( Index_t, Index_t, Index_t, Index_t, Numeric_t )));

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
         "Global output: \n"
         "   Tensor5 : The scaled tensor5. \n"
         "\n"
         "Global input: \n"
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
         "Creates a workspace tensor5 with the specified size and sets \n"
         "all values of the tensor5 to the specified value. \n"
         "\n"
         "Global output: \n"
         "   Tensor5 : The tensor5 to be created. \n"
         "\n"
         "Keywords:\n"
         "   nshelfs : The number of shelfs of the tensor5 to create. \n"
         "   nbooks  : The number of books of the tensor5 to create. \n"
         "   npages  : The number of pages of the tensor5 to create. \n"
         "   nrows   : The number of rows of the tensor5 to create. \n"
         "   ncols   : The number of columns of the tensor5 to create. \n"
         "   value   : The value of the tensor5 elements. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor5_ ),
        GINPUT(),
        KEYWORDS( "nshelfs", "nbooks", "npages", "nrows", "ncols", "value" ),
        TYPES( Index_t, Index_t, Index_t, Index_t, Index_t, Numeric_t )));

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
         "Global output: \n"
         "   Tensor6 : The scaled tensor6. \n"
         "\n"
         "Global input: \n"
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
         "Creates a workspace tensor6 with the specified size and sets \n"
         "all values of the tensor6 to the specified value. \n"
         "\n"
         "Global output: \n"
         "   Tensor6 : The tensor6 to be created. \n"
         "\n"
         "Keywords:\n"
         "   nvitrines : The number of vitrines of the tensor6 to create. \n"
         "   nshelfs   : The number of shelfs of the tensor6 to create. \n"
         "   nbooks    : The number of books of the tensor6 to create. \n"
         "   npages    : The number of pages of the tensor6 to create. \n"
         "   nrows     : The number of rows of the tensor6 to create. \n"
         "   ncols     : The number of columns of the tensor6 to create. \n"
         "   value     : The value of the tensor6 elements. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor6_ ),
        GINPUT(),
        KEYWORDS( "nvitrines", "nshelfs", "nbooks", "npages", "nrows",
                  "ncols", "value" ),
        TYPES( Index_t, Index_t, Index_t, Index_t, Index_t, Index_t,
               Numeric_t )));

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
         "Global output: \n"
         "   Tensor7 : The scaled tensor7. \n"
         "\n"
         "Global input: \n"
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
         "Creates a workspace tensor7 with the specified size and sets \n"
         "all values of the tensor7 to the specified value. \n"
         "\n"
         "Global output: \n"
         "   Tensor7 : The tensor7 to be created. \n"
         "\n"
         "Keywords:\n"
         "   nlibraries : The number of libraries of the tensor7 to create. \n"
         "   nvitrines  : The number of vitrines of the tensor7 to create. \n"
         "   nshelfs    : The number of shelfs of the tensor7 to create. \n"
         "   nbooks     : The number of books of the tensor7 to create. \n"
         "   npages     : The number of pages of the tensor7 to create. \n"
         "   nrows      : The number of rows of the tensor7 to create. \n"
         "   ncols      : The number of columns of the tensor7 to create. \n"
         "   value      : The value of the tensor7 elements. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor7_ ),
        GINPUT(),
        KEYWORDS( "nlibraries", "nvitrines", "nshelfs", "nbooks", "npages",
                  "nrows", "ncols", "value" ),
        TYPES( Index_t, Index_t, Index_t, Index_t, Index_t, Index_t,
               Index_t, Numeric_t )));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor6WriteIteration"),
        DESCRIPTION
        (
       "Write iterated fields.\n"
       "\n"
       "This function writes intermediate resultes, the iterations of \n"
       "fields to xml files. It can be used to check the solution method \n"
       "for the RTE with scattering integral, which is an iterative \n"
       "numerical method. It is useful to look how the radiation field \n"
       "*i_field* and the scattered field *scat_field* behave. \n"
       "\n"
       "The user can give an array containing the iterations which shall \n"
       "be written to files as a keyword to the method. E.g. if \n"
       "'iterations = [3, 6, 9]', the 3rd, 6th and 9th iterations are \n"
       "stored in the files 'iteration_field_3.xml', \n"
       "'iteration_field_6.xml' ...\n"
       "\n"
       "If you want to save all the iterations the array has to contain \n"
       "just one element set to 0: 'iterations = [0]'.\n"
       "\n"
       "Note: The workspace variable iteration_counter has to be set as 0 \n"
       "in the control file before using this method.\n"
       "\n"
        ),
        OUTPUT( ),
        INPUT(iteration_counter_),
        GOUTPUT( ),
        GINPUT(Tensor6_),
        KEYWORDS("iterations"),
        TYPES(Array_Index_t )));

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
      ( NAME("VectorPrint"),
        DESCRIPTION
        (
         "Prints a vector variable on the screen."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
        TYPES( )));

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
         "Creates a workspace vector with the specified length and sets \n"
         "all values of the vector to the specified value. \n"
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Keywords:\n"
         "   length : The length of the new vector. \n"
         "   value  : The value of the vector elements. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "length", "value"   ),
        TYPES(    Index_t,    Numeric_t )));

   md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorSetElement"),
        DESCRIPTION
        (
         "Sets the selected vector element to the given value.\n"
         "\n"
         "The vector must be initiated before calling the function. An error\n"
         "is issued if a position outside the range of the vector is\n"
         "selected.\n"
         "\n"
         "Note that the indexing is zero based. That is, the first element \n"
         "has index 0.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The vector to be modified. \n"
         "\n"
         "Keywords:\n"
         "   pos   : Column of the position element to set. \n"
         "   value : The value of the vector element. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "pos",   "value"   ),
        TYPES(    Index_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorSetTakingLengthFromVector"),
        DESCRIPTION
        (
         "Creates a workspace vector with the same length as another vector,\n"
         "and sets all values of the new vector to the specified value. \n"
         "\n"
         "A possible usage of the function is: \n"
         "  VectorSetLengthFromVector(vector1,f_grid){value=0.75} \n"
         "where *vector1* then can be used to set *e_ground*. \n"
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector specifying the length.. \n"
         "\n"
         "Keywords:\n"
         "   value  : The value of the vector elements. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "value"   ),
        TYPES(    Numeric_t )));

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
         "   Vector : A vector with radiance values. \n"
         "   Vector : A set of frequencies. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_, Vector_ ),
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
         "The conversion requieres that the frequency of each intensity \n"
         "value is known. The last general input vector is assumed to be a \n"
         "set of frequencies. It is further assumed that the frequencies of \n"
         "the first general input vector are obtained by repeating the last \n"
         "vector an integer number of times. \n"
         "\n"
         "If *y* shall be converted from radiances to brightness \n"
         "temperatures and no instrument response has been applied, the \n"
         "conversion is done as: \n"
         "   VectorToTbByRJ(y,y,f_grid){} \n"
         "\n"
         "Generic output: \n"
         "   Vector : A vector with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A vector with radiance values. \n"
         "   Vector : A set of frequencies. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_, Vector_ ),
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
      ( NAME( "yNoPolarisation" ),
        DESCRIPTION
        (
         "Converts *y_rte* to *y* assuming an unpolarised instrument.\n"
         "\n"
         "This function assumes that the instrument is equally sensitive for\n"
         "all polarisations. This corresponds to that *y* simply equals the \n"
         "first column of *y_rte*."
        ),
        OUTPUT( y_ ),
        INPUT( y_rte_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

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
               lon_grid_, z_field_, r_geoid_, z_ground_ ),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "z_recieve", "z_send", "t_sample", 
                  "z_scan_low", "z_scan_high" ),
        TYPES( Numeric_t, Numeric_t, Numeric_t,
               Numeric_t, Numeric_t )));








  //
  // Below this line you find methods not touched for ARTS 2. 
  // Please revise the documentation etc. before a methods is moved up.
  // Place functions in alphabetical order. 
  // Finally, all methods below the line will be removed.
  //--------------------------------------------------------------------

//======================================================================
//=== IO methods
//======================================================================

//=== Index ============================================================

  md_data_raw.push_back
    ( MdRecord
      ( NAME("IndexWriteAscii"),
        DESCRIPTION(
                    "Writes an index value to an ASCII file.\n"
                    "\n"
                    "The index value of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the index is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Index : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Index_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("IndexReadAscii"),
        DESCRIPTION(
                    "Reads a index value from an ASCII file.\n"
                    "\n"
                    "The index value is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the index is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   Index : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Index_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("IndexWriteBinary"),
        DESCRIPTION(
                    "Writes an index to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Index : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename: Name of the output file.\n"
                    ),
   
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Index_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("IndexReadBinary"),
        DESCRIPTION(
                    "Reads an index from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   Index : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Index_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT

//=== NUMERIC ==========================================================

  md_data_raw.push_back
    ( MdRecord
      ( NAME("NumericWriteAscii"),
        DESCRIPTION(
                    "Writes a numeric value to an ASCII file.\n"
                    "\n"
                    "The numeric value of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the numeric is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Numeric : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Numeric_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("NumericReadAscii"),
        DESCRIPTION(
                    "Reads a numeric value from an ASCII file.\n"
                    "\n"
                    "The numeric value is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the numeric is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   Numeric : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Numeric_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("NumericWriteBinary"),
        DESCRIPTION(
                    "Writes a numeric value to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Numeric : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Numeric_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("NumericReadBinary"),
        DESCRIPTION(
                    "Reads a numeric from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   Numeric : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Numeric_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT


//=== Vector ==========================================================



  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorWriteAscii"),
        DESCRIPTION(
                    "Writes a vector to an ASCII file.\n"
                    "\n"
                    "The vector of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the vector is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Vector : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Vector_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorReadAscii"),
        DESCRIPTION(
                    "Reads a vector from an ASCII file.\n"
                    "\n"
                    "The vector is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the vector is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   Vector : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorWriteBinary"),
        DESCRIPTION(
                    "Writes a vector to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Vector : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Vector_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorReadBinary"),
        DESCRIPTION(
                    "Reads a vector from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n" 
                    "Generic output: \n"
                    "   Vector : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT




//=== Matrix ==========================================================


  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixWriteAscii"),
        DESCRIPTION(
                    "Writes a matrix to an ASCII file.\n"
                    "\n"
                    "The matrix of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the matrix is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Matrix : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Matrix_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixReadAscii"),
        DESCRIPTION(
                    "Reads a matrix from an ASCII file.\n"
                    "\n"
                    "The matrix is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the matrix is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   Matrix : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixWriteBinary"),
        DESCRIPTION(
                    "Writes a matrix to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Matrix : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Matrix_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixReadBinary"),
        DESCRIPTION(
                    "Reads a matrix from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   Matrix : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT





//=== ArrayOfIndex =====================================================

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexWriteAscii"),
        DESCRIPTION(
                    "Writes a index array to an ASCII file.\n"
                    "\n"
                    "The index array of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the index array is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfIndex : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfIndex_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexReadAscii"),
        DESCRIPTION(
                    "Reads a index array from an ASCII file.\n"
                    "\n"
                    "The index array is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the index array is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   ArrayOfIndex : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"                    
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfIndex_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexWriteBinary"),
        DESCRIPTION(
                    "Writes an index array to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfIndex : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfIndex_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexReadBinary"),
        DESCRIPTION(
                    "Reads an index array from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   ArrayOfIndex : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfIndex_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT


//=== ArrayOfVector ====================================================

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfVectorWriteAscii"),
        DESCRIPTION(
                    "Writes an array of vectors to an ASCII file.\n"
                    "\n"
                    "The array of vectors of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the array of vectors is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfVector : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfVector_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfVectorReadAscii"),
        DESCRIPTION(
                    "Reads an array of vectors from an ASCII file.\n"
                    "\n"
                    "The array of vectors is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the array of vectors is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   ArrayOfVector : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfVector_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfVectorWriteBinary"),
        DESCRIPTION(
                    "Writes a vector array to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfVector : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfVector_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfVectorReadBinary"),
        DESCRIPTION(
                    "Reads a vector array from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output:  \n"
                    "   ArrayOfVector : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfVector_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT



//=== ArrayOfMatrix ====================================================

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixWriteAscii"),
        DESCRIPTION(
                    "Writes an array of matrices to an ASCII file.\n"
                    "\n"
                    "The array of matrices of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the array of matrices is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "The format is as follows:\n"
                    "\n"
                    "# <comments>\n"
                    "<n_array_elements>\n"
                    "<n_rows> <n_columns>\n"
                    "<elements>\n"
                    "<n_rows> <n_columns>\n"
                    "<elements>\n"
                    "...\n"
                    "\n"
                    "Example:\n"
                    "# Generated by arts-0.0.16, Apr 29 2000, 17:38:44\n"
                    "2\n"
                    "3 4\n"
                    "xx xx xx xx\n"
                    "xx xx xx xx\n"
                    "xx xx xx xx\n"
                    "2 2\n"
                    "yy yy\n"
                    "yy yy"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfMatrix : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfMatrix_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixReadAscii"),
        DESCRIPTION(
                    "Reads an array of matrices from an ASCII file.\n"
                    "\n"
                    "The array of matrices is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the array of matrices is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   ArrayOfMatrix : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfMatrix_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixWriteBinary"),
        DESCRIPTION(
                    "Writes a matrix array to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfMatrix : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfMatrix_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixReadBinary"),
        DESCRIPTION(
                    "Reads a matrix array from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   ArrayOfMatrix : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfMatrix_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT

//=== STRING ============================================================


  md_data_raw.push_back
    ( MdRecord
      ( NAME("StringWriteAscii"),
        DESCRIPTION(
                    "Writes a string to an ASCII file.\n"
                    "\n"
                    "The string of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the string is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfStringWriteAscii* for file format.\n"
                    "\n"
                    "Generic input: \n"
                    "   String : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( String_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("StringReadAscii"),
        DESCRIPTION(
                    "Reads a string from an ASCII file.\n"
                    "\n"
                    "The string is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the string is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfStringWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   String : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( String_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("StringWriteBinary"),
        DESCRIPTION(
                    "Writes a String to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   String : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( String_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("StringReadBinary"),
        DESCRIPTION(
                    "Reads a String from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   String : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( String_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT


//=== ArrayOfSTRING =========================================================


  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfStringWriteAscii"),
        DESCRIPTION(
                    "Writes an array of strings to an ASCII file.\n"
                    "\n"
                    "The array of strings of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the string is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "The format is as follows:\n"
                    "\n"
                    "# <comments>\n"
                    "<n_Strings>\n"
                    "<String 1>\n"
                    "<String 2>\n"
                    "...\n"
                    "\n"
                    "Example:\n"
                    "# Generated by arts-0.0.16, Apr 29 2000, 17:38:44\n"
                    "2\n"
                    "A String\n"
                    "Another String\n"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfString : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfString_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfStringReadAscii"),
        DESCRIPTION(
                    "Reads an array of strings from an ASCII file.\n"
                    "\n"
                    "The array of strings is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the string is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfStringWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   ArrayOfString : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfString_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfStringWriteBinary"),
        DESCRIPTION(
                    "Writes a String array to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfString : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfString_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfStringReadBinary"),
        DESCRIPTION(
                    "Reads a String array from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   ArrayOfString : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfString_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT



//======================================================================
//=== Absorption methods
//======================================================================

//=== Spectroscopic methods ============================================

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgSetEmpty"),
//      DESCRIPTION
//      (
//       "Sets lines_per_tg to empty line lists.\n"
//       "\n"
//       "You can use this method to set lines per tag if you do not reall want\n"
//       "to compute line spectra. Formally, absCalc will still require\n"
//       "lines_per_tg to be set.\n"
//       ),
//      OUTPUT(   lines_per_tg_      ),
//      INPUT(    tgs_        ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(  ),
//      TYPES(    )));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgReadFromCatalogues"),
//      DESCRIPTION(
//                  "This method can read lines from different line \n"
//                  "catalogues.\n"
//                  "\n"
//                  "For each tag group, you can specify which catalogue\n"
//                  "to use. Because the method creates lines_per_tg directly,\n"
//                  "it replaces for example thefollowing two method calls:\n"
//                  "  - linesReadFromHitran\n"
//                  "  - lines_per_tgCreateFromLines\n"
//                  "   This method needs as input WSVs the list of tag \n"
//                  "groups. Keyword parameters must specify the names of\n"
//                  "the catalogue files to use and the matching formats.\n"
//                  "Names can be anything, formats can currently be \n"
//                  "HITRAN96, MYTRAN2, JPL, or ARTS. Furthermore, keyword\n"
//                  "parameters have to specify minimum and maximum \n"
//                  "frequency for each tag group. To safe typing, if there\n"
//                  "are less elements in the keyword parameters than there\n"
//                  "are tag groups, the last parameters are applied to all\n"
//                  "following tag groups.\n"
//                  "\n"
//                  "Example usage:\n"
//                  "\n"
//                  "lines_per_tgReadFromCatalogues{\n"
//                  "  filenames = [ \"../data/cat1.dat\", \"../data/cat2.dat\" ]\n"
//                  "  formats   = [ \"MYTRAN2\",          \"HITRAN96\"         ]\n"
//                  "  fmin      = [ 0,                  0                  ]\n"
//                  "  fmax      = [ 2000e9,             100e9              ]\n"
//                  "}\n"
//                  "   In this example, lines for the first tag group will\n"
//                  "be taken from cat1, lines for all other tag groups \n"
//                  "will be taken from cat2.\n"
//                  "   This methods allows you for example to use a \n"
//                  "special line file just for water vapor lines. This\n"
//                  "could be the  improved water vapor line file \n"
//                  "generated by Thomas Kuhn.\n"
//                  "   Catalogues are only read once, even if several tag\n"
//                  "groups have the same catalogue. However, in that case\n"
//                  "the frequency ranges MUST be the same. (If you want \n"
//                  "to do fine-tuning of the frequency ranges, you can do \n"
//                  "this inside the tag definitions, e.g., \"H2O-*-0-2000e9\".)\n"
//                  "   This function uses the various reading routines\n"
//                  "(linesReadFromHitran, etc.), as well as\n"
//                  "lines_per_tgCreateFromLines.\n"
//                  "\n"
//                  "Keywords: \n"
//                  "   filenames = Name (and path) of the catalogue files.\n"
//                  "   formats   = allowed formats are HITRAN96,MYTRAN2,JPL,ARTS \n"
//                  "   fmin      = Minimum frequency for lines to read in Hz.\n"
//                  "   fmax      = Maximum frequency for lines to read in Hz.\n"),
//      OUTPUT(   lines_per_tg_      ),
//      INPUT(    tgs_        ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "filenames",    "formats",      "fmin",   "fmax" ),
//      TYPES(    Array_String_t, Array_String_t, Vector_t, Vector_t)));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesReadFromHitran"),
//      DESCRIPTION(
//                  "Read all the lines from a HITRAN catalogue file in the \n"
//                  "given frequency range. Otherwise a runtime error will be\n"
//                  "thrown\n"
//                  "\n"
//                  "Please note that all lines must correspond\n"
//                  "to the legal species / isotope combinations\n"
//                  "\n"
//                  "Keywords: \n"
//                  "   filename = Name (and path) of the catalogue file.\n"
//                  "   fmin     = Minimum frequency for lines to read in Hz.\n"
//                  "   fmax     = Maximum frequency for lines to read in Hz."),
//      OUTPUT(   lines_   ),
//      INPUT(),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "filename",  "fmin",    "fmax"),
//      TYPES(    String_t,    Numeric_t, Numeric_t)));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesReadFromMytran2"),
//      DESCRIPTION(
//                  "Read all the lines from a MYTRAN2 catalogue file in the \n"
//                  "given frequency range. Otherwise a runtime error will be\n"
//                  "thrown\n"
//                  "\n"
//                  "Please note that all lines must correspond\n"
//                  "to the legal species / isotope combinations\n"
//                  "\n"
//                  "Keywords: \n"
//                  "   filename = Name (and path) of the catalogue file.\n"
//                  "   fmin     = Minimum frequency for lines to read in Hz.\n"
//                  "   fmax     = Maximum frequency for lines to read in Hz."),
//      OUTPUT(   lines_   ),
//      INPUT(),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "filename",  "fmin",    "fmax"),
//      TYPES(    String_t,    Numeric_t, Numeric_t)));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesReadFromJpl"),
//      DESCRIPTION(
//                  "Read all the lines from a JPL catalogue file in the \n"
//                  "given frequency range. Otherwise a runtime error will be\n"
//                  "thrown\n"
//                  "\n"
//                  "Please note that all lines must correspond\n"
//                  "to the legal species / isotope combinations.\n"
//                  "\n"
//                  "Keywords: \n"
//                  "   filename = Name (and path) of the catalogue file.\n"
//                  "   fmin     = Minimum frequency for lines to read in Hz.\n"
//                  "   fmax     = Maximum frequency for lines to read in Hz."),
//      OUTPUT(   lines_   ),
//      INPUT(),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "filename",  "fmin",    "fmax"),
//      TYPES(    String_t,    Numeric_t, Numeric_t)));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesReadFromArts"),
//      DESCRIPTION(
//                  "Read all the lines from an Arts catalogue file in the \n"
//                  "given frequency range. Otherwise a runtime error will be\n"
//                  "thrown \n"
//                  "\n"
//                  "Please note that all lines must correspond\n"
//                  "to the legal species / isotope combinations\n"
//                  "\n"
//                  "Keywords: \n"
//                  "   filename = Name (and path) of the catalogue file.\n"
//                  "   fmin     = Minimum frequency for lines to read in Hz.\n"
//                  "   fmax     = Maximum frequency for lines to read in Hz."),
//      OUTPUT(   lines_   ),
//      INPUT(),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "filename",  "fmin",    "fmax"),
//      TYPES(    String_t,    Numeric_t, Numeric_t)));
  
//   // FIXME: Remove this one.
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesElowToJoule"),
//      DESCRIPTION(
//                  "Just a little helper to convert the lower state energy from cm^-1\n"
//                  "(ARTSCAT-2) to Joule (ARTSCAT-3). This should be removed soon\n"),
//      OUTPUT(   lines_   ),
//      INPUT(),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( ),
//      TYPES(    )));
      
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgCreateFromLines"),
//      DESCRIPTION(
//                  "Split lines up into the different tag groups.\n"
//                  "\n"
//                  "The tag groups are tested in the order in which they are\n" 
//                  "specified in the controlfile. The lines are assigned to \n"
//                  "the tag groups in the order as the groups  are specified.\n"
//                  "That means if you do [\"O3-666\",\"O3\"],the last group O3 \n"
//                  "gets assigned all the O3 lines that do not fit in the first group."),
//      OUTPUT(   lines_per_tg_      ),
//      INPUT(    lines_, tgs_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgAddMirrorLines"),
//      DESCRIPTION(
//                  "Adds mirror lines at negative frequencies to the *lines_per_tg*.\n"
//                  "\n"
//                  "For each line at frequency +f in *lines_per_tg* a corresponding\n"
//                  "entry at frequency -f is added to *lines_per_tg*.The mirror \n"
//                  "lines are appended to the line lists after the original lines."),
//      OUTPUT(   lines_per_tg_      ),
//      INPUT(    lines_per_tg_      ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgCompact"),
//      DESCRIPTION(
//                  "Removes all lines outside the defined lineshape cutoff frequency\n"
//                  "from the *lines_per_tg*. This can save computation time.\n"
//                  "It should be particularly useful to call this method after\n"
//                  "*lines_per_tgAddMirrorLines*."),
//      OUTPUT(   lines_per_tg_      ),
//      INPUT(    lines_per_tg_, lineshape_, f_mono_  ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesWriteAscii"),
//      DESCRIPTION(
//                     "Writes the workspace variable *lines* to an ASCII file.\n"
//                     "\n"
//                  "The content of the workspace variable 'lines`\n"
//                  "The content of the workspace variable *lines*\n"
//                  "is written in ARTS line format to the file with\n"
//                     "the specified name. If the filename is omitted, the\n"
//                     "lines are written to <basename>.lines.aa.\n"
//                     "\n"
//                     "Keywords: \n"
//                     "   filename : Name of the output file.\n"
//                     ), 
//      OUTPUT(),
//      INPUT( lines_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "filename" ),
//      TYPES(    String_t   )));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgWriteAscii"),
//      DESCRIPTION(
//                     "Writes the workspace variable *lines_per_tg* to an ASCII file.\n"
//                     "\n"
//                     "The content of the workspace variable *lines_per_tg*\n"
//                     "is written in ARTS line format to the file with\n"
//                     "the specified name. If the filename is omitted, the\n"
//                     "lines are written to <basename>.lines_per_tg.aa.\n"
//                     "\n"
//                     "The array dimension is handled in a similar way as by the\n"
//                     "array of vector and matrix output functions:\n"
//                     "First an integer stating the number of tag groups.\n"
//                     "Then an integer specifying the number of lines for the\n"
//                     "first group. Then the other groups in similar fashion."
//                     "\n"
//                     "Keywords: \n"
//                     "   filename : Name of the output file.\n"
//                     ),
//      OUTPUT(),
//      INPUT( lines_per_tg_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "filename" ),
//      TYPES(    String_t   )));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("tgsDefineAllInScenario"),
//      DESCRIPTION
//      (
//       "Define one tag group for each species known to ARTS and included in an\n"
//       "atmospheric scenario.\n"
//       "\n"
//       "You can use this as an alternative to tgsDefine if you want to make an\n"
//       "absorption calculation that is as complete as possible. The method\n"
//       "goes through all defined species and tries to open the VMR file. If\n"
//       "this works the tag is included, otherwise it is skipped.\n"
//       "\n"
//       "Keywords:\n"
//       "   basename : The name and path of a particular atmospheric scenario.\n"
//       "              For example: /pool/lookup2/arts-data/atmosphere/fascod/tropical"
//            ),
//      OUTPUT( tgs_ ),
//      INPUT(),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "basename" ),
//      TYPES(    String_t   )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lineshapeDefine"),
//      DESCRIPTION(
//           "Sets the lineshape for all calculated lines.\n\n"
//           "\n"
//           "   A general lineshape profile is specified, according to a given  \n"
//           "approximation. Alongside a normalization factor is to be set - a  \n"
//           "multiplicative forefactor through which the profile can be \n"
//           "modified. This factor is just the 0th or 1st, or 2nd power of the \n"
//           "ratio between the frequency of calculation f and the center frequency\n"
//           "for a specific line f0. A cutoff frequency must also be specified in\n"
//           "order to restrict the calculation within a desired frequency region or\n"
//           "not, when there's no such region.\n"
//           "   The general lineshape profile is given by the keyword shape,\n"
//           "while the normalization factor and the cutoff frequency by\n"
//           "normalizationfactor and cutoff respectively.\n"
//           "\n"
//           "   The available values for these keywords are given below.\n"
//           "shape - \"no_shape\" : no specified shape\n"
//           "        \"Doppler\" : Doppler lineshape\n"
//           "        \"Lorentz\" : Lorentz lineshape\n"
//           "        \"Voigt_Kuntz3\" : Kuntz approximation to the Voigt profile,\n"
//           "                         accuracy > 2x10^(-3)\n"
//           "        \"Voigt_Kuntz4\" : Kuntz approximation to the Voigt profile,\n"
//           "                         accuracy > 2x10^(-4)\n"
//           "        \"Voigt_Kuntz6\" : Kuntz approximation to the Voigt profile,\n"
//           "                         accuracy > 2x10^(-6)\n"   
//           "        \"Voigt_Drayson\" : Drayson approximation to the Voigt profile \n"
//           "        \"Rosenkranz_Voigt_Drayson\" : Rosenkrantz oxygen absortion with overlap correction\n" 
//           "                                     on the basis of Drayson routine\n"                                    
//           "        \"Rosenkranz_Voigt_Kuntz6\" : Rosenkrantz oxygen absortion with overlap correction\n"
//           "                                    on the basis of Kuntz routine, accuracy > 2x10^(-6)\n"
//           "normalizationfactor - \"no_norm\": 1\n"
//           "                      \"linear\": f/f0\n" 
//           "                      \"quadratic\": (f/f0)^2.\n"
//           "cutoff - \" -1\" : no cutoff\n"
//           "         \"Number\": positive cutoff frequency in Hz.\n"
//           "\n"
//           "Example usage:\n"
//        "shape=[\"Lorentz\"]\n"
//           "normalizationfactor=[\"linear\"]\n"
//           "cutoff= [650e9]"
//           "\n"
//           "Keywords:\n"
//           "   shape               : The general profile according to an approximation.\n"
//           "   normalizationfactor : The multiplicative forefactor for the general profile.\n"
//           "   cutoff              : The frequency at which a cutoff can be made.\n"),
//      OUTPUT( lineshape_ ),
//      INPUT( tgs_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(  "shape",    "normalizationfactor",  "cutoff" ),
//      TYPES(     String_t,        String_t,         Numeric_t )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lineshape_per_tgDefine"),
//      DESCRIPTION(
//           "Sets the lineshape per tag group for all calculated lines.\n\n"
//           "\n" 
//           "   A general lineshape profile is specified, according to a given  \n"
//           "approximation for each tag group. Alongside a normalization factor\n" 
//           "is to be set also for each tag group - a multiplicative forefactor through\n"
//           "which the profile can be modified. This factor is just the 0th or 1st,\n"
//           "or 2nd power of the ratio between the frequency of calculation f and\n"
//           "the center frequency for a specific line f0. A cutoff frequency must also be\n"
//           "specified for each of the tags in  order to restrict the calculation within\n" 
//           "a desired region or not, when there's no such region.\n"
//           "   The general lineshape profile is given by the keyword shape,\n"
//           "while the normalization factor and the cutoff frequency by\n"
//           "normalizationfactor and cutoff respectively.\n"
//           "\n"
//           "   The available values for these keywords are given below.\n"
//           "shape - \"no_shape\" : no specified shape\n"
//           "        \"Doppler\" : Doppler lineshape\n"
//           "        \"Lorentz\" : Lorentz lineshape\n"
//           "        \"Voigt_Kuntz3\" : Kuntz approximation to the Voigt profile,\n"
//           "                        accuracy > 2x10^(-3)\n"
//           "        \"Voigt_Kuntz4\" : Kuntz approximation to the Voigt profile,\n"
//           "                         accuracy > 2x10^(-4)\n"
//           "        \"Voigt_Kuntz6\" : Kuntz approximation to the Voigt profile,\n"
//           "                         accuracy > 2x10^(-6)\n"   
//           "        \"Voigt_Drayson\" : Drayson approximation to the Voigt profile \n"
//           "        \"Rosenkranz_Voigt_Drayson\" : Rosenkrantz oxygen absortion with overlap correction\n" 
//           "                                     on the basis of Drayson routine\n"                                    
//           "        \"Rosenkranz_Voigt_Kuntz6\" : Rosenkrantz oxygen absortion with overlap correction\n"
//           "                                    on the basis of Kuntz routine, accuracy > 2x10^(-6)\n"
//           "normalizationfactor - \"no_norm\": 1\n"
//           "                      \"linear\": f/f0\n" 
//           "                      \"quadratic\": (f/f0)^2.\n"
//           "cutoff - \" -1\" : no cutoff\n"
//           "           \"Number\": positive cutoff frequency in Hz.\n"
//           "\n"
//           "Example usage:\n"
//        "shape = [\"Lorentz\",\"Voigt_Kuntz6\"] \n"
//        "normalizationfactor= [\"linear\", \"quadratic\"] \n"
//        "cutoff = [ 650e9, -1 ]"
//           "\n"
//           "Keywords:\n"
//           "   shape               : The general profile according to an approximation.\n"
//           "   normalizationfactor : The multiplicative forefactor for the general profile.\n"
//           "   cutoff              : The frequency at which a cutoff can be made.\n"),
//      OUTPUT( lineshape_ ),
//      INPUT( tgs_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(  "shape",           "normalizationfactor",    "cutoff" ),
//      TYPES(   Array_String_t,         Array_String_t,        Vector_t )));


// //=== Continuum methods ============================================

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("cont_descriptionInit"),
//      DESCRIPTION
//      (
//       "Initializes the two workspace variables for the continuum description,\n"
//       "*cont_description_names* and *cont_description_parameters*.\n"
//       " \n"
//       "This method does not really do anything, except setting the two\n"
//       "variables to empty Arrays. It is just necessary because the method\n"
//       "*cont_descriptionAppend* wants to append to the variables.\n"
//       "   Formally, the continuum description workspace variables are required\n"
//       "by the absorption calculation methods (e.g., *absCalc*). Therefore you\n"
//       "always have to call at least *cont_descriptionInit*, even if you do\n"
//       "not want to use any continua."
//       ),
//      OUTPUT( cont_description_names_, 
//                 cont_description_models_,
//                 cont_description_parameters_ ),
//      INPUT(),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("cont_descriptionAppend"),
//      DESCRIPTION
//      (
//       "Appends the description of a continuum model or a complete absorption\n"
//       "model to *cont_description_names* and *cont_description_parameters*.\n"
//       "\n"
//       "See online documentation for *cont_description_names* for a list of\n"
//       "allowed models and for information what parameters they require. See\n"
//       "file cont.arts in the doc/examples directory for usage examples and\n"
//       "default parameters for the various models. \n"
//       "\n"
//       "Keywords:\n"
//       "   name       : The name of a continuum model. Must match one of the models\n"
//       "                implemented in ARTS. \n"
//          "   option     : give here the option of this continuum/full model.\n"
//       "   parameters : A Vector containing the required number of parameters\n"
//       "                for the model given. The meaning of the parameters and\n"
//       "                how many parameters are required depends on the model.\n"
//       ),
//      OUTPUT( cont_description_names_, 
//                 cont_description_models_,
//                 cont_description_parameters_ ),
//      INPUT(  cont_description_names_, 
//                 cont_description_models_,
//                 cont_description_parameters_),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "tagname",  "model",   "userparameters" ),
//      TYPES(    String_t,   String_t,   Vector_t         )));


//=== Input Atmosphere methods ===========================================

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("raw_vmrsReadFromFiles"),
//         DESCRIPTION(
//           "Reads the individual VMR profile for each TAGS from file.\n"
//           "\n"
//           "Using this function one can read VMRs of specific TAGS from\n"
//           "explicitly specified files and the remaing from a scenario.\n"
//           "The filenames and the base name of atmospheric scenario\n"
//           "should be specified as keywords. One file name must\n"
//           "be specified for each tag group(each element of *tgs*).\n"
//           "The name may include a path.\n"
//        "\n"
//        "Keywords:\n"
//        "   seltags   : Must be a sub group of tags which should be read from files.\n"
//        "   filenames : Names of the files containing VMR profiles of seltags.\n"
//        "   basename  : The name of a particular atmospheric scenario.\n"
//        "               See *raw_vmrsReadFromScenario* for details. Remaining\n"
//        "               VMRs will be read from the scenario.\n"
//        "\n"
//           ),
//         OUTPUT(   raw_vmrs_         ),
//         INPUT(    tgs_                 ),
//         GOUTPUT(                       ),
//         GINPUT(                        ),
//         KEYWORDS( "seltags",       "filenames",    "basename"),
//         TYPES(    Array_String_t,  Array_String_t, String_t)));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("raw_vmrsReadFromScenario"),
//      DESCRIPTION(
//        "Reads the individual VMR profile for each tag group from a standard\n"
//        "atmospheric scenario.\n" 
//        "\n"
//           "Five different atmospheric scenarios are available in arts data:\n"
//           "tropical, midlatitude-summer, midlatitude-winter, subartic-summer\n"
//           "and subartic-winter.\n"
//        "\n"
//        "   Files in the scenarios look like this: tropical.H2O.aa\n"
//        "\n"
//        "   The basename must include the path, i.e., the files can be anywhere,\n"
//        "but they must be all in the same directory.\n"
//        "   The profile is chosen by the species name. If you have more than one\n"
//        "tag group for the same species, the same profile will be used.\n"
//        "\n"
//        "Keywords:\n"
//        "   basename :The name and path of a particular atmospheric scenario.\n"
//        "   For example:\n"
//        "   /pool/lookup2/arts-data/atmosphere/fascod/tropical\n"
//        "\n"
//        ),
//      OUTPUT(   raw_vmrs_    ),
//      INPUT(    tgs_                 ),
//      GOUTPUT(                       ),
//      GINPUT(                        ),
//      KEYWORDS( "basename"           ),
//      TYPES(    String_t             )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("AtmFromRaw"),
//      DESCRIPTION(
//        "Interpolates temperature, altitude, and VMRs to the pressure grid\n"
//        "given by p_abs.\n" 
//        "\n"
//           "The altitude is not used by the absorption routines,\n"
//        "but later on by the RT routines.\n"
//        "\n"
//        "Interpolations used: \n"
//        "\n"
//        "Temperature      : Linear interpolation in ln(p)\n"
//        "Altitude         : Linear interpolation in ln(p)\n"
//        "VMRs             : Linear interpolation in ln(p)\n"
//        "Cloud Parameters : Linear interpolation in ln(p)\n"
//        "\n"
//        ),
//      OUTPUT(   t_abs_    , z_abs_   , vmrs_           ),
//      INPUT(    tgs_, p_abs_    , raw_ptz_ , raw_vmrs_ ),
//      GOUTPUT(                                         ),         
//      GINPUT(                                          ),
//      KEYWORDS(                                        ),
//      TYPES(                                           )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("WaterVaporSaturationInClouds"),
//      DESCRIPTION(
//        "Calculates the water vapor saturation volume mixing ratio (VMR) in the\n"
//        "vertical range where liquid or ice clouds are in the atmosphere.\n"
//        "At the pressure/altitude grid points where the liquid water content (LWC)\n"
//        "or ice water content (IWC) of the clouds (tags 'liquidcloud' and 'icecloud')\n"
//           "is larger than zero the H2O-VMR is set to liquid water/ice saturation VMR.\n"
//           "The saturation pressure is calculated according to Goff-Gratch equations.\n"
//        ),
//      OUTPUT(   vmrs_ , p_abs_                         ),
//      INPUT(    vmrs_ , p_abs_ , t_abs_ , tgs_         ),
//      GOUTPUT(                                         ),         
//      GINPUT(                                          ),
//      KEYWORDS(                                        ),
//      TYPES(                                           )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("vmrsScale"),
//      DESCRIPTION(
//           "Scales the vmr input of the tgs given in scaltgs by the\n"
//        "factors given in scalfac.\n"
//        "\n"
//        "Keywords:\n"
//        "   scaltgs : subgroup of tags which has to be scaled.\n"
//        "   scalfac : the factor with which vmr to be scaled.\n"
//        "\n"
//        ),
//      OUTPUT( vmrs_ ),
//      INPUT(  tgs_, vmrs_  ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "scaltgs", "scalfac"),
//      TYPES( Array_String_t, Vector_t)));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("h2o_absSet"),
//      DESCRIPTION(
//           "Sets h2o_abs to the profile of the first tag group containing\n"
//        "water.\n" 
//        "\n"
//           "This is necessary, because for example *absCalc* requires h2o_abs\n"
//        "to contain the water vapour profile(the reason for this is the\n"
//           "calculation of oxygen line brodening requires water vapour profile).\n"
//        "Then this function can be used to copy the profile of the first tag\n"
//           "group of water.\n"
//        "\n"
//        ),
//      OUTPUT( h2o_abs_ ),
//      INPUT(  tgs_, vmrs_  ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("n2_absSet"),
//      DESCRIPTION(
//           "Sets n2_abs to the profile of the first tag group containing\n"
//        "molecular nitrogen. See *h2o_absSet* for more details.\n"
//        "\n"
//        ),
//      OUTPUT(     n2_abs_ ),
//      INPUT(  tgs_, vmrs_  ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));




// //=== Absorption methods ===============================================

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME( "absCalc" ),
//      DESCRIPTION(
//         "Calculate absorption coefficients. \n"
//         "\n"
//         "This function calculates both, the total absorption (*abs*)\n"
//         "and the absorption per tag group (*abs_per_tg*).\n"
//             ) ,
//      OUTPUT(abs_  , abs_per_tg_ ),
//      INPUT(tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_, 
//               lines_per_tg_, lineshape_,
//            cont_description_names_, cont_description_models_, 
//               cont_description_parameters_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("absCalcFromXsec"),
//      DESCRIPTION(
//                  "Calculate absorption coefficients from cross sections.\n"
//                  "\n"
//                  "This calculates both the total absorption and the\n"
//                  "absorption per tag group. \n"
//                  "This method calls three other  methods:\n"
//                  "1. *xsec_per_tgInit* - initialize *xsec_per_tg* \n"
//                  "2. *xsec_per_tgAddLine* - calculate cross sections per \n"
//                  "                   tag group for line spectra.\n"
//                  "3. *xsec_per_tgAddConts* - calculate cross sections per \n"
//                  "                   tag group for continua.\n"
//                  "Then it calculates the absorption coefficient by multiplying\n"
//                  "the cross section by VMR.\n"
//                     "This is done once for each tag group (output: *abs_per_tg*)\n"
//                  "and for the sum of all tag group to get the total absorption\n"
//                  "coefficient (output: *abs*)\n"
//                  ),
//      OUTPUT(     abs_  , abs_per_tg_ ),
//      INPUT(      xsec_per_tg_, vmrs_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME( "xsec_per_tgInit" ),
//      DESCRIPTION(
//         "Initialize *xsec_per_tg*.\n"
//         "\n"
//         "The initialization is\n"
//         "necessary, because methods *xsec_per_tgAddLines*\n"
//         "and *xsec_per_tgAddConts* just add to *xsec_per_tg*.\n"
//         "The size is determined from *tgs*.\n"
//         ),
//      OUTPUT( xsec_per_tg_ ),
//      INPUT(tgs_, f_mono_, p_abs_),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("xsec_per_tgAddLines"),
//      DESCRIPTION(
//                  "Calculate cross sections per tag group for line spectra.\n"
//                 ),
//      OUTPUT(     xsec_per_tg_                             ),
//      INPUT(      tgs_, f_mono_, p_abs_, t_abs_, h2o_abs_, vmrs_, 
//                  lines_per_tg_, lineshape_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("xsec_per_tgAddConts"),
//      DESCRIPTION(
//                  "Calculate cross sections per tag group for continua.\n"
//                      ),
//      OUTPUT(     xsec_per_tg_                             ),
//      INPUT(      tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_,
//                  cont_description_names_, cont_description_parameters_,
//                     cont_description_models_),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));


// //=== Methods operating on absorption ========================================

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("abs_per_tgReduce"),
//      DESCRIPTION(
//                  "Reduces absorption coefficients. Only absorption\n"
//                  "coefficients for which weighting functions are\n"
//                  "calculated are kept in memory.\n"
//                  ),
//      OUTPUT(     abs_per_tg_ ),
//      INPUT(      abs_per_tg_, tgs_, wfs_tgs_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));



}

