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


#if HAVE_CONFIG_H
#include <config.h>
#endif

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

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_vec_partCalc"),
	DESCRIPTION
	(
	 "This function sums up the absorption vectors for all particle \n"
	 "types weighted with particle number density.\n"
	 "\n"
	 "The output of this method is *abs_vec_part* (stokes_dim).\n"
	 "The inputs are the absorption vector for the single particle type \n"
	 "*abs_vec_spt* (part_types, stokes_dim) and the local particle\n"
	 " number densities for all particle types namely the \n"
	 "*pnd_field* (part_types, p_grid, lat_grid, lon_grid, ) for given \n"
	 "*p_grid*, *lat_grid*, and *lon_grid*. The particle types required \n"
	 "are specified in the control file.  \n"
	 ),
	OUTPUT(abs_vec_part_),
	INPUT(abs_vec_spt_, pnd_field_, atmosphere_dim_, scat_p_index_, 
	      scat_lat_index_, scat_lon_index_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_vec_sptCalc"),
  	DESCRIPTION
	(
	 "This method calculates the absorption vector for a single particle \n"
	 "type.\n"
	 "\n"
	 "All the elements of absorption vector for a particle can be expressed\n"
	 "in terms of extinction matrix and phase matrix. So it is necessary \n"
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
	OUTPUT(abs_vec_spt_),
	INPUT(abs_vec_spt_,ext_mat_spt_,pha_mat_spt_,scat_za_grid_,
	      scat_aa_grid_),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_vecCalc"),
  	DESCRIPTION
	(
	 "This function sums up the absorption vectors of particle and gas\n"
	 "and gives the total absorption vector.\n"
	 "\n"
	 "The method *abs_vec_partCalc* calculates the particle absorption\n"
	 "vector *abs_vec_part*.  In the case of gases it can either be \n"
	 "calculated using a similar method *abs_vec_gasCalc* or can be \n"
	 "be read in from a file.  The second option can really save some\n"
	 "computation time.\n"
	 "\n"
	 "The output of this method is *abs_vec* (stokes_dim). The inputs\n"
	 "the particle absorption vector *abs_vec_part*( stokes_dim )\n"
	 "and the gaseous absorption vector *abs_vec_gas* ( stokes_dim )\n"
	 ),
	OUTPUT( abs_vec_  ),
        INPUT( abs_vec_part_, abs_vec_gas_  ),
	GOUTPUT( ),
	GINPUT( ),
	KEYWORDS( ),
	TYPES( )));

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
	 " *a_pos*, *lat* as the second element and *lon* as third element.\n"
	 "However, the length of *a_pos* is set to *atmosphere_dim* and\n"
	 " keyword arguments for dimensions not used are ignored.\n"
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
      ( NAME("CloudboxOff"),
	DESCRIPTION
        (
         "Deactivates the cloud box. \n"
         "\n"
         "The function sets *cloudbox_on* to 0, and *cloudbox_limits* to be\n"
         "a an empty vector."
        ),
	OUTPUT( cloudbox_on_, cloudbox_limits_ ),
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
	INPUT( atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, 
               blackbody_ground_ ),
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
	INPUT(i_field_, i_field_old_, cloudbox_limits_, scat_za_grid_, 
              scat_aa_grid_, stokes_dim_, atmosphere_dim_),
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
      ( NAME("e_groundSet"),
	DESCRIPTION
	(
	 "Sets the ground emissivity to a constant value.\n"
	 "\n"
	 "The size of *e_ground* is automatically determined by the\n"
	 "latitaude , longitude and frequency grids.\n"
	 "\n"
	 "If the ground shall be treated as a blackbody, use the\n"
	 "function *SetBlackbodyGround*.\n"
	 "\n"
         "Keywords: \n"
         "   value : The ground emissivity."
	),
	OUTPUT( e_ground_ ),
	INPUT( atmosphere_dim_, lat_grid_, lon_grid_, f_grid_ ),
	GOUTPUT( ),
	GINPUT( ),
	KEYWORDS( "value" ),
	TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ext_matCalc"),
  	DESCRIPTION
	(
	 "This function sums up the extinction matrices of particle and gas\n"
	 "and gives the total extinction matrix.\n"
	 "\n"
	 "The extinction due to particle is the sum of scattering and \n"
	 "absorption whereas the extinction due to gas is due only to \n"
	 "absorption. The method *ext_mat_partCalc* calculates the particle\n" 
	 "extinction matrix *ext_mat_part*. In the case of gases it can\n"
	 "either be calculated in a similar method *ext_mat_gasCalc* or\n"
	 "can be read in from a file.  The second option can really save\n"
	 "some computation time.\n"
	 "\n"
	 "The output of this method is *ext_mat* (stokes_dim, stokes_dim).\n"
	 "The inputs are the particle extinction matrix *ext_mat_part*\n"
	 "(stokes_dim,  stokes_dim) and the gaseous extinction matrix \n"
	 "*ext_mat_gas* (stokes_dim,  stokes_dim).\n"
	 ),
	OUTPUT( ext_mat_  ),
        INPUT( ext_mat_part_, ext_mat_gas_  ),
	GOUTPUT( ),
	GINPUT( ),
	KEYWORDS( ),
	TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ext_mat_partCalc"),
  	DESCRIPTION
	(
	 "This function sums up the extinction matrices for all particle \n"
	 "types weighted with particle number density.\n"
	 "\n"
	 "The output of this method is *ext_mat_part* (stokes_dim, stokes_dim).\n"
	 "The inputs are the extinction matrix for the single particle type \n"
	 "*ext_mat_spt* (part_types, stokes_dim, stokes_dim) and the local \n"
	 "particle number densities for all particle types namely the \n"
	 "*pnd_field* (part_types, p_grid, lat_grid, lon_grid ) for given \n"
	 "*p_grid*, *lat_grid*, and *lon_grid*. The particle types required \n"
	 "are specified in the control file.  \n"
	 ),
	OUTPUT( ext_mat_part_  ),
        INPUT( ext_mat_spt_, pnd_field_, atmosphere_dim_, scat_p_index_, 
	       scat_lat_index_, scat_lon_index_),
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
	 "*ext_mat_spt*(Tensor 3, size = [Npt,stokes_dim,stokes_dim])and the \n"
	 "inputs are *ext_mat_spt*,*amp_mat*(Tensor 6, Size=[Npt,Nza,Naa,Nza,Naa,8]), \n"
	 "*scat_za_index*,*scat_aa_index*,*scat_f_index* and *scat_f_grid*. \n"
	 "\n"
	 "The variables *scat_za_index* and *scat_aa_index picks the right \n"
	 "element of the Tensor *amp_mat*. *scat_f_grid* and *scat_f_index* picks \n"
	 "the right frequeny for calculation. Frequeny is needed because the\n"
	 "computation of extinction matrix from amplitude matrix involves \n"
	 "multiplication by wavelength.  Then this method calls the \n"
	 "function amp2ext which does the actual physics, that of computing\n"
	 "the elements of extinction matrix from the elements of amplitude \n"
	 "matrix.\n"
	 ),
	OUTPUT( ext_mat_spt_  ),
        INPUT( ext_mat_spt_,amp_mat_, scat_za_index_, scat_aa_index_,
	       scat_f_index_, f_grid_  ),
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
      ( NAME("GroundSetToBlackbody"),
	DESCRIPTION
        (
         "Sets *blackbody_ground* to 1 and creates *e_ground* as a tensor\n"
	 "with size [1,1,1] having the value 1."
        ),
	OUTPUT( blackbody_ground_, e_ground_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

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
	OUTPUT(i_field_, ppath_step_, i_field_old_, scat_field_, sca_vec_, 
               stokes_vec_, planck_function_, l_step_,
               convergence_flag_, pha_mat_part_, pha_mat_spt_, abs_vec_spt_,
               ext_mat_spt_, ext_mat_, abs_vec_, scat_p_index_, 
               scat_lat_index_, scat_lon_index_),
	INPUT(ext_mat_agenda_, abs_vec_agenda_, convergence_test_agenda_,
              ppath_step_agenda_, scat_rte_agenda_, 
              amp_mat_, cloudbox_limits_,
              scat_za_grid_, scat_aa_grid_, 
	      p_grid_, lat_grid_, lon_grid_, t_field_, z_field_,
	      r_geoid_, f_grid_, scat_f_index_, 
	      stokes_dim_, atmosphere_dim_, pnd_field_, part_types_),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "i_fieldUpdate1D" ),
	DESCRIPTION
        (
	 "Updates the i_field during the iteration. It performs the RT \n"
         "calculation using a fixed value for the scattering integral stored \n"
         "in *scat_field*.\n"
         "   " 
        ),
	OUTPUT(i_field_, ppath_step_, stokes_vec_, 
               sca_vec_, planck_function_, l_step_,
               abs_vec_spt_, ext_mat_spt_, ext_mat_, abs_vec_,
               scat_p_index_),
	INPUT(ext_mat_agenda_, abs_vec_agenda_, ppath_step_agenda_,
              scat_rte_agenda_, amp_mat_, scat_field_, cloudbox_limits_,
              scat_za_grid_, p_grid_, t_field_, z_field_, 
              r_geoid_, f_grid_, scat_f_index_, 
	      pnd_field_, stokes_dim_, atmosphere_dim_, part_types_,
              pha_mat_spt_),
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
      ( NAME("i_spaceCBR"),
	DESCRIPTION
	(
	 "Sets *i_space* to hold cosmic background radiation."
	),
	OUTPUT( i_space_ ),
	INPUT( f_grid_ ),
	GOUTPUT( ),
	GINPUT( ),
	KEYWORDS( ),
	TYPES( )));

  
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
	OUTPUT(),
        INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES(),
	AGENDAMETHOD( true )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixFillWithVector"),
	DESCRIPTION
        (
         "Forms a matrix by repeating a vector.\n"
         "\n"
         "The vector can either form the rows or the columns of the matrix.\n"
	 "The direction of the vector inside the matrix is selected by\n"
         "setting the size determined by the vector length to 0. For \n"
         "example, if the keyword *ncols* is set to 0, the vector will be\n"
         "put in as rows of the matrix and the number of rows will equal\n"
         "*nrows*.\n"
	 "\n"
	 "One, but only one, keyword argument must be 0.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied. \n"
         "\n"
         "Keyword: \n"
         "   nrows : Number of rows in the matrix.\n"
         "   ncols : Number of columns in the matrix. "
        ),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Matrix_ ),
	GINPUT( Vector_ ),
	KEYWORDS( "nrows", "ncols"   ),
	TYPES(    Index_t, Index_t )));

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
      ( NAME( "pha_mat_partCalc" ),
	DESCRIPTION
        (
	 "This function sums up the phase matrices for all particle\n"
	 "types weighted with particle number density.\n"
	 "\n"
	 "The output of this method is *pha_mat_part* (Nza, Naa, stokes_dim,\n"
	 "stokes_dim). The inputs are the phase matrix for the single particle\n"
	 "type *pha_mat_spt* (part_types, Nza, Naa, stokes_dim, stokes_dim)\n"
	 "and the local particle  number densities for all particle types namely \n"
	 "the *pnd_field* (part_types, p_grid, lat_grid, lon_grid ) for given\n"
	 "*p_grid*, *lat_grid*, and *lon_grid*. The particle types required \n"
	 "are specified in the control file.\n"
	 ),
	OUTPUT(pha_mat_part_),
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
	INPUT(pha_mat_spt_, amp_mat_, scat_za_index_, scat_aa_index_ ),
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
	OUTPUT( ppath_ ),
	INPUT( atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_, 
               r_geoid_, z_ground_, blackbody_ground_, 
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
         "points are included for crossings with the grids, tangent points\n"
	 "and points of ground intersections.\n"
	 "\n"
	 "Functions to determine propagation path steps, as this function,\n"
	 "are normally not called directly. They are instead normally part\n"
	 "of *ppath_step_agenda* and then called from *ppathCalc*.\n"
	 "\n"
	 "As functions of this kind should very seldom be called directly,\n"
	 "and that the functions can be called a high number of times, these\n"
	 "functions do not perform any checks of the input that give\n" 
	 "detailed error messages, but asserts are performed (if turned on).\n"
	 "\n"
         "For further information, type see the on-line information for\n"
	 "*ppath_step_agenda* (type \"arts -d ppath_step_agenda\")."
        ),
	OUTPUT( ppath_step_ ),
	INPUT( ppath_step_, atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, 
	       z_field_, r_geoid_, z_ground_, blackbody_ground_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ppath_stepGeometricWithLmax" ),
	DESCRIPTION
        (
	 "As *ppath_stepGeometric* but with a length criterion for the\n"
         "distance between the path points.\n"
	 "\n"
	 "This function works as *ppath_stepGeometric* but additional points\n"
	 "are included in the propagation path to ensure that the distance\n"
	 "along the path between the points does not exceed the selected\n"
	 "maximum length. The length criterion is set by the keyword\n"
	 "argument.\n"
	 "\n"
         "See further the on-line information for *ppath_stepGeometric*\n"
	 "(type \"arts -d ppath_stepGeometric\").\n"
	 "\n"
         "Keywords: \n"
         "   lmax : Maximum allowed length between path points.\n"
        ),
	OUTPUT( ppath_step_ ),
	INPUT( ppath_step_, atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, 
	       z_field_, r_geoid_, z_ground_, blackbody_ground_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "lmax"),
	TYPES(    Numeric_t )));

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
      ( NAME( "RteCalc" ),
	DESCRIPTION
        (
	 "Main function for calculation of spectra and WFs.\n"
         "\n"
         "Text will be written (PE)."
        ),
	OUTPUT( y_ ),
	INPUT( atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_,
               t_field_, r_geoid_, z_ground_, t_ground_, e_ground_, 
               blackbody_ground_, cloudbox_on_,  cloudbox_limits_, f_grid_, 
               i_space_, antenna_dim_, mblock_za_grid_, mblock_aa_grid_,
               stokes_dim_, sensor_pos_, sensor_los_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

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
      ( NAME( "scat_integralCalc" ),
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
	 "these two workspace variables we calculate *pha_mat_part* with\n"
	 "the method *pha_mat_partCalc*.  \n"
	 ),
	OUTPUT( scat_field_, pha_mat_part_ ),
	INPUT( i_field_, pha_mat_spt_, pnd_field_, scat_za_grid_, 
	       scat_aa_grid_, p_grid_, lat_grid_, lon_grid_, stokes_dim_,
	       atmosphere_dim_, cloudbox_limits_ ),
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
      ( NAME("stokes_vecGeneral"),
	DESCRIPTION
	(
	 "Calculate vector radiative transfer with fixed scattering integral."
         "\n"
         "This function computes the radiative transfer for a thin layer.\n"
         "It is a general function which works for both, the vector and the \n"
         "scalar RTE. But for the scalar equation it is more efficient to \n"
         "the method  *stokes_vecScalar*.\n"
         "All coefficients and the scattered field vector are assumed to be\n"
         "constant inside the grid cell/layer.\n"
         "Then an analytic solution can be found (see AUG for details).\n"
	),
	OUTPUT(stokes_vec_),
	INPUT(ext_mat_, abs_vec_, sca_vec_, l_step_, planck_function_,
              stokes_dim_),
	GOUTPUT( ),
	GINPUT( ),
	KEYWORDS( ),
	TYPES( )));

md_data_raw.push_back     
    ( MdRecord
      ( NAME("stokes_vecScalar"),
	DESCRIPTION
	(
	 "Calculate scalar radiative transfer with fixed scattering integral."
         "\n"
         "This function computes the radiative transfer for a thin layer. All\n"
         "coefficients and the scattered field vector  are assumed to be \n"
         "constant inside the grid cell/layer.\n"
         "Then an analytic solution can be found (see AUG for details).\n"
	),
	OUTPUT(stokes_vec_),
	INPUT(ext_mat_, abs_vec_, sca_vec_, l_step_, planck_function_,
              stokes_dim_),
	GOUTPUT( ),
	GINPUT( ),
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
      ( NAME("VectorSetTakingLengthFromVector"),
	DESCRIPTION
        (
         "Creates a workspace vector with the same length as another vector,\n"
         "and sets all values of the new vector to the specified value. \n"
         "\n"
         "A common usage of the function should be: \n"
         "  VectorSetLengthFromVector(e_ground,f_mono){value=0.75} \n"
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
	INPUT(),
	GOUTPUT( ),
	GINPUT(  Any_ ),
	KEYWORDS( "filename" ),
	TYPES(    String_t   ),
	AGENDAMETHOD(   false ),
	SUPPRESSHEADER( true  )));











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
//   	DESCRIPTION
// 	(
// 	 "Sets lines_per_tg to empty line lists.\n"
// 	 "\n"
// 	 "You can use this method to set lines per tag if you do not reall want\n"
// 	 "to compute line spectra. Formally, absCalc will still require\n"
// 	 "lines_per_tg to be set.\n"
// 	 ),
// 	OUTPUT(   lines_per_tg_      ),
// 	INPUT(    tgs_        ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(  ),
// 	TYPES(    )));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgReadFromCatalogues"),
//   	DESCRIPTION(
// 		    "This method can read lines from different line \n"
// 		    "catalogues.\n"
// 		    "\n"
// 		    "For each tag group, you can specify which catalogue\n"
// 		    "to use. Because the method creates lines_per_tg directly,\n"
// 		    "it replaces for example thefollowing two method calls:\n"
// 		    "  - linesReadFromHitran\n"
// 		    "  - lines_per_tgCreateFromLines\n"
// 		    "   This method needs as input WSVs the list of tag \n"
// 		    "groups. Keyword parameters must specify the names of\n"
// 		    "the catalogue files to use and the matching formats.\n"
// 		    "Names can be anything, formats can currently be \n"
// 		    "HITRAN96, MYTRAN2, JPL, or ARTS. Furthermore, keyword\n"
// 		    "parameters have to specify minimum and maximum \n"
// 		    "frequency for each tag group. To safe typing, if there\n"
// 		    "are less elements in the keyword parameters than there\n"
// 		    "are tag groups, the last parameters are applied to all\n"
// 		    "following tag groups.\n"
// 		    "\n"
// 		    "Example usage:\n"
// 		    "\n"
// 		    "lines_per_tgReadFromCatalogues{\n"
// 		    "  filenames = [ \"../data/cat1.dat\", \"../data/cat2.dat\" ]\n"
// 		    "  formats   = [ \"MYTRAN2\",          \"HITRAN96\"         ]\n"
// 		    "  fmin      = [ 0,                  0                  ]\n"
// 		    "  fmax      = [ 2000e9,             100e9              ]\n"
// 		    "}\n"
// 		    "   In this example, lines for the first tag group will\n"
// 		    "be taken from cat1, lines for all other tag groups \n"
// 		    "will be taken from cat2.\n"
// 		    "   This methods allows you for example to use a \n"
// 		    "special line file just for water vapor lines. This\n"
// 		    "could be the  improved water vapor line file \n"
// 		    "generated by Thomas Kuhn.\n"
// 		    "   Catalogues are only read once, even if several tag\n"
// 		    "groups have the same catalogue. However, in that case\n"
// 		    "the frequency ranges MUST be the same. (If you want \n"
// 		    "to do fine-tuning of the frequency ranges, you can do \n"
// 		    "this inside the tag definitions, e.g., \"H2O-*-0-2000e9\".)\n"
// 		    "   This function uses the various reading routines\n"
// 		    "(linesReadFromHitran, etc.), as well as\n"
// 		    "lines_per_tgCreateFromLines.\n"
// 		    "\n"
// 		    "Keywords: \n"
// 		    "   filenames = Name (and path) of the catalogue files.\n"
// 		    "   formats   = allowed formats are HITRAN96,MYTRAN2,JPL,ARTS \n"
// 		    "   fmin      = Minimum frequency for lines to read in Hz.\n"
// 		    "   fmax      = Maximum frequency for lines to read in Hz.\n"),
// 	OUTPUT(   lines_per_tg_      ),
// 	INPUT(    tgs_        ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS( "filenames",    "formats",      "fmin",   "fmax" ),
// 	TYPES(    Array_String_t, Array_String_t, Vector_t, Vector_t)));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesReadFromHitran"),
//   	DESCRIPTION(
// 		    "Read all the lines from a HITRAN catalogue file in the \n"
// 		    "given frequency range. Otherwise a runtime error will be\n"
// 		    "thrown\n"
// 		    "\n"
// 		    "Please note that all lines must correspond\n"
// 		    "to the legal species / isotope combinations\n"
// 		    "\n"
// 		    "Keywords: \n"
// 		    "   filename = Name (and path) of the catalogue file.\n"
// 		    "   fmin     = Minimum frequency for lines to read in Hz.\n"
// 		    "   fmax     = Maximum frequency for lines to read in Hz."),
// 	OUTPUT(   lines_   ),
// 	INPUT(),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS( "filename",  "fmin",    "fmax"),
// 	TYPES(    String_t,    Numeric_t, Numeric_t)));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesReadFromMytran2"),
//   	DESCRIPTION(
// 		    "Read all the lines from a MYTRAN2 catalogue file in the \n"
// 		    "given frequency range. Otherwise a runtime error will be\n"
// 		    "thrown\n"
// 		    "\n"
// 		    "Please note that all lines must correspond\n"
// 		    "to the legal species / isotope combinations\n"
// 		    "\n"
// 		    "Keywords: \n"
// 		    "   filename = Name (and path) of the catalogue file.\n"
// 		    "   fmin     = Minimum frequency for lines to read in Hz.\n"
// 		    "   fmax     = Maximum frequency for lines to read in Hz."),
// 	OUTPUT(   lines_   ),
// 	INPUT(),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS( "filename",  "fmin",    "fmax"),
// 	TYPES(    String_t,    Numeric_t, Numeric_t)));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesReadFromJpl"),
//   	DESCRIPTION(
// 		    "Read all the lines from a JPL catalogue file in the \n"
// 		    "given frequency range. Otherwise a runtime error will be\n"
// 		    "thrown\n"
// 		    "\n"
// 		    "Please note that all lines must correspond\n"
// 		    "to the legal species / isotope combinations.\n"
// 		    "\n"
// 		    "Keywords: \n"
// 		    "   filename = Name (and path) of the catalogue file.\n"
// 		    "   fmin     = Minimum frequency for lines to read in Hz.\n"
// 		    "   fmax     = Maximum frequency for lines to read in Hz."),
// 	OUTPUT(   lines_   ),
// 	INPUT(),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS( "filename",  "fmin",    "fmax"),
// 	TYPES(    String_t,    Numeric_t, Numeric_t)));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesReadFromArts"),
//   	DESCRIPTION(
// 		    "Read all the lines from an Arts catalogue file in the \n"
// 		    "given frequency range. Otherwise a runtime error will be\n"
// 		    "thrown \n"
// 		    "\n"
// 		    "Please note that all lines must correspond\n"
// 		    "to the legal species / isotope combinations\n"
// 		    "\n"
// 		    "Keywords: \n"
// 		    "   filename = Name (and path) of the catalogue file.\n"
// 		    "   fmin     = Minimum frequency for lines to read in Hz.\n"
// 		    "   fmax     = Maximum frequency for lines to read in Hz."),
// 	OUTPUT(   lines_   ),
// 	INPUT(),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS( "filename",  "fmin",    "fmax"),
// 	TYPES(    String_t,    Numeric_t, Numeric_t)));
  
//   // FIXME: Remove this one.
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesElowToJoule"),
//   	DESCRIPTION(
// 		    "Just a little helper to convert the lower state energy from cm^-1\n"
// 		    "(ARTSCAT-2) to Joule (ARTSCAT-3). This should be removed soon\n"),
// 	OUTPUT(   lines_   ),
// 	INPUT(),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS( ),
// 	TYPES(    )));
      
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgCreateFromLines"),
//   	DESCRIPTION(
// 		    "Split lines up into the different tag groups.\n"
// 		    "\n"
// 		    "The tag groups are tested in the order in which they are\n" 
// 		    "specified in the controlfile. The lines are assigned to \n"
// 		    "the tag groups in the order as the groups  are specified.\n"
// 		    "That means if you do [\"O3-666\",\"O3\"],the last group O3 \n"
// 		    "gets assigned all the O3 lines that do not fit in the first group."),
// 	OUTPUT(   lines_per_tg_      ),
// 	INPUT(    lines_, tgs_ ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(),
// 	TYPES()));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgAddMirrorLines"),
//   	DESCRIPTION(
// 		    "Adds mirror lines at negative frequencies to the *lines_per_tg*.\n"
// 		    "\n"
// 		    "For each line at frequency +f in *lines_per_tg* a corresponding\n"
// 		    "entry at frequency -f is added to *lines_per_tg*.The mirror \n"
// 		    "lines are appended to the line lists after the original lines."),
// 	OUTPUT(   lines_per_tg_      ),
// 	INPUT(    lines_per_tg_      ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(),
// 	TYPES()));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgCompact"),
//   	DESCRIPTION(
// 		    "Removes all lines outside the defined lineshape cutoff frequency\n"
// 		    "from the *lines_per_tg*. This can save computation time.\n"
// 		    "It should be particularly useful to call this method after\n"
// 		    "*lines_per_tgAddMirrorLines*."),
// 	OUTPUT(   lines_per_tg_      ),
// 	INPUT(    lines_per_tg_, lineshape_, f_mono_  ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(),
// 	TYPES()));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesWriteAscii"),
//   	DESCRIPTION(
//                     "Writes the workspace variable *lines* to an ASCII file.\n"
//                     "\n"
// 		    "The content of the workspace variable 'lines`\n"
// 		    "The content of the workspace variable *lines*\n"
// 		    "is written in ARTS line format to the file with\n"
//                     "the specified name. If the filename is omitted, the\n"
//                     "lines are written to <basename>.lines.aa.\n"
//                     "\n"
//                     "Keywords: \n"
//                     "   filename : Name of the output file.\n"
//                     ), 
// 	OUTPUT(),
// 	INPUT( lines_ ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS( "filename" ),
// 	TYPES(    String_t   )));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgWriteAscii"),
//   	DESCRIPTION(
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
// 	OUTPUT(),
// 	INPUT( lines_per_tg_ ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS( "filename" ),
// 	TYPES(    String_t   )));
  
//   md_data_raw.push_back
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
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("tgsDefineAllInScenario"),
//   	DESCRIPTION
// 	(
// 	 "Define one tag group for each species known to ARTS and included in an\n"
// 	 "atmospheric scenario.\n"
// 	 "\n"
// 	 "You can use this as an alternative to tgsDefine if you want to make an\n"
// 	 "absorption calculation that is as complete as possible. The method\n"
// 	 "goes through all defined species and tries to open the VMR file. If\n"
// 	 "this works the tag is included, otherwise it is skipped.\n"
// 	 "\n"
// 	 "Keywords:\n"
// 	 "   basename : The name and path of a particular atmospheric scenario.\n"
// 	 "              For example: /pool/lookup2/arts-data/atmosphere/fascod/tropical"
//            ),
// 	OUTPUT( tgs_ ),
// 	INPUT(),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS( "basename" ),
// 	TYPES(    String_t   )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lineshapeDefine"),
//   	DESCRIPTION(
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
//    	  "shape=[\"Lorentz\"]\n"
//           "normalizationfactor=[\"linear\"]\n"
//           "cutoff= [650e9]"
//           "\n"
//           "Keywords:\n"
//           "   shape               : The general profile according to an approximation.\n"
//           "   normalizationfactor : The multiplicative forefactor for the general profile.\n"
//           "   cutoff              : The frequency at which a cutoff can be made.\n"),
// 	OUTPUT( lineshape_ ),
// 	INPUT( tgs_ ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(  "shape",    "normalizationfactor",  "cutoff" ),
// 	TYPES(     String_t,        String_t,         Numeric_t )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lineshape_per_tgDefine"),
//   	DESCRIPTION(
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
// 	  "shape = [\"Lorentz\",\"Voigt_Kuntz6\"] \n"
// 	  "normalizationfactor= [\"linear\", \"quadratic\"] \n"
// 	  "cutoff = [ 650e9, -1 ]"
//           "\n"
//           "Keywords:\n"
//           "   shape               : The general profile according to an approximation.\n"
//           "   normalizationfactor : The multiplicative forefactor for the general profile.\n"
//           "   cutoff              : The frequency at which a cutoff can be made.\n"),
// 	OUTPUT( lineshape_ ),
// 	INPUT( tgs_ ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(  "shape",           "normalizationfactor",    "cutoff" ),
// 	TYPES(   Array_String_t,         Array_String_t,        Vector_t )));


// //=== Continuum methods ============================================

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("cont_descriptionInit"),
//   	DESCRIPTION
// 	(
// 	 "Initializes the two workspace variables for the continuum description,\n"
// 	 "*cont_description_names* and *cont_description_parameters*.\n"
// 	 " \n"
// 	 "This method does not really do anything, except setting the two\n"
// 	 "variables to empty Arrays. It is just necessary because the method\n"
// 	 "*cont_descriptionAppend* wants to append to the variables.\n"
// 	 "   Formally, the continuum description workspace variables are required\n"
// 	 "by the absorption calculation methods (e.g., *absCalc*). Therefore you\n"
// 	 "always have to call at least *cont_descriptionInit*, even if you do\n"
// 	 "not want to use any continua."
// 	 ),
// 	OUTPUT( cont_description_names_, 
//                 cont_description_models_,
//                 cont_description_parameters_ ),
// 	INPUT(),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(),
// 	TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("cont_descriptionAppend"),
//   	DESCRIPTION
// 	(
// 	 "Appends the description of a continuum model or a complete absorption\n"
// 	 "model to *cont_description_names* and *cont_description_parameters*.\n"
// 	 "\n"
// 	 "See online documentation for *cont_description_names* for a list of\n"
// 	 "allowed models and for information what parameters they require. See\n"
// 	 "file cont.arts in the doc/examples directory for usage examples and\n"
// 	 "default parameters for the various models. \n"
// 	 "\n"
// 	 "Keywords:\n"
// 	 "   name       : The name of a continuum model. Must match one of the models\n"
// 	 "                implemented in ARTS. \n"
//          "   option     : give here the option of this continuum/full model.\n"
// 	 "   parameters : A Vector containing the required number of parameters\n"
// 	 "                for the model given. The meaning of the parameters and\n"
// 	 "                how many parameters are required depends on the model.\n"
// 	 ),
// 	OUTPUT( cont_description_names_, 
//                 cont_description_models_,
//                 cont_description_parameters_ ),
// 	INPUT(  cont_description_names_, 
//                 cont_description_models_,
//                 cont_description_parameters_),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS( "tagname",  "model",   "userparameters" ),
// 	TYPES(    String_t,   String_t,   Vector_t         )));


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
// 	  "\n"
// 	  "Keywords:\n"
// 	  "   seltags   : Must be a sub group of tags which should be read from files.\n"
// 	  "   filenames : Names of the files containing VMR profiles of seltags.\n"
// 	  "   basename  : The name of a particular atmospheric scenario.\n"
// 	  "               See *raw_vmrsReadFromScenario* for details. Remaining\n"
// 	  "               VMRs will be read from the scenario.\n"
// 	  "\n"
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
//   	DESCRIPTION(
// 	  "Reads the individual VMR profile for each tag group from a standard\n"
// 	  "atmospheric scenario.\n" 
// 	  "\n"
//           "Five different atmospheric scenarios are available in arts data:\n"
//           "tropical, midlatitude-summer, midlatitude-winter, subartic-summer\n"
//           "and subartic-winter.\n"
// 	  "\n"
// 	  "   Files in the scenarios look like this: tropical.H2O.aa\n"
// 	  "\n"
// 	  "   The basename must include the path, i.e., the files can be anywhere,\n"
// 	  "but they must be all in the same directory.\n"
// 	  "   The profile is chosen by the species name. If you have more than one\n"
// 	  "tag group for the same species, the same profile will be used.\n"
// 	  "\n"
// 	  "Keywords:\n"
// 	  "   basename :The name and path of a particular atmospheric scenario.\n"
// 	  "   For example:\n"
// 	  "   /pool/lookup2/arts-data/atmosphere/fascod/tropical\n"
// 	  "\n"
// 	  ),
// 	OUTPUT(   raw_vmrs_    ),
// 	INPUT(    tgs_                 ),
// 	GOUTPUT(                       ),
// 	GINPUT(                        ),
// 	KEYWORDS( "basename"           ),
// 	TYPES(    String_t             )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("AtmFromRaw"),
//   	DESCRIPTION(
// 	  "Interpolates temperature, altitude, and VMRs to the pressure grid\n"
// 	  "given by p_abs.\n" 
// 	  "\n"
//           "The altitude is not used by the absorption routines,\n"
// 	  "but later on by the RT routines.\n"
// 	  "\n"
// 	  "Interpolations used: \n"
// 	  "\n"
// 	  "Temperature      : Linear interpolation in ln(p)\n"
// 	  "Altitude         : Linear interpolation in ln(p)\n"
// 	  "VMRs             : Linear interpolation in ln(p)\n"
// 	  "Cloud Parameters : Linear interpolation in ln(p)\n"
// 	  "\n"
// 	  ),
// 	OUTPUT(   t_abs_    , z_abs_   , vmrs_           ),
// 	INPUT(    tgs_, p_abs_    , raw_ptz_ , raw_vmrs_ ),
// 	GOUTPUT(                       			 ),         
// 	GINPUT(                        			 ),
// 	KEYWORDS(                             		 ),
// 	TYPES(                          		 )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("WaterVaporSaturationInClouds"),
//   	DESCRIPTION(
// 	  "Calculates the water vapor saturation volume mixing ratio (VMR) in the\n"
// 	  "vertical range where liquid or ice clouds are in the atmosphere.\n"
// 	  "At the pressure/altitude grid points where the liquid water content (LWC)\n"
// 	  "or ice water content (IWC) of the clouds (tags 'liquidcloud' and 'icecloud')\n"
//           "is larger than zero the H2O-VMR is set to liquid water/ice saturation VMR.\n"
//           "The saturation pressure is calculated according to Goff-Gratch equations.\n"
// 	  ),
// 	OUTPUT(   vmrs_ , p_abs_                         ),
// 	INPUT(    vmrs_ , p_abs_ , t_abs_ , tgs_         ),
// 	GOUTPUT(                       			 ),         
// 	GINPUT(                        			 ),
// 	KEYWORDS(                             		 ),
// 	TYPES(                          		 )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("vmrsScale"),
// 	DESCRIPTION(
//           "Scales the vmr input of the tgs given in scaltgs by the\n"
// 	  "factors given in scalfac.\n"
// 	  "\n"
// 	  "Keywords:\n"
// 	  "   scaltgs : subgroup of tags which has to be scaled.\n"
// 	  "   scalfac : the factor with which vmr to be scaled.\n"
// 	  "\n"
// 	  ),
// 	OUTPUT(	vmrs_ ),
// 	INPUT( 	tgs_, vmrs_  ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS( "scaltgs", "scalfac"),
// 	TYPES( Array_String_t, Vector_t)));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("h2o_absSet"),
// 	DESCRIPTION(
//           "Sets h2o_abs to the profile of the first tag group containing\n"
// 	  "water.\n" 
// 	  "\n"
//           "This is necessary, because for example *absCalc* requires h2o_abs\n"
// 	  "to contain the water vapour profile(the reason for this is the\n"
//           "calculation of oxygen line brodening requires water vapour profile).\n"
// 	  "Then this function can be used to copy the profile of the first tag\n"
//           "group of water.\n"
// 	  "\n"
// 	  ),
// 	OUTPUT(	h2o_abs_ ),
// 	INPUT( 	tgs_, vmrs_  ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(),
// 	TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("n2_absSet"),
// 	DESCRIPTION(
//           "Sets n2_abs to the profile of the first tag group containing\n"
// 	  "molecular nitrogen. See *h2o_absSet* for more details.\n"
// 	  "\n"
// 	  ),
// 	OUTPUT(	    n2_abs_ ),
// 	INPUT( 	tgs_, vmrs_  ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(),
// 	TYPES()));




// //=== Absorption methods ===============================================

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME( "absCalc" ),
// 	DESCRIPTION(
// 	   "Calculate absorption coefficients. \n"
// 	   "\n"
// 	   "This function calculates both, the total absorption (*abs*)\n"
// 	   "and the absorption per tag group (*abs_per_tg*).\n"
//             ) ,
// 	OUTPUT(abs_  , abs_per_tg_ ),
// 	INPUT(tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_, 
//               lines_per_tg_, lineshape_,
// 	      cont_description_names_, cont_description_models_, 
//               cont_description_parameters_ ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(),
// 	TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("absCalcFromXsec"),
// 	DESCRIPTION(
// 		    "Calculate absorption coefficients from cross sections.\n"
// 		    "\n"
// 		    "This calculates both the total absorption and the\n"
// 		    "absorption per tag group. \n"
// 		    "This method calls three other  methods:\n"
// 		    "1. *xsec_per_tgInit* - initialize *xsec_per_tg* \n"
// 		    "2. *xsec_per_tgAddLine* - calculate cross sections per \n"
// 		    "                   tag group for line spectra.\n"
// 		    "3. *xsec_per_tgAddConts* - calculate cross sections per \n"
// 		    "                   tag group for continua.\n"
// 		    "Then it calculates the absorption coefficient by multiplying\n"
// 		    "the cross section by VMR.\n"
//                     "This is done once for each tag group (output: *abs_per_tg*)\n"
// 		    "and for the sum of all tag group to get the total absorption\n"
// 		    "coefficient (output: *abs*)\n"
// 		    ),
// 	OUTPUT(	    abs_  , abs_per_tg_ ),
// 	INPUT( 	    xsec_per_tg_, vmrs_ ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(),
// 	TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME( "xsec_per_tgInit" ),
// 	DESCRIPTION(
// 	   "Initialize *xsec_per_tg*.\n"
// 	   "\n"
// 	   "The initialization is\n"
// 	   "necessary, because methods *xsec_per_tgAddLines*\n"
// 	   "and *xsec_per_tgAddConts* just add to *xsec_per_tg*.\n"
// 	   "The size is determined from *tgs*.\n"
// 	   ),
// 	OUTPUT( xsec_per_tg_ ),
// 	INPUT(tgs_, f_mono_, p_abs_),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(),
// 	TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("xsec_per_tgAddLines"),
// 	DESCRIPTION(
// 		    "Calculate cross sections per tag group for line spectra.\n"
// 		   ),
// 	OUTPUT(	    xsec_per_tg_                             ),
// 	INPUT( 	    tgs_, f_mono_, p_abs_, t_abs_, h2o_abs_, vmrs_, 
// 		    lines_per_tg_, lineshape_ ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(),
// 	TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("xsec_per_tgAddConts"),
// 	DESCRIPTION(
// 		    "Calculate cross sections per tag group for continua.\n"
//                      ),
// 	OUTPUT(	    xsec_per_tg_                             ),
// 	INPUT( 	    tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_,
// 		    cont_description_names_, cont_description_parameters_,
//                     cont_description_models_),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(),
// 	TYPES()));


// //=== Methods operating on absorption ========================================

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("abs_per_tgReduce"),
// 	DESCRIPTION(
// 		    "Reduces absorption coefficients. Only absorption\n"
// 		    "coefficients for which weighting functions are\n"
// 		    "calculated are kept in memory.\n"
// 		    ),
// 	OUTPUT(	    abs_per_tg_ ),
// 	INPUT( 	    abs_per_tg_, tgs_, wfs_tgs_ ),
// 	GOUTPUT(),
// 	GINPUT(),
// 	KEYWORDS(),
// 	TYPES()));



}

