/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>
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
#include "wsv.h"
#include "methods.h"


void define_md_data()
{
  // The variable md_data is defined in file methods_aux.cc.
  extern ARRAY<MdRecord> md_data;

  // Initialize to zero, just in case:
  md_data.clear();

  /* Here's an empty template record entry:

  md_data.push_back
    ( MdRecord
      ( NAME(""),
	DESCRIPTION(""),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(""),
	TYPES()
	));
  */


//======================================================================
//=== IO methods
//======================================================================

//=== INDEX ============================================================

  md_data.push_back
    ( MdRecord
      ( NAME("IntSet"),
	DESCRIPTION("Sets an integer workspace variable to the given value."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( int_t ),
	GINPUT(),
	KEYWORDS( "value" ),
	TYPES( int_t )));



//=== NUMERIC ==========================================================

  md_data.push_back
    ( MdRecord
      ( NAME("NumericSet"),
	DESCRIPTION("Sets a workspace variable of type Numeric to a value."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Numeric_ ),
	GINPUT(),
	KEYWORDS( "value" ),
	TYPES( Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("NumericWriteBinary"),
	DESCRIPTION("Writes a numeric value to a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "See ??? for details about the file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( Numeric_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("NumericReadBinary"),
	DESCRIPTION("Reads a numeric from a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "See ??? for details about the file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Numeric_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));



//=== VECTOR ==========================================================

  md_data.push_back
    ( MdRecord
      ( NAME("VectorSet"),
	DESCRIPTION("Creates a workspace vector with the specified length\n"
                    "and initializes the vector with the given value."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(VECTOR_),
	GINPUT(),
	KEYWORDS("length", "value"),
	TYPES(int_t, Numeric_t)));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorLinSpace"),
	DESCRIPTION("Creates a linearly spaced vector with defined spacing.\n"
                    "Format: VectorLinSpace(x){start,stop,step}\n"
		    "The first element of x is always start.\n"
		    "The next value is start+step etc.\n"
		    "Note that the last value can deviate from stop.\n"
		    "The step can be both positive and negative."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(VECTOR_),
	GINPUT(),
	KEYWORDS("start", "stop", "step"),
	TYPES(Numeric_t, Numeric_t, Numeric_t)));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorNLinSpace"),
	DESCRIPTION("Creates a vector with defined length, equally spaced\n"
                    "between the given values.\n"
		    "The length must be larger than 1."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(VECTOR_),
	GINPUT(),
	KEYWORDS("start", "stop", "n"),
	TYPES(Numeric_t, Numeric_t, int_t)));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorNLogSpace"),
	DESCRIPTION("Creates a vector with defined length, logarithmically\n"
                    "spaced between the given values.\n"
		    "The length must be larger than 1."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(VECTOR_),
	GINPUT(),
	KEYWORDS("start", "stop", "n"),
	TYPES(Numeric_t, Numeric_t, int_t)));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorCopy"),
	DESCRIPTION("Copies a vector."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( VECTOR_ ),
	GINPUT( VECTOR_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorWriteAscii"),
	DESCRIPTION("Writes a vector to an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( VECTOR_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorReadAscii"),
	DESCRIPTION("Reads a vector from an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( VECTOR_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorWriteBinary"),
	DESCRIPTION("Writes a vector to a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "See ??? for details about the file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( VECTOR_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorReadBinary"),
	DESCRIPTION("Reads a vector from a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "See ??? for details about the file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( VECTOR_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));



//=== MATRIX ==========================================================

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixCopy"),
	DESCRIPTION("Copies a matrix."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( MATRIX_ ),
	GINPUT( MATRIX_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixWriteAscii"),
	DESCRIPTION("Writes a matrix to an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( MATRIX_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixReadAscii"),
	DESCRIPTION("Reads a matrix from an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( MATRIX_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixWriteBinary"),
	DESCRIPTION("Writes a matrix to a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "See ??? for details about the file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( MATRIX_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixReadBinary"),
	DESCRIPTION("Reads a matrix from a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "See ??? for details about the file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( MATRIX_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));



//=== ARRAYofINDEX =====================================================

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexWriteBinary"),
	DESCRIPTION("Writes an index array to a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "See ??? for details about the file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( ARRAYofsizet_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexReadBinary"),
	DESCRIPTION("Reads an index array from a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "See ??? for details about the file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( ARRAYofsizet_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));



//=== ARRAYofVECTOR ====================================================

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfVectorWriteAscii"),
	DESCRIPTION("Writes an ARRAYofVECTOR to an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( ARRAYofVECTOR_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfVectorReadAscii"),
	DESCRIPTION("Reads an ARRAYofVECTOR from an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( ARRAYofVECTOR_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfVectorWriteBinary"),
	DESCRIPTION("Writes a vector array to a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "See ??? for details about the file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( ARRAYofVECTOR_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfVectorReadBinary"),
	DESCRIPTION("Reads a vector array from a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "See ??? for details about the file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( ARRAYofVECTOR_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));
  



//=== ARRAYofMATRIX ====================================================

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixWriteAscii"),
	DESCRIPTION("Writes an ARRAYofMATRIX to an ASCII file.\n"
		    "The filename can also be an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "The format is as follows:\n\n"
		    "# <comments>\n"
		    "<n_array_elements>\n"
		    "<n_rows> <n_columns>\n"
		    "<elements>\n"
		    "<n_rows> <n_columns>\n"
		    "<elements>\n"
		    "...\n\n"
		    "Example:\n"
		    "# Generated by arts-0.0.16, Apr 29 2000, 17:38:44\n"
		    "2\n"
		    "3 4\n"
		    "xx xx xx xx\n"
		    "xx xx xx xx\n"
		    "xx xx xx xx\n"
		    "2 2\n"
		    "yy yy\n"
		    "yy yy"),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( ARRAYofMATRIX_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixReadAscii"),
	DESCRIPTION("Reads an ARRAYofMATRIX from an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( ARRAYofMATRIX_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixWriteBinary"),
	DESCRIPTION("Writes a matrix array to a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "See ??? for details about the file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( ARRAYofMATRIX_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixReadBinary"),
	DESCRIPTION("Reads a matrix array from a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "See ??? for details about the file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( ARRAYofMATRIX_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));



//=== MAYBESPARSE ====================================================

  md_data.push_back
    ( MdRecord
      ( NAME("HmatrixReadAscii"),
	DESCRIPTION("Reads a H matrix from an ASCII file.\n"
		    "The filename can be specified or be an empty string\n"
		    "If empty, it is set to <basename>.<variable_name>.am."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Hmatrix_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));



//======================================================================
//=== Absorption methods
//======================================================================

//=== Spectroscopic methods ============================================

  md_data.push_back
    ( MdRecord
      ( NAME("linesReadFromHitran"),
  	DESCRIPTION(
          "Read all the lines from a HITRAN catalogue file that\n"
	  "correspond to legal species / isotope combinations.\n"
	  "\n"
	  "The output line array is not overwritten, but the new data\n"
	  "is appended!\n" 
	  "\n"
	  "This is mainly for testing, the method probably will become\n"
	  "more complicated later on.\n"
	  "\n"
	  "filename = Name (and path) of the catalogue file.\n"
	  "fmin     = Minimum frequency for lines to read in Hz.\n"
	  "fmax     = Maximum frequency for lines to read in Hz."),
	OUTPUT(   lines_   ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "filename",  "fmin",    "fmax"),
	TYPES(    string_t,    Numeric_t, Numeric_t)));

  md_data.push_back
    ( MdRecord
      ( NAME("lines_per_tgCreateFromLines"),
  	DESCRIPTION(
          "Split lines up into the different tag groups.\n"
	  "The tag groups are tested in the order in which they are\n" 
          "specified in the controlfile. The line is assigned to the\n" 
	  "first tag group that fits."),
	OUTPUT(   lines_per_tg_      ),
	INPUT(    lines_, tag_groups_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("linesWriteToFile"),
  	DESCRIPTION(
          "Write the content of the workspace variable lines to the\n"
	  "given file in ARTS line format."),
	OUTPUT(),
	INPUT( lines_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("lines_per_tgWriteToFile"),
  	DESCRIPTION(
          "Write the content of the workspace variable lines_per_tg to the\n"
	  "given file in ARTS line format.\n"
	  "The array dimension is handled in a similar way as by the\n"
	  "array of vector and matrix output function:\n"
	  "First an integer stating the number of tag groups.\n"
	  "Then an integer specifying the number of lines for the\n"
	  "first group. Then the other groups in similar fashion."),
	OUTPUT(),
	INPUT( lines_per_tg_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("tag_groupsDefine"),
  	DESCRIPTION(
          "Set up the list of tag groups.\n"
	  "Specify one string for each tag group that you want to create.\n"
	  "Inside the string, separate the tags by comma (plus optional blanks).\n"
	  "Example:\n"
	  "tag = [\"O3-666-500e9-501e9, O3-686\",\"H2O\",\"O2-*-*-*\"]"),
	OUTPUT( tag_groups_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "tags" ),
	TYPES(    ARRAY_string_t   )));



//=== Input Atmosphere methods ===========================================

  md_data.push_back
    ( MdRecord
      ( NAME("raw_vmrs_1dReadFromScenario"),
  	DESCRIPTION(
	  "Read the individual VMR profile for each tag group from a standard\n"
	  "atmospheric scenario. Files must look like this example:\n"
	  "<basename>.ClO.am\n"
	  "\n"
	  "The basename can include a path, i.e., the files can be anywhere,\n"
	  "but they must be all in the same directory.\n"
	  "\n"
	  "The profile is chosen by the species name. If you have more than one\n"
	  "tag group for the same species, the same profile will be used."
	  ),
	OUTPUT(   raw_vmrs_1d_    ),
	INPUT(    tag_groups_          ),
	GOUTPUT(                       ),
	GINPUT(                        ),
	KEYWORDS( "basename"           ),
	TYPES(    string_t             )));

  md_data.push_back
    ( MdRecord
      ( NAME("Atm2dFromRaw1D"),
  	DESCRIPTION(
	  "This method is not currently useful for anything, since\n"
	  "there is no method to calculate absorption from the 2D\n"
	  "parameters.\n"
	  "\n"
	  "Interpolates temperature, altitude, and VMRs to the pressure grid\n"
	  "given by p_abs. The altitude is not used by the absorption routines,\n"
	  "But later on by the RT routines."
	  "\n"
	  "Interpolations used: FIXME: Add these.f\n"
	  "Temperature [K]: \n"
	  "Altitude    [m]: \n"
	  "VMRs        [1]: \n"
	  "\n"
	  "Uses interp_lin(...)."
	  ),
	OUTPUT(   t_abs_2d_ , z_abs_2d_   , vmrs_2d_     ),
	INPUT(    p_abs_    , raw_ptz_1d_ , raw_vmrs_1d_ ),
	GOUTPUT(                       			 ),         
	GINPUT(                        			 ),
	KEYWORDS(                      			 ),
	TYPES(                         			 )));

  md_data.push_back
    ( MdRecord
      ( NAME("AtmFromRaw1D"),
  	DESCRIPTION(
	  "Interpolates temperature, altitude, and VMRs to the pressure grid\n"
	  "given by p_abs. The altitude is not used by the absorption routines,\n"
	  "But later on by the RT routines."
	  "\n"
	  "Interpolations used: FIXME: Add these.f\n"
	  "Temperature [K]: \n"
	  "Altitude    [m]: \n"
	  "VMRs        [1]: \n"
	  "\n"
	  "Uses interp_lin(...)."
	  ),
	OUTPUT(   t_abs_    , z_abs_   , vmrs_           ),
	INPUT(    p_abs_    , raw_ptz_1d_ , raw_vmrs_1d_ ),
	GOUTPUT(                       			 ),         
	GINPUT(                        			 ),
	KEYWORDS(                      			 ),
	TYPES(                         			 )));



//=== 1D absorption methods ===============================================

  md_data.push_back
    ( MdRecord
      ( NAME("absCalc"),
	DESCRIPTION("Calculate absorption coefficients. This\n"
		    "calculates both the total absorption and the\n"
		    "absorption per tag group"
		    "\n"
		    "Line shape function is hardwired, and temperature\n"
		    "dependence is not calculated. This method is\n"
		    "really only for demonstration."),
	OUTPUT(	    abs_  , abs_per_tg_                         ),
	INPUT( 	    f_mono_, p_abs_, t_abs_, vmrs_, lines_per_tg_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("refr_indexBoudourisDryAir"),
	DESCRIPTION("Calculates the refractive index for dry air at micro-\n"
		    "wave frequncies following Boudouris 1963.\n"
		    "The expression is also found in Chapter 5 of the\n"
		    "Janssen book."),
	OUTPUT(	    refr_index_ ),
	INPUT( 	    p_abs_, t_abs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));
  


//======================================================================
//=== LOS/RTE methods
//======================================================================

  md_data.push_back
    ( MdRecord
      ( NAME("losCalc"),
  	DESCRIPTION(
          "A general function to determine LOS for a 1D atmosphere.\n"
          "Refraction is selected by a flag and the refraction variables\n"
          "must be set when using this function. The ground altitude must\n"
          "also be specified."),
	OUTPUT( los_ ),
	INPUT( z_plat_ ,za_pencil_, l_step_, p_abs_, z_abs_, 
               refr_, l_step_refr_, refr_index_, z_ground_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("losNoRefraction"),
  	DESCRIPTION(
          "Determines the LOS for a 1D atmosphere without refraction.\n"
          "The ground altitude must be specified."),
	OUTPUT( los_ ),
	INPUT( z_plat_ ,za_pencil_, l_step_, p_abs_, z_abs_, z_ground_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("losUpward"),
  	DESCRIPTION(
          "Determines the LOS for a 1D atmosphere when neglecting refraction\n"
          "and there is no ground intersections. The typical case is upward\n"
          "observations, but the function could also be of interest for limb\n"
          "sounding observations strictly above about 20 km.\n"
          "Refraction and ground altitude variables are NOT needed."),
	OUTPUT( los_ ),
	INPUT( z_plat_ ,za_pencil_, l_step_, p_abs_, z_abs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("sourceCalc"),
  	DESCRIPTION(
          "Calculates source function values valid between the points of"
          "of a 1D LOS.\n" 
          "No scattering and local thermodynamic equilibrium are assumed,\n"
          "that is, the source function equals the Planck function.\n"
          "The source function is set to the mean of the Planck function at\n"
          "the two LOS points limiting the steps. The temperature at the LOS\n"
          "points is obtained by linear interpolation"),
	OUTPUT( source_ ),
	INPUT( los_, p_abs_, t_abs_, f_mono_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("transCalc"),
  	DESCRIPTION(
          "Calculates the transmission between the points of a 1D LOS.\n"
          "The absorption is assumed to vary linear between the LOS points."
          "The absorption at the LOS points is obtained by linear\n"
          "interpolation of abs."),
	OUTPUT( trans_ ),
	INPUT( los_, p_abs_, abs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("y_spaceStd"),
  	DESCRIPTION(
          "Standard selection for the radiation entering the atmosphere at\n"
          "the start of the LOS. The selections are:\n"
          "  0 no radiation\n"
          "  1 cosmic background radiation (planck for COSMIC_BG_TEMP)\n"
          "  2 solar radiation (planck for SUN_TEMP)\n"
          "COSMIC_BG_TEMP and SUN_TEMP are global variables, defined in\n"
          "constants.cc."),
	OUTPUT( y_space_ ),
	INPUT( f_mono_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "nr" ),
	TYPES( int_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("y_spacePlanck"),
  	DESCRIPTION(
          "Sets the radiation entering the atmosphere at the start of the\n"
          "LOS to the Planck function for the given temperature."),
	OUTPUT( y_space_ ),
	INPUT( f_mono_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS("temp"),
	TYPES(Numeric_t)));

  md_data.push_back
    ( MdRecord
      ( NAME("yRte"),
  	DESCRIPTION(
          "Solves the general radiative transfer equation (RTE) along the\n"
          "LOS. With other words, both absorption and emission are\n"
          "considered.\n"
          "This function requires that e_ground and t_ground are set."),
	OUTPUT( y_ ),
	INPUT( los_, f_mono_, y_space_, source_, trans_, e_ground_, 
                                                                  t_ground_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("yRteNoGround"),
  	DESCRIPTION(
          "This function can be used instead of yRte if there is no\n"
          "intersection with the ground.\n"
          "With other words, the ground emission and temperature are NOT\n"
          "needed when using this function."),
	OUTPUT( y_ ),
	INPUT( los_, f_mono_, y_space_, source_, trans_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("yBl"),
  	DESCRIPTION(
          "Calculates the total transmission throught the atmosphere,\n"
          "using the Beer-Lambert (BL) law.\n"
          "This function requires that e_ground is set."),
	OUTPUT( y_ ),
	INPUT( los_, trans_, e_ground_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("yBlNoGround"),
  	DESCRIPTION(
          "This function can be used instead of yBl if there is no\n"
          "intersection with the ground. The ground emission is NOT needed \n"
          "when using this function."),
	OUTPUT( y_ ),
	INPUT( los_, trans_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));



//======================================================================
//=== Weighting function (WF) methods
//======================================================================

  md_data.push_back
    ( MdRecord
      ( NAME("absloswfsCalc"),
  	DESCRIPTION(
          "Calculates absorption line of sight weighting functions (LOS WFs)\n"
          "for 1D.\n"
          "These WFs are the derivative of the monochromatic pencil beam\n"
          "intensity with respect to the absorption at the LOS points.\n"
          "See further the ARTS user guide.\n"
          "This function requires that e_ground and t_ground are set."),
	OUTPUT( absloswfs_ ),
	INPUT( los_, source_, trans_, y_, y_space_, f_mono_, e_ground_, 
                                                                   t_ground_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("absloswfsNoGround"),
  	DESCRIPTION(
          "As klosCalc but does not need any ground variables"),
	OUTPUT( absloswfs_ ),
	INPUT( los_, source_, trans_, y_, y_space_, f_mono_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kSpecies"),
  	DESCRIPTION(
          "Calculates species 1D weighting functions for a single tag.\n"
          "The tag is selected by the parameter nr, which is the position\n"
          "for the tag in abs_per_tg.\n"
          "The avaliable units are\n"
          "  1 fractions of linearisation state \n"
          "  2 volume mixing ratio \n"
          "  3 number density"),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( los_, absloswfs_, p_abs_, t_abs_, tag_groups_, abs_per_tg_, 
                                                              vmrs_, k_grid_),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "nr",   "unit"  ),
	TYPES(    int_t,  int_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kSpeciesAll"),
  	DESCRIPTION(
          "Calculates species 1D weighting functions for all tags that\n"
          "are included in abs_per_tg. Units as for kSpecies."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( los_, absloswfs_, p_abs_, t_abs_, tag_groups_, abs_per_tg_, 
                                                              vmrs_, k_grid_),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "unit"  ),
	TYPES(    int_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kContAbs"),
  	DESCRIPTION(
          "Calculates 1D weighting functions for fit of continuum absorption\n"
          "by polynomials with selectable order.\n"  
          "The continuum is fitted be determining an off-set at a number of\n"
          "points (order+1) that are evenly spread between the lowest and\n"
          "highest frequency of f_mono."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( los_, absloswfs_, f_mono_, k_grid_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "order" ),
	TYPES(    int_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("kTempNoHydro"),
  	DESCRIPTION(
          "Calculates temperature 1D weighting functions WITHOUT including\n"
          "hydrostatic equilibrium."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( los_, absloswfs_, f_mono_, p_abs_, t_abs_, vmrs_, 
                  lines_per_tg_, abs_, trans_, e_ground_, t_ground_, k_grid_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kTempNoHydroNoGround"),
  	DESCRIPTION(
          "As kTempNoHydro but does not need any ground variables"),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( los_, absloswfs_, f_mono_, p_abs_, t_abs_, vmrs_, 
                   lines_per_tg_, abs_, trans_, k_grid_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kManual"),
  	DESCRIPTION(
          "Calculates a weighting function using y and y0.\n"
          "The weighting function is calculated as: k = (y-y0)/delta\n"
          "That is, delta is the magnitude of the perturbation done.\n"
          "The other variables are:\n"
          "  name     name on retrieval identity\n"
          "  grid     grid point value\n"
          "  apriori  a priori value"),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( y0_, y_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "name", "delta", "grid", "apriori" ),
	TYPES( string_t, Numeric_t, Numeric_t, Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kDiffHSmall"),
  	DESCRIPTION(
          "Calculates a weighting function using y, h1 and h2.\n"
          "This function minimizes memory usage. For faster calculations,\n"
          "use kDiffHFast.\n"
          "The weighting function is calculated as: k = (h2*y-h1*y)/delta\n"
          "That is, delta is the magnitude of the perturbation done.\n"
          "The other variables are:\n"
          "  name     name on retrieval identity\n"
          "  grid     grid point value\n"
          "  apriori  a priori value\n"
          "Note that the obtained k matrix includes effects of h1."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( h1_, h2_, y_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "name", "delta", "grid", "apriori" ),
	TYPES( string_t, Numeric_t, Numeric_t, Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kDiffHFast"),
  	DESCRIPTION(
          "Calculates a weighting function using y, h1 and h2.\n"
          "This function tries to be as fast as possible. To save memory,\n"
          "use kDiffHSmall.\n"
          "The weighting function is calculated as: k = (h2-h1)*y/delta\n"
          "See further kDiffHSmall."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( h1_, h2_, y_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "name", "delta", "grid", "apriori" ),
	TYPES( string_t, Numeric_t, Numeric_t, Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kxInit"),
  	DESCRIPTION(
          "Initializes Kx weighting function matrix and help variables\n"
          "(kx_names, kx_index and kx_aux).\n"
          "Use this function before the WF calculations are started or\n"
          "restarted."),
	OUTPUT( kx_, kx_names_, kx_index_, kx_aux_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kbInit"),
  	DESCRIPTION(
          "Initializes Kb weighting function matrix and help variables\n"
          "(kb_names, kb_index and kb_aux).\n"
          "Use this function before the WF calculations are started or\n"
          "restarted."),
	OUTPUT( kb_, kb_names_, kb_index_, kb_aux_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kxAppend"),
  	DESCRIPTION(
          "Appends the K matrix to Kx and handles additional data\n"
          "correspondingly."),
	OUTPUT( kx_, kx_names_, kx_index_, kx_aux_ ),
        INPUT( kx_, kx_names_, kx_index_, kx_aux_, k_, k_names_, k_aux_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kbAppend"),
  	DESCRIPTION(
          "Appends the K matrix to Kb and handles additional data\n"
          "correspondingly."),
	OUTPUT( kb_, kb_names_, kb_index_, kb_aux_ ),
        INPUT( kb_, kb_names_, kb_index_, kb_aux_, k_, k_names_, k_aux_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kxAppendUsingH"),
  	DESCRIPTION(
          "Applies a H matrix on the K matrix and appends the result to Kx.\n"
          "Additional data are treated correspondingly."),
	OUTPUT( kx_, kx_names_, kx_index_, kx_aux_ ),
        INPUT( kx_, kx_names_, kx_index_, kx_aux_, k_, k_names_, k_aux_ ),
	GOUTPUT(),
	GINPUT( Hmatrix_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kbAppendUsingH"),
  	DESCRIPTION(
          "Applies a H matrix on the K matrix and appends the result to Kb.\n"
          "Additional data are treated correspondingly."),
	OUTPUT( kb_, kb_names_, kb_index_, kb_aux_ ),
        INPUT( kb_, kb_names_, kb_index_, kb_aux_, k_, k_names_, k_aux_ ),
	GOUTPUT(),
	GINPUT( Hmatrix_ ),
	KEYWORDS(),
	TYPES()));



//======================================================================
//=== H Matrix Methods
//======================================================================

  md_data.push_back
    ( MdRecord
      ( NAME("VectorApplyH"),
  	DESCRIPTION(
          "Applies a H matrix on a vector."),
	OUTPUT(),
        INPUT(),
	GOUTPUT( VECTOR_ ),
	GINPUT( Hmatrix_, VECTOR_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixApplyH"),
  	DESCRIPTION(
          "Applies a H matrix on a matrix."),
	OUTPUT(),
        INPUT(),
	GOUTPUT( MATRIX_ ),
	GINPUT( Hmatrix_, MATRIX_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("Test"),
  	DESCRIPTION(
          "xxx."),
	OUTPUT(),
        INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

}
