/* Copyright (C) 2000, 2001 Stefan Buehler <sbuehler@uni-bremen.de>
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

void define_md_data()
{
  // The variable md_data is defined in file methods_aux.cc.
  extern Array<MdRecord> md_data;

  // Initialize to zero, just in case:
  md_data.resize(0);

  /* Here's a template record entry:  (PE 2001-09-18)

  md_data.push_back
    ( MdRecord
      ( NAME( "FunctionName" ),
	DESCRIPTION
        (
	 "A summary of the function in one sentence.\n"
         "\n"
         "A detailed description of the function. Please, try to be as \n"
         "clear and detailed as possible, this will help both you and \n"
         "others in the long run. \n"
         "   Additional paragraphs are indented with three blanks, as \n"
         "exemplified here.\n"
         "   The names of workspace variables and other methods\n"
         "are marked by stars, for example *z_plat*.\n"
         "   Generic input and output, and keywords shall be described \n"
         "as exemplified below. If there is no variables of a group, \n"
         "(e.g. generic input) remove that part totally. Note that the \n"
         "on-line help just gives the type of generic input/output and the \n"
         "keyword names, and additional information is for sure needed.\n"
         "   Leave space and brake lines when listing input and output \n"
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

  md_data.push_back
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


//======================================================================
//=== Overall ARTS functions
//======================================================================

  md_data.push_back
    ( MdRecord
      ( NAME( "Echo" ),
	DESCRIPTION
        (
	 "Outputs a string.\n"
         "\n"
         "Output the given message string, level follows the same convention as the -r\n"
         "command line flag of arts.\n"
         "\n"
         "Keywords:\n"
         "   message   : Message for output on screen.\n"
         "   level     : Output level for the message.\n"
        ),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "message", "level" ),
	TYPES(    String_t, Index_t    )));


  md_data.push_back     
    ( MdRecord
      ( NAME("Exit"),
	DESCRIPTION
	(
	 "Stops the execution and exits ARTS.\n"
	 "\n"
	 "This method is handy if you want to debug one of your\n"
	 "controlfiles. You can insert it anywhere in the controlfile. When\n"
	 "it is reached, it will terminate the program."
	),
	OUTPUT( ),
	INPUT( ),
	GOUTPUT( ),
	GINPUT( ),
	KEYWORDS( ),
	TYPES( )));

  md_data.push_back     
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




//======================================================================
//=== IO methods
//======================================================================

//=== Index ============================================================

  md_data.push_back     
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

  md_data.push_back
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

  md_data.push_back
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
  md_data.push_back
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

  md_data.push_back
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

  md_data.push_back
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

  md_data.push_back
    ( MdRecord
      ( NAME("NumericCopyFirstOfVector"),
	DESCRIPTION
        (
         "Sets a numeric workspace variable to the value of the first \n"
         "element of a vector. \n"
         "\n"
         "Generic output: \n"
         "   Numeric : The numeric variable to be set. \n"
         "\n"
         "Generic input:\n"
         "   Vector : The vector from which the value shall be obtained." 
        ),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Numeric_ ),
	GINPUT(  Vector_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("NumericCopyLastOfVector"),
	DESCRIPTION
        (
         "Sets a numeric workspace variable to the value of the last \n"
         "element of a vector. \n"
         "\n"
         "Generic output: \n"
         "   Numeric : The numeric variable to be set. \n"
         "\n"
         "Generic input:\n"
         "   Vector : The vector from which the value shall be obtained." 
        ),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Numeric_ ),
	GINPUT(  Vector_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
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

  md_data.push_back
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
  md_data.push_back
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

  md_data.push_back
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

  md_data.push_back
    ( MdRecord
      ( NAME("VectorCopy"),
        DESCRIPTION
        (
         "Creates a copy of a vector. \n"
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied. "
        ),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Vector_ ),
	GINPUT( Vector_ ),
	KEYWORDS(),
	TYPES()));


  md_data.push_back
    ( MdRecord
      ( NAME("VectorCopyFromMatrix"),
        DESCRIPTION
        (
         "Copies a row or a column from a matrix to a vector"
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Generic input: \n"
         "   Matrix : The source matrix. \n"
         "\n"
         "Keywords:\n"
         "   orientation : Could be either \"col\" or \"row\". \n"
         "   index       : Row or column number to be copied. \n"
        ),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Vector_ ),
	GINPUT(  Matrix_ ),
	KEYWORDS( "orientation", "index"   ),
	TYPES(    String_t,      Index_t   )));


  md_data.push_back
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

  md_data.push_back
    ( MdRecord
      ( NAME("VectorSetLengthFromVector"),
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

  md_data.push_back
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

  md_data.push_back
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

  md_data.push_back
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

  md_data.push_back
    ( MdRecord
      ( NAME("VectorPressuresForLinAltitudes"),
	DESCRIPTION
        (
         "Calculates a set of pressures corresponding to a set of\n"
         "linearly spaced altitudes. \n"
         "\n"
         "The linear set of altitudes is defined by an altitude step and \n"
         "a start and stop pressure. \n"
         "   The conversions between pressures and altitudes are based on\n"
         "*p_abs* and *z_abs*. \n"
         "\n"
         "Generic output: \n"
         "   Vector : Return vector for the pressure grid created. \n"
         "\n"
         "Keywords:\n"
         "     delta_z : Altitude step.\n"
         "     p_start : Start pressure.\n"
         "     p_stop  : Stop pressure."  
        ),
	OUTPUT(),
	INPUT( p_abs_, z_abs_ ),
	GOUTPUT( Vector_ ),
	GINPUT(),
	KEYWORDS( "delta_z", "p_start", "p_stop"  ),
	TYPES(    Numeric_t, Numeric_t, Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorFlip"),
	DESCRIPTION
        (
         "Creates a copy of a vector in reversed order. \n"
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied. "
        ),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Vector_ ),
	GINPUT( Vector_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
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
	OUTPUT(),
	INPUT(),
	GOUTPUT( Vector_ ),
	GINPUT( Matrix_, Vector_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
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

  md_data.push_back
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
  md_data.push_back
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

  md_data.push_back
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


  md_data.push_back
    ( MdRecord
      ( NAME("VectorPlanck"),
  	DESCRIPTION
        (
         "Sets a vector to the Planck function for the given frequencies\n"
         "and temperature. \n"
         "\n"
         "An example:\n"
         "   VectorPlanck(y_space,f_mono){temp=2.7} \n"
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : Input frequencies. \n"
         "\n"
         "Keywords:\n"
         "   temp : The blackbody temperature."  
        ),
	OUTPUT( ),
	INPUT( ),
	GOUTPUT( Vector_ ),
	GINPUT( Vector_ ),
	KEYWORDS( "temp"    ),
	TYPES(    Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorCalcLog10"),
  	DESCRIPTION
        (
         "Calculates the base 10 logarithm of a vector. \n"
         "\n"
         "The result can either be stored in the same or another vector. \n"
         "\n"
         "Generic output: \n"
         "   Vector : Return vector. \n"
         "\n"
         "Generic input: \n"
         "   Vector : Input vector. "
        ),
	OUTPUT( ),
	INPUT( ),
	GOUTPUT( Vector_ ),
	GINPUT( Vector_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorAdd"),
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

  md_data.push_back
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


//=== Matrix ==========================================================

  md_data.push_back
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

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixCopy"),
	DESCRIPTION
        (
         "Creates a copy of a matrix. \n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Matrix : The matrix to be copied. "
        ),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Matrix_ ),
	GINPUT( Matrix_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixFillWithVector"),
	DESCRIPTION
        (
         "Forms a matrix with n columns, and put the given vector in \n"
         "each column. \n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied. \n"
         "Keyword: \n"
         "   n : Number of columns in the matrix. "
        ),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Matrix_ ),
	GINPUT( Vector_ ),
	KEYWORDS( "n"   ),
	TYPES(    Index_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixMatrixMultiply"),
	DESCRIPTION
        (
	 "Multiply a Matrix with another Matrix and store the result in the result\n"
	 "Matrix.\n"
	 "\n"
	 "This just computes the normal Matrix-Matrix product, Y=M*X. It is ok\n"
	 "if Y and X are the same Matrix. This function is handy for\n"
	 "multiplying the H Matrix to weighting functions.\n"
	 "\n"
	 "Generic output:\n"
	 "   Matrix : The result of the multiplication (dimension mxc).\n"
	 "\n"
	 "Generic input:\n"
	 "   Matrix : The Matrix to multiply (dimension mxn).\n"
	 "   Matrix : The original Matrix (dimension nxc).\n"
        ),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Matrix_ ),
	GINPUT( Matrix_, Matrix_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixMatrixAdd"),
	DESCRIPTION
        (
	 "Adds two matrices. \n"
	 "\n"
	 "The function makes an element-wise addition. The size of the two \n"
         "matrices to add must have the same size. \n"
	 "\n"
	 "Generic output:\n"
	 "   Matrix : The result of the addition (dimension m x n).\n"
	 "\n"
	 "Generic input:\n"
	 "   Matrix : A matrix (dimension m x n).\n"
	 "   Matrix : A matrix (dimension m x n)."
        ),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Matrix_ ),
	GINPUT( Matrix_, Matrix_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
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

  md_data.push_back
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
  md_data.push_back
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

  md_data.push_back
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

  md_data.push_back
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

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixDiagonal"),
	DESCRIPTION
        (
         "Creates a diagonal matrix. \n"
         "\n"
         "All diagonal elements are set to the same value.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Keywords: \n"
         "   nrows : The number of rows (and columns) of the matrix to \n"
         "           create. \n"
         "   value : The value of the diagonal matrix elements. " 
        ),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Matrix_ ),
	GINPUT(),
	KEYWORDS( "nrows", "value"   ),
	TYPES(    Index_t,   Numeric_t )));



//=== ArrayOfIndex =====================================================

  md_data.push_back
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

  md_data.push_back
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
  md_data.push_back
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

  md_data.push_back
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

  md_data.push_back
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

  md_data.push_back
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
  md_data.push_back
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

  md_data.push_back
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

  md_data.push_back
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

  md_data.push_back
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
  md_data.push_back
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

  md_data.push_back
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

  md_data.push_back
    ( MdRecord
      ( NAME("StringSet"),
	DESCRIPTION("Sets a String to the given text String."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( String_ ),
	GINPUT(),
	KEYWORDS( "text"   ),
	TYPES(    String_t )));

  md_data.push_back
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

  md_data.push_back
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
  md_data.push_back
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

  md_data.push_back
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

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfStringSet"),
	DESCRIPTION("Sets a String array according the given text.\n"
                    "The format is text = [\"String1\",\"String2\",...]"),
	OUTPUT(),
	INPUT(),
	GOUTPUT( ArrayOfString_ ),
	GINPUT(),
	KEYWORDS( "text"         ),
	TYPES(    Array_String_t )));

  md_data.push_back
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

  md_data.push_back
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
  md_data.push_back
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

  md_data.push_back
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



//=== LOS ==================================================================

#ifdef HDF_SUPPORT
  md_data.push_back
    ( MdRecord
      ( NAME("LosWriteBinary"),
	DESCRIPTION(
		    "Writes a LOS structure to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   LOS : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
		    ),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( Los_ ),
	KEYWORDS( "filename" ),
	TYPES(    String_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("LosReadBinary"),
	DESCRIPTION(
		    "Reads a LOS structure from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   LOS : Name of the workspace variable to read.\n"
                    "\n"
		    "Keywords:\n"
		    "   filename : Name of the input file.\n"
		    ),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Los_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    String_t   )));
#endif // HDF_SUPPPORT



//======================================================================
//=== Absorption methods
//======================================================================

//=== Spectroscopic methods ============================================

  md_data.push_back
    ( MdRecord
      ( NAME("lines_per_tgSetEmpty"),
  	DESCRIPTION
	(
	 "Sets lines_per_tg to empty line lists.\n"
	 "\n"
	 "You can use this method to set lines per tag if you do not reall want\n"
	 "to compute line spectra. Formally, absCalc will still require\n"
	 "lines_per_tg to be set.\n"
	 ),
	OUTPUT(   lines_per_tg_      ),
	INPUT(    tgs_        ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(  ),
	TYPES(    )));
  
  md_data.push_back
    ( MdRecord
      ( NAME("lines_per_tgReadFromCatalogues"),
  	DESCRIPTION(
		    "This method can read lines from different line \n"
		    "catalogues.\n"
		    "\n"
		    "For each tag group, you can specify which catalogue\n"
		    "to use. Because the method creates lines_per_tg directly,\n"
		    "it replaces for example thefollowing two method calls:\n"
		    "  - linesReadFromHitran\n"
		    "  - lines_per_tgCreateFromLines\n"
		    "   This method needs as input WSVs the list of tag \n"
		    "groups. Keyword parameters must specify the names of\n"
		    "the catalogue files to use and the matching formats.\n"
		    "Names can be anything, formats can currently be \n"
		    "HITRAN96, MYTRAN2, JPL, or ARTS. Furthermore, keyword\n"
		    "parameters have to specify minimum and maximum \n"
		    "frequency for each tag group. To safe typing, if there\n"
		    "are less elements in the keyword parameters than there\n"
		    "are tag groups, the last parameters are applied to all\n"
		    "following tag groups.\n"
		    "\n"
		    "Example usage:\n"
		    "\n"
		    "lines_per_tgReadFromCatalogues{\n"
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
		    "(linesReadFromHitran, etc.), as well as\n"
		    "lines_per_tgCreateFromLines.\n"
		    "\n"
		    "Keywords: \n"
		    "   filenames = Name (and path) of the catalogue files.\n"
		    "   formats   = allowed formats are HITRAN96,MYTRAN2,JPL,ARTS \n"
		    "   fmin      = Minimum frequency for lines to read in Hz.\n"
		    "   fmax      = Maximum frequency for lines to read in Hz.\n"),
	OUTPUT(   lines_per_tg_      ),
	INPUT(    tgs_        ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "filenames",    "formats",      "fmin",   "fmax" ),
	TYPES(    Array_String_t, Array_String_t, Vector_t, Vector_t)));
  
  md_data.push_back
    ( MdRecord
      ( NAME("linesReadFromHitran"),
  	DESCRIPTION(
		    "Read all the lines from a HITRAN catalogue file in the \n"
		    "given frequency range. Otherwise a runtime error will be\n"
		    "thrown\n"
		    "\n"
		    "Please note that all lines must correspond\n"
		    "to the legal species / isotope combinations\n"
		    "\n"
		    "Keywords: \n"
		    "   filename = Name (and path) of the catalogue file.\n"
		    "   fmin     = Minimum frequency for lines to read in Hz.\n"
		    "   fmax     = Maximum frequency for lines to read in Hz."),
	OUTPUT(   lines_   ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "filename",  "fmin",    "fmax"),
	TYPES(    String_t,    Numeric_t, Numeric_t)));
  
  md_data.push_back
    ( MdRecord
      ( NAME("linesReadFromMytran2"),
  	DESCRIPTION(
		    "Read all the lines from a MYTRAN2 catalogue file in the \n"
		    "given frequency range. Otherwise a runtime error will be\n"
		    "thrown\n"
		    "\n"
		    "Please note that all lines must correspond\n"
		    "to the legal species / isotope combinations\n"
		    "\n"
		    "Keywords: \n"
		    "   filename = Name (and path) of the catalogue file.\n"
		    "   fmin     = Minimum frequency for lines to read in Hz.\n"
		    "   fmax     = Maximum frequency for lines to read in Hz."),
	OUTPUT(   lines_   ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "filename",  "fmin",    "fmax"),
	TYPES(    String_t,    Numeric_t, Numeric_t)));
  
  md_data.push_back
    ( MdRecord
      ( NAME("linesReadFromJpl"),
  	DESCRIPTION(
		    "Read all the lines from a JPL catalogue file in the \n"
		    "given frequency range. Otherwise a runtime error will be\n"
		    "thrown\n"
		    "\n"
		    "Please note that all lines must correspond\n"
		    "to the legal species / isotope combinations.\n"
		    "\n"
		    "Keywords: \n"
		    "   filename = Name (and path) of the catalogue file.\n"
		    "   fmin     = Minimum frequency for lines to read in Hz.\n"
		    "   fmax     = Maximum frequency for lines to read in Hz."),
	OUTPUT(   lines_   ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "filename",  "fmin",    "fmax"),
	TYPES(    String_t,    Numeric_t, Numeric_t)));
  
  md_data.push_back
    ( MdRecord
      ( NAME("linesReadFromArts"),
  	DESCRIPTION(
		    "Read all the lines from an Arts catalogue file in the \n"
		    "given frequency range. Otherwise a runtime error will be\n"
		    "thrown \n"
		    "\n"
		    "Please note that all lines must correspond\n"
		    "to the legal species / isotope combinations\n"
		    "\n"
		    "Keywords: \n"
		    "   filename = Name (and path) of the catalogue file.\n"
		    "   fmin     = Minimum frequency for lines to read in Hz.\n"
		    "   fmax     = Maximum frequency for lines to read in Hz."),
	OUTPUT(   lines_   ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "filename",  "fmin",    "fmax"),
	TYPES(    String_t,    Numeric_t, Numeric_t)));
  
  // FIXME: Remove this one.
  md_data.push_back
    ( MdRecord
      ( NAME("linesElowToJoule"),
  	DESCRIPTION(
		    "Just a little helper to convert the lower state energy from cm^-1\n"
		    "(ARTSCAT-2) to Joule (ARTSCAT-3). This should be removed soon\n"),
	OUTPUT(   lines_   ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( ),
	TYPES(    )));
      
  md_data.push_back
    ( MdRecord
      ( NAME("lines_per_tgCreateFromLines"),
  	DESCRIPTION(
		    "Split lines up into the different tag groups.\n"
		    "\n"
		    "The tag groups are tested in the order in which they are\n" 
		    "specified in the controlfile. The lines are assigned to \n"
		    "the tag groups in the order as the groups  are specified.\n"
		    "That means if you do [\"O3-666\",\"O3\"],the last group O3 \n"
		    "gets assigned all the O3 lines that do not fit in the first group."),
	OUTPUT(   lines_per_tg_      ),
	INPUT(    lines_, tgs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));
  
  md_data.push_back
    ( MdRecord
      ( NAME("lines_per_tgAddMirrorLines"),
  	DESCRIPTION(
		    "Adds mirror lines at negative frequencies to the *lines_per_tg*.\n"
		    "\n"
		    "For each line at frequency +f in *lines_per_tg* a corresponding\n"
		    "entry at frequency -f is added to *lines_per_tg*.The mirror \n"
		    "lines are appended to the line lists after the original lines."),
	OUTPUT(   lines_per_tg_      ),
	INPUT(    lines_per_tg_      ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));
  
  md_data.push_back
    ( MdRecord
      ( NAME("lines_per_tgCompact"),
  	DESCRIPTION(
		    "Removes all lines outside the defined lineshape cutoff frequency\n"
		    "from the *lines_per_tg*. This can save computation time.\n"
		    "It should be particularly useful to call this method after\n"
		    "*lines_per_tgAddMirrorLines*."),
	OUTPUT(   lines_per_tg_      ),
	INPUT(    lines_per_tg_, lineshape_, f_mono_  ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));
  
  md_data.push_back
    ( MdRecord
      ( NAME("linesWriteAscii"),
  	DESCRIPTION(
                    "Writes the workspace variable *lines* to an ASCII file.\n"
                    "\n"
		    "The content of the workspace variable 'lines`\n"
		    "The content of the workspace variable *lines*\n"
		    "is written in ARTS line format to the file with\n"
                    "the specified name. If the filename is omitted, the\n"
                    "lines are written to <basename>.lines.aa.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ), 
	OUTPUT(),
	INPUT( lines_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    String_t   )));
  
  md_data.push_back
    ( MdRecord
      ( NAME("lines_per_tgWriteAscii"),
  	DESCRIPTION(
                    "Writes the workspace variable *lines_per_tg* to an ASCII file.\n"
                    "\n"
                    "The content of the workspace variable *lines_per_tg*\n"
                    "is written in ARTS line format to the file with\n"
                    "the specified name. If the filename is omitted, the\n"
                    "lines are written to <basename>.lines_per_tg.aa.\n"
                    "\n"
                    "The array dimension is handled in a similar way as by the\n"
                    "array of vector and matrix output functions:\n"
                    "First an integer stating the number of tag groups.\n"
                    "Then an integer specifying the number of lines for the\n"
                    "first group. Then the other groups in similar fashion."
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
	OUTPUT(),
	INPUT( lines_per_tg_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    String_t   )));
  
  md_data.push_back
    ( MdRecord
      ( NAME("tgsDefine"),
  	DESCRIPTION(
		    "Set up the list of tag groups.\n"
		    "\n"
		    "The workspace variable *tgs* contains several tag groups. Each \n"
		    "tag group contain one or more tags. This method converts \n"
		    "description of tag groups  given in the keyword to the internal \n"
		    "representation *tgs*. A tag group selects spectral features which \n"
		    "belong to the same species. \n"
		    "   A tag group can contain a mixture of general and special \n"
		    "tags.  All the continuum tags belong to the special tags and \n"
		    "the rest come under the general tags.\n"
		    "   A general tag is defined in terms of the name of the species,\n"
		    "isotope and a range of frequencies. Species are named after the \n"
		    "standard chemical names,e.g., \"O3\".  Isotopes are given by the \n"
		    "last digit of the atomic weight, i.e., \"O3-668\" for the \n"
		    "asymmetric ozone molecule including an oxygen 18 atom.  Groups\n"
		    "of transitions are specified by giving a lower and upper limit \n"
		    "of a frequency range,\"O3-666-500e9-501e9\".Moreover the symbol\n"
		    "'*' acts as a wild card. Furthermore, frequency range or frequency\n"
		    "range and isotope may be omitted.\n"
		    "Example for some tag groups containing only general tags:\n"
		    "tags = [\"O3-666-500e9-501e9, O3-686\",\"O3\"]\n"
		    "The first tag group consist of all O3-666 lines between 500 and\n"
		    "501 GHz plus all O3-686 lines.  The second tag group will contain\n"
		    "all remaining O3 transitions.\n"
		    "\n"
		    "Keywords:\n"
		    "   tags : Specify one String for each tag group that you want to create.\n"
		    "   Inside the String, separate the tags by comma (plus optional blanks).\n"
		    "   Example:\n"
		    "   tag = [\"O3-686\",\"H2O\"]"),
	OUTPUT( tgs_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "tags" ),
	TYPES(    Array_String_t   )));
  
  md_data.push_back
    ( MdRecord
      ( NAME("tgsDefineAllInScenario"),
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
	 "              For example: /pool/lookup2/arts-data/atmosphere/fascod/tropical"
           ),
	OUTPUT( tgs_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "basename" ),
	TYPES(    String_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("lineshapeDefine"),
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
	OUTPUT( lineshape_ ),
	INPUT( tgs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(  "shape",    "normalizationfactor",  "cutoff" ),
	TYPES(     String_t,        String_t,         Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("lineshape_per_tgDefine"),
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
	OUTPUT( lineshape_ ),
	INPUT( tgs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(  "shape",           "normalizationfactor",    "cutoff" ),
	TYPES(   Array_String_t,         Array_String_t,        Vector_t )));


//=== Continuum methods ============================================

  md_data.push_back
    ( MdRecord
      ( NAME("cont_descriptionInit"),
  	DESCRIPTION
	(
	 "Initializes the two workspace variables for the continuum description,\n"
	 "*cont_description_names* and *cont_description_parameters*.\n"
	 " \n"
	 "This method does not really do anything, except setting the two\n"
	 "variables to empty Arrays. It is just necessary because the method\n"
	 "*cont_descriptionAppend* wants to append to the variables.\n"
	 "   Formally, the continuum description workspace variables are required\n"
	 "by the absorption calculation methods (e.g., *absCalc*). Therefore you\n"
	 "always have to call at least *cont_descriptionInit*, even if you do\n"
	 "not want to use any continua."
	 ),
	OUTPUT( cont_description_names_, 
                cont_description_models_,
                cont_description_parameters_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("cont_descriptionAppend"),
  	DESCRIPTION
	(
	 "Appends the description of a continuum model or a complete absorption\n"
	 "model to *cont_description_names* and *cont_description_parameters*.\n"
	 "\n"
	 "See online documentation for *cont_description_names* for a list of\n"
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
	OUTPUT( cont_description_names_, 
                cont_description_models_,
                cont_description_parameters_ ),
	INPUT(  cont_description_names_, 
                cont_description_models_,
                cont_description_parameters_),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "tagname",  "model",   "userparameters" ),
	TYPES(    String_t,   String_t,   Vector_t         )));


//=== Input Atmosphere methods ===========================================

  md_data.push_back
    ( MdRecord
      ( NAME("raw_vmrsReadFromFiles"),
        DESCRIPTION(
          "Reads the individual VMR profile for each TAGS from file.\n"
          "\n"
          "Using this function one can read VMRs of specific TAGS from\n"
          "explicitly specified files and the remaing from a scenario.\n"
          "The filenames and the base name of atmospheric scenario\n"
          "should be specified as keywords. One file name must\n"
          "be specified for each tag group(each element of *tgs*).\n"
          "The name may include a path.\n"
	  "\n"
	  "Keywords:\n"
	  "   seltags   : Must be a sub group of tags which should be read from files.\n"
	  "   filenames : Names of the files containing VMR profiles of seltags.\n"
	  "   basename  : The name of a particular atmospheric scenario.\n"
	  "               See *raw_vmrsReadFromScenario* for details. Remaining\n"
	  "               VMRs will be read from the scenario.\n"
	  "\n"
          ),
        OUTPUT(   raw_vmrs_         ),
        INPUT(    tgs_                 ),
        GOUTPUT(                       ),
        GINPUT(                        ),
        KEYWORDS( "seltags",       "filenames",    "basename"),
        TYPES(    Array_String_t,  Array_String_t, String_t)));

  md_data.push_back
    ( MdRecord
      ( NAME("raw_vmrsReadFromScenario"),
  	DESCRIPTION(
	  "Reads the individual VMR profile for each tag group from a standard\n"
	  "atmospheric scenario.\n" 
	  "\n"
          "Five different atmospheric scenarios are available in arts data:\n"
          "tropical, midlatitude-summer, midlatitude-winter, subartic-summer\n"
          "and subartic-winter.\n"
	  "\n"
	  "   Files in the scenarios look like this: tropical.H2O.aa\n"
	  "\n"
	  "   The basename must include the path, i.e., the files can be anywhere,\n"
	  "but they must be all in the same directory.\n"
	  "   The profile is chosen by the species name. If you have more than one\n"
	  "tag group for the same species, the same profile will be used.\n"
	  "\n"
	  "Keywords:\n"
	  "   basename :The name and path of a particular atmospheric scenario.\n"
	  "   For example:\n"
	  "   /pool/lookup2/arts-data/atmosphere/fascod/tropical\n"
	  "\n"
	  ),
	OUTPUT(   raw_vmrs_    ),
	INPUT(    tgs_                 ),
	GOUTPUT(                       ),
	GINPUT(                        ),
	KEYWORDS( "basename"           ),
	TYPES(    String_t             )));

  md_data.push_back
    ( MdRecord
      ( NAME("AtmFromRaw"),
  	DESCRIPTION(
	  "Interpolates temperature, altitude, and VMRs to the pressure grid\n"
	  "given by p_abs.\n" 
	  "\n"
          "The altitude is not used by the absorption routines,\n"
	  "but later on by the RT routines.\n"
	  "\n"
	  "Interpolations used: \n"
	  "\n"
	  "Temperature      : Linear interpolation in ln(p)\n"
	  "Altitude         : Linear interpolation in ln(p)\n"
	  "VMRs             : Linear interpolation in ln(p)\n"
	  "Cloud Parameters : Linear interpolation in ln(p)\n"
	  "\n"
	  ),
	OUTPUT(   t_abs_    , z_abs_   , vmrs_           ),
	INPUT(    tgs_, p_abs_    , raw_ptz_ , raw_vmrs_ ),
	GOUTPUT(                       			 ),         
	GINPUT(                        			 ),
	KEYWORDS(                             		 ),
	TYPES(                          		 )));

  md_data.push_back
    ( MdRecord
      ( NAME("WaterVaporSaturationInClouds"),
  	DESCRIPTION(
	  "Calculates the water vapor saturation volume mixing ratio (VMR) in the\n"
	  "vertical range where liquid or ice clouds are in the atmosphere.\n"
	  "At the pressure/altitude grid points where the liquid water content (LWC)\n"
	  "or ice water content (IWC) of the clouds (tags 'liquidcloud' and 'icecloud')\n"
          "is larger than zero the H2O-VMR is set to liquid water/ice saturation VMR.\n"
          "The saturation pressure is calculated according to Goff-Gratch equations.\n"
	  ),
	OUTPUT(   vmrs_ , p_abs_                         ),
	INPUT(    vmrs_ , p_abs_ , t_abs_ , tgs_         ),
	GOUTPUT(                       			 ),         
	GINPUT(                        			 ),
	KEYWORDS(                             		 ),
	TYPES(                          		 )));

  md_data.push_back
    ( MdRecord
      ( NAME("vmrsScale"),
	DESCRIPTION(
          "Scales the vmr input of the tgs given in scaltgs by the\n"
	  "factors given in scalfac.\n"
	  "\n"
	  "Keywords:\n"
	  "   scaltgs : subgroup of tags which has to be scaled.\n"
	  "   scalfac : the factor with which vmr to be scaled.\n"
	  "\n"
	  ),
	OUTPUT(	vmrs_ ),
	INPUT( 	tgs_, vmrs_  ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "scaltgs", "scalfac"),
	TYPES( Array_String_t, Vector_t)));

  md_data.push_back
    ( MdRecord
      ( NAME("h2o_absSet"),
	DESCRIPTION(
          "Sets h2o_abs to the profile of the first tag group containing\n"
	  "water.\n" 
	  "\n"
          "This is necessary, because for example *absCalc* requires h2o_abs\n"
	  "to contain the water vapour profile(the reason for this is the\n"
          "calculation of oxygen line brodening requires water vapour profile).\n"
	  "Then this function can be used to copy the profile of the first tag\n"
          "group of water.\n"
	  "\n"
	  ),
	OUTPUT(	h2o_abs_ ),
	INPUT( 	tgs_, vmrs_  ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("n2_absSet"),
	DESCRIPTION(
          "Sets n2_abs to the profile of the first tag group containing\n"
	  "molecular nitrogen. See *h2o_absSet* for more details.\n"
	  "\n"
	  ),
	OUTPUT(	    n2_abs_ ),
	INPUT( 	tgs_, vmrs_  ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("hseSet"),
	DESCRIPTION(
          "Sets the vector of parameters for calculation of hydrostatic \n"
          "equilibrium (*hse*). The on/off flag is set to 1. \n"
          "\n"
          "Type \"arts -d hse\" for more information. \n"
          "\n"
          "Keywords \n"
          "  pref  : Pressure of the reference point. \n"
          "  zref  : The geometrical altitude at pref. \n"
          "  g0    : Gravitational acceleration at the geoid surface.\n"
          "  niter : Number of iterations (1-2 should suffice normally)."),
	OUTPUT( hse_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "pref",    "zref",    "g0",      "niter" ),
	TYPES(    Numeric_t, Numeric_t, Numeric_t, Index_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("hseSetFromLatitude"),
	DESCRIPTION(
          "Sets the vector of parameters for calculation of hydrostatic \n"
          "equilibrium (*hse*). The on/off flag is set to 1. The gravitational\n"
          "acceleration is calculated following the international gravity formula.\n"
          "\n"
          "Type \"arts -d hse\" for more information. \n"
          "\n"
          "Keywords \n"
          "  pref     : Pressure of the reference point. \n"
          "  zref     : The geometrical altitude at pref. \n"
          "  latitude : Geocentric latitude of observation point (-90 to 90 Degree).\n"
          "  niter    : Number of iterations (1-2 should suffice normally)."),
	OUTPUT( hse_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "pref",    "zref",    "latitude", "niter" ),
	TYPES(    Numeric_t, Numeric_t, Numeric_t,  Index_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("hseSetFromLatitudeIndex"),
	DESCRIPTION(
          "As hseSetFromLatitude, but sets pref, zref to the values given by index \n"
          "from vectors p_abs, z_abs (e.g. 0 = ground).\n"
          "\n"
          "Type \"arts -d hse\" for more information. \n"
          "\n"
          "Keywords \n"
          "  latitude : Geocentric latitude of observation point (-90 to 90 Degree).\n"
          "  index    : Reference index within p_abs, z_abs for setting pref, zref.\n"
          "  niter    : Number of iterations (1-2 should suffice normally)."),
	OUTPUT( hse_ ),
	INPUT(p_abs_, z_abs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "latitude", "index", "niter" ),
	TYPES(    Numeric_t,  Index_t, Index_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("hseFromBottom"),
	DESCRIPTION(
          "As hseSet but uses the first values of p_abs and z_abs for pref\n"
          "and zref, respectively.\n"
          "\n"
          "Keywords \n"
          "  g0    : Gravitational acceleration at the geoid surface.\n"
          "  niter : Number of iterations (1-2 should suffice normally)."),
	OUTPUT( hse_ ),
	INPUT( p_abs_, z_abs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "g0",      "niter" ),
	TYPES(    Numeric_t, Index_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("hseOff"),
	DESCRIPTION(
          "Turns off hydrostatic equilibrium. \n"
          "\n"
          "The on/off flag off *hse* is set to 0 and *hse* is set to be a \n"
          "vector of length 1."),
	OUTPUT( hse_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("hseCalc"),
	DESCRIPTION(
          "Ensures that 'z_abs' fulfills hydrostatic equilibrium. \n"
          "\n"
          "Nothing is done if the on/off flag of *hse* is set to 0. The \n"
          "reference point, g at the ground and number of iterations \n"
          "are taken from *hse*. \n"
          "   The given altitudes (*z_abs*) are used as a first guess when \n"
          "starting the calculations (to estimate g etc.). The altitude \n"
          "variation of the gravitational acceleration is considered. The \n"
          "average molecular weight is assumed to be 28.96 at all altitudes.\n"
          "The amount of water vapour is taken into account. \n"
          "    The calculations are repeated according to the number of \n"
          "iterations specified. A higher number of iterations \n" 
          "improves the accuracy, but one iteration should be normally \n"
          "enough if *z_abs* already has reasonable values. Two iterations \n"
          "should suffice for basically all applications."),
	OUTPUT( z_abs_ ),
	INPUT( z_abs_, p_abs_, t_abs_, h2o_abs_, r_geoid_, hse_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));



//=== Absorption methods ===============================================

  md_data.push_back
    ( MdRecord
      ( NAME( "absCalc" ),
	DESCRIPTION(
	   "Calculate absorption coefficients. \n"
	   "\n"
	   "This function calculates both, the total absorption (*abs*)\n"
	   "and the absorption per tag group (*abs_per_tg*).\n"
            ) ,
	OUTPUT(abs_  , abs_per_tg_ ),
	INPUT(tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_, 
              lines_per_tg_, lineshape_,
	      cont_description_names_, cont_description_models_, 
              cont_description_parameters_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME( "absCalcSaveMemory" ),
	DESCRIPTION(
	   "Calculate absorption coefficients, trying to conserve memory. \n"
	   "\n"
	   "This function calculates only the total absorption (*abs*),\n"
	   "NOT the absorption per tag group (*abs_per_tg*).\n"
           "\n"
           "This means you cannot use it if you want to calculate Jacobians\n"
           "later.\n"
           "\n"
           "The implementation follows absCalc."
            ) ,
	OUTPUT(abs_ ),
	INPUT(tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_, 
              lines_per_tg_, lineshape_,
	      cont_description_names_, cont_description_models_, 
              cont_description_parameters_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("absCalcFromXsec"),
	DESCRIPTION(
		    "Calculate absorption coefficients from cross sections.\n"
		    "\n"
		    "This calculates both the total absorption and the\n"
		    "absorption per tag group. \n"
		    "This method calls three other  methods:\n"
		    "1. *xsec_per_tgInit* - initialize *xsec_per_tg* \n"
		    "2. *xsec_per_tgAddLine* - calculate cross sections per \n"
		    "                   tag group for line spectra.\n"
		    "3. *xsec_per_tgAddConts* - calculate cross sections per \n"
		    "                   tag group for continua.\n"
		    "Then it calculates the absorption coefficient by multiplying\n"
		    "the cross section by VMR.\n"
                    "This is done once for each tag group (output: *abs_per_tg*)\n"
		    "and for the sum of all tag group to get the total absorption\n"
		    "coefficient (output: *abs*)\n"
		    ),
	OUTPUT(	    abs_  , abs_per_tg_ ),
	INPUT( 	    xsec_per_tg_, vmrs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME( "xsec_per_tgInit" ),
	DESCRIPTION(
	   "Initialize *xsec_per_tg*.\n"
	   "\n"
	   "The initialization is\n"
	   "necessary, because methods *xsec_per_tgAddLines*\n"
	   "and *xsec_per_tgAddConts* just add to *xsec_per_tg*.\n"
	   "The size is determined from *tgs*.\n"
	   ),
	OUTPUT( xsec_per_tg_ ),
	INPUT(tgs_, f_mono_, p_abs_),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("xsec_per_tgAddLines"),
	DESCRIPTION(
		    "Calculate cross sections per tag group for line spectra.\n"
		   ),
	OUTPUT(	    xsec_per_tg_                             ),
	INPUT( 	    tgs_, f_mono_, p_abs_, t_abs_, h2o_abs_, vmrs_, 
		    lines_per_tg_, lineshape_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("xsec_per_tgAddConts"),
	DESCRIPTION(
		    "Calculate cross sections per tag group for continua.\n"
                     ),
	OUTPUT(	    xsec_per_tg_                             ),
	INPUT( 	    tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_,
		    cont_description_names_, cont_description_parameters_,
                    cont_description_models_),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));


//=== Methods operating on absorption ========================================

  md_data.push_back
    ( MdRecord
      ( NAME("abs_per_tgReduce"),
	DESCRIPTION(
		    "Reduces absorption coefficients. Only absorption\n"
		    "coefficients for which weighting functions are\n"
		    "calculated are kept in memory.\n"
		    ),
	OUTPUT(	    abs_per_tg_ ),
	INPUT( 	    abs_per_tg_, tgs_, wfs_tgs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));


//=== Refraction ==========================================================

  md_data.push_back
    ( MdRecord
      ( NAME("refrSet"),
	DESCRIPTION(
           "Sets the refraction input arguments (refr, refr_model and \n"
           "refr_lfac) to the specified values. \n"
           "\n"
           "Type \"arts -d refr\" etc. for more information on the input \n"
           "arguments. See *refrCalc* for avaliable refraction models.\n"
           "\n"
           "Keywords:\n"
           "     on    : On/off boolean.\n"
           "     model : Name on parametization for the refractive index.\n"
           "     lfac  : Length factor for ray tracing."  ),
	OUTPUT( refr_, refr_lfac_, refr_model_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "on",  "model",  "lfac"    ),
	TYPES(    Index_t, String_t, Index_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("refrOff"),
  	DESCRIPTION(
          "Sets the refraction boolean (*refr*) to zero and gives the other \n"
          "refraction input arguments (*refr_lfac* and *refr_model*) dummy \n"
          "values (that will give error messages if used)."),
	OUTPUT( refr_, refr_lfac_, refr_model_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("refrCalc"),
	DESCRIPTION(
           "Calculates the refractive index using the parameterization\n"
           "specified by *refr_model*. \n"
           "\n"
           "If *refr* is set to zero, the refractive index is set to be an \n"
           "empty vector. \n"
           "\n"
           "Existing parameterizations are: \n"
           "\n"
           "   'Unity': \n"
           "      Sets the refractive index to 1 at all altitudes. \n"
           "\n"
           "   'Boudouris': \n"
           "      Refractive index at microwave frequencies following \n"
           "      Boudouris 1963. The k-parameter values were taken from \n"
           "      Section 5.1.1 of the Janssen book. The Z parameters are \n"
           "      set to 1. \n"
           "\n"
           "  'BoudourisDryAir': \n"
           "      As Boudouris but setting the water content to zero. "),
	OUTPUT(	refr_index_ ),
	INPUT( 	p_abs_, t_abs_, h2o_abs_, refr_, refr_model_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));



//======================================================================
//=== LOS/RTE methods
//======================================================================

  md_data.push_back
    ( MdRecord
      ( NAME("zaFromZtan"),
	DESCRIPTION(
           "Calculates the zenith angles corresponding to a set of tangent\n"
           "altitudes.\n"
           "\n"
           "Refraction is considered if it is turned on (refr=1)."),
	OUTPUT(),
	INPUT( z_tan_, z_plat_ , p_abs_, z_abs_, refr_, refr_index_, r_geoid_,
               z_ground_ ),
	GOUTPUT( Vector_ ),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME( "zaFromDeltat" ),
	DESCRIPTION(
	   "Calculates the zenith angles for a LEO-LEO cross-link.\n"
           "\n"
           "The function calculates the zenith angles corresponding to a \n"
           "LEO-LEO cross-link for an occultation between two altitudes and \n"
           "the given time increment. The LEOs are supposed to be moving \n"
           "in opposite directions in identical orbits (but not colliding!).\n"
           "   The time window where the LEOs are interacting is defined by \n"
           "a start and stop tangent altitudes.\n"
           "   The function uses REFR to determine if refraction shall be \n"
           "considered or not. \n"
           "\n"
           "Generic output: \n"
           "   Vector : Return vector for the zenith angles. The normal \n"
           "            options are ZA_PENCIL and ZA_SENSOR. \n"
           "\n"
           "Keywords:\n"
           "   delta_t   : Time increment between observations.\n"
           "   z_tan_lim : Vector with start and stop tangent altitudes." ),
	OUTPUT(),
	INPUT( z_plat_ , p_abs_, z_abs_, l_step_, refr_, refr_lfac_, 
               refr_index_, r_geoid_, z_ground_ ),
	GOUTPUT( Vector_ ),
	GINPUT(),
	KEYWORDS( "delta_t", "z_tan_lim" ),
	TYPES(    Numeric_t, Vector_t    )));

  md_data.push_back
    ( MdRecord
      ( NAME("r_geoidStd"),
  	DESCRIPTION(
          "Sets the geoid radius to the standard Earth radius defined in \n"
          "constants.cc."),
	OUTPUT( r_geoid_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("r_geoidWGS84"),
  	DESCRIPTION(
          "Sets the geoid radius according to WGS-84. \n"
          "\n"
          "The function is based on Section 9.4.1 in the Rodgers book. \n"
          "The observation direction is given as the angle to the meridian \n"           "plane (that is, S=N=0, W=E=90).\n"
          "\n"
          "Keywords \n"
          "      latitude : Latitude at the measurement.\n"
          "  obsdirection : Observation direction (see above)."),
	OUTPUT( r_geoid_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "latitude", "obsdirection" ),
	TYPES(    Numeric_t,  Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("groundOff"),
  	DESCRIPTION(
          "Sets dummy values to the ground variables. \n"
          "\n"
          "The ground altitude is set to the first element of 'z_abs'.\n"
          "The ground temperature (t_ground) is set to 0. \n"
          "The ground emission vector (e_ground) is set to be empty.\n"
          "   If there is a ground intersection and only this function is\n"
          "used to set the ground variables, there will be error messages."),
	OUTPUT( z_ground_, t_ground_, e_ground_ ),
	INPUT( z_abs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("groundSet"),
  	DESCRIPTION(
          "Sets the ground altitude and emission to the specified values,\n"
          "and selects a ground temperature.\n"
          "\n"
          "The emission is set to be identical for all frequencies. \n"
          "The ground temperature is obtained by interpolating *t_abs*.\n"
          "\n"
          "Keywords \n"
          "  z : Altitude above the geoid of the ground.\n"
          "  e : Ground emission factor."),
	OUTPUT( z_ground_, t_ground_, e_ground_ ),
	INPUT( p_abs_, t_abs_, z_abs_, f_mono_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "z",       "e"       ),
	TYPES(    Numeric_t, Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("groundAtBottom"),
  	DESCRIPTION(
          "Sets the ground emission to the specified value, and sets ground \n"
          "altitude and temperature to the first values of z_abs and t_abs.\n"
          "\n"
          "The emission is set to be identical for all frequencies. \n"
          "\n"
          "Keywords \n"
          "  e : Ground emission factor."),
	OUTPUT( z_ground_, t_ground_, e_ground_ ),
	INPUT( t_abs_, z_abs_, f_mono_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "e"       ),
	TYPES(    Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("groundFlatSea"),
  	DESCRIPTION(
          "Models the emission from a flat sea. \n"
          "\n"
          "The method sets the ground variables to match the properties of \n"
          "the sea, without wind effects. The emissivity is calculated from\n"
          "the dielectric constant by the Fresnel equations.\n"
          "\n"
          "The altitude is set to 0 m. The skin temperature is set by the\n"
          "keyword argument *t_skin*. If this argument is set to be <= 0,\n"
          "the temperature is obtained by interpolating *t_abs*. \n"
          "\n"
          "The incident angle (of reflection) is calculated for max of \n"
          "*za_pencil*. Refraction is considered or not, depending on value \n"
          "of *refr*. This means that *refrCalc* must be called before this \n"
          "method if refraction is considered. The emissivity depends on the\n"
          "selected polarisation. The refractive index of air is set to 1.\n"
          "\n"
          "The relative dielectric constant is calculated following\n"
          "Liebe et al. 1991 Int. J. IR+mm Waves 12(12), 659-675. \n"
          "The method does not consider salinity and is restricted to the\n"
          "range 5 - 1000 GHz (below 5 GHz salinity must be considered).\n"
          "\n"
          "Keywords \n"
          "  pol    : Polarisation. Calid options are \"v\" or \"h\".\n"
          "  t_skin : Skin temperature. "),
	OUTPUT( z_ground_, t_ground_, e_ground_ ),
	INPUT( p_abs_, t_abs_, z_abs_, f_mono_, za_pencil_, z_plat_, r_geoid_,
           refr_, refr_index_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "pol",    "t_skin"  ),
	TYPES(    String_t, Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("emissionOn"),
  	DESCRIPTION(
	  "Turns on emission by setting the emission flag to 1. \n"),
	OUTPUT( emission_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("emissionOff"),
  	DESCRIPTION(
	  "Turns off emission by setting the emission flag to 0."),
	OUTPUT( emission_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("losCalc"),
  	DESCRIPTION(
          "Calculates the line-of-sight (LOS).\n"
          "\n"
          "See AUG for details about the calculations."),
	OUTPUT( los_, z_tan_ ),
	INPUT( z_plat_ ,za_pencil_, l_step_, p_abs_, z_abs_, 
                refr_, refr_lfac_, refr_index_, z_ground_, r_geoid_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("sourceCalc"),
  	DESCRIPTION(
          "Calculates source function values valid between the points "
          "of the LOS.\n"
          "\n" 
          "No scattering and local thermodynamic equilibrium are assumed,\n"
          "that is, the source function equals the Planck function.\n"
          "The source function is set to the mean of the Planck function at\n"
          "the two LOS points limiting the steps. The temperature at the LOS\n"
          "points is obtained by linear interpolation.\n"
          "   If emission is neglected (emission=0), the WSV source is set \n"
          "to be empty."),
	OUTPUT( source_ ),
	INPUT( emission_, los_, p_abs_, t_abs_, f_mono_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("transCalc"),
  	DESCRIPTION(
          "Calculates the transmission between the points of the LOS.\n"
          "\n"
          "The absorption is assumed to vary linear between the LOS points."
          "The absorption at the LOS points is obtained by linear\n"
          "interpolation of *abs*."),
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
          "Standard choices for the radiation entering the atmosphere at\n"
          "the top of the atmosphere. \n"
          "\n"
          "The selections are:\n"
          "  zero : no radiation\n"
          "  cbgr : cosmic background radiation (planck for COSMIC_BG_TEMP)\n"
          "  sun  : solar radiation (planck for SUN_TEMP)\n"
          "\n"
          "COSMIC_BG_TEMP and SUN_TEMP are global variables, defined in\n"
          "constants.cc.\n"
          "\n"
          "Keywords \n"
          "  choice : Selection String (see above)."),
	OUTPUT( y_space_ ),
	INPUT( f_mono_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "nr"  ),
	TYPES(    String_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("yCalc"),
  	DESCRIPTION(
          "Performs the integration of the radiative transfer equation\n"
          "along the LOS, with or without emission.\n"
          "\n"
          "If emission is considered (emission=1) the outout unit is \n"
          "intensity, while without emission (emission=0) optical \n"
          "thicknesses are returned. "),
	OUTPUT( y_ ),
	INPUT( emission_, los_, f_mono_, y_space_, source_, trans_, 
                                                     e_ground_, t_ground_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("sourcetransyCalcSaveMemory"),
  	DESCRIPTION(
          "Combines:\nsourceCalc\ntransCalc\nyCalc\n\n"
          "Calculation is performed in frequency chunks thus allowing\n"
          "larger jobs to run. This means you cannot use it if you want \n"
           "to calculate Jacobians later on. "),
	OUTPUT( y_ ),
	INPUT( emission_, los_, p_abs_, t_abs_, f_mono_, abs_, y_space_,
                                          e_ground_, t_ground_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS("f_chunksize"),
	TYPES(Index_t)));

  md_data.push_back
    ( MdRecord
      ( NAME("CoolingRates"),
  	DESCRIPTION(
     "Calculates cooling rates due to exchange of longwave radiation.\n"
     "\n"
     "This function applies a straightforward algorith to obtain cooling \n"
     "rates. The algorithm is described in AUG. The basic idea is to \n"
     "calculate incoming radiation from all directions (by using yCalc), \n"
     "instead of follow the vertical flux through the atmosphere as usually \n"
     "done. No assumptions on a flat Earth is made.\n"
     "\n"
     "The atmosphere is described in usual way and absorption shall be pre- \n"
     "calculated. Emission is always activated (of course) and *y_space* is \n"
     "set to cosmic background radiation by calling *y_spaceStd*. Refraction\n"
     "follows corresponding WSV. \n" 
     "\n"
     "The zenith angle grid is set by *Za_pencil*, where the vector must \n"
     "start with 0 and end with 180. \n"
     "\n"
     "The WSV *l_step* is here treated to give the radiative step length for\n"
     "the zenith and nadir directions. The step length is scaled by \n"
     "abs( 1/cos(za) ), where za is the zenith angle, for other directions. \n"
     "The keyword argument *lstep_limit* sets an upper limit for *l_step*.  \n"
     "For example, for za=90, the expression above gives infinity for  \n"
     "*l_step*.\n"
     "\n"
     "The function returns the spectral cooling rate (by *coolrate*).\n"
     "\n"
     "Keywords \n"
     "  lstep_limit : Upper limit on l_step for off-zenith/nadir angles."),
	OUTPUT( coolrate_ ),
	INPUT( l_step_, p_abs_, z_abs_, t_abs_, f_mono_, abs_, za_pencil_,
               refr_, refr_lfac_, refr_index_, 
               r_geoid_, z_ground_, e_ground_, t_ground_, p_coolrate_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "lstep_limit" ),
	TYPES(    Numeric_t)));

  md_data.push_back
    ( MdRecord
      ( NAME("yTB"),
  	DESCRIPTION(
           "Converts a radiance spectrum to Planck brightness temperatures.\n"
           "\n"
           "The conversion is done by the Planck expression.\n"
	   "   The frequency of each value of *y* is determined by *f_mono* \n"
           "and za_pencil."),
	OUTPUT( y_ ),
	INPUT( y_, f_mono_, za_pencil_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixTB"),
	DESCRIPTION(
	   "Converts a radiance matrix to Planck brightness temperatures.\n"
           "\n"
           "Applies the function yTB on each column of the matrix. \n"
           "\n"
           "Generic input: \n"
           "   Matrix : Any matrix, but typically *ybatch* or a WF matrix.\n"
           "\n"
           "Generic output: \n"
           "   Matrix : Any matrix, but typically the same as the input \n"
           "            matrix."),
	OUTPUT(),
	INPUT(   f_mono_, za_pencil_ ),
	GOUTPUT( Matrix_ ),
	GINPUT(  Matrix_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("yTRJ"),
  	DESCRIPTION(
           "Converts a radiance spectrum to Rayleigh-Jean temperatures.\n"
           "\n"
           "The conversion is done by the Rayleigh-Jean approximation of the\n"
           "Planck expression.\n"
	   "   The frequency of each value of *y* is determined by *f_mono* \n"
           "and za_pencil."),
	OUTPUT( y_ ),
	INPUT( y_, f_mono_, za_pencil_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixTRJ"),
	DESCRIPTION(
	   "Converts a radiance matrix to Rayleigh-Jean temperatures.\n"
           "\n"
           "Applies the function yTRJ on each column of the matrix. \n"
           "\n"
           "Generic input: \n"
           "   Matrix : Any matrix, but typically *ybatch* or a WF matrix.\n"
           "\n"
           "Generic output: \n"
           "   Matrix : Any matrix, but typically the same as the input \n"
           "            matrix."),
	OUTPUT(),
	INPUT(   f_mono_, za_pencil_ ),
	GOUTPUT( Matrix_ ),
	GINPUT(  Matrix_ ),
	KEYWORDS(),
	TYPES()));



//======================================================================
//=== Weighting function (WF) methods
//======================================================================

  md_data.push_back
    ( MdRecord
      ( NAME("wfs_tgsDefine"),
  	DESCRIPTION(
          "Set up the list of tag groups for which weighting functions will \n"
	  "be calculated. \n"
          "\n"
          "The *wfs_tgs* are specified exactly as *tgs* (see tgsDefine). \n"
          "The selected tag groups must be a subgroup of the absorption \n"
          "tags (*tgs*). \n"
	  "   See the functions abs_per_tgReduce and kSpecies for some more \n"
          "information around *wfs_tgs*. \n"
          "\n"
          "Keywords \n"
          "  wfs_tgs : String with tag groups."),
	OUTPUT( wfs_tgs_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "wfs_tgs" ),
	TYPES(    Array_String_t   )));

md_data.push_back
    ( MdRecord
      ( NAME("wfss_tgsDefine"),
  	DESCRIPTION(
          "Set up the list of tag groups for which weighting functions will \n"
	  "be calculated. \n"
          "\n"
          "The *wfs_tgs* are specified exactly as *tgs* (see tgsDefine). \n"
          "The selected tag groups must be a subgroup of the absorption \n"
          "tags (*tgs*). \n"
	  "   See the functions abs_per_tgReduce and kSpecies for some more \n"
          "information around *wfs_tgs*. \n"
          "\n"
          "Keywords \n"
          "  wfss_tgs : String with tag groups."),
	OUTPUT( wfss_tgs_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "wfss_tgs" ),
	TYPES(    Array_String_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("absloswfsCalc"),
  	DESCRIPTION(
          "Calculates absorption line of sight weighting functions (LOS WFs)\n"
          "\n"
          "These WFs are the derivative of the spectra with respect to the \n"
          "absorption at the LOS points. See AUG for more detailed \n"
          "definition and details about the calculations."),
	OUTPUT( absloswfs_ ),
	INPUT( emission_, los_, source_, trans_, y_, y_space_, f_mono_, 
                                                        e_ground_, t_ground_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kSpecies"),
  	DESCRIPTION(
          "Calculates species weighting functions (WFs) for all *wfs_tgs*.\n"
          "\n"
          "This function is the simplest option if a single retrieval grid \n"
          "and a single retrieval unit are used for all species. If this is \n"
          "not the case, the function kSpeciesSingle must be used.\n"
          "   The WFs are calculated by (semi-)analytical expressions, where\n"
          "it is assumed that there is a linear relationship between the \n"
          "amount of the species and the absorption, and that the LOS is not\n"
          "affected by changes of the species. These assumtions should be\n"
          "valid generally for observations above the tropopause (as long \n"
          "LTE applies), but is not true for tropospheric water vapor. \n"
          "See AUG for details about the calculations. \n"
          "   The WFs for the different tag groups in *wfs_tgs* are appended\n"
          "to form a single matrix. The absorption array (*abs_per_tg*) must\n"
          "have been reduced to match *wfs_tags* (by using the function \n"
          "abs_per_tgReduce). The unit of the returned WFs are described \n"
          "below.\n"
          "\n"
          "The avaliable units are\n"
          "  frac : fractions of linearisation profile \n"
          "  vmr  : volume mixing ratio \n"
          "  nd   : number density\n"
          "\n"
          "Keywords \n"
          "  unit : Retrieval unit string (see above)."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( los_, absloswfs_, p_abs_, t_abs_, wfs_tgs_, abs_per_tg_, 
                                                              vmrs_, k_grid_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "unit"   ),
	TYPES(    String_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kSpeciesSingle"),
  	DESCRIPTION(
          "Calculates species weighting functions (WFs) for a single tag \n"
          "group.\n"
          "\n"
          "The tag group is selected by the giving the full name. This \n"
          "string must match exactly the string in *wfs_tgs* (and then the \n"
          "string in *tgs*). Otherwise as the function kSpecies (this \n"
          "including units). \n"
          "\n"
          "Keywords \n"
          "  tg   : Tag group string.\n"
          "  unit : Retrieval unit string (see kSpecies)."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( los_, absloswfs_, p_abs_, t_abs_, wfs_tgs_, abs_per_tg_, 
                                                              vmrs_, k_grid_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "tg",      "unit"  ),
	TYPES(    String_t,  String_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kContAbs"),
  	DESCRIPTION(
          "Calculates weighting functions (WFs) for polynomial fit of \n"
          "continuum  absorption. \n"
          "\n"  
          "The continuum is fitted by determining an off-set at a number of \n"
          "points (order+1) that are evenly spread between the lowest and \n"
          "upper frequency limit. See AUG for more details.\n"
          "   If the limits are set to be negative, *f_low* is set to the \n"
          "first value of *f_mono*, and *f_high* to the last value of \n"
          "*f_mono*. The frequency limits cannot be outside the range of \n"
          "*f_mono*. \n"
          "   The WFs can be calculated for different length units, selected\n"
          "by the keyword *l_unit*. For example, if *l_unit* is set to \n"
          "\"km\", the WFs corresponds to an absorpion with unit [1/km]. \n"
          "   The WFs for each frequency point are kept together, and the \n"
          "WF matrix for the different frequency points are appended. \n"
          "\n"
          "Keywords \n"
          "  order : Polynomial order (>=0). \n"
          "  f_low : Frequency of first fit point. \n"
          " f_high : Frequency of last fit point. \n"
          " l_unit : Length unit. Avaliable units are \"m\" and \"km\"." ),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( los_, absloswfs_, f_mono_, k_grid_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "order", "f_low",   "f_high",    "l_unit" ),
	TYPES(    Index_t,   Numeric_t, Numeric_t, String_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kTemp"),
  	DESCRIPTION(
          "Calculates temperature weighting functions (WFs).\n"
          "\n"
          "The calculations can be performed both with and without \n"
          "hydrostatic equilibrium (HSE). \n"
          "   If HSE is not considered (hse=0), the WFs are obtained by \n"
          "semi-analytical expressions and the calculations are relatively \n"
          "fast. See AUG for details. \n"
          "  If HSE is considered (hse=1), perturbation calculations are \n"
          "done. If the keyword *fast* is set to 0, the absorption is re-\n"
          "calculated for each temperature disturbance and the calculations\n"
          "are slow. With fast=1, it is assumed that the absorption is \n" 
          "linear with temparature between the present state and 1K higher \n"
          "temperature, and the new absorption is calculated once, that \n"
          "decreases the total calculation time considerbly. The accuracy of\n"
          "the latter option should suffice normally. \n"
          "   Note that the keyword *hse* here is not the workspace variable\n"
          "*hse*. If the keyword *hse* is set to 1, a constraint is that \n"
          "HSE is considered generally, that is, that the do-field of the \n"
          "variable 'hse* is turned on. The data to calculate assure HSE is\n"
          "of course taken from the workspace variable *hse*. \n"
          "   The fast keyword has no importance if the keyword hse is set \n"
          "to 0. \n"
          "\n"
          "Keywords \n"
          "  hse : Flag for hydrostatic eq. 0=no HSE, 1=HSE. \n"
          " fast : Flag to perform fast calculations with hse=1. " ),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_, abs_, 
          lines_per_tg_, lineshape_, e_ground_, emission_, k_grid_, 
	  cont_description_names_, cont_description_parameters_,
          cont_description_models_, los_, absloswfs_, trans_,
          z_plat_ ,za_pencil_, l_step_, z_abs_, refr_, refr_lfac_, refr_index_,
	  refr_model_, z_ground_, t_ground_, y_space_, r_geoid_, hse_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "hse",   "fast"  ),
	TYPES(    Index_t, Index_t )));
	
  md_data.push_back
    ( MdRecord
      ( NAME("kSpectro"),
  	DESCRIPTION(
	  "Calculates the spectroscopic parameters weighting functions (WFs).\n"
          "\n"
          "The calculation can be performed for the intensity,  line position, pressure \n"
           "broadening parameters  and pressure shift.\n"
           "For each parameter a do flag has to be specified.\n"
          "\n"
          "Keywords \n"
	  "do_intens: Flag for calculating the weighting function for"
           "the intensity do_intens=1. \n" 
	  "do_position: flag for line possition. \n"
	  "do_agam: flag for agam. \n"
	  "do_sgam: flag for sgam. \n" 
	  "do_nair: flag for temperature dependence of agam. \n"
	  "do_nself: flag for temperature dependence of sgam. \n" 
	  "do_pSift: flag for pressure shift .\n" ),

	OUTPUT(k_, k_names_, k_aux_ , S_S_),
	INPUT( wfss_tgs_, tgs_, f_mono_, p_abs_, t_abs_, z_abs_, h2o_abs_, vmrs_, 
          lines_per_tg_, lineshape_, los_, absloswfs_),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS("do_intens",  "do_position",  "do_agam", 
                  "do_sgam", "do_nair", "do_nself", "do_pSift"),
	TYPES(  Index_t, Index_t, Index_t, Index_t, Index_t, Index_t, Index_t)));
	
  md_data.push_back
    ( MdRecord
      ( NAME("kFrequencyOffSet"),
  	DESCRIPTION(
          "Calculates the weighting function (WF) for a frequency off-set.\n"
          "\n"
          "The Wf is simply the difference between *y* and the spectrum \n"
          "obtained when adding *delta* to *f_mono*, diveded by *delta*.\n"
          "That is, a pure perturbation calculation is performed. \n"
          "\n"
          "   The WF can be calculated for different frequency units,\n"
	  "selected by the keyword *f_unit*.\n"
          "\n"
          "Keywords \n"
          "  delta  : Size of frequency perturbation (in units of *l_unit*).\n"
          "  l_unit : Frequency unit. Avaliable units are \"Hz\", \"kHz\"\n"
          "           \"MHz\"." ),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_, 
          lines_per_tg_, lineshape_, e_ground_, emission_, 
	  cont_description_names_, cont_description_parameters_,
          cont_description_models_, los_, t_ground_, y_space_, y_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "delta",   "f_unit"  ),
	TYPES(    Numeric_t, String_t  )));

  md_data.push_back
    ( MdRecord
      ( NAME("kPointingOffSet"),
  	DESCRIPTION(
          "Calculates the WF for a pointing off-set.\n"
          "\n"
          "The Wf is simply the difference between *y* and the spectrum \n"
          "obtained when adding *delta* to *za_pencil*, diveded by *delta*.\n"
          "That is, a pure perturbation calculation is performed. \n"
          "\n"
          "Keywords \n"
          "  delta : Size of zenith angle perturbation."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( z_plat_, za_pencil_, l_step_, p_abs_, z_abs_, t_abs_, f_mono_,
               refr_, refr_lfac_, refr_index_, z_ground_, r_geoid_, 
               abs_, emission_, y_space_, e_ground_, t_ground_, y_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "delta"   ),
	TYPES(    Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME( "kEground" ),
	DESCRIPTION(
	   "Calculates the WF(s) for ground emission coefficent(s).\n"
           "\n"
           "The ground emission WF(s) are calculated by semi-analytical\n"
           "expressions (see AUG). With single_e=0, a WF is returned for \n"
           "the emission coefficient of each monochromatic frequency. \n"
           "On the other hand, when single_e=1, the ground emission is \n"
           "treated as a single varaible (that is, no frequency dependency) \n"
           "and there is only a single WF to be calculated. The latter \n"
           "option requieres that all elements of E_GROUND are set to the \n"
           "same value. \n"
           "\n"
           "Keywords:\n"
           "   single_e : Boolean to treat the ground emission as a single\n"
           "              variable. See further above." ),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( za_pencil_, f_mono_, emission_, y_space_, e_ground_, t_ground_,
               los_, source_, trans_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "single_e" ),
	TYPES(    Index_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kCalibration"),
  	DESCRIPTION(
          "Calculates the WF for a proportional calibration error. \n"
          "\n"
          "The WF is simply : k = y - y0 where y0 is the specified \n"
          "vector. The y0-vector shhould typically be the radiance (or TB) \n"
          "of the load used for load switching. For example: \n"
          "   VectorPlanck(y0,f_mono){temp=2.7} \n"
          "   kCalibration(y0){} \n"
          "\n"
          "Generic input: \n"
	  "   Vector : A vector with spectrum for calibration reference \n"
          "            point. This vector should typically be *y0*. "),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( y_, f_mono_ ),
	GOUTPUT(),
	GINPUT( Vector_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kManual"),
  	DESCRIPTION(
          "Calculates a weighting function using y and y0.\n"
          "\n"
          "The weighting function is calculated as: k = (y-y0)/delta\n"
          "That is, delta is the magnitude of the perturbation done.\n"
          "\n"
          "Keywords \n"
          "  name    : Name on retrieval/error identity.\n"
          "  delta   : Magnitude of perturbation.\n"
          "  grid    : Grid point value.\n"
          "  apriori : A priori value."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( y0_, y_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "name",   "delta",   "grid",    "apriori" ),
	TYPES(    String_t, Numeric_t, Numeric_t, Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kxInit"),
  	DESCRIPTION(
          "Initializes Kx weighting function matrix and help variables\n"
          "(kx_names, kx_lengths and kx_aux).\n"
          "\n"
          "Use this function before the WF calculations are started and\n"
          "together with kxAppend."),
	OUTPUT( kx_, kx_names_, kx_lengths_, kx_aux_ ),
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
          "(kb_names, kb_lengths and kb_aux).\n"
          "\n"
          "Use this function before the WF calculations are started and\n"
          "together with kbAppend."),
	OUTPUT( kb_, kb_names_, kb_lengths_, kb_aux_ ),
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
          "correspondingly. \n"
          "\n"
          "All the data are reallocated to make space for the new data.\n"
          "This function is accordingly slow for large data sizes,\n"
          "and it can be better to use kxAllocate and kxPutInK."),
	OUTPUT( kx_, kx_names_, kx_lengths_, kx_aux_ ),
        INPUT( kx_, kx_names_, kx_lengths_, kx_aux_, k_, k_names_, k_aux_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kbAppend"),
  	DESCRIPTION(
          "Appends the K matrix to Kb and handles additional data\n"
          "correspondingly. \n"
          "\n"
          "All the data are reallocated to make space for the new data.\n"
          "This function is accordingly slow for large data sizes,\n"
          "and it can be better to use kbAllocate and kbPutInK."),
	OUTPUT( kb_, kb_names_, kb_lengths_, kb_aux_ ),
        INPUT( kb_, kb_names_, kb_lengths_, kb_aux_, k_, k_names_, k_aux_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kxAllocate"),
  	DESCRIPTION(
          "Allocates memory for kx and help variables (kx_names, kx_lengths \n"
          "and kx_aux).\n" 
          "\n"
          "The total number of frequencies is taken from the length \n"
          "of the given vector (typically y).\n"
          "   Use this function before the WF calculations are started and\n"
          "together with kxPutInK.\n"
          "\n"
          "Generic input: \n"
	  "   Vector : A vector with same length as the appended spectra.\n"
          "            The typical choice is *y*.\n"
          "\n"
          "Keywords \n"
          "  ni : Number of retrieval identities (species profiles,\n"
          "              pointing off-set etc.).\n"
          "  nx : Final length of x."),
	OUTPUT( kx_, kx_names_, kx_lengths_, kx_aux_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT( Vector_ ),
	KEYWORDS( "ni",  "nx"  ),
	TYPES(    Index_t, Index_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kbAllocate"),
  	DESCRIPTION(
          "Allocates memory for kx and help variables (kb_names, kb_lengths \n"
          "and kb_aux). \n"
          "\n"
          "The total number of frequencies is taken from the length \n"
          "of the given vector (typically y).\n"
          "   Use this function before the WF calculations are started and\n"
          "together with kbPutInK.\n"
          "\n"
          "Generic input: \n"
	  "   Vector : A vector with same length as the appended spectra.\n"
          "            The typical choice is *y*.\n"
          "\n"
          "Keywords \n"
          "  ni : Number of retrieval identities (species profiles,\n"
          "              pointing off-set etc.).\n"
          "  nb : Final length of b."),
	OUTPUT( kb_, kb_names_, kb_lengths_, kb_aux_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT( Vector_ ),
	KEYWORDS( "ni",  "nb"  ),
	TYPES(    Index_t, Index_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kxPutInK"),
  	DESCRIPTION(
          "Puts K in Kx and handles additional data correspondingly.\n"
          "\n"
          "K is placed in the first free columns of Kx.\n"
          "   No reallocation is performed (in contrast to kxAppend) and an\n"
          "error message is given if k does not fit into kx. The kx-data are\n"
          "allocated by the function kxAllocate."),
	OUTPUT( kx_, kx_names_, kx_lengths_, kx_aux_ ),
        INPUT( kx_, kx_names_, kx_lengths_, kx_aux_, k_, k_names_, k_aux_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kbPutInK"),
  	DESCRIPTION(
          "Puts K in Kb and handles additional data correspondingly.\n"
          "\n"
          "K is placed in the first free columns of Kb.\n"
          "   No reallocation is performed (in contrast to kbAppend) and an\n"
          "error message is given if k does not fit into kb. The kb-data are\n"
          "allocated by the function kbAllocate."),
	OUTPUT( kb_, kb_names_, kb_lengths_, kb_aux_ ),
        INPUT( kb_, kb_names_, kb_lengths_, kb_aux_, k_, k_names_, k_aux_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));



//======================================================================
//=== Batch Calculation Methods
//======================================================================

#ifdef HDF_SUPPORT
  md_data.push_back
    ( MdRecord
      ( NAME("ybatchCalc"),
	DESCRIPTION(
          "Calculates a batch of spectra from a set of profiles and\n"
	  "frequency and viewing angle grids.\n"
          "The following workspace methods are used:\n"
          "   absCalc    \n"
          "   losCalc    \n"
          "   sourceCalc \n"
          "   transCalc  \n"
          "   yCalc  \n"
          "and all the workspace variables needed by these functions must be\n"
          "set before starting this method.\n"
          "The refractive index is kept constant for all spectra.\n"
          "The values of the workspace variables are used as defaults when\n"
          "appropiate. For example, if temperature profiles are not read\n"
          "from a file (do_t=0), then the values of t_abs are used for all\n"
          "spectra.\n"
          "All input files shall be readable by MatrixReadBinary.\n"
          "The profiles and grids are stored as columns in the file matrix.\n"
          "When a filename is empty (""), filenames are created as:\n"
          "   batchname.XXX.ab \n"
          "where XXX is\n"
          "   t_abs     : For temperature profiles.\n" 
          "   z_abs     : For vertical altitude grids.\n" 
          "   f_mono    : For frequency grids.\n" 
          "   za_pencil : For zenith angle grids.\n" 
          "For species XXX is the molecule name, e.g. H2O and O3.\n" 
          "If a filename is given, batchname is ignored.\n"
          "When a flag is 0, the corresponding filename is of no importance.\n"
          "The length of profiles (t, z and species) must match p_abs (no\n"
          "interpolation is performed).\n" 
          "\n"
          "Keywords \n"
          "  ncalc     : The number of spectra to calculate. The files can\n"
          "              contain data for more spectra (but not less).\n"
          "  do_t      : Temperature flag (0/1).\n"
          "  t_file    : Filename for temperature data.\n"
          "  do_z      : Vertical altitude flag (0/1).\n"
          "  z_file    : Filename for vertical grid data.\n"
          "  do_f      : Frequency flag (0/1).\n"
          "  f_file    : Filename for frequency data.\n"
          "  do_za     : Zenith angle flag (0/1).\n"
          "  za_file   : Filename for zenith angle data.\n"
          "  do_tags   : This String array gives the tags for which profiles\n"
          "              shall be read from a file, e.g. [\"H2O\",\"O3\"].\n"
          "              These tags must match some tag in tags.\n"
          "  tag_files : Filenames for species data."),
	OUTPUT( ybatch_ ),
	INPUT( // Variables needed for absCalc
               f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_, 
	       lines_per_tg_, lineshape_, 
               // Additional variables for losCalc
	       z_abs_, z_plat_ ,za_pencil_, l_step_, refr_, 
               refr_lfac_, refr_index_, z_ground_, r_geoid_,
               // Additional variables for yRte
	       emission_, y_space_, e_ground_, t_ground_,
               // Additional variables needed for this function
               batchname_, tgs_, cont_description_names_, cont_description_parameters_,
	       cont_description_models_),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS("ncalc", "do_t", "t_file", "do_z", "z_file",
                 "do_tags", "tag_files",
                 "do_f", "f_file", "do_za", "za_file"),
	TYPES(   Index_t,   Index_t,  String_t, Index_t,  String_t, 
                 Array_String_t, Array_String_t,
                 Index_t,  String_t, Index_t,   String_t  )));

#endif // HDF_SUPPORT
  
  md_data.push_back
    ( MdRecord
      ( NAME("ybatchFromRadiosonde"),
	DESCRIPTION
	("Calculate spectra for a batch of radiosonde data."
	 "\n"
	 "   We set the oxygen and nitrogen VMR to constant values of 0.209 and\n"
	 "0.782, respectively. Some other methods are called implicitly by this\n"
	 "method. Specifically:\n"
	 "- absCalc\n"
	 "- refrCalc\n"
	 "- losCalc\n"
	 "- sourceCalc\n"
	 "- transCalc\n"
	 "- yCalc\n"
	 "\n"
	 "Keywords: \n"
	 "\n"
         "  finegrid : Flag for a fine *p_abs* grid, 0 = Radiosonde levels, 1 = finer grid. \n"
         "\n"
         "If the keyword finegrid is set to 0 (finegrid = 0),\n"
         "absorption coeff. are calculated on the same grid as in the radiosonde launch.\n" 
         "It does not check whether the launch has reached up to a certain height. \n"
         "\n"
         "If finegrid = 1, the absorption is calculated on a very fine grid (about 15 m)\n"
         "which is finer than the high resolution radiosonde levels (about 60 m). In this\n"
         "case *p_abs* grid is only up to 100 hPa so that we assume there is no atmosphere\n"
         "above this pressure level. The RT calculation is not done for any profile which\n"
         "do not fly up to 100 hPa. In this case we put -1 as the Tbs for all the frequencies.\n"
         "Planck brightness temperature is calculated for all the launches those reach 100 hPa.\n"
         "Note that *l_step* which is given by the user in the control file should be less than\n"
	 "15 m, may be 5 m.\n"
         "\n"
         "  interp_rh : Flag for interpolation of H2O profile in RH\n"
         "  0 = Normal ARTS interpolation in VMRs\n"
         "  1 = Interpolation in RH."),
	OUTPUT( ybatch_ ),
	INPUT( // Variables needed for absCalc
               radiosonde_data_, f_mono_, lines_per_tg_, lineshape_, 
               // Additional variables for losCalc
	       z_plat_ ,za_pencil_, l_step_, refr_, refr_model_,
               refr_lfac_, r_geoid_,
               // Additional variables for yRte
	       emission_, y_space_, e_ground_, 
               // Additional variables needed for this function
               tgs_, 
               cont_description_names_, 
               cont_description_models_, 
               cont_description_parameters_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "finegrid" ,  "interp_rh"),
	TYPES(  Index_t, Index_t  )));

md_data.push_back
    ( MdRecord
      ( NAME("ybatchFromRadiosondeGlobal"),
	DESCRIPTION
	("Calculate spectra for a batch of Global radiosonde data."
	 "\n"
	 "This method is almost similar to #ybatchFromRadiosonde#. Since the pressure \n"
	 "grid in the global radiosonde data is coarse it is to be interpolated in a \n"
	 "fine grid.\n"
	 "\n"
	 "   We set the oxygen and nitrogen VMR to constant values of 0.209 and\n"
	 "0.782, respectively. Some other methods are called implicitly by this\n"
	 "method. Specifically:\n"
	 "- absCalc\n"
	 "- refrCalc\n"
	 "- losCalc\n"
	 "- sourceCalc\n"
	 "- transCalc\n"
	 "- yCalc"),
	OUTPUT( ybatch_ ),
	INPUT( // Variables needed for absCalc
               radiosonde_data_, f_mono_, lines_per_tg_, lineshape_, 
               // Additional variables for losCalc
	       z_plat_ ,za_pencil_, l_step_, refr_, refr_model_,
               refr_lfac_, r_geoid_,
               // Additional variables for yRte
	       emission_, y_space_, e_ground_, 
               // Additional variables needed for this function
               tgs_, 
               cont_description_names_, 
               cont_description_models_, 
               cont_description_parameters_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES(     )));

//======================================================================
//=== Methods as Workspace Variables
//======================================================================

  md_data.push_back
    ( MdRecord
      ( NAME("MethodListDefine"),
  	DESCRIPTION
	( "Set up a method list.\n"
	  "\n"
	  "A method list just contains indices (in md_data) of methods\n"
	  "intended for sequential execution. Only methods without keyword\n"
	  "arguments are allowed. It is the task of this method to\n"
	  "set this up. For example, it must be checked, whether the given\n"
	  "names really correspond to methods.\n"
	  "\n"
	  "Generic Output:\n"
	  "   ArrayOfIndex : The newly generated method list.\n"
	  "\n"
	  "Keywords:\n"
	  "   methods      : An array of names of methods." ),
	OUTPUT(  ),
        INPUT(  ),
	GOUTPUT( ArrayOfIndex_ ),
	GINPUT(),
	KEYWORDS( "methods" ),
	TYPES(    Array_String_t )));

}

