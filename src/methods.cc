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
	DESCRIPTION(
	   "A summary of the function in one sentence.\n"
           "\n"
           "A detailed description of the function. Please, try to be as \n"
           "clear and detailed as possible, this will help both you and \n"
           "others in the long run. \n"
           "   Additional paragraphs are indented with three blanks, as \n"
           "exemplified here.\n"
           "   The names of workspace variables and other methods\n"
           "are marked by stars, for example *z_plat*.\n"
           "   Global input and output, and keywords shall be described \n"
           "as exemplified below. If there is no variables of a group, \n"
           "(e.g. global input) remove that part totally. Note that the \n"
           "on-line help just gives the type of global input/output and the \n"
           "keyword names, and additional information is for sure needed.\n"
           "   Leave space and brake lines whem listing input and output \n"
           "variabales to make the code easier to read. See example below. \n"
           "\n"
           "Global input: \n"
           "   Vector : Vector giving some very important input. Don't \n"
           "            be too short. Use the type of indention used here. \n"
           "\n"
           "Global output: \n"
           "   Vector : Return vector for the zenith angles. The normal \n"
           "            options are ZA_PENCIL and ZA_SENSOR. \n"
           "\n"
           "Keywords:\n"
           "   delta_t   : Time increment between observations.\n"
           "   z_tan_lim : Vector with start and stop tangent altitudes." ),
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
	DESCRIPTION(
	   "\n"
           "\n"
           "Global input: \n"
           "   \n"
           "\n"
           "Global output: \n"
           "   \n"
           "\n"
           "Keywords:\n"
           "   " ),
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
      ( NAME("Exit"),
	DESCRIPTION
	(
	 "Stops the execution and exits ARTS.\n"
	 "\n"
	 "This method is handy if you want to debug one of your\n"
	 "controlfiles. You can insert it anywhere in the controlfile. When it\n"
	 "is reached, it will terminate the program."
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
	 "This method can be used by ARTS developers to quickly test stuff. The\n"
	 "implementation is in file m_io.cc. This just saves you the trouble of\n"
	 "adding a dummy method everytime you want to try something out quickly."
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

 // These functions should be changed to handle Index

  md_data.push_back     
    ( MdRecord
      ( NAME("IndexSet"),
	DESCRIPTION("Sets an integer workspace variable to the given value."),
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
                    "See `ArrayOfMatrixWriteAscii` for file format.\n"
                    "\n"
                    "Global input: \n"
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
                    "See `ArrayOfMatrixWriteAscii` for file format.\n"
                    "\n"
                    "Global output: \n"
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
		    "Global input: \n"
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
		    "Global output: \n"
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
	DESCRIPTION("Sets a workspace variable of type Numeric to a value."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Numeric_ ),
	GINPUT(),
	KEYWORDS( "value"   ),
	TYPES(    Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("NumericCopyFirstOfVector"),
	DESCRIPTION(
           "Sets a workspace variable of type Numeric to the value of the"
           "first element in a vector." ),
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
                    "See `ArrayOfMatrixWriteAscii` for file format.\n"
                    "\n"
                    "Global input: \n"
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
                    "See `ArrayOfMatrixWriteAscii` for file format.\n"
                    "\n"
                    "Global output: \n"
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
		    "Global input: \n"
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
		    "Global output: \n"
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
      ( NAME("VectorSet"),
	DESCRIPTION("Creates a workspace vector with the specified length\n"
                    "and initializes the vector with the given value."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Vector_ ),
	GINPUT(),
	KEYWORDS( "length", "value"   ),
	TYPES(    Index_t,    Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorSetLengthFromVector"),
	DESCRIPTION("Creates a workspace vector with the same length as the\n"
		    "given vector and initializes the new vector with the\n"
                    "given value. For example\n"
                    " VectorSetLengthFromVector(e_ground,f_mono){value=0.75}"),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Vector_ ),
	GINPUT( Vector_ ),
	KEYWORDS( "value"   ),
	TYPES(    Numeric_t )));

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
	GOUTPUT( Vector_ ),
	GINPUT(),
	KEYWORDS( "start",   "stop",    "step"    ),
	TYPES(    Numeric_t, Numeric_t, Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorNLinSpace"),
	DESCRIPTION("Creates a vector with defined length, equally spaced\n"
                    "between the given values.\n"
		    "The length must be larger than 1."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(Vector_),
	GINPUT(),
	KEYWORDS( "start",   "stop",    "n"   ),
	TYPES(    Numeric_t, Numeric_t, Index_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorNLogSpace"),
	DESCRIPTION("Creates a vector with defined length, logarithmically\n"
                    "spaced between the given values.\n"
		    "The length must be larger than 1."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(Vector_),
	GINPUT(),
	KEYWORDS( "start",   "stop",    "n"   ),
	TYPES(    Numeric_t, Numeric_t, Index_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorCopy"),
	DESCRIPTION("Copies a vector."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Vector_ ),
	GINPUT( Vector_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorFlip"),
	DESCRIPTION(
           "Flips a vector. The result is the vector in reversed order"),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Vector_ ),
	GINPUT( Vector_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorCopyFromArrayOfVector"),
	DESCRIPTION("Copies a vector from a vector array.\n"
          "\n"
          "Keywords \n"
          "  index : The index of the vector in the array to copy. " ),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Vector_ ),
	GINPUT( ArrayOfVector_),
	KEYWORDS( "index" ),
	TYPES( Index_t )));

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
                    "See `ArrayOfMatrixWriteAscii` for file format.\n"
                    "\n"
                    "Global input: \n"
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
                    "See `ArrayOfMatrixWriteAscii` for file format.\n"
                    "\n"
                    "Global output: \n"
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
                    "Global input: \n"
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
                    "Global output: \n"
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
  	DESCRIPTION(
          "Sets a vector to the Planck function for the given frequency\n"
          "vector and temperature. An example:\n"
          "   VectorPlanck(y_space,f_mono){temp=2.7}"),
	OUTPUT( ),
	INPUT( ),
	GOUTPUT( Vector_ ),
	GINPUT( Vector_ ),
	KEYWORDS( "temp"    ),
	TYPES(    Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorCalcLog10"),
  	DESCRIPTION(
          "Calculates the base 10 logarithm of a vector.\n"
          "The result can either be stored in the same or another vector."),
	OUTPUT( ),
	INPUT( ),
	GOUTPUT( Vector_ ),
	GINPUT( Vector_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorAdd"),
  	DESCRIPTION(
	  "Adds a scalar to all elements of a vector."),
	OUTPUT( ),
	INPUT( ),
	GOUTPUT( Vector_ ),
	GINPUT( Vector_ ),
	KEYWORDS( "value" ),
	TYPES( Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorScale"),
  	DESCRIPTION(
	  "Scales all elements of a vector with the same value."),
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
	DESCRIPTION("Creates a workspace matrix of the specified size\n"
                    "and initializes the matrix with the given value."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Matrix_ ),
	GINPUT(),
	KEYWORDS( "nrows", "ncols", "value"   ),
	TYPES(    Index_t,   Index_t,   Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixCopy"),
	DESCRIPTION("Copies a matrix."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Matrix_ ),
	GINPUT( Matrix_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixFillWithVector"),
	DESCRIPTION("Forms a vector with n columns, and put the given.\n"
                    "in each column."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Matrix_ ),
	GINPUT( Vector_ ),
	KEYWORDS( "n"   ),
	TYPES(    Index_t )));

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
                    "See `ArrayOfMatrixWriteAscii` for file format.\n"
                    "\n"
                    "Global input: \n"
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
                    "See `ArrayOfMatrixWriteAscii` for file format.\n"
                    "\n"
                    "Global output: \n"
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
                    "Global input: \n"
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
                    "Global output: \n"
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
	DESCRIPTION(
          "Scales all elements of a matrix with the same value."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Matrix_ ),
	GINPUT( Matrix_ ),
	KEYWORDS( "value" ),
	TYPES(    Numeric_t   )));



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
                    "See `ArrayOfMatrixWriteAscii` for file format.\n"
                    "\n"
                    "Global input: \n"
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
                    "See `ArrayOfMatrixWriteAscii` for file format.\n"
                    "\n"
                    "Global output: \n"
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
                    "Global input: \n"
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
                    "Global output: \n"
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
                    "See `ArrayOfMatrixWriteAscii` for file format.\n"
                    "\n"
                    "Global input: \n"
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
                    "See `ArrayOfMatrixWriteAscii` for file format.\n"
                    "\n"
                    "Global output: \n"
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
                    "Global input: \n"
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
                    "Global output:  \n"
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
                    "Global input: \n"
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
                    "See `ArrayOfMatrixWriteAscii` for file format.\n"
                    "\n"
                    "Global output: \n"
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
                    "Global input: \n"
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
                    "Global output: \n"
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
                    "See `ArrayOfStringWriteAscii` for file format.\n"
                    "\n"
                    "Global input: \n"
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
                    "See `ArrayOfStringWriteAscii` for file format.\n"
                    "\n"
                    "Global output: \n"
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
                    "Global input: \n"
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
                    "Global output: \n"
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
                    "Global input: \n"
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
                    "See `ArrayOfStringWriteAscii` for file format.\n"
                    "\n"
                    "Global output: \n"
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
                    "Global input: \n"
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
                    "Global output: \n"
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

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixDiagonal"),
	DESCRIPTION("Creates a diagonal matrix.\n"
                    "All diagonal elements are set to the same value."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Matrix_ ),
	GINPUT(),
	KEYWORDS( "nrows", "value"   ),
	TYPES(    Index_t,   Numeric_t )));



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
                    "Global input: \n"
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
                    "Global output: \n"
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
      ( NAME("lines_per_tgReadFromCatalogues"),
  	DESCRIPTION(
	"This method can read lines from different line catalogues.\n"
	"For each tag group, you can specify which catalogue to use. Because\n"
	"the method creates lines_per_tg directly, it replaces for example the\n"
	"following two method calls:\n"
	"  - linesReadFromHitran\n"
	"  - lines_per_tgCreateFromLines\n"
	"\n"
	"This method needs as input WSVs the list of tag groups. Keyword\n"
	"parameters must specify the names of the catalogue files to use and\n"
	"the matching formats. Names can be anything, formats can currently\n"
	"be HITRAN96, MYTRAN2, JPL, or ARTS. Furthermore, keyword parameters\n"
	"have to specify minimum and maximum frequency for each tag group. To\n"
	"safe typing, if there are less elements in the keyword parameters\n"
	"than there are tag groups, the last parameters are applied to all\n"
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
	"\n"
	"In this example, lines for the first tag group will be taken from\n"
	"cat1, lines for all other tag groups will be taken from cat2.\n"
	"\n"
	"This methods allows you for example to use a special line file just\n"
	"for water vapor lines. This could be the improved water vapor line\n"
	"file generated by Thomas Kuhn.\n"
	"\n"
	"Catalogues are only read once, even if several tag groups have the\n"
	"same catalogue. However, in that case the frequency ranges MUST be\n"
	"the same. (If you want to do fine-tuning of the frequency ranges,\n"
	"you can do this inside the tag definitions, e.g., \"H2O-*-0-2000e9\".)\n"
	"\n"
	"This function uses the various reading routines\n"
	"(linesReadFromHitran, etc.), as well as\n"
	"lines_per_tgCreateFromLines."),
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
	TYPES(    String_t,    Numeric_t, Numeric_t)));

  md_data.push_back
    ( MdRecord
      ( NAME("linesReadFromMytran2"),
  	DESCRIPTION(
          "Read all the lines from a MYTRAN2 catalogue file that\n"
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
	TYPES(    String_t,    Numeric_t, Numeric_t)));

  md_data.push_back
    ( MdRecord
      ( NAME("linesReadFromJpl"),
  	DESCRIPTION(
          "Read all the lines from a JPL catalogue file that\n"
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
	TYPES(    String_t,    Numeric_t, Numeric_t)));

  md_data.push_back
    ( MdRecord
      ( NAME("linesReadFromArts"),
  	DESCRIPTION(
          "Read all the lines from an Arts catalogue file that\n"
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
	TYPES(    String_t,    Numeric_t, Numeric_t)));

  // FIXME: Remove this one.
  md_data.push_back
    ( MdRecord
      ( NAME("linesElowToJoule"),
  	DESCRIPTION(
          "Just a little helper to convert the lower state energy from\n"
	  " cm^-1 (ARTSCAT-2) to Joule (ARTSCAT-3). This should be removed soon."),
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
	  "The tag groups are tested in the order in which they are\n" 
          "specified in the controlfile. The line is assigned to the\n" 
	  "first tag group that fits."),
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
          "Adds mirror lines at negative frequencies to the lines_per_tg.\n"
	  "For each line at frequency +f in lines_per_tg a corresponding\n"
	  "entry at frequency -f is added to lines_per_tg.\n"
	  "The mirror lines are appended to the line lists after the original lines."),
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
	  "from the lines_per_tg, in order to save computation time.\n"
	  "It should be particularly useful to call this method after\n"
	  "lines_per_tgAddMirrorLines."),
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
                    "Writes the workspace variable `lines` to an ASCII file.\n"
                    "\n"
                    "The content of the workspace variable `lines`\n"
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
                    "Writes the workspace variable `lines_per_tg` to an ASCII file.\n"
                    "\n"
                    "The content of the workspace variable `lines_per_tg`\n"
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
	  "Specify one String for each tag group that you want to create.\n"
	  "Inside the String, separate the tags by comma (plus optional blanks).\n"
	  "Example:\n"
	  "tag = [\"O3-666-500e9-501e9, O3-686\",\"H2O\",\"O2-*-*-*\"]"),
	OUTPUT( tgs_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "tags" ),
	TYPES(    Array_String_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("lineshapeDefine"),
  	DESCRIPTION(
          "Sets the lineshape for all calculated lines. Specify an available\n"
	  "lineshape, together with a normalization factor and a cutoff frequency.\n\n"
	  "Shape: no_shape Doppler Lorentz Drayson Voigt_Kuntz Voigt_Drayson3\n"
	  "       Voigt_Drayson4 Voigt_Drayson6 Rosenkranz_Voigt_Drayson\n"
	  "       Rosenkranz_Voigt_Kuntz6.\n"
	  "Normalization Factors:  no_norm: 1 linear: f/f0  quadratic: (f/f0)^2.\n"
	  "Cutoff: -1: no cutoff Number: positive cutoff frequency in Hz.\n"
	  "Example:\n\n"
	  "shape=\"Lorentz\" normalizationfactor=\"linear\" cutoff=650e9"),
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
          "Sets the lineshape per tag group for all calculated lines. Specify\n"
	  "an available lineshape, together with a normalization factor and a\n"
	  "cutoff frequency.\n\n"
	  "Shape: no_shape Doppler Lorentz Drayson Voigt_Kuntz Voigt_Drayson3\n"
	  "       Voigt_Drayson4 Voigt_Drayson6 Rosenkranz_Voigt_Drayson\n"
	  "       Rosenkranz_Voigt_Kuntz6.\n"
	  "Normalization Factors:  no_norm: 1 linear: f/f0  quadratic: (f/f0)^2.\n"
	  "Cutoff: -1: no cutoff Number: positive cutoff frequency in Hz.\n\n"
	  "Example:\n"
	  "shape=[\"Lorentz\",\"Voigt_Kuntz6\"] \n"
	  "normalizationfactor=[\"linear\", \"quadratic\"] \n"
	  "cutoff=[ 650e9, -1 ]"),
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
	OUTPUT( cont_description_names_, cont_description_parameters_ ),
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
	 "   parameters : A Vector containing the required number of parameters\n"
	 "                for the model given. The meaning of the parameters and\n"
	 "                how many parameters are required depends on the model."
	 ),
	OUTPUT( cont_description_names_, cont_description_parameters_ ),
	INPUT(  cont_description_names_, cont_description_parameters_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "name",   "parameters" ),
	TYPES(    String_t, Vector_t     )));


//=== Input Atmosphere methods ===========================================

  md_data.push_back
    ( MdRecord
      ( NAME("raw_vmrsReadFromFiles"),
        DESCRIPTION(
          "Read the individual VMR profile for each TAGS from the list of\n"
          "files given as keyword parameters. One file name must be specified for\n"
          "each TAGS. The name may include a path."
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
	  "Read the individual VMR profile for each tag group from a standard\n"
	  "atmospheric scenario. Files must look like this example:\n"
	  "<basename>.ClO.aa\n"
	  "\n"
	  "The basename can include a path, i.e., the files can be anywhere,\n"
	  "but they must be all in the same directory.\n"
	  "\n"
	  "The profile is chosen by the species name. If you have more than one\n"
	  "tag group for the same species, the same profile will be used."
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
	  "given by p_abs. The altitude is not used by the absorption routines,\n"
	  "But later on by the RT routines."
	  "\n"
	  "Interpolations used: FIXME: Add these.\n"
	  "Temperature [K]: \n"
	  "Altitude    [m]: \n"
	  "VMRs        [1]: \n"
	  "\n"
	  "Uses interp_lin(...)."
	  ),
	OUTPUT(   t_abs_    , z_abs_   , vmrs_           ),
	INPUT(    tgs_, p_abs_    , raw_ptz_ , raw_vmrs_ ),
	GOUTPUT(                       			 ),         
	GINPUT(                        			 ),
	KEYWORDS(  "CloudSatWV"             		 ),
	TYPES(     String_t                    		 )));

  md_data.push_back
    ( MdRecord
      ( NAME("h2o_absSet"),
	DESCRIPTION(
          "Sets h2o_abs to the profile of the first tag group containing\n"
	  "water."),
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
	  "molecular nitrogen."),
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
          "equilibrium (HSE). The on/off flag is set to 1. \n"
          "Type 'arts -d hse' for more information. \n"
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
          "The on/off flag off hse is set to 0 and hse is set to be a vector\n"
          "of length 1."),
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
          "Ensures that 'z_abs' fulfills hydrostatic equilibrium.\n"
          "Nothing is done if the on/off flag of 'hse' is set to 0.\n"
          "Calculates a vertical grid fulfilling hydrostatic equilibrium\n"
          "The given altitudes (i.e. z_abs) are used as a first guess when \n"
          "starting the calculations (to estimate g etc.).\n"
          "The altitude variation of the gravitational acceleration is\n"
          "considered.\n"
          "The average molecular weight is assumed to be 28.96 at all\n"
          "altitudes.\n"
          "The amount of water vapour is taken into account.\n"
          "The reference point, g at the ground and the number of iterations\n"
          "are taken from 'hse'. A higher number of iterations \n" 
          "improves the accuracy, but one iteration should be normally \n"
          "enough if z_abs already has reasonable values. Two iterations \n"
          "should suffice for basically all applications."),
	OUTPUT( z_abs_ ),
	INPUT( z_abs_, p_abs_, t_abs_, h2o_abs_, r_geoid_, hse_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));


  md_data.push_back
    ( MdRecord
      ( NAME("vmrsScale"),
	DESCRIPTION(
          "Scales the vmr input of the tgs given in scaltgs by the\n"
	  "factors given in scalfac."),
	OUTPUT(	vmrs_ ),
	INPUT( 	tgs_, vmrs_  ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "scaltgs", "scalfac"),
	TYPES( Array_String_t, Vector_t)));


//=== 1D absorption methods ===============================================

  md_data.push_back
    ( MdRecord
      ( NAME("absCalc"),
	DESCRIPTION("Calculate absorption coefficients. This\n"
		    "calculates both the total absorption and the\n"
		    "absorption per tag group."),
	OUTPUT(	    abs_  , abs_per_tg_                         ),
	INPUT( 	    tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_, 
                    lines_per_tg_, lineshape_,
		    cont_description_names_, cont_description_parameters_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("absCalcFromXsec"),
	DESCRIPTION("Calculate absorption coefficients from cross sections.\n"
		    "This calculates both the total absorption and the\n"
		    "absorption per tag group."),
	OUTPUT(	    abs_  , abs_per_tg_                         ),
	INPUT( 	    xsec_per_tg_, vmrs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("xsec_per_tgInit"),
	DESCRIPTION("Initialize xsec_per_tg. The initialization is\n"
		    "necessary, because methods `xsec_per_tgAddLines'\n"
		    "and `xsec_per_tgAddConts' just add to xsec_per_tg.\n"
		    "The size is determined from `tgs'."),
	OUTPUT(	    xsec_per_tg_                             ),
	INPUT(      tgs_, f_mono_, p_abs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("xsec_per_tgAddLines"),
	DESCRIPTION("Calculate cross sections per tag group for line spectra."),
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
	DESCRIPTION("Calculate cross sections per tag group for continua."),
	OUTPUT(	    xsec_per_tg_                             ),
	INPUT( 	    tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_,
		    cont_description_names_, cont_description_parameters_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));


//=== methods operating on absorption ========================================

  md_data.push_back
    ( MdRecord
      ( NAME("abs_per_tgReduce"),
	DESCRIPTION("Reduces absorption coefficients. Only absorption\n"
		    "coefficients for which weighting functions are\n"
		    "calculated are kept in memory."),
	OUTPUT(	    abs_per_tg_ ),
	INPUT( 	    abs_per_tg_, tgs_, wfs_tgs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));


//=== 1D refraction ==========================================================

  md_data.push_back
    ( MdRecord
      ( NAME("refrSet"),
	DESCRIPTION(
           "Sets the refraction input arguments (refr, refr_model and \n"
           "refr_lfac) to the specified values. Type e.g. 'arts -d refr'\n"
           "for more information on the input arguments.\n"
           "See refrCalc for avaliable refraction models.\n"
           "\n"
           "Keywords:\n"
           "     on    : On/off flagg.\n"
           "     model : Model/parametization for the refractive index.\n"
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
          "Sets the refraction flag (refr) to zero and gives other \n"
          "refraction input arguments (refr_lfac and refr_model) dummy \n"
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
           "Calculates the refractive index using the model/parameterization\n"
           "specified by 'refr_model'. \n"
           "If 'refr' is set to zero, the refractive index is set to be an \n"
           "empty vector. \n"
           "\n"
           "Available models are: \n"
           "\n"
           "   'Unity': \n"
           "      Sets the refractive to be 1 at all altitudes. \n"
           "\n"
           "   'Boudouris': \n"
           "      Refractive index at microwave frequencies following \n"
           "      Boudouris 1963. The parameter values were taken from \n"
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
            "altitudes. See the WSV z_tan for definitions."),
	OUTPUT(),
	INPUT( z_tan_, z_plat_ , p_abs_, z_abs_, refr_, refr_index_, r_geoid_,
               z_ground_ ),
	GOUTPUT(Vector_ ),
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
           "Global output: \n"
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
          "Sets the geoid radius to the Earth radius defined in\n"
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
          "Sets the geoid radius according to WGS-84. See the Rodgers book \n"
          "Sec. 9.4.1. The observation direction is given as the angle to\n"
          "the meridian plane (i.e. S=N=0, W=E=90).\n"
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
      ( NAME("groundSet"),
  	DESCRIPTION(
          "Sets the ground altitude and emission to the specified values,\n"
          "and selects a ground temperature.\n"
          "The emission is set to be identical for all frequencies. \n"
          "The ground temperature is obtained by interpolating 't_abs'.\n"
          "\n"
          "Keywords \n"
          "  z : Altitude above the geoid of the ground.\n"
          "  e : Emission factor."),
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
          "Sets the ground emission to the specified value, and sets the \n"
          "altitude and temperature to the first values of z_abs and t_abs.\n"
          "The emission is set to be identical for all frequencies. \n"
          "\n"
          "Keywords \n"
          "  e : Emission factor."),
	OUTPUT( z_ground_, t_ground_, e_ground_ ),
	INPUT( t_abs_, z_abs_, f_mono_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "e"       ),
	TYPES(    Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("groundOff"),
  	DESCRIPTION(
          "Sets the ground altitude and gives the other ground variables\n"
          "(t_ground and e_ground) dummy values. \n"
          "The ground altitude is set to the first element of 'z_abs'.\n"
          "The ground temperature (t_ground) is set to 0. \n"
          "The ground emission vector (e_ground) is set to be empty.\n"
          "If there is a ground intersection and only this function is\n"
          "used to set the ground variables, there will be error messages."),
	OUTPUT( z_ground_, t_ground_, e_ground_ ),
	INPUT( z_abs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("emissionOn"),
  	DESCRIPTION(
	  "Turns on emission by setting the emission flag to 1."),
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
          "Calculates the line-of-sight (LOS) for 1D atmospheres with and\n"
          "without refraction."),
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
          "Calculates source function values valid between the points of"
          "of a 1D LOS.\n" 
          "No scattering and local thermodynamic equilibrium are assumed,\n"
          "that is, the source function equals the Planck function.\n"
          "The source function is set to the mean of the Planck function at\n"
          "the two LOS points limiting the steps. The temperature at the LOS\n"
          "points is obtained by linear interpolation.\n"
          "If emission is neglected (emission=0), the WSV source is set to be"
          "empty."),
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
          "Standard choices for the radiation entering the atmosphere at\n"
          "the start of the LOS. The selections are:\n"
          "  zero : no radiation\n"
          "  cbgr : cosmic background radiation (planck for COSMIC_BG_TEMP)\n"
          "  sun  : solar radiation (planck for SUN_TEMP)\n"
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
          "along the LOS with or without emission.\n"
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
      ( NAME("yTau"),
  	DESCRIPTION(
          "As yCalc but to be used only when emission is neglected.\n"
          "The function returns the optical thicknesses (tau) along the LOS.\n"
          "Some variables needed for yCalc are not needed here (such as \n"
          "y_space, source and t_ground). The emission WSV emission must \n"
          "be set to 0. "),
	OUTPUT( y_ ),
	INPUT( emission_, los_, trans_, e_ground_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("yTB"),
  	DESCRIPTION(
          "Converts a spectrum from intensity to brightness temperature.\n"
          "The used frequency vector is f_sensor."),
	OUTPUT( y_ ),
	INPUT( y_, f_sensor_, za_sensor_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("yTRJ"),
  	DESCRIPTION(
          "Converts a spectrum from intensity to Rayleigh-Jean temperature.\n"
          "The used frequency vetor is f_sensor."),
	OUTPUT( y_ ),
	INPUT( y_, f_sensor_, za_sensor_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixTRJ"),
	DESCRIPTION("Converts a matrix of radiance to Rayleigh-Jeans Brightness\n"
		    "temperatures. This generic method does exatly the\n"
		    "same as the old specific method ybatchTRJ. It can\n"
		    "be used on ybatch or on k-matrices."),
	OUTPUT(),
	INPUT(   f_sensor_, za_sensor_ ),
	GOUTPUT( Matrix_ ),
	GINPUT(  Matrix_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixTB"),
	DESCRIPTION("Converts a matrix of radiance to Planck Brightness\n"
		    "temperatures. This generic method does exatly the\n"
		    "same as the old specific method ybatchTB. It can\n"
		    "be used on ybatch or on k-matrices."),
	OUTPUT(),
	INPUT(   f_sensor_, za_sensor_ ),
	GOUTPUT( Matrix_ ),
	GINPUT(  Matrix_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("yLoadCalibration"),
  	DESCRIPTION(
          "Simulates a load switch calibration as\n"
          "  y = i_cal1 + (i_cal2-i_cal1)*(y-y_cal1)/(y_cal2-y_cal1)\n"
          "The unit for i_cal1,2 can either be intensity or brightness\n"
          "temperature. An example:\n"
          "...\n"
          "VectorPlanck(y_cal1,f_sensor){temp=78}\n"
          "VectorPlanck(y_cal2,f_sensor){temp=300}\n"
          "VectorSet2(i_cal1,y_cal1){value=78}\n"
          "VectorSet2(i_cal2,y_cal2){value=300}\n"
          "yLoadCalibration{}"),
	OUTPUT( y_ ),
	INPUT( y_, i_cal1_, i_cal2_, y_cal1_, y_cal2_, za_sensor_ ),
	GOUTPUT(),
	GINPUT(),
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
	  "be calculated. The specified Strings must be a subgroup of the\n"
          "absorption tag groups (tgs).\n"
	  "Example:\n"
	  "wfs_tgs = [\"O3-666-500e9-501e9, O3-686\",\"H2O\",\"O2-*-*-*\"]"),
	OUTPUT( wfs_tgs_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "wfs_tgs" ),
	TYPES(    Array_String_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("absloswfsCalc"),
  	DESCRIPTION(
          "Calculates absorption line of sight weighting functions (LOS WFs)\n"
          "for 1D atmospheres with or without emission.\n"
          "These WFs are the derivative of the spectra with respect to the \n"
          "absorption at the LOS points. See further the ARTS user guide."),
	OUTPUT( absloswfs_ ),
	INPUT( emission_, los_, source_, trans_, y_, y_space_, f_mono_, 
                                                        e_ground_, t_ground_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("absloswfsTau"),
  	DESCRIPTION(
          "As absloswfsCalc to be used when neglecting emission and \n"
          "thus requires less input (e.g. y_space, source and t_ground are\n"
          "not needed). When using this function the WSV emission must be 0."),
	OUTPUT( absloswfs_ ),
	INPUT( emission_, los_, f_mono_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kSpecies"),
  	DESCRIPTION(
          "Calculates species 1D weighting functions for a single tag.\n"
          "The tag is selected by the giving the tag name. This String must \n"
          "match exactly the String in wfs_tgs. The original absorption\n"
          "array (abs_per_tg) must been reduced to match wfs_tags. \n"
          "The avaliable units are\n"
          "  frac : fractions of linearisation profile \n"
          "  vmr  : volume mixing ratio \n"
          "  nd   : number density\n"
          "\n"
          "Keywords \n"
          "  tag  : Tag String.\n"
          "  unit : Retrieval unit String (see above)."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( los_, absloswfs_, p_abs_, t_abs_, wfs_tgs_, abs_per_tg_, 
                                                              vmrs_, k_grid_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "tag",     "unit"  ),
	TYPES(    String_t,  String_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kSpeciesAll"),
  	DESCRIPTION(
          "Calculates species 1D weighting functions for all tags that\n"
          "are included in wfs_tags. Units as for kSpecies.\n"
          "The original absorption array (abs_per_tg) must been reduced to \n"
          "match wfs_tags. \n"
          "\n"
          "Keywords \n"
          "  unit : Retrieval unit String."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( los_, absloswfs_, p_abs_, t_abs_, wfs_tgs_, abs_per_tg_, 
                                                              vmrs_, k_grid_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "unit"   ),
	TYPES(    String_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kContAbs"),
  	DESCRIPTION(
          "Calculates 1D weighting functions for fit of continuum absorption\n"
          "by polynomials with selectable order.\n"  
          "The continuum is fitted be determining an off-set at a number of\n"
          "points (order+1) that are evenly spread between the lowest and\n"
          "highest frequency of f_mono.\n"
          "\n"
          "Keywords \n"
          "  order : Polynomial order (>=0)."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( los_, absloswfs_, f_mono_, k_grid_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "order" ),
	TYPES(    Index_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("kContAbsSpecifiedLimits"),
  	DESCRIPTION(
          "Calculates 1D weighting functions for fit of continuum absorption\n"
          "by polynomials with selectable order.\n"
          "The continuum is fitted be determining an off-set at a number of\n"
          "points (order+1) that are evenly spread between the given\n"
          "frequency limits.\n"
          "This functions can be used to make seperate fits in the primary\n"
          "and image bands.\n"
          "\n"
          "Keywords \n"
          "  order : Polynomial order (>=0)."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( los_, absloswfs_, f_mono_, k_grid_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "order", "f_low",   "f_high" ),
	TYPES(    Index_t,   Numeric_t, Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kTemp"),
  	DESCRIPTION(
          "Calculates temperature weighting functions with hydrostatic \n"
          "equilibrium. \n"
          "This function includes all effects of a temperature change in a \n"
          "a detailed manner, most notably absorption is recalculated for \n"
          "each point."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_, abs_, 
          lines_per_tg_, lineshape_, e_ground_, emission_, k_grid_, 
          cont_description_names_, cont_description_parameters_,
          z_plat_ ,za_pencil_, l_step_, z_abs_, refr_, refr_lfac_, refr_index_,
	  z_ground_, t_ground_, y_space_, r_geoid_, hse_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kTempFast"),
  	DESCRIPTION(
          "As kTemp but faster as the absorption is assumed to be perfectly\n"
          "linear between t_abs and t_abs+1K.\n"
          "The difference between this function and kTemp should in general\n"
          "be negliable."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_, abs_, 
          lines_per_tg_, lineshape_, e_ground_, emission_, k_grid_, 
          cont_description_names_, cont_description_parameters_,
          z_plat_ ,za_pencil_, l_step_, z_abs_, refr_, refr_lfac_, refr_index_,
	  z_ground_, t_ground_, y_space_, r_geoid_, hse_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kTempNoHydro"),
  	DESCRIPTION(
          "Calculates temperature weighting functions WITHOUT including\n"
          "hydrostatic equilibrium."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( tgs_, los_, absloswfs_, f_mono_, p_abs_, t_abs_, n2_abs_, 
          h2o_abs_, vmrs_, lines_per_tg_, lineshape_, abs_, trans_, e_ground_, 
          k_grid_, cont_description_names_, cont_description_parameters_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kPointingOffSet"),
  	DESCRIPTION(
          "Calculates the WF for a pointing off-set.\n"
          "The functions uses losCalc, transCalc, sourceCalc and yRte to\n"
          "calculate a new spectrum for a changed zenith angles grid.\n"
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
           "expressions (see AUG). With SINGLE_E=0, a WF is returned for \n"
           "the emission coefficient of each monochromatic frequency. \n"
           "On the other hand, when SINGLE_E=1, the ground emission is \n"
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
          "The WF is simply : k = y - y_ref where y_ref is the specified\n"
          "vector."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( y_ ),
	GOUTPUT(),
	GINPUT( Vector_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("kManual"),
  	DESCRIPTION(
          "Calculates a weighting function using y and y0.\n"
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
          "correspondingly. All the data are reallocated to make space for\n"
          "the new data. This function is accordingly slow for large data\n"
          "sizes, and it can be better to use kxAllocate and kxPutInK."),
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
          "correspondingly. All the data are reallocated to make space for\n"
          "the new data. This function is accordingly slow for large data\n"
          "sizes, and it can be better to use kbAllocate and kbPutInK."),
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
          "and kx_aux). The number of frequencies is taken from the length \n"
          "of the given vector (typically y)\n"
          "Use this function before the WF calculations are started and\n"
          "together with kxPutInK.\n"
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
          "Allocates memory for kb and help variables (kb_names, kb_lengths\n"
          "and kb_aux). The number of frequencies is taken from the length \n"
          "of the given vector (typically y)\n"
          "Use this function before the WF calculations are started and\n"
          "together with kbPutInK.\n"
          "\n"
          "Keywords \n"
          "  ni : Number of error identities (temperature profile etc.).\n"
          "  nx : Final length of x."),
	OUTPUT( kb_, kb_names_, kb_lengths_, kb_aux_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT( Vector_ ),
	KEYWORDS( "ni",  "nx"  ),
	TYPES(    Index_t, Index_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kxPutInK"),
  	DESCRIPTION(
          "Puts K in Kx and handles additional data correspondingly.\n"
          "K is placed in the first free columns of Kx.\n"
          "No reallocation is performed (in contrast to kxAppend) and an\n"
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
          "K is placed in the first free columns of Kb.\n"
          "No reallocation is performed (in contrast to kbAppend) and an\n"
          "error message is given if k does not fit into kb. The kb-data are\n"
          "allocated by the function kbAllocate."),
	OUTPUT( kb_, kb_names_, kb_lengths_, kb_aux_ ),
        INPUT( kb_, kb_names_, kb_lengths_, kb_aux_, k_, k_names_, k_aux_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("LinAltsFromPres"),
	DESCRIPTION(
           "Calculates a set of pressures corresponding to a linear set\n"
           "of altitudes, for instance, to set a k_grid with the altitudes\n"
           "equally spaced. The method uses p_abs and z_abs to go an\n"
           "forward and backwards from pressure to altitude. The linear \n"
           "set of altitudes is given with an altitude step and a starting\n"
           "and stopping pressure.\n"
           "Keywords:\n"
           "     delta_z   : altitude step\n"
           "     p_start   : starting pressure\n"
           "     p_stop    : stopping pressure."  ),
	OUTPUT(),
	INPUT( p_abs_, z_abs_ ),
	GOUTPUT(Vector_ ),
	GINPUT(),
	KEYWORDS("delta_z","p_start","p_stop"),
	TYPES(   Numeric_t, Numeric_t, Numeric_t )));




//======================================================================
//=== Methods To Generate Random Data
//======================================================================

  md_data.push_back
    ( MdRecord
      ( NAME("RandSetSeed"),
  	DESCRIPTION(
          "Sets the random seed in a \"random\" way (using the clock).\n"
          "If not this function is called, the same random sequence is\n"
          "obtained for each run."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorRandUniform"),
  	DESCRIPTION(
          "Fills the vector with random data uniformerly distributed between\n"
          "the lower and higher limit given. The length of the vector shall\n"
          "also be given. The random data is uncorrelated."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Vector_ ),
	GINPUT(),
	KEYWORDS( "low",     "high",    "n" ),
	TYPES(    Numeric_t, Numeric_t, Index_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorRandGaussian"),
  	DESCRIPTION(
          "Fills the vector with random data having a normal PDF, zero mean\n"
          "and the standard deviation given. The length of the vector shall\n"
          "also be given. The random data is uncorrelated."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Vector_ ),
	GINPUT(),
	KEYWORDS( "stddev",  "n" ),
	TYPES(    Numeric_t, Index_t )));



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
          "   yRte  \n"
          "and all the workspace variables needed by these functions must be\n"
          "set before starting this method.\n"
          "The refractive index is kept constant for all spectra.\n"
          "The values of the workspace variables are used as defaults when\n"
          "appropiate. For example, if temperature profiles not are read\n"
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
               batchname_, tgs_, cont_description_names_, cont_description_parameters_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS("ncalc", "do_t", "t_file", "do_z", "z_file",
                 "do_tags", "tag_files",
                 "do_f", "f_file", "do_za", "za_file"),
	TYPES(   Index_t,   Index_t,  String_t, Index_t,  String_t, 
                 Array_String_t, Array_String_t,
                 Index_t,  String_t, Index_t,   String_t  )));

#endif // HDF_SUPPORT
}

