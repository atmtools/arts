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
  resize(md_data,0);

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
//=== Overall ARTS functions
//======================================================================

  md_data.push_back     
    ( MdRecord
      ( NAME("Exit"),
	DESCRIPTION("Stops the execution and exits ARTS."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));




//======================================================================
//=== IO methods
//======================================================================

//=== INDEX ============================================================

 // These functions should be changed to handle INDEX

  md_data.push_back     
    ( MdRecord
      ( NAME("IntSet"),
	DESCRIPTION("Sets an integer workspace variable to the given value."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( int_ ),
	GINPUT(),
	KEYWORDS( "value" ),
	TYPES(    int_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("IndexWriteAscii"),
	DESCRIPTION("Writes an integer value to an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( int_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("IndexReadAscii"),
	DESCRIPTION("Reads an integer value from an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( int_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("IndexWriteBinary"),
	DESCRIPTION("Writes an index to a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"
		    "See ??? for details about the file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( int_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("IndexReadBinary"),
	DESCRIPTION("Reads an index from a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( int_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));


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
	GINPUT(  VECTOR_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("NumericWriteAscii"),
	DESCRIPTION("Writes a numeric value to an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( Numeric_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("NumericReadAscii"),
	DESCRIPTION("Reads a numeric value from an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( Numeric_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("NumericWriteBinary"),
	DESCRIPTION("Writes a numeric value to a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab."),
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
		    "If empty, it is set to <basename>.<variable_name>.ab."),
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
	GOUTPUT( VECTOR_ ),
	GINPUT(),
	KEYWORDS( "length", "value"   ),
	TYPES(    int_t,    Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorSetLengthFromVector"),
	DESCRIPTION("Creates a workspace vector with the same length as the\n"
                    "given vector and initializes the new vector with the\n"
                    "given value. For example\n"
                    " VectorSetLengthFromVector(e_ground,f_mono){value=0.75}"),
	OUTPUT(),
	INPUT(),
	GOUTPUT( VECTOR_ ),
	GINPUT( VECTOR_ ),
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
	GOUTPUT( VECTOR_ ),
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
	GOUTPUT(VECTOR_),
	GINPUT(),
	KEYWORDS( "start",   "stop",    "n"   ),
	TYPES(    Numeric_t, Numeric_t, int_t )));

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
	KEYWORDS( "start",   "stop",    "n"   ),
	TYPES(    Numeric_t, Numeric_t, int_t )));

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
      ( NAME("VectorCopyFromArrayOfVector"),
	DESCRIPTION("Copies a vector from a vector array.\n"
          "\n"
          "Keywords \n"
          "  index : The index of the vector in the array to copy. " ),
	OUTPUT(),
	INPUT(),
	GOUTPUT( VECTOR_ ),
	GINPUT( ARRAYofVECTOR_),
	KEYWORDS( "index" ),
	TYPES( int_t )));

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
		    "If empty, it is set to <basename>.<variable_name>.ab."),
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
		    "If empty, it is set to <basename>.<variable_name>.ab."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( VECTOR_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorPlanck"),
  	DESCRIPTION(
          "Sets a vector to the Planck function for the given frequency\n"
          "vector and temperature. An example:\n"
          "   VectorPlanck(y_space,f_mono){temp=2.7}"),
	OUTPUT( ),
	INPUT( ),
	GOUTPUT( VECTOR_ ),
	GINPUT( VECTOR_ ),
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
	GOUTPUT( VECTOR_ ),
	GINPUT( VECTOR_ ),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorAdd"),
  	DESCRIPTION(
	  "Adds a scalar to all elements of a vector."),
	OUTPUT( ),
	INPUT( ),
	GOUTPUT( VECTOR_ ),
	GINPUT( VECTOR_ ),
	KEYWORDS( "value" ),
	TYPES( Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorScale"),
  	DESCRIPTION(
	  "Scales all elements of a vector with the same value."),
	OUTPUT( ),
	INPUT( ),
	GOUTPUT( VECTOR_ ),
	GINPUT( VECTOR_ ),
	KEYWORDS( "value" ),
	TYPES( Numeric_t )));



//=== MATRIX ==========================================================

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixSet"),
	DESCRIPTION("Creates a workspace matrix of the specified size\n"
                    "and initializes the matrix with the given value."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( MATRIX_ ),
	GINPUT(),
	KEYWORDS( "nrows", "ncols", "value"   ),
	TYPES(    int_t,   int_t,   Numeric_t )));

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
      ( NAME("MatrixFillWithVector"),
	DESCRIPTION("Forms a vector with n columns, and put the given.\n"
                    "in each column."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( MATRIX_ ),
	GINPUT( VECTOR_ ),
	KEYWORDS( "n"   ),
	TYPES(    int_t )));

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
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"),
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
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"),
	OUTPUT(),
	INPUT(),
	GOUTPUT( MATRIX_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixScale"),
	DESCRIPTION(
          "Scales all elements of a matrix with the same value."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( MATRIX_ ),
	GINPUT( MATRIX_ ),
	KEYWORDS( "value" ),
	TYPES(    Numeric_t   )));



//=== ARRAYofINDEX =====================================================

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexWriteAscii"),
	DESCRIPTION("Writes an integer array to an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( ARRAYofsizet_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexReadAscii"),
	DESCRIPTION("Reads an integer array from an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( ARRAYofsizet_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexWriteBinary"),
	DESCRIPTION("Writes an index array to a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab."),
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
		    "If empty, it is set to <basename>.<variable_name>.ab."),
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
		    "If empty, it is set to <basename>.<variable_name>.ab."),
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
		    "If empty, it is set to <basename>.<variable_name>.ab."),
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
		    "If empty, it is set to <basename>.<variable_name>.ab."),
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
		    "If empty, it is set to <basename>.<variable_name>.ab."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( ARRAYofMATRIX_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));


//=== STRING ============================================================

  md_data.push_back
    ( MdRecord
      ( NAME("StringSet"),
	DESCRIPTION("Sets a string to the given text string."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( string_ ),
	GINPUT(),
	KEYWORDS( "text"   ),
	TYPES(    string_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("StringWriteAscii"),
	DESCRIPTION("Writes a string to an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.as.\n"
		    "See `ArrayOfStringWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( string_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("StringReadAscii"),
	DESCRIPTION("Reads a string from an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.as.\n"
		    "See `ArrayOfStringWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( string_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("StringWriteBinary"),
	DESCRIPTION("Writes a string to a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( string_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("StringReadBinary"),
	DESCRIPTION("Reads a string from a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( string_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));


//=== ARRAYofSTRING =========================================================

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfStringSet"),
	DESCRIPTION("Sets a string array according the given text.\n"
                    "The format is text = [\"string1\",\"string2\",...]"),
	OUTPUT(),
	INPUT(),
	GOUTPUT( ARRAYofstring_ ),
	GINPUT(),
	KEYWORDS( "text"         ),
	TYPES(    ARRAY_string_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfStringWriteAscii"),
	DESCRIPTION("Writes a string array to an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.as.\n"
		    "The format is as follows:\n\n"
		    "# <comments>\n"
		    "<n_strings>\n"
		    "<string 1>\n"
		    "<string 2>\n"
		    "...\n\n"
		    "Example:\n"
		    "# Generated by arts-0.0.16, Apr 29 2000, 17:38:44\n"
		    "2\n"
		    "A string\n"
		    "Another string"),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( ARRAYofstring_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfStringReadAscii"),
	DESCRIPTION("Reads a string array from an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.as.\n"
		    "See `ArrayOfStringWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( ARRAYofstring_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfStringWriteBinary"),
	DESCRIPTION("Writes a string array to a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( ARRAYofstring_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfStringReadBinary"),
	DESCRIPTION("Reads a string array from a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( ARRAYofstring_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));



//=== SYMMETRIC ==========================================================

  md_data.push_back
    ( MdRecord
      ( NAME("SymmetricWriteAscii"),
	DESCRIPTION("Writes a covariance matrix to an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( SYMMETRIC_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("SymmetricReadAscii"),
	DESCRIPTION("Reads a covariance matrix from an ASCII file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.am.\n"
		    "See `ArrayOfMatrixWriteAscii' for file format."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( SYMMETRIC_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("SymmetricWriteBinary"),
	DESCRIPTION("Writes a covariance matrix to a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( SYMMETRIC_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("SymmetricReadBinary"),
	DESCRIPTION("Reads a covariance matrix from a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab.\n"),
	OUTPUT(),
	INPUT(),
	GOUTPUT( SYMMETRIC_ ),
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



//=== LOS ==================================================================

  md_data.push_back
    ( MdRecord
      ( NAME("LosWriteBinary"),
	DESCRIPTION("Writes a LOS structure to a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT( LOS_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("LosReadBinary"),
	DESCRIPTION("Reads a LOS structure from a binary file.\n"
		    "The filename can be specified or an empty string.\n"
		    "If empty, it is set to <basename>.<variable_name>.ab."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( LOS_ ),
	GINPUT(),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));



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
	TYPES(    ARRAY_string_t, ARRAY_string_t, VECTOR_t, VECTOR_t)));

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
	TYPES(    string_t,    Numeric_t, Numeric_t)));

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
	TYPES(    string_t,    Numeric_t, Numeric_t)));

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
	  "entry at frequency -f is added to lines_per_tg."),
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
	  "from the lines_per_tg, in order to save computation time"),
	OUTPUT(   lines_per_tg_      ),
	INPUT(    lines_per_tg_, lineshape_, f_mono_  ),
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
      ( NAME("tgsDefine"),
  	DESCRIPTION(
          "Set up the list of tag groups.\n"
	  "Specify one string for each tag group that you want to create.\n"
	  "Inside the string, separate the tags by comma (plus optional blanks).\n"
	  "Example:\n"
	  "tag = [\"O3-666-500e9-501e9, O3-686\",\"H2O\",\"O2-*-*-*\"]"),
	OUTPUT( tgs_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "tags" ),
	TYPES(    ARRAY_string_t   )));

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
	TYPES(     string_t,        string_t,         Numeric_t )));

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
	TYPES(   ARRAY_string_t,         ARRAY_string_t,        VECTOR_t )));


//=== Continuum methods ============================================

  md_data.push_back
    ( MdRecord
      ( NAME("cont_descriptionInit"),
  	DESCRIPTION(
          "Initializes the two continuum description WSVs,\n"
	  "`cont_description_names' and `cont_description_parameters'.\n"  
	  "\n"
	  "This method does not really do anything, except setting the two\n"
	  "variables to empty ARRAYs. It is just necessary\n"
	  "because the method `cont_descriptionAppend' wants to append to the\n"
	  "variables.\n"
	  "\n"
	  "Formally, the continuum description WSVs are required by the\n"
	  "absorption calculation methods (e.g., `absCalc'). Therefore you\n"
	  "always have to call `cont_descriptionInit'.\n"
	  "\n"
	  "Usage example: cont_descriptionInit{}"),
	OUTPUT( cont_description_names_, cont_description_parameters_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("cont_descriptionAppend"),
  	DESCRIPTION(
          ""),
	OUTPUT( cont_description_names_, cont_description_parameters_ ),
	INPUT(  cont_description_names_, cont_description_parameters_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "name",   "parameters" ),
	TYPES(    string_t, VECTOR_t     )));


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
	INPUT(    tgs_          ),
	GOUTPUT(                       ),
	GINPUT(                        ),
	KEYWORDS( "basename"           ),
	TYPES(    string_t             )));

//   md_data.push_back
//     ( MdRecord
//       ( NAME("Atm2dFromRaw1D"),
//   	DESCRIPTION(
// 	  "This method is not currently useful for anything, since\n"
// 	  "there is no method to calculate absorption from the 2D\n"
// 	  "parameters.\n"
// 	  "\n"
// 	  "Interpolates temperature, altitude, and VMRs to the pressure grid\n"
// 	  "given by p_abs. The altitude is not used by the absorption routines,\n"
// 	  "But later on by the RT routines."
// 	  "\n"
// 	  "Interpolations used: FIXME: Add these.f\n"
// 	  "Temperature [K]: \n"
// 	  "Altitude    [m]: \n"
// 	  "VMRs        [1]: \n"
// 	  "\n"
// 	  "Uses interp_lin(...)."
// 	  ),
// 	OUTPUT(   t_abs_2d_ , z_abs_2d_   , vmrs_2d_     ),
// 	INPUT(    p_abs_    , raw_ptz_1d_ , raw_vmrs_1d_ ),
// 	GOUTPUT(                       			 ),         
// 	GINPUT(                        			 ),
// 	KEYWORDS(                      			 ),
// 	TYPES(                         			 )));

  md_data.push_back
    ( MdRecord
      ( NAME("AtmFromRaw1D"),
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
	INPUT(    p_abs_    , raw_ptz_1d_ , raw_vmrs_1d_ ),
	GOUTPUT(                       			 ),         
	GINPUT(                        			 ),
	KEYWORDS(                      			 ),
	TYPES(                         			 )));

  md_data.push_back
    ( MdRecord
      ( NAME("z_absHydrostatic"),
	DESCRIPTION(
          "Calculates a vertical grid fulfilling hydrostatic equilibrium\n"
          "taking effects from water vapour into account.\n"
          "The given altitudes are used as a first guess when starting the\n"
          "calculations (to estimate g etc.).\n"
          "The altitude of one pressure level must be given (the reference \n"
          "point).\n"
          "The altitude variation of the gravitational acceleration is\n"
          "considered.\n"
          "The average molecular weight is assumed to be 28.96 at all\n"
          "altitudes \n"
          "\n"
          "Keywords \n"
          "  g0    : Gravitational acceleration at the geoid surface.\n"
          "  pref  : Pressure reference point.\n"
          "  zref  : The geometrical altitude at pref.\n"
          "  niter : Number of iterations (1-2 should suffice normally)."),
	OUTPUT( z_abs_ ),
	INPUT( z_abs_, p_abs_, t_abs_, h2o_abs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "g0",      "pref",    "zref",    "niter" ),
	TYPES(    Numeric_t, Numeric_t, Numeric_t, int_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("h2o_absSet"),
	DESCRIPTION(
          "Sets h2o_abs to the profile of the first tag group containing\n"
	  "water."),
	OUTPUT(	    h2o_abs_ ),
	INPUT( 	tgs_, vmrs_  ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));




//=== 1D absorption methods ===============================================

  md_data.push_back
    ( MdRecord
      ( NAME("absCalc"),
	DESCRIPTION("Calculate absorption coefficients. This\n"
		    "calculates both the total absorption and the\n"
		    "absorption per tag group."),
	OUTPUT(	    abs_  , abs_per_tg_                         ),
	INPUT( 	    tgs_, f_mono_, p_abs_, t_abs_, h2o_abs_, vmrs_, lines_per_tg_, 
		    lineshape_,
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
	INPUT( 	    tgs_, f_mono_, p_abs_, t_abs_, h2o_abs_, vmrs_,
		    cont_description_names_, cont_description_parameters_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("refr_indexBoudourisDryAir"),
	DESCRIPTION(
           "Calculates the refractive index for dry air at micro-wave\n"
	   "frequencies following Boudouris 1963.\n"
           "The effect of water vapor is neglected (dry air).\n"
	   "The expression is also found in Chapter 5 of the Janssen book."),
	OUTPUT( refr_index_ ),
	INPUT( 	p_abs_, t_abs_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("refr_indexBoudouris"),
	DESCRIPTION(
           "Calculates the refractive index at microwave frequencies\n"
           "following Boudouris 1963.\n"
	   "The expression is also found in Chapter 5 of the Janssen book."),
	OUTPUT(	refr_index_ ),
	INPUT( 	p_abs_, t_abs_, h2o_abs_ ),
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



//======================================================================
//=== LOS/RTE methods
//======================================================================


  md_data.push_back
    ( MdRecord
      ( NAME("zaFromZtan"),
	DESCRIPTION("Calculates the zenith angles from a set of tangent\n"
                    "altitudes and a given LOS geometry."),
	OUTPUT(),
	INPUT( z_tan_, z_plat_ , p_abs_, z_abs_, refr_, refr_index_, r_geoid_, z_ground_ ),
	GOUTPUT(VECTOR_ ),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

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
      ( NAME("NoGround"),
  	DESCRIPTION(
          "Sets the ground altitude and gives the other ground variables\n"
          "(t_ground and e_ground) dummy values. \n"
          "The ground altitude is set to the specified value. Note that\n"
          "z_abs must cover the ground altitude.\n"
          "The ground temperature (t_ground) is set to 0. \n"
          "The ground emission vector (e_ground) is set to be empty.\n"
          "If there is a ground intersection and only this function is\n"
          "used to set the ground variables, there will be error messages."),
	OUTPUT( z_ground_, t_ground_, e_ground_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "z_ground" ),
	TYPES(    Numeric_t  )));

  md_data.push_back
    ( MdRecord
      ( NAME("NoRefraction"),
  	DESCRIPTION(
          "Sets the refraction flag (refr) to zero and gives other \n"
          "refraction variables (refr_index and refr_lfac) dummy values. \n"
          "If the refraction later is turned on, the variables refr_index \n"
          "and refr_lfac must be given propoer values, or there will be\n"
          "error messages."),
	OUTPUT( refr_, refr_index_, refr_lfac_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

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
          "Standard choices for the radiation entering the atmosphere at\n"
          "the start of the LOS. The selections are:\n"
          "  zero : no radiation\n"
          "  cbgr : cosmic background radiation (planck for COSMIC_BG_TEMP)\n"
          "  sun  : solar radiation (planck for SUN_TEMP)\n"
          "COSMIC_BG_TEMP and SUN_TEMP are global variables, defined in\n"
          "constants.cc.\n"
          "\n"
          "Keywords \n"
          "  choice : Selection string (see above)."),
	OUTPUT( y_space_ ),
	INPUT( f_mono_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "nr"  ),
	TYPES(    string_t )));

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
	  "be calculated. The specified strings must be a subgroup of the\n"
          "absorption tag groups (tgs).\n"
	  "Example:\n"
	  "wfs_tgs = [\"O3-666-500e9-501e9, O3-686\",\"H2O\",\"O2-*-*-*\"]"),
	OUTPUT( wfs_tgs_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "wfs_tgs" ),
	TYPES(    ARRAY_string_t   )));

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
      ( NAME("kSpecies"),
  	DESCRIPTION(
          "Calculates species 1D weighting functions for a single tag.\n"
          "The tag is selected by the giving the tag name. This string must \n"
          "match exactly the string in wfs_tgs. The original absorption\n"
          "array (abs_per_tg) must been reduced to match wfs_tags. \n"
          "The avaliable units are\n"
          "  frac : fractions of linearisation profile \n"
          "  vmr  : volume mixing ratio \n"
          "  nd   : number density\n"
          "\n"
          "Keywords \n"
          "  tag  : Tag string.\n"
          "  unit : Retrieval unit string (see above)."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( los_, absloswfs_, p_abs_, t_abs_, wfs_tgs_, abs_per_tg_, 
                                                              vmrs_, k_grid_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "tag",     "unit"  ),
	TYPES(    string_t,  string_t )));

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
          "  unit : Retrieval unit string."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( los_, absloswfs_, p_abs_, t_abs_, wfs_tgs_, abs_per_tg_, 
                                                              vmrs_, k_grid_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "unit"   ),
	TYPES(    string_t )));

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
	TYPES(    int_t   )));

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
	TYPES(    int_t,   Numeric_t, Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kTempNoHydro"),
  	DESCRIPTION(
          "Calculates temperature 1D weighting functions WITHOUT including\n"
          "hydrostatic equilibrium."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( tgs_, los_, absloswfs_, f_mono_, p_abs_, t_abs_, 
               h2o_abs_, vmrs_, lines_per_tg_, lineshape_, 
               abs_, trans_, e_ground_, k_grid_,
	       cont_description_names_, cont_description_parameters_ ),
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
               abs_, y_space_, e_ground_, t_ground_, y_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "delta"   ),
	TYPES(    Numeric_t )));

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
	GINPUT( VECTOR_ ),
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
	TYPES(    string_t, Numeric_t, Numeric_t, Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("kDiffHSmall"),
  	DESCRIPTION(
          "Calculates a weighting function using y, h1 and h2.\n"
          "This function minimizes memory usage. For faster calculations,\n"
          "use kDiffHFast.\n"
          "The weighting function is calculated as: k = (h2*y-h1*y)/delta\n"
          "That is, delta is the magnitude of the perturbation done.\n"
          "Note that the obtained k matrix includes effects of h1.\n"
          "\n"
          "Keywords \n"
          "  name    : Name on retrieval/error identity.\n"
          "  delta   : Magnitude of perturbation.\n"
          "  grid    : Grid point value.\n"
          "  apriori : A priori value."),
	OUTPUT( k_, k_names_, k_aux_ ),
	INPUT( h1_, h2_, y_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "name",   "delta",   "grid",   "apriori"  ),
	TYPES(    string_t, Numeric_t, Numeric_t, Numeric_t )));

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
	KEYWORDS( "name",   "delta",   "grid",   "apriori"  ),
	TYPES(    string_t, Numeric_t, Numeric_t, Numeric_t )));

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
	GINPUT( VECTOR_ ),
	KEYWORDS( "ni",  "nx"  ),
	TYPES(    int_t, int_t )));

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
	GINPUT( VECTOR_ ),
	KEYWORDS( "ni",  "nx"  ),
	TYPES(    int_t, int_t )));

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



//======================================================================
//=== Covariance Matrix Methods
//======================================================================

  md_data.push_back
    ( MdRecord
      ( NAME("sDiagonal"),
  	DESCRIPTION(
          "Creates a diagonal covariance matrix.\n"
          "\n"
          "Keywords \n"
          "       n : Size of covariance matrix.\n"
          "  stddev : Standard deviation."),
	OUTPUT( s_ ),
        INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "n",   "stddev"  ),
	TYPES(    int_t, Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("sDiagonalLengthFromVector"),
  	DESCRIPTION(
          "Creates a diagonal covariance matrix matching a vector.\n"
          "The size of s is determined of the length of the given vector.\n"
          "\n"
          "Keywords \n"
          "  stddev : Standard deviation."),
	OUTPUT( s_ ),
        INPUT(),
	GOUTPUT(),
	GINPUT( VECTOR_ ),
	KEYWORDS( "stddev"  ),
	TYPES(    Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("sDiagonalLengthFromVectors"),
  	DESCRIPTION(
          "Creates a diagonal covariance matrix matching 2 vectors.\n"
          "This function can be used to create a covariance matrix for all\n"
          "combinations between the values in the two vectors. For example,\n"
          "   sDiagonalLengthFromVectors(f_sensor,za_sensor){2}\n"
          "gives a covariance matrix for uncorrelated thermal noise with a\n"
          "magnitude of 2 (probably Kelvins).\n"
          "The size of s is the product of the length of the given vectors.\n"
          "\n"
          "Keywords \n"
          "  stddev : Standard deviation."),
	OUTPUT( s_ ),
        INPUT(),
	GOUTPUT(),
	GINPUT( VECTOR_, VECTOR_ ),
	KEYWORDS( "stddev"  ),
	TYPES(    Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("sSimple"),
  	DESCRIPTION(
          "Creates a covariance matrix where the standard deviation and\n"
          "correlation length are constant.\n"
          "The abscissa is set to 1,2,3,...,n. The correlation length shall\n"
          "accordingly be treated as a index distance.\n"
          "\n"
          "Keywords \n"
          "  n          : Size of covariance matrix.\n"
          "  stddev     : Standard deviation.\n"
          "  corrfun    : The correlation function: \n"
          "               0  no correlation, diagonal matrix. The corrlength\n"
          "                  and cutoff are of no importance here.\n"
          "               1  linearly decreasing to 0 (tenth function) \n"
          "               2  exponential \n"
          "               3  gaussian \n"
          "  cutoff     : Correlations below this value are set to 0.\n"
          "               This variable can be used to make s more sparse.\n"
          "  corrlength : Correlation length (in index). The length where\n"
          "               the correlation has decreased to exp(-1)."),
	OUTPUT( s_ ),
        INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "n",   "stddev",  "corrfun", "cutoff",  "corrlength" ),
	TYPES(    int_t, Numeric_t, int_t,     Numeric_t, Numeric_t  )));

  md_data.push_back
    ( MdRecord
      ( NAME("sSimpleLengthFromVector"),
  	DESCRIPTION(
          "Creates a covariance matrix where the standard deviation and\n"
          "correlation length are constant. The size of s is determined of\n"
          "the length of the given vector.\n"
          "The given vector is also used as abscissa when calulating \n"
          "correlation values. \n"
          "\n"
          "Keywords \n"
          "  stddev     : Standard deviation.\n"
          "  corrfun    : The correlation function. See sSimple. \n"
          "  cutoff     : Correlations below this value are set to 0.\n"
          "               This variable can be used to make s more sparse.\n"
          "  corrlength : Correlation length. The length where the corre-\n"
          "               lation has decreased to exp(-1)."),
	OUTPUT( s_ ),
        INPUT(),
	GOUTPUT(),
	GINPUT( VECTOR_ ),
	KEYWORDS( "stddev",  "corrfun", "cutoff",  "corrlength" ),
	TYPES(    Numeric_t, int_t,     Numeric_t, Numeric_t  )));

  md_data.push_back
    ( MdRecord
      ( NAME("sSimpleLengthFromVectors"),
  	DESCRIPTION(
          "Creates a covariance matrix where the standard deviation and\n"
          "correlation length are constant.\n"
          "The size of s is the product of the length of the two given\n"
          "vectors. See further sDiagonalLengthFromVectors.\n"
          "The vector corresponding to s is assumed to have the structure\n"
          "   [vector1,vector1,...,vector1]\n"
          "where vector1 is the first of the two input vectors. The number \n"
          "of repititions of vector1 equals the length of the second vector.\n"
          "The correlation between each repitition of vector 1 is set to 0.\n"
          "The total covariance matrix has accordingly the structure:\n"
          "   s1  0  0 \n"
          "    0 s1  0 \n"
          "    0  0 s1 \n"
          "where s1 is the matrix given by sSimpleLengthFromVector for\n"
          "vector 1.\n"
          "\n"
          "Keywords \n"
          "  stddev     : Standard deviation.\n"
          "  corrfun    : The correlation function. See sSimple. \n"
          "  cutoff     : Correlations below this value are set to 0.\n"
          "               This variable can be used to make s more sparse.\n"
          "  corrlength : Correlation length in units of the first of the\n"
          "               given vectors. This is the length where the corre-\n"
          "               lation has decreased to exp(-1)."),
	OUTPUT( s_ ),
        INPUT(),
	GOUTPUT(),
	GINPUT( VECTOR_, VECTOR_ ),
	KEYWORDS( "stddev",  "corrfun", "cutoff",  "corrlength" ),
	TYPES(    Numeric_t, int_t,  Numeric_t, Numeric_t  )));

  md_data.push_back
    ( MdRecord
      ( NAME("sFromFile"),
  	DESCRIPTION(
          "Creates a covariance matrix based on definition data in a file.\n"
          "The abscissa of the covariance matrix is determined by the given\n"
          "vector.\n"
          "Only the ASCII files are allowed.\n"
          "The covaraince matrix can be constructed as a sum of an arbitrary\n"
          "number of covariance matrices.\n"
          "The definition data shall be stored as a ARRAYofMATRIX.\n"
          "The first row of the first matrix gives the correlation function\n"
          "flag for each covariance matrix. See sSimple.\n"
          "The second row of the first matrix gives the correlation cut-off\n"
          "for each covariance matrix. See sSimple.\n"
          "The number of columns of the first matrix must equal the length\n"
          "of the matrix array - 1.\n"
          "Each covariance matrix is defined by a three column matrix, where\n"
          "column 1 is the abscissa (matching the given vector), column 2 \n"
          "standard deviations and column 3 correlation lengths.\n"
          "Linear interpolation is used to get values at intermediate\n"
          "points.\n"
          "An example, where a diagonal matrix and a gaussian matrix, both\n"
          "having a standard deviation of 1, are summed, the gaussian\n"
          "correlation length is 2 and k_grid has values between 0 and 10:\n"
          "3           \n"
          "2 2         \n"
          "0 2         \n"
          "0 0         \n"
          "2 3         \n"
          "0   1  0    \n"
          "10  1  0    \n"
          "2 3         \n"
          "0   1  2    \n"
          "10  1  2    \n"
          "\n"
          "Keywords \n"
          "  filename : Name on file. The file name must be specified, no\n"
          "             default name exists."),
	OUTPUT( s_ ),
        INPUT(),
	GOUTPUT(),
	GINPUT( VECTOR_ ),
	KEYWORDS( "filename" ),
	TYPES(    string_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("CovmatrixInit"),
  	DESCRIPTION(
          "Initializes a covariance matrix.\n"
          "The matrix is set to be empty (0 x 0)."),
	OUTPUT(),
        INPUT(),
	GOUTPUT( SYMMETRIC_ ),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("sxAppend"),
  	DESCRIPTION(
          "Appends s to sx assuming no correlation.\n"
          "CovmatrixInit(sx){} must be called before using this function for\n"
          "the first time."),
	OUTPUT( sx_ ),
        INPUT( sx_, s_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("sbAppend"),
  	DESCRIPTION(
          "Appends s to sb assuming no correlation.\n"
          "CovmatrixInit(sb){} must be called before using this function for\n"
          "the first time."),
	OUTPUT( sb_ ),
        INPUT( sb_, s_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));



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
	GOUTPUT( VECTOR_ ),
	GINPUT(),
	KEYWORDS( "low",     "high",    "n" ),
	TYPES(    Numeric_t, Numeric_t, int_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorRandGaussian"),
  	DESCRIPTION(
          "Fills the vector with random data having a normal PDF, zero mean\n"
          "and the standard deviation given. The length of the vector shall\n"
          "also be given. The random data is uncorrelated."),
	OUTPUT(),
	INPUT(),
	GOUTPUT( VECTOR_ ),
	GINPUT(),
	KEYWORDS( "stddev",  "n" ),
	TYPES(    Numeric_t, int_t )));



//======================================================================
//=== Batch Calculation Methods
//======================================================================

  md_data.push_back
    ( MdRecord
      ( NAME("BatchdataGaussianNoiseNoCorrelation"),
  	DESCRIPTION(
          "Creates a set of noise vectors suitable for batch calculations.\n"
          "The produced vectors have elements taken from a Gaussia\n"
          "distribution of zero mean and stddev, simulating the thermal\n"
          "noise of a sensor with no inter-channel correlation.\n"
          "The length of the vectors is the product of the sizes of the two\n"
          "given vectors, typically f_sensor, za_sensor.\n"
          "The vectors are written to the file batchname.noise.ab  \n"
          "\n"
          "Keywords \n"
          "  n       : Number of random vectors to produce.\n"
          "  stddev  : Standard deviation.\n"),
	OUTPUT(),
	INPUT( batchname_ ),
	GOUTPUT(),
	GINPUT( VECTOR_, VECTOR_ ),
	KEYWORDS( "n", "stddev" ),
	TYPES(    int_t, Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("BatchdataGaussianZeroMean"),
  	DESCRIPTION(
          "Creates a set of vectors suitable for batch calculations.\n"
          "The produced vectors have zero mean and its statistics match the\n"  
          "given covariance matrix.\n"
          "The vectors are written to the file\n"
          "   batchname.varname.ab  \n"
          "where varname is the given keyword string.\n"
          "\n"
          "Keywords \n"
          "  n       : Number of random vectors to produce.\n"
          "  varname : Variable name."),
	OUTPUT(),
	INPUT( batchname_, s_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "n",   "varname"  ),
	TYPES(    int_t, string_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("BatchdataGaussianTemperatureProfiles"),
  	DESCRIPTION(
          "Creates temperature profiles suitable for batch calculations.\n"
          "The profiles have gaussian statistics and are stored in a file.\n"
          "The mean of the data equals t_abs. Standard deviations and \n"
          "correlations follow the covariance matrix s.\n"
          "The corresponding altitudes are calculated by z_absHydrostatic\n"
          "and are stored in a seperate file.\n"
          "The data are stored to the files:\n"
          "   batchname.t_abs.ab  \n"
          "   batchname.z_abs.ab  \n"
          "\n"
          "Keywords \n"
          "  n        : Number of random vectors to produce.\n"
          "  g0       : Gravitational acceleration at the geoid surface.\n"
          "  pref     : Pressure reference point.\n"
          "  zref     : The geometrical altitude at pref.\n"
          "  niter    : Number of iterations (1-2 should suffice normally)."),
	OUTPUT(),
	INPUT( p_abs_, t_abs_, z_abs_, h2o_abs_, s_, batchname_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "n",   "g0",      "pref",    "zref",    "niter" ),
	TYPES(    int_t, Numeric_t, Numeric_t, Numeric_t, int_t   )));

  md_data.push_back
    ( MdRecord
      ( NAME("BatchdataGaussianTemperatureProfilesNoHydro"),
  	DESCRIPTION(
          "As BatchdataGaussianTemperatureProfiles but no altitudes are\n"
          "calculated.\n"
          "The data are stored to the file:\n"
          "   batchname.t_abs.ab  \n"
          "\n"
          "Keywords \n"
          "  n        : Number of random vectors to produce."),
	OUTPUT(),
	INPUT( t_abs_, s_, batchname_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "n"   ),
	TYPES(    int_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("BatchdataGaussianSpeciesProfiles"),
  	DESCRIPTION(
          "Creates species profiles suitable for batch calculations.\n"
          "The profiles have gaussian statistics and are stored in a file.\n"
          "The mean of the data equals the corresponding profile in vmrs.\n"
          "Standard dev. and correlations follow the covariance matrix s.\n"
          "Profiles for several species can be produced if s is valid for\n"
          "all the species.\n"
          "The data are stored to the files:\n"
          "   batchname.XXX.ab  \n"
          "where XX is the molecule name.\n"
          "\n"
          "Keywords \n"
          "  n       : Number of random vectors to produce.\n"
          "  do_tags : This string array gives the tags for which profiles\n"
          "            shall be generated, e.g. [\"H2O\",\"O3\"].\n"
          "            These tags must match some tag in tags.\n"
          "  unit    : Unit string for the given standard deviation. Unit\n"
          "            coding as for kSpecies."),
	OUTPUT(),
	INPUT( tgs_, vmrs_, p_abs_, t_abs_, s_, batchname_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "n",   "do_tags",      "unit"   ),
	TYPES(    int_t, ARRAY_string_t, string_t    )));

  md_data.push_back
    ( MdRecord
      ( NAME("BatchdataGaussianOffSets"),
  	DESCRIPTION(
          "Creates a set of vectors suitable for batch calculations.\n"
          "The produced vectors are shifted, having off-sets, compared\n"
          "to the given vector where the off-sets have gaussian PDF.\n"
          "There is no correlation between the off-sets.\n"
          "The vectors are written to the file\n"
          "   batchname.varname.ab  \n"
          "where varname is the workspace name of the reference vector.\n"
          "\n"
          "Keywords \n"
          "  n        : Number of random vectors to produce.\n"
          "  stddev   : Standard deviation for the off-sets."),
	OUTPUT(),
	INPUT( batchname_ ),
	GOUTPUT(),
	GINPUT( VECTOR_ ),
	KEYWORDS( "n",   "stddev"  ),
	TYPES(    int_t, Numeric_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("BatchdataUniformOffSets"),
  	DESCRIPTION(
          "Creates a set of vectors suitable for batch calculations.\n"
          "The produced vectors are shifted, having off-sets, compared\n"
          "to the given vector where the off-sets have uniform PDF\n"
          "between low and high.\n"
          "Note that the off-set range (low to high) can be strictly\n"
          "negative or positive.\n"
          "There is no correlation between the off-sets.\n"
          "The vectors are written to the file\n"
          "   batchname.varname.ab  \n"
          "where varname is the workspace name of the reference vector.\n"
          "\n"
          "Keywords \n"
          "  n    : Number of random vectors to produce.\n"
          "  low  : Lowest limit for the off-sets.\n"
          "  high : Upper limit for the off-sets." ),
	OUTPUT(),
	INPUT( batchname_ ),
	GOUTPUT(),
	GINPUT( VECTOR_ ),
	KEYWORDS( "n",   "low",    "high"     ),
	TYPES(    int_t, Numeric_t, Numeric_t )));


  md_data.push_back
    ( MdRecord
      ( NAME("BatchdataSinusoidalRippleNoCorrelations"),
  	DESCRIPTION(
          "Creates a batch of sinusoidal baseline ripples. \n"
          "The ripple is modelled as\n"
          "   a*sin(2*pi*f/period+phase)  \n"
          "The frequency period is kept constant to the selected value.\n"
          "The phase is given a uniform PDF between 0 and 2*pi.\n"
          "The amplitude can be modelled in three different ways by the\n"
          "keyword pdf. The choices for pdf are:\n"
          "  none     : The amplitude is kept constant to the selected value\n"
          "  gaussian : The selected value is used as standard deviation.\n"
          "  uniform  : Uniform PDF between -+amplitude.\n"
          "The phase and amplitude variations are uncorrelated, both\n"
          "mutually and between different spectra.\n"
          "The first input vector is used as frequency vector.\n"
          "The second input vector is treated as the zenith angle grid.\n"
          "The vectors are written to the file\n"
          "   batchname.varname.ab  \n"
          "where varname is the given keyword string.\n"
          "\n"
          "Keywords \n"
          "  n         : Number of random vectors to produce.\n"
          "  period    : The frequency period.\n"
          "  amplitude : Amplitude of the ripple (see above).\n"
          "  pdf       : Probability density function. Possible choices:\n"
          "                \"none\", \"gaussian\" and \"uniform\" \n"
          "  varname   : Variable name."),
	OUTPUT(),
	INPUT( batchname_ ),
	GOUTPUT(),
	GINPUT( VECTOR_, VECTOR_ ),
	KEYWORDS( "n",   "period", "amplitude", "pdf",    "varname" ),
	TYPES(    int_t, Numeric_t, Numeric_t,  string_t, string_t  )));

  md_data.push_back
    ( MdRecord
      ( NAME("ybatchAbsAndRte"),
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
          "  do_tags   : This string array gives the tags for which profiles\n"
          "              shall be read from a file, e.g. [\"H2O\",\"O3\"].\n"
          "              These tags must match some tag in tags.\n"
          "  tag_files : Filenames for species data."),
	OUTPUT( ybatch_ ),
	INPUT( // Variables needed for absCalc
               f_mono_, p_abs_, t_abs_, h2o_abs_, vmrs_, lines_per_tg_, 
               lineshape_, 
               // Additional variables for losCalc
	       z_abs_, z_plat_ ,za_pencil_, l_step_, refr_, 
               refr_lfac_, refr_index_, z_ground_, r_geoid_,
               // Additional variables for yRte
	       y_space_, e_ground_, t_ground_,
               // Additional variables needed for this function
               batchname_, tgs_, cont_description_names_, cont_description_parameters_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS("ncalc", "do_t", "t_file", "do_z", "z_file",
                 "do_tags", "tag_files",
                 "do_f", "f_file", "do_za", "za_file"),
	TYPES(   int_t,   int_t,  string_t, int_t,  string_t, 
                 ARRAY_string_t, ARRAY_string_t,
                 int_t,  string_t, int_t,   string_t  )));

  md_data.push_back
    ( MdRecord
      ( NAME("ybatchTB"),
  	DESCRIPTION(
          "Converts a batch of spectra from intensity to brightness\n"
          "temperature. The used frequency vector is f_sensor."),
	OUTPUT( ybatch_ ),
	INPUT( ybatch_, f_sensor_, za_sensor_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("ybatchTRJ"),
  	DESCRIPTION(
          "Converts a batch of spectra from intensity to Rayleigh-Jeans\n"
          "temperature. The used frequency vetor is f_sensor."),
	OUTPUT( ybatch_ ),
	INPUT( ybatch_, f_sensor_, za_sensor_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("ybatchLoadCalibration"),
  	DESCRIPTION(
          "Applies load switching calibration on a batch of spectra.\n"
          "See further yLoadCalibration."),
	OUTPUT( ybatch_ ),
	INPUT( ybatch_, i_cal1_, i_cal2_, y_cal1_, y_cal2_, za_sensor_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("ybatchAdd"),
  	DESCRIPTION(
          "Loads batch data from a file and adds the data to ybatch.\n"
          "The data to add are read from the file\n"
          "   batchname.varname.ab  \n"
          "where varname is the given keyword string.\n"
          "The file can contain data for more spectra than the number of\n"
          "spectra in ybatch (n). In such cases, the only first n columns of\n"
          "the file data is considered.\n"
          "\n"
          "Keywords \n"
          "  n       : Number of random vectors to produce.\n"
          "  varname : Variable name."),
	OUTPUT( ybatch_ ),
	INPUT( ybatch_, batchname_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "varname" ),
	TYPES(    string_t  )));

  md_data.push_back
    ( MdRecord
      ( NAME("ybatchAddScaled"),
  	DESCRIPTION(
          "Loads batch data from a file and adds the data to ybatch.\n"
          "The file data can be scaled by scalefac. By setting the scalfac\n"
          "to -1.0, for example, thermal noise can be removed from a batch\n"
          "of spectra.\n"
          "See also ybatchAdd.\n"
          "\n"
          "Keywords \n"
          "  varname  : Variable name.\n"
          "  scalefac : Scale factor for data. Can be negative."),
	OUTPUT( ybatch_ ),
	INPUT( ybatch_, batchname_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "varname", "scalefac" ),
	TYPES(    string_t,  Numeric_t)));


}
