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
         "\n"
         "Paragraphs are seperated with blank lines.\n"
         "\n"
         "The names of workspace variables and other methods are marked by\n"
         "stars, for example *z_plat*.\n"
         "\n"
         "Global input and output, and keywords shall be described \n"
         "as exemplified below. If there is no variables of a group, \n"
         "(e.g. global input) remove that part totally. Note that the \n"
         "on-line help just gives the type of global input/output and the \n"
         "keyword names, and additional information is for sure needed.\n"
         "\n"
         "Leave space and brake lines when listing input and output \n"
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
         "Global input: \n"
         "   \n"
         "\n"
         "Global output: \n"
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




void define_md_data()
{
  // The variable md_data is defined in file methods_aux.cc.
  extern Array<MdRecord> md_data;

  // Initialize to zero, just in case:
  md_data.resize(0);


  /////////////////////////////////////////////////////////////////////////////
  // Let's put in the functions in alphabetical order. This gives a clear rule
  // for where to place a new function and this gives a nicer results when
  // the functions are listed by "arts -m all".
  // No distinction is made between uppercase and lowercase letters. The sign
  // "_" comes after all letters.
  // Patrick Eriksson 2002-05-08
  /////////////////////////////////////////////////////////////////////////////


  md_data.push_back
    ( MdRecord
      ( NAME("AgendaDefine"),
  	DESCRIPTION
	( 
         "Set up an agenda.\n"
	 "\n"
	 "A method list just contains indices (in md_data) of methods\n"
	 "intended for sequential execution. Only methods without keyword\n"
	 "arguments are allowed. It is the task of this method to set this\n"
	 "up. For example, it must be checked, whether the given names\n"
	 "really correspond to methods.\n"
	 "\n"
	 "Global output:\n"
	 "   ArrayOfIndex : The newly generated method list.\n"
	 "\n"
	 "Keywords:\n"
	 "   methods      : An array of names of methods." 
	),
	OUTPUT(  ),
        INPUT(  ),
	GOUTPUT( Agenda_ ),
	GINPUT(),
	KEYWORDS(),
	TYPES(),
	AGENDAMETHOD( true )));

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexReadXML"),
	DESCRIPTION(
                    "Reads a index array from an XML file.\n"
                    "\n"
                    "The index array is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the index array is read\n"
                    "from <basename>.<variable_name>.xml.\n"
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

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexWriteXML"),
	DESCRIPTION(
                    "Writes a index array to an XML file.\n"
                    "\n"
                    "The index array of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the index array is written\n"
                    "to <basename>.<variable_name>.xml.\n"
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
      ( NAME("ArrayOfMatrixReadXML"),
	DESCRIPTION(
                    "Reads an array of matrices from an XML file.\n"
                    "\n"
                    "The array of matrices is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the array of matrices is read\n"
                    "from <basename>.<variable_name>.xml.\n"
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

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixWriteXML"),
	DESCRIPTION(
                    "Writes an array of matrices to an XML file.\n"
                    "\n"
                    "The array of matrices of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the array of matrices is written\n"
                    "to <basename>.<variable_name>.xml.\n"
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
      ( NAME("ArrayOfStringReadXML"),
	DESCRIPTION(
                    "Reads an array of strings from an XML file.\n"
                    "\n"
                    "The array of strings is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the string is read\n"
                    "from <basename>.<variable_name>.xml.\n"
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
      ( NAME("ArrayOfStringWriteXML"),
	DESCRIPTION(
                    "Writes an array of strings to an XML file.\n"
                    "\n"
                    "The array of strings of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the string is written\n"
                    "to <basename>.<variable_name>.xml.\n"
                    "\n"
                    "The format is as follows:\n"
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
      ( NAME("ArrayOfVectorReadXML"),
	DESCRIPTION(
                    "Reads an array of vectors from an XML file.\n"
                    "\n"
                    "The array of vectors is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the array of vectors is read\n"
                    "from <basename>.<variable_name>.xml.\n"
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

  md_data.push_back
    ( MdRecord
      ( NAME("ArrayOfVectorWriteXML"),
	DESCRIPTION(
                    "Writes an array of vectors to an XML file.\n"
                    "\n"
                    "The array of vectors of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the array of vectors is written\n"
                    "to <basename>.<variable_name>.xml.\n"
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

  md_data.push_back
    ( MdRecord
      ( NAME("IndexReadXML"),
	DESCRIPTION(
                    "Reads a index value from an XML file.\n"
                    "\n"
                    "The index value is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the index is read\n"
                    "from <basename>.<variable_name>.xml.\n"
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

  md_data.push_back
    ( MdRecord
      ( NAME("IndexWriteXML"),
	DESCRIPTION(
                    "Writes an index value to an XML file.\n"
                    "\n"
                    "The index value of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the index is written\n"
                    "to <basename>.<variable_name>.xml.\n"
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
      ( NAME("IndexSet"),
	DESCRIPTION
        (
         "Sets an index workspace variable to the given value. \n"
         "\n"
         "Global output: \n"
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

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixCopy"),
	DESCRIPTION
        (
         "Creates a copy of a matrix. \n"
         "\n"
         "Global output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Global input: \n"
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
         "Global output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Global input: \n"
         "   Vector : The vector to be copied. \n"
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

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixReadXML"),
	DESCRIPTION(
                    "Reads a matrix from an XML file.\n"
                    "\n"
                    "The matrix is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the matrix is read\n"
                    "from <basename>.<variable_name>.xml.\n"
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

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixScale"),
	DESCRIPTION
        (
         "Scales all elements of a matrix with the same value. \n"
         "\n"
         "The result can either be stored in the same or another matrix. \n"
         "\n"
         "Global output: \n"
         "   Matrix : Return matrix. \n"
         "\n"
         "Global input: \n"
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
      ( NAME("MatrixSet"),
	DESCRIPTION
        (
         "Creates a workspace matrix with the specified size and sets \n"
         "all values of the matrix to the specified value. \n"
         "\n"
         "Global output: \n"
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
      ( NAME("MatrixWriteXML"),
	DESCRIPTION(
                    "Writes a matrix to an XML file.\n"
                    "\n"
                    "The matrix of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the matrix is written\n"
                    "to <basename>.<variable_name>.xml.\n"
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
      ( NAME("NoCloudbox"),
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

  md_data.push_back
    ( MdRecord
      ( NAME("NumericReadXML"),
	DESCRIPTION(
                    "Reads a numeric value from an XML file.\n"
                    "\n"
                    "The numeric value is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the numeric is read\n"
                    "from <basename>.<variable_name>.xml.\n"
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

  md_data.push_back
    ( MdRecord
      ( NAME("NumericSet"),
	DESCRIPTION
        (
         "Sets a numeric workspace variable to the given value. \n"
         "\n"
         "Global output: \n"
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
      ( NAME("NumericWriteXML"),
	DESCRIPTION(
                    "Writes a numeric value to an XML file.\n"
                    "\n"
                    "The numeric value of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the numeric is written\n"
                    "to <basename>.<variable_name>.xml.\n"
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
      ( NAME( "ppathCalc" ),
	DESCRIPTION
        (
	 "To be written. The function is being implemented.\n"
         ""
        ),
	OUTPUT( ppath_ ),
	INPUT( atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_, 
               r_geoid_, z_ground_, blackbody_ground_, 
               cloudbox_on_, cloudbox_limits_, a_pos_, a_los_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("SetAtmosphericDimensionality"),
	DESCRIPTION
        (
         "Sets *atmosphere_dim* and gives some variables dummy values.\n"
         "\n"
         "The function sets *atmosphere_dim* to the selected value. The\n"
         "latitude and longitude grids are set to be empty."
        ),
	OUTPUT( atmosphere_dim_, lat_grid_, lon_grid_ ),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS( "dim" ),
	TYPES(    Index_t )));

  md_data.push_back
    ( MdRecord
      ( NAME("StringReadXML"),
	DESCRIPTION(
                    "Reads a string from an XML file.\n"
                    "\n"
                    "The string is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the string is read\n"
                    "from <basename>.<variable_name>.xml.\n"
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
      ( NAME("StringWriteXML"),
	DESCRIPTION(
                    "Writes a string to an XML file.\n"
                    "\n"
                    "The string of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the string is written\n"
                    "to <basename>.<variable_name>.xml.\n"
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
         "Global output: \n"
         "   Tensor3 : The tensor to be created. \n"
         "\n"
         "Global input: \n"
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

  md_data.push_back
    ( MdRecord
      ( NAME("VectorAddScalar"),
  	DESCRIPTION
        (
         "Adds a scalar to all elements of a vector. \n"
         "\n"
         "The result can either be stored in the same or another vector. \n"
         "\n"
         "Global output: \n"
         "   Vector : Return vector. \n"
         "\n"
         "Global input: \n"
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
      ( NAME("VectorCopy"),
	DESCRIPTION
        (
         "Creates a copy of a vector. \n"
         "\n"
         "Global output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Global input: \n"
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
         "Global output: \n"
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
         "Global output: \n"
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
         "Global output: \n"
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
      ( NAME("VectorReadXML"),
	DESCRIPTION(
                    "Reads a vector from an XML file.\n"
                    "\n"
                    "The vector is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the vector is read\n"
                    "from <basename>.<variable_name>.xml.\n"
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

  md_data.push_back
    ( MdRecord
      ( NAME("VectorScale"),
  	DESCRIPTION
        (
         "Scales all elements of a vector with the same value. \n"
         "\n"
         "The result can either be stored in the same or another vector. \n"
         "\n"
         "Global output: \n"
         "   Vector : Return vector. \n"
         "\n"
         "Global input: \n"
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

  md_data.push_back
    ( MdRecord
      ( NAME("VectorSet"),
	DESCRIPTION
        (
         "Creates a workspace vector with the specified length and sets \n"
         "all values of the vector to the specified value. \n"
         "\n"
         "Global output: \n"
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
      ( NAME("VectorSetTakingLengthFromVector"),
	DESCRIPTION
        (
         "Creates a workspace vector with the same length as another vector,\n"
         "and sets all values of the new vector to the specified value. \n"
         "\n"
         "A common usage of the function should be: \n"
         "  VectorSetLengthFromVector(e_ground,f_mono){value=0.75} \n"
         "\n"
         "Global output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Global input: \n"
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
      ( NAME("VectorWriteXML"),
	DESCRIPTION(
                    "Writes a vector to an XML file.\n"
                    "\n"
                    "The vector of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the vector is written\n"
                    "to <basename>.<variable_name>.xml.\n"
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
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
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
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
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
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
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




//=== Matrix ==========================================================


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
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
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
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
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
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
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
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
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
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
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
                    "See *ArrayOfStringWriteAscii* for file format.\n"
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
                    "See *ArrayOfStringWriteAscii* for file format.\n"
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
          "Sets the lineshape for all calculated lines.\n\n"
          "\n"
          "   A general lineshape profile is specified, according to a given  \n"
          "approximation. Alongside a normalization factor is to be set - a  \n"
          "multiplicative forefactor through which the profile can be \n"
          "modified. This factor is just the 0th or 1st, or 2nd power of the \n"
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
          "normalizationfactor - \"no_norm\": 1\n"
          "                      \"linear\": f/f0\n" 
          "                      \"quadratic\": (f/f0)^2.\n"
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



}

