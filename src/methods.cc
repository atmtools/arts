/* Implementation of MdRecord and of the compound to hold the
   records. See methods.h for more documentation. */

#include "arts.h"
#include "make_array.h"
#include "workspace.h"
#include "methods.h"

// Some #defines and typedefs to make the records better readable:
#define NAME(x) x
#define DESCRIPTION(x) x
#define OUTPUT   make_array<size_t> 
#define INPUT    make_array<size_t> 
#define GOUTPUT  make_array<size_t> 
#define GINPUT   make_array<size_t> 
#define KEYWORDS make_array<string>
#define TYPES    make_array<TokValType>

/** The lookup information for the workspace methods. */
ARRAY<MdRecord> md_data;

/** Initializes the method lookup data. */
void define_md_data()
{
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

  md_data.push_back
    ( MdRecord
      ( NAME("AllAbsExample"),
	DESCRIPTION("Reads all important absorption related variables from the\n"
		    "given files."),
	OUTPUT(f_abs_, p_abs_, t_abs_, abs_),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorWriteToFile"),
	DESCRIPTION("Writes a workspace variable that is a vector to a file.\n"
		    "The filename is <basename>.<variable_name>.a.\n"
		    "The format is exactly the same as for matrices, with the\n"
		    "second dimension given as 1. Explicitly:\n\n"
		    "# <comments>\n\n"
		    //		    "\"<variable_name>\"\n\n"
		    "<n_elements> 1\n\n"
		    "<elements>"),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT(VECTOR_),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixWriteToFile"),
	DESCRIPTION("Writes a workspace variable that is a matrix to a file.\n"
		    "The filename is <basename>.<variable_name>.a.\n"
		    "The format is as follows:\n\n"
		    "# <comments>\n\n"
		    //		    "\"<variable_name>\"\n\n"
		    "<n_rows> <n_columns>\n\n"
		    "<elements>"),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT(MATRIX_),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorWriteToNamedFile"),
	DESCRIPTION("Writes a workspace variable that is a vector to a file.\n"
		    "The filename has to be specified.\n"
		    "The format is exactly the same as for matrices, with the\n"
		    "second dimension given as 1. Explicitly:\n\n"
		    "# <comments>\n\n"
		    //		    "\"<variable_name>\"\n\n"
		    "<n_elements> 1\n\n"
		    "<elements>"),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT(VECTOR_),
	KEYWORDS("filename"),
	TYPES(str_)));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixWriteToNamedFile"),
	DESCRIPTION("Writes a workspace variable that is a matrix to a file.\n"
		    "The filename has to be specified.\n"
		    "The format is as follows:\n\n"
		    "# <comments>\n\n"
		    //		    "\"<variable_name>\"\n\n"
		    "<n_rows> <n_columns>\n\n"
		    "<elements>"),
	OUTPUT(),
	INPUT(),
	GOUTPUT(),
	GINPUT(MATRIX_),
	KEYWORDS("filename"),
	TYPES(str_)));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorReadFromFile"),
	DESCRIPTION("Reads a workspace variable that is a vector from a file.\n"
		    "The filename is <basename>.<variable_name>.a.\n"
		    "The format is exactly the same as for matrices, with the\n"
		    "second dimension given as 1. Explicitly:\n\n"
		    "# <comments>\n\n"
		    //		    "\"<variable_name>\"\n\n"
		    "<n_elements> 1\n\n"
		    "<elements>"),
	OUTPUT(),
	INPUT(),
	GOUTPUT(VECTOR_),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixReadFromFile"),
	DESCRIPTION("Reads a workspace variable that is a matrix from a file.\n"
		    "The filename is <basename>.<variable_name>.a.\n"
		    "The format is as follows:\n\n"
		    "# <comments>\n\n"
		    //		    "\"<variable_name>\"\n\n"
		    "<n_rows> <n_columns>\n\n"
		    "<elements>"),
	OUTPUT(),
	INPUT(),
	GOUTPUT(MATRIX_),
	GINPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("NumericSet"),
	DESCRIPTION("Sets a workspace variable of type Numeric to a value."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(Numeric_),
	GINPUT(),
	KEYWORDS("value"),
	TYPES(num_)));

  md_data.push_back
    ( MdRecord
      ( NAME("losTest"),
	DESCRIPTION("Just to see if los works."),
	OUTPUT(los_),
	INPUT(),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));
  

}
