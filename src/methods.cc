/* Implementation of MdRecord and of the compound to hold the
   records. See methods.h for more documentation. */

#include "arts.h"
#include "make_array.h"
#include "workspace.h"
#include "methods.h"

// Some #defines and typedefs to make the records better readable:
#define NAME(x) x
#define DESCRIPTION(x) x
#define GENERIC(x) x
#define OUTPUT make_array<size_t> 
#define INPUT make_array<size_t> 
#define KEYWORDS make_array<string>
#define TYPES make_array<TokValType>

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
	GENERIC(false),
	OUTPUT(),
	INPUT(),
	KEYWORDS(""),
	TYPES()
	));
  */

  md_data.push_back
    ( MdRecord
      ( NAME("AllAbsExample"),
	DESCRIPTION("Reads all important absorption related variables from the\n"
		    "given files."),
	GENERIC(false),
	OUTPUT(f_abs_, p_abs_, t_abs_, abs_),
	INPUT(),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorWriteToFile"),
	DESCRIPTION("Writes a workspace variable that is a vector to a file.\n"
		    "The filename is <basename>.<variable_name>.\n"
		    "The format is exactly the same as for matrices, with the\n"
		    "second dimension given as zero. Explicitly:\n\n"
		    "# <comments>\n\n"
		    //		    "\"<variable_name>\"\n\n"
		    "<n_elements> 0\n\n"
		    "<elements>"),
	GENERIC(true),
	OUTPUT(),
	INPUT(VECTOR_),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("MatrixWriteToFile"),
	DESCRIPTION("Writes a workspace variable that is a matrix to a file.\n"
		    "The filename is <basename>.<variable_name>.\n"
		    "The format is as follows:\n\n"
		    "# <comments>\n\n"
		    //		    "\"<variable_name>\"\n\n"
		    "<n_rows> <n_columns>\n\n"
		    "<elements>"),
	GENERIC(true),
	OUTPUT(),
	INPUT(MATRIX_),
	KEYWORDS(),
	TYPES()));

  md_data.push_back
    ( MdRecord
      ( NAME("VectorWriteToNamedFile"),
	DESCRIPTION("Writes a workspace variable that is a vector to a file.\n"
		    "The filename has to be specified.\n"
		    "The format is exactly the same as for matrices, with the\n"
		    "second dimension given as zero. Explicitly:\n\n"
		    "# <comments>\n\n"
		    //		    "\"<variable_name>\"\n\n"
		    "<n_elements> 0\n\n"
		    "<elements>"),
	GENERIC(true),
	OUTPUT(),
	INPUT(VECTOR_),
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
	GENERIC(true),
	OUTPUT(),
	INPUT(MATRIX_),
	KEYWORDS("filename"),
	TYPES(str_)));

  

}
