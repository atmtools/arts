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


//======================================================================
//=== IO methods
//======================================================================

//
//=== Scalars
//
  md_data.push_back
    ( MdRecord
      ( NAME("IntSet"),
	DESCRIPTION("Sets an integer workspace variable to the given value."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(int_t),
	GINPUT(),
	KEYWORDS("value"),
	TYPES(int_t)));

  md_data.push_back
    ( MdRecord
      ( NAME("NumericSet"),
	DESCRIPTION("Sets a workspace variable of type Numeric to a value."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(Numeric_),
	GINPUT(),
	KEYWORDS("value"),
	TYPES(Numeric_t)));



//
//=== Vector initialization
//

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
	DESCRIPTION("Creates a linearly spaced vector with defined length.\n"
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
	DESCRIPTION("Creates a logarithmically spaced vector with defined length.\n"
		    "The length must be larger than 1."),
	OUTPUT(),
	INPUT(),
	GOUTPUT(VECTOR_),
	GINPUT(),
	KEYWORDS("start", "stop", "n"),
	TYPES(Numeric_t, Numeric_t, int_t)));



//
//=== Matrix and Vector Write Methods
//
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
	TYPES(string_t)));

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
	TYPES(string_t)));



//
//=== Matrix and Vector Read Methods
//

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



//======================================================================
//=== Absorption methods
//======================================================================
  
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
  


//======================================================================
//=== LOS methods
//======================================================================

  md_data.push_back
    ( MdRecord
      ( NAME("losGeneral"),
  	DESCRIPTION("A general function to determine LOS for a 1D atmosphere\n"
                  "Refraction variables and ground altitude and emission\n"
                  "must be set when using this function."),
	OUTPUT(los_),
	INPUT(z_plat_ ,view1_, l_step_, p_abs_, z_abs_, refr_, l_step_refr_, 
              z_ground_, e_ground_ ),
	GOUTPUT(),
	GINPUT(),
	KEYWORDS(),
	TYPES()));


}
