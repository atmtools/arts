/* Implementation of MdRecord and of the compound to hold the
   records. See methods.h for more documentation. */

#include "arts.h"
#include "vecmat.h"
#include "token.h"
#include "wsv.h"
#include "make_vector.h"
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
	DESCRIPTION("Reads all important absorption related variables from the\n."
		    "given files."),
	GENERIC(false),
	OUTPUT(p_abs_, f_abs_, t_abs_, abs_),
	INPUT(),
	KEYWORDS(),
	TYPES()));

  
// MdRecord md_data[] = {
//   MdRecord,
//   MdRecord("SomeMethod",
// 	   "Do whatever it takes to fulfil the expectations that have\n"
// 	   "been put in you.",
// 	   false,
// 	   make_array<size_t>(abs_),
// 	   make_array<size_t>(p_abs_, t_abs_, f_abs_),
// 	   make_array<string>("sample","ints","numm"),
// 	   make_array<TokValType>(strvec_,intvec_,numvec_),
// 	   make_array<TokVal>(make_array<string>("ab","cd"),
// 			      make_array<int>(1,2),
// 			      make_array<Numeric>(1.2,2.3)
// 			     )
//           ),
//   MdRecord("WriteVectorToFile",
// 	   "This writes any vector to an ASCI File.\n"
// 	   "The filename is given by the basename plus the name\n"
// 	   "of the workspace variable as the extension.",
// 	   true,
// 	   make_array<size_t>(),
// 	   make_array<size_t>(NumVector_),
// 	   make_array<string>(),
// 	   make_array<TokValType>(),
// 	   make_array<TokVal>()),
//   MdRecord("End",
// 	   "This method terminates the program.",
// 	   false,
// 	   make_array<size_t>(),
// 	   make_array<size_t>(),
// 	   make_array<string>(),
// 	   make_array<TokValType>(),
// 	   make_array<TokVal>())  
// };

}
