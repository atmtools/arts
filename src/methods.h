/* This file contains the declaration of the class MdRecord, which
   contains all information for one workspace method.

   History:
   SAB 27.07.1999 Started.
*/

#ifndef methods_h
#define methods_h

#include "vecmat.h"
#include "token.h"
#include "wsv.h"

/** This class contains all information for one workspace method. */
class MdRecord {
public:
  /** The only one non-trivial constructor, which sets all the
      fields. */
  MdRecord(const char name[],
	   const char description[],
	   const ARRAY<size_t>&  output,
	   const ARRAY<size_t>&  input,   
	   const ARRAY<size_t>&  goutput,
	   const ARRAY<size_t>&  ginput,   
	   const ARRAY<string>&     keywords,
	   const ARRAY<TokValType>& types) :
    mname(name),
    mdescription(description),
    moutput(output),  
    minput(input),   
    mgoutput(goutput),  
    mginput(ginput),   
    mkeywords(keywords),
    mtypes(types)
    { 
      // Keywords and type must have the same number of
      // elements. (Types specifies the types associated with each
      // keyword.)
      assert( mkeywords.size() == mtypes.size() );
    }
  
  const string&            Name()         const { return mname;        }   
  const string&            Description()  const { return mdescription; }
  const ARRAY<size_t>&     Output()       const { return moutput;      }
  const ARRAY<size_t>&     Input()        const { return minput;       }
  const ARRAY<size_t>&     GOutput()      const { return mgoutput;      }
  const ARRAY<size_t>&     GInput()       const { return mginput;       }
  const ARRAY<string>&     Keywords()     const { return mkeywords;    }
  const ARRAY<TokValType>& Types()        const { return mtypes;       }

private:

  /** The name of this method. */
  string mname;

  /** A text string describing this method. */
  string mdescription;

  /** Workspace Output. */
  ARRAY<size_t> moutput;

  /** Workspace Input. */
  ARRAY<size_t> minput;

  /** Generic Workspace Output. */
  ARRAY<size_t> mgoutput;

  /** Generic Workspace Input. */
  ARRAY<size_t> mginput;

  /** Keywords. */
  ARRAY<string> mkeywords;

  /** Types associated with keywords. */
  ARRAY<TokValType> mtypes;

};

/** Define the lookup data for the workspace methods. The array
    md_data contains all that we need to know about each method. The
    lookup table is a global variable. It can be made visible anywhere
    with an extern declaration. */
void define_md_data();


#endif  // methods_h
