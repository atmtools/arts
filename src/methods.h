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
	   bool  generic,
	   const ARRAY<size_t>&  output,
	   const ARRAY<size_t>&  input,   
	   const ARRAY<string>&     keywords,
	   const ARRAY<TokValType>& types) :
    mname(name),
    mdescription(description),
    mgeneric(generic),
    moutput(output),  
    minput(input),   
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
  bool                     Generic()      const { return mgeneric;     }
  const ARRAY<size_t>&     Output()       const { return moutput;      }
  const ARRAY<size_t>&     Input()        const { return minput;       }
  const ARRAY<string>&     Keywords()     const { return mkeywords;    }
  const ARRAY<TokValType>& Types()        const { return mtypes;       }

private:

  /** The name of this method. */
  string mname;

  /** A text string describing this method. */
  string mdescription;

  /** Is this a generic method? The workspace output and input is
      speciefied either as lists of individual workspace variables
      (identified by their handle), or as lists of groups. It is not
      possible to merge these to concepts!  Generic==true indicates
      that this is a generic method, hence groups are
      used. Generic==false means that this is a normal method with
      explicit WsvHandles. */
  bool mgeneric;
  
  /** Workspace Output. */
  ARRAY<size_t> moutput;

  /** Workspace Input. */
  ARRAY<size_t> minput;

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
