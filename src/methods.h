/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>

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
  \file   methods.h
  \brief  Declaration of the class MdRecord.

  The class MdRecord contains all information for one workspace method. 

  \author Stefan Buehler
  \date   1999-07-27
*/

#ifndef methods_h
#define methods_h

#include "vecmat.h"
#include "token.h"
#include "auto_wsv.h"

/** This class contains all information for one workspace method. */
class MdRecord {
public:

  /** Default constructor. */
  MdRecord(){};

  /** The only one non-trivial constructor, which sets all the
      fields. */
  MdRecord(const char 		    name[],
	   const char 		    description[],
	   const Array<size_t>&     output,
	   const Array<size_t>&     input,   
	   const Array<size_t>&     goutput,
	   const Array<size_t>&     ginput,   
	   const Array<string>&     keywords,
	   const Array<TokValType>& types) :
    mname(        name            ),
    mdescription( description     ),
    moutput(      output.size()   ),  
    minput(       input.size()    ),   
    mgoutput(     goutput.size()  ),  
    mginput(      ginput.size()   ),   
    mkeywords(    keywords.size() ),
    mtypes(       types.size()    )
    { 
      // We need to use copy to initialize the Array members. If we use
      // the assignment operator they end up all pointing to the same
      // data!
      copy( output   , moutput   );
      copy( input    , minput    );
      copy( goutput  , mgoutput  );
      copy( ginput   , mginput   );
      copy( keywords , mkeywords );
      copy( types    , mtypes    );

      // Keywords and type must have the same number of
      // elements. (Types specifies the types associated with each
      // keyword.)
      assert( mkeywords.size() == mtypes.size() );
    }
  
  const string&            Name()         const { return mname;        }   
  const string&            Description()  const { return mdescription; }
  const Array<size_t>&     Output()       const { return moutput;      }
  const Array<size_t>&     Input()        const { return minput;       }
  const Array<size_t>&     GOutput()      const { return mgoutput;      }
  const Array<size_t>&     GInput()       const { return mginput;       }
  const Array<string>&     Keywords()     const { return mkeywords;    }
  const Array<TokValType>& Types()        const { return mtypes;       }

  /** Print method template for the control file. This prints the
      method data exactly in the same way how it can be included in
      the control file. The description string is also printed as a
      comment, but this can be turned off by setting show_comment to
      false.

      @param show_description Should the description string also be printed?   */
  ostream& PrintTemplate(ostream& os, bool show_description=true) const;
private:

  /** The name of this method. */
  string mname;

  /** A text string describing this method. */
  string mdescription;

  /** Workspace Output. */
  Array<size_t> moutput;

  /** Workspace Input. */
  Array<size_t> minput;

  /** Generic Workspace Output. */
  Array<size_t> mgoutput;

  /** Generic Workspace Input. */
  Array<size_t> mginput;

  /** Keywords. */
  Array<string> mkeywords;

  /** Types associated with keywords. */
  Array<TokValType> mtypes;

};


// Some #defines and typedefs to make the records better readable:
#define NAME(x) x
#define DESCRIPTION(x) x
#define OUTPUT   make_array<size_t> 
#define INPUT    make_array<size_t> 
#define GOUTPUT  make_array<size_t> 
#define GINPUT   make_array<size_t> 
#define KEYWORDS make_array<string>
#define TYPES    make_array<TokValType>


/** Define the lookup data for the workspace methods. The array
    md_data contains all that we need to know about each method. The
    lookup table is a global variable. It can be made visible anywhere
    with an extern declaration. */
void define_md_data();

/** Define MdMap. MdMap can be used to find method data by method
    name. */
void define_md_map();

/** Output operator for MdRecord.
    @author Stefan Buehler */
ostream& operator<<(ostream& os, const MdRecord& mdr);


#endif  // methods_h
