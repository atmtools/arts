/* Copyright (C) 2000, 2001 Stefan Buehler <sbuehler@uni-bremen.de>

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

#include "token.h"
#include "make_array.h"

/** This class contains all information for one workspace method. */
class MdRecord {
public:

  /** Default constructor. */
  MdRecord():
      mname(      ""              ),
    mdescription( ""              ),
    moutput(      0               ),  
    minput(       0               ),   
    mgoutput(     0               ),  
    mginput(      0               ),   
    mkeywords(    0               ),
    mtypes(       0               )
  {};

  /** The only non-trivial constructor, which sets all the
      fields. */
  MdRecord(const char                   name[],
           const char                   description[],
           const MakeArray<Index>&      output,
           const MakeArray<Index>&      input,   
           const MakeArray<Index>&      goutput,
           const MakeArray<Index>&      ginput,   
           const MakeArray<String>&     keywords,
           const MakeArray<TokValType>& types) :
    mname(        name            ),
    mdescription( description     ),
    moutput(      output          ),  
    minput(       input           ),   
    mgoutput(     goutput         ),  
    mginput(      ginput          ),   
    mkeywords(    keywords        ),
    mtypes(       types           )
    { 
      // Initializing the various arrays with input data should now
      // work correctly.  

      // Keywords and type must have the same number of
      // elements. (Types specifies the types associated with each
      // keyword.)
      assert( mkeywords.nelem() == mtypes.nelem() );
    }
  
  const String&            Name()         const { return mname;        }   
  const String&            Description()  const { return mdescription; }
  const ArrayOfIndex&      Output()       const { return moutput;      }
  const ArrayOfIndex&      Input()        const { return minput;       }
  const ArrayOfIndex&      GOutput()      const { return mgoutput;      }
  const ArrayOfIndex&      GInput()       const { return mginput;       }
  const Array<String>&     Keywords()     const { return mkeywords;    }
  const Array<TokValType>& Types()        const { return mtypes;       }

  /** Print method template for the control file. This prints the
      method data exactly in the same way how it can be included in
      the control file. The description string is also printed as a
      comment, but this can be turned off by setting show_comment to
      false.

      @param show_description Should the description string also be printed?   */
  ostream& PrintTemplate(ostream& os, bool show_description=true) const;

  /** To override the default assignment operator. MdRecords cannot be
      assigned! */
  MdRecord operator=(const MdRecord&){
    cout << "MdRecord cannot be assigned!\n";
    exit(1);
      }
private:

  /** The name of this method. */
  String mname;

  /** A text string describing this method. */
  String mdescription;

  /** Workspace Output. */
  ArrayOfIndex moutput;

  /** Workspace Input. */
  ArrayOfIndex minput;

  /** Generic Workspace Output. */
  ArrayOfIndex mgoutput;

  /** Generic Workspace Input. */
  ArrayOfIndex mginput;

  /** Keywords. */
  ArrayOfString mkeywords;

  /** Types associated with keywords. */
  Array<TokValType> mtypes;

};


// Some #defines and typedefs to make the records better readable:
#define NAME(x) x 
#define DESCRIPTION(x) x
#define OUTPUT   MakeArray<Index>
#define INPUT    MakeArray<Index>
#define GOUTPUT  MakeArray<Index>
#define GINPUT   MakeArray<Index>
#define KEYWORDS MakeArray<String>
#define TYPES    MakeArray<TokValType>

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
