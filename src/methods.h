/* Copyright (C) 2000-2007 Stefan Buehler <sbuehler@ltu.se>

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

#include <iosfwd>
#include "token.h"
#include "messages.h"
#include "make_array.h"

//! All information for one workspace method.
class MdRecord {
public:

  //! Default constructor.
  MdRecord():
    mname(        ""              ),
    mdescription( ""              ),
    mauthors(     0               ),
    moutput(      0               ),  
    minput(       0               ),   
    mgoutput(     0               ),  
    mginput(      0               ),   
    mkeywords(    0               ),
    mtypes(       0               ),
    magenda_method(false),
    msupergeneric(false),
    msuppress_header(false)
  {};

  // Initializing constructor. Implementation in methods_aux.cc.
  MdRecord(const char                   name[],
           const char                   description[],
           const MakeArray<String>&     authors,
           const MakeArray<Index>&      output,
           const MakeArray<Index>&      input,   
           const MakeArray<Index>&      goutput,
           const MakeArray<Index>&      ginput,   
           const MakeArray<String>&     keywords,
           const MakeArray<TokValType>& types,
           bool                         agenda_method   = false,
           bool                         suppress_header = false
           );

  // Methods returning the lookup information:
  const String&            Name()           const { return mname; }   
  const String&            Description()    const { return mdescription; }
  const ArrayOfString&     Authors()        const { return mauthors; }   
  const ArrayOfIndex&      Output()         const { return moutput; }
  const ArrayOfIndex&      Input()          const { return minput; }
  const ArrayOfIndex&      GOutput()        const { return mgoutput; }
  const ArrayOfIndex&      GInput()         const { return mginput; }
  const Array<String>&     Keywords()       const { return mkeywords; }
  const Array<TokValType>& Types()          const { return mtypes; }
  bool                     AgendaMethod()   const { return magenda_method; }
  bool                     Supergeneric()   const { return msupergeneric; }
  bool                     SuppressHeader() const { return msuppress_header; }
  Index                    ActualGroup()    const { return mactual_group; }

  // Expand supergeneric method record to an actual group
  // (documentation with implementation in method_aux.cc):
  void subst_any_with_group( Index g );

  //! Print method template for the control file. 
  /*!
    This prints the method data exactly in the same way how it can
    be included in the control file. The description string is also
    printed as a comment, but this can be turned off by setting
    show_comment to false.

    @param os Output stream
    @param show_description Should the description string also be printed?
  */
  ostream& PrintTemplate(ostream& os, bool show_description=true) const;

  //! To override the default assignment operator.
  /*! MdRecords cannot be assigned! */
  MdRecord operator=(const MdRecord& m){
    out0 << "MdRecord cannot be assigned!\n"
         << "You tried to assign: " << m << "\n";
    arts_exit ();
    return MdRecord();
  }

  // Needed by make_auto_md_h.cc. See documentation there.
  friend void subst_any_with_group( MdRecord& mdd, Index g );

private:

  //! The name of this method.
  String mname;

  //! A text string describing this method.
  String mdescription;

  //! Author(s) of this method.
  ArrayOfString mauthors;

  //! Workspace Output.
  ArrayOfIndex moutput;

  //! Workspace Input.
  ArrayOfIndex minput;

  //! Generic Workspace Output.
  ArrayOfIndex mgoutput;

  //! Generic Workspace Input.
  ArrayOfIndex mginput;

  //! Keywords.
  ArrayOfString mkeywords;

  //! Types associated with keywords.
  Array<TokValType> mtypes;

  //! Flag, whether this is an agenda method. 
  /*!
    Agenda methods do not expect keywords, but other method
    definitions inside their body in the controlfile.
  */
  bool magenda_method;

  //! Flag, whether this method is supergeneric.
  /*!
    This flag is set automatically if the goutput or ginput contains
    Any_.
  */ 
  bool msupergeneric;

  //! Flag, whether method header should be suppressed.
  /*!
    If we want to implement a supergeneric method by a template
    function, we must not include method headers in auto_md.h.
  */ 
  bool msuppress_header;

  //! The actual group of a supergeneric method.
  /*! 
    This holds the actual group after expansion of a supergeneric method.
  */
  Index mactual_group;
};

void define_md_data_raw();

void expand_md_data_raw_to_md_data();

void define_md_map();

void define_md_raw_map();

ostream& operator<<(ostream& os, const MdRecord& mdr);

#endif  // methods_h
