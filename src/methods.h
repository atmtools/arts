/* Copyright (C) 2000-2008 Stefan Buehler <sbuehler@ltu.se>

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
#include "messages.h"
#include "make_array.h"

#define NODEF "@@THIS_KEYWORD_HAS_NO_DEFAULT_VALUE@@"

//! All information for one workspace method.
class MdRecord {
public:

  //! Default constructor.
  MdRecord():
    mname(        ""              ),
    mdescription( ""              ),
    mauthors(     0               ),
    moutput(      0               ),  
    mgout(        0               ),  
    mgouttype(    0               ),  
    minput(       0               ),   
    mgin(         0               ),   
    mgintype(     0               ),   
    mgindefault(  0               ),
    mset_method(false),
    magenda_method(false),
    msupergeneric(false),
    msuppress_header(false),
    mpass_workspace(false),
    mpass_wsv_names(false),
    mactual_group(-1)
  {};

  // Initializing constructor. Implementation in methods_aux.cc.
  MdRecord(const char                   name[],
           const char                   description[],
           const MakeArray<String>&     authors,
           const MakeArray<String>&     output,
           const MakeArray<String>&     gout,
           const MakeArray<String>&     gouttype,
           const MakeArray<String>&     input,   
           const MakeArray<String>&     gin,
           const MakeArray<String>&     gintype,   
           const MakeArray<String>&     gindefault _U_,
           bool                         set_method      = false,
           bool                         agenda_method   = false,
           bool                         suppress_header = false,
           bool                         pass_workspace  = false,
           bool                         pass_wsv_names  = false
           );

  // Methods returning the lookup information:
  const String&            Name()           const { return mname; }   
  const String&            Description()    const { return mdescription; }
  const ArrayOfString&     Authors()        const { return mauthors; }   
  const ArrayOfIndex&      Out()            const { return moutput; }
  const ArrayOfString&     GOut()           const { return mgout; }
  const ArrayOfIndex&      GOutType()       const { return mgouttype; }
  const ArrayOfIndex&      In()             const { return minput; }
  const ArrayOfString&     GIn()            const { return mgin; }
  const ArrayOfIndex&      GInType()        const { return mgintype; }
  const Array<String>&     GInDefault()     const { return mgindefault; }
  bool                     SetMethod()      const { return mset_method; }
  bool                     AgendaMethod()   const { return magenda_method; }
  bool                     Supergeneric()   const { return msupergeneric; }
  bool                     SuppressHeader() const { return msuppress_header; }
  bool                     PassWorkspace()  const { return mpass_workspace; }
  bool                     PassWsvNames()   const { return mpass_wsv_names; }
  Index                    ActualGroup()    const { return mactual_group; }
  void                     SetPassWorkspace()  { mpass_workspace = true; }

  // Expand supergeneric method record to an actual group
  // (documentation with implementation in method_aux.cc):
  void subst_any_with_group( Index g );

  // Helper function returning a list of WSVs which are only input, not
  // output. (documentation with implementation in method_aux.cc):
  void input_only(ArrayOfIndex& inonly) const;

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
  MdRecord& operator=(const MdRecord& m){
    out0 << "MdRecord cannot be assigned!\n"
         << "You tried to assign: " << m << "\n";
    arts_exit ();
    return *this;
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

  //! Generic Workspace Output Names.
  ArrayOfString mgout;

  //! Generic Workspace Output Type.
  ArrayOfIndex mgouttype;

  //! Workspace Input.
  ArrayOfIndex minput;

  //! Generic Workspace Input Names.
  ArrayOfString mgin;

  //! Generic Workspace Input.
  ArrayOfIndex mgintype;

  //! Generic Workspace Input Defaults.
  ArrayOfString mgindefault;

  //! Flag, whether this is a set method. 
  /*!
    Set methods have exactly one generic input. Unlike other inputs these
    are not stored in a workspace variable.
  */
  bool mset_method;

  //! Flag, whether this is an agenda method. 
  /*!
    Agenda methods expect other method definitions inside their body
    in the controlfile.
  */
  bool magenda_method;

  //! Flag, whether this method is supergeneric.
  /*!
    This flag is set automatically if the gouttype or gintype contains
    Any_.
  */ 
  bool msupergeneric;

  //! Flag, whether method header should be suppressed.
  /*!
    If we want to implement a supergeneric method by a template
    function, we must not include method headers in auto_md.h.
  */ 
  bool msuppress_header;

  //! Flag, whether a workspace reference should be passed to the WSM.
  /*!
    Some WSMs (like Delete) need direct access to the workspace object. If
    this flag is set to true, the gateway function will take care of that.
  */ 
  bool mpass_workspace;

  //! Flag, whether WSV names should be passed to the WSM.
  /*!
    Some WSMs (like ReadXML, WriteXML) need to know the names of the WSVs that
    have been passed to them. If this flag is set to true, the gateway function
    will take care of that.
  */ 
  bool mpass_wsv_names;

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
