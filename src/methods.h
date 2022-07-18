/* Copyright (C) 2000-2012 Stefan Buehler <sbuehler@ltu.se>

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
#include "matpackI.h"
#include "messages.h"

#define NODEF "@@THIS_KEYWORD_HAS_NO_DEFAULT_VALUE@@"

//! All information for one workspace method.
class MdRecord {
 public:
  //! Default constructor.
  MdRecord()
      : mname(""),
        mdescription(""),
        mauthors(0),
        moutput(0),
        mgout(0),
        mgouttype(0),
        mgoutdesc(0),
        minput(0),
        mgin(0),
        mgintype(0),
        mgindefault(0),
        mgindesc(0),
        mset_method(false),
        magenda_method(false),
        msupergeneric(false),
        muses_templates(false),
        mpass_workspace(true),
        mpass_wsv_names(false),
        mactual_groups("") {}

  // Initializing constructor. Implementation in methods_aux.cc.
  MdRecord(const char* name,
           const char* description,
           const ArrayOfString& authors,
           const ArrayOfString& output,
           const ArrayOfString& gout,
           const ArrayOfString& gouttype,
           const ArrayOfString& goutdesc,
           const ArrayOfString& input,
           const ArrayOfString& gin,
           const ArrayOfString& gintype,
           const ArrayOfString& gindefault,
           const ArrayOfString& gindesc,
           bool set_method = false,
           bool agenda_method = false,
           bool uses_templates = false,
           bool pass_workspace = false,
           bool pass_wsv_names = false);

  MdRecord(const MdRecord&) = default;
  MdRecord(MdRecord&&) = default;

  // Methods returning the lookup information:
  [[nodiscard]] const String& Name() const { return mname; }
  [[nodiscard]] const String& Description() const { return mdescription; }
  [[nodiscard]] const ArrayOfString& Authors() const { return mauthors; }
  [[nodiscard]] const ArrayOfIndex& Out() const { return moutput; }
  [[nodiscard]] const ArrayOfString& GOut() const { return mgout; }
  [[nodiscard]] const ArrayOfIndex& GOutType() const { return mgouttype; }
  [[nodiscard]] const ArrayOfArrayOfIndex& GOutSpecType() const { return mgoutspectype; }
  [[nodiscard]] const Array<String>& GOutDescription() const { return mgoutdesc; }
  [[nodiscard]] const ArrayOfIndex& In() const { return minput; }
  [[nodiscard]] const ArrayOfString& GIn() const { return mgin; }
  [[nodiscard]] const ArrayOfIndex& GInType() const { return mgintype; }
  [[nodiscard]] const ArrayOfArrayOfIndex& GInSpecType() const { return mginspectype; }
  [[nodiscard]] const Array<String>& GInDefault() const { return mgindefault; }
  [[nodiscard]] const Array<String>& GInDescription() const { return mgindesc; }
  [[nodiscard]] const ArrayOfIndex& InOnly() const { return minonly; }
  [[nodiscard]] const ArrayOfIndex& InOut() const { return minout; }
  [[nodiscard]] const ArrayOfIndex& OutOnly() const { return moutonly; }
  [[nodiscard]] bool SetMethod() const { return mset_method; }
  [[nodiscard]] bool AgendaMethod() const { return magenda_method; }
  [[nodiscard]] bool Supergeneric() const { return msupergeneric; }
  [[nodiscard]] bool UsesTemplates() const { return muses_templates; }
  [[nodiscard]] bool PassWorkspace() const { return mpass_workspace; }
  [[nodiscard]] bool PassWsvNames() const { return mpass_wsv_names; }
  [[nodiscard]] const String& ActualGroups() const { return mactual_groups; }
  void SetPassWorkspace() { mpass_workspace = true; }

  // Expand supergeneric method record to an actual group
  // (documentation with implementation in method_aux.cc):
  void subst_any_with_group(Index g);
  void subst_any_with_specific_group(Index g);

  //! Print method template for the control file.
  /*!
    This prints the method data exactly in the same way how it can
    be included in the control file. The description string is also
    printed as a comment, but this can be turned off by setting
    show_comment to false.

    @param os Output stream
    @param show_description Should the description string also be printed?
  */
  ostream& PrintTemplate(ostream& os, bool show_description = true) const;

  //! To override the default assignment operator.
  /*! MdRecords cannot be assigned! */
  MdRecord& operator=(const MdRecord& m) {
    Verbosity verbosity;
    ArtsOut0 out0(verbosity);
    out0 << "MdRecord cannot be assigned!\n"
         << "You tried to assign: " << m << "\n";
    arts_exit();
    return *this;
  }

  // Needed by make_auto_md_h.cc. See documentation there.
  friend void subst_any_with_group(MdRecord& mdd, Index g);

  friend ostream& operator<<(ostream& os, const MdRecord& mdr);

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

  //! Generic Workspace Output Types (Contains the valid types if the method
  // is not completely supergeneric but only defined for these types).
  ArrayOfArrayOfIndex mgoutspectype;

  //! Generic Workspace Output Description.
  ArrayOfString mgoutdesc;

  //! Workspace Input.
  ArrayOfIndex minput;

  //! Generic Workspace Input Names.
  ArrayOfString mgin;

  //! Generic Workspace Input.
  ArrayOfIndex mgintype;

  //! Generic Workspace Input Types (Contains the valid types if the method
  // is not completely supergeneric but only defined for these types).
  ArrayOfArrayOfIndex mginspectype;

  //! Generic Workspace Input Defaults.
  ArrayOfString mgindefault;

  //! Generic Workspace Input Description.
  ArrayOfString mgindesc;

  //! Indexes of Input-only variables.
  ArrayOfIndex minonly;

  //! Indexes of Output-only variables.
  ArrayOfIndex moutonly;

  //! Indexes of Input-Output variables.
  ArrayOfIndex minout;

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

  //! Flag, whether method implementation relies on templates.
  /*!
    If we want to implement a supergeneric method by a template
    function, we must not generate explicit method headers in auto_md.h.
  */
  bool muses_templates;

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

  //! The actual groups of a supergeneric method.
  /*! 
    This holds the actual groups after expansion of a supergeneric method.
  */
  String mactual_groups;
};

void define_md_data_raw();

void expand_md_data_raw_to_md_data();

bool format_paragraph(String& s,
                      const String& indent,
                      const size_t linelen,
                      const size_t offset = 0);

void define_md_map();

void define_md_raw_map();

void get_short_wsv_description(String& s, const String& desc);

#endif  // methods_h
