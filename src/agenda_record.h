/*!
  \file   agenda_record.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Thu Mar 14 08:49:33 2002
  
  \brief  Declarations for AgRecord, storing lookup information for
  one agenda. 
*/

#ifndef agenda_record_h
#define agenda_record_h

#include <iostream>
#include "arts.h"
#include "mystring.h"

//! Lookup information for one agenda.
/*! 
  An object of this class contains the documentation for an agenda,
  plus the list of output and input variables.

  This cannot be set in the workspace lookup data, since we need to
  use the handles for the WSVs.
*/

class AgRecord {
 public:
  //! Default constructor.
  AgRecord() : mname(""), mdescription(""), moutput(0), minput(0) {
    // Nothing to do here.
  }

  // Initializing constructor. Implementation in .cc file.
  AgRecord(const char* name,
           const char* description,
           const ArrayOfString& output,
           const ArrayOfString& input);

  AgRecord(const AgRecord&) = default;
  AgRecord& operator=(const AgRecord&) = default;

  [[nodiscard]] const String& Name() const { return mname; }
  [[nodiscard]] const String& Description() const { return mdescription; }
  [[nodiscard]] const ArrayOfIndex& Out() const { return moutput; }
  [[nodiscard]] const ArrayOfIndex& In() const { return minput; }

  friend ostream& operator<<(ostream& os, const AgRecord& agr);

 private:
  //! The name of this agenda.
  String mname;

  //! A text string describing this agenda.
  String mdescription;

  //! Workspace Output.
  ArrayOfIndex moutput;

  //! Workspace Input.
  ArrayOfIndex minput;
};

void define_agenda_data();

void define_agenda_map();

bool check_agenda_data();

void write_agenda_wrapper_header(ofstream& ofs,
                                 const AgRecord& agr,
                                 bool is_agenda_array);

#endif
