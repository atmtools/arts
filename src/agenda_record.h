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
  AgRecord(const char name[],
           const char description[],
           const ArrayOfString& output,
           const ArrayOfString& input);

  const String& Name() const { return mname; }
  const String& Description() const { return mdescription; }
  const ArrayOfIndex& Out() const { return moutput; }
  const ArrayOfIndex& In() const { return minput; }

  //! Assignment operator.
  /*! To override the default assignment operator. AgRecords cannot be
      assigned! */
  AgRecord operator=(const AgRecord& /* m */) {
    cout << "AgRecord cannot be assigned!\n";
    arts_exit();
    return AgRecord();
  }

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

ostream& operator<<(ostream& os, const AgRecord& agr);

#endif
