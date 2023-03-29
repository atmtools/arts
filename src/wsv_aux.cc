/* Copyright (C) 2022 Oliver Lemke <oliver.lemke@uni-hamburg.de>

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
  \file   wsv_aux.cc
  \author Oliver Lemke <oliver.lemke@uni-hamburg.de>
  \date   Thu Jul 14 14:43:37 CEST 2022
  
  \brief  Implementation of WSV aux functions.
*/

#include "wsv_aux.h"

/** Default constructor. */
WsvRecord::WsvRecord() : mname(), mdescription(), defval() { /* Nothing to do here */
}

/** Initializing constructor.

    This is used by Workspace::define_wsv_data() to set the information for
    each workspace variable. */
WsvRecord::WsvRecord(const char* name,
          const char* description,
          const String& group,
          const TokVal& val)
    : mname(name), mdescription(description), defval(val) {
  // Map the group names to groups' indexes
  mgroup = get_wsv_group_id(group);
  if (mgroup == -1) {
    ostringstream os;

    os << "Unknown WSV Group " << group << " WSV " << mname;
    throw runtime_error(os.str());
  }
}

/** Initializing constructor.

    This is used by the parser to create automatically allocated variables */
WsvRecord::WsvRecord(const char* name,
          const char* description,
          const Index group,
          const TokVal& val)
    : mname(name), mdescription(description), mgroup(group), defval(val) {
  // Nothing to do here
}

/** Initializing constructor.

    This is used by Workspace::define_wsv_data() to set the information for
    each workspace variable. */
WsvRecord::WsvRecord(const char* name, const char* description, const String& group)
    : mname(name), mdescription(description), defval() {
  // Map the group names to groups' indexes
  mgroup = get_wsv_group_id(group);
  if (mgroup == -1) {
    ostringstream os;

    os << "Unknown WSV Group " << group << " WSV " << mname;
    throw runtime_error(os.str());
  }
}

/** Initializing constructor.

    This is used by the parser to create automatically allocated variables */
WsvRecord::WsvRecord(const char* name, const char* description, const Index group)
    : mname(name), mdescription(description), mgroup(group), defval() {
  // Nothing to do here
}

bool WsvRecord::has_defaults() const {
  return not defval.holdsAny();
}

std::shared_ptr<void> WsvRecord::get_copy() const {
  if (has_defaults())
    return defval.copy_value();
  return nullptr;
}
