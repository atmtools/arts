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
