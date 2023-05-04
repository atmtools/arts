/**
  \file   wsv_aux.h
  \brief  Auxiliary header stuff related to workspace variable
          groups. Normally you should not need to edit this file. 


  \author Stefan Buehler
  \date   2000-06-10
*/

#ifndef wsv_aux_h
#define wsv_aux_h

#include <tokval.h>

#include "array.h"
#include "arts.h"
#include "exceptions.h"

//! Returns list of ids of the given group names
void get_wsv_group_ids(ArrayOfIndex& ids, String name);

//! Check if group is an agenda group
bool is_agenda_group_id(const Index group_id);

//! Returns the id of the given group
Index get_wsv_group_id(const String& name);

//! Return string list of array types
String get_array_groups_as_string(bool basetype_is_group = false,
                                  bool return_basetype_only = false);

/** This class contains all static information for one workspace
    variable.

    The program make_wsv_h.cc uses these records to generate the file
    wsv.h, which contains both the declaration of the wsv handles and
    the declaration of the workspace itself.

    \author Stefan Buehler */
class WsvRecord {
 public:
  /** Default constructor. */
  WsvRecord();

  /** Initializing constructor.

    This is used by Workspace::define_wsv_data() to set the information for
    each workspace variable. */
  WsvRecord(const char* name,
            const char* description,
            const String& group,
            const TokVal& val);

  /** Initializing constructor.

    This is used by the parser to create automatically allocated variables */
  WsvRecord(const char* name,
            const char* description,
            const Index group,
            const TokVal& val);

  /** Initializing constructor.

    This is used by Workspace::define_wsv_data() to set the information for
    each workspace variable. */
  WsvRecord(const char* name, const char* description, const String& group);

  /** Initializing constructor.

    This is used by the parser to create automatically allocated variables */
  WsvRecord(const char* name, const char* description, const Index group);

  /** Name of this workspace variable. */
  [[nodiscard]] const String& Name() const { return mname; }

  /** A text describing this workspace variable. */
  [[nodiscard]] const String& Description() const { return mdescription; }

  /** The wsv group to which this variable belongs. */
  [[nodiscard]] Index Group() const { return mgroup; }

  [[nodiscard]] bool has_defaults() const;

  [[nodiscard]] std::shared_ptr<void> get_copy() const;

  [[nodiscard]] const TokVal& default_value() const { return defval; }
  void update_default_value(ArtsType auto&& v) {defval=std::forward<decltype(v)>(v);}

 private:
  String mname;

  String mdescription;

  Index mgroup{-1};

  TokVal defval;
};

#endif  // wsv_aux_h
