#pragma once

#include "mystring.h"

String unwrap_stars(String);

String get_agenda_io(const String&);

String short_doc(const String& x);

String method_docs(const String& name);

String to_defval_str(const String& x);

namespace Python {
String group_generics_inout(const String& group);

String group_workspace_types(const String& group);
}  // namespace Python
