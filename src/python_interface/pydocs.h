#pragma once

#include "mystring.h"

String unwrap_stars(String);

String add_type(String, const String&);

String get_agenda_io(const String&);

namespace Python {
String group_generics_inout(const String& group);

String group_workspace_types(const String& group);
}  // namespace Python
