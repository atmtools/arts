#pragma once

#include "mystring.h"

String unwrap_stars(String x);

namespace Python {
String group_generics_inout(const String& group);

String group_workspace_types(const String& group);
}  // namespace Python
