#pragma once

#include "mystring.h"

uint32_t hlist_num_cols(const std::vector<String>& v);

bool str_compare_nocase(const std::string& lhs, const std::string& rhs);

String unwrap_stars(String);

String get_agenda_io(const String&);

String short_doc(const String& x);

String method_docs(const String& name);

String to_defval_str(const String& x, const String& group);

namespace Python {
String group_generics_inout(const String& group);

String group_workspace_types(const String& group);
}  // namespace Python
