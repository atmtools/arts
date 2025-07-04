#pragma once

#include <workspace.h>

uint32_t hlist_num_cols(const std::vector<String>& v,
                        const std::vector<String>& v2 = {});

bool str_compare_nocase(const std::string& lhs, const std::string& rhs);

std::string fix_newlines(std::string x);

String unwrap_stars(const String&);

String get_agenda_io(const String&);

String short_doc(const String& x);

String method_docs(const String& name);

String variable_used_by(const String& name);

String to_defval_str(const Wsv& wsv, std::string_view x);

namespace Python {
String group_generics_inout(const String& group);

String group_workspace_types(const String& group);
}  // namespace Python
