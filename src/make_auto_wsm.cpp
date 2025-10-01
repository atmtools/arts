#include <auto_wsg.h>
#include <debug.h>
#include <mystring.h>

#include <algorithm>
#include <exception>
#include <format>
#include <iostream>
#include <ranges>
#include <stdexcept>
#include <string>

#include "workspace_group_friends.h"
#include "workspace_groups.h"
#include "workspace_meta_methods.h"
#include "workspace_methods.h"
#include "workspace_variables.h"

namespace {
void scan_wsmr_for_errors(ArrayOfString& errors,
                          const std::string& name,
                          const WorkspaceMethodInternalRecord& wsmr) {
  const auto& wsv = internal_workspace_variables();

  if (wsmr.desc.size() == 0) {
    errors.push_back(std::format("No description for \"{}\"", name));
  }

  if (wsmr.desc.back() not_eq '\n') {
    errors.push_back(
        std::format("Description for \"{}\" ends without a newline", name));
  }

  if (wsmr.author.size() == 0) {
    errors.push_back(std::format("No authors for \"{}\"", name));
  }

  if (wsmr.gout.size() not_eq wsmr.gout_type.size() or
      wsmr.gout.size() not_eq wsmr.gout_desc.size()) {
    errors.push_back(std::format(
        "Mismatch sizes of gout ({}), gout_type ({}), gout_des ({}), for \"{}\"",
        wsmr.gout.size(),
        wsmr.gout_type.size(),
        wsmr.gout_desc.size(),
        name));
  }

  for (auto& type : wsmr.gout_type) {
    if (not valid_wsg(type)) {
      errors.push_back(
          std::format(R"(Invalid gout_type for "{}": "{}")", name, type));
    }
  }

  for (auto& gvar : wsmr.gout) {
    if (wsv.find(gvar) not_eq wsv.end()) {
      errors.push_back(std::format(
          R"(Invalid gout for "{}": "{}" - it is already a WSV)", name, gvar));
    }
  }

  if (wsmr.gin.size() not_eq wsmr.gin_type.size() or
      wsmr.gin.size() not_eq wsmr.gin_value.size() or
      wsmr.gin.size() not_eq wsmr.gin_desc.size()) {
    errors.push_back(std::format(
        "Mismatch sizes of gin ({}), gin_type ({}), gin_value ({}), gin_desc ({}), for \"{}\"",
        wsmr.gin.size(),
        wsmr.gin_type.size(),
        wsmr.gin_value.size(),
        wsmr.gin_desc.size(),
        name));
  }

  for (auto& type : wsmr.gin_type) {
    if (not valid_wsg(type)) {
      errors.push_back(
          std::format(R"(Invalid gin_type for "{}": "{}")", '"', name, type));
    }
  }

  for (auto& gvar : wsmr.gin) {
    if (wsv.find(gvar) not_eq wsv.end()) {
      errors.push_back(std::format(
          R"(Invalid gin for "{}": "{}" - it is already a WSV)", name, gvar));
    }
  }

  for (auto& var : wsmr.out) {
    if (wsv.find(var) == wsv.end()) {
      errors.push_back(std::format(
          R"(Invalid out for "{}": "{}" - it is not a WSV)", name, var));
    }
  }

  for (auto& var : wsmr.in) {
    if (wsv.find(var) == wsv.end()) {
      errors.push_back(std::format(
          R"(Invalid in for "{}": "{}" - it is not a WSV)", name, var));
    }
  }
}

std::vector<std::vector<std::string>> overloads(
    const WorkspaceMethodInternalRecord& wsmr) {
  std::vector<std::vector<std::string>> out;

  for (auto& str : wsmr.gout_type) {
    if (std::any_of(str.begin(), str.end(), Cmp::eq(',')))
      out.push_back(split(str, ","));
  }

  for (auto& str : wsmr.gin_type) {
    if (std::any_of(str.begin(), str.end(), Cmp::eq(',')))
      out.push_back(split(str, ","));
  }

  for (auto& v : out) {
    for (auto& s : v) {
      trim(s);
    }
  }

  return out;
}

std::optional<WorkspaceMethodInternalRecord> make_overload(
    const WorkspaceMethodInternalRecord& wsmr, std::size_t i) {
  const auto ol = overloads(wsmr);
  if (ol.size() == 0) return std::nullopt;
  if (i >= ol.front().size()) return std::nullopt;

  WorkspaceMethodInternalRecord owsmr(wsmr);

  std::size_t ol_i = 0;
  for (std::size_t garg = 0; garg < owsmr.gout_type.size(); garg++) {
    if (std::any_of(wsmr.gout_type[garg].begin(),
                    wsmr.gout_type[garg].end(),
                    Cmp::eq(','))) {
      owsmr.gout_type[garg] = ol[ol_i][i];
      ol_i++;
    }
  }
  for (std::size_t garg = 0; garg < owsmr.gin_type.size(); garg++) {
    if (std::any_of(wsmr.gin_type[garg].begin(),
                    wsmr.gin_type[garg].end(),
                    Cmp::eq(','))) {
      owsmr.gin_type[garg] = ol[ol_i][i];
      ol_i++;
    }
  }
  return owsmr;
}

ArrayOfString scan_for_errors() {
  const auto& wsm = internal_workspace_methods();

  ArrayOfString errors{};

  for (auto& [name, wsmr] : wsm) {
    if (wsmr.has_overloads()) {
      std::size_t i = 0;
      auto owsmr    = make_overload(wsmr, i);
      while (owsmr) {
        scan_wsmr_for_errors(errors, name, *owsmr);
        owsmr = make_overload(wsmr, ++i);
      }
    } else {
      scan_wsmr_for_errors(errors, name, wsmr);
    }
  }

  return errors;
}

void signature(std::ostream& os,
               const std::string& name,
               const WorkspaceMethodInternalRecord& wsmr) try {
  const int n = wsmr.count_overloads();

  os << wsmr.docstring() << '\n';

  for (int i = 0; i < n; i++) os << wsmr.header(name, i) << ";\n";
} catch (std::exception& e) {
  throw std::runtime_error("Error in signature():\n\n" + std::string(e.what()));
}

void header(std::ostream& os) try {
  const auto& wsm = internal_workspace_methods();

  os << R"--(#pragma once

//! auto-generated by make_auto_wsm.cpp

#include <string>
#include <vector>

#include "auto_wsg.h"

struct Workspace;

struct WorkspaceMethodRecord {
  std::vector<std::string> out;
  std::vector<std::string> in;
  std::unordered_map<std::string, Wsv> defs;
  std::function<void(Workspace&, const std::vector<std::string>&, const std::vector<std::string>&)> func;
};

const std::unordered_map<std::string, WorkspaceMethodRecord>& workspace_methods();

)--";

  for (auto& [name, wsmr] : wsm) {
    signature(os, name, wsmr);
    os << '\n';
  }
} catch (std::exception& e) {
  throw std::runtime_error("Error in header():\n\n" + std::string(e.what()));
}

bool is_generic(std::string_view type) {
  return type == "Any"sv or type.contains(',');
}

std::string first_generic(const std::vector<std::string>& gout_type,
                          const std::vector<std::string>& gin_type) {
  for (std::size_t i = 0; i < gout_type.size(); i++) {
    if (is_generic(gout_type[i])) return std::format("out[{}]", i);
  }

  for (std::size_t i = 0; i < gin_type.size(); i++) {
    if (is_generic(gin_type[i])) return std::format("in[{}]", i);
  }

  return "INVALID";
}

std::vector<std::string> all_generics(std::string_view type) {
  std::vector<std::string> out{};

  if (type == "Any"sv) {
    for (const auto& group : internal_workspace_groups()) {
      out.push_back(group.first);
    }
  } else {
    out = split(std::string{type}, ",");
    for (auto& s : out) {
      trim(s);
    }
  }

  // trailing commas
  while (not out.empty() and out.back().empty()) out.pop_back();

  return out;
}

template <typename T>
std::size_t max_elem(const std::vector<T>& vecs) {
  return std::ranges::max_element(vecs, {}, &T::size)->size();
}

std::vector<std::vector<std::string>> all_gets(
    const WorkspaceMethodInternalRecord& wsmr) {
  std::vector<std::vector<std::string>> out;

  out.reserve(wsmr.gout_type.size() + wsmr.gin_type.size() + wsmr.out.size() +
              wsmr.in.size());

  if (wsmr.pass_workspace) out.emplace_back().emplace_back("ws");

  std::size_t COUNT = 0;
  for (const auto& wsm : wsmr.out) {
    out.emplace_back().emplace_back(std::format(
        "ws.get{1}<{2}>(out[{0}])",
        COUNT++,
        internal_workspace_variables().at(wsm).type,
        std::ranges::any_of(wsmr.in, Cmp::eq(wsm)) ? ""sv : "_or"sv));
  }

  for (const auto& type : wsmr.gout_type) {
    const auto name = std::format("out[{}]", COUNT++);
    auto vec        = all_generics(type);
    if (vec.empty()) {
      out.emplace_back().emplace_back(
          std::format("ws.get_or<{0}>({1})", type, name));
    } else {
      auto& v = out.emplace_back();
      for (const auto& s : vec)
        v.emplace_back(std::format("ws.get_or<{0}>({1})", s, name));
    }
  }

  COUNT = 0;
  for (const auto& wsm : wsmr.in) {
    const auto name = std::format("in[{}]", COUNT++);
    if (std::ranges::any_of(wsmr.out, Cmp::eq(wsm))) continue;
    out.emplace_back().emplace_back(std::format(
        "ws.get<{0}>({1})", internal_workspace_variables().at(wsm).type, name));
  }

  for (std::size_t i = 0; i < wsmr.gin_type.size(); i++) {
    const auto& type = wsmr.gin_type[i];
    const auto name  = std::format("in[{}]", COUNT++);
    if (std::ranges::any_of(wsmr.gout, Cmp::eq(wsmr.gin[i]))) continue;

    auto vec = all_generics(type);
    if (vec.empty()) {
      out.emplace_back().emplace_back(
          std::format("ws.get<{0}>({1})", type, name));
    } else {
      auto& v = out.emplace_back();
      for (const auto& s : vec)
        v.emplace_back(std::format("ws.get<{0}>({1})", s, name));
    }
  }

  const std::size_t max_size = max_elem(out);

  for (auto& vec : out) {
    while (vec.size() < max_size) {
      vec.push_back(vec.front());
    }
  }

  return out;
}

std::string error_signature(const std::string& name,
                            const WorkspaceMethodInternalRecord& wsmr) {
  const auto& wsvs = internal_workspace_variables();

  std::size_t type_spaces = 0;
  for (auto& wsv : wsmr.out) type_spaces = std::max(type_spaces, wsv.size());
  for (auto& wsv : wsmr.in) type_spaces = std::max(type_spaces, wsv.size());
  for (auto& wsv : wsmr.gout) type_spaces = std::max(type_spaces, wsv.size());
  for (auto& wsv : wsmr.gin) type_spaces = std::max(type_spaces, wsv.size());
  const std::string spaces(name.size() + 3, ' ');

  std::string out = "  " + name + "(";

  std::string sep       = "";
  std::string other_sep = ",\n" + spaces;

  for (auto& wsv : wsmr.out) {
    out += std::format("{}{}{} : {}",
                       std::exchange(sep, other_sep),
                       wsv,
                       std::string(type_spaces - wsv.size(), ' '),
                       wsvs.at(wsv).type);
  }

  for (auto& wsv : wsmr.in) {
    if (std::ranges::any_of(wsmr.out, Cmp::eq(wsv))) continue;
    out += std::format("{}{}{} : {}",
                       std::exchange(sep, other_sep),
                       wsv,
                       std::string(type_spaces - wsv.size(), ' '),
                       wsvs.at(wsv).type);
  }

  for (std::size_t i = 0; i < wsmr.gout.size(); i++) {
    out += std::format("{}{}{} : {}",
                       std::exchange(sep, other_sep),
                       wsmr.gout[i],
                       std::string(type_spaces - wsmr.gout[i].size(), ' '),
                       wsmr.gout_type[i]);
  }

  for (std::size_t i = 0; i < wsmr.gin.size(); i++) {
    if (std::ranges::any_of(wsmr.gout, Cmp::eq(wsmr.gin[i]))) continue;
    out += std::format("{}{}{} : {}",
                       std::exchange(sep, other_sep),
                       wsmr.gin[i],
                       std::string(type_spaces - wsmr.gin[i].size(), ' '),
                       wsmr.gin_type[i]);
  }

  return out + ')';
}

void call_function(std::ostream& os,
                   const std::string& name,
                   const WorkspaceMethodInternalRecord& wsmr) try {
  const auto& wsv = internal_workspace_variables();

  if (wsmr.has_any()) {
    const std::string first = first_generic(wsmr.gout_type, wsmr.gin_type);

    os << "[](Workspace& ws [[maybe_unused]], const std::vector<std::string>& out [[maybe_unused]], const std::vector<std::string>& in [[maybe_unused]]) {\n";
    os << "    try {\n";

    std::println(os, "      auto& _first = ws.share({0});", first);

    const auto gets = all_gets(wsmr);
    std::size_t i   = 0;

    // one method per group
    std::println(os, "      switch (_first.value_index()) {{");
    for (auto& group : internal_workspace_groups() | std::views::keys) {
      std::print(
          os,
          R"(        case WorkspaceGroupInfo<{0}>::index: return {1}({2})",
          group,
          name,
          gets.front()[i]);
      for (auto& v : gets | std::views::drop(1)) {
        std::print(os, ", {}", v[i]);
      }
      ++i;
      std::println(os, ");");
    }
    std::println(os,
                 R"(      }}
    }} catch (std::exception& e) {{
      throw std::runtime_error(std::format(R"-x-(Error in agenda call to generic method

{}

{{}})-x-", e.what()));
    }}
  }})",
                 error_signature(name, wsmr));
  } else if (wsmr.has_overloads()) {
    os << "[map = std::unordered_map<std::string_view, std::function<void(Workspace&, const std::vector<std::string>&, const std::vector<std::string>&)>>{\n";

    const auto ol = overloads(wsmr);

    bool first = true;
    for (std::size_t i = 0; i < ol.front().size(); i++) {
      WorkspaceMethodInternalRecord owsmr(wsmr);

      std::size_t ol_i = 0;
      for (std::size_t garg = 0; garg < owsmr.gout_type.size(); garg++) {
        if (std::any_of(wsmr.gout_type[garg].begin(),
                        wsmr.gout_type[garg].end(),
                        Cmp::eq(','))) {
          owsmr.gout_type[garg] = ol[ol_i][i];
          ol_i++;
        }
      }
      for (std::size_t garg = 0; garg < owsmr.gin_type.size(); garg++) {
        if (std::any_of(wsmr.gin_type[garg].begin(),
                        wsmr.gin_type[garg].end(),
                        Cmp::eq(','))) {
          owsmr.gin_type[garg] = ol[ol_i][i];
          ol_i++;
        }
      }
      os << comma(first, "     ");
      os << "{\"";
      bool other_first = true;
      for (const auto& garg : ol) {
        os << comma(other_first, "") << garg[i];
      }
      os << "\", ";
      call_function(os, name, owsmr);
      os << "}";
    }

    os << R"--(}] (Workspace& ws [[maybe_unused]],const std::vector<std::string>& out [[maybe_unused]], const std::vector<std::string>& in [[maybe_unused]]) {
      const auto& func = map.at()--";

    bool final_first = true;
    for (std::size_t garg = 0; garg < wsmr.gout_type.size(); garg++) {
      if (std::any_of(wsmr.gout_type[garg].begin(),
                      wsmr.gout_type[garg].end(),
                      Cmp::eq(','))) {
        if (not final_first) os << ", \", \", ";
        os << "ws.share(out[" << garg + wsmr.out.size() << "]).type_name()";
        final_first = false;
      }
    }
    for (std::size_t garg = 0; garg < wsmr.gin_type.size(); garg++) {
      if (std::any_of(wsmr.gin_type[garg].begin(),
                      wsmr.gin_type[garg].end(),
                      Cmp::eq(','))) {
        if (not final_first) os << ", \", \", ";
        os << "ws.share(in[" << garg + wsmr.in.size() << "]).type_name()";
        final_first = false;
      }
    }

    os << R"--();
      func(ws, out, in);
    }
)--";
  } else {
    os << "[](Workspace& ws [[maybe_unused]], const std::vector<std::string>& out [[maybe_unused]], const std::vector<std::string>& in [[maybe_unused]]) {\n";
    os << "    try {\n";

    bool first = true;
    os << "      " << name << "(";
    if (wsmr.pass_workspace) {
      os << "ws";
      first = false;
    }

    const String spaces(name.size() + 7, ' ');

    int out_count = 0;
    for (auto& str : wsmr.out) {
      os << comma(first, spaces) << "ws.get";
      if (std::count(wsmr.in.begin(), wsmr.in.end(), str) == 0) os << "_or";
      os << "<" << wsv.at(str).type << ">(out[" << out_count++
         << "]) /* out */";
    }

    for (std::size_t i = 0; i < wsmr.gout.size(); i++) {
      os << comma(first, spaces) << "ws.get_or<"
         << any_is_typename(wsmr.gout_type[i]) << ">(out[" << out_count++
         << "]) /* gout */";
    }

    int in_count = 0;
    for (auto& str : wsmr.in) {
      if (std::any_of(wsmr.out.begin(), wsmr.out.end(), [&str](auto& var) {
            return str == var;
          })) {
        in_count++;
        continue;
      }
      os << comma(first, spaces) << "ws.get<" << wsv.at(str).type << ">(in["
         << in_count++ << "]) /* in */";
    }

    for (std::size_t i = 0; i < wsmr.gin.size(); i++) {
      os << comma(first, spaces) << "ws.get<" << wsmr.gin_type[i] << ">(in["
         << in_count++ << "]) /* gin */";
    }

    os << "\n      );\n"
          "    } catch (std::exception& e) {\n"
          "      throw std::runtime_error(std::format(R\"-x-(Error in agenda call to specific method\n\n"
       << error_signature(name, wsmr)
       << "\n\n{})-x-\", e.what()));\n"
          "    }\n  }\n";
  }
} catch (std::exception& e) {
  throw std::runtime_error("Error in call_function():\n\n" +
                           std::string(e.what()));
}

void wsm_record(std::ostream& os,
                const std::string& name,
                const WorkspaceMethodInternalRecord& wsmr) try {
  os << "    .out={";
  bool first = true;
  for (auto& str : wsmr.out) {
    os << comma(first, "          ") << "\"" << str << "\"";
  }
  for (auto& str : wsmr.gout) {
    os << comma(first, "          ") << "\"_" << str << "\"";
  }
  os << "},\n";
  os << "    .in={";
  first = true;
  for (auto& str : wsmr.in) {
    os << comma(first, "         ") << "\"" << str << "\"";
  }
  for (auto& str : wsmr.gin) {
    os << comma(first, "         ") << "\"_" << str << "\"";
  }
  os << "},\n";
  os << "    .defs={";
  first = true;
  for (std::size_t i = 0; i < wsmr.gin.size(); i++) {
    if (wsmr.gin_value[i]) {
      os << comma(first, "           ") << "{\"_" << wsmr.gin[i]
         << "\", internal_workspace_methods().at(\"" << name << "\").gin_value["
         << i << "].value()}";
    }
  }
  os << "},\n";

  os << "    .func=";
  call_function(os, name, wsmr);

  os << "\n";
} catch (std::exception& e) {
  throw std::runtime_error("Error in wsm_record():\n\n" +
                           std::string(e.what()));
}

void meta_implementation(std::ostream& os) try {
  const auto& wsm = internal_workspace_methods();

  os << R"--(//! auto-generated by make_auto_wsm.cpp

#include <workspace.h>
#include <workspace_meta_methods.h>

static const auto& _wsmmeta = internal_meta_methods();

)--";

  for (auto& m : internal_meta_methods()) {
    os << m.call(wsm) << '\n' << '\n';
  }
} catch (std::exception& e) {
  throw std::runtime_error("Error in meta_implementation:\n\n" +
                           std::string(e.what()));
}

std::ofstream& select_ofstream(std::vector<std::ofstream>& ofs, int i) {
  return ofs[i % ofs.size()];
}

void implementation(std::ostream& os, const int n) try {
  const auto& wsm = internal_workspace_methods();

  std::vector<std::ofstream> ofs;
  ofs.reserve(n);
  for (int i = 0; i < n; i++) {
    ofs.emplace_back(std::format("auto_wsm_{}.cpp", i));
  }

  for (auto& of : ofs) {
    of << R"--(//! auto-generated by make_auto_wsm.cpp

#include <iomanip>

#include "auto_wsm.h"

#include "workspace_methods.h"

#include "workspace_agenda_class.h"
#include "workspace_method_class.h"

#include "workspace_class.h"

#include <m_ignore.h>
#include <m_xml.h>

)--";
  }

  int i = 0;
  for (auto& [name, wsmr] : wsm) {
    select_ofstream(ofs, i)
        << "WorkspaceMethodRecord _wsm_" << name << "() {\n  return {\n";
    wsm_record(select_ofstream(ofs, i), name, wsmr);
    select_ofstream(ofs, i) << "  };\n}\n\n";
    i++;
  }

  for (i = 0; i < n; i++) {
    select_ofstream(ofs, i)
        << "std::unordered_map<std::string, WorkspaceMethodRecord> get_workspace_methods"
        << i << R"--(() {
  return {
)--";
  }

  i = 0;
  for (auto& [name, wsmr] : wsm) {
    select_ofstream(ofs, i++)
        << "    { \"" << name << "\", _wsm_" << name << "()},\n";
  }

  for (auto& of : ofs)
    of << R"--(  };
}
)--";

  os << "#include \"auto_wsv.h\"\n\n";
  os << "#include <workspace.h>\n\n";
  for (i = 0; i < n; i++) {
    os << "std::unordered_map<std::string, WorkspaceMethodRecord> get_workspace_methods"
       << i << "();\n";
  }
  os << R"--(
const std::unordered_map<std::string, WorkspaceMethodRecord>& workspace_methods() {
  static std::unordered_map<std::string, WorkspaceMethodRecord> wsm;
  static bool first = true;
  if (first) {
    first = false;
)--";

  for (i = 0; i < n; i++) {
    os << "    const std::unordered_map<std::string, WorkspaceMethodRecord> m"
       << i << " = get_workspace_methods" << i << "();\n    const std::size_t n"
       << i << " = m" << i << ".size();\n";
  }
  os << "\n    wsm.reserve(";
  for (i = 0; i < n; i++) {
    os << "n" << i << " + ";
  }
  os << "0);\n\n";
  for (i = 0; i < n; i++) {
    os << R"--(    for (auto& [k, v] : m)--" << i << R"--() {  wsm[k] = v; }
)--";
  }

  os << R"--(
  }
  return wsm;
}
)--";
} catch (std::exception& e) {
  throw std::runtime_error("Error in implementation():\n\n" +
                           std::string(e.what()));
}

void various_checks_or_throw() {
  ArrayOfString errors{};
  for (auto& [name, wsmr] : internal_workspace_methods()) {
    bool is_void = wsmr.return_type == "void";
    if (not is_void) {
      bool is_workspace = wsmr.return_type == "Workspace";
      bool is_wsg = internal_workspace_groups().contains(wsmr.return_type);
      bool is_wsg_friend = workspace_group_friends().contains(wsmr.return_type);
      if (not is_workspace and not is_wsg and not is_wsg_friend) {
        errors.push_back(std::format(
            "Method {} has a return type.  Its return type \"{}\" is not void or Workspace, and it is not a workspace group (or friend of one).",
            name,
            wsmr.return_type));
      }

      if (wsmr.return_desc.empty()) {
        errors.push_back(std::format(
            "Method {} has a return type.  It lacks a return description.",
            name));
      }
    }
  }

  if (not errors.empty()) {
    throw std::runtime_error(
        std::format("Errors found in workspace methods:\n\n{:n}", errors));
  }
}
}  // namespace

int main(int argc, char** argv) try {
  if (argc != 2) throw std::runtime_error("usage: make_auto_wsm <num_methods>");

  various_checks_or_throw();

  const int num_methods = std::stoi(argv[1]);
  std::ofstream head("auto_wsm.h");
  std::ofstream impl("auto_wsm.cpp");
  std::ofstream of("auto_wsmmeta.cpp");

  ARTS_USER_ERROR_IF(const auto err = scan_for_errors();
                     not err.empty(), "{}", err)

  header(head);
  implementation(impl, num_methods);
  meta_implementation(of);
} catch (std::exception& e) {
  std::cerr << "Cannot create the automatic methods with error:\n\n"
            << e.what() << '\n';
  return 1;
}
