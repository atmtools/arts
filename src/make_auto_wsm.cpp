#include <auto_wsg.h>

#include <algorithm>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "auto_wsv.h"
#include "workspace_agendas.h"
#include "workspace_methods.h"

const static auto wsm = internal_workspace_methods();
const static auto wsv = workspace_variables();

bool scan_wsmr_for_errors(const std::string& name,
                          const WorkspaceMethodInternalRecord& wsmr) {
  bool any_errors = false;
  if (wsmr.desc.size() == 0) {
    any_errors = true;
    std::cerr
        << R"(static_assert(false, "\n\n    Methods cannot be compiled\n    No description for \")"
        << name << R"(\"\n");

)";
    return true;
  }

  if (wsmr.desc.back() not_eq '\n') {
    any_errors = true;
    std::cerr
        << R"(static_assert(false, "\n\n    Methods cannot be compiled\n    No final newline for description of \")"
        << name << R"(\"\n");

)";
  }

  if (wsmr.author.size() == 0) {
    any_errors = true;
    std::cerr
        << R"(static_assert(false, "\n\n    Methods cannot be compiled\n    No authors for \")"
        << name << R"(\"\n");

)";
  }

  if (wsmr.gout.size() not_eq wsmr.gout_type.size() or
      wsmr.gout.size() not_eq wsmr.gout_desc.size()) {
    any_errors = true;
    std::cerr
        << R"(static_assert(false, "\n\n    Methods cannot be compiled\n    gout, gout_type, and gout_desc size mismatch for \")"
        << name << R"(\"");

)";
  }

  for (auto& type : wsmr.gout_type) {
    if (not valid_wsg(type)) {
      any_errors = true;
      std::cerr
          << R"(static_assert(false, "\n\n    Methods cannot be compiled\n    Invalid gout_type for \")"
          << name << R"(\": \")" << type << R"(\"\n");
  
)";
    }
  }

  for (auto& gvar : wsmr.gout) {
    if (wsv.find(gvar) not_eq wsv.end()) {
      any_errors = true;
      std::cerr
          << R"(static_assert(false, "\n\n    Methods cannot be compiled\n    Invalid gout for \")"
          << name << R"(\": \")" << gvar
          << R"(\", the gout is already a WSV\n");
  
)";
    }
  }

  if (wsmr.gin.size() not_eq wsmr.gin_type.size() or
      wsmr.gin.size() not_eq wsmr.gin_value.size() or
      wsmr.gin.size() not_eq wsmr.gin_desc.size()) {
    any_errors = true;
    std::cerr
        << R"(static_assert(false, "\n\n    Methods cannot be compiled\n    gin, gin_type, gin_value, and gin_desc size mismatch for \")"
        << name << R"(\"\n");

)";
  }

  for (auto& type : wsmr.gin_type) {
    if (not valid_wsg(type)) {
      any_errors = true;
      std::cerr
          << R"(static_assert(false, "\n\n    Methods cannot be compiled\n    Invalid gin_type for \")"
          << name << R"(\": \")" << type << R"(\"\n");
  
)";
    }
  }

  for (auto& gvar : wsmr.gin) {
    if (wsv.find(gvar) not_eq wsv.end()) {
      any_errors = true;
      std::cerr
          << R"(static_assert(false, "\n\n    Methods cannot be compiled\n    Invalid gin for \")"
          << name << R"(\": \")" << gvar << R"(\", the gin is already a WSV\n");
  
)";
    }
  }

  for (auto& var : wsmr.out) {
    if (wsv.find(var) == wsv.end()) {
      any_errors = true;
      std::cerr
          << R"(static_assert(false, "\n\n    Methods cannot be compiled\n    Invalid out for \")"
          << name << R"(\": \")" << var << R"(\", the variable is not a WSV\n");

)";
    }
  }

  for (auto& var : wsmr.in) {
    if (wsv.find(var) == wsv.end()) {
      any_errors = true;
      std::cerr
          << R"(static_assert(false, "\n\n    Methods cannot be compiled\n    Invalid in for \")"
          << name << R"(\": \")" << var << R"(\", the variable is not a WSV\n");

)";
    }
  }
  return any_errors;
}

bool has_overloads(const WorkspaceMethodInternalRecord& wsmr) {
  for (auto& str : wsmr.gout_type) {
    if (std::any_of(str.begin(), str.end(), Cmp::eq(','))) return true;
  }

  for (auto& str : wsmr.gin_type) {
    if (std::any_of(str.begin(), str.end(), Cmp::eq(','))) return true;
  }

  return false;
}

std::vector<std::string> split(const std::string& str, char delim) {
  std::vector<std::string> out;
  std::stringstream ss(str);
  std::string item;
  while (std::getline(ss, item, delim)) {
    out.push_back(item);
  }

  for (auto& arg : out) {
    while (arg.front() == ' ') arg.erase(arg.begin());
    while (arg.back() == ' ') arg.pop_back();
  }

  return out;
}

std::vector<std::vector<std::string>> overloads(
    const WorkspaceMethodInternalRecord& wsmr) {
  std::vector<std::vector<std::string>> out;

  for (auto& str : wsmr.gout_type) {
    if (std::any_of(str.begin(), str.end(), Cmp::eq(',')))
      out.push_back(split(str, ','));
  }

  for (auto& str : wsmr.gin_type) {
    if (std::any_of(str.begin(), str.end(), Cmp::eq(',')))
      out.push_back(split(str, ','));
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

bool scan_for_errors() {
  bool any_errors = false;

  for (auto& [name, wsmr] : wsm) {
    if (has_overloads(wsmr)) {
      std::size_t i = 0;
      auto owsmr = make_overload(wsmr, i);
      while (owsmr) {
        any_errors |= scan_wsmr_for_errors(name, *owsmr);
        owsmr = make_overload(wsmr, ++i);
      }
    } else {
      any_errors |= scan_wsmr_for_errors(name, wsmr);
    }
  }

  return any_errors;
}

std::string comma(bool& first, const std::string& spaces = "") {
  if (first) {
    first = false;
    return "";
  }
  return var_string(',', (spaces.size() ? '\n' : ' '), spaces);
}

bool has_any(const WorkspaceMethodInternalRecord& wsmr) {
  return std::count(wsmr.gout_type.begin(), wsmr.gout_type.end(), "Any") +
         std::count(wsmr.gin_type.begin(), wsmr.gin_type.end(), "Any");
}

bool needs_workspace(const WorkspaceMethodInternalRecord& wsmr) {
  return wsmr.pass_workspace;
}

std::string_view any(const std::string& type) {
  if (type == "Any") {
    return "T";
  }
  return type;
}

void docstring(std::ostream& os, const WorkspaceMethodInternalRecord& wsmr) {
  os << "/** " << wsmr.desc << "\n";
  if (needs_workspace(wsmr)) os << "  @param[in] ws Workspace reference\n";

  for (auto& str : wsmr.out) {
    if (std::any_of(wsmr.in.begin(), wsmr.in.end(), [&str](auto& var) {
          return str == var;
        }))
      os << "  @param[inout] " << str << " As WSV" << '\n';
    else
      os << "  @param[out] " << str << " As WSV" << '\n';
  }

  for (std::size_t i = 0; i < wsmr.gout.size(); i++) {
    os << "  @param[out] " << wsmr.gout[i] << " " << wsmr.gout_desc[i] << '\n';
  }

  for (auto& str : wsmr.in) {
    if (std::any_of(wsmr.out.begin(), wsmr.out.end(), [&str](auto& var) {
          return str == var;
        }))
      continue;
    os << "  @param[in] " << str << " As WSV" << '\n';
  }

  for (std::size_t i = 0; i < wsmr.gin.size(); i++) {
    os << "  @param[in] " << wsmr.gin[i] << " " << wsmr.gin_desc[i] << '\n';
  }

  os << " */\n";
}

void signature(std::ostream& os,
               const std::string& name,
               const WorkspaceMethodInternalRecord& wsmr,
               int overload = -1) {
  const std::string spaces(name.size() + 6, ' ');

  if (overload == -1) docstring(os, wsmr);

  if (overload < 0 and has_overloads(wsmr)) {
    const auto ol = overloads(wsmr);

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

      signature(os, name, owsmr, static_cast<int>(i));
    }
    return;
  }

  if (has_any(wsmr)) {
    os << "template <WorkspaceGroup T>\n";
  }

  os << "void " << name << '(';

  bool first = true;
  if (needs_workspace(wsmr)) {
    os << "const Workspace& ws";
    first = false;
  }

  for (auto& str : wsmr.out) {
    os << comma(first, spaces) << wsv.at(str).type << "& " << str;
  }

  for (std::size_t i = 0; i < wsmr.gout.size(); i++) {
    os << comma(first, spaces) << any(wsmr.gout_type[i]) << "& "
       << wsmr.gout[i];
  }

  for (auto& str : wsmr.in) {
    if (std::any_of(wsmr.out.begin(), wsmr.out.end(), [&str](auto& var) {
          return str == var;
        }))
      continue;
    os << comma(first, spaces) << "const " << wsv.at(str).type << "& " << str;
  }

  for (std::size_t i = 0; i < wsmr.gin.size(); i++) {
    os << comma(first, spaces) << "const " << any(wsmr.gin_type[i]) << "& "
       << wsmr.gin[i];
  }

  os << ");\n";
}

void header(std::ostream& os) {
  os << R"--(#pragma once

//! auto-generated by make_auto_wsm.cpp

#include <string>
#include <vector>

#include <auto_wsg.h>

class Workspace;

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
}

void call_function(std::ostream& os,
                   const std::string& name,
                   const WorkspaceMethodInternalRecord& wsmr) {
  if (has_any(wsmr)) {
    const bool any_out = std::any_of(
        wsmr.gout_type.begin(), wsmr.gout_type.end(), Cmp::eq("Any"));

    const std::size_t first_any =
        any_out
            ? (wsmr.out.size() + std::distance(wsmr.gout_type.begin(),
                                               std::find(wsmr.gout_type.begin(),
                                                         wsmr.gout_type.end(),
                                                         "Any")))
            : (wsmr.in.size() + std::distance(wsmr.gin_type.begin(),
                                              std::find(wsmr.gin_type.begin(),
                                                        wsmr.gin_type.end(),
                                                        "Any")));

    // MOSTLY COPY-PASTA

    os << "[](Workspace& ws [[maybe_unused]], const std::vector<std::string>& out [[maybe_unused]], const std::vector<std::string>& in [[maybe_unused]]) {\n";

    bool first = true;

    os << "      std::visit([&](auto& first_any){\n        " << name << "(";
    if (needs_workspace(wsmr)) {
      os << "ws";
      first = false;
    }

    const String spaces(name.size() + 9, ' ');

    int any_count = 0;
    int out_count = 0;
    for (auto& str : wsmr.out) {
      os << comma(first, spaces) << "ws.get";
      if (std::count(wsmr.in.begin(), wsmr.in.end(), str) == 0) os << "_or";
      os << "<" << wsv.at(str).type << ">(out[" << out_count++
         << "]) /* out */";
    }

    for (std::size_t i = 0; i < wsmr.gout.size(); i++) {
      if (wsmr.gout_type[i] == "Any") {
        if (any_count == 0) {
          os << comma(first, spaces) << "*first_any /* gout */";
          out_count++;
        } else {
          os << comma(first, spaces)
             << "ws.get_or<std::remove_cvref_t<decltype(*first_any)>>(out["
             << out_count++ << "]) /* gout */";
        }
        any_count++;
      } else {
        os << comma(first, spaces) << "ws.get_or<" << any(wsmr.gout_type[i])
           << ">(out[" << out_count++ << "]) /* gout */";
      }
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
      if (wsmr.gin_type[i] == "Any") {
        if (any_count == 0) {
          os << comma(first, spaces) << "*first_any /* gin */";
          in_count++;
        } else {
          os << comma(first, spaces)
             << "ws.get<std::remove_cvref_t<decltype(*first_any)>>(in["
             << in_count++ << "]) /* gin */";
        }
        any_count++;
      } else {
        os << comma(first, spaces) << "ws.get<" << wsmr.gin_type[i] << ">(in["
           << in_count++ << "]) /* gin */";
      }
    }

    os << "\n        );\n      }, ";
    if (any_out) {
      os << "ws.share(out[" << first_any << "]) -> value);\n    }";
    } else {
      os << "ws.share(in[" << first_any << "]) -> value);\n    }";
    }

  } else if (has_overloads(wsmr)) {
    os << "[map = std::unordered_map<std::string, std::function<void(Workspace&, const std::vector<std::string>&, const std::vector<std::string>&)>>{\n";

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
      const auto& func = map.at(var_string()--";

    bool final_first = true;
    for (std::size_t garg = 0; garg < wsmr.gout_type.size(); garg++) {
      if (std::any_of(wsmr.gout_type[garg].begin(),
                      wsmr.gout_type[garg].end(),
                      Cmp::eq(','))) {
        if (not final_first) os << ", \", \", ";
        os << "ws.share(out[" << garg + wsmr.out.size() << "]) -> type_name()";
        final_first = false;
      }
    }
    for (std::size_t garg = 0; garg < wsmr.gin_type.size(); garg++) {
      if (std::any_of(wsmr.gin_type[garg].begin(),
                      wsmr.gin_type[garg].end(),
                      Cmp::eq(','))) {
        if (not final_first) os << ", \", \", ";
        os << "ws.share(in[" << garg + wsmr.in.size() << "]) -> type_name()";
        final_first = false;
      }
    }

    os << R"--());
      func(ws, out, in);
    }
)--";
  } else {
    os << "[](Workspace& ws [[maybe_unused]], const std::vector<std::string>& out [[maybe_unused]], const std::vector<std::string>& in [[maybe_unused]]) {\n";

    bool first = true;
    os << "      " << name << "(";
    if (needs_workspace(wsmr)) {
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
      os << comma(first, spaces) << "ws.get_or<" << any(wsmr.gout_type[i])
         << ">(out[" << out_count++ << "]) /* gout */";
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

    os << "\n      );\n    }";
  }
}

void wsm_record(std::ostream& os,
                const std::string& name,
                const WorkspaceMethodInternalRecord& wsmr) {
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
}

void wsm_record(std::ostream& os,
                const std::string& name,
                const WorkspaceAgendaInternalRecord& wsmr) try {
  os << "    .out={";
  bool first = true;
  for (auto& str : wsmr.output) {
    os << comma(first, "          ") << "\"" << str << "\"";
  }
  os << "},\n";
  os << "    .in={";
  first = true;
  for (auto& str : wsmr.input) {
    os << comma(first, "         ") << "\"" << str << "\"";
  }
  os << comma(first, "         ") << '\"' << name << '\"';
  os << "},\n";
  os << "    .defs={},\n";

  os << "    .func=";
  
  os << "[](Workspace& ws, const std::vector<std::string>& out [[maybe_unused]], const std::vector<std::string>& in [[maybe_unused]]) {\n";

  os << "      " << name << "Execute(";
  os << "ws";
  for (std::size_t i=0; i<wsmr.output.size(); i++) {
    os << ",\n" << std::string(14 + name.size(), ' ') << "ws.get_or<"<< wsv.at(wsmr.output[i]).type <<">(out[" << i << "])";
  }
  for (std::size_t i=0; i<wsmr.input.size(); i++) {
    if (std::ranges::any_of(wsmr.output, Cmp::eq(wsmr.input[i]))) continue;
    os << ",\n" << std::string(14 + name.size(), ' ') << "ws.get<"<< wsv.at(wsmr.input[i]).type <<">(out[" << i << "])";
  }
  os << ",\n" << std::string(14 + name.size(), ' ') << "ws.get<";
  if (wsmr.array) os << "ArrayOf";
  os << "Agenda>(in.back())";
  os << "\n      );";
  os << "\n    }";
  os << "\n";
} catch(...) {
  throw std::runtime_error(var_string("Failure in: ", std::quoted(name)));
}

std::ofstream& select_ofstream(std::vector<std::ofstream>& ofs, int i) {
  return ofs[i % ofs.size()];
}

void implementation(std::ostream& os, const int n) {
  std::vector<std::ofstream> ofs;
  ofs.reserve(n);
  for (int i = 0; i < n; i++) {
    ofs.emplace_back("auto_wsm_" + std::to_string(i) + ".cpp");
  }

  for (auto& of : ofs) {
    of << R"--(//! auto-generated by make_auto_wsm.cpp

#include <iomanip>

#include <auto_wsm.h>

#include "workspace_methods.h"

#include "workspace_agenda_class.h"
#include "workspace_method_class.h"

#include "workspace_class.h"

#include <m_copy.h>
#include <m_delete.h>
#include <m_general.h>
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

  i = 0;
  for (auto& [name, wsar]: internal_workspace_agendas()) {
    select_ofstream(ofs, i)
        << "WorkspaceMethodRecord _wsm_" << name << "Execute() {\n  return {\n";
    wsm_record(select_ofstream(ofs, i), name, wsar);
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

  i = 0;
  for (auto& [name, wsar]: internal_workspace_agendas()) {
    select_ofstream(ofs, i)
        << "    { \"" << name << "Execute\", _wsm_" << name << "Execute()},\n";
    i++;
  }

  for (auto& of : ofs)
    of << R"--(  };
}
)--";

  os << "#include <auto_wsv.h>\n\n";
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
}

int main(int argc, char** argv) {
  if (scan_for_errors()) return 0;

  if (argc != 2) throw std::runtime_error("usage: make_auto_wsm <num_methods>");

  const int num_methods = std::stoi(argv[1]);

  header(std::cout);
  implementation(std::cerr, num_methods);
}
