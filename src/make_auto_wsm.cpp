#include <auto_wsg.h>
#include <mystring.h>

#include <algorithm>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "auto_wsv.h"
#include "debug.h"
#include "workspace_meta_methods.h"
#include "workspace_methods.h"
#include "workspace_variables.h"

void scan_wsmr_for_errors(ArrayOfString& errors,
                          const std::string& name,
                          const WorkspaceMethodInternalRecord& wsmr) {
  const auto wsm = internal_workspace_methods();
  const auto& wsv = workspace_variables();

  if (wsmr.desc.size() == 0) {
    errors.push_back(var_string("No description for ", std::quoted(name)));
  }

  if (wsmr.desc.back() not_eq '\n') {
    errors.push_back(var_string(
        "Description for ", std::quoted(name), " does not end with a newline"));
  }

  if (wsmr.author.size() == 0) {
    errors.push_back(var_string("No authors for ", std::quoted(name)));
  }

  if (wsmr.gout.size() not_eq wsmr.gout_type.size() or
      wsmr.gout.size() not_eq wsmr.gout_desc.size()) {
    errors.push_back(var_string("Mismatch sizes of gout (",
                                wsmr.gout.size(),
                                "), gout_type (",
                                wsmr.gout_type.size(),
                                "), gout_desc (",
                                wsmr.gout_desc.size(),
                                ") for ",
                                std::quoted(name)));
  }

  for (auto& type : wsmr.gout_type) {
    if (not valid_wsg(type)) {
      errors.push_back(var_string("Invalid gout_type for ",
                                  std::quoted(name),
                                  ": ",
                                  std::quoted(type)));
    }
  }

  for (auto& gvar : wsmr.gout) {
    if (wsv.find(gvar) not_eq wsv.end()) {
      errors.push_back(var_string("Invalid gout for ",
                                  std::quoted(name),
                                  ": ",
                                  std::quoted(gvar),
                                  " - it is already a WSV"));
    }
  }

  if (wsmr.gin.size() not_eq wsmr.gin_type.size() or
      wsmr.gin.size() not_eq wsmr.gin_value.size() or
      wsmr.gin.size() not_eq wsmr.gin_desc.size()) {
    errors.push_back(var_string("Mismatch sizes of gin (",
                                wsmr.gin.size(),
                                ") gin_type (",
                                wsmr.gin_type.size(),
                                "), gin_value (",
                                wsmr.gin_value.size(),
                                "), gin_desc (",
                                wsmr.gin_desc.size(),
                                ") for ",
                                std::quoted(name)));
  }

  for (auto& type : wsmr.gin_type) {
    if (not valid_wsg(type)) {
      errors.push_back(var_string(
          "Invalid gin_type for ", std::quoted(name), ": ", std::quoted(type)));
    }
  }

  for (auto& gvar : wsmr.gin) {
    if (wsv.find(gvar) not_eq wsv.end()) {
      errors.push_back(var_string("Invalid gin for ",
                                  std::quoted(name),
                                  ": ",
                                  std::quoted(gvar),
                                  " - it is already a WSV"));
    }
  }

  for (auto& var : wsmr.out) {
    if (wsv.find(var) == wsv.end()) {
      errors.push_back(var_string("Invalid out for ",
                                  std::quoted(name),
                                  ": ",
                                  std::quoted(var),
                                  " - it is not a WSV"));
    }
  }

  for (auto& var : wsmr.in) {
    if (wsv.find(var) == wsv.end()) {
      errors.push_back(var_string("Invalid in for ",
                                  std::quoted(name),
                                  ": ",
                                  std::quoted(var),
                                  " - it is not a WSV"));
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
  const auto wsm = internal_workspace_methods();

  ArrayOfString errors{};

  for (auto& [name, wsmr] : wsm) {
    if (wsmr.has_overloads()) {
      std::size_t i = 0;
      auto owsmr = make_overload(wsmr, i);
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
  const auto wsm = internal_workspace_methods();

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
} catch (std::exception& e) {
  throw std::runtime_error("Error in header():\n\n" + std::string(e.what()));
}

void call_function(std::ostream& os,
                   const std::string& name,
                   const WorkspaceMethodInternalRecord& wsmr) try {
  const auto wsm = internal_workspace_methods();
  const auto& wsv = workspace_variables();

  if (wsmr.has_any()) {
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
    if (wsmr.pass_workspace) {
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

  } else if (wsmr.has_overloads()) {
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
  const auto wsm = internal_workspace_methods();

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
  const auto wsm = internal_workspace_methods();

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
} catch (std::exception& e) {
  throw std::runtime_error("Error in implementation():\n\n" +
                           std::string(e.what()));
}

int main(int argc, char** argv) try {
  if (argc != 2) throw std::runtime_error("usage: make_auto_wsm <num_methods>");

  const int num_methods = std::stoi(argv[1]);
  std::ofstream head("auto_wsm.h");
  std::ofstream impl("auto_wsm.cpp");
  std::ofstream of("auto_wsmmeta.cpp");

  ARTS_USER_ERROR_IF(const auto err = scan_for_errors(); not err.empty(), err)

  header(head);
  implementation(impl, num_methods);
  meta_implementation(of);
} catch (std::exception& e) {
  std::cerr << "Cannot create the automatic methods with error:\n\n"
            << e.what() << '\n';
  return 1;
}
