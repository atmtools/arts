#include "workspace_dimensions.h"

#include <algorithm>
#include <format>
#include <ranges>

#include "workspace_groups.h"
#include "workspace_variables.h"

namespace {
std::unordered_map<std::string, WorkspaceDimensionRecord> internal_workspace_dimensions_creator() {
  std::unordered_map<std::string, WorkspaceDimensionRecord> wsd_data;

  wsd_data["NF"] = {.desc = "number of frequency points"};
  wsd_data["NP"] = {.desc = "number of path points"};
  wsd_data["NT"] = {.desc = "number of Jacobian targets"};
  wsd_data["NX"] = {.desc = "size of the model state vector"};
  wsd_data["NL"] = {.desc = "number of Legendre coefficients"};

  return wsd_data;
}

std::string subst(const std::string& expr, const std::string& name) {
  const auto pos = expr.find("{}");
  if (pos == std::string::npos) return expr;
  return expr.substr(0, pos) + name + expr.substr(pos + 2);
}

const WorkspaceVariableInternalRecord* find_wsv(const std::string& name) {
  const auto& wsvs = internal_workspace_variables();
  const auto  ptr  = wsvs.find(name);
  return ptr == wsvs.end() ? nullptr : &ptr->second;
}

//! The expression template that reads dimension i of a variable, or empty if it cannot be read
std::string dim_size_template(const WorkspaceVariableInternalRecord& wsv, std::size_t i) {
  if (i < wsv.dim_size.size() and not wsv.dim_size[i].empty()) return wsv.dim_size[i];

  const auto& wsgs = internal_workspace_groups();
  const auto  ptr  = wsgs.find(wsv.type);
  if (ptr == wsgs.end() or i >= ptr->second.dim_size.size()) return {};

  return ptr->second.dim_size[i];
}

//! How many dimensions of this variable can actually be read
std::size_t readable_dims(const WorkspaceVariableInternalRecord& wsv) {
  std::size_t n = 0;
  while (n < wsv.dims.size() and not dim_size_template(wsv, n).empty()) n++;
  return n;
}

struct Ref {
  std::string var;
  std::string expr;
};

//! Records the first variable of the list that can be read for each of its dimensions
void collect_refs(std::unordered_map<std::string, Ref>& refs,
                  const std::vector<std::string>&       names,
                  const DimAccess&                      access) {
  for (const auto& name : names) {
    const auto* wsv = find_wsv(name);
    if (wsv == nullptr) continue;

    const auto n = readable_dims(*wsv);
    for (std::size_t i = 0; i < n; i++) {
      const auto& sym = wsv->dims[i];
      if (sym.empty() or refs.contains(sym)) continue;

      refs[sym] = {.var = name, .expr = subst(dim_size_template(*wsv, i), access(name))};
    }
  }
}

std::string dim_with_desc(const std::string& sym) {
  const auto& wsds = internal_workspace_dimensions();
  const auto  ptr  = wsds.find(sym);
  if (ptr == wsds.end()) return sym;
  return std::format("{} ({})", sym, ptr->second.desc);
}

std::string join(const std::vector<std::string>& v, std::string_view sep, std::string_view last_sep) {
  std::string out;
  for (std::size_t i = 0; i < v.size(); i++) {
    if (i != 0) out += (i + 1 == v.size()) ? last_sep : sep;
    out += v[i];
  }
  return out;
}

//! Removes the repeats while keeping the order
template <typename T> void keep_first(std::vector<T>& v) {
  std::vector<T> out;
  for (auto& x : v) {
    if (std::ranges::find(out, x) == out.end()) out.push_back(std::move(x));
  }
  v = std::move(out);
}
}  // namespace

const std::unordered_map<std::string, WorkspaceDimensionRecord>& internal_workspace_dimensions() {
  static const auto out = internal_workspace_dimensions_creator();
  return out;
}

std::string dim_access_by_name(const std::string& name) { return name; }

std::vector<SizeCheck> size_checks(const std::vector<std::vector<std::string>>& references,
                                   const std::vector<std::string>&              targets,
                                   const DimAccess&                             access) {
  std::unordered_map<std::string, Ref> refs;
  for (const auto& pass : references) collect_refs(refs, pass, access);

  std::vector<SizeCheck> out;

  for (const auto& name : targets) {
    const auto* wsv = find_wsv(name);
    if (wsv == nullptr) continue;

    const auto n = readable_dims(*wsv);

    std::vector<std::string>                         tests;
    std::vector<std::string>                         syms;
    std::vector<std::string>                         ref_exprs;
    std::vector<std::string>                         ref_vars;
    std::vector<std::pair<std::string, std::string>> printables;
    bool                                             all_covered = n > 0;

    for (std::size_t i = 0; i < n; i++) {
      const auto& sym = wsv->dims[i];
      const auto  ptr = refs.find(sym);

      // Nothing to check against, or this variable is itself what others are checked against
      if (sym.empty() or ptr == refs.end() or ptr->second.var == name) {
        all_covered = false;
        continue;
      }

      const auto tmpl  = dim_size_template(*wsv, i);
      const auto expr  = subst(tmpl, access(name));
      const auto label = subst(tmpl, name);

      tests.push_back(std::format("static_cast<Size>({}) == static_cast<Size>({})", expr, ptr->second.expr));
      syms.push_back(sym);
      ref_exprs.push_back(ptr->second.expr);
      ref_vars.push_back(ptr->second.var);
      printables.emplace_back(label, expr);

      const auto* ref_wsv = find_wsv(ptr->second.var);
      if (ref_wsv != nullptr) {
        for (std::size_t j = 0; j < readable_dims(*ref_wsv); j++) {
          if (ref_wsv->dims[j] == sym) {
            printables.emplace_back(subst(dim_size_template(*ref_wsv, j), ptr->second.var), ptr->second.expr);
            break;
          }
        }
      }
    }

    if (tests.empty()) continue;

    keep_first(ref_vars);
    keep_first(printables);

    const auto refs_text = join(ref_vars | std::views::transform([](const std::string& v) {
                                  return std::format("*{}*", v);
                                }) | std::ranges::to<std::vector<std::string>>(),
                                ", ",
                                " and ");

    SizeCheck check;

    // Whole-shape checks read better than one check per dimension, but only work
    // when every dimension of the variable takes part
    if (all_covered and tests.size() > 1) {
      check.test       = std::format("same_shape({{{}}}, {})", join(ref_exprs, ", ", ", "), access(name));
      check.constraint = std::format("*{}* has the shape ({}), matching {}.", name, join(syms, ", ", ", "), refs_text);
      printables.emplace_back(std::format("{}.shape()", name), std::format("{}.shape()", access(name)));
    } else {
      check.test = join(tests, " and ", " and ");
      check.constraint =
          std::format("*{}* matches {} in {}.",
                      name,
                      refs_text,
                      join(syms | std::views::transform(dim_with_desc) | std::ranges::to<std::vector<std::string>>(),
                           ", ",
                           " and "));
    }

    check.printables = std::move(printables);
    out.push_back(std::move(check));
  }

  return out;
}

std::vector<SizeCheck> agenda_output_size_checks(const WorkspaceAgendaInternalRecord& ag) {
  if (not ag.output_constraints) return {};

  // Prefer checking against what the agenda was given over what it produced
  std::vector<std::string> pure_input;
  for (const auto& name : ag.input) {
    if (std::ranges::find(ag.output, name) == ag.output.end()) pure_input.push_back(name);
  }

  auto out = size_checks({pure_input, ag.input, ag.output}, ag.output);
  for (auto& c : out) c.constraint = "On output, " + c.constraint;
  return out;
}

std::vector<SizeCheck> method_input_size_checks(const WorkspaceMethodInternalRecord& wsmr, const DimAccess& access) {
  // What the method only reads says what the sizes are.  What it also writes is
  // a buffer it accumulates into, so that is what gets checked against them.
  std::vector<std::string> pure_input;
  for (const auto& name : wsmr.in) {
    if (std::ranges::find(wsmr.out, name) == wsmr.out.end()) pure_input.push_back(name);
  }

  auto out = size_checks({pure_input, wsmr.in}, wsmr.in, access);
  for (auto& c : out) c.constraint = "On input, " + c.constraint;
  return out;
}

std::vector<SizeCheck> method_output_size_checks(const WorkspaceMethodInternalRecord& wsmr, const DimAccess& access) {
  // The outputs that were also inputs have been checked before the method ran
  std::vector<std::string> created;
  for (const auto& name : wsmr.out) {
    if (std::ranges::find(wsmr.in, name) == wsmr.in.end()) created.push_back(name);
  }

  auto out = size_checks({wsmr.in, wsmr.out}, created, access);
  for (auto& c : out) c.constraint = "On output, " + c.constraint;
  return out;
}

std::string size_check_code(const SizeCheck& check, std::string_view indent) {
  // The message is a format string, so anything that looks like a field must be escaped
  std::string msg;
  for (const auto c : check.constraint) {
    if (c == '{' or c == '}') msg.push_back(c);
    msg.push_back(c);
  }

  std::size_t width = 0;
  for (const auto& [label, expr] : check.printables) width = std::max(width, label.size());

  std::string out = std::format("{}if (not ({}))\n{}  throw std::runtime_error(std::format(R\"ERR({}\n",
                                indent,
                                check.test,
                                indent,
                                msg);

  for (const auto& [label, expr] : check.printables) {
    out += std::format("\n{}:{} {{}}\n", label, std::string(width - label.size(), ' '));
  }

  out += ")ERR\"";
  for (const auto& [label, expr] : check.printables) out += std::format(", {}", expr);
  out += "));\n";

  return out;
}

std::string size_check_docs(const std::vector<SizeCheck>& checks) {
  if (checks.empty()) return {};

  std::string out = std::format("\n.. rubric:: Constraint{}\n\n", checks.size() > 1 ? "s" : "");
  for (const auto& c : checks) out += std::format("#. {}\n", c.constraint);
  return out + "\n";
}

std::string variable_dimension_docs(const std::string& name) {
  const auto* wsv = find_wsv(name);
  if (wsv == nullptr or wsv->dims.empty()) return {};

  std::string out = std::format("Shape: ({})", join(wsv->dims, ", ", ", "));

  const auto& wsds = internal_workspace_dimensions();
  for (const auto& sym : wsv->dims) {
    const auto ptr = wsds.find(sym);
    if (ptr == wsds.end()) continue;
    out += std::format(", {} is the {}", sym, ptr->second.desc);
  }

  return out + ".\n";
}
