#include "workspace_dimensions.h"

#include <algorithm>
#include <format>

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

//! How many dimensions of this variable can actually be read
std::size_t readable_dims(const WorkspaceVariableInternalRecord& wsv) {
  const auto& wsgs = internal_workspace_groups();

  const auto ptr = wsgs.find(wsv.type);
  if (ptr == wsgs.end()) return 0;

  return std::min(wsv.dims.size(), std::max(ptr->second.dim_size.size(), wsv.dim_size.size()));
}

//! The expression that reads dimension i of the named variable
std::string dim_expr(const std::string& name, const WorkspaceVariableInternalRecord& wsv, std::size_t i) {
  if (i < wsv.dim_size.size() and not wsv.dim_size[i].empty()) return subst(wsv.dim_size[i], name);

  const auto& wsgs = internal_workspace_groups();
  const auto  ptr  = wsgs.find(wsv.type);
  if (ptr == wsgs.end() or i >= ptr->second.dim_size.size()) return {};

  return subst(ptr->second.dim_size[i], name);
}

const WorkspaceVariableInternalRecord* find_wsv(const std::string& name) {
  const auto& wsvs = internal_workspace_variables();
  const auto  ptr  = wsvs.find(name);
  return ptr == wsvs.end() ? nullptr : &ptr->second;
}

struct Ref {
  std::string var;
  std::string expr;
};

//! Records the first variable of the list that can be read for each of its dimensions
void collect_refs(std::unordered_map<std::string, Ref>& refs, const std::vector<std::string>& names) {
  for (const auto& name : names) {
    const auto* wsv = find_wsv(name);
    if (wsv == nullptr) continue;

    const auto n = readable_dims(*wsv);
    for (std::size_t i = 0; i < n; i++) {
      const auto& sym = wsv->dims[i];
      if (sym.empty() or refs.contains(sym)) continue;

      auto expr = dim_expr(name, *wsv, i);
      if (expr.empty()) continue;

      refs[sym] = {.var = name, .expr = std::move(expr)};
    }
  }
}

std::string dim_with_desc(const std::string& sym) {
  const auto& wsds = internal_workspace_dimensions();
  const auto  ptr  = wsds.find(sym);
  if (ptr == wsds.end()) return sym;
  return std::format("{} ({})", sym, ptr->second.desc);
}
}  // namespace

const std::unordered_map<std::string, WorkspaceDimensionRecord>& internal_workspace_dimensions() {
  static const auto out = internal_workspace_dimensions_creator();
  return out;
}

std::vector<StringVectorAgendaHelper> agenda_output_size_checks(const WorkspaceAgendaInternalRecord& ag) {
  if (not ag.output_constraints) return {};

  // Prefer checking against what the agenda was given over what it produced
  std::vector<std::string> pure_input;
  for (const auto& name : ag.input) {
    if (std::ranges::find(ag.output, name) == ag.output.end()) pure_input.push_back(name);
  }

  std::unordered_map<std::string, Ref> refs;
  collect_refs(refs, pure_input);
  collect_refs(refs, ag.input);
  collect_refs(refs, ag.output);

  std::vector<StringVectorAgendaHelper> out;

  for (const auto& name : ag.output) {
    const auto* wsv = find_wsv(name);
    if (wsv == nullptr) continue;

    const auto n = readable_dims(*wsv);

    std::vector<std::string> tests;
    std::vector<std::string> syms;
    std::vector<std::string> ref_exprs;
    std::vector<std::string> ref_vars;
    std::vector<std::string> printables;
    bool                     all_covered = n > 0;

    for (std::size_t i = 0; i < n; i++) {
      const auto& sym = wsv->dims[i];
      const auto  ptr = refs.find(sym);

      // Nothing to check against, or this variable is itself what others are checked against
      if (sym.empty() or ptr == refs.end() or ptr->second.var == name) {
        all_covered = false;
        continue;
      }

      const auto expr = dim_expr(name, *wsv, i);
      if (expr.empty()) {
        all_covered = false;
        continue;
      }

      tests.push_back(std::format("static_cast<Size>({}) == static_cast<Size>({})", expr, ptr->second.expr));
      syms.push_back(sym);
      ref_exprs.push_back(ptr->second.expr);
      if (std::ranges::find(ref_vars, ptr->second.var) == ref_vars.end()) ref_vars.push_back(ptr->second.var);
      printables.push_back(expr);
      printables.push_back(ptr->second.expr);
    }

    if (tests.empty()) continue;

    std::string refs_text;
    for (std::size_t i = 0; i < ref_vars.size(); i++) {
      if (i != 0) refs_text += (i + 1 == ref_vars.size()) ? " and " : ", ";
      refs_text += std::format("*{}*", ref_vars[i]);
    }

    std::string test;
    std::string constraint;

    // Whole-shape checks read better than one check per dimension, but only work
    // when every dimension of the variable takes part
    if (all_covered and tests.size() > 1) {
      std::string shape;
      for (std::size_t i = 0; i < ref_exprs.size(); i++) {
        if (i != 0) shape += ", ";
        shape += ref_exprs[i];
      }
      test = std::format("same_shape({{{}}}, {})", shape, name);

      std::string sym_text;
      for (std::size_t i = 0; i < syms.size(); i++) {
        if (i != 0) sym_text += ", ";
        sym_text += syms[i];
      }
      constraint = std::format("On output, *{}* has the shape ({}), matching {}.", name, sym_text, refs_text);

      printables.push_back(std::format("{}.shape()", name));
    } else {
      for (std::size_t i = 0; i < tests.size(); i++) {
        if (i != 0) test += " and ";
        test += tests[i];
      }

      std::string sym_text;
      for (std::size_t i = 0; i < syms.size(); i++) {
        if (i != 0) sym_text += (i + 1 == syms.size()) ? " and " : ", ";
        sym_text += dim_with_desc(syms[i]);
      }
      constraint = std::format("On output, *{}* matches {} in {}.", name, refs_text, sym_text);
    }

    std::vector<std::string> unique_printables;
    for (auto& p : printables) {
      if (std::ranges::find(unique_printables, p) == unique_printables.end()) unique_printables.push_back(p);
    }
    printables = std::move(unique_printables);

    StringVectorAgendaHelper helper;
    helper.test       = std::move(test);
    helper.constraint = std::move(constraint);
    helper.printables = std::move(printables);
    out.push_back(std::move(helper));
  }

  return out;
}

std::string variable_dimension_docs(const std::string& name) {
  const auto* wsv = find_wsv(name);
  if (wsv == nullptr or wsv->dims.empty()) return {};

  std::string shape;
  for (std::size_t i = 0; i < wsv->dims.size(); i++) {
    if (i != 0) shape += ", ";
    shape += wsv->dims[i];
  }

  std::string out = std::format("Shape: ({})", shape);

  const auto& wsds = internal_workspace_dimensions();
  for (const auto& sym : wsv->dims) {
    const auto ptr = wsds.find(sym);
    if (ptr == wsds.end()) continue;
    out += std::format(", {} is the {}", sym, ptr->second.desc);
  }

  return out + ".\n";
}
