#include "workspace_dimensions.h"

#include <algorithm>
#include <format>
#include <optional>
#include <stdexcept>
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

/*! The group of the elements of an array group, or empty if it is not one.
 *
 * An ArrayOf group only knows how to read its own length.  What is inside it is
 * shaped the way its element group is shaped.
 */
std::string element_group(const std::string& type) {
  constexpr std::string_view prefix = "ArrayOf";
  if (not type.starts_with(prefix)) return {};

  const auto elem = type.substr(prefix.size());

  // A nested array holds elements that have a length rather than a shape, and
  // the lengths are free to differ between them, as the paths of a path field
  // do.  There is nothing to compare them against.
  if (elem.starts_with(prefix)) return {};

  const auto& wsgs = internal_workspace_groups();
  const auto  ptr  = wsgs.find(elem);
  if (ptr == wsgs.end() or ptr->second.dim_size.empty() or ptr->second.map_type) return {};

  return elem;
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

/*! Checks that every element of an array is shaped the way its dimensions say.
 *
 * An array names the dimensions of its elements after its own length, but
 * reading them means visiting every element.  So all of the element dimensions
 * are checked together in a single pass, and the check can only say that the
 * elements do not agree, not which element or which dimension was wrong.
 *
 * Returns nothing unless every dimension inside the element is known, since a
 * partial shape cannot be compared against.
 */
std::optional<SizeCheck> inner_shape_check(const std::string&                          name,
                                           const WorkspaceVariableInternalRecord&      wsv,
                                           const std::unordered_map<std::string, Ref>& refs,
                                           const DimAccess&                            access) {
  if (wsv.inner_dims.empty()) return std::nullopt;

  const auto elem = element_group(wsv.type);
  if (elem.empty()) {
    throw std::runtime_error(std::format(
        R"(Workspace variable "{}" names inner dimensions but its group "{}" holds no elements that are shaped.

Inner dimensions say that every element of an array has the same shape, so they
need an ArrayOf whose element group declares dim_size.  They are not available for:

  - a nested array, e.g. ArrayOfArrayOf..., whose elements have a length each
    and are free to differ, as the paths of a path field do
  - a map, whose values would have to be reached through its keys
  - a group whose elements carry no sizes at all
)",
        name,
        wsv.type));
  }

  const auto& esize = internal_workspace_groups().at(elem).dim_size;
  if (wsv.inner_dims.size() != esize.size()) {
    throw std::runtime_error(std::format(
        R"(Workspace variable "{}" names {} inner dimension(s) but its elements of group "{}" have {}.

The whole shape of an element is compared at once, so every one of its dimensions must be named.
)",
        name,
        wsv.inner_dims.size(),
        elem,
        esize.size()));
  }

  std::vector<std::string> shape;
  std::vector<std::string> syms;
  std::vector<std::string> ref_vars;

  for (const auto& sym : wsv.inner_dims) {
    const auto ptr = refs.find(sym);

    // Nothing in scope says what this size is, and a partial shape cannot be compared
    if (sym.empty() or ptr == refs.end()) return std::nullopt;

    shape.push_back(ptr->second.expr);
    syms.push_back(sym);
    ref_vars.push_back(ptr->second.var);
  }

  keep_first(ref_vars);

  SizeCheck check;
  check.test = std::format("all_same_shape({{{}}}, {})", join(shape, ", ", ", "), access(name));
  check.constraint =
      std::format("every element of *{}* has the shape ({}), matching {}.",
                  name,
                  join(syms, ", ", ", "),
                  join(ref_vars | std::views::transform([](const std::string& v) { return std::format("*{}*", v); }) |
                           std::ranges::to<std::vector<std::string>>(),
                       ", ",
                       " and "));

  for (std::size_t i = 0; i < shape.size(); i++) { check.printables.emplace_back(syms[i], shape[i]); }

  return check;
}

}  // namespace

const std::unordered_map<std::string, WorkspaceDimensionRecord>& internal_workspace_dimensions() {
  static const auto out = internal_workspace_dimensions_creator();
  return out;
}

std::string dim_access_by_name(const std::string& name) { return name; }

void check_workspace_dimensions() {
  const auto& wsds = internal_workspace_dimensions();
  const auto& wsgs = internal_workspace_groups();

  for (const auto& [name, wsv] : internal_workspace_variables()) {
    for (const auto& sym : wsv.dims) {
      if (not wsds.contains(sym)) {
        throw std::runtime_error(
            std::format(R"(Workspace variable "{}" names the unknown dimension "{}".)", name, sym));
      }
    }

    const auto readable = readable_dims(wsv);
    if (readable != wsv.dims.size()) {
      throw std::runtime_error(std::format(
          R"(Workspace variable "{}" names {} dimension(s) but its group "{}" can read {}.

A dimension inside the elements of an array belongs in inner_dims, not in dims.
)",
          name,
          wsv.dims.size(),
          wsv.type,
          readable));
    }

    if (wsv.inner_dims.empty()) continue;

    for (const auto& sym : wsv.inner_dims) {
      if (not wsds.contains(sym)) {
        throw std::runtime_error(
            std::format(R"(Workspace variable "{}" names the unknown inner dimension "{}".)", name, sym));
      }
    }

    const auto elem = element_group(wsv.type);
    if (elem.empty()) {
      throw std::runtime_error(std::format(
          R"(Workspace variable "{}" names inner dimensions but its group "{}" holds no elements that are shaped.

Inner dimensions say that every element of an array has the same shape, so they
need an ArrayOf whose element group declares dim_size.  They are not available for:

  - a nested array, e.g. ArrayOfArrayOf..., whose elements have a length each
    and are free to differ, as the paths of a path field do
  - a map, whose values would have to be reached through its keys
  - a group whose elements carry no sizes at all
)",
          name,
          wsv.type));
    }

    const auto rank = wsgs.at(elem).dim_size.size();
    if (wsv.inner_dims.size() != rank) {
      throw std::runtime_error(std::format(
          R"(Workspace variable "{}" names {} inner dimension(s) but its elements of group "{}" have {}.

The whole shape of an element is compared at once, so every one of its dimensions must be named.
)",
          name,
          wsv.inner_dims.size(),
          elem,
          rank));
    }
  }
}

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

  for (const auto& name : targets) {
    const auto* wsv = find_wsv(name);
    if (wsv == nullptr) continue;

    const auto inner = inner_shape_check(name, *wsv, refs, access);
    if (inner) out.push_back(*inner);
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
  if (wsv == nullptr or (wsv->dims.empty() and wsv->inner_dims.empty())) return {};

  auto all = wsv->dims;
  all.insert(all.end(), wsv->inner_dims.begin(), wsv->inner_dims.end());

  std::string out = std::format("Shape: ({})", join(all, ", ", ", "));

  const auto& wsds = internal_workspace_dimensions();
  for (const auto& sym : all) {
    const auto ptr = wsds.find(sym);
    if (ptr == wsds.end()) continue;
    out += std::format(", {} is the {}", sym, ptr->second.desc);
  }

  return out + ".\n";
}
