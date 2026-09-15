#include "workspace_dimensions.h"

#include <algorithm>
#include <format>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <unordered_set>

#include "workspace_group_friends.h"
#include "workspace_groups.h"
#include "workspace_variables.h"

namespace {
std::unordered_map<std::string, WorkspaceDimensionRecord> internal_workspace_dimensions_creator() {
  std::unordered_map<std::string, WorkspaceDimensionRecord> wsd_data;

  wsd_data["NFREQ"]       = {.desc = "number of frequency points"};
  wsd_data["NPATH"]       = {.desc = "number of path points"};
  wsd_data["NTARGET"]     = {.desc = "number of Jacobian targets"};
  wsd_data["NSTATE"]      = {.desc = "size of the model state vector"};
  wsd_data["NLEGENDRE"]   = {.desc = "number of Legendre coefficients"};
  wsd_data["NQUADRATURE"] = {.desc = "number of quadrature points"};
  wsd_data["NALT"]        = {.desc = "number of altitude points"};
  wsd_data["NLAT"]        = {.desc = "number of latitude points"};
  wsd_data["NLON"]        = {.desc = "number of longitude points"};
  wsd_data["NZEN"]        = {.desc = "number of zenith angles"};
  wsd_data["NAZI"]        = {.desc = "number of azimuth angles"};
  wsd_data["NAUX"]        = {.desc = "number of auxiliary points"};
  wsd_data["NDEPTH"]      = {.desc = "number of subsurface depth points"};
  wsd_data["NMEAS"]       = {.desc = "number of measurements"};
  return wsd_data;
}

std::string subst(const std::string& expr, const std::string& name) {
  const auto pos = expr.find("{}");
  if (pos == std::string::npos) return expr;
  return expr.substr(0, pos) + name + expr.substr(pos + 2);
}

//! As subst, but for expressions that name the variable more than once
std::string subst_all(const std::string& expr, const std::string& name) {
  std::string out;
  std::size_t pos = 0;
  while (true) {
    const auto next = expr.find("{}", pos);
    if (next == std::string::npos) return out + expr.substr(pos);
    out += expr.substr(pos, next - pos) + name;
    pos  = next + 2;
  }
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

  auto elem = type.substr(prefix.size());

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

  //! Set when expr names a local that has to be declared before it is read
  std::string decl_name{};
  std::string decl_expr{};
};

//! The name of the local that holds a dimension read out of the elements of an array
std::string inner_ref_name(const std::string& sym) { return std::format("{}_DIM_", sym); }

/*! Attaches the declaration a reference needs, to the first check that reads it.
 *
 * The checks of one list are written in order into one scope, so declaring it on
 * the first reader puts it in scope for every later one.
 */
void declare_ref(SizeCheck& check, const Ref& ref, std::unordered_set<std::string>& declared) {
  if (ref.decl_name.empty() or not declared.insert(ref.decl_name).second) return;
  check.locals.emplace_back(ref.decl_name, ref.decl_expr);
}

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

/*! Fills the gaps by reading a dimension out of the elements of an array.
 *
 * Only for the dimensions that nothing names directly.  A variable that has the
 * dimension as its own is the better authority, because it says what the size is
 * even when it is zero, whereas an empty array says nothing about the shape of
 * the elements it does not have.
 *
 * Reading it means reaching into an element, so it is read once into a local and
 * every check of that dimension then shares it.  An array that is empty gives
 * zero, which is what the elements of an empty array agree on.
 */
void collect_inner_refs(std::unordered_map<std::string, Ref>& refs,
                        const std::vector<std::string>&       names,
                        const DimAccess&                      access) {
  for (const auto& name : names) {
    const auto* wsv = find_wsv(name);
    if (wsv == nullptr or wsv->inner_dims.empty()) continue;

    const auto elem = element_group(wsv->type);
    if (elem.empty()) continue;

    const auto& esize = internal_workspace_groups().at(elem).dim_size;
    if (wsv->inner_dims.size() != esize.size()) continue;

    for (std::size_t i = 0; i < wsv->inner_dims.size(); i++) {
      const auto& sym = wsv->inner_dims[i];
      if (sym.empty() or refs.contains(sym)) continue;

      const auto var  = access(name);
      const auto read = subst(esize[i], std::format("{}.front()", var));

      refs[sym] = {.var       = name,
                   .expr      = inner_ref_name(sym),
                   .decl_name = inner_ref_name(sym),
                   .decl_expr = std::format("{}.empty() ? Size{{0}} : static_cast<Size>({})", var, read)};
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
                                           std::unordered_set<std::string>&            declared,
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
  std::vector<std::string> needed;

  for (const auto& sym : wsv.inner_dims) {
    const auto ptr = refs.find(sym);

    // Nothing in scope says what this size is, and a partial shape cannot be compared
    if (sym.empty() or ptr == refs.end()) return std::nullopt;

    shape.push_back(ptr->second.expr);
    syms.push_back(sym);
    needed.push_back(sym);

    // An array compared against a size read from its own elements is being asked
    // whether its elements agree with each other, which is worth saying plainly
    if (ptr->second.var != name) ref_vars.push_back(ptr->second.var);
  }

  keep_first(ref_vars);

  SizeCheck check;
  // Qualified, because the generated code is compiled wherever the method lives
  // and cannot rely on the element group pulling matpack in by argument lookup
  check.test = std::format("matpack::all_same_shape({{{}}}, {})", join(shape, ", ", ", "), access(name));

  const auto shape_text = std::format("every element of *{}* has the shape ({})", name, join(syms, ", ", ", "));
  check.constraint = ref_vars.empty() ? shape_text + ", the same for all of them."
                                      : std::format("{}, matching {}.",
                                                    shape_text,
                                                    join(ref_vars | std::views::transform([](const std::string& v) {
                                                           return std::format("*{}*", v);
                                                         }) | std::ranges::to<std::vector<std::string>>(),
                                                         ", ",
                                                         " and "));

  for (const auto& sym : needed) declare_ref(check, refs.at(sym), declared);

  for (std::size_t i = 0; i < shape.size(); i++) { check.printables.emplace_back(syms[i], shape[i]); }

  return check;
}

}  // namespace

const std::unordered_map<std::string, WorkspaceDimensionRecord>& internal_workspace_dimensions() {
  static const auto out = internal_workspace_dimensions_creator();
  return out;
}

std::string dim_access_by_name(const std::string& name) { return name; }

namespace {
/*! Finds a group among the workspace groups or among the friends.
 *
 * The friends are documented next to the groups and may be promoted to groups
 * one day, so what is written about a group is written about them too.  Only the
 * documentation looks here: a workspace variable can only be of a real group, so
 * nothing that generates a check needs the friends.
 */
const WorkspaceGroupRecord* find_group(const std::string& group) {
  const auto& wsgs = internal_workspace_groups();
  if (const auto ptr = wsgs.find(group); ptr != wsgs.end()) return &ptr->second;

  const auto& friends = workspace_group_friends();
  if (const auto ptr = friends.find(group); ptr != friends.end()) return &ptr->second;

  return nullptr;
}
}  // namespace

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
  // A variable that has the dimension as its own is the better authority, so every
  // pass is asked for those before any of them is asked to reach into its elements
  // Each pass is asked for its own dimensions before it is asked to reach into
  // its elements, so an outer dimension wins over an inner one.  But a whole pass
  // is exhausted before the next is tried, because which variables are trusted to
  // carry the right size matters more: a buffer the method accumulates into names
  // its dimensions too, and it is what the checks are for rather than what they
  // are made against.
  std::unordered_map<std::string, Ref> refs;
  for (const auto& pass : references) {
    collect_refs(refs, pass, access);
    collect_inner_refs(refs, pass, access);
  }

  // Every size a check is made against is read once into a local named after its
  // dimension, so the generated code opens by saying where each size came from
  for (auto& [sym, ref] : refs) {
    if (not ref.decl_name.empty()) continue;
    ref.decl_name = inner_ref_name(sym);
    ref.decl_expr = std::format("static_cast<Size>({})", ref.expr);
    ref.expr      = ref.decl_name;
  }

  std::unordered_set<std::string> declared;
  std::vector<SizeCheck>          out;

  for (const auto& name : targets) {
    const auto* wsv = find_wsv(name);
    if (wsv == nullptr) continue;

    const auto n = readable_dims(*wsv);

    std::vector<std::string>                         tests;
    std::vector<std::string>                         syms;
    std::vector<std::string>                         ref_vars;
    std::vector<std::pair<std::string, std::string>> printables;
    std::vector<std::string>                         needed;
    std::vector<std::string>                         ref_exprs;
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

      tests.push_back(std::format("static_cast<Size>({}) == {}", expr, ptr->second.expr));
      syms.push_back(sym);
      ref_vars.push_back(ptr->second.var);
      printables.emplace_back(label, expr);
      needed.push_back(sym);
      ref_exprs.push_back(ptr->second.expr);

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

    // The dimensions are always compared one by one rather than as a whole shape.
    // Comparing whole shapes would mean calling matpack, and a group is free to
    // say how its dimensions are read without being a matpack tensor at all, as
    // the descriptor groups do.
    check.test = join(tests, " and ", " and ");

    // Naming the whole shape reads better than listing the dimensions, but only
    // says everything when every dimension of the variable takes part
    if (all_covered and tests.size() > 1) {
      check.constraint = std::format("*{}* has the shape ({}), matching {}.", name, join(syms, ", ", ", "), refs_text);

      // Reported whole rather than one dimension per line: what is wrong with a
      // shape is usually which axes were swapped, which only the whole shape shows
      const auto tuple = [](const std::vector<std::string>& v) {
        std::vector<std::string> fields(v.size(), "{}");
        return std::format(R"(std::format("[{}]", {}))", join(fields, ", ", ", "), join(v, ", ", ", "));
      };

      std::vector<std::string> actual;
      actual.reserve(n);
      for (std::size_t i = 0; i < n; i++) actual.push_back(subst(dim_size_template(*wsv, i), access(name)));

      printables.clear();
      printables.emplace_back("expected", tuple(ref_exprs));
      printables.emplace_back(std::format("{}.shape()", name), tuple(actual));
    } else {
      check.constraint =
          std::format("*{}* matches {} in {}.",
                      name,
                      refs_text,
                      join(syms | std::views::transform(dim_with_desc) | std::ranges::to<std::vector<std::string>>(),
                           ", ",
                           " and "));
    }

    check.printables = std::move(printables);
    for (const auto& sym : needed) declare_ref(check, refs.at(sym), declared);

    out.push_back(std::move(check));
  }

  for (const auto& name : targets) {
    const auto* wsv = find_wsv(name);
    if (wsv == nullptr) continue;

    const auto inner = inner_shape_check(name, *wsv, refs, declared, access);
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
  if (not wsmr.size_constraints) return {};

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

std::vector<SizeCheck> method_input_invariants(const WorkspaceMethodInternalRecord& wsmr, const DimAccess& access) {
  const auto& wsvs = internal_workspace_variables();
  const auto& wsgs = internal_workspace_groups();

  std::vector<SizeCheck> out;

  for (const auto& name : wsmr.in) {
    // A method that also writes the variable is building it, so it is allowed to
    // see it part-built.  Only what a method purely reads has to be usable.
    if (std::ranges::find(wsmr.out, name) != wsmr.out.end()) continue;

    const auto wsv = wsvs.find(name);
    if (wsv == wsvs.end()) continue;

    const auto wsg = wsgs.find(wsv->second.type);
    if (wsg == wsgs.end() or wsg->second.invariant.empty()) continue;

    const auto& rec = wsg->second;

    SizeCheck check;
    check.test       = subst_all(rec.invariant, access(name));
    check.constraint = std::format("On input, *{}* {}", name, rec.invariant_desc);
    for (const auto& expr : rec.invariant_printables) {
      check.printables.emplace_back(subst_all(expr, name), subst_all(expr, access(name)));
    }

    out.push_back(std::move(check));
  }

  return out;
}

std::vector<SizeCheck> method_output_size_checks(const WorkspaceMethodInternalRecord& wsmr, const DimAccess& access) {
  if (not wsmr.size_constraints) return {};

  // The outputs that were also inputs have been checked before the method ran
  std::vector<std::string> created;
  for (const auto& name : wsmr.out) {
    if (std::ranges::find(wsmr.in, name) == wsmr.in.end()) created.push_back(name);
  }

  auto out = size_checks({wsmr.in, wsmr.out}, created, access);
  for (auto& c : out) c.constraint = "On output, " + c.constraint;
  return out;
}

std::string size_check_code(const std::vector<SizeCheck>& checks, std::string_view indent) {
  if (checks.empty()) return {};

  // Scoped, so the locals a check declares cannot collide with those of the next
  // list written beside it, and so the accumulator lives no longer than it is used
  std::string out = std::format("{}{{\n", indent);

  for (const auto& check : checks) {
    for (const auto& [name, expr] : check.locals) {
      out += std::format("{}  const Size {} = {};\n", indent, name, expr);
    }
  }

  out += std::format("{}  std::string _err;\n", indent);

  for (const auto& check : checks) {
    // The message is a format string, so anything that looks like a field must be escaped
    std::string msg;
    for (const auto c : check.constraint) {
      if (c == '{' or c == '}') msg.push_back(c);
      msg.push_back(c);
    }

    std::size_t width = 0;
    for (const auto& [label, expr] : check.printables) width = std::max(width, label.size());

    out += std::format("{}  if (not ({}))\n{}    _err += std::format(R\"ERR(\n{}\n", indent, check.test, indent, msg);

    for (const auto& [label, expr] : check.printables) {
      out += std::format("\n{}:{} {{}}\n", label, std::string(width - label.size(), ' '));
    }

    out += ")ERR\"";
    for (const auto& [label, expr] : check.printables) out += std::format(", {}", expr);
    out += ");\n";
  }

  out += std::format("{}  if (not _err.empty()) throw std::runtime_error(_err);\n{}}}\n", indent, indent);

  return out;
}

std::string size_check_docs(const std::vector<SizeCheck>& checks) {
  if (checks.empty()) return {};

  std::string out = std::format("\n.. rubric:: Constraint{}\n\n", checks.size() > 1 ? "s" : "");
  for (const auto& c : checks) out += std::format("#. {}\n", c.constraint);
  return out + "\n";
}

std::string variable_dimension_docs(const std::string& name) {
  /* Links to the dimension page that gen_overview_list.py writes.  The label it
   * puts on each dimension is spelled the same way there, so the two have to be
   * changed together. */
  const auto dimension_link = [](const std::string& sym) {
    if (not internal_workspace_dimensions().contains(sym)) return sym;
    return std::format(":ref:`{} <wsd-{}>`", sym, sym);
  };

  const auto* wsv = find_wsv(name);
  if (wsv == nullptr or (wsv->dims.empty() and wsv->inner_dims.empty())) return {};

  // An array of equally shaped elements is as deep as its own length plus their
  // shape, so the effective shape is what an index into it has to provide
  auto all = wsv->dims;
  all.insert(all.end(), wsv->inner_dims.begin(), wsv->inner_dims.end());

  return std::format(
      "\n.. rubric:: Effective shape\n\n[{}]\n",
      join(all | std::views::transform(dimension_link) | std::ranges::to<std::vector<std::string>>(), ", ", ", "));
}

std::string group_invariant_docs(const std::string& group) {
  const auto* rec = find_group(group);
  if (rec == nullptr or rec->invariant.empty()) return {};

  return std::format("\n.. rubric:: Invariant\n\nA variable ``x`` of this group {}\n",
                     subst_all(rec->invariant_desc, "x"));
}

std::string group_dimension_docs(const std::string& group) {
  const auto* rec = find_group(group);
  if (rec == nullptr or rec->dim_size.empty()) return {};

  const auto& sizes = rec->dim_size;

  // Written against a variable called "x", since the expressions are templates
  // that a variable's own name is substituted into
  auto reads = sizes |
               std::views::transform([](const std::string& expr) { return std::format("``{}``", subst(expr, "x")); }) |
               std::ranges::to<std::vector<std::string>>();

  return std::format("\n.. rubric:: Sizes\n\nA variable ``x`` of this group may name {} dimension{}, read as {}.\n",
                     sizes.size(),
                     sizes.size() > 1 ? "s" : "",
                     join(reads, ", ", " and "));
}
