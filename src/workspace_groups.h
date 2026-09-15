#pragma once

#include <format_tags.h>

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct WorkspaceGroupRecord {
  std::string file;
  std::string desc;

  //! Set to true if python must treat this as a value
  bool value_type{false};

  //! Set to true if the type is a map of some type
  bool map_type{false};

  /*! How to read the size of each dimension of this group, outermost first.
   *
   * Each entry is a C++ expression where "{}" is replaced by the name of a
   * variable of this group, e.g. "{}.size()".  For tensor-like groups these
   * are the actual axes.  For descriptor types the entries are the implied
   * sizes the type carries, e.g. *JacobianTargets* provides both its target
   * count and its model state vector size.
   *
   * Leave empty for groups that carry no queryable sizes.  Sizes are only
   * checked for variables whose record names the dimensions, see
   * WorkspaceVariableInternalRecord::dims.
   */
  std::vector<std::string> dim_size{};

  /*! A C++ predicate that holds for every usable value of this group.
   *
   * "{}" is replaced by the name of a variable of this group, so
   * "not {}.bad_ellipsoid()" becomes "not surf_field.bad_ellipsoid()".  It is
   * verified wherever a variable of this group is given to a method, which is
   * where a user can supply one that is not usable.
   *
   * It runs on every such call, including inside the loops of an agenda, so it
   * must be cheap: reading a few sizes or comparing two shapes, never walking
   * the data.  Leave it empty for a group that is usable in any state it can be
   * constructed in, which is most of them.
   */
  std::string invariant{};

  //! What the invariant means, for the error message and the documentation
  std::string invariant_desc{};

  /*! What to report when the invariant fails, as expressions to read.
   *
   * Each is a C++ expression where "{}" is replaced by the name of a variable of
   * this group, e.g. "{}.ellipsoid".  The expression is also what labels the
   * value in the message, so it says where the number came from.
   */
  std::vector<std::string> invariant_printables{};
};

const std::unordered_map<std::string, WorkspaceGroupRecord>& internal_workspace_groups();

void add_arrays_of(std::unordered_map<std::string, WorkspaceGroupRecord>& wsg_data,
                   const std::vector<std::string>&                        types,
                   std::vector<std::string>                               extra_headers);

template <> struct std::formatter<WorkspaceGroupRecord> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const WorkspaceGroupRecord& wsg, FmtContext& ctx) const {
    return tags.format(ctx,
                       "WorkspaceGroupRecord{"sv,
                       "\n  .file="sv,
                       wsg.file,
                       "\n  .desc="sv,
                       wsg.desc,
                       "\n  .value_type="sv,
                       wsg.value_type ? "true"sv : "false"sv,
                       "\n  .map_type="sv,
                       wsg.map_type ? "true"sv : "false"sv,
                       "\n  .dim_size="sv,
                       wsg.dim_size,
                       "\n  .invariant="sv,
                       wsg.invariant,
                       "\n  .invariant_desc="sv,
                       wsg.invariant_desc,
                       "\n}"sv);
  }
};
