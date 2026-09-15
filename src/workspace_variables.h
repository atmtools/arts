#pragma once

#include <string>
#include <unordered_map>
#include <vector>

struct WorkspaceVariableInternalRecord {
  std::string desc;
  std::string type;
  std::string default_value{""};

  /*! Names the dimensions of this variable, outermost first.
   *
   * Each entry is a symbol from internal_workspace_dimensions().  The entries
   * line up with the dim_size of this variable's group, so a variable of a
   * group with N dimension expressions may name up to N dimensions.
   *
   * Variables that share a symbol must agree on that size.  This is what
   * agenda output size verification is generated from, and what the generated
   * documentation reports as the shape of the variable.
   */
  std::vector<std::string> dims{};

  /*! Names the dimensions inside each element of an array variable.
   *
   * An array grows along the dimension named in dims, e.g. a path gains points
   * to keep the optical depth per step small, while every element it holds is
   * shaped the same.  Those element dimensions are named here, in the order the
   * element's group reads them.
   *
   * Checking these means visiting every element, so they are verified together
   * in a single pass.  Such a check can only report that the elements disagree,
   * not which element or which dimension was wrong.
   */
  std::vector<std::string> inner_dims{};

  /*! Overrides the group's dim_size expressions for this variable.
   *
   * Only for variables whose size is not read the way its group reads it,
   * e.g. *legendre_degree* counts one less than the dimension it names.
   * Entries that are left empty use the group's expression.
   */
  std::vector<std::string> dim_size{};
};

const std::unordered_map<std::string, WorkspaceVariableInternalRecord>& internal_workspace_variables();

std::string_view any_is_typename(const std::string& type);
