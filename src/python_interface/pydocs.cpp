#include "pydocs.h"

#include <string>

#include "auto_wsm.h"
#include "auto_wsv.h"
#include "compare.h"
#include "debug.h"

namespace Python {
String group_generics_inout(const String& group) try {
  const auto& wsms = internal_workspace_methods();
  const auto pred  = Cmp::eq(group);

  std::pair<std::vector<String>, std::vector<String>> outdocs;
  for (auto& [name, wsm] : wsms) {
    if (stdr::any_of(wsm.gout_type, pred)) outdocs.first.emplace_back(name);
    if (stdr::any_of(wsm.gin_type, pred)) outdocs.second.emplace_back(name);
  }

  std::sort(outdocs.first.begin(), outdocs.first.end(), str_compare_nocase);
  std::sort(outdocs.second.begin(), outdocs.second.end(), str_compare_nocase);

  String out;

  if (outdocs.first.size()) {
    out += std::format(R"(
.. rubric:: Workspace methods that can generate {}

.. hlist::
    :columns: {}
)",
                       group,
                       hlist_num_cols(outdocs.first));
    for (auto& m : outdocs.first)
      out += std::format("\n    *  :func:`~pyarts3.workspace.Workspace.{}`", m);
  }
  out += '\n';

  if (outdocs.second.size()) {
    out += std::format(R"(
.. rubric:: Workspace methods that require {}

.. hlist::
    :columns: {}
)",
                       group,
                       hlist_num_cols(outdocs.second));
    for (auto& m : outdocs.second)
      out += std::format("\n    *  :func:`~pyarts3.workspace.Workspace.{}`", m);
  }
  out += '\n';

  return out;
} catch (const std::exception& e) {
  return std::format("Error in group_generics_inout: {}",
                     std::string_view(e.what()));
}

String group_workspace_types(const String& group) try {
  const auto& wsvs = workspace_variables();

  std::vector<String> vars;
  for (auto& [name, wsv] : wsvs) {
    if (group == wsv.type) vars.emplace_back(name);
  }

  std::sort(vars.begin(), vars.end(), str_compare_nocase);

  String out;
  if (vars.size()) {
    out += std::format(R"(
.. rubric:: Workspace variables of type {}

.. hlist::
    :columns: {}
)",
                       group,
                       hlist_num_cols(vars));
    for (auto& m : vars)
      out += std::format("\n    *  :attr:`~pyarts3.workspace.Workspace.{}`", m);
  }

  return out + "\n";
} catch (const std::exception& e) {
  return std::format("Error in group_workspace_types: {}",
                     std::string_view(e.what()));
}
}  // namespace Python
