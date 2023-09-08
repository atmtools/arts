#include "pydocs.h"

#include <string>

#include "auto_wsm.h"
#include "auto_wsv.h"
#include "compare.h"
#include "debug.h"

namespace Python {
String group_generics_inout(const String& group) {
  const auto wsms = internal_workspace_methods();

  std::pair<std::vector<String>, std::vector<String>> outdocs;
  for (auto& [name, wsm] : wsms) {
    if (std::ranges::any_of(wsm.gout_type, Cmp::eq(group)))
      outdocs.first.emplace_back(name);
    if (std::ranges::any_of(wsm.gin_type, Cmp::eq(group)))
      outdocs.second.emplace_back(name);
  }

  std::sort(outdocs.first.begin(), outdocs.first.end(), str_compare_nocase);
  std::sort(outdocs.second.begin(), outdocs.second.end(), str_compare_nocase);

  String out;

  if (outdocs.first.size()) {
    out += var_string("\nWorkspace methods that can generate ",
                      group,
                      "\n",
                      String(36 + group.size(), '-'),
                      "\n\n.. hlist::\n    :columns: ",
                      hlist_num_cols(outdocs.first),
                      "\n");
    for (auto& m : outdocs.first)
      out += var_string("\n    * :func:`~pyarts.workspace.Workspace.", m, '`');
  }
  out += '\n';

  if (outdocs.second.size()) {
    out += var_string("\nWorkspace methods that require ",
                      group,
                      "\n",
                      String(31 + group.size(), '-'),
                      "\n\n.. hlist::\n    :columns: ",
                      hlist_num_cols(outdocs.second),
                      "\n");
    for (auto& m : outdocs.second)
      out += var_string("\n    * :func:`~pyarts.workspace.Workspace.", m, '`');
  }
  out += '\n';

  return out;
}

String group_workspace_types(const String& group) {
  const auto& wsvs = workspace_variables();

  std::vector<String> vars;
  for (auto& [name, wsv] : wsvs) {
    if (group == wsv.type) vars.emplace_back(name);
  }

  std::sort(vars.begin(), vars.end(), str_compare_nocase);

  String out;
  if (vars.size()) {
    out += var_string("\nWorkspace variables of type ",
                      group,
                      "\n",
                      String(28 + group.size(), '-'),
                      "\n\n.. hlist::\n    :columns: ",
                      hlist_num_cols(vars),
                      "\n");
    for (auto& m : vars)
      out += var_string("\n    * :attr:`~pyarts.workspace.Workspace.", m, '`');
  }

  return out + "\n";
}
}  // namespace Python
