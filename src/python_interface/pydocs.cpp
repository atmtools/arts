#include "pydocs.h"

#include <string>

#include "compare.h"
#include "debug.h"
#include "global_data.h"
#include "workspace_global_data.h"

namespace Python {
String group_generics_inout(const String& group) {
  const Index gr = global_data::WsvGroupMap.at(group);

  std::pair<std::vector<String>, std::vector<String>> outdocs;
  for (auto& method : global_data::md_data_raw) {
    if (std::any_of(method.GOutType().cbegin(),
                    method.GOutType().cend(),
                    Cmp::eq(gr)))
      outdocs.first.push_back(method.Name());
    if (std::any_of(method.GInType().cbegin(),
                    method.GInType().cend(),
                    Cmp::eq(gr)))
      outdocs.second.push_back(method.Name());
  }

  String out;

  if (outdocs.first.size()) {
    out += var_string("\nGeneric methods that can generate ",
                      group,
                      "\n",
                      String(34 + group.size(), '-'),
                      "\n");
    for (auto& m : outdocs.first)
      out += var_string("\n\n    :func:`~pyarts.workspace.Workspace.", m, '`');
  }

  if (outdocs.second.size()) {
    out += var_string("\nGeneric methods that require ",
                      group,
                      "\n",
                      String(29 + group.size(), '-'),
                      "\n");
    for (auto& m : outdocs.second)
      out += var_string("\n\n    :func:`~pyarts.workspace.Workspace.", m, '`');
  }

  //! FIXME: Add Workspace Variables of type XYZ
  //! FIXME: Use RAW MD instead of generated MD

  return out;
}

String group_workspace_types(const String& group) {
  const Index gr = global_data::WsvGroupMap.at(group);

  std::vector<String> vars;
  for (auto& var: global_data::wsv_data) {
    if (gr == var.Group()) vars.push_back(var.Name());
  }

  String out;
  if (vars.size()) {
    out += var_string("\nWorkspace variables of type ",
                      group,
                      "\n",
                      String(28 + group.size(), '-'),
                      "\n");
    for (auto& m : vars)
      out += var_string("\n\n    :attr:`~pyarts.workspace.Workspace.", m, '`');
  }

  return out;
}
}  // namespace Python