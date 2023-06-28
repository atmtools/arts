#include <mystring.h>

#include <algorithm>
#include <iomanip>
#include <iterator>
#include <stdexcept>
#include <vector>

#include "debug.h"
#include "global_data.h"
#include "nonstd.h"
#include "workspace.h"
#include "workspace_global_data.h"

String as_pyarts(const String& x) {
  static bool once = false;
  if (not once) {
    define_wsv_groups();
    define_wsv_data();
    define_wsv_map();
    define_md_data_raw();
    define_md_raw_map();
    expand_md_data_raw_to_md_data();
    define_md_map();
    define_agenda_data();
    define_agenda_map();
    once = true;
  }

  const auto find = [&](auto& map) { return map.find(x) not_eq map.end(); };

  if (find(global_data::WsvGroupMap))
    return var_string(":class:`~pyarts.arts.", x, '`');
  if (find(global_data::MdRawMap))
    return var_string(":func:`~pyarts.workspace.Workspace.", x, '`');
  if (find(global_data::WsvMap))
    return var_string(":attr:`~pyarts.workspace.Workspace.", x, '`');

  throw std::invalid_argument(
      var_string(std::quoted(x),
                 " is not a valid group, method, or workspace variable.\n"
                 "If it is a GIN or GOUT, consider representing it as \"``",
                 x,
                 "``\" instead?\nIf it is an old or deleted method or "
                 "variable or group, please remove it from the documentation!\n"));
}

String unwrap_stars(String x) {
  const auto find = [&](auto p) { return std::find(p + 1, x.end(), '*'); };
  const auto idx = [&](auto p) { return std::distance(x.begin(), p); };
  const auto mv = [&](auto p) {
    while (p not_eq x.end() and not nonstd::isspace(*p) and (*p) not_eq '*')
      p++;
    return std::min(p, x.end());
  };

  std::vector<Index> start_pos;

  auto ptr = find(x.begin());
  while (ptr not_eq x.end()) {
    const Index pos = idx(ptr);
    ptr = mv(ptr + 1);
    if (ptr not_eq x.end() and *(ptr) == '*') start_pos.push_back(pos);
    ptr = find(ptr + 1);
  }

  while (not start_pos.empty()) {
    const Index pos = start_pos.back();
    start_pos.pop_back();
    const auto start_ptr = x.begin() + pos;
    const auto end_ptr = 1 + mv(start_ptr + 1);
    const auto replace_name =
        x.substr(pos, std::distance(start_ptr, end_ptr));
    x.replace(start_ptr,
              end_ptr,
              as_pyarts(replace_name.substr(1, replace_name.size() - 2)));
  }

  return x;
}

String add_type(String x, const String& type) {
  const auto first_newline = std::find(x.begin(), x.end(), '\n');

  x.replace(first_newline,
            first_newline,
            var_string("\n\nThis workspace variable holds the group: :class:`~pyarts.arts.", type, "`\n\n"));
  return x;
}