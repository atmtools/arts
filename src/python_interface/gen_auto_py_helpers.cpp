#include <mystring.h>

#include <algorithm>
#include <iomanip>
#include <iterator>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "compare.h"
#include "debug.h"
#include "global_data.h"
#include "matpack_constexpr.h"
#include "nonstd.h"
#include "python_interface/pydocs.h"
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

  const auto found_in = [&](auto& map) { return map.find(x) not_eq map.end(); };

  if (found_in(global_data::WsvGroupMap))
    return var_string(":class:`~pyarts.arts.", x, '`');
  if (found_in(global_data::MdRawMap))
    return var_string(":func:`~pyarts.workspace.Workspace.", x, '`');
  if (found_in(global_data::WsvMap))
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
  const auto find = [&](auto p) { p = std::min(p, x.end()); return std::find(p, x.end(), '*'); };
  const auto idx = [&](auto p) { return std::distance(x.begin(), p); };
  const auto mv = [&](auto p) {
    do 
      p++;
    while (p < x.end() and not nonstd::isspace(*p) and (*p) not_eq '*');
    return std::min(p, x.end());
  };

  std::vector<Index> start_pos;

  auto ptr = find(x.begin());
  while (ptr < x.end()) {
    const Index pos = idx(ptr);
    ptr = mv(ptr);
    if (ptr < x.end() and *(ptr) == '*') start_pos.push_back(pos);
    ptr = find(ptr + 1);
  }

  while (not start_pos.empty()) {
    const Index pos = start_pos.back();
    start_pos.pop_back();
    const auto start_ptr = x.begin() + pos;
    const auto end_ptr = 1 + mv(start_ptr);
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
            var_string("\n\nThis workspace variable holds the group: :class:`~pyarts.arts.", type, "`\n"));
  while (x.ends_with("\n\n")) x.pop_back();  // Renove extra spaces
  return x;
}

String get_agenda_io(const String & x) {
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

  String out{R"(
Parameters
----------
)"};

  struct AgendaIO {
    bool in;
    bool out;
    String group;
    String name;
  };

  auto& ag = global_data::agenda_data[global_data::AgendaMap.at(x)];
  auto& out_ind = ag.Out();
  auto& in_ind = ag.In();

  const auto output = [&](const Index ind) {
    return std::any_of(out_ind.begin(), out_ind.end(), Cmp::eq(ind));
  };

  const auto input = [&](const Index ind) {
    return std::any_of(in_ind.begin(), in_ind.end(), Cmp::eq(ind));
  };

  std::vector<AgendaIO> writer;
  for (auto& var : out_ind) {
    writer.emplace_back(
        input(var),
        true,
        global_data::wsv_groups[global_data::wsv_data[var].Group()].name,
        global_data::wsv_data[var].Name());
  }

  for (auto& var : in_ind) {
    if (not output(var))
      writer.emplace_back(
          true,
          false,
          global_data::wsv_groups[global_data::wsv_data[var].Group()].name,
          global_data::wsv_data[var].Name());
  }

  constexpr matpack::matpack_constant_data<std::string_view, 2, 2> inout{
      "(ERROR)", "(OUT)", "(IN)", "(INOUT)"};
  for (auto& var : writer) {
    out += var_string(var.name,
                      " : ~pyarts.arts.",
                      var.group,
                      '\n',
                      "    See :attr:`~pyarts.workspace.Workspace.",
                      var.name,
                      "` for more details ", inout[var.in][var.out], '\n');
  }

  return out;
}

String until_first_newline(const String& x) {
  std::string out{x.begin(), std::find(x.begin(), x.end(), '\n')};
  if (not out.ends_with('.')) out.push_back('.');
  return out;
}

String short_doc(const String& x) {
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

  const auto found_in = [&](auto& map) { return map.find(x) not_eq map.end(); };

  if (found_in(global_data::WsvGroupMap))
    return until_first_newline(global_data::wsv_groups[global_data::WsvGroupMap.at(x)].desc);
  if (found_in(global_data::MdRawMap))
    return until_first_newline(global_data::md_data_raw[global_data::MdRawMap.at(x)].Description());
  if (found_in(global_data::WsvMap))
    return until_first_newline(global_data::wsv_data[global_data::WsvMap.at(x)].Description());

  throw std::invalid_argument(
      var_string(std::quoted(x),
                 " is not a valid group, method, or workspace variable.\n"
                 "If it is a GIN or GOUT, consider representing it as \"``",
                 x,
                 "``\" instead?\nIf it is an old or deleted method or "
                 "variable or group, please remove it from the documentation!\n"));
}
