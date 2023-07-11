#include <mystring.h>

#include <algorithm>
#include <iomanip>
#include <iterator>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "array.h"
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

  throw std::invalid_argument(var_string(
      std::quoted(x),
      " is not a valid group, method, or workspace variable.\n"
      "If it is a GIN or GOUT, consider representing it as \"``",
      x,
      "``\" instead?\nIf it is an old or deleted method or "
      "variable or group, please remove it from the documentation!\n"));
}

uint32_t hlist_num_cols(const std::vector<String>& v) {
  return v.size() < 5 ? 1 : 2;
};

String unwrap_stars(String x) {
  const auto find = [&](auto p) {
    p = std::min(p, x.end());
    return std::find(p, x.end(), '*');
  };
  const auto idx = [&](auto p) { return std::distance(x.begin(), p); };
  const auto mv = [&](auto p) {
    do p++;
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
    const auto replace_name = x.substr(pos, std::distance(start_ptr, end_ptr));
    x.replace(start_ptr,
              end_ptr,
              as_pyarts(replace_name.substr(1, replace_name.size() - 2)));
  }

  return x;
}

String get_agenda_io(const String& x) {
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
    writer.emplace_back(AgendaIO{
        input(var),
        true,
        global_data::wsv_groups[global_data::wsv_data[var].Group()].name,
        global_data::wsv_data[var].Name()});
  }

  for (auto& var : in_ind) {
    if (not output(var))
      writer.emplace_back(AgendaIO{
          true,
          false,
          global_data::wsv_groups[global_data::wsv_data[var].Group()].name,
          global_data::wsv_data[var].Name()});
  }

  constexpr matpack::matpack_constant_data<std::string_view, 2, 2> inout{
      "[ERROR]", "[OUT]", "[IN]", "[INOUT]"};
  for (auto& var : writer) {
    out += var_string(var.name,
                      " : ~pyarts.arts.",
                      var.group,
                      '\n',
                      "    ", unwrap_stars(short_doc(var.name)),
                      " See :attr:`~pyarts.workspace.Workspace.",
                      var.name,
                      "` **",
                      inout[var.in][var.out],
                      "**\n");
  }

  return writer.size() ? out : "";
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
    return until_first_newline(
        global_data::wsv_groups[global_data::WsvGroupMap.at(x)].desc);
  if (found_in(global_data::MdRawMap))
    return until_first_newline(
        global_data::md_data_raw[global_data::MdRawMap.at(x)].Description());
  if (found_in(global_data::WsvMap))
    return until_first_newline(
        global_data::wsv_data[global_data::WsvMap.at(x)].Description());

  throw std::invalid_argument(var_string(
      std::quoted(x),
      " is not a valid group, method, or workspace variable.\n"
      "If it is a GIN or GOUT, consider representing it as \"``",
      x,
      "``\" instead?\nIf it is an old or deleted method or "
      "variable or group, please remove it from the documentation!\n"));
}

void remove_trailing(String& x, char n) {
  while (x.ends_with(n)) x.pop_back();
}

String compose_generic_groups(Index grp, const ArrayOfIndex& inds) {
  if (inds.nelem() == 0)
    return var_string("~pyarts.arts.", global_data::wsv_groups[grp].name);

  String out;
  for (auto i : inds)
    out +=
        var_string("~pyarts.arts.", global_data::wsv_groups[i].name, " or ");

  // Remove trailing " or "
  out.pop_back();
  out.pop_back();
  out.pop_back();
  out.pop_back();
  return out;
}

String to_defval_str(const String& x, const String& group) {
  std::string out = var_string(x);

  auto pos = out.find('\n');
  while (pos not_eq out.npos) {
    out.replace(pos, 1, " ");
    pos = out.find('\n');
  }
  
  while (out.front() == ' ') out.erase(out.begin());
  while (out.back() == ' ') out.pop_back();

  if (group == "String" and out.front() not_eq '"' and out.back() not_eq '"') {
    return var_string('"', out, '"');
  }

  if (out.size() == 0) {
    if (group.starts_with("Array") or group == "Vector" or group == "Matrix" or
        group == "Tensor3" or group == "Tensor4" or group == "Tensor5" or
        group == "Tensor6" or group == "Tensor7")
      return "[]";

    if (group == "PredefinedModelData")
      return "pyarts.arts.PredefinedModelData()";

    if (group == "Numeric" or group == "Index") return "0";

    throw std::runtime_error(var_string(
        "Cannot interpret empty default value for ",
        group,
        " to a good python type. Please add one.\n"
        "Ensure that the default value is usable in the constructor inside python!\n"
        "Even better, please create a formatter for all Arts's types that print them\n"
        "in agood pythonesque way.\n\nTHIS IS A DEVELOPER ERROR!\n"));
  }

  return out;
}

String method_docs(const String& name) {
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

  //! WARNING: Raw method
  const auto& method = global_data::md_data_raw[global_data::MdRawMap.at(name)];
  String out;

  const auto is_output = [&](auto ind) {
    return std::any_of(method.Out().begin(), method.Out().end(), Cmp::eq(ind));
  };
  const auto is_input = [&](auto ind) {
    return std::any_of(method.In().begin(), method.In().end(), Cmp::eq(ind));
  };
  const auto is_ginput = [&](auto& ind) {
    return std::any_of(method.GIn().begin(), method.GIn().end(), Cmp::eq(ind));
  };
  const auto is_goutput = [&](auto& ind) {
    return std::any_of(
        method.GOut().begin(), method.GOut().end(), Cmp::eq(ind));
  };
  const auto fix = [&] {
    remove_trailing(out, '\n');
    out += "\n";
  };

  out = "--(";
  out += unwrap_stars(method.Description());
  fix();

  out += "\nAuthor(s):";
  for (auto& author : method.Authors()) out += " " + author + ", ";
  out.pop_back();  // Remove ' '
  out.pop_back();  // Remove ','
  fix();

  out += "\nParameters\n----------";
  for (auto i : method.Out()) {
    const String io = is_input(i) ? "INOUT" : "OUT";
    const auto& wsv = global_data::wsv_data[i];
    const auto& varname = wsv.Name();
    const auto& grpname = global_data::wsv_groups[wsv.Group()];
    out += var_string('\n',
                      varname,
                      " : ~pyarts.arts.",
                      grpname,
                      ", optional\n    ",
                      unwrap_stars(short_doc(varname)),
                      " See :attr:`~pyarts.workspace.Workspace.",
                      varname,
                      "`, defaults to ``self.",
                      varname,
                      "`` **[",
                      io,
                      "]**");
  }

  for (Index i = 0; i < method.GOut().nelem(); i++) {
    const auto& varname = method.GOut()[i];
    const String io = is_ginput(varname) ? "INOUT" : "OUT";
    const auto& grpname =
        compose_generic_groups(method.GOutType()[i], method.GOutSpecType()[i]);
    out += var_string('\n',
                      varname,
                      " : ",
                      grpname,
                      "\n    ",
                      unwrap_stars(until_first_newline(method.GOutDescription()[i])),
                      " **[",
                      io,
                      "]**");
  }

  for (auto i : method.In()) {
    if (is_output(i)) continue;

    const auto& wsv = global_data::wsv_data[i];
    const auto& varname = wsv.Name();
    const auto& grpname = global_data::wsv_groups[wsv.Group()];
    out += var_string('\n',
                      varname,
                      " : ~pyarts.arts.",
                      grpname,
                      ", optional\n    ",
                      unwrap_stars(short_doc(varname)),
                      " See :attr:`~pyarts.workspace.Workspace.",
                      varname,
                      "`, defaults to ``self.",
                      varname,
                      "`` **[IN]**");
  }

  for (Index i = 0; i < method.GIn().nelem(); i++) {
    const auto& varname = method.GIn()[i];
    if (is_goutput(varname)) continue;
    const auto& grpname =
        compose_generic_groups(method.GInType()[i], method.GInSpecType()[i]);
    const auto& defval = method.GInDefault()[i];
    const bool has_defval = defval not_eq NODEF;
    const String opt{has_defval ? ", optional" : ""};
    const String optval{
        has_defval ? var_string(" Defaults to ``", to_defval_str(defval, global_data::wsv_groups[method.GInType()[i]].name), "``") : ""};
    out += var_string('\n',
                      varname,
                      " : ",
                      grpname,
                      opt,
                      "\n    ",
                      unwrap_stars(until_first_newline(method.GInDescription()[i])),
                      optval,
                      " **[IN]**");
  }

  static const Index verbpos = global_data::WsvMap.at("verbosity");
  if (not(is_input(verbpos) or is_output(verbpos))) {
    out +=
        "\nverbosity : ~pyarts.arts.Verbosity\n    ARTS verbosity. See :attr:`~pyarts.workspace.Workspace.verbosity`, defaults to ``self.verbosity`` **[IN]**";
  }
  fix();

  out += ")--";

  return out;
}
