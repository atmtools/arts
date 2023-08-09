#include <matpack.h>
#include <mystring.h>
#include <workspace.h>

#include <algorithm>
#include <iomanip>
#include <iterator>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "array.h"
#include "compare.h"
#include "debug.h"
#include "nonstd.h"
#include "python_interface/pydocs.h"
#include "workspace_agendas.h"

String as_pyarts(const String& x) {
  const static auto& wsgs = internal_workspace_groups();
  const static auto& wsvs = workspace_variables();
  const static auto& wsms = workspace_methods();

  const auto found_in = [&](auto& map) { return map.find(x) not_eq map.end(); };

  if (found_in(wsgs)) return var_string(":class:`~pyarts.arts.", x, '`');
  if (found_in(wsms))
    return var_string(":func:`~pyarts.workspace.Workspace.", x, '`');
  if (found_in(wsvs))
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

bool str_compare_nocase(const std::string& lhs, const std::string& rhs) {
  auto str_toupper = [](std::string s) {
    std::transform(s.begin(),
                   s.end(),
                   s.begin(),
                   [](unsigned char c) { return std::toupper(c); }  // correct
    );
    return s;
  };

  return str_toupper(lhs) < str_toupper(rhs);
};

std::string fix_newlines(std::string x) {
  constexpr auto old = "\n\n\n";
  constexpr auto new_ = "\n\n";
  auto ptr = x.find(old);
  while (ptr not_eq x.npos) {
    x.replace(ptr, 3, new_);
    ptr = x.find(old);
  }

  while (x.back() == '\n') x.pop_back();
  x.push_back('\n');
  return x;
}

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
  const static auto& wsas = internal_workspace_agendas();
  const static auto& wsvs = workspace_variables();

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

  auto& ag = wsas.at(x);
  auto& out_ind = ag.output;
  auto& in_ind = ag.input;

  const auto output = [&](const std::string& ind) {
    return std::ranges::any_of(out_ind, Cmp::eq(ind));
  };

  const auto input = [&](const std::string& ind) {
    return std::ranges::any_of(in_ind, Cmp::eq(ind));
  };

  std::vector<AgendaIO> writer;
  writer.reserve(out_ind.size());
  for (auto& var : out_ind) {
    writer.emplace_back(AgendaIO{input(var), true, wsvs.at(var).type, var});
  }

  for (auto& var : in_ind) {
    if (not output(var))
      writer.emplace_back(AgendaIO{true, false, wsvs.at(var).type, var});
  }

  constexpr matpack::matpack_constant_data<std::string_view, 2, 2> inout{
      "[ERROR]", "[OUT]", "[IN]", "[INOUT]"};
  for (auto& var : writer) {
    out += var_string(var.name,
                      " : ~pyarts.arts.",
                      var.group,
                      '\n',
                      "    ",
                      unwrap_stars(short_doc(var.name)),
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
  const static auto& wsgs = internal_workspace_groups();
  const static auto& wsvs = workspace_variables();
  const static auto& wsms = internal_workspace_methods();

  const auto found_in = [&](auto& map) { return map.find(x) not_eq map.end(); };

  if (found_in(wsgs)) return until_first_newline(wsgs.at(x).desc);
  if (found_in(wsms)) return until_first_newline(wsms.at(x).desc);
  if (found_in(wsvs)) return until_first_newline(wsvs.at(x).desc);

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

String compose_generic_groups(const String& grps) {
  return var_string("~pyarts.arts.", grps);
}

String to_defval_str(const Wsv& wsv, const String& group) {
  std::string out =
      std::visit([](auto& a) { return var_string(*a); }, wsv.value);

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

String method_docs(const String& name) try {
  const static auto& wsms = internal_workspace_methods();
  const static auto& wsvs = workspace_variables();

  //! WARNING: Raw method
  const auto& method = wsms.at(name);
  String out;

  const auto is_output = [&](auto ind) {
    return std::ranges::any_of(method.out, Cmp::eq(ind));
  };
  const auto is_input = [&](auto ind) {
    return std::ranges::any_of(method.in, Cmp::eq(ind));
  };
  const auto is_ginput = [&](auto& ind) {
    return std::ranges::any_of(method.gin, Cmp::eq(ind));
  };
  const auto is_goutput = [&](auto& ind) {
    return std::ranges::any_of(method.gout, Cmp::eq(ind));
  };
  const auto fix = [&] {
    remove_trailing(out, '\n');
    out += "\n";
  };

  out = "--(";
  out += unwrap_stars(method.desc);
  fix();

  out += "\nAuthor(s):";
  for (auto& author : method.author) out += " " + author + ", ";
  out.pop_back();  // Remove ' '
  out.pop_back();  // Remove ','
  fix();

  out += "\nParameters\n----------";
  for (auto& varname : method.out) {
    const String io = is_input(varname) ? "INOUT" : "OUT";
    const auto& wsv = wsvs.at(varname);
    const auto& grpname = wsv.type;
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

  for (std::size_t i = 0; i < method.gout.size(); i++) {
    const auto& varname = method.gout[i];
    const String io = is_ginput(varname) ? "INOUT" : "OUT";
    const auto& grpname = compose_generic_groups(method.gout_type[i]);
    out += var_string('\n',
                      varname,
                      " : ",
                      grpname,
                      "\n    ",
                      unwrap_stars(until_first_newline(method.gout_desc[i])),
                      " **[",
                      io,
                      "]**");
  }

  for (auto varname : method.in) {
    if (is_output(varname)) continue;

    const auto& wsv = wsvs.at(varname);
    const auto& grpname = wsv.type;
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

  for (std::size_t i = 0; i < method.gin.size(); i++) {
    const auto& varname = method.gin[i];
    if (is_goutput(varname)) continue;
    const auto& grpname = compose_generic_groups(method.gin_type[i]);
    const auto& defval = method.gin_value[i];
    const bool has_defval = bool(defval);
    const String opt{has_defval ? ", optional" : ""};
    const String optval{
        has_defval ? var_string(" Defaults to ``", to_defval_str(*defval, method.gin_type[i]), "``")
                   : ""};
    out += var_string('\n',
                      varname,
                      " : ",
                      grpname,
                      opt,
                      "\n    ",
                      unwrap_stars(until_first_newline(method.gin_desc[i])),
                      optval,
                      " **[IN]**");
  }

  fix();

  out += ")--";

  return out;
} catch (std::out_of_range& e) {
  throw std::runtime_error(var_string("Cannot find: ", std::quoted(name)));
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("Error in method_docs(", name, "): ", e.what()));
}

String variable_used_by(const String& name) {
  const static auto& wsms = internal_workspace_methods();
  const static auto& wsas = internal_workspace_agendas();

  struct wsv_io {
    std::vector<String> wsm_out;
    std::vector<String> wsm_in;
    std::vector<String> ag_out;
    std::vector<String> ag_in;
  };

  wsv_io usedocs;
  for (auto& [mname, method] : wsms) {
    if (std::ranges::any_of(method.out, Cmp::eq(name)))
      usedocs.wsm_out.emplace_back(mname);

    if (std::ranges::any_of(method.in, Cmp::eq(name)))
      usedocs.wsm_in.emplace_back(mname);
  }

  for (auto& [aname, agenda] : wsas) {
    if (std::ranges::any_of(agenda.output, Cmp::eq(name)))
      usedocs.ag_out.emplace_back(aname);

    if (std::ranges::any_of(agenda.input, Cmp::eq(name)))
      usedocs.ag_in.emplace_back(aname);
  }

  String val;
  if (usedocs.wsm_out.size()) {
    val += var_string("\n\nWorkspace methods that can generate ",
                      name,
                      "\n",
                      String(36 + name.size(), '-'),
                      "\n\n.. hlist::",
                      "\n    :columns: ",
                      hlist_num_cols(usedocs.wsm_out),
                      "\n");
    for (auto& m : usedocs.wsm_out)
      val += var_string("\n    * :func:`~pyarts.workspace.Workspace.", m, '`');
    val += "\n";
  }

  if (usedocs.wsm_in.size()) {
    val += var_string("\n\nWorkspace methods that require ",
                      name,
                      "\n",
                      String(31 + name.size(), '-'),
                      "\n\n.. hlist::",
                      "\n    :columns: ",
                      hlist_num_cols(usedocs.wsm_in),
                      "\n");
    for (auto& m : usedocs.wsm_in)
      val += var_string("\n    * :func:`~pyarts.workspace.Workspace.", m, '`');
    val += "\n";
  }

  if (usedocs.ag_out.size()) {
    val += var_string("\n\nWorkspace agendas that can generate ",
                      name,
                      "\n",
                      String(36 + name.size(), '-'),
                      "\n\n.. hlist::",
                      "\n    :columns: ",
                      hlist_num_cols(usedocs.ag_out),
                      "\n");
    for (auto& m : usedocs.ag_out)
      val += var_string("\n    * :attr:`~pyarts.workspace.Workspace.", m, '`');
    val += "\n";
  }

  if (usedocs.ag_in.size()) {
    val += var_string("\n\nWorkspace agendas that require ",
                      name,
                      "\n",
                      String(31 + name.size(), '-'),
                      "\n\n.. hlist::",
                      "\n    :columns: ",
                      hlist_num_cols(usedocs.ag_in),
                      "\n");
    for (auto& m : usedocs.ag_in)
      val += var_string("\n    * :attr:`~pyarts.workspace.Workspace.", m, '`');
    val += "\n";
  }

  return val;
}