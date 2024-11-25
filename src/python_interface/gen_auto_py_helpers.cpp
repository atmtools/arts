#include <matpack.h>
#include <mystring.h>
#include <workspace.h>

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "arts_options.h"
#include "compare.h"
#include "debug.h"
#include "nonstd.h"
#include "python_interface/pydocs.h"
#include "workspace_agendas.h"
#include "workspace_variables.h"

String as_pyarts(const String& x) try {
  const auto& wsgs = internal_workspace_groups();
  const auto& wsvs = workspace_variables();
  const auto& wsms = workspace_methods();

  const auto found_in = [&](auto& map) { return map.find(x) not_eq map.end(); };
  const auto found_in_options = [&](auto& key) {
    return std::ranges::any_of(
        internal_options(), Cmp::eq(key), &EnumeratedOption::name);
  };

  if (found_in_options(x) or found_in(wsgs))
    return std::format(":class:`~pyarts.arts.{}`", x);
  if (found_in(wsms))
    return std::format(":func:`~pyarts.workspace.Workspace.{}`", x);
  if (found_in(wsvs))
    return std::format(":attr:`~pyarts.workspace.Workspace.{}`", x);

  throw std::invalid_argument(
      std::format(R"("{0}"  is not a valid group, method, or workspace variable.

If it is a GIN or GOUT, consider representing it as ``{0}`` instead?
If it is an old or deleted method or variable or group, please remove it from the documentation! 
  )",
                  x));
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Could not convert \"{}\" to pyarts:\n{}",
                  x,
                  std::string_view(e.what())));
}

uint32_t hlist_num_cols(const std::vector<String>& v1,
                        const std::vector<String>& v2) {
  return (v1.size() + v2.size()) < 5 ? 1 : 2;
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
  while (not x.empty() and x.back() == '\n') x.pop_back();
  x.push_back('\n');
  return x;
}

String unwrap_stars(String x) try {
  const auto find = [&](auto p) {
    p = std::min(p, x.end());
    return std::find(p, x.end(), '*');
  };
  const auto idx = [&](auto p) { return std::distance(x.begin(), p); };
  const auto mv  = [&](auto p) {
    do p++;
    while (p < x.end() and not nonstd::isspace(*p) and (*p) not_eq '*');
    return std::min(p, x.end());
  };

  std::vector<Index> start_pos;

  auto ptr = find(x.begin());
  while (ptr < x.end()) {
    const Index pos = idx(ptr);
    ptr             = mv(ptr);
    if (ptr < x.end() and *(ptr) == '*') start_pos.push_back(pos);
    ptr = find(ptr + 1);
  }

  while (not start_pos.empty()) {
    const Index pos = start_pos.back();
    start_pos.pop_back();
    const auto start_ptr    = x.begin() + pos;
    const auto end_ptr      = 1 + mv(start_ptr);
    const auto replace_name = x.substr(pos, std::distance(start_ptr, end_ptr));
    x.replace(start_ptr,
              end_ptr,
              as_pyarts(replace_name.substr(1, replace_name.size() - 2)));
  }

  return x;
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Could not unwrap stars: {}", std::string_view(e.what())));
}

String get_agenda_io(const String& x) try {
  const auto& wsas = internal_workspace_agendas();
  const auto& wsvs = workspace_variables();

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

  auto& ag      = wsas.at(x);
  auto& out_ind = ag.output;
  auto& in_ind  = ag.input;

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
    out += std::format(R"({0} : ~pyarts.arts.{1}
    {2} See :attr: `~pyarts.workspace.Workspace.{0}` **{3}**
)",
                       var.name,
                       var.group,
                       unwrap_stars(short_doc(var.name)),
                       inout[var.in][var.out]);
  }

  return writer.size() ? out + "\n\n" : "";
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Could not get agenda IO for \"{}\":\n{}",
                  x,
                  std::string_view(e.what())));
}

String until_first_newline(const String& x) {
  std::string out{x.begin(), std::find(x.begin(), x.end(), '\n')};
  if (not out.ends_with('.')) out.push_back('.');
  return out;
}

String short_doc(const String& x) try {
  const auto& wsgs = internal_workspace_groups();
  const auto& wsvs = workspace_variables();
  const auto& wsms = internal_workspace_methods();

  const auto found_in = [&](auto& map) { return map.find(x) not_eq map.end(); };

  if (found_in(wsgs)) return until_first_newline(wsgs.at(x).desc);
  if (found_in(wsms)) return until_first_newline(wsms.at(x).desc);
  if (found_in(wsvs)) return until_first_newline(wsvs.at(x).desc);

  throw std::invalid_argument(
      std::format(R"("{0}" is not a valid group, method, or workspace variable.
If it If it is a GIN or GOUT, consider representing it as ``{0}`` instead?",
If it is an old or deleted method or variable or group, please remove it from the documentation!
)",
                  x));
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Could not get short doc for \"{}\":\n{}",
                  x,
                  std::string_view(e.what())));
}

void remove_trailing(String& x, char n) {
  while (x.ends_with(n)) x.pop_back();
}

String compose_generic_groups(const String& grps) {
  return std::format("~pyarts.arts.{}", grps);
}

String to_defval_str(const Wsv& wsv) try {
  const auto& group = wsv.type_name();

  std::string out =
      std::visit([](auto& a) { return std::format("{}", *a); }, wsv.value());

  while (not out.empty() and out.front() == ' ') out.erase(out.begin());
  while (not out.empty() and out.back() == ' ') out.pop_back();

  if (group == "String" and
      (out.empty() or (out.front() not_eq '"' and out.back() not_eq '"'))) {
    return std::format("\"{}\"", out);
  }

  if (out.size() == 0) {
    if (group.starts_with("Array") or group == "Vector" or group == "Matrix" or
        group == "Tensor3" or group == "Tensor4" or group == "Tensor5" or
        group == "Tensor6" or group == "Tensor7")
      return "[]";

    if (group == "Numeric" or group == "Index") return "0";

    return std::format("pyarts.arts.{}()", group);
  }

  return out;
} catch (std::bad_variant_access&) {
  throw std::runtime_error("Cannot convert to defval string");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Error in to_defval_str: {}\n", std::string_view(e.what())));
}

String method_docs(const String& name) try {
  const auto& wsms = internal_workspace_methods();
  const auto& wsvs = workspace_variables();

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

  out  = "--(";
  out += unwrap_stars(method.desc);
  fix();

  out += "\nAuthor(s):";
  for (auto& author : method.author) out += " " + author + ", ";
  out.pop_back();  // Remove ' '
  out.pop_back();  // Remove ','
  fix();

  out += "\nParameters\n----------";
  for (auto& varname : method.out) {
    const String io      = is_input(varname) ? "INOUT" : "OUT";
    const auto& wsv      = wsvs.at(varname);
    const auto& grpname  = wsv.type;
    out                 += std::format(R"(
{0} : ~pyarts.arts.{1}, optional
    {2} See :attr:`~pyarts.workspace.Workspace.{0}`, defaults to ``self.{0}`` **[{3}]**)",
                       varname,
                       grpname,
                       unwrap_stars(short_doc(varname)),
                       io);
  }

  for (std::size_t i = 0; i < method.gout.size(); i++) {
    const auto& varname  = method.gout[i];
    const String io      = is_ginput(varname) ? "INOUT" : "OUT";
    const auto& grpname  = compose_generic_groups(method.gout_type[i]);
    out                 += std::format(R"(
{0} : {1}
    {2}  Defaults to create and/or use ``self.{0}`` : :class:`{1}`. **[{3}]**)",
                       varname,
                       grpname,
                       unwrap_stars(until_first_newline(method.gout_desc[i])),
                       io);
  }

  for (auto&& varname : method.in) {
    if (is_output(varname)) continue;

    const auto& wsv      = wsvs.at(varname);
    const auto& grpname  = wsv.type;
    out                 += std::format(R"(
{0} : ~pyarts.arts.{1}, optional
    {2} See :attr:`~pyarts.workspace.Workspace.{0}`, defaults to ``self.{0}`` **[IN]**)",
                       varname,
                       grpname,
                       unwrap_stars(short_doc(varname)));
  }

  for (std::size_t i = 0; i < method.gin.size(); i++) {
    const auto& varname = method.gin[i];
    if (is_goutput(varname)) continue;
    const auto& grpname   = compose_generic_groups(method.gin_type[i]);
    const auto& defval    = method.gin_value[i];
    const bool has_defval = bool(defval);
    const String opt{has_defval ? ", optional" : ""};
    const String optval{
        has_defval ? std::format(" Defaults to ``{}``", to_defval_str(*defval))
                   : ""};
    out += std::format(R"({0} : {1}{2}
    {2}{3} **[IN]**)",
                      varname,
                      grpname,
                      opt,
                      unwrap_stars(until_first_newline(method.gin_desc[i])),
                      optval);
  }

  fix();

  out += ")--";

  return out;
} catch (std::out_of_range& e) {
  throw std::runtime_error(std::format("Cannot find: \"{}\"", name));
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Error in method_docs({}): {}", name, std::string_view(e.what())));
}

String variable_used_by(const String& name) {
  const auto& wsvs = internal_workspace_variables();
  const auto& wsms = internal_workspace_methods();
  const auto& wsas = internal_workspace_agendas();

  struct wsv_io {
    std::vector<String> wsm_out;
    std::vector<String> wsm_in;
    std::vector<String> ag_out;
    std::vector<String> ag_in;
    std::vector<String> wsvs;
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

  std::ranges::copy_if(wsvs | std::views::keys,
                       std::back_inserter(usedocs.wsvs),
                       [&name](auto& vname) {
                         return workspace_variables_keywords_match(vname, name);
                       });

  std::ranges::sort(usedocs.wsm_out);
  std::ranges::sort(usedocs.wsm_in);
  std::ranges::sort(usedocs.ag_out);
  std::ranges::sort(usedocs.ag_in);
  std::ranges::sort(usedocs.wsvs);

  String val;
  if (usedocs.wsm_out.size() or usedocs.wsm_out.size()) {
    std::array<std::vector<std::string>, 3> io;
    for (auto& m : usedocs.wsm_out) {
      if (std::ranges::any_of(usedocs.wsm_in, Cmp::eq(m))) {
        io[1].emplace_back(m);
      } else {
        io[2].emplace_back(m);
      }
    }
    for (auto& m : usedocs.wsm_in) {
      if (std::ranges::none_of(usedocs.wsm_out, Cmp::eq(m))) {
        io[0].emplace_back(m);
      }
    }

    if (io[0].size()) {
      val += std::format(R"(

.. rubric:: Input to workspace methods

.. hlist::
    :columns: {0}
)",
                        hlist_num_cols(io[0]));
      for (auto& m : io[0]) {
        val +=
            std::format("\n    * :func:`~pyarts.workspace.Workspace.{}`", m);
      }
    }

    if (io[1].size()) {
      val += std::format(R"(

.. rubric:: Modified by workspace methods

.. hlist::
    :columns: {0}
)",
                        hlist_num_cols(io[1]));
      for (auto& m : io[1]) {
        val +=
            std::format("\n    * :func:`~pyarts.workspace.Workspace.{}`", m);
      }
    }

    if (io[2].size()) {
      val += std::format(R"(

.. rubric:: Output from workspace methods

.. hlist::
    :columns: {0}
)",
                        hlist_num_cols(io[2]));
      for (auto& m : io[2]) {
        val +=
            std::format("\n    * :func:`~pyarts.workspace.Workspace.{}`", m);
      }
    }
    val += "\n";
  }

  if (usedocs.ag_out.size()) {
    std::array<std::vector<std::string>, 3> io;
    for (auto& m : usedocs.ag_out) {
      if (std::ranges::any_of(usedocs.ag_in, Cmp::eq(m))) {
        io[1].emplace_back(m);
      } else {
        io[2].emplace_back(m);
      }
    }
    for (auto& m : usedocs.ag_in) {
      if (std::ranges::none_of(usedocs.ag_out, Cmp::eq(m))) {
        io[0].emplace_back(m);
      }
    }

    if (io[0].size()) {
      val += std::format(R"(
.. rubric:: Input to workspace agendas

.. hlist::
    :columns: {0}
)",
                        hlist_num_cols(io[0]));
      for (auto& m : io[0]) {
        val +=
            std::format("\n    * :attr:`~pyarts.workspace.Workspace.{}`", m);
      }
    }

    if (io[1].size()) {
      val += std::format(R"(

.. rubric:: Modified by workspace agendas

.. hlist::
    :columns: {0}
)",
                        hlist_num_cols(io[1]));
      for (auto& m : io[1]) {
        val +=
            std::format("\n    * :attr:`~pyarts.workspace.Workspace.{}`", m);
      }
    }

    if (io[2].size()) {
      val += std::format(R"(
.. rubric:: Output from workspace agendas

.. hlist::
    :columns: {0}
)",
                        hlist_num_cols(io[2]));
      for (auto& m : io[2]) {
        val +=
            std::format("\n    * :attr:`~pyarts.workspace.Workspace.{}`", m);
      }
    }
    val += "\n";
  }

  if (usedocs.wsvs.size()) {
    val += std::format(R"(

.. rubric:: Related workspace variables

.. hlist::
    :columns: {0}
)",
                      hlist_num_cols(usedocs.wsvs));
    for (auto& m : usedocs.wsvs) {
      val += std::format("\n    *  :attr:`~pyarts.workspace.Workspace.{}`", m);
    }
    val += "\n";
  }

  return val;
}
