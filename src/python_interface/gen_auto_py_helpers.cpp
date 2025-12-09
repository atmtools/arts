#include <arts_options.h>
#include <compare.h>
#include <debug.h>
#include <nonstd.h>
#include <workspace.h>
#include <workspace_agendas.h>
#include <workspace_group_friends.h>
#include <workspace_meta_methods.h>
#include <workspace_method_extra_doc.h>
#include <workspace_variables.h>

#include <algorithm>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "pydocs.h"

namespace {
String rawify(const String& x) {
  std::stringstream os{x};

  std::string out{};

  for (std::string line; std::getline(os, line);) {
    out += std::format(R"(R"-x-({})-x-" "\n"
)",
                       line);
  }

  return out;
}

String as_pyarts(const String& x) try {
  const auto& wsgs          = internal_workspace_groups();
  const auto& wsvs          = workspace_variables();
  const auto& wsms          = workspace_methods();
  const auto& group_friends = workspace_group_friends();

  const auto found_in = [&](auto& map) { return map.find(x) not_eq map.end(); };
  const auto found_in_options = [&](auto& key) {
    return stdr::any_of(
        internal_options(), Cmp::eq(key), &EnumeratedOption::name);
  };

  if (found_in_options(x) or found_in(wsgs) or found_in(group_friends))
    return std::format(":class:`~pyarts3.arts.{}`", x);
  if (found_in(wsms))
    return std::format(":func:`~pyarts3.workspace.Workspace.{}`", x);
  if (found_in(wsvs))
    return std::format(":attr:`~pyarts3.workspace.Workspace.{}`", x);

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
}  // namespace

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
  while (not x.empty() and nonstd::isspace(x.back())) x.pop_back();
  x.push_back('\n');
  return x;
}

namespace {
bool isnotletter(char c) { return not nonstd::isabc(c); }
}  // namespace

String unwrap_stars(const String& x) try {
  auto ptr       = x.begin();
  const auto end = x.end();

  String out{};
  out.reserve(x.size() + 1000);

  bool maybe_first = true;
  while (ptr < end) {
    const auto c = *ptr;
    ++ptr;

    if (not maybe_first or c != '*' or ptr == end) {
      out.push_back(c);
      maybe_first = nonstd::isspace(c) or nonstd::isbracket(c);
      continue;
    }

    const auto first = ptr;

    if (isnotletter(*first)) {
      out.push_back('*');
      out.push_back(*first);
      ptr++;
      continue;
    }

    while (ptr < end and *ptr != '*') ++ptr;
    out += as_pyarts({first, ptr});
    ptr++;  // May be beyond end
  }

  return out;
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Could not unwrap stars: {}", std::string_view(e.what())));
}

String get_agenda_io(const String& x) try {
  const auto& wsas = internal_workspace_agendas();

  String out{};

  auto& ag = wsas.at(x);

  if (not ag.output.empty()) {
    out += std::format(R"(
.. rubric:: Agenda output

.. hlist::
    :columns: {}

)",
                       1 + (ag.output.size() > 3));

    for (auto& varname : ag.output) {
      out += std::format(R"(    * :attr:`~pyarts3.workspace.Workspace.{}`
)",
                         varname);
    }
  }

  if (not ag.input.empty()) {
    out += std::format(R"(
.. rubric:: Agenda input

.. hlist::
    :columns: {}

)",
                       1 + (ag.input.size() > 3));

    for (auto& varname : ag.input) {
      out += std::format(R"(    * :attr:`~pyarts3.workspace.Workspace.{}`
)",
                         varname);
    }
  }

  return out;
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Could not get agenda IO for \"{}\":\n{}",
                  x,
                  std::string_view(e.what())));
}

namespace {
String until_first_newline(const String& x) {
  std::string out{x.begin(), std::find(x.begin(), x.end(), '\n')};
  if (not out.ends_with('.')) out.push_back('.');
  return out;
}
}  // namespace

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

namespace {
void remove_trailing(String& x, char n) {
  while (x.ends_with(n)) x.pop_back();
}

String compose_generic_groups(const String& grps) {
  return std::format("~pyarts3.arts.{}", grps);
}
}  // namespace

String to_defval_str(const Wsv& wsv, const std::string_view x) try {
  const auto& group = wsv.type_name();

  std::string out = wsv.vformat("{}"sv);

  while (not out.empty() and out.front() == ' ') out.erase(out.begin());
  while (not out.empty() and out.back() == ' ') out.pop_back();

  if (group == "String" and
      (out.empty() or (out.front() not_eq '"' and out.back() not_eq '"'))) {
    return std::format("{1}\"{0}\"{1}", out, x);
  }

  if (group == "Agenda") {
    return unwrap_stars(wsv.get<Agenda>().sphinx_list("#. "sv));
  }

  if (out.size() == 0) {
    if (group.starts_with("Array") or group == "Vector" or group == "Matrix" or
        group == "Tensor3" or group == "Tensor4" or group == "Tensor5" or
        group == "Tensor6" or group == "Tensor7")
      return std::format("{0}[]{0}", x);

    if (group == "Numeric" or group == "Index")
      return std::format("{0}0{0}", x);

    return std::format("{1}pyarts3.arts.{0}(){1}", group, x);
  }

  return std::format("{1}{0}{1}", out, x);
} catch (std::bad_variant_access&) {
  throw std::runtime_error("Cannot convert to defval string");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Error in to_defval_str: {}\n", std::string_view(e.what())));
}

String method_docs(const String& name) try {
  const auto& wsms       = internal_workspace_methods();
  const auto& wsvs       = workspace_variables();
  const auto wsadoc      = get_agenda_enum_documentation();
  const auto& wsms_extra = workspace_method_extra_doc();
  const auto& wsms_meta  = internal_meta_methods();

  //! WARNING: Raw method
  const auto& method = wsms.at(name);
  String out;

  const auto is_output = [&](auto&& ind) {
    return stdr::any_of(method.out, Cmp::eq(ind));
  };
  const auto is_input = [&](auto&& ind) {
    return stdr::any_of(method.in, Cmp::eq(ind));
  };
  const auto is_ginput = [&](auto& ind) {
    return stdr::any_of(method.gin, Cmp::eq(ind));
  };
  const auto is_goutput = [&](auto&& ind) {
    return stdr::any_of(method.gout, Cmp::eq(ind));
  };
  const auto fix = [&] {
    remove_trailing(out, '\n');
    out += "\n";
  };

  out += unwrap_stars(method.desc);
  fix();

  out += std::format("\nAuthor{}: {:,}",
                     method.author.size() > 1 ? "s"sv : ""sv,
                     method.author);
  fix();

  std::vector<String> metamethods =
      wsms_meta | stdv::filter([&name](const auto& mm) {
        return stdr::any_of(mm.methods, Cmp::eq(name));
      }) |
      stdv::transform([](auto&& m) -> String { return m.name; }) |
      stdr::to<std::vector<String>>();
  stdr::sort(metamethods);
  out += metamethods.empty()
             ? ""s
             : std::format(
                   R"(
.. rubric:: Used by wrapper method{0}

.. hlist::
    :columns: {1}
{2}
)",
                   metamethods.size() > 1 ? "s"sv : ""sv,
                   hlist_num_cols(metamethods),
                   metamethods | stdv::transform([](const auto& m) {
                     return std::format(
                         "\n    * :func:`~pyarts3.workspace.Workspace.{}`", m);
                   }) | stdr::to<std::vector<String>>());
  fix();

  out += "\nParameters\n----------";
  for (auto& varname : method.out) {
    const String io      = is_input(varname) ? "INOUT" : "OUT";
    const auto& wsv      = wsvs.at(varname);
    const auto& grpname  = wsv.type;
    out                 += std::format(R"(
{0} : ~pyarts3.arts.{1}, optional
    {2} See :attr:`~pyarts3.workspace.Workspace.{0}`, defaults to ``self.{0}`` **[{3}]**)",
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
{0} : ~pyarts3.arts.{1}, optional
    {2} See :attr:`~pyarts3.workspace.Workspace.{0}`, defaults to ``self.{0}`` **[IN]**)",
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
    const String optval{has_defval ? std::format(" Defaults to {}",
                                                 to_defval_str(*defval, "``"sv))
                                   : ""};
    out += std::format(R"(
{0} : {1}{2}
    {3}{4} **[IN]**)",
                       varname,
                       grpname,
                       opt,
                       unwrap_stars(until_first_newline(method.gin_desc[i])),
                       optval);
  }

  if (method.return_type != "void") {
    out += std::format(
        R"(
Returns
-------
opt : {0}
    {1})",
        std::format("~pyarts3.arts.{}{}",
                    method.return_type == "Workspace" ? "Cxx"sv : ""sv,
                    method.return_type),
        unwrap_stars(until_first_newline(method.return_desc)));
  }

  fix();

  if (auto ptr = wsadoc.find(name); ptr != wsadoc.end()) {
    out += std::format(R"(

.. rubric:: Valid options

These are the valid options for the ``{}`` method.
The listed method calls describe the order of the agenda calls for each ``option``.

)",
                       name);

    for (auto& [opt, doc] : ptr->second) {
      out += std::format(R"(

------------------------------------------------------------

``{}(option="{}")``

{}

)",
                         name,
                         opt,
                         unwrap_stars(doc));
    }
    fix();
  }

  if (wsms_extra.find(name) != wsms_extra.end()) {
    out += std::format(R"(

.. rubric:: Extra

{}

)",
                       unwrap_stars(wsms_extra.at(name)));
  }

  return rawify(out);
} catch (std::out_of_range& e) {
  throw std::runtime_error(std::format("Cannot find: \"{}\"", name));
} catch (std::exception& e) {
  throw std::runtime_error(std::format(
      "Error in method_docs({}): {}", name, std::string_view(e.what())));
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

    wsv_io(std::vector<String>&& wsm_out_,
           std::vector<String>&& wsm_in_,
           std::vector<String>&& ag_out_,
           std::vector<String>&& ag_in_,
           std::vector<String>&& wsvs_)
        : wsm_out(std::move(wsm_out_)),
          wsm_in(std::move(wsm_in_)),
          ag_out(std::move(ag_out_)),
          ag_in(std::move(ag_in_)),
          wsvs(std::move(wsvs_)) {
      stdr::sort(wsm_out);
      stdr::sort(wsm_in);
      stdr::sort(ag_out);
      stdr::sort(ag_in);
      stdr::sort(wsvs);
    }
  };

  const auto to_vstring = stdr::to<std::vector<String>>();

  const auto filter_wsmout = stdv::filter([&name](auto& m) {
                               return stdr::any_of(m.second.out, Cmp::eq(name));
                             }) |
                             stdv::keys;

  const auto filter_wsmin = stdv::filter([&name](auto& m) {
                              return stdr::any_of(m.second.in, Cmp::eq(name));
                            }) |
                            stdv::keys;

  const auto filter_wsaout =
      stdv::filter([&name](auto& a) {
        return stdr::any_of(a.second.output, Cmp::eq(name));
      }) |
      stdv::keys;

  const auto filter_wsain =
      stdv::filter([&name](auto& a) {
        return stdr::any_of(a.second.input, Cmp::eq(name));
      }) |
      stdv::keys;

  const auto filter_wsvs = stdv::keys | stdv::filter([&name](auto& v) {
                             return workspace_variables_keywords_match(v, name);
                           });

  const wsv_io usedocs{wsms | filter_wsmout | to_vstring,
                       wsms | filter_wsmin | to_vstring,
                       wsas | filter_wsaout | to_vstring,
                       wsas | filter_wsain | to_vstring,
                       wsvs | filter_wsvs | to_vstring};

  const auto to_wsmout = stdv::filter([&usedocs](const String& x) {
    return stdr::none_of(usedocs.wsm_in, Cmp::eq(x));
  });

  const auto to_wsminout = stdv::filter([&usedocs](const String& x) {
    return stdr::any_of(usedocs.wsm_in, Cmp::eq(x));
  });

  const auto to_wsmin = stdv::filter([&usedocs](const String& x) {
    return stdr::none_of(usedocs.wsm_out, Cmp::eq(x));
  });

  const auto to_wsaout = stdv::filter([&usedocs](const String& x) {
    return stdr::none_of(usedocs.ag_in, Cmp::eq(x));
  });

  const auto to_wsainout = stdv::filter([&usedocs](const String& x) {
    return stdr::any_of(usedocs.ag_in, Cmp::eq(x));
  });

  const auto to_wsain = stdv::filter([&usedocs](const String& x) {
    return stdr::none_of(usedocs.ag_out, Cmp::eq(x));
  });

  const auto to_attr = stdv::transform([](const String& x) -> String {
    return std::format("    * :attr:`~pyarts3.workspace.Workspace.{}`", x);
  });

  const auto to_func = stdv::transform([](const String& x) -> String {
    return std::format("    * :func:`~pyarts3.workspace.Workspace.{}`", x);
  });

  const auto wsmout   = usedocs.wsm_out | to_wsmout | to_vstring;
  const auto wsminout = usedocs.wsm_out | to_wsminout | to_vstring;
  const auto wsmin    = usedocs.wsm_in | to_wsmin | to_vstring;
  const auto wsaout   = usedocs.ag_out | to_wsaout | to_vstring;
  const auto wsainout = usedocs.ag_out | to_wsainout | to_vstring;
  const auto wsain    = usedocs.ag_in | to_wsain | to_vstring;

  const auto to_str = [to_vstring](auto& arr,
                                   auto& f,
                                   const std::string_view& m) -> std::string {
    if (arr.empty()) return ""s;
    return std::format(R"(

.. rubric:: {3}{2}

.. hlist::
    :columns: {0}

{1:n}
)",
                       hlist_num_cols(arr),
                       arr | f | to_vstring,
                       arr.size() > 1 ? "s"sv : ""sv,
                       m);
  };

  return std::format(
      R"({}{}{}{}{}{}{}

)",
      to_str(wsmin, to_func, "Input to workspace method"),
      to_str(wsminout, to_func, "Modified by workspace method"),
      to_str(wsmout, to_func, "Output from workspace method"),
      to_str(wsain, to_attr, "Input to workspace agenda"),
      to_str(wsainout, to_attr, "Modified by workspace agenda"),
      to_str(wsaout, to_attr, "Output from workspace agenda"),
      to_str(usedocs.wsvs, to_attr, "Related workspace variable"));
}
