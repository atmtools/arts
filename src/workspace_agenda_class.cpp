#include "workspace_agenda_class.h"

#include <auto_wsa.h>
#include <auto_wsm.h>

#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "compare.h"
#include "debug.h"
#include "workspace_class.h"
#include "workspace_method_class.h"

void Agenda::add(const Method& method) {
  checked = false;
  method.add_defaults_to_agenda(*this);
  methods.push_back(method);
}

auto is_in(const std::vector<std::string>& seq) {
  return [&seq](auto& str) {
    return std::find(seq.begin(), seq.end(), str) != seq.end();
  };
}

auto is_not_in(const std::vector<std::string>& seq) {
  return [&seq](auto& str) {
    return std::find(seq.begin(), seq.end(), str) == seq.end();
  };
}

void Agenda::finalize(bool fix) try {
  const auto& wsa = workspace_agendas();

  auto ag_ptr                          = wsa.find(name);
  const std::vector<std::string> empty = {};
  const std::vector<std::string>& must_out =
      ag_ptr == wsa.end() ? empty : ag_ptr->second.output;
  const std::vector<std::string>& must_in =
      ag_ptr == wsa.end() ? empty : ag_ptr->second.input;

  std::vector<std::string> ins_first;
  std::vector<std::string> outs_first;
  std::vector<std::string> in_then_out;

  for (const Method& method : methods) {
    const auto& ins  = method.get_ins();
    const auto& outs = method.get_outs();

    std::ranges::copy_if(
        ins, std::back_inserter(ins_first), [&](const std::string& i) {
          const auto cmp = Cmp::eq(i);
          return not std::ranges::any_of(ins_first, cmp) and
                 not std::ranges::any_of(outs_first, cmp) and
                 not std::ranges::any_of(in_then_out, cmp);
        });

    std::ranges::copy_if(
        outs, std::back_inserter(in_then_out), [&](const std::string& i) {
          const auto cmp = Cmp::eq(i);
          return std::ranges::any_of(ins_first, cmp) and
                 not std::ranges::any_of(in_then_out, cmp);
        });

    std::ranges::copy_if(
        outs, std::back_inserter(outs_first), [&](const std::string& i) {
          const auto cmp = Cmp::eq(i);
          return not std::ranges::any_of(outs_first, cmp) and
                 not std::ranges::any_of(ins_first, cmp) and
                 not std::ranges::any_of(in_then_out, cmp);
        });
  }

  auto sort_and_erase_copies = [](std::vector<std::string>& vec) {
    std::ranges::sort(vec);
    auto [s, e] = std::ranges::unique(vec);
    vec.erase(s, e);
  };
  sort_and_erase_copies(ins_first);
  sort_and_erase_copies(outs_first);
  sort_and_erase_copies(in_then_out);

  for (const std::string& i : must_in) {
    if (std::ranges::binary_search(outs_first, i)) {
      throw std::runtime_error(std::format(
          R"(Agenda "{}" first uses "{}" as an input but it is an output)",
          name,
          i));
    }

    if (not std::ranges::binary_search(ins_first, i)) {
      if (fix) {
        methods.emplace_back("Ignore",
                             std::vector<std::string>{i},
                             std::unordered_map<std::string, std::string>{});
      } else {
        throw std::runtime_error(
            std::format(R"(Agenda "{}" does not use "{}")", name, i));
      }
    }
  }

  for (const std::string& o : must_out) {
    if (std::ranges::binary_search(ins_first, o) and
        not std::ranges::binary_search(in_then_out, o)) {
      throw std::runtime_error(std::format(
          R"(Agenda "{}" uses "{}" only as an input but it is Agenda inoutput

Agenda user input:      {:B,}
Agenda user output:     {:B,}
Agenda user inoutput:   {:B,}

Agenda required input:  {:B,}
Agenda required output: {:B,}
)",
          name,
          o,
          ins_first,
          outs_first,
          in_then_out,
          must_in,
          must_out));
    }

    if (not std::ranges::binary_search(outs_first, o)) {
      if (fix) {
        methods.emplace_back("Touch",
                             std::vector<std::string>{o},
                             std::unordered_map<std::string, std::string>{});
      } else {
        throw std::runtime_error(
            std::format(R"(Agenda "{}" does not set "{}")", name, o));
      }
    }
  }

  std::erase_if(ins_first, [&must_in](const auto& str) {
    return std::ranges::any_of(must_in, Cmp::eq(str));
  });

  std::erase_if(ins_first, [&in_then_out](const auto& str) {
    return std::ranges::binary_search(in_then_out, str);
  });

  copy  = in_then_out;
  share = ins_first;

  checked = true;
} catch (std::exception& e) {
  throw std::runtime_error(std::format(R"(Error finalizing agenda "{}"

{}
)",
                                       name,
                                       std::string_view(e.what())));
}

void agenda_add_inner_logic(Workspace& out,
                            const Workspace& in,
                            WorkspaceAgendaBoolHandler handle) {
startover:
  for (auto& var : out) {
    if (var.second.holds<Agenda>()) {
      if (not handle.has(var.first)) {
        auto& ag = var.second.get<Agenda>();
        handle.set(var.first);
        ag.copy_workspace(out, in, true);
        goto startover;
      }
    }
  }
}

void Agenda::copy_workspace(Workspace& out,
                            const Workspace& in,
                            bool share_only) const try {
  if (share_only) {
    for (auto& str : share) {
      if (not out.contains(str)) out.set(str, in.share(str));
    }
    for (auto& str : copy) {
      if (not out.contains(str)) out.set(str, in.share(str));
    }
  } else {
    for (auto& str : share) {
      out.set(str, in.share(str));
    }

    for (auto& str : copy) {
      if (out.contains(str)) {
        //! If copy and share are the same, copy will overwrite share (keep them unique!)
        //! Also if named-agenda call has set output variable, copy will take a copy of variable
        out.overwrite(str, out.copy(str));
      } else {
        out.set(str, in.copy(str));
      }
    }

    WorkspaceAgendaBoolHandler handle;
    handle.set(name);
    agenda_add_inner_logic(out, in, handle);
  }
} catch (std::exception& e) {
  throw std::runtime_error(std::format(
      R"(
Error with workspace copying in Agenda

Workspace contains:
{:s}

{})",
      in,
      e.what()));
}

Workspace Agenda::copy_workspace(const Workspace& in) const {
  Workspace out{WorkspaceInitialization::Empty};
  copy_workspace(out, in);
  return out;
}

void Agenda::execute(Workspace& ws) const try {
  for (auto& method : methods) {
    method(ws);
  }
} catch (std::exception& e) {
  throw std::runtime_error(std::format(R"(Error executing agenda "{}"

{}
)",
                                       name,
                                       std::string_view(e.what())));
}

bool Agenda::has_method(const std::string& method) const {
  for (auto& m : methods) {
    if (m.get_name() == method) {
      return true;
    }
  }
  return false;
}

Agenda::Agenda(std::string n,
               const std::vector<Method>& m,
               const std::vector<std::string>& s,
               const std::vector<std::string>& c,
               bool check)
    : name(std::move(n)), methods(m), share(s), copy(c), checked(check) {}

std::vector<std::string> split(const std::string& s, char c) {
  std::vector<std::string> out;
  std::string tmp;
  for (auto& ch : s) {
    if (ch == c) {
      out.push_back(tmp);
      tmp.clear();
    } else {
      tmp.push_back(ch);
    }
  }
  out.push_back(tmp);
  return out;
}

std::ostream& operator<<(std::ostream& os, const Agenda& a) {
  const auto& wsa = workspace_agendas();

  auto ptr           = wsa.find(a.name);
  const bool named   = ptr != wsa.end();
  const bool checked = a.checked;

  os << "Agenda " << a.name;

  if (a.is_checked())
    os << " (checked)";
  else
    os << " (unchecked)";

  if (checked or named) {
    os << ":\n";

    if (named) {
      os << "  Output : ";
      for (auto& o : ptr->second.output) {
        os << o << ", ";
      }
      os << '\n';

      os << "  Input  : ";
      for (auto& i : ptr->second.input) {
        os << i << ", ";
      }
    }

    if (checked) {
      if (named) {
        os << '\n';
      }

      os << "  Share  : ";
      for (auto& s : a.share) {
        os << s << ", ";
      }
      os << '\n';

      os << "  Copy   : ";
      for (auto& c : a.copy) {
        os << c << ", ";
      }
    }
  }

  return os;
}

std::string Agenda::sphinx_list(const std::string_view prep) const {
  std::string out{};

  for (auto& method : methods) {
    out += std::format("{}{}\n", prep, method.sphinx_list_item());
  }

  return out;
}

void xml_io_stream<Agenda>::write(std::ostream& os,
                                  const Agenda& x,
                                  bofstream* pbofs,
                                  std::string_view name) {
  std::println(R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.get_name(), pbofs);
  xml_write_to_stream(os, x.get_methods(), pbofs);
  xml_write_to_stream(os, x.get_share(), pbofs);
  xml_write_to_stream(os, x.get_copy(), pbofs);
  xml_write_to_stream(os, x.is_checked(), pbofs);

  std::println(R"(</{0}>)", type_name);
}

void xml_io_stream<Agenda>::read(std::istream& is,
                                 Agenda& x,
                                 bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  std::string name{};
  std::vector<Method> methods;
  std::vector<std::string> share{};
  std::vector<std::string> copy{};
  bool checked{false};
  xml_read_from_stream(is, name, pbifs);
  xml_read_from_stream(is, methods, pbifs);
  xml_read_from_stream(is, share, pbifs);
  xml_read_from_stream(is, copy, pbifs);
  xml_read_from_stream(is, checked, pbifs);
  x = Agenda{name, methods, share, copy, checked};

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}