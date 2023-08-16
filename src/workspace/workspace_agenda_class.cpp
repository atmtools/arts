#include "workspace_agenda_class.h"

#include <auto_wsa.h>
#include <auto_wsm.h>

#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "debug.h"
#include "workspace_class.h"
#include "workspace_method_class.h"

Agenda::Agenda(std::string n) : name(std::move(n)), methods() {} 

void Agenda::add(const Method& method) {
  checked = false;
  method.add_defaults_to_agenda(*this);
  methods.push_back(method);
}

auto is_in(const std::vector<std::string>& seq) {
  return [&seq](auto& str) { return std::ranges::find(seq, str) != seq.end(); };
}

auto is_not_in(const std::vector<std::string>& seq) {
  return [&seq](auto& str) { return std::ranges::find(seq, str) == seq.end(); };
}

void Agenda::finalize() try {
  static const auto& wsa = workspace_agendas();

  auto ag_ptr = wsa.find(name);
  if (ag_ptr == wsa.end()) {
    throw std::runtime_error("Only pre-compiled agendas may be finalized and thus live on the workspace");
  }

  auto must_out = ag_ptr->second.output;
  auto must_in = ag_ptr->second.input;

  for (auto& method : methods) {
    const auto& ins = method.get_ins();
    const auto& outs = method.get_outs();

    for (auto& i : ins) {
      if (auto ptr = std::ranges::find(must_in, i); ptr != must_in.end()) {
        must_in.erase(ptr);
      } else if (std::ranges::any_of(must_out, Cmp::eq(i))) {
        throw std::runtime_error(
            var_string("method:\n",
                       method,
                       "\nvariable: ",
                       std::quoted(i),
                       '\n',
                       "Despite being an output of the agenda, "
                       "the variable is first encountered as an input to a method"));
      }
    }

    for (auto& i : outs) {
      if (auto ptr = std::ranges::find(must_out, i); ptr != must_out.end()) {
        must_out.erase(ptr);
      }
    }

    std::ranges::copy_if(ins, std::back_inserter(copy), is_in(outs));
    std::ranges::copy_if(ins, std::back_inserter(share), is_not_in(outs));
    std::ranges::copy_if(outs, std::back_inserter(share), is_in(ins));
    std::ranges::copy_if(outs, std::back_inserter(copy), is_not_in(ins));
  }

  if (must_in.size() or must_out.size()) {
    std::ostringstream os;
    os << "Agenda has unused variables:\nRequired output : ";
    for (auto& o: ag_ptr->second.output) {
      os << o << ", ";
    }

    os << "\nRequired input  : ";
    for (auto& o: ag_ptr->second.input) {
      os << o << ", ";
    }

    os << "\nUnused output   : ";
    for (auto& o: must_out) {
      os << o << ", ";
    }

    os << "\nUnused input    : ";
    for (auto& o: must_in) {
      os << o << ", ";
    }
    os << "\nPlease Touch unused output and Ignore unused input if this is intentional.";
    throw std::runtime_error(os.str());
  }

  const auto remove_copies = [](std::vector<string>& seq) {
    std::ranges::sort(seq);
    auto [first, last] = std::ranges::unique(seq);
    seq.erase(first, last);
  };
  remove_copies(share);
  remove_copies(copy);

  if (auto ptr = wsa.find(name); ptr != wsa.end()) {
    std::erase_if(share, is_in(ptr->second.output));
    std::erase_if(copy, is_in(ptr->second.output));
    std::erase_if(share, is_in(ptr->second.input));
  }

  // Don't overlap
  for (auto& str : copy) {
    if (auto ptr = std::ranges::find(share, str); ptr != share.end()) {
      share.erase(ptr);
    }
  }

  std::erase_if(share, [](auto& str) { return str.front() == '_'; });
  std::erase_if(copy, [](auto& str) { return str.front() == '_'; });

  checked = true;
} catch (std::exception& e) {
  throw std::runtime_error(var_string(
      "Error finalizing agenda ", std::quoted(name), '\n', e.what()));
}

void agenda_add_inner_logic(Workspace& out, const Workspace& in, WorkspaceAgendaBoolHandler handle) {
startover:
  for (auto& var: out) {
    if (var.second->holds<Agenda>()) {
      if (not handle.has(var.first)) {
        auto& ag = var.second->get<Agenda>();
        handle.set(var.first);
        ag.copy_workspace(out, in, true);
        goto startover;
      }
    } else if (var.second->holds<ArrayOfAgenda>()) {
      if (not handle.has(var.first)) {
        auto& aag = var.second->get<ArrayOfAgenda>(); 
        handle.set(var.first);
        for (auto& ag: aag) {
          ag.copy_workspace(out, in, true);
        }
        goto startover;
      }
    }
  }
}

void Agenda::copy_workspace(Workspace& out, const Workspace& in, bool share_only) const try {
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
  throw std::runtime_error(var_string(
      "Cannot get value from\n\n", in, "\ninto\n\n", out, '\n', e.what()));
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
  throw std::runtime_error(
      var_string("Error executing agenda ", std::quoted(name), '\n', e.what()));
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
  static const auto& wsa = workspace_agendas();

  auto ptr = wsa.find(a.name);
  const bool named = ptr != wsa.end();
  const bool checked = a.checked;

  os << "Agenda " << a.name;

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
