#include "workspace_agenda_class.h"

#include <auto_wsa.h>
#include <auto_wsm.h>

#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

#include "workspace_class.h"
#include "workspace_method_class.h"

Agenda::Agenda(std::string n) : name(std::move(n)), methods() {} 

void Agenda::add(const Method& method) {
  checked = false;

  const std::vector<std::string>& outs = method.get_outs();
  const std::vector<std::string>& ins = method.get_ins();

  for (auto& var: outs) {
    if (var.front() == '<') {
      throw std::runtime_error(var_string(
          "Undefined output variable ", std::quoted(var), " in agenda ", std::quoted(name)));
    }
  }

  for (auto& var: ins) {
    if (var.front() == '<') {
      const auto& wsm = workspace_methods();

      const auto& wsmr = wsm.at(method.get_name());
      auto ptr = wsmr.defs.find(var);
      if (ptr == wsmr.defs.end()){
        throw std::runtime_error(var_string(
            "Undefined input variable ", std::quoted(var), " in agenda ", std::quoted(name)));
      }
      
      methods.emplace_back(var, ptr->second);
    }
  }

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

  for (auto& method : methods) {
    std::ranges::copy_if(
        method.get_ins(), std::back_inserter(copy), is_in(method.get_outs()));
    std::ranges::copy_if(method.get_ins(),
                         std::back_inserter(share),
                         is_not_in(method.get_outs()));

    // FIXME: Is this really necessary???
    std::ranges::copy_if(
        method.get_outs(), std::back_inserter(share), is_in(method.get_ins()));
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

  checked = true;
} catch (std::exception& e) {
  throw std::runtime_error(var_string(
      "Error finalizing agenda ", std::quoted(name), "\n", e.what()));
}

void Agenda::copy_workspace(Workspace& out, const Workspace& in) const try {
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
} catch (std::exception& e) {
  throw std::runtime_error(var_string(
      "Cannot get value from\n\n", in, "\ninto\n\n", out, "\n", e.what()));
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
      var_string("Error executing agenda ", std::quoted(name), "\n", e.what()));
}

std::ostream& operator<<(std::ostream& os, const Agenda& a) {
  os << "Agenda " << std::quoted(a.name) << ":\n";
  os << "  Methods:\n";
  for (auto& method : a.methods) {
    os << "    " << method.get_name() << "\n";
  }
  return os;
}

bool Agenda::has_method(const std::string& method) const {
  for (auto& m : methods) {
    if (m.get_name() == method) {
      return true;
    }
  }
  return false;
}
