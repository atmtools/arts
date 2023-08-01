#include "workspace_agenda_class.h"

#include <auto_wsa.h>
#include <auto_wsm.h>

#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include <string>

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

void Agenda::finalize() try {
  static const auto& wsa = workspace_agendas();

  checked = false;
  output.clear();
  pure_input.clear();

  for (auto& method : methods) {
    const std::vector<std::string>& m_out = method.get_outs();
    const std::vector<std::string>& m_in = method.get_ins();

    if (auto& v = method.get_setval(); v) {
      output.push_back(method.get_name());
    }

    output.insert(output.end(), m_out.begin(), m_out.end());
    pure_input.insert(pure_input.end(), m_in.begin(), m_in.end());
  }

  std::sort(output.begin(), output.end());
  output.erase(std::unique(output.begin(), output.end()), output.end());

  std::sort(pure_input.begin(), pure_input.end());
  pure_input.erase(std::unique(pure_input.begin(), pure_input.end()),
                   pure_input.end());

  /*
  Named agendas are different from unnamed agendas

  The name of an agenda indicates that a set of inputs
  and outputs are expected to be present in the workspace
  or will be added to the workspace after-the-fact in the
  corresponding agenda execute method.  They must thus be
  named by at least one of the methods that make up the
  agenda execution.  Here we find and remove them from the
  callstack as their presence would otherwise be an error.
  */
  if (auto ptr = wsa.find(name); ptr not_eq wsa.end()) {
    const std::vector<std::string>& ag_inout = ptr->second.inout;
    const std::vector<std::string>& ag_out = ptr->second.output;
    const std::vector<std::string>& ag_in = ptr->second.input;

    for (auto& item : ag_inout) {
      if ((std::erase(pure_input, item) + std::erase(output, item)) not_eq 2) {
        throw std::runtime_error(var_string(
            "Expected input and output ", std::quoted(item), " not found"));
      }
    }

    for (auto& item : ag_in) {
      if (std::erase(pure_input, item) not_eq 1) {
        throw std::runtime_error(
            var_string("Expected input ", std::quoted(item), " not found"));
      }
      std::erase(output, item);
    }

    for (auto& item : ag_out) {
      std::erase(pure_input, item);
      if (std::erase(output, item) not_eq 1) {
        throw std::runtime_error(
            var_string("Expected output ", std::quoted(item), " not found"));
      }
    }
  }

  // Finally, remove duplicates from output that are also in pure_input
  for (auto& item : output) {
    std::erase(pure_input, item);
  }

  output.shrink_to_fit();
  pure_input.shrink_to_fit();

  checked = true;
} catch (std::exception& e) {
  throw std::runtime_error(var_string(
      "Error finalizing agenda ", std::quoted(name), "\n", e.what()));
}

void Agenda::copy_workspace(Workspace& out, const Workspace& in) const try {
  for (auto& str : output) {
    out.set(str, in.copy_type(str));
  }

  for (auto& str : pure_input) {
    out.set(str, in.share(str));
  }
} catch (std::exception& e) {
  throw std::runtime_error(var_string(
      "Error setting up agenda ", std::quoted(name), "\n", e.what()));
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
