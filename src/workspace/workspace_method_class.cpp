#include "workspace_agenda_class.h"
#include "workspace_method_class.h"

#include "workspace_class.h"
#include <auto_wsm.h>

#include <algorithm>
#include <exception>
#include <iomanip>
#include <stdexcept>

const auto& wsms = workspace_methods();

std::ostream& operator<<(std::ostream& os, const Method& m) {
  if (m.setval) {
    os << "setting " << std::quoted(m.name);
  } else {
    os << "calling " << std::quoted(m.name);
  }
  return os;
}

Method::Method(const std::string& n,
               const std::vector<std::string>& a,
               const std::unordered_map<std::string, std::string>& kw) try
    : name(n), outargs(wsms.at(name).out), inargs(wsms.at(name).in), func(wsms.at(name).func) {
  const std::size_t nargout = outargs.size();
  const std::size_t nargin = [&]() {
    std::size_t num = 0;
    for (const auto& arg : inargs) {
      num += (std::find(outargs.begin(), outargs.end(), arg) == outargs.end());
    }
    return num;
  }();

  if ((a.size() + kw.size()) > (nargin + nargout)) {
    throw std::runtime_error(var_string("Too many arguments to method ",
                                        std::quoted(n),
                                        ".  At most ",
                                        nargin + nargout,
                                        " arguments are accepted, but got ",
                                        a.size() + kw.size(),
                                        " arguments"));
  }

  if (nargout < a.size()) {
    std::copy(a.begin(), a.begin() + nargout, outargs.begin());
    std::copy(a.begin() + nargout, a.end(), inargs.begin());
  } else {
    std::copy(a.begin(), a.end(), outargs.begin());
  }

  for (auto& [k, v] : kw) {
    auto out = std::find(outargs.begin(), outargs.end(), k);
    auto in = std::find(inargs.begin(), inargs.end(), k);

    if (out not_eq outargs.end()) {
      *out = v;
    }

    if (in not_eq inargs.end()) {
      *in = v;
    }

    if (out == outargs.end() and in == inargs.end()) {
      throw std::runtime_error(
          var_string("Keyword argument ", std::quoted(k), " not found."));
    }
  }
} catch (std::out_of_range&) {
  throw std::runtime_error(var_string("No method named ", std::quoted(n)));
} catch (std::exception& e) {
  throw std::runtime_error(var_string(
      "Error in method construction for ", std::quoted(n), "\n", e.what()));
}

Method::Method(std::string n, const Wsv& wsv)
    : name(std::move(n)), setval(wsv) {
  if (not setval) {
    throw std::runtime_error(var_string("Cannot set workspace variable ",
                                        std::quoted(name),
                                        " to empty value"));
  }
}

void Method::operator()(Workspace& ws) const try {
  if (setval) {
    ws.set(name, std::make_shared<Wsv>(setval.value()));
  } else {
    func(ws, outargs, inargs);
  }
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("Error in method ", *this, "\n", e.what()));
}


void Method::add_setvals(Agenda&) const {
  
}
