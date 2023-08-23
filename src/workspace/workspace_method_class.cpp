#include "workspace_agenda_class.h"
#include "workspace_method_class.h"

#include "workspace_class.h"

#include <auto_wsm.h>

#include <algorithm>
#include <exception>
#include <iomanip>
#include <stdexcept>
#include <string_view>

const auto& wsms = workspace_methods();

std::ostream& operator<<(std::ostream& os, const Method& m) {
  os << m.name;

  if (m.setval) {
    std::string var = std::visit([](auto& v) { return var_string(*v); },
                                 m.setval.value().value);
    constexpr std::size_t maxsize = 50;
    if (var.size() > maxsize) {
      var = std::string(var.begin(), var.begin() + maxsize) + "...";
    }
    std::replace(var.begin(), var.end(), '\n', ' ');

    os << " = " << var;
  } else {
    os << '(';

    os << "o : ";
    bool first = true;
    for (auto& o : m.outargs) {
      if (not first) os << ", ";
      first = false;
      os << o;
    }

    os << '\n' << std::string(m.name.size() + 1, ' ') << "i : ";
    first = true;
    for (auto& o : m.inargs) {
      if (not first) os << ", ";
      first = false;
      os << o;
    }
    os << ')';
  }

  return os;
}

Method::Method(const std::string& n,
               const std::vector<std::string>& a,
               const std::unordered_map<std::string, std::string>& kw) try
    : name(n), outargs(wsms.at(name).out), inargs(wsms.at(name).in) {
  const std::size_t nargout = outargs.size();
  const std::size_t nargin = inargs.size();

  // FIXME: IN C++23, USE ZIP HERE INSTEAD AS WE CAN REMOVE LATER CODE DOING THAT
  std::vector<std::pair<std::string, bool>> outargs_set(nargout);
  std::vector<std::pair<std::string, bool>> inargs_set(nargin);
  for (std::size_t i=0; i<nargout; ++i) outargs_set[i] = {outargs[i], false};
  for (std::size_t i=0; i<nargin; ++i) inargs_set[i] = {inargs[i], false};

  // Common filter
  const auto unset =
      std::views::filter([](const auto& p) { return not p.second; });
  
  // Common G-name
  const auto is_gname = [](const auto& str1, auto& str2) {
    return str1.front() == '_' and
           std::string_view(str1.begin() + 1, str1.end()) == str2;
  };

  // Positional arguments
  {
    const auto fuzzy_equals = [is_gname](auto& arg) {
      return std::views::filter([&arg, is_gname](const auto& p) {
        return p.first == arg or is_gname(p.first, arg);
      });
    };

    std::transform(a.begin(),
                   a.begin() + std::min(nargout, a.size()),
                   outargs.begin(),
                   outargs_set.begin(),
                   [&](auto& arg, auto& orig) {
                     if (auto x = inargs_set | unset | fuzzy_equals(orig); bool(x)) {
                       auto ptr = x.begin();
                       ptr->first = arg;
                       ptr->second = true;
                     }
                     return std::pair<std::string, bool>{arg, true};
                   });

    auto unset_filt = inargs_set | unset;
    std::transform(a.begin() + std::min(nargout, a.size()),
                   a.end(),
                   unset_filt.begin(),
                   [](auto& arg) {
                     return std::pair<std::string, bool>{arg, true};
                   });
  }

  // Named arguments
  {
    for (auto& [key, val] : kw) {
      bool any=false;

      for (auto& arg : outargs_set | unset) {
        if (arg.first == key) {
          arg.first = val;
          arg.second = true;
          any = true;
        }

        if (is_gname(arg.first, key)) {
          arg.first = val;
          arg.second = true;
          any = true;
        }
      }

      for (auto& arg : inargs_set | unset) {
        if (arg.first == key) {
          arg.first = val;
          arg.second = true;
          any = true;
        }

        if (is_gname(arg.first, key)) {
          arg.first = val;
          arg.second = true;
          any = true;
        }
      }

      if (not any) {
        throw std::runtime_error(var_string("No named argument ", std::quoted(key)));
      }
    }
  }

  // FIXME: REMOVE THESE TWO IN C++23 WITH ZIP
  for (std::size_t i=0; i<nargout; ++i) {
    if (outargs_set[i].second) {
      outargs[i] = outargs_set[i].first;
    }
  }
  for (std::size_t i=0; i<nargin; ++i) {
    if (inargs_set[i].second) {
      inargs[i] = inargs_set[i].first;
    }
  }

  // Check that all non-defaulted GINS are set
  for (std::size_t i=0; i<nargin; i++) {
    if (inargs[i].front() == '_' and not wsms.at(n).defs.contains(inargs[i])) {
      throw std::runtime_error(
          var_string("Missing required generic input argument ",
                     std::quoted(std::string_view(inargs[i].begin() + 1,
                                                  inargs[i].end()))));
    }
  }
} catch (std::out_of_range&) {
  throw std::runtime_error(var_string("No method named ", std::quoted(n)));
} catch (std::exception& e) {
  throw std::runtime_error(var_string(
      "Error in method construction for ", std::quoted(n), "\n", e.what()));
}

Method::Method(std::string n, const Wsv& wsv, bool overwrite)
    : name(std::move(n)), setval(wsv), overwrite_setval(overwrite) {
  if (not setval) {
    throw std::runtime_error(var_string("Cannot set workspace variable ",
                                        std::quoted(name),
                                        " to empty value"));
  }
}

void Method::operator()(Workspace& ws) const try {
  if (setval) {
    if (overwrite_setval)
      ws.overwrite(name, std::make_shared<Wsv>(setval.value()));
    else
      ws.set(name, std::make_shared<Wsv>(setval.value()));
  } else {
    wsms.at(name).func(ws, outargs, inargs);
  }
} catch (std::out_of_range&) {
  throw std::runtime_error(var_string("No method named ", std::quoted(name)));
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("Error in method ", *this, "\n", e.what()));
}

void Method::add_defaults_to_agenda(Agenda& agenda) const {
  if (not setval) {
    const auto& map = wsms.at(name).defs; 
    for (auto& arg : inargs) {
      if (arg.front() == '_' and map.contains(arg)) {
       agenda.add(Method{arg, map.at(arg), true});
      }
    }
  }
}

Method::Method(const std::string& n,
               const std::vector<std::string>& ins,
               const std::vector<std::string>& outs,
               const std::optional<Wsv>& wsv,
               bool overwrite)
    : name(n),
      outargs(outs),
      inargs(ins),
      setval(wsv),
      overwrite_setval(overwrite) {}
