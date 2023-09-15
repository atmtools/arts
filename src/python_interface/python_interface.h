#ifndef python_interface_h
#define python_interface_h

#include <py_auto_wsg.h>
#include <pybind11/attr.h>
#include <pybind11/cast.h>
#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/detail/common.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>
#include <pybind11/stl_bind.h>
#include <workspace.h>

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <variant>

#include "auto_wsg.h"

// Exposed types whose base type are not workspace variables
PYBIND11_MAKE_OPAQUE(Array<SpeciesTag>);
PYBIND11_MAKE_OPAQUE(Array<LagrangeInterpolation>);
PYBIND11_MAKE_OPAQUE(Array<Species::Species>);
PYBIND11_MAKE_OPAQUE(Array<Species::IsotopeRecord>);
PYBIND11_MAKE_OPAQUE(Array<AbsorptionSingleLine>);
PYBIND11_MAKE_OPAQUE(Array<SpeciesErrorCorrectedSuddenData>);
PYBIND11_MAKE_OPAQUE(Array<ErrorCorrectedSuddenData>);

//! Contains a bunch of helper functions to manipulate python objects inside C++
namespace Python {
namespace py = pybind11;

template <typename T, typename... Ts>
using artsclass = py::class_<T, Ts..., std::shared_ptr<T>>;

ENUMCLASS(ArrayOptions, char, index, nolist, noobjlist, nopickle, nocopy)

template <typename PythonListable, ArrayOptions... opts>
py::class_<PythonListable, std::shared_ptr<PythonListable>> artsarray(
    py::handle scope, const std::string& name, auto&&... args) {
  using base = std::remove_cvref_t<decltype(PythonListable{}.front())>;

  auto out = py::bind_vector<PythonListable, std::shared_ptr<PythonListable>>(
      scope, name, std::forward<decltype(args)>(args)...);

  if constexpr (not((opts == ArrayOptions::nocopy) or ...)) {
    out.def("__copy__", [](const PythonListable& v) {
      return std::make_shared<PythonListable>(v);
    });

    out.def("__deepcopy__", [](const PythonListable& v, const py::dict&) {
      return std::make_shared<PythonListable>(v);
    });
  }

  if constexpr (not((opts == ArrayOptions::noobjlist) or ...)) {
    out.def(py::init([](const Array<py::object>& x) {
      auto v = std::make_shared<PythonListable>(x.size());
      std::transform(x.begin(), x.end(), v->begin(), [](const auto& s) {
        return py::cast<base>(s);
      });
      return v;
    }));
    py::implicitly_convertible<Array<py::object>, PythonListable>();
  }

  if constexpr (not((opts == ArrayOptions::nolist) or ...)) {
    out.def(py::init([](const py::list& x) {
      auto v = std::make_shared<PythonListable>(x.size());
      std::transform(x.begin(), x.end(), v->begin(), [](const auto& s) {
        return py::cast<base>(s);
      });
      return v;
    }));
    py::implicitly_convertible<py::list, PythonListable>();
  }

  if constexpr (not((opts == ArrayOptions::nopickle) or ...)) {
    out.def(
        py::pickle([](const PythonListable& t) { return py::make_tuple(t); },
                   [](const py::tuple& t) {
                     ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
                     return t[0].cast<PythonListable>();
                   }));
  }

  if constexpr (((opts == ArrayOptions::index) or ...)) {
    out.def("index", [](const PythonListable& v, const base& f) {
      auto ptr = std::ranges::find(v, f);
      if (ptr == v.end())
        throw py::value_error(var_string(f, " is not in list"));
      return std::distance(v.begin(), ptr);
    });
  }

  return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T, PythonWorkspaceGroup U>
T& select_out(std::optional<std::shared_ptr<U>>& x,
              Workspace& ws,
              const char* const name) {
  return x ? *x.value() : ws.get_or<T>(name);
}

template <WorkspaceGroup T>
T& select_out(std::optional<ValueHolder<T>>& x,
              Workspace& ws,
              const char* const name) {
  return x ? static_cast<T&>(x.value()) : ws.get_or<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T, PythonWorkspaceGroup U>
T& select_gout(std::optional<std::shared_ptr<U>>& x, const char* const name) {
  return x ? *x.value()
           : throw std::runtime_error(
                 var_string("Unknown ouput: ", std::quoted(name)));
}

template <WorkspaceGroup T>
T& select_gout(std::optional<ValueHolder<T>>& x, const char* const name) {
  return x ? x.value()
           : throw std::runtime_error(
                 var_string("Unknown ouput: ", std::quoted(name)));
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T, PythonWorkspaceGroup U>
T& select_inout(std::optional<std::shared_ptr<U>>& x,
                const Workspace& ws,
                const char* const name) {
  return x ? *x.value() : ws.get<T>(name);
}

template <WorkspaceGroup T>
T& select_inout(std::optional<ValueHolder<T>>& x,
                const Workspace& ws,
                const char* const name) {
  return x ? static_cast<T&>(x.value()) : ws.get<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T, PythonWorkspaceGroup U>
const T& select_in(const std::optional<const std::shared_ptr<U>>& x,
                   const Workspace& ws,
                   const char* const name) {
  return x ? *x.value() : ws.get<T>(name);
}

template <WorkspaceGroup T>
const T& select_in(const std::optional<const ValueHolder<T>>& x,
                   const Workspace& ws,
                   const char* const name) {
  return x ? static_cast<const T&>(x.value()) : ws.get<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T, PythonWorkspaceGroup U>
const T& select_gin(const std::optional<const std::shared_ptr<U>>& x,
                    const char* const name) {
  return x ? *x.value()
           : throw std::runtime_error(
                 var_string("Unknown input: ", std::quoted(name)));
}

template <WorkspaceGroup T>
const T& select_gin(const std::optional<const ValueHolder<T>>& x,
                    const char* const name) {
  return x ? x.value()
           : throw std::runtime_error(
                 var_string("Unknown input: ", std::quoted(name)));
}

template <WorkspaceGroup T, PythonWorkspaceGroup U>
const T& select_gin(const std::optional<const std::shared_ptr<U>>& x,
                    const T& defval) {
  return x ? *x.value() : defval;
}

template <WorkspaceGroup T>
const T& select_gin(const std::optional<const ValueHolder<T>>& x,
                    const T& defval) {
  return x ? static_cast<const T&>(x.value()) : defval;
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
}  // namespace Python

#endif  // python_interface_h
