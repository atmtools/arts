#ifndef python_interface_h
#define python_interface_h

#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include <filesystem>
#include <iomanip>
#include <iostream>

#include "debug.h"
#include "auto_md.h"
#include "enums.h"

//! Contains a bunch of helper functions to manipulate python objects inside C++
namespace Python {
namespace py = pybind11;

struct Numeric_ {
  Numeric val;
  constexpr operator Numeric&() { return val; }
  constexpr operator const Numeric&() const { return val; }
};

struct Index_ {
  Index val;
  constexpr operator Index&() { return val; }
  constexpr operator const Index&() const { return val; };
};

template<typename OUT>
OUT& refcast(const py::object& in) {
  return *in.cast<OUT*>();
}

template<>
inline Index& refcast<Index>(const py::object& in) {
  return in.cast<Index_*>() -> val;
}

template<>
inline Numeric& refcast<Numeric>(const py::object& in) {
  return in.cast<Numeric_*>() -> val;
}

template <class T, class Tptr>
T& select(T& default_,
          Tptr* val,
          py::kwargs& kwargs,
          const char* const name) {
  if (val not_eq nullptr) return *val;
  if (kwargs.contains(name)) return refcast<T>(kwargs[name]);
  return default_;
}

template <class T, class Tptr>
T& select(Workspace& ws,
          size_t iws,
          Tptr* val,
          py::kwargs& kwargs,
          const char* const name) {
  if (val not_eq nullptr) return *val;
  if (kwargs.contains(name)) return refcast<T>(kwargs[name]);
  ARTS_USER_ERROR_IF(not ws.is_initialized(iws),
                     Workspace::wsv_data[iws].Name(),
                     " is not initialized")
  return *reinterpret_cast<T*>(ws[iws]);
}

template <class T, class Tptr>
T& select(Tptr* val,
          py::kwargs& kwargs,
          const char* const name) {
  if (val not_eq nullptr) return *val;
  if (kwargs.contains(name)) return refcast<T>(kwargs[name]);
  ARTS_USER_ERROR("User generated data is not defined ", name)
}
}  // namespace Python

#endif  // python_interface_h
