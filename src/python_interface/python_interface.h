#ifndef python_interface_h
#define python_interface_h

#include <auto_md.h>
#include <debug.h>
#include <enums.h>
#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>
#include <pybind11/stl_bind.h>

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <variant>

#include <auto_md.h>
#include <global_data.h>
#include <m_append.h>
#include <m_basic_types.h>
#include <m_conversion.h>
#include <m_copy.h>
#include <m_delete.h>
#include <m_extract.h>
#include <m_general.h>
#include <m_gridded_fields.h>
#include <m_ignore.h>
#include <m_nc.h>
#include <m_reduce.h>
#include <m_select.h>
#include <m_xml.h>
#include <parameters.h>
#include <supergeneric.h>
#include <xml_io.h>

#include <py_auto_interface.h>

#include "python_interface_value_type.h"

//! Contains a bunch of helper functions to manipulate python objects inside C++
namespace Python {
namespace py = pybind11;

//! Only for debugs (and to suppress warning about the line above)
template <class... Ts> void print(const Ts&... args) {py::print(args...);}

template <class T, class VariantT>
T& select_gout(VariantT& val) {
  return std::visit([](auto&& out) -> T& { return *out; }, val);
}

template <class T, class VariantT>
const T& select_gin(const VariantT& val) {
  return std::visit([](auto&& out) -> const T& { return *out; }, val);
}

template <class T, class VariantT>
const T& select_gin(const T& default_, const std::optional<VariantT>& val) {
  if (val)
    return std::visit([](auto&& out) -> const T& { return *out; }, val.value());
  return default_;
}

template <class T, class WorkspaceVariable, class VariantT>
T& select_out(WorkspaceVariable wsv, std::optional<VariantT>& val) {
  if (val)
    return std::visit([](auto&& out) -> T& { return *out; }, val.value());
  return wsv;
}

template <class T, class WorkspaceVariable, class VariantT>
T& select_inout(const WorkspaceVariable wsv, std::optional<VariantT>& val) {
  if (val)
    return std::visit([](auto&& out) -> T& { return *out; }, val.value());
  return wsv;
}

template <class T, class WorkspaceVariable, class VariantT>
const T& select_in(const WorkspaceVariable wsv,
                   const std::optional<VariantT>& val) {
  if (val)
    return std::visit([](auto&& out) -> const T& { return *out; }, val.value());
  return wsv;
}

template <std::size_t i, class SelectT, class VariantT>
SelectT select_internal(VariantT val) {
  constexpr std::size_t n = std::variant_size_v<VariantT>;

  if constexpr (i + 1 == n) {
    return std::get<i>(val);
  } else {
    if (val.index() == i) return std::get<i>(val);
    return select_internal<i + 1, SelectT>(val);
  }
}

template <class SelectT, class VariantT>
SelectT select_wvv(VariantT val) {
  constexpr std::size_t n = std::variant_size_v<VariantT>;
  static_assert(n > 1,
                "Cannot select a workspace variable variant from 0 arguments");

  if (val.index() == 0) return *std::get<0>(val);
  return select_internal<1, SelectT>(val);
}
}  // namespace Python

#endif  // python_interface_h
