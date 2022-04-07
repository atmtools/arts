#ifndef python_interface_h
#define python_interface_h

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <variant>

#include <pybind11/pybind11.h>
#include <py_auto_interface.h>

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>
#include <pybind11/stl_bind.h>

#include <auto_md.h>
#include <debug.h>
#include <enums.h>
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


//! Contains a bunch of helper functions to manipulate python objects inside C++
namespace Python {
namespace py = pybind11;

//! Only for debugs (and to suppress warning about the line above)
template <class... Ts> void print(const Ts&... args) {py::print(args...);}

template <class T, class U>
T& select_gout(std::variant<WorkspaceVariable *, U *>& val) {
  return std::visit([](auto&& out) -> T& { return *out; }, val);
}

template <class T, class U>
const T& select_gin(const std::variant<const WorkspaceVariable *, U *>& val) {
  return std::visit([](auto&& out) -> const T& { return *out; }, val);
}

template <class T, class U>
const T& select_gin(const T& default_, const std::optional<std::variant<const WorkspaceVariable *, U *>>& val) {
  if (val)
    return std::visit([](auto&& out) -> const T& { return *out; }, val.value());
  return default_;
}

template <class T, class U>
T& select_out(WorkspaceVariable wsv, std::optional<std::variant<WorkspaceVariable *, U *>>& val) {
  if (val)
    return std::visit([](auto&& out) -> T& { return *out; }, val.value());
  return wsv;
}

template <class T, class U>
T& select_inout(const WorkspaceVariable wsv, std::optional<std::variant<const WorkspaceVariable *, U *>>& val) {
  if (val)
    return std::visit([](auto&& out) -> T& { return *out; }, val.value());
  return wsv;
}

template <class T, class U>
const T& select_in(const WorkspaceVariable wsv,
                   const std::optional<std::variant<const WorkspaceVariable *, U *>>& val) {
  if (val)
    return std::visit([](auto&& out) -> const T& { return *out; }, val.value());
  return wsv;
}

template <class... Types>
WorkspaceVariablesVariant select_wvv(std::variant<WorkspaceVariable *, Types...>&  val) {
  return std::visit([](auto&& x) -> WorkspaceVariablesVariant {
    if constexpr (std::is_convertible_v<decltype(x), WorkspaceVariable *>) return *x;
    else return x;
    }, val);
}

template <class... Types>
WorkspaceVariablesVariant select_wvv(std::variant<const WorkspaceVariable *, Types...>&  val) {
  return std::visit([](auto&& x) -> WorkspaceVariablesVariant {
    if constexpr (std::is_convertible_v<decltype(x), const WorkspaceVariable *>) return *x;
    else return x;
    }, val);
}
}  // namespace Python

#endif  // python_interface_h
