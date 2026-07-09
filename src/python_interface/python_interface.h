#ifndef python_interface_h
#define python_interface_h

#include <nanobind/nanobind.h>
#include <py_auto_options.h>
#include <py_auto_wsg.h>
#include <workspace.h>

#include <memory>
#include <variant>

#include "hpy_opaque.h"

using ssize_t = Py_ssize_t;

//! Contains a bunch of helper functions to manipulate python objects inside C++
namespace Python {
namespace py = nanobind;
using namespace py::literals;

//! A ``std::variant<std::shared_ptr<T>...>`` used for supergeneric
//! workspace-method arguments, with a custom nanobind type caster (defined
//! below) that preserves the implicit-conversion (``convert``) flag when
//! probing each alternative.
//!
//! This works around a nanobind 2.13 regression in
//! ``type_caster<std::shared_ptr<T>>::from_python`` that strips the
//! ``convert`` flag, which silently disabled implicit conversions such as
//! ``"H2O"`` -> ``SpeciesEnum`` registered via
//! ``py::implicitly_convertible<std::string, SpeciesEnum>``.  As a result,
//! passing a ``str`` for a supergeneric ``key`` argument stopped matching any
//! variant alternative.
//!
//! ``Supergeneric<...>`` is layout-compatible with the underlying variant and
//! exposes it via ``operator const std::variant<...>&`` so the existing
//! ``select_gin`` / ``select_gout`` / ``has_selected_value`` overloads (which
//! take ``const std::variant<std::shared_ptr<T>...>&``) keep working without
//! modification.
template <typename... T>
struct Supergeneric {
  std::variant<std::shared_ptr<T>...> value;

  operator const std::variant<std::shared_ptr<T>...>&() const { return value; }
};

template <typename... T>
bool has_selected_value(const std::variant<std::shared_ptr<T>...>& x) {
  return std::visit([](const auto& ptr) { return static_cast<bool>(ptr); }, x);
}

template <typename... T>
bool has_selected_value(const Supergeneric<T...>& x) {
  return has_selected_value(x.value);
}

template <typename... T>
std::string type(const Supergeneric<T...>* const x) {
  return type(x ? &x->value : nullptr);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
T& select_out(T* const x, Workspace& ws, const char* const name) {
  return x ? *x : ws.get_or<T>(name);
}

template <WorkspaceGroup T>
T& select_out(ValueHolder<T>* const x, Workspace& ws, const char* const name) {
  return x ? static_cast<T&>(*x) : ws.get_or<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
T& select_gout(T* const x, Workspace& ws, const char* const name) {
  return x ? *x : ws.get_or<T>(name);
}

template <WorkspaceGroup T>
T& select_gout(ValueHolder<T>* const x, Workspace& ws, const char* const name) {
  return x ? static_cast<T&>(*x) : ws.get_or<T>(name);
}

template <typename... T>
std::variant<std::shared_ptr<T>...>& select_gout(
    std::variant<std::shared_ptr<T>...>* const x,
    Workspace&,
    const char* const name) {
  if (not x or not has_selected_value(*x)) {
    throw std::runtime_error(std::format("Unknown output: \"{}\"", name));
  }

  return *x;
}

template <typename... T>
std::variant<std::shared_ptr<T>...>& select_gout(Supergeneric<T...>* const x,
                                                 Workspace& ws,
                                                 const char* const name) {
  return select_gout(x ? &x->value : nullptr, ws, name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
T& select_inout(T* const x, const Workspace& ws, const char* const name) {
  return x ? *x : ws.get<T>(name);
}

template <WorkspaceGroup T>
T& select_inout(ValueHolder<T>* const x,
                const Workspace& ws,
                const char* const name) {
  return x ? static_cast<T&>(*x) : ws.get<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
const T& select_in(const T* const x,
                   const Workspace& ws,
                   const char* const name) {
  return x ? *x : ws.get<T>(name);
}

template <WorkspaceGroup T>
const T& select_in(const ValueHolder<T>* const x,
                   const Workspace& ws,
                   const char* const name) {
  return x ? static_cast<const T&>(*x) : ws.get<T>(name);
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

template <WorkspaceGroup T>
const T& select_gin(const T* const x, const char* const name) {
  return x ? *x
           : throw std::runtime_error(
                 std::format("Unknown input: \"{}\"", name));
}

template <WorkspaceGroup T>
const T& select_gin(const ValueHolder<T>* const x, const char* const name) {
  return x ? static_cast<const T&>(*x)
           : throw std::runtime_error(
                 std::format("Unknown input: \"{}\"", name));
}

template <typename... T>
const std::variant<std::shared_ptr<T>...>& select_gin(
    const std::variant<std::shared_ptr<T>...>* const x,
    const char* const name) {
  if (not x or not has_selected_value(*x)) {
    throw std::runtime_error(std::format("Unknown input: \"{}\"", name));
  }

  return *x;
}

template <typename... T>
const std::variant<std::shared_ptr<T>...>& select_gin(
    const Supergeneric<T...>* const x, const char* const name) {
  return select_gin(x ? &x->value : nullptr, name);
}

template <WorkspaceGroup T>
const T& select_gin(const T* const x, const T& defval) {
  return x ? *x : defval;
}

template <WorkspaceGroup T>
const T& select_gin(const ValueHolder<T>* const x, const T& defval) {
  return x ? static_cast<const T&>(*x) : defval;
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
}  // namespace Python

template <typename T>
struct std::hash<Python::ValueHolder<T>> {
  static std::size_t operator()(const Python::ValueHolder<T>& x) {
    return std::hash<T>{}(x->val);
  }
};

NAMESPACE_BEGIN(NB_NAMESPACE)
NAMESPACE_BEGIN(detail)

//! Custom type caster for ``Python::Supergeneric<T...>``.
//!
//! It mirrors nanobind's stock ``std::variant<...>`` caster (see
//! ``nanobind/stl/variant.h``) but, crucially, forwards the ``convert`` flag
//! intact to each alternative's caster.  nanobind 2.13 added
//! ``flags &= ~cast_flags::convert`` at the top of
//! ``type_caster<std::shared_ptr<T>>::from_python`` (to fix a use-after-free
//! with implicitly-converted temporaries), which as a side effect silently
//! disabled implicit conversions such as ``str`` -> ``SpeciesEnum`` whenever
//! the alternative is wrapped in ``std::shared_ptr<...>``.  ARTS exposes
//! supergeneric workspace-method arguments as
//! ``const std::variant<std::shared_ptr<pyT1>, ...> * const``, so passing
//! ``key="H2O"`` stopped matching any alternative after the 2.13 upgrade.
//!
//! By taking ownership of the conversion here, the ``convert`` flag reaches
//! the inner ``T`` caster (e.g. ``SpeciesEnum``), which then honors
//! ``py::implicitly_convertible<std::string, SpeciesEnum>`` as before.  On
//! success the produced ``T*`` (owned by the cleanup list for the duration of
//! the call) is wrapped in a ``std::shared_ptr<T>`` whose deleter decrefs the
//! Python source object, matching the lifetime semantics of the stock
//! shared_ptr caster.
template <typename... T>
struct type_caster<Python::Supergeneric<T...>> {
  using Value   = Python::Supergeneric<T...>;
  using Variant = std::variant<std::shared_ptr<T>...>;

  static constexpr auto Name =
      union_name(make_caster<std::shared_ptr<T>>::Name...);

  template <typename T_>
  using Cast = movable_cast_t<T_>;
  template <typename T_>
  static constexpr bool can_cast() {
    return true;
  }

  explicit operator Value*() { return &value; }
  explicit operator Value&() { return (Value&)value; }
  explicit operator Value&&() { return (Value&&)value; }

  bool from_python(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
    // Same two-pass strategy as nanobind's stock variant caster: first without
    // convert, then with convert.  The difference is that the per-alternative
    // probe below bypasses type_caster<std::shared_ptr<T>> (which strips the
    // convert flag on nanobind 2.13+) and invokes the inner T caster
    // directly, preserving convert so implicit conversions succeed.
    if ((try_alt<T>(src, flags & ~(uint8_t)cast_flags::convert, cleanup) ||
         ...) ||
        (try_alt<T>(src, flags, cleanup) || ...))
      return true;
    return false;
  }

 private:
  template <typename Alt>
  bool try_alt(handle src, uint8_t flags, cleanup_list* cleanup) noexcept {
    using CasterT = make_caster<Alt>;
    CasterT caster;
    if (!caster.from_python(src, flags, cleanup)) return false;
    Alt* ptr = caster.operator Alt*();
    if (!ptr) return false;
    // Keep the Python source alive for the lifetime of the shared_ptr, the
    // same way the stock shared_ptr caster does (see py_deleter in
    // nanobind/stl/shared_ptr.h).
    src.inc_ref();
    value.value.template emplace<std::shared_ptr<Alt>>(
        std::shared_ptr<Alt>(ptr, [src = src.ptr()](void*) noexcept {
          if (is_alive()) {
            gil_scoped_acquire guard;
            Py_DECREF(src);
          }
        }));
    return true;
  }

 public:
  static handle from_cpp(const Value& v,
                         rv_policy policy,
                         cleanup_list* cleanup) noexcept {
    return make_caster<Variant>::from_cpp(v.value, policy, cleanup);
  }

  static handle from_cpp(const Value* v,
                         rv_policy policy,
                         cleanup_list* cleanup) noexcept {
    if (!v) return none().release();
    return from_cpp(*v, policy, cleanup);
  }

  Value value;
};

NAMESPACE_END(detail)
NAMESPACE_END(NB_NAMESPACE)

#endif  // python_interface_h
