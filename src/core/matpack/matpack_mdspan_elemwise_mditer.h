#pragma once

#include <configtypes.h>

#include <ranges>
#include <type_traits>

#include "matpack_mdspan_common_types.h"

namespace matpack {
template <class T>
struct elemwise_mditer {
  T* data{nullptr};
  Index pos{0};

  using difference_type = Index;
  using value_type      = std::remove_cvref_t<decltype(data->elem_at(0))>;

  constexpr elemwise_mditer()                                      = default;
  constexpr elemwise_mditer(elemwise_mditer&&) noexcept            = default;
  constexpr elemwise_mditer(const elemwise_mditer&)                = default;
  constexpr elemwise_mditer& operator=(elemwise_mditer&&) noexcept = default;
  constexpr elemwise_mditer& operator=(const elemwise_mditer&)     = default;

  constexpr elemwise_mditer(T& x) : data(&x) {}

  constexpr elemwise_mditer& operator++() noexcept {
    pos++;
    return *this;
  }
  constexpr elemwise_mditer& operator--() noexcept {
    pos--;
    return *this;
  }
  constexpr elemwise_mditer operator++(int) noexcept {
    elemwise_mditer out(*this);
    ++pos;
    return out;
  }
  constexpr elemwise_mditer operator--(int) noexcept {
    elemwise_mditer out(*this);
    --pos;
    return out;
  }

  constexpr elemwise_mditer& operator+=(Index i) noexcept {
    pos += i;
    return *this;
  }
  constexpr elemwise_mditer& operator-=(Index i) noexcept {
    pos -= i;
    return *this;
  }
  [[nodiscard]] constexpr elemwise_mditer operator+(Index i) const noexcept {
    elemwise_mditer out(*this);
    out.pos += i;
    return out;
  }
  [[nodiscard]] constexpr elemwise_mditer operator-(Index i) const noexcept {
    elemwise_mditer out(*this);
    out.pos -= i;
    return out;
  }
  [[nodiscard]] constexpr friend elemwise_mditer operator+(
      Index i, const elemwise_mditer& m) noexcept {
    return m + i;
  }
  [[nodiscard]] constexpr friend elemwise_mditer operator-(
      Index i, const elemwise_mditer& m) noexcept {
    return m - i;
  }

  [[nodiscard]] constexpr difference_type operator-(
      const elemwise_mditer& other) const noexcept {
    return pos - other.pos;
  }
  [[nodiscard]] constexpr difference_type operator+(
      const elemwise_mditer& other) const noexcept {
    return pos + other.pos;
  }

  [[nodiscard]] constexpr auto operator<=>(
      const elemwise_mditer&) const noexcept = default;

  [[nodiscard]] constexpr decltype(auto) operator*() const {
    return data->elem_at(pos);
  }

  [[nodiscard]] constexpr decltype(auto) operator[](Index i) const {
    return data->elem_at(pos + i);
  }
};

template <any_md T>
auto elemwise_range(T&& x) {
  return stdr::subrange{x.elem_begin(), x.elem_end()};
};

// Catch for non-matpack types (might not work, be careful)
template <typename T>
auto elemwise_range(T&& x)
  requires(rank<T>() == 1 and not any_md<T>)
{
  return stdr::subrange{x.begin(), x.end()};
};
}  // namespace matpack
