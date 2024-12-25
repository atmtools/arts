#pragma once

#include <configtypes.h>

#include <concepts>
#include <iterator>
#include <ranges>
#include <type_traits>

namespace matpack {
template <class T>
struct elemwise_mditer {
  const T* data{nullptr};
  Index pos{0};

  using difference_type = Index;
  using value_type      = std::remove_cvref_t<decltype(data->elem_at(0))>;

  constexpr elemwise_mditer()                                      = default;
  constexpr elemwise_mditer(elemwise_mditer&&) noexcept            = default;
  constexpr elemwise_mditer(const elemwise_mditer&)                = default;
  constexpr elemwise_mditer& operator=(elemwise_mditer&&) noexcept = default;
  constexpr elemwise_mditer& operator=(const elemwise_mditer&)     = default;

  constexpr elemwise_mditer(const T& x) : data(&x) {}
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

template <typename T>
concept elemwise_mditerable = requires(T a) {
  a.elem_at(Index{});
  a.elem_begin();
  a.elem_end();
};

// Catch-all for non-iterable types -- fallback, if you are here you are probably doing something weird
template <typename T>
auto elemwise_range(T&& x) {
  return std::ranges::subrange{x.begin(), x.end()};
};

template <elemwise_mditerable T>
auto elemwise_range(T&& x) {
  return std::ranges::subrange{x.elem_begin(), x.elem_end()};
};
}  // namespace matpack
