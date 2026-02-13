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

  using iterator_concept  = std::random_access_iterator_tag;
  using iterator_category = std::random_access_iterator_tag;
  using difference_type   = Index;
  using value_type        = std::remove_cvref_t<decltype(data->elem_at(0))>;

  constexpr elemwise_mditer()                                      = default;
  constexpr elemwise_mditer(elemwise_mditer&&) noexcept            = default;
  constexpr elemwise_mditer(const elemwise_mditer&)                = default;
  constexpr elemwise_mditer& operator=(elemwise_mditer&&) noexcept = default;
  constexpr elemwise_mditer& operator=(const elemwise_mditer&)     = default;

  constexpr elemwise_mditer(T& x) : data(&x) {}
  constexpr elemwise_mditer(T* x, Index p) : data(x), pos(p) {}

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
    return {data, pos + i};
  }

  [[nodiscard]] constexpr elemwise_mditer operator-(Index i) const noexcept {
    return {data, pos - i};
  }

  [[nodiscard]] constexpr friend elemwise_mditer operator+(
      Index i, const elemwise_mditer& m) noexcept {
    return {m.data, m.pos + i};
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
      const elemwise_mditer& m) const noexcept {
    return pos <=> m.pos;
  }

  [[nodiscard]] constexpr bool operator==(
      const elemwise_mditer& m) const noexcept {
    return pos == m.pos;
  }

  [[nodiscard]] constexpr decltype(auto) operator*() const {
    return data->elem_at(pos);
  }

  [[nodiscard]] constexpr decltype(auto) operator[](Index i) const {
    return data->elem_at(pos + i);
  }

  friend constexpr void iter_swap(const elemwise_mditer& a,
                                  const elemwise_mditer& b) noexcept {
    using std::swap;
    swap(*a, *b);
  }
};

template <any_md T>
constexpr auto elemwise_range(T&& x) {
  if constexpr (requires {
                  std::span{x.view_as(
                      std::array<Index, 1>{static_cast<Index>(x.size())})};
                })
    return std::span{
        x.view_as(std::array<Index, 1>{static_cast<Index>(x.size())})};
  else
    return stdr::subrange{x.elem_begin(), x.elem_end()};
};

// Catch for non-matpack types (might not work, be careful)
template <ranked<1> T>
constexpr auto elemwise_range(T&& x)
  requires(not any_md<T>)
{
  if constexpr (requires { std::span{x}; })
    return std::span{x};
  else
    return stdr::subrange{x.begin(), x.end()};
};

struct elemwise_r {
  template <matpack::any_md T>
  friend constexpr auto operator|(T&& x, elemwise_r) {
    return elemwise_range(std::forward<T>(x));
  }
};
}  // namespace matpack

static constexpr matpack::elemwise_r by_elem{};
