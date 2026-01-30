#pragma once

#include <configtypes.h>

#include <algorithm>
#include <iterator>
#include <type_traits>

namespace matpack {
template <class T>
struct left_mditer {
  T* data{nullptr};
  Index pos{0};

  using iterator_concept  = std::random_access_iterator_tag;
  using iterator_category = std::random_access_iterator_tag;
  using difference_type   = Index;
  using value_type        = std::remove_cvref_t<decltype(data->operator[](0))>;

  constexpr left_mditer()                                  = default;
  constexpr left_mditer(left_mditer&&) noexcept            = default;
  constexpr left_mditer(const left_mditer&)                = default;
  constexpr left_mditer& operator=(left_mditer&&) noexcept = default;
  constexpr left_mditer& operator=(const left_mditer&)     = default;

  constexpr left_mditer(T& x) : data(&x) {}

  constexpr left_mditer& operator++() noexcept {
    pos++;
    return *this;
  }
  constexpr left_mditer& operator--() noexcept {
    pos--;
    return *this;
  }
  constexpr left_mditer operator++(int) noexcept {
    left_mditer out(*this);
    ++pos;
    return out;
  }
  constexpr left_mditer operator--(int) noexcept {
    left_mditer out(*this);
    --pos;
    return out;
  }

  constexpr left_mditer& operator+=(Index i) noexcept {
    pos += i;
    return *this;
  }
  constexpr left_mditer& operator-=(Index i) noexcept {
    pos -= i;
    return *this;
  }
  [[nodiscard]] constexpr left_mditer operator+(Index i) const noexcept {
    left_mditer out(*this);
    out.pos += i;
    return out;
  }
  [[nodiscard]] constexpr left_mditer operator-(Index i) const noexcept {
    left_mditer out(*this);
    out.pos -= i;
    return out;
  }
  [[nodiscard]] constexpr friend left_mditer operator+(
      Index i, const left_mditer& m) noexcept {
    return m + i;
  }
  [[nodiscard]] constexpr friend left_mditer operator-(
      Index i, const left_mditer& m) noexcept {
    return m - i;
  }

  [[nodiscard]] constexpr difference_type operator-(
      const left_mditer& other) const noexcept {
    return pos - other.pos;
  }
  [[nodiscard]] constexpr difference_type operator+(
      const left_mditer& other) const noexcept {
    return pos + other.pos;
  }

  [[nodiscard]] constexpr auto operator<=>(const left_mditer&) const noexcept =
      default;

  [[nodiscard]] constexpr bool operator==(const left_mditer&) const noexcept =
      default;

  [[nodiscard]] constexpr decltype(auto) operator*() const {
    return data->operator[](pos);
  }

  [[nodiscard]] constexpr decltype(auto) operator[](Index i) const {
    return data->operator[](pos + i);
  }

  friend constexpr auto iter_move(const left_mditer& x) noexcept { return *x; }

  friend constexpr void iter_swap(const left_mditer& a,
                                  const left_mditer& b) noexcept {
    using std::swap;
    swap(*a, *b);
  }
};
}  // namespace matpack