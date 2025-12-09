#pragma once

#include <configtypes.h>

#include <array>
#include <cstddef>
#include <tuple>
#include <utility>

#include "matpack_mdspan_flat_shape_pos.h"

namespace matpack {
template <class... Formats, Size N, Size... Is>
std::tuple<Formats...> as_tuple(const std::array<std::nullptr_t, N>& arr,
                                std::index_sequence<Is...>) {
  return std::make_tuple(Formats{arr[Is]}...);
}

template <class... Formats, size_t N>
std::tuple<Formats...> as_tuple(const std::array<std::nullptr_t, N>& arr)
  requires(N == sizeof...(Formats))
{
  return as_tuple<Formats...>(arr, std::make_index_sequence<N>{});
}

template <typename... iters>
  requires(((rank<iters>() == 1) and ...))
struct elemwise {
  struct elemwise_iteration {
    static constexpr Index N = sizeof...(iters);

    matpack::flat_shape_pos<N> pos{std::array<Index, N>{}};
    const std::tuple<const iters* const...> orig{
        as_tuple<const iters* const...>(std::array<std::nullptr_t, N>{})};

    constexpr elemwise_iteration()                              = default;
    constexpr elemwise_iteration(elemwise_iteration&&) noexcept = default;
    constexpr elemwise_iteration(const elemwise_iteration&)     = default;
    constexpr elemwise_iteration& operator=(
        elemwise_iteration&& other) noexcept {
      this->pos  = std::move(other.pos);
      this->orig = std::move(other.orig);
      return *this;
    }
    constexpr elemwise_iteration& operator=(const elemwise_iteration& other) {
      this->pos  = other.pos;
      this->orig = other.orig;
      return *this;
    }

    constexpr elemwise_iteration(const iters&... x)
        : pos(std::array{static_cast<Index>(x.size())...}), orig((&x)...) {}

    constexpr elemwise_iteration& operator++() noexcept {
      pos++;
      return *this;
    }
    constexpr elemwise_iteration& operator--() noexcept {
      pos--;
      return *this;
    }
    constexpr elemwise_iteration operator++(int) noexcept {
      elemwise_iteration out(*this);
      ++pos;
      return out;
    }
    constexpr elemwise_iteration operator--(int) noexcept {
      elemwise_iteration out(*this);
      --pos;
      return out;
    }

    constexpr elemwise_iteration& operator+=(Index i) noexcept {
      pos += i;
      return *this;
    }
    constexpr elemwise_iteration& operator-=(Index i) noexcept {
      pos -= i;
      return *this;
    }
    [[nodiscard]] constexpr elemwise_iteration operator+(
        Index i) const noexcept {
      elemwise_iteration out(*this);
      out.pos += i;
      return out;
    }
    [[nodiscard]] constexpr elemwise_iteration operator-(
        Index i) const noexcept {
      elemwise_iteration out(*this);
      out.pos -= i;
      return out;
    }
    [[nodiscard]] constexpr friend elemwise_iteration operator+(
        Index i, const elemwise_iteration& m) noexcept {
      return m + i;
    }
    [[nodiscard]] constexpr friend elemwise_iteration operator-(
        Index i, const elemwise_iteration& m) noexcept {
      return m - i;
    }

    [[nodiscard]] constexpr Index operator-(
        const elemwise_iteration& other) const noexcept {
      return pos - other.pos;
    }
    [[nodiscard]] constexpr Index operator+(
        const elemwise_iteration& other) const noexcept {
      return pos + other.pos;
    }

    [[nodiscard]] constexpr auto operator==(
        const elemwise_iteration& other) const noexcept {
      return pos == other.pos;
    }
    [[nodiscard]] constexpr auto operator!=(
        const elemwise_iteration& other) const noexcept {
      return pos != other.pos;
    }
    [[nodiscard]] constexpr auto operator<(
        const elemwise_iteration& other) const noexcept {
      return pos < other.pos;
    }
    [[nodiscard]] constexpr auto operator>(
        const elemwise_iteration& other) const noexcept {
      return pos > other.pos;
    }
    [[nodiscard]] constexpr auto operator<=(
        const elemwise_iteration& other) const noexcept {
      return pos <= other.pos;
    }
    [[nodiscard]] constexpr auto operator>=(
        const elemwise_iteration& other) const noexcept {
      return pos >= other.pos;
    }

    template <Index... ints>
    constexpr auto values(std::integer_sequence<Index, ints...>) const {
      return std::apply(
          [this](auto&&... i) {
            return std::tuple{(*std::get<ints>(orig))[i]...};
          },
          pos.pos);
    }

    [[nodiscard]] constexpr auto operator*() const {
      return values(std::make_integer_sequence<Index, N>{});
    }

    [[nodiscard]] constexpr auto operator[](Index i) const {
      return *(*this + i);
    }
  };

  elemwise_iteration d;
  Index len;
  constexpr elemwise(iters&... x)
      : d(x...), len(static_cast<Index>((... * x.size()))) {}
  constexpr elemwise(const iters&... x)
      : d(x...), len(static_cast<Index>((... * x.size()))) {}
  constexpr auto begin() { return d; }
  constexpr auto end() { return d + len; }
  constexpr auto begin() const { return d; }
  constexpr auto end() const { return d + len; }
};
}  // namespace matpack