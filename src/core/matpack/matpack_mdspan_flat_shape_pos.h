#pragma once

#include <configtypes.h>

#include <array>

namespace matpack {
template <Index N>
struct flat_shape_pos {
  std::array<Index, N> pos{};
  std::array<Index, N> shp{};

  constexpr flat_shape_pos(std::array<Index, N> shape) : shp(shape) {}
  constexpr flat_shape_pos(std::array<Index, N> shape, Index n) : shp(shape) {
    pos.back() += n;
    adapt();
  }

  constexpr flat_shape_pos& operator++() noexcept {
    ++pos.back();
    adapt();
    return *this;
  }
  constexpr flat_shape_pos& operator--() noexcept {
    --pos.back();
    adapt();
    return *this;
  }
  constexpr flat_shape_pos operator++(int) noexcept {
    flat_shape_pos out(*this);
    ++pos.back();
    adapt();
    return out;
  }
  constexpr flat_shape_pos operator--(int) noexcept {
    flat_shape_pos out(*this);
    --pos.back();
    adapt();
    return out;
  }

  constexpr flat_shape_pos& operator+=(Index i) noexcept {
    pos.back() += i;
    adapt();
    return *this;
  }
  constexpr flat_shape_pos& operator-=(Index i) noexcept {
    pos.back() -= i;
    adapt();
    return *this;
  }
  [[nodiscard]] constexpr flat_shape_pos operator+(Index i) const noexcept {
    flat_shape_pos out(*this);
    out.pos.back() += i;
    out.adapt();
    return out;
  }
  [[nodiscard]] constexpr flat_shape_pos operator-(Index i) const noexcept {
    flat_shape_pos out(*this);
    out.pos.back() -= i;
    out.adapt();
    return out;
  }
  [[nodiscard]] constexpr friend flat_shape_pos operator+(
      Index i, const flat_shape_pos& m) noexcept {
    return m + i;
  }
  [[nodiscard]] constexpr friend flat_shape_pos operator-(
      Index i, const flat_shape_pos& m) noexcept {
    return m - i;
  }

  [[nodiscard]] constexpr Index to_index() const {
    Index str = 1;
    Index out = pos.back();
    for (Index i = N - 2; i >= 0; i--) {
      str *= shp[i + 1];
      out += str * pos[i];
    }
    return out;
  }

  [[nodiscard]] constexpr Index operator-(
      const flat_shape_pos& other) const noexcept {
    return to_index() - other.to_index();
  }
  [[nodiscard]] constexpr Index operator+(
      const flat_shape_pos& other) const noexcept {
    return to_index() + other.to_index();
  }

  [[nodiscard]] constexpr auto operator==(
      const flat_shape_pos& other) const noexcept {
    return pos == other.pos;
  }
  [[nodiscard]] constexpr auto operator!=(
      const flat_shape_pos& other) const noexcept {
    return pos != other.pos;
  }
  [[nodiscard]] constexpr auto operator<(
      const flat_shape_pos& other) const noexcept {
    return pos < other.pos;
  }
  [[nodiscard]] constexpr auto operator>(
      const flat_shape_pos& other) const noexcept {
    return pos > other.pos;
  }
  [[nodiscard]] constexpr auto operator<=(
      const flat_shape_pos& other) const noexcept {
    return pos <= other.pos;
  }
  [[nodiscard]] constexpr auto operator>=(
      const flat_shape_pos& other) const noexcept {
    return pos >= other.pos;
  }

 private:
  constexpr void adapt() noexcept {
    for (Index i = N - 1; i > 0; i--) {
      while (shp[i] > 0 and pos[i] >= shp[i]) {
        pos[i] -= shp[i];
        pos[i - 1]++;
      };
    }
  }
};
}  // namespace matpack
