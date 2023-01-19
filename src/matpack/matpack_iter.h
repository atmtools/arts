#pragma once

#include "matpack_concepts.h"

namespace matpack {
//! Return an array of jokers
template <Index N> [[nodiscard]] constexpr std::array<Joker, N> jokers() {
  std::array<Joker, N> out;
  out.fill(joker);
  return out;
}

//! Create a tuple to access an index surrounded by jokers
template <Index M, Index N> [[nodiscard]] constexpr auto tup(Index i) {
  if constexpr (M > 0 and N - M - 1 > 0) {
    return std::tuple_cat(jokers<M>(), std::tuple{i}, jokers<N - M - 1>());
  } else if constexpr (M > 0) {
    return std::tuple_cat(jokers<M>(), std::tuple{i});
  } else if constexpr (N - M - 1 > 0) {
    return std::tuple_cat(std::tuple{i}, jokers<N - M - 1>());
  } else {
    return std::tuple{i};
  }
}

//! Helper to return the original span type when accessed by a joker
template <Index M, Index N, class mdspan_type>
[[nodiscard]] constexpr auto sub(mdspan_type &v, Joker) {
  return v;
}

//! Helper to return  smaller span type when accessed by a index
template <Index M, Index N, class mdspan_type>
[[nodiscard]] constexpr auto sub(mdspan_type &v, Index i) {
  using namespace stdx;
  return std::apply(
      [&v](auto &&...slices) {
        if constexpr (requires { v(slices...); }) {
          return v(slices...);
        } else {
          return submdspan(v, slices...);
        }
      },
      tup<M, N>(i));
}

//! The multidimensional iteration type --- FIXME: tested only for M == 0
template <Index M, bool constant, class mdspan_type,
          class submdspan_type>
class matpack_mditer {
  //! The rank of the original mdspan
  static constexpr auto N = mdspan_type::rank();

  //! The underlying type of the original mdspan
  using T = std::conditional_t<constant, const typename mdspan_type::value_type,
                               typename mdspan_type::value_type>;

  static_assert(N > 0, "Not for 0-rank mdspan");
  static_assert(M >= 0 and M < N, "Out-of-rank");

  //! The iterator represents this dimension along M
  Index pos{0};

  //! The original data is at this location
  std::conditional_t<constant, const mdspan_type *, mdspan_type *> orig{nullptr};

public:
  using difference_type = Index;
  using value_type = std::conditional_t<N == 1, T, submdspan_type>;

  constexpr matpack_mditer() = default;
  constexpr matpack_mditer(matpack_mditer&&) noexcept = default;
  constexpr matpack_mditer(const matpack_mditer&) = default;
  constexpr matpack_mditer& operator=(matpack_mditer&&) noexcept = default;
  constexpr matpack_mditer& operator=(const matpack_mditer&) = default;

  //! Construct the iterator from this type
  constexpr matpack_mditer(mdspan_type &x) requires(not constant) : orig(&x) {}

  //! Construct the iterator from this type
  constexpr matpack_mditer(const mdspan_type &x) requires(constant) : orig(&x) {}

  constexpr matpack_mditer& operator++() noexcept {pos++; return *this;}
  constexpr matpack_mditer& operator--() noexcept {pos--; return *this;}
  [[nodiscard]] constexpr matpack_mditer operator++(int) const noexcept {matpack_mditer out(*this); ++out.pos; return out;}
  [[nodiscard]] constexpr matpack_mditer operator--(int) const noexcept {matpack_mditer out(*this); --out.pos; return out;}

  constexpr matpack_mditer& operator+=(Index i) noexcept {pos+=i; return *this;}
  constexpr matpack_mditer& operator-=(Index i) noexcept {pos-=i; return *this;}
  [[nodiscard]] constexpr matpack_mditer operator+(Index i) const noexcept {matpack_mditer out(*this); out.pos+=i; return out;}
  [[nodiscard]] constexpr matpack_mditer operator-(Index i) const noexcept {matpack_mditer out(*this); out.pos-=i; return out;}
  [[nodiscard]] constexpr friend matpack_mditer operator+(Index i, const matpack_mditer& m) noexcept {return m + i;}
  [[nodiscard]] constexpr friend matpack_mditer operator-(Index i, const matpack_mditer& m) noexcept {return m - i;}

  [[nodiscard]] constexpr difference_type operator-(const matpack_mditer& other) const noexcept {return pos-other.pos;}
  [[nodiscard]] constexpr difference_type operator+(const matpack_mditer& other) const noexcept {return pos+other.pos;}

  [[nodiscard]] constexpr auto operator<=>(const matpack_mditer&) const noexcept = default;

  [[nodiscard]] constexpr auto operator*() const
      -> std::conditional_t<N == 1, std::conditional_t<constant, T, T&>,
                            value_type> {
    if constexpr (N == 1) {
      return orig->elem_at(pos);
    } else {
      return sub<M, N>(*orig, pos);
    }
  }

  [[nodiscard]] constexpr auto operator[](Index i) const
      -> std::conditional_t<N == 1, std::conditional_t<constant, T, T&>,
                            value_type> {
    if constexpr (N == 1) {
      return orig->elem_at(pos + i);
    } else {
      return sub<M, N>(*orig, pos + i);
    }
  }
};

//! An elementwise iterator for a mdspan type
template <typename T, bool constant, class mdspan_type>
class matpack_elemwise_mditer {
  //! The position of this object in the original data
  Index pos{0};

  //! A pointer to the original data
  std::conditional_t<constant, const mdspan_type *, mdspan_type *> orig{nullptr};

public:
  using difference_type = Index;
  using value_type = std::conditional_t<constant, const T, T>;

  constexpr matpack_elemwise_mditer() = default;
  constexpr matpack_elemwise_mditer(matpack_elemwise_mditer&&) noexcept = default;
  constexpr matpack_elemwise_mditer(const matpack_elemwise_mditer&) = default;
  constexpr matpack_elemwise_mditer& operator=(matpack_elemwise_mditer&&) noexcept = default;
  constexpr matpack_elemwise_mditer& operator=(const matpack_elemwise_mditer&) = default;

  //! Construct the iterator from this type
  constexpr matpack_elemwise_mditer(mdspan_type &x) requires(not constant) : orig(&x) {}

  //! Construct the iterator from this type
  constexpr matpack_elemwise_mditer(const mdspan_type &x) requires(constant) : orig(&x) {}

  constexpr matpack_elemwise_mditer& operator++() noexcept {pos++; return *this;}
  constexpr matpack_elemwise_mditer& operator--() noexcept {pos--; return *this;}
  [[nodiscard]] constexpr matpack_elemwise_mditer operator++(int) const noexcept {matpack_elemwise_mditer out(*this); ++out.pos; return out;}
  [[nodiscard]] constexpr matpack_elemwise_mditer operator--(int) const noexcept {matpack_elemwise_mditer out(*this); --out.pos; return out;}

  constexpr matpack_elemwise_mditer& operator+=(Index i) noexcept {pos+=i; return *this;}
  constexpr matpack_elemwise_mditer& operator-=(Index i) noexcept {pos-=i; return *this;}
  [[nodiscard]] constexpr matpack_elemwise_mditer operator+(Index i) const noexcept {matpack_elemwise_mditer out(*this); out.pos+=i; return out;}
  [[nodiscard]] constexpr matpack_elemwise_mditer operator-(Index i) const noexcept {matpack_elemwise_mditer out(*this); out.pos-=i; return out;}
  [[nodiscard]] constexpr friend matpack_elemwise_mditer operator+(Index i, const matpack_elemwise_mditer& m) noexcept {return m + i;}
  [[nodiscard]] constexpr friend matpack_elemwise_mditer operator-(Index i, const matpack_elemwise_mditer& m) noexcept {return m - i;}

  [[nodiscard]] constexpr difference_type operator-(const matpack_elemwise_mditer& other) const noexcept {return pos-other.pos;}
  [[nodiscard]] constexpr difference_type operator+(const matpack_elemwise_mditer& other) const noexcept {return pos+other.pos;}

  [[nodiscard]] constexpr auto operator<=>(const matpack_elemwise_mditer&) const noexcept = default;

  [[nodiscard]] constexpr auto operator*() const
      -> std::conditional_t<constant, T, T&> {
    return orig -> elem_at(pos);
  }

  [[nodiscard]] constexpr auto operator[](Index i) const
      -> std::conditional_t<constant, T, T&> {
    return orig -> elem_at(pos + i);
  }
};
}  // namespace matpack
