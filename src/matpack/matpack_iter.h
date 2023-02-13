#pragma once

#include "matpack_concepts.h"
#include <__concepts/same_as.h>
#include <tuple>

namespace matpack {
//! Return an array of jokers
template <Index N> [[nodiscard]] constexpr std::array<Joker, N> jokers() {
  std::array<Joker, N> out;
  out.fill(joker);
  return out;
}

//! Create a tuple to access an index surrounded by jokers
template <Index M, Index N, typename T> [[nodiscard]] constexpr auto tup(T&& i) {
  if constexpr (M > 0 and N - M - 1 > 0) {
    return std::tuple_cat(jokers<M>(), std::tuple{std::forward<T>(i)}, jokers<N - M - 1>());
  } else if constexpr (M > 0) {
    return std::tuple_cat(jokers<M>(), std::tuple{std::forward<T>(i)});
  } else if constexpr (N - M - 1 > 0) {
    return std::tuple_cat(std::tuple{std::forward<T>(i)}, jokers<N - M - 1>());
  } else {
    return std::tuple{std::forward<T>(i)};
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

//! Helper to return  smaller span type when accessed by a index
template <Index M, Index N, class mdspan_type>
[[nodiscard]] constexpr auto sub(mdspan_type &v, const stdx::strided_slice<Index, Index, Index>& i) {
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
  constexpr matpack_mditer operator++(int) noexcept {matpack_mditer out(*this); ++pos; return out;}
  constexpr matpack_mditer operator--(int) noexcept {matpack_mditer out(*this); --pos; return out;}

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
      return orig->operator[](pos);
    } else {
      return sub<M, N>(*orig, pos);
    }
  }

  [[nodiscard]] constexpr auto operator[](Index i) const
      -> std::conditional_t<N == 1, std::conditional_t<constant, T, T&>,
                            value_type> {
    if constexpr (N == 1) {
      return orig->operator[](pos + i);
    } else {
      return sub<M, N>(*orig, pos + i);
    }
  }
};

//! Helper struct to pretty-print the shape of the object
template <std::size_t N> struct shape_help {
std::array<Index, N> shape;
constexpr shape_help(const std::array<Index, N>& x) : shape(x) {}
friend std::ostream& operator<<(std::ostream& os, const shape_help& sh) {
  bool first=true;
  os << '(';
  for (auto& x: sh.shape) {
    if (not first) os << ", ";
    first = false;
    os << x;
  }
  return os << ')';
}
};

//! Return a constant Index array of size N with value v
template <Index N, Index v> 
constexpr std::array<Index, N> constant_array() {
  std::array<Index, N> x;
  x.fill(v);
  return x;
}

template <Index N>
struct flat_shape_pos {
  std::array<Index, N> pos{constant_array<N, 0>()};
  std::array<Index, N> shp{constant_array<N, 0>()};

  constexpr flat_shape_pos(std::array<Index, N> shape) : shp(shape) {}
  constexpr flat_shape_pos(std::array<Index, N> shape, Index n) : shp(shape) {pos.back() += n; adapt();}

  constexpr flat_shape_pos& operator++() noexcept {++pos.back(); adapt(); return *this;}
  constexpr flat_shape_pos& operator--() noexcept {--pos.back(); adapt(); return *this;}
  constexpr flat_shape_pos operator++(int) noexcept {flat_shape_pos out(*this); ++pos.back(); adapt(); return out;}
  constexpr flat_shape_pos operator--(int) noexcept {flat_shape_pos out(*this); --pos.back(); adapt(); return out;}

  constexpr flat_shape_pos& operator+=(Index i) noexcept {pos.back()+=i; adapt(); return *this;}
  constexpr flat_shape_pos& operator-=(Index i) noexcept {pos.back()-=i; adapt(); return *this;}
  [[nodiscard]] constexpr flat_shape_pos operator+(Index i) const noexcept {flat_shape_pos out(*this); out.pos.back()+=i; out.adapt(); return out;}
  [[nodiscard]] constexpr flat_shape_pos operator-(Index i) const noexcept {flat_shape_pos out(*this); out.pos.back()-=i; out.adapt(); return out;}
  [[nodiscard]] constexpr friend flat_shape_pos operator+(Index i, const flat_shape_pos& m) noexcept {return m + i;}
  [[nodiscard]] constexpr friend flat_shape_pos operator-(Index i, const flat_shape_pos& m) noexcept {return m - i;}

  [[nodiscard]] constexpr Index to_index() const {
    Index str=1;
    Index out=pos.back();
    for (Index i=N-2; i>=0; i--) {
      str *= shp[i+1];
      out += str * pos[i];
    }
    return out;
  }

  [[nodiscard]] constexpr Index operator-(const flat_shape_pos& other) const noexcept {return to_index()-other.to_index();}
  [[nodiscard]] constexpr Index operator+(const flat_shape_pos& other) const noexcept {return to_index()+other.to_index();}

  [[nodiscard]] constexpr auto operator==(const flat_shape_pos& other) const noexcept {return pos == other.pos;}
  [[nodiscard]] constexpr auto operator!=(const flat_shape_pos& other) const noexcept {return pos != other.pos;}
  [[nodiscard]] constexpr auto operator<(const flat_shape_pos& other) const noexcept {return pos < other.pos;}
  [[nodiscard]] constexpr auto operator>(const flat_shape_pos& other) const noexcept {return pos > other.pos;}
  [[nodiscard]] constexpr auto operator<=(const flat_shape_pos& other) const noexcept {return pos <= other.pos;}
  [[nodiscard]] constexpr auto operator>=(const flat_shape_pos& other) const noexcept {return pos >= other.pos;}

  friend std::ostream& operator<<(std::ostream& os, const flat_shape_pos& fsp) {
    return os << "POS: " << shape_help<N>(fsp.pos) << " SHP: " << shape_help<N>(fsp.shp);
  }

private:
  constexpr void adapt() noexcept {
    for (Index i=N-1; i>0; i--) {
      while (shp[i] > 0 and pos[i] >= shp[i]) {
        pos[i] -= shp[i];
        pos[i-1]++;
      };
    }
  }
};

//! An elementwise iterator for a mdspan type
template <typename T, bool constant, class mdspan_type>
class matpack_elemwise_mditer {
  static constexpr Index N = rank<mdspan_type>();

  //! The position of this object in the original data
  flat_shape_pos<N> pos{constant_array<N, 0>()};

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

  constexpr matpack_elemwise_mditer(mdspan_type &x) requires(not constant) : pos(x.shape()), orig(&x) {}
  constexpr matpack_elemwise_mditer(const mdspan_type &x) requires(constant) : pos(x.shape()), orig(&x) {}

  constexpr matpack_elemwise_mditer& operator++() noexcept {pos++; return *this;}
  constexpr matpack_elemwise_mditer& operator--() noexcept {pos--; return *this;}
  constexpr matpack_elemwise_mditer operator++(int) noexcept {matpack_elemwise_mditer out(*this); ++pos; return out;}
  constexpr matpack_elemwise_mditer operator--(int) noexcept {matpack_elemwise_mditer out(*this); --pos; return out;}

  constexpr matpack_elemwise_mditer& operator+=(Index i) noexcept {pos+=i; return *this;}
  constexpr matpack_elemwise_mditer& operator-=(Index i) noexcept {pos-=i; return *this;}
  [[nodiscard]] constexpr matpack_elemwise_mditer operator+(Index i) const noexcept {matpack_elemwise_mditer out(*this); out.pos+=i; return out;}
  [[nodiscard]] constexpr matpack_elemwise_mditer operator-(Index i) const noexcept {matpack_elemwise_mditer out(*this); out.pos-=i; return out;}
  [[nodiscard]] constexpr friend matpack_elemwise_mditer operator+(Index i, const matpack_elemwise_mditer& m) noexcept {return m + i;}
  [[nodiscard]] constexpr friend matpack_elemwise_mditer operator-(Index i, const matpack_elemwise_mditer& m) noexcept {return m - i;}

  [[nodiscard]] constexpr difference_type operator-(const matpack_elemwise_mditer& other) const noexcept {return pos-other.pos;}
  [[nodiscard]] constexpr difference_type operator+(const matpack_elemwise_mditer& other) const noexcept {return pos+other.pos;}

  [[nodiscard]] constexpr auto operator==(const matpack_elemwise_mditer& other) const noexcept {return pos == other.pos;}
  [[nodiscard]] constexpr auto operator!=(const matpack_elemwise_mditer& other) const noexcept {return pos != other.pos;}
  [[nodiscard]] constexpr auto operator<(const matpack_elemwise_mditer& other) const noexcept {return pos < other.pos;}
  [[nodiscard]] constexpr auto operator>(const matpack_elemwise_mditer& other) const noexcept {return pos > other.pos;}
  [[nodiscard]] constexpr auto operator<=(const matpack_elemwise_mditer& other) const noexcept {return pos <= other.pos;}
  [[nodiscard]] constexpr auto operator>=(const matpack_elemwise_mditer& other) const noexcept {return pos >= other.pos;}

  [[nodiscard]] constexpr auto operator*() const
      -> std::conditional_t<constant, T, T &> {
    if constexpr (N == 1)
      return orig->operator[](pos.pos[0]);
    else
      return std::apply(
          [this](auto... inds) -> std::conditional_t<constant, T, T &> {
            return orig->operator()(inds...);
          },
          pos.pos);
  }

  [[nodiscard]] constexpr auto operator[](Index i) const
      -> std::conditional_t<constant, T, T &> {
    if constexpr (N == 1)
      return orig->operator[](pos.pos[0] + i);
    else
      return std::apply(
          [this](auto... inds) -> std::conditional_t<constant, T, T &> {
            return orig->operator()(inds...);
          },
          (pos + i).pos);
  }
};

template <typename ... iters> requires(((rank<iters>() == 1) and ...))
struct elemwise {
  class elemwise_iteration {
    static constexpr Index N = sizeof...(iters); 

    matpack::flat_shape_pos<N> pos{matpack::constant_array<N, 0>()};
    const std::tuple<iters*...> orig{std::array<nullptr_t, N>{}};

  public:
    constexpr elemwise_iteration() = default;
    constexpr elemwise_iteration(elemwise_iteration&&) noexcept = default;
    constexpr elemwise_iteration(const elemwise_iteration&) = default;
    constexpr elemwise_iteration& operator=(elemwise_iteration&&) noexcept = default;
    constexpr elemwise_iteration& operator=(const elemwise_iteration&) = default;

    constexpr elemwise_iteration(iters&...x) : pos(std::array{static_cast<Index>(x.size())...}), orig((&x)...) {}
    constexpr elemwise_iteration(const iters &...x) : pos(std::array{static_cast<Index>(x.size())...}), orig((&x)...) {}

    constexpr elemwise_iteration& operator++() noexcept {pos++; return *this;}
    constexpr elemwise_iteration& operator--() noexcept {pos--; return *this;}
    constexpr elemwise_iteration operator++(int) noexcept {elemwise_iteration out(*this); ++pos; return out;}
    constexpr elemwise_iteration operator--(int) noexcept {elemwise_iteration out(*this); --pos; return out;}

    constexpr elemwise_iteration& operator+=(Index i) noexcept {pos+=i; return *this;}
    constexpr elemwise_iteration& operator-=(Index i) noexcept {pos-=i; return *this;}
    [[nodiscard]] constexpr elemwise_iteration operator+(Index i) const noexcept {elemwise_iteration out(*this); out.pos+=i; return out;}
    [[nodiscard]] constexpr elemwise_iteration operator-(Index i) const noexcept {elemwise_iteration out(*this); out.pos-=i; return out;}
    [[nodiscard]] constexpr friend elemwise_iteration operator+(Index i, const elemwise_iteration& m) noexcept {return m + i;}
    [[nodiscard]] constexpr friend elemwise_iteration operator-(Index i, const elemwise_iteration& m) noexcept {return m - i;}

    [[nodiscard]] constexpr Index operator-(const elemwise_iteration& other) const noexcept {return pos-other.pos;}
    [[nodiscard]] constexpr Index operator+(const elemwise_iteration& other) const noexcept {return pos+other.pos;}

    [[nodiscard]] constexpr auto operator==(const elemwise_iteration& other) const noexcept {return pos == other.pos;}
    [[nodiscard]] constexpr auto operator!=(const elemwise_iteration& other) const noexcept {return pos != other.pos;}
    [[nodiscard]] constexpr auto operator<(const elemwise_iteration& other) const noexcept {return pos < other.pos;}
    [[nodiscard]] constexpr auto operator>(const elemwise_iteration& other) const noexcept {return pos > other.pos;}
    [[nodiscard]] constexpr auto operator<=(const elemwise_iteration& other) const noexcept {return pos <= other.pos;}
    [[nodiscard]] constexpr auto operator>=(const elemwise_iteration& other) const noexcept {return pos >= other.pos;}

  private:
  template<Index ... ints>
    constexpr auto values(std::integer_sequence<Index, ints...>) const {
      return std::apply([this](auto&&... i){return std::tuple{(*std::get<ints>(orig))[i] ...}; }, pos.pos);
    }

  public:

    [[nodiscard]] constexpr auto operator*() const {
      return values(std::make_integer_sequence<Index, N>{});
    }

    [[nodiscard]] constexpr auto operator[](Index i) const {
      return *(*this + i);
    }
  };

  elemwise_iteration d;
  Index len;
  constexpr elemwise(iters&... x) : d(x...), len(static_cast<Index>((x.size() * ...))) {}
  constexpr elemwise(const iters&... x) : d(x...), len(static_cast<Index>((x.size() * ...))) {}
  constexpr auto begin() {return d;}
  constexpr auto end() {return d+len;}
  constexpr auto begin() const {return d;}
  constexpr auto end() const {return d+len;}
};
}  // namespace matpack
