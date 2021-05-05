#ifndef grids_h
#define grids_h

#include <array>
#include <iostream>
#include <vector>

/*! Compute the multiplication of all inds */
template <typename... Inds>
constexpr std::size_t mul(Inds... inds) noexcept {
  return (std::size_t(inds) * ...);
}

/*! Compute the multiplication of all inds in arr */
template <std::size_t N>
constexpr std::size_t mul(const std::array<std::size_t, N>& arr) noexcept {
  if constexpr (N == 0) {
    return 0;
  } else {
    std::size_t out = 1;
    for (auto i : arr) out *= i;
    return out;
  }
}

/*! Compute gridsizes from indices */
template <typename... Inds>
constexpr std::array<std::size_t, sizeof...(Inds)> gridsize_from_index(Inds... inds) noexcept {
  std::array<std::size_t, sizeof...(Inds)> out{};
  std::array<std::size_t, sizeof...(Inds)> arr{std::size_t(inds)...};
  std::size_t i = sizeof...(Inds)-1;
  std::size_t s = 1;
  while (i < sizeof...(Inds)) {
    out[i] = s;
    s *= arr[i];
    i--;
  }
  return out;
}

/*! Compute gridsizes from indices */
template <std::size_t N>
constexpr std::size_t index_from_gridsize(std::array<std::size_t, N> gridsize,
                                          std::array<std::size_t, N> inds) noexcept {
  std::size_t pos=0;
  for (std::size_t i=0; i<N; i++) pos += gridsize[i] * inds[i];
  return pos;
}

/** Row-major grid creation */
template <typename b, std::size_t n>
class Grid {
  /*! List of values */
  std::vector<b> ptr;
  
  /*! List of sizes of the grid so that the numbers represent the stepping */
  std::array<std::size_t, n> gridsize;
  
  /*! Number of elements in total */
  std::size_t size() const noexcept { return ptr.size(); }
  
public:
  /*! Number of dimensions */
  static constexpr std::size_t N = n;
  
  /*! Base unit of the Grid */
  using base = b;
  
  static_assert(N, "Must have size");
  
  /*! Initialization of grid with size of the input list of indices */
  template <typename... Inds>
  Grid(Inds... inds) noexcept : ptr(mul(inds...)), gridsize(gridsize_from_index(inds...)) {
    static_assert(sizeof...(Inds) == N,
                  "Must have same size for initialization");
  }
  
  /*! Move initialization */
  Grid(Grid&& g) noexcept : ptr(std::move(g.ptr)), gridsize(std::move(g.gridsize)) {}
  
  /*! Move and set operator */
  Grid& operator=(Grid&& g) noexcept {
    ptr = std::move(g.ptr);
    gridsize = std::move(g.gridsize);
    return *this;
  }
  
  /*! Access operator */
  template <typename... Inds>
  base& operator()(Inds... inds) noexcept {
    return ptr[index_from_gridsize(gridsize, std::array<std::size_t, N>{std::size_t(inds)...})];
  }
  
  /*! Access operator */
  template <typename... Inds>
  const base& operator()(Inds... inds) const noexcept {
    return ptr[index_from_gridsize(gridsize, std::array<std::size_t, N>{std::size_t(inds)...})];
  }
  
  /*! Access operator for VectorType */
  base& operator[](std::size_t ind) noexcept {
    static_assert(N == 1);
    return ptr[ind];
  }
  
  /*! Access operator for VectorType */
  const base& operator[](std::size_t ind) const noexcept {
    static_assert(N == 1);
    return ptr[ind];
  }
  
  /*! Friendly stream operator */
  friend std::ostream& operator<<(std::ostream& os, const Grid& g) {
    const std::size_t nel = g.size();
    for (std::size_t i = 0; i < nel; i++) {
      os << g.ptr[i];
      if constexpr (N > 1) {
        if (i not_eq 0 and i not_eq nel and i % g.gridsize[N-2] == 0)
          os << '\n';
        else if (i not_eq nel)
          os << ' ';
      } else if (i not_eq nel) {
        os << ' ';
      }
    }
    return os;
  }
};  // Grid

/** Row-major fixed grid creation */
template <typename b, std::size_t... Sizes>
class FixedGrid {
  /*! List of values */
  std::array<b, mul(Sizes...)> ptr;
  
public:
  /*! Number of dimensions */
  static constexpr std::size_t N = sizeof...(Sizes);
  
  /*! Base unit of the Grid */
  using base = b;
  
  static_assert(N, "Must have size");
  
  /*! throw-away constructor so that Grid and FixedGrid can be output in same templates */
  template <typename... Inds>
  constexpr FixedGrid(Inds... inds [[maybe_unused]]) noexcept : ptr({}) {
    static_assert(sizeof...(Inds) == N, "Must have same size for initialization");
  }
  
  /*! Standard constructor uses standard constructor of the base */
  constexpr FixedGrid() noexcept : ptr({}) {static_assert(mul(Sizes...), "Must have size");}
  
  /*! Move initialization */
  constexpr FixedGrid(FixedGrid&& g) noexcept : ptr(std::move(g.ptr)) {}
  
  /*! Move and set operator */
  constexpr FixedGrid& operator=(FixedGrid&& g) noexcept {
    ptr = std::move(g.ptr);
    return *this;
  }
  
  /*! Access operator */
  template <typename... Inds>
  constexpr base& operator()(Inds... inds) noexcept {
    return ptr[index_from_gridsize(gridsize_from_index(Sizes...), std::array<std::size_t, N>{std::size_t(inds)...})];
  }
  
  /*! Access operator */
  template <typename... Inds>
  constexpr const base& operator()(Inds... inds) const noexcept {
    return ptr[index_from_gridsize(gridsize_from_index(Sizes...), std::array<std::size_t, N>{std::size_t(inds)...})];
  }
  
  /*! Access operator for VectorType */
  constexpr base& operator[](std::size_t ind) noexcept {
    static_assert(N == 1);
    return ptr[ind];
  }
  
  /*! Access operator for VectorType */
  constexpr const base& operator[](std::size_t ind) const noexcept {
    static_assert(N == 1);
    return ptr[ind];
  }
  
  /*! Friendly stream operator */
  friend std::ostream& operator<<(std::ostream& os, const FixedGrid& g) {
    constexpr std::size_t nel = mul(Sizes...);
    constexpr std::size_t last_of = (std::array<std::size_t, N>{std::size_t(Sizes)...}).back();
    for (std::size_t i = 0; i < nel; i++) {
      os << g.ptr[i];
      if (i not_eq 0 and i not_eq nel and i % last_of == 0)
        os << '\n';
      else if (i not_eq nel)
        os << ' ';
    }
    return os;
  }
};  // FixedGrid

#endif  // grids_h
