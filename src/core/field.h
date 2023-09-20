/**
 * @file   field.h
 * @author Richard Larsson
 * @date   2019-02-26
 * 
 * @brief This file contains the definition of Field3D.
 */

#ifndef FIELD_HEADER
#define FIELD_HEADER
#include "interpolation.h"
#include <vector>

/** Creates a 3D field of a base unit */
template <class base>
class Field3D {
 private:
  size_t mpages, mrows, mcols;
  std::vector<base> data;

 public:
  /**
   * @brief Construct a new Field3D object
   * 
   * @param[in] g Another field
   */
  Field3D(const Field3D& g) = default;

  /**
   * @brief Default assignment operator
   * 
   * @param[in] g Another field
   * @return Field3D<base>& *this
   */
  Field3D& operator=(const Field3D& g) = default;

  /**
   * @brief Default move operator
   * 
   * @param[in] g Another field
   * @return Field3D<base>& *this
   */
  Field3D& operator=(Field3D&& g) noexcept = default;

  /**
   * @brief Construct a new Field3D object
   * 
   * @param[in] g Another field
   */
  Field3D(Field3D&& g) noexcept
      : mpages(std::move(g.mpages)),
        mrows(std::move(g.mrows)),
        mcols(std::move(g.mcols)),
        data(std::move(g.data)) {}

  /**
   * @brief Construct a new Field 3 D object
   * 
   * @param[in] pages Number of pages
   * @param[in] rows Number of rows
   * @param[in] cols Numeber of columns
   * @param[in] init Const value
   */
  Field3D(size_t pages = 0,
          size_t rows = 0,
          size_t cols = 0,
          const base& init = base())
      : mpages(pages),
        mrows(rows),
        mcols(cols),
        data(cols * rows * pages, init) {}

  /**
   * @brief Access operator
   * 
   * Returns a ref to the object that can be changed in place
   * 
   * @param[in] page Outer dim
   * @param[in] row Middle dim
   * @param[in] col Inner dim
   * @return base& 
   */
  base& operator()(size_t page = 0, size_t row = 0, size_t col = 0) {
    return data[col + row * mcols + page * mrows * mcols];
  }

  /**
   * @brief Access operator
   * 
   * Returns a ref to the object that cannot be changed
   * 
   * @param[in] page Outer dim
   * @param[in] row Middle dim
   * @param[in] col Inner dim
   * @return base& 
   */
  const base& operator()(size_t col = 0,
                         size_t row = 0,
                         size_t page = 0) const {
    return data[col + row * mcols + page * mrows * mcols];
  }


  /**
   * @brief Weighted access operator by GridPos
   * 
   * Returns a weighted new object.  The object
   * must support multiplication with a Numeric
   * and the addition of (object) * (Numeric).
   * 
   * @param[in] page Outer dim
   * @param[in] row Middle dim
   * @param[in] col Inner dim
   * @return base& 
   */
  base operator()(const GridPos& page = {0, {0, 1}},
                  const GridPos& row = {0, {0, 1}},
                  const GridPos& col = {0, {0, 1}}) const {
    const std::array pos{
        col.idx + 0 + (row.idx + 0) * mcols + (page.idx + 0) * mrows * mcols,
        col.idx + 1 + (row.idx + 0) * mcols + (page.idx + 0) * mrows * mcols,
        col.idx + 0 + (row.idx + 1) * mcols + (page.idx + 0) * mrows * mcols,
        col.idx + 1 + (row.idx + 1) * mcols + (page.idx + 0) * mrows * mcols,
        col.idx + 0 + (row.idx + 0) * mcols + (page.idx + 1) * mrows * mcols,
        col.idx + 1 + (row.idx + 0) * mcols + (page.idx + 1) * mrows * mcols,
        col.idx + 0 + (row.idx + 1) * mcols + (page.idx + 1) * mrows * mcols,
        col.idx + 1 + (row.idx + 1) * mcols + (page.idx + 1) * mrows * mcols};

    const std::array w{page.fd[1] * row.fd[1] * col.fd[1],
                       page.fd[1] * row.fd[1] * col.fd[0],
                       page.fd[1] * row.fd[0] * col.fd[1],
                       page.fd[1] * row.fd[0] * col.fd[0],
                       page.fd[0] * row.fd[1] * col.fd[1],
                       page.fd[0] * row.fd[1] * col.fd[0],
                       page.fd[0] * row.fd[0] * col.fd[1],
                       page.fd[0] * row.fd[0] * col.fd[0]};

    bool any_base = false;
    base out;
    for (size_t i = 0; i < 8; i++) {
      if (w[i] not_eq 0) {
        if (any_base)
          out += w[i] * data[pos[i]];
        else
          out = w[i] * data[pos[i]];
        any_base = true;
      }
    }
    return out;
  }

  /** Number of pages */
  [[nodiscard]] size_t npages() const { return mpages; }

  /** Number of rows */
  [[nodiscard]] size_t nrows() const { return mrows; }

  /** Number of columns */
  [[nodiscard]] size_t ncols() const { return mcols; }

  friend std::ostream& operator<<(std::ostream& os, const Field3D& v) {
  for (size_t i = 0; i < v.npages(); i++)
    for (size_t j = 0; j < v.nrows(); j++)
      for (size_t k = 0; k < v.ncols(); k++) os << v(i, j, k) << '\n';
  return os;
}
};

#endif  // FIELD_HEADER
