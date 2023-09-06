/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   m_basic_types.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2002-05-08 

  \brief  Workspace functions for straightforward operations on variables 
          of basic types.

  This file includes workspace functions for variables of basic types, 
  such as Matrix and ArrayOfIndex. The functions are mainly of two types: <br>
  1. Initiation of variables by keyword arguments, such as *StringSet*. <br>
  2. Basic math, such as *VectorMatrixMultiply*.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>

#include "sorting.h"

#include <matpack.h>
#include <workspace.h>

inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric NAT_LOG_2=Constant::ln_2;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfIndexSetConstant(ArrayOfIndex& aoi,
                             const Index& nelem,
                             const Index& value) {
  aoi.resize(nelem);
  for (Index i = 0; i < nelem; i++) aoi[i] = value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfIndexLinSpace(ArrayOfIndex& x,
                          const Index& start,
                          const Index& stop,
                          const Index& step) {
  Index n = (Index)floor((stop - start) / step) + 1;
  if (n < 1) n = 1;

  x.resize(n);

  for (Index i = 0; i < n; i++) x[i] = start + i * step;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void FlagOff(Index& x) { x = 0; }

/* Workspace method: Doxygen documentation will be auto-generated */
void FlagOn(Index& x) { x = 1; }

/* Workspace method: Doxygen documentation will be auto-generated */
void IndexAdd(Index& out,
              const Index& in,
              const Index& value) {
  out = in + value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void IndexDivide(Index& out,
                 const Index& in,
                 const Index& value) {
  out = in / value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void IndexMultiply(Index& out,
                   const Index& in,
                   const Index& value) {
  out = in * value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void IndexStepDown(Index& xout, const Index& xin) {
  xout = xin - 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void IndexStepUp(Index& xout, const Index& xin) {
  xout = xin + 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void IndexSubtract(Index& out,
                   const Index& in,
                   const Index& value) {
  out = in - value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixAdd(Matrix& out,
               const Matrix& in,
               const Numeric& value) {
  // Note that in and out can be the same vector
  if (&out == &in) {
    // Out and in are the same. Just add the scalar value.
    out += value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then add the scalar value.
    out.resize(in.nrows(), in.ncols());
    out = in;
    out += value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixDivide(Matrix& out,
                  const Matrix& in,
                  const Numeric& value) {
  // Note that in and out can be the same matrix
  if (&out == &in) {
    // Out and in are the same. Just multiply by the scalar value.
    out /= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize(in.nrows(), in.ncols());
    out = in;
    out /= value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixGaussian(
    Matrix& Y,
    const Vector& x_row,
    const Numeric& x0_row,
    const Numeric& si_row,
    const Numeric& fwhm_row,
    const Vector& x_col,
    const Numeric& x0_col,
    const Numeric& si_col,
    const Numeric& fwhm_col)
{
  ARTS_USER_ERROR_IF ((si_row<=0 && fwhm_row<=0) ||
                      (si_row>0 && fwhm_row>0),
     "One of the GINs *si_row* and *fwhm_row* shall be >0, but just one.");
  ARTS_USER_ERROR_IF ((si_col<=0 && fwhm_col<=0) ||
                      (si_col>0 && fwhm_col>0),
     "One of the GINs *si_col* and *fwhm_col* shall be >0, but just one.");

  const Index nrow = x_row.nelem();
  const Index ncol = x_col.nelem();
  Y.resize(nrow, ncol);

  const Numeric si4row = si_row > 0 ? si_row : fwhm_row / (2 * sqrt(2 * NAT_LOG_2));
  const Numeric si4col = si_col > 0 ? si_col : fwhm_col / (2 * sqrt(2 * NAT_LOG_2));
  Vector row_term(nrow);
  for (Index r=0; r<nrow; ++r)
    row_term[r] = pow((x_row[r] - x0_row) / si4row, 2.0);
  Vector col_term(ncol);
  for (Index c=0; c<ncol; ++c)
    col_term[c] = pow((x_col[c] - x0_col) / si4col, 2.0);
  const Numeric fac = 1 / (2 * PI * si4row * si4col);
  for (Index r=0; r<nrow; ++r) {
     for (Index c=0; c<ncol; ++c) {
       Y(r,c) = fac * exp(-0.5 * (row_term[r] + col_term[c]));
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixMultiply(Matrix& out,
                    const Matrix& in,
                    const Numeric& value) {
  // Note that in and out can be the same matrix
  if (&out == &in) {
    // Out and in are the same. Just multiply by the scalar value.
    out *= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize(in.nrows(), in.ncols());
    out = in;
    out *= value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixSubtract(Matrix& out,
                    const Matrix& in,
                    const Numeric& value) {
  // Note that in and out can be the same matrix
  if (&out == &in) {
    // Out and in are the same. Just multiply by the scalar value.
    out -= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize(in.nrows(), in.ncols());
    out = in;
    out -= value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixCopySparse(Matrix& out, const Sparse& in) {
  // There is probably a more efficient way to do this
  out.resize(in.nrows(), in.ncols());
  for (Index r = 0; r < in.nrows(); r++) {
    for (Index c = 0; c < in.ncols(); c++) {
      out(r, c) = in(r, c);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixExtractFromTensor3(
    // WS Generic Output:
    Matrix& m,
    // WS Input:
    // WS Generic Input:
    const Tensor3& t3,
    const Index& index,
    // Control Parameters:
    const String& direction) {
  if (direction == "page") {
    if (index >= t3.npages()) {
      std::ostringstream os;
      os << "The index " << index
         << " is outside the page range of the Matrix.";
      throw std::runtime_error(os.str());
    }

    m.resize(t3.nrows(), t3.ncols());
    m = t3(index, joker, joker);
  } else if (direction == "row") {
    if (index >= t3.nrows()) {
      std::ostringstream os;
      os << "The index " << index << " is outside the row range of the Matrix.";
      throw std::runtime_error(os.str());
    }

    m.resize(t3.npages(), t3.ncols());
    m = t3(joker, index, joker);
  } else if (direction == "column") {
    if (index >= t3.ncols()) {
      std::ostringstream os;
      os << "The index " << index
         << " is outside the column range of the Matrix.";
      throw std::runtime_error(os.str());
    }

    m.resize(t3.npages(), t3.nrows());
    m = t3(joker, joker, index);
  } else {
    std::ostringstream os;
    os << "Keyword *direction* must be either *page* or *row* or *column*,"
       << "but you gave: " << direction << ".";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixMatrixMultiply(  // WS Generic Output:
    Matrix& Y,
    // WS Generic Input:
    const Matrix& M,
    const Matrix& X) {
  // Check that dimensions are right, M.ncols() must match X.nrows():
  if (M.ncols() != X.nrows()) {
    std::ostringstream os;
    os << "Matrix dimensions must be consistent!\n"
       << "Matrix1.ncols() = " << M.ncols() << "\n"
       << "Matrix2.nrows() = " << X.nrows();
    throw std::runtime_error(os.str());
  }

  // Temporary for the result:
  Matrix dummy(M.nrows(), X.ncols());

  mult(dummy, M, X);

  // Copy result to Y:

  Y.resize(dummy.nrows(), dummy.ncols());

  Y = dummy;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix1ColFromVector(  // WS Generic Output:
    Matrix& m,
    // WS Generic Input:
    const Vector& v) {
  const Index nv = v.nelem();

  m.resize(nv, 1);
  m(joker, 0) = v;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix2ColFromVectors(  // WS Generic Output:
    Matrix& m,
    // WS Generic Input:
    const Vector& v1,
    const Vector& v2) {
  const Index nv = v1.nelem();

  if (v2.nelem() != nv)
    throw std::runtime_error("Vectors must be of the same size.");

  m.resize(nv, 2);
  m(joker, 0) = v1;
  m(joker, 1) = v2;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix3ColFromVectors(  // WS Generic Output:
    Matrix& m,
    // WS Generic Input:
    const Vector& v1,
    const Vector& v2,
    const Vector& v3) {
  const Index nv = v1.nelem();

  if (v3.nelem() != nv || v2.nelem() != nv)
    throw std::runtime_error("Vectors must be of the same size.");

  m.resize(nv, 3);
  m(joker, 0) = v1;
  m(joker, 1) = v2;
  m(joker, 2) = v3;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix1RowFromVector(  // WS Generic Output:
    Matrix& m,
    // WS Generic Input:
    const Vector& v) {
  const Index nv = v.nelem();

  m.resize(1, nv);
  m(0, joker) = v;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix2RowFromVectors(  // WS Generic Output:
    Matrix& m,
    // WS Generic Input:
    const Vector& v1,
    const Vector& v2) {
  const Index nv = v1.nelem();

  if (v2.nelem() != nv)
    throw std::runtime_error("Vectors must be of the same size.");

  m.resize(2, nv);
  m(0, joker) = v1;
  m(1, joker) = v2;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix3RowFromVectors(  // WS Generic Output:
    Matrix& m,
    // WS Generic Input:
    const Vector& v1,
    const Vector& v2,
    const Vector& v3) {
  const Index nv = v1.nelem();

  if (v3.nelem() != nv || v2.nelem() != nv)
    throw std::runtime_error("Vectors must be of the same size.");

  m.resize(3, nv);
  m(0, joker) = v1;
  m(1, joker) = v2;
  m(2, joker) = v3;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixIdentity(Matrix& out,
                    const Index& n,
                    const Numeric& value) {
  out.resize(n, n);
  id_mat(out);
  if (value != 1) {
    out *= value;
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixReshapeTensor3(Matrix& m,
                          const Tensor3& t)
{
  const Index npages = t.npages();
  const Index nrows = t.nrows();
  const Index ncols = t.ncols();

  m.resize(npages * nrows, ncols);

  Index i = 0;
  for (Index p=0; p<npages; ++p) { 
    for (Index r=0; r<nrows; ++r) {
        m(i++, joker) = t(p, r, joker);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixSetConstant(Matrix& x,
                       const Index& nrows,
                       const Index& ncols,
                       const Numeric& value) {
  x.resize(nrows, ncols);
  x = value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void NumericAdd(Numeric& out,
                const Numeric& in,
                const Numeric& value) {
  out = in + value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void NumericClip(Numeric& out,
                 const Numeric& in,
                 const Numeric& limit_low,
                 const Numeric& limit_high) {
  if (in < limit_low)
    out = limit_low;
  else if (in > limit_high)
    out = limit_high;
  else
    out = in;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void NumericDivide(Numeric& out,
                   const Numeric& in,
                   const Numeric& value) {
  out = in / value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void NumericFromVector(Numeric& out,
                       const Vector& in,
                       const String& op) {
  if (op == "first")
    out = in[0];
  else if (op == "last")
    out = in[in.nelem() - 1];
  else if (op == "max")
    out = max(in);
  else if (op == "min")
    out = min(in);
  else if (op == "mean")
    out = mean(in);
  else {
    std::ostringstream os;
    os << "Your choice, *op* = \"" << op << "\", is not recognised.\n"
       << R"(Valid options are: "first", "last", "max", "min" and "mean".)";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void NumericMultiply(Numeric& out,
                     const Numeric& in,
                     const Numeric& value) {
  out = in * value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void NumericSubtract(Numeric& out,
                     const Numeric& in,
                     const Numeric& value) {
  out = in - value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void RationalAdd(Rational& out,
                 const Rational& in,
                 const Rational& value) {
  out = in + value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void RationalDivide(Rational& out,
                    const Rational& in,
                    const Rational& value) {
  out = in / value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void RationalMultiply(Rational& out,
                      const Rational& in,
                      const Rational& value) {
  out = in * value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void RationalSubtract(Rational& out,
                      const Rational& in,
                      const Rational& value) {
  out = in - value;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void SparseSparseMultiply(  // WS Generic Output:
    Sparse& Y,
    // WS Generic Input:
    const Sparse& M,
    const Sparse& X) {
  // Check that dimensions are right, M.ncols() must match X.nrows():
  if (M.ncols() != X.nrows()) {
    std::ostringstream os;
    os << "Matrix dimensions must be consistent!\n"
       << "Matrix1.ncols() = " << M.ncols() << "\n"
       << "Matrix2.nrows() = " << X.nrows();
    throw std::runtime_error(os.str());
  }

  // Temporary for the result:
  Sparse dummy(M.nrows(), X.ncols());

  mult(dummy, M, X);

  // Copy result to Y:
  Y = dummy;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void SparseIdentity(Sparse& X,
                    const Index& n,
                    const Numeric& value) {
  X.resize(n, n);
  id_mat(X);

  if (value != 1.0) X *= value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void DiagonalMatrix(Matrix& X, const Vector& diag /*v*/) {
  Index n = diag.nelem();
  X.resize(n, n);
  X = 0.0;

  for (Index i = 0; i < n; ++i) {
    X(i, i) = diag[i];
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void DiagonalMatrix(Sparse& X, const Vector& diag /*v*/) {
  Index n = diag.nelem();
  X.resize(n, n);

  ArrayOfIndex indices(n);

  for (Index i = 0; i < n; ++i) {
    indices[i] = i;
  }

  X.insert_elements(n, indices, indices, diag);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor3Add(Tensor3& out,
                const Tensor3& in,
                const Numeric& value) {
  // Note that in and out can be the same vector
  if (&out == &in) {
    // Out and in are the same. Just multiply by the scalar value.
    out += value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize(in.npages(), in.nrows(), in.ncols());
    out = in;
    out += value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor3FromVector(  // WS Generic Output:
    Tensor3& m,
    // WS Generic Input:
    const Vector& v) {
  const Index nv = v.nelem();

  m.resize(nv, 1, 1);
  m(joker, 0, 0) = v;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor3Multiply(Tensor3& out,
                     const Tensor3& in,
                     const Numeric& value) {
  // Note that in and out can be the same vector
  if (&out == &in) {
    // Out and in are the same. Just multiply by the scalar value.
    out *= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize(in.npages(), in.nrows(), in.ncols());
    out = in;
    out *= value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor3SetConstant(Tensor3& x,
                        const Index& npages,
                        const Index& nrows,
                        const Index& ncols,
                        const Numeric& value) {
  x.resize(npages, nrows, ncols);
  x = value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor3ExtractFromTensor4(
    // WS Generic Output:
    Tensor3& t3,
    // WS Input:
    // WS Generic Input:
    const Tensor4& t4,
    const Index& index,
    // Control Parameters:
    const String& direction) {
  if (direction == "book") {
    if (index >= t4.nbooks()) {
      std::ostringstream os;
      os << "The index " << index
         << " is outside the book range of the Tensor4.";
      throw std::runtime_error(os.str());
    }

    t3.resize(t4.npages(), t4.nrows(), t4.ncols());
    t3 = t4(index, joker, joker, joker);
  } else if (direction == "page") {
    if (index >= t4.npages()) {
      std::ostringstream os;
      os << "The index " << index
         << " is outside the pages range of the Tensor4.";
      throw std::runtime_error(os.str());
    }

    t3.resize(t4.nbooks(), t4.nrows(), t4.ncols());
    t3 = t4(joker, index, joker, joker);
  } else if (direction == "row") {
    if (index >= t4.nrows()) {
      std::ostringstream os;
      os << "The index " << index
         << " is outside the row range of the Tensor4.";
      throw std::runtime_error(os.str());
    }

    t3.resize(t4.npages(), t4.nbooks(), t4.ncols());
    t3 = t4(joker, joker, index, joker);
  } else if (direction == "column") {
    if (index >= t4.ncols()) {
      std::ostringstream os;
      os << "The index " << index
         << " is outside the column range of the Tensor4.";
      throw std::runtime_error(os.str());
    }

    t3.resize(t4.npages(), t4.nbooks(), t4.nrows());
    t3 = t4(joker, joker, joker, index);
  } else {
    std::ostringstream os;
    os << "Keyword *direction* must be either *page*, *book*, *row* or *column*,"
       << "but you gave: " << direction << ".";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor4Add(Tensor4& out,
                const Tensor4& in,
                const Numeric& value) {
  // Note that in and out can be the same vector
  if (&out == &in) {
    // Out and in are the same. Just multiply by the scalar value.
    out += value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize(in.nbooks(), in.npages(), in.nrows(), in.ncols());
    out = in;
    out += value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor4Multiply(Tensor4& out,
                     const Tensor4& in,
                     const Numeric& value) {
  // Note that in and out can be the same vector
  if (&out == &in) {
    // Out and in are the same. Just multiply by the scalar value.
    out *= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize(in.nbooks(), in.npages(), in.nrows(), in.ncols());
    out = in;
    out *= value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor4SetConstant(Tensor4& x,
                        const Index& nbooks,
                        const Index& npages,
                        const Index& nrows,
                        const Index& ncols,
                        const Numeric& value) {
  x.resize(nbooks, npages, nrows, ncols);
  x = value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor5Multiply(Tensor5& out,
                     const Tensor5& in,
                     const Numeric& value) {
  // Note that in and out can be the same vector
  if (&out == &in) {
    // Out and in are the same. Just multiply by the scalar value.
    out *= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize(in.nshelves(), in.nbooks(), in.npages(), in.nrows(), in.ncols());
    out = in;
    out *= value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor5SetConstant(Tensor5& x,
                        const Index& nshelves,
                        const Index& nbooks,
                        const Index& npages,
                        const Index& nrows,
                        const Index& ncols,
                        const Numeric& value) {
  x.resize(nshelves, nbooks, npages, nrows, ncols);
  x = value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor6Multiply(Tensor6& out,
                     const Tensor6& in,
                     const Numeric& value) {
  // Note that in and out can be the same vector
  if (&out == &in) {
    // Out and in are the same. Just multiply by the scalar value.
    out *= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize(in.nvitrines(),
               in.nshelves(),
               in.nbooks(),
               in.npages(),
               in.nrows(),
               in.ncols());
    out = in;
    out *= value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor6SetConstant(Tensor6& x,
                        const Index& nvitrines,
                        const Index& nshelves,
                        const Index& nbooks,
                        const Index& npages,
                        const Index& nrows,
                        const Index& ncols,
                        const Numeric& value) {
  x.resize(nvitrines, nshelves, nbooks, npages, nrows, ncols);
  x = value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor7Multiply(Tensor7& out,
                     const Tensor7& in,
                     const Numeric& value) {
  // Note that in and out can be the same vector
  if (&out == &in) {
    // Out and in are the same. Just multiply by the scalar value.
    out *= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize(in.nlibraries(),
               in.nvitrines(),
               in.nshelves(),
               in.nbooks(),
               in.npages(),
               in.nrows(),
               in.ncols());
    out = in;
    out *= value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor7SetConstant(Tensor7& x,
                        const Index& nlibraries,
                        const Index& nvitrines,
                        const Index& nshelves,
                        const Index& nbooks,
                        const Index& npages,
                        const Index& nrows,
                        const Index& ncols,
                        const Numeric& value) {
  x.resize(nlibraries, nvitrines, nshelves, nbooks, npages, nrows, ncols);
  x = value;
}

/* Workspace method: Doxygen documentation will be auto-generated  */
void Trapz(
    Numeric& out,
    const Vector& x,
    const Vector& y) {
  const Index n = x.nelem();
  if (y.nelem() != n) 
    throw std::runtime_error("The vectors *x* and *y* must have the same length.");
    
  out = 0;
  for (Index i=1; i<n; i++)
    out += 0.5*(y[i-1]+y[i]) * (x[i]-x[i-1]);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorAdd(Vector& out,
               const Vector& in,
               const Numeric& value) {
  // Note that in and out can be the same vector
  if (&out == &in) {
    // Out and in are the same. Just add the scalar value.
    out += value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then add the scalar value.
    out.resize(in.nelem());
    out = in;
    out += value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorDivide(Vector& out,
                  const Vector& in,
                  const Numeric& value) {
  // Note that in and out can be the same vector
  if (&out == &in) {
    // Out and in are the same. Just add the scalar value.
    out /= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then add the scalar value.
    out.resize(in.nelem());
    out = in;
    out /= value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorMultiply(Vector& out,
                    const Vector& in,
                    const Numeric& value) {
  // Note that in and out can be the same vector
  if (&out == &in) {
    // Out and in are the same. Just multiply by the scalar value.
    out *= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize(in.nelem());
    out = in;
    out *= value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorSubtract(Vector& out,
                    const Vector& in,
                    const Numeric& value) {
  // Note that in and out can be the same vector
  if (&out == &in) {
    // Out and in are the same. Just add the scalar value.
    out -= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then add the scalar value.
    out.resize(in.nelem());
    out = in;
    out -= value;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorAddElementwise(Vector& c,
                          const Vector& a,
                          const Vector& b) {
  // b has length 1. Here we easily can avoid just adding 0.
  if (b.nelem() == 1) {
    // a and c are the same WSV
    if (&c == &a) {
      if (b[0] != 0) {
        c += b[0];
      }
    } else {
      c = a;
      if (b[0] != 0) {
        c += b[0];
      }
    }
  }

  // b is a vector
  else if (b.nelem() == a.nelem()) {
    // a and c are the same WSV
    if (&c == &a) {
      c += b;
    } else {
      c = a;
      c += b;
    }
  }

  else
    throw std::runtime_error(
        "The vector *b* must have length 1 or match *a* in length.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorSubtractElementwise(Vector& c,
                          const Vector& a,
                          const Vector& b) {
  // b has length 1. Here we easily can avoid just adding 0.
  if (b.nelem() == 1) {
    // a and c are the same WSV
    if (&c == &a) {
      if (b[0] != 0) {
        c -= b[0];
      }
    } else {
      c = a;
      if (b[0] != 0) {
        c -= b[0];
      }
    }
  }

  // b is a vector
  else if (b.nelem() == a.nelem()) {
    // a and c are the same WSV
    if (&c == &a) {
      c -= b;
    } else {
      c = a;
      c -= b;
    }
  }

  else
    throw std::runtime_error(
        "The vector *b* must have length 1 or match *a* in length.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorMultiplyElementwise(Vector& c,
                               const Vector& a,
                               const Vector& b) {
  // b has length 1. Here we easily can avoid just adding 0.
  if (b.nelem() == 1) {
    // a and c are the same WSV
    if (&c == &a) {
      if (b[0] != 0) {
        c *= b[0];
      }
    } else {
      c = a;
      if (b[0] != 0) {
        c *= b[0];
      }
    }
  }

  // b is a vector
  else if (b.nelem() == a.nelem()) {
    // a and c are the same WSV
    if (&c == &a) {
      c *= b;
    } else {
      c = a;
      c *= b;
    }
  }

  else
    throw std::runtime_error(
        "The vector *b* must have length 1 or match *a* in length.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorDivideElementwise(Vector& c,
                             const Vector& a,
                             const Vector& b) {
  // b has length 1. Here we easily can avoid just adding 0.
  if (b.nelem() == 1) {
    // a and c are the same WSV
    if (&c == &a) {
      if (b[0] != 0) {
        c /= b[0];
      }
    } else {
      c = a;
      if (b[0] != 0) {
        c /= b[0];
      }
    }
  }

  // b is a vector
  else if (b.nelem() == a.nelem()) {
    // a and c are the same WSV
    if (&c == &a) {
      c /= b;
    } else {
      c = a;
      c /= b;
    }
  }

  else
    throw std::runtime_error(
        "The vector *b* must have length 1 or match *a* in length.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorClip(Vector& out,
                const Vector& in,
                const Numeric& limit_low,
                const Numeric& limit_high) {
  const Index l = in.nelem();
  if (out.nelem() != l)
    out.resize(l);
  
  for (Index i=0; i<l; i++) {
    if (in[i] < limit_low)
      out[i] = limit_low;
   else if (in[i] > limit_high)
     out[i] = limit_high;
   else
     out[i] = in[i];
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorCrop(Vector& out,
                const Vector& in,
                const Numeric& min_value,
                const Numeric& max_value) {
  const Index nin = in.nelem();

  Index nout = 0;
  //
  for (Index i = 0; i < nin; i++) {
    if (in[i] >= min_value && in[i] <= max_value) {
      nout += 1;
    }
  }

  // Make copy if in-vector, as it also can be the out one
  Vector c(in);

  out.resize(nout);

  nout = 0;
  //
  for (Index i = 0; i < nin; i++) {
    if (c[i] >= min_value && c[i] <= max_value) {
      out[nout] = c[i];
      nout += 1;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated 

   2004-09-15 Patrick Eriksson

   Added keyword to control if row or column is extracted.

   2007-07-24 Stefan Buehler */
void VectorExtractFromMatrix(
    // WS Generic Output:
    Vector& v,
    // WS Input:
    // WS Generic Input:
    const Matrix& m,
    const Index& index,
    // Control Parameters:
    const String& direction) {
  if (direction == "row") {
    if (index >= m.nrows()) {
      std::ostringstream os;
      os << "The index " << index << " is outside the row range of the Matrix.";
      throw std::runtime_error(os.str());
    }

    v.resize(m.ncols());
    v = m(index, joker);
  } else if (direction == "column") {
    if (index >= m.ncols()) {
      std::ostringstream os;
      os << "The index " << index
         << " is outside the column range of the Matrix.";
      throw std::runtime_error(os.str());
    }

    v.resize(m.nrows());
    v = m(joker, index);
  } else {
    std::ostringstream os;
    os << "Keyword *direction* must be either *row* or *column*,"
       << "but you gave: " << direction << ".";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorFlip(Vector& out, const Vector& in) {
  const Index n = in.nelem();

  // Note that in and out can be the same vector
  if (&out == &in) {
    // Out and in are the same. A copy is needed
    const Vector v = in;
    for (Index i = 0; i < n; i++) out[i] = v[n - 1 - i];
  } else {
    // Out and in are different.
    out.resize(n);
    for (Index i = 0; i < n; i++) out[i] = in[n - 1 - i];
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorGaussian(
    Vector& y,
    const Vector& x,
    const Numeric& x0,
    const Numeric& si,
    const Numeric& fwhm)
{
  ARTS_USER_ERROR_IF ((si<=0 && fwhm<=0) || (si>0 && fwhm>0),
     "One of the GINs *si* and *fwhm* shall be >0, but just one.");

  const Index n = x.nelem();

  // Note that y and x can be the same vector
  if (&y != &x) {
    y.resize(n);
  }

  const Numeric si2use = si > 0 ? si : fwhm / (2 * sqrt(2 * NAT_LOG_2));
  const Numeric fac = 1 / (sqrt(2 * PI) * si2use);
  for (Index i=0; i<n; ++i) {
    y[i] = fac * exp(-0.5 * pow((x[i] - x0) / si2use, 2.0));
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorInsertGridPoints(  // WS Generic Output:
    Vector& og,               // Output grid
    // WS Generic Input:
    const Vector& ingrid,  // Input grid
    const Vector& points) {
  // First make duplicates of the input vectors, in case one of them
  // happens to be identical to the output vector. Also, we can fool
  // around with these, if we want.
  Vector ig(ingrid);
  Vector p(points);

  // Check how the input grid is sorted. If the grid is sorted in
  // descending order, we simply turn it around. (But don't
  // forget to turn it back at the end!)
  Index ascending;  // 1=ascending, 0=descending
  if (is_increasing(ig)) {
    ascending = 1;
  } else if (is_decreasing(ig)) {
    ascending = 0;

    // Turn grid round.

    // Copy ig to dummy vector in reverse order:
    const Vector dummy = reverse(ig);

    // Copy dummy back to ig vector:
    ig = dummy;
  } else {
    std::ostringstream os;
    os << "The input Vector must be either\n"
       << "strictly increasing or strictly decreasing,\n"
       << "but this is not the case.\n";
    os << "The vector contains:\n" << ig;
    throw std::runtime_error(os.str());
  }

  // Sort also the vector of points to insert in increasing order:
  {
    ArrayOfIndex si;            // Sorted indices
    get_sorted_indexes(si, p);  // Get sorted p indices
    const Vector dummy = p;     // Copy p to dummy
    // Copy back dummy to p in right order:
    for (Index j = 0; j < p.nelem(); j++) p[j] = dummy[si[j]];
  }

  // The idea is to step through both ig and p, and build up the
  // output in a temporary array.
  Array<Numeric> x;
  Index iig = 0, ip = 0;  // indices to ig and p
  while (iig < ig.nelem() && ip < p.nelem()) {
    if (p[ip] < ig[iig]) {
      x.push_back(p[ip]);
      ++ip;
    } else if (p[ip] > ig[iig]) {
      x.push_back(ig[iig]);
      ++iig;
    } else {
      ++ip;
    }
  }

  // Add remaining points of either p or ig, depending on which is
  // longer:
  if (ip == p.nelem()) {
    // p has reached its end.
    while (iig < ig.nelem()) {
      x.push_back(ig[iig]);
      ++iig;
    }
  } else if (iig == ig.nelem()) {
    // ig has reached its end
    while (ip < p.nelem()) {
      x.push_back(p[ip]);
      ++ip;
    }
  } else {
    // We should never be here.
    ARTS_ASSERT(false);
    std::exit(EXIT_FAILURE);
  }

  // Ok, x should now contain the new grid.

  og.resize(x.nelem());

  // Copy to result vector, turn around if necessary.
  if (ascending)
    for (Index i = 0; i < x.nelem(); ++i) og[i] = x[i];  // Just copy.
  else
    for (Index i = 0; i < x.nelem(); ++i)
      og[i] = x[x.nelem() - 1 - i];  // Copy in reverse order.
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorLinSpace(Vector& x,
                    const Numeric& start,
                    const Numeric& stop,
                    const Numeric& step) {
  linspace(x, start, stop, step);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorLogSpace(Vector& x,
                    const Numeric& start,
                    const Numeric& stop,
                    const Numeric& step) {
  linspace(x, log(start), log(stop), step);
  transform(x, exp, x);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorMatrixMultiply(  // WS Generic Output:
    Vector& y,
    // WS Generic Input:
    const Matrix& M,
    const Vector& x) {
  // Check that dimensions are right, x must match columns of M:
  if (M.ncols() != x.nelem()) {
    std::ostringstream os;
    os << "Matrix and vector dimensions must be consistent!\n"
       << "Matrix.ncols() = " << M.ncols() << "\n"
       << "Vector.nelem() = " << x.nelem();
    throw std::runtime_error(os.str());
  }

  // Temporary for the result:
  Vector dummy(M.nrows());

  mult(dummy, M, x);

  y.resize(dummy.nelem());

  y = dummy;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorSparseMultiply(  // WS Generic Output:
    Vector& y,
    // WS Generic Input:
    const Sparse& M,
    const Vector& x) {
  // Check that dimensions are right, x must match columns of M:
  if (M.ncols() != x.nelem()) {
    std::ostringstream os;
    os << "Sparse and vector dimensions must be consistent!\n"
       << "Sparse.ncols() = " << M.ncols() << "\n"
       << "Vector.nelem() = " << x.nelem();
    throw std::runtime_error(os.str());
  }

  // Temporary for the result:
  Vector dummy(M.nrows());

  mult(dummy, M, x);

  y.resize(dummy.nelem());

  y = dummy;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorNLinSpace(Vector& x,
                     const Index& n,
                     const Numeric& start,
                     const Numeric& stop) {
  if (n < 2) throw std::runtime_error("The number of points must be > 1.");
  nlinspace(x, start, stop, n);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorNLinSpaceVector(Vector& x,
                     const Index& n,
                     const Vector& y) {
  VectorNLinSpace(x, n, y[0], last(y));
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfTimeNLinSpace(ArrayOfTime& x,
                          const Index& n,
                          const String& start,
                          const String& stop) {
  Vector seconds;
  Time t0(start), t1(stop);

  ARTS_USER_ERROR_IF (n < 2, "The number of points must be > 1.");
  nlinspace(seconds, t0.Seconds(), t1.Seconds(), n);
  
  x = time_vector(seconds);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorNLogSpace(Vector& x,
                     const Index& n,
                     const Numeric& start,
                     const Numeric& stop) {
  if (n < 2) throw std::runtime_error("The number of points must be > 1.");
  if ((start <= 0) || (stop <= 0))
    throw std::runtime_error("Only positive numbers are allowed.");

  nlogspace(x, start, stop, n);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorPower(Vector& out,
                 const Vector& in,
                 const Numeric& value) {
  const Index n = in.nelem();
  // Note that in and out can be the same vector
  if (&out == &in) {
    // Out and in are the same. 
  } else {
    out.resize(n);
  }
  for (Index i=0; i<n; i++) {
    out[i] = pow( in[i], value );
  }  
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorReshapeMatrix(Vector& v,
                         const Matrix& m,
                         const String& direction) {
  const Index nrows = m.nrows();
  const Index ncols = m.ncols();

  v.resize(nrows * ncols);

  Index iv = 0;

  if (direction == "column") {
    for (Index col = 0; col < ncols; col++) {
      for (Index row = 0; row < nrows; row++) {
        v[iv] = m(row, col);
        iv++;
      }
    }
  } else if (direction == "row") {
    for (Index row = 0; row < nrows; row++) {
      for (Index col = 0; col < ncols; col++) {
        v[iv] = m(row, col);
        iv++;
      }
    }
  } else {
    std::ostringstream os;
    os << "Keyword *direction* must be either *row* or *column*,"
       << "but you gave: " << direction << ".";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void VectorSetConstant(Vector& x,
                       const Index& n,
                       const Numeric& value) {
  x.resize(n);
  x = value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfTimeSetConstant(ArrayOfTime& x,
                       const Index& n,
                       const Time& value) {
  x.resize(n);
  x = value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Numeric& var1,
             const Numeric& var2,
             const Numeric& maxabsdiff,
             const String& error_message) {
  Numeric maxdiff = var1 - var2;

  if (std::isnan(var1) || std::isnan(var2)) {
    if (std::isnan(var1) && std::isnan(var2)) {
      maxdiff = 0;
    } else if (std::isnan(var1)) {
      std::ostringstream os;
      os << "Nan found in var1, but there is no "
         << "NaN at same position in var2.\nThis "
         << "is not allowed.";
      throw std::runtime_error(os.str());
    } else {
      std::ostringstream os;
      os << "Nan found in var2, but there is no "
         << "NaN at same position in var1.\nThis "
         << "is not allowed.";
      throw std::runtime_error(os.str());
    }
  }

  if (abs(maxdiff) > maxabsdiff) {
    std::ostringstream os;
    os << "var1-var2 FAILED!\n";
    if (error_message.length()) os << error_message << "\n";
    os << "Max allowed deviation set to: " << maxabsdiff << std::endl
       << "but the value deviates with:  " << maxdiff << std::endl;
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Vector& var1,
             const Vector& var2,
             const Numeric& maxabsdiff,
             const String& error_message) {
  const Index n = var1.nelem();

  if (var2.nelem() != n) {
    std::ostringstream os;
    os << "var1 (" << n << ") and var2 (" << var2.nelem()
       << ") do not have the same size.";
    throw std::runtime_error(os.str());
  }

  Numeric maxdiff = 0.0;
  for (Index i = 0; i < n; i++) {
    Numeric diff = var1[i] - var2[i];

    if (std::isnan(var1[i]) || std::isnan(var2[i])) {
      if (std::isnan(var1[i]) && std::isnan(var2[i])) {
        diff = 0;
      } else if (std::isnan(var1[i])) {
        std::ostringstream os;
        os << "Nan found in var1, but there is no "
           << "NaN at same position in var2.\nThis "
           << "is not allowed.";
        throw std::runtime_error(os.str());
      } else {
        std::ostringstream os;
        os << "Nan found in var2, but there is no "
           << "NaN at same position in var1.\nThis "
           << "is not allowed.";
        throw std::runtime_error(os.str());
      }
    }

    if (abs(diff) > abs(maxdiff)) {
      maxdiff = diff;
    }
  }

  if (std::isnan(maxdiff) || abs(maxdiff) > maxabsdiff) {
    std::ostringstream os;
    os << "var1-var2 FAILED!\n";
    if (error_message.length()) os << error_message << "\n";
    os << "Max allowed deviation set to: " << maxabsdiff << std::endl
       << "but the vectors deviate with: " << maxdiff << std::endl;
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Matrix& var1,
             const Matrix& var2,
             const Numeric& maxabsdiff,
             const String& error_message) {
  const Index nrows = var1.nrows();
  const Index ncols = var1.ncols();

  if (var2.nrows() != nrows || var2.ncols() != ncols) {
    std::ostringstream os;
    os << "var1 (" << nrows << "," << ncols << ") and var2"
       << " (" << var2.nrows() << "," << var2.ncols()
       << ") do not have the same size.";
    throw std::runtime_error(os.str());
  }

  Numeric maxdiff = 0.0;

  for (Index r = 0; r < nrows; r++) {
    for (Index c = 0; c < ncols; c++) {
      Numeric diff = var1(r, c) - var2(r, c);

      if (std::isnan(var1(r, c)) || std::isnan(var2(r, c))) {
        if (std::isnan(var1(r, c)) && std::isnan(var2(r, c))) {
          diff = 0;
        } else if (std::isnan(var1(r, c))) {
          std::ostringstream os;
          os << "Nan found in var1, but there is no "
             << "NaN at same position in var2.\nThis "
             << "is not allowed.";
          throw std::runtime_error(os.str());
        } else {
          std::ostringstream os;
          os << "Nan found in var2, but there is no "
             << "NaN at same position in var1.\nThis "
             << "is not allowed.";
          throw std::runtime_error(os.str());
        }
      }

      if (abs(diff) > abs(maxdiff)) {
        maxdiff = diff;
      }
    }
  }

  if (abs(maxdiff) > maxabsdiff) {
    std::ostringstream os;
    os << "var1-var2 FAILED!\n";
    if (error_message.length()) os << error_message << "\n";
    os << "Max allowed deviation set to : " << maxabsdiff << std::endl
       << "but the matrices deviate with: " << maxdiff << std::endl;
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Tensor3& var1,
             const Tensor3& var2,
             const Numeric& maxabsdiff,
             const String& error_message) {
  const Index ncols = var1.ncols();
  const Index nrows = var1.nrows();
  const Index npages = var1.npages();

  if (var2.ncols() != ncols || var2.nrows() != nrows ||
      var2.npages() != npages) {
    std::ostringstream os;
    os << "var1 and var2 do not have the same size.";
    throw std::runtime_error(os.str());
  }

  Numeric maxdiff = 0.0;

  for (Index c = 0; c < ncols; c++)
    for (Index r = 0; r < nrows; r++)
      for (Index p = 0; p < npages; p++) {
        Numeric diff = var1(p, r, c) - var2(p, r, c);

        if (std::isnan(var1(p, r, c)) || std::isnan(var2(p, r, c))) {
          if (std::isnan(var1(p, r, c)) && std::isnan(var2(p, r, c))) {
            diff = 0;
          } else if (std::isnan(var1(p, r, c))) {
            std::ostringstream os;
            os << "Nan found in var1, but there is no "
               << "NaN at same position in var2.\nThis "
               << "is not allowed.";
            throw std::runtime_error(os.str());
          } else {
            std::ostringstream os;
            os << "Nan found in var2, but there is no "
               << "NaN at same position in var1.\nThis "
               << "is not allowed.";
            throw std::runtime_error(os.str());
          }
        }

        if (abs(diff) > abs(maxdiff)) {
          maxdiff = diff;
        }
      }

  if (abs(maxdiff) > maxabsdiff) {
    std::ostringstream os;
    os << "var1-var2 FAILED!\n";
    if (error_message.length()) os << error_message << "\n";
    os << "Max allowed deviation set to : " << maxabsdiff << std::endl
       << "but the tensors deviate with: " << maxdiff << std::endl;
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Tensor4& var1,
             const Tensor4& var2,
             const Numeric& maxabsdiff,
             const String& error_message) {
  const Index ncols = var1.ncols();
  const Index nrows = var1.nrows();
  const Index npages = var1.npages();
  const Index nbooks = var1.nbooks();

  if (var2.ncols() != ncols || var2.nrows() != nrows ||
      var2.npages() != npages || var2.nbooks() != nbooks) {
    std::ostringstream os;
    os << "var1 and var2 do not have the same size.";
    throw std::runtime_error(os.str());
  }

  Numeric maxdiff = 0.0;

  for (Index c = 0; c < ncols; c++)
    for (Index r = 0; r < nrows; r++)
      for (Index p = 0; p < npages; p++)
        for (Index b = 0; b < nbooks; b++) {
          Numeric diff = var1(b, p, r, c) - var2(b, p, r, c);

          if (std::isnan(var1(b, p, r, c)) || std::isnan(var2(b, p, r, c))) {
            if (std::isnan(var1(b, p, r, c)) && std::isnan(var2(b, p, r, c))) {
              diff = 0;
            } else if (std::isnan(var1(b, p, r, c))) {
              std::ostringstream os;
              os << "Nan found in var1, but there is no "
                 << "NaN at same position in var2.\nThis "
                 << "is not allowed.";
              throw std::runtime_error(os.str());
            } else {
              std::ostringstream os;
              os << "Nan found in var2, but there is no "
                 << "NaN at same position in var1.\nThis "
                 << "is not allowed.";
              throw std::runtime_error(os.str());
            }
          }

          if (abs(diff) > abs(maxdiff)) {
            maxdiff = diff;
          }
        }

  if (abs(maxdiff) > maxabsdiff) {
    std::ostringstream os;
    os << "var1-var2 FAILED!\n";
    if (error_message.length()) os << error_message << "\n";
    os << "Max allowed deviation set to : " << maxabsdiff << std::endl
       << "but the tensors deviate with: " << maxdiff << std::endl;
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Tensor5& var1,
             const Tensor5& var2,
             const Numeric& maxabsdiff,
             const String& error_message) {
  const Index ncols = var1.ncols();
  const Index nrows = var1.nrows();
  const Index npages = var1.npages();
  const Index nbooks = var1.nbooks();
  const Index nshelves = var1.nshelves();

  if (var2.ncols() != ncols || var2.nrows() != nrows ||
      var2.npages() != npages || var2.nbooks() != nbooks ||
      var2.nshelves() != nshelves) {
    std::ostringstream os;
    os << "var1 and var2 do not have the same size.";
    throw std::runtime_error(os.str());
  }

  Numeric maxdiff = 0.0;

  for (Index c = 0; c < ncols; c++)
    for (Index r = 0; r < nrows; r++)
      for (Index p = 0; p < npages; p++)
        for (Index b = 0; b < nbooks; b++)
          for (Index s = 0; s < nshelves; s++) {
            Numeric diff = var1(s, b, p, r, c) - var2(s, b, p, r, c);

            if (std::isnan(var1(s, b, p, r, c)) ||
                std::isnan(var2(s, b, p, r, c))) {
              if (std::isnan(var1(s, b, p, r, c)) &&
                  std::isnan(var2(s, b, p, r, c))) {
                diff = 0;
              } else if (std::isnan(var1(s, b, p, r, c))) {
                std::ostringstream os;
                os << "Nan found in var1, but there is no "
                   << "NaN at same position in var2.\nThis "
                   << "is not allowed.";
                throw std::runtime_error(os.str());
              } else {
                std::ostringstream os;
                os << "Nan found in var2, but there is no "
                   << "NaN at same position in var1.\nThis "
                   << "is not allowed.";
                throw std::runtime_error(os.str());
              }
            }

            if (abs(diff) > abs(maxdiff)) {
              maxdiff = diff;
            }
          }

  if (abs(maxdiff) > maxabsdiff) {
    std::ostringstream os;
    os << "var1-var2 FAILED!\n";
    if (error_message.length()) os << error_message << "\n";
    os << "Max allowed deviation set to : " << maxabsdiff << std::endl
       << "but the tensors deviate with: " << maxdiff << std::endl;
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Tensor7& var1,
             const Tensor7& var2,
             const Numeric& maxabsdiff,
             const String& error_message) {
  const Index ncols = var1.ncols();
  const Index nrows = var1.nrows();
  const Index npages = var1.npages();
  const Index nbooks = var1.nbooks();
  const Index nshelves = var1.nshelves();
  const Index nvitrines = var1.nvitrines();
  const Index nlibraries = var1.nlibraries();

  if (var2.ncols() != ncols || var2.nrows() != nrows ||
      var2.npages() != npages || var2.nbooks() != nbooks ||
      var2.nshelves() != nshelves || var2.nvitrines() != nvitrines ||
      var2.nlibraries() != nlibraries) {
    std::ostringstream os;
    os << "var1 and var2 do not have the same size.";
    throw std::runtime_error(os.str());
  }

  Numeric maxdiff = 0.0;

  for (Index c = 0; c < ncols; c++)
    for (Index r = 0; r < nrows; r++)
      for (Index p = 0; p < npages; p++)
        for (Index b = 0; b < nbooks; b++)
          for (Index s = 0; s < nshelves; s++)
            for (Index v = 0; v < nvitrines; v++)
              for (Index l = 0; l < nlibraries; l++) {
                Numeric diff =
                    var1(l, v, s, b, p, r, c) - var2(l, v, s, b, p, r, c);

                if (std::isnan(var1(l, v, s, b, p, r, c)) ||
                    std::isnan(var2(l, v, s, b, p, r, c))) {
                  if (std::isnan(var1(l, v, s, b, p, r, c)) &&
                      std::isnan(var2(l, v, s, b, p, r, c))) {
                    diff = 0;
                  } else if (std::isnan(var1(l, v, s, b, p, r, c))) {
                    std::ostringstream os;
                    os << "Nan found in var1, but there is no "
                       << "NaN at same position in var2.\nThis "
                       << "is not allowed.";
                    throw std::runtime_error(os.str());
                  } else {
                    std::ostringstream os;
                    os << "Nan found in var2, but there is no "
                       << "NaN at same position in var1.\nThis "
                       << "is not allowed.";
                    throw std::runtime_error(os.str());
                  }
                }

                if (abs(diff) > abs(maxdiff)) {
                  maxdiff = diff;
                }
              }

  if (abs(maxdiff) > maxabsdiff) {
    std::ostringstream os;
    os << "var1-var2 FAILED!\n";
    if (error_message.length()) os << error_message << "\n";
    os << "Max allowed deviation set to : " << maxabsdiff << std::endl
       << "but the tensors deviate with: " << maxdiff << std::endl;
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const ArrayOfVector& var1,
             const ArrayOfVector& var2,
             const Numeric& maxabsdiff,
             const String& error_message) {
  if (var1.nelem() != var2.nelem()) {
    std::ostringstream os;
    os << "The two arrays do not have the same size." << std::endl
       << "var1 nelem: " << var1.nelem() << std::endl
       << "var2" " nelem: " << var2.nelem() << std::endl;
    throw std::runtime_error(os.str());
  }

  bool failed = false;
  std::ostringstream fail_msg;
  for (Index i = 0; i < var1.nelem(); i++) {
    try {
      std::ostringstream vn1, vn2;
      vn1 << "var1[" << i << "]";
      vn2 << "var2" "[" << i << "]";
      Compare(var1[i],
              var2[i],
              maxabsdiff,
              error_message);
    } catch (const std::runtime_error& e) {
      failed = true;
      fail_msg << std::endl
               << e.what() << std::endl
               << "Mismatch at array index: " << i << std::endl;
    }
  }

  if (failed) throw std::runtime_error(fail_msg.str());
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const ArrayOfMatrix& var1,
             const ArrayOfMatrix& var2,
             const Numeric& maxabsdiff,
             const String& error_message) {
  if (var1.nelem() != var2.nelem()) {
    std::ostringstream os;
    os << "The two arrays do not have the same size." << std::endl
       << "var1 nelem: " << var1.nelem() << std::endl
       << "var2" " nelem: " << var2.nelem() << std::endl;
    throw std::runtime_error(os.str());
  }

  bool failed = false;
  std::ostringstream fail_msg;
  for (Index i = 0; i < var1.nelem(); i++) {
    try {
      std::ostringstream vn1, vn2;
      vn1 << "var1[" << i << "]";
      vn2 << "var2" "[" << i << "]";
      Compare(var1[i],
              var2[i],
              maxabsdiff,
              error_message);
    } catch (const std::runtime_error& e) {
      failed = true;
      fail_msg << std::endl
               << e.what() << std::endl
               << "Mismatch at array index: " << i << std::endl;
    }
  }

  if (failed) throw std::runtime_error(fail_msg.str());
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const ArrayOfTensor7& var1,
             const ArrayOfTensor7& var2,
             const Numeric& maxabsdiff,
             const String& error_message) {
  if (var1.nelem() != var2.nelem()) {
    std::ostringstream os;
    os << "The two arrays do not have the same size." << std::endl
       << "var1 nelem: " << var1.nelem() << std::endl
       << "var2" " nelem: " << var2.nelem() << std::endl;
    throw std::runtime_error(os.str());
  }

  bool failed = false;
  std::ostringstream fail_msg;
  for (Index i = 0; i < var1.nelem(); i++) {
    try {
      std::ostringstream vn1, vn2;
      vn1 << "var1[" << i << "]";
      vn2 << "var2" "[" << i << "]";
      Compare(var1[i],
              var2[i],
              maxabsdiff,
              error_message);
    } catch (const std::runtime_error& e) {
      failed = true;
      fail_msg << std::endl
               << e.what() << std::endl
               << "Mismatch at array index: " << i << std::endl;
    }
  }

  if (failed) throw std::runtime_error(fail_msg.str());
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const GriddedField3& var1,
             const GriddedField3& var2,
             const Numeric& maxabsdiff,
             const String& error_message) {
  for (Index i = 0; i < var1.get_dim(); i++) {
    if (var1.get_grid_size(i) != var2.get_grid_size(i)) {
      std::ostringstream os;
      os << "var1 and var2 grid " << i
         << " do not have the same size: " << var1.get_grid_size(i)
         << " != " << var2.get_grid_size(i);
      throw std::runtime_error(os.str());
    }
    if (var1.get_grid_name(i) != var2.get_grid_name(i)) {
      std::ostringstream os;
      os << "var1 and var2 grid " << i
         << " do not have the same name: " << var1.get_grid_name(i)
         << " != " << var2.get_grid_name(i);
      throw std::runtime_error(os.str());
    }
  }

  Compare(var1.data,
          var2.data,
          maxabsdiff,
          error_message);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Sparse& var1,
             const Sparse& var2,
             const Numeric& maxabsdiff,
             const String& error_message) {
  const Index nrows = var1.nrows();
  const Index ncols = var1.ncols();

  if (var2.nrows() != nrows || var2.ncols() != ncols) {
    std::ostringstream os;
    os << "var1 (" << nrows << "," << ncols << ") and var2"
       << " (" << var2.nrows() << "," << var2.ncols()
       << ") do not have the same size.";
    throw std::runtime_error(os.str());
  }

  Numeric maxdiff = 0.0;

  for (Index r = 0; r < nrows; r++) {
    for (Index c = 0; c < ncols; c++) {
      Numeric diff = var1(r, c) - var2(r, c);

      if (std::isnan(var1(r, c)) || std::isnan(var2(r, c))) {
        if (std::isnan(var1(r, c)) && std::isnan(var2(r, c))) {
          diff = 0;
        } else if (std::isnan(var1(r, c))) {
          std::ostringstream os;
          os << "Nan found in var1, but there is no "
             << "NaN at same position in var2.\nThis "
             << "is not allowed.";
          throw std::runtime_error(os.str());
        } else {
          std::ostringstream os;
          os << "Nan found in var2, but there is no "
             << "NaN at same position in var1.\nThis "
             << "is not allowed.";
          throw std::runtime_error(os.str());
        }
      }

      if (abs(diff) > abs(maxdiff)) {
        maxdiff = diff;
      }
    }
  }

  if (abs(maxdiff) > maxabsdiff) {
    std::ostringstream os;
    os << "var1-var2 FAILED!\n";
    if (error_message.length()) os << error_message << "\n";
    os << "Max allowed deviation set to : " << maxabsdiff << std::endl
       << "but the matrices deviate with: " << maxdiff << std::endl;
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const SingleScatteringData& var1,
             const SingleScatteringData& var2,
             const Numeric& maxabsdiff,
             const String& error_message) {
  if (var1.ptype != var2.ptype) {
    std::ostringstream os;
    os << "The particle types don't match: " << std::endl
       << "var1 = " << PTypeToString(var1.ptype) << ", var2"
       << " = " << PTypeToString(var2.ptype) << std::endl;
    throw std::runtime_error(os.str());
  }
  Compare(var1.f_grid,
          var2.f_grid,
          maxabsdiff,
          error_message);
  Compare(var1.T_grid,
          var2.T_grid,
          maxabsdiff,
          error_message);
  Compare(var1.za_grid,
          var2.za_grid,
          maxabsdiff,
          error_message);
  Compare(var1.aa_grid,
          var2.aa_grid,
          maxabsdiff,
          error_message);
  Compare(var1.pha_mat_data,
          var2.pha_mat_data,
          maxabsdiff,
          error_message);
  Compare(var1.ext_mat_data,
          var2.ext_mat_data,
          maxabsdiff,
          error_message);
  Compare(var1.abs_vec_data,
          var2.abs_vec_data,
          maxabsdiff,
          error_message);
}

inline void _cr_internal_(const Numeric& var1,
                          const Numeric& var2,
                          const Numeric& maxabsreldiff,
                          const String& error_message) {
  if (var1 not_eq 0. and var2 not_eq 0.) {
    const Numeric absreldiff = abs(var1 / var2 - 1);
    if (absreldiff > maxabsreldiff) {
      std::ostringstream os;
      os << "var1-var2 FAILED!\n";
      if (error_message.length()) os << error_message << "\n";
      os << "Max allowed deviation set to: " << maxabsreldiff * 100.0 << "%"
         << std::endl
         << "but the input deviate with: " << absreldiff * 100.0 << "%\n" 
         << "If you compare non-scalar variables, the reported deviation is\n"
         << "the first one found violating the criterion. The maximum\n"
         << "difference can be higher.\n";
      throw std::runtime_error(os.str());
    }
  }
}

inline void _cr_internal_(const ConstVectorView var1,
                          const ConstVectorView var2,
                          const Numeric& maxabsreldiff,
                          const String& error_message) {
  const Index n = var1.nelem();
  if (var2.nelem() not_eq n)
    throw std::runtime_error("Cannot compare variables of different size");
  for (Index i = 0; i < n; i++)
    _cr_internal_(var1[i],
                  var2[i],
                  maxabsreldiff,
                  error_message);
}

inline void _cr_internal_(const ConstMatrixView var1,
                          const ConstMatrixView var2,
                          const Numeric& maxabsreldiff,
                          const String& error_message) {
  const Index n = var1.nrows();
  if (var2.nrows() not_eq n)
    throw std::runtime_error("Cannot compare variables of different size");
  for (Index i = 0; i < n; i++)
    _cr_internal_(var1(i, joker),
                  var2(i, joker),
                  maxabsreldiff,
                  error_message);
}

inline void _cr_internal_(const ConstTensor3View var1,
                          const ConstTensor3View var2,
                          const Numeric& maxabsreldiff,
                          const String& error_message) {
  const Index n = var1.npages();
  if (var2.npages() not_eq n)
    throw std::runtime_error("Cannot compare variables of different size");
  for (Index i = 0; i < n; i++)
    _cr_internal_(var1(i, joker, joker),
                  var2(i, joker, joker),
                  maxabsreldiff,
                  error_message);
}

inline void _cr_internal_(const ConstTensor4View var1,
                          const ConstTensor4View var2,
                          const Numeric& maxabsreldiff,
                          const String& error_message) {
  const Index n = var1.nbooks();
  if (var2.nbooks() not_eq n)
    throw std::runtime_error("Cannot compare variables of different size");
  for (Index i = 0; i < n; i++)
    _cr_internal_(var1(i, joker, joker, joker),
                  var2(i, joker, joker, joker),
                  maxabsreldiff,
                  error_message);
}

inline void _cr_internal_(const ConstTensor5View var1,
                          const ConstTensor5View var2,
                          const Numeric& maxabsreldiff,
                          const String& error_message) {
  const Index n = var1.nshelves();
  if (var2.nshelves() not_eq n)
    throw std::runtime_error("Cannot compare variables of different size");
  for (Index i = 0; i < n; i++)
    _cr_internal_(var1(i, joker, joker, joker, joker),
                  var2(i, joker, joker, joker, joker),
                  maxabsreldiff,
                  error_message);
}

inline void _cr_internal_(const ConstTensor6View var1,
                          const ConstTensor6View var2,
                          const Numeric& maxabsreldiff,
                          const String& error_message) {
  const Index n = var1.nvitrines();
  if (var2.nvitrines() not_eq n)
    throw std::runtime_error("Cannot compare variables of different size");
  for (Index i = 0; i < n; i++)
    _cr_internal_(var1(i, joker, joker, joker, joker, joker),
                  var2(i, joker, joker, joker, joker, joker),
                  maxabsreldiff,
                  error_message);
}

inline void _cr_internal_(const ConstTensor7View var1,
                          const ConstTensor7View var2,
                          const Numeric& maxabsreldiff,
                          const String& error_message) {
  const Index n = var1.nlibraries();
  if (var2.nlibraries() not_eq n)
    throw std::runtime_error("Cannot compare variables of different size");
  for (Index i = 0; i < n; i++)
    _cr_internal_(var1(i, joker, joker, joker, joker, joker, joker),
                  var2(i, joker, joker, joker, joker, joker, joker),
                  maxabsreldiff,
                  error_message);
}

template <class T>
inline void _cr_internal_(const Array<T>& var1,
                          const Array<T>& var2,
                          const Numeric& maxabsreldiff,
                          const String& error_message) {
  const Index n = var1.nelem();
  if (var2.nelem() not_eq n)
    throw std::runtime_error("Cannot compare arrays of different length");
  for (Index i = 0; i < n; i++)
    _cr_internal_(var1[i],
                  var2[i],
                  maxabsreldiff,
                  error_message);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void CompareRelative(const Numeric& var1,
                     const Numeric& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const Vector& var1,
                     const Vector& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const Matrix& var1,
                     const Matrix& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const Tensor3& var1,
                     const Tensor3& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const Tensor4& var1,
                     const Tensor4& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const Tensor5& var1,
                     const Tensor5& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const Tensor6& var1,
                     const Tensor6& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const Tensor7& var1,
                     const Tensor7& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const ArrayOfVector& var1,
                     const ArrayOfVector& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const ArrayOfMatrix& var1,
                     const ArrayOfMatrix& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const ArrayOfTensor3& var1,
                     const ArrayOfTensor3& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const ArrayOfTensor4& var1,
                     const ArrayOfTensor4& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const ArrayOfTensor5& var1,
                     const ArrayOfTensor5& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const ArrayOfTensor6& var1,
                     const ArrayOfTensor6& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const ArrayOfTensor7& var1,
                     const ArrayOfTensor7& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const ArrayOfArrayOfVector& var1,
                     const ArrayOfArrayOfVector& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const ArrayOfArrayOfMatrix& var1,
                     const ArrayOfArrayOfMatrix& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const ArrayOfArrayOfTensor3& var1,
                     const ArrayOfArrayOfTensor3& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const ArrayOfArrayOfTensor4& var1,
                     const ArrayOfArrayOfTensor4& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const ArrayOfArrayOfTensor5& var1,
                     const ArrayOfArrayOfTensor5& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const ArrayOfArrayOfTensor6& var1,
                     const ArrayOfArrayOfTensor6& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}
void CompareRelative(const ArrayOfArrayOfTensor7& var1,
                     const ArrayOfArrayOfTensor7& var2,
                     const Numeric& maxabsreldiff,
                     const String& error_message) {
  _cr_internal_(var1,
                var2,
                maxabsreldiff,
                error_message);
}

void PrintPhysicalConstants() {
  std::cout << std::setprecision(15) << std::scientific;
  std::cout << "---------------------------------------------------------\n"
            << "Numerical const in ARTS \tValue\n"
            << "Avogadro's constant:    \t " << Constant::avogadro_constant
            << '\n'
            << "Bohr's magneton:        \t " << Constant::bohr_magneton << '\n'
            << "Boltzmann's constant:   \t " << Constant::boltzmann_constant
            << '\n'
            << "Elemental charge:       \t " << Constant::elementary_charge
            << '\n'
            << "Electron mass:          \t " << Constant::electron_mass << '\n'
            << "Ideal gas constant:     \t " << Constant::ideal_gas_constant
            << '\n'
            << "Planck's constant:      \t " << Constant::planck_constant
            << '\n'
            << "Speed of light:         \t " << Constant::speed_of_light << '\n'
            << "Vacuum permittivity:    \t " << Constant::vacuum_permittivity
            << '\n'
            << "Doppler constant:       \t "
            << std::sqrt(Constant::doppler_broadening_const_squared) << '\n'
            << "---------------------------------------------------------\n";
}
