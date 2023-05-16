/**
 * @file   propagationmatrix.h
 * @author Richard Larsson
 * @date   2017-06-23
 * 
 * @brief  Stuff related to the propagation matrix.
 * 
 * This implementation takes advantage of symmetries to lower memory and computational costs
 * 
 * Some standard functions applied on the propagation matrix have been included but far from all.
 * 
 * Present support is for stokes dim 1-4 however please take note that  circular polarization delay variable's
 * position in the internal mdata structure moves depending on if stokes_dim is 3 or 4.  It is not present elsewhere...
 */

#ifndef propagationmatrix_h
#define propagationmatrix_h

#include "array.h"
#include "matpack_complex.h"
#include "matpack_data.h"

/** Class to help with hidden temporary variables for operations of type Numeric times Class
 * 
 * This class only exists in public settings and other classes are required to implement its
 * understanding in their own operators to define the level of laziness
 */
template <class base>
class LazyScale {
 public:
 
  /** Construct a new Lazy Scale object
   * 
   * @param t Constant parameter
   * @param x Scale factor
   */
  LazyScale(const base& t, const Numeric& x) : bas(t), scale(x) {}

  //! A const reference to a value
  const base& bas;

  //! A const reference to a Numeric
  const Numeric& scale;
};

template<bool matrix>
constexpr Index need2stokes(Index nstokes_needed) {
  if constexpr (matrix) {
    if (nstokes_needed == 7) return 4;
    if (nstokes_needed == 4) return 3;
    if (nstokes_needed == 2) return 2;
    if (nstokes_needed == 1) return 1;
  }

  if (nstokes_needed < 5 and nstokes_needed > 0) return nstokes_needed;

  ARTS_ASSERT (false, "Cannot understand the input Stokes dimensions");
  return std::numeric_limits<Index>::max();
}

template<bool matrix>
constexpr Index stokes2need(Index nstokes) {
  if constexpr (matrix) {
      if (nstokes == 4) return 7;
      if (nstokes == 3) return 4;
      if (nstokes == 2) return 2;
      if (nstokes == 1) return 1;
  }
  
  if (nstokes < 5 and nstokes > 0) return nstokes;
  
  ARTS_ASSERT (false, "Cannot understand the input Stokes dimensions");
  return std::numeric_limits<Index>::max();
}

/*! Propagation Matrix Holder Class With Some Computational Capabilities
 * 
 * The idea comes from the fact that the propagation matrix has the looks
 * 
 * /            \
 * | a  b  c  d |
 * | b  a  u  v |
 * | c -u  a  w |
 * | d -v -w  a |
 * \            /
 * 
 * So we instead store the inner parts of the matrix linearly as [a b c d u v w],
 * and all computations happens on these variables instead.  Note that the variables
 * changes with the Stokes dimension
 * 
 * Stokes Dim 4: [a b c d u v w]
 * Stokes Dim 3: [a b c u]
 * Stokes Dim 2: [a b]
 * Stokes Dim 1: [a]
 * 
 * And for devs: the u-variable is at vector position of the Stokes Dim in case 4 and 3, 
 * and that the other variables never change their positions, which is why the switch
 * cases are all consistently fall-through-able in this class
 */
class PropagationMatrix {
 public:
  /** Initialize variable sizes
   * 
   * Will create a Tensor4 of the size and order (nr_aa, nr_za, 
   * NumberOfNeededVectors(), nr_frequencies)
   * 
   * @param[in] nr_frequencies Number of Dirac frequencies
   * @param[in] stokes_dim Stokes dimensionality
   * @param[in] nr_za Number of Dirac Zeniths
   * @param[in] nr_aa Number of Dirac Azimuths
   * @param[in] v Initial values of things in the created Tensor4
   */
  PropagationMatrix(const Index nr_frequencies = 0,
                    const Index stokes_dim = 1,
                    const Index nr_za = 1,
                    const Index nr_aa = 1,
                    const Numeric v = 0.0)
      : mfreqs(nr_frequencies),
        mstokes_dim(stokes_dim),
        mza(nr_za),
        maa(nr_aa) {
    ARTS_ASSERT(mstokes_dim < 5 and mstokes_dim > 0);
    mdata = Tensor4(maa, mza, mfreqs, NumberOfNeededVectors(), v);
  }

  /** Construct a new Propagation Matrix object
   * 
   * @param[in] pm Old propagation matrix object to copy
   */
  PropagationMatrix(const PropagationMatrix& pm) = default;

  /** Construct a new Propagation Matrix object
   * 
   * @param[in] pm old Propagation Matrix object to move from
   */
  PropagationMatrix(PropagationMatrix&& pm) noexcept
      : mfreqs(pm.mfreqs),
        mstokes_dim(pm.mstokes_dim),
        mza(pm.mza),
        maa(pm.maa),
        mdata(std::move(pm.mdata)),
        mvectortype(pm.mvectortype) {}

  /** Construct a new Propagation Matrix object
   * 
   * @param[in] x Tensor4 object to use to initialize from
   */
  explicit PropagationMatrix(Tensor4 x)
      : mfreqs(x.nrows()),
        mza(x.npages()),
        maa(x.nbooks()),
        mdata(std::move(x)) {
    switch (x.ncols()) {
      case 7:
        mstokes_dim = 4;
        break;
      case 4:
        mstokes_dim = 3;
        break;
      case 2:
        mstokes_dim = 2;
        break;
      case 1:
        mstokes_dim = 1;
        break;
      default:
        ARTS_ASSERT(false, "Tensor4 not representative of PropagationMatrix");
    }
  }

  /** Initialize from single stokes_dim-by-stokes_dim matrix
   * 
   * @param[in] x The matrix
   * @param[in] assume_fit Assume a correct fit?  Do not set this in manual interface
   */
  explicit PropagationMatrix(const ConstMatrixView& x, const bool& assume_fit = false)
      : mfreqs(1), mstokes_dim(x.ncols()), mza(1), maa(1) {
    ARTS_ASSERT(mstokes_dim < 5 and mstokes_dim > 0);
    mvectortype = false;

    if (not assume_fit) {
      ARTS_ASSERT(FittingShape(x),
        "Matrix not fit as propagation matrix");
    }

    mdata.resize(1, 1, 1, NumberOfNeededVectors());

    switch (mstokes_dim) {
      case 4:
        mdata(0, 0, 0, 5) = x(1, 3);
        mdata(0, 0, 0, 6) = x(2, 3);
        mdata(0, 0, 0, 3) = x(0, 3); /* FALLTHROUGH */
      case 3:
        mdata(0, 0, 0, mstokes_dim) = x(1, 2);
        mdata(0, 0, 0, 2) = x(0, 2); /* FALLTHROUGH */
      case 2:
        mdata(0, 0, 0, 1) = x(0, 1); /* FALLTHROUGH */
      case 1:
        mdata(0, 0, 0, 0) = x(0, 0); /* FALLTHROUGH */
    }
  };

  /** The stokes dimension of the propagation matrix */
  [[nodiscard]] Index StokesDimensions() const { return mstokes_dim; };

  /** The number of frequencies of the propagation matrix */
  [[nodiscard]] Index NumberOfFrequencies() const { return mfreqs; };

  /** The number of zenith angles of the propagation matrix */
  [[nodiscard]] Index NumberOfZenithAngles() const { return mza; };

  /** The number of azimuth angles of the propagation matrix */
  [[nodiscard]] Index NumberOfAzimuthAngles() const { return maa; };
  
  [[nodiscard]] bool OK() const {return mdata.ncols() == NumberOfNeededVectors() and mdata.nrows() == mfreqs and mdata.npages() == mza and mdata.nbooks() == maa;}

  /** Asks if the class is empty */
  [[nodiscard]] bool IsEmpty() const { return not mfreqs or not mza or not maa; };

  /** False if any non-zeroes in internal Matrix representation
   * 
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return False if any non-zeroes
   */
  [[nodiscard]] bool IsZero(const Index iv = 0,
                            const Index iz = 0,
                            const Index ia = 0) const {
    // FIXME: matpack does not do pointers in a clear manner
    // return std::any_of(mdata(ia, iz, iv, joker).begin(), mdata(ia, iz, iv, joker).end(), [](auto& x){return x not_eq 0;});
    for (auto& n : mdata(ia, iz, iv, joker))
      if (n not_eq 0) return false;
    return true;
  };

  /** False if diagonal element is non-zero in internal Matrix representation 
   * 
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return False if diagonal is non-zero
   */
  [[nodiscard]] bool IsRotational(const Index iv = 0,
                                  const Index iz = 0,
                                  const Index ia = 0) const {
    if (mdata(ia, iz, iv, 0) == 0.0)
      return true;
    return false;
  };

  /** The number of required vectors to fill this PropagationMatrix */
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wreturn-type"
  [[nodiscard]] Index NumberOfNeededVectors() const {
    if (not mvectortype) {
      return stokes2need<true>(mstokes_dim);
    }
   return stokes2need<false>(mstokes_dim);
  }
  #pragma GCC diagnostic pop

  /** access operator.  Please refrain from using this if possible since it copies */
  Numeric operator()(const Index iv = 0,
                     const Index is1 = 0,
                     const Index is2 = 0,
                     const Index iz = 0,
                     const Index ia = 0) const;

  /** Adds the Faraday rotation to the PropagationMatrix at required position
   * 
   * No vector function exists since rot is a function of frequency
   * 
   * @param[in] rot rotation
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void AddFaraday(const Numeric& rot,
                  const Index iv = 0,
                  const Index iz = 0,
                  const Index ia = 0) {
    mdata(ia, iz, iv, mstokes_dim) += rot;
  }

  /** Sets the Faraday rotation to the PropagationMatrix at required position
   * 
   * No vector function exists since rot is a function of frequency
   * 
   * @param[in] rot rotation
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void SetFaraday(const Numeric& rot,
                  const Index iv = 0,
                  const Index iz = 0,
                  const Index ia = 0) {
    mdata(ia, iz, iv, mstokes_dim) = rot;
  }

  /** Sets the dense matrix.  Avoid using if possible.
   * 
   * @param[in,out] Full propagation matrix
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void MatrixAtPosition(MatrixView ret,
                        const Index iv = 0,
                        const Index iz = 0,
                        const Index ia = 0) const;

  /** Move operator
   * 
   * @param[in] pm PropagationMatrix to move from
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator=(PropagationMatrix&& pm) noexcept {
    if (this != &pm) {
      mfreqs = pm.mfreqs;
      mstokes_dim = pm.mstokes_dim;
      mza = pm.mza;
      maa = pm.maa;
      mdata = std::move(pm.mdata);
      mvectortype = pm.mvectortype;
    }
    return *this;
  }

  /** Laze equal to opeartor
   * 
   * @param[in] lpms lazified propagation matrix
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator=(const LazyScale<PropagationMatrix>& lpms) {
    operator=(lpms.bas);
    mdata *= lpms.scale;
    return *this;
  }

  /** Copy operator
   * 
   * @param[in] other PropagationMatrix to copy
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator=(const PropagationMatrix& other) = default;

  /** Sets all data to constant
   * 
   * @param[in] x new constant value
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator=(const Numeric& x) {
    mdata = x;
    return *this;
  }

  /** Set the At Position object
   * 
   * @param[in] x Object to set
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void SetAtPosition(const PropagationMatrix& x,
                     const Index iv = 0,
                     const Index iz = 0,
                     const Index ia = 0) {
    mdata(ia, iz, iv, joker) = x.mdata(ia, iz, iv, joker);
  }

  /** Set the At Position object
   * 
   * @param[in] x Object to set
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void SetAtPosition(const ConstMatrixView& x,
                     const Index iv = 0,
                     const Index iz = 0,
                     const Index ia = 0);

  
  /** Set the At Position object
   * 
   * @param[in] x Object to set
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void SetAtPosition(const Numeric& x,
                     const Index iv = 0,
                     const Index iz = 0,
                     const Index ia = 0) {
    mdata(ia, iz, iv, joker) = x;
  }

  /** Divide operator
   * 
   * @param[in] other Divide by other
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator/=(const PropagationMatrix& other) {
    mdata /= other.mdata;
    return *this;
  }

  /** Divide operator
   * 
   * @param[in] other Divide by other
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator/=(const ConstVectorView& x) {
    for (Index i = 0; i < NumberOfNeededVectors(); i++) {
      for (Index j = 0; j < mza; j++) {
        for (Index k = 0; k < maa; k++) {
          mdata(k, j, joker, i) /= x;
        }
      }
    }
    return *this;
  }

  /** Divide operator
   * 
   * @param[in] other Divide by other
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator/=(const Numeric& x) {
    mdata /= x;
    return *this;
  }

  /** Divide at position
   * 
   * @param[in] other Divide by other
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void DivideAtPosition(const PropagationMatrix& x,
                        const Index iv = 0,
                        const Index iz = 0,
                        const Index ia = 0) {
    mdata(ia, iz, iv, joker) /= x.mdata(ia, iz, iv, joker);
  }

  /** Divide at position
   * 
   * @param[in] other Divide by other
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void DivideAtPosition(const ConstMatrixView& x,
                        const Index iv = 0,
                        const Index iz = 0,
                        const Index ia = 0);

  /** Divide at position
   * 
   * @param[in] other Divide by other
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void DivideAtPosition(const Numeric& x,
                        const Index iv = 0,
                        const Index iz = 0,
                        const Index ia = 0) {
    mdata(ia, iz, iv, joker) /= x;
  }

  /** Multiply operator
   * 
   * @param[in] other Multiply by other
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator*=(const PropagationMatrix& other) {
    mdata *= other.mdata;
    return *this;
  }

  /** Multiply operator
   * 
   * @param[in] other Multiply by other
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator*=(const ConstVectorView& x) {
    for (Index i = 0; i < NumberOfNeededVectors(); i++)
      for (Index j = 0; j < mza; j++)
        for (Index k = 0; k < maa; k++) mdata(k, j, joker, i) *= x;
    return *this;
  }

  /** Multiply operator
   * 
   * @param[in] other Multiply by other
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator*=(const Numeric& x) {
    mdata *= x;
    return *this;
  }

  /** Multiply operator at position
   * 
   * @param[in] other Multiply by other
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void MultiplyAtPosition(const PropagationMatrix& x,
                          const Index iv = 0,
                          const Index iz = 0,
                          const Index ia = 0) {
    mdata(ia, iz, iv, joker) *= x.mdata(ia, iz, iv, joker);
  }

  /** Multiply operator at position
   * 
   * @param[in] other Multiply by other
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void MultiplyAtPosition(const ConstMatrixView& x,
                          const Index iv = 0,
                          const Index iz = 0,
                          const Index ia = 0);

  /** Multiply operator at position
   * 
   * @param[in] other Multiply by other
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void MultiplyAtPosition(const Numeric& x,
                          const Index iv = 0,
                          const Index iz = 0,
                          const Index ia = 0) {
    mdata(ia, iz, iv, joker) *= x;
  }

  /** Addition operator
   * 
   * @param[in] other Addition by other
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator+=(const PropagationMatrix& other) {
    mdata += other.mdata;
    return *this;
  }

  /** Addition operator
   * 
   * @param[in] other Addition by other
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator+=(const LazyScale<PropagationMatrix>& lpms) {
    MultiplyAndAdd(lpms.scale, lpms.bas);
    return *this;
  }

  /** Addition operator
   * 
   * @param[in] other Addition by other
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator+=(const ConstVectorView& x) {
    for (Index i = 0; i < NumberOfNeededVectors(); i++)
      for (Index j = 0; j < mza; j++)
        for (Index k = 0; k < maa; k++) mdata(k, j, joker, i) += x;
    return *this;
  }

  /** Addition operator
   * 
   * @param[in] other Addition by other
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator+=(const Numeric& x) {
    mdata += x;
    return *this;
  }

  /** Addition at position operator
   * 
   * @param[in] other Addition by other
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void AddAtPosition(const PropagationMatrix& x,
                     const Index iv = 0,
                     const Index iz = 0,
                     const Index ia = 0) {
    mdata(ia, iz, iv, joker) += x.mdata(ia, iz, iv, joker);
  }

  /** Addition at position operator
   * 
   * @param[in] other Addition by other
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void AddAtPosition(const ConstMatrixView& x,
                     const Index iv = 0,
                     const Index iz = 0,
                     const Index ia = 0);

  /** Addition at position operator
   * 
   * @param[in] other Addition by other
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void AddAtPosition(const Numeric& x,
                     const Index iv = 0,
                     const Index iz = 0,
                     const Index ia = 0) {
    mdata(ia, iz, iv, joker) += x;
  }

  /** Subtraction operator
   * 
   * @param[in] other Subtract by other
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator-=(const PropagationMatrix& other) {
    mdata -= other.mdata;
    return *this;
  }

  /** Subtraction operator
   * 
   * @param[in] other Subtract by other
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator-=(const ConstVectorView& x) {
    for (Index i = 0; i < NumberOfNeededVectors(); i++) {
      for (Index j = 0; j < mza; j++) {
        for (Index k = 0; k < maa; k++) {
          mdata(k, j, joker, i) -= x;
        }
      }
    }
    return *this;
  }

  /** Subtraction operator
   * 
   * @param[in] other Subtract by other
   * @return PropagationMatrix& *this
   */
  PropagationMatrix& operator-=(const Numeric& x) {
    mdata -= x;
    return *this;
  }

  /** Subtraction at position
   * 
   * @param[in] other Subtract by other
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void RemoveAtPosition(const PropagationMatrix& x,
                        const Index iv = 0,
                        const Index iz = 0,
                        const Index ia = 0) {
    mdata(ia, iz, iv, joker) -= x.mdata(ia, iz, iv, joker);
  }

  /** Subtraction at position
   * 
   * @param[in] other Subtract by other
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void RemoveAtPosition(const ConstMatrixView& x,
                        const Index iv = 0,
                        const Index iz = 0,
                        const Index ia = 0);

  /** Subtraction at position
   * 
   * @param[in] other Subtract by other
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void RemoveAtPosition(const Numeric& x,
                        const Index iv = 0,
                        const Index iz = 0,
                        const Index ia = 0) {
    mdata(ia, iz, iv, joker) -= x;
  }

  /** Adds as a Stokes vector at position
   * 
   * @param[in] x 
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void AddAbsorptionVectorAtPosition(const ConstVectorView& x,
                                     const Index iv = 0,
                                     const Index iz = 0,
                                     const Index ia = 0) {
    for (Index i = 0; i < mstokes_dim; i++) mdata(ia, iz, iv, i) += x[i];
  }

  /** Add the average of the two input at position
   * 
   * @param[in] mat1 input 1
   * @param[in] mat2 input 2
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void AddAverageAtPosition(const ConstMatrixView& mat1,
                            const ConstMatrixView& mat2,
                            const Index iv = 0,
                            const Index iz = 0,
                            const Index ia = 0);

  /** Multiply input by scalar and add to this
   * 
   * @param x scalar
   * @param y input
   */
  void MultiplyAndAdd(const Numeric x, const PropagationMatrix& y);

  /** Return the matrix inverse at the position
   * 
   * @param[in, out] ret Inversed matrix
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void MatrixInverseAtPosition(MatrixView ret,
                               const Index iv = 0,
                               const Index iz = 0,
                               const Index ia = 0) const;

  /** Tests of the input matrix fits Propagation Matrix style
   * 
   * @param[in] x Input matrix
   */
  [[nodiscard]] bool FittingShape(const ConstMatrixView& x) const;

  /** Get a Tensor3 object from this
   * 
   * @param[in,out] tensor3 New tensor 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void GetTensor3(Tensor3View tensor3, const Index iz = 0, const Index ia = 0);

  /** Vector view to diagonal elements
   * 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView Diagonal elements
   */
  VectorView Kjj(const Index iz = 0, const Index ia = 0) {
    return mdata(ia, iz, joker, 0);
  }

  /** Vector view to K(0, 1) elements
   * 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView K(0, 1) elements
   */
  VectorView K12(const Index iz = 0, const Index ia = 0) {
    return mdata(ia, iz, joker, 1);
  }

  /** Vector view to K(0, 2) elements
   * 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView K(0, 2) elements
   */
  VectorView K13(const Index iz = 0, const Index ia = 0) {
    return mdata(ia, iz, joker, 2);
  }

  /** Vector view to K(0, 3) elements
   * 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView K(0, 3) elements
   */
  VectorView K14(const Index iz = 0, const Index ia = 0) {
    return mdata(ia, iz, joker, 3);
  }

  /** Vector view to K(1, 2) elements
   * 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView K(1, 2) elements
   */
  VectorView K23(const Index iz = 0, const Index ia = 0) {
    return mdata(ia, iz, joker, mstokes_dim);
  }

  /** Vector view to K(1, 3) elements
   * 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView K(1, 3) elements
   */
  VectorView K24(const Index iz = 0, const Index ia = 0) {
    return mdata(ia, iz, joker, 5);
  }

  /** Vector view to K(2, 3) elements
   * 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView K(2, 3) elements
   */
  VectorView K34(const Index iz = 0, const Index ia = 0) {
    return mdata(ia, iz, joker, 6);
  }

  /** Vector view to diagonal elements
   * 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView Diagonal elements
   */
  [[nodiscard]] ConstVectorView Kjj(const Index iz = 0, const Index ia = 0) const {
    return mdata(ia, iz, joker, 0);
  }

  /** Vector view to K(0, 1) elements
   * 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView K(0, 1) elements
   */
  [[nodiscard]] ConstVectorView K12(const Index iz = 0, const Index ia = 0) const {
    return mdata(ia, iz, joker, 1);
  }

  /** Vector view to K(0, 2) elements
   * 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView K(0, 2) elements
   */
  [[nodiscard]] ConstVectorView K13(const Index iz = 0, const Index ia = 0) const {
    return mdata(ia, iz, joker, 2);
  }

  /** Vector view to K(0, 3) elements
   * 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView K(0, 3) elements
   */
  [[nodiscard]] ConstVectorView K14(const Index iz = 0, const Index ia = 0) const {
    return mdata(ia, iz, joker, 3);
  }

  /** Vector view to K(1, 3) elements
   * 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView K(1, 3) elements
   */
  [[nodiscard]] ConstVectorView K23(const Index iz = 0, const Index ia = 0) const {
    return mdata(ia, iz, joker, mstokes_dim);
  }

  /** Vector view to K(1, 3) elements
   * 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView K(1, 3) elements
   */
  [[nodiscard]] ConstVectorView K24(const Index iz = 0, const Index ia = 0) const {
    return mdata(ia, iz, joker, 5);
  }

  /** Vector view to diagonal elements
   * 
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView Diagonal elements
   */
  [[nodiscard]] ConstVectorView K34(const Index iz = 0, const Index ia = 0) const {
    return mdata(ia, iz, joker, 6);
  }

  /** Sets all data to zero */
  void SetZero() { mdata = 0.0; }

  /** Get full view to data */
  Tensor4& Data() { return mdata; }

  /** Get full const view to data */
  [[nodiscard]] const Tensor4& Data() const { return mdata; }

  /** Multiply the matrix input from the left of this at position
   * 
   * @param[in,out] out Outgoing matrix
   * @param[in] in Incoming matrix
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void LeftMultiplyAtPosition(MatrixView out,
                              const ConstMatrixView& in,
                              const Index iv = 0,
                              const Index iz = 0,
                              const Index ia = 0) const;


  /** Multiply the matrix input from the right of this at position
   * 
   * @param[in,out] out Outgoing matrix
   * @param[in] in Incoming matrix
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void RightMultiplyAtPosition(MatrixView out,
                               const ConstMatrixView& in,
                               const Index iv = 0,
                               const Index iz = 0,
                               const Index ia = 0) const;

  friend std::ostream& operator<<(std::ostream& os, const PropagationMatrix& pm);

  void swap(PropagationMatrix&) noexcept;
  friend void swap(PropagationMatrix&, PropagationMatrix&) noexcept;

 protected:
  Index mfreqs, mstokes_dim;
  Index mza, maa;
  Tensor4 mdata;
  bool mvectortype{false};
};

using ArrayOfPropagationMatrix = Array<PropagationMatrix>;
using ArrayOfArrayOfPropagationMatrix = Array<ArrayOfPropagationMatrix>;

/** Compute the matrix exponent as the transmission matrix of this propagation matrix
 * 
 * The propagation matrix is multiplied by -r and level-averaged before exponent is applied.
 * 
 * upper_level and lower_level propagation matrices should thus be the level matrices and r the distance between these
 * levels.  The same is true for the derivative matrices.
 * 
 * Stokes dim 1 and 4 have been tested more.  Stokes dim 2 and 3 have been found to work but could still have hidden errors for uncommon cases
 * 
 * @param[in,out] T transmission tensor with outmost dimension being frequency
 * @param[in] r the distance over which the propagation matrix causes the transmission
 * @param[in] upper_level The upper level propagation matrix
 * @param[in] lower_level The lower level propagation matrix
 * @param[in] iz Zenith index
 * @param[in] ia Azimuth index
 */
void compute_transmission_matrix(Tensor3View T,
                                 const Numeric& r,
                                 const PropagationMatrix& upper_level,
                                 const PropagationMatrix& lower_level,
                                 const Index iz = 0,
                                 const Index ia = 0);


/** Compute the matrix exponent as the transmission matrix of this propagation matrix
 * 
 * The propagation matrix is multiplied by -r and level-averaged before exponent is applied.
 * 
 * Stokes dim 1 and 4 have been tested more.  Stokes dim 2 and 3 have been found to work but could still have hidden errors for uncommon cases
 * 
 * @param[in,out] T transmission matrix at the position
 * @param[in] r the distance over which the propagation matrix causes the transmission
 * @param[in] averaged_propagation_matrix The propagation matrix
 * @param[in] iv Frequency index
 * @param[in] iz Zenith index
 * @param[in] ia Azimuth index
 */
void compute_transmission_matrix_from_averaged_matrix_at_frequency(
    MatrixView T,
    const Numeric& r,
    const PropagationMatrix& averaged_propagation_matrix,
    const Index iv,
    const Index iz = 0,
    const Index ia = 0);

/*! Compute the matrix exponent as the transmission matrix of this propagation matrix
 * 
 * The propagation matrix is multiplied by -r and level-averaged before exponent is applied.
 * 
 * upper_level and lower_level propagation matrices should thus be the level matrices and r the distance between these
 * levels.  The same is true for the derivative matrices.
 * 
 * Stokes dim 1 and 4 have been tested more.  Stokes dim 2 and 3 have been found to work but could still have hidden errors for uncommon cases
 * 
 * @param[in,out] T transmission tensor with outmost dimension being frequency
 * @param[in,out] dT_upp transmission tensors derivative with respect to derivatives of the propagation matrix for upper level
 * @param[in,out] dT_low transmission tensors derivative with respect to derivatives of the propagation matrix for lower level
 * @param[in] r the distance over which the propagation matrix causes the transmission
 * @param[in] upper_level The upper level propagation matrix
 * @param[in] lower_level The lower level propagation matrix
 * @param[in] dprop_mat_upp derivatives of the upper propagation matrix with respect to some parameter (is multiplied by -0.5 r)
 * @param[in] dprop_mat_low derivatives of the lower propagation matrix with respect to some parameter (is multiplied by -0.5 r)
 * @param[in] dr_dTu Distance temperature derivative for upper level
 * @param[in] dr_dTl Distance temperature derivative for lower level
 * @param[in] it Position of temperature derivative (ignored at -1)
 * @param[in] iz Zenith index
 * @param[in] ia Azimuth index
 */
void compute_transmission_matrix_and_derivative(
    Tensor3View T,
    Tensor4View dT_upper_level,
    Tensor4View dT_lower_level,
    const Numeric& r,
    const PropagationMatrix& upper_level,
    const PropagationMatrix& lower_level,
    const Array<PropagationMatrix>& dprop_mat_upper_level,
    const Array<PropagationMatrix>& dprop_mat_lower_level,
    const Numeric& dr_dTu = 0.0,
    const Numeric& dr_dTl = 0.0,
    const Index it = -1,
    const Index iz = 0,
    const Index ia = 0);

/** Stokes vector is as Propagation matrix but only has 4 possible values */
class StokesVector final : public PropagationMatrix {
 public:

  /** Initialize variable sizes
   * 
   * Will create a Tensor4 of the size and order (nr_aa, nr_za, 
   * NumberOfNeededVectors(), nr_frequencies)
   * 
   * \param[in] nr_frequencies Number of Dirac frequencies
   * \param[in] stokes_dim Stokes dimensionality
   * \param[in] nr_za Number of Dirac Zeniths
   * \param[in] nr_aa Number of Dirac Azimuths
   * \param[in] v Initial values of things in the created Tensor4
   */
  StokesVector(const Index nr_frequencies = 0,
               const Index stokes_dim = 1,
               const Index nr_za = 1,
               const Index nr_aa = 1,
               const Numeric& v = 0.0) {
    mvectortype = true;
    mfreqs = nr_frequencies;
    mstokes_dim = stokes_dim;
    mza = nr_za;
    maa = nr_aa;
    ARTS_ASSERT(mstokes_dim < 5 and mstokes_dim > 0);
    mdata = Tensor4(maa, mza, mfreqs, mstokes_dim, v);
  };

  /** Construct a new Propagation Matrix object
   * 
   * @param[in] x Tensor4 object to use to initialize from
   */
  explicit StokesVector(Tensor4 x) {
    mfreqs = x.nrows();
    mstokes_dim = x.ncols();
    mza = x.npages();
    maa = x.nbooks();
    mdata = std::move(x);
    mvectortype = true;
    ARTS_ASSERT (not (mstokes_dim > 4 or mstokes_dim < 1),
      "Tensor4 is bad for StokesVector");
  }

  /** Construct a new Stokes Vector object
   * 
   * @param x Single Stokes vector
   */
  explicit StokesVector(const ConstVectorView& x) {
    mfreqs = 1;
    mstokes_dim = x.nelem();
    mza = 1;
    maa = 1;
    ARTS_ASSERT(mstokes_dim < 5 and mstokes_dim > 0);
    mvectortype = true;
    mdata.resize(1, 1, 1, mstokes_dim);
    for (Index i = 0; i < mstokes_dim; i++) mdata(0, 0, 0, i) = x[i];
  };

  /** Construct a new Stokes Vector object
   *
   * @param x Spectrum of Stokes vector
   */
  explicit StokesVector(ConstMatrixView x) {
    mfreqs = x.nrows();
    mstokes_dim = x.ncols();
    mza = 1;
    maa = 1;
    ARTS_ASSERT(mstokes_dim < 5 and mstokes_dim > 0);
    mvectortype = true;
    mdata.resize(1, 1, mfreqs, mstokes_dim);
    for (Index j = 0; j < mstokes_dim; j++) {
      for (Index i = 0; i < mstokes_dim; i++) {
        mdata(0, 0, j, i) = x(j, i);
      }
    }
  };

  /** Construct a new Stokes Vector as a scale between two others
   * 
   * @param a input 1
   * @param b input 2
   * @param scale Scale of inputs
   */
  StokesVector(const StokesVector& a,
               const StokesVector& b,
               const Numeric& scale = 0.5) {
    mfreqs = a.NumberOfFrequencies();
    mstokes_dim = a.StokesDimensions();
    mza = a.NumberOfZenithAngles();
    maa = a.NumberOfAzimuthAngles();
    mdata.resize(maa, mza, mfreqs, mstokes_dim);
    mvectortype = true;

    for (Index i = 0; i < maa; i++)
      for (Index j = 0; j < mza; j++)
        for (Index k = 0; k < mfreqs; k++) {
          switch (mstokes_dim) {
            case 4:
              mdata(i, j, k, 3) = (a.mdata(i, j, k, 3) + b.mdata(i, j, k, 3)) *
                                  scale; /* FALLTHROUGH */
            case 3:
              mdata(i, j, k, 2) = (a.mdata(i, j, k, 2) + b.mdata(i, j, k, 2)) *
                                  scale; /* FALLTHROUGH */
            case 2:
              mdata(i, j, k, 1) = (a.mdata(i, j, k, 1) + b.mdata(i, j, k, 1)) *
                                  scale; /* FALLTHROUGH */
            case 1:
              mdata(i, j, k, 0) = (a.mdata(i, j, k, 0) + b.mdata(i, j, k, 0)) *
                                  scale; /* FALLTHROUGH */
          }
        }
  };

  /** The number of required vectors to fill this StokesVector */
  [[nodiscard]] Index NumberOfNeededVectors() const { return mstokes_dim; }

  /** Addition operator
   * 
   * @param x Stokes Vector or Propagation Matrix to add
   * @return StokesVector& *this
   */
  StokesVector& operator+=(const PropagationMatrix& x) {
    mdata += x.Data()(joker, joker, joker, Range(0, mstokes_dim, 1));
    return *this;
  }

  /** Addition operator
   * 
   * @param x Lazy addition
   * @return StokesVector& *this
   */
  StokesVector& operator+=(const LazyScale<PropagationMatrix>& lpms) {
    MultiplyAndAdd(lpms.scale, lpms.bas);
    return *this;
  }

  /** Set *this from a Propagation matrix
   * 
   * @param x Propagation matrix
   * @return StokesVector& *this
   */
  StokesVector& operator=(const PropagationMatrix& x) {
    mstokes_dim = x.StokesDimensions();
    mfreqs = x.NumberOfFrequencies();
    mza = x.NumberOfZenithAngles();
    maa = x.NumberOfAzimuthAngles();
    mdata.resize(maa, mza, mfreqs, mstokes_dim);
    mdata(joker, joker, joker, Range(0, mstokes_dim, 1)) = x.Data()(joker, joker, joker, Range(0, mstokes_dim, 1));
    return *this;
  }

  /** Set this lazily
   * 
   * @param[in] lpms Lazy propagation matrix
   * @return StokesVector& *this
   */
  StokesVector& operator=(const LazyScale<PropagationMatrix>& lpms) {
    operator=(lpms.bas);
    mdata *= lpms.scale;
    return *this;
  }

  /** Set this to constant value
   * 
   * @param[in] x constant value
   * @return StokesVector& *this
   */
  StokesVector& operator=(const Numeric& x) {
    mdata = x;
    return *this;
  }

  /** Add a scaled version of the input
   * 
   * @param[in] x Scale
   * @param[in] y Input
   */
  void MultiplyAndAdd(const Numeric x, const PropagationMatrix& y) {
    ARTS_ASSERT(mstokes_dim == y.StokesDimensions());
    ARTS_ASSERT(mfreqs == y.NumberOfFrequencies());
    ARTS_ASSERT(mza == y.NumberOfZenithAngles());
    ARTS_ASSERT(maa == y.NumberOfAzimuthAngles());

    const Tensor4& data = y.Data();

    for (Index i = 0; i < maa; i++) {
      for (Index j = 0; j < mza; j++) {
        for (Index k = 0; k < mfreqs; k++) {
          switch (mstokes_dim) {
            case 4:
              mdata(i, j, k, 3) += x * data(i, j, k, 3); /* FALLTHROUGH */
            case 3:
              mdata(i, j, k, 2) += x * data(i, j, k, 2); /* FALLTHROUGH */
            case 2:
              mdata(i, j, k, 1) += x * data(i, j, k, 1); /* FALLTHROUGH */
            case 1:
              mdata(i, j, k, 0) += x * data(i, j, k, 0); /* FALLTHROUGH */
          }
        }
      }
    }
  }

  /** Get a vectorview to the position
   * 
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView To StokesVector
   */
  auto VectorAtPosition(const Index iv = 0,
                        const Index iz = 0,
                        const Index ia = 0) {
    return mdata(ia, iz, iv, joker);
  }

  /** Get a vectorview to the position
   * 
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return VectorView To StokesVector
   */
  [[nodiscard]] auto VectorAtPosition(const Index iv = 0,
                                      const Index iz = 0,
                                      const Index ia = 0) const {
    return mdata(ia, iz, iv, joker);
  }

  void SetAtPosition(const ConstVectorView& x,
                     const Index iv = 0,
                     const Index iz = 0,
                     const Index ia = 0) {
    mdata(ia, iz, iv, joker) = x;
  }

  /** Add the average of both inputs at position
   * 
   * @param[in] vec1 input 1
   * @param[in] vec2 input 2
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   */
  void AddAverageAtPosition(const ConstVectorView& vec1,
                            const ConstVectorView& vec2,
                            const Index iv = 0,
                            const Index iz = 0,
                            const Index ia = 0) {
    switch (mstokes_dim) {
      case 4:
        mdata(ia, iz, iv, 3) += (vec1[3] + vec2[3]) * 0.5; /* FALLTHROUGH */
      case 3:
        mdata(ia, iz, iv, 2) += (vec1[2] + vec2[2]) * 0.5; /* FALLTHROUGH */
      case 2:
        mdata(ia, iz, iv, 1) += (vec1[1] + vec2[1]) * 0.5; /* FALLTHROUGH */
      case 1:
        mdata(ia, iz, iv, 0) += (vec1[0] + vec2[0]) * 0.5; /* FALLTHROUGH */
    }
  }

  /** Checks if vector is polarized
   * 
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return true if polarized
   * @return false if not polarized
   */
  [[nodiscard]] bool IsPolarized(const Index iv = 0,
                                 const Index iz = 0,
                                 const Index ia = 0) const {
    switch (mstokes_dim) {
      case 4:
        if (K14(iz, ia)[iv] not_eq 0.0) return true; /* FALLTHROUGH */
      case 3:
        if (K13(iz, ia)[iv] not_eq 0.0) return true; /* FALLTHROUGH */
      case 2:
        if (K12(iz, ia)[iv] not_eq 0.0) return true; /* FALLTHROUGH */
    }
    return false;
  }

  /** Checks if vector is polarized
   * 
   * @param[in] iv Frequency index
   * @param[in] iz Zenith index
   * @param[in] ia Azimuth index
   * @return true if not polarized
   * @return false if polarized
   */
  [[nodiscard]] bool IsUnpolarized(const Index iv = 0,
                                   const Index iz = 0,
                                   const Index ia = 0) const {
    return not IsPolarized(iv, iz, ia);
  }
  
  [[nodiscard]] bool allZeroes() const {
    for (Index i=0; i<maa; i++)
      for (Index j=0; j<mza; j++)
        for (Index k=0; k<mfreqs; k++)
          for (Index m=0; m<mstokes_dim; m++)
            if (mdata(i, j, k, m) not_eq 0)
              return false;
    return true;
  }

  friend std::ostream& operator<<(std::ostream& os, const StokesVector& pm);
};

using ArrayOfStokesVector = Array<StokesVector>;
using ArrayOfArrayOfStokesVector = Array<ArrayOfStokesVector>;
using ArrayOfArrayOfArrayOfStokesVector = Array<ArrayOfArrayOfStokesVector>;

/** Returns a lazy multiplier
 * 
 * @param[in] pm Propagation matrix
 * @param[in] x Scale
 * @return LazyScale<PropagationMatrix> A lazy multiplier
 */
LazyScale<PropagationMatrix> operator*(const PropagationMatrix& pm,
                                       const Numeric& x);

/** Returns a lazy multiplier
 * 
 * @param[in] x Scale
 * @param[in] pm Propagation matrix
 * @return LazyScale<PropagationMatrix> A lazy multiplier
 */
LazyScale<PropagationMatrix> operator*(const Numeric& x,
                                       const PropagationMatrix& pm);

/** Checks if a Propagation Matrix or something similar has good grids */
template <class T>
bool bad_propmat(const Array<T>& main,
                 const Vector& f_grid,
                 const Index sd = 1) noexcept {
  const Index nf = f_grid.nelem();
  for (auto& var : main) {
    const bool bad_stokes = sd > var.StokesDimensions();
    const bool bad_freq = nf not_eq var.NumberOfFrequencies();
    if (bad_freq or bad_stokes) return true;
  }

  return false;
}

#endif  //propagationmatrix_h
