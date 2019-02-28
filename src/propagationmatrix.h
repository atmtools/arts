/* Copyright (C) 2017
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/*!
 * \file   propagationmatrix.h
 * \brief  Stuff related to the propagation matrix.
 * 
 * This implementation takes advantage of symmetries to lower memory and computational costs
 * 
 * Some standard functions applied on the propagation matrix have been included but far from all.
 * 
 * Present support is for stokes dim 1-4 however please take note that  circular polarization delay variable's
 * position in the internal mdata structure moves depending on if stokes_dim is 3 or 4.  It is not present elsewhere...
 * 
 * \author Richard Larsson
 * \date   2017-06-23
 */

#ifndef propagationmatrix_h
#define propagationmatrix_h

#include "matpackIV.h"
#include "complex.h"

class PropagationMatrix;

class LazyPropagationMatrixScale {
public:
  LazyPropagationMatrixScale(const PropagationMatrix& p, const Numeric& x) : pm(p), scale(x) {}
  const PropagationMatrix& pm;
  const Numeric& scale;
};

enum class CaseOfPropagationMatrix
{
  Diagonal=1,
  FullDimensional=100 // Must always be at end even if we find out faster calcs for other common matrices
};

typedef Array<CaseOfPropagationMatrix> ArrayOfCaseOfPropagationMatrix;
typedef Array<ArrayOfCaseOfPropagationMatrix> ArrayOfArrayOfCaseOfPropagationMatrix;

/*! Propagation Matrix Holder Class With Some Computational Capabilities*/
class PropagationMatrix
{

public:
  
  //! Initialize variable sizes
  PropagationMatrix(const Index nr_frequencies=0, const Index stokes_dim=1, const Index nr_za=1, const Index nr_aa=1, const Numeric v=0.0) :
  mfreqs(nr_frequencies), mstokes_dim(stokes_dim), mza(nr_za), maa(nr_aa), mvectortype(false) 
  {
    assert(mstokes_dim < 5 and mstokes_dim > 0);
    mdata = Tensor4(maa, mza, mfreqs, NumberOfNeededVectors(), v);
  }
    
  //! Initialize from a constant other
  PropagationMatrix(const PropagationMatrix& pm) :
  mfreqs(pm.mfreqs), mstokes_dim(pm.mstokes_dim), 
  mza(pm.mza), maa(pm.maa), 
  mdata(pm.mdata), mvectortype(pm.mvectortype) {}
    
  PropagationMatrix(PropagationMatrix&& pm) :
  mfreqs(std::move(pm.mfreqs)), mstokes_dim(std::move(pm.mstokes_dim)),
  mza(std::move(pm.mza)), maa(std::move(pm.maa)),
  mdata(), mvectortype(std::move(pm.mvectortype))
  { swap(mdata, pm.mdata); }
  
  explicit PropagationMatrix(ConstTensor4View x) : mfreqs(x.nrows()), mza(x.npages()), maa(x.nbooks()), mdata(x), mvectortype(false)
  {
    switch(x.ncols()) {
      case 7: mstokes_dim=4; break;
      case 4: mstokes_dim=3; break;
      case 2: mstokes_dim=2; break;
      case 1: mstokes_dim=1; break;
      default: throw std::runtime_error("Tensor4 not representative of PropagationMatrix");
    }
  }

  //! Initialize from matrix
  explicit PropagationMatrix(ConstMatrixView x, const bool& assume_fit=false) : mfreqs(1), mstokes_dim(x.ncols()), mza(1), maa(1)
  {
    assert(mstokes_dim < 5 and mstokes_dim > 0);
    mvectortype = false;
    
    if(not assume_fit)
    {
      if(not FittingShape(x))
      {
        throw std::runtime_error("Matrix not fit as propagation matrix");
      }
    }
    
    mdata.resize(1, 1, 1, NumberOfNeededVectors());
    
    switch(mstokes_dim)
    {
      case 4: mdata(0, 0, 0, 5) = x(1, 3); mdata(0, 0, 0, 5) = x(2, 3); mdata(0, 0, 0, 3) = x(0, 3); /* FALLTHROUGH */
      case 3: mdata(0, 0, 0, mstokes_dim) = x(1, 2); mdata(0, 0, 0, 2) = x(0, 2); /* FALLTHROUGH */
      case 2: mdata(0, 0, 0, 1) = x(0, 1); /* FALLTHROUGH */
      case 1: mdata(0, 0, 0, 0) = x(0, 0); /* FALLTHROUGH */
    }
  };
  
  PropagationMatrix(const PropagationMatrix& a, const PropagationMatrix& b, const Numeric& scale = 0.5) : 
  mfreqs(a.mfreqs),
  mstokes_dim(a.mstokes_dim),
  mza(a.mza),
  maa(a.maa),
  mdata(maa, mza, mfreqs, mstokes_dim),
  mvectortype(false)
  {
    for(Index i = 0; i < maa; i++)
    for(Index j = 0; j < mza; j++)
    for(Index k = 0; k < mfreqs; k++)
    {
      switch(mstokes_dim)
      {
        case 4: 
          mdata(i, j, k, 3) = (a.mdata(i, j, k, 3) + b.mdata(i, j, k, 3)) * scale;
          mdata(i, j, k, 5) = (a.mdata(i, j, k, 5) + b.mdata(i, j, k, 5)) * scale;
          mdata(i, j, k, 6) = (a.mdata(i, j, k, 6) + b.mdata(i, j, k, 6)) * scale; /* FALLTHROUGH */
        case 3: 
          mdata(i, j, k, 2) = (a.mdata(i, j, k, 2) + b.mdata(i, j, k, 2)) * scale;
          mdata(i, j, k, mstokes_dim) = (a.mdata(i, j, k, mstokes_dim) + b.mdata(i, j, k, mstokes_dim)) * scale; /* FALLTHROUGH */
        case 2: 
          mdata(i, j, k, 1) = (a.mdata(i, j, k, 1) + b.mdata(i, j, k, 1)) * scale; /* FALLTHROUGH */
        case 1: 
          mdata(i, j, k, 0) = (a.mdata(i, j, k, 0) + b.mdata(i, j, k, 0)) * scale; /* FALLTHROUGH */
      }
    }
  };
  
  
  /*! The stokes dimension of the propagation matrix */
  Index StokesDimensions() const {return mstokes_dim;};
  
  /*! The number of frequencies of the propagation matrix */
  Index NumberOfFrequencies() const {return mfreqs;};
  
  /*! The number of frequencies of the propagation matrix */
  Index NumberOfZenithAngles() const {return mza;};
  
  /*! The number of frequencies of the propagation matrix */
  Index NumberOfAzimuthAngles() const {return maa;};
  
  void SetVectorType(bool vectortype) {mvectortype = vectortype;}
  
  bool IsEmpty() const {return not mfreqs or not mza or not maa;};
  
  bool IsZero(const Index iv=0, const Index iz=0, const Index ia=0) const 
  {
    for(auto& n : mdata(ia, iz, iv, joker))
      if(n not_eq 0)
        return false;
    return true;
  };
  
  bool IsRotational(const Index iv=0, const Index iz=0, const Index ia=0) const 
  {
    if(mdata(ia, iz, iv, 0) == 0.0) 
      return true;
    else 
      return false;
  };
  
  /* The number of required vectors to fill this PropagationMatrix --- designed for GetVector(i) */
  Index NumberOfNeededVectors() const
  {
    if(not mvectortype)
    {
      switch(mstokes_dim)
      {
        case 1: return 1; break;
        case 2: return 2; break;
        case 3: return 4; break;
        case 4: return 7; break;
        default: throw std::runtime_error("Cannot understand the input in PropagationMatrix");
      }
    }
    else
    {
      return mstokes_dim;
    }
  }
  
  /*! access operator.  Please refrain from using this if possible */
  Numeric operator()(const Index iv=0, const Index is1=0, const Index is2=0, const Index iz=0, const Index ia=0) const ;
  
  /*! Adds the Faraday rotation to the PropagationMatrix at required ifreq
   * 
   * No vector function exists since rot is a function of frequency
   * 
   * \param rot: rotation
   * \param ifreq: frequency index
   */
  void AddFaraday(const Numeric& rot, const Index iv=0, const Index iz=0, const Index ia=0) {mdata(ia, iz, iv, mstokes_dim) += rot;}
  void SetFaraday(const Numeric& rot, const Index iv=0, const Index iz=0, const Index ia=0) {mdata(ia, iz, iv, mstokes_dim)  = rot;}
  
  /*! Sets the dense matrix.  Avoid using if possible. */
  void MatrixAtPosition(MatrixView ret,
                        const Index iv=0, const Index iz=0, const Index ia=0)
    const;
  
  PropagationMatrix& operator=(PropagationMatrix&& pm)
  {
    if (this != &pm)
    {
      mfreqs = std::move(pm.mfreqs);
      mstokes_dim = std::move(pm.mstokes_dim);
      mza = std::move(pm.mza);
      maa = std::move(pm.maa);
      swap(mdata, pm.mdata);
      mvectortype = std::move(pm.mvectortype);
    }
    return *this;
  }
  
  PropagationMatrix& operator=(const LazyPropagationMatrixScale& lpms)
  {
    operator=(lpms.pm);
    mdata *= lpms.scale;
    return *this;
  }
  PropagationMatrix&  operator=(const PropagationMatrix& other)
  {
    mvectortype = other.mvectortype;
    mstokes_dim = other.mstokes_dim;
    mfreqs = other.mfreqs;
    mza = other.mza;
    maa = other.maa;
    mdata = other.mdata;
    return *this;
  }
  PropagationMatrix&  operator=(ConstVectorView x) 
  { 
    for(Index i = 0; i < NumberOfNeededVectors(); i++){
      for(Index j = 0; j < mza; j++){
        for(Index k = 0; k < maa; k++){
          mdata(k, j, joker, i) = x;}}}
    return *this; 
  }
  PropagationMatrix&  operator=(const Numeric& x)
    { mdata = x; return *this; }
  
  void SetAtPosition(const PropagationMatrix& x,
                     const Index iv=0, const Index iz=0, const Index ia=0)
    { mdata(ia, iz, iv, joker)  = x.mdata(ia, iz, iv, joker); }
  void SetAtPosition(ConstMatrixView x,
                     const Index iv=0, const Index iz=0, const Index ia=0);
  void SetAtPosition(const Numeric& x,
                     const Index iv=0, const Index iz=0, const Index ia=0)
    { mdata(ia, iz, iv, joker)  = x; }
  
  PropagationMatrix& operator/=(const PropagationMatrix& other)
    { mdata /= other.mdata; return *this; }
  PropagationMatrix& operator/=(ConstVectorView x)
  { 
    for(Index i = 0; i < NumberOfNeededVectors(); i++){
      for(Index j = 0; j < mza; j++){
        for(Index k = 0; k < maa; k++){
          mdata(k, j, joker, i) /= x;}}}
          return *this; 
  }
  PropagationMatrix& operator/=(const Numeric& x)
    { mdata /= x; return *this; }
  
  void DivideAtPosition(const PropagationMatrix& x,
                        const Index iv=0, const Index iz=0, const Index ia=0)
    { mdata(ia, iz, iv, joker) /= x.mdata(ia, iz, iv, joker); }
  void DivideAtPosition(ConstMatrixView x,
                        const Index iv=0, const Index iz=0, const Index ia=0);
  void DivideAtPosition(const Numeric& x,
                        const Index iv=0, const Index iz=0, const Index ia=0)
    { mdata(ia, iz, iv, joker) /= x; }
  
  PropagationMatrix& operator*=(const PropagationMatrix& other)
    { mdata *= other.mdata; return *this; }
  PropagationMatrix& operator*=(ConstVectorView x) 
  { 
    for(Index i = 0; i < NumberOfNeededVectors(); i++)
      for(Index j = 0; j < mza; j++)
        for(Index k = 0; k < maa; k++)
          mdata(k, j, joker, i) *= x;
    return *this; 
  }
  PropagationMatrix& operator*=(const Numeric& x)
    { mdata *= x; return *this; }
  
  void MultiplyAtPosition(const PropagationMatrix& x,
                          const Index iv=0, const Index iz=0, const Index ia=0)
    { mdata(ia, iz, iv, joker) *= x.mdata(ia, iz, iv, joker); }
  void MultiplyAtPosition(ConstMatrixView x,
                          const Index iv=0, const Index iz=0, const Index ia=0);
  void MultiplyAtPosition(const Numeric& x,
                          const Index iv=0, const Index iz=0, const Index ia=0)
    { mdata(ia, iz, iv, joker) *= x; }
  
  PropagationMatrix& operator+=(const PropagationMatrix& other)
  { mdata += other.mdata; return *this; }
  PropagationMatrix& operator+=(const LazyPropagationMatrixScale& lpms)
  { MultiplyAndAdd(lpms.scale, lpms.pm); return *this; }
  PropagationMatrix& operator+=(ConstVectorView x) 
  { 
    for(Index i = 0; i < NumberOfNeededVectors(); i++)
      for(Index j = 0; j < mza; j++)
        for(Index k = 0; k < maa; k++)
          mdata(k, j, joker, i) += x;
    return *this; 
  }
  PropagationMatrix& operator+=(const Numeric& x)
    { mdata += x; return *this; }
  
  void AddAtPosition(const PropagationMatrix& x,
                     const Index iv=0, const Index iz=0, const Index ia=0)
    { mdata(ia, iz, iv, joker) += x.mdata(ia, iz, iv, joker); }
  void AddAtPosition(ConstMatrixView x,
                     const Index iv=0, const Index iz=0, const Index ia=0);
  void AddAtPosition(const Numeric& x,
                     const Index iv=0, const Index iz=0, const Index ia=0)
    { mdata(ia, iz, iv, joker) += x; }
  
  PropagationMatrix& operator-=(const PropagationMatrix& other)
    { mdata -= other.mdata; return *this; }
  PropagationMatrix& operator-=(ConstVectorView x) 
  { 
    for(Index i = 0; i < NumberOfNeededVectors(); i++){
      for(Index j = 0; j < mza; j++){
        for(Index k = 0; k < maa; k++){
          mdata(k, j, joker, i) -= x;}}}
          return *this; 
  }
  PropagationMatrix& operator-=(const Numeric& x)
    { mdata -= x; return *this; }
  
  void RemoveAtPosition(const PropagationMatrix& x,
                        const Index iv=0, const Index iz=0, const Index ia=0)
    { mdata(ia, iz, iv, joker) -= x.mdata(ia, iz, iv, joker); }
  void RemoveAtPosition(ConstMatrixView x,
                        const Index iv=0, const Index iz=0, const Index ia=0);
  void RemoveAtPosition(const Numeric& x,
                        const Index iv=0, const Index iz=0, const Index ia=0)
    { mdata(ia, iz, iv, joker) -= x; }
  
  void AddAbsorptionVectorAtPosition(ConstVectorView x,
                                     const Index iv=0, const Index iz=0, const Index ia=0)
    { for(Index i = 0; i < mstokes_dim; i++) mdata(ia, iz, iv, i) += x[i]; }
  
  void AddAverageAtPosition(ConstMatrixView mat1, ConstMatrixView mat2,
                            const Index iv=0, const Index iz=0, const Index ia=0);
  
  void MultiplyAndAdd(const Numeric x, const PropagationMatrix& y);
  
  void MatrixInverseAtPosition(MatrixView ret,
                               const Index iv=0, const Index iz=0, const Index ia=0)
    const;
  
  bool FittingShape(ConstMatrixView x) const;
  
  void GetTensor3(Tensor3View tensor3, const Index iz=0, const Index ia=0);
  
  VectorView Kjj(const Index iz=0, const Index ia=0){return mdata(ia, iz, joker, 0);}
  VectorView K12(const Index iz=0, const Index ia=0){return mdata(ia, iz, joker, 1);}
  VectorView K13(const Index iz=0, const Index ia=0){return mdata(ia, iz, joker, 2);}
  VectorView K14(const Index iz=0, const Index ia=0){return mdata(ia, iz, joker, 3);}
  VectorView K23(const Index iz=0, const Index ia=0){return mdata(ia, iz, joker, mstokes_dim);}
  VectorView K24(const Index iz=0, const Index ia=0){return mdata(ia, iz, joker, 5);}
  VectorView K34(const Index iz=0, const Index ia=0){return mdata(ia, iz, joker, 6);}
  
  ConstVectorView Kjj(const Index iz=0, const Index ia=0)const{return mdata(ia, iz, joker, 0);}
  ConstVectorView K12(const Index iz=0, const Index ia=0)const{return mdata(ia, iz, joker, 1);}
  ConstVectorView K13(const Index iz=0, const Index ia=0)const{return mdata(ia, iz, joker, 2);}
  ConstVectorView K14(const Index iz=0, const Index ia=0)const{return mdata(ia, iz, joker, 3);}
  ConstVectorView K23(const Index iz=0, const Index ia=0)const{return mdata(ia, iz, joker, mstokes_dim);}
  ConstVectorView K24(const Index iz=0, const Index ia=0)const{return mdata(ia, iz, joker, 5);}
  ConstVectorView K34(const Index iz=0, const Index ia=0)const{return mdata(ia, iz, joker, 6);}
  
  void SetZero() {mdata = 0.0;}
  
  Tensor4View GetData() { return mdata; }
  ConstTensor4View GetData() const { return mdata; }
  
  // Set level of case of calculations
  void CalculationCase(ArrayOfCaseOfPropagationMatrix& cases) const;
  
  // Increases level of case if too low
  void CalculationCaseMaximize(ArrayOfCaseOfPropagationMatrix& cases) const;
  
  void LeftMultiplyAtPosition(MatrixView out, ConstMatrixView in,
                              const Index iv=0, const Index iz=0, const Index ia=0)
    const;
  
  void RightMultiplyAtPosition(MatrixView out, ConstMatrixView in,
                               const Index iv=0, const Index iz=0, const Index ia=0)
    const;
    
protected:
  Index mfreqs, mstokes_dim;
  Index mza, maa;
  Tensor4 mdata;
  bool mvectortype;
};

typedef Array<PropagationMatrix> ArrayOfPropagationMatrix;
typedef Array<ArrayOfPropagationMatrix> ArrayOfArrayOfPropagationMatrix;


/*! Compute the matrix exponent as the transmission matrix of this propagation matrix
 * 
 * The propagation matrix is multiplied by -r and level-averaged before exponent is applied.
 * 
 * upper_level and lower_level propagation matrices should thus be the level matrices and r the distance between these
 * levels.  The same is true for the derivative matrices.
 * 
 * Stokes dim 1 and 4 have been tested more.  Stokes dim 2 and 3 have been found to work but could still have hidden errors for uncommon cases
 * 
 * \param T: transmission matrix with outmost dimension being frequency
 * \param r: the distance over which the propagation matrix causes the transmission
 */
void compute_transmission_matrix(Tensor3View T, 
                                 const Numeric& r, 
                                 const PropagationMatrix& upper_level, 
                                 const PropagationMatrix& lower_level,
                                 const Index iz=0,
                                 const Index ia=0);


void compute_transmission_matrix_from_averaged_matrix_at_frequency(MatrixView T, 
                                                                   const Numeric& r, 
                                                                   const PropagationMatrix& averaged_propagation_matrix,
                                                                   const Index iv,
                                                                   const Index iz=0,
                                                                   const Index ia=0);

/*! Compute the matrix exponent as the transmission matrix of this propagation matrix
 * 
 * The propagation matrix is multiplied by -r and level-averaged before exponent is applied.
 * 
 * upper_level and lower_level propagation matrices should thus be the level matrices and r the distance between these
 * levels.  The same is true for the derivative matrices.
 * 
 * Stokes dim 1 and 4 have been tested more.  Stokes dim 2 and 3 have been found to work but could still have hidden errors for uncommon cases
 * 
 * \param T: transmission matrix with outmost dimension being frequency
 * \param dT_upp: transmission matrix derivative with respect to derivatives of the propagation matrix for upper level
 * \param dT_low: transmission matrix derivative with respect to derivatives of the propagation matrix for lower level
 * \param r: the distance over which the propagation matrix causes the transmission
 * \param dprop_mat_upp: derivatives of the upper propagation matrix with respect to some parameter (is multiplied by -0.5 r)
 * \param dprop_mat_low: derivatives of the lower propagation matrix with respect to some parameter (is multiplied by -0.5 r)
 */
void compute_transmission_matrix_and_derivative(Tensor3View T, 
                                                Tensor4View dT_upper_level, 
                                                Tensor4View dT_lower_level, 
                                                const Numeric& r, 
                                                const PropagationMatrix& upper_level, 
                                                const PropagationMatrix& lower_level, 
                                                const Array<PropagationMatrix>& dprop_mat_upper_level, 
                                                const Array<PropagationMatrix>& dprop_mat_lower_level,
                                                const Numeric& dr_dTu=0.0,
                                                const Numeric& dr_dTl=0.0,
                                                const Index it=-1,
                                                const Index iz=0,
                                                const Index ia=0);

std::ostream& operator<<(std::ostream& os, const PropagationMatrix& pm);
std::ostream& operator<<(std::ostream& os, const ArrayOfPropagationMatrix& apm);
std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfPropagationMatrix& aapm);

class StokesVector: public PropagationMatrix
{
public:
  
  // Initialize variable size depending on requirements...
  StokesVector(const Index nr_frequencies=0, const Index stokes_dim=1, const Index nr_za=1, const Index nr_aa=1, const Numeric& v=0.0)
  {
    mvectortype = true;
    mfreqs = nr_frequencies;
    mstokes_dim = stokes_dim;
    mza = nr_za;
    maa = nr_aa;
    assert(mstokes_dim < 5 and mstokes_dim > 0);
    mdata = Tensor4(maa, mza, mfreqs, mstokes_dim, v);
  };
  
  explicit StokesVector(ConstTensor4View x) 
  {
    mfreqs = x.nrows();
    mstokes_dim = x.ncols();
    mza = x.npages();
    maa = x.nbooks();
    mdata = x;
    mvectortype = true;
    if(mstokes_dim > 4 or mstokes_dim < 1)
      throw std::runtime_error("Tensor4 is bad for StokesVector");
  }
  
  explicit StokesVector(ConstVectorView x)
  {
    mfreqs = 1;
    mstokes_dim = x.nelem();
    mza = 1;
    maa = 1;
    assert(mstokes_dim < 5 and mstokes_dim > 0);
    mvectortype = true;
    mdata.resize(1, 1, 1, mstokes_dim);
    for(Index i = 0; i < mstokes_dim; i++)
      mdata(0, 0, 0, i) = x[i];
  };
  
  StokesVector(const StokesVector& a, const StokesVector& b, const Numeric& scale = 0.5)
  {
    mfreqs = a.NumberOfFrequencies();
    mstokes_dim = a.StokesDimensions();
    mza = a.NumberOfZenithAngles();
    maa = a.NumberOfAzimuthAngles();
    mdata.resize(maa, mza, mfreqs, mstokes_dim);
    mvectortype = true;
    
    for(Index i = 0; i < maa; i++)
    for(Index j = 0; j < mza; j++)
    for(Index k = 0; k < mfreqs; k++)
    {
      switch(mstokes_dim)
      {
        case 4:  mdata(i,j,k, 3) = (a.mdata(i,j,k, 3) + b.mdata(i,j,k, 3)) * scale; /* FALLTHROUGH */
        case 3:  mdata(i,j,k, 2) = (a.mdata(i,j,k, 2) + b.mdata(i,j,k, 2)) * scale; /* FALLTHROUGH */
        case 2:  mdata(i,j,k, 1) = (a.mdata(i,j,k, 1) + b.mdata(i,j,k, 1)) * scale; /* FALLTHROUGH */
        case 1:  mdata(i,j,k, 0) = (a.mdata(i,j,k, 0) + b.mdata(i,j,k, 0)) * scale; /* FALLTHROUGH */
      }
    }
  };
  
  /* The number of required vectors to fill this StokesVector */
  Index NumberOfNeededVectors() const { return mstokes_dim; }
  
  StokesVector& operator+=(const PropagationMatrix& x)
  {
    mdata += x.GetData()(joker, joker, joker, Range(0, mstokes_dim, 1));
    return *this;
  }
  
  StokesVector& operator+=(const LazyPropagationMatrixScale& lpms)
  { MultiplyAndAdd(lpms.scale, lpms.pm); return *this; }
  
  StokesVector& operator=(const PropagationMatrix& x)
  {
    mstokes_dim = x.StokesDimensions();
    mfreqs = x.NumberOfFrequencies();
    mza = x.NumberOfZenithAngles();
    maa = x.NumberOfAzimuthAngles();
    mdata = x.GetData()(joker, joker, joker, Range(0, mstokes_dim, 1));
    return *this;
  }
  
  StokesVector& operator=(const LazyPropagationMatrixScale& lpms)
  {
    operator=(lpms.pm);
    mdata *= lpms.scale;
    return *this;
  }
  
  StokesVector& operator=(const Numeric& x)
  {
    mdata = x;
    return *this;
  }
  
  void MultiplyAndAdd(const Numeric x, const PropagationMatrix& y)
  {
    assert(mstokes_dim == y.StokesDimensions());
    assert(mfreqs == y.NumberOfFrequencies());
    assert(mza == y.NumberOfZenithAngles());
    assert(maa == y.NumberOfAzimuthAngles());
    
    const ConstTensor4View data = y.GetData();
    
    for(Index i = 0; i < maa; i++)
    for(Index j = 0; j < mza; j++)
    for(Index k = 0; k < mfreqs; k++) {
      switch(mstokes_dim) {
        case 4: mdata(i,j,k, 3) += x * data(i,j,k, 3); /* FALLTHROUGH */
        case 3: mdata(i,j,k, 2) += x * data(i,j,k, 2); /* FALLTHROUGH */
        case 2: mdata(i,j,k, 1) += x * data(i,j,k, 1); /* FALLTHROUGH */
        case 1: mdata(i,j,k, 0) += x * data(i,j,k, 0); /* FALLTHROUGH */
      }
    }
  }
  
  VectorView VectorAtPosition(const Index iv=0, const Index iz=0, const Index ia=0) { return mdata(ia, iz, iv, joker); }
  ConstVectorView VectorAtPosition(const Index iv=0, const Index iz=0, const Index ia=0) const { return mdata(ia, iz, iv, joker); }
  
  void VectorAtPosition(VectorView ret, const Index iv=0, const Index iz=0, const Index ia=0) { ret = mdata(ia, iz, iv, joker); }
  void VectorAtPosition(VectorView ret, const Index iv=0, const Index iz=0, const Index ia=0) const { ret = mdata(ia, iz, iv, joker); }
  
  void SetAtPosition(ConstVectorView x, const Index iv=0, const Index iz=0, const Index ia=0)  { mdata(ia, iz, iv, joker) = x; }
  
  void AddAverageAtPosition(ConstVectorView vec1, ConstVectorView vec2, const Index iv=0, const Index iz=0, const Index ia=0)
  {
    switch(mstokes_dim)
    {
      case 4: mdata(ia, iz , iv, 3) += (vec1[3] + vec2[3]) * 0.5; /* FALLTHROUGH */
      case 3: mdata(ia, iz , iv, 2) += (vec1[2] + vec2[2]) * 0.5; /* FALLTHROUGH */
      case 2: mdata(ia, iz , iv, 1) += (vec1[1] + vec2[1]) * 0.5; /* FALLTHROUGH */
      case 1: mdata(ia, iz , iv, 0) += (vec1[0] + vec2[0]) * 0.5; /* FALLTHROUGH */
    }
  }
  
  bool IsPolarized(const Index iv=0, const Index iz=0, const Index ia=0) const
  {
    switch(mstokes_dim)
    {
      case 4: if(K14(iz, ia)[iv] not_eq 0.0) return true; /* FALLTHROUGH */
      case 3: if(K13(iz, ia)[iv] not_eq 0.0) return true; /* FALLTHROUGH */
      case 2: if(K12(iz, ia)[iv] not_eq 0.0) return true; /* FALLTHROUGH */
    }
    return false;
  }
  
  bool IsUnpolarized(const Index iv=0, const Index iz=0, const Index ia=0) const
  {
    return not IsPolarized(iv, iz, ia);
  }
};

typedef Array<StokesVector> ArrayOfStokesVector;
typedef Array<ArrayOfStokesVector> ArrayOfArrayOfStokesVector;
typedef Array<ArrayOfArrayOfStokesVector> ArrayOfArrayOfArrayOfStokesVector;

std::ostream& operator<<(std::ostream& os, const StokesVector& pm);
std::ostream& operator<<(std::ostream& os, const ArrayOfStokesVector& apm);
std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfStokesVector& aapm);

inline LazyPropagationMatrixScale operator*(const PropagationMatrix& pm, const Numeric& x) {return LazyPropagationMatrixScale(pm, x);}
inline LazyPropagationMatrixScale operator*(const Numeric& x, const PropagationMatrix& pm) {return LazyPropagationMatrixScale(pm, x);}

#endif //propagationmatrix_h
