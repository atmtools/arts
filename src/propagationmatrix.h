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
  PropagationMatrix(const Index nr_frequencies=0, const Index stokes_dim=1) :
  mfreqs(nr_frequencies), mstokes_dim(stokes_dim),  mvectortype(false) 
  {
    assert(mstokes_dim < 5 and mstokes_dim > 0);
    mdata.resize(nr_frequencies, NumberOfNeededVectors());
  }
    
  //! Initialize from a constant other
  PropagationMatrix(const PropagationMatrix& pm) :
  mfreqs(pm.mfreqs), mstokes_dim(pm.mstokes_dim), 
  mdata(pm.mdata), mvectortype(pm.mvectortype) {}
    
  PropagationMatrix(PropagationMatrix&& pm) :
  mfreqs(std::move(pm.mfreqs)), mstokes_dim(std::move(pm.mstokes_dim)),
  mdata(), mvectortype(std::move(pm.mvectortype))
  { swap(mdata, pm.mdata); }

  //! Initialize from matrix
  explicit PropagationMatrix(ConstMatrixView x, const bool& assume_fit=false) : mfreqs(1), mstokes_dim(x.ncols())
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
    
    mdata.resize(1, NumberOfNeededVectors());
    
    switch(mstokes_dim)
    {
      case 4: mdata(0, 5) = x(1, 3); mdata(0, 5) = x(2, 3); mdata(0, 3) = x(0, 3);
      case 3: mdata(0, mstokes_dim) = x(1, 2); mdata(0, 2) = x(0, 2);
      case 2: mdata(0, 1) = x(0, 1);
      case 1: mdata(0, 0) = x(0, 0);
    }
  };
  
  PropagationMatrix(const PropagationMatrix& a, const PropagationMatrix& b, const Numeric& scale = 0.5) : 
  mfreqs(a.mfreqs),
  mstokes_dim(a.mstokes_dim),
  mdata(mfreqs, mstokes_dim),
  mvectortype(false)
  {
    for(Index i = 0; i < mfreqs; i++)
    {
      switch(mstokes_dim)
      {
        case 4: 
          mdata(i, 3) = (a.mdata(i, 3) + b.mdata(i, 3)) * scale;
          mdata(i, 5) = (a.mdata(i, 5) + b.mdata(i, 5)) * scale;
          mdata(i, 6) = (a.mdata(i, 6) + b.mdata(i, 6)) * scale;
        case 3: 
          mdata(i, 2) = (a.mdata(i, 2) + b.mdata(i, 2)) * scale;
          mdata(i, mstokes_dim) = (a.mdata(i, mstokes_dim) + b.mdata(i, mstokes_dim)) * scale;
        case 2: 
          mdata(i, 1) = (a.mdata(i, 1) + b.mdata(i, 1)) * scale;
        case 1: 
          mdata(i, 0) = (a.mdata(i, 0) + b.mdata(i, 0)) * scale;
      }
    }
  };
  
  
  /*! The stokes dimension of the propagation matrix */
  Index StokesDimensions() const {return mstokes_dim;};
  
  /*! The number of frequencies of the propagation matrix */
  Index NumberOfFrequencies() const {return mfreqs;};
  
  void SetVectorType(bool vectortype) {mvectortype = vectortype;}
  
  bool IsEmpty() const {return not ((bool) mfreqs);};
  
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
  Numeric operator()(const Index iv, const Index is1, const Index is2=0) const ;
  
  void AddToDiagonalAttenuation(                  ConstVectorView add) { mdata(joker, 0) += add; }
  void AddToMainAxisLinearPolarizationAttenuation(ConstVectorView add) { mdata(joker, 1) += add; }
  void AddToOffAxisLinearPolarizationAttenuation( ConstVectorView add) { mdata(joker, 2) += add; }
  void AddToCircularPolarizationAttenuation(      ConstVectorView add) { mdata(joker, 3) += add; }
  void AddToCircularPolarizationPhaseDelay(       ConstVectorView add) { mdata(joker, mstokes_dim) += add; } 
  void AddToOffAxisLinearPolarizationPhaseDelay(  ConstVectorView add) { mdata(joker, 5) += add; }
  void AddToMainAxisLinearPolarizationPhaseDelay( ConstVectorView add) { mdata(joker, 6) += add; }
  
  void AddToSinglePointOnDiagonalAttenuation(                  const Numeric& add, const Index ifreq) { mdata(ifreq, 0) += add; }
  void AddToSinglePointOnMainAxisLinearPolarizationAttenuation(const Numeric& add, const Index ifreq) { mdata(ifreq, 1) += add; }
  void AddToSinglePointOnOffAxisLinearPolarizationAttenuation( const Numeric& add, const Index ifreq) { mdata(ifreq, 2) += add; }
  void AddToSinglePointOnCircularPolarizationAttenuation(      const Numeric& add, const Index ifreq) { mdata(ifreq, 3) += add; }
  void AddToSinglePointOnCircularPolarizationPhaseDelay(       const Numeric& add, const Index ifreq) { mdata(ifreq, mstokes_dim) += add; }
  void AddToSinglePointOnOffAxisLinearPolarizationPhaseDelay(  const Numeric& add, const Index ifreq) { mdata(ifreq, 5) += add; }
  void AddToSinglePointOnMainAxisLinearPolarizationPhaseDelay( const Numeric& add, const Index ifreq) { mdata(ifreq, 6) += add; }
  
  /*! Adds the Faraday rotation to the PropagationMatrix at required ifreq
   * 
   * No vector function exists since rot is a function of frequency
   * 
   * \param rot: rotation
   * \param ifreq: frequency index
   */
  void AddFaraday(const Numeric& rot, const Index ifreq) {mdata(ifreq, mstokes_dim) += rot;}
  void SetFaraday(const Numeric& rot, const Index ifreq) {mdata(ifreq, mstokes_dim)  = rot;}
  
  /*! Adds the Zeeman effect to the PropagationMatrix
   * 
   * \param attenuation: attenuation component
   * \param phase: phase component
   * \param nd: number density
   * \param theta: magnetic angle of los to magnetic field
   * \param eta: magnetic angle of main-axis linear polarization rotation as seen by satellite
   */
  void AddZeemanPiComponent(ConstVectorView attenuation,
                            ConstVectorView phase,
                            const Numeric& nd,
                            const Numeric& theta,
                            const Numeric& eta, 
                            ConstVectorView extra=Vector(0));
  
  void AddZeemanPiComponentThetaDerivative(ConstVectorView attenuation,
                                           ConstVectorView phase,
                                           const Numeric& nd, 
                                           const Numeric& theta,
                                           const Numeric& eta, 
                                           ConstVectorView extra=Vector(0));
  
  void AddZeemanPiComponentEtaDerivative(ConstVectorView attenuation,
                                         ConstVectorView phase,
                                         const Numeric& nd,
                                         const Numeric& theta,
                                         const Numeric& eta, 
                                         ConstVectorView extra=Vector(0));
  
  void AddZeemanPiComponentDerivative(ConstVectorView attenuation,
                                      ConstVectorView dattenuation,
                                      ConstVectorView phase,
                                      ConstVectorView dphase,
                                      const Numeric& nd,
                                      const Numeric& theta,
                                      const Numeric& dtheta,
                                      const Numeric& eta,
                                      const Numeric& deta, 
                                      ConstVectorView extra=Vector(0));
  
  /*! Adds the Zeeman effect to the PropagationMatrix
   * 
   * \param attenuation: attenuation component
   * \param phase: phase component
   * \param nd: number density
   * \param theta: magnetic angle of los to magnetic field
   * \param eta: magnetic angle of main-axis linear polarization rotation as seen by satellite
   */
  void AddZeemanSigmaPlusComponent(ConstVectorView attenuation, 
                                   ConstVectorView phase, 
                                   const Numeric& nd, 
                                   const Numeric& theta, 
                                   const Numeric& eta, 
                                   ConstVectorView extra=Vector(0));
  
  void AddZeemanSigmaPlusComponentThetaDerivative(ConstVectorView attenuation, 
                                                  ConstVectorView phase, 
                                                  const Numeric& nd, 
                                                  const Numeric& theta, 
                                                  const Numeric& eta, 
                                                  ConstVectorView extra=Vector(0));
  
  void AddZeemanSigmaPlusComponentEtaDerivative(ConstVectorView attenuation, 
                                                ConstVectorView phase, 
                                                const Numeric& nd, 
                                                const Numeric& theta, 
                                                const Numeric& eta, 
                                                ConstVectorView extra=Vector(0));
  
  void AddZeemanSigmaPlusComponentDerivative(ConstVectorView attenuation, 
                                             ConstVectorView dattenuation, 
                                             ConstVectorView phase, 
                                             ConstVectorView dphase, 
                                             const Numeric& nd, 
                                             const Numeric& theta, 
                                             const Numeric& dtheta, 
                                             const Numeric& eta,
                                             const Numeric& deta, 
                                             ConstVectorView extra=Vector(0));
  
  /*! Adds the Zeeman effect to the PropagationMatrix
   * 
   * \param attenuation: attenuation component
   * \param phase: phase component
   * \param nd: number density
   * \param theta: magnetic angle of los to magnetic field
   * \param eta: magnetic angle of main-axis linear polarization rotation as seen by satellite
   */
  void AddZeemanSigmaMinusComponent(ConstVectorView attenuation, 
                                    ConstVectorView phase, 
                                    const Numeric& nd, 
                                    const Numeric& theta, 
                                    const Numeric& eta, 
                                    ConstVectorView extra=Vector(0));
  
  void AddZeemanSigmaMinusComponentThetaDerivative(ConstVectorView attenuation, 
                                                   ConstVectorView phase, 
                                                   const Numeric& nd, 
                                                   const Numeric& theta, 
                                                   const Numeric& eta, 
                                                   ConstVectorView extra=Vector(0));
  
  void AddZeemanSigmaMinusComponentEtaDerivative(ConstVectorView attenuation, 
                                                 ConstVectorView phase, 
                                                 const Numeric& nd, 
                                                 const Numeric& theta, 
                                                 const Numeric& eta, 
                                                 ConstVectorView extra=Vector(0));
  
  void AddZeemanSigmaMinusComponentDerivative(ConstVectorView attenuation,
                                              ConstVectorView dattenuation,
                                              ConstVectorView phase,
                                              ConstVectorView dphase,
                                              const Numeric& nd, 
                                              const Numeric& theta, 
                                              const Numeric& dtheta, 
                                              const Numeric& eta,
                                              const Numeric& deta, 
                                              ConstVectorView extra=Vector(0));
  
  /*! Sets the dense matrix.  Aboid using if possible. */
  void MatrixAtFrequency(MatrixView ret, const Index ifreq) const;
  
  PropagationMatrix& operator=(PropagationMatrix&& pm)
  {
    if (this != &pm)
    {
      mfreqs = std::move(pm.mfreqs);
      mstokes_dim = std::move(pm.mstokes_dim);
      swap(mdata, pm.mdata);
      mvectortype = std::move(pm.mvectortype);
    }
    return *this;
  }

  PropagationMatrix&  operator=(const PropagationMatrix& other) { mvectortype = other.mvectortype; mstokes_dim = other.mstokes_dim; mfreqs = other.mfreqs; mdata = other.mdata; return *this; }
  PropagationMatrix&  operator=(ConstVectorView x) { for(Index i = 0; i < NumberOfNeededVectors(); i++){for(Index j = 0; j < mfreqs; j++) {mdata(j, i) = x[j];}} return *this; }
  PropagationMatrix&  operator=(const Numeric& x) { mdata  = x; return *this; }
  
  void SetAtFrequency(const Index ifreq, const PropagationMatrix& x)      { mdata(ifreq, joker)  = x.mdata(ifreq, joker); }
  void SetAtFrequency(const Index ifreq, ConstMatrixView x);
  void SetAtFrequency(const Index ifreq, const Numeric& x)      { mdata(ifreq, joker)  = x; }
  
  PropagationMatrix& operator/=(const PropagationMatrix& other) { mdata /= other.mdata; return *this; }
  PropagationMatrix& operator/=(ConstVectorView x) { for(Index i = 0; i < NumberOfNeededVectors(); i++){mdata(joker, i) /= x;} return *this; }
  PropagationMatrix& operator/=(const Numeric& x) { mdata /= x; return *this; }
  
  void DivideAtFrequency(const Index ifreq, const PropagationMatrix& x)   { mdata(ifreq, joker) /= x.mdata(ifreq, joker); }
  void DivideAtFrequency(const Index ifreq, ConstMatrixView x);
  void DivideAtFrequency(const Index ifreq, const Numeric& x)   { mdata(ifreq, joker) /= x; }
  
  PropagationMatrix& operator*=(const PropagationMatrix& other) { mdata *= other.mdata; return *this; }
  PropagationMatrix& operator*=(ConstVectorView x) { for(Index i = 0; i < NumberOfNeededVectors(); i++){mdata(joker, i) *= x;} return *this; }
  PropagationMatrix& operator*=(const Numeric& x) { mdata *= x; return *this; }
  
  void MultiplyAtFrequency(const Index ifreq, const PropagationMatrix& x) { mdata(ifreq, joker) *= x.mdata(ifreq, joker); }
  void MultiplyAtFrequency(const Index ifreq, ConstMatrixView x);
  void MultiplyAtFrequency(const Index ifreq, const Numeric& x) { mdata(ifreq, joker) *= x; }
  
  PropagationMatrix& operator+=(const PropagationMatrix& other) { mdata += other.mdata; return *this; }
  PropagationMatrix& operator+=(ConstVectorView x) { for(Index i = 0; i < NumberOfNeededVectors(); i++){mdata(joker, i) += x;} return *this; }
  PropagationMatrix& operator+=(const Numeric& x) { mdata += x; return *this; }
  
  void AddAtFrequency(const Index ifreq, const PropagationMatrix& x)      { mdata(ifreq, joker) += x.mdata(ifreq, joker); }
  void AddAtFrequency(const Index ifreq, ConstMatrixView x);
  void AddAtFrequency(const Index ifreq, const Numeric& x)      { mdata(ifreq, joker) += x; }
  
  PropagationMatrix& operator-=(const PropagationMatrix& other) { mdata -= other.mdata; return *this; }
  PropagationMatrix& operator-=(ConstVectorView x) { for(Index i = 0; i < NumberOfNeededVectors(); i++){mdata(joker, i) -= x;} return *this; }
  PropagationMatrix& operator-=(const Numeric& x) { mdata -= x; return *this; }
  
  void RemoveAtFrequency(const Index ifreq, const PropagationMatrix& x)   { mdata(ifreq, joker) -= x.mdata(ifreq, joker); }
  void RemoveAtFrequency(const Index ifreq, ConstMatrixView x);
  void RemoveAtFrequency(const Index ifreq, const Numeric& x)   { mdata(ifreq, joker) -= x; }
  
  void AddAbsorptionVectorAtFrequency(const Index ifreq, ConstVectorView x) { for(Index i = 0; i < mstokes_dim; i++) mdata(ifreq, i) += x[i]; }
  
  void AddAverageAtFrequency(const Index ifreq, ConstMatrixView mat1, ConstMatrixView mat2);
  
  void MultiplyAndAdd(const Numeric x, const PropagationMatrix& y);
  
  void MatrixInverseAtFrequency(MatrixView ret, const Index j) const;
  
  bool FittingShape(ConstMatrixView x) const;
  
  void GetTensor3(Tensor3View tensor3);
  
  VectorView Kjj(){return mdata(joker, 0);}
  VectorView K12(){return mdata(joker, 1);}
  VectorView K13(){return mdata(joker, 2);}
  VectorView K14(){return mdata(joker, 3);}
  VectorView K23(){return mdata(joker, mstokes_dim);}
  VectorView K24(){return mdata(joker, 5);}
  VectorView K34(){return mdata(joker, 6);}
  
  ConstVectorView Kjj()const{return mdata(joker, 0);}
  ConstVectorView K12()const{return mdata(joker, 1);}
  ConstVectorView K13()const{return mdata(joker, 2);}
  ConstVectorView K14()const{return mdata(joker, 3);}
  ConstVectorView K23()const{return mdata(joker, mstokes_dim);}
  ConstVectorView K24()const{return mdata(joker, 5);}
  ConstVectorView K34()const{return mdata(joker, 6);}
  
  void SetZero() {mdata = 0.0;}
  
  MatrixView GetMatrix() { return mdata; }
  ConstMatrixView GetMatrix() const { return mdata; }
  
  // Set level of case of calculations
  void CalculationCase(ArrayOfCaseOfPropagationMatrix& cases) const;
  
  // Increases level of case if too low
  void CalculationCaseMaximize(ArrayOfCaseOfPropagationMatrix& cases) const;
  
  void LeftMultiplyAtFrequency(const Index ifreq, MatrixView out, ConstMatrixView in) const;
  
  void RightMultiplyAtFrequency(const Index ifreq, MatrixView out, ConstMatrixView in) const;
    
protected:
  Index mfreqs, mstokes_dim;
  Matrix mdata;
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
                                 const PropagationMatrix& lower_level);


void compute_transmission_matrix_from_averaged_matrix_at_frequency(MatrixView T, 
                                                                   const Numeric& r, 
                                                                   const PropagationMatrix& averaged_propagation_matrix,
                                                                   const Index ifreq);

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
                                                const Array<PropagationMatrix>& dprop_mat_lower_level);

std::ostream& operator<<(std::ostream& os, const PropagationMatrix& pm);
std::ostream& operator<<(std::ostream& os, const ArrayOfPropagationMatrix& apm);
std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfPropagationMatrix& aapm);

class StokesVector: public PropagationMatrix
{
public:
  
  // Initialize variable size depending on requirements...
  StokesVector(const Index nr_frequencies=0, const Index stokes_dim=1)
  {
    mvectortype = true;
    mfreqs = nr_frequencies;
    mstokes_dim = stokes_dim;
    assert(mstokes_dim < 5 and mstokes_dim > 0);
    mdata.resize(mfreqs, mstokes_dim);
  };
  
  explicit StokesVector(ConstVectorView x)
  {
    mfreqs = 1;
    mstokes_dim = x.nelem();
    assert(mstokes_dim < 5 and mstokes_dim > 0);
    mvectortype = true;
    mdata.resize(1, mstokes_dim);
    for(Index i = 0; i < mstokes_dim; i++)
      mdata(0, i) = x[i];
  };
  
  StokesVector(const StokesVector& a, const StokesVector& b, const Numeric& scale = 0.5)
  {
    mfreqs = a.NumberOfFrequencies();
    mstokes_dim = a.StokesDimensions();
    mdata.resize(mfreqs, mstokes_dim);
    mvectortype = true;
    
    for(Index i = 0; i < mfreqs; i++)
    {
      switch(mstokes_dim)
      {
        case 4:  mdata(i, 3) = (a.mdata(i, 3) + b.mdata(i, 3)) * scale;
        case 3:  mdata(i, 2) = (a.mdata(i, 2) + b.mdata(i, 2)) * scale;
        case 2:  mdata(i, 1) = (a.mdata(i, 1) + b.mdata(i, 1)) * scale;
        case 1:  mdata(i, 0) = (a.mdata(i, 0) + b.mdata(i, 0)) * scale;
      }
    }
  };
  
  StokesVector& operator+=(const PropagationMatrix& x)
  {
    mdata += x.GetMatrix()(joker, Range(0, mstokes_dim, 1));
    return *this;
  }
  
  StokesVector& operator=(const PropagationMatrix& x)
  {
    mstokes_dim = x.StokesDimensions();
    mfreqs = x.NumberOfFrequencies();
    mdata = x.GetMatrix()(joker, Range(0, mstokes_dim, 1));
    return *this;
  }
  
  StokesVector& operator=(const Numeric& x)
  {
    mdata = x;
    return *this;
  }
  
  void MultiplyAndAdd(const Numeric x, const StokesVector& y)
  {
    assert(mstokes_dim == y.mstokes_dim);
    assert(mfreqs == y.mfreqs);
    
    for(Index i = 0; i < mfreqs; i++)
    {
      switch(mstokes_dim)
      {
        case 4: mdata(i, 3) += x * y.mdata(i, 3);
        case 3: mdata(i, 2) += x * y.mdata(i, 2);
        case 2: mdata(i, 1) += x * y.mdata(i, 1);
        case 1: mdata(i, 0) += x * y.mdata(i, 0);
      }
    }
  }
  
  VectorView VectorAtFrequency(const Index ifreq) { return mdata(ifreq, joker); }
  ConstVectorView VectorAtFrequency(const Index ifreq) const { return mdata(ifreq, joker); }
  
  void VectorAtFrequency(VectorView ret, const Index ifreq) { ret = mdata(ifreq, joker); }
  
  void VectorAtFrequency(VectorView ret, const Index ifreq) const { ret = mdata(ifreq, joker); }
  
  void SetAtFrequency(const Index ifreq, ConstVectorView x)  { mdata(ifreq, joker) = x; }
  
  void AddAverageAtFrequency(const Index ifreq, ConstVectorView vec1, ConstVectorView vec2)
  {
    switch(mstokes_dim)
    {
      case 4: mdata(ifreq, 3) += (vec1[3] + vec2[3]) * 0.5;
      case 3: mdata(ifreq, 2) += (vec1[2] + vec2[2]) * 0.5;
      case 2: mdata(ifreq, 1) += (vec1[1] + vec2[1]) * 0.5;
      case 1: mdata(ifreq, 0) += (vec1[0] + vec2[0]) * 0.5; 
    }
  }
  
  bool IsPolarized(Index ifreq) const
  {
    switch(mstokes_dim)
    {
      case 4: if(K14()[ifreq] not_eq 0.0) return true;
      case 3: if(K13()[ifreq] not_eq 0.0) return true;
      case 2: if(K12()[ifreq] not_eq 0.0) return true;
    }
    return false;
  }
  
  bool IsUnpolarized(Index ifreq) const
  {
    return not IsPolarized(ifreq);
  }
};

typedef Array<StokesVector> ArrayOfStokesVector;
typedef Array<ArrayOfStokesVector> ArrayOfArrayOfStokesVector;
typedef Array<ArrayOfArrayOfStokesVector> ArrayOfArrayOfArrayOfStokesVector;

std::ostream& operator<<(std::ostream& os, const StokesVector& pm);
std::ostream& operator<<(std::ostream& os, const ArrayOfStokesVector& apm);
std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfStokesVector& aapm);

void get_diydx_replacement(MatrixView diydx_this,
                           MatrixView diydx_next,
                           ConstMatrixView iy,
                           ConstMatrixView sibi,
                           const StokesVector& nlte_this,
                           const StokesVector& nlte_next,
                           const StokesVector& dnltedx_this,
                           const StokesVector& dnltedx_next,
                           const PropagationMatrix& K_this,
                           const PropagationMatrix& K_next,
                           const PropagationMatrix& dKdx_this,
                           const PropagationMatrix& dKdx_next,
                           ConstTensor3View T_this,
                           ConstTensor3View dTdx_this,
                           ConstTensor3View dTdx_next,
                           ConstTensor3View PiT_this,
                           ConstTensor3View PiT_next,
                           const Numeric& temperature_this,
                           const Numeric& temperature_next,
                           const Numeric& dt,
                           ConstVectorView dBdx_this,
                           ConstVectorView dBdx_next,
                           const Numeric& r,
                           const bool& do_Bsource,
                           const bool& do_HSE,
                           const bool& do_nonlte );

#endif //propagationmatrix_h
