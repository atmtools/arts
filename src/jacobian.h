/* Copyright (C) 2004-2012 Mattias Ekstrom <ekstrom@rss.chalmers.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/** \file
    Declarations required for the calculation of jacobians.

    \author Mattias Ekstrom
*/

#ifndef jacobian_h
#define jacobian_h

#include <map>
#include <iostream>
#include <stdexcept>
#include "matpackI.h"
#include "array.h"
#include "mystring.h"
#include "bifstream.h"
#include "interpolation.h"
#include "logic.h"
#include "methods.h"
#include "ppath.h"
#include "agenda_class.h"
#include "abs_species_tags.h"

#include "quantum.h"

/** Contains the data for one retrieval quantity.
    \author Mattias Ekstrom */
class RetrievalQuantity {
public:

  /** Default constructor. Needed by make_array. */
  RetrievalQuantity() : mmaintag(),
                        msubtag(),
                        msubsubtag(),
                        mmode(),
                        manalytical(-1),
                        mperturbation(0.),
                        mgrids(),
                        mquantumidentifier(),
                        mintegration_flag(false)
  { /* Nothing to do here. */ }


  /** Constructor that sets the values. */
  RetrievalQuantity(const String&             maintag,
                    const String&             subtag,
                    const String&             subsubtag,
                    const String&             mode,
                    const Index&              analytical,
                    const Numeric&            perturbation,
                    const ArrayOfVector&      grids ) :
    mmaintag(maintag),
    msubtag(subtag),
    msubsubtag(subsubtag),
    mmode(mode),
    manalytical(analytical),
    mperturbation(perturbation),
    mgrids(grids),
    mquantumidentifier(),
    mintegration_flag(false)
  {
    // With Matpack, initialization of mgrids from grids should work correctly.
  }

  /** Main tag. */
  const String& MainTag() const { return mmaintag; }
  void MainTag( const String& mt ) { mmaintag = mt; }
  /** Subtag. Eg. for gas species: O3, ClO. */
  const String& Subtag() const { return msubtag; }
  void Subtag( const String& st ) { msubtag = st; }
  /** SubSubtag. Eg. for scat species fields: mass_density, mass_flux, ... */
  const String& SubSubtag() const { return msubsubtag; }
  void SubSubtag( const String& sst ) { msubsubtag = sst; }
  /** Calculation mode. Eg. "abs", "rel", "vmr", "nd", "From propagation matrix". 
       Note that the latter of these only supports "vmr" for abs species. */
  const String& Mode() const { return mmode; }
  void Mode( const String& m ) { mmode = m; }
  /** Boolean to make analytical calculations (if possible). */
  const Index& Analytical() const { return manalytical; }
  void Analytical( const Index& m ) { manalytical = m; }
  /** Size of perturbation used for perturbation calculations. */
  const Numeric& Perturbation() const { return mperturbation; }
  void Perturbation( const Numeric& p ) { mperturbation = p; }
  /** Grids. Definition grids for the jacobian, eg. p, lat and lon. */
  const ArrayOfVector& Grids() const { return mgrids; }
  void Grids( const ArrayOfVector& g ) { mgrids = g; }

  Index nelem() const {
      Index i = 1;
      for (Index j = 0; j < mgrids.nelem(); ++j) {
          i *= mgrids[j].nelem();
      }
      return i;
  }

  /** QuantumIdentifier as necessary for matching line specific parameters to jacobian grid */
  const QuantumIdentifier& QuantumIdentity() const { return mquantumidentifier; }
  void QuantumIdentity( const QuantumIdentifier& qi ) { mquantumidentifier = qi; }
  
  /** Do integrations? */
  void IntegrationOn() { mintegration_flag = true; }
  void IntegrationOff() { mintegration_flag = false; }
  const bool& Integration() const { return mintegration_flag; }

  /** Transformation **/
  void SetTransformationFunc(const String& s)   {transformation_func = s;}
  void SetTransformationMatrix(const Matrix& A) {transformation_matrix = A;}
  void SetOffsetVector(const Vector& b)         {offset_vector = b;}
  bool HasAffine()                     const    {return !transformation_matrix.empty();}
  const String& TransformationFunc()   const    {return transformation_func;}
  const Matrix& TransformationMatrix() const    {return transformation_matrix;}
  const Vector& OffsetVector()         const    {return offset_vector;}

  /** Constraints **/
  bool HasConstraints() const {return constraints.nelem() > 0;}
  void AddConstraint(const String& s, Numeric b)
  {
      constraints.push_back(s);
      boundaries.push_back(b);
  }
  const ArrayOfString&  GetConstraints() const {return constraints;}
  const ArrayOfNumeric& GetBoundaries()  const {return boundaries;}


private:

  String mmaintag;
  String msubtag;
  String msubsubtag;
  String mmode;
  Index manalytical;
  Numeric mperturbation;
  ArrayOfVector mgrids;
  QuantumIdentifier mquantumidentifier;
  bool mintegration_flag;

  String transformation_func;
  
  Matrix transformation_matrix;
  Vector offset_vector;

  ArrayOfString  constraints;
  ArrayOfNumeric boundaries;

};



/** Output operator for RetrievalQuantity.

    \author Mattias Ekstrom */
ostream& operator << (ostream& os, const RetrievalQuantity& ot);

typedef Array<RetrievalQuantity> ArrayOfRetrievalQuantity;



// A macro to loop analytical jacobian quantities
#define FOR_ANALYTICAL_JACOBIANS_DO(what_to_do) \
  for( Index iq=0; iq<jacobian_quantities.nelem(); iq++ ) \
    { \
      if( jacobian_quantities[iq].Analytical() ) \
        { what_to_do } \
    } 


//======================================================================
//             Index ranges and transformation functions
//======================================================================

void jac_ranges_indices(
          ArrayOfArrayOfIndex&      jis,  
          bool&                     any_affine,
    const ArrayOfRetrievalQuantity& jqs,
    const bool&                     before_affine = false );

void transform_jacobian(
    Matrix&                           jacobian,
    const Vector                      x,
    const ArrayOfRetrievalQuantity&   jqs );

void transform_x(
    Vector&                           x,
    const ArrayOfRetrievalQuantity&   jqs );

void transform_x_back(
    Vector&                           x_t,
    const ArrayOfRetrievalQuantity&   jqs );




//======================================================================
//             Functions related to calculation of Jacobian
//======================================================================

void calc_nd_field(       Tensor3View& nd,
                    const VectorView&  p,
                    const Tensor3View& t);

bool check_retrieval_grids(       ArrayOfVector& grids,
                                  ostringstream& os,
                            const Vector&        p_grid,
                            const Vector&        lat_grid,
                            const Vector&        lon_grid,
                            const Vector&        p_retr,
                            const Vector&        lat_retr,
                            const Vector&        lon_retr,
                            const String&        p_retr_name,
                            const String&        lat_retr_name,
                            const String&        lon_retr_name,
                            const Index&         dim);

void diy_from_path_to_rgrids(
         Tensor3View          diy_dx,
   const RetrievalQuantity&   jacobian_quantity,
   ConstTensor3View           diy_dpath,
   const Index&               atmosphere_dim,
   const Ppath&               ppath,
   ConstVectorView            ppath_p );

void get_perturbation_gridpos(      ArrayOfGridPos& gp,
                              const Vector&         atm_grid,
                              const Vector&         jac_grid,
                              const bool&           is_pressure);

void get_perturbation_limit(       ArrayOfIndex& limit,
                             const Vector&       pert_grid,
                             const Vector&       atm_limit);

void get_perturbation_range(       Range& range,
                             const Index& index,
                             const Index& length);

void get_pointers_for_analytical_jacobians( 
         ArrayOfIndex&               abs_species_i, 
         ArrayOfIndex&               scat_species_i, 
         ArrayOfIndex&               is_t,
         ArrayOfIndex&               wind_i,
         ArrayOfIndex&               magfield_i,
         ArrayOfIndex&               integrate_i,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const ArrayOfString&              scat_species );

void jacobian_type_extrapol( ArrayOfGridPos&   gp );

void perturbation_field_1d(       VectorView      field,
                            const ArrayOfGridPos& p_gp,
                            const Index&          p_pert_n,
                            const Range&          p_range,
                            const Numeric&        size,
                            const Index&          method);
                                
void perturbation_field_2d(       MatrixView      field,
                            const ArrayOfGridPos& p_gp,
                            const ArrayOfGridPos& lat_gp,
                            const Index&          p_pert_n,
                            const Index&          lat_pert_n,
                            const Range&          p_range,
                            const Range&          lat_range,
                            const Numeric&        size,
                            const Index&          method);
                                
void perturbation_field_3d(       Tensor3View     field,
                            const ArrayOfGridPos& p_gp,
                            const ArrayOfGridPos& lat_gp,
                            const ArrayOfGridPos& lon_gp,
                            const Index&          p_pert_n,
                            const Index&          lat_pert_n,
                            const Index&          lon_pert_n,
                            const Range&          p_range,
                            const Range&          lat_range,
                            const Range&          lon_range,
                            const Numeric&        size,
                            const Index&          method);

void polynomial_basis_func(
        Vector&   b,
  const Vector&   x,
  const Index&    poly_coeff );


//! Calculate baseline fit
/**
 * Computes the baseline fit from a given state vector.
 *
 * Given a retrieval quantitiy which is either a polynomial or a sine baseline fit
 * this function computes the baseline offset in y_baseline.
 *
 *  \param y_baseline (output) The computed baseline offset. Computed baseline offset are
 *  accumulated into this vector, so it must be initialized externally!
 *  \param x State vector consisten with given retrieval quantity
 *  \param mblock_index The index of the measurement block.
 *  \param sensor_response Must be consistent with size of y_baseline.
 *  \param sensor_response_pol_grid Must be consistent with size of y_baseline.
 *  \param sensor_response_f_grid Must be consistent with size of y_baseline.
 *  \param sensor_dlos_grid Must be consistent with size of y_baseline.
 *  \param rq The poly- or sinefit retrieval quantity
 *  \param rq_index The index of the retrieval quantity
 *  \param jacobian_indices
 */
void calcBaselineFit(
    Vector&                    y_baseline,
    const Vector&                    x,
    const Index&                     mblock_index,
    const Sparse&                    sensor_response,
    const ArrayOfIndex&              sensor_response_pol_grid,
    const Vector&                    sensor_response_f_grid,
    const Matrix&                    sensor_response_dlos_grid,
    const RetrievalQuantity&         rq,
    const Index                      rq_index,
    const ArrayOfArrayOfIndex&       jacobian_indices);

void vmrunitscf(  
        Numeric&   x, 
  const String&    unit, 
  const Numeric&   vmr,
  const Numeric&   p,
  const Numeric&   t );                                

void dxdvmrscf(  
        Numeric&   x, 
  const String&    unit, 
  const Numeric&   vmr,
  const Numeric&   p,
  const Numeric&   t );


// Enum for knowing what Jacobian scheme is in-play in the m_rte.cc methods.
enum {
    JAC_IS_NONE=0,             // Setting to nil means that (bool)0 and (bool)N still works.
    JAC_IS_T_SEMI_ANALYTIC,
    JAC_IS_T_FROM_PROPMAT,
    JAC_IS_WIND_U_SEMI_ANALYTIC,
    JAC_IS_WIND_V_SEMI_ANALYTIC,
    JAC_IS_WIND_W_SEMI_ANALYTIC,
    JAC_IS_WIND_U_FROM_PROPMAT,
    JAC_IS_WIND_V_FROM_PROPMAT,
    JAC_IS_WIND_W_FROM_PROPMAT,
    JAC_IS_WIND_ABS_FROM_PROPMAT,
    JAC_IS_MAG_U_SEMI_ANALYTIC,
    JAC_IS_MAG_V_SEMI_ANALYTIC,
    JAC_IS_MAG_W_SEMI_ANALYTIC,
    JAC_IS_MAG_V_FROM_PROPMAT,
    JAC_IS_MAG_U_FROM_PROPMAT,
    JAC_IS_MAG_W_FROM_PROPMAT,
    JAC_IS_FLUX,
    JAC_IS_INTEGRATION,
    JAC_IS_OTHER
};

void get_diydx( VectorView diydx_this,
                VectorView diydx_next,
                Index& extmat_case,
                ConstVectorView iy,
                ConstVectorView sibi,
                ConstVectorView nlte_this,
                ConstVectorView nlte_next,
                ConstVectorView dnltedx_this,
                ConstVectorView dnltedx_next,
                ConstMatrixView K_this,
                ConstMatrixView K_next,
                ConstMatrixView dKdx_this,
                ConstMatrixView dKdx_next,
                ConstMatrixView T_this,
                ConstMatrixView dTdx_this,
                ConstMatrixView dTdx_next,
                ConstMatrixView PiT_this,
                ConstMatrixView PiT_next,
                const Numeric& temperature_this,
                const Numeric& temperature_next,
                const Numeric& dt,
                const Numeric& dBdx_this,
                const Numeric& dBdx_next,
                const Numeric& r,
                const Index& stokes_dim,
                const bool& do_Bsource,
                const bool& do_HSE,
                const bool& do_nonlte );

void get_diydx(VectorView diy1,
               VectorView diy2,
               ConstMatrixView ImT,
               ConstMatrixView cumulative_transmission,
               ConstMatrixView dT1,
               ConstMatrixView dT2,
               ConstVectorView iYmJ,
               ConstVectorView dJ1,
               ConstVectorView dJ2,
               const Index stokes_dim);

#endif // jacobian_h

