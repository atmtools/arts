/* Copyright (C) 2004-2012 Mattias Ekstrom  <ekstrom@rss.chalmers.se>

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

/*!
  \file   jacobian.cc
  \author Mattias Ekstrom <ekstrom@rss.chalmers.se>
  \date   2004-09-14

  \brief  Routines for setting up the jacobian.
*/

#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "jacobian.h"
#include "special_interp.h"
#include "physics_funcs.h"
#include "lin_alg.h"
#include "rte.h"

extern const Numeric NAT_LOG_TEN;

extern const String  ABSSPECIES_MAINTAG;
extern const String  SCATSPECIES_MAINTAG;
extern const String  TEMPERATURE_MAINTAG;
extern const String  WIND_MAINTAG;
extern const String  MAGFIELD_MAINTAG;
extern const String  FLUX_MAINTAG;
extern const String  PROPMAT_SUBSUBTAG;


ostream& operator << (ostream& os, const RetrievalQuantity& ot)
{
  return os << "\n       Main tag = " << ot.MainTag() 
            << "\n       Sub  tag = " << ot.Subtag()
            << "\n           Mode = " << ot.Mode()
            << "\n     Analytical = " << ot.Analytical();
}


//! jac_ranges_indices
/*!
    Determines the index range inside x and the Jacobian for each retrieval quantity

    The ranges are given as an ArrayOfArrayOfIndex, where outermots dimension
    corresponds to retrieval quantity. The inner dimension has throughout size
    2, where element 0 is the first index and element 1 is the last index of
    the range.

    \param   jis             Out: Indices, as described above
    \param   any_affine      Out: True if at least one  quantity has affine
                             transformation. 
    \param   jqs             The WSV jacobian_quantities
    \param   before_affine   Set to true to get indices without any affine
                             transformation. Default is false.

    \author Simon Pfreundschuh and Patrick Eriksson 
    \date   2017-12-30
*/
void jac_ranges_indices(
          ArrayOfArrayOfIndex&      jis,
          bool&                     any_affine,
    const ArrayOfRetrievalQuantity& jqs,
    const bool&                     before_affine )
{
  jis.resize( jqs.nelem() );

  any_affine = false;
  
  // Indices before affine transformation
  if( before_affine )
    {
      for( Index i = 0; i < jqs.nelem(); ++i )
        {
          jis[i] = ArrayOfIndex(2);
          if (i > 0) {
            jis[i][0] = jis[i-1][1] + 1;
          } else {
            jis[i][0] = 0;
          }
          const RetrievalQuantity &jq = jqs[i];
          jis[i][1] = jis[i][0] + jq.nelem() - 1;
          if (jq.HasAffine()) {
            any_affine = true;
          }          
        }
    }
  
  // After affine transformation
  else
    {
      for( Index i = 0; i < jqs.nelem(); ++i )
        {
          jis[i] = ArrayOfIndex(2);
          if (i > 0) {
            jis[i][0] = jis[i-1][1] + 1;
          } else {
            jis[i][0] = 0;
          }
          const RetrievalQuantity &jq = jqs[i];
          if (jq.HasAffine()) {
            jis[i][1] = jis[i][0] + jq.TransformationMatrix().ncols() - 1;
            any_affine = true;
          } else {
            jis[i][1] = jis[i][0] + jq.nelem() - 1;
          }
        }
    }
}



//! Handles transformations of the Jacobian
/**
 * Applies both functional and affine transformations.
 *
 *  \param jacobian As the WSV jacobian
 *  \param jqs As the WSV jacobian_quantities
 *
 *  \author Simon Pfreundschuh and Patrick Eriksson 
 *  \date   2017-12-30
 */
void transform_jacobian(
    Matrix&                           jacobian,
    const Vector                      x,
    const ArrayOfRetrievalQuantity&   jqs )
{
  // Range indices without affine trans
  ArrayOfArrayOfIndex jis;
  bool any_affine;
  //
  jac_ranges_indices( jis, any_affine, jqs, true );

  // Apply functional transformations
  for (Index i = 0; i < jqs.nelem(); ++i) {
    const RetrievalQuantity &jq = jqs[i];
    const String tfun = jq.TransformationFunc();
    // Remember to add new functions also to transform_jacobian and transform_x_back
    if (tfun == "") {
      // Nothing to do
    }
    else if (tfun == "log") {
      for (Index c = jis[i][0]; c <= jis[i][1]; ++c) {
        jacobian(joker,c) *= exp( x[c] );
      }
    }
    else if (tfun == "log10") {
      for (Index c = jis[i][0]; c <= jis[i][1]; ++c) {
        jacobian(joker,c) *= NAT_LOG_TEN * pow(10.0,x[c]);
      }
    }
    else if (tfun == "atanh") {
      const Numeric xmax = jq.TFuncParameter();;
      for (Index c = jis[i][0]; c <= jis[i][1]; ++c) {
        jacobian(joker,c) *= xmax * 2 / pow(exp(-x[c])+exp(x[c]),2.0);
      }
    }
    else{
      assert(0);
    }
  }
  
  // Apply affine transformations
  if( any_affine )
    {
      ArrayOfArrayOfIndex jis_t;
      jac_ranges_indices( jis_t, any_affine, jqs );

      Matrix jacobian_t(jacobian.nrows(), jis_t.back()[1] + 1);

      for (Index i = 0; i < jqs.nelem(); ++i) {
        const RetrievalQuantity &jq = jqs[i];
        Index col_start  = jis[i][0];
        Index col_extent = jis[i][1] - jis[i][0] + 1;
        Range col_range(col_start, col_extent);
        Index col_start_t  = jis_t[i][0];
        Index col_extent_t = jis_t[i][1] - jis_t[i][0] + 1;
        Range col_range_t(col_start_t, col_extent_t);
        if (jq.HasAffine()) {
            mult(jacobian_t(joker, col_range_t),
                 jacobian(joker, col_range),
                 jq.TransformationMatrix());
        } else {
            jacobian_t(joker, col_range_t) = jacobian(joker, col_range);
        }
      }
      swap(jacobian_t, jacobian);
    }
}


//! Handles transformations of the state vector
/**
 * Applies both functional and affine transformations.
 *
 *  \param x As the WSV x
 *  \param jqs As the WSV jacobian_quantities
 *
 *  \author Simon Pfreundschuh and Patrick Eriksson 
 *  \date   2017-12-30
 */
void transform_x(
    Vector&                           x,
    const ArrayOfRetrievalQuantity&   jqs )
{
  // Range indices without affine trans
  ArrayOfArrayOfIndex jis;
  bool any_affine;
  //
  jac_ranges_indices( jis, any_affine, jqs, true );
  
  // Apply functional transformations
  for (Index i = 0; i < jqs.nelem(); ++i) {
    const RetrievalQuantity &jq = jqs[i];
    const String tfun = jq.TransformationFunc();
    // Remember to add new functions also to transform_jacobian and transform_x_back
    if (tfun == "") {
      // Nothing to do
    }
    else if (tfun == "log") {
      for (Index r = jis[i][0]; r <= jis[i][1]; ++r) {
        if( x[r] <= 0 )
          {
            ostringstream os;
            os << "log-transformation selected for retrieval quantity with\n"
               << "index " << i << " (0-based), but at least one value <= 0\n"
               << "found for this quantity. This is not allowed.";
            throw std::runtime_error(os.str());
          }
        x[r] = log( x[r] );
      }
    }
    else if (tfun == "log10") {
      for (Index r = jis[i][0]; r <= jis[i][1]; ++r) {
        if( x[r] <= 0 )
          {
            ostringstream os;
            os << "log10-transformation selected for retrieval quantity with\n"
               << "index " << i << " (0-based), but at least one value <= 0\n"
               << "found for this quantity. This is not allowed.";
            throw std::runtime_error(os.str());
          }
        x[r] = log10( x[r] );
      }
    }
    else if (tfun == "atanh") {
      const Numeric xmax = jq.TFuncParameter();;
      for (Index r = jis[i][0]; r <= jis[i][1]; ++r) {
        if( x[r] <= 0 )
          {
            ostringstream os;
            os << "atanh-transformation selected for retrieval quantity with\n"
               << "index " << i << " (0-based), but at least one value <= 0\n"
               << "found for this quantity. This is not allowed.";
            throw std::runtime_error(os.str());
          }
        if( x[r] >= xmax )
          {
            ostringstream os;
            os << "atanh-transformation selected for retrieval quantity with\n"
               << "index " << i << " (0-based), but at least one value is\n"
               << ">= tfunc_parameter. This is not allowed.";
            throw std::runtime_error(os.str());
          }
        x[r] = atanh( 2*x[r]/xmax - 1 );
      }
    }
    else{
      assert(0);
    }
  }
  
  // Apply affine transformations
  if( any_affine )
    {
      ArrayOfArrayOfIndex jis_t;
      jac_ranges_indices( jis_t, any_affine, jqs );

      Vector x_t(jis_t.back()[1] + 1);

      for (Index i = 0; i < jqs.nelem(); ++i) {
        const RetrievalQuantity &jq = jqs[i];
        Index col_start  = jis[i][0];
        Index col_extent = jis[i][1] - jis[i][0] + 1;
        Range col_range(col_start, col_extent);
        Index col_start_t  = jis_t[i][0];
        Index col_extent_t = jis_t[i][1] - jis_t[i][0] + 1;
        Range col_range_t(col_start_t, col_extent_t);
        if (jq.HasAffine()) {
            Vector t(x[col_range]);
            t -= jq.OffsetVector();
            mult(x_t[col_range_t], transpose(jq.TransformationMatrix()), t);
        } else {
            x_t[col_range_t] = x[col_range];
        }
      }
      swap(x, x_t);
    }
}


//! Handles back-transformations of the state vector
/**
 * Applies both functional and affine transformations.
 *
 *  \param x As the WSV x
 *  \param jqs As the WSV jacobian_quantities
 *
 *  \author Simon Pfreundschuh and Patrick Eriksson 
 *  \date   2017-12-30
 */
void transform_x_back(
    Vector&                           x_t,
    const ArrayOfRetrievalQuantity&   jqs )
{
  // Range indices without affine trans
  ArrayOfArrayOfIndex jis;
  bool any_affine;
  //
  jac_ranges_indices( jis, any_affine, jqs, true );

  // Revert affine transformations
  // Apply affine transformations
  if( any_affine )
    {
      ArrayOfArrayOfIndex jis_t;
      jac_ranges_indices( jis_t, any_affine, jqs );

      Vector x(jis.back()[1] + 1);

      for (Index i = 0; i < jqs.nelem(); ++i) {
        const RetrievalQuantity &jq = jqs[i];
        Index col_start  = jis[i][0];
        Index col_extent = jis[i][1] - jis[i][0] + 1;
        Range col_range(col_start, col_extent);
        Index col_start_t  = jis_t[i][0];
        Index col_extent_t = jis_t[i][1] - jis_t[i][0] + 1;
        Range col_range_t(col_start_t, col_extent_t);
        if (jq.HasAffine()) {
            mult(x[col_range], jq.TransformationMatrix(), x_t[col_range_t]);
            x[col_range] += jq.OffsetVector();
        } else {
            x[col_range] = x_t[col_range_t];
        }
      }
      swap(x_t, x);
    }

  // Revert functional transformations
  for (Index i = 0; i < jqs.nelem(); ++i) {
    const RetrievalQuantity &jq = jqs[i];
    const String tfun = jq.TransformationFunc();
    // Remember to add new functions also to transform_jacobian and transform_x_back
    if (tfun == "") {
      // Nothing to do
    }
    else if (tfun == "log") {
      for (Index r = jis[i][0]; r <= jis[i][1]; ++r) {
        x_t[r] = exp( x_t[r] );
      }
    }
    else if (tfun == "log10") {
      for (Index r = jis[i][0]; r <= jis[i][1]; ++r) {
        x_t[r] = pow( 10.0, x_t[r] );
      }
    }
    else if (tfun == "atanh") {
      const Numeric xmax = jq.TFuncParameter();
      for (Index r = jis[i][0]; r <= jis[i][1]; ++r) {
        x_t[r] = xmax/2 * ( 1 + tanh( x_t[r] ) );
      }
    }
    else{
      assert(0);
    }
  }  
}


/*===========================================================================
  === Help sub-functions to handle analytical jacobians (in alphabetical order)
  ===========================================================================*/

// Small help function, to make the code below cleaner
void from_dpath_to_dx(
        MatrixView   diy_dx,
   ConstMatrixView   diy_dq,
   const Numeric&    w )
{
  for( Index irow=0; irow<diy_dx.nrows(); irow++ )
    { 
      for( Index icol=0; icol<diy_dx.ncols(); icol++ )
        { diy_dx(irow,icol) += w * diy_dq(irow,icol); }
    }
}

//! diy_from_path_to_rgrids
/*!
    Maps jacobian data for points along the propagation path, to
    jacobian retrieval grid data.

    \param   diy_dx              Out: One lement of the WSV *diy_dx*.
    \param   jacobian_quantity   One element of of the WSV *jacobian_quantities*.
    \param   diy_dpath           Jacobians along the propagation path.
    \param   atmosphere_dim      As the WSV.
    \param   ppath               As the WSV.
    \param   ppath_p             The pressure at each ppath point.

    \author Patrick Eriksson 
    \date   2009-10-08
*/
void diy_from_path_to_rgrids(
         Tensor3View          diy_dx,
   const RetrievalQuantity&   jacobian_quantity,
   ConstTensor3View           diy_dpath,
   const Index&               atmosphere_dim,
   const Ppath&               ppath,
   ConstVectorView            ppath_p )
{
  // If this is an integration target then diy_dx is just the sum of all in diy_dpath
  if( jacobian_quantity.Integration() )
    {
      diy_dx(0, joker, joker) = diy_dpath(0, joker, joker);
      for(Index i = 1; i < diy_dpath.npages(); i++)
        diy_dx(0, joker, joker) += diy_dpath(i, joker, joker);
      return;
    }

  assert( jacobian_quantity.Grids().nelem() == atmosphere_dim );
   
  // We want here an extrapolation to infinity -> 
  //                                        extremly high extrapolation factor
  const Numeric   extpolfac = 1.0e99;

  if( ppath.np > 1 )  // Otherwise nothing to do here
    {
      // Pressure
      Index            nr1 = jacobian_quantity.Grids()[0].nelem();
      ArrayOfGridPos   gp_p(ppath.np);
      if( nr1 > 1 )
        {
          p2gridpos( gp_p, jacobian_quantity.Grids()[0], ppath_p, extpolfac );
          jacobian_type_extrapol( gp_p );
        }
      else
        { gp4length1grid( gp_p ); }        

      // Latitude
      Index            nr2 = 1;
      ArrayOfGridPos   gp_lat;
      if( atmosphere_dim > 1 )
        {          
          gp_lat.resize(ppath.np);
          nr2    = jacobian_quantity.Grids()[1].nelem();
          if( nr2 > 1 )
            {
              gridpos( gp_lat, jacobian_quantity.Grids()[1], 
                       ppath.pos(joker,1), extpolfac );
              jacobian_type_extrapol( gp_lat );
            }
          else
            { gp4length1grid( gp_lat ); }
        }

      // Longitude
      ArrayOfGridPos   gp_lon;
      if( atmosphere_dim > 2 )
        {
          gp_lon.resize(ppath.np);
          if( jacobian_quantity.Grids()[2].nelem() > 1 )
            {          
              gridpos( gp_lon, jacobian_quantity.Grids()[2], 
                       ppath.pos(joker,2), extpolfac );
              jacobian_type_extrapol( gp_lon );
            }
          else
            { gp4length1grid( gp_lon ); }
        }
      
      //- 1D
      if( atmosphere_dim == 1 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              if( gp_p[ip].fd[1] > 0 )
                {
                  from_dpath_to_dx( diy_dx(gp_p[ip].idx,joker,joker),
                                    diy_dpath(ip,joker,joker), gp_p[ip].fd[1] );
                }
              if( gp_p[ip].fd[0] > 0 )
                {
                  from_dpath_to_dx( diy_dx(gp_p[ip].idx+1,joker,joker),
                                    diy_dpath(ip,joker,joker), gp_p[ip].fd[0] );
                }
            }
        }

      //- 2D
      else if( atmosphere_dim == 2 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              Index ix = nr1*gp_lat[ip].idx + gp_p[ip].idx;
              // Low lat, low p
              if( gp_lat[ip].fd[1]>0 && gp_p[ip].fd[1]>0 )
                from_dpath_to_dx( diy_dx(ix,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                                  gp_lat[ip].fd[1]*gp_p[ip].fd[1] );
              // Low lat, high p
              if( gp_lat[ip].fd[1]>0 && gp_p[ip].fd[0]>0 )
                from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                                  gp_lat[ip].fd[1]*gp_p[ip].fd[0] );
              // High lat, low p
              if( gp_lat[ip].fd[0]>0 && gp_p[ip].fd[1]>0 )
                from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                                  gp_lat[ip].fd[0]*gp_p[ip].fd[1] );
              // High lat, high p
              if( gp_lat[ip].fd[0]>0 && gp_p[ip].fd[0]>0 )
                from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                                  gp_lat[ip].fd[0]*gp_p[ip].fd[0] );
            }
        }

      //- 3D
      else if( atmosphere_dim == 3 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              Index   ix = nr2*nr1*gp_lon[ip].idx +
                           nr1*gp_lat[ip].idx + gp_p[ip].idx;
              // Low lon, low lat, low p
              if( gp_lon[ip].fd[1]>0 && gp_lat[ip].fd[1]>0 && gp_p[ip].fd[1]>0 )
                from_dpath_to_dx( diy_dx(ix,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[1]*gp_p[ip].fd[1]);
              // Low lon, low lat, high p
              if( gp_lon[ip].fd[1]>0 && gp_lat[ip].fd[1]>0 && gp_p[ip].fd[0]>0 )
                from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[1]*gp_p[ip].fd[0]);
              // Low lon, high lat, low p
              if( gp_lon[ip].fd[1]>0 && gp_lat[ip].fd[0]>0 && gp_p[ip].fd[1]>0 )
                from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[0]*gp_p[ip].fd[1]);
              // Low lon, high lat, high p
              if( gp_lon[ip].fd[1]>0 && gp_lat[ip].fd[0]>0 && gp_p[ip].fd[0]>0 )
                from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[0]*gp_p[ip].fd[0]);

              // Increase *ix* (to be valid for high lon level)
              ix += nr2*nr1;

              // High lon, low lat, low p
              if( gp_lon[ip].fd[0]>0 && gp_lat[ip].fd[1]>0 && gp_p[ip].fd[1]>0 )
                from_dpath_to_dx( diy_dx(ix,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[1]*gp_p[ip].fd[1]);
              // High lon, low lat, high p
              if( gp_lon[ip].fd[0]>0 && gp_lat[ip].fd[1]>0 && gp_p[ip].fd[0]>0 )
                from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[1]*gp_p[ip].fd[0]);
              // High lon, high lat, low p
              if( gp_lon[ip].fd[0]>0 && gp_lat[ip].fd[0]>0 && gp_p[ip].fd[1]>0 )
                from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[0]*gp_p[ip].fd[1]);
              // High lon, high lat, high p
              if( gp_lon[ip].fd[0]>0 && gp_lat[ip].fd[0]>0 && gp_p[ip].fd[0]>0 )
                from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[0]*gp_p[ip].fd[0]);
            }
        }
    }
}



//! diy_from_pos_to_rgrids
/*!
    Maps jacobian data for a surface position, to jacobian retrieval grid data.

    \param   diy_dx              Out: One lement of the WSV *diy_dx*.
    \param   jacobian_quantity   One element of of the WSV *jacobian_quantities*.
    \param   diy_dpos            Jacobian for the position itself.
    \param   atmosphere_dim      As the WSV.
    \param   rtp_pos             As the WSV.

    \author Patrick Eriksson 
    \date   2018-04-10
*/
void diy_from_pos_to_rgrids(
         Tensor3View          diy_dx,
   const RetrievalQuantity&   jacobian_quantity,
   ConstMatrixView            diy_dpos,
   const Index&               atmosphere_dim,
   ConstVectorView            rtp_pos )
{
  assert( jacobian_quantity.Grids().nelem() == max(atmosphere_dim-1,Index(1)) );
  assert( rtp_pos.nelem() == atmosphere_dim );
  
  // We want here an extrapolation to infinity -> 
  //                                        extremly high extrapolation factor
  const Numeric   extpolfac = 1.0e99;

  // Handle 1D separately
  if( atmosphere_dim == 1 )
    {
      diy_dx(0,joker,joker) = diy_dpos;
      return;
    }

  // Latitude
  Index            nr1 = 1;
  ArrayOfGridPos   gp_lat;
  {          
    gp_lat.resize(1);
    nr1 = jacobian_quantity.Grids()[0].nelem();
    if( nr1 > 1 )
      {
        gridpos( gp_lat, jacobian_quantity.Grids()[0], 
                 Vector(1,rtp_pos[1]), extpolfac );
        jacobian_type_extrapol( gp_lat );
      }
    else
      { gp4length1grid( gp_lat ); }
  }

  // Longitude
  ArrayOfGridPos   gp_lon;
  if( atmosphere_dim > 2 )
    {
      gp_lon.resize(1);
      if( jacobian_quantity.Grids()[1].nelem() > 1 )
        {          
          gridpos( gp_lon, jacobian_quantity.Grids()[1], 
                   Vector(1,rtp_pos[2]), extpolfac );
          jacobian_type_extrapol( gp_lon );
        }
      else
        { gp4length1grid( gp_lon ); }
    }

  //- 2D
  if( atmosphere_dim == 2 )
    {
      if( gp_lat[0].fd[1] > 0 )
        {
          from_dpath_to_dx( diy_dx(gp_lat[0].idx,joker,joker),
                            diy_dpos(joker,joker), gp_lat[0].fd[1] );
        }
      if( gp_lat[0].fd[0] > 0 )
        {
          from_dpath_to_dx( diy_dx(gp_lat[0].idx+1,joker,joker),
                            diy_dpos(joker,joker), gp_lat[0].fd[0] );
        }
    }
  //- 3D
  else 
    {
      Index ix = nr1*gp_lon[0].idx + gp_lat[0].idx;
      // Low lon, low lat
      if( gp_lon[0].fd[1]>0 && gp_lat[0].fd[1]>0 )
        from_dpath_to_dx( diy_dx(ix,joker,joker),
                          diy_dpos(joker,joker), 
                          gp_lon[0].fd[1]*gp_lat[0].fd[1] );
      // Low lon, high lat
      if( gp_lon[0].fd[1]>0 && gp_lat[0].fd[0]>0 )
        from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                          diy_dpos(joker,joker), 
                          gp_lon[0].fd[1]*gp_lat[0].fd[0] );
      // High lon, low lat
      if( gp_lon[0].fd[0]>0 && gp_lat[0].fd[1]>0 )
        from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                          diy_dpos(joker,joker), 
                          gp_lon[0].fd[0]*gp_lat[0].fd[1] );
      // High lon, high lat
      if( gp_lon[0].fd[0]>0 && gp_lat[0].fd[0]>0 )
        from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                          diy_dpos(joker,joker), 
                          gp_lon[0].fd[0]*gp_lat[0].fd[0] );
    }
}



//! Help function for analytical jacobian calculations
/*!
    The function determines which terms in jacobian_quantities that are 
    analytical absorption species and temperature jacobians. 

    *abs_species_i* and *is_t* shall be sized to have the same length
    as *jacobian_quantities*. For analytical absorption species
    jacobians, *abs_species_i* is set to the matching index in
    *abs_species*. Otherwise, to -1. For analytical temperature
    jacobians, *is_t* is set to 1. Otherwise to 0.

    \param   abs_species_i         Out: Matching index in abs_species 
    \param   scat_species_i        Out: Matching index among scattering species 
    \param   is_t                  Out: Flag for: Is a temperature jacobian?
    \param   jacobian_quantities   As the WSV.
    \param   abs_species           As the WSV.


    \author Patrick Eriksson 
    \date   2009-10-07
*/
void get_pointers_for_analytical_jacobians( 
         ArrayOfIndex&               abs_species_i, 
         ArrayOfIndex&               scat_species_i, 
         ArrayOfIndex&               is_t,
         ArrayOfIndex&               wind_i,
         ArrayOfIndex&               magfield_i,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const ArrayOfString&              scat_species )
{
  FOR_ANALYTICAL_JACOBIANS_DO( 
    //
    if( jacobian_quantities[iq].MainTag()   == TEMPERATURE_MAINTAG  &&
        jacobian_quantities[iq].SubSubtag() == PROPMAT_SUBSUBTAG  )
    { is_t[iq] = Index(JacobianType::Temperature); }
    else
      { is_t[iq] = Index(JacobianType::None); }
    //
    if( jacobian_quantities[iq].MainTag() == ABSSPECIES_MAINTAG )
      {
        if( jacobian_quantities[iq].SubSubtag() == PROPMAT_SUBSUBTAG )
          {
            bool test_available=false;
            for(Index ii=0; ii<abs_species.nelem(); ii++)
              {
                if( abs_species[ii][0].Species() ==
                    SpeciesTag(jacobian_quantities[iq].Subtag()).Species() )
                  {
                    test_available=true;
                    abs_species_i[iq]=ii;
                    break;
                  }
              }
            if(!test_available)
              {
                ostringstream os;
                os << "Could not find " << jacobian_quantities[iq].Subtag() <<
                " in species of abs_species.\n";
                throw std::runtime_error(os.str());
              }
          }
        else
          {
            ArrayOfSpeciesTag  atag;
            array_species_tag_from_string( atag, jacobian_quantities[iq].Subtag() );
            abs_species_i[iq] = chk_contains( "abs_species", abs_species, atag );
          }
      }
    else if( jacobian_quantities[iq].MainTag() == PARTICULATES_MAINTAG ||
             jacobian_quantities[iq].MainTag() == ELECTRONS_MAINTAG)
      {
        abs_species_i[iq] = -9999;
      }
    else
      { abs_species_i[iq] = -1; }
    //
    if( jacobian_quantities[iq].MainTag() == SCATSPECIES_MAINTAG )
      {
        scat_species_i[iq] = find_first( scat_species,
                                             jacobian_quantities[iq].Subtag() );
        if( scat_species_i[iq] < 0 )
          {
            ostringstream os;
            os << "Jacobian quantity with index " << iq << " refers to\n"
               << "  " << jacobian_quantities[iq].Subtag()
               << "\nbut this species could not be found in *scat_species*.";
            throw runtime_error(os.str());
          }
      }
    else
      { scat_species_i[iq] = -1; }
    //
    if( jacobian_quantities[iq].MainTag() == WIND_MAINTAG  &&
        jacobian_quantities[iq].SubSubtag() == PROPMAT_SUBSUBTAG  )
      {
        // Map u, v and w to 1, 2 and 3, respectively
        char c = jacobian_quantities[iq].Subtag()[0];
        const Index test = Index( c ) - 116;
        if( test == 1 )
          wind_i[iq] = Index(JacobianType::WindFieldU);
        else if(test == 2 )
          wind_i[iq] = Index(JacobianType::WindFieldV);
        else if(test == 3 )
          wind_i[iq] = Index(JacobianType::WindFieldW);
        else if(test == (Index('s')-116) )
          wind_i[iq] = Index(JacobianType::AbsWind);
      }
    else
      { wind_i[iq] = Index(JacobianType::None); }
    //
    if( jacobian_quantities[iq].MainTag() == MAGFIELD_MAINTAG  &&
        jacobian_quantities[iq].SubSubtag() == PROPMAT_SUBSUBTAG )
      {
        // Map u, v and w to 1, 2 and 3, respectively
        char c = jacobian_quantities[iq].Subtag()[0];
        const Index test = Index( c ) - 116;
        if( test == 1 )
          magfield_i[iq] = Index(JacobianType::MagFieldU);
        else if(test == 2 )
          magfield_i[iq] = Index(JacobianType::MagFieldV);
        else if(test == 3 )
          magfield_i[iq] = Index(JacobianType::MagFieldW);
        else if(test == (Index('s')-116) )
          magfield_i[iq] = Index(JacobianType::AbsMag);
      }
    else
      { magfield_i[iq] = Index(JacobianType::None); }
      
   )
}





/*===========================================================================
  === Other functions, in alphabetical order
  ===========================================================================*/

//! Calculate the number density field
/*! 
   This function returns the number density for each grid point in the 
   Tensor3View.
   
   \param nd  The number density field
   \param p   The pressure grid
   \param t   The temperature field
   
   \author Mattias Ekstrom
   \date   2005-06-03
*/
void calc_nd_field(       Tensor3View& nd,
                    const VectorView&  p,
                    const Tensor3View& t)
{
  assert( nd.npages()==t.npages() );
  assert( nd.nrows()==t.nrows() );               
  assert( nd.ncols()==t.ncols() );
  assert( nd.npages()==p.nelem() );
  
  for (Index p_it=0; p_it<nd.npages(); p_it++)
  {
    for (Index lat_it=0; lat_it<nd.nrows(); lat_it++)
    {
      for (Index lon_it=0; lon_it<nd.ncols(); lon_it++)
      {
        nd(p_it,lat_it,lon_it) = number_density( p[p_it], t(p_it,lat_it,lon_it));
      }
    }
  }
}



//! Check that the retrieval grids are defined for each atmosphere dim
/*!
   Use this version for atmospheric fields.

   This function checks for the given atmosphere dimension that 
     I)  the retrieval grids are defined 
     II) and that they are covered by the corresponding atmospheric grid. 
   If not the return is false and an output string is created to print 
   the error to the user. If the grids are ok they are stored in an array 
   and true is  returned.
   
   \param grids         The array of retrieval grids.
   \param os            The output string stream.
   \param p_grid        The atmospheric pressure grid
   \param lat_grid      The atmospheric latitude grid
   \param lon_grid      The atmospheric longitude grid
   \param p_retr        The pressure retrieval grid.
   \param lat_retr      The latitude retrieval grid.
   \param lon_retr      The longitude retrieval grid.
   \param p_retr_name   The control file name used for the pressure retrieval grid.
   \param lat_retr_name The control file name for the latitude retrieval grid.
   \param lon_retr_name The control file name for the longitude retrieval grid.
   \param dim           The atmosphere dimension
   \return              Boolean for check.
   
   \author Mattias Ekstrom
   \date   2005-05-11
*/ 
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
                            const Index&         dim)
{
  assert( grids.nelem() == dim );

  if ( p_retr.nelem()==0 )
    {
      os << "The grid vector *" << p_retr_name << "* is empty,"
         << " at least one pressure level\n"
         << "should be specified.";
      return false;
    }
  else if( !is_decreasing( p_retr ) )
    {
      os << "The pressure grid vector *" << p_retr_name << "* is not a\n"
         << "strictly decreasing vector, which is required.";
      return false;      
    }
  else if( p_grid.nelem() == 1 and p_grid.nelem() == p_retr.nelem())
  {
    if(p_grid[0] not_eq p_retr[0])
    {
      os << "Mismatching 1-long grids for " << p_retr_name;
      return false;
    }
    
    // Necessary repeat but grids are OK
    grids[0] = p_retr; 
  }
  else if ( log(p_retr[0])> 1.5*log(p_grid[0])-0.5*log(p_grid[1]) || 
            log(p_retr[p_retr.nelem()-1])<1.5*log(p_grid[p_grid.nelem()-1])-
                                          0.5*log(p_grid[p_grid.nelem()-2])) 
    {
      os << "The grid vector *" << p_retr_name << "* is not covered by the\n"
         << "corresponding atmospheric grid.";
      return false;
    }
  else
    {
      // Pressure grid ok, add it to grids
      grids[0]=p_retr;  
    }

  if (dim>=2)
  {
    // If 2D and 3D atmosphere, check latitude grid
    if ( lat_retr.nelem()==0 )
    {
      os << "The grid vector *" << lat_retr_name << "* is empty,"
         << " at least one latitude\n"
         << "should be specified for a 2D/3D atmosphere.";
      return false;
    }
    else if( !is_increasing( lat_retr ) )
    {
      os << "The latitude grid vector *" << lat_retr_name << "* is not a\n"
         << "strictly increasing vector, which is required.";
      return false;      
    }
    else if( lat_grid.nelem() == 1 and lat_grid.nelem() == lat_retr.nelem())
    {
      if(lat_grid[0] not_eq lat_retr[0])
      {
        os << "Mismatching 1-long grids for " << lat_retr_name;
        return false;
      }
      
      // Necessary repeat but grids are OK
      grids[1] = lat_retr; 
    }
    else if ( lat_retr[0]<1.5*lat_grid[0]-0.5*lat_grid[1] || 
              lat_retr[lat_retr.nelem()-1]>1.5*lat_grid[lat_grid.nelem()-1]-
                                           0.5*lat_grid[lat_grid.nelem()-2] )
    {
      os << "The grid vector *" << lat_retr_name << "* is not covered by the\n"
         << "corresponding atmospheric grid.";
      return false;
    }
    else
    {
      // Latitude grid ok, add it to grids
      grids[1]=lat_retr;
    }
    if (dim==3)
    {
      // For 3D atmosphere check longitude grid
      if ( lon_retr.nelem()==0 )
      {
        os << "The grid vector *" << lon_retr_name << "* is empty,"
           << " at least one longitude\n"
           << "should be specified for a 3D atmosphere.";
        return false;
      }
      else if( !is_increasing( lon_retr ) )
      {
      os << "The longitude grid vector *" << lon_retr_name << "* is not a\n"
         << "strictly increasing vector, which is required.";
      return false;      
      }
      else if( lon_grid.nelem() == 1 and lon_grid.nelem() == lon_retr.nelem())
      {
        if(lon_grid[0] not_eq lon_retr[0])
        {
          os << "Mismatching 1-long grids for " << lon_retr_name;
          return false;
        }
        
        // Necessary repeat but grids are OK
        grids[2] = lon_retr; 
      }
      else if ( lon_retr[0]<1.5*lon_grid[0]-0.5*lon_grid[1] || 
                lon_retr[lon_retr.nelem()-1]>1.5*lon_grid[lon_grid.nelem()-1]-
                                             0.5*lon_grid[lon_grid.nelem()-2] )
      {
        os << "The grid vector *" << lon_retr_name << "* is not covered by the\n"
           << "corresponding atmospheric grid.";
        return false;
      }
      else
      {
        // Longitude grid ok, add it to grids      
        grids[2]=lon_retr;
      }
    }
  }
  return true;
}             



//! Check that the retrieval grids are defined for each atmosphere dim
/*!
   Use this version for surface variables

   This function checks for the given atmosphere dimension that 
     I)  the retrieval grids are defined 
     II) and that they are covered by the corresponding atmospheric grid. 
   If not the return is false and an output string is created to print 
   the error to the user. If the grids are ok they are stored in an array 
   and true is  returned.
   
   \param grids         The array of retrieval grids.
   \param os            The output string stream.
   \param lat_grid      The atmospheric latitude grid
   \param lon_grid      The atmospheric longitude grid
   \param lat_retr      The latitude retrieval grid.
   \param lon_retr      The longitude retrieval grid.
   \param lat_retr_name The control file name for the latitude retrieval grid.
   \param lon_retr_name The control file name for the longitude retrieval grid.
   \param dim           The atmosphere dimension
   \return              Boolean for check.
   
   \author Mattias Ekstrom
   \date   2005-05-11
*/ 
bool check_retrieval_grids(       ArrayOfVector& grids,
                                  ostringstream& os,
                            const Vector&        lat_grid,
                            const Vector&        lon_grid,
                            const Vector&        lat_retr,
                            const Vector&        lon_retr,
                            const String&        lat_retr_name,
                            const String&        lon_retr_name,
                            const Index&         dim)
{
  assert( grids.nelem() == max(dim-1,Index(1)) );

  if (dim==1)
  {
    // Here we only need to create a length 1 dummy grid
    grids[0].resize(1);
    grids[0][0] = 0;
  }
  
  if (dim>=2)
  {
    // If 2D and 3D atmosphere, check latitude grid
    if ( lat_retr.nelem()==0 )
    {
      os << "The grid vector *" << lat_retr_name << "* is empty,"
         << " at least one latitude\n"
         << "should be specified for a 2D/3D atmosphere.";
      return false;
    }
    else if( !is_increasing( lat_retr ) )
    {
      os << "The latitude grid vector *" << lat_retr_name << "* is not a\n"
         << "strictly increasing vector, which is required.";
      return false;      
    }
    else if( lat_grid.nelem() == 1 and lat_grid.nelem() == lat_retr.nelem())
    {
      if(lat_grid[0] not_eq lat_retr[0])
      {
        os << "Mismatching 1-long grids for " << lat_retr_name;
        return false;
      }
      
      // Necessary repeat but grids are OK
      grids[0] = lat_retr; 
    }
    else if ( lat_retr[0]<1.5*lat_grid[0]-0.5*lat_grid[1] || 
              lat_retr[lat_retr.nelem()-1]>1.5*lat_grid[lat_grid.nelem()-1]-
                                           0.5*lat_grid[lat_grid.nelem()-2] )
    {
      os << "The grid vector *" << lat_retr_name << "* is not covered by the\n"
         << "corresponding atmospheric grid.";
      return false;
    }
    else
    {
      // Latitude grid ok, add it to grids
      grids[0]=lat_retr;
    }
    
    if (dim==3)
    {
      // For 3D atmosphere check longitude grid
      if ( lon_retr.nelem()==0 )
      {
        os << "The grid vector *" << lon_retr_name << "* is empty,"
           << " at least one longitude\n"
           << "should be specified for a 3D atmosphere.";
        return false;
      }
      else if( !is_increasing( lon_retr ) )
      {
      os << "The longitude grid vector *" << lon_retr_name << "* is not a\n"
         << "strictly increasing vector, which is required.";
      return false;      
      }
      else if( lon_grid.nelem() == 1 and lon_grid.nelem() == lon_retr.nelem())
      {
        if(lon_grid[0] not_eq lon_retr[0])
        {
          os << "Mismatching 1-long grids for " << lon_retr_name;
          return false;
        }
        
        // Necessary repeat but grids are OK
        grids[1] = lon_retr; 
      }
      else if ( lon_retr[0]<1.5*lon_grid[0]-0.5*lon_grid[1] || 
                lon_retr[lon_retr.nelem()-1]>1.5*lon_grid[lon_grid.nelem()-1]-
                                             0.5*lon_grid[lon_grid.nelem()-2] )
      {
        os << "The grid vector *" << lon_retr_name << "* is not covered by the\n"
           << "corresponding atmospheric grid.";
        return false;
      }
      else
      {
        // Longitude grid ok, add it to grids      
        grids[1]=lon_retr;
      }
    }
  }
  return true;
}             





//! Calculate array of GridPos for perturbation interpolation
/*!
   This function constructs a perturbation grid which consists of the
   given retrieval grid with an extra endpoint added at each end.
   These endpoints lies outside the atmospheric grid. This enables the
   interpolation of an perturbation on the perturbation grid to be
   interpolated to the atmospheric grid. For this reason the function
   returns an ArrayOfGridPos. 
   
   If the atmospheric grid is a pressure grid, interpolation is made
   in logarithm of the atmospheric grid.
   
   \param gp          Array of GridPos for interpolation.
   \param atm_grid    Atmospheric grid.
   \param jac_grid    Retrieval grid.
   \param is_pressure True for pressure grid 
   
   \author Mattias Ekstrom
   \date   2005-05-12
*/   
void get_perturbation_gridpos(       ArrayOfGridPos& gp,
                               const Vector&         atm_grid,
                               const Vector&         jac_grid,
                               const bool&           is_pressure)
{
  Index nj = jac_grid.nelem();
  Index na = atm_grid.nelem();
  Vector pert(nj+2);
  
  // Create perturbation grid, with extension outside the atmospheric grid
  if ( is_pressure )
  {
    pert[0] = atm_grid[0]*10.0;
    pert[nj+1] = atm_grid[na-1]*0.1;
  }
  else
  {    
    pert[0] = atm_grid[0]-1.0;
    pert[nj+1] = atm_grid[na-1]+1.0;
  }
  pert[Range(1,nj)] = jac_grid;
  
  // Calculate the gridpos
  gp.resize(na);
  if( is_pressure ){
    p2gridpos( gp, pert, atm_grid);
  }
  else
  { 
    gridpos( gp, pert, atm_grid);
  }
}



//! Get limits for perturbation of a box
/*!
   This is a helper function that calculates the limits where the 
   perturbation should be added to the perturbation grid. 
   This is needed for example by the particle perturbation that only
   should be applied for the cloudbox. The limits are defined as the 
   outermost points lying within or just outside the box limits.
   
   The atmospheric limits should be given in the same unit as the
   perturbation grid. And only the first and last element will be 
   considered as limits. 
   
   Assertions are used to perform checks. The input grids are 
   checked so that the atmospheric limits are containg within 
   the perturbation grid. The limit indices are checked so 
   that they are ordered in increasing order before return.
   
   \param limit     The limit indices in the perturbation grid
   \param pert_grid The perturbation grid
   \param atm_limit The atmospheric limits of the box.

   \author Mattias Ekstrom
   \date   2005-02-25
*/   
void get_perturbation_limit(       ArrayOfIndex& limit,
                             const Vector&       pert_grid,
                             const Vector&       atm_limit)
{
  limit.resize(2);
//   Index np = pert_grid.nelem()-1;
  Index na = atm_limit.nelem()-1;
  
  // If the field is ordered in decreasing order set the
  // increment factor to -1
  Numeric inc = 1;
  if (is_decreasing(pert_grid))
    inc = -1;

  // Check that the pert_grid is encompassing atm_limit
//   assert( inc*pert_grid[0] < inc*atm_limit[0] &&
//           inc*pert_grid[np] > inc*atm_limit[na]);
      
  // Find first limit, check if following value is above lower limit
  limit[0]=0;
  while (inc*pert_grid[limit[0]+1] < inc*atm_limit[0]) 
  {
    limit[0]++;
  }
  
  // Find last limit, check if previous value is below upper limit
  limit[1]=pert_grid.nelem();
  while (inc*pert_grid[limit[1]-1] > inc*atm_limit[na]) 
  {
    limit[1]--;
  }
  // Check that the limits are ok
  assert(limit[1]>limit[0]);
  
}

                             

//! Get range for perturbation 
/*!
   This is a helper function that calculates the range in which the 
   perturbation should be added to the perturbation grid. This is needed
   to handle the edge effects. At the edges we want the perturbation to 
   continue outwards. 
   
   \param range     The range in the perturbation grid.
   \param index     The index of the perturbation in the retrieval grid.
   \param length    The length of retrieval grid
   
   \author Mattias Ekstrom
   \date   2004-10-14
*/   
void get_perturbation_range(       Range& range,
                             const Index& index,
                             const Index& length)
{
  if (index==0)
    range = Range(index,2);
  else if (index==length-1)
    range = Range(index+1,2);
  else 
    range = Range(index+1,1);

}



//! Adopts grid postions to extrapolation used for jacobians
/*!
  The standard interpolation scheme applies a linear extrapolation, while for
  the jacobians the extrapolation can be seen as a "closest" interpolation.
  That is, for points outisde the covered grid, the value at closest end point
  is taken. And here extrapolation to +-Inf is allowed.

  This function modifies grid positions to jacobaina extrapolation approach.
  For efficiency, the input grid positions are not asserted at all, and
  "extrapolation points" are identified simply  by a fd outside [0,1].

  \param[in/out] gp   Array of grid positions.

  \author Patrick Eriksson 
  \date   2015-09-10
*/
void jacobian_type_extrapol( ArrayOfGridPos&   gp )
{
  for( Index i=0; i<gp.nelem(); i++ )
    { 
      if( gp[i].fd[0] < 0 ) 
        {  
          gp[i].fd[0] = 0; 
          gp[i].fd[1] = 1; 
        }
      else if( gp[i].fd[0] > 1 ) 
        {  
          gp[i].fd[0] = 1; 
          gp[i].fd[1] = 0; 
        }
    }
}


//! Calculate the 1D perturbation for a relative perturbation.
/*!
   This is a helper function that interpolated the perturbation field for
   a 1D relative perturbation onto the atmospheric field. 
   
   \param field     The interpolated perturbation field.
   \param p_gp      The GridPos for interpolation.
   \param p_pert_n  The number of perturbations.
   \param p_range   The perturbation range in the perturbation grid.
   \param size      The size of the perturbation.
   \param method    Relative perturbation==0, absolute==1
   
   \author Mattias Ekstrom
   \date   2005-05-11
*/   
void perturbation_field_1d(       VectorView      field,
                            const ArrayOfGridPos& p_gp,
                            const Index&          p_pert_n,
                            const Range&          p_range,
                            const Numeric&        size,
                            const Index&          method)
{
  // Here we only perturb a vector
  Vector pert(field.nelem());
  Matrix itw(p_gp.nelem(),2);
  interpweights(itw,p_gp);
  // For relative pert_field should be 1.0 and for absolute 0.0
  Vector pert_field(p_pert_n,1.0-(Numeric)method);
  pert_field[p_range] += size;
  interp( pert, itw, pert_field, p_gp);
  if (method==0)
  {
    field *= pert;
  }
  else
  {
    field += pert;
  }
}            



//! Calculate the 2D perturbation for a relative perturbation.
/*!
   This is a helper function that interpolated the perturbation field for
   a 2D relative perturbation onto the atmospheric field. 
   
   \param field       The interpolated perturbation field.
   \param p_gp        The GridPos for interpolation in the 1st dim.
   \param lat_gp      The GridPos for interpolation in the 2nd dim.
   \param p_pert_n    The number of perturbations in the 1st dim.
   \param lat_pert_n  The number of perturbations in the 2nd dim.
   \param p_range     The perturbation range in the 1st dim.
   \param lat_range   The perturbation range in the 2nd dim.
   \param size        The size of the perturbation.
   \param method      Relative perturbation==0, absolute==1
   
   \author Mattias Ekstrom
   \date   2005-05-11
*/   
void perturbation_field_2d(       MatrixView      field,
                            const ArrayOfGridPos& p_gp,
                            const ArrayOfGridPos& lat_gp,
                            const Index&          p_pert_n,
                            const Index&          lat_pert_n,
                            const Range&          p_range,
                            const Range&          lat_range,
                            const Numeric&        size,
                            const Index&          method)
{
  // Here we perturb a matrix
  Matrix pert(field.nrows(),field.ncols());
  Tensor3 itw(p_gp.nelem(),lat_gp.nelem(),4);
  interpweights(itw,p_gp,lat_gp);
  // Init pert_field to 1.0 for relative and 0.0 for absolute
  Matrix pert_field(p_pert_n,lat_pert_n,1.0-(Numeric)method);
  pert_field(p_range,lat_range) += size;
  interp( pert, itw, pert_field, p_gp, lat_gp);
  if (method==0)
  {
    field *= pert;
  }
  else
  { 
    field += pert;
  }
}            



//! Calculate the 3D perturbation for a relative perturbation.
/*!
   This is a helper function that interpolated the perturbation field for
   a 3D relative perturbation onto the atmospheric field. 
   
   \param field       The interpolated perturbation field.
   \param p_gp        The GridPos for interpolation in the 1st dim.
   \param lat_gp      The GridPos for interpolation in the 2nd dim.
   \param lon_gp      The GridPos for interpolation in the 3rd dim.
   \param p_pert_n    The number of perturbations in the 1st dim.
   \param lat_pert_n  The number of perturbations in the 2nd dim.
   \param lon_pert_n  The number of perturbations in the 3rd dim.
   \param p_range     The perturbation range in the 1st dim.
   \param lat_range   The perturbation range in the 2nd dim.
   \param lon_range   The perturbation range in the 3rd dim.
   \param size        The size of the perturbation.
   \param method      Set to 0 for relative, and 1 for absolute.
   
   \author Mattias Ekstrom
   \date   2005-05-11
*/   
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
                            const Index&          method)
{
  // Here we need to perturb a tensor3
  Tensor3 pert(field.npages(),field.nrows(),field.ncols());
  Tensor4 itw(p_gp.nelem(),lat_gp.nelem(),lon_gp.nelem(),8);
  interpweights(itw,p_gp,lat_gp,lon_gp);
  // Init pert_field to 1.0 for relative and 0.0 for absolute
  Tensor3 pert_field(p_pert_n,lat_pert_n,lon_pert_n,1.0-(Numeric)method);
  pert_field(p_range,lat_range,lon_range) += size;
  interp( pert, itw, pert_field, p_gp, lat_gp, lon_gp);
  if (method==0)
  {
    field *= pert;
  }
  else
  {
    field += pert;
  }
}            



//! Calculates polynomial basis functions
/*!
   The basis function is b(x) = 1 for poly_coeff = 0. For higher
   coefficients, x^poly_coeff - m, where first the range covered by
   *x* is normalised to [-1,1] and m is selected in such way that
   sum(b) = 0.
   
   \param b            Calculated basis function.
   \param x            The grid over which the fit shall be performed.
   \param poly_coeff   Polynomial coefficient.
   
   \author Patrick Eriksson
   \date   2008-11-07
*/   
void polynomial_basis_func(
        Vector&   b,
  const Vector&   x,
  const Index&    poly_coeff )
{
  const Index l = x.nelem();
  
  assert( l > poly_coeff );

  if( b.nelem() != l )
    b.resize( l );

  if( poly_coeff == 0 )
    { b = 1.0; }
  else
    {
      const Numeric xmin = min( x );
      const Numeric dx = 0.5 * ( max( x ) - xmin );
      //
      for( Index i=0; i<l; i++ )
        {
          b[i] = ( x[i] - xmin ) / dx - 1.0;
          b[i] = pow( b[i], int(poly_coeff) );
        }
      //
      b -= mean( b );
    }  
}

extern const String POLYFIT_MAINTAG;
extern const String SINEFIT_MAINTAG;
extern const Numeric PI;

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
  const ArrayOfArrayOfIndex&       jacobian_indices)
{

    bool is_sine_fit = false;
    if (rq.MainTag() == POLYFIT_MAINTAG) {
        is_sine_fit = false;
    } else if (rq.MainTag() == SINEFIT_MAINTAG) {
        is_sine_fit = true;
    } else {
        throw runtime_error("Retrieval quantity is neither a polynomial or a sine "
                            " baseline fit." );
    }

    // Size and check of sensor_response
    //
    const Index nf     = sensor_response_f_grid.nelem();
    const Index npol   = sensor_response_pol_grid.nelem();
    const Index nlos    = sensor_response_dlos_grid.nrows();

    // Evaluate basis functions for fits.
    Vector w, s, c;
    if (is_sine_fit) {
        s.resize(nf);
        c.resize(nf);
        Numeric period = rq.Grids()[0][0];
        for( Index f=0; f<nf; f++ )
        {
            Numeric a = (sensor_response_f_grid[f]-sensor_response_f_grid[0]) * 
                2 * PI / period;
            s[f] = sin( a );
            c[f] = cos( a );
        }
    } else {
        Numeric poly_coeff = rq.Grids()[0][0];
        polynomial_basis_func( w, sensor_response_f_grid, static_cast<Index>(poly_coeff));
    }

    // Compute baseline
    ArrayOfVector jg   = rq.Grids();
    const Index n1     = jg[1].nelem();
    const Index n2     = jg[2].nelem();
    const Index n3     = jg[3].nelem();
    const Range rowind = get_rowindex_for_mblock( sensor_response, mblock_index );
    const Index row4   = rowind.get_start();
    Index col4   = jacobian_indices[rq_index][0];

    if( n3 > 1 ) {
        col4 += mblock_index*n2*n1;
    }

    for( Index l=0; l<nlos; l++ ) {

        const Index row3 = row4 + l*nf*npol;
        const Index col3 = col4 + l * n1 * (is_sine_fit ? 2 : 1);

        for( Index f=0; f<nf; f++ ) {

            const Index row2 = row3 + f*npol;

            for( Index p=0; p<npol; p++ ) {
                Index col1 = col3;
                if( n1 > 1 ) {
                    col1 += p;
                }
                if (is_sine_fit) {
                    y_baseline[row2+p] += x[col1] * s[f] + x[col1 + 1] * c[f];
                } else {
                    y_baseline[row2+p] += w[f] * x[col1];
                }
            }
        }
    }
}

//! vmrunitscf
/*!
    Scale factor for conversion between gas species units.

    The function finds the factor with which the total absorption of a
    gas species shall be multiplicated to match the selected
    (jacobian) unit. 

    \param   x      Out: scale factor
    \param   unit   Unit selected.
    \param   vmr    VMR value.
    \param   p      Pressure
    \param   t      Temperature.

    \author Patrick Eriksson 
    \date   2009-10-08
*/
void vmrunitscf(  
        Numeric&   x, 
  const String&    unit, 
  const Numeric&   vmr,
  const Numeric&   p,
  const Numeric&   t )
{  
  if( unit == "rel"  ||  unit == "logrel" )
    { x = 1; }
  else if( unit == "vmr" )
    {
        if(vmr==0)
        { 
            x = 0; 
            return;
        }
        x = 1 / vmr;
    }
  else if( unit == "nd" )
    {
        if(vmr==0)
        { 
            x = 0; 
            return;
        }
        x = 1 / ( vmr * number_density( p, t ) ); 
    }
  else
    {
        ostringstream os;
        os << "Allowed options for gas species jacobians are "
        "\"rel\", \"vmr\" and \"nd\".\nYou have selected: "<<unit<<std::endl; 
      throw std::runtime_error( os.str() );
    }
}




//! dxdvmrscf
/*!
    Scale factor for conversion of derivatives with respect to VMR.

    The function finds the factor with which a partial derivative with respect
    to gas species VMR shall be multiplicated to match the selected (jacobian)
    unit. The function was implemented for scaling of *dpropmat_clearsky_dx*
    but it could also be of use in other contexts.

    \param   x      Out: scale factor
    \param   unit   Unit selected.
    \param   vmr    VMR value.
    \param   p      Pressure
    \param   t      Temperature.

    \author Patrick Eriksson 
    \date   2015-12-11
*/
void dxdvmrscf(  
        Numeric&   x, 
  const String&    unit, 
  const Numeric&   vmr,
  const Numeric&   p,
  const Numeric&   t )
{
  if( unit == "rel"  ||  unit == "logrel" )
    { x = vmr; }
  else if( unit == "vmr" )
    { x = 1; }
  else if( unit == "nd" )
    { x = 1 / number_density( p, t ); }
  else
    {
        ostringstream os;
        os << "Allowed options for gas species jacobians are "
        "\"rel\", \"vmr\" and \"nd\".\nYou have selected: "<<unit<<std::endl; 
      throw std::runtime_error( os.str() );
    }
}



/*!
 *   The function helps to calculate the partial derivative of iy with respect
 *   to one input at one pressure.  The formalism here assumes that the radiation
 *   terms are averaged rather than the absorption parameters, thus this can be 
 *   solved per layer rather than for two layers at a time.  Still, the absorption
 *   parameters for the transmission needs to be considered by the two layer derivatives
 * 
 *   FIXME:  Add HSE support
 * 
 *   \param   diy1                      Out: diy for a level encountered the first time
 *   \param   diy2                      Out: diy for a level encountered the second time
 *   \param   ImT                       In: identity matrix minus tranmsission matrix
 *   \param   cumulative_transmission   In: cumulative transmission from level to sensor
 *   \param   dT1                       In: transmission matrix derivative for the first time
 *   \param   dT2                       In: transmission matrix derivative for the second time
 *   \param   iYmJ                      In: incoming radiation to layer minus source of layer
 *   \param   dJ1                       In: derivative of source term emitted for the first time
 *   \param   dJ2                       In: derivative of source term emitted for the second time
 *   \param   stokes_dim                In: essentially the size of the problem
 *   \param   transmission_only         In: remove all computations on source terms, making iYmJ pure incoming radiation
 * 
 *   \author Richard Larsson
 *   \date   2017-09-20
 */
void get_diydx(VectorView diy1,
               VectorView diy2,
               ConstMatrixView ImT,
               ConstMatrixView cumulative_transmission,
               ConstMatrixView dT1,
               ConstMatrixView dT2,
               ConstVectorView iYmJ,
               ConstVectorView dJ1,
               ConstVectorView dJ2,
               const Index stokes_dim,
               const bool transmission_only)
{
  /*
   * Solves 
   * 
   * diy1 = PiT [ dT1 iYmJ + (1-T) dJ1 ],
   * 
   * and
   * 
   * diy2 += PiT [ dT2 iYmJ + (1-T) dJ2 ],
   * 
   * where diy2 is diy1 from a prior layer
   * 
   * FIXME:  Needs HSE
  */
  
  // Computation vectors
  Vector a(stokes_dim), b(stokes_dim);
  
  // The first time a level is involved in a layer
  mult(a, dT1, iYmJ);
  if(not transmission_only)
  {
    mult(b, ImT, dJ1);
    a += b;
  }
  mult(diy1, cumulative_transmission, a);
  
  // The second time a level is involved in a layer
  mult(a, dT2, iYmJ);
  if(not transmission_only)
  {
    mult(b, ImT, dJ2);
    a += b;
  }
  mult(b, cumulative_transmission, a);
  diy2 += b; 
}

