/* Copyright (C) 2012 Stefan Buehler <sbuehler(at)ltu.se>

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
  \file   interpolation_poly.h
  \author Stefan Buehler <sbuehler(at)ltu.se>
  \date   2008-02-04
  
  \brief  Header file for interpolation_poly.cc.
*/

#ifndef interpolation_poly_h
#define interpolation_poly_h

#include "matpackI.h"
#include "interpolation.h"

//! Structure to store a grid position for higher order interpolation.
/*! 
  This serves the same purpose as GridPos for linear
  interpolation. The main difference is that we store not fractional
  distances (fd), but interpolation weights (w). For the linear case:
  w[0] = fd[1]
  w[1] = fd[0]

  Because the w are associated with each interpolation point, they
  work also for higher interpolation order, whereas the concept of
  fractional distance does not.

  There is some confusion between the meaning of fractional distances
  and interpolation weights. In fact the two are almost the same in the 1D
  case! (But the sorting and signs can be different.) 

  The w here correspond exactly to the terms in front of the yi in
  Numerical Recipes, 2nd edition, section 3.1, eq. 3.1.1.

  Only for 2D or higher dimensional interpolation are further
  calculations necessary, to multiply the w for the individual
  dimensions. 

  The size of w must equal the number of points used in the
  interpolation m. (m=2 For linear, m=3 for quadratic, etc..)

  Another important difference is that we store not only the index of
  the first interpolation point in the original grid, but the indices
  of all. We store these explicitly, to allow correct handling of the
  case of cyclic interpolation (in longitude or phi angle) by a future
  specialized gridpos function.
*/
struct GridPosPoly {
  /*! Indices of the interpolation point in the original grid. 
      (Dimension is the number of points in the interpolation, m.)*/
  ArrayOfIndex idx;
  /*! Interpolation weight for each grid point to use.
      (Dimension is the number of points in the interpolation, m.)  */
  Vector w;
};

//! An Array of grid positions.
/*! 
  See \ref GridPosPoly for details.
*/

typedef Array<GridPosPoly> ArrayOfGridPosPoly;

ostream& operator<<(ostream& os, const GridPosPoly& gp);

void gridpos_poly(ArrayOfGridPosPoly& gp,
                  ConstVectorView old_grid,
                  ConstVectorView new_grid,
                  const Index     order,
                  const Numeric&  extpolfac = 0.5);

void gridpos_poly(GridPosPoly& gp,
                  ConstVectorView old_grid,
                  const Numeric&  new_grid,
                  const Index     order,
                  const Numeric&  extpolfac = 0.5);

void gridpos_poly_longitudinal(const String&   error_msg,
                               ArrayOfGridPosPoly& gp,
                               ConstVectorView old_grid,
                               ConstVectorView new_grid,
                               const Index     order,
                               const Numeric&  extpolfac = 0.5);

void gridpos_poly_cyclic_longitudinal(ArrayOfGridPosPoly& gp,
                                      ConstVectorView old_grid,
                                      ConstVectorView new_grid,
                                      const Index     order,
                                      const Numeric&  extpolfac = 0.5);

void gridpos_poly_longitudinal( GridPosPoly&    gp,
                   ConstVectorView old_grid,
                   const Numeric&  new_grid,
                   const Index     order,
                   const Numeric&  extpolfac = 0.5 );

void gridpos_poly_cyclic_longitudinal( GridPosPoly&    gp,
                   ConstVectorView old_grid,
                   const Numeric&  new_grid,
                   const Index     order,
                   const Numeric&  extpolfac = 0.5 );




////////////////////////////////////////////////////////////////////////////
//                      Red Interpolation
////////////////////////////////////////////////////////////////////////////

void interpweights( VectorView itw,
                    const GridPosPoly& tc );

void interpweights( VectorView itw,
                    const GridPosPoly& tr,
                    const GridPosPoly& tc );

void interpweights( VectorView itw,
                    const GridPosPoly& tp,
                    const GridPosPoly& tr,
                    const GridPosPoly& tc );

void interpweights( VectorView itw,
                    const GridPosPoly& tb,
                    const GridPosPoly& tp,
                    const GridPosPoly& tr,
                    const GridPosPoly& tc );

void interpweights( VectorView itw,
                    const GridPosPoly& ts,
                    const GridPosPoly& tb,
                    const GridPosPoly& tp,
                    const GridPosPoly& tr,
                    const GridPosPoly& tc );

void interpweights( VectorView itw,
                    const GridPosPoly& tv,
                    const GridPosPoly& ts,
                    const GridPosPoly& tb,
                    const GridPosPoly& tp,
                    const GridPosPoly& tr,
                    const GridPosPoly& tc );

Numeric interp( ConstVectorView itw,
                ConstVectorView a,    
                const GridPosPoly&  tc );

Numeric interp( ConstVectorView  itw,
                ConstMatrixView  a,    
                const GridPosPoly&   tr,
                const GridPosPoly&   tc );

Numeric interp( ConstVectorView  itw,
                ConstTensor3View a,    
                const GridPosPoly&   tp,
                const GridPosPoly&   tr,
                const GridPosPoly&   tc );

Numeric interp( ConstVectorView  itw,
                ConstTensor4View a,    
                const GridPosPoly&   tb,
                const GridPosPoly&   tp,
                const GridPosPoly&   tr,
                const GridPosPoly&   tc );

Numeric interp( ConstVectorView  itw,
                ConstTensor5View a,    
                const GridPosPoly&   ts,
                const GridPosPoly&   tb,
                const GridPosPoly&   tp,
                const GridPosPoly&   tr,
                const GridPosPoly&   tc );

Numeric interp( ConstVectorView  itw,
                ConstTensor6View a,    
                const GridPosPoly&   tv,
                const GridPosPoly&   ts,
                const GridPosPoly&   tb,
                const GridPosPoly&   tp,
                const GridPosPoly&   tr,
                const GridPosPoly&   tc );





////////////////////////////////////////////////////////////////////////////
//                      Blue interpolation
////////////////////////////////////////////////////////////////////////////

void interpweights( MatrixView itw,
                    const ArrayOfGridPosPoly& cgp );

void interpweights( MatrixView itw,
                    const ArrayOfGridPosPoly& rgp,
                    const ArrayOfGridPosPoly& cgp );

void interpweights( MatrixView itw,
                    const ArrayOfGridPosPoly& pgp,
                    const ArrayOfGridPosPoly& rgp,
                    const ArrayOfGridPosPoly& cgp );

void interpweights( MatrixView itw,
                    const ArrayOfGridPosPoly& bgp,
                    const ArrayOfGridPosPoly& pgp,
                    const ArrayOfGridPosPoly& rgp,
                    const ArrayOfGridPosPoly& cgp );

void interpweights( MatrixView itw,
                    const ArrayOfGridPosPoly& sgp,
                    const ArrayOfGridPosPoly& bgp,
                    const ArrayOfGridPosPoly& pgp,
                    const ArrayOfGridPosPoly& rgp,
                    const ArrayOfGridPosPoly& cgp );

void interpweights( MatrixView itw,
                    const ArrayOfGridPosPoly& vgp,
                    const ArrayOfGridPosPoly& sgp,
                    const ArrayOfGridPosPoly& bgp,
                    const ArrayOfGridPosPoly& pgp,
                    const ArrayOfGridPosPoly& rgp,
                    const ArrayOfGridPosPoly& cgp );

void interp( VectorView            ia,
             ConstMatrixView       itw,
             ConstVectorView       a,    
             const ArrayOfGridPosPoly& cgp);

void interp( VectorView            ia,
             ConstMatrixView       itw,
             ConstMatrixView       a,    
             const ArrayOfGridPosPoly& rgp,
             const ArrayOfGridPosPoly& cgp);

void interp( VectorView            ia,
             ConstMatrixView       itw,
             ConstTensor3View      a,    
             const ArrayOfGridPosPoly& pgp,
             const ArrayOfGridPosPoly& rgp,
             const ArrayOfGridPosPoly& cgp);

void interp( VectorView            ia,
             ConstMatrixView       itw,
             ConstTensor4View      a,    
             const ArrayOfGridPosPoly& bgp,
             const ArrayOfGridPosPoly& pgp,
             const ArrayOfGridPosPoly& rgp,
             const ArrayOfGridPosPoly& cgp);

void interp( VectorView            ia,
             ConstMatrixView       itw,
             ConstTensor5View      a,    
             const ArrayOfGridPosPoly& sgp,
             const ArrayOfGridPosPoly& bgp,
             const ArrayOfGridPosPoly& pgp,
             const ArrayOfGridPosPoly& rgp,
             const ArrayOfGridPosPoly& cgp);

void interp( VectorView            ia,
             ConstMatrixView       itw,
             ConstTensor6View      a,    
             const ArrayOfGridPosPoly& vgp,
             const ArrayOfGridPosPoly& sgp,
             const ArrayOfGridPosPoly& bgp,
             const ArrayOfGridPosPoly& pgp,
             const ArrayOfGridPosPoly& rgp,
             const ArrayOfGridPosPoly& cgp);





////////////////////////////////////////////////////////////////////////////
//                      Green interpolation
////////////////////////////////////////////////////////////////////////////

void interpweights( Tensor3View itw,
                    const ArrayOfGridPosPoly& rgp,
                    const ArrayOfGridPosPoly& cgp );

void interpweights( Tensor4View itw,
                    const ArrayOfGridPosPoly& pgp,
                    const ArrayOfGridPosPoly& rgp,
                    const ArrayOfGridPosPoly& cgp );

void interpweights( Tensor5View itw,
                    const ArrayOfGridPosPoly& bgp,
                    const ArrayOfGridPosPoly& pgp,
                    const ArrayOfGridPosPoly& rgp,
                    const ArrayOfGridPosPoly& cgp );

void interpweights( Tensor6View itw,
                    const ArrayOfGridPosPoly& sgp,
                    const ArrayOfGridPosPoly& bgp,
                    const ArrayOfGridPosPoly& pgp,
                    const ArrayOfGridPosPoly& rgp,
                    const ArrayOfGridPosPoly& cgp );

void interpweights( Tensor7View itw,
                    const ArrayOfGridPosPoly& vgp,
                    const ArrayOfGridPosPoly& sgp,
                    const ArrayOfGridPosPoly& bgp,
                    const ArrayOfGridPosPoly& pgp,
                    const ArrayOfGridPosPoly& rgp,
                    const ArrayOfGridPosPoly& cgp );

void interp( MatrixView            ia,
             ConstTensor3View      itw,
             ConstMatrixView       a,   
             const ArrayOfGridPosPoly& rgp,
             const ArrayOfGridPosPoly& cgp);

void interp( Tensor3View           ia,
             ConstTensor4View      itw,
             ConstTensor3View      a,   
             const ArrayOfGridPosPoly& pgp,
             const ArrayOfGridPosPoly& rgp,
             const ArrayOfGridPosPoly& cgp);

void interp( Tensor4View           ia,
             ConstTensor5View      itw,
             ConstTensor4View      a,   
             const ArrayOfGridPosPoly& bgp,
             const ArrayOfGridPosPoly& pgp,
             const ArrayOfGridPosPoly& rgp,
             const ArrayOfGridPosPoly& cgp);

void interp( Tensor5View           ia,
             ConstTensor6View      itw,
             ConstTensor5View      a,   
             const ArrayOfGridPosPoly& sgp,
             const ArrayOfGridPosPoly& bgp,
             const ArrayOfGridPosPoly& pgp,
             const ArrayOfGridPosPoly& rgp,
             const ArrayOfGridPosPoly& cgp);

void interp( Tensor6View           ia,
             ConstTensor7View      itw,
             ConstTensor6View      a,   
             const ArrayOfGridPosPoly& vgp,
             const ArrayOfGridPosPoly& sgp,
             const ArrayOfGridPosPoly& bgp,
             const ArrayOfGridPosPoly& pgp,
             const ArrayOfGridPosPoly& rgp,
             const ArrayOfGridPosPoly& cgp);


#endif // interpolation_poly_h
