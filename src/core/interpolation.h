/*!
  \file   interpolation.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri May  3 08:54:45 2002
  
  \brief  Header file for interpolation.cc.
*/

#ifndef interpolation_h
#define interpolation_h

#include <matpack.h>
#include <mystring.h>

//! Structure to store a grid position.
/*! 
  A grid position specifies, where an interpolation point is, relative
  to the original grid. It consists of three parts, an Index giving the
  original grid index below the interpolation point, a Numeric
  giving the fractional distance to the next original grid point, and a
  Numeric giving 1 minus this number. Of course, the last element is
  redundant. However, it is efficient to store this, since it is used
  many times over. We store the two Numerics in a plain C array of
  dimension 2. (No need to use fancy Array or Vector for this, since
  the dimension is fixed.)

  For example, idx=3 and fd=.5 means that this interpolation point is
  half-way between index 3 and 4 of the original grid.

  Note, that below in the first paragraph means "with a lower
  index". If the original grid is sorted in descending order, the
  value at the grid point below the interpolation point will be
  numerically higher than the interpolation point.

  In other words, grid positions and fractional distances are defined
  relative to the order of the original grid. Examples:

  old grid = 2 3
  new grid = 2.25
  idx      = 0
  fd[0]    = 0.25

  old grid = 3 2
  new grid = 2.25
  idx      = 0
  fd[0]    = 0.75

  Note that fd[0] is different in the second case, because the old grid
  is sorted in descending order. Note also that idx is the same in
  both cases.

  Grid positions for a whole new grid are stored in an Array<GridPos>
  (called ArrayOfGridPos). 
*/
struct GridPos {
  Index idx; /*!< Original grid index below interpolation point. */
  std::array<Numeric, 2> fd; /*!< Fractional distance to next point
                                    (0<=fd[0]<=1), fd[1] = 1-fd[0]. */

  friend std::ostream& operator<<(std::ostream& os, const GridPos& gp);
};

//! An Array of grid positions.
/*! 
  See \ref GridPos for details.
*/

typedef Array<GridPos> ArrayOfGridPos;
typedef Array<Array<GridPos> > ArrayOfArrayOfGridPos;
typedef Array<Array<Array<GridPos> > > ArrayOfArrayOfArrayOfGridPos;
typedef Array<Array<Array<Array<GridPos> > > >
    ArrayOfArrayOfArrayOfArrayOfGridPos;

std::ostream& operator<<(std::ostream& os, const ArrayOfGridPos& a);

std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfGridPos& a);

std::ostream& operator<<(std::ostream& os,
                         const ArrayOfArrayOfArrayOfGridPos& a);

std::ostream& operator<<(std::ostream& os,
                         const ArrayOfArrayOfArrayOfArrayOfGridPos& a);

// Function headers (documentation is in .cc file):

void gridpos(ArrayOfGridPos& gp,
             StridedConstVectorView old_grid,
             StridedConstVectorView new_grid,
             const Numeric& extpolfac = 0.5);

void gridpos(GridPos& gp,
             StridedConstVectorView old_grid,
             const Numeric& new_grid,
             const Numeric& extpolfac = 0.5);

void gridpos_1to1(ArrayOfGridPos& gp, StridedConstVectorView grid);

void gridpos_copy(GridPos& gp_new, const GridPos& gp_old);

Numeric fractional_gp(const GridPos& gp);

void gridpos_check_fd(GridPos& gp);

void gridpos_force_end_fd(GridPos& gp, const Index& n);

void gridpos_upperend_check(GridPos& gp, const Index& ie);

void gridpos_upperend_check(ArrayOfGridPos& gp, const Index& ie);

void gp4length1grid(ArrayOfGridPos& gp);

bool is_gridpos_at_index_i(const GridPos& gp,
                           const Index& i,
                           const bool& strict = true);

Index gridpos2gridrange(const GridPos& gp, const bool& upwards);

////////////////////////////////////////////////////////////////////////////
//                      Red Interpolation
////////////////////////////////////////////////////////////////////////////

void interpweights(StridedVectorView itw, const GridPos& tc);

void interpweights(StridedVectorView itw, const GridPos& tr, const GridPos& tc);

void interpweights(StridedVectorView itw,
                   const GridPos& tp,
                   const GridPos& tr,
                   const GridPos& tc);

void interpweights(StridedVectorView itw,
                   const GridPos& tb,
                   const GridPos& tp,
                   const GridPos& tr,
                   const GridPos& tc);

void interpweights(StridedVectorView itw,
                   const GridPos& ts,
                   const GridPos& tb,
                   const GridPos& tp,
                   const GridPos& tr,
                   const GridPos& tc);

void interpweights(StridedVectorView itw,
                   const GridPos& tv,
                   const GridPos& ts,
                   const GridPos& tb,
                   const GridPos& tp,
                   const GridPos& tr,
                   const GridPos& tc);

Numeric interp(StridedConstVectorView itw,
               StridedConstVectorView a,
               const GridPos& tc);

Numeric interp(StridedConstVectorView itw,
               StridedConstMatrixView a,
               const GridPos& tr,
               const GridPos& tc);

Numeric interp(StridedConstVectorView itw,
               StridedConstTensor3View a,
               const GridPos& tp,
               const GridPos& tr,
               const GridPos& tc);

Numeric interp(StridedConstVectorView itw,
               StridedConstTensor4View a,
               const GridPos& tb,
               const GridPos& tp,
               const GridPos& tr,
               const GridPos& tc);

Numeric interp(StridedConstVectorView itw,
               StridedConstTensor5View a,
               const GridPos& ts,
               const GridPos& tb,
               const GridPos& tp,
               const GridPos& tr,
               const GridPos& tc);

Numeric interp(StridedConstVectorView itw,
               StridedConstTensor6View a,
               const GridPos& tv,
               const GridPos& ts,
               const GridPos& tb,
               const GridPos& tp,
               const GridPos& tr,
               const GridPos& tc);

////////////////////////////////////////////////////////////////////////////
//                      Blue interpolation
////////////////////////////////////////////////////////////////////////////

void interpweights(StridedMatrixView itw, const ArrayOfGridPos& cgp);

void interpweights(StridedMatrixView itw,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp);

void interpweights(StridedMatrixView itw,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp);

void interpweights(StridedMatrixView itw,
                   const ArrayOfGridPos& bgp,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp);

void interpweights(StridedMatrixView itw,
                   const ArrayOfGridPos& sgp,
                   const ArrayOfGridPos& bgp,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp);

void interpweights(StridedMatrixView itw,
                   const ArrayOfGridPos& vgp,
                   const ArrayOfGridPos& sgp,
                   const ArrayOfGridPos& bgp,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp);

void interp(StridedVectorView ia,
            StridedConstMatrixView itw,
            StridedConstVectorView a,
            const ArrayOfGridPos& cgp);

void interp(VectorView ia,
            ConstMatrixView itw,
            ConstMatrixView a,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp);

void interp(StridedVectorView ia,
            StridedConstMatrixView itw,
            StridedConstTensor3View a,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp);

void interp(StridedVectorView ia,
            StridedConstMatrixView itw,
            StridedConstTensor4View a,
            const ArrayOfGridPos& bgp,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp);

void interp(StridedVectorView ia,
            StridedConstMatrixView itw,
            StridedConstTensor5View a,
            const ArrayOfGridPos& sgp,
            const ArrayOfGridPos& bgp,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp);

void interp(StridedVectorView ia,
            StridedConstMatrixView itw,
            StridedConstTensor6View a,
            const ArrayOfGridPos& vgp,
            const ArrayOfGridPos& sgp,
            const ArrayOfGridPos& bgp,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp);

////////////////////////////////////////////////////////////////////////////
//                      Green interpolation
////////////////////////////////////////////////////////////////////////////

void interpweights(StridedTensor3View itw,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp);

void interpweights(StridedTensor4View itw,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp);

void interpweights(StridedTensor5View itw,
                   const ArrayOfGridPos& bgp,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp);

void interpweights(StridedTensor6View itw,
                   const ArrayOfGridPos& sgp,
                   const ArrayOfGridPos& bgp,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp);

void interpweights(StridedTensor7View itw,
                   const ArrayOfGridPos& vgp,
                   const ArrayOfGridPos& sgp,
                   const ArrayOfGridPos& bgp,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp);

void interp(StridedMatrixView ia,
            StridedConstTensor3View itw,
            StridedConstMatrixView a,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp);

void interp(StridedTensor3View ia,
            StridedConstTensor4View itw,
            StridedConstTensor3View a,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp);

void interp(StridedTensor4View ia,
            StridedConstTensor5View itw,
            StridedConstTensor4View a,
            const ArrayOfGridPos& bgp,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp);

void interp(StridedTensor5View ia,
            StridedConstTensor6View itw,
            StridedConstTensor5View a,
            const ArrayOfGridPos& sgp,
            const ArrayOfGridPos& bgp,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp);

void interp(Tensor6View ia,
            ConstTensor7View itw,
            ConstTensor6View a,
            const ArrayOfGridPos& vgp,
            const ArrayOfGridPos& sgp,
            const ArrayOfGridPos& bgp,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp);

Numeric interp_poly(ConstVectorView x,
                    StridedConstVectorView y,
                    const Numeric& x_i,
                    const GridPos& gp);

#endif  // interpolation_h
