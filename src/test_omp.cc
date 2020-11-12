/* Copyright (C) 2013 Oliver Lemke

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

#include <iostream>
#include "arts_omp.h"
#include "lin_alg.h"
#include "matpackIV.h"
#include "rte.h"

void ext_mat_case(Index& icase,
                  ConstMatrixView ext_mat,
                  const Index stokes_dim) {
  if (icase == 0) {
    icase = 1;  // Start guess is diagonal

    //--- Scalar case ----------------------------------------------------------
    if (stokes_dim == 1) {
    }

    //--- Vector RT ------------------------------------------------------------
    else {
      // Check symmetries and analyse structure of exp_mat:
      assert(ext_mat(1, 1) == ext_mat(0, 0));
      assert(ext_mat(1, 0) == ext_mat(0, 1));

      if (ext_mat(1, 0) != 0) {
        icase = 2;
      }

      if (stokes_dim >= 3) {
        assert(ext_mat(2, 2) == ext_mat(0, 0));
        assert(ext_mat(2, 1) == -ext_mat(1, 2));
        assert(ext_mat(2, 0) == ext_mat(0, 2));

        if (ext_mat(2, 0) != 0 || ext_mat(2, 1) != 0) {
          icase = 3;
        }

        if (stokes_dim > 3) {
          assert(ext_mat(3, 3) == ext_mat(0, 0));
          assert(ext_mat(3, 2) == -ext_mat(2, 3));
          assert(ext_mat(3, 1) == -ext_mat(1, 3));
          assert(ext_mat(3, 0) == ext_mat(0, 3));

          if (icase < 3)  // if icase==3, already at most complex case
          {
            if (ext_mat(3, 0) != 0 || ext_mat(3, 1) != 0) {
              icase = 3;
            } else if (ext_mat(3, 2) != 0) {
              icase = 2;
            }
          }
        }
      }
    }
  }
}

/** Converts an extinction matrix to a transmission matrix

    The function performs the calculations differently depending on the
    conditions, to improve the speed. There are three cases: <br>
       1. Scalar RT and/or the matrix ext_mat_av is diagonal. <br>
       2. Special expression for "azimuthally_random" case. <br>
       3. The total general case.

    If the structure of *ext_mat* is known, *icase* can be set to "case index"
    (1, 2 or 3) and some time is saved. This includes that no asserts are
    performed on *ext_mat*.

    Otherwise, *icase* must be set to 0. *ext_mat* is then analysed and *icase*
    is set by the function and is returned.

    trans_mat must be sized before calling the function.

    @param[out]   trans_mat      Transmission matrix of slab.
    @param[out]   icase          Index giving ext_mat case.
    @param[in]    ext_mat        Averaged extinction matrix.
    @param[in]    lstep          The length of the RTE step.

    @author Patrick Eriksson (based on earlier version started by Claudia)
    @date   2013-05-17 
 */
void ext2trans(MatrixView trans_mat,
               Index& icase,
               ConstMatrixView ext_mat,
               const Numeric& lstep) {
  const Index stokes_dim = ext_mat.ncols();

  assert(ext_mat.nrows() == stokes_dim);
  assert(trans_mat.nrows() == stokes_dim && trans_mat.ncols() == stokes_dim);

  // Theoretically ext_mat(0,0) >= 0, but to demand this can cause problems for
  // iterative retrievals, and the assert is skipped. Negative should be a
  // result of negative vmr, and this issue is checked in basics_checkedCalc.
  //assert( ext_mat(0,0) >= 0 );

  assert(icase >= 0 && icase <= 3);
  assert(!is_singular(ext_mat));
  assert(lstep >= 0);

  // Analyse ext_mat?
  ext_mat_case(icase, ext_mat, stokes_dim);

  // Calculation options:
  if (icase == 1) {
    trans_mat = 0;
    trans_mat(0, 0) = exp(-ext_mat(0, 0) * lstep);
    for (Index i = 1; i < stokes_dim; i++) {
      trans_mat(i, i) = trans_mat(0, 0);
    }
  }

  else if (icase == 2 && stokes_dim < 3) {
    // Expressions below are found in "Polarization in Spectral Lines" by
    // Landi Degl'Innocenti and Landolfi (2004).
    const Numeric tI = exp(-ext_mat(0, 0) * lstep);
    const Numeric HQ = ext_mat(0, 1) * lstep;
    trans_mat(0, 0) = tI * cosh(HQ);
    trans_mat(1, 1) = trans_mat(0, 0);
    trans_mat(1, 0) = -tI * sinh(HQ);
    trans_mat(0, 1) = trans_mat(1, 0);
    /* Does not work for stokes_dim==3, and commnted out 180502 by PE: 
      if( stokes_dim >= 3 )
        {
          trans_mat(2,0) = 0;
          trans_mat(2,1) = 0;
          trans_mat(0,2) = 0;
          trans_mat(1,2) = 0;
          const Numeric RQ = ext_mat(2,3) * lstep;
          trans_mat(2,2) = tI * cos( RQ );
          if( stokes_dim > 3 )
            {
              trans_mat(3,0) = 0;
              trans_mat(3,1) = 0;
              trans_mat(0,3) = 0;
              trans_mat(1,3) = 0;
              trans_mat(3,3) = trans_mat(2,2);
              trans_mat(3,2) = tI * sin( RQ );
              trans_mat(2,3) = -trans_mat(3,2); 
            }
        }
      */
  } else {
    Matrix ext_mat_ds = ext_mat;
    ext_mat_ds *= -lstep;
    //
    Index q = 10;  // index for the precision of the matrix exp function
    //
    switch (stokes_dim) {
      case 4:
        cayley_hamilton_fitted_method_4x4_propmat_to_transmat__eigen(
            trans_mat, ext_mat_ds);
        break;
      default:
        matrix_exp(trans_mat, ext_mat_ds, q);
    }
  }
}

int main() {
  Index nloop = 2000;
  Index nf = 115;
  Index np = 50;
  Index stokes_dim = 1;

  Tensor4 ppath_abs(nf, stokes_dim, stokes_dim, np, 0.);
  Tensor4 trans_cumulat(nf, stokes_dim, stokes_dim, np, 0.);
  Tensor4 trans_partial(nf, stokes_dim, stokes_dim, np, 0.);
  Vector scalar_tau(nf, 0.);
  Vector lstep(np, 1.);

  ArrayOfArrayOfIndex extmat_case(np);
  for (Index ip = 0; ip < np; ip++) {
    extmat_case[ip].resize(nf);
  }

  // Commenting in the next two lines shouldn't change anything, because
  // the outer loop will still run sequentially and the inner one parallel.
  // Nonetheless, performance grinds to a halt because OpenMP can not reuse
  // the initial threads. It has to create new ones for the inner loop on
  // every iteration.
  /*#pragma omp parallel for \
if (0)*/
  for (Index n = 0; n < nloop; n++) {
    for (Index ip = 1; ip < np; ip++) {
#pragma omp parallel for
      for (Index iv = 0; iv < nf; iv++) {
        // Transmission due to absorption
        Matrix ext_mat(stokes_dim, stokes_dim);
        for (Index is1 = 0; is1 < stokes_dim; is1++) {
          for (Index is2 = 0; is2 < stokes_dim; is2++) {
            ext_mat(is1, is2) = 0.5 * (ppath_abs(iv, is1, is2, ip - 1) +
                                       ppath_abs(iv, is1, is2, ip));
          }
        }
        scalar_tau[iv] += lstep[ip - 1] * ext_mat(0, 0);
        extmat_case[ip - 1][iv] = 0;
        ext2trans(trans_partial(iv, joker, joker, ip - 1),
                  extmat_case[ip - 1][iv],
                  ext_mat,
                  lstep[ip - 1]);

        // Cumulative transmission
        // (note that multiplication below depends on ppath loop order)
        mult(trans_cumulat(iv, joker, joker, ip),
             trans_cumulat(iv, joker, joker, ip - 1),
             trans_partial(iv, joker, joker, ip - 1));
      }
    }
  }
  return 0;
}
