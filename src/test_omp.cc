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
#include "matpackIV.h"
#include "rte.h"

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
