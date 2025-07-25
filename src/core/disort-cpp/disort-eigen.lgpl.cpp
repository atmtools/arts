#include "disort-eigen.h"

#include <lin_alg.h>

Index diagonalize_inplace(MatrixView P,
                          VectorView W,
                          MatrixView A,
                          real_diagonalize_workdata& workdata) {
  /* Wraps the standard lin-alg diagonalize method
   */

  diagonalize_inplace(P, W, workdata.imag, A, workdata.workdata);

  return 0;
}
