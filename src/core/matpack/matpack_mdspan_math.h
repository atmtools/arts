#pragma once

#include "matpack_mdspan_strided_view_t.h"

/** Makes A = alpha * B * C + beta * A
 * 
 * @param[out] A May not point at the same data as B or C
 * @param B Any matrix
 * @param C Any matrix
 */
void mult(StridedMatrixView A,
          const StridedConstMatrixView &B,
          const StridedConstMatrixView &C,
          Numeric alpha = 1.0,
          Numeric beta  = 0.0);

/** Makes A = B * C
 * 
 * @param[out] A May not point at the same data as B or C
 * @param B Any matrix
 * @param C Any matrix
 */
void mult(StridedComplexMatrixView A,
          const StridedConstComplexMatrixView &B,
          const StridedConstComplexMatrixView &C,
          Complex alpha = Complex{1.0},
          Complex beta  = Complex{0.0});

/** Makes y = alpha * M * x + beta * y
 * 
 * @param[out] y May not point at the same data as M or x
 * @param M Any matrix
 * @param x Any vector
 */
void mult(StridedVectorView y,
          const StridedConstMatrixView &M,
          const StridedConstVectorView &x,
          Numeric alpha = 1.0,
          Numeric beta  = 0.0);

/** Makes y = M * x
 * 
 * @param[out] y May not point at the same data as M or x
 * @param M Any matrix
 * @param x Any vector
 */
void mult(StridedComplexVectorView y,
          const StridedConstComplexMatrixView &M,
          const StridedConstComplexVectorView &x,
          Complex alpha = Complex{1.0},
          Complex beta  = Complex{0.0});
