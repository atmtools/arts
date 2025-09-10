#include "disort-eigen.h"

#include <arts_constexpr_math.h>
#include <matpack.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace {
/*============================= c_asymmetric_matrix() ===================*/

/*
  Solves eigenfunction problem for real asymmetric matrix for which it
  is known a priori that the eigenvalues are real. This is an adaptation
  of a subroutine EIGRF in the IMSL library to use real instead of complex
  arithmetic, accounting for the known fact that the eigenvalues and
  eigenvectors in the discrete ordinate solution are real.
  
  EIGRF is based primarily on EISPACK routines.  The matrix is first
  balanced using the Parlett-Reinsch algorithm.  Then the Martin-Wilkinson
  algorithm is applied. There is a statement 'j = wk(i)' that converts a
  Numeric precision variable to an Indexeger variable; this seems dangerous
  to us in principle, but seems to work fine in practice.
  
  References:

  Dongarra, J. and C. Moler, EISPACK -- A Package for Solving Matrix
      Eigenvalue Problems, in Cowell, ed., 1984: Sources and Development of
      Mathematical Software, Prentice-Hall, Englewood Cliffs, NJ
  Parlett and Reinsch, 1969: Balancing a Matrix for Calculation of
      Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
  Wilkinson, J., 1965: The Algebraic Eigenvalue Problem, Clarendon Press,
      Oxford

   I N P U T    V A R I A B L E S:

       aa    :  input asymmetric matrix, destroyed after solved
        m    :  order of aa
       ia    :  first dimension of aa
    ievec    :  first dimension of evec

   O U T P U T    V A R I A B L E S:

       evec  :  (unnormalized) eigenvectors of aa (column j corresponds to W[J-1])
       eval  :  (unordered) eigenvalues of aa (dimension m)
       ier   :  if != 0, signals that W[ier-1] failed to converge;
                   in that case eigenvalues ier+1,ier+2,...,m  are
                   correct but eigenvalues 1,...,ier are set to zero.

   S C R A T C H   V A R I A B L E S:

       wk    :  work area (dimension at least 2*m)
       
   Called by- c_solve_eigen
   Calls- c_errmsg
 -------------------------------------------------------------------*/

Index c_asymmetric_matrix(MatrixView P,
                          VectorView W,
                          MatrixView A,
                          VectorView work,
                          std::span<Index> iwork) {
  using Math::pow2;
  using nonstd::abs;
  using std::copysign;
  using std::min;
  using std::sqrt;

  constexpr Numeric c1 = 0.4375, c2 = 0.5, c3 = 0.75, c4 = 0.95, c5 = 16,
                    c6 = 256;

  constexpr Numeric tol = std::numeric_limits<Numeric>::epsilon();

  const Index m = W.size();

  assert(work.size() >= static_cast<Size>(2 * m));
  assert(iwork.size() >= static_cast<Size>(m));
  assert(m == A.nrows() and m == A.ncols());
  assert(m == P.nrows() and m == P.ncols());
  assert(P.shape() == A.shape());

  const Numeric rnorm = abssum(A);

  /*
   * Handle 1x1 and 2x2 special cases
   */
  if (m == 1) {
    W[0]    = A[0, 0];
    P[0, 0] = 1;

    return 0;
  }

  if (m == 2) {
    const Numeric discri = pow2(A[0, 0] - A[1, 1]) + 4 * A[0, 1] * A[1, 0];
    const Numeric sgn    = (A[0, 0] < A[1, 1]) ? -1 : 1;

    W[0]    = 0.5 * (A[0, 0] + A[1, 1] + sgn * sqrt(discri));
    W[1]    = 0.5 * (A[0, 0] + A[1, 1] - sgn * sqrt(discri));
    P[0, 0] = 1;
    P[1, 1] = 1;

    if (A[0, 0] == A[1, 1] && (A[1, 0] == 0 || A[0, 1] == 0)) {
      const Numeric w = tol * rnorm;
      P[1, 0]         = A[1, 0] / w;
      P[0, 1]         = -A[0, 1] / w;
    } else {
      P[1, 0] = A[1, 0] / (W[0] - A[1, 1]);
      P[0, 1] = A[0, 1] / (W[1] - A[0, 0]);
    }

    return 0;
  }

  /*
   * Initialize output variables
   */
  W           = 0;
  P           = 0;
  diagonal(P) = 1;

  // lowest and highest valid index - to be adjusted
  Index l = 0, k = m - 1;

  /*
   * Balance the input matrix and reduce its norm by diagonal similarity transformation stored in wk;
   * then search for rows isolating an eigenvalue and push them down.
   */

  for (Index j = k; j >= 0; j--) {
    const Range oj{0, j}, jm{j + 1, m - j - 1};
    const Numeric row = abssum(A[j, oj]) + abssum(A[j, jm]);

    if (row == 0) {
      iwork[k] = j;

      if (j != k) {
        stdr::swap_ranges(A[Range{0, k + 1}, j], A[Range{0, k + 1}, k]);
        stdr::swap_ranges(A[j], A[k]);
      }

      j = k--;
    }
  }

  const Range zk{0, k + 1};

  /*
   * Search for columns isolating an eigenvalue and push them left.
   */

  for (Index j = l; j <= k; j++) {
    const Range oj{0, j}, jm{j + 1, m - j - 1};
    const Numeric col = abssum(A[oj, j]) + abssum(A[jm, j]);

    if (col == 0) {
      iwork[l] = j;

      if (j != l) {
        stdr::swap_ranges(A[zk, j], A[zk, l]);
        stdr::swap_ranges(A[j, Range{l, m - l}], A[l, Range{l, m - l}]);
      }

      j = l++;
    }
  }

  const Range lm{l, m - l};
  const Range lk{l, k - l + 1};

  /*
   * Balance the submatrix in rows L through K
   */
  work[lk] = 1;

  for (bool noconv = true; noconv;) {
    noconv = false;
    for (Index i = l; i <= k; i++) {
      Numeric col = 0, row = 0;

      for (Index j = l; j <= k; j++) {
        if (j != i) {
          col += abs(A[j, i]);
          row += abs(A[i, j]);
        }
      }

      Numeric f       = 1;
      Numeric g       = row / c5;
      const Numeric h = col + row;

      while (col < g) {
        f   *= c5;
        col *= c6;
      }

      g = row * c5;

      while (col >= g) {
        f   /= c5;
        col /= c6;
      }

      /*
       * Now balance
       */
      if ((col + row) / f < c4 * h) {
        work[i] *= f;
        noconv   = true;

        A[i, lm] /= f;
        A[zk, i] *= f;
      }
    }
  }

  Index n1{}, n2{};

  if (k - 1 >= l + 1) {
    VectorView rwork = work[Range{m, m}];
    /*
     * Transfer A to a Hessenberg form.
     */
    for (Index n = l + 1; n <= k - 1; n++) {
      Numeric h = 0;
      rwork[n]  = 0;

      /*
       * Scale column
       */
      const Numeric scale = abssum(A[Range{n, k - n + 1}, n - 1]);

      if (scale != 0) {
        for (Index i = k; i >= n; i--) {
          rwork[i]  = A[i, n - 1] / scale;
          h        += pow2(rwork[i]);
        }

        const Numeric g  = -copysign(sqrt(h), rwork[n]);
        h               -= rwork[n] * g;
        rwork[n]        -= g;

        /*
         * Form (I-(U*UT)/H)*A
         */
        for (Index j = n; j < m; j++) {
          const Range nk{n, k - n + 1};
          const Numeric f = dot(rwork[nk], A[nk, j]) / h;

          for (Index i = n; i <= k; i++) A[i, j] -= rwork[i] * f;
        }

        /*
         * Form (i-(u*ut)/h)*a*(i-(u*ut)/h)
         */
        for (Index i = 0; i <= k; i++) {
          const Range nk{n, k - n + 1};
          const Numeric f = dot(rwork[nk], A[i, nk]) / h;

          for (Index j = n; j <= k; j++) A[i, j] -= rwork[j] * f;
        }

        rwork[n]    *= scale;
        A[n, n - 1]  = scale * g;
      }
    }

    for (Index n = k - 2; n >= l; n--) {
      n1        = n + 1;
      n2        = n + 2;
      Numeric f = A[n + 1, n];

      if (f != 0) {
        f *= rwork[n + 1];

        const Range rs{n + 2, k - n - 1};
        rwork[rs] = A[rs, n];

        if (n + 1 <= k) {
          for (Index j = 0; j < m; j++) {
            const Range nk{n + 1, k - n};
            const Numeric g = dot(rwork[nk], P[nk, j]) / f;

            for (Index i = n + 1; i <= k; i++) P[i, j] += g * rwork[i];
          }
        }
      }
    }
  }

  const Range zl{0, l};
  const Range km{k + 1, m - k - 1};
  auto dA = diagonal(A);
  W[zl]   = dA[zl];
  W[km]   = dA[km];

  Index i, in, j, ka, lb = 0;
  Numeric p = 0, q = 0, r = 0, s, uu, vv, w, x, y, z;

  Index n   = k;
  Numeric t = 0;
  /*
   * Search for next eigenvalues
   */

  do {
    if (n < l) goto S550;

    in = 0;
    n1 = n - 1;
    n2 = n - 2;

    /*
     * Look for single small sub-diagonal element
     */

  S410:

    for (i = l; i <= n; i++) {
      lb = n + l - i;

      if (lb == l) break;

      s = abs(dA[lb - 1]) + abs(dA[lb]);
      if (s == 0) s = rnorm;

      if (abs(A[lb, lb - 1]) <= tol * s) break;
    }

    x = dA[n];

    if (lb == n) {
      /*
       * One eigenvalue found
       */
      dA[n] = x + t;
      W[n]  = dA[n];
      n     = n1;
      continue;
    }

    y = dA[n1];
    w = A[n, n1] * A[n1, n];

    if (lb == n1) {
      /*
       * Two eigenvalues found
       */
      p      = (y - x) * c2;
      q      = p * p + w;
      z      = sqrt(abs(q));
      dA[n]  = x + t;
      x      = dA[n];
      dA[n1] = y + t;

      /*
       * Real pair
       */
      z     = p + copysign(z, p);
      W[n1] = x + z;
      W[n]  = W[n1];

      if (z != 0) W[n] = x - w / z;

      x = A[n, n1];

      /*
       * Employ scale factor in case X and Z are very small
       */
      r = sqrt(x * x + z * z);
      p = x / r;
      q = z / r;

      /*
       * Row modification
       */
      for (j = n1; j < m; j++) {
        z        = A[n1, j];
        A[n1, j] = q * z + p * A[n, j];
        A[n, j]  = -p * z + q * A[n, j];
      }

      /*
       * Column modification
       */
      for (i = 0; i <= n; i++) {
        z        = A[i, n1];
        A[i, n1] = q * z + p * A[i, n];
        A[i, n]  = -p * z + q * A[i, n];
      }

      /*
       * Accumulate transformations
       */
      for (i = l; i <= k; i++) {
        z        = P[i, n1];
        P[i, n1] = q * z + p * P[i, n];
        P[i, n]  = -p * z + q * P[i, n];
      }

      n = n2;
      continue;
    }

    break;
  } while (true);

  /*
   * No convergence after 30 iterations; set error indicator to
   * the index of the current eigenvalue, and return.
   */
  if (in == 30) return n;

  /*
   * Form shift
   */
  if (in == 10 || in == 20) {
    t                       += x;
    dA[Range{l, n - l + 1}] -= x;

    s = abs(A[n, n1]) + abs(A[n1, n2]);
    x = c3 * s;
    y = x;
    w = -c1 * s * s;
  }

  in++;

  /*
   * Look for two consecutive small sub-diagonal elements
   */
  for (j = lb; j <= n2; j++) {
    i  = n2 + lb - j;
    z  = dA[i];
    r  = x - z;
    s  = y - z;
    p  = (r * s - w) / A[i + 1, i] + A[i, i + 1];
    q  = dA[i + 1] - z - r - s;
    r  = A[i + 2, i + 1];
    s  = abs(p) + abs(q) + abs(r);
    p /= s;
    q /= s;
    r /= s;

    if (i == lb) break;

    uu = abs(A[i, i - 1]) * (abs(q) + abs(r));
    vv = abs(p) * (abs(dA[i - 1]) + abs(z) + abs(dA[i + 1]));

    if (uu <= tol * vv) break;
  }

  A[i + 2, i] = 0;
  for (j = i + 3; j <= n; j++) {
    A[j, j - 2] = 0;
    A[j, j - 3] = 0;
  }

  /*
   * Double QR step involving rows K to N and columns M to N
   */
  for (ka = i; ka <= n1; ka++) {
    const bool notlas = (ka != n1);
    if (ka == i) {
      s = copysign(sqrt(p * p + q * q + r * r), p);
      if (lb != i) A[ka, ka - 1] *= -1;

    } else {
      p = A[ka, ka - 1];
      q = A[ka + 1, ka - 1];
      r = 0;

      if (notlas) r = A[ka + 2, ka - 1];

      x = abs(p) + abs(q) + abs(r);
      if (x == 0) continue;

      p /= x;
      q /= x;
      r /= x;
      s  = copysign(sqrt(p * p + q * q + r * r), p);

      A[ka, ka - 1] = -s * x;
    }

    p += s;
    x  = p / s;
    y  = q / s;
    z  = r / s;
    q /= p;
    r /= p;

    /*
     * Row modification
     */
    for (j = ka; j < m; j++) {
      p = A[ka, j] + q * A[ka + 1, j];

      if (notlas) {
        p            += r * A[ka + 2, j];
        A[ka + 2, j] -= p * z;
      }

      A[ka + 1, j] -= p * y;
      A[ka, j]     -= p * x;
    }

    /*
     * Column modification
     */
    for (Index ii = 0; ii <= min(n, ka + 3); ii++) {
      p = x * A[ii, ka] + y * A[ii, ka + 1];

      if (notlas) {
        p             += z * A[ii, ka + 2];
        A[ii, ka + 2] -= p * r;
      }

      A[ii, ka + 1] -= p * q;
      A[ii, ka]     -= p;
    }

    /*
     * Accumulate transformations
     */
    for (Index ii = l; ii <= k; ii++) {
      p = x * P[ii, ka] + y * P[ii, ka + 1];

      if (notlas) {
        p             += z * P[ii, ka + 2];
        P[ii, ka + 2] -= p * r;
      }

      P[ii, ka + 1] -= p * q;
      P[ii, ka]     -= p;
    }
  }

  goto S410;

  /*
   * All evals found, now backsubstitute real vector
   */

S550:

  if (rnorm != 0) {
    for (n = m - 1; n >= 0; n--) {
      n2    = n;
      dA[n] = 1;

      for (i = n - 1; i >= 0; i--) {
        w = dA[i] - W[n];
        if (w == 0) w = tol * rnorm;

        r = A[i, n] + dot(A[i, Range(n2, n - n2)], A[Range(n2, n - n2), n]);

        A[i, n] = -r / w;
        n2      = i;
      }
    }
    /*
     * End backsubstitution vectors of isolated evals
     */
    for (i = 0; i < m; i++) {
      if (i < l || i > k) {
        const Range im{i, m - i};
        P[i, im] = A[i, im];
      }
    }
    /*
     * Multiply by transformation matrix
     */
    if (k != 0) {
      for (j = m - 1; j >= l; j--) {
        for (i = l; i <= k; i++) {
          const Range nk(l, min(j, k) - l + 1);
          //! NOTE: This can be turned into a mult(vec, mat, vec)
          P[i, j] = dot(P[i, nk], A[nk, j]);
        }
      }
    }
  }

  for (i = l; i <= k; i++) P[i] *= work[i];

  /*
   * Interchange rows if permutations occurred
   */
  for (i = l - 1; i >= 0; i--) {
    j = iwork[i];
    if (i != j) stdr::swap_ranges(P[i], P[j]);
  }

  for (i = k + 1; i < m; i++) {
    j = iwork[i];
    if (i != j) stdr::swap_ranges(P[i], P[j]);
  }

  return 0;
}
}  // namespace

Index diagonalize_inplace(MatrixView P,
                          VectorView W,
                          MatrixView A,
                          real_diagonalize_workdata& workdata) {
  /* Note that this is code ported from cdisort.

  If it turns out that porting the code has broken something,
  it is relatively easy to just use the original.

  The return of this function is the int-pointer that is part
  of the last variables.

  The evec matrix must be transposed again upon output.

  Please use static casts instead of implicit casts as well,
  because there are weird int-double casts going on.
*/

  return c_asymmetric_matrix(P, W, A, workdata.reals, workdata.ints);
}
