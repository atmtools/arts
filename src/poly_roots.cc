/* Copyright (C) 2002-2012 Oliver Lemke <olemke@core-dump.info>

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

////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
   \file   poly_roots.cc

   Contains the code to determine roots of polynomials.

   Code was taken from the GNU Scientific library.
   http://sources.redhat.com/gsl/

   \author Oliver Lemke
   \date 2002-03-06
*/

////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "poly_roots.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

/* C-style matrix elements */
#define MAT(m, i, j, n) ((m)[(i) * (n) + (j)])

/* Fortran-style matrix elements */
#define FMAT(m, i, j, n) ((m)[((i)-1) * (n) + ((j)-1)])

#define GSL_DBL_EPSILON 2.2204460492503131e-16

#define RADIX 2
#define RADIX2 (RADIX * RADIX)

#define GSL_SUCCESS 0
#define GSL_FAILURE -1
#define GSL_EINVAL 4
#define GSL_EFAILED 5

#define GSL_SET_COMPLEX_PACKED(zp, n, x, y) \
  do {                                      \
    *((zp) + 2 * (n)) = (x);                \
    *((zp) + (2 * (n) + 1)) = (y);          \
  } while (0)

typedef double *gsl_complex_packed_ptr;

typedef struct {
  size_t nc;
  double *matrix;
} gsl_poly_complex_workspace;

/* Begin Internal GSL function prototypes */

static gsl_poly_complex_workspace *gsl_poly_complex_workspace_alloc(size_t n);

static void gsl_poly_complex_workspace_free(gsl_poly_complex_workspace *w);

static void set_companion_matrix(const double *a, size_t n, double *m);

static int qr_companion(double *h, size_t nc, gsl_complex_packed_ptr z);

static void balance_companion_matrix(double *m, size_t n);

static int gsl_poly_complex_solve(const double *a,
                                  size_t n,
                                  gsl_poly_complex_workspace *w,
                                  gsl_complex_packed_ptr z);

/* End Internal GSL function prototypes */

int poly_root_solve(Matrix &roots, Vector &coeffs) {
  Index a;
  double *c;
  double *s;

  a = coeffs.nelem();

  ARTS_ASSERT(roots.nrows() == a - 1);
  ARTS_ASSERT(roots.ncols() == 2);
  ARTS_ASSERT(coeffs[a - 1] != 0);

#ifdef USE_DOUBLE
  c = coeffs.data_handle();
  s = roots.data_handle();
#else
  c = (double *)malloc(a * sizeof(double));
  for (Index i = 0; i < a; i++) c[i] = (double)coeffs[i];

  s = (double *)malloc((a - 1) * 2 * sizeof(double));
#endif

  gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(a);

  int status = gsl_poly_complex_solve(c, a, w, s);

#ifndef USE_DOUBLE
  if (status == GSL_SUCCESS) {
    for (Index i = 0; i < a - 1; i++) {
      roots(i, 0) = (Numeric)s[i * 2];
      roots(i, 1) = (Numeric)s[i * 2 + 1];
    }
  }

  free(c);
  free(s);
#endif

  gsl_poly_complex_workspace_free(w);

  return status;
}

static gsl_poly_complex_workspace *gsl_poly_complex_workspace_alloc(size_t n) {
  size_t nc;

  gsl_poly_complex_workspace *w;

  if (n == 0) {
    cerr << "matrix size n must be positive integer" << endl;

    return (NULL);
  }

  w = (gsl_poly_complex_workspace *)malloc(sizeof(gsl_poly_complex_workspace));

  if (w == 0) {
    cerr << "failed to allocate space for struct" << endl;

    return (NULL);
  }

  nc = n - 1;

  w->nc = nc;

  w->matrix = (double *)malloc(nc * nc * sizeof(double));

  if (w->matrix == 0) {
    free(w); /* error in constructor, avoid memory leak */

    cerr << "failed to allocate space for workspace matrix" << endl;

    return (NULL);
  }

  return w;
}

static void gsl_poly_complex_workspace_free(gsl_poly_complex_workspace *w) {
  free(w->matrix);
  free(w);
}

static void balance_companion_matrix(double *m, size_t nc) {
  int not_converged = 1;

  double row_norm = 0;
  double col_norm = 0;

  while (not_converged) {
    size_t i, j;
    double g, f, s;

    not_converged = 0;

    for (i = 0; i < nc; i++) {
      /* column norm, excluding the diagonal */

      if (i != nc - 1) {
        col_norm = fabs(MAT(m, i + 1, i, nc));
      } else {
        col_norm = 0;

        for (j = 0; j < nc - 1; j++) {
          col_norm += fabs(MAT(m, j, nc - 1, nc));
        }
      }

      /* row norm, excluding the diagonal */

      if (i == 0) {
        row_norm = fabs(MAT(m, 0, nc - 1, nc));
      } else if (i == nc - 1) {
        row_norm = fabs(MAT(m, i, i - 1, nc));
      } else {
        row_norm = (fabs(MAT(m, i, i - 1, nc)) + fabs(MAT(m, i, nc - 1, nc)));
      }

      if (col_norm == 0 || row_norm == 0) {
        continue;
      }

      g = row_norm / RADIX;
      f = 1;
      s = col_norm + row_norm;

      while (col_norm < g) {
        f *= RADIX;
        col_norm *= RADIX2;
      }

      g = row_norm * RADIX;

      while (col_norm > g) {
        f /= RADIX;
        col_norm /= RADIX2;
      }

      if ((row_norm + col_norm) < 0.95 * s * f) {
        not_converged = 1;

        g = 1 / f;

        if (i == 0) {
          MAT(m, 0, nc - 1, nc) *= g;
        } else {
          MAT(m, i, i - 1, nc) *= g;
          MAT(m, i, nc - 1, nc) *= g;
        }

        if (i == nc - 1) {
          for (j = 0; j < nc; j++) {
            MAT(m, j, i, nc) *= f;
          }
        } else {
          MAT(m, i + 1, i, nc) *= f;
        }
      }
    }
  }
}

static int qr_companion(double *h, size_t nc, gsl_complex_packed_ptr zroot) {
  double t = 0.0;

  size_t iterations, e, i, j, k, m;

  double w, x, y, s, z;

  double p = 0, q = 0, r = 0;

  /* FIXME: if p,q,r, are not set to zero then the compiler complains
     that they ``might be used uninitialized in this
     function''. Looking at the code this does seem possible, so this
     should be checked. */

  int notlast;

  size_t n = nc;

next_root:

  if (n == 0) return GSL_SUCCESS;

  iterations = 0;

next_iteration:

  for (e = n; e >= 2; e--) {
    double a1 = fabs(FMAT(h, e, e - 1, nc));
    double a2 = fabs(FMAT(h, e - 1, e - 1, nc));
    double a3 = fabs(FMAT(h, e, e, nc));

    if (a1 <= GSL_DBL_EPSILON * (a2 + a3)) break;
  }

  x = FMAT(h, n, n, nc);

  if (e == n) {
    GSL_SET_COMPLEX_PACKED(zroot, n - 1, x + t, 0); /* one real root */
    n--;
    goto next_root;
    /*continue;*/
  }

  y = FMAT(h, n - 1, n - 1, nc);
  w = FMAT(h, n - 1, n, nc) * FMAT(h, n, n - 1, nc);

  if (e == n - 1) {
    p = (y - x) / 2;
    q = p * p + w;
    y = sqrt(fabs(q));

    x += t;

    if (q > 0) /* two real roots */
    {
      if (p < 0) y = -y;
      y += p;

      GSL_SET_COMPLEX_PACKED(zroot, n - 1, x - w / y, 0);
      GSL_SET_COMPLEX_PACKED(zroot, n - 2, x + y, 0);
    } else {
      GSL_SET_COMPLEX_PACKED(zroot, n - 1, x + p, -y);
      GSL_SET_COMPLEX_PACKED(zroot, n - 2, x + p, y);
    }
    n -= 2;

    goto next_root;
    /*continue;*/
  }

  /* No more roots found yet, do another iteration */

  if (iterations == 60) /* increased from 30 to 60 */
  {
    /* too many iterations - give up! */

    return GSL_FAILURE;
  }

  if (iterations % 10 == 0 && iterations > 0) {
    /* use an exceptional shift */

    t += x;
    for (i = 1; i <= n; i++) {
      FMAT(h, i, i, nc) -= x;
    }

    s = fabs(FMAT(h, n, n - 1, nc)) + fabs(FMAT(h, n - 1, n - 2, nc));
    y = 0.75 * s;
    x = y;
    w = -0.4375 * s * s;
  }

  iterations++;

  for (m = n - 2; m >= e; m--) {
    double a1, a2, a3;

    z = FMAT(h, m, m, nc);
    r = x - z;
    s = y - z;
    p = FMAT(h, m, m + 1, nc) + (r * s - w) / FMAT(h, m + 1, m, nc);
    q = FMAT(h, m + 1, m + 1, nc) - z - r - s;
    r = FMAT(h, m + 2, m + 1, nc);
    s = fabs(p) + fabs(q) + fabs(r);
    p /= s;
    q /= s;
    r /= s;

    if (m == e) break;

    a1 = fabs(FMAT(h, m, m - 1, nc));
    a2 = fabs(FMAT(h, m - 1, m - 1, nc));
    a3 = fabs(FMAT(h, m + 1, m + 1, nc));

    if (a1 * (fabs(q) + fabs(r)) <= GSL_DBL_EPSILON * fabs(p) * (a2 + a3))
      break;
  }

  for (i = m + 2; i <= n; i++) {
    FMAT(h, i, i - 2, nc) = 0;
  }

  for (i = m + 3; i <= n; i++) {
    FMAT(h, i, i - 3, nc) = 0;
  }

  /* double QR step */

  for (k = m; k <= n - 1; k++) {
    notlast = (k != n - 1);

    if (k != m) {
      p = FMAT(h, k, k - 1, nc);
      q = FMAT(h, k + 1, k - 1, nc);
      r = notlast ? FMAT(h, k + 2, k - 1, nc) : 0.0;

      x = fabs(p) + fabs(q) + fabs(r);

      if (x == 0) continue; /* FIXME????? */

      p /= x;
      q /= x;
      r /= x;
    }

    s = sqrt(p * p + q * q + r * r);

    if (p < 0) s = -s;

    if (k != m) {
      FMAT(h, k, k - 1, nc) = -s * x;
    } else if (e != m) {
      FMAT(h, k, k - 1, nc) *= -1;
    }

    p += s;
    x = p / s;
    y = q / s;
    z = r / s;
    q /= p;
    r /= p;

    /* do row modifications */

    for (j = k; j <= n; j++) {
      p = FMAT(h, k, j, nc) + q * FMAT(h, k + 1, j, nc);

      if (notlast) {
        p += r * FMAT(h, k + 2, j, nc);
        FMAT(h, k + 2, j, nc) -= p * z;
      }

      FMAT(h, k + 1, j, nc) -= p * y;
      FMAT(h, k, j, nc) -= p * x;
    }

    j = (k + 3 < n) ? (k + 3) : n;

    /* do column modifications */

    for (i = e; i <= j; i++) {
      p = x * FMAT(h, i, k, nc) + y * FMAT(h, i, k + 1, nc);

      if (notlast) {
        p += z * FMAT(h, i, k + 2, nc);
        FMAT(h, i, k + 2, nc) -= p * r;
      }
      FMAT(h, i, k + 1, nc) -= p * q;
      FMAT(h, i, k, nc) -= p;
    }
  }

  goto next_iteration;
}

static void set_companion_matrix(const double *a, size_t nc, double *m) {
  size_t i, j;

  for (i = 0; i < nc; i++)
    for (j = 0; j < nc; j++) MAT(m, i, j, nc) = 0.0;

  for (i = 1; i < nc; i++) MAT(m, i, i - 1, nc) = 1.0;

  for (i = 0; i < nc; i++) MAT(m, i, nc - 1, nc) = -a[i] / a[nc];
}

static int gsl_poly_complex_solve(const double *a,
                                  size_t n,
                                  gsl_poly_complex_workspace *w,
                                  gsl_complex_packed_ptr z) {
  int status;
  double *m;

  if (n == 0) {
    cerr << "number of terms must be a positive integer" << endl;

    return (GSL_FAILURE);
  }

  if (n == 1) {
    cerr << "cannot solve for only one term" << endl;

    return (GSL_FAILURE);
  }

  if (a[n - 1] == 0) {
    cerr << "leading term of polynomial must be non-zero" << endl;

    return (GSL_FAILURE);
  }

  if (w->nc != n - 1) {
    cerr << "size of workspace does not match polynomial" << endl;

    return (GSL_FAILURE);
  }

  m = w->matrix;

  set_companion_matrix(a, n - 1, m);

  balance_companion_matrix(m, n - 1);

  status = qr_companion(m, n - 1, z);

  if (status) {
    //cerr << "root solving qr method failed to converge" << endl;

    return (GSL_FAILURE);
  }

  return GSL_SUCCESS;
}
