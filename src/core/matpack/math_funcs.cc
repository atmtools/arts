/*****************************************************************************
 ***  File description 
 *****************************************************************************/

/*!
   \file   math_funcs.cc
   \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   \date   2000-09-18 

   Contains basic mathematical functions.
*/

/*****************************************************************************
 *** External declarations
 *****************************************************************************/

#include "math_funcs.h"

#include <array.h>
#include <arts_constants.h>
#include <arts_conversions.h>
#include <debug.h>
#include <mystring.h>

#include <cmath>

#include "matpack_mdspan_helpers_reduce.h"

inline constexpr Numeric DEG2RAD = Conversion::deg2rad(1);
inline constexpr Numeric PI      = Constant::pi;

/*****************************************************************************
 *** The functions (in alphabetical order)
 *****************************************************************************/

//! fac
/*!
    Calculates the factorial.

    The function asserts that n must be >= 0

    \return      The factorial
    \param   n   Nominator

    \author Oliver Lemke
    \date   2003-08-15
*/
Numeric fac(const Index n) {
  Numeric sum;

  if (n == 0) return (1.0);

  sum = 1.0;
  for (Index i = 1; i <= n; i++) sum *= Numeric(i);

  return (sum);
}

//! integer_div
/*! 
    Performs an integer division.

    The function asserts that the reminder of the division x/y is 0.

    \return      The quotient
    \param   x   Nominator
    \param   y   Denominator

    \author Patrick Eriksson 
    \date   2002-08-11
*/
Index integer_div(const Index& x, const Index& y) {
  const auto retval = std::div(x, y);
  assert(retval.rem == 0);
  return retval.quot;
}

//! last
/*! 
    Returns the last value of a vector.

    \return      The last value of x.
    \param   x   A vector.

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric last(ConstVectorView x) {
  assert(x.size() > 0);
  return x[x.size() - 1];
}

//! last
/*! 
    Returns the last value of an index array.

    \return      The last value of x.
    \param   x   An index array.

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Index last(const ArrayOfIndex& x) {
  assert(x.size() > 0);
  return x[x.size() - 1];
}

//! linspace
/*! 
    Linearly spaced vector with specified spacing. 

    The first element of x is always start. The next value is start+step etc.
    Note that the last value can deviate from stop.
    The step can be both positive and negative. 
    (in Matlab notation: start:step:stop)

    Size of result is adjusted within this function!

    \param    x       Output: linearly spaced vector
    \param    start   first value in x
    \param    stop    last value of x <= stop
    \param    step    distance between values in x

    \author Patrick Eriksson
    \date   2000-06-27
*/
void linspace(Vector& x,
              const Numeric start,
              const Numeric stop,
              const Numeric step) {
  Index n = (Index)floor((stop - start) / step) + 1;
  if (n < 1) n = 1;
  x.resize(n);
  for (Index i = 0; i < n; i++) x[i] = start + (double)i * step;
}

Vector linspace(const Numeric start, const Numeric stop, const Numeric step) {
  Vector x;
  linspace(x, start, stop, step);
  return x;
}

//! nlinspace
/*! 
    Linearly spaced vector with specified length. 

    Returns a vector equally and linearly spaced between start and stop 
    of length n. (equals the Matlab function linspace)

    The length must be > 1.

    \param    x       Output: linearly spaced vector
    \param    start   first value in x
    \param    stop    last value of x <= stop
    \param    n       length of x

    \author Patrick Eriksson
    \date   2000-06-27
*/
void nlinspace(Vector& x,
               const Numeric start,
               const Numeric stop,
               const Index n) {
  x.resize(n);
  if (x.empty()) return;

  Numeric step = (stop - start) / ((double)n - 1);
  for (Index i = 0; i < n - 1; i++) x[i] = start + (double)i * step;
  x[n - 1] = stop;
}
void nlinspace(VectorView x,
               const Numeric start,
               const Numeric stop,
               const Index n) {
  Numeric step = (stop - start) / ((double)n - 1);
  for (Index i = 0; i < n - 1; i++) x[i] = start + (double)i * step;
  x[n - 1] = stop;
}

Vector nlinspace(const Numeric start, const Numeric stop, const Index step) {
  Vector x;
  nlinspace(x, start, stop, step);
  return x;
}

//! nlogspace
/*! 
    Logarithmically spaced vector with specified length. 

    Returns a vector logarithmically spaced vector between start and 
    stop of length n (equals the Matlab function logspace)

    The length must be > 1.

    \param    x       Output: logarithmically spaced vector
    \param    start   first value in x
    \param    stop    last value of x <= stop
    \param    n       length of x

    \author Patrick Eriksson
    \date   2000-06-27
*/
void nlogspace(Vector& x,
               const Numeric start,
               const Numeric stop,
               const Index n) {
  // Number of points must be greater than 1:
  assert(1 < n);
  // Only positive numbers are allowed for start and stop:
  assert(0 < start);
  assert(0 < stop);

  x.resize(n);
  Numeric a    = log(start);
  Numeric step = (log(stop) - a) / ((double)n - 1);
  x[0]         = start;
  for (Index i = 1; i < n - 1; i++) x[i] = std::exp(a + (double)i * step);
  x[n - 1] = stop;
}

Vector nlogspace(const Numeric start, const Numeric stop, const Index step) {
  Vector x;
  nlogspace(x, start, stop, step);
  return x;
}

Vector binary_grid(const Numeric x0, const Numeric xn, const Numeric dx) {
  assert(dx > 0);

  Numeric dX = xn - x0;
  Size N     = 2;
  Size dN    = 1;

  while (dX > dx) {
    dX *= 0.5;
    N  += dN;
    dN += dN;
  }

  return nlinspace(x0, xn, N);
}

//! trapz
/*! 
    Integration by the basic trapezoidal rule

    \param x   Grid values
    \param y   Integrand

    \return The resulting integral

    \author Patrick Eriksson
    \date 2022-03-06
*/
Numeric trapz(ConstVectorView x, ConstVectorView y) {
  const Size n = x.size();
  assert(y.size() == n);
  Numeric sum = 0.0;
  for (Size i = 1; i < n; ++i) sum += (x[i] - x[i - 1]) * (y[i] + y[i - 1]);
  return sum / 2.0;
}

//! cumsum
/*! 
    Cumulative sum of a vector

    csum[0] becomes x[0] and last value in csum equals the full sum of x

    \param csum  Vector with cumulative values
    \param x     Input vector

    \author Patrick Eriksson
    \date 2022-12-22
*/
void cumsum(VectorView csum, ConstVectorView x) {
  const Size n = x.size();
  assert(csum.size() == n);
  csum[0] = x[0];
  for (Size i = 1; i < n; ++i) csum[i] = csum[i - 1] + x[i];
}

//! AngIntegrate_trapezoid
/*! 
    Performs an integration of a matrix over all directions defined in angular
    grids using the trapezoidal integration method.

    \param Integrand The Matrix to be integrated
    \param zen_grid   The zenith angle grid 
    \param azi_grid   The azimuth angle grid 
    
    \return The resulting integral
*/
Numeric AngIntegrate_trapezoid(ConstMatrixView Integrand,
                               ConstVectorView zen_grid,
                               ConstVectorView azi_grid) {
  Index n = zen_grid.size();
  Index m = azi_grid.size();
  Vector res1(n);
  assert((Integrand.shape() == std::array{n, m}));

  for (Index i = 0; i < n; ++i) {
    res1[i] = 0.0;

    for (Index j = 0; j < m - 1; ++j) {
      res1[i] += 0.5 * DEG2RAD * (Integrand[i, j] + Integrand[i, j + 1]) *
                 (azi_grid[j + 1] - azi_grid[j]) *
                 std::sin(zen_grid[i] * DEG2RAD);
    }
  }
  Numeric res = 0.0;
  for (Index i = 0; i < n - 1; ++i) {
    res += 0.5 * DEG2RAD * (res1[i] + res1[i + 1]) *
           (zen_grid[i + 1] - zen_grid[i]);
  }

  return res;
}

//! AngIntegrate_trapezoid_opti
/*! 
    Performs an integration of a matrix over all directions defined in angular
    grids using the trapezoidal integration method.

    In addition to the "old fashined" integration method, it checks whether
    the stepsize is constant. If it is, it uses a faster method, if not, it
    uses the old one.

    \param Integrand Input : The Matrix to be integrated
    \param zen_grid Input : The zenith angle grid 
    \param azi_grid Input : The azimuth angle grid
    \param grid_stepsize Input : stepsize of the grid
    
    \return The resulting integral

    \author Claas Teichmann <claas@sat.physik.uni-bremen.de>
    \date 2003/05/28
*/
Numeric AngIntegrate_trapezoid_opti(ConstMatrixView Integrand,
                                    ConstVectorView zen_grid,
                                    ConstVectorView azi_grid,
                                    ConstVectorView grid_stepsize) {
  Numeric res = 0;
  if ((grid_stepsize[0] > 0) && (grid_stepsize[1] > 0)) {
    Index n             = zen_grid.size();
    Index m             = azi_grid.size();
    Numeric stepsize_za = grid_stepsize[0];
    Numeric stepsize_aa = grid_stepsize[1];
    Vector res1(n);
    assert((Integrand.shape() == std::array{n, m}));

    Numeric temp = 0.0;

    for (Index i = 0; i < n; ++i) {
      temp = Integrand[i, 0];
      for (Index j = 1; j < m - 1; j++) {
        temp += Integrand[i, j] * 2;
      }
      temp    += Integrand[i, m - 1];
      temp    *= 0.5 * DEG2RAD * stepsize_aa * std::sin(zen_grid[i] * DEG2RAD);
      res1[i]  = temp;
    }

    res = res1[0];
    for (Index i = 1; i < n - 1; i++) {
      res += res1[i] * 2;
    }
    res += res1[n - 1];
    res *= 0.5 * DEG2RAD * stepsize_za;
  } else {
    res = AngIntegrate_trapezoid(Integrand, zen_grid, azi_grid);
  }

  return res;
}

//! AngIntegrate_trapezoid
/*! 
    Performs an integration of a matrix over all directions defined in angular
    grids using the trapezoidal integration method.
    The integrand is independant of the azimuth angle. The integration over
    the azimuth angle gives a 2*PI

    \param Integrand Input : The vector to be integrated
    \param zen_grid Input : The zenith angle grid 

    \author Claas Teichmann
    \date   2003-05-13
    
    \return The resulting integral
*/
Numeric AngIntegrate_trapezoid(ConstVectorView Integrand,
                               ConstVectorView zen_grid) {
  Index n = zen_grid.size();
  assert((Integrand.shape() == std::array{n}));

  Numeric res = 0.0;
  for (Index i = 0; i < n - 1; ++i) {
    // in this place 0.5 * 2 * PI is calculated:
    res += PI * DEG2RAD *
           (Integrand[i] * std::sin(zen_grid[i] * DEG2RAD) +
            Integrand[i + 1] * std::sin(zen_grid[i + 1] * DEG2RAD)) *
           (zen_grid[i + 1] - zen_grid[i]);
  }

  return res;
}

//! sign
/*! 
    Returns the sign of a numeric value.

    The function returns 1 if the value is greater than zero, 0 if it 
    equals zero and -1 if it is less than zero.

    \return      The sign of x (see above).
    \param   x   A Numeric.

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric sign(const Numeric& x) {
  if (x < 0) return -1.0;
  if (x == 0) return 0.0;
  return 1.0;
}

//! sign
/*! 
    Returns the sign of an Index value.

    The function returns 1 if the value is greater than zero, 0 if it 
    equals zero and -1 if it is less than zero.

    \return      The sign of x (see above).
    \param   x   An Index.

    \author Patrick Eriksson 
    \date   2022-09-29
*/
Index sign(const Index& x) { return (x > 0) - (x < 0); }

Numeric min_geq(const Numeric n1, const Numeric n2, const Numeric limit) {
  if (n1 < limit) {
    if (n2 < limit) return limit - 1;
    return n2;
  }

  if (n2 < limit) {
    if (n1 < limit) return limit - 1;
    return n1;
  }

  return std::min(n1, n2);
}

Index n_int_between(const Numeric x, const Numeric y) {
  Index n;
  if (y > x) {
    n = Index(std::ceil(y)) - Index(std::floor(x)) - 1;
    return n > 0 ? n : 0;
  }
  n = Index(std::ceil(x)) - Index(std::floor(y)) - 1;
  return n > 0 ? -n : 0;
}

Index int_at_step(const Numeric gp, const Index step) {
  assert(step != 0);
  if (step > 0) return Index(std::floor(gp)) + step;
  return Index(std::ceil(gp)) + step;
}

/*! Modified gamma distribution
 *  
 *  Uses all four free parameters (n0, mu, la, ga) to calculate
 *    psd(D) = n0 * D^mu * std::exp( -la * x^ga )
 *  
 *  Reference: Eq 1 of Petty & Huang, JAS, (2011).

    \param psd       Particle number density per x-interval. Sizing of vector
                     should be done before calling the function.
    \param x       Mass or size.
    \param n0      See above.
    \param mu      See above.
    \param la      See above.
    \param ga      See above.
  
  \author Jana Mendrok, Patrick Eriksson
  \date 2017-06-07

*/
void mgd(VectorView psd,
         const Vector& x,
         const Numeric& n0,
         const Numeric& mu,
         const Numeric& la,
         const Numeric& ga) {
  const Size nx = x.size();

  assert(psd.size() == nx);

  if (ga == 1) {
    if (mu == 0) {
      // std::Exponential distribution
      for (Size ix = 0; ix < nx; ix++) {
        const Numeric eterm = std::exp(-la * x[ix]);
        psd[ix]             = n0 * eterm;
      }
    } else {
      ARTS_USER_ERROR_IF(mu > 10,
                         "Given mu is {}\n"
                         "Seems unreasonable. Have you mixed up the inputs?",
                         mu)
      // Gamma distribution
      for (Size ix = 0; ix < nx; ix++) {
        const Numeric eterm = std::exp(-la * x[ix]);
        const Numeric xterm = nonstd::pow(x[ix], mu);
        psd[ix]             = n0 * xterm * eterm;
        psd[ix] = n0 * nonstd::pow(x[ix], mu) * std::exp(-la * x[ix]);
      }
    }
  } else {
    // Complete MGD
    ARTS_USER_ERROR_IF(mu > 10,
                       "Given mu is {}\n"
                       "Seems unreasonable. Have you mixed up the inputs?",
                       mu)
    ARTS_USER_ERROR_IF(ga > 10,
                       "Given gamma is {}\n"
                       "Seems unreasonable. Have you mixed up the inputs?",
                       ga)
    for (Size ix = 0; ix < nx; ix++) {
      const Numeric pterm = nonstd::pow(x[ix], ga);
      const Numeric eterm = std::exp(-la * pterm);
      const Numeric xterm = nonstd::pow(x[ix], mu);
      psd[ix]             = n0 * xterm * eterm;
    }
  }
}

/*! Modified gamma distribution, and derivatives
 *  
 *  As mgd, but this version can also return the derivate of psd with respect 
 *  to the four parameters.   

    \param psd       Particle number density per x-interval. Sizing of vector
                     should be done before calling the function.
    \param jac_data  Container for returning jacobian data. Shall be a matrix
                     with four rows, where the rows match n0, mu, la and ga.  
                     Number of columns same as length of psd.
    \param x       Mass or size.
    \param n0      See above.
    \param mu      See above.
    \param la      See above.
    \param ga      See above.
    \param do_n0_jac  Flag to actually calculate d_psd/d_n0
    \param do_mu_jac  Flag to actually calculate d_psd/d_mu
    \param do_la_jac  Flag to actually calculate d_psd/d_la
    \param do_ga_jac  Flag to actually calculate d_psd/d_ga
  
  \author Patrick Eriksson
  \date 2017-06-07

*/
void mgd_with_derivatives(VectorView psd,
                          MatrixView jac_data,
                          const Vector& x,
                          const Numeric& n0,
                          const Numeric& mu,
                          const Numeric& la,
                          const Numeric& ga,
                          const bool& do_n0_jac,
                          const bool& do_mu_jac,
                          const bool& do_la_jac,
                          const bool& do_ga_jac) {
  const Size nx = x.size();

  assert(psd.size() == nx);
  assert(static_cast<Size>(jac_data.nrows()) == 4);
  assert(static_cast<Size>(jac_data.ncols()) == nx);

  if (ga == 1 && !do_ga_jac) {
    if (mu == 0 && !do_mu_jac) {
      // std::Exponential distribution
      for (Size ix = 0; ix < nx; ix++) {
        const Numeric eterm = std::exp(-la * x[ix]);
        psd[ix]             = n0 * eterm;
        if (do_n0_jac) {
          jac_data[0, ix] = eterm;
        }
        if (do_la_jac) {
          jac_data[2, ix] = -x[ix] * psd[ix];
        }
      }
    } else {
      ARTS_USER_ERROR_IF(mu > 10,
                         "Given mu is {}\n"
                         "Seems unreasonable. Have you mixed up the inputs?",
                         mu)
      // Gamma distribution
      for (Size ix = 0; ix < nx; ix++) {
        const Numeric eterm = std::exp(-la * x[ix]);
        const Numeric xterm = nonstd::pow(x[ix], mu);
        psd[ix]             = n0 * xterm * eterm;
        if (do_n0_jac) {
          jac_data[0, ix] = xterm * eterm;
        }
        if (do_mu_jac) {
          jac_data[1, ix] = std::log(x[ix]) * psd[ix];
        }
        if (do_la_jac) {
          jac_data[2, ix] = -x[ix] * psd[ix];
        }
        psd[ix] = n0 * nonstd::pow(x[ix], mu) * std::exp(-la * x[ix]);
      }
    }
  } else {
    // Complete MGD
    ARTS_USER_ERROR_IF(mu > 10,
                       "Given mu is {}\n"
                       "Seems unreasonable. Have you mixed up the inputs?",
                       mu)
    ARTS_USER_ERROR_IF(ga > 10,
                       "Given gamma is {}\n"
                       "Seems unreasonable. Have you mixed up the inputs?",
                       ga)
    for (Size ix = 0; ix < nx; ix++) {
      const Numeric pterm = nonstd::pow(x[ix], ga);
      const Numeric eterm = std::exp(-la * pterm);
      const Numeric xterm = nonstd::pow(x[ix], mu);
      psd[ix]             = n0 * xterm * eterm;
      if (do_n0_jac) {
        jac_data[0, ix] = xterm * eterm;
      }
      if (do_mu_jac) {
        jac_data[1, ix] = std::log(x[ix]) * psd[ix];
      }
      if (do_la_jac) {
        jac_data[2, ix] = -pterm * psd[ix];
      }
      if (do_ga_jac) {
        jac_data[3, ix] = -la * pterm * std::log(x[ix]) * psd[ix];
      }
    }
  }
}

void delanoe_shape_with_derivative(VectorView psd,
                                   MatrixView jac_data,
                                   const Vector& x,
                                   const Numeric& alpha,
                                   const Numeric& beta) {
  Numeric f_c  = std::tgamma(4.0) / 256.0;
  f_c         *= nonstd::pow(std::tgamma((alpha + 5.0) / beta), 4 + alpha);
  f_c         /= nonstd::pow(std::tgamma((alpha + 4.0) / beta), 5 + alpha);

  Numeric f_d  = std::tgamma((alpha + 5.0) / beta);
  f_d         /= std::tgamma((alpha + 4.0) / beta);

  for (Size i = 0; i < x.size(); ++i) {
    Numeric xi = x[i];
    psd[i]     = beta * f_c * nonstd::pow(xi, alpha) *
             std::exp(-nonstd::pow(f_d * xi, beta));
    jac_data[0, i] =
        psd[i] * (alpha / xi - beta * f_d * nonstd::pow(f_d * xi, beta - 1.0));
  }
}

//! Generalized Modified Gamma Distribution
/*! Returns number density per unit of 'x' as function of 'x'.
 
 \return  dN Number density as function of x.
 \param   x       Numeric
 \param   N0      Numeric, Scaling parameter
 \param   Lambda  Numeric, Shape parameter
 \param   mu      Numeric, Shape parameter
 \param   gamma   Numeric, Shape parameter
 
 \author Manfred Brath
 \date   2015-01-19
 */

Numeric mod_gamma_dist(
    Numeric x, Numeric N0, Numeric Lambda, Numeric mu, Numeric gamma) {
  Numeric dN;

  if (x > 0. && N0 > 0. && Lambda > 0. && (mu + 1) / gamma > 0.) {
    //Distribution function
    dN = N0 * nonstd::pow(x, mu) * std::exp(-Lambda * nonstd::pow(x, gamma));

    return dN;
  }

  ARTS_USER_ERROR(
      "At least one argument is zero or negative.\n"
      "Modified gamma distribution can not be calculated.\n"
      "x      = {}"
      "\n"
      "N0     = {}"
      "\n"
      "lambda = {}"
      "\n"
      "mu     = {}"
      "\n"
      "gamma  = {}\n",
      x,
      N0,
      Lambda,
      mu,
      gamma)
}

//! unitl
/*!
    Normalises a vector to have unit length.

    The standard Euclidean norm is used (2-norm).

    param    x   In/Out: A vector.

    \author Patrick Eriksson
    \date   2012-02-12
*/
void unitl(Vector& x) {
  assert(x.size() > 0);

  const Numeric l = std::sqrt(dot(x, x));
  for (Size i = 0; i < x.size(); i++) x[i] /= l;
}

//! flat
/*!
    Flattens a matrix to a vector

    The matrix is read from front, i.e. rows are looped first. 
    In Matlab this equals x=X(:).

    \param[out] x   The vector. Should already be sized
    \param[in]  X   The matrix.

    \author Patrick Eriksson
    \date   2015-09-09
*/
void flat(VectorView x, ConstMatrixView X) {
  assert(static_cast<Index>(x.size()) == X.nrows() * X.ncols());

  Index i = 0;

  for (Index c = 0; c < X.ncols(); c++) {
    for (Index r = 0; r < X.nrows(); r++) {
      x[i]  = X[r, c];
      i    += 1;
    }
  }
}

//! flat
/*!
    Converts Tensor3 to a vector

    The matrix is read from front, i.e. pages are looped first, followed by rows. 
    In Matlab this equals x=X(:).

    \param[out] x   The vector. Should already be sized
    \param[in]  X   The tensor.

    \author Patrick Eriksson
    \date   2015-09-09
*/
void flat(VectorView x, ConstTensor3View X) {
  assert(static_cast<Index>(x.size()) == X.nrows() * X.ncols() * X.npages());

  Index i = 0;

  for (Index c = 0; c < X.ncols(); c++) {
    for (Index r = 0; r < X.nrows(); r++) {
      for (Index p = 0; p < X.npages(); p++) {
        x[i]  = X[p, r, c];
        i    += 1;
      }
    }
  }
}

//! reshape
/*!
    Converts vector to Tensor3

    The tensor is filled from front, i.e. pages are looped first, followed by rows. 
    In Matlab this equals X = reshape( x, [ X.npages(), X.nrows(), X.ncols() ]

    \param[out] X   The tensor. Should already be sized
    \param[in]  x   The vector.

    \author Patrick Eriksson
    \date   2015-09-10
*/
void reshape(Tensor3View X, ConstVectorView x) {
  assert(static_cast<Index>(x.size()) == X.nrows() * X.ncols() * X.npages());

  Index i = 0;

  for (Index c = 0; c < X.ncols(); c++) {
    for (Index r = 0; r < X.nrows(); r++) {
      for (Index p = 0; p < X.npages(); p++) {
        X[p, r, c]  = x[i];
        i          += 1;
      }
    }
  }
}

//! reshape
/*!
    Converts vector to Matrix

    The matrix is filled from front, i.e. rows are looped first, followed by cols. 
    In Matlab this equals X = reshape( x, [ X.nrows(), X.ncols() ]

    \param[out] X   The matrix. Should already be sized
    \param[in]  x   The vector.

    \author Patrick Eriksson
    \date   2015-09-10
*/
void reshape(MatrixView X, ConstVectorView x) {
  assert(static_cast<Index>(x.size()) == X.nrows() * X.ncols());

  Index i = 0;

  for (Index c = 0; c < X.ncols(); c++) {
    for (Index r = 0; r < X.nrows(); r++) {
      X[r, c]  = x[i];
      i       += 1;
    }
  }
}

//! calculate_weights_linear
/*!
  Function to set the evaluation points and the corresponding weights
  for numerical integration on the domain from [-1,1]

  Parameters:

  \param[out] x evaluation points
  \paran[out] w integration weights
  \param[in] nph       number of evaluation points per hemisphere
*/
void calculate_weights_linear(Vector& x, Vector& w, const Index nph) {
  Index N = nph * 2;

  //directions in total
  nlinspace(x, -1, 1, N);

  //allocate
  w.resize(x.size());

  // calculate weights
  w[0] = (x[1] - x[0]) / 2.;

  for (Index i = 1; i < nph * 2 - 1; i++) {
    w[i] = (x[i + 1] - x[i - 1]) / 2.;
  }
  w[x.size() - 1] = (x[x.size() - 1] - x[x.size() - 2]) / 2.;
}

void calculate_int_weights_arbitrary_grid(Vector& w, const Vector& x) {
  ARTS_USER_ERROR_IF(x.size() < 1, "Grid needs at least 2 points.");

  Index N = x.size();

  w.resize(N);

  w[0] = (x[1] - x[0]) / 2;
  if (N > 2) {
    for (Index i = 1; i < N - 1; i++) {
      w[i] = (x[i + 1] - x[i - 1]) / 2;
    }
  }
  w[N - 1] = (x[N - 1] - x[N - 2]) / 2;
}
