/* Copyright (C) 2003 Oliver Lemke  <olemke@uni-bremen.de>

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
   USA.
*/



////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
  \file   legendre.cc

  Contains the code to calculate Legendre polynomials.

  \author Oliver Lemke
  \date 2003-08-14
  */

#include "legendre.h"

#ifdef HAVE_SSTREAM
#include <sstream>
#else
#include "sstream.h"
#endif
#include <stdexcept>

#include "math_funcs.h"


//! legendre_poly
/*!
    Returns the associated Legendre polynomial Plm(x).

    The input parameters must fulfill the following conditions:
    0 <= m <= l and |x| <= 1

    The code is based on the Numerical recipes. Results were compared
    to the Legendre calculations from the GNU Scientific library and found
    to be identical.

    \return      Plm
    \param   l   Index
    \param   m   Index
    \param   x   Value

    \author Oliver Lemke
    \date   2003-08-14
*/
Numeric
legendre_poly (Index l, Index m, Numeric x)
{
  Numeric pmm;
  Numeric result;

  if (m < 0 || m > l || abs (x) > 1.0)
    {
      ostringstream os;
      os << "legendre_poly: Condition 0 <= m <= l && -1 < x < 1 failed"
        << endl << "  l = " << l << "  m = " << m << "  x = " << x << endl;
      throw runtime_error (os.str ());
    }

  pmm = 1.0;
  if (m > 0)
    {
      Numeric fact, somx2;

      somx2 = sqrt ((1.0 - x) * (1.0 + x));
      fact = 1.0;
      for (Index i = 1; i <= m; i++)
        {
          pmm *= -fact * somx2;
          fact += 2.0;
        }
    }

  if (l == m)
    result = pmm;
  else
    {
      Numeric pmmp1;

      pmmp1 = x * (2*m + 1) * pmm;
      if (l == (m+1))
        result = pmmp1;
      else
        {
          for (Index ll = m+2; ll <= l; ll++)
            {
              result = (x * (2*ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
              pmm = pmmp1;
              pmmp1 = result;
            }
        }
    }

  return (result);
}


//! legendre_poly_norm
/*!
    Returns the normalized associated Legendre polynomial Plm(x).

    The input parameters must fulfill the following conditions:
    0 <= m <= l and |x| <= 1

    The code is based on the Numerical recipes. Results were compared
    to the Legendre calculations from the GNU Scientific library and found
    to be identical.

    \return      Plm
    \param   l   Index
    \param   m   Index
    \param   x   Value

    \author Oliver Lemke
    \date   2003-08-15
*/
Numeric
legendre_poly_norm (Index l, Index m, Numeric x)
{
  return (sqrt (2.0 * fac (1 - m) / fac (1 + m))
          * legendre_poly (l, m, x));
}


//! legendre_poly_deriv
/*!
    Returns the derivative of the associated Legendre polynomial Plm(x).

    The input parameters must fulfill the following conditions:
    0 <= m <= l and |x| < 1

    \return      dPlm
    \param   l   Index
    \param   m   Index
    \param   x   Value

    \author Oliver Lemke
    \date   2003-08-18
*/
Numeric
legendre_poly_deriv (Index l, Index m, Numeric x)
{
  assert (x != 1.);
  if (x == 1.)
    {
      ostringstream os;
      os << "legendre_poly_deriv: Condition x != 1 failed"
        << endl << "  x = " << x << endl;
      throw runtime_error (os.str ());
    }
  return (((l + 1) * (l + m) + legendre_poly (l, m-1, x)
           - l * (l - m + 1) * legendre_poly (l+1, m, x))
          / ((2*l + 1) * (1 - x * x)));
}

