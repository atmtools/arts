/* Copyright (C) 2003-2012 Oliver Lemke <olemke@core-dump.info>

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

#include <cmath>
#include <sstream>

#include "legendre.h"
#include "exceptions.h"
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
  Numeric result = 0.;

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

      pmmp1 = x * (Numeric)(2*m + 1) * pmm;
      if (l == (m+1))
        result = pmmp1;
      else
        {
          for (Index ll = m+2; ll <= l; ll++)
            {
              result = (x * (Numeric)(2*ll - 1) * pmmp1 - (Numeric)(ll + m - 1) * pmm) / (Numeric)(ll - m);
              pmm = pmmp1;
              pmmp1 = result;
            }
        }
    }

  return (result);
}


//! legendre_poly_norm_schmidt
/*!
    Returns the Schmidt quasi-normalized associated Legendre polynomial Plm(x).

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
legendre_poly_norm_schmidt (Index l, Index m, Numeric x)
{
  Numeric result;

  if (m != 0)

    {
      result = ((sqrt (2.0 * fac (l - m) / fac (l + m))
                 * legendre_poly (l, m, x)));
    }
  else
    {
      result = (legendre_poly (l, m, x));
    }

  return(result);

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
  Numeric result;

  if (l == 1)
    {
      if (m == 0)
        {
          result = 1;
        }
      else if (m == 1)
        {
          result = x / sqrt(1 - x*x);
        }
      else
        {
          ostringstream os;
          os << "legendre_poly_deriv: "
            << "Condition l == 1 && (m == 0 || m == 1) failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }
  else if ( m < l )
    {
      try
        {
          result = ((Numeric)(l + m) * legendre_poly (l-1, m, x) - 
                    (Numeric)l * x * legendre_poly (l, m, x)) /
            (1 - x * x);
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "legendre_poly_deriv: "
            << "Condition m < l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }
  else
    {
      try
        {
          result = (Numeric)m * x * legendre_poly (l, m, x) / (1 - x * x)
            + (Numeric)((l + m) * (l - m + 1)) * legendre_poly (l, m - 1, x)
            / sqrt (1 - x * x);
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "legendre_poly_norm_schmidt_deriv: "
            << "Condition m = l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }
 return (result);
}

//! legendre_poly_norm_schmidt_deriv
/*!
  Returns the derivative of the Schmidt quasi-normalized associated Legendre polynomial Plm(x).

  The input parameters must fulfill the following conditions:
  0 <= m <= l and |x| < 1

  \return      dPlm
  \param   l   Index
  \param   m   Index
  \param   x   Value

  \author Nikolay Koulev
  \date   2003-08-18
*/
Numeric
legendre_poly_norm_schmidt_deriv (Index l, Index m, Numeric x)
{
  assert (x != 1.);
  if (x == 1.)
    {
      ostringstream os;
      os << "legendre_poly_norm_schmidt_deriv: Condition x != 1 failed"
        << endl << "  x = " << x << endl;
      throw runtime_error (os.str ());
    }

  Numeric result;

  if (l == 1)
    {
      if (m == 0)
        {
          result = sqrt (2.0 * fac (1 - m) / fac (1 + m));
        }
      else if (m == 1)
        {
          result = sqrt (2.0 * fac (1 - m) / fac (1 + m)) * x / sqrt(1 - x*x);
        }
      else
        {
          ostringstream os;
          os << "legendre_poly_norm_schmidt_deriv: "
            << "Condition l == 1 && (m == 0 || m == 1) failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }

    }
  else if ( m < l )
    {
      try
        {
          result = ((Numeric)(l + m) * legendre_poly_norm_schmidt (l-1, m, x) - 
                    (Numeric)l * x * legendre_poly_norm_schmidt (l, m, x)) /
            (1 - x * x);
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "legendre_poly_norm_schmidt_deriv: "
            << "Condition m < l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }
  else
    {
      try
        {
          result = (Numeric)m * x * legendre_poly_norm_schmidt (l, m, x) / (1 - x * x)
            + (Numeric)((l + m) * (l - m + 1)) * legendre_poly_norm_schmidt (l, m - 1, x)
            / sqrt (1 - x * x);
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "legendre_poly_norm_schmidt_deriv: "
            << "Condition m = l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }

  return (result);
}

//! g_legendre_poly
/*!
  Returns the associated Legendre polynomial Plm(x) without 
  the factor (-1)^m.

  The input parameters must fulfill the following conditions:
  0 <= m <= l and |x| <= 1

  The code is based on the Numerical recipes. Results were compared
  to the Legendre calculations from the GNU Scientific library and found
  to be identical.

  \return      Plm
  \param   l   Index
  \param   m   Index
  \param   x   Value

  \author Nikolay Koulev
  \date   2003-09-02
  */
Numeric
g_legendre_poly (Index l, Index m, Numeric x)
{
  Numeric pmm;
  Numeric result = 0.;

  if (m < 0 || m > l || abs (x) > 1.0)
    {
      ostringstream os;
      os << "g_legendre_poly: Condition 0 <= m <= l && -1 < x < 1 failed"
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
          pmm *= fact * somx2;
          fact += 2.0;
        }
    }

  if (l == m)
    result = pmm;
  else
    {
      Numeric pmmp1;

      pmmp1 = x * (Numeric)(2*m + 1) * pmm;
      if (l == (m+1))
        result = pmmp1;
      else
        {
          for (Index ll = m+2; ll <= l; ll++)
            {
              result = (x * (Numeric)(2*ll - 1) * pmmp1 - (Numeric)(ll + m - 1) * pmm) / (Numeric)(ll - m);
              pmm = pmmp1;
              pmmp1 = result;
            }
        }
    }

  return (result);
}


//! g_legendre_poly_norm_schmidt
/*!
  Returns the Schmidt quasi-normalized associated Legendre polynomial Plm(x)) without 
  the factor (-1)^m..

  The input parameters must fulfill the following conditions:
  0 <= m <= l and |x| <= 1

  The code is based on the Numerical recipes. Results were compared
  to the Legendre calculations from the GNU Scientific library and found
  to be identical.

  \return      Plm
  \param   l   Index
  \param   m   Index
  \param   x   Value

  \author Nikolay koulev
  \date   2003-09-02
  */
Numeric
g_legendre_poly_norm_schmidt (Index l, Index m, Numeric x)
{
  Numeric result;

  if (m != 0)

    {
      result = ((sqrt (2.0 * fac (l - m) / fac (l + m))
                 * g_legendre_poly (l, m, x)));
    }
  else
    {
      result = (g_legendre_poly (l, m, x));
    }

  return(result);

}


//! g_legendre_poly_deriv
/*!
  Returns the derivative of the associated Legendre polynomial Plm(x)) without 
  the factor (-1)^m..

  The input parameters must fulfill the following conditions:
  0 <= m <= l and |x| < 1

  \return      dPlm
  \param   l   Index
  \param   m   Index
  \param   x   Value

  \author Nikolay Koulev
  \date   2003-09-02
  */
Numeric
g_legendre_poly_deriv (Index l, Index m, Numeric x)
{
  assert (x != 1.);
  if (x == 1.)
    {
      ostringstream os;
      os << "g_legendre_poly_deriv: Condition x != 1 failed"
        << endl << "  x = " << x << endl;
      throw runtime_error (os.str ());
    }
  Numeric result;

  if (l == 1)
    {
      if (m == 0)
        {
          result = 1;
        }
      else if (m == 1)
        {
          result = x / sqrt(1 - x*x);
        }
      else
        {
          ostringstream os;
          os << "g_legendre_poly_deriv: "
            << "Condition l == 1 && (m == 0 || m == 1) failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }
  else if ( m < l )
    {
      try
        {
          result = ((Numeric)(l + m) * g_legendre_poly (l-1, m, x) - 
                    (Numeric)l * x * g_legendre_poly (l, m, x)) /
            (1 - x * x);
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "g_legendre_poly_deriv: "
            << "Condition m < l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }
  else
    {
      try
        {
          result = - (Numeric)m * x * g_legendre_poly (l, m, x) / (1 - x * x) + 
            (Numeric)((l + m) * (l - m + 1)) * g_legendre_poly (l, m - 1, x) / 
            sqrt (1 - x * x);
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition m = l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }
 return (result);
}

//! g_legendre_poly_norm_schmidt_deriv
/*!
  Returns the derivative of the Schmidt quasi-normalized associated Legendre polynomial Plm(x)) without 
  the factor (-1)^m.

  The input parameters must fulfill the following conditions:
  0 <= m <= l and |x| < 1

  \return      dPlm
  \param   l   Index
  \param   m   Index
  \param   x   Value

  \author Nikolay Koulev
  \date   2003-09-02
*/
Numeric
g_legendre_poly_norm_schmidt_deriv (Index l, Index m, Numeric x)
{
  assert (x != 1.);
  if (x == 1.)
    {
      ostringstream os;
      os << "g_legendre_poly_norm_schmidt_deriv: Condition x != 1 failed"
        << endl << "  x = " << x << endl;
      throw runtime_error (os.str ());
    }

  Numeric result;

  if (l == 1)
    {
      if (m == 0)
        {
          result = 1;
        }
      else if (m == 1)
        {
          result = x / sqrt(1 - x * x);
        }
      else
        {
          ostringstream os;
          os << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition l == 1 && (m == 0 || m == 1) failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }

    }
  else if ( m < l )
    {
      try
        {
          result = sqrt (2.0 * fac (l - m) / fac (l + m)) * 
            ((Numeric)(l + m) * g_legendre_poly(l-1, m, x) - 
             (Numeric)l * x * g_legendre_poly (l, m, x)) / (1 - x * x);
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition m < l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }
  else
    {
      try
        {
          result = sqrt (2.0 * fac (l - m) / fac (l + m)) * 
            ( - (Numeric)m * x * g_legendre_poly (l, m, x) / (1 - x * x) + 
              (Numeric)((l + m) * (l - m + 1)) * g_legendre_poly (l, m - 1, x) / 
              sqrt (1 - x * x));
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition m = l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }

  return (result);
}

//! g_legendre_poly_norm_schmidt_deriv1
/*!
  Returns the derivative of the Schmidt quasi-normalized associated Legendre polynomial Plm(x)) without 
  the factor (-1)^m. Utilizes the simplest recurrence scheme.

  The input parameters must fulfill the following conditions:
  0 <= m <= l and |x| < 1

  \return      dPlm
  \param   l   Index
  \param   m   Index
  \param   x   Value

  \author Nikolay Koulev
  \date   2003-09-02
*/
Numeric
g_legendre_poly_norm_schmidt_deriv1 (Index l, Index m, Numeric x)
{
  assert (x != 1.);
  if (x == 1.)
    {
      ostringstream os;
      os << "g_legendre_poly_norm_schmidt_deriv: Condition x != 1 failed"
        << endl << "  x = " << x << endl;
      throw runtime_error (os.str ());
    }

  Numeric result;

  if (l == 1)
    {
      if (m == 0)
        {
          result = 1;
        }
      else if (m == 1)
        {
          result = x / sqrt(1 - x * x);
        }
      else
        {
          ostringstream os;
          os << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition l == 1 && (m == 0 || m == 1) failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }

    }
  else if ( m <= l - 1 )
    {
      try
        {
          // result = - m * x * g_legendre_poly_norm_schmidt (l, m, x) / (1 - x * x)
          //  +  sqrt((double)(l + m + 1 / l - m - 1)) * g_legendre_poly_norm_schmidt (l, m + 1, x)
          //  / sqrt (1 - x * x);

          result = sqrt (2.0 * fac (l - m) / fac (l + m)) * 
            ( - (Numeric)m * x * g_legendre_poly (l, m, x) / (1 - x * x)  + 
              g_legendre_poly (l, m + 1, x)
              / sqrt (1 - x * x));
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition m <= l - 1 failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }
  else
    {
      try
        {
          result = - sqrt (2.0 * fac (l - m) / fac (l + m)) * 
            (Numeric)m * x * g_legendre_poly (l, m, x) / (1 - x * x);
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition m = l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }

  return (result);
}

//! g_legendre_poly_norm_schmidt_deriv2
/*!
  Returns the derivative of the Schmidt quasi-normalized associated Legendre polynomial Plm(x)) without 
  the factor (-1)^m.

  The input parameters must fulfill the following conditions:
  0 <= m <= l and |x| < 1

  \return      dPlm
  \param   l   Index
  \param   m   Index
  \param   x   Value

  \author Nikolay Koulev
  \date   2003-09-02
*/
Numeric
g_legendre_poly_norm_schmidt_deriv2 (Index l, Index m, Numeric x)
{
  assert (x != 1.);
  if (x == 1.)
    {
      ostringstream os;
      os << "g_legendre_poly_norm_schmidt_deriv: Condition x != 1 failed"
        << endl << "  x = " << x << endl;
      throw runtime_error (os.str ());
    }

  Numeric result;

  if (l == 1)
    {
      if (m == 0)
        {
          result = 1;
        }
      else if (m == 1)
        {
          result = x / sqrt(1 - x * x);
        }
      else
        {
          ostringstream os;
          os << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition l == 1 && (m == 0 || m == 1) failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }

    }
  else if ( m < l )
    {
      try
        {
          result = - sqrt (2.0 * fac (l - m) / fac (l + m)) * 
            ((Numeric)(l + m) * g_legendre_poly (l-1, m, x) - 
             (Numeric)l * x * g_legendre_poly (l, m, x)) /
            (1 - x * x);
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition m < l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }
  else
    {
      try
        {
          result = - (Numeric)m * x * g_legendre_poly_norm_schmidt (l, m, x) / (1 - x * x)
            ;
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition m = l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }

  return (result);
}

//! g_legendre_poly_norm_schmidt_deriv3
/*!
  Returns the derivative of the Schmidt quasi-normalized associated Legendre polynomial Plm(x)) without 
  the factor (-1)^m.

  The input parameters must fulfill the following conditions:
  0 <= m <= l and |x| < 1

  \return      dPlm
  \param   l   Index
  \param   m   Index
  \param   x   Value

  \author Nikolay Koulev
  \date   2003-09-02
  */
Numeric
g_legendre_poly_norm_schmidt_deriv3 (Index l, Index m, Numeric x)
{
  assert (x != 1.);
  if (x == 1.)
    {
      ostringstream os;
      os << "g_legendre_poly_norm_schmidt_deriv: Condition x != 1 failed"
        << endl << "  x = " << x << endl;
      throw runtime_error (os.str ());
    }

  Numeric result;

  if (l == 1)
    {
      if (m == 0)
        {
          result = 1;
        }
      else if (m == 1)
        {
          result = x / sqrt(1 - x * x);
        }
      else
        {
          ostringstream os;
          os << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition l == 1 && (m == 0 || m == 1) failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }

    }
  else if ( m < l )
    {
      try
        {
          result = sqrt(2.0 * fac (l - m) / fac (l + m)) * 
            ((Numeric)l * g_legendre_poly (l - 1, m, x) + 
             (Numeric)(m - l) * x * g_legendre_poly (l, m, x)) /
            (1 - x * x);
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition m < l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }
  else
    {
      try
        {
          result = - sqrt(2.0 * fac (l - m) / fac (l + m)) * 
            (Numeric)m * x * g_legendre_poly (l, m, x) / (1 - x * x);
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition m = l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }

  return (result);
}

//! g_legendre_poly_norm_schmidt_deriv4
/*!
  Returns the derivative of the Schmidt quasi-normalized associated Legendre polynomial Plm(x)) without 
  the factor (-1)^m.

  The input parameters must fulfill the following conditions:
  0 <= m <= l and |x| < 1

  \return      dPlm
  \param   l   Index
  \param   m   Index
  \param   x   Value

  \author Nikolay Koulev
  \date   2003-09-02
*/
Numeric
g_legendre_poly_norm_schmidt_deriv4 (Index l, Index m, Numeric x)
{
  assert (x != 1.);
  if (x == 1.)
    {
      ostringstream os;
      os << "g_legendre_poly_norm_schmidt_deriv: Condition x != 1 failed"
        << endl << "  x = " << x << endl;
      throw runtime_error (os.str ());
    }

  Numeric result;

  if (l == 1)
    {
      if (m == 0)
        {
          result = 1;
        }
      else if (m == 1)
        {
          result = x / sqrt(1 - x * x);
        }
      else
        {
          ostringstream os;
          os << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition l == 1 && (m == 0 || m == 1) failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }

    }
  else if ( m < l )
    {
      try
        {
          result = sqrt (2.0 * fac (l - m) / fac (l + m)) * 
            ((Numeric)((l + m) * (l + 1)) * g_legendre_poly (l - 1, m, x)  - 
             (Numeric)((l + 2 * m) * (l - m + 1)) * g_legendre_poly (l + 1, m, x)
             / ((Numeric)(2 * l + 1 ) * (1 - x * x)));
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition m < l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }
  else
    {
      try
        {
          result = - sqrt (2.0 * fac (l - m) / fac (l + m)) * 
            (Numeric)m * x * g_legendre_poly (l, m, x) / (1 - x * x);
        }
      catch (runtime_error e)
        {
          ostringstream os;
          os << e.what () << "g_legendre_poly_norm_schmidt_deriv: "
            << "Condition m = l failed" << endl
            << "l = " << l << "  m = " << m << endl;
          throw runtime_error (os.str ());
        }
    }

  return (result);
}
