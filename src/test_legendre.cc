/* Copyright (C) 2003-2012 Oliver Lemke <olemke@core-dump.info>

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

#include "arts.h"

#if HAVE_UNISTD_H
#include <sys/types.h>
#include <unistd.h>
#endif
#ifdef _POSIX_VERSION
#include <sys/times.h>
#endif
#include "exceptions.h"
#include "legendre.h"

void test_gsl_int() {
  Vector x, w;

  std::cout << gsl_integration_glfixed_table_alloc(x, w, 20) << std::endl;
  std::cout << x << std::endl;
  std::cout << w << std::endl << std::endl;
}

int main(int argc, char *argv[]) {
  Index l, m;
  Numeric x;

  if (argc != 4) {
    cerr << "Usage: " << argv[0] << " l m x" << endl;
    exit(1);
  }

  l = atoi(argv[1]);
  m = atoi(argv[2]);
  x = strtod(argv[3], NULL);

  try {
    cout << "l = " << l << "  m = " << m << "  x = " << x << endl;
    cout << "Pml = " << legendre_poly(l, m, x) << endl;
    cout << "dPml = " << legendre_poly_deriv(l, m, x) << endl;
    cout << "Norm Pml = " << legendre_poly_norm_schmidt(l, m, x) << endl;
    cout << "Norm dPml = " << legendre_poly_norm_schmidt_deriv(l, m, x) << endl;
  } catch (const std::runtime_error &e) {
    cerr << e.what();
  }

  /*  struct tms cput_start, cput_end;
  const Index n = 1000000;
  Index clktck;
  Vector v2 (n), r2 (n);

  cout << endl << "Timing ("
    << n << " loops):" << endl;

  for (Index i = 0; i < n; i++)
    {
      v2[i] = -1.0 + Numeric (i) / n * 2.0;
    }

  if ((clktck = sysconf (_SC_CLK_TCK)) < 0)
    throw runtime_error ("Timer error: Unable to determine CPU clock ticks");

  if (times (&cput_start) == -1)
    throw runtime_error ("Timer error: Unable to get current CPU time");

  for (Index i = 0; i < n; i++)
    r2[i] = legendre_poly (l, m, v2[i]);

  if (times (&cput_end) == -1)
    throw runtime_error ("Timer error: Unable to get current CPU time");

  cout << "CPU time: " << setprecision (2)
    << ((cput_end.tms_stime - cput_start.tms_stime)
        + (cput_end.tms_utime - cput_start.tms_utime))
    / (Numeric)clktck << " s" << endl;
*/

  return (0);
}
