#include <iostream>
#include <sys/times.h>
#include <unistd.h>
#include "legendre.h"
#include "exceptions.h"

int
main (int argc, char *argv[])
{
  Index l, m;
  Numeric x;

  if (argc != 4)
    {
      cerr << "Usage: " << argv[0] << " l m x" << endl;
      exit (1);
    }

  l = atoi (argv[1]);
  m = atoi (argv[2]);
  x = strtod (argv[3], NULL);

  try
    {
      cout << "l = " << l << "  m = " << m << "  x = " << x << endl;
      cout << "Pml = " << legendre_poly (l, m, x) << endl;
      cout << "Norm Pml = " << legendre_poly_norm (l, m, x) << endl;
    }
  catch (runtime_error e)
    {
      cerr << e.what ();
    }

  struct tms cput_start, cput_end;
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

  return (0);
}

