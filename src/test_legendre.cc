#include <iostream>
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
      cout << "Pml = " << legendre_polynomial (l, m, x) << endl;
    }
  catch (runtime_error e)
    {
      cerr << e.what ();
    }

  return (0);
}

