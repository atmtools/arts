#include "poly_roots.h"

int
main (void)
{
  Vector v(3, 0);

  v[0] = 1;
  v[2] = 1;

  cout << poly_root_solve (v) << endl;

  return (0);
}

