#include "poly_roots.h"

int
main (void)
{
  Vector v(9, 0);

  v[0] = 1.5;
  v[4] = 1;
  v[8] = 1;

  cout << poly_root_solve (v) << endl;

  return (0);
}

