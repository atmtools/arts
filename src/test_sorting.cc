#include <iostream>
using namespace std;

#include "matpackI.h"
#include "sorting.h"

int
main (void)
{
  // Array for output of sorted indexes
  ArrayOfIndex i;

  /*
   * First sort a vector
   */
  Vector v (10);
  v [0] = 2.2;
  v [1] = 1.1;
  v [2] = 3.3;
  v [3] = 7.7;
  v [4] = 6.6;
  v [5] = 9.9;
  v [6] = 4.4;
  v [7] = 8.8;
  v [8] = 5.5;
  v [9] = 10.01;

  cout << "Vector before sort:     " << v << endl;
  get_sorted_indexes (i, v);
  cout << "Index array after sort: " << i << endl;
  cout << "Sorted Vector:         ";
  for (Index j = 0; j < v.nelem (); j++)
    cout << " " << setw (3) << v[i[j]];
  cout << endl << endl;

  /*
   * Now sorting an ArrayOfIndex
   */
  ArrayOfIndex a (10);
  a [0] = 2;
  a [1] = 1;
  a [2] = 3;
  a [3] = 7;
  a [4] = 6;
  a [5] = 9;
  a [6] = 4;
  a [7] = 8;
  a [8] = 5;
  a [9] = 10;

  cout << "Array before sort:      " << a << endl;
  get_sorted_indexes (i, a);
  cout << "Index array after sort: " << i << endl;
  cout << "Sorted Array:          ";
  for (Index j = 0; j < a.nelem (); j++)
    cout << " " << setw (3) << a[i[j]];
  cout << endl;

  return (0);
}

