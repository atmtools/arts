#include <cstdlib>
#include <ctime>
#include <iostream>

#if HAVE_UNISTD_H
#include <sys/types.h>
#include <unistd.h>
#endif

#include "matpack_data.h"
#include "sorting.h"

void testVector() {
  // Array for output of sorted indexes
  ArrayOfIndex i;

  Vector v(10);
  v[0] = 2.2;
  v[1] = 1.1;
  v[2] = 3.3;
  v[3] = 7.7;
  v[4] = 6.6;
  v[5] = 9.9;
  v[6] = 4.4;
  v[7] = 8.8;
  v[8] = 5.5;
  v[9] = 10.01;

  std::cout << "Vector before sort:     " << v << std::endl;
  get_sorted_indexes(i, v);
  std::cout << "Index array after sort: " << i << std::endl;
  std::cout << "Sorted Vector:         ";
  for (Index j = 0; j < v.nelem(); j++) std::cout << " " << std::setw(3) << v[i[j]];
  std::cout << std::endl << std::endl;
}

#ifdef _POSIX_VERSION
void testArray() {
  // Array for output of sorted indexes
  ArrayOfIndex i;

  ArrayOfIndex a(10);
  a[0] = 2;
  a[1] = 1;
  a[2] = 3;
  a[3] = 7;
  a[4] = 6;
  a[5] = 9;
  a[6] = 4;
  a[7] = 8;
  a[8] = 5;
  a[9] = 10;

  std::cout << "Array before sort:      " << a << std::endl;
  get_sorted_indexes(i, a);
  std::cout << "Index array after sort: " << i << std::endl;
  std::cout << "Sorted Array:          ";
  for (Index j = 0; j < a.nelem(); j++) std::cout << " " << std::setw(3) << a[i[j]];
  std::cout << std::endl << std::endl;
}

void profileVector(Index n) {
  std::cout << "Creating Vector with random numbers" << std::endl;

  srandom((unsigned int)time(NULL));
  Vector v(n);
  for (Index i = 0; i < n; i++) v[i] = Numeric(random());

  std::cout << "Now sorting" << std::endl;
  ArrayOfIndex i;
  get_sorted_indexes(i, v);
}
#endif

int main(void) {
#ifdef _POSIX_VERSION
  testVector();
  testArray();
#else
  std::cerr << "This test is only available when compiled with POSIX support."
       << std::endl;
#endif

  //  profileVector (100 * 100 * 20 * 20);

  return (0);
}
