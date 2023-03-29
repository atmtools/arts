#include <stdio.h>

/* Header must be included to have the functions declared as extern "C". */
#include "wigxjpf.h"

/* Wrap the calls in a class, to make explicit that we are compiling C++. */
class test
{
public:
  void testcalls();
};

void test::testcalls()
{
  double val3j, val6j, val9j;

  printf ("WIGXJPF C++ test program\n");

  wig_table_init(2*100,9);
  wig_temp_init(2*100);

  /* Note that the arguments to wig3jj, wig6jj and wig9jj are 2*j
   * and 2*m.  To be able to handle half-integer arguments.
   */

  val3j = wig3jj(2* 10 , 2* 15 , 2* 10 ,
		 2*(-3), 2* 12 , 2*(-9));

  printf ("3J(10  15  10; -3  12  -9):      %#25.15f\n", val3j);

  val6j = wig6jj(2* 10 , 2* 15 , 2* 10 ,
		 2*  7,  2*  7 , 2*  9 );

  printf ("6J{10  15  10;  7   7   9}:      %#25.15f\n", val6j);

  val9j = wig9jj(1,  2,  3,
		 4,  6,  8,
		 3,  6,  9);

  printf ("9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}: %#25.15f\n", val9j);

  wig_temp_free();
  wig_table_free();
}

int main()
{
  test t;
  t.testcalls();
  return 0;
}
