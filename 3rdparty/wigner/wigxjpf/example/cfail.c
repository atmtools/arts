#include <stdio.h>
#include <stdlib.h>

#include "wigxjpf.h"

/* This program is used to test (some) of the wigxjpf error handling.
 *
 * DO NOT USE as a template!
 */

int main(int argc, char *argv[])
{
  double val3j, val6j, val9j;
  int failure = 0;

  if (argc > 1)
    failure = atoi(argv[1]);

  printf ("WIGXJPF C error handling test (%d)\n", failure);
  fflush(stdout); /* Make sure error messages come in order... */

  if (failure != 1)
    wig_table_init((failure != 2) ? 2*100 : 2*1,
		   (failure != 3) ? 9 : -1);
  if (failure != 4)
    wig_temp_init((failure != 5) ? 2*100 : 1);

  val3j = wig3jj(2* 10 , 2* 15 , 2* 10,
		 2*(-3), 2* 12 , 2*(-9));

  if (failure == 6)
    val3j = wig3jj(2* 1000 , 2* 1500 , 2* 1000,
		   2*(-300), 2* 1200 , 2*(-900));

  printf ("3J(10  15  10; -3  12  -9):      %#25.15f\n", val3j);
  fflush(stdout);
  
  val6j = wig6jj(2* 10 , 2* 15 , 2* 10 ,
		 2*  7,  2*  7 , 2*  9 );

  if (failure == 7)
    val6j = wig6jj(2* 1000 , 2* 1500 , 2* 1000 ,
		   2*  700,  2*  700 , 2*  900 );

  printf ("6J{10  15  10;  7   7   9}:      %#25.15f\n", val6j);
  fflush(stdout);
  
  val9j = wig9jj(1,  2,  3,
		 4,  6,  8,
		 3,  6,  9);

  if (failure == 8)
    val9j = wig9jj(100,  200,  300,
		   400,  600,  800,
		   300,  600,  900);

  printf ("9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}: %#25.15f\n", val9j);
  fflush(stdout);
  
  wig_temp_free();
  wig_table_free();

  return 0;
}
