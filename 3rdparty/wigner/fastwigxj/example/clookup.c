#include <stdio.h>
#include <string.h>

#include "fastwigxj.h"

#include "wigxjpf.h"

int main(int argc, char *argv[])
{
  double val3j, val6j, val9j;

  printf ("FASTWIGXJ C test program\n");

  /* Load tables produced during build test. */
  
  fastwigxj_load("test_table_18.3j",         3, NULL);
  fastwigxj_load("test_table_8.6j",          6, NULL);
  fastwigxj_load("test_hashed_3_6_12.9j",    9, NULL);

  /* For fallback of 9j to sum-of-6j-products when
   * symbols too large for table.
   *
   * And a dirty hack to not test this when quadmath/float128
   * is not available.
   */

  if (!(argc > 1 && strcmp (argv[1],"--no-float128") == 0))
    fastwigxj_load("test_table_8_float128.6j", 7, NULL);

  /* For fallback to WIGXJPF when symbols too large for table. */

  wig_table_init(2*100,9);
  wig_temp_init(2*100);

  /* Note that the arguments to wig3jj, wig6jj and wig9jj are 2*j
   * and 2*m.  To be able to handle half-integer arguments.
   *
   * Also note that for 3j symbols, the third m is not given,
   * as it is fully redundant (unless trivial-0 symbols are wanted).
   *
   * (A version which takes the third m exists as fw3jja6).
   */

  val3j = fw3jja(2*  5 , 2*  7 , 2*  5 ,
		 2*(-3), 2*  5);

  printf ("3J( 5   7   5; -3   5  -2):      %#25.15f\n", val3j);

  val3j = fw3jja(2* 10 , 2* 15 , 2* 10 ,
		 2*(-3), 2* 12);

  printf ("3J(10  15  10; -3  12  -9):      %#25.15f\n", val3j);

  val3j = fw3jja6(2*  5 , 2*  7 , 2*  5 ,
		  2*(-3), 2*  5 , 2*(-2));

  printf ("3J( 5   7   5; -3   5  -2):      %#25.15f\n", val3j);

  val3j = fw3jja6(2* 10 , 2* 15 , 2* 10 ,
		  2*(-3), 2* 12 , 2*(-9));

  printf ("3J(10  15  10; -3  12  -9):      %#25.15f\n", val3j);

  val6j = fw6jja(2*  3 , 2*  4 , 2*  2 ,
		 2*  2 , 2*  2 , 2*  3 );

  printf ("6J{ 3   4   2;  2   2   3}:      %#25.15f\n", val6j);

  val6j = fw6jja(2* 10 , 2* 15 , 2* 10 ,
		 2*  7 , 2*  7 , 2*  9 );

  printf ("6J{10  15  10;  7   7   9}:      %#25.15f\n", val6j);

  val9j = fw9jja(1,  2,  3,
		 2,  3,  5,
		 3,  3,  6);

  printf ("9J{0.5 1 1.5; 1 1.5 2.5; 1.5 1.5 3}:%#22.15f\n", val9j);

  /* The following symbol is not in the 9j table, and will
   * automatically go via 6j fallback (through float128).
   */
  
  val9j = fw9jja(1,  2,  3,
		 4,  6,  8,
		 3,  6,  9);

  printf ("9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}: %#25.15f\n", val9j);

  {
    /* Call the lookup functions with argument lists instead.
     * The lists are not destroyed, so can be partially reused.
     */
    int arg3j[5], arg6j[6], arg9j[9];

    arg3j[0] = 2*  5;  arg3j[1] = 2*  7;  arg3j[2] = 2*  5;
    arg3j[3] = 2*(-3); arg3j[4] = 2*  5;

    val3j = fw3jjl(arg3j);

    printf ("3J( 5   7   5; -3   5  -2):      %#25.15f\n", val3j);

    arg3j[0] = 2* 10;  arg3j[1] = 2* 15;  arg3j[2] = 2* 10;
    arg3j[3] = 2*(-3); arg3j[4] = 2* 12;

    val3j = fw3jjl(arg3j);

    printf ("3J(10  15  10; -3  12  -9):      %#25.15f\n", val3j);

    arg3j[0] = 2* 11;      /* Partial reuse. */

    val3j = fw3jjl(arg3j);

    printf ("3J(11  15  10; -3  12  -9):      %#25.15f\n", val3j);

    arg6j[0] = 2*  3;  arg6j[1] = 2*  4;  arg6j[2] = 2*  2;
    arg6j[3] = 2*  2;  arg6j[4] = 2*  2;  arg6j[5] = 2*  3;

    val6j = fw6jjl(arg6j);

    printf ("6J{ 3   4   2;  2   2   3}:      %#25.15f\n", val6j);

    arg6j[0] = 2* 10;  arg6j[1] = 2* 15;  arg6j[2] = 2* 10;
    arg6j[3] = 2*  7;  arg6j[4] = 2*  7;  arg6j[5] = 2*  9;

    val6j = fw6jjl(arg6j);

    printf ("6J{10  15  10;  7   7   9}:      %#25.15f\n", val6j);

    arg9j[0] = 1;  arg9j[1] = 2;  arg9j[2] = 3;
    arg9j[3] = 4;  arg9j[4] = 6;  arg9j[5] = 8;
    arg9j[6] = 3;  arg9j[7] = 6;  arg9j[8] = 9;

    val9j = fw9jjl(arg9j);

    printf ("9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}: %#25.15f\n", val9j);
  }

  {
    /* By using separate canonicalisation (and prefetch) and value
     * retrieval, several lookups can be interleaved.
     */

    int arg3j_1[5];
    int arg3j_2[5];
    int arg3j_3[5];
    uint64_t x_1;
    uint64_t x_2;
    uint64_t x_3;

    arg3j_1[0] = 2*  5;  arg3j_1[1] = 2*  7;  arg3j_1[2] = 2*  5;
    arg3j_1[3] = 2*(-3); arg3j_1[4] = 2*  5;

    fw3jj_canon(arg3j_1, &x_1);
    fw3jj_prefetch(x_1);

    arg3j_2[0] = 2* 10;  arg3j_2[1] = 2* 15;  arg3j_2[2] = 2* 10;
    arg3j_2[3] = 2*(-3); arg3j_2[4] = 2* 12;

    fw3jj_canon(arg3j_2, &x_2);
    fw3jj_prefetch(x_2);

    arg3j_3[0] = 2* 11;  arg3j_3[1] = 2* 15;  arg3j_3[2] = 2* 10;
    arg3j_3[3] = 2*(-3); arg3j_3[4] = 2* 12;

    fw3jj_canon(arg3j_3, &x_3);
    fw3jj_prefetch(x_3);

    /* The list of original arguments must be provided at lookup
     * also, as they are used in case the value is not in the table
     * and a fallback calculation is required.
     */
    
    val3j = fw3jj_get(arg3j_1, x_1);

    printf ("3J( 5   7   5; -3   5  -2):      %#25.15f\n", val3j);

    val3j = fw3jj_get(arg3j_2, x_2);

    printf ("3J(10  15  10; -3  12  -9):      %#25.15f\n", val3j);

    val3j = fw3jj_get(arg3j_3, x_3);

    printf ("3J(11  15  10; -3  12  -9):      %#25.15f\n", val3j);
  } 

  /* Print statistics. */

  printf ("\n");
  fastwigxj_print_stats();

  /* Remove tables from memory. */

  fastwigxj_unload(3);
  fastwigxj_unload(6);
  fastwigxj_unload(9);

  /* Release WIGXJPF memory. */
  
  wig_temp_free();
  wig_table_free();

  return 0;
}
