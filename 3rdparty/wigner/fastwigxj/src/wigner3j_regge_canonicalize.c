
/* Copyright 2015 Haakan T. Johansson */

/*  This file is part of FASTWIGXJ.
 *
 *  FASTWIGXJ is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  FASTWIGXJ is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with FASTWIGXJ.  If not, see
 *  <http://www.gnu.org/licenses/>.
 */

#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <stdio.h>
#include "canonicalise_inc.h"

#include "wigner36j_regge_canonicalize.h"

/* From WIGXJPF: trivial_zero.c
 */

#define COLLECT_NEGATIVE(two_j1,two_j2,two_j3) do {   \
  collect_sign |= (two_j1) | (two_j2) | (two_j3);     \
  } while (0)

#define COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j1,two_j2,two_j3) do {      \
  collect_odd |= (two_j1) + (two_j2) + (two_j3);                      \
  collect_sign |= (((two_j2) + (two_j3)) - (two_j1));                 \
  collect_sign |= ((two_j1) - ((two_j2) - (two_j3)));                 \
  collect_sign |= ((two_j1) - ((two_j3) - (two_j2)));                 \
  } while (0)

#define COLLECT_ABS_M_WITHIN_J(two_m,two_j) do {      \
  collect_odd |= ((two_m) + (two_j));                 \
  collect_sign |= ((two_j) - (two_m));                \
  collect_sign |= ((two_j) + (two_m));                \
  } while (0)

/* End from WIGXJPF. */

#define DEBUG_WIG3J_R_C  0

uint64_t wigner3j_regge_canonicalise(const int *two_jv)
{
  /* Canonicalisation (index) according to:
   *
   * EFFICIENT STORAGE SCHEME FOR PRECALCULATED WIGNER
   * 3J, 6J AND GAUNT COEFFICIENTS
   *
   * by J. RASCH AND A. C. H. YU
   *
   * SIAM J. SCI. COMPUT.  Vol. 25, No. 4, pp. 1416--1428
   */

  /* Note: we will not use the index, but a hash table instead. */

  /* Organisation:
   *
   * Calculate:
   *
   * -j1 + j2 + j3   j1 - j2 + j3   j1 + j2 - j3
   *  j1 - m1        j2 - m2        j3 - m3
   *  j1 + m1        j2 + m2        j3 + m3
   *
   * Which at the end shall be transformed (by row & column
   * interchange, and transposition) into:
   *
   * S       L       B+X-T
   * X       B       S+L-T
   * B+L-T   S+X-T   T
   *
   * S <= B <= T <= X <= L
   *
   * That is done in three steps:
   *
   * - Find (one of) the smallest element(s), place it as S.
   *   (need to consider all 9).
   *   (row/col interchange, no transposition)
   *
   * - Find the largest element, place it as L.
   *   (need to consider 4 elements, same row/col as S).
   *   (transposition, col interchange).
   *
   * - Make R22 < R32, or if R22 = R32, then R23 < R33.
   *   (row interchange)
   */

  int ja =  two_jv[0];
  int jb =  two_jv[1];
  int jc =  two_jv[2];
  int ma =  two_jv[3];
  int mb =  two_jv[4];
  int mc =  - ma - mb;
  
#if DEBUG_WIG3J_R_C
  fprintf (stderr, "%3d %3d %3d %3d %3d %3d jabcmabc\n",
	   ja, jb, jc, ma, mb, -ma - mb);
#endif

  /* Note that mc is always ma + mb (for non-zero 3j-symbol), so not
   * given as an argument to this function.
   *
   * Also note that we never (for non-zero symbols) have any negative
   * numbers in the matrix.
   */

  uint64_t j1 = (uint64_t) ja;
  uint64_t j2 = (uint64_t) jb;
  uint64_t j3 = (uint64_t) jc;
  uint64_t m1 = (uint64_t) ma;
  uint64_t m2 = (uint64_t) mb;
  uint64_t m3 = (uint64_t) mc;

  {
#define two_j1 ja
#define two_j2 jb
#define two_j3 jc
#define two_m1 ma
#define two_m2 mb
#define two_m3 mc

    /* From WIGXJPF: trivial_zero.c
     */

    int collect_sign = 0;
    int collect_odd = 0;

    COLLECT_NEGATIVE(two_j1, two_j2, two_j3);
    COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j1, two_j2, two_j3);

    COLLECT_ABS_M_WITHIN_J(two_m1, two_j1);
    COLLECT_ABS_M_WITHIN_J(two_m2, two_j2);
    COLLECT_ABS_M_WITHIN_J(two_m3, two_j3);

    if ((two_m1 + two_m2 + two_m3) |
	(collect_sign & (1 << (sizeof (int) * 8 - 1))) | 
	(collect_odd & 1))
      {
	STATS_3J->_trivial0++;
	return 0;
      }
  }

#undef two_j1
#undef two_j2
#undef two_j3
#undef two_m1
#undef two_m2
#undef two_m3
#define two_j1 j1
#define two_j2 j2
#define two_j3 j3
#define two_m1 m1
#define two_m2 m2
#define two_m3 m3

  uint64_t two_js = j1 + j2 + j3;

  uint64_t two_js_js_js =
    (two_js << (32+4)) | (two_js << (16+4)) | (two_js << (4));
  /* By adding, instead of or'ing, the m values, any bits marking
   * negative values are superpositioned.  They will then cancel
   * properly when subtracted or added to the j values.
   * Adding or or'ing the j values does not matter.
   */
  uint64_t two_j1_j2_j3 =
    (two_j1 << (32+4)) + (two_j2 << (16+4)) + (two_j3 << (4));
  uint64_t two_m1_m2_m3 =
    (two_m1 << (32+4)) + (two_m2 << (16+4)) + (two_m3 << (4));

  uint64_t col_ind = (0ll << (32+2)) | (1ll << (16+2)) | (2ll << (2));
  uint64_t row_ind = (1ll << (32+0)) | (1ll << (16+0)) | (1ll << (0));

  uint64_t row1 = two_js_js_js - 2 * two_j1_j2_j3 + col_ind;
  uint64_t row2 = two_j1_j2_j3 - two_m1_m2_m3 + col_ind + row_ind;
  uint64_t row3 = two_j1_j2_j3 + two_m1_m2_m3 + col_ind + (row_ind << 1);

#if DEBUG_WIG3J_R_C
  fprintf (stderr, "%016" PRIx64 "  %016" PRIx64 "  %016" PRIx64 "\n",
	   two_js_js_js, two_j1_j2_j3, two_m1_m2_m3);
  fprintf (stderr, "%016" PRIx64 "  %016" PRIx64 "  %016" PRIx64 " rows\n",
	   row1, row2, row3);
  fprintf (stderr, "%3d %3d %3d  %3d %3d %3d  %3d %3d %3d\n",
	   (row1 >> 36) & 0xfff, (row1 >> 20) & 0xfff, (row1 >> 4) & 0xfff,
	   (row2 >> 36) & 0xfff, (row2 >> 20) & 0xfff, (row2 >> 4) & 0xfff,
	   (row3 >> 36) & 0xfff, (row3 >> 20) & 0xfff, (row3 >> 4) & 0xfff);
#endif

  /* There would now be two ways of finding and selecting the smallest
   * elements.  The issue is that we must not only find the element,
   * but do the same transformations to the entire matrix when placing
   * the smallest element at the upper left corner.
   *
   * - Find the element, keep its index, and then do one transformation.
   *
   * - Make the transformations to move the smallest and keep as such.
   *
   * Transforming is rather expensive, so rather find the element.
   * That has been prepared by putting the column and row indices just
   * after the values.  They will thus not affect the comparisons.
   */

  uint64_t smallest  = row1 >> 32;
  /* Working with two smallest variables and combing at the end
   * did not improve speed.
   */

  if (((row1 >> 16) & 0xffff) < smallest) smallest = (row1 >> 16) & 0xffff;
  if (((row1      ) & 0xffff) < smallest) smallest = (row1      ) & 0xffff;

  if (((row2 >> 32)         ) < smallest) smallest = (row2 >> 32);
  if (((row2 >> 16) & 0xffff) < smallest) smallest = (row2 >> 16) & 0xffff;
  if (((row2      ) & 0xffff) < smallest) smallest = (row2      ) & 0xffff;

  if (((row3 >> 32)         ) < smallest) smallest = (row3 >> 32);
  if (((row3 >> 16) & 0xffff) < smallest) smallest = (row3 >> 16) & 0xffff;
  if (((row3      ) & 0xffff) < smallest) smallest = (row3      ) & 0xffff;

#if DEBUG_WIG3J_R_C
  fprintf (stderr, "%016" PRIx64 " smallest\n", smallest);
#endif

  /* The column index of the smallest is now in bits 3..2 of smallest,
   * and the row index is in bits 1..0. */

  uint64_t small_col = (smallest & (3 << 2)) << 2; /* 0, 4, 8 -> 0, 16, 32 */

  /* Column permutations.  We also kill the index information.
   * We either change no columns, or interchange twice, so no sign change.
   */

  uint64_t small_col_rsl = 48 - small_col;

  row1 = ((row1 << small_col) | (row1 >> small_col_rsl)) & 0xfff0fff0fff0ll;
  row2 = ((row2 << small_col) | (row2 >> small_col_rsl)) & 0xfff0fff0fff0ll;
  row3 = ((row3 << small_col) | (row3 >> small_col_rsl)) & 0xfff0fff0fff0ll;

  /* The columns are now in order.  Exchange rows as needed.
   * If exchange is made, we have also changed the sign.
   */

  uint64_t tmprow1 = row1;
#if 0
  if (smallest & 1) { row1 = row2; row2 = tmprow1; }
  /* tmprow1 = row1; */ /* we never do both interchanges */
  if (smallest & 2) { row1 = row3; row3 = tmprow1; }
#else
  uint64_t tmprow2 = row2;
  uint64_t tmprow3 = row3;
#if 1
  /* This masking game saves ~2 ns compared to the if statements above. */
  /* all 0 when bit 1 set, all 1 when bit 1 not set */
  uint64_t noswapmask12 = (smallest & 1) - 1;
  uint64_t   swapmask12 = ~noswapmask12;
  /* all 0 when bit 1 set, all 1 when bit 1 not set */
  uint64_t noswapmask13 = ((smallest >> 1) & 1) - 1;
  uint64_t   swapmask13 = ~noswapmask13;
  row1 = (tmprow1 & noswapmask12 & noswapmask13) |
    (tmprow2 & swapmask12) | (tmprow3 & swapmask13);
  row2 = (tmprow2 & noswapmask12) | (tmprow1 & swapmask12);
  row3 = (tmprow3 & noswapmask13) | (tmprow1 & swapmask13);
#else
  /* This did not help since the compiler produced conditional jumps. */
  row1 = (smallest & 1) ? row2 : tmprow1;
  row2 = (smallest & 1) ? tmprow1 : row2;
  tmprow1 = row1;
  row1 = (smallest & 2) ? row3 : tmprow1;
  row3 = (smallest & 2) ? tmprow1 : row3;
#endif
#endif  

  uint64_t sign = smallest | (smallest >> 1); /* bit 0 matters */
  
#if DEBUG_WIG3J_R_C
  fprintf (stderr, "%016" PRIx64 "  %016" PRIx64 "  %016" PRIx64 " rows\n",
	   row1, row2, row3);
  fprintf (stderr, "%3d %3d %3d  %3d %3d %3d  %3d %3d %3d\n",
	   (row1 >> 36) & 0xfff, (row1 >> 20) & 0xfff, (row1 >> 4) & 0xfff,
	   (row2 >> 36) & 0xfff, (row2 >> 20) & 0xfff, (row2 >> 4) & 0xfff,
	   (row3 >> 36) & 0xfff, (row3 >> 20) & 0xfff, (row3 >> 4) & 0xfff);
#endif

  /* We now need to find the largest element.  It will be either in
   * the first row or the first column.  We need to remember both if
   * we have to perform a transposition, and a column interchange.
   */

  uint64_t largest =  (row2 >> 32) | (2    ); /* transpose */
  uint64_t lrow3   =  (row3 >> 32) | (2 | 1); /* transpose, col23 exchange */

  if (lrow3 > largest) { largest = lrow3; }
  uint64_t lcol2   =  (row1 >> 16) & 0xfff0;  /* none */
  uint64_t lcol3   = ((row1      ) & 0xfff0) | (    1); /* col23 exchange */

  if (lcol2 > largest) { largest = lcol2; }
  if (lcol3 > largest) { largest = lcol3; }

#if DEBUG_WIG3J_R_C
  fprintf (stderr, "%016" PRIx64 "\n", largest);
#endif


  /* Calculate transposition.  Use it if wanted. */
  /* Can be beautifully done with SSE permutation functions. */
  
  uint64_t trow1 = (((row1      ) & 0xfff000000000ll) |
		    ((row2 >> 16) & 0xfff00000ll) |
		    ((row3 >> 32)));
  uint64_t trow2 = (((row1 << 16) & 0xfff000000000ll) |
		    ((row2      ) & 0xfff00000ll) |
		    ((row3 >> 16) & 0xfff0ll));
  uint64_t trow3 = (((row1 << 32)) |
		    ((row2 << 16) & 0xfff00000ll) |
		    ((row3      ) & 0xfff0ll));

#if DEBUG_WIG3J_R_C
  fprintf (stderr, "%016" PRIx64 "  %016" PRIx64 "  %016" PRIx64 " trows\n",
	   trow1, trow2, trow3);
  fprintf (stderr, "%3d %3d %3d  %3d %3d %3d  %3d %3d %3d\n",
	   (row1 >> 36) & 0xfff, (row1 >> 20) & 0xfff, (row1 >> 4) & 0xfff,
	   (row2 >> 36) & 0xfff, (row2 >> 20) & 0xfff, (row2 >> 4) & 0xfff,
	   (row3 >> 36) & 0xfff, (row3 >> 20) & 0xfff, (row3 >> 4) & 0xfff);
#endif

  if (largest & 2) { row1 = trow1; row2 = trow2; row3 = trow3; }

#if DEBUG_WIG3J_R_C
  fprintf (stderr, "%016" PRIx64 "  %016" PRIx64 "  %016" PRIx64 " rows\n",
	   row1, row2, row3);
  fprintf (stderr, "%3d %3d %3d  %3d %3d %3d  %3d %3d %3d\n",
	   (row1 >> 36) & 0xfff, (row1 >> 20) & 0xfff, (row1 >> 4) & 0xfff,
	   (row2 >> 36) & 0xfff, (row2 >> 20) & 0xfff, (row2 >> 4) & 0xfff,
	   (row3 >> 36) & 0xfff, (row3 >> 20) & 0xfff, (row3 >> 4) & 0xfff);
#endif

  /* Then do the column exchanges. */

  uint64_t colshift = (largest & 1) << 4; /* 0 or 16 */

  row1 = (row1 & 0xfff000000000ll) | 
    ((row1 & 0xfff00000ll) >> colshift) | ((row1 & 0xfff0ll) << colshift);
  row2 = (row2 & 0xfff000000000ll) | 
    ((row2 & 0xfff00000ll) >> colshift) | ((row2 & 0xfff0ll) << colshift);
  row3 = (row3 & 0xfff000000000ll) | 
    ((row3 & 0xfff00000ll) >> colshift) | ((row3 & 0xfff0ll) << colshift);

#if DEBUG_WIG3J_R_C
  fprintf (stderr, "%016" PRIx64 "  %016" PRIx64 "  %016" PRIx64 " rows\n",
	   row1, row2, row3);
  fprintf (stderr, "%3d %3d %3d  %3d %3d %3d  %3d %3d %3d\n",
	   (row1 >> 36) & 0xfff, (row1 >> 20) & 0xfff, (row1 >> 4) & 0xfff,
	   (row2 >> 36) & 0xfff, (row2 >> 20) & 0xfff, (row2 >> 4) & 0xfff,
	   (row3 >> 36) & 0xfff, (row3 >> 20) & 0xfff, (row3 >> 4) & 0xfff);
#endif

  sign ^= largest; /* bit 0 matters */

  /* We have now placed S and L at their locations. */
  
  if ((row3 & 0xfff0fff0ll) < (row2 & 0xfff0fff0ll))
    {
      uint64_t tmp = row2; row2 = row3; row3 = tmp;
      sign ^= 1; /* bit 0 matters */
    }

#if DEBUG_WIG3J_R_C
  fprintf (stderr, "%016" PRIx64 "  %016" PRIx64 "  %016" PRIx64 " rows\n",
	   row1, row2, row3);
  fprintf (stderr, "%3d %3d %3d  %3d %3d %3d  %3d %3d %3d\n",
	   (row1 >> 36) & 0xfff, (row1 >> 20) & 0xfff, (row1 >> 4) & 0xfff,
	   (row2 >> 36) & 0xfff, (row2 >> 20) & 0xfff, (row2 >> 4) & 0xfff,
	   (row3 >> 36) & 0xfff, (row3 >> 20) & 0xfff, (row3 >> 4) & 0xfff);
#endif

  /* And the array is fully organised. */
  
  /* Pick out the items and place them in a 64-bit array.  Five items
   * in 64 bits give 12 bits per item.
   *
   * 36 52 S
   * 20 40 L
   * 36 28 X
   * 20 16 B
   *  4  4 T
   *  0  0 sign
   */

  uint64_t key =
    ((row1 & 0xfff000000000ll) << 16) | // S
    ((row1 & 0x0000fff00000ll) << 20) | // L
    ((row2 & 0xfff000000000ll) >>  8) | // X
    ((row2 & 0x0000fff00000ll) >>  4) | // B
    ((row3 & 0x00000000fff0ll)      ) | // T
    (sign & 1);

  /* + 4 to shift away the empty row/col locations, +1 to shift away
   * the odd sign (which never is).
   */

  uint64_t S = ((row1 & 0xfff000000000ll) >> (32 + 4 + 1));
  uint64_t L = ((row1 & 0x0000fff00000ll) >> (16 + 4 + 1));
  uint64_t X = ((row2 & 0xfff000000000ll) >> (32 + 4 + 1));
  uint64_t B = ((row2 & 0x0000fff00000ll) >> (16 + 4 + 1));
  uint64_t T = ((row3 & 0x00000000fff0ll) >> ( 0 + 4 + 1));

#if 1
  /* As we have the values 2*j, we multiply by this 2 for all other
   * contributions.  Maximum a power j^6 is encountered in Eind.
   */

  uint64_t Lind = L * (24 + L * (50 + L * (35 + L * (10 + L))));
  uint64_t Xind = X * (6 + X * (11 + X * (6 + X)));
  uint64_t Tind = T * (2 + T * (3 + T));
  uint64_t Bind = B * (1 + B);
  uint64_t Sind = S;

  uint64_t index = (1  * Lind +
		    5  * Xind +
		    20 * Tind +
		    60 * Bind +
		    120 * Sind) / (120) + 1;
  
  sign &= (two_js >> 1) & 1;

#if DEBUG_WIG3J_R_C
  fprintf (stderr,"%d %d %d %d %d , %d %d %d %d %d -> %d s%d\n",
	   L, X, T, B, S,
	   Lind, Xind, Tind, Bind, Sind, index, sign);
#endif

#ifdef KEEP_MAX_RASCH_YU_INDEX
  if (index > _max_RaschYu_index)
    _max_RaschYu_index = index;
#endif
#if 1
  return (index << 1) | sign;
#endif
#endif

  /* As we intend to use a hash table, we do not need the values ordered
   * in an index.
   */

#if DEBUG_WIG3J_R_C
  fprintf (stderr, "-> %016" PRIx64 " (%016" PRIx64 ")\n", key, key & ~1);
#endif

  return key;

  /* And when the need arises to retrieve a set of jn from the c14n
   * variables (e.g. to calculate the value of the symbol):
   *
   * -j1 + j2 + j3 == S
   *  j1 - j2 + j3 == L
   *  j1 - m1      == X
   *  j2 - m2      == B
   *  j3 + m3      == T
   *
   * and  m1 + m2 + m3 == 0
   *
   * gives:
   *
   * j1 = (B + L - T + X)/2
   * j2 = (B + S - T + X)/2
   * j3 = (L + S)/2
   * m1 = j1 - 2*X/2
   * m2 = j2 - 2*B/2
   * m3 = - m1 - m2
   */
}

#undef COLLECT_NEGATIVE
#undef COLLECT_TRIANGLE_TRIVIAL_ZERO
#undef COLLECT_ABS_M_WITHIN_J
