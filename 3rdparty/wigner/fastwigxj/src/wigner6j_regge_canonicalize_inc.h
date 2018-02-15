
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

/* From WIGXJPF: trivial_zero.c
 * (disabled collect_odd)
 */

#define COLLECT_NEGATIVE(two_j1,two_j2,two_j3) do {     \
    collect_sign |= (two_j1) | (two_j2) | (two_j3);     \
  } while (0)

#define COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j1,two_j2,two_j3) do {        \
    /* collect_odd |= (two_j1) + (two_j2) + (two_j3); */		\
    collect_sign |= (((two_j2) + (two_j3)) - (two_j1));                 \
    collect_sign |= ((two_j1) - ((two_j2) - (two_j3)));                 \
    collect_sign |= ((two_j1) - ((two_j3) - (two_j2)));                 \
  } while (0)

/* End from WIGXJPF. */

#if WIGNER6J_REGGE_CANONICALISE_ARRAY4
#define W6U32 v4su32
#define W6U64 v4su64
#else
#define W6U32 uint32_t
#define W6U64 uint64_t
#endif

#if WIGNER6J_REGGE_CANONICALISE_ARRAY4
void
#else
W6U64
#endif
WIGNER6J_REGGE_CANONICALISE_NAME(
#if WIGNER6J_REGGE_CANONICALISE_ARRAY4
				 v4su64 *key,
#endif
#if WIGNER6J_REGGE_CANONICALISE_IN_STRUCT
# if WIGNER6J_REGGE_CANONICALISE_ARRAY4
				 const wigner6j_4_symbol *js
#  define two_j1 js->jv.two_j1
#  define two_j2 js->jv.two_j2
#  define two_j3 js->jv.two_j3
#  define two_j4 js->jv.two_j4
#  define two_j5 js->jv.two_j5
#  define two_j6 js->jv.two_j6
# else
				 const wigner6j_symbol *js
#  define two_j1 js->j.two_j1
#  define two_j2 js->j.two_j2
#  define two_j3 js->j.two_j3
#  define two_j4 js->j.two_j4
#  define two_j5 js->j.two_j5
#  define two_j6 js->j.two_j6
# endif
#else
				 const int *two_jv
#endif
				 )
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

  /* Note that with jx having half-integer values (integer in 2j),
   * both ai and bj are by construction integers (thus even in 2j),
   * for non-trivially 0 input.
   * Thus dividing by two already here to have simpler indexing.
   */

#if !WIGNER6J_REGGE_CANONICALISE_IN_STRUCT
  uint32_t two_j1 =  (uint32_t) two_jv[0];
  uint32_t two_j2 =  (uint32_t) two_jv[1];
  uint32_t two_j3 =  (uint32_t) two_jv[2];
  uint32_t two_j4 =  (uint32_t) two_jv[3];
  uint32_t two_j5 =  (uint32_t) two_jv[4];
  uint32_t two_j6 =  (uint32_t) two_jv[5];
#endif

  W6U32 b1 = (W6U32) (two_j1 + two_j2 + two_j3);
  W6U32 b2 = (W6U32) (two_j1 + two_j5 + two_j6);
  W6U32 b3 = (W6U32) (two_j4 + two_j2 + two_j6);
  W6U32 b4 = (W6U32) (two_j4 + two_j5 + two_j3);

#if WIGNER6J_REGGE_CANONICALISE_TRIVIAL_0_CHECK
  W6U32 collect_sign = 0;
  W6U32 collect_odd = 0;

  collect_odd = b1 | b2 | b3 | b4;

  COLLECT_NEGATIVE(two_j1, two_j2, two_j3);
  COLLECT_NEGATIVE(two_j4, two_j5, two_j6);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j1, two_j2, two_j3);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j1, two_j5, two_j6);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j4, two_j2, two_j6);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j4, two_j5, two_j3);

  if ((collect_sign & (uint32_t) (1 << (sizeof (int) * 8 - 1))) | 
      (collect_odd & 1))
    {
      /* trivial 0 */
      STATS_6J->_trivial0++;
# if WIGNER6J_REGGE_CANONICALISE_RETURN_INDEX
      return 0;
# else
      return (uint64_t) -2;
# endif
    }
#endif

#if WIGNER6J_REGGE_CANONICALISE_ARRAY4
# define SHIFT_DOWN_1(a)  __builtin_ia32_psrldi128(a,1)
#else
# define SHIFT_DOWN_1(a)  ((W6U32) (a) >> 1)
#endif

  b1 = SHIFT_DOWN_1(b1); /* /= 2 */
  b2 = SHIFT_DOWN_1(b2);
  b3 = SHIFT_DOWN_1(b3);
  b4 = SHIFT_DOWN_1(b4);

  W6U32 a1 = SHIFT_DOWN_1(two_j1 + two_j2 + two_j4 + two_j5); /* /= 2 */
  W6U32 a2 = SHIFT_DOWN_1(two_j1 + two_j3 + two_j4 + two_j6);
  W6U32 a3 = SHIFT_DOWN_1(two_j2 + two_j3 + two_j5 + two_j6);

#if WIGNER6J_REGGE_CANONICALISE_ARRAY4
#define SWAP_TO_FIRST_LARGER(tmptype,a,b) do {		\
    W6C_KEEPMAXMIN(a,b);				\
  } while (0)
#else
  /* Looks a bit backwards, but generates code using conditional moves. */
#define SWAP_TO_FIRST_LARGER(tmptype,a,b) do {		\
    tmptype __tmp_a = a;				\
    tmptype __tmp_b = b;				\
    a = (__tmp_a > __tmp_b) ? __tmp_a : __tmp_b;	\
    b = (__tmp_a < __tmp_b) ? __tmp_a : __tmp_b;	\
  } while (0)
#endif

  /*  
#if WIGNER6J_REGGE_CANONICALISE_ARRAY4
  DEBUG_V4SI(a1);
  DEBUG_V4SI(a2);
  DEBUG_V4SI(a3);

  DEBUG_V4SI(b1);
  DEBUG_V4SI(b2);
  DEBUG_V4SI(b3);
  DEBUG_V4SI(b4);
#endif
  */

  SWAP_TO_FIRST_LARGER(W6U32, a1, a2);
  SWAP_TO_FIRST_LARGER(W6U32, a2, a3);
  SWAP_TO_FIRST_LARGER(W6U32, a1, a2);

  SWAP_TO_FIRST_LARGER(W6U32, b1, b2);
  SWAP_TO_FIRST_LARGER(W6U32, b2, b3);
  SWAP_TO_FIRST_LARGER(W6U32, b3, b4);
  SWAP_TO_FIRST_LARGER(W6U32, b1, b2);
  SWAP_TO_FIRST_LARGER(W6U32, b2, b3);
  SWAP_TO_FIRST_LARGER(W6U32, b1, b2);

  /*
#if WIGNER6J_REGGE_CANONICALISE_ARRAY4
  DEBUG_V4SI(a1);
  DEBUG_V4SI(a2);
  DEBUG_V4SI(a3);

  DEBUG_V4SI(b1);
  DEBUG_V4SI(b2);
  DEBUG_V4SI(b3);
  DEBUG_V4SI(b4);
#endif
  */

  // We now have a1 >= a2 >= a3, and b1 >= b2 >= b3 >= b4

  W6U32 S32 = a3 - b1;
  W6U32 B32 = a3 - b2;
  W6U32 T32 = a3 - b3;
  W6U32 X32 = a3 - b4;
  W6U32 L32 = a2 - b4;
  W6U32 E32 = a1 - b4;

  /* This is where one would have transformed to use 64-bit variables.
   * We just do aliasing.
   */

#if 0 /* Go to 64-bit variables */
#if WIGNER6J_REGGE_CANONICALISE_ARRAY4
  W6U64 S = __builtin_ia32_pmovzxdq256(S32);
  W6U64 B = __builtin_ia32_pmovzxdq256(B32);
  W6U64 T = __builtin_ia32_pmovzxdq256(T32);
  W6U64 X = __builtin_ia32_pmovzxdq256(X32);
  W6U64 L = __builtin_ia32_pmovzxdq256(L32);
  W6U64 E = __builtin_ia32_pmovzxdq256(E32);
#else
  W6U64 S = S32;
  W6U64 B = B32;
  W6U64 T = T32;
  W6U64 X = X32;
  W6U64 L = L32;
  W6U64 E = E32;
#endif
#endif

  /*
  if (S > B || B > T || T > X || X > L || L > E)
    fprintf (stderr,"*** %d %d %d  %d %d %d | %d %d %d ; %d %d %d %d ; %d %d %d %d %d %d\n",
	     ja, jb, jc, jd, je, jf,
	     a1, a2, a3,
	     b1, b2, b3, b4,
	     (int) S, (int) B, (int) T, (int) X, (int) L, (int) E);
  */
  /*
  fprintf (stderr,"*** %d %d %d  %d %d %d | "
	   "%d %d %d ; %d %d %d %d ; "
	   "%d %d %d %d %d %d\n",
	   ja, jb, jc, jd, je, jf,
	   a1, a2, a3,
	   b1, b2, b3, b4,
	   (int) E, (int) L, (int) X, (int) T, (int) B, (int) S);
  */
#if defined(KEEP_MAX_RASCH_YU_INDEX) || WIGNER6J_REGGE_CANONICALISE_RETURN_INDEX
  /* As we have the values 2*j, we multiply by this 2 for all other
   * contributions.  Maximum a power j^6 is encountered in Eind.
   */

  /* SSE2 / AVX(2) does not have a 64x64 multiplier, only 32x32->64 bits.
   * We therefore do as much as possible of the calculations directly
   * in 32 bits, and only when needed extend the results to 64 bits.
   *
   */
  /*
# if WIGNER6J_REGGE_CANONICALISE_ARRAY4
#  define INIT_W6U64(x,val)  W6U32 x = { val, val, val, val };
# else
#  define INIT_W6U64(x,val)  W6U32 x = val;
# endif

  INIT_W6U64(c15,15);
  INIT_W6U64(c85,85);
  INIT_W6U64(c225,225);
  INIT_W6U64(c274,274);
  INIT_W6U64(c120,120);

  INIT_W6U64(c10,10);
  INIT_W6U64(c35,35);
  INIT_W6U64(c50,50);
  INIT_W6U64(c24,24);

  INIT_W6U64(c6,6);
  INIT_W6U64(c11,11);

  INIT_W6U64(c3,3);
  INIT_W6U64(c2,2);
  INIT_W6U64(c1,1);

  INIT_W6U64(c30,30);
  INIT_W6U64(c360,360);
  INIT_W6U64(c720,720);
  */
  int i;
#if WIGNER6J_REGGE_CANONICALISE_ARRAY4
  uint32_t __E[4] __attribute__ ((aligned (32)));  *(v4si *) &__E[0] = E32;
  uint32_t __L[4] __attribute__ ((aligned (32)));  *(v4si *) &__L[0] = L32;
  uint32_t __X[4] __attribute__ ((aligned (32)));  *(v4si *) &__X[0] = X32;
  uint32_t __T[4] __attribute__ ((aligned (32)));  *(v4si *) &__T[0] = T32;
  uint32_t __B[4] __attribute__ ((aligned (32)));  *(v4si *) &__B[0] = B32;
  uint32_t __S[4] __attribute__ ((aligned (32)));  *(v4si *) &__S[0] = S32;

  uint64_t __indexm720[4] __attribute__ ((aligned (32)));

  for (i = 0; i < 4; i++)
#else
  uint32_t __E[1] = { E32 };
  uint32_t __L[1] = { L32 };
  uint32_t __X[1] = { X32 };
  uint32_t __T[1] = { T32 };
  uint32_t __B[1] = { B32 };
  uint32_t __S[1] = { S32 };
  
  uint64_t __indexm720[1];

  for (i = 0; i < 1; i++)
#endif
    {
      uint32_t S = __S[i];
      uint32_t B = __B[i];
      uint32_t T = __T[i];
      uint32_t X = __X[i];
      uint32_t L = __L[i];
      uint32_t E = __E[i];

      /*  
  printf("ELXTBS: %" PRIu32 "  %" PRIu32 "  %" PRIu32 ""
	 "  %" PRIu32 "  %" PRIu32 "  %" PRIu32 "  \n",
	 E, L, X, T, B, S);
      */
  uint64_t Eind = E * (120 + E * (274 + (uint64_t) (E * (225 + E * (85 + E * (15 + E))))));
  uint64_t Lind = L * (24 + (uint64_t) (L * (50 + L * (35 + L * (10 + L)))));
  uint64_t Xind = X * (6 + X * (11 + X * (6 + X)));
  uint64_t Tind = T * (2 + T * (3 + T));
  uint64_t Bind = B * (1 + B);
  uint64_t Sind = S;
  __indexm720[i] = (/* */  Eind +
		     6   * Lind +
		     30  * Xind +
		     120 * Tind +
		     360 * Bind +
		     720 * Sind);
  /*
  printf("ELXTBS: %" PRIu64 "  %" PRIu64 "  %" PRIu64 ""
	 "  %" PRIu32 "  %" PRIu32 "  %" PRIu32 "  \n",
	 Eind, Lind, Xind, Tind, Bind, Sind);
  */
    }
  
#if WIGNER6J_REGGE_CANONICALISE_ARRAY4
  union wigner6j_type_pun_v4di_uint64_t index;
  /*W6U64 index;*/
  /*
  int __indexm720[4] __attribute__ ((aligned (32)));
  *(v4si *) &__indexm720[0] = indexm720;
  */
  /*uint64_t __index[4] __attribute__ ((aligned (32)));*/

  for (i = 0; i < 4; i++)
    {
      //printf ("%lld\n",__index[i]);
      /*__index*/index._u64[i] = __indexm720[i] / 720 + 1;
      //printf ("%lld\n",__index[i]);
    }

  /*index = *(v4di *) &__index[0];*/
  
#else
  uint64_t index = __indexm720[0] / 720 + 1;
  //printf("idx: %" PRIu64 "  %" PRIu64 "\n",__indexm720[0],index);
#endif

  /*
  fprintf (stderr,"%d %d %d %d %d %d -> %d\n",
	   Eind, Lind, Xind, Tind, Bind, Sind, index);
  */
# if defined(KEEP_MAX_RASCH_YU_INDEX)
#  if !WIGNER6J_REGGE_CANONICALISE_ARRAY4
  if (index > _max_RaschYu_index)
    _max_RaschYu_index = index;
#  endif
# endif
# if WIGNER6J_REGGE_CANONICALISE_RETURN_INDEX
  /*
  printf ("%016lld\n", (long long) index);
  */
#  if WIGNER6J_REGGE_CANONICALISE_ARRAY4
  *key = index._v4di;
#  else
  return index;
#  endif
# endif
#endif

  /* When using a hash table, we do not need the values ordered in an
   * index.
   */

  /*
  fprintf (stderr,"-> %016llx\n",
	   (E) | (L << 8) | (X << 16) | (T << 24) |
	   (((uint64_t) B) << 32) | (((uint64_t) S) << 40));
  */

  /* When the need arises to retrieve a set of jn from the c14n
   * variables (e.g. to calculate the value of the symbol):
   *
   * -j3 + j4 + j5 == S
   *  j2 + j4 - j6 == B
   *  j1 + j5 - j6 == T
   *  j1 + j2 - j3 == X
   *  j1 - j5 + j6 == L
   *  j2 - j4 + j6 == E
   *
   * gives (do not do /2 when using 2j-notation):
   *
   * j1 = (L + T)/2
   * j2 = (B + E)/2
   * j3 = (B + E + L + T - 2*X)/2 = j1 + j2 - 2*X/2
   * j4 = (B + L + S - X)/2
   * j5 = (E + S + T - X)/2
   * j6 = (E + L + S - X)/2
   *
   * A potential problem is that this gives larger jx than necessary,
   * and therefore a sumbol which incurs more penalty (speed and error)
   * when directly calculated.  (Has not been investigated.)
   */
}

#undef two_j1
#undef two_j2
#undef two_j3
#undef two_j4
#undef two_j5
#undef two_j6

#undef W6U32
#undef W6U64

#undef SWAP_TO_FIRST_LARGER
#undef SHIFT_DOWN_1
#undef INIT_W6U64

#undef COLLECT_NEGATIVE
#undef COLLECT_TRIANGLE_TRIVIAL_ZERO
