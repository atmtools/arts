
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

#include <stdio.h>
#include <stdint.h>
#include "canonicalise_inc.h"
#include "fastwigxj_struct.h"

/* Note:
 *
 * v2di     is well defined in type and size
 * __m128i  is not well defined, and does not work with arithmetic
 *          operations.  (But seems to be helpful for the casts
 *          in calls to _mm_shuffle_pd).
 *
 * the __builtin_ia32 (gcc style??) seems to be associated with v2di
 * style, while _mm (intel style intrinsics) use the less descriptive
 * types.
 */

#if !FASTWIGXJ_HAVE_AVX2
#define V2DI v2di
#define V4DI v2di
#else
#define V2DI __m128i
#define V4DI v4di
#endif

inline uint64_t wigner9j_canonicalise(const int *two_jv)
{
  long long int reg_to_xmm[6] __attribute__ ((aligned (16)));

  V2DI j128_32_dcba = *(const V2DI *) &(two_jv[0]);
  V2DI j128_32_hgfe = *(const V2DI *) &(two_jv[4]);
  V2DI j128_32_xxxi = *(const V2DI *) &(two_jv[8]);
  
  DEBUG_NL;
  DEBUG_V2DI(j128_32_dcba);
  DEBUG_V2DI(j128_32_hgfe);
  DEBUG_V2DI(j128_32_xxxi);

  // a b c d e f g h i
  // a d g b e h c f i
  // a       e       i
  //       b       f
  //   d       h
  //     g
  //             c    

  int64_t j64_32_10 = (0xffffffffll << 32);
  int64_t j64_32_01 = (0xffffffffll <<  0);

  v2di j128_32_1010 = { j64_32_10, j64_32_10 };
  v2di j128_32_0101 = { j64_32_01, j64_32_01 };

  __m128i j128_32__c_a = j128_32_dcba & j128_32_0101;
  __m128i j128_32_d_b_ = j128_32_dcba & j128_32_1010;
  __m128i j128_32__g_e = j128_32_hgfe & j128_32_0101;
  __m128i j128_32_h_f_ = j128_32_hgfe & j128_32_1010;
  __m128i j128_32____i = j128_32_xxxi & j128_32_0101;
  
  DEBUG_NL;
  DEBUG_V2DI(j128_32__c_a);
  DEBUG_V2DI(j128_32_d_b_);
  DEBUG_V2DI(j128_32__g_e);
  DEBUG_V2DI(j128_32_h_f_);
  DEBUG_V2DI(j128_32____i);

  /* __builtin_ia32_punpcklqdq128 and __builtin_ia32_punpckhqdq128 */
#define UNP_LOLO(a,b) _mm_unpacklo_epi64(a,b)
#define UNP_HIHI(a,b) _mm_unpackhi_epi64(a,b)

#define SHL128(a,s) __builtin_ia32_psllqi128(a,s)
#define SHR128(a,s) __builtin_ia32_psrlqi128(a,s)

  V2DI j__a________a______ =
    SHL128(UNP_LOLO(j128_32__c_a,j128_32__c_a),43);
  V2DI j_____b________d___ = SHR128(j128_32_d_b_,32-22);
  V2DI j________c________g =
    SHL128(UNP_HIHI(j128_32__c_a,j128_32__g_e),1);
  V2DI j_____e________e___ =
    SHL128(UNP_LOLO(j128_32__g_e,j128_32__g_e),22);
  V2DI j________f________h = SHR128(j128_32_h_f_,32-1);
  V2DI j________i________i =
    SHL128(UNP_LOLO(j128_32____i,j128_32____i),1);
  
  DEBUG_NL;
  DEBUG_V2DI_9(j__a________a______);
  DEBUG_V2DI_9(j_____b________d___);
  DEBUG_V2DI_9(j________c________g);
  DEBUG_V2DI_9(j_____e________e___);
  DEBUG_V2DI_9(j________f________h);
  DEBUG_V2DI_9(j________i________i);
  
#define SWAP_HILO(a) (__m128i) _mm_shuffle_pd((__m128d) (__m128i) a,(__m128d) (__m128i) a,_MM_SHUFFLE2(0,1))
  //#define SWAP_HILO(a) _mm_shuffle_epi32(a,_MM_SHUFFLE2(0,1))

  V2DI j__d________b______ = SHL128(SWAP_HILO(j_____b________d___),21);
  V2DI j__g________c______ = SHL128(SWAP_HILO(j________c________g),42);
  V2DI j_____h________f___ = SHL128(SWAP_HILO(j________f________h),21);
  
  DEBUG_NL;
  DEBUG_V2DI_9(j__d________b______);
  DEBUG_V2DI_9(j__g________c______);
  DEBUG_V2DI_9(j_____h________f___);
  
  V2DI v128j__a__b__c/*__a__d__g*/ =
    j__a________a______ | j_____b________d___ | j________c________g;
  V2DI v128j__d__e__f/*__b__e__h*/ =
    j__d________b______ | j_____e________e___ | j________f________h;
  V2DI v128j__g__h__i/*__c__f__i*/ =
    j__g________c______ | j_____h________f___ | j________i________i;

#ifndef WIGNER9J_NO_OVERFLOW_CHECK
  /* Check for 2j overflowing the 7 bits available. */

  int64_t j110110110_ =
    (0x7fll << 57) | (0x7fll << 50) |
    (0x7fll << 36) | (0x7fll << 29) |
    (0x7fll << 15) | (0x7fll <<  8);

  V2DI j110110110 = { j110110110_, j110110110_ };

  V2DI joverflow =
    (v128j__a__b__c | v128j__d__e__f | v128j__g__h__i) & j110110110;

  *(v2di *) &reg_to_xmm[0] = joverflow;

  if (reg_to_xmm[0])
    {
      /* We would overflow while c14n. */
      /* Return a special non-existent  key. */
      return (uint64_t) -4;
      exit(1);
    }
#endif

#ifndef WIGNER9J_NO_TRIVIAL0_CHECK
  /* Trivial-0 check. */

  /* First see if any of the sums are odd.
   *
   * As the values are separated by 0s, overflows will not spill
   * into the next line values.
   */

  DEBUG_V2DI(v128j__a__b__c);
  DEBUG_V2DI(v128j__d__e__f);
  DEBUG_V2DI(v128j__g__h__i);

  v2di jsum_adg = v128j__a__b__c + v128j__d__e__f + v128j__g__h__i;

  DEBUG_V2DI(jsum_adg);

  /* If any sum was odd, bit 1 (, 22, 43) of jodd will be 1. */

  /* And then check the triangle conditions. */

  int64_t jguard_ = (1ll << 63) | (1ll << 42) | (1ll << 21);

  v2di jguard = { jguard_, jguard_ };

  v2di jsum_ad  = v128j__a__b__c + v128j__d__e__f;
  v2di jdiff_ad = (jguard + v128j__a__b__c) - v128j__d__e__f;
  v2di jdiff_da = (jguard + v128j__d__e__f) - v128j__a__b__c;

  DEBUG_V2DI(jsum_ad);
  DEBUG_V2DI(jdiff_ad);
  DEBUG_V2DI(jdiff_da);

  /* a > d => a-d > 0 , d-a < 0
   *
   * guard                .. 1000000 0000000 0000000
   * a                    .. 0000000 0000000 0011111
   * d                    .. 0000000 0000000 0000000
   * diff_ad              .. 1000000 0000000 0011111 valid to test
   * diff_da              .. 0111111 1111111 1100001 not valid test
   *
   * diff_ad 2compl       .. 0111111 1111111 1100001
   * diff_da 2compl       .. 1000000 0000000 0011111
   *
   * g > a-d :
   *
   * g                    .. 0000000 0000000 1000000
   * guard+g              .. 1000000 0000000 1000000
   * g_minus_abs_diff_ad  .. 0000000 0000000 0100001 valid, g ok
   * g_minus_abs_diff_da  .. 1000000 0000000 1011111 not valid test
   *
   * g < a-d :
   *
   * g                    .. 0000000 0000000 0000100
   * guard+g              .. 1000000 0000000 0000100
   * g_minus_abs_diff_ad  .. 1111111 1111111 1100101 valid, g bad
   * g_minus_abs_diff_da  .. 1000000 0000000 0100011 not valid test
   */

  /* We want to check that g <= a+d.  If g > a+d = sum_ad, then we are
   * trivially 0.  That is easily done by subtracting g from sum_ad.
   * If the subtraction starts to borrow bits, then sum_ad is smaller,
   * so the test is if a higher bit than we normally have numbers in
   * got set.  It is no problem if we start to borrow from a higher
   * line (like h > b+e), as it is enough that one fails the test, for
   * the symbol to be trivially 0.
   */

  v2di jsum_ad_minus_g = jsum_ad - v128j__g__h__i;

  DEBUG_V2DI(jsum_ad_minus_g);

  /* We want to check that g >= |a-d|.  This means that g >= a-d and
   * also g > -(a-d).  However, to not have the subtraction of a-d (or
   * d-a) borrow bits from the next value to check, we have introduced
   * guard bits before doing the first subtraction.
   */

  v2di jg_minus_abs_diff_ad = (jguard + v128j__g__h__i) - jdiff_ad;
  v2di jg_minus_abs_diff_da = (jguard + v128j__g__h__i) - jdiff_da;

  DEBUG_V2DI(jg_minus_abs_diff_ad);
  DEBUG_V2DI(jg_minus_abs_diff_da);

  /* One of the differences a-d and d-a will be negative, i.e.
   * have lots of 1s up till (not including) the guard bit.
   * But since we subtract it, we will just be dumping at lot of
   * 0s up to the guard bit.  The guard bit also becomes 0.
   */

  v2di jtriangle_check =
    jsum_ad_minus_g | jg_minus_abs_diff_ad | jg_minus_abs_diff_da;

  DEBUG_V2DI(jtriangle_check);
  DEBUG_V2DI_9(jtriangle_check);

  v2di jtriangle_check_shft = __builtin_ia32_psrlqi128(jtriangle_check,19);

  DEBUG_V2DI(jtriangle_check_shft);

  /* Shift the badness down from bits 62, 41, 20 to 43, 22 and 1,
   * and merge with the parity check.
   */

  v2di jcollect_trivial0 =
    jtriangle_check_shft | jsum_adg;

  DEBUG_V2DI(jcollect_trivial0);

  /* Collect the badness bits. */

  v2di jtrivial0_bit1 = jcollect_trivial0 |
    __builtin_ia32_psrlqi128(jcollect_trivial0,21) |
    __builtin_ia32_psrlqi128(jcollect_trivial0,42);

  DEBUG_V2DI(jtrivial0_bit1);

  // 2 values now in jtrivial0_bit1:lo and jtrivial0_bit1:hi
  // keep i jtrivial0_bit1:lo and get the other into jtrivial0_bit1_2:lo

  v2di jtrivial0_bit1_w2;

  W9C_HILO_SPLIT(jtrivial0_bit1, jtrivial0_bit1_w2);

  DEBUG_V2DI(jtrivial0_bit1_w2);

  jtrivial0_bit1 = jtrivial0_bit1 | jtrivial0_bit1_w2;

  DEBUG_V2DI(jtrivial0_bit1);

  *(v2di *) &reg_to_xmm[0] = jtrivial0_bit1;

  if (reg_to_xmm[0] & (1 << 1))
    {
      /* We are trivially zero. */
      STATS_9J->_trivial0++;
      /* Return the special trivial-0 key. */
      return (uint64_t) -2;
    }

  /* End of trivial-0 check. */
#endif

  DEBUG_NL;
  DEBUG_V2DI_9(v128j__a__b__c);
  DEBUG_V2DI_9(v128j__d__e__f);
  DEBUG_V2DI_9(v128j__g__h__i);
  DEBUG_NL;

  int64_t j111000000_ = (0x7fll << 57) | (0x7fll << 50) | (0x7fll << 43);
  int64_t j000111000_ = (0x7fll << 36) | (0x7fll << 29) | (0x7fll << 22);
  int64_t j000000111_ = (0x7fll << 15) | (0x7fll <<  8) | (0x7fll <<  1);

  V2DI j111000000 = { j111000000_, j111000000_ };
  V2DI j000111000 = { j000111000_, j000111000_ };
  V2DI j000000111 = { j000000111_, j000000111_ };

#if !FASTWIGXJ_HAVE_AVX2
#else
  v4di v256j111000000 = { j111000000_, j111000000_, j111000000_, j111000000_ };
  v4di v256j000111000 = { j000111000_, j000111000_, j000111000_, j000111000_ };
  v4di v256j000000111 = { j000000111_, j000000111_, j000000111_, j000000111_ };

#endif

#if !FASTWIGXJ_HAVE_AVX2
#define W9C_SHIFTUP_6 do {					\
    j_a__b__c_ = __builtin_ia32_psllqi128(j__a__b__c,7);	\
    ja__b__c__ = __builtin_ia32_psllqi128(j__a__b__c,14);	\
    j_d__e__f_ = __builtin_ia32_psllqi128(j__d__e__f,7);	\
    jd__e__f__ = __builtin_ia32_psllqi128(j__d__e__f,14);	\
    j_g__h__i_ = __builtin_ia32_psllqi128(j__g__h__i,7);	\
    jg__h__i__ = __builtin_ia32_psllqi128(j__g__h__i,14);	\
  } while (0)
#else
#define W9C_SHIFTUP_6 do {					\
    j_a__b__c_ = __builtin_ia32_psllqi256(j__a__b__c,7);	\
    ja__b__c__ = __builtin_ia32_psllqi256(j__a__b__c,14);	\
    j_d__e__f_ = __builtin_ia32_psllqi256(j__d__e__f,7);	\
    jd__e__f__ = __builtin_ia32_psllqi256(j__d__e__f,14);	\
    j_g__h__i_ = __builtin_ia32_psllqi256(j__g__h__i,7);	\
    jg__h__i__ = __builtin_ia32_psllqi256(j__g__h__i,14);	\
  } while (0)
#endif

#if !FASTWIGXJ_HAVE_AVX2
#define W9C_SWAPCOL12(jxxxyyyzzz) do {				\
    jxxxyyyzzz =						\
      __builtin_ia32_psrlqi128(jxxxyyyzzz & j111000000,21) |	\
      __builtin_ia32_psllqi128(jxxxyyyzzz & j000111000,21) |	\
      /* */                   (jxxxyyyzzz & j000000111);	\
  } while (0)

#define W9C_SWAPCOL12_3 do {			\
    W9C_SWAPCOL12(j__a__b__c);			\
    W9C_SWAPCOL12(j__d__e__f);			\
    W9C_SWAPCOL12(j__g__h__i);			\
  } while (0)

#define W9C_SWAPCOL23(jxxxyyyzzz) do {				\
    jxxxyyyzzz =						\
      __builtin_ia32_psrlqi128(jxxxyyyzzz & j000111000,21) |	\
      __builtin_ia32_psllqi128(jxxxyyyzzz & j000000111,21) |	\
      /* */                   (jxxxyyyzzz & j111000000);	\
  } while (0)

#define W9C_SWAPCOL23_3 do {			\
    W9C_SWAPCOL23(j__a__b__c);			\
    W9C_SWAPCOL23(j__d__e__f);			\
    W9C_SWAPCOL23(j__g__h__i);			\
  } while (0)
#else
#define W9C_SWAPCOL12_256(jxxxyyyzzz) do {				\
    jxxxyyyzzz =						\
      __builtin_ia32_psrlqi256(jxxxyyyzzz & v256j111000000,21) |	\
      __builtin_ia32_psllqi256(jxxxyyyzzz & v256j000111000,21) |	\
      /* */                   (jxxxyyyzzz & v256j000000111);	\
  } while (0)

#define W9C_SWAPCOL12_3 do {			\
    W9C_SWAPCOL12_256(j__a__b__c);			\
    W9C_SWAPCOL12_256(j__d__e__f);			\
    W9C_SWAPCOL12_256(j__g__h__i);			\
  } while (0)

#define W9C_SWAPCOL23_256(jxxxyyyzzz) do {				\
    jxxxyyyzzz =						\
      __builtin_ia32_psrlqi256(jxxxyyyzzz & v256j000111000,21) |	\
      __builtin_ia32_psllqi256(jxxxyyyzzz & v256j000000111,21) |	\
      /* */                   (jxxxyyyzzz & v256j111000000);	\
  } while (0)

#define W9C_SWAPCOL23_3 do {			\
    W9C_SWAPCOL23_256(j__a__b__c);			\
    W9C_SWAPCOL23_256(j__d__e__f);			\
    W9C_SWAPCOL23_256(j__g__h__i);			\
  } while (0)
#endif


#if !FASTWIGXJ_HAVE_AVX2
  V2DI j__a__b__c = v128j__a__b__c;
  V2DI j__d__e__f = v128j__d__e__f;
  V2DI j__g__h__i = v128j__g__h__i;

#else

// The above did it in the following order

// XYZ
// YXZ
// YZX
// ZYX
// ZXY
// XZY

// We now combine them again

// XYZ : XZY
// YXZ : ZXY
// YZX : ZYX

  /*
  __m128i q128 = *(__m128i *) &reg_to_xmm[0];
  __m256i q256 = _mm256_castsi128_si256(q128);
  */
  
  DEBUG_V2DI(v128j__a__b__c);
  DEBUG_V2DI(v128j__d__e__f);
  DEBUG_V2DI(v128j__g__h__i);
  
  /* Copies where we have shifted around parts 2 and 3. */

  __m128i w128j__a__b__c = v128j__a__b__c;
  __m128i w128j__d__e__f = v128j__d__e__f;
  __m128i w128j__g__h__i = v128j__g__h__i;

#define W9C_SWAPCOL23(jxxxyyyzzz) do {				\
    jxxxyyyzzz =						\
      __builtin_ia32_psrlqi128(jxxxyyyzzz & j000111000,21) |	\
      __builtin_ia32_psllqi128(jxxxyyyzzz & j000000111,21) |	\
      /* */                   (jxxxyyyzzz & j111000000);	\
  } while (0)

  W9C_SWAPCOL23(w128j__a__b__c);
  W9C_SWAPCOL23(w128j__d__e__f);
  W9C_SWAPCOL23(w128j__g__h__i);

  DEBUG_V2DI(w128j__a__b__c);
  DEBUG_V2DI(w128j__d__e__f);
  DEBUG_V2DI(w128j__g__h__i);

  /* Join the two 128-bit pieces into 256 bit pieces. */

#define W9C_JOIN128_256(res,lo,hi) do {			\
    __m256i qq;						\
    qq = _mm256_castsi128_si256(lo);			\
    res= _mm256_inserti128_si256(qq, hi, 1);		\
  } while(0)

  __m256i j__a__b__c;
  __m256i j__d__e__f;
  __m256i j__g__h__i;

  W9C_JOIN128_256(j__a__b__c, v128j__a__b__c, w128j__a__b__c);
  W9C_JOIN128_256(j__d__e__f, v128j__d__e__f, w128j__d__e__f);
  W9C_JOIN128_256(j__g__h__i, v128j__g__h__i, w128j__g__h__i);

  DEBUG_V4DI(j__a__b__c);
  DEBUG_V4DI(j__d__e__f);
  DEBUG_V4DI(j__g__h__i);
#endif

  V4DI j_a__b__c_;
  V4DI ja__b__c__;
  V4DI j_d__e__f_;
  V4DI jd__e__f__;
  V4DI j_g__h__i_;
  V4DI jg__h__i__;

#if !FASTWIGXJ_HAVE_AVX2
  v2di jsign = { 1, 1 };

#define OR_SIGN_A
#define OR_SIGN_B | jsign

#define W9C_KEEPMAX_v4di W9C_KEEPMAX
#else
  v4di jsignA = { 0, 0, 1, 1 };
  v4di jsignB = { 1, 1, 0, 0 };

#define OR_SIGN_A | jsignA
#define OR_SIGN_B | jsignB
#endif

  W9C_SHIFTUP_6;
  
  V4DI jadgbehcfi = ja__b__c__ | j_d__e__f_ | j__g__h__i OR_SIGN_A; // +
  V4DI jmin1 = jadgbehcfi;
  V4DI jagdbhecif = ja__b__c__ | j_g__h__i_ | j__d__e__f OR_SIGN_B; // -
  V4DI jmin2 = jagdbhecif;

  V4DI jdagebhfci = jd__e__f__ | j_a__b__c_ | j__g__h__i OR_SIGN_B; // -
  W9C_KEEPMAX_v4di(jmin1, jdagebhfci);
  V4DI jdgaehbfic = jd__e__f__ | j_g__h__i_ | j__a__b__c OR_SIGN_A; // +
  W9C_KEEPMAX_v4di(jmin2, jdgaehbfic);

  V4DI jgadhbeicf = jg__h__i__ | j_a__b__c_ | j__d__e__f OR_SIGN_A; // +
  W9C_KEEPMAX_v4di(jmin1, jgadhbeicf);
  V4DI jgdahebifc = jg__h__i__ | j_d__e__f_ | j__a__b__c OR_SIGN_B; // -
  W9C_KEEPMAX_v4di(jmin2, jgdahebifc);

#define W9C_GENERATE_KEEP(sign1,sign2) do {			\
    jadgbehcfi = ja__b__c__ | j_d__e__f_ | j__g__h__i sign1;	\
    W9C_KEEPMAX_v4di(jmin1, jadgbehcfi);			\
    jagdbhecif = ja__b__c__ | j_g__h__i_ | j__d__e__f sign2;	\
    W9C_KEEPMAX_v4di(jmin2, jagdbhecif);			\
    jdagebhfci = jd__e__f__ | j_a__b__c_ | j__g__h__i sign2;	\
    W9C_KEEPMAX_v4di(jmin1, jdagebhfci);			\
    jdgaehbfic = jd__e__f__ | j_g__h__i_ | j__a__b__c sign1;	\
    W9C_KEEPMAX_v4di(jmin2, jdgaehbfic);			\
    jgadhbeicf = jg__h__i__ | j_a__b__c_ | j__d__e__f sign1;	\
    W9C_KEEPMAX_v4di(jmin1, jgadhbeicf);			\
    jgdahebifc = jg__h__i__ | j_d__e__f_ | j__a__b__c sign2;	\
    W9C_KEEPMAX_v4di(jmin2, jgdahebifc);			\
  } while (0)

  // Then mix around xxxyyyzzz -> yyyxxxzzz  [ AVX2: xxxzzzyyy -> zzzxxxyyy ]

#define WIGNER9J_LOOP_VARIATIONS 0
  /* These variations had no effect on speed.  Main idea is to have
   * explicit loops in the code, to save instruction cache.
   *
   * What was 'worse': the compiler unrolled them! when the
   * conditional became too simple.
   *
   * And when not (avx2) it got worse...  c14n time, but overall look-up
   * time improved.  Caching effect?
   */

#if WIGNER9J_LOOP_VARIATIONS
#if !FASTWIGXJ_HAVE_AVX2
#undef OR_SIGN_A
#undef OR_SIGN_B
  v2di jsignB = jsign;
  v2di jsignA = { 0, 0 };

#define OR_SIGN_A | jsignA
#define OR_SIGN_B | jsignB
#endif
  int i;
  for (i = 0; i < 5; i++)
    {
      switch (i)
	{
	case 0:
	case 2:
	case 4:
	  W9C_SWAPCOL12_3;
	  break;
	case 1:
	case 3:
	  W9C_SWAPCOL23_3;
	  break;
	}

      W9C_SHIFTUP_6;
      /*
      if (i & 1)
	{
	  W9C_GENERATE_KEEP(OR_SIGN_A,OR_SIGN_B);
	}
      else
      {*/
	  W9C_GENERATE_KEEP(OR_SIGN_B,OR_SIGN_A);
      /*}*/
      
      V4DI jsign_tmp = jsignA;
      jsignA = jsignB;
      jsignB = jsign_tmp;

      // oups - loop got unrolled...      
    }
#else
  W9C_SWAPCOL12_3;
  W9C_SHIFTUP_6;
  W9C_GENERATE_KEEP(OR_SIGN_B,OR_SIGN_A);

  // yyyxxxzzz -> yyyzzzxxx  [ AVX2: zzzxxxyyy -> zzzyyyxxx ]

  W9C_SWAPCOL23_3;
  W9C_SHIFTUP_6;
  W9C_GENERATE_KEEP(OR_SIGN_A,OR_SIGN_B);

#if !FASTWIGXJ_HAVE_AVX2
  // yyyzzzxxx -> zzzyyyxxx

  W9C_SWAPCOL12_3;
  W9C_SHIFTUP_6;
  W9C_GENERATE_KEEP(OR_SIGN_B,OR_SIGN_A);

  // zzzyyyxxx -> zzzxxxyyy

  W9C_SWAPCOL23_3;
  W9C_SHIFTUP_6;
  W9C_GENERATE_KEEP(OR_SIGN_A,OR_SIGN_B);

  // zzzxxxyyy -> xxxzzzyyy

  W9C_SWAPCOL12_3;
  W9C_SHIFTUP_6;
  W9C_GENERATE_KEEP(OR_SIGN_B,OR_SIGN_A);
#endif
#endif

  W9C_KEEPMAX_v4di(jmin1, jmin2);

  V2DI v128jmin1;
  V2DI v128jmin2;

#if !FASTWIGXJ_HAVE_AVX2
  v128jmin1 = jmin1;

  // 2 min values is now in v128jmin1:lo and v128jmin1:hi
#else /* FASTWIGXJ_HAVE_AVX2 */

  DEBUG_V4DI(jmin1);

  // 4 min values are now in jmin1:0, jmin1:1, jmin1:2, and jmin1:3
  // keep in jmin1:0 and 1, and get the other into jmin2:0 and 1

  v128jmin2 = _mm256_extracti128_si256(jmin1, 1);
  v128jmin1 = _mm256_castsi256_si128(jmin1);

  DEBUG_V2DI(v128jmin1);
  DEBUG_V2DI(v128jmin2);

  W9C_KEEPMAX(v128jmin1, v128jmin2);

  DEBUG_V2DI(v128jmin1);
#endif/* FASTWIGXJ_HAVE_AVX2 */

  // 2 min values is now in jmin1:lo and jmin1:hi
  // keep i jmin1:lo and get the other into jmin2:lo

  W9C_HILO_SPLIT(v128jmin1, v128jmin2);

  DEBUG_V2DI(v128jmin1);
  DEBUG_V2DI(v128jmin2);

  W9C_KEEPMAX(v128jmin1, v128jmin2);

  DEBUG_V2DI(v128jmin1);

  // minimum now in jmin1:lo

  *(v2di *) &reg_to_xmm[0] = v128jmin1;

  return (uint64_t) reg_to_xmm[0];
}
