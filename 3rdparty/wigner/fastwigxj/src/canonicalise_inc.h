
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

#ifndef __CANONICALISE_INC_H__
#define __CANONICALISE_INC_H__

#include <x86intrin.h>

#include "fastwigxj_config.h"

#include "fastwigxj_vector.h"

#define DEBUG_CANONICALISE 0

typedef double  v2df __attribute__ ((vector_size (16)));
typedef float   v4sf __attribute__ ((vector_size (16)));

typedef union v2dif_t
{
  v2di _i;
  v2df _f;
} v2dif;

typedef union v4sif_t
{
  v4si _i;
  v4sf _f;
} v4sif;

typedef union v2d4si_t
{
  v2di _2d;
  v4si _4s;
} v2d4si;

#if !FASTWIGXJ_HAVE_SSE4_1
  // Without sse4 we have to do comparision by explicit subtraction,
  // and the selection (blend) by and/or-ing with the result.
  // The code below does not work when the highest bit is used for data,
  // so we are effectively limited to two_j <= 63 (instead of 127).
#define W9C_KEEPMAX(jx,jy) do {						\
    v2di __v2di_1 = { 1, 1 };						\
    v2di __y_sub_x = jy - jx; /* high bit set if jx > jy */		\
    /* shift the bit down to the lowest bit position */			\
    v2di __x_gt_y = __builtin_ia32_psrlqi128(__y_sub_x,63);		\
    /* and make a bitmask with the inverse condition. */		\
    v2di __x_le_y = __x_gt_y - __v2di_1;				\
    /* use jx, but jy if __x_le_y */					\
    jx = jx + (__x_le_y & (/*jy - jx*/__y_sub_x));			\
  } while (0)
#define W9C_KEEPMIN(jx,jy) do {						\
    v2di __v2di_1 = { 1, 1 };						\
    v2di __y_sub_x = jx - jy; /* high bit set if jx < jy */		\
    /* shift the bit down to the lowest bit position */			\
    v2di __x_gt_y = __builtin_ia32_psrlqi128(__y_sub_x,63);		\
    /* and make a bitmask with the inverse condition. */		\
    v2di __x_le_y = __x_gt_y - __v2di_1;				\
    /* use jx, but jy if __x_le_y */					\
    jx = jx + (__x_le_y & (/*jy - jx*/__y_sub_x));			\
  } while (0)
#define W6C_KEEPMAXMIN(jx,jy) do {					\
    v4si __v4si_1 = { 1, 1, 1, 1 };					\
    v4si __y_sub_x = jy - jx; /* high bit set if jx > jy */		\
    /* shift the bit down to the lowest bit position */			\
    v4si __x_gt_y = __builtin_ia32_psrldi128(__y_sub_x,31);		\
    /* and make a bitmask with the inverse condition. */		\
    v4si __x_le_y = __x_gt_y - __v4si_1;				\
    /* use jx, but jy if __x_le_y */					\
    jx = jx + (__x_le_y & (/*jy - jx*/__y_sub_x));			\
    jy = jy - (__x_le_y & (/*jy - jx*/__y_sub_x));			\
  } while (0)
#if 0
  // Attempt with a 32-bit arithmetic shift followed by a shuffle
  // turned out slower.  Also did not manage to remove the
  // memory-operations.
#define W9C_KEEPMAX(jx,jy) do {						\
    v2d4si __y_sub_xds;							\
    __y_sub_xds._2d = jy - jx; /* high bit set if jx > jy */		\
    /* shift the bit down to the lowest bit position */			\
    v4si __x_gt_y_hi32 = __builtin_ia32_psradi128(__y_sub_xds._4s,31);	\
    v2d4si __x_gt_y;							\
    __x_gt_y._4s =							\
      __builtin_ia32_pshufd(__x_gt_y_hi32,				\
			    (1 << 0) | (1 << 2) | (3 << 4) | (3 << 6));	\
    /* use jx, but jy if __x_le_y */					\
    jx = jy - (__x_gt_y._2d & (/*jy - jx*/__y_sub_xds._2d));		\
  } while (0)
#endif
#else // for sse4
#define W9C_KEEPMAX(jx,jy) do {						\
    v2dif __x_gt_y;							\
    /* __builtin_ia32_pcmpgtq */					\
    __x_gt_y._i = (v2di) _mm_cmpgt_epi64((__m128i) jx,(__m128i) jy);	\
    v2dif __x, __y, __z;						\
    __x._i = jx;							\
    __y._i = jy;							\
    __z._f = __builtin_ia32_blendvpd(__y._f,__x._f,__x_gt_y._f);	\
    jx = __z._i;							\
  } while (0)
#define W9C_KEEPMIN(jx,jy) do {						\
    v2dif __y_gt_x;							\
    /* __builtin_ia32_pcmpgtq */					\
    __y_gt_x._i = (v2di) _mm_cmpgt_epi64((__m128i) jy,(__m128i) jx);	\
    v2dif __x, __y, __z;						\
    __x._i = jx;							\
    __y._i = jy;							\
    __z._f = __builtin_ia32_blendvpd(__y._f,__x._f,__y_gt_x._f);	\
    jx = __z._i;							\
  } while (0)
#define W6C_KEEPMAXMIN(jx,jy) do {					\
    v4sif __x_gt_y;							\
    /* __builtin_ia32_pcmpgtd128 */					\
    __x_gt_y._i = (v4si) _mm_cmpgt_epi32((__m128i) jx,(__m128i) jy);	\
    v4sif __x, __y, __z, __w;						\
    __x._i = jx;							\
    __y._i = jy;							\
    __z._f = __builtin_ia32_blendvps(__y._f,__x._f,__x_gt_y._f);	\
    jx = __z._i;							\
    __w._f = __builtin_ia32_blendvps(__x._f,__y._f,__x_gt_y._f);	\
    jy = __w._i;							\
  } while (0)
#endif

/* split from x:lo and x:hi to x:lo and y:lo */
#define W9C_HILO_SPLIT(x,y) do {			\
    v2dif __x;						\
    v2dif __y;						\
    __x._i = x;						\
    /* __builtin_ia32_unpckhpd */			\
    __y._f = (v2df) _mm_unpackhi_pd(__x._f, __x._f);	\
    y = __y._i;						\
  } while ( 0 )

#if FASTWIGXJ_HAVE_AVX2

typedef double  v4df __attribute__ ((vector_size (32)));

typedef union v4dif_t
{
  v4di _i;
  v4df _f;
} v4dif;

#define W9C_KEEPMAX_v4di(jx,jy) do {					\
    v4dif __x_gt_y;							\
    /* __builtin_ia32_pcmpgtq256 */					\
    __x_gt_y._i = _mm256_cmpgt_epi64(jx,jy);				\
    v4dif __x, __y, __z;						\
    __x._i = jx;							\
    __y._i = jy;							\
    __z._f = __builtin_ia32_blendvpd256(__y._f,__x._f,__x_gt_y._f);	\
    jx = __z._i;							\
  } while (0)

typedef union v42di_t
{
  v4di _4d;
  v2di _2d;
} v42di;

#if DEBUG_CANONICALISE
#define DEBUG_V4DI(x) do {						\
    long long int __reg_to_ymm[4] __attribute__ ((aligned (32)));	\
    *(v4di *) &__reg_to_ymm[0] = x;					\
    printf ("%-20s %016llx %016llx %016llx %016llx\n",			\
	    #x,								\
	    __reg_to_ymm[3],__reg_to_ymm[2],				\
	    __reg_to_ymm[1],__reg_to_ymm[0]);				\
  } while ( 0 )
#else
#define DEBUG_V4DI(x)    (void) x
#endif

#endif/*FASTWIGXJ_HAVE_AVX2*/

#if DEBUG_CANONICALISE
#define DEBUG_V2DI(x) do {						\
    long long int __reg_to_xmm[2] __attribute__ ((aligned (16)));	\
    *(v2di *) &__reg_to_xmm[0] = x;					\
    printf ("%-20s %016llx %016llx\n",					\
	    #x,								\
	    __reg_to_xmm[1],__reg_to_xmm[0]);				\
  } while ( 0 )

#define DEBUG_V4SI(x) do {						\
    int __reg_to_xmm[4] __attribute__ ((aligned (16)));			\
    *(v4si *) &__reg_to_xmm[0] = x;					\
    printf ("%-20s %08x %08x %08x %08x\n",				\
	    #x,								\
	    __reg_to_xmm[3],__reg_to_xmm[2],				\
	    __reg_to_xmm[1],__reg_to_xmm[0]);				\
  } while ( 0 )

#define DEBUG_V2DI_9(x) do {						\
    long long int __reg_to_xmm[2] __attribute__ ((aligned (16)));	\
    *(v2di *) &__reg_to_xmm[0] = x;					\
    printf ("%-20s "							\
	    "%02llx %02llx %02llx %02llx %02llx %02llx %02llx %02llx %02llx " \
	    "%02llx %02llx %02llx %02llx %02llx %02llx %02llx %02llx %02llx\n", \
	    #x,								\
	    (__reg_to_xmm[1] >> 57) & 0x7f,				\
	    (__reg_to_xmm[1] >> 50) & 0x7f,				\
	    (__reg_to_xmm[1] >> 43) & 0x7f,				\
	    (__reg_to_xmm[1] >> 36) & 0x7f,				\
	    (__reg_to_xmm[1] >> 29) & 0x7f,				\
	    (__reg_to_xmm[1] >> 22) & 0x7f,				\
	    (__reg_to_xmm[1] >> 15) & 0x7f,				\
	    (__reg_to_xmm[1] >>  8) & 0x7f,				\
	    (__reg_to_xmm[1] >>  1) & 0x7f,				\
	    (__reg_to_xmm[0] >> 57) & 0x7f,				\
	    (__reg_to_xmm[0] >> 50) & 0x7f,				\
	    (__reg_to_xmm[0] >> 43) & 0x7f,				\
	    (__reg_to_xmm[0] >> 36) & 0x7f,				\
	    (__reg_to_xmm[0] >> 29) & 0x7f,				\
	    (__reg_to_xmm[0] >> 22) & 0x7f,				\
	    (__reg_to_xmm[0] >> 15) & 0x7f,				\
	    (__reg_to_xmm[0] >>  8) & 0x7f,				\
	    (__reg_to_xmm[0] >>  1) & 0x7f);				\
  } while ( 0 )

#define DEBUG_NL         printf ("\n")
#else
#define DEBUG_V2DI(x)    (void) x
#define DEBUG_V2DI_9(x)  (void) x
#define DEBUG_NL         do { } while (0)
#endif

#endif/*__CANONICALISE_INC_H__*/
