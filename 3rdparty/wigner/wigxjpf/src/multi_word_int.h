
/* Copyright 2015 Haakan T. Johansson */

/*  This file is part of WIGXJPF.
 *
 *  WIGXJPF is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  WIGXJPF is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with WIGXJPF.  If not, see
 *  <http://www.gnu.org/licenses/>.
 */

#ifndef __WIGXJPF_MULTI_WORD_INT_H__
#define __WIGXJPF_MULTI_WORD_INT_H__

#include "wigxjpf_config.h"
#include "wigxjpf_error.h"

#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <math.h>

#if MULTI_WORD_INT_SIZEOF_ITEM == 8
typedef __int128      int128_t;
typedef __uint128_t  uint128_t;

typedef uint64_t   mwi_u_word_t;
typedef uint128_t  mwi_u_dword_t;

typedef  int64_t   mwi_s_word_t;
typedef  int128_t  mwi_s_dword_t;

typedef uint64_t   mwi_u_mul_word_t;
# define MWI_MULW_LARGER 0

# define PRIxMWI    PRIx64
# define PRI0xMWI   "016" PRIx64
#elif MULTI_WORD_INT_SIZEOF_ITEM == 4
typedef uint32_t   mwi_u_word_t;
typedef uint64_t   mwi_u_dword_t;

typedef  int32_t   mwi_s_word_t;
typedef  int64_t   mwi_s_dword_t;

# if MULTI_WORD_INT_SIZEOF_MULW == 8
typedef uint64_t   mwi_u_mul_word_t; /* for prime factor up mul */
#  define MWI_MULW_LARGER 1 /* mwi_u_mul_word_t larger than mwi_u_word_t */
# elif MULTI_WORD_INT_SIZEOF_MULW == 4
typedef uint32_t   mwi_u_mul_word_t; /* for prime factor up mul */
#  define MWI_MULW_LARGER 0
# else
#  error MULTI_WORD_INT_SIZEOF_MULW must be 4 or 8.
# endif

# define PRIxMWI    PRIx32
# define PRI0xMWI   "08" PRIx32
#else
# error MULTI_WORD_INT_SIZEOF_ITEM must be 4 or 8.
#endif

#define MWI_SHIFT_BITS (sizeof (mwi_u_word_t) * 8)

#define MWI_MULW_SHIFT_BITS (sizeof (mwi_u_mul_word_t) * 8)

#define MWI_FULL_SIGN_WORD(x)					\
  ((mwi_u_word_t) (((mwi_s_word_t) (x)) >>			\
		   (mwi_s_word_t) (MWI_SHIFT_BITS - 1)))
#define MWI_SIGN_BIT (((mwi_u_word_t) 1) << (MWI_SHIFT_BITS - 1))

/* The sign of the multi-word integer is handled in the usual way: the
 * highest bit of the highest word gives the sign.  This means that a
 * positive result that sets the highest bit has to be extended with
 * an all-zero higher word.
 */

struct multi_word_int
{
  size_t nw;
  size_t nalloc;
  mwi_u_word_t *w;
};

static inline void mwi_realloc(struct multi_word_int *mwi,
			       size_t min_alloc)
{
  if (mwi->nalloc >= min_alloc)
    return;

  while (mwi->nalloc < min_alloc)
    mwi->nalloc *= 2;

  size_t alloc_size = mwi->nalloc * sizeof (mwi_u_word_t);

  mwi->w = (mwi_u_word_t*) realloc (mwi->w, alloc_size);

  if (mwi->w == NULL)
    {
      fprintf (stderr,
	       "wigxjpf: "
	       "Memory allocation error (multi-word int), %zd bytes.\n",
	       alloc_size);
      wigxjpf_error();
    }
}

static inline void mwi_alloc(struct multi_word_int *mwi)
{
  mwi->w = NULL;
  mwi->nalloc = 1; /* Non-0, so *2 gives a value */
  /* 8 is some sort of sane lower default.  6 required by mwi_to_float128. */
  mwi_realloc(mwi, 8); /* Will realloc, as nalloc < 8. */
}

static inline void mwi_free(struct multi_word_int *mwi)
{
  free(mwi->w);
}

static inline size_t mwi_alloc_size(struct multi_word_int *mwi)
{
  return mwi->nalloc * sizeof (mwi_u_word_t);
}

static inline void mwi_dump(const char *name,
			    const struct multi_word_int *mwi)
{
  printf ("%15s:", name);
  ssize_t i;
  for (i = (ssize_t) mwi->nw - 1; i >= 0; i--)
    printf (" 0x%" PRI0xMWI "", mwi->w[i]);
  printf ("\n");
}

static inline void mwi_copy(struct multi_word_int *mwi,
			    const struct multi_word_int *src)
{
  size_t i;
  mwi->nw = src->nw;
  mwi_realloc(mwi, mwi->nw);
  for (i = 0; i < mwi->nw; i++)
    mwi->w[i] = src->w[i];
}

static inline void mwi_set_one_word(struct multi_word_int *mwi,
				    mwi_u_word_t w)
{
  mwi->w[0] = w;
  mwi->nw = 1;
}

static inline void mwi_set_one_mul_word(struct multi_word_int *mwi,
					mwi_u_mul_word_t w)
{
#if !MWI_MULW_LARGER
  mwi_set_one_word(mwi, w);
#else
  mwi->w[0] = (mwi_u_word_t)  w;
  mwi->w[1] = (mwi_u_word_t) (w >> MWI_SHIFT_BITS);
  mwi->nw = (mwi->w[1] == MWI_FULL_SIGN_WORD(mwi->w[0]) ? 1 : 2);
#endif
}

#define MWI_ADD_KERNEL(dest,src1,src2) do {		\
    s = (src1);						\
    t = (src2);						\
    s = s + t + carry;					\
    dest = (mwi_u_word_t) s;				\
    carry = (mwi_u_word_t) (s >> MWI_SHIFT_BITS);	\
  } while ( 0 )

static inline void mwi_add_mwi(struct multi_word_int *mwi,
			       struct multi_word_int *term)
{
  size_t i;
  mwi_u_dword_t t, s;
  mwi_u_word_t carry = 0;

  mwi_realloc(mwi, mwi->nw+1);

  mwi_u_word_t mwi_sign_bits  = MWI_FULL_SIGN_WORD(mwi->w[mwi->nw-1]);
  mwi_u_word_t term_sign_bits = MWI_FULL_SIGN_WORD(term->w[term->nw-1]);

  if (term->nw <= mwi->nw)
    {
      for (i = 0; i < term->nw; i++)
	MWI_ADD_KERNEL(mwi->w[i], mwi->w[i], term->w[i]);
      for ( ; i < mwi->nw; i++)
	MWI_ADD_KERNEL(mwi->w[i], mwi->w[i], term_sign_bits);
    }
  else
    {
      mwi_realloc(mwi, term->nw+1);

      for (i = 0; i < mwi->nw; i++)
	MWI_ADD_KERNEL(mwi->w[i], mwi->w[i], term->w[i]);
      for ( ; i < term->nw; i++)
	MWI_ADD_KERNEL(mwi->w[i], mwi_sign_bits, term->w[i]);

      mwi->nw = term->nw;
    }
  /* We may have to add another word.  Only case when we
   * do not is when the new word has all bits the same as the sign bit,
   * and the sign bit of the lesser word matches.
   */
  {
    mwi_u_word_t next_word =
      mwi_sign_bits + term_sign_bits + carry;
    mwi_u_word_t next_sign_word =
      MWI_FULL_SIGN_WORD(next_word);

    if (next_word != next_sign_word ||
	((next_word ^ mwi->w[mwi->nw-1]) & MWI_SIGN_BIT))
      {
	mwi->w[mwi->nw] = next_word;
	mwi->nw++;
      }
  }
}

#define MWI_SUB_KERNEL(dest,src1,src2) do {		\
    s = (src1);						\
    t = (src2);						\
    s = s - t - carry;					\
    dest = (mwi_u_word_t) s;				\
    carry = ((mwi_u_word_t) (s >> MWI_SHIFT_BITS)) & 1;	\
  } while ( 0 )

static inline void mwi_sub_mwi(struct multi_word_int *mwi,
			       const struct multi_word_int *term)
{
  size_t i;
  mwi_u_dword_t t, s;
  mwi_u_word_t carry = 0;

  mwi_realloc(mwi, mwi->nw+1);

  mwi_u_word_t mwi_sign_bits  = MWI_FULL_SIGN_WORD(mwi->w[mwi->nw-1]);
  mwi_u_word_t term_sign_bits = MWI_FULL_SIGN_WORD(term->w[term->nw-1]);

  if (term->nw <= mwi->nw)
    {
      for (i = 0; i < term->nw; i++)
	MWI_SUB_KERNEL(mwi->w[i], mwi->w[i], term->w[i]);
      for ( ; i < mwi->nw; i++)
	MWI_SUB_KERNEL(mwi->w[i], mwi->w[i], term_sign_bits);
    }
  else
    {
      mwi_realloc(mwi, term->nw+1);

      for (i = 0; i < mwi->nw; i++)
	MWI_SUB_KERNEL(mwi->w[i], mwi->w[i], term->w[i]);
      for ( ; i < term->nw; i++)
	MWI_SUB_KERNEL(mwi->w[i], mwi_sign_bits, term->w[i]);

      mwi->nw = term->nw;
    }
  /* See discussion about extra word in add function above. */
  {
    mwi_u_word_t next_word =
      mwi_sign_bits - term_sign_bits - carry;
    mwi_u_word_t next_sign_word =
      MWI_FULL_SIGN_WORD(next_word);

    if (next_word != next_sign_word ||
	((next_word ^ mwi->w[mwi->nw-1]) & MWI_SIGN_BIT))
      {
	mwi->w[mwi->nw] = next_word;
	mwi->nw++;
      }
  }
}

#define MWI_MUL_KERNEL(dest,src,add_src) do {		\
    mwi_u_word_t s;					\
    s = (src);						\
    p = s * f;						\
    p += from_lower + (add_src);			\
    dest = ((mwi_u_word_t) p);				\
    from_lower = (mwi_u_word_t) (p >> MWI_SHIFT_BITS);	\
  } while ( 0 )

/* Only unsigned factors! */
static inline void mwi_mul_plain(struct multi_word_int *mwi,
				 mwi_u_word_t factor)
{
  size_t i;
  mwi_u_dword_t p, f;
  f = (mwi_u_dword_t) factor;
  mwi_u_word_t from_lower = 0;

  mwi_realloc(mwi, mwi->nw+1);

  for (i = 0; i < mwi->nw; i++)
    MWI_MUL_KERNEL(mwi->w[i],mwi->w[i],0);
  if (from_lower || (mwi->w[mwi->nw-1] & MWI_SIGN_BIT))
    {
      mwi->w[mwi->nw] = from_lower;
      mwi->nw++;
    }
}

static inline void mwi_mul_mwi(struct multi_word_int *dest,
			       const struct multi_word_int *src,
			       const struct multi_word_int *factor)
{
  size_t i, j;

  dest->nw = src->nw + factor->nw;
  mwi_realloc(dest, dest->nw);

  for (j = 0; j < dest->nw; j++)
    dest->w[j] = 0;

  mwi_u_word_t factor_sign_bits = MWI_FULL_SIGN_WORD(factor->w[factor->nw-1]);
  mwi_u_word_t src_sign_bits    = MWI_FULL_SIGN_WORD(src->w[src->nw-1]);

  for (j = 0; j < factor->nw; j++)
    {
      mwi_u_dword_t p, f;
      f = factor->w[j];
      mwi_u_word_t from_lower = 0;
      size_t lim_i = dest->nw - j;
      size_t lim_i2 = src->nw < lim_i ? src->nw : lim_i;
      for (i = 0; i < lim_i2; i++)
	MWI_MUL_KERNEL(dest->w[i+j],src->w[i],
		       (mwi_u_dword_t) dest->w[i+j]);
      if (src_sign_bits)
	for ( ; i < lim_i; i++)
	  MWI_MUL_KERNEL(dest->w[i+j],src_sign_bits,
			 (mwi_u_dword_t) dest->w[i+j]);
      else
	for ( ; from_lower && i < lim_i; i++)
	  MWI_MUL_KERNEL(dest->w[i+j],0,
			 (mwi_u_dword_t) dest->w[i+j]);
    }
  if (factor_sign_bits)
    for ( ; j < dest->nw; j++)
      {
	mwi_u_dword_t p, f;
	f = factor_sign_bits;
	mwi_u_word_t from_lower = 0;
	size_t lim_i = dest->nw - j;
	size_t lim_i2 = src->nw < lim_i ? src->nw : lim_i;
	for (i = 0; i < lim_i2; i++)
	  MWI_MUL_KERNEL(dest->w[i+j],src->w[i],
			 (mwi_u_dword_t) dest->w[i+j]);
	if (src_sign_bits)
	  for ( ; i < lim_i; i++)
	    MWI_MUL_KERNEL(dest->w[i+j],src_sign_bits,
			   (mwi_u_dword_t) dest->w[i+j]);
	else
	  for ( ; from_lower && i < lim_i; i++)
	    MWI_MUL_KERNEL(dest->w[i+j],0,
			   (mwi_u_dword_t) dest->w[i+j]);
      }
  while (dest->nw > 1 &&
	 dest->w[dest->nw-1] == MWI_FULL_SIGN_WORD(dest->w[dest->nw-2]))
    dest->nw--;
}

#define DOUBLE_TYPE         double
#define DBL_FCN_POSTFIX(x)  x##_double
#define DBL_MATH_FCN_LQ(x)  x
#include "multi_word_int_dbl.h"
#undef DBL_MATH_FCN_LQ
#undef DBL_FCN_POSTFIX
#undef DOUBLE_TYPE

#if WIGXJPF_IMPL_LONG_DOUBLE
# define DOUBLE_TYPE         long double
# define DBL_FCN_POSTFIX(x)  x##_long_double
# define DBL_MATH_FCN_LQ(x)  x##l
# include "multi_word_int_dbl.h"
# undef DBL_MATH_FCN_LQ
# undef DBL_FCN_POSTFIX
# undef DOUBLE_TYPE
#endif

#endif/*__WIGXJPF_MULTI_WORD_INT_H__*/
