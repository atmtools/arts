
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

/* This file is included several times from wigxjpf_multi_word_int.h
 * to implement for different floating point types.
 */

static inline void DBL_FCN_POSTFIX(mwi_to)(DOUBLE_TYPE *d, int *exponent,
					   struct multi_word_int *mwi)
{
  size_t high = mwi->nw - 1;
  mwi_u_word_t wi;
  DOUBLE_TYPE ds, di;
  size_t i;

  while (high > 0 &&
	 mwi->w[high] == MWI_FULL_SIGN_WORD(mwi->w[high]) &&
	 !((mwi->w[high] ^ mwi->w[high-1]) & MWI_SIGN_BIT))
    {
      high--;
    }

  /* A integer that is small enough to need no exponent can always be
   * exactly represented in a floating point variable.  This means
   * that the integer parts that fit into the mantissa will come
   * unharmed.  And smaller parts will just discard.  We need not go
   * via large integers to do the conversion from integer to floating
   * point.
   */

  ds = 0;

  for (i = 2; i >= 1; i--)
    {
      wi = (high >= i) ? mwi->w[high-i] : 0;
      di = (DOUBLE_TYPE)  wi;
      di = DBL_MATH_FCN_LQ(ldexp)(di, -((int) i * (int) MWI_SHIFT_BITS));
      ds += di;
    }
  
  wi = mwi->w[high];
  di = (DOUBLE_TYPE) ((mwi_s_word_t) wi);
  ds += di;

  /* TODO: e3, w4 for _float128. */

  *d = ds;
  
  /* Note, exponent is always divisible by two (for sqrt). */
  *exponent = (int) ((high) * MWI_SHIFT_BITS);
}
