
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

#ifndef __WIGXJPF_PRIME_FACTOR_H__
#define __WIGXJPF_PRIME_FACTOR_H__

#include "wigxjpf_config.h"

#include <assert.h>
#include <stdint.h>
#include "multi_word_int.h"

#if PRIME_LIST_SIZEOF_ITEM == 2
typedef int16_t       prime_exp_t;
typedef uint16_t    u_prime_exp_t;
#elif PRIME_LIST_SIZEOF_ITEM == 4
typedef int32_t       prime_exp_t;
typedef uint32_t    u_prime_exp_t;
#else
# error PRIME_LIST_SIZEOF_ITEM must be 2 or 4.
#endif

#if PRIME_LIST_USE_VECTOR
typedef prime_exp_t v_prime_exp_t 
__attribute__ ((vector_size (PRIME_LIST_VECT_SIZE)));
# define V_PRIME_EXP_T  v_prime_exp_t
# define V_EXPO         v_expo
# define PRIME_LIST_VECT_BLOCK_ITEMS \
  (PRIME_LIST_VECT_SIZE / PRIME_LIST_SIZEOF_ITEM)
# if PRIME_LIST_VECT_BLOCK_ITEMS == 2
#  define V_ITEM_INIT(x)   { x, x }
# elif PRIME_LIST_VECT_BLOCK_ITEMS == 4
#  define V_ITEM_INIT(x)   { x, x, x, x }
# elif PRIME_LIST_VECT_BLOCK_ITEMS == 8
#  define V_ITEM_INIT(x)   { x, x, x, x, x, x, x, x }
# elif PRIME_LIST_VECT_BLOCK_ITEMS == 16
#  define V_ITEM_INIT(x)   { x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x }
# endif
#else
# define V_PRIME_EXP_T    prime_exp_t
# define V_EXPO           expo
# define PRIME_LIST_VECT_BLOCK_ITEMS 1
# define V_ITEM_INIT(x)   x
#endif

struct prime_exponents
{
  int num_blocks; /* how many items (in units of vector blocks). */
  union
  {
    /* Note: the size of the arrays below are not included in memory
     * allocations.  They are dynamic and use the offset of expo in
     * the structure.
     */
    prime_exp_t      expo[PRIME_LIST_VECT_BLOCK_ITEMS * 4];
#if PRIME_LIST_USE_VECTOR
    v_prime_exp_t  v_expo[4];
#endif
  };
};

extern uint32_t *wigxjpf_prime_list;

extern size_t wigxjpf_prime_fact_stride;
extern int    wigxjpf_max_prime_decomp;
extern void *wigxjpf_prime_factors_base;
extern void *wigxjpf_prime_factors_1;
extern void *wigxjpf_prime_factors_2;

#define PRIME_FACTOR(x)							\
  ((struct prime_exponents*) (((char *) wigxjpf_prime_factors_1) +	\
			      (size_t) (x) * wigxjpf_prime_fact_stride))

#define FACTORIAL_PRIME_FACTOR(x)					\
  ((struct prime_exponents*) (((char *) wigxjpf_prime_factors_2) +	\
			      (size_t) (x) * wigxjpf_prime_fact_stride))

#define PRIME_FACTOR_UPDOWN(x, type, cast_type, direction_op) do {	\
    x = ((type *) (((cast_type *) x) direction_op			\
		   wigxjpf_prime_fact_stride));				\
  } while ( 0 )

#define PRIME_FACTOR_UP(x) \
  PRIME_FACTOR_UPDOWN(x, struct prime_exponents, char, +)
#define CONST_PRIME_FACTOR_UP(x) \
  PRIME_FACTOR_UPDOWN(x, const struct prime_exponents, const char, +)
#define CONST_PRIME_FACTOR_DOWN(x) \
  PRIME_FACTOR_UPDOWN(x, const struct prime_exponents, const char, -)

static inline void pexpo_set_zero(struct prime_exponents *fpf,
				  int num_blocks)
{
  int i;

  fpf->num_blocks = num_blocks;

  V_PRIME_EXP_T *p_dest = fpf->V_EXPO;
  V_PRIME_EXP_T  init = V_ITEM_INIT(0);

  for (i = 0; i < num_blocks; i++)
    *(p_dest++) = init;
}

static inline void pexpo_expand_blocks(struct prime_exponents *fpf,
				       int num_blocks)
{
  int i;

  if (fpf->num_blocks >= num_blocks)
    return;

  V_PRIME_EXP_T *p_dest = fpf->V_EXPO;
  V_PRIME_EXP_T  init = V_ITEM_INIT(0);

  for (i = fpf->num_blocks; i < num_blocks; i++)
    p_dest[i] = init;

  fpf->num_blocks = num_blocks;
}

/* To be able to handle also negative exponents with bit-tricks in
 * comparisons (for vector types), we can only use half the range,
 * i.e. 0x3fff.
 */

static inline void pexpo_set_max(struct prime_exponents *fpf,
				 int num_blocks)
{
  int i;
  const prime_exp_t max_exp = (prime_exp_t) (((u_prime_exp_t) -1) >> 2);

  fpf->num_blocks = num_blocks;

  V_PRIME_EXP_T *p_dest = fpf->V_EXPO;
  V_PRIME_EXP_T  init = V_ITEM_INIT(max_exp);

  for (i = 0; i < num_blocks; i++)
    *(p_dest++) = init;
}

static inline void pexpo_keep_min(struct prime_exponents *keep_fpf,
				  const struct prime_exponents *in_fpf)
{
  int i;

  assert(keep_fpf->num_blocks == in_fpf->num_blocks);

#if !PRIME_LIST_USE_VECTOR
  for (i = 0; i < keep_fpf->num_blocks; i++)
    {
      keep_fpf->expo[i] =
	(in_fpf->expo[i] < keep_fpf->expo[i]) ? 
	in_fpf->expo[i] : keep_fpf->expo[i];
    }
#else
  for (i = 0; i < keep_fpf->num_blocks; i++)
    {
      v_prime_exp_t in_smaller_sign_bit =
	in_fpf->v_expo[i] - keep_fpf->v_expo[i];

      v_prime_exp_t shift_sign =
	V_ITEM_INIT(sizeof (prime_exp_t) * 8 - 1);

      v_prime_exp_t in_smaller_all_bits   =
	in_smaller_sign_bit >> shift_sign;
      v_prime_exp_t keep_smaller_all_bits = ~in_smaller_all_bits;

      v_prime_exp_t smallest =
	(in_fpf->v_expo[i] & in_smaller_all_bits) |
	(keep_fpf->v_expo[i] & keep_smaller_all_bits);

      keep_fpf->v_expo[i] = smallest;
    }
#endif
}

static inline void pexpo_keep_min_in_as_diff(struct prime_exponents *keep_fpf,
					     struct prime_exponents *in_fpf)
{
  int i;

  assert(keep_fpf->num_blocks == in_fpf->num_blocks);

#if !PRIME_LIST_USE_VECTOR
  for (i = 0; i < keep_fpf->num_blocks; i++)
    {
      prime_exp_t tmp = in_fpf->expo[i] - keep_fpf->expo[i];

      keep_fpf->expo[i] =
	(in_fpf->expo[i] < keep_fpf->expo[i]) ? 
	in_fpf->expo[i] : keep_fpf->expo[i];

      in_fpf->expo[i] = tmp;
    }
#else
  for (i = 0; i < keep_fpf->num_blocks; i++)
    {
      v_prime_exp_t tmp = in_fpf->v_expo[i] - keep_fpf->v_expo[i];

      v_prime_exp_t in_smaller_sign_bit =
	in_fpf->v_expo[i] - keep_fpf->v_expo[i];

      v_prime_exp_t shift_sign =
	V_ITEM_INIT(sizeof (prime_exp_t) * 8 - 1);

      v_prime_exp_t in_smaller_all_bits   =
	in_smaller_sign_bit >> shift_sign;
      v_prime_exp_t keep_smaller_all_bits = ~in_smaller_all_bits;

      v_prime_exp_t smallest =
	(in_fpf->v_expo[i] & in_smaller_all_bits) |
	(keep_fpf->v_expo[i] & keep_smaller_all_bits);

      keep_fpf->v_expo[i] = smallest;

      in_fpf->v_expo[i] = tmp;
    }
#endif
}

struct pexpo_eval_temp
{
  struct multi_word_int prod_pos[2];
  struct multi_word_int prod_neg[2];
  struct multi_word_int factor[2];
  struct multi_word_int big_up[2];
};

static inline void pexpo_prime_factor(struct multi_word_int factor[2],
				      int *factor_active,
				      int64_t prime, prime_exp_t fpf,
				      struct pexpo_eval_temp *temp)
{
  mwi_u_mul_word_t fact = 1;
  mwi_u_mul_word_t up = (mwi_u_mul_word_t) prime;

  for ( ; ; )
    {
      mwi_u_mul_word_t mult;

      mult = (fpf & 1) ? up : 1;

      fact *= mult;

      up = up * up;

      fpf >>= 1;
	  
      if (!fpf)
	break;

      /* We must not hit the sign bit. */
      if (up & (((mwi_u_mul_word_t) -1) << (MWI_MULW_SHIFT_BITS/2-1)))
	{
	  goto full_mult;
	}
    }

  mwi_set_one_mul_word(&factor[0], fact);
  *factor_active = 0;
  return;
  
 full_mult:
  {
    int up_active = 0;
    int fact_active = 0;
 
    mwi_set_one_mul_word(&temp->big_up[up_active], up);
    mwi_set_one_mul_word(&factor[fact_active], fact);
  
    for ( ; ; )
      {
	if (fpf & 1)
	  {
	    mwi_mul_mwi(&factor[!fact_active],
			&factor[fact_active], &temp->big_up[up_active]);
	    fact_active = !fact_active;
	  }

	fpf >>= 1;
	if (!fpf)
	  break;

	mwi_mul_mwi(&temp->big_up[!up_active],
		    &temp->big_up[up_active], &temp->big_up[up_active]);
	up_active = !up_active;
      }
    *factor_active = fact_active;
  }
}

static inline void pexpo_evaluate(struct multi_word_int *big_prod,
				  const struct prime_exponents *in_fpf,
				  struct pexpo_eval_temp *temp)
{
  int i;
  int active = 0;
  mwi_set_one_word(&temp->prod_pos[active], 1);

  for (i = 0; i < in_fpf->num_blocks * PRIME_LIST_VECT_BLOCK_ITEMS; i++)
    {
      prime_exp_t fpf = in_fpf->expo[i];

      if (!fpf)
	continue;

      int64_t prime = wigxjpf_prime_list[i];
      int factor_active;

      pexpo_prime_factor(temp->factor, &factor_active, prime, fpf, temp);

      if (temp->factor[factor_active].nw == 1)
	mwi_mul_plain(&temp->prod_pos[active],
		      temp->factor[factor_active].w[0]);
      else
	{
	  mwi_mul_mwi(&temp->prod_pos[!active], &temp->prod_pos[active],
		      &temp->factor[factor_active]);	
	  active = !active;
	}
    }

  mwi_copy(big_prod,&temp->prod_pos[active]);
}

static inline void pexpo_evaluate2(struct multi_word_int *big_prod_pos,
				   struct multi_word_int *big_prod_neg,
				   const struct prime_exponents *in_fpf,
				   struct pexpo_eval_temp *temp)
{
  int i;
  int active_pos = 0, active_neg = 0;
  mwi_set_one_word(&temp->prod_pos[active_pos], 1);
  mwi_set_one_word(&temp->prod_neg[active_neg], 1);

  for (i = 0; i < in_fpf->num_blocks * PRIME_LIST_VECT_BLOCK_ITEMS; i++)
    {
      prime_exp_t fpf = in_fpf->expo[i];
      
      if (!fpf)
	continue;
      
      int64_t prime = wigxjpf_prime_list[i];

      if (fpf > 0)
	{
	  int factor_active;

	  pexpo_prime_factor(temp->factor, &factor_active, prime, fpf, temp);

	  if (temp->factor[factor_active].nw == 1)
	    mwi_mul_plain(&temp->prod_pos[active_pos],
			  temp->factor[factor_active].w[0]);
	  else
	    {
	      mwi_mul_mwi(&temp->prod_pos[!active_pos],
			  &temp->prod_pos[active_pos],
			  &temp->factor[factor_active]);	
	      active_pos = !active_pos;
	    }
	}
      else
	{
	  fpf = -fpf;

	  int factor_active;

	  pexpo_prime_factor(temp->factor, &factor_active, prime, fpf, temp);

	  if (temp->factor[factor_active].nw == 1)
	    mwi_mul_plain(&temp->prod_neg[active_neg],
			  temp->factor[factor_active].w[0]);
	  else
	    {
	      mwi_mul_mwi(&temp->prod_neg[!active_neg],
			  &temp->prod_neg[active_neg],
			  &temp->factor[factor_active]);	
	      active_neg = !active_neg;
	    }
	}
    }

  mwi_copy(big_prod_pos,&temp->prod_pos[active_pos]);
  mwi_copy(big_prod_neg,&temp->prod_neg[active_neg]);
}

static inline void pexpo_copy(struct prime_exponents *dest_fpf,
			      const struct prime_exponents *a_fpf)
{
  dest_fpf->num_blocks = a_fpf->num_blocks;

  V_PRIME_EXP_T *p_dest = dest_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_a = a_fpf->V_EXPO;

  int i;

  for (i = 0; i < dest_fpf->num_blocks; i++)
    *(p_dest++) = *(p_a++);
}

static inline void pexpo_expand_add(struct prime_exponents *dest_fpf,
				    const struct prime_exponents *a_fpf)
{
  pexpo_expand_blocks(dest_fpf, a_fpf->num_blocks);

  V_PRIME_EXP_T *p_dest = dest_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_a = a_fpf->V_EXPO;

  int i;

  for (i = 0; i < dest_fpf->num_blocks; i++)
    *(p_dest++) += *(p_a++);
}

static inline void pexpo_expand_sub(struct prime_exponents *dest_fpf,
				    const struct prime_exponents *a_fpf)
{
  pexpo_expand_blocks(dest_fpf, a_fpf->num_blocks);

  V_PRIME_EXP_T *p_dest = dest_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_a = a_fpf->V_EXPO;

  int i;

  for (i = 0; i < dest_fpf->num_blocks; i++)
    *(p_dest++) -= *(p_a++);
}

static inline void pexpo_expand_sum3(struct prime_exponents *dest_fpf,
				     struct prime_exponents *a_fpf,
				     struct prime_exponents *b_fpf,
				     struct prime_exponents *c_fpf)
{
  dest_fpf->num_blocks = a_fpf->num_blocks;
  if (b_fpf->num_blocks > dest_fpf->num_blocks)
    dest_fpf->num_blocks = b_fpf->num_blocks;
  if (c_fpf->num_blocks > dest_fpf->num_blocks)
    dest_fpf->num_blocks = c_fpf->num_blocks;

  pexpo_expand_blocks(a_fpf, dest_fpf->num_blocks);
  pexpo_expand_blocks(b_fpf, dest_fpf->num_blocks);
  pexpo_expand_blocks(c_fpf, dest_fpf->num_blocks);

  V_PRIME_EXP_T *p_dest = dest_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_a = a_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_b = b_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_c = c_fpf->V_EXPO;

  int i;

  for (i = 0; i < dest_fpf->num_blocks; i++)
    *(p_dest++) = *(p_a++) + *(p_b++) + *(p_c++);
}

static inline void pexpo_add3_sub(struct prime_exponents *dest_fpf,
				  const struct prime_exponents *a_fpf,
				  const struct prime_exponents *b_fpf,
				  const struct prime_exponents *c_fpf,
				  const struct prime_exponents *d_fpf)
{
  V_PRIME_EXP_T *p_dest = dest_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_a = a_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_b = b_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_c = c_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_d = d_fpf->V_EXPO;

  int i;

  for (i = 0; i < dest_fpf->num_blocks; i++)
    *(p_dest++) += *(p_a++) + *(p_b++) + *(p_c++) - *(p_d++);
}

static inline void pexpo_add6(struct prime_exponents *dest_fpf,
			      const struct prime_exponents *a_fpf,
			      const struct prime_exponents *b_fpf,
			      const struct prime_exponents *c_fpf,
			      const struct prime_exponents *d_fpf,
			      const struct prime_exponents *e_fpf,
			      const struct prime_exponents *f_fpf)
{
  V_PRIME_EXP_T *p_dest = dest_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_a = a_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_b = b_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_c = c_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_d = d_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_e = e_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_f = f_fpf->V_EXPO;

  int i;

  for (i = 0; i < dest_fpf->num_blocks; i++)
    *(p_dest++) += 
      *(p_a++) + *(p_b++) + *(p_c++) + *(p_d++) + *(p_e++) + *(p_f++);
}

static inline void pexpo_add_sub3(struct prime_exponents *dest_fpf,
				  const struct prime_exponents *a_fpf,
				  const struct prime_exponents *b_fpf,
				  const struct prime_exponents *c_fpf,
				  const struct prime_exponents *d_fpf)
{
  V_PRIME_EXP_T *p_dest = dest_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_a = a_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_b = b_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_c = c_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_d = d_fpf->V_EXPO;

  int i;

  for (i = 0; i < dest_fpf->num_blocks; i++)
    *(p_dest++) += *(p_a++) - *(p_b++) - *(p_c++) - *(p_d++);
}

static inline void pexpo_sum_sub7(struct prime_exponents *dest_fpf,
				  const struct prime_exponents *a_fpf,
				  const struct prime_exponents *b_fpf,
				  const struct prime_exponents *c_fpf,
				  const struct prime_exponents *d_fpf,
				  const struct prime_exponents *e_fpf,
				  const struct prime_exponents *f_fpf,
				  const struct prime_exponents *g_fpf,
				  const struct prime_exponents *h_fpf,
				  int num_blocks)
{
  dest_fpf->num_blocks = num_blocks;

  V_PRIME_EXP_T *p_dest = dest_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_a = a_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_b = b_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_c = c_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_d = d_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_e = e_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_f = f_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_g = g_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_h = h_fpf->V_EXPO;

  int i;

  for (i = 0; i < dest_fpf->num_blocks; i++)
    *(p_dest++) =
      *(p_a++) -
      *(p_b++) - *(p_c++) - *(p_d++) - *(p_e++) - 
      *(p_f++) - *(p_g++) - *(p_h++);
}

static inline void pexpo_sum0_sub6(struct prime_exponents *dest_fpf,
				   const struct prime_exponents *a_fpf,
				   const struct prime_exponents *b_fpf,
				   const struct prime_exponents *c_fpf,
				   const struct prime_exponents *d_fpf,
				   const struct prime_exponents *e_fpf,
				   const struct prime_exponents *f_fpf,
				   int num_blocks)
{
  dest_fpf->num_blocks = num_blocks;

  V_PRIME_EXP_T *p_dest = dest_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_a = a_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_b = b_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_c = c_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_d = d_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_e = e_fpf->V_EXPO;
  const V_PRIME_EXP_T *p_f = f_fpf->V_EXPO;

  int i;

  for (i = 0; i < dest_fpf->num_blocks; i++)
    *(p_dest++) = 
      - *(p_a++) - *(p_b++) - *(p_c++) - *(p_d++) - *(p_e++) - *(p_f++);
}

static inline void pexpo_split_sqrt_add(struct prime_exponents *src_dest_fpf,
					struct multi_word_int *big_sqrt,
					struct prime_exponents *add_fpf)
{
  int i;
  mwi_set_one_word(big_sqrt, 1);

  int num_blocks = src_dest_fpf->num_blocks > add_fpf->num_blocks ?
    src_dest_fpf->num_blocks : add_fpf->num_blocks;

  pexpo_expand_blocks(src_dest_fpf,num_blocks);
  pexpo_expand_blocks(add_fpf,num_blocks);

  for (i = 0; i < src_dest_fpf->num_blocks * PRIME_LIST_VECT_BLOCK_ITEMS; i++)
    {
      prime_exp_t sqrt_fpf = src_dest_fpf->expo[i] & 1;
      src_dest_fpf->expo[i] += sqrt_fpf;
      src_dest_fpf->expo[i] /= 2;

      src_dest_fpf->expo[i] = src_dest_fpf->expo[i] + add_fpf->expo[i];

      if (!sqrt_fpf)
	continue;

      mwi_mul_plain(big_sqrt, wigxjpf_prime_list[i]);
    }
}

#endif/*__WIGXJPF_PRIME_FACTOR_H__*/
