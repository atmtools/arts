
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

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <inttypes.h>

#include "wigxjpf.h"

#include "wigxjpf_config.h"
#include "calc_minmax.h"

#include "calc.h"

#if ACCOUNT_MAX_FACT_ITER
int totmaxfact = 0;
int totmaxiter = 0;

# define ACCOUNT_MAX_FACT(maxfact) do {				\
    if ((maxfact) > totmaxfact) { totmaxfact = (maxfact); }	\
  } while ( 0 )
# define ACCOUNT_MAX_ITER(iter) do {			        \
    if ((iter) > totmaxiter) { totmaxiter = (iter); }	        \
  } while ( 0 )

void print_max_factiter(int lim)
{
  printf ("MAXITER %3d %4d %4d\n", lim, totmaxfact, totmaxiter);
}
#else
# define ACCOUNT_MAX_FACT(maxfact) do { } while ( 0 )
# define ACCOUNT_MAX_ITER(iter)    do { } while ( 0 )
#endif

#define CHECK_MAX_PRIMEFACT_FACTORIAL(maxfact) do {		\
    if ((maxfact) > wigxjpf_max_prime_decomp) {			\
      if (wigxjpf_max_prime_decomp == -1)			\
	fprintf (stderr, "wigxjpf: "				\
		 "No factorials table.\n");			\
      else							\
	fprintf (stderr, "wigxjpf: "				\
		 "Too large factorial (%d!).  "			\
		 "Increase MAX_FACT.\n", (maxfact));		\
      wigxjpf_error();						\
    }								\
    ACCOUNT_MAX_FACT(maxfact);					\
  } while (0)

#define CHECK_TEMP do {						\
    if (csi == NULL) {						\
      fprintf (stderr, "wigxjpf: "				\
	       "Temp array not allocated.  "			\
	       "(For this thread?)\n");				\
      wigxjpf_error();						\
    }							        \
    if (csi->inuse) {						\
      fprintf (stderr, "wigxjpf: "				\
	       "Temp array already in use.  "			\
	       "Running multi-threaded, "			\
	       "but wigxjpf not compiled for that?\n");		\
      wigxjpf_error();						\
    }								\
    csi->inuse = 1;						\
  } while (0)

#define CHECK_MAX_ITER(iter) do {			        \
    if ((iter) > csi->max_iter) {			        \
      fprintf (stderr, "wigxjpf: "				\
	       "More iterations (%d) than allocated "		\
	       "in temp array.\n", (iter));			\
      wigxjpf_error();						\
    }							        \
    ACCOUNT_MAX_ITER(iter);				        \
  } while (0)

void WIGXJPF_NOINLINE delta_coeff(int two_a, int two_b, int two_c,
				  struct prime_exponents *prefact_fpf)
{
  /* Check maximum factorial. */

  int max_factorial = (two_a + two_b + two_c) / 2 + 1;

  CHECK_MAX_PRIMEFACT_FACTORIAL(max_factorial);

  /* Calculate. */

  const struct prime_exponents *p_n1 =
    FACTORIAL_PRIME_FACTOR((two_a + two_b - two_c) / 2);
  const struct prime_exponents *p_n2 =
    FACTORIAL_PRIME_FACTOR((two_a + two_c - two_b) / 2);
  const struct prime_exponents *p_n3 =
    FACTORIAL_PRIME_FACTOR((two_b + two_c - two_a) / 2);

  const struct prime_exponents *p_d1 =
    FACTORIAL_PRIME_FACTOR((two_a + two_b + two_c) / 2 + 1);

  /* We may have to expand the blocks in prefact_fpf to handle the
   * maximum new factorials.
   */

  pexpo_expand_blocks(prefact_fpf, p_d1->num_blocks);

  FPF_DUMP_FACT("n1_fpf", p_n1);
  FPF_DUMP_FACT("n2_fpf", p_n2);
  FPF_DUMP_FACT("n3_fpf", p_n3);
  FPF_DUMP_FACT("d1_fpf", p_d1);

  pexpo_add3_sub(prefact_fpf, p_n1, p_n2, p_n3, p_d1);
}

void calcsum_3j(struct wigxjpf_temp *csi,
		int two_j1, int two_j2, int two_j3,
		int two_m1, int two_m2, int two_m3)
{
  int kmin = max3(two_j1 + two_m2 - two_j3,
		  two_j2 - two_m1 - two_j3,
		  0) / 2;
  int kmax = min3(two_j2 + two_m2,
		  two_j1 - two_m1,
		  two_j1 + two_j2 - two_j3) / 2;

  CHECK_TEMP;

#if DEBUG_PRINT
  printf ("[3j] krange: [%d,%d] steps %d\n",
	  kmin, kmax,
	  kmax - kmin + 1);
#endif

  /* What is the maximum factorial? */

  int max_factorial = (two_j1 + two_j2 + two_j3) / 2 + 1;

  CHECK_MAX_PRIMEFACT_FACTORIAL(max_factorial);

  /* Limit: max fact = 3 * j + 1. */

  pexpo_set_max(CSI_MIN_NUME_FPF(csi),
		FACTORIAL_PRIME_FACTOR(max_factorial)->num_blocks);

  /* Then we do a first round over the loop.  Find out what factors
   * the numerator - denominator actually has, and what we thus have
   * to multiply the numerator with.  The important result is: what
   * is common to all numerators in the sum?
   */

  int k;

  int k_lim = kmax - kmin;

  CHECK_MAX_ITER(k_lim + 1);

  /* Limit: max iter = j + 1. */

  const struct prime_exponents *p_d1 =
    FACTORIAL_PRIME_FACTOR(kmin);

  const struct prime_exponents *p_d2 =
    FACTORIAL_PRIME_FACTOR(kmin + (two_j3 - two_j1 - two_m2) / 2);
  const struct prime_exponents *p_d3 =
    FACTORIAL_PRIME_FACTOR(kmin + (two_j3 - two_j2 + two_m1) / 2);

  const struct prime_exponents *p_d4 =
    FACTORIAL_PRIME_FACTOR((two_j2 + two_m2) / 2 - kmin);
  const struct prime_exponents *p_d5 =
    FACTORIAL_PRIME_FACTOR((two_j1 - two_m1) / 2 - kmin);
  const struct prime_exponents *p_d6 =
    FACTORIAL_PRIME_FACTOR((two_j1 + two_j2 - two_j3) / 2 - kmin);

  struct prime_exponents *nume_fpf;

  nume_fpf = CSI_K_ITER_FPF(csi);

  for (k = 0; k <= k_lim; k++)
    {
      pexpo_sum0_sub6(nume_fpf,
		      p_d1,
		      p_d2, p_d3, p_d4, p_d5, p_d6,
		      CSI_MIN_NUME_FPF(csi)->num_blocks);

      FPF_DUMP("nume_fpf", nume_fpf);

      pexpo_keep_min(CSI_MIN_NUME_FPF(csi), nume_fpf);

      PRIME_FACTOR_UP(nume_fpf);

      CONST_PRIME_FACTOR_UP(p_d1);
      CONST_PRIME_FACTOR_UP(p_d2);
      CONST_PRIME_FACTOR_UP(p_d3);
      CONST_PRIME_FACTOR_DOWN(p_d4);
      CONST_PRIME_FACTOR_DOWN(p_d5);
      CONST_PRIME_FACTOR_DOWN(p_d6);
    }

  FPF_DUMP("min_nume_fpf", CSI_MIN_NUME_FPF(csi));

  /* And then we do the second loop, where we multiply up and sum the
   * numerators.
   */

  mwi_set_one_word(&csi->sum_prod, 0);

  int sign = kmin ^ ((two_j1 - two_j2 - two_m3) / 2);

  nume_fpf = CSI_K_ITER_FPF(csi);

  for (k = 0; k <= k_lim; k++)
    {
      pexpo_expand_sub(nume_fpf, CSI_MIN_NUME_FPF(csi));

      FPF_DUMP("nume_fpf", nume_fpf);

      pexpo_evaluate(&csi->big_prod, nume_fpf, &csi->pexpo_eval);

      MWI_DUMP("big_prod", &csi->big_prod);

      if ((k ^ sign) & 1)
	mwi_sub_mwi(&csi->sum_prod, &csi->big_prod);
      else
	mwi_add_mwi(&csi->sum_prod, &csi->big_prod);

      PRIME_FACTOR_UP(nume_fpf);
   }

  MWI_DUMP("sum_prod", &csi->sum_prod);

  pexpo_set_zero(CSI_PREFACT_FPF(csi),
		 CSI_MIN_NUME_FPF(csi)->num_blocks);

  {
    const struct prime_exponents *p_n1 =
      FACTORIAL_PRIME_FACTOR((  two_j1 + two_j2 - two_j3) / 2);
    const struct prime_exponents *p_n2 =
      FACTORIAL_PRIME_FACTOR((  two_j1 - two_j2 + two_j3) / 2);
    const struct prime_exponents *p_n3 =
      FACTORIAL_PRIME_FACTOR((- two_j1 + two_j2 + two_j3) / 2);

    const struct prime_exponents *p_d1 =
      FACTORIAL_PRIME_FACTOR((two_j1 + two_j2 + two_j3) / 2 + 1);

    pexpo_add3_sub(CSI_PREFACT_FPF(csi), p_n1, p_n2, p_n3, p_d1);

    const struct prime_exponents *p_n4 =
      FACTORIAL_PRIME_FACTOR((two_j1 - two_m1) / 2);
    const struct prime_exponents *p_n5 =
      FACTORIAL_PRIME_FACTOR((two_j1 + two_m1) / 2);
    const struct prime_exponents *p_n6 =
      FACTORIAL_PRIME_FACTOR((two_j2 - two_m2) / 2);
    const struct prime_exponents *p_n7 =
      FACTORIAL_PRIME_FACTOR((two_j2 + two_m2) / 2);
    const struct prime_exponents *p_n8 =
      FACTORIAL_PRIME_FACTOR((two_j3 - two_m3) / 2);
    const struct prime_exponents *p_n9 =
      FACTORIAL_PRIME_FACTOR((two_j3 + two_m3) / 2);

    /* We do it in two steps (3-1 and 6) instead of 9-1 directly,
     * as the compiler otherwise is likely to run out of registers...
     */

    pexpo_add6(CSI_PREFACT_FPF(csi), p_n4, p_n5, p_n6, p_n7, p_n8, p_n9);
  }

  FPF_DUMP("prefact_fpf", CSI_PREFACT_FPF(csi));
}

void factor_6j(struct wigxjpf_temp *csi,
	       int two_j1, int two_j2, int two_j3,
	       int two_j4, int two_j5, int two_j6,
	       struct prime_exponents *min_nume_fpf,
	       struct multi_word_int *sum_prod)
{
  int two_a = two_j1;
  int two_b = two_j2;
  int two_c = two_j5;
  int two_d = two_j4;
  int two_e = two_j3;
  int two_f = two_j6;

  int kmin = max4(two_a + two_b + two_e,
		  two_c + two_d + two_e,
		  two_a + two_c + two_f,
		  two_b + two_d + two_f) / 2;
  int kmax = min3(two_a + two_b + two_c + two_d,
		  two_a + two_d + two_e + two_f,
		  two_b + two_c + two_e + two_f) / 2;

  /* Trying to use the Regge symmetries to find symbols that have
   * fewer terms in the sum leads nowhere: always same number of
   * iterations.
   */

#if DEBUG_PRINT
  printf ("[6j] krange: [%d,%d] steps %d\n",
	  kmin, kmax,
	  kmax - kmin + 1);
#endif

  /* What is the maximum factorial? */

  int max_factorial = max4(kmax + 1,
			   (two_a + two_b + two_c + two_d) / 2,
			   (two_a + two_d + two_e + two_f) / 2,
			   (two_b + two_c + two_e + two_f) / 2);

  CHECK_MAX_PRIMEFACT_FACTORIAL(max_factorial);

  /* Limit 6j: max fact = 4 * j + 1. */
  /* Limit 9j: max fact = 5 * j + 1. */

  pexpo_set_max(min_nume_fpf,
		FACTORIAL_PRIME_FACTOR(max_factorial)->num_blocks);

  /* Then we do a first round over the loop.  Find out what factors
   * the numerator - denominator actually has, and what we thus have
   * to multiply the numerator with.  The important result is: what
   * is common to all numerators in the sum?
   */

  int k;

  int k_lim = kmax - kmin;

  CHECK_MAX_ITER(k_lim + 1);

  /* Limit: max iter = ?. */

  const struct prime_exponents *p_d1 =
    FACTORIAL_PRIME_FACTOR(kmin - (two_a + two_b + two_e) / 2);
  const struct prime_exponents *p_d2 =
    FACTORIAL_PRIME_FACTOR(kmin - (two_c + two_d + two_e) / 2);
  const struct prime_exponents *p_d3 =
    FACTORIAL_PRIME_FACTOR(kmin - (two_a + two_c + two_f) / 2);
  const struct prime_exponents *p_d4 =
    FACTORIAL_PRIME_FACTOR(kmin - (two_b + two_d + two_f) / 2);

  const struct prime_exponents *p_d5 =
    FACTORIAL_PRIME_FACTOR((two_a + two_b + two_c + two_d) / 2 - kmin);
  const struct prime_exponents *p_d6 =
    FACTORIAL_PRIME_FACTOR((two_a + two_d + two_e + two_f) / 2 - kmin);
  const struct prime_exponents *p_d7 =
    FACTORIAL_PRIME_FACTOR((two_b + two_c + two_e + two_f) / 2 - kmin);

  const struct prime_exponents *p_n1 =
    FACTORIAL_PRIME_FACTOR(kmin + 1);

  struct prime_exponents *nume_fpf;

  nume_fpf = CSI_K_ITER_FPF(csi);

  for (k = 0; k <= k_lim; k++)
    {
      FPF_DUMP_FACT("n1_fpf", p_n1);
      FPF_DUMP_FACT("d1_fpf", p_d1);
      FPF_DUMP_FACT("d2_fpf", p_d2);
      FPF_DUMP_FACT("d3_fpf", p_d3);
      FPF_DUMP_FACT("d4_fpf", p_d4);
      FPF_DUMP_FACT("d5_fpf", p_d5);
      FPF_DUMP_FACT("d6_fpf", p_d6);
      FPF_DUMP_FACT("d7_fpf", p_d7);

      pexpo_sum_sub7(nume_fpf,
		     p_n1,
		     p_d1, p_d2, p_d3, p_d4, p_d5, p_d6, p_d7,
		     min_nume_fpf->num_blocks);

      FPF_DUMP("nume_fpf", nume_fpf);

      pexpo_keep_min(min_nume_fpf, nume_fpf);

      PRIME_FACTOR_UP(nume_fpf);

      CONST_PRIME_FACTOR_UP(p_n1);
      CONST_PRIME_FACTOR_UP(p_d1);
      CONST_PRIME_FACTOR_UP(p_d2);
      CONST_PRIME_FACTOR_UP(p_d3);
      CONST_PRIME_FACTOR_UP(p_d4);
      CONST_PRIME_FACTOR_DOWN(p_d5);
      CONST_PRIME_FACTOR_DOWN(p_d6);
      CONST_PRIME_FACTOR_DOWN(p_d7);
    }

  FPF_DUMP("(LCD) min_nume_fpf", min_nume_fpf);

  /* And then we do the second loop, where we multiply up and sum the
   * numerators.
   */

  mwi_set_one_word(sum_prod, 0);

  nume_fpf = CSI_K_ITER_FPF(csi);

  for (k = 0; k <= k_lim; k++)
    {
      FPF_DUMP("term_nume_fpf", nume_fpf);

      pexpo_expand_sub(nume_fpf, min_nume_fpf);

      FPF_DUMP("term-lcd_nume_fpf", nume_fpf);

      pexpo_evaluate(&csi->big_prod, nume_fpf, &csi->pexpo_eval);

      MWI_DUMP("big_prod", &csi->big_prod);

#if DEBUG_PRINT
      printf ("%15s: %d\n","sign",(k ^ kmin) & 1 ? -1 : 1);
#endif

      if ((k ^ kmin) & 1)
	mwi_sub_mwi(sum_prod, &csi->big_prod);
      else
	mwi_add_mwi(sum_prod, &csi->big_prod);

      PRIME_FACTOR_UP(nume_fpf);
   }

  MWI_DUMP("sum_prod", sum_prod);
}

void calcsum_6j(struct wigxjpf_temp *csi,
	       int two_j1, int two_j2, int two_j3,
	       int two_j4, int two_j5, int two_j6)
{
  int two_a = two_j1;
  int two_b = two_j2;
  int two_c = two_j5;
  int two_d = two_j4;
  int two_e = two_j3;
  int two_f = two_j6;

#if DEBUG_PRINT
  printf ("calc 6j : jj : %d %d %d  %d %d %d\n",
	  two_j1, two_j2, two_j3,
	  two_j4, two_j5, two_j6);
#endif

  CHECK_TEMP;

  factor_6j(csi,
	    two_j1, two_j2, two_j3,
	    two_j4, two_j5, two_j6,
	    CSI_MIN_NUME_FPF(csi), &csi->sum_prod);

  pexpo_set_zero(CSI_PREFACT_FPF(csi), 0);

  delta_coeff(two_a, two_b, two_e, CSI_PREFACT_FPF(csi));
  FPF_DUMP("cum_prefact_fpf", CSI_PREFACT_FPF(csi));
  delta_coeff(two_c, two_d, two_e, CSI_PREFACT_FPF(csi));
  FPF_DUMP("cum_prefact_fpf", CSI_PREFACT_FPF(csi));
  delta_coeff(two_a, two_c, two_f, CSI_PREFACT_FPF(csi));
  FPF_DUMP("cum_prefact_fpf", CSI_PREFACT_FPF(csi));
  delta_coeff(two_b, two_d, two_f, CSI_PREFACT_FPF(csi));
  FPF_DUMP("cum_prefact_fpf", CSI_PREFACT_FPF(csi));
}

void calcsum_9j(struct wigxjpf_temp *csi,
		int two_a, int two_b, int two_c,
		int two_d, int two_e, int two_f,
		int two_g, int two_h, int two_i)
{
  int two_kmin = max3(abs(two_h - two_d),
		      abs(two_b - two_f),
		      abs(two_a - two_i));
  int two_kmax = min3(two_h + two_d,
		      two_b + two_f,
		      two_a + two_i);
  int two_k;

  CHECK_TEMP;

#if DEBUG_PRINT
  printf ("[9j] krange: [%d,%d] %d steps\n",
	  two_kmin, two_kmax,
	  (two_kmax - two_kmin) / 2 + 1);
#endif

  /* One could think of doing a first round over the loop, to find out
   * the least common denominator.  That will however not help very
   * far, as we then just have to explicity multiply it in in each
   * term anyhow.  We may as well do just one loop, and calculate each
   * term (product of the non-prefactor part of three 6js).  For each
   * such term, we then scale it (and the previous sum) to the new
   * common denominator.
   */

  pexpo_set_zero(CSI_MIN_NUME_FPF(csi), 0);
  mwi_set_one_word(&csi->sum_prod, 0);

  for (two_k = two_kmin; two_k <= two_kmax; two_k += 2)
    {
      factor_6j(csi, two_a, two_b, two_c, two_f, two_i, two_k,
		CSI_TRIPROD_Fx_FPF(csi,0), &csi->triprod);

      factor_6j(csi, two_f, two_d, two_e, two_h, two_b, two_k,
		CSI_TRIPROD_Fx_FPF(csi,1), &csi->triprod_factor);

      mwi_mul_mwi(&csi->triprod_tmp, &csi->triprod, &csi->triprod_factor);
      MWI_DUMP("triprod_tmp", &csi->triprod_tmp);

      factor_6j(csi, two_h, two_i, two_g, two_a, two_d, two_k,
		CSI_TRIPROD_Fx_FPF(csi,2), &csi->triprod_factor);

      mwi_mul_mwi(&csi->triprod, &csi->triprod_tmp, &csi->triprod_factor);
      MWI_DUMP("triprod", &csi->triprod);

      FPF_DUMP("nume_triprod_f1_fpf", CSI_TRIPROD_Fx_FPF(csi,0));
      FPF_DUMP("nume_triprod_f2_fpf", CSI_TRIPROD_Fx_FPF(csi,1));
      FPF_DUMP("nume_triprod_f3_fpf", CSI_TRIPROD_Fx_FPF(csi,2));

      pexpo_expand_sum3(CSI_NUME_TRIPROD_FPF(csi),
			CSI_TRIPROD_Fx_FPF(csi,0),
			CSI_TRIPROD_Fx_FPF(csi,1),
			CSI_TRIPROD_Fx_FPF(csi,2));

      FPF_DUMP("nume_triprod_fpf", CSI_NUME_TRIPROD_FPF(csi));

      /* And multiply by the inner delta coefficients. */

      delta_coeff(two_a, two_i, two_k, CSI_NUME_TRIPROD_FPF(csi));
      FPF_DUMP("nume_triprod_fpf", CSI_NUME_TRIPROD_FPF(csi));
      delta_coeff(two_f, two_b, two_k, CSI_NUME_TRIPROD_FPF(csi));
      FPF_DUMP("nume_triprod_fpf", CSI_NUME_TRIPROD_FPF(csi));
      delta_coeff(two_h, two_d, two_k, CSI_NUME_TRIPROD_FPF(csi));
      FPF_DUMP("nume_triprod_fpf", CSI_NUME_TRIPROD_FPF(csi));

      /* And multiply by the coefficient 2*k+1. */

      const struct prime_exponents *p_f1 =
	PRIME_FACTOR(two_k + 1);

      pexpo_expand_add(CSI_NUME_TRIPROD_FPF(csi), p_f1);
      FPF_DUMP("nume_triprod_fpf", CSI_NUME_TRIPROD_FPF(csi));

      /* Before we can add triprod to sum_triprod,
       * we must normalise both to the least common denominator so far.
       */

      if (two_k == two_kmin)
	{
	  pexpo_copy(CSI_MIN_NUME_FPF(csi), CSI_NUME_TRIPROD_FPF(csi));
	  mwi_set_one_word(&csi->big_nume, 1);
	  mwi_set_one_word(&csi->big_div, 1);
	}
      else
	{
	  FPF_DUMP("min_nume_triprod_fpf", CSI_NUME_TRIPROD_FPF(csi));

	  pexpo_expand_blocks(CSI_MIN_NUME_FPF(csi),
			      CSI_NUME_TRIPROD_FPF(csi)->num_blocks);
	  pexpo_keep_min_in_as_diff(CSI_MIN_NUME_FPF(csi),
				    CSI_NUME_TRIPROD_FPF(csi));

	  FPF_DUMP("min_nume_triprod_fpf", CSI_NUME_TRIPROD_FPF(csi));
	  FPF_DUMP("nume_triprod_fpf", CSI_NUME_TRIPROD_FPF(csi));

	  pexpo_evaluate2(&csi->big_div, &csi->big_nume,
			  CSI_NUME_TRIPROD_FPF(csi),
			  &csi->pexpo_eval);
	}

      MWI_DUMP("triprod", &csi->triprod);
      MWI_DUMP("sum_triprod", &csi->sum_prod);

      MWI_DUMP("big_nume", &csi->big_nume);
      MWI_DUMP("big_div", &csi->big_div);

      /* We abuse as temporary: */
      struct multi_word_int *sum_triprod_tmp = &csi->triprod_factor;

      mwi_mul_mwi(&csi->triprod_tmp, &csi->triprod,  &csi->big_div);
      mwi_mul_mwi(sum_triprod_tmp,   &csi->sum_prod, &csi->big_nume);

      MWI_DUMP("triprod_tmp", &csi->triprod_tmp);
      MWI_DUMP("sum_triprod_tmp", sum_triprod_tmp);

      mwi_copy(&csi->sum_prod, sum_triprod_tmp);

      if ((two_k) & 1)
	mwi_sub_mwi(&csi->sum_prod, &csi->triprod_tmp);
      else
	mwi_add_mwi(&csi->sum_prod, &csi->triprod_tmp);

      MWI_DUMP("sum_triprod", &csi->sum_prod);
   }

  pexpo_set_zero(CSI_PREFACT_FPF(csi), 0);

  delta_coeff(two_a, two_b, two_c, CSI_PREFACT_FPF(csi));
  FPF_DUMP("cum_prefact_fpf", CSI_PREFACT_FPF(csi));
  delta_coeff(two_d, two_e, two_f, CSI_PREFACT_FPF(csi));
  FPF_DUMP("cum_prefact_fpf", CSI_PREFACT_FPF(csi));
  delta_coeff(two_g, two_h, two_i, CSI_PREFACT_FPF(csi));
  FPF_DUMP("cum_prefact_fpf", CSI_PREFACT_FPF(csi));
  delta_coeff(two_a, two_d, two_g, CSI_PREFACT_FPF(csi));
  FPF_DUMP("cum_prefact_fpf", CSI_PREFACT_FPF(csi));
  delta_coeff(two_b, two_e, two_h, CSI_PREFACT_FPF(csi));
  FPF_DUMP("cum_prefact_fpf", CSI_PREFACT_FPF(csi));
  delta_coeff(two_c, two_f, two_i, CSI_PREFACT_FPF(csi));
  FPF_DUMP("cum_prefact_fpf", CSI_PREFACT_FPF(csi));
}

#define DOUBLE_TYPE         double
#define DBL_FCN_POSTFIX(x)  x##_double
#define DBL_MATH_FCN_LQ(x)  x
#include "calc_dbl.h"
#undef DBL_MATH_FCN_LQ
#undef DBL_FCN_POSTFIX
#undef DOUBLE_TYPE

#if WIGXJPF_IMPL_LONG_DOUBLE
# define DOUBLE_TYPE         long double
# define DBL_FCN_POSTFIX(x)  x##_long_double
# define DBL_MATH_FCN_LQ(x)  x##l
# include "calc_dbl.h"
# undef DBL_MATH_FCN_LQ
# undef DBL_FCN_POSTFIX
# undef DOUBLE_TYPE
#endif

/* */

struct wigxjpf_temp *wigxjpf_temp_alloc(int max_iter)
{
  size_t cache_alignment = 64;
  size_t prime_factor_size;
  size_t total_size;

  if (max_iter < 1)
    max_iter = 1;

  prime_factor_size =
    (CSI_PF_K_ITER_START + (size_t) max_iter) * wigxjpf_prime_fact_stride;

  total_size = sizeof (struct wigxjpf_temp) +
    prime_factor_size + cache_alignment;

  struct wigxjpf_temp *temp = (struct wigxjpf_temp *) malloc (total_size);

  if (temp == NULL)
    {
      fprintf (stderr,
	       "wigxjpf: Memory allocation error (wigxjpf_temp), %zd bytes.\n",
	       total_size);
      wigxjpf_error();
    }

  size_t align_pf = (size_t) (temp + 1);
  align_pf += (cache_alignment - 1);
  align_pf -= (align_pf) % cache_alignment;

  temp->prime_exp_base = (void *) align_pf;
  temp->max_iter = max_iter;
  temp->size = total_size;

  mwi_alloc(&temp->sum_prod);

  mwi_alloc(&temp->big_prod);

  mwi_alloc(&temp->big_sqrt);
  mwi_alloc(&temp->big_nume);
  mwi_alloc(&temp->big_div);
  mwi_alloc(&temp->big_nume_prod);

  mwi_alloc(&temp->triprod);
  mwi_alloc(&temp->triprod_tmp);
  mwi_alloc(&temp->triprod_factor);

  int i;

  for (i = 0; i < 2; i++)
    {
      mwi_alloc(&temp->pexpo_eval.prod_pos[i]);
      mwi_alloc(&temp->pexpo_eval.prod_neg[i]);
      mwi_alloc(&temp->pexpo_eval.factor[i]);
      mwi_alloc(&temp->pexpo_eval.big_up[i]);
    }

  temp->inuse = 0;

  return temp;
}

void wigxjpf_temp_free(struct wigxjpf_temp *temp)
{
  if (temp == NULL)
    return;

  if (temp->inuse)
    {
      fprintf (stderr,
	       "wigxjpf: Freeing temp array while still in use.\n");
      wigxjpf_error();
    }

  mwi_free(&temp->sum_prod);

  mwi_free(&temp->big_prod);

  mwi_free(&temp->big_sqrt);
  mwi_free(&temp->big_nume);
  mwi_free(&temp->big_div);
  mwi_free(&temp->big_nume_prod);

  mwi_free(&temp->triprod);
  mwi_free(&temp->triprod_tmp);
  mwi_free(&temp->triprod_factor);

  int i;

  for (i = 0; i < 2; i++)
    {
      mwi_free(&temp->pexpo_eval.prod_pos[i]);
      mwi_free(&temp->pexpo_eval.prod_neg[i]);
      mwi_free(&temp->pexpo_eval.factor[i]);
      mwi_free(&temp->pexpo_eval.big_up[i]);
    }

  free(temp);
}

size_t wigxjpf_temp_size(struct wigxjpf_temp *temp)
{
  size_t size;

  size = temp->size;

  size += mwi_alloc_size(&temp->sum_prod);

  size += mwi_alloc_size(&temp->big_prod);

  size += mwi_alloc_size(&temp->big_sqrt);
  size += mwi_alloc_size(&temp->big_nume);
  size += mwi_alloc_size(&temp->big_div);
  size += mwi_alloc_size(&temp->big_nume_prod);

  size += mwi_alloc_size(&temp->triprod);
  size += mwi_alloc_size(&temp->triprod_tmp);
  size += mwi_alloc_size(&temp->triprod_factor);

  int i;

  for (i = 0; i < 2; i++)
    {
      size += mwi_alloc_size(&temp->pexpo_eval.prod_pos[i]);
      size += mwi_alloc_size(&temp->pexpo_eval.prod_neg[i]);
      size += mwi_alloc_size(&temp->pexpo_eval.factor[i]);
      size += mwi_alloc_size(&temp->pexpo_eval.big_up[i]);
    }

  return size;
}
