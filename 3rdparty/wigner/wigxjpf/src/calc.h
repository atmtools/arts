
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

#ifndef __WIGXJPF_CALC_H__
#define __WIGXJPF_CALC_H__

#include <stdlib.h>

#include "wigxjpf_config.h"

#include "multi_word_int.h"
#include "prime_factor.h"

#if DEBUG_PRINT
#define MWI_DUMP(name, mwi)  do { mwi_dump(name, mwi); } while (0)
#define FPF_DUMP(name, fpf)  do { dump_fpf(name, fpf, -1); } while (0)
#define FPF_DUMP_FACT(name, fpf)			\
  do {							\
    dump_fpf(name, fpf,					\
	     ((const char *) fpf -			\
	      (char *) wigxjpf_prime_factors_2) /	\
	     (ssize_t) wigxjpf_prime_fact_stride);	\
  } while (0)

static void dump_fpf(const char *name, const struct prime_exponents *fpf,
		     ssize_t fact)
{
  int j;

  printf ("%20s:", name);
  if (fact != -1)
    printf (" %zd!", fact);

  for (j = 0; j < fpf->num_blocks * PRIME_LIST_VECT_BLOCK_ITEMS; j++)
    {
      printf (" %d", fpf->expo[j]);
    }
  printf ("\n");
}
#else
#define FPF_DUMP(name, fpf)  do { } while (0)
#define FPF_DUMP_FACT(name, fpf)  do { } while (0)
#define MWI_DUMP(name, mwi)  do { } while (0)
#endif

struct wigxjpf_temp
{
  struct multi_word_int sum_prod;

  struct multi_word_int big_prod;

  struct multi_word_int big_sqrt;
  struct multi_word_int big_nume;
  struct multi_word_int big_div;
  struct multi_word_int big_nume_prod;

  struct multi_word_int triprod;
  struct multi_word_int triprod_tmp;
  struct multi_word_int triprod_factor;

  struct pexpo_eval_temp pexpo_eval;

  int inuse;

  int max_iter;
  size_t size;

  /* remaining space is for a list of
   * struct prime_exponents
   */

  struct prime_exponents *prime_exp_base;
};

#define CSI_PF_PREFACT_NO      0
#define CSI_PF_MIN_NUME_NO     1
#define CSI_PF_NUME_TRIPROD_NO 2
#define CSI_PF_TRIPROD_Fx_NO   3 /* ,4,5 */
#define CSI_PF_K_ITER_START    6

#define CSI_PRIME_FACTOR(csi,i)					    \
  ((struct prime_exponents*) (((char *) (csi)->prime_exp_base) +    \
                              (i) * wigxjpf_prime_fact_stride))

#define CSI_PREFACT_FPF(csi)       CSI_PRIME_FACTOR(csi,CSI_PF_PREFACT_NO)
#define CSI_MIN_NUME_FPF(csi)      CSI_PRIME_FACTOR(csi,CSI_PF_MIN_NUME_NO)
#define CSI_NUME_TRIPROD_FPF(csi)  CSI_PRIME_FACTOR(csi,CSI_PF_NUME_TRIPROD_NO)
#define CSI_TRIPROD_Fx_FPF(csi,i)  CSI_PRIME_FACTOR(csi,		\
						    CSI_PF_TRIPROD_Fx_NO+(i))
#define CSI_K_ITER_FPF(csi)        CSI_PRIME_FACTOR(csi,CSI_PF_K_ITER_START)

void calcsum_3j(struct wigxjpf_temp *csi,
		int two_j1, int two_j2, int two_j3,
		int two_m1, int two_m2, int two_m3);

void calcsum_6j(struct wigxjpf_temp *csi,
		int two_j1, int two_j2, int two_j3,
		int two_j4, int two_j5, int two_j6);

void calcsum_9j(struct wigxjpf_temp *csi,
		int two_a, int two_b, int two_c,
		int two_d, int two_e, int two_f,
		int two_g, int two_h, int two_i);

#endif/*__WIGXJPF_CALC_H__*/
