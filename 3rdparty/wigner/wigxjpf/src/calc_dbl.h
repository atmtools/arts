
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

/* This file is included several times from wigxjpf_calc.c
 * to implement for different floating point types.
 */

void DBL_FCN_POSTFIX(eval_calcsum_info)(DOUBLE_TYPE *result,
					struct wigxjpf_temp *csi)
{
  /* The prefactor shall suffer a sqrt().  So split it into parts
   * which are divisible by 2, and the remainders (that need a sqrt).
   * The remainders are kept as pure numerators (denominator better?,
   * due to more denominator outside needing cancelling?), and the
   * rest is moved out.
   */

  pexpo_split_sqrt_add(CSI_PREFACT_FPF(csi),
		       &csi->big_sqrt, CSI_MIN_NUME_FPF(csi));

  MWI_DUMP("big_sqrt", &csi->big_sqrt);
  FPF_DUMP("prefact_fpf", CSI_PREFACT_FPF(csi));

  pexpo_evaluate2(&csi->big_nume, &csi->big_div, CSI_PREFACT_FPF(csi),
		  &csi->pexpo_eval);

  MWI_DUMP("big_nume", &csi->big_nume);
  MWI_DUMP("big_div", &csi->big_div);

  mwi_mul_mwi(&csi->big_nume_prod, &csi->big_nume, &csi->sum_prod);

  MWI_DUMP("big_nume_prod", &csi->big_nume_prod);

  DOUBLE_TYPE d_nume_prod;
  DOUBLE_TYPE d_div;
  DOUBLE_TYPE d_sqrt;

  int exp_nume_prod;
  int exp_div;
  int exp_sqrt;

  DBL_FCN_POSTFIX(mwi_to)(&d_nume_prod, &exp_nume_prod, &csi->big_nume_prod);
  DBL_FCN_POSTFIX(mwi_to)(&d_div,       &exp_div,       &csi->big_div);
  DBL_FCN_POSTFIX(mwi_to)(&d_sqrt,      &exp_sqrt,      &csi->big_sqrt);

  /* Test with 3j shows a speed advantage for dividing with the sqrt
   * factor.  About 4.5 % at 2j ~ 100, and less for other cases.
   */
  DOUBLE_TYPE r =
    ((d_nume_prod) / d_div) / DBL_MATH_FCN_LQ(sqrt) (d_sqrt);

  /* Note: exp_sqrt is always en even number. */
  int res_exponent = exp_nume_prod - exp_div - exp_sqrt / 2;

  *result = DBL_MATH_FCN_LQ(ldexp)(r, res_exponent);

#if DEBUG_PRINT
  printf ("%20s: %20.15Lf * 2^%d\n", "d_nume_prod",
	  (long double) d_nume_prod, exp_nume_prod);
  printf ("%20s: %20.15Lf * 2^%d\n", "d_div",
	  (long double) d_div, exp_div);
  printf ("%20s: %20.15Lf * 2^%d\n", "d_sqrt(pre)",
	  (long double) d_sqrt, exp_sqrt);
  printf ("%20s: %20.15Lf * 2^%d\n", "d_sqrt",
	  (long double) DBL_MATH_FCN_LQ(sqrt) (d_sqrt), exp_sqrt / 2);
  printf ("%20s: %20.15Lf\n", "result", (long double) *result);
#endif
}

void DBL_FCN_POSTFIX(calc_3j)(DOUBLE_TYPE *result,
			      int two_j1, int two_j2, int two_j3,
			      int two_m1, int two_m2, int two_m3,
			      struct wigxjpf_temp *temp)
{
  if (trivial_zero_3j(two_j1, two_j2, two_j3,
		      two_m1, two_m2, two_m3))
    {
      *result = 0;
      return;
    }

  calcsum_3j(temp,
	     two_j1, two_j2, two_j3,
	     two_m1, two_m2, two_m3);

  DBL_FCN_POSTFIX(eval_calcsum_info)(result, temp);

  temp->inuse = 0;
}

void DBL_FCN_POSTFIX(calc_6j)(DOUBLE_TYPE *result,
			      int two_j1, int two_j2, int two_j3,
			      int two_j4, int two_j5, int two_j6,
			      struct wigxjpf_temp *temp)
{
  if (trivial_zero_6j(two_j1, two_j2, two_j3,
		      two_j4, two_j5, two_j6))
    {
      *result = 0;
      return;
    }

  calcsum_6j(temp,
	     two_j1, two_j2, two_j3,
	     two_j4, two_j5, two_j6);

  DBL_FCN_POSTFIX(eval_calcsum_info)(result, temp);

  temp->inuse = 0;
}

void DBL_FCN_POSTFIX(calc_9j)(DOUBLE_TYPE *result,
			      int two_j1, int two_j2, int two_j3,
			      int two_j4, int two_j5, int two_j6,
			      int two_j7, int two_j8, int two_j9,
			      struct wigxjpf_temp *temp)
{
  if (trivial_zero_9j(two_j1, two_j2, two_j3,
		      two_j4, two_j5, two_j6,
		      two_j7, two_j8, two_j9))
    {
      *result = 0;
      return;
    }
  calcsum_9j(temp,
	     two_j1, two_j2, two_j3,
	     two_j4, two_j5, two_j6,
	     two_j7, two_j8, two_j9);

  DBL_FCN_POSTFIX(eval_calcsum_info)(result, temp);

  temp->inuse = 0;
}
