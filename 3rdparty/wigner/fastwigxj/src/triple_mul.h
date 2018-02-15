#include <stdio.h>
#include <stdint.h>

typedef __uint128_t  uint128_t;
typedef __int128_t   int128_t;

#if 0
#define DEBUG(fmt,...)  printf(fmt,__VA_ARGS__);
#else
#define DEBUG(fmt,...)  do { } while (0)
#endif
    
#define MUL_2w64_2w64_GET_HI2w64(hi_a,lo_a,hi_b,lo_b,hi_res,lo_res) do { \
    uint128_t lalb = lo_a * (uint128_t) lo_b;				\
    uint128_t lahb = lo_a * (uint128_t) hi_b;				\
    uint128_t halb = hi_a * (uint128_t) lo_b;				\
    uint128_t hahb = hi_a * (uint128_t) hi_b;				\
									\
    /* uint64_t lalb_l = (uint64_t) lalb; */				\
    uint64_t lahb_l = (uint64_t) lahb;					\
    uint64_t halb_l = (uint64_t) halb;					\
    uint64_t hahb_l = (uint64_t) hahb;					\
									\
    uint64_t lalb_h = (uint64_t) (lalb >> 64);				\
    uint64_t lahb_h = (uint64_t) (lahb >> 64);				\
    uint64_t halb_h = (uint64_t) (halb >> 64);				\
    uint64_t hahb_h = (uint64_t) (hahb >> 64);				\
									\
    DEBUG("lalb: %016" PRIx64 ":----------------\n",lalb_h);		\
    DEBUG("lahb: %016" PRIx64 ":%016" PRIx64 "\n",lahb_h,lahb_l);	\
    DEBUG("halb: %016" PRIx64 ":%016" PRIx64 "\n",halb_h,halb_l);	\
    DEBUG("hahb: %016" PRIx64 ":%016" PRIx64 "\n",hahb_h,hahb_l);	\
									\
    /* uint64_t  ab_0 = lalb_l; */					\
    uint128_t ab_1 =							\
      ((uint128_t) lalb_h) + ((uint128_t) lahb_l) + ((uint128_t) halb_l); \
    uint64_t ab_1_carry = (uint64_t) (ab_1 >> 64);			\
    uint128_t ab_2 =							\
      ((uint128_t) hahb_l) + ((uint128_t) lahb_h) + ((uint128_t) halb_h) + \
      ((uint128_t) ab_1_carry);						\
    uint64_t ab_2_carry = (uint64_t) (ab_2 >> 64);			\
    uint128_t ab_3 = ((uint128_t) hahb_h) + ((uint128_t) ab_2_carry);	\
									\
    lo_res = (uint64_t) ab_2;						\
    hi_res = (uint64_t) ab_3;						\
  } while (0)

#define MUL_2w64_1w64_GET_LO2w64(hi_a,lo_a,lo_b,hi_res,lo_res) do { \
    uint128_t lalb = lo_a * (uint128_t) lo_b;				\
    uint128_t halb = hi_a * (uint128_t) lo_b;				\
									\
    uint64_t lalb_l = (uint64_t) lalb;					\
    uint64_t halb_l = (uint64_t) halb;					\
									\
    uint64_t lalb_h = (uint64_t) (lalb >> 64);				\
    /* uint64_t halb_h = (uint64_t) (halb >> 64); */			\
									\
    DEBUG("lalb: %016" PRIx64 ":%016" PRIx64 "\n",lalb_h,lalb_l);	\
    DEBUG("halb: ----------------:%016" PRIx64 "\n",halb_l);		\
									\
    uint64_t  ab_0 = lalb_l;						\
    uint128_t ab_1 =							\
      ((uint128_t) lalb_h) + ((uint128_t) halb_l);			\
    /*uint64_t ab_1_carry = (uint64_t) (ab_1 >> 64);*/			\
									\
    lo_res = (uint64_t) ab_0;						\
    hi_res = (uint64_t) ab_1;						\
  } while (0)

void triple_mul_accum(uint64_t *a, uint64_t *b, uint64_t *c,
		      uint64_t d,
		      int128_t *accum, int32_t *accum_exp)
{
  /* We have 112 bits of fraction in each of a, b, and c.  And would
   * like to end up having that after being done as well, since we
   * would like the higher bits to be useful for the accumulation.
   * Simplest is to shift a and b up, and do c without shifting.  Also
   * note that a, b and c have an implicit 1, which we need to keep
   * track of.  Another issue is the sign bit.  Since applying the
   * sign requires us to subtract from 0, it is easier to only do that
   * once, i.e. on the product.
   */

  uint32_t exp_a = (uint32_t) (a[1] >> 48);
  uint32_t exp_b = (uint32_t) (b[1] >> 48);
  uint32_t exp_c = (uint32_t) (c[1] >> 48);

  uint32_t sign_bit_15 = exp_a ^ exp_b ^ exp_c;

  int32_t exp_abcd = (int32_t)
    ((exp_a & 0x7fff) + (exp_b & 0x7fff) + (exp_c & 0x7fff));

  DEBUG("exps: %4x %4x %4x -> exp_abc %5x sign %d\n",
	exp_a, exp_b, exp_c, exp_abcd, sign_bit_15 >> 15);

  /* If any of a, b or c is zero, we cannot follow the normal path,
   * since that then inserts a wrong (implicit) leading one.  Also, it
   * would move around the exponent of the accumulated result, leaving
   * it harder to find out if we have had massive precision loss.
   *
   * When detecting a zero we ca cheat a bit however.  We know that no
   * value is greater than 1.  I.e. has exponent > 16383.  We also know
   * that no value (that is not exactly zero) is very small.  So we say
   * that if the product exponent is smaller than having two ones, one
   * of them was zero.
   */

  if (exp_abcd < 2*16383)
    return;  

  /* We shift b one less than a, since a*b may carry one bit that
   * now ends up at the highest bit (and does not get lost).
   */

  uint64_t hi_a = 0x8000000000000000ll | (a[1] << 15) | (a[0] >> 49);
  uint64_t hi_b = 0x8000000000000000ll | (b[1] << 15) | (b[0] >> 49);
  /*
  uint64_t hi_b = 0x4000000000000000ll | ((b[1] & 0x0000ffffffffffffll)
					  << 14) | (b[0] >> 50);
  */
  uint64_t lo_a = (a[0] << 15);
  /*uint64_t lo_b = (b[0] << 14);*/
  uint64_t lo_b = (b[0] << 15);

  uint64_t lo_ab, hi_ab;

  DEBUG("  a: %016" PRIx64 ":%016" PRIx64 "  "
	"  b: %016" PRIx64 ":%016" PRIx64 "\n",
	hi_a, lo_a, hi_b, lo_b);

  MUL_2w64_2w64_GET_HI2w64(hi_a, lo_a, hi_b, lo_b, hi_ab, lo_ab);

  uint64_t hi_c =
    0x0010000000000000ll | ((c[1] & 0x0000ffffffffffffll) << 4) | (c[0] >> 60);
  uint64_t lo_c = c[0] << 4;

  uint64_t lo_cd, hi_cd;

  DEBUG("  c: %016" PRIx64 ":%016" PRIx64 "  "
	"  d: %016" PRIx64 ":%016" PRIx64 "\n",
	hi_c, lo_c, 0, d);

  MUL_2w64_1w64_GET_LO2w64(hi_c, lo_c, d, hi_cd, lo_cd);

  DEBUG(" ab: %016" PRIx64 ":%016" PRIx64 "  "
	" cd: %016" PRIx64 ":%016" PRIx64 "\n",
	hi_ab, lo_ab, hi_cd, lo_cd);
  
  uint64_t lo_abcd, hi_abcd;

  MUL_2w64_2w64_GET_HI2w64(hi_ab, lo_ab, hi_cd, lo_cd, hi_abcd, lo_abcd);

  DEBUG("abcd: %016" PRIx64 ":%016" PRIx64 "\n",
	hi_abcd, lo_abcd);

  int128_t abcd =
    (int128_t) (((uint128_t) lo_abcd) | (((uint128_t) hi_abcd) << 64));

  if (sign_bit_15 & 0x8000) {
    abcd = -abcd;
  }

  /* d will have 'shifted' the value up a number of bits, but that
   * does not improve our number of significant bits.  Shift them out
   * again, as the detection of accuracy loss due to cancellation
   * otherwise is fooled.  Also, this 'ensures' that the high free
   * bits we have can handle the accumulations.  If d is one, it has
   * one bit, and we do not want to shift.
   */

  int d_bits = (int) (sizeof (int) * 8 - 1) - __builtin_clz((unsigned int) d);

  abcd >>= d_bits;

  exp_abcd += d_bits;

  int32_t exp_diff = *accum_exp - exp_abcd;

  DEBUG("exp: %d %d => %d\n", *accum_exp, exp_abcd, exp_diff);

  int128_t acc = accum[0]; // | (((uint128_t) accum[1]) << 64);

  if (exp_diff > 0)
    {
      /* Accumulated result has larger exponent.  Shift abc down. */
      if (exp_diff < 128)
	{
	  abcd >>= exp_diff;
	  acc += abcd;
	}
      else
	{ } /* Nothing to add, abc would be completely shifted out. */
    }
  else
    {
      /* Accumulated result has smaller exponent.  Shift accum down. */
      if (exp_diff > -128)
	{
	  acc >>= (-exp_diff);
	  acc += abcd;
	  *accum_exp = exp_abcd;
	}
      else
	{
	  /* Accumulated result completely shifted out. */
	  acc = abcd;
	  *accum_exp = exp_abcd;
	}
    }

  //accum[1] = (uint64_t) acc;
  //accum[0] = (uint64_t) (acc >> 64);
  accum[0] = acc;
}

