
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

#include "prime_factor.h"

#include <string.h>
#include <stddef.h>

uint32_t *wigxjpf_prime_list = NULL;

size_t wigxjpf_prime_fact_stride = 0;
int wigxjpf_max_prime_decomp = -1;
void *wigxjpf_prime_factors_base = NULL;
void *wigxjpf_prime_factors_1 = NULL;
void *wigxjpf_prime_factors_2 = NULL;

size_t wigxjpf_fill_factors(int max_factorial)
{
  int i;
  int p;

  if (max_factorial == 0)
    {
      /* Free memory. */

      free(wigxjpf_prime_list);
      wigxjpf_prime_list = NULL;
      wigxjpf_max_prime_decomp = -1;

      free(wigxjpf_prime_factors_base);
      wigxjpf_prime_factors_base = NULL;

      wigxjpf_prime_factors_1 = NULL;
      wigxjpf_prime_factors_2 = NULL;
      
      return 0;
    }

  if (max_factorial < 2) /* Otherwise routine will hang. */
    max_factorial = 2;

  /* Simple sieve of Eratosthenes to find all prime numbers up to
   * max_factorial.
   */
  
  char *isprime;
  size_t isprime_size = (size_t) (max_factorial+1) * sizeof (char);

  isprime = (char *) malloc (isprime_size);

  if (isprime == NULL)
    {
      fprintf (stderr,
	       "wigxjpf: "
	       "Memory allocation error (isprimearray), %zd bytes.\n",
	       isprime_size);
      exit(1);
    }

  memset(isprime, 1, isprime_size);

  for (p = 2; p * p <= max_factorial; p++)
    {
      int mark;

      if (!isprime[p])
	continue;

      for (mark = p*p; mark <= max_factorial; mark += p)
	isprime[mark] = 0;
    }

  /* We can now count the number of prime numbers, i.e. the size of
   * out basis.
   */

  int num_primes = 0;

  for (i = 2; i <= max_factorial; i++)
    num_primes += isprime[i];

  int num_blocks = (num_primes + (PRIME_LIST_VECT_BLOCK_ITEMS - 1)) /
    PRIME_LIST_VECT_BLOCK_ITEMS;

  /* List of primes. */

  wigxjpf_prime_list =
    (uint32_t *) realloc (wigxjpf_prime_list,
			  (size_t) num_primes * sizeof (uint32_t));

  if (wigxjpf_prime_list == NULL)
    {
      fprintf (stderr,
	       "wigxjpf: "
	       "Memory allocation error (prime list), %zd bytes.\n",
	       sizeof (uint32_t));
      exit(1);
    }

  p = 0;

  for (i = 2; i <= max_factorial; i++)
    if (isprime[i])
      wigxjpf_prime_list[p++] = (uint32_t) i;

  /* We will store the decomposition of all numbers (including 0) up
   * to and including @max_factorial, i.e. max_factorial+1 numbers.
   */

  wigxjpf_prime_fact_stride =
    offsetof (struct prime_exponents, expo) +
    (size_t) num_blocks * PRIME_LIST_VECT_BLOCK_ITEMS * sizeof (prime_exp_t);

  /* We want each decomposition to start on a new cache line.
   * (It is seldom that the last part of a decomposition is read).
   */

  size_t cache_alignment = 64;

  wigxjpf_prime_fact_stride += cache_alignment - 1;
  wigxjpf_prime_fact_stride -= (wigxjpf_prime_fact_stride) % cache_alignment;

  size_t prime_factor_size =
    (size_t) (max_factorial + 1) * wigxjpf_prime_fact_stride;

  size_t prime_factors_base_size =
    2 * prime_factor_size + cache_alignment;

  wigxjpf_prime_factors_base =
    realloc (wigxjpf_prime_factors_base,
	     prime_factors_base_size);

  if (wigxjpf_prime_factors_base == NULL)
    {
      fprintf (stderr,
	       "wigxjpf: "
	       "Memory allocation error (prime factors), %zd bytes.\n",
	       prime_factors_base_size);
      exit(1);
    }

  size_t align_pf = (size_t) wigxjpf_prime_factors_base;
  align_pf += (cache_alignment - 1);
  align_pf -= (align_pf) % cache_alignment;

  wigxjpf_prime_factors_1 = (void *)  align_pf;
  wigxjpf_prime_factors_2 = (void *) (align_pf + prime_factor_size);

  /* We now fill the prime_factors, by enumerating the base...
   * We play around with item 0 temporarily.
   */

  struct prime_exponents *tmp = PRIME_FACTOR(0);

  memset(tmp, 0, wigxjpf_prime_fact_stride);
  memset(PRIME_FACTOR(1), 0, wigxjpf_prime_fact_stride);

  uint64_t cur = 1;
  int max_p = 0;

  for ( ; ; )
    {
      /* We always try to increment the smallest exponent. */
      p = 0;

      for ( ; ; )
	{
	  /* If incrementing the current exponent by one more exceeds
	   * our destination table, we reset this exponent to zero,
	   * and try the next.
	   */
	  
	  if (cur * wigxjpf_prime_list[p] <= (uint64_t) max_factorial)
	    {
	      tmp->expo[p]++;
	      cur *= wigxjpf_prime_list[p];
	      break; /* successful */
	    }

	  /* This exponent did overflow.  Remove it. */
	  while (tmp->expo[p])
	    {
	      cur /= wigxjpf_prime_list[p];
	      tmp->expo[p]--;	      
	    }
	  /* We shall try the next exponent. */
	  p++;
	  if (p > max_p)
	    max_p = p;
	  if (p >= num_primes)
	    goto done_enum;
	}

      struct prime_exponents *dest = PRIME_FACTOR(cur);

      memcpy (dest, tmp, wigxjpf_prime_fact_stride);

      dest->num_blocks = max_p / PRIME_LIST_VECT_BLOCK_ITEMS + 1;
    }
 done_enum:
  /* Reset the temporary item. */
  memset (tmp, 0, wigxjpf_prime_fact_stride);

  memset (FACTORIAL_PRIME_FACTOR(0), 0, wigxjpf_prime_fact_stride);

  /* The factorial tables are set by accumulation. */
  for (i = 1; i <= max_factorial; i++)
    {
      struct prime_exponents *add  = PRIME_FACTOR(i);
      struct prime_exponents *src  = FACTORIAL_PRIME_FACTOR(i-1);
      struct prime_exponents *dest = FACTORIAL_PRIME_FACTOR(i);

      for (p = 0; p < num_primes; p++)
	dest->expo[p] = src->expo[p] + add->expo[p];

      dest->num_blocks = src->num_blocks;
      if (add->num_blocks > dest->num_blocks)
	dest->num_blocks = add->num_blocks;
    }

  free(isprime);

  /* We need a larger type for the exponents if there is any
   * danger of overflow.
   *
   * The worst formula is 9j, which contain 
   * 3*(factor6j:(1add+7sub)+delta:(3add+1sub)+coeff:1add) and
   * 6*(delta:(3add+1sub)).  I.e. 3*(5add+8sub)=15add+8sub
   * and 18add+6sub.  Worst case 35 added.
   *
   * Hmm, lets require a margin of 50.
   *
   * And we cannot use last bit,
   */

  uint32_t max_allowed = (uint32_t) (1 << (sizeof (prime_exp_t) * 8 - 1));

  /* We can never pick up more than @max_factorial factors of [2] during the
   * factorials.  I.e. FACTORIAL_PRIME_FACTOR(i)->expo[0] is always < i.
   */
  
  if ((uint32_t) max_factorial * /*50*/5 > max_allowed)
    {
      fprintf (stderr,
	       "wigxjpf: "
	       "Type prime_exp_t too small!  "
	       "Exponent for [2] could overflow.\n");
      exit(1); 
    }

  wigxjpf_max_prime_decomp = max_factorial;

  return prime_factors_base_size;
};
