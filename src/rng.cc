/* Copyright (C) 2003-2012 Cory Davis <cory@met.ed.ac.uk>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   rng.cc
  \author Cory Davis <cory@met.ed.ac.uk>
  \date   2003-06-26 

  \brief  member functions of the Rng class and gsl_rng code

  The Rng class is a simple class that uses the gsl_rng_mt_19937 random number 
generator from the GNU Scientific Library <http://www.gnu.org/software/gsl/>.
  
  The period of this generator is 2^{19937} - 1.
*/

#include "rng.h"
#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include "arts.h"
#include "messages.h"

/*!
Constructor creates instance of gsl_rng of type gsl_rng_mt19937
*/
Rng::Rng() { r = gsl_rng_alloc(gsl_rng_mt19937); }

/*!
Destructor frees memory allocated to gsl_rng 
*/
Rng::~Rng() { gsl_rng_free(r); }

/*!
 Seeds the Rng with the integer argument. 

 Every seed is only used once. The provided seed is increased by 1 until an
 unused seed is found.
*/
void Rng::seed(unsigned long int n, const Verbosity &verbosity) {
  // Static pool of previously used seeds.
  static vector<unsigned long int> seeds;

#pragma omp critical(Rng_seed)
  {
    unsigned long int n_orig = n;
    //cout << "Requested seed: " << n;
    while (find(seeds.begin(), seeds.end(), n) != seeds.end()) {
      if (n >= ULONG_MAX - 1)
        n = 0;
      else
        n++;

      // If all possible seeds were already used, we empty the pool and
      // start over.
      if (n == n_orig) {
        CREATE_OUT0;
        out0
            << "Rng Warning: Couldn't find an unused seed. Clearing seed pool.\n";
        seeds.empty();
        break;
      }
    }
    seeds.push_back(n);
    seed_no = n;
    //cout << " Got seed: " << seed_no << endl;
  }

  gsl_rng_set(r, seed_no);
}

/*!
 Seeds the Rng with the integer argument. 
 */
void Rng::force_seed(unsigned long int n) {
  seed_no = n;
  gsl_rng_set(r, seed_no);
}

/*!
Draws a double from the uniform distribution [0,1)
*/
double Rng::draw() { return gsl_rng_uniform(r); }

/*!
Returns the seed number
*/
unsigned long int Rng::showseed() const { return seed_no; }

/* 
   rng/mt.c
   
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.  You should have received
   a copy of the GNU General Public License along with this program;
   if not, write to the Free Foundation, Inc., 59 Temple Place, Suite
   330, Boston, MA 02111-1307 USA

   Original implementation was copyright (C) 1997 Makoto Matsumoto and
   Takuji Nishimura. Coded by Takuji Nishimura, considering the
   suggestions by Topher Cooper and Marc Rieffel in July-Aug. 1997, "A
   C-program for MT19937: Integer version (1998/4/6)"

   This implementation copyright (C) 1998 Brian Gough. I reorganized
   the code to use the module framework of GSL.  The license on this
   implementation was changed from LGPL to GPL, following paragraph 3
   of the LGPL, version 2.

   Update:

   The seeding procedure has been updated to match the 10/99 release
   of MT19937.

   Update:

   The seeding procedure has been updated again to match the 2002
   release of MT19937

   The original code included the comment: "When you use this, send an
   email to: matumoto@math.keio.ac.jp with an appropriate reference to
   your work".

   Makoto Matsumoto has a web page with more information about the
   generator, http://www.math.keio.ac.jp/~matumoto/emt.html. 

   The paper below has details of the algorithm.

   From: Makoto Matsumoto and Takuji Nishimura, "Mersenne Twister: A
   623-dimensionally equidistributerd uniform pseudorandom number
   generator". ACM Transactions on Modeling and Computer Simulation,
   Vol. 8, No. 1 (Jan. 1998), Pages 3-30

   You can obtain the paper directly from Makoto Matsumoto's web page.

   The period of this generator is 2^{19937} - 1.

*/

static inline unsigned long int mt_get(void *vstate);
static double mt_get_double(void *vstate);
static void mt_set(void *state, unsigned long int s);

#define N 624 /* Period parameters */
#define M 397

/* most significant w-r bits */
static const unsigned long UPPER_MASK = 0x80000000UL;

/* least significant r bits */
static const unsigned long LOWER_MASK = 0x7fffffffUL;

typedef struct {
  unsigned long mt[N];
  int mti;
} mt_state_t;

static inline unsigned long mt_get(void *vstate) {
  mt_state_t *state = (mt_state_t *)vstate;

  unsigned long k;
  unsigned long int *const mt = state->mt;

#define MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)

  if (state->mti >= N) { /* generate N words at one time */
    int kk;

    for (kk = 0; kk < N - M; kk++) {
      unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAGIC(y);
    }
    for (; kk < N - 1; kk++) {
      unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ MAGIC(y);
    }

    {
      unsigned long y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAGIC(y);
    }

    state->mti = 0;
  }

  /* Tempering */

  k = mt[state->mti];
  k ^= (k >> 11);
  k ^= (k << 7) & 0x9d2c5680UL;
  k ^= (k << 15) & 0xefc60000UL;
  k ^= (k >> 18);

  state->mti++;

  return k;
}

static double mt_get_double(void *vstate) {
  return (double)mt_get(vstate) / 4294967296.0;
}

static void mt_set(void *vstate, unsigned long int s) {
  mt_state_t *state = (mt_state_t *)vstate;
  int i;

  if (s == 0) s = 4357; /* the default seed is 4357 */

  state->mt[0] = s & 0xffffffffUL;

  for (i = 1; i < N; i++) {
    /* See Knuth's "Art of Computer Programming" Vol. 2, 3rd
         Ed. p.106 for multiplier. */

    state->mt[i] =
        (1812433253UL * (state->mt[i - 1] ^ (state->mt[i - 1] >> 30)) + i);

    state->mt[i] &= 0xffffffffUL;
  }

  state->mti = i;
}

static const gsl_rng_type mt_type = {"mt19937",    /* name */
                                     0xffffffffUL, /* RAND_MAX  */
                                     0,            /* RAND_MIN  */
                                     sizeof(mt_state_t),
                                     &mt_set,
                                     &mt_get,
                                     &mt_get_double};

const gsl_rng_type *gsl_rng_mt19937 = &mt_type;
unsigned long int gsl_rng_default_seed = 0;
/* rng/types.c
 * 
 * Copyright (C) 2001 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#define N1 100

const gsl_rng_type *gsl_rng_generator_types[N1];

#define ADD(t)                        \
  {                                   \
    if (i == N1) abort();             \
    gsl_rng_generator_types[i] = (t); \
    i++;                              \
  };

const gsl_rng_type **gsl_rng_types_setup(void) {
  int i = 0;

  /*  ADD(gsl_rng_borosh13);
  ADD(gsl_rng_cmrg);
  ADD(gsl_rng_coveyou);
  ADD(gsl_rng_fishman18);
  ADD(gsl_rng_fishman20);
  ADD(gsl_rng_fishman2x);
  ADD(gsl_rng_gfsr4);
  ADD(gsl_rng_knuthran);
  ADD(gsl_rng_knuthran2);
  ADD(gsl_rng_lecuyer21);
  ADD(gsl_rng_minstd);
  ADD(gsl_rng_mrg);*/
  ADD(gsl_rng_mt19937);
  /*  ADD(gsl_rng_mt19937_1999);
  ADD(gsl_rng_mt19937_1998);
  ADD(gsl_rng_r250);
  ADD(gsl_rng_ran0);
  ADD(gsl_rng_ran1);
  ADD(gsl_rng_ran2);
  ADD(gsl_rng_ran3);
  ADD(gsl_rng_rand);
  ADD(gsl_rng_rand48);
  ADD(gsl_rng_random128_bsd);
  ADD(gsl_rng_random128_glibc2);
  ADD(gsl_rng_random128_libc5);
  ADD(gsl_rng_random256_bsd);
  ADD(gsl_rng_random256_glibc2);
  ADD(gsl_rng_random256_libc5);
  ADD(gsl_rng_random32_bsd);
  ADD(gsl_rng_random32_glibc2);
  ADD(gsl_rng_random32_libc5);
  ADD(gsl_rng_random64_bsd);
  ADD(gsl_rng_random64_glibc2);
  ADD(gsl_rng_random64_libc5);
  ADD(gsl_rng_random8_bsd);
  ADD(gsl_rng_random8_glibc2);
  ADD(gsl_rng_random8_libc5);
  ADD(gsl_rng_random_bsd);
  ADD(gsl_rng_random_glibc2);
  ADD(gsl_rng_random_libc5);
  ADD(gsl_rng_randu);
  ADD(gsl_rng_ranf);
  ADD(gsl_rng_ranlux);
  ADD(gsl_rng_ranlux389);
  ADD(gsl_rng_ranlxd1);
  ADD(gsl_rng_ranlxd2);
  ADD(gsl_rng_ranlxs0);
  ADD(gsl_rng_ranlxs1);
  ADD(gsl_rng_ranlxs2);
  ADD(gsl_rng_ranmar);
  ADD(gsl_rng_slatec);
  ADD(gsl_rng_taus);
  ADD(gsl_rng_taus2);
  ADD(gsl_rng_taus113);
  ADD(gsl_rng_transputer);
  ADD(gsl_rng_tt800);
  ADD(gsl_rng_uni);
  ADD(gsl_rng_uni32);
  ADD(gsl_rng_vax);
  ADD(gsl_rng_waterman14);
  ADD(gsl_rng_zuf);*/
  ADD(0);

  return gsl_rng_generator_types;
}

/* err/stream.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

FILE *gsl_stream = NULL;
gsl_stream_handler_t *gsl_stream_handler = NULL;

void gsl_stream_printf(const char *label,
                       const char *file,
                       int line,
                       const char *reason) {
  if (gsl_stream == NULL) {
    gsl_stream = stderr;
  }
  if (gsl_stream_handler) {
    (*gsl_stream_handler)(label, file, line, reason);
    return;
  }
  fprintf(gsl_stream, "gsl: %s:%d: %s: %s\n", file, line, label, reason);
}

gsl_stream_handler_t *gsl_set_stream_handler(
    gsl_stream_handler_t *new_handler) {
  gsl_stream_handler_t *previous_handler = gsl_stream_handler;
  gsl_stream_handler = new_handler;
  return previous_handler;
}

FILE *gsl_set_stream(FILE *new_stream) {
  FILE *previous_stream;
  if (gsl_stream == NULL) {
    gsl_stream = stderr;
  }
  previous_stream = gsl_stream;
  gsl_stream = new_stream;
  return previous_stream;
}

/* err/error.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

gsl_error_handler_t *gsl_error_handler = NULL;

static void no_error_handler(const char *reason,
                             const char *file,
                             int line,
                             int gsl_errno);

void gsl_error(const char *reason, const char *file, int line, int gsl_errno) {
  if (gsl_error_handler) {
    (*gsl_error_handler)(reason, file, line, gsl_errno);
    return;
  }

  gsl_stream_printf("ERROR", file, line, reason);

  fprintf(stderr, "Default GSL error handler invoked.\n");
  abort();
}

gsl_error_handler_t *gsl_set_error_handler(gsl_error_handler_t *new_handler) {
  gsl_error_handler_t *previous_handler = gsl_error_handler;
  gsl_error_handler = new_handler;
  return previous_handler;
}

gsl_error_handler_t *gsl_set_error_handler_off(void) {
  gsl_error_handler_t *previous_handler = gsl_error_handler;
  gsl_error_handler = no_error_handler;
  return previous_handler;
}

static void no_error_handler(const char *reason _U_,
                             const char *file _U_,
                             int line _U_,
                             int gsl_errno _U_) {
  /* do nothing */
  return;
}

/* rng/rng.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

gsl_rng *gsl_rng_alloc(const gsl_rng_type *T) {
  gsl_rng *r = (gsl_rng *)malloc(sizeof(gsl_rng));

  if (r == 0) {
    GSL_ERROR_VAL("failed to allocate space for rng struct", GSL_ENOMEM, 0);
  };

  r->state = malloc(T->size);

  if (r->state == 0) {
    free(r); /* exception in constructor, avoid memory leak */

    GSL_ERROR_VAL("failed to allocate space for rng state", GSL_ENOMEM, 0);
  };

  r->type = T;

  gsl_rng_set(r, gsl_rng_default_seed); /* seed the generator */

  return r;
}

int gsl_rng_memcpy(gsl_rng *dest, const gsl_rng *src) {
  if (dest->type != src->type) {
    GSL_ERROR("generators must be of the same type", GSL_EINVAL);
  }

  memcpy(dest->state, src->state, src->type->size);

  return GSL_SUCCESS;
}

gsl_rng *gsl_rng_clone(const gsl_rng *q) {
  gsl_rng *r = (gsl_rng *)malloc(sizeof(gsl_rng));

  if (r == 0) {
    GSL_ERROR_VAL("failed to allocate space for rng struct", GSL_ENOMEM, 0);
  };

  r->state = malloc(q->type->size);

  if (r->state == 0) {
    free(r); /* exception in constructor, avoid memory leak */

    GSL_ERROR_VAL("failed to allocate space for rng state", GSL_ENOMEM, 0);
  };

  r->type = q->type;

  memcpy(r->state, q->state, q->type->size);

  return r;
}

void gsl_rng_set(const gsl_rng *r, unsigned long int seed) {
  (r->type->set)(r->state, seed);
}

#ifndef HIDE_INLINE_STATIC
unsigned long int gsl_rng_get(const gsl_rng *r) {
  return (r->type->get)(r->state);
}

double gsl_rng_uniform(const gsl_rng *r) {
  return (r->type->get_double)(r->state);
}

double gsl_rng_uniform_pos(const gsl_rng *r) {
  double x;
  do {
    x = (r->type->get_double)(r->state);
  } while (x == 0);

  return x;
}

unsigned long int gsl_rng_uniform_int(const gsl_rng *r, unsigned long int n) {
  unsigned long int offset = r->type->min;
  unsigned long int range = r->type->max - offset;
  unsigned long int scale = range / n;
  unsigned long int k;

  if (n > range) {
    GSL_ERROR_VAL("n exceeds maximum value of generator", GSL_EINVAL, 0);
  }

  do {
    k = (((r->type->get)(r->state)) - offset) / scale;
  } while (k >= n);

  return k;
}
#endif

unsigned long int gsl_rng_max(const gsl_rng *r) { return r->type->max; }

unsigned long int gsl_rng_min(const gsl_rng *r) { return r->type->min; }

const char *gsl_rng_name(const gsl_rng *r) { return r->type->name; }

size_t gsl_rng_size(const gsl_rng *r) { return r->type->size; }

void *gsl_rng_state(const gsl_rng *r) { return r->state; }

void gsl_rng_print_state(const gsl_rng *r) {
  size_t i;
  unsigned char *p = (unsigned char *)(r->state);
  const size_t n = r->type->size;

  for (i = 0; i < n; i++) {
    /* FIXME: we're assuming that a char is 8 bits */
    printf("%.2x", *(p + i));
  }
}

void gsl_rng_free(gsl_rng *r) {
  free(r->state);
  free(r);
}
