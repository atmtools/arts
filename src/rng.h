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

/**
 * @file   rng.h
 * @author Cory Davis <cory@met.ed.ac.uk>
 * @date   2003-06-26
 *
 * @brief  Defines the Rng random number generator class
 *
 * The Rng class is described at the very end of this file.  The rest of the file,
 * which describes the code that actually does the work, was obtained from the
 * GNU Scientific Library <http://www.gnu.org/software/gsl/>.
 *
 * The Rng class uses the gsl_rng_mt_19937 random number generator whose original
 * implementation was copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
 * Coded by Takuji Nishimura, considering the suggestions by Topher Cooper and
 * Marc Rieffel in July-Aug. 1997, "A C-program for MT19937: Integer version
 * (1998/4/6)".
 *
 * The period of this generator is 2^{19937} - 1.
 */

/* gsl_types.h
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
#ifndef rng_h
#define rng_h

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#ifdef _DLL
#ifdef DLL_EXPORT
#define GSL_VAR __declspec(dllexport)
#else
#define GSL_VAR __declspec(dllimport)
#endif
#else
#define GSL_VAR extern
#endif
#else
#define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* rng/gsl_rng.h
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

#ifndef __GSL_RNG_H__
#define __GSL_RNG_H__
#include <cstdlib>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS /* empty */
#define __END_DECLS   /* empty */
#endif

__BEGIN_DECLS

typedef struct {
  const char *name;
  unsigned long int max;
  unsigned long int min;
  size_t size;
  void (*set)(void *state, unsigned long int seed);
  unsigned long int (*get)(void *state);
  double (*get_double)(void *state);
} gsl_rng_type;

typedef struct {
  const gsl_rng_type *type;
  void *state;
} gsl_rng;

/* These structs also need to appear in default.c so you can select
   them via the environment variable GSL_RNG_TYPE */

GSL_VAR const gsl_rng_type *gsl_rng_borosh13;
GSL_VAR const gsl_rng_type *gsl_rng_coveyou;
GSL_VAR const gsl_rng_type *gsl_rng_cmrg;
GSL_VAR const gsl_rng_type *gsl_rng_fishman18;
GSL_VAR const gsl_rng_type *gsl_rng_fishman20;
GSL_VAR const gsl_rng_type *gsl_rng_fishman2x;
GSL_VAR const gsl_rng_type *gsl_rng_gfsr4;
GSL_VAR const gsl_rng_type *gsl_rng_knuthran;
GSL_VAR const gsl_rng_type *gsl_rng_knuthran2;
GSL_VAR const gsl_rng_type *gsl_rng_lecuyer21;
GSL_VAR const gsl_rng_type *gsl_rng_minstd;
GSL_VAR const gsl_rng_type *gsl_rng_mrg;
GSL_VAR const gsl_rng_type *gsl_rng_mt19937;
GSL_VAR const gsl_rng_type *gsl_rng_mt19937_1999;
GSL_VAR const gsl_rng_type *gsl_rng_mt19937_1998;
GSL_VAR const gsl_rng_type *gsl_rng_r250;
GSL_VAR const gsl_rng_type *gsl_rng_ran0;
GSL_VAR const gsl_rng_type *gsl_rng_ran1;
GSL_VAR const gsl_rng_type *gsl_rng_ran2;
GSL_VAR const gsl_rng_type *gsl_rng_ran3;
GSL_VAR const gsl_rng_type *gsl_rng_rand;
GSL_VAR const gsl_rng_type *gsl_rng_rand48;
GSL_VAR const gsl_rng_type *gsl_rng_random128_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random128_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random128_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random256_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random256_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random256_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random32_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random32_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random32_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random64_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random64_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random64_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random8_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random8_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random8_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_random_bsd;
GSL_VAR const gsl_rng_type *gsl_rng_random_glibc2;
GSL_VAR const gsl_rng_type *gsl_rng_random_libc5;
GSL_VAR const gsl_rng_type *gsl_rng_randu;
GSL_VAR const gsl_rng_type *gsl_rng_ranf;
GSL_VAR const gsl_rng_type *gsl_rng_ranlux;
GSL_VAR const gsl_rng_type *gsl_rng_ranlux389;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxd1;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxd2;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxs0;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxs1;
GSL_VAR const gsl_rng_type *gsl_rng_ranlxs2;
GSL_VAR const gsl_rng_type *gsl_rng_ranmar;
GSL_VAR const gsl_rng_type *gsl_rng_slatec;
GSL_VAR const gsl_rng_type *gsl_rng_taus;
GSL_VAR const gsl_rng_type *gsl_rng_taus2;
GSL_VAR const gsl_rng_type *gsl_rng_taus113;
GSL_VAR const gsl_rng_type *gsl_rng_transputer;
GSL_VAR const gsl_rng_type *gsl_rng_tt800;
GSL_VAR const gsl_rng_type *gsl_rng_uni;
GSL_VAR const gsl_rng_type *gsl_rng_uni32;
GSL_VAR const gsl_rng_type *gsl_rng_vax;
GSL_VAR const gsl_rng_type *gsl_rng_waterman14;
GSL_VAR const gsl_rng_type *gsl_rng_zuf;

const gsl_rng_type **gsl_rng_types_setup(void);

GSL_VAR const gsl_rng_type *gsl_rng_default;
GSL_VAR unsigned long int gsl_rng_default_seed;

gsl_rng *gsl_rng_alloc(const gsl_rng_type *T);
int gsl_rng_memcpy(gsl_rng *dest, const gsl_rng *src);
gsl_rng *gsl_rng_clone(const gsl_rng *r);

void gsl_rng_free(gsl_rng *r);

void gsl_rng_set(const gsl_rng *r, unsigned long int seed);
unsigned long int gsl_rng_max(const gsl_rng *r);
unsigned long int gsl_rng_min(const gsl_rng *r);
const char *gsl_rng_name(const gsl_rng *r);
size_t gsl_rng_size(const gsl_rng *r);
void *gsl_rng_state(const gsl_rng *r);

void gsl_rng_print_state(const gsl_rng *r);

const gsl_rng_type *gsl_rng_env_setup(void);

unsigned long int gsl_rng_get(const gsl_rng *r);
double gsl_rng_uniform(const gsl_rng *r);
double gsl_rng_uniform_pos(const gsl_rng *r);
unsigned long int gsl_rng_uniform_int(const gsl_rng *r, unsigned long int n);

#ifdef HAVE_INLINE
extern inline unsigned long int gsl_rng_get(const gsl_rng *r);

extern inline unsigned long int gsl_rng_get(const gsl_rng *r) {
  return (r->type->get)(r->state);
}

extern inline double gsl_rng_uniform(const gsl_rng *r);

extern inline double gsl_rng_uniform(const gsl_rng *r) {
  return (r->type->get_double)(r->state);
}

extern inline double gsl_rng_uniform_pos(const gsl_rng *r);

extern inline double gsl_rng_uniform_pos(const gsl_rng *r) {
  double x;
  do {
    x = (r->type->get_double)(r->state);
  } while (x == 0);

  return x;
}

extern inline unsigned long int gsl_rng_uniform_int(const gsl_rng *r,
                                                    unsigned long int n);

extern inline unsigned long int gsl_rng_uniform_int(const gsl_rng *r,
                                                    unsigned long int n) {
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
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_RNG_H__ */

/* err/gsl_errno.h
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

#ifndef __GSL_ERRNO_H__
#define __GSL_ERRNO_H__

#include <cerrno>
#include <cstdio>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS /* empty */
#define __END_DECLS   /* empty */
#endif

__BEGIN_DECLS

enum {
  GSL_SUCCESS = 0,
  GSL_FAILURE = -1,
  GSL_CONTINUE = -2, /* iteration has not converged */
  GSL_EDOM = 1,      /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE = 2,    /* output range error, e.g. exp(1e100) */
  GSL_EFAULT = 3,    /* invalid pointer */
  GSL_EINVAL = 4,    /* invalid argument supplied by user */
  GSL_EFAILED = 5,   /* generic failure */
  GSL_EFACTOR = 6,   /* factorization failed */
  GSL_ESANITY = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM = 8,    /* malloc failed */
  GSL_EBADFUNC = 9,  /* problem with user-supplied function */
  GSL_ERUNAWAY = 10, /* iterative process is out of control */
  GSL_EMAXITER = 11, /* exceeded max number of iterations */
  GSL_EZERODIV = 12, /* tried to divide by zero */
  GSL_EBADTOL = 13,  /* user specified an invalid tolerance */
  GSL_ETOL = 14,     /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15, /* underflow */
  GSL_EOVRFLW = 16,  /* overflow  */
  GSL_ELOSS = 17,    /* loss of accuracy */
  GSL_EROUND = 18,   /* failed because of roundoff error */
  GSL_EBADLEN = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR = 20,  /* matrix not square */
  GSL_ESING = 21,    /* apparent singularity detected */
  GSL_EDIVERGE = 22, /* integral or series is divergent */
  GSL_EUNSUP = 23,   /* requested feature is not supported by the hardware */
  GSL_EUNIMPL = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE = 25,   /* cache limit exceeded */
  GSL_ETABLE = 26,   /* table limit exceeded */
  GSL_ENOPROG = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28, /* jacobian evaluations are not improving the solution */
  GSL_ETOLF = 29,    /* cannot reach the specified tolerance in F */
  GSL_ETOLX = 30,    /* cannot reach the specified tolerance in X */
  GSL_ETOLG = 31,    /* cannot reach the specified tolerance in gradient */
  GSL_EOF = 32       /* end of file */
};

void gsl_error(const char *reason, const char *file, int line, int gsl_errno);

void gsl_warning(const char *reason, const char *file, int line, int gsl_errno);

void gsl_stream_printf(const char *label,
                       const char *file,
                       int line,
                       const char *reason);

const char *gsl_strerror(const int gsl_errno);

typedef void gsl_error_handler_t(const char *reason,
                                 const char *file,
                                 int line,
                                 int gsl_errno);

typedef void gsl_stream_handler_t(const char *label,
                                  const char *file,
                                  int line,
                                  const char *reason);

gsl_error_handler_t *gsl_set_error_handler(gsl_error_handler_t *new_handler);

gsl_error_handler_t *gsl_set_error_handler_off(void);

gsl_stream_handler_t *gsl_set_stream_handler(gsl_stream_handler_t *new_handler);

FILE *gsl_set_stream(FILE *new_stream);

/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno)                  \
  do {                                                \
    gsl_error(reason, __FILE__, __LINE__, gsl_errno); \
    return gsl_errno;                                 \
  } while (0)

/* GSL_ERROR_VAL: call the error handler, and return the given value */

#define GSL_ERROR_VAL(reason, gsl_errno, value)       \
  do {                                                \
    gsl_error(reason, __FILE__, __LINE__, gsl_errno); \
    return value;                                     \
  } while (0)

/* GSL_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define GSL_ERROR_VOID(reason, gsl_errno)             \
  do {                                                \
    gsl_error(reason, __FILE__, __LINE__, gsl_errno); \
    return;                                           \
  } while (0)

/* GSL_ERROR_NULL suitable for out-of-memory conditions */

#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)

/* GSL library code can occasionally generate warnings, which are not
   intended to be fatal. You can compile a version of the library with
   warnings turned off globally by defining the preprocessor constant
   GSL_WARNINGS_OFF. This turns off the warnings, but does not disable
   error handling in any way or turn off error messages.
 
   GSL_WARNING() is not intended for use in client code -- use
   GSL_MESSAGE() instead.  */

#ifdef GSL_WARNINGS_OFF /* throw away warnings */
#define GSL_WARNING(warning, gsl_errno) \
  do {                                  \
  } while (0)
#else /* output all warnings */
#define GSL_WARNING(warning, gsl_errno)                  \
  do {                                                   \
    gsl_warning(warning, __FILE__, __LINE__, gsl_errno); \
  } while (0)
#endif

/* Warnings can also be turned off at runtime by setting the variable
   gsl_warnings_off to a non-zero value */

GSL_VAR int gsl_warnings_off;

/* Sometimes you have several status results returned from
 * function calls and you want to combine them in some sensible
 * way. You cannot produce a "total" status condition, but you can
 * pick one from a set of conditions based on an implied hierarchy.
 *
 * In other words:
 *    you have: status_a, status_b, ...
 *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
 *
 * In this example you consider status_a to be more important and
 * it is checked first, followed by the others in the order specified.
 *
 * Here are some dumb macros to do this.
 */
#define GSL_ERROR_SELECT_2(a, b) \
  ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a, b, c) \
  ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b, c))
#define GSL_ERROR_SELECT_4(a, b, c, d) \
  ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b, c, d))
#define GSL_ERROR_SELECT_5(a, b, c, d, e) \
  ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b, c, d, e))

#define GSL_STATUS_UPDATE(sp, s)         \
  do {                                   \
    if ((s) != GSL_SUCCESS) *(sp) = (s); \
  } while (0)

__END_DECLS

#endif /* __GSL_ERRNO_H__ */

/* err/gsl_message.h
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

#ifndef __GSL_MESSAGE_H__
#define __GSL_MESSAGE_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS /* empty */
#define __END_DECLS   /* empty */
#endif

__BEGIN_DECLS

/* Provide a general messaging service for client use.  Messages can
 * be selectively turned off at compile time by defining an
 * appropriate message mask. Client code which uses the GSL_MESSAGE()
 * macro must provide a mask which is or'ed with the GSL_MESSAGE_MASK.
 *
 * The messaging service can be completely turned off
 * by defining GSL_MESSAGING_OFF.  */

void gsl_message(const char *message,
                 const char *file,
                 int line,
                 unsigned int mask);

#ifndef GSL_MESSAGE_MASK
#define GSL_MESSAGE_MASK 0xffffffffu /* default all messages allowed */
#endif

GSL_VAR unsigned int gsl_message_mask;

/* Provide some symolic masks for client ease of use. */

enum {
  GSL_MESSAGE_MASK_A = 1,
  GSL_MESSAGE_MASK_B = 2,
  GSL_MESSAGE_MASK_C = 4,
  GSL_MESSAGE_MASK_D = 8,
  GSL_MESSAGE_MASK_E = 16,
  GSL_MESSAGE_MASK_F = 32,
  GSL_MESSAGE_MASK_G = 64,
  GSL_MESSAGE_MASK_H = 128
};

#ifdef GSL_MESSAGING_OFF /* throw away messages */
#define GSL_MESSAGE(message, mask) \
  do {                             \
  } while (0)
#else /* output all messages */
#define GSL_MESSAGE(message, mask)                    \
  do {                                                \
    if (mask & GSL_MESSAGE_MASK)                      \
      gsl_message(message, __FILE__, __LINE__, mask); \
  } while (0)
#endif

__END_DECLS

#endif /* __GSL_MESSAGE_H__ */

///////////////////////////////////////////////////////////////

/*CPD: 26-06-02. Here is my contribution to this file: a simple 
random number generator class*/

#include <ctime>

class Rng {
  gsl_rng *r;  //the gsl random number generator instance

  unsigned long int seed_no;  //The integer used to seen the Rng

 public:
  Rng();  //constructor

  ~Rng();  //destructor

 /**
  * Seeds the Rng with the integer argument.
  *
  * Every seed is only used once. The provided seed is increased by 1 until an
  * unused seed is found.
  *
  * The default is to seed the Rng using the seconds elapsed since 1970.
  */
  void seed(unsigned long int n);

  void force_seed(unsigned long int n);

  double draw();  //draw a random number between [0,1)

  unsigned long int showseed() const;  //return the seed.
};

#endif /* rng_h  */
