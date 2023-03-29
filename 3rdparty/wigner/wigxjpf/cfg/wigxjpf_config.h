
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

#ifndef __WIGXJPF_CONFIG_H__
#define __WIGXJPF_CONFIG_H__

#include "wigxjpf_auto_config.h"

/* Generate code for 'long double' results. */

#ifndef WIGXJPF_IMPL_LONG_DOUBLE    /* normally from wigxjpf_auto_config.h */ 
#define WIGXJPF_IMPL_LONG_DOUBLE    0
#endif

/* Generate code for '__float128' results. */

#ifndef WIGXJPF_IMPL_FLOAT128       /* normally from wigxjpf_auto_config.h */
#define WIGXJPF_IMPL_FLOAT128       0
#endif

/* Generate code to include semi-factorials. */

#ifndef WIGXJPF_IMPL_DOUBLE_FACTORIAL
#define WIGXJPF_IMPL_DOUBLE_FACTORIAL  0
#endif

/* Size in bytes of each prime exponent.  2 or 4. */

#define PRIME_LIST_SIZEOF_ITEM      4

/* Use vector types / instructions, e.g. SSE / AVX, for the exponents. */

/* Using vector instructions give ~0 (or at most single-digit %)
 * improvements for large symbols, and some 20 % for small symbols.
 */

#define PRIME_LIST_USE_VECTOR       0

/* Size in bytes of the vector instances.  SSE: 16, AVX:32 */

/* Using vector_size(32), i.e. AVX instructions on x86-64 has been
 * tested, and gives ~0 difference.  We are memory/cache bandwidth
 * limited.
 */

#define PRIME_LIST_VECT_SIZE        16

/* Size in bytes of each word in multi-word integers.  4 or 8. */

/* We can efficiently only use words that have half the size of the
 * maximum multiplication product the machine can easily produce.
 *
 * For 64-bit x86_64 CPUs, there is 128-bit (16 byte) result multiply
 * instructions.  However, code generation of compilers does not
 * handle the carry handling nicely, such that performance currently
 * equals that of using 32-bit words.
 */

#ifndef MULTI_WORD_INT_SIZEOF_ITEM  /* normally from wigxjpf_auto_config.h */
#define MULTI_WORD_INT_SIZEOF_ITEM  4
#endif

/* Size in bytes of word used to multiply from first prime exponents.
 * Only used when item is smaller.
 */

#ifndef MULTI_WORD_INT_SIZEOF_MULW  /* normally from wigxjpf_auto_config.h */
#define MULTI_WORD_INT_SIZEOF_MULW  8
#endif

/* Array size for precalculated factorial used in fpsimple routines. */

#define FPSIMPLE_MAX_FACTORIAL      2500

/* Debug printing. */

#define DEBUG_PRINT                 0

/* Accounting to find maximum factorial and iteration count. */

#define ACCOUNT_MAX_FACT_ITER       0

/* Special rules for MSVC. */

#ifdef _MSC_VER

# define __thread         __declspec(thread)

# include <BaseTsd.h>
typedef SSIZE_T ssize_t;

#define WIGXJPF_NOINLINE  __declspec(noinline)

#endif

/* Avoid inlining of function. */

#ifndef WIGXJPF_NOINLINE
#define WIGXJPF_NOINLINE  __attribute__((noinline))
#endif

#endif/*__WIGXJPF_CONFIG_H__*/
