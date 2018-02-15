
/* Copyright 2015 Haakan T. Johansson */

/*  This file is part of FASTWIGXJ.
 *
 *  FASTWIGXJ is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  FASTWIGXJ is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with FASTWIGXJ.  If not, see
 *  <http://www.gnu.org/licenses/>.
 */

#ifndef __WIGXJPF_CONFIG_H__
#define __WIGXJPF_CONFIG_H__

#include "fastwigxj_auto_config.h"

/* Use WIGXJPF 'long double' code. */

#ifndef FASTWIGXJ_USE_LONG_DOUBLE  /* normally from fastwigxj_auto_config.h */ 
#define FASTWIGXJ_USE_LONG_DOUBLE  0
#endif

/* Use WIGXJPF '__float128' code for 9j fallback to 6j. */

#ifndef FASTWIGXJ_USE_FLOAT128     /* normally from fastwigxj_auto_config.h */
#define FASTWIGXJ_USE_FLOAT128     0
#endif

#ifndef FASTWIGXJ_HAVE_THREAD      /* normally from fastwigxj_auto_config.h */ 
#define FASTWIGXJ_HAVE_THREAD      0
#endif

#ifndef FASTWIGXJ_HAVE_SSE4_1      /* normally from fastwigxj_auto_config.h */ 
#define FASTWIGXJ_HAVE_SSE4_1      0
#endif

#ifndef FASTWIGXJ_HAVE_AVX2        /* normally from fastwigxj_auto_config.h */ 
#define FASTWIGXJ_HAVE_AVX2        0
#endif

#endif/*__FASTWIGXJ_CONFIG_H__*/
