
/* Copyright 2019 Haakan T. Johansson */

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

#ifndef __WIGXJPF_ERROR_H__
#define __WIGXJPF_ERROR_H__

#include "wigxjpf_config.h"

void wigxjpf_error(void);

#include <setjmp.h>

extern __thread jmp_buf _error_jmp_env;

#if PYWIGXJPF_ERROR_HANDLING || ERRNO_ERROR_HANDLING
#define NONABORT_ERROR_SETUP(x)           \
  do {                                    \
    if (setjmp(_error_jmp_env)) return x; \
  } while (0)
#else
#define NONABORT_ERROR_SETUP(x) do {} while(0)
#endif

void wigxjpf_drop_temp(void);

# define NONABORT_ERROR_SETUP_void  NONABORT_ERROR_SETUP()
# define NONABORT_ERROR_SETUP_NaN   NONABORT_ERROR_SETUP(strtof("NAN",NULL))

#endif/*__WIGXJPF_ERROR_H__*/
