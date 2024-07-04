
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

typedef void (*wigxjpf_error_handler_func_t)(void);

extern wigxjpf_error_handler_func_t wigxjpf_error_handler;

void wigxjpf_error(void);

#if PYWIGXJPF_ERROR_HANDLING

#include <setjmp.h>

void pywigxjpf_error_handler(void);

extern __thread jmp_buf pywigxjpf_jmp_env;

/* This is called before every evaluation in the python code, in order
 * to deal gracefully with invalid user input.
 *
 * It is not recommended to do something similar in other code - rather
 * fix the code to initialise the temporary arrays large enough for all
 * used symbols!
 *
 * Note: wigxjpf_error_handler is a not thread local global, which
 * does not matter since it will be set to the same every time.
 */
# define PYWIGXJPF_ERROR_SETUP(x) do {			\
    wigxjpf_error_handler = pywigxjpf_error_handler;	\
    if (setjmp(pywigxjpf_jmp_env))			\
      return x;						\
  } while (0)
#else
# define PYWIGXJPF_ERROR_SETUP(x) do { } while (0)
#endif

# define PYWIGXJPF_ERROR_SETUP_void  PYWIGXJPF_ERROR_SETUP()
# define PYWIGXJPF_ERROR_SETUP_NaN   PYWIGXJPF_ERROR_SETUP(strtof("NAN",NULL))

void wigxjpf_drop_temp(void);

#endif/*__WIGXJPF_ERROR_H__*/
