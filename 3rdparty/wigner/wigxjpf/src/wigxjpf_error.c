
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

#include "wigxjpf_error.h"

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>

__thread jmp_buf _error_jmp_env;


void wigxjpf_error_pywigxjpf(void)
{
  /* Allow reuse of the temp array. */
  wigxjpf_drop_temp();

  fprintf(stderr,
	  "\n"
	  "pywigxjpf: Error detected! "
	  "** Library misuse?  See documentation. **\n"
	  "\n");

  longjmp(_error_jmp_env, 1);
}


void wigxjpf_error_exit(void)
{
  fprintf (stderr, "wigxjpf: Abort.\n");
  exit(1);
}

void wigxjpf_error_set_errno(void)
{
  wigxjpf_drop_temp();
  errno = EDOM;
  longjmp(_error_jmp_env, 1);
}

void wigxjpf_error(void) {
  #ifdef PYWIGXJPF_ERROR_HANDLING
  wigxjpf_error_pywigxjpf();
  #else
  #if ERRNO_ERROR_HANDLING
  wigxjpf_error_set_errno();
  #else
  wigxjpf_error_exit();
  #endif
  #endif
}
