
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

#ifndef __WIGXJPF_C_WRAP_H__
#define __WIGXJPF_C_WRAP_H__

#include "wigxjpf_auto_config.h"

#include "wigxjpf.h"

#if WIGXJPF_HAVE_THREAD
extern __thread struct wigxjpf_temp *wigxjpf_global_temp;
#else
extern struct wigxjpf_temp *wigxjpf_global_temp;
#endif

#endif/*__WIGXJPF_C_WRAP_H__*/
