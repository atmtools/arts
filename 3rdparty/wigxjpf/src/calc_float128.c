
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

#include "multi_word_int_float128.h"

#include "calc.h"

#include "wigxjpf.h"

#if WIGXJPF_IMPL_FLOAT128
# define DOUBLE_TYPE         __float128
# define DBL_FCN_POSTFIX(x)  x##_float128
# define DBL_MATH_FCN_LQ(x)  x##q
# include "calc_dbl.h"
# undef DBL_MATH_FCN_LQ
# undef DBL_FCN_POSTFIX
# undef DOUBLE_TYPE
#endif
