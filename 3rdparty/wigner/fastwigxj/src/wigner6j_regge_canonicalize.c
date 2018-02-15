
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

#include <stdint.h>
#include <stdio.h>
#include "canonicalise_inc.h"
#include "fastwigxj_struct.h"

#include "wigner36j_regge_canonicalize.h"

#define WIGNER6J_REGGE_CANONICALISE_NAME wigner6j_regge_canonicalise_index
#define WIGNER6J_REGGE_CANONICALISE_RETURN_INDEX     1
#define WIGNER6J_REGGE_CANONICALISE_TRIVIAL_0_CHECK  1
#define WIGNER6J_REGGE_CANONICALISE_IN_STRUCT        0
#define WIGNER6J_REGGE_CANONICALISE_ARRAY4           0

#include "wigner6j_regge_canonicalize_inc.h"

#undef WIGNER6J_REGGE_CANONICALISE_NAME
#undef WIGNER6J_REGGE_CANONICALISE_RETURN_INDEX
#undef WIGNER6J_REGGE_CANONICALISE_TRIVIAL_0_CHECK
#undef WIGNER6J_REGGE_CANONICALISE_IN_STRUCT
#undef WIGNER6J_REGGE_CANONICALISE_ARRAY4

#define WIGNER6J_REGGE_CANONICALISE_NAME wigner6j_struct_array4_regge_canonicalise_index_no_t0_chk
#define WIGNER6J_REGGE_CANONICALISE_RETURN_INDEX     1
#define WIGNER6J_REGGE_CANONICALISE_TRIVIAL_0_CHECK  0
#define WIGNER6J_REGGE_CANONICALISE_IN_STRUCT        1
#define WIGNER6J_REGGE_CANONICALISE_ARRAY4           1

#include "wigner6j_regge_canonicalize_inc.h"

#undef WIGNER6J_REGGE_CANONICALISE_NAME
#undef WIGNER6J_REGGE_CANONICALISE_RETURN_INDEX
#undef WIGNER6J_REGGE_CANONICALISE_TRIVIAL_0_CHECK
#undef WIGNER6J_REGGE_CANONICALISE_IN_STRUCT
#undef WIGNER6J_REGGE_CANONICALISE_ARRAY4
