
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

#ifndef __WIGNER369J_HASH_FCN_H__
#define __WIGNER369J_HASH_FCN_H__

/* All interfaces in this header are internal.
 * Do not use directly.
 */

#define WIGNER369_HASH_FCN_XOR_XOR

#ifdef WIGNER369_HASH_FCN_XOR_XOR
#define WIGNER369_HASH_FCN_PRM_XOR1A  15
#define WIGNER369_HASH_FCN_PRM_XOR1B  5
#define WIGNER369_HASH_FCN_PRM_XOR1C  17
#define WIGNER369_HASH_FCN_PRM_XOR2A  21
#define WIGNER369_HASH_FCN_PRM_XOR2B  17
#define WIGNER369_HASH_FCN_PRM_XOR2C  24
#define WIGNER369_HASH_FCN(key_in)					\
  ({									\
    uint64_t x_out = (key_in);						\
    /* XOR1 */								\
    x_out = x_out ^ (x_out >> (WIGNER369_HASH_FCN_PRM_XOR1C));		\
    x_out = x_out ^ (x_out << (WIGNER369_HASH_FCN_PRM_XOR1B));		\
    x_out = x_out ^ (x_out >> (WIGNER369_HASH_FCN_PRM_XOR1A));		\
    /* XOR1 */								\
    x_out = x_out ^ (x_out >> (WIGNER369_HASH_FCN_PRM_XOR2A));		\
    x_out = x_out ^ (x_out << (WIGNER369_HASH_FCN_PRM_XOR2B));		\
    x_out = x_out ^ (x_out >> (WIGNER369_HASH_FCN_PRM_XOR2C));		\
    /* done */								\
    x_out;								\
  })
#define WIGNER369_HASH_FCN_DESCR "XOR_RLR_XOR_RLR"
#define WIGNER369_HASH_FCN_PARAMS {				\
    (WIGNER369_HASH_FCN_PRM_XOR1A), (WIGNER369_HASH_FCN_PRM_XOR1B), \
    (WIGNER369_HASH_FCN_PRM_XOR1C),				\
    (WIGNER369_HASH_FCN_PRM_XOR2A), (WIGNER369_HASH_FCN_PRM_XOR2B), \
    (WIGNER369_HASH_FCN_PRM_XOR2C),				\
    0, 0, 0, 0						        \
  }
#endif

#endif/*__WIGNER369J_HASH_FCN_H__*/
