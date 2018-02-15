
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

#ifndef __GEN_BUF_HH__
#define __GEN_BUF_HH__

#include <stdint.h>
#include <stdlib.h>

extern uint64_t buf[1024];
extern size_t nbuf;

void write_buf();

inline void insert_buf(uint64_t key)
{
  buf[nbuf++] = key;
  if (nbuf >= sizeof(buf)/sizeof(buf[0]))
    write_buf();
}

#endif//__GEN_BUF_HH__
