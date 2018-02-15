
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

#include "gen_buf.hh"

#include <stdio.h>

uint64_t buf[1024];
size_t nbuf = 0;

void write_buf()
{
  size_t n = fwrite(buf, sizeof (buf[0]), nbuf, stdout);
  if (n != nbuf)
    {
      fprintf(stderr,"Error writing symbols.\n");
      exit(1);
    }
  nbuf = 0;
}
