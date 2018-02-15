
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

/* TO BE REMOVED FROM TAR-BALL. */

#include "fastwigxj_cc.hh"

wigner369j_hash::wigner369j_hash()
{
}

wigner369j_hash::~wigner369j_hash()
{
  deinit();
#if DUMP_SYMBOLS
  flush_dump();
#endif
}

void wigner369j_hash::deinit()
{
}

void wigner369j_hash::init(const char *filename, int type,
			   fastwigxj_header *header)
{
  if (type == -1)
    {
      struct wigner369j_table dummy;
      
      wig369j_ht_clear(&dummy);
      wig369j_ht_init(&dummy, filename, type, -1, header);
      wig369j_ht_deinit(&dummy);
      return;
    }
  
  fastwigxj_load(filename, type, header);
}
