
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

#include <stdio.h>
#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "gen_header.hh"

#include <algorithm>

// We read data continously from stdin, in chunks of a few ten
// thousand elements we sort them, and inject them by merge-sort into
// the larger buffer.  Once the larger buffer is full, this is ejected
// to the next program for insertation into the big list.  By ejecting
// to a next program, we get the file-injection async with the input
// of this program for free.

#define LARGE_BUF_ENTRIES (128*1024*1024)

uint64_t *_buf = NULL;
size_t  _nbuf = 0;

template<typename T>
bool compare (T i,T j)
{
  return i < j;
}

#include "sort_deduplicate.h"

uint64_t ejected = 0;

void eject()
{
  ejected += _nbuf;
  
  if (fwrite (_buf, sizeof(_buf[0]), _nbuf, stdout) != _nbuf)
    {
      fprintf (stderr, "Error writing sorted chunk.\n");
      exit(1);
    }

  _nbuf = 0; 
}

int main()
{
  uint64_t got = 0;

  size_t sz = LARGE_BUF_ENTRIES * sizeof (_buf[0]);
  _buf = (uint64_t *) malloc (sz);
  if (!_buf)
    {
      fprintf (stderr, "Error allocating sort array (%zd MB).\n",
	       sz / 1000000);
      exit(1);
    }

  fastwigxj_header header;

#if 0 /* when combining from a raw 64-bit source */
  prepare_header(&header);
#else
  read_header(&header);
#endif
  write_header(&header);

  while (!feof(stdin))
    {
      size_t entries = LARGE_BUF_ENTRIES;
      size_t left = entries - _nbuf;

      if (!left)
	{
	  sort_deduplicate();
	  
	  if (left < entries / 10)
	    eject();
	}

      size_t n = fread(_buf + _nbuf,
		       sizeof(_buf[0]), left, 
		       stdin);
      
      if (!n)
	{
	  if (ferror(stdin))
	    {
	      fprintf (stderr, "Error reading chunk to sort.\n");
	      exit(1);
	    }
	  continue;
	}
      
      got += n;

      _nbuf += n;
   }

  sort_deduplicate();
  eject();
  
  fprintf (stderr,"Read symbols: %" PRIu64 ".\n",got);
  fprintf (stderr,"Wrote symbols: %" PRIu64 ".\n",ejected);

  return 0;
}
