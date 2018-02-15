
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

#include "gen_header.hh"

#include <stdio.h>
#include <string.h>

void prepare_header(fastwigxj_header *header)
{
  memset(header, 0, sizeof (*header));

  header->_magic       = FASTWIGXJ_MAGIC;
  header->_magic_end   = FASTWIGXJ_MAGIC;
  header->_version     = FASTWIGXJ_VERSION;
}

void write_header(fastwigxj_header *header, FILE *fid)
{
  size_t n = fwrite(header, sizeof (*header), 1, fid);

  if (n != 1)
    {
      fprintf(stderr,"Error writing header.\n");
      exit(1);
    }
}

void read_header(fastwigxj_header *header, FILE *fid)
{
  size_t n = fread(header, sizeof (*header), 1, fid);

  if (n != 1)
    {
      fprintf(stderr,"Error reading header.\n");
      exit(1);
    }

  if (header->_magic != FASTWIGXJ_MAGIC)
    {
      fprintf (stderr, "Bad magic in input header (%08x != %08x).\n",
               header->_magic, FASTWIGXJ_MAGIC);
      exit(1);
    }

  if (header->_version != FASTWIGXJ_VERSION)
    {
      fprintf (stderr, "Bad version in input header (%08x != %08x).\n",
               header->_version, FASTWIGXJ_MAGIC);
      exit(1);
    }
}
