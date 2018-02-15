
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
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <unistd.h>

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
size_t   _nbuf = 0;

template<typename T>
bool compare (T i,T j)
{
  return i < j;
}

#include "sort_deduplicate.h"

uint64_t ejected = 0;
uint64_t lastejected = 0;

fastwigxj_header _header;

const char *_filename = NULL;
char *_filename_tmp = NULL;
char *_filename_rd  = NULL;

void eject(FILE* fid, uint64_t *ptr, size_t len)
{
  lastejected += len;
  
  if (fwrite (ptr, sizeof(ptr[0]), len, fid) != len)
    {
      fprintf (stderr, "Error writing sorted chunk.\n");
      exit(1);
    }
}

void merge_with_file()
{
#define FBUF_SIZE (16*1024)

  uint64_t in[FBUF_SIZE];
  uint64_t out[FBUF_SIZE];
  uint64_t *pin, *pin_end;
  uint64_t *pout, *pout_end;

  pin = pin_end = in;
  pout = out;
  pout_end = out + FBUF_SIZE;

  uint64_t *p = _buf;
  uint64_t *p_end = _buf + _nbuf;

  rename(_filename_tmp, _filename_rd);

  FILE *fin  = fopen(_filename_rd,"r");
  FILE *fout = fopen(_filename_tmp,"w");

  write_header(&_header, fout);
  
  for ( ; ; )
    {
      if (pin >= pin_end)
	{
	  // We need more data from the file

	  size_t n = fread(in, sizeof(in[0]), FBUF_SIZE, fin);

	  //fprintf (stderr,"%zd\n",n);

	  if (!n)
	    {
	      if (ferror(fin))
		{
		  fprintf (stderr, "Error reading from tmp file.\n");
		  exit(1);
		}
	      if (feof(fin))
		break;
	      fprintf (stderr, "Got no entry while reading from tmp file.\n");
	      exit(1);
	    }

	  pin = in;
	  pin_end = in + n;
	}

      if (p >= p_end)
	break;

      // So, who is the smaller one?

      if (*pin < *p)
	*(pout++) = *(pin++);
      else if (*pin > *p)
	*(pout++) = *(p++);
      else // equal
	{
	  *(pout++) = *(pin++);
	  p++;
	}

      if (pout >= pout_end)
	{
	  eject(fout, out, pout - out);
	  pout = out;
	}
    }

  // Write out whatever is left in the output buffer.
  eject(fout, out, pout - out);

  // Stuff left over either from file or from buffer.
  // Just write the two of them, only one will have entries.

  eject(fout, p, p_end - p);

  eject(fout, pin, pin_end - pin);

  // We consumed all buffered data
  _nbuf = 0;

  // Is there data left in the input file?

  for ( ; ; )
    {
      size_t n = fread(in, sizeof(in[0]), FBUF_SIZE, fin);

      if (!n)
	{
	  if (ferror(fin))
	    {
	      fprintf (stderr, "Error reading from tmp file.\n");
	      exit(1);
	    }
	  if (feof(fin))
	    break;
	  fprintf (stderr, "Got no entry while reading from tmp file.\n");
	  exit(1);
	}

      eject(fout, in, n);
    }

  fclose (fin);
  fclose (fout);

  unlink(_filename_rd);
}

int main(int argc, char *argv[])
{
  uint64_t got = 0;
  uint64_t cycles = 0;

  bool first = true;

  if (argc < 2)
    {
      fprintf (stderr,"Usage:  %s outfile\n",argv[0]);
      exit(1);
    }

  _filename = argv[1];
  _filename_tmp = (char *) malloc(strlen(_filename)+strlen(".tmp")+1);
  sprintf(_filename_tmp, "%s%s",_filename,".tmp");
  _filename_rd  = (char *) malloc(strlen(_filename)+strlen(".rd")+1);
  sprintf(_filename_rd, "%s%s",_filename,".rd");

  size_t sz = LARGE_BUF_ENTRIES * sizeof (_buf[0]);
  _buf = (uint64_t *) malloc (sz);
  if (!_buf)
    {
      fprintf (stderr, "Error allocating sort array (%zd MB).\n",
	       sz / 1000000);
      exit(1);
    }

  read_header(&_header);

  while (!feof(stdin))
    {
      size_t entries = LARGE_BUF_ENTRIES;
      
      size_t n = fread(_buf,
		       sizeof(_buf[0]), entries, 
		       stdin);
      
      if (!n)
	{
	  if (ferror(stdin))
	    {
	      fprintf (stderr, "Error reading chunk to insert.\n");
	      exit(1);
	    }
	  fprintf (stderr, "Got no entry while reading chunk to insert.\n");
	  exit(1);
	}
      
      got += n;
      _nbuf += n;

      // Sort and remove duplicates

      sort_deduplicate();

      // Now, insert it into the file.

      lastejected = 0;
      if (first)
	{
	  FILE *fout = fopen(_filename_tmp,"w");
	  write_header(&_header, fout);
	  eject(fout, _buf, _nbuf);
	  fclose(fout);
	  _nbuf = 0;
	  first = false;
	}
      else
	merge_with_file();
      ejected += lastejected;
      cycles++;
   }
  
  fprintf (stderr,"Read chunks: %" PRIu64 ".\n",got);
  fprintf (stderr,"Wrote to file(s): %" PRIu64 ".\n",ejected);
  fprintf (stderr,"Last file entries: %" PRIu64 ".\n",lastejected);
  fprintf (stderr,"File cycles: %" PRIu64 ".\n",cycles);

  // All finished well, rename file to final name

  rename(_filename_tmp, _filename);

  return 0;
}
