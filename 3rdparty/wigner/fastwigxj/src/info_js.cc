
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

#include "fastwigxj_cc.hh"
#include "fastwigxj_header.h"

#include <string.h>
#include <math.h>

void usage(const char *cmd)
{
  printf ("Usage:  %s hashfile\n",cmd);
}

int main(int argc, char *argv[])
{
  const char *hashfile = NULL;

  for (int i = 1; i < argc; i++)
    {
      if (!hashfile)
	hashfile = argv[i];
      else
	{
	  fprintf (stderr,"Bad argument: %s\n", argv[i]);
	  
	  usage(argv[0]);
	  exit(1);
	}
   }

  if (!hashfile)
    {
      usage(argv[0]);
      exit(1);
    }

  wigner369j_hash hash;
  fastwigxj_header header;

  hash.init(hashfile, -1, &header);

  printf ("Hash file:     %s\n",  hashfile);
  printf ("Version:       %d\n",  header._version);
  printf ("Entries:       %" PRIu64 " (%.1f %%)\n", header._entries,
	  (100.0 * (double) header._entries) / (double) header._hashentries);
  printf ("Type:          %dj\n", header._type);
  printf ("c14n:          %d\n", header._c14n);
  printf ("Hashentries:   %" PRIu64 "\n", header._hashentries);
  printf ("Hash function: %s\n", header._hash_func_descr);
  printf ("Hash func prm:");
  for (size_t i = 0; i < (sizeof (header._hash_param) /
		       sizeof (header._hash_param[0])); i++)
    printf (" %" PRIu64 "", header._hash_param[i]);
  printf ("\n");    
  printf ("Max 2j:        %d\n",  header._max_two_j);
  printf ("Max abs err:   %.5e\n",header._max_abs_err);
  printf ("Max rel err:   %.5e (|val| >  1e-10)\n",
	  header._max_rel_err_gt_1em10);
  printf ("Max small val: %.5e (|val| <= 1e-10)\n",header._max_small_val);
  printf ("Checksums:     0x%016" PRIx64 ", 0x%016" PRIx64 "\n",
	  header._checksum_xor, header._checksum_xor_header);

  hash.deinit();

  return 0;
}
