
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

#include <sys/time.h>
#include <float.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "canonicalise.c"

#include "gen_buf.hh"
#include "gen_header.hh"

int _max_two_j_abde;
int _max_two_j_cfgh;
int _max_two_j_i;

int _raw = 0;
int _all = 0;
int _c14n = 1;

struct stats
{
  uint64_t count;
  uint64_t null_loop;
};

// ja jb jc
// jd je jf
// jg jh ji
// Mimimum values:
// 0  ja jb
// jb ja ja
// jb ja ja
// Minimum values like above does not work (miss some combinations
// that min values all zeros would find).  These however work
// (empirically):
// 0  ja 0
// jb ja 0
// 0  0  0

// Generation order:
// ja, jb, jc, jd  (exhaustive)
// jc, jf, jg, jh  (according to traingular conditions)
// ji              (according to traingular conditions)

#define START_MIN(jmin)      (_all ? 0 : (jmin))
#define START_MIN_NONE(jmin) 0

void calc_min_max(int j1, int j2, int &j_min, int &j_max,
		  int j_start, int j_end)
{
  j_min = abs(j1 - j2);
  j_max = j1 + j2;
  
  if (j_min < j_start)
    j_min = j_start;
  if ((j_min ^ j_max) & 1) // j_max still has correct odd/even
    j_min++;
  if (j_max > j_end)
    j_max = j_end;
}

int trivial_zero_9j(int ja, int jb, int jc,
		    int jd, int je, int jf,
		    int jg, int jh, int ji)
{
#define OUTSIDE_TRIANG(j1,j2,j3) (j3 < abs(j1 - j2) || j3 > (j1 + j2))

  if (OUTSIDE_TRIANG(ja, jb, jc) ||
      OUTSIDE_TRIANG(jd, je, jf) ||
      OUTSIDE_TRIANG(jg, jh, ji) ||
      OUTSIDE_TRIANG(ja, jd, jg) ||
      OUTSIDE_TRIANG(jb, je, jh) ||
      OUTSIDE_TRIANG(jc, jf, ji))
    return 1;

  if (((ja ^ jb ^ jc) |
       (jd ^ je ^ jf) |
       (jg ^ jh ^ ji) |
       (ja ^ jd ^ jg) |
       (jb ^ je ^ jh) |
       (jc ^ jf ^ ji)) & 1)
    return 1;

  return 0;
}

void loop_four(stats &st, 
	       int ja, int jb, int jd, int je,
	       int jc_min, int jc_max)
{
  int jf_min, jf_max;

  calc_min_max(jd, je, jf_min, jf_max, START_MIN_NONE(ja), _max_two_j_cfgh);

  if (jf_max < jf_min) // no cutting?
    {
      st.null_loop++;
      return;
    }

  int jg_min, jg_max;

  calc_min_max(ja, jd, jg_min, jg_max, START_MIN_NONE(jb), _max_two_j_cfgh);

  if (jg_max < jg_min) // no cutting?
    {
      st.null_loop++;
      return;
    }

  int jh_min, jh_max;

  calc_min_max(jb, je, jh_min, jh_max, START_MIN_NONE(ja), _max_two_j_cfgh);

  if (jh_max < jh_min) // no cutting?
    {
      st.null_loop++;
      return;
    }

  wigner9j_symbol js;

  js.j.two_j1 = ja;
  js.j.two_j2 = jb;
  js.j.two_j4 = jd;
  js.j.two_j5 = je;

  for (int jc = jc_min; jc <= jc_max; jc += 2)
    for (int jf = jf_min; jf <= jf_max; jf += 2)
      for (int jg = jg_min; jg <= jg_max; jg += 2)
	for (int jh = jh_min; jh <= jh_max; jh += 2)
	  {
	    int j_start = START_MIN_NONE(ja);
	    
	    int ji_min1 = abs(jc - jf);
	    int ji_min2 = abs(jg - jh);
	    int ji_max1 = jc + jf;
	    int ji_max2 = jg + jh;
	    
	    int ji_min = ji_min1 > ji_min2 ? ji_min1 : ji_min2;
	    int ji_max = ji_max1 < ji_max2 ? ji_max1 : ji_max2;
	    
	    if (ji_min < j_start)
	      ji_min = j_start;

	    if ((ji_min ^ ji_max) & 1) // j_max still has correct odd/even
	      ji_min++;
	    if (ji_max > _max_two_j_i)
	      ji_max = _max_two_j_i;

		/*
	    printf ("LIM: (%2d) %2d %2d [%2d] : [%2d] %2d %2d (%2d)\n",
		    ja,abs(jc-jf),abs(jg-jh),
		    ji_min,ji_max,
		    jc+jf,jg+jh,
		    MAX_J2);
		*/

	    js.j.two_j3 = jc;
	    js.j.two_j6 = jf;
	    js.j.two_j7 = jg;
	    js.j.two_j8 = jh;

	    for (int ji = ji_min; ji <= ji_max; ji += 2)
	    {
	      st.count++;
	      
	      uint64_t simple_key =
		(((uint64_t) ja) << 57) |
		(((uint64_t) jb) << 50) |
		(((uint64_t) jc) << 43) |
		(((uint64_t) jd) << 36) |
		(((uint64_t) je) << 29) |
		(((uint64_t) jf) << 22) |
		(((uint64_t) jg) << 15) |
		(((uint64_t) jh) <<  8) |
		(((uint64_t) ji) <<  1);
	      
	      js.j.two_j9 = ji;

	      if (_raw)
		printf ("%3d %3d %3d  %3d %3d %3d  %3d %3d %3d\n",
			ja, jb, jc, jd, je, jf, jg, jh, ji);
	      else
		{
		  uint64_t key;

		  if (_c14n)
		    key = wigner9j_canonicalise(js.v);
		  else
		    key = simple_key;

		  key &= ~1; // ignore the sign bit
	      
		  insert_buf(key);
		}
	    }
	  }
}

int main(int argc, char *argv[])
{
  _max_two_j_abde = 0;
  _max_two_j_cfgh = 0;
  _max_two_j_i    = 0;

  for (int i = 1; i < argc; i++)
    {
      if (strncmp(argv[i],"--lim=",6) == 0)
	{
	  _max_two_j_abde = atoi(argv[i] + 6);
	  _max_two_j_cfgh = _max_two_j_abde * 2;
	  _max_two_j_i    = _max_two_j_cfgh * 2;
	}
      else if (strncmp(argv[i],"--flat-lim=",11) == 0)
	{
	  _max_two_j_abde = atoi(argv[i] + 11);
	  _max_two_j_cfgh = _max_two_j_abde;
	  _max_two_j_i    = _max_two_j_abde;
	}
      else if (strcmp(argv[i],"--raw") == 0)
	{
	  _raw = 1;
	}
      else if (strcmp(argv[i],"--no-c14n") == 0)
	{
	  _c14n = 0;
	}
      else if (strcmp(argv[i],"--all") == 0)
	{
	  _all = 1;
	}
      else
	{
	  fprintf (stderr, "Bad option: %s.\n",argv[i]);
	  exit(1);
	}
    }

  if (!_raw)
    {
      fastwigxj_header header;

      prepare_header(&header);
      header._type = 9;
      header._c14n = _c14n;
      header._max_two_j = _max_two_j_i;
      write_header(&header);
    }

  /* int ja, jb, jc, jd, je, jf, jg, jh, ji; */

  stats st;

  memset(&st, 0, sizeof (st));

  fprintf (stderr,"Limits: %d %d %d.\n",
	   _max_two_j_abde, _max_two_j_cfgh, _max_two_j_i);
  
  for (int ja = /*START_MIN*/(0); ja <= _max_two_j_abde; ja++)
    for (int jb = START_MIN(ja); jb <= _max_two_j_abde; jb++)
      {
	int jc_min, jc_max;
	
	calc_min_max(ja, jb, jc_min, jc_max, START_MIN_NONE(jb), _max_two_j_cfgh);
	
	if (jc_max < jc_min)
	  {
	    st.null_loop++;
	    continue;
	  }	

	for (int jd = START_MIN(jb); jd <= _max_two_j_abde; jd++)
	  for (int je = START_MIN(ja); je <= _max_two_j_abde; je++)
	    loop_four(st, ja, jb, jd, je, jc_min, jc_max);
      }

  /*uint64_t exhaustive = 0;*/

#if 0
  printf ("--- Exhaustive check ---\n");

  for (int ja = 0; ja <= MAX_J2_abde; ja++)
    for (int jb = 0; jb <= MAX_J2_abde; jb++)
      for (int jc = 0; jc <= MAX_J2_cfghi; jc++)
	for (int jd = 0; jd <= MAX_J2_abde; jd++)
	  for (int je = 0; je <= MAX_J2_abde; je++)
	    for (int jf = 0; jf <= MAX_J2_cfghi; jf++)
	      for (int jg = 0; jg <= MAX_J2_cfghi; jg++)
		for (int jh = 0; jh <= MAX_J2_cfghi; jh++)
		  for (int ji = 0; ji <= MAX_J2_cfghi; ji++)
		    {		      
		      if (trivial_zero_9j(ja, jb, jc,
					  jd, je, jf,
					  jg, jh, ji))
			continue;

		      exhaustive++;

		      uint64_t simple_key =
			(((uint64_t) ja) << 57) |
			(((uint64_t) jb) << 50) |
			(((uint64_t) jc) << 43) |
			(((uint64_t) jd) << 36) |
			(((uint64_t) je) << 29) |
			(((uint64_t) jf) << 22) |
			(((uint64_t) jg) << 15) |
			(((uint64_t) jh) <<  8) |
			(((uint64_t) ji) <<  1);

		      uint64_t key =
			wigner9_canonicalise(ja, jb, jc,
					     jd, je, jf,
					     jg, jh, ji);
		      
		      if (_seen.find(key) == _seen.end())
			{
			  printf ("%2d %2d %2d  %2d %2d %2d  %2d %2d %2d "
				  "%016" PRIx64 "  "
				  "%016" PRIx64 "  "
				  "%2d %2d %2d  %2d %2d %2d  %2d %2d %2d "
				  "\n",
				  ja,jb,jc,jd,je,jf,jg,jh,ji,
				  simple_key,
				  key,
				  (key >> 57) & 0x7f,
				  (key >> 50) & 0x7f,
				  (key >> 43) & 0x7f,
				  (key >> 36) & 0x7f,
				  (key >> 29) & 0x7f,
				  (key >> 22) & 0x7f,
				  (key >> 15) & 0x7f,
				  (key >>  8) & 0x7f,
				  (key >>  1) & 0x7f);
			}
		    }
#endif

  write_buf();

  fprintf (stderr,"Symbols generated: %" PRIu64 ".\n",st.count);

  return 0;
}
