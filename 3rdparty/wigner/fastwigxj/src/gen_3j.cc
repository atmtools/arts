#include <sys/time.h>
#include <float.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#define KEEP_MAX_RASCH_YU_INDEX 1

#include "fastwigxj.h"
#include "wigner36j_regge_canonicalize.h"

#include "gen_buf.hh"
#include "gen_header.hh"

int _max_two_j_ab;
int _max_two_j_c;

bool _c14n_regge = false;

struct stats
{
  uint64_t count;
  uint64_t null_loop;
};

// ja jb jc
// ma mb mc

// Generation order:
// ja, jb          (exhaustive)
// ma, mb          (according to limits)
// mc              (given)
// jc              (according to triangular conditions)

//#define START_MIN(jmin) (jmin)
#define START_MIN(jmin) 0

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

void loop_four(stats &st, 
	       int ja, int jb, int ma, int mb,
	       int jc_min, int jc_max)
{
  for (int jc = jc_min; jc <= jc_max; jc++)
    {
      st.count++;

      int jv[5] = { ja, jb, jc, ma, mb };
	    
      uint64_t key =
	wigner3j_regge_canonicalise(jv);
      
      key &= ~1; // ignore the sign bit
      
      insert_buf(key);
    }
}

int main(int argc, char *argv[])
{
  _max_two_j_ab  = 0;
  _max_two_j_c   = 0;

  for (int i = 1; i < argc; i++)
    {
      if (strncmp(argv[i],"--lim=",6) == 0)
	{
	  _max_two_j_ab  = atoi(argv[i] + 6);
	  _max_two_j_c   = _max_two_j_ab * 2;
	}
      else if (strncmp(argv[i],"--flat-lim=",11) == 0)
	{
	  _max_two_j_ab  = atoi(argv[i] + 11);
	  _max_two_j_c   = _max_two_j_ab;
	}
      else
	{
	  fprintf (stderr, "Bad option: %s.\n",argv[i]);
	  exit(1);
	}
    }

  fastwigxj_header header;

  prepare_header(&header);
  header._type = 3;
  header._c14n = 1;
  header._max_two_j = _max_two_j_c;
  write_header(&header);

  /*int ja, jb, jc, jd, je, jf;*/

  stats st;

  memset(&st, 0, sizeof (st));

  fprintf (stderr,"Limits: %d %d.\n",
	   _max_two_j_ab, _max_two_j_c);
  
  for (int ja = START_MIN(0); ja <= _max_two_j_ab; ja++)
    for (int jb = START_MIN(ja); jb <= _max_two_j_ab; jb++)
      for (int ma = -ja; ma <= ja; ma += 2)
	for (int mb = -jb; mb <= jb; mb += 2)
	  {
	    int mc = - ma - mb;

	    int jc_min, jc_max;
	    
	    calc_min_max(ja, jb, jc_min, jc_max, abs(mc), _max_two_j_c);
	    
	    if (jc_max < jc_min)
	      {
		st.null_loop++;
		continue;
	      }
	    
	    loop_four(st, ja, jb, ma, mb, jc_min, jc_max);
	  }

  write_buf();

  fprintf (stderr,"Symbols generated: %" PRIu64 ".\n",st.count);

  if (_c14n_regge)
    fprintf (stderr, "Max Rasch-Yu index: %" PRIu64 ".\n", _max_RaschYu_index);

  return 0;
}
