
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
#include <math.h>
#include <float.h>
#include <sys/time.h>

#include <algorithm>

#include "fastwigxj_config.h"

#include "wigxjpf.h"
#if FASTWIGXJ_USE_FLOAT128
#include "wigxjpf_quadmath.h"
#include <quadmath.h>
#endif

#include "gen_header.hh"

#include "fastwigxj_header.h"
#include "fastwigxj_hash_fcn.h"

uint64_t _max_index = 0;

double _max_abs_err = 0;
double _max_rel_err_gt_1em10 = 0;
double _max_small_val = 0;

struct wigxjpf_temp *_wigxjpf_temp = NULL;

void update_error(double val, double err)
{
  if (err > _max_abs_err)
    {
      _max_abs_err = err;
      // printf ("%e %e\n",result.val,result.err);
    }

  if (fabs(val) > 1e-10)
    {
      double rel_err = err / fabs(val);
      
      if (rel_err > _max_rel_err_gt_1em10)
	_max_rel_err_gt_1em10 = rel_err;
    }
  else
    {
      double small_val = fabs(val);

      if (small_val > _max_small_val)
	_max_small_val = small_val;
    }
}

bool _calc_dummy = false;
bool _calc_sign_update_error = true;

void dummy_result_9(int a, int b, int c,
		    int d, int e, int f,
		    int g, int h, int i,
		    double *val, double *err)
{
  *val = a + b + c + d + e + f + g + h + i;
  *err = 0;
}

void dummy_result_6(int a, int b, int c,
		    int d, int e, int f,
		    double *val, double *err)
{
  *val = a + b + c + d + e + f;
  *err = 0;
}

#if FASTWIGXJ_USE_LONG_DOUBLE
#define WIGXJPF_DOUBLE        long double
// The rounding error is not given by wigxjpf itself, but
// our conversion to double.
#define WIGXJPF_DOUBLE_MAX_REL_ERR   (DBL_EPSILON*0.5)
#define WIGXJPF_CALC_FNC(x)   x##_long_double
#else
#define WIGXJPF_DOUBLE        double
// DBL_EPSILON from math.h is not the rounding error, but the minimum
// difference from 1.0 that is representable
#define WIGXJPF_DOUBLE_MAX_REL_ERR   ((DBL_EPSILON*0.5)*6)
#define WIGXJPF_CALC_FNC(x)   x##_double
#endif
#if FASTWIGXJ_USE_FLOAT128
#define WIGXJPF_FLOAT128      __float128
#define WIGXJPF_FLOAT128_MAX_REL_ERR   (FLT128_EPSILON*0.5)
#define WIGXJPF_CALC_FNC_FLOAT128(x)   x##_float128
#endif

double calc_9j_val_for_key(uint64_t key, uint64_t &sign)
{
  int ja = (int) (key >> 57) & 0x7f;
  int jb = (int) (key >> 50) & 0x7f;
  int jc = (int) (key >> 43) & 0x7f;
  int jd = (int) (key >> 36) & 0x7f;
  int je = (int) (key >> 29) & 0x7f;
  int jf = (int) (key >> 22) & 0x7f;
  int jg = (int) (key >> 15) & 0x7f;
  int jh = (int) (key >>  8) & 0x7f;
  int ji = (int) (key >>  1) & 0x7f;

  double val, err;
  
  if (_calc_dummy)
    {
      dummy_result_9(ja, jb, jc,
		     jd, je, jf,
		     jg, jh, ji,
		     &val, &err);
    }
  else
    {
      WIGXJPF_DOUBLE v;
      WIGXJPF_CALC_FNC(calc_9j)(&v,
				ja, jb, jc,
				jd, je, jf,
				jg, jh, ji,
				_wigxjpf_temp);
      val = (double) v;
      err = (double) (WIGXJPF_DOUBLE_MAX_REL_ERR * v);
    }

  if (_calc_sign_update_error)
    {
      sign = ((ja + jb + jc + jd + je + jf + jg + jh + ji) >> 1) & 1;

      update_error(val, err);
    }

  return val;
}

double calc_3j_val_for_key(uint64_t key, uint64_t &sign)
{
  int ja, jb, jc, ma, mb, mc;

  int S, B, T, X, L;

  S = (int) (key >> 52) & 0xfff;
  L = (int) (key >> 40) & 0xfff;
  X = (int) (key >> 28) & 0xfff;
  B = (int) (key >> 16) & 0xfff;
  T = (int) (key >>  4) & 0xfff;

  ja = (B + L - T + X) / 2;
  jb = (B + S - T + X) / 2;
  jc = (L + S) / 2;
  ma = ja - X;
  mb = jb - B;
  mc = - ma - mb;
  /*
  printf ("%016" PRIx64 "%3d %3d %3d %3d %3d : %3d %3d %3d  %3d %3d %3d  ",
	  key, S, L, X, B, T, ja, jb, jc, ma, mb, mc);
  */
  double val, err;

  if (_calc_dummy)
    {
      dummy_result_6(ja, jb, jc,
		     ma, mb, mc,
		     &val, &err);
    }
  else
    {
      WIGXJPF_DOUBLE v;
      WIGXJPF_CALC_FNC(calc_3j)(&v,
				ja, jb, jc,
				ma, mb, mc,
				_wigxjpf_temp);
      val = (double) v;
      err = (double) (WIGXJPF_DOUBLE_MAX_REL_ERR * v);
    }
  /*
  printf ("%15.10f\n",result.val);
  */
  if (_calc_sign_update_error)
    {
      sign = ((ja + jb + jc) >> 1) & 1;

      update_error(val, err);
    }

  return val;
}

typedef double (*calc_js_val_for_key)(uint64_t key, uint64_t &sign);

calc_js_val_for_key _calc_js_val_for_key = NULL;

struct hash_item
{
  uint64_t _key;
  double   _value;
};

uint64_t _collisions     = 0;
uint64_t _collisions_sqr = 0;

void insert_entry(hash_item *hash_table, uint64_t hash_table_mask,
		  uint64_t key, uint64_t sign, double value)
{
  uint64_t x = key;
  
  x = WIGNER369_HASH_FCN(x);
  
  x = x & hash_table_mask;
  
  uint64_t coll = 0;
  
  while (hash_table[x]._key != (uint64_t) -1)
    {
      x++;
      x &= hash_table_mask;
      coll++;
    }
  
  _collisions += coll;
  _collisions_sqr += coll * coll;
  
  hash_table[x]._key   = key | sign;
  hash_table[x]._value = value;
}

void fill_entries(hash_item *hash_table, uint64_t hash_table_mask,
		  uint64_t *buf, uint64_t entries)
{
  for (size_t i = 0; i < entries; i++)
    {
      uint64_t key = buf[i];

      uint64_t sign = 0;

      double value =
	_calc_js_val_for_key(key, sign);

      insert_entry(hash_table, hash_table_mask, key, sign, value);
    }
}

//#define FBUF_SIZE (16*1024*1024)
#define FBUF_SIZE (2*1024*1024)

uint64_t in[FBUF_SIZE];

void fill_entries(FILE *fin, long startpos,
		  uint64_t entries,
		  hash_item *hash, uint64_t hash_table_mask,
		  size_t *total_calc, double *elapsed)
{
  fseek(fin, startpos, SEEK_SET);

  /* uint64_t in[FBUF_SIZE]; */

  struct timeval t1;
  struct timeval t2;

  gettimeofday(&t1, NULL);

  uint64_t total = 0;
  time_t lastprint = t1.tv_sec;

  for ( ; ; )
    {
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

      (void) total_calc;

	fill_entries(hash, hash_table_mask, in, n);

      total += n;

      gettimeofday(&t2, NULL);

      if (t2.tv_sec != lastprint)
	{
	  lastprint = t2.tv_sec;
	  printf ("%" PRIu64 " (%.1f%%)\r",
		  total, 100. * (double) total / (double) entries);
	  fflush(stdout);
	}
    }

  if (total != entries)
    {
      fprintf (stderr, "Did not get the expected number of entries "
	       "(%" PRIu64 " != %" PRIu64 ").\n",
	       total, entries);
      exit(1);
    }

  *elapsed =
    /*   */ (double) (t2.tv_sec - t1.tv_sec) +
    1.e-6 * (double) (t2.tv_usec - t1.tv_usec);
}

void fill_entries(FILE *fin, long startpos,
                  uint64_t entries,
                  hash_item *hash, uint64_t hash_table_mask)
{
      double elapsed;

      fill_entries(fin, startpos,
		   entries,
		   hash, hash_table_mask,
		   NULL, &elapsed);
}

void clear_entries(hash_item *hash, uint64_t hashentries)
{
  for (size_t i = 0; i < hashentries; i++)
    {
      hash[i]._key = (uint64_t) -1;
      hash[i]._value = 0;
    }
}

void fill_table(FILE *fin, long startpos, 
		hash_item *hash, uint64_t hashentries,
		uint64_t entries)
{
  uint64_t hash_table_mask = hashentries - 1;

  clear_entries(hash, hashentries);

  _collisions     = 0;
  _collisions_sqr = 0;

  /* Insert the value 0.0 for trival-0 elements. */
  insert_entry(hash, hash_table_mask, (uint64_t) -2, 0, 0.0);
  /* The special key -4 shall not exist, it marks symbols that
   * did an overflow in/before c14n!
   */

  fill_entries(fin, startpos, entries, hash, hash_table_mask);

  double fill_factor = (double) entries / (double) hashentries;

  fprintf (stderr,"Collisions: %" PRIu64 ".\n", _collisions);

  double limit_coll_tot = -(log(1 - fill_factor) + fill_factor);
  double limit_coll     = limit_coll_tot / fill_factor;

  fprintf (stderr,"Coll factor: %.3f. (%.3f %.3f)\n",
	   (double) _collisions / (double) entries,
	   limit_coll, limit_coll_tot);
  fprintf (stderr,"Coll sqrfct: %.3f.\n",
	   (double) _collisions_sqr / (double) entries);

}

void try_calc_3j(int S, int B, int T, int X, int L,
		 double *val)
{
  int two_ja, two_jb, two_jc, two_ma, two_mb, two_mc;
  double err;

  two_ja = (B + L - T + X);
  two_jb = (B + S - T + X);
  two_jc = (L + S);
  two_ma = two_ja - 2*X;
  two_mb = two_jb - 2*B;
  two_mc = - two_ma - two_mb;

  /*
  printf ("%2d %2d %2d %2d %2d : %3d %3d %3d %3d %3d %3d : ",
	  L, X, T, B, S,
	  two_ja, two_jb, two_jc, two_ma, two_mb, two_mc);
  */
  
  WIGXJPF_DOUBLE v;
  WIGXJPF_CALC_FNC(calc_3j)(&v,
			    two_ja, two_jb, two_jc,
			    two_ma, two_mb, two_mc,
			    _wigxjpf_temp);
  *val = (double) v;

  err = WIGXJPF_DOUBLE_MAX_REL_ERR * (double) v;
  /*
  printf ("%.15f\n",*val);
  */

  update_error(*val, err);
}

void try_calc_6j(int S, int B, int T, int X, int L, int E,
		 double *val)
{
  int two_ja, two_jb, two_jc, two_jd, two_je, two_jf;
  double err;

  two_ja = (L + T);
  two_jb = (B + E);
  two_jc = (B + E + L + T - 2*X);
  two_jd = (B + L + S - X);
  two_je = (E + S + T - X);
  two_jf = (E + L + S - X);
  /*
  printf ("%d %d %d %d %d %d ", ja, jb, jc, jd, je, jf);
  */
  WIGXJPF_DOUBLE v;
  WIGXJPF_CALC_FNC(calc_6j)(&v,
			    two_ja, two_jb, two_jc,
			    two_jd, two_je, two_jf,
			    _wigxjpf_temp);
  *val = (double) v;
  err = WIGXJPF_DOUBLE_MAX_REL_ERR * (double) v;
  
  update_error(*val, err);
}

#if FASTWIGXJ_USE_FLOAT128
void try_calc_6j_float128(int S, int B, int T, int X, int L, int E,
			  __float128 *val)
{
  int two_ja, two_jb, two_jc, two_jd, two_je, two_jf;
  double err;

  two_ja = (L + T);
  two_jb = (B + E);
  two_jc = (B + E + L + T - 2*X);
  two_jd = (B + L + S - X);
  two_je = (E + S + T - X);
  two_jf = (E + L + S - X);
  /*
  printf ("%d %d %d %d %d %d ", ja, jb, jc, jd, je, jf);
  */
  WIGXJPF_FLOAT128 v;
  WIGXJPF_CALC_FNC_FLOAT128(calc_6j)(&v,
				     two_ja, two_jb, two_jc,
				     two_jd, two_je, two_jf,
				     _wigxjpf_temp);
  *val = v;
  err = (double) WIGXJPF_FLOAT128_MAX_REL_ERR * (double) v;
  
  update_error((double) *val, err);
}
#endif

void fill_table_3j(double *table, uint64_t entries)
{
  if (entries > 0)
    table[0] = 0.0;

#define DIVIDE_INDEX_3J (120)

  for (uint64_t L = 0; ; L++)
    {
      uint64_t Lind = L * (24 + L * (50 + L * (35 + L * (10 + L))));
      uint64_t indexL = (Lind);

      if (indexL / DIVIDE_INDEX_3J >= _max_index)
	break;

      for (uint64_t X = 0; X <= L; X++)
	{
	  uint64_t Xind = X * (6 + X * (11 + X * (6 + X)));
	  uint64_t indexX = indexL + (5  * Xind);

	  if (indexX / DIVIDE_INDEX_3J >= _max_index)
	    break;

	  for (uint64_t T = 0; T <= X; T++)
	    {
	      uint64_t Tind = T * (2 + T * (3 + T));
	      uint64_t indexT = indexX + (20 * Tind);

	      if (indexT / DIVIDE_INDEX_3J >= _max_index)
		break;
		  
	      for (uint64_t B = 0; B <= T; B++)
		{
		  uint64_t Bind = B * (1 + B);
		  uint64_t indexB = indexT + (60 * Bind);
		      
		  if (indexB / DIVIDE_INDEX_3J >= _max_index)
		    break;

		  for (uint64_t S = 0; S <= B; S++)
		    {
		      uint64_t Sind = S;
		      uint64_t indexS = indexB + (120 * Sind);
			  
		      if (indexS / DIVIDE_INDEX_3J + 1 >= _max_index)
			break;

		      uint64_t index = (indexS) / DIVIDE_INDEX_3J + 1;

		      try_calc_3j((int) S, (int) B, (int) T, (int) X, (int) L,
				  &table[index]);
		    }
		}
	    }
	}
    }
}

void fill_table_6j(double *table, uint64_t entries, int item_float128)
{
  if (entries > 0)
    table[0] = 0.0;

#define DIVIDE_INDEX_6J (720)

  for (uint64_t E = 0; ; E++)
    {
      uint64_t Eind = E * (120 + E * (274 + E * (225 + E * (85 + E * (15 + E)))));
      uint64_t indexE = (Eind);

      if (indexE / DIVIDE_INDEX_6J >= _max_index)
	break;

      for (uint64_t L = 0; L <= E; L++)
	{
	  uint64_t Lind = L * (24 + L * (50 + L * (35 + L * (10 + L))));
	  uint64_t indexL = indexE + (6   * Lind);

	  if (indexL / DIVIDE_INDEX_6J >= _max_index)
	    break;

	  for (uint64_t X = 0; X <= L; X++)
	    {
	      uint64_t Xind = X * (6 + X * (11 + X * (6 + X)));
	      uint64_t indexX = indexL + (30  * Xind);

	      if (indexX / DIVIDE_INDEX_6J >= _max_index)
		break;

	      for (uint64_t T = 0; T <= X; T++)
		{
		  uint64_t Tind = T * (2 + T * (3 + T));
		  uint64_t indexT = indexX + (120 * Tind);

		  if (indexT / DIVIDE_INDEX_6J >= _max_index)
		    break;
		  
		  for (uint64_t B = 0; B <= T; B++)
		    {
		      uint64_t Bind = B * (1 + B);
		      uint64_t indexB = indexT + (360 * Bind);
		      
		      if (indexB / DIVIDE_INDEX_6J >= _max_index)
			break;

		      for (uint64_t S = 0; S <= B; S++)
			{
			  uint64_t Sind = S;
			  uint64_t indexS = indexB + (720 * Sind);
			  
			  if (indexS / DIVIDE_INDEX_6J + 1 >= _max_index)
			    break;

			  uint64_t index = (indexS) / DIVIDE_INDEX_6J + 1;

			  if (item_float128)
#if FASTWIGXJ_USE_FLOAT128
			    try_calc_6j_float128((int) S, (int) B, (int) T,
						 (int) X, (int) L, (int) E,
						 &((__float128 *)table)[index]);
#else
			    /* This should have been caught earlier! */
			    { fprintf(stderr,"Bug!\n"); exit(1); }
#endif
			  else
			    try_calc_6j((int) S, (int) B, (int) T,
					(int) X, (int) L, (int) E,
					&table[index]);
			}
		    }
		}
	    }
	}
    }
}

void usage(const char *cmd)
{
  fprintf (stderr,"Usage:  %s [opt] infile outfile\n",cmd);
}

uint64_t array_xor(const void *array, size_t size)
{
  size_t n = size / sizeof (uint64_t);
  const uint64_t *p = (const uint64_t *) array;
  uint64_t s = 0;

  for (size_t i = 0; i < n; i++)
    s ^= p[i];

  return s;
}

int main(int argc, char *argv[])
{
  const char *infile  = NULL;
  const char *outfile = NULL;

  int max_two_j = -1;
  int table_type = -1;
  int item_float128 = 0;

  for (int i = 1; i < argc; i++)
    {
      if (strncmp(argv[i],"--max-index-3j=",15) == 0)
	{
	  _max_index = atoll(argv[i]+15);

	  uint64_t L = 0;

	  for ( ; ; )
	    {
	      uint64_t max_index = 
		(L * (274 + L * (225 + L * (85 + L * (15 + L))))
		 / 120) + 1;

	      L++;

	      if (max_index >= _max_index)
		break;
	    }

	  max_two_j = (int) (2 * L);
	  table_type = 3;
	}
      else if (strncmp(argv[i],"--max-E-3j=",11) == 0)
	{
	  uint64_t L = atoll(argv[i]+11);

	  _max_index = 
	    (L * (274 + L * (225 + L * (85 + L * (15 + L))))
	     / 120) + 1;

	  /* And add one more for the first entry for trivial 0 */

	  _max_index++;

	  max_two_j = (int) (2 * L);
	  table_type = 3;

	  printf ("3j: %2" PRIu64 " -> %2" PRIu64 " \n", L, _max_index);
	}
      else if (strncmp(argv[i],"--max-index-6j=",15) == 0)
	{
	  _max_index = atoll(argv[i]+15);

	  uint64_t E = 0;

	  for ( ; ; )
	    {
	      uint64_t max_index = 
		(E * (1764 + E * (1624 + E * (735 + E * (175 + E * (21 + E)))))
		 / 720) + 1;

	      E++;

	      if (max_index >= _max_index)
		break;
	    }

	  max_two_j = (int) (2 * E);
	  table_type = 6;
	}
      else if (strncmp(argv[i],"--max-E-6j=",11) == 0)
	{
	  uint64_t E = atoll(argv[i]+11);

	  _max_index = 
	    (E * (1764 + E * (1624 + E * (735 + E * (175 + E * (21 + E)))))
	     / 720) + 1;

	  /* And add one more for the first entry for trivial 0 */

	  _max_index++;

	  max_two_j = (int) (2 * E);
	  table_type = 6;

	  printf ("6j: %2" PRIu64 " -> %2" PRIu64 " \n", E, _max_index);
	}
      else if (strcmp(argv[i],"--float128") == 0)
	{
#if FASTWIGXJ_USE_FLOAT128
	  item_float128 = 1;
#else
	  fprintf (stderr, "Support for float128 not available/built.\n");
	  exit(1);
#endif
	}
      else if (strcmp(argv[i],"--help") == 0)
	{
	  usage(argv[0]);
	  exit(0);
	}
      else
	{
	  if (!infile)
	    infile = argv[i];
	  else if (!outfile)
	    outfile = argv[i];
	  else
	    {
	      printf ("Unhandled option '%s'\n", argv[i]);
	      usage(argv[0]);
	      exit(1);
	    }
	}
    }

  if (!infile || !outfile)
    {
      usage(argv[0]);
      exit(1);
    }
  
  FILE *fin = NULL;
  fastwigxj_header header;
  uint64_t hashentries = 0;
  uint64_t entries = 0;
  long startpos = 0;
  long endpos = 0;

  if (!_max_index)
    {
      fin = fopen(infile,"r");

      if (fin == NULL)
	{
	  fprintf (stderr, "Failure opening '%s' for reading.\n", infile);
	  exit(1);
	}

      read_header(&header, fin);

      if (header._type == 3)
	{
	  _calc_js_val_for_key = calc_3j_val_for_key;
	}
      else if (header._type == 9)
	{
	  _calc_js_val_for_key = calc_9j_val_for_key;
	}
      else
	{
	  fprintf (stderr, "Unknown type of js: %d.\n", header._type);
	  exit(1);
	}

      startpos = ftell(fin);
      fseek(fin, 0, SEEK_END);
      endpos = ftell(fin);
      fseek(fin, startpos, SEEK_SET);

      entries = (endpos - startpos) / sizeof (uint64_t);

      fprintf (stderr, "Entries: %" PRIu64 ".\n", entries);

      hashentries = 1;
      while (hashentries * 2 < entries * 3)
	hashentries *= 2;
      /* We must have at least one empty item, and one for the trivial 0. */
      if (hashentries < entries + 2)
	hashentries *= 2;

      fprintf (stderr, "Hash entries: %" PRIu64 ".\n", hashentries);

      double fill_factor = (double) entries / (double) hashentries;

      fprintf (stderr, "Fill factor: %.3f.\n", fill_factor);
    }
  else
    {
      hashentries = _max_index;

      fprintf (stderr, "Table entries: %" PRIu64 ".\n", hashentries);

      prepare_header(&header);
      header._type = table_type;
      header._c14n = (item_float128 ? 2 : 0);
      header._max_two_j = max_two_j;
    }

  size_t entry_size;

  switch (header._c14n)
    {
    case 0: entry_size = sizeof (double); break;
    case 1: entry_size = sizeof (hash_item); break;
    case 2:
#if FASTWIGXJ_USE_FLOAT128
      entry_size = sizeof (__float128);
#else
      fprintf (stderr, "Support for float128 not available/built.\n");
      exit(1);
#endif
      if (table_type != 6)
	{
	  /* This is just because we only use it for 9j fallback.
	   * WIGXJPF could do the calculations.
	   */
	  fprintf (stderr, "Support for float128 only for 6j.\n");
	  exit(1);  
	}
      break;
    default:
        fprintf (stderr, "Bad item type (c14n = %d).\n",
		 header._c14n);
      exit(1);    
    }

  size_t hashsz = hashentries * entry_size;

  fprintf (stderr, "Hash/table size: %.2f GB, %zd B/entry.\n",
	   (double) hashsz * 1.e-9, entry_size);

  hash_item *hash = (hash_item *) malloc (hashsz);

  if (max_two_j == -1)
    max_two_j = header._max_two_j;

  int max_factorial =
    (/*3->3, 6->4, 9->5*/(header._type / 3 + 2) * max_two_j) / 2 + 1;
  int max_iter = max_two_j / 2 + 1;

  wigxjpf_fill_factors(max_factorial);
  _wigxjpf_temp = wigxjpf_temp_alloc(max_iter);

  if (!_max_index)
    {
	fill_table(fin, startpos, hash, hashentries, entries);

      fclose(fin);
    }
  else if (table_type == 3)
    {
      fill_table_3j((double *) hash, entries);
    }
  else if (table_type == 6)
    {
      fill_table_6j((double *) hash, entries, item_float128);
    }

  wigxjpf_temp_free(_wigxjpf_temp);
  _wigxjpf_temp = NULL;
  wigxjpf_fill_factors(0);
  

  // header._magic       = WIGNER69J_MAGIC;   // set by gen_xj
  // header._version     = WIGNER69J_VERSION; // set by gen_xj
  header._entries     = entries + 1; /* +1 for the trivial-0 entry. */
  // header._type        = type;              // set by gen_xj

  header._hashentries = hashentries;
  if (!_max_index)
    {
      strcpy(header._hash_func_descr, WIGNER369_HASH_FCN_DESCR);
      uint64_t hash_params[10] = WIGNER369_HASH_FCN_PARAMS;
      if (sizeof (hash_params) != sizeof (header._hash_param))
	{
	  fprintf (stderr, "Internal error - "
		   "wrong number of parameters for hash function.\n");
	  exit(1);
	}
      memcpy(header._hash_param, hash_params, sizeof (hash_params));
    }
  header._max_abs_err = _max_abs_err;
  header._max_rel_err_gt_1em10 = _max_rel_err_gt_1em10;
  header._max_small_val = _max_small_val;
  header._checksum_xor = array_xor(hash, hashsz);
  header._checksum_xor_header = 0;
  header._checksum_xor_header = array_xor(&header, sizeof (header));

  FILE *fout = fopen(outfile,"w");

  if (fout == NULL)
    {
      fprintf (stderr, "Failure opening '%s' for writing.\n", outfile);
      exit(1);
    }

  if (fwrite(hash,
	     (header._c14n ? sizeof(hash_item) : sizeof (double)),
	     hashentries, fout) != hashentries ||
      fwrite(&header, sizeof (header), 1, fout) != 1)
    {
      fprintf (stderr, "Failure writing to hash file.\n");
      exit(1);
    }
  if (fclose(fout) != 0)
    {
      perror("fclose");
      fprintf (stderr, "Failure when closing hash file.\n");
      exit(1);      
    }

  return 0;
}
