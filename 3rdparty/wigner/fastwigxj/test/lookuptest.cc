
#include "fastwigxj_cc.hh"
#include "fastwigxj_header.h"

#include "wigxjpf.h"

#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <float.h>

void usage(const char *cmd)
{
  printf ("Usage:  %s "
	  "[--dyn-3j=N] "
	  "[--dyn-6j=N] "
	  "[--dyn-9j=N] "
	  "--lim=N hashfile [hashfile]\n",cmd);
}

#define WITHIN_TRIANGLE(ja, jb, jc) \
  (((jc) >= abs ((ja) - (jb))) && ((jc) <= (ja) + (jb)))

#define MAX_SPEED_TEST 10000000

wigner9j_symbol totest[MAX_SPEED_TEST];
uint64_t ntotest = 0;

volatile int _dummy_zero = 0;

int lookup_9j(int lim, const char *hashfile,
	      const char *tablefile_6j)
{
  wigner6j_table table_6j;
  wigner9j_hash hash_9j;
  fastwigxj_header header;

  hash_9j.init(hashfile, &header);

  if (tablefile_6j)
    {
      table_6j.init(tablefile_6j, &header);
    }

  /* Try some set of js... */

  int had_error = 0;

  for (int ja = 0; ja <= lim; ja++)
    for (int jb = 0; jb <= lim; jb++)
      for (int jc = 0; jc <= lim; jc++)
	for (int jd = 0; jd <= lim; jd++)
	  for (int je = 0; je <= lim; je++)
	    for (int jf = 0; jf <= lim; jf++)
	      for (int jg = 0; jg <= lim; jg++)
		for (int jh = 0; jh <= lim; jh++)
		  for (int ji = 0; ji <= lim; ji++)
		    {
		      double hash_result;
		      bool got_hash;

		      wigner9j_symbol js;

		      js.j.two_j1 = ja;
		      js.j.two_j2 = jb;
		      js.j.two_j3 = jc;
		      js.j.two_j4 = jd;
		      js.j.two_j5 = je;
		      js.j.two_j6 = jf;
		      js.j.two_j7 = jg;
		      js.j.two_j8 = jh;
		      js.j.two_j9 = ji;

		      bool trivial_0 =
			(((ja + jb + jc) |
			  (jd + je + jf) |
			  (jg + jh + ji) |
			  (ja + jd + jg) |
			  (jb + je + jh) |
			  (jc + jf + ji)) & 1) ||
			!WITHIN_TRIANGLE(ja, jb, jc) ||
			!WITHIN_TRIANGLE(jd, je, jf) ||
			!WITHIN_TRIANGLE(jg, jh, ji) ||
			!WITHIN_TRIANGLE(ja, jd, jg) ||
			!WITHIN_TRIANGLE(jb, je, jh) ||
			!WITHIN_TRIANGLE(jc, jf, ji);

		      hash_result = hash_9j.lookup(js.v);
		      got_hash = true;

		      if (got_hash)
			{
			  double result, result_err;
			  
			  result =
			    wig9jj(ja, jb, jc,
				   jd, je, jf,
				   jg, jh, ji);
			  result_err = fabs(result * 7.e-16);

			  double diff = fabs(result - hash_result);

			  if (diff > result_err &&
			      diff > header._max_abs_err)
			    {
			      fflush (stdout);
			      fprintf (stderr,
				       "%2d %2d %2d  "
				       "%2d %2d %2d  "
				       "%2d %2d %2d  "
				       "hashed:%15.20f "
				       "wigxjpf:%15.20f %15.20f (%4.1f)\n",
				       ja,jb,jc, jd,je,jf, jg,jh,ji,
				       hash_result, result, result_err,
				       log(diff / fabs(result)) / log(10));

			      if (++had_error > 10)
				goto report_errors;
			    }

			  if (!trivial_0)
			    {
			      totest[ntotest % MAX_SPEED_TEST] = js;
			      ntotest++;
			    }
			}
		    }

  while (ntotest < MAX_SPEED_TEST)
    {
      uint64_t moretotest = MAX_SPEED_TEST - ntotest;
      if (moretotest > ntotest)
	moretotest = ntotest;

      memcpy(&totest[ntotest], &totest[0], sizeof (totest[0]) * moretotest);

      ntotest += moretotest;
    }

  hash_9j.deinit();

  if (had_error)
    goto report_errors;

  return 0;

 report_errors:
  fprintf (stderr, "Had errors!\n");
  return 1;
}

template<typename lookup_t>
int lookup_6j(int lim, const char *hashfile)
{
  lookup_t hash_6j;
  fastwigxj_header header;

  hash_6j.init(hashfile, &header);

  /* Try some set of js... */

  int had_error = 0;

  for (int ja = 0; ja <= lim; ja++)
    for (int jb = 0; jb <= lim; jb++)
      for (int jc = 0; jc <= lim; jc++)
	for (int jd = 0; jd <= lim; jd++)
	  for (int je = 0; je <= lim; je++)
	    for (int jf = 0; jf <= lim; jf++)
	      {
		bool trivial_0 =
		  (((ja + jb + jc) |
		    (ja + je + jf) |
		    (jd + jb + jf) |
		    (jd + je + jc)) & 1) ||
		  !WITHIN_TRIANGLE(ja, jb, jc) ||
		  !WITHIN_TRIANGLE(ja, je, jf) ||
		  !WITHIN_TRIANGLE(jd, jb, jf) ||
		  !WITHIN_TRIANGLE(jd, je, jc);

		double hash_result;
		bool got_hash;

		hash_result = hash_6j.lookup(ja, jb, jc,
					     jd, je, jf);
		got_hash = true;

		if (got_hash)
		  {
		    double result, result_err;

		    result =
		      wig6jj(ja, jb, jc,
			     jd, je, jf);
		    result_err = fabs(result * 7.e-16);
		    
		    double diff = fabs(result - hash_result);

		    if (diff > result_err &&
			diff > header._max_abs_err)
		      {
			fprintf (stderr,
				 "%2d %2d %2d  "
				 "%2d %2d %2d  "
				 "hashed:%15.20f "
				 "wigxjpf:%15.20f %15.20f (%4.1f)\n",
				 ja,jb,jc, jd,je,jf,
				 hash_result, result, result_err,
				 log(diff / fabs(result)) / log(10));

			if (++had_error > 10)
			  goto report_errors;
		      }
		  }
		else if (trivial_0)
		  { }
		else
		  {
		    double result, result_err;

		    result =
		      wig6jj(ja, jb, jc,
			     jd, je, jf);
		    result_err = result * 7.e-17;
		    (void) result_err;

		  }
	      }

  hash_6j.deinit();

  if (had_error)
    goto report_errors;

  return 0;

 report_errors:
  fprintf (stderr, "Had errors!\n");
  return 1;
}

int lookup_3j(int lim, const char *hashfile)
{
  wigner3j_table hash_3j;
  fastwigxj_header header;

  hash_3j.init(hashfile, &header);

  /* Try some set of js... */

  int had_error = 0;

  for (int ja = 0; ja <= lim; ja++)
    for (int jb = 0; jb <= lim; jb++)
      for (int jc = (ja + jb) & 1; jc <= lim; jc += 2)
	for (int ma = -ja; ma <= ja; ma += 2)
	  for (int mb = -jb; mb <= jb; mb += 2)
	      {
		int mc = - ma - mb;

		if (mc > jc || -mc > jc)
		  continue;

		bool trivial_0 =
		  ((ja + jb + jc) & 1) ||
		  !WITHIN_TRIANGLE(ja, jb, jc);

		double hash_result;
		bool got_hash = false;

		hash_result = hash_3j.lookup(ja, jb, jc,
					     ma, mb /*, mc*/);
		got_hash = true;

		if (got_hash)
		  {
		    double result, result_err;

		    result =
		      wig3jj(ja, jb, jc,
			     ma, mb, mc);
		    result_err = fabs(result * 7.e-16);

		    double diff = fabs(result - hash_result);

		    if (diff > result_err ||
			diff > header._max_abs_err)
		      {
			fprintf (stderr,
				 "%2d %2d %2d  "
				 "%2d %2d %2d  "
				 "hashed:%15.20f "
				 "wigxjpf:%15.20f %15.20f (%15.20f %4.1f)\n",
				 ja,jb,jc,ma,mb,mc,
				 hash_result, result, result_err, diff,
				 log(diff / fabs(result)) / log(10));

			if (++had_error > 10)
			  goto report_errors;
		      }
		  }
		else if (trivial_0)
		  { }
		else
		  {
		    fprintf (stderr,
			     "%2d %2d %2d  "
			     "%2d %2d %2d  "
			     "miss\n",
			     ja,jb,jc,ma,mb,mc);

		  }
	      }

  hash_3j.deinit();

  if (had_error)
    goto report_errors;

  return 0;

 report_errors:
  fprintf (stderr, "Had errors!\n");
  return 1;
}

int main(int argc, char *argv[])
{
  const char *hashfile[2] = { NULL, NULL };
  int lim = 0;
  int ret;
  size_t dyn_3j = 0;
  size_t dyn_6j = 0;
  size_t dyn_9j = 0;

  for (int i = 1; i < argc; i++)
    {
      if (strncmp(argv[i],"--lim=",6) == 0)
	{
	  lim = atoi(argv[i]+6);
	}
      else if (strncmp(argv[i],"--dyn-3j=",9) == 0)
	{
	  dyn_3j = atoi(argv[i]+9);
	}
      else if (strncmp(argv[i],"--dyn-6j=",9) == 0)
	{
	  dyn_6j = atoi(argv[i]+9);
	}
      else if (strncmp(argv[i],"--dyn-9j=",9) == 0)
	{
	  dyn_9j = atoi(argv[i]+9);
	}
      else
	{
	  if (!hashfile[0])
	    hashfile[0] = argv[i];
	  else if (!hashfile[1])
	    hashfile[1] = argv[i];
	  else
	    {
	      fprintf (stderr, "Bad argument: %s\n", argv[i]);

	      usage(argv[0]);
	      exit(1);
	    }
	}
   }

  if (!lim || !hashfile[0])
    {
      usage(argv[0]);
      exit(1);
    }

  const char *hashfile_9j = NULL;
  const char *tablefile_3j = NULL;

  const char *tablefile_6j = NULL;
  const char *tablefile_6j_float128 = NULL;

  for (int i = 0; i < 2; i++)
    if (hashfile[i])
      {
	wigner369j_hash hash;
	fastwigxj_header header;

	hash.init(hashfile[i], -1, &header);

	hash.deinit();

	if (header._type == 9 && header._c14n == 1)
	  {
	    if (hashfile_9j)
	      {
		fprintf (stderr, "Cannot handle two 9j hash files.\n");
		exit(1);
	      }
	    hashfile_9j = hashfile[i];
	  }
	else if (header._type == 6 && header._c14n == 0)
	  {
	    if (tablefile_6j)
	      {
		fprintf (stderr, "Cannot handle two 6j table files.\n");
		exit(1);
	      }
	    tablefile_6j = hashfile[i];
	  }
	else if (header._type == 6 && header._c14n == 2)
	  {
	    if (tablefile_6j_float128)
	      {
		fprintf (stderr,"Cannot handle two 6j float128 table files.\n");
		exit(1);
	      }
	    tablefile_6j_float128 = hashfile[i];
	  }
	else if (header._type == 3 && header._c14n == 0)
	  {
	    if (tablefile_3j)
	      {
		fprintf (stderr, "Cannot handle two 3j hash files.\n");
		exit(1);
	      }
	    tablefile_3j = hashfile[i];
	  }
	else
	  {
	    fprintf (stderr, "Cannot test hash table of type/c14n: %d/%d\n",
		     header._type, header._c14n);
	    exit(1);
	  }
      }

  wig_table_init(2*100, 9);
  wig_temp_init(2*100);

  if (dyn_3j)
    fastwigxj_dyn_init(3, dyn_3j);
  if (dyn_6j)
    fastwigxj_dyn_init(6, dyn_6j);
  if (dyn_9j)
    fastwigxj_dyn_init(9, dyn_9j);

  ret = -1;

  if (hashfile_9j)
    {
      ret = lookup_9j(lim, hashfile_9j,
		      tablefile_6j);
    }
  else if (tablefile_6j)
    {
      ret = lookup_6j<wigner6j_table>(lim, tablefile_6j);
    }
  else if (tablefile_6j_float128)
    {
#if FASTWIGXJ_USE_FLOAT128
      ret = lookup_6j<wigner6j_table_float128>(lim, tablefile_6j_float128);
#else
      fprintf (stderr, "Support for float128 not available/built.\n");
      exit(1);
#endif
    }
  else if (tablefile_3j)
    {
      ret = lookup_3j(lim, tablefile_3j);
    }

  printf ("\n");
  fastwigxj_print_stats();

  fastwigxj_dyn_free(3);
  fastwigxj_dyn_free(6);
  fastwigxj_dyn_free(9);

  wig_temp_free();
  wig_table_free();

  return ret;
}
