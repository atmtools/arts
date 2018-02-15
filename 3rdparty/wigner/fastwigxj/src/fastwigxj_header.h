#ifndef __FASTWIGXJ_HEADER_H__
#define __FASTWIGXJ_HEADER_H__

#include <stdint.h>
#include <stdlib.h>

#define FASTWIGXJ_MAGIC   0x9ac65238
#define FASTWIGXJ_VERSION          6 /* change when changing format */

struct fastwigxj_header /* actually a footer */
{
  /* File/header version. */
  uint32_t _magic;
  uint32_t _version;

  /* Content description. */
  uint64_t _entries;
  int      _type;      // 3, 6 or 9
  int      _c14n;      /* 0 = none (=indexed),
                        * 1 = hash table,
			* 2 = quadmath (=indexed) */
  int      _max_two_j;

  /* Information about the hash table. */
  uint64_t _hashentries;
  /* These entries are for verification.  The actual values are
   * for performance compiled into the code.
   */
  char     _hash_func_descr[32];
  uint64_t _hash_param[10];

  /* Some statistics from the pre-calculation. */
  double   _max_abs_err;
  double   _max_rel_err_gt_1em10; // for |values| >  1e-10
  double   _max_small_val;        // for |values| <= 1e-10

  /* Cheap test against corrupted files. */
  uint64_t _checksum_xor;
  uint64_t _checksum_xor_header;

  /* Another magic at the end, since the header is a footer... */
  uint32_t _magic_end;
};

#endif/*__FASTWIGXJ_HEADER_H__*/
