/*
 * Copyright (c) 2010-2020 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
 * 
 * nathanael.schaeffer@univ-grenoble-alpes.fr
 * 
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * 
 */

#include "sht_private.h"

#ifndef SHTNS4MAGIC

	#define S2D_STORE(mem, idx, n, s) \
		vstor(((double*)mem), idx, n); \
		vstor(((double*)mem) + NLAT-VSIZE2, -(idx), vreverse(s));

  #ifdef _GCC_VEC_
	void inline static
	cstore_north_south(double* mem, double* mem_m, long idx, long nlat, rnd nr, rnd sr, rnd ni, rnd si) {
		ni = vxchg_even_odd(ni);	sr = vreverse(sr);
		rnd aa = nr + ni;		rnd bb = nr - ni;
		const long ridx = nlat - (idx+1)*VSIZE2;
		si = vreverse_pairs(si);
		rnd cc = sr + si;		rnd dd = sr - si;
		nr = vblend_even_odd( bb, aa );		ni = vblend_even_odd( aa, bb );
		vstor(mem, idx, nr);
		sr = vblend_even_odd( dd, cc );		si = vblend_even_odd( cc, dd );
		vstor(mem_m, idx, ni);
		vstor(mem + ridx, 0, sr);
		vstor(mem_m + ridx, 0, si);
	}
  #else
	#define S2D_CSTOREX(mem, idx, v, nr, sr, ni, si) { \
		double a0 = (nr[(v)])   + (ni[(v)+1]); \
		double b0 = (nr[(v)])   - (ni[(v)+1]); \
		double a1 = (nr[(v)+1]) + (ni[(v)]); \
		double b1 = (nr[(v)+1]) - (ni[(v)]); \
		((cplx*)mem)[(idx)] = b0 + I*a1; \
		((cplx*)mem)[(NPHI-2*im)*(shtns->nlat_padded>>1) + (idx)] = a0 + I*b1; \
		a1 = (sr[(v)])   + (si[(v)+1]); \
		b1 = (sr[(v)])   - (si[(v)+1]); \
		a0 = (sr[(v)+1]) + (si[(v)]); \
		b0 = (sr[(v)+1]) - (si[(v)]); \
		((cplx*)mem)[NLAT_2-1 -(idx)] = b0 + I*a1; \
		((cplx*)mem)[(NPHI-2*im)*(shtns->nlat_padded>>1) +NLAT_2-1 -(idx)] = a0 + I*b1; }
  #endif

/*	// "NATIVE" LAYOUT ?
	void inline static
	cstore_north_south(double* mem, double* mem_m, long idx, rnd nr, rnd sr, rnd ni, rnd si) {
		ni = vxchg_even_odd(ni);	si = vxchg_even_odd(si);
		rnd aa = nr + ni;		rnd bb = nr - ni;
		rnd cc = sr + si;		rnd dd = sr - si;
		nr = vblend_even_odd( bb, aa );		ni = vblend_even_odd( aa, bb );
		sr = vblend_even_odd( dd, cc );		si = vblend_even_odd( cc, dd );
		vstor(mem, idx, nr);
		vstor(mem, idx+1, sr);
		vstor(mem_m, idx, ni);
		vstor(mem_m, idx+1, si);
	}	*/
/*	void inline static
	cstore_north_south(double* mem, double* mem_m, long idx, long nlat, rnd nr, rnd sr, rnd ni, rnd si) {
		rnd neg_even = vneg_even_xor_cte;	ni = vxor(ni, neg_even);		si = vxor(si, neg_even);		// change sign of even values in vector
		ni = vxchg_even_odd(ni);	sr = vreverse(sr);		si = vreverse_pairs(si);
		const long ridx = nlat - (idx+1)*VSIZE2;
		vstor(mem, idx, nr - ni);
		vstor(mem_m, idx, nr + ni);
		vstor(mem + ridx, 0, sr + si);
		vstor(mem_m + ridx, 0, sr - si);
	}	*/

#else
	inline static void S2D_STORE_4MAGIC(double* mem, long idx, rnd n, rnd s) {
		vinterleave(n,s);
		vstor(mem, idx*2,   n);		vstor(mem, idx*2+1, s);
	}

	inline static void S2D_CSTORE_4MAGIC(double* mem, double* mem_m, long idx, rnd nr, rnd sr, rnd ni, rnd si) {
		rnd a0 = nr-si;		rnd a1 = sr+ni;
		rnd b0 = nr+si;		rnd b1 = sr-ni;
		vinterleave(a0,a1);
		vinterleave(b0,b1);
		vstor(mem,   idx*2, a0);	vstor(mem,   idx*2+1, a1);
		vstor(mem_m, idx*2, b0);	vstor(mem_m, idx*2+1, b1);
	}
#endif


#define MTR MMAX
#define SHT_VAR_LTR

#define GEN(name,sfx) GLUE2(name,sfx)
#define GEN3(name,nw,sfx) GLUE3(name,nw,sfx)

// genaral case, hi lmax
#undef SUFFIX
#define SUFFIX _l
#define HI_LLIM

  #ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
  #endif
	#define NWAY 2
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
  #if VSIZE2 <= 4
	#define NWAY 6
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
  #endif

#define SHT_GRAD
	#define NWAY 2
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#endif
	#define NWAY 4
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
#undef SHT_GRAD

#define SHT_3COMP
  #ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
  #endif
	#define NWAY 2
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
  #if VSIZE2 <= 4
	#define NWAY 6
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
  #endif
#undef SHT_3COMP

// genaral case, low lmax
#undef SUFFIX
#define SUFFIX _l
#undef HI_LLIM

  #ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
  #endif
	#define NWAY 2
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
  #if VSIZE2 <= 4
	#define NWAY 6
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
  #endif

#define SHT_GRAD
	#define NWAY 2
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#endif
	#define NWAY 4
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
#undef SHT_GRAD

#define SHT_3COMP
  #ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
  #endif
	#define NWAY 2
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
  #if VSIZE2 <= 4
	#define NWAY 6
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
  #endif
#undef SHT_3COMP


// axisymmetric
#define SHT_AXISYM
#undef SUFFIX
#define SUFFIX _m0l

  #ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
  #endif
	#define NWAY 2
	#include "SHT/SH_to_spat_kernel.c"
	#include "SHT/SHst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/SH_to_spat_kernel.c"
	#undef NWAY
  #if VSIZE2 <= 4
	#define NWAY 6
	#include "SHT/SH_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/SH_to_spat_kernel.c"
	#undef NWAY
  #endif

#define SHT_GRAD
	#define NWAY 2
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
	#endif
	#define NWAY 4
	#include "SHT/SHs_to_spat_kernel.c"
	#include "SHT/SHt_to_spat_kernel.c"
	#undef NWAY
#undef SHT_GRAD

#define SHT_3COMP
  #ifdef _GCC_VEC_
	#define NWAY 1
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
  #endif
	#define NWAY 2
	#include "SHT/SHqst_to_spat_kernel.c"
	#undef NWAY
#undef SHT_3COMP
