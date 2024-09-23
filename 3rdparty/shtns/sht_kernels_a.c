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
// Fm = F + im*m_inc
// Fm2 = F + (NPHI-im)*m_inc
static void split_north_south_real_imag(double* Fm, double* Fm2, double* eori, long k0v, unsigned nlat, int k_inc)
{
    unsigned nk = ((((nlat+1)>>1) +VSIZE2-1)/VSIZE2)*VSIZE2;
    k0v *= VSIZE2;
  #if VSIZE2 >= 2
	if LIKELY(k_inc == 1) {		// optimized and vectorized for k_inc==1
		const rnd neg_even = vneg_even_xor_cte;
		eori += k0v*4;
		for (long k=k0v; k<nk; k+=VSIZE2) {
			rnd a = vread(Fm+k, 0);
			rnd b = vread(Fm2+k, 0);
			rnd xm = a-b;		rnd xp = a+b;
			xm = vxchg_even_odd(xm);
			xm = vxor(xm, neg_even);	// change sign of even values in vector

			rnd c = vread(Fm+(nlat-VSIZE2-k), 0);
			rnd d = vread(Fm2+(nlat-VSIZE2-k), 0);
			rnd yp = vreverse(c+d);
			rnd ym = d-c;
			ym = vxor(ym, neg_even);	// change sign of even values in vector
			ym = vreverse_pairs(ym);

			vstor(eori, 0, xp);		// north real
			vstor(eori, 1, xm);		// north imag
			vstor(eori, 2, yp);		// south real
			vstor(eori, 3, ym);		// south imag
			eori += 4*VSIZE2;
		}
		return;
	}
    for (long k=k0v; k<nk; k+=2) {
        v2d a = vread2(Fm + k*k_inc, 0);
        v2d b = vread2(Fm2 + k*k_inc, 0);
        v2d xm = vxchg(a-b);		v2d xp = a+b;
		xm = vxor2(xm, *(v2d*)_neg0);	// change sign of even values in vector

		v2d c = vread2(Fm  + (nlat-2-k)*k_inc, 0);
		v2d d = vread2(Fm2 + (nlat-2-k)*k_inc, 0);
		v2d yp = vxchg(c+d);		v2d ym = d-c;
		ym = vxor2(ym, *(v2d*)_neg0); 	// change sign of even values in vector

		unsigned long kk = (((unsigned)k) % VSIZE2) + 4*VSIZE2*(((unsigned)k)/VSIZE2);
		*(v2d*)(eori+kk) = xp;
		*(v2d*)(eori+kk + VSIZE2) = xm;
		*(v2d*)(eori+kk + 2*VSIZE2) = yp;
		*(v2d*)(eori+kk + 3*VSIZE2) = ym;
    }
  #else
    for (long k=k0v; k<nk; k+=2) {
        double an, bn, ani, bni, bs, as, bsi, asi, t;
        ani = Fm[k*k_inc];	bni = Fm[k*k_inc+1];
        an  = Fm2[k*k_inc];	bn = Fm2[k*k_inc+1];
        t = ani-an;	an += ani;		ani = bn-bni;		bn += bni;		bni = t;
        bsi = Fm[(nlat-2-k)*k_inc];	asi = Fm[(nlat-2-k)*k_inc +1];
        bs = Fm2[(nlat-2-k)*k_inc];	as = Fm2[(nlat-2-k)*k_inc +1];
        t = bsi-bs;		bs += bsi;		bsi = as-asi;		as += asi;		asi = t;
		#if VSIZE2 >= 2
			unsigned long kk = (((unsigned)k) % VSIZE2) + 4*VSIZE2*(((unsigned)k)/VSIZE2);
			eori[kk] = an;			eori[kk+1] = bn;			// north real
			eori[kk+VSIZE2] = ani;	eori[kk+VSIZE2 +1] = bni;		// north imag
			eori[kk+2*VSIZE2] = as;		eori[kk+2*VSIZE2 +1] = bs;		// south real
			eori[kk+3*VSIZE2] = asi;	eori[kk+3*VSIZE2 +1] = bsi;		// south imag
		#else
			unsigned long kk = 4*k;
			eori[kk] = an;		eori[kk+4] = bn;		// north real
			eori[kk+1] = ani;	eori[kk+5] = bni;		// north imag
			eori[kk+2] = as;	eori[kk+6] = bs;		// south real
			eori[kk+3] = asi;	eori[kk+7] = bsi;		// south imag
		#endif
    }
  #endif
}

// compute symmetric and antisymmetric parts, and reorganize data.
static
void split_sym_asym_m0(double* F0, double* eo, unsigned nlat_2, int k_inc)
{
	unsigned nk = ((nlat_2 +(VSIZE2-1))/VSIZE2)*VSIZE2;
	long int k=0;
	#if VSIZE2 >= 2
	if LIKELY(k_inc == 1) {		// optimized and vectorized for k_inc==1
		do {
			rnd an = vread(F0 + k, 0);
			rnd as = vread(F0 + 2*nlat_2-k-VSIZE2, 0);
			as = vreverse(as);
			vstor(eo+2*k, 0, an+as);
			vstor(eo+2*k, 1, an-as);
			k+=VSIZE2;
		} while(k < nk);
		return;
	}
	do {
		v2d an[VSIZE2/2];	v2d as[VSIZE2/2];
		for (int j=0; j<VSIZE2/2; j++) {
			an[j] = vread2(F0 + (k+2*j)*k_inc, 0);
			as[j] = vread2(F0 + (2*nlat_2-2-(k+2*j))*k_inc, 0);
			as[j] = vxchg(as[j]);
			*(v2d*)(eo+2*(k+j)) = an[j]+as[j];
			*(v2d*)(eo+2*(k+j) + VSIZE2) = an[j]-as[j];
		}
		k+=VSIZE2;
	} while(k < nk);
	#else
	do {
		double an = F0[k*k_inc];				double bn = F0[k*k_inc +1];
		double as = F0[(2*nlat_2-2-k)*k_inc +1];	double bs = F0[(2*nlat_2-2-k)*k_inc];
		unsigned long kk = (((unsigned)k) % VSIZE2) + 2*VSIZE2*(((unsigned)k)/VSIZE2);
		eo[kk] = an+as;		eo[kk + VSIZE2] = an-as;
		eo[kk+1+(VSIZE2==1)] = bn+bs;		eo[kk+1 + VSIZE2 + (VSIZE2==1)] = bn-bs;
		k+=2;
	} while(k < nk);
	#endif 
}

static
double split_sym_asym_m0_accl0(double* F0, double* eo, unsigned nlat_2, int k_inc, double* wg)
{
	unsigned nk = ((nlat_2 +(VSIZE2-1))/VSIZE2)*VSIZE2;
	long int k=0;
  #if VSIZE2 >= 2
	if LIKELY(k_inc == 1) {		// optimized and vectorized for k_inc==1
		rnd r0 = vall(0.0);
		do {
			rnd an = vread(F0 + k, 0);
			rnd as = vread(F0 + 2*nlat_2-k-VSIZE2, 0);
			as = vreverse(as);
			rnd ev = an+as;
			rnd od = an-as;
			vstor(eo+2*k, 0, ev);
			vstor(eo+2*k, 1, od);
			r0 += ev * vread(wg+k, 0);
			k+=VSIZE2;
		} while(k < nk);
		return reduce_add(r0);
	}
	v2d r0[VSIZE2/2];
	for (int j=0; j<VSIZE2/2; j++) r0[j] = vdup(0.0);	// independent accumulators
	do {
		v2d an[VSIZE2/2];	v2d as[VSIZE2/2];
		for (int j=0; j<VSIZE2/2; j++) {
			an[j] = vread2(F0 + (k+2*j)*k_inc, 0);
			as[j] = vread2(F0 + (2*nlat_2-2-(k+2*j))*k_inc, 0);
			as[j] = vxchg(as[j]);
			*(v2d*)(eo+2*(k+j)) = an[j]+as[j];
			*(v2d*)(eo+2*(k+j) + VSIZE2) = an[j]-as[j];
			r0[j] += (an[j]+as[j]) * vread2(wg+k, j);
		}
		k+=VSIZE2;
	} while(k < nk);
	#if VSIZE2 >= 4
	for (int j=2; j<VSIZE2/2; j+=2) {	r0[0] += r0[j];		r0[1] += r0[j+1];	}
	r0[0] += r0[1];
	#endif
	return vlo_to_dbl(r0[0]) + vhi_to_dbl(r0[0]);
  #else
	double r0a = 0.0;	double r0b = 0.0;	// two independent accumulators
	do {
		double an = F0[k*k_inc];				double bn = F0[k*k_inc +1];
		double as = F0[(2*nlat_2-2-k)*k_inc +1];	double bs = F0[(2*nlat_2-2-k)*k_inc];
		unsigned long kk = (((unsigned)k) % VSIZE2) + 2*VSIZE2*(((unsigned)k)/VSIZE2);
		eo[kk] = an+as;						eo[kk + VSIZE2] = an-as;
		eo[kk+1+(VSIZE2==1)] = bn+bs;		eo[kk+1 + VSIZE2 + (VSIZE2==1)] = bn-bs;
		r0a += (an+as)*wg[k];	r0b += (bn+bs)*wg[k+1];
		k+=2;
	} while(k < nk);
	return r0a+r0b;
  #endif
}

#else /* SHTNS4MAGIC */

// Fm = F + im*m_inc
// Fm2 = F + (NPHI-im)*m_inc
static void split_north_south_real_imag(double* Fm, double* Fm2, double* eori, long k0v, unsigned nlat, int k_inc)
{
    const unsigned nk = (nlat+1)>>1;
	for (long k = ((k0v*VSIZE2)>>1)*2; k<nk; k++) {
		double ar,ai,br,bi, sr,si,nr,ni;
		br = Fm[2*k*k_inc];		bi = Fm[2*k*k_inc +1];
		ar = Fm2[2*k*k_inc];	ai = Fm2[2*k*k_inc +1];
		nr = ar + br;		ni = ai - bi;
		sr = ai + bi;		si = br - ar;
		unsigned long kk = (((unsigned)k) % VSIZE2) + 4*VSIZE2*(((unsigned)k)/VSIZE2);
		eori[kk] = nr;				eori[kk +VSIZE2] = ni;
		eori[kk +2*VSIZE2] = sr;	eori[kk +3*VSIZE2] = si;
	}
	// here, it is assumed nk is a multiple of vector size.
}

// compute symmetric and antisymmetric parts, and reorganize data.
static
void split_sym_asym_m0(double* F0, double* eo, unsigned nlat_2, int k_inc)
{
	long k=0; do {
		double an = F0[2*k*k_inc];		double as = F0[2*k*k_inc +1];
		unsigned long kk = (((unsigned)k) % VSIZE2) + 2*VSIZE2*(((unsigned)k)/VSIZE2);
		eo[kk] = an+as;			eo[kk +VSIZE2] = an-as;
		k+=1;
	} while(k < nlat_2);
	// here, it is assumed nlat_2 is a multiple of vector size.
}

static
double split_sym_asym_m0_accl0(double* F0, double* eo, unsigned nlat_2, int k_inc, double* wg)
{
	double acc0 = 0.0;
	long k=0; do {
		double an = F0[2*k*k_inc];	double as = F0[2*k*k_inc +1];
		unsigned long kk = (((unsigned)k) % VSIZE2) + 2*VSIZE2*(((unsigned)k)/VSIZE2);
		eo[kk] = (an+as);			eo[kk +VSIZE2] = (an-as);
		acc0 += (an+as)*wg[k];
		k+=1;
	} while(k < nlat_2);
	// here, it is assumed nlat_2 is a multiple of vector size.
	return acc0;
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

	#define NWAY 1
	#include "SHT/spat_to_SH_kernel.c"
	#include "SHT/spat_to_SHst_kernel.c"
	#undef NWAY
	#define NWAY 2
	#include "SHT/spat_to_SH_kernel.c"
	#include "SHT/spat_to_SHst_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/spat_to_SH_kernel.c"
	#include "SHT/spat_to_SHst_kernel.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/spat_to_SH_kernel.c"
	#include "SHT/spat_to_SHst_kernel.c"
	#undef NWAY
  #if VSIZE2 <= 4
	#define NWAY 6
	#include "SHT/spat_to_SH_kernel.c"
	#include "SHT/spat_to_SHst_kernel.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/spat_to_SH_kernel.c"
	#include "SHT/spat_to_SHst_kernel.c"
	#undef NWAY
  #endif

#define SHT_3COMP
	#define NWAY 1
	#include "SHT/spat_to_SHqst_kernel.c"
	#undef NWAY
	#define NWAY 2
	#include "SHT/spat_to_SHqst_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/spat_to_SHqst_kernel.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/spat_to_SHqst_kernel.c"
	#undef NWAY
  #if VSIZE2 <= 4
	#define NWAY 6
	#include "SHT/spat_to_SHqst_kernel.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/spat_to_SHqst_kernel.c"
	#undef NWAY
  #endif
#undef SHT_3COMP

// genaral case, low lmax
#undef SUFFIX
#define SUFFIX _l
#undef HI_LLIM

	#define NWAY 1
	#include "SHT/spat_to_SH_kernel.c"
	#include "SHT/spat_to_SHst_kernel.c"
	#undef NWAY
	#define NWAY 2
	#include "SHT/spat_to_SH_kernel.c"
	#include "SHT/spat_to_SHst_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/spat_to_SH_kernel.c"
	#include "SHT/spat_to_SHst_kernel.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/spat_to_SH_kernel.c"
	#include "SHT/spat_to_SHst_kernel.c"
	#undef NWAY
  #if VSIZE2 <= 4
	#define NWAY 6
	#include "SHT/spat_to_SH_kernel.c"
	#include "SHT/spat_to_SHst_kernel.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/spat_to_SH_kernel.c"
	#include "SHT/spat_to_SHst_kernel.c"
	#undef NWAY
  #endif

#define SHT_3COMP
	#define NWAY 1
	#include "SHT/spat_to_SHqst_kernel.c"
	#undef NWAY
	#define NWAY 2
	#include "SHT/spat_to_SHqst_kernel.c"
	#undef NWAY
	#define NWAY 3
	#include "SHT/spat_to_SHqst_kernel.c"
	#undef NWAY
	#define NWAY 4
	#include "SHT/spat_to_SHqst_kernel.c"
	#undef NWAY
  #if VSIZE2 <= 4
	#define NWAY 6
	#include "SHT/spat_to_SHqst_kernel.c"
	#undef NWAY
	#define NWAY 8
	#include "SHT/spat_to_SHqst_kernel.c"
	#undef NWAY
  #endif
#undef SHT_3COMP

