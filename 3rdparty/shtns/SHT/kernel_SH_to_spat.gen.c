/*
 * Copyright (c) 2010-2021 Centre National de la Recherche Scientifique.
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

# This file is meta-code for SHT.c (spherical harmonic transform).
# it is intended for "make" to generate C code for 3 similar SHT functions,
# (namely spat_to_SH [Q tag]), spat_to_SHsphtor [V tag], spat_to_SH3 [both Q&V tags])
# from one generic function + tags.
# Basically, there are tags at the beginning of lines (Q,V) that are information
# to keep or remove the line depending on the function to build. (Q for scalar, V for vector, # for comment)
#
//////////////////////////////////////////////////

	#ifdef HI_LLIM
QX	#define BASE _sy1_hi
3	#define BASE _sy3_hi
	#ifndef SHT_GRAD
VX	#define BASE _sy2_hi
	#else
S	#define BASE _sy1s_hi
T	#define BASE _sy1t_hi
	#endif
	#else
QX	#define BASE _sy1
3	#define BASE _sy3
	#ifndef SHT_GRAD
VX	#define BASE _sy2
	#else
S	#define BASE _sy1s
T	#define BASE _sy1t
	#endif
	#endif

3	void GEN3(BASE,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, v2d *BrF, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1) {
QX	void GEN3(BASE,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Qlm, v2d *BrF, long int llim, const unsigned im, int it0, int it1) {
  #ifndef SHT_GRAD
VX	void GEN3(BASE,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, cplx *Tlm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1) {
  #else
S	void GEN3(BASE,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1) {
T	void GEN3(BASE,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Tlm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1) {
  #endif

	#if !defined( _GCC_VEC_ ) && (NWAY & 1)
	#error "NWAY must be even when compiled without explicit vectorization."
	#endif
	#if VSIZE2*NWAY > 32
	#error "VSIZE2*NWAY must not exceed 32"
	#endif

  #ifndef SHT_AXISYM
   #ifndef SHTNS_ISHIOKA
Q	#define qr(l) vall(creal(Ql[l]))
Q	#define qi(l) vall(cimag(Ql[l]))
   #else
Q	#define qr(l) vall( ((double*) QQl)[2*(l)]   )
Q	#define qi(l) vall( ((double*) QQl)[2*(l)+1] )
   #endif
V	#define vr(l) vall( ((double*) VWl)[4*(l)]   )
V	#define vi(l) vall( ((double*) VWl)[4*(l)+1] )
V	#define wr(l) vall( ((double*) VWl)[4*(l)+2] )
V	#define wi(l) vall( ((double*) VWl)[4*(l)+3] )
  #endif
	long int nk,k,l,m;
	double *alm, *al;
	double *ct, *st;
V	int robert_form;
QX	double Ql0[llim+2];
V	v2d VWl[llim*2+4];
  #ifdef SHTNS_ISHIOKA
Q	v2d QQl[llim+2];
  #endif

	ct = shtns->ct;		st = shtns->st;
	nk = it1;	//NLAT_2;
	#if _GCC_VEC_
		nk = ((unsigned)(nk+VSIZE2-1)) / VSIZE2;
		it0 = ((unsigned)(it0+VSIZE2-1)) / VSIZE2;
	#endif
V	robert_form = shtns->robert_form;

	if (im == 0)
	{	//	im=0;
S		double* const Sl0 = (double*) VWl;
T		double* const Tl0 = (double*) VWl + llim+2;
3		double* const Ql0 = (double*) (VWl + llim+2);
		#ifdef SHT_GRAD
			// TODO: fix k,nk bounds
S			if (BpF != NULL) memset(BpF, 0, sizeof(v2d) * NLAT_2);
T			if (BtF != NULL) memset(BtF, 0, sizeof(v2d) * NLAT_2);
		#endif
 		l=1;
		alm = shtns->alm;
Q		Ql0[0] = (double) Qlm[0];		// l=0
		do {		// for m=0, compress the complex Q,S,T to double
Q			Ql0[l] = creal( Qlm[l] );	//	Ql[l+1] = (double) Qlm[l+1];
S			Sl0[l-1] = creal( Slm[l] );	//	Sl[l] = (double) Slm[l+1];
T			Tl0[l-1] = creal( Tlm[l] );	//	Tl[l] = (double) Tlm[l+1];
			++l;
		} while(l<=llim);
		k=it0;
		do {
			l=0;	al = alm;
			rnd cost[NWAY], y0[NWAY], y1[NWAY];
V			rnd sint[NWAY], dy0[NWAY], dy1[NWAY];
Q			rnd re[NWAY], ro[NWAY];
S			rnd te[NWAY], to[NWAY];
T			rnd pe[NWAY], po[NWAY];
			for (int j=0; j<NWAY; ++j) {
				cost[j] = vread(ct, j+k);
V				sint[j] = -vread(st, j+k);
				y0[j] = vall(al[0]);
V				dy0[j] = vall(0.0);
Q				re[j] = y0[j] * vall(Ql0[0]);
S				to[j] = dy0[j];
T				po[j] = dy0[j];
			}
V			if (robert_form) {
V				for (int j=0; j<NWAY; ++j) sint[j] *= -sint[j];
V			}
			for (int j=0; j<NWAY; ++j) {
				y1[j]  = vall(al[0]*al[1]) * cost[j];
V				dy1[j] = vall(al[0]*al[1]) * sint[j];
			}
			for (int j=0; j<NWAY; ++j) {
Q				ro[j] = y1[j] * vall(Ql0[1]);
S				te[j] = dy1[j] * vall(Sl0[0]);
T				pe[j] = -dy1[j] * vall(Tl0[0]);
			}
			al+=2;	l+=2;
			while(l<llim) {
				for (int j=0; j<NWAY; ++j) {
V					dy0[j] = vall(al[1])*(cost[j]*dy1[j] + y1[j]*sint[j]) + vall(al[0])*dy0[j];
					y0[j]  = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; ++j) {
Q					re[j] += y0[j] * vall(Ql0[l]);
S					to[j] += dy0[j] * vall(Sl0[l-1]);
T					po[j] -= dy0[j] * vall(Tl0[l-1]);
				}
				for (int j=0; j<NWAY; ++j) {
V					dy1[j] = vall(al[3])*(cost[j]*dy0[j] + y0[j]*sint[j]) + vall(al[2])*dy1[j];
					y1[j]  = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
				}
				for (int j=0; j<NWAY; ++j) {
Q					ro[j] += y1[j] * vall(Ql0[l+1]);
S					te[j] += dy1[j] * vall(Sl0[l]);
T					pe[j] -= dy1[j] * vall(Tl0[l]);
				}
				al+=4;	l+=2;
			}
			if (l==llim) {
				for (int j=0; j<NWAY; ++j) {
V					dy0[j] = vall(al[1])*(cost[j]*dy1[j] + y1[j]*sint[j]) + vall(al[0])*dy0[j];
					y0[j]  = vall(al[1])*cost[j]*y1[j] + vall(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; ++j) {
Q					re[j] += y0[j] * vall(Ql0[l]);
S					to[j] += dy0[j] * vall(Sl0[l-1]);
T					po[j] -= dy0[j] * vall(Tl0[l-1]);
				}
			}
			// combine even/odd into north/south
Q			for (int j=0; j<NWAY; ++j) {
Q				rnd s = re[j] - ro[j];		re[j] = re[j] + ro[j];
Q				ro[j] = s;
Q			}
V			for (int j=0; j<NWAY; ++j) {			
S				rnd ts = te[j] - to[j];	te[j] = te[j] + to[j];
T				rnd ps = pe[j] - po[j];	pe[j] = pe[j] + po[j];
S				to[j] = ts;
T				po[j] = ps;
V			}
		#ifndef SHTNS4MAGIC
			for (int j=0; j<NWAY; ++j) {
Q				S2D_STORE(BrF, j+k, re[j], ro[j])
S				S2D_STORE(BtF, j+k, te[j], to[j])
T				S2D_STORE(BpF, j+k, pe[j], po[j])
			}
		#else
			for (int j=0; j<NWAY; ++j) {
				if ((k+j)>=nk) break;
Q				S2D_STORE_4MAGIC((double*)BrF, j+k, re[j], ro[j]);
S				S2D_STORE_4MAGIC((double*)BtF, j+k, te[j], to[j]);
T				S2D_STORE_4MAGIC((double*)BpF, j+k, pe[j], po[j]);
			}
		#endif
			k+=NWAY;
		} while (k < nk);
	}
  #ifndef SHT_AXISYM
	else
	{		// im > 0
Q		BrF += im*(shtns->nlat_padded >>1);
V		BtF += im*(shtns->nlat_padded >>1);
V		BpF += im*(shtns->nlat_padded >>1);
		m = im*MRES;
		l = (im*(2*(LMAX+1)-(m+MRES)))>>1;		//l = LiM(shtns, 0,im);
		#ifndef SHTNS_ISHIOKA
		alm = shtns->alm + 2*(l+m);		// shtns->alm + im*(2*(LMAX+1) -m+MRES);
		#else
		alm = shtns->clm + (l+m);		// shtns->clm + im*(2*(LMAX+1) -m+MRES)/2;
		#endif

  #ifndef SHT_GRAD
V		SH_vect_to_2scal(shtns->mx_stdt + 2*l, llim, m, &Slm[l], &Tlm[l], (cplx*) VWl);
  #else
S		SHsph_to_2scal(shtns->mx_stdt + 2*l, llim, m, &Slm[l], (cplx*) VWl);
T		SHtor_to_2scal(shtns->mx_stdt + 2*l, llim, m, &Tlm[l], (cplx*) VWl);
  #endif

	#ifndef SHTNS_ISHIOKA
Q		cplx* Ql = &Qlm[l];	// virtual pointer for l=0 and im
	#else
		// pre-processing for recurrence relation of Ishioka
		const double* restrict xlm = shtns->xlm + 3*im*(2*(LMAX+4) -m+MRES)/4;
Q		v2d* Ql = (v2d*) &Qlm[l];	// virtual pointer for l=0 and im
Q		SH_to_ishioka(xlm, Ql+m, llim-m, QQl+m);
V		SH2_to_ishioka(xlm, VWl+2*m, llim-m+1);
	#endif

		// polar optimization
		k = shtns->tm[im];		// start index in theta (=0 if no polar optimization)
		#if _GCC_VEC_
		k = ((unsigned) k) / VSIZE2;	// in vector size units.
		#else
		k = (k>>1)*2;		// k must be even.
		#endif
		if (it0 < k) {
			const long ofsm = (NPHI-2*im)*(shtns->nlat_padded >>1);
		#ifndef SHTNS4MAGIC
			#if _GCC_VEC_
			const long ofs1 = NLAT_2 - k*(VSIZE2/2);
			#else
			const long ofs1 = NLAT_2 - k/2;
			#endif
Q			zero_poles4_vect(BrF+it0*(VSIZE2/2), ofsm, ofs1, k-it0);
V			zero_poles4_vect(BtF+it0*(VSIZE2/2), ofsm, ofs1, k-it0);
V			zero_poles4_vect(BpF+it0*(VSIZE2/2), ofsm, ofs1, k-it0);
		#else
Q			zero_poles2_vect(BrF+it0*VSIZE2, ofsm, 2*(k-it0));
V			zero_poles2_vect(BtF+it0*VSIZE2, ofsm, 2*(k-it0));
V			zero_poles2_vect(BpF+it0*VSIZE2, ofsm, 2*(k-it0));
		#endif
		} else k=it0;

		do {
			rnd cost[NWAY], y0[NWAY], y1[NWAY];
Q			rnd rer[NWAY], rei[NWAY], ror[NWAY], roi[NWAY];
V			rnd ter[NWAY], tei[NWAY], tor[NWAY], toi[NWAY];
V			rnd per[NWAY], pei[NWAY], por[NWAY], poi[NWAY];
			for (int j=0; j<NWAY; ++j) {
				cost[j] = vread(st, k+j);
				y0[j] = vall(1.0);
			}
			l=m;
V			if (robert_form == 0) l=m-1;
	#ifndef HI_LLIM
			const long int ny = 0;
			while(1) {		// sin(theta)^m
				if (l&1) for (int j=0; j<NWAY; ++j) y0[j] *= cost[j];
				l >>= 1;
				if (l==0) break;
				for (int j=0; j<NWAY; ++j) cost[j] *= cost[j];
			}
	#else
			long int ny = 0;
			for (long int nsint = 0;;) {		// sin(theta)^m		(use rescaling to avoid underflow)
				if (l&1) {
					for (int j=NWAY-1; j>=0; --j) y0[j] *= cost[j];
					ny += nsint;
					if (vlo(y0[NWAY-1]) < (SHT_ACCURACY+1.0/SHT_SCALE_FACTOR)) {
						ny--;
						for (int j=NWAY-1; j>=0; --j) y0[j] *= vall(SHT_SCALE_FACTOR);
					}
				}
				l >>= 1;
				if (l==0) break;
				for (int j=NWAY-1; j>=0; --j) cost[j] *= cost[j];
				nsint += nsint;
				if (vlo(cost[NWAY-1]) < 1.0/SHT_SCALE_FACTOR) {
					nsint--;
					for (int j=NWAY-1; j>=0; --j) cost[j] *= vall(SHT_SCALE_FACTOR);
				}
			}
	#endif
			al = alm;
			for (int j=0; j<NWAY; ++j) {
				cost[j] = vread(ct, j+k);
Q				ror[j] = vall(0.0);		roi[j] = vall(0.0);
Q				rer[j] = vall(0.0);		rei[j] = vall(0.0);
				#ifndef SHTNS_ISHIOKA
				y0[j] *= vall(al[0]);
				#else
				cost[j] *= cost[j];		// cos(theta)^2
				#endif
			}
			for (int j=0; j<NWAY; ++j) {
				#ifndef SHTNS_ISHIOKA
				y1[j]  = (vall(al[1])*y0[j]) *cost[j];		//	y1[j] = vall(al[1])*cost[j]*y0[j];
				#else
				y1[j] = (vall(al[1])*cost[j] + vall(al[0]))*y0[j];
				#endif
V				por[j] = vall(0.0);		tei[j] = vall(0.0);
V				tor[j] = vall(0.0);		pei[j] = vall(0.0);
V				poi[j] = vall(0.0);		ter[j] = vall(0.0);
V				toi[j] = vall(0.0);		per[j] = vall(0.0);
			}
			l=m;		al+=2;
	#ifdef HI_LLIM
		  if (ny < 0) {		// ylm treated as zero and ignored if ny < 0
			const rnd scale = vall(1.0/SHT_SCALE_FACTOR);
			while (l<llim) {
				#ifndef SHTNS_ISHIOKA
				for (int j=0; j<NWAY; ++j)	y0[j] = (vall(al[1])*cost[j])*y1[j] + vall(al[0])*y0[j];
				for (int j=0; j<NWAY; ++j)	y1[j] = (vall(al[3])*cost[j])*y0[j] + vall(al[2])*y1[j];
				al+=4;
				#else
				rnd a[NWAY];
				for (int j=0; j<NWAY; ++j)	a[j] = vall(al[1])*cost[j] + vall(al[0]);
				al+=2;
				for (int j=0; j<NWAY; ++j) {
					a[j] = a[j]*y1[j] + y0[j];
					y0[j] = y1[j];		y1[j] = a[j];
				}
				#endif
				l+=2;
				if (fabs(vlo(y0[NWAY-1])) > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0) {		// rescale when value is significant
					for (int j=0; j<NWAY; ++j) {
						y0[j] *= scale;		y1[j] *= scale;
					}
					if (++ny == 0) break;
				}
			}
		  }
	#endif
		  if LIKELY(ny == 0) {
		#ifndef SHTNS_ISHIOKA
			while (l<llim) {	// compute even and odd parts
Q				for (int j=0; j<NWAY; ++j) {	rer[j] += y0[j]  * qr(l);		rei[j] += y0[j] * qi(l);	}
V				for (int j=0; j<NWAY; ++j) {	ter[j] += y0[j]  * vr(l);		tei[j] += y0[j] * vi(l);	}
V				for (int j=0; j<NWAY; ++j) {	per[j] += y0[j]  * wr(l);		pei[j] += y0[j] * wi(l);	}
				for (int j=0; j<NWAY; ++j) {
					y0[j] = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
				}
Q				for (int j=0; j<NWAY; ++j) {	ror[j] += y1[j]  * qr(l+1);		roi[j] += y1[j] * qi(l+1);	}
V				for (int j=0; j<NWAY; ++j) {	tor[j] += y1[j]  * vr(l+1);		toi[j] += y1[j] * vi(l+1);	}
V				for (int j=0; j<NWAY; ++j) {	por[j] += y1[j]  * wr(l+1);		poi[j] += y1[j] * wi(l+1);	}
				for (int j=0; j<NWAY; ++j) {
					y1[j] = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
				}
				l+=2;	al+=4;
			}
V				for (int j=0; j<NWAY; ++j) {	ter[j] += y0[j]  * vr(l);		tei[j] += y0[j] * vi(l);	}
V				for (int j=0; j<NWAY; ++j) {	per[j] += y0[j]  * wr(l);		pei[j] += y0[j] * wi(l);	}
			if (l==llim) {
Q				for (int j=0; j<NWAY; ++j) {	rer[j] += y0[j]  * qr(l);		rei[j] += y0[j] * qi(l);	}
V				for (int j=0; j<NWAY; ++j) {	tor[j] += y1[j]  * vr(l+1);		toi[j] += y1[j] * vi(l+1);	}
V				for (int j=0; j<NWAY; ++j) {	por[j] += y1[j]  * wr(l+1);		poi[j] += y1[j] * wi(l+1);	}
			}
			// combine even/odd into north/south
Q			for (int j=0; j<NWAY; ++j) {
Q				rnd sr = rer[j] - ror[j];	rer[j] = rer[j] + ror[j];
Q			  	rnd si = rei[j] - roi[j];	rei[j] = rei[j] + roi[j];
Q				ror[j] = sr;		roi[j] = si;
Q			}
V			for (int j=0; j<NWAY; ++j) {
V				rnd sr = ter[j] - tor[j];	ter[j] = ter[j] + tor[j];
V			  	rnd si = tei[j] - toi[j];	tei[j] = tei[j] + toi[j];
V				tor[j] = sr;		toi[j] = si;
V				sr = per[j] - por[j];	per[j] = per[j] + por[j];
V			  	si = pei[j] - poi[j];	pei[j] = pei[j] + poi[j];
V				por[j] = sr;		poi[j] = si;
V			}
		#else
			while (l<llim) {	// compute even and odd parts
QX				for (int j=0; j<NWAY; ++j) {	rer[j] += y0[j]  * qr(l);		rei[j] += y0[j] * qi(l);	}
QX				for (int j=0; j<NWAY; ++j) {	ror[j] += y0[j]  * qr(l+1);		roi[j] += y0[j] * qi(l+1);	}
3				for (int j=0; j<NWAY; ++j)	rer[j] += y0[j] * qr(l);
3				for (int j=0; j<NWAY; ++j)	rei[j] += y0[j] * qi(l);
3				for (int j=0; j<NWAY; ++j)	ror[j] += y0[j] * qr(l+1);
3				for (int j=0; j<NWAY; ++j)	roi[j] += y0[j] * qi(l+1);
V				for (int j=0; j<NWAY; ++j)	ter[j] += y0[j] * vr(l);
V				for (int j=0; j<NWAY; ++j)	tei[j] += y0[j] * vi(l);
V				for (int j=0; j<NWAY; ++j)	per[j] += y0[j] * wr(l);
V				for (int j=0; j<NWAY; ++j)	pei[j] += y0[j] * wi(l);
V				for (int j=0; j<NWAY; ++j)	tor[j] += y0[j] * vr(l+1);
V				for (int j=0; j<NWAY; ++j)	toi[j] += y0[j] * vi(l+1);
V				for (int j=0; j<NWAY; ++j)	por[j] += y0[j] * wr(l+1);
V				for (int j=0; j<NWAY; ++j)	poi[j] += y0[j] * wi(l+1);
				for (int j=0; j<NWAY; ++j) {
					rnd tmp = (vall(al[1])*cost[j] + vall(al[0]))*y1[j] + y0[j];
					y0[j] = y1[j];
					y1[j] = tmp;
				}
				l+=2;	al+=2;
			}
			for (int j=0; j<NWAY; ++j) cost[j] = vread(ct, k+j);		// read ahead to correct the odd part below
V			for (int j=0; j<NWAY; ++j)	ter[j] += y0[j] * vr(l);
V			for (int j=0; j<NWAY; ++j)	tei[j] += y0[j] * vi(l);
V			for (int j=0; j<NWAY; ++j)	per[j] += y0[j] * wr(l);
V			for (int j=0; j<NWAY; ++j)	pei[j] += y0[j] * wi(l);
			if LIKELY(l==llim) {
V				for (int j=0; j<NWAY; ++j)	tor[j] += y0[j] * vr(l+1);
V				for (int j=0; j<NWAY; ++j)	toi[j] += y0[j] * vi(l+1);
V				for (int j=0; j<NWAY; ++j)	por[j] += y0[j] * wr(l+1);
V				for (int j=0; j<NWAY; ++j)	poi[j] += y0[j] * wi(l+1);
Q				for (int j=0; j<NWAY; ++j)	rer[j] += y0[j] * qr(l);
Q				for (int j=0; j<NWAY; ++j)	rei[j] += y0[j] * qi(l);
			}
			// correct the odd part:
Q		//	for (int j=0; j<NWAY; ++j) {  ror[j] *= cost[j];	roi[j] *= cost[j]; }
V		//	for (int j=0; j<NWAY; ++j) {  tor[j] *= cost[j];	toi[j] *= cost[j]; }
V		//	for (int j=0; j<NWAY; ++j) {  por[j] *= cost[j];	poi[j] *= cost[j]; }

			// combine even/odd into north/south, and correct the odd part for free with FMA
Q			for (int j=0; j<NWAY; ++j) {
Q				rnd sr = rer[j] - ror[j]*cost[j];	rer[j] = rer[j] + ror[j]*cost[j];
Q			  	rnd si = rei[j] - roi[j]*cost[j];	rei[j] = rei[j] + roi[j]*cost[j];
Q				ror[j] = sr;		roi[j] = si;
Q			}
V			for (int j=0; j<NWAY; ++j) {
V				rnd sr = ter[j] - tor[j]*cost[j];	ter[j] = ter[j] + tor[j]*cost[j];
V			  	rnd si = tei[j] - toi[j]*cost[j];	tei[j] = tei[j] + toi[j]*cost[j];
V				tor[j] = sr;		toi[j] = si;
V				sr = per[j] - por[j]*cost[j];	per[j] = per[j] + por[j]*cost[j];
V			  	si = pei[j] - poi[j]*cost[j];	pei[j] = pei[j] + poi[j]*cost[j];
V				por[j] = sr;		poi[j] = si;
V			}
		#endif
3			if LIKELY(robert_form == 0) {
3				for (int j=0; j<NWAY; ++j) cost[j]  = vread(st, k+j);
3				for (int j=0; j<NWAY; ++j) {  rer[j] *= cost[j];  ror[j] *= cost[j];	rei[j] *= cost[j];  roi[j] *= cost[j];  }
3			}
		  }
		#ifndef SHTNS4MAGIC
		  #ifdef _GCC_VEC_
			const long ofs = (NPHI-2*im)*shtns->nlat_padded;
V			for (int j=0; j<NWAY; ++j) cstore_north_south((double*) BtF, ((double*) (BtF)) +ofs, k+j, NLAT, ter[j], tor[j], tei[j], toi[j]);
V			for (int j=0; j<NWAY; ++j) cstore_north_south((double*) BpF, ((double*) (BpF)) +ofs, k+j, NLAT, per[j], por[j], pei[j], poi[j]);
Q			for (int j=0; j<NWAY; ++j) cstore_north_south((double*) BrF, ((double*) (BrF)) +ofs, k+j, NLAT, rer[j], ror[j], rei[j], roi[j]);
		  #else
		  	// NWAY is even when _GCC_VEC_ is not defined
V			for (int j=0; j<NWAY/2; ++j) {	S2D_CSTOREX(BtF, k/2+j, 2*j, ter, tor, tei, toi)  }
V			for (int j=0; j<NWAY/2; ++j) {	S2D_CSTOREX(BpF, k/2+j, 2*j, per, por, pei, poi)  }
Q			for (int j=0; j<NWAY/2; ++j) {	S2D_CSTOREX(BrF, k/2+j, 2*j, rer, ror, rei, roi)  }
		  #endif
		#else
			for (int j=0; j<NWAY; ++j) {
				if ((k+j)>=nk) break;
V				S2D_CSTORE_4MAGIC((double*) BtF, (double*) (BtF + (NPHI-2*im)*(shtns->nlat_padded>>1)), k+j, ter[j], tor[j], tei[j], toi[j]);
V				S2D_CSTORE_4MAGIC((double*) BpF, (double*) (BpF + (NPHI-2*im)*(shtns->nlat_padded>>1)), k+j, per[j], por[j], pei[j], poi[j]);
Q				S2D_CSTORE_4MAGIC((double*) BrF, (double*) (BrF + (NPHI-2*im)*(shtns->nlat_padded>>1)), k+j, rer[j], ror[j], rei[j], roi[j]);
			}
		#endif
			k+=NWAY;
		} while (k < nk);
	}
  #endif
}

Q	#undef qr
Q	#undef qi
S	#undef sr
S	#undef si
T	#undef tr
T	#undef ti

	#undef BASE
