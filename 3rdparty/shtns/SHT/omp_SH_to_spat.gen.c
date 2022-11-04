/*
 * Copyright (c) 2010-2019 Centre National de la Recherche Scientifique.
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
# it is intended for "make" to generate C code for similar SHT functions,
# from one generic function + tags.
# > See Makefile and SHT.c
# Basically, there are tags at the beginning of lines that are information
# to keep or remove the line depending on the function to build.
# tags :
# Q : line for scalar transform
# V : line for vector transform (both spheroidal and toroidal)
# S : line for vector transfrom, spheroidal component
# T : line for vector transform, toroidal component.

3	void GEN3(_sy3,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, v2d *BrF, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1);
QX	void GEN3(_sy1,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Qlm, v2d *BrF, long int llim, const unsigned im, int it0, int it1);
  #ifndef SHT_GRAD
VX	void GEN3(_sy2,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, cplx *Tlm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1);
  #else
S	void GEN3(_sy1s,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1);
T	void GEN3(_sy1t,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Tlm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1);
  #endif

3	void GEN3(_sy3_hi,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, v2d *BrF, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1);
QX	void GEN3(_sy1_hi,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Qlm, v2d *BrF, long int llim, const unsigned im, int it0, int it1);
  #ifndef SHT_GRAD
VX	void GEN3(_sy2_hi,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, cplx *Tlm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1);
  #else
S	void GEN3(_sy1s_hi,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1);
T	void GEN3(_sy1t_hi,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Tlm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1);
  #endif




3	static void GEN3(SHqst_to_spat_omp_a,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp, long int llim) {
QX	static void GEN3(SH_to_spat_omp_a,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Qlm, double *Vr, long int llim) {
  #ifndef SHT_GRAD
VX	static void GEN3(SHsphtor_to_spat_omp_a,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, cplx *Tlm, double *Vt, double *Vp, long int llim) {
  #else
S	static void GEN3(SHsph_to_spat_omp_a,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, double *Vt, double *Vp, long int llim) {
T	static void GEN3(SHtor_to_spat_omp_a,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Tlm, double *Vt, double *Vp, long int llim) {
  #endif

	unsigned imlim = 0;
Q	v2d* BrF = (v2d*) Vr;
V	v2d* BtF = (v2d*) Vt;	v2d* BpF = (v2d*) Vp;

  #ifndef SHT_AXISYM
	imlim = MTR;
	#ifdef SHT_VAR_LTR
		if (imlim*MRES > (unsigned) llim) imlim = ((unsigned) llim)/MRES;		// 32bit mul and div should be faster
	#endif
	if (shtns->fftc_mode > 0) {		// alloc memory for the FFT
		unsigned long nv = shtns->nspat;
QX		BrF = (v2d*) VMALLOC( nv * sizeof(double) );
VX		BtF = (v2d*) VMALLOC( 2*nv * sizeof(double) );
VX		BpF = BtF + nv/2;
3		BrF = (v2d*) VMALLOC( 3*nv * sizeof(double) );
3		BtF = BrF + nv/2;		BpF = BrF + nv;
	}
  #endif

  #pragma omp parallel num_threads(shtns->nthreads)
  {
	const int it0=0;
	const int it1=NLAT_2;
	#pragma omp for schedule(static,1) nowait
	for (int im=0; im<=imlim; im++)
	{
3		GEN3(_sy3_hi,NWAY,SUFFIX)(shtns, Qlm, Slm, Tlm, BrF, BtF, BpF, llim, im, it0, it1);
QX		GEN3(_sy1_hi,NWAY,SUFFIX)(shtns, Qlm, BrF, llim, im, it0, it1);
	#ifndef SHT_GRAD
VX		GEN3(_sy2_hi,NWAY,SUFFIX)(shtns, Slm, Tlm, BtF, BpF, llim, im, it0, it1);
	#else
S		GEN3(_sy1s_hi,NWAY,SUFFIX)(shtns, Slm, BtF, BpF, llim, im, it0, it1);
T		GEN3(_sy1t_hi,NWAY,SUFFIX)(shtns, Tlm, BtF, BpF, llim, im, it0, it1);
	#endif
	}

  #ifndef SHT_AXISYM
	// padding for high m's
	if (NPHI-1 > 2*imlim) {
		const int m_inc = shtns->nlat_padded >> 1;
		#pragma omp for schedule(static) nowait
		for (int im=imlim+1; im < NPHI-imlim; im++)  {
Q			memset(BrF + m_inc*im, 0, sizeof(cplx)* m_inc );
V			memset(BtF + m_inc*im, 0, sizeof(cplx)* m_inc );
V			memset(BpF + m_inc*im, 0, sizeof(cplx)* m_inc );
		}
	}
  #endif
  }

  #ifndef SHT_AXISYM
    // NPHI > 1 as SHT_AXISYM is not defined.
	if (shtns->fftc_mode >= 0) {
		if (shtns->fftc_mode != 1) {
Q			fftw_execute_dft(shtns->ifftc, ((cplx *) BrF), ((cplx *) Vr));
V			fftw_execute_dft(shtns->ifftc, ((cplx *) BtF), ((cplx *) Vt));
V			fftw_execute_dft(shtns->ifftc, ((cplx *) BpF), ((cplx *) Vp));
		} else {		// split dft
Q			fftw_execute_split_dft(shtns->ifftc,((double*)BrF)+1, ((double*)BrF), Vr+NPHI, Vr);
V			fftw_execute_split_dft(shtns->ifftc,((double*)BtF)+1, ((double*)BtF), Vt+NPHI, Vt);
V			fftw_execute_split_dft(shtns->ifftc,((double*)BpF)+1, ((double*)BpF), Vp+NPHI, Vp);
Q			VFREE(BrF);
VX			VFREE(BtF);		// this frees also BpF.
		}
	}
  #endif
  }


3	static void GEN3(SHqst_to_spat_omp_b,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp, long int llim) {
QX	static void GEN3(SH_to_spat_omp_b,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Qlm, double *Vr, long int llim) {
  #ifndef SHT_GRAD
VX	static void GEN3(SHsphtor_to_spat_omp_b,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, cplx *Tlm, double *Vt, double *Vp, long int llim) {
  #else
S	static void GEN3(SHsph_to_spat_omp_b,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, double *Vt, double *Vp, long int llim) {
T	static void GEN3(SHtor_to_spat_omp_b,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Tlm, double *Vt, double *Vp, long int llim) {
  #endif

	unsigned imlim = 0;
Q	v2d* BrF = (v2d*) Vr;
V	v2d* BtF = (v2d*) Vt;	v2d* BpF = (v2d*) Vp;

  #ifndef SHT_AXISYM
	imlim = MTR;
	#ifdef SHT_VAR_LTR
		if (imlim*MRES > (unsigned) llim) imlim = ((unsigned) llim)/MRES;		// 32bit mul and div should be faster
	#endif
	if (shtns->fftc_mode > 0) {		// alloc memory for the FFT
		unsigned long nv = shtns->nspat;
QX		BrF = (v2d*) VMALLOC( nv * sizeof(double) );
VX		BtF = (v2d*) VMALLOC( 2*nv * sizeof(double) );
VX		BpF = BtF + nv/2;
3		BrF = (v2d*) VMALLOC( 3*nv * sizeof(double) );
3		BtF = BrF + nv/2;		BpF = BrF + nv;
	}
  #endif
  
  #pragma omp parallel num_threads(shtns->nthreads)
  {
	  int it0 = 0;
	  int it1 = NLAT_2;
	  const int it_step = NLAT_2 / shtns->nthreads;
	#pragma omp for schedule(dynamic,1) collapse(1) nowait
	//for (int it0=0; it0<NLAT_2; it0 += it_step)
	for (int im=0; im <= imlim/2; im++)
	{
		//it1 = it0 + it_step;
3		GEN3(_sy3_hi,NWAY,SUFFIX)(shtns, Qlm, Slm, Tlm, BrF, BtF, BpF, llim, im, it0, it1);
QX		GEN3(_sy1_hi,NWAY,SUFFIX)(shtns, Qlm, BrF, llim, im, it0, it1);
Q		if (imlim-im > im) {
3			GEN3(_sy3_hi,NWAY,SUFFIX)(shtns, Qlm, Slm, Tlm, BrF, BtF, BpF, llim, imlim-im, it0, it1);
QX			GEN3(_sy1_hi,NWAY,SUFFIX)(shtns, Qlm, BrF, llim, imlim-im, it0, it1);
Q		}
V		#ifndef SHT_GRAD
VX		GEN3(_sy2_hi,NWAY,SUFFIX)(shtns, Slm, Tlm, BtF, BpF, llim, im, it0, it1);
VX		if (imlim-im > im) {
VX			GEN3(_sy2_hi,NWAY,SUFFIX)(shtns, Slm, Tlm, BtF, BpF, llim, imlim-im, it0, it1);
VX		}
V		#else
S		GEN3(_sy1s_hi,NWAY,SUFFIX)(shtns, Slm, BtF, BpF, llim, im, it0, it1);
T		GEN3(_sy1t_hi,NWAY,SUFFIX)(shtns, Tlm, BtF, BpF, llim, im, it0, it1);
V		if (imlim-im > im) {
S			GEN3(_sy1s_hi,NWAY,SUFFIX)(shtns, Slm, BtF, BpF, llim, imlim-im, it0, it1);
T			GEN3(_sy1t_hi,NWAY,SUFFIX)(shtns, Tlm, BtF, BpF, llim, imlim-im, it0, it1);		
V		}
V		#endif
	}

  #ifndef SHT_AXISYM
	// padding for high m's
	if (NPHI-1 > 2*imlim) {
		#pragma omp for schedule(dynamic) nowait
		for (int im=imlim+1; im < NPHI-imlim; im++)  {
Q			memset(BrF + NLAT_2*im, 0, sizeof(cplx)* NLAT_2 );
V			memset(BtF + NLAT_2*im, 0, sizeof(cplx)* NLAT_2 );
V			memset(BpF + NLAT_2*im, 0, sizeof(cplx)* NLAT_2 );
		}
	}

	if (shtns->fftc_mode >= 0) {
		const int nblk = (NLAT/2) / shtns->nthreads;
		#pragma omp barrier
		if (shtns->fftc_mode != 1) {
			for (int k=0; k<shtns->nthreads; k++) {
Q				#pragma omp single nowait
Q				fftw_execute_dft(shtns->ifftc_block, ((cplx *) BrF) + k*nblk, ((cplx *) Vr) + k*nblk);
V				#pragma omp single nowait
V				fftw_execute_dft(shtns->ifftc_block, ((cplx *) BtF) + k*nblk, ((cplx *) Vt) + k*nblk);
V				#pragma omp single nowait
V				fftw_execute_dft(shtns->ifftc_block, ((cplx *) BpF) + k*nblk, ((cplx *) Vp) + k*nblk);
			}
		} else {
			for (int k=0; k<shtns->nthreads; k++) {
Q				#pragma omp single nowait
Q				fftw_execute_split_dft(shtns->ifftc_block,((double*)BrF)+1 +2*k*nblk, ((double*)BrF) +2*k*nblk, Vr+NPHI*(1+2*k*nblk), Vr +NPHI*2*k*nblk);
V				#pragma omp single nowait
V				fftw_execute_split_dft(shtns->ifftc_block,((double*)BtF)+1 +2*k*nblk, ((double*)BtF) +2*k*nblk, Vt+NPHI*(1+2*k*nblk), Vt +NPHI*2*k*nblk);
V				#pragma omp single nowait
V				fftw_execute_split_dft(shtns->ifftc_block,((double*)BpF)+1 +2*k*nblk, ((double*)BpF) +2*k*nblk, Vp+NPHI*(1+2*k*nblk), Vp +NPHI*2*k*nblk);
			}
		}
	}
  #endif
  }

  #ifndef SHT_AXISYM
	if (shtns->fftc_mode > 0) {
Q			VFREE(BrF);
VX			VFREE(BtF);		// this frees also BpF.
	}
  #endif

  }
