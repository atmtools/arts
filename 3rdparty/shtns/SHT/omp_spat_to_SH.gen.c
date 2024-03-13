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

# This file is meta-code for SHT.c (spherical harmonic transform).
# it is intended for "make" to generate C code for 3 similar SHT functions,
# (namely spat_to_SH [Q tag]), spat_to_SHsphtor [V tag], spat_to_SH3 [both Q&V tags])
# from one generic function + tags.
# Basically, there are tags at the beginning of lines (Q,V) that are information
# to keep or remove the line depending on the function to build. (Q for scalar, V for vector, # for comment)
#
//////////////////////////////////////////////////

QX	void GEN3(_an1,NWAY,SUFFIX)(shtns_cfg shtns, double *BrF, cplx *Qlm, const long int llim, const int im);
VX	void GEN3(_an2,NWAY,SUFFIX)(shtns_cfg shtns, double *BtF, double *BpF, cplx *Slm, cplx *Tlm, const long int llim, const int im);
3	void GEN3(_an3,NWAY,SUFFIX)(shtns_cfg shtns, double *BrF, double *BtF, double *BpF, cplx *Qlm, cplx *Slm, cplx *Tlm, const long int llim, const int im);
QX	void GEN3(_an1_hi,NWAY,SUFFIX)(shtns_cfg shtns, double *BrF, cplx *Qlm, const long int llim, const int im);
VX	void GEN3(_an2_hi,NWAY,SUFFIX)(shtns_cfg shtns, double *BtF, double *BpF, cplx *Slm, cplx *Tlm, const long int llim, const int im);
3	void GEN3(_an3_hi,NWAY,SUFFIX)(shtns_cfg shtns, double *BrF, double *BtF, double *BpF, cplx *Qlm, cplx *Slm, cplx *Tlm, const long int llim, const int im);


	static
QX	void GEN3(spat_to_SH_omp_a,NWAY,SUFFIX)(shtns_cfg shtns, double *Vr, cplx *Qlm, long int llim) {
VX	void GEN3(spat_to_SHsphtor_omp_a,NWAY,SUFFIX)(shtns_cfg shtns, double *Vt, double *Vp, cplx *Slm, cplx *Tlm, long int llim) {
3	void GEN3(spat_to_SHqst_omp_a,NWAY,SUFFIX)(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, long int llim) {

Q	double *BrF;		// contains the Fourier transformed data
V	double *BtF, *BpF;	// contains the Fourier transformed data
	unsigned imlim=0;

Q	BrF = Vr;
V	BtF = Vt;	BpF = Vp;
  #ifndef SHT_AXISYM
	imlim = MTR;
	#ifdef SHT_VAR_LTR
		if (imlim*MRES > (unsigned) llim) imlim = ((unsigned) llim)/MRES;		// 32bit mul and div should be faster
	#endif

	if (shtns->fftc_mode >= 0) {
		if (shtns->fftc_mode > 0) {		// alloc memory for out-of-place FFT
			unsigned long nv = shtns->nspat;
QX			BrF = (double*) VMALLOC( nv * sizeof(double) );
VX			BtF = (double*) VMALLOC( 2*nv * sizeof(double) );
VX			BpF = BtF + nv;
3			BrF = (double*) VMALLOC( 3*nv * sizeof(double) );
3			BtF = BrF + nv;		BpF = BtF + nv;
		}
		if (shtns->fftc_mode != 1) {
Q			fftw_execute_dft(shtns->fftc, ((cplx *) Vr), ((cplx *) BrF));
V			fftw_execute_dft(shtns->fftc, ((cplx *) Vt), ((cplx *) BtF));
V			fftw_execute_dft(shtns->fftc, ((cplx *) Vp), ((cplx *) BpF));
		} else {
Q			fftw_execute_split_dft(shtns->fftc, Vr+NPHI, Vr, BrF+1, BrF);
V			fftw_execute_split_dft(shtns->fftc, Vt+NPHI, Vt, BtF+1, BtF);
V			fftw_execute_split_dft(shtns->fftc, Vp+NPHI, Vp, BpF+1, BpF);
		}
	}
  #endif

	#pragma omp parallel num_threads(shtns->nthreads)
	{
		if (llim < SHT_L_RESCALE_FLY) {
			#pragma omp for schedule(static,1) nowait
			for (int im=0; im<=imlim; im++) {
QX				GEN3(_an1,NWAY,SUFFIX)(shtns, BrF, Qlm, llim, im);
VX				GEN3(_an2,NWAY,SUFFIX)(shtns, BtF, BpF, Slm, Tlm, llim, im);
3				GEN3(_an3,NWAY,SUFFIX)(shtns, BrF, BtF, BpF, Qlm, Slm, Tlm, llim, im);
			}
		} else {
			#pragma omp for schedule(static,1) nowait
			for (int im=0; im<=imlim; im++) {
QX				GEN3(_an1_hi,NWAY,SUFFIX)(shtns, BrF, Qlm, llim, im);
VX				GEN3(_an2_hi,NWAY,SUFFIX)(shtns, BtF, BpF, Slm, Tlm, llim, im);
3				GEN3(_an3_hi,NWAY,SUFFIX)(shtns, BrF, BtF, BpF, Qlm, Slm, Tlm, llim, im);
			}
		}
		#ifndef SHT_AXISYM
			if (imlim < MMAX) {		// zero out m > imlim
				long l = LiM(shtns, (imlim+1)*MRES, imlim+1);
Q				#pragma omp single nowait
Q				memset(Qlm+l, 0, (shtns->nlm - l)*sizeof(cplx));
V				#pragma omp single nowait
V				memset(Slm+l, 0, (shtns->nlm - l)*sizeof(cplx));
V				#pragma omp single nowait
V				memset(Tlm+l, 0, (shtns->nlm - l)*sizeof(cplx));
			}
		#endif
	}

  #ifndef SHT_AXISYM
  	if (shtns->fftc_mode > 0) {		// free memory
Q	    VFREE(BrF);
VX	    VFREE(BtF);	// this frees also BpF.
	}
  #endif

  }

	static
QX	void GEN3(spat_to_SH_omp_b,NWAY,SUFFIX)(shtns_cfg shtns, double *Vr, cplx *Qlm, long int llim) {
VX	void GEN3(spat_to_SHsphtor_omp_b,NWAY,SUFFIX)(shtns_cfg shtns, double *Vt, double *Vp, cplx *Slm, cplx *Tlm, long int llim) {
3	void GEN3(spat_to_SHqst_omp_b,NWAY,SUFFIX)(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, long int llim) {

Q	double *BrF;		// contains the Fourier transformed data
V	double *BtF, *BpF;	// contains the Fourier transformed data
	unsigned imlim=0;

Q	BrF = Vr;
V	BtF = Vt;	BpF = Vp;
  #ifndef SHT_AXISYM
	imlim = MTR;
	#ifdef SHT_VAR_LTR
		if (imlim*MRES > (unsigned) llim) imlim = ((unsigned) llim)/MRES;		// 32bit mul and div should be faster
	#endif

	if (shtns->fftc_mode > 0) {		// alloc memory for out-of-place FFT
		unsigned long nv = shtns->nspat;
QX		BrF = (double*) VMALLOC( nv * sizeof(double) );
VX		BtF = (double*) VMALLOC( 2*nv * sizeof(double) );
VX		BpF = BtF + nv;
3		BrF = (double*) VMALLOC( 3*nv * sizeof(double) );
3		BtF = BrF + nv;		BpF = BtF + nv;
	}

	#pragma omp parallel num_threads(shtns->nthreads)
	{
		const int nblk = (NLAT/2) / shtns->nthreads;
		if (shtns->fftc_mode != 1) {
Q			#pragma omp for schedule(dynamic) nowait
Q			for (int k=0; k<shtns->nthreads; k++)
Q				fftw_execute_dft(shtns->fftc_block, ((cplx *) Vr) + k*nblk, ((cplx *) BrF) + k*nblk);
V			#pragma omp for schedule(dynamic) nowait
V			for (int k=0; k<shtns->nthreads; k++) 
V				fftw_execute_dft(shtns->fftc_block, ((cplx *) Vt) + k*nblk, ((cplx *) BtF) + k*nblk);
V			#pragma omp for schedule(dynamic) nowait
V			for (int k=0; k<shtns->nthreads; k++) 
V				fftw_execute_dft(shtns->fftc_block, ((cplx *) Vp) + k*nblk, ((cplx *) BpF) + k*nblk);
		} else {
Q			#pragma omp for schedule(dynamic) nowait
Q			for (int k=0; k<shtns->nthreads; k++) 
Q				fftw_execute_split_dft(shtns->fftc_block, Vr+NPHI*(1+2*k*nblk), Vr+NPHI*2*k*nblk, BrF+1+NPHI*2*k*nblk, BrF+NPHI*2*k*nblk);
V			#pragma omp for schedule(dynamic) nowait
V			for (int k=0; k<shtns->nthreads; k++) 
V				fftw_execute_split_dft(shtns->fftc_block, Vt+NPHI*(1+2*k*nblk), Vt+NPHI*2*k*nblk, BtF+1+NPHI*2*k*nblk, BtF+NPHI*2*k*nblk);
V			#pragma omp for schedule(dynamic) nowait
V			for (int k=0; k<shtns->nthreads; k++) 
V				fftw_execute_split_dft(shtns->fftc_block, Vp+NPHI*(1+2*k*nblk), Vp+NPHI*2*k*nblk, BpF+1+NPHI*2*k*nblk, BpF+NPHI*2*k*nblk);
		}
		#pragma omp barrier
  #else
	#pragma omp parallel num_threads(shtns->nthreads)
	{
  #endif
		if (llim < SHT_L_RESCALE_FLY) {
			#pragma omp for schedule(dynamic,1) nowait
			for (int im=0; im <= imlim/2; im++) {
QX				GEN3(_an1,NWAY,SUFFIX)(shtns, BrF, Qlm, llim, im);
VX				GEN3(_an2,NWAY,SUFFIX)(shtns, BtF, BpF, Slm, Tlm, llim, im);
3				GEN3(_an3,NWAY,SUFFIX)(shtns, BrF, BtF, BpF, Qlm, Slm, Tlm, llim, im);
				if (imlim-im > im) {
QX					GEN3(_an1,NWAY,SUFFIX)(shtns, BrF, Qlm, llim, imlim-im);
VX					GEN3(_an2,NWAY,SUFFIX)(shtns, BtF, BpF, Slm, Tlm, llim, imlim-im);
3					GEN3(_an3,NWAY,SUFFIX)(shtns, BrF, BtF, BpF, Qlm, Slm, Tlm, llim, imlim-im);
				}
			}
		} else {
			#pragma omp for schedule(dynamic,1) nowait
			for (int im=0; im<=imlim; im++) {
QX				GEN3(_an1_hi,NWAY,SUFFIX)(shtns, BrF, Qlm, llim, im);
VX				GEN3(_an2_hi,NWAY,SUFFIX)(shtns, BtF, BpF, Slm, Tlm, llim, im);
3				GEN3(_an3_hi,NWAY,SUFFIX)(shtns, BrF, BtF, BpF, Qlm, Slm, Tlm, llim, im);
				if (imlim-im > im) {
QX					GEN3(_an1_hi,NWAY,SUFFIX)(shtns, BrF, Qlm, llim, imlim-im);
VX					GEN3(_an2_hi,NWAY,SUFFIX)(shtns, BtF, BpF, Slm, Tlm, llim, imlim-im);
3					GEN3(_an3_hi,NWAY,SUFFIX)(shtns, BrF, BtF, BpF, Qlm, Slm, Tlm, llim, imlim-im);
				}
			}
		}
		#ifndef SHT_AXISYM
			if (imlim < MMAX) {		// zero out m > imlim
				long l = LiM(shtns, (imlim+1)*MRES, imlim+1);
Q				#pragma omp single nowait
Q				memset(Qlm+l, 0, (shtns->nlm - l)*sizeof(cplx));
V				#pragma omp single nowait
V				memset(Slm+l, 0, (shtns->nlm - l)*sizeof(cplx));
V				#pragma omp single nowait
V				memset(Tlm+l, 0, (shtns->nlm - l)*sizeof(cplx));
			}
		#endif
	}

  #ifndef SHT_AXISYM
  	if (shtns->fftc_mode > 0) {		// free memory
Q	    VFREE(BrF);
VX	    VFREE(BtF);	// this frees also BpF.
	}
  #endif
  }

	#undef LSPAN
