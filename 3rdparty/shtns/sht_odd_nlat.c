#include "sht_private.h"

void SH_to_spat_odd_nlat(shtns_cfg shtns, cplx* Qlm, double* q, const long llim)
{
	const int nlat = shtns->nlat;
	const int nphi = shtns->nphi;
	const int ncplx = nphi/2 +1;

	int mlim = shtns->mmax;
	if (mlim*MRES > (unsigned) llim) mlim = ((unsigned) llim) /MRES;
	cplx* qm = (cplx*) VMALLOC(sizeof(cplx) * ncplx * nlat);
	for (int im=0; im <= mlim; im++) {
		long lm = LiM(shtns, im*MRES, im);
		SH_to_spat_ml(shtns, im, Qlm + lm, qm + im*nlat, llim);
	}
	if LIKELY(nphi > 1) {
		memset(qm + (mlim+1)*nlat, 0, (ncplx - (mlim+1))*nlat * sizeof(cplx));
		fftw_execute_dft_c2r(shtns->ifftc, qm, q);
	} else {
		for (int i=0; i<nlat; i++) q[i] = creal( qm[i] );	// compress to real
	}
	VFREE(qm);
}


void SHsphtor_to_spat_odd_nlat(shtns_cfg shtns, cplx* Slm, cplx* Tlm, double* vt, double* vp, const long llim)
{
	const int nlat = shtns->nlat;
	const int nphi = shtns->nphi;
	const int ncplx = nphi/2 +1;
	
	int mlim = shtns->mmax;
	if (mlim*MRES > (unsigned) llim) mlim = ((unsigned) llim) /MRES;
	cplx* vtm = (cplx*) VMALLOC(sizeof(cplx) * ncplx * nlat);
	cplx* vpm = (cplx*) VMALLOC(sizeof(cplx) * ncplx * nlat);
	for (int im=0; im <= mlim; im++) {
		long lm = LiM(shtns, im*MRES, im);
		SHsphtor_to_spat_ml(shtns, im, Slm + lm, Tlm + lm, vtm + im*nlat, vpm + im*nlat, llim);
	}
	if LIKELY(nphi > 1) {
		memset(vtm + (mlim+1)*nlat, 0, (ncplx - (mlim+1))*nlat * sizeof(cplx));
		fftw_execute_dft_c2r(shtns->ifftc, vtm, vt);
		memset(vpm + (mlim+1)*nlat, 0, (ncplx - (mlim+1))*nlat * sizeof(cplx));
		fftw_execute_dft_c2r(shtns->ifftc, vpm, vp);
	} else {
		for (int i=0; i<nlat; i++) vt[i] = creal( vtm[i] );	// compress to real
		for (int i=0; i<nlat; i++) vp[i] = creal( vpm[i] );	// compress to real
	}
	VFREE(vpm);
	VFREE(vtm);
}

void SHqst_to_spat_odd_nlat(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double *Vr, double *Vt, double *Vp, int ltr)
{
	SH_to_spat_odd_nlat(shtns, Qlm, Vr, ltr);
	SHsphtor_to_spat_odd_nlat(shtns, Slm, Tlm, Vt, Vp, ltr);
}

void SHsph_to_spat_odd_nlat(shtns_cfg shtns, cplx* Slm, double* vt, double* vp, const long llim)
{
	const int nlat = shtns->nlat;
	const int nphi = shtns->nphi;
	const int ncplx = nphi/2 +1;
	
	int mlim = shtns->mmax;
	if (mlim*MRES > (unsigned) llim) mlim = ((unsigned) llim) /MRES;
	cplx* vtm = (cplx*) VMALLOC(sizeof(cplx) * ncplx * nlat);
	cplx* vpm = (cplx*) VMALLOC(sizeof(cplx) * ncplx * nlat);
	for (int im=0; im <= mlim; im++) {
		long lm = LiM(shtns, im*MRES, im);
		SHsph_to_spat_ml(shtns, im, Slm + lm, vtm + im*nlat, vpm + im*nlat, llim);
	}
	if LIKELY(nphi > 1) {
		memset(vtm + (mlim+1)*nlat, 0, (ncplx - (mlim+1))*nlat * sizeof(cplx));
		fftw_execute_dft_c2r(shtns->ifftc, vtm, vt);
		memset(vpm + (mlim+1)*nlat, 0, (ncplx - (mlim+1))*nlat * sizeof(cplx));
		fftw_execute_dft_c2r(shtns->ifftc, vpm, vp);
	} else {
		for (int i=0; i<nlat; i++) vt[i] = creal( vtm[i] );	// compress to real
		if (vp) for (int i=0; i<nlat; i++) vp[i] = 0.0;	// v_phi = 0 for m=0
	}
	VFREE(vpm);
	VFREE(vtm);
}

void SHtor_to_spat_odd_nlat(shtns_cfg shtns, cplx* Tlm, double* vt, double* vp, const long llim)
{
	const int nlat = shtns->nlat;
	const int nphi = shtns->nphi;
	const int ncplx = nphi/2 +1;
	
	int mlim = shtns->mmax;
	if (mlim*MRES > (unsigned) llim) mlim = ((unsigned) llim) /MRES;
	cplx* vtm = (cplx*) VMALLOC(sizeof(cplx) * ncplx * nlat);
	cplx* vpm = (cplx*) VMALLOC(sizeof(cplx) * ncplx * nlat);
	for (int im=0; im <= mlim; im++) {
		long lm = LiM(shtns, im*MRES, im);
		SHtor_to_spat_ml(shtns, im, Tlm + lm, vtm + im*nlat, vpm + im*nlat, llim);
	}
	if LIKELY(nphi > 1) {
		memset(vtm + (mlim+1)*nlat, 0, (ncplx - (mlim+1))*nlat * sizeof(cplx));
		fftw_execute_dft_c2r(shtns->ifftc, vtm, vt);
		memset(vpm + (mlim+1)*nlat, 0, (ncplx - (mlim+1))*nlat * sizeof(cplx));
		fftw_execute_dft_c2r(shtns->ifftc, vpm, vp);
	} else {
		for (int i=0; i<nlat; i++) vp[i] = creal( vpm[i] );	// compress to real
		if (vt) for (int i=0; i<nlat; i++) vt[i] = 0.0;		// v_theta = 0 for m=0
	}
	VFREE(vpm);
	VFREE(vtm);
}

void spat_to_SH_odd_nlat(shtns_cfg shtns, double* q, cplx* Qlm, const long llim)
{
	const int nlat = shtns->nlat;
	const int nphi = shtns->nphi;
	const int ncplx = nphi/2 +1;
	const double norm = 1.0 / nphi;

	cplx* qm = (cplx*) VMALLOC(sizeof(cplx) * ncplx * nlat);
	if LIKELY(nphi > 1) {
		fftw_execute_dft_r2c(shtns->fftc, q, qm);
	} else {
		for (int i=0; i<nlat; i++) qm[i] = q[i];	// real to complex
	}
	int mlim = shtns->mmax;
	if (mlim*MRES > llim) mlim = llim /MRES;
	for (int im=0; im <= mlim; im++) {
		long lm = LiM(shtns, im*MRES, im);
		spat_to_SH_ml(shtns, im, qm + im*nlat, Qlm + lm, llim);
		for (int k=0; k<=llim-im*MRES; k++) Qlm[lm+k] *= norm;
	}
	if (mlim < shtns->mmax) {
		long lm = LiM(shtns, (mlim+1)*MRES, mlim+1);
		memset(Qlm+lm, 0, (shtns->nlm - lm)*sizeof(cplx));
	}
	VFREE(qm);
}

void spat_to_SHsphtor_odd_nlat(shtns_cfg shtns, double* vt, double* vp, cplx* Slm, cplx* Tlm, const long llim)
{
	const int nlat = shtns->nlat;
	const int nphi = shtns->nphi;
	const int ncplx = nphi/2 +1;
	const double norm = 1.0 / nphi;

	cplx* vtm = (cplx*) VMALLOC(sizeof(cplx) * ncplx * nlat);
	cplx* vpm = (cplx*) VMALLOC(sizeof(cplx) * ncplx * nlat);
	if LIKELY(nphi > 1) {
		fftw_execute_dft_r2c(shtns->fftc, vt, vtm);
		fftw_execute_dft_r2c(shtns->fftc, vp, vpm);
	} else {
		for (int i=0; i<nlat; i++) vtm[i] = vt[i];	// real to complex
		for (int i=0; i<nlat; i++) vpm[i] = vp[i];	// real to complex
	}
	int mlim = shtns->mmax;
	if (mlim*MRES > llim) mlim = llim /MRES;
	for (int im=0; im <= mlim; im++) {
		long lm = LiM(shtns, im*MRES, im);
		spat_to_SHsphtor_ml(shtns, im, vtm + im*nlat, vpm + im*nlat, Slm + lm, Tlm +lm, llim);
		for (int k=0; k<=llim-im*MRES; k++) {	Slm[lm+k] *= norm;		Tlm[lm+k] *= norm;  }
	}
	if (mlim < shtns->mmax) {
		long lm = LiM(shtns, (mlim+1)*MRES, mlim+1);
		memset(Slm+lm, 0, (shtns->nlm - lm)*sizeof(cplx));
		memset(Tlm+lm, 0, (shtns->nlm - lm)*sizeof(cplx));
	}
	VFREE(vpm);
	VFREE(vtm);
}

void spat_to_SHqst_odd_nlat(shtns_cfg shtns, double *Vr, double *Vt, double *Vp, cplx *Qlm, cplx *Slm, cplx *Tlm, int ltr)
{
	spat_to_SH_odd_nlat(shtns, Vr, Qlm, ltr);
	spat_to_SHsphtor_odd_nlat(shtns, Vt,Vp, Slm,Tlm, ltr);
}

void* fodd[SHT_NTYP] =
	{ SH_to_spat_odd_nlat, spat_to_SH_odd_nlat, SHsphtor_to_spat_odd_nlat, spat_to_SHsphtor_odd_nlat,
		SHsph_to_spat_odd_nlat, SHtor_to_spat_odd_nlat, SHqst_to_spat_odd_nlat, spat_to_SHqst_odd_nlat };

