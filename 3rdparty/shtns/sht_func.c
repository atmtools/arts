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

/** \internal \file sht_func.c
 * \brief Rotation of Spherical Harmonics.
 */


/** \addtogroup gimbutas_rotation Pseudo-spectral rotations of Spherical Harmonic fields (deprecated).
 \b deprecated use \ref rotation instead.
Rotation around axis other than Z should be considered of beta quality (they have been tested but may still contain bugs).
They also require \c mmax = \c lmax. They use an Algorithm inspired by the pseudospectral rotation described in
Gimbutas Z. and Greengard L. 2009 "A fast and stable method for rotating spherical harmonic expansions" <i>Journal of Computational Physics</i>.
doi:<a href="http://dx.doi.org/10.1016/j.jcp.2009.05.014">10.1016/j.jcp.2009.05.014</a>

These functions do only require a call to \ref shtns_create, but not to \ref shtns_set_grid.
*/
///@{

/// Rotate a SH representation Qlm around the z-axis by angle alpha (in radians),
/// which is the same as rotating the reference frame by angle -alpha.
/// Result is stored in Rlm (which can be the same array as Qlm).
void SH_Zrotate(shtns_cfg shtns, cplx *Qlm, double alpha, cplx *Rlm)
{
	int im, l, lmax, mmax, mres;

	lmax = shtns->lmax;		mmax = shtns->mmax;		mres = shtns->mres;

	if (Rlm != Qlm) {		// copy m=0 which does not change.
		l=0;	do { Rlm[l] = Qlm[l]; } while(++l <= lmax);
	}
	for (int im=1; im<=mmax; im++) {
		cplx eima = cos(im*mres*alpha) - I*sin(im*mres*alpha);		// rotate reference frame by angle -alpha
		for (l=im*mres; l<=lmax; ++l)	Rlm[LiM(shtns, l, im)] = Qlm[LiM(shtns, l, im)] * eima;
	}
}

///@}

/// \internal initialize pseudo-spectral rotations
static void SH_rotK90_init(shtns_cfg shtns)
{
	cplx *q;
	double *q0;
	int nfft, nrembed, ncembed;
	
//	if ((shtns->mres != 1) || (shtns->mmax != shtns->lmax)) runerr("Arbitrary rotations require lmax=mmax and mres=1");

#define NWAY 4

	const int lmax = shtns->lmax;
	const int fac = 2*VSIZE2*NWAY;		// we need a multiple of 2*VSIZE2*NWAY ...
	const int ntheta = fft_int( ((lmax+fac)/fac) , 7) * fac;		// ... and also an fft-friendly value

	// generate the equispaced grid for synthesis
	shtns->ct_rot = malloc( sizeof(double)*ntheta );
	shtns->st_rot = shtns->ct_rot + (ntheta/2);
	for (int k=0; k<ntheta/2; ++k) {
		double cost = cos(((0.5*M_PI)*(2*k+1))/ntheta);
		double sint = sqrt((1.0-cost)*(1.0+cost));
		shtns->ct_rot[k] = cost;
		shtns->st_rot[k] = sint;
	}

	// plan FFT
	size_t sze = sizeof(double)*(2*ntheta+2)*lmax;
	q0 = VMALLOC(sze);		// alloc.
	#ifdef OMP_FFTW
		int k = (lmax < 63) ? 1 : shtns->nthreads;
		fftw_plan_with_nthreads(k);
	#endif
	q = (cplx*) q0;		// in-place FFT
	nfft = 2*ntheta;	nrembed = nfft+2;		ncembed = nrembed/2;
	shtns->fft_rot = fftw_plan_many_dft_r2c(1, &nfft, lmax, q0, &nrembed, lmax, 1, q, &ncembed, lmax, 1, FFTW_MEASURE);

	VFREE(q0);
	shtns->npts_rot = ntheta;		// save ntheta, and mark as initialized.
}

/** \internal rotation kernel used by SH_Yrotate90(), SH_Xrotate90() and SH_rotate().
 Algorithm based on the pseudospectral rotation[1] :
 - rotate around Z by angle dphi0.
 - synthetize for each l the spatial description for phi=0 and phi=pi on an equispaced latitudinal grid.
 - Fourier ananlyze as data on the equator to recover the m in the 90 degrees rotated frame.
 - rotate around new Z by angle dphi1.
 [1] Gimbutas Z. and Greengard L. 2009 "A fast and stable method for rotating spherical harmonic expansions" Journal of Computational Physics. **/
static void SH_rotK90(shtns_cfg shtns, cplx *Qlm, cplx *Rlm, double dphi0, double dphi1)
{
//	if (shtns->npts_rot == 0)	SH_rotK90_init(shtns);

//	ticks tik0, tik1, tik2, tik3;

	const int lmax = shtns->lmax;
	const int ntheta = shtns->npts_rot;
	size_t sze = sizeof(double)*(2*ntheta+2)*lmax;
	double* const q0 = VMALLOC(sze);		// alloc.

	// rotate around Z by dphi0,  and also pre-multiply imaginary parts by m
	if (Rlm != Qlm) {		// copy m=0 which does not change.
		long l=0;	do { Rlm[l] = Qlm[l]; } while(++l <= lmax);
	}
	for (int m=1; m<=lmax; m++) {
		cplx eima = cos(m*dphi0) - I*sin(m*dphi0);		// rotate reference frame by angle -dphi0
		long lm = LiM(shtns,m,m);
		double em = m;
		for (long l=m; l<=lmax; ++l) {
			cplx qrot = Qlm[lm] * eima;
			((double*)Rlm)[2*lm]   = creal(qrot);
			((double*)Rlm)[2*lm+1] = cimag(qrot) * em;		// multiply every imaginary part by m  (part of im/sin(theta)*Ylm)
			lm++;
		}
	}
	Qlm = Rlm;

//	tik0 = getticks();

		rnd* const qve = (rnd*) VMALLOC( sizeof(rnd)*NWAY*4*lmax );	// vector buffer
		rnd* const qvo = qve + NWAY*2*lmax;		// for odd m
		double* const ct = shtns->ct_rot;
		double* const st = shtns->st_rot;
		double* const alm = shtns->alm;
		const long nk = ntheta/(2*VSIZE2);		// ntheta is a multiple of (2*VSIZE2)
		long k = 0;
		do {
			rnd cost[NWAY], y0[NWAY], y1[NWAY];
			long l=0;
			// m=0
			double*	al = alm;
			for (int j=0; j<NWAY; ++j) {
				y0[j] = vall(al[0]) / vread(st, j+k);		// l=0  (discarded) DIVIDED by sin(theta) [to be consistent with m>0]
				cost[j] = vread(ct, j+k);
			}
			for (int j=0; j<NWAY; ++j) {
				y1[j]  = vall(al[1]) * y0[j] * cost[j];
			}
			al += 2;	l+=2;
			while(l<=lmax) {
				for (int j=0; j<NWAY; ++j) {
					qve[ (l-2)*2*NWAY + 2*j]   = y1[j] * vall(creal(Qlm[l-1]));	// l-1
					qve[ (l-2)*2*NWAY + 2*j+1] = vall(0.0);
					qvo[ (l-2)*2*NWAY + 2*j]   = vall(0.0);
					qvo[ (l-2)*2*NWAY + 2*j+1] = vall(0.0);
					y0[j]  = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; ++j) {
					qve[ (l-1)*2*NWAY + 2*j]   = y0[j] * vall(creal(Qlm[l]));	// l
					qve[ (l-1)*2*NWAY + 2*j+1] = vall(0.0);
					qvo[ (l-1)*2*NWAY + 2*j]   = vall(0.0);
					qvo[ (l-1)*2*NWAY + 2*j+1] = vall(0.0);
					y1[j]  = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
				}
				al+=4;	l+=2;
			}
			if (l==lmax+1) {
				for (int j=0; j<NWAY; ++j) {
					qve[ (l-2)*2*NWAY + 2*j]   = y1[j] * vall(creal(Qlm[l-1]));	// l-1
					qve[ (l-2)*2*NWAY + 2*j+1] = vall(0.0);
					qvo[ (l-2)*2*NWAY + 2*j]   = vall(0.0);
					qvo[ (l-2)*2*NWAY + 2*j+1] = vall(0.0);
				}
			}
			// m > 0
			for (long m=1; m<=lmax; ++m) {
				rnd* qv = qve;
				if (m&1) qv = qvo;		// store even and odd m separately.
				double*	al = shtns->alm + m*(2*(lmax+1) -m+1);
				cplx* Ql = &Qlm[LiM(shtns, 0,m)];	// virtual pointer for l=0 and m
				rnd cost[NWAY], y0[NWAY], y1[NWAY];
				for (int j=0; j<NWAY; ++j) {
					cost[j] = vread(st, k+j);
					y0[j] = vall(2.0);		// *2 for m>0
				}
				long l=m-1;
				long int ny = 0;
				  if ((int)lmax <= SHT_L_RESCALE_FLY) {
					do {		// sin(theta)^(m-1)
						if (l&1) for (int j=0; j<NWAY; ++j) y0[j] *= cost[j];
						for (int j=0; j<NWAY; ++j) cost[j] *= cost[j];
					} while(l >>= 1);
				  } else {
					long int nsint = 0;
					do {		// sin(theta)^(m-1)		(use rescaling to avoid underflow)
						if (l&1) {
							for (int j=0; j<NWAY; ++j) y0[j] *= cost[j];
							ny += nsint;
							if (vlo(y0[0]) < (SHT_ACCURACY+1.0/SHT_SCALE_FACTOR)) {
								ny--;
								for (int j=0; j<NWAY; ++j) y0[j] *= vall(SHT_SCALE_FACTOR);
							}
						}
						for (int j=0; j<NWAY; ++j) cost[j] *= cost[j];
						nsint += nsint;
						if (vlo(cost[0]) < 1.0/SHT_SCALE_FACTOR) {
							nsint--;
							for (int j=0; j<NWAY; ++j) cost[j] *= vall(SHT_SCALE_FACTOR);
						}
					} while(l >>= 1);
				  }
				for (int j=0; j<NWAY; ++j) {
					y0[j] *= vall(al[0]);
					cost[j] = vread(ct, j+k);
				}
				for (int j=0; j<NWAY; ++j) {
					y1[j]  = (vall(al[1])*y0[j]) *cost[j];
				}
				l=m;		al+=2;
				while ((ny<0) && (l<lmax)) {		// ylm treated as zero and ignored if ny < 0
					for (int j=0; j<NWAY; ++j) {
						y0[j] = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
					}
					for (int j=0; j<NWAY; ++j) {
						y1[j] = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
					}
					l+=2;	al+=4;
					if (fabs(vlo(y0[NWAY-1])) > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0) {		// rescale when value is significant
						++ny;
						for (int j=0; j<NWAY; ++j) {
							y0[j] *= vall(1.0/SHT_SCALE_FACTOR);		y1[j] *= vall(1.0/SHT_SCALE_FACTOR);
						}
					}
				}
			  if (ny == 0) {
				while (l<lmax) {
					rnd qr = vall(creal(Ql[l]));		rnd qi = vall(cimag(Ql[l]));
					for (int j=0; j<NWAY; ++j) {
						qv[ (l-1)*2*NWAY + 2*j]   += y0[j] * qr;	// l
						qv[ (l-1)*2*NWAY + 2*j+1] += y0[j] * qi;
					}
					qr = vall(creal(Ql[l+1]));		qi = vall(cimag(Ql[l+1]));
					for (int j=0; j<NWAY; ++j) {
						qv[ (l)*2*NWAY + 2*j]   += y1[j] * qr;	// l+1
						qv[ (l)*2*NWAY + 2*j+1] += y1[j] * qi;
					}
					for (int j=0; j<NWAY; ++j) {
						y0[j] = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];	// l+2
					}
					for (int j=0; j<NWAY; ++j) {
						y1[j] = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];	// l+3
					}
					l+=2;	al+=4;
				}
				if (l==lmax) {
					rnd qr = vall(creal(Ql[l]));		rnd qi = vall(cimag(Ql[l]));
					for (int j=0; j<NWAY; ++j) {
						qv[ (l-1)*2*NWAY + 2*j]   += y0[j] * qr;	// l
						qv[ (l-1)*2*NWAY + 2*j+1] += y0[j] * qi;
					}
				}
			  }
			}
			// construct ring using symmetries + transpose...
			double signl = -1.0;
			double* qse = (double*) qve;
			double* qso = (double*) qvo;
			for (long l=1; l<=lmax; ++l) {
				for (int j=0; j<NWAY; j++) {
					for (int i=0; i<VSIZE2; i++) {
						double qre = qse[(l-1)*2*NWAY*VSIZE2 + 2*j*VSIZE2 + i];		// m even
						double qie = qse[(l-1)*2*NWAY*VSIZE2 + (2*j+1)*VSIZE2 + i];
						double qro = qso[(l-1)*2*NWAY*VSIZE2 + 2*j*VSIZE2 + i];		// m odd
						double qio = qso[(l-1)*2*NWAY*VSIZE2 + (2*j+1)*VSIZE2 + i];
						long ijk = (k+j)*VSIZE2 + i;
						qre *= st[ijk];			qro *= st[ijk];		// multiply real part by sin(theta)  [to get Ylm from Ylm/sin(theta)]
						// because qr and qi map on different parities with respect to the future Fourier tranform, we can add them !!
						// note that this may result in leak between even and odd m's if their amplitude is widely different.
						q0[ijk*lmax +(l-1)]              =  (qre + qro) - (qie + qio);
						q0[(ntheta+ijk)*lmax +(l-1)]     = ((qre + qro) + (qie + qio)) * signl;				// * (-1)^l
						q0[(2*ntheta-1-ijk)*lmax +(l-1)] =  (qre - qro) + (qie - qio);
						q0[(ntheta-1-ijk)*lmax +(l-1)]   = ((qre - qro) - (qie - qio)) * signl;				// (-1)^(l-m)
					}
				}
				signl *= -1.0;
			}
			k += NWAY;
		} while (k<nk);
	VFREE(qve);
#undef NWAY

//	tik1 = getticks();

	// perform FFT
	cplx* q = (cplx*) q0;		// in-place FFT
	fftw_execute_dft_r2c(shtns->fft_rot, q0, q);

//	tik2 = getticks();

	const int nphi = 2*ntheta;
	double ydyl[lmax+1];
	long m=0;		long lm=1;		// start at l=1,m=0
	long l;
		legendre_sphPlm_deriv_array_equ(shtns, lmax, m, ydyl+m);
		for (l=1; l<lmax; l+=2) {
			Rlm[lm] =  -creal(q[m*lmax +(l-1)])/(ydyl[l]*nphi);
			Rlm[lm+1] =  creal(q[m*lmax +l])/(ydyl[l+1]*nphi);
			lm+=2;
		}
		if (l==lmax) {
			Rlm[lm] =  -creal(q[m*lmax +(l-1)])/(ydyl[l]*nphi);
			lm+=1;
		}
	dphi1 += M_PI/nphi;	// shift rotation angle by angle of first synthesis latitude.
	for (m=1; m<=lmax; ++m) {
		legendre_sphPlm_deriv_array_equ(shtns, lmax, m, ydyl+m);
		cplx eimdp = (cos(m*dphi1) - I*sin(m*dphi1))/nphi;
		for (l=m; l<lmax; l+=2) {
			Rlm[lm] =  eimdp*q[m*lmax +(l-1)]*(1./ydyl[l]);
			Rlm[lm+1] =  eimdp*q[m*lmax +l]*(-1./ydyl[l+1]);
			lm+=2;
		}
		if (l==lmax) {
			Rlm[lm] =  eimdp*q[m*lmax +(l-1)]*(1./ydyl[l]);
			lm++;
		}
	}
	VFREE(q0);

//	tik3 = getticks();
//	printf("    tick ratio : %.3f  %.3f  %.3f\n", elapsed(tik1,tik0)/elapsed(tik3,tik0), elapsed(tik2,tik1)/elapsed(tik3,tik0), elapsed(tik3,tik2)/elapsed(tik3,tik0));

}


/// \addtogroup gimbutas_rotation
///@{

/// rotate Qlm by 90 degrees around X axis and store the result in Rlm.
/// shtns->mres MUST be 1, and lmax=mmax.
void SH_Xrotate90(shtns_cfg shtns, cplx *Qlm, cplx *Rlm)
{
	int lmax= shtns->lmax;
	if ((shtns->mres != 1) || (shtns->mmax < lmax)) shtns_runerr("truncature makes rotation not closed.");

	if (lmax == 1) {
		Rlm[0] = Qlm[0];	// l=0 is invariant.
		int l=1;													// rotation matrix for rotX(90), l=1 : m=[0, 1r, 1i]
			double q0 = creal(Qlm[LiM(shtns, l, 0)]);
			Rlm[LiM(shtns, l, 0)] = sqrt(2.0) * cimag(Qlm[LiM(shtns, l, 1)]);			//[m=0]     0        0    sqrt(2)
			Rlm[LiM(shtns, l ,1)] = creal(Qlm[LiM(shtns, l, 1)]) - I*(sqrt(0.5)*q0);	//[m=1r]    0        1      0
		return;																			//[m=1i] -sqrt(2)/2  0      0
	}

	SH_rotK90(shtns, Qlm, Rlm, 0.0,  -M_PI/2);
}

/// rotate Qlm by 90 degrees around Y axis and store the result in Rlm.
/// shtns->mres MUST be 1, and lmax=mmax.
void SH_Yrotate90(shtns_cfg shtns, cplx *Qlm, cplx *Rlm)
{
	int lmax= shtns->lmax;
	if ((shtns->mres != 1) || (shtns->mmax < lmax)) shtns_runerr("truncature makes rotation not closed.");

	if (lmax == 1) {
		Rlm[0] = Qlm[0];	// l=0 is invariant.
		int l=1;											// rotation matrix for rotY(90), l=1 : m=[0, 1r, 1i]
			double q0 = creal(Qlm[LiM(shtns, l, 0)]);									//[m=0]       0     sqrt(2)  0
			Rlm[LiM(shtns, l, 0)] = sqrt(2.0) * creal(Qlm[LiM(shtns, l, 1)]);			//[m=1r] -sqrt(2)/2   0      0
			Rlm[LiM(shtns, l ,1)] = I*cimag(Qlm[LiM(shtns, l, 1)]) - sqrt(0.5) * q0;	//[m=1i]      0       0      1
		return;
	}

	SH_rotK90(shtns, Qlm, Rlm, -M_PI/2, 0.0);
}

/// rotate Qlm around Y axis by arbitrary angle, using composition of rotations. Store the result in Rlm.
void SH_Yrotate(shtns_cfg shtns, cplx *Qlm, double alpha, cplx *Rlm)
{
	if ((shtns->mres != 1) || (shtns->mmax < shtns->lmax)) shtns_runerr("truncature makes rotation not closed.");

	SH_rotK90(shtns, Qlm, Rlm, 0.0, M_PI/2 + alpha);	// Zrotate(pi/2) + Yrotate90 + Zrotate(pi+alpha)
	SH_rotK90(shtns, Rlm, Rlm, 0.0, M_PI/2);			// Yrotate90 + Zrotate(pi/2)
}

///@}



/** \addtogroup operators Special operators
 * Apply special operators in spectral space: multiplication by cos(theta), sin(theta).d/dtheta.
*/
///@{


/// \internal generates the cos(theta) matrix up to lmax+1
/// \param mx : an array of 2*NLM double that will be filled with the matrix coefficients.
/// xq[lm] = mx[2*lm-1] * q[lm-1] + mx[2*lm] * q[lm+1];			[note the shift in indices compared to the public functions]
static void mul_ct_matrix_shifted(shtns_cfg shtns, double* mx)
{
	long int im,l,lm;
	double a_1;

	if (SHT_NORM == sht_schmidt) {
		lm=0;
		for (im=0; im<=MMAX; im++) {
			double* al = alm_im(shtns,im);
			long int m=im*MRES;
			a_1 = 1.0 / al[1];
			l=m;
			while(++l <= LMAX) {
				al+=2;				
				mx[2*lm+1] = a_1;
				a_1 = 1.0 / al[1];
				mx[2*lm] = -a_1*al[0];        // = -al[2*(lm+1)] / al[2*(lm+1)+1];
				lm++;
			}
			if (l == LMAX+1) {	// the last one needs to be computed (used in vector to scalar transform)
				mx[2*lm+1] = a_1;
				mx[2*lm] = sqrt((l+m)*(l-m))/(2*l+1);		// LMAX+1
				lm++;
			}
		}
	} else {
		lm=0;
		for (im=0; im<=MMAX; im++) {
			double* al = alm_im(shtns, im);
			l=im*MRES;
			while(++l <= LMAX+1) {	// compute coeff up to LMAX+1, it fits into the 2*NLM bloc, and is needed for vector <> scalar conversions.
				a_1 = 1.0 / al[1];
				mx[2*lm] = a_1;		// specific to orthonormal.
				mx[2*lm+1] = a_1;
				lm++;	al+=2;
			}
		}
	}
}

static void st_dt_matrix_shifted(shtns_cfg shtns, double* mx)
{
	mul_ct_matrix_shifted(shtns, mx);
	for (int lm=0; lm<NLM; lm++) {
		mx[2*lm]   *= -(shtns->li[lm] + 2);		// coeff (l+1)
		mx[2*lm+1] *=   shtns->li[lm];			// coeff (l-1)
	}
}


/// fill mx with the coefficients for multiplication by cos(theta)
/// \param mx : an array of 2*NLM double that will be filled with the matrix coefficients.
/// xq[lm] = mx[2*lm] * q[lm-1] + mx[2*lm+1] * q[lm+1];
void mul_ct_matrix(shtns_cfg shtns, double* mx)
{
	long int im,l,lm;
	double a_1;

	mul_ct_matrix_shifted(shtns, mx);
	for (int lm=2*NLM-1; lm>0; lm--)	mx[lm] = mx[lm-1];		// shift back indices (copy, slow)
	mx[0] = 0.0;
	for (int im=1; im<=MMAX; im++) {				// remove the coeff for lmax+1 (for backward compatibility)
		int lm = LiM(shtns, im*MRES, im);
		mx[2*lm-1] = 0.0;		mx[2*lm] = 0.0;
	}
	mx[2*NLM-1] = 0.0;
}

/// fill mx with the coefficients of operator sin(theta).d/dtheta
/// \param mx : an array of 2*NLM double that will be filled with the matrix coefficients.
/// stdq[lm] = mx[2*lm] * q[lm-1] + mx[2*lm+1] * q[lm+1];
void st_dt_matrix(shtns_cfg shtns, double* mx)
{
	mul_ct_matrix(shtns, mx);
	for (int lm=0; lm<NLM; lm++) {
		mx[2*lm]   *=   shtns->li[lm] - 1;		// coeff (l-1)
		mx[2*lm+1] *= -(shtns->li[lm] + 2);		// coeff (l+1)
	}
}

/// Multiplication of Qlm by a matrix involving l+1 and l-1 only.
/// The result is stored in Rlm, which MUST be different from Qlm.
/// mx is an array of 2*NLM values as returned by \ref mul_ct_matrix or \ref st_dt_matrix
/// compute: Rlm[lm] = mx[2*lm] * Qlm[lm-1] + mx[2*lm+1] * Qlm[lm+1];
void SH_mul_mx(shtns_cfg shtns, double* mx, cplx *Qlm, cplx *Rlm)
{
	long int nlmlim, lm;
	v2d* vq = (v2d*) Qlm;
	v2d* vr = (v2d*) Rlm;
	nlmlim = NLM-1;
	lm = 0;
		s2d mxu = vdup(mx[1]);
		vr[0] = mxu*vq[1];
	for (lm=1; lm<nlmlim; lm++) {
		s2d mxl = vdup(mx[2*lm]);		s2d mxu = vdup(mx[2*lm+1]);
		vr[lm] = mxl*vq[lm-1] + mxu*vq[lm+1];
	}
	lm = nlmlim;
		s2d mxl = vdup(mx[2*lm]);
		vr[lm] = mxl*vq[lm-1];
}

///@}

// truncation at LMAX and MMAX
#define LTR LMAX
#define MTR MMAX

/** \addtogroup local Local and partial evaluation of SH fields.
 * These do only require a call to \ref shtns_create, but not to \ref shtns_set_grid.
 * These functions are not optimized and can be relatively slow, but they provide good
 * reference implemenation for the transforms.
*/
///@{

/// Evaluate scalar SH representation of complex field \b alm at physical point defined by \b cost = cos(theta) and \b phi
cplx SH_to_point_cplx(shtns_cfg shtns, cplx *alm, double cost, double phi)
{
	double yl[LMAX+1];
	long int l,m;
	cplx z = 0.0;

	v2d vr0 = vdup(0.0);		v2d vr1 = vdup(0.0);
	m=0;
		legendre_sphPlm_array(shtns, LTR, m, cost, &yl[m]);

		int ll = 0;
		for (l=0; l<LTR; l+=2) {
			ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
			vr0 += ((v2d*)alm)[ll] * vdup(yl[l]);
			ll += (l<MMAX) ? 2*l+2 : 2*MMAX+1;
			vr1 += ((v2d*)alm)[ll] * vdup(yl[l+1]);
		}
		if (l==LTR) {
			ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
			vr0 += ((v2d*)alm)[ll] * vdup(yl[l]);
		}
		vr0 += vr1;
		z = vcplx_real(vr0) + I*vcplx_imag(vr0);
	if (MTR>0) {
		cplx eip = cos(MRES*phi) + I*sin(MRES*phi);
		v2d vrc = vdup(0.0);
		v2d vrs = vdup(0.0);
		cplx eimp = eip;
		for (long im=1; im<=MTR; im++) {
			const long m = im*MRES;
			//memset(yl+m, 0xFF, sizeof(double)*(LMAX+1-m));
			long lnz = legendre_sphPlm_array(shtns, LTR, m, cost, &yl[m]);
			//if (lnz > m) printf("m=%d, lnz=%d  [ %g, %g, %g]\n", m, lnz, yl[lnz-1],yl[lnz],yl[lnz+1]);
			if (lnz > LTR) break;		// nothing else to do

			v2d vrm = vdup(0.0);		v2d vim = vdup(0.0);
			long ll = m*m;
			long l=m;
			cplx* almm = alm - 2*m;		// m<0

			//for (; l<lnz; l++) 	ll += (l<=MMAX) ? 2*l : 2*MMAX+1;	// skip zeros in yl, replaced by direct computation in next block:
			if (lnz > m) {	// skip zeros in yl:
				int lx = (lnz <= MMAX) ? lnz : MMAX+1;
				ll += ((m+lx-1)*(lx-m));					// for (l=m; l<min(lnz,MMAX+1); l++) ll += 2*l;
				if (lnz > MMAX+1) {
					ll += (2*MMAX+1)*(lnz-(MMAX+1));		// for (l=MMAX+1; l<lnz; l++) ll += 2*MMAX+1;
				}
				l=lnz;
			}
			for (; l<=LMAX; l++) {			// gather from cplx rep
				ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
				vim += vdup(yl[l]) * ((v2d*)almm)[ll];	// -m
				vrm += vdup(yl[l]) * ((v2d*)alm)[ll];	// +m
			}
			//cplx eimp = cos(m*phi) + I*sin(m*phi);		// we need something accurate here.
			if (m&1) vim = -vim;				// m<0, m odd
			//vrc += vdup(cos(m*phi)) * (vrm+vim);	// error @lmax=1023 = 2e-7
			//vrs += vdup(sin(m*phi)) * (vrm-vim);
			vrc += vdup(creal(eimp)) * (vrm+vim);	// error @lmax=1023 = 2e-9
			vrs += vdup(cimag(eimp)) * (vrm-vim);
			eimp *= eip;
		}
		z += vcplx_real(vrc) - vcplx_imag(vrs) + I*(vcplx_imag(vrc) + vcplx_real(vrs));
	}
	return z;
}


/// Evaluate scalar SH representation \b Qlm at physical point defined by \b cost = cos(theta) and \b phi
double SH_to_point(shtns_cfg shtns, cplx *Qlm, double cost, double phi)
{
	double yl[LMAX+1];
	double vr0, vr1;
	long int l,m,im;

	vr0 = 0.0;		vr1 = 0.0;
	m=0;	im=0;
		legendre_sphPlm_array(shtns, LTR, im, cost, &yl[m]);
		for (l=m; l<LTR; l+=2) {
			vr0 += yl[l]   * creal( Qlm[l] );
			vr1 += yl[l+1] * creal( Qlm[l+1] );
		}
		if (l==LTR) {
			vr0 += yl[l] * creal( Qlm[l] );
		}
		vr0 += vr1;
	
	for (im=1; im<=MTR; im++) {
		m = im*MRES;
		long lnz = legendre_sphPlm_array(shtns, LTR, im, cost, &yl[m]);
		if (lnz > LTR) break;		// nothing else to do

		v2d* Ql = (v2d*) &Qlm[LiM(shtns, 0,im)];	// virtual pointer for l=0 and im
		v2d vrm0 = vdup(0.0);		v2d vrm1 = vdup(0.0);
		for (l=lnz; l<LTR; l+=2) {
			vrm0 += vdup(yl[l])   * Ql[l];
			vrm1 += vdup(yl[l+1]) * Ql[l+1];
		}
		cplx eimp = 2.*(cos(m*phi) + I*sin(m*phi));		// we need something accurate here.
		vrm0 += vrm1;
		if (l==LTR) {
			vrm0 += vdup(yl[l]) * Ql[l];
		}
		vr0 += vcplx_real(vrm0)*creal(eimp) - vcplx_imag(vrm0)*cimag(eimp);
	}
	return vr0;
}

void SH_to_grad_point(shtns_cfg shtns, cplx *DrSlm, cplx *Slm, double cost, double phi,
					   double *gr, double *gt, double *gp)
{
	double yl[LMAX+1];
	double dtyl[LMAX+1];
	double vtt, vpp, vr0, vrm;
	long int l,m,im;

	const double sint = sqrt((1.-cost)*(1.+cost));
	vtt = 0.;  vpp = 0.;  vr0 = 0.;  vrm = 0.;
	m=0;	im=0;
		legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
		for (l=m; l<=LTR; ++l) {
			vr0 += yl[l] * creal( DrSlm[l] );
			vtt += dtyl[l] * creal( Slm[l] );
		}
	if (MTR>0) {
		im=1;  do {
			m = im*MRES;
			long lnz = legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
			if (lnz > LTR) break;		// nothing else to do

			cplx eimp = 2.*(cos(m*phi) + I*sin(m*phi));
			cplx imeimp = eimp*m*I;
			l = LiM(shtns, 0,im);
			v2d* Ql = (v2d*) &DrSlm[l];		v2d* Sl = (v2d*) &Slm[l];
			v2d qm = vdup(0.0);
			v2d dsdt = vdup(0.0);		v2d dsdp = vdup(0.0);
			for (l=lnz; l<=LTR; ++l) {
				qm += vdup(yl[l]) * Ql[l];
				dsdt += vdup(dtyl[l]) * Sl[l];
				dsdp += vdup(yl[l]) * Sl[l];
			}
			vrm += vcplx_real(qm)*creal(eimp) - vcplx_imag(qm)*cimag(eimp);			// dS/dr
			vtt += vcplx_real(dsdt)*creal(eimp) - vcplx_imag(dsdt)*cimag(eimp);		// dS/dt
			vpp += vcplx_real(dsdp)*creal(imeimp) - vcplx_imag(dsdp)*cimag(imeimp);	// + I.m/sint *S
		} while (++im <= MTR);
		vr0 += vrm*sint;
	}
	*gr = vr0;	// Gr = dS/dr
	*gt = vtt;	// Gt = dS/dt
	*gp = vpp;	// Gp = I.m/sint *S
}

/// Evaluate vector SH representation \b Qlm at physical point defined by \b cost = cos(theta) and \b phi
void SHqst_to_point(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double cost, double phi,
					   double *vr, double *vt, double *vp)
{
	double yl[LMAX+1];
	double dtyl[LMAX+1];
	double vtt, vpp, vr0, vrm;
	long int l,m,im;

	const double sint = sqrt((1.-cost)*(1.+cost));
	vtt = 0.;  vpp = 0.;  vr0 = 0.;  vrm = 0.;
	m=0;	im=0;
		legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
		for (l=m; l<=LTR; ++l) {
			vr0 += yl[l] * creal( Qlm[l] );
			vtt += dtyl[l] * creal( Slm[l] );
			vpp -= dtyl[l] * creal( Tlm[l] );
		}
	if (MTR>0) {
		im=1;  do {
			m = im*MRES;
			legendre_sphPlm_deriv_array(shtns, LTR, im, cost, sint, &yl[m], &dtyl[m]);
			cplx eimp = 2.*(cos(m*phi) + I*sin(m*phi));
			cplx imeimp = eimp*m*I;
			l = LiM(shtns, 0,im);
			v2d* Ql = (v2d*) &Qlm[l];	v2d* Sl = (v2d*) &Slm[l];	v2d* Tl = (v2d*) &Tlm[l];
			v2d qm = vdup(0.0);
			v2d dsdt = vdup(0.0);		v2d dtdt = vdup(0.0);
			v2d dsdp = vdup(0.0);		v2d dtdp = vdup(0.0);
			for (l=m; l<=LTR; ++l) {
				qm += vdup(yl[l]) * Ql[l];
				dsdt += vdup(dtyl[l]) * Sl[l];
				dtdt += vdup(dtyl[l]) * Tl[l];
				dsdp += vdup(yl[l]) * Sl[l];
				dtdp += vdup(yl[l]) * Tl[l];
			}
			vrm += vcplx_real(qm)*creal(eimp) - vcplx_imag(qm)*cimag(eimp);
			vtt += (vcplx_real(dtdp)*creal(imeimp) - vcplx_imag(dtdp)*cimag(imeimp))	// + I.m/sint *T
					+ (vcplx_real(dsdt)*creal(eimp) - vcplx_imag(dsdt)*cimag(eimp));	// + dS/dt
			vpp += (vcplx_real(dsdp)*creal(imeimp) - vcplx_imag(dsdp)*cimag(imeimp))	// + I.m/sint *S
					- (vcplx_real(dtdt)*creal(eimp) - vcplx_imag(dtdt)*cimag(eimp));	// - dT/dt
		} while (++im <= MTR);
		vr0 += vrm*sint;
	}
	*vr = vr0;
	*vt = vtt;	// Bt = I.m/sint *T  + dS/dt
	*vp = vpp;	// Bp = I.m/sint *S  - dT/dt
}
///@}
	
#undef LTR
#undef MTR


/*
	SYNTHESIS AT A GIVEN LATITUDE
	(does not require a previous call to shtns_set_grid)
*/

/// synthesis at a given latitude, on nphi equispaced longitude points.
/// vr, vt, and vp arrays must have nphi doubles allocated.
/// It does not require a previous call to shtns_set_grid, but it is NOT thread-safe, 
/// unless called with a different shtns_cfg
/// \ingroup local
void SHqst_to_lat(shtns_cfg shtns, cplx *Qlm, cplx *Slm, cplx *Tlm, double cost,
					double *vr, double *vt, double *vp, int nphi, int ltr, int mtr)
{
	cplx vst, vtt, vsp, vtp, vrr;
	cplx *vrc, *vtc, *vpc;
	double* ylm_lat;
	double* dylm_lat;
	double st_lat;

	if (ltr > LMAX) ltr=LMAX;
	if (mtr > MMAX) mtr=MMAX;
	if (mtr*MRES > ltr) mtr=ltr/MRES;
	if (mtr*2*MRES >= nphi) mtr = (nphi-1)/(2*MRES);

	ylm_lat = shtns->ylm_lat;
	if (ylm_lat == NULL) {		// alloc memory for Legendre functions ?
		ylm_lat = (double *) malloc(sizeof(double)* 2*NLM);
		shtns->ylm_lat = ylm_lat;
	}
	dylm_lat = ylm_lat + NLM;

	st_lat = sqrt((1.-cost)*(1.+cost));	// sin(theta)
	if (cost != shtns->ct_lat) {		// compute Legendre functions ?
		shtns->ct_lat = cost;
		for (int m=0,j=0; m<=mtr; ++m) {
			legendre_sphPlm_deriv_array(shtns, ltr, m, cost, st_lat, &ylm_lat[j], &dylm_lat[j]);
			j += LMAX -m*MRES +1;
		}
	}

	vrc = (cplx*) fftw_malloc(sizeof(double) * 3*(nphi+2));
	vtc = vrc + (nphi/2+1);
	vpc = vtc + (nphi/2+1);

	if (nphi != shtns->nphi_lat) {		// compute FFTW plan ?
		if (shtns->ifft_lat) fftw_destroy_plan(shtns->ifft_lat);
		#ifdef OMP_FFTW
			fftw_plan_with_nthreads(1);
		#endif
		shtns->ifft_lat = fftw_plan_dft_c2r_1d(nphi, vrc, vr, FFTW_ESTIMATE);
		shtns->nphi_lat = nphi;
	}

	for (int m = 0; m<nphi/2+1; ++m) {	// init with zeros
		vrc[m] = 0.0;	vtc[m] = 0.0;	vpc[m] = 0.0;
	}
	long j=0;
	int m=0;
		vrr=0;	vtt=0;	vst=0;
		for(int l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * creal(Qlm[j]);
			vst += dylm_lat[j] * creal(Slm[j]);
			vtt += dylm_lat[j] * creal(Tlm[j]);
		}
		j += (LMAX-ltr);
		vrc[m] = vrr;
		vtc[m] =  vst;	// Vt =   dS/dt
		vpc[m] = -vtt;	// Vp = - dT/dt
	for (int m=MRES; m<=mtr*MRES; m+=MRES) {
		vrr=0;	vtt=0;	vst=0;	vsp=0;	vtp=0;
		for(int l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * Qlm[j];
			vst += dylm_lat[j] * Slm[j];
			vtt += dylm_lat[j] * Tlm[j];
			vsp += ylm_lat[j] * Slm[j];
			vtp += ylm_lat[j] * Tlm[j];
		}
		j+=(LMAX-ltr);
		vrc[m] = vrr*st_lat;
		vtc[m] = I*m*vtp + vst;	// Vt = I.m/sint *T  + dS/dt
		vpc[m] = I*m*vsp - vtt;	// Vp = I.m/sint *S  - dT/dt
	}
	fftw_execute_dft_c2r(shtns->ifft_lat,vrc,vr);
	fftw_execute_dft_c2r(shtns->ifft_lat,vtc,vt);
	fftw_execute_dft_c2r(shtns->ifft_lat,vpc,vp);
	fftw_free(vrc);
}

/// synthesis at a given latitude, on nphi equispaced longitude points.
/// vr arrays must have nphi doubles allocated.
/// It does not require a previous call to shtns_set_grid, but it is NOT thread-safe,
/// unless called with a different shtns_cfg
/// \ingroup local
void SH_to_lat(shtns_cfg shtns, cplx *Qlm, double cost,
					double *vr, int nphi, int ltr, int mtr)
{
	cplx vrr;
	cplx *vrc;
	double* ylm_lat;
	double* dylm_lat;
	double st_lat;

	if (ltr > LMAX) ltr=LMAX;
	if (mtr > MMAX) mtr=MMAX;
	if (mtr*MRES > ltr) mtr=ltr/MRES;
	if (mtr*2*MRES >= nphi) mtr = (nphi-1)/(2*MRES);

	ylm_lat = shtns->ylm_lat;
	if (ylm_lat == NULL) {
		ylm_lat = (double *) malloc(sizeof(double)* 2*NLM);
		shtns->ylm_lat = ylm_lat;
	}
	dylm_lat = ylm_lat + NLM;

	st_lat = sqrt((1.-cost)*(1.+cost));	// sin(theta)
	if (cost != shtns->ct_lat) {
		shtns->ct_lat = cost;
		for (int m=0,j=0; m<=mtr; ++m) {
			legendre_sphPlm_deriv_array(shtns, ltr, m, cost, st_lat, &ylm_lat[j], &dylm_lat[j]);
			j += LMAX -m*MRES +1;
		}
	}

	vrc = (cplx*) fftw_malloc(sizeof(double) * (nphi+2));

	if (nphi != shtns->nphi_lat) {
		if (shtns->ifft_lat) fftw_destroy_plan(shtns->ifft_lat);
		#ifdef OMP_FFTW
			fftw_plan_with_nthreads(1);
		#endif
		shtns->ifft_lat = fftw_plan_dft_c2r_1d(nphi, vrc, vr, FFTW_ESTIMATE);
		shtns->nphi_lat = nphi;
	}

	for (int m = 0; m<nphi/2+1; ++m) {	// init with zeros
		vrc[m] = 0.0;
	}
	long j=0;
	int m=0;
		vrr=0;
		for(int l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * creal(Qlm[j]);
		}
		j += (LMAX-ltr);
		vrc[m] = vrr;
	for (int m=MRES; m<=mtr*MRES; m+=MRES) {
		vrr=0;
		for(int l=m; l<=ltr; ++l, ++j) {
			vrr += ylm_lat[j] * Qlm[j];
		}
		j+=(LMAX-ltr);
		vrc[m] = vrr*st_lat;
	}
	fftw_execute_dft_c2r(shtns->ifft_lat,vrc,vr);
	fftw_free(vrc);
}




// SPAT_CPLX transform indexing scheme:
// if (l<=MMAX) : l*(l+1) + m
// if (l>=MMAX) : l*(2*mmax+1) - mmax*mmax + m  = mmax*(2*l-mmax) + l+m
///\internal
void SH_2real_to_cplx(shtns_cfg shtns, cplx* Rlm, cplx* Ilm, cplx* Zlm)
{
	// combine into complex coefficients
	unsigned ll = 0;
	unsigned lm = 0;
	for (unsigned l=0; l<=LMAX; l++) {
		ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
		Zlm[ll] = creal(Rlm[lm]) + I*creal(Ilm[lm]);		// m=0
		lm++;
	}
	for (unsigned m=1; m<=MMAX; m++) {
		ll = (m-1)*m;
		for (unsigned l=m; l<=LMAX; l++) {
			ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
			cplx rr = Rlm[lm];
			cplx ii = Ilm[lm];
			Zlm[ll+m] = rr + I*ii;			// m>0
			rr = conj(rr) + I*conj(ii);		// m<0, m even
			if (m&1) rr = -rr;				// m<0, m odd
			Zlm[ll-m] = rr;
			lm++;
		}
	}
}

///\internal
void SH_cplx_to_2real(shtns_cfg shtns, cplx* Zlm, cplx* Rlm, cplx* Ilm)
{
	// extract complex coefficients corresponding to real and imag
	unsigned ll = 0;
	unsigned lm = 0;
	for (unsigned l=0; l<=LMAX; l++) {
		ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
		Rlm[lm] = creal(Zlm[ll]);		// m=0
		Ilm[lm] = cimag(Zlm[ll]);
		lm++;
	}
	double half_parity = 0.5;
	for (unsigned m=1; m<=MMAX; m++) {
		ll = (m-1)*m;
		half_parity = -half_parity;		// (-1)^m * 0.5
		for (unsigned l=m; l<=LMAX; l++) {
			ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
			cplx b = Zlm[ll-m] * half_parity;		// (-1)^m for m negative.
			cplx a = Zlm[ll+m] * 0.5;
			Rlm[lm] = (conj(b) + a);		// real part
			Ilm[lm] = (conj(b) - a)*I;		// imag part
			lm++;
		}
	}
}


/// complex scalar transform.
/// in: complex spatial field z.
/// out: alm[LM(shtns,l,m)] is the SH coefficients of order l and degree m (with -l <= m <= l)
/// for a total of nlm_cplx_calc(lmax,mmax,mres) coefficients.
void spat_cplx_to_SH(shtns_cfg shtns, cplx *z, cplx *alm)
{
	const long int nspat = shtns->nspat;
	cplx *rlm, *ilm, *Q, *mem;

	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	// alloc temporary fields
	mem = (cplx*) VMALLOC( (nspat+2*NLM)*sizeof(cplx) );
	rlm = mem + nspat;
	ilm = rlm + NLM;

	Q = z;
	if (NPHI>1) {
		if (shtns->fftc_mode != 0) Q = mem;			// out-of-place transform
		fftw_execute_dft(shtns->fft_cplx, z, Q);
	}

	const double norm = 1.0/NPHI;
	#pragma omp parallel for schedule(static,1) num_threads(shtns->nthreads)
	for (int m=0; m<=MMAX; m++) {
	if (m==0) {	// m=0
		spat_to_SH_ml(shtns, 0, Q,      rlm, LMAX);						// real
		spat_to_SH_ml(shtns, 0, (cplx*)(((double*)Q)+1), ilm, LMAX);	// imag
		int lm = 0;
		int ll = 0;
		for (int l=0; l<=LMAX; l++) {
			ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
			alm[ll] = (creal(rlm[lm]) + I*creal(ilm[lm]))*norm;		// m=0
			lm++;
		}
	} else {
		long lm = LM(shtns,m,m);
		spat_to_SH_ml(shtns, m, Q + (NPHI-m)*NLAT, rlm + lm, LMAX);		// m>0
		spat_to_SH_ml(shtns, m, Q + m*NLAT,        ilm + lm, LMAX);		// m<0
		int ll = (m-1)*m;
		for (int l=m; l<=LMAX; l++) {
			ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
			cplx rr = rlm[lm];	// +m
			cplx ii = ilm[lm];	// -m
			alm[ll+m] = rr*norm;
			if (m&1) ii = -ii;				// m<0, m odd
			alm[ll-m] = ii*norm;
			lm++;
		}
	}
	}

	VFREE(mem);
}

/// complex scalar transform.
/// in: alm[LM_cplx(shtns,l,m)] is the SH coefficients of order l and degree m (with -l <= m <= l)
/// for a total of nlm_cplx_calc(lmax,mmax,mres) coefficients.
/// out: complex spatial field z.
void SH_to_spat_cplx(shtns_cfg shtns, cplx *alm, cplx *z)
{
	const long int nspat = shtns->nspat;
	cplx *rlm, *ilm, *Q, *mem;

	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	// alloc temporary fields
	mem = (cplx*) VMALLOC( 2*(nspat + NLM*2)*sizeof(double) );
	rlm = mem + nspat;
	ilm = rlm + NLM;
	
	Q = z;
	if ((NPHI>1) && (shtns->fftc_mode != 0)) Q = mem;			// out-of-place transform

	#pragma omp parallel for schedule(static,1) num_threads(shtns->nthreads)
	for (int m=0; m<=MMAX; m++) {
	if (m==0) {	// m=0
		int lm = 0;
		int ll = 0;
		for (int l=0; l<=LMAX; l++) {
			ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
			rlm[lm] = creal(alm[ll]);
			ilm[lm] = cimag(alm[ll]);
			lm++;
		}
		cplx tmp[NLAT] SSE;
		SH_to_spat_ml(shtns, 0, rlm, Q,   LMAX);
		SH_to_spat_ml(shtns, 0, ilm, tmp, LMAX);
		for (int it=0; it<NLAT; it++) {
			((double*)Q)[2*it+1] = creal(tmp[it]);		// copy imaginary part to destination array.
		}
		// fill m>MMAX with zeros!
		for (long i=(MMAX+1)*NLAT; i<(NPHI-MMAX)*NLAT; i++)	 Q[i] = 0.0;
	} else {
		long lm = LM(shtns,m,m);
		int ll = m*m;
		cplx* almm = alm - 2*m;
		for (int l=m; l<=LMAX; l++) {			// gather from cplx rep
			ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
			cplx rr = alm[ll];	// +m
			cplx ii = almm[ll];	// -m
			if (m&1) ii = -ii;	// m<0, m odd
			rlm[lm] = rr;
			ilm[lm] = ii;
			lm++;
		}
		lm = LM(shtns,m,m);
		SH_to_spat_ml(shtns, m, rlm + lm, Q + m*NLAT,        LMAX);		// m>0
		SH_to_spat_ml(shtns, m, ilm + lm, Q + (NPHI-m)*NLAT, LMAX);		// m<0
	}
	}

	if (NPHI>1) fftw_execute_dft(shtns->ifft_cplx, Q, z);
	VFREE(mem);
}



/// complex vector transform (2D).
/// in: slm, tlm are the spheroidal/toroidal SH coefficients of order l and degree m (with -l <= m <= l)
/// out: zt, zp are respectively the theta and phi components of the complex spatial vector field.
void SHsphtor_to_spat_cplx(shtns_cfg shtns, cplx *slm, cplx *tlm, cplx *zt, cplx *zp)
{
	const long int nspat = shtns->nspat;
	cplx *zzt, *zzp, *mem, *stlm;

	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	// alloc temporary fields
	mem = (cplx*) VMALLOC( 4*(nspat + NLM*2)*sizeof(double) );
	stlm = mem + 2*nspat;

	zzt = zt;		zzp = zp;
	if ((NPHI>1) && (shtns->fftc_mode != 0)) {	zzt = mem;	zzp = mem + nspat;  }			// out-of-place transform

	#pragma omp parallel for schedule(static,1) num_threads(shtns->nthreads)
	for (int m=0; m<=MMAX; m++) {
		const long stride = LMAX+1-m;
		if (m==0) {	// m=0
			int lm = 0;
			int ll = 0;
			for (int l=0; l<=LMAX; l++) {			// gather from cplx rep
				ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
				stlm[lm]          = creal(slm[ll]);
				stlm[lm+2*stride] = cimag(slm[ll]);
				stlm[lm+stride]   = creal(tlm[ll]);
				stlm[lm+3*stride] = cimag(tlm[ll]);
				lm++;
			}
			cplx tt[NLAT] SSE;
			cplx pp[NLAT] SSE;
			SHsphtor_to_spat_ml(shtns, 0, stlm,          stlm+stride,   zzt, zzp, LMAX);
			SHsphtor_to_spat_ml(shtns, 0, stlm+2*stride, stlm+3*stride, tt,  pp,  LMAX);
			for (int it=0; it<NLAT; it++) {
				((double*)zzt)[2*it+1] = creal(tt[it]);		// copy imaginary part to destination array.
				((double*)zzp)[2*it+1] = creal(pp[it]);
			}
			// fill m>MMAX with zeros!
			for (long i=(MMAX+1)*NLAT; i<(NPHI-MMAX)*NLAT; i++)	{
				zzt[i] = 0.0;		zzp[i] = 0.0;
			}
		} else {
			long lm = 4 * LM(shtns,m,m);
			int ll = (m-1)*m;
			for (int l=m; l<=LMAX; l++) {			// gather from cplx rep
				ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
				cplx sr = slm[ll+m];	// +m
				cplx si = slm[ll-m];	// -m
				cplx tr = tlm[ll+m];	// +m
				cplx ti = tlm[ll-m];	// -m
				if (m&1) {				// m<0, m odd
					si = -si;
					ti = -ti;
				}
				stlm[lm] = sr;          stlm[lm+2*stride] = si;
				stlm[lm+stride] = tr;	stlm[lm+3*stride] = ti;
				lm++;
			}
			lm = 4 * LM(shtns,m,m);
			SHsphtor_to_spat_ml(shtns, m,  stlm+lm,          stlm+stride+lm,   zzt + m*NLAT,        zzp + m*NLAT,        LMAX);		// m>0
			SHsphtor_to_spat_ml(shtns, -m, stlm+2*stride+lm, stlm+3*stride+lm, zzt + (NPHI-m)*NLAT, zzp + (NPHI-m)*NLAT, LMAX);		// m<0
		}
	}

	if (NPHI>1) {
		fftw_execute_dft(shtns->ifft_cplx, zzt, zt);
		fftw_execute_dft(shtns->ifft_cplx, zzp, zp);
	}
	VFREE(mem);
}


/// complex vector transform (2D).
/// zt,zp: theta,phi components of the complex spatial vector field.
/// out: slm[LM_cplx(l,m)] and tlm[LM_cplx(l,m)] are the SH coefficients of order l and degree m (with -l <= m <= l)
/// for a total of shtns->nlm_cplx = nlm_cplx_calc(lmax, mmax, mres) coefficients.
void spat_cplx_to_SHsphtor(shtns_cfg shtns, cplx *zt, cplx *zp, cplx *slm, cplx *tlm)
{
	const long int nspat = shtns->nspat;
	cplx *zzt, *zzp, *mem, *stlm;

	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	// alloc temporary fields
	mem = (cplx*) VMALLOC( (2*nspat + NLM*4)*sizeof(cplx) );
	stlm = mem + 2*nspat;

	zzt = zt;		zzp = zp;
	if (NPHI>1) {
		if (shtns->fftc_mode != 0) {
			zzt = mem;		zzp = mem + nspat;	// out-of-place transform
		}
		fftw_execute_dft(shtns->fft_cplx, zt, zzt);
		fftw_execute_dft(shtns->fft_cplx, zp, zzp);
	}

	const double norm = 1.0/NPHI;
	#pragma omp parallel for schedule(static,1) num_threads(shtns->nthreads)
	for (int m=0; m<=MMAX; m++) {
		const long stride = LMAX+1-m;
		if (m==0) {	// m=0
			spat_to_SHsphtor_ml(shtns, 0, zzt, zzp,     stlm, stlm+stride, LMAX);	// real
			spat_to_SHsphtor_ml(shtns, 0, (cplx*)(((double*)zzt)+1), (cplx*)(((double*)zzp)+1), stlm+2*stride, stlm+3*stride, LMAX);	// imag
			int lm = 0;
			int ll = 0;
			for (int l=0; l<=LMAX; l++) {
				ll += (l<=MMAX) ? 2*l : 2*MMAX+1;
				slm[ll] = (creal(stlm[lm])        + I*creal(stlm[lm+2*stride]))*norm;		// m=0
				tlm[ll] = (creal(stlm[lm+stride]) + I*creal(stlm[lm+3*stride]))*norm;		// m=0
				lm++;
			}
		} else {
			long lm = 4 * LM(shtns,m,m);
			spat_to_SHsphtor_ml(shtns, m,  zzt + (NPHI-m)*NLAT, zzp + (NPHI-m)*NLAT, stlm+lm,          stlm+lm+stride,   LMAX);		// m>0
			spat_to_SHsphtor_ml(shtns, -m, zzt + m*NLAT,        zzp + m*NLAT,        stlm+lm+2*stride, stlm+lm+3*stride, LMAX);		// m<0
			int ll = (m-1)*m;
			for (int l=m; l<=LMAX; l++) {
				ll += (l<=MMAX) ? 2*l : 2*MMAX+1;	
				cplx sr = stlm[lm];				// +m
				cplx tr = stlm[lm+stride];		// +m
				cplx si = stlm[lm+2*stride];	// -m
				cplx ti = stlm[lm+3*stride];	// -m
				slm[ll+m] = sr*norm;
				tlm[ll+m] = tr*norm;
				if (m&1) {	si = -si;		ti = -ti;  }				// m<0, m odd
				slm[ll-m] = si*norm;
				tlm[ll-m] = ti*norm;
				lm++;
			}
		}
	}

	VFREE(mem);
}

/// complex vector transform (3D).
/// in: zr,zt,zp are the r,theta,phi components of the complex spatial vector field.
/// out: {qlm,slm,tlm}[LM_cplx(l,m)] are the SH coefficients of order l and degree m (with -l <= m <= l)
void spat_cplx_to_SHqst(shtns_cfg shtns, cplx *zr, cplx *zt, cplx *zp, cplx *qlm, cplx *slm, cplx *tlm)
{
	spat_cplx_to_SH(shtns, zr, qlm);
	spat_cplx_to_SHsphtor(shtns, zt,zp, slm,tlm);
}

/// complex vector transform (3D).
/// in: {qlm,slm,tlm}[LM_cplx(l,m)] are the SH coefficients of order l and degree m (with -l <= m <= l)
/// out: zr,zt,zp: r,theta,phi components of the complex spatial vector field.
void SHqst_to_spat_cplx(shtns_cfg shtns, cplx *qlm, cplx *slm, cplx *tlm, cplx *zr, cplx *zt, cplx *zp)
{
	SH_to_spat_cplx(shtns, qlm, zr);
	SHsphtor_to_spat_cplx(shtns, slm,tlm, zt,zp);
}



void SH_cplx_Xrotate90(shtns_cfg shtns, cplx *Qlm, cplx *Rlm)
{
	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	// alloc temporary fields
	cplx* rlm = (cplx*) VMALLOC( NLM*2*sizeof(cplx) );
	cplx* ilm = rlm + NLM;

	// extract complex coefficients corresponding to real and imag
	SH_cplx_to_2real(shtns, Qlm, rlm, ilm);

	// perform two real rotations:
	SH_Xrotate90(shtns, rlm, rlm);
	SH_Xrotate90(shtns, ilm, ilm);

	// combine back into complex coefficients
	SH_2real_to_cplx(shtns, rlm, ilm, Rlm);

	VFREE(rlm);
}

void SH_cplx_Yrotate90(shtns_cfg shtns, cplx *Qlm, cplx *Rlm)
{
	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	// alloc temporary fields
	cplx* rlm = (cplx*) VMALLOC( NLM*2*sizeof(cplx) );
	cplx* ilm = rlm + NLM;

	// extract complex coefficients corresponding to real and imag
	SH_cplx_to_2real(shtns, Qlm, rlm, ilm);

	// perform two real rotations:
	SH_Yrotate90(shtns, rlm, rlm);
	SH_Yrotate90(shtns, ilm, ilm);

	// combine back into complex coefficients
	SH_2real_to_cplx(shtns, rlm, ilm, Rlm);

	VFREE(rlm);
}

/// complex scalar rotation around Y
/// in: Qlm[l*(l+1)+m] is the SH coefficients of order l and degree m (with -l <= m <= l)
/// out: Qlm[l*(l+1)+m] is the rotated SH coefficients of order l and degree m (with -l <= m <= l)
void SH_cplx_Yrotate(shtns_cfg shtns, cplx *Qlm, double alpha, cplx *Rlm)
{
	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	// alloc temporary fields
	cplx* rlm = (cplx*) VMALLOC( NLM*2*sizeof(cplx) );
	cplx* ilm = rlm + NLM;

	// extract complex coefficients corresponding to real and imag
	SH_cplx_to_2real(shtns, Qlm, rlm, ilm);

	// perform two real rotations:
	SH_Yrotate(shtns, rlm, alpha, rlm);
	SH_Yrotate(shtns, ilm, alpha, ilm);

	// combine back into complex coefficients
	SH_2real_to_cplx(shtns, rlm, ilm, Rlm);

	VFREE(rlm);
}

/// complex scalar rotation around Z
/// in: Qlm[l*(l+1)+m] is the SH coefficients of order l and degree m (with -l <= m <= l)
/// out: Qlm[l*(l+1)+m] is the rotated SH coefficients of order l and degree m (with -l <= m <= l)
void SH_cplx_Zrotate(shtns_cfg shtns, cplx *Qlm, double alpha, cplx *Rlm)
{
	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	// alloc temporary fields
	cplx* rlm = (cplx*) VMALLOC( NLM*2*sizeof(cplx) );
	cplx* ilm = rlm + NLM;

	// extract complex coefficients corresponding to real and imag
	SH_cplx_to_2real(shtns, Qlm, rlm, ilm);

	// perform two real rotations:
	SH_Zrotate(shtns, rlm, alpha, rlm);
	SH_Zrotate(shtns, ilm, alpha, ilm);

	// combine back into complex coefficients
	SH_2real_to_cplx(shtns, rlm, ilm, Rlm);

	VFREE(rlm);
}

/*
/// Rotate a SH representation of complex field Qlm around the z-axis by angle alpha (in radians),
/// which is the same as rotating the reference frame by angle -alpha.
/// Result is stored in Rlm (which can be the same array as Qlm).
void SH_cplx_Zrotate(shtns_cfg shtns, cplx *Qlm, double alpha, cplx *Rlm)
{
	if (MRES != 1) shtns_runerr("complex SH requires mres=1.");

	cplx* eima = (cplx*) VMALLOC( (2*MMAX+1)*sizeof(cplx) );
	eima += MMAX;
	eima[0] = 1.0;
	for (int m=1; m<=MMAX; m++) {		// precompute the complex numbers
		double cma = cos(m*alpha);
		double sma = sin(m*alpha);
		eima[m]  = cma - I*sma;
		eima[-m] = cma + I*sma;
	}

	unsigned ll=0;
	for (unsigned l=0; l<=MMAX; l++) {
		for (int m=-l; m<=l; m++) {
			Rlm[ll] = Qlm[ll] * eima[m];
			ll++;
		}
	}
	for (unsigned l=MMAX+1; l<=LMAX; l++) {
		for (int m=-MMAX; m<=MMAX; m++) {
			Rlm[ll] = Qlm[ll] * eima[m];
			ll++;
		}
	}

	VFREE(eima);
}
*/

/*
void SH_to_spat_grad(shtns_cfg shtns, cplx *alm, double *gt, double *gp)
{
	double *mx;
	cplx *blm, *clm;
	
	blm = (cplx*) VMALLOC( 3*NLM*sizeof(cplx) );
	clm = blm + NLM;
	mx = (double*)(clm + NLM);

	st_dt_matrix(shtns, mx);
	SH_mul_mx(shtns, mx, alm, blm);
	int lm=0;
	for (int im=0; im<=MMAX; im++) {
		int m = im*MRES;
		for (int l=m; l<=LMAX; l++) {
			clm[lm] = alm[lm] * I*m;
			lm++;
		}
	}
	SH_to_spat(shtns, blm, gt);
	SH_to_spat(shtns, clm, gp);
	for (int ip=0; ip<NPHI; ip++) {
		for (int it=0; it<NLAT; it++) {
			gt[ip*NLAT+it] /= shtns->st[it];
			gp[ip*NLAT+it] /= shtns->st[it];
		}
	}
	VFREE(blm);
}
*/

#ifdef _GCC_VEC_
typedef double rndu __attribute__ ((vector_size (VSIZE2*8), aligned (8)));		///< \internal UNALIGNED vector that contains a complex number
typedef double v2du __attribute__ ((vector_size (16), aligned (8)));		///< \internal UNALIGNED vector that contains a complex number
#else
typedef rnd rndu;	///< \internal
typedef v2d v2du;	///< \internal
#endif

/** \addtogroup rotation Rotations of Spherical Harmonic fields.
Rotation of spherical harmonics, using an on-the-fly algorithm (does not store the rotation matrix) inspired by
the GUMEROV's algorithm to generate the Wigner-d matrices describing rotation of Spherical Harmonics.
See https://arxiv.org/abs/1403.7698  or  https://doi.org/10.1007/978-3-319-13230-3_5
Thanks to Alex J. Yuffa, ayuffa@gmail.com  for his suggestions and help.

These functions are more accurate and faster than the \link gimbutas_rotation pseudo-spectral rotations of Gimbutas\endlink and do not require mmax=lmax.
Note that if mmax<lmax, the rotations are apprixmate only, with the quality of the approximation decreasing with increasing
'beta' angle (rotation angle around Y-axis).
*/
///@{

/// Allocate memory and precompute some recurrence coefficients for rotation (independent of angle).
/// Setting mmax < lmax will result in approximate rotations if not aligned with the z-axis, while mmax=lmax leads to exact rotations.
shtns_rot shtns_rotation_create(const int lmax, const int mmax, int norm)
{
	shtns_rot r = (shtns_rot) malloc(sizeof(struct shtns_rot_));
	r->lmax = lmax;
	r->mmax = mmax;
	r->no_cs_phase = (norm & SHT_NO_CS_PHASE) ? -1. : 1.;	// adapt rotations to Condon-Shortley phase
	r->m0_renorm = (norm & SHT_REAL_NORM) ? sqrt(2.) : 1.;	// adapt to real norm
	r->plm_beta = (double*) malloc( sizeof(double) * nlm_calc(lmax+1, lmax+1, 1) );
	r->sht = shtns_create(lmax+1, lmax+1, 1, sht_for_rotations | SHT_NO_CS_PHASE);		// need SH up to lmax+1, with Schmidt semi-normalization.
	r->alpha = 0.0;
	r->beta = 0.0;
	r->gamma = 0.0;
	r->eia = 0;		r->eig = 0;
	r->flag_alpha_gamma = 0;		// mark alpha and gamma as zero
	return r;
}

void shtns_rotation_destroy(shtns_rot r)
{
	if (r) {
		shtns_destroy(r->sht);
		if (r->plm_beta) free(r->plm_beta);
		free(r);
	}
}

/// \internal retruns cos(phi) + I*sin(phi), with specialization for some particular values of phi.
static cplx special_eiphi(const double phi)
{
	cplx eip;
	if (phi == 0.0) {
		eip = 1.0;
	} else if (phi == M_PI) {
		eip = -1.0;
	} else if (phi == M_PI/2) {
		eip = I;
	} else if ((phi == 3*M_PI/2) || (phi == -M_PI/2)) {
		eip = -I;
	} else if (phi == M_PI/4) {
		eip = sqrt(0.5) + I*sqrt(0.5);
	} else if (phi == M_PI/3) {
		eip = 0.5 + I*sqrt(3)*0.5;
	} else {
		eip = cos(phi) + I*sin(phi);
	}
	return eip;
}

/// Set the rotation angles, and compute associated Legendre functions, given the 3 intrinsic Euler angles in ZYZ convention.
void shtns_rotation_set_angles_ZYZ(shtns_rot r, double alpha, double beta, double gamma)
{
	beta *= r->no_cs_phase;		// condon-shortley phase is the same thing as rotating by 180°, or changing sign of beta.
	if UNLIKELY(fabs(beta) > M_PI) {
		printf("ERROR: angle 'beta' must be between -pi and pi\n");
		exit(1);
	}
	if (beta < 0.0) {	// translate to beta>0 as beta<0 is not supported.
		alpha = (alpha>0) ? alpha-M_PI : alpha+M_PI;		// rotate by 180°
		beta = fabs(beta);
		gamma = (gamma>0) ? gamma-M_PI : gamma+M_PI;		// rotate by 180°
	} else if (beta == 0.0) {
		alpha += gamma;
		gamma = 0.0;
	}

	// step 0 : compute plm(beta)
	const cplx eib = special_eiphi(beta);
	r->cos_beta = creal(eib);
	r->sin_beta = cimag(eib);
	r->eia = special_eiphi(-alpha);
	r->eig = special_eiphi(-gamma);
	r->alpha = alpha;
	r->beta = beta;
	r->gamma = gamma;
	r->flag_alpha_gamma = (alpha != 0) + 2*(gamma != 0);
	if (beta != 0.0) {
		const int lmax = r->lmax + 1;			// need SH up to lmax+1.
		#pragma omp parallel for schedule(dynamic) firstprivate(lmax)
		for (int m=0; m<=lmax; m++) {
			const long ofs = m*(lmax+2) - (m*(m+1))/2;
			legendre_sphPlm_array(r->sht, lmax, m, r->cos_beta, r->plm_beta + ofs);
		}
	}
}

/// Set the rotation angles, and compute associated Legendre functions, given the 3 intrinsic Euler angles in ZXZ convention.
void shtns_rotation_set_angles_ZXZ(shtns_rot r, double alpha, double beta, double gamma)
{
	shtns_rotation_set_angles_ZYZ(r, alpha+M_PI/2, beta, gamma-M_PI/2);
}

/// Sets a rotation by angle theta around axis of cartesian coordinates (Vx,Vy,Vz), and compute associated Legendre functions.
void shtns_rotation_set_angle_axis(shtns_rot r, double theta, double Vx, double Vy, double Vz)
{
	if ((Vx==0) && (Vy==0)) {	// rotation along Z-axis
		if (Vz<0) theta = -theta;
		shtns_rotation_set_angles_ZYZ(r, theta, 0, 0);
	} else {
		// 1) convert to normalized quaternion representation (c, Vx, Vy, Vz). See https://en.wikipedia.org/wiki/Versor
		double s = sin(0.5*theta);
		double c = cos(0.5*theta);
		double n = s / sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
		Vx *= n;	Vy *= n;	Vz *= n;

		// 2) convert from quaternion to extrinsic Euler angles:
		double beta = acos( 1.0 - 2.0*(Vx*Vx + Vy*Vy) );  // = acos( c*c + Vz*Vz - (Vx*Vx + Vy*Vy) );
		double Vxz = Vx*Vz;
		double Vyz = Vy*Vz;
		double cVx = c*Vx;
		double cVy = c*Vy;
		double alpha = atan2( Vyz - cVx, cVy + Vxz );		// note: for ZXZ convention: switch x with y and alpha with gamma.
		double gamma = atan2( Vyz + cVx, cVy - Vxz );
		shtns_rotation_set_angles_ZYZ(r, gamma, beta, alpha);	// extrinsic rotation: swap gamma and alpha. See https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations
	}
}

/// \internal lw is the line-width. Use lw=2*l+1 for the full matrix, or l+1 for a compressed matrix.
/// It always has 2*l+1 lines.
/// \return the number of columns in the matrix (l+1 for compressed or 2*l+1 for the full matrix; 0 if l is out of the 0..lmax range.
int quarter_wigner_d_matrix(shtns_rot r, const int l, double* mx, const int compressed)
{
	if ((l > r->lmax) || (l<0)) {
		printf("ERROR: 0 <= l <= lmax not satified.\n");
		return 0;	// error, nothing written into mx.
	}

	const int lmax = r->lmax + 1;
	const double cos_beta = r->cos_beta;
	const double sin_beta = r->sin_beta;
	const double* const plm_beta = r->plm_beta;

	const int lw = (compressed) ? l+1 : 2*l+1;
	// shift mx to index by negative values:
	mx += l*lw + l*(compressed==0);

	mx[0] = plm_beta[l];		// d(m'=0, m=0)
	// step 2+3:  d(m'=0,m) and d(m'=1,m)
	double cb_p1 = (1. + cos_beta)*0.5;
	double cb_m1 = (1. - cos_beta)*0.5;
	const double d_1 = 1.0/((l)*(l+1.));
	for (int m=1; m<=l; m++) {
		const long ofs_m = m*(lmax+2) - (m*(m+1))/2;
		const long ofs_mm1 = (m-1)*(lmax+2) - ((m-1)*m)/2;
		const long ofs_mp1 = (m+1)*(lmax+2) - ((m+1)*(m+2))/2;
		mx[m] = plm_beta[ofs_m + (l-m)];			// m'=0 : copy plm values (eq 32)

		double a  = sqrt( ((l+1+m)*(l+1-m)) * d_1 ) * sin_beta;
		double b1 = sqrt( ((l-m+1)*(l-m+2)) * d_1 ) * cb_p1;
		double b2 = sqrt( ((l+m+1)*(l+m+2)) * d_1 ) * cb_m1;
		mx[lw + m] = b2*plm_beta[ofs_mp1 + (l+1)-(m+1)]  +  b1*plm_beta[ofs_mm1 + (l+1)-(m-1)]   + a*plm_beta[ofs_m + (l+1)-m];		// m'=1 (eq 105)
	}

	double clm[l+1];		// temporary array holding clm
	for (int m=0; m<l; m++)  clm[m] = sqrt( (l-m)*(l+m+1) );		// precompute clm
	clm[l] = 0.0;		// boundary condition handled with this

	// step 5:	recursively compute d(m',m),  m'=-1; m=-m'..l
	{ const int mp=-1;	// m'=-1
		const double c_1 = 1.0/clm[0];
		for (int m=-mp; m<l; m++) {
			double H = mx[lw*(mp+2) + m] + c_1*( -clm[m-1] * mx[lw*(mp+1) + (m-1)] + clm[m] * mx[lw*(mp+1) + (m+1)]);	// eq 108
			mx[lw*mp + m] = H;		// d(m',m)
		}
		const int m=l;
			double H = mx[lw*(mp+2) + m] - c_1*clm[m-1] * mx[lw*(mp+1) + (m-1)];	// eq 108
			mx[lw*mp + m] = H;		// d(m',m)
	}

	// step 4+5 merged:	recursively compute d(m',m),  m'=2..l; m=m'..l  AND m'=-2..-l; m=-m'..l
	for (int mp=2; mp<=l; mp++) {
		const double c_1 = 1.0/clm[mp-1];
		const double cmp2 = clm[mp-2];
		for (int m=mp; m<l; m++) {
			double H  = cmp2 * mx[lw*(mp-2) + m]  + clm[m-1] * mx[lw*(mp-1) + (m-1)]  - clm[m] * mx[lw*(mp-1) + (m+1)];		// eq 106
			double H_ = cmp2 * mx[-lw*(mp-2) + m] - clm[m-1] * mx[-lw*(mp-1) + (m-1)] + clm[m] * mx[-lw*(mp-1) + (m+1)];	// eq 108
			mx[lw*mp + m]  = H  * c_1;		// d(m',m)
			mx[-lw*mp + m] = H_ * c_1;		// d(-m',m)
		}
		const int m=l;
			double H  = cmp2 * mx[lw*(mp-2) + m]  + clm[m-1] * mx[lw*(mp-1) + (m-1)];	// eq 106
			double H_ = cmp2 * mx[-lw*(mp-2) + m] - clm[m-1] * mx[-lw*(mp-1) + (m-1)];	// eq 108
			mx[lw*mp + m]  = H  * c_1;		// d(m',m)
			mx[-lw*mp + m] = H_ * c_1;		// d(-m',m)
	}

	return lw;	// return the number of lines (or columns) in the matrix.
}


/// Generate spherical-harmonic rotation matrix for given degree l and orthonormal convention (Wigner-d matrix)
/// \param[out] mx is an (2*l+1)*(2*l+1) array that will be filled with the Wigner-d matrix elements (rotation matrix along Y-axis in orthonormal spherical harmonic space).
/// \return 0 if error, or 2*l+1 (size of the square matrix) otherwise.
/// \warning The returned rotation matrix only applies to **orthonormal** convention (see \ref norm)
int shtns_rotation_wigner_d_matrix(shtns_rot r, const int l, double* mx)
{
	// step 1:
	if (l==0) {
		mx[0] = 1;
		return 1;
	}

	const int lw = quarter_wigner_d_matrix(r, l, mx, 0);		// generate quarter of wigner-d matrix (but full storage)
	if (lw <= 0) return 0;	// error

	// shift mx to index by negative values:
	mx += l*(2*l+1) + l;

	// step 6: fill matrix using symmetries
	for (int m=1; m<=l; m++) {		// diagonals
		mx[(2*l+1)*m + -m]   = mx[(2*l+1)*(-m) + m];
		mx[(2*l+1)*(-m) - m] = mx[(2*l+1)*m + m];
	}
	for (int mp=-l+1; mp<l; mp++) {		// off-diagonals
		for (int m=abs(mp)+1; m<=l; m++) {
			double x = mx[(2*l+1)*mp + m];
			double parity = (1-2*((m-mp)&1));
			mx[(2*l+1)*(-m) - mp] = x;
			mx[(2*l+1)*m + mp]    = x * parity;
			mx[(2*l+1)*(-mp) - m] = x * parity;
		}
	}
	return lw;
}

/*
/// rotate Zlm, store the result in Rlm. Generates quarter of Wigner-d matrices internally.
/// Reference implementation.
void shtns_rotation_apply_cplx(shtns_rot r, cplx* Zlm, cplx* Rlm)
{
	const int lmax = r->lmax;

	Rlm[0] = Zlm[0];		// copy l=0

	#pragma omp for schedule(dynamic)
	for (int l=1; l<=lmax; l++) {
		const int lw = l+1;		// line width of compressed matrix.
		double* mx0 = (double*) malloc( sizeof(double) * (2*l+1)*lw );
		memset(mx0, 0, sizeof(double)*(2*l+1)*lw);

		quarter_wigner_d_matrix(r, l, mx0, 1);		// 1 : compressed matrix (saves memory and time)
		const double* mx = mx0 + l*lw;		// allow indexing by negative m and m'

		const cplx* zl = Zlm + l*(l+1);
		cplx* rl = Rlm + l*(l+1);
		// apply the matrix
		for (int m=-l; m<=l; m++) {
			cplx rm = 0.0;
			// only 1/4 of the matrix is read (4 times).
			for (int mp=l; mp>=abs(m); mp--) {			// third domain (right) : populated.
				rm += zl[mp] * mx[m*lw + mp];
			}
			if (m<0) {
				for (int mp=-m-1; mp>=m; mp--) {		// second domain (top)
					rm += zl[mp] * mx[-mp*lw - m];		// exchange m and m' and change signs  => no parity factor
				}
			} else {
				for (int mp=m-1; mp>=-m; mp-=2) {		// second domain (bottom) => always even number of elements
					rm -= zl[mp] * mx[mp*lw + m];		// exchange m and m' (negative parity)
					rm += zl[mp-1] * mx[(mp-1)*lw + m];	// exchange m and m' (positive parity)
				}
			}
			for (int mp=-abs(m)-1; mp >-l; mp-=2) {		// first domain (left)
				rm -= zl[mp] * mx[-m*lw - mp];			// change sign of m and m'
				rm += zl[mp-1] * mx[-m*lw - (mp-1)];	// change sign of m and m'
			}
			if ((l+m)&1) {	int mp = -l;	// boundary condition.
				rm -= zl[mp] * mx[-m*lw - mp];			// change sign of m and m'				
			}
			rl[m] = rm;
		}

		free(mx0);
	}
}
*/

void shtns_rotation_apply_real(shtns_rot r, cplx* Qlm, cplx* Rlm)
{
	const int lmax = r->lmax+1;
	const int mmax = r->mmax;

	if (r->beta == 0.0) {	// only rotation along Z-axis.
		if (Rlm != Qlm)		for (int l=0; l<lmax; l++) 	Rlm[l] = Qlm[l];		// copy m=0
		long lm = lmax;
		const cplx ei_alpha = r->eia;
		cplx eim_alpha = ei_alpha;
		for (int m=1; m<=mmax; m++) {
			for (int l=m; l<lmax; l++) {
				Rlm[lm] = Qlm[lm] * eim_alpha;
				lm++;
			}
			eim_alpha *= ei_alpha;
		}
		return;
	}

	Rlm[0] = Qlm[0];		// copy l=0 (invariant by rotation)

	#pragma omp parallel firstprivate(lmax,mmax)
	{
		const double cos_beta = r->cos_beta;
		const double sin_beta = r->sin_beta;
		const double* const plm_beta = r->plm_beta;
		const int flag_ag = r->flag_alpha_gamma;

		const double cb_p1 = (1. + cos_beta)*0.5;
		const double cb_m1 = (1. - cos_beta)*0.5;

		const double m0_renorm = r->m0_renorm;
		const double m0_renorm_1 = 1.0 / m0_renorm;

		const cplx eia = r->eia;
		const cplx eig = r->eig;

//		cplx* const rl = (cplx*) malloc(sizeof(double) * (4*(l+1) + 5*lw +2));	// 2 real temp array 2*(l+1), + 5 lines of storage (lw) + a bit
		cplx* const rl = (cplx*) malloc(sizeof(double) * (9*lmax +2));	// 2 real temp array 2*(l+1), + 5 lines of storage (lw) + a bit

		#pragma omp for schedule(dynamic,1)
		for (int l=lmax-1; l>0; l--) {
			const int lw = l+1;		// line width of matrix.
			const int mlim = (l <= mmax) ? l       : mmax;
			cplx* const ql = rl + (l+1);
			double* const m0 = (double*) (rl + 2*(l+1));
			memset(rl, 0, sizeof(cplx)*(2*l+1));		// zero the temp destination array.

			// gather all m's for given l of source array:
			long lm = l;
			ql[0] = creal(Qlm[lm]) * m0_renorm;			// m=0
			if (flag_ag & 1) {		// pre-rotate along Z-axis
				cplx eima = eia;
				for (int m=1; m<=mlim; m++) {
					lm += lmax-m;
					ql[m] = Qlm[lm] * eima;
					eima *= eia;
				}
			} else {		// copy
				for (int m=1; m<=mlim; m++) {
					lm += lmax-m;
					((v2d*)ql)[m] = ((v2d*)Qlm)[lm];
				}
			}
			for (int m=mlim+1; m<=l; m++) ql[m] = 0;		// zero out m that are absent.

			m0[0] = plm_beta[l];		// d(m'=0, m=0)
			m0[lw] = 0.0;				// required, as a boundary condition
			double r0 = 0.5*creal(ql[0])*m0[0];		// accumulator for rl[0]

			// step 2+3:  d(m'=0,m) and d(m'=1,m)
			const double d_1 = 1.0/((l)*(l+1.));
			cplx r1 = 0.0;		// accumulator for rl[1]
			double parity = -1;
			long ofs_m = (lmax+2) - 1;	// m*(lmax+2) - (m*(m+1))/2  for m=1
			long ofs_mm1 = 0;			// m*(lmax+2) - (m*(m+1))/2  for m=0
			for (int m=1; m<=l; m++) {
				double H = plm_beta[ofs_m + (l-m)];			// m'=0 : copy plm values (eq 32)
				m0[m] = H;
				r0 += creal(ql[m]) * H;		// right quadrant (m>0)
				const double a  = sqrt( ((l+1+m)*(l+1-m)) * d_1 ) * sin_beta;
				rl[m] += creal(ql[0]) * H*parity;			// bottom quadrant (m>0) : exchange m and m' => parity factor
				const double b1 = sqrt( ((l-m+1)*(l-m+2)) * d_1 ) * cb_p1;
				long ofs_mp1 = (m+1)*(lmax+2) - ((m+1)*(m+2))/2;
				const double b2 = sqrt( ((l+m+1)*(l+m+2)) * d_1 ) * cb_m1;
				parity = -parity;
				H = b2*plm_beta[ofs_mp1 + (l+1)-(m+1)]  +  b1*plm_beta[ofs_mm1 + (l+1)-(m-1)]   + a*plm_beta[ofs_m + (l+1)-m];		// m'=1 (eq 105)
				m0[lw+m] = H;		// m'=1
				r1 += ql[m] * H;		// right quadrant (m>0)
				ofs_mm1 = ofs_m;
				ofs_m = ofs_mp1;
				if (m>1) {	// avoid duplicates
					rl[m]  += ql[1] * H*parity;	// bottom quadrant (m>0) : exchange m and m' => parity factor
				}
			}
			rl[0] = r0+r0;
			rl[1] += r1;

			double* clm = m0 + 4*lw +1;		// array holding clm (size l+2, adressing by m=-1 to l)
			clm[-1] = 0.0;
			for (int m=0; m<l; m++)  clm[m] = sqrt( (l-m)*(l+m+1) );		// precompute clm
			clm[l] = 0.0;		// boundary condition handled with this

			// step 5:	recursively compute d(m',m),  for m'=-1;  m=1..l
			double* mx1_ = m0 + 2*lw;
			double* mx0_ = m0 + 3*lw;
			memcpy(mx1_, m0, sizeof(double)*2*lw);		// first copy the initial lines (m'=0 and m'=1).
			#if _GCC_VEC_
			s2d conj_parity = {-1.0, 1.0};		// change sign of real or imaginary part only.
			#endif
			{ const int mp = -1;
				double clm_1 = clm[0];
				const double c_1 = 1.0/clm_1;  // 1.0/clm[l+mp];
				double mx1_1 = mx1_[-mp-1];
				double mx1_0 = mx1_[-mp];
				v2d rmp = vdup(0.0);
				#if _GCC_VEC_
				v2d zlmp = ((v2d*)ql)[-mp] * conj_parity;
				#else
				cplx zlmp = conj(ql[-mp]) * (1-2*(mp&1));
				#endif
				int m = -mp;
				for (; m<l; m+=2) {
					double clm0 = clm[m];
					double clm1 = clm[m+1];
					double mx11 = mx1_[m+1];
					double mx12 = mx1_[m+2];
					double H0 = mx0_[m]   + c_1*( -clm_1 * mx1_1 + clm0 * mx11 );	// eq 108
					double H1 = mx0_[m+1] + c_1*(  -clm0 * mx1_0 + clm1 * mx12 );
					clm_1 = clm1;		mx1_1 = mx11;	mx1_0 = mx12;		// cycle coefficients and matrix elements.
					mx0_[m] = H0;			// d(m',m) -> overwrite d(m'+2,m)
					mx0_[m+1] = H1;			// d(m',m) -> overwrite d(m'+2,m)
					rmp += ((v2d*)ql)[m] * vdup(H0)  + ((v2d*)ql)[m+1] * vdup(H1);	// left quadrant, change signs => parity factor
					if (m>-mp) {	// avoid duplicates
						((v2d*)rl)[m]  += zlmp * vdup(H0);	// bottom quadrant (m>0) : exchange m and m' => parity factor
					}
					((v2d*)rl)[m+1]  -= zlmp * vdup(H1);	// bottom quadrant (m>0) : exchange m and m' => parity factor
				}
				if (m==l) {
					double H0 = mx0_[m]   - c_1 * clm_1 * mx1_1;	// eq 108
					mx0_[m] = H0;			// d(m',m) -> overwrite d(m'+2,m)
					rmp += ((v2d*)ql)[m] * vdup(H0);	// left quadrant, change signs => parity factor
					if (m>-mp) {	// avoid duplicates
						((v2d*)rl)[m]  += zlmp * vdup(H0);	// bottom quadrant (m>0) : exchange m and m' => parity factor
					}
				}
				#if _GCC_VEC_
				((v2d*)rl)[-mp] += rmp * conj_parity;
				conj_parity = - conj_parity;	// switch parity
				#else
				rl[-mp] += conj(rmp)*(1-2*(mp&1));		// -mp > 0
				#endif
				double* t = mx0_;		mx0_ = mx1_;	mx1_ = t;	// cycle the buffers for lines.
			}

			// step 4 + 5 merged:	recursively compute and apply d(m',m),  m'=2..l; m=m'..l  AND  m'=-2..-l; m=-m'..l
			double* mx0 = m0;
			double* mx1 = m0 + lw;
			for (int mp=2; mp<=l; mp++) {
				const double cmp2 = clm[mp-2];
				double clm_1 = clm[mp-1];
				const double c_1 = 1.0/clm_1;

				v2d rmp = vdup(0.0);
				v2d zlmp = ((v2d*)ql)[mp];
				v2d rmp_ = vdup(0.0);
				#if _GCC_VEC_
				v2d zlmp_ = zlmp * conj_parity;
				#else
				cplx zlmp_ = conj(zlmp) * (1-2*(mp&1));
				#endif

				int m = mp;
				for (; m<=l-(VSIZE2-1); m+=VSIZE2) {
					rnd clm_1 = *((rndu*)(clm+m-1));
					rnd clm0 =  *((rndu*)(clm+m));

					rnd mx1_1 =  *((rndu*)(mx1+m-1));
					rnd mx11 =  *((rndu*)(mx1+m+1));

					rnd mx1__1 =  *((rndu*)(mx1_+m-1));
					rnd mx11_  =  *((rndu*)(mx1_+m+1));

					rnd mx00  = *((rndu*)(mx0+m));
					rnd mx00_ = *((rndu*)(mx0_+m));

					rnd H  = (vall(cmp2) * mx00   + clm_1 * mx1_1  - clm0 * mx11)  * vall(c_1);		// eq 106
					rnd H_ = (vall(cmp2) * mx00_  - clm_1 * mx1__1 + clm0 * mx11_) * vall(c_1);		// eq 108
					*((rndu*)(mx0+m)) = H;
					*((rndu*)(mx0_+m)) = H_;
				}
		/*		v2d mx1_1 = *((v2d*)(mx1+mp-1));
				v2d mx1__1 = *((v2d*)(mx1_+mp-1));
				for (; m<l; m+=2) {
					v2d clm_1 = *((v2du*)(clm+m-1));
					v2d clm0 =  *((v2du*)(clm+m));

					v2d mx11 =  *((v2du*)(mx1+m+1));

					v2d mx11_  =  *((v2du*)(mx1_+m+1));

					v2d mx00  = *((v2du*)(mx0+m));
					v2d mx00_ = *((v2du*)(mx0_+m));

					v2d H = vdup(cmp2) * mx00   + clm_1 * mx1_1 - clm0 * mx11;		// eq 106
					v2d H_ = vdup(cmp2) * mx00_  - clm_1 * mx1__1 + clm0 * mx11_;	// eq 108
					mx1__1 = mx11_;		mx1_1 = mx11;
					H *= vdup(c_1);
					H_ *= vdup(c_1);
					*((v2du*)(mx0+m)) = H;
					*((v2du*)(mx0_+m)) = H_;
				}	*/
				for (; m<=l; m++) {
					double clm_1 = clm[m-1];
					double clm0 = clm[m];
					double mx1_1 = mx1[m-1];
					double mx1__1 = mx1_[m-1];
					double mx11 = mx1[m+1];
					double mx11_ = mx1_[m+1];
					double H0 = cmp2 * mx0[m]   + clm_1 * mx1_1  - clm0 * mx11;		// eq 106
					double H0_ = cmp2 * mx0_[m] - clm_1 * mx1__1 + clm0 * mx11_;	// eq 108
					H0 *= c_1;			//mx[lw*mp + m] = H * c_1;		// d(m',m)
					H0_ *= c_1;			//mx[lw*mp + m] = H * c_1;		// d(m',m)
					mx0[m] = H0;			// d(m',m) -> overwrite d(m'-2,m)
					mx0_[m] = H0_;			// d(m',m) -> overwrite d(m'+2,m)
				}

				m = mp;
				for (; m<l; m+=2) {
					double H0 = mx0[m];		double H1 = mx0[m+1];
					double H0_ = mx0_[m];	double H1_ = mx0_[m+1];
					rmp  += ((v2d*)ql)[m] * vdup(H0)  + ((v2d*)ql)[m+1] * vdup(H1);	//mx[mp*lw + m];		// right quadrant (m>0, mp>0)
					rmp_ += ((v2d*)ql)[m] * vdup(H0_) + ((v2d*)ql)[m+1] * vdup(H1_);	// left quadrant, change signs => parity factor
					if (m>mp) {	// avoid duplicates
						((v2d*)rl)[m]  += zlmp * vdup(H0)  +   zlmp_ * vdup(H0_);	// bottom quadrant (m>0) : exchange m and m' => parity factor
					}
					((v2d*)rl)[m+1] -= zlmp * vdup(H1) + zlmp_ * vdup(H1_);	// bottom quadrant (m>0) : exchange m and m' => parity factor
				}
				if (m==l) {
					double H0 = mx0[m];			// d(m',m) -> overwrite d(m'-2,m)
					double H0_ = mx0_[m];			// d(m',m) -> overwrite d(m'+2,m)
					rmp  += ((v2d*)ql)[m] * vdup(H0);	//mx[mp*lw + m];		// right quadrant (m>0, mp>0)
					rmp_ += ((v2d*)ql)[m] * vdup(H0_);	// left quadrant, change signs => parity factor
					if (m>mp) {	// avoid duplicates
						((v2d*)rl)[m]  += zlmp * vdup(H0) + zlmp_ * vdup(H0_);	// bottom quadrant (m>0) : exchange m and m' => parity factor
					}
				}
				#if _GCC_VEC_
				((v2d*)rl)[mp] += rmp + rmp_ * conj_parity;
				conj_parity = - conj_parity;	// switch parity
				#else
				rl[mp] += rmp + conj(rmp_)*(1-2*(mp&1));		// -mp > 0
				#endif
				double* t = mx0;		mx0 = mx1;		mx1 = t;	// cycle the buffers for lines.
				double* t_ = mx0_;		mx0_ = mx1_;	mx1_ = t_;	// cycle the buffers for lines.
			}

			// scatter all m's for current l into dest array:
			lm = l;
			Rlm[lm] = creal(rl[0]) * m0_renorm_1;			// m=0
			if (flag_ag & 2) {		// post-rotate along Z-axis
				cplx eimg = eig;
				for (int m=1; m<=mlim; m++) {
					lm += lmax-m;
					Rlm[lm] = rl[m] * eimg;
					eimg *= eig;
				}
			} else {		// copy values to dest
				for (int m=1; m<=mlim; m++) {
					lm += lmax-m;
					((v2d*)Rlm)[lm] = ((v2d*)rl)[m];
				}
			}
		}

		free(rl);
	}
}

/// rotate Zlm, store the result in Rlm, without ever storing the wigner-d matrices (on-the-fly operation)
void shtns_rotation_apply_cplx(shtns_rot r, cplx* Zlm, cplx* Rlm)
{
	const int lmax = r->lmax+1;
	const int mmax = r->mmax;

// indexing scheme:
// if (l<=MMAX) : l*(l+1) + m
// if (l>=MMAX) : l*(2*mmax+1) - mmax*mmax + m  = mmax*(2*l-mmax) + l+m

	if (r->beta == 0.0) {	// only rotation along Z-axis.
		long lm = 0;
		const cplx eia = r->eia;
		for (int l=0; l<lmax; l++) {
			long ll = (l<=mmax) ? l*(l+1) : mmax*(2*l-mmax) + l;
			cplx eim_alpha = eia;
			Rlm[ll + 0] = Zlm[ll + 0];	// copy m=0;
			for (int m=1; m<=l; m++) {
				Rlm[ll - m] = Zlm[ll - m] * conj(eim_alpha);
				Rlm[ll + m] = Zlm[ll + m] * eim_alpha;
				eim_alpha *= eia;
			}
		}
		return;
	}

	Rlm[0] = Zlm[0];		// copy l=0 (invariant by rotation)

	#pragma omp parallel firstprivate(lmax,mmax)
	{
		const double cos_beta = r->cos_beta;
		const double sin_beta = r->sin_beta;
		const double* const plm_beta = r->plm_beta;

		const double cb_p1 = (1. + cos_beta)*0.5;
		const double cb_m1 = (1. - cos_beta)*0.5;

		const int flag_ag = r->flag_alpha_gamma;
		const cplx eia = r->eia;
		const cplx eig = r->eig;
		const int ntmp = 1 + (flag_ag & 1);

//		v2d* rl = (v2d*) malloc(sizeof(double) * (2*ntmp*(2*l+1) + 5*lw +2));	// 1 cplx temp array (2l+1), + 5 lines of storage (lw) + a bit
		v2d* const buf = (v2d*) malloc(sizeof(double) * (2*ntmp*(2*lmax-1) + 5*lmax +2));	// 1 cplx temp array (2l+1), + 5 lines of storage (lw) + a bit

		#pragma omp for schedule(dynamic,1)
		for (int l=lmax-1; l>0; l--) {
			const int lw = l+1;		// line width of matrix.
			const int mlim = (l <= mmax) ? l       : mmax;
			const long ll  = (l <= mmax) ? l*(l+1) : mmax*(2*l-mmax) + l;
			v2d* rl = buf;
			double* const m0 = (double*) (buf + ntmp*(2*l+1));
			memset(rl, 0, sizeof(cplx)*(2*l+1));		// zero the temp dest array.
			rl += l;	// shift pointer to allow indexing by m = -l to l
			v2d* zl = (v2d*)Zlm + l*(l+1);		// source array, index by mp = -l to l, assumed aligned for sse2
			if (flag_ag & 1) {		// pre-rotate arount Z-axis
				zl = rl + 2*l+1;
				((cplx*)zl)[0] = Zlm[l*(l+1) + 0];
				cplx eim_alpha = eia;
				for (int m=1; m<=mlim; m++) {
					((cplx*)zl)[-m] = Zlm[ll - m] * conj(eim_alpha);
					((cplx*)zl)[m]  = Zlm[ll + m] * eim_alpha;
					eim_alpha *= eia;
				}
			} else if (mmax < l) {		// we must copy the source to pad it with zeros
				zl = rl + 2*l+1;
				for (int m=-mmax; m<=mmax; m++) ((cplx*)zl)[m] = Zlm[ll+m];
			}
			for (int m=mlim+1; m<=l; m++) {		// pad source with zeros
				((cplx*)zl)[-m] = 0;
				((cplx*)zl)[m]  = 0;
			}

			m0[0] = plm_beta[l];		// d(m'=0, m=0)
			m0[lw] = 0.0;				// required, as a boundary condition
			v2d r0 = zl[0]*vdup(m0[0]);		// accumulator for rl[0]

			// step 2+3:  d(m'=0,m) and d(m'=1,m)
			const double d_1 = 1.0/((l)*(l+1.));
			v2d r1 = vdup(0.0);		// accumulator for rl[1]
			v2d rm1 = vdup(0.0);		// accumulator for rl[-1]
			double parity = -1;
			long ofs_m = (lmax+2) - 1;	// m*(lmax+2) - (m*(m+1))/2  for m=1
			long ofs_mm1 = 0;			// m*(lmax+2) - (m*(m+1))/2  for m=0
			for (int m=1; m<=l; m++) {
				double H = plm_beta[ofs_m + (l-m)];			// m'=0 : copy plm values (eq 32)
				m0[m] = H;
				r0 += zl[m] * vdup(H);		// right quadrant (m>0)
				r0 += zl[-m] * vdup(H*parity);		// left quadrant (m<0) : change signs => parity factor
				const double a  = sqrt( ((l+1+m)*(l+1-m)) * d_1 ) * sin_beta;
				rl[-m] += zl[0] * vdup(H);		// top quadrant (m<0) : exchange m and m' and change signs  => no parity factor
				rl[m] += zl[0] * vdup(H*parity);			// bottom quadrant (m>0) : exchange m and m' => parity factor
				const double b1 = sqrt( ((l-m+1)*(l-m+2)) * d_1 ) * cb_p1;
				long ofs_mp1 = (m+1)*(lmax+2) - ((m+1)*(m+2))/2;
				const double b2 = sqrt( ((l+m+1)*(l+m+2)) * d_1 ) * cb_m1;
				parity = -parity;
				H = b2*plm_beta[ofs_mp1 + (l+1)-(m+1)]  +  b1*plm_beta[ofs_mm1 + (l+1)-(m-1)]   + a*plm_beta[ofs_m + (l+1)-m];		// m'=1 (eq 105)
				m0[lw+m] = H;
				r1 += zl[m] * vdup(H);		// right quadrant (m>0)
				ofs_mm1 = ofs_m;
				rm1 += zl[-m] * vdup(H*parity);		// left quadrant: change signs => parity factor
				ofs_m = ofs_mp1;
				if (m>1) {	// avoid duplicates
					rl[-m] += zl[-1] * vdup(H);	// top quadrant (m<0) : exchange m and m' and change signs => no parity factor
					rl[m]  += zl[1] * vdup(H*parity);	// bottom quadrant (m>0) : exchange m and m' => parity factor
				}
			}
			rl[-1] += rm1;
			rl[0] = r0;
			rl[1] += r1;

			double* clm = m0 + 4*lw +1;		// array holding clm (size l+2, adressing by m=-1 to l)
			clm[-1] = 0.0;		// also for BC.
			for (int m=0; m<l; m++)  clm[m] = sqrt( (l-m)*(l+m+1) );		// precompute clm
			clm[l] = 0.0;		// boundary condition handled with this

			// step 5:	recursively compute d(m',m),  for m'=-1;  m=1..l
			double* mx1_ = m0 + 2*lw;
			double* mx0_ = m0 + 3*lw;
			memcpy(mx1_, m0, sizeof(double)*2*lw);		// first copy the initial lines (m'=0 and m'=1).
			{ const int mp = -1;		// m'=-1
				double clm_1 = clm[-mp-1];
				const double c_1 = 1.0/clm_1;  // 1.0/clm[l+mp];
				double mx1_1 = mx1_[-mp-1];
				double mx1_0 = mx1_[-mp];
				v2d rmp = vdup(0.0);
				v2d rmmp = vdup(0.0);
				v2d zlmmp = zl[-mp];
				v2d zlmp = zl[mp];
				int m = -mp;
				for (; m<l; m+=2) {
					double clm0 = clm[m];
					double clm1 = clm[m+1];
					double mx11 = mx1_[m+1];
					double mx12 = mx1_[m+2];
					double H0 = mx0_[m]   + c_1*( - clm_1 * mx1_1 + clm0 * mx11);	// eq 108
					double H1 = mx0_[m+1] + c_1*( -  clm0 * mx1_0 + clm1 * mx12);
					clm_1 = clm1;		mx1_1 = mx11;	mx1_0 = mx12;		// cycle coefficients and matrix elements.
					mx0_[m] = H0;			// d(m',m) -> overwrite d(m'+2,m)
					mx0_[m+1] = H1;			// d(m',m) -> overwrite d(m'+2,m)
					rmp += zl[m] * vdup(H0)  + zl[m+1] * vdup(H1);	//mx[mp*lw + m];		// right quadrant (m>0, mp<0)
					rmmp += zl[-m] * vdup(H0)  - zl[-m-1] * vdup(H1);	// left quadrant, change signs => parity factor
					rl[-m-1] += zlmmp * vdup(H1);	// top quadrant (-m<0) : exchange m and m' and change signs => no parity factor
					if (m>-mp) {	// avoid duplicates
						rl[-m] += zlmmp * vdup(H0);	// top quadrant (-m<0) : exchange m and m' and change signs => no parity factor
						rl[m]  += zlmp * vdup(H0);	// bottom quadrant (m>0) : exchange m and m' => parity factor
					}
					rl[m+1]  -= zlmp * vdup(H1);	// bottom quadrant (m>0) : exchange m and m' => parity factor
				}
				if(m==l) {
					double H0 = mx0_[m]   - c_1 * clm_1 * mx1_1;	// eq 108
					mx0_[m] = H0;			// d(m',m) -> overwrite d(m'+2,m)
					rmp += zl[m] * vdup(H0);	//mx[mp*lw + m];		// right quadrant (m>0, mp<0)
					rmmp += zl[-m] * vdup(H0);	// left quadrant, change signs => parity factor
					if (m>-mp) {	// avoid duplicates
						rl[-m] += zlmmp * vdup(H0);	// top quadrant (-m<0) : exchange m and m' and change signs => no parity factor
						rl[m]  += zlmp * vdup(H0);	// bottom quadrant (m>0) : exchange m and m' => parity factor
					}
				}
				rl[-mp] += rmmp;
				rl[mp] += rmp;
				double* t = mx0_;		mx0_ = mx1_;	mx1_ = t;	// cycle the buffers for lines.
			}

			// step 4 + 5 merged:	recursively compute and apply d(m',m),  m'=2..l; m=m'..l  AND  m'=-2..-l; m=-m'..l
			double* mx0 = m0;
			double* mx1 = m0 + lw;
			for (int mp=2; mp<=l; mp++) {
				const double cmp2 = clm[mp-2];
				double clm_1 = clm[mp-1];
				const double c_1 = 1.0/clm_1;
				//v2d mx1_1 = *((v2d*)(mx1+mp-1));
				//v2d mx1__1 = *((v2d*)(mx1_+mp-1));
				v2d rmp = vdup(0.0);
				v2d rmmp = vdup(0.0);
				v2d zlmmp = zl[-mp];
				v2d zlmp = zl[mp];
				int m = mp;
				for (; m<=l-(VSIZE2-1); m+=VSIZE2) {
					rnd clm_1 = *((rndu*)(clm+m-1));
					rnd clm0 =  *((rndu*)(clm+m));

					rnd mx1_1 =  *((rndu*)(mx1+m-1));
					rnd mx11 =  *((rndu*)(mx1+m+1));

					rnd mx1__1 =  *((rndu*)(mx1_+m-1));
					rnd mx11_  =  *((rndu*)(mx1_+m+1));

					rnd mx00  = *((rndu*)(mx0+m));
					rnd mx00_ = *((rndu*)(mx0_+m));

					rnd H  = (vall(cmp2) * mx00   + clm_1 * mx1_1  - clm0 * mx11)  * vall(c_1);		// eq 106
					rnd H_ = (vall(cmp2) * mx00_  - clm_1 * mx1__1 + clm0 * mx11_) * vall(c_1);		// eq 108
					//mx1__1 = mx11_;		mx1_1 = mx11;
					*((rndu*)(mx0+m)) = H;			// d(m',m) -> overwrite d(m'-2,m)
					*((rndu*)(mx0_+m)) = H_;		// d(-m',m) -> overwrite d(-m'+2,m)
				}
		/*		v2d mx1_1 =  *((v2du*)(mx1+m-1));
				v2d mx1__1 =  *((v2du*)(mx1_+m-1));
				for (; m<l; m+=2) {
					v2d clm_1 = *((v2du*)(clm+m-1));
					v2d clm0 =  *((v2du*)(clm+m));

					v2d mx11 =  *((v2du*)(mx1+m+1));

					v2d mx11_  =  *((v2du*)(mx1_+m+1));

					v2d mx00  = *((v2du*)(mx0+m));
					v2d mx00_ = *((v2du*)(mx0_+m));

					v2d H = vdup(cmp2) * mx00   + clm_1 * mx1_1 - clm0 * mx11;		// eq 106
					v2d H_ = vdup(cmp2) * mx00_  - clm_1 * mx1__1 + clm0 * mx11_;	// eq 108
					mx1__1 = mx11_;		mx1_1 = mx11;
					H *= vdup(c_1);
					H_ *= vdup(c_1);
					*((v2du*)(mx0+m)) = H;
					*((v2du*)(mx0_+m)) = H_;
				}	*/
				for (; m<=l; m++) {
					double clm_1 = clm[m-1];
					double clm0 = clm[m];
					double mx1_1 = mx1[m-1];
					double mx1__1 = mx1_[m-1];
					double mx11 = mx1[m+1];
					double mx11_ = mx1_[m+1];
					double H0  = (cmp2 * mx0[m]  + clm_1 * mx1_1  - clm0 * mx11)  * c_1;	// eq 106
					double H0_ = (cmp2 * mx0_[m] - clm_1 * mx1__1 + clm0 * mx11_) * c_1;	// eq 108
					mx0[m] = H0;			// d(m',m) -> overwrite d(m'-2,m)
					mx0_[m] = H0_;			// d(-m',m) -> overwrite d(-m'+2,m)
				}

				m = mp;
				for (; m<l; m+=2) {
					double H0 = mx0[m];		double H1 = mx0[m+1];
					double H0_ = mx0_[m];	double H1_ = mx0_[m+1];
					rmp  += zl[m]  * vdup(H0)  + zl[m+1]  * vdup(H1);	//mx[mp*lw + m];		// right quadrant (m>0, mp>0)
					rmmp += zl[m]  * vdup(H0_) + zl[m+1]  * vdup(H1_);	//mx[mp*lw + m];		// right quadrant (m>0, mp<0)
					rmmp += zl[-m] * vdup(H0)  - zl[-m-1] * vdup(H1);	// left quadrant, change signs => parity factor
					rmp  += zl[-m] * vdup(H0_) - zl[-m-1] * vdup(H1_);	// left quadrant, change signs => parity factor
					rl[-m-1] += zlmmp * vdup(H1)  +   zlmp  * vdup(H1_);			// top quadrant (m<0) : exchange m and m' and change signs => no parity factor
					if (m>mp) {	// avoid duplicates
						rl[-m] += zlmmp * vdup(H0)  +  zlmp  * vdup(H0_);	// top quadrant (m<0) : exchange m and m' and change signs => no parity factor
						rl[m]  += zlmp  * vdup(H0)  +  zlmmp * vdup(H0_);	// bottom quadrant (m>0) : exchange m and m' => parity factor
					}
					rl[m+1] -= zlmp  * vdup(H1)  +  zlmmp * vdup(H1_);	// bottom quadrant (m>0) : exchange m and m' => parity factor
				}
				if (m==l) {
					double H0 = mx0[m];
					double H0_ = mx0_[m];
					rmp  += zl[m]  * vdup(H0);	//mx[mp*lw + m];		// right quadrant (m>0, mp>0)
					rmmp += zl[m]  * vdup(H0_);	//mx[mp*lw + m];		// right quadrant (m>0, mp<0)
					rmmp += zl[-m] * vdup(H0);	// left quadrant, change signs => parity factor
					rmp  += zl[-m] * vdup(H0_);	// left quadrant, change signs => parity factor
					if (m>mp) {	// avoid duplicates
						rl[-m] += zlmmp * vdup(H0)  +  zlmp  * vdup(H0_);	// top quadrant (m<0) : exchange m and m' and change signs => no parity factor
						rl[m]  += zlmp  * vdup(H0)  +  zlmmp * vdup(H0_);	// bottom quadrant (m>0) : exchange m and m' => parity factor
					}
				}
				rl[-mp] += rmmp;
				rl[mp] += rmp;
				double* t = mx0;		mx0 = mx1;		mx1 = t;	// cycle the buffers for lines.
				double* t_ = mx0_;		mx0_ = mx1_;	mx1_ = t_;	// cycle the buffers for lines.
			}

			if (flag_ag & 2) {		// apply post-rotation along Z-axis.
				Rlm[ll + 0] = ((cplx*)rl)[0];
				cplx eim_gamma = eig;
				for (int m=1; m<=mlim; m++) {
					Rlm[ll - m] = ((cplx*)rl)[-m] * conj(eim_gamma);
					Rlm[ll + m] = ((cplx*)rl)[m]  * eim_gamma;
					eim_gamma *= eig;
				}
			} else {
				memcpy(Rlm + ll-mlim, rl-mlim, sizeof(cplx)*(2*mlim+1));		// copy to destination (avoid false sharing when doing openmp).
			}
		}

		free(buf);
	}
}
