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


/** \example SHT_example.c
 \brief An example program that performs backward and forward Spherical Harmonic Transforms using SHTns.

  Compile using : \code make SHT_example \endcode
  \see \ref usage
**/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include <shtns.h>

/// a simple function that writes a vector to a file
void write_vect(char *fn, double *vec, int N);

/// a simple function that writes a matrix to a file.
void write_mx(char *fn, double *mx, int N1, int N2);


int main()
{
	shtns_cfg shtns;		// handle to a sht transform configuration
	long int NLM;
	complex double *Slm, *Tlm;	// spherical harmonics coefficients (l,m space): complex numbers.
	double *Sh, *Th;		// real space : theta,phi
	long int i,im,lm;
	double t;

	const int lmax = 5;		// maximum degree of spherical harmonics
	const int mmax = 3;		// maximum order of spherical harmonics
	const int mres = 1;		// periodicity in phi (1 for full-sphere, 2 for half the sphere, 3 for 1/3, etc...)
	const int nlat = 64;	// number of points in the latitude direction  (constraint: nlat >= lmax+1)
	const int nphi = 16;	// number of points in the longitude direction (constraint: nphi >= 2*mmax+1)

	shtns_verbose(1);			// displays informations during initialization.
	shtns_use_threads(0);		// enable multi-threaded transforms (if supported).
	shtns = shtns_init( sht_gauss, lmax, mmax, mres, nlat, nphi );
//	shtns = shtns_create(lmax, mmax, mres, sht_orthonormal | SHT_REAL_NORM);
//	shtns_set_grid(shtns, sht_gauss, 0.0, nlat, nphi);
	NLM = shtns->nlm;

// Memory allocation : the use of shtns_malloc is recommended for proper alignement, 
// or to use 'pinned' memory for faster GPU transfers (if CUDA is used).
// Use shtns_free() to free the memory when no more needed.

// allocate spatial fields.
	Sh = (double *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(double));
	Th = (double *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(double));

// allocate SH representations.
	Slm = (complex double *) shtns_malloc( NLM * sizeof(complex double));
	Tlm = (complex double *) shtns_malloc( NLM * sizeof(complex double));

// SH_to_spat
	LM_LOOP(shtns,  Slm[lm]=0.0;  Tlm[lm] = 0.0; )		/* this is the same as :
						for (lm=0; lm<shtns->nlm; lm++) {
							Slm[lm] = 0.0;	Tlm[lm] = 0.0;
						} */

//	Slm[LM(shtns, 1,1)] = sh11_st(shtns);				// access to SH coefficient
	Slm[LM(shtns, 2,0)] = 1.0;
// 	Slm[LiM(shtns, 1,0)] = sh10_ct(shtns);
//	Slm[LiM(shtns, 0,0)] = 0.5*sh00_1(shtns);
	SH_to_spat(shtns, Slm,Sh);
	SHtor_to_spat(shtns, Slm,Th,Sh);
	write_vect("ylm",(double *) Slm,NLM*2);
	write_mx("spat",Sh,nphi,nlat);

// compute value of SH expansion at a given physical point.
 	double t2;
	SHqst_to_point(shtns, Tlm, Tlm, Slm, shtns->ct[nlat/3], 2.*M_PI/(mres*nphi),&t2,&t2,&t);
	printf("check if SH_to_point coincides with SH_to_spat : %f = %f\n",t,Sh[nlat/3]);
	printf("ct*st = %f\n",shtns->ct[nlat/3]*shtns->st[nlat/3]);

// check non-linear behaviour
	for (im=0;im<nphi;im++) {
		for (i=0;i<nlat;i++) {
			Sh[im*nlat+i] *= Sh[im*nlat+i];
		}
	}
	spat_to_SH(shtns, Sh, Tlm);		//  /!\ WARNING! this destroys the spatial data in Sh
	write_vect("ylm_nl",(double *) Tlm, NLM*2);

// vector transform
	LM_LOOP(shtns,  Slm[lm]=0.0;  Tlm[lm] = 0.0; )
	SHsphtor_to_spat(shtns, Slm,Tlm, Sh,Th);		// vector transform
	write_mx("spatt",Sh,nphi,nlat);
	write_mx("spatp",Th,nphi,nlat);

// spat_to_SH
	for (im=0;im<nphi;im++) {
		for (i=0;i<nlat;i++) {
			Sh[im*nlat+i] = shtns->ct[i];			// use cos(theta) array
		}
	}
	spat_to_SH(shtns, Sh,Slm);		//  /!\ WARNING! this destroys the spatial data in Sh
	write_vect("ylm_v",(double *) Slm,NLM*2);
}


void write_vect(char *fn, double *vec, int N)
{
	FILE *fp;
	int i;
	
	fp = fopen(fn,"w");
	for (i=0;i<N;i++) {
		fprintf(fp,"%.6g ",vec[i]);
	}
	fclose(fp);
}

void write_mx(char *fn, double *mx, int N1, int N2)
{
	FILE *fp;
	int i,j;
	
	fp = fopen(fn,"w");
	for (i=0;i<N1;i++) {
		for(j=0;j<N2;j++) {
			fprintf(fp,"%.6g ",mx[i*N2+j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

