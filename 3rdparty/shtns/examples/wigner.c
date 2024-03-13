
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "shtns.h"

#include <time.h>
#include <sys/time.h>

void SH_2real_to_cplx(shtns_cfg shtns, cplx* Rlm, cplx* Ilm, cplx* Zlm);
void SH_cplx_to_2real(shtns_cfg shtns, cplx* Zlm, cplx* Rlm, cplx* Ilm);

void write_vect(char *fn, double *vec, int N)
{
	FILE *fp;
	int i;
	
	fp = fopen(fn,"w");
	for (i=0;i<N;i++) {
		fprintf(fp,"%.16g ",vec[i]);
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
			fprintf(fp,"%.16g ",mx[i*N2+j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

/// for real-time performance measurements, returns time in mili-seconds.
double tdiff(struct timeval *start, struct timeval *end)
{
	double sec = ((long) end->tv_sec - (long) start->tv_sec);
	sec += 1.e-6*((long) end->tv_usec - (long) start->tv_usec);
	return sec * (1.e3);	// time in ms.
}

double SHnorm(shtns_cfg shtns, cplx* zlm)
{
	double x=0;
	int lmax = shtns->lmax;
	for (int j=0;j<(lmax+1)*(lmax+1); j++) {
		double r = creal(zlm[j]);
		double i = cimag(zlm[j]);
		x += r*r + i*i;
	}
	return sqrt(x);
}

double SHnorm_real(shtns_cfg shtns, cplx* zlm)
{
	double x=0;
	int lmax = shtns->lmax;
	for (int l=0; l<=lmax; l++) {
		double r = creal(zlm[l]);
		x += r*r;
	}
	x *= 0.5;
	for (int lm=lmax+1; lm<shtns->nlm; lm++) {
		double r = creal(zlm[lm]);
		double i = cimag(zlm[lm]);
		x += r*r + i*i;
	}
	return sqrt(x+x);
}

int main()
{
	struct timeval t1, t2, t3;
	const int lmax = 300;
	const double beta = M_PI/2;
	shtns_rot rot;
	shtns_cfg sht;

	printf("Rotation using Wigner-d matrices, at lmax=%d\n",lmax);
	srand( time(NULL) );	// initialize random numbers.

	gettimeofday(&t1, NULL);
	rot = shtns_rotation_create(lmax,lmax,sht_orthonormal);	// third argument is normalization
	gettimeofday(&t2, NULL);
	shtns_rotation_set_angles_ZYZ(rot, 0.0, beta, 0.0);
	gettimeofday(&t3, NULL);
	printf("rotation: creation time=%.3g ms   set_angle time=%.3g ms\n", tdiff(&t1,&t2), tdiff(&t2,&t3));

	if (lmax>=2) {
		const int l = 2;
		double* mx = (double*) malloc( sizeof(double) * (2*l+1)*(2*l+1) );
		memset(mx, 0, sizeof(double)*(2*l+1)*(2*l+1));
		shtns_rotation_wigner_d_matrix(rot, l, mx);			
		write_mx("Wigner-d_l2", mx, 2*l+1, 2*l+1);

		printf("# rotation matrix for l=%d:\n#  ", l);
		for (int j=0; j<2*l+1; j++) printf("%10d ", j-l);
		for (int i=0; i<2*l+1; i++) {
			printf("\n%3d ", i-l);
			for (int j=0; j<2*l+1; j++) {
				double x = mx[i*(2*l+1) + j];
				printf("%10.5g ", (fabs(x) < 1e-14) ? 0.0 : x);
			}
		}
		printf("\n");
		free(mx);
	}

	if (lmax >=9) {
		const int l = 9;
		double* mx = (double*) malloc( sizeof(double) * (2*l+1)*(2*l+1) );
		memset(mx, 0, sizeof(double)*(2*l+1)*(2*l+1));
		shtns_rotation_wigner_d_matrix(rot, l, mx);			
		write_mx("Wigner-d_l9", mx, 2*l+1, 2*l+1);
		free(mx);
	}

	
	//// TEST ////
	
	sht = shtns_create(lmax,lmax,1,sht_orthonormal);

	long nlm = nlm_calc(lmax, lmax, 1);
	cplx* Qlm = malloc(sizeof(cplx) * nlm);
	memset(Qlm, 0, sizeof(cplx)*nlm);
	cplx* Slm = malloc(sizeof(cplx) * nlm);
	memset(Slm, 0, sizeof(cplx)*nlm);
	cplx* Tlm = malloc(sizeof(cplx) * nlm);
	memset(Tlm, 0, sizeof(cplx)*nlm);
	cplx* Zlm = malloc(sizeof(cplx) * (lmax+1)*(lmax+1));
	memset(Zlm, 0, sizeof(cplx)*(lmax+1)*(lmax+1));
	cplx* Rlm = malloc(sizeof(cplx) * (lmax+1)*(lmax+1));
	memset(Rlm, 0, sizeof(cplx)*(lmax+1)*(lmax+1));

	Qlm[0] = 1;
//	Qlm[LiM(sht,5,2)] = 1;
//	Qlm[LiM(sht,7,3)] = I;
//	write_vect("qlm",(double *)Qlm,nlm*2);

// test case...
	printf("generating random test case...\n");
	double t = 1.0 / (RAND_MAX/2);
	for (int i=0;i<nlm;i++) {
		Qlm[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		//Slm[i] = t*((double) (rand() - RAND_MAX/2)) + I*t*((double) (rand() - RAND_MAX/2));
		Slm[i] = 0.0;
	}
	for (int i=0;i<=lmax;i++)	Qlm[i] = creal(Qlm[i]);		// m=0 is REAL
	for (int i=0;i<=lmax;i++)	Slm[i] = 0.0;	//creal(Slm[i]);		// m=0 is REAL

	{
		printf("*** real field ***\n");
		printf("norm=%g\n", SHnorm_real(sht,Qlm));
		gettimeofday(&t1, NULL);
		shtns_rotation_apply_real(rot,Qlm,Zlm);
		gettimeofday(&t2, NULL);
		printf("norm=%g    time=%.3g ms\n", SHnorm_real(sht,Zlm), tdiff(&t1,&t2));
		shtns_rotation_apply_real(rot,Zlm,Rlm);		printf("norm=%g\n", SHnorm_real(sht,Rlm));
		shtns_rotation_apply_real(rot,Rlm,Zlm);		printf("norm=%g\n", SHnorm_real(sht,Zlm));
		gettimeofday(&t1, NULL);
		shtns_rotation_apply_real(rot,Zlm,Rlm);
		gettimeofday(&t2, NULL);
		printf("norm=%g    time=%.3g ms\n", SHnorm_real(sht,Rlm), tdiff(&t1,&t2));

		// evaluate error:
		double emax = 0.0;		int imax = 0;
		double esum = 0.0;
		for (int i=0;i<nlm;i++) {
			double qr = creal(Qlm[i]);		double qi = cimag(Qlm[i]);
			double sr = creal(Rlm[i]);		double si = cimag(Rlm[i]);
			double e = (sr-qr)*(sr-qr) + (si-qi)*(si-qi);
			if (e > emax) { emax = e;  imax=i; }
			esum += e;
		}
		printf("after 4 90° rotations along Y:    rms error = %.3g,   max error = %.3g (lm=%d)\n", sqrt(esum/nlm), sqrt(emax), imax);
	}
	
	{
		printf("\n*** complex field ***\n");
		SH_2real_to_cplx(sht,Qlm,Slm,Zlm);
		printf("norm=%g\n", SHnorm(sht,Zlm));

		gettimeofday(&t1, NULL);
		shtns_rotation_apply_cplx(rot,Zlm,Rlm);
		gettimeofday(&t2, NULL);
		printf("norm=%g    time=%.3g ms\n", SHnorm(sht,Rlm), tdiff(&t1,&t2));

		shtns_rotation_apply_cplx(rot,Rlm,Zlm);		printf("norm=%g\n", SHnorm(sht,Zlm));
		shtns_rotation_apply_cplx(rot,Zlm,Rlm);		printf("norm=%g\n", SHnorm(sht,Rlm));
		gettimeofday(&t1, NULL);
		shtns_rotation_apply_cplx(rot,Rlm,Zlm);
		gettimeofday(&t2, NULL);
		printf("norm=%g    time=%.3g ms\n", SHnorm(sht,Zlm), tdiff(&t1,&t2));
		SH_cplx_to_2real(sht,Zlm, Tlm,Slm);

		// evaluate error:
		double emax = 0.0;		int imax = 0;
		double esum = 0.0;
		for (int i=0;i<nlm;i++) {
			double qr = creal(Qlm[i]);		double qi = cimag(Qlm[i]);
			double sr = creal(Tlm[i]);		double si = cimag(Tlm[i]);
			double e = (sr-qr)*(sr-qr) + (si-qi)*(si-qi);
			if (e > emax) { emax = e;  imax=i; }
			esum += e;
		}
		printf("after 4 90° rotations along Y:    rms error = %.3g,   max error = %.3g (lm=%d)\n", sqrt(esum/nlm), sqrt(emax), imax);
	}


//	write_vect("rlm",(double *)Tlm,nlm*2);
}
