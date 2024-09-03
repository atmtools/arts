#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _GCC_VEC_
	// SIMD vectorization by default, turn off by compiling with -D_GCC_VEC=0
	#define _GCC_VEC_ 1
#endif
#include "shtns_simd.h"

#include "shtns.h"

const double ftol = 1e-14;	// tolerance for approximate tests.
long global_count = 0;		// count the number of tests
long global_fail = 0;		// count the number of failures

#define COLOR_OK  "\033[92m"
#define COLOR_WRN "\033[93m"
#define COLOR_ERR "\033[91m"
#define COLOR_END "\033[0m"

const char* res_string[2] = {COLOR_ERR "FAIL" COLOR_END, 	COLOR_OK "OK" COLOR_END};

int test_condition(int line, int test, int display) {
	global_count++;
	global_fail += (!test);
	if ((!test) || (display))
		printf("line #%d : %s\n", line, res_string[!!test]);
}

int test_approx(int line, double a, double b, int display) {
	global_count++;
	double err = fabs(a-b);
	if (err > ftol) {
		global_fail++;
		printf("line #%d : %s (|%.3g - %.3g| = %.3g)\n" , line, res_string[0], a,b,err);
	} else if (display)
		printf("line #%d : %s (%.3g)\n" , line, res_string[1], err);
}

#define PRINT_TEST(test)      test_condition(__LINE__, (test), 1)
#define SILENT_TEST(test)     test_condition(__LINE__, (test), 0)
#define PRINT_TEST_APPROX(a, b)	 test_approx(__LINE__, a, b, 1)
#define SILENT_TEST_APPROX(a, b) test_approx(__LINE__, a, b, 0)


// *** HERE START THE TESTS DEFINITIONS ****

void test_simd()
{
	const double a[2] = {4.0, -0.3};

	printf(COLOR_WRN "** SIMD TESTS **" COLOR_END " (VSIZE=%d)\n", VSIZE2);
	v2d va = *(v2d*)a;

	SILENT_TEST( vcplx_real(va) == a[0] );
	SILENT_TEST( vcplx_imag(va) == a[1] );

	va = vread2(a,0);	// check that vread2 is the same as direct vector load
	SILENT_TEST( vcplx_real(va) == a[0] );
	SILENT_TEST( vcplx_imag(va) == a[1] );

	#ifdef _GCC_VEC_
	SILENT_TEST( vcplx_real(vxchg(va)) == a[1] );
	SILENT_TEST( vcplx_imag(vxchg(va)) == a[0] );
	#endif
	SILENT_TEST( vcplx_real(IxKxZ(3,va)) == -3*a[1] );
	SILENT_TEST( vcplx_imag(IxKxZ(3,va)) == 3*a[0] );

	rnd v = vall(a[1]);
	#if VSIZE2 > 1
	for (int k=0; k<VSIZE2; k++)
		SILENT_TEST( v[k] == a[1] );
	#else
		SILENT_TEST( v == a[1] );
	#endif

	#if VSIZE2 >= 4
		const double b[4] = {1.0, 2.0, 3.0, 4.0};
		v4d vx = vread4(b, 0);
		v4d vx_r = vreverse4(vx);
		SILENT_TEST( vlo(vx_r) == 4.0 );
		SILENT_TEST( vx_r[1] == 3.0 );
		SILENT_TEST( vx_r[2] == 2.0 );
		SILENT_TEST( vx_r[3] == 1.0 );
	#endif
}

double spat_func_test(double theta, double phi) {	
	return 1.0 + 0.01*cos(theta)  + 0.1*(3.*cos(theta)*cos(theta) - 1.0)	// Y00, Y10, Y20
	+ (cos(phi) + 0.3*sin(phi)) * sin(theta)	// Y11
	+ (cos(2.*phi) + 0.1*sin(2.*phi)) * sin(theta)*sin(theta) * (7.0* cos(theta)*cos(theta) - 1.0) * 3./8.; 	// Y42
}

void test_analys()
{
	printf(COLOR_WRN "** ANALYSE TESTS **" COLOR_END "\n");

	const int nlat = 64;
	const int nphi = 8;
	const int lmax = 7;
	shtns_cfg sht = shtns_init(sht_orthonormal | sht_quick_init, lmax, 3, 1, nlat, nphi);
	
	double d[nphi][nlat];
	complex double q[sht->nlm];
	
	for (int it=0; it<nlat; it++) {
		double theta = acos( sht->ct[it] );
		for (int ip=0; ip<nphi; ip++) {
			double phi = (2.*M_PI*ip)/nphi;
			d[ip][it] = spat_func_test(theta,phi);
		}
	}
	
	spat_to_SH(sht, (double*) d, q);	// Note that this destroys the content of 'd'
	
	//for (int lm=0; lm < sht->nlm; lm++) 	printf("lm=%d [%g %g]\n", lm, creal(q[lm]), cimag(q[lm]));

	for (int i=0; i<=lmax; i++)
		SILENT_TEST( cimag(q[0]) == 0.0 );
	SILENT_TEST( LM(sht, 1, 1) == lmax+1 );

	PRINT_TEST_APPROX( creal(q[0]), sqrt(4.*M_PI) );
	PRINT_TEST_APPROX( creal(q[1]), 0.01*sqrt(4.*M_PI/3.) );
	PRINT_TEST_APPROX( creal(q[2]), 0.1*sqrt(16.*M_PI/5.) );

	PRINT_TEST_APPROX( creal(q[lmax+1]), -sqrt(2.*M_PI/3.) );
	PRINT_TEST_APPROX( cimag(q[lmax+1]), 0.3*sqrt(2.*M_PI/3.) );
	PRINT_TEST_APPROX( creal(q[LM(sht, 4, 2)]), 0.5*sqrt(2.*M_PI/5.) );
	PRINT_TEST_APPROX( cimag(q[LM(sht, 4, 2)]), -0.05*sqrt(2.*M_PI/5.) );
	
	for (int lm=3; lm<sht->nlm; lm++) {
		if ((lm != LM(sht,1,1)) && (lm != LM(sht,4,2))) {
			SILENT_TEST_APPROX( creal(q[lm]), 0);
			SILENT_TEST_APPROX( cimag(q[lm]), 0);
		}
	}
	
	SH_to_spat(sht, q, (double*) d);
	
	double err_max = 0.0;
	for (int it=0; it<nlat; it++) {
		double theta = acos( sht->ct[it] );
		for (int ip=0; ip<nphi; ip++) {
			double phi = (2.*M_PI*ip)/nphi;
			double err = fabs( spat_func_test(theta,phi) - d[ip][it] );
			if (err > err_max) err_max = err;
		}
	}
	PRINT_TEST_APPROX(err_max, 0);
}


int main()
{
	test_simd();
	test_analys();
	
	if (global_fail==0)	{
		printf("\n** %ld tests " COLOR_OK "ALL OK" COLOR_END " **\n", global_count);
	} else {
		printf("\n** "COLOR_ERR "ERROR" COLOR_END ": %ld tests out of %ld " COLOR_ERR "FAILED" COLOR_END "** \n", global_fail, global_count);
	}
	return global_fail;
}
