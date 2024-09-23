#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fftw3.h"
#include "shtns.h"
#include "shtns_cuda.h"

template <typename T1, typename T2> T1 maxrelerr(T1 *expected, T2 *actual, unsigned int n) {
    T1 mre = 0.0f;
    unsigned int imax = 0;
    for (unsigned int i=0; i<n; i++)
    {
        T1 re = abs(expected[i] - actual[i]) / abs(expected[i]);
        if (i < 10)
            printf("re[%d] %g <- (%g, %g)\n", i, re, expected[i], actual[i]);
        if (re > mre) { mre = re; imax=i; }
    }
    printf("max re[%d] %g <- (%g, %g)\n", imax, mre, expected[imax], actual[imax]);
    return mre;
}

double test_double_host(int nlon, int nlat, int lmax, double *init_x) {
    int nspat = nlon * nlat;
    double *x1, *x2;
    cplx *q;
    shtns_cfg sht = shtns_create(lmax, lmax, 1, sht_orthonormal);
    shtns_set_grid(sht, ((enum shtns_type)(sht_quick_init)), 0, nlat, nlon);
    shtns_print_cfg(sht);
    x1 = (double*) fftw_malloc(nspat * sizeof(double));
    x2 = (double*) fftw_malloc(nspat * sizeof(double));
    q = (cplx*) fftw_malloc(sht->nlm * sizeof(cplx));
    memcpy(x1, init_x, nspat * sizeof(double));
//    memset(x2, 0, nspat * sizeof(double));
    spat_to_SH(sht, x1, q);
    SH_to_spat(sht, q, x2);
	memcpy(init_x, x2, nspat * sizeof(double));
    spat_to_SH(sht, x2, q);
    SH_to_spat(sht, q, x1);
    //return maxrelerr(x2, x1, nspat);
    return maxrelerr(init_x, x1, nspat);
}

double test_double_gpu(int nlon, int nlat, int lmax, double *init_x) {
    int nspat = nlon * nlat;
    double *x1, *x2, *hx1, *hx2;
    cplx *q;
    shtns_cfg sht = shtns_create(lmax, lmax, 1, sht_orthonormal);
    shtns_set_grid(sht, ((enum shtns_type)(sht_quick_init|SHT_ALLOW_GPU)), 0, nlat, nlon);
    shtns_print_cfg(sht);
    cudaMalloc(&x1, nspat * sizeof(double));
    cudaMalloc(&x2, nspat * sizeof(double));
    cudaMalloc(&q, sht->nlm * sizeof(cplx));
    cudaMemcpy(x1, init_x, nspat * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(x2, 0, nspat * sizeof(double));

    hx1 = (double*) malloc(nspat * sizeof(double));
    hx2 = (double*) malloc(nspat * sizeof(double));

    cu_spat_to_SH(sht, x1, q, lmax);
    cu_SH_to_spat(sht, q, x2, lmax);
    	cudaMemcpy(hx2, x2, nspat * sizeof(double), cudaMemcpyDeviceToHost);
    cu_spat_to_SH(sht, x2, q, lmax);
    cu_SH_to_spat(sht, q, x1, lmax);

    cudaMemcpy(hx1, x1, nspat * sizeof(double), cudaMemcpyDeviceToHost);
    maxrelerr(hx2, hx1, nspat);               // error for back & forth within GPU
    return maxrelerr(init_x, hx1, nspat);     // error comparing GPU with CPU (assuming init_x is the CPU data)
}

//double test_float() {
//    cu_spat_to_SH_float(sht, c_x64, c_Qlm, lmax);
//    cu_SH_to_spat_float(sht, c_Qlm, c_x64, lmax);
//}

int main()
{
    int nlon=520, nlat=256, lmax=150;
    double *init_x = (double*) malloc(nlat*nlon*sizeof(double));

//    shtns_verbose(1);
//    shtns_print_version();
    shtns_use_gpu(0);

    for (int i=0; i<(nlat*nlon); i++)
        init_x[i] = (1.0 * rand()) / RAND_MAX - 0.5;

    printf("TEST HOST\n");
    test_double_host(nlon, nlat, lmax, init_x);
    printf("TEST GPU\n");
    test_double_gpu(nlon, nlat, lmax, init_x);
//    double float_max_relerr = test_float(nlon, nlat, lmax, init_x);

    return 0;
}
