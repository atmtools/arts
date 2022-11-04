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

// Various CUDA kernels for SHTns

// adjustment for cuda
#undef SHT_L_RESCALE_FLY
#undef SHT_ACCURACY
#define SHT_L_RESCALE_FLY 1800
#define SHT_L_RESCALE_FLY_FLOAT 128
#define SHT_ACCURACY 1.0e-40
#define SHT_ACCURACY_FLOAT 1.0e-15
#define SHT_SCALE_FACTOR_FLOAT 72057594037927936.0


#if (__CUDACC_VER_MAJOR__ < 8) || ( defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600 )
__device__ double atomicAdd(double* address, double val)
{
	unsigned long long int* address_as_ull =
							 (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do {
		assumed = old;
	old = atomicCAS(address_as_ull, assumed,
						__double_as_longlong(val +
							   __longlong_as_double(assumed)));
	} while (assumed != old);
	return __longlong_as_double(old);
}
#endif

// define our own suffle macros, to accomodate cuda<9 and cuda>=9
#if __CUDACC_VER_MAJOR__ < 9
	#define shfl_xor(...) __shfl_xor(__VA_ARGS__)
	#define shfl_down(...) __shfl_down(__VA_ARGS__)
	#define shfl(...) __shfl(__VA_ARGS__)
	#define _any(p) __any(p)
	#define _all(p) __all(p)
#else
	#define shfl_xor(...) __shfl_xor_sync(0xFFFFFFFF, __VA_ARGS__)
	#define shfl_down(...) __shfl_down_sync(0xFFFFFFFF, __VA_ARGS__)
	#define shfl(...) __shfl_sync(0xFFFFFFFF, __VA_ARGS__)
	#define _any(p) __any_sync(0xFFFFFFFF, p)
	#define _all(p) __all_sync(0xFFFFFFFF, p)
#endif

/*
__device__ __forceinline__ int getLaneId() {
  int laneId;
  asm("mov.s32 %0, %laneid;" : "=r"(laneId) );
  return laneId;
}

__device__ __forceinline__ void namedBarrierWait(int name, int numThreads) {
  asm volatile("bar.sync %0, %1;" : : "r"(name), "r"(numThreads) : "memory");
}

__device__ __forceinline__ void namedBarrierArrived(int name, int numThreads) {
  asm volatile("bar.arrive %0, %1;" : : "r"(name), "r"(numThreads) : "memory");
}
*/

__global__ void copy_kernel(const double *in, float *out, const int n) {
	for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x) 
		out[i] = in[i];
}

static void copy_convert(shtns_cfg shtns) {
	const int BLOCKSIZE = 256;		// good value
	int n = 2*shtns->nlm;
	copy_kernel<<<(n+BLOCKSIZE-1)/BLOCKSIZE, BLOCKSIZE, 0, shtns->comp_stream>>>(shtns->d_alm, shtns->d_alm_f, n);
	n = 2*shtns->nlat_2;
	copy_kernel<<<(n+BLOCKSIZE-1)/BLOCKSIZE, BLOCKSIZE, 0, shtns->comp_stream>>>(shtns->d_ct, shtns->d_ct_f, n);
}


/// dim0, dim1 : size in complex numbers !
/// BLOCK_DIM_Y must be between 1 and 16
template<int BLOCK_DIM_Y, typename real> __global__ void
transpose_cplx_kernel(const real* in, real* out, const int dim0, const int dim1)
{
	const int TILE_DIM = WARPSZE/2;		// 16 complex numbers per warp, read as 32 reals.
	__shared__ real shrdMem[TILE_DIM][TILE_DIM+1][2];		// avoid shared mem conflicts

	const int lx = threadIdx.x >> 1;
	const int ly = threadIdx.y;
	const int ri = threadIdx.x & 1;		// real/imag index

	const int bx = TILE_DIM * blockIdx.x;
	const int by = TILE_DIM * blockIdx.y;

	int gx = lx + bx;
	int gy = ly + by;
	#pragma unroll
	for (int repeat = 0; repeat < TILE_DIM; repeat += BLOCK_DIM_Y) {
		int gy_ = gy+repeat;
		shrdMem[ly + repeat][lx][ri] = in[2*(gy_ * dim0 + gx) + ri];
	}

	// transpose tiles:
	gx = lx + by;
	gy = ly + bx;

	__syncthreads();
	// transpose within tile:
	#pragma unroll
	for (unsigned repeat = 0; repeat < TILE_DIM; repeat += BLOCK_DIM_Y) {
		int gy_ = gy+repeat;
		out[2*(gy_ * dim1 + gx) + ri] = shrdMem[lx][ly + repeat][ri];
	}
}

/// dim0, dim1 : size in complex numbers !
/// BLOCK_DIM_Y must be a power of 2 between 1 and 16
template<int BLOCK_DIM_Y, typename real> __global__ void
transpose_cplx_zero_kernel(const real* in, real* out, const int dim0, const int dim1, const int mmax)
{
	const int TILE_DIM = WARPSZE/2;		// 16 complex numbers per warp, read as 32 reals.
	__shared__ real shrdMem[TILE_DIM][TILE_DIM+1][2];		// avoid shared mem conflicts

	const int ly = threadIdx.y;
	const int lx = threadIdx.x >> 1;
	const int ri = threadIdx.x & 1;		// real/imag index

	const int by = TILE_DIM * blockIdx.y;
	const int bx = TILE_DIM * blockIdx.x;

	int gy = ly + by;
	int gx = lx + bx;

	if ((gy+(TILE_DIM-BLOCK_DIM_Y) <= mmax) || (gy >= dim1 - mmax)) {		// SAFE, no zero to insert
		#pragma unroll
		for (int repeat = 0; repeat < TILE_DIM; repeat += BLOCK_DIM_Y) {
			int gy_ = gy+repeat;
			shrdMem[ly + repeat][lx][ri] = in[2*(gy_ * dim0 + gx) + ri];
		}
	} else {
		for (int repeat = 0; repeat < TILE_DIM; repeat += BLOCK_DIM_Y) {
			int gy_ = gy+repeat;
			if ((gy_ <= mmax) || (gy_ >= dim1 - mmax)) {
				shrdMem[ly + repeat][lx][ri] = in[2*(gy_ * dim0 + gx) + ri];
			} else {
				shrdMem[ly + repeat][lx][ri] = 0.0;
			}
		}
	}

	// transpose tiles:
	gy = ly + bx;
	gx = lx + by;

	__syncthreads();
	// transpose within tile:
	#pragma unroll
	for (unsigned repeat = 0; repeat < TILE_DIM; repeat += BLOCK_DIM_Y) {
		int gy_ = gy+repeat;
		out[2*(gy_ * dim1 + gx) + ri] = shrdMem[lx][ly + repeat][ri];
	}
}

/// dim0, dim1 : size in complex numbers !
/// BLOCK_DIM_Y must be a power of 2 between 1 and 16
template<int BLOCK_DIM_Y, typename real> __global__ void
transpose_cplx_skip_kernel(const real* in, real* out, const int dim0, const int dim1, const int mmax)
{
	const int TILE_DIM = WARPSZE/2;		// 16 complex per warp, read as 32 reals.
	__shared__ real shrdMem[TILE_DIM][TILE_DIM+1][2];		// avoid shared mem conflicts

	const int lx = threadIdx.x >> 1;
	const int ly = threadIdx.y;
	const int ri = threadIdx.x & 1;		// real/imag index

	const int bx = TILE_DIM * blockIdx.x;
	const int by = TILE_DIM * blockIdx.y;

	int gx = lx + bx;
	int gy = ly + by;

	if ((gx <= mmax) || (gx >= dim0 - mmax)) {		// read only data if m<=mmax
		#pragma unroll
		for (int repeat = 0; repeat < TILE_DIM; repeat += BLOCK_DIM_Y) {
			int gy_ = gy+repeat;
			shrdMem[ly + repeat][lx][ri] = in[2*(gy_ * dim0 + gx) + ri];
		}
	}

	// transpose tiles:
	gy = ly + bx;
	gx = lx + by;

	__syncthreads();
	// transpose within tile:
	if ((gy <= mmax) || (gy+(TILE_DIM-BLOCK_DIM_Y) >= dim0 - mmax)) {		// write all useful data (+a bit more), the rest is ignored anyway
		#pragma unroll
		for (unsigned repeat = 0; repeat < TILE_DIM; repeat += BLOCK_DIM_Y) {
			int gy_ = gy+repeat;
			out[2*(gy_ * dim1 + gx) + ri] = shrdMem[lx][ly + repeat][ri];
		}
	}
}

/// dim0, dim1 must be multiple of 16.
template <typename real> static void
transpose_cplx(cudaStream_t stream, const real* in, real* out, const int dim0, const int dim1)
{
	const int block_dim_y = 4;		// good performance with 4 (MUST be power of 2 between 1 and 16)
	dim3 blocks(dim0/16, dim1/16);
	dim3 threads(32, block_dim_y);
	transpose_cplx_kernel<block_dim_y,real> <<<blocks, threads, 0, stream>>>(in, out, dim0, dim1);
}

/// dim0, dim1 must be multiple of 16.
template <typename real> static void
transpose_cplx_zero(cudaStream_t stream, const real* in, real* out, const int dim0, const int dim1, const int mmax)
{
	const int block_dim_y = 4;		// good performance with 4 (MUST be power of 2 between 1 and 16)
	dim3 blocks(dim0/16, dim1/16);
	dim3 threads(32, block_dim_y);
	transpose_cplx_zero_kernel<block_dim_y,real> <<<blocks, threads, 0, stream>>>(in, out, dim0, dim1, mmax);
}

/// dim0, dim1 must be multiple of 16.
template <typename real> static void
transpose_cplx_skip(cudaStream_t stream, const real* in, real* out, const int dim0, const int dim1, const int mmax)
{
	const int block_dim_y = 4;		// good performance with 4 (MUST be power of 2 between 1 and 16)
	dim3 blocks(dim0/16, dim1/16);
	dim3 threads(32, block_dim_y);
	transpose_cplx_skip_kernel<block_dim_y,real> <<<blocks, threads, 0, stream>>>(in, out, dim0, dim1, mmax);
}



/// On KEPLER, This kernel is fastest with THREADS_PER_BLOCK=256 and NW=1
template<int BLOCKSIZE, int S, int NW, typename real=double> __global__ void
leg_m0_kernel(const real* __restrict__ al, const real* __restrict__ ct, const real* __restrict__ ql, real* q, const int llim, const int nlat_2)
{
	// im = 0
	const int it = blockDim.x * blockIdx.x + threadIdx.x;
	const int j = threadIdx.x;

	__shared__ real ak[BLOCKSIZE];		// size blockDim.x
	__shared__ real qk[BLOCKSIZE/2];	// size blockDim.x / 2

	ak[j] = al[j];
	if ((j <= llim)&&(j<blockDim.x/2)) qk[j] = ql[2*j];
	__syncthreads();

	int l = 0;
	int k = 0;	int kq = 0;
	real cost[NW];
	real y0[NW];    real y1[NW];
	real re[NW];    real ro[NW];

	for (int i=0; i<NW; i++) {
	cost[i] = (it+i<nlat_2) ? ct[it+i] : 0.0;
	y0[i] = ak[0];
	if (S==1) y0[i] *= rsqrt(1.0 - cost[i]*cost[i]);	// for vectors, divide by sin(theta)
	}
	for (int i=0; i<NW; i++) {
	re[i] = y0[i] * qk[0];
	y1[i] = y0[i] * ak[1] * cost[i];
	}
	for (int i=0; i<NW; i++) {
	ro[i] = y1[i] * qk[1];
	}
	al+=2;    l+=2;	k+=2;	kq+=2;
	while(l<llim) {
	if (k+6 >= blockDim.x) {
		__syncthreads();
		ak[j] = al[j];
		if ((j <= llim)&&(j<blockDim.x/2)) qk[j] = ql[2*(l+j)];
		k=0;	kq=0;
		__syncthreads();
	}
	for (int i=0; i<NW; i++)	y0[i]  = ak[k+1]*cost[i]*y1[i] + ak[k]*y0[i];
	for (int i=0; i<NW; i++)	re[i] += y0[i] * qk[kq];
	for (int i=0; i<NW; i++)	y1[i]  = ak[k+3]*cost[i]*y0[i] + ak[k+2]*y1[i];
	for (int i=0; i<NW; i++)	ro[i] += y1[i] * qk[kq+1];
	al+=4;	l+=2;	k+=4;	kq+=2;
	}
	if (l==llim) {
	for (int i=0; i<NW; i++)	y0[i]  = ak[k+1]*cost[i]*y1[i] + ak[k]*y0[i];
	for (int i=0; i<NW; i++)	re[i] += y0[i] * qk[kq];
	}

	for (int i=0; i<NW; i++) {
	if (it+i < nlat_2) {
		q[it+i] = re[i]+ro[i];
		q[nlat_2*2-1-(it+i)] = re[i]-ro[i];
	}
	}
/*
	if (it < nlat_2) {
		int l = 0;
		real cost = ct[it];
		real y0 = al[0];
		real re = y0 * ql[0];
		real y1 = y0 * al[1] * cost;
		real ro = y1 * ql[1];
		al+=2;    l+=2;
		while(l<llim) {
			y0  = al[1]*(cost*y1) + al[0]*y0;
			re += y0 * ql[l];
			y1  = al[3]*(cost*y0) + al[2]*y1;
			ro += y1 * ql[l+1];
			al+=4;	l+=2;
		}
		if (l==llim) {
			y0  = al[1]*cost*y1 + al[0]*y0;
			re += y0 * ql[l];
		}

		q[it] = re+ro;
		q[nlat_2*2-1-it] = re-ro;
	}
	*/
}

template<int S, int NFIELDS, typename real=double>
static void leg_m0(shtns_cfg shtns, const real *ql, real *q, const int llim, int spat_dist = 0)
{
	const int nlat_2 = shtns->nlat_2;
	real *d_alm = (sizeof(real) >= 8) ? (real*) shtns->d_alm : (real*) shtns->d_alm_f;
	real *d_ct = (sizeof(real) >= 8) ? (real*) shtns->d_ct : (real*) shtns->d_ct_f;
	cudaStream_t stream = shtns->comp_stream;

	const int BLOCKSIZE = 256;		// good value
	const int NW = 1;

	// Launch the Legendre CUDA Kernel
	const int threadsPerBlock = BLOCKSIZE;	// can be from 32 to 1024, we should try to measure the fastest !
	const int blocksPerGrid = (nlat_2 + BLOCKSIZE*NW - 1) / (BLOCKSIZE*NW);
	if (spat_dist == 0) spat_dist = shtns->spat_stride;
	for (int f=0; f<NFIELDS; f++) {
		leg_m0_kernel<BLOCKSIZE, S,1,real> <<<blocksPerGrid, threadsPerBlock, 0, stream>>>(d_alm, d_ct, ql + f*shtns->nlm_stride, q + f*spat_dist, llim, nlat_2);
	}
}


/*
__inline__ __device__
void warp_reduce_add_4(double& re, double& ro, double& ie, double& io) {
  for (int offset = warpSize/2; offset > 0; offset >>= 1) {
	re += shfl_down(re, offset);
	ro += shfl_down(ro, offset);
	ie += shfl_down(ie, offset);
	io += shfl_down(io, offset);
  }
}

__inline__ __device__
void warp_reduce_add_2(double& ev, double& od) {
  for (int offset = warpSize/2; offset > 0; offset >>= 1) {
	ev += shfl_down(ev, offset);
	od += shfl_down(od, offset);
  }
}

__inline__ __device__
void warp_reduce_add(double& ev) {
  for (int offset = warpSize/2; offset > 0; offset >>= 1) {
	ev += shfl_down(ev, offset);
  }
}
*/

template<int BLOCKSIZE, int LSPAN, int S, int NFIELDS, typename real=double> __global__ void
ileg_m0_kernel(const real* __restrict__ al, const real* __restrict__ ct, const real* __restrict__ q, real *ql, const int llim, const int nlat_2, const int lmax, const int q_dist=0, const int ql_dist=0)
{
	const int it = BLOCKSIZE * blockIdx.x + threadIdx.x;
	const int j = threadIdx.x;

	// re-assign each thread an l (transpose)
	const int ll = j / (BLOCKSIZE/LSPAN);

	__shared__ real ak[2*LSPAN+2];	// cache
	__shared__ real yl[LSPAN*BLOCKSIZE];		// yl is also used for even/odd computation. Ensure LSPAN >= 4.
	const int l_inc = BLOCKSIZE;
	const real cost = (it < nlat_2) ? ct[it] : 0.0;
	real y0, y1;

	if (LSPAN < 4) printf("ERROR: LSPAN<4\n");

	real my_reo[NFIELDS][LSPAN];			// in registers
	if (j < 2*LSPAN+2) ak[j] = al[j];

	#pragma unroll
	for (int f=0; f<NFIELDS; f++) {
		y0 = (it < nlat_2) ? q[it + f*q_dist] : 0.0;				// north
		y1 = (it < nlat_2) ? q[nlat_2*2-1 - it + f*q_dist] : 0.0;	// south

		if ((f>0) && (BLOCKSIZE > WARPSZE)) 	__syncthreads();
		yl[j] = y0+y1;					// even
		yl[BLOCKSIZE +j] = y0-y1;		// odd
		if (BLOCKSIZE > WARPSZE) 	__syncthreads();

		// transpose reo to my_reo
		#pragma unroll
		for (int i=0, k=0; i<BLOCKSIZE; i+= BLOCKSIZE/LSPAN, k++) {
			int it = j % (BLOCKSIZE/LSPAN) + i;
			my_reo[f][k] = yl[(ll&1)*BLOCKSIZE +it];
		}
	}

	int l = 0;
	y0 = (it < nlat_2) ? ct[it + nlat_2] : 0.0;		// weights are stored just after ct.
	if (S==1) y0 *= rsqrt(1.0 - cost*cost);
	y0 *= ak[0];
	y1 = y0 * ak[1] * cost;

	if (BLOCKSIZE > WARPSZE)	__syncthreads();
	
	yl[j] = y0;
	yl[l_inc +j] = y1;
	al+=2;
	while (l <= llim) {
		for (int k=0; k<LSPAN; k+=2) {		// compute a block of the matrix, write it in shared mem.
			yl[k*l_inc +j]     = y0;
			y0 = ak[2*k+3]*cost*y1 + ak[2*k+2]*y0;
			yl[(k+1)*l_inc +j] = y1;
			y1 = ak[2*k+5]*cost*y0 + ak[2*k+4]*y1;
			al += 4;
		}
		if(BLOCKSIZE > WARPSZE)	__syncthreads();

		real qll[NFIELDS];	// accumulator
		// now re-assign each thread an l (transpose)
		const int itl = ll*l_inc + j % (BLOCKSIZE/LSPAN);
		#pragma unroll
		for (int f=0; f<NFIELDS; f++) qll[f] = my_reo[f][0] * yl[itl];			// first element
		#pragma unroll
		for (int i=BLOCKSIZE/LSPAN, k=1; i<BLOCKSIZE; i+= BLOCKSIZE/LSPAN, k++) {		// accumulate
			#pragma unroll
			for (int f=0; f<NFIELDS; f++)	qll[f] += my_reo[f][k] * yl[itl+i];
		}

		if (BLOCKSIZE/LSPAN <= WARPSZE) {	// reduce_add within same l is in same warp too:
			if (WARPSZE % (BLOCKSIZE/LSPAN)) printf("ERROR\n");
			#pragma unroll
			for (int ofs = BLOCKSIZE/(LSPAN*2); ofs > 0; ofs>>=1) {
				#pragma unroll
				for (int f=0; f<NFIELDS; f++)	qll[f] += shfl_down(qll[f], ofs, BLOCKSIZE/LSPAN);
			}
			if ( ((j % (BLOCKSIZE/LSPAN)) == 0) && ((l+ll)<=llim) ) {	// write result
				if (nlat_2 <= BLOCKSIZE) {		// do we need atomic add or not ?
					#pragma unroll
					for (int f=0; f<NFIELDS; f++)	ql[2*(l+ll) + f*ql_dist] = qll[f];
				} else {
					#pragma unroll
					for (int f=0; f<NFIELDS; f++)	atomicAdd(ql+2*(l+ll) + f*ql_dist, qll[f]);		// VERY slow atomic add on Kepler.
				}
			}
		} else {	// only partial reduction possible, finish with atomicAdd():
			if ((BLOCKSIZE/LSPAN) % WARPSZE) printf("ERROR\n");
			#pragma unroll
			for (int ofs = WARPSZE/2; ofs > 0; ofs>>=1) {
				#pragma unroll
				for (int f=0; f<NFIELDS; f++)	qll[f] += shfl_down(qll[f], ofs, WARPSZE);
			}
			__syncthreads();
			const int nsum = (BLOCKSIZE/(LSPAN*WARPSZE));
			if ((j % WARPSZE) == 0) {
				for (int f=0; f<NFIELDS; f++)  yl[ll*nsum + ((j/WARPSZE) % nsum) + f*LSPAN*nsum] = qll[f];
			}
			__syncthreads();
			if ( ((j % (BLOCKSIZE/LSPAN)) == 0) && ((l+ll)<=llim) ) {	// write result
				for (int i=1; i<nsum; i++) {
					for (int f=0; f<NFIELDS; f++)	qll[f] += yl[ll*nsum + i + f*LSPAN*nsum];
				}
				if (nlat_2 <= BLOCKSIZE) {		// do we need atomic add or not ?
					#pragma unroll
					for (int f=0; f<NFIELDS; f++)	ql[2*(l+ll) + f*ql_dist] = qll[f];
				} else {
					#pragma unroll
					for (int f=0; f<NFIELDS; f++)	atomicAdd(ql+2*(l+ll) + f*ql_dist, qll[f]);		// VERY slow atomic add on Kepler.
				}
			}
		/*	if ( ((j % WARPSZE) == 0) && ((l+ll)<=llim) ) {	// write result
				#pragma unroll
				for (int f=0; f<NFIELDS; f++)	atomicAdd(ql+2*(l+ll) + f*ql_dist, qll[f]);		// VERY slow atomic add on Kepler.
			}*/
		}

		if (j<2*LSPAN) ak[j+2] = al[j];
		if (BLOCKSIZE > WARPSZE)	__syncthreads();
		l+=LSPAN;
	}
}

template<int S, int NFIELDS, typename real=double>
static void ileg_m0(shtns_cfg shtns, const real* q, real *ql, const int llim, int q_dist=0, int ql_dist=0)
{
	const int nlat_2 = shtns->nlat_2;
	real *d_alm = (sizeof(real) >= 8) ? (real*) shtns->d_alm : (real*) shtns->d_alm_f;
	real *d_ct = (sizeof(real) >= 8) ? (real*) shtns->d_ct : (real*) shtns->d_ct_f;
	cudaStream_t stream = shtns->comp_stream;

	const int BLOCKSIZE = 256/NFIELDS;
	const int LSPAN_ = 8/NFIELDS;
	const int NW = 1;

	const int threadsPerBlock = BLOCKSIZE;	// can be from 32 to 1024, we should try to measure the fastest !
	const int blocksPerGrid = (nlat_2 + BLOCKSIZE*NW - 1) / (BLOCKSIZE*NW);
	if (q_dist == 0) q_dist = shtns->spat_stride;
	if (ql_dist == 0) ql_dist = shtns->nlm_stride;
	ileg_m0_kernel<BLOCKSIZE, LSPAN_, S, NFIELDS, real><<<blocksPerGrid, threadsPerBlock, 0, stream>>>(d_alm, d_ct, q, ql, llim, nlat_2, q_dist, ql_dist);
}


/** \internal convert from vector SH to scalar SH
	Vlm =  st*d(Slm)/dtheta + I*m*Tlm
	Wlm = -st*d(Tlm)/dtheta + I*m*Slm
*/
template<int BLOCKSIZE, typename real=double> __global__ void
sphtor2scal_kernel(const double* __restrict__ mx, const real* __restrict__ slm, const real* __restrict__ tlm, real *vlm, real *wlm, const int llim, const int lmax, const int mres)
{
	// indices for overlapping blocks:
	const int ll = (blockDim.x-4) * blockIdx.x + threadIdx.x - 2;		// = 2*l + ((imag) ? 1 : 0)
	const int j = threadIdx.x;
	const int im = blockIdx.y;

	__shared__ real sl[BLOCKSIZE];
	__shared__ real tl[BLOCKSIZE];
	__shared__ real M[BLOCKSIZE];

	const int m = im*mres;
	//int ofs = im*(2*(lmax+1) -m + mres);
//    const int xchg = 1 - 2*(ll&1);	// +1 for real and -1 for imag
//	const int xchg = ll - (ll^1);	// -1 for real and +1 for imag
	const int ofs = im*(((lmax+1)<<1) -m + mres) + ll;

	if ( (ll >= 0) && (ll < 2*(llim+1-m)) ) {
		M[j] = mx[ofs];
		sl[j] = slm[ofs];
		tl[j] = tlm[ofs];
	} else {
		M[j] = 0.0;
		sl[j] = 0.0;
		tl[j] = 0.0;
	}
	const real mimag = im * mres * (ll - (ll^1));

	__syncthreads();

//    if ((j>=2) && (j<BLOCKSIZE-2) && (ll < 2*(llim+2-m))) {
	if ((j<BLOCKSIZE-4) && (ll < 2*(llim+1-m))) {
		real ml = M[2*(j>>1)+1];
		real mu = M[2*(j>>1)+2];
		real v = mimag*tl[(j+2)^1]  +  (ml*sl[j] + mu*sl[j+4]);
		real w = mimag*sl[(j+2)^1]  -  (ml*tl[j] + mu*tl[j+4]);
		vlm[ofs+2*im+2] = v;
		wlm[ofs+2*im+2] = w;
	}
}

/** \internal convert from 2 scalar SH to vector SH
	Slm = - (I*m*Wlm + MX*Vlm) / (l*(l+1))
	Tlm = - (I*m*Vlm - MX*Wlm) / (l*(l+1))
**/
template<int BLOCKSIZE, typename real=double> __global__ void
scal2sphtor_kernel(const double* __restrict__ mx, const real* __restrict__ vlm, const real* __restrict__ wlm, real *slm, real *tlm, const int llim, const int lmax, const int mres)
{
	// indices for overlapping blocks:
	const int ll = (blockDim.x-4) * blockIdx.x + threadIdx.x - 2;		// = 2*l + ((imag) ? 1 : 0)
	const int j = threadIdx.x;
	const int im = blockIdx.y;

	__shared__ real vl[BLOCKSIZE];
	__shared__ real wl[BLOCKSIZE];
	__shared__ real M[BLOCKSIZE];

	const int m = im * mres;
	//const int xchg = 1 - 2*(j&1);	// +1 for real and -1 for imag
	//const int xchg = (j^1) - j;		// +1 for real and -1 for imag
	int ofs = im*(2*(lmax+1) -m + mres)  + ll;

	if ( (ll >= 0) && (ll < 2*(llim+1-m)) ) {
		M[j] = mx[ofs];
	} else M[j] = 0.0;

	if ( (ll >= 0) && (ll < 2*(llim+2-m)) ) {
		vl[j] = vlm[ofs+2*im];
		wl[j] = wlm[ofs+2*im];
	} else {
		vl[j] = 0.0;
		wl[j] = 0.0;
	}

	int ell = (ll>>1) + m + 1;		// +1 because we shift below

	__syncthreads();

//    if ((j>=2) && (j<THREADS_PER_BLOCK-2) && (ll < 2*(llim+1-m))) {
	if (j<BLOCKSIZE-4) {
		if ((ell <= llim) && (ell>0)) {
			const real mimag = im * mres * ((j^1) -j);
			real ll_1 = 1.0 / (ell*(ell+1));
			real ml = M[2*(j>>1)+1];
			real mu = M[2*(j>>1)+2];
			real s = mimag*wl[(j+2)^1]  -  (ml*vl[j] + mu*vl[j+4]);
			real t = mimag*vl[(j+2)^1]  +  (ml*wl[j] + mu*wl[j+4]);
			slm[ofs+2] = s * ll_1;
			tlm[ofs+2] = t * ll_1;
		} else if (ell <= lmax) {	// fill with zeros up to lmax (and l=0 too).
			slm[ofs+2] = 0.0;
			tlm[ofs+2] = 0.0;
		}
	}
}

void sphtor2scal_gpu(shtns_cfg shtns, cplx* d_Slm, cplx* d_Tlm, cplx* d_Vlm, cplx* d_Wlm, int llim, int mmax)
{
	dim3 blocks((2*(shtns->lmax+2)+MAX_THREADS_PER_BLOCK-5)/(MAX_THREADS_PER_BLOCK-4), mmax+1);
	dim3 threads(MAX_THREADS_PER_BLOCK, 1);
	sphtor2scal_kernel<MAX_THREADS_PER_BLOCK> <<< blocks, threads,0, shtns->comp_stream >>>
		(shtns->d_mx_stdt, (double*) d_Slm, (double*) d_Tlm, (double*) d_Vlm, (double*) d_Wlm, llim, shtns->lmax, shtns->mres);
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess) { printf("sphtor2scal_gpu error : %s!\n", cudaGetErrorString(err));	return; }
}

void scal2sphtor_gpu(shtns_cfg shtns, cplx* d_Vlm, cplx* d_Wlm, cplx* d_Slm, cplx* d_Tlm, int llim)
{
	dim3 blocks((2*(shtns->lmax+2)+MAX_THREADS_PER_BLOCK-5)/(MAX_THREADS_PER_BLOCK-4), shtns->mmax+1);
	dim3 threads(MAX_THREADS_PER_BLOCK, 1);
	scal2sphtor_kernel<MAX_THREADS_PER_BLOCK> <<<blocks, threads, 0, shtns->comp_stream>>>
		(shtns->d_mx_van, (double*) d_Vlm, (double*) d_Wlm, (double*)d_Slm, (double*)d_Tlm, llim, shtns->lmax, shtns->mres);
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess) { printf("scal2sphtor_gpu error : %s!\n", cudaGetErrorString(err));	return; }
}



/// requirements : blockSize must be 1 in the y-direction and THREADS_PER_BLOCK in the x-direction.
/// llim MUST BE <= 1800
/// S can only be 0 (for scalar) or 1 (for spin 1 / vector)
template<int BLOCKSIZE, int S, int NFIELDS, int NW, typename real=double>
static __global__ void leg_m_lowllim_kernel(
	const real* __restrict__ al, const real* __restrict__ ct, const real* __restrict__ ql, real *q,
	const int llim, const int nlat_2, const int lmax, const int mres, const int nphi, const int ql_dist=0, const int q_dist=0)
{
	const int it = BLOCKSIZE*NW * blockIdx.x + threadIdx.x;
	const int im = blockIdx.y;
	const int j = threadIdx.x;
	const int m_inc = 2*nlat_2;
	const int k_inc = 1;

	__shared__ real ak[BLOCKSIZE];		// size blockDim.x
	__shared__ real qk[NFIELDS][BLOCKSIZE];	// size blockDim.x * NFIELDS

	real cost[NW];
	real y0[NW];
	real y1[NW];
	#pragma unroll
	for (int i=0; i<NW; i++) {
		const int iit = it+i*BLOCKSIZE;
		cost[i] = (iit < nlat_2) ? ct[iit] : 0.0;
	}

	if (im==0) {
		ak[j] = al[j+2];
		if (j<2*(llim+1)) {
			#pragma unroll
			for (int f=0; f<NFIELDS; f++) 	qk[f][j] = ql[j  + f*ql_dist];
		}
		real re[NFIELDS][NW], ro[NFIELDS][NW];
		#pragma unroll
		for (int f=0; f<NFIELDS; f++) {
			#pragma unroll
			for (int i=0; i<NW; i++) {
				re[f][i] = 0.0;
				ro[f][i] = 0.0;
			}
		}
		int l = 0;
		for (int i=0; i<NW; i++) y0[i] = al[0];
		if (S==1) for (int i=0; i<NW; i++) y0[i] *= rsqrt(1.0 - cost[i]*cost[i]);	// for vectors, divide by sin(theta)
		for (int i=0; i<NW; i++) y1[i] = y0[i] * al[1] * cost[i];
		al+=2;
		__syncthreads();
		while(l<=llim-BLOCKSIZE/2) {
			#pragma unroll
			for (int k=0; k<BLOCKSIZE; k+=4) {
				#pragma unroll
				for (int i=0; i<NW; i++) {
					#pragma unroll
					for (int f=0; f<NFIELDS; f++)	re[f][i] += y0[i] * qk[f][k];
				}
				#pragma unroll
				for (int i=0; i<NW; i++) 	y0[i] = ak[k+1]*cost[i]*y1[i] + ak[k]*y0[i];
				#pragma unroll
				for (int i=0; i<NW; i++) {
					#pragma unroll
					for (int f=0; f<NFIELDS; f++)	ro[f][i] += y1[i] * qk[f][k+2];
				}
				#pragma unroll
				for (int i=0; i<NW; i++)	y1[i] = ak[k+3]*cost[i]*y0[i] + ak[k+2]*y1[i];
			}
			al += BLOCKSIZE;
			l += BLOCKSIZE/2;
			__syncthreads();
			if (l+j/2 <= llim) {
				ak[j] = al[j];
				#pragma unroll
				for (int f=0; f<NFIELDS; f++)	qk[f][j] = ql[2*l+j + f*ql_dist];
			}
			__syncthreads();
		}
		int k=0;
		while(l<llim) {
			#pragma unroll
			for (int i=0; i<NW; i++) {
				#pragma unroll
				for (int f=0; f<NFIELDS; f++)	re[f][i] += y0[i] * qk[f][k];
			}
			#pragma unroll
			for (int i=0; i<NW; i++) 	y0[i]  = ak[k+1]*cost[i]*y1[i] + ak[k]*y0[i];
			#pragma unroll
			for (int i=0; i<NW; i++) {
				#pragma unroll
				for (int f=0; f<NFIELDS; f++)	ro[f][i] += y1[i] * qk[f][k+2];
			}
			#pragma unroll
			for (int i=0; i<NW; i++)	y1[i]  = ak[k+3]*cost[i]*y0[i] + ak[k+2]*y1[i];
			l+=2;	  k+=4;
		}
		if (l==llim) {
			#pragma unroll
			for (int i=0; i<NW; i++) {
				#pragma unroll
				for (int f=0; f<NFIELDS; f++)	re[f][i] += y0[i] * qk[f][k];
			}
		}
		#pragma unroll
		for (int i=0; i<NW; i++) {
			const int iit = it+i*BLOCKSIZE;
			if (iit < nlat_2) {
				// store mangled for complex fft
				const int iit = it+i*BLOCKSIZE;
				#pragma unroll
				for (int f=0; f<NFIELDS; f++) {
					q[iit*k_inc + f*q_dist] = re[f][i]+ro[f][i];
					q[(nlat_2*2-1-iit)*k_inc + f*q_dist] = re[f][i]-ro[f][i];
				}
			}
		}
	} else { 	// m>0
		real rer[NFIELDS][NW], ror[NFIELDS][NW], rei[NFIELDS][NW], roi[NFIELDS][NW];
		int m = im*mres;
		int l = (im*(2*(lmax+1)-(m+mres)))>>1;
		al += 2*(l+m);
		ql += 2*(l + S*im);	// allow vector transforms where llim = lmax+1

		#pragma unroll
		for (int i=0; i<NW; i++) 	y1[i] = sqrt(1.0 - cost[i]*cost[i]);		// y1 = sin(theta)
		ak[j] = al[j+2];
		#pragma unroll
		for (int f=0; f<NFIELDS; f++)	if (m+j/2 <= llim) qk[f][j] = ql[2*m+j + f*ql_dist];

		#pragma unroll
		for (int i=0; i<NW; i++) {
			#pragma unroll
			for (int f=0; f<NFIELDS; f++) {
				ror[f][i] = 0.0;		roi[f][i] = 0.0;
				rer[f][i] = 0.0;		rei[f][i] = 0.0;
			}
			y0[i] = 1.0;
		}
		l = m - S;
		do {		// sin(theta)^(m-S)
			if (l&1) {
				#pragma unroll
				for (int i=0; i<NW; i++) y0[i] *= y1[i];
			}
			#pragma unroll
			for (int i=0; i<NW; i++) y1[i] *= y1[i];
		} while(l >>= 1);

		#pragma unroll
		for (int i=0; i<NW; i++) y0[i] *= al[0];
		#pragma unroll
		for (int i=0; i<NW; i++) y1[i] = al[1]*y0[i]*cost[i];

		__syncthreads();
		l=m;		al+=2;
		while (l<=llim - BLOCKSIZE/2) {	// compute even and odd parts
			#pragma unroll
			for (int k = 0; k<BLOCKSIZE; k+=4) {
				#pragma unroll
				for (int f=0; f<NFIELDS; f++) {
					#pragma unroll
					for (int i=0; i<NW; i++) {
						rer[f][i] += y0[i] * qk[f][k];		// real
						rei[f][i] += y0[i] * qk[f][k+1];	// imag
					}
				}
				#pragma unroll
				for (int i=0; i<NW; i++) 	y0[i] = ak[k+1]*(cost[i]*y1[i]) + ak[k]*y0[i];
				#pragma unroll
				for (int f=0; f<NFIELDS; f++) {
					#pragma unroll
					for (int i=0; i<NW; i++) {
						ror[f][i] += y1[i] * qk[f][k+2];	// real
						roi[f][i] += y1[i] * qk[f][k+3];	// imag
					}
				}
				#pragma unroll
				for (int i=0; i<NW; i++) 	y1[i] = ak[k+3]*(cost[i]*y0[i]) + ak[k+2]*y1[i];
			}
			al += BLOCKSIZE;
			l += BLOCKSIZE/2;
			__syncthreads();
			if (l+j/2 <= llim) {
				ak[j] = al[j];
				#pragma unroll
				for (int f=0; f<NFIELDS; f++)	qk[f][j] = ql[2*l+j + f*ql_dist];
			}
			__syncthreads();
		}
		int k=0;
		while (l<llim) {	// compute even and odd parts
			#pragma unroll
			for (int f=0; f<NFIELDS; f++) {
				#pragma unroll
				for (int i=0; i<NW; i++) {
					rer[f][i] += y0[i] * qk[f][k];		// real
					rei[f][i] += y0[i] * qk[f][k+1];	// imag
				}
			}
			#pragma unroll
			for (int i=0; i<NW; i++) 	y0[i] = ak[k+1]*(cost[i]*y1[i]) + ak[k]*y0[i];
			#pragma unroll
			for (int f=0; f<NFIELDS; f++) {
				#pragma unroll
				for (int i=0; i<NW; i++) {
					ror[f][i] += y1[i] * qk[f][k+2];	// real
					roi[f][i] += y1[i] * qk[f][k+3];	// imag
				}
			}
			#pragma unroll
			for (int i=0; i<NW; i++) 	y1[i] = ak[k+3]*(cost[i]*y0[i]) + ak[k+2]*y1[i];
			l+=2;	k+=4;
		}
		if (l==llim) {
			#pragma unroll
			for (int f=0; f<NFIELDS; f++) {
				#pragma unroll
				for (int i=0; i<NW; i++) {
					rer[f][i] += y0[i] * qk[f][k];		// real
					rei[f][i] += y0[i] * qk[f][k+1];	// imag
				}
			}
		}

		/// store mangled for complex fft
		#pragma unroll
		for (int i=0; i<NW; i++) {
			#pragma unroll
			for (int f=0; f<NFIELDS; f++) {
				rei[f][i] = shfl_xor(rei[f][i], 1);
				roi[f][i] = shfl_xor(roi[f][i], 1);
			}
		}
		real nr[NFIELDS][NW];
		const real sgn = 1 - 2*(j&1);
		#pragma unroll
		for (int i=0; i<NW; i++) {
			const int iit = it+i*BLOCKSIZE;
			if (iit < nlat_2) {
				#pragma unroll
				for (int f=0; f<NFIELDS; f++) {
					nr[f][i] =  rer[f][i]+ror[f][i];
					rer[f][i] = rer[f][i]-ror[f][i];
					ror[f][i] = sgn*(rei[f][i]+roi[f][i]);
					rei[f][i] = sgn*(rei[f][i]-roi[f][i]);
				}
				#pragma unroll
				for (int f=0; f<NFIELDS; f++) {
					q[im*m_inc + iit*k_inc + f*q_dist]                     = nr[f][i]  - ror[f][i];
					q[(nphi-im)*m_inc + iit*k_inc + f*q_dist]              = nr[f][i]  + ror[f][i];
					q[im*m_inc + (nlat_2*2-1-iit)*k_inc + f*q_dist]        = rer[f][i] + rei[f][i];
					q[(nphi-im)*m_inc + (nlat_2*2-1-iit)*k_inc + f*q_dist] = rer[f][i] - rei[f][i];
				}
			}
		}
	}
}

template<int S, int NFIELDS, typename real=double>
static void leg_m_lowllim(shtns_cfg shtns, const real *ql, real *q, const int llim, const int mmax, int spat_dist=0)
{
	const int lmax = shtns->lmax;
	const int mres = shtns->mres;
	const int nlat_2 = shtns->nlat_2;
	const int nphi = shtns->nphi;
	real *d_alm = (sizeof(real) >= 8) ? (real*) shtns->d_alm : (real*) shtns->d_alm_f;
	real *d_ct = (sizeof(real) >= 8) ? (real*) shtns->d_ct : (real*) shtns->d_ct_f;
	cudaStream_t stream = shtns->comp_stream;

	const int BLOCKSIZE = 256;		// good value
	const int NW = 2;

	// Launch the Legendre CUDA Kernel
	const int threadsPerBlock = BLOCKSIZE;	// can be from 32 to 1024, we should try to measure the fastest !
	const int blocksPerGrid = (nlat_2 + BLOCKSIZE*NW - 1) / (BLOCKSIZE*NW);
	if (spat_dist == 0) spat_dist = shtns->spat_stride;
	dim3 blocks(blocksPerGrid, mmax+1);
	dim3 threads(threadsPerBlock, 1);
	leg_m_lowllim_kernel<BLOCKSIZE, S, NFIELDS, NW, real> <<<blocks, threads, 0, stream>>>(d_alm, d_ct, (real*) ql, (real*) q, llim, nlat_2, lmax,mres, nphi, shtns->nlm_stride, spat_dist);
}

/// requirements : blockSize must be 1 in the y-direction and THREADS_PER_BLOCK in the x-direction.
/// llim can be arbitrarily large (> 1800)
template<int BLOCKSIZE, int S, typename real=double> __global__ void
leg_m_highllim_kernel(const real *al, const real *ct, const real *ql, real *q, const int llim, const int nlat_2, const int lmax, const int mres, const int nphi)
{
	const int it = blockDim.x * blockIdx.x + threadIdx.x;
	const int im = blockIdx.y;
	const int j = threadIdx.x;
	const int m_inc = 2*nlat_2;
	const int k_inc = 1;
	const real accuracy = (sizeof(real) >= 8) ? SHT_ACCURACY : SHT_ACCURACY_FLOAT;
	const real scale_factor = (sizeof(real) >= 8) ? SHT_SCALE_FACTOR : SHT_SCALE_FACTOR_FLOAT;

	__shared__ real ak[BLOCKSIZE];	// cache
	__shared__ real qk[BLOCKSIZE];

	const real cost = (it < nlat_2) ? ct[it] : 0.0;

	if (im==0) {
	int l = 0;
	real y0 = al[0];
	if (S==1) y0 *= rsqrt(1.0 - cost*cost);
	real re = y0 * ql[0];
	real y1 = y0 * al[1] * cost;
	real ro = y1 * ql[2];
	al+=2;    l+=2;
	while(l<llim) {
		y0  = al[1]*(cost*y1) + al[0]*y0;
		re += y0 * ql[2*l];
		y1  = al[3]*(cost*y0) + al[2]*y1;
		ro += y1 * ql[2*l+2];
		al+=4;	l+=2;
	}
	if (l==llim) {
		y0  = al[1]*cost*y1 + al[0]*y0;
		re += y0 * ql[2*l];
	}
	if (it < nlat_2) {
		// store mangled for complex fft
		q[it*k_inc] = re+ro;
		q[(nlat_2*2-1-it)*k_inc] = re-ro;
	}
	} else { 	// m>0
	int m = im*mres;
	int l = (im*(2*(lmax+1)-(m+mres)))>>1;
	al += 2*(l+m);
	ql += 2*(l + S*im);
	real rer,ror, rei,roi, y0, y1;
	ror = 0.0;	roi = 0.0;
	rer = 0.0;	rei = 0.0;
	y1 = sqrt(1.0 - cost*cost);	// sin(theta)
	if (_any(m - llim*y1 <= max(80, llim>>7))) {		// polar optimization (see Reinecke 2013), avoiding warp divergence
		y0 = 1.0;	// y0
		l = m - S;
		int ny = 0;
		int nsint = 0;
		do {		// sin(theta)^(m-S)		(use rescaling to avoid underflow)
		if (l&1) {
			y0 *= y1;
			ny += nsint;
			if (y0 < (accuracy+1.0/scale_factor)) {		// possible warp divergence
			ny--;
			y0 *= scale_factor;
			}
		}
		y1 *= y1;
		nsint += nsint;
		if (y1 < 1.0/scale_factor) {		// possible warp divergence
			nsint--;
			y1 *= scale_factor;
		}
		} while(l >>= 1);
		y0 *= al[0];
		y1 = 0.0;
//	    y1 = al[1]*y0*cost;

		l=m;	int ka = WARPSZE;
		const int ofs = j & 0xFFE0;

		while ( _all(ny<0) && (l<llim) ) {
			if (ka+4 >= WARPSZE) {
				ak[j] = al[(j&31)];
				ka=0;
			}
			y1 = ak[ka+1+ofs]*cost*y0 + ak[ka+ofs]*y1;
			y0 = ak[ka+3+ofs]*cost*y1 + ak[ka+2+ofs]*y0;
			l+=2;	al+=4;	ka+=4;
			if (fabs(y1) > accuracy*scale_factor + 1.0)
			{	// rescale when value is significant
				++ny;
				y0 *= 1.0/scale_factor;
				y1 *= 1.0/scale_factor;
			}
		}

		ka = WARPSZE;
		while (l<llim) {
			if (ka+4 >= WARPSZE) {		// cache coefficients
				ak[j] = al[(j&31)];
				qk[j] = ql[2*l+(j&31)];
				ka = 0;
			}
			y1 = ak[ka+1+ofs]*cost*y0 + ak[ka+ofs]*y1;
			if (ny==0) {
				rer += y0 * qk[ka+ofs];	// real
				rei += y0 * qk[ka+1+ofs];	// imag
				ror += y1 * qk[ka+2+ofs];	// real
				roi += y1 * qk[ka+3+ofs];	// imag
			}
			else if (fabs(y0) > accuracy*scale_factor + 1.0)
			{	// rescale when value is significant
				++ny;
				y0 *= 1.0/scale_factor;
				y1 *= 1.0/scale_factor;
			}
			l+=2;	al+=4;
			y0 = ak[ka+3+ofs]*cost*y1 + ak[ka+2+ofs]*y0;
			ka+=4;
		}
		if ((l==llim) && (ny==0)) {
			rer += y0 * ql[2*l];
			rei += y0 * ql[2*l+1];
		}
	}

	/// store mangled for complex fft
	real nr = rer+ror;
	real sr = rer-ror;
	const real sgn = 1 - 2*(j&1);
	rei = shfl_xor(rei, 1);
	roi = shfl_xor(roi, 1);
	real nix = sgn*(rei+roi);
	real six = sgn*(rei-roi);
	if (it < nlat_2) {
		q[im*m_inc + it*k_inc]                     = nr - nix;
		q[(nphi-im)*m_inc + it*k_inc]              = nr + nix;
		q[im*m_inc + (nlat_2*2-1-it)*k_inc]        = sr + six;
		q[(nphi-im)*m_inc + (nlat_2*2-1-it)*k_inc] = sr - six;
	}
	}
}

template<int S, int NFIELDS, typename real=double>
static void leg_m_highllim(shtns_cfg shtns, const real *ql, real *q, const int llim, const int mmax, int spat_dist = 0)
{
	const int lmax = shtns->lmax;
	const int mres = shtns->mres;
	const int nlat_2 = shtns->nlat_2;
	const int nphi = shtns->nphi;
	real *d_alm = (sizeof(real) >= 8) ? (real*) shtns->d_alm : (real*) shtns->d_alm_f;
	real *d_ct = (sizeof(real) >= 8) ? (real*) shtns->d_ct : (real*) shtns->d_ct_f;
	cudaStream_t stream = shtns->comp_stream;

	const int BLOCKSIZE = 256;		// good value
	const int NW = 1;

	// Launch the Legendre CUDA Kernel
	const int threadsPerBlock = BLOCKSIZE;	// can be from 32 to 1024, we should try to measure the fastest !
	const int blocksPerGrid = (nlat_2 + BLOCKSIZE*NW - 1) / (BLOCKSIZE*NW);
	if (spat_dist == 0) spat_dist = shtns->spat_stride;
	dim3 blocks(blocksPerGrid, mmax+1);
	dim3 threads(threadsPerBlock, 1);
	for (int f=0; f<NFIELDS; f++) {
		leg_m_highllim_kernel<BLOCKSIZE, S,real> <<<blocks, threads, 0, stream>>>(d_alm, d_ct, ql + f*shtns->nlm_stride, q + f*spat_dist, llim, nlat_2, lmax,mres, nphi);
	}
}


template<int BLOCKSIZE, int LSPAN, int S, int NFIELDS, typename real=double> __global__ void
ileg_m_lowllim_kernel(const real* __restrict__ al, const real* __restrict__ ct, const real* __restrict__ q, real *ql, const int llim, const int nlat_2, const int lmax, const int mres, const int nphi, const real mpos_scale, const int q_dist=0, const int ql_dist=0)
{
	const int it = BLOCKSIZE * blockIdx.x + threadIdx.x;
	const int j = threadIdx.x;
	const int im = blockIdx.y;
	const int m_inc = 2*nlat_2;
//    const int k_inc = 1;

	// re-assign each thread an l (transpose)
	const int ll = j / (BLOCKSIZE/LSPAN);
	const int ri = j / (BLOCKSIZE/(2*LSPAN)) % 2;	// real (0) or imag (1)

	__shared__ real ak[2*LSPAN+2];	// cache
	__shared__ real yl[LSPAN*BLOCKSIZE];		// yl is also used for even/odd computation. Ensure LSPAN >= 4.
	const int l_inc = BLOCKSIZE;
	const real cost = (it < nlat_2) ? ct[it] : 0.0;
	real y0, y1;

	if (LSPAN < 4) printf("ERROR: LSPAN<4\n");

	if (im == 0) {
		real my_reo[NFIELDS][LSPAN];			// in registers
		if (j < 2*LSPAN+2) ak[j] = al[j];

		#pragma unroll
		for (int f=0; f<NFIELDS; f++) {
			y0 = (it < nlat_2) ? q[it + f*q_dist] : 0.0;				// north
			y1 = (it < nlat_2) ? q[nlat_2*2-1 - it + f*q_dist] : 0.0;	// south

			if ((f>0) && (BLOCKSIZE > WARPSZE)) 	__syncthreads();
			yl[j] = y0+y1;					// even
			yl[BLOCKSIZE +j] = y0-y1;		// odd
			if (BLOCKSIZE > WARPSZE) 	__syncthreads();

			// transpose reo to my_reo
			#pragma unroll
			for (int i=0, k=0; i<BLOCKSIZE; i+= BLOCKSIZE/LSPAN, k++) {
				int it = j % (BLOCKSIZE/LSPAN) + i;
				my_reo[f][k] = yl[(ll&1)*BLOCKSIZE +it];
			}
		}

		int l = 0;
		y0 = (it < nlat_2) ? ct[it + nlat_2] : 0.0;		// weights are stored just after ct.
		if (S==1) y0 *= rsqrt(1.0 - cost*cost);
		y0 *= ak[0];
		y1 = y0 * ak[1] * cost;

		if (BLOCKSIZE > WARPSZE)	__syncthreads();
		
		yl[j] = y0;
		yl[l_inc +j] = y1;
		al+=2;
		while (l <= llim) {
			for (int k=0; k<LSPAN; k+=2) {		// compute a block of the matrix, write it in shared mem.
				yl[k*l_inc +j]     = y0;
				y0 = ak[2*k+3]*cost*y1 + ak[2*k+2]*y0;
				yl[(k+1)*l_inc +j] = y1;
				y1 = ak[2*k+5]*cost*y0 + ak[2*k+4]*y1;
				al += 4;
			}
			if(BLOCKSIZE > WARPSZE)	__syncthreads();

			real qll[NFIELDS];	// accumulator
			// now re-assign each thread an l (transpose)
			const int itl = ll*l_inc + j % (BLOCKSIZE/LSPAN);
			#pragma unroll
			for (int f=0; f<NFIELDS; f++) qll[f] = my_reo[f][0] * yl[itl];			// first element
			#pragma unroll
			for (int i=BLOCKSIZE/LSPAN, k=1; i<BLOCKSIZE; i+= BLOCKSIZE/LSPAN, k++) {		// accumulate
				#pragma unroll
				for (int f=0; f<NFIELDS; f++)	qll[f] += my_reo[f][k] * yl[itl+i];
			}

			if (BLOCKSIZE/LSPAN <= WARPSZE) {	// reduce_add within same l is in same warp too:
				if (WARPSZE % (BLOCKSIZE/LSPAN)) printf("ERROR\n");
				#pragma unroll
				for (int ofs = BLOCKSIZE/(LSPAN*2); ofs > 0; ofs>>=1) {
					#pragma unroll
					for (int f=0; f<NFIELDS; f++)	qll[f] += shfl_down(qll[f], ofs, BLOCKSIZE/LSPAN);
				}
				if ( ((j % (BLOCKSIZE/LSPAN)) == 0) && ((l+ll)<=llim) ) {	// write result
					if (nlat_2 <= BLOCKSIZE) {		// do we need atomic add or not ?
						#pragma unroll
						for (int f=0; f<NFIELDS; f++)	ql[2*(l+ll) + f*ql_dist] = qll[f];
					} else {
						#pragma unroll
						for (int f=0; f<NFIELDS; f++)	atomicAdd(ql+2*(l+ll) + f*ql_dist, qll[f]);		// VERY slow atomic add on Kepler.
					}
				}
			} else {	// only partial reduction possible, finish with atomicAdd():
				if ((BLOCKSIZE/LSPAN) % WARPSZE) printf("ERROR\n");
				#pragma unroll
				for (int ofs = WARPSZE/2; ofs > 0; ofs>>=1) {
					#pragma unroll
					for (int f=0; f<NFIELDS; f++)	qll[f] += shfl_down(qll[f], ofs, WARPSZE);
				}
				__syncthreads();
				const int nsum = (BLOCKSIZE/(LSPAN*WARPSZE));
				if ((j % WARPSZE) == 0) {
					for (int f=0; f<NFIELDS; f++)  yl[ll*nsum + ((j/WARPSZE) % nsum) + f*LSPAN*nsum] = qll[f];
				}
				__syncthreads();
				if ( ((j % (BLOCKSIZE/LSPAN)) == 0) && ((l+ll)<=llim) ) {	// write result
					for (int i=1; i<nsum; i++) {
						for (int f=0; f<NFIELDS; f++)	qll[f] += yl[ll*nsum + i + f*LSPAN*nsum];
					}
					if (nlat_2 <= BLOCKSIZE) {		// do we need atomic add or not ?
						#pragma unroll
						for (int f=0; f<NFIELDS; f++)	ql[2*(l+ll) + f*ql_dist] = qll[f];
					} else {
						#pragma unroll
						for (int f=0; f<NFIELDS; f++)	atomicAdd(ql+2*(l+ll) + f*ql_dist, qll[f]);		// VERY slow atomic add on Kepler.
					}
				}
			/*	if ( ((j % WARPSZE) == 0) && ((l+ll)<=llim) ) {	// write result
					#pragma unroll
					for (int f=0; f<NFIELDS; f++)	atomicAdd(ql+2*(l+ll) + f*ql_dist, qll[f]);		// VERY slow atomic add on Kepler.
				}*/
			}

			if (j<2*LSPAN) ak[j+2] = al[j];
			if (BLOCKSIZE > WARPSZE)	__syncthreads();
			l+=LSPAN;
		}
	} else {	// im > 0
		real my_reo[NFIELDS][2*LSPAN];			// in registers
		int m = im*mres;
		int l = (im*(2*(lmax+1)-(m+mres)))>>1;
		al += 2*(l+m);
		ql += 2*(l + S*im);	// allow vector transforms where llim = lmax+1

		if (j < 2*LSPAN+2) ak[j] = al[j];
		const real sgn = 2*(j&1) - 1;	// -/+
		
		#pragma unroll
		for (int f=0; f<NFIELDS; f++) {
			y0         = (it < nlat_2) ? q[im*m_inc + it + f*q_dist] : 0.0;		// north imag (ani)
			real qer = (it < nlat_2) ? q[(nphi-im)*m_inc + it + f*q_dist] : 0.0;	// north real (an)
			y1         = (it < nlat_2) ? q[im*m_inc + nlat_2*2-1-it + f*q_dist] : 0.0;	// south imag (asi)
			real qor = (it < nlat_2) ? q[(nphi-im)*m_inc + nlat_2*2-1-it + f*q_dist] : 0.0;	// south real (as)
			real qei = y0-qer;		qer += y0;		// ani = -qei[lane+1],   bni = qei[lane-1]
			real qoi = y1-qor;		qor += y1;		// bsi = -qoi[lane-1],   asi = qoi[lane+1];
			y0 = shfl_xor(qei, 1);	// exchange between adjacent lanes.
			y1 = shfl_xor(qoi, 1);

			if ((f>0) && (BLOCKSIZE > WARPSZE)) 	__syncthreads();

			yl[j] 		       = qer + qor;	// rer
			yl[BLOCKSIZE +j]   = qer - qor;	// ror
			yl[2*BLOCKSIZE +j] = sgn*(y0 - y1);	// rei
			yl[3*BLOCKSIZE +j] = sgn*(y0 + y1);	// roi

			if (BLOCKSIZE > WARPSZE) 	__syncthreads();
			// transpose yl to my_reo
			#pragma unroll
			for (int i=0, k=0; i<BLOCKSIZE; i+= BLOCKSIZE/(2*LSPAN), k++) {
				int it = j % (BLOCKSIZE/(2*LSPAN)) + i;
				my_reo[f][k] = yl[((ll&1)+2*ri)*BLOCKSIZE +it];
			}
		}

		y1 = sqrt(1.0 - cost*cost);	// sin(theta)

		y0 = mpos_scale * ak[0];	// y0
		l = m - S;
		do {		// sin(theta)^(m-S)
		if (l&1) y0 *= y1;
		y1 *= y1;
		} while(l >>= 1);
		if (it < nlat_2)     y0 *= ct[it + nlat_2];		// include quadrature weights.
		y1 = ak[1]*y0*cost;

		l=m;		al+=2;
		while (l <= llim) {
			if (BLOCKSIZE > WARPSZE) 	__syncthreads();
			for (int k=0; k<LSPAN; k+=2) {		// compute a block of the matrix, write it in shared mem.
				yl[k*l_inc +j]     = y0;
				y0 = ak[2*k+3]*cost*y1 + ak[2*k+2]*y0;
				yl[(k+1)*l_inc +j] = y1;
				y1 = ak[2*k+5]*cost*y0 + ak[2*k+4]*y1;
				al += 4;
			}

			// transposed work:
			if (BLOCKSIZE > WARPSZE)	__syncthreads();
			real qlri[NFIELDS];	// accumulator
			const int itl = ll*l_inc + j % (BLOCKSIZE/(2*LSPAN));
			#pragma unroll
			for (int f=0; f<NFIELDS; f++)	qlri[f] = my_reo[f][0] * yl[itl];		// first element
			#pragma unroll
			for (int i=BLOCKSIZE/(2*LSPAN), k=1; i<BLOCKSIZE; i+= BLOCKSIZE/(2*LSPAN),k++) {		// accumulate
				#pragma unroll
				for (int f=0; f<NFIELDS; f++)	qlri[f] += my_reo[f][k] * yl[itl + i];
			}

			
			if (BLOCKSIZE/(2*LSPAN) <= WARPSZE) {		// reduce_add within same l is in same warp too:
				if (WARPSZE % (BLOCKSIZE/(2*LSPAN))) printf("ERROR\n");
				#pragma unroll
				for (int ofs = BLOCKSIZE/(LSPAN*4); ofs > 0; ofs>>=1) {
					#pragma unroll
					for (int f=0; f<NFIELDS; f++)	qlri[f] += shfl_down(qlri[f], ofs, BLOCKSIZE/(LSPAN*2));
				}
				if ( ((j % (BLOCKSIZE/(2*LSPAN))) == 0) && ((l+ll)<=llim) ) {	// write result
					if (nlat_2 <= BLOCKSIZE) {		// do we need atomic add or not ?
						#pragma unroll
						for (int f=0; f<NFIELDS; f++)	ql[2*(l+ll)+ri + f*ql_dist]   = qlri[f];
					} else {
						#pragma unroll
						for (int f=0; f<NFIELDS; f++)	atomicAdd(ql+2*(l+ll)+ri + f*ql_dist, qlri[f]);		// VERY slow atomic add on Kepler.
					}
				}
			} else {	// only partial reduction possible, finish with atomicAdd():
				if ((BLOCKSIZE/(2*LSPAN)) % WARPSZE) printf("ERROR\n");
				#pragma unroll
				for (int ofs = WARPSZE; ofs > 0; ofs>>=1) {
					#pragma unroll
					for (int f=0; f<NFIELDS; f++)	qlri[f] += shfl_down(qlri[f], ofs, WARPSZE);
				}
				if ( ((j % WARPSZE) == 0) && ((l+ll)<=llim) ) {	// write result
					#pragma unroll
					for (int f=0; f<NFIELDS; f++)	atomicAdd(ql+2*(l+ll)+ri + f*ql_dist, qlri[f]);		// VERY slow atomic add on Kepler.
				}
			}

			if (j<2*LSPAN) ak[j+2] = al[j];
			l+=LSPAN;
		}
	}
}

template<int S, int NFIELDS, typename real=double>
static void ileg_m_lowllim(shtns_cfg shtns, const real* q, real *ql, const int llim, int q_dist=0, int ql_dist=0)
{
	const int lmax = shtns->lmax;
	const int mres = shtns->mres;
	const int nlat_2 = shtns->nlat_2;
	const int nphi = shtns->nphi;
	int mmax = shtns->mmax;
	real *d_alm = (sizeof(real) >= 8) ? (real*) shtns->d_alm : (real*) shtns->d_alm_f;
	real *d_ct = (sizeof(real) >= 8) ? (real*) shtns->d_ct : (real*) shtns->d_ct_f;
	cudaStream_t stream = shtns->comp_stream;

	const int BLOCKSIZE = 256/NFIELDS;
	const int LSPAN_ = 8/NFIELDS;
	const int NW = 1;

	const int threadsPerBlock = BLOCKSIZE;	// can be from 32 to 1024, we should try to measure the fastest !
	const int blocksPerGrid = (nlat_2 + BLOCKSIZE*NW - 1) / (BLOCKSIZE*NW);
	if (q_dist == 0) q_dist = shtns->spat_stride;
	if (ql_dist == 0) ql_dist = shtns->nlm_stride;
	if (llim < mmax*mres) mmax = llim / mres;	// truncate mmax too !
	dim3 blocks(blocksPerGrid, mmax+1);
	dim3 threads(threadsPerBlock, 1);
	ileg_m_lowllim_kernel<BLOCKSIZE, LSPAN_, S, NFIELDS, real><<<blocks, threads, 0, stream>>>(d_alm, d_ct, (real*) q, (real*) ql, llim, nlat_2, lmax,mres, nphi, shtns->mpos_scale_analys, q_dist, ql_dist);
}


template<int BLOCKSIZE, int LSPAN, int S, typename real=double> __global__ void
ileg_m_highllim_kernel(const real *al, const real *ct, const real *q, real *ql, const int llim, const int nlat_2, const int lmax, const int mres, const int nphi, const real mpos_scale)
{
	const int it = BLOCKSIZE * blockIdx.x + threadIdx.x;
	const int j = threadIdx.x;
	const int im = blockIdx.y;
	const int m_inc = 2*nlat_2;
//    const int k_inc = 1;
	const real accuracy = (sizeof(real) >= 8) ? SHT_ACCURACY : SHT_ACCURACY_FLOAT;
	const real scale_factor = (sizeof(real) >= 8) ? SHT_SCALE_FACTOR : SHT_SCALE_FACTOR_FLOAT;

	__shared__ real ak[2*LSPAN+2];	// cache
	__shared__ real yl[LSPAN*BLOCKSIZE];
	__shared__ real reo[4*BLOCKSIZE];
	const int l_inc = BLOCKSIZE;
	const real cost = (it < nlat_2) ? ct[it] : 0.0;
	real y0, y1;


	if (im == 0) {
		if (j < 2*LSPAN+2) ak[j] = al[j];
		if (BLOCKSIZE > WARPSZE)	__syncthreads();
		y0 = (it < nlat_2) ? q[it] : 0.0;		// north
		y1 = (it < nlat_2) ? q[nlat_2*2-1 - it] : 0.0;	// south
		reo[j] = y0+y1;				// even
		reo[BLOCKSIZE +j] = y0-y1;		// odd

		int l = 0;
		y0 = (it < nlat_2) ? ct[it + nlat_2] : 0.0;		// weights are stored just after ct.
		if (S==1) y0 *= rsqrt(1.0 - cost*cost);
		y0 *= ak[0];
		y1 = y0 * ak[1] * cost;
		yl[j] = y0;
		yl[l_inc +j] = y1;
		al+=2;
		while (l <= llim) {
			for (int k=0; k<LSPAN; k+=2) {		// compute a block of the matrix, write it in shared mem.
				yl[k*l_inc +j]     = y0;
				y0 = ak[2*k+3]*cost*y1 + ak[2*k+2]*y0;
				yl[(k+1)*l_inc +j] = y1;
				y1 = ak[2*k+5]*cost*y0 + ak[2*k+4]*y1;
				al += 4;
			}
			if (BLOCKSIZE > WARPSZE)	__syncthreads();
			real qll = 0.0;	// accumulator
			// now re-assign each thread an l (transpose)
			const int ll = j / (BLOCKSIZE/LSPAN);
			for (int i=0; i<BLOCKSIZE; i+= BLOCKSIZE/LSPAN) {
				int it = j % (BLOCKSIZE/LSPAN) + i;
				qll += reo[(ll&1)*BLOCKSIZE +it] * yl[ll*l_inc +it];
			}

			// reduce_add within same l must be in same warp too:
			if (BLOCKSIZE/LSPAN > WARPSZE) printf("ERROR\n");

			for (int ofs = BLOCKSIZE/(LSPAN*2); ofs > 0; ofs>>=1) {
				qll += shfl_down(qll, ofs, BLOCKSIZE/LSPAN);
			}
			if ( ((j % (BLOCKSIZE/LSPAN)) == 0) && ((l+ll)<=llim) ) {	// write result
				if (nlat_2 <= BLOCKSIZE) {		// do we need atomic add or not ?
					ql[2*(l+ll)] = qll;
				} else {
					atomicAdd(ql+2*(l+ll), qll);		// VERY slow atomic add on Kepler.
				}
			}
			if (j<2*LSPAN) ak[j+2] = al[j];
			if (BLOCKSIZE > WARPSZE)	__syncthreads();
			l+=LSPAN;
		}
	} else {	// im > 0
		int m = im*mres;
		int l = (im*(2*(lmax+1)-(m+mres)))>>1;
		al += 2*(l+m);
		ql += 2*(l + S*im);	// allow vector transforms where llim = lmax+1

		if (j < 2*LSPAN+2) ak[j] = al[j];
		if (BLOCKSIZE > WARPSZE)	__syncthreads();
		const real sgn = 2*(j&1) - 1;	// -/+
		y0    = (it < nlat_2) ? q[im*m_inc + it] : 0.0;		// north imag (ani)
		real qer    = (it < nlat_2) ? q[(nphi-im)*m_inc + it] : 0.0;	// north real (an)
		y1    = (it < nlat_2) ? q[im*m_inc + nlat_2*2-1-it] : 0.0;	// south imag (asi)
		real qor    = (it < nlat_2) ? q[(nphi-im)*m_inc + nlat_2*2-1-it] : 0.0;	// south real (as)
		real qei = y0-qer;		qer += y0;		// ani = -qei[lane+1],   bni = qei[lane-1]
		real qoi = y1-qor;		qor += y1;		// bsi = -qoi[lane-1],   asi = qoi[lane+1];
		y0 = shfl_xor(qei, 1);	// exchange between adjacent lanes.
		y1 = shfl_xor(qoi, 1);
		reo[j] 			    = qer + qor;	// rer
		reo[BLOCKSIZE +j]   = qer - qor;	// ror
		reo[2*BLOCKSIZE +j] = sgn*(y0 - y1);	// rei
		reo[3*BLOCKSIZE +j] = sgn*(y0 + y1);	// roi

		y1 = sqrt(1.0 - cost*cost);	// sin(theta)

		y0 = mpos_scale;	// y0
		l = m - S;
		int ny = 0;
		int nsint = 0;
		do {		// sin(theta)^(m-S)		(use rescaling to avoid underflow)
			if (l&1) {
				y0 *= y1;
				ny += nsint;
				// the use of _any leads to wrong results. On KEPLER it is also slower.
				if (y0 < (accuracy+1.0/scale_factor)) {		// possible warp divergence
					ny--;
					y0 *= scale_factor;
				}
			}
			y1 *= y1;
			nsint += nsint;
			if (y1 < 1.0/scale_factor) {	// possible warp divergence
				nsint--;
				y1 *= scale_factor;
			}
		} while(l >>= 1);
		y0 *= ak[0];
		if (it < nlat_2)     y0 *= ct[it + nlat_2];		// include quadrature weights.
		y1 = ak[1]*y0*cost;


		l=m;		al+=2;
		while (l <= llim) {
			for (int k=0; k<LSPAN; k+=2) {		// compute a block of the matrix, write it in shared mem.
				yl[k*l_inc +j]     = (ny==0) ? y0 : 0.0;
				y0 = ak[2*k+3]*cost*y1 + ak[2*k+2]*y0;
				yl[(k+1)*l_inc +j] = (ny==0) ? y1 : 0.0;
				y1 = ak[2*k+5]*cost*y0 + ak[2*k+4]*y1;
				if (ny<0) {
					if (fabs(y0) > accuracy*scale_factor + 1.0)		// possible warp divergence
					{	// rescale when value is significant
						++ny;
						y0 *= 1.0/scale_factor;
						y1 *= 1.0/scale_factor;
					}
				}
				al += 4;
			}

			if (BLOCKSIZE > WARPSZE)	__syncthreads();
			real qlri = 0.0;	// accumulator
			// now re-assign each thread an l (transpose)
			const int ll = j / (BLOCKSIZE/LSPAN);
			const int ri = j / (BLOCKSIZE/(2*LSPAN)) % 2;	// real (0) or imag (1)
			if (ll+l <= llim) {
				for (int i=0; i<BLOCKSIZE; i+= BLOCKSIZE/(2*LSPAN)) {
				int it = j % (BLOCKSIZE/(2*LSPAN)) + i;
				qlri += reo[((ll&1)+2*ri)*BLOCKSIZE +it]   * yl[ll*l_inc +it];
				}
			}

			// reduce_add within same l must be in same warp too:
			if (BLOCKSIZE/(2*LSPAN) > WARPSZE) printf("ERROR\n");

			for (int ofs = BLOCKSIZE/(LSPAN*4); ofs > 0; ofs>>=1) {
				qlri += shfl_down(qlri, ofs, BLOCKSIZE/(LSPAN*2));
			}
			if ( ((j % (BLOCKSIZE/(2*LSPAN))) == 0) && ((l+ll)<=llim) ) {	// write result
				if (nlat_2 <= BLOCKSIZE) {		// do we need atomic add or not ?
				ql[2*(l+ll)+ri]   = qlri;
				} else {
				atomicAdd(ql+2*(l+ll)+ri, qlri);		// VERY slow atomic add on Kepler.
				}
			}
			if (j<2*LSPAN) ak[j+2] = al[j];
			if (BLOCKSIZE > WARPSZE)	__syncthreads();
			l+=LSPAN;
		}
	}
}

template<int S, int NFIELDS, typename real=double>
static void ileg_m_highllim(shtns_cfg shtns, const real* q, real *ql, const int llim, int q_dist=0, int ql_dist=0)
{
	const int lmax = shtns->lmax;
	const int mres = shtns->mres;
	const int nlat_2 = shtns->nlat_2;
	const int nphi = shtns->nphi;
	int mmax = shtns->mmax;
	real *d_alm = (sizeof(real) >= 8) ? (real*) shtns->d_alm : (real*) shtns->d_alm_f;
	real *d_ct = (sizeof(real) >= 8) ? (real*) shtns->d_ct : (real*) shtns->d_ct_f;
	cudaStream_t stream = shtns->comp_stream;

	const int BLOCKSIZE = 256/NFIELDS;
	const int LSPAN_ = 8/NFIELDS;
	const int NW = 1;

	const int threadsPerBlock = BLOCKSIZE;	// can be from 32 to 1024, we should try to measure the fastest !
	const int blocksPerGrid = (nlat_2 + BLOCKSIZE*NW - 1) / (BLOCKSIZE*NW);
	if (q_dist == 0) q_dist = shtns->spat_stride;
	if (ql_dist == 0) ql_dist = shtns->nlm_stride;
	if (llim < mmax*mres) mmax = llim / mres;	// truncate mmax too !
	dim3 blocks(blocksPerGrid, mmax+1);
	dim3 threads(threadsPerBlock, 1);
	for (int f=0; f<NFIELDS; f++) {
		ileg_m_highllim_kernel<BLOCKSIZE, LSPAN_, S, real><<<blocks, threads, 0, stream>>>(d_alm, d_ct, q + f*q_dist, ql + f*ql_dist, llim, nlat_2, lmax,mres, nphi, shtns->mpos_scale_analys);
	}
}


template<int S, int NFIELDS, typename real=double>
static void legendre(shtns_cfg shtns, const real *ql, real *q, const int llim, const int mmax, int spat_dist = 0)
{
	if (spat_dist == 0) spat_dist = shtns->spat_stride;
	if (mmax==0) {
		leg_m0<S,NFIELDS,real>(shtns, ql, q, llim);
	} else {
		const int limit = (sizeof(real) >= 8) ? SHT_L_RESCALE_FLY : SHT_L_RESCALE_FLY_FLOAT;
		if (llim <= limit) {
			leg_m_lowllim<S,NFIELDS,real>(shtns, ql, q, llim, mmax, spat_dist);
		} else {
			leg_m_highllim<S,NFIELDS,real>(shtns, ql, q, llim, mmax);
		}
	}
}

/// Perform SH transform on data that is already on the GPU. d_Qlm and d_Vr are pointers to GPU memory (obtained by cudaMalloc() for instance)
template<int S, int NFIELDS, typename real=double>
static void ilegendre(shtns_cfg shtns, const real *q, real* ql, const int llim, int spat_dist = 0)
{
	int mmax = shtns->mmax;
	const int mres = shtns->mres;

	if (spat_dist == 0) spat_dist = shtns->spat_stride;
	cudaMemsetAsync(ql, 0, sizeof(real) * NFIELDS * shtns->nlm_stride, shtns->comp_stream);		// set to zero before we start.
	if (llim < mmax*mres) mmax = llim / mres;	// truncate mmax too !
	if (mmax==0) {
		ileg_m0<S, NFIELDS, real>(shtns, q, ql, llim, spat_dist, shtns->nlm_stride);
	} else {
		const int limit = (sizeof(real) >= 8) ? SHT_L_RESCALE_FLY : SHT_L_RESCALE_FLY_FLOAT;
		if (llim <= limit) {
			ileg_m_lowllim<S, NFIELDS, real>(shtns, q, ql, llim, spat_dist, shtns->nlm_stride);
		} else {
			ileg_m_highllim<S, NFIELDS, real>(shtns, q, ql, llim, spat_dist, shtns->nlm_stride);
		}
	}
}
