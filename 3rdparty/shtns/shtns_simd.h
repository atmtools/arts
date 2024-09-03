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
 
/****************************************************************************
 * SIMD macros and functions for processor agnostic vectorization of SHTns. *
 *    A subset also adapts to various vector-length (2, 4 or 8 doubles).    *
 *    Written by Nathanael Schaeffer / CNRS                                 *
 ****************************************************************************/

/// define _GCC_VEC_ to activate SIMD using gcc extensions.
#if _GCC_VEC_ == 0
	#undef _GCC_VEC_
#endif

/* are there supported vector extensions available ? */
#if !(defined __SSE2__ || defined __VSX__ || __ARM_NEON_FP >= 8)
	#undef _GCC_VEC_
#endif
#ifdef __INTEL_COMPILER
	#if __INTEL_COMPILER < 1400
		#undef _GCC_VEC_
		#warning "no vector extensions available ! use gcc 4+ or icc 14+ for best performance."
	#endif
#endif
#ifdef __GNUC__
	#if __GNUC__ < 4
		#undef _GCC_VEC_
		#warning "no vector extensions available ! use gcc 4+ or icc 14+ for best performance."
	#endif
#endif


#if _GCC_VEC_ && ( __ARM_NEON_FP >= 8  )
	// support ARM NEON
	#include <arm_neon.h>
	#define MIN_ALIGNMENT 16
	#define VSIZE 2
	typedef float64x2_t s2d;
	typedef float64x2_t v2d;
	typedef float64x2_t rnd;
	#define VSIZE2 2
	#define _SIMD_NAME_ "neon"
	#warning "arm neon"
	#define vall(x) vdupq_n_f64(x)
	// _mm_unpacklo_pd => vzip1q_f64 	// _mm_unpackhi_pd => vzip2q_f64

	#define vread(mem, idx) ((s2d*)(mem))[idx]
	#define vstor(mem, idx, v) ((s2d*)(mem))[idx] = v
	#define vread2 vread
	#define vstor2 vstor
	#define vxor2(v,x) veorq_u64(v,x)
	inline static v2d v2d_reduce(v2d a, v2d b) { return vpaddq_f64(a,b); }
	inline static v2d vneg_even_precalc(v2d v) {		// don't use in an intensive loop.
		v[0] = -v[0];
		return v;
	}
	inline static v2d vreverse(v2d a) { return vextq_f64(a,a,1); }
	#define vxchg(a) vreverse(a)
	#define vreverse_pairs(v) (v)
	#define vdup_even(v) vdupq_laneq_f64(v,0)
	#define vdup_odd(v) vdupq_laneq_f64(v,1)
	#define vxchg_even_odd(v) vreverse(v)

	static const unsigned long long _neg0[2] __attribute__((aligned (16))) = {0x8000000000000000ULL, 0} ;		// a constant needed to change the sign of vectors

	#define vblend_even_odd(a,b) vsetq_lane_f64(vgetq_lane_f64(b,1),a,1)	//	__builtin_shufflevector(a,b,0,3)
	#define vneg_even_xor_cte (*(const v2d*)_neg0)
	#define vxor(v,x) vxor2(v,x)
	#define reduce_add(a) vaddvq_f64(a)
	#define vinterleave(a,b) {  rnd x = vzip1q_f64(a,b);	b = vzip2q_f64(a, b);	a = x; }
	#define vinterleave_reverse(a,b)	{  rnd x = vzip1q_f64(a,b);	a = vzip2q_f64(a, b);	b = x; }

	inline static void S2D_CSTORE2_4MAGIC(double* mem, long idx, rnd nr, rnd sr, rnd ni, rnd si) {
		((s2d*)mem)[idx*4]   = vzip1q_f64(nr, ni);	// aa = north_ri[0]
		((s2d*)mem)[idx*4+1] = vzip1q_f64(sr, si);	// cc = south_ri[0]
		((s2d*)mem)[idx*4+2] = vzip2q_f64(nr, ni);	// bb = north_ri[1]
		((s2d*)mem)[idx*4+3] = vzip2q_f64(sr, si);	// dd = south_ri[1]
	}
	inline static v2d IxKxZ(double k, v2d z) {		// I*k*z,  allowing to use FMA.
		v2d vk = {-k,k};
		return vk * vxchg(z);
	}

	// vset(lo, hi) takes two doubles and pack them in a vector
	inline static v2d vset(double lo, double hi) {
		v2d v = {lo, hi};
		return v;
	}
	// vdup(x) takes a double and duplicate it to a vector of 2 doubles.
	#define vdup(x) vall(x)
	// vxchg(a) exchange hi and lo component of vector a
	#define vlo_to_cplx(a) vsetq_lane_f64(0.0,a,1)
	#define vhi_to_cplx(a) vsetq_lane_f64(vgetq_lane_f64(a,1),vall(0.0),0)
	#define vcplx_real(a) vlo_to_dbl(a)
	#define vcplx_imag(a) vhi_to_dbl(a)
	#ifdef __clang__
		// allow to compile with clang (llvm)
		#define vlo(a) (a)[0]
		#define vlo_to_dbl(a) (a)[0]
		#define vhi_to_dbl(a) (a)[1]
	#else
		// gcc extensions
		#define vlo(a) vgetq_lane_f64(a,0)
		#define vlo_to_dbl(a) vgetq_lane_f64(a,0)
		#define vhi_to_dbl(a) vgetq_lane_f64(a,1)
	#endif
	#define v2d_lo(a) (a)
#endif

#if _GCC_VEC_ && __VSX__
	// support VSX (IBM Power)
	#include <altivec.h>
	#define MIN_ALIGNMENT 16
	#define VSIZE 2
	typedef __vector double s2d;
	typedef __vector double v2d;
	typedef __vector double rnd;
	#define VSIZE2 2
	#define _SIMD_NAME_ "vsx"
	#define vall(x) vec_splats((double)x)
	inline static s2d vreverse(s2d a) {
		const vector unsigned char perm = { 8,9,10,11,12,13,14,15, 0,1,2,3,4,5,6,7 };
		return vec_perm(a,a,perm);
	}
	//#define vreverse(a) vec_reve(a)
	#define vxchg(a) vreverse(a)
	#define vdup_even(v) vec_mergeh(v,v)
	#define vdup_odd(v)  vec_mergel(v,v)
	#define vxchg_even_odd(a) vreverse(a)
	#define vread(mem, idx) ((s2d*)(mem))[idx]
	#define vstor(mem, idx, v) ((s2d*)(mem))[idx] = v
	#define vread2 vread
	#define vstor2 vstor
	static const unsigned long long _neg0[2] __attribute__((aligned (16))) = {0x8000000000000000ULL, 0} ;		// a constant needed to change the sign of vectors
	#define vneg_even_xor_cte (*(const v2d*)_neg0)
	#define vxor(v,x) vec_xor(v,x)
	#define vxor2 vxor
	inline static double reduce_add(rnd v) { return vec_extract(v,0) + vec_extract(v,1); }
	inline static v2d v2d_reduce(rnd a, rnd b) {
		v2d c = vec_mergel(a, b);		b = vec_mergeh(a, b);
		return b + c;
	}
	inline static s2d vneg_even_precalc(s2d a) {
		const s2d ne = {-1.0, 1.0};
		return ne * a;
	}
	#define vinterleave(a,b) { rnd x = vec_mergeh(a,b);		b = vec_mergel(a,b);	a = x; }
	#define vinterleave_reverse(a,b)	{ rnd x = vec_mergeh(a,b);		a = vec_mergel(a,b);	b = x; }
/*	inline static v2d addi(v2d a, v2d b) {
		const s2d mp = {-1.0, 1.0};
		return a + vxchg(b)*mp;
	}
	inline static s2d subadd(s2d a, s2d b) {
		const s2d mp = {-1.0, 1.0};
		return a + b*mp;
	}
*/
	#define vdup(x) vec_splats((double)x)
	#define vlo(a) vec_extract(a, 0)
	#define vhi(a) vec_extract(a, 1)
	#define vcplx_real(a) vec_extract(a, 0)
	#define vcplx_imag(a) vec_extract(a, 1)
	#define vlo_to_dbl(a) vec_extract(a, 0)
	#define vhi_to_dbl(a) vec_extract(a, 1)
	#define vreverse_pairs(v) (v)
	#define v2d_lo(a) (a)
	// vset(lo, hi) takes two doubles and pack them in a vector
	//#define vset(lo, hi) _mm_set_pd(hi, lo)
	//#define vlo_to_cplx(a) vec_insert(0.0, a, 1)
	//#define vhi_to_cplx(a) vec_mergel(a, vdup(0.0))

	inline static s2d vblend_even_odd(s2d a, s2d b) {	// same as _mm_shuffle_pd(a,b,2)
		const vector unsigned char perm = {0,1,2,3,4,5,6,7, 24,25,26,27,28,29,30,31};
		return vec_perm(a,b,perm);
	}
	inline static v2d IxKxZ(double k, v2d z) {		// I*k*z,  allowing to use FMA.
		const s2d vk = {-k, k};
		return vk*vxchg(z);
	}

		inline static void S2D_CSTORE2_4MAGIC(double* mem, long idx, rnd nr, rnd sr, rnd ni, rnd si) {
			((s2d*)mem)[idx*4]   = vec_mergeh(nr, ni);
			((s2d*)mem)[idx*4+1] = vec_mergeh(sr, si);
			((s2d*)mem)[idx*4+2] = vec_mergel(nr, ni);
			((s2d*)mem)[idx*4+3] = vec_mergel(sr, si);
		}
#endif



#if _GCC_VEC_ && __SSE2__
	#define VSIZE 2
	typedef double s2d __attribute__ ((vector_size (8*VSIZE)));		// vector that should behave like a real scalar for complex number multiplication.
	typedef double v2d __attribute__ ((vector_size (8*VSIZE)));		// vector that contains a complex number
	#define vxchg(a) ((v2d)_mm_shuffle_pd(a,a,1))					// swap the two elements of a vector of 2 doubles.
	#define vxor2(v,x) ((v2d)_mm_xor_pd(v, x))
	#define vread2(mem, idx) ((v2d)_mm_loadu_pd( ((double*)(mem)) + (idx)*2 ))
	#define vstor2(mem, idx, v) _mm_storeu_pd( ((double*)(mem)) + (idx)*2 , v)
	static const unsigned long long _neg0[2] __attribute__((aligned (16))) = {0x8000000000000000ULL, 0} ;		// a constant needed to change the sign of vectors
	#ifdef __AVX__
		#include <immintrin.h>
		typedef double v4d __attribute__ ((vector_size (8*4)));		// vector that contains 2 complex numbers
		#define vall4(x) ((v4d) _mm256_set1_pd(x))
		#define vread4(mem, idx) ((v4d)_mm256_loadu_pd( ((double*)(mem)) + (idx)*4 ))
		#define vstor4(mem, idx, v) _mm256_storeu_pd( ((double*)(mem)) + (idx)*4 , v)
		inline static v4d v2d_x2_to_v4d(v2d a, v2d b) {
			return (v4d) _mm256_insertf128_pd( _mm256_castpd128_pd256( a ), b, 1);
		}
		inline static v4d vreverse4(v4d a) {		// reverse vector: [0,1,2,3] => [3,2,1,0]
			#if defined( __AVX2__ )
			return (v4d) _mm256_permute4x64_pd(a, 0x1B);	// 3 cycles on intel; 6 cycles on Zen2
			#else
			a = (v4d)_mm256_permute2f128_pd(a,a, 1);	// => [2,3,0,1]			// 2 cycles on SandyBridge, 3 cycles on Haswell+, 3 cycles on Zen2
			return (v4d)_mm256_shuffle_pd(a, a, 5);		// [2,3,0,1] => [3,2,1,0]	// 1 cycle on intel; 3 cycles on Zen2
			#endif
		}
		#define vdup_even4(v) ((v4d)_mm256_movedup_pd(v))
	#endif
	#ifdef __AVX512F__
		#define MIN_ALIGNMENT 64
		#define VSIZE2 8
		// Allocate memory aligned on 64 bytes for AVX-512
		#define _SIMD_NAME_ "avx512"
		typedef double rnd __attribute__ ((vector_size (VSIZE2*8)));		// vector of 8 doubles.
		#define vall(x) ((rnd) _mm512_set1_pd(x))
		#define vread(mem, idx) ((rnd)_mm512_loadu_pd( ((double*)(mem)) + (idx)*8 ))
		#define vstor(mem, idx, v) _mm512_storeu_pd( ((double*)(mem)) + (idx)*8 , v)
		inline static rnd vreverse(rnd a) {		// reverse vector: [0,1,2,3,4,5,6,7] => [7,6,5,4,3,2,1,0]
			a = _mm512_permute_pd(a,0x55);	// [1,0,3,2,5,4,7,6]
			return _mm512_shuffle_f64x2(a,a,0x1B);	// [7,6,5,4,3,2,1,0]
			//return (rnd) _mm512_permutexvar_pd(_mm512_set_epi64(0,1,2,3,4,5,6,7), a);		// same speed on KNL, but requires a constant to be loaded...
		}
		#define vreverse_pairs(v) ((rnd)_mm512_shuffle_f64x2(v,v, 0x1B))
		#define vdup_even(v) ((rnd)_mm512_movedup_pd(v))
		#define vdup_odd(v)	 ((rnd)_mm512_permute_pd(v,0xFF))
		#define vxchg_even_odd(v) ((rnd)_mm512_permute_pd(v,0x55))
		#define vblend_even_odd(a,b) ((rnd)_mm512_mask_blend_pd((__mmask8) 0xAA, a,b))
		inline static rnd vneg_even_precalc(rnd v) {		// don't use in an intesive loop.
			return _mm512_fmaddsub_pd(vall(0.0), vall(0.0), v);
		}
		#define vneg_even_xor_cte ((rnd)_mm512_castsi512_pd(_mm512_broadcast_i32x4(*(const __m128i*)_neg0)))
		#define vxor(v,x) ((rnd)_mm512_castsi512_pd( _mm512_xor_epi64(_mm512_castpd_si512(v), _mm512_castpd_si512(x))))
		inline static double reduce_add(rnd a) {
			return _mm512_reduce_add_pd(a);
		}
		/*	inline static v2d v2d_reduce(rnd a, rnd b) {	// KNL Latency=39
			rnd x = (rnd)_mm512_unpackhi_pd(a, b) + (rnd)_mm512_unpacklo_pd(a, b);		// Latency 15c
			v4d y = (v4d)_mm512_castpd512_pd256(x) + (v4d)_mm512_extractf64x4_pd(x,1);	// Latency 12c
			return (v2d)_mm256_castpd256_pd128(y) + (v2d)_mm256_extractf128_pd(y,1);	// Latency 12c
		}	*/
		inline static v2d v2d_reduce(rnd a, rnd b) {	// KNL Latency=37
			rnd x = (rnd)_mm512_shuffle_pd(a,b,0x55);	// a1,b0,a3,b2,...	// L=7
			a = vblend_even_odd(a,b);		// a0,b1,a2,b3,...  // (L=2, hidden)
			x+=a;		// L=13
			v4d y = (v4d)_mm512_castpd512_pd256(x) + (v4d)_mm512_extractf64x4_pd(x,1);	// Latency 12c
			return (v2d)_mm256_castpd256_pd128(y) + (v2d)_mm256_extractf128_pd(y,1);	// Latency 12c
		}
		inline static v4d v4d_reduce(rnd a, rnd b, rnd c, rnd d) {		// KNL Latency=40c
			rnd x = (rnd)_mm512_shuffle_pd(a,b,0x55);	// a1,b0,a3,b2,...	// L=7
			a = vblend_even_odd(a,b);		// a0,b1,a2,b3,...  // (L=2, hidden)
			rnd y = (rnd)_mm512_shuffle_pd(c,d,0x55);
			c = vblend_even_odd(c,d);
			x += a;		y += c;
			v4d a4 = (v4d)_mm512_castpd512_pd256(x) + (v4d)_mm512_extractf64x4_pd(x,1);
			v4d c4 = (v4d)_mm512_castpd512_pd256(y) + (v4d)_mm512_extractf64x4_pd(y,1);
			v4d y4 = _mm256_permute2f128_pd(a4,c4,0x21);	// _mm256_setr_pd(a[2],a[3],c[0],c[1]);
			v4d x4 = _mm256_blend_pd(a4,c4,0xc);			// _mm256_setr_pd(a[0],a[1],c[2],c[3]);
			return x4+y4;
		}

		// inplace operation: a=[a0,b0,a1,b1, a2,b2,a3,b3]  b=[a4,b4,a5,b5, a6,b6,a7,b7]
		/*#define vinterleave(a,b) { \
			a = (rnd)_mm512_permutex_pd(a, 0xD8);           b = (rnd)_mm512_permutex_pd(b, 0xD8);   \
			rnd ev = (rnd)_mm512_unpacklo_pd(a, b);         rnd od = (rnd)_mm512_unpackhi_pd(a, b); \
			a = _mm512_shuffle_f64x2(ev,od, 0x44);          b = _mm512_shuffle_f64x2(ev,od, 0xEE); } */
		#define vinterleave(a,b) { /*
		*/	a = (rnd)_mm512_permutex_pd(a, 0xD8);	 /* a0,a2,a1,a3, a4,a6,a5,a7
		*/	b = (rnd)_mm512_permutex_pd(b, 0x72);	 /* b2,b0,b3,b1, b6,b4,b7,b5
		*/	rnd od = (rnd)_mm512_shuffle_pd(a, b, 0x55);				/* a2,b2,a3,b3, a6,b6,a7,b7
		*/	rnd ev = (rnd)_mm512_mask_blend_pd((__mmask8) 0xAA, a,b);	/* a0,b0,a1,b1, a4,b4,a5,b5
		*/	a = (rnd)_mm512_shuffle_f64x2(ev,od, 0x44);	  /* a0,b0,a1,b1, a2,b2,a3,b3
		*/  b = (rnd)_mm512_shuffle_f64x2(ev,od, 0xEE); } /* a4,b4,a5,b5, a6,b6,a7,b7 */
		/*#define vinterleave(a,b) { \
			rnd x = _mm512_permutex2var_pd(a, _mm512_setr_epi64(0,8,1,9,2,10,3,11), b);	\
			b = _mm512_permutex2var_pd(a, _mm512_setr_epi64(4,12,5,13,6,14,7,15), b);	a = x;	}	*/
		/*#define vinterleave_reverse(a,b) { \
			a = (rnd)_mm512_permutex_pd(a, 0x8D);		b = (rnd)_mm512_permutex_pd(b, 0x8D);	\
			rnd ev = (rnd)_mm512_unpacklo_pd(a, b);		rnd od = (rnd)_mm512_unpackhi_pd(a, b);	\
			b = _mm512_shuffle_f64x2(od,ev, 0x44);		a = _mm512_shuffle_f64x2(od,ev, 0xEE); }	*/
		#define vinterleave_reverse(a,b) { /*
		*/	a = (rnd)_mm512_permutex_pd(a, 0x8D);	 /* a1,a3,a0,a2, a5,a7,a4,a6
		*/	b = (rnd)_mm512_permutex_pd(b, 0x27);	 /* b3,b1,b2,b0, b7,b5,b6,b4
		*/	rnd od = (rnd)_mm512_shuffle_pd(a, b, 0x55);				/* a3,b3,a2,b2, a7,b7,a6,b6
		*/	rnd ev = (rnd)_mm512_mask_blend_pd((__mmask8) 0xAA, a,b);	/* a1,b1,a0,b0, a5,b5,a4,b4
		*/	b = (rnd)_mm512_shuffle_f64x2(od,ev, 0x44);	  /* a3,b3,a2,b2, a1,b1,a0,b0
		*/  a = (rnd)_mm512_shuffle_f64x2(od,ev, 0xEE); } /* a7,b7,a6,b6, a5,b5,a4,b4 */

		/* Some AVX512 vector tricks:
		 * pair-wise exchange: _mm512_permute_pd(a, 0x55)	=> [1,0,3,2,5,4,7,6]
		 * exchange middle elements of quartets: _mm512_permutex_pd(a, 0xD8)  => [0,2,1,3,4,6,5,7]
		 * exchange middle elements of quartets + swap pairs: _mm512_permutex_pd(a, 0x8D)	=> [1,3,0,2,5,7,4,6]
		 * reverse vector of pairs: _mm512_shuffle_f64x2(a,a, 0x1B)	=> [6,7,4,5,2,3,0,1]
		 */

	/*	inline static void S2D_CSTORE2_4MAGIC_OLD(double* mem, long idx, rnd nr, rnd sr, rnd ni, rnd si) {
			rnd aa = (rnd)_mm512_unpacklo_pd(nr, sr);	rnd bb = (rnd)_mm512_unpackhi_pd(nr, sr);	// [n0,s0,n2,s2,n4,s4,n6,s6] and [n1,s1,n3,s3,n5,s5,n7,s7]
			rnd cc = (rnd)_mm512_unpacklo_pd(ni, si);	rnd dd = (rnd)_mm512_unpackhi_pd(ni, si);	// same with imaginary part
			aa = _mm512_permutex_pd(aa, 0xD8);		bb = _mm512_permutex_pd(bb, 0xD8);	// [n0,n2,s0,s2,n4,n6,s4,s6] and [n1,n3,s1,s3,n5,n7,s5,s7]
			cc = _mm512_permutex_pd(cc, 0xD8);		dd = _mm512_permutex_pd(dd, 0xD8);	// same with imaginary part
			nr = (rnd)_mm512_unpacklo_pd(aa, cc);	// [n0,s0,n4,s4] packed real/imag
			sr = (rnd)_mm512_unpackhi_pd(aa, cc);	// [n2,s2,n6,s6] packed real/imag
			ni = (rnd)_mm512_unpacklo_pd(bb, dd);	// [n1,s1,n5,s5] packed real/imag
			si = (rnd)_mm512_unpackhi_pd(bb, dd);	// [n3,s3,n7,s7] packed real/imag
			_mm512_storeu_pd(mem + idx*32,     _mm512_shuffle_f64x2(nr,ni, 0x44) ); //[n0,s0,n1,s1]
			_mm512_storeu_pd(mem + idx*32 +8,  _mm512_shuffle_f64x2(sr,si, 0x44) ); //[n2,s2,n3,s3]
			_mm512_storeu_pd(mem + idx*32 +16, _mm512_shuffle_f64x2(nr,ni, 0xEE) ); //[n4,s4,n5,s5]
			_mm512_storeu_pd(mem + idx*32 +24, _mm512_shuffle_f64x2(sr,si, 0xEE) ); //[n6,s6,n7,s7]
		}	*/
		// FASTER ALTERNATIVE (easier on shuffle port)
		inline static void S2D_CSTORE2_4MAGIC(double* mem, long idx, rnd nr, rnd sr, rnd ni, rnd si) {
			nr = (rnd)_mm512_permutex_pd(nr, 0x72);		sr = (rnd)_mm512_permutex_pd(sr, 0x27);		//	nr:[2,0,3,1, 6,4,7,5], sr:[3,1,2,0, 7,5,6,4]
			ni = (rnd)_mm512_permutex_pd(ni, 0xD8);		si = (rnd)_mm512_permutex_pd(si, 0x8D);		//  ni:[0,2,1,3, 4,6,5,7], si:[1,3,0,2, 5,7,4,6]
			rnd ar = (rnd)_mm512_mask_blend_pd((__mmask8) 0xCC, nr,sr);		// [n2,n0,s2,s0,  n6,n4,s6,s4]
			rnd br = (rnd)_mm512_mask_blend_pd((__mmask8) 0xCC, sr,nr);		// [s3,s1,n3,n1,  s7,s5,n7,n5]
			rnd ai = (rnd)_mm512_mask_blend_pd((__mmask8) 0xCC, ni,si);		// [n0,n2,s0,s2,  n4,n6,s4,s6]
			rnd bi = (rnd)_mm512_mask_blend_pd((__mmask8) 0xCC, si,ni);		// [s1,s3,n1,n3,  s5,s7,n5,n7]
			nr = (rnd)_mm512_shuffle_f64x2(ar,br, 0x14);		// n2,n0,s2,s0, n3,n1,s3,s1
			ni = (rnd)_mm512_shuffle_f64x2(ai,bi, 0x14);		// n0,n2,s0,s2, n1,n3,s1,s3
			sr = (rnd)_mm512_shuffle_f64x2(ar,br, 0xBE);		// n6,n4,s6,s4, n7,n5,s7,s5
			si = (rnd)_mm512_shuffle_f64x2(ai,bi, 0xBE);		// n4,n6,s4,s6, n5,n7,s5,s7
			vstor(mem, idx*4,   _mm512_shuffle_pd(nr,ni, 0x55) ); //[n0,s0,n1,s1]
			vstor(mem, idx*4+1, vblend_even_odd(nr,ni) );		  //[n2,s2,n3,s3]
			vstor(mem, idx*4+2, _mm512_shuffle_pd(sr,si, 0x55) ); //[n4,s4,n5,s5]
			vstor(mem, idx*4+3, vblend_even_odd(sr,si) );		  //[n6,s6,n7,s7]
		}

	#elif defined __AVX__
		#define MIN_ALIGNMENT 32
		#define VSIZE2 4
		#ifdef __AVX2__
			#define _SIMD_NAME_ "avx2"
		#else
			#define _SIMD_NAME_ "avx"
		#endif
		typedef double rnd __attribute__ ((vector_size (VSIZE2*8)));		// vector of 4 doubles.
		#define vall(x) ((rnd) _mm256_set1_pd(x))
		#define vread(mem, idx) ((rnd)_mm256_loadu_pd( ((double*)(mem)) + (idx)*4 ))
		#define vstor(mem, idx, v) _mm256_storeu_pd( ((double*)(mem)) + (idx)*4 , v)
		#define vreverse vreverse4
		#define vreverse_pairs(v) ((rnd)_mm256_permute2f128_pd(v,v,1))
		#define vdup_even(v) ((rnd)_mm256_movedup_pd(v))
		#define vdup_odd(v)  ((rnd)_mm256_unpackhi_pd(v,v))
		// intra-lane exchange, same as _mm256_permute_pd(xm, 0x5), but faster with gcc<10
		#define vxchg_even_odd(v) ((rnd)_mm256_shuffle_pd(v,v,0x5))
		#define vblend_even_odd(a,b) ((rnd)_mm256_blend_pd (a,b,0xA))
		//#define vblend_even_odd(a,b) ((rnd)_mm256_shuffle_pd (a,b,0xA))
		inline static rnd vneg_even_precalc(rnd v) {		// don't use in an intesive loop.
			return _mm256_addsub_pd(vall(0.0), v);
		}
		#define vneg_even_xor_cte ((rnd)_mm256_broadcast_pd((const __m128d*)_neg0))
		//#define vneg_even_xor_cte ((rnd)_mm256_castsi256_pd( _mm256_setr_epi32(0,0x80000000, 0,0, 0,0x80000000, 0,0)))	// BUGGY ON GCC! DON'T USE!
		#define vxor(v,x) ((rnd)_mm256_xor_pd(v, x))
		inline static double reduce_add(rnd a) {	// Latency=12c Skylake
			v2d t = (v2d)_mm256_castpd256_pd128(a) + (v2d)_mm256_extractf128_pd(a,1);
			return _mm_cvtsd_f64(t) + _mm_cvtsd_f64(_mm_unpackhi_pd(t,t));
		}
		inline static v2d v2d_reduce(rnd a, rnd b) {		// Latency=13c Skylake.
			a = _mm256_hadd_pd(a, b);
			return (v2d)_mm256_castpd256_pd128(a) + (v2d)_mm256_extractf128_pd(a,1);
		}
		inline static v4d v4d_reduce(rnd a, rnd b, rnd c, rnd d) {		// Latency=15c Skylake
			a = _mm256_hadd_pd(a, b);	// a0+a1, b0+b1, a2+a3, b2+b3
			c = _mm256_hadd_pd(c, d);	// c0+c1, d0+d1, c2+c3, d2+d3
		/*	rnd x4 = (rnd)_mm256_shuffle_pd(a,b,5);	// a1,b0,a3,b2	// L=1
			a = (rnd)_mm256_blend_pd(a, b, 0xA);	// a0,b1,a2,b3,...  // (L=1, paired, hidden)
			rnd y4 = (rnd)_mm256_shuffle_pd(c,d,5);
			c = (rnd)_mm256_blend_pd(c, d, 0xA);
			a += x4;		c += y4;	*/
			v4d y = _mm256_permute2f128_pd(a,c,0x21);	// _mm256_setr_pd(a[2],a[3],c[0],c[1]);
			v4d x = _mm256_blend_pd(a,c,0xc);			// _mm256_setr_pd(a[0],a[1],c[2],c[3]);
			return x+y;
		}

		// inplace operation: a=[a0,b0,a1,b1] b=[a2,b2,a3,b3]
		#define vinterleave(a,b) { \
			rnd ev = (rnd)_mm256_unpacklo_pd(a, b);		rnd od = (rnd)_mm256_unpackhi_pd(a, b);	\
			a = _mm256_permute2f128_pd(ev, od, 0x20);	b = _mm256_permute2f128_pd(ev, od, 0x31); }
		#define vinterleave_reverse(a,b) { \
			rnd ev = (rnd)_mm256_unpacklo_pd(a, b);		rnd od = (rnd)_mm256_unpackhi_pd(a, b);	\
			b = _mm256_permute2f128_pd(od, ev, 0x20);	a = _mm256_permute2f128_pd(od, ev, 0x31); }

		#define vinterleave_x4(nr, sr, ni, si) { \
			rnd ar = (rnd)_mm256_permute2f128_pd(nr,sr,0x20);	rnd ai = (rnd)_mm256_permute2f128_pd(ni,si,0x20);	\
			rnd br = (rnd)_mm256_permute2f128_pd(nr,sr,0x31);	rnd bi = (rnd)_mm256_permute2f128_pd(ni,si,0x31);	\
			nr = _mm256_unpacklo_pd(ar, ai);	ni = _mm256_unpackhi_pd(ar, ai);  \
			sr = _mm256_unpacklo_pd(br, bi);	si = _mm256_unpackhi_pd(br, bi); }

		inline static void S2D_CSTORE2_4MAGIC(double* mem, long idx, rnd nr, rnd sr, rnd ni, rnd si) {
			rnd ar = (rnd)_mm256_permute2f128_pd(nr,sr,0x20);	rnd ai = (rnd)_mm256_permute2f128_pd(ni,si,0x20);	// n0,n1,s0,s1
			rnd br = (rnd)_mm256_permute2f128_pd(nr,sr,0x31);	rnd bi = (rnd)_mm256_permute2f128_pd(ni,si,0x31);	// n2,n3,s2,s3
			_mm256_storeu_pd(mem + idx*16,     _mm256_unpacklo_pd(ar, ai));		// nr0,ni0, sr0,si0
			_mm256_storeu_pd(mem + idx*16 +4,  _mm256_unpackhi_pd(ar, ai));		// nr1,ni1, sr1,si1
			_mm256_storeu_pd(mem + idx*16 +8,  _mm256_unpacklo_pd(br, bi));		// nr2,ni2, sr2,si2
			_mm256_storeu_pd(mem + idx*16 +12, _mm256_unpackhi_pd(br, bi));		// nr3,ni3, sr3,si3
		}
	#else
		#define MIN_ALIGNMENT 16
		#define VSIZE2 2
		typedef double rnd __attribute__ ((vector_size (VSIZE2*8)));		// vector of 2 doubles.
		// Allocate memory aligned on 16 bytes for SSE2 (fftw_malloc works only if fftw was compiled with --enable-sse2)
		// in 64 bit systems, malloc should be 16 bytes aligned anyway.
		#define vall(x) ((rnd) _mm_set1_pd(x))
		#define vread(mem, idx) ((s2d*)(mem))[idx]
		#define vstor(mem, idx, v) ((s2d*)(mem))[idx] = v
		#ifdef __SSE3__
			#include <pmmintrin.h>
			#define _SIMD_NAME_ "sse3"
			inline static v2d v2d_reduce(v2d a, v2d b) {
				return _mm_hadd_pd(a,b);
			}
			inline static rnd vneg_even_precalc(rnd v) {		// don't use in an intesive loop.
				return _mm_addsub_pd(vall(0.0), v);
			}
		#else
			#include <emmintrin.h>
			#define _SIMD_NAME_ "sse2"
			inline static v2d v2d_reduce(v2d a, v2d b) {
				v2d c = _mm_unpacklo_pd(a, b);		b = _mm_unpackhi_pd(a, b);
				return b + c;
			}
			inline static rnd vneg_even_precalc(rnd v) {		// don't use in an intesive loop.
				rnd nv = vall(0.0) - v;
				return _mm_shuffle_pd(nv, v, 2);
			}
		#endif
		inline static rnd vreverse(rnd a) {	return (rnd)_mm_shuffle_pd(a,a,1);	}
		#define vreverse_pairs(v) (v)
		#define vdup_even(v) ((rnd)_mm_unpacklo_pd(v,v))
		#define vdup_odd(v)  ((rnd)_mm_unpackhi_pd(v,v))
		#define vxchg_even_odd(v) vxchg(v)
		#define vblend_even_odd(a,b) ((rnd)_mm_shuffle_pd(a,b,2))
		#define vneg_even_xor_cte (*(const v2d*)_neg0)
		#define vxor(v,x) vxor2(v,x)
		#define reduce_add(a) ( _mm_cvtsd_f64(a) + _mm_cvtsd_f64(_mm_unpackhi_pd(a,a)) )
		#define vinterleave(a,b) {  rnd x = _mm_unpacklo_pd(a,b);	b = _mm_unpackhi_pd(a, b);	a = x; }
		#define vinterleave_reverse(a,b)	{  rnd x = _mm_unpacklo_pd(a,b);	a = _mm_unpackhi_pd(a, b);	b = x; }

		inline static void S2D_CSTORE2_4MAGIC(double* mem, long idx, rnd nr, rnd sr, rnd ni, rnd si) {
			((s2d*)mem)[idx*4]   = _mm_unpacklo_pd(nr, ni);	// aa = north_ri[0]
			((s2d*)mem)[idx*4+1] = _mm_unpacklo_pd(sr, si);	// cc = south_ri[0]
			((s2d*)mem)[idx*4+2] = _mm_unpackhi_pd(nr, ni);	// bb = north_ri[1]
			((s2d*)mem)[idx*4+3] = _mm_unpackhi_pd(sr, si);	// dd = south_ri[1]
		}
	#endif
/*	#ifdef __SSE3__
		//#define addi(a,b) _mm_addsub_pd(a, _mm_shuffle_pd((b),(b),1))		// a + I*b
		//#define subadd(a,b) _mm_addsub_pd(a, b)		// [al-bl, ah+bh]
		//#define CMUL(a,b) _mm_addsub_pd(_mm_shuffle_pd(a,a,0)*b, _mm_shuffle_pd(a,a,3)*_mm_shuffle_pd(b,b,1))
	#else
		//#define addi(a,b) ( (a) + (_mm_shuffle_pd((b),(b),1) * _mm_set_pd(1.0, -1.0)) )		// a + I*b		[note: _mm_set_pd(imag, real)) ]
		//#define subadd(a,b) ( (a) + (b) * _mm_set_pd(1.0, -1.0) )		// [al-bl, ah+bh]
	#endif	*/
	inline static v2d IxKxZ(double k, v2d z) {		// I*k*z,  allowing to use FMA.
		return (v2d) _mm_setr_pd(-k,k) * vxchg(z);
	}

	// vset(lo, hi) takes two doubles and pack them in a vector
	#define vset(lo, hi) _mm_set_pd(hi, lo)
	// vdup(x) takes a double and duplicate it to a vector of 2 doubles.
	#define vdup(x) ((s2d)_mm_set1_pd(x))
	// vxchg(a) exchange hi and lo component of vector a
	#define vlo_to_cplx(a) _mm_unpacklo_pd(a, vdup(0.0))
	#define vhi_to_cplx(a) _mm_unpackhi_pd(a, vdup(0.0))
	#define vcplx_real(a) vlo_to_dbl(a)
	#define vcplx_imag(a) vhi_to_dbl(a)
	#ifdef __clang__
		// allow to compile with clang (llvm)
		#define vlo(a) (a)[0]
		#define vlo_to_dbl(a) (a)[0]
		#define vhi_to_dbl(a) (a)[1]
	#else
		// gcc extensions
		#ifdef __AVX512F__
			#define vlo(a) _mm_cvtsd_f64(_mm512_castpd512_pd128(a))
		#elif defined __AVX__
			#define vlo(a) _mm_cvtsd_f64(_mm256_castpd256_pd128(a))
		#else
			#define vlo(a) _mm_cvtsd_f64(a)
		#endif
		#define vlo_to_dbl(a) _mm_cvtsd_f64(a)
		#define vhi_to_dbl(a) _mm_cvtsd_f64(_mm_unpackhi_pd(a,a))
	#endif
	#ifdef __AVX512F__
		#define v2d_lo(a) (v2d)_mm512_castpd512_pd128(a)
	#elif defined __AVX__
		#define v2d_lo(a) (v2d)_mm256_castpd256_pd128(a)
	#else
		#define v2d_lo(a) (a)
	#endif
#endif



#ifndef _GCC_VEC_
	#define MIN_ALIGNMENT 16
	#define VSIZE 1
	#define VSIZE2 1
	#define _SIMD_NAME_ "scalar"
	typedef double s2d;
	#ifndef __cplusplus
	#include <complex.h>
	typedef complex double v2d;
	#else
	#include <complex>
	typedef std::complex<double> v2d;
	#define I v2d(0.,1.)
	#endif
	typedef double rnd;
	#define vread(mem, idx) ((double*)(mem))[idx]
	#define vstor(mem, idx, v) ((double*)(mem))[idx] = v;
	#define vread2(mem, idx) ((v2d*)(mem))[idx]
	#define vstor2(mem, idx, v) ((v2d*)(mem))[idx] = v;
	#define reduce_add(a) (a)
	#define v2d_reduce(a,b) ((a) +I*(b))	
	#define vlo(a) (a)
	#define vall(x) (x)
	#define vdup(x) (x)
	#define vxchg(x) (x)
	//#define addi(a,b) ((a) + I*(b))
	#define vlo_to_dbl(a) (a)
	#define vhi_to_dbl(a) (a)
	#define vcplx_real(a) creal(a)
	#define vcplx_imag(a) cimag(a)
	inline static v2d IxKxZ(double k, v2d z) {		// I*k*z
		return I*k*z;
	}
	#define vinterleave(a,b)
	#define vinterleave_reverse(a,b)
	#define vreverse(a) (a)
	//#define vreverse_pairs(a) (a)
	#define vxchg_even_odd(a) (a)

	inline static void S2D_CSTORE2_4MAGIC(double* mem, long idx, double nr, double sr, double ni, double si) {
		((v2d*)mem)[2*idx]   = (nr) + I*(ni);
		((v2d*)mem)[2*idx+1] = (sr) + I*(si);
	}
#endif

/// Aligned malloc on 64 bytes (cache-line) that fits any vector size up to 512 bits.
#if _GCC_VEC_ && __SSE2__
	#define VMALLOC(s)	_mm_malloc(s, 64)
	#define VFREE(s)	_mm_free(s)
#else
	inline static void* VMALLOC(size_t s) {
		void* ptr = 0;		// return value will be zero on failure.
		posix_memalign(&ptr, 64, s);
		return ptr;
	}
	#define VFREE(s)	free(s)
#endif

#define SSE __attribute__((aligned (MIN_ALIGNMENT)))

/// align pointer on MIN_ALIGNMENT (must be a power of 2)
#define PTR_ALIGN(p) ((((size_t)(p)) + (MIN_ALIGNMENT-1)) & (~((size_t)(MIN_ALIGNMENT-1))))

#ifdef __GNUC__
#define LIKELY(x)    (__builtin_expect (!!(x), 1))
#define UNLIKELY(x)  (__builtin_expect (!!(x), 0))
#else
#define LIKELY(x)    (x)
#define UNLIKELY(x)  (x)
#endif
