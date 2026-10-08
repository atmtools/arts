#pragma once

#include <matpack.h>

/* Evans' real FFT of 3rdparty/polradtran/radscat3.f (FFT1DR, and the FFTC,
   FIXREAL and MAKEPHASE it calls), ported to C++.  It is RT3's FFT, and
   the default.  It calls no Fortran and keeps its state in fft_workdata
   only.

   It is in a file pair of its own because a build may use FFTW instead, as
   a compile-time option only (FFTW's license keeps it out of the default
   build).  The rest of RT3 uses only fft1dr, fft_direction and
   fft_workdata.  An FFTW fft1dr must give the packed format and scaling
   below (FFTW's r2c followed by a complex conjugate, and a conjugate
   followed by c2r) and may keep its plans in fft_workdata.
   cpp.fast.rt3-radtran-test checks the format against direct sums. */
namespace rt3 {
//! FFT1DR's ISIGN
enum class fft_direction {
  forward,  //!< +1: real to complex conjugate
  inverse,  //!< -1: complex conjugate to real
};

/** FFT1DR's SAVEd state, kept between calls: MN, the largest length so
 * far, and PHASE, MAKEPHASE's table of exp(+-i pi k / m) for every power of
 * two m < MN, the + half from 0 and the - half from 2 MN.  Default
 * constructed it is empty; fft1dr builds the table as it needs it.
 */
struct fft_workdata {
  Index  mn{0};
  Vector phase{};
};

/** MAKEPHASE: the phase table of FFT1DR for lengths up to nmax, a power of
 * two.  For each power of two n < nmax (n = 1 at least), the n pairs
 * cos(pi k / n), sin(pi k / n), k = 0 to n - 1, follow each other from
 * phase[0], and the same for -pi from phase[2 nmax].  Entries in between
 * are not written.  Throws for another nmax, for which the halves would
 * overlap and run past 4 nmax.
 *
 *   phase  [4 nmax]  PHASE(4*NMAX), output
 */
void makephase(VectorView phase);

/** FFTC: complex FFT in place of the n complex values in data (interleaved
 * real and imaginary parts), n a power of two, with the twiddle factors of
 * one half of MAKEPHASE's table (exp(+i pi k / m) for the + half, which is
 * the transform with exp(+2 pi i j k / n), exp(-i pi k / m) for the -
 * half).  Bit-reversal reordering, then radix-2 butterflies.  No
 * normalization.
 *
 *   data   [2 n]          DATA(2*N), input and output
 *   phase  [2 n - 2 or more]  one half of PHASE
 */
void fftc(VectorView data, ConstVectorView phase);

/** FIXREAL: the step between the complex FFT of n points (fftc) and the
 * real FFT of 2 n points that FFT1DR packs into them.  forward: from the
 * transform of the 2 n reals taken as n complex values to the packed real
 * transform, the Nyquist value going to nyquist (nyquist[1] = 0) and
 * data[1] set to 0; inverse: the reverse, taking the Nyquist value from
 * nyquist[0].  phase is the half of MAKEPHASE's table of the direction
 * (FIXREAL reads its block for n, from phase[2 n]).
 *
 *   data   [2 n]                   DATA(2*N), input and output, n a power of two
 *   phase  [3 n or more, n > 1]    one half of PHASE
 */
void fixreal(VectorView data, Vector2& nyquist, fft_direction isign, ConstVectorView phase);

/** FFT1DR: real 1D FFT in place of data, of n = data.size() values, a
 * power of two from 2 to 512 (512 is FFT1DR's MAXN).  No normalization.
 *
 * forward: the n reals x_j become X_k = sum_j x_j exp(+2 pi i j k / n),
 *          packed as data = [X_0, X_{n/2}, Re X_1, Im X_1, ...,
 *          Re X_{n/2-1}, Im X_{n/2-1}] (X_0 and X_{n/2} are real).
 * inverse: the packed X become
 *          x_j = X_0 + (-1)^j X_{n/2} + 2 sum_{k=1}^{n/2-1} Re(X_k exp(-2 pi i j k / n)),
 *          so that inverse(forward(x)) = n x.
 */
void fft1dr(VectorView data, fft_direction isign, fft_workdata& work);
}  // namespace rt3
