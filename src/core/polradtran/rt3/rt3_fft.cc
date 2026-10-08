#include "rt3_fft.h"

#include <arts_constants.h>
#include <debug.h>

#include <cmath>
#include <utility>

namespace polradtran::rt3 {
void makephase(VectorView phase) {
  using Constant::pi;

  // Only for a power of two do the two halves fit in 4 nmax values; for
  // another nmax MAKEPHASE writes past them (FFT1DR's MN is a power of two)
  const Index size = static_cast<Index>(phase.size());
  const Index nmax = size / 4;
  ARTS_USER_ERROR_IF(nmax < 1 or (nmax & (nmax - 1)) != 0 or size != 4 * nmax,
                     "MAKEPHASE needs phase [4 nmax] with nmax a power of two, got {} values",
                     size);

  // exp(+i pi k / n) from 0, exp(-i pi k / n) from 2 nmax
  for (const Numeric p : {pi, -pi}) {
    Index j = p > 0.0 ? 0 : 2 * nmax;
    Index n = 1;
    do {
      const Numeric f = p / static_cast<Numeric>(n);
      for (Index i = 0; i < n; i++) {
        phase[j]      = std::cos(f * static_cast<Numeric>(i));
        phase[j + 1]  = std::sin(f * static_cast<Numeric>(i));
        j            += 2;
      }
      n *= 2;
    } while (n < nmax);
  }
}

void fftc(VectorView data, ConstVectorView phase) {
  const Index size = static_cast<Index>(data.size());
  const Index n    = size / 2;
  ARTS_USER_ERROR_IF(n < 1 or size != 2 * n or (n & (n - 1)) != 0,
                     "FFTC transforms a power of two of complex values (2 n reals), got {} values",
                     size);
  ARTS_USER_ERROR_IF(static_cast<Index>(phase.size()) < 2 * n - 2,
                     "FFTC of {} complex values needs a phase table of at least {} values, got {}",
                     n,
                     2 * n - 2,
                     phase.size());

  if (n <= 1) return;

  // Reorder the complex values to bit-reversed order
  Index irev = 0;
  for (Index i = 0; i < n; i++) {
    if (i > irev) {
      std::swap(data[2 * i], data[2 * irev]);
      std::swap(data[2 * i + 1], data[2 * irev + 1]);
    }
    Index m = n;
    do {
      m /= 2;
      if (irev < m) break;
      irev -= m;
    } while (m > 1);
    irev += m;
  }

  // Combine the transforms of length power into ones of length 2 power.
  // m0 and m1 are the real parts of the two values of a butterfly, iph the
  // cosine of its twiddle factor (PHASE(IPH-1)), all 0-based.
  Index jmax  = n;
  Index power = 1;
  do {
    jmax     /= 2;
    Index m0  = 0;
    Index m1  = power * 2;
    for (Index j = 0; j < jmax; j++) {
      Index iph = 2 * power - 2;
      for (Index k = 0; k < power; k++) {
        const Numeric phr   = phase[iph];
        const Numeric phi   = phase[iph + 1];
        iph                += 2;
        const Numeric tmpr  = phr * data[m1] - phi * data[m1 + 1];
        const Numeric tmpi  = phi * data[m1] + phr * data[m1 + 1];
        data[m1]            = data[m0] - tmpr;
        data[m1 + 1]        = data[m0 + 1] - tmpi;
        data[m0]            = data[m0] + tmpr;
        data[m0 + 1]        = data[m0 + 1] + tmpi;
        m0                 += 2;
        m1                 += 2;
      }
      m0 += power * 2;
      m1 += power * 2;
    }
    power *= 2;
  } while (jmax > 1);
}

void fixreal(VectorView data, Vector2& nyquist, fft_direction isign, ConstVectorView phase) {
  const Index size = static_cast<Index>(data.size());
  const Index n    = size / 2;
  ARTS_USER_ERROR_IF(n < 1 or size != 2 * n or (n & (n - 1)) != 0,
                     "FIXREAL works on a power of two of complex values (2 n reals), got {} values",
                     size);
  ARTS_USER_ERROR_IF(n > 1 and static_cast<Index>(phase.size()) < 3 * n,
                     "FIXREAL of {} complex values needs a phase table of at least {} values, got {}",
                     n,
                     3 * n,
                     phase.size());

  // The real parts of the pairs k and n - k, from k = 1, and the cosine of
  // the twiddle factor of k (PHASE(2*N+1) for k = 1), all 0-based.  The
  // middle pair, k = n / 2, has m = mc.
  Index iph = 2 * n;
  Index m   = 2;
  Index mc  = 2 * n - 2;
  if (isign == fft_direction::forward) {
    nyquist[0] = data[0] - data[1];
    nyquist[1] = 0.0;
    data[0]    = data[0] + data[1];
    data[1]    = 0.0;
    for (Index i = 2; i <= n / 2 + 1; i++) {
      const Numeric phr    = phase[iph];
      const Numeric phi    = phase[iph + 1];
      iph                 += 2;
      const Numeric tmp0r  = data[m] + data[mc];
      const Numeric tmp0i  = data[m + 1] - data[mc + 1];
      const Numeric dr     = data[m] - data[mc];
      const Numeric si     = data[m + 1] + data[mc + 1];
      const Numeric tmp1r  = -phi * dr - phr * si;
      const Numeric tmp1i  = phr * dr - phi * si;
      data[m]              = 0.5 * (tmp0r - tmp1r);
      data[m + 1]          = 0.5 * (tmp0i - tmp1i);
      data[mc]             = 0.5 * (tmp0r + tmp1r);
      data[mc + 1]         = -0.5 * (tmp0i + tmp1i);
      m                   += 2;
      mc                  -= 2;
    }
  } else {
    data[1] = data[0] - nyquist[0];
    data[0] = data[0] + nyquist[0];
    for (Index i = 2; i <= n / 2 + 1; i++) {
      const Numeric phr    = phase[iph];
      const Numeric phi    = phase[iph + 1];
      iph                 += 2;
      const Numeric tmp0r  = data[m] + data[mc];
      const Numeric tmp0i  = data[m + 1] - data[mc + 1];
      const Numeric dr     = data[m] - data[mc];
      const Numeric si     = data[m + 1] + data[mc + 1];
      const Numeric tmp1r  = phi * dr + phr * si;
      const Numeric tmp1i  = -phr * dr + phi * si;
      data[m]              = tmp0r - tmp1r;
      data[m + 1]          = tmp0i - tmp1i;
      data[mc]             = tmp0r + tmp1r;
      data[mc + 1]         = -(tmp0i + tmp1i);
      m                   += 2;
      mc                  -= 2;
    }
  }
}

void fft1dr(VectorView data, fft_direction isign, fft_workdata& work) {
  //! FFT1DR's MAXN, the size of its phase table
  constexpr Index maxn = 512;

  const Index n = data.size();
  ARTS_USER_ERROR_IF(n < 2 or (n & (n - 1)) != 0, "FFT1DR transforms a power of two of at least 2 values, got {}", n);
  ARTS_USER_ERROR_IF(n > maxn, "Phase array too small: FFT1DR transforms at most {} values, got {}", maxn, n);

  if (work.mn < n) {
    work.mn = n;
    work.phase.resize(4 * work.mn);
    makephase(work.phase);
  }

  Vector2 nyquist{};
  if (isign == fft_direction::forward) {
    // Forward transform:  real to complex-conjugate
    const auto phase = work.phase[Range{0, 2 * work.mn}];
    fftc(data, phase);
    fixreal(data, nyquist, fft_direction::forward, phase);
    data[1] = nyquist[0];
  } else {
    // Inverse transform:  complex-conjugate to real
    const auto phase = work.phase[Range{2 * work.mn, 2 * work.mn}];
    nyquist[0]       = data[1];
    fixreal(data, nyquist, fft_direction::inverse, phase);
    fftc(data, phase);
  }
}
}  // namespace polradtran::rt3
