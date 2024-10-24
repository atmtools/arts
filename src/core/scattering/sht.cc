#include "sht.h"

#ifndef ARTS_NO_SHTNS
#include <fftw3.h>
#include <shtns.h>
#endif

namespace scattering {
namespace sht {
ComplexVector SHT::transform(const ConstMatrixView &view [[maybe_unused]]) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (is_trivial_) {
    ComplexVector result(1);
    result[0] = view(0, 0);
    return result;
  }
  set_spatial_coeffs(view);
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_aa_, n_za_);
  spat_to_SH(shtns, spatial_coeffs_, spectral_coeffs_);
  return static_cast<ComplexVector>(get_spectral_coeffs());
#endif
}

ComplexVector SHT::transform_cmplx(const ConstComplexMatrixView &view
                                   [[maybe_unused]]) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (is_trivial_) {
    ComplexVector result(1);
    result[0] = view(0, 0);
    return result;
  }
  set_spatial_coeffs(view);
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_aa_, n_za_);
  spat_cplx_to_SH(shtns, spatial_coeffs_cmplx_, spectral_coeffs_cmplx_);
  return static_cast<ComplexVector>(get_spectral_coeffs_cmplx());
#endif
}

Matrix SHT::synthesize(const ConstComplexVectorView &view [[maybe_unused]]) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (is_trivial_) {
    Matrix result(1, 1);
    result = view[0].real();
    return result;
  }
  set_spectral_coeffs(view);
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_aa_, n_za_);
  SH_to_spat(shtns, spectral_coeffs_, spatial_coeffs_);
  return static_cast<Matrix>(get_spatial_coeffs());
#endif
}

ComplexMatrix SHT::synthesize_cmplx(const ConstComplexVectorView &view
                                    [[maybe_unused]]) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (is_trivial_) {
    ComplexMatrix result(1, 1);
    result(0, 0) = view[0];
    return result;
  }
  set_spectral_coeffs_cmplx(view);
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_aa_, n_za_);
  SH_to_spat_cplx(shtns, spectral_coeffs_cmplx_, spatial_coeffs_cmplx_);
  return static_cast<ComplexMatrix>(get_spatial_coeffs_cmplx());
#endif
}

Numeric SHT::evaluate(const ConstComplexVectorView &view [[maybe_unused]],
                      Numeric phi [[maybe_unused]],
                      Numeric theta [[maybe_unused]]) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (is_trivial_) {
    return view[0].real();
  }
  set_spectral_coeffs(view);
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_aa_, n_za_);
  return SH_to_point(shtns, spectral_coeffs_, cos(theta), phi);
#endif
}

Vector SHT::evaluate(const ComplexVectorView &view [[maybe_unused]],
                     const MatrixView &points [[maybe_unused]]) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (is_trivial_) {
    Vector results(points.nrows());
    results = view[0].real();
    return results;
  }
  set_spectral_coeffs(view);
  auto n_points = points.nrows();
  Vector result(n_points);
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_aa_, n_za_);
  for (auto i = 0; i < n_points; ++i) {
    result[i] =
        SH_to_point(shtns, spectral_coeffs_, cos(points(i, 1)), points(i, 0));
  }
  return result;
#endif
}

Vector SHT::evaluate(const ConstComplexVectorView &view [[maybe_unused]],
                     const Vector &thetas [[maybe_unused]]) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (is_trivial_) {
    Vector results(view.size());
    for (Index i = 0; i < view.size(); ++i) {
      results[i] = view[i].real();
    }
    return results;
  }
  ARTS_ASSERT(m_max_ == 0);
  set_spectral_coeffs(view);
  auto n_points = thetas.size();
  Vector result(n_points);
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_aa_, n_za_);
  for (auto i = 0; i < n_points; ++i) {
    result[i] = SH_to_point(shtns, spectral_coeffs_, cos(thetas[i]), 0.0);
  }
  return result;
#endif
}

/// Deleter function for use with smart pointers.
struct FFTWDeleter {
  template <typename T>
  void operator()(T *t [[maybe_unused]]) {
#ifdef ARTS_NO_SHTNS
    ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
    if (t) {
      fftw_free(t);
      t = nullptr;
    }
#endif
  }
};

template <typename Numeric>
FFTWArray<Numeric>::FFTWArray(Index n [[maybe_unused]]) : ptr_(nullptr) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (n > 0) {
    ptr_ = std::shared_ptr<Numeric>(
        reinterpret_cast<Numeric *>(fftw_malloc(n * sizeof(Numeric))),
        FFTWDeleter());
  }
#endif
}

shtns_cfg ShtnsHandle::get(Index l_max [[maybe_unused]],
                           Index m_max [[maybe_unused]],
                           Index n_aa [[maybe_unused]],
                           Index n_za [[maybe_unused]]) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  std::array<Index, 4> config = {l_max, m_max, n_aa, n_za};
  if (config == current_config_) {
    return shtns_;
  } else {
    shtns_reset();
    shtns_          = shtns_init(sht_reg_fast,
                        static_cast<int>(l_max),
                        static_cast<int>(m_max),
                        1,
                        static_cast<int>(n_za),
                        static_cast<int>(n_aa));
    current_config_ = config;
  }
  return shtns_;
#endif
}

shtns_cfg ShtnsHandle::shtns_                     = nullptr;
std::array<Index, 4> ShtnsHandle::current_config_ = {-1, -1, -1, -1};

////////////////////////////////////////////////////////////////////////////////
// SHT
////////////////////////////////////////////////////////////////////////////////

SHT::SHT(Index l_max, Index m_max, Index n_aa, Index n_za)
    : l_max_(l_max), m_max_(m_max), n_aa_(n_aa), n_za_(n_za) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (l_max == 0) {
    is_trivial_              = true;
    n_spectral_coeffs_       = 1;
    n_spectral_coeffs_cmplx_ = 1;
  } else {
    is_trivial_ = false;
    shtns_verbose(1);
    shtns_use_threads(0);
    n_spectral_coeffs_       = calc_n_spectral_coeffs(l_max, m_max);
    n_spectral_coeffs_cmplx_ = calc_n_spectral_coeffs_cmplx(l_max, m_max);
    spectral_coeffs_ =
        sht::FFTWArray<std::complex<double> >(n_spectral_coeffs_);
    spectral_coeffs_cmplx_ =
        sht::FFTWArray<std::complex<double> >(n_spectral_coeffs_cmplx_);
    spatial_coeffs_ = sht::FFTWArray<double>(n_aa * n_za);
    spatial_coeffs_cmplx_ =
        sht::FFTWArray<std::complex<double> >(n_aa * n_za);
    za_grid_ = std::make_shared<ZenithAngleGrid>(get_zenith_angle_grid());
    aa_grid_ = std::make_shared<Vector>(get_azimuth_angle_grid());
  }
#endif
}

SHT::SHT(Index l_max, Index m_max)
    : SHT(l_max,
          m_max,
          (m_max > 0) ? 2 * m_max + 2 : 1,
          (l_max > 0) ? 2 * l_max + 2 : 1) {}

SHT::SHT(Index l_max) : SHT(l_max, l_max) {}

Vector SHT::get_cos_za_grid() {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (is_trivial_) {
    return Vector(1, 0.0);
  }
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_aa_, n_za_);
  Vector result(n_za_);
  std::copy_n(shtns->ct, n_za_, result.begin());
  return result;
#endif
}

ArrayOfIndex SHT::get_l_indices() {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (is_trivial_) {
    return {0};
  }
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_aa_, n_za_);
  ArrayOfIndex result(n_spectral_coeffs_);
  for (Index i = 0; i < n_spectral_coeffs_; ++i) {
    result[i] = shtns->li[i];
  }
  return result;
#endif
}

ArrayOfIndex SHT::get_m_indices() {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (is_trivial_) {
    return {0};
  }
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_aa_, n_za_);
  ArrayOfIndex result(n_spectral_coeffs_);
  for (Index i = 0; i < n_spectral_coeffs_; ++i) {
    result[i] = shtns->mi[i];
  }
  return result;
#endif
}

FejerGrid SHT::get_zenith_angle_grid(Index n_za, bool radians) {
  auto result = FejerGrid(n_za);
  if (radians) {
    result *= Conversion::deg2rad(1.0);
  }
  return result;
};

////////////////////////////////////////////////////////////////////////////////
// SHTProvider
////////////////////////////////////////////////////////////////////////////////

SHTProvider provider{};

}  // namespace sht
}  // namespace scattering
