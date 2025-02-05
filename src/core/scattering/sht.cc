#include "sht.h"

#include <mutex>

#ifndef ARTS_NO_SHTNS
#include <fftw3.h>
#include <shtns.h>
#endif

namespace scattering {
namespace sht {
static std::mutex shtns_mutex{};

Vector SHT::get_azimuth_angle_grid(Index n_aa, bool radians) {
  if (n_aa == 1) {
    Vector result(1);
    result = 0.0;
    return result;
  }
  Vector result(n_aa);
  double dx = 360.0 / static_cast<double>(n_aa);
  for (Index i = 0; i < n_aa; ++i) {
    result[i] = dx * static_cast<double>(i);
  }
  if (radians) {
    result *= Conversion::deg2rad(1.0);
  }
  return result;
}

FejerGrid SHT::get_zenith_angle_grid(bool radians) const {
  return get_zenith_angle_grid(n_za_, radians);
}

std::shared_ptr<const ZenithAngleGrid> SHT::get_za_grid_ptr() const {
  return za_grid_;
}

Index SHT::calc_n_spectral_coeffs(Index l_max, Index m_max) {
  return (l_max + 1) * (m_max + 1) - (m_max * (m_max + 1)) / 2;
}

Index SHT::calc_n_spectral_coeffs_cmplx(Index l_max, Index m_max) {
  return (2 * m_max + 1) * (l_max + 1) - m_max * (m_max + 1);
}

Index SHT::calc_l_max(Index n_spectral_coeffs) {
  return static_cast<Index>(
      sqrt(2.0 * static_cast<double>(n_spectral_coeffs) + 0.25) - 1.5);
}

std::array<Index, 4> SHT::get_config_lonlat(Index n_aa, Index n_za) {
  if (n_aa > 1) {
    n_aa -= n_aa % 2;
  }
  n_za -= n_za % 2;

  Index l_max = (n_za > 2) ? (n_za / 2) - 1 : 0;
  Index m_max = (n_aa > 2) ? (n_aa / 2) - 1 : 0;
  m_max       = std::min(l_max, m_max);
  return {l_max, m_max, n_aa, n_za};
}

std::array<Index, 4> SHT::get_config_lm(Index l_max, Index m_max) {
  return {l_max,
          m_max,
          (m_max > 0) ? 2 * m_max + 2 : 1,
          (l_max > 0) ? 2 * l_max + 2 : 1};
}

std::ostream &SHT::serialize(std::ostream &output) const {
  output.write(reinterpret_cast<const char *>(&l_max_), sizeof(Index));
  output.write(reinterpret_cast<const char *>(&m_max_), sizeof(Index));
  output.write(reinterpret_cast<const char *>(&n_aa_), sizeof(Index));
  output.write(reinterpret_cast<const char *>(&n_za_), sizeof(Index));
  return output;
}

SHT SHT::deserialize(std::istream &input) {
  Index l_max, m_max, n_aa, n_za;
  input.read(reinterpret_cast<char *>(&l_max), sizeof(Index));
  input.read(reinterpret_cast<char *>(&m_max), sizeof(Index));
  input.read(reinterpret_cast<char *>(&n_aa), sizeof(Index));
  input.read(reinterpret_cast<char *>(&n_za), sizeof(Index));
  return SHT(l_max, m_max, n_aa, n_za);
}

Vector SHT::get_azimuth_angle_grid(bool radians) {
  return SHT::get_azimuth_angle_grid(n_aa_, radians);
}

Index SHT::get_coeff_index(Index l, Index m) {
  // l parameter varies fastest.
  m = std::abs(m);
  return m * (l_max_ + 1) - (m * (m - 1)) / 2 + l - m;
}

void SHT::set_spectral_coeffs(
    const matpack::strided_view_t<const Complex, 1> &view) const {
  // Input size must match number of spectral coefficients of SHT.
  ARTS_ASSERT(view.size() == static_cast<Size>(n_spectral_coeffs_));
  Index index = 0;
  for (auto &x : view) {
    spectral_coeffs_[index] = x;
    ++index;
  }
}

void SHT::set_spectral_coeffs_cmplx(
    const matpack::strided_view_t<const Complex, 1> &view) const {
  // Input size must match number of spectral coefficients of SHT.
  ARTS_ASSERT(view.size() == static_cast<Size>(n_spectral_coeffs_cmplx_));
  Index index = 0;
  for (auto &x : view) {
    spectral_coeffs_cmplx_[index] = x;
    ++index;
  }
}

StridedConstMatrixView SHT::get_spatial_coeffs() const {
  return StridedConstMatrixView{
      matpack::mdview_t<Numeric, 2>(spatial_coeffs_, std::array{n_aa_, n_za_})};
}

StridedConstComplexMatrixView SHT::get_spatial_coeffs_cmplx() const {
  return StridedConstComplexMatrixView{matpack::mdview_t<Complex, 2>(
      spatial_coeffs_cmplx_, std::array{n_aa_, n_za_})};
}

StridedConstComplexVectorView SHT::get_spectral_coeffs_cmplx() const {
  return StridedConstComplexVectorView{matpack::mdview_t<Complex, 1>(
      spectral_coeffs_cmplx_, std::array{n_spectral_coeffs_cmplx_})};
}

Index SHT::get_n_zenith_angles() const { return n_za_; }

Index SHT::get_n_azimuth_angles() const { return n_aa_; }

Index SHT::get_n_spectral_coeffs() const { return n_spectral_coeffs_; }

Index SHT::get_n_spectral_coeffs_cmplx() const {
  return n_spectral_coeffs_cmplx_;
}

Index SHT::get_l_max() const { return l_max_; }

Index SHT::get_m_max() const { return m_max_; }

StridedConstComplexVectorView SHT::get_spectral_coeffs() const {
  return StridedConstComplexVectorView{matpack::mdview_t<Complex, 1>(
      spectral_coeffs_, std::array{n_spectral_coeffs_})};
}

ComplexVector SHT::transform(const StridedConstMatrixView &view
                             [[maybe_unused]]) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (is_trivial_) {
    ComplexVector result(1);
    result[0] = view[0, 0];
    return result;
  }
  set_spatial_coeffs(view);
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_aa_, n_za_);

  shtns_mutex.lock();
  spat_to_SH(shtns, spatial_coeffs_, spectral_coeffs_);
  shtns_mutex.unlock();

  return static_cast<ComplexVector>(get_spectral_coeffs());
#endif
}

ComplexVector SHT::transform_cmplx(const StridedConstComplexMatrixView &view
                                   [[maybe_unused]]) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (is_trivial_) {
    ComplexVector result(1);
    result[0] = view[0, 0];
    return result;
  }
  set_spatial_coeffs(view);
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_aa_, n_za_);

  shtns_mutex.lock();
  spat_cplx_to_SH(shtns, spatial_coeffs_cmplx_, spectral_coeffs_cmplx_);
  shtns_mutex.unlock();

  return static_cast<ComplexVector>(get_spectral_coeffs_cmplx());
#endif
}

Matrix SHT::synthesize(const StridedConstComplexVectorView &view
                       [[maybe_unused]]) {
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

  shtns_mutex.lock();
  SH_to_spat(shtns, spectral_coeffs_, spatial_coeffs_);
  shtns_mutex.unlock();

  return static_cast<Matrix>(get_spatial_coeffs());
#endif
}

ComplexMatrix SHT::synthesize_cmplx(const StridedConstComplexVectorView &view
                                    [[maybe_unused]]) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (is_trivial_) {
    ComplexMatrix result(1, 1);
    result[0, 0] = view[0];
    return result;
  }
  set_spectral_coeffs_cmplx(view);
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_aa_, n_za_);

  shtns_mutex.lock();
  SH_to_spat_cplx(shtns, spectral_coeffs_cmplx_, spatial_coeffs_cmplx_);
  shtns_mutex.unlock();

  return static_cast<ComplexMatrix>(get_spatial_coeffs_cmplx());
#endif
}

Numeric SHT::evaluate(const StridedConstComplexVectorView &view
                      [[maybe_unused]],
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

  const std::scoped_lock lock(shtns_mutex);
  return SH_to_point(shtns, spectral_coeffs_, std::cos(theta), phi);
#endif
}

Vector SHT::evaluate(const StridedComplexVectorView &view [[maybe_unused]],
                     const StridedMatrixView &points [[maybe_unused]]) {
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

  shtns_mutex.lock();
  for (auto i = 0; i < n_points; ++i) {
    result[i] = SH_to_point(
        shtns, spectral_coeffs_, std::cos(points[i, 1]), points[i, 0]);
  }
  shtns_mutex.lock();

  return result;
#endif
}

Vector SHT::evaluate(const StridedConstComplexVectorView &view [[maybe_unused]],
                     const Vector &thetas [[maybe_unused]]) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
#else
  if (is_trivial_) {
    Vector results(view.size());
    for (Size i = 0; i < view.size(); ++i) {
      results[i] = view[i].real();
    }
    return results;
  }
  ARTS_ASSERT(m_max_ == 0);
  set_spectral_coeffs(view);
  Size n_points = thetas.size();
  Vector result(n_points);
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_aa_, n_za_);

  shtns_mutex.lock();
  for (Size i = 0; i < n_points; ++i) {
    result[i] = SH_to_point(shtns, spectral_coeffs_, cos(thetas[i]), 0.0);
  }
  shtns_mutex.unlock();

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
      shtns_mutex.lock();
      fftw_free(t);
      shtns_mutex.unlock();
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
    shtns_mutex.lock();
    ptr_ = std::shared_ptr<Numeric>(
        reinterpret_cast<Numeric *>(fftw_malloc(n * sizeof(Numeric))),
        FFTWDeleter());
    shtns_mutex.unlock();
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
    shtns_mutex.lock();
    shtns_reset();
    shtns_ = shtns_init(sht_reg_fast,
                        static_cast<int>(l_max),
                        static_cast<int>(m_max),
                        1,
                        static_cast<int>(n_za),
                        static_cast<int>(n_aa));
    shtns_mutex.unlock();

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

    shtns_mutex.lock();
    shtns_verbose(1);
    shtns_use_threads(0);
    shtns_mutex.unlock();

    n_spectral_coeffs_       = calc_n_spectral_coeffs(l_max, m_max);
    n_spectral_coeffs_cmplx_ = calc_n_spectral_coeffs_cmplx(l_max, m_max);
    spectral_coeffs_ =
        sht::FFTWArray<std::complex<double> >(n_spectral_coeffs_);
    spectral_coeffs_cmplx_ =
        sht::FFTWArray<std::complex<double> >(n_spectral_coeffs_cmplx_);
    spatial_coeffs_       = sht::FFTWArray<double>(n_aa * n_za);
    spatial_coeffs_cmplx_ = sht::FFTWArray<std::complex<double> >(n_aa * n_za);
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
    result.data *= Conversion::deg2rad(1.0);
  }
  return result;
};

////////////////////////////////////////////////////////////////////////////////
// SHTProvider
////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<SHT> SHTProvider::get_instance(SHTParams params
                                               [[maybe_unused]]) {
#ifdef ARTS_NO_SHTNS
  ARTS_USER_ERROR("Not compiled with SHTNS or FFTW support.");
  std::unreachable();
#else
  if (sht_instances_.count(params) == 0) {
    sht_instances_[params] =
        std::make_shared<SHT>(params[0], params[1], params[2], params[3]);
  }
  return sht_instances_[params];
#endif
}

std::shared_ptr<SHT> SHTProvider::get_instance(Index n_aa, Index n_za) {
  return get_instance(SHT::get_config_lonlat(n_aa, n_za));
}

std::shared_ptr<SHT> SHTProvider::get_instance_lm(Index l_max, Index m_max) {
  return get_instance(SHT::get_config_lm(l_max, m_max));
}

SHTProvider provider{};

}  // namespace sht
}  // namespace scattering
