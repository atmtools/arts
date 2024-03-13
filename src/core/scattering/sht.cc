#include "sht.h"

namespace scattering {
namespace sht {

/// Deleter function for use with smart pointers.
struct FFTWDeleter {
  template <typename T>
  void operator()(T *t) {
    if (t) {
      fftw_free(t);
      t = nullptr;
    }
  }
};

template <typename Numeric>
FFTWArray<Numeric>::FFTWArray(Index n) : ptr_(nullptr) {
  if (n > 0) {
    ptr_ = std::shared_ptr<Numeric>(
        reinterpret_cast<Numeric *>(fftw_malloc(n * sizeof(Numeric))),
        FFTWDeleter());
  }
}

shtns_cfg ShtnsHandle::get(Index l_max, Index m_max, Index n_lon, Index n_lat) {
  std::array<Index, 4> config = {l_max, m_max, n_lon, n_lat};
  if (config == current_config_) {
    return shtns_;
  } else {
    shtns_reset();
    shtns_ = shtns_init(sht_reg_fast,
                        static_cast<int>(l_max),
                        static_cast<int>(m_max),
                        1,
                        static_cast<int>(n_lat),
                        static_cast<int>(n_lon));
    current_config_ = config;
  }
  return shtns_;
}

shtns_cfg ShtnsHandle::shtns_ = nullptr;
std::array<Index, 4> ShtnsHandle::current_config_ = {-1, -1, -1, -1};

////////////////////////////////////////////////////////////////////////////////
// SHT
////////////////////////////////////////////////////////////////////////////////

SHT::SHT(Index l_max, Index m_max, Index n_lon, Index n_lat)
    : l_max_(l_max), m_max_(m_max), n_lon_(n_lon), n_lat_(n_lat) {
  if (l_max == 0) {
    is_trivial_ = true;
    n_spectral_coeffs_ = 1;
    n_spectral_coeffs_cmplx_ = 1;
  } else {
    is_trivial_ = false;
    shtns_verbose(1);
    shtns_use_threads(0);
    n_spectral_coeffs_ = calc_n_spectral_coeffs(l_max, m_max);
    n_spectral_coeffs_cmplx_ = calc_n_spectral_coeffs_cmplx(l_max, m_max);
    spectral_coeffs_ = sht::FFTWArray<std::complex<double>>(n_spectral_coeffs_);
    spectral_coeffs_cmplx_ =
        sht::FFTWArray<std::complex<double>>(n_spectral_coeffs_cmplx_);
    spatial_coeffs_ = sht::FFTWArray<double>(n_lon * n_lat);
    spatial_coeffs_cmplx_ = sht::FFTWArray<std::complex<double>>(n_lon * n_lat);
    za_grid_ = std::make_shared<LatGrid>(get_latitude_grid());
    aa_grid_ = std::make_shared<Vector>(get_longitude_grid());
  }
}

SHT::SHT(Index l_max, Index m_max)
    : SHT(l_max,
          m_max,
          (m_max > 0) ? 2 * m_max + 2 : 1,
          (l_max > 0) ? 2 * l_max + 2 : 1) {}

SHT::SHT(Index l_max) : SHT(l_max, l_max) {}

Vector SHT::get_colatitude_grid() {
  if (is_trivial_) {
    return Vector(1, 0.0);
  }
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_lon_, n_lat_);
  Vector result(n_lat_);
  std::copy_n(shtns->ct, n_lat_, result.begin());
  return result;
}

ArrayOfIndex SHT::get_l_indices() {
  if (is_trivial_) {
    return {0};
  }
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_lon_, n_lat_);
  ArrayOfIndex result(n_spectral_coeffs_);
  for (Index i = 0; i < n_spectral_coeffs_; ++i) {
    result[i] = shtns->li[i];
  }
  return result;
}

ArrayOfIndex SHT::get_m_indices() {
  if (is_trivial_) {
    return {0};
  }
  auto shtns = ShtnsHandle::get(l_max_, m_max_, n_lon_, n_lat_);
  ArrayOfIndex result(n_spectral_coeffs_);
  for (Index i = 0; i < n_spectral_coeffs_; ++i) {
    result[i] = shtns->mi[i];
  }
  return result;
}

////////////////////////////////////////////////////////////////////////////////
// SHTProvider
////////////////////////////////////////////////////////////////////////////////

SHTProvider provider{};

}  // namespace sht
}  // namespace scattering
