#pragma once

#include <arts_constants.h>
#include <disort.h>

#include <functional>

namespace disort::brdf {

using RawFunction = std::function<Numeric(Numeric outgoing_mu, Numeric incoming_mu, Numeric relative_azimuth)>;

struct Hapke {
  Numeric opposition_amplitude{1.0};
  Numeric opposition_width{0.06};
  Numeric single_scattering_albedo{0.6};

  Numeric operator()(Numeric outgoing_mu, Numeric incoming_mu, Numeric relative_azimuth) const;
};

struct CoxMunk {
  Numeric wind_speed{12.0};
  Numeric refractive_index{1.34};
  bool    shadowing{false};

  Numeric operator()(Numeric outgoing_mu, Numeric incoming_mu, Numeric relative_azimuth) const;
};

struct RPV {
  Numeric rho0{0.027};
  Numeric kappa{0.647};
  Numeric asymmetry{-0.169};
  Numeric hotspot{0.1};

  Numeric operator()(Numeric outgoing_mu, Numeric incoming_mu, Numeric relative_azimuth) const;
};

struct RossLi {
  Numeric isotropic{0.091};
  Numeric volumetric{0.02};
  Numeric geometric{0.01};
  Numeric hotspot_angle{1.5 * Constant::pi / 180.0};

  Numeric operator()(Numeric outgoing_mu, Numeric incoming_mu, Numeric relative_azimuth) const;
};

/** Convert an azimuth-resolved scalar BRDF into the cosine Fourier modes
 * expected by disort::main_data.
 *
 * The coefficient convention is
 *   rho_m = (2 - delta_m0) / pi * integral_0^pi rho(dphi) cos(m dphi) dphi.
 * This is the convention used by DISORT's SURFAC/DISOBRDF routines and by
 * disort::BDRF.  incoming_mu may be signed when the resulting callbacks are
 * evaluated; raw functions always receive its positive magnitude.
 */
std::vector<BDRF> fourier_modes(RawFunction raw, Index number_of_modes, Index azimuth_quadrature_points = 100);

}  // namespace disort::brdf
