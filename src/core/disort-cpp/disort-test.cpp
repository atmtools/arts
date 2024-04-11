#include <disort.h>
#include <iomanip>
#include <iostream>

void test1a() {
  const Vector tau_arr{0.03125};
  const Vector omega_arr{0.2};
  const Index NQuad = 16;
  Matrix Leg_coeffs_all(1, 17, 0);
  Leg_coeffs_all(0, 0) = 1;

  const Numeric mu0 = 0.1;
  const Numeric I0 = Constant::pi / mu0;
  const Numeric phi0 = 0;

  // Optional (unused)
  const Index NLeg = NQuad;
  const Index NFourier = NQuad;
  const Matrix b_pos(1, 1, 0);
  const Matrix b_neg(1, 1, 0);
  const bool only_flux = false;
  const Vector f_arr{0};
  const bool NT_cor = false;
  const std::vector<disort::BDRF> BDRF_Fourier_modes{};
  const Matrix s_poly_coeffs(0, 0);

  const disort::main_data dis(NQuad,
                              tau_arr,
                              omega_arr,
                              Leg_coeffs_all,
                              b_pos,
                              b_neg,
                              f_arr,
                              s_poly_coeffs,
                              BDRF_Fourier_modes,
                              mu0,
                              I0,
                              phi0);

const Numeric tau=0.01;
const Numeric phi=1;

disort::u0_data u0data;
dis.u0(u0data, tau);
std::cout << u0data.u0 << '\n';

disort::u_data udata;
dis.u(udata, tau, phi, false);
std::cout << udata.intensities << '\n';

disort::flux_data fluxdata;
std::cout << dis.flux_up(fluxdata, tau) << '\n';
const auto [dif, dir] = dis.flux_down(fluxdata, tau);
std::cout << dif << ' ' << dir << '\n';
}

int main() { test1a(); }
