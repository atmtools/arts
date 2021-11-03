/* Copyright (C) 2020
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/*!
 * @file   predefined_absorption_models.cc
 * @author Richard Larsson
 * @date   2020-01-29
 * 
 * @brief  Full absorption models of various kinds
 */

#include "predefined_absorption_models.h"

#include <Faddeeva/Faddeeva.hh>
#include <algorithm>

#include "debug.h"
#include "jacobian.h"
#include "lin_alg.h"
#include "linescaling.h"
#include "matpack.h"
#include "propagationmatrix.h"
#include "quantum.h"
#include "species.h"

namespace Absorption::PredefinedModel::Makarov2020etal {
constexpr Index num = 38;

constexpr Numeric sum_lines(Numeric f,
                            const std::array<Numeric, num>& c,
                            const std::array<Numeric, num>& ga,
                            const std::array<Numeric, num>& g,
                            const std::array<Numeric, num>& y,
                            const std::array<Numeric, num>& dv,
                            const std::array<Numeric, num>& f0) noexcept {
  using Constant::pow2;

  Numeric a = 0;
  for (Index i = 0; i < num; i++)
    a += c[i] * ((ga[i] * (1 + g[i]) + y[i] * (f - f0[i] - dv[i])) /
                     (pow2(ga[i]) + pow2(f - f0[i] - dv[i])) +
                 (ga[i] * (1 + g[i]) - y[i] * (f + f0[i] + dv[i])) /
                     (pow2(ga[i]) + pow2(f + f0[i] + dv[i])));
  return a;
}

void impl(PropagationMatrix& propmat_clearsky,
          const Vector& f_grid,
          const Numeric& p_pa,
          const Numeric& t,
          const Numeric& oxygen_vmr) noexcept {
  using Constant::pow2, Constant::pow3, Constant::log10_euler;
  using Conversion::hz2ghz, Conversion::pa2bar;

  // constexpr std::array ids = identifiers();
  std::array<Numeric, num> c{
      940.3,  543.4,  1503.0, 1442.1, 2103.4, 2090.7, 2379.9, 2438.0,
      2363.7, 2479.5, 2120.1, 2275.9, 1746.6, 1915.4, 1331.8, 1490.2,
      945.3,  1078.0, 627.1,  728.7,  389.7,  461.3,  227.3,  274.0,
      124.6,  153.0,  64.29,  80.40,  31.24,  39.80,  14.32,  18.56,
      6.193,  8.172,  2.529,  3.397,  0.975,  1.334};
  constexpr std::array<Numeric, num> a2{
      0.01,  0.014, 0.083, 0.083, 0.207, 0.207, 0.387, 0.386, 0.621, 0.621,
      0.910, 0.910, 1.255, 1.255, 1.654, 1.654, 2.109, 2.108, 2.618, 2.617,
      3.182, 3.181, 3.800, 3.800, 4.474, 4.473, 5.201, 5.200, 5.983, 5.982,
      6.819, 6.818, 7.709, 7.708, 8.653, 8.652, 9.651, 9.650};
  std::array<Numeric, num> ga{
      1.685, 1.703, 1.513, 1.495, 1.433, 1.408, 1.353, 1.353, 1.303, 1.319,
      1.262, 1.265, 1.238, 1.217, 1.207, 1.207, 1.137, 1.137, 1.101, 1.101,
      1.037, 1.038, 0.996, 0.996, 0.955, 0.955, 0.906, 0.906, 0.858, 0.858,
      0.811, 0.811, 0.764, 0.764, 0.717, 0.717, 0.669, 0.669};
  std::array<Numeric, num> y0{
      -0.041, 0.277,  -0.372, 0.559,  -0.573, 0.618,  -0.366, 0.278,
      -0.089, -0.021, 0.060,  -0.152, 0.216,  -0.293, 0.373,  -0.436,
      0.491,  -0.542, 0.571,  -0.613, 0.636,  -0.670, 0.690,  -0.718,
      0.740,  -0.763, 0.788,  -0.807, 0.834,  -0.849, 0.876,  -0.887,
      0.915,  -0.922, 0.950,  -0.955, 0.987,  -0.988};
  constexpr std::array<Numeric, num> y1{
      0.0,   0.124,  -0.002, 0.008,  0.045, -0.093, 0.264, -0.351,
      0.359, -0.416, 0.326,  -0.353, 0.484, -0.503, 0.579, -0.590,
      0.616, -0.619, 0.611,  -0.609, 0.574, -0.568, 0.574, -0.566,
      0.60,  -0.59,  0.63,   -0.62,  0.64,  -0.63,  0.65,  -0.64,
      0.65,  -0.64,  0.65,   -0.64,  0.64,  -0.62};
  std::array<Numeric, num> g0{
      -0.000695, -0.090, -0.103, -0.239, -0.172, -0.171, 0.028,  0.150,
      0.132,     0.170,  0.087,  0.069,  0.083,  0.067,  0.007,  0.016,
      -0.021,    -0.066, -0.095, -0.115, -0.118, -0.140, -0.173, -0.186,
      -0.217,    -0.227, -0.234, -0.242, -0.266, -0.272, -0.301, -0.304,
      -0.334,    -0.333, -0.361, -0.358, -0.348, -0.344};
  constexpr std::array<Numeric, num> g1{
      0.,     -0.045, 0.007,  0.033,  0.081,  0.162,  0.179,  0.225,
      0.054,  0.003,  0.0004, -0.047, -0.034, -0.071, -0.180, -0.210,
      -0.285, -0.323, -0.363, -0.380, -0.378, -0.387, -0.392, -0.394,
      -0.424, -0.422, -0.465, -0.46,  -0.51,  -0.50,  -0.55,  -0.54,
      -0.58,  -0.56,  -0.62,  -0.59,  -0.68,  -0.65};
  std::array<Numeric, num> dv0{
      -0.00028, 0.00597, -0.0195, 0.032,   -0.0475, 0.0541,  -0.0232, 0.0154,
      0.0007,   -0.0084, -0.0025, -0.0014, -0.0004, -0.0020, 0.005,   -0.0066,
      0.0072,   -0.008,  0.0064,  -0.0070, 0.0056,  -0.0060, 0.0047,  -0.0049,
      0.0040,   -0.0041, 0.0036,  -0.0037, 0.0033,  -0.0034, 0.0032,  -0.0032,
      0.0030,   -0.0030, 0.0028,  -0.0029, 0.0029,  -0.0029};
  constexpr std::array<Numeric, num> dv1{
      -0.00039, 0.009,   -0.012, 0.016,   -0.027, 0.029,   0.006,  -0.015,
      0.010,    -0.014,  -0.013, 0.013,   0.004,  -0.005,  0.010,  -0.010,
      0.010,    -0.011,  0.008,  -0.009,  0.003,  -0.003,  0.0009, -0.0009,
      0.0017,   -0.0016, 0.0024, -0.0023, 0.0024, -0.0024, 0.0024, -0.0020,
      0.0017,   -0.0016, 0.0013, -0.0012, 0.0005, -0.0004};
  constexpr std::array<Numeric, num> f0 = {
      118.750334, 56.264774, 62.486253, 58.446588, 60.306056, 59.590983,
      59.164204,  60.434778, 58.323877, 61.150562, 57.612486, 61.800158,
      56.968211,  62.411220, 56.363399, 62.997984, 55.783815, 63.568526,
      55.221384,  64.127775, 54.671180, 64.678910, 54.130025, 65.224078,
      53.595775,  65.764779, 53.066934, 66.302096, 52.542418, 66.836834,
      52.021429,  67.369601, 51.503360, 67.900868, 50.987745, 68.431006,
      50.474214,  68.960312};

  // Conversion factor to ARTS units
  constexpr Numeric conv = 0.1820 * 1e-7 / (2.0946 * log10_euler);

  {
    // Pressure and temperature adaptation constants
    constexpr Numeric x = 0.754;
    const Numeric p = pa2bar(p_pa);
    const Numeric theta = 300. / t;
    const Numeric dt = theta - 1;
    const Numeric tadapt = std::pow(theta, x);

    // Temperature adaptation formulas
    constexpr auto div = [](auto& a, auto& b) { return a / b; };
    const auto linear = [dt, ta = tadapt * p](auto& Z0, auto& Z1) {
      return (Z0 + Z1 * dt) * ta;
    };
    const auto square = [dt, ta = pow2(tadapt * p)](auto& Z0, auto& Z1) {
      return (Z0 + Z1 * dt) * ta;
    };
    const auto gamma = [ta = tadapt * p](auto& Z) { return Z * ta; };
    const auto strength = [dt, tp = pow3(theta) * p](auto& a, auto& b) {
      return a * tp * std::exp(-b * dt);
    };

    // Apply transformation
    std::transform(y0.begin(), y0.end(), y1.begin(), y0.begin(), linear);
    std::transform(g0.begin(), g0.end(), g1.begin(), g0.begin(), square);
    std::transform(dv0.begin(), dv0.end(), dv1.begin(), dv0.begin(), square);
    std::transform(ga.begin(), ga.end(), ga.begin(), gamma);
    std::transform(c.begin(), c.end(), f0.begin(), c.begin(), div);
    std::transform(c.begin(), c.end(), a2.begin(), c.begin(), strength);
  }

  // Sum up positive absorption
  const Index nf = f_grid.nelem();
  for (Index iv = 0; iv < nf; iv++) {
    const Numeric f = hz2ghz(f_grid[iv]);
    if (const Numeric a = sum_lines(f, c, ga, g0, y0, dv0, f0); a > 0)
      propmat_clearsky.Kjj()[iv] += conv * oxygen_vmr * pow2(f) * a;
  }
}

void compute(PropagationMatrix& propmat_clearsky,
             ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
             const Vector& f_grid,
             const Numeric& p,
             const Numeric& t,
             const Numeric& oxygen_vmr,
             const ArrayOfRetrievalQuantity& jacs) ARTS_NOEXCEPT {
  ARTS_ASSERT(t > 0)

  if (const Index nq = dpropmat_clearsky_dx.nelem(); nq > 0) {
    // Local copies are required
    PropagationMatrix propmat(propmat_clearsky.NumberOfFrequencies(),
                              propmat_clearsky.StokesDimensions());
    PropagationMatrix dpropmat(propmat_clearsky.NumberOfFrequencies(),
                               propmat_clearsky.StokesDimensions());

    impl(propmat, f_grid, p, t, oxygen_vmr);
    propmat_clearsky.Kjj() += propmat.Kjj();

    for (Index iq = 0; iq < nq; iq++) {
      const auto& deriv = jacs[iq];

      if (not deriv.propmattype()) continue;

      if (deriv == Jacobian::Atm::Temperature) {
        ARTS_ASSERT(deriv.Target().Perturbation() not_eq 0)

        dpropmat.SetZero();

        impl(
            dpropmat, f_grid, p, t + deriv.Target().Perturbation(), oxygen_vmr);
        dpropmat.Kjj() -= propmat.Kjj();
        dpropmat.Kjj() /= deriv.Target().Perturbation();
        dpropmat_clearsky_dx[iq].Kjj() += dpropmat.Kjj();
      } else if (is_frequency_parameter(deriv)) {
        ARTS_ASSERT(deriv.Target().Perturbation() not_eq 0)

        dpropmat.SetZero();

        Vector f_pert{f_grid};
        f_pert += deriv.Target().Perturbation();

        impl(dpropmat, f_pert, p, t, oxygen_vmr);
        dpropmat.Kjj() -= propmat.Kjj();
        dpropmat.Kjj() /= deriv.Target().Perturbation();
        dpropmat_clearsky_dx[iq].Kjj() += dpropmat.Kjj();
      } else if (constexpr auto& o266 =
                     Species::select(Species::Species::Oxygen, "66");
                 deriv.Target().isSpeciesVMR() and
                 deriv.Target().QuantumIdentity().Isotopologue() == o266 and
                 oxygen_vmr > 0) {
        dpropmat = propmat;
        dpropmat /= oxygen_vmr;
        dpropmat_clearsky_dx[iq].Kjj() += dpropmat.Kjj();
      }
    }
  } else {
    impl(propmat_clearsky, f_grid, p, t, oxygen_vmr);
  }
}
}  // namespace Absorption::PredefinedModel::Makarov2020etal
