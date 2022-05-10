#include <interpolation_lagrange.h>
#include <propagationmatrix.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>

#include "constants.h"
#include "debug.h"

#include "predef_data.h"

/*!  This code is ported from Fortran90 code that contains the following statement:

!  --------------------------------------------------------------------------
! |  Copyright �, Atmospheric and Environmental Research, Inc., 2022         |
! |                                                                          |
! |  All rights reserved. This source code was developed as part of the      |
! |  LBLRTM software and is designed for scientific and research purposes.   |
! |  Atmospheric and Environmental Research Inc. (AER) grants USER the right |
! |  to download, install, use and copy this software for scientific and     |
! |  research purposes only. This software may be redistributed as long as   |
! |  this copyright notice is reproduced on any copy made and appropriate    |
! |  acknowledgment is given to AER. This software or any modified version   |
! |  of this software may not be incorporated into proprietary software or   |
! |  commercial software offered for sale without the express written        |
! |  consent of AER.                                                         |
! |                                                                          |
! |  This software is provided as is without any express or implied          |
! |  warranties.                                                             |
!  --------------------------------------------------------------------------
*/

namespace Absorption::PredefinedModel::Hitran::MTCKD {
Numeric RADFN_FUN(const Numeric XVI, const Numeric XKT) noexcept {
  // ---------------------------------------------------------------------- B18060
  //              LAST MODIFICATION:    12 AUGUST 1991                      B17940
  //                                                                        B17950
  //                 IMPLEMENTATION:    R.D. WORSHAM                        B17960
  //                                                                        B17970
  //            ALGORITHM REVISIONS:    S.A. CLOUGH                         B17980
  //                                    R.D. WORSHAM                        B17990
  //                                    J.L. MONCET                         B18000
  //                                                                        B18010
  //                                                                        B18020
  //                    ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         B18030
  //                    840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          B18040
  //                                                                        B18050
  //                                                                        B18070
  //              WORK SUPPORTED BY:    THE ARM PROGRAM                     B18080
  //                                    OFFICE OF ENERGY RESEARCH           B18090
  //                                    DEPARTMENT OF ENERGY                B18100
  //                                                                        B18110
  //                                                                        B18120
  //     SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             B18130
  //                                                                        B18140
  //                                            FASCOD3                     B18150
  //                                                                        B18160
  // ---------------------------------------------------------------------- B18060
  //                                                                        B18170
  //  IN THE SMALL XVIOKT REGION 0.5 IS REQUIRED

  if (XKT > 0.0) {
    const Numeric XVIOKT = XVI / XKT;

    if (XVIOKT <= 0.01) return 0.5 * XVIOKT * XVI;

    if (XVIOKT <= 10) {
      Numeric EXPVKT = std::expm1(-XVIOKT);
      return -XVI * EXPVKT / (2 + EXPVKT);
    }

    return XVI;
  }
  return XVI;
}

/*!  Single XINT function to compute the interpolation of some data to a point X

  @param[in] P (X1 - X) / (X1 - X0) if A is on X0, X1, X2, X3 and X is the target of the interpolation
  @param[in] A The value of the original function at four positions (at X0, X1, X2, X3)
*/
constexpr Numeric XINT_FUN(const Numeric P,
                           const std::array<Numeric, 4>& A) noexcept {
  const Numeric C = (3 - 2 * P) * P * P;
  const Numeric B = 0.5 * P * (1 - P);
  const Numeric B1 = B * (1 - P);
  const Numeric B2 = B * P;

  return -A[0] * B1 + A[1] * (1 - C + B2) + A[2] * (C + B1) - A[3] * B2;
}

void check(const WaterData& data) {
  bool c = data.self_absco_ref.size() == 0 or data.for_absco_ref.size() == 0 or
           data.wavenumbers.size() == 0 or data.self_texp.size() == 0;
  ARTS_USER_ERROR_IF(c, "No data")
}

void compute_foreign_h2o(PropagationMatrix& propmat_clearsky,
                         const Vector& f_grid,
                         const Numeric& P,
                         const Numeric& T,
                         const Numeric& vmrh2o,
                         const WaterData& data) {
  using Conversion::freq2kaycm;

  // Perform checks to ensure the calculation data is good
  check(data);

  const Index n{f_grid.nelem()};
  if (n == 0) return;
  if (freq2kaycm(f_grid[0]) > data.wavenumbers.back()) return;

  // Constants
  constexpr Numeric RADCN2 = 1.4387752;
  constexpr Numeric P0 = 101300;
  constexpr Numeric T0 = 296;

  // Data constants
  const Numeric last_wavenumber = data.wavenumbers.back();
  const Numeric dvc = data.wavenumbers[1] - data.wavenumbers[0];
  const Numeric recdvc = 1 / dvc;
  const auto data_size = Index(data.wavenumbers.size());

  // Level constants
  const Numeric xkt = T / RADCN2;
  const Numeric rho_rat = (P / P0) * (T0 / T);
  const Numeric num_den_cm2 = 1e-6 * vmrh2o * P / (Constant::k * T);

  // Scaling logic
  auto scl = [vmrh2o, rho_rat, xkt](auto& input_value, auto& input_wavenumber) {
    return input_value * (1.0 - vmrh2o) * rho_rat *
           RADFN_FUN(input_wavenumber, xkt);
  };

  // Compute data
  Index cur = std::distance(data.wavenumbers.begin(),
                            std::lower_bound(data.wavenumbers.begin(),
                                             data.wavenumbers.end(),
                                             freq2kaycm(f_grid[0]) - 2 * dvc));
  const Numeric* v = data.wavenumbers.data() + cur;
  const Numeric* y = data.for_absco_ref.data() + cur;
  std::array<Numeric, 4> k{0, 0, 0, 0};
  for (auto i = -1; i < 3 and cur + i < data_size; i++) {
    if (i < 0 and cur == 0)
      k[i + 1] = scl(*(y + i + 2), *(v + i + 2));
    else
      k[i + 1] = scl(*(y + i), *(v + i));
  }

  // Compute loop
  for (Index s = 0; s < n; ++s) {
    if (f_grid[s] < 0) continue;
    auto x = freq2kaycm(f_grid[s]);
    if (x > last_wavenumber) return;

    while (x > *(v + 1)) {
      std::rotate(k.begin(), k.begin() + 1, k.end());

      k.back() = data_size > cur + 3 ? scl(*(y + 3), *(v + 3)) : 0;

      cur++;
      v++;
      y++;
    }

    auto out = 1e2 * num_den_cm2 * XINT_FUN(recdvc * (x - *v), k);
    propmat_clearsky.Kjj()[s] += out >= 0 ? out : 0;
  }
}

void compute_self_h2o(PropagationMatrix& propmat_clearsky,
                      const Vector& f_grid,
                      const Numeric& P,
                      const Numeric& T,
                      const Numeric& vmrh2o,
                      const WaterData& data) {
  using Conversion::freq2kaycm;

  // Perform checks to ensure the calculation data is good
  check(data);

  const Index n{f_grid.nelem()};
  if (n == 0) return;
  if (freq2kaycm(f_grid[0]) > data.wavenumbers.back()) return;

  // Constants
  constexpr Numeric RADCN2 = 1.4387752;
  constexpr Numeric P0 = 101300;
  constexpr Numeric T0 = 296;

  // Data constants
  const Numeric last_wavenumber = data.wavenumbers.back();
  const Numeric dvc = data.wavenumbers[1] - data.wavenumbers[0];
  const Numeric recdvc = 1 / dvc;
  const auto data_size = Index(data.wavenumbers.size());

  // Level constants
  const Numeric xkt = T / RADCN2;
  const Numeric rho_rat = (P / P0) * (T0 / T);
  const Numeric num_den_cm2 = 1e-6 * vmrh2o * P / (Constant::k * T);

  // Scaling logic
  auto scl = [vmrh2o, rho_rat, xkt, r = T0 / T](auto& input_value,
                                                auto& input_exponent,
                                                auto& input_wavenumber) {
    return input_value * vmrh2o * rho_rat * std::pow(r, input_exponent) *
           RADFN_FUN(input_wavenumber, xkt);
  };

  // Compute data
  Index cur = std::distance(data.wavenumbers.begin(),
                            std::lower_bound(data.wavenumbers.begin(),
                                             data.wavenumbers.end(),
                                             freq2kaycm(f_grid[0]) - 2 * dvc));
  const Numeric* v = data.wavenumbers.data() + cur;
  const Numeric* y = data.self_absco_ref.data() + cur;
  const Numeric* e = data.self_texp.data() + cur;
  std::array<Numeric, 4> k{0, 0, 0, 0};
  for (auto i = -1; i < 3 and cur + i < data_size; i++) {
    if (i < 0 and cur == 0)
      k[i + 1] = scl(*(y + i + 2), *(e + i + 2), *(v + i + 2));
    else
      k[i + 1] = scl(*(y + i), *(e + i), *(v + i));
  }

  // Compute loop
  for (Index s = 0; s < n; ++s) {
    if (f_grid[s] < 0) continue;
    auto x = freq2kaycm(f_grid[s]);
    if (x > last_wavenumber) return;

    while (x > *(v + 1)) {
      std::rotate(k.begin(), k.begin() + 1, k.end());

      k.back() = data_size > cur + 3 ? scl(*(y + 3), *(e + 3), *(v + 3)) : 0;

      cur++;
      v++;
      y++;
      e++;
    }

    auto out = 1e2 * num_den_cm2 * XINT_FUN(recdvc * (x - *v), k);
    propmat_clearsky.Kjj()[s] += out >= 0 ? out : 0;
  }
}
}  // namespace Absorption::PredefinedModel::Hitran::MTCKD