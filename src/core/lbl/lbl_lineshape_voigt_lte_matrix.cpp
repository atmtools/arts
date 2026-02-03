#include "lbl_lineshape_voigt_lte_matrix.h"

#include <arts_omp.h>
#include <jacobian.h>
#include <matpack.h>
#include <partfun.h>
#include <physics_funcs.h>

#include <Faddeeva.hh>
#include <algorithm>

#include "matpack_mdspan_common_select.h"

namespace lbl::voigt::lte::matrix {
namespace {
Complex z_(const ConstVectorView& v, const Numeric& f) {
  return {(f - v[0]) * v[1], v[2]};
}

Complex s_(const ConstVectorView& v) { return {v[3], v[4]}; }

Complex dz_(const ConstVectorView& v,
            const ConstVectorView& dv,
            const Numeric& f) {
  return {(f - v[0]) * dv[1] - dv[0] * v[1], dv[2]};
}

Complex ds_(const ConstVectorView& dv) { return {dv[3], dv[4]}; }

Complex dw_(const Complex& z, const Complex& w) {
  /*! FIXME: We should use a proper algorithm here.  This produces
   *         no errors in tests, but the actual derivative form is
   *         analytically known.  Its numerical instability, however,
   *         makes it completely useless.
   *
   *         The analytical form is:
   *
   *           dF = -2 * z * F(z) + 2 * i / sqrt(pi)
   *
   * Tests show that for y < 1e7 it works until x > 1e7, but for
   * y > 1e7, it always fails.  This is about the analytical form
   * above using the latest version of the MIT Faddeeva package.
  */
  const Complex dz{std::max(1e-4 * nonstd::abs(z.real()), 1e-4),
                   std::max(1e-4 * nonstd::abs(z.imag()), 1e-4)};
  const Complex w_2 = Faddeeva::w(z + dz);
  return (w_2 - w) / dz;
}

constexpr Complex dexpr(const Complex& w,
                        const Complex& dw,
                        const Complex& s,
                        const Complex& ds,
                        const Complex& dz) {
  return w * ds + dw * dz * s;
}

Complex dexpr_df(const Complex& dw,
                 const Complex& s,
                 const ConstVectorView& v) {
  return dw * v[1] * s;
}

constexpr Complex expr(const Complex& w, const Complex& s) { return w * s; }

Numeric scaled_gd_divf0(const Numeric& T, const Numeric& mass) {
  constexpr auto c = Constant::doppler_broadening_const_squared;
  return std::sqrt(c * T / mass);
}

void prepare_wo_jac(MatrixView& mat,
                    const AtmPoint& atm,
                    const std::span<const flat_band_data>& bands,
                    const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;
  const Numeric H = hypot(atm.mag);

#pragma omp parallel for if (arts_omp_parallel(bands.size()))
  for (const auto& flat : bands) {
    Size idx = flat.prev_size;

    const auto& ir   = flat.band_key.isot;
    const Numeric Q  = PartitionFunctions::Q(T, ir);
    const Numeric rx = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];
    const Numeric GD = scaled_gd_divf0(T, ir.mass);

    for (auto& line : flat.band) {
      const Size nz     = line.z.size(line.qn, pol);
      const Numeric s_  = line.s(T, Q);
      const Numeric G   = line.ls.G(atm);
      const Numeric Y   = line.ls.Y(atm);
      const Numeric G0  = line.ls.G0(atm);
      const Numeric D0  = line.ls.D0(atm);
      const Numeric DV  = line.ls.DV(atm);
      const Numeric lmr = 1 + G;
      const Numeric lmi = -Y;

      for (Size iz = 0; iz < nz; ++iz) {
        auto v           = mat[idx + iz];
        const Numeric Sz = line.z.Strength(line.qn, pol, iz);
        const Numeric dH = line.z.Splitting(line.qn, pol, iz);
        const Numeric s  = Sz * s_;
        v[0]             = line.f0 + D0 + DV + H * dH;
        v[1]             = 1.0 / (GD * v[0]);
        v[2]             = G0 * v[1];
        const Numeric c  = rx * s * v[1];
        v[3]             = c * lmr;
        v[4]             = c * lmi;
      }
      idx += nz;
    }
  }
}

void prepare_t_deriv(MatrixView& mat,
                     const AtmPoint& atm,
                     const std::span<const flat_band_data>& bands,
                     const Range& range,
                     const ZeemanPolarization& pol) {
  const Numeric T     = atm.temperature;
  const Numeric inv2T = 1.0 / (2.0 * T);

#pragma omp parallel for if (arts_omp_parallel(bands.size()))
  for (const auto& flat : bands) {
    Size idx = flat.prev_size;

    const auto& ir     = flat.band_key.isot;
    const Numeric Q    = PartitionFunctions::Q(T, ir);
    const Numeric dQdT = PartitionFunctions::dQdT(T, ir);
    const Numeric rx   = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];

    for (auto& line : flat.band) {
      const Size nz       = line.z.size(line.qn, pol);
      const Numeric s_    = line.s(T, Q);
      const Numeric ds_   = line.ds_dT(T, Q, dQdT);
      const Numeric G     = line.ls.G(atm);
      const Numeric Y     = line.ls.Y(atm);
      const Numeric G0    = line.ls.G0(atm);
      const Numeric dGdT  = line.ls.dG_dT(atm);
      const Numeric dYdT  = line.ls.dY_dT(atm);
      const Numeric dG0dT = line.ls.dG0_dT(atm);
      const Numeric dD0dT = line.ls.dD0_dT(atm);
      const Numeric dDVdT = line.ls.dDV_dT(atm);
      const Numeric lmr   = 1 + G;
      const Numeric lmi   = -Y;
      const Numeric dlmr  = dGdT;
      const Numeric dlmi  = -dYdT;

      for (Size iz = 0; iz < nz; ++iz) {
        auto v           = mat[idx + iz];
        auto dv          = v[range];
        const Numeric Sz = line.z.Strength(line.qn, pol, iz);
        const Numeric s  = Sz * s_;
        const Numeric ds = Sz * ds_;
        dv[0]            = dD0dT + dDVdT;
        dv[1]            = -v[1] * (inv2T + dv[0] / v[0]);
        dv[2]            = dG0dT * v[1] + G0 * dv[1];
        const Numeric c  = rx * s * v[1];
        const Numeric dc = rx * (ds * v[1] + s * dv[1]);
        dv[3]            = dc * lmr + c * dlmr;
        dv[4]            = dc * lmi + c * dlmi;
      }
      idx += nz;
    }
  }
}

void prepare_p_deriv(MatrixView& mat,
                     const AtmPoint& atm,
                     const std::span<const flat_band_data>& bands,
                     const Range& range,
                     const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;

#pragma omp parallel for if (arts_omp_parallel(bands.size()))
  for (const auto& flat : bands) {
    Size idx = flat.prev_size;

    const auto& ir   = flat.band_key.isot;
    const Numeric Q  = PartitionFunctions::Q(T, ir);
    const Numeric rx = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];

    for (auto& line : flat.band) {
      const Size nz      = line.z.size(line.qn, pol);
      const auto s_      = line.s(T, Q);
      const auto G       = line.ls.G(atm);
      const auto Y       = line.ls.Y(atm);
      const auto G0      = line.ls.G0(atm);
      const auto dGdP    = line.ls.dG_dP(atm);
      const auto dYdP    = line.ls.dY_dP(atm);
      const auto dG0dP   = line.ls.dG0_dP(atm);
      const auto dD0dP   = line.ls.dD0_dP(atm);
      const auto dDVdP   = line.ls.dDV_dP(atm);
      const Numeric lmr  = 1 + G;
      const Numeric lmi  = -Y;
      const Numeric dlmr = dGdP;
      const Numeric dlmi = -dYdP;

      for (Size iz = 0; iz < nz; ++iz) {
        auto v           = mat[idx + iz];
        auto dv          = v[range];
        const Numeric Sz = line.z.Strength(line.qn, pol, iz);
        const Numeric s  = Sz * s_;
        dv[0]            = dD0dP + dDVdP;
        dv[1]            = -v[1] * dv[0] / v[0];
        dv[2]            = dG0dP * v[1] + G0 * dv[1];
        const Numeric c  = rx * s * v[1];
        const Numeric dc = rx * s * dv[1];
        dv[3]            = dc * lmr + c * dlmr;
        dv[4]            = dc * lmi + c * dlmi;
      }
      idx += nz;
    }
  }
}

void prepare_mag_u_deriv(MatrixView& mat,
                         const AtmPoint& atm,
                         const std::span<const flat_band_data>& bands,
                         const Range& range,
                         const ZeemanPolarization& pol) {
  const Numeric dHdu = atm.mag[0] / hypot(atm.mag);

#pragma omp parallel for if (arts_omp_parallel(bands.size()))
  for (const auto& flat : bands) {
    Size idx = flat.prev_size;

    for (auto& line : flat.band) {
      const Size nz = line.z.size(line.qn, pol);
      const auto G0 = line.ls.G0(atm);

      for (Size iz = 0; iz < nz; ++iz) {
        auto v  = mat[idx + iz];
        auto dv = v[range];
        dv[0]   = line.z.Splitting(line.qn, pol, iz) * dHdu;
        dv[1]   = -v[1] * dv[0] / v[0];
        dv[2]   = G0 * dv[1];
        dv[3]   = v[3] * dv[1] / v[1];
        dv[4]   = v[4] * dv[1] / v[1];
      }
      idx += nz;
    }
  }
}

void prepare_mag_v_deriv(MatrixView& mat,
                         const AtmPoint& atm,
                         const std::span<const flat_band_data>& bands,
                         const Range& range,
                         const ZeemanPolarization& pol) {
  const Numeric dHdv = atm.mag[1] / hypot(atm.mag);

#pragma omp parallel for if (arts_omp_parallel(bands.size()))
  for (const auto& flat : bands) {
    Size idx = flat.prev_size;

    for (auto& line : flat.band) {
      const Size nz = line.z.size(line.qn, pol);
      const auto G0 = line.ls.G0(atm);

      for (Size iz = 0; iz < nz; ++iz) {
        auto v  = mat[idx + iz];
        auto dv = v[range];
        dv[0]   = line.z.Splitting(line.qn, pol, iz) * dHdv;
        dv[1]   = -v[1] * dv[0] / v[0];
        dv[2]   = G0 * dv[1];
        dv[3]   = v[3] * dv[1] / v[1];
        dv[4]   = v[4] * dv[1] / v[1];
      }
      idx += nz;
    }
  }
}

void prepare_mag_w_deriv(MatrixView& mat,
                         const AtmPoint& atm,
                         const std::span<const flat_band_data>& bands,
                         const Range& range,
                         const ZeemanPolarization& pol) {
  const Numeric dHdw = atm.mag[2] / hypot(atm.mag);

#pragma omp parallel for if (arts_omp_parallel(bands.size()))
  for (const auto& flat : bands) {
    Size idx = flat.prev_size;

    for (auto& line : flat.band) {
      const Size nz = line.z.size(line.qn, pol);
      const auto G0 = line.ls.G0(atm);

      for (Size iz = 0; iz < nz; ++iz) {
        auto v  = mat[idx + iz];
        auto dv = v[range];
        dv[0]   = line.z.Splitting(line.qn, pol, iz) * dHdw;
        dv[1]   = -v[1] * dv[0] / v[0];
        dv[2]   = G0 * dv[1];
        dv[3]   = v[3] * dv[1] / v[1];
        dv[4]   = v[4] * dv[1] / v[1];
      }
      idx += nz;
    }
  }
}

void prepare_deriv(MatrixView& mat,
                   const AtmPoint& atm,
                   const std::span<const flat_band_data>& bands,
                   const SpeciesEnum& target_spec,
                   const Range& range,
                   const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;

#pragma omp parallel for if (arts_omp_parallel(bands.size()))
  for (const auto& flat : bands) {
    Size idx              = flat.prev_size;
    const band_data& band = flat.band;

    const auto& ir = flat.band_key.isot;

    if (ir.spec != target_spec) {
      idx += band.count_zeeman_lines(pol);
      continue;
    }

    const Numeric Q   = PartitionFunctions::Q(T, ir);
    const Numeric drx = Constant::inv_sqrt_pi * atm[ir];
    const Numeric rx  = drx * atm[ir.spec];

    for (auto& line : band) {
      const Size nz         = line.z.size(line.qn, pol);
      const Numeric s_      = line.s(T, Q);
      const Numeric G       = line.ls.G(atm);
      const Numeric Y       = line.ls.Y(atm);
      const Numeric G0      = line.ls.G0(atm);
      const Numeric dGdVMR  = line.ls.dG_dVMR(atm, target_spec);
      const Numeric dYdVMR  = line.ls.dY_dVMR(atm, target_spec);
      const Numeric dG0dVMR = line.ls.dG0_dVMR(atm, target_spec);
      const Numeric dD0dVMR = line.ls.dD0_dVMR(atm, target_spec);
      const Numeric dDVdVMR = line.ls.dDV_dVMR(atm, target_spec);
      const Numeric lmr     = 1 + G;
      const Numeric lmi     = -Y;
      const Numeric dlmr    = dGdVMR;
      const Numeric dlmi    = -dYdVMR;

      for (Size iz = 0; iz < nz; ++iz) {
        auto v           = mat[idx + iz];
        auto dv          = v[range];
        const Numeric Sz = line.z.Strength(line.qn, pol, iz);
        const Numeric s  = Sz * s_;
        dv[0]            = dD0dVMR + dDVdVMR;
        dv[1]            = -v[1] * dv[0] / v[0];
        dv[2]            = dG0dVMR * v[1] + G0 * dv[1];
        const Numeric c  = rx * s * v[1];
        const Numeric dc = rx * s * dv[1] + drx * s * v[1];
        dv[3]            = dc * lmr + c * dlmr;
        dv[4]            = dc * lmi + c * dlmi;
      }
      idx += nz;
    }
  }
}

void prepare_deriv(MatrixView& mat,
                   const AtmPoint& atm,
                   const std::span<const flat_band_data>& bands,
                   const SpeciesIsotope& target_spec,
                   const Range& range,
                   const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;

#pragma omp parallel for if (arts_omp_parallel(bands.size()))
  for (const auto& flat : bands) {
    Size idx              = flat.prev_size;
    const band_data& band = flat.band;

    const auto& ir = flat.band_key.isot;

    if (ir != target_spec) {
      idx += band.count_zeeman_lines(pol);
      continue;
    }

    const Numeric Q   = PartitionFunctions::Q(T, ir);
    const Numeric drx = Constant::inv_sqrt_pi * atm[ir.spec];

    for (auto& line : band) {
      const Size nz     = line.z.size(line.qn, pol);
      const Numeric s_  = line.s(T, Q);
      const Numeric G   = line.ls.G(atm);
      const Numeric Y   = line.ls.Y(atm);
      const Numeric lmr = 1 + G;
      const Numeric lmi = -Y;

      for (Size iz = 0; iz < nz; ++iz) {
        auto v           = mat[idx + iz];
        auto dv          = v[range];
        const Numeric Sz = line.z.Strength(line.qn, pol, iz);
        const Numeric s  = Sz * s_;
        const Numeric dc = drx * s * v[1];
        dv[3]            = dc * lmr;
        dv[4]            = dc * lmi;
      }
      idx += nz;
    }
  }
}

void prepare_f0_deriv(MatrixView& mat,
                      const AtmPoint& atm,
                      const std::span<const flat_band_data>& bands,
                      const Jacobian::LineTarget& target,
                      const Range& range,
                      const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;

  for (const auto& flat : bands) {
    Size idx    = flat.prev_size;
    auto& key   = flat.band_key;
    auto& lines = flat.band.lines;

    if (target.type.band == key) {
      for (Size i = 0; i < target.type.line; ++i) {
        idx += lines[i].z.size(lines[i].qn, pol);
      }

      const auto& ir    = key.isot;
      auto& line        = lines[target.type.line];
      const Size nz     = line.z.size(line.qn, pol);
      const Numeric Q   = PartitionFunctions::Q(T, ir);
      const Numeric rx  = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];
      const Numeric s_  = line.s(T, Q);
      const Numeric ds_ = line.ds_df0(T, Q);
      const Numeric G   = line.ls.G(atm);
      const Numeric Y   = line.ls.Y(atm);
      const Numeric G0  = line.ls.G0(atm);
      const Numeric lmr = 1 + G;
      const Numeric lmi = -Y;

      for (Size iz = 0; iz < nz; ++iz) {
        auto v           = mat[idx + iz];
        auto dv          = v[range];
        const Numeric Sz = line.z.Strength(line.qn, pol, iz);
        const Numeric s  = Sz * s_;
        const Numeric ds = Sz * ds_;
        dv[0]            = 1;
        dv[1]            = -v[1] * dv[0] / v[0];
        dv[2]            = G0 * dv[1];
        const Numeric dc = rx * (ds * v[1] + s * dv[1]);
        dv[3]            = dc * lmr;
        dv[4]            = dc * lmi;
      }

      return;
    }
  }
}

void prepare_e0_deriv(MatrixView& mat,
                      const AtmPoint& atm,
                      const std::span<const flat_band_data>& bands,
                      const Jacobian::LineTarget& target,
                      const Range& range,
                      const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;

  for (const auto& flat : bands) {
    Size idx    = flat.prev_size;
    auto& key   = flat.band_key;
    auto& lines = flat.band.lines;

    if (target.type.band == key) {
      for (Size i = 0; i < target.type.line; ++i) {
        idx += lines[i].z.size(lines[i].qn, pol);
      }

      const auto& ir    = key.isot;
      const Numeric rx  = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];
      const Numeric Q   = PartitionFunctions::Q(T, ir);
      auto& line        = lines[target.type.line];
      const Size nz     = line.z.size(line.qn, pol);
      const Numeric ds_ = line.ds_de0(T, Q);
      const Numeric G   = line.ls.G(atm);
      const Numeric Y   = line.ls.Y(atm);
      const Numeric lmr = 1 + G;
      const Numeric lmi = -Y;

      for (Size iz = 0; iz < nz; ++iz) {
        auto v           = mat[idx + iz];
        auto dv          = v[range];
        const Numeric Sz = line.z.Strength(line.qn, pol, iz);
        const Numeric ds = Sz * ds_;
        const Numeric dc = rx * ds * v[1];
        dv[3]            = dc * lmr;
        dv[4]            = dc * lmi;
      }

      return;
    }
  }
}

void prepare_a_deriv(MatrixView& mat,
                     const AtmPoint& atm,
                     const std::span<const flat_band_data>& bands,
                     const Jacobian::LineTarget& target,
                     const Range& range,
                     const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;

  for (const auto& flat : bands) {
    Size idx    = flat.prev_size;
    auto& key   = flat.band_key;
    auto& lines = flat.band.lines;

    if (target.type.band == key) {
      for (Size i = 0; i < target.type.line; ++i) {
        idx += lines[i].z.size(lines[i].qn, pol);
      }

      const auto& ir    = key.isot;
      const Numeric Q   = PartitionFunctions::Q(T, ir);
      const Numeric rx  = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];
      auto& line        = lines[target.type.line];
      const Size nz     = line.z.size(line.qn, pol);
      const Numeric ds_ = line.ds_da(T, Q);
      const Numeric G   = line.ls.G(atm);
      const Numeric Y   = line.ls.Y(atm);
      const Numeric lmr = 1 + G;
      const Numeric lmi = -Y;

      for (Size iz = 0; iz < nz; ++iz) {
        auto v           = mat[idx + iz];
        auto dv          = v[range];
        const Numeric Sz = line.z.Strength(line.qn, pol, iz);
        const Numeric ds = Sz * ds_;
        const Numeric dc = rx * ds * v[1];
        dv[3]            = dc * lmr;
        dv[4]            = dc * lmi;
      }

      return;
    }
  }
}

void prepare_g0_deriv(MatrixView& mat,
                      const AtmPoint& atm,
                      const std::span<const flat_band_data>& bands,
                      const Jacobian::LineTarget& target,
                      const Range& range,
                      const ZeemanPolarization& pol) {
  for (const auto& flat : bands) {
    Size idx    = flat.prev_size;
    auto& key   = flat.band_key;
    auto& lines = flat.band.lines;

    if (target.type.band == key) {
      for (Size i = 0; i < target.type.line; ++i) {
        idx += lines[i].z.size(lines[i].qn, pol);
      }

      auto& line    = lines[target.type.line];
      const Size nz = line.z.size(line.qn, pol);
      const Numeric dG0 =
          line.ls.dG0_dX(atm, target.type.band.isot.spec, target.type.ls_coeff);

      for (Size iz = 0; iz < nz; ++iz) {
        auto v  = mat[idx + iz];
        auto dv = v[range];
        dv[2]   = dG0 * v[1];
      }

      return;
    }
  }
}

void prepare_d0_deriv(MatrixView& mat,
                      const AtmPoint& atm,
                      const std::span<const flat_band_data>& bands,
                      const Jacobian::LineTarget& target,
                      const Range& range,
                      const ZeemanPolarization& pol) {
  for (const auto& flat : bands) {
    Size idx    = flat.prev_size;
    auto& key   = flat.band_key;
    auto& lines = flat.band.lines;

    if (target.type.band == key) {
      for (Size i = 0; i < target.type.line; ++i) {
        idx += lines[i].z.size(lines[i].qn, pol);
      }

      auto& line       = lines[target.type.line];
      const Size nz    = line.z.size(line.qn, pol);
      const Numeric G0 = line.ls.G0(atm);
      const Numeric dD0 =
          line.ls.dD0_dX(atm, target.type.band.isot.spec, target.type.ls_coeff);

      for (Size iz = 0; iz < nz; ++iz) {
        auto v  = mat[idx + iz];
        auto dv = v[range];
        dv[0]   = dD0;
        dv[1]   = -v[1] * dv[0] / v[0];
        dv[2]   = G0 * dv[1];
      }

      return;
    }
  }
}

void prepare_dv_deriv(MatrixView& mat,
                      const AtmPoint& atm,
                      const std::span<const flat_band_data>& bands,
                      const Jacobian::LineTarget& target,
                      const Range& range,
                      const ZeemanPolarization& pol) {
  for (const auto& flat : bands) {
    Size idx    = flat.prev_size;
    auto& key   = flat.band_key;
    auto& lines = flat.band.lines;

    if (target.type.band == key) {
      for (Size i = 0; i < target.type.line; ++i) {
        idx += lines[i].z.size(lines[i].qn, pol);
      }

      auto& line       = lines[target.type.line];
      const Size nz    = line.z.size(line.qn, pol);
      const Numeric G0 = line.ls.G0(atm);
      const Numeric dDV =
          line.ls.dDV_dX(atm, target.type.band.isot.spec, target.type.ls_coeff);

      for (Size iz = 0; iz < nz; ++iz) {
        auto v  = mat[idx + iz];
        auto dv = v[range];
        dv[0]   = dDV;
        dv[1]   = -v[1] * dv[0] / v[0];
        dv[2]   = G0 * dv[1];
      }

      return;
    }
  }
}

void prepare_y_deriv(MatrixView& mat,
                     const AtmPoint& atm,
                     const std::span<const flat_band_data>& bands,
                     const Jacobian::LineTarget& target,
                     const Range& range,
                     const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;

  for (const auto& flat : bands) {
    Size idx    = flat.prev_size;
    auto& key   = flat.band_key;
    auto& lines = flat.band.lines;

    if (target.type.band == key) {
      for (Size i = 0; i < target.type.line; ++i) {
        idx += lines[i].z.size(lines[i].qn, pol);
      }

      const auto& ir   = key.isot;
      const Numeric Q  = PartitionFunctions::Q(T, ir);
      const Numeric rx = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];
      auto& line       = lines[target.type.line];
      const Size nz    = line.z.size(line.qn, pol);
      const Numeric s_ = line.s(T, Q);

      for (Size iz = 0; iz < nz; ++iz) {
        auto dv          = mat[idx + iz, range];
        const Numeric Sz = line.z.Strength(line.qn, pol, iz);
        const Numeric s  = Sz * s_;
        dv[4]            = -rx * s;
      }

      return;
    }
  }
}

void prepare_g_deriv(MatrixView& mat,
                     const AtmPoint& atm,
                     const std::span<const flat_band_data>& bands,
                     const Jacobian::LineTarget& target,
                     const Range& range,
                     const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;

  for (const auto& flat : bands) {
    Size idx        = flat.prev_size;
    const auto& key = flat.band_key;
    auto& lines     = flat.band.lines;

    if (target.type.band == key) {
      for (Size i = 0; i < target.type.line; ++i) {
        idx += lines[i].z.size(lines[i].qn, pol);
      }

      const auto& ir   = key.isot;
      const Numeric Q  = PartitionFunctions::Q(T, ir);
      const Numeric rx = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];
      auto& line       = lines[target.type.line];
      const Size nz    = line.z.size(line.qn, pol);
      const Numeric s_ = line.s(T, Q);

      for (Size iz = 0; iz < nz; ++iz) {
        auto dv          = mat[idx + iz, range];
        const Numeric Sz = line.z.Strength(line.qn, pol, iz);
        const Numeric s  = Sz * s_;
        dv[3]            = rx * s;
      }

      return;
    }
  }
}

void prepare_deriv(MatrixView& mat,
                   const AtmPoint& atm,
                   const std::span<const flat_band_data>& bands,
                   const AtmKey& key,
                   const Range& range,
                   const ZeemanPolarization& pol) {
  switch (key) {
    using enum AtmKey;
    case t:      prepare_t_deriv(mat, atm, bands, range, pol); break;
    case p:      prepare_p_deriv(mat, atm, bands, range, pol); break;
    case mag_u:  prepare_mag_u_deriv(mat, atm, bands, range, pol); break;
    case mag_v:  prepare_mag_v_deriv(mat, atm, bands, range, pol); break;
    case mag_w:  prepare_mag_w_deriv(mat, atm, bands, range, pol); break;
    case wind_u: [[fallthrough]];
    case wind_v: [[fallthrough]];
    case wind_w: break;
  }
}

void prepare_deriv(MatrixView&,
                   const AtmPoint&,
                   const std::span<const flat_band_data>&,
                   const ScatteringSpeciesProperty&,
                   const Range&,
                   const ZeemanPolarization&) {}

void prepare_deriv(MatrixView&,
                   const AtmPoint&,
                   const std::span<const flat_band_data>&,
                   const QuantumLevelIdentifier&,
                   const Range&,
                   const ZeemanPolarization&) {}

void prepare_deriv(MatrixView& mat,
                   const AtmPoint& atm,
                   const std::span<const flat_band_data>& bands,
                   const Jacobian::LineTarget& key,
                   const Range& range,
                   const ZeemanPolarization& pol) {
  switch (key.type.var) {
    using enum LineByLineVariable;
    case f0:     prepare_f0_deriv(mat, atm, bands, key, range, pol); break;
    case e0:     prepare_e0_deriv(mat, atm, bands, key, range, pol); break;
    case a:      prepare_a_deriv(mat, atm, bands, key, range, pol); break;
    case unused: break;
  }

  switch (key.type.ls_var) {
    using enum LineShapeModelVariable;
    case G0:     prepare_g0_deriv(mat, atm, bands, key, range, pol); break;
    case D0:     prepare_d0_deriv(mat, atm, bands, key, range, pol); break;
    case Y:      prepare_y_deriv(mat, atm, bands, key, range, pol); break;
    case G:      prepare_g_deriv(mat, atm, bands, key, range, pol); break;
    case DV:     prepare_dv_deriv(mat, atm, bands, key, range, pol); break;
    case G2:     [[fallthrough]];
    case D2:     [[fallthrough]];
    case FVC:    [[fallthrough]];
    case ETA:    return;
    case unused: break;
  }
}

void prepare_with_jac(MatrixView mat,
                      const AtmPoint& atm,
                      const std::span<const flat_band_data>& bands,
                      const JacobianTargets& jac_targets,
                      const ZeemanPolarization& pol) {
  assert(static_cast<Size>(static_cast<Size>(mat.ncols())) ==
         jac_targets.target_count() * 5 + 5);

  prepare_wo_jac(mat, atm, bands, pol);

  for (auto& target : jac_targets.atm) {
    std::visit(
        [&](const auto& x) {
          prepare_deriv(
              mat, atm, bands, x, {target.target_pos * 5 + 5, 5}, pol);
        },
        target.type);
  }

  for (auto& target : jac_targets.line) {
    prepare_deriv(mat, atm, bands, target, {target.target_pos * 5 + 5, 5}, pol);
  }
}
}  // namespace

void str_scale(ComplexVectorView a,
               const AtmPoint& atm,
               const ConstVectorView& fs) {
  using Constant::pi, Constant::c, Constant::h, Constant::k;
  constexpr Numeric sc = c * c / (8 * pi);

  const Size n = static_cast<Size>(a.ncols());

  assert(static_cast<Size>(a.nrows()) == n);

  const Numeric T = atm.temperature;
  const Numeric P = atm.pressure;
  const Numeric N = number_density(P, T);

#pragma omp parallel for if (arts_omp_parallel(n))
  for (Size i = 0; i < n; i++) {
    const auto f       = fs[i];
    const Numeric r    = (h * f) / (k * T);
    const Numeric e    = std::expm1(-r);
    const Numeric scl  = -N * f * e * sc;
    a[i]              *= scl;
  }
}

void str_scale(ComplexMatrixView a,
               const AtmPoint& atm,
               const ConstVectorView& fs,
               const std::vector<bool>& df,
               const Size it) {
  using Constant::pi, Constant::c, Constant::h, Constant::k;
  constexpr Numeric sc = c * c / (8 * pi);

  const Size nf = fs.ncols();
  const Size nq = df.size();

  if (nq == 0) return str_scale(a.view_as(nf), atm, fs);

  assert(static_cast<Size>(a.nrows()) == nf);
  assert(static_cast<Size>(a.ncols()) == (1 + nq));

  const Numeric T = atm.temperature;
  const Numeric P = atm.pressure;
  const Numeric N = number_density(P, T);

#pragma omp parallel for if (arts_omp_parallel(nf))
  for (Size iv = 0; iv < nf; iv++) {
    const auto f    = fs[iv];
    const Numeric r = (h * f) / (k * T);

    const Numeric e   = std::expm1(-r);
    const Numeric scl = -N * f * e * sc;

    a[iv, Range{1, nq}] *= scl;

    if (it > 0) {
      a[iv, it] -= a[iv, 0] * f * N * (r * (e + 1) - e) * sc / T;
    }

    for (Size iq = 0; iq < nq; iq++) {
      if (df[iq]) {
        a[iv, iq + 1] += a[iv, 0] * N * (r * (e + 1) - e) * sc;
      }
    }

    a[iv, 0] *= scl;
  }
}

void sumup(ComplexVectorView a,
           const ConstMatrixView& mat,
           const ConstVectorView& f) {
  const Size n = mat.nrows();
  const Size m = f.ncols();
  assert(mat.ncols() >= 5);
  assert(static_cast<Size>(a.ncols()) == m);

#pragma omp parallel for if (arts_omp_parallel(m))
  for (Size j = 0; j < m; ++j) {
    const auto& f_ = f[j];
    auto& aj       = a[j];
    for (Size i = 0; i < n; ++i) {
      const auto&& v  = mat[i];
      aj             += expr(Faddeeva::w(z_(v, f_)), s_(v));
    }
  }
}

void sumup(ComplexVectorView a,
           const ConstMatrixView& mat,
           const ConstVectorView& f,
           const Numeric cutoff) {
  const Size m = f.ncols();
  assert(mat.ncols() >= 5);
  assert(static_cast<Size>(a.ncols()) == m);

  if (m == 0) return;

  constexpr auto fproj = [](const auto& v) { return v[0]; };

  const auto lower = stdr::lower_bound(mat, f.front() - cutoff, {}, fproj);
  const auto upper = stdr::lower_bound(mat, f.back() + cutoff, {}, fproj);

  Size i0       = lower.pos;
  const Size i1 = upper.pos;

  const ComplexVector cutoffs = [i0, i1, cutoff, &mat]() -> ComplexVector {
    ComplexVector res(mat.nrows());

#pragma omp parallel for if (arts_omp_parallel(i1-i0))
    for (Size i = i0; i < i1; i++) {
      const auto&& v = mat[i];
      res[i]         = expr(Faddeeva::w(z_(v, v[0] + cutoff)), s_(v));
    }

    return res;
  }();

#pragma omp parallel for if (arts_omp_parallel(m)) firstprivate(i0)
  for (Size j = 0; j < m; ++j) {
    const auto& f_ = f[j];

    while (i0 < i1 and (mat[i0, 0] < (f_ - cutoff))) ++i0;

    auto& aj = a[j];

    for (Size i = i0; i < i1; ++i) {
      const auto&& v = mat[i];

      if ((v[0] > (f_ + cutoff))) break;

      aj += expr(Faddeeva::w(z_(v, f_)), s_(v)) - cutoffs[i];
    }
  }
}

void sumup(ComplexMatrixView res,
           const ConstMatrixView& mat,
           const ConstVectorView& fs,
           const std::vector<bool>& df) {
  const Size nf = fs.size();
  const Size nl = mat.nrows();
  const Size nq = df.size();

  assert(nf == static_cast<Size>(res.nrows()));
  assert(nq == static_cast<Size>(res.ncols() - 1));
  assert(5 + 5 * nq == static_cast<Size>(mat.ncols()));

#pragma omp parallel for if (arts_omp_parallel(nf))
  for (Size iv = 0; iv < nf; ++iv) {
    const auto& f = fs[iv];
    auto aj       = res[iv];
    for (Size il = 0; il < nl; ++il) {
      const auto v = mat[il];

      const Complex z  = z_(v, f);
      const Complex w  = Faddeeva::w(z);
      const Complex dw = dw_(z, w);
      const Complex s  = s_(v);

      aj[0] += w * s;

      for (Size iq = 0; iq < nq; ++iq) {
        const auto dv  = v[Range(5 + iq * 5, 5)];
        aj[iq + 1]    += df[iq] ? dexpr_df(dw, s, v)
                                : dexpr(w, dw, s, ds_(dv), dz_(v, dv, f));
      }
    }
  }
}

void sumup(ComplexMatrixView res,
           const ConstMatrixView& mat,
           const ConstVectorView& fs,
           const Numeric cutoff,
           const std::vector<bool>& df) {
  const Size nf = fs.size();
  const Size nq = df.size();

  assert(nf == static_cast<Size>(res.nrows()));
  assert(nq == static_cast<Size>(res.ncols() - 1));
  assert((5 + 5 * nq) == static_cast<Size>(mat.ncols()));

  if (nf == 0) return;

  constexpr auto fproj = [](const auto& v) { return v[0]; };

  const auto lower = stdr::lower_bound(mat, fs.front() - cutoff, {}, fproj);
  const auto upper = stdr::lower_bound(mat, fs.back() + cutoff, {}, fproj);

  Size i0       = lower.pos;
  const Size i1 = upper.pos;

  const ComplexMatrix cutoffs =
      [nq, i0, i1, cutoff, &mat, &df]() -> ComplexMatrix {
    ComplexMatrix res(mat.nrows(), 1 + nq);

#pragma omp parallel for if (arts_omp_parallel(i1-i0))
    for (Size i = i0; i < i1; i++) {
      const auto&& v = mat[i];

      const Numeric f  = v[0] + cutoff;
      const Complex z  = z_(v, f);
      const Complex w  = Faddeeva::w(z);
      const Complex dw = dw_(z, w);
      const Complex s  = s_(v);

      res[i, 0] += w * s;

      for (Size iq = 0; iq < nq; ++iq) {
        const auto dv   = v[Range(5 + iq * 5, 5)];
        res[i, iq + 1] += df[iq] ? dexpr_df(dw, s, v)
                                 : dexpr(w, dw, s, ds_(dv), dz_(v, dv, f));
      }
    }

    return res;
  }();

#pragma omp parallel for if (arts_omp_parallel(nf)) firstprivate(i0)
  for (Size iv = 0; iv < nf; ++iv) {
    const auto& f = fs[iv];

    while (i0 < i1 and (mat[i0, 0] < (f - cutoff))) ++i0;

    auto aj = res[iv];
    for (Size i = i0; i < i1; ++i) {
      const auto v = mat[i];

      if ((v[0] > (f + cutoff))) break;

      const Complex z  = z_(v, f);
      const Complex w  = Faddeeva::w(z);
      const Complex dw = dw_(z, w);
      const Complex s  = s_(v);

      aj[0] += w * s;

      for (Size iq = 0; iq < nq; ++iq) {
        const auto dv  = v[Range(5 + iq * 5, 5)];
        aj[iq + 1]    += df[iq] ? dexpr_df(dw, s, v)
                                : dexpr(w, dw, s, ds_(dv), dz_(v, dv, f));
      }

      aj -= cutoffs[i];
    }
  }
}

void prepare(Matrix& mat,
             const AtmPoint& atm,
             const std::span<const flat_band_data>& bands,
             const ZeemanPolarization& pol) {
  if (bands.empty()) {
    mat.resize(0, 5);
    return;
  }

  auto& last = bands.back();
  auto n     = last.prev_size + last.band.count_zeeman_lines(pol);
  mat.resize(n, 5);
  prepare_wo_jac(mat, atm, bands, pol);
}

void prepare(Matrix& mat,
             const AtmPoint& atm,
             const AbsorptionBands& bands,
             const ZeemanPolarization& pol) {
  auto flat = flatter_view(bands, pol, [](auto&...) { return true; });
  prepare(mat, atm, flat, pol);
}

void prepare(Matrix& mat,
             const AtmPoint& atm,
             const std::span<const flat_band_data>& bands,
             const JacobianTargets& jac_targets,
             const ZeemanPolarization& pol) {
  const Size m = jac_targets.target_count();

  if (bands.empty()) {
    mat.resize(0, 5 * (1 + m));
    return;
  }

  auto& last = bands.back();
  auto n     = last.prev_size + last.band.count_zeeman_lines(pol);
  mat.resize(n, 5 * (1 + m));
  mat = 0;  // reset for derivatives
  prepare_with_jac(mat, atm, bands, jac_targets, pol);
}

void prepare(Matrix& mat,
             const AtmPoint& atm,
             const AbsorptionBands& bands,
             const JacobianTargets& jac_targets,
             const ZeemanPolarization& pol) {
  auto flat = flatter_view(bands, pol, [](auto&...) { return true; });
  prepare(mat, atm, flat, jac_targets, pol);
}

void sort(MatrixView mat) {
  matpack::sort(mat, {}, [](auto&& a) { return a[0]; });
}
}  // namespace lbl::voigt::lte::matrix
