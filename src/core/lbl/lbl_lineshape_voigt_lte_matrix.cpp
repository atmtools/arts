#include "lbl_lineshape_voigt_lte_matrix.h"

#include <arts_omp.h>
#include <jacobian.h>
#include <matpack.h>
#include <partfun.h>
#include <physics_funcs.h>

#include <Faddeeva.hh>

#include "matpack_mdspan_common_select.h"

namespace lbl::voigt::lte::matrix {
namespace {
Complex z_(const ConstVectorView& v, const Numeric& f) {
  return {(v[0] - f) * v[1], v[2]};
}

Complex s_(const ConstVectorView& v) { return {v[3], v[4]}; }

Complex dz_(const ConstVectorView& v,
            const ConstVectorView& dv,
            const Numeric& f) {
  return {(v[0] - f) * dv[1] + dv[0] * v[1], dv[2]};
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
  return -dw * v[1] * s;
}

constexpr Complex expr(const Complex& w, const Complex& s) { return w * s; }

Numeric scaled_gd_divf0(const Numeric& T, const Numeric& mass) {
  constexpr auto c = Constant::doppler_broadening_const_squared;
  return std::sqrt(c * T / mass);
}

void preparex(MatrixView& mat,
              const AtmPoint& atm,
              const AbsorptionBands& bands,
              const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;
  const Numeric H = hypot(atm.mag);
  Size idx        = 0;

  for (const auto& [key, band] : bands) {
    const auto& ir   = key.isot;
    const Numeric Q  = PartitionFunctions::Q(T, ir);
    const Numeric rx = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];
    const Numeric GD = scaled_gd_divf0(T, ir.mass);

    for (auto& line : band) {
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
                     const AbsorptionBands& bands,
                     const Range& range,
                     const ZeemanPolarization& pol) {
  const Numeric T     = atm.temperature;
  const Numeric inv2T = 1.0 / (2.0 * T);
  Size idx            = 0;

  for (const auto& [key, band] : bands) {
    const auto& ir     = key.isot;
    const Numeric Q    = PartitionFunctions::Q(T, ir);
    const Numeric dQdT = PartitionFunctions::dQdT(T, ir);
    const Numeric rx   = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];

    for (auto& line : band.lines) {
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
                     const AbsorptionBands& bands,
                     const Range& range,
                     const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;
  Size idx        = 0;

  for (const auto& [key, band] : bands) {
    const auto& ir   = key.isot;
    const Numeric Q  = PartitionFunctions::Q(T, ir);
    const Numeric rx = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];

    for (auto& line : band.lines) {
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

void prepare_mu_deriv(MatrixView& mat,
                      const AtmPoint& atm,
                      const AbsorptionBands& bands,
                      const Range& range,
                      const ZeemanPolarization& pol) {
  Size idx           = 0;
  const Numeric dHdu = atm.mag[0] / hypot(atm.mag);

  for (const auto& [key, band] : bands) {
    for (auto& line : band.lines) {
      const Size nz = line.z.size(line.qn, pol);
      const auto G0 = line.ls.G0(atm);

      for (Size iz = 0; iz < nz; ++iz) {
        auto v  = mat[idx + iz];
        auto dv = v[range];
        dv[0]   = line.z.Splitting(line.qn, pol, iz) * dHdu;
        dv[1]   = -v[1] * dv[0] / v[0];
        dv[2]   = G0 * dv[1];
      }
      idx += nz;
    }
  }
}

void prepare_mv_deriv(MatrixView& mat,
                      const AtmPoint& atm,
                      const AbsorptionBands& bands,
                      const Range& range,
                      const ZeemanPolarization& pol) {
  Size idx           = 0;
  const Numeric dHdu = atm.mag[1] / hypot(atm.mag);

  for (const auto& [key, band] : bands) {
    for (auto& line : band.lines) {
      const Size nz = line.z.size(line.qn, pol);
      const auto G0 = line.ls.G0(atm);

      for (Size iz = 0; iz < nz; ++iz) {
        auto v  = mat[idx + iz];
        auto dv = v[range];
        dv[0]   = line.z.Splitting(line.qn, pol, iz) * dHdu;
        dv[1]   = -v[1] * dv[0] / v[0];
        dv[2]   = G0 * dv[1];
      }
      idx += nz;
    }
  }
}

void prepare_mw_deriv(MatrixView& mat,
                      const AtmPoint& atm,
                      const AbsorptionBands& bands,
                      const Range& range,
                      const ZeemanPolarization& pol) {
  Size idx           = 0;
  const Numeric dHdu = atm.mag[1] / hypot(atm.mag);

  for (const auto& [key, band] : bands) {
    for (auto& line : band.lines) {
      const Size nz = line.z.size(line.qn, pol);
      const auto G0 = line.ls.G0(atm);

      for (Size iz = 0; iz < nz; ++iz) {
        auto v  = mat[idx + iz];
        auto dv = v[range];
        dv[0]   = line.z.Splitting(line.qn, pol, iz) * dHdu;
        dv[1]   = -v[1] * dv[0] / v[0];
        dv[2]   = G0 * dv[1];
      }
      idx += nz;
    }
  }
}

void prepare_deriv(MatrixView& mat,
                   const AtmPoint& atm,
                   const AbsorptionBands& bands,
                   const SpeciesEnum& target_spec,
                   const Range& range,
                   const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;
  Size idx        = 0;

  for (const auto& [key, band] : bands) {
    const auto& ir = key.isot;

    if (ir.spec != target_spec) {
      idx += band.lines.size();
      continue;
    }

    const Numeric Q   = PartitionFunctions::Q(T, ir);
    const Numeric drx = Constant::inv_sqrt_pi * atm[ir];
    const Numeric rx  = drx * atm[ir.spec];

    for (auto& line : band.lines) {
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
                   const AbsorptionBands& bands,
                   const SpeciesIsotope& target_spec,
                   const Range& range,
                   const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;
  Size idx        = 0;

  for (const auto& [key, band] : bands) {
    const auto& ir = key.isot;

    if (ir != target_spec) {
      idx += band.lines.size();
      continue;
    }

    const Numeric Q   = PartitionFunctions::Q(T, ir);
    const Numeric drx = Constant::inv_sqrt_pi * atm[ir.spec];

    for (auto& line : band.lines) {
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
                      const AbsorptionBands& bands,
                      const Jacobian::LineTarget& target,
                      const Range& range,
                      const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;
  Size idx        = 0;

  for (const auto& [key, band] : bands) {
    if (target.type.band == key) {
      const auto& ir     = key.isot;
      auto& line         = band.lines[target.type.line];
      idx               += target.type.line;
      const Size nz      = line.z.size(line.qn, pol);
      const Numeric Q    = PartitionFunctions::Q(T, ir);
      const Numeric rx   = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];
      const Numeric s_   = line.s(T, Q);
      const Numeric ds_  = line.ds_df0(T, Q);
      const Numeric G    = line.ls.G(atm);
      const Numeric Y    = line.ls.Y(atm);
      const Numeric G0   = line.ls.G0(atm);
      const Numeric lmr  = 1 + G;
      const Numeric lmi  = -Y;

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
    idx += band.lines.size();
  }
}

void prepare_e0_deriv(MatrixView& mat,
                      const AtmPoint& atm,
                      const AbsorptionBands& bands,
                      const Jacobian::LineTarget& target,
                      const Range& range,
                      const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;
  Size idx        = 0;

  for (const auto& [key, band] : bands) {
    if (target.type.band == key) {
      const auto& ir     = key.isot;
      const Numeric rx   = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];
      const Numeric Q    = PartitionFunctions::Q(T, ir);
      auto& line         = band.lines[target.type.line];
      idx               += target.type.line;
      const Size nz      = line.z.size(line.qn, pol);
      const Numeric ds_  = line.ds_de0(T, Q);
      const Numeric G    = line.ls.G(atm);
      const Numeric Y    = line.ls.Y(atm);
      const Numeric lmr  = 1 + G;
      const Numeric lmi  = -Y;

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
    idx += band.lines.size();
  }
}

void prepare_a_deriv(MatrixView& mat,
                     const AtmPoint& atm,
                     const AbsorptionBands& bands,
                     const Jacobian::LineTarget& target,
                     const Range& range,
                     const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;
  Size idx        = 0;

  for (const auto& [key, band] : bands) {
    if (target.type.band == key) {
      const auto& ir     = key.isot;
      const Numeric Q    = PartitionFunctions::Q(T, ir);
      const Numeric rx   = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];
      auto& line         = band.lines[target.type.line];
      idx               += target.type.line;
      const Size nz      = line.z.size(line.qn, pol);
      const Numeric ds_  = line.ds_da(T, Q);
      const Numeric G    = line.ls.G(atm);
      const Numeric Y    = line.ls.Y(atm);
      const Numeric lmr  = 1 + G;
      const Numeric lmi  = -Y;

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
    idx += band.lines.size();
  }
}

void prepare_g0_deriv(MatrixView& mat,
                      const AtmPoint& atm,
                      const AbsorptionBands& bands,
                      const Jacobian::LineTarget& target,
                      const Range& range,
                      const ZeemanPolarization& pol) {
  Size idx = 0;
  for (const auto& [key, band] : bands) {
    if (target.type.band == key) {
      auto& line     = band.lines[target.type.line];
      idx           += target.type.line;
      const Size nz  = line.z.size(line.qn, pol);
      const Numeric dG0 =
          line.ls.dG0_dX(atm, target.type.band.isot.spec, target.type.ls_coeff);

      for (Size iz = 0; iz < nz; ++iz) {
        auto v  = mat[idx + iz];
        auto dv = v[range];
        dv[2]   = dG0 * v[1];
      }

      return;
    }
    idx += band.lines.size();
  }
}

void prepare_d0_deriv(MatrixView& mat,
                      const AtmPoint& atm,
                      const AbsorptionBands& bands,
                      const Jacobian::LineTarget& target,
                      const Range& range,
                      const ZeemanPolarization& pol) {
  Size idx = 0;
  for (const auto& [key, band] : bands) {
    if (target.type.band == key) {
      auto& line        = band.lines[target.type.line];
      idx              += target.type.line;
      const Size nz     = line.z.size(line.qn, pol);
      const Numeric G0  = line.ls.G0(atm);
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
    idx += band.lines.size();
  }
}

void prepare_dv_deriv(MatrixView& mat,
                      const AtmPoint& atm,
                      const AbsorptionBands& bands,
                      const Jacobian::LineTarget& target,
                      const Range& range,
                      const ZeemanPolarization& pol) {
  Size idx = 0;
  for (const auto& [key, band] : bands) {
    if (target.type.band == key) {
      auto& line        = band.lines[target.type.line];
      idx              += target.type.line;
      const Size nz     = line.z.size(line.qn, pol);
      const Numeric G0  = line.ls.G0(atm);
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
    idx += band.lines.size();
  }
}

void prepare_y_deriv(MatrixView& mat,
                     const AtmPoint& atm,
                     const AbsorptionBands& bands,
                     const Jacobian::LineTarget& target,
                     const Range& range,
                     const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;
  Size idx        = 0;

  for (const auto& [key, band] : bands) {
    if (target.type.band == key) {
      const auto& ir    = key.isot;
      const Numeric Q   = PartitionFunctions::Q(T, ir);
      const Numeric rx  = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];
      auto& line        = band.lines[target.type.line];
      idx              += target.type.line;
      const Size nz     = line.z.size(line.qn, pol);
      const Numeric s_  = line.s(T, Q);

      for (Size iz = 0; iz < nz; ++iz) {
        auto dv          = mat[idx + iz, range];
        const Numeric Sz = line.z.Strength(line.qn, pol, iz);
        const Numeric s  = Sz * s_;
        dv[4]            = -rx * s;
      }

      return;
    }
    idx += band.lines.size();
  }
}

void prepare_g_deriv(MatrixView& mat,
                     const AtmPoint& atm,
                     const AbsorptionBands& bands,
                     const Jacobian::LineTarget& target,
                     const Range& range,
                     const ZeemanPolarization& pol) {
  const Numeric T = atm.temperature;
  Size idx        = 0;

  for (const auto& [key, band] : bands) {
    if (target.type.band == key) {
      const auto& ir    = key.isot;
      const Numeric Q   = PartitionFunctions::Q(T, ir);
      const Numeric rx  = Constant::inv_sqrt_pi * atm[ir] * atm[ir.spec];
      auto& line        = band.lines[target.type.line];
      idx              += target.type.line;
      const Size nz     = line.z.size(line.qn, pol);
      const Numeric s_  = line.s(T, Q);

      for (Size iz = 0; iz < nz; ++iz) {
        auto dv          = mat[idx + iz, range];
        const Numeric Sz = line.z.Strength(line.qn, pol, iz);
        const Numeric s  = Sz * s_;
        dv[3]            = rx * s;
      }

      return;
    }
    idx += band.lines.size();
  }
}

void prepare_deriv(MatrixView& mat,
                   const AtmPoint& atm,
                   const AbsorptionBands& bands,
                   const AtmKey& key,
                   const Range& range,
                   const ZeemanPolarization& pol) {
  switch (key) {
    using enum AtmKey;
    case t:      prepare_t_deriv(mat, atm, bands, range, pol); break;
    case p:      prepare_p_deriv(mat, atm, bands, range, pol); break;
    case mag_u:  prepare_mu_deriv(mat, atm, bands, range, pol); break;
    case mag_v:  prepare_mv_deriv(mat, atm, bands, range, pol); break;
    case mag_w:  prepare_mw_deriv(mat, atm, bands, range, pol); break;
    case wind_u: [[fallthrough]];
    case wind_v: [[fallthrough]];
    case wind_w: break;
  }
}

void prepare_deriv(MatrixView&,
                   const AtmPoint&,
                   const AbsorptionBands&,
                   const ScatteringSpeciesProperty&,
                   const Range&,
                   const ZeemanPolarization&) {}

void prepare_deriv(MatrixView&,
                   const AtmPoint&,
                   const AbsorptionBands&,
                   const QuantumLevelIdentifier&,
                   const Range&,
                   const ZeemanPolarization&) {}

void prepare_deriv(MatrixView& mat,
                   const AtmPoint& atm,
                   const AbsorptionBands& bands,
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

void preparex(MatrixView mat,
              const AtmPoint& atm,
              const AbsorptionBands& bands,
              const JacobianTargets& jac_targets,
              const ZeemanPolarization& pol) {
  assert(static_cast<Size>(static_cast<Size>(mat.ncols())) ==
         jac_targets.target_count() * 5 + 5);

  preparex(mat, atm, bands, pol);

#pragma omp parallel for if (arts_omp_parallel())
  for (auto& target : jac_targets.atm) {
    std::visit(
        [&](const auto& x) {
          prepare_deriv(
              mat, atm, bands, x, {target.target_pos * 5 + 5, 5}, pol);
        },
        target.type);
  }

#pragma omp parallel for if (arts_omp_parallel())
  for (auto& target : jac_targets.line) {
    prepare_deriv(mat, atm, bands, target, {target.target_pos * 5 + 5, 5}, pol);
  }
}
}  // namespace

void str_scale(ComplexVectorView a,
               const AtmPoint& atm,
               const ConstVectorView& fs) {
  const Size n        = static_cast<Size>(a.ncols());
  constexpr Numeric c = Constant::c * Constant::c / (8 * Constant::pi);

  const Numeric T = atm.temperature;
  const Numeric P = atm.pressure;
  const Numeric N = number_density(P, T);

#pragma omp parallel for if (arts_omp_parallel())
  for (Size i = 0; i < n; i++) {
    const Numeric f  = fs[i];
    const Numeric r  = (Constant::h * f) / (Constant::k * T);
    a[i]            *= -N * f * std::expm1(-r) * c;
  }
}

void str_scale(ComplexMatrixView a,
               const AtmPoint& atm,
               const ConstVectorView& fs,
               const std::vector<bool>& df,
               const Size it) {
  const Size nf = fs.ncols();
  const Size nq = df.size();

  assert(static_cast<Size>(a.nrows()) == nf);
  assert(static_cast<Size>(a.ncols()) == (1 + nq));

  constexpr Numeric c = Constant::c * Constant::c / (8 * Constant::pi);

  const Numeric T  = atm.temperature;
  const Numeric P  = atm.pressure;
  const Numeric N  = number_density(P, T);
  const Numeric dN = dnumber_density_dt(P, T);

#pragma omp parallel for if (arts_omp_parallel())
  for (Size iv = 0; iv < nf; iv++) {
    const Numeric f = fs[iv];
    const Numeric r = (Constant::h * f) / (Constant::k * T);

    const Numeric em1 = std::expm1(-r);
    const Numeric emr = std::exp(-r);
    const Numeric n   = -N * f * em1 * c;

    a[iv, 0] *= n;

    if (it > 0) {
      a[iv, it] =
          a[iv, it] * n + a[iv, 0] * -f * (N * r * emr / T + dN * em1) * c;
    }

    for (Size iq = 0; iq < nq; iq++) {
      if (df[iq]) {
        a[iv, iq + 1] = a[iv, iq + 1] * n + a[iv, 0] * N * (r * emr - em1) * c;
      }
    }
  }
}

void sumup(ComplexVectorView a,
           const ConstMatrixView& mat,
           const ConstVectorView& f) {
  const Size n = mat.nrows();
  const Size m = f.ncols();
  assert(mat.ncols() >= 5);
  assert(static_cast<Size>(a.ncols()) == m);

#pragma omp parallel for if (arts_omp_parallel())
  for (Size j = 0; j < m; ++j) {
    const auto& f_ = f[j];
    auto& aj       = a[j];
    for (Size i = 0; i < n; ++i) {
      const auto&& v  = mat[i];
      aj             += expr(Faddeeva::w(z_(v, f_)), s_(v));
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

#pragma omp parallel for if (arts_omp_parallel())
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

void prepare(Matrix& mat,
             const AtmPoint& atm,
             const AbsorptionBands& bands,
             const ZeemanPolarization& pol) {
  mat.resize(lbl::count_zeeman_lines(bands, pol), 5);
  preparex(mat, atm, bands, pol);
}

void prepare(Matrix& mat,
             const AtmPoint& atm,
             const AbsorptionBands& bands,
             const JacobianTargets& jac_targets,
             const ZeemanPolarization& pol) {
  const Size n = lbl::count_zeeman_lines(bands, pol);
  mat.resize(n, 5 + jac_targets.target_count() * 5);
  mat = 0;
  preparex(mat, atm, bands, jac_targets, pol);
}
}  // namespace lbl::voigt::lte::matrix
