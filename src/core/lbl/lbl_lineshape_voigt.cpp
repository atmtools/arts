#include "lbl_lineshape_voigt.h"

#include <new_jacobian.h>

#include <cmath>
#include <limits>

#include "arts_constants.h"
#include "atm.h"
#include "debug.h"
#include "isotopologues.h"
#include "lbl_data.h"
#include "lbl_lineshape_model.h"
#include "lbl_temperature_model.h"
#include "lbl_zeeman.h"
#include "matpack_view.h"
#include "partfun.h"
#include "physics_funcs.h"
#include "quantum_numbers.h"

namespace lbl::voigt::lte {
Complex line_strength_calc(const SpeciesIsotopeRecord& spec,
                           const line& line,
                           const AtmPoint& atm) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const Complex lm{1 + line.ls.G(atm), -line.ls.Y(atm)};

  const auto n = atm[spec] * atm[spec.spec];

  return n * lm * s;
}

Complex dline_strength_calc_dVMR(const SpeciesIsotopeRecord& spec,
                                 const Species::Species target_spec,
                                 const line& line,
                                 const AtmPoint& atm) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const Complex lm{1 + line.ls.G(atm), -line.ls.Y(atm)};
  const Complex dlm{line.ls.dG_dVMR(atm, target_spec),
                    -line.ls.dY_dVMR(atm, target_spec)};

  const auto n = atm[spec] * atm[spec.spec];
  const auto dn = target_spec == spec.spec ? atm[spec] : 0.0;

  return (dn * lm + n * dlm) * s;
}

Complex dline_strength_calc_dT(const SpeciesIsotopeRecord& spec,
                               const line& line,
                               const AtmPoint& atm) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));
  const auto ds = line.ds_dT(atm.temperature,
                             PartitionFunctions::Q(atm.temperature, spec),
                             PartitionFunctions::dQdT(atm.temperature, spec));

  const Complex lm{1 + line.ls.G(atm), -line.ls.Y(atm)};
  const Complex dlm{line.ls.dG_dT(atm), -line.ls.dY_dT(atm)};

  const auto n = atm[spec] * atm[spec.spec];

  return n * (dlm * s + lm * ds);
}

Complex dline_strength_calc_disot(const SpeciesIsotopeRecord& spec,
                                  const SpeciesIsotopeRecord& target_isot,
                                  const line& line,
                                  const AtmPoint& atm) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const Complex lm{1 + line.ls.G(atm), -line.ls.Y(atm)};

  const auto dn = spec == target_isot ? atm[spec.spec] : 0.0;

  return dn * lm * s;
}

Complex line_strength_calc(const SpeciesIsotopeRecord& spec,
                           const line& line,
                           const AtmPoint& atm,
                           const Size ispec) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const auto& ls = line.ls.single_models[ispec];

  const Complex lm{1 + ls.G(line.ls.T0, atm.temperature, atm.pressure),
                   -ls.Y(line.ls.T0, atm.temperature, atm.pressure)};

  const auto n = atm[spec] * atm[spec.spec] * atm[ls.species];

  return n * lm * s;
}

Complex dline_strength_calc_dVMR(const SpeciesIsotopeRecord& spec,
                                 const Species::Species target_spec,
                                 const line& line,
                                 const AtmPoint& atm,
                                 const Size ispec) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const auto& ls = line.ls.single_models[ispec];

  const Complex lm{1 + ls.G(line.ls.T0, atm.temperature, atm.pressure),
                   -ls.Y(line.ls.T0, atm.temperature, atm.pressure)};

  const auto dn =
      (spec.spec == target_spec ? atm[spec.spec] * atm[ls.species] : 0.0) +
      (spec.spec == ls.species ? atm[spec] * atm[spec.spec] : 0.0);

  return dn * lm * s;
}

Complex dline_strength_calc_dT(const SpeciesIsotopeRecord& spec,
                               const line& line,
                               const AtmPoint& atm,
                               const Size ispec) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));
  const auto ds = line.ds_dT(atm.temperature,
                             PartitionFunctions::Q(atm.temperature, spec),
                             PartitionFunctions::dQdT(atm.temperature, spec));

  const auto& ls = line.ls.single_models[ispec];

  const Complex lm{1 + ls.G(line.ls.T0, atm.temperature, atm.pressure),
                   -ls.Y(line.ls.T0, atm.temperature, atm.pressure)};
  const Complex dlm{ls.dG_dT(line.ls.T0, atm.temperature, atm.pressure),
                    -ls.dY_dT(line.ls.T0, atm.temperature, atm.pressure)};

  const auto n = atm[spec] * atm[spec.spec] * atm[ls.species];

  return n * (dlm * s + lm * ds);
}

Complex dline_strength_calc_disot(const SpeciesIsotopeRecord& spec,
                                  const SpeciesIsotopeRecord& target_isot,
                                  const line& line,
                                  const AtmPoint& atm,
                                  const Size ispec) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const auto& ls = line.ls.single_models[ispec];

  const Complex lm{1 + ls.G(line.ls.T0, atm.temperature, atm.pressure),
                   -ls.Y(line.ls.T0, atm.temperature, atm.pressure)};

  const auto dn =
      (target_isot == spec) ? atm[spec.spec] * atm[ls.species] : 0.0;

  return dn * lm * s;
}

Numeric line_center_calc(const line& line, const AtmPoint& atm) {
  return line.f0 + line.ls.D0(atm) + line.ls.DV(atm);
}

Numeric dline_center_calc_dT(const line& line, const AtmPoint& atm) {
  return line.ls.dD0_dT(atm) + line.ls.dDV_dT(atm);
}

Numeric dline_center_calc_dVMR(const line& line,
                               const Species::Species spec,
                               const AtmPoint& atm) {
  return line.ls.dD0_dVMR(atm, spec) + line.ls.dDV_dVMR(atm, spec);
}

Numeric line_center_calc(const line& line, const AtmPoint& atm, Size ispec) {
  const auto& ls = line.ls.single_models[ispec];
  return line.f0 + ls.D0(line.ls.T0, atm.temperature, atm.pressure) +
         ls.DV(line.ls.T0, atm.temperature, atm.pressure);
}

Numeric dline_center_calc_dT(const line& line,
                             const AtmPoint& atm,
                             Size ispec) {
  const auto& ls = line.ls.single_models[ispec];
  return ls.dD0_dT(line.ls.T0, atm.temperature, atm.pressure) +
         ls.dDV_dT(line.ls.T0, atm.temperature, atm.pressure);
}

Numeric dline_center_calc_dVMR(const line& line,
                               const Species::Species spec,
                               const AtmPoint& atm,
                               Size ispec) {
  const auto& ls = line.ls.single_models[ispec];
  return ls.species == spec
             ? ls.D0(line.ls.T0, atm.temperature, atm.pressure) +
                   ls.DV(line.ls.T0, atm.temperature, atm.pressure)
             : 0;
}

Numeric scaled_gd(const Numeric T, const Numeric mass, const Numeric f0) {
  constexpr auto c = Constant::doppler_broadening_const_squared;
  return std::sqrt(c * T / mass) * f0;
}

single_shape::single_shape(const SpeciesIsotopeRecord& spec,
                           const line& line,
                           const AtmPoint& atm)
    : f0(line_center_calc(line, atm)),
      inv_gd(1.0 / scaled_gd(atm.temperature, spec.mass, f0)),
      z_imag(line.ls.G0(atm) * inv_gd),
      s(Constant::inv_sqrt_pi * inv_gd * line_strength_calc(spec, line, atm)) {}

single_shape::single_shape(const SpeciesIsotopeRecord& spec,
                           const line& line,
                           const AtmPoint& atm,
                           const Size ispec)
    : f0(line_center_calc(line, atm, ispec)),
      inv_gd(1.0 / scaled_gd(atm.temperature, spec.mass, f0)),
      z_imag(line.ls.single_models[ispec].G0(
                 line.ls.T0, atm.temperature, atm.pressure) *
             inv_gd),
      s(Constant::inv_sqrt_pi * inv_gd *
        line_strength_calc(spec, line, atm, ispec)) {}

void single_shape::as_zeeman(const line& line,
                             const Numeric H,
                             zeeman::pol pol,
                             Index iz) {
  s *= line.z.Strength(line.qn.val, pol, iz);
  f0 += H * line.z.Splitting(line.qn.val, pol, iz);
}

Numeric get_line_shape_dx(const Numeric T,
                          const Numeric P,
                          const line_shape::model& ls,
                          const ::Jacobian::LineTarget& line_target) {
  const auto& type = line_target.type;
  const auto& model = ls.single_models[type.spec];

#define TEMPCOEFF(COEFF, VAR)           \
  case temperature::coefficient::COEFF: \
    return model.d##VAR##_d##COEFF(ls.T0, T, P)

#define TEMPCOEFFS(VAR)                    \
  switch (type.ls_coeff) {                 \
    TEMPCOEFF(X0, VAR);                    \
    TEMPCOEFF(X1, VAR);                    \
    TEMPCOEFF(X2, VAR);                    \
    TEMPCOEFF(X3, VAR);                    \
    case temperature::coefficient::FINAL:; \
  }                                        \
  break

#define LSVAR(VAR)                \
  case line_shape::variable::VAR: \
    TEMPCOEFFS(VAR)

  switch (type.ls_var) {
    LSVAR(G0);
    LSVAR(D0);
    LSVAR(G2);
    LSVAR(D2);
    LSVAR(FVC);
    LSVAR(ETA);
    LSVAR(Y);
    LSVAR(G);
    LSVAR(DV);
    case line_shape::variable::FINAL:;
  }

#undef LSVAR
#undef TEMPCOEFFS
#undef TEMPCOEFF

  return 0;
}

struct ComputeData {
  std::vector<single_shape>
      lines{};  //! Line shapes; save for reuse, assume moved from
  std::vector<line_pos> pos{};  //! Save for reuse, size of line shapes

  ComplexVector cut{};   //! Size of line shapes
  ComplexVector dz{};    //! Size of line shapes
  ComplexVector ds{};    //! Size of line shapes
  ComplexVector dcut{};  //! Size of line shapes

  Vector scl{};            //! Size of frequency
  Vector dscl{};           //! Size of frequency
  ComplexVector shape{};   //! Size of frequency
  ComplexVector dshape{};  //! Size of frequency

  Propmat npm{};      //! The orientation of the polarization
  Propmat dnpm_du{};  //! The orientation of the polarization
  Propmat dnpm_dv{};  //! The orientation of the polarization
  Propmat dnpm_dw{};  //! The orientation of the polarization

  //! Sizes scl, dscl, shape, dshape.  Sets scl, npm, dnpm_du, dnpm_dv, dnpm_dw
  ComputeData(const Vector& f_grid,
              const AtmPoint& atm,
              const Vector2 los,
              const zeeman::pol pol)
      : scl(f_grid.size()),
        dscl(f_grid.size()),
        shape(f_grid.size()),
        dshape(f_grid.size()) {
    using Constant::h, Constant::k;
    std::transform(f_grid.begin(),
                   f_grid.end(),
                   scl.begin(),
                   [N = number_density(atm.pressure, atm.temperature),
                    T = atm.temperature](auto f) {
                     return -N * f * std::expm1(-h * f / (k * T));
                   });

    npm = zeeman::norm_view(pol, atm.mag, los);
    if (pol != zeeman::pol::no) {
      dnpm_du = zeeman::dnorm_view_du(pol, atm.mag, los);
      dnpm_dv = zeeman::dnorm_view_dv(pol, atm.mag, los);
      dnpm_dw = zeeman::dnorm_view_dw(pol, atm.mag, los);
    }
  }

  //! Sizes cut, dcut, dz, ds; sets shape
  void core_calc(const band_shape& shp, const band& bnd, const Vector& f_grid) {
    cut.resize(shp.size());
    dz.resize(shp.size());
    ds.resize(shp.size());
    dcut.resize(shp.size());

    if (bnd.cutoff != CutoffType::None) {
      shp(cut);
      std::transform(
          f_grid.begin(), f_grid.end(), shape.begin(), [this, &shp](Numeric f) {
            return shp(cut, f);
          });
    } else {
      std::transform(
          f_grid.begin(), f_grid.end(), shape.begin(), [&shp](Numeric f) {
            return shp(f);
          });
    }
  }

  //! Sets dshape and dscl and ds and dz
  void dt_core_calc(const SpeciesIsotopeRecord& spec,
                    const band_shape& shp,
                    const band& bnd,
                    const Vector& f_grid,
                    const AtmPoint& atm,
                    const zeeman::pol pol) {
    using Constant::h, Constant::k;

    std::transform(f_grid.begin(),
                   f_grid.end(),
                   dscl.begin(),
                   [N = number_density(atm.pressure, atm.temperature),
                    dN = dnumber_density_dt(atm.pressure, atm.temperature),
                    T = atm.temperature](auto f) {
                     return -f *
                            (N * f * h * exp(-h * f / (k * T)) / (T * T * k) +
                             dN * std::expm1(-h * f / (k * T)));
                   });

    const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
    for (Size i = 0; i < pos.size(); i++) {
      const bool single_line = pos[i].spec == std::numeric_limits<Size>::max();
      const Numeric& inv_gd = shp.lines[i].inv_gd;
      const auto& line = bnd.lines[pos[i].line];
      const Numeric dinv_gd_dT_ratio = -0.5 / atm.temperature;
      ds[i] = dinv_gd_dT_ratio * shp.lines[i].s +
              Constant::inv_sqrt_pi * inv_gd *
                  line.z.Strength(line.qn.val, pol, pos[i].iz) *
                  (single_line
                       ? dline_strength_calc_dT(spec, line, atm)
                       : dline_strength_calc_dT(spec, line, atm, pos[i].spec));
      dz[i] =
          dinv_gd_dT_ratio *
              Complex{-inv_gd * shp.lines[i].f0, shp.lines[i].z_imag} +
          inv_gd *
              Complex{-(single_line
                            ? dline_center_calc_dT(line, atm)
                            : dline_center_calc_dT(line, atm, pos[i].spec)) -
                          H * line.z.Splitting(line.qn.val, pol, pos[i].iz),
                      (single_line
                           ? line.ls.dG0_dT(atm)
                           : line.ls.single_models[pos[i].spec].dG0_dT(
                                 line.ls.T0, atm.temperature, atm.pressure))};
    }

    if (bnd.cutoff != CutoffType::None) {
      shp.dT(dcut, ds, dz);
      std::transform(
          f_grid.begin(),
          f_grid.end(),
          dshape.begin(),
          [this, &shp](Numeric f) { return shp.dT(dcut, ds, dz, f); });
    } else {
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dshape.begin(),
                     [this, &shp](Numeric f) { return shp.dT(ds, dz, f); });
    }
  }

  //! Sets dshape and dscl
  void df_core_calc(const band_shape& shp,
                    const band& bnd,
                    const Vector& f_grid,
                    const AtmPoint& atm) {
    using Constant::h, Constant::k;

    std::transform(f_grid.begin(),
                   f_grid.end(),
                   dscl.begin(),
                   [N = number_density(atm.pressure, atm.temperature),
                    T = atm.temperature](auto f) {
                     return N * (f * h * std::exp(-f * h / (T * k)) / (T * k) -
                                 std::expm1(-f * h / (T * k)));
                   });

    if (bnd.cutoff != CutoffType::None) {
      shp.df(dcut);
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dshape.begin(),
                     [this, &shp](Numeric f) { return shp.df(dcut, f); });
    } else {
      std::transform(
          f_grid.begin(), f_grid.end(), dshape.begin(), [&shp](Numeric f) {
            return shp.df(f);
          });
    }
  }

  //! Sets dshape and dz
  void dmag_u_core_calc(const band_shape& shp,
                        const band& bnd,
                        const Vector& f_grid,
                        const AtmPoint& atm,
                        const zeeman::pol pol) {
    const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
    const Numeric dH_dmag_u = atm.mag[0] / H;

    for (Size i = 0; i < pos.size(); i++) {
      const auto& line = bnd.lines[pos[i].line];
      dz[i] = -shp.lines[pos[i].line].inv_gd * dH_dmag_u *
              line.z.Splitting(line.qn.val, pol, pos[i].iz);
    }

    if (bnd.cutoff != CutoffType::None) {
      shp.dH(dcut, dz);
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dshape.begin(),
                     [this, &shp](Numeric f) { return shp.dH(dcut, dz, f); });
    } else {
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dshape.begin(),
                     [this, &shp](Numeric f) { return shp.dH(dz, f); });
    }
  }

  //! Sets dshape and dz
  void dmag_v_core_calc(const band_shape& shp,
                        const band& bnd,
                        const Vector& f_grid,
                        const AtmPoint& atm,
                        const zeeman::pol pol) {
    const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
    const Numeric dH_dmag_v = atm.mag[1] / H;

    for (Size i = 0; i < pos.size(); i++) {
      const auto& line = bnd.lines[pos[i].line];
      dz[i] = -shp.lines[pos[i].line].inv_gd * dH_dmag_v *
              line.z.Splitting(line.qn.val, pol, pos[i].iz);
    }

    if (bnd.cutoff != CutoffType::None) {
      shp.dH(dcut, dz);
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dshape.begin(),
                     [this, &shp](Numeric f) { return shp.dH(dcut, dz, f); });
    } else {
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dshape.begin(),
                     [this, &shp](Numeric f) { return shp.dH(dz, f); });
    }
  }

  //! Sets dshape and dz
  void dmag_w_core_calc(const band_shape& shp,
                        const band& bnd,
                        const Vector& f_grid,
                        const AtmPoint& atm,
                        const zeeman::pol pol) {
    const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
    const Numeric dH_dmag_w = atm.mag[2] / H;

    for (Size i = 0; i < pos.size(); i++) {
      const auto& line = bnd.lines[pos[i].line];
      dz[i] = -shp.lines[pos[i].line].inv_gd * dH_dmag_w *
              line.z.Splitting(line.qn.val, pol, pos[i].iz);
    }

    if (bnd.cutoff != CutoffType::None) {
      shp.dH(dcut, dz);
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dshape.begin(),
                     [this, &shp](Numeric f) { return shp.dH(dcut, dz, f); });
    } else {
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dshape.begin(),
                     [this, &shp](Numeric f) { return shp.dH(dz, f); });
    }
  }

  //! Sets ds and dz and dcut and dshape
  void dVMR_core_calc(const SpeciesIsotopeRecord& spec,
                      const band_shape& shp,
                      const band& bnd,
                      const Vector& f_grid,
                      const AtmPoint& atm,
                      const zeeman::pol pol,
                      const Species::Species target_spec) {
    using Constant::h, Constant::k;

    const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
    for (Size i = 0; i < pos.size(); i++) {
      const bool single_line = pos[i].spec == std::numeric_limits<Size>::max();
      const Numeric& inv_gd = shp.lines[i].inv_gd;
      const auto& line = bnd.lines[pos[i].line];
      ds[i] =
          Constant::inv_sqrt_pi * inv_gd *
          line.z.Strength(line.qn.val, pol, pos[i].iz) *
          (single_line ? dline_strength_calc_dVMR(spec, target_spec, line, atm)
                       : dline_strength_calc_dVMR(
                             spec, target_spec, line, atm, pos[i].spec));
      dz[i] =
          inv_gd *
          Complex{
              -(single_line ? dline_center_calc_dVMR(line, target_spec, atm)
                            : dline_center_calc_dVMR(
                                  line, target_spec, atm, pos[i].spec)) -
                  H * line.z.Splitting(line.qn.val, pol, pos[i].iz),
              (single_line
                   ? line.ls.dG0_dVMR(atm, target_spec)
                   : (line.ls.single_models[pos[i].spec].species == target_spec
                          ? line.ls.single_models[pos[i].spec].G0(
                                line.ls.T0, atm.temperature, atm.pressure)
                          : 0.0))};
    }

    if (bnd.cutoff != CutoffType::None) {
      shp.dVMR(dcut, ds, dz);
      std::transform(
          f_grid.begin(),
          f_grid.end(),
          dshape.begin(),
          [this, &shp](Numeric f) { return shp.dVMR(dcut, ds, dz, f); });
    } else {
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dshape.begin(),
                     [this, &shp](Numeric f) { return shp.dVMR(ds, dz, f); });
    }
  }
};

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const Vector& f_grid,
                        const SpeciesIsotopeRecord& spec,
                        const band_shape& shape,
                        const band& bnd,
                        const AtmPoint& atm,
                        const zeeman::pol pol,
                        const Atm::Key& key) {
  using enum Atm::Key;
  switch (key) {
    case t:
      com_data.dt_core_calc(spec, shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dscl[i] * com_data.shape[i] +
                                    com_data.scl[i] * com_data.dshape[i]);
      }
      break;
    case p:
      ARTS_USER_ERROR("Not implemented, pressure derivative");
      break;
    case mag_u:
      com_data.dmag_u_core_calc(shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dnpm_du,
                                com_data.scl[i] * com_data.shape[i],
                                com_data.scl[i] * com_data.dshape[i]);
      }
      break;
    case mag_v:
      com_data.dmag_v_core_calc(shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dnpm_dv,
                                com_data.scl[i] * com_data.shape[i],
                                com_data.scl[i] * com_data.dshape[i]);
      }
      break;
    case mag_w:
      com_data.dmag_w_core_calc(shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dnpm_dw,
                                com_data.scl[i] * com_data.shape[i],
                                com_data.scl[i] * com_data.dshape[i]);
      }
      break;
    case wind_u:
    case wind_v:
    case wind_w:
      com_data.df_core_calc(shape, bnd, f_grid, atm);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dscl[i] * com_data.shape[i] +
                                    com_data.scl[i] * com_data.dshape[i]);
      }
      break;
    case FINAL:
      break;
  }
}

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const Vector& f_grid,
                        const SpeciesIsotopeRecord& spec,
                        const band_shape&,
                        const band&,
                        const AtmPoint& atm,
                        const zeeman::pol,
                        const SpeciesIsotopeRecord& deriv_spec) {
  if (deriv_spec != spec) return;

  const Numeric isorat = atm[spec];

  ARTS_USER_ERROR_IF(
      isorat == 0,
      "Does not support 0 for isotopologue ratios (may be added upon request)")

  for (Index i = 0; i < f_grid.size(); i++) {
    dpm[i] += zeeman::scale(com_data.npm,
                            com_data.scl[i] * com_data.shape[i] / isorat);
  }
}

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const Vector& f_grid,
                        const SpeciesIsotopeRecord& spec,
                        const band_shape& shape,
                        const band& bnd,
                        const AtmPoint& atm,
                        const zeeman::pol pol,
                        const Species::Species& deriv_spec) {
  com_data.dVMR_core_calc(spec, shape, bnd, f_grid, atm, pol, deriv_spec);
  for (Index i = 0; i < f_grid.size(); i++) {
    dpm[i] += zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
  }
}

void compute_derivative(PropmatVectorView,
                        ComputeData&,
                        const Vector&,
                        const SpeciesIsotopeRecord&,
                        const band_shape&,
                        const band&,
                        const AtmPoint,
                        const zeeman::pol,
                        const line_key&) {}

void compute_derivative(PropmatVectorView,
                        ComputeData&,
                        const Vector&,
                        const SpeciesIsotopeRecord&,
                        const band_shape&,
                        const band&,
                        const AtmPoint,
                        const zeeman::pol,
                        const auto&) {}

void calculate(PropmatVector& pm,
               PropmatMatrix& dpm,
               const Vector& f_grid,
               const JacobianTargets& jacobian_targets,
               const band_key& bnd_qid,
               const band& bnd,
               const AtmPoint& atm,
               const Vector2 los,
               const zeeman::pol pol) {
  const SpeciesIsotopeRecord spec = bnd_qid.Isotopologue();
  const Numeric fmin = f_grid.front();
  const Numeric fmax = f_grid.back();

  const Index nf = f_grid.size();
  const Size njac = jacobian_targets.target_count();

  ARTS_ASSERT(njac == static_cast<Size>(dpm.nrows()) and nf == dpm.ncols())
  ARTS_ASSERT(nf == pm.nelem())

  //! FIXME: should be input so data sizes can be "stored"
  ComputeData com_data(f_grid, atm, los, pol);

  //! Reuse lines and positions if possible
  band_shape_helper(
      com_data.lines, com_data.pos, spec, bnd, atm, fmin, fmax, pol);

  //! Not const to save lines for reuse
  band_shape shape{std::move(com_data.lines), bnd.get_cutoff_frequency()};

  com_data.core_calc(shape, bnd, f_grid);

  for (Index i = 0; i < nf; i++) {
    pm[i] += zeeman::scale(com_data.npm, com_data.scl[i] * com_data.shape[i]);
  }

  for (auto& atm_target : jacobian_targets.atm()) {
    std::visit(
        [&](auto& target) {
          compute_derivative(dpm[atm_target.target_pos],
                             com_data,
                             f_grid,
                             spec,
                             shape,
                             bnd,
                             atm,
                             pol,
                             target);
        },
        atm_target.type);
  }
  
  for (auto& line_target : jacobian_targets.line()) {
      compute_derivative(dpm[line_target.target_pos],
                          com_data,
                          f_grid,
                          spec,
                          shape,
                          bnd,
                          atm,
                          pol,
                          line_target.type);
  }

  com_data.lines = std::move(shape.lines);
}
}  // namespace lbl::voigt::lte
