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

Complex single_shape::F(const Complex z_) const noexcept {
  return Faddeeva::w(z_);  // FIXME: Should factor be part of s?
}

Complex single_shape::F(const Numeric f) const noexcept { return F(z(f)); }

Complex single_shape::operator()(const Numeric f) const noexcept {
  return s * F(f);
}

Complex single_shape::dF(const Numeric f) const noexcept {
  Complex z_ = z(f);
  return dF(z_, F(z_));
}

single_shape::zFdF single_shape::all(const Numeric f) const noexcept {
  zFdF out;
  out.z = z(f);
  out.F = F(out.z);
  out.dF = dF(out.z, out.F);
  return out;
}

Complex single_shape::df(const Numeric f) const noexcept { return s * dF(f); }

Complex single_shape::df0(const Complex ds_df0,
                          const Complex dz_df0,
                          const Numeric f) const noexcept {
  const auto [z_, F_, dF_] = all(f);
  return ds_df0 * F_ + s * dz_df0 * dF_;
}

Complex single_shape::dDV(const Complex dz_dDV,
                          const Numeric f) const noexcept {
  return s * dz_dDV * dF(f);
}

Complex single_shape::dD0(const Complex dz_dD0,
                          const Numeric f) const noexcept {
  return s * dz_dD0 * dF(f);
}

Complex single_shape::dG0(const Complex dz_dG0,
                          const Numeric f) const noexcept {
  return s * dz_dG0 * dF(f);
}

Complex single_shape::dH(const Complex dz_dH, const Numeric f) const noexcept {
  return s * dz_dH * dF(f);
}

Complex single_shape::dVMR(const Complex ds_dVMR,
                           const Complex dz_dVMR,
                           const Numeric f) const noexcept {
  const auto [z_, F_, dF_] = all(f);
  return ds_dVMR * F_ + dz_dVMR * dF_;
}

Complex single_shape::dT(const Complex ds_dT,
                         const Complex dz_dT,
                         const Numeric f) const noexcept {
  const auto [z_, F_, dF_] = all(f);
  return ds_dT * F_ + dz_dT * dF_;
  //FIXME: invGD factor of F_ is missing
}

Complex single_shape::da(const Complex ds_da, const Numeric f) const noexcept {
  return ds_da * F(f);
}

Complex single_shape::de0(const Complex ds_de0,
                          const Numeric f) const noexcept {
  return ds_de0 * F(f);
}

Complex single_shape::dG(const Complex ds_dG, const Numeric f) const noexcept {
  return ds_dG * F(f);
}

Complex single_shape::dY(const Complex ds_dY, const Numeric f) const noexcept {
  return ds_dY * F(f);
}

Size count_lines(const band& bnd, const zeeman::pol type) {
  return std::transform_reduce(
      bnd.begin(), bnd.end(), Index{}, std::plus<>{}, [type](auto& line) {
        const Index factor =
            line.ls.one_by_one ? line.ls.single_models.size() : 1;
        return factor * line.z.size(line.qn.val, type);
      });
}

void zeeman_set_back(std::vector<single_shape>& lines,
                            std::vector<line_pos>& pos,
                            const single_shape& s,
                            const line& line,
                            const Numeric H,
                            const Size spec,
                            const zeeman::pol pol,
                            Size& last_single_shape_pos) {
  const auto line_nr = static_cast<Size>(pos.size() ? pos.back().line + 1 : 0);

  if (pol == zeeman::pol::no) {
    lines[last_single_shape_pos] = s;
    pos[last_single_shape_pos] = {line_nr, spec};
    ++last_single_shape_pos;
  } else {
    const auto nz = line.z.size(line.qn.val, pol);
    for (Index iz = 0; iz < nz; iz++) {
      lines[last_single_shape_pos] = s;
      lines[last_single_shape_pos].as_zeeman(line, H, pol, iz);
      pos[last_single_shape_pos] = {line_nr, spec, static_cast<Size>(iz)};
      ++last_single_shape_pos;
    }
  }
}

void lines_set(std::vector<single_shape>& lines,
                      std::vector<line_pos>& pos,
                      const SpeciesIsotopeRecord& spec,
                      const line& line,
                      const AtmPoint& atm,
                      const zeeman::pol pol,
                      Size& last_single_shape_pos) {
  const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  if (line.ls.one_by_one) {
    for (Size i = 0; i < line.ls.single_models.size(); ++i) {
      zeeman_set_back(lines,
                      pos,
                      {spec, line, atm, i},
                      line,
                      H,
                      i,
                      pol,
                      last_single_shape_pos);
    }
  } else {
    zeeman_set_back(lines,
                    pos,
                    {spec, line, atm},
                    line,
                    H,
                    std::numeric_limits<Size>::max(),
                    pol,
                    last_single_shape_pos);
  }
}

 void band_shape_helper(std::vector<single_shape>& lines,
                              std::vector<line_pos>& pos,
                              const SpeciesIsotopeRecord& spec,
                              const band& bnd,
                              const AtmPoint& atm,
                              const Numeric fmin,
                              const Numeric fmax,
                              const zeeman::pol pol) {
  lines.resize(count_lines(bnd, pol));
  pos.resize(lines.size());

  Size i{0};

  using enum CutoffType;
  switch (bnd.cutoff) {
    case None:
      for (auto& line : bnd) {
        lines_set(lines, pos, spec, line, atm, pol, i);
      }
      break;
    case Freq: {
      const auto by_range =
          std::views::filter([fmin, fmax, c = bnd.cutoff_value](auto& line) {
            return (line.f0 - c) > fmin and (line.f0 + c) < fmax;
          });
      for (auto& line : bnd | by_range) {
        lines_set(lines, pos, spec, line, atm, pol, i);
      }
    } break;
    case FINAL:
      ARTS_USER_ERROR("Bad state")
  }

  bubble_sort_by(
      [&](const Size l1, const Size l2) { return lines[l1].f0 < lines[l2].f0; },
      lines,
      pos);
}

struct ComputeData {
  std::vector<single_shape>
      lines{};  //! Line shapes; save for reuse, assume moved from
  std::vector<line_pos> pos{};  //! Save for reuse, size of line shapes
  std::vector<Size>
      filter;  //! Filter for line parameters; resized all the time but reserves size of line shapes

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
    filter.reserve(shp.size());

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

  void set_filter(const line_key& key, bool check_spec) {
    filter.resize(0);

    if (not check_spec) {
      for (Size i = 0; i < pos.size(); i++) {
        if (pos[i].line == key.line) filter.push_back(i);
      }
    } else {
      for (Size i = 0; i < pos.size(); i++) {
        if (pos[i].line == key.line and pos[i].spec == key.spec)
          filter.push_back(i);
      }
    }
  }

  //! Sets dshape and ds and dz and dcut and dshape
  void df0_core_calc(const band_shape& shp,
                     const band& bnd,
                     const Vector& f_grid,
                     const line_key& key) {
    using Constant::h, Constant::k;

    set_filter(key, false);

    for (Size i : filter) {
      const Numeric& inv_gd = shp.lines[i].inv_gd;
      const auto& line = bnd.lines[pos[i].line];
      const Numeric dinv_gd_df0_ratio = -1 / shp.lines[i].f0;
      const Numeric ds_df0_ratio = line.ds_df0_s_ratio();
      ds[i] =
          dinv_gd_df0_ratio * shp.lines[i].s + ds_df0_ratio * shp.lines[i].s;
      dz[i] = dinv_gd_df0_ratio *
                  Complex{-inv_gd * shp.lines[i].f0, shp.lines[i].z_imag} -
              inv_gd;
    }

    if (bnd.cutoff != CutoffType::None) {
      shp.df0(dcut, ds, dz, filter);
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.df0(dcut, ds, dz, f_grid[i], filter);
      }
    } else {
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.df0(ds, dz, f_grid[i], filter);
      }
    }
  }

  //! Sets dshape and ds and dcut and dshape
  void de0_core_calc(const band_shape& shp,
                     const band& bnd,
                     const Vector& f_grid,
                     const AtmPoint& atm,
                     const line_key& key) {
    using Constant::h, Constant::k;

    set_filter(key, false);

    for (Size i : filter) {
      const Numeric ds_de0_ratio =
          bnd.lines[pos[i].line].ds_de0_s_ratio(atm.temperature);
      ds[i] = ds_de0_ratio * shp.lines[i].s;
    }

    if (bnd.cutoff != CutoffType::None) {
      shp.de0(dcut, ds, filter);
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.de0(dcut, ds, f_grid[i], filter);
      }
    } else {
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.de0(ds, f_grid[i], filter);
      }
    }
  }

  //! Sets dshape and ds and dcut and dshape
  void da_core_calc(const band_shape& shp,
                    const band& bnd,
                    const Vector& f_grid,
                    const line_key& key) {
    using Constant::h, Constant::k;

    set_filter(key, false);

    for (Size i : filter) {
      const Numeric ds_da_ratio = 1.0 / bnd.lines[pos[i].line].a;
      ds[i] = ds_da_ratio * shp.lines[i].s;
    }

    if (bnd.cutoff != CutoffType::None) {
      shp.da(dcut, ds, filter);
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.da(dcut, ds, f_grid[i], filter);
      }
    } else {
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.da(ds, f_grid[i], filter);
      }
    }
  }

  //! Sets dshape and dz and dcut and dshape
  void dG0_core_calc(const band_shape& shp,
                     const band& bnd,
                     const Vector& f_grid,
                     const AtmPoint& atm,
                     const line_key& key) {
    using Constant::h, Constant::k;

    set_filter(key, true);

    using enum temperature::coefficient;
    switch (key.ls_coeff) {
      case X0:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          dz[i] = Complex(
              0,
              shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dG0_dX0(
                                        ls.T0, atm.temperature, atm.pressure));
        }
        break;
      case X1:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          dz[i] = Complex(
              0,
              shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dG0_dX1(
                                        ls.T0, atm.temperature, atm.pressure));
        }
        break;
      case X2:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          dz[i] = Complex(
              0,
              shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dG0_dX2(
                                        ls.T0, atm.temperature, atm.pressure));
        }
        break;
      case X3:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          dz[i] = Complex(
              0,
              shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dG0_dX3(
                                        ls.T0, atm.temperature, atm.pressure));
        }
        break;
      case FINAL:;
    }

    if (bnd.cutoff != CutoffType::None) {
      shp.dG0(dcut, dz, filter);
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.dG0(dcut, dz, f_grid[i], filter);
      }
    } else {
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.dG0(dz, f_grid[i], filter);
      }
    }
  }

  //! Sets dshape and dz and dcut and dshape
  void dD0_core_calc(const band_shape& shp,
                     const band& bnd,
                     const Vector& f_grid,
                     const AtmPoint& atm,
                     const line_key& key) {
    using Constant::h, Constant::k;

    set_filter(key, true);

    using enum temperature::coefficient;
    switch (key.ls_coeff) {
      case X0:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          dz[i] = -lines[i].inv_gd * ls.single_models[pos[i].spec].dD0_dX0(
                                         ls.T0, atm.temperature, atm.pressure);
        }
        break;
      case X1:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          dz[i] =
              -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dD0_dX1(
                                         ls.T0, atm.temperature, atm.pressure);
        }
        break;
      case X2:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          dz[i] =
              -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dD0_dX2(
                                         ls.T0, atm.temperature, atm.pressure);
        }
        break;
      case X3:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          dz[i] =
              -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dD0_dX3(
                                         ls.T0, atm.temperature, atm.pressure);
        }
        break;
      case FINAL:;
    }

    if (bnd.cutoff != CutoffType::None) {
      shp.dD0(dcut, dz, filter);
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.dD0(dcut, dz, f_grid[i], filter);
      }
    } else {
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.dD0(dz, f_grid[i], filter);
      }
    }
  }

  //! Sets dshape and ds and dcut and dshape
  void dY_core_calc(const band_shape& shp,
                    const band& bnd,
                    const Vector& f_grid,
                    const AtmPoint& atm,
                    const line_key& key) {
    using Constant::h, Constant::k;

    set_filter(key, true);

    using enum temperature::coefficient;
    switch (key.ls_coeff) {
      case X0:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          ds[i] = Complex{
              0,
              -lines[i].inv_gd * ls.single_models[pos[i].spec].dY_dX0(
                                     ls.T0, atm.temperature, atm.pressure)};
        }
        break;
      case X1:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          ds[i] = Complex{
              0,
              -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dY_dX1(
                                         ls.T0, atm.temperature, atm.pressure)};
        }
        break;
      case X2:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          ds[i] = Complex{
              0,
              -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dY_dX2(
                                         ls.T0, atm.temperature, atm.pressure)};
        }
        break;
      case X3:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          ds[i] = Complex{
              0,
              -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dY_dX3(
                                         ls.T0, atm.temperature, atm.pressure)};
        }
        break;
      case FINAL:;
    }

    if (bnd.cutoff != CutoffType::None) {
      shp.dY(dcut, dz, filter);
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.dY(dcut, ds, f_grid[i], filter);
      }
    } else {
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.dY(ds, f_grid[i], filter);
      }
    }
  }

  //! Sets dshape and ds and dcut and dshape
  void dG_core_calc(const band_shape& shp,
                    const band& bnd,
                    const Vector& f_grid,
                    const AtmPoint& atm,
                    const line_key& key) {
    using Constant::h, Constant::k;

    set_filter(key, true);

    using enum temperature::coefficient;
    switch (key.ls_coeff) {
      case X0:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          ds[i] = lines[i].inv_gd * ls.single_models[pos[i].spec].dG_dX0(
                                        ls.T0, atm.temperature, atm.pressure);
        }
        break;
      case X1:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          ds[i] =
              shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dG_dX1(
                                        ls.T0, atm.temperature, atm.pressure);
        }
        break;
      case X2:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          ds[i] =
              shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dG_dX2(
                                        ls.T0, atm.temperature, atm.pressure);
        }
        break;
      case X3:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          ds[i] =
              shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dG_dX3(
                                        ls.T0, atm.temperature, atm.pressure);
        }
        break;
      case FINAL:;
    }

    if (bnd.cutoff != CutoffType::None) {
      shp.dG(dcut, dz, filter);
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.dG(dcut, ds, f_grid[i], filter);
      }
    } else {
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.dG(ds, f_grid[i], filter);
      }
    }
  }

  //! Sets dshape and dz and dcut and dshape
  void dDV_core_calc(const band_shape& shp,
                     const band& bnd,
                     const Vector& f_grid,
                     const AtmPoint& atm,
                     const line_key& key) {
    using Constant::h, Constant::k;

    set_filter(key, true);

    using enum temperature::coefficient;
    switch (key.ls_coeff) {
      case X0:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          dz[i] = -lines[i].inv_gd * ls.single_models[pos[i].spec].dDV_dX0(
                                         ls.T0, atm.temperature, atm.pressure);
        }
        break;
      case X1:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          dz[i] =
              -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dDV_dX1(
                                         ls.T0, atm.temperature, atm.pressure);
        }
        break;
      case X2:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          dz[i] =
              -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dDV_dX2(
                                         ls.T0, atm.temperature, atm.pressure);
        }
        break;
      case X3:
        for (Size i : filter) {
          const auto& ls = bnd.lines[pos[i].line].ls;
          dz[i] =
              -shp.lines[i].inv_gd * ls.single_models[pos[i].spec].dDV_dX3(
                                         ls.T0, atm.temperature, atm.pressure);
        }
        break;
      case FINAL:;
    }

    if (bnd.cutoff != CutoffType::None) {
      shp.dDV(dcut, dz, filter);
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.dDV(dcut, dz, f_grid[i], filter);
      }
    } else {
      for (Index i = 0; i < f_grid.size(); i++) {
        dshape[i] = shp.dDV(dz, f_grid[i], filter);
      }
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

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const Vector& f_grid,
                        const band_shape& shape,
                        const band& bnd,
                        const AtmPoint& atm,
                        const zeeman::pol,
                        const line_key& deriv) {
  switch (deriv.var) {
    case variable::f0:
      com_data.df0_core_calc(shape, bnd, f_grid, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case variable::e0:
      com_data.de0_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case variable::a:
      com_data.da_core_calc(shape, bnd, f_grid, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case variable::FINAL:
      break;
  }

  switch (deriv.ls_var) {
    case line_shape::variable::G0:
      com_data.dG0_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case line_shape::variable::D0:
      com_data.dD0_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case line_shape::variable::G2:
      return;
    case line_shape::variable::D2:
      return;
    case line_shape::variable::FVC:
      return;
    case line_shape::variable::ETA:
      return;
    case line_shape::variable::Y:
      com_data.dY_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case line_shape::variable::G:
      com_data.dG_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case line_shape::variable::DV:
      com_data.dDV_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case line_shape::variable::FINAL:
      break;
  }
}

void compute_derivative(PropmatVectorView,
                        ComputeData&,
                        const Vector&,
                        const SpeciesIsotopeRecord&,
                        const band_shape&,
                        const band&,
                        const AtmPoint&,
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
    if (line_target.type.band == bnd_qid) {
      compute_derivative(dpm[line_target.target_pos],
                         com_data,
                         f_grid,
                         shape,
                         bnd,
                         atm,
                         pol,
                         line_target.type);
    }
  }

  com_data.lines = std::move(shape.lines);
}
}  // namespace lbl::voigt::lte
