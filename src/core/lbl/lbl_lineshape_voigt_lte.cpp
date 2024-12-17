#include "lbl_lineshape_voigt_lte.h"

#include <atm.h>
#include <jacobian.h>
#include <partfun.h>
#include <physics_funcs.h>
#include <sorting.h>

#include <Faddeeva/Faddeeva.hh>
#include <cmath>
#include <limits>
#include <numeric>

#include "lbl_data.h"
#include "lbl_zeeman.h"

namespace lbl::voigt::lte {
Complex line_strength_calc(const Numeric inv_gd,
                           const SpeciesIsotope& spec,
                           const line& line,
                           const AtmPoint& atm) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));
  const Numeric G = line.ls.G(atm);
  const Numeric Y = line.ls.Y(atm);

  const Complex lm{1 + G, -Y};
  const Numeric r = atm[spec];
  const Numeric x = atm[spec.spec];

  return Constant::inv_sqrt_pi * inv_gd * r * x * lm * s;
}

Complex dline_strength_calc_dY(const Numeric dY,
                               const Numeric inv_gd,
                               const SpeciesIsotope& spec,
                               const line& line,
                               const AtmPoint& atm) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const Numeric r = atm[spec];
  const Numeric x = atm[spec.spec];

  return Constant::inv_sqrt_pi * inv_gd * r * x * Complex(0, -dY) * s;
}

Complex dline_strength_calc_dG(const Numeric dG,
                               const Numeric inv_gd,
                               const SpeciesIsotope& spec,
                               const line& line,
                               const AtmPoint& atm) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const Numeric r = atm[spec];
  const Numeric x = atm[spec.spec];

  return Constant::inv_sqrt_pi * inv_gd * r * x * dG * s;
}

Complex dline_strength_calc_df0(const Numeric f0,
                                const Numeric inv_gd,
                                const SpeciesIsotope& spec,
                                const line& line,
                                const AtmPoint& atm) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));
  const auto ds = line.ds_df0_s_ratio() * s;

  const Numeric G = line.ls.G(atm);
  const Numeric Y = line.ls.Y(atm);

  const Complex lm{1 + G, -Y};

  const Numeric r = atm[spec];
  const Numeric x = atm[spec.spec];

  return Constant::inv_sqrt_pi * inv_gd * r * x * (f0 * ds - s) * lm / f0;
}

Complex dline_strength_calc_dVMR(const Numeric inv_gd,
                                 const Numeric f0,
                                 const SpeciesIsotope& spec,
                                 const SpeciesEnum target_spec,
                                 const line& line,
                                 const AtmPoint& atm) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const Numeric G   = line.ls.G(atm);
  const Numeric Y   = line.ls.Y(atm);
  const Numeric dG  = line.ls.dG_dVMR(atm, target_spec);
  const Numeric dY  = line.ls.dY_dVMR(atm, target_spec);
  const Numeric dD0 = line.ls.dD0_dVMR(atm, target_spec);
  const Numeric dDV = line.ls.dDV_dVMR(atm, target_spec);

  const Numeric df0 = dD0 + dDV;
  const Complex lm{1 + G, -Y};
  const Complex dlm = {dG, -dY};
  const Numeric r   = atm[spec];
  const Numeric x   = atm[spec.spec];

  if (target_spec == spec.spec) {
    return -Constant::inv_sqrt_pi * inv_gd * r * s *
           (x * (df0 / f0) * lm - (x * dlm + lm));
  }

  return -Constant::inv_sqrt_pi * inv_gd * r * s * x * ((df0 / f0) * lm - dlm);
}

Complex dline_strength_calc_dT(const Numeric inv_gd,
                               const Numeric f0,
                               const SpeciesIsotope& spec,
                               const line& line,
                               const AtmPoint& atm) {
  const Numeric T = atm.temperature;
  const auto s    = line.s(T, PartitionFunctions::Q(T, spec));
  const auto ds   = line.ds_dT(
      T, PartitionFunctions::Q(T, spec), PartitionFunctions::dQdT(T, spec));

  const Numeric G   = line.ls.G(atm);
  const Numeric Y   = line.ls.Y(atm);
  const Numeric dG  = line.ls.dG_dT(atm);
  const Numeric dY  = line.ls.dY_dT(atm);
  const Numeric dD0 = line.ls.dD0_dT(atm);
  const Numeric dDV = line.ls.dDV_dT(atm);

  const Numeric df0 = dD0 + dDV;
  const Complex lm{1 + G, -Y};
  const Complex dlm = {dG, -dY};
  const Numeric r   = atm[spec];
  const Numeric x   = atm[spec.spec];

  return Constant::inv_sqrt_pi * inv_gd * r * x *
         (2 * T * (dlm * s + lm * ds) * f0 - 2 * T * df0 * lm * s -
          f0 * lm * s) /
         (2 * T * f0);
}

Complex line_strength_calc(const Numeric inv_gd,
                           const SpeciesIsotope& spec,
                           const line& line,
                           const AtmPoint& atm,
                           const Size ispec) {
  const auto& ls   = line.ls.single_models[ispec];
  const Numeric T0 = line.ls.T0;
  const Numeric T  = atm.temperature;
  const Numeric P  = atm.pressure;
  const Numeric x  = atm[spec.spec];
  const Numeric r  = atm[spec];
  const Numeric v  = ls.species == SpeciesEnum::Bath
                         ? 1 - std::transform_reduce(
                                  line.ls.single_models.begin(),
                                  line.ls.single_models.end() - 1,
                                  0.0,
                                  std::plus<>{},
                                  [&atm](auto& s) { return atm[s.species]; })
                         : atm[ls.species];

  const auto s = line.s(T, PartitionFunctions::Q(T, spec));

  const auto G = ls.G(T0, T, P);
  const auto Y = ls.Y(T0, T, P);
  const Complex lm{1 + G, -Y};

  return Constant::inv_sqrt_pi * inv_gd * x * r * v * lm * s;
}

Complex dline_strength_calc_dG(const Numeric dG,
                               const Numeric inv_gd,
                               const SpeciesIsotope& spec,
                               const line& line,
                               const AtmPoint& atm,
                               const Size ispec) {
  const auto& ls  = line.ls.single_models[ispec];
  const Numeric T = atm.temperature;
  const Numeric x = atm[spec.spec];
  const Numeric r = atm[spec];
  const Numeric v = ls.species == SpeciesEnum::Bath
                        ? 1 - std::transform_reduce(
                                  line.ls.single_models.begin(),
                                  line.ls.single_models.end() - 1,
                                  0.0,
                                  std::plus<>{},
                                  [&atm](auto& s) { return atm[s.species]; })
                        : atm[ls.species];

  const auto s = line.s(T, PartitionFunctions::Q(T, spec));

  const Numeric dlm{dG};

  return Constant::inv_sqrt_pi * inv_gd * x * r * v * dlm * s;
}

Complex dline_strength_calc_dY(const Numeric dY,
                               const Numeric inv_gd,
                               const SpeciesIsotope& spec,
                               const line& line,
                               const AtmPoint& atm,
                               const Size ispec) {
  const auto& ls  = line.ls.single_models[ispec];
  const Numeric T = atm.temperature;
  const Numeric x = atm[spec.spec];
  const Numeric r = atm[spec];
  const Numeric v = ls.species == SpeciesEnum::Bath
                        ? 1 - std::transform_reduce(
                                  line.ls.single_models.begin(),
                                  line.ls.single_models.end() - 1,
                                  0.0,
                                  std::plus<>{},
                                  [&atm](auto& s) { return atm[s.species]; })
                        : atm[ls.species];

  const auto s = line.s(T, PartitionFunctions::Q(T, spec));

  const Complex dlm{0, -dY};

  return Constant::inv_sqrt_pi * inv_gd * x * r * v * dlm * s;
}

Complex dline_strength_calc_df0(const Numeric f0,
                                const Numeric inv_gd,
                                const SpeciesIsotope& spec,
                                const line& line,
                                const AtmPoint& atm,
                                const Size ispec) {
  const auto& ls   = line.ls.single_models[ispec];
  const Numeric T0 = line.ls.T0;
  const Numeric T  = atm.temperature;
  const Numeric P  = atm.pressure;
  const Numeric x  = atm[spec.spec];
  const Numeric r  = atm[spec];
  const Numeric v  = ls.species == SpeciesEnum::Bath
                         ? 1 - std::transform_reduce(
                                  line.ls.single_models.begin(),
                                  line.ls.single_models.end() - 1,
                                  0.0,
                                  std::plus<>{},
                                  [&atm](auto& s) { return atm[s.species]; })
                         : atm[ls.species];
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));
  const auto ds = line.ds_df0_s_ratio() * s;

  const auto G = ls.G(T0, T, P);
  const auto Y = ls.Y(T0, T, P);
  const Complex lm{1 + G, -Y};

  return Constant::inv_sqrt_pi * inv_gd * r * x * v * (f0 * ds - s) * lm / f0;
}

Complex dline_strength_calc_dT(const Numeric f0,
                               const Numeric inv_gd,
                               const SpeciesIsotope& spec,
                               const line& line,
                               const AtmPoint& atm,
                               const Size ispec) {
  const auto& ls   = line.ls.single_models[ispec];
  const Numeric T0 = line.ls.T0;
  const Numeric T  = atm.temperature;
  const Numeric P  = atm.pressure;
  const Numeric x  = atm[spec.spec];
  const Numeric r  = atm[spec];
  const Numeric v  = ls.species == SpeciesEnum::Bath
                         ? 1 - std::transform_reduce(
                                  line.ls.single_models.begin(),
                                  line.ls.single_models.end() - 1,
                                  0.0,
                                  std::plus<>{},
                                  [&atm](auto& s) { return atm[s.species]; })
                         : atm[ls.species];

  const auto s  = line.s(T, PartitionFunctions::Q(T, spec));
  const auto ds = line.ds_dT(
      T, PartitionFunctions::Q(T, spec), PartitionFunctions::dQdT(T, spec));

  const Numeric G   = ls.G(T0, T, P);
  const Numeric Y   = ls.Y(T0, T, P);
  const Numeric dG  = ls.dG_dT(T0, T, P);
  const Numeric dY  = ls.dY_dT(T0, T, P);
  const Numeric dD0 = ls.dD0_dT(T0, T, P);
  const Numeric dDV = ls.dDV_dT(T0, T, P);

  const Numeric df0 = dD0 + dDV;
  const Complex lm{1 + G, -Y};
  const Complex dlm = {dG, -dY};

  return Constant::inv_sqrt_pi * inv_gd * r * x * v *
         (2 * T * (dlm * s + lm * ds) * f0 - 2 * T * df0 * lm * s -
          f0 * lm * s) /
         (2 * T * f0);
}

Numeric line_center_calc(const line& line, const AtmPoint& atm) {
  return line.f0 + line.ls.D0(atm) + line.ls.DV(atm);
}

Numeric dline_center_calc_dT(const line& line, const AtmPoint& atm) {
  return line.ls.dD0_dT(atm) + line.ls.dDV_dT(atm);
}

Numeric dline_center_calc_dVMR(const line& line,
                               const SpeciesEnum spec,
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
                               const SpeciesEnum spec,
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

//! Should only live in CC-file since it holds references
struct single_shape_builder {
  const SpeciesIsotope& spec;
  const line& ln;
  const AtmPoint& atm;
  Numeric f0;
  Numeric scaled_gd_part;
  Numeric G0;
  Size ispec{std::numeric_limits<Size>::max()};

  single_shape_builder(const SpeciesIsotope& s,
                       const line& l,
                       const AtmPoint& a)
      : spec(s),
        ln(l),
        atm(a),
        f0(line_center_calc(ln, atm)),
        scaled_gd_part(std::sqrt(Constant::doppler_broadening_const_squared *
                                 atm.temperature / s.mass)),
        G0(ln.ls.G0(atm)) {}

  single_shape_builder(const SpeciesIsotope& s,
                       const line& l,
                       const AtmPoint& a,
                       const Size is)
      : spec(s),
        ln(l),
        atm(a),
        f0(line_center_calc(ln, atm, is)),
        scaled_gd_part(std::sqrt(Constant::doppler_broadening_const_squared *
                                 atm.temperature / s.mass)),
        G0(ln.ls.single_models[is].G0(ln.ls.T0, atm.temperature, atm.pressure)),
        ispec(is) {}

  [[nodiscard]] single_shape as_zeeman(const Numeric H,
                                       const zeeman::pol pol,
                                       const Size iz) const {
    single_shape s;
    s.f0     = f0 + H * ln.z.Splitting(ln.qn.val, pol, iz);
    s.inv_gd = 1.0 / (scaled_gd_part * f0);
    s.z_imag = G0 * s.inv_gd;
    s.s      = ln.z.Strength(ln.qn.val, pol, iz) *
          (ispec == std::numeric_limits<Size>::max()
               ? line_strength_calc(s.inv_gd, spec, ln, atm)
               : line_strength_calc(s.inv_gd, spec, ln, atm, ispec));
    return s;
  }

  operator single_shape() const {
    single_shape s;
    s.f0     = f0;
    s.inv_gd = 1.0 / (scaled_gd_part * f0);
    s.z_imag = G0 * s.inv_gd;
    s.s      = (ispec == std::numeric_limits<Size>::max()
                    ? line_strength_calc(s.inv_gd, spec, ln, atm)
                    : line_strength_calc(s.inv_gd, spec, ln, atm, ispec));
    return s;
  }
};

single_shape::single_shape(const SpeciesIsotope& spec,
                           const line& line,
                           const AtmPoint& atm,
                           const zeeman::pol pol,
                           const Index iz)
    : f0(line_center_calc(line, atm) +
         std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]) *
             line.z.Splitting(line.qn.val, pol, iz)),
      inv_gd(1.0 / scaled_gd(atm.temperature, spec.mass, f0)),
      z_imag(line.ls.G0(atm) * inv_gd),
      s(line.z.Strength(line.qn.val, pol, iz) *
        line_strength_calc(inv_gd, spec, line, atm)) {}

single_shape::single_shape(const SpeciesIsotope& spec,
                           const line& line,
                           const AtmPoint& atm,
                           const zeeman::pol pol,
                           const Index iz,
                           const Size ispec)
    : f0(line_center_calc(line, atm, ispec) +
         std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]) *
             line.z.Splitting(line.qn.val, pol, iz)),
      inv_gd(1.0 / scaled_gd(atm.temperature, spec.mass, f0)),
      z_imag(line.ls.single_models[ispec].G0(
                 line.ls.T0, atm.temperature, atm.pressure) *
             inv_gd),
      s(line.z.Strength(line.qn.val, pol, iz) *
        line_strength_calc(inv_gd, spec, line, atm, ispec)) {}

Complex single_shape::F(const Complex z_) { return Faddeeva::w(z_); }

Complex single_shape::F(const Numeric f) const { return F(z(f)); }

Complex single_shape::operator()(const Numeric f) const { return s * F(f); }

Complex single_shape::dF(const Numeric f) const {
  const Complex z_ = z(f);
  return dF(z_, F(z_));
}

Complex single_shape::dF(const Complex z_, const Complex F_) {
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
  const Complex dz{std::max(1e-4 * z_.real(), 1e-4),
                   std::max(1e-4 * z_.imag(), 1e-4)};
  const Complex F_2 = Faddeeva::w(z_ + dz);
  return (F_2 - F_) / dz;
}

single_shape::zFdF::zFdF(const Complex z_)
    : z{z_}, F{single_shape::F(z_)}, dF{single_shape::dF(z_, F)} {}

single_shape::zFdF single_shape::all(const Numeric f) const { return z(f); }

Complex single_shape::df(const Numeric f) const { return s * inv_gd * dF(f); }

Complex single_shape::df0(const Complex ds_df0,
                          const Complex dz_df0,
                          const Numeric dz_df0_fac,
                          const Numeric f) const {
  const auto [z_, F_, dF_] = all(f);
  return ds_df0 * F_ + s * (dz_df0 + dz_df0_fac * z_) * dF_;
}

Complex single_shape::dDV(const Complex ds_dDV,
                          const Complex dz_dDV,
                          const Numeric dz_dDV_fac,
                          const Numeric f) const {
  const auto [z_, F_, dF_] = all(f);
  return ds_dDV * F_ + s * (dz_dDV + dz_dDV_fac * z_) * dF_;
}

Complex single_shape::dD0(const Complex ds_dD0,
                          const Complex dz_dD0,
                          const Numeric dz_dD0_fac,
                          const Numeric f) const {
  const auto [z_, F_, dF_] = all(f);
  return ds_dD0 * F_ + s * (dz_dD0 + dz_dD0_fac * z_) * dF_;
}

Complex single_shape::dG0(const Complex dz_dG0, const Numeric f) const {
  return s * dz_dG0 * dF(f);
}

Complex single_shape::dH(const Complex dz_dH, const Numeric f) const {
  return s * dz_dH * dF(f);
}

Complex single_shape::dVMR(const Complex ds_dVMR,
                           const Complex dz_dVMR,
                           const Numeric dz_dVMR_fac,
                           const Numeric f) const {
  const auto [z_, F_, dF_] = all(f);
  return ds_dVMR * F_ + s * (dz_dVMR + dz_dVMR_fac * z_) * dF_;
}

Complex single_shape::dT(const Complex ds_dT,
                         const Complex dz_dT,
                         const Numeric dz_dT_fac,
                         const Numeric f) const {
  const auto [z_, F_, dF_] = all(f);
  return ds_dT * F_ + s * (dz_dT + dz_dT_fac * z_) * dF_;
}

Complex single_shape::da(const Complex ds_da, const Numeric f) const {
  return ds_da * F(f);
}

Complex single_shape::de0(const Complex ds_de0, const Numeric f) const {
  return ds_de0 * F(f);
}

Complex single_shape::dG(const Complex ds_dG, const Numeric f) const {
  return ds_dG * F(f);
}

Complex single_shape::dY(const Complex ds_dY, const Numeric f) const {
  return ds_dY * F(f);
}

Size count_lines(const band_data& bnd, const zeeman::pol type) {
  return std::transform_reduce(
      bnd.begin(), bnd.end(), Index{}, std::plus<>{}, [type](auto& line) {
        const Index factor =
            line.ls.one_by_one ? line.ls.single_models.size() : 1;
        return factor * line.z.size(line.qn.val, type);
      });
}

void zeeman_push_back(std::vector<single_shape>& lines,
                      std::vector<line_pos>& pos,
                      const single_shape_builder& s,
                      const line& line,
                      const AtmPoint& atm,
                      const zeeman::pol pol,
                      const Size ispec,
                      const Size iline) {
  if (pol == zeeman::pol::no) {
    lines.emplace_back(s);
    pos.emplace_back(line_pos{.line = iline, .spec = ispec});
  } else {
    const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
    const auto nz   = static_cast<Size>(line.z.size(line.qn.val, pol));
    for (Size iz = 0; iz < nz; iz++) {
      lines.emplace_back(s.as_zeeman(H, pol, iz));
      pos.emplace_back(line_pos{.line = iline, .spec = ispec, .iz = iz});

      if (lines.back().s == 0.0) {
        lines.pop_back();
        pos.pop_back();
      }
    }
  }
}

void lines_push_back(std::vector<single_shape>& lines,
                     std::vector<line_pos>& pos,
                     const SpeciesIsotope& spec,
                     const line& line,
                     const AtmPoint& atm,
                     const zeeman::pol pol,
                     const Size iline) {
  if (line.ls.one_by_one) {
    for (Size i = 0; i < line.ls.single_models.size(); ++i) {
      if ((line.z.on and pol != zeeman::pol::no) or
          (not line.z.on and pol == zeeman::pol::no)) {
        zeeman_push_back(lines,
                         pos,
                         single_shape_builder{spec, line, atm, i},
                         line,
                         atm,
                         pol,
                         i,
                         iline);
      }
    }
  } else {
    if ((line.z.on and pol != zeeman::pol::no) or
        (not line.z.on and pol == zeeman::pol::no)) {
      zeeman_push_back(lines,
                       pos,
                       single_shape_builder{spec, line, atm},
                       line,
                       atm,
                       pol,
                       std::numeric_limits<Size>::max(),
                       iline);
    }
  }
}

void band_shape_helper(std::vector<single_shape>& lines,
                       std::vector<line_pos>& pos,
                       const SpeciesIsotope& spec,
                       const band_data& bnd,
                       const AtmPoint& atm,
                       const Numeric fmin,
                       const Numeric fmax,
                       const zeeman::pol pol) {
  lines.resize(0);
  pos.resize(0);

  lines.reserve(count_lines(bnd, pol));
  pos.reserve(lines.capacity());

  using enum LineByLineCutoffType;
  switch (bnd.cutoff) {
    case None:
      for (Size iline = 0; iline < bnd.size(); iline++) {
        lines_push_back(lines, pos, spec, bnd.lines[iline], atm, pol, iline);
      }
      break;
    case ByLine: {
      auto [iline, active_lines] = bnd.active_lines(fmin, fmax);
      for (auto& line : active_lines) {
        lines_push_back(lines, pos, spec, line, atm, pol, iline++);
      }
    } break;
  }

  bubble_sort_by(
      [&](const Size l1, const Size l2) { return lines[l1].f0 > lines[l2].f0; },
      lines,
      pos);
}

band_shape::band_shape(std::vector<single_shape>&& ls, const Numeric cut)
    : lines(std::move(ls)), cutoff(cut) {}

Complex band_shape::operator()(const Numeric f) const {
  return std::transform_reduce(
      lines.begin(), lines.end(), Complex{}, std::plus<>{}, [f](auto& ls) {
        return ls(f);
      });
}

Complex band_shape::df(const Numeric f) const {
  return std::transform_reduce(
      lines.begin(), lines.end(), Complex{}, std::plus<>{}, [f](auto& ls) {
        return ls.df(f);
      });
}

Complex band_shape::dH(const ExhaustiveConstComplexVectorView& dz_dH,
                       const Numeric f) const {
  ARTS_ASSERT(static_cast<Size>(dz_dH.size()) == lines.size())

  return std::transform_reduce(lines.begin(),
                               lines.end(),
                               dz_dH.begin(),
                               Complex{},
                               std::plus<>{},
                               [f](auto& ls, auto& d) { return ls.dH(d, f); });
}

Complex band_shape::dT(const ExhaustiveConstComplexVectorView& ds_dT,
                       const ExhaustiveConstComplexVectorView& dz_dT,
                       const ExhaustiveConstVectorView& dz_dT_fac,
                       const Numeric f) const {
  ARTS_ASSERT(ds_dT.size() == dz_dT.size())
  ARTS_ASSERT(static_cast<Size>(ds_dT.size()) == lines.size())

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i = 0; i < lines.size(); ++i) {
    out += lines[i].dT(ds_dT[i], dz_dT[i], dz_dT_fac[i], f);
  }

  return out;
}

Complex band_shape::dVMR(const ExhaustiveConstComplexVectorView& ds_dVMR,
                         const ExhaustiveConstComplexVectorView& dz_dVMR,
                         const ExhaustiveConstVectorView& dz_dVMR_fac,
                         const Numeric f) const {
  ARTS_ASSERT(ds_dVMR.size() == dz_dVMR.size())
  ARTS_ASSERT(static_cast<Size>(ds_dVMR.size()) == lines.size())

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i = 0; i < lines.size(); ++i) {
    out += lines[i].dVMR(ds_dVMR[i], dz_dVMR[i], dz_dVMR_fac[i], f);
  }

  return out;
}

Complex band_shape::df0(const ExhaustiveConstComplexVectorView ds_df0,
                        const ExhaustiveConstComplexVectorView dz_df0,
                        const ExhaustiveConstVectorView dz_df0_fac,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].df0(ds_df0[i], dz_df0[i], dz_df0_fac[i], f);
  }

  return out;
}

Complex band_shape::da(const ExhaustiveConstComplexVectorView ds_da,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].da(ds_da[i], f);
  }

  return out;
}

Complex band_shape::de0(const ExhaustiveConstComplexVectorView ds_de0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].de0(ds_de0[i], f);
  }

  return out;
}

Complex band_shape::dDV(const ExhaustiveConstComplexVectorView ds_dDV,
                        const ExhaustiveConstComplexVectorView dz_dDV,
                        const ExhaustiveConstVectorView dz_dDV_fac,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dDV(ds_dDV[i], dz_dDV[i], dz_dDV_fac[i], f);
  }

  return out;
}

Complex band_shape::dD0(const ExhaustiveConstComplexVectorView ds_dD0,
                        const ExhaustiveConstComplexVectorView dz_dD0,
                        const ExhaustiveConstVectorView dz_dD0_fac,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dD0(ds_dD0[i], dz_dD0[i], dz_dD0_fac[i], f);
  }

  return out;
}

Complex band_shape::dG0(const ExhaustiveConstComplexVectorView dz_dG0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dG0(dz_dG0[i], f);
  }

  return out;
}

Complex band_shape::dY(const ExhaustiveConstComplexVectorView ds_dY,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dY(ds_dY[i], f);
  }

  return out;
}

Complex band_shape::dG(const ExhaustiveConstComplexVectorView ds_dG,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dG(ds_dG[i], f);
  }

  return out;
}

Complex band_shape::operator()(const ExhaustiveConstComplexVectorView& cut,
                               const Numeric f) const {
  const auto [s, cs] = frequency_spans(cutoff, f, lines, cut);
  return std::transform_reduce(s.begin(),
                               s.end(),
                               cs.begin(),
                               Complex{},
                               std::plus<>{},
                               [f](auto& ls, auto& c) { return ls(f) - c; });
}

void band_shape::operator()(ExhaustiveComplexVectorView cut) const {
  std::transform(
      lines.begin(),
      lines.end(),
      cut.begin(),
      [cutoff_freq = cutoff](auto& ls) { return ls(ls.f0 + cutoff_freq); });
}

Complex band_shape::df(const ExhaustiveConstComplexVectorView& cut,
                       const Numeric f) const {
  const auto [s, cs] = frequency_spans(cutoff, f, lines, cut);
  return std::transform_reduce(s.begin(),
                               s.end(),
                               cs.begin(),
                               Complex{},
                               std::plus<>{},
                               [f](auto& ls, auto& c) { return ls.df(f) - c; });
}

void band_shape::df(ExhaustiveComplexVectorView cut) const {
  std::transform(
      lines.begin(),
      lines.end(),
      cut.begin(),
      [cutoff_freq = cutoff](auto& ls) { return ls.df(ls.f0 + cutoff_freq); });
}

Complex band_shape::dH(const ExhaustiveConstComplexVectorView& cut,
                       const ExhaustiveConstComplexVectorView& dz_dH,
                       const Numeric f) const {
  ARTS_ASSERT(static_cast<Size>(dz_dH.size()) == lines.size())

  const auto [s, cs, dH] = frequency_spans(cutoff, f, lines, cut, dz_dH);

  Complex out{};  //! Fixme, use zip in C++ 23...
  for (Size i = 0; i < s.size(); ++i) {
    out += s[i].dH(dH[i], f) - cs[i];
  }

  return out;
}

void band_shape::dH(ExhaustiveComplexVectorView cut,
                    const ExhaustiveConstComplexVectorView& df0_dH) const {
  ARTS_ASSERT(static_cast<Size>(df0_dH.size()) == lines.size())

  std::transform(lines.begin(),
                 lines.end(),
                 df0_dH.begin(),
                 cut.begin(),
                 [cutoff_freq = cutoff](auto& ls, auto& d) {
                   return ls.dH(d, ls.f0 + cutoff_freq);
                 });
}

Complex band_shape::dT(const ExhaustiveConstComplexVectorView& cut,
                       const ExhaustiveConstComplexVectorView& ds_dT,
                       const ExhaustiveConstComplexVectorView& dz_dT,
                       const ExhaustiveConstVectorView& dz_dT_fac,
                       const Numeric f) const {
  ARTS_ASSERT(ds_dT.size() == dz_dT.size())
  ARTS_ASSERT(static_cast<Size>(ds_dT.size()) == lines.size())

  Complex out{};  //! Fixme, use zip in C++ 23...

  const auto [s, cs, ds, dz, dzf] =
      frequency_spans(cutoff, f, lines, cut, ds_dT, dz_dT, dz_dT_fac);

  for (Size i = 0; i < s.size(); ++i) {
    out += s[i].dT(ds[i], dz[i], dzf[i], f) - cs[i];
  }

  return out;
}

void band_shape::dT(ExhaustiveComplexVectorView cut,
                    const ExhaustiveConstComplexVectorView& ds_dT,
                    const ExhaustiveConstComplexVectorView& dz_dT,
                    const ExhaustiveConstVectorView& dz_dT_fac) const {
  ARTS_ASSERT(ds_dT.size() == dz_dT.size())
  ARTS_ASSERT(static_cast<Size>(ds_dT.size()) == lines.size())

  for (Size i = 0; i < lines.size(); ++i) {
    cut[i] =
        lines[i].dT(ds_dT[i], dz_dT[i], dz_dT_fac[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dVMR(const ExhaustiveConstComplexVectorView& cut,
                         const ExhaustiveConstComplexVectorView& ds_dVMR,
                         const ExhaustiveConstComplexVectorView& dz_dVMR,
                         const ExhaustiveConstVectorView& dz_dVMR_fac,
                         const Numeric f) const {
  ARTS_ASSERT(ds_dVMR.size() == dz_dVMR.size())
  ARTS_ASSERT(static_cast<Size>(ds_dVMR.size()) == lines.size())

  Complex out{};  //! Fixme, use zip in C++ 23...

  const auto [s, cs, ds, dz, dzf] =
      frequency_spans(cutoff, f, lines, cut, ds_dVMR, dz_dVMR, dz_dVMR_fac);

  for (Size i = 0; i < s.size(); ++i) {
    out += s[i].dVMR(ds[i], dz[i], dzf[i], f) - cs[i];
  }

  return out;
}

void band_shape::dVMR(ExhaustiveComplexVectorView cut,
                      const ExhaustiveConstComplexVectorView& ds_dVMR,
                      const ExhaustiveConstComplexVectorView& dz_dVMR,
                      const ExhaustiveConstVectorView& dz_dVMR_fac) const {
  ARTS_ASSERT(ds_dVMR.size() == dz_dVMR.size())
  ARTS_ASSERT(static_cast<Size>(ds_dVMR.size()) == lines.size())

  for (Size i = 0; i < lines.size(); ++i) {
    cut[i] = lines[i].dVMR(
        ds_dVMR[i], dz_dVMR[i], dz_dVMR_fac[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::df0(const ExhaustiveConstComplexVectorView& cut,
                        const ExhaustiveConstComplexVectorView ds_df0,
                        const ExhaustiveConstComplexVectorView dz_df0,
                        const ExhaustiveConstVectorView dz_df0_fac,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  const auto [s, cs, ds, dz, dzf] =
      frequency_spans(cutoff, f, lines, cut, ds_df0, dz_df0, dz_df0_fac);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].df0(ds[i], dz[i], dzf[i], f) - cs[i];
  }

  return out;
}

void band_shape::df0(ExhaustiveComplexVectorView cut,
                     const ExhaustiveConstComplexVectorView ds_df0,
                     const ExhaustiveConstComplexVectorView dz_df0,
                     const ExhaustiveConstVectorView dz_df0_fac,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] =
        lines[i].df0(ds_df0[i], dz_df0[i], dz_df0_fac[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::da(const ExhaustiveConstComplexVectorView& cut,
                       const ExhaustiveConstComplexVectorView ds_da,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_da);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].da(ds[i], f) - cs[i];
  }

  return out;
}

void band_shape::da(ExhaustiveComplexVectorView cut,
                    const ExhaustiveConstComplexVectorView ds_da,
                    const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].da(ds_da[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::de0(const ExhaustiveConstComplexVectorView& cut,
                        const ExhaustiveConstComplexVectorView ds_de0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_de0);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].de0(ds[i], f) - cs[i];
  }

  return out;
}

void band_shape::de0(ExhaustiveComplexVectorView cut,
                     const ExhaustiveConstComplexVectorView ds_de0,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].de0(ds_de0[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dDV(const ExhaustiveConstComplexVectorView& cut,
                        const ExhaustiveConstComplexVectorView ds_dDV,
                        const ExhaustiveConstComplexVectorView dz_dDV,
                        const ExhaustiveConstVectorView dz_dDV_fac,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  const auto [s, cs, ds, dz, dzf] =
      frequency_spans(cutoff, f, lines, cut, ds_dDV, dz_dDV, dz_dDV_fac);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].dDV(ds[i], dz[i], dzf[i], f) - cs[i];
  }

  return out;
}

void band_shape::dDV(ExhaustiveComplexVectorView cut,
                     const ExhaustiveConstComplexVectorView ds_dDV,
                     const ExhaustiveConstComplexVectorView dz_dDV,
                     const ExhaustiveConstVectorView dz_dDV_fac,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] =
        lines[i].dDV(ds_dDV[i], dz_dDV[i], dz_dDV_fac[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dD0(const ExhaustiveConstComplexVectorView& cut,
                        const ExhaustiveConstComplexVectorView ds_dD0,
                        const ExhaustiveConstComplexVectorView dz_dD0,
                        const ExhaustiveConstVectorView dz_dD0_fac,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  const auto [s, cs, ds, dz, dzf] =
      frequency_spans(cutoff, f, lines, cut, ds_dD0, dz_dD0, dz_dD0_fac);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].dD0(ds[i], dz[i], dzf[i], f) - cs[i];
  }

  return out;
}

void band_shape::dD0(ExhaustiveComplexVectorView cut,
                     const ExhaustiveConstComplexVectorView ds_dD0,
                     const ExhaustiveConstComplexVectorView dz_dD0,
                     const ExhaustiveConstVectorView dz_dD0_fac,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] =
        lines[i].dD0(ds_dD0[i], dz_dD0[i], dz_dD0_fac[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dG0(const ExhaustiveConstComplexVectorView& cut,
                        const ExhaustiveConstComplexVectorView dz_dG0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  const auto [s, cs, dz] = frequency_spans(cutoff, f, lines, cut, dz_dG0);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].dG0(dz[i], f) - cs[i];
  }

  return out;
}

void band_shape::dG0(ExhaustiveComplexVectorView cut,
                     const ExhaustiveConstComplexVectorView dz_dG0,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].dG0(dz_dG0[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dY(const ExhaustiveConstComplexVectorView& cut,
                       const ExhaustiveConstComplexVectorView ds_dY,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_dY);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].dY(ds[i], f) - cs[i];
  }

  return out;
}

void band_shape::dY(ExhaustiveComplexVectorView cut,
                    const ExhaustiveConstComplexVectorView ds_dY,
                    const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].dY(ds_dY[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dG(const ExhaustiveConstComplexVectorView& cut,
                       const ExhaustiveConstComplexVectorView ds_dG,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_dG);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].dG(ds[i], f) - cs[i];
  }

  return out;
}

void band_shape::dG(ExhaustiveComplexVectorView cut,
                    const ExhaustiveConstComplexVectorView ds_dG,
                    const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].dG(ds_dG[i], lines[i].f0 + cutoff);
  }
}

void ComputeData::update_zeeman(const Vector2& los,
                                const Vector3& mag,
                                const zeeman::pol pol) {
  npm = zeeman::norm_view(pol, mag, los);
  if (pol != zeeman::pol::no) {
    dnpm_du = zeeman::dnorm_view_du(pol, mag, los);
    dnpm_dv = zeeman::dnorm_view_dv(pol, mag, los);
    dnpm_dw = zeeman::dnorm_view_dw(pol, mag, los);
  }
}

ComputeData::ComputeData(const ExhaustiveConstVectorView& f_grid,
                         const AtmPoint& atm,
                         const Vector2& los,
                         const zeeman::pol pol)
    : scl(f_grid.size()),
      dscl(f_grid.size()),
      shape(f_grid.size()),
      dshape(f_grid.size()) {
  std::transform(f_grid.begin(),
                 f_grid.end(),
                 scl.begin(),
                 [N = number_density(atm.pressure, atm.temperature),
                  T = atm.temperature](auto f) {
                   constexpr Numeric c =
                       Constant::c * Constant::c / (8 * Constant::pi);
                   const Numeric r = (Constant::h * f) / (Constant::k * T);
                   return -N * f * std::expm1(-r) * c;
                 });

  update_zeeman(los, atm.mag, pol);
}

//! Sizes cut, dcut, dz, ds; sets shape
void ComputeData::core_calc(const band_shape& shp,
                            const band_data& bnd,
                            const ExhaustiveConstVectorView& f_grid) {
  cut.resize(shp.size());
  dz.resize(shp.size());
  dz_fac.resize(shp.size());
  ds.resize(shp.size());
  dcut.resize(shp.size());
  filter.reserve(shp.size());

  if (bnd.cutoff != LineByLineCutoffType::None) {
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
void ComputeData::dt_core_calc(const SpeciesIsotope& spec,
                               const band_shape& shp,
                               const band_data& bnd,
                               const ExhaustiveConstVectorView& f_grid,
                               const AtmPoint& atm,
                               const zeeman::pol pol) {
  std::transform(f_grid.begin(),
                 f_grid.end(),
                 dscl.begin(),
                 [N  = number_density(atm.pressure, atm.temperature),
                  dN = dnumber_density_dt(atm.pressure, atm.temperature),
                  T  = atm.temperature](auto f) {
                   constexpr Numeric c =
                       Constant::c * Constant::c / (8 * Constant::pi);
                   const Numeric r = (Constant::h * f) / (Constant::k * T);
                   return -f * (N * r * exp(-r) / T + dN * std::expm1(-r)) * c;
                 });

  const Numeric T = atm.temperature;
  for (Size i = 0; i < pos.size(); i++) {
    const auto& line = bnd.lines[pos[i].line];
    const auto& lshp = shp.lines[i];

    const Numeric& inv_gd = lshp.inv_gd;
    const Numeric& f0     = lshp.f0;

    if (pos[i].spec == std::numeric_limits<Size>::max()) {
      dz_fac[i] =
          (-2 * T * line.ls.dD0_dT(atm) - 2 * T * line.ls.dDV_dT(atm) - f0) /
          (2 * T * f0);

      ds[i] = line.z.Strength(line.qn.val, pol, pos[i].iz) *
              dline_strength_calc_dT(inv_gd, f0, spec, line, atm);

      dz[i] = inv_gd *
              Complex{-dline_center_calc_dT(line, atm), line.ls.dG0_dT(atm)};
    } else {
      const auto& ls = line.ls.single_models[pos[i].spec];

      dz_fac[i] = (-2 * T * ls.dD0_dT(line.ls.T0, T, atm.pressure) -
                   2 * T * ls.dDV_dT(line.ls.T0, T, atm.pressure) - f0) /
                  (2 * T * f0);

      ds[i] = line.z.Strength(line.qn.val, pol, pos[i].iz) *
              dline_strength_calc_dT(f0, inv_gd, spec, line, atm, pos[i].spec);

      dz[i] = inv_gd * Complex{-ls.dD0_dT(line.ls.T0, T, atm.pressure) -
                                   ls.dDV_dT(line.ls.T0, T, atm.pressure),
                               ls.dG0_dT(line.ls.T0, T, atm.pressure)};
    }
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dT(dcut, ds, dz, dz_fac);
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dT(dcut, ds, dz, dz_fac, f);
        });
  } else {
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dT(ds, dz, dz_fac, f);
        });
  }
}

//! Sets dshape and dscl
void ComputeData::df_core_calc(const band_shape& shp,
                               const band_data& bnd,
                               const ExhaustiveConstVectorView& f_grid,
                               const AtmPoint& atm) {
  std::transform(f_grid.begin(),
                 f_grid.end(),
                 dscl.begin(),
                 [N = number_density(atm.pressure, atm.temperature),
                  T = atm.temperature](auto f) {
                   constexpr Numeric c =
                       Constant::c * Constant::c / (8 * Constant::pi);
                   const Numeric r = (Constant::h * f) / (Constant::k * T);
                   return N * (r * std::exp(-r) - std::expm1(-r)) * c;
                 });

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.df(dcut);
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.df(dcut, f);
        });
  } else {
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [&shp](Numeric f) {
          return shp.df(f);
        });
  }
}

//! Sets dshape and dz
void ComputeData::dmag_u_core_calc(const band_shape& shp,
                                   const band_data& bnd,
                                   const ExhaustiveConstVectorView& f_grid,
                                   const AtmPoint& atm,
                                   const zeeman::pol pol) {
  const Numeric H         = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  const Numeric dH_dmag_u = atm.mag[0] / H;

  for (Size i = 0; i < pos.size(); i++) {
    const auto& line = bnd.lines[pos[i].line];
    dz[i]            = -shp.lines[i].inv_gd * dH_dmag_u *
            line.z.Splitting(line.qn.val, pol, pos[i].iz);
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dH(dcut, dz);
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dH(dcut, dz, f);
        });
  } else {
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dH(dz, f);
        });
  }
}

//! Sets dshape and dz
void ComputeData::dmag_v_core_calc(const band_shape& shp,
                                   const band_data& bnd,
                                   const ExhaustiveConstVectorView& f_grid,
                                   const AtmPoint& atm,
                                   const zeeman::pol pol) {
  const Numeric H         = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  const Numeric dH_dmag_v = atm.mag[1] / H;

  for (Size i = 0; i < pos.size(); i++) {
    const auto& line = bnd.lines[pos[i].line];
    dz[i]            = -shp.lines[i].inv_gd * dH_dmag_v *
            line.z.Splitting(line.qn.val, pol, pos[i].iz);
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dH(dcut, dz);
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dH(dcut, dz, f);
        });
  } else {
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dH(dz, f);
        });
  }
}

//! Sets dshape and dz
void ComputeData::dmag_w_core_calc(const band_shape& shp,
                                   const band_data& bnd,
                                   const ExhaustiveConstVectorView& f_grid,
                                   const AtmPoint& atm,
                                   const zeeman::pol pol) {
  const Numeric H         = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  const Numeric dH_dmag_w = atm.mag[2] / H;

  for (Size i = 0; i < pos.size(); i++) {
    const auto& line = bnd.lines[pos[i].line];
    dz[i]            = -shp.lines[i].inv_gd * dH_dmag_w *
            line.z.Splitting(line.qn.val, pol, pos[i].iz);
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dH(dcut, dz);
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dH(dcut, dz, f);
        });
  } else {
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dH(dz, f);
        });
  }
}

//! Sets ds and dz and dcut and dshape
void ComputeData::dVMR_core_calc(const SpeciesIsotope& spec,
                                 const band_shape& shp,
                                 const band_data& bnd,
                                 const ExhaustiveConstVectorView& f_grid,
                                 const AtmPoint& atm,
                                 const zeeman::pol pol,
                                 const SpeciesEnum target_spec) {
  const Numeric x = atm[target_spec];

  for (Size i = 0; i < pos.size(); i++) {
    const auto& line      = bnd.lines[pos[i].line];
    const auto& lshp      = shp.lines[i];
    const Numeric& inv_gd = lshp.inv_gd;
    const Numeric& f0     = lshp.f0;

    if (pos[i].spec == std::numeric_limits<Size>::max()) {
      dz_fac[i] = -(line.ls.dD0_dVMR(atm, target_spec) +
                    line.ls.dDV_dVMR(atm, target_spec)) /
                  f0;

      ds[i] =
          line.z.Strength(line.qn.val, pol, pos[i].iz) *
          dline_strength_calc_dVMR(inv_gd, f0, spec, target_spec, line, atm);

      dz[i] = inv_gd * Complex{-dline_center_calc_dVMR(line, target_spec, atm),
                               line.ls.dG0_dVMR(atm, target_spec)};
    } else {
      const auto ls_spec = line.ls.single_models[pos[i].spec].species;

      dz_fac[i] = 0;

      if (target_spec == ls_spec) {
        ds[i] = lshp.s * (1 + (target_spec == spec.spec)) / x;
      } else if (ls_spec == SpeciesEnum::Bath) {
        const Numeric v = 1.0 - std::transform_reduce(
                                    line.ls.single_models.begin(),
                                    line.ls.single_models.end() - 1,
                                    0.0,
                                    std::plus<>{},
                                    [&atm](auto& s) { return atm[s.species]; });
        ds[i] = lshp.s * (v - x) / (x * v);
      } else {
        ds[i] = 0;
      }

      dz[i] = 0;
    }
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dVMR(dcut, ds, dz, dz_fac);
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dVMR(dcut, ds, dz, dz_fac, f);
        });
  } else {
    std::transform(
        f_grid.begin(), f_grid.end(), dshape.begin(), [this, &shp](Numeric f) {
          return shp.dVMR(ds, dz, dz_fac, f);
        });
  }
}

void ComputeData::set_filter(const line_key& key) {
  if (key.line == filtered_line and key.spec == filtered_spec) return;

  filtered_spec = key.spec;
  filtered_line = key.line;
  filter.resize(0);

  if (not good_enum(key.var)) {
    for (Size i = 0; i < pos.size(); i++) {
      if (pos[i].line == key.line and
          (pos[i].spec == key.spec or
           pos[i].spec == std::numeric_limits<Size>::max()))
        filter.push_back(i);
    }
  } else {
    for (Size i = 0; i < pos.size(); i++) {
      if (pos[i].line == key.line) filter.push_back(i);
    }
  }
}

//! Sets dshape and ds and dz and dcut and dshape
void ComputeData::df0_core_calc(const SpeciesIsotope& spec,
                                const band_shape& shp,
                                const band_data& bnd,
                                const ExhaustiveConstVectorView& f_grid,
                                const AtmPoint& atm,
                                const zeeman::pol pol,
                                const line_key& key) {
  set_filter(key);

  for (Size i : filter) {
    const auto& lshp = shp.lines[i];
    const auto& line = bnd.lines[pos[i].line];

    const Numeric& inv_gd = lshp.inv_gd;
    const Numeric& f0     = lshp.f0;

    if (pos[i].spec == std::numeric_limits<Size>::max()) {
      dz_fac[i] = -1.0 / f0;

      ds[i] = line.z.Strength(line.qn.val, pol, pos[i].iz) *
              dline_strength_calc_df0(f0, inv_gd, spec, line, atm);

      dz[i] = -inv_gd;
    } else {
      dz_fac[i] = -1.0 / f0;

      ds[i] = line.z.Strength(line.qn.val, pol, pos[i].iz) *
              dline_strength_calc_df0(f0, inv_gd, spec, line, atm, pos[i].spec);

      dz[i] = -inv_gd;
    }
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.df0(dcut, ds, dz, dz_fac, filter);
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.df0(dcut, ds, dz, dz_fac, f_grid[i], filter);
    }
  } else {
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.df0(ds, dz, dz_fac, f_grid[i], filter);
    }
  }
}

//! Sets dshape and ds and dcut and dshape
void ComputeData::de0_core_calc(const band_shape& shp,
                                const band_data& bnd,
                                const ExhaustiveConstVectorView& f_grid,
                                const AtmPoint& atm,
                                const line_key& key) {
  using Constant::h, Constant::k;

  set_filter(key);

  for (Size i : filter) {
    const Numeric ds_de0_ratio =
        bnd.lines[pos[i].line].ds_de0_s_ratio(atm.temperature);
    ds[i] = ds_de0_ratio * shp.lines[i].s;
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
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
void ComputeData::da_core_calc(const band_shape& shp,
                               const band_data& bnd,
                               const ExhaustiveConstVectorView& f_grid,
                               const line_key& key) {
  using Constant::h, Constant::k;

  set_filter(key);

  for (Size i : filter) {
    const Numeric ds_da_ratio = 1.0 / bnd.lines[pos[i].line].a;
    ds[i]                     = ds_da_ratio * shp.lines[i].s;
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
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
void ComputeData::dG0_core_calc(const band_shape& shp,
                                const band_data& bnd,
                                const ExhaustiveConstVectorView& f_grid,
                                const AtmPoint& atm,
                                const line_key& key) {
  set_filter(key);

  for (Size i : filter) {
    const auto& ls = bnd.lines[pos[i].line].ls;

    if (pos[i].spec == std::numeric_limits<Size>::max()) {
      dz[i] = Complex(
          0, shp.lines[i].inv_gd * ls.dG0_dX(atm, key.spec, key.ls_coeff));
    } else {
      dz[i] =
          Complex(0,
                  shp.lines[i].inv_gd *
                      ls.single_models[pos[i].spec].dG0_dX(
                          ls.T0, atm.temperature, atm.pressure, key.ls_coeff));
    }
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
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
void ComputeData::dD0_core_calc(const band_shape& shp,
                                const band_data& bnd,
                                const ExhaustiveConstVectorView& f_grid,
                                const AtmPoint& atm,
                                const line_key& key) {
  set_filter(key);

  for (Size i : filter) {
    const auto& lshp = shp.lines[i];
    const auto& ls   = bnd.lines[pos[i].line].ls;

    const Numeric& inv_gd = lshp.inv_gd;
    const Numeric& f0     = lshp.f0;

    if (pos[i].spec == std::numeric_limits<Size>::max()) {
      const Numeric d = ls.dD0_dX(atm, key.spec, key.ls_coeff);

      dz_fac[i] = -d / f0;

      ds[i] = -d * lshp.s / f0;

      dz[i] = -d * inv_gd;
    } else {
      const Numeric d = ls.single_models[pos[i].spec].dD0_dX(
          ls.T0, atm.temperature, atm.pressure, key.ls_coeff);

      dz_fac[i] = -d / f0;

      ds[i] = -d * lshp.s / f0;

      dz[i] = -d * inv_gd;
    }
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dD0(dcut, ds, dz, dz_fac, filter);
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dD0(dcut, ds, dz, dz_fac, f_grid[i], filter);
    }
  } else {
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dD0(ds, dz, dz_fac, f_grid[i], filter);
    }
  }
}

//! Sets dshape and ds and dcut and dshape
void ComputeData::dY_core_calc(const SpeciesIsotope& spec,
                               const band_shape& shp,
                               const band_data& bnd,
                               const ExhaustiveConstVectorView& f_grid,
                               const AtmPoint& atm,
                               const zeeman::pol pol,
                               const line_key& key) {
  set_filter(key);

  for (Size i : filter) {
    const auto& line = bnd.lines[pos[i].line];
    const auto& lshp = shp.lines[i];

    if (pos[i].spec == std::numeric_limits<Size>::max()) {
      ds[i] = line.z.Strength(line.qn.val, pol, pos[i].iz) *
              dline_strength_calc_dY(line.ls.dY_dX(atm, key.spec, key.ls_coeff),
                                     lshp.inv_gd,
                                     spec,
                                     line,
                                     atm);
    } else {
      ds[i] = line.z.Strength(line.qn.val, pol, pos[i].iz) *
              dline_strength_calc_dY(
                  line.ls.single_models[pos[i].spec].dY_dX(
                      line.ls.T0, atm.temperature, atm.pressure, key.ls_coeff),
                  lshp.inv_gd,
                  spec,
                  line,
                  atm,
                  pos[i].spec);
    }
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
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
void ComputeData::dG_core_calc(const SpeciesIsotope& spec,
                               const band_shape& shp,
                               const band_data& bnd,
                               const ExhaustiveConstVectorView& f_grid,
                               const AtmPoint& atm,
                               const zeeman::pol pol,
                               const line_key& key) {
  set_filter(key);

  for (Size i : filter) {
    const auto& line = bnd.lines[pos[i].line];
    const auto& lshp = shp.lines[i];

    if (pos[i].spec == std::numeric_limits<Size>::max()) {
      ds[i] = line.z.Strength(line.qn.val, pol, pos[i].iz) *
              dline_strength_calc_dG(line.ls.dG_dX(atm, key.spec, key.ls_coeff),
                                     lshp.inv_gd,
                                     spec,
                                     line,
                                     atm);
    } else {
      ds[i] = line.z.Strength(line.qn.val, pol, pos[i].iz) *
              dline_strength_calc_dG(
                  line.ls.single_models[pos[i].spec].dG_dX(
                      line.ls.T0, atm.temperature, atm.pressure, key.ls_coeff),
                  lshp.inv_gd,
                  spec,
                  line,
                  atm,
                  pos[i].spec);
    }
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
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
void ComputeData::dDV_core_calc(const band_shape& shp,
                                const band_data& bnd,
                                const ExhaustiveConstVectorView& f_grid,
                                const AtmPoint& atm,
                                const line_key& key) {
  using Constant::h, Constant::k;

  set_filter(key);

  for (Size i : filter) {
    const auto& lshp = shp.lines[i];
    const auto& ls   = bnd.lines[pos[i].line].ls;

    const Numeric& inv_gd = lshp.inv_gd;
    const Numeric& f0     = lshp.f0;

    if (pos[i].spec == std::numeric_limits<Size>::max()) {
      const Numeric d = ls.dDV_dX(atm, key.spec, key.ls_coeff);

      dz_fac[i] = -d / f0;

      ds[i] = lshp.s * dz_fac[i];

      dz[i] = -d * inv_gd;
    } else {
      const Numeric d = ls.single_models[pos[i].spec].dDV_dX(
          ls.T0, atm.temperature, atm.pressure, key.ls_coeff);

      dz_fac[i] = -d / f0;

      ds[i] = lshp.s * dz_fac[i];

      dz[i] = -d * inv_gd;
    }
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dDV(dcut, ds, dz, dz_fac, filter);
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dDV(dcut, ds, dz, dz_fac, f_grid[i], filter);
    }
  } else {
    for (Index i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dDV(ds, dz, dz_fac, f_grid[i], filter);
    }
  }
}

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const ExhaustiveConstVectorView& f_grid,
                        const SpeciesIsotope& spec,
                        const band_shape& shape,
                        const band_data& bnd,
                        const AtmPoint& atm,
                        const zeeman::pol pol,
                        const AtmKey& key) {
  using enum AtmKey;
  switch (key) {
    case t:
      com_data.dt_core_calc(spec, shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dscl[i] * com_data.shape[i] +
                                    com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case p: ARTS_USER_ERROR("Not implemented, pressure derivative"); break;
    case mag_u:
      if (pol == zeeman::pol::no) return;
      com_data.dmag_u_core_calc(shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dnpm_du,
                                com_data.scl[i] * com_data.shape[i],
                                com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case mag_v:
      if (pol == zeeman::pol::no) return;
      com_data.dmag_v_core_calc(shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dnpm_dv,
                                com_data.scl[i] * com_data.shape[i],
                                com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case mag_w:
      if (pol == zeeman::pol::no) return;
      com_data.dmag_w_core_calc(shape, bnd, f_grid, atm, pol);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dnpm_dw,
                                com_data.scl[i] * com_data.shape[i],
                                com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case wind_u:
    case wind_v:
    case wind_w:
      com_data.df_core_calc(shape, bnd, f_grid, atm);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dscl[i] * com_data.shape[i] +
                                    com_data.scl[i] * com_data.dshape[i]);
      }
      return;
  }
}

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const ExhaustiveConstVectorView& f_grid,
                        const SpeciesIsotope& spec,
                        const band_shape&,
                        const band_data&,
                        const AtmPoint& atm,
                        const zeeman::pol,
                        const SpeciesIsotope& deriv_spec) {
  if (deriv_spec != spec) return;

  const Numeric isorat = atm[spec];

  ARTS_USER_ERROR_IF(
      isorat == 0,
      "Does not support 0 for isotopologue ratios (may be added upon request)")

  for (Index i = 0; i < f_grid.size(); i++) {
    const auto dF  = com_data.scl[i] * com_data.shape[i] / isorat;
    dpm[i]        += zeeman::scale(com_data.npm, dF);
  }
}

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const ExhaustiveConstVectorView& f_grid,
                        const SpeciesIsotope& spec,
                        const band_shape& shape,
                        const band_data& bnd,
                        const AtmPoint& atm,
                        const zeeman::pol pol,
                        const SpeciesEnum& deriv_spec) {
  com_data.dVMR_core_calc(spec, shape, bnd, f_grid, atm, pol, deriv_spec);
  for (Index i = 0; i < f_grid.size(); i++) {
    dpm[i] += zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
  }
}

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const ExhaustiveConstVectorView& f_grid,
                        const SpeciesIsotope& spec,
                        const band_shape& shape,
                        const band_data& bnd,
                        const AtmPoint& atm,
                        const zeeman::pol pol,
                        const line_key& deriv) {
  switch (deriv.var) {
    case LineByLineVariable::f0:
      com_data.df0_core_calc(spec, shape, bnd, f_grid, atm, pol, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case LineByLineVariable::e0:
      com_data.de0_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case LineByLineVariable::a:
      com_data.da_core_calc(shape, bnd, f_grid, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
  }

  switch (deriv.ls_var) {
    case LineShapeModelVariable::G0:
      com_data.dG0_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case LineShapeModelVariable::D0:
      com_data.dD0_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case LineShapeModelVariable::G2:  return;
    case LineShapeModelVariable::D2:  return;
    case LineShapeModelVariable::FVC: return;
    case LineShapeModelVariable::ETA: return;
    case LineShapeModelVariable::Y:
      com_data.dY_core_calc(spec, shape, bnd, f_grid, atm, pol, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case LineShapeModelVariable::G:
      com_data.dG_core_calc(spec, shape, bnd, f_grid, atm, pol, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case LineShapeModelVariable::DV:
      com_data.dDV_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Index i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
  }
}

void compute_derivative(PropmatVectorView,
                        ComputeData&,
                        const ExhaustiveConstVectorView&,
                        const SpeciesIsotope&,
                        const band_shape&,
                        const band_data&,
                        const AtmPoint&,
                        const zeeman::pol,
                        const auto&) {}

std::ostream& operator<<(std::ostream& os, const ComputeData& cd) {
  for (auto line : cd.lines) {
    os << line.f0 << ' ' << line.s << ' ' << line.inv_gd << ' ' << line.z_imag
       << '\n';
  }
  return os;
}

void calculate(PropmatVectorView pm,
               matpack::matpack_view<Propmat, 2, false, true> dpm,
               ComputeData& com_data,
               const ExhaustiveConstVectorView& f_grid,
               const JacobianTargets& jacobian_targets,
               const QuantumIdentifier& bnd_qid,
               const band_data& bnd,
               const AtmPoint& atm,
               const zeeman::pol pol,
               const bool no_negative_absorption) {
  if (std::ranges::all_of(com_data.npm, [](auto& n) { return n == 0; })) return;

  const Index nf = f_grid.size();
  if (nf == 0) return;

  const SpeciesIsotope spec = bnd_qid.Isotopologue();
  const Numeric fmin        = f_grid.front();
  const Numeric fmax        = f_grid.back();

  ARTS_ASSERT(jacobian_targets.target_count() ==
                  static_cast<Size>(dpm.nrows()) and
              nf == dpm.ncols())
  ARTS_ASSERT(nf == pm.size())

  band_shape_helper(
      com_data.lines, com_data.pos, spec, bnd, atm, fmin, fmax, pol);
  if (com_data.lines.empty()) return;

  //! Not const to save lines for reuse
  band_shape shape{std::move(com_data.lines), bnd.get_cutoff_frequency()};

  com_data.core_calc(shape, bnd, f_grid);

  for (Index i = 0; i < nf; i++) {
    const auto F = com_data.scl[i] * com_data.shape[i];
    if (no_negative_absorption and F.real() < 0) continue;
    pm[i] += zeeman::scale(com_data.npm, F);
  }

  for (auto& atm_target : jacobian_targets.atm()) {
    std::visit(
        [&](auto& target) {
          compute_derivative(dpm.as_slice(atm_target.target_pos),
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
      compute_derivative(dpm.as_slice(line_target.target_pos),
                         com_data,
                         f_grid,
                         spec,
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
