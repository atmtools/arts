#include "lbl_lineshape_voigt_lte.h"

#include <arts_omp.h>
#include <atm.h>
#include <enumsSpeciesEnum.h>
#include <jacobian.h>
#include <partfun.h>
#include <physics_funcs.h>
#include <quantum.h>
#include <sorting.h>

#include <Faddeeva/Faddeeva.hh>
#include <algorithm>
#include <cmath>
#include <numeric>

#include "lbl_data.h"
#include "lbl_zeeman.h"

namespace lbl::voigt::lte {
namespace {
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
                           const SpeciesEnum ispec) {
  const auto& ls   = line.ls.single_models.at(ispec);
  const Numeric T0 = line.ls.T0;
  const Numeric T  = atm.temperature;
  const Numeric P  = atm.pressure;
  const Numeric x  = atm[spec.spec];
  const Numeric r  = atm[spec];
  const Numeric v =
      ispec == SpeciesEnum::Bath
          ? 1 - std::transform_reduce(
                    line.ls.single_models.begin(),
                    line.ls.single_models.end(),
                    0.0,
                    std::plus<>{},
                    [&atm](auto& s) {
                      return s.first == SpeciesEnum::Bath ? 0.0 : atm[s.first];
                    })
          : atm[ispec];

  const auto s = line.s(T, PartitionFunctions::Q(T, spec));

  const auto G = ls.G(T0, T, P);
  const auto Y = ls.Y(T0, T, P);
  const Complex lm{1 + G, -Y};

  return Constant::inv_sqrt_pi * inv_gd * x * r * v * lm * s;
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

Numeric line_center_calc(const line& line,
                         const AtmPoint& atm,
                         SpeciesEnum ispec) {
  const auto& ls = line.ls.single_models.at(ispec);
  return line.f0 + ls.D0(line.ls.T0, atm.temperature, atm.pressure) +
         ls.DV(line.ls.T0, atm.temperature, atm.pressure);
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

  single_shape_builder(const SpeciesIsotope& s,
                       const line& l,
                       const AtmPoint& a)
      : spec(s),
        ln(l),
        atm(a),
        f0(line_center_calc(l, atm)),
        scaled_gd_part(std::sqrt(Constant::doppler_broadening_const_squared *
                                 atm.temperature / s.mass)),
        G0(l.ls.G0(atm)) {}


  [[nodiscard]] single_shape as_zeeman(const Numeric H,
                                       const ZeemanPolarization pol,
                                       const Size iz) const {
    single_shape s;
    s.f0     = f0 + H * ln.z.Splitting(ln.qn, pol, iz);
    s.inv_gd = 1.0 / (scaled_gd_part * f0);
    s.z_imag = G0 * s.inv_gd;
    s.s      = ln.z.Strength(ln.qn, pol, iz) *
          line_strength_calc(s.inv_gd, spec, ln, atm);
    return s;
  }

  operator single_shape() const {
    single_shape s;
    s.f0     = f0;
    s.inv_gd = 1.0 / (scaled_gd_part * f0);
    s.z_imag = G0 * s.inv_gd;
    s.s      = line_strength_calc(s.inv_gd, spec, ln, atm);
    return s;
  }
};
}  // namespace

[[nodiscard]] single_shape single_shape::update_iz(const SpeciesIsotope& spec,
                                                   const line& line,
                                                   const AtmPoint& atm,
                                                   const ZeemanPolarization pol,
                                                   const Index iz) const {
  single_shape out{*this};

  if (iz != 0) {
    out.f0 = line_center_calc(line, atm) +
             std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]) *
                 line.z.Splitting(line.qn, pol, iz);
    out.inv_gd = 1.0 / scaled_gd(atm.temperature, spec.mass, out.f0);
    out.s      = line.z.Strength(line.qn, pol, iz) *
            line_strength_calc(out.inv_gd, spec, line, atm);
  }

  return out;
}

single_shape::single_shape(const SpeciesIsotope& spec,
                           const line& line,
                           const AtmPoint& atm,
                           const ZeemanPolarization pol,
                           const Index iz)
    : f0(line_center_calc(line, atm) +
         std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]) *
             line.z.Splitting(line.qn, pol, iz)),
      inv_gd(1.0 / scaled_gd(atm.temperature, spec.mass, f0)),
      z_imag(line.ls.G0(atm) * inv_gd),
      s(line.z.Strength(line.qn, pol, iz) *
        line_strength_calc(inv_gd, spec, line, atm)) {}

single_shape::single_shape(const SpeciesIsotope& spec,
                           const line& line,
                           const AtmPoint& atm,
                           const ZeemanPolarization pol,
                           const Index iz,
                           const SpeciesEnum ispec)
    : f0(line_center_calc(line, atm, ispec) +
         std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]) *
             line.z.Splitting(line.qn, pol, iz)),
      inv_gd(1.0 / scaled_gd(atm.temperature, spec.mass, f0)),
      z_imag(line.ls.single_models.at(ispec).G0(
                 line.ls.T0, atm.temperature, atm.pressure) *
             inv_gd),
      s(line.z.Strength(line.qn, pol, iz) *
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

Size count_lines(const band_data& bnd, const ZeemanPolarization type) {
  return std::transform_reduce(
      bnd.begin(), bnd.end(), Index{}, std::plus<>{}, [type](auto& line) {
        const Index factor = 1;
        return factor * line.z.size(line.qn, type);
      });
}

namespace {
void zeeman_push_back(std::vector<single_shape>& lines,
                      std::vector<line_pos>& pos,
                      const single_shape_builder& s,
                      const line& line,
                      const AtmPoint& atm,
                      const ZeemanPolarization pol,
                      const Size iline) {
  if (pol == ZeemanPolarization::no) {
    lines.emplace_back(s);
    pos.emplace_back(line_pos{.line = iline});
  } else {
    const Numeric H = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
    const auto nz   = static_cast<Size>(line.z.size(line.qn, pol));
    for (Size iz = 0; iz < nz; iz++) {
      lines.emplace_back(s.as_zeeman(H, pol, iz));
      pos.emplace_back(line_pos{.line = iline, .iz = iz});

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
                     const ZeemanPolarization pol,
                     const Size iline) {
  if ((line.z.on and pol != ZeemanPolarization::no) or
      (not line.z.on and pol == ZeemanPolarization::no)) {
    zeeman_push_back(lines,
                     pos,
                     single_shape_builder{spec, line, atm},
                     line,
                     atm,
                     pol,
                     iline);
  }
}
}  // namespace

void band_shape_helper(std::vector<single_shape>& lines,
                       std::vector<line_pos>& pos,
                       const SpeciesIsotope& spec,
                       const band_data& bnd,
                       const AtmPoint& atm,
                       const Numeric fmin,
                       const Numeric fmax,
                       const ZeemanPolarization pol) {
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

Complex band_shape::dH(const ConstComplexVectorView& dz_dH,
                       const Numeric f) const {
  assert(static_cast<Size>(dz_dH.size()) == lines.size());

  return std::transform_reduce(lines.begin(),
                               lines.end(),
                               dz_dH.begin(),
                               Complex{},
                               std::plus<>{},
                               [f](auto& ls, auto& d) { return ls.dH(d, f); });
}

Complex band_shape::dT(const ConstComplexVectorView& ds_dT,
                       const ConstComplexVectorView& dz_dT,
                       const ConstVectorView& dz_dT_fac,
                       const Numeric f) const {
  assert(ds_dT.size() == dz_dT.size());
  assert(static_cast<Size>(ds_dT.size()) == lines.size());

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i = 0; i < lines.size(); ++i) {
    out += lines[i].dT(ds_dT[i], dz_dT[i], dz_dT_fac[i], f);
  }

  return out;
}

Complex band_shape::dVMR(const ConstComplexVectorView& ds_dVMR,
                         const ConstComplexVectorView& dz_dVMR,
                         const ConstVectorView& dz_dVMR_fac,
                         const Numeric f) const {
  assert(ds_dVMR.size() == dz_dVMR.size());
  assert(static_cast<Size>(ds_dVMR.size()) == lines.size());

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i = 0; i < lines.size(); ++i) {
    out += lines[i].dVMR(ds_dVMR[i], dz_dVMR[i], dz_dVMR_fac[i], f);
  }

  return out;
}

Complex band_shape::df0(const ConstComplexVectorView ds_df0,
                        const ConstComplexVectorView dz_df0,
                        const ConstVectorView dz_df0_fac,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].df0(ds_df0[i], dz_df0[i], dz_df0_fac[i], f);
  }

  return out;
}

Complex band_shape::da(const ConstComplexVectorView ds_da,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].da(ds_da[i], f);
  }

  return out;
}

Complex band_shape::de0(const ConstComplexVectorView ds_de0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].de0(ds_de0[i], f);
  }

  return out;
}

Complex band_shape::dDV(const ConstComplexVectorView ds_dDV,
                        const ConstComplexVectorView dz_dDV,
                        const ConstVectorView dz_dDV_fac,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dDV(ds_dDV[i], dz_dDV[i], dz_dDV_fac[i], f);
  }

  return out;
}

Complex band_shape::dD0(const ConstComplexVectorView ds_dD0,
                        const ConstComplexVectorView dz_dD0,
                        const ConstVectorView dz_dD0_fac,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dD0(ds_dD0[i], dz_dD0[i], dz_dD0_fac[i], f);
  }

  return out;
}

Complex band_shape::dG0(const ConstComplexVectorView dz_dG0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dG0(dz_dG0[i], f);
  }

  return out;
}

Complex band_shape::dY(const ConstComplexVectorView ds_dY,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dY(ds_dY[i], f);
  }

  return out;
}

Complex band_shape::dG(const ConstComplexVectorView ds_dG,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += lines[i].dG(ds_dG[i], f);
  }

  return out;
}

Complex band_shape::operator()(const ConstComplexVectorView& cut,
                               const Numeric f) const {
  const auto [s, cs] = frequency_spans(cutoff, f, lines, cut);
  return std::transform_reduce(s.begin(),
                               s.end(),
                               cs.begin(),
                               Complex{},
                               std::plus<>{},
                               [f](auto& ls, auto& c) { return ls(f) - c; });
}

void band_shape::operator()(ComplexVectorView cut) const {
  std::transform(
      lines.begin(),
      lines.end(),
      cut.begin(),
      [cutoff_freq = cutoff](auto& ls) { return ls(ls.f0 + cutoff_freq); });
}

Complex band_shape::df(const ConstComplexVectorView& cut,
                       const Numeric f) const {
  const auto [s, cs] = frequency_spans(cutoff, f, lines, cut);
  return std::transform_reduce(s.begin(),
                               s.end(),
                               cs.begin(),
                               Complex{},
                               std::plus<>{},
                               [f](auto& ls, auto& c) { return ls.df(f) - c; });
}

void band_shape::df(ComplexVectorView cut) const {
  std::transform(
      lines.begin(),
      lines.end(),
      cut.begin(),
      [cutoff_freq = cutoff](auto& ls) { return ls.df(ls.f0 + cutoff_freq); });
}

Complex band_shape::dH(const ConstComplexVectorView& cut,
                       const ConstComplexVectorView& dz_dH,
                       const Numeric f) const {
  assert(static_cast<Size>(dz_dH.size()) == lines.size());

  const auto [s, cs, dH] = frequency_spans(cutoff, f, lines, cut, dz_dH);

  Complex out{};  //! Fixme, use zip in C++ 23...
  for (Size i = 0; i < s.size(); ++i) {
    out += s[i].dH(dH[i], f) - cs[i];
  }

  return out;
}

void band_shape::dH(ComplexVectorView cut,
                    const ConstComplexVectorView& df0_dH) const {
  assert(static_cast<Size>(df0_dH.size()) == lines.size());

  std::transform(lines.begin(),
                 lines.end(),
                 df0_dH.begin(),
                 cut.begin(),
                 [cutoff_freq = cutoff](auto& ls, auto& d) {
                   return ls.dH(d, ls.f0 + cutoff_freq);
                 });
}

Complex band_shape::dT(const ConstComplexVectorView& cut,
                       const ConstComplexVectorView& ds_dT,
                       const ConstComplexVectorView& dz_dT,
                       const ConstVectorView& dz_dT_fac,
                       const Numeric f) const {
  assert(ds_dT.size() == dz_dT.size());
  assert(static_cast<Size>(ds_dT.size()) == lines.size());

  Complex out{};  //! Fixme, use zip in C++ 23...

  const auto [s, cs, ds, dz, dzf] =
      frequency_spans(cutoff, f, lines, cut, ds_dT, dz_dT, dz_dT_fac);

  for (Size i = 0; i < s.size(); ++i) {
    out += s[i].dT(ds[i], dz[i], dzf[i], f) - cs[i];
  }

  return out;
}

void band_shape::dT(ComplexVectorView cut,
                    const ConstComplexVectorView& ds_dT,
                    const ConstComplexVectorView& dz_dT,
                    const ConstVectorView& dz_dT_fac) const {
  assert(ds_dT.size() == dz_dT.size());
  assert(static_cast<Size>(ds_dT.size()) == lines.size());

  for (Size i = 0; i < lines.size(); ++i) {
    cut[i] =
        lines[i].dT(ds_dT[i], dz_dT[i], dz_dT_fac[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dVMR(const ConstComplexVectorView& cut,
                         const ConstComplexVectorView& ds_dVMR,
                         const ConstComplexVectorView& dz_dVMR,
                         const ConstVectorView& dz_dVMR_fac,
                         const Numeric f) const {
  assert(ds_dVMR.size() == dz_dVMR.size());
  assert(static_cast<Size>(ds_dVMR.size()) == lines.size());

  Complex out{};  //! Fixme, use zip in C++ 23...

  const auto [s, cs, ds, dz, dzf] =
      frequency_spans(cutoff, f, lines, cut, ds_dVMR, dz_dVMR, dz_dVMR_fac);

  for (Size i = 0; i < s.size(); ++i) {
    out += s[i].dVMR(ds[i], dz[i], dzf[i], f) - cs[i];
  }

  return out;
}

void band_shape::dVMR(ComplexVectorView cut,
                      const ConstComplexVectorView& ds_dVMR,
                      const ConstComplexVectorView& dz_dVMR,
                      const ConstVectorView& dz_dVMR_fac) const {
  assert(ds_dVMR.size() == dz_dVMR.size());
  assert(static_cast<Size>(ds_dVMR.size()) == lines.size());

  for (Size i = 0; i < lines.size(); ++i) {
    cut[i] = lines[i].dVMR(
        ds_dVMR[i], dz_dVMR[i], dz_dVMR_fac[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::df0(const ConstComplexVectorView& cut,
                        const ConstComplexVectorView ds_df0,
                        const ConstComplexVectorView dz_df0,
                        const ConstVectorView dz_df0_fac,
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

void band_shape::df0(ComplexVectorView cut,
                     const ConstComplexVectorView ds_df0,
                     const ConstComplexVectorView dz_df0,
                     const ConstVectorView dz_df0_fac,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] =
        lines[i].df0(ds_df0[i], dz_df0[i], dz_df0_fac[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::da(const ConstComplexVectorView& cut,
                       const ConstComplexVectorView ds_da,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_da);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].da(ds[i], f) - cs[i];
  }

  return out;
}

void band_shape::da(ComplexVectorView cut,
                    const ConstComplexVectorView ds_da,
                    const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].da(ds_da[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::de0(const ConstComplexVectorView& cut,
                        const ConstComplexVectorView ds_de0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_de0);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].de0(ds[i], f) - cs[i];
  }

  return out;
}

void band_shape::de0(ComplexVectorView cut,
                     const ConstComplexVectorView ds_de0,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].de0(ds_de0[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dDV(const ConstComplexVectorView& cut,
                        const ConstComplexVectorView ds_dDV,
                        const ConstComplexVectorView dz_dDV,
                        const ConstVectorView dz_dDV_fac,
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

void band_shape::dDV(ComplexVectorView cut,
                     const ConstComplexVectorView ds_dDV,
                     const ConstComplexVectorView dz_dDV,
                     const ConstVectorView dz_dDV_fac,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] =
        lines[i].dDV(ds_dDV[i], dz_dDV[i], dz_dDV_fac[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dD0(const ConstComplexVectorView& cut,
                        const ConstComplexVectorView ds_dD0,
                        const ConstComplexVectorView dz_dD0,
                        const ConstVectorView dz_dD0_fac,
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

void band_shape::dD0(ComplexVectorView cut,
                     const ConstComplexVectorView ds_dD0,
                     const ConstComplexVectorView dz_dD0,
                     const ConstVectorView dz_dD0_fac,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] =
        lines[i].dD0(ds_dD0[i], dz_dD0[i], dz_dD0_fac[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dG0(const ConstComplexVectorView& cut,
                        const ConstComplexVectorView dz_dG0,
                        const Numeric f,
                        const std::vector<Size>& filter) const {
  const auto [s, cs, dz] = frequency_spans(cutoff, f, lines, cut, dz_dG0);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].dG0(dz[i], f) - cs[i];
  }

  return out;
}

void band_shape::dG0(ComplexVectorView cut,
                     const ConstComplexVectorView dz_dG0,
                     const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].dG0(dz_dG0[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dY(const ConstComplexVectorView& cut,
                       const ConstComplexVectorView ds_dY,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_dY);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].dY(ds[i], f) - cs[i];
  }

  return out;
}

void band_shape::dY(ComplexVectorView cut,
                    const ConstComplexVectorView ds_dY,
                    const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].dY(ds_dY[i], lines[i].f0 + cutoff);
  }
}

Complex band_shape::dG(const ConstComplexVectorView& cut,
                       const ConstComplexVectorView ds_dG,
                       const Numeric f,
                       const std::vector<Size>& filter) const {
  const auto [s, cs, ds] = frequency_spans(cutoff, f, lines, cut, ds_dG);

  Complex out{};  //! Fixme, use zip in C++ 23...

  for (Size i : filter) {
    out += s[i].dG(ds[i], f) - cs[i];
  }

  return out;
}

void band_shape::dG(ComplexVectorView cut,
                    const ConstComplexVectorView ds_dG,
                    const std::vector<Size>& filter) const {
  for (Size i : filter) {
    cut[i] = lines[i].dG(ds_dG[i], lines[i].f0 + cutoff);
  }
}

void ComputeData::update_zeeman(const Vector2& los,
                                const Vector3& mag,
                                const ZeemanPolarization pol) {
  npm = zeeman::norm_view(pol, mag, los);
  if (pol != ZeemanPolarization::no) {
    dnpm_du = zeeman::dnorm_view_du(pol, mag, los);
    dnpm_dv = zeeman::dnorm_view_dv(pol, mag, los);
    dnpm_dw = zeeman::dnorm_view_dw(pol, mag, los);
  }
}

ComputeData::ComputeData(const ConstVectorView& f_grid,
                         const AtmPoint& atm,
                         const Vector2& los,
                         const ZeemanPolarization pol)
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
                            const ConstVectorView& f_grid) {
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
                               const ConstVectorView& f_grid,
                               const AtmPoint& atm,
                               const ZeemanPolarization pol) {
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

    dz_fac[i] =
        (-2 * T * line.ls.dD0_dT(atm) - 2 * T * line.ls.dDV_dT(atm) - f0) /
        (2 * T * f0);

    ds[i] = line.z.Strength(line.qn, pol, pos[i].iz) *
            dline_strength_calc_dT(inv_gd, f0, spec, line, atm);

    dz[i] = inv_gd *
            Complex{-dline_center_calc_dT(line, atm), line.ls.dG0_dT(atm)};
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
                               const ConstVectorView& f_grid,
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
                                   const ConstVectorView& f_grid,
                                   const AtmPoint& atm,
                                   const ZeemanPolarization pol) {
  const Numeric H         = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  const Numeric dH_dmag_u = atm.mag[0] / H;

  for (Size i = 0; i < pos.size(); i++) {
    const auto& line = bnd.lines[pos[i].line];
    dz[i]            = -shp.lines[i].inv_gd * dH_dmag_u *
            line.z.Splitting(line.qn, pol, pos[i].iz);
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
                                   const ConstVectorView& f_grid,
                                   const AtmPoint& atm,
                                   const ZeemanPolarization pol) {
  const Numeric H         = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  const Numeric dH_dmag_v = atm.mag[1] / H;

  for (Size i = 0; i < pos.size(); i++) {
    const auto& line = bnd.lines[pos[i].line];
    dz[i]            = -shp.lines[i].inv_gd * dH_dmag_v *
            line.z.Splitting(line.qn, pol, pos[i].iz);
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
                                   const ConstVectorView& f_grid,
                                   const AtmPoint& atm,
                                   const ZeemanPolarization pol) {
  const Numeric H         = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
  const Numeric dH_dmag_w = atm.mag[2] / H;

  for (Size i = 0; i < pos.size(); i++) {
    const auto& line = bnd.lines[pos[i].line];
    dz[i]            = -shp.lines[i].inv_gd * dH_dmag_w *
            line.z.Splitting(line.qn, pol, pos[i].iz);
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
                                 const ConstVectorView& f_grid,
                                 const AtmPoint& atm,
                                 const ZeemanPolarization pol,
                                 const SpeciesEnum target_spec) {
  for (Size i = 0; i < pos.size(); i++) {
    const auto& line      = bnd.lines[pos[i].line];
    const auto& lshp      = shp.lines[i];
    const Numeric& inv_gd = lshp.inv_gd;
    const Numeric& f0     = lshp.f0;

    dz_fac[i] = -(line.ls.dD0_dVMR(atm, target_spec) +
                  line.ls.dDV_dVMR(atm, target_spec)) /
                f0;

    ds[i] =
        line.z.Strength(line.qn, pol, pos[i].iz) *
        dline_strength_calc_dVMR(inv_gd, f0, spec, target_spec, line, atm);

    dz[i] = inv_gd * Complex{-dline_center_calc_dVMR(line, target_spec, atm),
                             line.ls.dG0_dVMR(atm, target_spec)};
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

  for (Size i = 0; i < pos.size(); i++) {
    if (pos[i].line == key.line) filter.push_back(i);
  }
}

//! Sets dshape and ds and dz and dcut and dshape
void ComputeData::df0_core_calc(const SpeciesIsotope& spec,
                                const band_shape& shp,
                                const band_data& bnd,
                                const ConstVectorView& f_grid,
                                const AtmPoint& atm,
                                const ZeemanPolarization pol,
                                const line_key& key) {
  set_filter(key);

  for (Size i : filter) {
    const auto& lshp = shp.lines[i];
    const auto& line = bnd.lines[pos[i].line];

    const Numeric& inv_gd = lshp.inv_gd;
    const Numeric& f0     = lshp.f0;

    dz_fac[i] = -1.0 / f0;

    ds[i] = line.z.Strength(line.qn, pol, pos[i].iz) *
            dline_strength_calc_df0(f0, inv_gd, spec, line, atm);

    dz[i] = -inv_gd;
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.df0(dcut, ds, dz, dz_fac, filter);
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.df0(dcut, ds, dz, dz_fac, f_grid[i], filter);
    }
  } else {
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.df0(ds, dz, dz_fac, f_grid[i], filter);
    }
  }
}

//! Sets dshape and ds and dcut and dshape
void ComputeData::de0_core_calc(const band_shape& shp,
                                const band_data& bnd,
                                const ConstVectorView& f_grid,
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
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.de0(dcut, ds, f_grid[i], filter);
    }
  } else {
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.de0(ds, f_grid[i], filter);
    }
  }
}

//! Sets dshape and ds and dcut and dshape
void ComputeData::da_core_calc(const band_shape& shp,
                               const band_data& bnd,
                               const ConstVectorView& f_grid,
                               const line_key& key) {
  using Constant::h, Constant::k;

  set_filter(key);

  for (Size i : filter) {
    const Numeric ds_da_ratio = 1.0 / bnd.lines[pos[i].line].a;
    ds[i]                     = ds_da_ratio * shp.lines[i].s;
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.da(dcut, ds, filter);
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.da(dcut, ds, f_grid[i], filter);
    }
  } else {
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.da(ds, f_grid[i], filter);
    }
  }
}

//! Sets dshape and dz and dcut and dshape
void ComputeData::dG0_core_calc(const band_shape& shp,
                                const band_data& bnd,
                                const ConstVectorView& f_grid,
                                const AtmPoint& atm,
                                const line_key& key) {
  set_filter(key);

  for (Size i : filter) {
    const auto& ls = bnd.lines[pos[i].line].ls;

    dz[i] = Complex(
        0, shp.lines[i].inv_gd * ls.dG0_dX(atm, key.spec, key.ls_coeff));
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dG0(dcut, dz, filter);
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dG0(dcut, dz, f_grid[i], filter);
    }
  } else {
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dG0(dz, f_grid[i], filter);
    }
  }
}

//! Sets dshape and dz and dcut and dshape
void ComputeData::dD0_core_calc(const band_shape& shp,
                                const band_data& bnd,
                                const ConstVectorView& f_grid,
                                const AtmPoint& atm,
                                const line_key& key) {
  set_filter(key);

  for (Size i : filter) {
    const auto& lshp = shp.lines[i];
    const auto& ls   = bnd.lines[pos[i].line].ls;

    const Numeric& inv_gd = lshp.inv_gd;
    const Numeric& f0     = lshp.f0;

    const Numeric d = ls.dD0_dX(atm, key.spec, key.ls_coeff);

    dz_fac[i] = -d / f0;

    ds[i] = -d * lshp.s / f0;

    dz[i] = -d * inv_gd;
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dD0(dcut, ds, dz, dz_fac, filter);
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dD0(dcut, ds, dz, dz_fac, f_grid[i], filter);
    }
  } else {
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dD0(ds, dz, dz_fac, f_grid[i], filter);
    }
  }
}

//! Sets dshape and ds and dcut and dshape
void ComputeData::dY_core_calc(const SpeciesIsotope& spec,
                               const band_shape& shp,
                               const band_data& bnd,
                               const ConstVectorView& f_grid,
                               const AtmPoint& atm,
                               const ZeemanPolarization pol,
                               const line_key& key) {
  set_filter(key);

  for (Size i : filter) {
    const auto& line = bnd.lines[pos[i].line];
    const auto& lshp = shp.lines[i];

    ds[i] = line.z.Strength(line.qn, pol, pos[i].iz) *
            dline_strength_calc_dY(line.ls.dY_dX(atm, key.spec, key.ls_coeff),
                                   lshp.inv_gd,
                                   spec,
                                   line,
                                   atm);
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dY(dcut, dz, filter);
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dY(dcut, ds, f_grid[i], filter);
    }
  } else {
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dY(ds, f_grid[i], filter);
    }
  }
}

//! Sets dshape and ds and dcut and dshape
void ComputeData::dG_core_calc(const SpeciesIsotope& spec,
                               const band_shape& shp,
                               const band_data& bnd,
                               const ConstVectorView& f_grid,
                               const AtmPoint& atm,
                               const ZeemanPolarization pol,
                               const line_key& key) {
  set_filter(key);

  for (Size i : filter) {
    const auto& line = bnd.lines[pos[i].line];
    const auto& lshp = shp.lines[i];

    ds[i] = line.z.Strength(line.qn, pol, pos[i].iz) *
            dline_strength_calc_dG(line.ls.dG_dX(atm, key.spec, key.ls_coeff),
                                   lshp.inv_gd,
                                   spec,
                                   line,
                                   atm);
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dG(dcut, dz, filter);
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dG(dcut, ds, f_grid[i], filter);
    }
  } else {
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dG(ds, f_grid[i], filter);
    }
  }
}

//! Sets dshape and dz and dcut and dshape
void ComputeData::dDV_core_calc(const band_shape& shp,
                                const band_data& bnd,
                                const ConstVectorView& f_grid,
                                const AtmPoint& atm,
                                const line_key& key) {
  using Constant::h, Constant::k;

  set_filter(key);

  for (Size i : filter) {
    const auto& lshp = shp.lines[i];
    const auto& ls   = bnd.lines[pos[i].line].ls;

    const Numeric& inv_gd = lshp.inv_gd;
    const Numeric& f0     = lshp.f0;

    const Numeric d = ls.dDV_dX(atm, key.spec, key.ls_coeff);

    dz_fac[i] = -d / f0;

    ds[i] = lshp.s * dz_fac[i];

    dz[i] = -d * inv_gd;
  }

  if (bnd.cutoff != LineByLineCutoffType::None) {
    shp.dDV(dcut, ds, dz, dz_fac, filter);
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dDV(dcut, ds, dz, dz_fac, f_grid[i], filter);
    }
  } else {
    for (Size i = 0; i < f_grid.size(); i++) {
      dshape[i] = shp.dDV(ds, dz, dz_fac, f_grid[i], filter);
    }
  }
}

namespace {
void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const ConstVectorView& f_grid,
                        const SpeciesIsotope& spec,
                        const band_shape& shape,
                        const band_data& bnd,
                        const AtmPoint& atm,
                        const ZeemanPolarization pol,
                        const AtmKey& key) {
  using enum AtmKey;
  switch (key) {
    case t:
      com_data.dt_core_calc(spec, shape, bnd, f_grid, atm, pol);
      for (Size i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dscl[i] * com_data.shape[i] +
                                    com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case p: ARTS_USER_ERROR("Not implemented, pressure derivative"); break;
    case mag_u:
      if (pol == ZeemanPolarization::no) return;
      com_data.dmag_u_core_calc(shape, bnd, f_grid, atm, pol);
      for (Size i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dnpm_du,
                                com_data.scl[i] * com_data.shape[i],
                                com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case mag_v:
      if (pol == ZeemanPolarization::no) return;
      com_data.dmag_v_core_calc(shape, bnd, f_grid, atm, pol);
      for (Size i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dnpm_dv,
                                com_data.scl[i] * com_data.shape[i],
                                com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case mag_w:
      if (pol == ZeemanPolarization::no) return;
      com_data.dmag_w_core_calc(shape, bnd, f_grid, atm, pol);
      for (Size i = 0; i < f_grid.size(); i++) {
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
      for (Size i = 0; i < f_grid.size(); i++) {
        dpm[i] += zeeman::scale(com_data.npm,
                                com_data.dscl[i] * com_data.shape[i] +
                                    com_data.scl[i] * com_data.dshape[i]);
      }
      return;
  }
}

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const ConstVectorView& f_grid,
                        const SpeciesIsotope& spec,
                        const band_shape&,
                        const band_data&,
                        const AtmPoint& atm,
                        const ZeemanPolarization,
                        const SpeciesIsotope& deriv_spec) {
  if (deriv_spec != spec) return;

  const Numeric isorat = atm[spec];

  ARTS_USER_ERROR_IF(isorat == 0, "Does not support 0 for isotopologue ratios")

  for (Size i = 0; i < f_grid.size(); i++) {
    const auto dF  = com_data.scl[i] * com_data.shape[i] / isorat;
    dpm[i]        += zeeman::scale(com_data.npm, dF);
  }
}

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const ConstVectorView& f_grid,
                        const SpeciesIsotope& spec,
                        const band_shape& shape,
                        const band_data& bnd,
                        const AtmPoint& atm,
                        const ZeemanPolarization pol,
                        const SpeciesEnum& deriv_spec) {
  com_data.dVMR_core_calc(spec, shape, bnd, f_grid, atm, pol, deriv_spec);
  for (Size i = 0; i < f_grid.size(); i++) {
    dpm[i] += zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
  }
}

void compute_derivative(PropmatVectorView dpm,
                        ComputeData& com_data,
                        const ConstVectorView& f_grid,
                        const SpeciesIsotope& spec,
                        const band_shape& shape,
                        const band_data& bnd,
                        const AtmPoint& atm,
                        const ZeemanPolarization pol,
                        const line_key& deriv) {
  switch (deriv.var) {
    case LineByLineVariable::f0:
      com_data.df0_core_calc(spec, shape, bnd, f_grid, atm, pol, deriv);
      for (Size i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case LineByLineVariable::e0:
      com_data.de0_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Size i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case LineByLineVariable::a:
      com_data.da_core_calc(shape, bnd, f_grid, deriv);
      for (Size i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case LineByLineVariable::unused: return;
  }

  switch (deriv.ls_var) {
    case LineShapeModelVariable::G0:
      com_data.dG0_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Size i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case LineShapeModelVariable::D0:
      com_data.dD0_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Size i = 0; i < f_grid.size(); i++) {
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
      for (Size i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case LineShapeModelVariable::G:
      com_data.dG_core_calc(spec, shape, bnd, f_grid, atm, pol, deriv);
      for (Size i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case LineShapeModelVariable::DV:
      com_data.dDV_core_calc(shape, bnd, f_grid, atm, deriv);
      for (Size i = 0; i < f_grid.size(); i++) {
        dpm[i] +=
            zeeman::scale(com_data.npm, com_data.scl[i] * com_data.dshape[i]);
      }
      return;
    case LineShapeModelVariable::unused: return;
  }
}

void compute_derivative(PropmatVectorView,
                        ComputeData&,
                        const ConstVectorView&,
                        const SpeciesIsotope&,
                        const band_shape&,
                        const band_data&,
                        const AtmPoint&,
                        const ZeemanPolarization,
                        const auto&) {}
}  // namespace

void calculate(PropmatVectorView pm_,
               PropmatMatrixView dpm,
               ComputeData& com_data,
               const ConstVectorView f_grid_,
               const Range& f_range,
               const JacobianTargets& jac_targets,
               const QuantumIdentifier& bnd_qid,
               const band_data& bnd,
               const AtmPoint& atm,
               const ZeemanPolarization pol,
               const bool no_negative_absorption) {
  if (std::ranges::all_of(com_data.npm, [](auto& n) { return n == 0; })) return;

  PropmatVectorView pm         = pm_[f_range];
  const ConstVectorView f_grid = f_grid_[f_range];

  const Size nf = f_grid.size();
  if (nf == 0) return;

  const SpeciesIsotope spec = bnd_qid.isot;
  const Numeric fmin        = f_grid.front();
  const Numeric fmax        = f_grid.back();

  assert(jac_targets.target_count() == static_cast<Size>(dpm.nrows()) and
         f_grid_.size() == static_cast<Size>(dpm.ncols()));
  assert(nf == pm.size());

  band_shape_helper(
      com_data.lines, com_data.pos, spec, bnd, atm, fmin, fmax, pol);
  if (com_data.lines.empty()) return;

  //! Not const to save lines for reuse
  band_shape shape{std::move(com_data.lines), bnd.get_cutoff_frequency()};

  com_data.core_calc(shape, bnd, f_grid);

  for (Size i = 0; i < nf; i++) {
    const auto F = com_data.scl[i] * com_data.shape[i];
    if (no_negative_absorption and F.real() < 0) continue;
    pm[i] += zeeman::scale(com_data.npm, F);
  }

  for (auto& atm_target : jac_targets.atm) {
    std::visit(
        [&](auto& target) {
          compute_derivative(dpm[atm_target.target_pos, f_range],
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

  for (auto& line_target : jac_targets.line) {
    if (line_target.type.band == bnd_qid) {
      compute_derivative(dpm[line_target.target_pos, f_range],
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

namespace {
void compute_derivative(ComplexVectorView dp,
                        const ConstComplexVectorView&,
                        const ConstVectorView& f_grid,
                        const single_shape& shp,
                        const QuantumIdentifier& qid,
                        const line& line,
                        const AtmPoint& atm,
                        const ZeemanPolarization pol,
                        const AtmKey& key,
                        const Size,
                        const Index iz) {
  using enum AtmKey;
  switch (key) {
    case t: {
      std::transform(
          f_grid.begin(),
          f_grid.end(),
          dp.begin(),
          dp.begin(),
          [&shp,
           dz_fac = (-2 * atm.temperature * line.ls.dD0_dT(atm) -
                     2 * atm.temperature * line.ls.dDV_dT(atm) - shp.f0) /
                    (2 * atm.temperature * shp.f0),
           ds = line.z.Strength(line.qn, pol, iz) *
                dline_strength_calc_dT(shp.inv_gd, shp.f0, qid.isot, line, atm),
           dz = shp.inv_gd *
                Complex{-dline_center_calc_dT(line, atm), line.ls.dG0_dT(atm)}](
              auto f, auto d) { return shp.dT(ds, dz, dz_fac, f) + d; });

    } break;
    case p:     ARTS_USER_ERROR("Not implemented, pressure derivative"); break;
    case mag_u: {
      if (pol == ZeemanPolarization::no) return;
      const Numeric H         = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
      const Numeric dH_dmag_u = atm.mag[0] / H;
      std::transform(
          f_grid.begin(),
          f_grid.end(),
          dp.begin(),
          dp.begin(),
          [&shp,
           dz = -shp.inv_gd * dH_dmag_u * line.z.Splitting(line.qn, pol, iz)](
              auto f, auto d) { return shp.dH(dz, f) + d; });
    } break;
    case mag_v: {
      if (pol == ZeemanPolarization::no) return;
      const Numeric H         = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
      const Numeric dH_dmag_v = atm.mag[1] / H;
      std::transform(
          f_grid.begin(),
          f_grid.end(),
          dp.begin(),
          dp.begin(),
          [&shp,
           dz = -shp.inv_gd * dH_dmag_v * line.z.Splitting(line.qn, pol, iz)](
              auto f, auto d) { return shp.dH(dz, f) + d; });
    } break;
    case mag_w: {
      if (pol == ZeemanPolarization::no) return;
      const Numeric H         = std::hypot(atm.mag[0], atm.mag[1], atm.mag[2]);
      const Numeric dH_dmag_w = atm.mag[2] / H;
      std::transform(
          f_grid.begin(),
          f_grid.end(),
          dp.begin(),
          dp.begin(),
          [&shp,
           dz = -shp.inv_gd * dH_dmag_w * line.z.Splitting(line.qn, pol, iz)](
              auto f, auto d) { return shp.dH(dz, f) + d; });
    } break;
    case wind_u:
    case wind_v:
    case wind_w: {
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dp.begin(),
                     dp.begin(),
                     [&shp](auto f, auto d) { return shp.df(f) + d; });
    } break;
  }
}

void compute_derivative(ComplexVectorView dp,
                        const ConstComplexVectorView&,
                        const ConstVectorView& f_grid,
                        const single_shape& shp,
                        const QuantumIdentifier& qid,
                        const line& line,
                        const AtmPoint& atm,
                        const ZeemanPolarization pol,
                        const SpeciesEnum& key,
                        const Size,
                        const Index iz) {
  const Numeric dz_fac =
      -(line.ls.dD0_dVMR(atm, key) + line.ls.dDV_dVMR(atm, key)) / shp.f0;
  const Complex ds =
      line.z.Strength(line.qn, pol, iz) *
      dline_strength_calc_dVMR(shp.inv_gd, shp.f0, qid.isot, key, line, atm);
  const Complex dz =
      shp.inv_gd * Complex{-dline_center_calc_dVMR(line, key, atm),
                           line.ls.dG0_dVMR(atm, key)};
  std::transform(f_grid.begin(),
                 f_grid.end(),
                 dp.begin(),
                 dp.begin(),
                 [&shp, dz, dz_fac, ds](auto f, auto d) {
                   return shp.dVMR(ds, dz, dz_fac, f) + d;
                 });
}

void compute_derivative(ComplexVectorView dp,
                        const ConstVectorView& f_grid,
                        const single_shape& shp,
                        const QuantumIdentifier& qid,
                        const line& line,
                        const AtmPoint& atm,
                        const ZeemanPolarization pol,
                        const line_key& key,
                        const Size il,
                        const Index iz) {
  if (key.line != il or key.band != qid) return;

  switch (key.var) {
    case LineByLineVariable::f0: {
      const Numeric dz_fac = -1.0 / shp.f0;
      const Complex ds =
          line.z.Strength(line.qn, pol, iz) *
          dline_strength_calc_df0(shp.f0, shp.inv_gd, qid.isot, line, atm);
      const Complex dz = -shp.inv_gd;
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dp.begin(),
                     dp.begin(),
                     [&shp, dz, dz_fac, ds](auto f, auto d) {
                       return shp.df0(ds, dz, dz_fac, f) + d;
                     });
    } break;
    case LineByLineVariable::e0: {
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dp.begin(),
                     dp.begin(),
                     [&shp, ds = line.ds_de0_s_ratio(atm.temperature) * shp.s](
                         auto f, auto d) { return shp.de0(ds, f) + d; });
    } break;
    case LineByLineVariable::a: {
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dp.begin(),
                     dp.begin(),
                     [&shp, ds = shp.s / line.a](auto f, auto d) {
                       return shp.da(ds, f) + d;
                     });
    } break;
    case LineByLineVariable::unused: break;
  }

  switch (key.ls_var) {
    case LineShapeModelVariable::G0: {
      std::transform(
          f_grid.begin(),
          f_grid.end(),
          dp.begin(),
          dp.begin(),
          [&shp,
           dz = Complex(
               0, shp.inv_gd * line.ls.dG0_dX(atm, key.spec, key.ls_coeff))](
              auto f, auto d) { return shp.dG0(dz, f) + d; });
    } break;
    case LineShapeModelVariable::D0: {
      const Numeric dD0 = line.ls.dD0_dX(atm, key.spec, key.ls_coeff);
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dp.begin(),
                     dp.begin(),
                     [&shp,
                      dz_fac = -dD0 / shp.f0,
                      ds     = -dD0 * shp.s / shp.f0,
                      dz     = -dD0 * shp.inv_gd](auto f, auto d) {
                       return shp.dD0(ds, dz, dz_fac, f) + d;
                     });
    } break;
    case LineShapeModelVariable::G2:  break;
    case LineShapeModelVariable::D2:  break;
    case LineShapeModelVariable::FVC: break;
    case LineShapeModelVariable::ETA: break;
    case LineShapeModelVariable::Y:   {
      const Complex ds =
          line.z.Strength(line.qn, pol, iz) *
          dline_strength_calc_dY(line.ls.dY_dX(atm, key.spec, key.ls_coeff),
                                 shp.inv_gd,
                                 qid.isot,
                                 line,
                                 atm);
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dp.begin(),
                     dp.begin(),
                     [&shp, ds](auto f, auto d) { return shp.dY(ds, f) + d; });
    } break;
    case LineShapeModelVariable::G: {
      const Complex ds =
          line.z.Strength(line.qn, pol, iz) *
          dline_strength_calc_dG(line.ls.dG_dX(atm, key.spec, key.ls_coeff),
                                 shp.inv_gd,
                                 qid.isot,
                                 line,
                                 atm);
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dp.begin(),
                     dp.begin(),
                     [&shp, ds](auto f, auto d) { return shp.dG(ds, f) + d; });
    } break;
    case LineShapeModelVariable::DV: {
      const Numeric d      = line.ls.dDV_dX(atm, key.spec, key.ls_coeff);
      const Numeric dz_fac = -d / shp.f0;
      const Complex ds     = shp.s * dz_fac;
      const Complex dz     = -d * shp.inv_gd;
      std::transform(f_grid.begin(),
                     f_grid.end(),
                     dp.begin(),
                     dp.begin(),
                     [&shp, ds, dz, dz_fac](auto f, auto d) {
                       return shp.dDV(ds, dz, dz_fac, f) + d;
                     });
    } break;
    case LineShapeModelVariable::unused: return;
  }
}

void compute_derivative(ComplexVectorView dp,
                        const ConstComplexVectorView& p,
                        const ConstVectorView&,
                        const single_shape&,
                        const QuantumIdentifier& qid,
                        const line&,
                        const AtmPoint& atm,
                        const ZeemanPolarization,
                        const SpeciesIsotope& key,
                        const Size,
                        const Index) {
  if (key != qid.isot) return;

  const Numeric isorat = atm[key];

  ARTS_USER_ERROR_IF(isorat == 0, "Does not support 0 for isotopologue ratios")

  einsum<"i", "", "i">(dp, 1.0 / isorat, p);
}

void compute_derivative(ComplexVectorView,
                        const ConstComplexVectorView&,
                        const ConstVectorView&,
                        const single_shape&,
                        const QuantumIdentifier&,
                        const line&,
                        const AtmPoint&,
                        const ZeemanPolarization,
                        const auto&,
                        const Size,
                        const Index) {}

bool calculate_shape(ComplexVectorView pm_,
                     ComplexMatrixView dpm_,
                     const ConstVectorView f_grid_,
                     const Range& f_range,
                     const JacobianTargets& jac_targets,
                     const QuantumIdentifier& qid,
                     const band_data& bnd,
                     const AtmPoint& atm,
                     const ZeemanPolarization pol,
                     const bool no_negative_absorption) {
  assert(pm_.size() == f_grid_.size());
  assert(dpm_.ncols() == static_cast<Index>(f_grid_.size()));
  assert(dpm_.nrows() == static_cast<Index>(jac_targets.target_count()));

  const ConstVectorView f_grid = f_grid_[f_range];
  const Size nf                = f_grid.size();
  const Size nt                = dpm_.nrows();

  //! Must have to remove negative absorption later (optimization when not?)
  ComplexVector pm(nf, 0);
  ComplexMatrix dpm(nt, nf, 0);

  const auto kernel = [&](ComplexVectorView res,
                          StridedComplexMatrixView dres,
                          const Size il,
                          const ConstVectorView freqs) {
    const auto& line = bnd.lines[il];
    const Size nz    = line.z.size(line.qn, pol);

    if (nz == 0) return true;

    const single_shape shp_(qid.isot, line, atm, pol, 0);
    for (Size iz = 0; iz < nz; iz++) {
      const single_shape shp = shp_.update_iz(qid.isot, line, atm, pol, iz);
      if (shp.s.real() == 0) continue;

      std::transform(freqs.begin(),
                     freqs.end(),
                     res.begin(),
                     res.begin(),
                     [&shp](auto f, auto p) { return shp(f) + p; });

      for (auto& atm_target : jac_targets.atm) {
        assert(atm_target.target_pos < nt);
        std::visit(
            [&](auto& target) {
              compute_derivative(
                  dres[atm_target.target_pos].to_exhaustive_view(),
                  res,
                  freqs,
                  shp,
                  qid,
                  line,
                  atm,
                  pol,
                  target,
                  il,
                  iz);
            },
            atm_target.type);
      }

      for (auto& line_target : jac_targets.line) {
        assert(line_target.target_pos < nt);
        compute_derivative(dres[line_target.target_pos].to_exhaustive_view(),
                           freqs,
                           shp,
                           qid,
                           line,
                           atm,
                           pol,
                           line_target.type,
                           il,
                           iz);
      }
    }

    return false;
  };

  bool has_pol = false;
  switch (bnd.cutoff) {
    case LineByLineCutoffType::None: {
      for (Size il = 0; il < bnd.lines.size(); il++) {
        has_pol = kernel(pm, dpm, il, f_grid) or has_pol;
      }
    } break;
    case LineByLineCutoffType::ByLine: {
      Complex cutoff = 0.0;
      ComplexVector dcutoffs(nt, 0.0);
      ComplexVectorView cut{cutoff};
      ComplexMatrixView dcut = dcutoffs.view_as(nt, 1);

      for (Size il = 0; il < bnd.lines.size(); il++) {
        const Numeric l = bnd.lines[il].f0 - bnd.cutoff_value;
        const Numeric u = bnd.lines[il].f0 + bnd.cutoff_value;
        const ConstVectorView freq{u};
        const auto in = sorted_range(f_grid, l, u);
        has_pol = kernel(pm[in], dpm[joker, in], il, f_grid[in]) or has_pol;
        has_pol = kernel(cut, dcut, il, freq) or has_pol;
      }

      pm -= cutoff;
      for (Size i = 0; i < nt; i++) dpm[i] -= dcutoffs[i];
    } break;
  }

  if (no_negative_absorption) {
    for (Size i = 0; i < nf; i++) {
      if (pm[i].real() < 0) {
        pm[i]         = Complex{0, 0};
        dpm[joker, i] = 0.0;
      }
    }
  }

  pm_[f_range]         += pm;
  dpm_[joker, f_range] += dpm;

  return has_pol;
}

void multiply_scale(ComplexVectorView pm,
                    ComplexMatrixView dpm,
                    const ConstVectorView f_grid,
                    const JacobianTargets& jac_targets,
                    const AtmPoint& atm) {
  using Constant::pi, Constant::c, Constant::h, Constant::k;
  constexpr Numeric sc = c * c / (8 * pi);

  const Numeric T = atm.temperature;
  const Numeric P = atm.pressure;
  const Numeric N = number_density(P, T);

  for (Size i = 0; i < pm.size(); i++) {
    const auto f    = f_grid[i];
    const Numeric r = (h * f) / (k * T);

    const Numeric e   = std::expm1(-r);
    const Numeric scl = -N * f * e * sc;

    auto& p = pm[i];
    for (auto& atm_target : jac_targets.atm) {
      auto& dp  = dpm[atm_target.target_pos, i];
      dp       *= scl;

      std::visit(
          [&]<typename U>(const U& key) {
            if constexpr (std::same_as<U, AtmKey>) {
              switch (key) {
                case AtmKey::p: dp += p * scl / P; break;
                case AtmKey::t:
                  dp -= p * f * N * (r * (e + 1) - e) * sc / T;
                  break;
                case AtmKey::wind_u:
                case AtmKey::wind_v:
                case AtmKey::wind_w:
                  dp += p * N * (r * (e + 1) - e) * sc;
                  break;
                default:;
              }
            }
          },
          atm_target.type);
    }

    p *= scl;
  }
}
}  // namespace

bool calculate(ComplexVectorView pm,
               ComplexMatrixView dpm,
               const AbsorptionBands& bands,
               const ConstVectorView f_grid,
               const JacobianTargets& jac_targets,
               const AtmPoint& atm,
               const ZeemanPolarization pol,
               const SpeciesEnum spec,
               const bool no_negative_absorption) {
  assert(pm.size() == f_grid.size());
  assert(dpm.ncols() == static_cast<Index>(f_grid.size()));
  assert(dpm.nrows() == static_cast<Index>(jac_targets.target_count()));

  const Size n = arts_omp_get_max_threads();
  bool has_pol = false;
  if (arts_omp_parallel(pm.size() > n)) {
    const auto f_ranges = matpack::omp_offset_count(f_grid.size(), n);
    std::vector<Index> thread_has_pol(n, 0);
    std::string error;
#pragma omp parallel for
    for (Size i = 0; i < n; i++) {
      try {
        for (auto& [qid, band] : bands) {
          thread_has_pol[i] += calculate_shape(pm,
                                               dpm,
                                               f_grid,
                                               f_ranges[i],
                                               jac_targets,
                                               qid,
                                               band,
                                               atm,
                                               pol,
                                               no_negative_absorption);
        }
      } catch (std::exception& e) {
#pragma omp critical
        if (error.empty()) error = e.what();
      }
    }

    if (not error.empty()) throw std::runtime_error(error);
    has_pol = std::ranges::any_of(thread_has_pol, [](auto v) { return v; });
  } else {
    for (auto& [qid, band] : bands) {
      if (band.lineshape != LineByLineLineshape::VP_LTE) continue;
      if (spec != "AIR"_spec and spec != qid.isot.spec) continue;

      has_pol = calculate_shape(pm,
                                dpm,
                                f_grid,
                                Range{0, f_grid.size()},
                                jac_targets,
                                qid,
                                band,
                                atm,
                                pol,
                                no_negative_absorption) or
                has_pol;
    }
  }

  multiply_scale(pm, dpm, f_grid, jac_targets, atm);
  return has_pol;
}
}  // namespace lbl::voigt::lte
