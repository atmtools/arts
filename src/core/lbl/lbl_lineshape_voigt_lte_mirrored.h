#pragma once

#include <arts_constants.h>
#include <atm.h>
#include <matpack.h>
#include <rtepack.h>

#include <limits>
#include <vector>

#include "lbl_data.h"

//! FIXME: These functions should be elsewhere?
namespace Jacobian {
struct Targets;
}  // namespace Jacobian

namespace lbl::voigt::lte_mirror {
struct single_shape {
  //! Line center after all adjustments
  Numeric f0{};

  //! Inverse of the Doppler broadening factor (missing an ln2)
  Numeric inv_gd{};

  //! The imaginary part of the complex argument that goes into the Faddeeva function
  Numeric z_imag{};

  //! Linestrength but lacks N * f * (1 - exp(-hf/kt)) * (c^2 / 8pi) factor; also has Constant::inv_sqrt_pi * inv_gd factor
  Complex s{};

  constexpr single_shape() = default;

  single_shape(const SpeciesIsotope&,
               const line&,
               const AtmPoint&,
               const ZeemanPolarization,
               const Index);

  [[nodiscard]] constexpr Complex z(Numeric f) const {
    return Complex{inv_gd * (f - f0), z_imag};
  }

  [[nodiscard]] constexpr Complex zm(Numeric f) const {
    return Complex{inv_gd * (f + f0), z_imag};
  }

  [[nodiscard]] static Complex F(const Complex z_);

  [[nodiscard]] Complex F(const Numeric f) const;

  [[nodiscard]] Complex operator()(const Numeric f) const;

  [[nodiscard]] static Complex dF(const Complex z_, const Complex F_);

  [[nodiscard]] Complex dF(const Numeric f) const;

 private:
  struct zFdF {
    Complex zp, zm;
    Complex Fp, Fm;
    Complex dFp, dFm;
    zFdF(const Complex zp_, const Complex zm_);
  };

  [[nodiscard]] zFdF all(const Numeric f) const;

 public:
  [[nodiscard]] Complex df(const Numeric f) const;

  [[nodiscard]] Complex df0(const Complex ds_df0,
                            const Complex dz_df0,
                            const Numeric dz_df0_fac,
                            const Numeric f) const;

  [[nodiscard]] Complex dDV(const Complex ds_dDV,
                            const Complex dz_dDV,
                            const Numeric dz_dDV_fac,
                            const Numeric f) const;

  [[nodiscard]] Complex dD0(const Complex ds_dD0,
                            const Complex dz_dD0,
                            const Numeric dz_dD0_fac,
                            const Numeric f) const;

  [[nodiscard]] Complex dG0(const Complex dz_dG0, const Numeric f) const;

  [[nodiscard]] Complex dH(const Complex dz_dH, const Numeric f) const;

  [[nodiscard]] Complex dVMR(const Complex ds_dVMR,
                             const Complex dz_dVMR,
                             const Numeric dz_dVMR_fac,
                             const Numeric f) const;

  [[nodiscard]] Complex dT(const Complex ds_dT,
                           const Complex dz_dT,
                           const Numeric dz_dT_fac,
                           const Numeric f) const;

  [[nodiscard]] Complex da(const Complex ds_da, const Numeric f) const;

  [[nodiscard]] Complex de0(const Complex ds_de0, const Numeric f) const;

  [[nodiscard]] Complex dG(const Complex ds_dG, const Numeric f) const;

  [[nodiscard]] Complex dY(const Complex ds_dY, const Numeric f) const;
};

Size count_lines(const band_data& bnd, const ZeemanPolarization type);

//! Helper for initializing the band_shape
void band_shape_helper(std::vector<single_shape>& lines,
                       std::vector<line_pos>& pos,
                       const SpeciesIsotope& spec,
                       const band_data& bnd,
                       const AtmPoint& atm,
                       const Numeric fmin,
                       const Numeric fmax,
                       const ZeemanPolarization pol);

constexpr std::pair<Index, Index> find_offset_and_count_of_frequency_range(
    const std::span<const single_shape> lines, Numeric f, Numeric cutoff) {
  if (cutoff < std::numeric_limits<Numeric>::infinity()) {
    auto low =
        std::ranges::lower_bound(lines, f - cutoff, {}, &single_shape::f0);
    auto upp =
        std::ranges::upper_bound(lines, f + cutoff, {}, &single_shape::f0);

    return {std::distance(lines.begin(), low), std::distance(low, upp)};
  }

  return {0, lines.size()};
}

namespace detail {
constexpr auto frequency_span(const auto& list,
                              const Size start,
                              const Size count) {
  std::span out{list};
  return out.subspan(start, count);
}
}  // namespace detail

template <typename... Ts>
constexpr auto frequency_spans(const Numeric cutoff,
                               const Numeric f,
                               const std::span<const single_shape>& lines,
                               const Ts&... lists) {
  assert(lines.size() == (static_cast<Size>(lists.size()) and ...));

  const auto [start, count] =
      find_offset_and_count_of_frequency_range(lines, f, cutoff);
  return std::tuple{detail::frequency_span(lines, start, count),
                    detail::frequency_span(lists, start, count)...};
}

//! A band shape is a collection of single shapes.  The shapes are sorted by frequency.
struct band_shape {
  //! Line absorption shapes (lacking the f * (1 - exp(-hf/kt)) factor)
  std::vector<single_shape> lines{};
  Numeric cutoff{-1};

  [[nodiscard]] Size size() const { return lines.size(); }

  band_shape() = default;

  band_shape(std::vector<single_shape>&& ls, const Numeric cut);

  [[nodiscard]] Complex operator()(const Numeric f) const;

  [[nodiscard]] Complex df(const Numeric f) const;

  [[nodiscard]] Complex dH(const ConstComplexVectorView& dz_dH,
                           const Numeric f) const;

  [[nodiscard]] Complex dT(const ConstComplexVectorView& ds_dT,
                           const ConstComplexVectorView& dz_dT,
                           const ConstVectorView& dz_dT_fac,
                           const Numeric f) const;

  [[nodiscard]] Complex dVMR(const ConstComplexVectorView& ds_dVMR,
                             const ConstComplexVectorView& dz_dVMR,
                             const ConstVectorView& dz_dVMR_fac,
                             const Numeric f) const;

  [[nodiscard]] Complex df0(const ConstComplexVectorView ds_df0,
                            const ConstComplexVectorView dz_df0,
                            const ConstVectorView dz_df0_fac,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  [[nodiscard]] Complex da(const ConstComplexVectorView ds_da,
                           const Numeric f,
                           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex de0(const ConstComplexVectorView ds_de0,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dDV(const ConstComplexVectorView ds_dDV,
                            const ConstComplexVectorView dz_dDV,
                            const ConstVectorView dz_dDV_fac,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dD0(const ConstComplexVectorView ds_dD0,
                            const ConstComplexVectorView dz_dD0,
                            const ConstVectorView dz_dD0_fac,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dG0(const ConstComplexVectorView dz_dG0,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dY(const ConstComplexVectorView ds_dY,
                           const Numeric f,
                           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dG(const ConstComplexVectorView ds_dG,
                           const Numeric f,
                           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex operator()(const ConstComplexVectorView& cut,
                                   const Numeric f) const;

  void operator()(ComplexVectorView cut) const;

  [[nodiscard]] Complex df(const ConstComplexVectorView& cut,
                           const Numeric f) const;

  void df(ComplexVectorView cut) const;

  [[nodiscard]] Complex dH(const ConstComplexVectorView& cut,
                           const ConstComplexVectorView& dz_dH,
                           const Numeric f) const;

  void dH(ComplexVectorView cut, const ConstComplexVectorView& df0_dH) const;

  [[nodiscard]] Complex dT(const ConstComplexVectorView& cut,
                           const ConstComplexVectorView& ds_dT,
                           const ConstComplexVectorView& dz_dT,
                           const ConstVectorView& dz_dT_fac,
                           const Numeric f) const;

  void dT(ComplexVectorView cut,
          const ConstComplexVectorView& ds_dT,
          const ConstComplexVectorView& dz_dT,
          const ConstVectorView& dz_dT_fac) const;

  [[nodiscard]] Complex dVMR(const ConstComplexVectorView& cut,
                             const ConstComplexVectorView& ds_dVMR,
                             const ConstComplexVectorView& dz_dVMR,
                             const ConstVectorView& dz_dVMR_fac,
                             const Numeric f) const;

  void dVMR(ComplexVectorView cut,
            const ConstComplexVectorView& ds_dVMR,
            const ConstComplexVectorView& dz_dVMR,
            const ConstVectorView& dz_dVMR_fac) const;

  [[nodiscard]] Complex df0(const ConstComplexVectorView& cut,
                            const ConstComplexVectorView ds_df0,
                            const ConstComplexVectorView dz_df0,
                            const ConstVectorView dz_df0_fac,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  void df0(ComplexVectorView cut,
           const ConstComplexVectorView ds_df0,
           const ConstComplexVectorView dz_df0,
           const ConstVectorView dz_df0_fac,
           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex da(const ConstComplexVectorView& cut,
                           const ConstComplexVectorView ds_da,
                           const Numeric f,
                           const std::vector<Size>& filter) const;

  void da(ComplexVectorView cut,
          const ConstComplexVectorView ds_da,
          const std::vector<Size>& filter) const;

  [[nodiscard]] Complex de0(const ConstComplexVectorView& cut,
                            const ConstComplexVectorView ds_de0,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  void de0(ComplexVectorView cut,
           const ConstComplexVectorView ds_de0,
           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dDV(const ConstComplexVectorView& cut,
                            const ConstComplexVectorView ds_dDV,
                            const ConstComplexVectorView dz_dDV,
                            const ConstVectorView dz_dDV_fac,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  void dDV(ComplexVectorView cut,
           const ConstComplexVectorView ds_dDV,
           const ConstComplexVectorView dz_dDV,
           const ConstVectorView dz_dDV_fac,
           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dD0(const ConstComplexVectorView& cut,
                            const ConstComplexVectorView ds_dD0,
                            const ConstComplexVectorView dz_dD0,
                            const ConstVectorView dz_dD0_fac,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  void dD0(ComplexVectorView cut,
           const ConstComplexVectorView ds_dD0,
           const ConstComplexVectorView dz_dD0,
           const ConstVectorView dz_dD0_fac,
           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dG0(const ConstComplexVectorView& cut,
                            const ConstComplexVectorView dz_dG0,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  void dG0(ComplexVectorView cut,
           const ConstComplexVectorView dz_dG0,
           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dY(const ConstComplexVectorView& cut,
                           const ConstComplexVectorView ds_dY,
                           const Numeric f,
                           const std::vector<Size>& filter) const;

  void dY(ComplexVectorView cut,
          const ConstComplexVectorView ds_dY,
          const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dG(const ConstComplexVectorView& cut,
                           const ConstComplexVectorView ds_dG,
                           const Numeric f,
                           const std::vector<Size>& filter) const;

  void dG(ComplexVectorView cut,
          const ConstComplexVectorView ds_dG,
          const std::vector<Size>& filter) const;
};

struct ComputeData {
  std::vector<single_shape>
      lines{};  //! Line shapes; save for reuse, assume moved from
  std::vector<line_pos> pos{};  //! Save for reuse, size of line shapes

  Size filtered_line{std::numeric_limits<
      Size>::max()};  //! filter is for this and filtered_spec
  SpeciesEnum filtered_spec{
      SpeciesEnum::unused};  //! filter is for this and filtered_spec
  std::vector<Size>
      filter;  //! Filter for line parameters; resized all the time but reserves size of line shapes

  ComplexVector cut{};   //! Size of line shapes
  ComplexVector dz{};    //! Size of line shapes
  Vector dz_fac{};       //! Size of line shapes
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
  ComputeData(const ConstVectorView& f_grid,
              const AtmPoint& atm,
              const Vector2& los           = {},
              const ZeemanPolarization pol = ZeemanPolarization::no);

  void update_zeeman(const Vector2& los,
                     const Vector3& mag,
                     const ZeemanPolarization pol);

  //! Sizes cut, dcut, dz, ds; sets shape
  void core_calc(const band_shape& shp,
                 const band_data& bnd,
                 const ConstVectorView& f_grid);

  //! Sets dshape and dscl and ds and dz
  void dt_core_calc(const SpeciesIsotope& spec,
                    const band_shape& shp,
                    const band_data& bnd,
                    const ConstVectorView& f_grid,
                    const AtmPoint& atm,
                    const ZeemanPolarization pol);

  //! Sets dshape and dscl
  void df_core_calc(const band_shape& shp,
                    const band_data& bnd,
                    const ConstVectorView& f_grid,
                    const AtmPoint& atm);

  //! Sets dshape and dz
  void dmag_u_core_calc(const band_shape& shp,
                        const band_data& bnd,
                        const ConstVectorView& f_grid,
                        const AtmPoint& atm,
                        const ZeemanPolarization pol);

  //! Sets dshape and dz
  void dmag_v_core_calc(const band_shape& shp,
                        const band_data& bnd,
                        const ConstVectorView& f_grid,
                        const AtmPoint& atm,
                        const ZeemanPolarization pol);

  //! Sets dshape and dz
  void dmag_w_core_calc(const band_shape& shp,
                        const band_data& bnd,
                        const ConstVectorView& f_grid,
                        const AtmPoint& atm,
                        const ZeemanPolarization pol);

  //! Sets ds and dz and dcut and dshape
  void dVMR_core_calc(const SpeciesIsotope& spec,
                      const band_shape& shp,
                      const band_data& bnd,
                      const ConstVectorView& f_grid,
                      const AtmPoint& atm,
                      const ZeemanPolarization pol,
                      const SpeciesEnum target_spec);

  void set_filter(const line_key& key);

  //! Sets dshape and ds and dz and dcut and dshape
  void df0_core_calc(const SpeciesIsotope& spec,
                     const band_shape& shp,
                     const band_data& bnd,
                     const ConstVectorView& f_grid,
                     const AtmPoint& atm,
                     const ZeemanPolarization pol,
                     const line_key& key);

  //! Sets dshape and ds and dcut and dshape
  void de0_core_calc(const band_shape& shp,
                     const band_data& bnd,
                     const ConstVectorView& f_grid,
                     const AtmPoint& atm,
                     const line_key& key);

  //! Sets dshape and ds and dcut and dshape
  void da_core_calc(const band_shape& shp,
                    const band_data& bnd,
                    const ConstVectorView& f_grid,
                    const line_key& key);

  //! Sets dshape and dz and dcut and dshape
  void dG0_core_calc(const band_shape& shp,
                     const band_data& bnd,
                     const ConstVectorView& f_grid,
                     const AtmPoint& atm,
                     const line_key& key);

  //! Sets dshape and dz and dcut and dshape
  void dD0_core_calc(const band_shape& shp,
                     const band_data& bnd,
                     const ConstVectorView& f_grid,
                     const AtmPoint& atm,
                     const line_key& key);

  //! Sets dshape and ds and dcut and dshape
  void dY_core_calc(const SpeciesIsotope& spec,
                    const band_shape& shp,
                    const band_data& bnd,
                    const ConstVectorView& f_grid,
                    const AtmPoint& atm,
                    const ZeemanPolarization pol,
                    const line_key& key);

  //! Sets dshape and ds and dcut and dshape
  void dG_core_calc(const SpeciesIsotope& spec,
                    const band_shape& shp,
                    const band_data& bnd,
                    const ConstVectorView& f_grid,
                    const AtmPoint& atm,
                    const ZeemanPolarization pol,
                    const line_key& key);

  //! Sets dshape and dz and dcut and dshape
  void dDV_core_calc(const band_shape& shp,
                     const band_data& bnd,
                     const ConstVectorView& f_grid,
                     const AtmPoint& atm,
                     const line_key& key);
};

void calculate(PropmatVectorView pm,
               PropmatMatrixView dpm,
               ComputeData& com_data,
               const ConstVectorView f_grid,
               const Range& f_range,
               const Jacobian::Targets& jac_targets,
               const QuantumIdentifier& bnd_qid,
               const band_data& bnd,
               const AtmPoint& atm,
               const ZeemanPolarization pol,
               const bool no_negative_absorption);
}  // namespace lbl::voigt::lte_mirror
