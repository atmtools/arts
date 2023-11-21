#pragma once

#include <arts_constants.h>
#include <atm.h>
#include <matpack.h>
#include <rtepack.h>

#include <algorithm>
#include <limits>
#include <numeric>
#include <type_traits>
#include <vector>

#include "lbl_data.h"
#include "lbl_lineshape_model.h"
#include "lbl_zeeman.h"

//! FIXME: These functions should be elsewhere?
namespace Jacobian {
struct Targets;
}  // namespace Jacobian

namespace lbl::voigt::lte {
struct single_shape {
  //! Line center after all adjustments, must be as_zeeman if zeeman effect is intended
  Numeric f0{};

  //! Inverse of the Doppler broadening factor (missing an ln2)
  Numeric inv_gd{};

  //! The imaginary part of the complex argument that goes into the Faddeeva function
  Numeric z_imag{};

  //! Linestrength but lacks N * f * (1 - exp(-hf/kt)) factor, must be as_zeeman if zeeman effect is intended, also has Constant::inv_sqrt_pi * inv_gd factor
  Complex s{};

  constexpr single_shape() = default;

  single_shape(const SpeciesIsotopeRecord& spec, const line&, const AtmPoint&);

  single_shape(const SpeciesIsotopeRecord& spec,
               const line&,
               const AtmPoint&,
               const Size);

  void as_zeeman(const line& line,
                 const Numeric H,
                 const zeeman::pol type,
                 const Index iz);

  [[nodiscard]] constexpr Complex z(Numeric f) const noexcept {
    return Complex{inv_gd * (f - f0), z_imag};
  }

  [[nodiscard]] static Complex F(const Complex z_) noexcept;

  [[nodiscard]] Complex F(const Numeric f) const noexcept;

  [[nodiscard]] Complex operator()(const Numeric f) const noexcept;

  [[nodiscard]] static Complex dF(const Complex z_, const Complex F_) noexcept;

  [[nodiscard]] Complex dF(const Numeric f) const noexcept;

 private:
  struct zFdF {
    Complex z, F, dF;
    zFdF(const Complex z_) noexcept;
  };

  [[nodiscard]] zFdF all(const Numeric f) const noexcept;

 public:
  [[nodiscard]] Complex df(const Numeric f) const noexcept;

  [[nodiscard]] Complex df0(const Complex ds_df0,
                            const Complex dz_df0,
                            const Numeric dz_df0_fac,
                            const Numeric f) const noexcept;

  [[nodiscard]] Complex dDV(const Complex dz_dDV,
                            const Numeric f) const noexcept;

  [[nodiscard]] Complex dD0(const Complex dz_dD0,
                            const Numeric f) const noexcept;

  [[nodiscard]] Complex dG0(const Complex dz_dG0,
                            const Numeric f) const noexcept;

  [[nodiscard]] Complex dH(const Complex dz_dH, const Numeric f) const noexcept;

  [[nodiscard]] Complex dVMR(const Complex ds_dVMR,
                             const Complex dz_dVMR,
                             const Numeric dz_dVMR_fac,
                             const Numeric f) const noexcept;

  [[nodiscard]] Complex dT(const Complex ds_dT,
                           const Complex dz_dT,
                           const Numeric dz_dT_fac,
                           const Numeric f) const noexcept;

  [[nodiscard]] Complex da(const Complex ds_da, const Numeric f) const noexcept;

  [[nodiscard]] Complex de0(const Complex ds_de0,
                            const Numeric f) const noexcept;

  [[nodiscard]] Complex dG(const Complex ds_dG, const Numeric f) const noexcept;

  [[nodiscard]] Complex dY(const Complex ds_dY, const Numeric f) const noexcept;
};

struct line_pos {
  Size line;
  Size spec{std::numeric_limits<Size>::max()};
  Size iz{std::numeric_limits<Size>::max()};
};

Size count_lines(const band_data& bnd, const zeeman::pol type);

void zeeman_set_back(std::vector<single_shape>& lines,
                     std::vector<line_pos>& pos,
                     const single_shape& s,
                     const line& line,
                     const Numeric H,
                     const Size spec,
                     const zeeman::pol pol,
                     Size& last_single_shape_pos);

void lines_set(std::vector<single_shape>& lines,
               std::vector<line_pos>& pos,
               const SpeciesIsotopeRecord& spec,
               const line& line,
               const AtmPoint& atm,
               const zeeman::pol pol,
               Size& last_single_shape_pos);

//! Helper for initializing the band_shape
void band_shape_helper(std::vector<single_shape>& lines,
                       std::vector<line_pos>& pos,
                       const SpeciesIsotopeRecord& spec,
                       const band_data& bnd,
                       const AtmPoint& atm,
                       const Numeric fmin,
                       const Numeric fmax,
                       const zeeman::pol pol);

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
  ARTS_ASSERT(lines.size() == (static_cast<Size>(lists.size()) and ...))

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

  [[nodiscard]] Size size() const noexcept { return lines.size(); }

  band_shape() = default;

  band_shape(std::vector<single_shape>&& ls, const Numeric cut) noexcept;

  [[nodiscard]] Complex operator()(const Numeric f) const noexcept;

  [[nodiscard]] Complex df(const Numeric f) const;

  [[nodiscard]] Complex dH(const ExhaustiveConstComplexVectorView& dz_dH,
                           const Numeric f) const;

  [[nodiscard]] Complex dT(const ExhaustiveConstComplexVectorView& ds_dT,
                           const ExhaustiveConstComplexVectorView& dz_dT,
                           const ExhaustiveConstVectorView& dz_dT_fac,
                           const Numeric f) const;

  [[nodiscard]] Complex dVMR(const ExhaustiveConstComplexVectorView& ds_dVMR,
                             const ExhaustiveConstComplexVectorView& dz_dVMR,
                             const ExhaustiveConstVectorView& dz_dVMR_fac,
                             const Numeric f) const;

  [[nodiscard]] Complex df0(const ExhaustiveConstComplexVectorView ds_df0,
                            const ExhaustiveConstComplexVectorView dz_df0,
                            const ExhaustiveConstVectorView dz_df0_fac,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  [[nodiscard]] Complex da(const ExhaustiveConstComplexVectorView ds_da,
                           const Numeric f,
                           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex de0(const ExhaustiveConstComplexVectorView ds_de0,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dDV(const ExhaustiveConstComplexVectorView dz_dDV,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dD0(const ExhaustiveConstComplexVectorView dz_dD0,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dG0(const ExhaustiveConstComplexVectorView dz_dG0,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dY(const ExhaustiveConstComplexVectorView ds_dY,
                           const Numeric f,
                           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dG(const ExhaustiveConstComplexVectorView ds_dG,
                           const Numeric f,
                           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex operator()(const ExhaustiveConstComplexVectorView& cut,
                                   const Numeric f) const;

  void operator()(ExhaustiveComplexVectorView cut) const;

  [[nodiscard]] Complex df(const ExhaustiveConstComplexVectorView& cut,
                           const Numeric f) const;

  void df(ExhaustiveComplexVectorView cut) const;

  [[nodiscard]] Complex dH(const ExhaustiveConstComplexVectorView& cut,
                           const ExhaustiveConstComplexVectorView& dz_dH,
                           const Numeric f) const;

  void dH(ExhaustiveComplexVectorView cut,
          const ExhaustiveConstComplexVectorView& df0_dH) const;

  [[nodiscard]] Complex dT(const ExhaustiveConstComplexVectorView& cut,
                           const ExhaustiveConstComplexVectorView& ds_dT,
                           const ExhaustiveConstComplexVectorView& dz_dT,
                           const ExhaustiveConstVectorView& dz_dT_fac,
                           const Numeric f) const;

  void dT(ExhaustiveComplexVectorView cut,
          const ExhaustiveConstComplexVectorView& ds_dT,
          const ExhaustiveConstComplexVectorView& dz_dT,
          const ExhaustiveConstVectorView& dz_dT_fac) const;

  [[nodiscard]] Complex dVMR(const ExhaustiveConstComplexVectorView& cut,
                             const ExhaustiveConstComplexVectorView& ds_dVMR,
                             const ExhaustiveConstComplexVectorView& dz_dVMR,
                             const ExhaustiveConstVectorView& dz_dVMR_fac,
                             const Numeric f) const;

  void dVMR(ExhaustiveComplexVectorView cut,
            const ExhaustiveConstComplexVectorView& ds_dVMR,
            const ExhaustiveConstComplexVectorView& dz_dVMR,
            const ExhaustiveConstVectorView& dz_dVMR_fac) const;

  [[nodiscard]] Complex df0(const ExhaustiveConstComplexVectorView& cut,
                            const ExhaustiveConstComplexVectorView ds_df0,
                            const ExhaustiveConstComplexVectorView dz_df0,
                            const ExhaustiveConstVectorView dz_df0_fac,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  void df0(ExhaustiveComplexVectorView cut,
           const ExhaustiveConstComplexVectorView ds_df0,
           const ExhaustiveConstComplexVectorView dz_df0,
           const ExhaustiveConstVectorView dz_df0_fac,
           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex da(const ExhaustiveConstComplexVectorView& cut,
                           const ExhaustiveConstComplexVectorView ds_da,
                           const Numeric f,
                           const std::vector<Size>& filter) const;

  void da(ExhaustiveComplexVectorView cut,
          const ExhaustiveComplexVectorView ds_da,
          const std::vector<Size>& filter) const;

  [[nodiscard]] Complex de0(const ExhaustiveConstComplexVectorView& cut,
                            const ExhaustiveConstComplexVectorView ds_de0,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  void de0(ExhaustiveComplexVectorView cut,
           const ExhaustiveComplexVectorView ds_de0,
           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dDV(const ExhaustiveConstComplexVectorView& cut,
                            const ExhaustiveConstComplexVectorView dz_dDV,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  void dDV(ExhaustiveComplexVectorView cut,
           const ExhaustiveComplexVectorView dz_dDV,
           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dD0(const ExhaustiveConstComplexVectorView& cut,
                            const ExhaustiveConstComplexVectorView dz_dD0,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  void dD0(ExhaustiveComplexVectorView cut,
           const ExhaustiveComplexVectorView dz_dD0,
           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dG0(const ExhaustiveConstComplexVectorView& cut,
                            const ExhaustiveConstComplexVectorView dz_dG0,
                            const Numeric f,
                            const std::vector<Size>& filter) const;

  void dG0(ExhaustiveComplexVectorView cut,
           const ExhaustiveComplexVectorView dz_dG0,
           const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dY(const ExhaustiveConstComplexVectorView& cut,
                           const ExhaustiveConstComplexVectorView ds_dY,
                           const Numeric f,
                           const std::vector<Size>& filter) const;

  void dY(ExhaustiveComplexVectorView cut,
          const ExhaustiveComplexVectorView ds_dY,
          const std::vector<Size>& filter) const;

  [[nodiscard]] Complex dG(const ExhaustiveConstComplexVectorView& cut,
                           const ExhaustiveConstComplexVectorView ds_dG,
                           const Numeric f,
                           const std::vector<Size>& filter) const;

  void dG(ExhaustiveComplexVectorView cut,
          const ExhaustiveComplexVectorView ds_dG,
          const std::vector<Size>& filter) const;
};

struct ComputeData {
  std::vector<single_shape>
      lines{};  //! Line shapes; save for reuse, assume moved from
  std::vector<line_pos> pos{};  //! Save for reuse, size of line shapes

  Size filtered_line{std::numeric_limits<
      Size>::max()};  //! filter is for this and filtered_spec
  Size filtered_spec{std::numeric_limits<
      Size>::max()};  //! filter is for this and filtered_spec
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
  ComputeData(const ExhaustiveConstVectorView& f_grid,
              const AtmPoint& atm,
              const Vector2& los = {},
              const zeeman::pol pol = zeeman::pol::no);

  void update_zeeman(const Vector2& los,
                     const Vector3& mag,
                     const zeeman::pol pol);

  //! Sizes cut, dcut, dz, ds; sets shape
  void core_calc(const band_shape& shp,
                 const band_data& bnd,
                 const ExhaustiveConstVectorView& f_grid);

  //! Sets dshape and dscl and ds and dz
  void dt_core_calc(const SpeciesIsotopeRecord& spec,
                    const band_shape& shp,
                    const band_data& bnd,
                    const ExhaustiveConstVectorView& f_grid,
                    const AtmPoint& atm,
                    const zeeman::pol pol);

  //! Sets dshape and dscl
  void df_core_calc(const band_shape& shp,
                    const band_data& bnd,
                    const ExhaustiveConstVectorView& f_grid,
                    const AtmPoint& atm);

  //! Sets dshape and dz
  void dmag_u_core_calc(const band_shape& shp,
                        const band_data& bnd,
                        const ExhaustiveConstVectorView& f_grid,
                        const AtmPoint& atm,
                        const zeeman::pol pol);

  //! Sets dshape and dz
  void dmag_v_core_calc(const band_shape& shp,
                        const band_data& bnd,
                        const ExhaustiveConstVectorView& f_grid,
                        const AtmPoint& atm,
                        const zeeman::pol pol);

  //! Sets dshape and dz
  void dmag_w_core_calc(const band_shape& shp,
                        const band_data& bnd,
                        const ExhaustiveConstVectorView& f_grid,
                        const AtmPoint& atm,
                        const zeeman::pol pol);

  //! Sets ds and dz and dcut and dshape
  void dVMR_core_calc(const SpeciesIsotopeRecord& spec,
                      const band_shape& shp,
                      const band_data& bnd,
                      const ExhaustiveConstVectorView& f_grid,
                      const AtmPoint& atm,
                      const zeeman::pol pol,
                      const Species::Species target_spec);

  void set_filter(const line_key& key);

  //! Sets dshape and ds and dz and dcut and dshape
  void df0_core_calc(const SpeciesIsotopeRecord& spec,const band_shape& shp,
                     const band_data& bnd,
                     const ExhaustiveConstVectorView& f_grid,
                     const AtmPoint& atm,
                     const zeeman::pol pol,
                     const line_key& key);

  //! Sets dshape and ds and dcut and dshape
  void de0_core_calc(const band_shape& shp,
                     const band_data& bnd,
                     const ExhaustiveConstVectorView& f_grid,
                     const AtmPoint& atm,
                     const line_key& key);

  //! Sets dshape and ds and dcut and dshape
  void da_core_calc(const band_shape& shp,
                    const band_data& bnd,
                    const ExhaustiveConstVectorView& f_grid,
                    const line_key& key);

  //! Sets dshape and dz and dcut and dshape
  void dG0_core_calc(const band_shape& shp,
                     const band_data& bnd,
                     const ExhaustiveConstVectorView& f_grid,
                     const AtmPoint& atm,
                     const line_key& key);

  //! Sets dshape and dz and dcut and dshape
  void dD0_core_calc(const band_shape& shp,
                     const band_data& bnd,
                     const ExhaustiveConstVectorView& f_grid,
                     const AtmPoint& atm,
                     const line_key& key);

  //! Sets dshape and ds and dcut and dshape
  void dY_core_calc(const band_shape& shp,
                    const band_data& bnd,
                    const ExhaustiveConstVectorView& f_grid,
                    const AtmPoint& atm,
                    const line_key& key);

  //! Sets dshape and ds and dcut and dshape
  void dG_core_calc(const band_shape& shp,
                    const band_data& bnd,
                    const ExhaustiveConstVectorView& f_grid,
                    const AtmPoint& atm,
                    const line_key& key);

  //! Sets dshape and dz and dcut and dshape
  void dDV_core_calc(const band_shape& shp,
                     const band_data& bnd,
                     const ExhaustiveConstVectorView& f_grid,
                     const AtmPoint& atm,
                     const line_key& key);
};

void calculate(PropmatVectorView pm,
               matpack::matpack_view<Propmat, 2, false, true> dpm,
               ComputeData& com_data,
               const ExhaustiveConstVectorView& f_grid,
               const Jacobian::Targets& jacobian_targets,
               const QuantumIdentifier& bnd_qid,
               const band_data& bnd,
               const AtmPoint& atm,
               const zeeman::pol pol = zeeman::pol::no);
}  // namespace lbl::voigt::lte
