#pragma once

#include <arts_constants.h>
#include <atm.h>
#include <matpack.h>
#include <rtepack.h>

#include <limits>
#include <iosfwd>
#include <vector>

#include "lbl_data.h"
#include "lbl_zeeman.h"
#include "quantum_numbers.h"

//! FIXME: These functions should be elsewhere?
namespace Jacobian {
struct Targets;
}  // namespace Jacobian

namespace lbl::voigt::nlte {
struct single_shape {
  //! Line center after all adjustments
  Numeric f0{};

  //! Inverse of the Doppler broadening factor (missing an ln2)
  Numeric inv_gd{};

  //! The imaginary part of the complex argument that goes into the Faddeeva function
  Numeric z_imag{};

  Numeric k{};

  Numeric e_ratio{};

  constexpr single_shape() = default;

  single_shape(const QuantumIdentifier&,
               const line&,
               const AtmPoint&,
               const zeeman::pol,
               const Index);

  single_shape(const QuantumIdentifier&,
               const line&,
               const AtmPoint&,
               const zeeman::pol,
               const Index,
               const Size);

  [[nodiscard]] constexpr Complex z(Numeric f) const {
    return Complex{inv_gd * (f - f0), z_imag};
  }

  [[nodiscard]] static Complex F(const Complex z_);

  [[nodiscard]] Complex F(const Numeric f) const;

  [[nodiscard]] std::pair<Complex, Complex> operator()(const Numeric f) const;

  [[nodiscard]] static Complex dF(const Complex z_, const Complex F_);

  [[nodiscard]] Complex dF(const Numeric f) const;

 private:
  struct zFdF {
    Complex z, F, dF;
    zFdF(const Complex z_);
  };

  [[nodiscard]] zFdF all(const Numeric f) const;

 public:
  [[nodiscard]] std::pair<Complex, Complex> dru(const Numeric dk_dru,
                                                const Numeric de_ratio_dru,
                                                const Numeric f) const;

  [[nodiscard]] std::pair<Complex, Complex> drl(const Numeric dk_drl,
                                                const Numeric de_ratio_drl,
                                                const Numeric f) const;

  [[nodiscard]] std::pair<Complex, Complex> df(const Numeric f) const;

  [[nodiscard]] std::pair<Complex, Complex> dH(const Complex dz_dH,
                                               const Numeric f) const;

  [[nodiscard]] std::pair<Complex, Complex> dT(const Numeric dk_dT,
                                               const Numeric de_ratio_dT,
                                               const Complex dz_dT,
                                               const Numeric dz_dT_fac,
                                               const Numeric f) const;
};

//! A band shape is a collection of single shapes.  The shapes are sorted by frequency.
struct band_shape {
  //! Line absorption shapes (lacking the f * (1 - exp(-hf/kt)) factor)
  std::vector<single_shape> lines{};
  Numeric cutoff{-1};

  [[nodiscard]] Size size() const { return lines.size(); }

  band_shape() = default;

  band_shape(std::vector<single_shape>&& ls, const Numeric cut);

  [[nodiscard]] std::pair<Complex, Complex> operator()(const Numeric f) const;

  [[nodiscard]] std::pair<Complex, Complex> df(const Numeric f) const;

  [[nodiscard]] std::pair<Complex, Complex> dH(
      const ExhaustiveConstComplexVectorView& dz_dH, const Numeric f) const;

  [[nodiscard]] std::pair<Complex, Complex> dT(
      const ExhaustiveConstVectorView& dk_dT,
      const ExhaustiveConstVectorView& de_ratio_dT,
      const ExhaustiveConstComplexVectorView& dz_dT,
      const ExhaustiveConstVectorView& dz_dT_fac,
      const Numeric f) const;

  using CutView =
      matpack::matpack_view<std::pair<Complex, Complex>, 1, false, false>;
  using CutViewConst =
      matpack::matpack_view<std::pair<Complex, Complex>, 1, true, false>;

  [[nodiscard]] std::pair<Complex, Complex> operator()(const CutViewConst& cut,
                                                       const Numeric f) const;

  void operator()(CutView cut) const;

  [[nodiscard]] std::pair<Complex, Complex> df(const CutViewConst& cut,
                                               const Numeric f) const;

  void df(CutView cut) const;

  [[nodiscard]] std::pair<Complex, Complex> dH(
      const CutViewConst& cut,
      const ExhaustiveConstComplexVectorView& dz_dH,
      const Numeric f) const;

  void dH(CutView cut, const ExhaustiveConstComplexVectorView& df0_dH) const;

  [[nodiscard]] std::pair<Complex, Complex> dT(
      const CutViewConst& cut,
      const ExhaustiveConstVectorView& dk_dT,
      const ExhaustiveConstVectorView& de_ratio_dT,
      const ExhaustiveConstComplexVectorView& dz_dT,
      const ExhaustiveConstVectorView& dz_dT_fac,
      const Numeric f) const;

  void dT(CutView cut,
          const ExhaustiveConstVectorView& dk_dT,
          const ExhaustiveConstVectorView& de_ratio_dT,
          const ExhaustiveConstComplexVectorView& dz_dT,
          const ExhaustiveConstVectorView& dz_dT_fac) const;
};

void band_shape_helper(std::vector<single_shape>& lines,
                       std::vector<line_pos>& pos,
                       const QuantumIdentifier& qid,
                       const band_data& bnd,
                       const AtmPoint& atm,
                       const Numeric fmin,
                       const Numeric fmax,
                       const zeeman::pol pol);

struct ComputeData {
  std::vector<single_shape>
      lines{};  //! Line shapes; save for reuse, assume moved from
  std::vector<line_pos> pos{};  //! Save for reuse, size of line shapes

  using PairDataC = matpack::matpack_data<std::pair<Complex, Complex>, 1>;

  PairDataC cut{};     //! Size of line shapes
  ComplexVector dz{};  //! Size of line shapes
  Vector dz_fac{};     //! Size of line shapes
  Vector dk{};         //! Size of line shapes
  Vector de_ratio{};   //! Size of line shapes
  PairDataC dcut{};    //! Size of line shapes

  Vector scl{};        //! Size of frequency
  Vector dscl{};       //! Size of frequency
  PairDataC shape{};   //! Size of frequency
  PairDataC dshape{};  //! Size of frequency

  Propmat npm{};      //! The orientation of the polarization
  Propmat dnpm_du{};  //! The orientation of the polarization
  Propmat dnpm_dv{};  //! The orientation of the polarization
  Propmat dnpm_dw{};  //! The orientation of the polarization

  //! Sizes scl, dscl, shape, dshape.  Sets scl, npm, dnpm_du, dnpm_dv, dnpm_dw
  ComputeData(const ExhaustiveConstVectorView& f_grid,
              const AtmPoint& atm,
              const Vector2& los    = {},
              const zeeman::pol pol = zeeman::pol::no);

  void update_zeeman(const Vector2& los,
                     const Vector3& mag,
                     const zeeman::pol pol);

  //! Sizes cut, dcut, dz, ds; sets shape
  void core_calc(const band_shape& shp,
                 const band_data& bnd,
                 const ExhaustiveConstVectorView& f_grid);

  //! Sets dshape and dscl and ds and dz
  void dt_core_calc(const QuantumIdentifier& qid,
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

  //! Pure debug print, will never be the same
  friend std::ostream& operator<<(std::ostream& os, const ComputeData& cd);
};

void calculate(PropmatVectorView pm,
               StokvecVectorView sv,
               matpack::matpack_view<Propmat, 2, false, true> dpm,
               matpack::matpack_view<Stokvec, 2, false, true> dsv,
               ComputeData& com_data,
               const ExhaustiveConstVectorView& f_grid,
               const Jacobian::Targets& jacobian_targets,
               const QuantumIdentifier& bnd_qid,
               const band_data& bnd,
               const AtmPoint& atm,
               const zeeman::pol pol,
               const bool no_negative_absorption);
}  // namespace lbl::voigt::nlte
