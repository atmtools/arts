#pragma once

#include <rtepack.h>

#include "array.h"
#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"

//! FIXME: These functions should be elsewhere?
namespace Jacobian {
struct Targets;
}  // namespace Jacobian

namespace lbl::voigt::ecs {
struct ComputeData {
  Numeric freq_offset{};  //! Frequency offset factor
  Numeric gd_fac{};       //! Doppler broadening factor of a band

  Vector pop{};         //! Size of line shapes
  Vector dip{};         //! Size of line shapes
  Vector dipr{};        //! Size of line shapes
  ArrayOfIndex sort{};  //! Size of line shapes
  Matrix Wimag{};       //! Size of line shapes

  ArrayOfNumeric vmrs;  //! [1, or broadening species]

  ArrayOfComplexVector
      eqv_strs{};  //! [1, or broadening species] x size of line shapes
  ArrayOfComplexVector
      eqv_vals{};  //! [1, or broadening species] x size of line shapes
  ArrayOfComplexMatrix
      Ws{};  //! [1, or broadening species] x size of line shapes x line shapes
  ArrayOfComplexMatrix
      Vs{};  //! [1, or broadening species] x size of line shapes x line shapes

  Vector scl{};           //! Size of frequency
  ComplexVector shape{};  //! Size of frequency

  Propmat npm{};  //! The orientation of the polarization

  //! Sizes scl, dscl, shape, dshape.  Sets scl, npm, dnpm_du, dnpm_dv, dnpm_dw
  ComputeData(const ExhaustiveConstVectorView& f_grid,
              const AtmPoint& atm,
              const Vector2& los = {},
              const zeeman::pol pol = zeeman::pol::no);

  void update_zeeman(const Vector2& los,
                     const Vector3& mag,
                     const zeeman::pol pol);

  void core_calc_eqv();
  void core_calc(const ExhaustiveConstVectorView& f_grid);
  void adapt(const QuantumIdentifier& bnd_qid,
             const band_data& bnd,
             const linemixing::species_data_map& rovib_data,
             const AtmPoint& atm);
};

void calculate(PropmatVectorView pm,
               matpack::matpack_view<Propmat, 2, false, true> dpm,
               ComputeData& com_data,
               const ExhaustiveConstVectorView& f_grid,
               const Jacobian::Targets& jacobian_targets,
               const QuantumIdentifier& bnd_qid,
               const band_data& bnd,
               const linemixing::species_data_map& rovib_data,
               const AtmPoint& atm,
               const zeeman::pol pol = zeeman::pol::no);
}  // namespace lbl::voigt::ecs
