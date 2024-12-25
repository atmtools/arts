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
  Numeric gd_fac{};  //! Doppler broadening factor of a band

  //! Size of line shapes
  Vector pop{};
  Vector dip{};
  Vector dipr{};
  ArrayOfIndex sort{};

  //! Size of line shapes x size of line shapes
  Matrix Wimag{};

  //! [1, or broadening species]
  Vector vmrs;

  //! [1, or broadening species] x size of line shapes
  ComplexMatrix eqv_strs{};
  ComplexMatrix eqv_vals{};

  //! [1, or broadening species] x size of line shapes x size of line shapes
  ComplexTensor3 Ws{};
  ComplexTensor3 Vs{};

  //! Size of frequency
  Vector scl{};
  ComplexVector shape{};

  //! The orientation of the polarization
  Propmat npm{};

  //! Sizes scl, dscl, shape, dshape.  Sets scl, npm, dnpm_du, dnpm_dv, dnpm_dw
  ComputeData(const ConstVectorView& f_grid,
              const AtmPoint& atm,
              const Vector2& los = {},
              const zeeman::pol pol = zeeman::pol::no);

  void update_zeeman(const Vector2& los,
                     const Vector3& mag,
                     const zeeman::pol pol);

  void core_calc_eqv();
  void core_calc(const ConstVectorView& f_grid);
  void adapt_single(const QuantumIdentifier& bnd_qid,
                    const band_data& bnd,
                    const linemixing::species_data_map& rovib_data,
                    const AtmPoint& atm,
                    const bool presorted = false);
  void adapt_multi(const QuantumIdentifier& bnd_qid,
                   const band_data& bnd,
                   const linemixing::species_data_map& rovib_data,
                   const AtmPoint& atm,
                   const bool presorted = false);
};

void calculate(PropmatVectorView pm,
               matpack::strided_view_t<Propmat, 2> dpm,
               ComputeData& com_data,
               const ConstVectorView& f_grid,
               const Jacobian::Targets& jacobian_targets,
               const QuantumIdentifier& bnd_qid,
               const band_data& bnd,
               const linemixing::species_data_map& rovib_data,
               const AtmPoint& atm,
               const zeeman::pol pol,
               const bool no_negative_absorption);

void equivalent_values(ComplexTensor3View eqv_str,
                       ComplexTensor3View eqv_val,
                       ComputeData& com_data,
                       const QuantumIdentifier& bnd_qid,
                       const band_data& bnd,
                       const linemixing::species_data_map& rovib_data,
                       const AtmPoint& atm,
                       const Vector& T);
}  // namespace lbl::voigt::ecs
