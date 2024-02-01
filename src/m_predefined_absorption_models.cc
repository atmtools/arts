/*!
 * @file   m_fullmodel.cc
 * @author Richard Larsson
 * @date   2020-01-29
 * 
 * @brief  Full absorption models of various kinds
 */

#include <predef_data.h>

#include <algorithm>

#include "atm.h"
#include "debug.h"
#include "isotopologues.h"
#include "jacobian.h"
#include "logic.h"
#include "predefined_absorption_models.h"
#include "species.h"

void propagation_matrix_predefined_model_dataInit(
    PredefinedModelData& propagation_matrix_predefined_model_data) {
  propagation_matrix_predefined_model_data = PredefinedModelData{};
}

void propagation_matrix_predefined_model_dataAddWaterMTCKD400(
    PredefinedModelData& propagation_matrix_predefined_model_data,
    const Numeric& ref_temp,
    const Numeric& ref_press,
    const Numeric& ref_h2o_vmr,
    const Vector& self_absco_ref,
    const Vector& for_absco_ref,
    const Vector& wavenumbers,
    const Vector& self_texp) {
  const auto sz = self_absco_ref.size();

  ARTS_USER_ERROR_IF(
      sz not_eq for_absco_ref.size() or sz not_eq wavenumbers.size() or
          sz not_eq self_texp.size(),
      "Mismatching size, all vector inputs must be the same length")
  ARTS_USER_ERROR_IF(sz < 4, "It makes no sense to have input shorter than 4")
  ARTS_USER_ERROR_IF(not is_regularly_increasing_within_epsilon(wavenumbers),
                     "The wavenumbers must be increasing in a regular manner")

  using Model = Absorption::PredefinedModel::MT_CKD400::WaterData;
  Model x;
  x.ref_temp = ref_temp;
  x.ref_press = ref_press;
  x.ref_h2o_vmr = ref_h2o_vmr;
  x.self_absco_ref.resize(sz);
  x.for_absco_ref.resize(sz);
  x.wavenumbers.resize(sz);
  x.self_texp.resize(sz);

  std::copy(
      self_absco_ref.begin(), self_absco_ref.end(), x.self_absco_ref.begin());
  std::copy(
      for_absco_ref.begin(), for_absco_ref.end(), x.for_absco_ref.begin());
  std::copy(wavenumbers.begin(), wavenumbers.end(), x.wavenumbers.begin());
  std::copy(self_texp.begin(), self_texp.end(), x.self_texp.begin());

  propagation_matrix_predefined_model_data
      .set<Model,
           find_species_index(Species::Species::Water, "ForeignContCKDMT400")>(
          x);
  propagation_matrix_predefined_model_data
      .set<Model,
           find_species_index(Species::Species::Water, "SelfContCKDMT400")>(x);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propagation_matrixAddPredefined(
    PropmatVector& propagation_matrix,
    PropmatMatrix& propagation_matrix_jacobian,
    const PredefinedModelData& propagation_matrix_predefined_model_data,
    const SpeciesEnum& select_species,
    const JacobianTargets& jacobian_targets,
    const AscendingGrid& f_grid,
    const AtmPoint& atm_point) {
  ARTS_USER_ERROR_IF(
      propagation_matrix.size() not_eq f_grid.size(),
      "Mismatch dimensions on internal matrices of xsec and frequency");

  // Derivatives and their error handling
  if (propagation_matrix_jacobian.nrows()) {
    ARTS_USER_ERROR_IF(
        static_cast<Size>(propagation_matrix_jacobian.nrows()) not_eq
            jacobian_targets.target_count(),
        "Mismatch dimensions on xsec derivatives and Jacobian grids");
    ARTS_USER_ERROR_IF(
        propagation_matrix_jacobian.ncols() not_eq f_grid.size(),
        "Mismatch dimensions on internal matrices of xsec derivatives and frequency");
  }

  const Absorption::PredefinedModel::VMRS vmr(atm_point);
  for (auto& [isot, data] : propagation_matrix_predefined_model_data) {
    if (select_species != SpeciesEnum::Bath and isot.spec != select_species)
      continue;
    Absorption::PredefinedModel::compute(propagation_matrix,
                                         propagation_matrix_jacobian,
                                         isot,
                                         f_grid,
                                         atm_point.pressure,
                                         atm_point.temperature,
                                         vmr,
                                         jacobian_targets,
                                         data);
  }
}
