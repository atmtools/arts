#include "general_tro_spectral.h"

#include <memory>

#include "sht.h"

ScatteringTroSpectralVector
ScatteringGeneralSpectralTRO::get_bulk_scattering_properties_tro_spectral(
    const AtmPoint& atm_point, const Vector& f_grid, Index degree) const {
  return f(atm_point, f_grid, degree);
}

ScatteringTroSpectralVector& ScatteringTroSpectralVector::operator+=(
    const ScatteringTroSpectralVector& other) {
  if (phase_matrix.has_value() and other.phase_matrix.has_value()) {
    phase_matrix.value() += other.phase_matrix.value();
  }

  extinction_matrix += other.extinction_matrix;
  absorption_vector += other.absorption_vector;

  return *this;
}

scattering::PhaseMatrixData<Numeric,
                            scattering::Format::TRO,
                            scattering::Representation::Spectral>
ScatteringTroSpectralVector::to_general(
    const SpecmatMatrix& phase_matrix,
    const std::shared_ptr<Vector>& f_grid_ptr) {
  ARTS_USER_ERROR_IF(not f_grid_ptr, "f grid must be provided")

  ARTS_USER_ERROR_IF(phase_matrix.nrows() != f_grid_ptr->size(),
                     "Phase matrix and f grid must have the same size")

  const auto t_grid = std::make_shared<Vector>(Vector{0.0});

  const Index l = phase_matrix.ncols() - 1;
  ARTS_USER_ERROR_IF(l < 0, "Legendre degree must be non-negative, is {}", l);

  scattering::PhaseMatrixData<Numeric,
                              scattering::Format::TRO,
                              scattering::Representation::Spectral>
      pm{t_grid, f_grid_ptr, scattering::sht::provider.get_instance_lm(l, 0)};

  for (Index f_ind = 0; f_ind < pm.npages(); ++f_ind) {
    for (Index ind = 0; ind < pm.nrows(); ++ind) {
      for (Index is = 0; is < pm.ncols(); is++) {
        pm(0, f_ind, ind, is) = phase_matrix(f_ind, ind).data[is];
      }
    }
  }

  return pm;
}

scattering::ExtinctionMatrixData<Numeric,
                                 scattering::Format::TRO,
                                 scattering::Representation::Spectral>
ScatteringTroSpectralVector::to_general(
    const PropmatVector& extinction_matrix,
    const std::shared_ptr<Vector>& f_grid_ptr) {
  ARTS_USER_ERROR_IF(not f_grid_ptr, "f grid must be provided")

  ARTS_USER_ERROR_IF(extinction_matrix.size() != f_grid_ptr->size(),
                     "Extinction matrix and f grid must have the same size")

  const auto t_grid = std::make_shared<Vector>(Vector{0.0});

  scattering::ExtinctionMatrixData<Numeric,
                                   scattering::Format::TRO,
                                   scattering::Representation::Spectral>
      emd{t_grid, f_grid_ptr};

  for (Index f_ind = 0; f_ind < emd.nrows(); ++f_ind) {
    for (Index i = 0; i < emd.ncols(); i++) {
      emd(0, f_ind, i) = extinction_matrix[f_ind][i];
    }
  }

  return emd;
}

scattering::AbsorptionVectorData<Numeric,
                                 scattering::Format::TRO,
                                 scattering::Representation::Spectral>
ScatteringTroSpectralVector::to_general(
    const StokvecVector& absorption_vector,
    const std::shared_ptr<Vector>& f_grid_ptr) {
  ARTS_USER_ERROR_IF(not f_grid_ptr, "f grid must be provided")

  ARTS_USER_ERROR_IF(absorption_vector.size() != f_grid_ptr->size(),
                     "Absorption vector and f grid must have the same size")

  const auto t_grid = std::make_shared<Vector>(Vector{0.0});

  scattering::AbsorptionVectorData<Numeric,
                                   scattering::Format::TRO,
                                   scattering::Representation::Spectral>
      av{t_grid, f_grid_ptr};

  for (Index f_ind = 0; av.nrows(); ++f_ind) {
    for (Index i = 0; i < av.ncols(); i++) {
      av(0, f_ind, i) = absorption_vector[f_ind][i];
    }
  }
  return av;
}

ScatteringTroSpectralVector::general_t ScatteringTroSpectralVector::to_general(
    const std::shared_ptr<Vector>& f) const {
  ARTS_USER_ERROR_IF(not phase_matrix.has_value(), "No phase matrix")

  if (f) {
    return {.phase_matrix      = to_general(phase_matrix.value(), f),
            .extinction_matrix = to_general(extinction_matrix, f),
            .absorption_vector = to_general(absorption_vector, f)};
  }

  const Index nf = phase_matrix.value().nrows();
  const auto fs =
      std::make_shared<Vector>(nlinspace(1, static_cast<Numeric>(nf), nf));
  return {.phase_matrix      = to_general(phase_matrix.value(), fs),
          .extinction_matrix = to_general(extinction_matrix, fs),
          .absorption_vector = to_general(absorption_vector, fs)};
}

ScatteringTroSpectralVector::gridded_t
ScatteringTroSpectralVector::to_lab_frame(
    std::shared_ptr<const Vector> /*za_inc_grid*/,
    std::shared_ptr<const Vector> /*delta_aa_grid*/,
    std::shared_ptr<const scattering::ZenithAngleGrid> /*za_scat_grid_new*/)
    const {
  // return to_general().to_lab_frame(std::move(za_inc_grid),
  //                                  std::move(delta_aa_grid),
  //                                  std::move(za_scat_grid_new));

  throw std::runtime_error("Not implemented.");
  std::unreachable();
}
