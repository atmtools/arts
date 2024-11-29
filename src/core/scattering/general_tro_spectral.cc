#include "general_tro_spectral.h"

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
