#include "scattering_species_retval.h"

ScatteringTroSpectralVector& ScatteringTroSpectralVector::operator+=(
    const ScatteringTroSpectralVector& other) {
  if (phase_matrix.has_value() and other.phase_matrix.has_value()) {
    phase_matrix.value() += other.phase_matrix.value();
  }

  extinction_matrix += other.extinction_matrix;
  absorption_vector += other.absorption_vector;

  return *this;
}
