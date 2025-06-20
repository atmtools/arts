#include "bulk_scattering_properties.h"

namespace scattering {

BulkScatteringPropertiesTROGridded
BulkScatteringPropertiesTROGridded::operator+(
    const BulkScatteringPropertiesTROGridded& other) const {
  std::optional<PhaseMatrixVector> new_phase_matrix = phase_matrix;
  if (phase_matrix.has_value()) {
    if (other.phase_matrix.has_value()) {
      *new_phase_matrix += other.phase_matrix.value();
    }
  }
  auto new_extinction_matrix  = extinction_matrix;
  new_extinction_matrix      += other.extinction_matrix;
  auto new_absorption_vector  = absorption_vector;
  new_absorption_vector      += other.absorption_vector;
  return BulkScatteringPropertiesTROGridded(
      new_phase_matrix, new_extinction_matrix, new_absorption_vector);
}

BulkScatteringPropertiesTROGridded&
BulkScatteringPropertiesTROGridded::operator+=(
    const BulkScatteringPropertiesTROGridded& other) {
  if (phase_matrix.has_value()) {
    if (other.phase_matrix.has_value()) {
      *phase_matrix += other.phase_matrix.value();
    }
  }
  extinction_matrix += other.extinction_matrix;
  absorption_vector += other.absorption_vector;
  return *this;
}

BulkScatteringPropertiesTROGridded&
BulkScatteringPropertiesTROGridded::operator*=(Numeric fac) {
  if (phase_matrix.has_value()) {
    *phase_matrix *= fac;
  }
  extinction_matrix *= fac;
  absorption_vector *= fac;
  return *this;
}
BulkScatteringPropertiesTROSpectral
BulkScatteringPropertiesTROSpectral::operator+(
    const BulkScatteringPropertiesTROSpectral& other) const {
  std::optional<PhaseMatrixVector> new_phase_matrix = phase_matrix;
  if (phase_matrix.has_value()) {
    if (other.phase_matrix.has_value()) {
      *new_phase_matrix += other.phase_matrix.value();
    }
  }
  auto new_extinction_matrix  = extinction_matrix;
  new_extinction_matrix      += other.extinction_matrix;
  auto new_absorption_vector  = absorption_vector;
  new_absorption_vector      += other.absorption_vector;
  return BulkScatteringPropertiesTROSpectral(
      new_phase_matrix, new_extinction_matrix, new_absorption_vector);
}

BulkScatteringPropertiesTROSpectral&
BulkScatteringPropertiesTROSpectral::operator+=(
    const BulkScatteringPropertiesTROSpectral& other) {
  if (phase_matrix.has_value()) {
    if (other.phase_matrix.has_value()) {
      *phase_matrix += other.phase_matrix.value();
    }
  }
  extinction_matrix += other.extinction_matrix;
  absorption_vector += other.absorption_vector;
  return *this;
}

BulkScatteringPropertiesTROSpectral&
BulkScatteringPropertiesTROSpectral::operator*=(Numeric fac) {
  if (phase_matrix.has_value()) {
    *phase_matrix *= fac;
  }
  extinction_matrix *= fac;
  absorption_vector *= fac;
  return *this;
}
}  // namespace scattering
