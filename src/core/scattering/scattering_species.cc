#include "scattering_species.h"

#include <stdexcept>
#include <utility>

void ArrayOfScatteringSpecies::add(const scattering::Species& species_) { species.push_back(species_); }

void ArrayOfScatteringSpecies::prepare_scattering_data(scattering::ScatteringDataSpec) {}

BulkScatteringProperties<scattering::Format::TRO, scattering::Representation::Gridded>
ArrayOfScatteringSpecies::get_bulk_scattering_properties_tro_gridded(
    const AtmPoint& atm_point, const Vector& f_grid, std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const {
  if (species.size() == 0) return {std::nullopt, {}, {}};

  const auto visitor =
      [&](const auto& spec) -> BulkScatteringProperties<scattering::Format::TRO, scattering::Representation::Gridded> {
    if constexpr (requires { spec.get_bulk_scattering_properties_tro_gridded(atm_point, f_grid, za_scat_grid); }) {
      return spec.get_bulk_scattering_properties_tro_gridded(atm_point, f_grid, za_scat_grid);
    } else {
      throw std::runtime_error(std::format("Method not implemented for TRO Gridded for species:\n{:N}", spec));
    }

    std::unreachable();
  };

  auto& scat_spec = species[0];
  auto  bsp       = std::visit(visitor, scat_spec);
  for (Size ind = 1; ind < species.size(); ++ind) {
    auto& scat_spec  = species[ind];
    bsp             += std::visit(visitor, scat_spec);
  }
  return bsp;
}

BulkScatteringProperties<scattering::Format::TRO, scattering::Representation::Gridded>
ArrayOfScatteringSpecies::get_bulk_scattering_properties_tro_gridded_derivative(
    const AtmPoint&                              atm_point,
    const Vector&                                f_grid,
    std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid,
    const AtmKeyVal&                             target) const {
  ARTS_USER_ERROR_IF(species.empty(), "Cannot differentiate an empty scattering-species array")
  const auto visitor =
      [&](const auto& spec) -> BulkScatteringProperties<scattering::Format::TRO, scattering::Representation::Gridded> {
    if constexpr (requires {
                    spec.get_bulk_scattering_properties_tro_gridded_derivative(atm_point, f_grid, za_scat_grid, target);
                  }) {
      return spec.get_bulk_scattering_properties_tro_gridded_derivative(atm_point, f_grid, za_scat_grid, target);
    } else {
      throw std::runtime_error(std::format("TRO-gridded derivatives are not implemented for species:\n{:N}", spec));
    }
    std::unreachable();
  };
  auto out = std::visit(visitor, species.front());
  for (Size i = 1; i < species.size(); ++i) out += std::visit(visitor, species[i]);
  return out;
}

ScatteringTroSpectralVector ArrayOfScatteringSpecies::get_bulk_scattering_properties_tro_spectral(
    const AtmPoint& atm_point, const Vector& f_grid, Index degree) const {
  if (species.size() == 0) return {std::nullopt, {}, {}};

  ARTS_USER_ERROR_IF(degree < 0, "The Legendre degree must be >= 0, got {}", degree)
  const auto visitor = [&](const auto& spec) -> ScatteringTroSpectralVector {
    if constexpr (requires { spec.get_bulk_scattering_properties_tro_spectral(atm_point, f_grid, degree); }) {
      return spec.get_bulk_scattering_properties_tro_spectral(atm_point, f_grid, degree);
    } else {
      throw std::runtime_error(std::format("Method not implemented for TRO Spectral for species:\n{:N}", spec));
    }

    std::unreachable();
  };

  // Every species, a user's function too, must give exactly what was asked
  const Size nf    = f_grid.size();
  const auto check = [&](const ScatteringTroSpectralVector& v, Size ind) {
    ARTS_USER_ERROR_IF(not v.phase_matrix.has_value(), "Scattering species {} gives no Legendre series", ind)
    ARTS_USER_ERROR_IF(v.phase_matrix->nrows() != static_cast<Index>(nf) or v.phase_matrix->ncols() != degree + 1,
                       "Scattering species {} gives Legendre series of shape [{}, {}], but [{} frequencies, degree "
                       "{} + 1] were asked for",
                       ind,
                       v.phase_matrix->nrows(),
                       v.phase_matrix->ncols(),
                       nf,
                       degree)
    ARTS_USER_ERROR_IF(v.extinction_matrix.size() != nf or v.absorption_vector.size() != nf,
                       "Scattering species {} gives {} extinction matrices and {} absorption vectors for {} "
                       "frequencies",
                       ind,
                       v.extinction_matrix.size(),
                       v.absorption_vector.size(),
                       nf)
  };

  auto bsp = std::visit(visitor, species[0]);
  check(bsp, 0);
  for (Size ind = 1; ind < species.size(); ++ind) {
    auto next = std::visit(visitor, species[ind]);
    check(next, ind);
    bsp += next;
  }
  return bsp;
}

BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded>
ArrayOfScatteringSpecies::get_bulk_scattering_properties_aro_gridded(
    const AtmPoint&                              atm_point,
    const Vector&                                f_grid,
    const Vector&                                za_inc_grid,
    const Vector&                                delta_aa_grid,
    std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const {
  if (species.size() == 0) return {std::nullopt, {}, {}};

  const auto visitor =
      [&](const auto& spec) -> BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded> {
    if constexpr (requires {
                    spec.get_bulk_scattering_properties_aro_gridded(
                        atm_point, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);
                  }) {
      return spec.get_bulk_scattering_properties_aro_gridded(
          atm_point, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);
    } else {
      throw std::runtime_error(std::format("Method not implemented for ARO Gridded for species:\n{:N}", spec));
    }

    std::unreachable();
  };

  auto& scat_spec = species[0];
  auto  bsp       = std::visit(visitor, scat_spec);
  for (Size ind = 1; ind < species.size(); ++ind) {
    auto& scat_spec  = species[ind];
    bsp             += std::visit(visitor, scat_spec);
  }
  return bsp;
}

BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded>
ArrayOfScatteringSpecies::get_bulk_scattering_properties_aro_gridded_derivative(
    const AtmPoint&                              atm_point,
    const Vector&                                f_grid,
    const Vector&                                za_inc_grid,
    const Vector&                                delta_aa_grid,
    std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid,
    const AtmKeyVal&                             target) const {
  ARTS_USER_ERROR_IF(species.empty(), "Cannot differentiate an empty scattering-species array")
  const auto visitor =
      [&](const auto& spec) -> BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded> {
    if constexpr (requires {
                    spec.get_bulk_scattering_properties_aro_gridded_derivative(
                        atm_point, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid, target);
                  }) {
      return spec.get_bulk_scattering_properties_aro_gridded_derivative(
          atm_point, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid, target);
    } else {
      throw std::runtime_error(std::format("ARO-gridded derivatives are not implemented for species:\n{:N}", spec));
    }
    std::unreachable();
  };
  auto out = std::visit(visitor, species.front());
  for (Size i = 1; i < species.size(); ++i) out += std::visit(visitor, species[i]);
  return out;
}

BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Spectral>
ArrayOfScatteringSpecies::get_bulk_scattering_properties_aro_spectral(const AtmPoint& atm_point,
                                                                      const Vector&   f_grid,
                                                                      const Vector&   za_inc_grid,
                                                                      const Vector&   za_scat_grid,
                                                                      Index           max_mode) const {
  if (species.size() == 0) return {.phase_matrix = std::nullopt, .extinction_matrix = {}, .absorption_vector = {}};

  const auto visitor =
      [&](const auto& spec) -> BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Spectral> {
    if constexpr (requires {
                    spec.get_bulk_scattering_properties_aro_spectral(
                        atm_point, f_grid, za_inc_grid, za_scat_grid, max_mode);
                  }) {
      return spec.get_bulk_scattering_properties_aro_spectral(atm_point, f_grid, za_inc_grid, za_scat_grid, max_mode);
    } else {
      throw std::runtime_error(std::format("Method not implemented for ARO Spectral for species:\n{:N}", spec));
    }

    std::unreachable();
  };

  auto& scat_spec = species[0];
  auto  bsp       = std::visit(visitor, scat_spec);
  for (Size ind = 1; ind < species.size(); ++ind) {
    auto& scat_spec  = species[ind];
    bsp             += std::visit(visitor, scat_spec);
  }
  return bsp;
}
