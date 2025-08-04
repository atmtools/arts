#include "scattering_species.h"

#include <stdexcept>
#include <utility>

void ArrayOfScatteringSpecies::add(const scattering::Species& species_) {
  species.push_back(species_);
}

void ArrayOfScatteringSpecies::prepare_scattering_data(
    scattering::ScatteringDataSpec) {}

BulkScatteringProperties<scattering::Format::TRO,
                         scattering::Representation::Gridded>
ArrayOfScatteringSpecies::get_bulk_scattering_properties_tro_gridded(
    const AtmPoint& atm_point,
    const Vector& f_grid,
    std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const {
  if (species.size() == 0) return {std::nullopt, {}, {}};

  const auto visitor = [&](const auto& spec)
      -> BulkScatteringProperties<scattering::Format::TRO,
                                  scattering::Representation::Gridded> {
    if constexpr (requires {
                    spec.get_bulk_scattering_properties_tro_gridded(
                        atm_point, f_grid, za_scat_grid);
                  }) {
      return spec.get_bulk_scattering_properties_tro_gridded(
          atm_point, f_grid, za_scat_grid);
    } else {
      throw std::runtime_error(std::format(
          "Method not implemented for TRO Gridded for species:\n{:N}", spec));
    }

    std::unreachable();
  };

  auto& scat_spec = species[0];
  auto bsp        = std::visit(visitor, scat_spec);
  for (Size ind = 1; ind < species.size(); ++ind) {
    auto& scat_spec  = species[ind];
    bsp             += std::visit(visitor, scat_spec);
  }
  return bsp;
}

ScatteringTroSpectralVector
ArrayOfScatteringSpecies::get_bulk_scattering_properties_tro_spectral(
    const AtmPoint& atm_point, const Vector& f_grid, Index degree) const {
  if (species.size() == 0) return {std::nullopt, {}, {}};

  const auto visitor = [&](const auto& spec) -> ScatteringTroSpectralVector {
    if constexpr (requires {
                    spec.get_bulk_scattering_properties_tro_spectral(
                        atm_point, f_grid, degree);
                  }) {
      return spec.get_bulk_scattering_properties_tro_spectral(
          atm_point, f_grid, degree);
    } else {
      throw std::runtime_error(std::format(
          "Method not implemented for TRO Spectral for species:\n{:N}", spec));
    }

    std::unreachable();
  };

  auto& scat_spec = species[0];
  auto bsp        = std::visit(visitor, scat_spec);
  for (Size ind = 1; ind < species.size(); ++ind) {
    auto& scat_spec  = species[ind];
    bsp             += std::visit(visitor, scat_spec);
  }
  return bsp;
}

BulkScatteringProperties<scattering::Format::ARO,
                         scattering::Representation::Gridded>
ArrayOfScatteringSpecies::get_bulk_scattering_properties_aro_gridded(
    const AtmPoint& atm_point,
    const Vector& f_grid,
    const Vector& za_inc_grid,
    const Vector& delta_aa_grid,
    std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const {
  if (species.size() == 0) return {std::nullopt, {}, {}};

  const auto visitor = [&](const auto& spec)
      -> BulkScatteringProperties<scattering::Format::ARO,
                                  scattering::Representation::Gridded> {
    if constexpr (requires {
                    spec.get_bulk_scattering_properties_aro_gridded(
                        atm_point,
                        f_grid,
                        za_inc_grid,
                        delta_aa_grid,
                        za_scat_grid);
                  }) {
      return spec.get_bulk_scattering_properties_aro_gridded(
          atm_point, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);
    } else {
      throw std::runtime_error(std::format(
          "Method not implemented for ARO Gridded for species:\n{:N}", spec));
    }

    std::unreachable();
  };

  auto& scat_spec = species[0];
  auto bsp        = std::visit(visitor, scat_spec);
  for (Size ind = 1; ind < species.size(); ++ind) {
    auto& scat_spec  = species[ind];
    bsp             += std::visit(visitor, scat_spec);
  }
  return bsp;
}

BulkScatteringProperties<scattering::Format::ARO,
                         scattering::Representation::Spectral>
ArrayOfScatteringSpecies::get_bulk_scattering_properties_aro_spectral(
    const AtmPoint& atm_point,
    const Vector& f_grid,
    const Vector& za_inc_grid,
    Index degree,
    Index order) const {
  if (species.size() == 0) return {std::nullopt, {}, {}};

  const auto visitor = [&](const auto& spec)
      -> BulkScatteringProperties<scattering::Format::ARO,
                                  scattering::Representation::Spectral> {
    if constexpr (requires {
                    spec.get_bulk_scattering_properties_aro_spectral(
                        atm_point, f_grid, za_inc_grid, degree, order);
                  }) {
      return spec.get_bulk_scattering_properties_aro_spectral(
          atm_point, f_grid, za_inc_grid, degree, order);
    } else {
      throw std::runtime_error(std::format(
          "Method not implemented for ARO Spectral for species:\n{:N}", spec));
    }

    std::unreachable();
  };

  auto& scat_spec = species[0];
  auto bsp        = std::visit(visitor, scat_spec);
  for (Size ind = 1; ind < species.size(); ++ind) {
    auto& scat_spec  = species[ind];
    bsp             += std::visit(visitor, scat_spec);
  }
  return bsp;
}
