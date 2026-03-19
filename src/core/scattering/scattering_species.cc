#include "scattering_species.h"

#include <arts_conversions.h>
#include <stdexcept>
#include <utility>

template <typename T>
concept GriddedBulkPropertiesProviderTRO =
    requires(const T& t,
             const AtmPoint& atm_point,
             const Vector& f_grid,
             std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) {
      {
        t.get_bulk_scattering_properties_tro_gridded(
            atm_point, f_grid, za_scat_grid)
      } -> std::same_as<
          BulkScatteringProperties<scattering::Format::TRO,
                                   scattering::Representation::Gridded>>;
    };

template <typename T>
concept SpectralBulkPropertiesProviderTRO = requires(
    const T& t, const AtmPoint& atm_point, const Vector& f_grid, Index degree) {
  {
    t.get_bulk_scattering_properties_tro_spectral(atm_point, f_grid, degree)
  } -> std::same_as<ScatteringTroSpectralVector>;
};

template <typename T>
concept GriddedBulkPropertiesProviderARO =
    requires(const T& t,
             const AtmPoint& atm_point,
             const Vector& f_grid,
             const Vector& za_inc_grid,
             const Vector& delta_aa_grid,
             std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) {
      {
        t.get_bulk_scattering_properties_aro_gridded(
            atm_point, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid)
      } -> std::same_as<
          BulkScatteringProperties<scattering::Format::ARO,
                                   scattering::Representation::Gridded>>;
    };

template <typename T>
concept SpectralBulkPropertiesProviderARO = requires(const T& t,
                                                     const AtmPoint& atm_point,
                                                     const Vector& f_grid,
                                                     const Vector& za_inc_grid,
                                                     Index degree,
                                                     Index order) {
  {
    t.get_bulk_scattering_properties_aro_spectral(
        atm_point, f_grid, za_inc_grid, degree, order)
  } -> std::same_as<
      BulkScatteringProperties<scattering::Format::ARO,
                               scattering::Representation::Spectral>>;
};

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
  if (species.size() == 0) {
    return {.phase_matrix      = std::nullopt,
            .extinction_matrix = {},
            .absorption_vector = {}};
  }

  const auto visitor = [&]<typename T>(const T& spec)
      -> BulkScatteringProperties<scattering::Format::TRO,
                                  scattering::Representation::Gridded> {
    if constexpr (GriddedBulkPropertiesProviderTRO<T>) {
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

MuelmatVector ArrayOfScatteringSpecies::get_phase_matrix_at_angle(
    const AtmPoint& atm_point,
    const Vector& f_grid,
    const Vector2& los_in,
    const Vector2& los_out) const {
  const Size nf = f_grid.size();
  MuelmatVector Z(nf, rtepack::muelmat{0.0});

  if (species.empty()) return Z;

  const Numeric cos_theta = rtepack::cos_scat_angle(los_in, los_out);
  const Numeric scat_angle_deg =
      Conversion::rad2deg(std::acos(cos_theta));

  auto za_grid = std::make_shared<scattering::ZenithAngleGrid>(
      scattering::IrregularZenithAngleGrid(Vector{scat_angle_deg}));

  auto tro = get_bulk_scattering_properties_tro_gridded(
      atm_point, f_grid, za_grid);

  if (tro.phase_matrix.has_value()) {
    const auto& pm = tro.phase_matrix.value();
    for (Size iv = 0; iv < nf; iv++) {
      const auto v = pm[0, iv, 0, joker];
      Z[iv] = rtepack::muelmat(
          scattering::detail::f11(v), scattering::detail::f12(v), 0.0, 0.0,
          scattering::detail::f12(v), scattering::detail::f22(v), 0.0, 0.0,
          0.0, 0.0, scattering::detail::f33(v), scattering::detail::f34(v),
          0.0, 0.0, scattering::detail::f34(v), scattering::detail::f33(v));
    }
  }

  return Z;
}

ScatteringTroSpectralVector
ArrayOfScatteringSpecies::get_bulk_scattering_properties_tro_spectral(
    const AtmPoint& atm_point, const Vector& f_grid, Index degree) const {
  if (species.size() == 0) {
    return {.phase_matrix      = std::nullopt,
            .extinction_matrix = {},
            .absorption_vector = {}};
  }

  const auto visitor = [&](const auto& spec) -> ScatteringTroSpectralVector {
    if constexpr (SpectralBulkPropertiesProviderTRO<decltype(spec)>) {
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
  if (species.size() == 0) {
    return {.phase_matrix      = std::nullopt,
            .extinction_matrix = {},
            .absorption_vector = {}};
  }

  const auto visitor = [&](const auto& spec)
      -> BulkScatteringProperties<scattering::Format::ARO,
                                  scattering::Representation::Gridded> {
    if constexpr (GriddedBulkPropertiesProviderARO<decltype(spec)>) {
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
    if constexpr (SpectralBulkPropertiesProviderARO<decltype(spec)>) {
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
