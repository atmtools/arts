#include "scattering_species.h"

#include <stdexcept>
#include <utility>

void ArrayOfScatteringSpecies::add(const scattering::Species& species_) { species.push_back(species_); }

void ArrayOfScatteringSpecies::prepare_scattering_data(scattering::ScatteringDataSpec) {}

template <scattering::medium_point Point, scattering::spectral_medium_point<Point> SpectralPoint>
BulkScatteringProperties<scattering::Format::TRO, scattering::Representation::Gridded>
ArrayOfScatteringSpecies::tro_gridded(const Point&                              point,
    const SpectralPoint& spectral_point,
                                      const Vector&                                f_grid,
                                      std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const {
  if (species.size() == 0) return {std::nullopt, {}, {}};

  const auto visitor =
      [&](const auto& spec) -> BulkScatteringProperties<scattering::Format::TRO, scattering::Representation::Gridded> {
    if constexpr (requires { spec.tro_gridded(point, spectral_point, f_grid, za_scat_grid); }) {
      return spec.tro_gridded(point, spectral_point, f_grid, za_scat_grid);
    } else {
      throw std::runtime_error(arts::format("Method not implemented for TRO Gridded for species:\n{:N}", spec));
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

template <scattering::medium_point Point, scattering::spectral_medium_point<Point> SpectralPoint>
BulkScatteringProperties<scattering::Format::TRO, scattering::Representation::Gridded>
ArrayOfScatteringSpecies::tro_gridded_derivative(const Point&                              point,
    const SpectralPoint& spectral_point,
                                                 const Vector&                                f_grid,
                                                 std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid,
                                                 const scattering::medium_key<Point>&                             target) const {
  if (species.empty()) arts::user_error("Cannot differentiate an empty scattering-species array");
  const auto visitor =
      [&](const auto& spec) -> BulkScatteringProperties<scattering::Format::TRO, scattering::Representation::Gridded> {
    if constexpr (requires { spec.tro_gridded_derivative(point, spectral_point, f_grid, za_scat_grid, target); }) {
      return spec.tro_gridded_derivative(point, spectral_point, f_grid, za_scat_grid, target);
    } else {
      throw std::runtime_error(arts::format("TRO-gridded derivatives are not implemented for species:\n{:N}", spec));
    }
    std::unreachable();
  };
  auto out = std::visit(visitor, species.front());
  for (Size i = 1; i < species.size(); ++i) out += std::visit(visitor, species[i]);
  return out;
}

template <scattering::medium_point Point, scattering::spectral_medium_point<Point> SpectralPoint>
ScatteringTroSpectralVector ArrayOfScatteringSpecies::tro_spectral(const Point& point,
    const SpectralPoint& spectral_point,
                                                                   const Vector&   f_grid,
                                                                   Index           degree) const {
  if (species.size() == 0) return {.phase_matrix = std::nullopt, .extinction_matrix = {}, .absorption_vector = {}};

  const auto visitor = [&](const auto& spec) -> ScatteringTroSpectralVector {
    if constexpr (requires { spec.tro_spectral(point, spectral_point, f_grid, degree); }) {
      return spec.tro_spectral(point, spectral_point, f_grid, degree);
    } else {
      throw std::runtime_error(arts::format("Method not implemented for TRO Spectral for species:\n{:N}", spec));
    }

    std::unreachable();
  };

  auto bsp = std::visit(visitor, species[0]);
  for (Size ind = 1; ind < species.size(); ++ind) bsp += std::visit(visitor, species[ind]);
  return bsp;
}

template <scattering::medium_point Point, scattering::spectral_medium_point<Point> SpectralPoint>
BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded>
ArrayOfScatteringSpecies::aro_gridded(const Point&                              point,
    const SpectralPoint& spectral_point,
                                      const Vector&                                f_grid,
                                      const Vector&                                za_inc_grid,
                                      const Vector&                                delta_aa_grid,
                                      std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const {
  if (species.size() == 0) return {.phase_matrix = std::nullopt, .extinction_matrix = {}, .absorption_vector = {}};

  const auto visitor =
      [&](const auto& spec) -> BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded> {
    if constexpr (requires { spec.aro_gridded(point, spectral_point, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid); }) {
      return spec.aro_gridded(point, spectral_point, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);
    } else {
      throw std::runtime_error(arts::format("Method not implemented for ARO Gridded for species:\n{:N}", spec));
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

template <scattering::medium_point Point, scattering::spectral_medium_point<Point> SpectralPoint>
BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded>
ArrayOfScatteringSpecies::aro_gridded_derivative(const Point&                              point,
    const SpectralPoint& spectral_point,
                                                 const Vector&                                f_grid,
                                                 const Vector&                                za_inc_grid,
                                                 const Vector&                                delta_aa_grid,
                                                 std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid,
                                                 const scattering::medium_key<Point>&                             target) const {
  if (species.empty()) arts::user_error("Cannot differentiate an empty scattering-species array");
  const auto visitor =
      [&](const auto& spec) -> BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded> {
    if constexpr (requires {
                    spec.aro_gridded_derivative(point, spectral_point, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid, target);
                  }) {
      return spec.aro_gridded_derivative(point, spectral_point, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid, target);
    } else {
      throw std::runtime_error(arts::format("ARO-gridded derivatives are not implemented for species:\n{:N}", spec));
    }
    std::unreachable();
  };
  auto out = std::visit(visitor, species.front());
  for (Size i = 1; i < species.size(); ++i) out += std::visit(visitor, species[i]);
  return out;
}

template <scattering::medium_point Point, scattering::spectral_medium_point<Point> SpectralPoint>
BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Spectral>
ArrayOfScatteringSpecies::aro_spectral(const Point&         point,
                                       const SpectralPoint& spectral_point,
                                       const Vector&        f_grid,
                                       const Vector&        za_inc_grid,
                                       Index                degree,
                                       Index                order) const {
  if (species.size() == 0) return {.phase_matrix = std::nullopt, .extinction_matrix = {}, .absorption_vector = {}};

  const auto visitor =
      [&](const auto& spec) -> BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Spectral> {
    if constexpr (requires { spec.aro_spectral(point, spectral_point, f_grid, za_inc_grid, degree, order); }) {
      return spec.aro_spectral(point, spectral_point, f_grid, za_inc_grid, degree, order);
    } else {
      throw std::runtime_error(arts::format("Method not implemented for ARO Spectral for species:\n{:N}", spec));
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

template <scattering::medium_point Point, scattering::spectral_medium_point<Point> SpectralPoint>
BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Fourier>
ArrayOfScatteringSpecies::aro_fourier(const Point& point,
    const SpectralPoint& spectral_point,
                                      const Vector&   f_grid,
                                      const Vector&   za_inc_grid,
                                      const Vector&   za_scat_grid,
                                      Index           max_mode) const {
  if (species.size() == 0) return {.phase_matrix = std::nullopt, .extinction_matrix = {}, .absorption_vector = {}};

  const auto visitor =
      [&](const auto& spec) -> BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Fourier> {
    if constexpr (requires { spec.aro_fourier(point, spectral_point, f_grid, za_inc_grid, za_scat_grid, max_mode); }) {
      return spec.aro_fourier(point, spectral_point, f_grid, za_inc_grid, za_scat_grid, max_mode);
    } else {
      throw std::runtime_error(arts::format("Method not implemented for ARO Fourier modes for species:\n{:N}", spec));
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

#define ARTS_INSTANTIATE_ARRAY_OF_SCATTERING_SPECIES(Point, SpectralPoint)                                          \
  template BulkScatteringProperties<scattering::Format::TRO, scattering::Representation::Gridded>                    \
  ArrayOfScatteringSpecies::tro_gridded(                                                                             \
      const Point&, const SpectralPoint&, const Vector&, std::shared_ptr<scattering::ZenithAngleGrid>) const;        \
  template BulkScatteringProperties<scattering::Format::TRO, scattering::Representation::Gridded>                    \
  ArrayOfScatteringSpecies::tro_gridded_derivative(const Point&,                                                     \
                                                   const SpectralPoint&,                                             \
                                                   const Vector&,                                                    \
                                                   std::shared_ptr<scattering::ZenithAngleGrid>,                     \
                                                   const scattering::medium_key<Point>&) const;                      \
  template ScatteringTroSpectralVector ArrayOfScatteringSpecies::tro_spectral(                                       \
      const Point&, const SpectralPoint&, const Vector&, Index) const;                                               \
  template BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded>                    \
  ArrayOfScatteringSpecies::aro_gridded(const Point&,                                                                \
                                        const SpectralPoint&,                                                        \
                                        const Vector&,                                                               \
                                        const Vector&,                                                               \
                                        const Vector&,                                                               \
                                        std::shared_ptr<scattering::ZenithAngleGrid>) const;                         \
  template BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded>                    \
  ArrayOfScatteringSpecies::aro_gridded_derivative(const Point&,                                                     \
                                                   const SpectralPoint&,                                             \
                                                   const Vector&,                                                    \
                                                   const Vector&,                                                    \
                                                   const Vector&,                                                    \
                                                   std::shared_ptr<scattering::ZenithAngleGrid>,                     \
                                                   const scattering::medium_key<Point>&) const;                      \
  template BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Spectral>                   \
  ArrayOfScatteringSpecies::aro_spectral(                                                                            \
      const Point&, const SpectralPoint&, const Vector&, const Vector&, Index, Index) const;                          \
  template BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Fourier>                    \
  ArrayOfScatteringSpecies::aro_fourier(                                                                             \
      const Point&, const SpectralPoint&, const Vector&, const Vector&, const Vector&, Index) const;

ARTS_INSTANTIATE_ARRAY_OF_SCATTERING_SPECIES(AtmPoint, AtmSpectralPoint)
ARTS_INSTANTIATE_ARRAY_OF_SCATTERING_SPECIES(SubsurfacePoint, SubsurfaceSpectralPoint)

#undef ARTS_INSTANTIATE_ARRAY_OF_SCATTERING_SPECIES
