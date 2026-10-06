#include "scattering_habit.h"

#include "configtypes.h"
#include "interpolation.h"

namespace scattering {

ScatteringHabit::ScatteringHabit(const ParticleHabit& particle_habit_,
                                 const PSD&           psd_,
                                 Numeric              mass_size_rel_a_,
                                 Numeric              mass_size_rel_b_)
    : particle_habit(particle_habit_), mass_size_rel_a(mass_size_rel_a_), mass_size_rel_b(mass_size_rel_b_), psd(psd_) {
  if ((mass_size_rel_a < 0.0) || (mass_size_rel_b < 0.0)) {
    auto size_param = std::visit([](const auto& psd) { return psd.get_size_parameter(); }, psd);
    auto [sizes, mass_size_rel_a_, mass_size_rel_b_] = particle_habit.get_size_mass_info(size_param);
    mass_size_rel_a                                  = mass_size_rel_a_;
    mass_size_rel_b                                  = mass_size_rel_b_;
  }
}

namespace detail {

Numeric max_relative_difference(const Vector& v1, const Vector& v2) {
  if (v1.size() != v2.size()) { return 1.0; }
  Numeric max_diff = 0.0;

  for (size_t i = 0; i < v1.size(); ++i) {
    Numeric denom = std::max(std::abs(v1[i]), std::abs(v2[i]));
    if (denom > std::numeric_limits<Numeric>::epsilon()) {  // Avoid division by zero
      Numeric rel_diff = std::abs(v1[i] - v2[i]) / denom;
      max_diff         = std::max(max_diff, rel_diff);
    }
  }
  return max_diff;
}

}  // namespace detail

////////////////////////////////////////////////////////////////////////////////
// Bulk properties TRO gridded
////////////////////////////////////////////////////////////////////////////////

/** Calculate bulk scattering properties for an amtospheric point
 *
 * This funciont evaluates the PSD at the given point, interpolates the scattering data along the temperature
 * dimension, and sum up the particle scattering data to calculate the bulk scattering properties.
 *
 * @param point: The AtmPoint for which to calculate the bulk properties.
 * @param f_grid: The frequencies for which to calculate the bulk scattering properties.
 * @param f_tol: Maximum relative difference between frequency vectors to trigger re-interpolation of
 *     along the frequency grid.
 * @return A struct containing the calculated bulk-scattering properties.
 */
BulkScatteringPropertiesTROGridded ScatteringHabit::get_bulk_scattering_properties_tro_gridded(const AtmPoint& point,
                                                                                               const Vector&   f_grid,
                                                                                               const Numeric) const {
  const auto  pnd         = number_densities(point);
  const Index n_particles = pnd.size();

  if (!particle_habit.grids.has_value()) {
    ARTS_USER_ERROR("Particle habit must be brought on a shared grid before buld properties can be computed.")
  }

  auto    grids  = particle_habit.grids.value();
  GridPos interp = find_interp_weights(*grids.t_grid, point[AtmKey::t]);

  Index n_freqs = grids.f_grid->size();
  Index n_angs  = grid_size(*grids.za_scat_grid);

  Tensor4 phase_matrix(n_freqs, n_angs, 4, 4);
  Tensor3 extinction_matrix(n_freqs, 4, 4);
  Matrix  absorption_vector(n_freqs, 4);

  for (Index part_ind = 0; part_ind < n_particles; ++part_ind) {
    try {
      auto ssd =
          std::get<SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>>(particle_habit[part_ind]);
      for (Index f_ind = 0; f_ind < n_freqs; ++f_ind) {
        // The scattering matrices at the two temperatures, weighted and expanded to the scattering-plane 4 x 4
        const std::array<std::pair<Numeric, Index>, 2> nodes{
            {{pnd[part_ind] * interp.fd[1], interp.idx}, {pnd[part_ind] * interp.fd[0], interp.idx + 1}}};
        for (Index ang_ind = 0; ang_ind < n_angs; ++ang_ind) {
          for (const auto& [w, it] : nodes) {
            const auto Z =
                (w * rtepack::compact_planar_muelmat{ssd.phase_matrix.value()[it, f_ind, ang_ind, joker]}).expand();
            for (Index i = 0; i < 4; ++i)
              for (Index j = 0; j < 4; ++j) phase_matrix[f_ind, ang_ind, i, j] += Z[i, j];
          }
        }
        // TRO extinction is K11 times the identity, and TRO absorption [a1, 0, 0, 0]
        for (const auto& [w, it] : nodes) {
          for (Index stokes_ind = 0; stokes_ind < 4; ++stokes_ind)
            extinction_matrix[f_ind, stokes_ind, stokes_ind] += w * ssd.extinction_matrix[it, f_ind, 0];
          absorption_vector[f_ind, 0] += w * ssd.absorption_vector[it, f_ind, 0];
        }
      }
    } catch (const std::bad_variant_access& e) {
      ARTS_USER_ERROR("Scattering habit must be in TRO gridded format to extract bulk scattering properties.");
    }
  }

  auto f_diff = detail::max_relative_difference(f_grid, *grids.f_grid);

  // If frequency grids are the same, return calculated bulk properties.
  if (f_diff < 1e-3) { return BulkScatteringPropertiesTROGridded(phase_matrix, extinction_matrix, absorption_vector); }

  Index n_freqs_new = f_grid.size();

  // Otherwise perform frequency interpolation.
  ArrayOfGridPos interp_weights;
  gridpos(interp_weights, *grids.f_grid, f_grid);
  Tensor4 phase_matrix_new(n_freqs_new, n_angs, 4, 4);
  Tensor3 extinction_matrix_new(n_freqs_new, 4, 4);
  Matrix  absorption_vector_new(n_freqs_new, 4);

  for (Size f_ind = 0; f_ind < f_grid.size(); ++f_ind) {
    auto weights = interp_weights[f_ind];
    for (Index za_scat_ind = 0; za_scat_ind < n_angs; ++za_scat_ind) {
      phase_matrix_new[f_ind, za_scat_ind, 0, 0] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 0, 0];
      phase_matrix_new[f_ind, za_scat_ind, 0, 0] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 0, 0];
      phase_matrix_new[f_ind, za_scat_ind, 0, 1] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 0, 1];
      phase_matrix_new[f_ind, za_scat_ind, 0, 1] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 0, 1];
      phase_matrix_new[f_ind, za_scat_ind, 1, 0] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 1, 0];
      phase_matrix_new[f_ind, za_scat_ind, 1, 0] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 1, 0];
      phase_matrix_new[f_ind, za_scat_ind, 1, 1] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 1, 1];
      phase_matrix_new[f_ind, za_scat_ind, 1, 1] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 1, 1];
      phase_matrix_new[f_ind, za_scat_ind, 2, 2] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 2, 2];
      phase_matrix_new[f_ind, za_scat_ind, 2, 2] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 2, 2];
      phase_matrix_new[f_ind, za_scat_ind, 2, 3] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 2, 3];
      phase_matrix_new[f_ind, za_scat_ind, 2, 3] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 2, 3];
      phase_matrix_new[f_ind, za_scat_ind, 3, 2] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 3, 2];
      phase_matrix_new[f_ind, za_scat_ind, 3, 2] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 3, 2];
      phase_matrix_new[f_ind, za_scat_ind, 3, 3] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind, 3, 3];
      phase_matrix_new[f_ind, za_scat_ind, 3, 3] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind, 3, 3];
    }

    for (Index stokes_ind = 0; stokes_ind < 4; ++stokes_ind) {
      extinction_matrix_new[f_ind, stokes_ind, stokes_ind] +=
          weights.fd[1] * extinction_matrix[weights.idx, stokes_ind, stokes_ind];
      extinction_matrix_new[f_ind, stokes_ind, stokes_ind] +=
          weights.fd[0] * extinction_matrix[weights.idx + 1, stokes_ind, stokes_ind];
      absorption_vector_new[f_ind, stokes_ind] += weights.fd[1] * absorption_vector[weights.idx, stokes_ind];
      absorption_vector_new[f_ind, stokes_ind] += weights.fd[0] * absorption_vector[weights.idx + 1, stokes_ind];
    }
  }
  return BulkScatteringPropertiesTROGridded(phase_matrix_new, extinction_matrix_new, absorption_vector_new);
}

BulkScatteringProperties<Format::TRO, Representation::Gridded>
ScatteringHabit::get_bulk_scattering_properties_tro_gridded(const AtmPoint&                  point,
                                                            const Vector&                    f_grid,
                                                            std::shared_ptr<ZenithAngleGrid> za_scat_grid) const {
  const auto pnd = number_densities(point);
  ARTS_USER_ERROR_IF(pnd.size() != static_cast<Size>(particle_habit.size()), "PSD and particle-habit sizes differ.")

  auto grids = ScatteringDataGrids(
      std::make_shared<Vector>(Vector{point.temperature}), std::make_shared<Vector>(f_grid), std::move(za_scat_grid));
  using Bulk = BulkScatteringProperties<Format::TRO, Representation::Gridded>;
  std::optional<Bulk> result;
  for (Index i = 0; i < particle_habit.size(); ++i) {
    auto data = std::visit([&grids](const auto& ssd) { return ssd_to_tro_gridded(grids, ssd); }, particle_habit[i]);
    if (pnd[i] == 0.0) {
      if (data.phase_matrix) std::fill_n(data.phase_matrix->data_handle(), data.phase_matrix->size(), Numeric{0.0});
      std::fill_n(data.extinction_matrix.data_handle(), data.extinction_matrix.size(), Numeric{0.0});
      std::fill_n(data.absorption_vector.data_handle(), data.absorption_vector.size(), Numeric{0.0});
    } else {
      if (data.phase_matrix) *data.phase_matrix *= pnd[i];
      data.extinction_matrix *= pnd[i];
      data.absorption_vector *= pnd[i];
    }
    Bulk bulk{std::move(data.phase_matrix), std::move(data.extinction_matrix), std::move(data.absorption_vector)};
    if (result) {
      *result += bulk;
    } else {
      result = std::move(bulk);
    }
  }
  ARTS_USER_ERROR_IF(not result, "Cannot calculate bulk properties for an empty particle habit.")
  return std::move(*result);
}

BulkScatteringProperties<Format::TRO, Representation::Gridded>
ScatteringHabit::get_bulk_scattering_properties_tro_gridded_derivative(const AtmPoint&                  point,
                                                                       const Vector&                    f_grid,
                                                                       std::shared_ptr<ZenithAngleGrid> za_scat_grid,
                                                                       const AtmKeyVal&                 target) const {
  const auto pnd = number_densities_with_derivatives(point);
  ARTS_USER_ERROR_IF(pnd.values.size() != static_cast<Size>(particle_habit.size()),
                     "PSD and particle-habit sizes differ.")

  const auto*   property   = std::get_if<ScatteringSpeciesProperty>(&target);
  const Vector* derivative = nullptr;
  if (property) {
    const auto iter = pnd.derivatives.find(*property);
    if (iter != pnd.derivatives.end()) derivative = &iter->second;
  }

  auto grids = ScatteringDataGrids(
      std::make_shared<Vector>(Vector{point.temperature}), std::make_shared<Vector>(f_grid), za_scat_grid);
  using Bulk = BulkScatteringProperties<Format::TRO, Representation::Gridded>;
  std::optional<Bulk> result;
  for (Index i = 0; i < particle_habit.size(); ++i) {
    auto data = std::visit([&grids](const auto& ssd) { return ssd_to_tro_gridded(grids, ssd); }, particle_habit[i]);
    const Numeric weight = derivative ? (*derivative)[i] : 0.0;
    if (weight == 0.0) {
      if (data.phase_matrix) std::fill_n(data.phase_matrix->data_handle(), data.phase_matrix->size(), Numeric{0.0});
      std::fill_n(data.extinction_matrix.data_handle(), data.extinction_matrix.size(), Numeric{0.0});
      std::fill_n(data.absorption_vector.data_handle(), data.absorption_vector.size(), Numeric{0.0});
    } else {
      if (data.phase_matrix) *data.phase_matrix *= weight;
      data.extinction_matrix *= weight;
      data.absorption_vector *= weight;
    }
    Bulk bulk{std::move(data.phase_matrix), std::move(data.extinction_matrix), std::move(data.absorption_vector)};

    if (target == AtmKeyVal{AtmKey::t} and pnd.values[i] != 0.0) {
      Bulk temperature_derivative = std::visit(
          [&](const auto& ssd) -> Bulk {
            if constexpr (std::remove_cvref_t<decltype(ssd)>::get_format() == Format::ARO) {
              ARTS_USER_ERROR("TRO radar derivatives cannot use ARO particle data")
            } else {
              const auto gridded = [&]() {
                if constexpr (std::remove_cvref_t<decltype(ssd)>::get_representation() == Representation::Gridded)
                  return ssd;
                else
                  return ssd.to_gridded();
              }();
              ARTS_USER_ERROR_IF(not gridded.phase_matrix, "Temperature derivatives require particle phase matrices")
              const Vector& temperatures = *gridded.phase_matrix->get_t_grid();
              ARTS_USER_ERROR_IF(temperatures.size() < 2,
                                 "Temperature-dependent particle derivatives require at least two temperatures")

              auto value_at = [&](Numeric temperature) {
                ScatteringDataGrids local_grids(
                    std::make_shared<Vector>(Vector{temperature}), std::make_shared<Vector>(f_grid), za_scat_grid);
                auto value = ssd_to_tro_gridded(local_grids, ssd);
                return Bulk{std::move(value.phase_matrix),
                            std::move(value.extinction_matrix),
                            std::move(value.absorption_vector)};
              };
              auto       upper   = std::ranges::lower_bound(temperatures, point.temperature);
              const Size located = static_cast<Size>(upper - temperatures.begin());
              Size       right   = std::clamp<Size>(located, 1, temperatures.size() - 1);
              Size       left    = right - 1;

              if (upper != temperatures.end() and *upper == point.temperature and located == right and right > 0 and
                  right + 1 < temperatures.size()) {
                auto below   = value_at(temperatures[right - 1]);
                auto middle  = value_at(temperatures[right]);
                auto above   = value_at(temperatures[right + 1]);
                below       *= -0.5 / (temperatures[right] - temperatures[right - 1]);
                middle      *= 0.5 / (temperatures[right] - temperatures[right - 1]) -
                               0.5 / (temperatures[right + 1] - temperatures[right]);
                above       *= 0.5 / (temperatures[right + 1] - temperatures[right]);
                below       += middle;
                below       += above;
                return below;
              }

              auto          below            = value_at(temperatures[left]);
              auto          above            = value_at(temperatures[right]);
              const Numeric inverse_spacing  = 1.0 / (temperatures[right] - temperatures[left]);
              below                         *= -inverse_spacing;
              above                         *= inverse_spacing;
              below                         += above;
              return below;
            }
            std::unreachable();
          },
          particle_habit[i]);
      temperature_derivative *= pnd.values[i];
      bulk                   += temperature_derivative;
    }
    if (result)
      *result += bulk;
    else
      result = std::move(bulk);
  }
  ARTS_USER_ERROR_IF(not result, "Cannot calculate bulk properties for an empty particle habit.")
  return std::move(*result);
}

BulkScatteringProperties<Format::ARO, Representation::Gridded>
ScatteringHabit::get_bulk_scattering_properties_aro_gridded(const AtmPoint&                  point,
                                                            const Vector&                    f_grid,
                                                            const Vector&                    za_inc_grid,
                                                            const Vector&                    delta_aa_grid,
                                                            std::shared_ptr<ZenithAngleGrid> za_scat_grid) const {
  const auto pnd = number_densities(point);
  ARTS_USER_ERROR_IF(pnd.size() != static_cast<Size>(particle_habit.size()), "PSD and particle-habit sizes differ.")

  auto grids = ScatteringDataGrids(std::make_shared<Vector>(Vector{point.temperature}),
                                   std::make_shared<Vector>(f_grid),
                                   std::make_shared<Vector>(za_inc_grid),
                                   std::make_shared<Vector>(delta_aa_grid),
                                   std::move(za_scat_grid));
  using SSD  = SingleScatteringData<Numeric, Format::ARO, Representation::Gridded>;
  std::optional<BulkScatteringProperties<Format::ARO, Representation::Gridded>> result;
  for (Index i = 0; i < particle_habit.size(); ++i) {
    SSD data = std::visit([&grids](const auto& ssd) { return ssd_to_aro_gridded(grids, ssd); }, particle_habit[i]);
    if (pnd[i] == 0.0) {
      // A laboratory-frame conversion can contain an indeterminate 0*NaN at
      // angular-coordinate poles.  An absent particle population is exactly
      // zero regardless of its single-particle coordinate representation.
      if (data.phase_matrix) std::fill_n(data.phase_matrix->data_handle(), data.phase_matrix->size(), Numeric{0.0});
      std::fill_n(data.extinction_matrix.data_handle(), data.extinction_matrix.size(), Numeric{0.0});
      std::fill_n(data.absorption_vector.data_handle(), data.absorption_vector.size(), Numeric{0.0});
    } else {
      if (data.phase_matrix) *data.phase_matrix *= pnd[i];
      data.extinction_matrix *= pnd[i];
      data.absorption_vector *= pnd[i];
    }
    BulkScatteringProperties<Format::ARO, Representation::Gridded> bulk{
        std::move(data.phase_matrix), std::move(data.extinction_matrix), std::move(data.absorption_vector)};
    if (result) {
      *result += bulk;
    } else {
      result = std::move(bulk);
    }
  }
  ARTS_USER_ERROR_IF(not result, "Cannot calculate bulk properties for an empty particle habit.")
  return std::move(*result);
}

BulkScatteringProperties<Format::ARO, Representation::Gridded>
ScatteringHabit::get_bulk_scattering_properties_aro_gridded_derivative(const AtmPoint&                  point,
                                                                       const Vector&                    f_grid,
                                                                       const Vector&                    za_inc_grid,
                                                                       const Vector&                    delta_aa_grid,
                                                                       std::shared_ptr<ZenithAngleGrid> za_scat_grid,
                                                                       const AtmKeyVal&                 target) const {
  const auto pnd = number_densities_with_derivatives(point);
  ARTS_USER_ERROR_IF(pnd.values.size() != static_cast<Size>(particle_habit.size()),
                     "PSD and particle-habit sizes differ.")

  const auto*   property   = std::get_if<ScatteringSpeciesProperty>(&target);
  const Vector* derivative = nullptr;
  if (property) {
    const auto iter = pnd.derivatives.find(*property);
    if (iter != pnd.derivatives.end()) derivative = &iter->second;
  }

  auto grids = ScatteringDataGrids(std::make_shared<Vector>(Vector{point.temperature}),
                                   std::make_shared<Vector>(f_grid),
                                   std::make_shared<Vector>(za_inc_grid),
                                   std::make_shared<Vector>(delta_aa_grid),
                                   za_scat_grid);
  using Bulk = BulkScatteringProperties<Format::ARO, Representation::Gridded>;
  std::optional<Bulk> result;
  for (Index i = 0; i < particle_habit.size(); ++i) {
    auto data = std::visit([&grids](const auto& ssd) { return ssd_to_aro_gridded(grids, ssd); }, particle_habit[i]);
    const Numeric weight = derivative ? (*derivative)[i] : 0.0;
    if (weight == 0.0) {
      if (data.phase_matrix) std::fill_n(data.phase_matrix->data_handle(), data.phase_matrix->size(), Numeric{0.0});
      std::fill_n(data.extinction_matrix.data_handle(), data.extinction_matrix.size(), Numeric{0.0});
      std::fill_n(data.absorption_vector.data_handle(), data.absorption_vector.size(), Numeric{0.0});
    } else {
      if (data.phase_matrix) *data.phase_matrix *= weight;
      data.extinction_matrix *= weight;
      data.absorption_vector *= weight;
    }
    Bulk bulk{std::move(data.phase_matrix), std::move(data.extinction_matrix), std::move(data.absorption_vector)};

    if (target == AtmKeyVal{AtmKey::t} and pnd.values[i] != 0.0) {
      Bulk temperature_derivative = std::visit(
          [&](const auto& ssd) -> Bulk {
            const auto gridded = [&]() {
              if constexpr (std::remove_cvref_t<decltype(ssd)>::get_representation() == Representation::Gridded)
                return ssd;
              else
                return ssd.to_gridded();
            }();
            ARTS_USER_ERROR_IF(not gridded.phase_matrix, "Temperature derivatives require particle phase matrices")
            const Vector& temperatures = *gridded.phase_matrix->get_t_grid();
            ARTS_USER_ERROR_IF(temperatures.size() < 2,
                               "Temperature-dependent particle derivatives require at least two temperatures")

            auto value_at = [&](Numeric temperature) {
              ScatteringDataGrids local_grids(std::make_shared<Vector>(Vector{temperature}),
                                              std::make_shared<Vector>(f_grid),
                                              std::make_shared<Vector>(za_inc_grid),
                                              std::make_shared<Vector>(delta_aa_grid),
                                              za_scat_grid);
              auto                value = ssd_to_aro_gridded(local_grids, ssd);
              return Bulk{std::move(value.phase_matrix),
                          std::move(value.extinction_matrix),
                          std::move(value.absorption_vector)};
            };
            auto       upper   = std::ranges::lower_bound(temperatures, point.temperature);
            const Size located = static_cast<Size>(upper - temperatures.begin());
            const Size right   = std::clamp<Size>(located, 1, temperatures.size() - 1);
            const Size left    = right - 1;

            if (upper != temperatures.end() and *upper == point.temperature and located == right and right > 0 and
                right + 1 < temperatures.size()) {
              auto below   = value_at(temperatures[right - 1]);
              auto middle  = value_at(temperatures[right]);
              auto above   = value_at(temperatures[right + 1]);
              below       *= -0.5 / (temperatures[right] - temperatures[right - 1]);
              middle      *= 0.5 / (temperatures[right] - temperatures[right - 1]) -
                             0.5 / (temperatures[right + 1] - temperatures[right]);
              above       *= 0.5 / (temperatures[right + 1] - temperatures[right]);
              below       += middle;
              below       += above;
              return below;
            }

            auto          below            = value_at(temperatures[left]);
            auto          above            = value_at(temperatures[right]);
            const Numeric inverse_spacing  = 1.0 / (temperatures[right] - temperatures[left]);
            below                         *= -inverse_spacing;
            above                         *= inverse_spacing;
            below                         += above;
            return below;
          },
          particle_habit[i]);
      temperature_derivative *= pnd.values[i];
      bulk                   += temperature_derivative;
    }
    if (result)
      *result += bulk;
    else
      result = std::move(bulk);
  }
  ARTS_USER_ERROR_IF(not result, "Cannot calculate bulk properties for an empty particle habit.")
  return std::move(*result);
}

namespace {
/** Why particle i of a habit has no Legendre series, for the error message */
std::string_view no_legendre_series(const ParticleData& data) {
  if (std::holds_alternative<SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>>(data))
    return "it holds gridded TRO data; convert the habit with ParticleHabit.to_tro_spectral (or "
           "to_tro_spectral_with_report, which also reports how well the series represents the data)";
  return "it holds ARO data, whose phase matrix depends on more than the scattering angle";
}

using TroSpectralSSD = SingleScatteringData<Numeric, Format::TRO, Representation::Spectral>;

/** The TRO Legendre series of particle i at the temperature and frequencies of grids, checked */
TroSpectralSSD tro_series(const ParticleData& data, Index i, Index degree, const ScatteringDataGrids& grids) {
  const auto* ssd = std::get_if<TroSpectralSSD>(&data);
  ARTS_USER_ERROR_IF(
      not ssd, "Particle {} of the scattering habit has no Legendre series: {}", i, no_legendre_series(data))
  ARTS_USER_ERROR_IF(not ssd->phase_matrix, "Particle {} of the scattering habit has no phase matrix", i)
  ARTS_USER_ERROR_IF(ssd->phase_matrix->get_degree() < degree,
                     "Particle {} of the scattering habit has a Legendre series to degree {}, which does not give "
                     "the coefficients to degree {}; convert its gridded data to that degree",
                     i,
                     ssd->phase_matrix->get_degree(),
                     degree)
  return ssd->regrid(grids);
}
}  // namespace

ScatteringTroSpectralVector ScatteringHabit::get_bulk_scattering_properties_tro_spectral(const AtmPoint& point,
                                                                                         const Vector&   f_grid,
                                                                                         const Index     degree) const {
  ARTS_USER_ERROR_IF(degree < 0, "The Legendre degree must be >= 0, got {}", degree)
  ARTS_USER_ERROR_IF(particle_habit.size() == 0, "Cannot calculate bulk properties for an empty particle habit.")
  const auto pnd = number_densities(point);
  const auto grids =
      ScatteringDataGrids(std::make_shared<Vector>(Vector{point.temperature}), std::make_shared<Vector>(f_grid));

  const Index                 nf = f_grid.size();
  ScatteringTroSpectralVector out{.phase_matrix      = SpecmatMatrix(nf, degree + 1, Specmat{0.0}),
                                  .extinction_matrix = PropmatVector(nf, Propmat{}),
                                  .absorption_vector = StokvecVector(nf, Stokvec{})};
  for (Index i = 0; i < particle_habit.size(); ++i) {
    const auto ssd = tro_series(particle_habit[i], i, degree, grids);
    if (pnd[i] == 0.0) continue;
    for (Index iv = 0; iv < nf; ++iv) {
      for (Index l = 0; l <= degree; ++l) {
        // The scattering-plane Mueller matrix of each (real) coefficient
        const auto coeffs = (*ssd.phase_matrix)[0, iv, l, joker];
        const auto F =
            pnd[i] * rtepack::compact_planar_muelmat{coeffs | std::views::transform([](const Complex& x) { return x.real(); })};
        (*out.phase_matrix)[iv, l] += Specmat{F.expand().data};
      }
      out.extinction_matrix[iv].A() += pnd[i] * ssd.extinction_matrix[0, iv, 0];
      out.absorption_vector[iv][0]  += pnd[i] * ssd.absorption_vector[0, iv, 0];
    }
  }
  return out;
}

BulkScatteringProperties<Format::ARO, Representation::Fourier>
ScatteringHabit::get_bulk_scattering_properties_aro_fourier(const AtmPoint& point,
                                                             const Vector&   f_grid,
                                                             const Vector&   za_inc_grid,
                                                             const Vector&   za_scat_grid,
                                                             Index           max_mode) const {
  ARTS_USER_ERROR_IF(particle_habit.size() == 0, "Cannot calculate bulk properties for an empty particle habit.")
  const auto pnd     = number_densities(point);
  auto       t_grid  = std::make_shared<Vector>(Vector{point.temperature});
  auto       f_ptr   = std::make_shared<Vector>(f_grid);
  auto       za_inc  = std::make_shared<const Vector>(za_inc_grid);
  auto       za_scat = std::make_shared<const ZenithAngleGrid>(IrregularZenithAngleGrid(za_scat_grid));
  const auto grids   = ScatteringDataGrids(t_grid, f_ptr, za_inc, nullptr, za_scat);
  const auto tf      = ScatteringDataGrids(t_grid, f_ptr);

  using Bulk = BulkScatteringProperties<Format::ARO, Representation::Fourier>;
  std::optional<Bulk> result;
  const auto          add = [&result](Bulk&& b) {
    if (result)
      *result += b;
    else
      result = std::move(b);
  };

  // The Legendre series of TRO particles add up before the one conversion to Fourier modes, which is linear
  std::optional<TroSpectralSSD> tro;
  for (Index i = 0; i < particle_habit.size(); ++i) {
    if (const auto* ssd = std::get_if<TroSpectralSSD>(&particle_habit[i])) {
      ARTS_USER_ERROR_IF(not ssd->phase_matrix, "Particle {} of the scattering habit has no phase matrix", i)
      if (pnd[i] == 0.0) continue;
      auto local               = ssd->regrid(tf);
      *local.phase_matrix     *= pnd[i];
      local.extinction_matrix *= pnd[i];
      local.absorption_vector *= pnd[i];
      if (not tro) {
        tro = std::move(local);
        continue;
      }
      // Series of different degrees add as the longer one
      if (local.phase_matrix->get_degree() > tro->phase_matrix->get_degree()) std::swap(*tro, local);
      for (Index iv = 0; iv < f_ptr->size(); ++iv)
        for (Index l = 0; l <= local.phase_matrix->get_degree(); ++l)
          for (Index e = 0; e < 6; ++e) (*tro->phase_matrix)[0, iv, l, e] += (*local.phase_matrix)[0, iv, l, e];
      tro->extinction_matrix += local.extinction_matrix;
      tro->absorption_vector += local.absorption_vector;
    } else {
      auto data = std::visit([&](const auto& s) { return ssd_to_aro_fourier(grids, max_mode, s); }, particle_habit[i]);
      if (pnd[i] == 0.0) continue;
      ARTS_USER_ERROR_IF(not data.phase_matrix, "Particle {} of the scattering habit has no phase matrix", i)
      Bulk bulk{.phase_matrix      = std::move(data.phase_matrix),
                .extinction_matrix = std::move(data.extinction_matrix),
                .absorption_vector = std::move(data.absorption_vector)};
      bulk *= pnd[i];
      add(std::move(bulk));
    }
  }
  if (tro) {
    auto data = tro->to_lab_frame_fourier_modes(grids, max_mode);
    add(Bulk{std::move(data.phase_matrix), std::move(data.extinction_matrix), std::move(data.absorption_vector)});
  }
  if (not result) {
    // Every particle is absent: zero optics on the requested grids
    Bulk zero{PhaseMatrixData<Numeric, Format::ARO, Representation::Fourier>(t_grid, f_ptr, za_inc, za_scat, max_mode),
              ExtinctionMatrixData<Numeric, Format::ARO, Representation::Fourier>(t_grid, f_ptr, za_inc),
              AbsorptionVectorData<Numeric, Format::ARO, Representation::Fourier>(t_grid, f_ptr, za_inc)};
    return zero;
  }
  return std::move(*result);
}

Vector ScatteringHabit::number_densities(const AtmPoint& point) const {
  const auto sizes = particle_habit.get_sizes(std::visit([](const auto& p) { return p.get_size_parameter(); }, psd));
  auto       pnd   = std::visit(
      [&](const auto& p) { return scattering::number_densities(p, point, sizes, mass_size_rel_a, mass_size_rel_b); },
      psd);
  ARTS_USER_ERROR_IF(pnd.size() != static_cast<Size>(particle_habit.size()), "PSD and particle-habit sizes differ.")
  return pnd;
}

PSDData ScatteringHabit::number_densities_with_derivatives(const AtmPoint& point) const {
  const auto sizes = particle_habit.get_sizes(std::visit([](const auto& p) { return p.get_size_parameter(); }, psd));
  auto       pnd   = std::visit(
      [&](const auto& p) {
        return scattering::number_densities_with_derivatives(p, point, sizes, mass_size_rel_a, mass_size_rel_b);
      },
      psd);
  ARTS_USER_ERROR_IF(pnd.values.size() != static_cast<Size>(particle_habit.size()),
                     "PSD and particle-habit sizes differ.")
  return pnd;
}

}  // namespace scattering

void xml_io_stream<scattering::ScatteringHabit>::write(std::ostream&,
                                                       const scattering::ScatteringHabit&,
                                                       bofstream*,
                                                       std::string_view) {
  throw std::runtime_error("private data not readable");
}

void xml_io_stream<scattering::ScatteringHabit>::read(std::istream&, scattering::ScatteringHabit&, bifstream*) {
  throw std::runtime_error("private data not writeable");
}
