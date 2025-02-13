#include "interpolation.h"
#include "scattering_habit.h"


namespace scattering {

  namespace detail {

    Numeric max_relative_difference(const Vector& v1, const Vector& v2) {
      if (v1.size() != v2.size()) {
        return 1.0;
      }
      Numeric max_diff = 0.0;

      for (size_t i = 0; i < v1.size(); ++i) {
        Numeric denom = std::max(std::abs(v1[i]), std::abs(v2[i]));
        if (denom > std::numeric_limits<Numeric>::epsilon()) {  // Avoid division by zero
          Numeric rel_diff = std::abs(v1[i] - v2[i]) / denom;
          max_diff = std::max(max_diff, rel_diff);
        }
      }
      return max_diff;
    }

  }

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
BulkScatteringPropertiesTROGridded
  ScatteringHabit::get_bulk_scattering_properties_tro_gridded(
      const AtmPoint& point,
      const Vector& f_grid,
      const Numeric f_tol) const {

  auto sizes = particle_habit.get_sizes(SizeParameter::DVeq);
  Index n_particles = sizes.size();
  auto pnd = std::visit([&point, &sizes, this](const auto& psd) {return psd.evaluate(point, sizes, mass_size_rel_a, mass_size_rel_b);}, psd);

  if (!particle_habit.grids.has_value()) {
    ARTS_USER_ERROR("Particle habit must be brought on a shared grid before buld properties can be computed.")
  }

  auto grids = particle_habit.grids.value();
  GridPos interp = find_interp_weights(*grids.t_grid, point[AtmKey::t]);

  Index n_freqs = grids.f_grid->size();
  Index n_angs = grid_size(*grids.za_scat_grid);

  Tensor4 phase_matrix(n_freqs, n_angs, 4, 4);
  Tensor3 extinction_matrix(n_freqs, 4, 4);
  Matrix absorption_vector(n_freqs, 4);

  for (Index part_ind = 0; part_ind < n_particles; ++part_ind) {
    try {
      auto ssd = std::get<SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>>(particle_habit[part_ind]);
      for (Index f_ind = 0; f_ind < n_freqs; ++f_ind) {
        for (Index ang_ind = 0; ang_ind < n_angs; ++ang_ind) {
          phase_matrix[f_ind, ang_ind, 0, 0] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 0];
          phase_matrix[f_ind, ang_ind, 0, 0] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 0];
          phase_matrix[f_ind, ang_ind, 0, 1] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 1];
          phase_matrix[f_ind, ang_ind, 0, 1] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 1];
          phase_matrix[f_ind, ang_ind, 1, 0] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 1];
          phase_matrix[f_ind, ang_ind, 1, 0] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 1];
          phase_matrix[f_ind, ang_ind, 1, 1] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 2];
          phase_matrix[f_ind, ang_ind, 1, 1] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 2];
          phase_matrix[f_ind, ang_ind, 2, 2] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 3];
          phase_matrix[f_ind, ang_ind, 2, 2] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 3];
          phase_matrix[f_ind, ang_ind, 2, 3] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 4];
          phase_matrix[f_ind, ang_ind, 2, 3] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 4];
          phase_matrix[f_ind, ang_ind, 3, 2] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 4];
          phase_matrix[f_ind, ang_ind, 3, 2] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 4];
          phase_matrix[f_ind, ang_ind, 3, 3] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, ang_ind, 5];
          phase_matrix[f_ind, ang_ind, 3, 3] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, ang_ind, 5];
        }
        for (Index stokes_ind = 0; stokes_ind < 4; ++stokes_ind) {
          extinction_matrix[f_ind, stokes_ind, stokes_ind] += pnd[part_ind] * interp.fd[1] * ssd.extinction_matrix[interp.idx, f_ind, stokes_ind];
          extinction_matrix[f_ind, stokes_ind, stokes_ind] += pnd[part_ind] * interp.fd[0] * ssd.extinction_matrix[interp.idx + 1, f_ind, stokes_ind];
          absorption_vector[f_ind, stokes_ind] += pnd[part_ind] * interp.fd[1] * ssd.absorption_vector[interp.idx, f_ind, stokes_ind];
          absorption_vector[f_ind, stokes_ind] += pnd[part_ind] * interp.fd[0] * ssd.absorption_vector[interp.idx + 1, f_ind, stokes_ind];
        }
      }
    } catch (const std::bad_variant_access& e) {
      ARTS_USER_ERROR("Scattering habit must be in TRO gridded format to extract bulk scattering properties.");
    }
  }

  auto f_diff = detail::max_relative_difference(f_grid, *grids.f_grid);

  // If frequency grids are the same, return calculated bulk properties.
  if (f_diff < 1e-3) {
    return BulkScatteringPropertiesTROGridded(phase_matrix,
                                              extinction_matrix,
                                              absorption_vector);
  }

  Index n_freqs_new = f_grid.size();

  // Otherwise perform frequency interpolation.
  ArrayOfGridPos interp_weights;
  gridpos(interp_weights, *grids.f_grid, f_grid);
  Tensor4 phase_matrix_new(n_freqs_new, n_angs, 4, 4);
  Tensor3 extinction_matrix_new(n_freqs_new, 4, 4);
  Matrix absorption_vector_new(n_freqs_new, 4);

  for (Index f_ind = 0; f_ind < f_grid.size(); ++f_ind) {
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
      extinction_matrix_new[f_ind, stokes_ind, stokes_ind] += weights.fd[1] * extinction_matrix[weights.idx, stokes_ind, stokes_ind];
      extinction_matrix_new[f_ind, stokes_ind, stokes_ind] += weights.fd[0] * extinction_matrix[weights.idx + 1, stokes_ind, stokes_ind];
      absorption_vector_new[f_ind, stokes_ind] += weights.fd[1] * absorption_vector[weights.idx, stokes_ind];
      absorption_vector_new[f_ind, stokes_ind] += weights.fd[0] * absorption_vector[weights.idx + 1, stokes_ind];
    }
  }
  return BulkScatteringPropertiesTROGridded(phase_matrix_new,
                                            extinction_matrix_new,
                                            absorption_vector_new);
}


ScatteringTroSpectralVector
  ScatteringHabit::get_bulk_scattering_properties_tro_spectral(
      const AtmPoint& point,
      const Vector& f_grid,
      const Numeric f_tol) const {

  auto sizes = particle_habit.get_sizes(SizeParameter::DVeq);
  Index n_particles = sizes.size();
  auto pnd = std::visit([&point, &sizes, this](const auto& psd) {return psd.evaluate(point, sizes, mass_size_rel_a, mass_size_rel_b);}, psd);

  if (!particle_habit.grids.has_value()) {
    ARTS_USER_ERROR("Particle habit must be brought on a shared grid before buld properties can be computed.")
  }

  auto grids = particle_habit.grids.value();
  GridPos interp = find_interp_weights(*grids.t_grid, point[AtmKey::t]);

  Index n_freqs = grids.f_grid->size();
  if (particle_habit.size() == 0) {
    ARTS_USER_ERROR("Encountered empty scattering habit without particles.");
  }
  auto ssd = std::get<SingleScatteringData<Numeric, Format::TRO, Representation::Spectral>>(particle_habit[0]);
  Index n_coeffs = ssd.phase_matrix.value().extent(2);

  SpecmatMatrix phase_matrix(n_freqs, n_coeffs, Specmat(0.0));
  PropmatVector extinction_matrix(n_freqs);
  StokvecVector absorption_vector(n_freqs);

  for (Index part_ind = 0; part_ind < n_particles; ++part_ind) {
    try {
      auto ssd = std::get<SingleScatteringData<Numeric, Format::TRO, Representation::Spectral>>(particle_habit[part_ind]);
      for (Index f_ind = 0; f_ind < n_freqs; ++f_ind) {
        for (Index coeff_ind = 0; coeff_ind < n_coeffs; ++coeff_ind) {
          phase_matrix[f_ind, coeff_ind][0, 0] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 0];
          phase_matrix[f_ind, coeff_ind][0, 0] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 0];
          phase_matrix[f_ind, coeff_ind][0, 1] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 1];
          phase_matrix[f_ind, coeff_ind][0, 1] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 1];
          phase_matrix[f_ind, coeff_ind][1, 0] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 1];
          phase_matrix[f_ind, coeff_ind][1, 0] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 1];
          phase_matrix[f_ind, coeff_ind][1, 1] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 2];
          phase_matrix[f_ind, coeff_ind][1, 1] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 2];
          phase_matrix[f_ind, coeff_ind][2, 2] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 3];
          phase_matrix[f_ind, coeff_ind][2, 2] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 3];
          phase_matrix[f_ind, coeff_ind][2, 3] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 4];
          phase_matrix[f_ind, coeff_ind][2, 3] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 4];
          phase_matrix[f_ind, coeff_ind][3, 2] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 4];
          phase_matrix[f_ind, coeff_ind][3, 2] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 4];
          phase_matrix[f_ind, coeff_ind][3, 3] += pnd[part_ind] * interp.fd[1] * ssd.phase_matrix.value()[interp.idx, f_ind, coeff_ind, 5];
          phase_matrix[f_ind, coeff_ind][3, 3] += pnd[part_ind] * interp.fd[0] * ssd.phase_matrix.value()[interp.idx + 1, f_ind, coeff_ind, 5];
        }

          extinction_matrix[f_ind].A() += pnd[part_ind] * interp.fd[1] * ssd.extinction_matrix[interp.idx, f_ind, 0];
          extinction_matrix[f_ind].A() += pnd[part_ind] * interp.fd[0] * ssd.extinction_matrix[interp.idx + 1, f_ind, 0];
          extinction_matrix[f_ind].B() += pnd[part_ind] * interp.fd[1] * ssd.extinction_matrix[interp.idx, f_ind, 1];
          extinction_matrix[f_ind].B() += pnd[part_ind] * interp.fd[0] * ssd.extinction_matrix[interp.idx + 1, f_ind, 1];
          extinction_matrix[f_ind].C() += pnd[part_ind] * interp.fd[1] * ssd.extinction_matrix[interp.idx, f_ind, 2];
          extinction_matrix[f_ind].C() += pnd[part_ind] * interp.fd[0] * ssd.extinction_matrix[interp.idx + 1, f_ind, 2];
          extinction_matrix[f_ind].D() += pnd[part_ind] * interp.fd[1] * ssd.extinction_matrix[interp.idx, f_ind, 3];
          extinction_matrix[f_ind].D() += pnd[part_ind] * interp.fd[0] * ssd.extinction_matrix[interp.idx + 1, f_ind, 3];

        for (Index stokes_ind = 0; stokes_ind < 4; ++stokes_ind) {
          absorption_vector[f_ind][stokes_ind] += pnd[part_ind] * interp.fd[1] * ssd.absorption_vector[interp.idx, f_ind, stokes_ind];
          absorption_vector[f_ind][stokes_ind] += pnd[part_ind] * interp.fd[0] * ssd.absorption_vector[interp.idx + 1, f_ind, stokes_ind];
        }
      }
    } catch (const std::bad_variant_access& e) {
      ARTS_USER_ERROR("Scattering habit must be in TRO gridded format to extract bulk scattering properties.");
    }
  }

  auto f_diff = detail::max_relative_difference(f_grid, *grids.f_grid);

  // If frequency grids are the same, return calculated bulk properties.
  if (f_diff < 1e-3) {
    return ScatteringTroSpectralVector(phase_matrix,
                                       extinction_matrix,
                                       absorption_vector);
  }

  Index n_freqs_new = f_grid.size();

  // Otherwise perform frequency interpolation.
  ArrayOfGridPos interp_weights{f_grid.size()};
  gridpos(interp_weights, *grids.f_grid, f_grid);
  SpecmatMatrix phase_matrix_new(n_freqs_new, n_coeffs, Specmat(0.0));
  PropmatVector extinction_matrix_new(n_freqs_new);
  StokvecVector absorption_vector_new(n_freqs_new);

  for (Index f_ind = 0; f_ind < f_grid.size(); ++f_ind) {
    auto weights = interp_weights[f_ind];
    for (Index za_scat_ind = 0; za_scat_ind < n_coeffs; ++za_scat_ind) {
        phase_matrix_new[f_ind, za_scat_ind] += weights.fd[1] * phase_matrix[weights.idx, za_scat_ind];
        phase_matrix_new[f_ind, za_scat_ind] += weights.fd[0] * phase_matrix[weights.idx + 1, za_scat_ind];
    }

    for (Index stokes_ind = 0; stokes_ind < 4; ++stokes_ind) {
      extinction_matrix_new[f_ind] += weights.fd[1] * extinction_matrix[weights.idx];
      extinction_matrix_new[f_ind] += weights.fd[0] * extinction_matrix[weights.idx + 1];
      absorption_vector_new[f_ind] += weights.fd[1] * absorption_vector[weights.idx];
      absorption_vector_new[f_ind] += weights.fd[0] * absorption_vector[weights.idx + 1];
    }
  }
  return ScatteringTroSpectralVector(phase_matrix_new,
                                     extinction_matrix_new,
                                     absorption_vector_new);
}



}
