#include <iostream>

#include "integration.h"
#include "particle_habit.h"
#include "scattering_habit.h"
#include "test_utils.h"

// Test
bool test_calculate_bulk_properties_tro_gridded() {

  Vector t_grid{280.0, 290.0, 300.0};
  Vector f_grid{10e9, 89e9, 183e9};
  Vector diameters{50e-6, 500e-6, 1e-3};
  scattering::IrregularZenithAngleGrid za_scat_grid = Vector{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0,
                                            90.0, 100., 110.0, 120.0, 130.0, 140.0, 150.0, 160., 170.0, 180.0};
  auto particle_habit = scattering::ParticleHabit::liquid_sphere(t_grid, f_grid, diameters, za_scat_grid);
  auto sizes = particle_habit.get_sizes(scattering::SizeParameter::VolumeEqDiameter);
  particle_habit = particle_habit.to_tro_gridded(t_grid, f_grid, za_scat_grid);
  Vector size_bins{1e-6, 100e-6, 800e-6, 1500e-6};
  Vector counts{1.0, 0.0, 0.0};
  auto psd = scattering::BinnedPSD(size_bins, counts);
  auto scattering_habit = scattering::ScatteringHabit(particle_habit, 1.0, 1.0, psd);

  ///
  /// Test temperature interpolation.
  ///
  auto point = AtmPoint{1e4, 290.0};
  auto bulk_props = scattering_habit.get_bulk_scattering_properties_tro_gridded(point, f_grid);

  auto pm = bulk_props.phase_matrix.value();
  using SSD = scattering::SingleScatteringData<Numeric, scattering::Format::TRO, scattering::Representation::Gridded>;
  auto pm_ref = std::get<SSD>(particle_habit[0]).phase_matrix.value();
  //pm_ref += 0.5 * std::get<SSD>(particle_habit[0]).phase_matrix.value()[2];
  auto em_ref = std::get<SSD>(particle_habit[0]).extinction_matrix;
  //em_ref += std::get<SSD>(particle_habit[0]).extinction_matrix[2];
  auto av_ref = std::get<SSD>(particle_habit[0]).absorption_vector;
  //av_ref += std::get<SSD>(particle_habit[0]).absorption_vector[2];
  //auto pm_ref = 0.5 * std::get<SSD>(particle_habit[0]).phase_matrix.value()[1];
  //pm_ref += 0.5 * std::get<SSD>(particle_habit[0]).phase_matrix.value()[2];
  //auto em_ref = std::get<SSD>(particle_habit[0]).extinction_matrix[1];
  //em_ref += std::get<SSD>(particle_habit[0]).extinction_matrix[2];
  //auto av_ref = std::get<SSD>(particle_habit[0]).absorption_vector[1];
  //av_ref += std::get<SSD>(particle_habit[0]).absorption_vector[2];

  /// Ensure relative errors are small.
  for (Index f_ind = 0; f_ind < f_grid.size(); ++f_ind) {
    for (Index ang_ind = 0; ang_ind < za_scat_grid.size(); ++ang_ind) {
      Numeric err = max_rel_error(pm[f_ind, ang_ind], scattering::expand_phase_matrix(pm_ref[1, f_ind, ang_ind]));
      if (err > 1e-3) return false;
    }

    for (Index stokes_ind = 0; stokes_ind < 4; ++stokes_ind) {
      Numeric ref = em_ref[1, f_ind, stokes_ind];
      Numeric rel_diff = std::abs((bulk_props.extinction_matrix[f_ind, stokes_ind, stokes_ind] - ref) / ref);
      if (rel_diff > 1e-3) return false;
    }

    Numeric err = max_rel_error(bulk_props.absorption_vector[f_ind], av_ref[1, f_ind]);
    if (err > 1e-3) return false;
  }
  return true;

  ///
  /// Test frequency interpolation
  ///
  Vector new_f_grid = {f_grid[1]};
  bulk_props = scattering_habit.get_bulk_scattering_properties_tro_gridded(point, new_f_grid);

  /// Ensure relative errors are small.
  Index f_ind = 1;
  for (Index ang_ind = 0; ang_ind < za_scat_grid.size(); ++ang_ind) {
    Numeric err = max_rel_error(pm[f_ind, ang_ind], scattering::expand_phase_matrix(pm_ref[1, f_ind, ang_ind]));
    if (err > 1e-3) return false;
  }

  for (Index stokes_ind = 0; stokes_ind < 4; ++stokes_ind) {
    Numeric ref = em_ref[1, f_ind, stokes_ind];
    Numeric rel_diff = std::abs((bulk_props.extinction_matrix[f_ind, stokes_ind, stokes_ind] - ref) / ref);
    if (rel_diff > 1e-3) return false;
  }

  Numeric err = max_rel_error(bulk_props.absorption_vector[f_ind], av_ref[1, f_ind]);
  if (err > 1e-3) return false;
  return true;

}

int main() {
  bool passed = false;
  std::cout << "Test calculate bulk properties TRO gridded: \t";
  passed = test_calculate_bulk_properties_tro_gridded();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }
  return 0;
}
