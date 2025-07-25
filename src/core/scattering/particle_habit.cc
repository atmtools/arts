#include "particle_habit.h"

#include "matpack.h"
#include "util/sorting.h"

namespace scattering {

std::pair<Numeric, Numeric> derive_scat_species_a_and_b(
    const Vector& sizes,
    const Vector& masses,
    const Numeric& fit_start,
    const Numeric& fit_end) {
  const Index nse = sizes.size();
  assert(nse > 1);

  ArrayOfIndex intarr_sort, intarr_unsort(0);
  Vector x_unsorted(nse), m_unsorted(nse);
  Vector q;
  Index nsev = 0;

  for (Index i = 0; i < nse; i++) {
    if (std::isnan(sizes[i]))
      ARTS_USER_ERROR("NaN found in selected size grid data.");
    if (std::isnan(masses[i]))
      ARTS_USER_ERROR("NaN found among particle mass data.");

    if (sizes[i] >= fit_start && sizes[i] <= fit_end) {
      x_unsorted[nsev]  = sizes[i];
      m_unsorted[nsev]  = masses[i];
      nsev             += 1;
    }
  }

  if (nsev < 2)
    ARTS_USER_ERROR(
        "Less than two size points found in the range "
        "[fit_start, fit_end]. It is then not possible "
        "to determine the a and b parameters.");

  get_sorted_indexes(intarr_sort, x_unsorted[Range(0, nsev)]);
  Vector log_x(nsev), log_m(nsev);

  for (Index i = 0; i < nsev; i++) {
    log_x[i] = log(x_unsorted[intarr_sort[i]]);
    log_m[i] = log(m_unsorted[intarr_sort[i]]);
  }

  linreg(q, log_x, log_m);
  return std::pair<Numeric, Numeric>(exp(q[0]), q[1]);
}

std::tuple<Vector, Numeric, Numeric> ParticleHabit::get_size_mass_info(
    SizeParameter size_parameter,
    const Numeric& fit_start,
    const Numeric& fit_end) {
  Index n_particles = scattering_data.size();

  Vector sizes(n_particles);
  Vector masses(n_particles);

  for (Index ind = 0; ind < n_particles; ++ind) {
    auto mass = std::visit([](const auto& ssd) { return ssd.get_mass(); },
                           scattering_data[ind]);
    auto size = std::visit(
        [&size_parameter](const auto& ssd) {
          return ssd.get_size(size_parameter);
        },
        scattering_data[ind]);
    if (mass.has_value()) {
      masses[ind] = mass.value();
      sizes[ind]  = size.value();
    } else {
      ARTS_USER_ERROR("Encountered particle without size information.");
    }
  }

  auto [a, b] = derive_scat_species_a_and_b(sizes, masses, fit_start, fit_end);
  return std::make_tuple(sizes, a, b);
}

ParticleHabit ParticleHabit::to_tro_spectral(const Vector& t_grid,
                                             const Vector& f_grid,
                                             Index l) {
  std::vector<
      SingleScatteringData<Numeric, Format::TRO, Representation::Spectral>>
      new_scat_data{};
  auto new_grids = ScatteringDataGrids(std::make_shared<Vector>(t_grid),
                                       std::make_shared<Vector>(f_grid));
  auto transform = [&new_grids, &l](const auto& ssd) {
    return ssd_to_tro_spectral(new_grids, l, ssd);
  };
  for (size_t p_ind = 0; p_ind < scattering_data.size(); ++p_ind) {
    new_scat_data.push_back(std::visit(transform, scattering_data[p_ind]));
  }
  return ParticleHabit(new_scat_data, new_grids);
}

ParticleHabit ParticleHabit::to_tro_gridded(
    const Vector& t_grid,
    const Vector& f_grid,
    const ZenithAngleGrid& za_scat_grid) {
  auto new_grids = ScatteringDataGrids(
      std::make_shared<const Vector>(t_grid),
      std::make_shared<const Vector>(f_grid),
      std::make_shared<const ZenithAngleGrid>(za_scat_grid));

  auto transform = [&new_grids](const auto& ssd) {
    return ssd_to_tro_gridded(new_grids, ssd);
  };
  std::vector<
      SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>>
      new_scattering_data;
  new_scattering_data.reserve(scattering_data.size());
  for (const ParticleData& pd : scattering_data) {
    new_scattering_data.push_back(std::visit(transform, pd));
  }
  return ParticleHabit(new_scattering_data, new_grids);
}

ParticleHabit ParticleHabit::to_aro_spectral(const Vector& t_grid,
                                             const Vector& f_grid,
                                             const Vector& za_inc_grid,
                                             Index l,
                                             Index m) {
  auto sht_ptr      = sht::provider.get_instance_lm(l, m);
  auto aa_scat_grid = sht_ptr->get_azimuth_angle_grid();
  auto za_scat_grid = sht_ptr->get_zenith_angle_grid();
  auto new_grids    = ScatteringDataGrids(
      std::make_shared<const Vector>(t_grid),
      std::make_shared<const Vector>(f_grid),
      std::make_shared<const Vector>(za_inc_grid),
      std::make_shared<const Vector>(aa_scat_grid),
      std::make_shared<const ZenithAngleGrid>(za_scat_grid));
  auto transform = [&new_grids, &l, &m](const auto& ssd) {
    return ssd_to_aro_spectral(new_grids, l, m, ssd);
  };
  std::vector<
      SingleScatteringData<Numeric, Format::ARO, Representation::Spectral>>
      new_scattering_data;
  new_scattering_data.reserve(scattering_data.size());
  for (const ParticleData& pd : scattering_data) {
    new_scattering_data.push_back(std::visit(transform, pd));
  }
  return ParticleHabit(new_scattering_data, new_grids);
}

ParticleHabit ParticleHabit::to_aro_gridded(const Vector& t_grid,
                                            const Vector& f_grid,
                                            const Vector&,
                                            const Vector& aa_scat_grid,
                                            const Vector& za_scat_grid) {
  auto new_grids = ScatteringDataGrids(
      std::make_shared<const Vector>(t_grid),
      std::make_shared<const Vector>(f_grid),
      std::make_shared<const Vector>(aa_scat_grid),
      std::make_shared<const Vector>(za_scat_grid),
      std::make_shared<const ZenithAngleGrid>(za_scat_grid));
  auto transform = [&new_grids](const auto& ssd) {
    return ssd_to_aro_gridded(new_grids, ssd);
  };
  std::vector<
      SingleScatteringData<Numeric, Format::ARO, Representation::Gridded>>
      new_scattering_data;
  new_scattering_data.reserve(scattering_data.size());
  for (const ParticleData& pd : scattering_data) {
    new_scattering_data.push_back(std::visit(transform, pd));
  }
  return ParticleHabit(new_scattering_data, new_grids);
}

Vector ParticleHabit::get_sizes(SizeParameter param) const {
  Index n_particles = scattering_data.size();
  Vector sizes(n_particles);
  for (Index ind = 0; ind < n_particles; ++ind) {
    auto size =
        std::visit([&param](const auto& ssd) { return ssd.get_size(param); },
                   scattering_data[ind]);
    if (size.has_value()) {
      sizes[ind] = size.value();
    } else {
      ARTS_USER_ERROR("Encountered particle without size information.");
    }
  }
  return sizes;
}

ParticleHabit ParticleHabit::from_legacy_tro(
    std::vector<::SingleScatteringData> ssd_,
    std::vector<::ScatteringMetaData> meta_) {
  std::vector<
      SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>>
      ssd;
  ssd.reserve(ssd_.size());
  for (auto ind = 0; ind < Index(ssd_.size()); ++ind) {
    ssd.push_back(
        SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>::
            from_legacy_tro(ssd_[ind], meta_[ind]));
  }
  return ParticleHabit(ssd);
}

ParticleHabit ParticleHabit::liquid_sphere(
    const StridedVectorView& t_grid,
    const StridedVectorView& f_grid,
    const StridedVectorView& diameters,
    const ZenithAngleGrid& za_scat_grid) {
  std::vector<
      SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>>
      ssd;
  ssd.reserve(diameters.size());
  for (auto ind = 0; ind < Index(diameters.size()); ++ind) {
    ssd.push_back(
        SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>::
            liquid_sphere(t_grid, f_grid, diameters[ind], za_scat_grid));
  }
  return ParticleHabit(ssd);
}
}  // namespace scattering
