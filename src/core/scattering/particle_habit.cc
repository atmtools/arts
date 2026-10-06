#include "particle_habit.h"

#include <matpack.h>
#include <sorting.h>

namespace scattering {

std::pair<Numeric, Numeric> derive_scat_species_a_and_b(const Vector&  sizes,
                                                        const Vector&  masses,
                                                        const Numeric& fit_start,
                                                        const Numeric& fit_end) {
  const Index nse = sizes.size();
  assert(nse > 1);

  ArrayOfIndex intarr_sort, intarr_unsort(0);
  Vector       x_unsorted(nse), m_unsorted(nse);
  Vector       q;
  Index        nsev = 0;

  for (Index i = 0; i < nse; i++) {
    if (std::isnan(sizes[i])) ARTS_USER_ERROR("NaN found in selected size grid data.");
    if (std::isnan(masses[i])) ARTS_USER_ERROR("NaN found among particle mass data.");

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

std::tuple<Vector, Numeric, Numeric> ParticleHabit::get_size_mass_info(SizeParameter  size_parameter,
                                                                       const Numeric& fit_start,
                                                                       const Numeric& fit_end) {
  Index n_particles = scattering_data.size();

  Vector sizes(n_particles);
  Vector masses(n_particles);

  for (Index ind = 0; ind < n_particles; ++ind) {
    auto mass = std::visit([](const auto& ssd) { return ssd.get_mass(); }, scattering_data[ind]);
    auto size =
        std::visit([&size_parameter](const auto& ssd) { return ssd.get_size(size_parameter); }, scattering_data[ind]);
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

ParticleHabit ParticleHabit::to_tro_spectral(const Vector& t_grid, const Vector& f_grid, Index l) {
  std::vector<SingleScatteringData<Numeric, Format::TRO, Representation::Spectral>> new_scat_data{};
  auto new_grids = ScatteringDataGrids(std::make_shared<Vector>(t_grid), std::make_shared<Vector>(f_grid));
  auto transform = [&new_grids, &l](const auto& ssd) { return ssd_to_tro_spectral(new_grids, l, ssd); };
  for (size_t p_ind = 0; p_ind < scattering_data.size(); ++p_ind) {
    new_scat_data.push_back(std::visit(transform, scattering_data[p_ind]));
  }
  return ParticleHabit(new_scat_data, new_grids);
}

ParticleHabit ParticleHabit::to_tro_gridded(const Vector&          t_grid,
                                            const Vector&          f_grid,
                                            const ZenithAngleGrid& za_scat_grid) {
  auto new_grids = ScatteringDataGrids(std::make_shared<const Vector>(t_grid),
                                       std::make_shared<const Vector>(f_grid),
                                       std::make_shared<const ZenithAngleGrid>(za_scat_grid));

  auto transform = [&new_grids](const auto& ssd) { return ssd_to_tro_gridded(new_grids, ssd); };
  std::vector<SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>> new_scattering_data;
  new_scattering_data.reserve(scattering_data.size());
  for (const ParticleData& pd : scattering_data) { new_scattering_data.push_back(std::visit(transform, pd)); }
  return ParticleHabit(new_scattering_data, new_grids);
}

std::pair<ParticleHabit, std::vector<LegendreReport>> ParticleHabit::to_tro_spectral_with_report(const Vector& t_grid,
                                                                                                 const Vector& f_grid,
                                                                                                 Index l) const {
  using Spectral = SingleScatteringData<Numeric, Format::TRO, Representation::Spectral>;
  using Gridded  = SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>;
  auto new_grids = ScatteringDataGrids(std::make_shared<Vector>(t_grid), std::make_shared<Vector>(f_grid));
  std::vector<Spectral>       new_scat_data;
  std::vector<LegendreReport> reports;
  for (const ParticleData& pd : scattering_data) {
    const auto* gridded = std::get_if<Gridded>(&pd);
    ARTS_USER_ERROR_IF(not gridded,
                       "A report on the Legendre conversion needs gridded TRO data, but the habit holds other data")
    auto [spectral, report] = gridded->to_spectral_with_report(l);
    new_scat_data.push_back(spectral.regrid(new_grids));
    reports.push_back(std::move(report));
  }
  return {ParticleHabit(new_scat_data, new_grids), std::move(reports)};
}

ParticleHabit ParticleHabit::to_aro_spectral(const Vector&          t_grid,
                                             const Vector&          f_grid,
                                             const Vector&          za_inc_grid,
                                             const ZenithAngleGrid& za_scat_grid,
                                             Index                  max_mode) const {
  auto new_grids = ScatteringDataGrids(std::make_shared<const Vector>(t_grid),
                                       std::make_shared<const Vector>(f_grid),
                                       std::make_shared<const Vector>(za_inc_grid),
                                       nullptr,
                                       std::make_shared<const ZenithAngleGrid>(za_scat_grid));
  auto transform = [&new_grids, max_mode](const auto& ssd) { return ssd_to_aro_spectral(new_grids, max_mode, ssd); };
  std::vector<SingleScatteringData<Numeric, Format::ARO, Representation::Spectral>> new_scattering_data;
  new_scattering_data.reserve(scattering_data.size());
  for (const ParticleData& pd : scattering_data) { new_scattering_data.push_back(std::visit(transform, pd)); }
  return ParticleHabit(new_scattering_data, new_grids);
}

ParticleHabit ParticleHabit::to_aro_gridded(
    const Vector& t_grid, const Vector& f_grid, const Vector&, const Vector& aa_scat_grid, const Vector& za_scat_grid) {
  auto new_grids = ScatteringDataGrids(std::make_shared<const Vector>(t_grid),
                                       std::make_shared<const Vector>(f_grid),
                                       std::make_shared<const Vector>(aa_scat_grid),
                                       std::make_shared<const Vector>(za_scat_grid),
                                       std::make_shared<const ZenithAngleGrid>(za_scat_grid));
  auto transform = [&new_grids](const auto& ssd) { return ssd_to_aro_gridded(new_grids, ssd); };
  std::vector<SingleScatteringData<Numeric, Format::ARO, Representation::Gridded>> new_scattering_data;
  new_scattering_data.reserve(scattering_data.size());
  for (const ParticleData& pd : scattering_data) { new_scattering_data.push_back(std::visit(transform, pd)); }
  return ParticleHabit(new_scattering_data, new_grids);
}

Vector ParticleHabit::get_sizes(SizeParameter param) const {
  Index  n_particles = scattering_data.size();
  Vector sizes(n_particles);
  for (Index ind = 0; ind < n_particles; ++ind) {
    auto size = std::visit([&param](const auto& ssd) { return ssd.get_size(param); }, scattering_data[ind]);
    if (size.has_value()) {
      sizes[ind] = size.value();
    } else {
      ARTS_USER_ERROR("Encountered particle without size information.");
    }
  }
  return sizes;
}

ParticleHabit ParticleHabit::from_legacy_tro(std::vector<::SingleScatteringData> ssd_,
                                             std::vector<::ScatteringMetaData>   meta_) {
  std::vector<SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>> ssd;
  ssd.reserve(ssd_.size());
  for (auto ind = 0; ind < Index(ssd_.size()); ++ind) {
    ssd.push_back(
        SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>::from_legacy_tro(ssd_[ind], meta_[ind]));
  }
  return ParticleHabit(ssd);
}

ParticleHabit ParticleHabit::from_legacy_aro(std::vector<::SingleScatteringData> ssd_,
                                             std::vector<::ScatteringMetaData>   meta_) {
  ARTS_USER_ERROR_IF(ssd_.size() != meta_.size(), "Scattering data and metadata sizes differ.")
  ARTS_USER_ERROR_IF(ssd_.empty(), "Cannot construct an empty particle habit.")
  std::vector<SingleScatteringData<Numeric, Format::ARO, Representation::Gridded>> converted;
  converted.reserve(ssd_.size());
  for (Index ind = 0; ind < static_cast<Index>(ssd_.size()); ++ind) {
    converted.push_back(
        SingleScatteringData<Numeric, Format::ARO, Representation::Gridded>::from_legacy_aro(ssd_[ind], meta_[ind]));
  }
  Vector signed_delta_aa(2 * ssd_[0].aa_grid.size() - 1);
  for (Size i = 1; i < ssd_[0].aa_grid.size(); ++i) {
    signed_delta_aa[ssd_[0].aa_grid.size() - 1 - i] = -ssd_[0].aa_grid[i];
  }
  for (Size i = 0; i < ssd_[0].aa_grid.size(); ++i) {
    signed_delta_aa[ssd_[0].aa_grid.size() - 1 + i] = ssd_[0].aa_grid[i];
  }
  auto grids = ScatteringDataGrids(std::make_shared<Vector>(ssd_[0].T_grid),
                                   std::make_shared<Vector>(ssd_[0].f_grid),
                                   std::make_shared<Vector>(ssd_[0].za_grid),
                                   std::make_shared<Vector>(std::move(signed_delta_aa)),
                                   std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(ssd_[0].za_grid)));
  return ParticleHabit(converted, grids);
}

ParticleHabit ParticleHabit::sphere(const StridedVectorView& t_grid,
                                    const StridedVectorView& f_grid,
                                    const StridedVectorView& diameters,
                                    const ZenithAngleGrid&   za_scat_grid,
                                    const ComplexMatrix&     refractive_index,
                                    Numeric                  density) {
  std::vector<SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>> ssd;
  ssd.reserve(diameters.size());
  for (auto diameter : diameters)
    ssd.push_back(SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>::sphere(
        t_grid, f_grid, diameter, za_scat_grid, refractive_index, density));
  return ParticleHabit(ssd);
}

ParticleHabit ParticleHabit::liquid_sphere(const StridedVectorView& t_grid,
                                           const StridedVectorView& f_grid,
                                           const StridedVectorView& diameters,
                                           const ZenithAngleGrid&   za_scat_grid) {
  std::vector<SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>> ssd;
  ssd.reserve(diameters.size());
  for (auto ind = 0; ind < Index(diameters.size()); ++ind) {
    ssd.push_back(SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>::liquid_sphere(
        t_grid, f_grid, diameters[ind], za_scat_grid));
  }
  return ParticleHabit(ssd);
}
}  // namespace scattering
