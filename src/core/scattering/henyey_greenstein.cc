#include "henyey_greenstein.h"

#include <cmath>
#include "math_funcs.h"

namespace scattering {

  HenyeyGreensteinScatterer::HenyeyGreensteinScatterer(ExtSSACallback ext_ssa_callback_,
                                                       const Numeric& g_)
    : ext_ssa_callback(ext_ssa_callback_),
      g(g_){};


  template <Index stokes_dim>
  BulkScatteringProperties<Format::TRO, Representation::Gridded, stokes_dim>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_gridded(const AtmPoint& atm_point,
                                                                        const Vector& f_grid,
                                                                        std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const {
    auto t_grid = std::make_shared<Vector>(Vector{0.0});
    auto f_grid_ptr = std::make_shared<Vector>(f_grid);

    PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded, stokes_dim> pm{t_grid,
                                                                                  f_grid_ptr,
                                                                                  zenith_angle_grid};
    ExtinctionMatrixData<Numeric, Format::TRO, Representation::Gridded, stokes_dim> emd{t_grid, f_grid_ptr};

    AbsorptionVectorData<Numeric, Format::TRO, Representation::Gridded, stokes_dim> av{t_grid, f_grid_ptr};

    auto zenith_angles = grid_vector(*zenith_angle_grid);
    for (Index f_ind = 0; f_ind < f_grid.size(); ++f_ind) {

      float extinction, ssa;
      std::tie(extinction, ssa) = ext_ssa_callback(f_grid[f_ind], atm_point);
      float scattering_xsec = extinction * ssa;

      emd(0, f_ind, 0) = extinction;
      av(0, f_ind, 0) = extinction - scattering_xsec;
      Numeric g2 = g * g;
      for (Index ind = 0; ind < zenith_angles.size(); ++ind) {
        pm(0, f_ind, ind, 0) = (1.0 - g2);
        pm(0, f_ind, ind, 0) /= std::pow(1.0 + g2 - 2.0 * g * cos(Conversion::deg2rad(zenith_angles[ind])), 3.0/2.0);
      }
      pm *= scattering_xsec / (4.0 * Constant::pi);

    }

    return BulkScatteringProperties<Format::TRO, Representation::Gridded, stokes_dim>{pm, emd, av};
  }

  template <Index stokes_dim>
  BulkScatteringProperties<Format::TRO, Representation::Spectral, stokes_dim>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_spectral(const AtmPoint& atm_point,
                                                                         const Vector& f_grid,
                                                                         Index l) const {
    auto t_grid = std::make_shared<Vector>(Vector{0.0});
    auto f_grid_ptr = std::make_shared<Vector>(f_grid);

    PhaseMatrixData<Numeric, Format::TRO, Representation::Spectral, stokes_dim> pm{t_grid,
                                                                                   f_grid_ptr,
                                                                                   sht::provider.get_instance_lm(l, 0)};
    ExtinctionMatrixData<Numeric, Format::TRO, Representation::Spectral, stokes_dim> emd{t_grid, f_grid_ptr};

    AbsorptionVectorData<Numeric, Format::TRO, Representation::Spectral, stokes_dim> av{t_grid, f_grid_ptr};

    for (Index f_ind = 0; f_ind < f_grid.size(); ++f_ind) {

      float extinction, ssa;
      std::tie(extinction, ssa) = ext_ssa_callback(f_grid[f_ind], atm_point);
      float scattering_xsec = extinction * ssa;

      for (Index ind = 0; ind <= l; ++ind) {
        pm(0, f_ind, ind, 0) = sqrt((2.0 * static_cast<Numeric>(ind) + 1.0) * 4.0 * Constant::pi) * std::pow(g, static_cast<Numeric>(ind));
      }

      pm *= scattering_xsec / (4.0 * Constant::pi);

      emd(0, f_ind, 0) = extinction;
      av(0, f_ind, 0) = extinction - scattering_xsec;
    }
    return BulkScatteringProperties<Format::TRO, Representation::Spectral, stokes_dim>{pm, emd, av};
  }

  template <Index stokes_dim>
  BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded, stokes_dim>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_gridded(const AtmPoint& atm_point,
                                                                        const Vector& f_grid,
                                                                        const Vector& za_inc_grid,
                                                                        const Vector& delta_aa_grid,
                                                                        std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const {
    auto bsp_tro = get_bulk_scattering_properties_tro_gridded<stokes_dim>(atm_point, f_grid, za_scat_grid);
    return bsp_tro.to_lab_frame(std::make_shared<Vector>(za_inc_grid),
                                std::make_shared<Vector>(delta_aa_grid),
                                za_scat_grid);
  }

  template <Index stokes_dim>
  BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Spectral, stokes_dim>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_spectral(const AtmPoint& atm_point,
                                                                         const Vector& f_grid,
                                                                         const Vector& za_inc_grid,
                                                                         Index degree,
                                                                         Index order) const {
    auto sht_ptr = sht::provider.get_instance(degree, order);
    auto aa_scat_grid_ptr = sht_ptr->get_aa_grid_ptr();
    auto za_scat_grid_ptr = std::make_shared<ZenithAngleGrid>(sht_ptr->get_zenith_angle_grid());
    auto bsp_tro = get_bulk_scattering_properties_tro_gridded<stokes_dim>(atm_point, f_grid, za_scat_grid_ptr);
    auto bsp_aro = bsp_tro.to_lab_frame(std::make_shared<Vector>(za_inc_grid), aa_scat_grid_ptr, za_scat_grid_ptr);
    return bsp_aro.to_spectral(degree, order);
  }

  std::ostream& operator<<(std::ostream& os,
                           const HenyeyGreensteinScatterer& scatterer) {
    return os << "HenyeyGreensteinScatterer(g = " << scatterer.g << ")";
  }

  template BulkScatteringProperties<Format::TRO, Representation::Gridded, 1>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_gridded(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                        std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const;
  template BulkScatteringProperties<Format::TRO, Representation::Gridded, 2>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_gridded(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                        std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const;
  template BulkScatteringProperties<Format::TRO, Representation::Gridded, 3>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_gridded(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                        std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const;
  template BulkScatteringProperties<Format::TRO, Representation::Gridded, 4>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_gridded(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                        std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const;

  template BulkScatteringProperties<Format::TRO, Representation::Spectral, 1>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_spectral(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                         Index l) const;
  template BulkScatteringProperties<Format::TRO, Representation::Spectral, 2>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_spectral(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                         Index l) const;
  template BulkScatteringProperties<Format::TRO, Representation::Spectral, 3>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_spectral(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                         Index l) const;
  template BulkScatteringProperties<Format::TRO, Representation::Spectral, 4>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_spectral(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                         Index l) const;

  template BulkScatteringProperties<Format::ARO, Representation::Gridded, 1>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_gridded(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                        const Vector& za_inc_grid,
                                                                        const Vector& delta_aa_grid,
                                                                        std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const;
  template BulkScatteringProperties<Format::ARO, Representation::Gridded, 2>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_gridded(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                        const Vector& za_inc_grid,
                                                                        const Vector& delta_aa_grid,
                                                                        std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const;
  template BulkScatteringProperties<Format::ARO, Representation::Gridded, 3>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_gridded(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                        const Vector& za_inc_grid,
                                                                        const Vector& delta_aa_grid,
                                                                        std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const;
  template BulkScatteringProperties<Format::ARO, Representation::Gridded, 4>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_gridded(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                        const Vector& za_inc_grid,
                                                                        const Vector& delta_aa_grid,
                                                                        std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const;

  template BulkScatteringProperties<Format::ARO, Representation::Spectral, 1>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_spectral(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                        const Vector& za_inc_grid,
                                                                        Index degree,
                                                                        Index order) const;
  template BulkScatteringProperties<Format::ARO, Representation::Spectral, 2>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_spectral(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                        const Vector& za_inc_grid,
                                                                        Index degree,
                                                                        Index order) const;
  template BulkScatteringProperties<Format::ARO, Representation::Spectral, 3>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_spectral(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                        const Vector& za_inc_grid,
                                                                        Index degree,
                                                                        Index order) const;
  template BulkScatteringProperties<Format::ARO, Representation::Spectral, 4>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_spectral(const AtmPoint& point,
                                                                        const Vector& f_grid,
                                                                        const Vector& za_inc_grid,
                                                                        Index degree,
                                                                        Index order) const;
}
