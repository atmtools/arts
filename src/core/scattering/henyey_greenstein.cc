#include "henyey_greenstein.h"

#include <cmath>
#include "math_funcs.h"

namespace scattering {

  HenyeyGreensteinScatterer::HenyeyGreensteinScatterer(ScatteringSpeciesProperty extinction,
                                     ScatteringSpeciesProperty single_scattering_albedo,
                                     const Numeric& g_) : extinction_property(extinction),
                                                          single_scattering_albedo_property(single_scattering_albedo),
                                                          g(g_){};


  template <Index stokes_dim>
  BulkScatteringProperties<Format::TRO, Representation::Gridded, stokes_dim>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_gridded(const AtmPoint& point,
                                                               std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) {
    auto t_grid = std::make_shared<Vector>(Vector{0.0});
    auto f_grid = std::make_shared<Vector>(Vector{0.0});

    PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded, stokes_dim> pm{t_grid,
                                                                                  f_grid,
                                                                                  zenith_angle_grid};
    pm *= 0.0;
    auto zenith_angles = grid_vector(*zenith_angle_grid);
    Numeric g2 = g * g;
    for (Index ind = 0; ind < zenith_angles.size(); ++ind) {
      pm(0, 0, ind, 0) = (1.0 - g2);
      pm(0, 0, ind, 0) /= std::pow(1.0 + g2 - 2.0 * g * cos(Conversion::deg2rad(zenith_angles[ind])), 3.0/2.0);
    }

    auto single_scattering_albedo = point[single_scattering_albedo_property];
    auto extinction = point[extinction_property];
    auto scattering_xsec = extinction * single_scattering_albedo;
    pm *= scattering_xsec / (4.0 * Constant::pi);

    ExtinctionMatrixData<Numeric, Format::TRO, Representation::Gridded, stokes_dim> emd{t_grid, f_grid};
    emd(0, 0, 0) = extinction;

    AbsorptionVectorData<Numeric, Format::TRO, Representation::Gridded, stokes_dim> av{t_grid, f_grid};
    av(0, 0, 0) = extinction - scattering_xsec;

    return BulkScatteringProperties<Format::TRO, Representation::Gridded, stokes_dim>{pm, emd, av};
  }

  template <Index stokes_dim>
  BulkScatteringProperties<Format::TRO, Representation::Spectral, stokes_dim>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_spectral(const AtmPoint& point, Index l) {
    auto t_grid = std::make_shared<Vector>(Vector{0.0});
    auto f_grid = std::make_shared<Vector>(Vector{0.0});

    PhaseMatrixData<Numeric, Format::TRO, Representation::Spectral, stokes_dim> pm(t_grid, f_grid, sht::provider.get_instance_lm(l, 0));
    pm *= 0.0;
    for (Index ind = 0; ind <= l; ++ind) {
      pm(0, 0, ind, 0) = sqrt((2.0 * static_cast<Numeric>(ind) + 1.0) * 4.0 * Constant::pi) * std::pow(g, static_cast<Numeric>(ind));
    }
    auto single_scattering_albedo = point[single_scattering_albedo_property];
    auto extinction = point[extinction_property];
    auto scattering_xsec = extinction * single_scattering_albedo;
    pm *= scattering_xsec / (4.0 * Constant::pi);

    ExtinctionMatrixData<Numeric, Format::TRO, Representation::Spectral, stokes_dim> emd{t_grid, f_grid};
    emd(0, 0, 0) = extinction;

    AbsorptionVectorData<Numeric, Format::TRO, Representation::Spectral, stokes_dim> av{t_grid, f_grid};
    av(0, 0, 0) = extinction - scattering_xsec;

    return BulkScatteringProperties<Format::TRO, Representation::Spectral, stokes_dim>{pm, emd, av};
  }

  std::ostream& operator<<(std::ostream& os,
                           const HenyeyGreensteinScatterer& scatterer) {
    return os << "HenyeyGreensteinScatterer(g = " << scatterer.g << ")";
  }

  template BulkScatteringProperties<Format::TRO, Representation::Gridded, 1>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_gridded(const AtmPoint& point,
                                                               std::shared_ptr<ZenithAngleGrid> zenith_angle_grid);
  template BulkScatteringProperties<Format::TRO, Representation::Gridded, 2>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_gridded(const AtmPoint& point,
                                                               std::shared_ptr<ZenithAngleGrid> zenith_angle_grid);
  template BulkScatteringProperties<Format::TRO, Representation::Gridded, 3>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_gridded(const AtmPoint& point,
                                                               std::shared_ptr<ZenithAngleGrid> zenith_angle_grid);
  template BulkScatteringProperties<Format::TRO, Representation::Gridded, 4>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_gridded(const AtmPoint& point,
                                                               std::shared_ptr<ZenithAngleGrid> zenith_angle_grid);

  template BulkScatteringProperties<Format::TRO, Representation::Spectral, 1>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_spectral(const AtmPoint& point,
                                                                Index l);
  template BulkScatteringProperties<Format::TRO, Representation::Spectral, 2>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_spectral(const AtmPoint& point,
                                                                Index l);
  template BulkScatteringProperties<Format::TRO, Representation::Spectral, 3>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_spectral(const AtmPoint& point,
                                                                Index l);
  template BulkScatteringProperties<Format::TRO, Representation::Spectral, 4>
  HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_spectral(const AtmPoint& point,
                                                                Index l);

}
