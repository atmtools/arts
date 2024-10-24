#ifndef HENYEY_GREENSTEIN_H_
#define HENYEY_GREENSTEIN_H_

#include <configtypes.h>
#include <functional>
#include <tuple>

#include "atm.h"
#include "properties.h"
#include "bulk_scattering_properties.h"


namespace scattering {


struct ExtinctionSSALookup {
  ScatteringSpeciesProperty extinction_field;
  ScatteringSpeciesProperty ssa_field;
  ExtinctionSSALookup(ScatteringSpeciesProperty extinction_field_,
                      ScatteringSpeciesProperty ssa_field_) :
    extinction_field(extinction_field_),
    ssa_field(ssa_field_)
    {}

  std::tuple<Numeric, Numeric> operator()(Numeric, const AtmPoint& atm_point) {
    auto extinction = atm_point[extinction_field];
    auto ssa = atm_point[ssa_field];
    return std::make_tuple(extinction, ssa);
  }

};


class HenyeyGreensteinScatterer {

  using ExtSSACallback = std::function<std::tuple<Numeric, Numeric>(Numeric, const AtmPoint&)>;
  ExtSSACallback ext_ssa_callback;

  Numeric g = 0.0;

 public:

  HenyeyGreensteinScatterer() {}
  HenyeyGreensteinScatterer(ExtSSACallback ext_ssa_callback, const Numeric &g);
  HenyeyGreensteinScatterer(ScatteringSpeciesProperty extinction_field, ScatteringSpeciesProperty ssa_field, const Numeric &g_) :
    ext_ssa_callback(ExtinctionSSALookup(extinction_field, ssa_field)),
    g(g_)
    {}

  template <Index stokes_dim>
  BulkScatteringProperties<Format::TRO, Representation::Gridded, stokes_dim>
  get_bulk_scattering_properties_tro_gridded(const AtmPoint&,
                                             const Vector& f_grid,
                                             std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const;

  template<Index stokes_dim>
  BulkScatteringProperties<Format::TRO, Representation::Spectral, stokes_dim>
  get_bulk_scattering_properties_tro_spectral(const AtmPoint&,
                                              const Vector& f_grid,
                                              Index l) const;

  template <Index stokes_dim>
  BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded, stokes_dim>
  get_bulk_scattering_properties_aro_gridded(const AtmPoint&,
                                             const Vector& f_grid,
                                             const Vector& za_inc_grid,
                                             const Vector& delta_aa_grid,
                                             std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const;
  template <Index stokes_dim>
  BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Spectral, stokes_dim>
  get_bulk_scattering_properties_aro_spectral(const AtmPoint&,
                                              const Vector& f_grid,
                                              const Vector& za_inc_grid,
                                              Index degree,
                                              Index order) const;

  Numeric get_g() const { return g; };
  void set_g(const Numeric& g_) { g = g_; };

  friend std::ostream& operator<<(std::ostream& os,
                                  const HenyeyGreensteinScatterer& scatterer);
};


}




#endif // HENYEY_GREENSTEIN_H_
