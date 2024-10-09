#ifndef HENYEY_GREENSTEIN_H_
#define HENYEY_GREENSTEIN_H_

#include <configtypes.h>

#include "atm.h"
#include "properties.h"
#include "bulk_scattering_properties.h"


namespace scattering {


class HenyeyGreensteinScatterer {

  ScatteringSpeciesProperty extinction_property;
  ScatteringSpeciesProperty single_scattering_albedo_property;

  Numeric g = 0.0;

 public:

  HenyeyGreensteinScatterer() {}
  HenyeyGreensteinScatterer(ScatteringSpeciesProperty extinction, ScatteringSpeciesProperty single_scattering_albedo, const Numeric &g);

  template <Index stokes_dim>
  BulkScatteringProperties<Format::TRO, Representation::Gridded, stokes_dim>
  get_bulk_scattering_properties_tro_gridded(const AtmPoint&,
                                             std::shared_ptr<ZenithAngleGrid> zenith_angle_grid);

  template<Index stokes_dim>
  BulkScatteringProperties<Format::TRO, Representation::Spectral, stokes_dim>
  get_bulk_scattering_properties_tro_spectral(const AtmPoint&, Index l);

  Numeric get_g() const { return g; };
  void set_g(const Numeric& g_) { g = g_; };

  friend std::ostream& operator<<(std::ostream& os,
                                  const HenyeyGreensteinScatterer& scatterer);
};


}




#endif // HENYEY_GREENSTEIN_H_
