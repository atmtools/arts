#ifndef HENYEY_GREENSTEIN_H_
#define HENYEY_GREENSTEIN_H_

#include <configtypes.h>
#include <operators.h>

#include <functional>
#include <tuple>

#include "atm.h"
#include "bulk_scattering_properties.h"
#include "properties.h"

using ExtSSACallback =
    CustomOperator<std::pair<Numeric, Numeric>, Numeric, const AtmPoint&>;

namespace scattering {

struct ExtinctionSSALookup {
  ScatteringSpeciesProperty extinction_field{};
  ScatteringSpeciesProperty ssa_field{};
  ExtinctionSSALookup(ScatteringSpeciesProperty extinction_field_,
                      ScatteringSpeciesProperty ssa_field_)
      : extinction_field(std::move(extinction_field_)),
        ssa_field(std::move(ssa_field_)) {}

  std::pair<Numeric, Numeric> operator()(Numeric, const AtmPoint& atm_point) {
    return {atm_point[extinction_field], atm_point[ssa_field]};
  }
};

class HenyeyGreensteinScatterer {
  ExtSSACallback ext_ssa_callback{};

  Numeric g = 0.0;

 public:
  HenyeyGreensteinScatterer() = default;
  HenyeyGreensteinScatterer(ExtSSACallback ext_ssa_callback, const Numeric& g);
  HenyeyGreensteinScatterer(ScatteringSpeciesProperty extinction_field,
                            ScatteringSpeciesProperty ssa_field,
                            const Numeric& g_);

  HenyeyGreensteinScatterer(const HenyeyGreensteinScatterer&) = default;
  HenyeyGreensteinScatterer(HenyeyGreensteinScatterer&&) noexcept= default;
  HenyeyGreensteinScatterer& operator=(const HenyeyGreensteinScatterer&) =
      default;
  HenyeyGreensteinScatterer& operator=(HenyeyGreensteinScatterer&&) noexcept=
      default;

  template <Index stokes_dim>
  BulkScatteringProperties<Format::TRO, Representation::Gridded, stokes_dim>
  get_bulk_scattering_properties_tro_gridded(
      const AtmPoint&,
      const Vector& f_grid,
      std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const;

  template <Index stokes_dim>
  BulkScatteringProperties<Format::TRO, Representation::Spectral, stokes_dim>
  get_bulk_scattering_properties_tro_spectral(const AtmPoint&,
                                              const Vector& f_grid,
                                              Index l) const;

  template <Index stokes_dim>
  BulkScatteringProperties<scattering::Format::ARO,
                           scattering::Representation::Gridded,
                           stokes_dim>
  get_bulk_scattering_properties_aro_gridded(
      const AtmPoint&,
      const Vector& f_grid,
      const Vector& za_inc_grid,
      const Vector& delta_aa_grid,
      std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const;
  template <Index stokes_dim>
  BulkScatteringProperties<scattering::Format::ARO,
                           scattering::Representation::Spectral,
                           stokes_dim>
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

}  // namespace scattering

#endif  // HENYEY_GREENSTEIN_H_
