#ifndef BULK_SCATTERING_PROPERTIES_H_
#define BULK_SCATTERING_PROPERTIES_H_

#include "phase_matrix.h"
#include "extinction_matrix.h"
#include "absorption_vector.h"


namespace scattering {

  template<Format format, Representation repr, Index stokes_dim>
  using PhaseMatrix = PhaseMatrixData<Numeric, format, repr, stokes_dim>;

  template<Format format, Representation repr, Index stokes_dim>
  using ExtinctionMatrix = ExtinctionMatrixData<Numeric, format, repr, stokes_dim>;

  template<Format format, Representation repr, Index stokes_dim>
  using AbsorptionVector = AbsorptionVectorData<Numeric, format, repr, stokes_dim>;

  template <Format format, Representation repr, Index stokes_dim>
    struct BulkScatteringProperties {
      std::optional<PhaseMatrix<format, repr, stokes_dim>> phase_matrix;
      ExtinctionMatrix<format, repr, stokes_dim> extinction_matrix;
      AbsorptionVector<format, repr, stokes_dim> absorption_vector;
    };

  

}


#endif // BULK_SCATTERING_PROPERTIES_H_
