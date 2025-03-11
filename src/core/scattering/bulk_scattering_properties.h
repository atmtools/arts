#ifndef BULK_SCATTERING_PROPERTIES_H_
#define BULK_SCATTERING_PROPERTIES_H_

#include "phase_matrix.h"
#include "extinction_matrix.h"
#include "absorption_vector.h"
#include <rtepack.h>



namespace scattering {


  template<Format format, Representation repr>
  using PhaseMatrix = PhaseMatrixData<Numeric, format, repr>;

  template<Format format, Representation repr>
  using ExtinctionMatrix = ExtinctionMatrixData<Numeric, format, repr>;

  template<Format format, Representation repr>
  using AbsorptionVector = AbsorptionVectorData<Numeric, format, repr>;

  template <Format format, Representation repr>
    struct BulkScatteringProperties {
      std::optional<PhaseMatrix<format, repr>> phase_matrix;
      ExtinctionMatrix<format, repr> extinction_matrix;
      AbsorptionVector<format, repr> absorption_vector;

    BulkScatteringProperties<Format::ARO, Representation::Gridded>
    to_lab_frame(std::shared_ptr<const Vector> za_inc_grid,
                                          std::shared_ptr<const Vector> delta_aa_grid,
                                          std::shared_ptr<const ZenithAngleGrid> za_scat_grid_new) {
      return BulkScatteringProperties<Format::ARO, Representation::Gridded>{
        phase_matrix.transform([&](const PhaseMatrix<format, repr>& pm) {return pm.to_lab_frame(za_inc_grid, delta_aa_grid, za_scat_grid_new);}),
        extinction_matrix.to_lab_frame(delta_aa_grid),
        absorption_vector.to_lab_frame(delta_aa_grid)
      };
    }

    BulkScatteringProperties<format, Representation::Spectral>
    to_spectral() {
      return BulkScatteringProperties<format, Representation::Spectral>{
        phase_matrix.transform([&](const PhaseMatrix<format, repr>& pm) {return pm.to_spectral();}),
        extinction_matrix,
        absorption_vector
      };
    }

    BulkScatteringProperties<format, Representation::Spectral>
    to_spectral(Index degree, Index order) {
      return BulkScatteringProperties<format, Representation::Spectral>{
        phase_matrix.transform([&](const PhaseMatrix<format, repr>& pm) {return pm.to_spectral(degree, order);}),
        extinction_matrix.to_spectral(),
        absorption_vector.to_spectral()
      };
    }

    BulkScatteringProperties& operator+=(const BulkScatteringProperties& other) {
      if (phase_matrix.has_value()) {
        if (!other.phase_matrix.has_value()) {
          ARTS_USER_ERROR("Phase matrix missing in calculation of bulk scattering properties.");
        }
        *phase_matrix += *other.phase_matrix;
      }
      extinction_matrix += other.extinction_matrix;
      absorption_vector += other.absorption_vector;
      return *this;
    }
  };

  //////////////////////////////////////////////////////////////////////////////
  //  TRO bulk scattering properties in gridded format
  //////////////////////////////////////////////////////////////////////////////

  struct BulkScatteringPropertiesTROGridded {

    using PhaseMatrixVector = Tensor4;
    using ExtinctionMatrixVector = Tensor3;
    using AbsorptionVectorVector = Matrix;

    std::optional<PhaseMatrixVector> phase_matrix;
    ExtinctionMatrixVector extinction_matrix;
    AbsorptionVectorVector absorption_vector;

    inline BulkScatteringPropertiesTROGridded operator+(const BulkScatteringPropertiesTROGridded &other) const {
      std::optional<PhaseMatrixVector> new_phase_matrix = phase_matrix;
      if (phase_matrix.has_value()) {
        if (other.phase_matrix.has_value()) {
          *new_phase_matrix += other.phase_matrix.value();
        }
      }
      auto new_extinction_matrix = extinction_matrix;
      new_extinction_matrix += other.extinction_matrix;
      auto new_absorption_vector = absorption_vector;
      new_absorption_vector += other.absorption_vector;
      return BulkScatteringPropertiesTROGridded(new_phase_matrix,
                                                new_extinction_matrix,
                                                new_absorption_vector);
    }

    inline BulkScatteringPropertiesTROGridded& operator+=(const BulkScatteringPropertiesTROGridded &other)  {
      if (phase_matrix.has_value()) {
        if (other.phase_matrix.has_value()) {
          *phase_matrix += other.phase_matrix.value();
        }
      }
      extinction_matrix += other.extinction_matrix;
      absorption_vector += other.absorption_vector;
      return *this;
    }

    inline BulkScatteringPropertiesTROGridded& operator*=(Numeric fac)  {
      if (phase_matrix.has_value()) {
        *phase_matrix *= fac;
      }
      extinction_matrix *= fac;
      absorption_vector *= fac;
      return *this;
    }
  };

  //////////////////////////////////////////////////////////////////////////////
  //  TRO bulk scattering properties in spectral format
  //////////////////////////////////////////////////////////////////////////////

  struct BulkScatteringPropertiesTROSpectral {

      using PhaseMatrixVector = ComplexTensor4;
      using ExtinctionMatrixVector = Tensor3;
      using AbsorptionVectorVector = Matrix;

      std::optional<PhaseMatrixVector> phase_matrix;
      ExtinctionMatrixVector extinction_matrix;
      AbsorptionVectorVector absorption_vector;

      inline BulkScatteringPropertiesTROSpectral operator+(const BulkScatteringPropertiesTROSpectral &other) const {
        std::optional<PhaseMatrixVector> new_phase_matrix = phase_matrix;
        if (phase_matrix.has_value()) {
          if (other.phase_matrix.has_value()) {
            *new_phase_matrix += other.phase_matrix.value();
          }
        }
        auto new_extinction_matrix = extinction_matrix;
        new_extinction_matrix += other.extinction_matrix;
        auto new_absorption_vector = absorption_vector;
        new_absorption_vector += other.absorption_vector;
        return BulkScatteringPropertiesTROSpectral(new_phase_matrix,
                                                  new_extinction_matrix,
                                                  new_absorption_vector);
      }

      inline BulkScatteringPropertiesTROSpectral& operator+=(const BulkScatteringPropertiesTROSpectral &other)  {
        if (phase_matrix.has_value()) {
          if (other.phase_matrix.has_value()) {
            *phase_matrix += other.phase_matrix.value();
          }
        }
        extinction_matrix += other.extinction_matrix;
        absorption_vector += other.absorption_vector;
        return *this;
      }

      inline BulkScatteringPropertiesTROSpectral& operator*=(Numeric fac)  {
        if (phase_matrix.has_value()) {
          *phase_matrix *= fac;
        }
        extinction_matrix *= fac;
        absorption_vector *= fac;
        return *this;
      }
  };

}


#endif // BULK_SCATTERING_PROPERTIES_H_
