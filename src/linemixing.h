#ifndef linemixing_h
#define linemixing_h

#include "absorption.h"
#include "complex.h"
#include "constants.h"
#include "gridded_fields.h"
#include "linescaling.h"

namespace Absorption::LineMixing {
struct EquivalentLines {
  ComplexVector val;
  ComplexVector str;
  
  explicit EquivalentLines(Index n=0) noexcept : val(n, 0), str(n, 0) {}
  EquivalentLines(const ComplexMatrix& W,
                  const Vector& pop,
                  const Vector& dip) noexcept;
  EquivalentLines(const EquivalentLines&) = delete;
  EquivalentLines(EquivalentLines&&) = default;
  EquivalentLines& operator=(const EquivalentLines&) = delete;
  EquivalentLines& operator=(EquivalentLines&&) = default;
};  // EquivalentLines

struct PopulationAndDipole {
  Vector pop;
  Vector dip;
  
  PopulationAndDipole(const Numeric T,
                      const AbsorptionLines& band,
                      const SpeciesAuxData::AuxType& partition_type,
                      const ArrayOfGriddedField1& partition_data) noexcept : pop(band.NumLines()), dip(band.NumLines())
  {
    const Index N = band.NumLines();
    
    const Numeric QT = single_partition_function(T, partition_type, partition_data);
    const Numeric QT0 = single_partition_function(band.T0(), partition_type, partition_data);
    const Numeric ratiopart = QT0 / QT;
    
    for (Index i=0; i<N; i++) {
      const Numeric pop0 = (band.g_upp(i) / QT0) * boltzman_factor(band.T0(), band.E0(i));
      pop[i] = pop0 * ratiopart * boltzman_ratio(T, band.T0(), band.E0(i));
      dip[i] = std::sqrt(band.I0(i)/(pop0 * band.F0(i) * (1-stimulated_emission(band.T0(), band.F0(i)))));
    }
  }
  
  //! Sort self by f0*pop*dip^2 and returns positions of sorted values in the original */
  ArrayOfIndex sort(const AbsorptionLines& band) noexcept;
};  // PopulationAndDipole

ComplexVector linemixing_ecs_absorption(const Numeric T,
                                        const Numeric P,
                                        const Numeric this_vmr,
                                        const Vector& vmrs,
                                        const Vector& mass,
                                        const Vector& f_grid,
                                        const AbsorptionLines& band,
                                        const SpeciesAuxData::AuxType& partition_type,
                                        const ArrayOfGriddedField1& partition_data);
}  // Absorption::LineMixing 

#endif  // linemixing_h
