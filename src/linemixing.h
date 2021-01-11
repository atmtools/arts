#ifndef linemixing_h
#define linemixing_h

#include "absorption.h"
#include "complex.h"
#include "constants.h"
#include "gridded_fields.h"
#include "linescaling.h"
#include "zeemandata.h"

namespace Absorption::LineMixing {
//! Contains recomputed equivalent lines (sorting is unknown)
struct EquivalentLines {
  ComplexVector val;
  ComplexVector str;
  
  /*! Construct from known size
   * 
   * @param[in] n The size of the problem
   */
  explicit EquivalentLines(Index n=0) noexcept : val(n, 0), str(n, 0) {}
  
  /*! Construct from known parameters
   * 
   * Note that W can be renormalized in frequency
   * 
   * @param[in] W The relaxation matrix
   * @param[in] pop The population distributions
   * @param[in] dip The dipoles
   */
  EquivalentLines(const ComplexMatrix& W, const Vector& pop, const Vector& dip) noexcept;
  
  //! Explicitly deleted
  EquivalentLines(const EquivalentLines&) = delete;
  
  //! Explicitly defaulted
  EquivalentLines(EquivalentLines&&) = default;
  
  //! Explicitly deleted
  EquivalentLines& operator=(const EquivalentLines&) = delete;
  
  //! Explicitly defaulted
  EquivalentLines& operator=(EquivalentLines&&) = default;
};  // EquivalentLines

//! Contains the population distribution and dipole
struct PopulationAndDipole {
  Vector pop;
  Vector dip;
  
  /*! Construct from known parameters
   * 
   * @param[in] T The temperature
   * @param[in] band The absorption band
   * @param[in] partition_type The type of partition function data
   * @param[in] partition_data The partition function data
   */
  PopulationAndDipole(const Numeric T,
                      const AbsorptionLines& band,
                      const SpeciesAuxData::AuxType& partition_type,
                      const ArrayOfGriddedField1& partition_data) noexcept;
  
  /*! Sort self by f0*pop*dip^2 and returns positions of sorted values in the original
   * 
   * @param[in] band The absorption band
   * @return Reshuffling of positions from [0 ... N-1] to new positions
   */
  ArrayOfIndex sort(const AbsorptionLines& band) noexcept;
  
  /*! Sort self by pre-sorted sorting mechanism
   * 
   * @param[in] presorting Changes positions from [0 ... N-1] to presorting positions
   */
  void sort(const ArrayOfIndex& presorting) noexcept;
};  // PopulationAndDipole


/*! Computed the Error Corrected Sudden Complex absorption
 * 
 * @param[in] T The temperature
 * @param[in] P The pressure
 * @param[in] this_vmr The VMR of this species
 * @param[in] vmrs The VMRs of all broadeners of the absorption band
 * @param[in] mass The mass of all broadeners of the absorption band
 * @param[in] f_grid The grid of frequencies
 * @param[in] band The absorption band
 * @param[in] partition_type The type of partition function data
 * @param[in] partition_data The partition function data
 * @return Complex absorption
 */
ComplexVector ecs_absorption(const Numeric T,
                             const Numeric P,
                             const Numeric this_vmr,
                             const Vector& vmrs,
                             const Vector& mass,
                             const Vector& f_grid,
                             const AbsorptionLines& band,
                             const SpeciesAuxData::AuxType& partition_type,
                             const ArrayOfGriddedField1& partition_data);


/*! Computed the Error Corrected Sudden Complex absorption with Zeeman effect perturbations
 * 
 * Note that Zeeman perturbations are only applied after the ECS computations
 * 
 * @param[in] T The temperature
 * @param[in] P The pressure
 * @param[in] this_vmr The VMR of this species
 * @param[in] vmrs The VMRs of all broadeners of the absorption band
 * @param[in] mass The mass of all broadeners of the absorption band
 * @param[in] f_grid The grid of frequencies
 * @param[in] zeeman_polarization The Zeeman polarization to consider
 * @param[in] band The absorption band
 * @param[in] partition_type The type of partition function data
 * @param[in] partition_data The partition function data
 * @return Complex absorption of the Zeeman component
 */
ComplexVector ecs_absorption_with_zeeman_perturbations(const Numeric T,
                                                       const Numeric H,
                                                       const Numeric P,
                                                       const Numeric this_vmr,
                                                       const Vector& vmrs,
                                                       const Vector& mass,
                                                       const Vector& f_grid,
                                                       const Zeeman::Polarization zeeman_polarization,
                                                       const AbsorptionLines& band,
                                                       const SpeciesAuxData::AuxType& partition_type,
                                                       const ArrayOfGriddedField1& partition_data);


/*! Adapts the band to use Rosenkranz parameters
 * 
 * This function does not work properly and using it will result in
 * bad parameters
 * 
 * @param[in] band The absorption band
 * @param[in] temperatures The temperature grid for fitting parameters upon
 * @param[in] mass The mass of all broadeners of the absorption band
 * @param[in] partition_type The type of partition function data
 * @param[in] partition_data The partition function data
 * @return EXIT_FAILURE when some parameterization fit fails
 * @return EXIT_SUCCESS if all algorithms worked (independent of if the absorption will be reasonable)
 */
Index ecs_rosenkranz_adaptation(AbsorptionLines& band,
                                const Vector& temperatures,
                                const Vector& mass,
                                const SpeciesAuxData::AuxType& partition_type,
                                const ArrayOfGriddedField1& partition_data);
}  // LineMixing

#endif  // linemixing_h
