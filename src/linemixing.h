#ifndef linemixing_h
#define linemixing_h

#include <functional>
#include <utility>

#include "absorption.h"
#include "arts_conversions.h"
#include "gridded_fields.h"
#include "linescaling.h"
#include "matpack_complex.h"
#include "species.h"
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
  
  void sort_by_frequency(Vector& f, const ArrayOfIndex& sorting);
  
  friend std::ostream& operator<<(std::ostream& os, const EquivalentLines& eqv);
};  // EquivalentLines

//! Contains the population distribution and dipole
struct PopulationAndDipole {
  Vector pop;
  Vector dip;
  
  /*! Construct from known parameters
   * 
   * @param[in] T The temperature
   * @param[in] band The absorption band
   */
  PopulationAndDipole(const Numeric T,
                      const AbsorptionLines& band);
  
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

using EnergyFunction = std::function<Numeric (const Rational)>;

struct SpeciesErrorCorrectedSuddenData {
  Species::Species spec{Species::Species::Bath};
  LineShapeModelParameters scaling;
  LineShapeModelParameters beta;
  LineShapeModelParameters lambda;
  LineShapeModelParameters collisional_distance;
  Numeric mass{1};
  
  constexpr SpeciesErrorCorrectedSuddenData(Species::Species inspec=Species::Species::Bath) noexcept :
    spec(inspec), scaling(LineShapeTemperatureModel::T0, 0, 0, 0, 0),
    beta(LineShapeTemperatureModel::T0, 0, 0, 0, 0),
    lambda(LineShapeTemperatureModel::T0, 0, 0, 0, 0),
    collisional_distance(LineShapeTemperatureModel::T0, 0, 0, 0, 0) {}
  
  [[nodiscard]] Numeric Q(const Rational J,
            const Numeric T,
            const Numeric T0,
            const Numeric energy) const;
  
  [[nodiscard]] Numeric Omega(const Numeric T,
                const Numeric T0,
                const Numeric other_mass,
                const Numeric energy_x,
                const Numeric energy_xm2) const;
  
  friend std::ostream& operator<<(std::ostream& os, const SpeciesErrorCorrectedSuddenData& srbd);
  
  friend std::istream& operator>>(std::istream& is, SpeciesErrorCorrectedSuddenData& srbd);
  
  constexpr bool operator==(Species::Species species) const noexcept {
    return species == spec;
  }
};  // SpeciesErrorCorrectedSuddenData

using ArrayOfSpeciesErrorCorrectedSuddenData = Array<SpeciesErrorCorrectedSuddenData>;

/** Rovibrational line mixing data following
 * the ideas of Collisional Effects On
 * Molecular Spectra by J.-M. Hartmann, C. Boulet,
 * and D. Robert, 1st edition, 2008.
 * 
 * This is for now a purely "data" struct, without
 * any logic.  If necessary, logic for the selection
 * of rotational energy calculations, Wigner symbols,
 * adiabatic factors and so on should be added here
 */
struct ErrorCorrectedSuddenData {
  QuantumIdentifier id;
  
  /** Data of species data
   * The program is considered ill-formed if data does not
   * contain a Bath-species, either the default one or modified */
  ArrayOfSpeciesErrorCorrectedSuddenData data;
  
  explicit ErrorCorrectedSuddenData(QuantumIdentifier qid={}) noexcept :
  id(std::move(qid)), data({SpeciesErrorCorrectedSuddenData()}) {}
  
  friend std::ostream& operator<<(std::ostream& os, const ErrorCorrectedSuddenData& rbd);
  
  friend std::istream& operator>>(std::istream& is, ErrorCorrectedSuddenData& rbd);
  
  bool operator==(const QuantumIdentifier& band_id) const noexcept {
    return band_id.part_of(id);
  }

  [[nodiscard]] Index pos(Species::Species spec) const;

  [[nodiscard]] Index size() const noexcept { return data.size(); }

  const SpeciesErrorCorrectedSuddenData& operator[](Species::Species spec) const;
  
  SpeciesErrorCorrectedSuddenData& operator[](Species::Species spec);
};  // ErrorCorrectedSuddenData

struct MapOfErrorCorrectedSuddenData : public Array<ErrorCorrectedSuddenData> {
  ErrorCorrectedSuddenData& operator[](const QuantumIdentifier& id);
  
  const ErrorCorrectedSuddenData& operator[](const QuantumIdentifier& id) const;
  
  const ErrorCorrectedSuddenData& operator[](Index i) const ARTS_NOEXCEPT {
    ARTS_ASSERT(i >= 0 and i < nelem())
    return * (begin() + i);
  }
  
  ErrorCorrectedSuddenData& operator[](Index i) ARTS_NOEXCEPT {
    ARTS_ASSERT(i >= 0 and i < nelem())
    return * (begin() + i);
  }
  
  friend std::ostream& operator<<(std::ostream& os, const MapOfErrorCorrectedSuddenData& m);
};  // MapOfErrorCorrectedSuddenData

// Return struct from calculations
struct EcsReturn {
  ComplexVector abs;
  ArrayOfComplexVector dabs;
  bool error;
  EcsReturn(EcsReturn&&) noexcept = default;
  EcsReturn(ComplexVector&& abs_, ArrayOfComplexVector&& dabs_, bool error_) noexcept :
  abs(std::move(abs_)), dabs(std::move(dabs_)), error(error_) {}
};


/*! Computed the Error Corrected Sudden Complex absorption with Zeeman effect perturbations
 * 
 * Note that Zeeman perturbations are only applied after the ECS computations, and only
 * if zeeman_polarization is not None
 * 
 * @param[in] T The temperature
 * @param[in] P The pressure
 * @param[in] this_vmr The VMR of this species
 * @param[in] vmrs The VMRs of all broadeners of the absorption band
 * @param[in] mass The mass of all broadeners of the absorption band
 * @param[in] f_grid The grid of frequencies
 * @param[in] zeeman_polarization The Zeeman polarization to consider
 * @param[in] band The absorption band
 * @return Complex absorption of the Zeeman component
 */
EcsReturn ecs_absorption(const Numeric T,
                         const Numeric H,
                         const Numeric P,
                         const Numeric this_vmr,
                          const Vector& vmrs,
                          const ErrorCorrectedSuddenData& ecs_data,
                          const Vector& f_grid,
                          const Zeeman::Polarization zeeman_polarization,
                          const AbsorptionLines& band,
                          const ArrayOfRetrievalQuantity& jacobian_quantities={});

/**  Adapts the band to the temperature data
 * 
 * @return EXIT_FAILURE on failure
 * @return EXIT_SUCCESS on success
 */
Index band_eigenvalue_adaptation(AbsorptionLines& band,
                                 const Tensor4& tempdata,
                                 const Vector& temperatures,
                                 const Numeric P0,
                                 const Index ord);


/** Adapts the relaxation matrix eigenvalues to a form
 * where they represent additions towards the three
 * Rosenkranz parameters and the 3rd order pressure
 * broadening term.
 * 
 * The EquivalentLines are sorted by the real part of the
 * val-component.  At a sufficiently small pressure, this
 * is equivalent to sorting by line frequency.  At higher
 * pressures, the sorting might be wrong
 * 
 * The output is adapted to fit into standard LBL calculations
 * assuming that the pressure-sorting works.  This happens via
 * an adaptation so that
 * 
 * X = val - Complex(F0 + D0(T), G0(T))
 * Y  = str / I(T) - 1
 * 
 * are the two return values, where val and str are the output of
 * the standard constructor of EquivalentLines.
 */
EquivalentLines eigenvalue_adaptation_of_relmat(const ComplexMatrix& W,
                                                const Vector& pop,
                                                const Vector& dip,
                                                const AbsorptionLines& band,
                                                const Numeric frenorm,
                                                const Numeric T,
                                                const Numeric P,
                                                const Numeric QT,
                                                const Numeric QT0,
                                                const Index broadener);


/*! Adapts the band to use Rosenkranz parameters from Eigenvalue decomposition
 * 
 * Fits will be of form:
 * 
 * Y  ~ P      * polyfit(x, y,  LineShape::ModelParameters::N - 1)
 * DV ~ P ** 2 * polyfit(x, dv, LineShape::ModelParameters::N - 1)
 * G  ~ P ** 2 * polyfit(x, g,  LineShape::ModelParameters::N - 1)
 * DG ~ P ** 3 * polyfit(x, dg, LineShape::ModelParameters::N - 1)
 * 
 * where all the values are computed at P0 but re-normalized by the order above
 * 
 * Note that these parameters will fail at high pressures
 * 
 * @param[in] band The absorption band
 * @param[in] temperatures The temperature grid for fitting parameters upon
 * @param[in] mass The mass of all broadeners of the absorption band
 * @param[in] P0 The pressure at which temperature dependencies are computed at
 * @param[in] ord The order of the parameters [1: Y; 2: Y, DF, G; 3: Y, DF, G, DG], the last is still not supported fully
 * @param[in] robust Doesn't throw on failure if true
 * @param[in] rosenkranz_adaptation Makes the explicit computation of Rosenkranz parameters
 * @return EXIT_FAILURE when some parameterization fit fails
 * @return EXIT_SUCCESS if all algorithms worked (independent of if the absorption will be reasonable)
 */
void ecs_eigenvalue_adaptation(AbsorptionLines& band,
                               const Vector& temperatures,
                               const ErrorCorrectedSuddenData& ecs_data,
                               const Numeric P0,
                               const Index ord,
                               const bool robust,
                               const bool rosenkranz_adaptation);

/*! Outputs the adaptation values used for ecs_eigenvalue_adaptation but as 
 * a function of pressure.  ecs_eigenvalue_adaptation makes strong assumptions
 * about the pressure order of the outputs used, and if this work for a given
 * band can be found out using this data
 * 
 * @param[in] band The absorption band [N lines]
 * @param[in] temperatures The temperature grid for fitting parameters upon [K temperatures]
 * @param[in] mass The mass of all broadeners of the absorption band [M broadeners]
 * @param[in] pressures The pressures for testing [L pressures]
 * @return Tensor with size [4, N, M, K, L]
 */
Tensor5 ecs_eigenvalue_adaptation_test(const AbsorptionLines& band,
                                       const Vector& temperatures,
                                       const ErrorCorrectedSuddenData& ecs_data,
                                       const Vector& pressures);
} // namespace Absorption::LineMixing

using ErrorCorrectedSuddenData = Absorption::LineMixing::ErrorCorrectedSuddenData;
using MapOfErrorCorrectedSuddenData = Absorption::LineMixing::MapOfErrorCorrectedSuddenData;
using SpeciesErrorCorrectedSuddenData = Absorption::LineMixing::SpeciesErrorCorrectedSuddenData;
using Absorption::LineMixing::ArrayOfSpeciesErrorCorrectedSuddenData;

#endif  // linemixing_h
