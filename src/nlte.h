/**
 * @file nlte.h
 * @author Richard Larsson
 * @date 2018-03-07
 * 
 * @brief Deep calculations for NLTE
 */

#ifndef NLTE_H
#define NLTE_H

#include "absorption.h"
#include "gridded_fields.h"
#include "matpack_data.h"
#include "quantum_numbers.h"

#include <unordered_map>

/** The data type for additioanl NLTE vibrational energy levels
 *
 * Note that this should just implement some of unordered_map
 * but it cannot be an unordered_map because that disagrees
 * with some functionality further down the chain (2023-01-31, RL)
 */
struct VibrationalEnergyLevels {
  using Key = QuantumIdentifier;
  using T = Numeric;
  using size_type = std::size_t;

  std::unordered_map<Key, T, Quantum::Number::GlobalStateHash> data;

  [[nodiscard]] auto begin() {return data.begin();}
  [[nodiscard]] auto begin() const {return data.begin();}
  [[nodiscard]] auto cbegin() const {return data.cbegin();}
  [[nodiscard]] auto end() {return data.end();}
  [[nodiscard]] auto end() const {return data.end();}
  [[nodiscard]] auto cend() const {return data.cend();}
  [[nodiscard]] T& at( const Key& key ) {return data.at(key);}
  [[nodiscard]] const T& at( const Key& key ) const {return data.at(key);}
  T& operator[]( const Key& key ) {return data[key];}
  T& operator[]( Key&& key ) {return data[std::move(key)];}
  [[nodiscard]] bool empty() const noexcept {return data.empty();}
  [[nodiscard]] size_type size() const noexcept {return data.size(); }
  [[nodiscard]] Index nelem() const noexcept {return static_cast<Index>(size());}

  friend std::ostream& operator<<(std::ostream& os, const VibrationalEnergyLevels& vib);
};

/** Sets up the solution matrix for linear statistical equilibrium equation
 * 
 * @param[in,out] A Matrix to solve SEE
 * @param[in] Aij Einstein coefficient for spontaneuos emission of all lines
 * @param[in] Bij Einstein coefficient for induced emission of all lines
 * @param[in] Bji Einstein coefficient for induced absorption of all lines
 * @param[in] Cij Collisional rate of change from upper to lower state level
 * @param[in] Cji Collisional rate of change from lower to upper state level
 * @param[in] Jij Radiation field for the upper to lower transition for each line
 * @param[in] upper Index list for upper state levels for each line
 * @param[in] lower Index list for lower state levels for each line
 */
void statistical_equilibrium_equation(MatrixView A,
                                      const ConstVectorView& Aij,
                                      const ConstVectorView& Bij,
                                      const ConstVectorView& Bji,
                                      const ConstVectorView& Cij,
                                      const ConstVectorView& Cji,
                                      const ConstVectorView& Jij,
                                      const ArrayOfIndex& upper,
                                      const ArrayOfIndex& lower);

/** Sets up the solution matrix for linear dampened statistical equilibrium equation
 * 
 * @param[in,out] A Matrix to solve SEE
 * @param[in] x Ratio of molecules for each state
 * @param[in] Aij Einstein coefficient for spontaneuos emission of all lines
 * @param[in] Bij Einstein coefficient for induced emission of all lines
 * @param[in] Bji Einstein coefficient for induced absorption of all lines
 * @param[in] Cij Collisional rate of change from upper to lower state level
 * @param[in] Cji Collisional rate of change from lower to upper state level
 * @param[in] Jij Radiation field for the upper to lower transition for each line
 * @param[in] Lambda Transmission for the upper to lower transition for each line
 * @param[in] upper Index list for upper state levels for each line
 * @param[in] lower Index list for lower state levels for each line
 */
void dampened_statistical_equilibrium_equation(
    MatrixView A,
    const ConstVectorView& x,
    const ConstVectorView& Aij,
    const ConstVectorView& Bij,
    const ConstVectorView& Bji,
    const ConstVectorView& Cij,
    const ConstVectorView& Cji,
    const ConstVectorView& Jij,
    const ConstVectorView& Lambda,
    const ArrayOfIndex& upper,
    const ArrayOfIndex& lower,
    const Numeric& total_number_count = 1.0);

/** Set a row of the SEE matrix and level distribution vector to constant
 * 
 * @param[in] A Matrix to solve SEE
 * @param[in] x Ratio of molecules for each state
 * @param[in] sem_ratio Largest ratio possible
 * @param[in] row Row of SEE Matrix being zero
 */
void set_constant_statistical_equilibrium_matrix(MatrixView A,
                                                 VectorView x,
                                                 const Numeric& sem_ratio,
                                                 const Index row);

/** Create a Aij object
 * 
 * @param[in] abs_lines All lines of interest
 * @return Vector Einstein coefficient for spontaneuos emission of all lines
 */
Vector createAij(const ArrayOfArrayOfAbsorptionLines& abs_lines);

/** Create a Bij object
 * 
 * @param[in] abs_lines All lines of interest
 * @return Vector Einstein coefficient for induced emission of all lines
 */
Vector createBij(const ArrayOfArrayOfAbsorptionLines& abs_lines);

/** Create a Bji object
 * 
 * @param[in] Bij Einstein coefficient for induced emission of all lines
 * @param[in] abs_lines All lines of interest
 * @return Vector Einstein coefficient for induced absorption of all lines
 */
Vector createBji(const Vector& Bij, const ArrayOfArrayOfAbsorptionLines& abs_lines);

/** Create a Cji object
 * 
 * @param[in] Cij Collisional rate of change from upper to lower state level
 * @param[in] abs_lines All lines of interest
 * @param[in] T Temperature
 * @return Vector Collisional rate of change from lower to upper state level
 */
Vector createCji(const Vector& Cij,
                 const ArrayOfArrayOfAbsorptionLines& abs_lines,
                 const Numeric& T);

/** Set the Cji object
 * 
 * @param[in,out] Cji Collisional rate of change from lower to upper state level
 * @param[in] Cij Collisional rate of change from upper to lower state level
 * @param[in] abs_lines All lines of interest
 * @param[in] T Temperature
 * @param[in] n Size of Cij
 */
void setCji(Vector& Cji,
            const Vector& Cij,
            const ArrayOfArrayOfAbsorptionLines& abs_lines,
            const Numeric& T);

/** Gets collisional factors from coefficients
 * 
 * @param[in,out] Cij Collisional rate of change from upper to lower state level
 * @param[in,out] Cji Collisional rate of change from lower to upper state level
 * @param[in] abs_lines All lines of interest
 * @param[in] abs_species All absorption species
 * @param[in] collision_coefficients As WSV
 * @param[in] collision_line_identifiers As WSV
 * @param[in] isotopologue_ratios As WSV
 * @param[in] vmr Volume mixing ratios of absoprtion species
 * @param[in] T Temperature
 * @param[in] P Pressure
 */
void nlte_collision_factorsCalcFromCoeffs(
  Vector& Cij,
  Vector& Cji,
  const ArrayOfArrayOfAbsorptionLines& abs_lines,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const ArrayOfArrayOfGriddedField1& collision_coefficients,
  const ArrayOfQuantumIdentifier& collision_line_identifiers,
  const SpeciesIsotopologueRatios& isotopologue_ratios,
  const ConstVectorView& vmr,
  const Numeric& T,
  const Numeric& P);

/** Finds upper and lower states in SEE Matrix
 * 
 * @param[in,out] upper Index list for upper state levels for each line
 * @param[in,out] lower Index list for lower state levels for each line
 * @param[in] abs_lines All lines of interest
 * @param[in] nlte_quantum_identifiers As WSV
 */
void nlte_positions_in_statistical_equilibrium_matrix(
  ArrayOfIndex& upper,
  ArrayOfIndex& lower,
  const ArrayOfArrayOfAbsorptionLines& abs_lines,
  const EnergyLevelMap& nlte_field);

/** Finds a unique lower state if one exists or returns index to last element
 * 
 * @param[in] upper Index list for upper state levels for each line
 * @param[in] lower Index list for lower state levels for each line
 * @return Index Pos of unique element or len-1
 */
Index find_first_unique_in_lower(const ArrayOfIndex& upper,
                                 const ArrayOfIndex& lower) ARTS_NOEXCEPT;

/** Checks that a WSV is OK or throws a run-time error
 * 
 * @param[in] collision_line_identifiers As WSV
 */
void check_collision_line_identifiers(
    const ArrayOfQuantumIdentifier& collision_line_identifiers);

#endif /* NLTE_H */

