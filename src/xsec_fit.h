/*!
  \file   xsec_fit.h
  \author Oliver Lemke <oliver.lemke@uni-hamburg.de>
  \date   2018-01-08

  \brief  Methods and classes for HITRAN absorption cross section data.
*/

#ifndef HITRAN_XSEC_H
#define HITRAN_XSEC_H

#include "array.h"
#include "arts.h"
#include "bifstream.h"
#include "gridded_fields.h"
#include "matpack_data.h"
#include "messages.h"
#include "mystring.h"
#include "species.h"

#include <memory>

/** Hitran crosssection class.
 *
 * Stores the coefficients from our model for hitran crosssection data and
 * applies them to calculate the crossections.
 */
class XsecRecord {
 public:
  /** Return species index */
  [[nodiscard]] const Species::Species& Species() const { return mspecies; };

  /** Return species name */
  [[nodiscard]] String SpeciesName() const;

  /** Set species name */
  void SetSpecies(const Species::Species species) { mspecies = species; };

  /** Return species index */
  [[nodiscard]] Index Version() const { return mversion; };

  /** Set species name */
  void SetVersion(Index version);

  /** Calculate hitran cross section data.

     Calculate crosssections at each frequency for given pressure and
     temperature.

     \param[out] result     Crosssections for given frequency grid.
     \param[in] f_grid      Frequency grid.
     \param[in] pressure    Scalar pressure.
     \param[in] temperature Scalar temperature.
     \param[in] verbosity   Verbosity.
     */
  void Extract(VectorView result,
               const Vector& f_grid,
               Numeric pressure,
               Numeric temperature,
               const Verbosity& verbosity) const;

  /************ VERSION 2 *************/
  /** Get mininum pressures from fit */
  [[nodiscard]] const Vector& FitMinPressures() const {
    return mfitminpressures;
  };

  /** Get maximum pressures from fit */
  [[nodiscard]] const Vector& FitMaxPressures() const {
    return mfitmaxpressures;
  };

  /** Get mininum temperatures from fit */
  [[nodiscard]] const Vector& FitMinTemperatures() const {
    return mfitmintemperatures;
  };

  /** Get maximum temperatures */
  [[nodiscard]] const Vector& FitMaxTemperatures() const {
    return mfitmaxtemperatures;
  };

  /** Get coefficients */
  [[nodiscard]] const ArrayOfGriddedField2& FitCoeffs() const {
    return mfitcoeffs;
  };

  /** Get mininum pressures from fit */
  [[nodiscard]] Vector& FitMinPressures() { return mfitminpressures; };

  /** Get maximum pressures from fit */
  [[nodiscard]] Vector& FitMaxPressures() { return mfitmaxpressures; };

  /** Get mininum temperatures from fit */
  [[nodiscard]] Vector& FitMinTemperatures() { return mfitmintemperatures; };

  /** Get maximum temperatures */
  [[nodiscard]] Vector& FitMaxTemperatures() { return mfitmaxtemperatures; };

  /** Get coefficients */
  [[nodiscard]] ArrayOfGriddedField2& FitCoeffs() { return mfitcoeffs; };

  friend std::ostream& operator<<(std::ostream& os, const XsecRecord& xd);

 private:
  /** Calculate crosssections */
  void CalcXsec(VectorView xsec,
                const Index dataset,
                const Numeric pressure,
                const Numeric temperature) const;

  // /** Calculate temperature derivative of crosssections */
  // void CalcDT(VectorView xsec_dt,
  //             Index dataset,
  //             Numeric temperature) const;

  // /** Calculate pressure derivative of crosssections */
  // void CalcDP(VectorView xsec_dp,
  //             Index dataset,
  //             Numeric pressure) const;

  static constexpr Index P00 = 0;
  static constexpr Index P10 = 1;
  static constexpr Index P01 = 2;
  static constexpr Index P20 = 3;

  static constexpr Index mversion{2};
  Species::Species mspecies;
  /* VERSION 2 */
  Vector mfitminpressures;
  Vector mfitmaxpressures;
  Vector mfitmintemperatures;
  Vector mfitmaxtemperatures;
  ArrayOfGriddedField2 mfitcoeffs;
};

using ArrayOfXsecRecord = Array<XsecRecord>;

Index hitran_xsec_get_index(const ArrayOfXsecRecord& xsec_data,
                            Species::Species species);

std::shared_ptr<XsecRecord> hitran_xsec_get_data(
    const std::vector<std::shared_ptr<XsecRecord>>& xsec_data,
    const Species::Species species);

#endif  // HITRAN_XSEC_H
