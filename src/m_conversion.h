/*!
  \file   m_conversion.h
  \author Claudia Emde <claudia.emde@lmu.de>
  \date   2010-07-21
  
  \brief  Implementation of unit conversion functions
  
*/

#ifndef m_conversion_h
#define m_conversion_h

#include "arts_constants.h"
#include "matpack_data.h"
#include "messages.h"

inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;
inline constexpr Numeric PI=Constant::pi;

/* Workspace method: Doxygen documentation will be auto-generated */
inline void FrequencyFromWavelength(  // WS Generic Output
    Numeric& frequency,
    // WS Generic Input
    const Numeric& wavelength,
    const Verbosity&) {
  // Convert from wavelength to frequency
  frequency = SPEED_OF_LIGHT / wavelength;
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void FrequencyFromWavelength(  // WS Generic Output
    Vector& frequency,
    // WS Generic Input
    const Vector& wavelength,
    const Verbosity&) {
  frequency.resize(wavelength.nelem());
  // Convert from wavelength to frequency
  for (Index i = 0; i < wavelength.nelem(); i++)
    frequency[i] = SPEED_OF_LIGHT / wavelength[i];
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void FrequencyFromCGSAngularWavenumber(  // WS Generic Output
    Numeric& frequency,
    // WS Generic Input
    const Numeric& angular_wavenumber,
    const Verbosity&) {
  frequency = SPEED_OF_LIGHT * angular_wavenumber / (2 * PI) * 100;
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void FrequencyFromCGSAngularWavenumber(  // WS Generic Output
    Vector& frequency,
    // WS Generic Input
    const Vector& angular_wavenumber,
    const Verbosity&) {
  frequency.resize(angular_wavenumber.nelem());
  // Convert from angular wavenumber to frequency
  for (Index i = 0; i < angular_wavenumber.nelem(); i++)
    frequency[i] = SPEED_OF_LIGHT * angular_wavenumber[i] / (2 * PI) * 100;
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void FrequencyFromCGSKayserWavenumber(  // WS Generic Output
    Numeric& frequency,
    // WS Generic Input
    const Numeric& kayser_wavenumber,
    const Verbosity&) {
  frequency = SPEED_OF_LIGHT * kayser_wavenumber * 100;
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void FrequencyFromCGSKayserWavenumber(  // WS Generic Output
    Vector& frequency,
    // WS Generic Input
    const Vector& kayser_wavenumber,
    const Verbosity&) {
  frequency.resize(kayser_wavenumber.nelem());
  // Convert from Kayser wavenumber to frequency
  for (Index i = 0; i < kayser_wavenumber.nelem(); i++)
    frequency[i] = SPEED_OF_LIGHT * kayser_wavenumber[i] * 100;
}

#endif /* m_conversion_h */
