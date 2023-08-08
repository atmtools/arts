#include "arts_constants.h"
#include "matpack_data.h"

constexpr Numeric SPEED_OF_LIGHT = Constant::speed_of_light;
constexpr Numeric PI = Constant::pi;

/* Workspace method: Doxygen documentation will be auto-generated */
void FrequencyFromWavelength(  // WS Generic Output
    Numeric& frequency,
    // WS Generic Input
    const Numeric& wavelength) {
  // Convert from wavelength to frequency
  frequency = SPEED_OF_LIGHT / wavelength;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void FrequencyFromWavelength(  // WS Generic Output
    Vector& frequency,
    // WS Generic Input
    const Vector& wavelength) {
  frequency.resize(wavelength.nelem());
  // Convert from wavelength to frequency
  for (Index i = 0; i < wavelength.nelem(); i++)
    frequency[i] = SPEED_OF_LIGHT / wavelength[i];
}

/* Workspace method: Doxygen documentation will be auto-generated */
void FrequencyFromCGSAngularWavenumber(  // WS Generic Output
    Numeric& frequency,
    // WS Generic Input
    const Numeric& angular_wavenumber) {
  frequency = SPEED_OF_LIGHT * angular_wavenumber / (2 * PI) * 100;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void FrequencyFromCGSAngularWavenumber(  // WS Generic Output
    Vector& frequency,
    // WS Generic Input
    const Vector& angular_wavenumber) {
  frequency.resize(angular_wavenumber.nelem());
  // Convert from angular wavenumber to frequency
  for (Index i = 0; i < angular_wavenumber.nelem(); i++)
    frequency[i] = SPEED_OF_LIGHT * angular_wavenumber[i] / (2 * PI) * 100;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void FrequencyFromCGSKayserWavenumber(  // WS Generic Output
    Numeric& frequency,
    // WS Generic Input
    const Numeric& kayser_wavenumber) {
  frequency = SPEED_OF_LIGHT * kayser_wavenumber * 100;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void FrequencyFromCGSKayserWavenumber(  // WS Generic Output
    Vector& frequency,
    // WS Generic Input
    const Vector& kayser_wavenumber) {
  frequency.resize(kayser_wavenumber.nelem());
  // Convert from Kayser wavenumber to frequency
  for (Index i = 0; i < kayser_wavenumber.nelem(); i++)
    frequency[i] = SPEED_OF_LIGHT * kayser_wavenumber[i] * 100;
}
