/**
 * @file radiation_field.h
 * @author Richard Larsson
 * @date 2019-09-04
 * 
 * @brief Radiation field calculations
 */

#ifndef radiation_field_h
#define radiation_field_h

#include "matpack_data.h"
#include "mystring.h"

#include <rtepack.h>

/** Throws an error if integration values are bad
 * 
 * @param[in] error_msg Error message to print
 * @param[in] value_that_should_be_unity Check value
 */
void error_in_integrate(const String& error_msg,
                        const Numeric& value_that_should_be_unity);

/** Integrate the line shape
 * 
 * f must be sorted and in same order as F
 * 
 * Return should be 1.0 for a full line shape
 * 
 * @param[in] F Line shape
 * @param[in] f Frequency grod
 * @return Numeric Integrated line shape
 */
Numeric test_integrate_convolved(const Eigen::Ref<Eigen::VectorXcd>& F,
                                 const Vector& f);

/** Integrate cos(za) over the angles
 * 
 * sorted_index must be in same order as cosza
 * 
 * Return should be 1.0 for a full set of zenith angles
 * 
 * @param[in] cosza cos of zenith angle
 * @param[in] sorted_index Order of zenith angles
 * @return Numeric Integration
 */
Numeric test_integrate_zenith(const Vector& cosza,
                              const Array<Index>& sorted_index);

/** Convolve intensity and line shape and integrate
 * 
 * f must be sorted and in same order as F and I
 * 
 * F must be normalized
 * 
 * @param[in] I Intensity vector
 * @param[in] F Line shape normalized
 * @param[in] f Frequency grid
 * @return Numeric Integrated absorption
 */
Numeric integrate_convolved(const StokvecVector& I,
                            const Eigen::VectorXcd& F,
                            const Vector& f);

/** Convolve transmission and line shape and integrate
 * 
 * f must be sorted and in same order as F and T
 * 
 * F must be normalized
 * 
 * Only consider [0, 0] position of T
 * 
 * @param[in] T Transmission matrix
 * @param[in] F Line shape normalized
 * @param[in] f Frequency grid
 * @return Numeric Integrated absorption
 */
Numeric integrate_convolved(const MuelmatVector& T,
                            const Eigen::VectorXcd& F,
                            const Vector& f);

/** Convolve source function with 1D sphere and integrate
 * 
 * @param[in] j Source vector
 * @param[in] cosza cos of zenith angle
 * @param[in] sorted_index Order of zenith angles
 * @return Numeric Integrated source function
 */
Numeric integrate_zenith(const ConstVectorView& j,
                         const Vector& cosza,
                         const Array<Index>& sorted_index);
#endif  // radiation_field_h
