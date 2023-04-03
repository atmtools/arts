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
#include "rte.h"

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
Numeric integrate_convolved(const RadiationVector& I,
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
Numeric integrate_convolved(const TransmissionMatrix& T,
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

/** Get a discrete position from grid pos
 * 
 * Assumes the grid pos has been set to 
 * extend over the path
 * 
 * @param[in] gp Grid position
 * @return Index Position
 */
Index grid_index_from_gp(const GridPos& gp);

/** Get sorting of zenith angles in field of ppath
 * 
 * @param[out] sorted_index Order of zenith angles
 * @param[out] cosza cos of zenith angle
 * @param[in] ppath_field As WSV
 */
void sorted_index_of_ppath_field(ArrayOfArrayOfIndex& sorted_index,
                                 ArrayOfVector& cosza,
                                 const ArrayOfPpath& ppath_field);

#endif  // radiation_field_h
