#pragma once

#include "grids.h"
#include "matpack_data.h"
namespace Legendre {

//! Stores the up (U), south (S), east (E) values of a field relative to the sphere
struct SphericalField {
  Numeric U{0};
  Numeric S{0};
  Numeric E{0};
  
  //! Returns the total strength of the field
  [[nodiscard]] Numeric total() const noexcept {return std::hypot(U, S, E);}
  
  //! Returns the total strength of the field
  [[nodiscard]] Numeric total_horizontal() const noexcept {return std::hypot(S, E);}
  
  //! Always construct to zeroes explicitly
  constexpr SphericalField() noexcept = default;
};


//! Holds a SphericalField for multiple dimensions (radius times longitudes)
using MatrixOfSphericalField = Grid<SphericalField, 2>;


/** Computes the spherical field
 * 
 * If latitude is beyond a limit, the longitude is set to zero
 * and the latitude at the limit away from the pole is used. For these
 * position, the absolute of the horizontal component is kept around in
 * the south-facing component and the east-facing component is set to zero
 * 
 * The latitude limit is defined in the ColatitudeConversion struct.
 * 
 * The longitude is constrained to the range [-180, 180) to have consistent
 * behavior at the daytime border.
 * 
 * This function is taken from the implementation by Isabela de Oliveira Martins
 * at https://github.com/de-oliveira/IsabelaFunctions/blob/master/IsabelaFunctions/fieldmodel.py
 * (2021-05-06).
 * 
 * @param[in] g A N x N matrix of g-coefficients
 * @param[in] h A N x N matrix of h-coefficients
 * @param[in] r0 The reference radius (spherical)
 * @param[in] r The actual radius (spherical)
 * @param[in] lat The latitude (spherical)
 * @param[in] lon The longitude (spherical)
 * @return A spherical field
 */
SphericalField schmidt_fieldcalc(const Matrix& g, const Matrix& h, const Numeric r0, const Numeric r, const Numeric lat, const Numeric lon) ARTS_NOEXCEPT;


/** Computes the spherical field for many radius and longitudes
 * 
 * See purely Numeric implementation for more information.
 * 
 * @param[in] g A N x N matrix of g-coefficients
 * @param[in] h A N x N matrix of h-coefficients
 * @param[in] r0 The reference radius (spherical)
 * @param[in] r The actual radius (spherical)
 * @param[in] lat The latitude (spherical)
 * @param[in] lon The longitude (spherical)
 * @return A spherical field
 */
MatrixOfSphericalField schmidt_fieldcalc(const Matrix& g, const Matrix& h, const Numeric r0, const Vector& r, const Numeric lat, const Vector& lon) ARTS_NOEXCEPT;
} // namespace Legendre
