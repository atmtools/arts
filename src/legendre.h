#pragma once

#include "grids.h"
#include "matpack_data.h"
namespace Legendre {
using Vector3 = std::array<Numeric, 3>;

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
 * @param[in] pos The position [r, lat, lon] (spherical)
 * @return A spherical field {Br, Btheta, Bphi}
 */
Vector3 schmidt_fieldcalc(const Matrix& g,
                          const Matrix& h,
                          const Numeric r0,
                          const Vector3 pos);
} // namespace Legendre
