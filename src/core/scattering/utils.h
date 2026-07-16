#ifndef UTILS_H_
#define UTILS_H_
#include <interpolation.h>

#include <algorithm>
#include <iostream>
#include <vector>

namespace scattering {
GridPos find_interp_weights(StridedConstVectorView grid, Numeric x_new);
Index   digitize(const Vector& boundaries, Numeric value);

}  // namespace scattering

#endif  // UTILS_H_
