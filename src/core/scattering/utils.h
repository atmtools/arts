#ifndef UTILS_H_
#define UTILS_H_
#include <iostream>
#include <vector>
#include <algorithm>

#include "interpolation.h"


namespace scattering {
GridPos find_interp_weights(StridedConstVectorView grid, Numeric x_new);
Index digitize(const Vector& boundaries, Numeric value);

}

#endif // UTILS_H_
