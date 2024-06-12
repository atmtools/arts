#include "igrf13.h"

#include <algorithm>

#include "arts_conversions.h"
#include "arts_omp.h"
#include "debug.h"
#include "legendre.h"
#include "matpack_constexpr.h"

/** International Geomagnetic Reference Field version 13
 *
 * The data from the years 2000-2020 inclusive is available below
 *
 * The IGRF model provides model data for \f$g\f$ and \f$h\f$ in the spherical
 * harmonics
 *
 * \f[
 * V(r, \theta, \phi, t) = a
 * \sum_{n=1}^N\sum_{m=0}^n\left(\frac{a}{r}\right)^{n+1}\left[g_n^m\left(t\right)\cos\left(m\phi\right)
 * + h_n^m\left(t\right)\sin\left(m\phi\right)\right]
 * P_n^m\left(\cos\theta\right), \f]
 *
 * where \f$r\f$ is the radius, \f$\theta\f$ is the colatitude, \f$\phi\f$ is
 * the longitude, \f$t\f$ is the time, \f$a\f$ is the reference radius at 6371.2
 * km, and \f$P_n^m\left(\cos\theta\right)\f$ are the Schmidth normalized
 * Legendre polynominal.
 *
 * The magnetic field itself is computed from the gradients of \f$ V(r, \theta,
 * \phi, t) \f$
 */
namespace IGRF {
//! g-coefficients for 2020 (14x14 matrix)
constexpr matpack::matpack_constant_data<Numeric, 14, 14> g2020{
    0.0,   0.0,     0.0,    0.0,    0.0,  0.0,      0.0,     0.0,     0.0,
    0.0,   0.0,     0.0,    0.0,    0.0,  -29404.8, -1450.9, 0.0,     0.0,
    0.0,   0.0,     0.0,    0.0,    0.0,  0.0,      0.0,     0.0,     0.0,
    0.0,   -2499.6, 2982.0, 1677.0, 0.0,  0.0,      0.0,     0.0,     0.0,
    0.0,   0.0,     0.0,    0.0,    0.0,  0.0,      1363.2,  -2381.2, 1236.2,
    525.7, 0.0,     0.0,    0.0,    0.0,  0.0,      0.0,     0.0,     0.0,
    0.0,   0.0,     903.0,  809.5,  86.3, -309.4,   48.0,    0.0,     0.0,
    0.0,   0.0,     0.0,    0.0,    0.0,  0.0,      0.0,     -234.3,  363.2,
    187.8, -140.7,  -151.2, 13.5,   0.0,  0.0,      0.0,     0.0,     0.0,
    0.0,   0.0,     0.0,    66.0,   65.5, 72.9,     -121.5,  -36.2,   13.5,
    -64.7, 0.0,     0.0,    0.0,    0.0,  0.0,      0.0,     0.0,     80.6,
    -76.7, -8.2,    56.5,   15.8,   6.4,  -7.2,     9.8,     0.0,     0.0,
    0.0,   0.0,     0.0,    0.0,    23.7, 9.7,      -17.6,   -0.5,    -21.1,
    15.3,  13.7,    -16.5,  -0.3,   0.0,  0.0,      0.0,     0.0,     0.0,
    5.0,   8.4,     2.9,    -1.5,   -1.1, -13.2,    1.1,     8.8,     -9.3,
    -11.9, 0.0,     0.0,    0.0,    0.0,  -1.9,     -6.2,    -0.1,    1.7,
    -0.9,  0.7,     -0.9,   1.9,    1.4,  -2.4,     -3.8,    0.0,     0.0,
    0.0,   3.0,     -1.4,   -2.5,   2.3,  -0.9,     0.3,     -0.7,    -0.1,
    1.4,   -0.6,    0.2,    3.1,    0.0,  0.0,      -2.0,    -0.1,    0.5,
    1.3,   -1.2,    0.7,    0.3,    0.5,  -0.3,     -0.5,    0.1,     -1.1,
    -0.3,  0.0,     0.1,    -0.9,   0.5,  0.7,      -0.3,    0.8,     0.0,
    0.8,   0.0,     0.4,    0.1,    0.5,  -0.5,     -0.4,
};

//! h-coefficients for 2020 (14x14 matrix)
constexpr matpack::matpack_constant_data<Numeric, 14, 14> h2020{
    0.0,    0.0,  0.0,   0.0,    0.0,   0.0,    0.0,   0.0,   0.0,    0.0,
    0.0,    0.0,  0.0,   0.0,    0.0,   4652.5, 0.0,   0.0,   0.0,    0.0,
    0.0,    0.0,  0.0,   0.0,    0.0,   0.0,    0.0,   0.0,   0.0,    -2991.6,
    -734.6, 0.0,  0.0,   0.0,    0.0,   0.0,    0.0,   0.0,   0.0,    0.0,
    0.0,    0.0,  0.0,   -82.1,  241.9, -543.4, 0.0,   0.0,   0.0,    0.0,
    0.0,    0.0,  0.0,   0.0,    0.0,   0.0,    0.0,   281.9, -158.4, 199.7,
    -349.7, 0.0,  0.0,   0.0,    0.0,   0.0,    0.0,   0.0,   0.0,    0.0,
    0.0,    47.7, 208.3, -121.2, 32.3,  98.9,   0.0,   0.0,   0.0,    0.0,
    0.0,    0.0,  0.0,   0.0,    0.0,   -19.1,  25.1,  52.8,  -64.5,  8.9,
    68.1,   0.0,  0.0,   0.0,    0.0,   0.0,    0.0,   0.0,   0.0,    -51.5,
    -16.9,  2.2,  23.5,  -2.2,   -27.2, -1.8,   0.0,   0.0,   0.0,    0.0,
    0.0,    0.0,  0.0,   8.4,    -15.3, 12.8,   -11.7, 14.9,  3.6,    -6.9,
    2.8,    0.0,  0.0,   0.0,    0.0,   0.0,    0.0,   -23.4, 11.0,   9.8,
    -5.1,   -6.3, 7.8,   0.4,    -1.4,  9.6,    0.0,   0.0,   0.0,    0.0,
    0.0,    3.4,  -0.2,  3.6,    4.8,   -8.6,   -0.1,  -4.3,  -3.4,   -0.1,
    -8.8,   0.0,  0.0,   0.0,    0.0,   0.0,    2.5,   -0.6,  -0.4,   0.6,
    -0.2,   -1.7, -1.6,  -3.0,   -2.0,  -2.6,   0.0,   0.0,   0.0,    -1.2,
    0.5,    1.4,  -1.8,  0.1,    0.8,   -0.2,   0.6,   0.2,   -0.9,   0.0,
    0.5,    0.0,  0.0,   -0.9,   0.6,   1.4,    -0.4,  -1.3,  -0.1,   0.3,
    -0.1,   0.5,  0.5,   -0.4,   -0.4,  -0.6,
};

//! g-coefficients for 2015 (14x14 matrix)
constexpr matpack::matpack_constant_data<Numeric, 14, 14> g2015{
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,    0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,    -29441.46, -1501.77,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,    0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      -2445.88, 3012.2, 1676.35,   0.0,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,    0.0,       0.0,
    0.0,    0.0,     1350.33, -2352.26, 1225.85,  581.69, 0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,    0.0,       0.0,
    907.42, 813.68,  120.49,  -334.85,  70.38,    0.0,    0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,    -232.91,   360.14,
    192.35, -140.94, -157.4,  4.3,      0.0,      0.0,    0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      69.55,    67.57,  72.79,     -129.85,
    -28.93, 13.14,   -70.85,  0.0,      0.0,      0.0,    0.0,       0.0,
    0.0,    0.0,     81.29,   -75.99,   -6.79,    51.82,  15.07,     9.32,
    -2.88,  6.61,    0.0,     0.0,      0.0,      0.0,    0.0,       0.0,
    23.98,  8.89,    -16.78,  -3.16,    -20.56,   13.33,  11.76,     -15.98,
    -2.02,  0.0,     0.0,     0.0,      0.0,      0.0,    5.33,      8.83,
    3.02,   -3.22,   0.67,    -13.2,    -0.1,     8.68,   -9.06,     -10.54,
    0.0,    0.0,     0.0,     0.0,      -2.01,    -6.26,  0.17,      0.55,
    -0.55,  1.7,     -0.67,   2.13,     2.33,     -1.8,   -3.59,     0.0,
    0.0,    0.0,     3.0,     -1.4,     -2.3,     2.08,   -0.79,     0.58,
    -0.7,   0.14,    1.7,     -0.22,    0.44,     3.49,   0.0,       0.0,
    -2.09,  -0.16,   0.46,    1.23,     -0.89,    0.85,   0.1,       0.54,
    -0.37,  -0.43,   0.22,    -0.94,    -0.03,    0.0,    -0.02,     -0.92,
    0.42,   0.63,    -0.42,   0.96,     -0.19,    0.81,   -0.13,     0.38,
    0.08,   0.46,    -0.35,   -0.36,
};

//! h-coefficients for 2015 (14x14 matrix)
constexpr matpack::matpack_constant_data<Numeric, 14, 14> h2015{
    0.0,    0.0,     0.0,     0.0,     0.0,     0.0,      0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,     0.0,     0.0,      0.0,     4795.99,
    0.0,    0.0,     0.0,     0.0,     0.0,     0.0,      0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,     0.0,     -2845.41, -642.17, 0.0,
    0.0,    0.0,     0.0,     0.0,     0.0,     0.0,      0.0,     0.0,
    0.0,    0.0,     0.0,     -115.29, 245.04,  -538.7,   0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,     0.0,     0.0,      0.0,     0.0,
    0.0,    283.54,  -188.43, 180.95,  -329.23, 0.0,      0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,     0.0,     0.0,      0.0,     46.98,
    196.98, -119.14, 15.98,   100.12,  0.0,     0.0,      0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,     0.0,     -20.61,   33.3,    58.74,
    -66.64, 7.35,    62.41,   0.0,     0.0,     0.0,      0.0,     0.0,
    0.0,    0.0,     0.0,     -54.27,  -19.53,  5.59,     24.45,   3.27,
    -27.5,  -2.32,   0.0,     0.0,     0.0,     0.0,      0.0,     0.0,
    0.0,    10.04,   -18.26,  13.18,   -14.6,   16.16,    5.69,    -9.1,
    2.26,   0.0,     0.0,     0.0,     0.0,     0.0,      0.0,     -21.77,
    10.76,  11.74,   -6.74,   -6.88,   7.79,    1.04,     -3.89,   8.44,
    0.0,    0.0,     0.0,     0.0,     0.0,     3.28,     -0.4,    4.55,
    4.4,    -7.92,   -0.61,   -4.16,   -2.85,   -1.12,    -8.72,   0.0,
    0.0,    0.0,     0.0,     0.0,     2.11,    -0.6,     -1.05,   0.76,
    -0.2,   -2.12,   -1.44,   -2.57,   -2.01,   -2.34,    0.0,     0.0,
    0.0,    -1.08,   0.37,    1.75,    -2.19,   0.27,     0.72,    -0.09,
    0.29,   0.23,    -0.89,   -0.16,   0.72,    0.0,      0.0,     -0.88,
    0.49,   1.56,    -0.5,    -1.24,   -0.1,    0.42,     -0.04,   0.48,
    0.48,   -0.3,    -0.43,   -0.71,
};

//! g-coefficients for 2010 (14x14 matrix)
constexpr matpack::matpack_constant_data<Numeric, 14, 14> g2010{
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,     -29496.57, -1586.42,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      -2396.06, 3026.34, 1668.17,   0.0,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,
    0.0,    0.0,     1339.85, -2326.54, 1232.1,   633.73,  0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,
    912.66, 808.97,  166.58,  -356.83,  89.4,     0.0,     0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,     -230.87,   357.29,
    200.26, -141.05, -163.17, -8.03,    0.0,      0.0,     0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      72.78,    68.69,   75.92,     -141.4,
    -22.83, 13.1,    -78.09,  0.0,      0.0,      0.0,     0.0,       0.0,
    0.0,    0.0,     80.44,   -75.0,    -4.55,    45.24,   14.0,      10.46,
    1.64,   4.92,    0.0,     0.0,      0.0,      0.0,     0.0,       0.0,
    24.41,  8.21,    -14.5,   -5.59,    -19.34,   11.61,   10.85,     -14.05,
    -3.54,  0.0,     0.0,     0.0,      0.0,      0.0,     5.5,       9.45,
    3.45,   -5.27,   3.13,    -12.38,   -0.76,    8.43,    -8.42,     -10.08,
    0.0,    0.0,     0.0,     0.0,      -1.94,    -6.24,   0.89,      -1.07,
    -0.16,  2.45,    -0.33,   2.13,     3.09,     -1.03,   -2.8,      0.0,
    0.0,    0.0,     3.05,    -1.48,    -2.03,    1.65,    -0.51,     0.54,
    -0.79,  0.37,    1.79,    0.12,     0.75,     3.75,    0.0,       0.0,
    -2.12,  -0.21,   0.3,     1.04,     -0.63,    0.95,    -0.11,     0.52,
    -0.39,  -0.37,   0.21,    -0.77,    0.04,     0.0,     -0.09,     -0.89,
    0.31,   0.42,    -0.45,   1.08,     -0.31,    0.78,    -0.18,     0.38,
    0.02,   0.42,    -0.26,   -0.26,
};

//! h-coefficients for 2010 (14x14 matrix)
constexpr matpack::matpack_constant_data<Numeric, 14, 14> h2010{
    0.0,    0.0,     0.0,     0.0,    0.0,     0.0,      0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,    0.0,     0.0,      0.0,     4944.26,
    0.0,    0.0,     0.0,     0.0,    0.0,     0.0,      0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,    0.0,     -2708.54, -575.73, 0.0,
    0.0,    0.0,     0.0,     0.0,    0.0,     0.0,      0.0,     0.0,
    0.0,    0.0,     0.0,     -160.4, 251.75,  -537.03,  0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,    0.0,     0.0,      0.0,     0.0,
    0.0,    286.48,  -211.03, 164.46, -309.72, 0.0,      0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,    0.0,     0.0,      0.0,     44.58,
    189.01, -118.06, -0.01,   101.04, 0.0,     0.0,      0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,    0.0,     -20.9,    44.18,   61.54,
    -66.26, 3.02,    55.4,    0.0,    0.0,     0.0,      0.0,     0.0,
    0.0,    0.0,     0.0,     -57.8,  -21.2,   6.54,     24.96,   7.03,
    -27.61, -3.28,   0.0,     0.0,    0.0,     0.0,      0.0,     0.0,
    0.0,    10.84,   -20.03,  11.83,  -17.41,  16.71,    6.96,    -10.74,
    1.64,   0.0,     0.0,     0.0,    0.0,     0.0,      0.0,     -20.54,
    11.51,  12.75,   -7.14,   -7.42,  7.97,    2.14,     -6.08,   7.01,
    0.0,    0.0,     0.0,     0.0,    0.0,     2.73,     -0.1,    4.71,
    4.44,   -7.22,   -0.96,   -3.95,  -1.99,   -1.97,    -8.31,   0.0,
    0.0,    0.0,     0.0,     0.13,   1.67,    -0.66,    -1.76,   0.85,
    -0.39,  -2.51,   -1.27,   -2.11,  -1.94,   -1.86,    0.0,     0.0,
    0.0,    -0.87,   0.27,    2.13,   -2.49,   0.49,     0.59,    0.0,
    0.13,   0.27,    -0.86,   -0.23,  0.87,    0.0,      0.0,     -0.87,
    0.3,    1.66,    -0.59,   -1.14,  -0.07,   0.54,     0.1,     0.49,
    0.44,   -0.25,   -0.53,   -0.79,
};

//! g-coefficients for 2005 (14x14 matrix)
constexpr matpack::matpack_constant_data<Numeric, 14, 14> g2005{
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,     -29554.63, -1669.05,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      -2337.24, 3047.69, 1657.76,   0.0,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,
    0.0,    0.0,     1336.3,  -2305.83, 1246.39,  672.51,  0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,     0.0,       0.0,
    920.55, 797.96,  210.65,  -379.86,  100.0,    0.0,     0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      0.0,      0.0,     -227.0,    354.41,
    208.95, -136.54, -168.05, -13.55,   0.0,      0.0,     0.0,       0.0,
    0.0,    0.0,     0.0,     0.0,      73.6,     69.56,   76.74,     -151.34,
    -14.58, 14.58,   -86.36,  0.0,      0.0,      0.0,     0.0,       0.0,
    0.0,    0.0,     79.88,   -74.46,   -1.65,    38.73,   12.3,      9.37,
    5.42,   1.94,    0.0,     0.0,      0.0,      0.0,     0.0,       0.0,
    24.8,   7.62,    -11.73,  -6.88,    -18.11,   10.17,   9.36,      -11.25,
    -4.87,  0.0,     0.0,     0.0,      0.0,      0.0,     5.58,      9.76,
    3.58,   -6.94,   5.01,    -10.76,   -1.25,    8.76,    -6.66,     -9.22,
    0.0,    0.0,     0.0,     0.0,      -2.17,    -6.12,   1.42,      -2.35,
    -0.15,  3.06,    0.29,    2.06,     3.77,     -0.21,   -2.09,     0.0,
    0.0,    0.0,     2.95,    -1.6,     -1.88,    1.44,    -0.31,     0.29,
    -0.79,  0.53,    1.8,     0.16,     0.96,     3.99,    0.0,       0.0,
    -2.15,  -0.29,   0.21,    0.89,     -0.38,    0.96,    -0.3,      0.46,
    -0.35,  -0.36,   0.08,    -0.49,    -0.08,    0.0,     -0.16,     -0.88,
    0.3,    0.28,    -0.43,   1.18,     -0.37,    0.75,    -0.26,     0.35,
    -0.05,  0.41,    -0.1,    -0.18,
};

//! h-coefficients for 2005 (14x14 matrix)
constexpr matpack::matpack_constant_data<Numeric, 14, 14> h2005{
    0.0,    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     5077.99,
    0.0,    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,     0.0,     -2594.5, -515.43, 0.0,
    0.0,    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
    0.0,    0.0,     0.0,     -198.86, 269.72,  -524.72, 0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
    0.0,    282.07,  -225.23, 145.15,  -305.36, 0.0,     0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     42.72,
    180.25, -123.45, -19.57,  103.85,  0.0,     0.0,     0.0,     0.0,
    0.0,    0.0,     0.0,     0.0,     0.0,     -20.33,  54.75,   63.63,
    -63.53, 0.24,    50.94,   0.0,     0.0,     0.0,     0.0,     0.0,
    0.0,    0.0,     0.0,     -61.14,  -22.57,  6.82,    25.35,   10.93,
    -26.32, -4.64,   0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
    0.0,    11.2,    -20.88,  9.83,    -19.71,  16.22,   7.61,    -12.76,
    -0.06,  0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     -20.11,
    12.69,  12.67,   -6.72,   -8.16,   8.1,     2.92,    -7.73,   6.01,
    0.0,    0.0,     0.0,     0.0,     0.0,     2.19,    0.1,     4.46,
    4.76,   -6.58,   -1.01,   -3.47,   -0.86,   -2.31,   -7.93,   0.0,
    0.0,    0.0,     0.0,     0.26,    1.44,    -0.77,   -2.27,   0.9,
    -0.58,  -2.69,   -1.08,   -1.58,   -1.9,    -1.39,   0.0,     0.0,
    0.0,    -0.55,   0.23,    2.38,    -2.63,   0.61,    0.4,     0.01,
    0.02,   0.28,    -0.87,   -0.34,   0.88,    0.0,     0.0,     -0.76,
    0.33,   1.72,    -0.54,   -1.07,   -0.04,   0.63,    0.21,    0.53,
    0.38,   -0.22,   -0.57,   -0.82,
};

//! g-coefficients for 2000 (14x14 matrix)
constexpr matpack::matpack_constant_data<Numeric, 14, 14> g2000{
    0.0,   0.0,     0.0,    0.0,    0.0,   0.0,      0.0,     0.0,     0.0,
    0.0,   0.0,     0.0,    0.0,    0.0,   -29619.4, -1728.2, 0.0,     0.0,
    0.0,   0.0,     0.0,    0.0,    0.0,   0.0,      0.0,     0.0,     0.0,
    0.0,   -2267.7, 3068.4, 1670.9, 0.0,   0.0,      0.0,     0.0,     0.0,
    0.0,   0.0,     0.0,    0.0,    0.0,   0.0,      1339.6,  -2288.0, 1252.1,
    714.5, 0.0,     0.0,    0.0,    0.0,   0.0,      0.0,     0.0,     0.0,
    0.0,   0.0,     932.3,  786.8,  250.0, -403.0,   111.3,   0.0,     0.0,
    0.0,   0.0,     0.0,    0.0,    0.0,   0.0,      0.0,     -218.8,  351.4,
    222.3, -130.4,  -168.6, -12.9,  0.0,   0.0,      0.0,     0.0,     0.0,
    0.0,   0.0,     0.0,    72.3,   68.2,  74.2,     -160.9,  -5.9,    16.9,
    -90.4, 0.0,     0.0,    0.0,    0.0,   0.0,      0.0,     0.0,     79.0,
    -74.0, 0.0,     33.3,   9.1,    6.9,   7.3,      -1.2,    0.0,     0.0,
    0.0,   0.0,     0.0,    0.0,    24.4,  6.6,      -9.2,    -7.9,    -16.6,
    9.1,   7.0,     -7.9,   -7.0,   0.0,   0.0,      0.0,     0.0,     0.0,
    5.0,   9.4,     3.0,    -8.4,   6.3,   -8.9,     -1.5,    9.3,     -4.3,
    -8.2,  0.0,     0.0,    0.0,    0.0,   -2.6,     -6.0,    1.7,     -3.1,
    -0.5,  3.7,     1.0,    2.0,    4.2,   0.3,      -1.1,    0.0,     0.0,
    0.0,   2.7,     -1.7,   -1.9,   1.5,   -0.1,     0.1,     -0.7,    0.7,
    1.7,   0.1,     1.2,    4.0,    0.0,   0.0,      -2.2,    -0.3,    0.2,
    0.9,   -0.2,    0.9,    -0.5,   0.3,   -0.3,     -0.4,    -0.1,    -0.2,
    -0.4,  0.0,     -0.2,   -0.9,   0.3,   0.1,      -0.4,    1.3,     -0.4,
    0.7,   -0.4,    0.3,    -0.1,   0.4,   0.0,      0.1,
};

//! h-coefficients for 2000 (14x14 matrix)
constexpr matpack::matpack_constant_data<Numeric, 14, 14> h2000{
    0.0,    0.0,  0.0,   0.0,    0.0,   0.0,    0.0,   0.0,   0.0,    0.0,
    0.0,    0.0,  0.0,   0.0,    0.0,   5186.1, 0.0,   0.0,   0.0,    0.0,
    0.0,    0.0,  0.0,   0.0,    0.0,   0.0,    0.0,   0.0,   0.0,    -2481.6,
    -458.0, 0.0,  0.0,   0.0,    0.0,   0.0,    0.0,   0.0,   0.0,    0.0,
    0.0,    0.0,  0.0,   -227.6, 293.4, -491.1, 0.0,   0.0,   0.0,    0.0,
    0.0,    0.0,  0.0,   0.0,    0.0,   0.0,    0.0,   272.6, -231.9, 119.8,
    -303.8, 0.0,  0.0,   0.0,    0.0,   0.0,    0.0,   0.0,   0.0,    0.0,
    0.0,    43.8, 171.9, -133.1, -39.3, 106.3,  0.0,   0.0,   0.0,    0.0,
    0.0,    0.0,  0.0,   0.0,    0.0,   -17.4,  63.7,  65.1,  -61.2,  0.7,
    43.8,   0.0,  0.0,   0.0,    0.0,   0.0,    0.0,   0.0,   0.0,    -64.6,
    -24.2,  6.2,  24.0,  14.8,   -25.4, -5.8,   0.0,   0.0,   0.0,    0.0,
    0.0,    0.0,  0.0,   11.9,   -21.5, 8.5,    -21.5, 15.5,  8.9,    -14.9,
    -2.1,   0.0,  0.0,   0.0,    0.0,   0.0,    0.0,   -19.7, 13.4,   12.5,
    -6.2,   -8.4, 8.4,   3.8,    -8.2,  4.8,    0.0,   0.0,   0.0,    0.0,
    0.0,    1.7,  0.0,   4.0,    4.9,   -5.9,   -1.2,  -2.9,  0.2,    -2.2,
    -7.4,   0.0,  0.0,   0.0,    0.0,   0.1,    1.3,   -0.9,  -2.6,   0.9,
    -0.7,   -2.8, -0.9,  -1.2,   -1.9,  -0.9,   0.0,   0.0,   0.0,    -0.4,
    0.3,    2.5,  -2.6,  0.7,    0.3,   0.0,    0.0,   0.3,   -0.9,   -0.4,
    0.8,    0.0,  0.0,   -0.9,   0.2,   1.8,    -0.4,  -1.0,  -0.1,   0.7,
    0.3,    0.6,  0.3,   -0.2,   -0.5,  -0.9,
};

//! The reference radius in IGRF13
constexpr Numeric r0{6371.2e3};

using Legendre::Vector3;
using Vector2 = std::array<Numeric, 2>;

Vector3 geocentric2ecef(const Vector3 pos) {
  const Numeric latrad = Conversion::deg2rad(pos[1]);
  const Numeric lonrad = Conversion::deg2rad(pos[2]);
  Vector3 ecef;
  ecef[0] = pos[0] * std::cos(latrad);  // Common term for x and z
  ecef[1] = ecef[0] * std::sin(lonrad);
  ecef[0] = ecef[0] * std::cos(lonrad);
  ecef[2] = pos[0] * std::sin(latrad);
  return ecef;
}

constexpr Numeric ellipsoid_radii_threshold = 1e-3;
constexpr bool is_ellipsoid_spherical(const Vector2 ellipsoid) {
  return (nonstd::abs(ellipsoid[0] - ellipsoid[1]) < ellipsoid_radii_threshold);
}

Vector3 geodetic2ecef(const Vector3 pos, const Vector2 refellipsoid) {
  Vector3 ecef;

  // Use geocentric function if geoid is spherical
  if (is_ellipsoid_spherical(refellipsoid)) {
    ecef = geocentric2ecef({pos[0] + refellipsoid[0], pos[1], pos[2]});
  } else {
    // See https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_geodetic_to_ECEF_coordinates
    const Numeric latrad = Conversion::deg2rad(pos[1]);
    const Numeric lonrad = Conversion::deg2rad(pos[2]);
    const Numeric sinlat = std::sin(latrad);
    const Numeric coslat = std::cos(latrad);
    const Numeric a2 = refellipsoid[0] * refellipsoid[0];
    const Numeric b2 = refellipsoid[1] * refellipsoid[1];
    const Numeric N =
        a2 / std::sqrt(a2 * coslat * coslat + b2 * sinlat * sinlat);
    const Numeric nhcos = (N + pos[0]) * coslat;
    ecef[0] = nhcos * std::cos(lonrad);
    ecef[1] = nhcos * std::sin(lonrad);
    ecef[2] = ((b2 / a2) * N + pos[0]) * sinlat;
  }

  return ecef;
}

Vector3 ecef2geocentric(const Vector3 ecef) {
  const Numeric r = std::hypot(ecef[0], ecef[1], ecef[2]);
  return {
      r, Conversion::asind(ecef[2] / r), Conversion::atan2d(ecef[1], ecef[0])};
}

Vector3 geodetic2geocentric(const Vector3 pos, const Vector2 ell) {
  return ecef2geocentric(geodetic2ecef(pos, ell));
}

/** Perform all computations on pre-allocated local data
 *
 * \param[in,out] out Size-initialized output data
 * \param[in] g The g-coefficients for the Legendre calculations
 * \param[in] h The h-coefficients for the Legendre calculations
 * \param[in] z_field As WSV
 * \param[in] lat_grid As WSV
 * \param[in] lon_grid As WSV
 * \param[in] ell As WSV called refellipsoid
 * \param[in] scale the addition scales with this, set to 1.0 for complete
 * calculations
 */
void compute_impl(MagneticField &out,
                  const Matrix &g,
                  const Matrix &h,
                  const Tensor3 &z_field,
                  const Vector &lat_grid,
                  const Vector &lon_grid,
                  const Vector &ell_,
                  const Numeric scale) ARTS_NOEXCEPT {
  using Conversion::cosd, Conversion::sind;

  // Constant sizes
  const Index nz = z_field.npages();
  const Index nlat = z_field.nrows();
  const Index nlon = z_field.ncols();

  const Vector2 ell{ell_[0], ell_[0] * std::sqrt(1.0 - ell_[1] * ell_[1])};

#pragma omp parallel for collapse(3) if (!arts_omp_in_parallel())
  for (Index iz = 0; iz < nz; iz++) {
    for (Index ilat = 0; ilat < nlat; ilat++) {
      for (Index ilon = 0; ilon < nlon; ilon++) {
        const Vector3 pos{
            z_field(iz, ilat, ilon), lat_grid[ilat], lon_grid[ilon]};
        const Vector3 geoc = geodetic2geocentric(pos, ell);
        const Vector3 mag = Legendre::schmidt_fieldcalc(g, h, r0, geoc);
        const Numeric ang = sind(lat_grid[ilat]) * sind(90.0 - geoc[1]) -
                            cosd(lat_grid[ilat]) * cosd(90.0 - geoc[1]);
        const Numeric ca = std::cos(ang);
        const Numeric sa = std::sin(ang);

        out.u(iz, ilat, ilon) += scale * mag[2];
        out.v(iz, ilat, ilon) += scale * (-ca * mag[1] - sa * mag[0]);
        out.w(iz, ilat, ilon) += scale * (-sa * mag[1] + ca * mag[0]);
      }
    }
  }
}

MagneticField compute(const Tensor3 &z_field,
                      const Vector &lat_grid,
                      const Vector &lon_grid,
                      const Time &time,
                      const Vector &ell) {
  ARTS_USER_ERROR_IF(
      ell.size() != 2, "ellipsoid must have two elements [a, e], is: ", ell)

  // Constant sizes
  const Index nz = z_field.npages();
  const Index nlat = z_field.nrows();
  const Index nlon = z_field.ncols();

  ARTS_USER_ERROR_IF(nlat != lat_grid.size(),
                     "z_field and lat_grid must have matching size: ",
                     nlat)
  ARTS_USER_ERROR_IF(nlon != lon_grid.size(),
                     "z_field and lon_grid must have matching size: ",
                     nlon)

  // Constant times
  static const Time y2020("2020-01-01 00:00:00");
  static const Time y2015("2015-01-01 00:00:00");
  static const Time y2010("2010-01-01 00:00:00");
  static const Time y2005("2005-01-01 00:00:00");
  static const Time y2000("2000-01-01 00:00:00");
  static const Matrix mg2020{g2020};
  static const Matrix mh2020{h2020};
  static const Matrix mg2015{g2015};
  static const Matrix mh2015{h2015};
  static const Matrix mg2010{g2010};
  static const Matrix mh2010{h2010};
  static const Matrix mg2005{g2005};
  static const Matrix mh2005{h2005};
  static const Matrix mg2000{g2000};
  static const Matrix mh2000{h2000};

  // Output
  MagneticField out(nz, nlat, nlon);  // Inits to zeroes

  // Select the correct time
  if (time >= y2020) {
    compute_impl(out, mg2020, mh2020, z_field, lat_grid, lon_grid, ell, 1.0);
  } else if (time >= y2015) {
    const Numeric scale = (time.Seconds() - y2015.Seconds()) /
                          (y2020.Seconds() - y2015.Seconds());
    ARTS_ASSERT(scale >= 0 and scale < 1)

    compute_impl(out, mg2020, mh2020, z_field, lat_grid, lon_grid, ell, scale);
    compute_impl(
        out, mg2015, mh2015, z_field, lat_grid, lon_grid, ell, 1.0 - scale);
  } else if (time >= y2010) {
    const Numeric scale = (time.Seconds() - y2010.Seconds()) /
                          (y2015.Seconds() - y2010.Seconds());
    ARTS_ASSERT(scale >= 0 and scale < 1)

    compute_impl(out, mg2015, mh2015, z_field, lat_grid, lon_grid, ell, scale);
    compute_impl(
        out, mg2010, mh2010, z_field, lat_grid, lon_grid, ell, 1.0 - scale);
  } else if (time >= y2005) {
    const Numeric scale = (time.Seconds() - y2005.Seconds()) /
                          (y2010.Seconds() - y2005.Seconds());
    ARTS_ASSERT(scale >= 0 and scale < 1)

    compute_impl(out, mg2010, mh2010, z_field, lat_grid, lon_grid, ell, scale);
    compute_impl(
        out, mg2005, mh2005, z_field, lat_grid, lon_grid, ell, 1.0 - scale);
  } else if (time >= y2000) {
    const Numeric scale = (time.Seconds() - y2000.Seconds()) /
                          (y2005.Seconds() - y2000.Seconds());
    ARTS_ASSERT(scale >= 0 and scale < 1)

    compute_impl(out, mg2005, mh2005, z_field, lat_grid, lon_grid, ell, scale);
    compute_impl(
        out, mg2000, mh2000, z_field, lat_grid, lon_grid, ell, 1.0 - scale);
  } else {
    compute_impl(out, mg2000, mh2000, z_field, lat_grid, lon_grid, ell, 1.0);
  }

  // Conversion from nano-Tesla to Tesla
  out.u *= 1e-9;
  out.v *= 1e-9;
  out.w *= 1e-9;

  return out;
}
}  // namespace IGRF
