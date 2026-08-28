#pragma once

#include <arts_constants.h>
#include <matpack.h>

#include <algorithm>
#include <array>
#include <initializer_list>

#include "reference-aerosol-moments.h"
#include "reference-cloud-moments.h"
#include "reference-problem-17.h"

namespace disort_test::reference {

inline Matrix matrix_2x6(std::initializer_list<Numeric> values) {
  Vector data(values.size());
  std::ranges::copy(values, data.begin());
  return Matrix{std::move(data).reshape(2, 6)};
}

inline Tensor3 tensor_3x6(const Index nphi, std::initializer_list<Numeric> values) {
  Vector data(values.size());
  std::ranges::copy(values, data.begin());
  return Tensor3{std::move(data).reshape(nphi, 3, 6)};
}

struct single_layer_case {
  const char* name;
  Numeric     depth;
  Numeric     omega;
  bool        beam;
  Vector      direct;
  Vector      diffuse_down;
  Vector      up;
  Matrix      radiance;
};

inline const Vector      problem_1_user_mu{-1.0, -0.5, -0.1, 0.1, 0.5, 1.0};
inline constexpr Index   problem_1_streams = 16;
inline constexpr Numeric problem_1_beam_mu = 0.1;
inline const std::array  problem_1{
    single_layer_case{"1a",
                      0.03125,
                      0.2,
                      true,
                      {3.14159, 2.29844},
                      {0, 0.0794108},
                      {0.0799451, 0},
                      matrix_2x6({0, 0, 0, 0.117771, 0.0264170, 0.0134041, 0.0133826, 0.0263324, 0.115898, 0, 0, 0})},
    single_layer_case{"1b",
                      0.03125,
                      1.0,
                      true,
                      {3.14159, 2.29844},
                      {0, 0.420233},
                      {0.422922, 0},
                      matrix_2x6({0, 0, 0, 0.622884, 0.139763, 0.0709192, 0.0708109, 0.139337, 0.613458, 0, 0, 0})},
    single_layer_case{"1c",
                      0.03125,
                      0.99,
                      false,
                      {0, 0},
                      {3.14159, 3.04897},
                      {0.0906556, 0},
                      matrix_2x6({1, 1, 1, 0.133177, 0.0299879, 0.0152233, 0.984447, 0.969363, 0.863946, 0, 0, 0})},
    single_layer_case{
        "1d",
        32.0,
        0.2,
        true,
        {3.14159, 0},
        {0, 0},
        {0.259686, 0},
        matrix_2x6({0, 0, 0, 0.262972, 0.0906967, 0.0502853, 1.22980e-15, 1.30698e-17, 6.88840e-18, 0, 0, 0})},
    single_layer_case{"1e",
                      32.0,
                      1.0,
                      true,
                      {3.14159, 0},
                      {0, 0.0676954},
                      {3.07390, 0},
                      matrix_2x6({0, 0, 0, 1.93321, 1.02732, 0.797199, 0.0271316, 0.0187805, 0.0116385, 0, 0, 0})},
    single_layer_case{
        "1f",
        32.0,
        0.99,
        false,
        {0, 0},
        {3.14159, 0.00460048},
        {2.49618, 0},
        matrix_2x6({1, 1, 1, 0.877510, 0.815136, 0.752715, 0.00186840, 0.00126492, 0.000779280, 0, 0, 0})},
};

inline const Vector      problem_2_user_mu{-0.981986, -0.538263, -0.018014, 0.018014, 0.538263, 0.981986};
inline constexpr Index   problem_2_streams = 16;
inline constexpr Numeric problem_2_beam_mu = 0.080442;
inline const std::array  problem_2{
    single_layer_case{
        "2a",
        0.2,
        0.5,
        true,
        {0.252716, 0.0210311},
        {0, 0.0441791},
        {0.0535063, 0},
        matrix_2x6({0, 0, 0, 0.161796, 0.0211501, 0.00786713, 0.00771897, 0.0200778, 0.0257685, 0, 0, 0})},
    single_layer_case{"2b",
                      0.2,
                      1.0,
                      true,
                      {0.252716, 0.0210311},
                      {0, 0.106123},
                      {0.125561, 0},
                      matrix_2x6({0, 0, 0, 0.347678, 0.0487120, 0.0189387, 0.0186027, 0.0464061, 0.0677603, 0, 0, 0})},
    single_layer_case{
        "2c",
        5.0,
        0.5,
        true,
        {0.252716, 2.56077e-28},
        {0, 0.000251683},
        {0.0624730, 0},
        matrix_2x6({0, 0, 0, 0.162566, 0.0245786, 0.0101498, 0.000170004, 0.0000397168, 0.0000132472, 0, 0, 0})},
    single_layer_case{
        "2d",
        5.0,
        1.0,
        true,
        {0.252716, 0},
        {0, 0.0268008},
        {0.225915, 0},
        matrix_2x6({0, 0, 0, 0.364010, 0.0826993, 0.0492370, 0.0105950, 0.00769337, 0.00379276, 0, 0, 0})},
};

inline const Vector      problem_3_user_mu   = problem_1_user_mu;
inline constexpr Index   problem_3_streams   = 16;
inline constexpr Index   problem_3_moments   = 33;
inline constexpr Numeric problem_3_asymmetry = 0.75;
inline const std::array  problem_3{
    single_layer_case{"3a",
                      1.0,
                      1.0,
                      true,
                      {3.14159, 1.15573},
                      {0, 1.73849},
                      {0.247374, 0},
                      matrix_2x6({0, 0, 0, 0.151159, 0.101103, 0.0395460, 3.05855, 0.266648, 0.213750, 0, 0, 0})},
    single_layer_case{"3b",
                      8.0,
                      1.0,
                      true,
                      {3.14159, 0.00105389},
                      {0, 1.54958},
                      {1.59096, 0},
                      matrix_2x6({0, 0, 0, 0.379740, 0.519598, 0.493302, 0.669581, 0.422350, 0.236362, 0, 0, 0})},
};

inline const Vector haze_l_weighted_coefficients{
    2.41260, 3.23047, 3.37296, 3.23150, 2.89350, 2.49594, 2.11361, 1.74812, 1.44692, 1.17714, 0.96643,
    0.78237, 0.64114, 0.51966, 0.42563, 0.34688, 0.28351, 0.23317, 0.18963, 0.15788, 0.12739, 0.10762,
    0.08597, 0.07381, 0.05828, 0.05089, 0.03971, 0.03524, 0.02720, 0.02451, 0.01874, 0.01711};

inline Matrix haze_l_moments() {
  Matrix moments(1, 33, 0.0);
  moments[0, 0] = 1.0;
  for (Index k = 1; k <= 32; ++k) moments[0, k] = haze_l_weighted_coefficients[k - 1] / static_cast<Numeric>(2 * k + 1);
  return moments;
}

struct haze_l_case {
  const char* name;
  Numeric     omega;
  Numeric     beam_mu;
  Vector      azimuth;
  Vector      direct;
  Vector      diffuse_down;
  Vector      up;
  Tensor3     radiance;
};

inline const Vector      problem_4_user_mu = problem_1_user_mu;
inline const Vector      problem_4_output_tau{0.0, 0.5, 1.0};
inline constexpr Index   problem_4_streams   = 32;
inline constexpr Numeric problem_4_total_tau = 1.0;
inline const std::array  problem_4{
    haze_l_case{"4a",
                1.0,
                1.0,
                {0},
                {3.14159, 1.90547, 1.15573},
                {0, 1.17401, 1.81264},
                {0.173223, 0.111113, 0},
                tensor_3x6(1,
                           {0,
                            0,
                            0,
                            0.0926837,
                            0.0659569,
                            0.0364755,
                            2.51608,
                            0.119287,
                            0.134962,
                            0.123887,
                            0.0402058,
                            0.0177746,
                            3.37302,
                            0.219835,
                            0.156893,
                            0,
                            0,
                            0})},
    haze_l_case{"4b",
                0.9,
                1.0,
                {0},
                {3.14159, 1.90547, 1.15573},
                {0, 1.01517, 1.51554},
                {0.123665, 0.0788690, 0},
                tensor_3x6(1,
                           {0,
                            0,
                            0,
                            0.0653056,
                            0.0455144,
                            0.0282693,
                            2.24258,
                            0.0966049,
                            0.0961335,
                            0.0843278,
                            0.0279473,
                            0.0138835,
                            2.97057,
                            0.167698,
                            0.108115,
                            0,
                            0,
                            0})},
    haze_l_case{
        "4c",
        0.9,
        0.5,
        {0, Constant::pi / 2, Constant::pi},
        {1.57080, 0.577864, 0.212584},
        {0, 0.702764, 0.803294},
        {0.225487, 0.123848, 0},
        tensor_3x6(3,
                   {0,         0,         0,          0.870812,  0.224960,  0.0227572, 0.0477016, 3.02631,   1.41195,
                    0.697692,  0.109130,  0.00932861, 0.0838488, 2.70538,   0.876523,  0,         0,         0,
                    0,         0,         0,          0.0888117, 0.0577411, 0.0227572, 0.0477016, 0.0580971, 0.104502,
                    0.0916071, 0.0295842, 0.00932861, 0.0838488, 0.0942187, 0.0895457, 0,         0,         0,
                    0,         0,         0,          0.0698247, 0.0502877, 0.0227572, 0.0477016, 0.0258544, 0.0625954,
                    0.0591273, 0.0247702, 0.00932861, 0.0838488, 0.0399383, 0.0467155, 0,         0,         0})},
};

inline const Vector cloud_c1_weighted_coefficients{
    1,      2.544,  3.883,  4.568,  5.235,  5.887,  6.457,  7.177,  7.859,  8.494,  9.286,  9.856,  10.615, 11.229,
    11.851, 12.503, 13.058, 13.626, 14.209, 14.660, 15.231, 15.641, 16.126, 16.539, 16.934, 17.325, 17.673, 17.999,
    18.329, 18.588, 18.885, 19.103, 19.345, 19.537, 19.721, 19.884, 20.024, 20.145, 20.251, 20.330, 20.401, 20.444,
    20.477, 20.489, 20.483, 20.467, 20.427, 20.382, 20.310, 20.236, 20.136, 20.036, 19.909, 19.785, 19.632, 19.486,
    19.311, 19.145, 18.949, 18.764, 18.551, 18.348, 18.119, 17.901, 17.659, 17.428, 17.174, 16.931, 16.668, 16.415,
    16.144, 15.883, 15.606, 15.338, 15.058, 14.784, 14.501, 14.225, 13.941, 13.662, 13.378, 13.098, 12.816, 12.536,
    12.257, 11.978, 11.703, 11.427, 11.156, 10.884, 10.618, 10.350, 10.090, 9.827,  9.574,  9.318,  9.072,  8.822,
    8.584,  8.340,  8.110,  7.874,  7.652,  7.424,  7.211,  6.990,  6.785,  6.573,  6.377,  6.173,  5.986,  5.790,
    5.612,  5.424,  5.255,  5.075,  4.915,  4.744,  4.592,  4.429,  4.285,  4.130,  3.994,  3.847,  3.719,  3.580,
    3.459,  3.327,  3.214,  3.090,  2.983,  2.866,  2.766,  2.656,  2.562,  2.459,  2.372,  2.274,  2.193,  2.102,
    2.025,  1.940,  1.869,  1.790,  1.723,  1.649,  1.588,  1.518,  1.461,  1.397,  1.344,  1.284,  1.235,  1.179,
    1.134,  1.082,  1.040,  0.992,  0.954,  0.909,  0.873,  0.832,  0.799,  0.762,  0.731,  0.696,  0.668,  0.636,
    0.610,  0.581,  0.557,  0.530,  0.508,  0.483,  0.463,  0.440,  0.422,  0.401,  0.384,  0.364,  0.349,  0.331,
    0.317,  0.301,  0.288,  0.273,  0.262,  0.248,  0.238,  0.225,  0.215,  0.204,  0.195,  0.185,  0.177,  0.167,
    0.160,  0.151,  0.145,  0.137,  0.131,  0.124,  0.118,  0.112,  0.107,  0.101,  0.097,  0.091,  0.087,  0.082,
    0.079,  0.074,  0.071,  0.067,  0.064,  0.060,  0.057,  0.054,  0.052,  0.049,  0.047,  0.044,  0.042,  0.039,
    0.038,  0.035,  0.034,  0.032,  0.030,  0.029,  0.027,  0.026,  0.024,  0.023,  0.022,  0.021,  0.020,  0.018,
    0.018,  0.017,  0.016,  0.015,  0.014,  0.013,  0.013,  0.012,  0.011,  0.011,  0.010,  0.009,  0.009,  0.008,
    0.008,  0.008,  0.007,  0.007,  0.006,  0.006,  0.006,  0.005,  0.005,  0.005,  0.005,  0.004,  0.004,  0.004,
    0.004,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.002,  0.002,  0.002,  0.002,  0.002,  0.002,  0.002,
    0.002,  0.002,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,
    0.001,  0.001,  0.001,  0.001,  0.001,  0.001};

inline Matrix cloud_c1_moments() {
  Matrix moments(1, cloud_c1_weighted_coefficients.size(), 1.0);
  for (Size k = 1; k < cloud_c1_weighted_coefficients.size(); ++k)
    moments[0, k] = cloud_c1_weighted_coefficients[k] / static_cast<Numeric>(2 * k + 1);
  return moments;
}

struct scalar_case {
  const char* name;
  Numeric     omega;
  Vector      tau;
  Vector      direct;
  Vector      diffuse_down;
  Vector      up;
  Matrix      radiance;
};

inline Matrix radiance(std::initializer_list<Numeric> values) {
  Vector data(values.size());
  std::ranges::copy(values, data.begin());
  return Matrix{std::move(data).reshape(3, 6)};
}

inline const scalar_case problem_5a{"5a",
                                    1.0,
                                    {0.0, 32.0, 64.0},
                                    {3.14159, 3.97856e-14, 5.03852e-28},
                                    {0.0, 2.24768, 0.479851},
                                    {2.66174, 1.76783, 0.0},
                                    radiance({0,
                                              0,
                                              0,
                                              0.458927,
                                              0.772983,
                                              1.07196,
                                              0.753662,
                                              0.696362,
                                              0.650541,
                                              0.627631,
                                              0.581809,
                                              0.524532,
                                              0.195230,
                                              0.131990,
                                              0.0720655,
                                              0,
                                              0,
                                              0})};

inline const scalar_case problem_5b{"5b",
                                    0.9,
                                    {3.2, 12.8, 48.0},
                                    {0.128058, 8.67322e-6, 4.47729e-21},
                                    {1.74767, 0.233975, 6.38345e-5},
                                    {0.270485, 0.0374252, 1.02904e-5},
                                    radiance({67.9623,
                                              0.221027,
                                              0.136619,
                                              0.114084,
                                              0.0873870,
                                              0.0881626,
                                              0.205706,
                                              0.0492736,
                                              0.0265449,
                                              0.0202154,
                                              0.0129661,
                                              0.00951334,
                                              3.41286e-5,
                                              1.39916e-5,
                                              7.47039e-6,
                                              5.65602e-6,
                                              3.58245e-6,
                                              2.57858e-6})};

inline constexpr Index   problem_5_streams   = 48;
inline constexpr Numeric problem_5_total_tau = 64.0;
inline const Vector      problem_5_user_mu   = problem_1_user_mu;

enum class surface_type { black, lambertian, hapke };

/** Solver-neutral physical inputs for original DISORT Problem 6.
 *
 * beam is FBEAM (irradiance), whereas top_isotropic is FISOT (radiance).
 * Temperatures describe Planck radiances over [0, 50000] cm^-1.  A Hapke
 * surface derives its directional emissivity from the BRDF.
 */
struct thermal_source_case {
  const char*            name;
  Numeric                optical_depth;
  Vector                 output_tau;
  Numeric                single_scattering_albedo;
  Numeric                beam_mu;
  Numeric                beam;
  Numeric                top_isotropic;
  Numeric                top_temperature;
  Numeric                top_emissivity;
  Numeric                bottom_temperature;
  Vector                 interface_temperature;
  surface_type           surface;
  Numeric                lambertian_albedo;
  std::array<Numeric, 3> hapke_parameters;
};

inline const Vector      problem_6_user_mu{-1.0, -0.1, 0.1, 1.0};
inline constexpr Numeric problem_6_azimuth         = Constant::pi / 2.0;
inline constexpr Numeric problem_6_wavenumber_low  = 0.0;
inline constexpr Numeric problem_6_wavenumber_high = 50000.0;

inline const std::array problem_6{
    thermal_source_case{"6a", 0.0, {0.0, 0.0}, 0.0, 0.5, 200.0, 0.0, 0.0, 0.0, 0.0, {}, surface_type::black, 0.0, {}},
    thermal_source_case{
        "6b", 1.0, {0.0, 0.5, 1.0}, 0.0, 0.5, 200.0, 0.0, 0.0, 0.0, 0.0, {}, surface_type::black, 0.0, {}},
    thermal_source_case{
        "6c", 1.0, {0.0, 0.5, 1.0}, 0.0, 0.5, 200.0, 0.0, 0.0, 0.0, 0.0, {}, surface_type::lambertian, 0.5, {}},
    thermal_source_case{"6d",
                        1.0,
                        {0.0, 0.5, 1.0},
                        0.0,
                        0.5,
                        200.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        {},
                        surface_type::hapke,
                        0.0,
                        {1.0, 0.06, 0.6}},
    thermal_source_case{"6e",
                        1.0,
                        {0.0, 0.5, 1.0},
                        0.0,
                        0.5,
                        200.0,
                        0.0,
                        0.0,
                        1.0,
                        300.0,
                        {},
                        surface_type::hapke,
                        0.0,
                        {1.0, 0.06, 0.6}},
    thermal_source_case{"6f",
                        1.0,
                        {0.0, 0.5, 1.0},
                        0.0,
                        0.5,
                        200.0,
                        100.0 / Constant::pi,
                        250.0,
                        1.0,
                        300.0,
                        {},
                        surface_type::hapke,
                        0.0,
                        {1.0, 0.06, 0.6}},
    thermal_source_case{"6g",
                        1.0,
                        {0.0, 0.5, 1.0},
                        0.0,
                        0.5,
                        200.0,
                        100.0 / Constant::pi,
                        250.0,
                        1.0,
                        300.0,
                        {250.0, 300.0},
                        surface_type::hapke,
                        0.0,
                        {1.0, 0.06, 0.6}},
    thermal_source_case{"6h",
                        10.0,
                        {0.0, 1.0, 10.0},
                        0.0,
                        0.5,
                        200.0,
                        100.0 / Constant::pi,
                        250.0,
                        1.0,
                        300.0,
                        {250.0, 300.0},
                        surface_type::hapke,
                        0.0,
                        {1.0, 0.06, 0.6}},
};

/** Active DISORT 4.0.99 reference fluxes for Problem 6.  These deliberately
 * follow the uncommented DATA statements in DISOTESTAUX.f; cDISORT still
 * prints the older thermal-emission values for cases 6e--6h. */
struct flux_reference {
  Vector direct;
  Vector diffuse_down;
  Vector up;
  Vector dfdt;
};

inline const std::array problem_6_flux{
    flux_reference{{100.0, 100.0}, {0.0, 0.0}, {0.0, 0.0}, {200.0, 200.0}},
    flux_reference{{100.0, 36.7879, 13.5335}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {200.0, 73.5759, 27.0671}},
    flux_reference{
        {100.0, 36.7879, 13.5335}, {0.0, 0.0, 0.0}, {1.48450, 2.99914, 6.76676}, {202.010, 77.9962, 40.6006}},
    flux_reference{
        {100.0, 36.7879, 13.5335}, {0.0, 0.0, 0.0}, {0.670783, 1.39084, 3.31655}, {200.936, 75.7187, 34.5317}},
    flux_reference{
        {100.0, 36.7879, 13.5335}, {0.0, 0.0, 0.0}, {94.4651, 190.616, 428.806}, {327.710, 353.897, 878.110}},
    flux_reference{{100.0, 36.7879, 13.5335},
                   {321.497, 142.493, 70.5305},
                   {97.7110, 197.247, 441.139},
                   {975.154, 573.874, 1005.98}},
    flux_reference{{100.0, 36.7879, 13.5335},
                   {321.497, 304.775, 363.632},
                   {350.211, 443.254, 513.521},
                   {601.025, 173.906, -7.36023}},
    flux_reference{{100.0, 13.5335, 2.06115e-7},
                   {321.497, 255.455, 443.444},
                   {237.351, 261.131, 528.258},
                   {423.780, 61.9842, 129.477}},
};

struct scattering_thermal_case {
  const char*  name;
  Index        streams;
  Numeric      optical_depth;
  Numeric      single_scattering_albedo;
  Numeric      asymmetry_parameter;
  Numeric      atmosphere_top_temperature;
  Numeric      atmosphere_bottom_temperature;
  Numeric      top_boundary_temperature;
  Numeric      bottom_boundary_temperature;
  Numeric      wavenumber_low;
  Numeric      wavenumber_high;
  Numeric      beam;
  Numeric      top_isotropic;
  surface_type surface;
  Numeric      lambertian_albedo;
};

inline const std::array problem_7{
    scattering_thermal_case{
        "7a", 16, 1.0, 0.10, 0.05, 200.0, 300.0, 0.0, 0.0, 300.0, 800.0, 0.0, 0.0, surface_type::black, 0.0},
    scattering_thermal_case{
        "7b", 16, 100.0, 0.95, 0.75, 200.0, 300.0, 0.0, 0.0, 2702.99, 2703.01, 0.0, 0.0, surface_type::black, 0.0},
    scattering_thermal_case{
        "7c", 12, 1.0, 0.50, 0.80, 300.0, 200.0, 100.0, 320.0, 0.0, 80000.0, 200.0, 100.0, surface_type::black, 0.0},
    scattering_thermal_case{"7d",
                            12,
                            1.0,
                            0.50,
                            0.80,
                            300.0,
                            200.0,
                            100.0,
                            320.0,
                            0.0,
                            80000.0,
                            200.0,
                            100.0,
                            surface_type::lambertian,
                            1.0},
    scattering_thermal_case{
        "7e", 12, 1.0, 0.50, 0.80, 300.0, 200.0, 100.0, 320.0, 0.0, 80000.0, 200.0, 100.0, surface_type::hapke, 0.0},
};

/** Active DISORT 4.0.99 reference fluxes for Problem 7. */
inline const std::array problem_7_flux{
    flux_reference{{0.0, 0.0}, {0.0, 121.204}, {86.2936, 0.0}, {-51.3731, -541.036}},
    flux_reference{{0.0, 0.0}, {0.0, 2.07786e-5}, {1.10949e-6, 0.0}, {8.23219e-8, -5.06461e-6}},
    flux_reference{{100.0, 36.7879, 13.5335},
                   {319.830, 354.099, 301.334},
                   {429.572, 447.018, 594.576},
                   {-80.4270, 251.589, 715.964}},
    flux_reference{{100.0, 36.7879, 13.5335},
                   {319.830, 350.555, 292.063},
                   {312.563, 268.126, 305.596},
                   {-168.356, 101.251, 409.326}},
    flux_reference{{100.0, 36.7879, 13.5335},
                   {319.830, 354.468, 302.366},
                   {440.940, 464.624, 623.842},
                   {-71.6977, 266.918, 750.170}},
};

/** Solver-neutral inputs and DISORT 4.0.99 reference output for Problem 8. */
struct layered_isotropic_case {
  const char* name;
  Vector      cumulative_tau;
  Vector      single_scattering_albedo;
  Vector      output_tau;
  Vector      direct;
  Vector      diffuse_down;
  Vector      up;
  Matrix      radiance;
};

inline Matrix radiance_3x4(std::initializer_list<Numeric> values) {
  Vector data(values.size());
  std::ranges::copy(values, data.begin());
  return Matrix{std::move(data).reshape(3, 4)};
}

inline const Vector      problem_8_user_mu{-1.0, -0.2, 0.2, 1.0};
inline constexpr Index   problem_8_streams = 8;
inline constexpr Numeric problem_8_azimuth = Constant::pi / 3.0;

inline const std::array problem_8{
    layered_isotropic_case{"8a",
                           {0.25, 0.5},
                           {0.5, 0.3},
                           {0.0, 0.25, 0.5},
                           {0.0, 0.0, 0.0},
                           {1.00000, 0.722235, 0.513132},
                           {0.0929633, 0.0278952, 0.0},
                           radiance_3x4({0.318310,
                                         0.318310,
                                         0.0562566,
                                         0.0194423,
                                         0.262711,
                                         0.136952,
                                         0.0184909,
                                         0.00552188,
                                         0.210014,
                                         0.0560376,
                                         0.0,
                                         0.0})},
    layered_isotropic_case{"8b",
                           {0.25, 0.5},
                           {0.8, 0.95},
                           {0.0, 0.25, 0.5},
                           {0.0, 0.0, 0.0},
                           {1.00000, 0.795332, 0.650417},
                           {0.225136, 0.126349, 0.0},
                           radiance_3x4({0.318310,
                                         0.318310,
                                         0.123687,
                                         0.0495581,
                                         0.277499,
                                         0.183950,
                                         0.0835695,
                                         0.0250575,
                                         0.240731,
                                         0.129291,
                                         0.0,
                                         0.0})},
    layered_isotropic_case{"8c",
                           {1.0, 3.0},
                           {0.8, 0.95},
                           {0.0, 1.0, 3.0},
                           {0.0, 0.0, 0.0},
                           {1.00000, 0.486157, 0.159984},
                           {0.378578, 0.243397, 0.0},
                           radiance_3x4({0.318310,
                                         0.318310,
                                         0.149335,
                                         0.104766,
                                         0.189020,
                                         0.0988158,
                                         0.0965192,
                                         0.0654445,
                                         0.0684762,
                                         0.0296698,
                                         0.0,
                                         0.0})},
};

enum class phase_type { isotropic, problem_9b, henyey_greenstein };

/** Solver-neutral physical inputs and DISORT 4.0.99 output for Problem 9. */
struct general_multilayer_case {
  const char* name;
  phase_type  phase;
  Vector      cumulative_tau;
  Vector      single_scattering_albedo;
  Vector      output_tau;
  Vector      azimuth;
  Numeric     beam_mu;
  Numeric     beam;
  Numeric     top_isotropic;
  bool        thermal;
  Numeric     top_temperature;
  Numeric     bottom_temperature;
  Vector      interface_temperature;
  Numeric     wavenumber_low;
  Numeric     wavenumber_high;
  Numeric     surface_albedo;
  Vector      direct;
  Vector      diffuse_down;
  Vector      up;
  Tensor3     radiance;
};

inline Tensor3 radiance_azimuth_tau_mu(const Index nazimuth, std::initializer_list<Numeric> values) {
  Vector data(values.size());
  std::ranges::copy(values, data.begin());
  return Tensor3{std::move(data).reshape(nazimuth, 5, 4)};
}

inline const Vector    problem_9_user_mu{-1.0, -0.2, 0.2, 1.0};
inline constexpr Index problem_9_streams = 8;
inline const Vector    problem_9_cumulative_tau{1.0, 3.0, 6.0, 10.0, 15.0, 21.0};
inline const Vector    problem_9_omega{0.65, 0.70, 0.75, 0.80, 0.85, 0.90};
inline const Vector    problem_9_output_tau{0.0, 1.05, 2.1, 6.0, 21.0};
inline const Vector    problem_9b_moments{1.0,
                                          2.00916 / 3.0,
                                          1.56339 / 5.0,
                                          0.67407 / 7.0,
                                          0.22215 / 9.0,
                                          0.04725 / 11.0,
                                          0.00671 / 13.0,
                                          0.00068 / 15.0,
                                          0.00005 / 17.0};
inline const Vector    problem_9c_asymmetry{1.0 / 7.0, 2.0 / 7.0, 3.0 / 7.0, 4.0 / 7.0, 5.0 / 7.0, 6.0 / 7.0};

inline const std::array problem_9{
    general_multilayer_case{"9a",
                            phase_type::isotropic,
                            problem_9_cumulative_tau,
                            problem_9_omega,
                            problem_9_output_tau,
                            {Constant::pi / 3.0},
                            0.5,
                            0.0,
                            Constant::inv_pi,
                            false,
                            0.0,
                            0.0,
                            {},
                            0.0,
                            0.0,
                            0.0,
                            {0, 0, 0, 0, 0},
                            {1.00000, 0.355151, 0.144265, 0.00671445, 6.16968e-7},
                            {0.227973, 0.0875098, 0.0361819, 0.00219291, 0.0},
                            radiance_azimuth_tau_mu(1, {0.318310,    0.318310,   0.0998915,  0.0591345,  0.153507,
                                                        0.0509531,   0.0367006,  0.0231903,  0.0706614,  0.0209119,
                                                        0.0148545,   0.00972307, 0.00372784, 0.00108815, 0.000883316,
                                                        0.000594743, 2.87656e-7, 1.05921e-7, 0.0,        0.0})},
    general_multilayer_case{
        "9b",
        phase_type::problem_9b,
        problem_9_cumulative_tau,
        problem_9_omega,
        problem_9_output_tau,
        {Constant::pi / 3.0},
        0.5,
        0.0,
        Constant::inv_pi,
        false,
        0.0,
        0.0,
        {},
        0.0,
        0.0,
        0.0,
        {0, 0, 0, 0, 0},
        {1.00000, 0.452357, 0.236473, 0.0276475, 7.41853e-5},
        {0.100079, 0.0452014, 0.0241941, 0.00416016, 0.0},
        radiance_azimuth_tau_mu(1, {0.318310,   0.318310,    0.0739198,  0.0132768,  0.196609,   0.0592369, 0.0300230,
                                    0.00705566, 0.115478,    0.0301809,  0.0152672,  0.00406932, 0.0146177, 0.00385590,
                                    0.00238301, 0.000777890, 3.37742e-5, 1.20858e-5, 0.0,        0.0})},
    general_multilayer_case{
        "9c",
        phase_type::henyey_greenstein,
        problem_9_cumulative_tau,
        problem_9_omega,
        problem_9_output_tau,
        {Constant::pi / 3.0, 2.0 * Constant::pi / 3.0, Constant::pi},
        0.5,
        Constant::pi,
        1.0,
        true,
        550.0,
        700.0,
        {600.0, 610.0, 620.0, 630.0, 640.0, 650.0, 660.0},
        999.0,
        1000.0,
        0.5,
        {1.57080, 0.192354, 0.0235550, 9.65131e-6, 9.03133e-19},
        {6.09217, 4.97279, 4.46616, 4.22731, 4.73767},
        {4.68414, 4.24381, 4.16941, 4.30667, 5.11524},
        radiance_azimuth_tau_mu(
            3, {1.93920, 1.93920, 1.61855, 1.43872, 1.66764, 1.44453, 1.38339, 1.33890, 1.48511, 1.35009,
                1.33079, 1.32794, 1.34514, 1.35131, 1.35980, 1.37918, 1.48927, 1.54270, 1.62823, 1.62823,
                1.93920, 1.93920, 1.57895, 1.43872, 1.66764, 1.42925, 1.37317, 1.33890, 1.48511, 1.34587,
                1.32921, 1.32794, 1.34514, 1.35129, 1.35979, 1.37918, 1.48927, 1.54270, 1.62823, 1.62823,
                1.93920, 1.93920, 1.56559, 1.43872, 1.66764, 1.42444, 1.37034, 1.33890, 1.48511, 1.34469,
                1.32873, 1.32794, 1.34514, 1.35128, 1.35979, 1.37918, 1.48927, 1.54270, 1.62823, 1.62823})},
};

inline const Vector    problem_10_output_tau{0.0, 2.1, 21.0};
inline const Vector    problem_10_azimuth{Constant::pi / 3.0, 2.0 * Constant::pi / 3.0};
inline const Vector    problem_10_user_mu{-0.788675129, -0.211324871, 0.211324871, 0.788675129};
inline constexpr Index problem_10_streams = 4;

inline const Vector      problem_11_output_tau{0.0, 0.05, 0.5, 1.0};
inline const Vector      problem_11_user_mu{-1.0, -0.1, 0.1, 1.0};
inline const Vector      problem_11_azimuth{0.0, Constant::pi / 2.0};
inline const Vector      problem_11_subdivided_tau{0.05, 0.5, 1.0};
inline constexpr Index   problem_11_streams        = 16;
inline constexpr Numeric problem_11_omega          = 0.9;
inline constexpr Numeric problem_11_beam_mu        = 0.5;
inline constexpr Numeric problem_11_beam           = 1.0;
inline constexpr Numeric problem_11_top_isotropic  = 0.5 * Constant::inv_pi;
inline constexpr Numeric problem_11_surface_albedo = 0.5;

inline const Vector      problem_12_output_tau{0.0, 10.0, 19.9, 20.1};
inline const Vector      problem_12_user_mu{-1.0, -0.1, 0.1, 1.0};
inline constexpr Numeric problem_12_azimuth = 0.0;
inline const Vector      problem_12_subdivided_tau{10.0, 19.9, 20.1};
inline constexpr Index   problem_12_streams        = 20;
inline constexpr Numeric problem_12_omega          = 0.5;
inline constexpr Numeric problem_12_asymmetry      = 0.9;
inline constexpr Numeric problem_12_beam_mu        = 1.0;
inline constexpr Numeric problem_12_beam           = 1.0;
inline constexpr Numeric problem_12_surface_albedo = 1.0;

/** Problem 13 values.  The shortcut cases 13a/13c return the same two
 * observables as the corresponding regular calculations 13b/13d. */
struct albedo_transmission_reference {
  const char* name;
  Vector      cumulative_tau;
  Vector      single_scattering_albedo;
  Numeric     albedo;
  Numeric     transmission;
};

inline constexpr Index   problem_13_streams        = 16;
inline constexpr Numeric problem_13_beam_mu        = 0.5;
inline constexpr Numeric problem_13_beam           = 2.0;
inline constexpr Numeric problem_13_surface_albedo = 0.5;
inline constexpr Numeric problem_13_asymmetry      = 0.8;
inline const std::array  problem_13{
    albedo_transmission_reference{"13a/13b", {1.0}, {0.99}, 0.54526, 0.84500},
    albedo_transmission_reference{"13c/13d", {0.5, 1.0}, {0.99, 0.50}, 0.27620, 0.50332},
};

enum class brdf_type { hapke, cox_munk, rpv, ross_li };

struct transparent_brdf_case {
  const char* name;
  brdf_type   type;
  Vector      azimuth;
  Numeric     direct;
  Numeric     diffuse_down;
  Numeric     up;
  Numeric     dfdt;
  Matrix      radiance;
};

inline Matrix radiance_3x4_phi_mu(std::initializer_list<Numeric> values) {
  Vector data(values.size());
  std::ranges::copy(values, data.begin());
  return Matrix{std::move(data).reshape(3, 4)};
}

inline constexpr Index   problem_14_streams = 32;
inline constexpr Numeric problem_14_beam_mu = 0.8660254037844386;
inline constexpr Numeric problem_14_beam    = 1.0;
inline const Vector      problem_14_user_mu{0.1, 0.2, 0.5, 1.0};

/** DISORT 4.0.99 Problem 14.  Radiance is stored as [azimuth, user_mu]. */
inline const std::array problem_14{
    transparent_brdf_case{"14a",
                          brdf_type::hapke,
                          {0.0, Constant::pi / 2.0, Constant::pi},
                          0.866025,
                          0.0,
                          0.180998,
                          1.37365,
                          radiance_3x4_phi_mu({0.0519252,
                                               0.0517173,
                                               0.0500654,
                                               0.0536730,
                                               0.0640461,
                                               0.0626774,
                                               0.0581123,
                                               0.0536730,
                                               0.0777477,
                                               0.0754426,
                                               0.0693954,
                                               0.0536730})},
    transparent_brdf_case{"14b",
                          brdf_type::cox_munk,
                          {0.0, Constant::pi / 4.0, Constant::pi / 2.0},
                          0.866025,
                          0.0,
                          0.0203449,
                          1.03199,
                          radiance_3x4_phi_mu({0.0169062,
                                               0.0164195,
                                               0.0268003,
                                               0.00985518,
                                               0.000160887,
                                               0.000387836,
                                               0.00374155,
                                               0.00985518,
                                               2.40499e-9,
                                               4.45562e-8,
                                               3.22034e-5,
                                               0.00985518})},
    transparent_brdf_case{"14c",
                          brdf_type::rpv,
                          {0.0, Constant::pi / 2.0, Constant::pi},
                          0.866025,
                          0.0,
                          0.148584,
                          1.32427,
                          radiance_3x4_phi_mu({0.0466076,
                                               0.0388121,
                                               0.0333376,
                                               0.0477098,
                                               0.0590715,
                                               0.0497124,
                                               0.0435047,
                                               0.0477099,
                                               0.0782148,
                                               0.0668223,
                                               0.0614060,
                                               0.0477098})},
    transparent_brdf_case{"14d",
                          brdf_type::ross_li,
                          {0.0, Constant::pi / 2.0, Constant::pi},
                          0.866025,
                          0.0,
                          0.213718,
                          1.38147,
                          radiance_3x4_phi_mu({0.0142562,
                                               0.0441596,
                                               0.0611839,
                                               0.0727943,
                                               0.0385329,
                                               0.0561919,
                                               0.0660729,
                                               0.0727943,
                                               0.0648269,
                                               0.0700422,
                                               0.0744868,
                                               0.0727943})},
};

struct layered_brdf_case {
  const char* name;
  brdf_type   type;
  Vector      direct;
  Vector      diffuse_down;
  Vector      up;
  Vector      dfdt;
  Tensor3     radiance;
};

inline Tensor3 radiance_3x3x4(std::initializer_list<Numeric> values) {
  Vector data(values.size());
  std::ranges::copy(values, data.begin());
  return Tensor3{std::move(data).reshape(3, 3, 4)};
}

inline constexpr Index   problem_15_streams = 32;
inline constexpr Numeric problem_15_beam_mu = problem_14_beam_mu;
inline constexpr Numeric problem_15_beam    = 1.0;
inline const Vector      problem_15_tau{0.32, 0.64};
inline const Vector      problem_15_output_tau{0.0, 0.4, 0.64};
inline const Vector      problem_15_user_mu = problem_14_user_mu;
inline const Vector      problem_15_azimuth{0.0, Constant::pi / 2.0, Constant::pi};

/** DISORT 4.0.99 Problem 15.  Radiance is [azimuth, output_tau, user_mu]. */
inline const std::array problem_15{
    layered_brdf_case{
        "15a",
        brdf_type::hapke,
        {0.866025, 0.545681, 0.413603},
        {0.0, 0.221186, 0.338138},
        {0.277760, 0.178593, 0.163465},
        {1.91452e-6, 1.80279e-6, 1.74169e-6},
        radiance_3x3x4({0.123554,  0.112727,  0.0833376, 0.0788616, 0.0978420, 0.0764028, 0.0540164, 0.0513113,
                        0.0555814, 0.0530299, 0.0479406, 0.0469991, 0.113583,  0.107286,  0.0890993, 0.0788616,
                        0.0758666, 0.0687194, 0.0570572, 0.0513113, 0.0641686, 0.0606886, 0.0534627, 0.0469992,
                        0.131842,  0.126113,  0.110240,  0.0788616, 0.0767841, 0.0736885, 0.0689669, 0.0513113,
                        0.0739453, 0.0696550, 0.0610954, 0.0469991})},
    layered_brdf_case{
        "15b",
        brdf_type::cox_munk,
        {0.866025, 0.545681, 0.413603},
        {0.0, 0.190697, 0.302160},
        {0.173025, 0.0433759, 0.0227533},
        {1.68943e-6, 1.38905e-6, 1.27981e-6},
        radiance_3x3x4({0.105039,  0.0919789,  0.0596677,  0.0442653,  0.0728008, 0.0502279,  0.0272793,  0.0114839,
                        0.0348211, 0.0279014,  0.0221136,  0.00709119, 0.0944399, 0.0832705,  0.0547476,  0.0442653,
                        0.0386983, 0.0255378,  0.00951889, 0.0114839,  0.0137107, 0.00944994, 0.00346688, 0.00709119,
                        0.112580,  0.100789,   0.0718390,  0.0442653,  0.0343420, 0.0233428,  0.0137850,  0.0114839,
                        0.0117920, 0.00806508, 0.00279583, 0.00709119})},
    layered_brdf_case{
        "15c",
        brdf_type::rpv,
        {0.866025, 0.545681, 0.413603},
        {0.0, 0.214201, 0.330068},
        {0.254555, 0.148380, 0.132172},
        {1.86327e-6, 1.70778e-6, 1.65841e-6},
        radiance_3x3x4({0.119166,  0.107082,  0.0741418, 0.0726226, 0.0900393, 0.0659806, 0.0405345, 0.0445512,
                        0.0494497, 0.0405951, 0.0333681, 0.0404333, 0.109205,  0.101727,  0.0806139, 0.0726226,
                        0.0684651, 0.0586094, 0.0448932, 0.0445511, 0.0582113, 0.0482087, 0.0403103, 0.0404332,
                        0.127564,  0.120988,  0.103628,  0.0726226, 0.0710388, 0.0656529, 0.0604316, 0.0445511,
                        0.0713241, 0.0598879, 0.0521536, 0.0404332})},
    layered_brdf_case{
        "15d",
        brdf_type::ross_li,
        {0.866025, 0.545681, 0.413603},
        {0.0, 0.225389, 0.342092},
        {0.297775, 0.202805, 0.187433},
        {1.95071e-6, 1.84998e-6, 1.77781e-6},
        radiance_3x3x4({0.126062,  0.114608,  0.0879932, 0.0899278, 0.0934865, 0.0733495, 0.0602902, 0.0654998,
                        0.0278309, 0.0454440, 0.0549975, 0.0619027, 0.116114,  0.109216,  0.0929265, 0.0899278,
                        0.0730116, 0.0661578, 0.0617029, 0.0654998, 0.0454126, 0.0541580, 0.0585154, 0.0619027,
                        0.134402,  0.128103,  0.113285,  0.0899278, 0.0757499, 0.0717037, 0.0720841, 0.0654998,
                        0.0643650, 0.0640986, 0.0642596, 0.0619027})},
};
}  // namespace disort_test::reference
