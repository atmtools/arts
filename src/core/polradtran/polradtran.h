#pragma once

#include <matpack.h>

/**
 * What the interfaces of Evans' RT3 and RT4 (polradtran::rt3 in rt3.h and
 * polradtran::rt4 in rt4.h) share: the quadrature of their streams and the
 * grounds both have.  The conventions of the streams and the grounds are
 * those of the two solvers (see rt3.h and rt4.h).
 */
namespace polradtran {
//! The quadrature rules of the streams on one hemisphere (QUAD_TYPE of
//! RADTRAN and RADTRANO)
enum class quadrature_type {
  gauss,         //!< 'G' (RT3's 'E' with extra angles): positive half of a 2*nmu-point Gauss-Legendre rule on [-1, 1]
  double_gauss,  //!< 'D': nmu-point Gauss-Legendre rule on [0, 1]
  lobatto,       //!< 'L': positive half of a 2*nmu-point Lobatto rule on [-1, 1]; includes mu = 1
};

//! One hemisphere's streams: ascending mu in (0, 1]; weights for the
//! integral over mu in [0, 1], summing to 1.  The 2 pi azimuth factor is not
//! included.
struct quadrature {
  Vector mu;
  Vector weights;
};

/** The streams of RT3 and RT4, nmu >= 1.
 *
 * ARTS's quadratures (scattering/integration.h) in place of Evans' own: the
 * positive half of scattering::DoubleGaussQuadrature,
 * GaussLegendreQuadrature or LobattoQuadrature of degree 2 nmu.  They are
 * Evans' rules: the nodes agree with RT4's DOUBLE_GAUSS_QUADRATURE,
 * GAUSS_LEGENDRE_QUADRATURE and LOBATTO_QUADRATURE to 4.4e-16 and the
 * weights to 2.4e-12 relative (RT4's Gauss weights are the less accurate),
 * and with RT3's RT3_ versions to rounding.  rt3::radtran and
 * rt4::radtrano use the same.
 */
quadrature get_quadrature(Index nmu, quadrature_type type);

/** 'L': a Lambertian ground of albedo A.  It reflects 2 A mu_j w_j from
 * stream j into every stream in the azimuth mode 0, I to I only, and emits
 * (1 - A) B(surface_temperature) in I.  RT3 also reflects the direct beam,
 * A F_direct / pi, and adds the emission only with its thermal source.
 * Energy is conserved on the streams only for double_gauss quadrature
 * (2 sum mu w = 1); gauss and lobatto are off by about 3e-3 A for 8
 * streams.
 */
struct lambertian_surface {
  Numeric albedo{0.0};
};

/** 'F': specular Fresnel reflection under a medium of index 1, in every
 * azimuth mode: with r_v, r_h the amplitude reflection coefficients,
 * R = [[R1, R2, 0, 0], [R2, R1, 0, 0], [0, 0, R3, -R4], [0, 0, R4, R3]]
 * (RT4 has its [I, Q] block), R1 = (|r_v|^2 + |r_h|^2) / 2,
 * R2 = (|r_v|^2 - |r_h|^2) / 2, R3 = Re(r_v r_h*), R4 = Im(r_v r_h*), and
 * the emission [(1 - R1) B, -R2 B, 0, 0] with B = B(surface_temperature),
 * always (also without RT3's thermal source).  RT3 does not allow a Fresnel
 * surface with a direct beam.
 */
struct fresnel_surface {
  Complex refractive_index{1.0, 0.0};
};
}  // namespace polradtran
