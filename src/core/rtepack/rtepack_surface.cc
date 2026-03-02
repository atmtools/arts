#include "rtepack_surface.h"

#include <arts_constexpr_math.h>

#include "rtepack_multitype.h"

namespace rtepack {
muelmat fresnel_reflectance(Complex Rv, Complex Rh) {
  const Numeric rv    = std::norm(Rv);
  const Numeric rh    = std::norm(Rh);
  const Numeric rmean = 0.5 * (rv + rh);
  const Numeric rdiff = 0.5 * (rv - rh);
  const Complex a     = Rh * std::conj(Rv);
  const Complex b     = Rv * std::conj(Rh);
  const Numeric c     = 0.5 * std::real(a + b);
  const Numeric d     = 0.5 * std::imag(a - b);

  muelmat out{};

  out[0, 0] = rmean;
  out[1, 0] = rdiff;
  out[0, 1] = rdiff;
  out[1, 1] = rmean;
  out[2, 2] = c;
  out[2, 3] = d;
  out[3, 2] = -d;
  out[3, 3] = c;

  return out;
}

namespace {
muelmat stokes_rotation(Numeric cos2psi, Numeric sin2psi) {
  muelmat L{};
  L[0, 0] = 1.0;
  L[1, 1] = cos2psi;
  L[1, 2] = sin2psi;
  L[2, 1] = -sin2psi;
  L[2, 2] = cos2psi;
  L[3, 3] = 1.0;
  return L;
}

muelmat stokes_rotation_refl(Numeric cos2psi, Numeric sin2psi) {
  muelmat L{};
  L[0, 0] = 1.0;
  L[1, 1] = cos2psi;
  L[1, 2] = -sin2psi;
  L[2, 1] = sin2psi;
  L[2, 2] = -cos2psi;
  L[3, 3] = -1.0;
  return L;
}

/** Compute polarization basis vectors (v, h) for a propagation direction k.
 *
 * Uses the local vertical z = (0, 0, 1) as reference.
 * h = k x z / |k x z|,  v = h x k
 */
static std::pair<Vector3, Vector3> pol_basis(const Vector3& k) {
  const Vector3 z{0.0, 0.0, 1.0};
  Vector3 h = cross(k, z);
  const Numeric norm_h = hypot(h);

  if (norm_h < 1e-12) {
    h = Vector3{1.0, 0.0, 0.0};
  } else {
    h[0] /= norm_h;
    h[1] /= norm_h;
    h[2] /= norm_h;
  }

  const Vector3 v = cross(h, k);
  return {v, h};
}
}  // namespace

/** Fresnel reflection Mueller matrix for specular reflection off a
 *  surface with arbitrary tilt.
 *
 *  @param Rv        Complex Fresnel coefficient for vertical polarization
 *  @param Rh        Complex Fresnel coefficient for horizontal polarization
 *  @param k_inc     Unit propagation vector of incident beam (toward surface)
 *  @param n_surface Outward unit surface normal
 *  @return          4x4 Mueller matrix in the incident beam's polarization frame
 *
 *  For a flat horizontal surface (n_surface = (0,0,1)) this reduces to
 *  fresnel_reflectance(Rv, Rh) with U and V sign flips.
 */
muelmat fresnel_reflectance_specular(Complex Rv,
                                     Complex Rh,
                                     const Vector3& k_inc,
                                     const Vector3& n_surface) {
  const muelmat M_fresnel = fresnel_reflectance(Rv, Rh);

  Vector3 m = cross(k_inc, n_surface);
  const Numeric norm_m = hypot(m);

  if (norm_m < 1e-12) {
    muelmat F = muelmat::id();
    F[2, 2] = -1.0;
    F[3, 3] = -1.0;
    return F * M_fresnel;
  }

  m[0] /= norm_m;
  m[1] /= norm_m;
  m[2] /= norm_m;

  const auto [v_i, h_i] = pol_basis(k_inc);

  const Numeric cos_psi1 = dot(h_i, m);
  const Numeric sin_psi1 = dot(v_i, m);
  const Numeric cos2psi1 = 2.0 * cos_psi1 * cos_psi1 - 1.0;
  const Numeric sin2psi1 = 2.0 * sin_psi1 * cos_psi1;

  const muelmat L1 = stokes_rotation(cos2psi1, sin2psi1);
  const muelmat L2 = stokes_rotation_refl(cos2psi1, -sin2psi1);

  return L2 * M_fresnel * L1;
}

/** Fresnel reflection Mueller matrix for non-specular (e.g. Lambertian,
 *  BRDF, or tilted-facet) reflection where incident and outgoing directions
 *  are independent.
 *
 *  @param Rv        Complex Fresnel coefficient for vertical polarization
 *  @param Rh        Complex Fresnel coefficient for horizontal polarization
 *  @param k_inc     Unit propagation vector of incident beam (toward surface)
 *  @param k_out     Unit propagation vector of outgoing beam (toward observer)
 *  @param n_surface Outward unit surface normal
 *  @return          4x4 Mueller matrix mapping incident Stokes to outgoing Stokes
 */
muelmat fresnel_reflectance_nonspecular(Complex Rv,
                                        Complex Rh,
                                        const Vector3& k_inc,
                                        const Vector3& k_out,
                                        const Vector3& n_surface) {
  // Fresnel matrix in the plane-of-incidence frame
  const muelmat M_fresnel = fresnel_reflectance(Rv, Rh);

  // Plane-of-incidence normal: m = k_inc x n / |k_inc x n|
  Vector3 m = cross(k_inc, n_surface);
  const Numeric norm_m = hypot(m);

  if (norm_m < 1e-12) {
    // Normal incidence: plane of incidence undefined
    muelmat F = muelmat::id();
    F[2, 2] = -1.0;
    F[3, 3] = -1.0;
    return F * M_fresnel;
  }

  m[0] /= norm_m;
  m[1] /= norm_m;
  m[2] /= norm_m;

  // Incident beam polarization basis and rotation angle psi_1
  const auto [v_i, h_i] = pol_basis(k_inc);
  const Numeric cos_psi1 = dot(h_i, m);
  const Numeric sin_psi1 = dot(v_i, m);
  const Numeric cos2psi1 = 2.0 * cos_psi1 * cos_psi1 - 1.0;
  const Numeric sin2psi1 = 2.0 * sin_psi1 * cos_psi1;

  const muelmat L1 = stokes_rotation(cos2psi1, sin2psi1);

  // Outgoing beam polarization basis and rotation angle psi_2
  const auto [v_r, h_r] = pol_basis(k_out);
  const Numeric cos_psi2 = dot(m, h_r);
  const Numeric sin_psi2 = dot(m, v_r);
  const Numeric cos2psi2 = 2.0 * cos_psi2 * cos_psi2 - 1.0;
  const Numeric sin2psi2 = 2.0 * sin_psi2 * cos_psi2;

  // Reflected-frame rotation: includes handedness change
  const muelmat L2 = stokes_rotation_refl(cos2psi2, sin2psi2);

  return L2 * M_fresnel * L1;
}

stokvec flat_scalar_reflection(stokvec I, const Numeric R, const stokvec B) {
  I      = I * R;
  I.V() *= -1.0;  //! NOTE: The elementwise multiplication is [R, R, R, -R]
  I     += B * (1.0 - R);
  return I;
}

stokvec dflat_scalar_reflection_dr(stokvec I, const stokvec B) {
  I.V() *= -1.0;  //! NOTE: The elementwise multiplication is [R, R, R, -R]
  I     -= B;
  return I;
}

stokvec reflection(stokvec I, const muelmat R, const stokvec B) {
  I      = R * I;
  I.V() *= -1.0;
  I     += (1.0 - R) * B;
  return I;
}

stokvec dreflection(stokvec I, const muelmat dR, const stokvec B) {
  I      = dR * I;
  I.V() *= -1.0;
  I     -= dR * B;
  return I;
}
}  // namespace rtepack
