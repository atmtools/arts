/* The polarization geometry of rtepack against electric fields.

   Every expected Stokes vector here is computed from a field by the
   definition I = |E_v|^2 + |E_h|^2, Q = |E_v|^2 - |E_h|^2,
   U = 2 Re(E_v E_h*), V = 2 Im(E_v E_h*), with the field components
   projected on the basis vectors.  No Mueller formula is shared with the
   code under test.  Random pure states span the Stokes space, so M S = S'
   for several of them tests all of M. */

#include <rtepack.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <format>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

namespace {
using rtepack::polarization_basis;

void require(bool ok, const std::string& what) {
  if (not ok) throw std::runtime_error(what);
}

using Field = std::array<Complex, 3>;

std::mt19937_64                         rng(42);
std::uniform_real_distribution<Numeric> uniform(-1.0, 1.0);

Complex random_complex() { return {uniform(rng), uniform(rng)}; }

Vector3 direction(Numeric mu, Numeric phi) {
  const Numeric s = std::sqrt(1 - mu * mu);
  return {s * std::cos(phi), s * std::sin(phi), mu};
}

Field field(const polarization_basis& b, Complex Ev, Complex Eh) {
  return {Ev * b.first[0] + Eh * b.second[0], Ev * b.first[1] + Eh * b.second[1], Ev * b.first[2] + Eh * b.second[2]};
}

Complex project(const Field& E, const Vector3& e) { return E[0] * e[0] + E[1] * e[1] + E[2] * e[2]; }

rtepack::stokvec stokes(const Field& E, const polarization_basis& b) {
  const Complex Ev = project(E, b.first), Eh = project(E, b.second);
  return {std::norm(Ev) + std::norm(Eh),
          std::norm(Ev) - std::norm(Eh),
          2 * std::real(Ev * std::conj(Eh)),
          2 * std::imag(Ev * std::conj(Eh))};
}

Numeric max_abs(const rtepack::stokvec& s) {
  return std::max({std::abs(s[0]), std::abs(s[1]), std::abs(s[2]), std::abs(s[3])});
}

//! max |M S_in - S_out| / I_in over random pure fields, S_out from the field map out = map(in)
Numeric deviation(const rtepack::muelmat&            M,
                  const polarization_basis&          in,
                  const polarization_basis&          out,
                  const std::function<Field(Field)>& map) {
  Numeric worst = 0.0;
  for (int i = 0; i < 8; i++) {
    const Field E   = field(in, random_complex(), random_complex());
    const auto  Sin = stokes(E, in);
    worst           = std::max(worst, max_abs(M * Sin - stokes(map(E), out)) / Sin[0]);
  }
  return worst;
}

Field scaled_sum(Complex a, const Vector3& x, Complex b, const Vector3& y) {
  return {a * x[0] + b * y[0], a * x[1] + b * y[1], a * x[2] + b * y[2]};
}

constexpr Numeric tol = 1e-13;

//! Orthonormal, right-handed with k, v . z >= 0, h horizontal; along +-z h = h_pole
void test_meridional_basis() {
  const Vector3 x{1, 0, 0}, z{0, 0, 1};
  for (const Numeric mu : {-0.999, -0.6, -0.1, 0.0, 0.3, 0.8, 0.9999}) {
    for (const Numeric phi : {0.0, 0.7, 2.0, 3.5, 5.9}) {
      const Vector3 k      = direction(mu, phi);
      const auto [v, h]    = rtepack::meridional_basis(k, x);
      const Vector3 vxh    = cross(v, h);
      const Numeric errors = std::max(
          {std::abs(dot(v, v) - 1), std::abs(dot(h, h) - 1), std::abs(dot(v, h)), hypot(vxh - k), std::abs(h[2])});
      require(errors < tol and dot(v, z) >= 0, std::format("meridional basis at mu {} phi {}: {}", mu, phi, errors));
    }
  }
  const Vector3 h_pole{std::sin(0.4), -std::cos(0.4), 0};
  for (const Numeric mu : {-1.0, 1.0}) {
    const auto [v, h] = rtepack::meridional_basis(Vector3{0, 0, mu}, h_pole);
    const auto limit  = rtepack::meridional_basis(direction(mu * (1 - 1e-14), 0.4), x);
    require(hypot(h - h_pole) < tol and hypot(h - limit.second) < 1e-6 and hypot(v - limit.first) < 1e-6,
            std::format("the pole at mu {} must be the limit along the meridian of h_pole", mu));
  }
}

//! Two bases of the same ray see the same field
void test_basis_rotation() {
  const Vector3 x{1, 0, 0};
  for (const Numeric mu : {-0.7, 0.2, 0.95}) {
    for (const Numeric psi : {0.0, 0.3, 1.2, 2.5, -0.8}) {
      const polarization_basis from = rtepack::meridional_basis(direction(mu, 1.1), x);
      const auto& [v, h]            = from;
      const polarization_basis to{std::cos(psi) * v + std::sin(psi) * h, std::cos(psi) * h - std::sin(psi) * v};
      const auto               M = rtepack::basis_rotation(from, to);
      const Numeric            d = deviation(M, from, to, [](Field E) { return E; });
      require(d < tol, std::format("basis rotation by {} at mu {}: {}", psi, mu, d));
    }
  }

  bool threw = false;
  try {
    const auto a = rtepack::meridional_basis(direction(0.3, 0.0), x);
    static_cast<void>(rtepack::basis_rotation(a, polarization_basis{a.second, a.first}));
  } catch (const std::exception&) { threw = true; }
  require(threw, "a basis of the opposite handedness is not a rotation and must throw");
}

//! E_out = J E_in in the bases' components, and the Jones matrix of that Mueller matrix
void test_jones() {
  const Vector3 x{1, 0, 0};
  for (int trial = 0; trial < 20; trial++) {
    const rtepack::jones     J{random_complex(), random_complex(), random_complex(), random_complex()};
    const polarization_basis in  = rtepack::meridional_basis(direction(uniform(rng), 3 * uniform(rng)), x);
    const polarization_basis out = rtepack::meridional_basis(direction(uniform(rng), 3 * uniform(rng)), x);
    const auto               map = [&](Field E) {
      const Complex Ev = project(E, in.first), Eh = project(E, in.second);
      return field(out, J[0, 0] * Ev + J[0, 1] * Eh, J[1, 0] * Ev + J[1, 1] * Eh);
    };
    const rtepack::muelmat M = rtepack::to_muelmat(J);
    const Numeric          d = deviation(M, in, out, map);
    require(d < tol, std::format("Jones matrix {}: {}", trial, d));

    // Back to J up to a global phase r: J_back = r J with |r| = 1
    const rtepack::jones J_back{M};
    Complex              overlap{};
    Numeric              norm = 0.0;
    for (Index i = 0; i < 2; i++)
      for (Index j = 0; j < 2; j++) {
        overlap += J_back[i, j] * std::conj(J[i, j]);
        norm    += std::norm(J[i, j]);
      }
    const Complex r     = overlap / norm;
    Numeric       worst = std::abs(std::abs(r) - 1.0);
    for (Index i = 0; i < 2; i++)
      for (Index j = 0; j < 2; j++) worst = std::max(worst, std::abs(J_back[i, j] - r * J[i, j]) / std::sqrt(norm));
    require(worst < 1e-12, std::format("Jones matrix {} from its Mueller matrix: {}", trial, worst));

    // The product of Jones matrices is the product of their Mueller matrices
    const rtepack::jones   K{random_complex(), random_complex(), random_complex(), random_complex()};
    const rtepack::muelmat MK      = rtepack::to_muelmat(J * K) - rtepack::to_muelmat(J) * rtepack::to_muelmat(K);
    Numeric                product = 0.0;
    for (const Numeric e : MK.data) product = std::max(product, std::abs(e));
    require(product < tol * 10, std::format("Jones product {}: {}", trial, product));
  }

  // The basis rotation and the diagonal Jones matrices
  for (const Numeric psi : {0.0, 0.4, 1.3, -2.2}) {
    const rtepack::muelmat a = rtepack::to_muelmat(rtepack::jones::rotation(std::cos(psi), std::sin(psi)));
    const rtepack::muelmat b = rtepack::stokes_rotation(std::cos(2 * psi), std::sin(2 * psi));
    for (Size n = 0; n < 16; n++) require(std::abs(a.data[n] - b.data[n]) < tol, "Jones basis rotation");
  }
  const Complex          Rv{0.6, 0.2}, Rh{-0.7, 0.1};
  const rtepack::muelmat a = rtepack::to_muelmat(rtepack::jones::diagonal(Rv, Rh));
  const rtepack::muelmat b = rtepack::fresnel_reflectance(Rv, Rh);
  for (Size n = 0; n < 16; n++) require(std::abs(a.data[n] - b.data[n]) < tol, "Jones of Fresnel reflection");

  // A depolarizing Mueller matrix has no Jones matrix
  const rtepack::muelmat mixed =
      0.5 * (rtepack::to_muelmat(rtepack::jones::diagonal(1.0, -1.0)) + rtepack::to_muelmat(rtepack::jones::id()));
  bool threw = false;
  try {
    static_cast<void>(rtepack::jones{mixed});
  } catch (const std::exception&) { threw = true; }
  require(threw, "A depolarizing Mueller matrix must not give a Jones matrix");
}

// The fixed-size arithmetic of rtepack is evaluated at compile time
#if __cpp_lib_constexpr_complex >= 201711L
static_assert(rtepack::to_muelmat(rtepack::jones::diagonal(1.0, -1.0))[3, 3] == -1.0);
static_assert(rtepack::fresnel_reflectance(1.0, -1.0)[2, 2] == -1.0);
#endif
static_assert(rtepack::power_reflectance(0.5)[3, 3] == -0.5);
static_assert(rtepack::reflection(rtepack::stokvec{1.0, 0.0, 0.0, 1.0}, rtepack::power_reflectance(1.0), {})[3] ==
              -1.0);

/* A scattering-plane amplitude matrix diag(S2, S1) on (par, perp), the
   field E_out = S2 (E . par_in) par_out + S1 (E . N) N, against
   lab_frame_phase_matrix of its scattering-plane matrix.  Rayleigh has
   S2 = cos(Theta) and S1 = 1; complex amplitudes give F34. */
void test_lab_frame_phase_matrix() {
  const Vector3 x{1, 0, 0};
  const auto    check =
      [&](const polarization_basis& in, const polarization_basis& out, Complex S2, Complex S1, bool mirrored = false) {
        const Vector3 k_in = cross(in.first, in.second), k_out = cross(out.first, out.second);
        Vector3       N      = cross(k_in, k_out);
        N                    = hypot(N) < 1e-12 ? in.second : (1 / hypot(N)) * N;
        const Vector3 par_in = cross(N, k_in), par_out = cross(N, k_out);
        const auto map = [&](Field E) { return scaled_sum(S2 * project(E, par_in), par_out, S1 * project(E, N), N); };

        // The scattering-plane matrix of diag(S2, S1) from its definition, Q = I_par - I_perp
        const Numeric                         a = std::norm(S2), b = std::norm(S1);
        const Complex                         c = S2 * std::conj(S1);
        const rtepack::compact_planar_muelmat F{
            0.5 * (a + b), 0.5 * (a - b), 0.5 * (a + b), std::real(c), -std::imag(c), std::real(c)};
        const auto Z = rtepack::lab_frame_phase_matrix(F, in, out);
        return deviation(mirrored ? rtepack::mirror(Z) : Z, in, out, map);
      };

  for (const Numeric mu_in : {-0.9, -0.3, 0.4, 0.85}) {
    for (const Numeric mu_out : {-0.95, -0.2, 0.5, 0.99}) {
      for (const Numeric phi : {0.0, 0.4, 1.9, 3.14159, 4.4}) {
        const polarization_basis in        = rtepack::meridional_basis(direction(mu_in, 0.3), x);
        const polarization_basis out       = rtepack::meridional_basis(direction(mu_out, 0.3 + phi), x);
        const Numeric            cos_theta = dot(cross(in.first, in.second), cross(out.first, out.second));
        const Numeric            rayleigh  = check(in, out, cos_theta, 1.0);
        const Numeric            mie       = check(in, out, Complex{0.3, -0.8}, Complex{-0.5, 0.6});
        require(rayleigh < tol and mie < tol,
                std::format("phase matrix mu_in {} mu_out {} phi {}: {} {}", mu_in, mu_out, phi, rayleigh, mie));
      }
    }
  }

  // Not blind to the side of the principal plane: the mirror image fails
  const polarization_basis in  = rtepack::meridional_basis(direction(-0.3, 0.3), x);
  const polarization_basis out = rtepack::meridional_basis(direction(0.5, 2.2), x);
  require(check(in, out, Complex{0.3, -0.8}, Complex{-0.5, 0.6}, true) > 1e-2,
          "The mirrored phase matrix must fail the field comparison");

  // Exactly forward and backward, also at the poles: any plane serves a regular matrix
  for (const Numeric mu : {-1.0, -0.6, 0.0, 0.7, 1.0}) {
    const Vector3            h_pole{std::sin(1.3), -std::cos(1.3), 0};
    const polarization_basis in      = rtepack::meridional_basis(direction(mu, 1.3), h_pole);
    const polarization_basis forward = in;
    const polarization_basis backward =
        rtepack::meridional_basis(-1.0 * cross(in.first, in.second), Vector3{-h_pole[0], -h_pole[1], 0});
    const Numeric f = check(in, forward, Complex{0.7, 0.2}, Complex{0.7, 0.2});
    const Numeric r = check(in, backward, Complex{-0.4, 0.1}, Complex{0.4, -0.1});
    require(f < tol and r < tol, std::format("forward and backward at mu {}: {} {}", mu, f, r));
  }
}

/* Fresnel reflection by a facet whose normal bisects k_in and k_out: the
   field E_r = Rh (E . s) s + Rv (E . p_in) p_out, s across the plane of
   incidence and p = s x k, and a perfect mirror E_r = 2 (E . n) n - E,
   which reverses circular polarization. */
void test_facet_reflection() {
  const Vector3 x{1, 0, 0};
  for (const Numeric mu_in : {-0.95, -0.6, -0.2}) {
    for (const Numeric mu_out : {0.15, 0.5, 0.9}) {
      for (const Numeric phi : {0.0, 0.8, 2.2, 3.14159265358979, 4.0}) {
        const Vector3 k_in = direction(mu_in, 0.5), k_out = direction(mu_out, 0.5 + phi);
        Vector3       n         = k_out - k_in;
        n                       = (1 / hypot(n)) * n;
        const Numeric cos_local = -dot(k_in, n);
        const Complex m{1.8, 0.3};
        const Complex root = std::sqrt(m * m - (1.0 - cos_local * cos_local));
        const Complex Rv   = (m * m * cos_local - root) / (m * m * cos_local + root);
        const Complex Rh   = (cos_local - root) / (cos_local + root);

        Vector3 s          = cross(k_in, n);
        s                  = (1 / hypot(s)) * s;
        const Vector3 p_in = cross(s, k_in), p_out = cross(s, k_out);
        const auto    in = rtepack::meridional_basis(k_in, x), out = rtepack::meridional_basis(k_out, x);

        const auto dielectric = [&](Field E) {
          return scaled_sum(Rh * project(E, s), s, Rv * project(E, p_in), p_out);
        };
        const auto mirror = [&](Field E) {
          const Complex En = 2.0 * project(E, n);
          return Field{En * n[0] - E[0], En * n[1] - E[1], En * n[2] - E[2]};
        };
        const Numeric d =
            deviation(rtepack::fresnel_reflectance_nonspecular(Rv, Rh, k_in, k_out, n), in, out, dielectric);
        const Numeric c =
            deviation(rtepack::fresnel_reflectance_nonspecular(1.0, -1.0, k_in, k_out, n), in, out, mirror);
        require(d < tol and c < tol, std::format("facet mu_in {} mu_out {} phi {}: {} {}", mu_in, mu_out, phi, d, c));

        // Not blind to V: the mirror with V negated once more, as ARTS's surface RT once did, fails
        const rtepack::muelmat flip{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1};
        const Numeric          f =
            deviation(flip * rtepack::fresnel_reflectance_nonspecular(1.0, -1.0, k_in, k_out, n), in, out, mirror);
        require(f > 1e-2, std::format("facet mu_in {} mu_out {} phi {}: a negated V must fail", mu_in, mu_out, phi));
      }
    }
  }

  // rtepack::reflection keeps the reflected V of the matrix: a mirror at nadir reverses it
  const rtepack::stokvec I{1.0, 0.1, 0.2, 0.3};
  const auto             R = rtepack::fresnel_reflectance(1.0, -1.0);
  const auto             r = rtepack::reflection(I, R, rtepack::stokvec{0.0, 0.0, 0.0, 0.0});
  require(r[0] == 1.0 and r[1] == 0.1 and r[2] == -0.2 and r[3] == -0.3,
          std::format("A mirror at normal incidence reverses U and V, got {}", r));
}
}  // namespace

int main() try {
  test_meridional_basis();
  test_basis_rotation();
  test_jones();
  test_lab_frame_phase_matrix();
  test_facet_reflection();
  std::cout << "rtepack polarization test passed\n";
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
