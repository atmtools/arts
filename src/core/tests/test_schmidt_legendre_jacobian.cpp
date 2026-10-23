#include <debug.h>
#include <functional_atm_mag_field.h>
#include <planet_data.h>

#include <algorithm>
#include <cmath>
#include <print>

/** Compares the derivatives of the Legendre field with respect to its coefficients with finite differences
 *
 * The field is linear in the coefficients, so the central differences are exact to rounding, and every
 * kept coefficient must change the field.
 */
int main() {
  // Up to degree 13, where the east component, u, does not depend on the m = 0 coefficients
  constexpr Size N = 14;

  for (const Vector3 pos : {Vector3{90e3, 35.0, 139.0}, Vector3{0.0, -60.0, -20.0}, Vector3{500e3, 10.0, 75.0}}) {
    for (const auto component : {FieldComponent::u, FieldComponent::v, FieldComponent::w}) {
      const Size ncoef = component == FieldComponent::u ? N * (N - 1) : N * N - 1;

      // Arbitrary non-zero coefficients [T]
      Vector gh(ncoef);
      for (Size i = 0; i < gh.size(); i++) gh[i] = (i % 2 ? -1e-5 : 1e-5) / static_cast<Numeric>(1 + i);

      Atm::SchmidthLegendre f{
          .gh = gh, .N = N, .r0 = 6371.2e3, .ell = {Body::Earth::a, Body::Earth::b}, .component = component};

      const auto w = f.w(pos[0], pos[1], pos[2]);

      if (w.size() != ncoef or static_cast<Size>(f.x().size()) != ncoef)
        arts::user_error("Expected {} derivatives and coefficients, got {} and {}", ncoef, w.size(), f.x().size());

      Numeric max_w = 0.0, max_error = 0.0;
      for (const auto& [i, dfdx] : w) {
        const Numeric x0 = f.x()[i];
        const Numeric dx = 1e-9;

        f.x()[i]         = x0 + dx;
        const Numeric fp = f(pos[0], pos[1], pos[2]);
        f.x()[i]         = x0 - dx;
        const Numeric fm = f(pos[0], pos[1], pos[2]);
        f.x()[i]         = x0;

        if (dfdx == 0.0) arts::user_error("The derivative of coefficient {} is zero, so it should not be kept", i);

        max_error = std::max(max_error, std::abs((fp - fm) / (2 * dx) - dfdx));
        max_w     = std::max(max_w, std::abs(dfdx));
      }

      std::println("pos {:B,}, component {}: {} derivatives, max |derivative| {:.3e}, max |error| {:.3e}",
                   pos,
                   component,
                   w.size(),
                   max_w,
                   max_error);

      if (max_error > 1e-9 * max_w) arts::user_error("The derivatives do not match the finite differences");
    }
  }

  return 0;
}
