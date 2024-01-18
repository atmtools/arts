/**
 * @file radiation_field.cc
 * @author Richard Larsson
 * @date 2019-09-04
 * 
 * @brief Radiation field calculations
 */

#include "radiation_field.h"
#include "arts_conversions.h"
#include "sorting.h"

void error_in_integrate(const String& error_msg,
                        const Numeric& value_that_should_be_unity) {
  ARTS_USER_ERROR_IF (std::abs(value_that_should_be_unity - 1.0) > 1e-4,
    "Failure in normalization:\n", error_msg, "\n")
}

Numeric test_integrate_convolved(const Eigen::Ref<Eigen::VectorXcd> &F,
                                 const Vector& f) {
  Numeric val = 0.0;

  const Index n = f.size();
  for (Index i = 0; i < n - 1; i++)
    val += 0.5 * (f[i + 1] - f[i]) * (F[i].real() + F[i + 1].real());

  return val;  // Should return 1.0 for normalized F
}

Numeric test_integrate_zenith(const Vector& cosza,
                              const Array<Index>& sorted_index) {
  Numeric val = 0.0;

  const Index n = cosza.size();
  for (Index i = 0; i < n - 1; i++)
    val += 0.5 * (cosza[sorted_index[i]] - cosza[sorted_index[i + 1]]);

  return val;  // Should return 1.0 for normalized cosza
}

Numeric integrate_convolved(const StokvecVector& I,
                            const Eigen::VectorXcd& F,
                            const Vector& f) {
  Numeric val = 0.0;

  const Index n = f.size();
  for (Index i = 0; i < n - 1; i++)
    val += 0.5 * (f[i + 1] - f[i]) *
           (I[i].I() * F[i].real() + I[i + 1].I() * F[i + 1].real());

  return val;
}

Numeric integrate_convolved(const MuelmatVector& T,
                            const Eigen::VectorXcd& F,
                            const Vector& f) {
  Numeric val = 0.0;

  const Index n = f.size();
  for (Index i = 0; i < n - 1; i++)
    val += 0.5 * (f[i + 1] - f[i]) *
           (T[i](0, 0) * F[i].real() + T[i + 1](0, 0) * F[i + 1].real());

  return 1.0 - val;
}

Numeric integrate_zenith(const ConstVectorView& j,
                         const Vector& cosza,
                         const Array<Index>& sorted_index) {
  Numeric val = 0.0;

  const Index n = cosza.size();
  for (Index i = 0; i < n - 1; i++)
    val += 0.25 * (cosza[sorted_index[i]] - cosza[sorted_index[i + 1]]) *
           (j[sorted_index[i]] + j[sorted_index[i + 1]]);

  return val;
}
