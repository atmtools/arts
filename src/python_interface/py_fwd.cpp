#include <arts_omp.h>
#include <debug.h>
#include <fwd.h>
#include <python_interface.h>

#include "hpy_arts.h"

namespace Python {
void py_fwd(py::module_& m) try {
  py::class_<SpectralRadianceOperator> sro(m, "SpectralRadianceOperator");
  generic_interface(sro);
  sro.def(
         "geometric_planar",
         [](const SpectralRadianceOperator& srad_op,
            const Numeric frequency,
            const Vector3 pos,
            const Vector2 los) {
           const auto path = srad_op.geometric_planar(pos, los);
           return srad_op(frequency, path);
         },
         "freq"_a,
         "pos"_a,
         "los"_a,
         "Geometric planar spectral radiance")
      .def(
          "geometric_planar",
          [](const SpectralRadianceOperator& srad_op,
             const Vector& frequency,
             const Vector3 pos,
             const Vector2 los) {
            const auto path = srad_op.geometric_planar(pos, los);

            StokvecVector out(frequency.size());

            if (arts_omp_in_parallel() or arts_omp_get_max_threads() == 1 or
                static_cast<Index>(frequency.size()) < arts_omp_get_max_threads()) {
              std::transform(
                  frequency.begin(),
                  frequency.end(),
                  out.begin(),
                  [&](const Numeric& f) { return srad_op(f, path); });
            } else {
              String error{};
#pragma omp parallel for
              for (Size i = 0; i < out.size(); ++i) {
                try {
                  out[i] = srad_op(frequency[i], path);
                } catch (std::exception& e) {
#pragma omp critical
                  error += e.what() + String{"\n"};
                }
              }

              ARTS_USER_ERROR_IF(not error.empty(), "{}", error)
            }
            return out;
          },
          "freq"_a,
          "pos"_a,
          "los"_a,
          "Geometric planar spectral radiance")
      .def_prop_ro("altitude",
                   &SpectralRadianceOperator::altitude,
                   "The altitude of the top of the atmosphere [m]");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize fwd\n{}", e.what()));
}
}  // namespace Python