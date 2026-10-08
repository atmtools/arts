#include "radutil.h"

#include <arts_conversions.h>
#include <debug.h>
#include <physics_funcs.h>
#include <rtepack.h>

#include <array>

namespace polradtran {
void lambert_surface_layer(Index           mode,
                           ConstVectorView mu_values,
                           ConstVectorView quad_weights,
                           Numeric         ground_albedo,
                           Tensor5View     reflect,
                           Tensor5View     trans,
                           Tensor3View     source) {
  const Index nummu   = mu_values.extent(0);
  const Index nstokes = source.extent(2);
  ARTS_USER_ERROR_IF(quad_weights.extent(0) != nummu or
                         reflect.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         trans.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         source.shape() != (std::array<Index, 3>{2, nummu, nstokes}),
                     "LAMBERT_SURFACE with {} mu_values needs quad_weights [nummu], reflect and trans [2, nummu, "
                     "nstokes, nummu, nstokes] and source [2, nummu, nstokes]; got {}, {:B,}, {:B,} and {:B,}",
                     nummu,
                     quad_weights.extent(0),
                     reflect.shape(),
                     trans.shape(),
                     source.shape());

  reflect = 0.0;
  source  = 0.0;
  for (Index h = 0; h < 2; h++) identity(trans[h].view_as(nummu * nstokes, nummu * nstokes));
  // The Lambertian ground reflects the flux equally in all direction
  // and completely unpolarizes the radiation
  if (mode == 0) {
    for (Index j2 = 0; j2 < nummu; j2++)
      reflect[1, j2, 0, joker, 0] = 2.0 * ground_albedo * mu_values[j2] * quad_weights[j2];
  }
}

void fresnel_surface_layer(
    ConstVectorView mu_values, Complex index, Tensor5View reflect, Tensor5View trans, Tensor3View source) {
  const Index nummu   = mu_values.extent(0);
  const Index nstokes = source.extent(2);
  ARTS_USER_ERROR_IF(nstokes > 4 or reflect.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         trans.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         source.shape() != (std::array<Index, 3>{2, nummu, nstokes}),
                     "FRESNEL_SURFACE with {} mu_values needs at most 4 Stokes components and reflect and trans "
                     "[2, nummu, nstokes, nummu, nstokes] and source [2, nummu, nstokes]; got {:B,}, {:B,} and {:B,}",
                     nummu,
                     reflect.shape(),
                     trans.shape(),
                     source.shape());

  reflect = 0.0;
  source  = 0.0;
  for (Index h = 0; h < 2; h++) identity(trans[h].view_as(nummu * nstokes, nummu * nstokes));

  // REFLECT(I, J, K, J, 2) = R(I, K): R1 and R2 in the [I, Q] block, R3 and
  // R4 in the [U, V] block
  const Range stokes{0, nstokes};
  for (Index j = 0; j < nummu; j++) {
    const auto [rv, rh]            = fresnel(1.0, index, Conversion::acosd(mu_values[j]));
    const Muelmat r                = rtepack::fresnel_reflectance(rv, rh);
    reflect[1, j, joker, j, joker] = transpose(r.view()[stokes, stokes]);
  }
}

void external_surface_layer(ConstTensor4View surf_reflect, Tensor5View reflect, Tensor5View trans, Tensor3View source) {
  const Index nummu   = surf_reflect.extent(0);
  const Index nstokes = surf_reflect.extent(1);
  ARTS_USER_ERROR_IF(surf_reflect.shape() != (std::array<Index, 4>{nummu, nstokes, nummu, nstokes}) or
                         reflect.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         trans.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         source.shape() != (std::array<Index, 3>{2, nummu, nstokes}),
                     "EXTERNAL_SURFACE needs surf_reflect [nummu, nstokes, nummu, nstokes], reflect and trans [2, "
                     "nummu, nstokes, nummu, nstokes] and source [2, nummu, nstokes]; got {:B,}, {:B,}, {:B,} and "
                     "{:B,}",
                     surf_reflect.shape(),
                     reflect.shape(),
                     trans.shape(),
                     source.shape());

  reflect = 0.0;
  source  = 0.0;
  for (Index h = 0; h < 2; h++) identity(trans[h].view_as(nummu * nstokes, nummu * nstokes));

  // REFLECT(I1, J1, I2, J2, 2) = SURF_REFL(I1, J1, I2, J2)
  reflect[1] = surf_reflect;
}
}  // namespace polradtran
