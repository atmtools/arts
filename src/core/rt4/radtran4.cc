#include "radtran4.h"

#include <debug.h>
#include <physics_funcs.h>

#include <algorithm>
#include <array>
#include <cmath>

#include "radintg4.h"
#include "radutil4.h"

namespace rt4 {
void radtrano(Numeric          max_delta_tau,
              quadrature_type  quad_type,
              ConstTensor4View surf_reflect,
              ConstMatrixView  gnd_radiance,
              Numeric          sky_temp,
              Numeric          frequency,
              ConstVectorView  height,
              ConstVectorView  temperatures,
              VectorView       gas_extinct,
              ConstVectorView  scatlayers,
              ConstTensor5View extinct_matrix,
              ConstTensor4View emis_vector,
              ConstTensor6View scatter_matrix,
              ConstVectorView  extra_mu,
              VectorView       mu_values,
              Tensor3View      up_rad,
              Tensor3View      down_rad,
              rt4_workdata&    work) {
  // NSTOKES, NUMMU, NUUMMU, NUM_LAYERS and NSL
  const Index nstokes    = up_rad.extent(2);
  const Index nummu      = mu_values.extent(0);
  const Index nuummu     = extra_mu.extent(0);
  const Index num_layers = height.extent(0) - 1;
  const Index nsl        = extinct_matrix.extent(0);

  // The array sizes of radtran4.f
  constexpr Index   maxv = 64, maxm = 4096, maxlay = 400, maxlm = 301 * 4096;
  constexpr Numeric zero = 0.0;

  // The Fortran trusts its declared extents; here they are checked
  ARTS_USER_ERROR_IF(nstokes < 1 or nuummu >= nummu,
                     "RADTRANO needs at least one Stokes parameter (up_rad's last extent) and one quadrature node "
                     "(mu_values longer than extra_mu); got {} Stokes parameters, {} mu_values and {} extra_mu",
                     nstokes,
                     nummu,
                     nuummu);
  ARTS_USER_ERROR_IF(num_layers < 0, "RADTRANO needs at least one height");
  ARTS_USER_ERROR_IF(temperatures.extent(0) != num_layers + 1 or gas_extinct.extent(0) != num_layers or
                         scatlayers.extent(0) != num_layers,
                     "RADTRANO with {} heights needs as many temperatures and one gas_extinct and scatlayers value "
                     "per layer; got {}, {} and {}",
                     height.size(),
                     temperatures.size(),
                     gas_extinct.size(),
                     scatlayers.size());
  ARTS_USER_ERROR_IF(extinct_matrix.shape() != (std::array<Index, 5>{nsl, 2, nummu, nstokes, nstokes}) or
                         emis_vector.shape() != (std::array<Index, 4>{nsl, 2, nummu, nstokes}) or
                         scatter_matrix.shape() != (std::array<Index, 6>{nsl, 4, nummu, nstokes, nummu, nstokes}),
                     "RADTRANO needs extinct_matrix [nsl, 2, nummu, nstokes, nstokes], emis_vector "
                     "[nsl, 2, nummu, nstokes] and scatter_matrix [nsl, 4, nummu, nstokes, nummu, nstokes]; got "
                     "{:B,}, {:B,} and {:B,}",
                     extinct_matrix.shape(),
                     emis_vector.shape(),
                     scatter_matrix.shape());
  ARTS_USER_ERROR_IF(gnd_radiance.shape() != (std::array<Index, 2>{nummu, nstokes}) or mu_values.extent(0) != nummu or
                         up_rad.shape() != (std::array<Index, 3>{num_layers + 1, nummu, nstokes}) or
                         down_rad.shape() != (std::array<Index, 3>{num_layers + 1, nummu, nstokes}),
                     "RADTRANO needs gnd_radiance [nummu, nstokes], mu_values [nummu] and up_rad and down_rad "
                     "[num_layers + 1, nummu, nstokes]; got {:B,}, {}, {:B,} and {:B,}",
                     gnd_radiance.shape(),
                     mu_values.size(),
                     up_rad.shape(),
                     down_rad.shape());
  ARTS_USER_ERROR_IF(surf_reflect.shape() != (std::array<Index, 4>{nummu, nstokes, nummu, nstokes}),
                     "RADTRANO needs surf_reflect [nummu, nstokes, nummu, nstokes]; got {:B,}",
                     surf_reflect.shape());

  // The layers' Planck function is ARTS's planck(), which is negative
  // below 0 K; PLANCK_FUNCTION gave 0 there
  ARTS_USER_ERROR_IF(not(frequency > 0.0), "RADTRANO needs a positive frequency, got {} Hz", frequency);
  ARTS_USER_ERROR_IF(stdr::any_of(temperatures, [](Numeric t) { return not(t >= 0.0); }),
                     "RADTRANO needs temperatures >= 0 K");

  // this is dangerous to do, as we use nstokes also as array-size
  // determining parameter. resetting it will cause problems later on
  // shaping the arrays and correct value extraction. so, don't do
  // this. just check whether (and fail if) condition is not fulfilled.
  //   NSTOKES = MIN(NSTOKES,2)
  ARTS_USER_ERROR_IF(nstokes > 2, "Number of Stokes parameters exceeded.  Maximum size : {}.  Yours is {}", 2, nstokes);

  const bool  symmetric = true;
  const Index n         = nstokes * nummu;
  ARTS_USER_ERROR_IF(n > maxv, "Vector size exceeded.  Maximum size : {}.  Yours is {}", maxv, n);
  ARTS_USER_ERROR_IF(n * n > maxm, "Matrix size exceeded.  Maximum size : {}.  Yours is {}*{} = {}", maxm, n, n, n * n);
  ARTS_USER_ERROR_IF(
      num_layers > maxlay, "Number of layers exceeded.  Maximum number : {}.  Yours is {}", maxlay, num_layers);
  ARTS_USER_ERROR_IF((num_layers + 1) * n * n > maxlm,
                     "Matrix layer size exceeded.  Maximum number (num_layers+1)*(nstokes*nummu)^2: {}.  Yours is "
                     "({}+1)*({}*{})^2 = {}",
                     maxlm,
                     num_layers,
                     nstokes,
                     nummu,
                     (num_layers + 1) * n * n);

  /* RADTRANO's work arrays, as the subroutines read them: those of the
     work data (rt4_workdata), sized for this problem.  The reflection and
     transmission of a slab are [2, n, n], the column-major n x n matrices
     of the + and - directions, a source or a radiance is [2, n].  The
     layers' REFLECT(KRT), TRANS(KRT) and SOURCE(KS), KRT = 1 + 2*N*N*(L-1)
     and KS = 1 + 2*N*(L-1), are reflect[L-1], trans[L-1] and source[L-1];
     L = NUM_LAYERS+1 is the surface. */
  work.resize(nstokes, nummu, num_layers);
  Vector&  quad_weights = work.quad_weights;
  Matrix&  lin_source   = work.lin_source;
  Tensor3& reflect1     = work.reflect1;
  Tensor3& upreflect    = work.upreflect;
  Tensor3& downreflect  = work.downreflect;
  Tensor3& trans1       = work.trans1;
  Tensor3& uptrans      = work.uptrans;
  Tensor3& downtrans    = work.downtrans;
  Matrix&  source1      = work.source1;
  Matrix&  upsource     = work.upsource;
  Matrix&  downsource   = work.downsource;
  Tensor4& reflect      = work.reflect;
  Tensor4& trans        = work.trans;
  Tensor3& source       = work.source;
  Matrix&  sky_radiance = work.sky_radiance;

  // The radiances are in SI, W m-2 Hz-1 sr-1: the Planck function is
  // ARTS's planck() at the frequency

  // Make the desired quadrature abscissas and weights (ARTS's quadratures,
  // rt4::get_quadrature); the extra angles follow them
  const Index j               = nummu - nuummu;
  const auto  q               = get_quadrature(j, quad_type);
  mu_values[Range{0, j}]      = q.mu;
  quad_weights[Range{0, j}]   = q.weights;
  mu_values[Range{j, nuummu}] = extra_mu;

  for (Index i = nummu - 1; i >= j; i--) quad_weights[i] = 0.0;

  // ------------------------------------------------------
  // Loop through the layers
  //   Do doubling to make the reflection and transmission matrices
  //   and soure vectors for each layer, which are stored.

  // jm: skip this check as it has no consequences at all anyways (at least
  //     since after we switched off the info message in CHECK_NORM).
  // !!! the calling code, i.e. ARTS has to handle this issue reliably. !!!
  //   IF (NSL .GT. 0) CALL CHECK_NORM (...)

  for (Index layer = 0; layer < num_layers; layer++) {
    // Calculate the layer thickness
    const Numeric zdiff = std::abs(height[layer] - height[layer + 1]);
    gas_extinct[layer]  = std::max(gas_extinct[layer], 0.0);

    // Do the stuff for thermal source in layer
    // Calculate the thermal source for end of layer
    const Numeric planck1 = planck(frequency, temperatures[layer + 1]);
    // Calculate the thermal source for beginning of layer
    const Numeric planck0 = planck(frequency, temperatures[layer]);

    const Index tsl = std::lround(scatlayers[layer]);
    if (tsl < 1) {
      // If the layer is purely absorbing then quickly
      // make the reflection and transmission matrices
      // and source vector instead of doubling.
      nonscatter_layer(zdiff * gas_extinct[layer],
                       mu_values,
                       planck0,
                       planck1,
                       reflect[layer].view_as(2, nummu, nstokes, nummu, nstokes),
                       trans[layer].view_as(2, nummu, nstokes, nummu, nstokes),
                       source[layer].view_as(2, nummu, nstokes));
    } else {
      ARTS_USER_ERROR_IF(tsl > nsl, "RADTRANO: scatlayers[{}] = {} but there are {} optics sets", layer, tsl, nsl);

      // Find initial thickness of sublayer and
      // the number of times to double
      const Numeric extinct     = extinct_matrix[tsl - 1, 0, 0, 0, 0] + gas_extinct[layer];
      const Numeric f           = std::log(std::max(extinct * zdiff, 1.0e-7) / max_delta_tau) / std::log(2.0);
      Index         num_doubles = 0;
      if (f > 0.0) num_doubles = static_cast<Index>(f) + 1;
      const Numeric num_sub_layers = std::pow(2.0, num_doubles);
      const Numeric delta_z        = zdiff / num_sub_layers;

      // Initialize the source vector
      initial_source(
          delta_z, mu_values, planck0, emis_vector[tsl - 1], gas_extinct[layer], lin_source.view_as(2, nummu, nstokes));
      Numeric linfactor;
      if (planck0 == 0.0) {
        linfactor = 0.0;
      } else {
        linfactor = (planck1 / planck0 - 1.0) / num_sub_layers;
      }

      // Generate the local reflection and transmission matrices
      initialize(delta_z,
                 mu_values,
                 quad_weights,
                 gas_extinct[layer],
                 extinct_matrix[tsl - 1],
                 scatter_matrix[tsl - 1],
                 reflect1.view_as(2, nummu, nstokes, nummu, nstokes),
                 trans1.view_as(2, nummu, nstokes, nummu, nstokes));

      // Double up to the thickness of the layer
      doubling_integration(num_doubles,
                           symmetric,
                           reflect1,
                           trans1,
                           lin_source,
                           linfactor,
                           reflect[layer],
                           trans[layer],
                           source[layer],
                           work);
    }
  }
  // End of layer loop

  // Get the surface reflection and transmission matrices.  The ground's
  // reflection and radiance are external data (rt4::ground_surface makes
  // them for each kind of ground), so only EXTERNAL_SURFACE is left of
  // the ground types; the radiance goes to INTERNAL_RADIANCE.
  const Index ground = num_layers;
  external_surface_layer(surf_reflect,
                         reflect[ground].view_as(2, nummu, nstokes, nummu, nstokes),
                         trans[ground].view_as(2, nummu, nstokes, nummu, nstokes),
                         source[ground].view_as(2, nummu, nstokes));

  // Assume the radiation coming from above is blackbody radiation
  thermal_radiance(sky_temp, zero, frequency, sky_radiance.view_as(2, nummu, nstokes));

  // For each desired output level (1 thru NL+2) add layers
  // above and below level and compute internal radiance.
  // OUTLEVELS gives the desired output levels.
  for (Index i = 0; i < num_layers + 1; i++) {
    const Index layer = std::min(std::max(i, Index{0}), num_layers + 1);
    upreflect         = 0.0;
    downreflect       = 0.0;
    matpack::identity(uptrans[0]);
    matpack::identity(uptrans[1]);
    matpack::identity(downtrans[0]);
    matpack::identity(downtrans[1]);
    upsource   = 0.0;
    downsource = 0.0;
    for (Index l = 0; l < layer; l++) {
      if (l == 0) {
        upreflect = reflect[l];
        uptrans   = trans[l];
        upsource  = source[l];
      } else {
        reflect1 = upreflect;
        trans1   = uptrans;
        source1  = upsource;
        combine_layers(reflect1, trans1, source1, reflect[l], trans[l], source[l], upreflect, uptrans, upsource, work);
      }
    }
    for (Index l = layer; l < num_layers + 1; l++) {
      if (l == layer) {
        downreflect = reflect[l];
        downtrans   = trans[l];
        downsource  = source[l];
      } else {
        reflect1 = downreflect;
        trans1   = downtrans;
        source1  = downsource;
        combine_layers(
            reflect1, trans1, source1, reflect[l], trans[l], source[l], downreflect, downtrans, downsource, work);
      }
    }
    internal_radiance(upreflect,
                      uptrans,
                      upsource,
                      downreflect,
                      downtrans,
                      downsource,
                      sky_radiance[0],
                      gnd_radiance.view_as(n),
                      up_rad[i].view_as(n),
                      down_rad[i].view_as(n),
                      work);
  }

  // Integrate the mu times the radiance to find the fluxes
  // jm: we don't care about the fluxes so far. so, just skip this.
}
}  // namespace rt4
