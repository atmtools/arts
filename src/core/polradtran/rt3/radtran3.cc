#include "radtran3.h"

#include <arts_constants.h>
#include <debug.h>
#include <radintg.h>
#include <radutil.h>

#include <algorithm>
#include <array>
#include <cmath>

#include "radintg3.h"
#include "radscat3.h"
#include "radutil3.h"

namespace polradtran::rt3 {
void radtran(Numeric                             max_delta_tau,
             Index                               src_code,
             quadrature_type                     quad_type,
             bool                                delta_m,
             Numeric                             direct_flux,
             Numeric                             direct_mu,
             ConstTensor5View                    surf_reflect,
             ConstTensor3View                    gnd_radiance,
             ConstTensor3View                    direct_reflect,
             Numeric                             sky_temp,
             Numeric                             frequency,
             ConstVectorView                     height,
             ConstVectorView                     temperatures,
             ConstVectorView                     gas_extinct,
             ConstVectorView                     scat_extinct,
             ConstVectorView                     scat_scatter,
             const ArrayOfIndex&                 scat_nlegen,
             CompactPlanarMuelmatConstMatrixView scat_coef,
             const ArrayOfIndex&                 scatlayers,
             const ArrayOfIndex&                 outlevels,
             ConstVectorView                     extra_mu,
             VectorView                          mu_values,
             MatrixView                          up_flux,
             MatrixView                          down_flux,
             Tensor4View                         up_rad,
             Tensor4View                         down_rad,
             rt3_workdata&                       work) {
  // NSTOKES, NUMMU, AZIORDER, NUM_LAYERS, NSL, LDCOEF and NOUTLEVELS
  const Index nstokes    = up_rad.extent(3);
  const Index nummu      = mu_values.extent(0);
  const Index nuummu     = extra_mu.extent(0);
  const Index aziorder   = up_rad.extent(1) - 1;
  const Index num_layers = height.extent(0) - 1;
  const Index nsl        = scat_coef.extent(0);
  const Index ldcoef     = scat_coef.extent(1);
  const Index noutlevels = static_cast<Index>(outlevels.size());

  constexpr Numeric pi = Constant::pi, twopi = Constant::two_pi, zero = 0.0;

  // The Fortran trusts its declared extents; here they are checked
  ARTS_USER_ERROR_IF(nstokes < 1 or nstokes > 4 or aziorder < 0 or nuummu >= nummu,
                     "RADTRAN needs 1 to 4 Stokes parameters and aziorder >= 0 (up_rad's last extents) and at least "
                     "one quadrature node (mu_values longer than extra_mu); got {} Stokes parameters, aziorder {}, "
                     "{} mu_values and {} extra_mu",
                     nstokes,
                     aziorder,
                     nummu,
                     nuummu);
  ARTS_USER_ERROR_IF(num_layers < 0, "RADTRAN needs at least one height");
  ARTS_USER_ERROR_IF(temperatures.extent(0) != num_layers + 1 or gas_extinct.extent(0) != num_layers or
                         static_cast<Index>(scatlayers.size()) != num_layers,
                     "RADTRAN with {} heights needs as many temperatures and one gas_extinct and scatlayers value per "
                     "layer; got {}, {} and {}",
                     height.size(),
                     temperatures.size(),
                     gas_extinct.size(),
                     scatlayers.size());
  ARTS_USER_ERROR_IF(
      scat_extinct.extent(0) != nsl or scat_scatter.extent(0) != nsl or static_cast<Index>(scat_nlegen.size()) != nsl,
      "RADTRAN needs scat_extinct, scat_scatter and scat_nlegen [nsl] and scat_coef [nsl, ldcoef]; "
      "got {}, {}, {} and {:B,}",
      scat_extinct.size(),
      scat_scatter.size(),
      scat_nlegen.size(),
      scat_coef.shape());
  ARTS_USER_ERROR_IF(stdr::any_of(scat_nlegen, [ldcoef](Index l) { return l < 0 or l + 1 > ldcoef; }),
                     "RADTRAN needs every SCAT_NLEGEN in [0, LDCOEF - 1] = [0, {}]",
                     ldcoef - 1);
  ARTS_USER_ERROR_IF(stdr::any_of(scatlayers, [nsl](Index s) { return s < 0 or s > nsl; }),
                     "RADTRAN needs every SCATLAYERS in [0, NSL] = [0, {}]",
                     nsl);
  ARTS_USER_ERROR_IF(stdr::any_of(outlevels, [num_layers](Index l) { return l < 1 or l > num_layers + 1; }),
                     "RADTRAN needs every OUTLEVELS in [1, NUM_LAYERS + 1] = [1, {}]",
                     num_layers + 1);
  ARTS_USER_ERROR_IF(src_code < 0 or src_code > 3, "RADTRAN needs SRC_CODE 0 to 3, got {}", src_code);
  ARTS_USER_ERROR_IF(
      up_rad.shape() != (std::array<Index, 4>{noutlevels, aziorder + 1, nummu, nstokes}) or
          down_rad.shape() != up_rad.shape() or up_flux.shape() != (std::array<Index, 2>{noutlevels, nstokes}) or
          down_flux.shape() != up_flux.shape(),
      "RADTRAN needs up_rad and down_rad [noutlevels, aziorder + 1, nummu, nstokes] and up_flux and down_flux "
      "[noutlevels, nstokes] (noutlevels = {}, nummu = {}); got {:B,}, {:B,}, {:B,} and {:B,}",
      noutlevels,
      nummu,
      up_rad.shape(),
      down_rad.shape(),
      up_flux.shape(),
      down_flux.shape());
  ARTS_USER_ERROR_IF(
      surf_reflect.shape() != (std::array<Index, 5>{aziorder + 1, nummu, nstokes, nummu, nstokes}) or
          gnd_radiance.shape() != (std::array<Index, 3>{aziorder + 1, nummu, nstokes}) or
          direct_reflect.shape() != gnd_radiance.shape(),
      "RADTRAN needs surf_reflect [aziorder + 1, nummu, nstokes, nummu, nstokes] and gnd_radiance and "
      "direct_reflect [aziorder + 1, nummu, nstokes] (aziorder = {}, nummu = {}, nstokes = {}); got {:B,}, {:B,} "
      "and {:B,}",
      aziorder,
      nummu,
      nstokes,
      surf_reflect.shape(),
      gnd_radiance.shape(),
      direct_reflect.shape());

  const bool solar   = src_code == 1 or src_code == 3;
  const bool thermal = src_code == 2 or src_code == 3;

  // The radiances are in SI, W m-2 Hz-1 sr-1, at the frequency
  ARTS_USER_ERROR_IF(not(frequency > 0.0), "RADTRAN needs a positive frequency, got {} Hz", frequency);

  const bool  symmetric = nstokes <= 2;
  const Index n         = nstokes * nummu;

  /* RADTRAN's work arrays, as the subroutines read them: those of the work
     data (rt3_workdata), sized for this problem.  The reflection and
     transmission of a slab are [2, n, n], the column-major n x n matrices
     of the + and - directions, a source or a radiance is [2, n].  The
     layers' REFLECT(KRT), TRANS(KRT) and SOURCE(KS), KRT = 1 +
     2*N*N*(L-1) and KS = 1 + 2*N*(L-1), are reflect[L-1], trans[L-1] and
     source[L-1]; L = NUM_LAYERS+1 is the surface.  SCATBUF and DIRECTBUF
     hold the scattering matrices (Muelmat) and direct vectors (Stokvec) of
     every set and mode,
     as SCATTERING and DIRECT_SCATTERING lay them out; scatbuf[SCAT_NUM-1]
     and directbuf[SCAT_NUM-1] are the set's parts (see rt3::scattering and
     rt3::direct_scattering); SCATTER_MATRIX (PHASE_FUNCTION to INITIALIZE)
     is [4, nummu, nstokes, nummu, nstokes]. */
  Index num_legendre = 2 * nummu;
  for (Index l : scat_nlegen) num_legendre = std::max(num_legendre, l + 1);
  work.resize(nstokes, nummu, aziorder, num_layers, nsl, num_legendre);
  Vector&                     quad_weights      = work.quad_weights;
  CompactPlanarMuelmatVector& legendre_coef     = work.legendre_coef;
  Vector&                     set_extinct       = work.set_extinct;
  Vector&                     set_scatter       = work.set_scatter;
  MuelmatTensor5&             scatbuf           = work.scatbuf;
  StokvecTensor4&             directbuf         = work.directbuf;
  ArrayOfIndex&               scat_nums         = work.scat_nums;
  Vector&                     extinctions       = work.extinctions;
  Vector&                     albedos           = work.albedos;
  Vector&                     direct_level_flux = work.direct_level_flux;
  Tensor5&                    scatter_matrix    = work.scatter_matrix;
  Matrix&                     direct_vector     = work.direct_vector;
  Matrix&                     exp_source        = work.exp_source;
  Matrix&                     thermal_vector    = work.thermal_vector;
  Matrix&                     lin_source        = work.lin_source;
  Tensor3&                    reflect1          = work.reflect1;
  Tensor3&                    trans1            = work.trans1;
  Matrix&                     source1           = work.source1;
  Tensor4&                    reflect           = work.reflect;
  Tensor4&                    trans             = work.trans;
  Tensor3&                    source            = work.source;
  Matrix&                     ground_radiance   = work.ground_radiance;
  Matrix&                     direct_radiance   = work.direct_radiance;
  Matrix&                     sky_radiance      = work.sky_radiance;

  // Make the desired quadrature abscissas and weights (ARTS's quadratures,
  // polradtran::get_quadrature); with QUAD_TYPE 'E' the extra angles follow them.
  // NLEGLIM, the highest Legendre degree kept, is max_legendre_degree.
  const Index j                  = nummu - nuummu;
  const auto  q                  = get_quadrature(j, quad_type);
  mu_values[Range{0, j}]         = q.mu;
  quad_weights[Range{0, j}]      = q.weights;
  mu_values[Range{j, nuummu}]    = extra_mu;
  quad_weights[Range{j, nuummu}] = 0.0;
  const Index nleglim            = max_legendre_degree(j, quad_type);

  // Make all of the scattering matrices ahead of time
  // and store them in memory
  for (Index scat_num = 1; scat_num <= nsl; scat_num++) {
    const Index s        = scat_num - 1;
    Index       numlegen = 0;
    get_scat_set(delta_m,
                 nummu,
                 scat_coef[s, Range{0, scat_nlegen[s] + 1}],
                 scat_extinct[s],
                 scat_scatter[s],
                 numlegen,
                 legendre_coef,
                 set_extinct[s],
                 set_scatter[s]);
    // Truncate the Legendre series to enforce normalization (RADTRAN also
    // prints that it does)
    if (numlegen > nleglim) numlegen = nleglim;
    // Make the scattering matrix
    scattering(mu_values, quad_weights, legendre_coef[Range{0, numlegen + 1}], nstokes, scatbuf[s], work);
    // Make the direct (solar) pseudo source
    if (solar)
      direct_scattering(mu_values, legendre_coef[Range{0, numlegen + 1}], direct_mu, nstokes, directbuf[s], work);
  }

  // SCATLAYERS(LAYER) is the set of each layer.  A non-scattering layer (0)
  // keeps the set number of the layer above, or 0 at the top.
  {
    Index scat_num = 0;
    for (Index layer = 0; layer < num_layers; layer++) {
      // Special case for a non-scattering layer
      Numeric extinct, scatter;
      if (scatlayers[layer] == 0) {
        extinct = 0.0;
        scatter = 0.0;
      } else {
        scat_num = scatlayers[layer];
        extinct  = set_extinct[scat_num - 1];
        scatter  = set_scatter[scat_num - 1];
      }

      scat_nums[layer]   = scat_num;
      extinctions[layer] = extinct + std::max(gas_extinct[layer], 0.0);
      if (extinctions[layer] > 0.0) {
        albedos[layer] = scatter / extinctions[layer];
      } else {
        albedos[layer] = 0.0;
      }
    }
  }

  // Compute the direct beam flux at each level
  if (solar) {
    Numeric tau          = 0.0;
    direct_level_flux[0] = direct_flux;
    for (Index layer = 0; layer < num_layers; layer++) {
      tau                          = tau + extinctions[layer] / direct_mu * std::abs(height[layer] - height[layer + 1]);
      direct_level_flux[layer + 1] = direct_flux * std::exp(-tau);
    }
  }

  // Loop through each azimuth mode
  for (Index mode = 0; mode <= aziorder; mode++) {
    Index scat_num = 0;
    // ------------------------------------------------------
    // Loop through the layers
    for (Index layer = 0; layer < num_layers; layer++) {
      // Calculate the layer thickness
      const Numeric zdiff      = std::abs(height[layer] - height[layer + 1]);
      const Numeric extinction = extinctions[layer];
      const Numeric albedo     = albedos[layer];

      if (scat_nums[layer] != scat_num) {
        scat_num = scat_nums[layer];
        // Get the scattering matrix from the buffer
        get_scattering(mode, scatbuf[scat_num - 1], scatter_matrix);
        // Check the normalization of the scattering matrix
        if (mode == 0) check_norm(quad_weights, scatter_matrix);
        // Get the direct (solar) vector from the buffer
        if (solar) get_direct(mode, directbuf[scat_num - 1], direct_vector.view_as(2, nummu, nstokes));
      }

      // Compute the thermal emission at top and bottom of layer
      Numeric planck0, planck1;
      if (thermal) {
        // Calculate the thermal source for end of layer
        thermal_radiance(mode, temperatures[layer + 1], albedo, frequency, thermal_vector.view_as(2, nummu, nstokes));
        planck1 = thermal_vector[0, 0];
        // Calculate the thermal source for beginning of layer
        thermal_radiance(mode, temperatures[layer], albedo, frequency, thermal_vector.view_as(2, nummu, nstokes));
        planck0 = thermal_vector[0, 0];
      } else {
        planck0 = 0.0;
        planck1 = 0.0;
      }

      if (albedo == zero) {
        // If the layer is purely absorbing then quickly
        // make the reflection and transmission matrices
        // and source vector instead of doubling.
        nonscatter_layer(mode,
                         zdiff * extinction,
                         mu_values,
                         planck0,
                         planck1,
                         reflect[layer].view_as(2, nummu, nstokes, nummu, nstokes),
                         trans[layer].view_as(2, nummu, nstokes, nummu, nstokes),
                         source[layer].view_as(2, nummu, nstokes));
      } else {
        // Find initial thickness of sublayer and
        // the number of times to double
        const auto [num_doubles, num_sub_layers, delta_z] = initial_sublayer(zdiff, extinction, max_delta_tau);

        // For a solar source make the pseudo source vector
        // and initialize it
        Numeric expfactor = 0.0, linfactor = 0.0;
        if (solar) {
          const Numeric tmp  = direct_level_flux[layer] * albedo / (4.0 * pi * direct_mu);
          source1            = direct_vector;
          source1           *= tmp;
          initial_source(delta_z,
                         mu_values,
                         extinction,
                         source1.view_as(2, nummu, nstokes),
                         exp_source.view_as(2, nummu, nstokes));
          expfactor = std::exp(-extinction * delta_z / direct_mu);
        }

        // Initialize the thermal source vector
        if (thermal) {
          initial_source(delta_z,
                         mu_values,
                         extinction,
                         thermal_vector.view_as(2, nummu, nstokes),
                         lin_source.view_as(2, nummu, nstokes));
          if (planck0 == 0.0) {
            linfactor = 0.0;
          } else {
            linfactor = (planck1 / planck0 - 1.0) / num_sub_layers;
          }
        }

        // Generate the local reflection and transmission matrices
        initialize(delta_z,
                   mu_values,
                   extinction,
                   albedo,
                   scatter_matrix,
                   reflect1.view_as(2, nummu, nstokes, nummu, nstokes),
                   trans1.view_as(2, nummu, nstokes, nummu, nstokes));

        // Double up to the thickness of the layer
        doubling_integration(num_doubles,
                             src_code,
                             symmetric,
                             reflect1,
                             trans1,
                             exp_source,
                             expfactor,
                             lin_source,
                             linfactor,
                             reflect[layer],
                             trans[layer],
                             source[layer],
                             work);
      }
    }
    // End of layer loop

    // Get the surface reflection and transmission matrices and the surface
    // radiance.  The ground is external data (rt3::ground_surface makes it
    // for each kind of ground): its reflection makes the surface layer, and
    // its radiance, with the direct beam it reflects, goes to
    // INTERNAL_RADIANCE as GND_RADIANCE.
    const Index ground = num_layers;
    external_surface_layer(surf_reflect[mode],
                           reflect[ground].view_as(2, nummu, nstokes, nummu, nstokes),
                           trans[ground].view_as(2, nummu, nstokes, nummu, nstokes),
                           source[ground].view_as(2, nummu, nstokes));
    ground_radiance = gnd_radiance[mode];
    if (solar) {
      direct_radiance  = direct_reflect[mode];
      direct_radiance *= direct_level_flux[num_layers];
      ground_radiance += direct_radiance;
    }

    // Assume the radiation coming from above is blackbody radiation
    thermal_radiance(mode, sky_temp, 0.0, frequency, sky_radiance.view_as(2, nummu, nstokes));

    // For each desired output level (1 thru NL+2) add layers
    // above and below level and compute internal radiance
    for (Index i = 0; i < noutlevels; i++) {
      // The 0-based layer just below the level, the number of layers above it
      const Index layer = std::min(std::max(outlevels[i], Index{1}), num_layers + 2) - 1;
      level_radiance(layer,
                     reflect,
                     trans,
                     source,
                     sky_radiance[0],
                     ground_radiance.view_as(n),
                     up_rad[i, mode].view_as(n),
                     down_rad[i, mode].view_as(n),
                     work);
    }
  }
  // End of azimuth mode loop

  // Integrate mu times the radiance to find the fluxes
  for (Index l = 0; l < noutlevels; l++) {
    for (Index i = 0; i < nstokes; i++) {
      up_flux[l, i]   = 0.0;
      down_flux[l, i] = 0.0;
      for (Index jj = 0; jj < nummu; jj++) {
        up_flux[l, i]   = up_flux[l, i] + twopi * quad_weights[jj] * mu_values[jj] * up_rad[l, 0, jj, i];
        down_flux[l, i] = down_flux[l, i] + twopi * quad_weights[jj] * mu_values[jj] * down_rad[l, 0, jj, i];
      }
    }
    // Add in direct beam fluxes
    if (solar) down_flux[l, 0] = down_flux[l, 0] + direct_level_flux[outlevels[l] - 1];
  }
}
}  // namespace polradtran::rt3
