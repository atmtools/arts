/**
 * @file    fluxes.h
 * @author Manfred Brath <manfred.brath@uni-hamburg.de>
 * @date   Tue Sep  26 14:30:0 2023,
 *
 * @brief  Functions for flux calculations.
 * 
 */

#include "matpack_data.h"



void FluxDivergenceFromIrradiance(Tensor3& flux_divergence,
                                 const Vector& p_grid,
                                 const Tensor4& irradiance_field) {
  //allocate
  flux_divergence.resize(irradiance_field.nbooks(),
                       irradiance_field.npages(),
                       irradiance_field.nrows());
  flux_divergence = 0;

  // allocate some auxiliary variables
  Numeric net_flux_b;  //net_flux bottom
  Numeric net_flux_c;  //net_flux center
  Numeric net_flux_t;  //net_flux top
  Index idx;

  // calculate flux divergence, we skip the upper and lower boundary here, because to achieve the same
  // second order accuracy for the derivation of the net flux at the edged, we use
  // a differentiation based on polynomial interpolation
  for (Index b = 1; b < irradiance_field.nbooks() - 1; b++) {
    for (Index p = 0; p < irradiance_field.npages(); p++) {
      for (Index r = 0; r < irradiance_field.nrows(); r++) {
        net_flux_b = (irradiance_field(b - 1, p, r, 0) +
                      irradiance_field(b - 1, p, r, 1));
        net_flux_t = (irradiance_field(b + 1, p, r, 0) +
                      irradiance_field(b + 1, p, r, 1));

        flux_divergence(b, p, r) = (net_flux_t - net_flux_b) /
                                 (p_grid[b + 1] - p_grid[b - 1]);
      }
    }
  }

  idx = irradiance_field.nbooks();

  // now calculate the heating rates for the upper and lower boundary
  for (Index p = 0; p < irradiance_field.npages(); p++) {
    for (Index r = 0; r < irradiance_field.nrows(); r++) {
      // lower boundary
      net_flux_b =
          (irradiance_field(0, p, r, 0) + irradiance_field(0, p, r, 1));
      net_flux_c =
          (irradiance_field(1, p, r, 0) + irradiance_field(1, p, r, 1));
      net_flux_t =
          (irradiance_field(2, p, r, 0) + irradiance_field(0, p, r, 1));

      flux_divergence(0, p, r) = (-3 * net_flux_b + 4 * net_flux_c - net_flux_t) /
                               (p_grid[2] - p_grid[0]);

      // upper boundary
      net_flux_t = (irradiance_field(idx - 1, p, r, 0) +
                    irradiance_field(idx - 1, p, r, 1));
      net_flux_c = (irradiance_field(idx - 2, p, r, 0) +
                    irradiance_field(idx - 2, p, r, 1));
      net_flux_b = (irradiance_field(idx - 3, p, r, 0) +
                    irradiance_field(idx - 3, p, r, 1));

      flux_divergence(idx - 1, p, r) =
          -(-3 * net_flux_t + 4 * net_flux_c - net_flux_b) /
          (p_grid[idx-1] - p_grid[idx-3]);
    }
  }
}