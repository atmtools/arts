/**
 * @file    fluxes.h
 * @author Manfred Brath <manfred.brath@uni-hamburg.de>
 * @date   Tue Sep  26 14:30:0 2023,
 *
 * @brief  Functions for flux calculations.
 * 
 */

#pragma once

#include "matpack_data.h"

void FluxDivergenceFromIrradiance(Tensor3& flux_divergence,
                                 const Vector& p_grid,
                                 const Tensor4& irradiance_field);