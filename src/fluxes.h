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

/** Calculates the derivative of the irradiance_field with respect to pressure.
 *
 * @param[out] flux_divergence The output tensor containing the calculated flux divergence.
 * @param[in] p_grid The vector of grid points at which to calculate the derivative.
 * @param[in] irradiance_field The input tensor containing the irradiance field.
 */
void FluxDivergenceFromIrradiance(Tensor3& flux_divergence,
                                 const Vector& p_grid,
                                 const Tensor4& irradiance_field);