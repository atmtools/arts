/* Copyright (C) 2015
   Richard Larsson <ric.larsson@gmail.com>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

#include "linescaling.h"
#include "interpolation_poly.h"

/*!
 *  Calculates the line strength scaling parameters for cross section calculations.
 * 
 *  The only atmospheric input is the temperature.  The line knows its energy and
 *  its reference temperature.  If a custom PF tag was applied, take that route,
 *  otherwise the partition function defaults to inbuilt partition function data.
 *  If atm_tv* are non-negative, then the non-LTE population levels are calculated.
 *  (This only works for rotational LTE at this point in time.)  The non-LTE implemnation
 *  has heritage from the FUTBOLIN implementation.
 *  
 * 
 *  \param  partition_ratio      Out:    The partition ratio to atmospheric temperature (LTE)
 *  \param  boltzmann_ratio      Out:    The boltzmann ratio to atmospheric temperature (LTE)
 *  \param  abs_nlte_ratio       Out:    The relative extra absorption due to NLTE effects
 *  \param  src_nlte_ratio       Out:    The relative extra source due to NLTE effects
 *  \param  partition_type       In:     Switch for partition type of line
 *  \param  partition_data       In:     Switch for partition data of line
 *  \param  atm_t                In:     The path point atmospheric temperature
 *  \param  line_t               In:     The line reference temperature
 *  \param  line_f               In:     The line central frequency
 *  \param  line_elow            In:     The line lower energy level
 *  \param  do_nlte              In:     Bool for "We need to to NLTE calculations"
 *  \param  line_evlow           In:     The line lower vibrational energy level
 *  \param  line_evupp           In:     The line upper vibrational energy level
 *  \param  line_evlow_index     In:     The line lower vibrational energy level index
 *  \param  line_evupp_index     In:     The line upper vibrational energy level index
 *  \param  atm_t_nlte           In:     Vector of NLTE temperatures.  The line knows which ones belong to it.
 * 
 *  \author Richard Larsson
 *  \date   2015-05-28
 */
void GetLineScalingData(Numeric& q_t,
                        Numeric& q_ref,
                        Numeric& partition_ratio,
                        Numeric& K1,
                        Numeric& K2,
                        Numeric& abs_nlte_ratio,
                        Numeric& src_nlte_ratio,
                        const SpeciesAuxData::AuxType& partition_type,
                        const ArrayOfGriddedField1& partition_data,
                        const Numeric& atm_t,
                        const Numeric& line_t,
                        const Numeric& line_f,
                        const Numeric& line_elow,
                        const bool& do_nlte,
                        const Numeric& line_evlow,
                        const Numeric& line_evupp,
                        const Index& line_evlow_index,
                        const Index& line_evupp_index,
                        ConstVectorView atm_t_nlte) {
  // Physical constants
  extern const Numeric PLANCK_CONST;
  extern const Numeric BOLTZMAN_CONST;

  if (q_t < 0 || q_ref < 0) {
    partition_function(
        q_ref, q_t, line_t, atm_t, partition_type, partition_data);

    partition_ratio = q_ref / q_t;
  }

  // Following Futbolin's division into two parts for the Boltzmann ratio because
  // gamma is also used for the NLTE part later on
  const Numeric gamma = exp(-PLANCK_CONST * line_f / (BOLTZMAN_CONST * atm_t));
  const Numeric gamma_ref =
      exp(-PLANCK_CONST * line_f / (BOLTZMAN_CONST * line_t));

  // Stimulated emission
  K2 = (1. - gamma) / (1. - gamma_ref);

  // Boltzmann level
  K1 = exp(line_elow / BOLTZMAN_CONST * (atm_t - line_t) / (atm_t * line_t));

  if (do_nlte) {
    // Test the NLTE of the line and find if we should use atmospheric temperatures or not
    const Numeric& atm_tv_low =
        line_evlow_index < 0 ? -1.0 : atm_t_nlte[line_evlow_index];
    const Numeric& atm_tv_upp =
        line_evupp_index < 0 ? -1.0 : atm_t_nlte[line_evupp_index];

    //r_low and r_upp are ratios for the population level compared to LTE conditions
    Numeric r_low, r_upp;
    if (atm_tv_low >
        0.0)  // where 1e-4 is considered a small number so that the multiplication in the denominator does not reach zero
      r_low = exp(-line_evlow / BOLTZMAN_CONST * (atm_t - atm_tv_low) /
                  (atm_t * atm_tv_low));
    else if (atm_tv_low == 0.0)
      throw std::runtime_error(
          "A line has been defined with zero vibrational temperature.\nThis is not physical.\n");
    else
      r_low = 1.0;

    if (atm_tv_upp >
        0.0)  // where 1e-4 is considered a small number so that the multiplication in the denominator does not reach zero
      r_upp = exp(-line_evupp / BOLTZMAN_CONST * (atm_t - atm_tv_upp) /
                  (atm_t * atm_tv_upp));
    else if (atm_tv_upp == 0.0)
      throw std::runtime_error(
          "A line has been defined with zero vibrational temperature.\nThis is not physical.\n");
    else
      r_upp = 1.0;

    // Both are unity when in LTE
    abs_nlte_ratio = (r_low - r_upp * gamma) / (1 - gamma);
    src_nlte_ratio = r_upp;
  }
}
void GetLineScalingData_dT(Numeric& dq_t_dT,
                           Numeric& dK2_dT,
                           Numeric& dpartition_ratio_dT,
                           Numeric& dabs_nlte_ratio_dT,
                           Numeric& atm_tv_low,
                           Numeric& atm_tv_upp,
                           const Numeric& q_t,
                           const Numeric& abs_nlte_ratio,
                           const SpeciesAuxData::AuxType& partition_type,
                           const ArrayOfGriddedField1& partition_data,
                           const Numeric& atm_t,
                           const Numeric& line_t,
                           const Numeric& dt,
                           const Numeric& line_f,
                           const bool& do_nlte,
                           const Numeric& line_evlow,
                           const Numeric& line_evupp,
                           const Index& line_evlow_index,
                           const Index& line_evupp_index,
                           ConstVectorView atm_t_nlte) {
  // Physical constants
  extern const Numeric PLANCK_CONST;
  extern const Numeric BOLTZMAN_CONST;

  if (dq_t_dT < 0) {
    switch (partition_type) {
      case SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF:
        CalculatePartitionFctFromCoeff_dT(
            dq_t_dT, atm_t, partition_data[0].data);
        break;
      case SpeciesAuxData::AT_PARTITIONFUNCTION_TFIELD:
        CalculatePartitionFctFromData_perturbed(
            dq_t_dT,
            atm_t,
            dt,
            q_t,
            partition_data[0].get_numeric_grid(0),
            partition_data[0].data,
            1);
        break;
      default:
        throw std::runtime_error("Unknown partition type requested.\n");
        break;
    }

    // Note that this should scale with q_tref, but we do not need that parameter here...
    dpartition_ratio_dT = -dq_t_dT / q_t;
  }

  // Following Futbolin's division into two parts for the Boltzmann ratio because
  // gamma is also used for the NLTE part later on
  const Numeric gamma = exp(-PLANCK_CONST * line_f / (BOLTZMAN_CONST * atm_t));
  const Numeric gamma_ref =
      exp(-PLANCK_CONST * line_f / (BOLTZMAN_CONST * line_t));
  dK2_dT = -line_f * PLANCK_CONST / BOLTZMAN_CONST / atm_t / atm_t *
           (gamma / (1.0 - gamma_ref));

  if (do_nlte) {
    atm_tv_low = line_evlow_index < 0 ? -1.0 : atm_t_nlte[line_evlow_index];
    atm_tv_upp = line_evupp_index < 0 ? -1.0 : atm_t_nlte[line_evupp_index];
    const Numeric gamma_p = 1 / gamma;

    //r_low and r_upp are ratios for the population level compared to LTE conditions
    Numeric dr_low, dr_upp;
    if (atm_tv_low >
        1e-4 *
            atm_t)  // where 1e-4 is considered a small number so that the multiplication in the denominator does not reach zero
    {
      const Numeric r_low = exp(-line_evlow / BOLTZMAN_CONST *
                                (atm_t - atm_tv_low) / (atm_t * atm_tv_low));
      dr_low = r_low * line_evlow;
    } else if (atm_tv_low >= 0.0) {
      dr_low = 0.0;
    } else {
      constexpr Numeric r_low = 1.0;
      dr_low = r_low * line_evlow;
    }

    if (atm_tv_upp >
        1e-4 *
            atm_t)  // where 1e-4 is considered a small number so that the multiplication in the denominator does not reach zero
    {
      const Numeric r_upp = exp(-line_evupp / BOLTZMAN_CONST *
                                (atm_t - atm_tv_upp) / (atm_t * atm_tv_upp));

      dr_upp = -r_upp * (line_evupp - PLANCK_CONST * line_f) * gamma;
    } else if (atm_tv_upp >= 0.0) {
      dr_upp = 0.0;
    } else {
      constexpr Numeric r_upp = 1.0;
      dr_upp = -r_upp * (line_evupp - PLANCK_CONST * line_f) * gamma;
    }

    // Both are unity when in LTE
    dabs_nlte_ratio_dT =
        ((dr_upp + dr_low) / (gamma - 1.0) +
         abs_nlte_ratio * PLANCK_CONST * line_f / (gamma_p - 1.0)) /
        BOLTZMAN_CONST / atm_t / atm_t;
  }
}
void GetLineScalingData_dF0(Numeric& dK2_dF0,
                            Numeric& dabs_nlte_ratio_dF0,
                            const Numeric& atm_t,
                            const Numeric& line_t,
                            const Numeric& atm_tv_low,
                            const Numeric& atm_tv_upp,
                            const Numeric& line_evlow,
                            const Numeric& line_evupp,
                            const Numeric& line_f) {
  // Physical constants
  extern const Numeric PLANCK_CONST;
  extern const Numeric BOLTZMAN_CONST;

  // Following Futbolin's division into two parts for the Boltzmann ratio because
  // gamma is also used for the NLTE part later on
  const Numeric gamma = exp(-PLANCK_CONST * line_f / (BOLTZMAN_CONST * atm_t));
  const Numeric gamma_ref =
      exp(-PLANCK_CONST * line_f / (BOLTZMAN_CONST * line_t));

  // Note lack of division with K2
  dK2_dF0 = (PLANCK_CONST * gamma_ref * (gamma - 1)) /
                (line_t * BOLTZMAN_CONST * (gamma_ref - 1) * (gamma_ref - 1)) -
            (PLANCK_CONST * gamma) / (atm_t * BOLTZMAN_CONST * (gamma_ref - 1));

  if (atm_tv_low > 0 || atm_tv_upp > 0) {
    //r_low and r_upp are ratios for the population level compared to LTE conditions
    Numeric r_low, r_upp;
    if (atm_tv_low >
        0.0)  // where 1e-4 is considered a small number so that the multiplication in the denominator does not reach zero
      r_low = exp(-line_evlow / BOLTZMAN_CONST * (atm_t - atm_tv_low) /
                  (atm_t * atm_tv_low));
    else if (atm_tv_low == 0.0)
      throw std::runtime_error(
          "A line has been defined with zero vibrational temperature.\nThis is not physical.\n");
    else
      r_low = 1.0;

    if (atm_tv_upp >
        0.0)  // where 1e-4 is considered a small number so that the multiplication in the denominator does not reach zero
      r_upp = exp(-line_evupp / BOLTZMAN_CONST * (atm_t - atm_tv_upp) /
                  (atm_t * atm_tv_upp));
    else if (atm_tv_upp == 0.0)
      throw std::runtime_error(
          "A line has been defined with zero vibrational temperature.\nThis is not physical.\n");
    else
      r_upp = 1.0;

    dabs_nlte_ratio_dF0 =
        -PLANCK_CONST * gamma * (r_low - r_upp) /
        (atm_t * BOLTZMAN_CONST * (gamma - 1.0) * (gamma - 1.0));
  }
}

void partition_function(Numeric& q_ref,
                        Numeric& q_t,
                        const Numeric& line_t,
                        const Numeric& atm_t,
                        const SpeciesAuxData::AuxType& partition_type,
                        const ArrayOfGriddedField1& partition_data) {
  switch (partition_type) {
    case SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF:
      CalculatePartitionFctFromCoeff(
          q_ref, q_t, line_t, atm_t, partition_data[0].data);
      break;
    case SpeciesAuxData::AT_PARTITIONFUNCTION_TFIELD:
      CalculatePartitionFctFromData(q_ref,
                                    q_t,
                                    line_t,
                                    atm_t,
                                    partition_data[0].get_numeric_grid(0),
                                    partition_data[0].data,
                                    1);
      break;
    default:
      throw std::runtime_error(
          "Unknown or deprecated partition type requested.\n");
      break;
  }
}

void dpartition_function_dT(Numeric& dq_t_dT,
                            const Numeric& q_t,
                            const Numeric& atm_t,
                            const Numeric& dT,
                            const SpeciesAuxData::AuxType& partition_type,
                            const ArrayOfGriddedField1& partition_data) {
  switch (partition_type) {
    case SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF:
      CalculatePartitionFctFromCoeff_dT(dq_t_dT, atm_t, partition_data[0].data);
      break;
    case SpeciesAuxData::AT_PARTITIONFUNCTION_TFIELD:
      CalculatePartitionFctFromData_perturbed(
          dq_t_dT,
          atm_t,
          dT,
          q_t,
          partition_data[0].get_numeric_grid(0),
          partition_data[0].data,
          1);
      break;
    default:
      throw std::runtime_error("Unknown partition type requested.\n");
      break;
  }
}

void CalculatePartitionFctFromData(Numeric& q_ref,
                                   Numeric& q_t,
                                   const Numeric& ref,
                                   const Numeric& t,
                                   ConstVectorView t_grid,
                                   ConstVectorView q_grid,
                                   const Index& interp_order) {
  GridPosPoly gp_t, gp_ref;
  gridpos_poly(gp_t, t_grid, t, interp_order);
  gridpos_poly(gp_ref, t_grid, ref, interp_order);
  Vector itw_t(gp_t.idx.nelem()), itw_ref(gp_ref.idx.nelem());
  interpweights(itw_t, gp_t);
  interpweights(itw_ref, gp_ref);
  q_t = interp(itw_t, q_grid, gp_t);
  q_ref = interp(itw_ref, q_grid, gp_ref);
}

void CalculatePartitionFctFromData_perturbed(Numeric& dQ_dT,
                                             const Numeric& t,
                                             const Numeric& dT,
                                             const Numeric& q_t,
                                             ConstVectorView t_grid,
                                             ConstVectorView q_grid,
                                             const Index& interp_order) {
  GridPosPoly gp_t;
  gridpos_poly(gp_t, t_grid, t + dT, interp_order);
  Vector itw_t(gp_t.idx.nelem());
  interpweights(itw_t, gp_t);
  const Numeric q_t2 = interp(itw_t, q_grid, gp_t);

  // FIXME:  Is there a way to have interp return the derivative instead?  This way always undershoots the curve...  Should we have other point-derivatives, e.g., 2-point?
  dQ_dT = (q_t - q_t2) / dT;
}

void CalculatePartitionFctFromCoeff(Numeric& q_ref,
                                    Numeric& q_t,
                                    const Numeric& ref,
                                    const Numeric& t,
                                    ConstVectorView q_grid) {
  Numeric result_t = 0.;
  Numeric result_ref = 0.;
  Numeric exponent_t = 1.;
  Numeric exponent_ref = 1.;

  Vector::const_iterator it;

  for (it = q_grid.begin(); it != q_grid.end(); ++it) {
    result_t += *it * exponent_t;
    result_ref += *it * exponent_ref;

    exponent_t *= t;
    exponent_ref *= ref;
  }

  q_t = result_t;
  q_ref = result_ref;
}

void CalculatePartitionFctFromCoeff_dT(Numeric& dQ_dT,
                                       const Numeric& t,
                                       ConstVectorView q_grid) {
  Numeric result_t = 0.;
  Numeric exponent_t = 1.;

  Vector::const_iterator it;

  for (Index ii = 1; ii < q_grid.nelem(); ii++) {
    result_t += q_grid[ii] * exponent_t * (Numeric)ii;

    exponent_t *= t;
  }

  dQ_dT = result_t;
}

void CalculatePartitionFctFromVibrotCoeff_dT(Numeric& dQ_dT,
                                             const Numeric& t_vib,
                                             const Numeric& t_rot,
                                             ConstVectorView qvib_grid,
                                             ConstVectorView qrot_grid) {
  throw std::runtime_error(
      "Vibrot does not yet work with propmat partial derivatives.\n");

  Numeric dQvibT = 0.;
  Numeric dQrotT = 0.;
  Numeric exponent_t_vib = 1.;
  Numeric exponent_t_rot = 1.;

  for (Index ii = 1; ii < qvib_grid.nelem(); ii++) {
    dQvibT += qvib_grid[ii] * exponent_t_vib * (Numeric)ii;
    dQrotT += qrot_grid[ii] * exponent_t_rot * (Numeric)ii;

    exponent_t_rot *= t_rot;
    exponent_t_vib *= t_vib;
  }
  //FIXME:  This is wrong...
  dQ_dT = dQvibT * dQrotT;
}

Numeric SingleCalculatePartitionFctFromCoeff(const Numeric& T,
                                             ConstVectorView q_grid) {
  Numeric result = 0.;
  Numeric TN = 1;

  for (Index i = 0; i < q_grid.nelem(); i++) {
    result += TN * q_grid[i];
    TN *= T;
  }

  return result;
}

Numeric SingleCalculatePartitionFctFromCoeff_dT(const Numeric& T,
                                                ConstVectorView q_grid) {
  Numeric result = 0;
  Numeric TN = 1;

  for (Index i = 1; i < q_grid.nelem(); i++) {
    result += Numeric(i) * TN * q_grid[i];
    TN *= T;
  }

  return result;
}

Numeric SingleCalculatePartitionFctFromData(const Numeric& T,
                                            ConstVectorView t_grid,
                                            ConstVectorView q_grid,
                                            const Index& interp_order) {
  GridPosPoly gp;
  gridpos_poly(gp, t_grid, T, interp_order);
  Vector itw(gp.idx.nelem());
  interpweights(itw, gp);
  return interp(itw, q_grid, gp);
}

Numeric SingleCalculatePartitionFctFromData_dT(const Numeric& QT,
                                               const Numeric& T,
                                               const Numeric& dT,
                                               ConstVectorView t_grid,
                                               ConstVectorView q_grid,
                                               const Index& interp_order) {
  GridPosPoly gp;
  gridpos_poly(gp, t_grid, T + dT, interp_order);
  Vector itw(gp.idx.nelem());
  interpweights(itw, gp);
  return (interp(itw, q_grid, gp) - QT) / dT;
}

Numeric single_partition_function(const Numeric& T,
                                  const SpeciesAuxData::AuxType& partition_type,
                                  const ArrayOfGriddedField1& partition_data) {
  switch (partition_type) {
    case SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF:
      return SingleCalculatePartitionFctFromCoeff(T, partition_data[0].data);
    case SpeciesAuxData::AT_PARTITIONFUNCTION_TFIELD:
      return SingleCalculatePartitionFctFromData(
          T, partition_data[0].get_numeric_grid(0), partition_data[0].data, 1);
    default:
      throw std::runtime_error(
          "Unknown or deprecated partition type requested.\n");
  }
}

Numeric dsingle_partition_function_dT(
    const Numeric& QT,
    const Numeric& T,
    const Numeric& dT,
    const SpeciesAuxData::AuxType& partition_type,
    const ArrayOfGriddedField1& partition_data) {
  switch (partition_type) {
    case SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF:
      return SingleCalculatePartitionFctFromCoeff_dT(T, partition_data[0].data);
    case SpeciesAuxData::AT_PARTITIONFUNCTION_TFIELD:
      return SingleCalculatePartitionFctFromData_dT(
          QT,
          T,
          dT,
          partition_data[0].get_numeric_grid(0),
          partition_data[0].data,
          1);
    default:
      throw std::runtime_error(
          "Unknown or deprecated partition type requested.\n");
  }
}

Numeric stimulated_emission(Numeric T, Numeric F0) {
  using namespace Constant;
  static constexpr Numeric c1 = -h / k;
  return std::exp(c1 * F0 / T);
}

Numeric dstimulated_emissiondT(Numeric T, Numeric F0) {
  using namespace Constant;
  static constexpr Numeric c1 = -h / k;
  return -F0 * c1 * std::exp(F0 * c1 / T) / pow2(T);
}

Numeric dstimulated_emissiondF0(Numeric T, Numeric F0) {
  using namespace Constant;
  static constexpr Numeric c1 = -h / k;
  return c1 * std::exp(F0 * c1 / T) / T;
}

Numeric stimulated_relative_emission(const Numeric& gamma,
                                     const Numeric& gamma_ref) {
  return (1. - gamma) / (1. - gamma_ref);
}

Numeric dstimulated_relative_emission_dT(const Numeric& gamma,
                                         const Numeric& gamma_ref,
                                         const Numeric& F0,
                                         const Numeric& T) {
  extern const Numeric PLANCK_CONST;
  extern const Numeric BOLTZMAN_CONST;
  static const Numeric c = -PLANCK_CONST / BOLTZMAN_CONST;

  return c * F0 * gamma / (T * T * (1. - gamma_ref));
}

Numeric dstimulated_relative_emission_dF0(const Numeric& gamma,
                                          const Numeric& gamma_ref,
                                          const Numeric& T,
                                          const Numeric& T0) {
  extern const Numeric PLANCK_CONST;
  extern const Numeric BOLTZMAN_CONST;
  static const Numeric c = -PLANCK_CONST / BOLTZMAN_CONST;

  const Numeric g0 = 1 - gamma_ref;
  const Numeric g = 1 - gamma;

  return c * (g * gamma_ref / (T0 * g0 * g0) - gamma / (T * g0));
}

// Ratio of boltzman emission at T and T0
Numeric boltzman_ratio(const Numeric& T, const Numeric& T0, const Numeric& E0) {
  extern const Numeric BOLTZMAN_CONST;
  static const Numeric c = 1 / BOLTZMAN_CONST;

  return exp(E0 * c * (T - T0) / (T * T0));
}

Numeric dboltzman_ratio_dT(const Numeric& boltzmann_ratio,
                           const Numeric& T,
                           const Numeric& E0) {
  extern const Numeric BOLTZMAN_CONST;
  static const Numeric c = 1 / BOLTZMAN_CONST;

  return E0 * c * boltzmann_ratio / (T * T);
}

// Boltzmann factor at T
Numeric boltzman_factor(Numeric T, Numeric E0) {
  using namespace Constant;
  static constexpr Numeric c1 = -1 / k;
  return std::exp(c1 * E0 / T);
}

// Boltzmann factor at T
Numeric dboltzman_factordT(Numeric T, Numeric E0) {
  using namespace Constant;
  static constexpr Numeric c1 = -1 / k;
  return -E0 * c1 * std::exp(E0 * c1 / T) / pow2(T);
}

// Boltzmann factor at T
Numeric dboltzman_factordE0(Numeric T, Numeric E0) {
  using namespace Constant;
  static constexpr Numeric c1 = -1 / k;
  return c1 * std::exp(E0 * c1 / T) / T;
}

Numeric absorption_nlte_ratio(const Numeric& gamma,
                              const Numeric& r_upp,
                              const Numeric& r_low) {
  return (r_low - r_upp * gamma) / (1 - gamma);
}

Numeric dabsorption_nlte_rate_dT(const Numeric& gamma,
                                 const Numeric& T,
                                 const Numeric& F0,
                                 const Numeric& El,
                                 const Numeric& Eu,
                                 const Numeric& r_upp,
                                 const Numeric& r_low) {
  extern const Numeric PLANCK_CONST;
  extern const Numeric BOLTZMAN_CONST;
  static const Numeric c = 1 / BOLTZMAN_CONST;

  if (El < 0 or Eu < 0) {
    std::ostringstream os;
    os << "It is considered undefined behavior to NLTE and "
       << "temperature Jacobian without defining all "
       << "vibrational energy states";
    throw std::runtime_error(os.str());
  }

  const Numeric x = 1 / (T * (gamma - 1));
  const Numeric hf = F0 * PLANCK_CONST;

  return x * x * c *
         ((gamma - 1) * (El * r_low - Eu * gamma * r_upp) -
          hf * gamma * (r_low - r_upp));
}

Numeric dabsorption_nlte_rate_dF0(const Numeric& gamma,
                                  const Numeric& T,
                                  const Numeric& r_upp,
                                  const Numeric& r_low) {
  extern const Numeric PLANCK_CONST;
  extern const Numeric BOLTZMAN_CONST;
  static const Numeric c = -PLANCK_CONST / BOLTZMAN_CONST;

  return c * gamma * (r_low - r_upp) / (T * (gamma * gamma - 2 * gamma + 1));
}

Numeric dabsorption_nlte_rate_dTl(const Numeric& gamma,
                                  const Numeric& T,
                                  const Numeric& Tl,
                                  const Numeric& El,
                                  const Numeric& r_low) {
  extern const Numeric BOLTZMAN_CONST;

  const Numeric x = 1 / (BOLTZMAN_CONST * T);
  const Numeric y = 1 / Tl;

  return El * x * y * y * T * r_low / (gamma - 1);
}

Numeric dabsorption_nlte_rate_dTu(const Numeric& gamma,
                                  const Numeric& T,
                                  const Numeric& Tu,
                                  const Numeric& Eu,
                                  const Numeric& r_upp) {
  extern const Numeric BOLTZMAN_CONST;

  const Numeric x = 1 / (BOLTZMAN_CONST * T);
  const Numeric y = 1 / Tu;

  return Eu * x * y * y * T * gamma * r_upp / (gamma - 1);
}
