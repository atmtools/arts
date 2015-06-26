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
                        Numeric& boltzmann_ratio, 
                        Numeric& abs_nlte_ratio, 
                        Numeric& src_nlte_ratio, 
                        const SpeciesAuxData::AuxType& partition_type,
                        const ArrayOfGriddedField1& partition_data,
                        const Numeric& atm_t,
                        const Numeric& line_t,
                        const Numeric& line_f,
                        const Numeric& line_elow,
                        const bool&    do_nlte,
                        const Numeric& line_evlow,
                        const Numeric& line_evupp,
                        const Index& line_evlow_index,
                        const Index& line_evupp_index,
                        ConstVectorView atm_t_nlte)
{
    // Physical constants
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    
    // 1:  partition_ratio is not calculated and therefore follows the old path
    //if( 1 == mpartitionfunctiondata.GetPartitionFunctionDataParams(partition_ratio, line_t, atm_t) )
    //partition_ratio = IsotopologueData().CalculatePartitionFctRatio(line_t, atm_t);
    partition_ratio=1;
    
    if(q_t<0 || q_ref<0)
    {
      switch(partition_type)
      {
        case SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF:
          CalculatePartitionFctFromCoeff(q_ref, q_t, line_t, atm_t,
                                         partition_data[0].data);
          break;
        case SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF_VIBROT:
          throw std::runtime_error("The vib-rotational partition functions specified by your user input is not yet implemented in ARTS.\n");
          break;
        case SpeciesAuxData::AT_PARTITIONFUNCTION_TFIELD:
          CalculatePartitionFctFromData(q_ref, q_t, line_t, atm_t,
                                        partition_data[0].get_numeric_grid(0),
                                        partition_data[0].data, 
                                        1);
          break;
        default:
          throw std::runtime_error("Unknown partition type requested.\n");
          break;
      }
    }
    
    partition_ratio = q_ref/q_t;
    
    // Following Futbolin's division into two parts for the Boltzmann ratio because
    // gamma is also used for the NLTE part later on
    const Numeric gamma = exp( - PLANCK_CONST * line_f / ( BOLTZMAN_CONST * atm_t ) );
    const Numeric gamma_ref = exp( - PLANCK_CONST * line_f / ( BOLTZMAN_CONST * line_t ) );
    
    // Stimulated emission
    const Numeric se = (1.-gamma)/(1.-gamma_ref);
    
    // Boltzmann level
    const Numeric sb = exp( line_elow / BOLTZMAN_CONST * (atm_t-line_t)/(atm_t*line_t) );

    boltzmann_ratio = sb*se;
    
    if(do_nlte)
    {
      // Test the NLTE of the line
      const Numeric& atm_tv_low = line_evlow_index<0?-1.0:atm_t_nlte[line_evlow_index];
      const Numeric& atm_tv_upp = line_evupp_index<0?-1.0:atm_t_nlte[line_evupp_index];
      
      //r_low and r_upp are ratios for the population level compared to LTE conditions
      Numeric r_low, r_upp;
      if( atm_tv_low > 1e-4 *atm_t ) // where 1e-4 is considered a small number so that the multiplication in the denominator does not reach zero
        r_low = exp( - line_evlow / BOLTZMAN_CONST * (atm_t-atm_tv_low) / (atm_t*atm_tv_low) );
      else if( atm_tv_low >= 0.0 )
        r_low = 0.0;
      else
        r_low = 1.0;

      if( atm_tv_upp > 1e-4 *atm_t ) // where 1e-4 is considered a small number so that the multiplication in the denominator does not reach zero
        r_upp = exp( - line_evupp / BOLTZMAN_CONST * (atm_t-atm_tv_upp) / (atm_t*atm_tv_upp) );
      else if( atm_tv_upp >= 0.0 )
        r_upp = 0.0;
      else
        r_upp = 1.0;

      abs_nlte_ratio = (r_low - r_upp * gamma ) / ( 1 - gamma );
      src_nlte_ratio = r_upp;
    }
}

void CalculatePartitionFctFromData( Numeric& q_ref, 
                                    Numeric& q_t, 
                                    const Numeric& ref, 
                                    const Numeric& t,
                                    ConstVectorView t_grid, 
                                    ConstVectorView q_grid, 
                                    const Index& interp_order)
{
  GridPosPoly gp_t, gp_ref;
  gridpos_poly( gp_t,   t_grid, t,   interp_order);
  gridpos_poly( gp_ref, t_grid, ref, interp_order);
  Vector itw_t(gp_t.idx.nelem()), itw_ref(gp_ref.idx.nelem());
  interpweights(itw_t,   gp_t );
  interpweights(itw_ref, gp_ref );
  q_t   = interp(  itw_t,   q_grid, gp_t );
  q_ref = interp(  itw_ref, q_grid, gp_ref );
}

void CalculatePartitionFctFromCoeff(Numeric& q_ref, 
                                    Numeric& q_t, 
                                    const Numeric& ref, 
                                    const Numeric& t,
                                    ConstVectorView q_grid)
{
  Numeric result_t   = 0.;
  Numeric result_ref = 0.;
  Numeric exponent_t   = 1.;
  Numeric exponent_ref = 1.;

  Vector::const_iterator it;

  for (it=q_grid.begin(); it != q_grid.end(); ++it)
  {
    result_t   += *it * exponent_t;
    result_ref += *it * exponent_ref;
    
    exponent_t   *= t;
    exponent_ref *= ref;
  }
  
  q_t   = result_t;
  q_ref = result_ref;
}