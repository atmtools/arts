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
    
    // This is for a future update possibility
    const bool do_rotational = false;
    
    if(q_t<0 || q_ref<0)
    {
      partition_function( q_ref, q_t,
                          line_t, atm_t,
                          partition_type, partition_data, do_rotational);
      
      partition_ratio = q_ref/q_t;
    }
    
    // Following Futbolin's division into two parts for the Boltzmann ratio because
    // gamma is also used for the NLTE part later on
    const Numeric gamma = exp( - PLANCK_CONST * line_f / ( BOLTZMAN_CONST * atm_t ) );
    const Numeric gamma_ref = exp( - PLANCK_CONST * line_f / ( BOLTZMAN_CONST * line_t ) );
    
    // Stimulated emission
    K2 = (1.-gamma)/(1.-gamma_ref);
    
    // Boltzmann level
    K1 = exp( line_elow / BOLTZMAN_CONST * (atm_t-line_t)/(atm_t*line_t) );
    
    if(do_nlte)
    {
        // Test the NLTE of the line and find if we should use atmospheric temperatures or not
        const Numeric& atm_tv_low = line_evlow_index<0?-1.0:atm_t_nlte[line_evlow_index];
        const Numeric& atm_tv_upp = line_evupp_index<0?-1.0:atm_t_nlte[line_evupp_index];
      
      //r_low and r_upp are ratios for the population level compared to LTE conditions
      Numeric r_low, r_upp;
      if( atm_tv_low > 0.0 ) // where 1e-4 is considered a small number so that the multiplication in the denominator does not reach zero
        r_low = exp( - line_evlow / BOLTZMAN_CONST * (atm_t-atm_tv_low) / (atm_t*atm_tv_low) );
      else if(atm_tv_low == 0.0)
        throw std::runtime_error("A line has been defined with zero vibrational temperature.\nThis is not physical.\n");
      else
        r_low = 1.0;

      if( atm_tv_upp > 0.0 ) // where 1e-4 is considered a small number so that the multiplication in the denominator does not reach zero
        r_upp = exp( - line_evupp / BOLTZMAN_CONST * (atm_t-atm_tv_upp) / (atm_t*atm_tv_upp) );
      else if(atm_tv_upp  == 0.0)
        throw std::runtime_error("A line has been defined with zero vibrational temperature.\nThis is not physical.\n");
      else
        r_upp = 1.0;

      // Both are unity when in LTE
      abs_nlte_ratio = (r_low - r_upp * gamma ) / ( 1 - gamma );
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
    
    const bool do_rotational = false;
    
    // Test the NLTE of the line and find if we should use atmospheric temperatures or not
    atm_tv_low = line_evlow_index<0?-1.0:atm_t_nlte[line_evlow_index];
    atm_tv_upp = line_evupp_index<0?-1.0:atm_t_nlte[line_evupp_index];
    
    if(dq_t_dT<0)
    {
        switch(partition_type)
        {
            case SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF:
                CalculatePartitionFctFromCoeff_dT(dq_t_dT, atm_t,
                                               partition_data[0].data);
                break;
            case SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF_VIBROT:
                if(!do_rotational)
                    CalculatePartitionFctFromVibrotCoeff_dT(dq_t_dT, atm_t, atm_t,
                                                         partition_data[0].data,partition_data[1].data);
                else
                    CalculatePartitionFctFromVibrotCoeff_dT(dq_t_dT, atm_t, /*t_rot*/ atm_t, //This must be implemented!  How...???
                                                            partition_data[0].data,partition_data[1].data);
                break;
            case SpeciesAuxData::AT_PARTITIONFUNCTION_TFIELD:
                CalculatePartitionFctFromData_perturbed(dq_t_dT, atm_t, dt, q_t,
                                              partition_data[0].get_numeric_grid(0),
                                              partition_data[0].data, 
                                              1);
                break;
            default:
                throw std::runtime_error("Unknown partition type requested.\n");
                break;
        }
        
        // Note that this should scale with q_tref, but we do not need that parameter here...
        dpartition_ratio_dT = -dq_t_dT/q_t; 
    }
    
    // Following Futbolin's division into two parts for the Boltzmann ratio because
    // gamma is also used for the NLTE part later on
    const Numeric gamma = exp( - PLANCK_CONST * line_f / ( BOLTZMAN_CONST * atm_t ) );
    const Numeric gamma_ref = exp( - PLANCK_CONST * line_f / ( BOLTZMAN_CONST * line_t ) );
    dK2_dT = -line_f*PLANCK_CONST/BOLTZMAN_CONST/atm_t/atm_t * (gamma/(1.0-gamma_ref));
    
    if(do_nlte)
    {   
        const Numeric gamma_p = 1/gamma;
        
        //r_low and r_upp are ratios for the population level compared to LTE conditions
        Numeric r_low, r_upp;
        Numeric dr_low, dr_upp;
        if( atm_tv_low > 1e-4 *atm_t ) // where 1e-4 is considered a small number so that the multiplication in the denominator does not reach zero
        {
            r_low = exp( - line_evlow / BOLTZMAN_CONST * (atm_t-atm_tv_low) / (atm_t*atm_tv_low) );
            dr_low = r_low * line_evlow ;
        }
        else if( atm_tv_low >= 0.0 )
        {
            r_low = 0.0;
            dr_low = 0.0;
        }
        else
        {
            r_low = 1.0;
            dr_low = r_low * line_evlow  ;
        }
        
        if( atm_tv_upp > 1e-4 *atm_t ) // where 1e-4 is considered a small number so that the multiplication in the denominator does not reach zero
        {
            r_upp = exp( - line_evupp / BOLTZMAN_CONST * (atm_t-atm_tv_upp) / (atm_t*atm_tv_upp) );
            dr_upp = - r_upp *(line_evupp - PLANCK_CONST*line_f)*gamma;
        }
        else if( atm_tv_upp >= 0.0 )
        {
            r_upp = 0.0;
            dr_upp = 0.0;
        }
        else
        {
            r_upp = 1.0;
            dr_upp = - r_upp *(line_evupp - PLANCK_CONST*line_f)*gamma;
        }
        
        // Both are unity when in LTE
        dabs_nlte_ratio_dT = ((dr_upp+dr_low) / ( gamma - 1.0 ) + abs_nlte_ratio* PLANCK_CONST*line_f/(gamma_p - 1.0))/BOLTZMAN_CONST / atm_t / atm_t;
    }
}
void GetLineScalingData_dF0(
                        Numeric& dK2_dF0, 
                        Numeric& dabs_nlte_ratio_dF0,
                        const Numeric& atm_t,
                        const Numeric& line_t,
                        const Numeric& atm_tv_low,
                        const Numeric& atm_tv_upp,
                        const Numeric& line_evlow,
                        const Numeric& line_evupp,
                        const Numeric& line_f)
{
    // Physical constants
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
        
    // Following Futbolin's division into two parts for the Boltzmann ratio because
    // gamma is also used for the NLTE part later on
    const Numeric gamma = exp( - PLANCK_CONST * line_f / ( BOLTZMAN_CONST * atm_t ) );
    const Numeric gamma_ref = exp( - PLANCK_CONST * line_f / ( BOLTZMAN_CONST * line_t ) );
    
    // Note lack of division with K2
    dK2_dF0 = (PLANCK_CONST*gamma_ref*(gamma - 1))/(line_t*BOLTZMAN_CONST*(gamma_ref - 1)*(gamma_ref - 1)) - (PLANCK_CONST*gamma)/(atm_t*BOLTZMAN_CONST*(gamma_ref - 1));
    
    
    if(atm_tv_low>0||atm_tv_upp>0)
    {
        
        //r_low and r_upp are ratios for the population level compared to LTE conditions
        Numeric r_low, r_upp;
        if( atm_tv_low > 0.0 ) // where 1e-4 is considered a small number so that the multiplication in the denominator does not reach zero
            r_low = exp( - line_evlow / BOLTZMAN_CONST * (atm_t-atm_tv_low) / (atm_t*atm_tv_low) );
        else if(atm_tv_low == 0.0)
            throw std::runtime_error("A line has been defined with zero vibrational temperature.\nThis is not physical.\n");
        else
            r_low = 1.0;
        
        if( atm_tv_upp > 0.0 ) // where 1e-4 is considered a small number so that the multiplication in the denominator does not reach zero
            r_upp = exp( - line_evupp / BOLTZMAN_CONST * (atm_t-atm_tv_upp) / (atm_t*atm_tv_upp) );
        else if(atm_tv_upp  == 0.0)
            throw std::runtime_error("A line has been defined with zero vibrational temperature.\nThis is not physical.\n");
        else
            r_upp = 1.0;
        
        dabs_nlte_ratio_dF0 = -PLANCK_CONST*gamma*(r_low-r_upp)/(atm_t*BOLTZMAN_CONST*(gamma-1.0)*(gamma-1.0));
    }
}


void partition_function( Numeric& q_ref,
                         Numeric& q_t,
                         const Numeric& line_t,
                         const Numeric& atm_t,
                         const SpeciesAuxData::AuxType& partition_type,
                         const ArrayOfGriddedField1& partition_data,
                         const bool& do_rotational)
{
  switch(partition_type)
  {
    case SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF:
      CalculatePartitionFctFromCoeff(q_ref, q_t, line_t, atm_t,
                                    partition_data[0].data);
      break;
    case SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF_VIBROT:
      if(!do_rotational)
        CalculatePartitionFctFromVibrotCoeff(q_ref, q_t, line_t, atm_t, atm_t,
                                            partition_data[0].data,partition_data[1].data);
        else
          CalculatePartitionFctFromVibrotCoeff(q_ref, q_t, line_t, atm_t, /*t_rot*/ atm_t, //This must be implemented!
                                              partition_data[0].data,partition_data[1].data);
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

void CalculatePartitionFctFromData_perturbed( Numeric& dQ_dT, 
                                              const Numeric& t,
                                              const Numeric& dT,
                                              const Numeric& q_t,
                                              ConstVectorView t_grid, 
                                              ConstVectorView q_grid, 
                                              const Index& interp_order)
{
    GridPosPoly gp_t;
    gridpos_poly( gp_t,   t_grid, t+dT,   interp_order);
    Vector itw_t(gp_t.idx.nelem());
    interpweights(itw_t,   gp_t );
    const Numeric q_t2   = interp(  itw_t,   q_grid, gp_t );
    
    // FIXME:  Is there a way to have interp return the derivative instead?  This way always undershoots the curve...  Should we have other point-derivatives, e.g., 2-point?
    dQ_dT = (q_t - q_t2)/dT;
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

void CalculatePartitionFctFromVibrotCoeff(Numeric& q_ref, 
                                          Numeric& q_t, 
                                          const Numeric& ref, 
                                          const Numeric& t_vib,
                                          const Numeric& t_rot,
                                          ConstVectorView qvib_grid,
                                          ConstVectorView qrot_grid)
{
    Numeric QvibT   = 0.;
    Numeric QvibT_ref = 0.;
    Numeric QrotT   = 0.;
    Numeric QrotT_ref = 0.;
    Numeric exponent_t_vib   = 1.;
    Numeric exponent_ref = 1.;
    Numeric exponent_t_rot   = 1.;
    
    for (Index ii = 0; ii<qvib_grid.nelem();ii++)
    {
        QvibT   +=   qvib_grid[ii] * exponent_t_vib;
        QrotT   +=   qrot_grid[ii] * exponent_t_rot;
        QvibT_ref += qvib_grid[ii] * exponent_ref;
        QrotT_ref += qrot_grid[ii] * exponent_ref;
        
        exponent_t_rot   *= t_rot;
        exponent_t_vib   *= t_vib;
        exponent_ref *= ref;
    }
    
    q_t   = QvibT*QrotT;
    q_ref = QvibT_ref*QrotT_ref;
}

void CalculatePartitionFctFromCoeff_dT(Numeric& dQ_dT, 
                                       const Numeric& t,
                                       ConstVectorView q_grid)
{
    Numeric result_t   = 0.;
    Numeric exponent_t = 1.;
    
    Vector::const_iterator it;
    
    for(Index ii = 1; ii<q_grid.nelem();ii++)
    {
        result_t   += q_grid[ii] * exponent_t * (Numeric)ii;
        
        exponent_t *= t;
    }
    
    dQ_dT   = result_t;
}

void CalculatePartitionFctFromVibrotCoeff_dT(Numeric& dQ_dT, 
                                             const Numeric& t_vib,
                                             const Numeric& t_rot,
                                             ConstVectorView qvib_grid,
                                             ConstVectorView qrot_grid)
{
    Numeric dQvibT   = 0.;
    Numeric dQrotT   = 0.;
    Numeric exponent_t_vib   = 1.;
    Numeric exponent_t_rot   = 1.;
    
    for(Index ii = 1; ii<qvib_grid.nelem();ii++)
    {
        dQvibT   +=   qvib_grid[ii] * exponent_t_vib*(Numeric)ii;
        dQrotT   +=   qrot_grid[ii] * exponent_t_rot*(Numeric)ii;
        
        exponent_t_rot   *= t_rot;
        exponent_t_vib   *= t_vib;
    }
    //FIXME:  This is wrong...
    dQ_dT   = dQvibT*dQrotT;
    throw std::runtime_error("Vibrot does not yet work with propmat partial derivatives.\n");
}