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

#include "partial_derivatives.h"
#include "absorption.h"
#include "linescaling.h"
#include "global_data.h"
#include "quantum.h"
#include "arts.h"

extern const String PROPMAT_SUBSUBTAG;


/* Helper function.
 * 
 * Function to calculate the partial deivative in a standardized manner.
 * We allow the partial derivative to propagate in three ways:
 * 
 *      1) The derivation of the line shape and the chain rule.
 *         This is the most theoretical method and is from direct derivations.
 *         A perfect example of this is the frequency derivations for wind calculations,
 *         where it is straightforward to calculate the derivative from the lineshape.
 * 
 * Input parameters:
 * 
 *      2) As a ratio on the output cross-section of the line.
 *         Some derivations are not necessary to perform because they only result in
 *         very simple ratios.  For instance, the derivation of Doppler broadening w.r.t.
 *         temperature is simple one over the temperature.  It is therefore convenient to
 *         use this property rather than to perform the long derivation of the line shape.
 * 
 *      3) As a ratio on the phase or attenuation of unmixed lines.
 *         This is necessary for line mixing.
 */
inline void calc_derivative(Numeric& dxsec_dtarget,
                            Numeric& dphase_dtarget,
                            Numeric& dsrc_dtarget,
                            const Numeric& dFA_dx,
                            const Numeric& dFB_dx,
                            const Numeric& dFA_dy,
                            const Numeric& dFB_dy,
                            const Numeric& dx_dtarget,
                            const Numeric& dy_dtarget,
                            const Numeric& FA, 
                            const Numeric& FB,
                            const Numeric& FA_ratio_to_dtarget,
                            const Numeric& FB_ratio_to_dtarget,
                            const Numeric& lma,
                            const Numeric& lmb,
                            const Numeric& lma_real_ratio_to_dtarget,
                            const Numeric& lma_imag_ratio_to_dtarget,
                            const Numeric& lmb_real_ratio_to_dtarget,
                            const Numeric& lmb_imag_ratio_to_dtarget,
                            const Numeric& nlte,
                            const Numeric& dnlte_dtarget,
                            const bool do_phase,
                            const bool do_src) 
{
    if(do_src)
    {
        const Numeric dtarget = dFA_dx * dx_dtarget +  // Lineshape attenuation derivative with frequency term
                                dFA_dy * dy_dtarget +   // Lineshape attenuation derivative with pressure term
                                FA * FA_ratio_to_dtarget +  // When derivation is just a ratio away
                                lma * lma_real_ratio_to_dtarget +  // line mixing attenuation contribution
                                lmb * lmb_real_ratio_to_dtarget;  // line mixing phase contribution
        dxsec_dtarget += dtarget;
        dsrc_dtarget += dtarget * nlte + FA * dnlte_dtarget;
    }
    else 
        dxsec_dtarget +=  dFA_dx * dx_dtarget +  // Lineshape attenuation derivative with frequency term
                          dFA_dy * dy_dtarget +   // Lineshape attenuation derivative with pressure term
                          FA * FA_ratio_to_dtarget +  // When derivation is just a ratio away
                          lma * lma_real_ratio_to_dtarget +  // line mixing attenuation contribution
                          lmb * lmb_real_ratio_to_dtarget;  // line mixing phase contribution
    
    if(do_phase)
        dphase_dtarget += dFB_dx * dx_dtarget +  // Lineshape phase derivative with frequency term
                          dFB_dy * dy_dtarget +   // Lineshape phase derivative with pressure term
                          FB * FB_ratio_to_dtarget +  // When derivation is just a ratio away
                          lma * lma_imag_ratio_to_dtarget +  // line mixing attenuation contribution
                          lmb * lmb_imag_ratio_to_dtarget;  // line mixing phase contribution
    
}


//This function will require much more inputs in the future
void partial_derivatives_lineshape_dependency(ArrayOfMatrix&  partials_attenuation,
                                              ArrayOfMatrix&  partials_phase, 
                                              ArrayOfMatrix&  partials_src, 
                                              const ArrayOfRetrievalQuantity&  flag_partials, 
                                              const ArrayOfIndex&  flag_partials_position,
                                              ConstVectorView CF_A,//no linemixing except DV!
                                              ConstVectorView CF_B,//no linemixing except DV!
                                              ConstVectorView C,
                                              ConstVectorView dFa_dx,
                                              ConstVectorView dFb_dx, 
                                              ConstVectorView dFa_dy, 
                                              ConstVectorView dFb_dy, 
                                              ConstVectorView f_grid,
                                              const Range&    this_f_grid,
                                              const Numeric&  temperature,
                                              const Numeric&  sigma,
                                              const Numeric&  K2,
                                              const Numeric&  dK2_dT,
                                              const Numeric&  K3,
                                              const Numeric&  dK3_dT,
                                              const Numeric&  K4,
                                              // Line parameters
                                              const Numeric&  line_frequency,
                                              const Numeric&  line_strength,
                                              const Numeric&  line_temperature,
                                              const Numeric&  line_E_low,
                                              const Numeric&  line_E_v_low,
                                              const Numeric&  line_E_v_upp,
                                              const Numeric&  line_T_v_low,
                                              const Numeric&  line_T_v_upp,
                                              const Numeric&  Y_LM,
                                              const Numeric&  dY_LM_dT,
                                              const Numeric&  dY_LM0,
                                              const Numeric&  dY_LM1,
                                              const Numeric&  dY_LMexp,
                                              const Numeric&  G_LM,
                                              const Numeric&  dG_LM_dT,
                                              const Numeric&  dG_LM0,
                                              const Numeric&  dG_LM1,
                                              const Numeric&  dG_LMexp,
                                              const Numeric&  DF_LM,
                                              const Numeric&  dDF_LM_dT,
                                              const Numeric&  dDF_LM0,
                                              const Numeric&  dDF_LM1,
                                              const Numeric&  dDF_LMexp,
                                              const QuantumIdentifier&  qi,
                                              // LINE SHAPE
                                              const Index& ind_ls,
                                              const Index& ind_lsn,
                                              const Numeric& df_0,
                                              const Numeric& ddf_dT,
                                              const Numeric& gamma,
                                              const Numeric& dgamma_dT,
                                              const Numeric& dgamma_dSelf,
                                              const Numeric& dgamma_dForeign,
                                              const Numeric& dgamma_dWater,
                                              const Numeric& dpsf_dSelf,
                                              const Numeric& dpsf_dForeign,
                                              const Numeric& dpsf_dWater,
                                              const Numeric& dgamma_dSelfExponent,
                                              const Numeric& dgamma_dForeignExponent,
                                              const Numeric& dgamma_dWaterExponent,
                                              const Numeric& dpsf_dSelfExponent,
                                              const Numeric& dpsf_dForeignExponent,
                                              const Numeric& dpsf_dWaterExponent,
                                              // Partition data parameters
                                              const Numeric&  dQ_dT,
                                              // Magnetic variables
                                              const Numeric&  DF_Zeeman,
                                              const Numeric&  H_mag_Zeeman,
                                              const bool      do_zeeman,
                                              // Programming variables
                                              const Index&    pressure_level_index,
                                              const bool      do_partials_phase,
                                              const bool      do_src)
{
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    
    if(this_f_grid.get_extent()==0)
        return;
    
    ConstVectorView this_f = f_grid[this_f_grid];
    const Index nv = this_f.nelem();
    const Numeric nlte = do_src? K4/K3-1.0 : 0.0;
    const Numeric f0 = line_frequency + df_0 + DF_LM + DF_Zeeman*H_mag_Zeeman;
    
    Vector LM_Fa(nv), LM_Fb(nv), empty_vector(nv);
    for(Index iv=0;iv<nv;iv++)
    {
        LM_Fa[iv] = ( (1.0 + G_LM)*CF_A[iv] + Y_LM*CF_B[iv]);
        if(do_partials_phase)
            LM_Fb[iv] = ( (1.0 + G_LM)*CF_B[iv] - Y_LM*CF_A[iv]);
    }
    
    // Loop over all jacobian_quantities, if a matching quantity is found, then apply the necessary steps to make the jacobian matrix
    for(Index ii=0; ii<flag_partials_position.nelem(); ii++)
    {
      if( (flag_partials[flag_partials_position[ii]] == JacPropMatType::MagneticMagnitude) )
        {
            if(!do_zeeman)
                continue;
            
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = partials_phase[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_src         = do_src?partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            const Numeric& dF_dH = DF_Zeeman;
            
            Numeric dx_dH;
            Vector  dfn_dH_div_dF_dH(nv);
            
            // Calculate the line shape derivative:
            global_data::lineshape_data[ind_ls].dInput_dH()(dx_dH,sigma,dF_dH);
            global_data::lineshape_norm_data[ind_lsn].dFunction_dF0()(dfn_dH_div_dF_dH, f0, this_f, temperature);
            for(Index iv=0;iv<nv;iv++)
            {
                const Numeric ls_A= ( (1.0 + G_LM)*CF_A[iv] + Y_LM*CF_B[iv]), 
                              ls_B= ( (1.0 + G_LM)*CF_B[iv] - Y_LM*CF_A[iv]);
                              
                this_partial_attenuation[iv] += dFa_dx[iv]*dx_dH + ls_A*dfn_dH_div_dF_dH[iv]*dF_dH/C[iv];
                this_partial_phase[iv]       += dFb_dx[iv]*dx_dH + ls_B*dfn_dH_div_dF_dH[iv]*dF_dH/C[iv];
                if(do_src)
                    this_partial_src[iv]     +=(dFa_dx[iv]*dx_dH + ls_A*dfn_dH_div_dF_dH[iv]*dF_dH/C[iv])*nlte;
            }
        }
        else if(flag_partials[flag_partials_position[ii]] == JacPropMatType::MagneticTheta)
        {/* Pass on this.  Can be done easier in later parts of the code execution.  
            Will not require line shape partials. */}
        else if(flag_partials[flag_partials_position[ii]] == JacPropMatType::MagneticEta)
        {/* Pass on this.  Can be done easier in later parts of the code execution.  
            Will not require line shape partials. */}
        else if(flag_partials[flag_partials_position[ii]] == JacPropMatType::MagneticU or
                flag_partials[flag_partials_position[ii]] == JacPropMatType::MagneticV or 
                flag_partials[flag_partials_position[ii]] == JacPropMatType::MagneticW)
        {/* Pass on these.  These are done by perturbation in later parts of the code execution since I found it too complicated for now.  FIXME: Richard */}
        else if(is_frequency_parameter(flag_partials[flag_partials_position[ii]]))
        {
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = do_partials_phase?
            partials_phase[ii](this_f_grid, pressure_level_index):empty_vector;
            VectorView this_partial_src         = do_src?
            partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            Numeric dF_dF;
            Vector dfn_dF(nv);
            
            // Calculate the line shape derivative:
            global_data::lineshape_data[ind_ls].dInput_dF()(dF_dF,sigma);
            global_data::lineshape_norm_data[ind_lsn].dFunction_dF()(dfn_dF, f0, this_f, temperature);
            
            for(Index iv=0;iv<nv;iv++)
            {
                const Numeric ls_A= ( (1.0 + G_LM)*CF_A[iv] + Y_LM*CF_B[iv]), 
                              ls_B= ( (1.0 + G_LM)*CF_B[iv] - Y_LM*CF_A[iv]);
                
                this_partial_attenuation[iv]+=dfn_dF[iv]/C[iv]*ls_A+dFa_dx[iv]*dF_dF;
                if(do_partials_phase)
                    this_partial_phase[iv]      +=dfn_dF[iv]/C[iv]*ls_B+dFb_dx[iv]*dF_dF;
                if(do_src)
                    this_partial_src[iv]+=(dfn_dF[iv]/C[iv]*ls_A+dFa_dx[iv]*dF_dF)*nlte;
                
                //NOTE:  Still missing wind term in this derivative.  Must be multiplied by cf/(c+W)^2 at some point to get this correct... As it stands though, this is the frequency derivative
                //NOTE:  Another missing aspect is that the dW/d{u,v,w} term must also be added for those components
            }
        }
        else if(flag_partials[flag_partials_position[ii]] == JacPropMatType::Temperature)
        {
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = do_partials_phase?partials_phase[ii](this_f_grid, pressure_level_index):empty_vector;
            VectorView this_partial_src         = do_src?partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            const Numeric kT2 = BOLTZMAN_CONST*temperature*temperature;
            
            // Line strength partials... NOTE: All but dK4 are divided by original values to fit ls_A/B
            const Numeric dK1 = line_E_low/kT2; 
            const Numeric dK2 = dK2_dT/K2; 
            const Numeric dK3 = do_src?(dK3_dT/K3):0.0;
            const Numeric dK4 = line_E_v_upp>=0?-K4*line_E_v_upp/kT2:0.0;
            const Numeric dS_dT = dK1+dK2+dQ_dT+dK3; // Note that missing dn/dT is handled later
            
            // Derivative of sigma with regards to temperature
            const Numeric dsigma_dT = 0.5*sigma / temperature;
            
            // Setting up for the partials of the inner loop
            Numeric dP_dT, dFu_dT;
            Vector  dfn_dT(nv), 
                    dF_dT(nv);
            
            // Calculate the line shape derivative:
            global_data::lineshape_data[ind_ls].dInput_dT()(dF_dT, dP_dT, dFu_dT,
                                                            this_f, f0, sigma, ddf_dT, dDF_LM_dT,
                                                            dsigma_dT, gamma, dgamma_dT);
            global_data::lineshape_norm_data[ind_lsn].dFunction_dT()(dfn_dT, f0,
                                                                     this_f, temperature);
            
            for(Index iv=0;iv<nv;iv++)
            {
                const Numeric ls_A= ( (1.0 + G_LM)*CF_A[iv] + Y_LM*CF_B[iv]), 
                              ls_B= ( (1.0 + G_LM)*CF_B[iv] - Y_LM*CF_A[iv]);
                
                this_partial_attenuation[iv] += 
                (dS_dT + dfn_dT[iv]/C[iv] + dFu_dT) * ls_A + //Line strength and factors
                dG_LM_dT  * CF_A[iv]  + dY_LM_dT * CF_B[iv]  + //Line Mixing (absolute)
                dF_dT[iv] * dFa_dx[iv] +                        //Frequency line shape
                dP_dT     * dFa_dy[iv];                         //Pressure line shape
                
                if(do_partials_phase)// Minus signs should be here due to iFb, though this must be tested!
                    this_partial_phase[iv]   += 
                    (dS_dT + dfn_dT[iv]/C[iv] + dFu_dT) * ls_B + //Line strength
                    dG_LM_dT  * CF_B[iv]  - dY_LM_dT * CF_A[iv]  + //Line Mixing (absolute)
                    dF_dT[iv] * dFb_dx[iv] +                       //Frequency line shape
                    dP_dT     * dFb_dy[iv];                        //Pressure line shape
                
                if(do_src)
                    this_partial_src[iv] += nlte * /*partial attenuation*/
                    ((dS_dT + dfn_dT[iv]/C[iv] + dFu_dT) * ls_A  + //Line strength
                    dG_LM_dT  * CF_A[iv]   + dY_LM_dT * CF_B[iv]   + //Line Mixing (absolute)
                    dF_dT[iv] * dFa_dx[iv]  +                         //Frequency line shape 
                    dP_dT     * dFa_dy[iv]) +                         //Pressure line shape
                    ls_A/K3*(dK4-K4*dK3);                            //Source term ratio
                    
            }
            
            // Ready and done!  So complicated that I need plenty of testing!
            
        } 
        else if(flag_partials[flag_partials_position[ii]] == JacPropMatType::VMR)
        {
            // Line shape cross section does not depend on VMR.  NOTE:  Ignoring self-pressure broadening.
        }
        else if(flag_partials[flag_partials_position[ii]] == JacPropMatType::LineCenter)
        {
          if(!line_match_line(flag_partials[flag_partials_position[ii]].QuantumIdentity(),
            qi.Species(), qi.Isotopologue(), qi.LowerQuantumNumbers(),qi.UpperQuantumNumbers()))
                continue;
            
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = do_partials_phase?partials_phase[ii](this_f_grid, pressure_level_index):empty_vector;
            VectorView this_partial_src         = do_src?partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            Numeric dK2_dF0, dK3_dF0=0.0;
            GetLineScalingData_dF0(dK2_dF0, dK3_dF0,temperature,line_temperature, line_T_v_low, line_T_v_upp, line_E_v_low, line_E_v_upp, line_frequency);
            dK2_dF0 /= K2; // to just multiply with ls_{A,B}
            if(do_src)
                dK3_dF0/=K3;
            const Numeric dS_dF0 = dK2_dF0+dK3_dF0;
            
            Numeric dF_dF0;
            global_data::lineshape_data[ind_ls].dInput_dF0()(dF_dF0, sigma);
            
            Vector dfn_dF0(nv);
            global_data::lineshape_norm_data[ind_lsn].dFunction_dF0()(dfn_dF0, f0,
                                                                      this_f, temperature);
            
            const Numeric dx_dF0_part = dF_dF0 / f0, dy_dF0 = gamma * dF_dF0 / f0;
            
            for(Index iv=0;iv<nv;iv++)
            {
                const Numeric ratio = (dS_dF0 + dfn_dF0[iv]/C[iv] - 1.0 / f0);
                calc_derivative(this_partial_attenuation[iv],
                                this_partial_phase[iv],
                                this_partial_src[iv],
                                dFa_dx[iv], dFb_dx[iv], dFa_dy[iv], dFb_dy[iv],
                                dx_dF0_part * this_f[iv], dy_dF0, // dx_dtarget, dy_target
                                LM_Fa[iv], LM_Fb[iv], // Full calculations
                                ratio, ratio, // ratios
                                CF_A[iv], CF_B[iv], 0.0, 0.0, 0.0, 0.0, //  When attenuation and phase matters
                                nlte, 0.0,
                                do_partials_phase, do_src);
            }
            
            // That's it!  Note that the output should be strongly correlated to Temperature and Pressure.
            // Also note that to do the fit for the catalog gamma is somewhat different than this
        }
        else if(flag_partials[flag_partials_position[ii]] == JacPropMatType::LineStrength)
        {
          if(!line_match_line(flag_partials[flag_partials_position[ii]].QuantumIdentity(),
            qi.Species(), qi.Isotopologue(), qi.LowerQuantumNumbers(),qi.UpperQuantumNumbers()))
                continue;
            
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = do_partials_phase?partials_phase[ii](this_f_grid, pressure_level_index):empty_vector;
            VectorView this_partial_src         = do_src?partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            const Numeric ratio = 1.0/line_strength;
            
            for(Index iv=0;iv<nv;iv++)
            {
                calc_derivative(this_partial_attenuation[iv],
                                this_partial_phase[iv],
                                this_partial_src[iv],
                                dFa_dx[iv], dFb_dx[iv], dFa_dy[iv], dFb_dy[iv],
                                0.0, 0.0, // dx_dtarget, dy_target
                                LM_Fa[iv], LM_Fb[iv], // Full calculations
                                ratio, ratio, // ratios
                                CF_A[iv], CF_B[iv], 0.0, 0.0, 0.0, 0.0, //  When attenuation and phase matters
                                nlte, 0.0,
                                do_partials_phase, do_src);
            }
            
            // That's it!  Now to wonder if this will grow extremely large, creating an unrealistic jacobian...
        }
        else if(flag_partials[flag_partials_position[ii]] == JacPropMatType::LineGammaSelf    or
                flag_partials[flag_partials_position[ii]] == JacPropMatType::LineGammaForeign or
                flag_partials[flag_partials_position[ii]] == JacPropMatType::LineGammaWater)
        {
          if(!line_match_line(flag_partials[flag_partials_position[ii]].QuantumIdentity(),
            qi.Species(), qi.Isotopologue(), qi.LowerQuantumNumbers(),qi.UpperQuantumNumbers()))
                continue;
            
            Numeric y_constant = 1.0;
            switch(flag_partials[flag_partials_position[ii]].PropMatType())
            {
              case JacPropMatType::LineGammaSelf: y_constant = dgamma_dSelf; break;
              case JacPropMatType::LineGammaForeign: y_constant = dgamma_dForeign; break;
              case JacPropMatType::LineGammaWater: y_constant = dgamma_dWater; break;
              default: throw std::runtime_error("This cannot happen.  A developer has made a mistake.\n");
            }
            
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = do_partials_phase?partials_phase[ii](this_f_grid, pressure_level_index):empty_vector;
            VectorView this_partial_src         = do_src?partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            Numeric dP_dgamma;
            // Calculate the line shape derivative:
            global_data::lineshape_data[ind_ls].dInput_dgamma()(dP_dgamma,sigma);
            
            dP_dgamma *= y_constant;
            
            for(Index iv=0;iv<nv;iv++)
            {
                calc_derivative(this_partial_attenuation[iv],
                                this_partial_phase[iv],
                                this_partial_src[iv],
                                dFa_dx[iv], dFb_dx[iv], dFa_dy[iv], dFb_dy[iv],
                                0.0, dP_dgamma, // dx_dtarget, dy_target
                                LM_Fa[iv], LM_Fb[iv], // Full calculations
                                0.0, 0.0, // ratios
                                CF_A[iv], CF_B[iv], 0.0, 0.0, 0.0, 0.0, //  When attenuation and phase matters
                                nlte, 0.0,
                                do_partials_phase, do_src);
            }   
            // That's it!
        }
        else if(flag_partials[flag_partials_position[ii]] == JacPropMatType::LineGammaSelfExp    or
                flag_partials[flag_partials_position[ii]] == JacPropMatType::LineGammaForeignExp or
                flag_partials[flag_partials_position[ii]] == JacPropMatType::LineGammaWaterExp   or
                flag_partials[flag_partials_position[ii]] == JacPropMatType::LineShiftSelf       or
                flag_partials[flag_partials_position[ii]] == JacPropMatType::LineShiftForeign    or
                flag_partials[flag_partials_position[ii]] == JacPropMatType::LineShiftWater)
        {
            if(!line_match_line(flag_partials[flag_partials_position[ii]].QuantumIdentity(),
              qi.Species(), qi.Isotopologue(), qi.LowerQuantumNumbers(),qi.UpperQuantumNumbers()))
              continue;
          
            Numeric x_constant = 0.0, y_constant = 0.0;
            switch(flag_partials[flag_partials_position[ii]].PropMatType())
            {
              case JacPropMatType::LineGammaSelfExp: 
                    x_constant = dpsf_dSelfExponent;
                    y_constant = dgamma_dSelfExponent; 
                    break;
              case JacPropMatType::LineGammaForeignExp:
                    x_constant = dpsf_dForeignExponent;
                    y_constant = dgamma_dForeignExponent; 
                    break;
              case JacPropMatType::LineGammaWaterExp: 
                    x_constant = dpsf_dWaterExponent;
                    y_constant = dgamma_dWaterExponent; 
                    break;
              case JacPropMatType::LineShiftSelf:
                    x_constant = dpsf_dSelf;
                    break;
              case JacPropMatType::LineShiftForeign:
                    x_constant = dpsf_dForeign;
                    break;
              case JacPropMatType::LineShiftWater:
                    x_constant = dpsf_dWater;
                    break;
              default: throw std::runtime_error("This cannot happen.  A developer has made a mistake.\n");
            }
            
            
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = do_partials_phase?partials_phase[ii](this_f_grid, pressure_level_index):empty_vector;
            VectorView this_partial_src         = do_src?partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            Numeric dP_dgamma, dF_dpsf;
            Vector dfn_dpsf(nv);
            // Calculate the line shape derivative:
            global_data::lineshape_data[ind_ls].dInput_dgamma()(dP_dgamma,sigma);
            global_data::lineshape_data[ind_ls].dInput_dF0()(dF_dpsf,sigma);
            global_data::lineshape_norm_data[ind_lsn].dFunction_dF0()(dfn_dpsf, f0, this_f, temperature);
            
            dP_dgamma *= y_constant;
            dF_dpsf *= x_constant;
            
            for(Index iv=0;iv<nv;iv++)
            {   
                calc_derivative(this_partial_attenuation[iv],
                                this_partial_phase[iv],
                                this_partial_src[iv],
                                dFa_dx[iv], dFb_dx[iv], dFa_dy[iv], dFb_dy[iv],
                                dF_dpsf, dP_dgamma, // dx_dtarget, dy_target
                                LM_Fa[iv], LM_Fb[iv], // Full calculations
                                dfn_dpsf[iv]/C[iv]*x_constant, 0.0, // ratios
                                CF_A[iv], CF_B[iv], 0.0, 0.0, 0.0, 0.0, //  When attenuation and phase matters
                                nlte, 0.0,
                                do_partials_phase, do_src);
            }   
            // That's it!
        }
        else if(is_line_mixing_line_strength_parameter(flag_partials[flag_partials_position[ii]]))
        {
          if(!line_match_line(flag_partials[flag_partials_position[ii]].QuantumIdentity(),
            qi.Species(), qi.Isotopologue(), qi.LowerQuantumNumbers(), qi.UpperQuantumNumbers()))
                continue;
            
            Numeric lma_real = 0.0, lmb_real = 0.0, lma_imag = 0.0, lmb_imag = 0.0;
            switch(flag_partials[flag_partials_position[ii]].PropMatType())
            {
                case JacPropMatType::LineMixingY0: 
                    lmb_real = dY_LM0;
                    lma_imag = dY_LM0; 
                    break;
                case JacPropMatType::LineMixingG0: 
                    lmb_imag = dG_LM0;
                    lma_real = dG_LM0; 
                    break;
                case JacPropMatType::LineMixingY1: 
                    lmb_real = dY_LM1;
                    lma_imag = dY_LM1;
                    break;
                case JacPropMatType::LineMixingG1: 
                    lmb_real = dG_LM1;
                    lma_imag = dG_LM1;
                    break;
                case JacPropMatType::LineMixingYExp: 
                    lmb_real = dY_LMexp;
                    lma_imag = dY_LMexp;
                    break;
                case JacPropMatType::LineMixingGExp: 
                    lmb_real = dG_LMexp;
                    lma_imag = dG_LMexp;
                    break;
                default: throw std::runtime_error("This cannot happen.  A developer has made a mistake.\n");
            }
            
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = do_partials_phase?partials_phase[ii](this_f_grid, pressure_level_index):empty_vector;
            VectorView this_partial_src         = do_src?partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            for(Index iv=0;iv<nv;iv++)
            {
                calc_derivative(this_partial_attenuation[iv],
                                this_partial_phase[iv],
                                this_partial_src[iv],
                                dFa_dx[iv], dFb_dx[iv], dFa_dy[iv], dFb_dy[iv],
                                0.0, 0.0, // dx_dtarget, dy_target
                                LM_Fa[iv], LM_Fb[iv], // Full calculations
                                0.0, 0.0, // ratios
                                CF_A[iv], CF_B[iv], 
                                lma_real, lma_imag, lmb_real, lmb_imag, //  When attenuation and phase matters
                                nlte, 0.0,
                                do_partials_phase, do_src);
            }
            
            // That's it!
        }
        else if(is_line_mixing_DF_parameter(flag_partials[flag_partials_position[ii]]))
        {
          if(!line_match_line(flag_partials[flag_partials_position[ii]].QuantumIdentity(),
            qi.Species(), qi.Isotopologue(),qi.LowerQuantumNumbers(),qi.UpperQuantumNumbers()))
                continue;
            
            Numeric constant = 1.0;
            switch(flag_partials[flag_partials_position[ii]].PropMatType())
            {
              case JacPropMatType::LineMixingDF0: constant = dDF_LM0; break;
              case JacPropMatType::LineMixingDF1: constant = dDF_LM1; break;
              case JacPropMatType::LineMixingDFExp: constant = dDF_LMexp; break;
              default: throw std::runtime_error("This cannot happen.  A developer has made a mistake.\n");
            }
            
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = do_partials_phase?partials_phase[ii](this_f_grid, pressure_level_index):empty_vector;
            VectorView this_partial_src         = do_src?partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            Numeric dF_dDF;
            
            // Calculate the line shape derivative:
            global_data::lineshape_data[ind_ls].dInput_dDF()(dF_dDF,sigma);
            
            for(Index iv=0;iv<nv;iv++)
            {
                
                this_partial_attenuation[iv] += dFa_dx[iv] * dF_dDF * constant;
                if(do_partials_phase)
                    this_partial_phase[iv]   += dFb_dx[iv] * dF_dDF * constant;
                if(do_src)
                    this_partial_src[iv]     += dFa_dx[iv] * dF_dDF * nlte * constant;
            }
            
            // That's it!  Note that the output should be strongly correlated to Temperature and Pressure.
            // Also note that to do the fit for the catalog gamma is somewhat different than this
        }
        else if(flag_partials[flag_partials_position[ii]] == JacPropMatType::NLTE)
        {
            
            /* 
             * WARNING:  This part will be in accordance with simplified formalism used in the transfer code
             * so that d/dTnlte[exp(f(Tnlte))/Q(T)], where Q(T) is independent of Tnlte.  In practice, Q(T)
             * will be Q(T)-exp(f(T))+exp(f(Tnlte)), and so there should be an additional term involved 
             * in the calculation of this partial derivative.  However, as this is not required in the simplified
             * formalism we presently include in ARTS, this is also not accounted for below.  I am not sure what
             * this implies for usability of these partial derivatives.
             */
            
            bool lower, upper;
            line_match_level(lower, upper, 
                             flag_partials[flag_partials_position[ii]].QuantumIdentity(),
                             qi.Species(), qi.Isotopologue(),
                             qi.LowerQuantumNumbers(), qi.UpperQuantumNumbers());
            if(!(lower||upper))
                continue;
            
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = do_partials_phase?partials_phase[ii](this_f_grid, pressure_level_index):empty_vector;
            VectorView this_partial_src         = do_src?partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            const Numeric Gamma = exp(-PLANCK_CONST * line_frequency / BOLTZMAN_CONST / temperature);
            
            Numeric dK4_dTu=0., dK3_dTu=0., dK3_dTl=0.;
            
            if(upper)
            {
                dK4_dTu =   K4*line_E_v_upp/line_T_v_upp/line_T_v_upp;
                dK3_dTu = - K4*line_E_v_upp/line_T_v_upp/line_T_v_upp/BOLTZMAN_CONST*Gamma/(Gamma - 1.0);
            }
            
            if(lower)
            {
                dK3_dTl = - K4*line_E_v_low/line_T_v_low/line_T_v_low/BOLTZMAN_CONST*Gamma/(Gamma - 1.0);
                //dK4_dTl = 0.0;
            }
            
            const Numeric dK_dTl = nlte*dK3_dTl/K3; // div K3 since there is a K3 in ls_A later that should not be there...
            const Numeric dK_dTu = nlte*dK3_dTu/K3;
            
            const Numeric dsourceC_dTl = - dK3_dTl/K3 ; // Simplifying some terms
            const Numeric dsourceC_dTu = ( dK4_dTu - dK3_dTu ) / K3 ;
            
            for(Index iv=0;iv<nv;iv++)
            {
                
                const Numeric ls_A  = ( (1.0 + G_LM)*CF_A[iv] + Y_LM*CF_B[iv]), 
                              ls_B  = ( (1.0 + G_LM)*CF_B[iv] - Y_LM*CF_A[iv]);
                
                this_partial_src[iv] += (dsourceC_dTl+dsourceC_dTu) * ls_A; //Note that this works since
                this_partial_attenuation[iv] += (dK_dTl+dK_dTu) * ls_A;     //Tu and Tl are independently
                if(do_partials_phase)                                       //changing absorption...
                    this_partial_phase[iv] += (dK_dTl+dK_dTu) * ls_B;       //Also, they are 0 when inactive...
            }
        }
    }
}



bool line_match_line(const QuantumIdentifier& from_jac,
                     const Index& species,
                     const Index& isotopologue,
                     const QuantumNumbers& lower_qn, 
                     const QuantumNumbers& upper_qn)
{
  if(species not_eq from_jac.Species() or isotopologue not_eq from_jac.Isotopologue())
  {
    return false;
  }
  else
  {
    QuantumMatchInfoEnum lower, upper;
    if(from_jac.Type() == QuantumIdentifier::TRANSITION)
    {
      lower_qn.CompareDetailed(lower, from_jac.QuantumMatch()[from_jac.TRANSITION_LOWER_INDEX]);
      upper_qn.CompareDetailed(upper, from_jac.QuantumMatch()[from_jac.TRANSITION_UPPER_INDEX]);
      
      return (lower==QMI_FULL&&upper==QMI_FULL); //Must be perfect match for this to be true.
    }
    else if(from_jac.Type() == QuantumIdentifier::ALL)
    {
      return true;
    }
    else 
    {
      return false;
    }
  }
}

void line_match_level(bool& lower_energy_level,
                      bool& upper_energy_level,
                      const QuantumIdentifier& from_jac,
                      const Index& species,
                      const Index& isotopologue,
                      const QuantumNumbers& lower_qn, 
                      const QuantumNumbers& upper_qn)
{
  if(species not_eq from_jac.Species() or isotopologue not_eq from_jac.Isotopologue())
  {
    lower_energy_level = false;
    upper_energy_level = false;
  } 
  else
  {
    // These can be partial, because vibrational energy levels are partial matches, so using simpler function
    lower_energy_level = lower_qn.Compare(from_jac.QuantumMatch()[from_jac.TRANSITION_UPPER_INDEX]);
    upper_energy_level = upper_qn.Compare(from_jac.QuantumMatch()[from_jac.TRANSITION_UPPER_INDEX]);
  }
}
