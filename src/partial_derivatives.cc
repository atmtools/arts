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


//This function will require much more inputs in the future
void partial_derivatives_lineshape_dependency(ArrayOfMatrix&  partials_attenuation,
                                              ArrayOfMatrix&  partials_phase, 
                                              ArrayOfMatrix&  partials_src, 
                                              const PropmatPartialsData&  flag_partials, 
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
                                              const QuantumNumberRecord&  qnr,
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
                                              const Numeric& dgamma_dSelfExponent,
                                              const Numeric& dgamma_dForeignExponent,
                                              const Numeric& dgamma_dWaterExponent,
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
    
    Vector empty_vector;
    
    if(this_f_grid.get_extent()==0)
        return;
    
    ConstVectorView this_f = f_grid[this_f_grid];
    const Index nv = this_f.nelem();
    const Numeric nlte = do_src? K4/K3-1.0 : 0.0;
    const Numeric f0 = line_frequency + df_0 + DF_LM + DF_Zeeman;
    
    // Loop over all jacobian_quantities, if a matching quantity is found, then apply the necessary steps to make the jacobian matrix
    for(Index ii=0; ii<flag_partials.nelem(); ii++)
    {
        if( (flag_partials(ii)==JQT_magnetic_magntitude) )
        {
            if(!do_zeeman)
                continue;
            else if( H_mag_Zeeman==0.0 )
                throw std::runtime_error("Sorry, but it is not supported to get Zeeman derivative"
                "when magnetic field is 0 T.\nIf required, let the devs know.\n");
            
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = partials_phase[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_src         = do_src?partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            const Numeric dF_dH = DF_Zeeman/H_mag_Zeeman;
            
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
        else if(flag_partials(ii)==JQT_magnetic_theta)
        {/* Pass on this.  Can be done easier in later parts of the code execution.  
            Will not require line shape partials. */}
        else if(flag_partials(ii)==JQT_magnetic_eta)
        {/* Pass on this.  Can be done easier in later parts of the code execution.  
            Will not require line shape partials. */}
        else if((flag_partials(ii)==JQT_magnetic_u)            ||
                (flag_partials(ii)==JQT_magnetic_v)            || 
                (flag_partials(ii)==JQT_magnetic_w))
        {/* Pass on these.  These are done by perturbation in later parts of the code execution since I found it too complicated for now.  FIXME: Richard */}
        else if((flag_partials(ii)==JQT_frequency)              || // This is ready
                (flag_partials(ii)==JQT_wind_magnitude)         || // Note: this requires one more step! dF/dW is missing
                (flag_partials(ii)==JQT_wind_u)                 || // Note: these three
                (flag_partials(ii)==JQT_wind_v)                 || // terms still need
                (flag_partials(ii)==JQT_wind_w))                   // two more steps! dF/dW and dW/d{u,v,w} are missing
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
        else if(flag_partials(ii)==JQT_temperature)
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
        else if(flag_partials(ii) == JQT_VMR)
        {
            // Line shape cross section does not depend on VMR.  NOTE:  Ignoring self-pressure broadening.
        }
        else if(flag_partials(ii)==JQT_line_center)
        {
            if(!line_match_line(flag_partials.jac()[ii].QuantumIdentity(),qnr.Lower(),qnr.Upper()))
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
            global_data::lineshape_data[ind_ls].dInput_dF0()(dF_dF0,sigma);
            
            Vector dfn_dF0(nv);
            global_data::lineshape_norm_data[ind_lsn].dFunction_dF0()(dfn_dF0, f0,
                                                                      this_f, temperature);
            
            for(Index iv=0;iv<nv;iv++)
            {
                const Numeric ls_A= ( (1.0 + G_LM)*CF_A[iv] + Y_LM*CF_B[iv]), 
                              ls_B= ( (1.0 + G_LM)*CF_B[iv] - Y_LM*CF_A[iv]);
                
                this_partial_attenuation[iv] += (dS_dF0+dfn_dF0[iv])*ls_A + dFa_dx[iv] * dF_dF0;
                if(do_partials_phase)
                    this_partial_phase[iv]   += (dS_dF0+dfn_dF0[iv])*ls_B + dFb_dx[iv] * dF_dF0;
                if(do_src)
                    this_partial_src[iv]     += ((dS_dF0+dfn_dF0[iv])*ls_A + dFa_dx[iv] * dF_dF0) * nlte - ls_A*K4/K3/K3*dK3_dF0;
            }
            
            // That's it!  Note that the output should be strongly correlated to Temperature and Pressure.
            // Also note that to do the fit for the catalog gamma is somewhat different than this
        }
        else if(flag_partials(ii) == JQT_line_strength)
        {
            if(!line_match_line(flag_partials.jac()[ii].QuantumIdentity(),qnr.Lower(),qnr.Upper()))
                continue;
            
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = do_partials_phase?partials_phase[ii](this_f_grid, pressure_level_index):empty_vector;
            VectorView this_partial_src         = do_src?partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            for(Index iv=0;iv<nv;iv++)
            {
                const Numeric ls_A= ( (1.0 + G_LM)*CF_A[iv] + Y_LM*CF_B[iv]), 
                              ls_B= ( (1.0 + G_LM)*CF_B[iv] - Y_LM*CF_A[iv]);
                            
                this_partial_attenuation[iv] += ls_A/line_strength;
                if(do_partials_phase)
                    this_partial_phase[iv]   -= ls_B/line_strength;
                if(do_src)
                    this_partial_src[iv]     += ls_A/line_strength*nlte;
                    
            }
            
            // That's it!  Now to wonder if this will grow extremely large, creating an unrealistic jacobian...
        }
        else if(flag_partials(ii) == JQT_line_gamma                 ||
                flag_partials(ii) == JQT_line_gamma_self            ||
                flag_partials(ii) == JQT_line_gamma_foreign         ||
                flag_partials(ii) == JQT_line_gamma_water           ||
                flag_partials(ii) == JQT_line_gamma_selfexponent    ||
                flag_partials(ii) == JQT_line_gamma_foreignexponent ||
                flag_partials(ii) == JQT_line_gamma_waterexponent)
        {
            if(!line_match_line(flag_partials.jac()[ii].QuantumIdentity(),qnr.Lower(),qnr.Upper()))
                continue;
            
            Numeric constant = 1.0;
            switch(flag_partials(ii))
            {
                case JQT_line_gamma: break;
                case JQT_line_gamma_self: constant = dgamma_dSelf; break;
                case JQT_line_gamma_foreign: constant = dgamma_dForeign; break;
                case JQT_line_gamma_water: constant = dgamma_dWater; break;
                case JQT_line_gamma_selfexponent: constant = dgamma_dSelfExponent; break;
                case JQT_line_gamma_foreignexponent: constant = dgamma_dForeignExponent; break;
                case JQT_line_gamma_waterexponent: constant = dgamma_dWaterExponent; break;
                default: throw std::runtime_error("This cannot happen.  A developer has made a mistake.\n");
            }
            
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = do_partials_phase?partials_phase[ii](this_f_grid, pressure_level_index):empty_vector;
            VectorView this_partial_src         = do_src?partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            Numeric dP_dgamma;
            // Calculate the line shape derivative:
            global_data::lineshape_data[ind_ls].dInput_dgamma()(dP_dgamma,sigma);
            
            for(Index iv=0;iv<nv;iv++)
            {
                
                this_partial_attenuation[iv] += dFa_dy[iv]*dP_dgamma * constant;
                if(do_partials_phase)
                    this_partial_phase[iv]   += dFb_dy[iv]*dP_dgamma * constant;
                if(do_src)
                    this_partial_src[iv]     += dFa_dy[iv]*dP_dgamma*nlte * constant;
                    
            }
            
            // That's it!
        }
        else if(flag_partials(ii)==JQT_line_mixing_Y  ||
                flag_partials(ii)==JQT_line_mixing_Y0 ||
                flag_partials(ii)==JQT_line_mixing_Y1 ||
                flag_partials(ii)==JQT_line_mixing_Yexp)
        {
            if(!line_match_line(flag_partials.jac()[ii].QuantumIdentity(),qnr.Lower(),qnr.Upper()))
                continue;
            
            Numeric constant = 1.0;
            switch(flag_partials(ii))
            {
                case JQT_line_mixing_Y: break;
                case JQT_line_mixing_Y0: constant = dY_LM0; break;
                case JQT_line_mixing_Y1: constant = dY_LM1; break;
                case JQT_line_mixing_Yexp: constant = dY_LMexp; break;
                default: throw std::runtime_error("This cannot happen.  A developer has made a mistake.\n");
            }
            
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = do_partials_phase?partials_phase[ii](this_f_grid, pressure_level_index):empty_vector;
            VectorView this_partial_src         = do_src?partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            for(Index iv=0;iv<nv;iv++)
            {
                
                this_partial_attenuation[iv] += CF_B[iv] * constant;
                if(do_partials_phase)
                    this_partial_phase[iv]   -= CF_A[iv] * constant;
                if(do_src)
                    this_partial_src[iv]     += CF_B[iv]*nlte * constant;
            }
            
            // That's it!  Note that this should be strongly correlated to Temperature and Pressure.
            // Also note that to do the fit for the catalog gamma is somewhat different than this
        }
        else if(flag_partials(ii)==JQT_line_mixing_G  ||
                flag_partials(ii)==JQT_line_mixing_G0 ||
                flag_partials(ii)==JQT_line_mixing_G1 ||
                flag_partials(ii)==JQT_line_mixing_Gexp)
        {
            if(!line_match_line(flag_partials.jac()[ii].QuantumIdentity(),qnr.Lower(),qnr.Upper()))
                continue;
            
            Numeric constant = 1.0;
            switch(flag_partials(ii))
            {
                case JQT_line_mixing_G: break;
                case JQT_line_mixing_G0: constant = dG_LM0; break;
                case JQT_line_mixing_G1: constant = dG_LM1; break;
                case JQT_line_mixing_Gexp: constant = dG_LMexp; break;
                default: throw std::runtime_error("This cannot happen.  A developer has made a mistake.\n");
            }
            
            VectorView this_partial_attenuation = partials_attenuation[ii](this_f_grid, pressure_level_index);
            VectorView this_partial_phase       = do_partials_phase?partials_phase[ii](this_f_grid, pressure_level_index):empty_vector;
            VectorView this_partial_src         = do_src?partials_src[ii](this_f_grid, pressure_level_index):empty_vector;
            
            for(Index iv=0;iv<nv;iv++)
            {
                
                this_partial_attenuation[iv] += CF_A[iv] * constant;
                if(do_partials_phase)
                    this_partial_phase[iv]   += CF_B[iv] * constant;
                if(do_src)
                    this_partial_src[iv]     += CF_A[iv]*nlte * constant;
            }
            
            // That's it!  Note that the output should be strongly correlated to Temperature and Pressure.
            // Also note that to do the fit for the catalog gamma is somewhat different than this
        }
        else if(flag_partials(ii)==JQT_line_mixing_DF  ||
                flag_partials(ii)==JQT_line_mixing_DF0 ||
                flag_partials(ii)==JQT_line_mixing_DF1 ||
                flag_partials(ii)==JQT_line_mixing_DFexp)
        {
            if(!line_match_line(flag_partials.jac()[ii].QuantumIdentity(),qnr.Lower(),qnr.Upper()))
                continue;
            
            Numeric constant = 1.0;
            switch(flag_partials(ii))
            {
                case JQT_line_mixing_DF: break;
                case JQT_line_mixing_DF0: constant = dDF_LM0; break;
                case JQT_line_mixing_DF1: constant = dDF_LM1; break;
                case JQT_line_mixing_DFexp: constant = dDF_LMexp; break;
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
        else if(flag_partials(ii)==JQT_nlte_temperature)
        {
            bool lower, upper;
            line_match_level(lower, upper, flag_partials.jac()[ii].QuantumIdentity(), qnr.Lower(), qnr.Upper());
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
                     const QuantumNumbers& lower_qn, 
                     const QuantumNumbers& upper_qn)
{
    QuantumMatchInfoEnum lower, upper;
    
    lower_qn.CompareDetailed(lower, from_jac.QuantumMatch()[from_jac.TRANSITION_LOWER_INDEX]);
    upper_qn.CompareDetailed(upper, from_jac.QuantumMatch()[from_jac.TRANSITION_UPPER_INDEX]);
    
    return (lower==QMI_FULL&&upper==QMI_FULL); //Must be perfect match for this to be true.
}

void line_match_level(bool& lower_energy_level,
                      bool& upper_energy_level,
                      const QuantumIdentifier& from_jac,
                      const QuantumNumbers& lower_qn, 
                      const QuantumNumbers& upper_qn)
{
    // These can be partial, because vibrational energy levels are partial matches, so using simpler function
    lower_energy_level = lower_qn.Compare(from_jac.QuantumMatch()[from_jac.TRANSITION_UPPER_INDEX]);
    upper_energy_level = upper_qn.Compare(from_jac.QuantumMatch()[from_jac.TRANSITION_UPPER_INDEX]);
}