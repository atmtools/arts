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


#ifndef propmatpartials_h
#define propmatpartials_h

#include "jacobian.h"

// Generic parameters
extern const String TEMPERATURE_MAINTAG;
extern const String WIND_MAINTAG;
extern const String MAGFIELD_MAINTAG;
extern const String FREQUENCY_MAINTAG;
extern const String ABSSPECIES_MAINTAG;

// Main tag for this type of Jacobian calculation
extern const String PROPMAT_SUBSUBTAG;

// Catalog parameters
extern const String CATALOGPARAMETER_MAINTAG;
extern const String PRESSUREBROADENINGGAMMA_MODE;
extern const String LINESTRENGTH_MODE;
extern const String LINECENTER_MODE;

// Specific catalog's parameters
extern const String SELFBROADENING_MODE;
extern const String FOREIGNBROADENING_MODE;
extern const String WATERBROADENING_MODE;
extern const String SELFBROADENINGEXPONENT_MODE;
extern const String FOREIGNBROADENINGEXPONENT_MODE;
extern const String WATERBROADENINGEXPONENT_MODE;


typedef enum 
{
    JQT_VMR,
    JQT_temperature,
    JQT_magnetic_magntitude,
    JQT_magnetic_eta,
    JQT_magnetic_theta,
    JQT_magnetic_u,
    JQT_magnetic_v,
    JQT_magnetic_w,
    JQT_wind_magnitude, // towards sensor?
    JQT_wind_u,
    JQT_wind_v,
    JQT_wind_w,
    JQT_frequency,
    JQT_line_strength,
    JQT_line_center,
    JQT_line_gamma,
    JQT_line_gamma_self,
    JQT_line_gamma_foreign,
    JQT_line_gamma_water,
    JQT_line_gamma_selfexponent,
    JQT_line_gamma_foreignexponent,
    JQT_line_gamma_waterexponent,
    JQT_line_mixing_Y,
    JQT_line_mixing_G,
    JQT_line_mixing_DF,
    JQT_level_vibrational_temperature,
    JQT_NOT_JQT
} JacobianQuantityType;



typedef Array<JacobianQuantityType> ArrayOfJacobianQuantityType;

// Helper internal class that is only used to define and know what propagation matrix partials can be calculated by propmat_clearsky agenda
class PropmatPartialsData
{
public:
    
    // Reduces jacobian_quantities down to only relevant quantities for the propmat_clearsky
    PropmatPartialsData(const ArrayOfRetrievalQuantity& jacobian_quantities)
    {
        const Index nq=jacobian_quantities.nelem();
        mreal_nelem = nq;
        
        for(Index iq=0; iq<nq; iq++)
        {
            if(jacobian_quantities[iq].SubSubtag()!=PROPMAT_SUBSUBTAG)
                mreal_nelem--;
        }
        
        mcontains_temperature = false;
        
        mmag_perturbation  = 0.0;
        mtemp_perturbation = 0.0;
        mfreq_perturbation = 0.0;
        
        mjacobian_quantities = jacobian_quantities;
        
        mqtype.resize(mreal_nelem);
        mjacobian_pos.resize(mreal_nelem);
        mspecies.resize(mreal_nelem);
        
        mcontains_pressure_term.resize(6);
        for(Index ii=0;ii<6;ii++)
            mcontains_pressure_term[ii]=0;
        
        Index ippdq = 0;
        for(Index iq=0; iq<nq; iq++)
        {
            if(jacobian_quantities[iq].SubSubtag()==PROPMAT_SUBSUBTAG)
            {
                if(jacobian_quantities[iq].MainTag() == ABSSPECIES_MAINTAG )
                {
                    mqtype[ippdq] = JQT_VMR;
                    mjacobian_pos[ippdq] = iq;
                    mspecies[ippdq]  = SpeciesTag(jacobian_quantities[iq].Subtag()).Species();
                    ippdq++;
                    mmag_perturbation = jacobian_quantities[iq].Perturbation();
                }
                else if(jacobian_quantities[iq].MainTag() == MAGFIELD_MAINTAG && 
                    jacobian_quantities[iq].Subtag() == "strength")
                {
                    mqtype[ippdq] = JQT_magnetic_magntitude;
                    mjacobian_pos[ippdq] = iq;
                    mspecies[ippdq] = -9999;//Flag for not a species...
                    ippdq++;
                }
                else if(jacobian_quantities[iq].MainTag() == MAGFIELD_MAINTAG && 
                    jacobian_quantities[iq].Subtag() == "eta")
                {
                    mqtype[ippdq] = JQT_magnetic_eta;
                    mjacobian_pos[ippdq] = iq;
                    mspecies[ippdq] = -9999;//Flag for not a species...
                    ippdq++;
                }
                else if(jacobian_quantities[iq].MainTag() == MAGFIELD_MAINTAG && 
                    jacobian_quantities[iq].Subtag() == "theta")
                {
                    mqtype[ippdq] = JQT_magnetic_theta;
                    mjacobian_pos[ippdq] = iq;
                    mspecies[ippdq] = -9999;//Flag for not a species...
                    ippdq++;
                }
                else if(jacobian_quantities[iq].MainTag() == MAGFIELD_MAINTAG && 
                    jacobian_quantities[iq].Subtag() == "u")
                {
                    mqtype[ippdq] = JQT_magnetic_u;  
                    mjacobian_pos[ippdq] = iq; 
                    mspecies[ippdq] = -9999;//Flag for not a species...
                    ippdq++;
                    mmag_perturbation = jacobian_quantities[iq].Perturbation();
                }
                else if(jacobian_quantities[iq].MainTag() == MAGFIELD_MAINTAG && 
                    jacobian_quantities[iq].Subtag() == "v")
                {
                    mqtype[ippdq] = JQT_magnetic_v;
                    mjacobian_pos[ippdq] = iq;
                    mspecies[ippdq] = -9999;//Flag for not a species...
                    ippdq++;
                    mmag_perturbation = jacobian_quantities[iq].Perturbation();
                }
                else if(jacobian_quantities[iq].MainTag() == MAGFIELD_MAINTAG && 
                    jacobian_quantities[iq].Subtag() == "w")
                {
                    mqtype[ippdq] = JQT_magnetic_w;
                    mjacobian_pos[ippdq] = iq;
                    mspecies[ippdq] = -9999;//Flag for not a species...
                    ippdq++;
                    mmag_perturbation = jacobian_quantities[iq].Perturbation();
                }
                else if(jacobian_quantities[iq].MainTag() == WIND_MAINTAG &&  
                    jacobian_quantities[iq].Subtag() == "strength")
                {
                    mqtype[ippdq] = JQT_wind_magnitude;
                    mjacobian_pos[ippdq] = iq;
                    mspecies[ippdq] = -9999;//Flag for not a species...
                    ippdq++;
                    mfreq_perturbation = jacobian_quantities[iq].Perturbation();
                }
                else if(jacobian_quantities[iq].MainTag() == WIND_MAINTAG &&  
                    jacobian_quantities[iq].Subtag() == "u")
                {
                    mqtype[ippdq] = JQT_wind_u;
                    mjacobian_pos[ippdq] = iq;
                    mspecies[ippdq] = -9999;//Flag for not a species...
                    ippdq++;
                    mfreq_perturbation = jacobian_quantities[iq].Perturbation();
                }
                else if(jacobian_quantities[iq].MainTag() == WIND_MAINTAG &&  
                    jacobian_quantities[iq].Subtag() == "v")
                {
                    mqtype[ippdq] = JQT_wind_v;
                    mjacobian_pos[ippdq] = iq;
                    mspecies[ippdq] = -9999;//Flag for not a species...
                    ippdq++;
                    mfreq_perturbation = jacobian_quantities[iq].Perturbation();
                }
                else if(jacobian_quantities[iq].MainTag() == WIND_MAINTAG &&  
                    jacobian_quantities[iq].Subtag() == "w")
                {
                    mqtype[ippdq] = JQT_wind_w;
                    mjacobian_pos[ippdq] = iq;
                    mspecies[ippdq] = -9999;//Flag for not a species...
                    ippdq++;
                    mfreq_perturbation = jacobian_quantities[iq].Perturbation();
                }
                else if(jacobian_quantities[iq].MainTag() == TEMPERATURE_MAINTAG)
                {
                    mqtype[ippdq] = JQT_temperature;
                    mjacobian_pos[ippdq] = iq;
                    mspecies[ippdq] = -9999;//Flag for not a species...
                    ippdq++;
                    mcontains_temperature = true;
                    mtemp_perturbation = jacobian_quantities[iq].Perturbation();
                }
                else if(jacobian_quantities[iq].MainTag() == FREQUENCY_MAINTAG)
                {
                    mqtype[ippdq] = JQT_frequency;
                    mjacobian_pos[ippdq] = iq;
                    mspecies[ippdq] = -9999;//Flag for not a species...
                    ippdq++;
                    mfreq_perturbation = jacobian_quantities[iq].Perturbation();
                }
                else if(jacobian_quantities[iq].MainTag() == CATALOGPARAMETER_MAINTAG)
                {
                    if(jacobian_quantities[iq].Mode() == LINESTRENGTH_MODE)
                    {
                        mqtype[ippdq] = JQT_line_strength;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == LINECENTER_MODE)
                    {
                        mqtype[ippdq] = JQT_line_center;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == PRESSUREBROADENINGGAMMA_MODE)
                    {
                        mqtype[ippdq] = JQT_line_gamma;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == SELFBROADENING_MODE)
                    {
                        mqtype[ippdq] = JQT_line_gamma_self;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        mcontains_pressure_term[0]=1;
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == FOREIGNBROADENING_MODE)
                    {
                        mqtype[ippdq] = JQT_line_gamma_foreign;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        mcontains_pressure_term[1]=1;
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == WATERBROADENING_MODE)
                    {
                        mqtype[ippdq] = JQT_line_gamma_water;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        mcontains_pressure_term[2]=1;
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == SELFBROADENINGEXPONENT_MODE)
                    {
                        mqtype[ippdq] = JQT_line_gamma_selfexponent;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        mcontains_pressure_term[3]=1;
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == FOREIGNBROADENINGEXPONENT_MODE)
                    {
                        mqtype[ippdq] = JQT_line_gamma_foreignexponent;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        mcontains_pressure_term[4]=1;
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == WATERBROADENINGEXPONENT_MODE)
                    {
                        mqtype[ippdq] = JQT_line_gamma_waterexponent;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        mcontains_pressure_term[5]=1;
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == "Line Mixing Y")
                    {
                        mqtype[ippdq] = JQT_line_mixing_Y;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == "Line Mixing G")
                    {
                        mqtype[ippdq] = JQT_line_mixing_G;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == "Line Mixing DF")
                    {
                        mqtype[ippdq] = JQT_line_mixing_DF;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == "Vibrational Temperature")
                    {
                        mqtype[ippdq] = JQT_level_vibrational_temperature;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                }
                else
                    throw std::runtime_error("Unknown Jacobian quantity for propagation calculations.\n");
            }
        }
    }
    
    Index nelem() const {return mreal_nelem;}
    bool do_temperature() const {return mcontains_temperature;}
    
    const JacobianQuantityType&  operator()(const Index& iq) const {return mqtype[iq];}
    const ArrayOfRetrievalQuantity& jac() const {return mjacobian_quantities;}
    
    Index species(Index iq) const { return mspecies[iq]; }
    
    bool supportsCIA() const 
    {
        bool testvar = false;
        
        for(Index iq=0; iq<nelem(); iq++)
        {
            if( mqtype[iq]==JQT_temperature             || //Supported types by CIA
                mqtype[iq]==JQT_frequency               || 
                mqtype[iq]==JQT_wind_magnitude          || 
                mqtype[iq]==JQT_wind_u                  || 
                mqtype[iq]==JQT_wind_v                  || 
                mqtype[iq]==JQT_wind_w                  || 
                mqtype[iq]==JQT_VMR)
                testvar=true;
            //else if(false) //There are no unsupported CIA variables.  Thus no runtime_error to throw
        }
        return testvar;
    };
    
    bool supportsContinuum() const 
    {
        bool testvar = false;
        
        for(Index iq=0; iq<nelem(); iq++)
        {
            if( mqtype[iq]==JQT_temperature             || //Supported types by Continuum
                mqtype[iq]==JQT_frequency               || 
                mqtype[iq]==JQT_wind_magnitude          || 
                mqtype[iq]==JQT_wind_u                  || 
                mqtype[iq]==JQT_wind_v                  || 
                mqtype[iq]==JQT_wind_w                  || 
                mqtype[iq]==JQT_VMR)
                testvar=true;
            else if(mqtype[iq]==JQT_line_center    || // We cannot know if any particular
                mqtype[iq]==JQT_line_gamma     || // line is in the Continuum?
                mqtype[iq]==JQT_line_gamma_self     || // line is in the Continuum?
                mqtype[iq]==JQT_line_gamma_foreign     || // line is in the Continuum?
                mqtype[iq]==JQT_line_gamma_water     || // line is in the Continuum?
                mqtype[iq]==JQT_line_gamma_selfexponent     || // line is in the Continuum?
                mqtype[iq]==JQT_line_gamma_foreignexponent     || // line is in the Continuum?
                mqtype[iq]==JQT_line_gamma_waterexponent     || // line is in the Continuum?
                    mqtype[iq]==JQT_line_mixing_DF || // Did not add JQT_level_vibrational_temperature,
                    mqtype[iq]==JQT_line_mixing_G  || // since no continuum model takes NLTE into account.
                    mqtype[iq]==JQT_line_mixing_Y  || // (if they do, add this to list of unsupported Jacobians.)
                    mqtype[iq]==JQT_line_strength)
                throw std::runtime_error("Line specific parameters are not supported while using continuum"
                " tags.\nWe do not track what lines are in the continuum.\n");
        }
        return testvar;
    };
    
    bool supportsLBLwithoutPhase() const 
    {
        
        for(Index iq=0; iq<nelem(); iq++)
        {
            if( mqtype[iq]!=JQT_VMR || 
                mqtype[iq]!=JQT_line_strength ||
                mqtype[iq]!=JQT_level_vibrational_temperature)
                return false;
        }
        
        return true;
    }
    
    bool supportsLookup() const 
    {
        bool testvar = false;
        
        for(Index iq=0; iq<nelem(); iq++)
        {
            if( mqtype[iq]==JQT_temperature             || //Supported types by Lookup
                mqtype[iq]==JQT_frequency               || 
                mqtype[iq]==JQT_wind_magnitude          || 
                mqtype[iq]==JQT_wind_u                  || 
                mqtype[iq]==JQT_wind_v                  || 
                mqtype[iq]==JQT_wind_w                  || 
                mqtype[iq]==JQT_VMR)
                testvar=true;
            else if(mqtype[iq]==JQT_line_center    || // We cannot know if any particular
                mqtype[iq]==JQT_line_gamma     || // line is in the Continuum?
                mqtype[iq]==JQT_line_gamma_self     || // line is in the Continuum?
                mqtype[iq]==JQT_line_gamma_foreign     || // line is in the Continuum?
                mqtype[iq]==JQT_line_gamma_water     || // line is in the Continuum?
                mqtype[iq]==JQT_line_gamma_selfexponent     || // line is in the Continuum?
                mqtype[iq]==JQT_line_gamma_foreignexponent     || // line is in the Continuum?
                mqtype[iq]==JQT_line_gamma_waterexponent     || // line is in the Continuum?
                mqtype[iq]==JQT_line_mixing_DF || // Did not add JQT_level_vibrational_temperature,
                mqtype[iq]==JQT_line_mixing_G  || // since no continuum model takes NLTE into account.
                mqtype[iq]==JQT_line_mixing_Y  || // (if they do, add this to list of unsupported Jacobians.)
                mqtype[iq]==JQT_line_strength)
                throw std::runtime_error("Line specific parameters are not supported while using Lookup "
                "table.\nWe do not track lines in the Lookup.\n");
        }
        return testvar;
    };
    
    bool supportsZeeman() const 
    {
        bool testvar = false;
        
        for(Index iq=0; iq<nelem(); iq++)
        {
            if( mqtype[iq]==JQT_magnetic_eta        || // Not Supported types by this Zeeman function
                mqtype[iq]==JQT_magnetic_magntitude || 
                mqtype[iq]==JQT_magnetic_theta      || 
                mqtype[iq]==JQT_magnetic_u          || 
                mqtype[iq]==JQT_magnetic_v          || 
                mqtype[iq]==JQT_magnetic_w)
                throw std::runtime_error("This method does not support magnetic Jacobian calculations.\n"
                    "Please use the Precalc method instead.\n");
            else if(mqtype[iq]!=JQT_NOT_JQT)
                testvar = true;
        }
        return testvar;
    };
    
    bool supportsZeemanPrecalc() const 
    {
        for(Index iq=0; iq<nelem(); iq++) 
            if(mqtype[iq]!=JQT_NOT_JQT)  // All is supported
                return true; 
        return false;
    }
    
    bool supportsFaraday() const
    {
        bool testvar = false;
        for(Index iq=0; iq<nelem(); iq++)
        {
            if( mqtype[iq]==JQT_magnetic_eta        || // Not Supported types by Faraday
                mqtype[iq]==JQT_magnetic_magntitude || //←could be supported?
                mqtype[iq]==JQT_magnetic_theta)
            {
                throw std::runtime_error("This method does not yet support Zeeman-style magnetic Jacobian calculations.\n"
                                         "Please use u, v, and w Jacobians instead.\n");
            }
            else if(mqtype[iq]==JQT_magnetic_u             || // Supported by Faraday
                    mqtype[iq]==JQT_magnetic_v             || 
                    mqtype[iq]==JQT_magnetic_w             ||
                    mqtype[iq]==JQT_frequency              ||
                    mqtype[iq]==JQT_wind_magnitude         ||
                    mqtype[iq]==JQT_wind_u                 || 
                    mqtype[iq]==JQT_wind_v                 || 
                    mqtype[iq]==JQT_wind_w)
            {
                testvar=true;
            }
        }
        return testvar;
    }
    
    bool supportsParticles() const
    {
        for(Index iq=0; iq<nelem(); iq++) 
            if(mqtype[iq]!=JQT_NOT_JQT)  // None is supported
                throw std::runtime_error("You cannot combine propmat_clearsky Jacobians with this method.\n"); 
        return false;
    }
    
    bool supportsPropmatClearsky(const Index this_species) const
    {   
        for(Index iq=0; iq<nelem(); iq++) 
            if( species(iq)!=this_species )  // Species flag is just for speed
                return true;
        return false;
    }
    
    Index this_jq_index(Index& jqt_index) const
    {
        Index counter = -1;
        for(Index iq=0; iq<=jqt_index; iq++)
            if(mjacobian_quantities[iq].SubSubtag()==PROPMAT_SUBSUBTAG)
                counter++;
        return counter;
    };
    
    bool is_this_linetype(Index& jqt_index) const
    {
        return mjacobian_quantities[jqt_index].MainTag()==CATALOGPARAMETER_MAINTAG;
    }
    
    // Note that wind and frequency means the same inside propmat_clearsky agenda.
    bool do_frequency() const
    {
        for(Index iq=0; iq<nelem(); iq++)
            if( mqtype[iq]==JQT_frequency || mqtype[iq]==JQT_wind_magnitude || mqtype[iq]==JQT_wind_u || mqtype[iq]==JQT_wind_v || mqtype[iq]==JQT_wind_w)
                return true;
        return false;
    }
    
    bool do_magnetic_field() const
    {
        for(Index iq=0; iq<nelem(); iq++)
            if(mqtype[iq]==JQT_magnetic_eta || mqtype[iq]==JQT_magnetic_magntitude || mqtype[iq]==JQT_magnetic_theta || mqtype[iq]==JQT_magnetic_u || mqtype[iq]==JQT_magnetic_v || mqtype[iq]==JQT_magnetic_w)
                return true;
        return false;
    }
    
    // Setting default values here so that the same perturbations are used everywhere where perturbations are required
    Numeric Temperature_Perturbation() const {return mtemp_perturbation;};   
    Numeric Magnetic_Field_Perturbation() const {return mmag_perturbation;}; 
    Numeric Frequency_Perturbation() const {return mfreq_perturbation;};  

    // Helper for pressure broadening terms
    Index PressureBroadeningTerm(Index this_index) const {return (bool)mcontains_pressure_term[this_index];};
    
private:
    ArrayOfJacobianQuantityType mqtype;
    ArrayOfIndex                mjacobian_pos;
    ArrayOfIndex                mspecies;
    ArrayOfRetrievalQuantity    mjacobian_quantities;
    Numeric                     mtemp_perturbation;
    Numeric                     mmag_perturbation;
    Numeric                     mfreq_perturbation;
    Index                       mreal_nelem;
    bool mcontains_temperature;
    bool mcontains_frequency_term;
    ArrayOfIndex mcontains_pressure_term;
};

// Generic function that will calculate all partial derivatives of the line-shape that depends on input per line
void partial_derivatives_lineshape_dependency(ArrayOfMatrix& partials_attenuation,
                                              ArrayOfMatrix& partials_phase, 
                                              ArrayOfMatrix&  partials_src,  
                                              const PropmatPartialsData&  flag_partials, 
                                              ConstVectorView CF_A,
                                              ConstVectorView CF_B,
                                              ConstVectorView C,
                                              ConstVectorView dFa_dF,
                                              ConstVectorView dFb_dF, 
                                              ConstVectorView dFa_dP, 
                                              ConstVectorView dFb_dP, 
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
                                              const Numeric&  G_LM,
                                              const Numeric&  dG_LM_dT,
                                              const Numeric&  DF_LM,
                                              const Numeric&  dDF_LM_dT,
                                              const QuantumNumberRecord&  qnr,
                                              // LINE SHAPE
                                              const Index& ind_ls,
                                              const Index& ind_lsn,
                                              const Numeric& df_0,
                                              const Numeric& ddf_0_dT,
                                              const Numeric& gamma,
                                              const Numeric& dgamma_dT,
                                              const Numeric& dgamma_dSelf,
                                              const Numeric& dgamma_dForeign,
                                              const Numeric& dgamma_dWater,
                                              const Numeric& dgamma_dSelfExponent,
                                              const Numeric& dgamma_dForeignExponent,
                                              const Numeric& dgamma_dWaterExponent,
                                              // Partition data parameters
                                              const Numeric& dQ_dT,
                                              // Magnetic variables
                                              const Numeric&  DF_Zeeman,
                                              const Numeric&  H_mag_Zeeman,
                                              // Programming variables
                                              const Index&    pressure_level_index,
                                              const bool      do_partials_phase,
                                              const bool      do_src);

void partial_derivatives_of_stokes_along_path(ArrayOfMatrix& dIn_dt,  //Size is [n_retrieval_points][f_grid,stokes_dim]
                                              const ArrayOfRetrievalQuantity& jacobian_quantities,
                                              ConstTensor4View                extinction_matrices,
                                              const ArrayOfTensor4&           partial_propmat_clearsky,
                                              const Numeric&                  path_length,
                                              const Numeric&                  atmospheric_temperature);

bool line_match_line(const QuantumIdentifier& from_jac,
                     const QuantumNumbers& lower_qn, 
                     const QuantumNumbers& upper_qn);

void line_match_level(bool& lower_energy_level,
                      bool& upper_energy_level,
                      const QuantumIdentifier& from_jac,
                      const QuantumNumbers& lower_qn, 
                      const QuantumNumbers& upper_qn);

#endif