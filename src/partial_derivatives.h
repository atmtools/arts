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
extern const String NLTE_MAINTAG;
extern const String WIND_MAINTAG;
extern const String MAGFIELD_MAINTAG;
extern const String FREQUENCY_MAINTAG;
extern const String ABSSPECIES_MAINTAG;
extern const String ELECTRONS_MAINTAG;
extern const String PARTICULATES_MAINTAG;
extern const String SCATSPECIES_MAINTAG;
//extern const String PRESSURE_MAINTAG

// Main tag for this type of Jacobian calculation
extern const String PROPMAT_SUBSUBTAG;

// Catalog parameters
extern const String CATALOGPARAMETER_MAINTAG;

// Generic modes
extern const String LINESTRENGTH_MODE;
extern const String LINECENTER_MODE;

//  Pressure Broadening
extern const String SELFBROADENING_MODE;
extern const String FOREIGNBROADENING_MODE;
extern const String WATERBROADENING_MODE;
extern const String SELFBROADENINGEXPONENT_MODE;
extern const String FOREIGNBROADENINGEXPONENT_MODE;
extern const String WATERBROADENINGEXPONENT_MODE;
extern const String SELFPRESSURESHIFT_MODE;
extern const String FOREIGNPRESSURESHIFT_MODE;
extern const String WATERPRESSURESHIFT_MODE;

//  Line Mixing
extern const String LINEMIXINGY0_MODE;
extern const String LINEMIXINGG0_MODE;
extern const String LINEMIXINGDF0_MODE;
extern const String LINEMIXINGY1_MODE;
extern const String LINEMIXINGG1_MODE;
extern const String LINEMIXINGDF1_MODE;
extern const String LINEMIXINGYEXPONENT_MODE;
extern const String LINEMIXINGGEXPONENT_MODE;
extern const String LINEMIXINGDFEXPONENT_MODE;


typedef enum 
{
    JQT_VMR,
    JQT_electrons,
    JQT_particulates,
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
    JQT_line_pressureshift_self,
    JQT_line_pressureshift_foreign,
    JQT_line_pressureshift_water,
    JQT_line_mixing_Y,
    JQT_line_mixing_G,
    JQT_line_mixing_DF,
    JQT_line_mixing_Y0,
    JQT_line_mixing_G0,
    JQT_line_mixing_DF0,
    JQT_line_mixing_Y1,
    JQT_line_mixing_G1,
    JQT_line_mixing_DF1,
    JQT_line_mixing_Yexp,
    JQT_line_mixing_Gexp,
    JQT_line_mixing_DFexp,
    JQT_nlte_temperature,
    JQT_population_level_ratio,
    JQT_NOT_JQT
} JacobianQuantityType;



typedef Array<JacobianQuantityType> ArrayOfJacobianQuantityType;

// Helper internal class that is only used to define and know what propagation matrix partials can be calculated by propmat_clearsky agenda
class PropmatPartialsData
{
public:
  
    // Empty init
    PropmatPartialsData()
    {
      mreal_nelem = 0;
      mcontains_temperature = false;
      mcontains_frequency_term = false;
      mcontains_linemixing_0_term = false;
      mcontains_linemixing_1_term = false;
      mcontains_linemixing_exponent_term = false;
      mmagnetic_u = false;
      mmagnetic_v = false; 
      mmagnetic_w = false;
      mmagnetic_abs = false;
      mmagnetic_theta = false;
      mmagnetic_eta = false;
    }
    
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
        mcontains_frequency_term = false;
        mcontains_linemixing_0_term = false;
        mcontains_linemixing_1_term = false;
        mcontains_linemixing_exponent_term = false;
        mmagnetic_u = false;
        mmagnetic_v = false; 
        mmagnetic_w = false;
        mmagnetic_abs = false;
        mmagnetic_theta = false;
        mmagnetic_eta = false;
        
        mmag_perturbation  = 0.0;
        mtemp_perturbation = 0.0;
        mfreq_perturbation = 0.0;
        
        mjacobian_quantities = jacobian_quantities;
        
        mqtype.resize(mreal_nelem);
        mjacobian_pos.resize(mreal_nelem);
        mspecies.resize(mreal_nelem);
        
        mcontains_pressure_term.resize(9);
        for(Index ii=0;ii<9;ii++)
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
                else if(jacobian_quantities[iq].MainTag() == ELECTRONS_MAINTAG )
                {
                  mqtype[ippdq] = JQT_electrons;
                  mjacobian_pos[ippdq] = iq;
                  mspecies[ippdq]  = -9999;//Flag for not a species...
                  ippdq++;
                }
                else if(jacobian_quantities[iq].MainTag() == PARTICULATES_MAINTAG )
                {
                  mqtype[ippdq] = JQT_particulates;
                  mjacobian_pos[ippdq] = iq;
                  mspecies[ippdq]  = -9999;//Flag for not a species...
                  ippdq++;
                }
                else if(jacobian_quantities[iq].MainTag() == MAGFIELD_MAINTAG)
                {
                  if(jacobian_quantities[iq].Subtag() == "strength")
                  {
                      mqtype[ippdq] = JQT_magnetic_magntitude;
                      mmagnetic_abs = true;
                      mjacobian_pos[ippdq] = iq;
                      mspecies[ippdq] = -9999;//Flag for not a species...
                      ippdq++;
                  }
                  else if(jacobian_quantities[iq].Subtag() == "eta")
                  {
                      mqtype[ippdq] = JQT_magnetic_eta;
                      mmagnetic_eta = true;
                      mjacobian_pos[ippdq] = iq;
                      mspecies[ippdq] = -9999;//Flag for not a species...
                      ippdq++;
                  }
                  else if(jacobian_quantities[iq].Subtag() == "theta")
                  {
                      mqtype[ippdq] = JQT_magnetic_theta;
                      mmagnetic_theta = true;
                      mjacobian_pos[ippdq] = iq;
                      mspecies[ippdq] = -9999;//Flag for not a species...
                      ippdq++;
                  }
                  else if(jacobian_quantities[iq].Subtag() == "u")
                  {
                      mqtype[ippdq] = JQT_magnetic_u;  
                      mmagnetic_u = true;
                      mjacobian_pos[ippdq] = iq; 
                      mspecies[ippdq] = -9999;//Flag for not a species...
                      ippdq++;
                      mmag_perturbation = jacobian_quantities[iq].Perturbation();
                  }
                  else if(jacobian_quantities[iq].Subtag() == "v")
                  {
                      mqtype[ippdq] = JQT_magnetic_v;
                      mmagnetic_v = true;
                      mjacobian_pos[ippdq] = iq;
                      mspecies[ippdq] = -9999;//Flag for not a species...
                      ippdq++;
                      mmag_perturbation = jacobian_quantities[iq].Perturbation();
                  }
                  else if(jacobian_quantities[iq].Subtag() == "w")
                  {
                      mqtype[ippdq] = JQT_magnetic_w;
                      mmagnetic_w = true;
                      mjacobian_pos[ippdq] = iq;
                      mspecies[ippdq] = -9999;//Flag for not a species...
                      ippdq++;
                      mmag_perturbation = jacobian_quantities[iq].Perturbation();
                  }
                }
                else if(jacobian_quantities[iq].MainTag() == WIND_MAINTAG)
                {
                  if(jacobian_quantities[iq].Subtag() == "strength")
                  {
                      mqtype[ippdq] = JQT_wind_magnitude;
                      mcontains_frequency_term = true;
                      mjacobian_pos[ippdq] = iq;
                      mspecies[ippdq] = -9999;//Flag for not a species...
                      ippdq++;
                      mfreq_perturbation = jacobian_quantities[iq].Perturbation();
                  }
                  else if(jacobian_quantities[iq].Subtag() == "u")
                  {
                      mqtype[ippdq] = JQT_wind_u;
                      mcontains_frequency_term = true;
                      mjacobian_pos[ippdq] = iq;
                      mspecies[ippdq] = -9999;//Flag for not a species...
                      ippdq++;
                      mfreq_perturbation = jacobian_quantities[iq].Perturbation();
                  }
                  else if(jacobian_quantities[iq].Subtag() == "v")
                  {
                      mqtype[ippdq] = JQT_wind_v;
                      mcontains_frequency_term = true;
                      mjacobian_pos[ippdq] = iq;
                      mspecies[ippdq] = -9999;//Flag for not a species...
                      ippdq++;
                      mfreq_perturbation = jacobian_quantities[iq].Perturbation();
                  }
                  else if(jacobian_quantities[iq].Subtag() == "w")
                  {
                      mqtype[ippdq] = JQT_wind_w;
                      mcontains_frequency_term = true;
                      mjacobian_pos[ippdq] = iq;
                      mspecies[ippdq] = -9999;//Flag for not a species...
                      ippdq++;
                      mfreq_perturbation = jacobian_quantities[iq].Perturbation();
                  }
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
                    mcontains_frequency_term = true;
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
                    /* else if(jacobian_quantities[iq].Mode() == PRESSUREBROADENINGGAMMA_MODE)
                    {
                        mqtype[ippdq] = JQT_line_gamma;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    } */
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
                    else if(jacobian_quantities[iq].Mode() == SELFPRESSURESHIFT_MODE)
                    {
                        mqtype[ippdq] = JQT_line_pressureshift_self;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        mcontains_pressure_term[6]=1;
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == FOREIGNPRESSURESHIFT_MODE)
                    {
                        mqtype[ippdq] = JQT_line_pressureshift_foreign;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        mcontains_pressure_term[7]=1;
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == WATERPRESSURESHIFT_MODE)
                    {
                        mqtype[ippdq] = JQT_line_pressureshift_water;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        mcontains_pressure_term[8]=1;
                        ippdq++;
                    }
                    /*else if(jacobian_quantities[iq].Mode() == LINEMIXINGY_MODE)
                    {
                        mqtype[ippdq] = JQT_line_mixing_Y;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == LINEMIXINGG_MODE)
                    {
                        mqtype[ippdq] = JQT_line_mixing_G;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == LINEMIXINGDF_MODE)
                    {
                        mqtype[ippdq] = JQT_line_mixing_DF;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }*/
                    else if(jacobian_quantities[iq].Mode() == LINEMIXINGY0_MODE)
                    {
                        mqtype[ippdq] = JQT_line_mixing_Y0;
                        mcontains_linemixing_0_term = true;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == LINEMIXINGG0_MODE)
                    {
                        mqtype[ippdq] = JQT_line_mixing_G0;
                        mcontains_linemixing_0_term = true;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == LINEMIXINGDF0_MODE)
                    {
                        mqtype[ippdq] = JQT_line_mixing_DF0;
                        mcontains_linemixing_0_term = true;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == LINEMIXINGY1_MODE)
                    {
                        mqtype[ippdq] = JQT_line_mixing_Y1;
                        mcontains_linemixing_1_term = true;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == LINEMIXINGG1_MODE)
                    {
                        mqtype[ippdq] = JQT_line_mixing_G1;
                        mcontains_linemixing_1_term = true;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == LINEMIXINGDF1_MODE)
                    {
                        mqtype[ippdq] = JQT_line_mixing_DF1;
                        mcontains_linemixing_1_term = true;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == LINEMIXINGYEXPONENT_MODE)
                    {
                        mqtype[ippdq] = JQT_line_mixing_Yexp;
                        mcontains_linemixing_exponent_term  = true;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == LINEMIXINGGEXPONENT_MODE)
                    {
                        mqtype[ippdq] = JQT_line_mixing_Gexp;
                        mcontains_linemixing_exponent_term  = true;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == LINEMIXINGDFEXPONENT_MODE)
                    {
                        mqtype[ippdq] = JQT_line_mixing_DFexp;
                        mcontains_linemixing_exponent_term = true;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                    else if(jacobian_quantities[iq].Mode() == NLTE_MAINTAG)
                    {
                        mqtype[ippdq] = JQT_nlte_temperature;
                        mjacobian_pos[ippdq] = iq;
                        mspecies[ippdq] = -9999;//Flag for not a species...
                        ippdq++;
                    }
                }
                else
                    throw std::runtime_error("Unknown Jacobian quantity for propagation calculations.\n");
            }
        }
        set_first_frequency();
        set_first_pressure_term();
    }
    
    void SetOnlyTemperatureTrue()
    {
      // Only meant for testing purposes
      mreal_nelem = 1;
      mqtype.resize(1);
      mqtype[0] = JQT_temperature;
      mjacobian_pos.resize(1);
      mjacobian_pos[0] = 0;
      mspecies.resize(1);
      mspecies[0] = -9999;//Flag for not a species...
      mcontains_temperature = true;
      mtemp_perturbation = 0.1;
    }
    
    // Returns the length of Jacobians that are computed inside the propagation agenda
    Index nelem() const {return mreal_nelem;}
    
    // Returns the type of the jacobian at the index (index < nelem())
    const JacobianQuantityType&  operator()(const Index& iq) const {return mqtype[iq];}
    
    // Returns the full jacobian ArrayOfRetrievalQuantity
    const ArrayOfRetrievalQuantity& jac() const {return mjacobian_quantities;}
    
    // Returns the jacobian RetrievalQuantity at the index (index < nelem())
    const RetrievalQuantity& jac(const Index& iq) const {return mjacobian_quantities[mjacobian_pos[iq]];}
    
    // Returns the species at the index
    Index species(Index iq) const { return mspecies[iq]; }
    
    // Run this check in CIA calculations to find if there is anything wrong with it
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
    
    // Run this check in Continuum calculations to find if there is anything wrong with it
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
    
    // Run this check in line-by-line calculations without phase to find if there is anything wrong with it
    bool supportsLBLwithoutPhase() const 
    {
        
        for(Index iq=0; iq<nelem(); iq++)
        {
            if( mqtype[iq]!=JQT_VMR || 
                mqtype[iq]!=JQT_line_strength ||
                mqtype[iq]!=JQT_nlte_temperature)
                return false;
        }
        
        return true;
    }
    
    // Run this check in Relaxation Matrix calculations to find if there is anything wrong with it
    bool supportsRelaxationMatrix() const 
    {
        bool testvar = false;
        
        for(Index iq=0; iq<nelem(); iq++)
        {
            if( mqtype[iq]==JQT_temperature             || //Supported types by relmat
                mqtype[iq]==JQT_frequency               || 
                mqtype[iq]==JQT_wind_magnitude          || 
                mqtype[iq]==JQT_wind_u                  || 
                mqtype[iq]==JQT_wind_v                  || 
                mqtype[iq]==JQT_wind_w)
                testvar=true;
            else if(mqtype[iq]==JQT_line_center    || // Not yet supported
                mqtype[iq]==JQT_line_gamma     ||
                mqtype[iq]==JQT_line_gamma_self     ||
                mqtype[iq]==JQT_line_gamma_foreign     || 
                mqtype[iq]==JQT_line_gamma_water     || 
                mqtype[iq]==JQT_line_gamma_selfexponent     || 
                mqtype[iq]==JQT_line_gamma_foreignexponent     ||
                mqtype[iq]==JQT_line_gamma_waterexponent     ||
                mqtype[iq]==JQT_line_mixing_DF || 
                mqtype[iq]==JQT_line_mixing_G  ||
                mqtype[iq]==JQT_line_mixing_Y  || 
                mqtype[iq]==JQT_line_strength)
                throw std::runtime_error("Line specific parameters are not supported while\n"
                "using the relaxation matrix line mixing routine.\n"
                "We do not yet track individual lines in the relaxation matrix calculations.\n");
        }
        return testvar;
    };
    
    // Run this check in Lookup calculations to find if there is anything wrong with it
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
                mqtype[iq]==JQT_line_gamma     || // line is in the Lookup?
                mqtype[iq]==JQT_line_gamma_self     || // line is in the Lookup?
                mqtype[iq]==JQT_line_gamma_foreign     || // line is in the Lookup?
                mqtype[iq]==JQT_line_gamma_water     || // line is in the Lookup?
                mqtype[iq]==JQT_line_gamma_selfexponent     || // line is in the Lookup?
                mqtype[iq]==JQT_line_gamma_foreignexponent     || // line is in the Lookup?
                mqtype[iq]==JQT_line_gamma_waterexponent     || // line is in the Lookup?
                mqtype[iq]==JQT_line_mixing_DF || // Did not add JQT_level_vibrational_temperature,
                mqtype[iq]==JQT_line_mixing_G  || // since no continuum model takes NLTE into account.
                mqtype[iq]==JQT_line_mixing_Y  || // (if they do, add this to list of unsupported Jacobians.)
                mqtype[iq]==JQT_line_strength)
                throw std::runtime_error("Line specific parameters are not supported while using Lookup "
                "table.\nWe do not track lines in the Lookup.\n");
        }
        return testvar;
    };
    
    // Run this check in Zeeman calculations to find if there is anything wrong with it
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
    
    // Run this check in Zeeman precalc calculations to find if there is anything wrong with it
    bool supportsZeemanPrecalc() const 
    {
        for(Index iq=0; iq<nelem(); iq++) 
            if(mqtype[iq]!=JQT_NOT_JQT)  // All is supported
                return true; 
        return false;
    }
    
    // Run this check in Faraday calculations to find if there is anything wrong with it
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
                    mqtype[iq]==JQT_wind_w                 ||
                    mqtype[iq]==JQT_electrons)
            {
                testvar=true;
            }
        }
        return testvar;
    }
    
    // Run this check in Particulates calculations to find if there is anything wrong with it
    bool supportsParticles() const
    {
        for(Index iq=0; iq<nelem(); iq++)
          if(mqtype[iq]!=JQT_NOT_JQT)
            return true;
        return false;
    }
    
    // Run this check to see if species is in jacobian list calculations to find if there is anything wrong with it
    bool supportsPropmatClearsky(const Index this_species) const
    {   
        for(Index iq=0; iq<nelem(); iq++) 
            if( species(iq)!=this_species )  // Species flag is just for speed
                return true;
        return false;
    }
    
    // Returns the counter to the jacobian position
    Index this_jq_index(Index& jqt_index) const
    {
        Index counter = -1;
        for(Index iq=0; iq<=jqt_index; iq++)
            if(mjacobian_quantities[iq].SubSubtag()==PROPMAT_SUBSUBTAG)
                counter++;
        return counter;
    };
    
    // Returns if this type should be computed
    bool is_this_propmattype(Index& jqt_index) const
    {
        const Index this_index = this_jq_index(jqt_index);
        if(this_index==-1)
            return mjacobian_quantities[jqt_index].MainTag()==CATALOGPARAMETER_MAINTAG;
        return mjacobian_quantities[jqt_index].MainTag()==CATALOGPARAMETER_MAINTAG ||
        mqtype[this_index]==JQT_magnetic_magntitude||mqtype[this_index]==JQT_magnetic_theta||mqtype[this_index]==JQT_magnetic_eta;
    }
    
    // Note that wind and frequency means the same inside propmat_clearsky agenda.
    bool do_frequency() const
    {
        for(Index iq=0; iq<nelem(); iq++)
        {
          if( mqtype[iq]==JQT_frequency || mqtype[iq]==JQT_wind_magnitude || mqtype[iq]==JQT_wind_u || mqtype[iq]==JQT_wind_v || mqtype[iq]==JQT_wind_w)
          {
            return true;
          }
        }
        return false;
    }
    
    bool do_temperature() const {return mcontains_temperature;}
    
    // Returns if any line center jacobains are needed
    bool do_line_center() const
    {
      for(Index iq=0; iq<nelem(); iq++)
      {
        if( mqtype[iq] == JQT_line_center)
        {
          return true;
        }
      }
      return false;
    }
    
    // Returns if any magnetic jacobian calculations are needed
    bool do_magnetic_field() const
    {
        for(Index iq=0; iq<nelem(); iq++)
            if(mqtype[iq]==JQT_magnetic_eta || mqtype[iq]==JQT_magnetic_magntitude || mqtype[iq]==JQT_magnetic_theta || mqtype[iq]==JQT_magnetic_u || mqtype[iq]==JQT_magnetic_v || mqtype[iq]==JQT_magnetic_w)
                return true;
        return false;
    }
    
    // Zeeman tests
    bool do_zeeman_u() const {return mmagnetic_u;};
    bool do_zeeman_v() const {return mmagnetic_v;};
    bool do_zeeman_w() const {return mmagnetic_w;};
    bool do_zeeman_abs() const {return mmagnetic_abs;};
    bool do_zeeman_theta() const {return mmagnetic_theta;};
    bool do_zeeman_eta() const {return mmagnetic_eta;};
    
    // Setting default values here so that the same perturbations are used everywhere where perturbations are required
    Numeric Temperature_Perturbation() const {return mtemp_perturbation;};   
    Numeric Magnetic_Field_Perturbation() const {return mmag_perturbation;}; 
    Numeric Frequency_Perturbation() const {return mfreq_perturbation;};  

    // Helper for pressure broadening terms
    Index PressureBroadeningTerm(Index this_index) const {return (bool)mcontains_pressure_term[this_index];};
    bool ZerothTermLM() const {return mcontains_linemixing_0_term;};
    bool FirstTermLM() const {return mcontains_linemixing_1_term;};
    bool ExponentLM() const {return mcontains_linemixing_exponent_term;};
    
    void set_first_frequency()
    {
      mfirst_frequency = -1;
      for(Index iq=0; iq<nelem(); iq++)
      {
        if(mqtype[iq] == JQT_frequency or
          mqtype[iq] == JQT_wind_magnitude or
          mqtype[iq] == JQT_wind_u or 
          mqtype[iq] == JQT_wind_v or 
          mqtype[iq]==JQT_wind_w)
        {
          mfirst_frequency = iq;
          return;
        }
      }
    }
    
    void set_first_pressure_term()
    {
      mfirst_pressure = -1;
      for(Index iq=0; iq<nelem(); iq++)
      {
        if(mqtype[iq] == JQT_line_gamma_self or 
          mqtype[iq] == JQT_line_gamma_selfexponent or
          mqtype[iq] == JQT_line_pressureshift_self or
          mqtype[iq] == JQT_line_gamma_foreign or
          mqtype[iq] == JQT_line_gamma_foreignexponent or
          mqtype[iq] == JQT_line_pressureshift_foreign or
          mqtype[iq] == JQT_line_gamma_water or
          mqtype[iq] == JQT_line_gamma_waterexponent or 
          mqtype[iq] == JQT_line_pressureshift_water)
        {
          mfirst_pressure = iq;
          return;
        }
      }
    }
    
    Index get_first_frequency() const {return mfirst_frequency;}
    Index get_first_pressure_term() const {return mfirst_pressure;}
    
    // Only for debug purposes
    String StringTypeAtIndex(Index ii) const
    {
      switch(mqtype[ii])
      {
        case JQT_VMR:
          return "VMR";
        case JQT_electrons:
          return "Electrons-VMR";
        case JQT_particulates:
          return "Particulate-VMR";
        case JQT_temperature:
          return "Temperature";
        case JQT_magnetic_magntitude:
          return "Magnetic-Strength";
        case JQT_magnetic_eta:
          return "Magnetic-Eta";
        case JQT_magnetic_theta:
          return "Magnetic-Theta";
        case JQT_magnetic_u:
          return "Magnetic-u";
        case JQT_magnetic_v:
          return "Magnetic-v";
        case JQT_magnetic_w:
          return "Magnetic-w";
        case JQT_wind_magnitude:
          return "Wind-Strength";
        case JQT_wind_u:
          return "Wind-u";
        case JQT_wind_v:
          return "Wind-v";
        case JQT_wind_w:
          return "Wind-w";
        case JQT_frequency:
          return "Frequency";
        case JQT_line_strength:
          return "Line-Strength";
        case JQT_line_center:
          return "Line-Center";
        case JQT_line_gamma:
          return "Line-Gamma";
        case JQT_line_gamma_self:
          return "Line-Gamma-Self";
        case JQT_line_gamma_foreign:
          return "Line-Gamma-Foreign";
        case JQT_line_gamma_water:
          return "Line-Gamma-Water";
        case JQT_line_gamma_selfexponent:
          return "Line-Gamma-SelfExponent";
        case JQT_line_gamma_foreignexponent:
          return "Line-Gamma-ForeignExponent";
        case JQT_line_gamma_waterexponent:
          return "Line-Gamma-WaterExponent";
        case JQT_line_pressureshift_self:
          return "Line-PressureShift-Self";
        case JQT_line_pressureshift_foreign:
          return "Line-PressureShift-Foreign";
        case JQT_line_pressureshift_water:
          return "Line-PressureShift-Water";
        case JQT_line_mixing_Y:
          return "Line-Mixing-Y";
        case JQT_line_mixing_G:
          return "Line-Mixing-G";
        case JQT_line_mixing_DF:
          return "Line-Mixing-DF";
        case JQT_line_mixing_Y0:
          return "Line-Mixing-Y0";
        case JQT_line_mixing_G0:
          return "Line-Mixing-G0";
        case JQT_line_mixing_DF0:
          return "Line-Mixing-DF0";
        case JQT_line_mixing_Y1:
          return "Line-Mixing-Y1";
        case JQT_line_mixing_G1:
          return "Line-Mixing-G1";
        case JQT_line_mixing_DF1:
          return "Line-Mixing-DF1";
        case JQT_line_mixing_Yexp:
          return "Line-Mixing-YExp";
        case JQT_line_mixing_Gexp:
          return "Line-Mixing-GExp";
        case JQT_line_mixing_DFexp:
          return "Line-Mixing-DFExp";
        case JQT_nlte_temperature:
          return "NLTE-Temperatures";
        case JQT_NOT_JQT:
          return "Not-A-PropMat-Variable";
        default:
          throw std::runtime_error("Dev forgot to add this...");
      }
    }
    
private:
    ArrayOfJacobianQuantityType mqtype;
    ArrayOfIndex                mjacobian_pos;
    ArrayOfIndex                mspecies;
    ArrayOfRetrievalQuantity    mjacobian_quantities;
    Numeric                     mtemp_perturbation;
    Numeric                     mmag_perturbation;
    Numeric                     mfreq_perturbation;
    Index                       mreal_nelem;
    Index                       mfirst_frequency;
    Index                       mfirst_pressure;
    bool mcontains_temperature;
    bool mcontains_frequency_term;
    bool mcontains_linemixing_0_term;
    bool mcontains_linemixing_1_term;
    bool mcontains_linemixing_exponent_term;
    bool mmagnetic_u;
    bool mmagnetic_v;
    bool mmagnetic_w;
    bool mmagnetic_abs;
    bool mmagnetic_theta;
    bool mmagnetic_eta;
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
                                              const Index& species,
                                              const Index& isotopologue,
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
                                              const Numeric& dQ_dT,
                                              // Magnetic variables
                                              const Numeric&  DF_Zeeman,
                                              const Numeric&  H_mag_Zeeman,
                                              const bool      do_zeeman,
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
                     const Index& species,
                     const Index& isotopologue,
                     const QuantumNumbers& lower_qn, 
                     const QuantumNumbers& upper_qn);

void line_match_level(bool& lower_energy_level,
                      bool& upper_energy_level,
                      const QuantumIdentifier& from_jac,
                      const Index& species,
                      const Index& isotopologue,
                      const QuantumNumbers& lower_qn, 
                      const QuantumNumbers& upper_qn);

#endif
