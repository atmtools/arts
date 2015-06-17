/* Copyright (C) 2014 
 R *ichard Larsson <ric.larsson@gmail.com>
 
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

/** Contains additional functions for the pressure broadening data class
 * \file   pressurebroadeningdata.cc
 * 
 * \author Richard Larsson
 * \date   2014-11-06
 **/
#include "pressurebroadeningdata.h"
#include "linerecord.h"


///////////////////////////////////////////
//  Broadening calculations below here
//////////////////////////////////////////

// Get broadening
void PressureBroadeningData::GetPressureBroadeningParams(Numeric& gamma,
                                                         Numeric& deltaf,
                                                         const Numeric& theta,
                                                         const Numeric& pressure,
                                                         const Numeric& self_pressure,
                                                         const Index    this_species,
                                                         const Index    h2o_species,
                                                         const ArrayOfIndex& broad_spec_locations,
                                                         ConstVectorView vmrs,
                                                         const Verbosity& verbosity) const
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            GetAirBroadening(gamma, deltaf, theta, pressure, self_pressure);
            break;
        case PB_AIR_AND_WATER_BROADENING:
            GetAirAndWaterBroadening(gamma, deltaf, theta, pressure, self_pressure, 
                                     this_species, h2o_species, vmrs, verbosity);
            break;
        case PB_PERRIN_BROADENING:
            GetPerrinBroadening(gamma, deltaf, theta, pressure, self_pressure,
                                this_species, broad_spec_locations, vmrs, verbosity);
            break;
        default:
            throw std::runtime_error("You have defined an unknown broadening mechanism.\n");
    }
}

// Get air broadening
void PressureBroadeningData::GetAirBroadening(Numeric& gamma,
                                              Numeric& deltaf,
                                              const Numeric& theta,
                                              const Numeric& pressure,
                                              const Numeric& self_pressure) const
{
    gamma   =
      mdata[2][0] * pow(theta,mdata[3][0]) * (pressure-self_pressure)
      + mdata[0][0] * pow(theta,mdata[1][0]) *         self_pressure;
    deltaf  =
      mdata[4][0] * pressure
      * pow (theta,(Numeric)0.25+(Numeric)1.5*mdata[3][0]);
}


// Get air broadening
void PressureBroadeningData::GetAirAndWaterBroadening(Numeric& gamma,
                                                      Numeric& deltaf,
                                                      const Numeric& theta,
                                                      const Numeric& pressure,
                                                      const Numeric& self_pressure,
                                                      const Index    this_species,
                                                      const Index    h2o_species,
                                                      ConstVectorView vmrs,
                                                      const Verbosity& verbosity) const
{
    CREATE_OUT2;
    
    if(this_species==h2o_species)    
    {
        gamma   =
        mdata[1][0] * pow(theta,mdata[1][1]) * (pressure-self_pressure)
        + mdata[0][0] * pow(theta,mdata[0][1]) *         self_pressure;
        
        deltaf  =
        mdata[1][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[1][1]) * (pressure-self_pressure)
        + mdata[0][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[0][1]) *         self_pressure;
    }
    else if(h2o_species==-1)
    {
        gamma   =
        mdata[1][0] * pow(theta,mdata[1][1]) * (pressure-self_pressure)
        + mdata[0][0] * pow(theta,mdata[0][1]) *         self_pressure;
        
        deltaf  =
        mdata[1][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[1][1]) * (pressure-self_pressure)
        + mdata[0][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[0][1]) *         self_pressure;
        
        out2 << "You have no H2O in species but you use water-broadened line shape.\n";
    }
    else
    {
        gamma   =
        mdata[1][0] * pow(theta,mdata[1][1]) * (pressure-self_pressure-vmrs[h2o_species]*pressure)
        + mdata[0][0] * pow(theta,mdata[0][1]) * self_pressure
        + mdata[2][0] * pow(theta,mdata[2][1]) * vmrs[h2o_species]*pressure;
        
        deltaf  =
        mdata[1][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[1][1]) * (pressure-self_pressure-vmrs[h2o_species]*pressure)
        + mdata[0][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[0][1]) * self_pressure
        + mdata[2][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[2][1]) * vmrs[h2o_species]*pressure;
    }
}


// Get Perrin broadening
void PressureBroadeningData::GetPerrinBroadening(Numeric& gamma,
                                                 Numeric& deltaf,
                                                 const Numeric& theta,
                                                 const Numeric& pressure,
                                                 const Numeric& self_pressure,
                                                 const Index    this_species,
                                                 const ArrayOfIndex& broad_spec_locations,
                                                 ConstVectorView vmrs,
                                                 const Verbosity& verbosity) const
{
    CREATE_OUT2;
    
    // Number of broadening species:
    const Index nbs = LineRecord::NBroadSpec();
    assert(nbs==broad_spec_locations.nelem());
    
    // Split total pressure in self and foreign part:
    const Numeric foreign_pressure = pressure - self_pressure;
    
    // Calculate sum of VMRs of all available foreign broadening species (we need this
    // for normalization). The species "Self" will not be included in the sum!
    Numeric broad_spec_vmr_sum = 0;
    
    // Gamma is the line width. We first initialize gamma with the self width
    gamma = mdata[0][0] * std::pow(theta, mdata[1][0]) * self_pressure;
    
    // Set deltaf to 0
    deltaf = 0;
    
    // and treat foreign width separately:
    Numeric foreign_gamma = 0;
    
    // Add up foreign broadening species, where available:
    for (Index i=0; i<nbs; ++i) {
        if ( broad_spec_locations[i] < -1 ) {
            // -2 means that this broadening species is identical to Self.
            // Throw runtime errors if the parameters are not identical.
            if (mdata[2][i]!=mdata[0][0] ||
                mdata[3][i]!=mdata[1][0])
            {
                std::ostringstream os;
                os << "Inconsistency in LineRecord, self broadening and line "
                << "broadening for " << LineRecord::BroadSpecName(i) << "\n"
                << "should be identical.\n";
                throw std::runtime_error(os.str());
            }
        } else if ( broad_spec_locations[i] >= 0 ) {
            
            // Add to VMR sum:
            broad_spec_vmr_sum += vmrs[broad_spec_locations[i]];
            
            // foreign broadening:
            foreign_gamma +=  mdata[2][i] * std::pow(theta, mdata[3][i])
                * vmrs[broad_spec_locations[i]];
            
            // Delta f (not .25+1.5*foreign_broadening)
            deltaf += mdata[4][i]
            * std::pow( theta , (Numeric).25 + (Numeric)1.5*mdata[3][i] )
            * vmrs[broad_spec_locations[i]];
        }
    }
    
    // Check that sum of self and all foreign VMRs is not too far from 1:
    if ( abs(vmrs[this_species]+broad_spec_vmr_sum-1) > 0.1
        && out2.sufficient_priority() )
    {
        std::ostringstream os;
        os << "Warning: The total VMR of all your defined broadening\n"
        << "species (including \"self\") is "
        << vmrs[this_species]+broad_spec_vmr_sum
        << ", more than 10% " << "different from 1.\n";
        out2 << os.str();
    }
    
    // Normalize foreign gamma and deltaf with the foreign VMR sum (but only if
    // we have any foreign broadening species):
    if (broad_spec_vmr_sum != 0.)
    {
        foreign_gamma /= broad_spec_vmr_sum;
        deltaf        /= broad_spec_vmr_sum;
    }
    else if (self_pressure > 0.)
        // If there are no foreign broadening species present, the best assumption
        // we can make is to use gamma_self in place of foreign_gamma. for deltaf
        // there is no equivalent solution, as we don't have a Delta_self and don't
        // know which other Delta we should apply (in this case delta_f gets 0,
        // which should be okayish):
    {
        foreign_gamma = gamma/self_pressure;
    }
    // It can happen that broad_spec_vmr_sum==0 AND p_self==0 (e.g., when p_grid
    // exceeds the given atmosphere and zero-padding is applied). In this case,
    // both gamma_foreign and deltaf are 0 and we leave it like that.
    
    // Multiply by pressure. For the width we take only the foreign pressure.
    // This is consistent with that we have scaled with the sum of all foreign
    // broadening VMRs. In this way we make sure that the total foreign broadening
    // scales with the total foreign pressure.
    foreign_gamma  *= foreign_pressure;
    
    // For the width, add foreign parts:
    gamma += foreign_gamma;
    
    // For the shift we simply take the total pressure, since there is no self part.
    deltaf *= pressure;
    
    // That's it, we're done.
}


///////////////////////////////////////////
//  Catalog interactions here
///////////////////////////////////////////

// Use these to insert the data in the required format from catalog readings
void PressureBroadeningData::SetAirBroadeningFromCatalog(const Numeric& sgam, 
                                                         const Numeric& nself,
                                                         const Numeric& agam,
                                                         const Numeric& nair,
                                                         const Numeric& air_pressure_DF,
                                                         const Numeric& dsgam, 
                                                         const Numeric& dnself,
                                                         const Numeric& dagam,
                                                         const Numeric& dnair,
                                                         const Numeric& dair_pressure_DF) 
{
    mtype = PB_AIR_BROADENING;
    mdata.resize(5);
    mdataerror.resize(5);
    for(Index ii=0;ii<5;ii++)
    {
        mdata[ii].resize(1);
        mdataerror[ii].resize(1);
    }
    mdata[0][0] = sgam;             // Self broadening gamma parameter
    mdata[1][0] = nself;            // Self broadening n parameter
    mdata[2][0] = agam;             // Air broadening gamma parameter
    mdata[3][0] = nair;             // Air broadening n parameter
    mdata[4][0] = air_pressure_DF;  // Pressure shift parameter
    
    mdataerror[0][0] = dsgam;             // Self broadening gamma parameter
    mdataerror[1][0] = dnself;            // Self broadening n parameter
    mdataerror[2][0] = dagam;             // Air broadening gamma parameter
    mdataerror[3][0] = dnair;             // Air broadening n parameter
    mdataerror[4][0] = dair_pressure_DF;  // Pressure shift parameter
}


// Use these to insert the data in the required format from catalog readings
void PressureBroadeningData::SetAirAndWaterBroadeningFromCatalog(const Numeric& sgam, 
                                                                 const Numeric& sn, 
                                                                 const Numeric& sdelta, 
                                                                 const Numeric& agam,
                                                                 const Numeric& an,
                                                                 const Numeric& adelta,
                                                                 const Numeric& wgam,
                                                                 const Numeric& wn,
                                                                 const Numeric& wdelta) 
{
    mtype = PB_AIR_AND_WATER_BROADENING;
    mdata.resize(3);
    mdataerror.resize(0);
    
    mdata[0].resize(3);
    mdata[0][0] = sgam;             // Self broadening gamma parameter
    mdata[0][1] = sn;               // Self broadening n parameter
    mdata[0][2] = sdelta;           // Self broadening shift parameter
    
    mdata[1].resize(3);
    mdata[1][0] = agam;             // Air broadening gamma parameter
    mdata[1][1] = an;               // Air broadening n parameter
    mdata[1][2] = adelta;           // Air broadening shift parameter
    
    mdata[2].resize(3);
    mdata[2][0] = wgam;             // Water broadening gamma parameter
    mdata[2][1] = wn;               // Water broadening n parameter
    mdata[2][2] = wdelta;           // Water broadening shift parameter
}


// Use these to insert the data in the required format from catalog readings
void PressureBroadeningData::SetPerrinBroadeningFromCatalog(const Numeric& sgam, 
                                                            const Numeric& nself,
                                                            const Vector&  foreign_gamma,
                                                            const Vector&  n_foreign,
                                                            const Vector&  foreign_pressure_DF) 
{
    // All vectors must have the same length
    assert(n_foreign.nelem()==foreign_gamma.nelem());
    assert(foreign_pressure_DF.nelem()==foreign_gamma.nelem());
    
    mtype = PB_PERRIN_BROADENING;
    mdata.resize(5);
    mdataerror.resize(0);
    for(Index ii=0;ii<2;ii++) 
    {
        mdata[ii].resize(1);
    }
    mdata[0][0] = sgam;         // Self broadening gamma parameter
    mdata[1][0] = nself;        // Self broadening n parameter
    mdata[2] = foreign_gamma;   // Gas broadening gamma parameter per species
    mdata[3] = n_foreign;       // Gas broadening n parameter per species
    mdata[4] = foreign_pressure_DF;      // Pressure shift parameter per species
}

///////////////////////////////////////////
// Formating and readings here
///////////////////////////////////////////

Index PressureBroadeningData::ExpectedVectorLengthFromType() const
{
    if(mtype == PB_NONE) // The none case
        return 0;
    else if(mtype == PB_AIR_BROADENING) // 10 Numerics
        return 10;
    else if(mtype == PB_AIR_AND_WATER_BROADENING) // 10 Numerics
        return 9;
    else if(mtype == PB_PERRIN_BROADENING) // 2 Numerics and 3 Vectors of 6-length
        return 20;
    else
        throw std::runtime_error("You are trying to store to a pressure broadening type that is unknown to ARTS.\n");
    return 0;
}

void PressureBroadeningData::SetDataFromVectorWithKnownType(const Vector & input)
{
    if(input.nelem()!=ExpectedVectorLengthFromType())
        throw std::runtime_error("Input pressure broadening is of wrong length.\n");
    
    if(mtype == PB_NONE) 
    {
        // The none case
        mdata.resize(0);
        mdataerror.resize(0);
    } 
    else if(mtype == PB_AIR_BROADENING) // 10 Numerics
    {
        SetAirBroadeningFromCatalog(input[0],
                                    input[1],
                                    input[2],
                                    input[3],
                                    input[4],
                                    input[5],
                                    input[6],
                                    input[7],
                                    input[8],
                                    input[9]);
    }
    else if(mtype == PB_AIR_AND_WATER_BROADENING) // 9 Numerics
    {
        SetAirAndWaterBroadeningFromCatalog(input[0],
                                            input[1],
                                            input[2],
                                            input[3],
                                            input[4],
                                            input[5],
                                            input[6],
                                            input[7],
                                            input[8]);
    }
    else if(mtype == PB_PERRIN_BROADENING) // 2 Numerics and 3 Vectors of 6-length
        SetPerrinBroadeningFromCatalog(input[0],
                                       input[7],
                                       input[Range(1,6)],
                                       input[Range(8,6)],
                                       input[Range(14,6)]);
}

void PressureBroadeningData::StorageTag2SetType(const String & input)
{
    if(input == "NA") // The none case
        mtype=PB_NONE;
    else if(input == "N2") // Air Broadening is N2 broadening mostly...
        mtype=PB_AIR_BROADENING;
    else if(input == "WA") // Water and Air Broadening
        mtype=PB_AIR_AND_WATER_BROADENING;
    else if(input == "AP") // Perrin broadening
        mtype=PB_PERRIN_BROADENING;
    else
        throw std::runtime_error("You are trying to set pressure broadening type that is unknown to ARTS.\n");
}

void PressureBroadeningData::GetVectorFromData(Vector& output) const 
{
    output.resize(ExpectedVectorLengthFromType());
    
    if(mtype == PB_NONE) {} //Nothing
    else if(mtype == PB_AIR_BROADENING) // 10 Numerics
    {
        output[0]=mdata[0][0];
        output[1]=mdata[1][0];
        output[2]=mdata[2][0];
        output[3]=mdata[3][0];
        output[4]=mdata[4][0];
        output[5]=mdataerror[0][0];
        output[6]=mdataerror[1][0];
        output[7]=mdataerror[2][0];
        output[8]=mdataerror[3][0];
        output[9]=mdataerror[4][0];
        
    }
    else if(mtype == PB_AIR_AND_WATER_BROADENING) // 9 Numerics
    {
        output[0]=mdata[0][0];
        output[1]=mdata[0][1];
        output[2]=mdata[0][2];
        
        output[3]=mdata[1][0];
        output[4]=mdata[1][1];
        output[5]=mdata[1][2];
        
        output[6]=mdata[2][0];
        output[7]=mdata[2][1];
        output[8]=mdata[2][2];
    }
    else if(mtype == PB_PERRIN_BROADENING) // 2 Numerics and 3 Vectors of 6-length
    {
        output[0]=mdata[0][0];
        output[7]=mdata[1][0];
        output[Range(1,6)]=mdata[2];
        output[Range(8,6)]=mdata[3];
        output[Range(14,6)]=mdata[4];
    }
}


String PressureBroadeningData::Type2StorageTag() const
{
    String output;

    if(mtype==PB_NONE) // The none case
        output = "NA";
    else if(mtype==PB_AIR_BROADENING) // Air Broadening is mostly N2 ...
        output = "N2";
    else if(mtype==PB_AIR_AND_WATER_BROADENING) // Water and Air Broadening
        output = "WA";
    else if(mtype==PB_PERRIN_BROADENING) // Perring broadening
        output="AP";
    else
        throw std::runtime_error("You are trying to set pressure broadening type that is unknown to ARTS.\n");

    return output;
}
