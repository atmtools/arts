/* Copyright (C) 2014
 R ichard Larsson <ric.larsson@gmail.com>
 
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

/** Contains the pressure broadening data class
 * \file   pressurebroadeningdata.h
 * 
 * \author Richard Larsson
 * \date   2014-11-06
 **/

#ifndef pressurebroadeningdata_h
#define pressurebroadeningdata_h

#include <cmath>
#include "array.h"
#include "messages.h"
#include "matpackI.h"



/**
 This class should contain the data and metadata associated with pressure broadening.
 
 The stored quantities are:
 \param mdata  containing the structured data
 \param mdataerror  containing the structured errors of the mdata input
 \param mtype  containing a tag for the type of data/the structure of the data
 */
class PressureBroadeningData
{
public:
    
    enum PB_Type {
        PB_NONE,                          // No pressure broadening
        PB_AIR_BROADENING,                // Air broadening and self broadening only
        PB_AIR_AND_WATER_BROADENING,      // Air, water, and self broadening
        PB_PERRIN_BROADENING              // Gas broadening as done by A. Perrin
    };
    
    // Defining an object with no data and no broadening
    PressureBroadeningData() : mtype(PB_NONE), mdata(), mdataerror() {}
    
    // Use these to get the raw data from this class
    const PB_Type& Type() const {return mtype;}
    const ArrayOfVector& Data() const {return mdata;}
    const ArrayOfVector& DataError() const {return mdataerror;}
    
    //  The functions below this point sets up pressure broadening schemes from external catalogs and old ARTSCAT
    
    /** 
     The HITRAN air broadening scheme was used exclusively for a long time by ARTS.
     This function follows that scheme and ARTSCAT-3
     
     The quantities are:
     \param sgam  self broadening gamma
     \param nself  self broadening n
     \param agam  air broadening gamma
     \param nair  air broadening n
     \param air_pressure_DF  air broadening pressure-based frequency shift
     \param d*  variable associated errors.  A tag of -1 is used when the errors are unknown.
     */
    void SetAirBroadeningFromCatalog(const Numeric& sgam, 
                                     const Numeric& nself,
                                     const Numeric& agam,
                                     const Numeric& nair,
                                     const Numeric& air_pressure_DF,
                                     const Numeric& dsgam, 
                                     const Numeric& dnself,
                                     const Numeric& dagam,
                                     const Numeric& dnair,
                                     const Numeric& dair_pressure_DF);
    
    /** 
     The HITRAN air broadening scheme was used exclusively for a long time by ARTS.
     This function follows that scheme and ARTSCAT-3, but also adds water broadening.
     Since it is new, we skip mdataerror usage for now. 
     
     The quantities are:
     \param sgam  self broadening gamma
     \param nself  self broadening n
     \param agam  air broadening gamma
     \param nair  air broadening n
     \param air_pressure_DF  air broadening pressure-based frequency shift
     \param d*  variable associated errors.  A tag of -1 is used when the errors are unknown.
     */
    void SetAirAndWaterBroadeningFromCatalog(const Numeric& sgam, 
                                             const Numeric& sn, 
                                             const Numeric& sdelta, 
                                             const Numeric& agam,
                                             const Numeric& an,
                                             const Numeric& adelta,
                                             const Numeric& wgam,
                                             const Numeric& wn,
                                             const Numeric& wdelta);
    
    
    /**
     The ESA planetary study introduced species-dependent pressure broadenings.
     This function can read ARTSCAT-4 line data
     
     The quantities are:
     \param sgam  self broadening gamma
     \param nself  self broadening n
     \param foreign_gamma  foreign broadening gamma
     \param n_foreign  foreign broadening n
     \param foreign_pressure_DF  foreign broadening pressure-based frequency shift
     */ 
    void SetPerrinBroadeningFromCatalog(const Numeric& sgam, 
                                        const Numeric& nself,
                                        const Vector&  foreign_gamma,
                                        const Vector&  n_foreign,
                                        const Vector&  foreign_pressure_DF);
    
    
    // The return functions below are self-explanatory and are all around for backwards compatibility
    ConstVectorView PerrinGammaForeign() const { assert(PB_PERRIN_BROADENING==mtype); return mdata[2]; }
    ConstVectorView PerrinNForeign() const { assert(PB_PERRIN_BROADENING==mtype); return mdata[3]; }
    ConstVectorView PerrinDeltaForeign() const { assert(PB_PERRIN_BROADENING==mtype); return mdata[4]; }
    Numeric PerrinGammaForeign(Index ii) const { assert(PB_PERRIN_BROADENING==mtype); return mdata[2][ii]; }
    Numeric PerrinNForeign(Index ii) const { assert(PB_PERRIN_BROADENING==mtype); return mdata[3][ii]; }
    Numeric PerrinDeltaForeign(Index ii) const {  assert(PB_PERRIN_BROADENING==mtype); return mdata[4][ii]; }
    Numeric AirBroadeningAgam() const { assert(PB_AIR_BROADENING==mtype);  return mdata[2][0]; }
    Numeric AirBroadeningNair() const { assert(PB_AIR_BROADENING==mtype);  return mdata[3][0]; }
    Numeric AirBroadeningPsf()  const { assert(PB_AIR_BROADENING==mtype);  return mdata[4][0]; }
    Numeric Sgam() const { assert(mdata.nelem() != 0); assert(mdata[0].nelem() != 0); return mdata[0][0]; }
    Numeric Nself() const { assert(mdata.nelem() > 0); assert(mdata[1].nelem() != 0); return mdata[1][0]; }
    Numeric AirBroadeningDAgam() const { assert(PB_AIR_BROADENING==mtype);  return mdataerror[2][0]; }
    Numeric AirBroadeningDNair() const { assert(PB_AIR_BROADENING==mtype);  return mdataerror[3][0]; }
    Numeric AirBroadeningDPsf()  const { assert(PB_AIR_BROADENING==mtype);  return mdataerror[4][0]; }
    Numeric dSgam()  const { assert(PB_AIR_BROADENING==mtype);  return mdataerror[0][0]; }
    Numeric dNself() const { assert(PB_AIR_BROADENING==mtype);  return mdataerror[1][0]; }
    
    
    // Use these to return data in the format required by the line shape calculator
    
    /**
     Air broadening calculations
     
     The quantities are:
     \param gamma  the pressure broadening in Hz
     \param deltaf  the pressure shift in Hz
     \param theta  the scaled temperature (T0/T)
     \param pressure  All gasses, in Pa
     \param self_pressure  pressure of the molecule the line belongs to
     */
    void GetAirBroadening(Numeric& gamma,
                          Numeric& deltaf,
                          const Numeric& theta,
                          const Numeric& pressure,
                          const Numeric& self_pressure) const;

    /**
    Air and water broadening calculations                              *
    
    The quantities are:
    \param gamma  the pressure broadening in Hz
    \param deltaf  the pressure shift in Hz
    \param theta  the scaled temperature (T0/T)
    \param pressure  All gasses, in Pa
    \param self_pressure  pressure of the molecule the line belongs to
    */
    void GetAirAndWaterBroadening(Numeric& gamma,
                                  Numeric& deltaf,
                                  const Numeric& theta,
                                  const Numeric& pressure,
                                  const Numeric& self_pressure,
                                  const Index    this_species,
                                  const Index    h2o_species,
                                  ConstVectorView vmrs,
                                  const Verbosity& verbosity) const;
    
    /**
     Perrin broadening calculations
     
     The quantities are:
     \param gamma  the pressure broadening in Hz
     \param deltaf  the pressure shift in Hz
     \param theta  the scaled temperature (T0/T)
     \param pressure  All gasses, in Pa
     \param self_pressure  pressure of the molecule the line belongs to
     \param this_species  location of this species in broad_spec_locations and vmrs
     \param broad_spec_locations  location of species that broadens the line in vmrs
     \param vmrs  volume mixing ratio of all species
     \param verbosity  output command used for warnings.
     */
    void GetPerrinBroadening(Numeric& gamma,
                             Numeric& deltaf,
                             const Numeric& theta,
                             const Numeric& pressure,
                             const Numeric& self_pressure,
                             const Index    this_species,
                             const ArrayOfIndex& broad_spec_locations,
                             ConstVectorView vmrs,
                             const Verbosity& verbosity) const;
    
     /**
       All broadening calculations
       
       The quantities are:
      \param gamma  the pressure broadening in Hz
      \param deltaf  the pressure shift in Hz
      \param theta  the scaled temperature (T0/T)
      \param pressure  All gasses, in Pa
      \param self_pressure  pressure of the molecule the line belongs to
      \param this_species  location of this species in broad_spec_locations and vmrs
      \param broad_spec_locations  location of species that broadens the line in vmrs
      \param vmrs  volume mixing ratio of all species
      \param verbosity  output command used for warnings.
      */
     void GetPressureBroadeningParams(Numeric& gamma_0,
                                      Numeric& gamma_2,
                                      Numeric& eta,
                                      Numeric& df_0,
                                      Numeric& df_2,
                                      Numeric& f_VC,
                                      const Numeric& theta,
                                      const Numeric& pressure,
                                      const Numeric& self_pressure,
                                      const Index    this_species,
                                      const Index    h2o_species,
                                      const ArrayOfIndex& broad_spec_locations,
                                      ConstVectorView vmrs,
                                      const Verbosity& verbosity) const;

    
    
   // Use these to read data from XML-formats and to return vectors for writing
   void StorageTag2SetType(const String& input);
   Index ExpectedVectorLengthFromType() const;
   void SetDataFromVectorWithKnownType(const Vector& input);
   void GetVectorFromData(Vector& output) const;
   String Type2StorageTag() const;
   
private:
    // mtype identifies the type of of pressure broadening and the other variables
    // are containers
    PB_Type       mtype;
    ArrayOfVector mdata;
    ArrayOfVector mdataerror;
};

#endif //pressurebroadeningdata_h