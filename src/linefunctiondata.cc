/* Copyright (C) 2018
 Richard Larsson <larsson@mps.mpg.de>
 
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
 * \file   linefunctiondata.cc
 * 
 * \author Richard Larsson
 * \date   2018-09-19
 **/


#include "linefunctiondata.h"
#include "file.h"


/*
 Main functions of this file is to compute a 
 variable based on some temperature fit.  The
 implemented temperature fits are:
  
  t0: a constant
  t1: standard HITRAN, x0 (T0 / T) ^ x1
  t2: line shifts, x0 (T0 / T) ^ x1 / (1 + x2 ln(T / T0))
  t3: speed-dependent parameters, x0 + x1 (T - T0)
  t4: second order line mixing, (x0 + x1 (T0 / T - 1)) (T0 / T) ^ x2
  t5: ARTS pressure shift, x0 * (T0 / T) ^ (0.25 + 1.5*x1)
  
 Note that each of these temperature fits need
 to have not only a main function implemented,
 but a derivative based on how many parameters
 are required for ALL OTHER temperature fits.
 
 At the moment, we have three coefficients, the 
 atmospheric temperature, and the temperature
 at which the line parameters were derived.
 This means that each main function must implement
 derivative functionality with regards to:
 
  x0
  x1,
  x2,
  T, and
  T0,
  
 where x0-x2 are coefficients often derived from lab
 work, T is the atmospheric temperature, and T0 is the
 temperature at which x0-x2 have been derived.  Note
 that the latter has to be the same not just for line
 broadening and mixing parameters, but also for other
 line parameters, notably, the line strength.
 
*/


/**********************************
 Main functionality, TH = T0/T 
**********************************/

//! Returns x0
inline Numeric main_t0(const Numeric& x0) noexcept { return x0; }

//! Returns x0 * pow(TH, x1)
inline Numeric main_t1(const Numeric& TH, const Numeric& x0, const Numeric& x1) noexcept  { return x0 * pow(TH, x1); }

//! Returns x0 * pow(TH, x1) * (1 + x2 * log(1/TH))
inline Numeric main_t2(const Numeric& TH, const Numeric& x0, const Numeric& x1, const Numeric& x2) noexcept { return x0 * pow(TH, x1) * (1 + x2 * log(1/TH)); }

//! Returns x0 + x1 * (T - T0)
inline Numeric main_t3(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1) noexcept { return x0 + x1 * (T - T0); }

//! Returns (x0 + x1 * (TH - 1)) * pow(TH, x2)
inline Numeric main_t4(const Numeric& TH, const Numeric& x0, const Numeric& x1, const Numeric& x2) noexcept { return (x0 + x1 * (TH - 1)) * pow(TH, x2); }

//! Returns x0 * pow(TH, 0.25 + 1.5*x1)
inline Numeric main_t5(const Numeric& TH, const Numeric& x0, const Numeric& x1) noexcept { return x0 * pow(TH, 0.25 + 1.5*x1); }


/**********************************
 Derivatives with regards to x0  
**********************************/

//! Returns 1
inline Numeric dmain_dx0_t0() noexcept { return 1; }

//! Returns pow(T0/T, x1)
inline Numeric dmain_dx0_t1(const Numeric& T, const Numeric& T0, const Numeric& x1) noexcept { return pow(T0/T, x1); }

//! Returns pow(T0/T, x1)*(x2*log(T/T0) + 1)
inline Numeric dmain_dx0_t2(const Numeric& T, const Numeric& T0, const Numeric& x1, const Numeric& x2) noexcept { return pow(T0/T, x1)*(x2*log(T/T0) + 1); }

//! Returns 1
inline Numeric dmain_dx0_t3() noexcept { return 1; }

//! Returns pow(T0/T, x2)
inline Numeric dmain_dx0_t4(const Numeric& T, const Numeric& T0, const Numeric& x2) noexcept { return pow(T0/T, x2); }

//! Returns pow(T0/T, 1.5*x1 + 0.25)
inline Numeric dmain_dx0_t5(const Numeric& T, const Numeric& T0, const Numeric& x1) noexcept { return pow(T0/T, 1.5*x1 + 0.25); }


/**********************************
 Derivatives with regards to x1  
**********************************/

//! Returns 0
inline Numeric dmain_dx1_t0() noexcept { return 0; }

//! Returns x0*pow(T0/T, x1)*log(T0/T)
inline Numeric dmain_dx1_t1(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1) noexcept { return x0*pow(T0/T, x1)*log(T0/T); }

//! Returns x0*pow(T0/T, x1)*(x2*log(T/T0) + 1)*log(T0/T)
inline Numeric dmain_dx1_t2(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1, const Numeric& x2) noexcept { return x0*pow(T0/T,x1)*(x2*log(T/T0)+1)*log(T0/T); }

//! Returns (T - T0)
inline Numeric dmain_dx1_t3(const Numeric& T, const Numeric& T0) noexcept { return (T - T0); }

//! Returns pow(T0/T, x2)*(T0/T - 1)
inline Numeric dmain_dx1_t4(const Numeric& T, const Numeric& T0, const Numeric& x2) noexcept { return pow(T0/T, x2)*(T0/T - 1); }

//! Returns 1.5*x0*pow(T0/T, 1.5*x1 + 0.25)*log(T0/T)
inline Numeric dmain_dx1_t5(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1) noexcept { return 1.5*x0*pow(T0/T, 1.5*x1 + 0.25)*log(T0/T); }


/**********************************
 Derivatives with regards to x2  
**********************************/

//! Returns 0
inline Numeric dmain_dx2_t0() noexcept { return 0; }

//! Returns 0
inline Numeric dmain_dx2_t1() noexcept { return 0; }

//! Returns x0*pow(T0/T, x1)*log(T/T0)
inline Numeric dmain_dx2_t2(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1) noexcept { return x0*pow(T0/T, x1)*log(T/T0); }

//! Returns 0
inline Numeric dmain_dx2_t3() noexcept { return 0; }

//! Returns pow(T0/T, x2)*(x0 + x1*(T0/T-1))*log(T0/T)
inline Numeric dmain_dx2_t4(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1, const Numeric& x2) noexcept { return pow(T0/T,x2)*(x0+x1*(T0/T-1))*log(T0/T); }

//! Returns 0
inline Numeric dmain_dx2_t5() noexcept { return 0; }


/**********************************
 Derivatives with regards to T 
**********************************/

//! Returns 0
inline Numeric dmain_dT_t0() noexcept { return 0; }

//! Returns -x0*x1*pow(T0/T, x1)/T
inline Numeric dmain_dT_t1(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1) noexcept { return -x0*x1*pow(T0/T, x1)/T; }

//! Returns -x0*x1*pow(T0/T, x1)*(x2*log(T/T0) + 1)/T + x0*x2*pow(T0/T, x1)/T
inline Numeric dmain_dT_t2(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1, const Numeric& x2) noexcept { return -x0*x1*pow(T0/T,x1)*(x2*log(T/T0)+1)/T+x0*x2*pow(T0/T,x1)/T; }

//! Returns x1
inline Numeric dmain_dT_t3(const Numeric& x1) noexcept { return x1; }

//! Returns -x2*pow(T0/T, x2)*(x0 + x1*(T0/T-1))/T - T0*x1*pow(T0/T, x2)/pow(T, 2)
inline Numeric dmain_dT_t4(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1, const Numeric& x2) noexcept { return -x2*pow(T0/T,x2)*(x0+x1*(T0/T-1))/T-T0*x1*pow(T0/T,x2)/pow(T,2); }

//! Returns -x0*pow(T0/T, 1.5*x1 + 0.25)*(1.5*x1 + 0.25)/T
inline Numeric dmain_dT_t5(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1) noexcept { return -x0*pow(T0/T,1.5*x1 + 0.25)*(1.5*x1 + 0.25)/T; }


/**********************************
 Derivatives with regards to T0  
**********************************/

//! Returns 0
inline Numeric dmain_dT0_t0() noexcept { return 0; }

//! Returns x0*x1*pow(T0/T, x1)/T0
inline Numeric dmain_dT0_t1(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1) noexcept { return x0*x1*pow(T0/T, x1)/T0; }

//! Returns x0*x1*pow(T0/T, x1)*(x2*log(T/T0) + 1)/T0 - x0*x2*pow(T0/T, x1)/T0
inline Numeric dmain_dT0_t2(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1, const Numeric& x2) noexcept { return x0*x1*pow(T0/T,x1)*(x2*log(T/T0)+1)/T0-x0*x2*pow(T0/T,x1)/T0; }

//! Returns -x1
inline Numeric dmain_dT0_t3(const Numeric& x1) noexcept { return -x1; }

//! Returns x2*pow(T0/T, x2)*(x0 + x1*(T0/T - 1))/T0 + x1*pow(T0/T, x2)/T
inline Numeric dmain_dT0_t4(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1, const Numeric& x2) noexcept { return x2*pow(T0/T,x2)*(x0+x1*(T0/T-1))/T0+x1*pow(T0/T,x2)/T; }

//! Returns x0*pow(T0/T, 1.5*x1 + 0.25)*(1.5*x1 + 0.25)/T0
inline Numeric dmain_dT0_t5(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1) noexcept { return x0*pow(T0/T,1.5*x1 + 0.25)*(1.5*x1 + 0.25)/T0; }


//! Select line parameter based on parameter order
/*! 
 * This helper code simply selects which of the input the
 * param-index is pointing at as a function of the type of
 * line shape under consideration
 * 
 * \param G0    The speed-independent pressure broadening term
 * \param D0    The speed-independent pressure shift term
 * \param G2    The speed-dependent pressure broadening term
 * \param D2    The speed-dependent pressure shift term
 * \param FVC   The velocity-changing collisional parameter
 * \param ETA   The correlation between velocity and rotational changes due to collisions
 * \param param Index enumerating the parameter
 * \param type  Enum class for knowing which parameter we are interested in
 * 
 * \return A reference pointing at G0, D0, G2, D2, FVC, or ETA.
 * 
 * \author Richard Larsson
 * \date   2018-07-16
 */
inline Numeric& select_line_shape_param(Numeric& G0,
                                        Numeric& D0,
                                        Numeric& G2,
                                        Numeric& D2,
                                        Numeric& FVC,
                                        Numeric& ETA,
                                        const Index param,
                                        const LineFunctionData::LineShapeType type) noexcept
{
  switch(type) {
    case LineFunctionData::LineShapeType::DP:
      // Do not let this happen!
    case LineFunctionData::LineShapeType::LP:
      switch(param) {
        case Index(LineFunctionData::LorentzParam::G0): return G0;
        case Index(LineFunctionData::LorentzParam::D0): return D0;
      }
    case LineFunctionData::LineShapeType::VP:
      switch(param) {
        case Index(LineFunctionData::VoigtParam::G0): return G0;
        case Index(LineFunctionData::VoigtParam::D0): return D0;
      }
    case LineFunctionData::LineShapeType::SDVP:
      switch(param) {
        case Index(LineFunctionData::SpeedVoigtParam::G0): return G0;
        case Index(LineFunctionData::SpeedVoigtParam::D0): return D0;
        case Index(LineFunctionData::SpeedVoigtParam::G2): return G2;
        case Index(LineFunctionData::SpeedVoigtParam::D2): return D2;
      }
    case LineFunctionData::LineShapeType::HTP:
      switch(param) {
        case Index(LineFunctionData::HTPParam::G0):  return G0;
        case Index(LineFunctionData::HTPParam::D0):  return D0;
        case Index(LineFunctionData::HTPParam::G2):  return G2;
        case Index(LineFunctionData::HTPParam::D2):  return D2;
        case Index(LineFunctionData::HTPParam::FVC): return FVC;
        case Index(LineFunctionData::HTPParam::ETA): return ETA;
      }
  }
  return G0;  // Only to suppress warnings!
}


//! Select line parameter based on parameter order
/*! 
 * This helper code simply selects which of the input the
 * param-index is pointing at as a function of the type of
 * line shape under consideration
 * 
 * \param Y     The imaginary part of the line shape due to line mixing
 * \param G     The strength-altering due to line mixing
 * \param DV    The frequency shift due to line mixing
 * \param param Index enumerating the parameter
 * \param type  Enum class for knowing which parameter we are interested in
 * 
 * \return A reference pointing at Y, G, or DV.
 * 
 * \author Richard Larsson
 * \date   2018-07-16
 */
inline Numeric& select_line_mixing_param(Numeric& Y,
                                         Numeric& G,
                                         Numeric& DV,
                                         const Index param,
                                         const LineFunctionData::LineMixingOrderType type) noexcept
{
  switch(type) {
    case LineFunctionData::LineMixingOrderType::None:
      // Do not let this happen
    case LineFunctionData::LineMixingOrderType::Interp:
      // This does not matter;
    case LineFunctionData::LineMixingOrderType::LM1:
      switch(param) {
        case Index(LineFunctionData::FirstOrderParam::Y): return Y;
      }
    case LineFunctionData::LineMixingOrderType::LM2:
      switch(param) {
        case Index(LineFunctionData::SecondOrderParam::Y):  return Y;
        case Index(LineFunctionData::SecondOrderParam::G):  return G;
        case Index(LineFunctionData::SecondOrderParam::DV): return DV;
      }
  }
  return Y;  // Only to suppress warnings!
}


//! Special function for line mixing of LBLRTM type
/*! 
 * LBLRTM interpolates linearly a set of variables
 * 
 * \param Y           The imaginary part of the line shape due to line mixing
 * \param G           The strength-altering due to line mixing
 * \param partial_vmr The VMR of the species in question
 * \param T           The atmospheric temperature
 * \param data        The associated data struction [T1, T2, T3, T4, Y1, Y2, Y3, Y4, G1, G2, G3, G4]
 * 
 * \author Richard Larsson
 * \date   2018-07-16
 */
inline void special_line_mixing_aer(Numeric& Y,
                                    Numeric& G,
                                    const Numeric& partial_vmr,
                                    const Numeric& T,
                                    ConstVectorView data) {
  const Index n=data.nelem();
  if(n not_eq 12)
    throw std::runtime_error("Data for line mixing not matching AER LM style despite marked as such...\n must have this structure: [T1, T2, T3, T4, Y1, Y2, Y3, Y4, G1, G2, G3, G4]");
  
  // Perform linear interpolation
  if(T < data[1]) {
    Y += partial_vmr * (data[4] + (T - data[0]) * (data[5] - data[4]) / (data[1] - data[0]));
    G += partial_vmr * (data[8] + (T - data[0]) * (data[9] - data[8]) / (data[1] - data[0]));
  }
  else if(T > data[2]) {
    Y += partial_vmr * (data[6 ] + (T - data[2]) * (data[7 ] - data[6 ]) / (data[3] - data[2]));
    G += partial_vmr * (data[10] + (T - data[2]) * (data[11] - data[10]) / (data[3] - data[2]));
  }
  else {
    Y = partial_vmr * (data[5] + (T - data[1]) * (data[6 ] - data[5]) / (data[2] - data[1]));
    G = partial_vmr * (data[9] + (T - data[1]) * (data[10] - data[9]) / (data[2] - data[1]));
  }
}


//! Compute the pressure broadening and line mixing parameters
/*! 
 * This helper code simply selects which of the input the
 * param-index is pointing at as a function of the type of
 * line shape under consideration
 * 
 * \param G0          The speed-independent pressure broadening term
 * \param D0          The speed-independent pressure shift term
 * \param G2          The speed-dependent pressure broadening term
 * \param D2          The speed-dependent pressure shift term
 * \param FVC         The velocity-changing collisional parameter
 * \param ETA         The correlation between velocity and rotational changes due to collisions
 * \param Y           The imaginary part of the line shape due to line mixing
 * \param G           The strength-altering due to line mixing
 * \param DV          The frequency shift due to line mixing
 * \param T0          The line reference temperature
 * \param T           The atmospheric temperature
 * \param P           The atmospheric pressure
 * \param self_vmr    The VMR of the species pressure
 * \param rtp_vmr     The VMR of all species in the atmosphere at a specific radiative transfer point
 * \param abs_species The list of all species in the atmosphere
 * 
 * \return A reference pointing at Y, G, or DV.
 * 
 * \author Richard Larsson
 * \date   2018-07-16
 */
void LineFunctionData::GetParams(Numeric& G0, Numeric& D0, 
                                 Numeric& G2, Numeric& D2, 
                                 Numeric& FVC, Numeric& ETA,
                                 Numeric& Y, Numeric& G, Numeric& DV,
                                 const Numeric& T0, const Numeric& T,
                                 const Numeric& P, const Numeric& self_vmr,
                                 const ConstVectorView rtp_vmr, 
                                 const ArrayOfSpeciesTag& abs_species,
                                 const bool normalization) const
{
  // Set to zero
  G0 = D0 = G2 = D2 = FVC = ETA = Y = G = DV = 0;
  
  // Doppler broadening has no types...
  if(mp == LineShapeType::DP) return;
  
  // Theta is useful in many models so it gets pre-computed
  const Numeric theta = T0/T;
  
  // Set up holders for partial and total VMR
  Numeric total_vmr=0, partial_vmr;
  
  // Add all species up
  for(Index i=0; i<mspecies.nelem(); i++) {
    if(i == 0 and mself)  // The first value might be self-broadening (use self_vmr)
      partial_vmr = self_vmr;
    else if(i == mspecies.nelem()-1 and mbath)  // The last value might be air-broadening (set total_vmr to 1)
      partial_vmr = 1 - total_vmr;
    else {  // Otherwise we have to find the species in the list of species tags
      Index this_species=-1;
      for(Index j=0; j<abs_species.nelem(); j++)
        if(abs_species[j].Species() == mspecies[i].Species())
          this_species = j;
        
      // If we cannot find the species, it does not exist and we are in a different atmosphere than covered by the line data
      if(this_species == -1) continue;  // Species does not exist
      partial_vmr = rtp_vmr[this_species];
    }
    // Sum up VMR
    total_vmr += partial_vmr;
    
    // Set up a current counter to keep track of the data
    Index current = 0;
    
    // Make reading easier
    #define x0 mdata[i][current]
    #define x1 mdata[i][current+1]
    #define x2 mdata[i][current+2]
    
    // Do the line shape parameters
    for(Index j=0; j<LineShapeTypeNelem(); j++) {
      // Selects one of G0, D0, G2, D2, FVC, ETA based on the order of the data for this line shape type
      Numeric& param = select_line_shape_param(G0, D0, G2, D2, FVC, ETA, j, mp);
      
      switch(mtypes[i][j]) {
        case TemperatureType::None:
          break;
        case TemperatureType::T0:
          param += partial_vmr * main_t0(x0);
          break;
        case TemperatureType::T1:
          param += partial_vmr * main_t1(theta, x0, x1);
          break;
        case TemperatureType::T2:
          param += partial_vmr * main_t2(theta, x0, x1, x2);
          break;
        case TemperatureType::T3:
          param += partial_vmr * main_t3(T, T0, x0, x1);
          break;
        case TemperatureType::T4:
          param += partial_vmr * main_t4(theta, x0, x1, x2);
          break;
        case TemperatureType::T5:
          param += partial_vmr * main_t5(theta, x0, x1);
          break;
        case TemperatureType::LM_AER:
          /* throw std::runtime_error("Not allowed for line shape parameters"); */
          break;
      }
      current += TemperatureTypeNelem(mtypes[i][j]);
    }
    
    // Do the line mixing parameters
    for(Index j=0; j<LineMixingTypeNelem(); j++) {
      Numeric& param = select_line_mixing_param(Y, G, DV, j, mlm);
      switch(mtypes[i][j]) {
        case TemperatureType::None:
          break;
        case TemperatureType::T0:
          param += partial_vmr * main_t0(x0);
          break;
        case TemperatureType::T1:
          param += partial_vmr * main_t1(theta, x0, x1);
          break;
        case TemperatureType::T2:
          param += partial_vmr * main_t2(theta, x0, x1, x2);
          break;
        case TemperatureType::T3:
          param += partial_vmr * main_t3(T, T0, x0, x1);
          break;
        case TemperatureType::T4:
          param += partial_vmr * main_t4(theta, x0, x1, x2);
          break;
        case TemperatureType::T5:
          param += partial_vmr * main_t5(theta, x0, x1);
          break;
        case TemperatureType::LM_AER:
          // Special case adding to both G and Y and nothing else!
          special_line_mixing_aer(Y, G, partial_vmr, T, mdata[i][Range(current, joker)]);
          break;
      }
      current += TemperatureTypeNelem(mtypes[i][j+LineShapeTypeNelem()]);
    }
    
    // Stop destroying names
    #undef x0
    #undef x1
    #undef x2
  }
  
  // If there was no species at all, then we cannot compute the pressure broadening
  if(total_vmr == 0) return;
  
  // Rescale to pressure with or without normalization
  const Numeric scale = normalization ? P / total_vmr : P;
  G0  *= scale;
  D0  *= scale;
  G2  *= scale;
  D2  *= scale;
  FVC *= scale;
  if(normalization) ETA /= total_vmr;
  Y   *= scale;
  G   *= scale * P;
  DV  *= scale * P;
}


//! Compute the pressure broadening and line mixing parameters temperature derivatives
/*! 
 * This helper code simply selects which of the input the
 * param-index is pointing at as a function of the type of
 * line shape under consideration
 * 
 * \param dG0          The speed-independent pressure broadening term
 * \param dD0          The speed-independent pressure shift term
 * \param dG2          The speed-dependent pressure broadening term
 * \param dD2          The speed-dependent pressure shift term
 * \param dFVC         The velocity-changing collisional parameter
 * \param dETA         The correlation between velocity and rotational changes due to collisions
 * \param dY           The imaginary part of the line shape due to line mixing
 * \param dG           The strength-altering due to line mixing
 * \param dDV          The frequency shift due to line mixing
 * \param dT0          The line reference temperature
 * \param T           The atmospheric temperature
 * \param P           The atmospheric pressure
 * \param self_vmr    The VMR of the species pressure
 * \param rtp_vmr     The VMR of all species in the atmosphere at a specific radiative transfer point
 * \param abs_species The list of all species in the atmosphere
 * 
 * \return A reference pointing at Y, G, or DV.
 * 
 * \author Richard Larsson
 * \date   2018-07-16
 */
void LineFunctionData::GetTemperatureDerivs(Numeric& dG0, Numeric& dD0, 
                                            Numeric& dG2, Numeric& dD2, 
                                            Numeric& dFVC, Numeric& dETA,
                                            Numeric& dY, Numeric& dG, Numeric& dDV,
                                            const Numeric& T0, const Numeric& T, const Numeric& dT,
                                            const Numeric& P, const Numeric& self_vmr,
                                            const ConstVectorView rtp_vmr, 
                                            const ArrayOfSpeciesTag& abs_species,
                                            const bool normalization) const
{
  // Set to zero
  dG0 = dD0 = dG2 = dD2 = dFVC = dETA = dY = dG = dDV = 0;
  
  // Doppler broadening has no types...
  if(mp == LineShapeType::DP) return;
  
  // Set up holders for partial and total VMR
  Numeric total_vmr=0, partial_vmr;
  
  // Add all species up
  for(Index i=0; i<mspecies.nelem(); i++) {
    if(i == 0 and mself)  // The first value might be self-broadening (use self_vmr)
      partial_vmr = self_vmr;
    else if(i == mspecies.nelem()-1 and mbath)  // The last value might be air-broadening (set total_vmr to 1)
      partial_vmr = 1 - total_vmr;
    else {  // Otherwise we have to find the species in the list of species tags
      Index this_species=-1;
      for(Index j=0; j<abs_species.nelem(); j++)
        if(abs_species[j].Species() == mspecies[i].Species())
          this_species = j;
        
      // If we cannot find the species, it does not exist and we are in a different atmosphere than covered by the line data
      if(this_species == -1) continue;  // Species does not exist
      partial_vmr = rtp_vmr[this_species];
    }
    // Sum up VMR
    total_vmr += partial_vmr;
    
    // Set up a current counter to keep track of the data
    Index current = 0;
    
    // Make reading easier
    #define x0 mdata[i][current]
    #define x1 mdata[i][current+1]
    #define x2 mdata[i][current+2]
    
    // Do the line shape parameters
    for(Index j=0; j<LineShapeTypeNelem(); j++) {
      // Selects one of G0, D0, G2, D2, FVC, ETA based on the order of the data for this line shape type
      Numeric& param = select_line_shape_param(dG0, dD0, dG2, dD2, dFVC, dETA, j, mp);
      
      switch(mtypes[i][j]) {
        case TemperatureType::None:
          break;
        case TemperatureType::T0:
          param += partial_vmr * dmain_dT_t0();
          break;
        case TemperatureType::T1:
          param += partial_vmr * dmain_dT_t1(T, T0, x0, x1);
          break;
        case TemperatureType::T2:
          param += partial_vmr * dmain_dT_t2(T, T0, x0, x1, x2);
          break;
        case TemperatureType::T3:
          param += partial_vmr * dmain_dT_t3(x1);
          break;
        case TemperatureType::T4:
          param += partial_vmr * dmain_dT_t4(T, T0, x0, x1, x2);
          break;
        case TemperatureType::T5:
          param += partial_vmr * dmain_dT_t5(T, T0, x0, x1);
          break;
        case TemperatureType::LM_AER:
          /* throw std::runtime_error("Not allowed for line shape parameters");*/
          break;
      }
      current += TemperatureTypeNelem(mtypes[i][j]);
    }
    
    // Do the line mixing parameters
    for(Index j=0; j<LineMixingTypeNelem(); j++) {
      Numeric& param = select_line_mixing_param(dY, dG, dDV, j, mlm);
      switch(mtypes[i][j]) {
        case TemperatureType::None:
          break;
        case TemperatureType::T0:
          param += partial_vmr * dmain_dT_t0();
          break;
        case TemperatureType::T1:
          param += partial_vmr * dmain_dT_t1(T, T0, x0, x1);
          break;
        case TemperatureType::T2:
          param += partial_vmr * dmain_dT_t2(T, T0, x0, x1, x2);
          break;
        case TemperatureType::T3:
          param += partial_vmr * dmain_dT_t3(x1);
          break;
        case TemperatureType::T4:
          param += partial_vmr * dmain_dT_t4(T, T0, x0, x1, x2);
          break;
        case TemperatureType::T5:
          param += partial_vmr * dmain_dT_t5(T, T0, x0, x1);
          break;
        case TemperatureType::LM_AER:
          // Special case adding to both G and Y and nothing else!
          {
              Numeric dYp=0, dGp=0;
              special_line_mixing_aer(dY, dG, partial_vmr, T+dT, mdata[i][Range(current, joker)]);
              special_line_mixing_aer(dYp, dGp, partial_vmr, T, mdata[i][Range(current, joker)]);
              dY -= dYp;
              dY /= dT;
              dG -= dGp;
              dG /= dT;
          }
          break;
      }
      current += TemperatureTypeNelem(mtypes[i][j+LineShapeTypeNelem()]);
    }
    
    // Stop destroying names
    #undef x0
    #undef x1
    #undef x2
  }
  
  // If there was no species at all, then we cannot compute the pressure broadening
  if(total_vmr == 0) return;
  
  // Rescale to pressure with or without normalization
  const Numeric scale = normalization ? P / total_vmr : P;
  dG0  *= scale;
  dD0  *= scale;
  dG2  *= scale;
  dD2  *= scale;
  dFVC *= scale;
  if(normalization) dETA /= total_vmr;
  dY   *= scale;
  dG   *= scale * P;
  dDV  *= scale * P;
}


//! Compute the pressure broadening and line mixing parameters reference temperature derivatives for one broadening species
/*! 
 * This helper code simply selects which of the input the
 * param-index is pointing at as a function of the type of
 * line shape under consideration
 * 
 * \param dG0          The speed-independent pressure broadening term
 * \param dD0          The speed-independent pressure shift term
 * \param dG2          The speed-dependent pressure broadening term
 * \param dD2          The speed-dependent pressure shift term
 * \param dFVC         The velocity-changing collisional parameter
 * \param dETA         The correlation between velocity and rotational changes due to collisions
 * \param dY           The imaginary part of the line shape due to line mixing
 * \param dG           The strength-altering due to line mixing
 * \param dDV          The frequency shift due to line mixing
 * \param dT0          The line reference temperature
 * \param T           The atmospheric temperature
 * \param P           The atmospheric pressure
 * \param self_vmr    The VMR of the species pressure
 * \param rtp_vmr     The VMR of all species in the atmosphere at a specific radiative transfer point
 * \param abs_species The list of all species in the atmosphere
 * 
 * \return A reference pointing at Y, G, or DV.
 * 
 * \author Richard Larsson
 * \date   2018-07-16
 */
void LineFunctionData::GetReferenceT0Derivs(Numeric& dG0, Numeric& dD0, 
                                            Numeric& dG2, Numeric& dD2, 
                                            Numeric& dFVC, Numeric& dETA,
                                            Numeric& dY, Numeric& dG, Numeric& dDV,
                                            const Numeric& T0, const Numeric& T,
                                            const Numeric& P, const Numeric& self_vmr,
                                            const ConstVectorView rtp_vmr, 
                                            const ArrayOfSpeciesTag& abs_species,
                                            const RetrievalQuantity& rt, 
                                            const QuantumIdentifier& line_qi,
                                            const bool normalization) const
{
  // Set to zero
  dG0 = dD0 = dG2 = dD2 = dFVC = dETA = dY = dG = dDV = 0;
  
  // Doppler broadening has no types...
  if(mp == LineShapeType::DP) return;
  
  // Set up holders for partial and total VMR
  Numeric total_vmr=0, partial_vmr=0;
  // Find calculation details
  Numeric this_vmr=0.0;
  Index this_derivative=-1;
  for(Index i=0; i<mspecies.nelem(); i++) {
    if(i == 0 and mself) { // The first value might be self-broadening (use self_vmr)
      if(rt.QuantumIdentity().In(line_qi)) {
        this_vmr = self_vmr;
        this_derivative = i;
      }
      partial_vmr = self_vmr;
    }
    else if(i == mspecies.nelem()-1 and mbath) { // The last value might be air-broadening (set total_vmr to 1)
      if(this_derivative == -1) { // if it still is not set by this point, it is going to
        this_vmr = 1 - total_vmr;
        this_derivative = i;
      }
      partial_vmr = 1 - total_vmr;
    }
    else {  // Otherwise we have to find the species in the list of species tags
      Index this_species=-1;
      for(Index j=0; j<abs_species.nelem(); j++) {
        if(abs_species[j].Species() == mspecies[i].Species())
          this_species = j;
      }
      
      if(rt.QuantumIdentity().Species() == mspecies[i].Species() and 
         rt.QuantumIdentity().Isotopologue() == mspecies[i].Isotopologue()) {
        this_vmr = rtp_vmr[this_species];
        this_derivative = i;
      }
      
      // If we cannot find the species, it does not exist and we are in a different atmosphere than covered by the line data
      if(this_species == -1) continue;  // Species does not exist
      partial_vmr = rtp_vmr[this_species];
    }
    
    // Sum up VMR
    total_vmr += partial_vmr;
  }
    
  if(this_derivative == -1) return; 
  
  // Set up a current counter to keep track of the data
  Index current = 0;
  
  // Make reading easier
  #define x0 mdata[this_derivative][current]
  #define x1 mdata[this_derivative][current+1]
  #define x2 mdata[this_derivative][current+2]
  
  // Do the line shape parameters
  for(Index j=0; j<LineShapeTypeNelem(); j++) {
    // Selects one of G0, D0, G2, D2, FVC, ETA based on the order of the data for this line shape type
    Numeric& param = select_line_shape_param(dG0, dD0, dG2, dD2, dFVC, dETA, j, mp);
    
    switch(mtypes[this_derivative][j]) {
      case TemperatureType::None:
        break;
      case TemperatureType::T0:
        param += partial_vmr * dmain_dT0_t0();
        break;
      case TemperatureType::T1:
        param += partial_vmr * dmain_dT0_t1(T, T0, x0, x1);
        break;
      case TemperatureType::T2:
        param += partial_vmr * dmain_dT0_t2(T, T0, x0, x1, x2);
        break;
      case TemperatureType::T3:
        param += partial_vmr * dmain_dT0_t3(x1);
        break;
      case TemperatureType::T4:
        param += partial_vmr * dmain_dT0_t4(T, T0, x0, x1, x2);
        break;
      case TemperatureType::T5:
        param += partial_vmr * dmain_dT0_t5(T, T0,x0, x1);
        break;
      case TemperatureType::LM_AER:
        /* throw std::runtime_error("Not allowed for line shape parameters"); */
        break;
    }
    current += TemperatureTypeNelem(mtypes[this_derivative][j]);
  }
  
  // Do the line mixing parameters
  for(Index j=0; j<LineMixingTypeNelem(); j++) {
    Numeric& param = select_line_mixing_param(dY, dG, dDV, j, mlm);
    switch(mtypes[this_derivative][j]) {
      case TemperatureType::None:
        break;
      case TemperatureType::T0:
        param += partial_vmr * dmain_dT0_t0();
        break;
      case TemperatureType::T1:
        param += partial_vmr * dmain_dT0_t1(T, T0, x0, x1);
        break;
      case TemperatureType::T2:
        param += partial_vmr * dmain_dT0_t2(T, T0, x0, x1, x2);
        break;
      case TemperatureType::T3:
        param += partial_vmr * dmain_dT0_t3(x1);
        break;
      case TemperatureType::T4:
        param += partial_vmr * dmain_dT0_t4(T, T0, x0, x1, x2);
        break;
      case TemperatureType::T5:
        param += partial_vmr * dmain_dT0_t5(T, T0, x0, x1);
        break;
      case TemperatureType::LM_AER:
        // No derivatives for it depends on T not T0
        break;
    }
    current += TemperatureTypeNelem(mtypes[this_derivative][j+LineShapeTypeNelem()]);
  }
  
  // Stop destroying names
  #undef x0
  #undef x1
  #undef x2
  
  // If there was no species at all, then we cannot compute the pressure broadening
  if(total_vmr == 0) return;
  
  // Rescale to pressure with or without normalization
  const Numeric scale = normalization ? this_vmr * P / total_vmr : this_vmr * P;
  dG0  *= scale;
  dD0  *= scale;
  dG2  *= scale;
  dD2  *= scale;
  dFVC *= scale;
  if(normalization) dETA /= total_vmr;
  dY   *= scale;
  dG   *= scale * this_vmr * P;
  dDV  *= scale * this_vmr * P;
}


bool LineFunctionData::ComputesParam(const String& type) const
{
  if(type == "G0")
    return mp == LineShapeType::LP or mp == LineShapeType::VP or
           mp == LineShapeType::SDVP or mp == LineShapeType::HTP;
  else if(type == "D0")
    return mp == LineShapeType::LP or mp == LineShapeType::VP or
           mp == LineShapeType::SDVP or mp == LineShapeType::HTP;
  else if(type == "G2")
    return mp == LineShapeType::SDVP or mp == LineShapeType::HTP;
  else if(type == "D2")
    return mp == LineShapeType::SDVP or mp == LineShapeType::HTP;
  else if(type == "FVC")
    return mp == LineShapeType::HTP;
  else if(type == "ETA")
    return mp == LineShapeType::HTP;
  else if(type == "Y")
    return mlm == LineMixingOrderType::LM1 or mlm == LineMixingOrderType::LM2 or
           mlm == LineMixingOrderType::Interp;
  else if(type == "G")
    return mlm == LineMixingOrderType::LM2 or mlm == LineMixingOrderType::Interp;
  else if(type == "DV")
    return mlm == LineMixingOrderType::LM2;
  return false;
}



//! Compute the pressure broadening and line mixing parameters temperature derivatives
/*! 
 * This helper code simply selects which of the input the
 * param-index is pointing at as a function of the type of
 * line shape under consideration
 * 
 * \param dG0          The speed-independent pressure broadening term
 * \param dD0          The speed-independent pressure shift term
 * \param dG2          The speed-dependent pressure broadening term
 * \param dD2          The speed-dependent pressure shift term
 * \param dFVC         The velocity-changing collisional parameter
 * \param dETA         The correlation between velocity and rotational changes due to collisions
 * \param dY           The imaginary part of the line shape due to line mixing
 * \param dG           The strength-altering due to line mixing
 * \param dDV          The frequency shift due to line mixing
 * \param dT0          The line reference temperature
 * \param T           The atmospheric temperature
 * \param P           The atmospheric pressure
 * \param self_vmr    The VMR of the species pressure
 * \param rtp_vmr     The VMR of all species in the atmosphere at a specific radiative transfer point
 * \param abs_species The list of all species in the atmosphere
 * 
 * \return A reference pointing at Y, G, or DV.
 * 
 * \author Richard Larsson
 * \date   2018-09-14
 */
Numeric LineFunctionData::GetLineParamDeriv(const Numeric& T0, const Numeric& T,
                                            const Numeric& P, const Numeric& self_vmr,
                                            const ConstVectorView rtp_vmr, 
                                            const ArrayOfSpeciesTag& abs_species,
                                            const RetrievalQuantity& rt, 
                                            const QuantumIdentifier& line_qi,
                                            const bool normalization) const
{
  // Doppler broadening has no types...
  if(mp == LineShapeType::DP) return 0.0;
  
  // Return nothing if this is not supposed to be called by this function
  if(not is_lineshape_lineparam(rt)) return 0.0;
  
  if(not rt.QuantumIdentity().In(line_qi)) return 0.0;
  
  // Find calculation details
  Numeric this_vmr=0.0, total_vmr=0, partial_vmr;
  Index this_derivative=-1;
  for(Index i=0; i<mspecies.nelem(); i++) {
    if(i == 0 and mself) { // The first value might be self-broadening (use self_vmr)
      if(rt.QuantumIdentity().In(line_qi)) {
        this_vmr = self_vmr;
        this_derivative = i;
      }
      partial_vmr = self_vmr;
    }
    else if(i == mspecies.nelem()-1 and mbath) { // The last value might be air-broadening (set total_vmr to 1)
      if(this_derivative == -1) { // if it still is not set by this point, it is going to
        this_vmr = 1 - total_vmr;
        this_derivative = i;
      }
      partial_vmr = 1 - total_vmr;
    }
    else {  // Otherwise we have to find the species in the list of species tags
      Index this_species=-1;
      for(Index j=0; j<abs_species.nelem(); j++) {
        if(abs_species[j].Species() == mspecies[i].Species())
          this_species = j;
      }
      
      if(rt.QuantumIdentity().Species() == mspecies[i].Species() and 
         rt.QuantumIdentity().Isotopologue() == mspecies[i].Isotopologue()) {
        this_vmr = rtp_vmr[this_species];
        this_derivative = i;
      }
      
      // If we cannot find the species, it does not exist and we are in a different atmosphere than covered by the line data
      if(this_species == -1) continue;  // Species does not exist
      partial_vmr = rtp_vmr[this_species];
    }
    
    // Sum up VMR
    total_vmr += partial_vmr;
  }
  
  if(this_derivative == -1) return 0.0;  // Species not in list
  
  Numeric val=0.0;
  Index current=0;
  
  // Make reading easier
  #define x0 mdata[this_derivative][current]
  #define x1 mdata[this_derivative][current+1]
  #define x2 mdata[this_derivative][current+2]
  
  Index skip;
  if(rt.PropMatType() == JacPropMatType::LineFunctionDataG0X0 or
      rt.PropMatType() == JacPropMatType::LineFunctionDataG0X1 or
      rt.PropMatType() == JacPropMatType::LineFunctionDataG0X2) {
    skip = 0;
    if(not ComputesParam("G0")) return 0;
  }
  else if(rt.PropMatType() == JacPropMatType::LineFunctionDataD0X0 or
          rt.PropMatType() == JacPropMatType::LineFunctionDataD0X1 or
          rt.PropMatType() == JacPropMatType::LineFunctionDataD0X2) {
    skip = 1;
    if(not ComputesParam("D0")) return 0;
  }
  else if(rt.PropMatType() == JacPropMatType::LineFunctionDataG2X0 or
          rt.PropMatType() == JacPropMatType::LineFunctionDataG2X1 or
          rt.PropMatType() == JacPropMatType::LineFunctionDataG2X2) {
    skip = 2;
    if(not ComputesParam("G2")) return 0;
  }
  else if(rt.PropMatType() == JacPropMatType::LineFunctionDataD2X0 or
          rt.PropMatType() == JacPropMatType::LineFunctionDataD2X1 or
          rt.PropMatType() == JacPropMatType::LineFunctionDataD2X2) {
    skip = 3;
    if(not ComputesParam("D2")) return 0;
  }
  else if(rt.PropMatType() == JacPropMatType::LineFunctionDataFVCX0 or
          rt.PropMatType() == JacPropMatType::LineFunctionDataFVCX1 or
          rt.PropMatType() == JacPropMatType::LineFunctionDataFVCX2) {
    skip = 4;
    if(not ComputesParam("FVC")) return 0;
  }
  else if(rt.PropMatType() == JacPropMatType::LineFunctionDataETAX0 or
          rt.PropMatType() == JacPropMatType::LineFunctionDataETAX1 or
          rt.PropMatType() == JacPropMatType::LineFunctionDataETAX2) {
    skip = 5;
    if(not ComputesParam("ETA")) return 0;
  }
  else
    skip=-1;
  
  // Skip to the right current position
  if(skip not_eq -1) {
    for(Index j=0; j<skip; j++)
      current += TemperatureTypeNelem(mtypes[this_derivative][j]);
  }
  else { // The derivative is not a HTP parameter but a line mixing parameter
    for(Index j=0; j<LineShapeTypeNelem(); j++)  // skip all of these
      current += TemperatureTypeNelem(mtypes[this_derivative][j]);
  
    // If it is a line mixing parameter we should maybe skip something
    if(rt.PropMatType() == JacPropMatType::LineFunctionDataYX0 or
       rt.PropMatType() == JacPropMatType::LineFunctionDataYX1 or
       rt.PropMatType() == JacPropMatType::LineFunctionDataYX2) {
      skip = 0;
      if(not ComputesParam("Y")) return 0;
    }
    else if(rt.PropMatType() == JacPropMatType::LineFunctionDataGX0 or
            rt.PropMatType() == JacPropMatType::LineFunctionDataGX1 or
            rt.PropMatType() == JacPropMatType::LineFunctionDataGX2) {
      skip = 1;
      if(not ComputesParam("G")) return 0;
    }
    else if(rt.PropMatType() == JacPropMatType::LineFunctionDataDVX0 or
            rt.PropMatType() == JacPropMatType::LineFunctionDataDVX1 or
            rt.PropMatType() == JacPropMatType::LineFunctionDataDVX2) {
      skip = 2;
      if(not ComputesParam("DV")) return 0;
    }
    else 
      throw std::runtime_error("Developer error:  This should not happen");
    
    if(skip not_eq -1) {
      for(Index j=0; j<skip; j++)
        current += TemperatureTypeNelem(mtypes[this_derivative][j+LineShapeTypeNelem()]);
    }
    else 
      throw std::runtime_error("Developer error:  This should not happen");
  }
  
  // Now perform the calculations, we know the type is available and in the below top level switch-statement 
  switch(rt.PropMatType()) {
    case JacPropMatType::LineFunctionDataG0X0:  
    case JacPropMatType::LineFunctionDataD0X0:
    case JacPropMatType::LineFunctionDataG2X0:
    case JacPropMatType::LineFunctionDataD2X0:
    case JacPropMatType::LineFunctionDataETAX0:
    case JacPropMatType::LineFunctionDataFVCX0:
    case JacPropMatType::LineFunctionDataYX0:
    case JacPropMatType::LineFunctionDataGX0:
    case JacPropMatType::LineFunctionDataDVX0:
      switch(mtypes[this_derivative][skip]) {
        case TemperatureType::None:
        case TemperatureType::LM_AER:
          val = 0;
          break;
        case TemperatureType::T0:
          val = dmain_dx0_t0();
          break;
        case TemperatureType::T1:
          val = dmain_dx0_t1(T, T0, x1);
          break;
        case TemperatureType::T2:
          val = dmain_dx0_t2(T, T0, x1, x2);
          break;
        case TemperatureType::T3:
          val = dmain_dx0_t3();
          break;
        case TemperatureType::T4:
          val = dmain_dx0_t4(T, T0, x2);
          break;
        case TemperatureType::T5:
          val = dmain_dx0_t5(T, T0, x1);
          break;
      }
      break;
    case JacPropMatType::LineFunctionDataG0X1:
    case JacPropMatType::LineFunctionDataD0X1:
    case JacPropMatType::LineFunctionDataG2X1:
    case JacPropMatType::LineFunctionDataD2X1:
    case JacPropMatType::LineFunctionDataETAX1:
    case JacPropMatType::LineFunctionDataFVCX1:
    case JacPropMatType::LineFunctionDataYX1:
    case JacPropMatType::LineFunctionDataGX1:
    case JacPropMatType::LineFunctionDataDVX1:
      switch(mtypes[this_derivative][skip]) {
        case TemperatureType::None:
        case TemperatureType::LM_AER:
          val = 0;
          break;
        case TemperatureType::T0:
          val = dmain_dx1_t0();
          break;
        case TemperatureType::T1:
          val = dmain_dx1_t1(T, T0, x0, x1);
          break;
        case TemperatureType::T2:
          val = dmain_dx1_t2(T, T0, x0, x1, x2);
          break;
        case TemperatureType::T3:
          val = dmain_dx1_t3(T, T0);
          break;
        case TemperatureType::T4:
          val = dmain_dx1_t4(T, T0, x2);
          break;
        case TemperatureType::T5:
          val = dmain_dx1_t5(T, T0, x0, x1);
          break;
      }
      break;
    case JacPropMatType::LineFunctionDataG0X2:
    case JacPropMatType::LineFunctionDataD0X2:
    case JacPropMatType::LineFunctionDataG2X2:
    case JacPropMatType::LineFunctionDataD2X2:
    case JacPropMatType::LineFunctionDataETAX2:
    case JacPropMatType::LineFunctionDataFVCX2:
    case JacPropMatType::LineFunctionDataYX2:
    case JacPropMatType::LineFunctionDataGX2:
    case JacPropMatType::LineFunctionDataDVX2:
      switch(mtypes[this_derivative][skip]) {
        case TemperatureType::None:
        case TemperatureType::LM_AER:
          val = 0;
          break;
        case TemperatureType::T0:
          val = dmain_dx2_t0();
          break;
        case TemperatureType::T1:
          val = dmain_dx2_t1();
          break;
        case TemperatureType::T2:
          val = dmain_dx2_t2(T, T0, x0, x1);
          break;
        case TemperatureType::T3:
          val = dmain_dx2_t3();
          break;
        case TemperatureType::T4:
          val = dmain_dx2_t4(T, T0, x0, x1, x2);
          break;
        case TemperatureType::T5:
          val = dmain_dx2_t5();
          break;
      }
      break;
    default: 
      throw std::runtime_error("Developer error:  This should not happen");
  }
  
  #undef x0
  #undef x1
  #undef x2
  
  switch(rt.PropMatType()){
    case JacPropMatType::LineFunctionDataGX0:
    case JacPropMatType::LineFunctionDataDVX0:
    case JacPropMatType::LineFunctionDataGX1:
    case JacPropMatType::LineFunctionDataDVX1:
    case JacPropMatType::LineFunctionDataGX2:
    case JacPropMatType::LineFunctionDataDVX2:
      val *= this_vmr * P;
      break;
    default: {/*pass*/}
  }
  
  if(normalization)
    return val * this_vmr / total_vmr * P;
  else 
    return val * this_vmr * P;
}


//! Prints data so as operator>> can read it
std::ostream& operator<<(std::ostream& os, const LineFunctionData& lfd) {
  // Keep track of the size of the problem
  const Index nshapes=lfd.LineShapeTypeNelem(), nmixing=lfd.LineMixingTypeNelem();
  
  // Init the problem
  os << lfd.LineShapeType2String() << " " << lfd.LineMixingType2String() << " " << lfd.mspecies.nelem() << " ";
  
  // For all species we should now do mostly the same
  for(Index i=0; i<lfd.mspecies.nelem(); i++) {
    if(i==0 and lfd.mself)
      os << SELF << " ";  // SELF can be the first species
    else if(i==lfd.mspecies.nelem()-1 and lfd.mbath)
      os << BATH << " ";  // BATH can be the last species
    else
      os << lfd.mspecies[i].SpeciesNameMain() << " ";  // Otherwise we have a species defined in the assoc. SpeciesTag
    
    // Now we must count the data
    Index counter=0;
    
    // For everything that relates to shapes, do the same thing
    for(Index j=0; j<nshapes; j++) {
      os << lfd.TemperatureType2String(lfd.mtypes[i][j]) << " ";  // Name the temperature type
      for(Index k=0; k<lfd.TemperatureTypeNelem(lfd.mtypes[i][j]); k++)
        os << lfd.mdata[i][counter+k] << " ";  // Print the assoc. data
      counter += lfd.TemperatureTypeNelem(lfd.mtypes[i][j]);  // Count how much data has been printed
    }
    
    // For mixing we must take care with the interp. values but otherwise we act the same as prev. loop
    for(Index j=nshapes; j<nmixing+nshapes; j++) {
      os << lfd.TemperatureType2String(lfd.mtypes[i][j]) << " ";
      if(lfd.mlm == LineFunctionData::LineMixingOrderType::Interp) {
        counter++;
        os << lfd.mdata[i].nelem() - counter << " ";
        for(; counter < lfd.mdata[i].nelem(); counter++)
          os << lfd.mdata[i][counter] << " ";
      }
      else {
        for(Index k=0; k<lfd.TemperatureTypeNelem(lfd.mtypes[i][j]); k++)
          os << lfd.mdata[i][counter+k] << " ";
        counter += lfd.TemperatureTypeNelem(lfd.mtypes[i][j]);
      }
    }
  }
  
  return os;
}


//! Reads data as created by operator<< can read it
std::istream& operator>>(std::istream& data, LineFunctionData& lfd) {
  lfd.mself = lfd.mbath = false;
  Index specs, c;
  Numeric n;
  String s;
  
  // The first tag should give the line shape scheme
  data >> s;
  lfd.StringSetLineShapeType(s);
  
  // The second tag should give the line mixing scheme
  data >> s;
  lfd.StringSetLineMixingType(s);
  
  // From the line shape and line mixing types, we know how many parameters are needed
  const Index count = lfd.LineShapeTypeNelem() + lfd.LineMixingTypeNelem();
  
  // The third tag should contain the number of species
  data >> specs;
  lfd.mspecies.resize(specs);
  lfd.mtypes = Array<Array<LineFunctionData::TemperatureType>>(specs, Array<LineFunctionData::TemperatureType>(count));
  lfd.mdata.resize(specs);
  
  if(lfd.mp not_eq LineFunctionData::LineShapeType::DP and not specs)
    throw std::runtime_error("Need at least one species for non-Doppler line shapes");
  
  // For all species, we need to set the methods to compute them
  for(Index i=0; i<specs; i++) {
    
    // This should be a species tag or one of the specials, SELF or BATH
    data >> s;
    if(s == SELF) {
      // If the species is self, then  we need to flag this
      lfd.mself = true;
      if(i not_eq 0)  // but self has to be first for consistent behavior
        throw std::runtime_error("Self broadening must be first, it is not\n");
    }
    else if(s == BATH) {
      // If the species is air, then we need to flag this
      lfd.mbath = true;
      if(i not_eq specs - 1)  // but air has to be last because it needs the rest's VMR
        throw std::runtime_error("Air/bath broadening must be last, it is not\n");
    }
    else {
      // Otherwise, we hope we find a species
      try {
        lfd.mspecies[i] = SpeciesTag(s);
      }
      catch(const std::runtime_error& e) {
        ostringstream os;
        os << "Encountered " << s << " in a position where a species should have been ";
        os << "defined.\nPlease check your pressure broadening data structure and ensure ";
        os << "that it follows the correct conventions.\n";
        os << "SpeciesTag error reads:  " << e.what();
        throw std::runtime_error(os.str());
      }
    }
    
    ArrayOfNumeric nums; nums.reserve(20); // buffers
    
    // For all parameters
    for(Index j=0; j<count; j++) {
      data >> s;  // Should contain a temperature tag
      lfd.StringSetTemperatureType(i, j, s);
      
      // Find the count of species, negative numbers means the next data value has this number
      c = lfd.TemperatureTypeNelem(lfd.mtypes[i][j]);
      if(c < 0)
        data >> c;
      
      // Add all new numbers for this line to the numbers
      for(Index k=0; k<c; k++) {
        data >> double_imanip() >> n;
        nums.push_back(n);
      }
    }
    
    // Set the data now that we know how many counts are required
    lfd.mdata[i].resize(nums.nelem());
    for(Index j=0; j<nums.nelem(); j++)
      lfd.mdata[i][j] = nums[j];
  }
  
  // Set the variable that determines whether or not 
  // the line has to be used under non-standard calculations,
  // such as full line mixing.  In case this line needs to be
  // computed separately, a special catalog that references
  // this line needs to be created, and this catalog
  // should ignore the variable below
  lfd.do_line_in_standard_calculations = true;
  
  return data;
}


LineFunctionData::LineFunctionData(const PressureBroadeningData& pb,
                                   const LineMixingData& lm,
                                   const String& species,
                                   const Numeric& T0)
{
  const PressureBroadeningData::PB_Type pbType = pb.Type();
  const LineMixingData::LM_Type lmType = lm.Type();
  
  Vector pbData, lmData, lmTranslated, pbTranslated;
  pb.GetVectorFromData(pbData);
  lm.GetVectorFromData(lmData);
  
  mself = false;
  mp = LineShapeType::VP;
  
  // Fox LM type and set the line mixing data
  switch(lmType) {
    case LineMixingData::LM_LBLRTM:
      do_line_in_standard_calculations = true;
      mlm = LineMixingOrderType::Interp;
      lmTranslated = std::move(lmData);
      break;
    case LineMixingData::LM_1STORDER:
      do_line_in_standard_calculations = true;
      mlm = LineMixingOrderType::LM1;
      lmTranslated = {lmData[1], lmData[2]};
      if(T0 not_eq lmData[0])
        throw std::runtime_error("Cannot translate data of line since it has bad temperature\n"
        "information.  The line would produce poor absorption profiles anyways, so please\n"
        "reform it in the new format of line catalogs.");
      break;
    case LineMixingData::LM_2NDORDER:
      do_line_in_standard_calculations = true;
      mlm = LineMixingOrderType::LM2;
      lmTranslated = {lmData[0], lmData[1], lmData[7], 
                      lmData[2], lmData[3], lmData[8],
                      lmData[4], lmData[5], lmData[9]};
      if(T0 not_eq lmData[6])
        throw std::runtime_error("Cannot translate data of line since it has bad temperature\n"
        "information.  The line would produce poor absorption profiles anyways, so please\n"
        "reform it in the new format of line catalogs.");
      break;
    case LineMixingData::LM_NONE:
      do_line_in_standard_calculations = true;
      mlm = LineMixingOrderType::None;
      lmTranslated = {};
      break;
    case LineMixingData::LM_BYBAND:
      do_line_in_standard_calculations = false;
      mlm = LineMixingOrderType::None;
      lmTranslated = {};
      break;
    default:
      throw std::runtime_error("Error, unsupported conversion.  Please update to\n"
      "new line catalog format manually.  Cannot support non-resonant LM as\n"
      "it was a hack to begin with.");
  }
  
  // Give an Index number to self if in the list of species for ARTSCAT-4
  const Index speciesID = SpeciesTag("N2" ).IsSpecies(species) ? 0 :
                          SpeciesTag("O2" ).IsSpecies(species) ? 1 :
                          SpeciesTag("H2O").IsSpecies(species) ? 2 :
                          SpeciesTag("CO2").IsSpecies(species) ? 3 :
                          SpeciesTag("H2" ).IsSpecies(species) ? 4 :
                          SpeciesTag("He" ).IsSpecies(species) ? 5 : - 1;

  // Fix mspecies and mdata for pressure broadening
  Index inx=0;
  switch(pbType) {
    case PressureBroadeningData::PB_AIR_BROADENING:
      mbath = true;
      mdata = ArrayOfVector(2, Vector(4+lmTranslated.nelem()));
      mspecies = {SpeciesTag(species), SpeciesTag()};
      
      // For SELF broadening
      pbTranslated = {pbData[0], pbData[1], 0, pbData[1]};
      mdata[inx][Range(0, 4)] = pbTranslated; inx++;
      
      // For AIR broadening
      pbTranslated = {pbData[2], pbData[3], pbData[4], pbData[3]};
      mdata[inx][Range(0, 4)] = pbTranslated; inx++;
      break;
    case PressureBroadeningData::PB_AIR_AND_WATER_BROADENING:
      mbath = true;
      
      if(speciesID == 2) { // Water
        mdata = ArrayOfVector(2, Vector(4+lmTranslated.nelem()));
        mspecies = {SpeciesTag("H2O"), SpeciesTag()};
        
        // For Water broadening
        pbTranslated = {pbData[6], pbData[7], pbData[8], pbData[7]};
        mdata[inx][Range(0, 4)] = pbTranslated; inx++;
        
        // For AIR broadening
        pbTranslated = {pbData[3], pbData[4], pbData[5], pbData[4]};
        mdata[inx][Range(0, 4)] = pbTranslated; inx++;
      }
      else { // not Water
        mdata = ArrayOfVector(3, Vector(4+lmTranslated.nelem()));
        mspecies = {SpeciesTag(species), SpeciesTag("H2O"), SpeciesTag()};
        
        // For SELF broadening
        pbTranslated = {pbData[0], pbData[1], pbData[2], pbData[1]};
        mdata[inx][Range(0, 4)] = pbTranslated; inx++;
        
        // For Water broadening
        pbTranslated = {pbData[6], pbData[7], pbData[8], pbData[7]};
        mdata[inx][Range(0, 4)] = pbTranslated; inx++;
        
        // For AIR broadening
        pbTranslated = {pbData[3], pbData[4], pbData[5], pbData[4]};
        mdata[inx][Range(0, 4)] = pbTranslated; inx++;
      }
      break;
    case PressureBroadeningData::PB_PLANETARY_BROADENING:
      mbath = false;
      if(speciesID not_eq -1) {
        mdata = ArrayOfVector(6, Vector(4+lmTranslated.nelem()));
        mspecies = {SpeciesTag("N2" ), SpeciesTag("O2" ), SpeciesTag("H2O"), SpeciesTag("CO2"), SpeciesTag("H2" ), SpeciesTag("He" )};
      }
      else {
        mdata = ArrayOfVector(7, Vector(4+lmTranslated.nelem()));
        mspecies = {SpeciesTag(species), SpeciesTag("N2"   ), SpeciesTag("O2"   ), SpeciesTag("H2O"  ), SpeciesTag("CO2"  ), SpeciesTag("H2"   ), SpeciesTag("He"   )};
        
        // For SELF broadening
        pbTranslated = {pbData[0], pbData[7], 0, pbData[7]};
        mdata[inx][Range(0, 4)] = pbTranslated; inx++;
      }
        
      // For N2 broadening
      pbTranslated = {pbData[1], pbData[8], pbData[14], pbData[8]};
      mdata[inx][Range(0, 4)] = pbTranslated; inx++;
      
      // For O2 broadening
      pbTranslated = {pbData[2], pbData[9], pbData[15], pbData[9]};
      mdata[inx][Range(0, 4)] = pbTranslated; inx++;
      
      // For H2O broadening
      pbTranslated = {pbData[3], pbData[10], pbData[16], pbData[10]};
      mdata[inx][Range(0, 4)] = pbTranslated; inx++;
      
      // For CO2 broadening
      pbTranslated = {pbData[4], pbData[11], pbData[17], pbData[11]};
      mdata[inx][Range(0, 4)] = pbTranslated; inx++;
      
      // For H2 broadening
      pbTranslated = {pbData[5], pbData[12], pbData[18], pbData[12]};
      mdata[inx][Range(0, 4)] = pbTranslated; inx++;
      
      // For He broadening
      pbTranslated = {pbData[6], pbData[13], pbData[19], pbData[13]};
      mdata[inx][Range(0, 4)] = pbTranslated; inx++;
      break;
    default:
      throw std::runtime_error("Error, unsupported conversion.  Please update to\n"
      "new line catalog format manually.  Only air, water+air, and planetary\n"
      "conversions are accepted as the rest where experimental at time of\n"
      "implementing the new code.");
  }
  
  // Apply line mixing data equally everywhere since we did not have line mixing data for different species before
  if(lmTranslated.nelem()) for(auto& d: mdata) d[Range(4, joker)] = lmTranslated;
  
  // Fix mtypes
  mtypes.resize(mdata.nelem());
  switch(mlm) {
    case LineMixingOrderType::None:
      for(auto& d: mtypes) d = {TemperatureType::T1, TemperatureType::T5};
      break;
    case LineMixingOrderType::LM1:
      for(auto& d: mtypes) d = {TemperatureType::T1, TemperatureType::T5, TemperatureType::T1};
      break;
    case LineMixingOrderType::LM2:
      for(auto& d: mtypes) d = {TemperatureType::T1, TemperatureType::T5, TemperatureType::T4, TemperatureType::T4, TemperatureType::T4};
      break;
    case LineMixingOrderType::Interp:
      for(auto& d: mtypes) d = {TemperatureType::T1, TemperatureType::T5, TemperatureType::LM_AER};
      break;
  }
}
