/* Copyright (C) 2014
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


#include "zeeman.h"
#include "linefunctions.h"
#include "linescaling.h"
#include "species_info.h"


void zeeman_on_the_fly(ArrayOfPropagationMatrix& propmat_clearsky, 
                       ArrayOfStokesVector& nlte_source,
                       ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                       ArrayOfStokesVector& dnlte_dx_source,
                       ArrayOfStokesVector& nlte_dsource_dx,
                       const ArrayOfArrayOfSpeciesTag& abs_species, 
                       const ArrayOfRetrievalQuantity& flag_partials,
                       const ArrayOfArrayOfLineRecord& zeeman_linerecord_precalc,
                       const SpeciesAuxData& isotopologue_ratios, 
                       const SpeciesAuxData& partition_functions,
                       const ConstVectorView f_grid,
                       const ConstVectorView rtp_vmrs, 
                       const ConstVectorView rtp_nlte, 
                       const ConstVectorView rtp_mag,
                       const ConstVectorView rtp_los,
                       const Numeric& rtp_pressure,
                       const Numeric& rtp_temperature,
                       const Index& manual_zeeman_tag,
                       const Numeric& manual_zeeman_magnetic_field_strength,
                       const Numeric& manual_zeeman_theta,
                       const Numeric& manual_zeeman_eta,
                       const Verbosity& verbosity)
{
  // Find relevant derivatives in retrieval quantities positions
  const ArrayOfIndex flag_partials_positions = equivlent_propmattype_indexes(flag_partials);
  
  // Size of problem
  const Index nf = f_grid.nelem();
  const Index nq = flag_partials_positions.nelem();
  const Index ns = abs_species.nelem();
  const Index nn = rtp_nlte.nelem();
  
  // Pressure information
  const Numeric dnumdens_dmvr = number_density(rtp_pressure, rtp_temperature);
  const Numeric dnumdens_dt_dmvr = dnumber_density_dt(rtp_pressure, rtp_temperature);
  
  // Main compute vectors
  ComplexVector F(nf), N(nn ? nf : 0);
  ComplexMatrix dF(nq, nf), dN(nn ? nq : 0, nf);
  Range frange(joker);
  
  // Magnetic field variables
  const Numeric u = rtp_mag[0];
  const Numeric v = rtp_mag[1];
  const Numeric w = rtp_mag[2];
  const Numeric z = DEG2RAD * rtp_los[0];
  const Numeric a = DEG2RAD * rtp_los[1];
  const Numeric H =     manual_zeeman_tag ? manual_zeeman_magnetic_field_strength : zeeman_magnetic_magnitude(u, v, w      );
  const Numeric eta =   manual_zeeman_tag ? manual_zeeman_eta                     : zeeman_magnetic_eta(      u, v, w, z, a);
  const Numeric theta = manual_zeeman_tag ? manual_zeeman_theta                   : zeeman_magnetic_theta(    u, v, w, z, a);
  
  // Magnetic field derivatives... FIXME:  Deal with asymptotes! (e.g., for theta == 90)
  const bool do_mag_jacs = do_magnetic_jacobian(flag_partials) and not manual_zeeman_tag;
  const Numeric dH_du     = do_mag_jacs ? zeeman_magnetic_dmagnitude_du(u, v, w      ) : 0;
  const Numeric dH_dv     = do_mag_jacs ? zeeman_magnetic_dmagnitude_dv(u, v, w      ) : 0;
  const Numeric dH_dw     = do_mag_jacs ? zeeman_magnetic_dmagnitude_dw(u, v, w      ) : 0;
  const Numeric deta_du   = do_mag_jacs ? zeeman_magnetic_deta_du(      u, v, w, z, a) : 0;
  const Numeric deta_dv   = do_mag_jacs ? zeeman_magnetic_deta_dv(      u, v, w, z, a) : 0;
  const Numeric deta_dw   = do_mag_jacs ? zeeman_magnetic_deta_dw(      u, v, w, z, a) : 0;
  const Numeric dtheta_du = do_mag_jacs ? zeeman_magnetic_dtheta_du(    u, v, w, z, a) : 0;
  const Numeric dtheta_dv = do_mag_jacs ? zeeman_magnetic_dtheta_dv(    u, v, w, z, a) : 0;
  const Numeric dtheta_dw = do_mag_jacs ? zeeman_magnetic_dtheta_dw(    u, v, w, z, a) : 0;
  
  for(Index ispecies=0; ispecies<ns; ispecies++) {
    for(const ArrayOfLineRecord& lines: zeeman_linerecord_precalc) {
      if(not lines.nelem()) continue;
      else if(lines[0].Species()      not_eq abs_species[ispecies][0].Species() or
              lines[0].Isotopologue() not_eq abs_species[ispecies][0].Isotopologue()) continue;
      
      // Polarization
      const Vector polarization_scale = lines[0].ZeemanEffect().Polarization(theta, eta);
      const Vector dpol_deta = do_mag_jacs ? lines[0].ZeemanEffect().dPolarization_deta(theta, eta) : Vector(0);
      const Vector dpol_dtheta = do_mag_jacs ? lines[0].ZeemanEffect().dPolarization_dtheta(theta, eta) : Vector(0);
      
      // Temperature constants
      Numeric t0=-1.0, qt, qt0, dqt_dT;
      const Numeric numdens = rtp_vmrs[ispecies] * dnumdens_dmvr;
      const Numeric dnumdens_dT = rtp_vmrs[ispecies] * dnumdens_dt_dmvr;
      const Numeric dc = Linefunctions::DopplerConstant(rtp_temperature, lines[0].IsotopologueData().Mass());
      const Numeric ddc_dT = Linefunctions::dDopplerConstant_dT(rtp_temperature, lines[0].IsotopologueData().Mass());
      const Numeric isotop_ratio = isotopologue_ratios.getParam(lines[0].Species(), lines[0].Isotopologue())[0].data[0];
      const Numeric partial_pressure = rtp_pressure * rtp_vmrs[ispecies];
      
      for(const LineRecord& line: lines) {
        if(line.Ti0() not_eq t0) {
          t0 = line.Ti0();
          
          partition_function(qt0, qt, t0, rtp_temperature,
                             partition_functions.getParamType(line.Species(), line.Isotopologue()),
                             partition_functions.getParam(line.Species(), line.Isotopologue()));
          
          if(do_temperature_jacobian(flag_partials))
            dpartition_function_dT(dqt_dT, qt, rtp_temperature, temperature_perturbation(flag_partials),
                                   partition_functions.getParamType(line.Species(), line.Isotopologue()),
                                   partition_functions.getParam(line.Species(), line.Isotopologue()));
        }
        
        for(Index iz=0; iz<line.ZeemanEffect().nelem(); iz++) {
          const Numeric B     = nn ? planck(    line.F(), rtp_temperature) : 0;
          const Numeric dB_dT = nn ? dplanck_dt(line.F(), rtp_temperature) : 0;
          
          Linefunctions::set_cross_section_for_single_line(F, dF, N, dN, frange,
                                                           flag_partials, flag_partials_positions, line, f_grid, rtp_vmrs, 
                                                           rtp_nlte, rtp_pressure, rtp_temperature, dc, partial_pressure, 
                                                           isotop_ratio, H, ddc_dT, qt, dqt_dT, qt0,
                                                           abs_species, ispecies, iz, verbosity);
          
          // range-based arguments that need be made to work for both complex and numeric
          const Index extent = (frange.get_extent() < 0) ? (nf - frange.get_start()) : frange.get_extent();
          const Range this_out_range(frange.get_start(), extent);
          
          const ComplexVectorView F_range_view  =       F[frange];
          const ComplexVectorView N_range_view  = nn ?  N[frange] : F;
          const ComplexMatrixView dF_range_view =      dF(joker, frange);
          const ComplexMatrixView dN_range_view = nn ? dN(joker, frange) : dF;
          
          for(Index i = 0; i < extent; i++) {
            const Index iv = frange(i);
            const Complex& CF = F_range_view[i];
            const Complex& CN = N_range_view[i];
            
            propmat_clearsky[ispecies].AddPolarized(polarization_scale, iv, CF*numdens);
            if(nn)
              nlte_source[ispecies].AddPolarized(polarization_scale, iv, CN*numdens*B);
            
            for(Index j=0; j<nq; j++) {
              const Complex& dCF = dF_range_view(j, i);
              const Complex& dCN = dN_range_view(j, i);
              
              if(flag_partials[flag_partials_positions[j]] == JacPropMatType::Temperature) {
                dpropmat_clearsky_dx[j].AddPolarized(polarization_scale, iv, dCF*numdens + CF*dnumdens_dT);
                if(nn) {
                  dnlte_dx_source[j].AddPolarized(polarization_scale, iv, (CN*dnumdens_dT + dCN*numdens)*B);
                  nlte_dsource_dx[j].AddPolarized(polarization_scale, iv, CN*numdens*dB_dT);
                }
              }
              else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::VMR) {
                if(flag_partials[flag_partials_positions[j]].QuantumIdentity() < line.QuantumIdentity()) {
                  dpropmat_clearsky_dx[j].AddPolarized(polarization_scale, iv, dCF*numdens + CF*dnumdens_dmvr);
                  if(nn)
                    dnlte_dx_source[j].AddPolarized(polarization_scale, iv, (dCN*numdens + CN*dnumdens_dmvr)*B);
                }
              }
              else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticEta) {
                dpropmat_clearsky_dx[j].AddPolarized(dpol_deta, iv, CF*numdens);
                if(nn)
                  dnlte_dx_source[j].AddPolarized(dpol_deta, iv, CN*numdens*B);
              }
              else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticTheta) {
                dpropmat_clearsky_dx[j].AddPolarized(dpol_dtheta, iv, CF*numdens);
                if(nn)
                  dnlte_dx_source[j].AddPolarized(dpol_dtheta, iv, CN*numdens*B);
              }
              else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticU) {
                dpropmat_clearsky_dx[j].AddPolarized(polarization_scale, iv, dCF*numdens*dH_du);
                dpropmat_clearsky_dx[j].AddPolarized(dpol_deta,          iv,  CF*numdens*deta_du);
                dpropmat_clearsky_dx[j].AddPolarized(dpol_dtheta,        iv,  CF*numdens*dtheta_du);
                if(nn) {
                  dnlte_dx_source[j].AddPolarized(polarization_scale, iv, dCN*numdens*B*dH_du);
                  dnlte_dx_source[j].AddPolarized(dpol_deta,          iv,  CN*numdens*B*deta_du);
                  dnlte_dx_source[j].AddPolarized(dpol_dtheta,        iv,  CN*numdens*B*dtheta_du);
                }
              }
              else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticV) {
                dpropmat_clearsky_dx[j].AddPolarized(polarization_scale, iv, dCF*numdens*dH_dv);
                dpropmat_clearsky_dx[j].AddPolarized(dpol_deta,          iv,  CF*numdens*deta_dv);
                dpropmat_clearsky_dx[j].AddPolarized(dpol_dtheta,        iv,  CF*numdens*dtheta_dv);
                if(nn) {
                  dnlte_dx_source[j].AddPolarized(polarization_scale, iv, dCN*numdens*B*dH_dv);
                  dnlte_dx_source[j].AddPolarized(dpol_deta,          iv,  CN*numdens*B*deta_dv);
                  dnlte_dx_source[j].AddPolarized(dpol_dtheta,        iv,  CN*numdens*B*dtheta_dv);
                }
              }
              else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticW) {
                dpropmat_clearsky_dx[j].AddPolarized(polarization_scale, iv, dCF*numdens*dH_dw);
                dpropmat_clearsky_dx[j].AddPolarized(dpol_deta,          iv,  CF*numdens*deta_dw);
                dpropmat_clearsky_dx[j].AddPolarized(dpol_dtheta,        iv,  CF*numdens*dtheta_dw);
                if(nn) {
                  dnlte_dx_source[j].AddPolarized(polarization_scale, iv, dCN*numdens*B*dH_dw);
                  dnlte_dx_source[j].AddPolarized(dpol_deta,          iv,  CN*numdens*B*deta_dw);
                  dnlte_dx_source[j].AddPolarized(dpol_dtheta,        iv,  CN*numdens*B*dtheta_dw);
                }
              }
              else {
                dpropmat_clearsky_dx[j].AddPolarized(polarization_scale, iv, dCF*numdens);
                if(nn)
                  dnlte_dx_source[j].AddPolarized(polarization_scale, iv, dCN*numdens*B);
              }
            }
          }
        }
      }
    }
  }
}


void set_magnetic_parameters(Numeric& H_mag,
                             Numeric& eta,
                             Numeric& theta,
                             const Index manual_zeeman_tag,
                             const Numeric& manual_zeeman_eta,
                             const Numeric& manual_zeeman_theta,
                             const Numeric& manual_zeeman_magnetic_field_strength,
                             ConstVectorView rtp_mag,
                             ConstVectorView r_path_los)
{
//Get the magnitude of the magnetic field and store a local unit Vector for simplified angle calculations.
  H_mag = manual_zeeman_tag != 0?manual_zeeman_magnetic_field_strength:sqrt( rtp_mag * rtp_mag );

  if(manual_zeeman_tag!=0)
  { // Leaving it up to the user to manually tag the angles for simplified magnetic fields.
    eta   = manual_zeeman_eta;
    theta = manual_zeeman_theta;
  }
  else if(H_mag==0.0)
  {
      eta = 0.0;
      theta = 0.0;
  }
  else
  { 
    const Numeric 
    aa=DEG2RAD*r_path_los[1], 
    za=DEG2RAD*r_path_los[0], 
    Bu = rtp_mag[0], 
    Bv = rtp_mag[1], 
    Bw = rtp_mag[2],
    cosaa=cos(aa),
    cosza=cos(za),
    sinaa=sin(aa),
    sinza=sin(za),
    cosaa2 = cosaa*cosaa,
    //cosza2 = cosza*cosza,
    sinza2 = sinza*sinza,
    sinaa2 = sinaa*sinaa,
    Bu2 = Bu*Bu,
    Bv2 = Bv*Bv,
    Bw2 = Bw*Bw,
    H2  = H_mag*H_mag,
    term1=(Bu*(cosaa2*sinza2 - 1.) + Bw*cosaa*cosza*sinza + Bv*cosaa*sinaa*sinza2),
    term2=(Bv*(sinaa2*sinza2 - 1.) + Bw*sinaa*cosza*sinza + Bu*cosaa*sinaa*sinza2),
    term3=( - Bw*sinza2            + Bu*cosaa*cosza*sinza + Bv*sinaa*cosza*sinza),
    eta_test=(PI - acos(((Bu*cosaa*cosza - Bw*sinza + Bv*sinaa*cosza)*sqrt(term1*term1 + term2*term2 +term3*term3))/(H2)))*RAD2DEG,
    x1 = (Bv*cosaa - Bu*sinaa),
    x2 = -((Bu2*cosaa-Bv2*cosaa+2.0*Bu*Bv*sinaa)*cosaa*sinza2 - Bu2 - Bw2 + 2.0*Bu*Bw*cosaa*cosza*sinza + (Bw2*cosza - Bv2*cosza+ 2.0*Bv*Bw*sinaa*sinza)*cosza),
    fx2=sqrt(x2),
    x3 = 1.0/H2,
    x=x1*fx2*x3;
    
    theta = acos((Bw*cosza + Bu*cosaa*sinza + Bv*sinaa*sinza)/H_mag) * RAD2DEG;
    
    if ((abs(x)-(Numeric)1.0)>0.0) // Numerical drifts can cause this...
      eta = 0.0;
    else 
    {
      eta=acos(x)*RAD2DEG;
      
      if(eta_test>90.0) eta*=-1.0;
    }
    
//      eta = RAD2DEG * zeeman_magnetic_eta(Bu, Bv, Bw, za, aa);
//     theta = RAD2DEG * zeeman_magnetic_theta(Bu, Bv, Bw, za, aa);
  }
}


void create_Zeeman_linerecordarrays(ArrayOfArrayOfLineRecord& aoaol,
                                    const ArrayOfArrayOfSpeciesTag& abs_species,
                                    const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                    const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  const static Numeric margin = 1e-4;
  
  // For all species
  for(Index II = 0; II<abs_species.nelem(); II++) {
    // If the species isn't Zeeman, look at the next species
    if(!is_zeeman(abs_species[II])) continue;
    
    // If there are no lines give up on this species
    const Index nlines = abs_lines_per_species[II].nelem();
    if(not nlines) continue;
    
    // One line record array per type of polarizer is created
    aoaol.push_back(abs_lines_per_species[II]); // First is negative
    aoaol.push_back(abs_lines_per_species[II]); // Second is 0
    aoaol.push_back(abs_lines_per_species[II]); // Third is positive
    
    ArrayOfLineRecord& temp_abs_lines_sm = aoaol[aoaol.nelem()-3]; // sigma minus
    ArrayOfLineRecord& temp_abs_lines_pi = aoaol[aoaol.nelem()-2]; // pi
    ArrayOfLineRecord& temp_abs_lines_sp = aoaol[aoaol.nelem()-1]; // sigma plus
    
    // Else loop over all the lines in the species.
    for (Index ii = 0; ii < nlines; ii++)
    {
      // local LineRecord
      temp_abs_lines_sm[ii].SetZeemanEffectData(ZeemanEffectData(temp_abs_lines_sm[ii].QuantumIdentity(), ZeemanPolarizationType::SigmaMinus));
      temp_abs_lines_pi[ii].SetZeemanEffectData(ZeemanEffectData(temp_abs_lines_pi[ii].QuantumIdentity(), ZeemanPolarizationType::Pi));
      temp_abs_lines_sp[ii].SetZeemanEffectData(ZeemanEffectData(temp_abs_lines_sp[ii].QuantumIdentity(), ZeemanPolarizationType::SigmaPlus));
      
      const Numeric rel_str = temp_abs_lines_sp[ii].ZeemanEffect().SumStrengthScale() + temp_abs_lines_sm[ii].ZeemanEffect().SumStrengthScale() + temp_abs_lines_pi[ii].ZeemanEffect().SumStrengthScale();
      if(abs(rel_str - 1) > margin) {
        ostringstream os;
        os << "The lines\n" << temp_abs_lines_sm[ii] << temp_abs_lines_pi[ii] << temp_abs_lines_sp[ii]
           << "The relative strength should be 1.0 but is " << rel_str << "\n";
        throw std::runtime_error(os.str());
      }
    }
  }
}

Index part_mag_strength(const ArrayOfRetrievalQuantity& flag_partials)
{
    for(Index ii=0; ii<flag_partials.nelem(); ii++)
        if(flag_partials[ii].MainTag() == "Zeeman" &&  flag_partials[ii].Subtag() == "Magnetic Strength" && flag_partials[ii].SubSubtag() == "From Propagation")
            return ii;
    return -1;
}

Index part_mag_theta(const ArrayOfRetrievalQuantity& flag_partials)
{
    for(Index ii=0; ii<flag_partials.nelem(); ii++)
        if(flag_partials[ii].MainTag() == "Zeeman" &&  flag_partials[ii].Subtag() == "Magnetic Theta" && flag_partials[ii].SubSubtag() == "From Propagation")
            return ii;
    return -1;
}

Numeric zeeman_magnetic_magnitude(const Numeric& u, const Numeric& v, const Numeric& w)
{
  /*
    H = ||{u, v, w}||
  */
  return sqrt(u*u + v*v + w*w);
}

Numeric zeeman_magnetic_dmagnitude_du(const Numeric& u, const Numeric& v, const Numeric& w)
{
  return u / sqrt(u*u + v*v + w*w);
}

Numeric zeeman_magnetic_dmagnitude_dv(const Numeric& u, const Numeric& v, const Numeric& w)
{
  return v / sqrt(u*u + v*v + w*w);
}

Numeric zeeman_magnetic_dmagnitude_dw(const Numeric& u, const Numeric& v, const Numeric& w)
{
  return w / sqrt(u*u + v*v + w*w);
}

Numeric zeeman_magnetic_theta(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a)
{
  /*
               ( / cos(a) sin(z) \   / u \   // || / u \ || )
   theta = acos( | sin(a) sin(z) | o | v |  //  || | v | || )
               ( \        cos(z) /   \ w / //   || \ w / || )
  */
  return acos((u*sin(z)*cos(a) + v*sin(a)*sin(z) + w*cos(z))/sqrt(u*u + v*v + w*w));
}

Numeric zeeman_magnetic_dtheta_du(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a)
{
  // autogenerated with sympy from theta
  return (u*v*sin(a)*sin(z) + u*w*cos(z) - v*v*sin(z)*cos(a) - w*w*sin(z)*cos(a))/(sqrt((u*u + v*v + w*w - pow(u*sin(z)*cos(a) + v*sin(a)*sin(z) + w*cos(z), 2))/(u*u + v*v + w*w))*pow(u*u + v*v + w*w, 1.5));
}

Numeric zeeman_magnetic_dtheta_dv(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a)
{
  // autogenerated with sympy from theta
  return (-u*u*sin(a)*sin(z) + u*v*sin(z)*cos(a) + v*w*cos(z) - w*w*sin(a)*sin(z))/(sqrt((u*u + v*v + w*w - pow(u*sin(z)*cos(a) + v*sin(a)*sin(z) + w*cos(z), 2))/(u*u + v*v + w*w))*pow(u*u + v*v + w*w, 1.5));
}

Numeric zeeman_magnetic_dtheta_dw(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a)
{
  // autogenerated with sympy from theta
  return (-u*u*cos(z) + u*w*sin(z)*cos(a) - v*v*cos(z) + v*w*sin(a)*sin(z))/(sqrt((u*u + v*v + w*w - pow(u*sin(z)*cos(a) + v*sin(a)*sin(z) + w*cos(z), 2))/(u*u + v*v + w*w))*pow(u*u + v*v + w*w, 1.5));
}

Numeric zeeman_magnetic_eta(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a)
{
  /*
               ( / -sin(a) \   [ / u \   / cos(a) sin(z) \   / u \   / u \ ]   / cos(a) sin(z) \   // / -sin(a) \   [ / u \   / cos(a) sin(z) \   / u \   / u \ ] )
     eta = atan( |  cos(a) | x [ | v | - | sin(a) sin(z) | o | v | * | v | ] o | sin(a) sin(z) |  //  |  cos(a) | o [ | v | - | sin(a) sin(z) | o | v | * | v | ] )
               ( \    0    /   [ \ w /   \        cos(z) /   \ w /   \ w / ]   \        cos(z) / //   \    0    /   [ \ w /   \        cos(z) /   \ w /   \ w / ] )
  */
  return atan2((u*cos(a)*cos(z) + v*sin(a)*cos(z) - w*sin(z))*(u*sin(z)*cos(a) + v*sin(a)*sin(z) + w*cos(z) - 1), (u*sin(a) - v*cos(a))*(u*sin(z)*cos(a) + v*sin(a)*sin(z) + w*cos(z) - 1));
}

Numeric zeeman_magnetic_deta_du(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a)
{
  // autogenerated with sympy from eta
  return (-v*cos(z) + w*sin(a)*sin(z))/(u*u*sin(a)*sin(a)*sin(z)*sin(z) - u*u*sin(z)*sin(z) + u*u - 2*u*v*sin(a)*sin(z)*sin(z)*cos(a) - 2*u*w*sin(z)*cos(a)*cos(z) - v*v*sin(a)*sin(a)*sin(z)*sin(z) + v*v - 2*v*w*sin(a)*sin(z)*cos(z) + w*w*sin(z)*sin(z));
}

Numeric zeeman_magnetic_deta_dv(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a)
{
  // autogenerated with sympy from eta
  return (u*cos(z) - w*sin(z)*cos(a))/(u*u*sin(a)*sin(a)*sin(z)*sin(z) - u*u*sin(z)*sin(z) + u*u - 2*u*v*sin(a)*sin(z)*sin(z)*cos(a) - 2*u*w*sin(z)*cos(a)*cos(z) - v*v*sin(a)*sin(a)*sin(z)*sin(z) + v*v - 2*v*w*sin(a)*sin(z)*cos(z) + w*w*sin(z)*sin(z));
}

Numeric zeeman_magnetic_deta_dw(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a)
{
  // autogenerated with sympy from eta
  return -(u*sin(a) - v*cos(a))*sin(z)/(pow(u*sin(a) - v*cos(a), 2) + pow(u*cos(a)*cos(z) + v*sin(a)*cos(z) - w*sin(z), 2));
}


Vector proj3D(const ConstVectorView a, const ConstVectorView b)
{
  return Vector({a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]});
}


Vector cross3D(const ConstVectorView a, const ConstVectorView b)
{
  return Vector({a[0] * b * a, a[1] * b * a, a[2] * b * a});
}
