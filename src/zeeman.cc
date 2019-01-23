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

typedef struct {Numeric H, theta, eta, dH_du, dH_dv, dH_dw, dtheta_du, dtheta_dv, dtheta_dw, deta_du, deta_dv, deta_dw;} ZeemanDerived;

ZeemanDerived zeeman_internal_variables(Numeric u, Numeric v, Numeric w, Numeric z, Numeric a) noexcept
{
  // Constants evaluated once for both eta and theta
  auto cz = std::cos(z), ca = std::cos(a), sz = std::sin(z), sa = std::sin(a);
  
  /*
   *   H = ||{u, v, w}||
   */
  auto H = std::hypot(std::hypot(u, v), w);
  auto dH_du = H>0 ? u / H : 0;
  auto dH_dv = H>0 ? v / H : 0;
  auto dH_dw = H>0 ? w / H : 0;
  
  /*
   *              ( / cos(a) sin(z) \   / u \   // || / u \ || )
   *  theta = acos( | sin(a) sin(z) | o | v |  //  || | v | || )
   *              ( \        cos(z) /   \ w / //   || \ w / || )
   */
  auto x = u*sz*ca + v*sa*sz + w*cz, d = std::sqrt(1 - std::pow(x/H, 2)) * H*H*H;
  
  auto theta = H>0 ? std::acos(x / H) : std::acos(0);
  auto dtheta_du = d not_eq 0 ? (u*x - H*H*sz*ca) / d : 0;
  auto dtheta_dv = d not_eq 0 ? (v*x - H*H*sa*sz) / d : 0;
  auto dtheta_dw = d not_eq 0 ? (w*x - H*H*cz)    / d : 0;
  
  /*
   *              ( / -sin(a) \   [ / u \   / cos(a) sin(z) \   / u \   / u \ ]   / cos(a) sin(z) \   // / -sin(a) \   [ / u \   / cos(a) sin(z) \   / u \   / u \ ] )
   *    eta = atan( |  cos(a) | x [ | v | - | sin(a) sin(z) | o | v | * | v | ] o | sin(a) sin(z) |  //  |  cos(a) | o [ | v | - | sin(a) sin(z) | o | v | * | v | ] )
   *              ( \    0    /   [ \ w /   \        cos(z) /   \ w /   \ w / ]   \        cos(z) / //   \    0    /   [ \ w /   \        cos(z) /   \ w /   \ w / ] )
   */
  auto p = std::pow(u*sa - v*ca, 2) + std::pow(u*ca*cz + v*sa*cz - w*sz, 2);
  
  auto eta = std::atan2(u*ca*cz + v*sa*cz - w*sz, u*sa - v*ca);
  auto deta_du = p not_eq 0 ? (-v*cz + w*sa*sz) / p : 0;
  auto deta_dv = p not_eq 0 ? ( u*cz - w*sz*ca) / p : 0;
  auto deta_dw = p not_eq 0 ? -(u*sa - v*ca)*sz / p : 0;
  
  return {H, theta, eta,  dH_du, dH_dv, dH_dw,  dtheta_du, dtheta_dv, dtheta_dw,  deta_du, deta_dv, deta_dw};
}

ZeemanDerived zeeman_internal_variables_manual(Numeric H, Numeric theta, Numeric eta) noexcept
{
  return {H, theta, eta,  0, 0, 0,  0, 0, 0,  0, 0, 0};
}

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
                       const Vector& f_grid,
                       const Vector& rtp_vmrs, 
                       const Vector& rtp_nlte, 
                       const Vector& rtp_mag,
                       const Vector& rtp_los,
                       const Numeric& rtp_pressure,
                       const Numeric& rtp_temperature,
                       const Index& manual_tag,
                       const Numeric& H0,
                       const Numeric& theta0,
                       const Numeric& eta0)
{
  extern const Numeric DEG2RAD;
  
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
  Eigen::VectorXcd F(nf);
  Eigen::MatrixXcd dF(nf, nq);
  Eigen::VectorXcd N(nf);
  Eigen::MatrixXcd dN(nf, nq);
  Eigen::MatrixXcd data(nf, Linefunctions::ExpectedDataSize());
  const auto f_grid_eigen = MapToEigen(f_grid);
  Index start, nelem;
  
  // Magnetic field variables
  const Numeric u = rtp_mag[0];
  const Numeric v = rtp_mag[1];
  const Numeric w = rtp_mag[2];
  const Numeric z = DEG2RAD * rtp_los[0];
  const Numeric a = DEG2RAD * rtp_los[1];
  
  // Magnetic field internals and derivatives...
  const auto X = manual_tag ? zeeman_internal_variables_manual(H0, theta0, eta0) :
                              zeeman_internal_variables(u, v, w, z, a);
  
  // Polarization
  const auto polarization_scale_data        = zeeman_polarization(        X.theta, X.eta);
  const auto polarization_scale_dtheta_data = zeeman_dpolarization_dtheta(X.theta, X.eta);
  const auto polarization_scale_deta_data   = zeeman_dpolarization_deta(  X.theta, X.eta);
  
  for(Index ispecies=0; ispecies<ns; ispecies++) {
    for(const ArrayOfLineRecord& lines: zeeman_linerecord_precalc) {
      if(not lines.nelem()) continue;
      else if(lines[0].Species()      not_eq abs_species[ispecies][0].Species() or
              lines[0].Isotopologue() not_eq abs_species[ispecies][0].Isotopologue()) continue;
      
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
        
        // Select polarization for this line
        auto& polarization_scale = select_zeeman_polarization(polarization_scale_data,        line.ZeemanPolarization());
        auto& dpol_dtheta        = select_zeeman_polarization(polarization_scale_dtheta_data, line.ZeemanPolarization());
        auto& dpol_deta          = select_zeeman_polarization(polarization_scale_deta_data,   line.ZeemanPolarization());
        
        for(Index iz=0; iz<line.ZeemanEffect().nelem(); iz++) {
          Linefunctions::set_cross_section_for_single_line(F, dF, N, dN, data, start, nelem, f_grid_eigen, line,
                                                           flag_partials, flag_partials_positions, rtp_vmrs, 
                                                           rtp_nlte, rtp_pressure, rtp_temperature, dc, partial_pressure, 
                                                           isotop_ratio, X.H, ddc_dT, qt, dqt_dT, qt0,
                                                           abs_species, ispecies, iz);
          
          auto F_seg =      F.segment(start, nelem);
          auto pol_real = polarization_scale.head<4>();
          auto pol_imag = polarization_scale.tail<3>();
          auto abs = propmat_clearsky[ispecies].GetData()(0, 0, joker, joker);
          
          // Propagation matrix calculations
          MapToEigen(abs).leftCols<4>().middleRows(start, nelem).noalias() += numdens * F_seg.real() * pol_real;
          MapToEigen(abs).rightCols<3>().middleRows(start, nelem).noalias() += numdens * F_seg.imag() * pol_imag;
          for(Index j=0; j<nq; j++) {
            auto dF_seg =      dF.col(j).segment(start, nelem);
            auto dabs = dpropmat_clearsky_dx[j].GetData()(0, 0, joker, joker);
            
            if(flag_partials[flag_partials_positions[j]] == JacPropMatType::Temperature) {
              MapToEigen(dabs).leftCols<4>().middleRows(start, nelem).noalias() += numdens * dF_seg.real() * pol_real + dnumdens_dT * F_seg.real() * pol_real;
              MapToEigen(dabs).rightCols<3>().middleRows(start, nelem).noalias() += numdens * dF_seg.imag() * pol_imag + dnumdens_dT * F_seg.imag() * pol_imag;
            }
            else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticEta) {
              MapToEigen(dabs).leftCols<4>().middleRows(start, nelem).noalias() += numdens * F_seg.real() * dpol_deta.head<4>();
              MapToEigen(dabs).rightCols<3>().middleRows(start, nelem).noalias() += numdens * F_seg.imag() * dpol_deta.tail<3>();
            }
            else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticTheta) {
              MapToEigen(dabs).leftCols<4>().middleRows(start, nelem).noalias() += numdens * F_seg.real() * dpol_dtheta.head<4>();
              MapToEigen(dabs).rightCols<3>().middleRows(start, nelem).noalias() += numdens * F_seg.imag() * dpol_dtheta.tail<3>();
            }
            else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticU) {
              MapToEigen(dabs).leftCols<4>().middleRows(start, nelem).noalias() += numdens * X.dH_du * dF_seg.real() * pol_real + numdens * X.deta_du * F_seg.real() * dpol_deta.head<4>() + numdens * X.dtheta_du * F_seg.real() * dpol_dtheta.head<4>();
              MapToEigen(dabs).rightCols<3>().middleRows(start, nelem).noalias() += numdens * X.dH_du * dF_seg.imag() * pol_imag + numdens * X.deta_du * F_seg.imag() * dpol_deta.tail<3>() + numdens * X.dtheta_du * F_seg.imag() * dpol_dtheta.tail<3>();
            }
            else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticV) {
              MapToEigen(dabs).leftCols<4>().middleRows(start, nelem).noalias() += numdens * X.dH_dv * dF_seg.real() * pol_real + numdens * X.deta_dv * F_seg.real() * dpol_deta.head<4>() + numdens * X.dtheta_dv * F_seg.real() * dpol_dtheta.head<4>();
              MapToEigen(dabs).rightCols<3>().middleRows(start, nelem).noalias() += numdens * X.dH_dv * dF_seg.imag() * pol_imag + numdens * X.deta_dv * F_seg.imag() * dpol_deta.tail<3>() + numdens * X.dtheta_dv * F_seg.imag() * dpol_dtheta.tail<3>();
            }
            else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticW) {
              MapToEigen(dabs).leftCols<4>().middleRows(start, nelem).noalias() += numdens * X.dH_dw * dF_seg.real() * pol_real + numdens * X.deta_dw * F_seg.real() * dpol_deta.head<4>() + numdens * X.dtheta_dw * F_seg.real() * dpol_dtheta.head<4>();
              MapToEigen(dabs).rightCols<3>().middleRows(start, nelem).noalias() += numdens * X.dH_dw * dF_seg.imag() * pol_imag + numdens * X.deta_dw * F_seg.imag() * dpol_deta.tail<3>() + numdens * X.dtheta_dw * F_seg.imag() * dpol_dtheta.tail<3>();
            }
            else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::VMR) {
              if(flag_partials[flag_partials_positions[j]].QuantumIdentity() < line.QuantumIdentity()) {
                MapToEigen(dabs).leftCols<4>().middleRows(start, nelem).noalias() += numdens * dF_seg.real() * pol_real + dnumdens_dmvr * F_seg.real() * pol_real;
                MapToEigen(dabs).rightCols<3>().middleRows(start, nelem).noalias() += numdens * dF_seg.imag() * pol_imag + dnumdens_dmvr * F_seg.imag() * pol_imag;
              }
            }
            else {
              MapToEigen(dabs).leftCols<4>().middleRows(start, nelem).noalias() += numdens * dF_seg.real() * pol_real;
              MapToEigen(dabs).rightCols<3>().middleRows(start, nelem).noalias() += numdens * dF_seg.imag() * pol_imag;
            }
          }
        
          // Source vector calculations
          if(nn) {
            auto B     =  planck(   line.F(), rtp_temperature);
            auto dB_dT = dplanck_dt(line.F(), rtp_temperature);
            
            auto N_seg = N.segment(start, nelem);
            auto nlte_src = nlte_source[ispecies].GetData()(0, 0, joker, joker);
            
            MapToEigen(nlte_src).leftCols<4>().middleRows(start, nelem).noalias() += B * numdens * N_seg.real() * pol_real;
            for(Index j=0; j<nq; j++) {
              auto dN_seg = dN.col(j).segment(start, nelem);
              auto dnlte_dx_src = dnlte_dx_source[j].GetData()(0, 0, joker, joker);
              auto nlte_dsrc_dx = nlte_dsource_dx[j].GetData()(0, 0, joker, joker);
              
              if(flag_partials[flag_partials_positions[j]] == JacPropMatType::Temperature) {
                MapToEigen(dnlte_dx_src).leftCols<4>().middleRows(start, nelem).noalias() += B * dnumdens_dT * N_seg.real() * pol_real + B * numdens * dN_seg.real() * pol_real;
                
                MapToEigen(nlte_dsrc_dx).leftCols<4>().middleRows(start, nelem).noalias() += numdens * dB_dT * N_seg.real() * pol_real;
              }
              else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticEta)
                MapToEigen(dnlte_dx_src).leftCols<4>().middleRows(start, nelem).noalias() += B * numdens * N_seg.real() * dpol_deta.head<4>();
              else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticTheta)
                MapToEigen(dnlte_dx_src).leftCols<4>().middleRows(start, nelem).noalias() += B * numdens * N_seg.real() * dpol_dtheta.head<4>();
              else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticU)
                MapToEigen(dnlte_dx_src).leftCols<4>().middleRows(start, nelem).noalias() += B * numdens * X.dH_du * dN_seg.real() * pol_real + B * numdens * X.deta_du * N_seg.real() * dpol_deta.head<4>() + B * numdens * X.dtheta_du * N_seg.real() * dpol_dtheta.head<4>();
              else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticV)
                MapToEigen(dnlte_dx_src).leftCols<4>().middleRows(start, nelem).noalias() += B * numdens * X.dH_dv * dN_seg.real() * pol_real + B * numdens * X.deta_dv * N_seg.real() * dpol_deta.head<4>() + B * numdens * X.dtheta_dv * N_seg.real() * dpol_dtheta.head<4>();
              else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::MagneticW)
                MapToEigen(dnlte_dx_src).leftCols<4>().middleRows(start, nelem).noalias() += B * numdens * X.dH_dw * dN_seg.real() * pol_real + B * numdens * X.deta_dw * N_seg.real() * dpol_deta.head<4>() + B * numdens * X.dtheta_dw * N_seg.real() * dpol_dtheta.head<4>();
              else if(flag_partials[flag_partials_positions[j]] == JacPropMatType::VMR) {
                if(flag_partials[flag_partials_positions[j]].QuantumIdentity() < line.QuantumIdentity())
                  MapToEigen(dnlte_dx_src).leftCols<4>().middleRows(start, nelem).noalias() += B * dnumdens_dmvr * N_seg.real() * pol_real + B * numdens * dN_seg.real() * pol_real;
              }
              else
                MapToEigen(dnlte_dx_src).leftCols<4>().middleRows(start, nelem).noalias() += B * numdens * dN_seg.real() * pol_real;
            }
          }
        }
      }
    }
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
