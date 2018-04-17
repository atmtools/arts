/* Copyright (C) 2018
   Richard Larsson
                            
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


#include "nlte.h"
#include "absorption.h"
#include "arts.h"
#include "auto_md.h"
#include "lin_alg.h"


void nlte_fieldForSingleSpeciesNonOverlappingLines(Workspace&                      ws,
                                                   Tensor4&                        nlte_field,
                                                   const ArrayOfArrayOfSpeciesTag& abs_species,
                                                   const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                                   const ArrayOfQuantumIdentifier& nlte_quantum_identifiers,
                                                   const Agenda&                   iy_space_agenda,
                                                   const Agenda&                   iy_surface_agenda,
                                                   const Agenda&                   iy_cloudbox_agenda,
                                                   const Agenda&                   propmat_clearsky_agenda,
                                                   const Agenda&                   water_psat_agenda,
                                                   const Tensor4&                  vmr_field,
                                                   const Tensor3&                  t_field,
                                                   const Tensor3&                  z_field,
                                                   const Vector&                   p_grid,
                                                   const Index&                    atmosphere_dim,
                                                   const Tensor3&                  surface_props_data,
                                                   const Index&                    nlte_do,
                                                   const Numeric&                  df,
                                                   const Index&                    nz,
                                                   const Index&                    nf,
                                                   const Index&                    dampened,
                                                   const Verbosity&                verbosity)
{
  if(not nlte_do)
    throw std::runtime_error("Must be set to do NLTE");
  if(nlte_field.empty())
    throw std::runtime_error("Error in NLTE field, it is empty");
  
  Matrix iy;
  Tensor3 iy_transmission;
  
  const Index nlevels = nlte_quantum_identifiers.nelem(), np = p_grid.nelem(), nq = nlte_quantum_identifiers.nelem();
  if(nlevels < 5)
    throw std::runtime_error("Must have more than a four levels");
  
  if(atmosphere_dim not_eq 1)
    throw std::runtime_error("Only for 1D atmosphere");
  
  ArrayOfLineRecord lines; lines.reserve(nlevels * 2);
  for(const auto& aolr : abs_lines_per_species)
    for(const auto& lr : aolr)
      if(lr.NLTELowerIndex() >= 0)
        lines.push_back(lr);
  
  Index nlines = lines.nelem();
  if(nlevels >= nlines)
    throw std::runtime_error("Bad number of lines... overlapping lines in nlte_quantum_identifiers?");
  
  // Create basic compute vectors
  const Vector Aij = createAij(lines);
  const Vector Bij = createBij(lines);
  const Vector Bji = createBji(Bij, lines);
  Vector Cij(nlines), Cji(nlines);
  
  // 
  const Index this_species  = find_first_species_tg(abs_species, lines[0].Species());
  const Index water_species = find_first_species_tg(abs_species, species_index_from_species_name("H2O"));
  ArrayOfIndex broad_spec_locations; find_broad_spec_locations(broad_spec_locations, abs_species, this_species);
  
  ArrayOfIndex upper, lower;
  nlte_positions_in_statistical_equilibrium_matrix(upper, lower, lines, nlte_quantum_identifiers);
  const Index lower_unique = find_first_unique_in_lower(upper, lower);
  
  // Compute arrays
  Matrix SEE(nq, nq, 0.0);
  Vector r(nq, 0.0), x(nq, 0.0);
  
  // Presently loop a fixed number of times.  FIXME:  Add relative error instead
  Index i = 0;
  while(i < 20) {
    // NB.  This function should become an Agenda at some point so that iy_transmission and iy are output as required by the functions below.  
    radiation_fieldCalcForSingleSpeciesNonOverlappingLines(ws, iy, iy_transmission, abs_species, abs_lines_per_species, 
                                                           nlte_field, vmr_field, t_field, z_field,
                                                           p_grid, atmosphere_dim, surface_props_data, iy_space_agenda, iy_surface_agenda,
                                                           iy_cloudbox_agenda, propmat_clearsky_agenda, water_psat_agenda, df, nz, nf, verbosity);
    
    if(dampened == 0) {
      for(Index ip = 0; ip < np; ip++) {
        r = nlte_field(joker, ip, 0, 0);
        setCijFromPressureBroadening(Cij, lines, vmr_field(joker, ip, 0, 0), broad_spec_locations, t_field(ip, 0, 0), this_species, water_species, nlines, verbosity);
        setCji(Cji, Cij, lines, t_field(ip, 0, 0), nlines);
        statistical_equilibrium_equation(SEE, Aij, Bij, Bji, Cij, Cji, iy(joker, ip), upper, lower);
        use_total_number_count_statistical_equilibrium_matrix(SEE, x, r, lower_unique);
        solve(r, SEE, x);
        
        std::cout<<r<<"\n\nSEE "<<i+1<<" (p "<<ip+1<<"/"<<np<<"):\n"<<MapToEigen(SEE)<<"\n\nX "<<i+1<<" (p "<<ip+1<<"/"<<np<<"):\n"<<x<<"\n\n";
        
        // Assume 1D
        nlte_field(joker, ip, 0, 0) = r;
      }
    }
    else {
      for(Index ip = 0; ip < np; ip++) {
        r = nlte_field(joker, ip, 0, 0);
        setCijFromPressureBroadening(Cij, lines, vmr_field(joker, ip, 0, 0), broad_spec_locations, t_field(ip, 0, 0), this_species, water_species, nlines, verbosity);
        setCji(Cji, Cij, lines, t_field(ip, 0, 0), nlines);
        dampened_statistical_equilibrium_equation(SEE, r, Aij, Bij, Bji, Cij, Cji, iy(joker, ip), iy_transmission(0, joker, ip), upper, lower);
        use_total_number_count_statistical_equilibrium_matrix(SEE, x, r, lower_unique);
        solve(nlte_field(joker, ip, 0, 0), SEE, x);
        solve(r, SEE, x);
        
        std::cout<<r<<"\n\nSEE "<<i+1<<" (p "<<ip+1<<"/"<<np<<"):\n"<<MapToEigen(SEE)<<"\n\nX "<<i+1<<" (p "<<ip+1<<"/"<<np<<"):\n"<<x<<"\n\n";
        
        // Assume 1D
        nlte_field(joker, ip, 0, 0) = r;
      }
    }
    i++;
  }
}
