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


void nlte_collision_factorsCalcFromCoeffs(Vector& Cij,
                                          Vector& Cji,
                                          const ArrayOfLineRecord& abs_lines,
                                          const ArrayOfGriddedField1& nlte_collision_coefficients,
                                          const ArrayOfQuantumIdentifier& nlte_collision_identifiers,
                                          const Numeric& T,
                                          const Numeric& P)
{
  if(nlte_collision_coefficients.nelem() not_eq nlte_collision_identifiers.nelem())
    throw std::runtime_error("Bad length of nlte_collision_* parameters.");
  
  extern const Numeric BOLTZMAN_CONST;
  
  const Numeric n = P / (BOLTZMAN_CONST * T);
  
  for(Index i=0; i<abs_lines.nelem(); i++) {
    const LineRecord& line = abs_lines[i];
    for(Index j=0; j<nlte_collision_coefficients.nelem(); j++) {
      if(line.QuantumIdentity() > nlte_collision_identifiers[j]) {
        const GriddedField1& gf1 = nlte_collision_coefficients[0];
        
        GridPosPoly gp;
        gridpos_poly(gp, gf1.get_numeric_grid(0), T, 1, 0.5);
        
        Vector itw(gp.idx.nelem());
        interpweights(itw,   gp);
        
        Cij[i] = interp(itw, gf1.data, gp) * n;
        
        break;
      }
    }
  }
  
  setCji(Cji, Cij, abs_lines, T, Cij.nelem());
}


void nlte_fieldRescalePopulationLevels(Tensor4& nlte_field, const Numeric& scale, const Verbosity&)
{ 
  nlte_field *= scale;
}


void nlte_fieldForSingleSpeciesNonOverlappingLines(Workspace&                      ws,
                                                   Tensor4&                        nlte_field,
                                                   const ArrayOfArrayOfSpeciesTag& abs_species,
                                                   const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                                   const ArrayOfQuantumIdentifier& nlte_levels,
                                                   const ArrayOfGriddedField1& nlte_collision_coefficients,
                                                   const ArrayOfQuantumIdentifier& nlte_collision_identifiers,
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
                                                   const Numeric&                  convergence_limit,
                                                   const Index&                    nz,
                                                   const Index&                    nf,
                                                   const Index&                    dampened,
                                                   const Index&                    iteration_limit,
                                                   const Verbosity&                verbosity)
{
  CREATE_OUT2;
  
  if(not nlte_do)
    throw std::runtime_error("Must be set to do NLTE");
  if(nlte_field.empty())
    throw std::runtime_error("Error in NLTE field, it is empty");
  
  Matrix iy;
  Tensor3 iy_transmission;
  
  const Index nlevels = nlte_levels.nelem(), np = p_grid.nelem();
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
    throw std::runtime_error("Bad number of lines... overlapping lines in nlte_levels?");
  
  // Create basic compute vectors
  const Vector Aij = createAij(lines);
  const Vector Bij = createBij(lines);
  const Vector Bji = createBji(Bij, lines);
  Vector Cij(nlines), Cji(nlines);
  
  ArrayOfIndex upper, lower;
  nlte_positions_in_statistical_equilibrium_matrix(upper, lower, lines, nlte_levels);
  const Index unique = find_first_unique_in_lower(upper, lower);
  
  // Compute arrays
  Matrix SEE(nlevels, nlevels, 0.0);
  Vector r(nlevels, 0.0), x(nlevels, 0.0);
  Numeric max_change=100;
  
  Index i = 0;
  while(i < iteration_limit and max_change > convergence_limit) {
    
    // Reset change
    max_change=0.0;
    
    //Compute radiation and transmission
    radiation_fieldCalcForSingleSpeciesNonOverlappingLines(ws, iy, iy_transmission, abs_species, abs_lines_per_species, 
                                                           nlte_field, vmr_field, t_field, z_field,
                                                           p_grid, atmosphere_dim, surface_props_data, iy_space_agenda, iy_surface_agenda,
                                                           iy_cloudbox_agenda, propmat_clearsky_agenda, water_psat_agenda, df, nz, nf, verbosity);
    
    for(Index ip = 0; ip < np; ip++) {
      r = nlte_field(joker, ip, 0, 0);
      nlte_collision_factorsCalcFromCoeffs(Cij, Cji, lines, nlte_collision_coefficients, nlte_collision_identifiers, t_field(ip, 0, 0), p_grid[ip]);
      
      if(dampened == 0)
        statistical_equilibrium_equation(SEE, Aij, Bij, Bji, Cij, Cji, iy(joker, ip), upper, lower);
      else
        dampened_statistical_equilibrium_equation(SEE, r, Aij, Bij, Bji, Cij, Cji, iy(joker, ip), iy_transmission(0, joker, ip), upper, lower);
      
      use_total_number_count_statistical_equilibrium_matrix(SEE, x, r, unique);
      solve(nlte_field(joker, ip, 0, 0), SEE, x);
      
      for(Index il=0; il<nlevels; il++) {
        max_change = max(abs(nlte_field(il, ip, 0, 0) - r[il])/r[il], max_change);
      }
    }
    i++;
  }
  
  if(i < iteration_limit)
    out2 << "Converged NLTE ratios (within convergence_limit) returned after " << i << " iterations\n";
  else
    out2 << "No convergence of NLTE ratios (within convergence_limit) returned even after " << iteration_limit << " iterations\n";
}
