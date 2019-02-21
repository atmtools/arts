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


void statistical_equilibrium_equation(MatrixView A,
                                      ConstVectorView Aij,
                                      ConstVectorView Bij,
                                      ConstVectorView Bji,
                                      ConstVectorView Cij,
                                      ConstVectorView Cji,
                                      ConstVectorView Jij,
                                      const ArrayOfIndex& upper,
                                      const ArrayOfIndex& lower)
{
  const Index nlines = Aij.nelem();
  
  A = 0.0;
  for(Index iline = 0; iline < nlines; iline++) {
    const Index i = upper[iline];
    const Index j = lower[iline];
    
    A(j, j) -=              Bji[iline] * Jij[iline] + Cji[iline];
    A(i, i) -= Aij[iline] + Bij[iline] * Jij[iline] + Cij[iline];
    
    A(j, i) += Aij[iline] + Bij[iline] * Jij[iline] + Cij[iline];
    A(i, j) +=              Bji[iline] * Jij[iline] + Cji[iline];
  }
}


void dampened_statistical_equilibrium_equation(MatrixView A,
                                               ConstVectorView x,
                                               ConstVectorView Aij,
                                               ConstVectorView Bij,
                                               ConstVectorView Bji,
                                               ConstVectorView Cij,
                                               ConstVectorView Cji,
                                               ConstVectorView Jij,
                                               ConstVectorView Lambda,
                                               const ArrayOfIndex& upper,
                                               const ArrayOfIndex& lower,
                                               const Numeric& total_number_count)
{
  const Index nlines = Aij.nelem();
  
  A = 0.0;
  for(Index iline = 0; iline < nlines; iline++) {
    const Index i = upper[iline];
    const Index j = lower[iline];
    
    const Numeric Source = total_number_count * (x[i] * Aij[iline] / (x[j] * Bji[iline] - x[i] * Bij[iline]));
    
    A(j, j) -=                                      Bji[iline] * (Jij[iline] - Lambda[iline] * Source) + Cji[iline];
    A(i, i) -= Aij[iline] * (1.0 - Lambda[iline]) + Bij[iline] * (Jij[iline] - Lambda[iline] * Source) + Cij[iline];
    
    A(j, i) += Aij[iline] * (1.0 - Lambda[iline]) + Bij[iline] * (Jij[iline] - Lambda[iline] * Source) + Cij[iline];
    A(i, j) +=                                      Bji[iline] * (Jij[iline] - Lambda[iline] * Source) + Cji[iline];
  }
}


void set_constant_statistical_equilibrium_matrix(MatrixView A, VectorView x, const Numeric& sem_ratio, const Index row)
{
  A(row, joker) = 1.0;
  x[row] = sem_ratio;
}


Vector createAij(const ArrayOfLineRecord& abs_lines)
{
  // Size of problem
  const Index n = abs_lines.nelem();
  Vector Aij(n);
  
  for(Index i = 0; i < n; i++)
    Aij[i] = abs_lines[i].A();
  return Aij;
}


Vector createBij(const ArrayOfLineRecord& abs_lines)
{
  extern const Numeric PLANCK_CONST, SPEED_OF_LIGHT;
  const static Numeric c0 = 2.0 * PLANCK_CONST / SPEED_OF_LIGHT / SPEED_OF_LIGHT;
  
  // Size of problem
  const Index n = abs_lines.nelem();
  Vector Bij(n);
  
  // Base equation for single state:  B21 = A21 c^2 / 2 h f^3  (nb. SI, don't use this without checking your need)
  for(Index i = 0; i < n; i++)
    Bij[i] =  abs_lines[i].A() / (c0 * abs_lines[i].F() * abs_lines[i].F() * abs_lines[i].F());
  return Bij;
}


Vector createBji(ConstVectorView Bij, const ArrayOfLineRecord& abs_lines)
{ 
  // Size of problem
  const Index n = Bij.nelem();
  Vector Bji(n);
  
  // Base equation for single state:  B12 = B21 g2 / g1
  for(Index i = 0; i < n; i++)
    Bji[i] =  Bij[i] * abs_lines[i].G_upper() / abs_lines[i].G_lower();
  return Bji;
}


Vector createCji(ConstVectorView Cij, const ArrayOfLineRecord& abs_lines, const Numeric& T)
{ 
  const Index n = abs_lines.nelem();
  Vector Cji(n);
  setCji(Cji, Cij, abs_lines, T, n);
  return Cji;
}


void setCji(VectorView Cji, ConstVectorView Cij, const ArrayOfLineRecord& abs_lines, const Numeric& T, const Index n)
{
  extern const Numeric PLANCK_CONST, BOLTZMAN_CONST;
  const static Numeric c0 = - PLANCK_CONST / BOLTZMAN_CONST;
  const Numeric constant = c0 / T;
  
  // Base equation for single state:  C12 = C21 exp(-hf / kT) g2 / g1
  for(Index i = 0; i < n; i++) 
    Cji[i] = Cij[i] * exp(constant * abs_lines[i].F()) * abs_lines[i].G_upper() / abs_lines[i].G_lower();
}


void nlte_collision_factorsCalcFromCoeffs(Vector& Cij,
                                          Vector& Cji,
                                          const ArrayOfLineRecord& abs_lines,
                                          const ArrayOfArrayOfSpeciesTag& abs_species,
                                          const ArrayOfArrayOfGriddedField1& collision_coefficients,
                                          const ArrayOfQuantumIdentifier& collision_line_identifiers,
                                          const SpeciesAuxData& isotopologue_ratios,
                                          const ConstVectorView vmr,
                                          const Numeric& T,
                                          const Numeric& P)
{
  extern const Numeric BOLTZMAN_CONST;
  
  // size of problem
  const Index nspec=abs_species.nelem();
  const Index ntrans=collision_line_identifiers.nelem();
  const Index nlines=abs_lines.nelem();
  
  // reset Cij for summing later
  Cij = 0;
  
  // For all species
  for(Index i=0; i<nspec; i++) {
    // Compute the number density noting that free_electrons will behave differently
    const Numeric numden = vmr[i] * (abs_species[i][0].SpeciesNameMain() == "free_electrons" ? 1.0 : P / (BOLTZMAN_CONST * T));
    for(Index k=0; k<nlines; k++) {
      const auto& line = abs_lines[k];
      
      const Numeric isot_ratio = isotopologue_ratios.getParam(line.Species(), line.Isotopologue())[0].data[0];
      
      for(Index j=0; j<ntrans; j++) {
        const auto& transition = collision_line_identifiers[j];
        const auto& gf1 = collision_coefficients[i][j];
        
        if(transition.In(line.QuantumIdentity())) {
          
          // Standard linear ARTS interpolation
          GridPosPoly gp;
          gridpos_poly(gp, gf1.get_numeric_grid(0), T, 1, 0.5);
          Vector itw(gp.idx.nelem());
          interpweights(itw, gp);
          
          Cij[k] += interp(itw, gf1.data, gp) * numden * isot_ratio;
          
          break; // A transition can only match one line
        }
      }
    }
  }
  
  // Compute the reverse
  setCji(Cji, Cij, abs_lines, T, Cij.nelem());
}


void nlte_positions_in_statistical_equilibrium_matrix(ArrayOfIndex& upper, ArrayOfIndex& lower, 
                                                      const ArrayOfLineRecord& abs_lines, 
                                                      const ArrayOfQuantumIdentifier& nlte_level_identifiers)
{
  const Index nl = abs_lines.nelem(), nq = nlte_level_identifiers.nelem();
  
  upper = ArrayOfIndex(nl, -1);
  lower = ArrayOfIndex(nl, -1);
  
  for(Index il = 0; il < nl; il++) {
    for(Index iq = 0; iq < nq; iq++) {
      if(nlte_level_identifiers[iq].InLower(abs_lines[il].QuantumIdentity()))
        lower[il] = iq;
      else  if(nlte_level_identifiers[iq].InUpper(abs_lines[il].QuantumIdentity()))
        upper[il] = iq;
    }
  }
  
  Index i = 0;
  for(Index il = 0; il < nl; il++)
    if(upper[il] < 0 or lower[il] < 0)
      i++;
  if(i > 1)
    throw std::runtime_error("Must set upper and lower levels completely for all but one level");
}


Index find_first_unique_in_lower(const ArrayOfIndex& upper, const ArrayOfIndex& lower) noexcept
{
  for(const Index& l : lower) {
    bool test = false;
    for(const Index& u : upper)
      if(l == u)
        test = true;
    if(not test)
      return l;
  }
  return upper.nelem() - 1;
}


void check_collision_line_identifiers(const ArrayOfQuantumIdentifier& collision_line_identifiers)
{
  if(collision_line_identifiers.nelem()) {
    const Index spec = collision_line_identifiers[0].Species();
    const Index isot = collision_line_identifiers[0].Isotopologue();
    for(const auto& x: collision_line_identifiers) {
      if(spec not_eq x.Species() or isot not_eq x.Isotopologue() or x.Type() not_eq QuantumIdentifier::TRANSITION) {
        std::ostringstream os;
        os << x << "\n" << "does not match the requirements for a line identifier\n"
           << "Your list of species is:\n" << collision_line_identifiers << "\n"
           << "This contains more than one isotopologue or it contains some non-transition type identifiers.\n"
           << "It will therefore fail in current code.  You can only input transitions, and a single isotopologue.\n";
        
        throw std::runtime_error(os.str());
      }
    }
  }
}
