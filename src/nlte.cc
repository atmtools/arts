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
  DEBUG_ONLY(const Index nlevels = A.ncols();)
  assert(Bij  .nelem() == nlines and Bji  .nelem() == nlines and 
         Cij  .nelem() == nlines and Cji  .nelem() == nlines and 
         upper.nelem() == nlines and lower.nelem() == nlines);
  assert(A.nrows() == nlevels);
  
  A = 0.0;
  for(Index iline = 0; iline < nlines; iline++) {
    const Index i = upper[iline];
    const Index j = lower[iline];
    assert(i < nlevels and j < nlevels);
    
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
  DEBUG_ONLY(const Index nlevels = x.nelem();)
  assert(Bij   .nelem() == nlines and Bji  .nelem() == nlines and 
         Cij   .nelem() == nlines and Cji  .nelem() == nlines and 
         upper .nelem() == nlines and lower.nelem() == nlines and
         Lambda.nelem() == nlines);
  assert(A.nrows() == nlevels and A.ncols() == nlevels);
  
  A = 0.0;
  for(Index iline = 0; iline < nlines; iline++) {
    const Index i = upper[iline];
    const Index j = lower[iline];
    assert(i < nlevels and j < nlevels);
    
    const Numeric Source = total_number_count * (x[i] * Aij[iline] / (x[j] * Bji[iline] - x[i] * Bij[iline]));
    
    A(j, j) -=                                      Bji[iline] * (Jij[iline] - Lambda[iline] * Source) + Cji[iline];
    A(i, i) -= Aij[iline] * (1.0 - Lambda[iline]) + Bij[iline] * (Jij[iline] - Lambda[iline] * Source) + Cij[iline];
    
    A(j, i) += Aij[iline] * (1.0 - Lambda[iline]) + Bij[iline] * (Jij[iline] - Lambda[iline] * Source) + Cij[iline];
    A(i, j) +=                                      Bji[iline] * (Jij[iline] - Lambda[iline] * Source) + Cji[iline];
  }
}


void use_total_number_count_statistical_equilibrium_matrix(MatrixView A, VectorView x, ConstVectorView r, const Index row)
{
  assert(A.ncols() == A.nrows() and A.ncols() == x.nelem() and row < A.ncols());
  A(row, joker) = 1;
  x[row] = r.sum();
}


Vector createAij(const ArrayOfLineRecord& abs_lines)
{
  // Size of problem
  const Index n = abs_lines.nelem();
  Vector Aij(n);
  
  // All must be defined
  for(Index i = 0; i < n; i++) {
    Aij[i] = abs_lines[i].A();
    if(Aij[i] <= 0.0)
      throw std::runtime_error("Undefined Einstein Coefficient");
  }
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
  assert(n == abs_lines.nelem());
  Vector Bji(n);
  
  // Base equation for single state:  B12 = B21 g2 / g1
  for(Index i = 0; i < n; i++)
    Bji[i] =  Bij[i] * abs_lines[i].G_upper() / abs_lines[i].G_lower();
  return Bji;
}


// See Pack 1979...
Vector createCijFromPressureBroadening(const ArrayOfLineRecord& abs_lines, 
                                       ConstVectorView vmrs,
                                       const ArrayOfIndex& broad_spec_locations,
                                       const Numeric& T,
                                       const Numeric& P,
                                       const Index this_species,
                                       const Index water_species)
{
  const Index n = abs_lines.nelem();
  Vector Cij(n);
  setCijFromPressureBroadening(Cij, abs_lines, vmrs, broad_spec_locations, T, P, this_species, water_species, n);
  return Cij;
}


void setCijFromPressureBroadening(VectorView Cij, 
                                  const ArrayOfLineRecord& abs_lines, 
                                  ConstVectorView vmrs,
                                  const ArrayOfIndex& broad_spec_locations,
                                  const Numeric& T,
                                  const Numeric& P,
                                  const Index this_species,
                                  const Index water_species,
                                  const Index n)
{
//FIXME: Fix this function!
  assert(n == Cij.nelem() and n == abs_lines.nelem());
  
  extern const Numeric PI, SPEED_OF_LIGHT, BOLTZMAN_CONST;
  const static Numeric c0 = 2.0 * PI * BOLTZMAN_CONST * SPEED_OF_LIGHT;
  const Numeric constant = T * c0;
  
  // Base equation for single state:  C21 = kT 2PI c G21 (where G21 is the pressure broadening at 1 Pascal)
  for(Index i = 0; i < n; i++) {
    const LineRecord& line = abs_lines[i];
    Numeric g0, g2, eta, df_0, df_2, f_VC;
    line.SetPressureBroadeningParameters(g0, g2, eta, df_0, df_2, f_VC, T, P, this_species, water_species, broad_spec_locations, vmrs);
    //Cij[i] = (g0 + 1.5*g2) / P * BOLTZMAN_CONST * T * 2 * PI * SPEED_OF_LIGHT * SPEED_OF_LIGHT / 100;
    Cij[i] = 2e12/3 * constant * g0;  // FIXME: MANY ERRORS IN THIS CODE!
  }
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
  assert(n == Cij.nelem() and n == Cij.nelem() and n == abs_lines.nelem());
  
  extern const Numeric PLANCK_CONST, BOLTZMAN_CONST;
  const static Numeric c0 = - PLANCK_CONST / BOLTZMAN_CONST;
  const Numeric constant = c0 / T;
  
  // Base equation for single state:  C12 = C21 exp(-hf / kT) g2 / g1
  for(Index i = 0; i < n; i++) 
    Cji[i] = Cij[i] * exp(constant * abs_lines[i].F()) * abs_lines[i].G_upper() / abs_lines[i].G_lower();
}


void nlte_positions_in_statistical_equilibrium_matrix(ArrayOfIndex& upper, ArrayOfIndex& lower, 
                                                      const ArrayOfLineRecord& abs_lines, 
                                                      const ArrayOfQuantumIdentifier& nlte_quantum_identifiers)
{
  const Index nl = abs_lines.nelem(), nq = nlte_quantum_identifiers.nelem();
  assert(nl > nq);
  DEBUG_ONLY(for(const auto& qi : nlte_quantum_identifiers) assert(qi.Type() == QuantumIdentifier::ENERGY_LEVEL);)
  
  upper = ArrayOfIndex(nl, -1);
  lower = ArrayOfIndex(nl, -1);
  
  for(Index il = 0; il < nl; il++) {
    for(Index iq = 0; iq < nq; iq++) {
      if(nlte_quantum_identifiers[iq].QuantumMatch()[0] > abs_lines[il].LowerQuantumNumbers())
        lower[il] = iq;
      else  if(nlte_quantum_identifiers[iq].QuantumMatch()[0] > abs_lines[il].UpperQuantumNumbers())
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
