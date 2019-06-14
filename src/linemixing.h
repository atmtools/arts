/* Copyright 2018, Richard Larsson
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/


// NOTE: This file contains only experimental code to test the F90-routines.
// It might evolve to replace it at some point but that is beyond present intent


#ifndef linemixing_h
#define linemixing_h

#include "rational.h"
#include "linerecord.h"
#include "absorption.h"
#include "constants.h"


template<class T> inline Numeric rotational_energy_hund_b_molecule(
  T J, T N, T J2,
  Numeric B, Numeric D, Numeric H, 
  Numeric gamma, Numeric gamma_D, Numeric gamma_H, 
  Numeric lambda, Numeric lambda_D, Numeric lambda_H)
{
  
  if(N == J and N == J2) {
    const Numeric JJ1 = J * (J + 1);
    const Numeric JJ2 = JJ1 * JJ1;
    const Numeric JJ3 = JJ1 * JJ1 * JJ1;
    return
    + (B * JJ1 - D * JJ2 + H * JJ3)
    - (gamma + gamma_D * JJ1 + gamma_H * JJ2)
    + 2/3 * (lambda + lambda_D * JJ1 + lambda_H * JJ2);
  }
  else if(N == (J - 1) and N == (J2 - 1)) {
    const Numeric JJ1 = J * (J - 1);
    const Numeric JJ2 = JJ1 * JJ1;
    const Numeric JJ3 = JJ1 * JJ1 * JJ1;
    return
    + (B * JJ1 - D * JJ2 + H * JJ3)
    + (gamma + gamma_D * JJ1 + gamma_H * JJ2) * (J - 1)
    + (lambda + lambda_D * JJ1 + lambda_H * JJ2) * (2/3 - 2*J /(2*J + 1));
  }
  else if(N == (J + 1) and N == (J2 + 1)) {
    const Numeric JJ1 = (J + 2) * (J + 1);
    const Numeric JJ2 = JJ1 * JJ1;
    const Numeric JJ3 = JJ1 * JJ1 * JJ1;
    return
    + (B * JJ1 - D * JJ2 + H * JJ3)
    - (gamma + gamma_D * JJ1 + gamma_H * JJ2) * (J + 2)
    + (lambda + lambda_D * JJ1 + lambda_H * JJ2) * (2/3 - 2*(J+1) /(2*J + 1));
    
  }
  else if((N == (J + 1) and N == (J2 - 1)) or (N == (J - 1) and N == (J2 + 1))) {
    const Numeric JJ1 = J * J + J + 1;
    const Numeric JJ2 = JJ1 * JJ1;
    return (lambda + lambda_D * JJ1 + lambda_H * JJ2) * 2 * std::sqrt(J*J+J)/(2*J + 1);
  }
  else {
    return 0;
  }
}

namespace Molecule {
  enum class Name {
    O2_66,
    CO2_626
  };
  
  namespace O2_66 {
    constexpr Numeric g_S = 2.002084;
    constexpr Numeric ge_l = 2.77e-3;
    constexpr Numeric g_r = -1.16e-4;
    
    constexpr Numeric B = 43100.44276e6;
    constexpr Numeric D = 145.1271e3;
    constexpr Numeric H = 49e-3;
    
    constexpr Numeric lambda = 59501.3438e6;
    constexpr Numeric lambda_D = 58.3680e3;
    constexpr Numeric lambda_H = 290.8e-3;
    
    constexpr Numeric gamma = -252.58634e6;
    constexpr Numeric gamma_D = -243.42;
    constexpr Numeric gamma_H = -1.46e-3;
    
    constexpr Numeric mass = 31.989830;
    
    /*! Compute Hamiltonian frequency
     
     Use dcol and drow to offset column from coupling of (N=J, N=J) to (N=J+dcol;N=J+drow).
     
     See Larsson, Lankhaar, Eriksson; Feb 2019; JQSRT for details of the matrix of these offsets
     */
    template<class T> Numeric hamiltonian_freq(T J, int dcol=0, int drow=0) {
      /*
      Matrix:
        (-1,-1), (0,-1), (1,-1)
        (-1, 0), (0, 0), (1, 0)
        (-1, 1), (0, 1), (1, 1)
      */
      return 
      rotational_energy_hund_b_molecule(
        J+dcol, J, J+drow, B, D, H, 
        gamma, gamma_D, gamma_H, 
        lambda, lambda_D, lambda_H);
    }
  };
  
  namespace CO2_626 {
    constexpr Numeric B = Conversion::kaycm2freq(0.39021);
    template<class T> constexpr Numeric hamiltonian_freq(T J) {
      return B * J * (J + 1);
    }
  }
};


class AdiabaticFactor {
public:
  enum class Type {Hartmann};
  enum class HartmannPos : Index {dc, SIZE};
  
  AdiabaticFactor(const Vector& v, Type t) : mtype(t), mdata(v) {
    bool error = false;
    switch(mtype) {
      case Type::Hartmann:
        if(Index(HartmannPos::SIZE) not_eq mdata.nelem())
          error = true;
        break;
    }
    
    if(error)
      throw std::runtime_error("Bad initializaton of BasisRate, type and size disagree...");
  }
  
  
  Numeric mol_X(const int L,
                const Numeric& B0,
                const Numeric& T,
                const Numeric& main_mass,
                const Numeric& collider_mass) const;
  
  Numeric get(const int L,
              const Numeric& B0,
              const Numeric& T,
              const Numeric& main_mass,
              const Numeric& collider_mass) const {
    switch(mtype) {
      case Type::Hartmann:
        return mol_X(L, B0, T, main_mass, collider_mass);
    }
    return 0;
  }
  
private:
  Type mtype;
  Vector mdata;
};


class BasisRate {
public:
  enum class Type {Hartmann};
  enum class HartmannPos : Index {a1, a2, a3, SIZE};
  
  BasisRate(const Vector& v, Type t) : mtype(t), mdata(v) {
    bool error = false;
    switch(mtype) {
      case Type::Hartmann:
        if(Index(HartmannPos::SIZE) not_eq mdata.nelem())
          error = true;
        break;
    }
    
    if(error)
      throw std::runtime_error("Bad initializaton of BasisRate, type and size disagree...");
  }
  
  Numeric mol_X(const int L, const Numeric& B0, const Numeric& T) const;
  
  Numeric get(const int L, const Numeric& B0, const Numeric& T) const {
    switch(mtype) {
      case Type::Hartmann:
        return mol_X(L, B0, T);
    }
    return 0;
  }
  
private:
  Type mtype;
  Vector mdata;
};


struct OffDiagonalElementOutput {Numeric ij, ji;};


namespace OffDiagonalElement {
  enum class Type {CO2_IR, O2_66_MW};
  
  OffDiagonalElementOutput CO2_IR(const LineRecord& j_line,
                                  const LineRecord& k_line,
                                  const Numeric& j_rho,
                                  const Numeric& k_rho,
                                  const BasisRate& br,
                                  const AdiabaticFactor& af,
                                  const Numeric& T,
                                  const Numeric& B0,
                                  const Numeric& main_mass,
                                  const Numeric& collider_mass);
  
  OffDiagonalElementOutput O2_66_MW(const LineRecord& line1,
                                    const LineRecord& line2,
                                    const Numeric& rho1,
                                    const Numeric& rho2,
                                    const Numeric& T,
                                    const Numeric& collider_mass);
};


Matrix hartmann_ecs_interface(const ArrayOfLineRecord& abs_lines,
                              const ArrayOfSpeciesTag& main_species,
                              const ArrayOfSpeciesTag& collider_species,
                              const Vector& collider_species_vmr,
                              const SpeciesAuxData& partition_functions,
                              const Numeric& T,
                              const Index& size);

Vector population_density_vector(const ArrayOfLineRecord& abs_lines,
                                 const SpeciesAuxData& partition_functions,
                                 const Numeric& T);

Vector dipole_vector(const ArrayOfLineRecord& abs_lines,
                     const SpeciesAuxData& partition_functions);

enum class RedPoleType {
  ElectricRoVibDipole,
  MagneticQuadrapole
};

Vector reduced_dipole_vector(const ArrayOfLineRecord& abs_lines,
                             const RedPoleType type);

Vector rosenkranz_scaling_second_order(const ArrayOfLineRecord& abs_lines,
                                       const Matrix& W,
                                       const Vector& d0);

Vector rosenkranz_shifting_second_order(const ArrayOfLineRecord& abs_lines,
                                        const Matrix& W);

Vector rosenkranz_first_order(const ArrayOfLineRecord& abs_lines,
                              const Matrix& W,
                              const Vector& d0);


struct SecondOrderLineMixingCoeffs {Numeric y0, y1;};

SecondOrderLineMixingCoeffs compute_2nd_order_lm_coeff(ConstVectorView y, ConstVectorView x, const Numeric exp, const Numeric x0);

ComplexVector equivalent_linestrengths(const Vector& population,
                                       const Vector& dipole,
                                       const Eigen::ComplexEigenSolver<Eigen::MatrixXcd>& M);

Numeric total_linestrengths(const Vector& population,
                            const Vector& dipole);

#endif // linemixing_h
