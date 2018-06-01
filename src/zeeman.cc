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

/*!
 *  Defines the phase of the propagation matrix.
 *  Read Larsson et al. (2013) for an explanation.
 * 
 *  \param  K       Out:    The rotation extinction matrix.
 *  \param  theta   In:     Angle between the magnetic field and the
 *                          propagation path. In radians.
 *  \param  eta     In:     Angle to rotate planar polarization clockwise to
 *                          fit the general coordinate system. In radians.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 *  \author Richard Larsson
 *  \date   2012-08-03
 */
void phase_matrix(MatrixView K, const Numeric& theta, const Numeric& eta, const Index& DM)
{
    assert(K.nrows() == 4 );
    assert(K.ncols() == 4 );
    
    const Numeric S2T = sin(theta)*sin(theta), CT = cos(theta), CE2 = cos(2*eta), SE2 = sin(2*eta);
    
    switch( DM )
    {
        case -1: // Transitions anti-parallel to the magnetic field
            K(0,0) =   0;  K(0,1) =           0;  K(0,2) =           0;  K(0,3) =           0;
            K(1,0) =   0;  K(1,1) =           0;  K(1,2) =      2 * CT;  K(1,3) =   S2T * SE2;
            K(2,0) =   0;  K(2,1) =    - 2 * CT;  K(2,2) =           0;  K(2,3) = - S2T * CE2;
            K(3,0) =   0;  K(3,1) = - S2T * SE2;  K(3,2) =   S2T * CE2;  K(3,3) =           0;
            break;
        case  1: // Transitions parallel to the magnetic field
            K(0,0) =   0;  K(0,1) =           0;  K(0,2) =           0;  K(0,3) =           0;
            K(1,0) =   0;  K(1,1) =           0;  K(1,2) =    - 2 * CT;  K(1,3) =   S2T * SE2;
            K(2,0) =   0;  K(2,1) =      2 * CT;  K(2,2) =           0;  K(2,3) = - S2T * CE2;
            K(3,0) =   0;  K(3,1) = - S2T * SE2;  K(3,2) =   S2T * CE2;  K(3,3) =           0;
            break;
        case  0:// Transitions perpendicular to the magnetic field
            K(0,0) =   0;  K(0,1) =           0;  K(0,2) =           0;  K(0,3) =           0;
            K(1,0) =   0;  K(1,1) =           0;  K(1,2) =           0;  K(1,3) = - S2T * SE2;
            K(2,0) =   0;  K(2,1) =           0;  K(2,2) =           0;  K(2,3) =   S2T * CE2;
            K(3,0) =   0;  K(3,1) =   S2T * SE2;  K(3,2) = - S2T * CE2;  K(3,3) =           0;
            break;
        default:
            throw std::runtime_error("Impossible Delta M to phase matrix");
    };
};


/*!
 *  Defines the attenuation of the propagation matrix.
 *  Read Larsson et al. (2013) for an explanation.
 * 
 *  \param  K       Out:    The rotation extinction matrix.
 *  \param  theta   In:     Angle between the magnetic field and the
 *                          propagation path. In radians.
 *  \param  eta     In:     Angle to rotate planar polarization clockwise to
 *                          fit the general coordinate system. In radians.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 *  \author Richard Larsson
 *  \date   2012-08-03
 */
void attenuation_matrix(MatrixView K, const Numeric& theta, const Numeric& eta, const Index& DM)
{
    assert(K.nrows() == 4 );
    assert(K.ncols() == 4 );
    
    const Numeric S2T = sin(theta)*sin(theta), CT = cos(theta), C2T = CT*CT, CE2 = cos(2*eta), SE2 = sin(2*eta);
    
    switch( DM )
    {
        case -1: // Transitions anti-parallel to the magnetic field
            K(0,0) =     1 + C2T;  K(0,1) =   S2T * CE2;  K(0,2) =   S2T * SE2;  K(0,3) =    2 * CT;
            K(1,0) =   S2T * CE2;  K(1,1) =     1 + C2T;  K(1,2) =           0;  K(1,3) =         0;
            K(2,0) =   S2T * SE2;  K(2,1) =           0;  K(2,2) =     1 + C2T;  K(2,3) =         0;
            K(3,0) =      2 * CT;  K(3,1) =           0;  K(3,2) =           0;  K(3,3) =   1 + C2T;
            break;
        case  1: // Transitions parallel to the magnetic field
            K(0,0) =     1 + C2T;  K(0,1) =   S2T * CE2;  K(0,2) =   S2T * SE2;  K(0,3) =  - 2 * CT;
            K(1,0) =   S2T * CE2;  K(1,1) =     1 + C2T;  K(1,2) =           0;  K(1,3) =         0;
            K(2,0) =   S2T * SE2;  K(2,1) =           0;  K(2,2) =     1 + C2T;  K(2,3) =         0;
            K(3,0) =    - 2 * CT;  K(3,1) =           0;  K(3,2) =           0;  K(3,3) =   1 + C2T;
            break;
        case  0: // Transitions perpendicular to the magnetic field
            K(0,0) =         S2T;  K(0,1) = - S2T * CE2;  K(0,2) = - S2T * SE2;  K(0,3) =         0;
            K(1,0) = - S2T * CE2;  K(1,1) =         S2T;  K(1,2) =           0;  K(1,3) =         0;
            K(2,0) = - S2T * SE2;  K(2,1) =           0;  K(2,2) =         S2T;  K(2,3) =         0;
            K(3,0) =           0;  K(3,1) =           0;  K(3,2) =           0;  K(3,3) =       S2T;
            break;
        default:
            throw std::runtime_error("Impossible Delta M to attenuation matrix");
    };
};


/*!
 *  Defines the derivative of the phase of the propagation matrix with regards to theta
 * 
 *  \param  dK      Out:    The rotation extinction matrix derivative with regards to theta.
 *  \param  theta   In:     Angle between the magnetic field and the
 *                          propagation path. In radians.
 *  \param  eta     In:     Angle to rotate planar polarization clockwise to
 *                          fit the general coordinate system. In radians.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 *  \author Richard Larsson
 *  \date   2015-08-14
 */
void dphase_matrix_dtheta(MatrixView dK, const Numeric& theta, const Numeric& eta, const Index& DM)
{
    assert(dK.nrows() == 4 );
    assert(dK.ncols() == 4 );
    
    const Numeric ST = sin(theta), CTST = sin(theta)*cos(theta), CE2 = cos(2*eta), SE2 = sin(2*eta);
    
    switch( DM )
    {
        case -1: // Transitions anti-parallel to the magnetic field
            dK(0,0) =   0;  dK(0,1) =                0;  dK(0,2) =                0;  dK(0,3) =                0;
            dK(1,0) =   0;  dK(1,1) =                0;  dK(1,2) =         - 2 * ST;  dK(1,3) =   2 * CTST * SE2;
            dK(2,0) =   0;  dK(2,1) =           2 * ST;  dK(2,2) =                0;  dK(2,3) = - 2 * CTST * CE2;
            dK(3,0) =   0;  dK(3,1) = - 2 * CTST * SE2;  dK(3,2) =   2 * CTST * CE2;  dK(3,3) =                0;
            break;
        case  1: // Transitions parallel to the magnetic field
            dK(0,0) =   0;  dK(0,1) =                0;  dK(0,2) =                0;  dK(0,3) =                0;
            dK(1,0) =   0;  dK(1,1) =                0;  dK(1,2) =           2 * ST;  dK(1,3) =   2 * CTST * SE2;
            dK(2,0) =   0;  dK(2,1) =        -  2 * ST;  dK(2,2) =                0;  dK(2,3) = - 2 * CTST * CE2;
            dK(3,0) =   0;  dK(3,1) = - 2 * CTST * SE2;  dK(3,2) =   2 * CTST * CE2;  dK(3,3) =                0;
            break;
        case  0:// Transitions perpendicular to the magnetic field
            dK(0,0) =   0;  dK(0,1) =                0;  dK(0,2) =                0;  dK(0,3) =                0;
            dK(1,0) =   0;  dK(1,1) =                0;  dK(1,2) =                0;  dK(1,3) = - 2 * CTST * SE2;
            dK(2,0) =   0;  dK(2,1) =                0;  dK(2,2) =                0;  dK(2,3) =   2 * CTST * CE2;
            dK(3,0) =   0;  dK(3,1) =   2 * CTST * SE2;  dK(3,2) = - 2 * CTST * CE2;  dK(3,3) =                0;
            break;
        default:
            throw std::runtime_error("Impossible Delta M to phase matrix");
    };
};


/*!
 *  Defines the derivative of the phase of the propagation matrix with regards to eta
 * 
 *  \param  dK      Out:    The rotation extinction matrix derivative with regards to eta.
 *  \param  theta   In:     Angle between the magnetic field and the
 *                          propagation path. In radians.
 *  \param  eta     In:     Angle to rotate planar polarization clockwise to
 *                          fit the general coordinate system. In radians.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 *  \author Richard Larsson
 *  \date   2015-08-14
 */
void dphase_matrix_deta(MatrixView dK, const Numeric& theta, const Numeric& eta, const Index& DM)
{
    assert(dK.nrows() == 4 );
    assert(dK.ncols() == 4 );
    
    const Numeric S2T = sin(theta)*sin(theta), CE2 = cos(2*eta), SE2 = sin(2*eta);
    
    switch( DM )
    {
        case -1: // Transitions anti-parallel to the magnetic field
            dK(0,0) =   0;  dK(0,1) =               0;  dK(0,2) =               0;  dK(0,3) =               0;
            dK(1,0) =   0;  dK(1,1) =               0;  dK(1,2) =               0;  dK(1,3) =   2 * S2T * CE2;
            dK(2,0) =   0;  dK(2,1) =               0;  dK(2,2) =               0;  dK(2,3) =   2 * S2T * SE2;
            dK(3,0) =   0;  dK(3,1) = - 2 * S2T * CE2;  dK(3,2) = - 2 * S2T * SE2;  dK(3,3) =               0;
            break;
        case  1: // Transitions parallel to the magnetic field
            dK(0,0) =   0;  dK(0,1) =               0;  dK(0,2) =               0;  dK(0,3) =               0;
            dK(1,0) =   0;  dK(1,1) =               0;  dK(1,2) =               0;  dK(1,3) =   2 * S2T * CE2;
            dK(2,0) =   0;  dK(2,1) =               0;  dK(2,2) =               0;  dK(2,3) =   2 * S2T * SE2;
            dK(3,0) =   0;  dK(3,1) = - 2 * S2T * CE2;  dK(3,2) = - 2 * S2T * SE2;  dK(3,3) =               0;
            break;
        case  0:// Transitions perpendicular to the magnetic field
            dK(0,0) =   0;  dK(0,1) =               0;  dK(0,2) =               0;  dK(0,3) =               0;
            dK(1,0) =   0;  dK(1,1) =               0;  dK(1,2) =               0;  dK(1,3) = - 2 * S2T * CE2;
            dK(2,0) =   0;  dK(2,1) =               0;  dK(2,2) =               0;  dK(2,3) = - 2 * S2T * SE2;
            dK(3,0) =   0;  dK(3,1) =   2 * S2T * CE2;  dK(3,2) =   2 * S2T * SE2;  dK(3,3) =               0;
            break;
        default:
            throw std::runtime_error("Impossible Delta M to phase matrix");
    };
};


/*!
 *  Defines the derivative of the attenuation of the propagation matrix with regards to theta.
 * 
 *  \param  dK      Out:    The rotation extinction matrix derivative with regards to theta.
 *  \param  theta   In:     Angle between the magnetic field and the
 *                          propagation path. In radians.
 *  \param  eta     In:     Angle to rotate planar polarization clockwise to
 *                          fit the general coordinate system. In radians.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 *  \author Richard Larsson
 *  \date   2015-08-14
 */
void dattenuation_matrix_dtheta(MatrixView dK, const Numeric& theta, const Numeric& eta, const Index& DM)
{
    assert(dK.nrows() == 4 );
    assert(dK.ncols() == 4 );
    
    const Numeric CTST = cos(theta)*sin(theta), ST = sin(theta), CE2 = cos(2*eta),  SE2 = sin(2*eta);
    
    switch( DM )
    {
        case -1: // Transitions anti-parallel to the magnetic field
            dK(0,0) =       - 2 * CTST;  dK(0,1) =   2 * CTST * CE2;  dK(0,2) =   2 * CTST * SE2;  dK(0,3) =   - 2 * ST;
            dK(1,0) =   2 * CTST * CE2;  dK(1,1) =       - 2 * CTST;  dK(1,2) =                0;  dK(1,3) =          0;
            dK(2,0) =   2 * CTST * SE2;  dK(2,1) =                0;  dK(2,2) =       - 2 * CTST;  dK(2,3) =          0;
            dK(3,0) =         - 2 * ST;  dK(3,1) =                0;  dK(3,2) =                0;  dK(3,3) = - 2 * CTST;
            break;
        case  1: // Transitions parallel to the magnetic field
            dK(0,0) =       - 2 * CTST;  dK(0,1) =   2 * CTST * CE2;  dK(0,2) =   2 * CTST * SE2;  dK(0,3) =     2 * ST;
            dK(1,0) =   2 * CTST * CE2;  dK(1,1) =       - 2 * CTST;  dK(1,2) =                0;  dK(1,3) =          0;
            dK(2,0) =   2 * CTST * SE2;  dK(2,1) =                0;  dK(2,2) =       - 2 * CTST;  dK(2,3) =          0;
            dK(3,0) =           2 * ST;  dK(3,1) =                0;  dK(3,2) =                0;  dK(3,3) = - 2 * CTST;
            break;
        case  0: // Transitions perpendicular to the magnetic field
            dK(0,0) =         2 * CTST;  dK(0,1) = - 2 * CTST * CE2;  dK(0,2) = - 2 * CTST * SE2;  dK(0,3) =          0;
            dK(1,0) = - 2 * CTST * CE2;  dK(1,1) =         2 * CTST;  dK(1,2) =                0;  dK(1,3) =          0;
            dK(2,0) = - 2 * CTST * SE2;  dK(2,1) =                0;  dK(2,2) =         2 * CTST;  dK(2,3) =          0;
            dK(3,0) =                0;  dK(3,1) =                0;  dK(3,2) =                0;  dK(3,3) =   2 * CTST;
            break;
        default:
            throw std::runtime_error("Impossible Delta M to attenuation matrix");
    };
};


/*!
 *  Defines the derivative of the attenuation of the propagation matrix with regards to eta.
 * 
 *  \param  dK      Out:    The rotation extinction matrix derivative with regards to eta.
 *  \param  theta   In:     Angle between the magnetic field and the
 *                          propagation path. In radians.
 *  \param  eta     In:     Angle to rotate planar polarization clockwise to
 *                          fit the general coordinate system. In radians.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 *  \author Richard Larsson
 *  \date   2015-08-14
 */
void dattenuation_matrix_deta(MatrixView dK, const Numeric& theta, const Numeric& eta, const Index& DM)
{
    assert(dK.nrows() == 4 );
    assert(dK.ncols() == 4 );
    
    const Numeric SE2 = sin(2*eta), CE2 = cos(2*eta), ST2 = sin(theta)*sin(theta);
    
    switch( DM )
    {
        case -1: // Transitions anti-parallel to the magnetic field
            dK(0,0) =               0;  dK(0,1) = - 2 * ST2 * SE2;  dK(0,2) =   2 * ST2 * CE2;  dK(0,3) =   0;
            dK(1,0) = - 2 * ST2 * SE2;  dK(1,1) =               0;  dK(1,2) =               0;  dK(1,3) =   0;
            dK(2,0) =   2 * ST2 * CE2;  dK(2,1) =               0;  dK(2,2) =               0;  dK(2,3) =   0;
            dK(3,0) =               0;  dK(3,1) =               0;  dK(3,2) =               0;  dK(3,3) =   0;
            break;
        case  1: // Transitions parallel to the magnetic field
            dK(0,0) =               0;  dK(0,1) = - 2 * ST2 * SE2;  dK(0,2) =   2 * ST2 * CE2;  dK(0,3) =   0;
            dK(1,0) = - 2 * ST2 * SE2;  dK(1,1) =               0;  dK(1,2) =               0;  dK(1,3) =   0;
            dK(2,0) =   2 * ST2 * CE2;  dK(2,1) =               0;  dK(2,2) =               0;  dK(2,3) =   0;
            dK(3,0) =               0;  dK(3,1) =               0;  dK(3,2) =               0;  dK(3,3) =   0;
            break;
        case  0: // Transitions perpendicular to the magnetic field
            dK(0,0) =               0;  dK(0,1) =   2 * ST2 * SE2;  dK(0,2) = - 2 * ST2 * CE2;  dK(0,3) =   0;
            dK(1,0) =   2 * ST2 * SE2;  dK(1,1) =               0;  dK(1,2) =               0;  dK(1,3) =   0;
            dK(2,0) = - 2 * ST2 * CE2;  dK(2,1) =               0;  dK(2,2) =               0;  dK(2,3) =   0;
            dK(3,0) =               0;  dK(3,1) =               0;  dK(3,2) =               0;  dK(3,3) =   0;
            break;
        default:
            throw std::runtime_error("Impossible Delta M to attenuation matrix");
    };
};


/*!
 * Return the relative strength of the split Zeeman line parts as found in
 * Berdyugina and Solnaki (2002). Note that this is the same as the general case
 * of Schadee (1978).
 * 
 * \param  __U__   Void.
 * \param  m       In:     Secondary rotational quantum number.
 * \param  j       In:     Spin-Orbit Coupling number.
 * \param  DJ      In:     Change in the main rotational quantum number.
 * \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 * \author Richard Larsson
 * \date   2012-10-26
 */
Rational relative_strength(const Rational& M, const Rational& J, const Index& dj, const Index& dm) {
  switch ( dj ) {
    case -1:
      switch ( dm ) {
        case -1: return (Rational(3, 8)*(J+M)*(J-1+M) / (J*(2*J-1)*(2*J+1))); // Transitions anti-parallel to the magnetic field
        case  0: return (Rational(3, 2)*(J*J-M*M)     / (J*(2*J-1)*(2*J+1))); // Transitions perpendicular to the magnetic field
        case +1: return (Rational(3, 8)*(J-M)*(J-1-M) / (J*(2*J-1)*(2*J+1))); // Transitions parallel to the magnetic field
      }
    case  0:
      switch ( dm ) {
        case -1: return (Rational(3, 8)*(J+M)*(J+1-M) / (J*(J+1)*(2*J+1))); // Transitions anti-parallel to the magnetic field
        case  0: return (Rational(3, 2)*M*M           / (J*(J+1)*(2*J+1))); // Transitions perpendicular to the magnetic field
        case +1: return (Rational(3, 8)*(J-M)*(J+1+M) / (J*(J+1)*(2*J+1))); // Transitions parallel to the magnetic field
      }
    case +1:
      switch ( dm ) {
        case -1: return (Rational(3, 8)*(J+1-M)*(J+2-M)   / ((J+1)*(2*J+1)*(2*J+3))); // Transitions anti-parallel to the magnetic field
        case  0: return (Rational(3, 2)*((J+1)*(J+1)-M*M) / ((J+1)*(2*J+1)*(2*J+3))); // Transitions perpendicular to the magnetic field
        case +1: return (Rational(3, 8)*(J+1+M)*(J+2+M)   / ((J+1)*(2*J+1)*(2*J+3))); // Transitions parallel to the magnetic field
      }
    default: throw std::runtime_error("Something is extremely wrong.");
  }
}


void xsec_species_line_mixing_wrapper_with_zeeman(  ArrayOfPropagationMatrix& propmat_clearsky, 
                                                    ArrayOfStokesVector& nlte_source,
                                                    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                                                    ArrayOfStokesVector& dnlte_dx_source,
                                                    ArrayOfStokesVector& nlte_dsource_dx,
                                                    const ArrayOfArrayOfSpeciesTag& abs_species, 
                                                    const ArrayOfRetrievalQuantity& flag_partials,
                                                    const ArrayOfIndex& flag_partials_positions,
                                                    const Index& abs_lineshape_ls, 
                                                    const Index& abs_lineshape_lsn, 
                                                    const Numeric& abs_lineshape_cutoff, 
                                                    const ArrayOfLineRecord& lr, 
                                                    const Vector&  planck_BT,
                                                    const Matrix&  dplanck_BT,
                                                    const SpeciesAuxData& isotopologue_ratios, 
                                                    const SpeciesAuxData& partition_functions,
                                                    const Matrix& abs_nlte, 
                                                    const Matrix& abs_vmrs, 
                                                    const Vector& abs_p,
                                                    const Vector& abs_t, 
                                                    const Vector& f_grid,
                                                    const Vector& rtp_mag,
                                                    const Vector& r_path_los,
                                                    const Numeric& lm_p_lim,
                                                    const Numeric& theta, 
                                                    const Numeric& eta, 
                                                    const Numeric& H_mag, 
                                                    const Index& DM, 
                                                    const Index& this_species,
                                                    const Verbosity& verbosity )
                                                
                                                
{
  const bool do_src  =  !nlte_source.empty();
  const Index nq = flag_partials_positions.nelem();
  const Index nf = f_grid.nelem();
  const Numeric n = abs_vmrs(this_species, 0)*number_density( abs_p[0],abs_t[0]);
  const Numeric dn_dT = abs_vmrs(this_species, 0)*dnumber_density_dt( abs_p[0],abs_t[0]);
  
  // Setting up variables
  Matrix attenuation(nf, 1),phase(nf, 1),source(do_src?nf:0,do_src?1:0),
  attenuation_du(nf, 1),phase_du(nf, 1),attenuation_dv(nf, 1),phase_dv(nf, 1),attenuation_dw(nf, 1),phase_dw(nf, 1),
  source_du(do_src?nf:0,do_src?1:0),source_dv(do_src?nf:0,do_src?1:0),source_dw(do_src?nf:0,do_src?1:0);
  
  
  ArrayOfMatrix partial_attenuation(nq, attenuation), partial_phase(nq, phase), partial_source(do_src?nq:0, source);
  
  // FIXME: Have to perturb for magnetic u, v, and w since it is too complicated otherwise
  Numeric dB=0.0;
  if(nq)
    dB = magnetic_field_perturbation(flag_partials);
  
  Numeric H_dummy, deta_du, dtheta_du, deta_dv, dtheta_dv, deta_dw, dtheta_dw;
  Vector dmag = rtp_mag;
  
  bool do_u=false, do_v=false, do_w=false;
  for(const auto& rq : flag_partials) {
    if(rq == JacPropMatType::MagneticU) {
      dmag[0]+=dB;
      set_magnetic_parameters(H_dummy,deta_du,dtheta_du,0,0,0,0,dmag,r_path_los);
      deta_du -= eta;
      deta_du /= dB;
      dtheta_du -= theta;
      dtheta_du /= dB;
      dmag[0]-=dB;
      
      do_u = true;
    }
    else if(rq == JacPropMatType::MagneticV) {
      dmag[1]+=dB;
      set_magnetic_parameters(H_dummy,deta_dv,dtheta_dv,0,0,0,0,dmag,r_path_los);
      deta_dv -= eta;
      deta_dv /= dB;
      dtheta_dv -= theta;
      dtheta_dv /= dB;
      dmag[1]-=dB;
      
      do_v = true;
    }
    else if(rq == JacPropMatType::MagneticW) {
      dmag[2]+=dB;
      set_magnetic_parameters(H_dummy,deta_dw,dtheta_dw,0,0,0,0,dmag,r_path_los);
      deta_dw -= eta;
      deta_dw /= dB;
      dtheta_dw -= theta;
      dtheta_dw /= dB;
      dmag[2]-=dB;
      
      do_w = true;
    }
  }
  // JACOBIAN SETUP END
  
  attenuation(joker, 0) = 0.;
  phase(joker, 0) = 0.;
  if(do_u)
  {
    attenuation_du(joker,0)=0.;
    phase_du(joker,0)=0.;
  }
  if(do_v)
  {
    attenuation_dv(joker,0)=0.;
    phase_dv(joker,0)=0.;
  }
  if(do_w)
  {
    attenuation_dw(joker,0)=0.;
    phase_dw(joker,0)=0.;
  }
  if(do_src)
    source(joker, 0)=0.;
  for(Index iq = 0; iq < nq; iq++)
  {
    partial_attenuation[iq](joker,0)=0.;
    partial_phase[iq](joker,0)=0.;
    if(do_src)
      partial_source[iq](joker,0)=0.;
  }
  
  xsec_species_line_mixing_wrapper(   attenuation,         source,         phase, 
                                      partial_attenuation, partial_source, partial_phase, flag_partials, flag_partials_positions,
                                      f_grid, abs_p, abs_t, abs_nlte, abs_vmrs, abs_species, 
                                      this_species, lr, H_mag,
                                      abs_lineshape_ls,abs_lineshape_lsn,lm_p_lim,abs_lineshape_cutoff,
                                      isotopologue_ratios, partition_functions, verbosity ); // Now in cross section
  if(do_u)
  {
    dmag[0]+=dB;
    
    xsec_species_line_mixing_wrapper(         attenuation_du,         source_du,         phase_du, 
                                              partial_attenuation, partial_source, partial_phase, ArrayOfRetrievalQuantity(0), ArrayOfIndex(0),
                                              f_grid, abs_p, abs_t, abs_nlte, abs_vmrs, abs_species, 
                                              this_species, lr, sqrt(dmag*dmag),
                                              abs_lineshape_ls,abs_lineshape_lsn,lm_p_lim,abs_lineshape_cutoff,
                                              isotopologue_ratios, partition_functions, verbosity ); // Now in cross section
    dmag[0]-=dB;
  }
  if(do_v)
  {
    dmag[1]+=dB;
    
    xsec_species_line_mixing_wrapper(         attenuation_dv,         source_dv,         phase_dv, 
                                              partial_attenuation, partial_source, partial_phase, ArrayOfRetrievalQuantity(0), ArrayOfIndex(0),
                                              f_grid, abs_p, abs_t, abs_nlte, abs_vmrs, abs_species, 
                                              this_species, lr, sqrt(dmag*dmag),
                                              abs_lineshape_ls,abs_lineshape_lsn,lm_p_lim,abs_lineshape_cutoff,
                                              isotopologue_ratios, partition_functions, verbosity ); // Now in cross section
    dmag[1]-=dB;
  }
  if(do_w)
  {
    dmag[2]+=dB;
    
    xsec_species_line_mixing_wrapper(         attenuation_dw,         source_dw,         phase_dw, 
                                              partial_attenuation, partial_source, partial_phase, ArrayOfRetrievalQuantity(0), ArrayOfIndex(0),
                                              f_grid, abs_p, abs_t, abs_nlte, abs_vmrs, abs_species, 
                                              this_species, lr, sqrt(dmag*dmag),
                                              abs_lineshape_ls,abs_lineshape_lsn,lm_p_lim,abs_lineshape_cutoff,
                                              isotopologue_ratios, partition_functions, verbosity ); // Now in cross section
    dmag[2]-=dB;
  }
  
  if(DM == 0)
  {
    propmat_clearsky[this_species].AddZeemanPiComponent(attenuation(joker, 0), phase(joker, 0), n, theta*DEG2RAD, 
                                                        eta*DEG2RAD);
    
    if(do_src)
      nlte_source[this_species].AddZeemanPiComponent(source(joker,0), phase(joker, 0), n, theta*DEG2RAD, 
                                                     eta*DEG2RAD, planck_BT);
    
    for(Index iq = 0; iq < nq; iq++)
    {
      if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticU)
      {
        attenuation_du -= attenuation;
        attenuation_du /= dB;
        phase_du -= phase;
        phase_du /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanPiComponentDerivative(attenuation(joker, 0), attenuation_du(joker, 0), 
                                                                phase(joker, 0), phase_du(joker, 0), n, theta*DEG2RAD, 
                                                                dtheta_du*DEG2RAD, eta*DEG2RAD, deta_du*DEG2RAD);
        
        if(do_src)
        {
          source_du -= source;
          source_du /= dB;
          dnlte_dx_source[iq].AddZeemanPiComponentDerivative(source(joker, 0), source_du(joker, 0), phase(joker, 0), 
                                                             phase_du(joker, 0), n, theta*DEG2RAD, dtheta_du*DEG2RAD, 
                                                             eta*DEG2RAD, deta_du*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticV)
      {
        attenuation_dv -= attenuation;
        attenuation_dv /= dB;
        phase_dv -= phase;
        phase_dv /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanPiComponentDerivative(attenuation(joker, 0), attenuation_dv(joker, 0), 
                                                                phase(joker, 0), phase_dv(joker, 0), n, theta*DEG2RAD, 
                                                                dtheta_dv*DEG2RAD, eta*DEG2RAD, deta_dv*DEG2RAD);
        
        if(do_src)
        {
          source_dv -= source;
          source_dv /= dB;
          dnlte_dx_source[iq].AddZeemanPiComponentDerivative(source(joker, 0), source_dv(joker, 0), phase(joker, 0),
                                                             phase_dv(joker, 0), n, theta*DEG2RAD, dtheta_dv*DEG2RAD, 
                                                             eta*DEG2RAD, deta_dv*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticW)
      {
        attenuation_dw -= attenuation;
        attenuation_dw /= dB;
        phase_dw -= phase;
        phase_dw /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanPiComponentDerivative(attenuation(joker, 0), attenuation_dw(joker, 0), 
                                                                phase(joker, 0), phase_dw(joker, 0), n, theta*DEG2RAD,
                                                                dtheta_dw*DEG2RAD, eta*DEG2RAD, deta_dw*DEG2RAD);
        
        if(do_src)
        {
          source_dw -= source;
          source_dw /= dB;
          dnlte_dx_source[iq].AddZeemanPiComponentDerivative(source(joker, 0), source_dw(joker, 0), phase(joker, 0),
                                                             phase_dw(joker, 0), n, theta*DEG2RAD, dtheta_dw*DEG2RAD,
                                                             eta*DEG2RAD, deta_dw*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticTheta)
      {
        dpropmat_clearsky_dx[iq].AddZeemanPiComponentThetaDerivative(attenuation(joker, 0), phase(joker, 0), n, 
                                                                     theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanPiComponentThetaDerivative(source(joker, 0), phase(joker, 0), n, theta*DEG2RAD,
                                                                  eta*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticEta)
      {
        dpropmat_clearsky_dx[iq].AddZeemanPiComponentEtaDerivative(attenuation(joker, 0), phase(joker, 0), n,
                                                                   theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dpropmat_clearsky_dx[iq].AddZeemanPiComponentEtaDerivative(source(joker, 0), phase(joker, 0), n, 
                                                                     theta*DEG2RAD, eta*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::Temperature)
      {
        dpropmat_clearsky_dx[iq].AddZeemanPiComponent(attenuation(joker, 0), phase(joker, 0), dn_dT, theta*DEG2RAD, 
                                                      eta*DEG2RAD);
        dpropmat_clearsky_dx[iq].AddZeemanPiComponent(partial_attenuation[iq](joker, 0), partial_phase[iq](joker, 0), n, 
                                                      theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanPiComponent(source(joker, 0), phase(joker, 0), dn_dT, theta*DEG2RAD, eta*DEG2RAD, 
                                                   planck_BT);
          dnlte_dx_source[iq].AddZeemanPiComponent(partial_source[iq](joker, 0), phase(joker, 0), n, theta*DEG2RAD, 
                                                   eta*DEG2RAD, planck_BT);
          
          nlte_dsource_dx[iq].AddZeemanPiComponent(source(joker, 0), phase(joker, 0), n, theta*DEG2RAD, eta*DEG2RAD, 
                                                   dplanck_BT(0, joker));
        }
      }
      else if(is_frequency_parameter(flag_partials[flag_partials_positions[iq]]))
      {
        dpropmat_clearsky_dx[iq].AddZeemanPiComponent(partial_attenuation[iq](joker, 0), partial_phase[iq](joker, 0), n, 
                                                      theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanPiComponent(partial_source[iq](joker, 0), phase(joker, 0), n, theta*DEG2RAD,
                                                   eta*DEG2RAD, planck_BT);
          
          nlte_dsource_dx[iq].AddZeemanPiComponent(source(joker, 0), phase(joker, 0), n, theta*DEG2RAD, eta*DEG2RAD,
                                                   dplanck_BT(1, joker));
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::VMR)
      {  
        if(species_match(flag_partials[flag_partials_positions[iq]], abs_species[this_species])) {
          dpropmat_clearsky_dx[iq].AddZeemanPiComponent(attenuation(joker, 0), phase(joker, 0), 
                                                        n/abs_vmrs(this_species, 0), theta*DEG2RAD, eta*DEG2RAD);
          
          if(do_src) {
            dnlte_dx_source[iq].AddZeemanPiComponent(source(joker, 0), phase(joker, 0), n/abs_vmrs(this_species, 0),
                                                    theta*DEG2RAD, eta*DEG2RAD, planck_BT);
          }
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] not_eq JacPropMatType::NotPropagationMatrixType)
      {
        dpropmat_clearsky_dx[iq].AddZeemanPiComponent(partial_attenuation[iq](joker, 0), partial_phase[iq](joker, 0), n,
                                                      theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanPiComponent(partial_source[iq](joker, 0), phase(joker, 0), n, theta*DEG2RAD,
                                                   eta*DEG2RAD, planck_BT);
        }
      }
    }
  }
  else if(DM == -1)
  {
    propmat_clearsky[this_species].AddZeemanSigmaMinusComponent(attenuation(joker, 0), phase(joker, 0), n,
                                                                theta*DEG2RAD, eta*DEG2RAD);
    
    if(do_src)
      nlte_source[this_species].AddZeemanSigmaMinusComponent(source(joker,0), phase(joker, 0), n, theta*DEG2RAD,
                                                             eta*DEG2RAD, planck_BT);
    
    for(Index iq = 0; iq < nq; iq++)
    {
      if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticU)
      {
        attenuation_du -= attenuation;
        attenuation_du /= dB;
        phase_du -= phase;
        phase_du /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponentDerivative(attenuation(joker, 0), attenuation_du(joker, 0),
                                                                        phase(joker, 0), phase_du(joker, 0), n, 
                                                                        theta*DEG2RAD, dtheta_du*DEG2RAD, eta*DEG2RAD,
                                                                        deta_du*DEG2RAD);
        
        if(do_src)
        {
          source_du -= source;
          source_du /= dB;
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponentDerivative(source(joker, 0), source_du(joker, 0),
                                                                     phase(joker, 0), phase_du(joker, 0), n, 
                                                                     theta*DEG2RAD, dtheta_du*DEG2RAD, eta*DEG2RAD,
                                                                     deta_du*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticV)
      {
        attenuation_dv -= attenuation;
        attenuation_dv /= dB;
        phase_dv -= phase;
        phase_dv /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponentDerivative(attenuation(joker, 0), attenuation_dv(joker, 0),
                                                                        phase(joker, 0), phase_dv(joker, 0), n,
                                                                        theta*DEG2RAD, dtheta_dv*DEG2RAD, eta*DEG2RAD,
                                                                        deta_dv*DEG2RAD);
        
        if(do_src)
        {
          source_dv -= source;
          source_dv /= dB;
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponentDerivative(source(joker, 0), source_dv(joker, 0), 
                                                                     phase(joker, 0), phase_dv(joker, 0), n, 
                                                                     theta*DEG2RAD, dtheta_dv*DEG2RAD, eta*DEG2RAD,
                                                                     deta_dv*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticW)
      {
        attenuation_dw -= attenuation;
        attenuation_dw /= dB;
        phase_dw -= phase;
        phase_dw /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponentDerivative(attenuation(joker, 0), attenuation_dw(joker, 0),
                                                                        phase(joker, 0), phase_dw(joker, 0), n,
                                                                        theta*DEG2RAD, dtheta_dw*DEG2RAD, eta*DEG2RAD,
                                                                        deta_dw*DEG2RAD);
        
        if(do_src)
        {
          source_dw -= source;
          source_dw /= dB;
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponentDerivative(source(joker, 0), source_dw(joker, 0),
                                                                     phase(joker, 0), phase_dw(joker, 0), n,
                                                                     theta*DEG2RAD, dtheta_dw*DEG2RAD, eta*DEG2RAD,
                                                                     deta_dw*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticTheta)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponentThetaDerivative(attenuation(joker, 0), phase(joker, 0), n,
                                                                             theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponentThetaDerivative(source(joker, 0), phase(joker, 0), n,
                                                                          theta*DEG2RAD, eta*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticEta)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponentEtaDerivative(attenuation(joker, 0), phase(joker, 0), n,
                                                                           theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponentEtaDerivative(source(joker, 0), phase(joker, 0), n,
                                                                             theta*DEG2RAD, eta*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::Temperature)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponent(attenuation(joker, 0), phase(joker, 0), dn_dT,
                                                              theta*DEG2RAD, eta*DEG2RAD);
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponent(partial_attenuation[iq](joker, 0), 
                                                              partial_phase[iq](joker, 0), n, theta*DEG2RAD, 
                                                              eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponent(source(joker, 0), phase(joker, 0), dn_dT, theta*DEG2RAD,
                                                           eta*DEG2RAD, planck_BT);
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponent(partial_source[iq](joker, 0), phase(joker, 0), n,
                                                           theta*DEG2RAD, eta*DEG2RAD, planck_BT);
          
          nlte_dsource_dx[iq].AddZeemanSigmaMinusComponent(source(joker, 0), phase(joker, 0), n, theta*DEG2RAD,
                                                           eta*DEG2RAD, dplanck_BT(0, joker));
        }
      }
      else if(is_frequency_parameter(flag_partials[flag_partials_positions[iq]]))
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponent(partial_attenuation[iq](joker, 0), 
                                                              partial_phase[iq](joker, 0), n, theta*DEG2RAD,
                                                              eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponent(partial_source[iq](joker, 0), phase(joker, 0), n,
                                                           theta*DEG2RAD, eta*DEG2RAD, planck_BT);
          
          nlte_dsource_dx[iq].AddZeemanSigmaMinusComponent(source(joker, 0), phase(joker, 0), n, theta*DEG2RAD,
                                                           eta*DEG2RAD, dplanck_BT(1, joker));
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::VMR)
      {  
        if(species_match(flag_partials[flag_partials_positions[iq]], abs_species[this_species])) {
          dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponent(attenuation(joker, 0), phase(joker, 0),
                                                                n/abs_vmrs(this_species, 0), theta*DEG2RAD, eta*DEG2RAD);
          
          if(do_src) {
            dnlte_dx_source[iq].AddZeemanSigmaMinusComponent(source(joker, 0), phase(joker, 0),
                                                            n/abs_vmrs(this_species, 0), theta*DEG2RAD, eta*DEG2RAD,
                                                            planck_BT);
          }
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] not_eq JacPropMatType::NotPropagationMatrixType)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaMinusComponent(partial_attenuation[iq](joker, 0),
                                                              partial_phase[iq](joker, 0), n, theta*DEG2RAD,
                                                              eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaMinusComponent(partial_source[iq](joker, 0), phase(joker, 0), n,
                                                           theta*DEG2RAD, eta*DEG2RAD, planck_BT);
        }
      }
    }
  }
  else if(DM == 1)
  {
    propmat_clearsky[this_species].AddZeemanSigmaPlusComponent(attenuation(joker, 0), phase(joker, 0), n,
                                                               theta*DEG2RAD, eta*DEG2RAD);
    
    if(do_src)
      nlte_source[this_species].AddZeemanSigmaPlusComponent(source(joker,0), phase(joker, 0), n, theta*DEG2RAD,
                                                            eta*DEG2RAD, planck_BT);
    
    for(Index iq = 0; iq < nq; iq++)
    {
      if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticU)
      {
        attenuation_du -= attenuation;
        attenuation_du /= dB;
        phase_du -= phase;
        phase_du /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponentDerivative(attenuation(joker, 0), attenuation_du(joker, 0),
                                                                       phase(joker, 0), phase_du(joker, 0), n, 
                                                                       theta*DEG2RAD, dtheta_du*DEG2RAD, eta*DEG2RAD, 
                                                                       deta_du*DEG2RAD);
        
        if(do_src)
        {
          source_du -= source;
          source_du /= dB;
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponentDerivative(source(joker, 0), source_du(joker, 0), 
                                                                    phase(joker, 0), phase_du(joker, 0), n, 
                                                                    theta*DEG2RAD, dtheta_du*DEG2RAD, eta*DEG2RAD, 
                                                                    deta_du*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticV)
      {
        attenuation_dv -= attenuation;
        attenuation_dv /= dB;
        phase_dv -= phase;
        phase_dv /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponentDerivative(attenuation(joker, 0), attenuation_dv(joker, 0),
                                                                       phase(joker, 0), phase_dv(joker, 0), n,
                                                                       theta*DEG2RAD, dtheta_dv*DEG2RAD, eta*DEG2RAD,
                                                                       deta_dv*DEG2RAD);
        
        if(do_src)
        {
          source_dv -= source;
          source_dv /= dB;
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponentDerivative(source(joker, 0), source_dv(joker, 0),
                                                                    phase(joker, 0), phase_dv(joker, 0), n,
                                                                    theta*DEG2RAD, dtheta_dv*DEG2RAD, eta*DEG2RAD,
                                                                    deta_dv*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticW)
      {
        attenuation_dw -= attenuation;
        attenuation_dw /= dB;
        phase_dw -= phase;
        phase_dw /= dB;
        
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponentDerivative(attenuation(joker, 0), attenuation_dw(joker, 0),
                                                                       phase(joker, 0), phase_dw(joker, 0), n, 
                                                                       theta*DEG2RAD, dtheta_dw*DEG2RAD, eta*DEG2RAD,
                                                                       deta_dw*DEG2RAD);
        
        if(do_src)
        {
          source_dw -= source;
          source_dw /= dB;
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponentDerivative(source(joker, 0), source_dw(joker, 0), 
                                                                    phase(joker, 0), phase_dw(joker, 0), n,
                                                                    theta*DEG2RAD, dtheta_dw*DEG2RAD, eta*DEG2RAD,
                                                                    deta_dw*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticTheta)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponentThetaDerivative(attenuation(joker, 0), phase(joker, 0), n,
                                                                            theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponentThetaDerivative(source(joker, 0), phase(joker, 0), n,
                                                                         theta*DEG2RAD, eta*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::MagneticEta)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponentEtaDerivative(attenuation(joker, 0), phase(joker, 0), n,
                                                                          theta*DEG2RAD, eta*DEG2RAD);
        
        if(do_src)
        {
          dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponentEtaDerivative(source(joker, 0), phase(joker, 0), n,
                                                                            theta*DEG2RAD, eta*DEG2RAD, planck_BT);
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::Temperature)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponent(attenuation(joker, 0), phase(joker, 0), dn_dT,
                                                             theta*DEG2RAD, eta*DEG2RAD);
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponent(partial_attenuation[iq](joker, 0),
                                                             partial_phase[iq](joker, 0), n, theta*DEG2RAD,
                                                             eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponent(source(joker, 0), phase(joker, 0), dn_dT, theta*DEG2RAD,
                                                          eta*DEG2RAD, planck_BT);
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponent(partial_source[iq](joker, 0), phase(joker, 0), n,
                                                          theta*DEG2RAD, eta*DEG2RAD, planck_BT);
          
          nlte_dsource_dx[iq].AddZeemanSigmaPlusComponent(source(joker, 0), phase(joker, 0), n, theta*DEG2RAD,
                                                          eta*DEG2RAD, dplanck_BT(0, joker));
        }
      }
      else if(is_frequency_parameter(flag_partials[flag_partials_positions[iq]]))
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponent(partial_attenuation[iq](joker, 0),
                                                             partial_phase[iq](joker, 0), n, theta*DEG2RAD,
                                                             eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponent(partial_source[iq](joker, 0), phase(joker, 0), n,
                                                          theta*DEG2RAD, eta*DEG2RAD, planck_BT);
          
          nlte_dsource_dx[iq].AddZeemanSigmaPlusComponent(source(joker, 0), phase(joker, 0), n, theta*DEG2RAD,
                                                          eta*DEG2RAD, dplanck_BT(1, joker));
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] == JacPropMatType::VMR)
      {  
        if(species_match(flag_partials[flag_partials_positions[iq]], abs_species[this_species])) {
          dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponent(attenuation(joker, 0), phase(joker, 0),
                                                              n/abs_vmrs(this_species, 0), theta*DEG2RAD, eta*DEG2RAD);
          
          if(do_src) {
            dnlte_dx_source[iq].AddZeemanSigmaPlusComponent(source(joker, 0), phase(joker, 0),
                                                            n/abs_vmrs(this_species, 0), theta*DEG2RAD, eta*DEG2RAD,
                                                            planck_BT);
          }
        }
      }
      else if(flag_partials[flag_partials_positions[iq]] not_eq JacPropMatType::NotPropagationMatrixType)
      {
        dpropmat_clearsky_dx[iq].AddZeemanSigmaPlusComponent(partial_attenuation[iq](joker, 0),
                                                             partial_phase[iq](joker, 0), n, theta*DEG2RAD,
                                                             eta*DEG2RAD);
        
        if(do_src)
        {
          dnlte_dx_source[iq].AddZeemanSigmaPlusComponent(partial_source[iq](joker, 0), phase(joker, 0), n,
                                                          theta*DEG2RAD, eta*DEG2RAD, planck_BT);
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
  }
}


void set_magnetic_parameters_derivative(
    Numeric& dH_du,
    Numeric& dH_dv,
    Numeric& dH_dw,
    Numeric& deta_du,
    Numeric& deta_dv,
    Numeric& deta_dw,
    Numeric& dtheta_du,
    Numeric& dtheta_dv,
    Numeric& dtheta_dw,
    ConstVectorView rtp_mag,
    ConstVectorView r_path_los)
{
    const Numeric& Bu = rtp_mag[0],Bv = rtp_mag[1],Bw = rtp_mag[2];
    
    const Numeric
    Bu2=Bu*Bu, Bv2 = Bv*Bv, Bw2 = Bw*Bw, 
    H_mag = sqrt(Bu2+Bv2+Bw2), H2 = H_mag*H_mag,
    aa=DEG2RAD*r_path_los[1], za=DEG2RAD*r_path_los[0], 
    cosaa=cos(aa), cosza=cos(za), sinaa=sin(aa), sinza=sin(za),
    sinza2=sinza*sinza,cosza2=cosza*cosza,cosaa2=cosaa*cosaa,sinaa2=sinaa*sinaa,
    a=(Bw*cos(za) + Bu*cos(aa)*sin(za) + Bv*sin(aa)*sin(za)),
    b=-cosaa*sinza*Bv2 + Bu*sinaa*sinza*Bv - cosaa*sinza*Bw2 + Bu*cosza*Bw,
    term1=(Bu*(cosaa2*sinza2 - 1.) + Bw*cosaa*cosza*sinza + Bv*cosaa*sinaa*sinza2),
    term2=(Bv*(sinaa2*sinza2 - 1.) + Bw*sinaa*cosza*sinza + Bu*cosaa*sinaa*sinza2),
    term3=(Bw*(cosza2 - 1.) + Bu*cosaa*cosza*sinza + Bv*sinaa*cosza*sinza),
    eta_test=(PI - acos(((Bu*cosaa*cosza - Bw*sinza + 
    Bv*sinaa*cosza)*sqrt(term1*term1/(H_mag*H_mag) + term2*term2/(H_mag*H_mag) +
    term3*term3/(H_mag*H_mag)))/(H_mag)))*RAD2DEG; 
    
    const Numeric
    x1 = (Bv*cosaa - Bu*sinaa),
    x2 = -((Bu2*cosaa-Bv2*cosaa+2.0*Bu*Bv*sinaa)*cosaa*sinza2
    - Bu2 - Bw2 + 2.0*Bu*Bw*cosaa*cosza*sinza +
    (Bw2*cosza - Bv2*cosza+ 2.0*Bv*Bw*sinaa*sinza)*cosza),
    fx2=sqrt(x2),
    x3 = 1.0/H2,
    x=x1*fx2*x3, // eta = acos(x) * RAD2DEG;
    dx1_dBu = -sinaa,
    dx1_dBv = cosaa,
    dx2_dBu = 2.0*Bu - cosaa*sinza2*(2.0*Bu*cosaa + 2.0*Bv*sinaa) - 2.0*Bw*cosaa*cosza*sinza,
    dx2_dBv = cosza*(2.0*Bv*cosza - 2.0*Bw*sinaa*sinza) + cosaa*sinza2*(2.0*Bv*cosaa - 2.0*Bu*sinaa),
    dx2_dBw = 2.0*Bw - cosza*(2.0*Bw*cosza + 2.0*Bv*sinaa*sinza) - 2.0*Bu*cosaa*cosza*sinza,
    dx3_dBu = -2.0*Bu*x3*x3,
    dx3_dBv = -2.0*Bv*x3*x3,
    dx3_dBw = -2.0*Bw*x3*x3,
    c1 = fx2*x3, c2 = x1*0.5*x3/fx2,c3=x1*fx2,
    d_acos  = -1/sqrt(1-x*x) * (eta_test>90.0?-RAD2DEG:RAD2DEG); //RAD2DEG is here for convenience...
    
    dH_du = Bu/H_mag;
    dH_dv = Bv/H_mag;
    dH_dw = Bw/H_mag;
    
    dtheta_du = b/sqrt(H_mag*H_mag-a*a)*H_mag*H_mag * RAD2DEG;
    dtheta_dv = dtheta_du * (- sinaa*sinza*Bu2 + Bv*cosaa*sinza*Bu - sinaa*sinza*Bw2 + Bv*cosza*Bw)/b;
    dtheta_dw = dtheta_du * (-(cosza*Bu2 - Bw*cosaa*sinza*Bu + cosza*Bv2 - Bw*sinaa*sinza*Bv)/b);
    
    deta_du = d_acos * (dx1_dBu*c1 + dx2_dBu*c2 +  dx3_dBu*c3);
    deta_dv = d_acos * (dx1_dBv*c1 + dx2_dBv*c2 +  dx3_dBv*c3);
    deta_dw = d_acos * (/* empty */  dx2_dBw*c2 +  dx3_dBw*c3);
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
