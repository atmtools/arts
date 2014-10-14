/* Copyright (C) 2012
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

#include <cmath>
#include <stdexcept>
#include <string>
#include "auto_md.h"
#include "ppath.h"
#include "messages.h"
#include "math_funcs.h"
#include "absorption.h"
#include "abs_species_tags.h"
#include "physics_funcs.h"
#include "matpackIII.h"
#include "rte.h"
#include "rational.h"

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PLANCK_CONST;
extern const Numeric BOHR_MAGNETON;
extern const Numeric LANDE_GS;


/*!
 *  Defines the phase of the propagation matrix as derived by the author.
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
    
    const Numeric
    S2T = sin(theta)*sin(theta),  
    CT  = cos(theta), 
    CE2 = cos(2*eta), 
    SE2 = sin(2*eta);
    
    
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
        default: // Nil matrix since this should not be called.
            K=0;
            break;
    };
};


/*!
 *  Defines the attenuation of the propagation matrix as derived by the author.
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
    
    
    const Numeric 
    S2T  = sin(theta)*sin(theta),  
    C2T  = cos(theta)*cos(theta),  
    CT   = cos(theta), 
    CE2  = cos(2*eta), 
    SE2  = sin(2*eta);
    
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
        default: // Unity matrix in attenuation
            K(0,0) = 1;            K(0,1) =           0;  K(0,2) =           0;  K(0,3) =         0;
            K(1,0) = 0;            K(1,1) =           1;  K(1,2) =           0;  K(1,3) =         0;
            K(2,0) = 0;            K(2,1) =           0;  K(2,2) =           1;  K(2,3) =         0;
            K(3,0) = 0;            K(3,1) =           0;  K(3,2) =           0;  K(3,3) =         1;
            break;
    };
};


Numeric gs_caseb(const Numeric& N, const Numeric& J, const Numeric& S, const Numeric& GS) { return (GS*(J*(J+1.)+S*(S+1.)-N*(N+1.))/(J*(J+1.))/2.0); }
Numeric gs_casea(const Numeric& Omega, const Numeric& J, const Numeric& Sigma, const Numeric& GS) { return GS/2.0/J/(J+1)*Omega*(Omega+Sigma); }


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
Numeric relative_strength(const Rational& m, const Rational& j, const Index& dj, const Index& dm)
{
    // Numeric conversions are necessary.
    const Numeric J = (j-dj).toNumeric(), M = (m).toNumeric();
    
    // Variable to be returned.
    Numeric ret_val;
    
    switch ( dj )
    {
        case -1:
            switch ( dm )
            {
                case -1: // Transitions anti-parallel to the magnetic field
                    ret_val = (0.75)*(J+M)*(J-1+M)/(2.*J*(2.*J-1)*(2.*J+1));
                    break;
                case  0: // Transitions perpendicular to the magnetic field
                    ret_val = (1.50)*(J*J-M*M)/(J*(2.*J-1.)*(2.*J+1.));
                    break;
                case +1: // Transitions parallel to the magnetic field
                    ret_val = (0.75)*(J-M)*(J-1.-M)/(2.*J*(2.*J-1.)*(2.*J+1.));
                    break;
                default:
                    throw std::runtime_error("Something is extremely wrong.");
                    break;
            }
            break;
        case  0:
            switch ( dm )
            {
                case -1: // Transitions anti-parallel to the magnetic field
                    ret_val = (0.75)*(J+M)*(J+1.-M)/(2.*J*(J+1.)*(2.*J+1.));
                    break;
                case  0: // Transitions perpendicular to the magnetic field
                    ret_val = (1.50)*M*M/(J*(J+1.)*(2.*J+1.));
                    break;
                case +1: // Transitions parallel to the magnetic field
                    ret_val = (0.75)*(J-M)*(J+1.+M)/(2.*J*(J+1.)*(2.*J+1.));
                    break;
                default:
                    throw std::runtime_error("Something is extremely wrong.");
                    break;
            }
            break;
        case +1:
            switch ( dm )
            {
                case -1: // Transitions anti-parallel to the magnetic field
                    ret_val = (0.75)*(J+1.-M)*(J+2.-M)/(2.*(J+1.)*(2.*J+1.)*(2.*J+3.));
                    break;
                case  0: // Transitions perpendicular to the magnetic field
                    ret_val = (1.5)*((J+1.)*(J+1.)-M*M)/((J+1.)*(2.*J+1.)*(2.*J+3.));
                    break;
                case +1: // Transitions parallel to the magnetic field
                    ret_val = (0.75)*(J+1.+M)*(J+2.+M)/(2.*(J+1.)*(2.*J+1.)*(2.*J+3.));
                    break;
                default:
                    throw std::runtime_error("Something is extremely wrong.");
                    break;
            }
            break;
        default:
            throw std::runtime_error("Something is extremely wrong.");
            break;
    }

    return ret_val;
}


Numeric (*frequency_change)(const Rational&, const  Rational&, const Rational&, const Numeric&, const Index&, const Index&, const Index&, const Numeric&, const Numeric&);


/*!
 * Return the frequency change of the split Zeeman line parts as found from
 * g_s x M - g_s' x M'
 * 
 * Takes into account non-free electron GS constant by GS/2 * g_s
 * 
 * \param  n       In:     Main quantum number.
 * \param  m       In:     Secondary rotational quantum number.
 * \param  j       In:     Coupled rotational quantum number.
 * \param  s       In:     Electron total spin quantum number.
 * \param  DJ      In:     Change in the coupled rotational quantum number.
 * \param  DM      In:     Change in the secondary rotational quantum number.
 * \param  DN      In:     Change in the main rotational quantum number.
 * \param  H_mag   In:     Magnitude of the magnetic field in Tesla.
 * \param  GS      In:     G-constant from data for the molecule.
 * 
 * \author Richard Larsson
 * \date   2012-11-13
 */
Numeric frequency_change_caseb(const Rational& n, const Rational& m, const Rational& j, const Numeric& S, const Index& DJ, const Index& DM, const Index& DN, const Numeric& H_mag, const Numeric& GS)
{
    const Numeric       N_up = n.toNumeric();
    const Numeric       N_lo = N_up - (Numeric)DN;
    const Numeric       M_lo = m.toNumeric();
    const Numeric       M_up = M_lo + (Numeric)DM; // This will probably confuse people even though it is correct. Should fix so outer loop is over M_up.
    const Numeric       J_up = j.toNumeric();
    const Numeric       J_lo = J_up - (Numeric)DJ;
    
    // The special case is because g_s is undefined for J=0.
    return (!(j == 0 && DJ == -1)) ? 
    - H_mag * (M_lo) * gs_caseb(N_lo,J_lo,S,GS) / PLANCK_CONST * BOHR_MAGNETON 
    + H_mag * (M_up) * gs_caseb(N_up,J_up,S,GS) / PLANCK_CONST * BOHR_MAGNETON : 
      H_mag * (M_lo) * gs_caseb(N_lo,J_lo,S,GS) / PLANCK_CONST * BOHR_MAGNETON ;
}


/*!
 * Return the frequency change of the split Zeeman line parts as found from
 * g_s x M - g_s' x M'
 * 
 * Takes into account non-free electron GS constant by GS/2 * g_s
 * 
 * \param  omega   In:     Main quantum number.
 * \param  m       In:     Secondary rotational quantum number.
 * \param  j       In:     Coupled rotational quantum number.
 * \param  s       In:     Electron total spin quantum number.
 * \param  DJ      In:     Change in the coupled rotational quantum number.
 * \param  DM      In:     Change in the secondary rotational quantum number.
 * \param  DN      In:     Change in the main rotational quantum number.
 * \param  H_mag   In:     Magnitude of the magnetic field in Tesla.
 * \param  GS      In:     G-constant from data for the molecule.
 * 
 * \author Richard Larsson
 * \date   2013-08-09
 */
Numeric frequency_change_casea(const Rational& omega, const Rational& m, const Rational& j, const Numeric& Sigma, const Index& DJ, const Index& DM, const Index& Domega, const Numeric& H_mag, const Numeric& GS)
{
    const Numeric       Omega_up = omega.toNumeric();
    const Numeric       Omega_lo = Omega_up - (Numeric)Domega;
    const Numeric       M_lo = m.toNumeric();
    const Numeric       M_up = M_lo + (Numeric)DM; // This will probably confuse people even though it is correct. Should fix so outer loop is over M_up.
    const Numeric       J_up = j.toNumeric();
    const Numeric       J_lo = J_up - (Numeric)DJ;
    
    // This follows Berdyugina and Solnaki
    return
        H_mag * (M_up) * gs_casea(Omega_up,J_up,Sigma,GS) / PLANCK_CONST * BOHR_MAGNETON -
        H_mag * (M_lo) * gs_casea(Omega_lo,J_lo,Sigma,GS) / PLANCK_CONST * BOHR_MAGNETON ;
}


/*!
  Helper function. This is the only place where m_zeeman interacts with other absorption protocols.
 */
void xsec_species_line_mixing_wrapper_with_zeeman(  
        Tensor3View part_abs_mat, 
        const ArrayOfArrayOfSpeciesTag& abs_species, 
        const ArrayOfLineshapeSpec& abs_lineshape, 
        const ArrayOfLineRecord& lr, 
        const SpeciesAuxData& isotopologue_ratios, 
        const Matrix& abs_vmrs, 
        const Vector& abs_p,
        const Vector& abs_t, 
        const Vector& f_grid, 
        const Numeric& theta, 
        const Numeric& eta, 
        const Index& DM, 
        const Index& this_species,
        const Verbosity& verbosity )
{
    assert( part_abs_mat.npages() == f_grid.nelem() && part_abs_mat.ncols() == 4 && part_abs_mat.nrows() == 4 );
    
    Matrix A(f_grid.nelem(), 1, 0.);
    Matrix B(f_grid.nelem(), 1, 0.);
    
    
    for ( Index i=0; i<abs_species[this_species].nelem(); ++i )
    {
        Matrix attenuation(f_grid.nelem(), 1, 0.), phase(f_grid.nelem(), 1, 0.);;

        xsec_species_line_mixing_wrapper( attenuation, phase,
                f_grid, abs_p, abs_t, abs_vmrs, abs_species, this_species, lr,
                abs_lineshape[i].Ind_ls(), abs_lineshape[i].Ind_lsn(), abs_lineshape[i].Cutoff(),
                isotopologue_ratios, verbosity ); // Now in cross section

        attenuation *= abs_vmrs(this_species, 0) * number_density( abs_p[0],abs_t[0]); // Now in absorption coef.
        phase *= abs_vmrs(this_species, 0) * number_density( abs_p[0],abs_t[0]); // Now in absorption coef.
        phase *= 2; // phase matrix is twice as large according to sources.

        A += attenuation;
        B += phase;
    }
    Matrix  K_a(4,4), K_b(4,4);
    attenuation_matrix(K_a, theta*DEG2RAD, eta*DEG2RAD, DM);
    phase_matrix(K_b, theta*DEG2RAD, eta*DEG2RAD, DM);

    Tensor3 temp_part_abs_mat=part_abs_mat;
    mult(part_abs_mat,A(joker,0), K_a);// Resets part_abs_mat
    mult(temp_part_abs_mat,B(joker,0), K_b);

    part_abs_mat+=temp_part_abs_mat;

}


/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddZeeman(Tensor4& propmat_clearsky,
        const Vector& f_grid,
        const ArrayOfArrayOfSpeciesTag& abs_species,
        const ArrayOfArrayOfLineRecord& abs_lines_per_species,
        const ArrayOfLineshapeSpec& abs_lineshape,
        const SpeciesAuxData& isotopologue_ratios,
        const SpeciesAuxData& isotopologue_quantum,
        const Numeric& rtp_pressure,
        const Numeric& rtp_temperature,
        const Vector& rtp_vmr,
        const Vector& rtp_mag,
        const Vector& ppath_los,
        const Index& atmosphere_dim,
        const Index& manual_zeeman_tag,
        const Numeric& manual_zeeman_magnetic_field_strength,
        const Numeric& manual_zeeman_theta,
        const Numeric& manual_zeeman_eta,
        const Verbosity& verbosity)
{
    CREATE_OUT3;

    const Numeric margin    = 1e-4; // This margin is for relative strength and can perhaps be lowered by returning RS as Rational?
    bool          do_zeeman = false;

    // Check that correct isotopologue ratios are defined for the species
    // we want to calculate
    checkIsotopologueRatios(abs_species, isotopologue_ratios);

    {// Begin TEST(s)
    if (abs_species.nelem() != abs_lines_per_species.nelem())
        throw std::runtime_error("Dimension of *abs_species* and *abs_lines_per_species* don't match.");
    for(Index II = 0; II<abs_species.nelem(); II++)
        if(is_zeeman(abs_species[II])) { do_zeeman = true; break; } // If any species is Zeeman, do it.
    if( propmat_clearsky.ncols()  != 4 )
        throw std::runtime_error("Zeeman Effect is only implemented for Stokes dimension 4.");
    if( propmat_clearsky.nrows()  != 4 )
        throw std::runtime_error("Zeeman Effect is only implemented for Stokes dimension 4.");
    if( propmat_clearsky.npages() != f_grid.nelem() )
        throw std::runtime_error("Frequency dimension of *propmat_clearsky* not equal to length of *f_grid*.");
    if( propmat_clearsky.nbooks() != abs_species.nelem() )
        throw std::runtime_error("Species dimension of *propmat_clearsky* not equal to length of *abs_species*.");
    if( rtp_mag.nelem() != 3 )
        throw std::runtime_error("*rtp_mag* must have length 3.");
    if( atmosphere_dim != 3 )
        throw std::runtime_error("*atmosphere_dim* must be 3.  Zeeman Effect is only implemented for 3D geometry.");
    if( ppath_los.nelem() != 2 )
        throw std::runtime_error("*ppath_los* is not set correctly.");
    }// End   TEST(s)

    Vector R_path_los;
    mirror_los(R_path_los, ppath_los, atmosphere_dim);

    // Using the standard scalar absorption functions to get physics parameters,
    Vector abs_p, abs_t; Matrix abs_vmrs;
    AbsInputFromRteScalars( abs_p, abs_t, abs_vmrs,                        // Output
            rtp_pressure, rtp_temperature, rtp_vmr,  //Input
            verbosity);                                  // Verbose!
    if( do_zeeman  && (( rtp_mag[0]!=0 || rtp_mag[1]!=0 || rtp_mag[2]!=0 ) || manual_zeeman_tag != 0 ))
    {
        //Get the magnitude of the magnetic field and store a local unit Vector for simplified angle calculations.
        const Numeric H_mag = manual_zeeman_tag != 0?manual_zeeman_magnetic_field_strength:sqrt( rtp_mag * rtp_mag );

        Numeric theta, eta;
        if(manual_zeeman_tag!=0)
        { // Leaving it up to the user to manually tag the angles for simplified magnetic fields.
            eta   = manual_zeeman_eta;
            theta = manual_zeeman_theta;
        }
        else
        { // Getting angles from coordinate system
            Vector H(3);
            H  = rtp_mag;
            H /= H_mag;
            // Direction vector of radiation
            Numeric dx, dy, dz;
            zaaa2cart(dx,dy,dz,R_path_los[0],R_path_los[1]);
            // Radiation path direction as per Mishchenko.
            Vector R_path(3);
            R_path[0] = dx;
            R_path[1] = dy;
            R_path[2] = dz;
            // Vertical polarization direction as per Mishchenko.
            Vector e_v(3);
            e_v[0] =  cos(R_path_los[0] * DEG2RAD) * cos(R_path_los[1] * DEG2RAD);
            e_v[1] =  cos(R_path_los[0] * DEG2RAD) * sin(R_path_los[1] * DEG2RAD);
            e_v[2] = -sin(R_path_los[0] * DEG2RAD);
            // Horizontal polarization direction as per Mishchenko.
            Vector e_h(3);
            e_h[0] = -sin(R_path_los[1] * DEG2RAD);
            e_h[1] =  cos(R_path_los[1] * DEG2RAD);
            e_h[2] =  0;
            // Get the coordinate system used by Lenoir.
            Vector tmp(3);
            proj(tmp, R_path, H);
            Vector R11(3);
            R11  = H;
            R11 -= tmp;
            R11 *= sqrt(R11*R11);
            Vector R22(3);
            cross3(R22, R11, R_path);
            // Test if the rotation is clockwise or counterclockwise.
            const Numeric eta_test = vector_angle(R22, e_h);
            // Find the angle between Mishchenko vertical/horizontal and Lenoir vertical/horizontal
            eta      = (eta_test > 90.0)?-vector_angle(R22, e_v):vector_angle(R22, e_v);

            // Angle between radiation propagation and magnetic field for determining how the radiation is polarized..
            theta = vector_angle(H, R_path);
        }

        Numeric DF, RS; // Delta Frequency and Relative Strength

            // Reinitialize every loop to empty the set.
        ArrayOfLineRecord temp_abs_lines_sm, temp_abs_lines_sp, //sigma minus, sigma plus
                          temp_abs_lines_pi; // pi

        // For all species
        for(Index II = 0; II<abs_species.nelem(); II++)
        {
            temp_abs_lines_sm.resize(0);
            temp_abs_lines_sp.resize(0);
            temp_abs_lines_pi.resize(0);

            temp_abs_lines_sm.reserve(25000);
            temp_abs_lines_sp.reserve(25000);
            temp_abs_lines_pi.reserve(25000);

            // If the species isn't Zeeman, look at the next species
            if(!is_zeeman(abs_species[II])) continue;

            // Else loop over all the lines in the species.
            for (Index ii = 0; ii< abs_lines_per_species[II].nelem(); ii++)
            {
                    // local LineRecord
		    LineRecord temp_LR = abs_lines_per_species[II][ii];
		    const Rational J   = temp_LR.QuantumNumbers().Lower(QN_J);
		    const Numeric hund = isotopologue_quantum.getParam(temp_LR.Species(), temp_LR.Isotopologue(), 2);
		    Numeric S;
		    const Numeric GS   = isotopologue_quantum.getParam(temp_LR.Species(), temp_LR.Isotopologue(), 0);
		    const Index DJ     = (J - temp_LR.QuantumNumbers().Upper(QN_J)).toIndex();
		    Numeric RS_sum     = 0; //Sum relative strength (which ought be close to one by the end)
		    Index number_M_sm  = 0, number_M_sp=0, number_M_pi=0;
                    // Only look at lines which have no change in the main rotational number
                    // std::cout << "RSpi=[];RSsp=[];RSsm=[];DFpi=[];DFsm=[];DFsp=[];\n";
		    
		    Rational Main;
		    Index DMain;
		    
		    if( hund ==0 )//Case a
		    {
			Main  = temp_LR.QuantumNumbers().Lower(QN_Omega);
			DMain = (Main - temp_LR.QuantumNumbers().Upper(QN_Omega)).toIndex();
			frequency_change=frequency_change_casea;
			S = 1-Main.toNumeric();
		    }
		    else if( hund == 1 )// Case b
		    {
			Main  = temp_LR.QuantumNumbers().Lower(QN_N);
			DMain = (Main - temp_LR.QuantumNumbers().Upper(QN_N)).toIndex();
			frequency_change=frequency_change_caseb;
			S = isotopologue_quantum.getParam(temp_LR.Species(), temp_LR.Isotopologue(), 1);
		    }
		    else
		    {
			ostringstream os;
			os << "There are undefined Hund cases: " << abs_lines_per_species[II][ii] << "\nThe case is: "<<hund<<", allowed are (a): "<<0<<" and (b): " << 1<<"\n";
			throw std::runtime_error(os.str());
		    }
		    
                    if (!J.isUndefined() != 0 && !Main.isUndefined() != 0 ) // This means the lines are considered erroneous.
                    {

                        for ( Rational M = -J+DJ; M<=J-DJ; M++ )
                        {
                            /*
                               Note that:
				sp := sigma plus,  which means DM =  1
				sm := sigma minus, which means DM = -1
				pi := planar,      which means DM =  0
                             */
                            if ( DJ ==  1 )
                            { // Then all DM transitions possible for all M
                                DF =  frequency_change(Main, M, J, S, DJ, -1, DMain, H_mag, GS);
                                RS = relative_strength(M, J, DJ, -1);
                                RS_sum += RS;
                                temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                temp_abs_lines_sm.push_back(temp_LR);
                                number_M_sm++;
                                // std::cout << "RSsm=[RSsm " << RS << "];DFsm=[DFsm " << DF << "];\n";

                                DF =  frequency_change(Main, M, J, S, DJ,  0, DMain, H_mag, GS);
                                RS = relative_strength(M, J, DJ,  0);
                                RS_sum += RS;
                                temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                temp_abs_lines_pi.push_back(temp_LR);
                                number_M_pi++;
                                // std::cout << "RSpi=[RSpi " << RS << "];DFpi=[DFpi " << DF << "];\n";

                                DF =  frequency_change(Main, M, J, S, DJ, +1, DMain, H_mag, GS);
                                RS = relative_strength(M, J, DJ, +1);
                                RS_sum += RS;
                                temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                temp_abs_lines_sp.push_back(temp_LR);
                                number_M_sp++;
                                // std::cout << "RSsp=[RSsp " << RS << "];DFsp=[DFsp " << DF << "];\n";
                            }
                            else if ( DJ ==  0 )
                            { // Then all DM transitions possible for all M
                                DF =  frequency_change(Main, M, J, S, DJ, -1, DMain, H_mag, GS);
                                RS = relative_strength(M, J, DJ, -1);
                                RS_sum += RS;
                                temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                temp_abs_lines_sm.push_back(temp_LR);
                                number_M_sm++;
                                // std::cout << "RSsm=[RSsm " << RS << "];DFsm=[DFsm " << DF << "];\n";
                                if( ! (M == 0) )
                                {
                                    DF =  frequency_change(Main, M, J, S, DJ,  0, DMain, H_mag, GS);
                                    RS = relative_strength(M, J, DJ,  0);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_pi.push_back(temp_LR);
                                    number_M_pi++;
                                    // std::cout << "RSpi=[RSpi " << RS << "];DFpi=[DFpi " << DF << "];\n";
                                }

                                DF =  frequency_change(Main, M, J, S, DJ, +1, DMain, H_mag, GS);
                                RS = relative_strength(M, J, DJ, +1);
                                RS_sum += RS;
                                temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                temp_abs_lines_sp.push_back(temp_LR);
                                number_M_sp++;
                                // std::cout << "RSsp=[RSsp " << RS << "];DFsp=[DFsp " << DF << "];\n";
                            }
                            else if ( DJ == -1 )
                            { // Then certain M results in blocked DM transitions
                                if ( M == -J + DJ && M!=0 )
                                { // Lower limit M only allows DM = 1
                                    DF =  frequency_change(Main, M, J, S, DJ, +1, DMain, H_mag, GS);
                                    RS = relative_strength(M, J, DJ, +1);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sp.push_back(temp_LR);
                                    number_M_sp++;
                                    // std::cout << "RSsp=[RSsp " << RS << "];DFsp=[DFsp " << DF << "];\n";

                                }
                                else if ( M == -J + DJ + 1 && M!=0 )
                                { // Next to lower limit M can only allow DM = 1, 0
                                    DF =  frequency_change(Main, M, J, S, DJ, +1, DMain, H_mag, GS);
                                    RS = relative_strength(M, J, DJ, +1);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sp.push_back(temp_LR);
                                    number_M_sp++;
                                    // std::cout << "RSsp=[RSsp " << RS << "];DFsp=[DFsp " << DF << "];\n";

                                    DF =  frequency_change(Main, M, J, S, DJ,  0, DMain, H_mag, GS);
                                    RS = relative_strength(M, J, DJ,  0);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_pi.push_back(temp_LR);
                                    number_M_pi++;
                                    // std::cout << "RSpi=[RSpi " << RS << "];DFpi=[DFpi " << DF << "];\n";
                                }
                                else if ( M ==  J - DJ - 1 && M!=0 )
                                { // Next to upper limit M can only allow DM = 0, -1
                                    DF =  frequency_change(Main, M, J, S, DJ,  0, DMain, H_mag, GS);
                                    RS = relative_strength(M, J, DJ,  0);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_pi.push_back(temp_LR);
                                    number_M_pi++;
                                    // std::cout << "RSpi=[RSpi " << RS << "];DFpi=[DFpi " << DF << "];\n";

                                    DF =  frequency_change(Main, M, J, S, DJ, -1, DMain, H_mag, GS);
                                    RS = relative_strength(M, J, DJ, -1);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sm.push_back(temp_LR);
                                    number_M_sm++;
                                    // std::cout << "RSsm=[RSsm " << RS << "];DFsm=[DFsm " << DF << "];\n";
                                }
                                else if ( M == J - DJ && M!=0 )
                                { // Upper limit M only allow DM = -1
                                    DF =  frequency_change(Main, M, J, S, DJ, -1, DMain, H_mag, GS);
                                    RS = relative_strength(M, J, DJ, -1);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sm.push_back(temp_LR);
                                    number_M_sm++;
                                    // std::cout << "RSsm=[RSsm " << RS << "];DFsm=[DFsm " << DF << "];\n";
                                }
                                else if( (-J + DJ + 1) ==  (J - DJ - 1) && M == 0)
                                { // Special case for N=1, J=0, M=0. Only allows DM = 0
                                    DF =  frequency_change(Main, M, J, S, DJ,  0, DMain, H_mag, GS);
                                    RS = relative_strength(M, J, DJ,  0);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_pi.push_back(temp_LR);
                                    number_M_pi++;
                                    // std::cout << "RSpi=[RSpi " << RS << "];DFpi=[DFpi " << DF << "];\n";
                                }
                                else
                                { // All DM transitions are possible for these M(s)
                                    DF =  frequency_change(Main, M, J, S, DJ, +1, DMain, H_mag, GS);
                                    RS = relative_strength(M, J, DJ, +1);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sp.push_back(temp_LR);
                                    number_M_sp++;
                                    // std::cout << "RSsp=[RSsp " << RS << "];DFsp=[DFsp " << DF << "];\n";

                                    DF =  frequency_change(Main, M, J, S, DJ,  0, DMain, H_mag, GS);
                                    RS = relative_strength(M, J, DJ,  0);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_pi.push_back(temp_LR);
                                    number_M_pi++;
                                    // std::cout << "RSpi=[RSpi " << RS << "];DFpi=[DFpi " << DF << "];\n";

                                    DF =  frequency_change(Main, M, J, S, DJ, -1, DMain, H_mag, GS);
                                    RS = relative_strength(M, J, DJ, -1);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sm.push_back(temp_LR);
                                    number_M_sm++;
                                    // std::cout << "RSsm=[RSsm " << RS << "];DFsm=[DFsm " << DF << "];\n";
                                }
                            }
                            else
                            { // The tests above failed and catastrophe follows
                                ostringstream os;
                                os << "There seems to be something wrong with the quantum numbers of at least one line in your *abs_lines*. " <<
                                    "Make sure this is a Zeeman line.\nThe upper quantum numbers are: " << 
                                    abs_lines_per_species[II][ii].QuantumNumbers().Upper() <<
                                    "\nThe lower quantum numbers are: " <<
                                    abs_lines_per_species[II][ii].QuantumNumbers().Lower() <<
                                    "\nThe entire line information: " << abs_lines_per_species[II][ii];
                                throw std::runtime_error(os.str());
                            }
                        }

                        if (abs(RS_sum-1.)>margin) //Reasonable confidence?
                        {
                            ostringstream os;
                            os << "The sum of relative strengths is not close to one. This is severly problematic and "
                                "you should look into why this happens.\nIt is currently " << RS_sum << " with DJ: "<<DJ<<", DMain: "<<DMain<<" for line: "<<
                                abs_lines_per_species[II][ii] <<"\n";
                            throw std::runtime_error(os.str());
                        }
                    }
                    else
                    {
                        ostringstream os;
                        os << "There are undefined quantum numbers in the line: " << abs_lines_per_species[II][ii] << "\nJ is "<<J<<" and Main is "<<Main<<std::endl;
                        throw std::runtime_error(os.str());
                    }
            }
            Tensor3 part_abs_mat(f_grid.nelem(), 4, 4);

            // Add Pi contribution to final propmat_clearsky
            xsec_species_line_mixing_wrapper_with_zeeman( part_abs_mat, abs_species, abs_lineshape,
                    temp_abs_lines_pi, isotopologue_ratios,
                    abs_vmrs, abs_p, abs_t, f_grid,
                    theta, eta, 0, II, 
                    verbosity );
            propmat_clearsky(II, joker, joker, joker) += part_abs_mat;

            // Add Sigma minus contribution to final propmat_clearsky
            xsec_species_line_mixing_wrapper_with_zeeman( part_abs_mat, abs_species, abs_lineshape,
                    temp_abs_lines_sm, isotopologue_ratios,
                    abs_vmrs, abs_p, abs_t, f_grid,
                    theta, eta, -1, II,
                    verbosity );
            propmat_clearsky(II, joker, joker, joker) += part_abs_mat;

            // Add Sigma plus contribution to final propmat_clearsky
            xsec_species_line_mixing_wrapper_with_zeeman( part_abs_mat, abs_species, abs_lineshape,
                    temp_abs_lines_sp, isotopologue_ratios,
                    abs_vmrs, abs_p, abs_t, f_grid,
                    theta, eta, 1, II, 
                    verbosity );
            propmat_clearsky(II, joker, joker, joker) += part_abs_mat;

        }
    }
    else // if the magnetic field is ignored
    {
        for(Index II = 0; II<abs_species.nelem(); II++)
        {
            // If the species isn't Zeeman, look at the next species.
            if(!is_zeeman(abs_species[II])) continue;
            
            Tensor3 part_abs_mat(f_grid.nelem(), 4, 4);
            xsec_species_line_mixing_wrapper_with_zeeman(  part_abs_mat, abs_species, abs_lineshape,
                                                            abs_lines_per_species[II], isotopologue_ratios,
                                                            abs_vmrs, abs_p, abs_t, f_grid,
                                                            0,0,1023, II, 
                                                            verbosity );
            propmat_clearsky(II, joker, joker, joker) += part_abs_mat;
        }
    }
}
