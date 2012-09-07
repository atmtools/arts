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
#include "arts.h"
#include "ppath.h"
#include "messages.h"
#include "math_funcs.h"
#include "absorption.h"
#include "abs_species_tags.h"
#include "matpackIII.h"
#include "rte.h"

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;



/*!
    Helper function for sort. Used to sort the LineRecord(s) in an
    ArrayOfLineRecord by increasing frequency. Usage example:

    sort(ArrayOfLineRecord_A.begin(), ArrayOfLineRecord_A.end(), sortF)
    returns ArrayOfLineRecord_A sorted by frequency

    \param  i   In:    LineRecord helper 1.
    \param  j   In:    LineRecord helper 2.

    \author Richard Larsson
    \date   2012-08-03
*/
bool sortF(LineRecord i, LineRecord j) { return ( i.F() < j.F() ); };

/*!
    Defines the rotation extinction matrix as derived by the author from
    the definitions of Lenoir (1967) and Mishchenko (2002).

    \param  K       Out:    The rotation extinction matrix.
    \param  theta   In:     Angle between the magnetic field and the
                            propagation path. In radians.
    \param  eta     In:     Angle to rotate planar polarization clockwise to
                            fit the general coordinate system. In radians.
    \param  DM      In:     Change in the secondary rotational quantum number.

    \author Richard Larsson
    \date   2012-08-03
*/
void K_mat(MatrixView K, const Numeric theta, const Numeric eta, const Index DM)
{
    assert(K.nrows() == 4 );
    assert(K.ncols() == 4 );
    
    switch( DM )
    {
        case -1:
            K(0,0) =   1 + cos(theta)*cos(theta);
            K(0,1) =   sin(theta)*sin(theta) * cos( 2 * eta );
            K(0,2) =   sin(theta)*sin(theta) * sin( 2 * eta );
            K(0,3) =   2 * cos(theta);

            K(1,0) = K(0,1);
            K(1,1) = K(0,0);
            K(1,2) = 0;
            K(1,3) = 0;

            K(2,0) = K(0,2);
            K(2,1) = K(1,2);
            K(2,2) = K(0,0);
            K(2,3) = 0;

            K(3,0) = K(0,3);
            K(3,1) = K(1,3);
            K(3,2) = K(2,3);
            K(3,3) = K(0,0);
            break;
        case  1:
            K(0,0) =   1 + cos(theta)*cos(theta);
            K(0,1) =   sin(theta)*sin(theta) * cos( 2 * eta );
            K(0,2) =   sin(theta)*sin(theta) * sin( 2 * eta );
            K(0,3) =  -2*cos(theta);

            K(1,0) = K(0,1);
            K(1,1) = K(0,0);
            K(1,2) = 0;
            K(1,3) = 0;

            K(2,0) = K(0,2);
            K(2,1) = K(1,2);
            K(2,2) = K(0,0);
            K(2,3) = 0;

            K(3,0) = K(0,3);
            K(3,1) = K(1,3);
            K(3,2) = K(2,3);
            K(3,3) = K(0,0);
            break;
        case  0:
            K(0,0) =   sin(theta)*sin(theta);
            K(0,1) = - sin(theta)*sin(theta) * cos( 2 * eta );
            K(0,2) = - sin(theta)*sin(theta) * sin( 2 * eta );
            K(0,3) = 0;

            K(1,0) = K(0,1);
            K(1,1) = K(0,0);
            K(1,2) = 0;
            K(1,3) = 0;

            K(2,0) = K(0,2);
            K(2,1) = K(1,2);
            K(2,2) = K(0,0);
            K(2,3) = 0;

            K(3,0) = K(0,3);
            K(3,1) = K(1,3);
            K(3,2) = K(2,3);
            K(3,3) = K(0,0);
            break;
        default: // Unity matrix
            K(0,0) = 1;
            K(0,1) = 0;
            K(0,2) = 0;
            K(0,3) = 0;

            K(1,0) = K(0,1);
            K(1,1) = K(0,0);
            K(1,2) = 0;
            K(1,3) = 0;

            K(2,0) = K(0,2);
            K(2,1) = K(1,2);
            K(2,2) = K(0,0);
            K(2,3) = 0;

            K(3,0) = K(0,3);
            K(3,1) = K(1,3);
            K(3,2) = K(2,3);
            K(3,3) = K(0,0);
            break;
    };
};

/*!
    Return the relative strength of the split Zeeman line parts as found in
    Table 2 of Liebe and Hufford (1989). Renormalized to one instead of two.

    \param  n       In:     Main quantum number.
    \param  m       In:     Secondary rotational quantum number.
    \param  DJ      In:     Change in the main rotational quantum number.
    \param  DM      In:     Change in the secondary rotational quantum number.

    \author Richard Larsson
    \date   2012-08-03
*/
Numeric RelativeStrength(Index n, Index m, Index DJ, Index DM)
{

    /*
    The following switch case is from table 2 of Liebe and Hufford, 1989.
    */
    
    Numeric N = (Numeric)n, M = (Numeric)m, relstr;
    
    switch ( DJ )
    {//* Switch over DJ, if DJ != 1 or -1, DJ is wrong.
        case 1:
            switch ( DM )
             {//* Switch over DM, if DM != 1, 0 or -1, DM is wrong.
                case  1:
                    relstr = 3.0 * ( N+M+1 ) * ( N+M+2 ) / ( 4.0 * ( N+1 ) * ( 2*N+1 ) * ( 2*N+3 ) ) ;
                    break;
                case  0:
                    relstr = 3.0 * ( (N+1) * (N+1) - M * M) / ( ( N+1 ) * ( 2*N+1 ) * ( 2*N+3 ) );
                    break;
                case -1:
                    relstr = 3.0 * ( N-M+1 ) * ( N-M+2 ) / ( 4.0 * ( N+1 ) * ( 2*N+1 ) * ( 2*N+3 ) ) ;
                    break;
                default:
                    relstr = -1e6;
                    break;
            };
            break;
        case -1:
            switch ( DM )
            {//* Switch over DM, if DM != 1, 0 or -1, DM is wrong.
                case  1:
                    relstr = 3.0 * ( N-M-1 ) * ( N-M ) / ( 4.0 * N * ( 2*N+1 ) * ( 2*N-1 ) );
                    break;
                case  0:
                    relstr = 3.0 * ( N * N - M * M ) / ( N * ( 2*N+1 ) * ( 2*N-1 ) );
                    break;
                case -1:
                    relstr = 3.0 * ( N+M-1 ) * ( N+M ) / ( 4.0 * N * ( 2*N+1 ) * ( 2*N-1 ) );
                    break;
                default:
                    relstr = -1e6;
                    break;
            };
            break;
        default:
            relstr = -1e6;
            break;
    };
    return relstr / 2.0; //FIXME? It seems the normalization used by Hufford and Liebe is to 2 and not 1... this factor might be erroneous?
};

/*!
    Return the frequency change of the split Zeeman line parts as found in
    Table 2 of Liebe and Hufford (1989).

    \param  n       In:     Main quantum number.
    \param  m       In:     Secondary rotational quantum number.
    \param  DJ      In:     Change in the main rotational quantum number.
    \param  DM      In:     Change in the secondary rotational quantum number.
    \param  H_mag   In:     Magnitude of the magnetic field in Tesla.

    \author Richard Larsson
    \date   2012-08-03
*/
Numeric FrequencyChange(Index n, Index m, Index DJ, Index DM, Numeric H_mag)
{
    /*
    The following switch case is from table 2 of Liebe and Hufford, 1989.
    */

    Numeric N = (Numeric)n, M = (Numeric)m, fcc;
    Numeric KH = 2.8026 * 1000000 * H_mag * 10000; // KH from Lenoir 1968. Factors are for conversion from Lenoir's' to SI units.

    switch ( DJ )
    {//* Switch over DJ, if DJ != 1 or -1, DJ is wrong.
        case 1:
            switch ( DM )
            {//* Switch over DM, if DM != 1, 0 or -1, DM is wrong.
                case 1:
                    fcc = 1.0 * ( M * ( N-1 ) - N ) / ( N * ( N+1 ) ) ;
                    break;
                case  0:
                    fcc = 1.0 * ( M * ( N-1 ) ) / ( N * ( N+1 ) ) ;
                    break;
                case  -1:
                    fcc = 1.0 * ( M * ( N-1 ) + N ) / ( N * ( N+1 ) ) ;
                    break;
                default:
                    fcc = -1e10;
                    break;
            };
            break;
        case -1:
            switch ( DM )
            {//* Switch over DM, if DM != 1, 0 or -1, DM is wrong.
                case 1:
                    fcc = 1.0 * ( M * ( N+2 ) - 1 ) / ( N * ( N+1 ) ) ;
                    break;
                case  0:
                    fcc = 1.0 * ( M * ( N+2 ) ) / ( N * ( N+1 ) ) ;
                    break;
                case  -1:
                    fcc = 1.0 * ( M * ( N+2 ) + 1 ) / ( N * ( N+1 ) ) ;
                    break;
                default:
                    fcc = -1e10;
                    break;
            };
            break;
        default:
            fcc = -1e10;
            break;
    };
    return KH*fcc;
};

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_mat_per_speciesAddZeemanLBL(Tensor4& abs_mat_per_species,
                                  const Vector& f_grid,
                                  const ArrayOfArrayOfSpeciesTag& abs_species,
                                  const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                  const ArrayOfLineshapeSpec& abs_lineshape,
                                  const ArrayOfVector& isotopologue_ratios,
                                  const Numeric& rte_pressure,
                                  const Numeric& rte_temperature,
                                  const Vector& rte_vmr_list,
                                  const Numeric& rte_doppler,
                                  const Vector& rte_mag,
				                  const Vector& ppath_los,
				                  const Index& atmosphere_dim,
                                  const Verbosity& verbosity)
{
    CREATE_OUT3;
    Vector R_path_los;
    mirror_los(R_path_los, ppath_los, atmosphere_dim);
    
    const Numeric margin    = 1e-3;
    bool          do_zeeman = false;

    /*
        This function will, for each Zeeman species, make a local
        ArrayOfLineRecord(s) for the various transition types with Zeeman
        altered LineRecord(s).  These are then composed into a single
        ArrayOfArrayOfLineRecord which is processed as per the scalar case.
        
        The line broadened absorption coefficients are finally multiplied with
        the transition type rotation matrix and the new variable is returned.
        
        Note that between 55 GHz and 65 GHz there is usually ~700 O_2 lines,
        however, when this Zeeman splitting method is used, the number of
        lines is increased to about 45,000. This is a time consuming method...
    */

    // Begin TEST(s)
    for(Index II = 0; II<abs_species.nelem(); II++)
        if(abs_species[II][0].Zeeman()) { do_zeeman = true; break; } // If any species is Zeeman, do it.
    if( abs_mat_per_species.ncols()  != 4 )
        throw runtime_error("Zeeman Effect is only implemented for Stokes dimension 4.");
    if( abs_mat_per_species.nrows()  != 4 )
        throw runtime_error("Zeeman Effect is only implemented for Stokes dimension 4.");
    if( abs_mat_per_species.npages() != f_grid.nelem() )
        throw runtime_error("Frequency dimension of abs_mat_per_species not equal to length of f_grid.");
    if( abs_mat_per_species.nbooks() != abs_species.nelem() )
        throw runtime_error("Species dimension of abs_mat_per_species not equal to length of abs_species.");
    if( rte_mag.nelem() != 3 )
    {
        if( rte_mag.nelem() == 1 )
            if(abs(rte_mag[0]+1) < margin)
                do_zeeman = false; //If no magnetic field, do not do Zeeman
            else
                throw runtime_error("rte_mag have the length of a flag"
                " but not the content of such.");
        else
            throw runtime_error("rte_mag does not have the length of"
                    " a flag nor the length of a magnetic field vector.");
    }
    // End   TEST(s)
    Vector local_f_grid;
    // Make pointer point to original.
    const Vector* f_grid_pointer = &f_grid;

    /*
        Doppler treatment, do this only if there is a non-zero Doppler
        shift. We do this after the frequency selection, so in the case
        that we have only a single frequency, we have to shift only that!

        Unfortunately, we need yet another local copy of f_grid. In
        contrast to the frequency selection, we here want to modify the
        actual frequency values inside!
    */
    Vector local_doppler_f_grid;
    if (rte_doppler==0)
    {
        out3 << "  Doppler shift: None\n";
    }
    else
    {
        ostringstream os;
        os << "  Doppler shift: " << rte_doppler << " Hz\n";
        out3 << os.str();

        Numeric local_doppler;
        NumericScale( local_doppler, rte_doppler, -1, verbosity );
        // I could just have multiplied by -1 directly, but I like using
        // the WSM here.

        VectorAddScalar( local_doppler_f_grid,  *f_grid_pointer, local_doppler, verbosity );

        // Make pointer point to the doppler shifted frequency grid.
        f_grid_pointer = &local_doppler_f_grid;
    }

    if(do_zeeman == 1)
    {
        //Get the magnitude of the magnetic field and store a local unit Vector for simplified angle calculations.
        const Numeric H_mag = sqrt( rte_mag * rte_mag );
        Vector H(3);
            H  = rte_mag;
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
        Vector temp(3);
            proj(temp, R_path, H);
        Vector R11(3);
            R11  = H;
            R11 -= temp;
            R11 *= sqrt(R11*R11);
        Vector R22(3);
            cross3(R22, R11, R_path);

        // Test if the rotation is clockwise or counterclockwise.
        const Numeric eta_test = vector_angle(R22, e_h);
        // Find the angle between Mishchenko vertical/horizontal and Lenoir vertical/horizontal
        const Numeric eta      = (eta_test > 90.0)?-vector_angle(R22, e_v):vector_angle(R22, e_v);
        const Numeric eta_same = (eta_test > 90.0)?-vector_angle(R11, e_h):vector_angle(R11, e_h);

        // Angle between radiation propagation and magnetic field for determining how the radiation is polarized..
        const Numeric theta = vector_angle(H, R_path);

        Numeric DF, RS; // Delta Frequency and Relative Strength

        // For all species
        for(Index II = 0; II<abs_species.nelem(); II++)
        {
            // Reinitialize every loop to empty the set.
            ArrayOfLineRecord temp_abs_lines_sm, temp_abs_lines_sp, //sigma minus, sigma plus
                              temp_abs_lines_pi, temp_abs_lines_dt; // pi, default

            // If the species isn't Zeeman, look at the next species
            if(!abs_species[II][0].Zeeman()) continue;
            // Else loop over all the lines in the species.
            for (Index ii = 0; ii< abs_lines_per_species[II].nelem(); ii++)
            {
                // local LineRecord
                LineRecord temp_LR = abs_lines_per_species[II][ii];
                const Index J  = temp_LR.Lower_J();
                const Index N  = temp_LR.Lower_N();
                const Index DJ = J - temp_LR.Upper_J();

                // Only look at lines which have no change in the main rotational number
                if (J != -1 && N != -1)
                {
                    if ( temp_LR.Upper_N() == N )
                    {
                        // Begin TEST(s)
                        if ( N <= 0 )
                            throw runtime_error("The main quantum number cannot be nil!?");
                        if (abs(DJ) != 1){ temp_abs_lines_dt.push_back(temp_LR); continue; }
                        if ( J-DJ != N ) //FIXME: default K??? throw runtime_error?
                        { // Since Hitran12 beta version, this does no longer apply
                            out3 << "Physics ERROR: J-DJ = " << J - DJ << " and N = "
                                << N << " for ii = " << ii << ". The string is:\n"
                                << temp_LR.Lower_LQuanta() <<"\nThe total line is:\n"
                                <<temp_LR<<"\n";
                            continue;
                        }
                        // End   TEST(s)

                        // For every upper molecular magnetic moment, M
                        for ( Index M = -J+DJ; M<=J-DJ; M++ )
                        {
                            /*
                                Note that:
                                    sp := sigma plus,  which means DM =  1
                                    sm := sigma minus, which means DM = -1
                                    pi := planar,      which means DM =  0
                            */
                            if ( DJ ==  1 )
                            { // Then all DM transitions possible for all M
                                DF =  FrequencyChange(N, M, DJ, -1, H_mag);
                                RS = RelativeStrength(N, M, DJ, -1);
                                temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                temp_abs_lines_sm.push_back(temp_LR);

                                DF =  FrequencyChange(N, M, DJ,  0, H_mag);
                                RS = RelativeStrength(N, M, DJ,  0);
                                temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                temp_abs_lines_pi.push_back(temp_LR);

                                DF =  FrequencyChange(N, M, DJ, +1, H_mag);
                                RS = RelativeStrength(N, M, DJ, +1);
                                temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                temp_abs_lines_sp.push_back(temp_LR);
                            }
                            else if ( DJ == -1 )
                            { // Then certain M results in blocked DM transitions
                                if ( M == -J + DJ && M!=0 )
                                { // Lower limit M only allows DM = 1
                                    DF =  FrequencyChange(N, M, DJ, +1, H_mag);
                                    RS = RelativeStrength(N, M, DJ, +1);
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sp.push_back(temp_LR);
                                }
                                else if ( M == -J + DJ + 1 && M!=0 )
                                { // Next to lower limit M can only allow DM = 1, 0
                                    DF =  FrequencyChange(N, M, DJ, +1, H_mag);
                                    RS = RelativeStrength(N, M, DJ, +1);
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sp.push_back(temp_LR);

                                    DF =  FrequencyChange(N, M, DJ,  0, H_mag);
                                    RS = RelativeStrength(N, M, DJ,  0);
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_pi.push_back(temp_LR);
                                }
                                else if ( M ==  J - DJ - 1 && M!=0 )
                                { // Next to upper limit M can only allow DM = 0, -1
                                    DF =  FrequencyChange(N, M, DJ,  0, H_mag);
                                    RS = RelativeStrength(N, M, DJ,  0);
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_pi.push_back(temp_LR);

                                    DF =  FrequencyChange(N, M, DJ, -1, H_mag);
                                    RS = RelativeStrength(N, M, DJ, -1);
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sm.push_back(temp_LR);
                                }
                                else if ( M == J - DJ && M!=0 )
                                { // Upper limit M only allow DM = -1
                                    DF =  FrequencyChange(N, M, DJ, -1, H_mag);
                                    RS = RelativeStrength(N, M, DJ, -1);
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sm.push_back(temp_LR);
                                }
                                else if( (-J + DJ + 1) ==  (J - DJ - 1) && M == 0)
                                { // Special case for N=1, J=0, M=0. Only allows DM = 0
                                    DF =  FrequencyChange(N, M, DJ,  0, H_mag);
                                    RS = RelativeStrength(N, M, DJ,  0);
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_pi.push_back(temp_LR);
                                }
                                else
                                { // All DM transitions are possible for these M(s)
                                    DF =  FrequencyChange(N, M, DJ, +1, H_mag);
                                    RS = RelativeStrength(N, M, DJ, +1);
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sp.push_back(temp_LR);

                                    DF =  FrequencyChange(N, M, DJ,  0, H_mag);
                                    RS = RelativeStrength(N, M, DJ,  0);
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_pi.push_back(temp_LR);

                                    DF =  FrequencyChange(N, M, DJ, -1, H_mag);
                                    RS = RelativeStrength(N, M, DJ, -1);
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sm.push_back(temp_LR);
                                }
                            }
                            else
                            { // The tests above failed and catastrophe follows
                                throw runtime_error("If this happens, something is horribly wrong... did some of the tests above fail?");
                            }
                        }
                    }
                    else
                        temp_abs_lines_dt.push_back(temp_LR);
                }
                else
                    temp_abs_lines_dt.push_back(temp_LR);
            }

            // Sort ArrayOfLineRecord(s) by frequency [low → high].
            sort(temp_abs_lines_pi.begin(), temp_abs_lines_pi.end(), sortF);
            sort(temp_abs_lines_sp.begin(), temp_abs_lines_sp.end(), sortF);
            sort(temp_abs_lines_sm.begin(), temp_abs_lines_sm.end(), sortF);
            sort(temp_abs_lines_dt.begin(), temp_abs_lines_dt.end(), sortF);
            // Treat the DM transitions as separate species,
            ArrayOfArrayOfLineRecord temp_abs_lines_per_species(4);
                temp_abs_lines_per_species[0] = temp_abs_lines_pi;
                temp_abs_lines_per_species[1] = temp_abs_lines_sp;
                temp_abs_lines_per_species[2] = temp_abs_lines_sm;
                temp_abs_lines_per_species[3] = temp_abs_lines_dt;

            out3 << "From initially having: " << abs_lines_per_species[II].nelem() << " lines, \nit is now split up into: "
                 << temp_abs_lines_pi.nelem()+temp_abs_lines_sp.nelem()+temp_abs_lines_sm.nelem()+temp_abs_lines_dt.nelem()
                 << " lines.\nOut of these, \nPI stands for: " << temp_abs_lines_pi.nelem()
                 << ",\nSM stands for: " << temp_abs_lines_sm.nelem()
                 << ",\nSP stands for: " << temp_abs_lines_sp.nelem()
                 << " and \nDT stands for: " << temp_abs_lines_dt.nelem() << "\n";
                 
            // with the same species information,
            ArrayOfArrayOfSpeciesTag temp_abs_species(4);
            for (Index i=0; i<4; ++i) {
                temp_abs_species[i].resize(1);
                temp_abs_species[i][0] = SpeciesTag("O2");
                // We are setting the tag explicitly to O2 here. The function
                // anyway currently works only for this molecule. By setting it
                // like this (without Zeeman flag) we ensure that we can use the
                // standard LBL function to calculate absorption.
            }
            
            // and the same volume mixing ratios.
            Vector temp_vmrs(4);
                temp_vmrs[0] = rte_vmr_list[II]; temp_vmrs[1] = rte_vmr_list[II];
                temp_vmrs[2] = rte_vmr_list[II]; temp_vmrs[3] = rte_vmr_list[II];

            // Using the standard scalar absorption functions to get physics parameters,
            Vector abs_p, abs_t; Matrix abs_vmrs;
                AbsInputFromRteScalars(abs_p, abs_t, abs_vmrs, // Output
                rte_pressure, rte_temperature, temp_vmrs,     //Input
                verbosity);                                  // Verbose!
            // and then initialize the cross section per species variable.
            ArrayOfMatrix abs_xsec_per_species;
                abs_xsec_per_speciesInit(abs_xsec_per_species, //Output
                temp_abs_species, *f_grid_pointer, abs_p,     //Input
                verbosity);                                  //Verbose!

            // then calculate the cross section per species
            abs_xsec_per_speciesAddLines(abs_xsec_per_species, //Output
                                         temp_abs_species, *f_grid_pointer, abs_p, abs_t, abs_vmrs, //Input
                                         temp_abs_lines_per_species, abs_lineshape,
                                         isotopologue_ratios,
                                         verbosity);

            // and take continua into account.
            // 2012-9-4 Stefan: I don't think this function should add
            //          continua, since they are already added elsewhere.
            //          This function should just do the Zeeman LBL part, as
            //          The name says.
            //
            //            abs_xsec_per_speciesAddConts(abs_xsec_per_species,                                    //Output
            //                                         temp_abs_species, *f_grid_pointer, abs_p, abs_t,        //Input
            //                                         abs_n2, abs_h2o, abs_vmrs,                             //Input
            //                                         abs_cont_names, abs_cont_parameters, abs_cont_models, //Input
            //                                         verbosity);                                          //Verbose!

            // We can finally obtain the coefficients per species!
            Matrix abs_coef; ArrayOfMatrix abs_coef_per_species;
                abs_coefCalcFromXsec(abs_coef, abs_coef_per_species, //Output
                abs_xsec_per_species, abs_vmrs, abs_p, abs_t,       //Input
                verbosity);                                        //Verbose!

            Tensor3 part_abs_mat_ZeemanO2((*f_grid_pointer).nelem(), 4, 4), abs_mat_ZeemanO2((*f_grid_pointer).nelem(), 4, 4, 0.0);
            Matrix  K(4,4); Index DM;

            // For pi, sp, sm and defaulting Zeeman transitions
            for(Index ii = 0; ii< 4; ii++)
            {
                // Assign the DM transition type,
                if(ii == 0){DM = 0;}
                else if(ii == 1){DM = 1;}
                else if(ii == 2){DM = -1;}
                else{DM = 1023;} // 1023 is a default, not physical, value.

                // then get the rotation extinction matrix for this DM transition
                K_mat(K, theta*DEG2RAD, eta*DEG2RAD, DM);
                // and multiplied it per frequency dependent absorption coefficient
                mult(part_abs_mat_ZeemanO2, abs_coef_per_species[ii](joker,0), K);
                // to build up the species extinction matrix as a function of frequency.
                abs_mat_ZeemanO2+=part_abs_mat_ZeemanO2;
            };

            // When done, add to the return Tensor4.
            abs_mat_per_species(II, joker, joker, joker) += abs_mat_ZeemanO2;

        }
    }
    else // if the magnetic field is ignored
    {
        for(Index II = 0; II<abs_species.nelem(); II++)
        {
            // If the species isn't Zeeman, look at the next species.
            if(!abs_species[II][0].Zeeman()) continue;

            //Assign temporary variables:
            ArrayOfArrayOfSpeciesTag temp_abs_species(1);
                temp_abs_species[0] = abs_species[II];
             Vector temp_vmrs(1);
                temp_vmrs[0] = rte_vmr_list[II];
             ArrayOfArrayOfLineRecord temp_abs_lines_per_species;
                temp_abs_lines_per_species.push_back(abs_lines_per_species[II]);
                
            // Using the standard scalar absorption functions to get physics parameters,
            Vector abs_p, abs_t; Matrix abs_vmrs;
            AbsInputFromRteScalars(abs_p, abs_t, abs_vmrs, // Output
            rte_pressure, rte_temperature, temp_vmrs,     //Input
            verbosity);                                  // Verbose!
            // and then initialize the cross section per species variable.
            ArrayOfMatrix abs_xsec_per_species;
            abs_xsec_per_speciesInit(abs_xsec_per_species, //Output
            temp_abs_species, *f_grid_pointer, abs_p,     //Input
            verbosity);                                  //Verbose!

            // Cheat to remove h2o dependency,
            Vector abs_h2o(1,-1.0);
            // then calculate the cross section per species
            abs_xsec_per_speciesAddLines(abs_xsec_per_species, //Output
                                         temp_abs_species, *f_grid_pointer, abs_p, abs_t, abs_vmrs, //Input
                                         temp_abs_lines_per_species, abs_lineshape,
                                         isotopologue_ratios, verbosity);

            // and take continua into account.
            // 2012-9-4 Stefan: I don't think this function should add
            //          continua, since they are already added elsewhere.
            //          This function should just do the Zeeman LBL part, as
            //          The name says.
            //
            //            abs_xsec_per_speciesAddConts(abs_xsec_per_species,                                    //Output
            //                                         temp_abs_species, *f_grid_pointer, abs_p, abs_t,        //Input
            //                                         abs_n2, abs_h2o, abs_vmrs,                             //Input
            //                                         abs_cont_names, abs_cont_parameters, abs_cont_models, //Input
            //                                         verbosity);                                          //Verbose!

            // We can finally obtain the coefficients per species!
            Matrix abs_coef; ArrayOfMatrix abs_coef_per_species;
                abs_coefCalcFromXsec(abs_coef, abs_coef_per_species, //Output
                abs_xsec_per_species, abs_vmrs, abs_p, abs_t,       //Input
                verbosity);                                        //Verbose!

            Tensor3 abs_mat_ZeemanO2(abs_mat_per_species.npages(), 4, 4, 0.);
                abs_mat_ZeemanO2(joker,0,0) = abs_coef_per_species[0](joker,0);
                abs_mat_ZeemanO2(joker,1,1) = abs_coef_per_species[0](joker,0);
                abs_mat_ZeemanO2(joker,2,2) = abs_coef_per_species[0](joker,0);
                abs_mat_ZeemanO2(joker,3,3) = abs_coef_per_species[0](joker,0);
            
            abs_mat_per_species(II, joker, joker, joker) = abs_mat_ZeemanO2;
        }
    }
}
