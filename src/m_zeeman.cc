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

#include "auto_md.h"
#include "zeeman.h"

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
    if( do_zeeman       && (
      ( rtp_mag[0]!=0   || 
        rtp_mag[1]!=0   || 
        rtp_mag[2]!=0 ) || 
        manual_zeeman_tag != 0 ))
    {
      Numeric H_mag,eta,theta;
      set_magnetic_parameters(H_mag,
                              eta,
                              theta,
                              manual_zeeman_tag,
                              manual_zeeman_eta,
                              manual_zeeman_theta,
                              manual_zeeman_magnetic_field_strength,
                              rtp_mag,
                              R_path_los);

            // Reinitialize every loop to empty the set.
        ArrayOfLineRecord temp_abs_lines_sm, temp_abs_lines_sp, //sigma minus, sigma plus
                          temp_abs_lines_pi; // pi
        Numeric (*frequency_change)(const Rational&, const  Rational&, const Rational&, 
                                    const Numeric&, const Index&, const Index&, 
                                    const Index&, const Numeric&, const Numeric&);

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
		    const Index hund = (Index) isotopologue_quantum.getParam(temp_LR.Species(), temp_LR.Isotopologue(), 2);
		    const Numeric GS   = isotopologue_quantum.getParam(temp_LR.Species(), temp_LR.Isotopologue(), 0);
		    Numeric RS_sum     = 0; //Sum relative strength (which ought be close to one by the end)
                    // Only look at lines which have no change in the main rotational number
		    
                    // Separate setting of the frequency_change function...
                    if( hund ==0 )//Case a
                        frequency_change=frequency_change_casea;
                    else if( hund == 1 )// Case b
                        frequency_change=frequency_change_caseb;
                    else
                    {
                        std::ostringstream os;
                        os << "There are undefined Hund cases: " << temp_LR << 
                        "\nThe case is: "<<hund<<", allowed are (a): "<<0<<" and (b): " << 1<<"\n";
                        throw std::runtime_error(os.str());
                    }
                    
                    // Quantum numbers
                    Rational Main,J,NA;
                    Index DMain,DJ,DNA;
                    Numeric S;
                    
                    set_quantum_numbers( Main, DMain, J, DJ, NA, DNA, S,
                                         temp_LR, hund, isotopologue_quantum,
                                         1,//DO_Main,
                                         1,//DO_J,
                                         0);//DO_M)
		    
                    if (!J.isUndefined() != 0 && !Main.isUndefined() != 0 ) // This means the lines are considered erroneous if they fail.
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
                                alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, -1,//DM
                                                  DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                temp_abs_lines_sm.push_back(temp_LR);

                                alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, 0,//DM
                                                  DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                temp_abs_lines_pi.push_back(temp_LR);

                                alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, +1,//DM
                                                  DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                temp_abs_lines_sp.push_back(temp_LR);
                            }
                            else if ( DJ ==  0 )
                            { // Then all DM transitions possible for all M
                                alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, -1,//DM
                                                  DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                temp_abs_lines_sm.push_back(temp_LR);
                                if( ! (M == 0) )
                                {
                                    alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, 0,//DM
                                                      DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                    temp_abs_lines_pi.push_back(temp_LR);
                                }

                                alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, +1,//DM
                                                  DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                temp_abs_lines_sp.push_back(temp_LR);
                            }
                            else if ( DJ == -1 )
                            { // Then certain M results in blocked DM transitions
                                if ( M == -J + DJ && M!=0 )
                                { // Lower limit M only allows DM = 1
                                    alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, +1,//DM
                                                      DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                    temp_abs_lines_sp.push_back(temp_LR);

                                }
                                else if ( M == -J + DJ + 1 && M!=0 )
                                { // Next to lower limit M can only allow DM = 1, 0
                                    alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, +1,//DM
                                                      DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                    temp_abs_lines_sp.push_back(temp_LR);

                                    alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, 0,//DM
                                                  DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                    temp_abs_lines_pi.push_back(temp_LR);
                                }
                                else if ( M ==  J - DJ - 1 && M!=0 )
                                { // Next to upper limit M can only allow DM = 0, -1
                                    alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, 0,//DM
                                                  DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                    temp_abs_lines_pi.push_back(temp_LR);

                                    alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, -1,//DM
                                                      DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                    temp_abs_lines_sm.push_back(temp_LR);
                                }
                                else if ( M == J - DJ && M!=0 )
                                { // Upper limit M only allow DM = -1
                                    alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, -1,//DM
                                                      DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                    temp_abs_lines_sm.push_back(temp_LR);
                                }
                                else if( (-J + DJ + 1) ==  (J - DJ - 1) && M == 0)
                                { // Special case for N=1, J=0, M=0. Only allows DM = 0
                                    alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, 0,//DM
                                                      DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                    temp_abs_lines_pi.push_back(temp_LR);
                                }
                                else
                                { // All DM transitions are possible for these M(s)
                                    alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, +1,//DM
                                                      DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                    temp_abs_lines_sp.push_back(temp_LR);

                                    alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, 0,//DM
                                                  DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                    temp_abs_lines_pi.push_back(temp_LR);

                                    alter_linerecord( temp_LR, RS_sum, abs_lines_per_species[II][ii], Main, M, J, S, DJ, -1,//DM
                                                  DMain, H_mag, GS, frequency_change, 1,1,0);//DO RS,DO_DF,DO_QR
                                    temp_abs_lines_sm.push_back(temp_LR);
                                }
                            }
                            else
                            { // The tests above failed and catastrophe follows
                                std::ostringstream os;
                                os << "There seems to be something wrong with the quantum numbers of at least one line in your *abs_lines*. " <<
                                    "Make sure this is a Zeeman line.\nThe upper quantum numbers are: " << 
                                    temp_LR.QuantumNumbers().Upper() <<
                                    "\nThe lower quantum numbers are: " <<
                                    temp_LR.QuantumNumbers().Lower() <<
                                    "\nThe entire line information: " << temp_LR;
                                throw std::runtime_error(os.str());
                            }
                        }

                        if (abs(RS_sum-1.)>margin) //Reasonable confidence?
                        {
                            std::ostringstream os;
                            os << "The sum of relative strengths is not close to one. This is severly problematic and "
                                "you should look into why this happens.\nIt is currently " << RS_sum 
                                << " with DJ: "<<DJ<<", DMain: "<<DMain<<" for line: "<<
                                temp_LR <<"\n";
                            throw std::runtime_error(os.str());
                        }
                    }
                    else
                    {
                        std::ostringstream os;
                        os << "There are undefined quantum numbers in the line: " << temp_LR 
                        << "\nJ is "<<J<<" and Main is "<<Main<<std::endl;
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
