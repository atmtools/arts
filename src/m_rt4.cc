/* Copyright (C) 2016 Jana Mendrok <jana.mendrok@gmail.com>

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
   USA. 
*/
  
/*!
  \file   m_rt4.cc
  \author Jana Mendrok <jana.mendrok@gmail.com>
  \date   2016-05-24
  
  \brief  Workspace functions related to application of scattering solver RT4.
  
  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include "auto_md.h"
#include <complex.h>
#include "messages.h"
#include "m_xml.h"
#include "rt4.h"


#ifdef ENABLE_RT4
/* Workspace method: Doxygen documentation will be auto-generated */
void RT4Test( Tensor4& out_rad,
              const Verbosity& verbosity)
{
/*
    rt4_test( z_file, T_file, abs_gas_file,
              ext_par_file, abs_par_file, sca_par_file,
              verbosity );
*/
    //emissivity.resize(4);
    //reflectivity.resize(4);

    Index nstokes=2;
    Index nummu=8;
    Numeric max_delta_tau=1.0E-6;
    String quad_type="L";
    Numeric ground_temp=300.;
    String ground_type="L";
    Numeric ground_albedo=0.05;
    Complex ground_index;
    Numeric sky_temp=0.;
    Numeric wavelength=880.;
    //Index noutlevels=1;
    //ArrayOfIndex outlevels(1);
    //outlevels[0]=1;

    Vector height, temperatures, gas_extinct;
    Tensor5 sca_data;
    Tensor4 ext_data;
    Tensor3 abs_data;
    ReadXML( height, "height", "testdata/z.xml", "", verbosity );
    ReadXML( temperatures, "temperatures", "testdata/T.xml", "", verbosity );
    ReadXML( gas_extinct, "gas_extinct", "testdata/abs_gas.xml", "", verbosity );
    ReadXML( abs_data, "abs_data", "testdata/abs_par.xml", "", verbosity );
    ReadXML( ext_data, "ext_data", "testdata/ext_par.xml", "", verbosity );
    ReadXML( sca_data, "sca_data", "testdata/sca_par.xml", "", verbosity );
    Index num_layers=height.nelem()-1;
    Index num_scatlayers=3;
    Vector scatlayers(num_layers,0.);
    scatlayers[3]=1.;
    scatlayers[4]=2.;
    scatlayers[5]=3.;

    // the read in sca/ext/abs_data is the complete set (and it's in the wrong
    // order for passing it directly to radtrano). before handing over to
    // fortran, we need to reduce it to the number of stokes elements to be
    // used. we can't use views here as all data needs to be continuous in
    // memory; that is, we have to explicitly copy the data we need.
    Tensor6 scatter_matrix(num_scatlayers,4,nummu,nstokes,nummu,nstokes);
    for( Index ii=0; ii<4; ii++ )
      for( Index ij=0; ij<nummu; ij++ )
        for( Index ik=0; ik<nstokes; ik++ )
          for( Index il=0; il<nummu; il++ )
            for( Index im=0; im<nstokes; im++ )
              for( Index in=0; in<num_scatlayers; in++ )
                scatter_matrix(in,ii,ij,ik,il,im) = sca_data(im,il,ik,ij,ii);
    Tensor5 extinct_matrix(num_scatlayers,2,nummu,nstokes,nstokes);
    for( Index ii=0; ii<2; ii++ )
      for( Index ij=0; ij<nummu; ij++ )
        for( Index ik=0; ik<nstokes; ik++ )
          for( Index il=0; il<nstokes; il++ )
            for( Index im=0; im<num_scatlayers; im++ )
              extinct_matrix(im,ii,ij,ik,il) = ext_data(il,ik,ij,ii);
    Tensor4 emis_vector(num_scatlayers,2,nummu,nstokes);
    for( Index ii=0; ii<2; ii++ )
      for( Index ij=0; ij<nummu; ij++ )
        for( Index ik=0; ik<nstokes; ik++ )
          for( Index il=0; il<num_scatlayers; il++ )
            emis_vector(il,ii,ij,ik) = abs_data(ik,ij,ii);

    // Declaration of Output variables
    Vector mu_values(nummu);
    Tensor3 up_rad(num_layers+1,nummu,nstokes, 0.);
    Tensor3 down_rad(num_layers+1,nummu,nstokes, 0.);

    radtrano_( nstokes,
               nummu,
               max_delta_tau,
               quad_type.c_str(),
               ground_temp,
               ground_type.c_str(),
               ground_albedo,
               ground_index,
               sky_temp,
               wavelength,
               num_layers,
               height.get_c_array(),
               temperatures.get_c_array(),
               gas_extinct.get_c_array(),
               num_scatlayers,
               scatlayers.get_c_array(),
               extinct_matrix.get_c_array(),
               emis_vector.get_c_array(),
               scatter_matrix.get_c_array(),
               //noutlevels,
               //outlevels.get_c_array(),
               mu_values.get_c_array(),
               up_rad.get_c_array(),
               down_rad.get_c_array()
             );

    //so far, output is in
    //    units W/m^2 um sr
    //    dimensions [numlayers+1,nummu,nstokes]
    //    sorted from high to low (altitudes) and 0 to |1| (mu)
    //WriteXML( "ascii", up_rad, "up_rad.xml", 0, "up_rad", "", "", verbosity );
    //WriteXML( "ascii", down_rad, "down_rad.xml", 0, "down_rad", "", "", verbosity );
    
    //to be able to compare with RT4 reference results, reshape output into
    //RT4-output type table (specifically, resort up_rad such that it runs from
    //zenith welling to horizontal, thus forms a continuous angle grid with
    //down_rad. if later changing up_rad/down_rad sorting such that it is in
    //line with doit_i_field, then this has to be adapted as well...
    out_rad.resize(num_layers+1,2,nummu,nstokes);
    for( Index ii=0; ii<nummu; ii++ )
      out_rad(joker,0,ii,joker) = up_rad(joker,nummu-1-ii,joker);
    //out_rad(joker,0,joker,joker) = up_rad;
    out_rad(joker,1,joker,joker) = down_rad;
    //WriteXML( "ascii", out_rad, "out_rad.xml", 0, "out_rad", "", "", verbosity );
}

#else /* ENABLE_RT4 */

void RT4Test( Tensor4&,
              const Verbosity& )
{
  throw runtime_error ("This version of ARTS was compiled without RT4 support.");
}

#endif /* ENABLE_RT4 */

