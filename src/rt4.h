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
  \file   rt4.h
  \author Jana Mendrok <jana.mendrok@gmail.com>
  \date   2016-05-24
  
  \brief  Contains functions related to application of scattering solver RT4.
  
*/

#ifndef rt4_h
#define rt4_h

#ifdef ENABLE_RT4
extern "C" {
#endif

    void radtrano_( const Index&   nstokes,
                    const Index&   nummu,
                    const Numeric& max_delta_tau,
                    const char*    quad_type,
                    const Numeric& ground_temp,
                    const char*    ground_type,
                    const Numeric& ground_albedo,
                    const Complex& ground_index,
                    const Numeric& sky_temp,
                    const Numeric& wavelength,
                    const Index&   num_layers,
                    const Numeric* height,
                    const Numeric* temperatures,
                    const Numeric* gas_extinct,
                    const Index&   num_scatlayers,
                    const Numeric* scatlayers,
                    const Numeric* ext_data,
                    const Numeric* abs_data,
                    const Numeric* sca_data,
                    //const Index&   noutlevels,
                    //const Index*   outlevels,
                    Numeric* mu_values,
                    Numeric* up_rad,
                    Numeric* down_rad
                    );

#ifdef ENABLE_RT4
}
#endif


// Define dummy function that throws a runtime error if ARTS is compiled without
// RT4 support.
#ifndef ENABLE_RT4

void radtrano_( const Index&,
                const Index&,
                const Numeric&,
                const char*,
                const Numeric&,
                const char*,
                const Numeric&,
                const Complex&,
                const Numeric&,
                const Numeric&,
                const Index&,
                const Numeric*,
                const Numeric*,
                const Numeric*,
                const Index&,
                const Numeric*,
                const Numeric*,
                const Numeric*,
                const Numeric*,
                //const Index&,
                //const Index*,
                Numeric*,
                Numeric*,
                Numeric*)
{
    throw std::runtime_error("This version of ARTS was compiled without RT4 support.");
}

#endif

/*
void rt4_test(// Output:
              // Input:
              const String& z_file,
              const String& T_file,
              const String& abs_gas_file,
              const String& ext_par_file,
              const String& abs_par_file,
              const String& sca_par_file,
              const Verbosity& verbosity
              );
*/

#endif /* rt4_h */

