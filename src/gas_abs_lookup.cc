/* Copyright (C) 2002,2003 Stefan Buehler <sbuehler@uni-bremen.de>

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

/*!
  \file   gas_abs_lookup.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Thu Sep 19 17:25:07 2002
  
  \brief  Implementation of scalar gas absorption lookup table functions.
*/

#include <cmath>
#include "gas_abs_lookup.h"
#include "interpolation.h"
#include "make_vector.h"
#include "logic.h"
#include "check_input.h"
#include "messages.h"
#include "physics_funcs.h"

//! Find positions of new grid points in old grid.
/*! 
  Uses gridpos to do most of the work.

  Comparison of Numerics is a bit tricky.

  \retval pos      Positions of new grid points in old grid.
  \param  old_grid The old grid.
  \param  new_grid The new grid.
*/
void find_new_grid_in_old_grid( ArrayOfIndex& pos,
                                ConstVectorView old_grid,
                                ConstVectorView new_grid )
{
  const Index n_new_grid = new_grid.nelem();

  // Make sure that pos has the right size:
  assert( n_new_grid == pos.nelem() );
  
  // We can use gridpos to do most of the work!
  ArrayOfGridPos gp( n_new_grid );
  gridpos( gp, old_grid, new_grid );

  // Convert to approximate numerical indices
  Vector approx_pos( n_new_grid );
  for ( Index i=0; i<n_new_grid; ++i )
    {
      const GridPos& tgp = gp[i];
      approx_pos[i] = tgp.idx + tgp.fd[0];
      // The last term is necessary, because we could in fact be
      // almost on the next grid point.
    }

  // We now have to check if the grid positions stored in gp are
  // sufficiently close to the grid points for our taste.
  for ( Index i=0; i<n_new_grid; ++i )
    {
      out3 << "  " << new_grid[i] << ": ";
      // This is the crucial if statement for the comparison of two
      // numerics!
      Numeric diff = approx_pos[i] - rint(approx_pos[i]);
      if ( 0 != diff )
        {
          ostringstream os;
          os << "Found no match for element [" << i << "] of the new grid.\n"
             << "Value: " << new_grid[i] << "\n"
             << "Diff:  " << diff;

          throw runtime_error( os.str() );          
        }
      else
        {
          // Assign to the output array:
          pos[i] = (Index) rint(approx_pos[i]);
          out3 << "found, index = " << pos[i] << ".\n";
        }
    }
}

//! Adapt lookup table to current calculation.
/*!
  This method has the following tasks:

  1. Find and remember the indices of the current species in the
  lookup table. At the same time verify that each species is included
  in the table exactly once.

  2. Find and remember the frequencies of the current calculation in
  the lookup table. At the same time verify that all frequencies are
  included and that no frequency occurs twice.

  3. Use the species and frequency index lists to build the new lookup
  table.

  4. Replace original table by the new one.

  5. Initialize log_p_grid.

  The method is intended to be called only once per ARTS job, more or
  less directly from a corresponding workspace method. Therefore,
  runtime errors are thrown, rather than assertions, if something is
  wrong. 

  \param current_species The list of species for the current calculation.
  \param current_f_grid  The list of frequencies for the current calculation.

  \date 2002-12-12
*/
void GasAbsLookup::Adapt( const ArrayOfArrayOfSpeciesTag& current_species,
                          ConstVectorView current_f_grid )
{
  // Some constants we will need:
  Index n_current_species = current_species.nelem();
  Index n_current_f_grid  = current_f_grid.nelem();

  Index n_species         = species.nelem();
  Index n_f_grid          = f_grid.nelem();
  Index n_p_grid          = p_grid.nelem();
  
  out2 << "  Original table: " << n_species << " species, "
       << n_f_grid << " frequencies.\n"
       << "  Adapt to:       " << n_current_species << " species, "
       << n_current_f_grid << " frequencies.\n";

  if ( 0 == nonlinear_species.nelem() )
    {
      out2 << "  Table contains no nonlinear species.\n";
    }

  if ( 0 == t_pert.nelem() )
    {
      out2 << "  Table contains no temperature perturbations.\n";
    }

  // We are constructing a new lookup table, containing just the
  // species and frequencies that are necessary for the current
  // calculation. We will build it in this local variable, then copy
  // it back to *this in the end.
  GasAbsLookup new_table;

  // First some checks on the lookup table itself:

  // Species:
  if ( 0 == n_species )
    {
      ostringstream os;
      os << "The lookup table should have at least one species.";
      throw runtime_error( os.str() );
    }
  
  // Nonlinear species:
  // They should be monotonically increasing and pointing at valid
  // species. 
  {
    Index max = -1;
    for ( Index i=0; i<nonlinear_species.nelem(); ++i )
      {
        ostringstream os;
        os << "nonlinear_species[" << i << "]";
        chk_if_in_range( os.str(),
                         nonlinear_species[i],
                         0,
                         n_species-1 );

        if ( max >= nonlinear_species[i] )
        {
          os << "The array of indices *nonlinear_species* should\n"
             << "be monotonically increasing.";
            throw runtime_error( os.str() );
        }

        max = nonlinear_species[i];
      }
  }

  // Frequency grid:
  chk_if_increasing( "f_grid", f_grid );

  // Pressure grid:
  chk_if_decreasing( "p_grid", p_grid );

  // Reference VMRs:
  chk_matrix_nrows( "vmrs_ref", vmrs_ref, n_species );
  chk_matrix_ncols( "vmrs_ref", vmrs_ref, n_p_grid );

  // Reference temperatur:
  chk_vector_length( "t_ref", t_ref, n_p_grid );

  // Temperature perturbations:
  // Nothing to check for t_pert, it seems.

  // Perturbations for nonlinear species:
  // Check that nls_pert is empty if and only if nonlinear_species is
  // empty:
  if ( 0 == nonlinear_species.nelem() )
    {
      chk_vector_length( "nls_pert", nls_pert, 0 );
    }
  else
    {
      if ( 0 == nls_pert.nelem() )
        {
          ostringstream os;
          os << "The vector nls_pert should contain the perturbations\n"
             << "for the nonlinear species, but it is empty.";
          throw runtime_error( os.str() );
        }
    }

  // The table itself, xsec:
  //
  // We have to separtely consider the 3 cases described in the
  // documentation of GasAbsLookup.
  //
  //     Dimension: [ a, b, c, d ]
  //
  if ( 0==nonlinear_species.nelem() )
    {
      if ( 0==t_pert.nelem() )
        {
          //     Simplest case (no temperature perturbations,
          //     no vmr perturbations):
          //     a = 1
          //     b = n_species
          //     c = n_f_grid 
          //     d = n_p_grid
          chk_size( "xsec", xsec,
                    1,
                    n_species,
                    n_f_grid,
                    n_p_grid );
        }
      else
        {
          //     Standard case (temperature perturbations,
          //     but no vmr perturbations):
          //     a = n_t_pert
          //     b = n_species
          //     c = n_f_grid
          //     d = n_p_grid
          chk_size( "xsec", xsec,
                    t_pert.nelem(),
                    n_species,
                    n_f_grid,
                    n_p_grid );
        }
    }
  else
    {
      //     Full case (with temperature perturbations and
      //     vmr perturbations):
      //     a = n_t_pert
      //     b = n_species + n_nonlinear_species * ( n_nls_pert - 1 )
      //     c = n_f_grid
      //     d = n_p_grid
      Index a = t_pert.nelem();
      Index b = n_species
        + nonlinear_species.nelem()
        * ( nls_pert.nelem() - 1 );
      Index c = n_f_grid;
      Index d = n_p_grid;

      chk_size( "xsec", xsec, a, b, c, d );
    }


  // Now some checks on the input data:

  // The list of current species should not be empty:
  if ( 0==n_current_species )
    {
      ostringstream os;
      os << "The list of current species should not be empty.";
      throw runtime_error( os.str() );
    }

  // The grid of current frequencies should be monotonically sorted:
  chk_if_increasing( "current_f_grid", current_f_grid );


  // 1. Find and remember the indices of the current species in the
  //    lookup table. At the same time verify that each species is
  //    included in the table exactly once.
  ArrayOfIndex i_current_species(n_current_species);
  out3 << "  Looking for species in lookup table:\n";
  for ( Index i=0; i<n_current_species; ++i )
    {
      out3 << "  " << get_tag_group_name( current_species[i] ) << ": ";
      // We need no error checking for the next statement, since the
      // chk_contains function throws a runtime error if the species
      // is not found exactly once.
      i_current_species[i] =
        chk_contains( "species", species, current_species[i] );
      out3 << "found.\n";
    }


  // 2. Find and remember the frequencies of the current calculation in
  //    the lookup table. At the same time verify that all frequencies are
  //    included and that no frequency occurs twice.

  // FIXME: This is a bit tricky, because we are comparing
  // Numerics. Let's see how well this works in practice.

  ArrayOfIndex i_current_f_grid(n_current_f_grid);
  out3 << "  Looking for Frequencies in lookup table:\n";

  // We need no error checking for the next statement, since the
  // function called throws a runtime error if a frequency
  // is not found, or if the grids are not ok.
  find_new_grid_in_old_grid( i_current_f_grid,
                             f_grid,
                             current_f_grid );


  // 3. Use the species and frequency index lists to build the new lookup
  // table.

  // Species:
  new_table.species.resize( n_current_species );
  for ( Index i=0; i<n_current_species; ++i )
    {
      new_table.species[i] = species[i_current_species[i]];

      // Is this a nonlinear species?
      if ( 0 <= find_first( nonlinear_species,
                            i_current_species[i] ) )
        {
          new_table.nonlinear_species.push_back( i );
        }
    }

  // Frequency grid:
  new_table.f_grid.resize( n_current_f_grid );
  for ( Index i=0; i<n_current_f_grid; ++i )
    {
      new_table.f_grid[i] = f_grid[i_current_f_grid[i]];
    }

  // Pressure grid:
  new_table.p_grid.resize( n_p_grid );
  new_table.p_grid = p_grid;

  // Reference VMR profiles:
  new_table.vmrs_ref.resize( n_current_species,
                             n_p_grid     );
  for ( Index i=0; i<n_current_species; ++i )
    {
      new_table.vmrs_ref( i,
                          Range(joker) )
        = vmrs_ref( i_current_species[i],
                    Range(joker)          );
    }

  // Reference temperature profile:
  new_table.t_ref.resize( t_ref.nelem() );
  new_table.t_ref = t_ref;

  // Vector of temperature perturbations:
  new_table.t_pert.resize( t_pert.nelem() );
  new_table.t_pert = t_pert;

  // Vector of perturbations for the VMRs of the nonlinear species: 
  // (Should stay empty if we have no nonlinear species)
  if ( 0 != new_table.nonlinear_species.nelem() )
    {
      new_table.nls_pert.resize( nls_pert.nelem() );
      new_table.nls_pert = nls_pert;
    }

  // Absorption coefficients:
  new_table.xsec.resize( xsec.nbooks(),
                         n_current_species,
                         n_current_f_grid,
                         xsec.ncols()        );

  // We have to copy the right species and frequencies from the old to
  // the new table. Temperature perturbations and pressure grid remain
  // the same.

  // Do species:
  for ( Index i_s=0; i_s<n_current_species; ++i_s )
    {
      // Do frequencies:
      for ( Index i_f=0; i_f<n_current_f_grid; ++i_f )
        {
          new_table.xsec( Range(joker),
                         i_s,
                         i_f,
                         Range(joker) )
            =
            xsec( Range(joker),
                 i_current_species[i_s],
                 i_current_f_grid[i_f],
                 Range(joker) );
        }
    }


  // 4. Replace original table by the new one.
  *this = new_table;

  // Obsolete!
  //   // 5. Initialize log_p_grid.
  //   log_p_grid.resize( n_p_grid );
  //   transform( log_p_grid,
  //              log10,
  //              p_grid );
}

//! Extract scalar gas absorption coefficients from the lookup table. 
/*!  
  This carries out a simple interpolation in temperature and
  pressure. The interpolated value is then scaled by the ratio between
  actual VMR and reference VMR. In the case of nonlinear species the
  interpolation goes also over VMR.

  All input parameters (f_index, p, T) must be in the range coverd by
  the table.

  In this case pressure is not an altitude coordinate, so we are free
  to choose the type of interpolation that gives lowest interpolation
  errors or is easiest. I tested both linear and log p interpolation
  with the result that it makes no difference. Therefore, linear
  interpolation is used.

  \retval sga A Matrix with scalar gas absorption coefficients
  [1/m]. Dimension is adjusted automatically to either
  [1,n_species] or [n_f_grid,n_species]!

  \param f_index The frequency index. If this is >=0, it means that
  absorption for this frequency will be extracted. (The leading
  dimension of sga will be 1.) If this is <0, it means that absorption
  for ALL frequencies is extracted. (The leading dimension of sga will
  be n_f_grid.)

  \param p The pressures [Pa].

  \param T The temperature [K].

  \param vmrs The VMRs [absolute number]. Dimension: [species].  

  \date 2002-09-20, 2003-02-22
*/
void GasAbsLookup::Extract( Matrix&         sga,
                            const Index&    f_index,
                            const Numeric&  p,
                            const Numeric&  T,
                            ConstVectorView vmrs ) const
{
  // Obtain some properties of the lookup table:
  
  // Number of gas species in the table:
  const Index n_species = species.nelem();

  // Number of nonlinear species:
  const Index n_nonlinear_species = nonlinear_species.nelem();

  // Number of frequencies in the table:
  const Index n_f_grid = f_grid.nelem();

  // Number of pressure grid points in the table:
  const Index n_p_grid = p_grid.nelem();

  // Number of temperature perturbations:
  const Index n_t_pert = t_pert.nelem();

  // Number of nonlinear species perturbations:
  //  const Index n_nls_pert = nls_pert.nelem();


  // First some checks on the lookup table itself:

  // Obsolete! FIXME: Do I stick with this?
  // Assert that log_p_grid has been initialized:
  //  assert( is_size( log_p_grid, n_p_grid) );

  // Check that the dimension of vmrs_ref is consistent with species and p_grid:
  assert( is_size( vmrs_ref, n_species, n_p_grid) );

  // Check dimension of t_ref:
  assert( is_size( t_ref, n_p_grid) );


  // Following are some checks on the input variables:

  // Assert that vmrs has the right dimension:
  assert( is_size( vmrs, n_species ) );

  // We need the log10 of the pressure:
  //  const Numeric log_p = log10( p );

  // We also set the start and extent for the frequency loop.
  Index f_start, f_extent;

  if ( f_index < 0 )
    {
      // This means we should extract for all frequencies.

      f_start  = 0;
      f_extent = n_f_grid;
    }
  else
    {
      // This means we should extract only for one frequency.

      // Check that f_index is inside f_grid:
      assert( f_index < n_f_grid );

      f_start  = f_index;
      f_extent = 1;
    }

  // Adjust size of sga accordingly:
  sga.resize(f_extent, n_species);


  // Now we will start to do some business.

  // Calculate the number density for the given pressure and
  // temperature: 
  // n = n0*T0/p0 * p/T or n = p/kB/t, ideal gas law
  const Numeric n = number_density( p, T );

  // For sure, we need to store the pressure grid position. 
  GridPos pgp;
  gridpos( pgp,
           p_grid,
           p );

  // This bit is obsolete, we do linear interpolation in p!
  //   Vector log_p_grid( n_p_grid );
  //   transform( log_p_grid,
  //              log10,
  //              p_grid );
  //   gridpos( pgp,
  //            log_p_grid,
  //            log10( p ) );

  // Pressure interpolation weights:
  Vector pitw(2);
  interpweights(pitw,pgp);


  // Let's first treat the case with no nonlinear species.
  // This case is marked by n_nonlinear_species == 0.
  if ( 0 == n_nonlinear_species )
    {
      // In this case, n_nls_pert must also be zero:
      assert( is_size( nls_pert, 0 ) );
      
      // The simplest case is that the table contains neither temperature
      // nor VMR perturbations. This means it is not really a
      // lookup table, just an absorption profile. In this case we ignore
      // the temperature, and interpolate in pressure only.
      // This case is marked by t_pert being an empty vector, since we
      // already know that there are no non-linear species.
      if ( 0 == n_t_pert )
        {
          // Verify, that abs has the right dimensions for this case:
          assert( is_size( xsec,
                           1,          // temperature perturbations
                           n_species,  // species
                           n_f_grid,   // frequency grid points
                           n_p_grid    // pressure levels
                           ) );

          // Loop over frequency:
          for ( Index s=0; s<f_extent; ++s )
            {
              // Loop over species:
              for ( Index i=0; i<n_species; ++i )
                {
                  // Get the right view on xsec. (Only a vector of
                  // pressure for this particular species and
                  // frequency):
                  ConstVectorView this_xsec = xsec( 0,         // T
                                                    i,         // species
                                                    f_start+s, // frequency
                                                    joker      // p
                                                    );

                  // Get the right view on our result variable,
                  // sga. (Only a scalar):
                  Numeric& this_sga = sga(s,i);

                  // Do the interpolation:
                  this_sga = interp( pitw,
                                     this_xsec,
                                     pgp );

                  // Watch out, this is not yet the final result, we
                  // need to multiply with the number density n:
                  this_sga *= ( n * vmrs[i] );
                }
            }
        }
      else
        {
          // This is the case with temperature variations, but without
          // nonlinear species.
          //
          // This means we have to do a simultaneous interpolation in
          // pressure and temperature.
          
          // Verify, that xsec has the right dimensions for this case:
          assert( is_size( xsec,
                           n_t_pert,   // temperature perturbations
                           n_species,  // species  
                           n_f_grid,   // frequency grid points
                           n_p_grid    // pressure levels      
                           ) );

          // The tricky bit is that we do not have a uniform grid in
          // temperature, but a different temperature grid for each
          // pressure in the table. Therefore, we have to do the
          // bi-linear interpolation in two steps:
          // 1. Interpolate to the correct temperature for the two
          //    neighboring pressure levels.
          // 2. Interpolate in pressure.

          // We need to create the temperature grid for each pressure
          // level, it is not stored directly in the table:
          Vector t_grid(n_t_pert);

          // To store interpolation weights:
          Vector titw(2);

          // To store the temperature interpolated result for the 2
          // pressure levels:
          Tensor3 xsec_T_interpolated( 2, f_extent, n_species );

          // The two lookup table pressure levels in question are
          // pgp.idx and pgp.idx+1
          for ( Index pi=0; pi<2; ++pi )
            {
              t_grid = t_pert;
              // Now t_grid contains just the perturbations. We
              // need to add the reference value for the given
              // pressure:
              t_grid += t_ref[pi+pgp.idx];

              // Temperature grid position:
              GridPos tgp;       // only a scalar
              gridpos( tgp, t_grid, T );

              // Temperature interpolation weights:
              interpweights(titw,tgp);

              // Loop over frequency:
              for ( Index s=0; s<f_extent; ++s )
                {
                  // Loop over species:
                  for ( Index i=0; i<n_species; ++i )
                    {
                      // Get the right view on xsec. (Only a vector of
                      // temperature for this particular pressure
                      // level, species and frequency):
                      ConstVectorView this_xsec = xsec( joker,      // T
                                                        i,          // species
                                                        f_start+s,  // frequency
                                                        pi+pgp.idx  // p
                                                        );

                      // Get the right view on our result variable,
                      // xsec_T_interpolated. (Only a scalar):
                      Numeric& this_xsec_T_interpolated =
                        xsec_T_interpolated( pi, s, i );

                      // Do the interpolation:
                      this_xsec_T_interpolated = interp( titw,
                                                         this_xsec,
                                                         tgp );
                    }
                }
            }

          // Now we have to interpolate between the two pressure levels

          // We can reuse the pressure grid position, but we have to
          // adjust the base index, since we have only two pressure
          // levels now.
          pgp.idx = 0;

          // Loop over frequency:
          for ( Index s=0; s<f_extent; ++s )
            {
              // Loop over species:
              for ( Index i=0; i<n_species; ++i )
                {
                  // Get the right view on
                  // xsec_T_interpolated (a vector for the
                  // two pressures):
                  VectorView this_xsec =
                    xsec_T_interpolated( joker, s, i );

                  // Get the right view on our result variable,
                  // sga. (Only a scalar):
                  Numeric& this_sga = sga(s,i);

                  // Do the interpolation:
                  this_sga = interp( pitw,
                                     this_xsec,
                                     pgp );

                  // Watch out, this is not yet the final result, we
                  // need to multiply with the number density n:
                  this_sga *= ( n * vmrs[i] );
                }
            }
        }
    }
  else
    {
      // So, we *do* have nonlinear species.
      throw(runtime_error("This case is not yet implemented"));
    }
}

// Below is the version of extract for the whole atmospheric field
// that I started, but never finished. 

// //! Extract scalar gas absorption coefficients from the lookup table. 
// /*! 
//   This carries out a simple interpolation in temperature and
//   pressure. The interpolated value is then scaled by the ratio between
//   actual VMR and reference VMR. In the case of nonlinear species the
//   interpolation goes also over VMR.

//   All input parameters (f_index, p, T) must be in the range coverd by
//   the table.

//   FIXME: Should we interpolate linearly in pressure or in
//   log-pressure? In this case pressure is not an altitude coordinate!
//   Anyway, I'll use log-pressure for now.

//   This can extract for a whole bunch of atmospheric conditions
//   simultaneously, which is much more efficient than doing it one by
//   one. Input variables p, T, and vmrs must have consistent dimensions.

//   \retval sga A Tensor5 with scalar gas absorption coefficients
//   [1/m]. Dimension must be either [1, n_species, p_grid, lat_grid,
//   lon_grid] or [n_f_grid, n_species, p_grid, lat_grid, lon_grid]!

//   \param f_index The frequency index. If this is >=0, it means that
//   absorption for this frequency will be extracted. (The relevant
//   dimension of sga will be 1.) If this is <0, it means that absorption
//   for ALL frequencies is extracted. (The relevant dimension of sga
//   must then be n_f_grid.)

//   \param p The pressures [Pa]. Dimension: [p_grid].

//   \param T The temperaturea [K]. Dimension: [p_grid, lat_grid, lon_grid].

//   \param vmrs The VMRs [absolute number]. Dimension: [species, p_grid,
//   lat_grid, lon_grid].  

//   \date 2002-09-20, 2003-02-22
// */
// void GasAbsLookup::Extract( Tensor5View      sga,
//                             const Index&     f_index,
//                             ConstVectorView  p,
//                             ConstTensor3View T,
//                             ConstTensor4iew  vmrs ) const
// {

//   // Obtain some properties of the lookup table:
  
//   // Number of gas species in the table:
//   const Index n_species = species.nelem();

//   // Number of nonlinear species:
//   const Index n_nonlinear_species = nonlinear_species.nelem();

//   // Number of frequencies in the table:
//   const Index n_f_grid = f_grid.nelem();

//   // Number of pressure grid points in the table:
//   const Index n_p_grid = p_grid.nelem();

//   // Number of temperature perturbations:
//   const Index n_t_pert = t_pert.nelem();

//   // Number of nonlinear species perturbations:
//   //  const Index n_nls_pert = nls_pert.nelem();


//   // First some checks on the lookup table itself:

//   // Assert that log_p_grid has been initialized:
//   assert( is_size( log_p_grid, n_p_grid) );

//   // Check that the dimension of vmrs_ref is consistent with species and p_grid:
//   assert( is_size( vmrs_ref, n_species, n_p_grid) );

//   // Check dimension of t_ref:
//   assert( is_size( t_ref, n_p_grid) );


//   // Following are some checks on the input variables:

//   // Number of pressures for which we want to extract:
//   const Index n_p = p.nelem();

//   // Number of latitudes for which we want to extract:
//   const Index n_lat = T.nrows();

//   // Number of longitudes for which we want to extract:
//   const Index n_lon = T.ncols();

//   // Assert, that first dimension of T is consistent with p:
//   assert( n_p == T.npages() );

//   // Assert that vmrs has the right dimension:
//   assert( is_size( vmrs, n_species, n_p, n_lat, n_lon ) );

//   // We need the log10 of the pressure:
//   Vector log_p(n_p);
//   transform( log_p, log10, p);


//   // We also set the start and extent for the frequency loop.
//   Index f_start, f_extent;

//   if ( f_index < 0 )
//     {
//       // This means we should extract for all frequencies.

//       // Assert size of sga accordingly:
//       assert(is_size(sga, n_p, n_f_grid, n_species));
      
//       f_start  = 0;
//       f_extent = n_f_grid;
//     }
//   else
//     {
//       // This means we should extract only for one frequency.

//       // Assert size of sga accordingly:
//       assert(is_size(sga, n_p, 1, n_species));

//       // Check that f_index is inside f_grid:
//       assert( f_index < n_f_grid );

//       f_start  = f_index;
//       f_extent = 1;
//     }


//   // Now we will start to do some business.

//   // For sure, we need to store the pressure grid positions. Remember
//   // that we use log interpolation here. (FIXME: Discuss with Patrick.)
//   ArrayOfGridPos pgp(n_p);
//   gridpos( pgp,
//            log_p_grid,
//            log_p );

//   // Let's first treat the case with no nonlinear species.
//   // This case is marked by n_nonlinear_species == 0.
//   if ( 0 == n_nonlinear_species )
//     {
//       // In this case, n_nls_pert must also be zero:
//       assert( is_size( nls_pert, 0 ) );
      
//       // The simplest case is that the table contains neither temperature
//       // nor VMR perturbations. This means it is not really a
//       // lookup table, just an absorption profile. In this case we ignore
//       // the temperature, and interpolate in pressure only.
//       // This case is marked by t_pert being an empty vector, since we
//       // already know that there are no non-linear species.
//       if ( 0 == n_t_pert )
//         {
//           // Verify, that abs has the right dimensions for this case:
//           assert( is_size( abs,
//                            1,          // temperature perturbations
//                            n_species,  // species  
//                            n_f_grid,   // frequency grid points
//                            n_p_grid    // pressure levels      
//                            ) );
          
//           // To store interpolation weights:
//           Matrix itw(n_p,2);
//           interpweights(itw,pgp);

//           // To temporarily store interpolation result:
//           Vector sga_interp(n_p);

//           // VMR scaling factors:
//           Vector vmr_scaling(n_p);

//           // Loop over frequency:
//           for ( Index s=f_start; s<f_extent; ++s )
//             {
//               // Loop over species:
//               for ( Index i=0; i<n_species; ++i )
//                 {
//                   // Get the right view on abs. (Only a vector of
//                   // pressure for this particular species and
//                   // frequency):
//                   ConstVectorView this_abs = abs( 0,      // T
//                                                   i,      // species
//                                                   s,      // frequency
//                                                   joker   // p
//                                                   );

//                   // Do the interpolation:
//                   interp( sga_interp,
//                           itw,
//                           this_abs,
//                           pgp );

//                   // Copy this result to all latitudes and
//                   // longitudes. This is ok, since we have no T
//                   // dependence in the table. Don't forget to scale
//                   // with the VMR!
//                   for ( j=0; j<n_lat; ++j )
//                     for ( k=0; k<n_lon; ++k )
//                       {
//                         // Get the right view on our result variable,
//                         // sga. (Only a vector of pressure.)
//                         VectorView this_sga = sga( s,         // frequency
//                                                    i,         // species
//                                                    joker,     // p
//                                                    j,         // lat
//                                                    k          // lon
//                                                    );
//                         this_sga = sga_interp;
//                         // Watch out, this is not yet the final
//                         // result, we need to scale with the actual
//                         // VMR values.

//                         // Determine VMR scaling factors:
//                         vmr_scaling = vmrs( i,            // species
//                                             joker,        // p
//                                             j,            // lat
//                                             k             // lon
//                                             );
//                         vmr_scaling /= vmrs_ref( i,       // species
//                                                  joker    // p
//                                                  );

//                         // Scale this_vmr with the VMR scaling factors:
//                         this_sga *= vmr_scaling;
//                       }

//                 }
//             }
//         }
//       else
//         {
//           // This is the case with temperature variations, but without nonlinear species. 
          
//           // This means we have to do a simultaneous interpolation in pressure and temperature. 

//           // Verify, that abs has the right dimensions for this case:
//           assert( is_size( abs,
//                            n_t_pert,   // temperature perturbations
//                            n_species,  // species  
//                            n_f_grid,   // frequency grid points
//                            n_p_grid    // pressure levels      
//                            ) );

//           // We need to create the temperature grid for each pressure
//           // level, it is not stored directly in the table:
//           Vector t_grid(n_t_pert);

//           // To store interpolation weights:
//           Vector itw(4);

//           // Loop pressures:
//           for ( m=0; m<n_p; ++m )
//             {
//               t_grid = t_pert;
//               // Now t_grid contains just the perturbations. We
//               // need to add the reference value for the given
//               // pressure:
//               t_grid += t_ref[m];

//               // Loop latitudes:
//               for ( j=0; j<n_lat; ++j )
//                 {
//                   // Loop longitudes:
//                   for ( k=0; k<n_lon; ++k )
//                     {
//                       // We already have the pressure grip positions in pgp, but
//                       // now we also need the temperature grid position:
//                       GridPos tgp;       // only a scalar 
          
//                       // Calculate T grid position:
//                       gridpos( tgp,
//                                t_grid,
//                                T( m,     // pressure
//                                   j,     // latitude
//                                   k      // longitude
//                                   )
//                                );

//                       interpweights( itw, tgp, pgp[m] );

//                       // Loop over frequency:
//                       for ( Index s=f_start; s<f_extent; ++s )
//                         {
//                           // Loop over species:
//                           for ( Index i=0; i<n_species; ++i )
//                             {
//                               // Get view on abs:
//                               ConstMatrixView this_abs = abs( joker,  // T
//                                                               i,      // species
//                                                               s,      // frequency
//                                                               joker   // p
//                                                               );
//                               // Get view on sga:
//                               Numeric& this_sga = sga( s,         // frequency
//                                                        i,         // species
//                                                        m,         // p
//                                                        j,         // lat
//                                                        k          // lon
//                                                        );

//                               // Interpolate:
//                               interp( this_sga,
//                                       itw,
//                                       this_abs,
//                                       tgp,
//                                       pgp );

//                               // Scale with VMR:
//                               this_sga *= vmrs( i,            // species
//                                                 m,            // p
//                                                 j,            // lat
//                                                 k             // lon
//                                                 );
//                               this_sga /= vmrs_ref( i,        // species
//                                                     m         // p
//                                                     );
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }
//   else
//     {
//       // So, we *do* have nonlinear species.
//       throw(runtime_error("This case is not yet implemented"));
//     }
// }

