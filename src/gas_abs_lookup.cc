/* Copyright (C) 2002 Stefan Buehler <sbuehler@uni-bremen.de>

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


//! Find positions of new grid points in old grid.
/*! 
  Uses gridpos to do most of the work.

  Comparison of Numerics is a bit tricky.

  \retval pos      Positions of new grid points in old grid.
  \param  old_grid The old grid.
  \param  new_grid The new grid.
*/
void find_new_grid_in_old_grid( ArrayOfIndex pos,
                                ConstVectorView old_grid,
                                ConstVectorView new_grid )
{
  Index n_new_grid = new_grid.nelem();
  
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
      // This is the crucial if statement for the comparison of two
      // numerics!
      Numeric diff = approx_pos[i] - floor(approx_pos[i]);
      if ( 0 != diff )
        {
          ostringstream os;
          os << "Found no match for element [" << i << "] of the new grid.\n"
             << "Value: " << new_grid[i] << "\n"
             << "Diff:  " << diff;

          throw runtime_error( os.str() );          
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

  // The table itself, abs:
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
          chk_size( "abs", abs,
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
          chk_size( "abs", abs,
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

      chk_size( "abs", abs, a, b, c, d );
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
      out3 << "  " << get_tag_group_name( current_species[i] );
      // We need no error checking for the next statement, since the
      // chk_contains function throws a runtime error if the species
      // is not found exactly once.
      i_current_species[i] =
        chk_contains( "species", species, current_species[i] );
      out3 << ": found.\n";
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
  new_table.abs.resize( abs.nbooks(),
                        n_current_species,
                        n_current_f_grid,
                        abs.ncols()        );

  // We have to copy the right species and frequencies from the old to
  // the new table. Temperature perturbations and pressure grid remain
  // the same.

  // Do species:
  for ( Index i_s=0; i_s<n_current_species; ++i_s )
    {
      // Do frequencies:
      for ( Index i_f=0; i_f<n_current_f_grid; ++i_f )
        {
          new_table.abs( Range(joker),
                         i_s,
                         i_f,
                         Range(joker) )
            =
            abs( Range(joker),
                 i_current_species[i_s],
                 i_current_f_grid[i_f],
                 Range(joker) );
        }
    }


  // 4. Replace original table by the new one.
  *this = new_table;


  // 5. Initialize log_p_grid.
  log_p_grid.resize( n_p_grid );
  transform( log_p_grid,
             log10,
             p_grid );
}

//! Extract scalar gas absorption coefficients from the lookup table. 
/*! 
  This carries out a simple interpolation in temperature and
  pressure. The interpolated value is then scaled by the ratio between
  actual VMR and reference VMR. In the case of nonlinear species the
  interpolation goes also over VMR.

  All input parameters (f_index, p, T) must be in the range coverd by
  the table.

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

  \param vmrs The VMRs [absolute number]. Dimension: [n_species]. 
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

  // Assert that log_p_grid has been initialized:
  assert( is_size( log_p_grid, n_p_grid) );

  // Check that the dimension of vmrs_ref is consistent with species and p_grid:
  assert( is_size( vmrs_ref, n_species, n_p_grid) );

  // Check dimension of t_ref:
  assert( is_size( t_ref, n_p_grid) );


  // Following are some checks on the input variables:

  // We also set the start and extent for the frequency loop.
  Index f_start, f_extent;

  if ( f_index < 0 )
    {
      // This means we should extract for all frequencies.

      // Adjust size of sga accordingly:
      sga.resize(n_f_grid, n_species);
      
      f_start  = 0;
      f_extent = n_f_grid;
    }
  else
    {
      // This means we should extract only for one frequency.

      // Adjust size of sga accordingly:
      sga.resize(1, n_species);

      // Check that f_index is inside f_grid:
      assert( f_index < n_f_grid );

      f_start  = f_index;
      f_extent = 1;
    }

  // Assert that vmrs has the right dimension:
  assert( is_size( vmrs, n_species ) );

  // We need the log10 of the pressure:
  const Numeric log_p = log10( p );


  // Now we will start to do some business.

  // For sure, we need to store the vertical grid positions:
  // We need only one, because we want to extract only for a single
  // atmospheric condition.
  ArrayOfGridPos vgp(1);
  gridpos( vgp,
           log_p_grid,
           log_p );

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
          
          assert( is_size( abs,
                           1,          // temperature perturbations
                           n_species,  // species  
                           n_f_grid,   // frequency grid points
                           n_p_grid    // pressure levels      
                           ) );
          
          // To store interpolation weights:
          Matrix itw(1,2);
          interpweights(itw,vgp);

          // Loop over frequency:
          for ( Index s=f_start; s<f_extent; ++s )
            {
              // Loop over species:
              for ( Index i=0; i<n_species; ++i )
                {
                  // Get the right view on abs. (Only a vector of
                  // pressure for this particular species and
                  // frequency):
                  ConstVectorView this_abs = abs( 0,
                                                  i,
                                                  s, 
                                                  Range(joker) );

                  // Get the right view on our result variable,
                  // sga. (Only a scalar):
                  Numeric& this_sga = sga(s,i);

                  // Do the interpolation:
                  interp( this_sga,
                          itw,
                          this_abs,
                          vgp );
                }
            }
        }
    }
}

