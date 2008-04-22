/* Copyright (C) 2002-2007 Stefan Buehler <sbuehler@ltu.se>
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

/*!
  \file   m_abs_lookup.cc
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Wed Nov 20 18:04:20 2002
  
  \brief  Methods related to absorption, lookup table, etc.
*/

#include <algorithm> 
#include <map>
#include <limits>

#include "auto_md.h"
#include "arts.h"
#include "messages.h"
#include "gas_abs_lookup.h"
#include "agenda_class.h"
#include "check_input.h"
#include "matpackV.h"
#include "physics_funcs.h"
#include "math_funcs.h"
#include "make_vector.h"
#include "arts_omp.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lookupInit(GasAbsLookup& /* x */)
{
  // Nothing to do here.
  // That means, we rely on the default constructor.

  out2 << "  Created an empty gas absorption lookup table.\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lookupCreate(// WS Output:
                      GasAbsLookup& abs_lookup,
                      Index& abs_lookup_is_adapted,
                      // WS Input:
                      const ArrayOfArrayOfSpeciesTag& abs_species,
                      const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                      const ArrayOfLineshapeSpec&     abs_lineshape,
                      const ArrayOfArrayOfSpeciesTag& abs_nls,
                      const Vector&                   f_grid,
                      const Vector&                   abs_p,
                      const Matrix&                   abs_vmrs,
                      const Vector&                   abs_t,
                      const Vector&                   abs_t_pert,
                      const Vector&                   abs_nls_pert, 
                      const Vector&                   abs_n2,            
                      const ArrayOfString&            abs_cont_names,    
                      const ArrayOfString&            abs_cont_models,   
                      const ArrayOfVector&            abs_cont_parameters )
{
  // We will be calling an absorption agenda one species at a
  // time. This is better than doing all simultaneously, because is
  // saves memory and allows for consistent treatment of nonlinear
  // species. But it means we need local copies of species, line list,
  // and line shapes for agenda communication.
  
  // 1. Output of absorption calculations:

  // Absorption coefficients:
  Matrix these_abs_coef;

  // Absorption cross sections per tag group. 
  ArrayOfMatrix abs_xsec_per_species;


  // 2. Determine various important sizes:
  const Index n_species = abs_species.nelem();   // Number of abs species
  const Index n_nls = abs_nls.nelem();           // Number of nonlinear species
  const Index n_f_grid = f_grid.nelem();         // Number of frequency grid points
  const Index n_p_grid = abs_p.nelem();          // Number of presure grid points
  const Index n_t_pert = abs_t_pert.nelem();     // Number of temp. perturbations
  const Index n_nls_pert = abs_nls_pert.nelem(); // Number of VMR pert. for NLS

  // 3. Input to absorption calculations:

  // Absorption vmrs and temperature:
  Matrix this_vmr(1,n_p_grid);
  Vector abs_h2o(n_p_grid);
  Vector this_t(n_p_grid);

  // Species list, lines, and line shapes, all with only 1 element:
  ArrayOfArrayOfSpeciesTag this_species(1);
  ArrayOfArrayOfLineRecord these_lines(1);
  ArrayOfLineshapeSpec this_lineshape(1);

  // Local copy of nls_pert and t_pert:
  Vector these_nls_pert;        // Is resized every time when used
  Vector these_t_pert;          // Is resized later on

  // 4. Checks of input parameter correctness:

  const Index h2o_index 
    = find_first_species_tg( abs_species,
                             species_index_from_species_name("H2O") );

  if ( h2o_index < 0 )
    {
      // If there are nonlinear species, then at least one species must be
      // H2O. We will use that to set h2o_abs, and to perturb in the case
      // of nonlinear species.
      if (n_nls>0)
        {
          ostringstream os;
          os << "If you have nonlinear species, at least one\n"
             << "species must be a H2O species.";
          throw runtime_error( os.str() );
        }
      else
        {
          out2 << "  You have no H2O species. Absorption continua will not work.\n"
               << "  You should get a runtime error if you try them anyway.\n";
        }
    }

  // abs_species, f_grid, and p_grid should not be empty:
  if ( 0==n_species ||
       0==n_f_grid ||
       0==n_p_grid )
    {
      ostringstream os;
      os << "One of the required input variables is empty:\n"
         << "abs_species.nelem() = " << n_species << ",\n"
         << "f_grid.nelem() = " << n_f_grid << ",\n"
         << "abs_p.nelem() = " << n_p_grid << ".";
      throw runtime_error( os.str() );
    }

  // Set up the index array abs_nls from the tag array
  // abs_nls. Give an error message if these
  // tags are not included in abs_species. 
  ArrayOfIndex abs_nls_idx;  
  for (Index i=0; i<n_nls; ++i)
    {
      Index s;
      for (s=0; s<n_species; ++s)
        {
          if (abs_nls[i]==abs_species[s])
            {
              abs_nls_idx.push_back(s);
              break;
            }
        }
      if (s==n_species)
        {
          ostringstream os;
          os << "Did not find *abs_nls* tag group \""
             << get_tag_group_name(abs_nls[i])
             << "\" in *abs_species*.";
          throw runtime_error( os.str() );
        }
    }

  // Furthermore abs_nls_idx must not contain duplicate values:
  if ( !is_unique(abs_nls_idx) )
    {
      ostringstream os;
      os << "The variable *abs_nls* must not have duplicate species.\n"
         << "Value of *abs_nls*: " << abs_nls_idx;
      throw runtime_error( os.str() );
    }

  // VMR matrix must match species list and pressure levels:
  chk_size( "abs_vmrs",
            abs_vmrs,
            n_species,
            n_p_grid );

  // Temperature vector must match number of pressure levels:
  chk_size( "abs_t",
            abs_t,
            n_p_grid ); 

  // abs_nls_pert should only be not-empty if we have nonlinear species:
  if ( ( 0==n_nls && 0 != n_nls_pert ) ||
       ( 0!=n_nls && 0 == n_nls_pert ))
    {
      ostringstream os;
      os << "You have to set both abs_nls and abs_nls_pert, or none.";
      throw runtime_error( os.str() );
    }


  // 4.a Set up a logical array for the nonlinear species.
  ArrayOfIndex non_linear(n_species,0);
  for ( Index s=0; s<n_nls; ++s )
    {
      non_linear[abs_nls_idx[s]] = 1;
    }


  // 5. Set general lookup table properties:
  abs_lookup.species = abs_species;    // Species list
  abs_lookup.nonlinear_species = abs_nls_idx;  // Nonlinear species   (e.g., H2O, O2)
  abs_lookup.f_grid = f_grid;           // Frequency grid
  abs_lookup.p_grid = abs_p;          // Pressure grid
  abs_lookup.vmrs_ref = abs_vmrs;
  abs_lookup.t_ref = abs_t;
  abs_lookup.t_pert = abs_t_pert;
  abs_lookup.nls_pert = abs_nls_pert;

  // 6. Create abs_lookup.xsec with the right dimensions:
  {
    Index a,b,c,d;

    if ( 0 == n_t_pert ) a = 1;
    else a = n_t_pert;

    b = n_species + n_nls * ( n_nls_pert - 1 );

    c = n_f_grid;

    d = n_p_grid;

    abs_lookup.xsec.resize( a, b, c, d );
  }


  // 6.a. Set up these_t_pert. This is done so that we can use the
  // same loop over the perturbations, independent of
  // whether we have temperature perturbations or not.
  if ( 0!=n_t_pert)
    {
      out2 << "  With temperature perturbations.\n";
      these_t_pert.resize(n_t_pert);
      these_t_pert = abs_t_pert;
    }
  else
    {
      out2 << "  No temperature perturbations.\n";
      these_t_pert.resize(1);
      these_t_pert = 0;
    }

  const Index these_t_pert_nelem = these_t_pert.nelem();
  

  // 7. Now we have to fill abs_lookup.xsec with the right values!

  // Loop species:
  for ( Index i=0,spec=0; i<n_species; ++i )
    {
      // spec is the index for the second dimension of abs_lookup.xsec.
      
      // Prepare absorption agenda input for this species:
      out2 << "  Doing species " << i+1 << " of " << n_species << ": "
           << abs_species[i] << ".\n";

      // Get a dummy list of tag groups with only the current element:
      this_species[0].resize(abs_species[i].nelem());
      this_species[0] = abs_species[i];

      // List of lines:
      these_lines[0].resize(abs_lines_per_species[i].nelem());
      these_lines[0] = abs_lines_per_species[i];
      
      // List of lineshapes. This requires special treatment: If there
      // is only 1 lineshape given, the same lineshape should be used
      // for all species.
      if (1==abs_lineshape.nelem())
        this_lineshape[0] = abs_lineshape[0];
      else
        this_lineshape[0] = abs_lineshape[i];

      // Set up these_nls_pert. This is done so that we can use the
      // same loop over the perturbations, independent of
      // whether we have nonlinear species or not.
      if ( non_linear[i] )
        {
          out2 << "  This is a species with H2O VMR perturbations.\n";
          these_nls_pert.resize(n_nls_pert);
          these_nls_pert = abs_nls_pert;
        }
      else
        {
          these_nls_pert.resize(1);
          these_nls_pert = 1;
        }
      
      // Loop these_nls_pert:
      for ( Index s=0; s<these_nls_pert.nelem(); ++s,++spec )
        {
          // Remember, spec is the index for the second dimension of
          // abs_lookup.xsec
          
          if ( non_linear[i] )
            {
              out2 << "  Doing H2O VMR variant " << s+1 << " of " << n_nls_pert << ": "
                   << abs_nls_pert[s] << ".\n";
            }

          // VMR for this species:
          this_vmr(0,joker) = abs_vmrs(i,joker);  
          if ( i==h2o_index )
            {
              //              out3 << "  Species is main H2O species.\n";
              this_vmr(0,joker) *= these_nls_pert[s]; // Add perturbation
            }

          // For abs_h2o, we can always add the perturbations (it will
          // not make a difference if the species itself is also H2O).
          // Attention, we need to treat here also the case that there
          // is no H2O species. We will then set abs_h2o to
          // -1. Absorption routines that do not really need abs_h2o
          // will still run.
          if ( h2o_index == -1 )
            {
              // The case without H2O species.
              abs_h2o.resize(1);
              abs_h2o = -1;
            }
          else
            {
              // The normal case.
              abs_h2o = abs_vmrs(h2o_index, joker);   
              abs_h2o *= these_nls_pert[s]; // Add perturbation
            }

          // Loop temperature perturbations:
#ifdef _OPENMP
#pragma omp parallel private(this_t, abs_xsec_per_species)
#pragma omp for 
#endif
          for ( Index j=0; j<these_t_pert_nelem; ++j )
            {
              if ( 0!=n_t_pert )
                {
                  out3 << "  Doing temperature variant " << j+1
                       << " of " << n_t_pert << ": "
                       << these_t_pert[j] << ".\n";
                }
              
              // Create perturbed temperature profile:
              this_t = abs_lookup.t_ref;
              this_t += these_t_pert[j];
      
              // The sequence of function calls here is inspired from
              // abs_coefCalcSaveMemory. 

              abs_xsec_per_speciesInit( abs_xsec_per_species, this_species,
                                        f_grid, abs_p );

              abs_xsec_per_speciesAddLines( abs_xsec_per_species,
                                            this_species,
                                            f_grid,
                                            abs_p,
                                            this_t,
                                            abs_h2o,
                                            this_vmr,
                                            these_lines,
                                            this_lineshape );

              abs_xsec_per_speciesAddConts( abs_xsec_per_species,
                                            this_species,
                                            f_grid,
                                            abs_p,
                                            this_t,
                                            abs_n2,
                                            abs_h2o,
                                            this_vmr,
                                            abs_cont_names,
                                            abs_cont_parameters,
                                            abs_cont_models);

              // Store in the right place:
              // Loop through all altitudes
              for ( Index p=0; p<n_p_grid; ++p )
                {
                  abs_lookup.xsec( j, spec, Range(joker), p )
                    = abs_xsec_per_species[0](Range(joker),p);

                  // There used to be a division by the number density
                  // n here. This is no longer necessary, since
                  // abs_xsec_per_species now contains true absorption
                  // cross sections.

                  // IMPORTANT: There was a bug in my old Matlab
                  // function "create_lookup.m" to generate the lookup
                  // table. (The number density was always the
                  // reference one, and did not change for different
                  // temperatures.) Patricks Atmlab function
                  // "arts_abstable_from_arts1.m" did *not* have this bug.

                  // Calculate the number density for the given pressure and
                  // temperature: 
                  // n = n0*T0/p0 * p/T or n = p/kB/t, ideal gas law
                  //                  const Numeric n = number_density( abs_lookup.p_grid[p],
                  //                                                    this_t[p]   );
                  //                  abs_lookup.xsec( j, spec, Range(joker), p ) /= n;
                }
            }
        }
    }

  // Set the abs_lookup_is_adapted flag. After all, the table fits the
  // current frequency grid and species selection.
  abs_lookup_is_adapted = 1;
}


//! Find continuum species in abs_species.
/*! 
  Returns an index array with indexes of those species in abs_species
  that have continuum tags that require h2o_abs, and hence require
  nonlinear treatment in the absorption lookup table.

  H2O itself is ignored here since that is treated separately.

  We are a bit conservative here, it is possible that some of the
  continua do not really require H2O. Check yourself, if you want, and
  improve the guessing here.

  \retval   cont         indices of those species with continua.
  \param    abs_species  list of absorption species.
  
  \author Stefan Buehler
  \date   2007-11-16
*/
void find_nonlinear_continua(ArrayOfIndex& cont,
                             const ArrayOfArrayOfSpeciesTag& abs_species)
{
  cont.resize(0);
  
  // This is quite complicated, unfortunately. The approach here
  // is copied from abs_xsec_per_speciesAddConts. For explanation,
  // see there.

  // Loop tag groups:
  for ( Index i=0; i<abs_species.nelem(); ++i )
    {
      extern const Array<SpeciesRecord> species_data; 

      // Loop tags in tag group
      for ( Index s=0; s<abs_species[i].nelem(); ++s )
        {

          // Check if this is a specific isotope, or "all". ("all" is
          // not a continuum tag.)
          if ( abs_species[i][s].Isotope() <
               species_data[abs_species[i][s].Species()].Isotope().nelem() )
            {
              // Continuum tags have isotope ratio -1
              if ( 0 >
                   species_data[abs_species[i][s].Species()].Isotope()[abs_species[i][s].Isotope()].Abundance() )
                {
                  const String thisname = abs_species[i][s].Name();
                  // Ok, now we know this is a continuum tag.
                  out3 << "  Continuum tag: " << thisname;

                  // Check whether we want nonlinear treatment for
                  // this or not. We have three classes of continua:
                  // 1. Those that we know do not require it
                  // 2. Those that we know require h2o_abs
                  // 3. Those for which we do not know

                  // The list here is not at all perfect, just a first
                  // guess. If you need particular models, you better
                  // check that the right thing is done with your model.

                  // 1. Continua known to not use h2o_abs
                  // We take also H2O itself here, since this is
                  // handled separately
                  if ( species_index_from_species_name("H2O")==abs_species[i][s].Species() ||
                       "N2-"==thisname.substr(0,3) ||
                       "CO2-"==thisname.substr(0,4) ||
                       "O2-CIA"==thisname.substr(0,6) ||
                       "O2-v0v"==thisname.substr(0,6) ||
                       "O2-v1v"==thisname.substr(0,6) )
                    {
                      out3 << " --> not added.\n";
                      break;
                    }

                  // 2. Continua known to use h2o_abs
                  if ( "O2-"==thisname.substr(0,3) )                       
                    {
                      cont.push_back(i);
                      out3 << " --> added to abs_nls.\n";
                      break;                      
                    }

                  // If we get here, then the tag was neither in the
                  // posivitive nor in the negative list. We through a
                  // runtime error.
                  out3 << " --> unknown.\n";
                  throw runtime_error("I don't know whether this tag uses h2o_abs or not.\n"
                                      "Cannot set abs_nls automatically.");
                }
            }              
        }
    }
}


//! Choose species for abs_nls
/*!
  Make an intelligent choice for abs_nls, based on abs_species.

  \author Stefan Buehler

  \param[out] abs_nls     The list of nonlinear species.
  \param[in]  abs_species Absorption species.
*/
void choose_abs_nls(ArrayOfArrayOfSpeciesTag& abs_nls,
                    const ArrayOfArrayOfSpeciesTag& abs_species)
{
  abs_nls.resize(0);

  // Add all H2O species as non-linear:
  Index next_h2o = 0;
  while (-1 != (next_h2o =
                find_next_species_tg(abs_species,
                                     species_index_from_species_name("H2O"),
                                     next_h2o)))
    {
      abs_nls.push_back(abs_species[next_h2o]);
      ++next_h2o;
    }

  // Certain continuum models also depend on abs_h2o. There is a
  // helper function that contains a list of these.
  ArrayOfIndex cont;
  find_nonlinear_continua(cont, abs_species);

  // Add these to abs_nls:
  for (Index i=0; i<cont.nelem(); ++i)
    {
      abs_nls.push_back(abs_species[cont[i]]);
    }

  out2 << "  Species marked for nonlinear treatment:\n";
  for (Index i=0; i<abs_nls.nelem(); ++i)
    {
      out2 << "  ";
      for (Index j=0; j<abs_nls[i].nelem(); ++j)
        {
          if (j!=0) out2 << ", ";
          out2 << abs_nls[i][j].Name();
        }
      out2 << "\n";
    }
}


//! Chose the temperature perturbations abs_t_pert
/*!  
  This simple function creates a vector of temperature
  perturbations, relative to the reference temperature profile, that
  covers the minimum and maximum temperature profile. The value 0 is
  always included.
  
  \author Stefan Buehler

  \param[out] abs_t_pert Temperature perturbations
  \param[in] abs_t       Reference temperature profile
  \param[in] tmin        Minimum temperature profile
  \param[in] tmax        Maximum temperature profile
  \param[in] t_step      Temperature perturbation step
*/
void choose_abs_t_pert(Vector&         abs_t_pert,
                       ConstVectorView abs_t,
                       ConstVectorView tmin,
                       ConstVectorView tmax,
                       const Numeric&  t_step)
{
  Vector tmindist = tmin;
  Vector tmaxdist = tmax;

  tmindist -= abs_t;
  tmaxdist -= abs_t;

  const Numeric mindev = min(tmindist);
  const Numeric maxdev = max(tmaxdist);

  out3 << "  mindev/maxdev : " << mindev << " / " << maxdev << "\n";
  Numeric start = 0;
  // n is used to make sure that there are at least 5 points in total,
  // as required for 4th order interpolation.
  // In other words, we add at least two points above and below
  // the reference profile.
  Index   n_points = 0;                
  while ( (mindev<start) || (n_points<2) ) 
    {
      start -= t_step;
      ++n_points;
    }
  
  Numeric stop  = 0;
  n_points = 0;                
  while ( (maxdev > stop) || (n_points<2) )
    {
      stop += t_step;
      ++n_points;
    }

  linspace(abs_t_pert,start,stop,t_step);

  // Special treatment for the case that there actually are no
  // temperature variations present (1-D case blown up to 3D):
  if (1==abs_t_pert.nelem())
    {
      out2 << "  No temperature variations, choosing default.\n";
      abs_t_pert = MakeVector(-t_step, 0, t_step);
    }

  out2 << "  abs_t_pert: " << abs_t_pert[0] << " K to " << abs_t_pert[abs_t_pert.nelem()-1]
       << " K in steps of " << t_step
       << " K (" << abs_t_pert.nelem() << " grid points)\n";

}


//! Chose the H2O perturbations abs_nls_pert
/*!  
  This simple function creates a vector of fractional H2O VMR
  perturbations, relative to the reference H2O profile, that
  covers the minimum and maximum profile. The value 1 is
  always included.
  
  \author Stefan Buehler

  \param[out] abs_nls_pert H2O VMR perturbations
  \param[in] refprof       Reference profile
  \param[in] minprof       Minimum profile
  \param[in] maxprof       Maximum profile
  \param[in] step          Fractional perturbation step
*/
void choose_abs_nls_pert(Vector&         abs_nls_pert,
                       ConstVectorView refprof,
                       ConstVectorView minprof,
                       ConstVectorView maxprof,
                       const Numeric&  step)
{
  Vector minprofdist = minprof;
  Vector maxprofdist = maxprof;

  minprofdist /= refprof;
  maxprofdist /= refprof;

  const Numeric mindev = min(minprofdist);
  const Numeric maxdev = max(maxprofdist);

  if (mindev<0) 
    out2 << "  Warning: I am getting a negative fractional distance to the H2O\n"
         << "  reference profile. Some of your H2O fields contain negative values.\n";

  out3 << "  mindev/maxdev : " << mindev << " / " << maxdev << "\n";
  Numeric start = 1;
  // n is used to make sure that there are at least 5 points in total,
  // as required for 4th order interpolation.
  // In other words, we add at least two points above and below
  // the reference profile.
  Index   n_points = 0;                
  while ( (mindev < start) || (n_points<2) )
    {
      start -= step;
      ++n_points;
    }

  Numeric stop  = 1;
  n_points = 0;                
  while ( (maxdev > stop) || (n_points<2) )
    {
      stop += step;
      ++n_points;
    }

  linspace(abs_nls_pert,start,stop,step);

  // Special treatment for the case that there actually are no
  // H2O variations present (1-D case blown up to 3D):
  if (1==abs_nls_pert.nelem())
    {
      out2 << "  No H2O variations, choosing default.\n";
      abs_nls_pert = MakeVector(1-step, 1, 1+step);
    }

  out2 << "  abs_nls_pert: " << abs_nls_pert[0] << " to " << abs_nls_pert[abs_nls_pert.nelem()-1]
       << " (fractional units) in steps of " << step
       << " (" << abs_nls_pert.nelem() << " grid points)\n";

}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lookupSetup(// WS Output:
                     Vector& abs_p,
                     Vector& abs_t,
                     Vector& abs_t_pert,
                     Matrix& abs_vmrs,
                     ArrayOfArrayOfSpeciesTag& abs_nls,
                     Vector& abs_nls_pert,
                     // WS Input:
                     const Index& atmosphere_dim,
                     const Vector& p_grid,
                     const Vector& lat_grid,
                     const Vector& lon_grid,
                     const Tensor3& t_field,
                     const Tensor4& vmr_field,
                     const ArrayOfArrayOfSpeciesTag& abs_species,
                     // Control Parameters:
                     const Numeric& p_step,
                     const Numeric& t_step,
                     const Numeric& h2o_step)
{
  // Checks on input parameters:
  
  // We don't actually need lat_grid and lon_grid, but we have them as
  // input variables, so that we can use the standard functions to
  // check atmospheric fields and grids. A bit cheesy, but I don't
  // want to program all the checks explicitly.

  // Check grids:
  chk_atm_grids(atmosphere_dim, p_grid, lat_grid, lon_grid);

  // Check T field:
  chk_atm_field("t_field", t_field, atmosphere_dim,
                p_grid, lat_grid, lon_grid);
 
  // Check VMR field (and abs_species):
  chk_atm_field("vmr_field", vmr_field, atmosphere_dim,
                abs_species.nelem(), p_grid, lat_grid, lon_grid);

  // Check the keyword arguments:
  if ( p_step <= 0 || t_step <= 0 || h2o_step <= 0 )
    {
      ostringstream os;
      os << "The keyword arguments p_step, t_step, and h2o_step must be >0.";
      throw runtime_error( os.str() );
    }

  // Ok, all input parameters seem to be reasonable.


  // We will need the log of the pressure grid:
  Vector log_p_grid(p_grid.nelem());
  transform(log_p_grid, log, p_grid);

  //  const Numeric epsilon = 0.01 * p_step; // This is the epsilon that
  //                                         // we use for comparing p grid spacings.

  // Construct abs_p
  // ---------------

  ArrayOfNumeric log_abs_p_a;  // We take log_abs_p_a as an array of
                             // Numeric, so that we can easily 
                             // build it up by appending new elements to the end. 

  // Check whether there are pressure levels that are further apart
  // (in log(p)) than p_step, and insert additional levels if
  // necessary:

  log_abs_p_a.push_back(log_p_grid[0]);

  for (Index i=1; i<log_p_grid.nelem(); ++i)
    {
      const Numeric dp =  log_p_grid[i-1] - log_p_grid[i]; // The grid is descending.

      const Numeric dp_by_p_step = dp/p_step;
      //          cout << "dp_by_p_step: " << dp_by_p_step << "\n";

      // How many times does p_step fit into dp?
      const Index n = (Index) ceil(dp_by_p_step); 
      // n is the number of intervals that we want to have in the
      // new grid. The number of additional points to insert is
      // n-1. But we have to insert the original point as well.
      //          cout << n << "\n";

      const Numeric ddp = dp/n;
      //          cout << "ddp: " << ddp << "\n";

      for (Index j=1; j<=n; ++j)
        log_abs_p_a.push_back(log_p_grid[i-1] - j*ddp);          
    }

  // Copy to a proper vector, we need this also later for
  // interpolation: 
  Vector log_abs_p(log_abs_p_a.nelem());
  for (Index i=0; i<log_abs_p_a.nelem(); ++i)
    log_abs_p[i] = log_abs_p_a[i];

  // Copy the new grid to abs_p, removing the log:
  abs_p.resize(log_abs_p.nelem());
  transform(abs_p, exp, log_abs_p);


  // We will also have to interpolate T and VMR profiles to the new
  // pressure grid. We interpolate in log(p), as usual in ARTS.

  // Grid positions:
  ArrayOfGridPos gp(log_abs_p.nelem());
  gridpos(gp, log_p_grid, log_abs_p);

  // Interpolation weights:
  Matrix itw(gp.nelem(),2);
  interpweights(itw,gp);


  // In the 1D case the lookup table is just a lookup table in
  // pressure. We treat this simple case first.
  if (1==atmosphere_dim)
    {
      // Reference temperature,
      // interpolate abs_t from t_field:
      abs_t.resize(log_abs_p.nelem());
      interp(abs_t,
             itw,
             t_field(joker,0,0),
             gp);

      // Temperature perturbations:
      abs_t_pert.resize(0);

      // Reference VMR profiles,
      // interpolate abs_vmrs from vmr_field:
      abs_vmrs.resize(abs_species.nelem(), log_abs_p.nelem());
      for (Index i=0; i<abs_species.nelem(); ++i)
        interp(abs_vmrs(i,joker),
               itw,
               vmr_field(i,joker,0,0),
               gp);

      // Species for which H2O VMR is perturbed:
      abs_nls.resize(0);

      // H2O VMR perturbations:
      abs_nls_pert.resize(0);
    }
  else
    {
      // 2D or 3D case. We have to set up T and nonlinear species variations.

      // Make an intelligent choice for the nonlinear species.
      choose_abs_nls(abs_nls, abs_species);

      // Now comes a part where we analyse the atmospheric fields.
      // We need to find the max, min, and mean profile for
      // temperature and VMRs.
      // We do this on the original p grid, not on the refined p
      // grid, to be more efficient.

      // Temperature:
      Vector tmin(p_grid.nelem());
      Vector tmax(p_grid.nelem());
      Vector tmean(p_grid.nelem()); 

      for (Index i=0; i<p_grid.nelem(); ++i)
        {
          tmin[i]  = min(t_field(i,joker,joker));
          tmax[i]  = max(t_field(i,joker,joker));
          tmean[i] = mean(t_field(i,joker,joker));
        }
      
//       cout << "Tmin: " << tmin << "\n";
//       cout << "Tmax: " << tmax << "\n";
//       cout << "Tmean: " << tmean << "\n";

      // Calculate mean profiles of all species. (We need all for abs_vmrs
      // not only H2O.)
      Matrix vmrmean(abs_species.nelem(), p_grid.nelem()); 
      for (Index s=0; s<abs_species.nelem(); ++s)
        for (Index i=0; i<p_grid.nelem(); ++i)
          {
            vmrmean(s,i) = mean(vmr_field(s,i,joker,joker));
          }      

      // If there are NLS, determine H2O statistics:

      // We have to define these here, outside the if block, because
      // we need the values later.
      Vector h2omin(p_grid.nelem());
      Vector h2omax(p_grid.nelem());
      const Index h2o_index 
        = find_first_species_tg( abs_species,
                                 species_index_from_species_name("H2O") );
      // We need this inside the if clauses for nonlinear species
      // treatment. The function returns "-1" if there is no H2O
      // species. There is a check for that in the next if block, with
      // an appropriate runtime error.

      if (0<abs_nls.nelem())
        {
          // Check if there really is a H2O species.
          if (h2o_index<0)
            {
              ostringstream os;
              os << "Some of your species require nonlinear treatment,\n"
                 << "but you have no H2O species.";
              throw runtime_error( os.str() );
            }

          for (Index i=0; i<p_grid.nelem(); ++i)
            {
              h2omin[i]  = min(vmr_field(h2o_index,i,joker,joker));
              h2omax[i]  = max(vmr_field(h2o_index,i,joker,joker));
            }
      
//           cout << "H2Omin: " << h2omin << "\n";
//           cout << "H2Omax: " << h2omax << "\n";
//           cout << "H2Omean: " << vmrmean(h2o_index,joker) << "\n";

        }


      // Interpolate in pressure, set abs_t, abs_vmr...

      // Reference temperature,
      // interpolate abs_t from tmean:
      abs_t.resize(log_abs_p.nelem());
      interp(abs_t,
             itw,
             tmean,
             gp);

      // Temperature perturbations:
      choose_abs_t_pert(abs_t_pert, tmean, tmin, tmax, t_step);
//       cout << "abs_t_pert: " << abs_t_pert << "\n";

      // Reference VMR profiles,
      // interpolate abs_vmrs from vmrmean:
      abs_vmrs.resize(abs_species.nelem(), log_abs_p.nelem());
      for (Index i=0; i<abs_species.nelem(); ++i)
        interp(abs_vmrs(i,joker),
               itw,
               vmrmean(i,joker),
               gp);

      if (0<abs_nls.nelem())
        {
          // Construct abs_nls_pert:
          choose_abs_nls_pert(abs_nls_pert,
                              vmrmean(h2o_index, joker),
                              h2omin,
                              h2omax,
                              h2o_step);
        }
      else
        {
          // Empty abs_nls_pert:
          abs_nls_pert.resize(0);
        }
//       cout << "abs_nls_pert: " << abs_nls_pert << "\n";

    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lookupSetupBatch(// WS Output:
                          Vector& abs_p,
                          Vector& abs_t,
                          Vector& abs_t_pert,
                          Matrix& abs_vmrs,
                          ArrayOfArrayOfSpeciesTag& abs_nls,
                          Vector& abs_nls_pert,
                          // WS Input:
                          const ArrayOfArrayOfSpeciesTag& abs_species,
                          const ArrayOfGriddedField4& batch_fields,
                          // Control Parameters:
                          const Numeric& p_step,
                          const Numeric& t_step,
                          const Numeric& h2o_step,
                          const Vector&  extremes)
{
  // FIXME: Some runtime input variable checks here, e.g.:
  // Field names for first batch case, is T, z, in the right place? Do
  // the species match?

  // Make an intelligent choice for the nonlinear species.
  choose_abs_nls(abs_nls, abs_species);

  // Find out maximum and minimum pressure.
  Numeric maxp=batch_fields[0].p_grid[0];
  Numeric minp=batch_fields[0].p_grid[batch_fields[0].p_grid.nelem()-1];
  for (Index i=0; i<batch_fields.nelem(); ++i)
    {
      if (maxp < batch_fields[i].p_grid[0])
        maxp = batch_fields[i].p_grid[0];
      if (minp > batch_fields[i].p_grid[batch_fields[i].p_grid.nelem()-1])
        minp = batch_fields[i].p_grid[batch_fields[i].p_grid.nelem()-1];
    }
  out3 << "  maxp/minp: " << maxp << " / " << minp << "\n";

  // We construct the pressure grid as follows:
  // - Everything is done in log(p).
  // - Start with maxp and go down in steps of p_step until we are <= minp.
  // - Adjust the final pressure value to be exactly minp, otherwise
  //   we have problems in getting min, max, and mean values for this
  //   point later.
  const Index np = (Index) ceil((log(maxp)-log(minp))/p_step)+1;
  // The +1 above has to be there because we must include also both grid end points. 
  Vector log_abs_p(log(maxp), np, -p_step);
  log_abs_p[np-1] = log(minp);

  abs_p.resize(np);
  transform(abs_p, exp, log_abs_p);  
  out2 << "  abs_p: " << abs_p[0] << " Pa to " << abs_p[np-1]
       << " Pa in log(p[hPa]) in steps of " << p_step
       << " (" << np << " grid points)\n";

  // Now we have to determine the statistics of T and VMRs, we need
  // profiles of min, max, and mean of these, on the abs_p grid.

  Index n_variables = batch_fields[0].data.nbooks();

  // The first dimension of datamin, datamax, and datamean is the
  // variable (T,Z,H2O,O3,...). The second dimension is pressure. We
  // assume all elements of the batch have the same variables.

  Matrix datamin  (n_variables, np, numeric_limits<Numeric>::max());
  Matrix datamax  (n_variables, np, numeric_limits<Numeric>::min());
  // The limits here are from the header file <limits>
  Matrix datamean (n_variables, np, 0);
  Vector mean_norm(np,0);          // Divide by this to get mean.

  for (Index i=0; i<batch_fields.nelem(); ++i)
    {
      // Check that really each case has the same variables (see
      // comment above.)
      if (batch_fields[i].field_names != batch_fields[0].field_names)
        throw runtime_error("All fields in the batch must have the same field names.");

      // Check that the first field is indeed T and the third field is
      // indeed H2O:
      if ("T" != batch_fields[i].field_names[0])
        {
          ostringstream os;
          os << "The first variable in the compact atmospheric state must be \"T\",\n"
             << "but yours is: " << batch_fields[i].field_names[0];
          throw runtime_error( os.str() );
        }

      if ("H2O" != batch_fields[i].field_names[2])
        {
          ostringstream os;
          os << "The third variable in the compact atmospheric state must be \"H2O\",\n"
             << "but yours is: " << batch_fields[i].field_names[0];
          throw runtime_error( os.str() );
        }

      // Create convenient handles:
      ConstVectorView  p_grid = batch_fields[i].p_grid;
      ConstTensor4View data   = batch_fields[i].data;
 

      // Interpolate the current batch fields to the abs_p grid. We
      // have to do this for each batch case, since the grids could
      // all be different.
      
      Vector log_p_grid(p_grid.nelem());
      transform(log_p_grid, log, p_grid);

      // There is a catch here: We can only use the part of abs_p that
      // is inside the current pressure grid p_grid, otherwise we
      // would have to extrapolate.
      const Numeric eps = (log_p_grid[0]-log_p_grid[1])/2.1;
      Index first_p=0;
      while (log_abs_p[first_p] > log_p_grid[0]+eps) ++first_p;
      Index last_p=log_abs_p.nelem()-1;
      while (log_abs_p[last_p] < log_p_grid[log_p_grid.nelem()-1]-eps) --last_p;
      const Index extent_p = last_p - first_p + 1;
      
      // This was too complicated to get right:
      //      const Index first_p   = (Index) round ( (log_abs_p[0]       - log_p_grid[0])                    / p_step);
      //      const Index extent_p  = (Index) round ( (log_abs_p[first_p] - log_p_grid[log_p_grid.nelem()-1]) / p_step) + 1;

      
      ConstVectorView active_log_abs_p = log_abs_p[Range(first_p, extent_p)];

//       cout << "first_p / last_p / extent_p : " << first_p << " / " << last_p << " / " << extent_p << "\n";
//       cout << "log_p_grid: "        << log_p_grid << "\n";
//       cout << "log_abs_p:  "        << log_abs_p << "\n";
//       cout << "active_log_abs_p:  " << active_log_abs_p << "\n";
//       cout << "=============================================================\n";
//       arts_exit();

      // Grid positions:
      ArrayOfGridPos gp(active_log_abs_p.nelem());
      gridpos(gp, log_p_grid, active_log_abs_p);
      //      gridpos(gp, log_p_grid, active_log_abs_p, 100);
      // We allow much more extrapolation here than normal (0.5 is
      // normal). If we do not do this, then we get problems for
      // p_grids that are much finer than abs_p.

      // Interpolation weights:
      Matrix itw(gp.nelem(),2);
      interpweights(itw,gp);

      // We have to loop over fields, latitudes, and longitudes here, doing the
      // pressure interpolation for all. The dimensions of data are: 
      // data[N_fields][N_p][N_lat][N_lon]
      Tensor4 data_interp(data.nbooks(),
                          active_log_abs_p.nelem(),
                          data.nrows(),
                          data.ncols());

      for (Index lo=0; lo<data.ncols(); ++lo)
        for (Index la=0; la<data.nrows(); ++la)
          for (Index fi=0; fi<data.nbooks(); ++fi)
              interp(data_interp(fi, joker, la, lo),
                     itw,
                     data(fi, joker, la, lo),
                     gp);

      // Now update our datamin, datamax, and datamean variables.
      // We need the min and max only for the T and H2O profile, 
      // not for others. But we need the mean for all. We are just
      // hopping over the case that we do not need below. This is not
      // very clean, but efficient. And it avoids handling all the
      // different cases explicitly.
      for (Index lo=0; lo<data_interp.ncols(); ++lo)
        for (Index la=0; la<data_interp.nrows(); ++la)
          {
            for (Index fi=0; fi<data_interp.nbooks(); ++fi)
              {
                if (1!=fi)      // We skip the z field, which we do not need
                  for (Index pr=0; pr<data_interp.npages(); ++pr)
                    {
                      if (0==fi || 2==fi) // Min and max only needed
                                          // for T and H2o
                        {
                          if (data_interp(fi, pr, la, lo) < datamin(fi, first_p+pr))
                            datamin(fi, first_p+pr) = data_interp(fi, pr, la, lo);
                          if (data_interp(fi, pr, la, lo) > datamax(fi, first_p+pr))
                            datamax(fi, first_p+pr) = data_interp(fi, pr, la, lo);
                        }
                      
                      datamean(fi, first_p+pr) += data_interp(fi, pr, la, lo);
                    }
              }

            // The mean_norm is actually a bit tricky. It depends on
            // pressure, since different numbers of cases contribute
            // to the mean for different pressures. At the very eges
            // of the grid, typically only a single case contributes.

            mean_norm[Range(first_p, extent_p)] += 1;
          }

    }  

  // Divide mean by mean_norm to get the mean:
  assert(np==mean_norm.nelem());
  for (Index fi=0; fi<datamean.nrows(); ++fi)
    if (1!=fi)      // We skip the z field, which we do not need
      for (Index pi=0; pi<np; ++pi)
        {
          // We do this in an explicit loop here, since we have to
          // check whether there really were data points to calculate
          // the mean at each level.
          if (0<mean_norm[pi])
            datamean(fi,pi) /= mean_norm[pi];
          else
            {
              ostringstream os;
              os << "No data at pressure level " << pi+1
                 << " of "<< np << " (" << abs_p[pi] << " hPa).";
              throw runtime_error( os.str() );
            }
        }
  // Set abs_t:
  abs_t.resize(np);
  abs_t = datamean(0,joker);
  //   cout << "abs_t: " << abs_t << "\n\n";
  //   cout << "tmin:  " << datamin(0,joker) << "\n\n";
  //   cout << "tmax:  " << datamax(0,joker) << "\n";
  
  // Set abs_vmrs:
  assert(abs_species.nelem()==datamean.nrows()-2);
  abs_vmrs.resize(abs_species.nelem(), np);
  abs_vmrs = datamean(Range(2,abs_species.nelem()),joker);
  //  cout << "\n\nabs_vmrs: " << abs_vmrs << "\n\n";
  
  // Construct abs_t_pert:
  ConstVectorView tmin = datamin(0,joker);
  ConstVectorView tmax = datamax(0,joker);
  choose_abs_t_pert(abs_t_pert, abs_t, tmin, tmax, t_step);
  //  cout << "abs_t_pert: " << abs_t_pert << "\n";

  // Construct abs_nls_pert:
  ConstVectorView h2omin = datamin(2,joker);
  ConstVectorView h2omax = datamax(2,joker);
  choose_abs_nls_pert(abs_nls_pert, abs_vmrs(0,joker), h2omin, h2omax, h2o_step);
  //  cout << "abs_nls_pert: " << abs_nls_pert << "\n";

  // Append the explicitly given user extreme values, if necessary:
  if (0!=extremes.nelem())
    {
      // There must be 4 values in this case: t_min, t_max, h2o_min, h2o_max
      if (4!=extremes.nelem())
        {
          ostringstream os;
          os << "There must be exactly 4 elements in extremes:\n"
             << "min(abs_t_pert), max(abs_t_pert), min(abs_nls_pert), max(abs_nls_pert)";
          throw runtime_error( os.str() );
        }
      
      // t_min:
      if (extremes[0] < abs_t_pert[0])
        {
          Vector dummy = abs_t_pert;
          abs_t_pert.resize(abs_t_pert.nelem()+1);
          abs_t_pert[0] = extremes[0];
          abs_t_pert[Range(1,dummy.nelem())] = dummy;
          out2 << "  Added min extreme value for abs_t_pert: "
               << abs_t_pert[0] << "\n";
        }

      // t_max:
      if (extremes[1] > abs_t_pert[abs_t_pert.nelem()-1])
        {
          Vector dummy = abs_t_pert;
          abs_t_pert.resize(abs_t_pert.nelem()+1);
          abs_t_pert[Range(0,dummy.nelem())] = dummy;
          abs_t_pert[abs_t_pert.nelem()-1] = extremes[1];
          out2 << "  Added max extreme value for abs_t_pert: "
               << abs_t_pert[abs_t_pert.nelem()-1] << "\n";
        }

      // nls_min:
      if (extremes[2] < abs_nls_pert[0])
        {
          Vector dummy = abs_nls_pert;
          abs_nls_pert.resize(abs_nls_pert.nelem()+1);
          abs_nls_pert[0] = extremes[2];
          abs_nls_pert[Range(1,dummy.nelem())] = dummy;
          out2 << "  Added min extreme value for abs_nls_pert: "
               << abs_nls_pert[0] << "\n";
        }

      // nls_max:
      if (extremes[3] > abs_nls_pert[abs_nls_pert.nelem()-1])
        {
          Vector dummy = abs_nls_pert;
          abs_nls_pert.resize(abs_nls_pert.nelem()+1);
          abs_nls_pert[Range(0,dummy.nelem())] = dummy;
          abs_nls_pert[abs_nls_pert.nelem()-1] = extremes[3];
          out2 << "  Added max extreme value for abs_nls_pert: "
               << abs_nls_pert[abs_nls_pert.nelem()-1] << "\n";
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesAdd(// WS Output:
                    ArrayOfArrayOfSpeciesTag& abs_species,
                    // Control Parameters:
                    const ArrayOfString& names)
{
  // Size of initial array
  Index n_gs = abs_species.nelem();
  
  // Temporary ArrayOfSpeciesTag
  ArrayOfSpeciesTag temp;
    
  // Each element of the array of Strings names defines one tag
  // group. Let's work through them one by one.
  for ( Index i=0; i<names.nelem(); ++i )
    {
      array_species_tag_from_string( temp, names[i] );  
      abs_species.push_back(temp);
    }

  // Print list of tag groups to the most verbose output stream:
  out3 << "  Added tag groups:";
  for ( Index i=n_gs; i<abs_species.nelem(); ++i )
    {
      out3 << "\n  " << i << ":";
      for ( Index s=0; s<abs_species[i].nelem(); ++s )
        {
          out3 << " " << abs_species[i][s].Name();
        }
    }
  out3 << '\n';
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesAdd2(// WS Output:
                    ArrayOfArrayOfSpeciesTag& abs_species,
                    ArrayOfRetrievalQuantity& jq,
                    Agenda&                   jacobian_agenda,
                    // WS Input:
                    const Matrix&             jac,
                    const Index&              atmosphere_dim,
                    const Vector&             p_grid,
                    const Vector&             lat_grid,
                    const Vector&             lon_grid,
                    // WS Generic Input:
                    const Vector&             rq_p_grid,
                    const Vector&             rq_lat_grid,
                    const Vector&             rq_lon_grid,
                    // Control Parameters:
                    const String&             species,
                    const String&             method,
                    const String&             mode,
                    const Numeric&            dx)
{
  // Add species to *abs_species*
  ArrayOfSpeciesTag tags;
  array_species_tag_from_string( tags, species );
  abs_species.push_back( tags );

  // Print list of added tag group to the most verbose output stream:
  out3 << "  Appended tag group:";
  out3 << "\n  " << abs_species.nelem()-1 << ":";
  for ( Index s=0; s<tags.nelem(); ++s )
  {
    out3 << " " << tags[s].Name();
  }
  out3 << '\n';

  // Do retrieval part
  jacobianAddAbsSpecies( jq, jacobian_agenda, jac, atmosphere_dim, 
                         p_grid, lat_grid, lon_grid, rq_p_grid, rq_lat_grid, 
                         rq_lon_grid, species, method, mode, dx);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesInit( ArrayOfArrayOfSpeciesTag& abs_species )
{
  abs_species.resize(0);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void SpeciesSet(// WS Generic Output:
                ArrayOfArrayOfSpeciesTag& abs_species,
                // Control Parameters:
                const ArrayOfString& names)
{
  abs_species.resize(names.nelem());

  //cout << "Names: " << names << "\n";

  // Each element of the array of Strings names defines one tag
  // group. Let's work through them one by one.
  for ( Index i=0; i<names.nelem(); ++i )
    {
      // This part has now been moved to array_species_tag_from_string.
      // Call this function.
      array_species_tag_from_string( abs_species[i], names[i] );  
    }

  // Print list of tag groups to the most verbose output stream:
  out3 << "  Defined tag groups: ";
  for ( Index i=0; i<abs_species.nelem(); ++i )
    {
      out3 << "\n  " << i << ":";
      for ( Index s=0; s<abs_species[i].nelem(); ++s )
        {
          out3 << " " << abs_species[i][s].Name();
        }
    }
  out3 << '\n';
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lookupAdapt( GasAbsLookup&                   abs_lookup,
                          Index&                          abs_lookup_is_adapted,
                          const ArrayOfArrayOfSpeciesTag& abs_species,
                          const Vector&                   f_grid)
{
  abs_lookup.Adapt( abs_species, f_grid );
  abs_lookup_is_adapted = 1;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_scalar_gasExtractFromLookup( Matrix&             abs_scalar_gas,
                                      const GasAbsLookup& abs_lookup,
                                      const Index&        abs_lookup_is_adapted, 
                                      const Index&        f_index,
                                      const Numeric&      a_pressure,
                                      const Numeric&      a_temperature,
                                      const Vector&       a_vmr_list)
{
  // Check if the table has been adapted:
  if ( 1!=abs_lookup_is_adapted )
    throw runtime_error("Gas absorption lookup table must be adapted,\n"
                        "use method abs_lookupAdapt.");

  // The function we are going to call here is one of the few helper
  // functions that adjust the size of their output argument
  // automatically. 
  abs_lookup.Extract( abs_scalar_gas,
                          f_index,
                          a_pressure,
                          a_temperature,
                          a_vmr_list );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_fieldCalc(// WS Output:
                   Tensor5& asg_field,
                   // WS Input:
                   const Agenda&  sga_agenda,
                   const Index&   f_index,
                   const Vector&  f_grid,
                   const Index&   atmosphere_dim,
                   const Vector&  p_grid,
                   const Vector&  lat_grid,
                   const Vector&  lon_grid,
                   const Tensor3& t_field,
                   const Tensor4& vmr_field )
{
  Matrix  asg;
  Vector  a_vmr_list;

  // Get the number of species from the leading dimension of vmr_field:
  const Index n_species = vmr_field.nbooks();

  // Number of frequencies:
  const Index n_frequencies = f_grid.nelem();

  // Number of pressure levels:
  const Index n_pressures = p_grid.nelem();

  // Number of latitude grid points (must be at least one):
  const Index n_latitudes = max( Index(1), lat_grid.nelem() );

  // Number of longitude grid points (must be at least one):
  const Index n_longitudes = max( Index(1), lon_grid.nelem() );
  
  // Check grids:
  chk_atm_grids( atmosphere_dim,
                 p_grid,
                 lat_grid,
                 lon_grid );
  
  // Check if t_field is ok:
  chk_atm_field( "t_field",
                 t_field,
                 atmosphere_dim,
                 p_grid,
                 lat_grid,
                 lon_grid );

  // Check if vmr_field is ok.
  // (Actually, we are not checking the first dimension, since
  // n_species has been set from this.)
  chk_atm_field( "vmr_field",
                 vmr_field,
                 atmosphere_dim,
                 n_species,
                 p_grid,
                 lat_grid,
                 lon_grid );

  // We also set the start and extent for the frequency loop.
  Index f_extent;

  if ( f_index < 0 )
    {
      // This means we should extract for all frequencies.

      f_extent = n_frequencies;
    }
  else
    {
      // This means we should extract only for one frequency.

      // Make sure that f_index is inside f_grid:
      if ( f_index >= n_frequencies )
        {
          ostringstream os;
          os << "The frequency index f_index points to a frequency outside"
             << "the frequency grid. (f_index = " << f_index
             << ", n_frequencies = " << n_frequencies << ")";
          throw runtime_error( os.str() );
        }

      f_extent = 1;
    }

  // Resize output field.
  // The dimension in lat and lon must be at least one, even if these
  // grids are empty.
  out2 << "  Creating field with dimensions:\n"
       << "    " << n_species << "    gas species,\n"
       << "    " << f_extent << "     frequencies,\n"
       << "    " << n_pressures << "  pressures,\n"
       << "    " << n_latitudes << "  latitudes,\n"
       << "    " << n_longitudes << " longitudes.\n";

  asg_field.resize( n_species,
                    f_extent,
                    n_pressures,
                    n_latitudes,
                    n_longitudes );


  out2 << "  Agenda output is suppressed, use reporting\n"
       <<"   level 4 if you want to see it.\n";

  // Now we have to loop all points in the atmosphere:
#ifdef _OPENMP
#pragma omp parallel private(asg, a_vmr_list)
#pragma omp for 
#endif
  for ( Index ipr=0; ipr<n_pressures; ++ipr )         // Pressure:  ipr
    {
      Numeric a_pressure = p_grid[ipr];

      out3 << "  p_grid[" << ipr << "] = " << a_pressure << "\n";

      for ( Index ila=0; ila<n_latitudes; ++ila )   // Latitude:  ila
        for ( Index ilo=0; ilo<n_longitudes; ++ilo ) // Longitude: ilo
          {
            Numeric a_temperature = t_field( ipr, ila, ilo );
            a_vmr_list    = vmr_field( Range(joker),
                                       ipr, ila, ilo );

            // Execute agenda to calculate local absorption.
            // Agenda input:  f_index, a_pressure, a_temperature, a_vmr_list
            // Agenda output: asg
            abs_scalar_gas_agendaExecute(asg,
                                         f_index, a_pressure,
                                         a_temperature, a_vmr_list,
                                         sga_agenda, true, true);

            // Verify, that the number of species in asg is
            // constistent with vmr_field:
            if ( n_species != asg.ncols() )
              {
                ostringstream os;
                os << "The number of gas species in vmr_field is "
                   << n_species << ",\n"
                   << "but the number of species returned by the agenda is "
                   << asg.ncols() << ".";
                throw runtime_error( os.str() );
              }

            // Verify, that the number of frequencies in asg is
            // constistent with f_extent:
            if ( f_extent != asg.nrows() )
              {
                ostringstream os;
                os << "The number of frequencies desired is "
                   << n_frequencies << ",\n"
                   << "but the number of frequencies returned by the agenda is "
                   << asg.nrows() << ".";
                throw runtime_error( os.str() );
              }

            // Store the result in output field.
            // We have to transpose asg, because the dimensions of the
            // two variables are:
            // asg_field: [ abs_species, f_grid, p_grid, lat_grid, lon_grid]
            // asg:       [ f_grid, abs_species ]
            asg_field( Range(joker),
                       Range(joker),
                       ipr, ila, ilo ) = transpose( asg );
            
          }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void f_gridFromGasAbsLookup(
             Vector&         f_grid,
       const GasAbsLookup&   abs_lookup )
{
  abs_lookup.GetFgrid( f_grid );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void p_gridFromGasAbsLookup(
             Vector&         p_grid,
       const GasAbsLookup&   abs_lookup )
{
  abs_lookup.GetPgrid( p_grid );
}

