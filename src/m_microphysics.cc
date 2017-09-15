/* Copyright (C) 2011-2017
   Jana Mendrok     <jana.mendrok@gmail.com>
   Daniel Kreyling  <daniel.kreyling@nict.go.jp>
   Manfred Brath    <manfred.brath@uni-hamburg.de>
   Patrick Eriksson <patrick.eriksson@chalmers.se>
                         
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



/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   m_microphysics.cc
  \author Jana Mendrok, Daniel Kreyling, Manfred Brath, Patrick Eriksson
  \date   2017-07-10 

  \brief  Workspace functions related to particle micophysics (e.g. size
          distributions).

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.*/



/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <cmath>
#include <stdexcept>
#include <cstdlib>

#include "arts.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "cloudbox.h"
#include "file.h"
#include "interpolation.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "microphysics.h"
#include "optproperties.h"
#include "parameters.h"
#include "physics_funcs.h"
#include "rte.h"
#include "sorting.h"
#include "special_interp.h"
#include "xml_io.h"

extern const String SCATSPECIES_MAINTAG;


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void particle_massesFromMetaDataSingleCategory(
        Matrix&                    particle_masses,
  const ArrayOfArrayOfScatteringMetaData& scat_meta,
  const Verbosity&)
{
  const Index np_total = TotalNumberOfElements(scat_meta);

  particle_masses.resize(np_total,1);

  Index i_se_flat = 0;
  for( Index i_ss=0; i_ss<scat_meta.nelem(); i_ss++ )
  {
      for( Index i_se=0; i_se<scat_meta[i_ss].nelem(); i_se++ )
      {
          if( isnan(scat_meta[i_ss][i_se].mass) ||
              scat_meta[i_ss][i_se].mass <= 0 ||
              scat_meta[i_ss][i_se].mass > 1. )
          {
              ostringstream os;
              os << "A presumably incorrect value found for "
              << "scat_meta[" << i_ss << "][" << i_se << "].mass.\n"
              << "The value is " << scat_meta[i_ss][i_se].mass;
              throw std::runtime_error(os.str());
          }

          if( scat_meta[i_ss][i_se].diameter_volume_equ <= 0 ||
             scat_meta[i_ss][i_se].diameter_volume_equ > 0.5 )
          {
              ostringstream os;
              os << "A presumably incorrect value found for "
              << "scat_meta[" << i_ss << "][" << i_se << "].diameter_volume_equ.\n"
              << "The value is " << scat_meta[i_ss][i_se].diameter_volume_equ;
              throw std::runtime_error(os.str());
          }

          particle_masses(i_se_flat,0) = scat_meta[i_ss][i_se].diameter_volume_equ;

          i_se_flat++;
      }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void particle_massesFromMetaData
                        (//WS Output:
                         Matrix& particle_masses,
                         // WS Input:
                         const ArrayOfArrayOfScatteringMetaData& scat_meta,
                         const Verbosity& )
{
  // resize particle_masses to required diemsions and properly initialize values
  particle_masses.resize ( TotalNumberOfElements(scat_meta), scat_meta.nelem() );
  particle_masses = 0.;

  // calculate and set particle_masses
  Index i_se_flat = 0;
  for ( Index i_ss=0; i_ss<scat_meta.nelem(); i_ss++ )
  {
    for ( Index i_se=0; i_se<scat_meta[i_ss].nelem(); i_se++ )
    {
      particle_masses (i_se_flat, i_ss) = scat_meta[i_ss][i_se].diameter_volume_equ;
      i_se_flat++;
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pndFromdNdD (//WS Output:
              Vector& pnd,
              //WS Input:
              const Vector& dNdD,
              const Vector& diameter,
              const Numeric& total_content,
              const Vector& scatelem_content,
              const Verbosity& verbosity)
{
  pnd.resize( diameter.nelem() );

  // scale dNdD by size bin width
  if (diameter.nelem() > 1)
      bin_integral( pnd, diameter, dNdD ); //[# m^-3]
  else
      pnd = dNdD;

  // scaling pnd to real mass/number density/flux (some PSDs have implicit
  // scaling - then this is only a check -, others don't
  chk_pndsum ( pnd, total_content, scatelem_content,
               0, 0, 0, "your field", verbosity );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_size_gridFromScatMeta(
          Vector&                             pnd_size_grid,
    const ArrayOfArrayOfScatteringMetaData&   scat_meta,
    const Index&                              scat_index,
    const String&                             unit,
    const Verbosity& )
{
  // Sizes
  const Index nss = scat_meta.nelem();

  // Checks
  if( nss == 0 )
    throw runtime_error( "*scat_meta* is empty!" );
  if( nss < scat_index+1 )
    {
      ostringstream os;
      os << "Selected scattering species index is " << scat_index << " but this "
         << "is not allowed since *scat_meta* has only" << scat_meta.nelem()
         << " elements.";
      throw runtime_error(os.str());
    }

  // Create size grid
  //
  const Index nse = scat_meta[scat_index].nelem();
  //
  pnd_size_grid.resize( nse );
  //
  if( unit == "dveq" )
    {
      for ( Index i=0; i<nse; i++ )
        { pnd_size_grid[i] = scat_meta[scat_index][i].diameter_volume_equ; }
    }
  else if( unit == "dmax" )
    {
      for ( Index i=0; i<nse; i++ )
        { pnd_size_grid[i] = scat_meta[scat_index][i].diameter_max; }
    }
  else if( unit == "mass" )
    {
      for ( Index i=0; i<nse; i++ )
        { pnd_size_grid[i] = scat_meta[scat_index][i].mass; }
    }
  else if( unit == "area" )
    {
      for ( Index i=0; i<nse; i++ )
        { pnd_size_grid[i] = scat_meta[scat_index][i].diameter_area_equ_aerodynamical; }
    }
  else
    {
      ostringstream os;
      os << "You have selected the unit: " << unit 
         << "while accepted choices are: \"dveq\", \"dmax\", \"mass\" and \"area\"";
      throw runtime_error(os.str());
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pndFromPsdBasic(
         Matrix&    pnd_data,
         Tensor3&   dpnd_data_dx,
   const Vector&    pnd_size_grid,
   const Matrix&    psd_data,
   const Vector&    psd_size_grid,
   const Tensor3&   dpsd_data_dx,
   const Index&     quad_order,
   const Verbosity& )
{
  // Some sizes 
  const Index np  = psd_data.nrows();
  const Index ng  = psd_size_grid.nelem();
        Index ndx = 0;
  const bool  do_dx = !dpsd_data_dx.empty();

  // Checks
  if( ng < 2 )
    throw runtime_error( "The method requires that length of *psd_size_grid* is >= 2." );
  if( ng != pnd_size_grid.nelem() )
    throw runtime_error( "So far, the method requires that *psd_size_grid* and "
                         "*pnd_size_grid* have same length." );
  for( Index i=0; i<ng; i++ )
    {
      if( psd_size_grid[i] != pnd_size_grid[i] )
        throw runtime_error( "So far, the method requires that *psd_size_grid* and "
                             "*pnd_size_grid* are identical." );
    }
  if( psd_data.ncols() != ng )
    throw runtime_error( "Number of columns in *psd_data* and length of "
                         "*psd_size_grid* must match." );

  pnd_data.resize( np, ng );
  if( do_dx )
    {
      if( dpsd_data_dx.ncols() != ng )
        throw runtime_error( "Number of columns in *dpsd_data_dx* and length of "
                             "*psd_size_grid* must match." );
      ndx = dpsd_data_dx.npages();
      dpnd_data_dx.resize( ndx, np, ng );
    }
  else
    { dpnd_data_dx.resize( 0, 0, 0 ); }


  // Get sorted version of psd_size_grid (and, since pnd_size_grid so far is
  // identical, of this as well implicitly)
  ArrayOfIndex intarr;
  Vector psd_size_grid_sorted(ng);
  get_sorted_indexes(intarr, psd_size_grid);
  for( Index i=0; i<ng; i++ )
    psd_size_grid_sorted[i] = psd_size_grid[intarr[i]];
    
      
  // Calculate pnd by intrgation of psd for given nodes/bins
  Vector quadweights( ng );
  bin_quadweights( quadweights, psd_size_grid_sorted, quad_order );

  for ( Index i=0; i<ng; i++ )
    {
      for( Index ip=0; ip<np; ip++ )
        { pnd_data(ip,intarr[i]) = quadweights[i] * psd_data(ip,intarr[i]); }

      if( do_dx )
        {
          for( Index ip=0; ip<np; ip++ )
            {
              for ( Index ix=0; ix<ndx; ix++ )
                { dpnd_data_dx(ix,ip,intarr[i]) = quadweights[i] *
                                          dpsd_data_dx(ix,ip,intarr[i]); }
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pndFromPsd(
         Matrix&    pnd_data,
         Tensor3&   dpnd_data_dx,
   const Vector&    pnd_size_grid,
   const Matrix&    psd_data,
   const Vector&    psd_size_grid,
   const Tensor3&   dpsd_data_dx,
   const ArrayOfArrayOfSingleScatteringData& scat_data,
   const Vector&    f_grid,
   const Index&     scat_data_checked,
   const Index&     quad_order,
   const Index&     scat_index,
   const Numeric&   threshold_rsec,
   const Numeric&   threshold_bext,
   const Numeric&   threshold_rpnd,
   const Verbosity& )
{
  // Some sizes 
  const Index np  = psd_data.nrows();
  const Index ng  = psd_size_grid.nelem();
        Index ndx = 0;
  const bool  do_dx = !dpsd_data_dx.empty();

  // Checks
  if( ng < 3 )
    throw runtime_error( "The method requires that length of *psd_size_grid* is >= 3." );
  if( ng != pnd_size_grid.nelem() )
    throw runtime_error( "So far, the method requires that *psd_size_grid* and"
                         " *pnd_size_grid* have the same length." );
  for( Index i=0; i<ng; i++ )
    {
      if( psd_size_grid[i] != pnd_size_grid[i] )
        throw runtime_error( "So far, the method requires that *psd_size_grid* and"
                             " *pnd_size_grid* are identical." );
    }
  if( psd_data.ncols() != ng )
    throw runtime_error( "Number of columns in *psd_data* and length of"
                         " *psd_size_grid* must match." );

  pnd_data.resize( np, ng );
  if( do_dx )
    {
      if( dpsd_data_dx.ncols() != ng )
        throw runtime_error( "Number of columns in *dpsd_data_dx* and length of"
                             " *psd_size_grid* must match." );
      ndx = dpsd_data_dx.npages();
      dpnd_data_dx.resize( ndx, np, ng );
    }
  else
    { dpnd_data_dx.resize( 0, 0, 0 ); }

  if( !scat_data_checked )
    throw runtime_error( "*scat_data* must have passed a consistency check"
                         " (scat_data_checked=1).\n"
                         "Alternatively, use *pndFromPsdBasic*." );
  if( scat_index >= scat_data.nelem() )
    throw runtime_error( "*scat_index* exceeds the number of available"
                         " scattering species." );
  if( scat_data[scat_index].nelem() != ng )
    throw runtime_error( "Number of scattering elements in this scattering"
                         " species (*scat_index*) inconsistent with length of"
                         " *pnd_size_grid*." );

  // Get sorted version of psd_size_grid (and, since pnd_size_grid so far is
  // identical, of this as well implicitly)
  ArrayOfIndex intarr;
  Vector psd_size_grid_sorted(ng);
  get_sorted_indexes(intarr, psd_size_grid);
  for( Index i=0; i<ng; i++ )
    psd_size_grid_sorted[i] = psd_size_grid[intarr[i]];

  // Calculate pnd by integration of psd for given nodes/bins
  Vector quadweights( ng );
  bin_quadweights( quadweights, psd_size_grid_sorted, quad_order );

  for ( Index i=0; i<ng; i++ ) //loop over pnd_size_grid aka scattering elements
    {
      for( Index ip=0; ip<np; ip++ ) //loop over pressure levels
        { pnd_data(ip,intarr[i]) = quadweights[i] * psd_data(ip,intarr[i]); }

      if( do_dx )
        {
          for( Index ip=0; ip<np; ip++ )
            {
              for ( Index ix=0; ix<ndx; ix++ )
                { dpnd_data_dx(ix,ip,intarr[i]) = quadweights[i] *
                                          dpsd_data_dx(ix,ip,intarr[i]); }
            }
        }
    }

  ArrayOfSingleScatteringData  sds = scat_data[scat_index];
  Index fstart = 0;
  Index nf = f_grid.nelem();
  Matrix bulkext(np, nf, 0.);
  Vector ext(nf), ext_s0(nf), ext_s1(nf), ext_l0(nf), ext_l1(nf);

  // check that pnd_size_grid and scat_data cover (scalar) bulk extinction properly.
  // (formerly we used mass. but heavy particles might not necessarily
  // contribute much to optprops. same is valid for total particle number.)
  //
  // We do that by 
  // (a) checking that psd*ext decreases at the edges, i.e. to increase the
  //     probability that the main extinction peak is covered and in order to
  //     avoid passing check (b) through making the edge bins tiny
  // (b) checking that the extinction contributions of the edge bins is small
  //     compared to the total bulk bin.
  // The tests are somewhat over-sensitive on the small particle side in the
  // sense that the checks are triggered easily by low mass density cases - which
  // have their main extinction contributed by small particles, but for which
  // the total extinction is quite low. Hence, there is also a total extinction
  // threshold, below which the checks are not applied. Furthermore, such low
  // densities do typically not occur throughout the whole atmospheric profile,
  // but higher density layers might exists, meaning the low density layers
  // contribute very little to the total optical properties, hence should not
  // dominate the check criteria (no need to be very strict with them if they
  // hardly contribute anyways). Therefore, checks are also not applied if the
  // edge bin number density at a given atmospheric level is low compared to the
  // maximum number of this scattering element in the profile.
  //
  // Technically, we need to cover all freqs separately (as optprops vary
  // strongly with freq). We indeed do that here.
  // FIXME: Shall a single freq or reduced freq-number option be implemented?
  // If so, the don't-apply-test criteria need to be checked again though (so it
  // does not happen that effectively all of the reduced-freq testing is skipped
  // due to them.
  //
  // Technically, we also need to use identical temperatures (and preferably the
  // ones used for deriving psd). But all scat elements can be on different
  // T-grids, ie we would need to interpolate. Instead we just use the first in
  // line.
  // FIXME: Maybe refine in future allowing to select between low, high, median
  // grid point or even specify one temp to use.
  //
  // And, technically ext might be directional dependent (for oriented
  // particles). Also this might be different for each scat element. However, at
  // least for now, we only check the 0-direction.
  // FIXME: Check further directions?

  for( Index ise=0; ise<ng; ise++ ) //loop over pnd_size_grid aka scattering elements
    {
      if( sds[intarr[ise]].ext_mat_data.nshelves()>1 )
        ext = sds[intarr[ise]].ext_mat_data(joker,0,0,0,0);
      else
        ext = sds[intarr[ise]].ext_mat_data(0,0,0,0,0);

      // keep the ext data of the edge bins 
      if( ise==0 )
        ext_s0 = ext;
      else if( ise==1 )
        ext_s1 = ext;
      else if( ise==ng-2 )
        ext_l1 = ext;
      else if( ise==ng-1 )
        ext_l0 = ext;
    
      for( Index ip=0; ip<np; ip++ ) //loop over pressure levels
        if( abs(pnd_data(ip,joker).sum()) > 0. )
          for( Index f=fstart; f<(fstart+nf); f++ )
            bulkext(ip,f) += pnd_data(ip,intarr[ise]) * ext[f];
    }

  Numeric max0=0, max1=0;
  for( Index ip=0; ip<np; ip++ ) //loop over pressure levels
    {
      max0 = max(abs(pnd_data(ip,intarr[0])),max0);
      max1 = max(abs(pnd_data(ip,intarr[ng-1])),max1);
    }

  Numeric contrib;
  for( Index ip=0; ip<np; ip++ ) //loop over pressure levels
  {
    if( abs(pnd_data(ip,joker).sum()) > 0. )
    {
      for( Index f=fstart; f<(fstart+nf); f++ )
      {
/*        for( Index ise=0; ise<ng; ise++ )
        {
          if( sds[ise].ext_mat_data.nshelves()>1 )
            ext = sds[ise].ext_mat_data(joker,0,0,0,0);
          else
            ext = sds[ise].ext_mat_data(0,0,0,0,0);
          cout << "    se #" << ise << ": contrib = pnd*ext/bext = "
               << abs(pnd_data(ip,ise)) << "*" << ext[f] << "/"
               << abs(bulkext(ip,f)) << " = "
               << 1e2*abs(pnd_data(ip,ise))*ext[f]/abs(bulkext(ip,f))
               << "%\n";
        }*/

        // check that bin-width normalized extinction (or ext*psd) is
        // decreasing.
        if( abs(bulkext(ip,f)) > 1e-2*threshold_bext )
        {
          if( abs(psd_data(ip,intarr[0])) > 0. and
              ext_s0[f]*abs(psd_data(ip,intarr[0])) >=
              ext_s1[f]*abs(psd_data(ip,intarr[1])) )
          {
            ostringstream os;
            os << "  Bin-width normalized extinction (ext*psd) not decreasing"
               << " at small size edge\n"
               << "  at atm level #" << ip
               << " and freq point #" << f << ".\n"
               << "  ext_s0=" << ext_s0[f]
               << ", psd_s0=" << abs(psd_data(ip,intarr[0]))
               << ", ext_s0*psd_s0=" << ext_s0[f]*abs(psd_data(ip,intarr[0]))
               << "\n    LARGER EQUAL\n"
               << "  ext_s1=" << ext_s1[f]
               << ", psd_s1=" << abs(psd_data(ip,intarr[1]))
               << ", ext_s1*psd_s1=" << ext_s1[f]*abs(psd_data(ip,intarr[1])) << "\n"
               << "    Total bulk ext = " << abs(bulkext(ip,f)) << "\n"
               << "  Need to add smaller sized particles!\n";
            throw runtime_error(os.str());
          }

          if( abs(psd_data(ip,intarr[ng-1])) > 0. and
              ext_l0[f]*abs(psd_data(ip,intarr[ng-1])) >=
              ext_l1[f]*abs(psd_data(ip,intarr[ng-2])) )
          {
            ostringstream os;
            os << "Bin-width normalized extinction (ext*psd) not decreasing"
              << " at large size edge\n"
              << "at atm level #" << ip
              << " and freq point #" << f << ".\n"
              << "  ext_l0=" << ext_l0[f]
              << ", psd_l0=" << abs(psd_data(ip,intarr[ng-1]))
              << ", ext_l0*psd_l0=" << ext_l0[f]*abs(psd_data(ip,intarr[ng-1]))
              << "\n    LARGER EQUAL\n"
              << "  ext_l1=" << ext_l1[f]
              << ", psd_l1=" << abs(psd_data(ip,intarr[ng-2]))
              << ", ext_l1*psd_l1=" << ext_l1[f]*abs(psd_data(ip,intarr[ng-2])) << "\n"
               << "    Total bulk ext = " << abs(bulkext(ip,f)) << "\n"
              << "  Need to add larger sized particles!\n";
            throw runtime_error(os.str());
          }
        }

        // check that contribution of edge bins to total extinction is
        // sufficiently small
        if( abs(bulkext(ip,f)) > threshold_bext )
        {
          if( abs(pnd_data(ip,intarr[0])) > threshold_rpnd*max0 )
          {
            contrib = abs(pnd_data(ip,intarr[0]))*ext_s0[f]/abs(bulkext(ip,f));
            //cout << "    small edge contrib = pnd*ext/bext = "
            //     << abs(pnd_data(ip,intarr[0])) << "*" << ext_s0[f] << "/"
            //     << abs(bulkext(ip,f)) << " = " << contrib << "\n";
            if( abs(pnd_data(ip,intarr[0]))*ext_s0[f] > threshold_rsec * abs(bulkext(ip,f)) )
            {
              ostringstream os;
              os << "Contribution of edge bin to total extinction too high"
                 << " (" << contrib*1e2 << "% of " << abs(bulkext(ip,f))
                 << ") at small size edge\n"
                 << "at atm level #" << ip
                 << " and freq point #" << f << ".\n"
                 << "  Need to add smaller sized particles or refine the size"
                 << " grid on the small size edge!\n";
              throw runtime_error(os.str());
            }
          }
          if( abs(pnd_data(ip,intarr[ng-1])) > threshold_rpnd*max1 )
          {
            contrib = abs(pnd_data(ip,intarr[ng-1]))*ext_l0[f]/abs(bulkext(ip,f));
            //cout << "    large edge contrib = pnd*ext/bext = "
            //     << abs(pnd_data(ip,ng-1)) << "*" << ext_l0[f] << "/"
            //     << abs(bulkext(ip,f)) << " = " << contrib << "\n";
            if( abs(pnd_data(ip,intarr[ng-1]))*ext_l0[f] > threshold_rsec * abs(bulkext(ip,f)) )
            {
              ostringstream os;
              os << "Contribution of edge bin to total extinction too high"
                 << " (" << contrib*1e2 << "% of " << abs(bulkext(ip,f))
                 << ") at large size edge\n"
                 << "at atm level #" << ip
                 << " and freq point #" << f << ".\n"
                 << "  Need to add larger sized particles or refine the size"
                 << " grid on the large size edge!\n";
              throw runtime_error(os.str());
            }
          }
        }
      }
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void pndAdjustFromScatMeta(
          Matrix&                             pnd_data,
//          Tensor3&                            dpnd_data_dx,
    const Matrix&                             pnd_agenda_input,
    const ArrayOfString&                      pnd_agenda_input_names,
    const ArrayOfArrayOfScatteringMetaData&   scat_meta,
    const String&                             pnd_agenda_input_tag,
    const Index&                              scat_index,
    const Verbosity&                          verbosity )
{
  const Index np  = pnd_data.nrows();
  if( pnd_agenda_input.nrows() != np )
    throw runtime_error( "Number of rows in *pnd_data* and *pnd_agenda_input*"
                         " must be equal." );

  const Index nin = pnd_agenda_input_names.nelem();
  if( pnd_agenda_input.ncols() != nin )
    throw runtime_error( "Length of *pnd_agenda_input_names* and number of "
                         "columns in *pnd_agenda_input* must be equal." );

  const Index nse  = pnd_data.ncols();
  if( scat_meta[scat_index].nelem() != nse )
    {
      ostringstream os;
      os << "Number of columns in *pnd_data* and elements in *scat_meta*\n"
         << "of species " << scat_index << " must be equal.";
      throw runtime_error(os.str());
    }

  Index not_found = 1;
  Index fieldID = 0;
  while( not_found && (fieldID<nin) )
    {
      if( pnd_agenda_input_names[fieldID]==pnd_agenda_input_tag )
        not_found = 0;
      else
        fieldID += 1;
    }
  if( not_found )
    {
      ostringstream os;
      os << "Requested adjustment field \"" << pnd_agenda_input_tag
         << "\" not found in *pnd_agenda_input_names*.";
      throw runtime_error(os.str());
    }

  Vector mass( nse );
  // we have to assume correspondence of scattering element order in pnd_data
  // and scat_meta[scat_index]. (any way to ensure that?)
  for ( Index i=0; i<nse; i++ )
    {
      if ( isnan(scat_meta[scat_index][i].mass) )
        {
          ostringstream os;
          os << "No mass data available for scattering element #"
             << i << " of scattering species with index "
             << scat_index << ".";
          throw runtime_error( os.str() );
        }
      mass[i] = scat_meta[scat_index][i].mass;
    }

  for( Index ip=0; ip<np; ip++ )
    {
      // scaling pnd to real mass/number density/flux (some PSDs have implicit
      // scaling - then this should only a minor adjustment, rather a check -,
      // but others don't
      Vector pnd = pnd_data(ip,joker);
      chk_pndsum ( pnd, pnd_agenda_input(ip,fieldID), mass,
                   ip, -1, -1, pnd_agenda_input_tag, verbosity );
      pnd_data(ip,joker) = pnd;
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MassSizeParamsFromScatMeta(
          Numeric&                            alpha,
          Numeric&                            beta,
    const ArrayOfArrayOfScatteringMetaData&   scat_meta,
    const Index&                              scat_index,
    const Numeric&                            D_min,
    const Numeric&                            D_max,
    const Numeric&                            beta_default,
    const Verbosity& )
{
  // Sizes
  const Index nss = scat_meta.nelem();

  // Checks
  if( nss == 0 )
    throw runtime_error( "*scat_meta* is empty!" );
  if( nss < scat_index+1 )
    {
      ostringstream os;
      os << "Selected scattering species index is " << scat_index
         << " but *scat_meta* contains only " << scat_meta.nelem()
         << " species.";
      throw runtime_error(os.str());
    }

  // Extract size and mass data from scat_meta
  const Index nse = scat_meta[scat_index].nelem();

  if( nse==0 )
    {
      ostringstream os;
      os << "*scat_meta* of selected scattering species with index "
         << scat_index << " contains has no scattering elements.";
      throw runtime_error(os.str());
    }

  if( nse>1 )
    {
      ArrayOfIndex intarr_sort, intarr_unsort( 0 );
      Vector dmax_unsorted( nse );
      Vector q;
      Index nsev=0;

      for ( Index i=0; i<nse; i++ )
      {
        if ( isnan(scat_meta[scat_index][i].diameter_max) )
        {
            ostringstream os;
            os << "No maximum diameter data available for scattering element #"
               << i << " of scattering species with index " << scat_index
               << ".";
            throw runtime_error( os.str() );
        }
        if( scat_meta[scat_index][i].diameter_max>=D_min &&
            scat_meta[scat_index][i].diameter_max<=D_max )
          {
            dmax_unsorted[nsev] = scat_meta[scat_index][i].diameter_max;
            intarr_unsort.push_back( i );
            nsev += 1;
          }
      }

      get_sorted_indexes(intarr_sort, dmax_unsorted[Range(0,nsev)]);
      Vector dmax( nsev ), log_D( nsev );
      Vector mass( nsev ), log_m( nsev );
    
      // extract scattering meta data
      for ( Index i=0; i<nsev; i++ )
      {
        dmax[i] = dmax_unsorted[intarr_sort[i]]; // [m]
        
        if ( isnan(scat_meta[scat_index][intarr_unsort[intarr_sort[i]]].mass) )
        {
            ostringstream os;
            os << "No mass data available for scattering element #"
               << intarr_unsort[intarr_sort[i]]
               << " of scattering species with index " << scat_index << ".";
            throw runtime_error( os.str() );
        }
        mass[i] = scat_meta[scat_index][intarr_unsort[intarr_sort[i]]].mass; // [kg]
        
        // logarithm of Dmax, needed for estimating mass-dimension-relationship
        log_D[i] = log(dmax[i]);
        
        // logarithm of the mass, even though it is  little weird to have
        // a logarithm of something with a unit...
        log_m[i] = log(mass[i]);
      }
    
      // estimate mass-dimension relationship from meta data by linear regression
      // Assumption of a power law for the mass dimension relationship
      // Approach: log(m) = log(alpha)+beta*log(dmax/D0)
      linreg(q,log_D, log_m);
      alpha = exp(q[0]);
      beta = q[1];
    }
  else
    {
      // for a monodispersion we can't estimate the m-D relation (1 relation, 2
      // unknowns), hence we fix one of the parameters and calculate the other
      // such that we have them consistent. but shouldn't make any difference on
      // the end result whatever we choose here (all ice has to end up with this
      // scattering anyways)
      beta = beta_default;
      if ( isnan(scat_meta[scat_index][0].diameter_max) ||
           isnan(scat_meta[scat_index][0].mass) )
        {
          ostringstream os;
          os << "Maximum diameter or mass data missing for scattering element #"
               << 0 << " of scattering species with index " << scat_index
               << ".";
          throw runtime_error( os.str() );
        }
      alpha = scat_meta[scat_index][0].mass / 
              pow(scat_meta[scat_index][0].diameter_max,beta);
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void psdF07 (
          Matrix&                             psd_data,
          Tensor3&                            dpsd_data_dx,
    const Vector&                             psd_size_grid,
    const Vector&                             pnd_agenda_input_t,
    const Matrix&                             pnd_agenda_input,
    const ArrayOfString&                      pnd_agenda_input_names,
    const ArrayOfString&                      dpnd_data_dx_names,
    const String&                             regime,
    const Numeric&                            alpha,
    const Numeric&                            beta,
    const Numeric&                            t_min,
    const Numeric&                            t_max,
    const Numeric&                            t_min_psd,
    const Numeric&                            t_max_psd,
    const Numeric&                            beta_min, 
    const Numeric&                            beta_max, 
    //const Numeric&                            density_deviation,
    const Index&                              picky, 
    const Verbosity&)
{
  // Some sizes 
  const Index nin = pnd_agenda_input_names.nelem();
  const Index ndx = dpnd_data_dx_names.nelem();
  const Index np  = pnd_agenda_input.nrows();
  const Index nsi = psd_size_grid.nelem();
  
  // Checks
  if( pnd_agenda_input.ncols() != nin )
    throw runtime_error( "Length of *pnd_agenda_input_names* and number of "
                         "columns in *pnd_agenda_input* must be equal." );
  if( pnd_agenda_input.ncols() != 1 )
    throw runtime_error( "*pnd_agenda_input* must have two columns (SWC)." );
  if( pnd_agenda_input_names[0] != "SWC" )
    throw runtime_error( "Element 0 of *pnd_agenda_input_names* must be \"SWC\"." );
  if( ndx )
    {
      if( ndx != 1 )
        throw runtime_error( "*dpnd_data_dx_names* must have length 0 or 1." );
      if( dpnd_data_dx_names[0] != "SWC" )
        throw runtime_error( "With F07, the only valid option for "
                             "*dpnd_data_dx_names* is: \"SWC\"." );
    }
  if( regime!="TR" && regime!="ML" )
    throw runtime_error( "regime must either be \"TR\" or \"ML\"." );
  
  // Init psd_data and dpsd_data_dx with zeros
  psd_data.resize( np, nsi );
  psd_data = 0.0;
  if( ndx )
    {
      dpsd_data_dx.resize( 1, np, nsi );   // SWC only possible retrieval quantity
      dpsd_data_dx = 0.0;
    }
  else
    { dpsd_data_dx.resize( 0, 0, 0  ); }  

  // Some sanity checks. If not picky, we just ignore the insanities.
  if( picky )
    {
      // Sane beta?
      if( beta < beta_min  ||  beta > beta_max )
        {
          ostringstream os;
          os << "Method called with a mass-dimension-relation exponent beta of "
             << beta << ".\n"
             << "This is outside the specified allowed range: ["
             << beta_min << "," << beta_max << "]";
          throw runtime_error(os.str());
        }

      // Do alpha and beta result in sane mass for ice particles (i.e. density
      // shouldn't be > DENSITY_ICE? we check the smallest and largest
      // particles.

      /* that seems unpractical. for beta<3., small particles will frequently
         exceed ice density quite significantly. hence, don't apply this.

      // first, find the smallest and largest particles D_small and D_large (as
      // psd_size grid is not necessarily sorted).
      Vector D_lims{1e3,-1.}; // 1km should be large enough to find
                                         // one that is smaller...
      for( Index iD=0; iD<nsi; iD++ )
        {
          if( D_lims[0] > psd_size_grid[iD] ) D_lims[0]=psd_size_grid[iD];
          if( D_lims[1] < psd_size_grid[iD] ) D_lims[1]=psd_size_grid[iD];
        }

      // now calc resulting densities for solid spheres D_small and D_large
      // diameters
      for( Index iD=0; iD<2; iD++ )
        {
          Numeric rho = 6.*alpha/PI * pow(D_lims[iD],(beta-3.));
          if( rho*(1.+density_deviation) > DENSITY_OF_ICE )
            {
              ostringstream os;
              os << "Method called with a mass-dimension-relation parameters"
                 << " alpha=" << alpha << " and beta=" << beta << "\n"
                 << "results in minimum density of " << rho << ".\n"
                 << "This is unphysical for ice particles and beyond the"
                 << " allowed unphysicality threshold: "
                 << rho*(1.+density_deviation) << ".";
              throw runtime_error(os.str());
            }
        }
      */
    }

  for( Index ip=0; ip<np; ip++ )
    {
      
      // Extract the input variables
      Numeric swc = pnd_agenda_input(ip,0);
      Numeric   t = pnd_agenda_input_t[ip];

      // No calc needed if swc==0 and no jacobians requested.
      if( (swc==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!

      // Outside of [t_min,tmax]?
      if( t < t_min  ||  t > t_max )
        {
          if( picky )
            {
              ostringstream os;
              os << "Method called with a temperature of " << t << " K.\n"
                 << "This is outside the specified allowed range: [ max(0.,"
                 << t_min << "), " << t_max << " ]";
              throw runtime_error(os.str());
            }
          else  
            { continue; }   // If here, we are ready with this point!
        }

      // PSD assumed to be constant outside [*t_min_psd*,*t_max_psd*]
      if( t < t_min_psd )
        { t = t_min_psd; }
      else if( t > t_max_psd )
        { t = t_max_psd; }
    
      // Negative swc?
      Numeric psd_weight = 1.0;
      if( swc < 0 )
        {
          psd_weight = -1.0;
          swc       *= -1.0;
        }
      
      // Calculate PSD
      Vector psd_1p(nsi);
      if( swc != 0 )
        {
          psd_snow_F07 ( psd_1p, psd_size_grid, swc, t, alpha, beta, regime );
          for ( Index i=0; i<nsi; i++ )
            { psd_data(ip,i) = psd_weight * psd_1p[i]; }
        }

      // Calculate derivative with respect to IWC
      if( ndx )
        {
          // Obtain derivative by perturbation of 0.1%, but not less than 0.1 mg/m3.
          // Note that the latter value becomes the perturbation for IWC=0.
          const Numeric dswc = max( 0.001*swc, 1e-7 );
          const Numeric swcp = swc + dswc;
          psd_snow_F07 ( psd_1p, psd_size_grid, swcp, t, alpha, beta, regime );
          for ( Index i=0; i<nsi; i++ )
            { dpsd_data_dx(0,ip,i) = ( psd_1p[i] - psd_data(ip,i) ) / dswc; }
        }   
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void psdMH97 (
          Matrix&                             psd_data,
          Tensor3&                            dpsd_data_dx,
    const Vector&                             psd_size_grid,
    const Vector&                             pnd_agenda_input_t,
    const Matrix&                             pnd_agenda_input,
    const ArrayOfString&                      pnd_agenda_input_names,
    const ArrayOfString&                      dpnd_data_dx_names,
    const Numeric&                            t_min, 
    const Numeric&                            t_max, 
    const Numeric&                            t_min_psd, 
    const Numeric&                            t_max_psd, 
    const Index&                              picky, 
    const Index&                              noisy,
    const Verbosity&)
{
  // Some sizes 
  const Index nin = pnd_agenda_input_names.nelem();
  const Index ndx = dpnd_data_dx_names.nelem();
  const Index np  = pnd_agenda_input.nrows();
  const Index nsi = psd_size_grid.nelem();
  
  // Checks
  if( pnd_agenda_input.ncols() != nin )
    throw runtime_error( "Length of *pnd_agenda_input_names* and number of "
                         "columns in *pnd_agenda_input* must be equal." );
  if( pnd_agenda_input.ncols() != 1 )
    throw runtime_error( "*pnd_agenda_input* must have one column (IWC)." );
  if( pnd_agenda_input_names[0] != "IWC" )
    throw runtime_error( "Element 0 of *pnd_agenda_input_names* must be \"IWC\"." );

  if( ndx )
    {
      if( ndx != 1 )
        throw runtime_error( "*dpnd_data_dx_names* must have length 0 or 1." );
      if( dpnd_data_dx_names[0] != "IWC" )
        throw runtime_error( "With MH97, the only valid option for "
                             "*dpnd_data_dx_names* is: \"IWC\"." );
    }        
  if( noisy  &&   ndx )
    throw runtime_error( "Jacobian calculations and \"noisy\" can not be "
                         "combined." );

  
  // Init psd_data and dpsd_data_dx with zeros
  psd_data.resize( np, nsi );
  psd_data = 0.0;
  if( ndx )
    {
      dpsd_data_dx.resize( 1, np, nsi );   // IWC only possible retrieval quantity
      dpsd_data_dx = 0.0;
    }
  else
    { dpsd_data_dx.resize( 0, 0, 0  ); }  

  for( Index ip=0; ip<np; ip++ )
    {
      
      // Extract the input variables
      Numeric iwc = pnd_agenda_input(ip,0);
      Numeric   t = pnd_agenda_input_t[ip];

      // No calc needed if iwc==0 and no jacobians requested.
      if( (iwc==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!

      // Outside of [t_min,tmax]?
      if( t < t_min  ||  t > t_max )
        {
          if( picky )
            {
              ostringstream os;
              os << "Method called with a temperature of " << t << " K.\n"
                 << "This is outside the specified allowed range: [ max(0.,"
                 << t_min << "), " << t_max << " ]";
              throw runtime_error(os.str());
            }
          else  
            { continue; }   // If here, we are ready with this point!
        }

      // PSD assumed to be constant outside [*t_min_psd*,*t_max_psd*]
      if( t < t_min_psd )
        { t = t_min_psd; }
      else if( t > t_max_psd )
        { t = t_max_psd; }
  
      // Negative iwc?
      Numeric psd_weight = 1.0;
      if( iwc < 0 )
        {
          psd_weight = -1.0;
          iwc       *= -1.0;
        }
      
      // Calculate PSD
      Vector psd_1p(nsi);
      if( iwc != 0 )
        {
          psd_cloudice_MH97 ( psd_1p, psd_size_grid, iwc, t, noisy );
          for ( Index i=0; i<nsi; i++ )
            { psd_data(ip,i) = psd_weight * psd_1p[i]; }
        }

      // Calculate derivative with respect to IWC
      if( ndx )
        {
          // Obtain derivative by perturbation of 0.1%, but not less than 0.1 mg/m3.
          // Note that the latter value becomes the perturbation for IWC=0.
          const Numeric diwc = max( 0.001*iwc, 1e-7 );
          const Numeric iwcp = iwc + diwc;
          psd_cloudice_MH97 ( psd_1p, psd_size_grid, iwcp, t, noisy );
          for ( Index i=0; i<nsi; i++ )
            { dpsd_data_dx(0,ip,i) = ( psd_1p[i] - psd_data(ip,i) ) / diwc; }
        }   
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void psdS2M (
             Matrix&                             psd_data,
             Tensor3&                            dpsd_data_dx,
             const Vector&                             psd_size_grid,
             const Vector&                             pnd_agenda_input_t,
             const Matrix&                             pnd_agenda_input,
             const ArrayOfString&                      pnd_agenda_input_names,
             const ArrayOfString&                      dpnd_data_dx_names,
             const String&                             hydrometeor_type,
             const Numeric&                            t_min,
             const Numeric&                            t_max,
             const Index&                              picky,
             const Verbosity&)
{
    // Some sizes
    const Index nin = pnd_agenda_input_names.nelem();
    const Index ndx = dpnd_data_dx_names.nelem();
    const Index np  = pnd_agenda_input.nrows();
    const Index nsi = psd_size_grid.nelem();
    
    // Checks
    if( pnd_agenda_input.ncols() != nin )
    {
        throw runtime_error( "Length of *pnd_agenda_input_names* and number of "
                            "columns in *pnd_agenda_input* must be equal." );
    }
    if( pnd_agenda_input.ncols() != 2 )
    {
        throw runtime_error( "*pnd_agenda_input* must have two columns"
                            "(mass density and number density)." );
    }
    
    if( ndx > 2 )
    {
        throw runtime_error( "*dpnd_data_dx_names* must have length <=2." );
    }
    
    
    // check name tags
    ArrayOfIndex input_idx={-1,-1};
    
    for (Index i = 0; i<nin; i++)
    {
        if ((Index) pnd_agenda_input_names[i].find("mass_density")!=String::npos)
        {
            input_idx[0]=i; //mass density index
        }
        else if ((Index) pnd_agenda_input_names[i].find("number_density")!=String::npos)
        {
            input_idx[1]=i; //number density index
        }
    }
    
    
    if (input_idx[0]==-1)
    {
        throw runtime_error( "mass_density-tag not found " );
    }
    if (input_idx[1]==-1)
    {
        throw runtime_error( "number_density-tag not found " );
    }
    
    
    // look after jacobian tags
    ArrayOfIndex dpnd_data_dx_idx={-1,-1};
    
    for (Index i = 0; i<ndx; i++)
    {
        if ((Index) dpnd_data_dx_names[i].find("mass_density")!=String::npos)
        {
            dpnd_data_dx_idx[0]=i; //mass density index
        }
        else if ((Index) dpnd_data_dx_names[i].find("number_density")!=String::npos)
        {
            dpnd_data_dx_idx[1]=i; //number density index
        }
    }
    
    
    
    
    
    //        if( dpnd_data_dx_names[0] != "SWC" )
    //            throw runtime_error( "With F07, the only valid option for "
    //                                "*dpnd_data_dx_names* is: \"SWC\"." );
    //    }
    
    // Init psd_data and dpsd_data_dx with zeros
    psd_data.resize( np, nsi );
    psd_data = 0.0;
    if( ndx!=0 )
    {
        dpsd_data_dx.resize( ndx, np, nsi );
        dpsd_data_dx = 0.0;
    }
    else
    { dpsd_data_dx.resize( 0, 0, 0  ); }
    
    
    for( Index ip=0; ip<np; ip++ )
    {
        
        // Extract the input variables
        Numeric WC = pnd_agenda_input(ip,input_idx[0]);
        Numeric N_tot = pnd_agenda_input(ip,input_idx[1]);
        Numeric   t = pnd_agenda_input_t[ip];
        
        // No calc needed if swc==0 and no jacobians requested.
        if( (WC==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!
        
        // Outside of [t_min,tmax]?
        if( t < t_min  ||  t > t_max )
        {
            if( picky )
            {
                ostringstream os;
                os << "Method called with a temperature of " << t << " K.\n"
                << "This is outside the specified allowed range: [ max(0.,"
                << t_min << "), " << t_max << " ]";
                throw runtime_error(os.str());
            }
            else
            { continue; }   // If here, we are ready with this point!
        }
        
        
        // Negative swc?
        Numeric psd_weight = 1.0;
        if( WC < 0 )
        {
            psd_weight = -1.0;
            WC       *= -1.0;
        }
        
        // Calculate PSD and derivatives
        Vector psd_1p(nsi);
        Matrix dpsd_1p(nsi,2);
        if( WC>0  )
        {
            psd_S2M ( psd_1p,dpsd_1p, psd_size_grid, N_tot, WC, hydrometeor_type );
            
            for ( Index i=0; i<nsi; i++ )
            {
                psd_data(ip,i) = psd_weight * psd_1p[i];
                
                
                for (Index idx=0; idx<dpnd_data_dx_idx.nelem(); idx++)
                {
                    // with respect to WC
                    
                    if (dpnd_data_dx_idx[idx]!=-1)
                    {
                        dpsd_data_dx(dpnd_data_dx_idx[idx],ip,i )=psd_weight *dpsd_1p(i,idx);
                    }
                    
                }
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void psdMY2 (
             Matrix&                             psd_data,
             Tensor3&                            dpsd_data_dx,
             const Vector&                             psd_size_grid,
             const Vector&                             pnd_agenda_input_t,
             const Matrix&                             pnd_agenda_input,
             const ArrayOfString&                      pnd_agenda_input_names,
             const ArrayOfString&                      dpnd_data_dx_names,
             const String&                             hydrometeor_type,
             const Numeric&                            t_min,
             const Numeric&                            t_max,
             const Index&                              picky,
             const Verbosity&)
{
    // Some sizes
    const Index nin = pnd_agenda_input_names.nelem();
    const Index ndx = dpnd_data_dx_names.nelem();
    const Index np  = pnd_agenda_input.nrows();
    const Index nsi = psd_size_grid.nelem();
    
    // Checks
    if( pnd_agenda_input.ncols() != nin )
    {
        throw runtime_error( "Length of *pnd_agenda_input_names* and number of "
                            "columns in *pnd_agenda_input* must be equal." );
    }
    if( pnd_agenda_input.ncols() != 2 )
    {
        throw runtime_error( "*pnd_agenda_input* must have two columns"
                            "(mass density and number density)." );
    }
    
    if( ndx > 2 )
    {
        throw runtime_error( "*dpnd_data_dx_names* must have length <=2." );
    }
    
    
    // check name tags
    ArrayOfIndex input_idx={-1,-1};
    
    for (Index i = 0; i<nin; i++)
    {
        if ((Index) pnd_agenda_input_names[i].find("mass_density")!=String::npos)
        {
            input_idx[0]=i; //mass density index
        }
        else if ((Index) pnd_agenda_input_names[i].find("number_density")!=String::npos)
        {
            input_idx[1]=i; //number density index
        }
    }
    
    
    if (input_idx[0]==-1)
    {
        throw runtime_error( "mass_density-tag not found " );
    }
    if (input_idx[1]==-1)
    {
        throw runtime_error( "number_density-tag not found " );
    }
    
    
    // look after jacobian tags
    ArrayOfIndex dpnd_data_dx_idx={-1,-1};
    
    for (Index i = 0; i<ndx; i++)
    {
        if ((Index) dpnd_data_dx_names[i].find("mass_density")!=String::npos)
        {
            dpnd_data_dx_idx[0]=i; //mass density index
        }
        else if ((Index) dpnd_data_dx_names[i].find("number_density")!=String::npos)
        {
            dpnd_data_dx_idx[1]=i; //number density index
        }
    }
    
    
    
    
    
    //        if( dpnd_data_dx_names[0] != "SWC" )
    //            throw runtime_error( "With F07, the only valid option for "
    //                                "*dpnd_data_dx_names* is: \"SWC\"." );
    //    }
    
    // Init psd_data and dpsd_data_dx with zeros
    psd_data.resize( np, nsi );
    psd_data = 0.0;
    if( ndx!=0 )
    {
        dpsd_data_dx.resize( ndx, np, nsi );
        dpsd_data_dx = 0.0;
    }
    else
    { dpsd_data_dx.resize( 0, 0, 0  ); }
    
    
    for( Index ip=0; ip<np; ip++ )
    {
        
        // Extract the input variables
        Numeric WC = pnd_agenda_input(ip,input_idx[0]);
        Numeric N_tot = pnd_agenda_input(ip,input_idx[1]);
        Numeric   t = pnd_agenda_input_t[ip];
        
        // No calc needed if wc==0 and no jacobians requested.
        if( (WC==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!
        
        // Outside of [t_min,tmax]?
        if( t < t_min  ||  t > t_max )
        {
            if( picky )
            {
                ostringstream os;
                os << "Method called with a temperature of " << t << " K.\n"
                << "This is outside the specified allowed range: [ max(0.,"
                << t_min << "), " << t_max << " ]";
                throw runtime_error(os.str());
            }
            else
            { continue; }   // If here, we are ready with this point!
        }
        
        
        // Negative wc?
        Numeric psd_weight = 1.0;
        if( WC < 0 )
        {
            psd_weight = -1.0;
            WC       *= -1.0;
        }
        
        
        // Calculate PSD and derivatives
        Vector psd_1p(nsi);
        Matrix dpsd_1p(nsi,2);
        if( WC>0  )
        {
            psd_MY2 ( psd_1p,dpsd_1p, psd_size_grid, N_tot, WC, hydrometeor_type );
            
            for ( Index i=0; i<nsi; i++ )
            {
                psd_data(ip,i) = psd_weight * psd_1p[i];
                
                
                for (Index idx=0; idx<dpnd_data_dx_idx.nelem(); idx++)
                {
                    // with respect to WC
                    
                    if (dpnd_data_dx_idx[idx]!=-1)
                    {
                        dpsd_data_dx(dpnd_data_dx_idx[idx],ip,i )=psd_weight *dpsd_1p(i,idx);
                    }
                    
                }
            }
        }
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void psdW16 (
          Matrix&                             psd_data,
          Tensor3&                            dpsd_data_dx,
    const Vector&                             psd_size_grid,
    const Vector&                             pnd_agenda_input_t,
    const Matrix&                             pnd_agenda_input,
    const ArrayOfString&                      pnd_agenda_input_names,
    const ArrayOfString&                      dpnd_data_dx_names,
    const Numeric&                            t_min, 
    const Numeric&                            t_max, 
    const Index&                              picky, 
    const Verbosity&)
{
  // Some sizes 
  const Index nin = pnd_agenda_input_names.nelem();
  const Index ndx = dpnd_data_dx_names.nelem();
  const Index np  = pnd_agenda_input.nrows();
  const Index nsi = psd_size_grid.nelem();
  
  // Checks
  if( pnd_agenda_input.ncols() != nin )
    throw runtime_error( "Length of *pnd_agenda_input_names* and number of "
                         "columns in *pnd_agenda_input* must be equal." );
  if( pnd_agenda_input.ncols() != 1 )
    throw runtime_error( "*pnd_agenda_input* must have one column." );
  if( pnd_agenda_input_names[0] != "RWC" )
    throw runtime_error( "Element 0 of *pnd_agenda_input_names* must be \"RWC\"." );
  if( ndx )
    {
      if( ndx != 1 )
        throw runtime_error( "*dpnd_data_dx_names* must have length 0 or 1." );
      if( dpnd_data_dx_names[0] != "RWC" )
        throw runtime_error( "With W16, the only valid option for "
                             "*dpnd_data_dx_names* is: \"RWC\"." );
    }        

  
  // Init psd_data and dpsd_data_dx with zeros
  psd_data.resize( np, nsi );
  psd_data = 0.0;
  if( ndx )
    {
      dpsd_data_dx.resize( 1, np, nsi );   // RWC only possible retrieval quantity
      dpsd_data_dx = 0.0;
    }
  else
    { dpsd_data_dx.resize( 0, 0, 0  ); }  

  for( Index ip=0; ip<np; ip++ )
    {
      // Extract the input variables
      Numeric rwc = pnd_agenda_input(ip,0);
      Numeric   t = pnd_agenda_input_t[ip];

      // No calc needed if swc==0 and no jacobians requested.
      if( (rwc==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!

      // Outside of [t_min,tmax]?
      if( t < t_min  ||  t > t_max )
        {
          if( picky )
            {
              ostringstream os;
              os << "Method called with a temperature of " << t << " K.\n"
                 << "This is outside the specified allowed range: [ max(0.,"
                 << t_min << "), " << t_max << " ]";
              throw runtime_error(os.str());
            }
          else  
            { continue; }   // If here, we are ready with this point!
        }
      
      // Negative rwc?
      Numeric psd_weight = 1.0;
      if( rwc < 0 )
        {
          psd_weight = -1.0;
          rwc       *= -1.0;
        }

      // Calculate PSD
      Vector psd_1p(nsi);
      if( rwc != 0 )
        {
          psd_rain_W16 ( psd_1p, psd_size_grid, rwc );
          for ( Index i=0; i<nsi; i++ )
            { psd_data(ip,i) = psd_weight * psd_1p[i]; }
        }

      // Calculate derivative with respect to IWC
      if( ndx )
        {
          // Obtain derivative by perturbation of 0.1%, but not less than 0.1 mg/m3.
          // Note that the latter value becomes the perturbation for RWC=0.
          const Numeric drwc = max( 0.001*rwc, 1e-7 );
          const Numeric rwcp = rwc + drwc;
          psd_rain_W16 ( psd_1p, psd_size_grid, rwcp );
          for ( Index i=0; i<nsi; i++ )
            { dpsd_data_dx(0,ip,i) = ( psd_1p[i] - psd_data(ip,i) ) / drwc; }
        }   
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldCalcFromParticleBulkProps(
         Workspace&                   ws,
         Tensor4&                     pnd_field,
         ArrayOfTensor4&              dpnd_field_dx,
   const Index&                       atmosphere_dim,
   const Vector&                      p_grid,
   const Vector&                      lat_grid,
   const Vector&                      lon_grid,
   const Tensor3&                     t_field,
   const Index&                       cloudbox_on,
   const ArrayOfIndex&                cloudbox_limits,
   const ArrayOfString&               scat_species,
   const ArrayOfArrayOfSingleScatteringData& scat_data,
   const ArrayOfArrayOfScatteringMetaData&   scat_meta,
   const Tensor4&                     particle_bulkprop_field,
   const ArrayOfString&               particle_bulkprop_names,
   const ArrayOfAgenda&               pnd_agenda_array,
   const ArrayOfArrayOfString&        pnd_agenda_array_input_names,
   const Index&                       jacobian_do,
   const ArrayOfRetrievalQuantity&    jacobian_quantities,
   const Verbosity& )
{
  // As we allow this method to be called without cloudbox_checkedCalc, it must
  // contain quite a number of checks.
  
  // Number of scattering species
  const Index nss = scat_data.nelem();

  // Checks (not totally complete, but should cover most mistakes)
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  chk_atm_field( "t_field", t_field, atmosphere_dim, 
                 p_grid, lat_grid, lon_grid );
  chk_atm_field( "particle_bulkprop_field", particle_bulkprop_field,
                 atmosphere_dim, particle_bulkprop_names.nelem(),
                 p_grid, lat_grid, lon_grid );

  // Further checks of *particle_bulkprop_field* below
  if( !cloudbox_on )
  {
    if( jacobian_do )
    {
      // FIXME: we might be able to avoid error throwing, but need to fill
      // dpnd_field_dx properly then, don't we?
      throw runtime_error( "*cloudbox_on* must be true to derive jacobians"
                           " using this method." );
    }
    else
      return;
  }

  if( cloudbox_limits.nelem() != 2*atmosphere_dim )
    throw runtime_error( "Length of *cloudbox_limits* incorrect with respect "
                         "to *atmosphere_dim*." );
  if( cloudbox_limits[1]<=cloudbox_limits[0] || cloudbox_limits[0]<0 ||
                                           cloudbox_limits[1]>=p_grid.nelem() )
    throw runtime_error( "Invalid data in pressure part of *cloudbox_limits*." );
  if( atmosphere_dim > 1 )
    {
      if( cloudbox_limits[3]<=cloudbox_limits[2] || cloudbox_limits[2]<0 ||
                                           cloudbox_limits[3]>=lat_grid.nelem() )
        throw runtime_error( "Invalid data in latitude part of *cloudbox_limits*." );
      if( atmosphere_dim > 2 )
        {
          if( cloudbox_limits[5]<=cloudbox_limits[4] || cloudbox_limits[4]<0 ||
                                           cloudbox_limits[5]>=lon_grid.nelem() )
            throw runtime_error( "Invalid data in longitude part of *cloudbox_limits*." );
        }
    }
  if( nss < 1 )
    throw runtime_error( "*scat_data* is empty!." );
  if( scat_species.nelem() != nss )
    throw runtime_error( "*scat_data* and *scat_species* are inconsistent in size." );
  if( scat_meta.nelem() != nss )
    throw runtime_error( "*scat_data* and *scat_meta* are inconsistent in size." );
  if( pnd_agenda_array.nelem() != nss )
    throw runtime_error( "*scat_data* and *pnd_agenda_array* are inconsistent "
                         "in size." );
  if( pnd_agenda_array_input_names.nelem() != nss )
    throw runtime_error( "*scat_data* and *pnd_agenda_array_input_names* are "
                         "inconsistent in size." );
  // Further checks of scat_data vs. scat_meta below  

  
  // Effective lengths of cloudbox
  const Index np = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index ip_offset = cloudbox_limits[0];
  Index nlat = 1;
  Index ilat_offset = 0;
  if( atmosphere_dim > 1 )
    {
      nlat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
      ilat_offset = cloudbox_limits[2];
    }
  Index nlon = 1;
  Index ilon_offset = 0;
  if( atmosphere_dim > 2 )
    {
      nlat = cloudbox_limits[5] - cloudbox_limits[4] + 1;
      ilon_offset = cloudbox_limits[4];
    }

  // Check that *particle_bulkprop_field* contains zeros outside and at
  // cloudbox boundaries
  const String estring = "*particle_bulkprop_field* can only contain non-zero "
    "values inside the cloudbox.";
  // Pressure end ranges
  for( Index ilon=0; ilon<nlon; ilon++ )
  {
    for( Index ilat=0; ilat<nlat; ilat++ )
    {
        if (cloudbox_limits[0])
        {
            for( Index ip=0; ip<=cloudbox_limits[0]; ip++ )
            {
                if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
                {
                    throw runtime_error( estring );
                }
            }
        }
        for( Index ip=cloudbox_limits[1]; ip<p_grid.nelem(); ip++ )
        {
            if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
            {
                throw runtime_error( estring );
            }
        }
    }
  }
  if( atmosphere_dim > 1 )
    {
      // Latitude end ranges
      for( Index ilon=0; ilon<nlon; ilon++ ) {
        for( Index ip=cloudbox_limits[0]+1; ip<cloudbox_limits[1]-1; ip++ ) {
          for( Index ilat=0; ilat<=cloudbox_limits[2]; ilat++ ) {
            if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
              throw runtime_error( estring ); }
          for( Index ilat=cloudbox_limits[3]; ilat<lat_grid.nelem(); ilat++ ) {
            if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
              throw runtime_error( estring ); } } }
      if( atmosphere_dim > 2 )
        {
          // Longitude end ranges
          for( Index ip=cloudbox_limits[0]+1; ip<cloudbox_limits[1]-1; ip++ ) {
            for( Index ilat=cloudbox_limits[2]+1; ilat<cloudbox_limits[3]-1; ilat++ ) {
              for( Index ilon=0; ilon<=cloudbox_limits[4]; ilon++ ) {
                if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
                  throw runtime_error( estring ); }
              for( Index ilon=cloudbox_limits[5]; ilon<lon_grid.nelem(); ilon++ ) {
                if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
                  throw runtime_error( estring ); } } }
        }
    }
  
  // Cumulative number of scattering elements
  ArrayOfIndex ncumse(nss+1);
  ncumse[0] = 0;
  for( Index i=0; i<nss; i++ )
    {
      if( scat_data[i].nelem() != scat_meta[i].nelem() )
        throw runtime_error( "*scat_data* and *scat_meta* have inconsistent sizes." );
      ncumse[i+1] = ncumse[i] + scat_data[i].nelem();
    }

  // Allocate output variables
  //
  pnd_field.resize( ncumse[nss], np, nlat, nlon );
  pnd_field = 0.0;  // To set all end values to zero
  //
  // Help variables for partial derivatives
  Index                nq = 0;
  ArrayOfArrayOfIndex  scatspecies_to_jq;
  //
  if( !jacobian_do )
    { dpnd_field_dx.resize(0); }
  else
    {
      nq = jacobian_quantities.nelem();
      dpnd_field_dx.resize( nq );
      scatspecies_to_jq.resize( nss );
      //
      for( Index iq=0; iq<nq; iq++ )
        {
          if( jacobian_quantities[iq].MainTag() == SCATSPECIES_MAINTAG  )
            {
              const Index ihit = find_first( scat_species,
                                             jacobian_quantities[iq].Subtag() );
              if( ihit < 0 )
                {
                  ostringstream os;
                  os << "Jacobian quantity with index " << iq << " refers to\n"
                     << "  " << jacobian_quantities[iq].Subtag()
                     << "\nbut this species could not be found in *scat_species*.";
                  throw runtime_error(os.str());
                }
              scatspecies_to_jq[ihit].push_back( iq );
              dpnd_field_dx[iq].resize( ncumse[nss], np, nlat, nlon );
              dpnd_field_dx[iq] = 0.0;  // To set all end values to zero
            }
        }
    }

  // Extract data from pnd-agenda array
  for( Index is=0; is<nss; is++ )
    {
      // Index range with respect to pnd_field
      Range se_range( ncumse[is], ncumse[is+1]-ncumse[is] );

      // Determine how pnd_agenda_array_input_names are related to input fields
      //
      const Index nin = pnd_agenda_array_input_names[is].nelem();
      ArrayOfIndex i_pbulkprop(nin);
      //
      for( Index i=0; i<nin; i++ )
        {
          i_pbulkprop[i] = find_first( particle_bulkprop_names,
                                       pnd_agenda_array_input_names[is][i] );
          if( i_pbulkprop[i] < 0 )
            {
              ostringstream os;
              os << "Pnd-agenda with index " << is << " is set to require \""
                 << pnd_agenda_array_input_names[is][i] << "\",\nbut this quantity "
                 << "could not found in *particle_bulkprop_names*.\n"
                 << "(Note that temperature must be written as \"Temperature\")";
                  throw runtime_error(os.str());
            }
        }

      // Set *dpnd_data_dx_names*
      //
      Index ndx = 0;
      ArrayOfString dpnd_data_dx_names( 0 );
      //
      if( jacobian_do )
        {
          ndx = scatspecies_to_jq[is].nelem();
          dpnd_data_dx_names.resize( ndx );
          for( Index ix=0; ix<ndx; ix++ )
            { dpnd_data_dx_names[ix] =
                jacobian_quantities[scatspecies_to_jq[is][ix]].SubSubtag(); }
        }
      
      // Loop lat/lon positions and call *pnd_agenda*
      for( Index ilon=0; ilon<nlon; ilon++ )
        { 
          for( Index ilat=0; ilat<nlat; ilat++ )
            {
              // Note that we don't need any calculations for end points

              // Here we consider this for lat and lon
              if( ( nlat > 1  &&  ( ilat == 0  ||  ilat == nlat-1 ) )  || 
                  ( nlon > 1  &&  ( ilon == 0  ||  ilon == nlon-1 ) ) )
                { continue; }
                  
              // Pressure handled here, by not including end points in loops

              Matrix pnd_agenda_input( np-2, nin );
              //
              for( Index i=0; i<nin; i++ )
                {
                  for( Index ip=1; ip<np-1; ip++ )
                    { pnd_agenda_input(ip-1,i) = particle_bulkprop_field(
                                                         i_pbulkprop[i],
                                                         ip_offset   + ip,
                                                         ilat_offset + ilat,
                                                         ilon_offset + ilon );
                    }
                }

              Vector pnd_agenda_input_t( np-2 );
              //
              for( Index ip=1; ip<np-1; ip++ )
                { pnd_agenda_input_t[ip-1] = t_field( ip_offset   + ip,
                                                      ilat_offset + ilat,
                                                      ilon_offset + ilon ); }
              
              // Call pnd-agenda array
              Matrix pnd_data;
              Tensor3 dpnd_data_dx;
              //
              pnd_agenda_arrayExecute( ws, pnd_data, dpnd_data_dx, is,
                                       pnd_agenda_input_t, pnd_agenda_input,
                                       pnd_agenda_array_input_names[is],
                                       dpnd_data_dx_names, pnd_agenda_array );

              // Copy to output variables
              for( Index ip=1; ip<np-1; ip++ )
                { pnd_field(se_range,ip,ilat,ilon) = pnd_data(ip-1,joker); }
              for( Index ix=0; ix<ndx; ix++ )
                {
                  for( Index ip=1; ip<np-1; ip++ )
                    { dpnd_field_dx[scatspecies_to_jq[is][ix]]
                            (se_range,ip,ilat,ilon) = dpnd_data_dx(ix,ip-1,joker); }
                }
            }
        }
    }
}
   
                           

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_F07 (//WS Output:
                Vector& dNdD,
                //WS Input:
                const Vector& diameter_max,
                const Numeric& SWC,
                const Numeric& T,
                const String& regime,
                const Numeric& alpha,
                const Numeric& beta,
                const Index& robust,
                const Verbosity&)
{
  Index n_se = diameter_max.nelem();
  dNdD.resize(n_se);
    
  // abort if SWC is negative
  if ( !robust && SWC<0. )
    {
      ostringstream os;
      os << "Snow water content can not be negative."
         << " Yours is " << SWC << "kg/m3.";
      throw runtime_error ( os.str() );
    }
  Numeric swc = max( SWC, 0. );

  // abort if T is negative
  if ( T<0. )
    {
      ostringstream os;
      os << "Negative temperatures not allowed.\n"
         << "Yours is " << T << "K.";
      throw runtime_error ( os.str() );
    }

  // check regime
  if( regime!="TR" && regime!="ML" )
    {
      ostringstream os;
      os << "Only regimes \"TR\" (tropical) and \"ML\" (midlatitude) known.\n"
         << "Yours is \"" << regime << "\".";
      throw runtime_error ( os.str() );
    }

  // calculate particle size distribution with F07TR
  // [# m^-3 m^-1]
  psd_snow_F07(dNdD, diameter_max, swc, T, alpha, beta, regime);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_H11 (//WS Output:
             Vector& dNdD,
             //WS Input:
             const Vector& Dmax,
             const Numeric& t,
             const Verbosity&)
{
  Index n_se = Dmax.nelem();
  dNdD.resize(n_se);

  for ( Index i=0; i<n_se; i++ ) //loop over number of scattering elementss
    {
      // calculate particle size distribution for H11
      // [# m^-3 m^-1]
      dNdD[i] = IWCtopnd_H11 ( Dmax[i], t );
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_H13_Ar (//WS Output:
                Vector& dNdD,
                Vector& Ar,
                //WS Input:
                const Vector& Dmax,
                const Numeric& t,
                const Verbosity&)
{
  Index n_se = Dmax.nelem();
  dNdD.resize(n_se);
  Ar.resize(n_se);

  for ( Index i=0; i<n_se; i++ ) //loop over number of scattering elementss
    {
      // calculate particle size distribution for H13Shape
      // [# m^-3 m^-1]
      dNdD[i] = IWCtopnd_H13Shape ( Dmax[i], t );
      // calculate Area ratio distribution for H13Shape
      Ar[i] = area_ratioH13 ( Dmax[i], t );
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_H98 (//WS Output:
             Vector& dNdD,
             //WS Input:
             const Vector& diameter_volume_equivalent,
             const Numeric& LWC,
             const Verbosity&)
{
  Index n_se = diameter_volume_equivalent.nelem();
  dNdD.resize(n_se);
  const Numeric dDdR = 2.; 

  for ( Index i=0; i<n_se; i++ )
    {
      // calculate particle size distribution for liquid
      // and compensate for LWCtopnd providing dNdR
      // [# m^-3 m^-1]
      dNdD[i] = LWCtopnd ( LWC, diameter_volume_equivalent[i]/2. ) / dDdR;
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_MGD_IWC (//WS Output:
                   Vector& dNdD,
                   //WS Input:
                   const Vector& diameter_volume_equ,
                   const Numeric& rho,
                   const Numeric& IWC,
                   const Verbosity&)
{
    Index n_se = diameter_volume_equ.nelem();
    dNdD.resize(n_se);
    
    for ( Index i=0; i<n_se; i++ )
    {
        // calculate particle size distribution with ModGamma for ice
        // [# m^-3 m^-1]
        dNdD[i] = IWCtopnd_MGD_IWC( diameter_volume_equ[i],rho,IWC );
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_MGD_LWC (//WS Output:
                 Vector& dNdD,
                 //WS Input:
                 const Vector& diameter_volume_equ,
                 const Numeric& rho,
                 const Numeric& LWC,
                 const Verbosity&)
{
    Index n_se = diameter_volume_equ.nelem();
    dNdD.resize(n_se);
    
    for ( Index i=0; i<n_se; i++ )
    {
        // calculate particle size distribution with ModGamma for liquid
        // [# m^-3 m^-1]
        dNdD[i] = LWCtopnd_MGD_LWC( diameter_volume_equ[i],rho ,LWC );
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_MH97 (//WS Output:
              Vector& dNdD,
              //WS Input:
              const Vector& diameter_mass_equivalent,
              const Numeric& IWC,
              const Numeric& T,
              const Index& noisy,
              const Index& robust,
              const Verbosity&)
{
  Index n_se = diameter_mass_equivalent.nelem();
  dNdD.resize(n_se);

  // abort if IWC is negative
  if ( !robust && IWC<0. )
    {
      ostringstream os;
      os << "Ice water content can not be negative."
         << " Yours is " << IWC << "kg/m3.";
      throw runtime_error ( os.str() );
    }
  Numeric iwc = max( IWC, 0. );

  // abort if T is negative
  if ( T<0. )
    {
      ostringstream os;
      os << "Negative temperatures not allowed.\n"
         << "Yours is " << T << "K.";
      throw runtime_error ( os.str() );
    }
  // abort if T is too high
  if ( !robust && T>280. )
    {
      ostringstream os;
      os << "Temperatures above 280K not allowed by MH97"
         << " (to allow: run with robust option).\n"
         << "Yours is " << T << "K.";
      throw runtime_error ( os.str() );
    }
  // allow some margin on T>0C (but use T=0C for PSD calc)
  Numeric t = min( T, 273.15 );

  // calculate particle size distribution with MH97
  // [# m^-3 m^-1]
  psd_cloudice_MH97 ( dNdD, diameter_mass_equivalent, iwc, t, noisy );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_MP48 (//WS Output:
              Vector& dNdD,
              //WS Input:
              const Vector& diameter_melted_equivalent,
              const Numeric& PR,
              const String& PRunit,
              const Numeric& rho,
              const Verbosity&)
{
  Numeric tPR;
  if (PRunit == "mm/h")
    tPR = PR;
  else if ( (PRunit == "SI") || (PRunit == "kg/m2/s") )
    {
      if (rho<=0.)
        {
          ostringstream os;
          os << "Precipitation unit " << PRunit
             << " requires valid material density (rho>0).\n"
             << "Yours is rho=" << rho << "kg/m3.\n";
          throw runtime_error ( os.str() );
        }
      tPR = PR * (3.6e6/rho);
    }
  else
    {
      ostringstream os;
      os << "Precipitation unit '" << PRunit << "' unknown.\n";
      throw runtime_error ( os.str() );
    }

  Index n_se = diameter_melted_equivalent.nelem();
  dNdD.resize(n_se);

  // derive particle number density for all given sizes
  for ( Index i=0; i<n_se; i++ )
    {
      // calculate particle size distribution with MP48
      // output: [# m^-3 m^-1]
      dNdD[i] = PRtopnd_MP48 ( tPR, diameter_melted_equivalent[i]);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_S2M (//WS Output:
                 Vector& dNdD,
                 //WS Input:
               const Vector& mass,
               const Numeric& N_tot,
               const Numeric& M,
               const String& psd_type,
               const Verbosity&)
{
    Matrix dummy;
    String hydrometeor_type;
    
    //Get the right hydrometeor type
    if (psd_type=="S2M_LWC")
    {
        hydrometeor_type="cloud_water";
    }
    else if (psd_type=="S2M_IWC")
    {
        hydrometeor_type="cloud_ice";
    }
    else if (psd_type=="S2M_RWC")
    {
        hydrometeor_type="rain";
    }
    else if (psd_type=="S2M_SWC")
    {
        hydrometeor_type="snow";
    }
    else if (psd_type=="S2M_GWC")
    {
        hydrometeor_type="graupel";
    }
    else if (psd_type=="S2M_HWC")
    {
        hydrometeor_type="hail";
    }
    else
    {
        ostringstream os;
        os << "You use a wrong tag! ";
        throw runtime_error( os.str() );
    }
    
    Numeric M1 = max( M, 0. );
    
    // calculate particle size distribution with S2M
    // [# m^-3 kg^-1]
    psd_S2M(dNdD, dummy,mass, N_tot, M1, hydrometeor_type) ;
    
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_S2M_M (//WS Output:
               Vector& dNdD,
               //WS Input:
               const Vector& mass,
               const Numeric& mean_mass,
               const Numeric& M,
               const String& psd_type,
               const Verbosity&)
{
    Numeric N_tot;
    Matrix dummy;
    String hydrometeor_type;
    
    // Calculate the total number density from mass density M and the
    // mean particle mass
    N_tot=M/mean_mass;
    
    //Get the right hydrometeor type
    if (psd_type=="S2M_LWC")
    {
        hydrometeor_type="cloud_water";
    }
    else if (psd_type=="S2M_IWC")
    {
        hydrometeor_type="cloud_ice";
    }
    else if (psd_type=="S2M_RWC")
    {
        hydrometeor_type="rain";
    }
    else if (psd_type=="S2M_SWC")
    {
        hydrometeor_type="snow";
    }
    else if (psd_type=="S2M_GWC")
    {
        hydrometeor_type="graupel";
    }
    else if (psd_type=="S2M_HWC")
    {
        hydrometeor_type="hail";
    }
    else
    {
        ostringstream os;
        os << "You use a wrong tag! ";
        throw runtime_error( os.str() );
    }
    
    Numeric M1 = max( M, 0. );
    
    // calculate particle size distribution with S2M
    // [# m^-3 kg^-1]
    psd_S2M(dNdD, dummy, mass, N_tot, M1, psd_type) ;

}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_MY2 (//WS Output:
               Vector& dNdD,
               //WS Input:
               const Vector& diameter_max,
               const Numeric& N_tot,
               const Numeric& M,
               const String& psd_type,
               const Verbosity&)
{
    Matrix dummy;
    String hydrometeor_type;
    
    //Get the right hydrometeor type
    if (psd_type=="MY2_LWC")
    {
        hydrometeor_type="cloud_water";
    }
    else if (psd_type=="MY2_IWC")
    {
        hydrometeor_type="cloud_ice";
    }
    else if (psd_type=="MY2_RWC")
    {
        hydrometeor_type="rain";
    }
    else if (psd_type=="MY2_SWC")
    {
        hydrometeor_type="snow";
    }
    else if (psd_type=="MY2_GWC")
    {
        hydrometeor_type="graupel";
    }
    else if (psd_type=="MY2_HWC")
    {
        hydrometeor_type="hail";
    }
    else
    {
        ostringstream os;
        os << "You use a wrong tag! ";
        throw runtime_error( os.str() );
    }
    
    Numeric M1 = max( M, 0. );
    
    // calculate particle size distribution with S2M
    // [# m^-3 kg^-1]
    psd_MY2(dNdD, dummy,diameter_max, N_tot, M1, hydrometeor_type);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_MY2_M (//WS Output:
                 Vector& dNdD,
                 //WS Input:
                 const Vector& diameter_max,
                 const Numeric& mean_mass,
                 const Numeric& M,
                 const String& psd_type,
                 const Verbosity&)
{
    Numeric N_tot;
    Matrix dummy;
    String hydrometeor_type;
    
    // Calculate the total number density from mass density M and the
    // mean particle mass
    N_tot=M/mean_mass;
    
    //Get the right hydrometeor type
    if (psd_type=="MY2_LWC")
    {
        hydrometeor_type="cloud_water";
    }
    else if (psd_type=="MY2_IWC")
    {
        hydrometeor_type="cloud_ice";
    }
    else if (psd_type=="MY2_RWC")
    {
        hydrometeor_type="rain";
    }
    else if (psd_type=="MY2_SWC")
    {
        hydrometeor_type="snow";
    }
    else if (psd_type=="MY2_GWC")
    {
        hydrometeor_type="graupel";
    }
    else if (psd_type=="MY2_HWC")
    {
        hydrometeor_type="hail";
    }
    else
    {
        ostringstream os;
        os << "You use a wrong tag! ";
        throw runtime_error( os.str() );
    }
    
    Numeric M1 = max( M, 0. );
    
    // calculate particle size distribution with S2M
    // [# m^-3 kg^-1]
    psd_MY2(dNdD, dummy,diameter_max, N_tot, M1, hydrometeor_type);

}



/* Workspace method: Doxygen documentation will be auto-generated */
void dNdD_W16 (//WS Output:
              Vector& dNdD,
              //WS Input:
              const Vector& diameter_mass_equivalent,
              const Numeric& RWC,
              const Index& robust,
              const Verbosity&)
{
  Index n_se = diameter_mass_equivalent.nelem();
  dNdD.resize(n_se);

  // abort if RWC is negative
  if ( !robust && RWC<0. )
    {
      ostringstream os;
      os << "Rain water content can not be negative."
         << " Yours is " << RWC << "kg/m3.";
      throw runtime_error ( os.str() );
    }
  Numeric rwc = max( RWC, 0. );

  // calculate particle size distribution with W16
  // [# m^-3 m^-1]
  psd_rain_W16 ( dNdD, diameter_mass_equivalent, rwc );
}

