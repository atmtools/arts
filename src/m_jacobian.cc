/* Copyright (C) 2004-2008 Mattias Ekstrom <ekstrom@rss.chalmers.se>

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



/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   m_jacobian.cc
  \author Mattias Ekstrom <ekstrom@rss.chalmers.se>
  \date   2004-09-14

  \brief  Workspace functions related to the jacobian.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <string>
#include "absorption.h"
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "math_funcs.h"
#include "messages.h"
#include "interpolation_poly.h"
#include "jacobian.h"
#include "physics_funcs.h"
#include "rte.h"

extern const String ABSSPECIES_MAINTAG;
extern const String FREQUENCY_MAINTAG;
extern const String FREQUENCY_SUBTAG_A;
extern const String FREQUENCY_CALCMODE_A;
extern const String POINTING_MAINTAG;
extern const String POINTING_SUBTAG_A;
extern const String POINTING_CALCMODE_A;
extern const String POLYFIT_MAINTAG;
extern const String TEMPERATURE_MAINTAG;



/*===========================================================================
  === The methods, with general methods first followed by the Add/Calc method
  === pairs for each retrieval quantity.
  ===========================================================================*/


//----------------------------------------------------------------------------
// General methods:
//----------------------------------------------------------------------------


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianClose(
        Index&                     jacobian_do,
        Matrix&                    jacobian,
        ArrayOfArrayOfIndex&       jacobian_indices,
  const ArrayOfRetrievalQuantity&  jacobian_quantities,
  const Matrix&                    sensor_pos,
  const Sparse&                    sensor_response )
{
  // Make sure that the array is not empty
  if( jacobian_quantities.nelem() == 0 )
    throw runtime_error(
          "No retrieval quantities has been added to *jacobian_quantities*." );

  // Check that sensor_pol and sensor_response has been initialised
  if( sensor_pos.nrows() == 0 )
    {
      ostringstream os;
      os << "The number of rows in *sensor_pos* is zero, i.e. no measurement\n"
         << "blocks has been defined. This has to be done before calling\n"
         << "jacobianClose.";
      throw runtime_error(os.str());
    }
  if( sensor_response.nrows() == 0 )
    {
      ostringstream os;
      os << "The sensor has either to be defined or turned off before calling\n"
         << "jacobianClose.";
      throw runtime_error(os.str());
    }

  // Loop over retrieval quantities, set JacobianIndices
  Index nrows = sensor_pos.nrows() * sensor_response.nrows();
  Index ncols = 0;
  //
  for( Index it=0; it<jacobian_quantities.nelem(); it++ )
    {
      // Store start jacobian index
      ArrayOfIndex indices(2);
      indices[0] = ncols;

      // Count total number of field points, i.e. product of grid lengths
      Index cols = 1;
      ArrayOfVector grids = jacobian_quantities[it].Grids();
      for( Index jt=0; jt<grids.nelem(); jt++ )
        { cols *= grids[jt].nelem(); }

      // Store stop index
      indices[1] = ncols + cols - 1;
      jacobian_indices.push_back( indices );

      ncols += cols;
    }
  
  // Resize *jacobian*
  jacobian.resize( nrows, ncols );
  jacobian = 0;

  jacobian_do = 1;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianInit(
        ArrayOfRetrievalQuantity&  jacobian_quantities,
        Agenda&                    jacobian_agenda )
{
  jacobian_quantities.resize(0);
  jacobian_agenda = Agenda();
  jacobian_agenda.set_name( "jacobian_agenda" );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianOff(
        Index&                     jacobian_do,
        Agenda&                    jacobian_agenda,
        ArrayOfRetrievalQuantity&  jacobian_quantities,
        ArrayOfArrayOfIndex&       jacobian_indices )
{
  jacobianInit( jacobian_quantities, jacobian_agenda );
  jacobian_do   = 0;
  jacobian_indices.resize(0);
}





//----------------------------------------------------------------------------
// Absorption species:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddAbsSpecies(
  Workspace&                  ws,
  ArrayOfRetrievalQuantity&   jq,
  Agenda&                     jacobian_agenda,
  const Index&                atmosphere_dim,
  const Vector&               p_grid,
  const Vector&               lat_grid,
  const Vector&               lon_grid,
  const Vector&               rq_p_grid,
  const Vector&               rq_lat_grid,
  const Vector&               rq_lon_grid,
  const String&               species,
  const String&               method,
  const String&               mode,
  const Numeric&              dx )
{
  // Check retrieval grids, here we just check the length of the grids
  // vs. the atmosphere dimension
  ArrayOfVector grids(atmosphere_dim);
  {
    ostringstream os;
    if( !check_retrieval_grids( grids, os, p_grid, lat_grid, lon_grid,
                                rq_p_grid, rq_lat_grid, rq_lon_grid,
                                "retrieval pressure grid", 
                                "retrieval latitude grid", 
                                "retrievallongitude_grid", 
                                atmosphere_dim ) )
    throw runtime_error(os.str());
  }
  
  // Check that method is either "analytical" or "perturbation"
  bool analytical;
  if( method == "perturbation" )
    { analytical = 0; }
  else if( method == "analytical" )
    { analytical = 1; }
  else
    {
      ostringstream os;
      os << "The method for absorption species retrieval can only be "
         << "\"analytical\"\n or \"perturbation\".";
      throw runtime_error(os.str());
    }
  
  // Check that mode is either "vmr", "nd" or "rel" with or without prefix log
  if( mode != "vmr" && mode != "nd" && mode != "rel" && mode != "logrel" )
    {
      throw runtime_error( "The retrieval mode can only be \"vmr\", \"nd\" "
                                                    "\"rel\" or \"logrel\"." );
    }

  // Check that this species is not already included in the jacobian.
  for( Index it=0; it<jq.nelem(); it++ )
    {
      if( jq[it].MainTag() == ABSSPECIES_MAINTAG  && jq[it].Subtag() == species )
        {
          ostringstream os;
          os << "The gas species:\n" << species << "\nis already included in "
             << "*jacobian_quantities*.";
          throw runtime_error(os.str());
        }
    }

  // Create the new retrieval quantity
  RetrievalQuantity rq = RetrievalQuantity();
  rq.MainTag( ABSSPECIES_MAINTAG );
  rq.Subtag( species );
  rq.Mode( mode );
  rq.Analytical( analytical );
  rq.Perturbation( dx );
  rq.Grids( grids );

  // Add it to the *jacobian_quantities*
  jq.push_back( rq );
  
  // Add gas species method to the jacobian agenda, if perturbation
  if( analytical )
    {
      out3 << "  Calculations done by semi-analytical expressions.\n"; 
    }
  else
    {
      out2 << "  Adding absorption species: " << species 
           << " to *jacobian_quantities*\n" << "  and *jacobian_agenda*\n";
      out3 << "  Calculations done by perturbation, size " << dx 
           << " " << mode << ".\n"; 

      jacobian_agenda.append( ws, "jacobianCalcAbsSpecies", species );
    }
}                    



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcAbsSpecies(
        Workspace&                  ws,
        Matrix&                     jacobian,
  const Index&                      imblock,
  const Vector&                     iyb _U_,
  const Vector&                     yb,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Tensor3&                    t_field,
  const Tensor4&                    vmr_field,
  const ArrayOfArrayOfSpeciesTag&   abs_species,
  const Index&                      cloudbox_on,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const Matrix&                     sensor_pos,
  const Matrix&                     sensor_los,
  const Vector&                     mblock_za_grid,
  const Vector&                     mblock_aa_grid,
  const Index&                      antenna_dim,
  const Sparse&                     sensor_response,
  const Agenda&                     iy_clearsky_agenda,
  const String&                     y_unit,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices,
  const String&                     species )
{
  // Set some useful variables. 
  RetrievalQuantity rq;
  ArrayOfIndex      ji;
  Index             it, pertmode;

  // Find the retrieval quantity related to this method, i.e. Abs. species -
  // species. This works since the combined MainTag and Subtag is individual.
  bool found = false;
  for( Index n=0; n<jacobian_quantities.nelem() && !found; n++ )
    {
      if( jacobian_quantities[n].MainTag() == ABSSPECIES_MAINTAG  &&  
          jacobian_quantities[n].Subtag()  == species )
        {
          found = true;
          rq    = jacobian_quantities[n];
          ji    = jacobian_indices[n];
        }
    }
  if( !found )
    {
      ostringstream os;
      os << "There is no gas species retrieval quantities defined for:\n"
         << species;
      throw runtime_error(os.str());
    }

  if( rq.Analytical() )
    {
      ostringstream os;
      os << "This WSM handles only perturbation calculations.\n"
         << "Are you using the method manually?";
      throw runtime_error(os.str());
    }
  
  // Store the start JacobianIndices and the Grids for this quantity
  it = ji[0];
  ArrayOfVector jg = rq.Grids();

  // Check if a relative pertubation is used or not, this information is needed
  // by the methods 'perturbation_field_?d'.
  // Note: both 'vmr' and 'nd' are absolute perturbations
  if( rq.Mode()=="rel" )
    pertmode = 0;
  else 
    pertmode = 1;

  // For each atmospheric dimension option calculate a ArrayOfGridPos, these
  // are the base functions for interpolating the perturbations into the
  // atmospheric grids.
  ArrayOfGridPos p_gp, lat_gp, lon_gp;
  Index j_p   = jg[0].nelem();
  Index j_lat = 1;
  Index j_lon = 1;
  //
  get_perturbation_gridpos( p_gp, p_grid, jg[0], true );
  //
  if( atmosphere_dim >= 2 ) 
    {
      j_lat = jg[1].nelem();
      get_perturbation_gridpos( lat_gp, lat_grid, jg[1], false );
      if( atmosphere_dim == 3 ) 
        {
          j_lon = jg[2].nelem();
          get_perturbation_gridpos( lon_gp, lon_grid, jg[2], false );
        }
    }

  // Find VMR field for this species. 
  ArrayOfSpeciesTag tags;
  array_species_tag_from_string( tags, species );
  Index si = chk_contains( "species", abs_species, tags );

  // Variables for vmr field perturbation unit conversion
  Tensor3 nd_field(0,0,0);
  if( rq.Mode()=="nd" )
    {
      nd_field.resize( t_field.npages(), t_field.nrows(), t_field.ncols() );
      calc_nd_field( nd_field, p_grid, t_field );
    }


  // Loop through the retrieval grid and calculate perturbation effect
  //
  const Index    n1y = sensor_response.nrows();
        Vector   dy( n1y ); 
  const Range    rowind = get_rowindex_for_mblock( sensor_response, imblock ); 
  //
  for( Index lon_it=0; lon_it<j_lon; lon_it++ )
    {
      for( Index lat_it=0; lat_it<j_lat; lat_it++ )
        {
          for (Index p_it=0; p_it<j_p; p_it++)
            {
              // Here we calculate the ranges of the perturbation. We want the
              // perturbation to continue outside the atmospheric grids for the
              // edge values.
              Range p_range   = Range(0,0);
              Range lat_range = Range(0,0);
              Range lon_range = Range(0,0);

              get_perturbation_range( p_range, p_it, j_p );

              if( atmosphere_dim>=2 )
                {
                  get_perturbation_range( lat_range, lat_it, j_lat );
                  if( atmosphere_dim == 3 )
                    {
                      get_perturbation_range( lon_range, lon_it, j_lon );
                    }
                }

              // Create VMR field to perturb
              Tensor4 vmr_p = vmr_field;
                              
              // If perturbation given in ND convert the vmr-field to ND before
              // the perturbation is added          
              if( rq.Mode() == "nd" )
                vmr_p(si,joker,joker,joker) *= nd_field;
        
              // Calculate the perturbed field according to atmosphere_dim, 
              // the number of perturbations is the length of the retrieval 
              // grid +2 (for the end points)
              switch (atmosphere_dim)
                {
                case 1:
                  {
                    // Here we perturb a vector
                    perturbation_field_1d( vmr_p(si,joker,lat_it,lon_it), 
                                           p_gp, jg[0].nelem()+2, p_range, 
                                           rq.Perturbation(), pertmode );
                    break;
                  }
                case 2:
                  {
                    // Here we perturb a matrix
                    perturbation_field_2d( vmr_p(si,joker,joker,lon_it),
                                           p_gp, lat_gp, jg[0].nelem()+2, 
                                           jg[1].nelem()+2, p_range, lat_range, 
                                           rq.Perturbation(), pertmode );
                    break;
                  }    
                case 3:
                  {  
                    // Here we need to perturb a tensor3
                    perturbation_field_3d( vmr_p(si,joker,joker,joker), 
                                           p_gp, lat_gp, lon_gp, 
                                           jg[0].nelem()+2,
                                           jg[1].nelem()+2, jg[2].nelem()+2, 
                                           p_range, lat_range, lon_range, 
                                           rq.Perturbation(), pertmode );
                    break;
                  }
                }

              // If perturbation given in ND convert back to VMR          
              if (rq.Mode()=="nd")
                vmr_p(si,joker,joker,joker) /= nd_field;
        
              // Calculate the perturbed spectrum  
              //
              Index         dummy2;
              Vector        iybp, dummy1;
              Matrix        dummy3;
              ArrayOfMatrix dummy4;      
              //
              iyb_calc( ws, iybp, dummy1, dummy2, dummy3, dummy4, imblock, 
                        atmosphere_dim, p_grid, lat_grid, lon_grid, 
                        t_field, vmr_p, cloudbox_on, stokes_dim, 
                        f_grid, sensor_pos, sensor_los, mblock_za_grid, 
                        mblock_aa_grid, antenna_dim, iy_clearsky_agenda, 0, 
                        y_unit, 0, ArrayOfRetrievalQuantity(), 
                        ArrayOfArrayOfIndex() );
              //
              mult( dy, sensor_response, iybp );

              // Difference spectrum
              for( Index i=0; i<n1y; i++ )
                { dy[i] = ( dy[i]- yb[i] ) / rq.Perturbation(); }

              // Put into jacobian
              jacobian(rowind,it) = dy;     

              // Result from next loop shall go into next column of J
              it++;
            }
        }
    }
}





//----------------------------------------------------------------------------
// Frequency:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddFreqShiftAndStretch(
  Workspace&                 ws,
  ArrayOfRetrievalQuantity&  jacobian_quantities,
  Agenda&                    jacobian_agenda,
  const String&              calcmode,
  const Numeric&             df,
  const Index&               do_stretch )
{
  // Check that do_stretch is 0 or 1
  if( do_stretch!=0 && do_stretch!=1 )
    throw runtime_error( "The argument *do_stretch* must be 0 or 1." );

  // Check that this type of frequency fit is not already included in the
  // jacobian.
  for( Index it=0; it<jacobian_quantities.nelem(); it++ )
    {
      if (jacobian_quantities[it].MainTag()== FREQUENCY_MAINTAG  &&  
          jacobian_quantities[it].Subtag() == FREQUENCY_SUBTAG_A )
        {
          ostringstream os;
          os << "This type of frequency fit is already included in\n"
             << "*jacobian_quantities*.";
          throw runtime_error(os.str());
        }
    }

  // Create the new retrieval quantity
  RetrievalQuantity rq = RetrievalQuantity();
  rq.MainTag( FREQUENCY_MAINTAG );
  rq.Subtag( FREQUENCY_SUBTAG_A );
  rq.Analytical( 0 );
  rq.Perturbation( df );
  // Shift and stretch are treated as order 0 and 1
  Vector grid(0,1+do_stretch,1);
  ArrayOfVector grids(1,grid);
  rq.Grids(grids);

  // Add pointing method to the jacobian agenda
  if( calcmode == "iybinterp" )
    { 
      rq.Mode( FREQUENCY_CALCMODE_A );
      jacobian_agenda.append( ws, "jacobianCalcFreqShiftAndStretchIybinterp", 
                                                                           "" );
   }
  else
    throw runtime_error( "Possible choices for *calcmode* are \"iybinterp\"." ); 

  // Add it to the *jacobian_quantities*
  jacobian_quantities.push_back( rq );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcFreqShiftAndStretchIybinterp(
        Matrix&                    jacobian,
  const Index&                     imblock,
  const Vector&                    iyb,
  const Vector&                    yb,
  const Index&                     stokes_dim,
  const Vector&                    f_grid,
  const Vector&                    mblock_za_grid,
  const Vector&                    mblock_aa_grid,
  const Index&                     antenna_dim,
  const Sparse&                    sensor_response,
  const ArrayOfIndex&              sensor_response_pol_grid,
  const Vector&                    sensor_response_f_grid,
  const Vector&                    sensor_response_za_grid,
  const ArrayOfRetrievalQuantity&  jacobian_quantities,
  const ArrayOfArrayOfIndex&       jacobian_indices )
{
  // Set some useful (and needed) variables.  
  RetrievalQuantity rq;
  ArrayOfIndex ji;

  // Find the retrieval quantity related to this method, i.e. Pointing
  // za offset. This works since the combined MainTag and Subtag is individual.
  bool found = false;
  for( Index n=0; n<jacobian_quantities.nelem() && !found; n++ )
    {
      if( jacobian_quantities[n].MainTag() == FREQUENCY_MAINTAG   && 
          jacobian_quantities[n].Subtag()  == FREQUENCY_SUBTAG_A  &&
          jacobian_quantities[n].Mode()    == FREQUENCY_CALCMODE_A )
        {
          found = true;
          rq = jacobian_quantities[n];
          ji = jacobian_indices[n];
        }
    }
  if( !found )
    {
      throw runtime_error(
                   "There is no such frequency retrieval quantity defined.\n" );
    }

  // Check that sensor_response is consistent with yb and iyb
  //
  if( sensor_response.nrows() != yb.nelem() )
    throw runtime_error( 
                       "Mismatch in size between *sensor_response* and *yb*." );
  if( sensor_response.ncols() != iyb.nelem() )
    throw runtime_error( 
                      "Mismatch in size between *sensor_response* and *iyb*." );


  // Size and another check of sensor_response
  //
  const Index nf     = sensor_response_f_grid.nelem();
  const Index npol   = sensor_response_pol_grid.nelem();
  const Index nza    = sensor_response_za_grid.nelem();


  // Get disturbed (part of) y
  //
  const Index    n1y = sensor_response.nrows(); 
        Vector   dy( n1y );
  {
    const Index   nf2      = f_grid.nelem();
    const Index   nza2     = mblock_za_grid.nelem();
          Index   naa2     = mblock_aa_grid.nelem();   
    if( antenna_dim == 1 )  
      { naa2 = 1; }
    const Index   niyb    = nf2 * nza2 * naa2 * stokes_dim;

    // Interpolation weights
    //
    const Index   porder = 3;
    //
    ArrayOfGridPosPoly   gp( nf2 );
                Matrix   itw( nf2, porder+1) ;
                Vector   fg_new = f_grid, iyb2(niyb);
    //
    fg_new += rq.Perturbation();
    gridpos_poly( gp, f_grid, fg_new, porder );
    interpweights( itw, gp );

    // Do interpolation
    for( Index iza=0; iza<nza2; iza++ )
      {
        for( Index iaa=0; iaa<naa2; iaa++ )
          {
            const Index row0 =( iza*naa2 + iaa ) * nf2 * stokes_dim;
            
            for( Index is=0; is<stokes_dim; is++ )
              { 
                interp( iyb2[Range(row0+is,nf2,stokes_dim)], itw, 
                         iyb[Range(row0+is,nf2,stokes_dim)], gp );
              }
          }
      }

    // Determine differnce
    //
    mult( dy, sensor_response, iyb2 );
    //
    for( Index i=0; i<n1y; i++ )
      { dy[i] = ( dy[i]- yb[i] ) / rq.Perturbation(); }
  }

  //--- Create jacobians ---
  //
  const Range rowind = get_rowindex_for_mblock( sensor_response, imblock );
  const Index lg = rq.Grids()[0].nelem();
        Index it = ji[0];

  // Shift jacobian
  jacobian(rowind,it) = dy;

  // Stretch jacobian
  if( lg > 1 )
    {
      const Index  row0 = rowind.get_start();
            Vector w;
           
      assert( Numeric(1) == rq.Grids()[0][1] );
      //
      polynomial_basis_func( w, sensor_response_f_grid, 1 );
      //
      it += 1;
      //
      for( Index l=0; l<nza; l++ )
        {
          for( Index f=0; f<nf; f++ )
            {
              const Index row1 = (l*nf + f)*npol;
              for( Index p=0; p<npol; p++ )
                { jacobian(row0+row1+p,it) = w[f] * dy[row1+p]; }
            }
        }
    }
}



//----------------------------------------------------------------------------
// Pointing:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddPointingZa(
  Workspace&                 ws,
  ArrayOfRetrievalQuantity&  jacobian_quantities,
  Agenda&                    jacobian_agenda,
  const Matrix&              sensor_pos,
  const Vector&              sensor_time,
  const Index&               poly_order,
  const String&              calcmode,
  const Numeric&             dza )
{
  // Check that poly_order is -1 or positive
  if( poly_order < -1 )
    throw runtime_error(
                  "The polynomial order has to be positive or -1 for gitter." );
 
  // Check that this jacobian type is not already included.
  for( Index it=0; it<jacobian_quantities.nelem(); it++ )
    {
      if (jacobian_quantities[it].MainTag()== POINTING_MAINTAG  &&  
          jacobian_quantities[it].Subtag() == POINTING_SUBTAG_A )
        {
          ostringstream os;
          os << "Fit of zenith angle pointing off-set is already included in\n"
             << "*jacobian_quantities*.";
          throw runtime_error(os.str());
        }
    }

  // Check that sensor_time is consistent with sensor_pos
  if( sensor_time.nelem() != sensor_pos.nrows() )
    {
      ostringstream os;
      os << "The WSV *sensor_time* must be defined for every "
         << "measurement block.\n";
      throw runtime_error(os.str());
    }

  // Do not allow that *poly_order* is not too large compared to *sensor_time*
  if( poly_order > sensor_time.nelem()-1 )
    { throw runtime_error( 
             "The polynomial order can not be >= length of *sensor_time*." ); }

  // Create the new retrieval quantity
  RetrievalQuantity rq = RetrievalQuantity();
  rq.MainTag( POINTING_MAINTAG );
  rq.Subtag( POINTING_SUBTAG_A );
  rq.Analytical( 0 );
  rq.Perturbation( dza );


  // To store the value or the polynomial order, create a vector with length
  // poly_order+1, in case of gitter set the size of the grid vector to be the
  // number of measurement blocks, all elements set to -1.
  Vector grid(0,poly_order+1,1);
  if( poly_order == -1 )
    {
      grid.resize(sensor_pos.nrows());
      grid = -1.0;
    }
  ArrayOfVector grids(1,grid);
  rq.Grids(grids);

  if( calcmode == "iybrecalc" )
    { 
      rq.Mode( POINTING_CALCMODE_A );  
      jacobian_agenda.append( ws, "jacobianCalcPointingZaIybrecalc", "" );
   }
  else
    throw runtime_error( "Possible choices for *calcmode* are \"iybrecalc\"." ); 

  // Add it to the *jacobian_quantities*
  jacobian_quantities.push_back( rq );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcPointingZaIybrecalc(
        Workspace&                  ws,
        Matrix&                     jacobian,
  const Index&                      imblock,
  const Vector&                     iyb _U_,
  const Vector&                     yb,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Tensor3&                    t_field,
  const Tensor4&                    vmr_field,
  const Index&                      cloudbox_on,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const Matrix&                     sensor_pos,
  const Matrix&                     sensor_los,
  const Vector&                     mblock_za_grid,
  const Vector&                     mblock_aa_grid,
  const Index&                      antenna_dim,
  const Sparse&                     sensor_response,
  const Vector&                     sensor_time,
  const Agenda&                     iy_clearsky_agenda,
  const String&                     y_unit,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices )
{
  // Set some useful variables.  
  RetrievalQuantity rq;
  ArrayOfIndex ji;

  // Find the retrieval quantity related to this method.
  // This works since the combined MainTag and Subtag is individual.
  bool found = false;
  for( Index n=0; n<jacobian_quantities.nelem() && !found; n++ )
    {
      if( jacobian_quantities[n].MainTag() == POINTING_MAINTAG    && 
          jacobian_quantities[n].Subtag()  == POINTING_SUBTAG_A   &&
          jacobian_quantities[n].Mode()    == POINTING_CALCMODE_A )
        {
          found = true;
          rq = jacobian_quantities[n];
          ji = jacobian_indices[n];
        }
    }
  if( !found )
    { throw runtime_error(
                "There is no such pointing retrieval quantity defined.\n" );
    }


  // Get disturbed (part of) y
  //
  const Index    n1y = sensor_response.nrows();
        Vector   dy( n1y );
  {
        Index         iyet;
        Vector        iyb2, iye;
        Matrix        iyb_aux, los = sensor_los;
        ArrayOfMatrix diyb_dx;      

    los(joker,0) += rq.Perturbation();

    iyb_calc( ws, iyb2, iye, iyet, iyb_aux, diyb_dx, imblock, 
              atmosphere_dim, p_grid, lat_grid, lon_grid, 
              t_field, vmr_field, cloudbox_on, stokes_dim, 
              f_grid, sensor_pos, los, mblock_za_grid, mblock_aa_grid, 
              antenna_dim, iy_clearsky_agenda, 0, y_unit, 
              0, ArrayOfRetrievalQuantity(), ArrayOfArrayOfIndex() );

    mult( dy, sensor_response, iyb2 );

    for( Index i=0; i<n1y; i++ )
      { dy[i] = ( dy[i]- yb[i] ) / rq.Perturbation(); }
  }

  //--- Create jacobians ---

  const Index lg = rq.Grids()[0].nelem();
  const Index it = ji[0];
  const Range rowind = get_rowindex_for_mblock( sensor_response, imblock );
  const Index row0 = rowind.get_start();

  // Handle gitter seperately
  if( rq.Grids()[0][0] == -1 )                  // Not all values are set here,
    {                                           // but should already have been 
      assert( lg == sensor_los.nrows() );       // set to 0
      assert( rq.Grids()[0][imblock] == -1 );
      jacobian(rowind,it+imblock) = dy;     
    }                                

  // Polynomial representation
  else
    {
      Vector w;
      for( Index c=0; c<lg; c++ )
        {
          assert( Numeric(c) == rq.Grids()[0][c] );
          //
          polynomial_basis_func( w, sensor_time, c );
          //
          for( Index i=0; i<n1y; i++ )
            { jacobian(row0+i,it+c) = w[imblock] * dy[i]; }
        }
    }
}




//----------------------------------------------------------------------------
// Polynomial baseline fits:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddPolyfit(
  Workspace&                 ws,
  ArrayOfRetrievalQuantity&  jq,
  Agenda&                    jacobian_agenda,
  const ArrayOfIndex&        sensor_response_pol_grid,
  const Vector&              sensor_response_f_grid,
  const Vector&              sensor_response_za_grid,
  const Matrix&              sensor_pos,
  const Index&               poly_order,
  const Index&               no_pol_variation,
  const Index&               no_los_variation,
  const Index&               no_mblock_variation )
{
  // Check that poly_order is >= 0
  if( poly_order < 0 )
    throw runtime_error(
      "The polynomial order has to be >= 0.");

  // Compare poly_order with number of frequency points
  if( poly_order > sensor_response_f_grid.nelem() - 1 )
    {
      ostringstream os;
      os << "The polynomial order can not exceed the length of\n"
         << "*sensor_response_f_grid* -1.";
      throw runtime_error(os.str());
    }
 
  // Check that polyfit is not already included in the jacobian.
  for( Index it=0; it<jq.nelem(); it++ )
    {
      if( jq[it].MainTag() == POLYFIT_MAINTAG )
        {
          ostringstream os;
          os << "Polynomial baseline fit is already included in\n"
             << "*jacobian_quantities*.";
          throw runtime_error(os.str());
        }
    }

  // "Grids"
  //
  // Grid dimensions correspond here to 
  //   1: polynomial order
  //   2: polarisation
  //   3: viewing direction
  //   4: measurement block
  //
  ArrayOfVector grids(4);
  //
  if( no_pol_variation )
    grids[1] = Vector(1,1);
  else
    grids[1] = Vector(0,sensor_response_pol_grid.nelem(),1);
  if( no_los_variation )
    grids[2] = Vector(1,1);
  else
    grids[2] = Vector(0,sensor_response_za_grid.nelem(),1); 
  if( no_mblock_variation )
    grids[3] = Vector(1,1);
  else
    grids[3] = Vector(0,sensor_pos.nrows(),1);

  // Create the new retrieval quantity
  RetrievalQuantity rq = RetrievalQuantity();
  rq.MainTag( POLYFIT_MAINTAG );
  rq.Mode( "" );
  rq.Analytical( 0 );
  rq.Perturbation( 0 );

  // Each polynomial coeff. is treated as a retrieval quantity
  //
  for( Index i=0; i<=poly_order; i++ )
    {
      ostringstream sstr;
      sstr << "Coefficient " << i;
      rq.Subtag( sstr.str() ); 

      // Grid is a scalar, use polynomial coeff.
      grids[0] = Vector(1,(Numeric)i);
      rq.Grids( grids );

      // Add it to the *jacobian_quantities*
      jq.push_back( rq );

      // Add pointing method to the jacobian agenda
      jacobian_agenda.append( ws, "jacobianCalcPolyfit", i );
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcPolyfit(
        Matrix&                   jacobian,
  const Index&                    imblock,
  const Vector&                   iyb _U_,
  const Vector&                   yb _U_,
  const Sparse&                   sensor_response,
  const ArrayOfIndex&             sensor_response_pol_grid,
  const Vector&                   sensor_response_f_grid,
  const Vector&                   sensor_response_za_grid,
  const ArrayOfRetrievalQuantity& jacobian_quantities,
  const ArrayOfArrayOfIndex&      jacobian_indices,
  const Index&                    poly_coeff )
{  
  // Find the retrieval quantity related to this method
  RetrievalQuantity rq;
  ArrayOfIndex ji;
  bool found = false;
  Index iq;
  ostringstream sstr;
  sstr << "Coefficient " << poly_coeff;
  for( iq=0; iq<jacobian_quantities.nelem() && !found; iq++ )
    {
      if( jacobian_quantities[iq].MainTag() == POLYFIT_MAINTAG  && 
          jacobian_quantities[iq].Subtag() == sstr.str() )
        {
          found = true;
          break;
        }
    }
  if( !found )
    {
      throw runtime_error( "There is no such polynomial baseline fit retrieval " 
                           "quantity defined.\n");
    }

  // Size and check of sensor_response
  //
  const Index nf     = sensor_response_f_grid.nelem();
  const Index npol   = sensor_response_pol_grid.nelem();
  const Index nza    = sensor_response_za_grid.nelem();

  // Make a vector with values to distribute over *jacobian*
  //
  Vector w; 
  //
  polynomial_basis_func( w, sensor_response_f_grid, poly_coeff );
  
  // Fill J
  //
  ArrayOfVector jg   = jacobian_quantities[iq].Grids();
  const Index n1     = jg[1].nelem();
  const Index n2     = jg[2].nelem();
  const Index n3     = jg[3].nelem();
  const Range rowind = get_rowindex_for_mblock( sensor_response, imblock );
  const Index row4   = rowind.get_start();
        Index col4   = jacobian_indices[iq][0];

  if( n3 > 1 )
    { col4 += imblock*n2*n1; }
      
  for( Index l=0; l<nza; l++ )
    {
      const Index row3 = row4 + l*nf*npol;
      const Index col3 = col4 + l*n1;

      for( Index f=0; f<nf; f++ )
        {
          const Index row2 = row3 + f*npol;

          for( Index p=0; p<npol; p++ )
            {
              Index col1 = col3;
              if( n1 > 1 )
                { col1 += p; }

              jacobian(row2+p,col1) = w[f];
            }
        }
    }
}



//----------------------------------------------------------------------------
// Temperatures (atmospheric):
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddTemperature(
  Workspace&                  ws,
  ArrayOfRetrievalQuantity&   jq,
  Agenda&                     jacobian_agenda,
  const Index&                atmosphere_dim,
  const Vector&               p_grid,
  const Vector&               lat_grid,
  const Vector&               lon_grid,
  const Vector&               rq_p_grid,
  const Vector&               rq_lat_grid,
  const Vector&               rq_lon_grid,
  const String&               hse,
  const String&               method,
  const Numeric&              dx )
{
  // Check retrieval grids, here we just check the length of the grids
  // vs. the atmosphere dimension
  ArrayOfVector grids(atmosphere_dim);
  {
    ostringstream os;
    if( !check_retrieval_grids( grids, os, p_grid, lat_grid, lon_grid,
                                rq_p_grid, rq_lat_grid, rq_lon_grid,
                                "retrieval pressure grid", 
                                "retrieval latitude grid", 
                                "retrievallongitude_grid", 
                                atmosphere_dim ) )
    throw runtime_error(os.str());
  }
  
  // Check that method is either "analytic" or "perturbation"
  bool analytical;
  if( method == "perturbation" )
    { analytical = 0; }
  else if( method == "analytical" )
    { analytical = 1; }
  else
    {
      ostringstream os;
      os << "The method for atmospheric temperature retrieval can only be "
         << "\"analytical\"\n or \"perturbation\".";
      throw runtime_error(os.str());
    }

  // Set subtag 
  String subtag;
  if( hse == "on" )
    {
      subtag = "HSE on";
      ostringstream os;
      os << "The keyword for hydrostatic equilibrium can so far only be set to\n"
         << "\"off\"\n";
      throw runtime_error(os.str());
    }
  else if( hse == "off" )
    {
      subtag = "HSE off";
    }
  else
    {
      ostringstream os;
      os << "The keyword for hydrostatic equilibrium can only be set to\n"
         << "\"on\" or \"off\"\n";
      throw runtime_error(os.str());
    }

  // HSE can not be combined with analytical 
  if( method == "analytical"  &&  hse == "on" )
    {
      ostringstream os;
      os << "Hydrostatic equilibrium can only be included for perturbation\n"
         << "calculations.";
      throw runtime_error(os.str());
    }

  // Check that temperature is not already included in the jacobian.
  // We only check the main tag.
  for (Index it=0; it<jq.nelem(); it++)
    {
      if( jq[it].MainTag() == "Temperature" )
        {
          ostringstream os;
          os << "Temperature is already included in *jacobian_quantities*.";
          throw runtime_error(os.str());
        }
    }

  // Create the new retrieval quantity
  RetrievalQuantity rq = RetrievalQuantity();
  rq.MainTag( TEMPERATURE_MAINTAG );
  rq.Subtag( subtag );
  rq.Mode( "abs" );
  rq.Analytical( analytical );
  rq.Perturbation( dx );
  rq.Grids( grids );

  // Add it to the *jacobian_quantities*
  jq.push_back( rq );

  if( analytical ) 
    {
      out3 << "  Calculations done by semi-analytical expression.\n"; 
    }
  else
    { 
      out3 << "  Calculations done by perturbations, of size " << dx << ".\n"; 

      jacobian_agenda.append( ws, "jacobianCalcTemperature", "" );
    }
}                    



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcTemperature(
        Workspace&                  ws,
        Matrix&                     jacobian,
  const Index&                      imblock,
  const Vector&                     iyb _U_,
  const Vector&                     yb,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Tensor3&                    t_field,
  const Tensor4&                    vmr_field,
  const Index&                      cloudbox_on,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const Matrix&                     sensor_pos,
  const Matrix&                     sensor_los,
  const Vector&                     mblock_za_grid,
  const Vector&                     mblock_aa_grid,
  const Index&                      antenna_dim,
  const Sparse&                     sensor_response,
  const Agenda&                     iy_clearsky_agenda,
  const String&                     y_unit,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices )
{
  // Set some useful variables. 
  RetrievalQuantity rq;
  ArrayOfIndex      ji;
  Index             it;

  // Find the retrieval quantity related to this method, i.e. Temperature.
  // For temperature only the main tag is checked.
  bool found = false;
  for( Index n=0; n<jacobian_quantities.nelem() && !found; n++ )
    {
      if( jacobian_quantities[n].MainTag() == TEMPERATURE_MAINTAG )
        {
          found = true;
          rq = jacobian_quantities[n];
          ji = jacobian_indices[n];
        }
    }
  if( !found )
    {
      ostringstream os;
      os << "There is no temperature retrieval quantities defined.\n";
      throw runtime_error(os.str());
    }

  if( rq.Analytical() )
    {
      ostringstream os;
      os << "This WSM handles only perturbation calculations.\n"
         << "Are you using the method manually?";
      throw runtime_error(os.str());
    }
  
  // FIXME: Only HSE off is implemented
  if( rq.Subtag() == "HSE on" )
      throw runtime_error( "HSE is so far not handled. Sorry!" );

  // Store the start JacobianIndices and the Grids for this quantity
  it = ji[0];
  ArrayOfVector jg = rq.Grids();

  // "Perturbation mode". 1 means absolute perturbations
  const Index pertmode = 1;
   
  // For each atmospheric dimension option calculate a ArrayOfGridPos, 
  // these will be used to interpolate a perturbation into the atmospheric 
  // grids.
  ArrayOfGridPos p_gp, lat_gp, lon_gp;
  Index j_p   = jg[0].nelem();
  Index j_lat = 1;
  Index j_lon = 1;
  //
  get_perturbation_gridpos( p_gp, p_grid, jg[0], true );
  //
  if( atmosphere_dim >= 2 ) 
    {
      j_lat = jg[1].nelem();
      get_perturbation_gridpos( lat_gp, lat_grid, jg[1], false );
      if( atmosphere_dim == 3 ) 
        {
          j_lon = jg[2].nelem();
          get_perturbation_gridpos( lon_gp, lon_grid, jg[2], false );
        }
    }


  // Loop through the retrieval grid and calculate perturbation effect
  //
  const Index    n1y = sensor_response.nrows();
        Vector   dy( n1y ); 
  const Range    rowind = get_rowindex_for_mblock( sensor_response, imblock ); 
  //
  for( Index lon_it=0; lon_it<j_lon; lon_it++ )
    {
      for( Index lat_it=0; lat_it<j_lat; lat_it++ )
        {
          for( Index p_it=0; p_it<j_p; p_it++ )
            {
              // Perturbed temperature field
              Tensor3 t_p = t_field;

              // Here we calculate the ranges of the perturbation. We want the
              // perturbation to continue outside the atmospheric grids for the
              // edge values.
              Range p_range   = Range(0,0);
              Range lat_range = Range(0,0);
              Range lon_range = Range(0,0);
              get_perturbation_range( p_range, p_it, j_p );
              if( atmosphere_dim >= 2 )
                {
                  get_perturbation_range( lat_range, lat_it, j_lat );
                  if( atmosphere_dim == 3 )
                    {
                      get_perturbation_range( lon_range, lon_it, j_lon );
                    }
                }
                           
              // Calculate the perturbed field according to atmosphere_dim, 
              // the number of perturbations is the length of the retrieval 
              // grid +2 (for the end points)
              switch (atmosphere_dim)
                {
                case 1:
                  {
                    // Here we perturb a vector
                    perturbation_field_1d( t_p(joker,lat_it,lon_it), 
                                           p_gp, jg[0].nelem()+2, p_range, 
                                           rq.Perturbation(), pertmode );
                    break;
                  }
                case 2:
                  {
                    // Here we perturb a matrix
                    perturbation_field_2d( t_p(joker,joker,lon_it), 
                                           p_gp, lat_gp, jg[0].nelem()+2, 
                                           jg[1].nelem()+2, p_range, lat_range, 
                                           rq.Perturbation(), pertmode );
                    break;
                  }    
                case 3:
                  {  
                    // Here we need to perturb a tensor3
                    perturbation_field_3d( t_p(joker,joker,joker), 
                                           p_gp, lat_gp, lon_gp, jg[0].nelem()+2,
                                           jg[1].nelem()+2, jg[2].nelem()+2, 
                                           p_range, lat_range, lon_range, 
                                           rq.Perturbation(), pertmode );
                    break;
                  }
                }
       
              // Calculate the perturbed spectrum  
              //
              Index         dummy2;
              Vector        iybp, dummy1;
              Matrix        dummy3;
              ArrayOfMatrix dummy4;      
              //
              iyb_calc( ws, iybp, dummy1, dummy2, dummy3, dummy4, imblock, 
                        atmosphere_dim, p_grid, lat_grid, lon_grid, 
                        t_p, vmr_field, cloudbox_on, stokes_dim, 
                        f_grid, sensor_pos, sensor_los, mblock_za_grid, 
                        mblock_aa_grid, antenna_dim, iy_clearsky_agenda, 0, 
                        y_unit, 0, ArrayOfRetrievalQuantity(), 
                        ArrayOfArrayOfIndex() );
              //
              mult( dy, sensor_response, iybp );

              // Difference spectrum
              for( Index i=0; i<n1y; i++ )
                { dy[i] = ( dy[i]- yb[i] ) / rq.Perturbation(); }

              // Put into jacobian
              jacobian(rowind,it) = dy;     

              // Result from next loop shall go into next column of J
              it++;
            }
        }
    }
}











//----------------------------------------------------------------------------
// Old:
//----------------------------------------------------------------------------


// /* Workspace method: Doxygen documentation will be auto-generated */
// void jacobianAddParticle(// WS Output:
//                          ArrayOfRetrievalQuantity& jq,
//                          Agenda&                   jacobian_agenda,
//                          // WS Input:
//                          const Matrix&             jac,
//                          const Index&              atmosphere_dim,
//                          const Vector&             p_grid,
//                          const Vector&             lat_grid,
//                          const Vector&             lon_grid,
//                          const Tensor4&            pnd_field,
//                          const Tensor5&            pnd_perturb,
//                          const ArrayOfIndex&       cloudbox_limits,
//                          // WS Generic Input:
//                          const Vector&             rq_p_grid,
//                          const Vector&             rq_lat_grid,
//                          const Vector&             rq_lon_grid)
// {
//   throw runtime_error("Particle jacobians not yet handled correctly.");

//   // Check that the jacobian matrix is empty. Otherwise it is either
//   // not initialised or it is closed.
//   if( jac.nrows()!=0 && jac.ncols()!=0 )
//     {
//       ostringstream os;
//       os << "The Jacobian matrix is not initialised correctly or closed.\n"
//          << "New retrieval quantities can not be added at this point.";
//       throw runtime_error(os.str());
//     }
  
//   // Check that pnd_perturb is consistent with pnd_field
//   if( pnd_perturb.nbooks()!=pnd_field.nbooks() ||
//       pnd_perturb.npages()!=pnd_field.npages() ||
//       pnd_perturb.nrows()!=pnd_field.nrows() ||
//       pnd_perturb.ncols()!=pnd_field.ncols() )
//     {
//       ostringstream os;
//       os << "The perturbation field *pnd_field_perturb* is not consistent with"
//          << "*pnd_field*,\none or several dimensions do not match.";
//       throw runtime_error(os.str());
//     }
  
//   // Check that particles are not already included in the jacobian.
//   for( Index it=0; it<jq.nelem(); it++ )
//     {
//       if( jq[it].MainTag()=="Particles" )
//         {
//           ostringstream os;
//           os << "The particles number densities are already included in "
//              << "*jacobian_quantities*.";
//           throw runtime_error(os.str());
//         }
//     }
  
//   // Particle Jacobian only defined for 1D and 3D atmosphere. check the 
//   // retrieval grids, here we just check the length of the grids vs. the 
//   // atmosphere dimension
//   if (atmosphere_dim==2) 
//   {
//     ostringstream os;
//     os << "Atmosphere dimension not equal to 1 or 3. " 
//        << "Jacobians for particle number\n"
//        << "density only available for 1D and 3D atmosphere.";
//     throw runtime_error(os.str());
//   }

//   ArrayOfVector grids(atmosphere_dim);
//   // The retrieval grids should only consists of gridpoints within
//   // the cloudbox. Setup local atmospheric fields inside the cloudbox
//   {
//     Vector p_cbox = p_grid;
//     Vector lat_cbox = lat_grid;
//     Vector lon_cbox = lon_grid;
//     switch (atmosphere_dim)
//       {
//       case 3:
//         {
//           lon_cbox = lon_grid[Range(cloudbox_limits[4], 
//                                     cloudbox_limits[5]-cloudbox_limits[4]+1)];
//         }
//       case 2:
//         {
//           lat_cbox = lat_grid[Range(cloudbox_limits[2], 
//                                     cloudbox_limits[3]-cloudbox_limits[2]+1)];
//         }    
//       case 1:
//         {  
//           p_cbox = p_grid[Range(cloudbox_limits[0], 
//                                 cloudbox_limits[1]-cloudbox_limits[0]+1)];
//         }
//       }
//     ostringstream os;
//     if( !check_retrieval_grids( grids, os, p_cbox, lat_cbox, lon_cbox,
//                                 rq_p_grid, rq_lat_grid, rq_lon_grid, 
//         // FIXMEOLE: These strings have to replaced later with the proper
//         //           names from the WSM documentation in methods.cc
//           "rq_p_grid", "rq_lat_grid", "rq_lon_grid", atmosphere_dim ))
//       throw runtime_error(os.str());
//   }

//   // Common part for all particle variables
//   RetrievalQuantity rq = RetrievalQuantity();
//   rq.MainTag("Particles");
//   rq.Grids(grids);
//   rq.Analytical(0);
//   rq.Perturbation(-999.999);
//   rq.Mode("Fields *mode* and *perturbation* are not defined");

//   // Set info for each particle variable
//   for( Index ipt=0; ipt<pnd_perturb.nshelves(); ipt++ )
//     {
//       out2 << "  Adding particle variable " << ipt +1 
//            << " to *jacobian_quantities / agenda*.\n";

//       ostringstream os;
//       os << "Variable " << ipt+1;
//       rq.Subtag(os.str());
      
//       jq.push_back(rq);
//     }

//   // Add gas species method to the jacobian agenda
//   String methodname = "jacobianCalcParticle";
//   String kwv = "";
//   jacobian_agenda.append (methodname, kwv);
// }                    












// /* Workspace method: Doxygen documentation will be auto-generated */
// void jacobianCalcParticle(
//            Workspace&                  ws,
//      // WS Output:
//            Matrix&                     jacobian,
//      // WS Input:
//      const Vector&                     y,
//      const ArrayOfRetrievalQuantity&   jq,
//      const ArrayOfArrayOfIndex&        jacobian_indices,
//      const Tensor5&                    pnd_field_perturb,
//      const Agenda&                     jacobian_particle_update_agenda,
//      const Agenda&                     ppath_step_agenda,
//      const Agenda&                     rte_agenda,
//      const Agenda&                     iy_space_agenda,
//      const Agenda&                     surface_prop_agenda,
//      const Agenda&                     iy_cloudbox_agenda,
//      const Index&                      atmosphere_dim,
//      const Vector&                     p_grid,
//      const Vector&                     lat_grid,
//      const Vector&                     lon_grid,
//      const Tensor3&                    z_field,
//      const Tensor3&                    t_field,
//      const Tensor4&                    vmr_field,
//      const Matrix&                     r_geoid,
//      const Matrix&                     z_surface,
//      const Index&                      cloudbox_on,
//      const ArrayOfIndex&               cloudbox_limits,
//      const Tensor4&                    pnd_field,
//      const Sparse&                     sensor_response,
//      const Matrix&                     sensor_pos,
//      const Matrix&                     sensor_los,
//      const Vector&                     f_grid,
//      const Index&                      stokes_dim,
//      const Index&                      antenna_dim,
//      const Vector&                     mblock_za_grid,
//      const Vector&                     mblock_aa_grid )
// {
//   // Set some useful (and needed) variables. 
//   Index n_jq = jq.nelem();
//   RetrievalQuantity rq;
//   ArrayOfIndex ji;
  
//   // Setup local atmospheric fields inside the cloudbox
//   Vector p_cbox = p_grid;
//   Vector lat_cbox = lat_grid;
//   Vector lon_cbox = lon_grid;
//   switch (atmosphere_dim)
//     {
//     case 3:
//       {
//         lon_cbox = lon_grid[Range(cloudbox_limits[4], 
//                                   cloudbox_limits[5]-cloudbox_limits[4]+1)];
//       }
//     case 2:
//       {
//         lat_cbox = lat_grid[Range(cloudbox_limits[2], 
//                                   cloudbox_limits[3]-cloudbox_limits[2]+1)];
//       }    
//     case 1:
//       {  
//         p_cbox = p_grid[Range(cloudbox_limits[0], 
//                               cloudbox_limits[1]-cloudbox_limits[0]+1)];
//       }
//     }


//   // Variables to handle and store perturbations
//   Vector yp;
//   Tensor4  pnd_p, base_pert = pnd_field;


//   // Loop particle variables (indexed by *ipt*, where *ipt* is zero based)
//   //
//   Index ipt      = -1;
//   bool not_ready = true;
  
//   while( not_ready )
//     {
//       // Step *ipt*
//       ipt++;

//       // Define sub-tag string
//       ostringstream os;
//       os << "Variable " << ipt+1;

//       // Find the retrieval quantity related to this particle type
//       //
//       bool  found = false;
//       //
//       for( Index n=0; n<n_jq; n++ )
//         {
//           if( jq[n].MainTag()=="Particles" && jq[n].Subtag()== os.str() )
//             {
//               found = true;
//               rq = jq[n];
//               ji = jacobian_indices[n];
//               n  = n_jq;                   // To jump out of for-loop
//             }
//         }

//       // At least one particle type must be found
//       assert( !( ipt==0  &&  !found ) );

//       // Ready or something to do?
//       if( !found )
//         { 
//           not_ready = false;
//         }
//       else
//         {
//           // Counters for report string
//           Index   it  = 0;
//           Index   nit = ji[1] -ji[0] + 1;
          
//           // Counter for column in *jacobian*
//           Index   icol = ji[0];

//           // Retrieval grid positions
//           ArrayOfVector jg = rq.Grids();
//           ArrayOfGridPos p_gp, lat_gp, lon_gp;
//           Index j_p = jg[0].nelem();
//           Index j_lat = 1;
//           Index j_lon = 1;
//           get_perturbation_gridpos( p_gp, p_cbox, jg[0], true );
//           if (atmosphere_dim==3) 
//             {
//               j_lat = jg[1].nelem();
//               get_perturbation_gridpos( lat_gp, lat_cbox, jg[1], false );
              
//               j_lon = jg[2].nelem();
//               get_perturbation_gridpos( lon_gp, lon_cbox, jg[2], false );
//             }

//           // Give verbose output
//           out1 << "  Calculating retrieval quantity:" << rq << "\n";
  

//           // Loop through the retrieval grid and calculate perturbation effect
//           for (Index lon_it=0; lon_it<j_lon; lon_it++)
//             {
//               for (Index lat_it=0; lat_it<j_lat; lat_it++)
//                 {
//                   for (Index p_it=0; p_it<j_p; p_it++)
//                     {
//                       // Update the perturbation field
//                       pnd_p = 
//                            pnd_field_perturb( ipt, joker, joker, joker, joker);

//                       it++;
//                       out1 << "  Calculating perturbed spectra no. " << it
//                            << " of " << nit << "\n";

//                       // Here we calculate the ranges of the perturbation. 
//                       // We want the perturbation to continue outside the 
//                       // atmospheric grids for the edge values.
//                       Range p_range   = Range(0,0);
//                       Range lat_range = Range(0,0);
//                       Range lon_range = Range(0,0);
//                       get_perturbation_range( p_range, p_it, j_p );
//                       if (atmosphere_dim==3)
//                         {
//                           get_perturbation_range( lat_range, lat_it, j_lat);
//                           get_perturbation_range( lon_range, lon_it, j_lon);
//                         }
                          
//                       // Make empty copy of pnd_pert for base functions
//                       base_pert *= 0;
            
//                       // Calculate the perturbed field according to atm_dim, 
//                       switch (atmosphere_dim)
//                         {
//                         case 1:
//                           {
//                             for( Index typ_it=0; typ_it<pnd_field.nbooks(); 
//                                                                      typ_it++ )
//                               {
//                                 perturbation_field_1d( 
//                                       base_pert(typ_it,joker,lat_it,lon_it),
//                                       p_gp, jg[0].nelem()+2, p_range, 1.0, 1 );
//                               }
//                             break;
//                           }
//                         case 3:
//                           {  
//                             for( Index typ_it=0; typ_it<pnd_field.nrows(); 
//                                                                      typ_it++ )
//                               {
//                                 perturbation_field_3d( 
//                                       base_pert(typ_it,joker,joker,joker),
//                                       p_gp, lat_gp, lon_gp, jg[0].nelem()+2, 
//                                       jg[1].nelem()+2, jg[2].nelem()+2, 
//                                       p_range, lat_range, lon_range, 1.0, 1);
//                               }
//                             break;
//                           }
//                         }
          
//                       // Now add the weighted perturbation field to the 
//                       // reference field and recalculate the scattered field
//                       pnd_p *= base_pert;
//                       pnd_p += pnd_field;
//                       jacobian_particle_update_agendaExecute( ws, pnd_p, 
//                                       jacobian_particle_update_agenda );
            
//                       // Calculate the perturbed spectrum  
//                       yCalc( ws, yp, ppath_step_agenda, rte_agenda, 
//                              iy_space_agenda, surface_prop_agenda, 
//                              iy_cloudbox_agenda, atmosphere_dim,
//                              p_grid, lat_grid, lon_grid, z_field, t_field, 
//                              vmr_field, r_geoid, z_surface, cloudbox_on, 
//                              cloudbox_limits, sensor_response, sensor_pos, 
//                              sensor_los, f_grid, stokes_dim, antenna_dim, 
//                              mblock_za_grid, mblock_aa_grid);
    
//                       // Add dy as column in jacobian. Note that we just return
//                       // the difference between the two spectra.
//                       for( Index y_it=0; y_it<yp.nelem(); y_it++ )
//                         {
//                           jacobian(y_it,icol) = yp[y_it]-y[y_it];
//                         }

//                       // Step *icol*
//                       icol++;
//                     }
//                 }
//             }
//         }
//     }
// }


                     











