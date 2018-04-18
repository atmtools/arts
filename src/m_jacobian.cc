/* Copyright (C) 2004-2012 Mattias Ekstrom <ekstrom@rss.chalmers.se>

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
#include "cloudbox.h"
#include "math_funcs.h"
#include "messages.h"
#include "interpolation_poly.h"
#include "jacobian.h"
#include "physics_funcs.h"
#include "rte.h"
#include "m_xml.h"

extern const Numeric PI;

extern const String ABSSPECIES_MAINTAG;
extern const String FREQUENCY_MAINTAG;
extern const String FREQUENCY_SUBTAG_0;
extern const String FREQUENCY_SUBTAG_1;
extern const String POINTING_MAINTAG;
extern const String POINTING_SUBTAG_A;
extern const String POINTING_CALCMODE_A;
extern const String POINTING_CALCMODE_B;
extern const String POLYFIT_MAINTAG;
extern const String SCATSPECIES_MAINTAG;
extern const String SINEFIT_MAINTAG;
extern const String TEMPERATURE_MAINTAG;
extern const String NLTE_MAINTAG;
extern const String WIND_MAINTAG;
extern const String MAGFIELD_MAINTAG;
extern const String FLUX_MAINTAG;
extern const String PROPMAT_SUBSUBTAG;

extern const String SURFACE_MAINTAG;


// Generic modes
extern const String PRESSUREBROADENINGGAMMA_MODE;
extern const String LINESTRENGTH_MODE;
extern const String LINECENTER_MODE;
extern const String LINEMIXINGY_MODE;
extern const String LINEMIXINGG_MODE;
extern const String LINEMIXINGDF_MODE;

// Modes for "some" catalogs
//  Pressure Broadening
extern const String SELFBROADENING_MODE;
extern const String FOREIGNBROADENING_MODE;
extern const String WATERBROADENING_MODE;
extern const String SELFBROADENINGEXPONENT_MODE;
extern const String FOREIGNBROADENINGEXPONENT_MODE;
extern const String WATERBROADENINGEXPONENT_MODE;
//  Line Mixing
extern const String LINEMIXINGY0_MODE;
extern const String LINEMIXINGG0_MODE;
extern const String LINEMIXINGDF0_MODE;
extern const String LINEMIXINGY1_MODE;
extern const String LINEMIXINGG1_MODE;
extern const String LINEMIXINGDF1_MODE;
extern const String LINEMIXINGYEXPONENT_MODE;
extern const String LINEMIXINGGEXPONENT_MODE;
extern const String LINEMIXINGDFEXPONENT_MODE;



/*===========================================================================
  === The methods, with general methods first followed by the Add/Calc method
  === pairs for each retrieval quantity.
  ===========================================================================*/


//----------------------------------------------------------------------------
// Basic methods:
//----------------------------------------------------------------------------


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcDoNothing(
        Matrix&     jacobian _U_,
  const Index&      mblock_index _U_,
  const Vector&     iyb _U_,
  const Vector&     yb _U_,
  const Verbosity& )
{
  /* Nothing to do here for the analytical case, this function just exists
   to satisfy the required inputs and outputs of the jacobian_agenda */
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianClose(
        Workspace&                 ws,
        Index&                     jacobian_do,
        Agenda&                    jacobian_agenda,
  const ArrayOfRetrievalQuantity&  jacobian_quantities,
  const Matrix&                    sensor_pos,
  const Sparse&                    sensor_response,
  const Verbosity&                 verbosity )
{
  // Make sure that the array is not empty
  if( jacobian_quantities.empty() )
    throw runtime_error(
          "No retrieval quantities has been added to *jacobian_quantities*." );

  // Check that sensor_pos and sensor_response has been initialised
  if( sensor_pos.empty() )
    {
      ostringstream os;
      os << "*sensor_pos* is empty, i.e. no measurement blocks has been "
         << "defined.\nThis has to be done before calling jacobianClose.";
      throw runtime_error(os.str());
    }
  if( sensor_response.empty() )
    {
      ostringstream os;
      os << "The sensor has either to be defined or turned off before calling\n"
         << "jacobianClose.";
      throw runtime_error(os.str());
    }

  jacobian_agenda.check( ws, verbosity );
  
  jacobian_do = 1;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianInit(
        ArrayOfRetrievalQuantity&  jacobian_quantities,
        Agenda&                    jacobian_agenda,
  const Verbosity& )
{
  jacobian_quantities.resize(0);
  jacobian_agenda = Agenda();
  jacobian_agenda.set_name( "jacobian_agenda" );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianOff(
        Index&                     jacobian_do,
        Index&                     jacobianDoit_do,
        Agenda&                    jacobian_agenda,
        ArrayOfRetrievalQuantity&  jacobian_quantities, 
   const Verbosity&                 verbosity )
{
  jacobian_do = 0;
  jacobianDoit_do = 0;
  jacobianInit( jacobian_quantities,
                jacobian_agenda, verbosity );
}





//----------------------------------------------------------------------------
// Absorption species:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddAbsSpecies(
        Workspace&,
        ArrayOfRetrievalQuantity&   jq,
        Agenda&                     jacobian_agenda,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Vector&                     rq_p_grid,
  const Vector&                     rq_lat_grid,
  const Vector&                     rq_lon_grid,
  const String&                     species,
  const String&                     method,
  const String&                     mode,
  const Index&                      for_species_tag,
  const Numeric&                    dx,
  const Verbosity&                  verbosity )
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  QuantumIdentifier qi;
  if(not for_species_tag)
  {
      ArrayOfSpeciesTag test;
      array_species_tag_from_string(test,species);
      if( test.nelem() not_eq 1 )
          throw std::runtime_error("Trying to add a species as a species tag of multiple species.\n"
          "This is not supported.  Please give just a single species instead.\n"
          "Otherwise consider if you intended for_species_tag to be evaluated true.\n");
      qi.SetAll();
      qi.SetIsotopologue(test[0].Isotopologue());
      qi.SetSpecies(test[0].Species());
  }
  
  // Check that this species is not already included in the jacobian.
  for( Index it=0; it<jq.nelem(); it++ )
    {
        if( jq[it].MainTag() == ABSSPECIES_MAINTAG  && jq[it].SubSubtag() != PROPMAT_SUBSUBTAG &&
          jq[it].Subtag()  == species )
        {
          ostringstream os;
          os << "The gas species:\n" << species << "\nis already included in "
             << "*jacobian_quantities*.";
          throw runtime_error(os.str());
        }
        else if( jq[it].MainTag() == ABSSPECIES_MAINTAG  && jq[it].SubSubtag() == PROPMAT_SUBSUBTAG )
        {
            if(SpeciesTag(jq[it].Subtag()) == SpeciesTag(species))
            {
                ostringstream os;
                os << "The atmospheric species of:\n" << species << "\nis already included in "
                << "*jacobian_quantities*.";
                throw runtime_error(os.str());
            }
        }
    }

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
  Index analytical;
  if( method == "perturbation" )
    { analytical = 0; }
  else if( method == "analytical")
    { analytical = 1; }
  else
    {
      ostringstream os;
      os << "The method for absorption species retrieval can only be "
         << "\"analytical\"\n or \"perturbation\".";
      throw runtime_error(os.str());
    }
  
  // Check that mode is correct
  if( mode != "vmr" && mode != "nd" && mode != "rel" && mode != "rh" && mode != "q" )
    {
      throw runtime_error( "The retrieval mode can only be \"vmr\", \"nd\", "
                           "\"rel\", \"rh\" or \"q\"." );
    }
  if( ( mode == "rh" || mode == "q" ) && species.substr(0,3) != "H2O" )
    {
      throw runtime_error( "Retrieval modes \"rh\" and \"q\" can only be applied "
                           "on species starting with H2O." );
    }
  
  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag( ABSSPECIES_MAINTAG );
  rq.Subtag( species );
  rq.Mode( mode );
  rq.Analytical( analytical );
  rq.Perturbation( dx );
  rq.Grids( grids );
  if(analytical and not for_species_tag) {
    rq.SubSubtag(species);
    rq.PropType(JacPropMatType::VMR);
  }
  else if((not analytical) and (not for_species_tag))
    throw std::runtime_error("perturbation only support for_species_tag true/\n ");
  rq.QuantumIdentity(qi);
  
  // Add it to the *jacobian_quantities*
  jq.push_back( rq );
  
  // Add gas species method to the jacobian agenda
  if( analytical )
    {
      out3 << "  Calculations done by semi-analytical expressions.\n"; 
      jacobian_agenda.append( "jacobianCalcDoNothing", TokVal() );
    }
  else
    {
      out2 << "  Adding absorption species: " << species 
           << " to *jacobian_quantities*\n" << "  and *jacobian_agenda*\n";
      out3 << "  Calculations done by perturbation, size " << dx 
           << " " << mode << ".\n"; 

      jacobian_agenda.append( "jacobianCalcAbsSpeciesPerturbations", species );
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddConstantVMRAbsSpecies(
  Workspace&,
  ArrayOfRetrievalQuantity&   jq,
  Agenda&                     jacobian_agenda,
  const String&               species,
  const String&               mode,
  const Index&                for_species_tag,
  const Numeric&              dx,
  const Verbosity&            verbosity )
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  if(not for_species_tag)
  {
    ArrayOfSpeciesTag test;
    array_species_tag_from_string(test,species);
    if( test.nelem()!=1 )
      throw std::runtime_error("Trying to add a species as a species tag of multiple species.\n"
      "This is not supported.  Please give just a single species instead.\n"
      "Otherwise consider if you intended for_species_tag to be evaluated true.\n");
  }
  
  // Check that this species is not already included in the jacobian.
  for( Index it=0; it<jq.nelem(); it++ )
  {
    RetrievalQuantity& rt = jq[it];
    if( rt.MainTag() == ABSSPECIES_MAINTAG)
    {
      if(rt.Subtag() == species)
      {
        ostringstream os;
        os << "jacobian already set for " << species << "\nThis is not allowed.";
        throw std::runtime_error(os.str());
      }
      
      if(not for_species_tag)
      {
        if(SpeciesTag(rt.Subtag()).Species() == SpeciesTag(species).Species())
        {
          ostringstream os;
          os << "for_species_tag is set to indicate full in one jacobian for species=\"" << species << "\""
             << "\nThis is the same species as exist in another jacobian species=\"" <<  rt.Subtag() << "\""
             << "\nSince this duplicates calculations, it is not allowed.";
          throw std::runtime_error(os.str());
        }
      }
    }
  }
  
  // Check that mode is either "vmr", "nd" or "rel" 
  if( mode != "vmr" && mode != "rel" && mode != "logrel" )
  {
    throw runtime_error( "The retrieval mode can only be \"vmr\", "
    "\"rel\" or \"logrel\"." );
  }
  
  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag( ABSSPECIES_MAINTAG );
  rq.Subtag( species );
  rq.Mode( mode);
  rq.Analytical( 1 );
  rq.Perturbation( dx );
  rq.SubSubtag(PROPMAT_SUBSUBTAG);
  rq.PropType(JacPropMatType::VMR);
  rq.IntegrationOn();
  
  // Add it to the *jacobian_quantities*
  jq.push_back( rq );
  
  // Add dummy
  jacobian_agenda.append( "jacobianCalcDoNothing", TokVal() );  
}      



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcAbsSpeciesPerturbations(
        Workspace&                  ws,
        Matrix&                     jacobian,
  const Index&                      mblock_index,
  const Vector&                     iyb _U_,
  const Vector&                     yb,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Tensor3&                    t_field,
  const Tensor3&                    z_field,
  const Tensor4&                    vmr_field,
  const Tensor4&                    nlte_field, 
  const ArrayOfArrayOfSpeciesTag&   abs_species,
  const Index&                      cloudbox_on,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const Matrix&                     sensor_pos,
  const Matrix&                     sensor_los,
  const Matrix&                     transmitter_pos,
  const Matrix&                     mblock_dlos_grid,
  const Sparse&                     sensor_response,
  const String&                     iy_unit,  
  const Agenda&                     iy_main_agenda,
  const Agenda&                     geo_pos_agenda,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const String&                     species,
  const Verbosity&                  verbosity)
{
  // Some useful variables. 
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
          bool any_affine;
          ArrayOfArrayOfIndex jacobian_indices;
          jac_ranges_indices( jacobian_indices, any_affine,
                              jacobian_quantities, true );
          //
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

  // Check if a relative perturbation is used or not, this information is needed
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
  const Range    rowind = get_rowindex_for_mblock( sensor_response, mblock_index ); 
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
              Vector        iybp;
              ArrayOfVector dummy3;      
              ArrayOfMatrix dummy4;
              Matrix        dummy5;
              //
              iyb_calc( ws, iybp, dummy3, dummy4, dummy5, mblock_index, 
                        atmosphere_dim, t_field, z_field,
                        vmr_p, nlte_field, cloudbox_on, 
                        stokes_dim, f_grid, sensor_pos, sensor_los, 
                        transmitter_pos, mblock_dlos_grid, 
                        iy_unit, iy_main_agenda, geo_pos_agenda,
                        0, ArrayOfRetrievalQuantity(), 
                        ArrayOfArrayOfIndex(), ArrayOfString(), verbosity );
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
// Frequency shift
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddFreqShift(
        Workspace&                 ws _U_,
        ArrayOfRetrievalQuantity&  jacobian_quantities,
        Agenda&                    jacobian_agenda,
  const Vector&                    f_grid,
  const Matrix&                    sensor_pos,
  const Vector&                    sensor_time,
  const Index&                     poly_order,
  const Numeric&                   df,
  const Verbosity& )
{
  // Check that poly_order is -1 or positive
  if( poly_order < -1 )
    throw runtime_error(
                  "The polynomial order has to be positive or -1 for gitter." );
 
  // Check that this jacobian type is not already included.
  for( Index it=0; it<jacobian_quantities.nelem(); it++ )
    {
      if (jacobian_quantities[it].MainTag()== FREQUENCY_MAINTAG  &&  
          jacobian_quantities[it].Subtag() == FREQUENCY_SUBTAG_0 )
        {
          ostringstream os;
          os << "Fit of frequency shift is already included in\n"
             << "*jacobian_quantities*.";
          throw runtime_error(os.str());
        }
    }

  // Checks of frequencies
  if( df <= 0 )
    throw runtime_error( "The argument *df* must be > 0." );
  if( df > 1e6 )
    throw runtime_error( "The argument *df* is not allowed to exceed 1 MHz." );
  const Index nf    = f_grid.nelem();
  if( nf < 2 )
    throw runtime_error( "Frequency shifts and *f_grid* of length 1 can "
                         "not be combined."  );
  const Numeric maxdf = f_grid[nf-1] - f_grid[nf-2]; 
  if( df > maxdf )
    {
      ostringstream os;
      os << "The value of *df* is too big with respect to spacing of "
         << "*f_grid*. The maximum\nallowed value of *df* is the spacing "
         << "between the two last elements of *f_grid*.\n"
         << "This spacing is   : " <<maxdf/1e3 << " kHz\n"
         << "The value of df is: " << df/1e3   << " kHz";
      throw runtime_error(os.str());
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
  RetrievalQuantity rq;
  rq.MainTag( FREQUENCY_MAINTAG );
  rq.Subtag( FREQUENCY_SUBTAG_0 );
  rq.Mode( "" );
  rq.Analytical( 0 );
  rq.Perturbation( df );

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

  // Add it to the *jacobian_quantities*
  jacobian_quantities.push_back( rq );

  // Add corresponding calculation method to the jacobian agenda
  jacobian_agenda.append( "jacobianCalcFreqShift", "" );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcFreqShift(
        Matrix&                    jacobian,
  const Index&                     mblock_index,
  const Vector&                    iyb,
  const Vector&                    yb,
  const Index&                     stokes_dim,
  const Vector&                    f_grid,
  const Matrix&                    DEBUG_ONLY(sensor_los),
  const Matrix&                    mblock_dlos_grid,
  const Sparse&                    sensor_response,
  const Vector&                    sensor_time,
  const ArrayOfRetrievalQuantity&  jacobian_quantities,
  const Verbosity& )
{
  // Set some useful (and needed) variables.  
  RetrievalQuantity rq;
  ArrayOfIndex ji;

  // Find the retrieval quantity related to this method.
  // This works since the combined MainTag and Subtag is individual.
  bool found = false;
  for( Index n=0; n<jacobian_quantities.nelem() && !found; n++ )
    {
      if( jacobian_quantities[n].MainTag() == FREQUENCY_MAINTAG   && 
          jacobian_quantities[n].Subtag()  == FREQUENCY_SUBTAG_0 )
        {
          bool any_affine;
          ArrayOfArrayOfIndex jacobian_indices;
          jac_ranges_indices( jacobian_indices, any_affine,
                              jacobian_quantities, true );
          //
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

  // Get disturbed (part of) y
  //
  const Index    n1y = sensor_response.nrows(); 
        Vector   dy( n1y );
  {
    const Index   nf2   = f_grid.nelem();
    const Index   nlos2 = mblock_dlos_grid.nrows();
    const Index   niyb  = nf2 * nlos2 * stokes_dim;

    // Interpolation weights
    //
    const Index   porder = 3;
    //
    ArrayOfGridPosPoly   gp( nf2 );
                Matrix   itw( nf2, porder+1) ;
                Vector   fg_new = f_grid, iyb2(niyb);
    //
    fg_new += rq.Perturbation();
    gridpos_poly( gp, f_grid, fg_new, porder, 1.0 );
    interpweights( itw, gp );

    // Do interpolation
    for( Index ilos=0; ilos<nlos2; ilos++ )
      {
        const Index row0 = ilos * nf2 * stokes_dim;
            
        for( Index is=0; is<stokes_dim; is++ )
          { 
            interp( iyb2[Range(row0+is,nf2,stokes_dim)], itw, 
                    iyb[Range(row0+is,nf2,stokes_dim)], gp );
          }
      }

    // Determine difference
    //
    mult( dy, sensor_response, iyb2 );
    //
    for( Index i=0; i<n1y; i++ )
      { dy[i] = ( dy[i]- yb[i] ) / rq.Perturbation(); }
  }

 //--- Create jacobians ---

  const Index lg = rq.Grids()[0].nelem();
  const Index it = ji[0];
  const Range rowind = get_rowindex_for_mblock( sensor_response, mblock_index );
  const Index row0 = rowind.get_start();

  // Handle gitter seperately
  if( rq.Grids()[0][0] == -1 )                  // Not all values are set here,
    {                                           // but should already have been 
      assert( lg == sensor_los.nrows() );       // set to 0
      assert( rq.Grids()[0][mblock_index] == -1 );
      jacobian(rowind,it+mblock_index) = dy;     
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
            { jacobian(row0+i,it+c) = w[mblock_index] * dy[i]; }
        }
    }
}




//----------------------------------------------------------------------------
// Frequency stretch
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddFreqStretch(
        Workspace&                 ws _U_,
        ArrayOfRetrievalQuantity&  jacobian_quantities,
        Agenda&                    jacobian_agenda,
  const Vector&                    f_grid,
  const Matrix&                    sensor_pos,
  const Vector&                    sensor_time,
  const Index&                     poly_order,
  const Numeric&                   df,
  const Verbosity& )
{
  // Check that poly_order is -1 or positive
  if( poly_order < -1 )
    throw runtime_error(
                  "The polynomial order has to be positive or -1 for gitter." );
 
  // Check that this jacobian type is not already included.
  for( Index it=0; it<jacobian_quantities.nelem(); it++ )
    {
      if (jacobian_quantities[it].MainTag()== FREQUENCY_MAINTAG  &&  
          jacobian_quantities[it].Subtag() == FREQUENCY_SUBTAG_1 )
        {
          ostringstream os;
          os << "Fit of frequency stretch is already included in\n"
             << "*jacobian_quantities*.";
          throw runtime_error(os.str());
        }
    }

  // Checks of df
  if( df <= 0 )
    throw runtime_error( "The argument *df* must be > 0." );
  if( df > 1e6 )
    throw runtime_error( "The argument *df* is not allowed to exceed 1 MHz." );
  const Index   nf    = f_grid.nelem();
  const Numeric maxdf = f_grid[nf-1] - f_grid[nf-2]; 
  if( df > maxdf )
    {
      ostringstream os;
      os << "The value of *df* is too big with respect to spacing of "
         << "*f_grid*. The maximum\nallowed value of *df* is the spacing "
         << "between the two last elements of *f_grid*.\n"
         << "This spacing is   : " <<maxdf/1e3 << " kHz\n"
         << "The value of df is: " << df/1e3   << " kHz";
      throw runtime_error(os.str());
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
  RetrievalQuantity rq;
  rq.MainTag( FREQUENCY_MAINTAG );
  rq.Subtag( FREQUENCY_SUBTAG_1 );
  rq.Mode( "" );
  rq.Analytical( 0 );
  rq.Perturbation( df );

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

  // Add it to the *jacobian_quantities*
  jacobian_quantities.push_back( rq );

  // Add corresponding calculation method to the jacobian agenda
  jacobian_agenda.append( "jacobianCalcFreqStretch", "" );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcFreqStretch(
        Matrix&                    jacobian,
  const Index&                     mblock_index,
  const Vector&                    iyb,
  const Vector&                    yb,
  const Index&                     stokes_dim,
  const Vector&                    f_grid,
  const Matrix&                    DEBUG_ONLY(sensor_los),
  const Matrix&                    mblock_dlos_grid,
  const Sparse&                    sensor_response,
  const ArrayOfIndex&              sensor_response_pol_grid,
  const Vector&                    sensor_response_f_grid,
  const Matrix&                    sensor_response_dlos_grid,
  const Vector&                    sensor_time,
  const ArrayOfRetrievalQuantity&  jacobian_quantities,
  const Verbosity& )
{
  // The code here is close to identical to the one for Shift. The main
  // difference is that dy is weighted with poly_order 1 basis function.

  // Set some useful (and needed) variables.  
  RetrievalQuantity rq;
  ArrayOfIndex ji;

  // Find the retrieval quantity related to this method.
  // This works since the combined MainTag and Subtag is individual.
  bool found = false;
  for( Index n=0; n<jacobian_quantities.nelem() && !found; n++ )
    {
      if( jacobian_quantities[n].MainTag() == FREQUENCY_MAINTAG   && 
          jacobian_quantities[n].Subtag()  == FREQUENCY_SUBTAG_1 )
        {
          bool any_affine;
          ArrayOfArrayOfIndex jacobian_indices;
          jac_ranges_indices( jacobian_indices, any_affine,
                              jacobian_quantities, true );
          //
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

  // Get disturbed (part of) y
  //
  const Index    n1y = sensor_response.nrows(); 
        Vector   dy( n1y );
  {
    const Index   nf2   = f_grid.nelem();
    const Index   nlos2 = mblock_dlos_grid.nrows();
    const Index   niyb  = nf2 * nlos2 * stokes_dim;

    // Interpolation weights
    //
    const Index   porder = 3;
    //
    ArrayOfGridPosPoly   gp( nf2 );
                Matrix   itw( nf2, porder+1) ;
                Vector   fg_new = f_grid, iyb2(niyb);
    //
    fg_new += rq.Perturbation();
    gridpos_poly( gp, f_grid, fg_new, porder, 1.0 );
    interpweights( itw, gp );

    // Do interpolation
    for( Index ilos=0; ilos<nlos2; ilos++ )
      {
        const Index row0 = ilos * nf2 * stokes_dim;
            
        for( Index is=0; is<stokes_dim; is++ )
          { 
            interp( iyb2[Range(row0+is,nf2,stokes_dim)], itw, 
                    iyb[Range(row0+is,nf2,stokes_dim)], gp );
          }
      }

    // Determine difference
    //
    mult( dy, sensor_response, iyb2 );
    //
    for( Index i=0; i<n1y; i++ )
      { dy[i] = ( dy[i]- yb[i] ) / rq.Perturbation(); }

    // dy above corresponds now to shift. Convert to stretch:
    //
    Vector w;
    polynomial_basis_func( w, sensor_response_f_grid, 1 );
    //
    const Index nf     = sensor_response_f_grid.nelem();
    const Index npol   = sensor_response_pol_grid.nelem();
    const Index nlos   = sensor_response_dlos_grid.nrows();
    //
    for( Index l=0; l<nlos; l++ )
      {    
        for( Index f=0; f<nf; f++ )
          {
            const Index row1 = (l*nf + f)*npol;
            for( Index p=0; p<npol; p++ )
              { dy[row1+p] *= w[f]; }
          }
      }
  }

 //--- Create jacobians ---

  const Index lg = rq.Grids()[0].nelem();
  const Index it = ji[0];
  const Range rowind = get_rowindex_for_mblock( sensor_response, mblock_index );
  const Index row0 = rowind.get_start();

  // Handle gitter seperately
  if( rq.Grids()[0][0] == -1 )                  // Not all values are set here,
    {                                           // but should already have been 
      assert( lg == sensor_los.nrows() );       // set to 0
      assert( rq.Grids()[0][mblock_index] == -1 );
      jacobian(rowind,it+mblock_index) = dy;     
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
            { jacobian(row0+i,it+c) = w[mblock_index] * dy[i]; }
        }
    }
}





//----------------------------------------------------------------------------
// Pointing:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddPointingZa(
        Workspace&                 ws _U_,
        ArrayOfRetrievalQuantity&  jacobian_quantities,
        Agenda&                    jacobian_agenda,
  const Matrix&                    sensor_pos,
  const Vector&                    sensor_time,
  const Index&                     poly_order,
  const String&                    calcmode,
  const Numeric&                   dza,
  const Verbosity& )
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

  // Checks of dza
  if( dza <= 0 )
    throw runtime_error( "The argument *dza* must be > 0." );
  if( dza > 0.1 )
    throw runtime_error( 
                     "The argument *dza* is not allowed to exceed 0.1 deg." );

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
  RetrievalQuantity rq;
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

  if( calcmode == "recalc" )
    { 
      rq.Mode( POINTING_CALCMODE_A );  
      jacobian_agenda.append( "jacobianCalcPointingZaRecalc", "" );
   }
  else if( calcmode == "interp" )
    { 
      rq.Mode( POINTING_CALCMODE_B );  
      jacobian_agenda.append( "jacobianCalcPointingZaInterp", "" );
   }
  else
    throw runtime_error( 
            "Possible choices for *calcmode* are \"recalc\" and \"interp\"." );

  // Add it to the *jacobian_quantities*
  jacobian_quantities.push_back( rq );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcPointingZaInterp(
        Matrix&                    jacobian,
  const Index&                     mblock_index,
  const Vector&                    iyb,
  const Vector&                    yb _U_,
  const Index&                     stokes_dim,
  const Vector&                    f_grid,
  const Matrix&                    DEBUG_ONLY(sensor_los),
  const Matrix&                    mblock_dlos_grid,
  const Sparse&                    sensor_response,
  const Vector&                    sensor_time,
  const ArrayOfRetrievalQuantity&  jacobian_quantities,
  const Verbosity& )
{
  if( mblock_dlos_grid.nrows() < 2 )
    throw runtime_error( "The method demands that *mblock_dlos_grid* has "
                         "more than one row." );

  if( !( is_increasing( mblock_dlos_grid(joker,0) )  || 
         is_decreasing( mblock_dlos_grid(joker,0) ) ) )
    throw runtime_error( "The method demands that the zenith angles in "
             "*mblock_dlos_grid* are sorted (increasing or decreasing)." );

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
          jacobian_quantities[n].Mode()    == POINTING_CALCMODE_B )
        {
          bool any_affine;
          ArrayOfArrayOfIndex jacobian_indices;
          jac_ranges_indices( jacobian_indices, any_affine,
                              jacobian_quantities, true );
          //
          found = true;
          rq = jacobian_quantities[n];
          ji = jacobian_indices[n];
        }
    }
  if( !found )
    { throw runtime_error(
                "There is no such pointing retrieval quantity defined.\n" );
    }


  // Get "dy", by inter/extra-polation of existing iyb
  //
  const Index    n1y = sensor_response.nrows();
        Vector   dy( n1y );
  {
    // Sizes
    const Index   nf  = f_grid.nelem();
    const Index   nza = mblock_dlos_grid.nrows();

    // Shifted zenith angles
    Vector za1 = mblock_dlos_grid(joker,0); za1 -= rq.Perturbation();
    Vector za2 = mblock_dlos_grid(joker,0); za2 += rq.Perturbation();

    // Find interpolation weights
    ArrayOfGridPos gp1(nza), gp2(nza);
    gridpos( gp1, mblock_dlos_grid(joker,0), za1, 1e6 );  // Note huge extrapolation!
    gridpos( gp2, mblock_dlos_grid(joker,0), za2, 1e6 );  // Note huge extrapolation!
    Matrix itw1(nza,2), itw2(nza,2);
    interpweights( itw1, gp1 );
    interpweights( itw2, gp2 );

    // Make interpolation (for all azimuth angles, frequencies and Stokes)
    //
    Vector  iyb1(iyb.nelem()), iyb2(iyb.nelem());
    //
    for( Index iza=0; iza<nza; iza++ )
      {
        for( Index iv=0; iv<nf; iv++ )
          {
            for( Index is=0; is<stokes_dim; is++ )
              {
                const Range r( iv*stokes_dim+is, nza, nf*stokes_dim );
                interp( iyb1[r], itw1, iyb[r], gp1 );
                interp( iyb2[r], itw2, iyb[r], gp2 );
              }
          }
      }

    // Apply sensor and take difference
    //
    Vector y1(n1y), y2(n1y);
    mult( y1, sensor_response, iyb1 );
    mult( y2, sensor_response, iyb2 );
    //
    for( Index i=0; i<n1y; i++ )
      { dy[i] = ( y2[i]- y1[i] ) / ( 2* rq.Perturbation() ); }
  }

  //--- Create jacobians ---

  const Index lg = rq.Grids()[0].nelem();
  const Index it = ji[0];
  const Range rowind = get_rowindex_for_mblock( sensor_response, mblock_index );
  const Index row0 = rowind.get_start();

  // Handle pointing "jitter" seperately
  if( rq.Grids()[0][0] == -1 )                  // Not all values are set here,
    {                                           // but should already have been 
      assert( lg == sensor_los.nrows() );       // set to 0
      assert( rq.Grids()[0][mblock_index] == -1 );
      jacobian(rowind,it+mblock_index) = dy;     
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
            { jacobian(row0+i,it+c) = w[mblock_index] * dy[i]; }
        }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcPointingZaRecalc(
        Workspace&                 ws,
        Matrix&                    jacobian,
  const Index&                     mblock_index,
  const Vector&                    iyb _U_,
  const Vector&                    yb,
  const Index&                     atmosphere_dim,
  const Tensor3&                   t_field,
  const Tensor3&                   z_field,
  const Tensor4&                   vmr_field,
  const Tensor4&                   nlte_field, 
  const Index&                     cloudbox_on,
  const Index&                     stokes_dim,
  const Vector&                    f_grid,
  const Matrix&                    sensor_pos,
  const Matrix&                    sensor_los,
  const Matrix&                    transmitter_pos,
  const Matrix&                    mblock_dlos_grid,
  const Sparse&                    sensor_response,
  const Vector&                    sensor_time,
  const String&                    iy_unit,  
  const Agenda&                    iy_main_agenda,
  const Agenda&                    geo_pos_agenda,
  const ArrayOfRetrievalQuantity&  jacobian_quantities,
  const Verbosity&                 verbosity )
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
          bool any_affine;
          ArrayOfArrayOfIndex jacobian_indices;
          jac_ranges_indices( jacobian_indices, any_affine,
                              jacobian_quantities, true );
          //
          found = true;
          rq = jacobian_quantities[n];
          ji = jacobian_indices[n];
        }
    }
  if( !found )
    { throw runtime_error(
                "There is no such pointing retrieval quantity defined.\n" );
    }

  // Get "dy", by calling iyb_calc with shifted sensor_los.
  //
  const Index    n1y = sensor_response.nrows();
        Vector   dy( n1y );
  {
        Vector        iyb2;
        Matrix        los = sensor_los;
        Matrix        geo_pos;
        ArrayOfVector iyb_aux;      
        ArrayOfMatrix diyb_dx;      

    los(joker,0) += rq.Perturbation();

    iyb_calc( ws, iyb2, iyb_aux, diyb_dx, geo_pos, 
              mblock_index, atmosphere_dim, 
              t_field, z_field, vmr_field, nlte_field,  cloudbox_on, stokes_dim, 
              f_grid, sensor_pos, los, transmitter_pos, mblock_dlos_grid, 
              iy_unit, iy_main_agenda, geo_pos_agenda,
              0, ArrayOfRetrievalQuantity(), ArrayOfArrayOfIndex(),
              ArrayOfString(), verbosity );

    // Apply sensor and take difference
    //
    mult( dy, sensor_response, iyb2 );
    //
    for( Index i=0; i<n1y; i++ )
      { dy[i] = ( dy[i]- yb[i] ) / rq.Perturbation(); }
  }

  //--- Create jacobians ---

  const Index lg = rq.Grids()[0].nelem();
  const Index it = ji[0];
  const Range rowind = get_rowindex_for_mblock( sensor_response, mblock_index );
  const Index row0 = rowind.get_start();

  // Handle "jitter" seperately
  if( rq.Grids()[0][0] == -1 )                  // Not all values are set here,
    {                                           // but should already have been 
      assert( lg == sensor_los.nrows() );       // set to 0
      assert( rq.Grids()[0][mblock_index] == -1 );
      jacobian(rowind,it+mblock_index) = dy;     
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
            { jacobian(row0+i,it+c) = w[mblock_index] * dy[i]; }
        }
    }
}





//----------------------------------------------------------------------------
// Polynomial baseline fits:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddPolyfit(
        Workspace&                 ws _U_,
        ArrayOfRetrievalQuantity&  jq,
        Agenda&                    jacobian_agenda,
  const ArrayOfIndex&              sensor_response_pol_grid,
  const Matrix&                    sensor_response_dlos_grid,
  const Matrix&                    sensor_pos,
  const Index&                     poly_order,
  const Index&                     no_pol_variation,
  const Index&                     no_los_variation,
  const Index&                     no_mblock_variation,
  const Verbosity& )
{
  // Check that poly_order is >= 0
  if( poly_order < 0 )
    throw runtime_error( "The polynomial order has to be >= 0.");

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
    grids[2] = Vector(0,sensor_response_dlos_grid.nrows(),1); 
  if( no_mblock_variation )
    grids[3] = Vector(1,1);
  else
    grids[3] = Vector(0,sensor_pos.nrows(),1);

  // Create the new retrieval quantity
  RetrievalQuantity rq;
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
      jacobian_agenda.append( "jacobianCalcPolyfit", i );
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcPolyfit(
        Matrix&                    jacobian,
  const Index&                     mblock_index,
  const Vector&                    iyb _U_,
  const Vector&                    yb _U_,
  const Sparse&                    sensor_response,
  const ArrayOfIndex&              sensor_response_pol_grid,
  const Vector&                    sensor_response_f_grid,
  const Matrix&                    sensor_response_dlos_grid,
  const ArrayOfRetrievalQuantity&  jacobian_quantities,
  const Index&                     poly_coeff,
  const Verbosity& )
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
      throw runtime_error( "There is no Polyfit jacobian defined, in general " 
                           "or for the selected polynomial coefficient.\n");
    }

  // Size and check of sensor_response
  //
  const Index nf     = sensor_response_f_grid.nelem();
  const Index npol   = sensor_response_pol_grid.nelem();
  const Index nlos   = sensor_response_dlos_grid.nrows();

  // Make a vector with values to distribute over *jacobian*
  //
  Vector w; 
  //
  polynomial_basis_func( w, sensor_response_f_grid, poly_coeff );
  
  // Fill J
  //
  ArrayOfArrayOfIndex jacobian_indices;
  {
    bool any_affine;
    jac_ranges_indices( jacobian_indices, any_affine,
                        jacobian_quantities, true );
  }
  //
  ArrayOfVector jg   = jacobian_quantities[iq].Grids();
  const Index n1     = jg[1].nelem();
  const Index n2     = jg[2].nelem();
  const Index n3     = jg[3].nelem();
  const Range rowind = get_rowindex_for_mblock( sensor_response, mblock_index );
  const Index row4   = rowind.get_start();
        Index col4   = jacobian_indices[iq][0];

  if( n3 > 1 )
    { col4 += mblock_index*n2*n1; }
      
  for( Index l=0; l<nlos; l++ )
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
// Scattering species:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddScatSpecies(
        Workspace&,
        ArrayOfRetrievalQuantity&   jq,
        Agenda&                     jacobian_agenda,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Vector&                     rq_p_grid,
  const Vector&                     rq_lat_grid,
  const Vector&                     rq_lon_grid,
  const String&                     species,
  const String&                     quantity,
  const Verbosity&                  verbosity )
{
  CREATE_OUT2;
  CREATE_OUT3;

  // Check that this species+quantity combination is not already included in
  // the jacobian.
  for( Index it=0; it<jq.nelem(); it++ )
    {
      if( jq[it].MainTag()   == SCATSPECIES_MAINTAG  &&
          jq[it].Subtag()    == species              &&
          jq[it].SubSubtag() == quantity  )
        {
          ostringstream os;
          os << "The combintaion of\n   scattering species: " << species
             << "\n   retrieval quantity: " <<quantity 
             << "\nis already included in *jacobian_quantities*.";
          throw runtime_error(os.str());
        }
    }

  
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
  

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag( SCATSPECIES_MAINTAG );
  rq.Subtag( species );
  rq.SubSubtag( quantity );
  rq.Analytical( 1 );
  rq.Grids( grids );
  
  // Add it to the *jacobian_quantities*
  jq.push_back( rq );
  
  // Add gas species method to the jacobian agenda
  jacobian_agenda.append( "jacobianCalcDoNothing", TokVal() );
}                    



//----------------------------------------------------------------------------
// Sinusoidal baseline fits:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddSinefit(
        Workspace&                 ws _U_,
        ArrayOfRetrievalQuantity&  jq,
        Agenda&                    jacobian_agenda,
  const ArrayOfIndex&              sensor_response_pol_grid,
  const Matrix&                    sensor_response_dlos_grid,
  const Matrix&                    sensor_pos,
  const Vector&                    period_lengths,
  const Index&                     no_pol_variation,
  const Index&                     no_los_variation,
  const Index&                     no_mblock_variation,
  const Verbosity& )
{
  const Index np = period_lengths.nelem();

  // Check that poly_order is >= 0
  if( np == 0 )
    throw runtime_error( "No sinusoidal periods has benn given.");

  // Check that polyfit is not already included in the jacobian.
  for( Index it=0; it<jq.nelem(); it++ )
    {
      if( jq[it].MainTag() == SINEFIT_MAINTAG )
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
    grids[2] = Vector(0,sensor_response_dlos_grid.nrows(),1); 
  if( no_mblock_variation )
    grids[3] = Vector(1,1);
  else
    grids[3] = Vector(0,sensor_pos.nrows(),1);

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag( SINEFIT_MAINTAG );
  rq.Mode( "" );
  rq.Analytical( 0 );
  rq.Perturbation( 0 );

  // Each sinefit coeff. pair is treated as a retrieval quantity
  //
  for( Index i=0; i<np; i++ )
    {
      ostringstream sstr;
      sstr << "Period " << i;
      rq.Subtag( sstr.str() ); 

      // "Grid" has length 2, set to period length
      grids[0] = Vector( 2, period_lengths[i] );
      rq.Grids( grids );

      // Add it to the *jacobian_quantities*
      jq.push_back( rq );

      // Add pointing method to the jacobian agenda
      jacobian_agenda.append( "jacobianCalcSinefit", i );
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcSinefit(
        Matrix&                    jacobian,
  const Index&                     mblock_index,
  const Vector&                    iyb _U_,
  const Vector&                    yb _U_,
  const Sparse&                    sensor_response,
  const ArrayOfIndex&              sensor_response_pol_grid,
  const Vector&                    sensor_response_f_grid,
  const Matrix&                    sensor_response_dlos_grid,
  const ArrayOfRetrievalQuantity&  jacobian_quantities,
  const Index&                     period_index,
  const Verbosity& )
{  
  // Find the retrieval quantity related to this method
  RetrievalQuantity rq;
  ArrayOfIndex ji;
  bool found = false;
  Index iq;
  ostringstream sstr;
  sstr << "Period " << period_index;
  for( iq=0; iq<jacobian_quantities.nelem() && !found; iq++ )
    {
      if( jacobian_quantities[iq].MainTag() == SINEFIT_MAINTAG  && 
          jacobian_quantities[iq].Subtag() == sstr.str() )
        {
          found = true;
          break;
        }
    }
  if( !found )
    {
      throw runtime_error( "There is no Sinefit jacobian defined, in general " 
                           "or for the selected period length.\n");
    }

  // Size and check of sensor_response
  //
  const Index nf     = sensor_response_f_grid.nelem();
  const Index npol   = sensor_response_pol_grid.nelem();
  const Index nlos   = sensor_response_dlos_grid.nrows();

  // Make vectors with values to distribute over *jacobian*
  //
  // (period length stored in grid 0)
  ArrayOfVector jg   = jacobian_quantities[iq].Grids();
  //
  Vector s(nf), c(nf); 
  //
  for( Index f=0; f<nf; f++ )
    {
      Numeric a = (sensor_response_f_grid[f]-sensor_response_f_grid[0]) * 
                                                             2 * PI / jg[0][0];
      s[f] = sin( a );
      c[f] = cos( a );
    }

  
  // Fill J
  //
  ArrayOfArrayOfIndex jacobian_indices;
  {
    bool any_affine;
    jac_ranges_indices( jacobian_indices, any_affine,
                        jacobian_quantities, true );
  }
  //  
  const Index n1     = jg[1].nelem();
  const Index n2     = jg[2].nelem();
  const Index n3     = jg[3].nelem();
  const Range rowind = get_rowindex_for_mblock( sensor_response, mblock_index );
  const Index row4   = rowind.get_start();
        Index col4   = jacobian_indices[iq][0];

  if( n3 > 1 )
    { col4 += mblock_index*n2*n1*2; }
      
  for( Index l=0; l<nlos; l++ )
    {
      const Index row3 = row4 + l*nf*npol;
      const Index col3 = col4 + l*n1*2;

      for( Index f=0; f<nf; f++ )
        {
          const Index row2 = row3 + f*npol;

          for( Index p=0; p<npol; p++ )
            {
              Index col1 = col3;
              if( n1 > 1 )
                { col1 += p*2; }

              jacobian(row2+p,col1)   = s[f];
              jacobian(row2+p,col1+1) = c[f];
            }
        }
    }
}



//----------------------------------------------------------------------------
// Surface quantities
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddSurfaceQuantity(
        Workspace&,
        ArrayOfRetrievalQuantity&   jq,
        Agenda&                     jacobian_agenda,
  const Index&                      atmosphere_dim,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Vector&                     rq_lat_grid,
  const Vector&                     rq_lon_grid,
  const String&                     quantity,
  const Verbosity&                  verbosity )
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  // Check that this species is not already included in the jacobian.
  for( Index it=0; it<jq.nelem(); it++ )
    {
      if( jq[it].MainTag() == SURFACE_MAINTAG  && 
          jq[it].Subtag()  == quantity )
        {
          ostringstream os;
          os << quantity << " is already included as a surface variable "
             << "in *jacobian_quantities*.";
          throw runtime_error(os.str());
        }
    }
    
  // Check retrieval grids, here we just check the length of the grids
  // vs. the atmosphere dimension
  ArrayOfVector grids( max(atmosphere_dim-1,Index(1)) );
  {
    ostringstream os;
    if( !check_retrieval_grids( grids, os, lat_grid, lon_grid,
                                rq_lat_grid, rq_lon_grid,
                                "retrieval latitude grid", 
                                "retrievallongitude_grid", 
                                atmosphere_dim ) )
    throw runtime_error(os.str());
  }
  
  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag( SURFACE_MAINTAG );
  rq.Subtag( quantity );
  rq.Analytical( 0 );
  rq.Grids( grids );

  // Add it to the *jacobian_quantities*
  jq.push_back( rq );
  
  // Add dummy
  jacobian_agenda.append( "jacobianCalcDoNothing", TokVal() );
}                    




//----------------------------------------------------------------------------
// Temperatures (atmospheric):
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddTemperature(
        Workspace&,
        ArrayOfRetrievalQuantity& jq,
        Agenda&                   jacobian_agenda,
  const Index&                    atmosphere_dim,
  const Vector&                   p_grid,
  const Vector&                   lat_grid,
  const Vector&                   lon_grid,
  const Vector&                   rq_p_grid,
  const Vector&                   rq_lat_grid,
  const Vector&                   rq_lon_grid,
  const String&                   hse,
  const String&                   method,
  const Numeric&                  dx,
  const Verbosity&                verbosity )
{
  CREATE_OUT3;
  
  // Check that temperature is not already included in the jacobian.
  // We only check the main tag.
  for (Index it=0; it<jq.nelem(); it++)
    {
      if( jq[it].MainTag() == TEMPERATURE_MAINTAG )
        {
          ostringstream os;
          os << "Temperature is already included in *jacobian_quantities*.";
          throw runtime_error(os.str());
        }
    }

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
  Index analytical;
  if( method == "perturbation" )
    { analytical = 0; }
  else if( method == "analytical" )
  { analytical = 1; }
  else
    {
      ostringstream os;
      os << "The method for atmospheric temperature retrieval can only be "
      << "\"analytical\"\n or \"perturbation\"\n.";
      throw runtime_error(os.str());
    }

  // Set subtag 
  String subtag;
  if( hse == "on" )
    { subtag = "HSE on"; }
  else if( hse == "off" )
    { subtag = "HSE off"; }
  else
    {
      ostringstream os;
      os << "The keyword for hydrostatic equilibrium can only be set to\n"
         << "\"on\" or \"off\"\n";
      throw runtime_error(os.str());
    }

  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.MainTag( TEMPERATURE_MAINTAG );
  rq.Subtag( subtag );
  rq.Mode( "abs" );
  rq.Analytical( analytical );
  rq.Perturbation( dx );
  rq.Grids( grids );
  if(analytical) {
    rq.SubSubtag(PROPMAT_SUBSUBTAG);
    rq.PropType(JacPropMatType::Temperature);
  }

  // Add it to the *jacobian_quantities*
  jq.push_back( rq );

  if( analytical ) 
    {
      out3 << "  Calculations done by semi-analytical expression.\n"; 
      jacobian_agenda.append( "jacobianCalcDoNothing", TokVal() );
    }
  else
    { 
      out3 << "  Calculations done by perturbations, of size " << dx << ".\n"; 

      jacobian_agenda.append( "jacobianCalcTemperaturePerturbations", "" );
    }
}




/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcTemperaturePerturbations(
        Workspace&                 ws,
        Matrix&                    jacobian,
  const Index&                      mblock_index,
  const Vector&                     iyb _U_,
  const Vector&                     yb,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Vector&                     lat_true,
  const Vector&                     lon_true,
  const Tensor3&                    t_field,
  const Tensor3&                    z_field,
  const Tensor4&                    vmr_field,
  const Tensor4&                    nlte_field, 
  const ArrayOfArrayOfSpeciesTag&   abs_species,
  const Vector&                     refellipsoid,
  const Matrix&                     z_surface,
  const Index&                      cloudbox_on,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const Matrix&                     sensor_pos,
  const Matrix&                     sensor_los,
  const Matrix&                     transmitter_pos,
  const Matrix&                     mblock_dlos_grid,
  const Sparse&                     sensor_response,
  const String&                     iy_unit,  
  const Agenda&                     iy_main_agenda,
  const Agenda&                     geo_pos_agenda,
  const Agenda&                     g0_agenda,
  const Numeric&                    molarmass_dry_air,
  const Numeric&                    p_hse,
  const Numeric&                    z_hse_accuracy,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const Verbosity&                  verbosity )
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
          bool any_affine;
          ArrayOfArrayOfIndex jacobian_indices;
          jac_ranges_indices( jacobian_indices, any_affine,
                              jacobian_quantities, true );
          //
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

  // Local copy of z_field. 
  Tensor3 z = z_field;

  // Loop through the retrieval grid and calculate perturbation effect
  //
  const Index    n1y = sensor_response.nrows();
        Vector   dy( n1y ); 
  const Range    rowind = get_rowindex_for_mblock( sensor_response, mblock_index ); 
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
                    perturbation_field_3d( t_p(joker,joker,joker), p_gp, 
                                           lat_gp, lon_gp, jg[0].nelem()+2,
                                           jg[1].nelem()+2, jg[2].nelem()+2, 
                                           p_range, lat_range, lon_range, 
                                           rq.Perturbation(), pertmode );
                    break;
                  }
                }

              // Apply HSE, if selected
              if( rq.Subtag() == "HSE on" )
                {
                  z_fieldFromHSE( ws,z, atmosphere_dim, p_grid, lat_grid, 
                                  lon_grid, lat_true, lon_true, abs_species, 
                                  t_p, vmr_field, refellipsoid, z_surface, 1,
                                  g0_agenda, molarmass_dry_air, 
                                  p_hse, z_hse_accuracy, verbosity );
                }
       
              // Calculate the perturbed spectrum  
              Vector        iybp;
              ArrayOfVector dummy3;      
              ArrayOfMatrix dummy4;
              Matrix        dummy5;
              //
              iyb_calc( ws, iybp, dummy3, dummy4, dummy5, mblock_index, 
                        atmosphere_dim, t_p, z,
                        vmr_field, nlte_field, cloudbox_on, 
                        stokes_dim, f_grid, sensor_pos, sensor_los, 
                        transmitter_pos, mblock_dlos_grid, 
                        iy_unit, iy_main_agenda, geo_pos_agenda,
                        0, ArrayOfRetrievalQuantity(), 
                        ArrayOfArrayOfIndex(), ArrayOfString(), verbosity );
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
// Winds:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddWind(
        Workspace&,
        ArrayOfRetrievalQuantity&   jq,
        Agenda&                     jacobian_agenda,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Vector&                     rq_p_grid,
  const Vector&                     rq_lat_grid,
  const Vector&                     rq_lon_grid,
  const String&                     component,
  const Numeric&                    dfrequency,
  const Verbosity&                  verbosity )
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  // Check that this species is not already included in the jacobian.
  for( Index it=0; it<jq.nelem(); it++ )
    {
      if( jq[it].MainTag() == WIND_MAINTAG  && 
          jq[it].Subtag()  == component )
        {
          ostringstream os;
          os << "The wind component:\n" << component << "\nis already included "
             << "in *jacobian_quantities*.";
          throw runtime_error(os.str());
        }
    }
    
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
  
  
  // Create the new retrieval quantity
  RetrievalQuantity rq;
  
  if(component == "u")
    rq.PropType(JacPropMatType::WindU);
  else if(component == "v")
    rq.PropType(JacPropMatType::WindV);
  else if(component == "w")
    rq.PropType(JacPropMatType::WindW);
  else if(component == "strength")
    rq.PropType(JacPropMatType::WindMagnitude);
  else
    throw std::runtime_error("The selection for *component* can only be \"u\", \"v\", \"w\" or \"strength\"." );
  
  rq.MainTag( WIND_MAINTAG );
  rq.Subtag( component );
  rq.Analytical( 1 );
  rq.Grids( grids );
  rq.SubSubtag(PROPMAT_SUBSUBTAG);
  rq.Perturbation( dfrequency );

  // Add it to the *jacobian_quantities*
  jq.push_back( rq );
  
  out3 << "  Calculations done by propagation matrix expression.\n"; 
  jacobian_agenda.append( "jacobianCalcDoNothing", TokVal() );
  
}                    




//----------------------------------------------------------------------------
// Magnetic field:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddMagField(
        Workspace&,
        ArrayOfRetrievalQuantity&   jq,
        Agenda&                     jacobian_agenda,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Vector&                     rq_p_grid,
  const Vector&                     rq_lat_grid,
  const Vector&                     rq_lon_grid,
  const String&                     component,
  const Numeric&                    dB,
  const Verbosity&                  verbosity )
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  // Check that this species is not already included in the jacobian.
  for( Index it=0; it<jq.nelem(); it++ )
    {
      if( jq[it].MainTag() == MAGFIELD_MAINTAG  && 
          jq[it].Subtag()  == component )
        {
          ostringstream os;
          os << "The magnetic field component:\n" << component << "\nis already "
             << "included in *jacobian_quantities*.";
          throw runtime_error(os.str());
        }
    }

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
  
  // Create the new retrieval quantity
  RetrievalQuantity rq;
  if(component == "u")
    rq.PropType(JacPropMatType::MagneticU);
  else if(component == "v")
    rq.PropType(JacPropMatType::MagneticV);
  else if(component == "w")
    rq.PropType(JacPropMatType::MagneticW);
  else if(component == "strength")
    rq.PropType(JacPropMatType::MagneticMagnitude);
  else if(component == "eta")
    rq.PropType(JacPropMatType::MagneticEta);
  else if(component == "theta")
    rq.PropType(JacPropMatType::MagneticTheta);
  else
    throw runtime_error("The selection for *component* can only be \"u\", \"v\", \"w\", \"strength\", \"eta\", or \"theta\"." );
  
  rq.MainTag( MAGFIELD_MAINTAG );
  rq.Subtag( component );
  rq.Analytical( 1 );
  rq.Grids( grids );
  
  rq.SubSubtag(PROPMAT_SUBSUBTAG);
  rq.Perturbation(dB);
  
  // Add it to the *jacobian_quantities*
  jq.push_back( rq );
  
  // Add gas species method to the jacobian agenda
  jacobian_agenda.append( "jacobianCalcDoNothing", TokVal() );
}                    




//----------------------------------------------------------------------------
// in presence of scattering:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianDoit(//WS Output:
                  Workspace& ws,
                  Vector& y0,
                  Matrix& jacobian,
                  // data that is passed to WSMs via the workspace
                  Tensor7& doit_i_field,
                  Tensor4& scat_species_mass_density_field,
                  Tensor4& scat_species_mass_flux_field,
                  Tensor4& scat_species_number_density_field,
                  Tensor4& scat_species_mean_mass_field,
                  Tensor4& pnd_field,
                  Tensor4& vmr_field,
                  Tensor3& t_field,
                  ArrayOfArrayOfSingleScatteringData& scat_data,
                  ArrayOfArrayOfScatteringMetaData& scat_meta,
                  ArrayOfString& scat_species,
                  // WS Input:
                  const ArrayOfRetrievalQuantity& jacobian_quantities,
                  const ArrayOfArrayOfSpeciesTag& abs_species,
                  const Index& atmosphere_dim,
                  const ArrayOfIndex& cloudbox_limits,
                  // input required for pnd_fieldCalcFromscat_speciesFields
                  // input required for ScatSpeciesMerge and cloudbox_checkedCalc
                  const Matrix& z_surface,
                  const Vector& p_grid,
                  // input required for DoitCalc
                  const Index& atmfields_checked,
                  const Index& atmgeom_checked,
                  const Index& cloudbox_checked,
                  const Index& scat_data_checked,
                  const Index& cloudbox_on,
                  const Vector& f_grid,
                  const Agenda& doit_mono_agenda,
                  const Index& doit_is_initialized,
                  // input required for yCalc
                  const Tensor3& z_field,
                  const Index& sensor_checked,
                  const Index& stokes_dim,
                  const Matrix& sensor_pos,
                  const Matrix& sensor_los,
                  const Matrix& transmitter_pos,
                  const Matrix& mblock_dlos_grid,
                  const Sparse& sensor_response,
                  const Vector& sensor_response_f,
                  const ArrayOfIndex& sensor_response_pol,
                  const Matrix& sensor_response_dlos,
                  const String& iy_unit,   
                  const Agenda& iy_main_agenda,
                  const Agenda& geo_pos_agenda,
                  const Agenda& jacobian_agenda,
                  const Index& jacobianDoit_do,
                  const ArrayOfString& iy_aux_vars,
                  // Keywords:
                  const Index& robust,
                  const Index& ScatSpeciesMerge_do,
                  const Index& debug,
                  const String& delim,
                  const Verbosity& verbosity)
{
/*
  // FIXME:
   1)- use jacobian_quantities.grids in perturbation level loops (currently
        relies on identiy between cloudbox_limits herein and when
        jacobianDoitAddSpecies was called)
   2)- check functionality for non-compact cases (use
        doit/TestDOITFromIndividualFields.arts as start)
       generally, improve handling of non-compact cases, which is currently done
        by an ugly workaround. not a JacobianDoit issues, though.
   3)- test (and adapt if necessary) for no-inARTS calculated pnd_fields
     - extend output: y_aux(?)
     - add some diagnostics: is input doit_i_field (and corresponding y0_1st)
        suffiently converged? check whether |y0_1st-y0_niter| << |ypert_niter-y0_niter|
     - check & disallow perturbations on scat_speciesXXfield that are not used
   4)- allow (a) smaller perturbation range in p
             (b) extend allowed perturbation range to outside cloudbox
             (c) arbitrary (in-cloudbox) p-levels
     - allow clearsky and cloudy jacobians in parallel
     - ...
*/
  CREATE_OUTS;

  if( jacobianDoit_do == 0 )
      throw std::runtime_error(
            "Doit Jacobians not switched on (use jacobianDoitClose).");

  if( jacobian_quantities.empty() )
      throw std::runtime_error(
            "No Jacobian quantities specified for DOIT Jacobian calculation.");

  if( atmosphere_dim != 1 )
      throw std::runtime_error(
            "DOIT Jacobians currently only work with a 1D atmosphere.");

  ArrayOfString fail_msg;
  bool do_abort = false;


  // yCalc will use doit_i_field via the workspace (by iyInterpCloudboxField
  // usually used in iy_cloudbox_agenda, passed to iyEmissionStandard). that is,
  // we will need to set doit_i_field, too. so we need to store the first guess
  // field in a container to pas it back as first guess for each DoitCalc
  Tensor7 doit_i_field_ref = doit_i_field;

  // for pnd_field recalculations, we need the original scat_data, not the
  // possibly merged one!
  // let's first check, whether scat_data is the original by comparing extend to
  // scat_meta extend (or likely. we can't be 100% sure).
  if( scat_data.nelem()!=scat_meta.nelem() ||
      scat_data[0].nelem()!=scat_meta[0].nelem() ||
      TotalNumberOfElements(scat_data)!=TotalNumberOfElements(scat_meta)  )
    {
      ostringstream os;
      os << "Size of scat_data and scat_meta not consistent.\n"
         << "Pass unmerged scat_data into JacobianDoit!";
      throw runtime_error(os.str());
    }
  // but, we need the pnd_fields corresponding to scat_data (i.e., if scat_data
  // is unmerged, we also need unmerged pnd_field for proper DoitCalc.
  // alternatively, we could calculate the original pnd_field again. that would
  // maybe be the better, because safer option...).
  // check, that we have that.
  if( TotalNumberOfElements(scat_data)!=pnd_field.nbooks() )
    {
      ostringstream os;
      os << "Size of scat_data and pnd_field not consistent.\n"
         << "Pass unmerged pnd_field into JacobianDoit!";
      throw runtime_error(os.str());
    }

  if( debug )
    {
      WriteXML( "ascii", pnd_field, "pnd_field_ref", 0, "pnd_field", "",
                "", verbosity );
      WriteXML( "ascii", scat_data, "scat_data_ref", 0, "scat_data", "",
                "", verbosity );
      WriteXML( "ascii", scat_meta, "scat_meta_ref", 0, "scat_meta", "",
                "", verbosity );
      WriteXML( "ascii", scat_species, "scat_species_ref", 0, "scat_species", "",
                "", verbosity );
    }


  ////// Calculation of reference case
  // We do a fixed number of iterations on top of the first-guess field from
  // outside to be consistent with the perturbation calculations.

  // if we are going to merge (i.e. to modify the scat_data), we need to keep
  // the original one. also, if we merge in the perturbation calculations, then
  // we merge here as well.
  ArrayOfArrayOfSingleScatteringData scat_data_ref;
  ArrayOfArrayOfScatteringMetaData scat_meta_ref;
  ArrayOfString scat_species_ref;
  ArrayOfTensor4 dummy_dpnd_field_dx( jacobian_quantities.nelem() );
  if( ScatSpeciesMerge_do )
    {
      Index cb_chk_internal=cloudbox_checked;
      scat_data_ref=scat_data;
      scat_meta_ref=scat_meta;
      scat_species_ref=scat_species;
      Vector latlon_dummy(0);
      Matrix part_mass_dummy(0,0);
      Tensor3 wind_dummy(0,0,0);
      ScatSpeciesMerge(	pnd_field, scat_data, scat_meta, scat_species,
                        cb_chk_internal,
                        atmosphere_dim, cloudbox_on, cloudbox_limits,
                        t_field, z_field,
                        z_surface, verbosity );
      // check that the merged scat_data still fulfills the requirements
      cloudbox_checkedCalc( cb_chk_internal, atmfields_checked,
                            atmosphere_dim, p_grid, latlon_dummy, latlon_dummy,
                            z_field, z_surface,
                            wind_dummy, wind_dummy, wind_dummy,
                            cloudbox_on, cloudbox_limits,
                            pnd_field, dummy_dpnd_field_dx, jacobian_quantities,
                            scat_data, scat_species, part_mass_dummy,
                            abs_species,
                            0, verbosity );
      if( debug )
        {
          WriteXML( "ascii", pnd_field, "pnd_field_refmerged", 0, "pnd_field",
                    "", "", verbosity );
          WriteXML( "ascii", scat_data, "scat_data_refmerged", 0, "scat_data",
                    "", "", verbosity );
          WriteXML( "ascii", scat_meta, "scat_meta_refmerged", 0, "scat_meta",
                    "", "", verbosity );
          WriteXML( "ascii", scat_species, "scat_species_refmerged", 0, "scat_species",
                    "", "", verbosity );
        }
    }
  DoitCalc( ws, doit_i_field,
            atmfields_checked, atmgeom_checked,
            cloudbox_checked, scat_data_checked,
            cloudbox_on, f_grid, doit_mono_agenda, doit_is_initialized,
            verbosity );

  Vector vec_dummy;
  ArrayOfIndex aoi_dummy;
  Matrix mat_dummy1, mat_dummy2, mat_dummy3, mat_dummy4;
  ArrayOfVector aov_dummy;
  Tensor4 nlte_field(0,0,0,0);
  yCalc( ws, y0,
         vec_dummy, aoi_dummy, mat_dummy1, mat_dummy2, aov_dummy,
         mat_dummy3, mat_dummy4,
         atmgeom_checked, atmfields_checked, atmosphere_dim,
         t_field, z_field, vmr_field, nlte_field, cloudbox_on,
         cloudbox_checked, scat_data_checked, sensor_checked,
         stokes_dim, f_grid,
         sensor_pos, sensor_los, transmitter_pos, mblock_dlos_grid,
         sensor_response, sensor_response_f, sensor_response_pol,
         sensor_response_dlos, iy_unit, iy_main_agenda, geo_pos_agenda,
         jacobian_agenda, 0, jacobian_quantities,
         iy_aux_vars, verbosity );


  ////// Now we start with the perturbations runs
  // per perturbation species and perturbation level we need to:
  // 1) perturb atmo
  // 1a) recalculated pnd_field for perturbed atmosphere
  // 2) DoitCalc
  // 3) yCalc
  // 4) Calculate jacobian = (y-y0)/perturbation and append to jacobian matrix

  // We need to safe the reference fields. This because we use DoitCalc, which
  // accesses the atmospheric state through the workspace. I.e., we
  // have to overwrite the workspace fields with the perturbed fields, then
  // restore them to the original state for the next perturbation level and/or
  // species.
  // This sounds like a lot of back-and-forth copying. Can we do this in a
  // better way? (Oliver says: not really)
  Tensor4 vmr_field_ref = vmr_field;
  Tensor3 t_field_ref = t_field;
  Tensor4 scat_species_mass_density_field_ref = scat_species_mass_density_field;
  Tensor4 scat_species_mass_flux_field_ref = scat_species_mass_flux_field;
  Tensor4 scat_species_number_density_field_ref = scat_species_number_density_field;
  Tensor4 scat_species_mean_mass_field_ref = scat_species_mean_mass_field;


  // Set some useful variables. 
  RetrievalQuantity jq;
  Index it=0;
  Index lstart, lend; //, pertmode;
  String speciesname;

  // as long as we limit the perturbation to within cloudbox and provide a
  // ready-made doit_i_field (were we do not modify the clear incoming part
  // again), we can only perturb until the last level below the upper
  // cloudbox limit (cause else we also modify the outside-cloudbox state,
  // i.e., we'd need a new clear incoming calc) unless cloudbox reaches till
  // TOA.
  if( cloudbox_limits[1]==p_grid.nelem() )
    lend = cloudbox_limits[1]+1;
  else
    lend = cloudbox_limits[1];
  // this also applies for lower limit, IF lower limit is not on the
  // surface.
  if( cloudbox_limits[0]==0 )
    lstart = 0;
  else
    lstart = cloudbox_limits[0]+1;
  Index np = lend-lstart;

  // for now we do the perturbations for ALL species on ALL the p-levels within
  // the cloudbox.
  // however, for making it easier adaptable to specific retrieval grids, we use
  // the same basic strategy using ArrayOfGridPos as for clearsky perturbation
  // jacobians (see jacobianCalcAbsSpeciesPerturbations in m_jacobians.cc)
  // ok, didn't work on first try. using the easy version (loop over
  // p_grid-levels) now. adapt to clearsky-equivalent use later on (when we
  // allow retrieval/perturbation grids different from the given p_grid)...
  
  //ArrayOfGridPos p_gp;
  //Vector jg_p = p_grid[Range(lstart,lend-lstart)]; // using correct extend?
                                                   // lend should NOT be included
  //Index np   = jg_p.nelem();
  //get_perturbation_gridpos( p_gp, p_grid, jg_p, true );

  // Initialize jacobian matrix
  //jacobian.resize(y0.nelem(), jacobian_quantities.nelem()*np);
  //
  ArrayOfArrayOfIndex jacobian_indices;
  {
    bool any_affine;
    jac_ranges_indices( jacobian_indices, any_affine,
                        jacobian_quantities, true );
  }
  //
  jacobian.resize(y0.nelem(),
                  jacobian_indices[jacobian_indices.nelem()-1][1]+1);
  jacobian = NAN;

  // loop over all perturbation species (aka jacobian quantities)
  for( Index iq=0; iq<jacobian_quantities.nelem(); iq++ )
    {
      if (do_abort) continue;

      // before start perturbing we need to put the original scat_* data back in
      // place. at least scat_species (rest comes later).
      if( ScatSpeciesMerge_do )
        scat_species=scat_species_ref;

      jq = jacobian_quantities[iq];
      Index si=-1;


      // check if iterator 'it' is consistent with jacobian_indices entry of the
      // species
      //assert( it == jacobian_indices[iq][0] );

      // Check if a relative perturbation is used or not, this information is needed
      //by the methods 'perturbation_field_?d'.
      //if( jq.Mode()=="rel" )
      //  pertmode = 0;
      //else 
      //  pertmode = 1;

      // check if perturbation species is valid. and extract the actual field we
      // are going to perturb.
      // first determine, which species type:
      if( jq.MainTag() == ABSSPECIES_MAINTAG )
        {
          // Find VMR field for this species. 
          ArrayOfSpeciesTag tags;
          array_species_tag_from_string( tags, jq.Subtag() );
          si = chk_contains( "abs_species", abs_species, tags );
          speciesname = jq.MainTag()+'.'+jq.Subtag();
        }
      else if( jq.MainTag() == SCATSPECIES_MAINTAG )
        {
          // we know, it's a scat_species. so, next we need to check, which
          // scat_species (or hydrometeor type) it is. for that, we need to
          // compare to scat_species entries.
          Index i=0;
          while( i<scat_species.nelem() && si<0 )
            {
              String scat_species_name;
              parse_partfield_name(scat_species_name, scat_species[i], delim);
              if( scat_species_name == jq.Subtag() )
                  si = i;
              i++;
            }
          if( si<0 )
            {
              ostringstream os;
              os << "scat_species does not contain " << jq.Subtag();
              throw runtime_error(os.str());
            }
          // whether a generally valid field has been selected was checked in
          // jacobianDoitAddSpecies. it's left to check, whether this field
          // contains valid values. 0 is perfectly valid. NaN is not.
          // Since this check is done within pnd_fieldCalcFromscat_speciesFields
          // by the pnd_fieldX methods, we skip that here.
          // Would also be good to check, whether the perturbed field is actually
          // applied by the PSD selected for this scat_species. However,
          // currently there is no way to check this (as we have no info here,
          // which PSD requires which field).
          // FIXME: check that to-be-perturbed field is of non-zero size (if no
          // scat_species has any entry, we don't set the scat_speciesXXfield at
          // all).
          // FIXME: we could add that check in the pnd_fieldX methods and give a
          // warning there (we can't throw error there as we currently can't
          // circumvent the extraction of fields from compact data. as for doing
          // this, we'd also need the info what fields are required for selected
          // PSD). or to allow for an error here (and/or to exclude extraction
          // of not needed fields in AtmFieldsFromCompact), we could built up an
          // internal variable that holds that information.
          // FIXME: when using basic atm fields instead of compact atmos, we
          // currently use dummy empty profiles in place of non-existing data,
          // i.e. set the fields to 0. when fixing the above (throwing error,
          // when non-NaN data provided for un-used fields), we need an
          // alternative way to buffer the non-existing data (if NaN is
          // interpolable (check!), then we can read NaN instead of 0 fields. but
          // maybe a proper WSM for setting the fields from non-compact data is
          // nicer...
          speciesname = jq.MainTag()+'.'+jq.Subtag()+'.'+jq.SubSubtag();
        }
      else if( jq.MainTag() == TEMPERATURE_MAINTAG )
        {
          speciesname = jq.MainTag();
        }
      else
        {
          ostringstream os;
          os << jq.MainTag() << " is an unknown Doit jacobian species.\n"
             << "Use jacobianDoitAddSpecies to properly set up cloudy-sky "
             << "jacobians.";
          throw runtime_error(os.str());
        }

      // loop over all perturbation levels
      for( Index il=lstart; il<lend; il++ )
        {
          out2 << "handling level #" << il << " of pert species " << speciesname
               << "\n";

          if (do_abort) continue;
          bool do_doit = true;

          // before start perturbing we need to put the original scat_* data
          // back in place.
          if( ScatSpeciesMerge_do )
            {
              scat_data=scat_data_ref;
              scat_meta=scat_meta_ref;
              scat_species=scat_species_ref;
            }

/*
          //use this if we once allow arbitrary perturbation grids. if so:
          //- correct for last-point outdrag
          //- check for proper use of pertmode (passed to perturbation_field_1d)
          Range p_range   = Range(0,0);
          get_perturbation_range( p_range, il, j_p );
          Tensor3 pertfield;

          if( jq.MainTag() == ABSSPECIES_MAINTAG )
            {
              pertfield = vmr_field(si,joker,joker,joker);
            }
          else if( jq.MainTag() == SCATSPECIES_MAINTAG )
            {
              ostringstream os;
              os << "Oops. Scatterig species perturbations not yet available.";
              throw runtime_error(os.str());
            }
          else //temperature
            {
              pertfield = t_field;
              WriteXMLIndexed( "ascii", il, t_field, "tfield",
                          "tfield", "", verbosity );
            }

          perturbation_field_1d( pertfield(joker, 0, 0), 
                                 p_gp, jg_p.nelem()+2, p_range, //extend correct?
                                 jq.Perturbation(), pertmode );

          // pasting pertured field back into the calculation field
          if( jq.MainTag() == ABSSPECIES_MAINTAG )
            {
              vmr_field(si,joker,joker,joker) = pertfield;
            }
          else if( jq.MainTag() == SCATSPECIES_MAINTAG )
            {
              ostringstream os;
              os << "Oops. Scattering species perturbations not yet available.";
              throw runtime_error(os.str());
            }
          else //temperature
            {
              t_field = pertfield;
              WriteXMLIndexed( "ascii", il, t_field, "pert_tfield",
                          "tfield", "", verbosity );
            }
*/
          if( jq.MainTag() == ABSSPECIES_MAINTAG )
            {
              if( jq.Mode() == "abs" )
                vmr_field(si,il,joker,joker) += jq.Perturbation();
              else if ( jq.Mode() == "rel" )
                if (vmr_field(si,il,0,0)==0.)
                  do_doit = false;
                else
                  vmr_field(si,il,joker,joker) *= (1.+jq.Perturbation());
              else
                // we shouldn't end up here. if we do, checks for allowed
                // retrieval modes above are incomplete.
                assert( 0 );
            }
          else if( jq.MainTag() == SCATSPECIES_MAINTAG )
            {
              if( jq.SubSubtag() == "mass_density" )
                {
                  if( jq.Mode() == "abs" )
                    {
                      scat_species_mass_density_field(si,il,joker,joker)
                        += jq.Perturbation();
                    }
                  else if ( jq.Mode() == "rel" )
                    {
                      // ATTENTION: in case ever allowing other than 1D, all
                      // these checks (in the rel branch) needs to be updated to
                      // check whether ALL lat/lon entries are 0.
                      if (scat_species_mass_density_field(si,il,0,0)==0.)
                        do_doit = false;
                      else
                        scat_species_mass_density_field(si,il,joker,joker)
                          *= (1.+jq.Perturbation());
                    }
                  else
                    // we shouldn't end up here. if we do, checks for allowed
                    // retrieval modes above are incomplete.
                    assert( 0 );
                }
              else if( jq.SubSubtag() == "mass_flux" )
                {
                  if( jq.Mode() == "abs" )
                    {
                      scat_species_mass_flux_field(si,il,joker,joker)
                        += jq.Perturbation();
                    }
                  else if ( jq.Mode() == "rel" )
                    {
                      if (scat_species_mass_flux_field(si,il,0,0)==0.)
                        do_doit = false;
                      else
                        scat_species_mass_flux_field(si,il,joker,joker)
                          *= (1.+jq.Perturbation());
                    }
                  else
                    // we shouldn't end up here. if we do, checks for allowed
                    // retrieval modes above are incomplete.
                    assert( 0 );
                }
              else if( jq.SubSubtag() == "number_density" )
                {
                  if( jq.Mode() == "abs" )
                    {
                      scat_species_number_density_field(si,il,joker,joker)
                        += jq.Perturbation();
                    }
                  else if ( jq.Mode() == "rel" )
                    {
                      if (scat_species_number_density_field(si,il,0,0)==0.)
                        do_doit = false;
                      else
                        scat_species_number_density_field(si,il,joker,joker)
                          *= (1.+jq.Perturbation());
                    }
                  else
                    // we shouldn't end up here. if we do, checks for allowed
                    // retrieval modes above are incomplete.
                    assert( 0 );
                }
              else if( jq.SubSubtag() == "mean_mass" )
                {
                  if( jq.Mode() == "abs" )
                    {
                      scat_species_mean_mass_field(si,il,joker,joker)
                        += jq.Perturbation();
                    }
                  else if ( jq.Mode() == "rel" )
                    {
                      if (scat_species_mean_mass_field(si,il,0,0)==0.)
                        do_doit = false;
                      else
                        scat_species_mean_mass_field(si,il,joker,joker)
                          *= (1.+jq.Perturbation());
                    }
                  else
                    // we shouldn't end up here. if we do, checks for allowed
                    // retrieval modes above are incomplete.
                    assert( 0 );
                }
            }
          else //temperature
            {
              t_field(il,joker,joker) += jq.Perturbation();
            }

          if ( do_doit )
          {
            try
            {
              // unless pnd_field is NOT calculated inside ARTS, we have to
              // recalculate it. not only for scat_speciesXXfield perturbances.
              // this since also other parameters could effect the pnd_field,
              // e.g., atmospheric temperature.
              // not straight forward, how we can check for whether pnd_field is
              // from external. but a good guess is that then scat_speciesXXfields
              // are not required, i.e. are likely empty. so, if not empty (here:
              // sized 0!), we try to recalculate them.
              if( scat_species_mass_density_field.npages()!=0 ||
                  scat_species_mass_flux_field.npages()!=0 ||
                  scat_species_number_density_field.npages()!=0 ||
                  scat_species_mean_mass_field.npages()!=0 )
                {
                  pnd_fieldCalcFromscat_speciesFields(
                    pnd_field, dummy_dpnd_field_dx,
                    atmosphere_dim, cloudbox_on, cloudbox_limits,
                    scat_species_mass_density_field, scat_species_mass_flux_field,
                    scat_species_number_density_field, scat_species_mean_mass_field,
                    t_field, scat_meta, scat_species, jacobian_quantities,
                    delim, verbosity );
                  if( debug )
                    {
                      WriteXMLIndexed( "ascii", iq*np+il, pnd_field,
                                       "pnd_field_perturbed", 0, "pnd_field", "", "",
                                       verbosity );
                    }
                  if( ScatSpeciesMerge_do )
                    {
                      //scat_data=scat_data_ref;
                      //scat_meta=scat_meta_ref;
                      //scat_species=scat_species_ref;
                      Index cb_chk_internal=cloudbox_checked;
                      Vector latlon_dummy(0);
                      Matrix part_mass_dummy(0,0);
                      Tensor3 wind_dummy(0,0,0);
                      ScatSpeciesMerge( pnd_field,
                                        scat_data, scat_meta, scat_species,
                                        cb_chk_internal,
                                        atmosphere_dim,
                                        cloudbox_on, cloudbox_limits,
                                        t_field, z_field,
                                        z_surface, verbosity );
                      cloudbox_checkedCalc( cb_chk_internal,
                                            atmfields_checked, atmosphere_dim,
                                            p_grid, latlon_dummy, latlon_dummy,
                                            z_field, z_surface,
                                            wind_dummy, wind_dummy, wind_dummy,
                                            cloudbox_on, cloudbox_limits,
                                            pnd_field, dummy_dpnd_field_dx,
                                            jacobian_quantities,
                                            scat_data, scat_species, 
                                            part_mass_dummy, abs_species,
                                            0, verbosity );
                      if( debug )
                        {
                          WriteXMLIndexed( "ascii", iq*np+il, scat_data,
                                           "scat_data_mergeperturbed", 0, "scat_data", "", "",
                                           verbosity );
                          WriteXMLIndexed( "ascii", iq*np+il, scat_meta,
                                           "scat_meta_mergeperturbed", 0, "scat_meta", "", "",
                                           verbosity );
                          WriteXMLIndexed( "ascii", iq*np+il, scat_species,
                                           "scat_species_mergeperturbed", 0, "scat_species", "", "",
                                           verbosity );
                          WriteXMLIndexed( "ascii", iq*np+il, pnd_field,
                                           "pnd_field_mergeperturbed", 0, "pnd_field", "", "",
                                           verbosity );
                        }
                    }
                }

              if( debug )
                {
                  WriteXMLIndexed( "ascii", iq*np+il, scat_data,
                                   "scat_data_final", 0, "scat_data", "", "",
                                    verbosity );
                  WriteXMLIndexed( "ascii", iq*np+il, scat_meta,
                                   "scat_meta_final", 0, "scat_meta", "", "",
                                    verbosity );
                  WriteXMLIndexed( "ascii", iq*np+il, scat_species,
                                   "scat_species_final", 0, "scat_species", "", "",
                                    verbosity );
                  WriteXMLIndexed( "ascii", iq*np+il, pnd_field,
                                   "pnd_field_final", 0, "pnd_field", "", "",
                                   verbosity );
                }
              doit_i_field = doit_i_field_ref;
              DoitCalc( ws, doit_i_field,
                        atmfields_checked, atmgeom_checked,
                        cloudbox_checked, scat_data_checked,
                        cloudbox_on, f_grid, doit_mono_agenda, doit_is_initialized,
                        verbosity );
              if( debug )
                {
                  WriteXMLIndexed( "ascii", iq*np+il, doit_i_field,
                                   "ifield_perturbed", 0, "doit_i_field", "", "",
                                    verbosity );
                }

              Vector y;
              yCalc( ws, y,
                     vec_dummy, aoi_dummy, mat_dummy1, mat_dummy2, aov_dummy,
                     mat_dummy3, mat_dummy4,
                     atmgeom_checked, atmfields_checked, atmosphere_dim,
                     t_field, z_field, vmr_field, nlte_field, cloudbox_on,
                     cloudbox_checked, scat_data_checked, sensor_checked,
                     stokes_dim, f_grid,
                     sensor_pos, sensor_los, transmitter_pos, mblock_dlos_grid,
                     sensor_response, sensor_response_f, sensor_response_pol,
                     sensor_response_dlos, iy_unit, iy_main_agenda, geo_pos_agenda,
                     jacobian_agenda, 0, jacobian_quantities, 
                     iy_aux_vars, verbosity );

              if( debug )
                {
                  WriteXMLIndexed( "ascii", iq*np+il, y, "y", 0, "y", "", "", verbosity );
                }

              Vector dydx(y0.nelem());
              for( Index i=0; i<y0.nelem(); i++ )
                {
                  dydx[i] = (y[i]-y0[i]) / jq.Perturbation();
                }

              jacobian(joker,it) = dydx;
            }
            catch (runtime_error e)
            {
              if( robust )
              {
                // Don't fail full calc if one of the perturbation calcs went
                // wrong.
                ostringstream os;
                os << "WARNING! Jacobian calculation for " << speciesname
                   << " at level " << il << " failed.\n"
                   << "jacobian matrix will contain NaN for this job.\n"
                   << "The runtime error produced was:\n"
                   << e.what() << "\n";
                out0 << os.str();
              }
              else
              {
                // The user wants the batch job to fail if one of the
                // jobs goes wrong.
                do_abort = true;
                ostringstream os;
                os << "Jacobian calculation for " << speciesname
                   << " at level " << il << " failed. Aborting...\n";
                out1 << os.str();
              }
              ostringstream os;
              os << "Run-time error at jacobianDoit species " << speciesname
                 << ", level " << il << ": \n" << e.what();
              fail_msg.push_back(os.str());
            }

            it++;

            // we need to restore the original atm fields again for to start from
            // original field again for next perturbation level (or species)
            if( jq.MainTag() == ABSSPECIES_MAINTAG )
              {
                vmr_field(si,joker,joker,joker) = vmr_field_ref(si,joker,joker,joker);
              }
            else if( jq.MainTag() == SCATSPECIES_MAINTAG )
              {
                scat_species_mass_density_field(si,joker,joker,joker) =
                  scat_species_mass_density_field_ref(si,joker,joker,joker);
                scat_species_mass_flux_field(si,joker,joker,joker) =
                  scat_species_mass_flux_field_ref(si,joker,joker,joker);
                scat_species_number_density_field(si,joker,joker,joker) =
                  scat_species_number_density_field_ref(si,joker,joker,joker);
                scat_species_mean_mass_field(si,joker,joker,joker) =
                  scat_species_mean_mass_field_ref(si,joker,joker,joker);
              }
            else //temperature
              {
                t_field = t_field_ref;
              }
          }
          else
          {
            Vector dydx(y0.nelem());
            for( Index i=0; i<y0.nelem(); i++ )
              {
                dydx[i] = 0.;
              }
            jacobian(joker,it) = dydx;
            it++;
          }

          if( debug )
            {
              // to have info on the jacobians already along the way...
              WriteXML( "ascii", jacobian, "jacobian", 0, "jacobian", "", "",
                         verbosity );
            }
        }
    }

  if (fail_msg.nelem())
    {
      ostringstream os;

      if (!do_abort) os << "\nError messages from failed jacobianDoit cases:\n";
      for (ArrayOfString::const_iterator cit = fail_msg.begin();
           cit != fail_msg.end(); cit++)
          os << *cit << '\n';

      if (do_abort)
          throw runtime_error(os.str());
      else
          out0 << os.str();
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianDoitClose(
        Index&                     jacobianDoit_do,
  const Index&                     jacobian_do,
  const ArrayOfRetrievalQuantity&  jacobian_quantities,
  const Verbosity&                 /* verbosity */ )
{
  // Make sure we're not trying to do both clearsky and cloudy Jacobians
  if( jacobian_do != 0 )
    throw runtime_error(
          "Currently not possible to combine clearksy and DOIT Jacobians.");

  // Make sure that the array is not empty
  if( jacobian_quantities.empty() )
    throw runtime_error(
          "No retrieval quantities has been added to *jacobian_quantities*." );

  jacobianDoit_do = 1;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianDoitAddSpecies(//WS Output:
                            ArrayOfRetrievalQuantity& jacobian_quantities,
                            // WS Input:
                            const Index& jacobian_do,
                            const Index& atmosphere_dim,
                            const Vector& p_grid,
                            const Vector& lat_grid,
                            const Vector& lon_grid,
                            const ArrayOfIndex& cloudbox_limits,
                            // Keywords:
                            const String& species,
                            const String& mode, //abs or rel
                            const Numeric& dx,
                            const Verbosity& /* verbosity */)
{
  // Make sure we're not trying to do both clearsky and cloudy Jacobians
  if( jacobian_do != 0 )
    throw runtime_error(
          "Currently not possible to combine clearksy and DOIT Jacobians.");

  // Create the new retrieval quantity
  RetrievalQuantity rq;

  // We don't allow zero or negative perturbations.
  if( dx<=0. )
    {
      ostringstream os;
      os << "Perturbations must be >0.\n"
         << "For " << species << " yours is " << dx << ".";
      throw runtime_error(os.str());
    }

  ArrayOfString strarr;
  // would rather like '-' as a delimiter.but for abs_species we need everything
  // except strarr[0] to go into the Subtag (e.g. for H20-PWR98 in abs_species,
  // Subtag must contain the PWR98 part. just H2O is not enough unless
  // abs_species contains H2O or H2O-*-*-*.).
  species.split( strarr, "." );

  if( strarr.size()>0 )
    {
      //first entry need to be "t", "abs_species" or "scat_species"
      if( strarr[0]=="T" || strarr[0]=="t" )
        {
          // Check that temperature is not already included in the jacobian.
          for( Index it=0; it<jacobian_quantities.nelem(); it++ )
            {
              if( jacobian_quantities[it].MainTag() == TEMPERATURE_MAINTAG )
                {
                  ostringstream os;
                  os << "Temperature is already included in *jacobian_quantities*.";
                  throw runtime_error(os.str());
                }
            }

          // Only abs perturbances allowed for temperature
          if( mode != "abs" )
            {
              ostringstream os;
              os << "Only absolute perturbances (mode='abs') allowed for "
                 << "temperature jacobians.";
              throw runtime_error(os.str());
            }

          // consider HSE on/off here? if so, do by subtag (see
          // jacobianAddTemperature)
          rq.MainTag( TEMPERATURE_MAINTAG );
          rq.Mode( "abs" );
          rq.Perturbation( dx );
        }

      else if( strarr[0]=="abs_species" )
        {
          if( strarr.size()>1 )
            {
              // Check that this species is not already included in the jacobian.
              for( Index it=0; it<jacobian_quantities.nelem(); it++ )
                {
                  if( jacobian_quantities[it].MainTag() == ABSSPECIES_MAINTAG  && 
                      jacobian_quantities[it].Subtag() == strarr[1] )
                    {
                      ostringstream os;
                      os << "The gas species:\n" << strarr[1] << " is already "
                         << "included in *jacobian_quantities*.";
                      throw runtime_error(os.str());
                    }
                }
            }
          else
            {
              ostringstream os;
              os << "No species tag for absorption given.";
              throw runtime_error(os.str());
            }

          if( mode != "abs" && mode != "rel" )
            {
              ostringstream os;
              os << mode << " is not a valid perturbation mode. "
                 << "Only 'abs' and 'rel' allowed.";
              throw runtime_error(os.str());
            }

          rq.MainTag( ABSSPECIES_MAINTAG );
          rq.Subtag( strarr[1] );
          rq.Mode( mode );
          rq.Perturbation( dx );
        }

      else if( strarr[0]=="scat_species" )
        {
          if( strarr.size()>1 )
            {
              if( strarr.size()>2 )
                {
                  // FIXME: we should also allow pnd_field (all? single scatt
                  // elements?)
                  if( strarr[2] == "mass_density" || strarr[2] == "mass_flux" ||
                      strarr[2] == "mean_mass" || strarr[2] == "number_density" )
                    {
                      // Check that this species&field combi is not already
                      // included in the jacobian.
                      for( Index it=0; it<jacobian_quantities.nelem(); it++ )
                        {
                          if( jacobian_quantities[it].MainTag() == SCATSPECIES_MAINTAG  && 
                              jacobian_quantities[it].Subtag() == strarr[1] &&
                              jacobian_quantities[it].SubSubtag()  == strarr[2])
                            {
                              ostringstream os;
                              os << "The " << strarr[2] << " field of "
                                 << "scattering species " << strarr[1] << "\n"
                                 << "is already included in "
                                 << "*jacobian_quantities*.";
                              throw runtime_error(os.str());
                            }
                        }
                    }
                  else
                    {
                      ostringstream os;
                      os << strarr[2] << " is not a valid scattering species "
                         << "field tag";
                      throw runtime_error(os.str());
                    }
                }
              else
                {
                  ostringstream os;
                  os << "No field tag for scattering species given.";
                  throw runtime_error(os.str());
                }
            }
          else
            {
              ostringstream os;
              os << "No species tag for scattering species given.";
              throw runtime_error(os.str());
            }

          if( mode != "abs" && mode != "rel" )
            {
              ostringstream os;
              os << mode << " is not a valid perturbation mode. "
                 << "Only 'abs' and 'rel' allowed.";
              throw runtime_error(os.str());
            }

          rq.MainTag( SCATSPECIES_MAINTAG );
          rq.Subtag( strarr[1] );
          rq.SubSubtag( strarr[2] );
          rq.Mode( mode );
          rq.Perturbation( dx );
        }
      else
        {
          ostringstream os;
          os << strarr << " is not a valid jacobianDOIT species.";
          throw runtime_error(os.str());
        }
    }
  else
    {
      ostringstream os;
      os << "No species string given.";
      throw runtime_error(os.str());
    }

  // Perturbations are done over the cloudbox. But not on the outermost grid
  // points (as this would also imply value changes outside the box, through the
  // linear interpolation between gridpoints. and we don't allow this since we
  // don't yet allow incoming field modifications!) unless outermost is at
  // surface or TOA.
  ArrayOfVector grids(atmosphere_dim);
  Vector rq_p_grid(0), rq_lat_grid(0), rq_lon_grid(0);

  Index range_start, // surface-closest point with perturb
        range_end,   // space-closest point with perturb
        range_extent;

  if( cloudbox_limits[0]==0 )
    range_start = cloudbox_limits[0];
  else
    range_start = cloudbox_limits[0]+1;
  if( cloudbox_limits[1]==p_grid.nelem() )
    range_end = cloudbox_limits[1];
  else
    range_end = cloudbox_limits[1]-1;
  range_extent = range_end - range_start + 1;

  rq_p_grid = p_grid[ Range(range_start,range_extent) ];

  if( atmosphere_dim>1 )
    {
      rq_lat_grid = lat_grid[ Range(cloudbox_limits[2]+1,
                                    cloudbox_limits[3]-cloudbox_limits[2] ) ];
      if( atmosphere_dim>2 )
          rq_lon_grid = lon_grid[ Range(cloudbox_limits[4]+1,
                                        cloudbox_limits[5]-cloudbox_limits[4]) ];
    }

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
  rq.Grids( grids );

  // Add it to the *jacobian_quantities*
  jacobian_quantities.push_back( rq );
}




//----------------------------------------------------------------------------
// Catalog parameters:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddCatalogParameter(
    Workspace&,
    ArrayOfRetrievalQuantity&   jq,
    Agenda&                     jacobian_agenda,
    const QuantumIdentifier&    catalog_identity,
    const String&               catalog_parameter,
    const Verbosity&            verbosity )
{
    CREATE_OUT3;
    
    // Check that this is not already included in the jacobian.
    for( Index it=0; it<jq.nelem(); it++ )
    {
        if( jq[it].MainTag() == CATALOGPARAMETER_MAINTAG  && 
            jq[it].QuantumIdentity()  == catalog_identity && 
            jq[it].Mode()  == catalog_parameter )
        {
            ostringstream os;
            os << "The catalog identifier:\n" << catalog_identity<< "\nis already included in "
            << "*jacobian_quantities*.";
            throw std::runtime_error(os.str());
        }
    }
    
    // Create the new retrieval quantity
    RetrievalQuantity rq;
    
    // Check catalog_parameter here
    if(LINESTRENGTH_MODE                   == catalog_parameter) rq.PropType(JacPropMatType::LineStrength);
    else if(LINECENTER_MODE                == catalog_parameter) rq.PropType(JacPropMatType::LineCenter);
    else if(SELFBROADENING_MODE            == catalog_parameter) rq.PropType(JacPropMatType::LineGammaSelf);
    else if(SELFPRESSURESHIFT_MODE         == catalog_parameter) rq.PropType(JacPropMatType::LineShiftSelf);
    else if(SELFBROADENINGEXPONENT_MODE    == catalog_parameter) rq.PropType(JacPropMatType::LineGammaSelfExp);
    else if(FOREIGNBROADENING_MODE         == catalog_parameter) rq.PropType(JacPropMatType::LineGammaForeign);
    else if(FOREIGNPRESSURESHIFT_MODE      == catalog_parameter) rq.PropType(JacPropMatType::LineShiftForeign);
    else if(FOREIGNBROADENINGEXPONENT_MODE == catalog_parameter) rq.PropType(JacPropMatType::LineGammaForeignExp);
    else if(WATERBROADENING_MODE           == catalog_parameter) rq.PropType(JacPropMatType::LineGammaWater);
    else if(WATERPRESSURESHIFT_MODE        == catalog_parameter) rq.PropType(JacPropMatType::LineShiftWater);
    else if(WATERBROADENINGEXPONENT_MODE   == catalog_parameter) rq.PropType(JacPropMatType::LineGammaWaterExp);
    else if(LINEMIXINGY0_MODE              == catalog_parameter) rq.PropType(JacPropMatType::LineMixingY0);
    else if(LINEMIXINGY1_MODE              == catalog_parameter) rq.PropType(JacPropMatType::LineMixingY1);
    else if(LINEMIXINGYEXPONENT_MODE       == catalog_parameter) rq.PropType(JacPropMatType::LineMixingYExp);
    else if(LINEMIXINGG0_MODE              == catalog_parameter) rq.PropType(JacPropMatType::LineMixingG0);
    else if(LINEMIXINGG1_MODE              == catalog_parameter) rq.PropType(JacPropMatType::LineMixingG1);
    else if(LINEMIXINGGEXPONENT_MODE       == catalog_parameter) rq.PropType(JacPropMatType::LineMixingGExp);
    else if(LINEMIXINGDF0_MODE             == catalog_parameter) rq.PropType(JacPropMatType::LineMixingG0);
    else if(LINEMIXINGDF1_MODE             == catalog_parameter) rq.PropType(JacPropMatType::LineMixingG1);
    else if(LINEMIXINGDFEXPONENT_MODE      == catalog_parameter) rq.PropType(JacPropMatType::LineMixingDFExp);
    else {
      ostringstream os;
      os << "You have selected:\n" << catalog_parameter << "\nas your catalog parameter. This is not supported.\n" 
          << "Please see user guide for supported parameters.\n";
          throw std::runtime_error(os.str());
    }
    
    
    rq.MainTag( CATALOGPARAMETER_MAINTAG );
    rq.Mode( catalog_parameter );
    rq.QuantumIdentity(catalog_identity);
    rq.Analytical(1);
    rq.SubSubtag(PROPMAT_SUBSUBTAG);
    rq.IntegrationOn();
    
    // Add it to the *jacobian_quantities*
    jq.push_back( rq );
    
    out3 << "  Calculations done by propagation matrix expressions.\n"; 
    
    jacobian_agenda.append( "jacobianCalcDoNothing", TokVal() );
}    



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddCatalogParameters(
    Workspace&                  ws,
    ArrayOfRetrievalQuantity&   jq,
    Agenda&                     jacobian_agenda,
    const ArrayOfQuantumIdentifier&    catalog_identities,
    const ArrayOfString&               catalog_parameters,
    const Verbosity&            verbosity )
{
    CREATE_OUT2;
    
    out2 << " Adding "<<catalog_identities.nelem()*catalog_parameters.nelem()
    <<" expression(s) to the Jacobian calculations.\n";
    
    for (Index ici = 0; ici<catalog_identities.nelem(); ici++)
    {
        for(Index icp = 0; icp<catalog_parameters.nelem(); icp++)
        {
            out2<<"    type: "<<catalog_parameters[icp]<<"; identifier: "<<catalog_identities[ici]<<"\n";
            jacobianAddCatalogParameter(ws, jq, jacobian_agenda,
                catalog_identities[ici], catalog_parameters[icp],
                verbosity );
        }
    }
}  



//----------------------------------------------------------------------------
// NLTE temperature:
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddNLTE(
    Workspace&,
    ArrayOfRetrievalQuantity&   jq,
    Agenda&                     jacobian_agenda,
    const Index&                atmosphere_dim,
    const Vector&               p_grid,
    const Vector&               lat_grid,
    const Vector&               lon_grid,
    const Vector&               rq_p_grid,
    const Vector&               rq_lat_grid,
    const Vector&               rq_lon_grid,
    const QuantumIdentifier&    energy_level_identity,
    const Numeric&              dx,
    const String&               mode,
    const Verbosity&            verbosity )
{
    CREATE_OUT3;
    
    // Check that this species is not already included in the jacobian.
    for( Index it=0; it<jq.nelem(); it++ )
    {
        if( jq[it].MainTag() == NLTE_MAINTAG and jq[it].QuantumIdentity()  == energy_level_identity )
        {
            ostringstream os;
            os << "The NLTE identifier:\n" << energy_level_identity<< "\nis already included in "
            << "*jacobian_quantities*.";
            throw std::runtime_error(os.str());
        }
    }
    
    // Check retrieval grids, here we just check the length of the grids
    // vs. the atmosphere dimension
    ArrayOfVector grids(atmosphere_dim);
    {
        ostringstream os;
        if(not check_retrieval_grids(grids, os, p_grid, lat_grid, lon_grid, rq_p_grid, rq_lat_grid, rq_lon_grid,
            "retrieval pressure grid",  "retrieval latitude grid",  "retrievallongitude_grid",  atmosphere_dim ) )
            throw runtime_error(os.str());
    }
    
    
    // Create the new retrieval quantity
    RetrievalQuantity rq;
    
    // Set the mode
    if(mode == "Tv") rq.PropType(JacPropMatType::VibrationalTemperature);
    else if(mode == "R") rq.PropType(JacPropMatType::PopulationRatio);
    else throw std::runtime_error("Mode must be either \"Tv\" or \"R\".  See function description");
    
    rq.MainTag( NLTE_MAINTAG );
    rq.QuantumIdentity(energy_level_identity);
    rq.Perturbation(dx);
    rq.Grids( grids );
    rq.Analytical(1);
    rq.SubSubtag(PROPMAT_SUBSUBTAG);
    
    // Add it to the *jacobian_quantities*
    jq.push_back( rq );
    
    out3 << "  Calculations done by propagation matrix expressions.\n"; 
    
    jacobian_agenda.append( "jacobianCalcDoNothing", TokVal() );
} 


void jacobianAddNLTEs(
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
  const ArrayOfQuantumIdentifier&    energy_level_identities,
  const Numeric&              dx,
  const String&               mode,
  const Verbosity&            verbosity )
{
  for(const auto& qi : energy_level_identities)
    jacobianAddNLTE(ws, jq, jacobian_agenda, atmosphere_dim, p_grid,lat_grid, lon_grid, rq_p_grid, rq_lat_grid, rq_lon_grid, qi, dx, mode, verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddSpecialSpecies(
  Workspace&,
  ArrayOfRetrievalQuantity&   jq,
  Agenda&                     jacobian_agenda,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Vector&                     rq_p_grid,
  const Vector&                     rq_lat_grid,
  const Vector&                     rq_lon_grid,
  const String&                     species,
  const Verbosity&                  verbosity )
{
  CREATE_OUT2;
  CREATE_OUT3;
  
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
  
  // Create the new retrieval quantity
  RetrievalQuantity rq;
  rq.Grids( grids );
  rq.Analytical( 1 );
  rq.SubSubtag(PROPMAT_SUBSUBTAG);
  
  // Make sure modes are valid and complain if they are repeated
  if( species == "electrons" )
  {
    for( Index it=0; it<jq.nelem(); it++ )
    {
      if( jq[it].MainTag() == ELECTRONS_MAINTAG )
      {
        ostringstream os;
        os << "Electrons are already included in *jacobian_quantities*.";
        throw std::runtime_error(os.str());
      }
    }
    rq.MainTag( ELECTRONS_MAINTAG );
    rq.PropType(JacPropMatType::Electrons);
  }
  else if( species == "particulates" )
  {
    for( Index it=0; it<jq.nelem(); it++ )
    {
      if( jq[it].MainTag() == PARTICULATES_MAINTAG )
      {
        ostringstream os;
        os << "Particulates are already included in *jacobian_quantities*.";
        throw std::runtime_error(os.str());
      }
    }
    rq.MainTag( PARTICULATES_MAINTAG );
    rq.PropType(JacPropMatType::Particulates);
  }
  else
  {
    ostringstream os;
    os << "Unknown special species jacobian: \""  << species <<
          "\"\nPlease see *jacobianAddSpecialSpecies* for viable options.";
    throw std::runtime_error(os.str());
  }
  
  // Add it to the *jacobian_quantities*
  jq.push_back( rq );
  
  jacobian_agenda.append( "jacobianCalcDoNothing", TokVal() );
}                    




//----------------------------------------------------------------------------
// Adjustments and transformations
//----------------------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAdjustAndTransform(
        Matrix&                    jacobian,
  const ArrayOfRetrievalQuantity&  jacobian_quantities,
  const Vector&                    x,
  const Verbosity& )
{
  // For flexibility inside inversion_iteration_agenda, we should accept empty
  // Jacobian
  if( jacobian.empty() )
    { return; }

  // Adjustments
  //
  // Unfortunately the adjustment requires both range indices and the
  // untransformed x, which makes things a bit messy
  bool vars_init = false;
  ArrayOfArrayOfIndex jis0;
  Vector x0;
  //
  for( Index q=0; q<jacobian_quantities.nelem(); q++ )
    {
      if( jacobian_quantities[q].MainTag() == ABSSPECIES_MAINTAG  &&
          jacobian_quantities[q].Mode()    == "rel")
        {
          if( !vars_init )
            {
              bool any_affine;
              jac_ranges_indices( jis0, any_affine, jacobian_quantities, true );
              x0 = x;
              transform_x_back( x0, jacobian_quantities );
              vars_init = true;
            }
          for( Index i=jis0[q][0]; i<=jis0[q][1]; i++ )
            {
              if( x[i] != 1 )
                { jacobian(joker,i) /= x[i]; }
            }
        }
    }

  // Transformations
  transform_jacobian( jacobian, x, jacobian_quantities );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianSetAffineTransformation(
    ArrayOfRetrievalQuantity& jqs,
    const Matrix& transformation_matrix,
    const Vector& offset_vector,
    const Verbosity& /*v*/
    )
{
    if (jqs.empty()) {
      runtime_error("Jacobian quantities is empty, so there is nothing to add the "
                    "transformation to.");
    }

    Index nelem = jqs.back().Grids().nelem();

    if (!(nelem == transformation_matrix.nrows())) {
      runtime_error("Dimension of transformation matrix incompatible with retrieval grids.");
    }
    if (!(nelem == offset_vector.nelem())) {
      runtime_error("Dimension of offset vector incompatible with retrieval grids.");
    }

    jqs.back().SetTransformationMatrix(transformation_matrix);
    jqs.back().SetOffsetVector(offset_vector);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianSetFuncTransformation(
    ArrayOfRetrievalQuantity& jqs,
    const String& transformation_func,
    const Numeric& z_min,
    const Numeric& z_max,
    const Verbosity& /*v*/
    )
{
  if( jqs.empty() )
    throw runtime_error( "Jacobian quantities is empty, so there is nothing to add the "
                   "transformation to." );

  if( transformation_func == "none" )
    {
      jqs.back().SetTransformationFunc( "" );
      return;
    }

  Vector pars;

  if( transformation_func == "atanh" )
    {
      if( z_max <= z_min )
        throw runtime_error(
          "For option atanh, the GIN *z_max* must be set and be > z_min." );
      pars.resize(2);
      pars[0] = z_min;
      pars[1] = z_max;
    }
  else if( transformation_func == "log"  ||  transformation_func == "log10"  )
    {
      pars.resize(1);
      pars[0] = z_min;
    }
  else
    {
      ostringstream os;
      os << "Valid options for *transformation_func* are:\n"
         << "\"none\", \"log\", \"log10\" and \"atanh\"\n"
         << "But found: \"" << transformation_func << "\""; 
      throw runtime_error( os.str() );
    }

  jqs.back().SetTransformationFunc( transformation_func );
  jqs.back().SetTFuncParameters( pars );
}



