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
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "math_funcs.h"
#include "messages.h"
#include "jacobian.h"
#include "rte.h"
#include "absorption.h"
#include "physics_funcs.h"


// Temporary solution

void yCalc(
         Vector&                     y,
   const Agenda&                     ppath_step_agenda,
   const Agenda&                     rte_agenda,
   const Agenda&                     iy_space_agenda,
   const Agenda&                     surface_prop_agenda,
   const Agenda&                     iy_cloudbox_agenda,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Vector&                     lat_grid,
   const Vector&                     lon_grid,
   const Tensor3&                    z_field,
   const Tensor3&                    t_field,
   const Tensor4&                    vmr_field,
   const Matrix&                     r_geoid,
   const Matrix&                     z_surface,
   const Index&                      cloudbox_on, 
   const ArrayOfIndex&               cloudbox_limits,
   const Sparse&                     sensor_response,
   const Matrix&                     sensor_pos,
   const Matrix&                     sensor_los,
   const Vector&                     f_grid,
   const Index&                      stokes_dim,
   const Index&                      antenna_dim,
   const Vector&                     mblock_za_grid,
   const Vector&                     mblock_aa_grid )
{
  Vector               y_f, y_za, y_aa;
  ArrayOfIndex         y_pol;
  const Index          n = sensor_response.nrows();
  const Vector         sensor_response_f(n);
  const ArrayOfIndex   sensor_response_pol(n);
  const Vector         sensor_response_za(n);
  const Vector         sensor_response_aa(n);

  RteCalcNoJacobian( y, y_f, y_pol, y_za, y_aa, 
                     ppath_step_agenda, rte_agenda,
                     iy_space_agenda, surface_prop_agenda, iy_cloudbox_agenda, 
                     atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, 
                     t_field, vmr_field, 
                     r_geoid, z_surface, cloudbox_on, cloudbox_limits, 
                     sensor_response, sensor_response_f,
               sensor_response_pol, sensor_response_za, sensor_response_aa,
                     sensor_pos, sensor_los, f_grid, 
                     stokes_dim, antenna_dim, mblock_za_grid, mblock_aa_grid, 
                     "1" );
}

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddAbsSpecies(// WS Output:
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
  // Check that the jacobian matrix is empty. Otherwise it is either
  // not initialised or it is closed.
  if (jac.nrows()!=0 && jac.ncols()!=0)
  {
    ostringstream os;
    os << "The Jacobian matrix is not initialised correctly or closed.\n"
       << "New retrieval quantities can not be added at this point.";
    throw runtime_error(os.str());
  }
  
  // Check retrieval grids, here we just check the length of the grids
  // vs. the atmosphere dimension
  ArrayOfVector grids(atmosphere_dim);
  {
  ostringstream os;
  if (!check_retrieval_grids( grids, os, p_grid, lat_grid, lon_grid,
        rq_p_grid, rq_lat_grid, rq_lon_grid,
        // FIXMEOLE: These strings have to replaced later with the proper
        //           names from the WSM documentation in methods.cc
        "rq_p_grid", "rq_lat_grid", "rq_lon_grid", atmosphere_dim))
    throw runtime_error(os.str());
  }
  
  // Check that method is either "analytic" or "perturbation"
  bool analytical;
  if (method=="perturbation")
    { analytical = 0; }
  else if (method=="analytical")
    { analytical = 1; }
  else
    {
      ostringstream os;
      os << "The method for absorption species retrieval can only be "
         << "\"analytical\"\n or \"perturbation\".";
      throw runtime_error(os.str());
    }
  
  // Check that mode is either "vmr", "nd" or "rel"
  if( mode!="vmr" && mode!="nd" && mode!="rel" )
    {
      throw runtime_error(
                "The retrieval mode can only be \"vmr\", \"nd\" or \"rel\"." );
    }

  // Check that this species is not already included in the jacobian.
  for (Index it=0; it<jq.nelem(); it++)
  {
    if (jq[it].MainTag()=="Abs. species" && jq[it].Subtag()==species)
    {
      ostringstream os;
      os << "The gas species:\n" << species << "\nis already included in "
         << "*jacobian_quantities*.";
      throw runtime_error(os.str());
    }
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq = RetrievalQuantity();
  rq.MainTag("Abs. species");
  rq.Subtag(species);
  rq.Mode(mode);
  rq.Analytical(analytical);
  rq.Perturbation(dx);
  rq.Grids(grids);

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);
  
  // Add gas species method to the jacobian agenda, if perturbation
  if( analytical )
    {
      out3 << "  Calculations done by analytical expression.\n"; 
    }
  else
    {
      out2 << "  Adding absorption species: " << species 
           << " to *jacobian_quantities*\n" << "  and *jacobian_agenda*\n";
      out3 << "  Calculations done by perturbation, size " << dx 
           << " " << mode << ".\n"; 
      String methodname = "jacobianCalcAbsSpecies";
      String kwv = species;
      jacobian_agenda.append (methodname, kwv);
    }
}                    


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddParticle(// WS Output:
                         ArrayOfRetrievalQuantity& jq,
                         Agenda&                   jacobian_agenda,
                         // WS Input:
                         const Matrix&             jac,
                         const Index&              atmosphere_dim,
                         const Vector&             p_grid,
                         const Vector&             lat_grid,
                         const Vector&             lon_grid,
                         const Tensor4&            pnd_field,
                         const Tensor5&            pnd_perturb,
                         const ArrayOfIndex&       cloudbox_limits,
                         // WS Generic Input:
                         const Vector&             rq_p_grid,
                         const Vector&             rq_lat_grid,
                         const Vector&             rq_lon_grid)
{
  throw runtime_error("Particle jacobians not yet handled correctly.");

  // Check that the jacobian matrix is empty. Otherwise it is either
  // not initialised or it is closed.
  if( jac.nrows()!=0 && jac.ncols()!=0 )
    {
      ostringstream os;
      os << "The Jacobian matrix is not initialised correctly or closed.\n"
         << "New retrieval quantities can not be added at this point.";
      throw runtime_error(os.str());
    }
  
  // Check that pnd_perturb is consistent with pnd_field
  if( pnd_perturb.nbooks()!=pnd_field.nbooks() ||
      pnd_perturb.npages()!=pnd_field.npages() ||
      pnd_perturb.nrows()!=pnd_field.nrows() ||
      pnd_perturb.ncols()!=pnd_field.ncols() )
    {
      ostringstream os;
      os << "The perturbation field *pnd_field_perturb* is not consistent with"
         << "*pnd_field*,\none or several dimensions do not match.";
      throw runtime_error(os.str());
    }
  
  // Check that particles are not already included in the jacobian.
  for( Index it=0; it<jq.nelem(); it++ )
    {
      if( jq[it].MainTag()=="Particles" )
        {
          ostringstream os;
          os << "The particles number densities are already included in "
             << "*jacobian_quantities*.";
          throw runtime_error(os.str());
        }
    }
  
  // Particle Jacobian only defined for 1D and 3D atmosphere. check the 
  // retrieval grids, here we just check the length of the grids vs. the 
  // atmosphere dimension
  if (atmosphere_dim==2) 
  {
    ostringstream os;
    os << "Atmosphere dimension not equal to 1 or 3. " 
       << "Jacobians for particle number\n"
       << "density only available for 1D and 3D atmosphere.";
    throw runtime_error(os.str());
  }

  ArrayOfVector grids(atmosphere_dim);
  // The retrieval grids should only consists of gridpoints within
  // the cloudbox. Setup local atmospheric fields inside the cloudbox
  {
    Vector p_cbox = p_grid;
    Vector lat_cbox = lat_grid;
    Vector lon_cbox = lon_grid;
    switch (atmosphere_dim)
      {
      case 3:
        {
          lon_cbox = lon_grid[Range(cloudbox_limits[4], 
                                    cloudbox_limits[5]-cloudbox_limits[4]+1)];
        }
      case 2:
        {
          lat_cbox = lat_grid[Range(cloudbox_limits[2], 
                                    cloudbox_limits[3]-cloudbox_limits[2]+1)];
        }    
      case 1:
        {  
          p_cbox = p_grid[Range(cloudbox_limits[0], 
                                cloudbox_limits[1]-cloudbox_limits[0]+1)];
        }
      }
    ostringstream os;
    if( !check_retrieval_grids( grids, os, p_cbox, lat_cbox, lon_cbox,
                                rq_p_grid, rq_lat_grid, rq_lon_grid, 
        // FIXMEOLE: These strings have to replaced later with the proper
        //           names from the WSM documentation in methods.cc
          "rq_p_grid", "rq_lat_grid", "rq_lon_grid", atmosphere_dim ))
      throw runtime_error(os.str());
  }

  // Common part for all particle variables
  RetrievalQuantity rq = RetrievalQuantity();
  rq.MainTag("Particles");
  rq.Grids(grids);
  rq.Analytical(0);
  rq.Perturbation(-999.999);
  rq.Mode("Fields *mode* and *perturbation* are not defined");

  // Set info for each particle variable
  for( Index ipt=0; ipt<pnd_perturb.nshelves(); ipt++ )
    {
      out2 << "  Adding particle variable " << ipt +1 
           << " to *jacobian_quantities / agenda*.\n";

      ostringstream os;
      os << "Variable " << ipt+1;
      rq.Subtag(os.str());
      
      jq.push_back(rq);
    }

  // Add gas species method to the jacobian agenda
  String methodname = "jacobianCalcParticle";
  String kwv = "";
  jacobian_agenda.append (methodname, kwv);
}                    


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddPointing(// WS Output:
                         ArrayOfRetrievalQuantity&  jq,
                         Agenda&                    jacobian_agenda,
                         // WS Input:
                         const Matrix&              jac,
                         const Matrix&              sensor_pos,
                         const Vector&              sensor_time,
                         // Control Parameters:
                         const Numeric&             dza,
                         const Index&               poly_order)
{
  // Check that the jacobian matrix is empty. Otherwise it is either
  // not initialised or it is closed.
  if (jac.nrows()!=0 && jac.ncols()!=0)
  {
    ostringstream os;
    os << "The Jacobian matrix is not initialised correctly or closed.\n"
       << "New retrieval quantities can not be added at this point.";
    throw runtime_error(os.str());
  }
  
  // Check that poly_order is -1 or positive
  if (poly_order<-1)
    throw runtime_error(
      "The polynomial order has to be positive or -1 for gitter.");
    
  // Define subtag here to easily expand function.
  String subtag="za offset";
  String mode="abs";

  // Check that this type of pointing is not already included in the
  // jacobian.
  for (Index it=0; it<jq.nelem(); it++)
  {
    if (jq[it].MainTag()=="Pointing" && jq[it].Subtag()==subtag)
    {
      ostringstream os;
      os << "A zenith angle pointing offset is already included in\n"
         << "*jacobian_quantities*.";
      throw runtime_error(os.str());
    }
  }

  // Check that sensor_time is consistent with sensor_pos
  if (sensor_time.nelem()!=sensor_pos.nrows())
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
  rq.MainTag("Pointing");
  rq.Subtag(subtag);
  rq.Mode(mode);
  rq.Analytical(0);
  rq.Perturbation(dza);
  // To store the value or the polynomial order, create a vector with length
  // poly_order+1, in case of gitter set the size of the grid vector to be the
  // number of measurement blocks, all elements set to -1.
  Vector grid(0,poly_order+1,1);
  if (poly_order==-1)
  {
    grid.resize(sensor_pos.nrows());
    grid = -1.0;
  }
  ArrayOfVector grids(1, grid);
  rq.Grids(grids);

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  // Add pointing method to the jacobian agenda
  String methodname = "jacobianCalcPointing";
  String kwv = "";
  jacobian_agenda.append (methodname, kwv);
  
  out2 << "  Adding zenith angle pointing offset to *jacobian_quantities*\n"
       << "  and *jacobian_agenda* with perturbation size " << dza << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianAddTemperature(// WS Output:
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
                    const String&             hse,
                    const String&             method,
                    const Numeric&            dx)
{
  throw runtime_error("Temperature jacobians not yet handled correctly.");


  // Check that the jacobian matrix is empty. Otherwise it is either
  // not initialised or it is closed.
  if (jac.nrows()!=0 && jac.ncols()!=0)
  {
    ostringstream os;
    os << "The Jacobian matrix is not initialised correctly or closed.\n"
       << "New retrieval quantities can not be added at this point.";
    throw runtime_error(os.str());
  }
  
  // Check retrieval grids, here we just check the length of the grids
  // vs. the atmosphere dimension
  ArrayOfVector grids(atmosphere_dim);
  {
  ostringstream os;
  if (!check_retrieval_grids( grids, os, p_grid, lat_grid, lon_grid,
        rq_p_grid, rq_lat_grid, rq_lon_grid,
        // FIXMEOLE: These strings have to replaced later with the proper
        //           names from the WSM documentation in methods.cc
        "rq_p_grid", "rq_lat_grid", "rq_lon_grid", atmosphere_dim))
    throw runtime_error(os.str());
  }
  
  // Set subtag to "HSE off"
  String subtag;
  if (hse=="on")
  {
    subtag = "HSE on";
    //FIXME: when implemented, remove this
    ostringstream os;
    os << "Temperature jacobian with HSE on is not implemented yet.";
    throw runtime_error(os.str());
  }
  else if (hse=="off")
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
    
  // Check that method is either "analytic" or "perturbation"
  bool analytical;
  if (method=="perturbation")
    {
      analytical = 0;
    }
  else if (method=="analytic")
    {
    //FIXME: Only "perturbation" implemented so far
    throw runtime_error(
      "Only perturbation method implemented for temperature retrieval.");
    
      analytical = 1;
    }
  else
  {
    ostringstream os;
    os << "The method for temperature retrieval can only be \"analytic\"\n"
       << "or \"perturbation\".";
    throw runtime_error(os.str());
  }
  
  // Check that temperature is not already included in the jacobian.
  // We only check the main tag.
  for (Index it=0; it<jq.nelem(); it++)
  {
    if (jq[it].MainTag()=="Temperature")
    {
      ostringstream os;
      os << "Temperature is already included in *jacobian_quantities*.";
      throw runtime_error(os.str());
    }
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq = RetrievalQuantity();
  rq.MainTag("Temperature");
  rq.Subtag(subtag);
  rq.Mode("abs");
  rq.Analytical(analytical);
  rq.Perturbation(dx);
  rq.Grids(grids);

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);
  
  // Add pointing method to the jacobian agenda
  String methodname = "jacobianCalcTemperature";
  String kwv = "";
  jacobian_agenda.append (methodname, kwv);
  
  out2 << "  Adding temperature to *jacobian_quantities* and "
       << "*jacobian_agenda*.\n";
  if( analytical ) 
  {
    out3 << "  Calculations done by analytical expression.\n"; 
  }
  else
  { 
    out3 << "  Calculations done by perturbations, of size " << dx << ".\n"; 
  }
}                    


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalc(// WS Output:
                        Matrix&                    jacobian,
                  // WS Input:
                  const Agenda&                    jacobian_agenda,
                  const ArrayOfArrayOfIndex&       jacobian_indices )
{
  // Check that the jacobian has been initialised. This is covered by the
  // next test, but gives a better output for this specific case.
  if (jacobian.ncols()==0)
  {
    ostringstream os;
    os << "The Jacobian matrix has not been properly initialised.\n"
       << "The WSM jacobianClose has to be called prior to jacobianCalc.\n";
    throw runtime_error(os.str());
  }

  // Check that *jacobian_indices* and *jacobian* are consistent
  ArrayOfIndex last_ind = jacobian_indices[jacobian_indices.nelem()-1];
  if (jacobian.ncols()-1!=last_ind[1])
  {
    ostringstream os;
    os << "There are more retrieval quantities in *jacobian_quantities*\n"
       << "than in *jacobian*. After calling jacobianClose no more\n"
       << "quantities can be added.";
    throw runtime_error(os.str());
  }

  // Output message
  out2 << "  Calculating *jacobian*.\n";
  
  // Run jacobian_agenda
  jacobian_agendaExecute( jacobian, jacobian_agenda, false );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcAbsSpecies(
     // WS Output:
           Matrix&                   jacobian,
     // WS Input:
     const Vector&                   y,
     const ArrayOfRetrievalQuantity& jq,
     const ArrayOfArrayOfIndex&      jacobian_indices,
     const ArrayOfArrayOfSpeciesTag& abs_species,
     const Agenda&                   ppath_step_agenda,
     const Agenda&                   rte_agenda,
     const Agenda&                   iy_space_agenda,
     const Agenda&                   surface_prop_agenda,
     const Agenda&                   iy_cloudbox_agenda,
     const Index&                    atmosphere_dim,
     const Vector&                   p_grid,
     const Vector&                   lat_grid,
     const Vector&                   lon_grid,
     const Tensor3&                  z_field,
     const Tensor3&                  t_field,
     const Tensor4&                  vmr_field,
     const Matrix&                   r_geoid,
     const Matrix&                   z_surface,
     const Index&                    cloudbox_on,
     const ArrayOfIndex&             cloudbox_limits,
     const Sparse&                   sensor_response,
     const Matrix&                   sensor_pos,
     const Matrix&                   sensor_los,
     const Vector&                   f_grid,
     const Index&                    stokes_dim,
     const Index&                    antenna_dim,
     const Vector&                   mblock_za_grid,
     const Vector&                   mblock_aa_grid,
     // Control Parameters:
     const String&                   species)
{
  // Set some useful (and needed) variables. 
  Index n_jq = jq.nelem();
  RetrievalQuantity rq;
  ArrayOfIndex ji;
  Index it, method;
  
  // Find the retrieval quantity related to this method, i.e. Abs. species -
  // species. This works since the combined MainTag and Subtag is individual.
  bool check_rq = false;
  for( Index n=0; n<n_jq && !check_rq; n++ )
    {
      if (jq[n].MainTag()=="Abs. species" && jq[n].Subtag()==species)
        {
          check_rq = true;
          rq = jq[n];
          ji = jacobian_indices[n];
        }
    }
  if( !check_rq )
    {
      ostringstream os;
      os << "There is no gas species retrieval quantities defined for:\n"
         << species;
      throw runtime_error(os.str());
    }

  // Should not be analytical
  assert( !rq.Analytical() );
   
  
  // Store the start JacobianIndices and the Grids for this quantity
  it = ji[0];
  ArrayOfVector jg = rq.Grids();

  // Check if a relative pertubation is used or not, this information is needed
  // by the methods 'perturbation_field_?d'.
  // Note: both 'vmr' and 'nd' are absolute perturbations
  if (rq.Mode()=="rel")
    method = 0;
  else 
    method = 1;
  
  // For each atmospheric dimension option calculate a ArrayOfGridPos, these
  // are the base functions for interpolating the perturbations into the
  // atmospheric grids.
  ArrayOfGridPos p_gp,lat_gp,lon_gp;
  Index j_p = jg[0].nelem();
  Index j_lat = 1;
  Index j_lon = 1;
  get_perturbation_gridpos(p_gp, p_grid, jg[0], true);
  if (atmosphere_dim>=2) 
  {
    j_lat = jg[1].nelem();
    get_perturbation_gridpos( lat_gp, lat_grid, jg[1], false);
    if (atmosphere_dim==3) 
    {
      j_lon = jg[2].nelem();
      get_perturbation_gridpos( lon_gp, lon_grid, jg[2], false);
    }
  }

  // Give verbose output
  out2 << "  Calculating retrieval quantity " << rq << "\n";
  
  // Find VMR field for these species. 
  ArrayOfSpeciesTag tags;
  array_species_tag_from_string( tags, species );
  Index si = chk_contains( "species", abs_species, tags );

  // Variables for vmr field perturbation unit conversion
  Tensor3 nd_field(t_field.npages(),t_field.nrows(),t_field.ncols(), 1.0);
  if (rq.Mode()=="nd")
    calc_nd_field(nd_field, p_grid, t_field);
  
  // Vector for perturbed measurement vector
  Vector yp;;
    
  // Loop through the retrieval grid and calculate perturbation effect
  for (Index lon_it=0; lon_it<j_lon; lon_it++)
  {
    for (Index lat_it=0; lat_it<j_lat; lat_it++)
    {
      for (Index p_it=0; p_it<j_p; p_it++)
      {
        // Here we calculate the ranges of the perturbation. We want the
        // perturbation to continue outside the atmospheric grids for the
        // edge values.
        Range p_range = Range(0,0);
        Range lat_range = Range(0,0);
        Range lon_range = Range(0,0);
        get_perturbation_range( p_range, p_it, j_p);
        if (atmosphere_dim>=2)
        {
          get_perturbation_range( lat_range, lat_it, j_lat);
          if (atmosphere_dim==3)
          {
            get_perturbation_range( lon_range, lon_it, j_lon);
          }
        }

        // Create VMR field to perturb
        Tensor4 vmr_p = vmr_field;
                              
        // If perturbation given in ND convert the vmr-field to ND before
        // the perturbation is added          
        if (rq.Mode()=="nd")
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
              p_gp, jg[0].nelem()+2, p_range, rq.Perturbation(), method);
            break;
          }
          case 2:
          {
            // Here we perturb a matrix
            perturbation_field_2d( vmr_p(si,joker,joker,lon_it),
              p_gp, lat_gp, jg[0].nelem()+2, jg[1].nelem()+2, p_range, 
              lat_range, rq.Perturbation(), method);
            break;
          }    
          case 3:
          {  
            // Here we need to perturb a tensor3
            perturbation_field_3d( vmr_p(si,joker,joker,joker), 
              p_gp, lat_gp, lon_gp, jg[0].nelem()+2, jg[1].nelem()+2, 
              jg[2].nelem()+2, p_range, lat_range, lon_range, 
              rq.Perturbation(), method);
            break;
          }
        }

        // If perturbation given in ND convert back to VMR          
        if (rq.Mode()=="nd")
          vmr_p(si,joker,joker,joker) /= nd_field;
        
        // Calculate the perturbed spectrum  
        out2 << "  Calculating perturbed spectra no. " << it+1 << " of "
             << ji[1]+1 << "\n";
        yCalc( yp, ppath_step_agenda, rte_agenda,
                 iy_space_agenda, surface_prop_agenda, iy_cloudbox_agenda, 
                 atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, 
                 t_field, vmr_p,
                 r_geoid, z_surface, cloudbox_on, cloudbox_limits, 
                 sensor_response, sensor_pos, sensor_los, f_grid, 
                 stokes_dim, antenna_dim, mblock_za_grid, mblock_aa_grid);
    
        // Add dy/dx as column in jacobian
        for (Index y_it=0; y_it<y.nelem(); y_it++)
          {
            jacobian(y_it,it) = (yp[y_it]-y[y_it])/rq.Perturbation();
          }

        it++;
      }
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcParticle(
     // WS Output:
           Matrix&                     jacobian,
     // WS Input:
     const Vector&                     y,
     const ArrayOfRetrievalQuantity&   jq,
     const ArrayOfArrayOfIndex&        jacobian_indices,
     const Tensor5&                    pnd_field_perturb,
     const Agenda&                     jacobian_particle_update_agenda,
     const Agenda&                     ppath_step_agenda,
     const Agenda&                     rte_agenda,
     const Agenda&                     iy_space_agenda,
     const Agenda&                     surface_prop_agenda,
     const Agenda&                     iy_cloudbox_agenda,
     const Index&                      atmosphere_dim,
     const Vector&                     p_grid,
     const Vector&                     lat_grid,
     const Vector&                     lon_grid,
     const Tensor3&                    z_field,
     const Tensor3&                    t_field,
     const Tensor4&                    vmr_field,
     const Matrix&                     r_geoid,
     const Matrix&                     z_surface,
     const Index&                      cloudbox_on,
     const ArrayOfIndex&               cloudbox_limits,
     const Tensor4&                    pnd_field,
     const Sparse&                     sensor_response,
     const Matrix&                     sensor_pos,
     const Matrix&                     sensor_los,
     const Vector&                     f_grid,
     const Index&                      stokes_dim,
     const Index&                      antenna_dim,
     const Vector&                     mblock_za_grid,
     const Vector&                     mblock_aa_grid )
{
  // Set some useful (and needed) variables. 
  Index n_jq = jq.nelem();
  RetrievalQuantity rq;
  ArrayOfIndex ji;
  
  // Setup local atmospheric fields inside the cloudbox
  Vector p_cbox = p_grid;
  Vector lat_cbox = lat_grid;
  Vector lon_cbox = lon_grid;
  switch (atmosphere_dim)
    {
    case 3:
      {
        lon_cbox = lon_grid[Range(cloudbox_limits[4], 
                                  cloudbox_limits[5]-cloudbox_limits[4]+1)];
      }
    case 2:
      {
        lat_cbox = lat_grid[Range(cloudbox_limits[2], 
                                  cloudbox_limits[3]-cloudbox_limits[2]+1)];
      }    
    case 1:
      {  
        p_cbox = p_grid[Range(cloudbox_limits[0], 
                              cloudbox_limits[1]-cloudbox_limits[0]+1)];
      }
    }


  // Variables to handle and store perturbations
  Vector yp;
  Tensor4  pnd_p, base_pert = pnd_field;


  // Loop particle variables (indexed by *ipt*, where *ipt* is zero based)
  //
  Index ipt      = -1;
  bool not_ready = true;
  
  while( not_ready )
    {
      // Step *ipt*
      ipt++;

      // Define sub-tag string
      ostringstream os;
      os << "Variable " << ipt+1;

      // Find the retrieval quantity related to this particle type
      //
      bool  found = false;
      //
      for( Index n=0; n<n_jq; n++ )
        {
          if( jq[n].MainTag()=="Particles" && jq[n].Subtag()== os.str() )
            {
              found = true;
              rq = jq[n];
              ji = jacobian_indices[n];
              n  = n_jq;                   // To jump out of for-loop
            }
        }

      // At least one particle type must be found
      assert( !( ipt==0  &&  !found ) );

      // Ready or something to do?
      if( !found )
        { 
          not_ready = false;
        }
      else
        {
          // Counters for report string
          Index   it  = 0;
          Index   nit = ji[1] -ji[0] + 1;
          
          // Counter for column in *jacobian*
          Index   icol = ji[0];

          // Retrieval grid positions
          ArrayOfVector jg = rq.Grids();
          ArrayOfGridPos p_gp, lat_gp, lon_gp;
          Index j_p = jg[0].nelem();
          Index j_lat = 1;
          Index j_lon = 1;
          get_perturbation_gridpos( p_gp, p_cbox, jg[0], true );
          if (atmosphere_dim==3) 
            {
              j_lat = jg[1].nelem();
              get_perturbation_gridpos( lat_gp, lat_cbox, jg[1], false );
              
              j_lon = jg[2].nelem();
              get_perturbation_gridpos( lon_gp, lon_cbox, jg[2], false );
            }

          // Give verbose output
          out1 << "  Calculating retrieval quantity:" << rq << "\n";
  

          // Loop through the retrieval grid and calculate perturbation effect
          for (Index lon_it=0; lon_it<j_lon; lon_it++)
            {
              for (Index lat_it=0; lat_it<j_lat; lat_it++)
                {
                  for (Index p_it=0; p_it<j_p; p_it++)
                    {
                      // Update the perturbation field
                      pnd_p = 
                           pnd_field_perturb( ipt, joker, joker, joker, joker);

                      it++;
                      out1 << "  Calculating perturbed spectra no. " << it
                           << " of " << nit << "\n";

                      // Here we calculate the ranges of the perturbation. 
                      // We want the perturbation to continue outside the 
                      // atmospheric grids for the edge values.
                      Range p_range   = Range(0,0);
                      Range lat_range = Range(0,0);
                      Range lon_range = Range(0,0);
                      get_perturbation_range( p_range, p_it, j_p );
                      if (atmosphere_dim==3)
                        {
                          get_perturbation_range( lat_range, lat_it, j_lat);
                          get_perturbation_range( lon_range, lon_it, j_lon);
                        }
                          
                      // Make empty copy of pnd_pert for base functions
                      base_pert *= 0;
            
                      // Calculate the perturbed field according to atm_dim, 
                      switch (atmosphere_dim)
                        {
                        case 1:
                          {
                            for( Index typ_it=0; typ_it<pnd_field.nbooks(); 
                                                                     typ_it++ )
                              {
                                perturbation_field_1d( 
                                      base_pert(typ_it,joker,lat_it,lon_it),
                                      p_gp, jg[0].nelem()+2, p_range, 1.0, 1 );
                              }
                            break;
                          }
                        case 3:
                          {  
                            for( Index typ_it=0; typ_it<pnd_field.nrows(); 
                                                                     typ_it++ )
                              {
                                perturbation_field_3d( 
                                      base_pert(typ_it,joker,joker,joker),
                                      p_gp, lat_gp, lon_gp, jg[0].nelem()+2, 
                                      jg[1].nelem()+2, jg[2].nelem()+2, 
                                      p_range, lat_range, lon_range, 1.0, 1);
                              }
                            break;
                          }
                        }
          
                      // Now add the weighted perturbation field to the 
                      // reference field and recalculate the scattered field
                      pnd_p *= base_pert;
                      pnd_p += pnd_field;
                      jacobian_particle_update_agendaExecute( pnd_p, 
                                      jacobian_particle_update_agenda, false );
            
                      // Calculate the perturbed spectrum  
                      yCalc( yp, ppath_step_agenda, rte_agenda, 
                             iy_space_agenda, surface_prop_agenda, 
                             iy_cloudbox_agenda, atmosphere_dim,
                             p_grid, lat_grid, lon_grid, z_field, t_field, 
                             vmr_field, r_geoid, z_surface, cloudbox_on, 
                             cloudbox_limits, sensor_response, sensor_pos, 
                             sensor_los, f_grid, stokes_dim, antenna_dim, 
                             mblock_za_grid, mblock_aa_grid);
    
                      // Add dy as column in jacobian. Note that we just return
                      // the difference between the two spectra.
                      for( Index y_it=0; y_it<yp.nelem(); y_it++ )
                        {
                          jacobian(y_it,icol) = yp[y_it]-y[y_it];
                        }

                      // Step *icol*
                      icol++;
                    }
                }
            }
        }
    }
}

                     
/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcPointing(
     // WS Output:
           Matrix&                   jacobian,
     // WS Input:
     const Vector&                   y,
     const ArrayOfRetrievalQuantity& jq,
     const ArrayOfArrayOfIndex&      jacobian_indices,
     const Vector&                   sensor_time,
     const Agenda&                   ppath_step_agenda,
     const Agenda&                   rte_agenda,
     const Agenda&                   iy_space_agenda,
     const Agenda&                   surface_prop_agenda,
     const Agenda&                   iy_cloudbox_agenda,
     const Index&                    atmosphere_dim,
     const Vector&                   p_grid,
     const Vector&                   lat_grid,
     const Vector&                   lon_grid,
     const Tensor3&                  z_field,
     const Tensor3&                  t_field,
     const Tensor4&                  vmr_field,
     const Matrix&                   r_geoid,
     const Matrix&                   z_surface,
     const Index&                    cloudbox_on,
     const ArrayOfIndex&             cloudbox_limits,
     const Sparse&                   sensor_response,
     const Matrix&                   sensor_pos,
     const Matrix&                   sensor_los,
     const Vector&                   f_grid,
     const Index&                    stokes_dim,
     const Index&                    antenna_dim,
     const Vector&                   mblock_za_grid,
     const Vector&                   mblock_aa_grid)
{
  // Set some useful (and needed) variables.  
  //  Index n_los = sensor_los.nrows();
  Matrix sensor_los_pert = sensor_los;
  RetrievalQuantity rq;
  ArrayOfIndex ji;
  bool gitter = false;
  
  // Check that sensor_time is consistent with sensor_pos
  assert(sensor_time.nelem()==sensor_pos.nrows());

  // Find the retrieval quantity related to this method, i.e. Pointing
  // za offset. This works since the combined MainTag and Subtag is individual.
  bool check_rq = false;
  for (Index n=0; n<jq.nelem(); n++)
  {
    if (jq[n].MainTag()=="Pointing" && jq[n].Subtag()=="za offset")
    {
      check_rq = true;
      rq = jq[n];
      ji = jacobian_indices[n];
    }
  }
  if (!check_rq)
  {
    throw runtime_error(
      "There is no pointing offset retrieval quantities defined.\n");
  }

  // This method only handles absolute perturbations
  assert( rq.Mode()=="abs" );
  
  // Check if pointing is of type gitter
  if (rq.Grids()[0][0]==-1)
  {
    gitter = true;
  }
  
  // Calculate the weight vector. We set sensor_time[0] to correspond to
  // -1 and sensor_time[end] to 1.
  chk_if_increasing("sensor_time", sensor_time);
  Vector weight(sensor_time.nelem());
  Numeric t0 = sensor_time[0];
  Numeric t1 = sensor_time[sensor_time.nelem()-1];
  for (Index tt=0; tt<weight.nelem(); tt++)
  {
    weight[tt] = 2*(sensor_time[tt]-t0)/(t1-t0)-1;
  }

  // Give verbose output
  out2 << "  Calculating retrieval quantity " << rq << "\n";

  // Declare variables for perturbed and difference spectra
  Vector yp;
  Vector dydx(y.nelem());
  
  // Add the pointing offset. 
  sensor_los_pert(joker,0) += rq.Perturbation();
     
  // Calculate the perturbed spectrum for the zeroth order polynomial
  yCalc( yp, ppath_step_agenda, rte_agenda,
         iy_space_agenda, surface_prop_agenda, iy_cloudbox_agenda, 
         atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, t_field, 
         vmr_field, r_geoid, 
         z_surface, cloudbox_on, cloudbox_limits, sensor_response, 
         sensor_pos, sensor_los_pert, f_grid, stokes_dim, antenna_dim, 
         mblock_za_grid, mblock_aa_grid);
    
  // Calculate difference in spectrum and divide by perturbation,
  // here we want the whole dy/dx vector so that we can use it later
  dydx = yp;
  dydx -= y;
  dydx /= rq.Perturbation();
    
  // Add the weighted dy/dx as column in jacobian
  Index ny = y.nelem()/sensor_pos.nrows();
  Index it = ji[0];
  Numeric exponent;
  while (it<=ji[1])
  {
    // For gitter the exponent is zero for all columns
    if (!gitter)
    {
      exponent = (Numeric) it-ji[0];
    }
    else
    {  
      exponent = 0.0;
    } 
    out2 << "  Calculating perturbed spectra no. " << it+1 << " of "
         << ji[1]+1 << "\n";
    Index y_it = 0;
    for (Index ns=0; ns<sensor_pos.nrows(); ns++)
    {
      for (Index dummy=0; dummy<ny; dummy++)
      {
        jacobian(y_it,it) = dydx[y_it]*pow(weight[ns], exponent);
        y_it++;
      }
      // If gitter then change column for each row in sensor_pos
      if (gitter)
        it++;
    }
    // If not gitter then change column for each order of polynomial
    if (!gitter)
      it++;
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianCalcTemperature(
     // WS Output:
           Matrix&                   jacobian,
     // WS Input:
     const Vector&                   y,
     const ArrayOfRetrievalQuantity& jq,
     const ArrayOfArrayOfIndex&      jacobian_indices,
     const Agenda&                   ppath_step_agenda,
     const Agenda&                   rte_agenda,
     const Agenda&                   iy_space_agenda,
     const Agenda&                   surface_prop_agenda,
     const Agenda&                   iy_cloudbox_agenda,
     const Index&                    atmosphere_dim,
     const Vector&                   p_grid,
     const Vector&                   lat_grid,
     const Vector&                   lon_grid,
     const Tensor3&                  z_field,
     const Tensor3&                  t_field,
     const Tensor4&                  vmr_field,
     const Matrix&                   r_geoid,
     const Matrix&                   z_surface,
     const Index&                    cloudbox_on,
     const ArrayOfIndex&             cloudbox_limits,
     const Sparse&                   sensor_response,
     const Matrix&                   sensor_pos,
     const Matrix&                   sensor_los,
     const Vector&                   f_grid,
     const Index&                    stokes_dim,
     const Index&                    antenna_dim,
     const Vector&                   mblock_za_grid,
     const Vector&                   mblock_aa_grid)
{
  // Set some useful (and needed) variables. 
  RetrievalQuantity rq;
  ArrayOfIndex ji;
  Index it;
  
  // Find the retrieval quantity related to this method, i.e. Temperature.
  // For temperature only the main tag is checked.
  bool check_rq = false;
  for (Index n=0; n<jq.nelem(); n++)
  {
    if (jq[n].MainTag()=="Temperature")
    {
      check_rq = true;
      rq = jq[n];
      ji = jacobian_indices[n];
    }
  }
  if (!check_rq)
  {
    ostringstream os;
    os << "There is no temperature retrieval quantities defined.\n";
    throw runtime_error(os.str());
  }
  
  // FIXME: Only HSE off is implemented
  assert( rq.Subtag()=="HSE off" );
  
  // Store the start JacobianIndices and the Grids for this quantity
  it = ji[0];
  ArrayOfVector jg = rq.Grids();

  // "Perturbation method". 1 means absolute perturbations
  const Index method = 1;
      
  // For each atmospheric dimension option calculate a ArrayOfGridPos, 
  // these will be used to interpolate a perturbation into the atmospheric 
  // grids.
  ArrayOfGridPos p_gp,lat_gp,lon_gp;
  Index j_p = jg[0].nelem();
  Index j_lat = 1;
  Index j_lon = 1;
  get_perturbation_gridpos(p_gp, p_grid, jg[0], true);
  if (atmosphere_dim>=2) 
  {
    j_lat = jg[1].nelem();
    get_perturbation_gridpos( lat_gp, lat_grid, jg[1], false);
    if (atmosphere_dim==3) 
    {
      j_lon = jg[2].nelem();
      get_perturbation_gridpos( lon_gp, lon_grid, jg[2], false);
    }
  }
  
  // Give verbose output
  out2 << "  Calculating retrieval quantity " << rq << "\n";

  // Loop through the retrieval grid and calculate perturbation effect
  for (Index lon_it=0; lon_it<j_lon; lon_it++)
  {
    for (Index lat_it=0; lat_it<j_lat; lat_it++)
    {
      for (Index p_it=0; p_it<j_p; p_it++)
      {
        // Perturbed spectrum and temperature field
        Vector yp;
        Tensor3 t_p = t_field;
  
        // Here we calculate the ranges of the perturbation. We want the
        // perturbation to continue outside the atmospheric grids for the
        // edge values.
        Range p_range = Range(0,0);
        Range lat_range = Range(0,0);
        Range lon_range = Range(0,0);
        get_perturbation_range( p_range, p_it, j_p);
        if (atmosphere_dim>=2)
        {
          get_perturbation_range( lat_range, lat_it, j_lat);
          if (atmosphere_dim==3)
          {
            get_perturbation_range( lon_range, lon_it, j_lon);
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
              p_gp, jg[0].nelem()+2, p_range, rq.Perturbation(), method);
            break;
          }
          case 2:
          {
            // Here we perturb a matrix
            perturbation_field_2d( t_p(joker,joker,lon_it), 
              p_gp, lat_gp, jg[0].nelem()+2, jg[1].nelem()+2, p_range, lat_range, 
              rq.Perturbation(), method);
            break;
          }    
          case 3:
          {  
            // Here we need to perturb a tensor3
            perturbation_field_3d( t_p(joker,joker,joker), 
              p_gp, lat_gp, lon_gp, jg[0].nelem()+2, jg[1].nelem()+2, 
              jg[2].nelem()+2, p_range, lat_range, lon_range, 
              rq.Perturbation(), method);
            break;
          }
        }
          
        // Calculate the perturbed spectrum  
        out2 << "  Calculating perturbed spectra no. " << it+1 << " of "
             << ji[1]+1 << "\n";

        yCalc( yp, ppath_step_agenda, rte_agenda, iy_space_agenda, 
               surface_prop_agenda, iy_cloudbox_agenda, atmosphere_dim, p_grid,
               lat_grid, lon_grid, z_field, t_p, vmr_field, 
               r_geoid, z_surface, cloudbox_on, cloudbox_limits, 
               sensor_response, sensor_pos, sensor_los, f_grid, 
               stokes_dim, antenna_dim, mblock_za_grid, mblock_aa_grid);
    
        // Add dy/dx as column in jacobian
        for (Index y_it=0; y_it<y.nelem(); y_it++)
          {
            jacobian(y_it,it) = (yp[y_it]-y[y_it])/rq.Perturbation();
          }

        it++;
      }
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianClose(// WS Output:
                   Matrix&                          jacobian,
                   ArrayOfArrayOfIndex&             jacobian_indices,
                   // WS Input:
                   const ArrayOfRetrievalQuantity&  jacobian_quantities,
                   const Matrix&                    sensor_pos,
                   const Sparse&                    sensor_response)
{
  // Check that *jacobian* has been initialised
  if (jacobian.nrows()!=0 && jacobian.ncols()!=0)
    throw runtime_error("The Jacobian matrix has not been initialised.");

  // Make sure that the array is not empty
  if (jacobian_quantities.nelem()==0)
    throw runtime_error(
      "No retrieval quantities has been added to *jacobian_quantities*");

  // Check that sensor_pol and sensor_response has been initialised
  if (sensor_pos.nrows()==0)
  {
    ostringstream os;
    os << "The number of rows in *sensor_pos* is zero, i.e. no measurement\n"
       << "blocks has been defined. This has to be done before calling\n"
       << "jacobianClose.";
    throw runtime_error(os.str());
  }
  if (sensor_response.nrows()==0)
  {
    ostringstream os;
    os << "The sensor has either to be defined or turned off before calling\n"
       << "jacobianClose.";
    throw runtime_error(os.str());
  }

  // Loop over retrieval quantities, set JacobianIndices
  Index nrows = sensor_pos.nrows()*sensor_response.nrows();
  Index ncols = 0;
  for (Index it=0; it<jacobian_quantities.nelem(); it++)
  {
    // Store start jacobian index
    ArrayOfIndex indices(2);
    indices[0] = ncols;

    // Count total number of field points, i.e. product of grid lengths
    Index cols = 1;
    ArrayOfVector grids = jacobian_quantities[it].Grids();
    for (Index jt=0; jt<grids.nelem(); jt++)
      cols *= grids[jt].nelem();
    ncols += cols;

    // Store stop index
    indices[1] = ncols-1;
    jacobian_indices.push_back(indices);
  }
  
  // Resize *jacobian*
  jacobian.resize( nrows, ncols);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianInit(
      Matrix&                    jacobian,
      ArrayOfRetrievalQuantity&  jacobian_quantities,
      ArrayOfArrayOfIndex&       jacobian_indices )
{
  jacobian.resize(0,0);
  jacobian_quantities.resize(0);
  jacobian_indices.resize(0);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianOff(
      Matrix&                    jacobian,
      ArrayOfRetrievalQuantity&  jacobian_quantities,
      ArrayOfArrayOfIndex&       jacobian_indices,
      String&                    jacobian_unit )
{
  jacobianInit( jacobian, jacobian_quantities, jacobian_indices );
  
  jacobian_unit = "-";
}



/* Workspace method: Doxygen documentation will be auto-generated */
void jacobianUnit(
              Matrix&   jacobian,
        const String&   jacobian_unit,
        const String&   y_unit,
        const Vector&   sensor_response_f )
{
  String j_unit = jacobian_unit;
  //
  if ( jacobian_unit == "-" )
    { j_unit = y_unit; }
  
  try
    {
      ybatchUnit( jacobian, j_unit, sensor_response_f );
    }
  catch( runtime_error e ) 
    {
      ostringstream os;
      os << "Unknown option: jacobian_unit = \"" << j_unit << "\"\n" 
         << "Recognised choices are: \"-\", \"1\", \"RJBT\" and \"PlanckBT\"";
      throw runtime_error( os.str() );      
    }
}
