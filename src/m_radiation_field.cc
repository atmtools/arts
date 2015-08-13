/* Copyright (C) 2015
   Richard Larsson
                            
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
   
#include "arts.h"
#include "arts_omp.h"
#include "auto_md.h"
#include "rte.h"


/* Workspace method: Doxygen documentation will be auto-generated */
void radiation_fieldCalcFromiyCalc(Workspace&              ws,
                                   //OUT:
                                   Matrix&                 iy,
                                   GriddedField4&          radiation_field,
                                   //IN:
                                   const Index&            atmfields_checked,
                                   const Index&            atmgeom_checked,
                                   const Vector&           f_grid,
                                   const Tensor3&          t_field,
                                   const Tensor3&          z_field,
                                   const Tensor4&          vmr_field,
                                   const Index&            cloudbox_on,
                                   const Index&            cloudbox_checked,
                                   const Index&            stokes_dim,
                                   const Vector&           rte_pos,
                                   const String&           iy_unit,  
                                   const Agenda&           iy_main_agenda,
                                   //GIN:
                                   const Vector&           za_coords,
                                   const Vector&           aa_coords,
                                   const Verbosity&        verbosity)
{
  extern const Numeric DEG2RAD, PI;
  
  // Test input
  if(za_coords.nelem()<2)
    throw std::runtime_error("za_coords must contain at least two elements.\n");
  if(aa_coords.nelem()<1)
    throw std::runtime_error("aa_coords must contain at least one element.\n");
  if(!(za_coords[0]==0. && za_coords[za_coords.nelem()-1]==180.))
    throw std::runtime_error("First and last coords in za_coords must be 0 and 180, respectively.\n");
  if(aa_coords.nelem()>1)
    if(!(aa_coords[0]==-180. && aa_coords[aa_coords.nelem()-1]==180.))
      throw std::runtime_error("First and last coords in aa_coords must be -180 and 180, respectively, if aa_coords.nelem()>1.\n");
  {
    Numeric tmp=-1.;   
    for(Index ii=0;ii<za_coords.nelem();ii++)
    {
      if(za_coords[ii]<tmp)
        throw std::runtime_error("za_coords must be strictly increasing.\n");
      else
        tmp=za_coords[ii];
    }
    tmp=-181.;  
    for(Index ii=0;ii<aa_coords.nelem();ii++)
    {
      if(aa_coords[ii]<tmp)
        throw std::runtime_error("aa_coords must be strictly increasing.\n");
      else
        tmp=aa_coords[ii];
    }
  }
  if(f_grid.nelem()==0)
    throw std::runtime_error("f_grid must contain at least one value.\n");
  
  // Necessary input dummy
  const ArrayOfString    iy_aux_vars;
  
  // Prepare a stokes vector so grid setting is possible
  Vector stokes_dim_vector(stokes_dim);
  for(Index ii=0;ii<stokes_dim;ii++)
    stokes_dim_vector[ii]=(Numeric)ii;
  
  // Prepare radiation_field
  radiation_field.set_grid(0, za_coords);
  radiation_field.set_grid_name(0, "Zenith Angle");
  radiation_field.set_grid(1, aa_coords);
  radiation_field.set_grid_name(1, "Azimuth Angle");
  radiation_field.set_grid(2, f_grid);
  radiation_field.set_grid_name(2, "Frequency");
  radiation_field.set_grid(3, stokes_dim_vector);
  radiation_field.set_grid_name(3, "Stokes Component");
  radiation_field.data.resize(za_coords.nelem(),aa_coords.nelem(),
                              f_grid.nelem(),stokes_dim);
  radiation_field.set_name("Radiation Field");
  radiation_field.checksize_strict();
  
  // Loop over coords
  for(Index ii=0; ii<za_coords.nelem();ii++)
    for(Index jj=0; jj<aa_coords.nelem();jj++)
    {
      // Special case of unpolarized zenith angles 0 and 180.  Only need one azimuth angle.
      if(jj>0&&(ii==0||ii==za_coords.nelem()-1)&&stokes_dim==1)
        break;
      
      // Necessary output dummys
      ArrayOfTensor4   iy_aux_dummy;
      Ppath            ppath_dummy;
      
      // Setup of local los vector from za_coords.
      Vector rte_los(2);
      rte_los[0]=za_coords[ii];
      rte_los[1]=aa_coords[jj];
      
      // Calculate iy for a particular rte_los.
      iyCalc(ws, iy, iy_aux_dummy, ppath_dummy,
            atmfields_checked, atmgeom_checked, iy_aux_vars,
            f_grid, t_field, z_field, vmr_field, cloudbox_on,
            cloudbox_checked, rte_pos, rte_los, rte_pos,
            iy_unit, iy_main_agenda, verbosity);
      
      // Add to radiation field
      radiation_field.data(ii,jj,joker,joker) = iy;
      
      // Special case of unpolarized zenith angles 0 and 180.  Only need one azimuth angle.
      if(jj>0&&(ii==0||ii==za_coords.nelem()-1)&&stokes_dim==1)
        for(Index kk=1; kk<aa_coords.nelem();kk++)
          radiation_field.data(ii,kk,joker,joker) = iy;
        
    }
    
    // Get iy from radiation_field
    iy.resize(f_grid.nelem(),stokes_dim);
    iy=0;
    
    // Area elements of the spherical observation geometry
    Matrix dA(za_coords.nelem()-1,aa_coords.nelem());
    for(Index ii=0; ii<za_coords.nelem()-1;ii++)
    {
      for(Index jj=0; jj<aa_coords.nelem()-1;jj++)
      {
          dA(ii,jj) = DEG2RAD * (za_coords[ii+1]-za_coords[ii]) * 
            sin(DEG2RAD * (za_coords[ii] + (za_coords[ii+1]-za_coords[ii])/2.)) * 
            DEG2RAD * (aa_coords[jj+1]-aa_coords[jj]);
      }
      if(aa_coords.nelem()==1)
      {
        dA(ii,0) = DEG2RAD * (za_coords[ii+1]-za_coords[ii]) * 
              sin(DEG2RAD * (za_coords[ii] + (za_coords[ii+1]-za_coords[ii])/2.)) * 
              DEG2RAD * (360.);
      }
    }
    
    // Sum the contributions by relative area of the sphere
    for(Index ii=0; ii<za_coords.nelem()-1;ii++)
    {
      if(aa_coords.nelem()>1)
      {
        for(Index jj=0; jj<aa_coords.nelem()-1;jj++)
        {
          Matrix iy_tmp(f_grid.nelem(),stokes_dim);
          iy_tmp  = radiation_field.data(ii,  jj,  joker,joker);
          iy_tmp += radiation_field.data(ii+1,jj,  joker,joker);
          iy_tmp += radiation_field.data(ii,  jj+1,joker,joker);
          iy_tmp += radiation_field.data(ii+1,jj+1,joker,joker);
          iy_tmp *= dA(ii,jj)/4.0 ;
          iy+=iy_tmp;
        }
      }
      else
      {
        Matrix iy_tmp(f_grid.nelem(),stokes_dim);
        iy_tmp  = radiation_field.data(ii,  0,joker,joker);
        iy_tmp += radiation_field.data(ii+1,0,joker,joker);
        iy_tmp *= dA(ii,0)/2.0 ;
        iy+=iy_tmp;
      }
      
    }
    
    // Normalize for surface area of the sphere
    iy /= 4*PI;
}