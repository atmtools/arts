/* Copyright (C) 2002 Claudia Emde <claudia@sat.physik.uni-bremen.de>
                      
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
  \file   cloudbox.cc
  \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
  \date   Thu May  23 10:59:55 2002
  
  \brief  Internal functions for scattering calculations.
*/

#include "cloudbox.h"

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <stdexcept>
#include <cmath>

#include "arts.h"
#include "auto_md.h"
#include "make_vector.h"
#include "array.h"
#include "logic.h"
#include "ppath.h"
#include "physics_funcs.h"
#include "array.h"


//! Get scattered radiance on 1D cloudbox boundary (cubic interpolation)
/* 
   This method interpolates radiances from the cloudbox boundary on a requested 
   point.
   This method is used to get the result of a scattering calculation as 
   radiative background for a clear-sky calculation. 

   \param i_out Outgoing radiance.
   \param scat_i_p Radiance on pressure surfaces.
   \param scat_i_lat Radiance on latitude surfaces.
   \param scat_i_lon Radiance on longitude surfaces.
   \param rte_gp_p Grid position (pressure).
   \param rte_gp_lat Grid position (latitude).
   \param rte_gp_lon Grid position (longitude). 
   \param rte_los Direction of LOS.
   \param cloudbox limits Cloudbox limits. 
   \param stokes_dim Stokes dimension.
   \param scat_za_grid Zenith angle grid.
   \param scat_aa_grid Azimuth angle grid.
   \param f_grid Frequency grid. 
   
   \author Claudia Emde
   \date 2004-03-20
*/
void cloudbox_getOutgoing1D(// Output:
                            MatrixView i_out,
                            //WS Input:
                            ConstTensor7View scat_i_p,
                            const GridPos& rte_gp_p,
                            ConstVectorView rte_los,
                            const ArrayOfIndex& cloudbox_limits,
                            const Index& stokes_dim,
                            ConstVectorView scat_za_grid,
                            ConstVectorView f_grid
                            )
{
  assert ( is_size( scat_i_p,
                    f_grid.nelem(), 2, 1, 1,
                    scat_za_grid.nelem(), 1,
                    stokes_dim ));
  
  //Check, if grid_positions correspond to cloudbox boundary
  if (rte_gp_p.idx != cloudbox_limits[0] &&
      rte_gp_p.idx != cloudbox_limits[1])
    {
      throw runtime_error(
                          "Gridpositions have to be on the boundary of the"
                          "cloudbox defined by *cloudbox_limits*."
                          );
    }
  
  //Define a vector to interpolate the outgoing radiance which is
  //defined on scat_za_grid on the requested zenith angle in
  //*cloudbox_los*.
      //        Vector zenith_angle(1);
      //  zenith_angle[0] = rte_los[0];
      Numeric zenith_angle = rte_los[0];
  
  //Array to store grid positions
  GridPos gp;
  gridpos(gp, scat_za_grid, zenith_angle);
  
  //Matrix to store interpolation weights
  //Matrix itw(scat_za_grid.nelem(),2);
  Vector itw(2);
  interpweights(itw, gp);
  
  for(Index i = 0; i < stokes_dim; i++)
    {
      for(Index f_index = 0; f_index < f_grid.nelem(); f_index ++)
        {
          //This variable holds the radiation for a specified frequency.
          //It is neccessairy because the interpolation is done for 
          //each frequency separately.
          Vector i_out_f(scat_za_grid.nelem());
	     
          //lower boundary
          if(rte_gp_p.idx == cloudbox_limits[0])
            {
              i_out_f = scat_i_p(f_index, 0, 0, 0, 
                                             Range(joker), 0, i);
            }
          //upper boundary
          else if(rte_gp_p.idx == cloudbox_limits[1])
            {
              i_out_f = scat_i_p(f_index, 1, 0, 0,
                                             Range(joker), 0, i);
            }
          //Define vector for the interpolated radiance.
          Numeric i_out_los;
	       
          //Do the interpolation:
          i_out_los = interp(itw, i_out_f, gp);
	       
          //Put the value into the matrix:
          i_out(f_index, i) = i_out_los;
        }//end frequency loop
    }//end stokes_dim loop
}

//! Get scattered radiance on 1D cloudbox boundary (cubic interpolation)
/* 
   This method interpolates radiances from the cloudbox boundary on a requested 
   point. For the zenith angle dimension cubic interpolation is used. This is usually 
   more accurate, especially for non-limb cases. 
   This method is used to get the result of a scattering calculation as 
   radiative background for a clear-sky calculation. 

   \param i_out Outgoing radiance.
   \param scat_i_p Radiance on pressure surfaces.
   \param scat_i_lat Radiance on latitude surfaces.
   \param scat_i_lon Radiance on longitude surfaces.
   \param rte_gp_p Grid position (pressure).
   \param rte_gp_lat Grid position (latitude).
   \param rte_gp_lon Grid position (longitude). 
   \param rte_los Direction of LOS.
   \param cloudbox limits Cloudbox limits. 
   \param stokes_dim Stokes dimension.
   \param scat_za_grid Zenith angle grid.
   \param scat_aa_grid Azimuth angle grid.
   \param f_grid Frequency grid. 
   
   \author Claudia Emde
   \date 2004-03-20
*/
void cloudbox_getOutgoingCubic1D(// Output:
                            MatrixView i_out,
                            //WS Input:
                            ConstTensor7View scat_i_p,
                            const GridPos& rte_gp_p,
                            ConstVectorView rte_los,
                            const ArrayOfIndex& cloudbox_limits,
                            const Index& stokes_dim,
                            ConstVectorView scat_za_grid,
                            ConstVectorView f_grid
                            )

{
  assert ( is_size( scat_i_p,
                    f_grid.nelem(), 2, 1, 1,
                    scat_za_grid.nelem(), 1,
                    stokes_dim ));
  
  //Check, if grid_positions correspond to cloudbox boundary
  if (rte_gp_p.idx != cloudbox_limits[0] &&
      rte_gp_p.idx != cloudbox_limits[1])
    {
      throw runtime_error(
                          "Gridpositions have to be on the boundary of the"
                          "cloudbox defined by *cloudbox_limits*."
                          );
    }
  
  //Define a vector to interpolate the outgoing radiance which is
  //defined on scat_za_grid on the requested zenith angle in
  //*cloudbox_los*.
      //        Vector zenith_angle(1);
      //  zenith_angle[0] = rte_los[0];
      Numeric zenith_angle = rte_los[0];
  
  //Array to store grid positions
  GridPos gp;
  gridpos(gp, scat_za_grid, zenith_angle);
  
  //Matrix to store interpolation weights
  //Matrix itw(scat_za_grid.nelem(),2);
  Vector itw(2);
  interpweights(itw, gp);
  
  for(Index i = 0; i < stokes_dim; i++)
    {
      for(Index f_index = 0; f_index < f_grid.nelem(); f_index ++)
        {
          //This variable holds the radiation for a specified frequency.
          //It is neccessairy because the interpolation is done for 
          //each frequency separately.
          Vector i_out_f(scat_za_grid.nelem());
	     
          //lower boundary
          if(rte_gp_p.idx == cloudbox_limits[0])
            {
              i_out_f = scat_i_p(f_index, 0, 0, 0, 
                                             Range(joker), 0, i);
            }
          //upper boundary
          else if(rte_gp_p.idx == cloudbox_limits[1])
            {
              i_out_f = scat_i_p(f_index, 1, 0, 0,
                                             Range(joker), 0, i);
            }
          //Define vector for the interpolated radiance.
          Numeric i_out_los;
	       
          //Do the interpolation:
          i_out_los = interp_cubic(scat_za_grid, i_out_f, zenith_angle, gp);

          //Put the value into the matrix:
          i_out(f_index, i) = i_out_los;
        }//end frequency loop
    }//end stokes_dim loop
}


//! Get scattered radiance on 3D cloudbox boundary.
/* 
   This method interpolates radiances from the cloudbox boundary on a requested 
   point. This method is used to get the result of a scattering calculation as 
   radiative background for a clear-sky calculation. 

   \param i_out Outgoing radiance.
   \param scat_i_p Radiance on pressure surfaces.
   \param scat_i_lat Radiance on latitude surfaces.
   \param scat_i_lon Radiance on longitude surfaces.
   \param rte_gp_p Grid position (pressure).
   \param rte_gp_lat Grid position (latitude).
   \param rte_gp_lon Grid position (longitude). 
   \param rte_los Direction of LOS.
   \param cloudbox limits Cloudbox limits. 
   \param stokes_dim Stokes dimension.
   \param scat_za_grid Zenith angle grid.
   \param scat_aa_grid Azimuth angle grid.
   \param f_grid Frequency grid. 
   
   \author Claudia Emde
   \date 2004-03-20
*/
void cloudbox_getOutgoing3D(// Output:
                            MatrixView   i_out,
                            // Input:
                            ConstTensor7View scat_i_p,
                            ConstTensor7View scat_i_lat,
                            ConstTensor7View scat_i_lon,
                            const GridPos& rte_gp_p,
                            const GridPos& rte_gp_lat,
                            const GridPos& rte_gp_lon,
                            ConstVectorView rte_los,
                            const ArrayOfIndex& cloudbox_limits,
                            const Index& stokes_dim,
                            ConstVectorView scat_za_grid,
                            ConstVectorView scat_aa_grid,
                            ConstVectorView f_grid)
{
  //Check consistency of input.

  assert ( is_size( scat_i_p,
                    f_grid.nelem(), 2, scat_i_p.nshelves(), 
                    scat_i_p.nbooks(), 
                    scat_za_grid.nelem(), scat_aa_grid.nelem(),
                    stokes_dim ));

  assert ( is_size( scat_i_lat,
                    f_grid.nelem(), scat_i_lat.nvitrines(), 2, 
                    scat_i_p.nbooks(), 
                    scat_za_grid.nelem(), scat_aa_grid.nelem(),
                    stokes_dim ));
  assert ( is_size( scat_i_lon,
                    f_grid.nelem(), scat_i_lat.nvitrines(), 
                    scat_i_p.nshelves(), 2, 
                    scat_za_grid.nelem(), scat_aa_grid.nelem(),
                    stokes_dim ));

  out3 << "\n" << "Get outgoing field from cloudbox boundary \n ";
  out3 << "    zenith_angle: " << rte_los[0] 
       << "    azimuth_angle: " << rte_los[1]+180 << "\n";
      
  // Boolians for pressure, latitude and longitude boundaries.
  bool on_p_bd = false;
  bool on_lat_bd = false;
  bool on_lon_bd = false;

  GridPos cloud_gp_p, cloud_gp_lat, cloud_gp_lon;
  cloud_gp_p = rte_gp_p;
  cloud_gp_lat = rte_gp_lat;
  cloud_gp_lon = rte_gp_lon;
  cloud_gp_p.idx -= cloudbox_limits[0];  
  cloud_gp_lat.idx -= cloudbox_limits[2];
  cloud_gp_lon.idx -= cloudbox_limits[4]; 

  cloudbox_boundary_check(on_p_bd, on_lat_bd, on_lon_bd, cloud_gp_p, 
                          cloud_gp_lat, cloud_gp_lon, cloudbox_limits);
      
  //Arrays to store grid positions
  GridPos gp_za, gp_aa;
  gridpos(gp_za, scat_za_grid, rte_los[0]);
  gridpos(gp_aa, scat_aa_grid, rte_los[1]+180);
     
  //Matrices to store interpolation weights
  Vector itw_angle(4);
  interpweights(itw_angle, gp_za, gp_aa);
      
  Vector itw_p_bd(4);
  interpweights(itw_p_bd, cloud_gp_lat, cloud_gp_lon);
      
  Vector itw_lat_bd(4);
  interpweights(itw_lat_bd, cloud_gp_p, cloud_gp_lon);
      
  Vector itw_lon_bd(4);
  interpweights(itw_lon_bd, cloud_gp_p, cloud_gp_lat);
     
  // Interpolation of the radiation field on the location
  // of the intersection point
  Matrix i_out_f(scat_za_grid.nelem(), scat_aa_grid.nelem());
      
  for(Index i = 0; i < stokes_dim; i++)
    {
      for(Index f_index = 0; f_index < f_grid.nelem(); f_index ++)
        {
          for(Index za_index = 0; za_index < scat_za_grid.nelem();
              za_index ++)
            {
              for(Index aa_index = 0; aa_index < scat_aa_grid.nelem();
                  aa_index ++)
                {
                  // The requested point lies on lower pressure boundary
                  if ( on_p_bd && cloud_gp_p.idx == 0)
                    {
                      // Interpolate on the right latitude longitude 
                      // position
                      i_out_f(za_index, aa_index) = 
                        interp( itw_p_bd, 
                                scat_i_p(f_index, 0, joker, joker,
                                         za_index, aa_index, i), 
                                cloud_gp_lat, cloud_gp_lon);
                    }
                  // The requested point lies on upper pressure boundary
                  else if ( on_p_bd && 
                            cloud_gp_p.idx == 
                            cloudbox_limits[1] - cloudbox_limits[0] - 1)
                    {
                      i_out_f(za_index, aa_index) = 
                        interp( itw_p_bd, 
                                scat_i_p(f_index, 1, joker, joker,
                                         za_index, aa_index, i),
                                cloud_gp_lat, cloud_gp_lon);
                    }
                  // The requested point lies on lower latitude boundary
                  else if (on_lat_bd && 
                           cloud_gp_lat.idx == 0)
                    {
                      i_out_f(za_index, aa_index) =
                        interp( itw_lat_bd, 
                                scat_i_lat(f_index, joker, 0, joker,
                                           za_index, aa_index, i),
                                cloud_gp_p, cloud_gp_lon);
                    }
                  // The requested point lies on upper latitude boundary
                  else if (on_lat_bd && 
                           cloud_gp_lat.idx == 
                           cloudbox_limits[3] - cloudbox_limits[2] - 1)
                    {
                      i_out_f(za_index, aa_index) = 
                        interp( itw_lat_bd, 
                                scat_i_lat(f_index, joker, 1, joker,
                                           za_index, aa_index, i),
                                cloud_gp_p, cloud_gp_lon);
                    }
                  // The requested point lies on lower longitude boundary
                  else if(on_lon_bd && 
                          cloud_gp_lon.idx == 0)
                    {
                      i_out_f(za_index, aa_index) = 
                        interp( itw_lon_bd, 
                                scat_i_lon(f_index, joker, joker, 0,
                                           za_index, aa_index, i),
                                cloud_gp_p, cloud_gp_lat);
                    }
                  // The requested point lies on upper longitude boundary
                  else if (on_lon_bd && 
                           cloud_gp_lon.idx == 
                           cloudbox_limits[5] - cloudbox_limits[4] - 1)
                    {
                      i_out_f(za_index, aa_index) = 
                        interp( itw_lon_bd, 
                                scat_i_lon(f_index, joker, joker, 1,
                                           za_index, aa_index, i),
                                cloud_gp_p, cloud_gp_lat);
                    }
                  else
                    { 
                      throw runtime_error(
                                          "Error in CloudboxGetOutgoing:"
                                          "The point where you want to \n"
                                          "calculcate the radiation "
                                          "coming out of the cloud is \n"
                                          "not on the cloudbox "
                                          "boundary."
                                          );
                    }
                }
            }
     
          // Outgoing radiances has to be interpolated on the
          // required direction.
             
          Numeric i_out_los;
             
          // Do the interpolation.
          i_out_los = interp(itw_angle, i_out_f, gp_za, gp_aa);
             
          //If radiance is greater than 1e-12, there must be a numerical 
          //error.
          assert(i_out_los < 1e-12);
             
          //Put the value into the matrix:
          i_out(f_index, i) = i_out_los;
             
        }//end frequency loop
    }//end stokes_dim loop
}

   
//! Check, which cloudbox surface is intersected and adjust grid positions.
/* 
   This method checks on which cloudbox boundary the propagation point,
   for which the outgoing radiation is required, lies. 
   One has to compare the grid-positions with cloudbox limits.
   It is not sufficient to take just the index as it is  
   ambiguous. If for example a point lies on the lower pressure 
   boundary of the cloudbox, the index of the grid-position for 
   the pressure coordinate can correspond to cloudbox_limits[0]
   or to cloudbox_limits[0]-1 depending on the interpolation 
   weights. 

   \param op_p_bd  Boolian. True if point lies on pressure surface.
   \param on_lat_bd True if point lies on latitude surface.
   \param on_lon_bd True if point lies on longitude surface.
   \param cloud_gp_p Grid position in cloudbox (pressure).
   \param cloud_gp_p Grid position in cloudbox (latitude).
   \param cloud_gp_p Grid position in cloudbox (longitude).
   \param cloudbox_limits Cloudbox limits. 
  
   \author Claudia Emde
   \date 2004-03-20
*/
void cloudbox_boundary_check(// Output:
                             bool& on_p_bd,
                             bool& on_lat_bd,
                             bool& on_lon_bd,
                             GridPos& cloud_gp_p,
                             GridPos& cloud_gp_lat,
                             GridPos& cloud_gp_lon,
                             // Input:
                             const ArrayOfIndex& cloudbox_limits
                             )
{
   
  //
     
  const Numeric TOL = 1e-6;
  // If the intersection points lies exactly on a 
  // lower boundary the gridposition index is 
  // increases by one and the first interpolation weight 
  // is set to 0.

  if (cloud_gp_p.idx == -1 &&
      abs(cloud_gp_p.fd[0]-1.) < TOL)
    {
      cloud_gp_p.idx = 0;
      cloud_gp_p.fd[0] = 0.;
      cloud_gp_p.fd[1] = 1.;
    }
                            
  if (cloud_gp_lat.idx == -1 &&
      abs(cloud_gp_lat.fd[0]-1.) < TOL)
    {
      cloud_gp_lat.idx = 0;
      cloud_gp_lat.fd[0] = 0.;
      cloud_gp_lat.fd[1] = 1.;
    }
      
  if (cloud_gp_lon.idx == -1 &&
      abs(cloud_gp_lon.fd[0]-1.) < TOL)
    {
      cloud_gp_lon.idx = 0;
      cloud_gp_lon.fd[0] = 0.;
      cloud_gp_lon.fd[1] = 1.;
    }
      
  // If the intersection points lies exactly on a 
  // upper boundary the gridposition index is 
  // reduced by one and the first interpolation weight 
  // is set to 1. This modification is necessary because otherwise
  // the interpolation does not work.
      
  if (cloud_gp_p.idx == (cloudbox_limits[1] - cloudbox_limits[0]) &&
      abs(cloud_gp_p.fd[0]) < TOL)
    {
      cloud_gp_p.idx -= 1;
      cloud_gp_p.fd[0] = 1.;
      cloud_gp_p.fd[1] = 0.;
    }
      
  if (cloud_gp_lat.idx == (cloudbox_limits[3] - cloudbox_limits[2]) &&
      abs(cloud_gp_lat.fd[0]) < TOL)
    {
      cloud_gp_lat.idx -= 1;
      cloud_gp_lat.fd[0] = 1.;
      cloud_gp_lat.fd[1] = 0.;
    }
      
  if (cloud_gp_lon.idx == (cloudbox_limits[5] - cloudbox_limits[4]) &&
      abs(cloud_gp_lon.fd[0]) < TOL)
    {
      cloud_gp_lon.idx -= 1;
      cloud_gp_lon.fd[0] = 1.;
      cloud_gp_lon.fd[1] = 0.;
    }
      
      
  //Check, on which boundary(s) the intersection point lies
  if(cloud_gp_p.idx == 0 ||
     cloud_gp_p.idx == cloudbox_limits[1] - cloudbox_limits[0] - 1)
    {
      on_p_bd = true;
      out3 << "(pressure surface, indx: "<< cloud_gp_p.idx << ")\n";
    }

  if(cloud_gp_lat.idx == 0 ||
     cloud_gp_lat.idx == cloudbox_limits[3] - cloudbox_limits[2] -1)
    {
      on_lat_bd = true;
      out3 << "(latitude surface), indx: " << cloud_gp_lat.idx << ") \n";
    }
      
  if(cloud_gp_lon.idx == 0 ||
     cloud_gp_lon.idx == cloudbox_limits[5] - cloudbox_limits[4] -1)
    {
      on_lon_bd = true;
      out3 << "(longitude surface), indx: " << cloud_gp_lon.idx << ") \n";
    }
      
}
