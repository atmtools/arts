/* Copyright (C) 2004-2012 Mattias Ekstrom  <ekstrom@rss.chalmers.se>

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
  \file   jacobian.cc
  \author Mattias Ekstrom <ekstrom@rss.chalmers.se>
  \date   2004-09-14

  \brief  Routines for setting up the jacobian.
*/

#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "jacobian.h"
#include "special_interp.h"
#include "physics_funcs.h"
#include "lin_alg.h"
#include "rte.h"

extern const String  ABSSPECIES_MAINTAG;
extern const String  TEMPERATURE_MAINTAG;
extern const String  WIND_MAINTAG;
extern const String  MAGFIELD_MAINTAG;
extern const String  FLUX_MAINTAG;
extern const String  PROPMAT_SUBSUBTAG;


ostream& operator << (ostream& os, const RetrievalQuantity& ot)
{
  return os << "\n       Main tag = " << ot.MainTag() 
            << "\n       Sub  tag = " << ot.Subtag()
            << "\n           Mode = " << ot.Mode()
            << "\n     Analytical = " << ot.Analytical();
}




/*===========================================================================
  === Help sub-functions to handle analytical jacobians (in alphabetical order)
  ===========================================================================*/

//! diy_from_path_to_rgrids
/*!
    Maps jacobian data for points along the propagation path, to
    jacobian retrieval grid data.

    \param   diy_dx              Out: Jacobians for selected retrieval grids.
    \param   jacobian_quantity   As the WSV.
    \param   diy_dpath           Jacobians along the propagation path.
    \param   atmosphere_dim      As the WSV.
    \param   ppath               As the WSV.
    \param   ppath_p             The pressure at each ppath point.

    \author Patrick Eriksson 
    \date   2009-10-08
*/
// Small help function, to make the code below cleaner
void from_dpath_to_dx(
        MatrixView   diy_dx,
   ConstMatrixView   diy_dq,
   const Numeric&    w )
{
  for( Index irow=0; irow<diy_dx.nrows(); irow++ )
    { 
      for( Index icol=0; icol<diy_dx.ncols(); icol++ )
        { diy_dx(irow,icol) += w * diy_dq(irow,icol); }
    }
}
//
void diy_from_path_to_rgrids(
         Tensor3View          diy_dx,
   const RetrievalQuantity&   jacobian_quantity,
   ConstTensor3View           diy_dpath,
   const Index&               atmosphere_dim,
   const Ppath&               ppath,
   ConstVectorView            ppath_p )
{
   if( jacobian_quantity.Integration() )
   {
       diy_dx = diy_dpath;
       return;
   }
    
  // We want here an extrapolation to infinity -> 
  //                                        extremly high extrapolation factor
  const Numeric   extpolfac = 1.0e99;

  if( ppath.np > 1 )  // Otherwise nothing to do here
    {
      // Pressure
      Index            nr1 = jacobian_quantity.Grids()[0].nelem();
      ArrayOfGridPos   gp_p(ppath.np);
      if( nr1 > 1 )
        {
          p2gridpos( gp_p, jacobian_quantity.Grids()[0], ppath_p, extpolfac );
          jacobian_type_extrapol( gp_p );
        }
      else
        { gp4length1grid( gp_p ); }        

      // Latitude
      Index            nr2 = 1;
      ArrayOfGridPos   gp_lat;
      if( atmosphere_dim > 1 )
        {          
          gp_lat.resize(ppath.np);
          nr2    = jacobian_quantity.Grids()[1].nelem();
          if( nr2 > 1 )
            {
              gridpos( gp_lat, jacobian_quantity.Grids()[1], 
                       ppath.pos(joker,1), extpolfac );
              jacobian_type_extrapol( gp_lat );
            }
          else
            { gp4length1grid( gp_lat ); }
        }

      // Longitude
      ArrayOfGridPos   gp_lon;
      if( atmosphere_dim > 2 )
        {
          gp_lon.resize(ppath.np);
          if( jacobian_quantity.Grids()[2].nelem() > 1 )
            {          
              gridpos( gp_lon, jacobian_quantity.Grids()[2], 
                       ppath.pos(joker,2), extpolfac );
              jacobian_type_extrapol( gp_lon );
            }
          else
            { gp4length1grid( gp_lon ); }
        }
      
      //- 1D
      if( atmosphere_dim == 1 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              if( gp_p[ip].fd[1] > 0 )
                {
                  from_dpath_to_dx( diy_dx(gp_p[ip].idx,joker,joker),
                                    diy_dpath(ip,joker,joker), gp_p[ip].fd[1] );
                }
              if( gp_p[ip].fd[0] > 0 )
                {
                  from_dpath_to_dx( diy_dx(gp_p[ip].idx+1,joker,joker),
                                    diy_dpath(ip,joker,joker), gp_p[ip].fd[0] );
                }
            }
        }

      //- 2D
      else if( atmosphere_dim == 2 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              Index   ix = nr1*gp_lat[ip].idx + gp_p[ip].idx;
              // Low lat, low p
              if( gp_lat[ip].fd[1]>0 && gp_p[ip].fd[1]>0 )
                from_dpath_to_dx( diy_dx(ix,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                                  gp_lat[ip].fd[1]*gp_p[ip].fd[1] );
              // Low lat, high p
              if( gp_lat[ip].fd[1]>0 && gp_p[ip].fd[0]>0 )
                from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                                  gp_lat[ip].fd[1]*gp_p[ip].fd[0] );
              // High lat, low p
              if( gp_lat[ip].fd[0]>0 && gp_p[ip].fd[1]>0 )
                from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                                  gp_lat[ip].fd[0]*gp_p[ip].fd[1] );
              // High lat, high p
              if( gp_lat[ip].fd[0]>0 && gp_p[ip].fd[0]>0 )
                from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                                  gp_lat[ip].fd[0]*gp_p[ip].fd[0] );
            }
        }

      //- 3D
      else if( atmosphere_dim == 3 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              Index   ix = nr2*nr1*gp_lon[ip].idx +
                           nr1*gp_lat[ip].idx + gp_p[ip].idx;
              // Low lon, low lat, low p
              if( gp_lon[ip].fd[1]>0 && gp_lat[ip].fd[1]>0 && gp_p[ip].fd[1]>0 )
                from_dpath_to_dx( diy_dx(ix,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[1]*gp_p[ip].fd[1]);
              // Low lon, low lat, high p
              if( gp_lon[ip].fd[1]>0 && gp_lat[ip].fd[1]>0 && gp_p[ip].fd[0]>0 )
                from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[1]*gp_p[ip].fd[0]);
              // Low lon, high lat, low p
              if( gp_lon[ip].fd[1]>0 && gp_lat[ip].fd[0]>0 && gp_p[ip].fd[1]>0 )
                from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[0]*gp_p[ip].fd[1]);
              // Low lon, high lat, high p
              if( gp_lon[ip].fd[1]>0 && gp_lat[ip].fd[0]>0 && gp_p[ip].fd[0]>0 )
                from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[0]*gp_p[ip].fd[0]);

              // Increase *ix* (to be valid for high lon level)
              ix += nr2*nr1;

              // High lon, low lat, low p
              if( gp_lon[ip].fd[0]>0 && gp_lat[ip].fd[1]>0 && gp_p[ip].fd[1]>0 )
                from_dpath_to_dx( diy_dx(ix,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[1]*gp_p[ip].fd[1]);
              // High lon, low lat, high p
              if( gp_lon[ip].fd[0]>0 && gp_lat[ip].fd[1]>0 && gp_p[ip].fd[0]>0 )
                from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[1]*gp_p[ip].fd[0]);
              // High lon, high lat, low p
              if( gp_lon[ip].fd[0]>0 && gp_lat[ip].fd[0]>0 && gp_p[ip].fd[1]>0 )
                from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[0]*gp_p[ip].fd[1]);
              // High lon, high lat, high p
              if( gp_lon[ip].fd[0]>0 && gp_lat[ip].fd[0]>0 && gp_p[ip].fd[0]>0 )
                from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                  diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[0]*gp_p[ip].fd[0]);
            }
        }
    }
}



//! Help function for analytical jacobian calculations
/*!
    The function determines which terms in jacobian_quantities that are 
    analytical absorption species and temperature jacobians. 

    *abs_species_i* and *is_t* shall be sized to have the same length
    as *jacobian_quantities*. For analytical absorption species
    jacobians, *abs_species_i* is set to the matching index in
    *abs_species*. Otherwise, to -1. For analytical temperature
    jacobians, *is_t* is set to 1. Otherwise to 0.

    \param   abs_species_i         Out: Matching index in abs_species 
    \param   is_t                  Out: Flag for: Is a temperature jacobian?
    \param   jacobian_quantities   As the WSV.
    \param   abs_species           As the WSV.


    \author Patrick Eriksson 
    \date   2009-10-07
*/
void get_pointers_for_analytical_jacobians( 
         ArrayOfIndex&               abs_species_i, 
         ArrayOfIndex&               is_t,
         ArrayOfIndex&               wind_i,
         ArrayOfIndex&               magfield_i,
         ArrayOfIndex&               flux_i,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfSpeciesTag&   abs_species )
{

  FOR_ANALYTICAL_JACOBIANS_DO( 
    //
    if( jacobian_quantities[iq].MainTag() == TEMPERATURE_MAINTAG  && jacobian_quantities[iq].SubSubtag() == PROPMAT_SUBSUBTAG  )
      { is_t[iq] = JAC_IS_T_FROM_PROPMAT; }
    else if( jacobian_quantities[iq].MainTag() == TEMPERATURE_MAINTAG )
      { is_t[iq] = JAC_IS_T_SEMI_ANALYTIC; }
    else
      { is_t[iq] = JAC_IS_NONE; }
    //
    if( jacobian_quantities[iq].MainTag() == ABSSPECIES_MAINTAG )
      {
        if( jacobian_quantities[iq].SubSubtag() == PROPMAT_SUBSUBTAG )
        {
            bool test_available=false;
            for(Index ii=0; ii<abs_species.nelem(); ii++)
            {
                if( abs_species[ii][0].Species()==SpeciesTag(jacobian_quantities[iq].Subtag()).Species() )
                {
                    test_available=true;
                    abs_species_i[iq]=ii;
                    break;
                }
            }
            if(!test_available)
            {
                ostringstream os;
                os << "Could not find " << jacobian_quantities[iq].Subtag() <<
                "in species of abs_species.\n";
                throw std::runtime_error(os.str());
            }
        }
        else
        {
            ArrayOfSpeciesTag  atag;
            array_species_tag_from_string( atag, jacobian_quantities[iq].Subtag() );
            abs_species_i[iq] = chk_contains( "abs_species", abs_species, atag );
        }
      }
    else
      { abs_species_i[iq] = -1; }
    //
    if( jacobian_quantities[iq].MainTag() == WIND_MAINTAG && jacobian_quantities[iq].SubSubtag() == PROPMAT_SUBSUBTAG  )
    {
        // Map u, v and w to 1, 2 and 3, respectively
        char c = jacobian_quantities[iq].Subtag()[0];
        const Index test = Index( c ) - 116;
        if( test == 1 )
            wind_i[iq] = JAC_IS_WIND_U_FROM_PROPMAT;
        else if(test == 2 )
            wind_i[iq] = JAC_IS_WIND_V_FROM_PROPMAT;
        else if(test == 3 )
            wind_i[iq] = JAC_IS_WIND_W_FROM_PROPMAT;
        else if(test == (Index('s')-116) )
            wind_i[iq] = JAC_IS_WIND_ABS_FROM_PROPMAT;
    }
    else if( jacobian_quantities[iq].MainTag() == WIND_MAINTAG )
      {
        // Map u, v and w to 1, 2 and 3, respectively
          char c = jacobian_quantities[iq].Subtag()[0];
          const Index test = Index( c ) - 116;
          if( test == 1 )
              wind_i[iq] = JAC_IS_WIND_U_SEMI_ANALYTIC;
          else if(test == 2 )
              wind_i[iq] = JAC_IS_WIND_V_SEMI_ANALYTIC;
          else if(test == 3 )
              wind_i[iq] = JAC_IS_WIND_W_SEMI_ANALYTIC;
      }
    else
      { wind_i[iq] = JAC_IS_NONE; }
    //
    if( jacobian_quantities[iq].MainTag() == MAGFIELD_MAINTAG && jacobian_quantities[iq].SubSubtag() == PROPMAT_SUBSUBTAG )
    {
        // Map u, v and w to 1, 2 and 3, respectively
        char c = jacobian_quantities[iq].Subtag()[0];
        const Index test = Index( c ) - 116;
        if( test == 1 )
            magfield_i[iq] = JAC_IS_MAG_U_FROM_PROPMAT;
        else if(test == 2 )
            magfield_i[iq] = JAC_IS_MAG_V_FROM_PROPMAT;
        else if(test == 3 )
            magfield_i[iq] = JAC_IS_MAG_W_FROM_PROPMAT;
    }
    else if( jacobian_quantities[iq].MainTag() == MAGFIELD_MAINTAG )
      {
        // Map u, v and w to 1, 2 and 3, respectively
        char c = jacobian_quantities[iq].Subtag()[0];
        const Index test = Index( c ) - 116;
        if( test == 1 )
            magfield_i[iq] = JAC_IS_MAG_U_SEMI_ANALYTIC;
        else if(test == 2 )
            magfield_i[iq] = JAC_IS_MAG_V_SEMI_ANALYTIC;
        else if(test == 3 )
            magfield_i[iq] = JAC_IS_MAG_W_SEMI_ANALYTIC;
      }
    else
      { magfield_i[iq] = JAC_IS_NONE; }
    //
    if( jacobian_quantities[iq].MainTag() == FLUX_MAINTAG ) { flux_i[iq] = 1; }
    else  { flux_i[iq] = 0; }
  )
}





/*===========================================================================
  === Other functions, in alphabetical order
  ===========================================================================*/

//! Calculate the number density field
/*! 
   This function returns the number density for each grid point in the 
   Tensor3View.
   
   \param nd  The number density field
   \param p   The pressure grid
   \param t   The temperature field
   
   \author Mattias Ekstrom
   \date   2005-06-03
*/
void calc_nd_field(       Tensor3View& nd,
                    const VectorView&  p,
                    const Tensor3View& t)
{
  assert( nd.npages()==t.npages() );
  assert( nd.nrows()==t.nrows() );               
  assert( nd.ncols()==t.ncols() );
  assert( nd.npages()==p.nelem() );
  
  for (Index p_it=0; p_it<nd.npages(); p_it++)
  {
    for (Index lat_it=0; lat_it<nd.nrows(); lat_it++)
    {
      for (Index lon_it=0; lon_it<nd.ncols(); lon_it++)
      {
        nd(p_it,lat_it,lon_it) = number_density( p[p_it], t(p_it,lat_it,lon_it));
      }
    }
  }
}



//! Check that the retrieval grids are defined for each atmosphere dim
/*!
   This function checks for the given atmosphere dimension that 
     I)  the retrieval grids are defined 
     II) and that they are covered by the corresponding atmospheric grid. 
   If not the return is false and an output string is created to print 
   the error to the user. If the grids are ok they are stored in an array 
   and true is  returned.
   
   \param grids         The array of retrieval grids.
   \param os            The output string stream.
   \param p_grid        The atmospheric pressure grid
   \param lat_grid      The atmospheric latitude grid
   \param lon_grid      The atmospheric longitude grid
   \param p_retr        The pressure retrieval grid.
   \param lat_retr      The latitude retrieval grid.
   \param lon_retr      The longitude retrieval grid.
   \param p_retr_name   The control file name used for the pressure retrieval grid.
   \param lat_retr_name The control file name for the latitude retrieval grid.
   \param lon_retr_name The control file name for the longitude retrieval grid.
   \param dim           The atmosphere dimension
   \return              Boolean for check.
   
   \author Mattias Ekstrom
   \date   2005-05-11
*/ 
bool check_retrieval_grids(       ArrayOfVector& grids,
                                  ostringstream& os,
                            const Vector&        p_grid,
                            const Vector&        lat_grid,
                            const Vector&        lon_grid,
                            const Vector&        p_retr,
                            const Vector&        lat_retr,
                            const Vector&        lon_retr,
                            const String&        p_retr_name,
                            const String&        lat_retr_name,
                            const String&        lon_retr_name,
                            const Index&         dim)
{
  if ( p_retr.nelem()==0 )
    {
      os << "The grid vector *" << p_retr_name << "* is empty,"
         << " at least one pressure level\n"
         << "should be specified.";
      return false;
    }
  else if( !is_decreasing( p_retr ) )
    {
      os << "The pressure grid vector *" << p_retr_name << "* is not a\n"
         << "strictly decreasing vector, which is required.";
      return false;      
    }
  else if ( log(p_retr[0])> 1.5*log(p_grid[0])-0.5*log(p_grid[1]) || 
            log(p_retr[p_retr.nelem()-1])<1.5*log(p_grid[p_grid.nelem()-1])-
                                          0.5*log(p_grid[p_grid.nelem()-2])) 
    {
      os << "The grid vector *" << p_retr_name << "* is not covered by the\n"
         << "corresponding atmospheric grid.";
      return false;
    }
  else
    {
      // Pressure grid ok, add it to grids
      grids[0]=p_retr;  
    }

  if (dim>=2)
  {
    // If 2D and 3D atmosphere, check latitude grid
    if ( lat_retr.nelem()==0 )
    {
      os << "The grid vector *" << lat_retr_name << "* is empty,"
         << " at least one latitude\n"
         << "should be specified for a 2D/3D atmosphere.";
      return false;
    }
    else if( !is_increasing( lat_retr ) )
    {
      os << "The latitude grid vector *" << lat_retr_name << "* is not a\n"
         << "strictly increasing vector, which is required.";
      return false;      
    }
    else if ( lat_retr[0]<1.5*lat_grid[0]-0.5*lat_grid[1] || 
              lat_retr[lat_retr.nelem()-1]>1.5*lat_grid[lat_grid.nelem()-1]-
                                           0.5*lat_grid[lat_grid.nelem()-2] )
    {
      os << "The grid vector *" << lat_retr_name << "* is not covered by the\n"
         << "corresponding atmospheric grid.";
      return false;
    }
    else
    {
      // Latitude grid ok, add it to grids
      grids[1]=lat_retr;
    }
    if (dim==3)
    {
      // For 3D atmosphere check longitude grid
      if ( lon_retr.nelem()==0 )
      {
        os << "The grid vector *" << lon_retr_name << "* is empty,"
           << " at least one longitude\n"
           << "should be specified for a 3D atmosphere.";
        return false;
      }
      else if( !is_increasing( lon_retr ) )
      {
      os << "The longitude grid vector *" << lon_retr_name << "* is not a\n"
         << "strictly increasing vector, which is required.";
      return false;      
      }
      else if ( lon_retr[0]<1.5*lon_grid[0]-0.5*lon_grid[1] || 
                lon_retr[lon_retr.nelem()-1]>1.5*lon_grid[lon_grid.nelem()-1]-
                                             0.5*lon_grid[lon_grid.nelem()-2] )
      {
        os << "The grid vector *" << lon_retr_name << "* is not covered by the\n"
           << "corresponding atmospheric grid.";
        return false;
      }
      else
      {
        // Longitude grid ok, add it to grids      
        grids[2]=lon_retr;
      }
    }
  }
  return true;
}             



//! Calculate array of GridPos for perturbation interpolation
/*!
   This function constructs a perturbation grid which consists of the
   given retrieval grid with an extra endpoint added at each end.
   These endpoints lies outside the atmospheric grid. This enables the
   interpolation of an perturbation on the perturbation grid to be
   interpolated to the atmospheric grid. For this reason the function
   returns an ArrayOfGridPos. 
   
   If the atmospheric grid is a pressure grid, interpolation is made
   in logarithm of the atmospheric grid.
   
   \param gp          Array of GridPos for interpolation.
   \param atm_grid    Atmospheric grid.
   \param jac_grid    Retrieval grid.
   \param is_pressure True for pressure grid 
   
   \author Mattias Ekstrom
   \date   2005-05-12
*/   
void get_perturbation_gridpos(       ArrayOfGridPos& gp,
                               const Vector&         atm_grid,
                               const Vector&         jac_grid,
                               const bool&           is_pressure)
{
  Index nj = jac_grid.nelem();
  Index na = atm_grid.nelem();
  Vector pert(nj+2);
  
  // Create perturbation grid, with extension outside the atmospheric grid
  if ( is_pressure )
  {
    pert[0] = atm_grid[0]*10.0;
    pert[nj+1] = atm_grid[na-1]*0.1;
  }
  else
  {    
    pert[0] = atm_grid[0]-1.0;
    pert[nj+1] = atm_grid[na-1]+1.0;
  }
  pert[Range(1,nj)] = jac_grid;
  
  // Calculate the gridpos
  gp.resize(na);
  if( is_pressure ){
    p2gridpos( gp, pert, atm_grid);
  }
  else
  { 
    gridpos( gp, pert, atm_grid);
  }
}



//! Get limits for perturbation of a box
/*!
   This is a helper function that calculates the limits where the 
   perturbation should be added to the perturbation grid. 
   This is needed for example by the particle perturbation that only
   should be applied for the cloudbox. The limits are defined as the 
   outermost points lying within or just outside the box limits.
   
   The atmospheric limits should be given in the same unit as the
   perturbation grid. And only the first and last element will be 
   considered as limits. 
   
   Assertions are used to perform checks. The input grids are 
   checked so that the atmospheric limits are containg within 
   the perturbation grid. The limit indices are checked so 
   that they are ordered in increasing order before return.
   
   \param limit     The limit indices in the perturbation grid
   \param pert_grid The perturbation grid
   \param atm_limit The atmospheric limits of the box.

   \author Mattias Ekstrom
   \date   2005-02-25
*/   
void get_perturbation_limit(       ArrayOfIndex& limit,
                             const Vector&       pert_grid,
                             const Vector&       atm_limit)
{
  limit.resize(2);
//   Index np = pert_grid.nelem()-1;
  Index na = atm_limit.nelem()-1;
  
  // If the field is ordered in decreasing order set the
  // increment factor to -1
  Numeric inc = 1;
  if (is_decreasing(pert_grid))
    inc = -1;

  // Check that the pert_grid is encompassing atm_limit
//   assert( inc*pert_grid[0] < inc*atm_limit[0] &&
//           inc*pert_grid[np] > inc*atm_limit[na]);
      
  // Find first limit, check if following value is above lower limit
  limit[0]=0;
  while (inc*pert_grid[limit[0]+1] < inc*atm_limit[0]) 
  {
    limit[0]++;
  }
  
  // Find last limit, check if previous value is below upper limit
  limit[1]=pert_grid.nelem();
  while (inc*pert_grid[limit[1]-1] > inc*atm_limit[na]) 
  {
    limit[1]--;
  }
  // Check that the limits are ok
  assert(limit[1]>limit[0]);
  
}

                             

//! Get range for perturbation 
/*!
   This is a helper function that calculates the range in which the 
   perturbation should be added to the perturbation grid. This is needed
   to handle the edge effects. At the edges we want the perturbation to 
   continue outwards. 
   
   \param range     The range in the perturbation grid.
   \param index     The index of the perturbation in the retrieval grid.
   \param length    The length of retrieval grid
   
   \author Mattias Ekstrom
   \date   2004-10-14
*/   
void get_perturbation_range(       Range& range,
                             const Index& index,
                             const Index& length)
{
  if (index==0)
    range = Range(index,2);
  else if (index==length-1)
    range = Range(index+1,2);
  else 
    range = Range(index+1,1);

}



//! Adopts grid postions to extrapolation used for jacobians
/*!
  The standard interpolation scheme applies a linear extrapolation, while for
  the jacobians the extrapolation can be seen as a "closest" interpolation.
  That is, for points outisde the covered grid, the value at closest end point
  is taken. And here extrapolation to +-Inf is allowed.

  This function modifies grid positions to jacobaina extrapolation approach.
  For efficiency, the input grid positions are not asserted at all, and
  "extrapolation points" are identified simply  by a fd outside [0,1].

  \param[in/out] gp   Array of grid positions.

  \author Patrick Eriksson 
  \date   2015-09-10
*/
void jacobian_type_extrapol( ArrayOfGridPos&   gp )
{
  for( Index i=0; i<gp.nelem(); i++ )
    { 
      if( gp[i].fd[0] < 0 ) 
        {  
          gp[i].fd[0] = 0; 
          gp[i].fd[1] = 1; 
        }
      else if( gp[i].fd[0] > 1 ) 
        {  
          gp[i].fd[0] = 1; 
          gp[i].fd[1] = 0; 
        }
    }
}


//! Calculate the 1D perturbation for a relative perturbation.
/*!
   This is a helper function that interpolated the perturbation field for
   a 1D relative perturbation onto the atmospheric field. 
   
   \param field     The interpolated perturbation field.
   \param p_gp      The GridPos for interpolation.
   \param p_pert_n  The number of perturbations.
   \param p_range   The perturbation range in the perturbation grid.
   \param size      The size of the perturbation.
   \param method    Relative perturbation==0, absolute==1
   
   \author Mattias Ekstrom
   \date   2005-05-11
*/   
void perturbation_field_1d(       VectorView      field,
                            const ArrayOfGridPos& p_gp,
                            const Index&          p_pert_n,
                            const Range&          p_range,
                            const Numeric&        size,
                            const Index&          method)
{
  // Here we only perturb a vector
  Vector pert(field.nelem());
  Matrix itw(p_gp.nelem(),2);
  interpweights(itw,p_gp);
  // For relative pert_field should be 1.0 and for absolute 0.0
  Vector pert_field(p_pert_n,1.0-(Numeric)method);
  pert_field[p_range] += size;
  interp( pert, itw, pert_field, p_gp);
  if (method==0)
  {
    field *= pert;
  }
  else
  {
    field += pert;
  }
}            



//! Calculate the 2D perturbation for a relative perturbation.
/*!
   This is a helper function that interpolated the perturbation field for
   a 2D relative perturbation onto the atmospheric field. 
   
   \param field       The interpolated perturbation field.
   \param p_gp        The GridPos for interpolation in the 1st dim.
   \param lat_gp      The GridPos for interpolation in the 2nd dim.
   \param p_pert_n    The number of perturbations in the 1st dim.
   \param lat_pert_n  The number of perturbations in the 2nd dim.
   \param p_range     The perturbation range in the 1st dim.
   \param lat_range   The perturbation range in the 2nd dim.
   \param size        The size of the perturbation.
   \param method      Relative perturbation==0, absolute==1
   
   \author Mattias Ekstrom
   \date   2005-05-11
*/   
void perturbation_field_2d(       MatrixView      field,
                            const ArrayOfGridPos& p_gp,
                            const ArrayOfGridPos& lat_gp,
                            const Index&          p_pert_n,
                            const Index&          lat_pert_n,
                            const Range&          p_range,
                            const Range&          lat_range,
                            const Numeric&        size,
                            const Index&          method)
{
  // Here we perturb a matrix
  Matrix pert(field.nrows(),field.ncols());
  Tensor3 itw(p_gp.nelem(),lat_gp.nelem(),4);
  interpweights(itw,p_gp,lat_gp);
  // Init pert_field to 1.0 for relative and 0.0 for absolute
  Matrix pert_field(p_pert_n,lat_pert_n,1.0-(Numeric)method);
  pert_field(p_range,lat_range) += size;
  interp( pert, itw, pert_field, p_gp, lat_gp);
  if (method==0)
  {
    field *= pert;
  }
  else
  { 
    field += pert;
  }
}            



//! Calculate the 3D perturbation for a relative perturbation.
/*!
   This is a helper function that interpolated the perturbation field for
   a 3D relative perturbation onto the atmospheric field. 
   
   \param field       The interpolated perturbation field.
   \param p_gp        The GridPos for interpolation in the 1st dim.
   \param lat_gp      The GridPos for interpolation in the 2nd dim.
   \param lon_gp      The GridPos for interpolation in the 3rd dim.
   \param p_pert_n    The number of perturbations in the 1st dim.
   \param lat_pert_n  The number of perturbations in the 2nd dim.
   \param lon_pert_n  The number of perturbations in the 3rd dim.
   \param p_range     The perturbation range in the 1st dim.
   \param lat_range   The perturbation range in the 2nd dim.
   \param lon_range   The perturbation range in the 3rd dim.
   \param size        The size of the perturbation.
   \param method      Set to 0 for relative, and 1 for absolute.
   
   \author Mattias Ekstrom
   \date   2005-05-11
*/   
void perturbation_field_3d(       Tensor3View     field,
                            const ArrayOfGridPos& p_gp,
                            const ArrayOfGridPos& lat_gp,
                            const ArrayOfGridPos& lon_gp,
                            const Index&          p_pert_n,
                            const Index&          lat_pert_n,
                            const Index&          lon_pert_n,
                            const Range&          p_range,
                            const Range&          lat_range,
                            const Range&          lon_range,
                            const Numeric&        size,
                            const Index&          method)
{
  // Here we need to perturb a tensor3
  Tensor3 pert(field.npages(),field.nrows(),field.ncols());
  Tensor4 itw(p_gp.nelem(),lat_gp.nelem(),lon_gp.nelem(),8);
  interpweights(itw,p_gp,lat_gp,lon_gp);
  // Init pert_field to 1.0 for relative and 0.0 for absolute
  Tensor3 pert_field(p_pert_n,lat_pert_n,lon_pert_n,1.0-(Numeric)method);
  pert_field(p_range,lat_range,lon_range) += size;
  interp( pert, itw, pert_field, p_gp, lat_gp, lon_gp);
  if (method==0)
  {
    field *= pert;
  }
  else
  {
    field += pert;
  }
}            



//! Calculates polynomial basis functions
/*!
   The basis function is b(x) = 1 for poly_coeff = 0. For higher
   coefficients, x^poly_coeff - m, where first the range covered by
   *x* is normalised to [-1,1] and m is selected in such way that
   sum(b) = 0.
   
   \param b            Calculated basis function.
   \param x            The grid over which the fit shall be performed.
   \param poly_coeff   Polynomial coefficient.
   
   \author Patrick Eriksson
   \date   2008-11-07
*/   
void polynomial_basis_func(
        Vector&   b,
  const Vector&   x,
  const Index&    poly_coeff )
{
  const Index l = x.nelem();
  
  assert( l > poly_coeff );

  if( b.nelem() != l )
    b.resize( l );

  if( poly_coeff == 0 )
    { b = 1.0; }
  else
    {
      const Numeric xmin = min( x );
      const Numeric dx = 0.5 * ( max( x ) - xmin );
      //
      for( Index i=0; i<l; i++ )
        {
          b[i] = ( x[i] - xmin ) / dx - 1.0;
          b[i] = pow( b[i], int(poly_coeff) );
        }
      //
      b -= mean( b );
    }  
}



//! vmrunitscf
/*!
    Scale factor for conversion between gas species units.

    The function finds the factor with which the total absorption of a
    gas species shall be multiplicated to match the selected
    (jacobian) unit. 

    \param   x      Out: scale factor
    \param   unit   Unit selected.
    \param   vmr    VMR value.
    \param   p      Pressure
    \param   t      Temperature.

    \author Patrick Eriksson 
    \date   2009-10-08
*/
void vmrunitscf(  
        Numeric&   x, 
  const String&    unit, 
  const Numeric&   vmr,
  const Numeric&   p,
  const Numeric&   t )
{  
  if( unit == "rel"  ||  unit == "logrel" )
    { x = 1; }
  else if( unit == "vmr" )
    {
        if(vmr==0)
        { 
            x = 0; 
            return;
        }
        x = 1 / vmr;
    }
  else if( unit == "nd" )
    {
        if(vmr==0)
        { 
            x = 0; 
            return;
        }
        x = 1 / ( vmr * number_density( p, t ) ); 
    }
  else
    {
        ostringstream os;
        os << "Allowed options for gas species jacobians are "
        "\"rel\", \"vmr\" and \"nd\".\nYou have selected: "<<unit<<std::endl; 
      throw std::runtime_error( os.str() );
    }
}




//! dxdvmrscf
/*!
    Scale factor for conversion of derivatives with respect to VMR.

    The function finds the factor with which a partial derivative with respect
    to gas species VMR shall be multiplicated to match the selected (jacobian)
    unit. The function was implemented for scaling of *dpropmat_clearsky_dx*
    but it could also be of use in other contexts.

    \param   x      Out: scale factor
    \param   unit   Unit selected.
    \param   vmr    VMR value.
    \param   p      Pressure
    \param   t      Temperature.

    \author Patrick Eriksson 
    \date   2015-12-11
*/
void dxdvmrscf(  
        Numeric&   x, 
  const String&    unit, 
  const Numeric&   vmr,
  const Numeric&   p,
  const Numeric&   t )
{
  if( unit == "rel"  ||  unit == "logrel" )
    { x = vmr; }
  else if( unit == "vmr" )
    { x = 1; }
  else if( unit == "nd" )
    { x = 1 / number_density( p, t ); }
  else
    {
        ostringstream os;
        os << "Allowed options for gas species jacobians are "
        "\"rel\", \"vmr\" and \"nd\".\nYou have selected: "<<unit<<std::endl; 
      throw std::runtime_error( os.str() );
    }
}


/*!
 *   The function helps to calculate the partial derivative of iy with respect
 *   to one input at one pressure.  Since one layer is influenced by two levels,
 *   the inputs for "this" and "next" levels are treated in parallel.  Note that
 *   "this" refers to the "ip" level in iyEmissionStandard and "next" refers to
 *   "ip+1" level in iyEmissionStandard.  Note that all vectors and matrices are
 *   for one frequency and one pressure level.  Parallelization over frequency
 *   should be possible. 
 * 
 *   \param   diydx_this                Out: diydx for this level
 *   \param   diydx_next                Out: diydx for next level
 *   \param   iy                        In:  iy for one frequency and this level
 *   \param   sibi                      In:  sibi for one frequency and this level
 *   \param   nlte_this                 In:  nlte source term for this level
 *   \param   nlte_next                 In:  nlte source term for next level
 *   \param   dnltedx_this              In:  partial derivative of nlte source term for this level
 *   \param   dnltedx_next              In:  partial derivative of nlte source term for next level
 *   \param   K_this                    In:  propagation matrix for this level
 *   \param   K_next                    In:  propagation matrix for next level
 *   \param   dKdx_this                 In:  partial derivative of propagation matrix for this level
 *   \param   dKdx_next                 In:  partial derivative of propagation matrix for next level
 *   \param   T_this                    In:  transmission matrix for this level
 *   \param   dTdx_this                 In:  partial derivative of transmission matrix for this level
 *   \param   dTdx_next                 In:  partial derivative of transmission matrix for next level
 *   \param   PiT_this                  In:  cumulative transmission matrix for this level
 *   \param   PiT_next                  In:  cumulative transmission matrix for next level
 *   \param   temperature_this          In:  temperature of this level
 *   \param   temperature_next          In:  temperature of next level
 *   \param   dt                        In:  Delta temperature in polarized HSE application
 *   \param   dBdx_this                 In:  partial derivative of Planck function for this level
 *   \param   dBdx_next                 In:  partial derivative of Planck function for next level
 *   \param   r                         In:  Layer thickness between this and next
 *   \param   stokes_dim                In:  Stokes dimension of inputs
 *   \param   extmat_case               In:  polarized or unpolarized transmission matrix?
 *   \param   do_Bsource                In:  calculate partial derivative due to Planck function?
 *   \param   do_HSE                    In:  calculate partial derivative due to hydrostatic equilibrium?
 *   \param   do_nonlte                 In:  calculate partial derivative due source in non-lte?
 * 
 * 
 *   \author Richard Larsson (adapted from most general jacobian in old iyEmissionStandard)
 *   \date   2015-12-10
 */
void get_diydx( VectorView diydx_this,
                VectorView diydx_next,
                Index& extmat_case,
                ConstVectorView iy,
                ConstVectorView sibi,
                ConstVectorView nlte_this,
                ConstVectorView nlte_next,
                ConstVectorView dnltedx_this,
                ConstVectorView dnltedx_next,
                ConstMatrixView K_this,
                ConstMatrixView K_next,
                ConstMatrixView dKdx_this,
                ConstMatrixView dKdx_next,
                ConstMatrixView T_this,
                ConstMatrixView dTdx_this,
                ConstMatrixView dTdx_next,
                ConstMatrixView PiT_this,
                ConstMatrixView PiT_next,
                const Numeric& temperature_this,
                const Numeric& temperature_next,
                const Numeric& dt,
                const Numeric& dBdx_this,
                const Numeric& dBdx_next,
                const Numeric& r,
                const Index& stokes_dim,
                const bool& do_Bsource,
                const bool& do_HSE,
                const bool& do_nonlte )
{
    // Diagonal transmission matrix
    if( extmat_case == 1 )
    {
        
        /*
         * Solves ds/dx = T [ dTi/dx * (si-bi-Ki^(-1)ji) + 
         *        (1 - Ti) * ( dbi/dx + K^(-2)dKi/dxji + K^(-1)dji/dx ) ]
         * 
         * In extmat_case 1, dTi/dx = T*dKi/dx.
         * 
         * (si-bi-Ki^(-1)ji) is sibi
         * T                 is PiT_this                   (NB PiT_next is PiT_this*T_this)
         * dTi/dx            is PiT_next*dKi/dx*r/2        (from T being scalar/diagonal exponent function)
         * dbi/dx            is dBdx_this/2, dBdx_next/2   (note how outer loop takes care of this)
         * Ki                is K_this/2, K_next/2         (note how outer loop takes care of this)
         * dKi/dx            is dKdx_this, dKdx_next       (note how outer loop takes care of this)
         * ji                is nlte_this/2, nlte_next/2   (note how outer loop takes care of this)
         * dji/dx            is dnltedx_this, dnltedx_next (note how outer loop takes care of this)
         * 
         * Begin solution
         */
        
        /*
         * T * dTi/dx * (si-bi-Ki^(-1)ji)
         */
        const Numeric x = -0.5 * r * PiT_next(0,0);
        const Numeric y = x * sibi[0];
        // Stokes 1:
        diydx_this[0] += y * dKdx_this(0,0);
        diydx_next[0] += y * dKdx_next(0,0);
        // Higher Stokes
        for( Index is=1; is<stokes_dim; is++ )
        { 
            const Numeric z = x * sibi[is]; //NB: sibi not iy since ji can be polarized as well?
            diydx_this[is] += z * dKdx_this(is,is);
            diydx_next[is] += z * dKdx_next(is,is);
        }
        
        /*
         * T * (1 - Ti) * dbi/dx
         */
        if(do_Bsource)
        {
            // The terms associated with B-bar:  For extmat_case==1 only works on first term
            const Numeric v = PiT_this(0,0) * ( 1.0 - T_this(0,0));
            diydx_this[0] += 0.5*v * dBdx_this;
            diydx_next[0] += 0.5*v * dBdx_next;
        }
        
        /*
         * T * (1 - Ti) *(K^(-2)dKi/dxji + K^(-1)dji/dx)
         */
        if(do_nonlte)
        {
            const Numeric v = 0.5 * PiT_this(0,0) * ( 1.0 - T_this(0,0));
            
            Numeric invK = 1.0/K_this(0,0);
            
            diydx_this[0] += -v*invK*
            (dKdx_this(0,0)*invK*nlte_this[0] + invK*dnltedx_this[0]);
            
            invK = 1.0/K_next(0,0);
            diydx_next[0] += -v*invK*
            (dKdx_next(0,0)*invK*nlte_next[0] + invK*dnltedx_next[0]);
        }
        
        /*
         * End of solution
         */
        
        // Zero for higher Stokes
        //
        // The terms associated with Delta-s:
        if( do_HSE )
        {
            if(do_nonlte)
                throw std::runtime_error("HSE and non-LTE not yet compatible\n");
                
            // Stokes 1:
            const Numeric kbar = 0.5 * ( K_this(0,0) + K_next(0,0) );
            
            diydx_this[0] += y * kbar / temperature_this;
            diydx_next[0] += y * kbar / temperature_next;
            // Higher Stokes
            for( Index is=1; is<stokes_dim; is++ )
            { 
                const Numeric z = x * iy[is];
                diydx_this[is] += z * kbar / temperature_this;
                diydx_next[is] += z * kbar / temperature_next;
            }
        } //hse
    }
    // General case
    else
    { 
        /*
         * Solves ds/dx = T [ dTi/dx * (si-bi-Ki^(-1)ji) + 
         *        (1 - Ti) * ( dbi/dx + K^(-1)dKi/dxK^(-1)ji + K^(-1)dji/dx ) ]
         * 
         * (si-bi-Ki^(-1)ji) is sibi
         * T                 is PiT_this
         * dTi/dx            is dTdx_this, dTdx_next       (note how outer loop takes care of this)
         * Ti                is T_this
         * dbi/dx            is dBdx_this/2, dBdx_next/2   (note how outer loop takes care of this)
         * Ki                is K_this/2, K_next/2         (note how outer loop takes care of this)
         * dKi/dx            is dKdx_this, dKdx_next       (note how outer loop takes care of this)
         * ji                is nlte_this/2, nlte_next/2   (note how outer loop takes care of this)
         * dji/dx            is dnltedx_this, dnltedx_next (note how outer loop takes care of this)
         * 
         * Begin solution
         */
        
        
        //Helpers
        Vector x(stokes_dim), y(stokes_dim);
        
        /*
         * This is the T * dTi/dx * (si-bi-Ki^(-1)ji) part
         */
        // Disturb for this
        mult( x, dTdx_this, sibi );
        mult( y, PiT_this, x );
        diydx_this += y;
        
        // Disturb for next
        mult( x, dTdx_next, sibi );
        mult( y, PiT_this, x );
        diydx_next += y; 
        
        
        
        // The terms associated with B-bar:
        /*
         * This is the T * (1-Ti)*dbi/dx part
         */
        if(do_Bsource)
        {
            // Do this
            const Numeric v = ( 1.0 - T_this(0,0) );
            x[0] = 0.5*v * dBdx_this;
            for( Index is=1; is<stokes_dim; is++ ) 
            { x[is] = -0.5*T_this(is,0)*dBdx_this; }
            mult( y, PiT_this, x );
            diydx_this += y; 
            
            // Same for next
            x[0] = 0.5*v * dBdx_next;
            for( Index is=1; is<stokes_dim; is++ ) 
            { x[is] = -0.5*T_this(is,0)*dBdx_next; }
            mult( y, PiT_this, x );
            diydx_next += y;
        }
        
        /*
         * This is the  T * (1-Ti)* ( K^(-1)dKi/dxK^(-1)ji + K^(-1)dji/dx ) part
         */
        if(do_nonlte)
        {
            Matrix invK(stokes_dim,stokes_dim),
            mat1(stokes_dim,stokes_dim),
            unit(stokes_dim,stokes_dim);
            Vector vec1(stokes_dim),
            vec2(stokes_dim);
            
            id_mat(unit);
            unit-=T_this;
            
            // for this
            inv(invK,K_this);
            mult(mat1,dKdx_this,invK);
            
            mult(vec1,mat1,nlte_this);
            vec1-=dnltedx_this;
            vec1*=-0.5;
            
            mult(vec2,unit,vec1);
            mult(vec1,PiT_this,vec2);
            
            diydx_this += vec1;
            
            // For next
            inv(invK,K_next);
            mult(mat1,dKdx_next,invK);
            
            mult(vec1,mat1,nlte_next);
            vec1-=dnltedx_next;
            vec1*=-0.5;
            
            mult(vec2,unit,vec1);
            mult(vec1,PiT_this,vec2);
            
            diydx_next += vec1;
        }
        
        /*
         *  End of solution
         * 
         * Below follows hydrostatic equilibrium calculations
         * FIXME:  Add theoretical description
         */
        
        //
        // The terms associated with Delta-s:
        if( do_HSE )
        {
            if(do_nonlte)
                throw std::runtime_error("HSE and non-LTE not yet compatible\n");
            
            Matrix ext_mat(stokes_dim,stokes_dim), dtdx(stokes_dim,stokes_dim);
            for( Index is1=0; is1<stokes_dim; is1++ )
                for( Index is2=0; is2<stokes_dim; is2++ )
                    ext_mat(is1,is2) = 0.5 * ( K_this(is1,is2) + K_next(is1,is2) );
            // dl for disturbed tbar
            const Numeric tbar = 0.5 * ( temperature_this + temperature_next );
            const Numeric dl = r * ( 1.0 + dt/tbar );
            
            ext2trans( dtdx, extmat_case,  ext_mat, dl ); 
            
            for( Index is1=0; is1<stokes_dim; is1++ ) {
                for( Index is2=0; is2<stokes_dim; is2++ ) {
                    dtdx(is1,is2) = (1/dt) * ( dtdx(is1,is2) - T_this(is1,is2) ); }}
                
            mult( x, dtdx, sibi );
            mult( y, PiT_this, x );
            
            // Contribution shared between the two
            // points  and is proportional to 1/t
            // See also AUG.
            for( Index is=0; is<stokes_dim; is++ ) 
            {
                diydx_this[is] += y[is] * 0.5 * tbar / temperature_this;
                diydx_next[is] += y[is] * 0.5 * tbar / temperature_next;
            }
        } // HSE
    } // General case
}
