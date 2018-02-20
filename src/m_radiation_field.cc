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
#include "physics_funcs.h"


void integrate_over_the_sphere(MatrixView iy, const GriddedField4& radiation_field, 
                               ConstVectorView za, ConstVectorView aa, 
                               ConstTensor3View weights=Tensor3(0, 0, 0)) noexcept
{
  extern const Numeric DEG2RAD, PI;
  const Index nf = iy.nrows(), ns = iy.ncols(), nz = za.nelem(), na = aa.nelem();
  const Index nwnz = weights.npages(), nwna = weights.nrows(), nwnf = weights.ncols();
  
  const bool use_weights = nwnf or nwnz or nwna;
  
  // Area elements of the spherical observation geometry
  Matrix dA(nz - 1, na);
  for(Index iz = 0; iz < nz-1; iz++)
  {
    for(Index ia = 0; ia<na-1; ia++)
    {
      dA(iz, ia) = DEG2RAD * (za[iz+1]-za[iz]) * 
               sin(DEG2RAD * (za[iz] + (za[iz+1] - za[iz])/2.)) * 
                   DEG2RAD * (aa[ia+1]-aa[ia]);
    }
    if(na == 1)
    {
      dA(iz, 0) = DEG2RAD * (za[iz+1]-za[iz]) *
              sin(DEG2RAD * (za[iz] + (za[iz+1] - za[iz])/2.)) *
                  DEG2RAD * (360.);
    }
  }
  
  
  
  // Sum the contributions by relative area of the sphere
  for(Index iz = 0; iz < nz - 1; iz++)
  {
    if(na > 1)
    {
      for(Index ia = 0; ia < na-1; ia++)
      {
        Matrix iy_tmp(nf, ns);
        if(use_weights)
        {
          for(Index iv = 0; iv < nf; iv++)
          {
            iy_tmp(iv, joker)  = radiation_field.data(iz,   ia,   iv, joker) * weights(iz,   ia,   iv);
            iy_tmp(iv, joker) += radiation_field.data(iz+1, ia,   iv, joker) * weights(iz+1, ia,   iv);
            iy_tmp(iv, joker) += radiation_field.data(iz,   ia+1, iv, joker) * weights(iz,   ia+1, iv);
            iy_tmp(iv, joker) += radiation_field.data(iz+1, ia+1, iv, joker) * weights(iz+1, ia+1, iv);
          }
        }
        else
        {
          iy_tmp  = radiation_field.data(iz,   ia,   joker, joker);
          iy_tmp += radiation_field.data(iz+1, ia,   joker, joker);
          iy_tmp += radiation_field.data(iz,   ia+1, joker, joker);
          iy_tmp += radiation_field.data(iz+1, ia+1, joker, joker);
        }
        iy_tmp *= dA(iz, ia) / 4.0 ;
        iy += iy_tmp;
      }
    }
    else
    {
      Matrix iy_tmp(nf, ns);
      if(use_weights)
      {
        for(Index iv = 0; iv < nf; iv++)
        {
          iy_tmp(iv, joker)  = radiation_field.data(iz,   0, iv, joker) * weights(iz,   0, iv);
          iy_tmp(iv, joker) += radiation_field.data(iz+1, 0, iv, joker) * weights(iz+1, 0, iv);
        }
      }
      else
      {
        iy_tmp  = radiation_field.data(iz,   0, joker, joker);
        iy_tmp += radiation_field.data(iz+1, 0, joker, joker);
      }
      iy_tmp *= dA(iz,0) / 2.0 ;
      iy+=iy_tmp;
    }
  }
  
  // Normalize for surface area of the sphere
  iy /= 4*PI;
}


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
                   atmfields_checked, atmgeom_checked, iy_aux_vars, 0,
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
    integrate_over_the_sphere(iy, radiation_field, za_coords, aa_coords);
}


bool compare_frequency(const Numeric& a, const Numeric& b) {return a < b;}


void VectorSort(VectorView V) noexcept
{
  // Fix if there were lines overlapping... (FIXME:  std::sort not working...)
  //   std::sort(f_grid.begin(), f_grid.end(), compare_frequency);
  
  // Gnome sort from wikipedia...  (aka, stupid sort...)
  const Index len = V.nelem();
  Index pos = 0;
  while(pos < len)
  {
    if(pos == 0)
      pos++;
    else if(V[pos] >= V[pos-1])
      pos++;
    else 
    {
      swap(V[pos], V[pos-1]);
      pos--;
    }
  }
}


void ppath_windSet(Matrix& ppath_wind, const Ppath& ppath, 
                   ConstTensor3View wind_u_field, ConstTensor3View wind_v_field, ConstTensor3View wind_w_field)
{
  // Winds:  COPY FROM get_ppath_atmvars
  ppath_wind.resize(3, ppath.np);
  ppath_wind = 0;
  Matrix itw_field;
  interp_atmfield_gp2itw(itw_field, 3, ppath.gp_p, ppath.gp_lat, ppath.gp_lon);
  if( wind_u_field.npages() > 0 ) 
    interp_atmfield_by_itw(ppath_wind(0,joker), 3, wind_u_field, ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field);
  if( wind_v_field.npages() > 0 )
    interp_atmfield_by_itw(ppath_wind(1,joker), 3, wind_v_field, ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field);
  if( wind_w_field.npages() > 0 )
    interp_atmfield_by_itw(ppath_wind(2,joker), 3, wind_w_field, ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field);
}


void radiation_fieldCalcForRotationalNLTE(Workspace&                      ws,
                                          Vector&                         y,
                                          const ArrayOfArrayOfSpeciesTag& abs_species,
                                          const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                          const SpeciesAuxData&           isotopologue_ratios,
                                          const SpeciesAuxData&           partition_functions,
                                          const Agenda&                   iy_main_agenda,
                                          const Tensor4&                  nlte_field,
                                          const Tensor4&                  vmr_field,
                                          const Tensor3&                  t_field,
                                          const Tensor3&                  z_field,
                                          const Tensor3&                  wind_u_field,
                                          const Tensor3&                  wind_v_field,
                                          const Tensor3&                  wind_w_field,
                                          const Vector&                   p_grid,
                                          const Numeric&                  df,
                                          const Index&                    nz,
                                          const Index&                    na,
                                          const Index&                    nf,
                                          const Verbosity&                verbosity)
{
  extern const Numeric SPEED_OF_LIGHT;
  bool use_wind=false;
  
  // Number of zenith angles and a zheck that there is enough to make simple integrations
  if(nz < 3)
    throw std::runtime_error("Must have at least three zenith angles");
  Vector za;
  nlinspace(za, 0., 180., nz);
  
  Vector aa;
  if(na < 1)
    throw std::runtime_error("Need at least one azimuth angle");
  else if(na == 1)
    aa = Vector(1, 0.0);
  else
    nlinspace(aa, -180.0, 180., na);
  
  // Number of species and a check that there is only 1 in all input
  const Index ns = abs_species.nelem();
  if(ns not_eq 1) 
    throw std::runtime_error("Only one species allowed in the present implementation");
  if(abs_species[0].nelem() not_eq 1)
    throw std::runtime_error("Only one species tag allowed in the present implementation");
  if(abs_lines_per_species.nelem() not_eq ns or vmr_field.nbooks() not_eq ns) 
    throw std::runtime_error("Bad number of species to absorption lines.");
  
  // Number of pressure grid-points and a check that the atmosphere is 1D
  const Index np = p_grid.nelem();
  if(t_field.npages() not_eq np or z_field.npages() not_eq np or vmr_field.npages() not_eq np)  
    throw std::runtime_error("Pressure grid is bad");
  
  // Number of lines and a check that there are lines
  const Index nl = abs_lines_per_species[0].nelem();
  if(not nl)
    throw std::runtime_error("Need atleast one line");
  
  if(nf < 2)
    throw std::runtime_error("Must have atleast two frequency points...");
  
  // Initialization of internal variables that reverts the order of the lines and species
  ArrayOfMatrix pconst_xsec(nl, Matrix(nl*nf, np, 0.0)), pconst_xsrc(nl, Matrix(nl*nf, np, 0.0));
  Vector f_grid(nl*nf), abs_t = t_field(joker, 0, 0);
  Matrix abs_vmrs(nl, np), abs_nlte;
  if(nlte_field.nbooks() and nlte_field.npages() and nlte_field.nrows() and nlte_field.ncols())
    abs_nlte = nlte_field(joker, joker, 0, 0);
  ArrayOfIndex aoi(nl);
  const ArrayOfArrayOfSpeciesTag aoaost(nl, abs_species[0]);
  ArrayOfArrayOfLineRecord aoaolr(nl, ArrayOfLineRecord(1));
  for(Index il = 0; il < nl; il++)
  {
    aoaolr[il][0] = abs_lines_per_species[0][il];
    aoi[il] = il;
    
    nlinspace(f_grid[Range((nl-1)*nf, nf)], aoaolr[il][0].F() - df, aoaolr[il][0].F() + df, nf);
    abs_vmrs(il, joker) = vmr_field(0, joker, 0, 0);
  }
  
  VectorSort(f_grid);
  
  // Set dummy variables 
  ArrayOfArrayOfMatrix dxsec(0), dxsrc(0);
  const static ArrayOfRetrievalQuantity jacs(0);
  const static Numeric lm = 0.;
  const static Index xsp = 0;
  
  Matrix ppath_wind, ppath_f;
  
  // Perform computations of cross-section and normalized source function...  Note that this is without wind at this point
  abs_xsec_per_speciesAddLines2(pconst_xsec, pconst_xsrc, dxsec, dxsrc, aoaost, jacs, aoi, f_grid, p_grid, abs_t, abs_nlte, 
                                lm, xsp, abs_vmrs, aoaolr, isotopologue_ratios, partition_functions, verbosity);
  
  // Compute the field with the slow method
  ArrayOfGriddedField4 data(np);
  for(Index ip = 0; ip < np; ip++)
  {
    Vector rte_pos(3, 0.); 
    rte_pos[0] = z_field(ip, 0, 0);
    if(ip == np - 1  or ip == 0)
      rte_pos[0] += 0.1;
    Matrix iy;
    radiation_fieldCalcFromiyCalc(ws, iy, data[ip], 1, 1, f_grid, t_field, z_field, vmr_field,
                                  0, 1, 1, rte_pos, "1", iy_main_agenda, za, aa, verbosity);
  }
  
  // Generate a weighting tensor 
  Tensor5 weights(nl, np, nz, na, nf*nl);
  for(Index iz = 0; iz < nz; iz++) 
  {
    for(Index ia = 0; ia < na; ia++)
    {
      Vector los(2); los[0] = za[iz]; los[1] = aa[ia];
      for(Index ip = 0; ip < np; ip++)
      {
        if(use_wind)
        {
          const Numeric vd = 1. - dotprod_with_los(los, wind_u_field(ip, 0, 0), wind_v_field(ip, 0, 0), 
                                                  wind_w_field(ip, 0, 0), 3)/SPEED_OF_LIGHT;
          Vector tmp(f_grid);
          tmp *= vd;
          ArrayOfGridPos gp(nf*nl);
          Matrix itw(2, nf*nl);
          ArrayOfMatrix xsec(nl, Matrix(nl*nf, np));
          gridpos(gp, f_grid, tmp, 0.5); // WARNING:  This part might fail...
          interpweights(itw, gp);
          for(Index il = 0; il < nl; il++)
          {
            interp(xsec[il](joker, ip), itw, pconst_xsec[il](joker, ip), gp);
            xsec[il](joker, ip) /= xsec[il](joker, ip).sum();
            weights(il, ip, iz, ia, joker) = xsec[il](joker, ip);
          }
        }
        else
        {
          for(Index il = 0; il < nl; il++)
          {
            weights(il, ip, iz, ia, joker) = pconst_xsec[il](joker, ip);
            weights(il, ip, iz, ia, joker) /= pconst_xsec[il](joker, ip).sum();
          }
        }
      }
    }
  }
  
  Matrix iy(nf*nl, 1);
  y = Vector(np*nl, 0.0);
  Index i = 0;
  for(Index ip = 0; ip < np; ip++)
  {
    for(Index il = 0; il < nl; il++)
    {
      integrate_over_the_sphere(iy, data[ip],  za, aa, weights(il, ip, joker, joker, joker));
      for(Index iv = 0; iv < nl*nf; iv++)
        y[i] += iy(iv, 0);
      i++;
    }
  }
}
