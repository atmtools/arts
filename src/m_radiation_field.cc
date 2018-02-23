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
#include "absorption.h"
#include "linefunctions.h"


// Integrates over the unit-area sphere for a za length over one, but only for evenly spaced aa
void integrate_over_the_sphere(MatrixView iy, ConstTensor4View radiation_field, 
                               ConstVectorView za, ConstVectorView aa) noexcept
{
  extern const Numeric DEG2RAD, PI;
  const Index nf = iy.nrows(), ns = iy.ncols(), nz = za.nelem(), na = aa.nelem();
  iy = 0.0;
  
  // Area elements of the spherical observation geometry
  Vector dA(nz - 1);
  for(Index iz = 0; iz < nz-1; iz++)
      dA[iz] = DEG2RAD * (za[iz+1]-za[iz]) *  sin(DEG2RAD * (za[iz] + 0.5*(za[iz+1] - za[iz])));
  
  // Normalize over the zenith, azimuth, unit sphere, and averaging of two points
  dA /= dA.sum() * Numeric(na) * 4.0 * PI * 2.0;
  
  for(Index iz = 0; iz < nz-1; iz++)
    for(Index ia = 0; ia < na; ia++)
        for(Index iv = 0; iv < nf; iv++)
          for(Index is = 0; is < ns; is++)
            iy(iv, is) += dA[iz] * (radiation_field(iz, ia, iv, is) + radiation_field(iz+1, ia, iv, is));
}


void total_line_source_and_transmission(Vector& J, 
                                        Vector& T, 
                                        const GriddedField4& radiation_field, 
                                        const GriddedField5& transmission_field, 
                                        const ArrayOfArrayOfLineRecord& lines, 
                                        const ArrayOfArrayOfSpeciesTag& abs_species,
                                        const ConstVectorView frequency,
                                        const ConstVectorView vmrs,
                                        const ConstVectorView wind,
                                        const ConstVectorView za,
                                        const ConstVectorView aa,
                                        const Numeric pressure, 
                                        const Numeric temperature, 
                                        const Verbosity& verbosity)
{
  extern const Numeric SPEED_OF_LIGHT;
  
  const Index nl = lines.nelem(), nf = frequency.nelem(), nz = za.nelem(), na = aa.nelem(), nw = wind.nelem();
  J = Vector(nl, 0.0);
  T = Vector(nl, 0.0);
  
  for(Index il = 0; il < nl; il++)
  {
    ArrayOfIndex bsl;
    Matrix iy(1, 1, 0.0);
    Matrix it(1, 1, 0.0);
    Tensor4 rad_data(nz, na, 1, 1, 0.0), tra_data(nz, na, 1, 1, 0.0);
    find_broad_spec_locations(bsl, abs_species, il);
    
    ComplexVector F(nf);
    Linefunctions::set_lineshape(F, lines[il][0], frequency, vmrs, temperature, pressure, 0.0, 0.0, 0.0, bsl, il, -1, verbosity);
    const Vector X = F.real();
    Vector x = X;
    
    for(Index iz = 0; iz < nz; iz++)
    {
      for(Index ia = 0; ia < na; ia++)
      {
        if(nw)
        {
          Vector f = frequency;
          Vector los(2); los[0] = za[iz]; los[1] = aa[ia];
          const Numeric vd = 1. - dotprod_with_los(los, wind[0], wind[1], wind[2], 3)/SPEED_OF_LIGHT;
          f *= vd;
          ArrayOfGridPos gp(nf);
          Matrix itw(2, nf);
          gridpos(gp, frequency, f, 0.5); // WARNING:  This part might fail...
          interpweights(itw, gp);
          interp(x, itw, X, gp);
        }
        
        for(Index iv = 0; iv < nf-1; iv++)
        {
          rad_data(iz, ia, 0, 0) += (frequency[iv+1]-frequency[iv]) * (x[iv+1] * radiation_field.data(iz, ia, iv+1, 0) + x[iv]  * radiation_field.data(iz, ia, iv, 0)) * 0.5;
          
          // FIXME:  How to compute factor r?
          tra_data(iz, ia, 0, 0) += (frequency[iv+1]-frequency[iv]) *  (x[iv+1] * transmission_field.data(iz, ia, iv+1, 0, 0) + x[iv] * transmission_field.data(iz, ia, iv, 0, 0)) * 0.5;
        }
      }
    }
    
    integrate_over_the_sphere(iy, rad_data, za, aa);
    J[il] = iy(0, 0);
    integrate_over_the_sphere(it, tra_data, za, aa);
    T[il] = it(0, 0);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void radiation_fieldCalcFromiyCalc(Workspace&              ws,
                                   //OUT:
                                   Matrix&                 iy,
                                   GriddedField4&          radiation_field,
                                   GriddedField5&          transmission_field,
                                   //IN:
                                   const Index&            atmfields_checked,
                                   const Index&            atmgeom_checked,
                                   const Vector&           f_grid,
                                   const Tensor3&          t_field,
                                   const Tensor3&          z_field,
                                   const Tensor4&          vmr_field,
                                   const Tensor4&          nlte_field,
                                   const Index&            cloudbox_on,
                                   const Index&            cloudbox_checked,
                                   const Index&            scat_data_checked,
                                   const Index&            stokes_dim,
                                   const Vector&           rte_pos,
                                   const String&           iy_unit,  
                                   const Agenda&           iy_main_agenda,
                                   //GIN:
                                   const Vector&           za_coords,
                                   const Vector&           aa_coords,
                                   const Index&            for_nlte,
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
    ArrayOfString    iy_aux_vars(0);

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
    
    if(for_nlte)
    {
      iy_aux_vars.resize(1);
      iy_aux_vars[0] = "Transmission";
      transmission_field.set_grid(0, za_coords);
      transmission_field.set_grid_name(0, "Zenith Angle");
      transmission_field.set_grid(1, aa_coords);
      transmission_field.set_grid_name(1, "Azimuth Angle");
      transmission_field.set_grid(2, f_grid);
      transmission_field.set_grid_name(2, "Frequency");
      transmission_field.set_grid(3, stokes_dim_vector);
      transmission_field.set_grid_name(3, "Stokes Component");
      transmission_field.set_grid(4, stokes_dim_vector);
      transmission_field.set_grid_name(4, "Stokes Component");
      transmission_field.data.resize(za_coords.nelem(),aa_coords.nelem(),
                                     f_grid.nelem(),stokes_dim,stokes_dim);
      transmission_field.set_name("Transmission Field");
      transmission_field.checksize_strict();
    }

    // Loop over coords
    for(Index ii=0; ii<za_coords.nelem();ii++)
    {
        for(Index jj=0; jj<aa_coords.nelem();jj++)
        {
            // Necessary output dummys
            ArrayOfTensor4   iy_aux;
            Ppath            ppath_dummy;

            // Setup of local los vector from za_coords.
            Vector rte_los(2);
            rte_los[0]=za_coords[ii];
            rte_los[1]=aa_coords[jj];

            // Calculate iy for a particular rte_los.
            iyCalc(ws, iy, iy_aux, ppath_dummy,
                   atmfields_checked, atmgeom_checked, iy_aux_vars, 0,
                   f_grid, t_field, z_field, vmr_field, nlte_field,
                   cloudbox_on, cloudbox_checked, scat_data_checked,
                   rte_pos, rte_los, Vector(0),
                   iy_unit, iy_main_agenda, verbosity);

            // Add to radiation field
            radiation_field.data(ii,jj,joker,joker) = iy;
            
            if(for_nlte)
              transmission_field.data(ii, jj, joker, joker, joker) = iy_aux[0](joker, joker, joker, iy_aux[0].ncols()-1);
        }
    }

    // Get iy from radiation_field
    if(not for_nlte)
    {
      iy.resize(f_grid.nelem(),stokes_dim);
      iy=0;
      integrate_over_the_sphere(iy, radiation_field.data, za_coords, aa_coords);
    }
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
  
  Matrix ppath_wind;
  if(use_wind)
  {
    ppath_wind.resize(3, np);
    for(Index ip = 0; ip < np; ip++)
    {
      ppath_wind(0, ip) = wind_u_field(ip, 0, 0);
      ppath_wind(1, ip) = wind_v_field(ip, 0, 0);
      ppath_wind(2, ip) = wind_w_field(ip, 0, 0);
    }
  }
  
  Vector t = Vector(nl*np);
  y = Vector(nl*np);
  Index i = 0;
  
  // Compute the field with the slow method
  GriddedField4 rady;
  GriddedField5 tran;
  for(Index ip = 0; ip < np; ip++)
  {
    Vector rte_pos(3, 0.); 
    rte_pos[0] = z_field(ip, 0, 0) + 0.01;  // The value at 1 cm above because of bug in ppath calculation
    Matrix iy;
    radiation_fieldCalcFromiyCalc(ws, iy, rady, tran, 1, 1, f_grid, t_field, z_field, vmr_field, nlte_field,
                                  0, 1, 0, 1, rte_pos, "1", iy_main_agenda, za, aa, 1, verbosity);
    
    Vector J, T;
    
    total_line_source_and_transmission(J, T, rady, tran, 
                                       aoaolr, abs_species, f_grid, abs_vmrs(joker, ip),
                                       use_wind?ppath_wind(joker, ip):Vector(0), za, aa, 
                                       p_grid[ip], abs_t[ip], verbosity);
    for(Index il = 0; il < nl; il++)
    {
      y[i] = J[il];
      t[i] = T[il];
      i++;
    }
  }
}
