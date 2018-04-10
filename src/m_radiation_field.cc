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
#include "ppath.h"


// Integrates over the unit-area sphere for a za length over one, but only for evenly spaced aa
void integrate_over_the_sphere(MatrixView iy, ConstTensor4View field, 
                               ConstVectorView za, ConstVectorView aa)
{
  extern const Numeric DEG2RAD;
  const Index nf = iy.nrows(), ns = iy.ncols(), nz = za.nelem(), na = aa.nelem();
  iy = 0.0;
  
  // Area elements of the spherical observation geometry
  Vector dA(nz - 1);
  for(Index iz = 0; iz < nz-1; iz++)
    dA[iz] = DEG2RAD * (za[iz+1]-za[iz]) * sin(DEG2RAD * (za[iz] + 0.5*(za[iz+1] - za[iz])));
  
  // Normalize to make a unit area sphere with dA.sum() times all the various angles
  dA /= dA.sum() * Numeric(na);
  
  for(Index iz = 0; iz < nz-1; iz++)
    for(Index ia = 0; ia < na; ia++)
        for(Index iv = 0; iv < nf; iv++)
          for(Index is = 0; is < ns; is++)
            iy(iv, is) += dA[iz] * (field(iz, ia, iv, is) + field(iz+1, ia, iv, is));
}


void total_line_source_and_transmission(Vector& J, 
                                        Vector& T, 
                                        const GriddedField4& radiation_field, 
                                        const GriddedField4& transmission_field, 
                                        const ArrayOfArrayOfLineRecord& lines, 
                                        const ArrayOfArrayOfSpeciesTag& abs_species,
                                        const ConstVectorView frequency,
                                        const ArrayOfArrayOfIndex& range_frequency,
                                        const ConstVectorView vmrs,
                                        const ConstVectorView wind,
                                        const ConstVectorView za,
                                        const ConstVectorView aa,
                                        const Numeric pressure, 
                                        const Numeric temperature, 
                                        const Verbosity& verbosity)
{
  extern const Numeric SPEED_OF_LIGHT;
  
  const Index nl = lines.nelem(), nz = za.nelem(), na = aa.nelem(), nw = wind.nelem();
  J = Vector(nl, 0.0);
  T = Vector(nl, 0.0);
  
  for(Index il = 0; il < nl; il++) {
    ArrayOfIndex bsl;
    Matrix iy(1, 1, 0.0);
    Matrix it(1, 1, 0.0);
    Tensor4 rad_data(nz, na, 1, 1, 0.0), tra_data(nz, na, 1, 1, 0.0);
    find_broad_spec_locations(bsl, abs_species, lines[il][0].Species());
    const ConstVectorView f = frequency[Range(range_frequency[il][0], range_frequency[il][1])];
    const Index nf = f.nelem();
    ComplexVector F(nf);
    Linefunctions::set_lineshape(F, lines[il][0], f, vmrs, temperature, pressure, 0.0, 0.0, 0.0, bsl, il, -1, verbosity);
    Vector X = F.real();
    
    Numeric integral_of_X = 0;
    for(Index iv = 0; iv < nf-1; iv++)
      integral_of_X += (f[iv+1]-f[iv]) * (X[iv]+X[iv+1]);
    X /= 2.0*integral_of_X; // Should renormalize X to 0.5 incase the user selects extremely strange numbers (and to reduce work in integral further down)
    
    Vector x = X; // Creates a temporary that can be changed by winds
    
    for(Index iz = 0; iz < nz; iz++) {
      for(Index ia = 0; ia < na; ia++) {
        if(nw) {
          Vector f_tmp = f;
          Vector los(2); los[0] = za[iz]; los[1] = aa[ia];
          const Numeric vd = 1. - dotprod_with_los(los, wind[0], wind[1], wind[2], 3)/SPEED_OF_LIGHT;
          f_tmp *= vd;
          ArrayOfGridPos gp(nf);
          Matrix itw(2, nf);
          gridpos(gp, f_tmp, f, 1e99); // WARNING:  This part might fail...  
          interpweights(itw, gp);
          interp(x, itw, X, gp);
        }
        for(Index iv = 0; iv < nf-1; iv++) {
          rad_data(iz, ia, 0, 0) += (f[iv+1]-f[iv]) * (x[iv+1] * radiation_field   .data(iz, ia, iv+1, 0) + x[iv] * radiation_field   .data(iz, ia, iv, 0   ));
          tra_data(iz, ia, 0, 0) += (f[iv+1]-f[iv]) * (x[iv+1] * (1-transmission_field.data(iz, ia, iv+1, 0)) + x[iv] * (1-transmission_field.data(iz, ia, iv, 0)));
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
                                   GriddedField4&          transmission_field,
                                   //IN:
                                   const Index&            atmosphere_dim,
                                   const Index&            cloudbox_on,
                                   const Index&            stokes_dim,
                                   const Numeric&          ppath_lmax,
                                   const Numeric&          ppath_lraytrace,
                                   const Vector&           f_grid,
                                   const Vector&           p_grid,
                                   const Vector&           rte_pos,
                                   const Tensor3&          t_field,
                                   const Tensor3&          z_field,
                                   const Tensor3&          wind_u_field,
                                   const Tensor3&          wind_v_field,
                                   const Tensor3&          wind_w_field,
                                   const Tensor3&          mag_u_field,
                                   const Tensor3&          mag_v_field,
                                   const Tensor3&          mag_w_field,
                                   const Tensor4&          vmr_field,
                                   const Tensor4&          nlte_field,
                                   const ArrayOfArrayOfSpeciesTag& abs_species,
                                   const String&           iy_unit,
                                   const Tensor3&          surface_props_data,
                                   const Agenda&           ppath_agenda,
                                   const Agenda&           iy_main_agenda,
                                   const Agenda&           iy_space_agenda,
                                   const Agenda&           iy_surface_agenda,
                                   const Agenda&           iy_cloudbox_agenda,
                                   const Agenda&           propmat_clearsky_agenda,
                                   const Agenda&           water_psat_agenda,   
                                   //GIN:
                                   const Vector&           za_coords,
                                   const Vector&           aa_coords,
                                   const Index&            for_nlte,
                                   const Verbosity&        verbosity)
{
  const Index nz = za_coords.nelem(), na = aa_coords.nelem(), nf = f_grid.nelem();
  
  // Test input
  if(nz < 2)
    throw std::runtime_error("za_coords must contain at least two elements.\n");
  if(na < 1)
    throw std::runtime_error("aa_coords must contain at least one element.\n");
  if(!(za_coords[0] == 0. && za_coords[nz-1] == 180.))
    throw std::runtime_error("First and last coords in za_coords must be 0 and 180, respectively.\n");
  if(na > 1)
    if(!(aa_coords[0] == -180. && aa_coords[na-1]==180.))
      throw std::runtime_error("First and last coords in aa_coords must be -180 and 180, respectively, if aa_coords.nelem()>1.\n");
    
  Numeric tmp = -1.;
  for(Index iz = 0; iz < nz; iz++) {
    if(za_coords[iz] < tmp)
      throw std::runtime_error("za_coords must be strictly increasing.\n");
    else
      tmp = za_coords[iz];
  }
  
  tmp = -181.;
  for(Index ia = 0; ia < na; ia++) {
    if(aa_coords[ia] < tmp)
      throw std::runtime_error("aa_coords must be strictly increasing.\n");
    else
      tmp = aa_coords[ia];
  }
  
  if(nf == 0)
    throw std::runtime_error("f_grid must contain at least one value.\n");
  
  // Necessary input dummy
  ArrayOfString iy_aux_vars(0);
  
  // Prepare a stokes vector so grid setting is possible
  Vector stokes_dim_vector(stokes_dim);
  for(Index is = 0; is < stokes_dim; is++)
    stokes_dim_vector[is] = Numeric(is);
  
  // Prepare radiation_field
  radiation_field.set_grid(0, za_coords);
  radiation_field.set_grid_name(0, "Zenith Angle");
  radiation_field.set_grid(1, aa_coords);
  radiation_field.set_grid_name(1, "Azimuth Angle");
  radiation_field.set_grid(2, f_grid);
  radiation_field.set_grid_name(2, "Frequency");
  radiation_field.set_grid(3, stokes_dim_vector);
  radiation_field.set_grid_name(3, "Stokes Component");
  radiation_field.data.resize(nz, na, nf, stokes_dim);
  radiation_field.set_name("Radiation Field");
  radiation_field.checksize_strict();
  
  if(for_nlte) {
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
    transmission_field.data.resize(nz, na, nf, stokes_dim);
    transmission_field.set_name("Transmission Field");
    transmission_field.checksize_strict();
  }
  
  const Tensor3 _tmp_transmission(nf, stokes_dim, stokes_dim, 1);
  
  // Loop over coords
  for(Index iz = 0; iz < nz; iz++) {
    for(Index ia = 0; ia < na; ia++) {
      
      // Setup of local los vector from za_coords.
      Vector rte_los(2);
      rte_los[0] = za_coords[iz];
      rte_los[1] = aa_coords[ia];
      
      Ppath ppath;
      ppath_agendaExecute(ws, ppath, ppath_lmax, ppath_lraytrace, rte_pos, rte_los, Vector(0), 
                          cloudbox_on, 0, t_field, z_field, vmr_field, f_grid, ppath_agenda );
      // Necessary output
      ArrayOfMatrix iy_aux;
      ArrayOfTensor3 diy_dx;
      Vector ppvar_t, ppvar_p;
      Matrix ppvar_nlte, ppvar_vmr, ppvar_wind, ppvar_mag, ppvar_f;
      Tensor3 ppvar_iy;
      Tensor4 ppvar_trans_cumulat;
      
      iyEmissionStandard(ws, iy, iy_aux, diy_dx, ppvar_p, ppvar_t,
                         ppvar_nlte, ppvar_vmr, ppvar_wind, ppvar_mag, ppvar_f,
                         ppvar_iy, ppvar_trans_cumulat,  
                         0, stokes_dim, f_grid, atmosphere_dim, p_grid,
                         z_field, t_field, nlte_field, vmr_field, abs_species,
                         wind_u_field, wind_v_field, wind_w_field, mag_u_field, mag_v_field, mag_w_field,
                         cloudbox_on, iy_unit, iy_aux_vars, 0, ArrayOfRetrievalQuantity(0),
                         ppath, Vector(0), propmat_clearsky_agenda, water_psat_agenda,
                         iy_main_agenda, iy_space_agenda,
                         iy_surface_agenda, iy_cloudbox_agenda, 0, _tmp_transmission,
                         0.0, surface_props_data, verbosity);
      
      // Add to radiation field
      radiation_field.data(iz, ia, joker, joker) = iy;
      
      if(for_nlte) {
        transmission_field.data(iz, ia, joker, joker) = iy_aux[0];
      }
    }
  }
  
  // Get iy from radiation_field
  if(not for_nlte) {
    iy.resize(f_grid.nelem(),stokes_dim);
    iy=0;
    integrate_over_the_sphere(iy, radiation_field.data, za_coords, aa_coords);
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
                                          Matrix&                         iy,
                                          Tensor3&                        iy_transmission,
                                          const ArrayOfArrayOfSpeciesTag& abs_species,
                                          const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                          const Tensor4&                  nlte_field,
                                          const Tensor4&                  vmr_field,
                                          const Tensor3&                  t_field,
                                          const Tensor3&                  z_field,
                                          const Tensor3&                  wind_u_field,
                                          const Tensor3&                  wind_v_field,
                                          const Tensor3&                  wind_w_field,
                                          const Vector&                   p_grid,
                                          const Index&                    atmosphere_dim,
                                          const Tensor3&                  surface_props_data,
                                          const Agenda&                   ppath_agenda,
                                          const Agenda&                   iy_main_agenda,
                                          const Agenda&                   iy_space_agenda,
                                          const Agenda&                   iy_surface_agenda,
                                          const Agenda&                   iy_cloudbox_agenda,
                                          const Agenda&                   propmat_clearsky_agenda,
                                          const Agenda&                   water_psat_agenda,   
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
  Vector f_grid(nf*nl), abs_t = t_field(joker, 0, 0);
  ArrayOfArrayOfIndex pos_f_grid(nl, ArrayOfIndex(2));
  Matrix abs_vmrs(nl, np), abs_nlte;
  if(nlte_field.nbooks() and nlte_field.npages() and nlte_field.nrows() and nlte_field.ncols())
    abs_nlte = nlte_field(joker, joker, 0, 0);
  ArrayOfIndex aoi(nl);
  const ArrayOfArrayOfSpeciesTag aoaost(nl, abs_species[0]);
  ArrayOfArrayOfLineRecord aoaolr(nl, ArrayOfLineRecord(1));
  Index itot = 0;
  for(Index il = 0; il < nl; il++) {
    aoaolr[il][0] = abs_lines_per_species[0][il];
    aoi[il] = il;
    Vector tmp_f(nf);
    nlinspace(tmp_f, aoaolr[il][0].F() - df, aoaolr[il][0].F() + df, nf);
    for(Index iv = 0; iv < nf; iv++) {
      f_grid[itot] = tmp_f[iv];
      itot++;
    }
    pos_f_grid[il][0] = il*nf;
    pos_f_grid[il][1] = nf;
    abs_vmrs(il, joker) = vmr_field(0, joker, 0, 0);
  }
  
  Matrix ppath_wind;
  if(use_wind) {
    ppath_wind.resize(3, np);
    for(Index ip = 0; ip < np; ip++) {
      ppath_wind(0, ip) = wind_u_field(ip, 0, 0);
      ppath_wind(1, ip) = wind_v_field(ip, 0, 0);
      ppath_wind(2, ip) = wind_w_field(ip, 0, 0);
    }
  }
  
  iy = Matrix(nl, np);
  iy_transmission = Tensor3(1, nl, np);
  
  // Compute the field with the slow method
//   #pragma omp parallel for
  for(Index ip = 0; ip < np; ip++) {
    GriddedField4 rady;
    GriddedField4 tran;
    Vector rte_pos(3, 0.); 
    rte_pos[0] = z_field(ip, 0, 0) + 0.01;  // The value at 1 cm above because of bug in ppath calculation
    Matrix _tmp;
    
    radiation_fieldCalcFromiyCalc(ws, _tmp, rady, tran, atmosphere_dim, 0, 1, 1e10, 1e10, f_grid, p_grid,
                                  rte_pos, t_field, z_field, wind_u_field, wind_v_field, wind_w_field,
                                  Tensor3(0, 0, 0), Tensor3(0, 0, 0), Tensor3(0, 0, 0), vmr_field, nlte_field,
                                  abs_species, "1", surface_props_data,
                                  ppath_agenda, iy_main_agenda, iy_space_agenda,
                                  iy_surface_agenda, iy_cloudbox_agenda,
                                  propmat_clearsky_agenda, water_psat_agenda,
                                  za, aa, ip+1, verbosity);
    
    Vector J(nl, 0), G(nl, 0);
    total_line_source_and_transmission(J, G, rady, tran, 
                                       aoaolr, abs_species, f_grid, pos_f_grid, 
                                       abs_vmrs(joker, ip), use_wind?ppath_wind(joker, ip):Vector(0), 
                                       za, aa, p_grid[ip], abs_t[ip], verbosity);
    
    iy(joker, ip) = J;
    iy_transmission(0, joker, ip) = G;
  }
}




/* Workspace method: Doxygen documentation will be auto-generated */
void doit_i_fieldClearskyPlaneParallel(
        Workspace&                  ws,
        Tensor7&                    doit_i_field,
        Tensor3&                    trans_field,
  const Agenda&                     propmat_clearsky_agenda,
  const Agenda&                     water_psat_agenda,   
  const Agenda&                     iy_space_agenda,
  const Agenda&                     iy_surface_agenda, 
  const Agenda&                     iy_cloudbox_agenda,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Tensor3&                    z_field,
  const Tensor3&                    t_field,
  const Tensor4&                    nlte_field,
  const Tensor4&                    vmr_field,
  const ArrayOfArrayOfSpeciesTag&   abs_species,
  const Tensor3&                    wind_u_field,
  const Tensor3&                    wind_v_field,
  const Tensor3&                    wind_w_field,
  const Tensor3&                    mag_u_field,
  const Tensor3&                    mag_v_field,
  const Tensor3&                    mag_w_field,
  const Matrix&                     z_surface,
  const Numeric&                    ppath_lmax,
  const Numeric&                    rte_alonglos_v,
  const Tensor3&                    surface_props_data,
  const Vector&                     scat_za_grid, 
  const Verbosity&                  verbosity )
{
  // Check input
  if( atmosphere_dim != 1 )
    throw runtime_error( "This method only works for atmosphere_dim = 1." );
  chk_if_increasing ( "scat_za_grid", scat_za_grid );
  if( abs( z_surface(0,0) - z_field(0,0,0) ) > 1e-3 )
    throw runtime_error(
        "The surface must be placed exactly at the first pressure level." );
  
  // Sizes
  const Index nl  = p_grid.nelem();
  const Index nf  = f_grid.nelem();
  const Index nza = scat_za_grid.nelem();
  
  // Init doit_i_field and trans_field
  doit_i_field.resize( nf, nl, 1, 1, nza, 1, stokes_dim );
  trans_field.resize( nf, nl, nza );
  
  // De-activate cloudbox 
  const Index cloudbox_on = 0, ppath_inside_cloudbox_do = 0;
  const ArrayOfIndex cloudbox_limits( 0 );

  // Various input variables
  const String                     iy_unit = "1";
  const ArrayOfString              iy_aux_vars(0);
  const Vector                     rte_pos2(0);
  const Index                      iy_agenda_call1 = 1;
  const Tensor3                    iy_transmission(0,0,0);
  const Index                      jacobian_do = 0;
  const ArrayOfRetrievalQuantity   jacobian_quantities(0);
  // Create one altitde just above TOA
  const Numeric z_space = z_field(nl-1,0,0) + 10;

  Workspace l_ws(ws);
  ArrayOfString fail_msg;
  bool failed = false;

  // Define iy_main_agenda to be consistent with the methods applied inside
  // this method. This definition of iy_main_agenda will be used to when
  // calculating the the radiation reflected by the surface
  Agenda iy_main_agenda;
  iy_main_agenda.append("ppathPlaneParallel", TokVal());
  iy_main_agenda.append("iyEmissionStandard", TokVal());
  iy_main_agenda.set_name("iy_main_agenda");
  iy_main_agenda.check(ws, verbosity);
  
  // Loop zenith angles
  //
  if (nza)
#pragma omp parallel for       \
  if (!arts_omp_in_parallel()  \
      && nza > 1)         \
  firstprivate(l_ws)
  for( Index i=0; i<nza; i++ )
    {
      if (failed)
        continue;
      try
        {
          // Define output variables
          Ppath          ppath;
          Vector         ppvar_p, ppvar_t;
          Matrix         iy, ppvar_nlte, ppvar_vmr, ppvar_wind, ppvar_mag, ppvar_f;
          Tensor3        ppvar_iy;
          Tensor4        ppvar_trans_cumulat;
          ArrayOfMatrix  iy_aux;
          ArrayOfTensor3 diy_dx;

          Index  iy_id = i;
          Vector rte_los( 1, scat_za_grid[i] );
          Vector rte_pos( 1, scat_za_grid[i] < 90 ? z_surface(0,0) : z_space );

          ppathPlaneParallel( ppath, atmosphere_dim, z_field, z_surface,cloudbox_on,
                              cloudbox_limits, ppath_inside_cloudbox_do,
                              rte_pos, rte_los, ppath_lmax, verbosity );

          iyEmissionStandard( l_ws, iy, iy_aux, diy_dx, ppvar_p, ppvar_t,ppvar_nlte,
                              ppvar_vmr, ppvar_wind, ppvar_mag, ppvar_f, ppvar_iy,
                              ppvar_trans_cumulat,iy_id, stokes_dim, f_grid,
                              atmosphere_dim,p_grid, z_field, t_field, nlte_field,
                              vmr_field, abs_species,wind_u_field, wind_v_field,
                              wind_w_field, mag_u_field, mag_v_field, mag_w_field,
                              cloudbox_on, iy_unit, iy_aux_vars, jacobian_do,
                              jacobian_quantities, ppath, rte_pos2,
                              propmat_clearsky_agenda, water_psat_agenda,
                              iy_main_agenda, iy_space_agenda,
                              iy_surface_agenda, iy_cloudbox_agenda,
                              iy_agenda_call1, iy_transmission, rte_alonglos_v,
                              surface_props_data, verbosity );
          assert( iy.nrows() == nf );
          assert( iy.ncols() == stokes_dim );

          // First and last points are most easily handled separately
          if( scat_za_grid[i] < 90 )
            {
              doit_i_field(joker,0,0,0,i,0,joker)    = ppvar_iy(joker,joker,0);
              doit_i_field(joker,nl-1,0,0,i,0,joker) = ppvar_iy(joker,joker,ppath.np-1);
              trans_field(joker,0,i)    = ppvar_trans_cumulat(0,joker,0,0);
              trans_field(joker,nl-1,i) = ppvar_trans_cumulat(ppath.np-1,joker,0,0);
            }
          else
            {
              doit_i_field(joker,nl-1,0,0,i,0,joker) = ppvar_iy(joker,joker,0);
              doit_i_field(joker,0,0,0,i,0,joker)    = ppvar_iy(joker,joker,ppath.np-1);
              trans_field(joker,nl-1,i) = ppvar_trans_cumulat(0,joker,0,0);
              trans_field(joker,0,i)    = ppvar_trans_cumulat(ppath.np-1,joker,0,0);
            }

          // Remaining points
          for( Index p=1; p<ppath.np-1; p++ )
            {
              // We just store values at pressure levels
              if( ppath.gp_p[p].fd[0] < 1e-2 )
                {
                  doit_i_field(joker,ppath.gp_p[p].idx,0,0,i,0,joker) =
                    ppvar_iy(joker,joker,p);
                  trans_field(joker,ppath.gp_p[p].idx,i) = ppvar_trans_cumulat(p,joker,0,0);
                }
            }
        }
      catch (runtime_error e)
        {
#pragma omp critical (planep_setabort)
          failed = true;

          ostringstream os;
          os << "Run-time error at nza #" << i << ": \n" << e.what();
#pragma omp critical (planep_push_fail_msg)
          fail_msg.push_back(os.str());
        }
    }

  if (fail_msg.nelem())
  {
    ostringstream os;
    for (auto& msg : fail_msg)
      os << msg << '\n';
    throw runtime_error(os.str());
  }
}
