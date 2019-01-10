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
                                        const Numeric temperature)
{
  extern const Numeric SPEED_OF_LIGHT;
  
  const Index nl = lines.nelem(), nz = za.nelem(), na = aa.nelem(), nw = wind.nelem();
  J = Vector(nl, 0.0);
  T = Vector(nl, 0.0);
  
  for(Index il = 0; il < nl; il++) {
    Matrix iy(1, 1, 0.0);
    Matrix it(1, 1, 0.0);
    Tensor4 rad_data(nz, na, 1, 1, 0.0), tra_data(nz, na, 1, 1, 0.0);
    const ConstVectorView f = frequency[Range(range_frequency[il][0], range_frequency[il][1])];
    const Index nf = f.nelem();
    Eigen::VectorXcd F(nf);
    Linefunctions::set_lineshape(F, MapToEigen(f), lines[il][0], vmrs, temperature, pressure, 0.0, abs_species, il, 0);
    Vector X(F.size()); for(int iv=0; iv<F.size(); iv++) X[iv] = F.real()[iv];
    
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
                                   const Agenda&           water_p_eq_agenda,   
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
    iy_aux_vars[0] = "Optical depth";
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
      Tensor4 ppvar_trans_cumulat, dummy;
      
      iyEmissionStandard(ws, iy, iy_aux, diy_dx, ppvar_p, ppvar_t,
                         ppvar_nlte, ppvar_vmr, ppvar_wind, ppvar_mag, ppvar_f,
                         ppvar_iy, ppvar_trans_cumulat,  dummy,
                         0, stokes_dim, f_grid, atmosphere_dim, p_grid,
                         z_field, t_field, nlte_field, vmr_field, abs_species,
                         wind_u_field, wind_v_field, wind_w_field, mag_u_field, mag_v_field, mag_w_field,
                         cloudbox_on, iy_unit, iy_aux_vars, 0, ArrayOfRetrievalQuantity(0),
                         ppath, Vector(0), propmat_clearsky_agenda, water_p_eq_agenda,
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


void radiation_fieldCalcForSingleSpeciesNonOverlappingLines(Workspace&                      ws,
                                                            Matrix&                         iy,
                                                            Tensor3&                        iy_transmission,
                                                            const ArrayOfArrayOfSpeciesTag& abs_species,
                                                            const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                                            const Tensor4&                  nlte_field,
                                                            const Tensor4&                  vmr_field,
                                                            const Tensor3&                  t_field,
                                                            const Tensor3&                  z_field,
                                                            const Vector&                   p_grid,
                                                            const Index&                    atmosphere_dim,
                                                            const Tensor3&                  surface_props_data,
                                                            const Agenda&                   iy_space_agenda,
                                                            const Agenda&                   iy_surface_agenda,
                                                            const Agenda&                   iy_cloudbox_agenda,
                                                            const Agenda&                   propmat_clearsky_agenda,
                                                            const Agenda&                   water_p_eq_agenda,   
                                                            const Numeric&                  df,
                                                            const Index&                    nz,
                                                            const Index&                    nf,
                                                            const Verbosity&                verbosity)
{
  extern const Numeric DEG2RAD;
  const Index ns = abs_species.nelem(); if(ns not_eq 1) throw std::runtime_error("Only for a single species in this test version");
  const Index nl = abs_lines_per_species[0].nelem();
  const Index np = p_grid.nelem();
  if(nz % 2) throw std::runtime_error("Cannot hit the tangent point in this test version; nz must be even");
  if(nf % 2 not_eq 1) throw std::runtime_error("Must include central frequency to test validity; nf must be uneven.");
  
  // Find Zenith angles and weigh them by their area
  Vector za_grid; nlinspace(za_grid, 0.0, 180.0, nz);
  Vector wzad(nz-1); 
  for(Index iz=0; iz<nz-1; iz++) 
    wzad[iz] = cos(DEG2RAD*za_grid[iz]) - cos(DEG2RAD*za_grid[iz+1]);
  
  iy = Matrix(nl, np, 0.0);
  iy_transmission = Tensor3(1, nl, np, 0.0);
  
  Tensor7 doit_i_field;
  Tensor3 trans_field;
  for(Index il=0; il<nl; il++) {
    const Numeric d = abs_lines_per_species[0][il].F() * df;
    Vector f_grid; nlinspace(f_grid, abs_lines_per_species[0][il].F()-d, abs_lines_per_species[0][il].F()+d, nf);
    doit_i_fieldClearskyPlaneParallel(ws, doit_i_field, trans_field, propmat_clearsky_agenda,
                                      water_p_eq_agenda, iy_space_agenda, iy_surface_agenda, 
                                      iy_cloudbox_agenda, 1, f_grid, atmosphere_dim,
                                      p_grid, z_field, t_field, nlte_field, vmr_field, abs_species,
                                      Tensor3(0, 0, 0), Tensor3(0, 0, 0), Tensor3(0, 0, 0),
                                      Tensor3(0, 0, 0), Tensor3(0, 0, 0), Tensor3(0, 0, 0),
                                      z_field(0, joker, joker), 1e99, 0.0, surface_props_data,
                                      za_grid, verbosity);
    
    // Integrate over the sphere
    for(Index ip=0; ip<np; ip++) {
      Eigen::VectorXcd F(nf);
      Linefunctions::set_lineshape(F, MapToEigen(f_grid), abs_lines_per_species[0][il], 
                                   Vector(1, vmr_field(0, ip, 0, 0)),  t_field(ip, 0, 0), p_grid[ip], 0.0,
                                   abs_species, 0, 0);
      Numeric sx = 0;
      for(Index iv=0; iv<nf-1; iv++) {
        const Numeric intF = (F[iv].real() + F[iv+1].real()) * (f_grid[iv+1] - f_grid[iv]);
        for(Index iz=0; iz<nz-1; iz++) {
          const Numeric x = wzad[iz] * 0.25 * intF;
          sx += x;
          iy(il, ip) += (doit_i_field(iv,   ip, 0, 0, iz,   0, 0) + 
                         doit_i_field(iv,   ip, 0, 0, iz+1, 0, 0) +
                         doit_i_field(iv+1, ip, 0, 0, iz,   0, 0) + 
                         doit_i_field(iv+1, ip, 0, 0, iz+1, 0, 0)) * 0.25 * x;
                         
          const Numeric exp_mtau = 0.25 * (trans_field(iv,   ip, iz  ) + 
                                           trans_field(iv,   ip, iz+1) + 
                                           trans_field(iv+1, ip, iz  ) + 
                                           trans_field(iv+1, ip, iz+1));
          if(abs(1-exp_mtau) < 1e-6)
          { /* do nothing */ }
          else 
            iy_transmission(0, il, ip) += (1 + (1 - exp_mtau)/log(exp_mtau)) * x;
        }
      }
      
      if(abs(sx-1) > 1e-3) {
        ostringstream os;
        os << "Integrated iy normalizes to " << sx << " instead of 1\n";
        os << "This means your frequency grid spanning " << f_grid[0] << " to " << f_grid[nf-1] << " is not good enough\n";
        if(sx < 1)
          os << "Please consider increasing df by about the inverse of the missing normalizations (a factor about " << 1/sx<<" is approximately enough)\n";
        else
          os << "Please consider increasing nf so that the frequency grid better represents the lineshape\n";
        throw std::runtime_error(os.str());
      }
    }
  }
}




/* Workspace method: Doxygen documentation will be auto-generated */
void doit_i_fieldClearskyPlaneParallel(
        Workspace&                  ws,
        Tensor7&                    doit_i_field,
        Tensor3&                    trans_field,
  const Agenda&                     propmat_clearsky_agenda,
  const Agenda&                     water_p_eq_agenda,   
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
          Tensor4        ppvar_trans_cumulat, ppvar_trans_partial;
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
                              ppvar_trans_cumulat,ppvar_trans_partial,iy_id, stokes_dim, f_grid,
                              atmosphere_dim,p_grid, z_field, t_field, nlte_field,
                              vmr_field, abs_species,wind_u_field, wind_v_field,
                              wind_w_field, mag_u_field, mag_v_field, mag_w_field,
                              cloudbox_on, iy_unit, iy_aux_vars, jacobian_do,
                              jacobian_quantities, ppath, rte_pos2,
                              propmat_clearsky_agenda, water_p_eq_agenda,
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
              trans_field(joker,0,i)    = ppvar_trans_partial(0,joker,0,0);
              trans_field(joker,nl-1,i) = ppvar_trans_partial(ppath.np-1,joker,0,0);
            }
          else
            {
              doit_i_field(joker,nl-1,0,0,i,0,joker) = ppvar_iy(joker,joker,0);
              doit_i_field(joker,0,0,0,i,0,joker)    = ppvar_iy(joker,joker,ppath.np-1);
              trans_field(joker,nl-1,i) = ppvar_trans_partial(0,joker,0,0);
              trans_field(joker,0,i)    = ppvar_trans_partial(ppath.np-1,joker,0,0);
            }

          // Remaining points
          for( Index p=1; p<ppath.np-1; p++ )
            {
              // We just store values at pressure levels
              if( ppath.gp_p[p].fd[0] < 1e-2 )
                {
                  doit_i_field(joker,ppath.gp_p[p].idx,0,0,i,0,joker) =
                    ppvar_iy(joker,joker,p);
                  trans_field(joker,ppath.gp_p[p].idx,i) = ppvar_trans_partial(p,joker,0,0);
                }
            }
        }
      catch (const std::runtime_error &e)
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
