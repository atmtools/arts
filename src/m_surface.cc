/*===========================================================================
  === File description
  ===========================================================================*/

/*!
  \file   m_surface.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2008-09-17

  \brief  Workspace functions associated wih the surface and its properties.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <iomanip>
#include "array.h"
#include "arts_constants.h"
#include <workspace.h>
#include "atm.h"
#include "check_input.h"
#include "debug.h"
#include "matpack_complex.h"
#include "fastem.h"
#include "geodetic.h"
#include "geodetic.h"
#include "interpolation.h"
#include "math_funcs.h"
#include "matpack_data.h"
#include "matpack_math.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"
#include "surf.h"
#include "surface.h"
#include "tessem.h"
#include "arts_conversions.h"
#include "gas_scattering.h"

inline constexpr Numeric EARTH_RADIUS=Constant::earth_radius;
inline constexpr Numeric DEG2RAD = Conversion::deg2rad(1);


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void FastemStandAlone(Matrix& emissivity,
                      Matrix& reflectivity,
                      const Vector& f_grid,
                      const Numeric& surface_skin_t,
                      const Numeric& za,
                      const Numeric& salinity,
                      const Numeric& wind_speed,
                      const Numeric& rel_aa,
                      const Vector& transmittance,
                      const Index& fastem_version) {
  const Index nf = f_grid.size();

  chk_if_in_range("zenith angle", za, 90, 180);
  chk_if_in_range_exclude("surface skin temperature", surface_skin_t, 260, 373);
  chk_if_in_range_exclude_high("salinity", salinity, 0, 1);
  chk_if_in_range_exclude_high("wind speed", wind_speed, 0, 100);
  chk_if_in_range("azimuth angle", rel_aa, -180, 180);
  chk_vector_length("transmittance", "f_grid", transmittance, f_grid);
  if (fastem_version < 3 || fastem_version > 6)
    throw std::runtime_error(
        "Invalid fastem version: 3 <= fastem_version <= 6");

  emissivity.resize(nf, 4);
  reflectivity.resize(nf, 4);

  const Numeric t = std::max(surface_skin_t, Numeric(270));

  for (Index i = 0; i < nf; i++) {
    if (f_grid[i] > 250e9)
      throw std::runtime_error("Only frequency <= 250 GHz are allowed");
    chk_if_in_range("transmittance", transmittance[i], 0, 1);

    Vector e, r;

    fastem(e,
           r,
           f_grid[i],
           za,
           t,
           salinity,
           wind_speed,
           transmittance[i],
           rel_aa,
           fastem_version);

    emissivity(i, joker) = e;
    reflectivity(i, joker) = r;
  }

  // FASTEM does not work close to the horizon (at least v6). Make sure values
  // are inside [0,1]. Then seems best to make sure that e+r=1.
  for (Index i = 0; i < nf; i++) {
    for (Index s = 0; s < 2; s++) {
      if (emissivity(i, s) > 1) {
        emissivity(i, s) = 1;
        reflectivity(i, s) = 0;
      }
      if (emissivity(i, s) < 0) {
        emissivity(i, s) = 0;
        reflectivity(i, s) = 1;
      }
      if (reflectivity(i, s) > 1) {
        emissivity(i, s) = 0;
        reflectivity(i, s) = 1;
      }
      if (reflectivity(i, s) < 0) {
        emissivity(i, s) = 1;
        reflectivity(i, s) = 0;
      }
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void InterpSurfaceFieldToPosition(SurfacePoint& surface_point,
                                  const Vector& rtp_pos,
                                  const SurfaceField& surface_field,
                                  const Numeric& surface_search_accuracy) {
  ARTS_USER_ERROR_IF(not surface_field.has(Surf::Key::h),
      "No elevation in the surface field, required by method")

  ARTS_USER_ERROR_IF(surface_search_accuracy < 0.0,
                     "Cannot have negative surface search accuracy")

  const Numeric alt = rtp_pos[0];
  const Numeric lat = rtp_pos[1];
  const Numeric lon = rtp_pos[2];
  surface_point = surface_field.at(lat, lon);

  const bool cmp = std::abs(alt - surface_point.elevation) > surface_search_accuracy;
  ARTS_USER_ERROR_IF(cmp,
                     "rtp_pos is not close enough to the surface.\nThe surface accuracy is: ",
                     surface_search_accuracy,
                     " m.\nThe surface altitude is: ",
                     surface_point.elevation,
                     " m.\nThe rtp_pos altitude is: ",
                     alt, " m.\n")
}

void surface_pointFromAtm(SurfacePoint& surface_point,
                          const AtmField& atm_field,
                          const Vector& rtp_pos) {
  surface_point = SurfacePoint{};
  surface_point.elevation = rtp_pos[0];
  surface_point.normal = {0, 0};
  
  ARTS_USER_ERROR_IF(not atm_field.has(Atm::Key::t),
                     "\"atm_field\" has no temperature field")
  surface_point.temperature =
      atm_field[Atm::Key::t].at(rtp_pos[0], rtp_pos[1], rtp_pos[2]);
  
  if (atm_field.has(Atm::Key::wind_u) and atm_field.has(Atm::Key::wind_v) and
      atm_field.has(Atm::Key::wind_w)) {
    surface_point.wind = {
        atm_field[Atm::Key::wind_u].at(rtp_pos[0], rtp_pos[1], rtp_pos[2]),
        atm_field[Atm::Key::wind_v].at(rtp_pos[0], rtp_pos[1], rtp_pos[2]),
        atm_field[Atm::Key::wind_w].at(rtp_pos[0], rtp_pos[1], rtp_pos[2])};
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iySurfaceInit(Matrix& iy,
                   const Vector& f_grid) {
  iy.resize(f_grid.size(), 4);
  iy = 0.;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iySurfaceRtpropAgenda(const Workspace& ws,
                           Matrix& iy,
                           ArrayOfTensor3& diy_dx,
                           SurfacePoint& surface_point,
                           Matrix& surface_los,
                           Tensor4& surface_rmatrix,
                           Matrix& surface_emission,
                           const Tensor3& iy_transmittance,
                           const Index& iy_id,
                           const Index& jacobian_do,
                           const Index& suns_do,
                           const AtmField& atm_field,
                           const Index& cloudbox_on,
                           const Vector& f_grid,
                           const Vector& rtp_pos,
                           const Vector& rtp_los,
                           const Vector& rte_pos2,
                           const String& iy_unit,
                           const Agenda& iy_main_agenda,
                           const Agenda& surface_rtprop_agenda) {
  // Input checks
  chk_rte_pos(rtp_pos);
  chk_rte_los(rtp_los);

  // Call *surface_rtprop_agenda*
  surface_rtprop_agendaExecute(ws,
                               surface_point,
                               surface_emission,
                               surface_los,
                               surface_rmatrix,
                               f_grid,
                               rtp_pos,
                               rtp_los,
                               surface_rtprop_agenda);

  // Check output of *surface_rtprop_agenda*
  const Index nlos = surface_los.nrows();
  const Index nf = f_grid.size();
  //
  if (nlos)  // if 0, blackbody ground and not all checks are needed
  {
    if (surface_los.ncols() != rtp_los.size())
      throw std::runtime_error("Number of columns in *surface_los* is not correct.");
    if (nlos != surface_rmatrix.nbooks())
      throw std::runtime_error(
          "Mismatch in size of *surface_los* and *surface_rmatrix*.");
    if (surface_rmatrix.npages() != nf)
      throw std::runtime_error(
          "Mismatch in size of *surface_rmatrix* and *f_grid*.");
    if (surface_rmatrix.nrows() != 4 ||
        surface_rmatrix.ncols() != 4)
      throw std::runtime_error(
          "Wrong size *surface_rmatrix*.");
  }
  if (surface_emission.ncols() != 4)
    throw std::runtime_error(
        "Wrong size *surface_emission*.");
  if (surface_emission.nrows() != nf)
    throw std::runtime_error("Mismatch in size of *surface_emission* and f_grid*.");

  // Variable to hold down-welling radiation
  Tensor3 I(nlos, nf, 4);

  ArrayOfString iy_aux_var(0);
  if (suns_do) iy_aux_var.emplace_back("Direct radiation");

  // Loop *surface_los*-es. If no such LOS, we are ready.
  if (nlos > 0) {
    for (Index ilos = 0; ilos < nlos; ilos++) {
      Vector los{surface_los(ilos, joker)};

      // Include surface reflection matrix in *iy_transmittance*
      // If iy_transmittance is empty, this is interpreted as the
      // variable is not needed.
      //
      Tensor3 iy_trans_new;
      //
      if (iy_transmittance.npages()) {
        iy_transmittance_mult(iy_trans_new,
                             iy_transmittance,
                             surface_rmatrix(ilos, joker, joker, joker));
      }

      // Calculate downwelling radiation for LOS ilos
      //
      {
        ArrayOfMatrix iy_aux;
        Ppath ppath;
        Vector geo_pos;
        Index iy_id_new = iy_id + ilos + 1;
        iy_main_agendaExecute(ws,
                              iy,
                              iy_aux,
                              ppath,
                              diy_dx,
                              geo_pos,
                              0,
                              iy_trans_new,
                              iy_aux_var,
                              iy_id_new,
                              iy_unit,
                              cloudbox_on,
                              jacobian_do,
                              f_grid,
                              atm_field,
                              rtp_pos,
                              los,
                              rte_pos2,
                              iy_main_agenda);

        //For the case that a sun is present and the los is towards a sun, we
        //subtract the direct radiation, so that only the diffuse radiation is considered here.
        //If sun is within los iy_aux[0] is the incoming and attenuated direct (sun)
        //radiation at the surface. Otherwise it is zero.
        if (suns_do){
          iy-=iy_aux[0];
        }

      }

      if (iy.ncols() != 4 || iy.nrows() != nf) {
        std::ostringstream os;
        os << "The size of *iy* returned from *" << iy_main_agenda.get_name()
           << "* is\n"
           << "not correct:\n"
           << "  expected size = [" << nf << "," << 4 << "]\n"
           << "  size of iy    = [" << iy.nrows() << "," << iy.ncols() << "]\n";
        throw std::runtime_error(os.str());
      }

      I(ilos, joker, joker) = iy;
    }
  }

  // Add up
  surface_calc(iy, I, surface_los, surface_rmatrix, surface_emission);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iySurfaceRtpropCalc(const Workspace& ws,
                         Matrix& iy,
                         ArrayOfTensor3& diy_dx,
                         const Matrix& surface_los,
                         const Tensor4& surface_rmatrix,
                         const Matrix& surface_emission,
                         const ArrayOfString& dsurface_names,
                         const ArrayOfTensor4& dsurface_rmatrix_dx,
                         const ArrayOfMatrix& dsurface_emission_dx,
                         const Tensor3& iy_transmittance,
                         const Index& iy_id,
                         const Index& jacobian_do,
                         const Index& suns_do,
                         const ArrayOfRetrievalQuantity& jacobian_quantities,
                         const AtmField& atm_field,                         
                         const Index& cloudbox_on,
                         
                         const Vector& f_grid,
                         const Vector& rtp_pos,
                         const Vector& rtp_los,
                         const Vector& rte_pos2,
                         const String& iy_unit,
                         const Agenda& iy_main_agenda) {
  // Input checks
  chk_rte_pos(rtp_pos);
  chk_rte_los(rtp_los);

  // Check provided surface rtprop variables
  const Index nlos = surface_los.nrows();
  const Index nf = f_grid.size();
  //
  if (nlos)  // if 0, blackbody ground and not all checks are needed
  {
    if (surface_los.ncols() != rtp_los.size())
      throw std::runtime_error("Number of columns in *surface_los* is not correct.");
    if (nlos != surface_rmatrix.nbooks())
      throw std::runtime_error(
          "Mismatch in size of *surface_los* and *surface_rmatrix*.");
    if (surface_rmatrix.npages() != nf)
      throw std::runtime_error(
          "Mismatch in size of *surface_rmatrix* and *f_grid*.");
    if (surface_rmatrix.nrows() != 4 ||
        surface_rmatrix.ncols() != 4)
      throw std::runtime_error(
          "Wrong size *surface_rmatrix*.");
  }
  if (surface_emission.ncols() != 4)
    throw std::runtime_error(
        "Wrong size *surface_emission*.");
  if (surface_emission.nrows() != nf)
    throw std::runtime_error("Mismatch in size of *surface_emission* and f_grid*.");

  // Variable to hold down-welling radiation
  Tensor3 I(nlos, nf, 4);

  ArrayOfString iy_aux_var(0);
  if (suns_do) iy_aux_var.emplace_back("Direct radiation");

  // Loop *surface_los*-es.
  if (nlos > 0) {
    for (Index ilos = 0; ilos < nlos; ilos++) {
      Vector los{surface_los(ilos, joker)};

      // Include surface reflection matrix in *iy_transmittance*
      // If iy_transmittance is empty, this is interpreted as the
      // variable is not needed.
      //
      Tensor3 iy_trans_new;
      //
      if (iy_transmittance.npages()) {
        iy_transmittance_mult(iy_trans_new,
                             iy_transmittance,
                             surface_rmatrix(ilos, joker, joker, joker));
      }

      // Calculate downwelling radiation for LOS ilos
      //
      {
        ArrayOfMatrix iy_aux;
        Ppath ppath;
        Vector geo_pos;
        iy_main_agendaExecute(ws,
                              iy,
                              iy_aux,
                              ppath,
                              diy_dx,
                              geo_pos,
                              0,
                              iy_trans_new,
                              iy_aux_var,
                              iy_id,
                              iy_unit,
                              cloudbox_on,
                              jacobian_do,
                              f_grid,
                              atm_field,
                              rtp_pos,
                              los,
                              rte_pos2,
                              iy_main_agenda);

        //For the case that a sun is present and the los is towards a sun, we
        //subtract the direct radiation, so that only the diffuse radiation is considered here.
        //If sun is within los iy_aux[0] is the incoming and attenuated direct (sun)
        //radiation at the surface. Otherwise it is zero.
        if (suns_do){
          iy-=iy_aux[0];
        }

      }

      if (iy.ncols() != 4 || iy.nrows() != nf) {
        std::ostringstream os;
        os << "The size of *iy* returned from *" << iy_main_agenda.get_name()
           << "* is\n"
           << "not correct:\n"
           << "  expected size = [" << nf << "," << 4 << "]\n"
           << "  size of iy    = [" << iy.nrows() << "," << iy.ncols() << "]\n";
        throw std::runtime_error(os.str());
      }

      I(ilos, joker, joker) = iy;
    }
  }

  // Add up
  surface_calc(iy, I, surface_los, surface_rmatrix, surface_emission);

  // Surface Jacobians
  if (jacobian_do && dsurface_names.size()) {
    // Loop dsurface_names
    for (Size i = 0; i < dsurface_names.size(); i++) {
      // Error if derivatives not calculated
      // Or should we accept this?
      if (dsurface_emission_dx[i].empty() || dsurface_rmatrix_dx[i].empty()) {
        std::ostringstream os;
        os << "The derivatives for surface quantity: " << dsurface_names[i]
           << "\nwere not calculated by *iy_surface_agenda*.\n"
           << "That is, *dsurface_emission_dx* and/or *dsurface_rmatrix_dx*\n"
           << "are empty.";
        throw std::runtime_error(os.str());
      } else {
        // Find index among jacobian quantities
        Index ihit = -1;
        for (Size j = 0; j < jacobian_quantities.size() && ihit < 0; j++) {
          if (dsurface_names[i] == jacobian_quantities[j].Subtag()) {
            ihit = j;
          }
        }
        ARTS_ASSERT(ihit >= 0);
        // Derivative, as observed at the surface
        Matrix diy_dpos0, diy_dpos;
        surface_calc(diy_dpos0,
                     I,
                     surface_los,
                     dsurface_rmatrix_dx[i],
                     dsurface_emission_dx[i]);
        // Weight with transmission to sensor
        iy_transmittance_mult(diy_dpos, iy_transmittance, diy_dpos0);
        // Put into diy_dx
        diy_from_pos_to_rgrids(diy_dx[ihit],
                               jacobian_quantities[ihit],
                               diy_dpos,
                               rtp_pos);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void specular_losCalcOldNoTopography(Vector &specular_los,
                                     Vector &surface_normal,
                                     const Vector &rtp_los) {
  chk_rte_los(rtp_los);

  specular_los = {180 - rtp_los[0], rtp_los[1]};
  surface_normal = {0, 0};
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceBlackbody(Matrix& surface_los,
                      Tensor4& surface_rmatrix,
                      Matrix& surface_emission,
                      const Vector& f_grid,
                      
                      const Vector& rtp_pos,
                      const Vector& rtp_los,
                      const SurfacePoint& surface_point) {
  chk_rte_pos(rtp_pos);
  chk_rte_los(rtp_los);
  chk_not_negative("surface_skin_t", surface_point.temperature);

  surface_los.resize(0, 0);
  surface_rmatrix.resize(0, 0, 0, 0);

  const Index nf = f_grid.size();

  Vector b(nf);
  planck(b, f_grid, surface_point.temperature);

  surface_emission.resize(nf, 4);
  surface_emission = 0.0;

  for (Index iv = 0; iv < nf; iv++) {
    surface_emission(iv, 0) = b[iv];
    for (Index is = 1; is < 4; is++) {
      surface_emission(iv, is) = 0;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceTessem(Matrix& surface_los,
                   Tensor4& surface_rmatrix,
                   Matrix& surface_emission,
                   
                   const Vector& f_grid,
                   const Vector& rtp_pos,
                   const Vector& rtp_los,
                   const Numeric& surface_skin_t,
                   const TessemNN& net_h,
                   const TessemNN& net_v,
                   const Numeric& salinity,
                   const Numeric& wind_speed) {
  // Input checks
  chk_rte_pos(rtp_pos);
  chk_rte_los(rtp_los);
  chk_if_in_range_exclude(
      "surface skin temperature", surface_skin_t, 260.0, 373.0);
  chk_if_in_range_exclude_high("salinity", salinity, 0, 1);
  chk_if_in_range_exclude_high("wind speed", wind_speed, 0, 100);

  // Determine specular direction
  Vector specular_los, surface_normal;
  specular_losCalcOldNoTopography(specular_los,
                               surface_normal,
                               rtp_los);

  // TESSEM in and out
  //
  Vector out(2);
  VectorView e_h = out[Range(0, 1)];
  VectorView e_v = out[Range(1, 1)];
  //
  Vector in(5);
  in[1] = 180.0 - abs(rtp_los[0]);
  in[2] = wind_speed;
  in[3] = surface_skin_t;
  in[4] = salinity;

  // Get Rv and Rh
  //
  const Index nf = f_grid.size();
  Matrix surface_rv_rh(nf, 2);
  //
  for (Index i = 0; i < nf; ++i) {
    if (f_grid[i] < 5e9)
      throw std::runtime_error("Only frequency >= 5 GHz are allowed");
    if (f_grid[i] > 900e9)
      throw std::runtime_error("Only frequency <= 900 GHz are allowed");

    in[0] = f_grid[i];

    tessem_prop_nn(e_h, net_h, in);
    tessem_prop_nn(e_v, net_v, in);

    surface_rv_rh(i, 0) = std::min(std::max(1 - e_v[0], (Numeric)0), (Numeric)1);
    surface_rv_rh(i, 1) = std::min(std::max(1 - e_h[0], (Numeric)0), (Numeric)1);
  }

  surfaceFlatRvRh(surface_los,
                  surface_rmatrix,
                  surface_emission,
                  f_grid,
                  rtp_pos,
                  rtp_los,
                  specular_los,
                  surface_skin_t,
                  surface_rv_rh);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatRefractiveIndex(Matrix& surface_los,
                                Tensor4& surface_rmatrix,
                                Matrix& surface_emission,
                                const Vector& f_grid,
                                     const Vector& rtp_pos,
                                const Vector& rtp_los,
                                const Vector& specular_los,
                                const Numeric& surface_skin_t,
                                const GriddedField3& surface_complex_refr_index) {
  chk_rte_pos(rtp_pos);
  chk_rte_los(rtp_los);
  chk_rte_los(specular_los);
  chk_not_negative("surface_skin_t", surface_skin_t);

  // Interpolate *surface_complex_refr_index*
  //
  const Index nf = f_grid.size();
  //
  Matrix n_real(nf, 1), n_imag(nf, 1);
  //
  complex_n_interp(n_real,
                   n_imag,
                   surface_complex_refr_index,
                   "surface_complex_refr_index",
                   f_grid,
                   Vector(1, surface_skin_t));

  surface_los.resize(1, specular_los.size());
  surface_los(0, joker) = specular_los;

  surface_emission.resize(nf, 4);
  surface_rmatrix.resize(1, nf, 4, 4);

  // Incidence angle
  const Numeric incang = calc_incang(rtp_los, specular_los);
  ARTS_ASSERT(incang <= 90);

  // Complex (amplitude) reflection coefficients
  Complex Rv, Rh;

  for (Index iv = 0; iv < nf; iv++) {
    // Set n2 (refractive index of surface medium)
    Complex n2(n_real(iv, 0), n_imag(iv, 0));

    // Amplitude reflection coefficients
    fresnel(Rv, Rh, Numeric(1.0), n2, incang);

    // Fill reflection matrix and emission vector
    surface_specular_R_and_b(surface_rmatrix(0, iv, joker, joker),
                             surface_emission(iv, joker),
                             Rv,
                             Rh,
                             f_grid[iv],
                             surface_skin_t);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatReflectivity(Matrix& surface_los,
                             Tensor4& surface_rmatrix,
                             Matrix& surface_emission,
                             const Vector& f_grid,
                               const Vector& rtp_pos,
                             const Vector& rtp_los,
                             const Vector& specular_los,
                             const Numeric& surface_skin_t,
                             const Tensor3& surface_reflectivity) {
  chk_rte_pos(rtp_pos);
  chk_rte_los(rtp_los);
  chk_rte_los(specular_los);
  chk_not_negative("surface_skin_t", surface_skin_t);

  const Index nf = f_grid.size();

  ARTS_USER_ERROR_IF(
      surface_reflectivity.nrows() != 4 && surface_reflectivity.ncols() != 4,
      "The number of rows and columns in *surface_reflectivity* must\n"
      "be 4.\nThe number of rows in *surface_reflectivity* : ",
      surface_reflectivity.nrows(), "\n")

  ARTS_USER_ERROR_IF(
      surface_reflectivity.npages() != nf && surface_reflectivity.npages() != 1,
      "The number of pages in *surface_reflectivity* should\n",
      "match length of *f_grid* or be 1.", "\n length of *f_grid* : ", nf,
      "\n dimension of *surface_reflectivity* : ",
      surface_reflectivity.npages(), "\n")

  surface_los.resize(1, specular_los.size());
  surface_los(0, joker) = specular_los;

  surface_emission.resize(nf, 4);
  surface_rmatrix.resize(1, nf, 4, 4);

  Matrix R, IR(4, 4);

  Vector b(nf);
  planck(b, f_grid, surface_skin_t);

  Vector B(4, 0);

  for (Index iv = 0; iv < nf; iv++) {
    if (iv == 0 || surface_reflectivity.npages() > 1) {
      R = surface_reflectivity(iv, joker, joker);
      for (Index i = 0; i < 4; i++) {
        for (Index j = 0; j < 4; j++) {
          if (i == j) {
            IR(i, j) = 1 - R(i, j);
          } else {
            IR(i, j) = -R(i, j);
          }
        }
      }
    }

    surface_rmatrix(0, iv, joker, joker) = R;

    B[0] = b[iv];
    mult(surface_emission(iv, joker), IR, B);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatRvRh(Matrix& surface_los,
                     Tensor4& surface_rmatrix,
                     Matrix& surface_emission,
                     const Vector& f_grid,
                     
                     const Vector& rtp_pos,
                     const Vector& rtp_los,
                     const Vector& specular_los,
                     const Numeric& surface_skin_t,
                     const Matrix& surface_rv_rh) {
  chk_rte_pos(rtp_pos);
  chk_rte_los(rtp_los);
  chk_rte_los(specular_los);
  chk_not_negative("surface_skin_t", surface_skin_t);

  const Index nf = f_grid.size();

  if (surface_rv_rh.ncols() != 2) {
    std::ostringstream os;
    os << "The number of columns in *surface_rv_rh* must be two,\n"
       << "but the actual number of columns is " << surface_rv_rh.ncols()
       << "\n";
    throw std::runtime_error(os.str());
  }

  if (surface_rv_rh.nrows() != nf && surface_rv_rh.nrows() != 1) {
    std::ostringstream os;
    os << "The number of rows in *surface_rv_rh* should\n"
       << "match length of *f_grid* or be 1."
       << "\n length of *f_grid* : " << nf
       << "\n rows in *surface_rv_rh* : " << surface_rv_rh.nrows() << "\n";
    throw std::runtime_error(os.str());
  }

  if (min(surface_rv_rh) < 0 || max(surface_rv_rh) > 1) {
    throw std::runtime_error("All values in *surface_rv_rh* must be inside [0,1].");
  }

  surface_los.resize(1, specular_los.size());
  surface_los(0, joker) = specular_los;

  surface_emission.resize(nf, 4);
  surface_rmatrix.resize(1, nf, 4, 4);

  surface_emission = 0;
  surface_rmatrix = 0;

  Vector b(nf);
  planck(b, f_grid, surface_skin_t);

  Numeric rmean = 0.0, rdiff = 0.0;

  for (Index iv = 0; iv < nf; iv++) {
    if (iv == 0 || surface_rv_rh.nrows() > 1) {
      rmean = 0.5 * (surface_rv_rh(iv, 0) + surface_rv_rh(iv, 1));
      rdiff = 0.5 * (surface_rv_rh(iv, 0) - surface_rv_rh(iv, 1));
    }

    surface_emission(iv, 0) = (1.0 - rmean) * b[iv];
    surface_rmatrix(0, iv, 0, 0) = rmean;

    surface_emission(iv, 1) = -rdiff * b[iv];

    surface_rmatrix(0, iv, 0, 1) = rdiff;
    surface_rmatrix(0, iv, 1, 0) = rdiff;
    surface_rmatrix(0, iv, 1, 1) = rmean;

    for (Index i = 2; i < 4; i++) {
      surface_rmatrix(0, iv, i, i) = rmean;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatScalarReflectivity(Matrix& surface_los,
                                   Tensor4& surface_rmatrix,
                                   Matrix& surface_emission,
                                   const Vector& f_grid,
                                           const Vector& rtp_pos,
                                   const Vector& rtp_los,
                                   const Vector& specular_los,
                                   const Numeric& surface_skin_t,
                                   const Vector& surface_scalar_reflectivity) {
  chk_rte_pos(rtp_pos);
  chk_rte_los(rtp_los);
  chk_rte_los(specular_los);
  chk_not_negative("surface_skin_t", surface_skin_t);

  const Index nf = f_grid.size();

  if (surface_scalar_reflectivity.size() != nf &&
      surface_scalar_reflectivity.size() != 1) {
    std::ostringstream os;
    os << "The number of elements in *surface_scalar_reflectivity* should\n"
       << "match length of *f_grid* or be 1."
       << "\n length of *f_grid* : " << nf
       << "\n length of *surface_scalar_reflectivity* : "
       << surface_scalar_reflectivity.size() << "\n";
    throw std::runtime_error(os.str());
  }

  if (min(surface_scalar_reflectivity) < 0 ||
      max(surface_scalar_reflectivity) > 1) {
    throw std::runtime_error(
        "All values in *surface_scalar_reflectivity* must be inside [0,1].");
  }

  surface_los.resize(1, specular_los.size());
  surface_los(0, joker) = specular_los;

  surface_emission.resize(nf, 4);
  surface_rmatrix.resize(1, nf, 4, 4);

  surface_emission = 0;
  surface_rmatrix = 0;

  Vector b(nf);
  planck(b, f_grid, surface_skin_t);

  Numeric r = 0.0;

  for (Index iv = 0; iv < nf; iv++) {
    if (iv == 0 || surface_scalar_reflectivity.size() > 1) {
      r = surface_scalar_reflectivity[iv];
    }

    surface_emission(iv, 0) = (1.0 - r) * b[iv];
    surface_rmatrix(0, iv, 0, 0) = r;
    for (Index i = 1; i < 4; i++) {
      surface_rmatrix(0, iv, i, i) = r;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceLambertianSimple(Matrix& surface_los,
                             Tensor4& surface_rmatrix,
                             Matrix& surface_emission,
                             const Vector& f_grid,
                               const Vector& rtp_pos,
                             const Vector& rtp_los,
                             const Vector& surface_normal,
                             const Numeric& surface_skin_t,
                             const Vector& surface_scalar_reflectivity,
                             const Index& lambertian_nza,
                             const Numeric& za_pos) {
  const Index nf = f_grid.size();

  chk_rte_pos(rtp_pos);
  chk_rte_los(rtp_los);
  chk_not_negative("surface_skin_t", surface_skin_t);
  chk_if_in_range("za_pos", za_pos, 0, 1);

  if (surface_scalar_reflectivity.size() != nf &&
      surface_scalar_reflectivity.size() != 1) {
    std::ostringstream os;
    os << "The number of elements in *surface_scalar_reflectivity* should\n"
       << "match length of *f_grid* or be 1."
       << "\n length of *f_grid* : " << nf
       << "\n length of *surface_scalar_reflectivity* : "
       << surface_scalar_reflectivity.size() << "\n";
    throw std::runtime_error(os.str());
  }

  if (min(surface_scalar_reflectivity) < 0 ||
      max(surface_scalar_reflectivity) > 1) {
    throw std::runtime_error(
        "All values in *surface_scalar_reflectivity* must be inside [0,1].");
  }

  // Allocate and init everything to zero
  //
  surface_los.resize(lambertian_nza, rtp_los.size());
  surface_rmatrix.resize(lambertian_nza, nf, 4, 4);
  surface_emission.resize(nf, 4);
  //
  surface_los = 0.0;
  surface_rmatrix = 0.0;
  surface_emission = 0.0;

  // Help variables
  //
  const Numeric dza = (90.0 - abs(surface_normal[0])) / (Numeric)lambertian_nza;
  const Vector za_lims=uniform_grid(0.0, lambertian_nza + 1, dza);

  // surface_los
  for (Index ip = 0; ip < lambertian_nza; ip++) {
    surface_los(ip, 0) = za_lims[ip] + za_pos * dza;
    surface_los(ip, 1) = rtp_los[1];
  }

  Vector b(nf);
  planck(b, f_grid, surface_skin_t);

  // Loop frequencies and set remaining values
  //
  Numeric r = 0.0;
  //
  for (Index iv = 0; iv < nf; iv++) {
    // Get reflectivity
    if (iv == 0 || surface_scalar_reflectivity.size() > 1) {
      r = surface_scalar_reflectivity[iv];
    }

    // surface_rmatrix:
    // Only element (0,0) is set to be non-zero. This follows VDISORT
    // that refers to: K. L. Coulson, Polarization and Intensity of Light in
    // the Atmosphere (1989), page 229
    // (Thanks to Michael Kahnert for providing this information!)
    // Update: Is the above for a later edition? We have not found a copy of
    // that edition. In a 1988 version of the book, the relevant page seems
    // to be 232.
    for (Index ip = 0; ip < lambertian_nza; ip++) {
      const Numeric w =
          r * 0.5 *
          (cos(2 * DEG2RAD * za_lims[ip]) - cos(2 * DEG2RAD * za_lims[ip + 1]));
      surface_rmatrix(ip, iv, 0, 0) = w;
    }

    // surface_emission
    surface_emission(iv, 0) = (1 - r) * b[iv];
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_scalar_reflectivityFromSurface_rmatrix(
    Vector& surface_scalar_reflectivity,
    const Tensor4& surface_rmatrix) {
  const Index nf = surface_rmatrix.npages();
  const Index nlos = surface_rmatrix.nbooks();

  surface_scalar_reflectivity.resize(nf);
  surface_scalar_reflectivity = 0;

  for (Index i = 0; i < nf; i++) {
    for (Index l = 0; l < nlos; l++) {
      surface_scalar_reflectivity[i] += surface_rmatrix(l, i, 0, 0);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_rtpropFromTypesManual(const Workspace& ws,
                                   Vector& surface_type_mix,
                                   Numeric& surface_skin_t,
                                   Matrix& surface_los,
                                   Tensor4& surface_rmatrix,
                                   Matrix& surface_emission,
                                   const Vector& f_grid,
                                   const Vector& rtp_pos,
                                   const Vector& rtp_los,
                                   const ArrayOfAgenda& surface_rtprop_agenda_array,
                                   const Index& surface_type)
{
  ARTS_USER_ERROR_IF(surface_type < 0 or
     static_cast<Size>(surface_type) >= surface_rtprop_agenda_array.size(),
     "Provided surface type index invalid (<0 or too high w.r.t. "
     "length of *surface_rtprop_agenda_array*).");
  
  surface_type_mix.resize(surface_rtprop_agenda_array.size());
  surface_type_mix = 0.0;
  surface_type_mix[surface_type] = 1.0;

  surface_rtprop_agenda_arrayExecute(ws,
                                     surface_skin_t,
                                     surface_emission,
                                     surface_los,
                                     surface_rmatrix,
                                     surface_type,
                                     f_grid,
                                     rtp_pos,
                                     rtp_los,
                                     surface_rtprop_agenda_array);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_rtpropInterpFreq(Vector& f_grid,
                              Tensor4& surface_rmatrix,
                              Matrix& surface_emission,
                              const Vector& f_new) {
  const Index nf = f_grid.size();
  const Index ns = surface_emission.ncols();
  const Index nlos = surface_rmatrix.nbooks();
  const Index nnew = f_new.size();

  // Checks
  ARTS_USER_ERROR_IF(surface_emission.nrows() not_eq nf,
     "Different number of frequencies in *f_grid* and *surface_emission*.");
  ARTS_USER_ERROR_IF(surface_rmatrix.npages() not_eq nf,
     "Different number of frequencies in *f_grid* and *surface_rmatrix*.");
  ARTS_USER_ERROR_IF(surface_rmatrix.ncols() not_eq ns,
     "Different number of Stokes elements in *surface_emission* and *surface_rmatrix*.");

  // Set up interpolation
  chk_interpolation_grids("Frequency interpolation", f_grid, f_new);
  ArrayOfGridPos gp(nnew);
  Matrix itw(nnew, 2);
  gridpos(gp, f_grid, f_new);
  interpweights(itw, gp);

  // Interpolate
  Tensor4 rmatrix = surface_rmatrix;
  Matrix emission = surface_emission;
  surface_rmatrix.resize(nlos, nnew, ns, ns);
  surface_emission.resize(nnew, ns);
  //
  for (Index is = 0; is < ns; ++is) {
    interp(surface_emission(joker, is), itw, emission(joker, is), gp);
    for (Index il = 0; il < nlos; ++il) {
      for (Index is2 = 0; is2 < ns; ++is2) {
        interp(surface_rmatrix(il, joker, is, is2),
               itw,
               rmatrix(il, joker, is, is2),
               gp);
      }
    }
  }
  f_grid = f_new;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void SurfaceDummy(ArrayOfTensor4& dsurface_rmatrix_dx,
                  ArrayOfMatrix& dsurface_emission_dx,
                  const ArrayOfString& surface_props_names,
                  const ArrayOfString& dsurface_names,
                  const Index& jacobian_do) {
  ARTS_USER_ERROR_IF (surface_props_names.size(),
        "When calling this method, *surface_props_names* should be empty.")

  if (jacobian_do) {
    dsurface_check(surface_props_names,
                   dsurface_names,
                   dsurface_rmatrix_dx,
                   dsurface_emission_dx);
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void specular_losCalc(Vector& specular_los,
                      const SurfaceField& surface_field,
                      const Vector& rtp_pos,
                      const Vector& rtp_los,
                      const Index& ignore_topography)
{
  chk_rte_pos("rtp_pos", rtp_pos);
  chk_rte_los("rtp_los", rtp_los);
  
  specular_los.resize(2);
  specular_los_calc(specular_los,
                    surface_field,
                    rtp_pos[Range(1, 2)],
                    rtp_los,
                    ignore_topography);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void surface_normalCalc(Vector& surface_normal,
                        const SurfaceField& surface_field,
                        const Vector& rtp_pos,
                        const Index& ignore_topography)
{
  chk_rte_pos("rtp_pos", rtp_pos);

  surface_normal.resize(2);
  
  // No surface tilt if told so or surface_elevation.data has size (1,1)
  if (ignore_topography || surface_field.constant_value(Surf::Key::h)) {
    surface_normal = 0;

  } else {
    Vector pos(3), ecef(3), decef(3);
    surface_normal_calc(pos,
                        ecef,
                        decef,
                        surface_field,
                        rtp_pos[Range(1, 2)]);
  
    ecef2geodetic_los(pos, surface_normal, ecef, decef, surface_field.ellipsoid);
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void transmittanceFromIy_aux(Vector& transmittance,
                             const ArrayOfString& iy_aux_vars,
                             const ArrayOfMatrix& iy_aux)
{
  Index ihit = -1;

  for (Size i = 0; i < iy_aux_vars.size(); i++) {
    if (iy_aux_vars[i] == "Optical depth") {
      ihit = i;
      break;
    }
  }

  if (ihit < 0)
    throw std::runtime_error("No element in *iy_aux* holds optical depths.");

  const Index n = iy_aux[ihit].nrows();

  transmittance.resize(n);

  for (Index i = 0; i < n; i++) {
    transmittance[i] = exp(-iy_aux[ihit](i, 0));
  }
}
