/*===========================================================================
  === File description 
  ===========================================================================*/

/**
   @file   surface.cc
   @author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   @date   2012-02-06 

   This file contains internal functions associated with the surface.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "surface.h"

#include <workspace.h>

#include <cmath>
#include <iomanip>

#include "debug.h"
#include "geodetic.h"
#include "physics_funcs.h"
#include "surf.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

// Expression double-checked 210330 (PE)
Numeric calc_incang(ConstVectorView rte_los, ConstVectorView specular_los) {
  return (180 - abs(rte_los[0]) + abs(specular_los[0])) / 2;
}

Index index_of_zsurface(const Numeric& z_surface,
                        ConstVectorView z_profile) {
  Index ip = 0;
  while (z_surface >= z_profile[ip+1]) {
    ip++;
  }
  return ip;
}

void surface_calc(Matrix& iy,
                  ConstTensor3View I,
                  ConstMatrixView surface_los,
                  ConstTensor4View surface_rmatrix,
                  ConstMatrixView surface_emission) {
  // Some sizes
  const Index nf = I.nrows();
  const Index nlos = surface_los.nrows();

  iy = surface_emission;

  // Loop *surface_los*-es. If no such LOS, we are ready.
  if (nlos > 0) {
    for (Index ilos = 0; ilos < nlos; ilos++) {
      Vector rtmp(4);  // Reflected Stokes vector for 1 frequency

      for (Index iv = 0; iv < nf; iv++) {
        mult(rtmp, surface_rmatrix(ilos, iv, joker, joker), I(ilos, iv, joker));
        iy(iv, joker) += rtmp;
      }
    }
  }
}

void surface_specular_R_and_b(MatrixView surface_rmatrix,
                              VectorView surface_emission,
                              const Complex& Rv,
                              const Complex& Rh,
                              const Numeric& f,
                              const Numeric& surface_skin_t) {
  ARTS_ASSERT(surface_rmatrix.nrows() == 4);
  ARTS_ASSERT(surface_rmatrix.ncols() == 4);
  ARTS_ASSERT(surface_emission.size() == 4);

  // Expressions are derived in the surface chapter in the user guide

  surface_rmatrix = 0.0;
  surface_emission = 0.0;

  Numeric B = planck(f, surface_skin_t);

  const Numeric rv = pow(abs(Rv), 2.0);
  const Numeric rh = pow(abs(Rh), 2.0);
  const Numeric rmean = (rv + rh) / 2;

  surface_rmatrix(0, 0) = rmean;
  surface_emission[0] = B * (1 - rmean);

    const Numeric rdiff = (rv - rh) / 2;

    surface_rmatrix(1, 0) = rdiff;
    surface_rmatrix(0, 1) = rdiff;
    surface_rmatrix(1, 1) = rmean;
    surface_emission[1] = -B * rdiff;

      const Complex a = Rh * conj(Rv);
      const Complex b = Rv * conj(Rh);
      const Numeric c = real(a + b) / 2.0;

      surface_rmatrix(2, 2) = c;

        const Numeric d = imag(a - b) / 2.0;

        surface_rmatrix(2, 3) = d;
        surface_rmatrix(3, 2) = -d;
        surface_rmatrix(3, 3) = c;
}

void surface_props_check(const SurfaceField& surface_field,
                         const ArrayOfString& surface_props_names) {
  // Check sizes
  ARTS_USER_ERROR_IF (surface_field.size<SurfacePropertyTag>() != surface_props_names.size(),
        "The number of pages in *surface_props_data* and "
        "length of *surface_props_names* differ.");
  // If no surface properties, then we are ready
  if (surface_props_names.size() == 0) {
    return;
  }

  for (auto &name : surface_props_names)
    ARTS_USER_ERROR_IF(not surface_field.contains(SurfacePropertyTag{name}),
                       "No ", std::quoted(name), " field in surface_field")

  for (Size i = 0; i < surface_props_names.size(); i++) {
    ARTS_USER_ERROR_IF (surface_props_names[i].size() == 0,
      "Element ", i, " (0-based) of *surface_props_names* is empty.")
    for (Size j = i + 1; j < surface_props_names.size(); j++) {
      ARTS_USER_ERROR_IF (surface_props_names[j] == surface_props_names[i],
        "Two surface properties with same name found!\n"
        "This found for these two properties\n"
        "   index: ", i, '\n',
        "   index: ", j, '\n',
        "    name: ", surface_props_names[i])
    }
  }
}

void dsurface_check(const ArrayOfString& surface_props_names,
                    const ArrayOfString& dsurface_names,
                    const ArrayOfTensor4 dsurface_rmatrix_dx,
                    const ArrayOfMatrix& dsurface_emission_dx) {
  const Size nq = dsurface_names.size();

  ARTS_USER_ERROR_IF (dsurface_rmatrix_dx.size() != nq,
        "The lengths of *dsurface_names* and *dsurface_rmatrix_dx* differ.");
  ARTS_USER_ERROR_IF (dsurface_emission_dx.size() != nq,
        "The lengths of *dsurface_names* and *dsurface_emission_dx* differ.");

  for (Size i = 0; i < nq; i++) {
    bool found = false;
    for (Size j = 0; j < surface_props_names.size() && !found; j++) {
      if (dsurface_names[i] == surface_props_names[j]) {
        found = true;
      }
    }
    ARTS_USER_ERROR_IF (!found,
        "String ", i, " (0-based) of *dsurface_names* is \"",
        dsurface_names[i], "\"\n"
        "but this string could not be found in *surface_props_names*.\n"
        "This is likely due to incorrect choice of quantity when\n"
        " calling *jacobianAddSurfaceQuantity*.")
  }
}


void surface_normal_calc(VectorView pos,
                         VectorView ecef,
                         VectorView decef,
                         const SurfaceField& surface_field,
                         ConstVectorView pos2D)
{
  ARTS_ASSERT(pos.size() == 3); 
  ARTS_ASSERT(ecef.size() == 3); 
  ARTS_ASSERT(decef.size() == 3); 
  ARTS_ASSERT(pos2D.size() == 2);
  
  // We need two orthogonal vectors inside the surface plane. We
  // follow ENU and the first one should be towards E and the second
  // towards N. We obtain the vectors quite easily by dl shifts,
  // except when we are are very close to the North or South pole. To
  // be sure that a dl shift does not pass any of the poles, we
  // consider here that all positions inside a distance 5*dl on the
  // side. These points are shifted to be at the pole.
  //
  const Numeric dl = 1.0;  
  const Numeric lat_limit = 90.0 - RAD2DEG * 5 * dl / surface_field.ellipsoid[1];

  // Determine 3D pos at pos
  Numeric lat = pos2D[0];
  if (lat > lat_limit)
    lat = 90.0;
  else if (lat < -lat_limit)
    lat = -90.0;
  //
  pos[1] = lat;
  pos[2] = pos2D[1];
  pos[0] = surface_field.single_value(Surf::Key::h, pos[0], pos[2]);

  // Radius at pos0
  const Numeric r = pos[0] + prime_vertical_radius(surface_field.ellipsoid, lat);

  // Shifted positions
  Vector posWE{pos};
  Vector posSN{pos};
  //
  // North pole case
  if (lat > lat_limit) {
    posSN[1] -= RAD2DEG * dl / r;  
    posSN[2] = 90; 
    posWE[1] = posSN[1];
    posWE[2] = 0; 
  // South pole case
  } else if (lat < -lat_limit) {
    posSN[1] += RAD2DEG * dl / r;  
    posSN[2] = 0; 
    posWE[1] = posSN[1];
    posWE[2] = 90; 
  // A general case 
  } else {
    posSN[1] += RAD2DEG * dl / r;  
    posWE[2] += RAD2DEG * dl / (r * cos(DEG2RAD * posWE[1]));
    if (posWE[2] >= 180)
      posWE[2] -= 360;
  }
  //
  posSN[0] = surface_field.single_value(Surf::Key::h, posSN[1], posSN[2]);
  posWE[0] = surface_field.single_value(Surf::Key::h, posWE[1], posWE[2]);
  
  // Convert all three positions to ECEF
  Vector ecefSN(3), ecefWE(3);
  geodetic2ecef(ecef, pos, surface_field.ellipsoid);
  geodetic2ecef(ecefSN, posSN, surface_field.ellipsoid);
  geodetic2ecef(ecefWE, posWE, surface_field.ellipsoid);

  // Directional vectors to shifted positions
  Vector decefSN(3), decefWE(3);
  ecef_vector_distance(decefSN, ecef, ecefSN);
  ecef_vector_distance(decefWE, ecef, ecefWE);

  // Normal is cross product of the two decef (
  cross3(decef, decefWE, decefSN);
  decef /= norm2(decef);
}
