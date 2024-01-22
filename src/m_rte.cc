/**
  @file   m_rte.cc
  @author Patrick Eriksson <patrick.eriksson@chalmers.se>
  @date   2002-05-11

  @brief  Workspace methods for solving clear sky radiative transfer.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <workspace.h>

#include <cmath>

#include "arts_constants.h"
#include "check_input.h"
#include "debug.h"
#include "rte.h"
#include "surf.h"

inline constexpr Numeric PI = Constant::pi;
inline constexpr Numeric SPEED_OF_LIGHT = Constant::speed_of_light;

/*===========================================================================
  === Workspace methods
  ===========================================================================*/

void spectral_radiance_backgroundFromMatrix(StokvecVector &background_rad,
                                            const Matrix &iy) {
  ARTS_USER_ERROR_IF(iy.ncols() not_eq 4, "Only for stokes dimensions 4.")
  background_rad = rtepack::to_stokvec_vector(iy);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ppvar_optical_depthFromPpvar_trans_cumulat(
    Matrix &ppvar_optical_depth, const Tensor4 &ppvar_trans_cumulat) {
  ppvar_optical_depth = ppvar_trans_cumulat(joker, joker, 0, 0);
  transform(ppvar_optical_depth, log, ppvar_optical_depth);
  ppvar_optical_depth *= -1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void yApplyUnit(Vector &y,
                Matrix &jacobian,
                const Vector &y_f,
                const ArrayOfIndex &y_pol,
                const String &iy_unit) {
  ARTS_USER_ERROR_IF(iy_unit == "1",
                     "No need to use this method with *iy_unit* = \"1\".");

  ARTS_USER_ERROR_IF(
      max(y) > 1e-3,
      "The spectrum vector *y* is required to have original radiance\n"
      "unit, but this seems not to be the case. This as a value above\n"
      "1e-3 is found in *y*.")

  // Is jacobian set?
  //
  const Index ny = y.size();
  //
  const bool do_j = jacobian.nrows() == ny;

  // Some jacobian quantities can not be handled
  ARTS_USER_ERROR_IF(
      do_j && max(jacobian) > 1e-3,
      "The method can not be used with jacobian quantities that are not\n"
      "obtained through radiative transfer calculations. One example on\n"
      "quantity that can not be handled is *jacobianAddPolyfit*.\n"
      "The maximum value of *jacobian* indicates that one or several\n"
      "such jacobian quantities are included.")

  // Planck-Tb
  //--------------------------------------------------------------------------
  if (iy_unit == "PlanckBT") {
    // Hard to use telescoping here as the data are sorted differently in y
    // and jacobian, than what is expected apply_iy_unit. Copy to temporary
    // variables instead.

    // Handle the elements in "frequency chunks"

    Index i0 = 0;  // Index of first element for present chunk
    //
    while (i0 < ny) {
      // Find number of values for this chunk
      Index n = 1;
      //
      while (i0 + n < ny && y_f[i0] == y_f[i0 + n]) {
        n++;
      }

      Matrix yv(1, n);
      ArrayOfIndex i_pol(n);
      bool any_quv = false;
      //
      for (Index i = 0; i < n; i++) {
        const Index ix = i0 + i;
        yv(0, i) = y[ix];
        i_pol[i] = y_pol[ix];
        if (i_pol[i] > 1 && i_pol[i] < 5) {
          any_quv = true;
        }
      }

      // Index of elements to convert
      Range ii(i0, n);

      if (do_j) {
        ARTS_USER_ERROR_IF(
            any_quv && i_pol[0] != 1,
            "The conversion to PlanckBT, of the Jacobian and "
            "errors for Q, U and V, requires that I (first Stokes "
            "element) is at hand and that the data are sorted in "
            "such way that I comes first for each frequency.")

        // Jacobian
        Tensor3 J(jacobian.ncols(), 1, n);
        J(joker, 0, joker) = transpose(jacobian(ii, joker));
        apply_iy_unit2(
            J, yv, iy_unit, ExhaustiveConstVectorView{y_f[i0]}, 1, i_pol);
        jacobian(ii, joker) = transpose(J(joker, 0, joker));
      }

      // y (must be done last)
      apply_iy_unit(yv, iy_unit, ExhaustiveConstVectorView{y_f[i0]}, 1, i_pol);
      y[ii] = yv(0, joker);

      i0 += n;
    }
  }

  // Other conversions
  //--------------------------------------------------------------------------
  else {
    // Here we take each element of y separately.

    Matrix yv(1, 1);
    ArrayOfIndex i_pol(1);

    for (Index i = 0; i < ny; i++) {
      yv(0, 0) = y[i];
      i_pol[0] = y_pol[i];

      // Jacobian
      if (do_j) {
        apply_iy_unit2(
            ExhaustiveTensor3View{ExhaustiveMatrixView{jacobian(i, joker)}},
            yv,
            iy_unit,
            ExhaustiveConstVectorView{y_f[i]},
            1,
            i_pol);
      }

      // y (must be done last)
      apply_iy_unit(yv, iy_unit, ExhaustiveConstVectorView{y_f[i]}, 1, i_pol);
      y[i] = yv(0, 0);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void y_geo_seriesFromY_geo(Matrix &y_geo_series,
                           const Matrix &y_geo,
                           const Vector &sensor_response_f_grid) {
  // Sizes
  const Index ly = y_geo.nrows();
  const Index nchannel = sensor_response_f_grid.size();
  const Index lseries = ly / nchannel;

  ARTS_USER_ERROR_IF(
      nchannel * lseries != ly,
      "Row size of *y_geo* not an even multiple of length of *sensor_response_f_grid*.")

  y_geo_series.resize(lseries, y_geo.ncols());

  Index i = 0;
  for (Index s = 0; s < lseries; ++s) {
    y_geo_series(s, joker) = y_geo(i, joker);
    i += nchannel;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void y_geo_swathFromY_geo(Tensor3 &y_geo_swath,
                          const Matrix &y_geo,
                          const Vector &sensor_response_f_grid,
                          const Index &npixel) {
  // Sizes
  const Index ly = y_geo.nrows();
  const Index nchannel = sensor_response_f_grid.size();
  const Index nswath = ly / (nchannel * npixel);

  ARTS_USER_ERROR_IF(
      nchannel * npixel * nswath != ly,
      "Row size of *y_geo* does not match given *npixel* and *sensor_response_f_grid*.")

  y_geo_swath.resize(nswath, npixel, y_geo.ncols());

  Index i = 0;
  for (Index s = 0; s < nswath; ++s) {
    for (Index p = 0; p < npixel; ++p) {
      y_geo_swath(s, p, joker) = y_geo(i, joker);
      i += nchannel;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void y_seriesFromY(Matrix &y_series,
                   const Vector &y,
                   const Vector &y_f,
                   const Vector &sensor_response_f_grid,
                   const Index &safe) {
  // Sizes
  const Index ly = y.size();
  const Index nchannel = sensor_response_f_grid.size();
  const Index lseries = ly / nchannel;

  ARTS_USER_ERROR_IF(
      nchannel * lseries != ly,
      "Length of *y* not an even multiple of length of *sensor_response_f_grid*.")

  y_series.resize(lseries, nchannel);

  Index i = 0;
  for (Index s = 0; s < lseries; ++s) {
    for (Index c = 0; c < nchannel; ++c) {
      if (safe && s > 0) {
        ARTS_USER_ERROR_IF(fabs(y_f[i] - y_f[i - nchannel]) > 1,
                           "At least one channel varies in frequency.")
      }
      y_series(s, c) = y[i++];
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void y_swathFromY(Tensor3 &y_swath,
                  const Vector &y,
                  const Vector &y_f,
                  const Vector &sensor_response_f_grid,
                  const Index &npixel,
                  const Index &safe) {
  // Sizes
  const Index ly = y.size();
  const Index nchannel = sensor_response_f_grid.size();
  const Index nswath = ly / (nchannel * npixel);

  ARTS_USER_ERROR_IF(
      nchannel * npixel * nswath != ly,
      "Length of *y* does not match given *npixel* and *sensor_response_f_grid*.")

  y_swath.resize(nswath, npixel, nchannel);

  Index i = 0;
  for (Index s = 0; s < nswath; ++s) {
    for (Index p = 0; p < npixel; ++p) {
      for (Index c = 0; c < nchannel; ++c) {
        if (safe && (p > 0 || s > 0)) {
          ARTS_USER_ERROR_IF(fabs(y_f[i] - y_f[i - nchannel]) > 1,
                             "At least one channel varies in frequency.")
        }
        y_swath(s, p, c) = y[i++];
      }
    }
  }
}
