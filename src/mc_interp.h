/**
 * @file   mc_interp.h
 * @author Cory Davis <cory@met.ed.ac.uk>
 * @date   2005-02-28
 *
 * @brief  Interpolation classes and functions created for use within Monte
 *         Carlo scattering simulations.
 *
 *** FIXMEDOC ***: interp_scat_angle_temperature (check)
 */
/*===========================================================================
  === External declarations
  ===========================================================================*/

#ifndef mc_interp_h
#define mc_interp_h
#include "array.h"
#include "arts.h"
#include "interpolation.h"
#include "matpack_data.h"
#include "optproperties.h"
#include "ppath.h"

/** A 2D sequential linear interpolation (SLI) lookup table
 * This class holds the gridded for 2D SLI as well as the
 * interpolate member function for retrieving interpolated values.
 *
 * @param[in] x1  desired x1 value
 * @param[in] x2  desired x2 value
 * @return    y   interpolated y value at x1,x2
 *
 * @author    Cory Davis <cdavis@staffmail.ed.ac.uk>
 * @date      2005-02-28
 */
class SLIData2 {
 public:
  //grid of x1 values where y is known
  Vector x1a;
  //A vector of x2 values for every x1a
  ArrayOfVector x2a;
  //y values for every x1a, x2a
  ArrayOfVector ya;
  //performs SLI.
  Numeric interpolate(Numeric x1, Numeric x2) const;
  //checks that it is not empty
  //void check() const;

  friend ostream& operator<<(ostream& os, const SLIData2& sli);
};

/** interp.
 *
 * Red 1D Interpolate.
 *
 * This is a slight modifiaction of Stefan's code to do 1_D interpolation
 * to get a Matrix from an array of Matrices.
 *
 * The dimension of itw must be consistent with the dimension of the
 * interpolation (2^n).
 *
 * @param[out] tia  Interpolated value.
 * @param[in]  itw  Interpolation weights.
 * @param[in]  a    The field to interpolate.(ArrayOfMatrix).
 * @param[in]  tc   The grid position for the column dimension.
 *
 * @author Cory Davis (modified original code by Stefan Buehler)
 * @date   2003-06-19
 */
void interp(MatrixView tia,
            ConstVectorView itw,
            const ArrayOfMatrix& a,
            const GridPos& tc);

/** interp.
 *
 * Red 1D Interpolate.
 *
 * This is a slight modifiaction of Stefan's code to do 1_D interpolation
 * to get a Vector from an array of Vectors.
 *
 * The dimension of itw must be consistent with the dimension of the
 * interpolation (2^n).
 *
 * @param[out] tia  Interpolated value.
 * @param[in]  itw  Interpolation weights.
 * @param[in]  a    The field to interpolate.(ArrayOfVector).
 * @param[in]  tc   The grid position for the column dimension.
 *
 * @author Cory Davis (modified original code by Stefan Buehler)
 * @date   2003-06-19
 */
void interp(VectorView tia,
            ConstVectorView itw,
            const ArrayOfVector& a,
            const GridPos& tc);

/** interp_scat_angle_temperature.
 *
 * Returns the interpolated phase matrix for given incident and scattered directions.
 *
 * @param[out] pha_mat_int       Interpolated phase matrix
 * @param[out] theta_rad         Scattering angle defined by inc and sca direct.
 * @param[in]  scat_data_single  A monochromatic SingleScatteringData object.
 * @param[in]  za_sca            Zenith angle of scattering direction.
 * @param[in]  aa_sca            Azimuth angle of scattering direction.
 * @param[in]  za_inc            Zenith angle of incident direction.
 * @param[in]  aa_inc            Azimuth angle of incident direction.
 * @param[in]  rtp_temperature   As the WSV.
 *
 * @author     Cory Davis
 * @date       2005-02-28
 */
void interp_scat_angle_temperature(  //Output:
    VectorView pha_mat_int,
    Numeric& theta_rad,
    //Input:
    const SingleScatteringData& scat_data_single,
    const Numeric& za_sca,
    const Numeric& aa_sca,
    const Numeric& za_inc,
    const Numeric& aa_inc,
    const Numeric& rtp_temperature);

#endif  // mc_interp_h
