/* Copyright (C) 2005-2012 Cory Davis <cdavis@staffmail.ed.ac.uk>
                            
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

/**
 * @file   mc_antenna.h
 * @author Cory Davis <cdavis@staffmail.ed.ac.uk>
 * @date   2005-12-02
 *
 * @brief  Workspace functions for the solution of cloud-box radiative transfer
 * by Monte Carlo methods.  All of these functions refer to 3D calculations.
 *
 *** FIXMEDOC *** : set_lookup; ran_uniform
 */
/*===========================================================================
  === External declarations
  ===========================================================================*/
#ifndef mc_antenna_h
#define mc_antenna_h

#include "arts.h"
#include "matpack_data.h"
#include "rng.h"
#include <cmath>
#include <stdexcept>

enum AntennaType {
  ANTENNA_TYPE_PENCIL_BEAM = 1,
  ANTENNA_TYPE_GAUSSIAN = 2,
  ANTENNA_TYPE_LOOKUP = 3
};

/** An Antenna object used by MCGeneral.
 * This class provides the means of sampling various types of 2D antenna
 * functions.
 */
struct MCAntenna {
  AntennaType atype{};
  Numeric sigma_aa{}, sigma_za{};
  Vector aa_grid, za_grid;
  Matrix G_lookup;

  MCAntenna()
      : aa_grid(),
        za_grid(),
        G_lookup() { /* Nothing to do here */
  }

  /** set_pencil_beam
   *
   * Makes the antenna pattern a pencil beam.
   */
  void set_pencil_beam();
  
  /** set_gaussian.
   *
   * Makes the antenna pattern a 2D gaussian specified by za and aa standard deviations.
   * Gives the MCAntenna object a 2D gaussian response function.
   *
   * @param[in]  za_sigma   The std. dev. parameter for zenith angle.
   * @param[in]  aa_sigma   The std. dev. parameter for azimuthal angle.
   *
   * @author     Cory Davis
   * @date       2005-12-02
   */
  void set_gaussian(const Numeric& za_sigma, const Numeric& aa_sigma);
    
  /** set_gaussian_fwhm.
   *
   * Makes the antenna pattern a 2D gaussian specified by za and aa FWHM.
   * Gives the MCAntenna object a 2D gaussian response function.
   *
   * @param[in]    za_fwhm    The full width half maximum zenith angle.
   * @param[in]    aa_fwhm    The full width half maximum azimuthal angle.
   *
   * @author       Cory Davis
   * @date         2005-12-02
   */
  void set_gaussian_fwhm(const Numeric& za_fwhm, const Numeric& aa_fwhm);

  /** set_lookup.
   *
   * Makes the antenna pattern use a 2D lookup table to define the antenna response.
   * The lookup antenna type is not yet implemented. *** FIXMEDOC ***
   *
   * @param[in]  za_grid_ zenith angle grid for the antenna response lookup table.
   * @param[in]  aa_grid_ azimuthal angle grid for the antenna response lookup table.
   * @param[in]  G_lookup_ the lookup table data.
   *
   * @author     Cory Davis
   * @date       2005-12-02
   */
  void set_lookup(ConstVectorView za_grid,
                  ConstVectorView aa_grid,
                  ConstMatrixView G_lookup);
    
  /** return_los
   *
   * Returns the normalized Gaussian weight for a photon line of sight
   * relative to the boresight.
   *
   * Modified 2016-09-07 by ISA to take a rotation matrix instead of
   * boresight los for reasons of computational efficiency.
   *
   * @param[out] wgt            Line-of-sight propagation vector in ENU frame.
   * @param[in]  rte_los        The line-of-sight of incoming photon.
   * @param[in]  bore_sight_los the bore sight LOS.
   *
   * @param[in]  R_enu2ant      Rotation matrix from ENU frame to antenna frame.
   *
   * @author     Ian S. Adams.
   * @date       2015-09-09.
   */
  void return_los(Numeric& wgt,
                  ConstMatrixView R_return,
                  ConstMatrixView R_enu2ant) const;
    
  /** draw_los.
   *
   * Draws a line of sight by sampling the antenna response function.
   *
   * @param[out]  sampled_rte_los  The sampled line of sight.
   * @param[out]  R_los            Line-of-sight propagation vector in ENU frame.
   * @param[in]   rng              A random number generator.
   * @param[in]   R_ant2enu        Rotation matrix from antenna frame to ENU frame.
   * @param[in]   bore_sight_los   The bore sight LOS.
   *
   * @author      Cory Davis
   * @date        2005-12-02
   */
  void draw_los(VectorView sampled_rte_los,
                MatrixView R_los,
                Rng& rng,
                ConstMatrixView R_ant2enu,
                ConstVectorView bore_sight_los) const;

  friend ostream& operator<<(ostream& os, const MCAntenna& mca);
};

/** ran_gaussian.
 *
 * Returns the gaussian random deviate.
 * Draw a random normal (Gaussian) deviate.
 * Polar (Box-Mueller) method; See Knuth v2, 3rd ed, p122.
 *
 * @param[in]   rng    Rng random number generator instance.
 * @param[in]   sigma  Standard deviation parameter for gaussian distribution.
 *
 * @author      Cory Davis
 * @date        2003-12-01
 */
Numeric ran_gaussian(Rng& rng, const Numeric sigma);

/** ran_uniform.  *** FIXMEDOC ***
 *
 */
Numeric ran_uniform(Rng& rng);

/** rotmat_enu.
 *
 * Calculates rotation matrix from antenna frame to ENU frame.
 * The columns of the rotation matrix are the v, h, and k
 * components of the propagating wave in the ENU frame.
 *
 * @param[out]   R_ant2enu  rotation matrix from antenna frame to ENU frame.
 * @param[in]    prop_los   los (zenith and azimuth).
 *
 * @author       Ian S. Adams
 * @date         2016-09-07
 */
void rotmat_enu(MatrixView R_ant2enu, ConstVectorView prop_los);

/** rotmat_stokes.
 *
 * Calculates the PRA matrix for the stokes vector
 * to account for polarization rotation from ENU
 * frame to antenna frame. Designed to handle sign
 * properly for radiometer and radar (both tx and rx)
 * using the bs_dir argument which (1 = away from sensor,
 * -1 = into sensor), based on Mishchenko's convention
 * for third Stokes. The assumption is that the
 * polarization basis vectors have magnitude of one;
 * therefore, a check is not made for the purpose of
 * computational efficiency.
 *
 * @param[out]    R_pra   rotation matrix.
 * @param[in]     stokes_dim  number of stokes vector elements to consider.
 * @param[in]     f1_dir  propagation direction of polarization basis 1 (-1.0 or 1.0).
 * @param[in]     f2_dir  propgation direction of polarization basis 2 (-1.0 or 1.0).
 * @param[in]     R_f1    rotation matrix (into ENU) for basis f1.
 * @param[in]     R_f2    photon rotation (into ENU) for basis f2.
 *
 * @author        Ian S. Adams
 * @date          2016-09-07
 */
void rotmat_stokes(MatrixView R_pra,
                   const Index& stokes_dim,
                   const Numeric& f1_dir,
                   const Numeric& f2_dir,
                   ConstMatrixView R_f1,
                   ConstMatrixView R_f2);

#endif  // mc_antenna_h
