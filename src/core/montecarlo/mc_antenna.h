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
#pragma once

#include <enumsAntennaType.h>
#include <matpack.h>
#include <rng.h>
#include <rtepack.h>

/** An Antenna object used by MCGeneral.
 * This class provides the means of sampling various types of 2D antenna
 * functions.
 */
struct MCAntenna {
  AntennaType atype{};
  Numeric sigma_aa{}, sigma_za{};
  Vector aa_grid{}, za_grid{};
  Matrix G_lookup{};

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
  [[nodiscard]] Numeric return_los(const Matrix33& R_return,
                                   const Matrix33& R_enu2ant) const;

  /** draw_los.
   *
   * Draws a line of sight by sampling the antenna response function.
   *
   * @param[in]   rng              A random number generator.
   * @param[in]   R_ant2enu        Rotation matrix from antenna frame to ENU frame.
   * @param[in]   bore_sight_los   The bore sight LOS.
   * @return      sampled_rte_los  The sampled line of sight.
   * @return      R_los            Line-of-sight propagation vector in ENU frame.
   *
   * @author      Cory Davis
   * @date        2005-12-02
   */
  [[nodiscard]] std::pair<Vector2, Matrix33> draw_los(
      RandomNumberGenerator<>& rng,
      const Matrix33& R_ant2enu,
      const Vector2& bore_sight_los) const;
};

template <>
struct xml_io_stream_name<MCAntenna> {
  static constexpr std::string_view name = "MCAntenna";
};

template <>
struct xml_io_stream_aggregate<MCAntenna> {
  static constexpr bool value = true;
};

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
[[nodiscard]] Matrix33 rotmat_enu(const Vector2& prop_los);

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
[[nodiscard]] Muelmat rotmat_stokes(Numeric f1_dir,
                                    Numeric f2_dir,
                                    const Matrix33& R_f1,
                                    const Matrix33& R_f2);
