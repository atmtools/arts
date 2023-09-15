#include <sun.h>

/*===========================================================================
  === Functions in sun.h
  ===========================================================================*/

/** Calculates the radiance spectrum of sun which is scattered by the atmospheric gases.
 *
 * @param[in,out] ws ARTS workspace.
 * @param[out] scattered_sunlight StokvecVector scattered monochromatic radiance
 *             spectrum of sun.
 * @param[in] f_grid Vector frequency grid.
 * @param[in] p Numeric pressure at location of scattering.
 * @param[in] T Numeric temperature at location of scattering.
 * @param[in] vmr Vector volume mixing ratios of absorption species at location
 *            of scattering.
 * @param[in] transmitted_sunlight Matrix transmitted monochromatic irradiance
 *             spectrum of sun at location of scattering.
 * @param[in] gas_scattering_los_in Vector incoming direction of the transmitted sun irradiance
 *            spectrum.
 * @param[in] gas_scattering_los_out outgoing direction of the transmitted sun irradiance
 *            spectrum.
 * @param[in] gas_scattering_agenda Agenda agenda calculating the gas scattering
 *            cross sectionand matrix.
 */
void get_scattered_sunsource(const Workspace& ws,
                              StokvecVector& scattered_sunlight,
                              const Vector& f_grid,
                              const AtmPoint& atm_point,
                              const Matrix& transmitted_sunlight,
                              const Vector& gas_scattering_los_in,
                              const Vector& gas_scattering_los_out,
                              const Agenda& gas_scattering_agenda);
