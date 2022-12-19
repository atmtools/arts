/* Copyright (C) 2019
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/*!
 * \file   arts_options.h
 * \brief  Options for ARTS from enumeration (including error handling)
 * 
 * \author Richard Larsson
 * \date   2019-04-01
 */

#ifndef OPTIONS_IN_ARTS_H
#define OPTIONS_IN_ARTS_H

#include "enums.h"

namespace Options {
/** Keep time options available to switch over them */
ENUMCLASS(
    TimeStep, char, hour, hours, h, minute, minutes, min, second, seconds, s)

/** Possible line shape coefficients */
ENUMCLASS(LineShapeCoeff, char, X0, X1, X2, X3)

/** Possible Jacobian for Wind and Magnetic field */
ENUMCLASS(WindMagJacobian, char, u, v, w, strength)

/** Possible Jacobian for basic line parameters */
ENUMCLASS(BasicCatParamJacobian, char, LineStrength, LineCenter)

/** Possible HITRAN types of catalogs */
ENUMCLASS(
    HitranType,
    char,
    Pre2004,   // 2004 version changed the .par-length
    Post2004,  // New par length
    Online  // Onine expects a modern .par line followed by Upper then Lower quantum numbers
)

/** Possible AddLines Speedups */
ENUMCLASS(LblSpeedup, char, None, QuadraticIndependent, LinearIndependent)

ENUMCLASS(SortingOption, char, ByFrequency, ByEinstein)

/** Options for setting iy_main_agenda */
ENUMCLASS(iy_main_agendaDefaultOptions,
          char,
          Emission,
          EmissionPlaneParallel,
          Clearsky,
          Transmission,
          TransmissionUnitUnpolIntensity,
          TransmissionUnitPolIntensity,
          Freqloop,
          ScattMC)

/** Options for setting iy_loop_freqs_agenda */
ENUMCLASS(iy_loop_freqs_agendaDefaultOptions,
          char,
          Emission,
          Transmission)

/** Options for setting iy_space_agenda */
ENUMCLASS(iy_space_agendaDefaultOptions,
          char,
          CosmicBackground)

/** Options for setting iy_surface_agenda */
ENUMCLASS(iy_surface_agendaDefaultOptions,
          char,
          UseSurfaceRtprop)

/** Options for setting iy_cloudbox_agenda */
ENUMCLASS(iy_cloudbox_agendaDefaultOptions,
          char,
          LinInterpField,
          QuarticInterpField)

/** Options for setting ppath_agenda */
ENUMCLASS(ppath_agendaDefaultOptions,
          char,
          FollowSensorLosPath,
          PlaneParallel,
          TransmitterReceiverPath)

/** Options for setting ppath_step_agenda */
ENUMCLASS(ppath_step_agendaDefaultOptions,
          char,
          GeometricPath,
          RefractedPath)

/** Options for setting refr_index_air_agenda */
ENUMCLASS(refr_index_air_agendaDefaultOptions,
          char,
          NoRefrac,
          GasMicrowavesEarth,
          GasInfraredEarth,
          GasMicrowavesGeneral,
          FreeElectrons,
          GasMicrowavesGeneralAndElectrons,
          GasMicrowavesEarthAndElectrons)

/** Options for setting water_p_eq_agenda */
ENUMCLASS(water_p_eq_agendaDefaultOptions,
          char,
          MK05)

/** Options for setting gas_scattering_agenda */
ENUMCLASS(gas_scattering_agendaDefaultOptions,
          char,
          Dummy)

/** Options for setting surface_rtprop_agenda */
ENUMCLASS(surface_rtprop_agendaDefaultOptions,
          char,
          Blackbody_SurfTFromt_surface,
          Blackbody_SurfTFromt_field,
          Specular_NoPol_ReflFix_SurfTFromt_surface,
          Specular_NoPol_ReflFix_SurfTFromt_field,
          Specular_WithPol_ReflFix_SurfTFromt_surface,
          lambertian_ReflFix_SurfTFromt_surface,
          lambertian_ReflFix_SurfTFromt_field)

/** Options for setting g0_agenda */
ENUMCLASS(g0_agendaDefaultOptions,
          char,
          Earth,
          Io,
          Jupiter,
          Mars,
          Venus)

/** Options for setting dobatch_calc_agenda */
ENUMCLASS(dobatch_calc_agendaDefaultOptions, char)

/** Options for setting ybatch_calc_agenda */
ENUMCLASS(ybatch_calc_agendaDefaultOptions, char)

/** Options for setting test_agenda */
ENUMCLASS(test_agendaDefaultOptions, char)

/** Options for setting surface_rtprop_sub_agenda */
ENUMCLASS(surface_rtprop_sub_agendaDefaultOptions, char)

/** Options for setting spt_calc_agenda */
ENUMCLASS(spt_calc_agendaDefaultOptions, char)

/** Options for setting sensor_response_agenda */
ENUMCLASS(sensor_response_agendaDefaultOptions, char)

/** Options for setting propmat_clearsky_agenda */
ENUMCLASS(propmat_clearsky_agendaDefaultOptions, char, Empty)

/** Options for setting pha_mat_spt_agenda */
ENUMCLASS(pha_mat_spt_agendaDefaultOptions, char)

/** Options for setting met_profile_calc_agenda */
ENUMCLASS(met_profile_calc_agendaDefaultOptions, char)

/** Options for setting main_agenda */
ENUMCLASS(main_agendaDefaultOptions, char)

/** Options for setting jacobian_agenda */
ENUMCLASS(jacobian_agendaDefaultOptions, char)

/** Options for setting iy_radar_agenda */
ENUMCLASS(iy_radar_agendaDefaultOptions, char)

/** Options for setting iy_independent_beam_approx_agenda */
ENUMCLASS(iy_independent_beam_approx_agendaDefaultOptions, char)

/** Options for setting inversion_iterate_agenda */
ENUMCLASS(inversion_iterate_agendaDefaultOptions, char)

/** Options for setting forloop_agenda */
ENUMCLASS(forloop_agendaDefaultOptions, char)

/** Options for setting doit_scat_field_agenda */
ENUMCLASS(doit_scat_field_agendaDefaultOptions, char)

/** Options for setting doit_rte_agenda */
ENUMCLASS(doit_rte_agendaDefaultOptions, char)

/** Options for setting doit_mono_agenda */
ENUMCLASS(doit_mono_agendaDefaultOptions, char)

/** Options for setting doit_conv_test_agenda */
ENUMCLASS(doit_conv_test_agendaDefaultOptions, char)

/** Options for setting planets */
ENUMCLASS(planetDefaultOptions,
          char,
          Earth,
          Io,
          Jupiter,
          Mars,
          Venus)
/** Options for setting pnd_agenda --- CHANGE TO ENUMCLASS WHEN ADDING ANY OPTIONS AND REMOVE THIS PART OF THE COMMENT */
    ENUMCLASS_EMPTY(pnd_agendaDefaultOptions, char)
}  // namespace Options


#endif
