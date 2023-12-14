!> SUBROUTINE for ARTS surfem_ocean
SUBROUTINE surfem_ocean(f_GHz, theta, windspeed, SST, SSS, phi, transmittance, emissivity, reflectivity)

USE parkind1, ONLY : jprb
USE rttov_surfem_ocean_mod, ONLY : rttov_surfem_ocean
IMPLICIT NONE
REAL(jprb), INTENT(IN)  :: f_GHz, theta, windspeed, SST, SSS, phi, transmittance
REAL(jprb), INTENT(OUT) :: emissivity(4), reflectivity(4)

CALL rttov_surfem_ocean(f_GHz, theta, windspeed, SST, SSS, phi, transmittance, emissivity, reflectivity)
RETURN

END SUBROUTINE surfem_ocean




