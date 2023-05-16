/*===========================================================================
  ===  File description
  ===========================================================================*/

#include "matpack_data.h"

/*!
  \file   gas_scattering.h
  \author Jon Petersen  <jon.petersen@studium.uni-hamburg.de>,
          Manfred Brath  <manfred.brath@.uni-hamburg.de>
  \date   2021-03-15

  \brief  Header file for functions related to gas scattering.
*/


#ifndef gas_scattering_h
#define gas_scattering_h

Vector calc_rayleighPhaMat(const Numeric& theta_rad,
                           const Index& stokes_dim);

#endif /* gas_scattering_h */