/* Copyright (C) 2016
   Oliver Lemke <olemke@core-dump.info>

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
   USA.
*/

/*!
  \file   tessem.cc

  \brief  This file contains functions that are adapted from TESSEM2
  code which is used to calculate surface emissivity.

  The implementation is based on the Fortran code v1.0 developed by Catherine
  Prigent and Filipe Aires in the EUMETSAT Study on surface emissivity at
  microwave and sub-millimeter frequencies project EUM/CO/14/4600001473/CJA
*/

#include "tessem.h"
#include "file.h"
#include "matpackI.h"
#include "mystring.h"

/*! Read TESSEM2 neural network parameters

  \param[in,out] is Input file stream
  \param[out] net Neural network parameters
*/
void tessem_read_ascii(std::ifstream& is, TessemNN& net) {
  is >> net.nb_inputs >> net.nb_cache >> net.nb_outputs;

  net.b1.resize(net.nb_cache);
  for (Index i = 0; i < net.nb_cache; i++) is >> double_imanip() >> net.b1[i];

  net.b2.resize(net.nb_outputs);
  for (Index i = 0; i < net.nb_outputs; i++) is >> double_imanip() >> net.b2[i];

  net.w1.resize(net.nb_cache, net.nb_inputs);
  for (Index i = 0; i < net.nb_cache; i++)
    for (Index j = 0; j < net.nb_inputs; j++) is >> double_imanip() >> net.w1(i, j);

  net.w2.resize(net.nb_outputs, net.nb_cache);
  for (Index i = 0; i < net.nb_outputs; i++)
    for (Index j = 0; j < net.nb_cache; j++) is >> double_imanip() >> net.w2(i, j);

  net.x_min.resize(net.nb_inputs);
  for (Index i = 0; i < net.nb_inputs; i++) is >> double_imanip() >> net.x_min[i];

  net.x_max.resize(net.nb_inputs);
  for (Index i = 0; i < net.nb_inputs; i++) is >> double_imanip() >> net.x_max[i];

  net.y_min.resize(net.nb_outputs);
  for (Index i = 0; i < net.nb_outputs; i++) is >> double_imanip() >> net.y_min[i];

  net.y_max.resize(net.nb_outputs);
  for (Index i = 0; i < net.nb_outputs; i++) is >> double_imanip() >> net.y_max[i];
}

/*! Tessem emissivity calculation

  When using the default neural network parameter files
  from the Tessem 2 distribution, the input Vector should contain
  5 elements:
    - Frequency (10e9-700e9) in Hz
    - Theta (0-90) Incidence angle in degrees
    - Windspeed (0-25) at 10m in m/s.
      Higher wind speed can be used but without garantee
    - Surface skin temperature (270-310) in K
    - Salinity (0.0-0.04) in kg/kg

  \param[out] ny  Calculated emissivity.
  \param[in] net  Neural network parameters.
  \param[in] nx  Input data.
*/
void tessem_prop_nn(VectorView& ny, const TessemNN& net, ConstVectorView nx) {
  ARTS_USER_ERROR_IF (nx.nelem() != net.nb_inputs,
    "Tessem NN requires ", net.nb_inputs,
    " values, but input vector has ", nx.nelem(), " element.")

  ARTS_USER_ERROR_IF (ny.nelem() != net.nb_outputs,
    "Tessem NN generates ", net.nb_outputs,
    " values, but output vector has ", ny.nelem(), " element.")

  // preprocessing
  Vector new_x(nx);
  new_x[0] *= 1e-9;
  new_x[4] *= 1e3;
  for (Index i = 0; i < net.nb_inputs; i++)
    new_x[i] =
        -1. + (new_x[i] - net.x_min[i]) / (net.x_max[i] - net.x_min[i]) * 2;

  // propagation
  Vector trans = net.b1;
  for (Index i = 0; i < net.nb_cache; i++) {
    for (Index j = 0; j < net.nb_inputs; j++)
      trans[i] += net.w1(i, j) * new_x[j];

    trans[i] = 2. / (1. + exp(-2. * trans[i])) - 1.;
  }

  Vector new_y = net.b2;
  for (Index i = 0; i < net.nb_outputs; i++) {
    for (Index j = 0; j < net.nb_cache; j++)
      new_y[i] += net.w2(i, j) * trans[j];
  }

  // postprocessing
  for (Index i = 0; i < net.nb_outputs; i++)
    ny[i] = net.y_min[i] + (new_y[i] + 1.) / 2. * (net.y_max[i] - net.y_min[i]);
}
