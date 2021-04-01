/* Copyright (C) 2021
   Jon Petersen <jon.petersen@studium.uni-hamburg.de>
   Manfred Brath  <manfred.brath@uni-hamburg.de>

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
  \file   star.cc
  \author Jon Petersen <jon.petersen@studium.uni-hamburg.de>
          Manfred Brath  <manfred.brath@uni-hamburg.de>
  \date   2021-02-22

  \brief  Functions needed for radiative transfer with direct sources.
*/

#include "star.h"
#include "auto_md.h"
#include "agenda_class.h"
#include "propagationmatrix.h"

void get_scattered_starsource(Workspace& ws,
                              StokesVector& scattered_starlight,
                              const Vector& f_grid,
                              const Numeric& p,
                              const Numeric& T,
                              const Vector& vmr,
                              const Matrix& transmitted_starlight,
                              const Vector& in_los,
                              const Vector& out_los,
                              const Agenda& gas_scattering_agenda) {
  PropagationMatrix K_sca;
  TransmissionMatrix sca_mat;


  // calculate gas scattering properties
  gas_scattering_agendaExecute(ws,
                               K_sca,
                               sca_mat,
                               f_grid,
                               p,
                               T,
                               vmr,
                               in_los,
                               out_los,
                               gas_scattering_agenda);

  //some basic quantities
  Index ns = transmitted_starlight.ncols();
  Index nf = f_grid.nelem();

  //allocate and resize
  Vector temp(ns);
//  Matrix Z(ns, ns);
//  Vector temp(ns);

  // Calculate the scattered radiation
  RadiationVector scattered_starlight_temp(nf,ns);
  scattered_starlight_temp=transmitted_starlight;
  scattered_starlight_temp.leftMul(sca_mat);


  // Richard will change the type of S to RadiationVector in iyClearsky
  //but for now we have to convert it
  for (Index i_f = 0; i_f < nf; i_f++) {

      temp=scattered_starlight_temp.Vec(i_f);
      temp*=K_sca.Kjj(0,0)[i_f];
      scattered_starlight.SetAtPosition(temp,i_f);
  }


}