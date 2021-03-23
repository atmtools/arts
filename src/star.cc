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
                              const StokesVector& transmitted_starlight,
                              const Vector& in_los,
                              const Vector& out_los,
                              const Agenda& gas_scattering_agenda) {
  PropagationMatrix K_sca;
  PropagationMatrix sca_mat;

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


  Index ns = transmitted_starlight.StokesDimensions();
  Index nf = f_grid.nelem();

  Vector I_in(ns, 0);
  Matrix Z(ns, ns, 0);

  Vector temp(ns, 0);

  for (Index i_f = 0; i_f < nf; i_f++) {
    I_in = transmitted_starlight.VectorAtPosition(i_f);

    sca_mat.MatrixAtPosition(Z, i_f);

    temp *= 0;

    for (Index i_s = 0; i_s < ns; i_s++) {
      for (Index j_s = 0; j_s < ns; j_s++) {
        temp[i_s] += Z(i_s, j_s) * I_in[j_s];
      }
    }
    scattered_starlight.SetAtPosition(temp, i_f);
  }
}