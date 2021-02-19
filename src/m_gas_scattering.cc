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
   USA. */

/*===========================================================================
  ===  File description
  ===========================================================================*/

#include "agenda_class.h"
#include "arts.h"
#include "rte.h"


/*!
  \file   m_gas_scattering.cc
  \author Jon Petersen  <jon.petersen@studium.uni-hamburg.de>,
          Manfred Brath  <manfred.brath@.uni-hamburg.de>
  \date   2021-02-08

  \brief  Workspace functions related to gas scattering.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric BOLTZMAN_CONST;

/*===========================================================================
  === The functions
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void gas_scatteringOff(Index& gas_scattering_do,
                       Agenda& gas_scattering_agenda,
                       const Verbosity&) {
  // set flag to False (default)
  gas_scattering_do = 0;

  gas_scattering_agenda = Agenda();
  gas_scattering_agenda.set_name("gas_scattering_agenda");
}

void gas_scatteringXsecConst(ArrayOfPropagationMatrix& sca_xsec,
                             const Vector& f_grid,
                             const Index& stokes_dim,
                             const Vector& ppvar_p,
                             const Vector& ppvar_t,
                             const Numeric& ConstXsec,
                             const Verbosity&) {
  // Some basic sizes
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppvar_p.nelem();

  ArrayOfPropagationMatrix dummy(np, PropagationMatrix(nf, ns));


  Vector Xsec(nf,ConstXsec);

  // Number density
  Numeric N;


  for (Index ip = 0; ip < np; ++ip) {

    N=ppvar_p[ip] / ppvar_t[ip] /  BOLTZMAN_CONST;

    dummy[ip].Kjj() += Xsec;
    dummy[ip].Kjj() *= N;
  }

  sca_xsec = dummy;
}