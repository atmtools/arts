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

/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   mc_antenna.h
  \author Cory Davis <cdavis@staffmail.ed.ac.uk>
  \date   2005-12-02 

  \brief  Workspace functions for the solution of cloud-box radiative transfer 
by Monte Carlo methods.  All of these functions refer to 3D calculations

  Modified 2015-09-09 Ian S. Adams
  Added class members to handle radar returns.
  Fixed gaussian line-of-sight calculation.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/
/*===========================================================================
  === External declarations
  ===========================================================================*/
#ifndef mc_antenna_h
#define mc_antenna_h

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "matpackI.h"
#include "ppath.h"
#include "rng.h"

enum AntennaType {
  ANTENNA_TYPE_PENCIL_BEAM = 1,
  ANTENNA_TYPE_GAUSSIAN = 2,
  ANTENNA_TYPE_LOOKUP = 3
};

//! An Antenna object used by MCGeneral
/*! This class provides the means of sampling various types of 2D antenna
functions.. */
class MCAntenna {
  AntennaType atype;
  Numeric sigma_aa, sigma_za;
  Vector aa_grid, za_grid;
  Matrix G_lookup;

 public:
  MCAntenna()
      : atype(),
        sigma_aa(0.),
        sigma_za(0.),
        aa_grid(),
        za_grid(),
        G_lookup() { /* Nothing to do here */
  }

  void set_pencil_beam(void);
  void set_gaussian(const Numeric& za_sigma, const Numeric& aa_sigma);
  void set_gaussian_fwhm(const Numeric& za_fwhm, const Numeric& aa_fwhm);
  void set_lookup(ConstVectorView za_grid,
                  ConstVectorView aa_grid,
                  ConstMatrixView G_lookup);
  AntennaType get_type(void) const;
  void return_los(Numeric& wgt,
                  ConstMatrixView R_return,
                  ConstMatrixView R_enu2ant) const;
  void draw_los(VectorView sampled_rte_los,
                MatrixView R_los,
                Rng& rng,
                ConstMatrixView R_ant2enu,
                ConstVectorView bore_sight_los) const;
};

ostream& operator<<(ostream& os, const MCAntenna& mca);

Numeric ran_gaussian(Rng& rng, const Numeric sigma);

Numeric ran_uniform(Rng& rng);

void rotmat_enu(MatrixView R_ant2enu, ConstVectorView prop_los);

void rotmat_stokes(MatrixView R_pra,
                   const Index& stokes_dim,
                   const Numeric& bs_dir,
                   const Numeric& prop_dir,
                   ConstMatrixView R_bs,
                   ConstMatrixView R_prop);

#endif  // mc_antenna_h
