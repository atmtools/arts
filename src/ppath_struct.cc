/* Copyright (C) 2002-2012 Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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

/**
 * @file   ppath_struct.cc
 * @author Manfred Brath <manfred.brath@uni-hamburg.de>
 * @date   2022-09-028
 *
 * @brief  Propagation path structure functions.
 *
 * This file contains the definition of the Ppath structure function.
 */

#include "ppath_struct.h"

std::ostream& operator<<(std::ostream& os, const Ppath& x) {
  os << "dim: " << x.dim << "\n";
  os << "np: " << x.np << "\n";
  os << "constant: " << x.constant << "\n";
  os << "background: " << x.background << "\n";
  os << "start_pos: " << x.start_pos << "\n";
  os << "start_los: " << x.start_los << "\n";
  os << "start_lstep: " << x.start_lstep << "\n";
  os << "pos: " << x.pos << "\n";
  os << "los: " << x.los << "\n";
  os << "r: " << x.r << "\n";
  os << "lstep: " << x.lstep << "\n";
  os << "end_pos: " << x.end_pos << "\n";
  os << "end_los: " << x.end_los << "\n";
  os << "end_lstep: " << x.end_lstep << "\n";
  os << "nreal: " << x.nreal << "\n";
  os << "ngroup: " << x.ngroup << "\n";
  os << "gp_p: " << x.gp_p << "\n";
  if (x.dim >= 2) {
    os << "gp_lat: " << x.gp_lat << "\n";
  }
  if (x.dim == 3) {
    os << "gp_lon: " << x.gp_lon << "\n";
  }
  return os;
}
