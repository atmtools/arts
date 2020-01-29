/* Copyright (C) 2020
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
 * @file   m_fullmodel.cc
 * @author Richard Larsson
 * @date   2020-01-29
 * 
 * @brief  Full absorption models of various kinds
 */


#include "absorption.h"


void abs_xsec_per_speciesAddO2Lines(ArrayOfMatrix& abs_xsec_per_species,
                                    ArrayOfArrayOfMatrix& dabs_xsec_per_species_dx,
                                    const ArrayOfArrayOfSpeciesTag& abs_species,
                                    const ArrayOfRetrievalQuantity& jacobian_quantities,
                                    const Vector& f_grid,
                                    const Vector& abs_p,
                                    const Vector& abs_t,
                                    const Matrix& abs_vmrs,
                                    const Verbosity& verbosity)
{
}
