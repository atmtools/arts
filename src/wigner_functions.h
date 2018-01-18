/* Copyright (C) 2012
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

#include "rational.h"
#include <cmath>

Numeric wigner3j(const Rational j1,const Rational j2,const Rational j3,
                 const Rational m1,const Rational m2,const Rational m3);

Numeric wigner6j(const Rational j1,const Rational j2,const Rational j3,
                 const Rational l1,const Rational l2,const Rational l3);

Numeric ECS_wigner(Rational L, Rational Nl, Rational Nk, 
                   Rational Jk_lower, Rational Jl_lower, 
                   Rational Jk_upper, Rational Jl_upper);

Numeric factorials(const ArrayOfIndex& NomFac, const ArrayOfIndex& DenomFac);

void primes(ArrayOfIndex& output, const Index input);
void powers(ArrayOfIndex& output, const ArrayOfIndex primes, const Index input);

Numeric triangle_coefficient(const Rational a, const Rational b, const Rational c);
bool triangular_inequality(const Rational x, const Rational y, const Rational z);
