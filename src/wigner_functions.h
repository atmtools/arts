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
#include "../3rdparty/wigner/wigxjpf/inc/wigxjpf.h"
#ifdef FAST_WIGNER_PATH_3J
#define DO_FAST_WIGNER 1
#include "../3rdparty/wigner/fastwigxj/inc/fastwigxj.h"
#else
#define DO_FAST_WIGNER 0
#endif

Numeric wigner3j(const Rational j1, const Rational j2, const Rational j3,
                 const Rational m1, const Rational m2, const Rational m3);

Numeric wigner6j(const Rational j1,const Rational j2,const Rational j3,
                 const Rational l1,const Rational l2,const Rational l3);

Numeric wigner9j(const Rational j11,const Rational j12,const Rational j13,
                 const Rational j21,const Rational j22,const Rational j23,
                 const Rational j31,const Rational j32,const Rational j33);

void ECS_wigner_CO2(Matrix& M, 
                    const ArrayOfRational& Jl, 
                    const ArrayOfRational& Ju, 
                    const Rational& ll, 
                    const Rational& lu, 
                    ConstVectorView G0, 
                    ConstVectorView population);
