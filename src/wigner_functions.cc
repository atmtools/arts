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

// #include "wigxjpf/inc/wigxjpf.h"
#include "wigner_functions.h"
#include <stdexcept>
#include <sstream>
#include <algorithm>

/*!
  Run wigxjpf wig3jj for Rational symbol
  
  /                \
  |  j1   j2   j3  |
  |                |
  |  m1   m2   m3  |
  \                /
  
  See for definition: http://dlmf.nist.gov/34.2
*/
Numeric wigner3j(const Rational j1,const Rational j2,const Rational j3,
                 const Rational m1,const Rational m2,const Rational m3)
{
  int a = int((j1.Denom() == 2)? j1.Nom() : 2 * j1.Nom()),
      b = int((j2.Denom() == 2)? j2.Nom() : 2 * j2.Nom()), 
      c = int((j3.Denom() == 2)? j3.Nom() : 2 * j3.Nom()),
      d = int((m1.Denom() == 2)? m1.Nom() : 2 * m1.Nom()),
      e = int((m2.Denom() == 2)? m2.Nom() : 2 * m2.Nom()),
      f = int((m3.Denom() == 2)? m3.Nom() : 2 * m3.Nom());
  int g = std::max(std::max(std::max(std::max(std::max(a, b), c), d), e), f);
  double h;
  
  wig_table_init(g, 3);
  wig_temp_init(g);
            
  h = wig3jj(a, b, c, d, e, f);
  
  wig_temp_free();
  wig_table_free();
  
  return Numeric(h);
}


/*!
 R un wigxjpf wig3jj* for Rational symbol
 
 /                \
 |  j1   j2   j3  |
 <                >
 |  l1   l2   l3  |
 \                /
 
 See for definition: http://dlmf.nist.gov/34.4
*/
Numeric wigner6j(const Rational j1,const Rational j2,const Rational j3,
                 const Rational l1,const Rational l2,const Rational l3)
{
  int a = int((j1.Denom() == 2)? j1.Nom() : 2 * j1.Nom()),
      b = int((j2.Denom() == 2)? j2.Nom() : 2 * j2.Nom()), 
      c = int((j3.Denom() == 2)? j3.Nom() : 2 * j3.Nom()),
      d = int((l1.Denom() == 2)? l1.Nom() : 2 * l1.Nom()),
      e = int((l2.Denom() == 2)? l2.Nom() : 2 * l2.Nom()),
      f = int((l3.Denom() == 2)? l3.Nom() : 2 * l3.Nom());
  int g = std::max(std::max(std::max(std::max(std::max(a, b), c), d), e), f);
  double h;
  
  wig_table_init(g, 6);
  wig_temp_init(g);
  
  h = wig6jj(a, b, c, d, e, f);
  
  wig_temp_free();
  wig_table_free();
  
  return Numeric(h);
}


Numeric wigner9j(const Rational j11,const Rational j12,const Rational j13,
                 const Rational j21,const Rational j22,const Rational j23,
                 const Rational j31,const Rational j32,const Rational j33)
{
  int a = int((j11.Denom() == 2)? j11.Nom() : 2 * j11.Nom()),
      b = int((j12.Denom() == 2)? j12.Nom() : 2 * j12.Nom()), 
      c = int((j13.Denom() == 2)? j13.Nom() : 2 * j13.Nom()),
      d = int((j21.Denom() == 2)? j21.Nom() : 2 * j21.Nom()),
      e = int((j22.Denom() == 2)? j22.Nom() : 2 * j22.Nom()),
      f = int((j23.Denom() == 2)? j23.Nom() : 2 * j23.Nom()),
      g = int((j31.Denom() == 2)? j31.Nom() : 2 * j31.Nom()),
      h = int((j32.Denom() == 2)? j32.Nom() : 2 * j32.Nom()),
      i = int((j33.Denom() == 2)? j33.Nom() : 2 * j33.Nom());
  int j = std::max(std::max(std::max(std::max(std::max(std::max(std::max(std::max(a, b), c), d), e), f), g), h), i);
  double k;
  
  wig_table_init(j, 9);
  wig_temp_init(j);
  
  k = wig9jj(a, b, c, d, e, f, g, h, i);
  
  wig_temp_free();
  wig_table_free();
  
  return Numeric(k);
}


Numeric ECS_wigner(Rational L, Rational Nl, Rational Nk, 
                   Rational Jk_lower, Rational Jl_lower, 
                   Rational Jk_upper, Rational Jl_upper)
{
    const Numeric A1 = wigner3j(Nl,Nk,L,0,0,0);
    if( A1 == 0)
        return 0;
    
    const Numeric A2 = wigner6j(L,Jk_upper,Jl_upper,1,Nl,Nk);
    if( A2 == 0)
        return 0;
    
    const Numeric A3 = wigner6j(L,Jk_lower,Jl_lower,1,Nl,Nk);
    if( A3 == 0)
        return 0;
    
    const Numeric A4 = wigner6j(L,Jk_upper,Jl_upper,1,Jl_lower,Jk_lower);
    if( A4 == 0)
        return 0;
    
    return A1*A2*A3*A4;
}
