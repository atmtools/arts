/* Copyright (C) 2002-2007 Stefan Buehler <sbuehler@ltu.se>

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

/*!
  \file   complex.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2002-12-19
  
  \brief  A class implementing complex numbers for ARTS.
*/

#include "complex.h"

complex<float> operator+ (const double &d, const complex<float> &c)
{
  return (float(d) + c);
}

complex<float> operator* (const double &d, const complex<float> &c)
{
  return (float(d) * c);
}


complex<float> operator+ (const complex<float> &c, const double &d)
{
  return (c + float(d));
}

complex<float> operator* (const complex<float> &c, const double &d)
{
  return (c * float(d));
}



complex<double> operator+ (const float &f, const complex<double> &c)
{
  return (double(f) + c);
}

complex<double> operator* (const float &f, const complex<double> &c)
{
  return (double(f) * c);
}

complex<double> operator+ (const complex<double> &c, const float &f)
{
  return (c + double(f));
}

complex<double> operator* (const complex<double> &c, const float &f)
{
  return (c * double(f));
}

