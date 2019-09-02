/* Copyright (C) 2017

   Jana Mendrok <jana.mendrok@gmail.com>
   Patrick Eriksson <patrick.eriksson@chalmers.se>
                      
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
  \file   psd.h
  \author Jana Mendrok Patrick Eriksson
  \date   2017-11-05 
  
  \brief  Internal functions of PSD type
*/

#ifndef psd_h
#define psd_h

#include "array.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "matpackVII.h"
#include "messages.h"
#include "optproperties.h"
#include "ppath.h"

void psd_cloudice_MH97(Vector& psd,
                       const Vector& diameter,
                       const Numeric& iwc,
                       const Numeric& t,
                       const bool noisy);

void psd_rain_A12(Vector& psd, const Vector& diameter, const Numeric& rwc);

void psd_rain_W16(Vector& psd, const Vector& diameter, const Numeric& rwc);

void psd_snow_F07(Vector& psd,
                  const Vector& diameter,
                  const Numeric& swc,
                  const Numeric& t,
                  const Numeric alpha,
                  const Numeric beta,
                  const String& regime);

#endif  //psd_h
