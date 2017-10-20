/* Copyright (C) 2017
   Oliver Lemke <olemke@core-dump.info>

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
  \file   telsem.cc

  \brief  This file contains functions to handle TELSEM 2 atlas data.
*/

#include <cmath>
#include "telsem.h"

void telsem_calc_correspondence(TelsemAtlas& ta)
{
    ta.correspondence.resize(660066);
    ta.correspondence = NAN;
    for (Index j = 0; j < ta.ndat; j++)
    {
        ta.correspondence[(Index) ta.cellnum[j]] = (Numeric)j;
    }
}


void telsem_read_ascii(std::istream& is, TelsemAtlas& ta)
{
    ta.name = "ssmi_mean_emis_climato";
    ta.nchan = 7;
    ta.dlat = 0.25;
    is >> ta.ndat;
    ta.emis.resize(ta.ndat, ta.nchan);
    ta.emis = NAN;
    ta.emis_err.resize(ta.ndat, ta.nchan);
    ta.emis_err = NAN;
    ta.class1.resize(ta.ndat);
    ta.class1 = NAN;
    ta.class2.resize(ta.ndat);
    ta.class2 = NAN;
    ta.cellnum.resize(ta.ndat);
    ta.cellnum = NAN;

    Index cellnum;
    Index cur_class1;
    Index cur_class2;
    Vector ssmi(2 * ta.nchan);
    Index ipos = -1;
    for (Index j = 0; j < ta.ndat; j++)
    {
        is >> cellnum;
        if (is.fail())
            throw std::runtime_error("Error reading cellnum.");
        for (Index nssmi = 0; nssmi < 2 * ta.nchan; nssmi++)
        {
            is >> ssmi[nssmi];
            if (is.fail())
                throw std::runtime_error("Error reading emissivity.");
        }

        is >> cur_class1 >> cur_class2;
        if (is.fail())
            throw std::runtime_error("Error reading classes.");
        if (cur_class1 > 0 && cur_class2 > 0 && ipos < ta.ndat)
        {
            ipos++;
            for (Index i = 0; i < ta.nchan; i++)
            {
                ta.emis(ipos, i) = ssmi[i];
                ta.emis_err(ipos, i) = std::sqrt(ssmi[ta.nchan + i]);
            }
            ta.cellnum[ipos] = (Numeric)cellnum;
            ta.class1[ipos] = (Numeric)cur_class1;
            ta.class2[ipos] = (Numeric)cur_class2;
        }
    }
    telsem_calc_correspondence(ta);
}


std::ostream& operator<<(std::ostream& os, const TelsemAtlas& ta)
{
    os << ta.name << std::endl;

    return os;
}


