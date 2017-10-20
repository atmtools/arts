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
  \file   m_telsem.cc

  \brief  This file contains functions to read TELSEM atlases.
*/

#include "matpackI.h"
#include "mystring.h"
#include "file.h"
#include "telsem.h"


/* Workspace method: Doxygen documentation will be auto-generated */
void telsem_atlasesReadAscii(ArrayOfTelsemAtlas& telsem_atlases,
                             const String& directory,
                             const String& filename_pattern,
                             const Verbosity& verbosity)
{
    CREATE_OUT2;
    const Index imonth = filename_pattern.find("@MM@");
    if (imonth < 0)
    {
        ostringstream os;
        os << "Substring '@MM@' not found in filename_pattern for" << std::endl
           << "month number replacement: "
           << filename_pattern;

    }

    telsem_atlases.resize(12);
    for (Index i = 1; i <= 12; i++)
    {
        std::ifstream is;
        ostringstream month;
        if (i < 10) month << 0;
        month << i;
        String this_filename = filename_pattern;
        this_filename.replace(imonth, 4, month.str());
        this_filename = directory + '/' + this_filename;

        out2 << "Reading TELSEM atlas: " << this_filename << '\n';
        open_input_file(is, this_filename);
        telsem_read_ascii(is, telsem_atlases[i - 1]);
        telsem_atlases[i - 1].month = i;
    }

    std::ifstream is;
    String corr_filename = directory + '/' + "correlations";
    out2 << "Reading correlations: " << corr_filename << '\n';
    open_input_file(is, corr_filename);
    Tensor3 correlation(10, 7, 7);
    String s;
    for (Index i = 0; i < 10; i++)
    {
        std::getline(is, s);
        for (Index j = 0; j < 7; j++)
        {
            for (Index k = 0; k < 7; k++)
            {
                is >> correlation(i, j, k);
                if (is.fail())
                    throw std::runtime_error("Error reading correlation.");
            }
            std::getline(is, s);
        }
    }

    for (Index i = 0; i < 12; i++)
    {
        telsem_atlases[i].correl = correlation;
    }
}

