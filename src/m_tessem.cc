/* Copyright (C) 2016
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
  \file   m_tessem.cc

  \brief  This file contains functions that are adapted from TESSEM
  code which is used to calculate surface emissivity.
*/

#include "file.h"
#include "matpack_data.h"
#include "mystring.h"
#include "tessem.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void TessemNNReadAscii(TessemNN& net,
                       const String& net_file) {
  std::ifstream net_is;

  open_input_file(net_is, net_file);
  tessem_read_ascii(net_is, net);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void TestTessem(Vector& outvalues,
                const TessemNN& net,
                const Vector& invalues) {
  outvalues.resize(net.nb_outputs);
  tessem_prop_nn(outvalues, net, invalues);
}
