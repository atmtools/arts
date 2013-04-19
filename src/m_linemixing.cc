/* Copyright (C) 2013
   Oliver Lemke  <olemke@core-dump.info>

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

#include "arts.h"
#include "absorption.h"

void line_mixing_o2Read(// WS Output:
                        ArrayOfVector& line_mixing_o2 _U_,
                        ArrayOfArrayOfIndex& line_mixing_o2_lut _U_,
                        // WS Input:
                        const ArrayOfArrayOfLineRecord& abs_lines_per_species _U_,
                        const ArrayOfArrayOfSpeciesTag& abs_species _U_,
                        const String& filename _U_,
                        const Verbosity& verbosity _U_)
{

}

