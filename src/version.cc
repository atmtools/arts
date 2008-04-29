/* Copyright (C) 2000-2008
   Stefan Buehler <sbuehler@ltu.se>
   Patrick Eriksson <patrick@rss.chalmers.se>

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

/**
   The ARTS running version number.
   Please always increase this before you do a `cvs
   commit', however minor your change may be.
 */

#include "arts.h"
#include "mystring.h"

// You now have to set the full version number (including subversion)
// in file configure.in in the top-level ARTS directory. Don't forget
// to run `autogen.sh' after editing configure.in!

String full_name  = static_cast<String>(PACKAGE)
  + "-" + VERSION;

