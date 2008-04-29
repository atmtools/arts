/* Copyright (C) 2003-2008 Stefan Buehler <sbuehler@ltu.se>

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
   \file   arts.cc

   This file contains global functions.

   \author Oliver Lemke
   \date   2003-05-07
*/

#include <cstdlib>
#include "arts.h"


/** This is the exit function of ARTS. Whenever arts has to be terminated
  at some point, call this function.

  You can call without any parameters, since the exit status then
  defaults to EXIT_FAILURE.  

  \param  status  Exit code. EXIT_FAILURE if omitted.
*/
void arts_exit (int status)
{
  exit (status);
}

