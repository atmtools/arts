/* Copyright (C) 2003 Oliver Lemke
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */

/**
  \file mpi.h

  This is header contains declarations needed for MPI support in ARTS.
  
  \author Oliver Lemke
  \date 2003-05-07
*/

#include "arts.h"

/* Define a macro to put around code that should only be used when
 * mpi support is active.
 */
#ifdef HAVE_MPI
#define MPI_ONLY(x) x
#else
#define MPI_ONLY(x)
#endif

#ifdef HAVE_MPI
#ifndef ARTS_MPI_H_INCLUDED
#define ARTS_MPI_H_INCLUDED

#include <mpi.h>

void mpi_startup (int &argc, char **&argv);

void mpi_shutdown ();

#endif
#endif

