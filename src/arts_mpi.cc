/* Copyright (C) 2003 Oliver Lemke
 * * This program is free software; you can redistribute it and/or
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
  \file arts_mpi.cc

  Implementation of the MPI interface for arts.

  \author Oliver Lemke <olemke@uni-bremen.de>
  \date 2003-05-14
  */

#include "arts.h"

#ifdef HAVE_MPI

#include "arts_mpi.h"

/** Global MPI Manager */
MpiManager mpi_manager;


MpiManager::MpiManager () : initialized (false)
{
}


MpiManager::~MpiManager ()
{
  if (initialized)
    MPI::Finalize ();
}


void MpiManager::startup (int &argc, char **&argv)
{
  MPI::Init (argc, argv);
  initialized = true;
}

#endif

