/* Copyright (C) 2003-2007 Oliver Lemke <olemke@core-dump.info>
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

/**
  \file arts_mpi.cc

  Implementation of the MPI interface for arts.

  \author Oliver Lemke <olemke@core-dump.info>
  \date 2003-05-14
  */

#include "arts.h"

#ifdef HAVE_MPI

#include "arts_mpi.h"
#include <iostream>

/** Global MPI Manager */
MpiManager mpi_manager;


MpiManager::MpiManager () : initialized (false), nprocs (0), rank (0)
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

  nprocs = MPI::COMM_WORLD.Get_size ();
  rank = MPI::COMM_WORLD.Get_rank ();
}

#endif

