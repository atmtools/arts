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
  \file arts_mpi.h

  This header contains declarations needed for MPI support in ARTS.

  If arts is configured without MPI support, all this code is omitted.
  
  \author Oliver Lemke
  \date 2003-05-07
*/

#include "arts.h"

/* Define a macro that should be put around code which is only used when
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

/**
  The MPI Manager class. This class is initializes the MPI interface
  and provides access to information about the current MPI status.

  \author Oliver Lemke <olemke@core-dump.info>
  \date 2003-05-14
  */
class MpiManager
{
public:
  /** Default constructor */
  MpiManager ();

  /** Destructor. Calls mpi finalize if the MPI interface was initialized. */
  virtual ~MpiManager ();

  /** Starts up the MPI interface. */
  void startup (int &argc, char **&argv);

  /** Returns the number of processes used for the current arts run.

    \return Number of processes involved.
   */
  int get_nprocs () { return nprocs; }

  /** Returns the MPI rank of the current process.

    \return Rank of the current process.
   */
  int get_rank () { return rank; }

private:
  /** Status of MPI. True if MPI initialization was successful. */
  bool initialized;

  /** Number of processes used for the current arts run. */
  int  nprocs;

  /** MPI rank of the current process. */
  int  rank;
};


extern MpiManager mpi_manager;

#endif // ARTS_MPI_H_INCLUDED

#endif // HAVE_MPI

