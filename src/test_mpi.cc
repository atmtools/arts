/* Copyright (C) 2003 Oliver Lemke <olemke@uni-bremen.de>

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

/*!
  \file   test_mpi.cc
  \author <olemke@uni-bremen.de>
  \date   2003-07-02

  \brief  Test the mpi interface.
*/

#include <cstdlib>
#include <iostream>
using namespace std;

#include "arts.h"
#include "arts_mpi.h"
#include "matpackI.h"

//#define NCOLS 10 * 10 * 20 * 100
#define NCOLS 2000
#define NROWS 100

#ifdef HAVE_MPI
void
mpi_broadcast (Matrix &m, Index source);

void
mpi_send (Matrix &m, Index source, Index dest);
#endif

#ifndef HAVE_MPI
int main()
{
  cerr << "For this example to work, you have to enable mpi support." << endl
    << "Rerun configure or autogen.sh with --with-mpi." << endl;

  return 0;
}

#else

int
main(int argc, char *argv[])
{
  double starttime;

  mpi_manager.startup (argc, argv);

  if (mpi_manager.get_nprocs () < 2)
    {
      cerr << "This example must be run on at least 2 processors/nodes."
        << endl;
      return 1;
    }

  system ("hostname");

  Matrix a;

  // Initialize data on master
  if (mpi_manager.get_rank () == 0)
    {
      const Index c = NCOLS, r = NROWS;

      cout << "nProcs: " << mpi_manager.get_nprocs () << endl;

      a.resize (r, c);

      Numeric k = 1.;
      for (Index i = 0; i < r; i++)
        for (Index j = 0; j < c; j++, k++)
          a (i, j) = k;

    }

  if (mpi_manager.get_rank () == 0)
    starttime = MPI::Wtime ();

  // Send matrix from master (0) to first client (1)
  mpi_send (a, 0, 1);

  // Output timing info on master
  if (mpi_manager.get_rank () == 0)
    cout << "Time: " << MPI::Wtime () - starttime << " secs" << endl;

  cout << "Rank: " << mpi_manager.get_rank ();
  cout << " - nCols: " << a.ncols ();
  cout << " - nRows: " << a.nrows ();
  cout << " - a (0, 50) = " << a (0, 50) << endl;

  if (mpi_manager.get_rank () == 0) starttime = MPI::Wtime ();

  // Broadcast matrix to all clients
  mpi_broadcast (a, 0);

  if (mpi_manager.get_rank () == 0)
    cout << "Time: " << MPI::Wtime () - starttime << " secs" << endl;

  cout << "Rank: " << mpi_manager.get_rank ();
  cout << " - nCols: " << a.ncols ();
  cout << " - nRows: " << a.nrows ();
  cout << " - a (0, 50) = " << a (0, 50) << endl;

  return 0;
}


void
mpi_broadcast (Matrix &m, Index source)
{
  Index n[2];
  n[0] = m.ncols ();
  n[1] = m.nrows ();

  if (mpi_manager.get_rank () == source)
    cout << "Broadcast matrix dimensions" << endl;

  MPI::COMM_WORLD.Bcast (n, 2, MPI::INT, source);

  // Resize the matrix on client nodes
  if (mpi_manager.get_rank () != source)
    {
      m.resize (n[1], n[0]);
    }

  if (mpi_manager.get_rank () == source)
    cout << "Broadcast matrix data" << endl;
  MPI::COMM_WORLD.Bcast (m.get_raw_data (), n[1] * n[0], MPI::DOUBLE, source);

}


void
mpi_send (Matrix &m, Index source, Index dest)
{
  Index n[2];
  n[0] = m.ncols ();
  n[1] = m.nrows ();

  if (mpi_manager.get_rank () == source)
    {
      cout << "Send matrix dimensions" << endl;
      MPI::COMM_WORLD.Send (n, 2, MPI::INT, dest, 0);
      cout << "Broadcast matrix data" << endl;
      MPI::COMM_WORLD.Send (m.get_raw_data (), n[1] * n[0], MPI::DOUBLE,
                            dest, 0);
    }
  else if (mpi_manager.get_rank () == dest)
    {
      MPI::Status status;
      MPI::COMM_WORLD.Recv (n, 2, MPI::INT, source, 0, status);
      m.resize (n[1], n[0]);
      MPI::COMM_WORLD.Recv (m.get_raw_data (), n[1] * n[0], MPI::DOUBLE,
                            source, 0, status);
    }

}

#endif // HAVE_MPI

