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

#include <iostream>
using namespace std;

#include "arts.h"
#include "arts_mpi.h"
#include "matpackI.h"


#ifndef HAVE_MPI
int main()
{
  cerr << "For this example to work, you have to enable mpi support." << endl
    << "Rerun configure or autogen.sh with --with-mpi." << endl;

  return 0;
}
#else
int main(int argc, char *argv[])
{
  mpi_manager.startup (argc, argv);

  if (mpi_manager.get_nprocs () < 2)
    {
      cerr << "This example must be run on at least 2 processors/nodes."
        << endl;
      return 1;
    }

  cout << "nProcs: " << mpi_manager.get_nprocs () << endl;

  // Create data on master and send it to all nodes
  if (mpi_manager.get_rank () == 0)
    {
      const Index n = 100;
      Matrix a (n, n);

      Numeric k = 1.;
      for (Index i = 0; i < n; i++)
        for (Index j = 0; j < n; j++, k++)
          a (i, j) = k;

      // TODO: Broadcast data
    }
  else
    {
      // TODO: Receive data
    }

  return 0;
}
#endif // HAVE_MPI
