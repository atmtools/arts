/* Copyright (C) 2004 Mattias Ekström <ekstrom@rss.chalmers.se>

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



/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   m_jacobian.cc
  \author Mattias Ekström <ekstrom@rss.chalmers.se>
  \date   2004-09-14

  \brief  Workspace functions related to the jacobian.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <string>
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "math_funcs.h"
#include "messages.h"
#include "jacobian.h"


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! jacobianAddPointing
/*!
   See the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date   2004-09-14
*/
void jacobianAddPointing(// WS Output:
                         ArrayOfRetrievalQuantity&  jq,
                         // WS Input:
                         const Sparse&              jac,
                         // Control Parameters:
                         const Numeric&             dza,
                         const Index&               poly_order)
{
  // Check that the jacobian matrix is empty. Otherwise it is either
  // not initialised or it is closed.
  if (jac.nrows()!=0 && jac.ncols()!=0)
  {
    ostringstream os;
    os << "The Jacobian matrix is not initialised correctly or closed.\n"
       << "New retrieval quantities can not be added at this point.\n";
    throw runtime_error(os.str());
  }

  // Check that poly_order is positive
  if (poly_order<0)
    throw runtime_error("The polynomial order has to be positive.");

  // Define subtag here to easily expand function.
  String subtag = "za offset";

  // Check that this type of pointing is not already included in the
  // jacobian.
  for (Index it=0; it<jq.nelem(); it++)
  {
    if (jq[it].MainTag()=="Pointing" && jq[it].Subtag()==subtag)
    {
      ostringstream os;
      os << "A zenith angle pointing offset is already included in\n"
         << "*jacobian_quantities*.";
      throw runtime_error(os.str());
    }
  }

  // Create the new retrieval quantity
  RetrievalQuantity rq = RetrievalQuantity();
  rq.MainTag("Pointing");
  rq.Subtag(subtag);
  // Unit?
  rq.Method(0);
  rq.Perturbation(dza);
  Vector grid(0,poly_order,1);
  ArrayOfVector grids(1, grid);
  rq.Grids(grids);

  // Add it to the *jacobian_quantities*
  jq.push_back(rq);

  out2 << "  Adding zenith angle pointing offset to *jacobian_quantities*\n"
       << "  with perturbation size " << dza << "\n";
}


//! jacobianClose
/*!
   See the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date   2004-09-19
*/
void jacobianClose(// WS Output:
                   Sparse&                    jacobian,
                   ArrayOfRetrievalQuantity&  jacobian_quantities,
                   // WS Input:
                   const Matrix&              sensor_pos,
                   const Sparse&              sensor_response)
{
  // Check that *jacobian* has been initialised
  if (jacobian.nrows()!=0 && jacobian.ncols()!=0)
    throw runtime_error("The Jacobian matrix has not been initialised.");

  // Make sure that the array is not empty
  if (jacobian_quantities.nelem()==0)
    throw runtime_error(
      "No retrieval quantities has been added to *jacobian_quantities*");

  // Check that sensor_pol and sensor_response has been initialised
  if (sensor_pos.nrows()==0)
  {
    ostringstream os;
    os << "The number of rows in *sensor_pos* is zero, i.e. no measurement\n"
       << "blocks has been defined. This has to be done before calling\n"
       << "jacobianClose.";
    throw runtime_error(os.str());
  }
  if (sensor_response.nrows()==0)
  {
    ostringstream os;
    os << "The sensor has either to be defined or turned off before calling\n"
       << "jacobianClose.";
    throw runtime_error(os.str());
  }

  // Loop over retrieval quantities, set JacobianIndices
  Index nrows = sensor_pos.nrows()*sensor_response.nrows();
  Index ncols = 0;
  for (Index it=0; it<jacobian_quantities.nelem(); it++)
  {
    // Store start jacobian index
    Vector indices(2);
    indices[0] = ncols;

    // Count total number of field points, i.e. product of grid lengths
    Index cols = 1;
    ArrayOfVector grids = jacobian_quantities[it].Grids();
    for (Index jt=0; jt<grids.nelem(); jt++)
      cols *= grids[jt].nelem();
    ncols += cols;

    // Store stop index
    indices[1] = ncols-1;
    jacobian_quantities[it].JacobianIndices(indices);
  }
  
  // Resize *jacobian*
  jacobian.resize( nrows, ncols);
}


//! jacobianInit
/*!
   See the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date   2004-09-14
*/
void jacobianInit(// WS Output:
                  Sparse&                    jacobian,
                  ArrayOfRetrievalQuantity&  jacobian_quantities)
{
  // FIXME: Check if array and sparse has already been initialised?
  // If so, what to do? Error/warning

  // Resize arrays and sparse to zero.
  jacobian_quantities.resize(0);
  jacobian.resize(0,0);

  out2 <<
    "  Initialising *jacobian* and *jacobian_quantities*.\n";
}

