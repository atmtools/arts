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

/** \file
    Declarations required for the calculation of jacobians.

    \author Mattias Ekström
*/

#ifndef jacobian_h
#define jacobian_h

#include <stdexcept>
#include "matpackI.h"
#include "array.h"
#include "mystring.h"
#include "make_array.h"
#include "bifstream.h"

/** Contains the data for one retrieval quantity.
    \author Mattias Ekström */
class RetrievalQuantity {
public:

  /** Default constructor. Needed by make_array. */
  RetrievalQuantity() { /* Nothing to do here */ }

  /** Copy constructor. We need this, since operator= does not work
      correctly for Arrays. (Target Array has to be resized first.) */
  RetrievalQuantity(const RetrievalQuantity& x) :
    mmaintag(x.mmaintag),
    msubtag(x.msubtag),
    munit(x.munit),
    mmethod(x.mmethod),
    mperturbation(x.mperturbation),
    mgrids(x.mgrids),
    mspeciesindex(x.mspeciesindex),
    mjacobianindices(x.mjacobianindices)
  { /* Nothing left to do here. */ }

  /** Constructor that sets the values. */
  RetrievalQuantity(const String&             maintag,
                    const String&             subtag,
                    const String&             unit,
                    const Index&              method,
                    const Numeric&            perturbation,
                    const MakeArray<Vector>&  grids,
                    const Index&              speciesindex,
                    const Matrix&             jacobianindices) :
    mmaintag(maintag),
    msubtag(subtag),
    munit(unit),
    mmethod(method),
    mperturbation(perturbation),
    mgrids(grids),
    mspeciesindex(speciesindex),
    mjacobianindices(jacobianindices)
  {
    // With Matpack, initialization of mgrids from grids should work correctly.

    // FIXME: Do we want consistency checks? What kind of checks do we want??
  }

  /** Main tag. */
  const String&        MainTag()         const { return mmaintag; }
  /** Subtag. Eg. for gas species: O3, ClO. */
  const String&        Subtag()          const { return msubtag; }
  /** Unit of retrieval quantity. Eg. "rel", "vmr" and "nd". */
  const String&        Unit()            const { return munit; }
  /** Method of calculation. Perturbation (=0) or analytical (=1). */
  const Index&         Method()          const { return mmethod; }
  /** Size of perturbation used for perturbation calculations. */
  const Numeric&       Perturbation()    const { return mperturbation; }
  /** Grids. Definition grids for the jacobian, eg. p, lat and lon. */
  const ArrayOfVector& Grids()           const { return mgrids; }
  /** Species index (= index of species in *gas_species*).
      Should be = -1 for other quantities, eg. Temp, Pointing. */
  const Index&         SpeciesIndex()    const { return mspeciesindex; }
  /** Jacobian indices (= start and stop columns in jacobian). */
  const Matrix&        JacobianIndices() const { return mjacobianindices; }

private:

  String mmaintag;
  String msubtag;
  String munit;
  Index mmethod;
  Numeric mperturbation;
  ArrayOfVector mgrids;
  Index mspeciesindex;
  Matrix mjacobianindices;
};

/** Output operator for RetrievalQuantity.

    \author Mattias Ekström */
ostream& operator << (ostream& os, const RetrievalQuantity& ot);

typedef Array<RetrievalQuantity> ArrayOfRetrievalQuantity;

#endif // jacobian_h
