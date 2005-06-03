/* Copyright (C) 2004 Mattias Ekstrom <ekstrom@rss.chalmers.se>

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

    \author Mattias Ekstrom
*/

#ifndef jacobian_h
#define jacobian_h

#include <map>
#include <iostream>
#include <stdexcept>
#include "matpackI.h"
#include "array.h"
#include "mystring.h"
#include "make_array.h"
#include "bifstream.h"
#include "interpolation.h"
#include "logic.h"
#include "methods.h"
#include "ppath.h"
#include "agenda_class.h"

/** Contains the data for one retrieval quantity.
    \author Mattias Ekstrom */
class RetrievalQuantity {
public:

  /** Default constructor. Needed by make_array. */
  RetrievalQuantity() { /* Nothing to do here */ }

  /** Copy constructor. We need this, since operator= does not work
      correctly for Arrays. (Target Array has to be resized first.) */
  RetrievalQuantity(const RetrievalQuantity& x) :
    mmaintag(x.mmaintag),
    msubtag(x.msubtag),
    mmode(x.mmode),
    manalytical(x.manalytical),
    mperturbation(x.mperturbation),
    mgrids(x.mgrids)
  { /* Nothing left to do here. */ }

  /** Constructor that sets the values. */
  RetrievalQuantity(const String&             maintag,
                    const String&             subtag,
                    const String&             mode,
                    const Index&              analytical,
                    const Numeric&            perturbation,
                    const MakeArray<Vector>&  grids ) :
    mmaintag(maintag),
    msubtag(subtag),
    mmode(mode),
    manalytical(analytical),
    mperturbation(perturbation),
    mgrids(grids)
  {
    // With Matpack, initialization of mgrids from grids should work correctly.
  }

  /** Main tag. */
  const String& MainTag() const { return mmaintag; }
  void MainTag( const String& mt ) { mmaintag = mt; }
  /** Subtag. Eg. for gas species: O3, ClO. */
  const String& Subtag() const { return msubtag; }
  void Subtag( const String& st ) { msubtag = st; }
  /** Mode of retrieval quantity. Eg. "abs", "rel", "vmr" and "nd". */
  const String& Mode() const { return mmode; }
  void Mode( const String& m ) { mmode = m; }
  /** Boolean to make analytical calculations (if possible). */
  const Index& Analytical() const { return manalytical; }
  void Analytical( const Index& m ) { manalytical = m; }
  /** Size of perturbation used for perturbation calculations. */
  const Numeric& Perturbation() const { return mperturbation; }
  void Perturbation( const Numeric& p ) { mperturbation = p; }
  /** Grids. Definition grids for the jacobian, eg. p, lat and lon. */
  const ArrayOfVector& Grids() const { return mgrids; }
  void Grids( const ArrayOfVector& g ) { mgrids = g; }

private:

  String mmaintag;
  String msubtag;
  String mmode;
  Index manalytical;
  Numeric mperturbation;
  ArrayOfVector mgrids;
};

/** Output operator for RetrievalQuantity.

    \author Mattias Ekstrom */
ostream& operator << (ostream& os, const RetrievalQuantity& ot);

typedef Array<RetrievalQuantity> ArrayOfRetrievalQuantity;

//======================================================================
//             Functions related to calculation of Jacobian
//======================================================================

void agenda_append(       Agenda& agenda,
                    const String& methodname,
                    const String& keywordvalue);

bool check_retrieval_grids(       ArrayOfVector& grids,
                                  ostringstream& os,
                            const Vector&        p_grid,
                            const Vector&        lat_grid,
                            const Vector&        lon_grid,
                            const Vector&        p_retr,
                            const Vector&        lat_retr,
                            const Vector&        lon_retr,
                            const String&        p_retr_name,
                            const String&        lat_retr_name,
                            const String&        lon_retr_name,
                            const Index&         dim);

void get_perturbation_gridpos(      ArrayOfGridPos& gp,
                              const Vector&         atm_grid,
                              const Vector&         jac_grid,
                              const bool&           is_pressure);

void get_perturbation_limit(       ArrayOfIndex& limit,
                             const Vector&       pert_grid,
                             const Vector&       atm_limit);

void get_perturbation_range(       Range& range,
                             const Index& index,
                             const Index& length);

void perturbation_field_1d(       VectorView      field,
                            const ArrayOfGridPos& p_gp,
                            const Index&          p_pert_n,
                            const Range&          p_range,
                            const Numeric&        size,
                            const Index&          method);
                                
void perturbation_field_2d(       MatrixView      field,
                            const ArrayOfGridPos& p_gp,
                            const ArrayOfGridPos& lat_gp,
                            const Index&          p_pert_n,
                            const Index&          lat_pert_n,
                            const Range&          p_range,
                            const Range&          lat_range,
                            const Numeric&        size,
                            const Index&          method);
                                
void perturbation_field_3d(       Tensor3View     field,
                            const ArrayOfGridPos& p_gp,
                            const ArrayOfGridPos& lat_gp,
                            const ArrayOfGridPos& lon_gp,
                            const Index&          p_pert_n,
                            const Index&          lat_pert_n,
                            const Index&          lon_pert_n,
                            const Range&          p_range,
                            const Range&          lat_range,
                            const Range&          lon_range,
                            const Numeric&        size,
                            const Index&          method);
                                
void jacobian_from_path_to_rgrids(
         Tensor3&             diy_dx,
   ConstTensor3View           diy_dq,
   const Index&               atmosphere_dim,
   const Ppath&               ppath,
   ConstVectorView            ppath_p,
   const RetrievalQuantity&   jacobian_quantity );
                           
#endif // jacobian_h
