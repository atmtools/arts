/* Copyright (C) 2013 Oliver Lemke
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*!
  \file   tmatrix.cc
  \author Oliver Lemke
  \date   2013-06-25
  
  \brief  Implementation of the T-Matrix interface.
*/

#include "tmatrix.h"
#include <stdexcept>
#include <cstring>
#include "messages.h"
#include "matpackI.h"


#ifdef ENABLE_TMATRIX

extern "C" {
    // T-matrix code for randomly oriented nonspherical particles
    //
    void tmd_(const Numeric& rat,
              const Index&   ndistr,
              const Numeric& axmax,
              const Index&   npnax,
              const Numeric& b,
              const Numeric& gam,
              const Index&   nkmax,
              const Numeric& eps,
              const Index&   np,
              const Numeric& lam,
              const Numeric& mrr,
              const Numeric& mri,
              const Numeric& ddelt,
              const Index&   npna,
              const Index&   ndgs,
              const Numeric& r1rat,
              const Numeric& r2rat,
              const Index&   quiet,
              Numeric& reff,  // OUT
              Numeric& veff,  // OUT
              Numeric& cext,  // OUT
              Numeric& csca,  // OUT
              Numeric& walb,  // OUT
              Numeric& asymm, // OUT
              Numeric* f11,   // Array npna
              Numeric* f22,   // Array npna
              Numeric* f33,   // Array npna
              Numeric* f44,   // Array npna
              Numeric* f12,   // Array npna
              Numeric* f34,   // Array npna
              char*    errmsg);
}


void tmatrix_tmd_test()
{
//    Numeric rat = 1.;
//    Index   ndistr = 4;
//    Numeric axmax = 10.;
//    Index   npnax = 1;
//    Numeric b = 0.1;
//    Numeric gam = 0.5;
//    Index   nkmax = -1;
//    Numeric eps = 0.33;
//    Index   np = -1;
//    Numeric lam = 1208.840556451613;
//    Numeric mrr = 1.78145992756;
//    Numeric mri = 0.00512840878218;
//    Numeric ddelt = 0.001;
//    Index   npna = 19;
//    Index   ndgs = 2;
//    Numeric r1rat = 0.9999999;
//    Numeric r2rat = 1.0000001;

    // Same inputs as in example included in original tmd.lp.f
    Numeric rat = 0.5;
    Index   ndistr = 3;
    Numeric axmax = 1.;
    Index   npnax = 2;
    Numeric b = 0.1;
    Numeric gam = 0.5;
    Index   nkmax = 5;
    Numeric eps = 2;
    Index   np = -1;
    Numeric lam = 0.5;
    Numeric mrr = 1.53;
    Numeric mri = 0.008;
    Numeric ddelt = 0.001;
    Index   npna = 19;
    Index   ndgs = 2;
    Numeric r1rat = 0.89031;
    Numeric r2rat = 1.56538;

    Index   quiet = 0;

    Numeric reff;
    Numeric veff;
    Numeric cext;
    Numeric csca;
    Numeric walb;
    Numeric asymm;
    Vector f11(npna, 0.);
    Vector f22(npna, 0.);
    Vector f33(npna, 0.);
    Vector f44(npna, 0.);
    Vector f12(npna, 0.);
    Vector f34(npna, 0.);
    char errmsg[1024] = "";

    tmd_(rat,
         ndistr,
         axmax,
         npnax,
         b,
         gam,
         nkmax,
         eps,
         np,
         lam,
         mrr,
         mri,
         ddelt,
         npna,
         ndgs,
         r1rat,
         r2rat,
         quiet,
         reff,
         veff,
         cext,
         csca,
         walb,
         asymm,
         f11.get_c_array(),
         f22.get_c_array(),
         f33.get_c_array(),
         f44.get_c_array(),
         f12.get_c_array(),
         f34.get_c_array(),
         errmsg);

    Verbosity verbosity(2,2,2);
    CREATE_OUT0;
    out0 << "reff: " << reff << "\n";
    out0 << "veff: " << veff << "\n";
    out0 << "cext: " << cext << "\n";
    out0 << "csca: " << csca << "\n";
    out0 << "walb: " << walb << "\n";
    out0 << "asymm: " << asymm << "\n";
    out0 << "f11: " << f11 << "\n";
    out0 << "f22: " << f22 << "\n";
    out0 << "f33: " << f33 << "\n";
    out0 << "f44: " << f44 << "\n";
    out0 << "f12: " << f12 << "\n";
    out0 << "f34: " << f34 << "\n";

    out0 << "Error message: " << (strlen(errmsg)?errmsg:"None") << "\n";
}

#else

void tmatrix_tmd_test()
{
    throw std::runtime_error("This version of ARTS was compiled without T-Matrix support.");
}

#endif
