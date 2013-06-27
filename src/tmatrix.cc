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
#include <cmath>
#include "complex.h"
#include "messages.h"
#include "matpackI.h"


void phasmat(Matrix& z,
             const Complex& s11,
             const Complex& s12,
             const Complex& s21,
             const Complex& s22);

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

    void tmatrix_(const Numeric& rat,
                  const Numeric& axi,
                  const Index&   np,
                  const Numeric& lam,
                  const Numeric& eps,
                  const Numeric& mrr,
                  const Numeric& mri,
                  const Numeric& ddelt,
                  Index&         nmax,
                  Numeric&       csca,
                  Numeric&       cext,
                  const Index&   quiet,
                  char*          errmsg);

    void ampl_(const Index& nmax,
              const Numeric& lam,
              const Numeric& thet0,
              const Numeric& thet,
              const Numeric& phi0,
              const Numeric& phi,
              const Numeric& alpha,
              const Numeric& beta,
              Complex&       s11,
              Complex&       s12,
              Complex&       s21,
              Complex&       s22);
}


void tmatrix_ampld_test()
{
    Verbosity verbosity(2, 2, 2);
    CREATE_OUT0;

    out0 << "======================================================\n";
    out0 << "Test for nonspherical particles in a fixed orientation\n";
    out0 << "Output should match 3rdparty/tmatrix/tmatrix_ampld.ref\n";
    out0 << "======================================================\n";

    // Same inputs as in example included in original ampld.lp.f
    Numeric rat = 1.;
    Numeric axi = 10.;
    Index   np = -1;
    Numeric lam = acos(-1.)*2.;
    Numeric eps = 0.5;
    Numeric mrr = 1.5;
    Numeric mri = 0.02;
    Numeric ddelt = 0.001;

    Index   quiet = 0;

    // Output variables
    Index   nmax;
    Numeric csca;
    Numeric cext;
    char errmsg[1024] = "";

    tmatrix_(rat, axi, np, lam, eps, mrr, mri, ddelt, nmax, csca, cext,
             quiet, errmsg);

    out0 << "nmax: " << nmax << "\n";
    out0 << "csca: " << csca << "\n";
    out0 << "cext: " << cext << "\n";

    out0 << "Error message: " << (strlen(errmsg)?errmsg:"None") << "\n";

    // Same inputs as in example included in original ampld.lp.f
    Numeric alpha = 145.;
    Numeric beta = 52.;
    Numeric thet0 = 56.;
    Numeric thet = 65.;
    Numeric phi0 = 114.;
    Numeric phi = 128.;

    // Output variables
    Complex s11;
    Complex s12;
    Complex s21;
    Complex s22;
    ampl_(nmax, lam, thet0, thet, phi0, phi, alpha, beta,
          s11, s12, s21, s22);

    out0 << "AMPLITUDE MATRIX: \n";
    out0 << "s11: " << s11 << "\n";
    out0 << "s12: " << s12 << "\n";
    out0 << "s21: " << s21 << "\n";
    out0 << "s22: " << s22 << "\n";

    Matrix z;
    phasmat(z, s11, s12, s21, s22);

    out0 << "PHASE MATRIX: \n" << z << "\n";
}


void tmatrix_tmd_test()
{
    Verbosity verbosity(2, 2, 2);
    CREATE_OUT0;

    out0 << "======================================================\n";
    out0 << "Test for randomly oriented nonspherical particles\n";
    out0 << "Output should match 3rdparty/tmatrix/tmatrix_tmd.ref\n";
    out0 << "======================================================\n";

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

    // Output variables
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

void tmatrix_ampld_test()
{
    throw std::runtime_error("This version of ARTS was compiled without T-Matrix support.");
}

void tmatrix_tmd_test()
{
    throw std::runtime_error("This version of ARTS was compiled without T-Matrix support.");
}

#endif


/** Calculate Phase Matrix.

 Ported from the T-Matrix Fortran code in 3rdparty/tmatrix/ampld.lp.f

 \param[out] z    Phase Matrix
 \param[in]  s11  Calculated by ampl_
 \param[in]  s12  Calculated by ampl_
 \param[in]  s21  Calculated by ampl_
 \param[in]  s22  Calculated by ampl_
 \return

 \author Oliver Lemke
 */
void phasmat(Matrix& z,
             const Complex& s11,
             const Complex& s12,
             const Complex& s21,
             const Complex& s22)
{
    z.resize(4, 4);
    z(0, 0) =               0.5 * (  s11 * conj(s11) + s12 * conj(s12)
                                   + s21 * conj(s21) + s22 * conj(s22)).real();
    z(0, 1) =               0.5 * (  s11 * conj(s11) - s12 * conj(s12)
                                   + s21 * conj(s21) - s22 * conj(s22)).real();
    z(0, 2) =                     ( -s11 * conj(s12) - s22 * conj(s21)).real();
    z(0, 3) = (Complex(0.,  1.) * (  s11 * conj(s12) - s22 * conj(s21))).real();

    z(1, 0) =               0.5 * (  s11 * conj(s11) + s12 * conj(s12)
                                   - s21 * conj(s21) - s22 * conj(s22)).real();
    z(1, 1) =               0.5 * (  s11 * conj(s11) - s12 * conj(s12)
                                   - s21 * conj(s21) + s22 * conj(s22)).real();
    z(1, 2) =                     ( -s11 * conj(s12) + s22 * conj(s21)).real();
    z(1, 3) = (Complex(0.,  1.) * (  s11 * conj(s12) + s22 * conj(s21))).real();

    z(2, 0) =                     ( -s11 * conj(s21) - s22 * conj(s12)).real();
    z(2, 1) =                     ( -s11 * conj(s21) + s22 * conj(s12)).real();
    z(2, 2) =                     (  s11 * conj(s22) + s12 * conj(s21)).real();
    z(2, 3) = (Complex(0., -1.) * (  s11 * conj(s22) + s21 * conj(s12))).real();

    z(3, 0) = (Complex(0.,  1.) * (  s21 * conj(s11) + s22 * conj(s12))).real();
    z(3, 1) = (Complex(0.,  1.) * (  s21 * conj(s11) - s22 * conj(s12))).real();
    z(3, 2) = (Complex(0., -1.) * (  s22 * conj(s11) - s12 * conj(s21))).real();
    z(3, 3) =                     (  s22 * conj(s11) - s12 * conj(s21)).real();
}

