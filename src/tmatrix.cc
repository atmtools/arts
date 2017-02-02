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
#include "make_vector.h"
#include "math_funcs.h"
#include "optproperties.h"


void calc_phamat(Matrix& z,
                 const Index& nmax,
                 const Numeric& lam,
                 const Numeric& thet0,
                 const Numeric& thet,
                 const Numeric& phi0,
                 const Numeric& phi,
                 const Numeric& beta,
                 const Numeric& alpha);

void ampmat_to_phamat(Matrix& z,
                      const Complex& s11,
                      const Complex& s12,
                      const Complex& s21,
                      const Complex& s22);

void integrate_phamat_alpha10(Matrix& phamat,
                              const Index& nmax,
                              const Numeric& lam,
                              const Numeric& thet0,
                              const Numeric& thet,
                              const Numeric& phi0,
                              const Numeric& phi,
                              const Numeric& beta,
                              const Numeric& alpha_1,
                              const Numeric& alpha_2);

void integrate_phamat_alpha6(Matrix& phamat,
                             const Index& nmax,
                             const Numeric& lam,
                             const Numeric& thet0,
                             const Numeric& thet,
                             const Numeric& phi0,
                             const Numeric& phi,
                             const Numeric& beta,
                             const Numeric& alpha_1,
                             const Numeric& alpha_2);

void integrate_phamat_theta0_phi10(Matrix& phamat,
                                   const Index& nmax,
                                   const Numeric& lam,
                                   const Numeric& thet0_1,
                                   const Numeric& thet0_2,
                                   const Numeric& thet,
                                   const Numeric& phi0,
                                   const Numeric& phi_1,
                                   const Numeric& phi_2,
                                   const Numeric& beta,
                                   const Numeric& alpha);

void integrate_phamat_theta0_phi_alpha6(Matrix& phamat,
                                        const Index& nmax,
                                        const Numeric& lam,
                                        const Numeric& thet0_1,
                                        const Numeric& thet0_2,
                                        const Numeric& thet,
                                        const Numeric& phi0,
                                        const Numeric& phi_1,
                                        const Numeric& phi_2,
                                        const Numeric& beta,
                                        const Numeric& alpha_1,
                                        const Numeric& alpha_2);


#ifdef ENABLE_TMATRIX
extern "C" {
#endif
    /** T-matrix code for randomly oriented nonspherical particles.

     This is the interface to the T-Matrix tmd Fortran subroutine.
     
     See 3rdparty/tmatrix/tmd.lp.f for the complete documentation of the
     T-Matrix codes.

     \param[in] rat    1 - particle size is specified in terms of the
                           equal-volume-sphere radius<br>
                       != 1 - particle size is specified in terms of the
                              equal-surface-area-sphere radius
     \param[in] ndistr Specifies the distribution of equivalent-sphere radii<br>
                       NDISTR = 1 - modified gamma distribution
                       [Eq. (40) of Ref. 7]<br>
                       AXI=alpha<br>
                       B=r_c<br>
                       GAM=gamma<br>
                       NDISTR = 2 - log-normal distribution
                       [Eq. 41) of Ref. 7]<br>
                       AXI=r_g<br>
                       B=[ln(sigma_g)]**2<br>
                       NDISTR = 3 - power law distribution
                       [Eq. (42) of Ref. 7]<br>
                       AXI=r_eff (effective radius)<br>
                       B=v_eff (effective variance)<br>
                       Parameters R1 and R2 (see below) are calculated
                       automatically for given AXI and B<br>
                       NDISTR = 4 - gamma distribution
                       [Eq. (39) of Ref. 7]<br>
                       AXI=a<br>
                       B=b<br>
                       NDISTR = 5 - modified power law distribution
                       [Eq. (24) in M. I. Mishchenko et al.,
                       Bidirectional reflectance of flat,
                       optically thick particulate laters: an efficient radiative
                       transfer solution and applications to snow and soil surfaces,
                       J. Quant. Spectrosc. Radiat. Transfer, Vol. 63, 409-432 (1999)].<br>
                       B=alpha
     \param[in] axmax  The code computes NPNAX size distributions of the same type
                       and with the same values of B and GAM in one run.
                       The parameter AXI varies from AXMAX to AXMAX/NPNAX in steps of
                       AXMAX/NPNAX.  To compute a single size distribution, use
                       NPNAX=1 and AXMAX equal to AXI of this size distribution.
     \param[in] npnax  See axmax above
     \param[in] b      See ndistr above
     \param[in] gam    See ndistr above
     \param[in] nkmax  NKMAX<=988 is such that NKMAX+2 is the
                       number of Gaussian quadrature points used in
                       integrating over the size distribution for particles with
                       AXI=AXMAX.  For particles with AXI=AXMAX-AXMAX/NPNAX,
                       AXMAX-2*AXMAX/NPNAX, etc. the number of Gaussian points
                       linearly decreases.
                       For the modified power law distribution, the number
                       of integration points on the interval [0,R1] is also
                       equal to NKMAX.
     \param[in] eps    Shape of the particles<br>
                       For spheroids NP=-1 and EPS is the ratio of the
                       horizontal to rotational axes.  EPS is larger than
                       1 for oblate spheroids and smaller than 1 for
                       prolate spheroids.<br>
                       For cylinders NP=-2 and EPS is the ratio of the
                       diameter to the length.<br>
                       For Chebyshev particles NP must be positive and
                       is the degree of the Chebyshev polynomial, while
                       EPS is the deformation parameter.
     \param[in] np     Shape of the particles (see eps above)
     \param[in] lam    Wavelength of light
     \param[in] mrr    Vector with real parts of refractive index
     \param[in] mri    Vector with imaginary parts of refractive index
     \param[in] ddelt  Accuracy of the computations
     \param[in] npna   Number of equidistant scattering angles (from 0
                       to 180 deg) for which the scattering matrix is
                       calculated
     \param[in] ndgs   Parameter controlling the number of division points
                       in computing integrals over the particle surface.
                       For compact particles, the recommended value is 2.
                       For highly aspherical particles larger values (3, 4,...)
                       may be necessary to obtain convergence.
                       The code does not check convergence over this parameter.
                       Therefore, control comparisons of results obtained with
                       different NDGS-values are recommended.
     \param[in] r1rat
     \param[in] r2rat
     \param[in] quiet  0 = Verbose output from Fortran code, 1 = silent
     \param[out] reff  Effective radius of the size distribution
     \param[out] veff  Effective variance of the size distribution
     \param[out] cext  Extinction cross section per particle
     \param[out] csca  Scattering cross section per particle
     \param[out] walb  Single scattering albedo
     \param[out] asymm Asymmetry parameter of the phase function
     \param[out] f11   Elements of the normalized scattering matrix
                       versus scattering angle
     \param[out] f22   Elements of the normalized scattering matrix
                       versus scattering angle
     \param[out] f33   Elements of the normalized scattering matrix
                       versus scattering angle
     \param[out] f44   Elements of the normalized scattering matrix
                       versus scattering angle
     \param[out] f12   Elements of the normalized scattering matrix
                       versus scattering angle
     \param[out] f34   Elements of the normalized scattering matrix
                       versus scattering angle
     \param[out] errmsg  Error message string from Fortran code
     \return

     \author Oliver Lemke
     */
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

    /** T-matrix code for nonspherical particles in a fixed orientation

     This is the interface to the T-Matrix tmatrix Fortran subroutine.
     It calculates extinction and scattering cross section per particle.
     The T-Matrix is calculated internally and used accessed later by
     ampl_() via common blocks.

     See 3rdparty/tmatrix/ampld.lp.f for the complete documentation of the
     T-Matrix codes.

     \param[in] rat    1 - particle size is specified in terms of the
                           equal-volume-sphere radius<br>
                       != 1 - particle size is specified in terms of the
                          equal-surface-area-sphere radius
     \param[in] axi    Equivalent-sphere radius [microns]
     \param[in] np     Shape of the particles<br>
                       For spheroids NP=-1 and EPS is the ratio of the
                       horizontal to rotational axes.  EPS is larger than
                       1 for oblate spheroids and smaller than 1 for
                       prolate spheroids.<br>
                       For cylinders NP=-2 and EPS is the ratio of the
                       diameter to the length.<br>
     \param[in] lam    Wavelength of light
     \param[in] eps    Shape of the particles (see np above)
     \param[in] mrr    Vector with real parts of refractive index
     \param[in] mri    Vector with imaginary parts of refractive index
     \param[in] ddelt  Accuracy of the computations
     \param[in] quiet  0 = Verbose output from Fortran code, 1 = silent
     \param[out] nmax  Iteration count
     \param[out] cext  Extinction cross section per particle
     \param[out] csca  Scattering cross section per particle
     \param[out] errmsg  Error message string from Fortran code
     */
    void tmatrix_(const Numeric& rat,
                  const Numeric& axi,
                  const Index&   np,
                  const Numeric& lam,
                  const Numeric& eps,
                  const Numeric& mrr,
                  const Numeric& mri,
                  const Numeric& ddelt,
                  const Index&   quiet,
                  Index&         nmax,
                  Numeric&       csca,
                  Numeric&       cext,
                  char*          errmsg);

    /** T-matrix code for nonspherical particles in a fixed orientation

     This is the interface to the T-Matrix ampl Fortran subroutine.
     It calculates the amplitude matrix.
     The T-Matrix is passed from tmatrix_() internally via common blocks.

     See 3rdparty/tmatrix/tmd.lp.f for the complete documentation of the
     T-Matrix codes.
     
     \param[in]  nmax   Iteration count. Calculated by tmatrix_()
     \param[in]  lam    Wavelength of incident light
     \param[in]  thet0  Zenith angle of the incident beam [degrees]
     \param[in]  thet   zenith angle of the scattered beam in degrees
     \param[in]  phi0   azimuth angle of the incident beam in degrees
     \param[in]  phi    azimuth angle of the scattered beam in degrees
     \param[in]  alpha  Euler angle (in degrees) specifying the orientation
                        of the scattering particle relative to the laboratory
                        reference frame
     \param[in]  beta   Euler angle (in degrees) specifying the orientation
                        of the scattering particle relative to the laboratory
                        reference frame
     \param[out] s11    
     */
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

    /** Perform orientation averaging

     This should be called after tmatrix_() for prolate particles.
     Data is passed from the tmatrix_() subroutine internally in the Fortran
     code via common blocks.

    \param[in]  nmax   Iteration count. Calculated by tmatrix_()
     */
    void avgtmatrix_(const Index& nmax);
#ifdef ENABLE_TMATRIX
}
#endif


// Define dummy functions that throw runtime errors if ARTS is
// compiled without T-Matrix support.
#ifndef ENABLE_TMATRIX

// T-matrix code for randomly oriented nonspherical particles.
void tmd_(const Numeric&,
          const Index&,
          const Numeric&,
          const Index&,
          const Numeric&,
          const Numeric&,
          const Index&,
          const Numeric&,
          const Index&,
          const Numeric&,
          const Numeric&,
          const Numeric&,
          const Numeric&,
          const Index&,
          const Index&,
          const Numeric&,
          const Numeric&,
          const Index&,
          Numeric&,
          Numeric&,
          Numeric&,
          Numeric&,
          Numeric&,
          Numeric&,
          Numeric*,
          Numeric*,
          Numeric*,
          Numeric*,
          Numeric*,
          Numeric*,
          char*)
{
    throw std::runtime_error("This version of ARTS was compiled without T-Matrix support.");
}


// T-matrix code for nonspherical particles in a fixed orientation
void tmatrix_(const Numeric&,
              const Numeric&,
              const Index&,
              const Numeric&,
              const Numeric&,
              const Numeric&,
              const Numeric&,
              const Numeric&,
              const Index&,
              Index&,
              Numeric&,
              Numeric&,
              char*)
{
    throw std::runtime_error("This version of ARTS was compiled without T-Matrix support.");
}


// T-matrix code for nonspherical particles in a fixed orientation
void ampl_(const Index&,
           const Numeric&,
           const Numeric&,
           const Numeric&,
           const Numeric&,
           const Numeric&,
           const Numeric&,
           const Numeric&,
           Complex&,
           Complex&,
           Complex&,
           Complex&)
{
    throw std::runtime_error("This version of ARTS was compiled without T-Matrix support.");
}

void avgtmatrix_(const Index&)
{
    throw std::runtime_error("This version of ARTS was compiled without T-Matrix support.");
}

#endif


void calc_phamat(Matrix& z,
                 const Index& nmax,
                 const Numeric& lam,
                 const Numeric& thet0,
                 const Numeric& thet,
                 const Numeric& phi0,
                 const Numeric& phi,
                 const Numeric& beta,
                 const Numeric& alpha)
{
    Complex s11;
    Complex s12;
    Complex s21;
    Complex s22;
    ampl_(nmax, lam, thet0, thet, phi0, phi, alpha, beta, s11, s12, s21, s22);

    ampmat_to_phamat(z, s11, s12, s21, s22);
}


/** Calculate phase matrix from amplitude matrix.

 Ported from the T-Matrix Fortran code in 3rdparty/tmatrix/ampld.lp.f

 \param[out] z    Phase Matrix
 \param[in]  s11  Calculated by ampl_()
 \param[in]  s12  Calculated by ampl_()
 \param[in]  s21  Calculated by ampl_()
 \param[in]  s22  Calculated by ampl_()
 \return

 \author Oliver Lemke
 */
void ampmat_to_phamat(Matrix& z,
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


static const Numeric GaussLeg6[][3] = {
    {0.23861918, 0.66120939, 0.93246951},
    {0.46791393, 0.36076157, 0.17132449}
};

static const Numeric GaussLeg10[][5] = {
    {0.14887434, 0.43339539, 0.67940957, 0.86506337, 0.97390653},
    {0.29552422, 0.26926672, 0.21908636, 0.14945135, 0.06667134}
};


/** Integrate phase matrix over particle orientation angle.

 Performs a ten point Gauss–Legendre integration over the orientation
 angle (first Euler angle) of the particles from alpha1 to alpha2.

 \param[out] phamat  Integrated phase matrix
 \param[in]  nmax    FIXME OLE
 \param[in]  lam     See ampl_()
 \param[in]  thet0   See ampl_()
 \param[in]  thet    See ampl_()
 \param[in]  phi0    See ampl_()
 \param[in]  phi     See ampl_()
 \param[in]  beta    See ampl_()
 \param[in]  alpha_1 See alpha in ampl_()
 \param[in]  alpha_2 See alpha in ampl_()

 \author Oliver Lemke
 */
void integrate_phamat_alpha10(Matrix& phamat,
                              const Index& nmax,
                              const Numeric& lam,
                              const Numeric& thet0,
                              const Numeric& thet,
                              const Numeric& phi0,
                              const Numeric& phi,
                              const Numeric& beta,
                              const Numeric& alpha_1,
                              const Numeric& alpha_2)
{
    const Numeric c = 0.5 * (alpha_2+alpha_1);
    const Numeric m = 0.5 * (alpha_2-alpha_1);

    phamat.resize(4, 4);
    phamat = 0.;
    Matrix z;

    for (Index i = 0; i < 5; ++i)
    {
        const Numeric abscissa = GaussLeg10[0][i];
        const Numeric weight = GaussLeg10[1][i];

        calc_phamat(z, nmax, lam, thet0, thet, phi0, phi, beta, c + m*abscissa);
        z *= weight;
        phamat += z;
        
        calc_phamat(z, nmax, lam, thet0, thet, phi0, phi, beta, c - m*abscissa);
        z *= weight;
        phamat += z;
    }
    phamat *= m;
}


/** Integrate phase matrix over particle orientation angle.

 Performs a six point Gauss–Legendre integration over the orientation
 angle (first Euler angle) of the particles from alpha1 to alpha2.

 \param[out] phamat  Integrated phase matrix
 \param[in]  nmax    FIXME OLE
 \param[in]  lam     See ampl_()
 \param[in]  thet0   See ampl_()
 \param[in]  thet    See ampl_()
 \param[in]  phi0    See ampl_()
 \param[in]  phi     See ampl_()
 \param[in]  beta    See ampl_()
 \param[in]  alpha_1 See alpha in ampl_()
 \param[in]  alpha_2 See alpha in ampl_()

 \author Oliver Lemke
 */
void integrate_phamat_alpha6(Matrix& phamat,
                             const Index& nmax,
                             const Numeric& lam,
                             const Numeric& thet0,
                             const Numeric& thet,
                             const Numeric& phi0,
                             const Numeric& phi,
                             const Numeric& beta,
                             const Numeric& alpha_1,
                             const Numeric& alpha_2)
{
    const Numeric c = 0.5 * (alpha_2+alpha_1);
    const Numeric m = 0.5 * (alpha_2-alpha_1);

    phamat.resize(4, 4);
    phamat = 0.;
    Matrix z;

    for (Index i = 0; i < 3; ++i)
    {
        const Numeric abscissa = GaussLeg6[0][i];
        const Numeric weight = GaussLeg6[1][i];

        calc_phamat(z, nmax, lam, thet0, thet, phi0, phi, beta, c + m*abscissa);
        z *= weight;
        phamat += z;
        
        calc_phamat(z, nmax, lam, thet0, thet, phi0, phi, beta, c - m*abscissa);
        z *= weight;
        phamat += z;
    }
    phamat *= m;
}


/** Integrate phase matrix over angles thet0 and phi.

 Performs a ten point Gauss–Legendre integration over the
 zenith angle of the incident beam (thet0) and the azimuth
 angle of the scattered beam (phi).

 \param[out] phamat  Integrated phase matrix
 \param[in]  nmax    FIXME OLE
 \param[in]  lam     See ampl_()
 \param[in]  thet0_1 See thet0 in ampl_()
 \param[in]  thet0_2 See thet0 in ampl_()
 \param[in]  thet    See ampl_()
 \param[in]  phi0    See ampl_()
 \param[in]  phi_1   See phi in ampl_()
 \param[in]  phi_2   See phi in ampl_()
 \param[in]  beta    See ampl_()
 \param[in]  alpha   See in ampl_()

 \author Oliver Lemke
 */
void integrate_phamat_theta0_phi10(Matrix& phamat,
                                   const Index& nmax,
                                   const Numeric& lam,
                                   const Numeric& thet0_1,
                                   const Numeric& thet0_2,
                                   const Numeric& thet,
                                   const Numeric& phi0,
                                   const Numeric& phi_1,
                                   const Numeric& phi_2,
                                   const Numeric& beta,
                                   const Numeric& alpha)
{
    extern const Numeric PI;

    phamat.resize(4, 4);
    phamat = 0.;
    Matrix z;

    const Numeric c_thet0 = 0.5 * (thet0_2+thet0_1);
    const Numeric m_thet0 = 0.5 * (thet0_2-thet0_1);
    const Numeric c_phi = 0.5 * (phi_2+phi_1);
    const Numeric m_phi = 0.5 * (phi_2-phi_1);

    for (Index t = 0; t < 5; ++t)
    {
        const Numeric abscissa_thet0 = GaussLeg10[0][t];
        const Numeric weight_thet0 = GaussLeg10[1][t];

        Matrix phamat_phi(4, 4, 0.);

        for (Index p = 0; p < 5; ++p)
        {
            const Numeric abscissa_phi = GaussLeg10[0][p];
            const Numeric weight_phi = GaussLeg10[1][p];

            Numeric this_thet0 = c_thet0 + m_thet0 * abscissa_thet0;
            calc_phamat(z, nmax, lam, this_thet0, thet, phi0, c_phi + m_phi * abscissa_phi, beta, alpha);
            z *= weight_phi * sin(this_thet0*PI/180.);
            phamat_phi += z;

            calc_phamat(z, nmax, lam, this_thet0, thet, phi0, c_phi - m_phi * abscissa_phi, beta, alpha);
            z *= weight_phi * sin(this_thet0*PI/180.);
            phamat_phi += z;

            this_thet0 = c_thet0 - m_thet0 * abscissa_thet0;
            calc_phamat(z, nmax, lam, this_thet0, thet, phi0, c_phi + m_phi * abscissa_phi, beta, alpha);
            z *= weight_phi * sin(this_thet0*PI/180.);
            phamat_phi += z;

            calc_phamat(z, nmax, lam, this_thet0, thet, phi0, c_phi - m_phi * abscissa_phi, beta, alpha);
            z *= weight_phi * sin(this_thet0*PI/180.);
            phamat_phi += z;
        }
        phamat_phi *= m_phi * weight_thet0;
        phamat += phamat_phi;
    }
    
    phamat *= m_thet0;
}


/** Integrate phase matrix over angles thet0, phi and alpha.

 Performs a six point Gauss–Legendre integration over the
 zenith angle of the incident beam (thet0), the azimuth
 angle of the scattered beam (phi) and the orientation
 angle (first Euler angle) of the particles.

 \param[out] phamat  Integrated phase matrix
 \param[in]  nmax    FIXME OLE
 \param[in]  lam     See ampl_()
 \param[in]  thet0_1 See thet0 in ampl_()
 \param[in]  thet0_2 See thet0 in ampl_()
 \param[in]  thet    See ampl_()
 \param[in]  phi0    See ampl_()
 \param[in]  phi_1   See phi in ampl_()
 \param[in]  phi_2   See phi in ampl_()
 \param[in]  beta    See ampl_()
 \param[in]  alpha_1 See alpha in ampl_()
 \param[in]  alpha_2 See alpha in ampl_()

 \author Oliver Lemke
 */
void integrate_phamat_theta0_phi_alpha6(Matrix& phamat,
                                        const Index& nmax,
                                        const Numeric& lam,
                                        const Numeric& thet0_1,
                                        const Numeric& thet0_2,
                                        const Numeric& thet,
                                        const Numeric& phi0,
                                        const Numeric& phi_1,
                                        const Numeric& phi_2,
                                        const Numeric& beta,
                                        const Numeric& alpha_1,
                                        const Numeric& alpha_2)
{
    extern const Numeric PI;

    Matrix z;
    Matrix phamat_phi(4, 4);


    const Numeric c_thet0 = 0.5 * (thet0_2+thet0_1);
    const Numeric m_thet0 = 0.5 * (thet0_2-thet0_1);
    const Numeric c_phi = 0.5 * (phi_2+phi_1);
    const Numeric m_phi = 0.5 * (phi_2-phi_1);

    phamat.resize(4, 4);
    phamat = 0.;

    for (Index t = 0; t < 3; ++t)
    {
        const Numeric abscissa_thet0 = GaussLeg6[0][t];
        const Numeric weight_thet0 = GaussLeg6[1][t];

        phamat_phi = 0.;

        for (Index p = 0; p < 3; ++p)
        {
            const Numeric abscissa_phi = GaussLeg6[0][p];
            const Numeric weight_phi = GaussLeg6[1][p];

            Numeric this_thet0 = c_thet0 + m_thet0 * abscissa_thet0;
            integrate_phamat_alpha6(z, nmax, lam, this_thet0, thet, phi0, c_phi + m_phi * abscissa_phi, beta, alpha_1, alpha_2);
            z *= weight_phi * sin(this_thet0*PI/180.);
            phamat_phi += z;

            integrate_phamat_alpha6(z, nmax, lam, this_thet0, thet, phi0, c_phi - m_phi * abscissa_phi, beta, alpha_1, alpha_2);
            z *= weight_phi * sin(this_thet0*PI/180.);
            phamat_phi += z;

            this_thet0 = c_thet0 - m_thet0 * abscissa_thet0;
            integrate_phamat_alpha6(z, nmax, lam, this_thet0, thet, phi0, c_phi + m_phi * abscissa_phi, beta, alpha_1, alpha_2);
            z *= weight_phi * sin(this_thet0*PI/180.);
            phamat_phi += z;

            integrate_phamat_alpha6(z, nmax, lam, this_thet0, thet, phi0, c_phi - m_phi * abscissa_phi, beta, alpha_1, alpha_2);
            z *= weight_phi * sin(this_thet0*PI/180.);
            phamat_phi += z;
        }
        phamat_phi *= m_phi * weight_thet0;
        phamat += phamat_phi;
    }
    
    phamat *= m_thet0;
}


/** Calculate properties for randomly oriented particles.
 
 This is a simplified interface to the tmd_() function for randomly oriented
 particles based on the PyARTS function tmat_rnd
 
 \param[out] cext  Extinction cross section per particle
 \param[out] csca  Scattering cross section per particle
 \param[out] f11   See tmd_()
 \param[out] f22   See tmd_()
 \param[out] f33   See tmd_()
 \param[out] f44   See tmd_()
 \param[out] f12   See tmd_()
 \param[out] f34   See tmd_()
 \param[in] equiv_radius    See parameter axmax in tmd_()
 \param[in] aspect_ratio    See parameter eps in tmd_()
 \param[in] np              See tmd_()
 \param[in] lam             See tmd_()
 \param[in] ref_index_real  See parameter mrr in tmd_()
 \param[in] ref_index_imag  See parameter mri in tmd_()
 \param[in] precision       See parameter ddelt in tmd_()
 \param[in] nza             See parameter npna in tmd_()
 \param[in] ndgs            See parameter ndgs in tmd_()
 \param[in] quiet           See tmd_()
*/
void tmatrix_random_orientation(Numeric& cext,
                                Numeric& csca,
                                Vector& f11,
                                Vector& f22,
                                Vector& f33,
                                Vector& f44,
                                Vector& f12,
                                Vector& f34,
                                const Numeric equiv_radius,
                                const Numeric aspect_ratio,
                                const Index np,
                                const Numeric lam,
                                const Numeric ref_index_real,
                                const Numeric ref_index_imag,
                                const Numeric precision,
                                const Index nza,
                                const Index ndgs,
                                const Index quiet = 1)
{
    Numeric reff;
    Numeric veff;
    Numeric walb;
    Numeric asymm;

    char errmsg[1024] = "";

    f11.resize(nza); f11 = NAN;
    f22.resize(nza); f22 = NAN;
    f33.resize(nza); f33 = NAN;
    f44.resize(nza); f44 = NAN;
    f12.resize(nza); f12 = NAN;
    f34.resize(nza); f34 = NAN;

    // It is necessary to make sure that the Fortran code is not
    // called from different threads at the same time. The common
    // blocks are not threadsafe.
#pragma omp critical(tmatrix_code)
    tmd_(1.0,
         4,
         equiv_radius,
         1,
         0.1,
         1.0,
         -1,
         aspect_ratio,
         np,
         lam,
         ref_index_real,
         ref_index_imag,
         precision,
         nza,
         ndgs,
         0.9999999,
         1.0000001,
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

    if (strlen(errmsg))
    {
        std::ostringstream os;
        os << "T-Matrix code failed: " <<  errmsg;
        throw std::runtime_error(os.str());
    }
}


/** Calculate properties for particles in a fixed orientation.

 This is a simplified interface to the tmatrix_() function for randomly oriented
 particles based on the PyARTS function tmat_fxd

 \param[out] cext  Extinction cross section per particle
 \param[out] csca  Scattering cross section per particle
 \param[out] nmax  FIXME OLE
 \param[in] equiv_radius    See parameter axmax in tmd_()
 \param[in] aspect_ratio    See parameter eps in tmd_()
 \param[in] np              See tmd_()
 \param[in] lam             See tmd_()
 \param[in] ref_index_real  See parameter mrr in tmd_()
 \param[in] ref_index_imag  See parameter mri in tmd_()
 \param[in] precision       See parameter ddelt in tmd_()
 \param[in] quiet           See tmd_()
 */
void tmatrix_fixed_orientation(Numeric& cext,
                               Numeric& csca,
                               Index& nmax,
                               const Numeric equiv_radius,
                               const Numeric aspect_ratio,
                               const Index np,
                               const Numeric lam,
                               const Numeric ref_index_real,
                               const Numeric ref_index_imag,
                               const Numeric precision,
                               const Index quiet = 1)
{
    char errmsg[1024] = "";

    // It is necessary to make sure that the Fortran code is not
    // called from different threads at the same time. The common
    // blocks are not threadsafe.
#pragma omp critical(tmatrix_code)
    tmatrix_(1.,
             equiv_radius,
             np, 
             lam, 
             aspect_ratio,
             ref_index_real, 
             ref_index_imag, 
             precision,
             quiet,
             nmax, 
             csca, 
             cext, 
             errmsg);

    if (strlen(errmsg))
    {
        std::ostringstream os;
        os << "T-Matrix code failed: " <<  errmsg;
        throw std::runtime_error(os.str());
    }
}


void calcSingleScatteringDataProperties(SingleScatteringData& ssd,
                                        ConstMatrixView ref_index_real,
                                        ConstMatrixView ref_index_imag,
                                        const Numeric equiv_radius,
                                        const Index np,
                                        const Numeric aspect_ratio,
                                        const Numeric precision,
                                        const Index ndgs,
                                        const Index robust,
                                        const Index quiet )
{
    const Index nf = ssd.f_grid.nelem();
    const Index nT = ssd.T_grid.nelem();

    extern const Numeric PI;
    extern const Numeric SPEED_OF_LIGHT;

    if (ref_index_real.nrows() != nf)
        throw std::runtime_error(
          "Number of rows of refractive index real part must match ssd.f_grid.");
    if (ref_index_real.ncols() != nT)
        throw std::runtime_error(
          "Number of cols of refractive index real part must match ssd.T_grid.");
    if (ref_index_imag.nrows() != nf)
        throw std::runtime_error(
          "Number of rows of refractive index imaginary part must match ssd.f_grid.");
    if (ref_index_imag.ncols() != nT)
        throw std::runtime_error(
          "Number of cols of refractive index imaginary part must match ssd.T_grid.");

    // The T-Matrix code needs wavelength
    Vector lam(nf, SPEED_OF_LIGHT);
    lam /= ssd.f_grid;

    switch (ssd.ptype) {
        case PTYPE_TOTAL_RND:
        {
            const Index nza = ssd.za_grid.nelem();

            ssd.pha_mat_data.resize(nf, nT, nza, 1, 1, 1, 6);
            ssd.ext_mat_data.resize(nf, nT, 1, 1, 1);
            ssd.abs_vec_data.resize(nf, nT, 1, 1, 1);

            ssd.pha_mat_data = NAN;
            ssd.ext_mat_data = NAN;
            ssd.abs_vec_data = NAN;

            // Output variables
            Numeric cext = NAN;
            Numeric csca = NAN;
            Vector f11;
            Vector f22;
            Vector f33;
            Vector f44;
            Vector f12;
            Vector f34;
            Matrix mono_pha_mat_data(nza, 6, NAN);

            ostringstream os;
            os << "Calculation of SingleScatteringData properties failed for\n\n";
            bool anyfailed = false;
#pragma omp critical(tmatrix_ssp)
            for (Index f_index = 0; f_index < nf; ++f_index)
                for (Index T_index = 0; T_index < nT; ++T_index)
                {
                    bool thisfailed = false;
                    try {
                        tmatrix_random_orientation
                        (cext, csca,
                         f11, f22, f33, f44, f12, f34,
                         equiv_radius, aspect_ratio, np, lam[f_index],
                         ref_index_real(f_index, T_index),
                         ref_index_imag(f_index, T_index),
                         precision, nza, ndgs, quiet );
                    } catch (std::runtime_error e) {
                        //ostringstream os;
                        //os << "Calculation of SingleScatteringData properties failed for\n"
                        os << "f_grid[" << f_index << "] = " << ssd.f_grid[f_index] << "\n"
                        << "T_grid[" << T_index << "] = " << ssd.T_grid[T_index] << "\n\n";
                        //<< e.what();
                        //throw std::runtime_error(os.str());
                        thisfailed = true;
                        anyfailed = true;
                        cout << "\n\n";
                        
                    }

                    if( !thisfailed )
                    {
                    mono_pha_mat_data(joker, 0) = f11;
                    mono_pha_mat_data(joker, 1) = f12;
                    mono_pha_mat_data(joker, 2) = f22;
                    mono_pha_mat_data(joker, 3) = f33;
                    mono_pha_mat_data(joker, 4) = f34;
                    mono_pha_mat_data(joker, 5) = f44;

                    mono_pha_mat_data *= csca/4./PI;
                    ssd.pha_mat_data(f_index, T_index, joker, 0, 0, 0, joker) = mono_pha_mat_data;

                    ssd.ext_mat_data(f_index, T_index, 0, 0, 0) = cext;
                    ssd.abs_vec_data(f_index, T_index, 0, 0, 0) = cext - csca;
                    }
                }
            if( anyfailed )
              if( robust )
                cout << os.str();
              else
                throw std::runtime_error(os.str());
            else
              {
                os << "None\n";
              }
            
            break;
        }
        case PTYPE_AZIMUTH_RND:
        {
            const Index nza = ssd.za_grid.nelem();
            const Index naa = ssd.aa_grid.nelem();

            Index nza_inc = (nza-1)/2 + 1;

            ssd.pha_mat_data.resize(nf, nT, nza, naa, nza_inc, 1, 16);
            ssd.ext_mat_data.resize(nf, nT, nza_inc, 1, 3);
            ssd.abs_vec_data.resize(nf, nT, nza_inc, 1, 2);

            ssd.ext_mat_data = NAN;
            ssd.pha_mat_data = NAN;
            ssd.abs_vec_data = NAN;

            // Output variables
            Numeric cext = NAN;
            Numeric csca = NAN;
            Index nmax = -1;

            Tensor5 csca_data(nf, nT, nza_inc, 1, 2);

#pragma omp critical(tmatrix_ssp)
            for (Index f_index = 0; f_index < nf; ++f_index)
            {
                const Numeric lam_f = lam[f_index];

                for (Index T_index = 0; T_index < nT; ++T_index)
                {
                    try {
                        tmatrix_fixed_orientation
                        (cext, csca, nmax,
                         equiv_radius, aspect_ratio, np, lam_f,
                         ref_index_real(f_index, T_index),
                         ref_index_imag(f_index, T_index),
                         precision);
                    } catch (std::runtime_error e) {
                        ostringstream os;
                        os << "Calculation of SingleScatteringData properties failed for\n"
                        << "f_grid[" << f_index << "] = " << ssd.f_grid[f_index] << "\n"
                        << "T_grid[" << T_index << "] = " << ssd.T_grid[T_index] << "\n"
                        << e.what();
                        throw std::runtime_error(os.str());
                    }

                    
                    Matrix phamat;
                    for (Index za_scat_index = 0; za_scat_index < nza; ++za_scat_index)
                        for (Index aa_index = 0; aa_index < naa; ++aa_index)
                            for (Index za_inc_index = 0; za_inc_index < nza_inc; ++za_inc_index)
                            {
                                if (aspect_ratio < 1.0)
                                {
                                    // Phase matrix for prolate particles
                                    integrate_phamat_alpha10(phamat,
                                                             nmax, lam_f,
                                                             ssd.za_grid[za_inc_index],
                                                             ssd.za_grid[za_scat_index],
                                                             0.0,
                                                             ssd.aa_grid[aa_index],
                                                             90.0, 0.0, 180.0);
                                    phamat /= 180.;
                                }
                                else
                                {
                                    // Phase matrix for oblate particles
                                    calc_phamat(phamat,
                                                nmax, lam_f,
                                                ssd.za_grid[za_inc_index],
                                                ssd.za_grid[za_scat_index],
                                                0.0,
                                                ssd.aa_grid[aa_index],
                                                0.0, 0.0);
                                }

                                ssd.pha_mat_data(f_index, T_index,
                                                 za_scat_index, aa_index, za_inc_index,
                                                 0, Range(0, 4)) = phamat(0, joker);
                                ssd.pha_mat_data(f_index, T_index,
                                                 za_scat_index, aa_index, za_inc_index,
                                                 0, Range(4, 4)) = phamat(1, joker);
                                ssd.pha_mat_data(f_index, T_index,
                                                 za_scat_index, aa_index, za_inc_index,
                                                 0, Range(8, 4)) = phamat(2, joker);
                                ssd.pha_mat_data(f_index, T_index,
                                                 za_scat_index, aa_index, za_inc_index,
                                                 0, Range(12, 4)) = phamat(3, joker);
                            }

                    // Csca integral
                    for (Index za_scat_index = 0; za_scat_index < nza_inc; ++za_scat_index)
                    {
                        Matrix csca_integral;
                        if (aspect_ratio < 1.0)
                        {
                            // Csca for prolate particles
                            integrate_phamat_theta0_phi_alpha6(csca_integral, nmax, lam_f,
                                                               0, 180,
                                                               ssd.za_grid[za_scat_index],
                                                               0.,
                                                               0., 180.,
                                                               90.,
                                                               0., 180.);
                            csca_integral /= 180.;
                        }
                        else
                        {
                            // Csca for oblate particles
                            integrate_phamat_theta0_phi10(csca_integral, nmax, lam_f,
                                                          0, 180,
                                                          ssd.za_grid[za_scat_index],
                                                          0.,
                                                          0., 180,
                                                          0., 0.);
                        }
                        csca_data(f_index, T_index, za_scat_index, 0, joker) = csca_integral(Range(0, 2), 0);
                    }

                    // Extinction matrix
                    if (aspect_ratio < 1.0)
                    {
                        // Average T-Matrix for prolate particles
                        avgtmatrix_(nmax);
                    }

                    for (Index za_inc_index = 0; za_inc_index < nza_inc; ++za_inc_index)
                    {
                        Complex s11;
                        Complex s12;
                        Complex s21;
                        Complex s22;
                        VectorView K = ssd.ext_mat_data(f_index, T_index, za_inc_index, 0, joker);

                        const Numeric beta = 0.;
                        const Numeric alpha = 0.;
                        ampl_(nmax, lam_f,
                              ssd.za_grid[za_inc_index],
                              ssd.za_grid[za_inc_index],
                              0., 0.,
                              alpha, beta,
                              s11, s12, s21, s22);


                        K[0] = (Complex(0., -1.) * (s11+s22)).real();
                        K[1] = (Complex(0.,  1.) * (s22-s11)).real();
                        K[2] =                     (s22-s11).real();

                        K *= lam_f;

                    }
                }
            }

            csca_data *= 2.*PI*PI/32400.;
            ssd.abs_vec_data = ssd.ext_mat_data(joker, joker, joker, joker, Range(0, 2));
            ssd.abs_vec_data -= csca_data;
            
            break;
        }
        default:
        {
            std::ostringstream os;
            os << "Only particle types 20 and 30 are currently supported: "
            << ssd.ptype;
            throw std::runtime_error(os.str());
            break;
        }
    }
}


// Documentation in header file.
void tmatrix_ampld_test(const Verbosity& verbosity)
{
    CREATE_OUT0;

    out0 << "======================================================\n";
    out0 << "Test for nonspherical particles in a fixed orientation\n";
    out0 << "Output should match 3rdparty/tmatrix/tmatrix_ampld.ref\n";
    out0 << "======================================================\n";

    // Same inputs as in example included in original ampld.lp.f
    Numeric rat = 1.;
    Numeric axi = 10.; //[um]
    Index   np = -1;
    Numeric lam = acos(-1.)*2.; //[um]
    Numeric eps = 0.5;
    Numeric mrr = 1.5;
    Numeric mri = 0.02;
    Numeric ddelt = 0.001;

    Index   quiet = 1;

    // Output variables
    Index   nmax;
    Numeric csca;
    Numeric cext;
    char errmsg[1024] = "";

    tmatrix_(rat, axi, np, lam, eps, mrr, mri, ddelt, quiet,
             nmax, csca, cext, errmsg);

    out0 << "nmax: " << nmax << "\n";
    out0 << "csca: " << csca << " um2\n";
    out0 << "cext: " << cext << " um2\n";

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

    out0 << "AMPLITUDE MATRIX (all in [um]): \n";
    out0 << "s11: " << s11 << "\n";
    out0 << "s12: " << s12 << "\n";
    out0 << "s21: " << s21 << "\n";
    out0 << "s22: " << s22 << "\n";

    Matrix z;
    ampmat_to_phamat(z, s11, s12, s21, s22);
    //z *= 1e12; // meter^2 to micron^2 for comparison with original results

    out0 << "PHASE MATRIX (all un [um2]): \n" << z << "\n";
}


// Documentation in header file.
void tmatrix_tmd_test(const Verbosity& verbosity)
{
    CREATE_OUT0;

    out0 << "======================================================\n";
    out0 << "Test for randomly oriented nonspherical particles\n";
    out0 << "Output should match 3rdparty/tmatrix/tmatrix_tmd.ref\n";
    out0 << "======================================================\n";

    // Same inputs as in example included in original tmd.lp.f
    Numeric rat = 0.5;
    Index   ndistr = 3;
    Numeric axmax = 1.; //[um]
    Index   npnax = 2;
    Numeric b = 0.1;
    Numeric gam = 0.5;
    Index   nkmax = 5;
    Numeric eps = 2;
    Index   np = -1;
    Numeric lam = 0.5; //[um]
    Numeric mrr = 1.53;
    Numeric mri = 0.008;
    Numeric ddelt = 0.001;
    Index   npna = 19;
    Index   ndgs = 2;
    Numeric r1rat = 0.89031;
    Numeric r2rat = 1.56538;

    Index   quiet = 1;

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

    out0 << "reff: " << reff << " um\n";
    out0 << "veff: " << veff << "\n";
    out0 << "cext: " << cext << " um2\n";
    out0 << "csca: " << csca << " um2\n";
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


// Documentation in header file.
void calc_ssp_random_test(const Verbosity& verbosity)
{
    CREATE_OUT0;
    out0 << "======================================================\n";
    out0 << "Test calculation of single scattering data\n";
    out0 << "for randomly oriented, oblate particles\n";
    out0 << "======================================================\n";


    SingleScatteringData ssd;

    ssd.ptype = PTYPE_TOTAL_RND;
    ssd.f_grid = MakeVector(230e9, 240e9);
    ssd.T_grid = MakeVector(220, 250);
    nlinspace(ssd.za_grid, 0, 180, 19);
    nlinspace(ssd.aa_grid, 0, 180, 19);

    // Refractive index real and imagenary parts
    // Dimensions: [nf, nT];
    Matrix mrr(ssd.f_grid.nelem(), ssd.T_grid.nelem(), 1.78031135);
    Matrix mri(ssd.f_grid.nelem(), ssd.T_grid.nelem(), 0.00278706);

    mrr(0, 0) = 1.78031135; mrr(0, 1) = 1.78150475;
    mrr(1, 0) = 1.78037238; mrr(1, 1) = 1.78147686;

    mri(0, 0) = 0.00278706; mri(0, 1) = 0.00507565;
    mri(1, 0) = 0.00287245; mri(1, 1) = 0.00523012;

    calcSingleScatteringDataProperties(ssd, mrr, mri,
                                       200.e-6,
                                       -1,
                                       1.5);

    out0 << "ssd.pha_mat_data(0, 0, joker, 0, 0, joker, joker):\n"
    << ssd.pha_mat_data(0, 0, joker, 0, 0, joker, joker) << "\n\n";

    out0 << "ssd.ext_mat_data:\n" << ssd.ext_mat_data << "\n\n";
    out0 << "ssd.abs_vec_data:\n" << ssd.abs_vec_data << "\n\n";

    out0 << "======================================================\n";
    out0 << "Test calculation of single scattering data\n";
    out0 << "for randomly oriented, prolate particles\n";
    out0 << "======================================================\n";

    calcSingleScatteringDataProperties(ssd, mrr, mri,
                                       200.e-6,
                                       -1,
                                       0.7);

    out0 << "ssd.pha_mat_data(0, 0, joker, 0, 0, joker, joker):\n"
    << ssd.pha_mat_data(0, 0, joker, 0, 0, joker, joker) << "\n\n";

    out0 << "ssd.ext_mat_data:\n" << ssd.ext_mat_data << "\n\n";
    out0 << "ssd.abs_vec_data:\n" << ssd.abs_vec_data << "\n\n";
}


// Documentation in header file.
void calc_ssp_fixed_test(const Verbosity& verbosity)
{
    CREATE_OUT0;
    out0 << "======================================================\n";
    out0 << "Test calculation of single scattering data\n";
    out0 << "for oblate particles with fixed orientation\n";
    out0 << "======================================================\n";


    SingleScatteringData ssd;

    ssd.ptype = PTYPE_AZIMUTH_RND;
    ssd.f_grid = MakeVector(230e9, 240e9);
    ssd.T_grid = MakeVector(220, 250);
    nlinspace(ssd.za_grid, 0, 180, 19);
    nlinspace(ssd.aa_grid, 0, 180, 19);

    // Refractive index real and imagenary parts
    // Dimensions: [nf, nT];
    Matrix mrr(ssd.f_grid.nelem(), ssd.T_grid.nelem(), 1.78031135);
    Matrix mri(ssd.f_grid.nelem(), ssd.T_grid.nelem(), 0.00278706);

    mrr(0, 0) = 1.78031135; mrr(0, 1) = 1.78150475;
    mrr(1, 0) = 1.78037238; mrr(1, 1) = 1.78147686;

    mri(0, 0) = 0.00278706; mri(0, 1) = 0.00507565;
    mri(1, 0) = 0.00287245; mri(1, 1) = 0.00523012;

    calcSingleScatteringDataProperties(ssd, mrr, mri,
                                       200.e-6,
                                       -1,
                                       1.5);

    out0 << "ssd.pha_mat_data(0, 0, joker, 0, 0, joker, joker):\n"
    << ssd.pha_mat_data(0, 0, joker, 0, 0, joker, joker) << "\n\n";

    out0 << "ssd.ext_mat_data(0, 0, joker, joker, joker):\n"
    << ssd.ext_mat_data(0, 0, joker, joker, joker) << "\n\n";
    
    out0 << "abs_vec_data:\n" << ssd.abs_vec_data << "\n\n";

    out0 << "======================================================\n";
    out0 << "Test calculation of single scattering data\n";
    out0 << "for prolate particles with fixed orientation\n";
    out0 << "======================================================\n";
    calcSingleScatteringDataProperties(ssd, mrr, mri,
                                       200.e-6,
                                       -1,
                                       0.7);
    
    out0 << "ssd.pha_mat_data(0, 0, joker, 0, 0, joker, joker):\n"
    << ssd.pha_mat_data(0, 0, joker, 0, 0, joker, joker) << "\n\n";

    out0 << "ssd.ext_mat_data(0, 0, joker, joker, joker):\n"
    << ssd.ext_mat_data(0, 0, joker, joker, joker) << "\n\n";

    out0 << "abs_vec_data:\n" << ssd.abs_vec_data << "\n\n";    
}

