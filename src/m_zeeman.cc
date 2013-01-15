/* Copyright (C) 2012
Richard Larsson <ric.larsson@gmail.com>

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

#include <cmath>
#include <stdexcept>
#include <string>
#include "auto_md.h"
#include "arts.h"
#include "ppath.h"
#include "messages.h"
#include "math_funcs.h"
#include "absorption.h"
#include "abs_species_tags.h"
#include "physics_funcs.h"
#include "matpackIII.h"
#include "rte.h"
#include "rational.h"

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PLANCK_CONST;
extern const Numeric BOHR_MAGNETON;
extern const Numeric LANDE_GS;


/* It is very common that people use Lenoir notations.
 * It may be better to use some direct equation.
 * This is an easy way to switch between the schemes. */
#define FREQUENCY_CHANGE FrequencyChangeDirectWithGS // Direct is obviously best if you understand Lande.
#define EFFECTIVE_LANDE_G Lande_gs // Lande is not well understood though. gs agrees most with other models. gs seems best for now.
#define K_MATRIX KMatrixAttenuationLenoir // Lenoir is best understood. The diffraction/refraction terms are poorly understood in Rees, though it seems to mirror things.
#define RELATIVE_STRENGTH RelativeStrengthSMMBerdyugina // SMMBerdyugina seems best and agrees most with Pardo O2


/*!
    Helper function for sort. Used to sort the LineRecord(s) in an
    ArrayOfLineRecord by increasing frequency. Usage example:

    sort(ArrayOfLineRecord_A.begin(), ArrayOfLineRecord_A.end(), sortF)
    returns ArrayOfLineRecord_A sorted by frequency

    \param  i   In:    LineRecord helper 1.
    \param  j   In:    LineRecord helper 2.

    \author Richard Larsson
    \date   2012-08-03
*/
bool sortF(LineRecord i, LineRecord j) { return ( i.F() < j.F() ); };


/*!
    Defines the rotation extinction matrix as derived by the author from
    the definitions of Lenoir (1967) and Mishchenko (2002).

    \param  K       Out:    The rotation extinction matrix.
    \param  theta   In:     Angle between the magnetic field and the
                            propagation path. In radians.
    \param  eta     In:     Angle to rotate planar polarization clockwise to
                            fit the general coordinate system. In radians.
    \param  DM      In:     Change in the secondary rotational quantum number.

    \author Richard Larsson
    \date   2012-08-03
*/
void KMatrixAttenuationLenoir(MatrixView K, const Numeric theta, const Numeric eta, const Index DM)
{
    assert(K.nrows() == 4 );
    assert(K.ncols() == 4 );
    
    switch( DM )
    {
        case -1:
            K(0,0) =   1 + cos(theta)*cos(theta);
            K(0,1) =       sin(theta)*sin(theta) * cos( 2 * eta );
            K(0,2) =       sin(theta)*sin(theta) * sin( 2 * eta );
            K(0,3) =   2 * cos(theta);

            K(1,0) = K(0,1);
            K(1,1) = K(0,0);
            K(1,2) = 0;
            K(1,3) = 0;

            K(2,0) = K(0,2);
            K(2,1) = K(1,2);
            K(2,2) = K(0,0);
            K(2,3) = 0;

            K(3,0) = K(0,3);
            K(3,1) = K(1,3);
            K(3,2) = K(2,3);
            K(3,3) = K(0,0);
            break;
        case  1:
            K(0,0) =   1 + cos(theta)*cos(theta);
            K(0,1) =       sin(theta)*sin(theta) * cos( 2 * eta );
            K(0,2) =       sin(theta)*sin(theta) * sin( 2 * eta );
            K(0,3) = - 2 * cos(theta);

            K(1,0) = K(0,1);
            K(1,1) = K(0,0);
            K(1,2) = 0;
            K(1,3) = 0;

            K(2,0) = K(0,2);
            K(2,1) = K(1,2);
            K(2,2) = K(0,0);
            K(2,3) = 0;

            K(3,0) = K(0,3);
            K(3,1) = K(1,3);
            K(3,2) = K(2,3);
            K(3,3) = K(0,0);
            break;
        case  0:
            K(0,0) =   sin(theta)*sin(theta);
            K(0,1) = - sin(theta)*sin(theta) * cos( 2 * eta );
            K(0,2) = - sin(theta)*sin(theta) * sin( 2 * eta );
            K(0,3) =   0;

            K(1,0) = K(0,1);
            K(1,1) = K(0,0);
            K(1,2) = 0;
            K(1,3) = 0;

            K(2,0) = K(0,2);
            K(2,1) = K(1,2);
            K(2,2) = K(0,0);
            K(2,3) = 0;

            K(3,0) = K(0,3);
            K(3,1) = K(1,3);
            K(3,2) = K(2,3);
            K(3,3) = K(0,0);
            break;
        default: // Unity matrix
            K(0,0) = 1;
            K(0,1) = 0;
            K(0,2) = 0;
            K(0,3) = 0;

            K(1,0) = K(0,1);
            K(1,1) = K(0,0);
            K(1,2) = 0;
            K(1,3) = 0;

            K(2,0) = K(0,2);
            K(2,1) = K(1,2);
            K(2,2) = K(0,0);
            K(2,3) = 0;

            K(3,0) = K(0,3);
            K(3,1) = K(1,3);
            K(3,2) = K(2,3);
            K(3,3) = K(0,0);
            break;
    };
};


/*!
 *  Defines the rotation extinction matrix as derived by the author from
 *  the definitions of Lenoir (1967) and Mishchenko (2002).
 * 
 *  \param  K       Out:    The rotation extinction matrix.
 *  \param  theta   In:     Angle between the magnetic field and the
 *                          propagation path. In radians.
 *  \param  eta     In:     Angle to rotate planar polarization clockwise to
 *                          fit the general coordinate system. In radians.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 * 
 *  \author Richard Larsson
 *  \date   2012-08-03
 */
void KMatrixRees(MatrixView K, const Numeric theta, const Numeric eta, const Index DM)
{
    assert(K.nrows() == 4 );
    assert(K.ncols() == 4 );
    
    switch( DM )
    {
        case -1:
            K(0,0) =   1 + cos(theta)*cos(theta);
            K(0,1) = - sin(theta)*sin(theta) * cos(2*eta);
            K(0,2) = - sin(theta)*sin(theta) * sin(2*eta);
            K(0,3) =   2 * cos(theta);
            
            K(1,0) =   K(0,1);
            K(1,1) =   K(0,0);
            K(1,2) =   2 * cos(theta);
            K(1,3) =   sin(theta)*sin(theta) * sin(2*eta);
            
            K(2,0) =   K(0,2);
            K(2,1) = - K(1,2);
            K(2,2) =   K(0,0);
            K(2,3) = - sin(theta)*sin(theta) * cos(2*eta);
            
            K(3,0) =   K(0,3);
            K(3,1) = - K(1,3);
            K(3,2) = - K(2,3);
            K(3,3) =   K(0,0);
            break;
        case  1:
            K(0,0) =   1 + cos(theta)*cos(theta);
            K(0,1) = - sin(theta)*sin(theta) * cos(2*eta);
            K(0,2) = - sin(theta)*sin(theta) * sin(2*eta);
            K(0,3) = - 2 * cos(theta);
            
            K(1,0) =   K(0,1);
            K(1,1) =   K(0,0);
            K(1,2) = - 2 * cos(theta);
            K(1,3) =   sin(theta)*sin(theta) * sin(2*eta);
            
            K(2,0) =   K(0,2);
            K(2,1) = - K(1,2);
            K(2,2) =   K(0,0);
            K(2,3) = - sin(theta)*sin(theta) * cos(2*eta);
            
            K(3,0) =   K(0,3);
            K(3,1) = - K(1,3);
            K(3,2) = - K(2,3);
            K(3,3) =   K(0,0);
            break;
        case  0:
            K(0,0) =   sin(theta)*sin(theta);
            K(0,1) =   sin(theta)*sin(theta) * cos(2*eta);
            K(0,2) =   sin(theta)*sin(theta) * sin(2*eta);
            K(0,3) =   0;
            
            K(1,0) =   K(0,1);
            K(1,1) =   K(0,0);
            K(1,2) =   0;
            K(1,3) = - sin(theta)*sin(theta) * sin(2*eta);
            
            K(2,0) =   K(0,2);
            K(2,1) = - K(1,2);
            K(2,2) =   K(0,0);
            K(2,3) =   sin(theta)*sin(theta) * cos(2*eta);
            
            K(3,0) =   K(0,3);
            K(3,1) = - K(1,3);
            K(3,2) = - K(2,3);
            K(3,3) =   K(0,0);
            break;
        default: // Unity matrix
            K(0,0) = 1;
            K(0,1) = 0;
            K(0,2) = 0;
            K(0,3) = 0;
            
            K(1,0) = K(0,1);
            K(1,1) = K(0,0);
            K(1,2) = 0;
            K(1,3) = 0;
            
            K(2,0) = K(0,2);
            K(2,1) = K(1,2);
            K(2,2) = K(0,0);
            K(2,3) = 0;
            
            K(3,0) = K(0,3);
            K(3,1) = K(1,3);
            K(3,2) = K(2,3);
            K(3,3) = K(0,0);
            break;
    };
};


Numeric Lande_gs(const Numeric N, const Numeric J, const Numeric S, const Numeric GS) { return (GS*(J*(J+1.)+S*(S+1.)-N*(N+1.))/(J*(J+1.))/2.0); }
Numeric Lande_RusselSaunders(const Numeric N, const Numeric J, const Numeric S, const Numeric) { return (1.  +    (J*(J+1.)+S*(S+1.)-N*(N+1.))/(J*(J+1.))/2.0); }


/*! Calculate the Wigner 3j function value
     From Wolfram-Alpha definition at http://mathworld.wolfram.com/Wigner3j-Symbol.html
   
     Also found in Landau and Lifshitz Volume 3 as equation (106.14). This source is more descriptive.
     
     This function may be incredible slow depending on input as there is an loop over several factorials involved.
     
     The matrix is:
     | jf  dj  ji |
     | mf  dm -mi |
     
      \param  jf
      \param  dj
      \param  ji
      \param  mf
      \param  dm
      \param  mi
     
      \author Richard Larsson
      \date   2012-10-12
 */
Numeric Wigner_3j_Racah_Formula(const Rational jf, const Index dj, const Rational ji, const Rational mf, const Index dm, const Rational mi)
{
    
    Numeric sum_object=0;
    Numeric a;
    Index t=0;
    
    while( ji - dj + t + mf >= 0 && ji - jf + t - dm >= 0 && jf + dj - ji -t>= 0 && jf - t - mf >= 0 && dj - t + dm >= 0)
    {
        a = t%2==0?1.:-1.;
        
        sum_object += a / ( fac(t) * fac(ji - dj + t + mf) * 
        fac(ji - jf + t - dm) * fac(jf + dj - ji -t) * fac(jf - t - mf) * fac(dj - t + dm) );
        
        t++;
    }
    
    a = (jf-dj-mi)%2==0?1.:-1.;

    return sum_object * a *    
    sqrt(fac(jf+dj-ji) * fac(jf-dj+ji) * fac(-jf+dj+ji) / fac(jf+dj+ji+1)) *
    sqrt(fac(jf+mf) * fac(jf-mf) * fac(dj+dm) * fac(dj-dm) * fac(ji+mi) * fac(ji-mi));
}


/*!
    Return the relative strength of the split Zeeman line parts as found in
    Landau and Lifshitz QED Volume 4 Equation (51.8).
    
    Renormalized (somewhat) to one instead of two.  Note that this 
    renormalization may cause some problem in the future because it
    is done within the equation, which is quite complex. See the description
    of Wigner 3j for more information and mathematics.
    
    Goal:
    PI consist of half the probability and SIGMA the other half.
    
    Only tested for Lenoir cases.
    
    \param  __U__   Void.
    \param  m       In:     Secondary rotational quantum number.
    \param  j       In:     Spin-Orbit Coupling number.
    \param  DJ      In:     Change in the main rotational quantum number.
    \param  DM      In:     Change in the secondary rotational quantum number.
    
    \author Richard Larsson
    \date   2012-10-15
 */
Numeric RelativeStrengthLandauAndLifshitz(Rational, Rational m, Rational j, Index dj, Index dm)
{
    Numeric ret_val = 0;
    
    if( dm == 1 )
    {
        if( dj!=0 )
            ret_val = Wigner_3j_Racah_Formula(j, abs(dj), j-dj, m+dm, -dm, -m) ;
        else
            ret_val = Wigner_3j_Racah_Formula(j, abs(dj), j-dj, -(m+dm), -(-dm), -(-m)) ;
    }
    else if( dm == -1 )
    {
        if(dj!=0)
            ret_val = Wigner_3j_Racah_Formula(j, abs(dj), j-dj, -(m+dm), -(-dm), -(-m)) ;
        else
            ret_val = Wigner_3j_Racah_Formula(j, 0, j, m+dm, -dm, -m) ;
    }
    else if( dm == 0)
    {
        if( dj < 0 )
            ret_val = Wigner_3j_Racah_Formula(j, abs(dj), j-dj, -m, 0 , m);
        else if( dj > 0 )
            ret_val = Wigner_3j_Racah_Formula(j-dj, abs(dj), j, -m, 0 , m);
        else if( dj == 0 && dm == 0 && m == 0 )
            ret_val = 0;
        else
            ret_val = Wigner_3j_Racah_Formula(j, 0, j, m, 0 , -m);
    }
    else
    {
        throw runtime_error("This cannot happen. Please notify someone or fix the problem if you see this error in Wigner_3j_Racah_Formula.");
    }
    return (dm==0 ? 1.5: 0.75) * ret_val * ret_val; // Again, div by 2 so to keep sum approx 1.
}


/*!
   Return the relative strength of the split Zeeman line parts as found in
   Berdyugina and Solnaki (2002). Note that this is the same as Hund
   case (a)-(a) transition as in Schadee (1978).
   
   Goal:
   PI consist of half the probability and SIGMA the other half.*
   
   Only tested for Lenoir cases.
   
   \param  __U__   Void.
   \param  m       In:     Secondary rotational quantum number.
   \param  j       In:     Spin-Orbit Coupling number.
   \param  DJ      In:     Change in the main rotational quantum number.
   \param  DM      In:     Change in the secondary rotational quantum number.
   
   \author Richard Larsson
   \date   2012-10-26
 */
Numeric RelativeStrengthSMMBerdyugina(Rational, Rational m, Rational j, Index dj, Index dm)
{
    const Numeric J = (Numeric)(j-dj), M = (Numeric)(m);
    
    Numeric ret_val;
    
    switch ( dj )
    {
        case -1:
            switch ( dm )
            {
                case -1:
                    ret_val = (0.75)*(J+M)*(J-1+M)/(2.*J*(2.*J-1)*(2.*J+1));
                    break;
                case  0:
                    ret_val = (1.50)*(J*J-M*M)/(J*(2.*J-1.)*(2.*J+1.));
                    break;
                case +1:
                    ret_val = (0.75)*(J-M)*(J-1.-M)/(2.*J*(2.*J-1.)*(2.*J+1.));
                    break;
                default:
                    throw runtime_error("Something is extremely wrong in RelativeStrengthSMMBerdyugina or the use thereof.");
                    break;
                
            }
            break;
        case  0:
            switch ( dm )
            {
                case -1:
                    ret_val = (0.75)*(J+M)*(J+1.-M)/(2.*J*(J+1.)*(2.*J+1.));
                    break;
                case  0:
                    ret_val = (1.50)*M*M/(J*(J+1.)*(2.*J+1.));
                    break;
                case +1:
                    ret_val = (0.75)*(J-M)*(J+1.+M)/(2.*J*(J+1.)*(2.*J+1.));
                    break;
                default:
                    throw runtime_error("Something is extremely wrong in RelativeStrengthSMMBerdyugina or the use thereof.");
                    break;
            }
            break;
        case +1:
            switch ( dm )
            {
                case -1:
                    ret_val = (0.75)*(J+1.-M)*(J+2.-M)/(2.*(J+1.)*(2.*J+1.)*(2.*J+3.));
                    break;
                case  0:
                    ret_val = (1.5)*((J+1.)*(J+1.)-M*M)/((J+1.)*(2.*J+1.)*(2.*J+3.));
                    break;
                case +1:
                    ret_val = (0.75)*(J+1.+M)*(J+2.+M)/(2.*(J+1.)*(2.*J+1.)*(2.*J+3.));
                    break;
                default:
                    throw runtime_error("Something is extremely wrong in RelativeStrengthSMMBerdyugina or the use thereof.");
                    break;
            }
            break;
        default:
            throw runtime_error("Something is extremely wrong in RelativeStrengthSMMBerdyugina or the use thereof.");
            break;
    }
    
    return ret_val;
}


/*!
    Return the relative strength of the split Zeeman line parts as found in
    Table 2 of Liebe and Hufford (1989). Renormalized to one instead of two.
   
    \param  n       In:     Main quantum number.
    \param  m       In:     Secondary rotational quantum number.
    \param  __U__   Void.
    \param  DJ      In:     Change in the main rotational quantum number.
    \param  DM      In:     Change in the secondary rotational quantum number.
   
    \author Richard Larsson
    \date   2012-08-03
 */
Numeric RelativeStrengthLenoir(Rational n, Rational m, Rational , Index DJ, Index DM)
{
    
    Numeric N = (Numeric)n, M = (Numeric)m, relstr;
    
    switch ( DJ )
    {//* Switch over DJ, if DJ != 1 or -1, DJ is wrong in Lenoir.
        case 1:
            switch ( DM )
            {//* Switch over DM, if DM != 1, 0 or -1, DM is wrong.
                case  1:
                    relstr = 3.0 * ( N+M+1 ) * ( N+M+2 ) / ( 4.0 * ( N+1 ) * ( 2*N+1 ) * ( 2*N+3 ) ) ;
                    break;
                case  0:
                    relstr = 3.0 * ( (N+1) * (N+1) - M * M) / ( ( N+1 ) * ( 2*N+1 ) * ( 2*N+3 ) );
                    break;
                case -1:
                    relstr = 3.0 * ( N-M+1 ) * ( N-M+2 ) / ( 4.0 * ( N+1 ) * ( 2*N+1 ) * ( 2*N+3 ) ) ;
                    break;
                default:
                    throw runtime_error("Something is extremely wrong in RelativeStrengthLenoir or the use thereof.");
                    break;
            }
            break;
            case -1:
                switch ( DM )
                {//* Switch over DM, if DM != 1, 0 or -1, DM is wrong.
                    case  1:
                        relstr = 3.0 * ( N-M-1 ) * ( N-M ) / ( 4.0 * N * ( 2*N+1 ) * ( 2*N-1 ) );
                        break;
                    case  0:
                        relstr = 3.0 * ( N * N - M * M ) / ( N * ( 2*N+1 ) * ( 2*N-1 ) );
                        break;
                    case -1:
                        relstr = 3.0 * ( N+M-1 ) * ( N+M ) / ( 4.0 * N * ( 2*N+1 ) * ( 2*N-1 ) );
                        break;
                    default:
                        throw runtime_error("Something is extremely wrong in RelativeStrengthLenoir or the use thereof.");
                        break;
                }
                break;
        default:
            throw runtime_error("Something is extremely wrong in RelativeStrengthLenoir or the use thereof. " 
            "Are you sure you checked that there are no DJ = 0 particles in the Zeeman mix?");
            break;
    };
return relstr / 2.0; //FIXME? It seems the normalization used by Hufford and Liebe is to 2 and not 1... this factor might be erroneous?
}


/*!
 *  Return the frequency change of the split Zeeman line parts as found in
 *  Table 2 of Liebe and Hufford (1989).
 *
 *  \param  n       In:     Main quantum number.
 *  \param  m       In:     Secondary rotational quantum number.
 *  \param  DJ      In:     Change in the coupled rotational quantum number.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 *  \param  H_mag   In:     Magnitude of the magnetic field in Tesla.
 *
 *  \author Richard Larsson
 *  \date   2012-08-03
*/
Numeric FrequencyChangeLenoir(Rational n, Rational m, Rational, Numeric, Index DJ, Index DM, Index, Numeric H_mag, Numeric)
{
    /*
    The following switch case is from table 2 of Lenoir, 1968.
    */

    const Numeric N = (Numeric)n, M = (Numeric)m;
    Numeric fcc;
    const Numeric KH = 2.8026 * 1000000 * H_mag * 10000; // KH from Lenoir 1968. Factors are for conversion from Lenoir's' to SI units.

    switch ( DJ )
    {//* Switch over DJ, if DJ != 1 or -1, DJ is wrong.
        case 1:
            switch ( DM )
            {//* Switch over DM, if DM != 1, 0 or -1, DM is wrong.
                case 1:
                    fcc = 1 / ( N + 1 ) * ( 1 + M * ( N - 1 ) / N );
                    break;
                case  0:
                    fcc = 1 / ( N + 1 ) * M * ( N - 1 ) / N;
                    break;
                case  -1:
                    fcc = 1 / ( N + 1 ) * ( M * ( N - 1 ) / N - 1 );
                    break;
                default:
                    fcc = 0;
                    break;
            };
            break;
        case -1:
            switch ( DM )
            {//* Switch over DM, if DM != 1, 0 or -1, DM is wrong.
                case 1:
                    fcc = - 1 / N * ( 1 + M * ( N +2 ) / ( N + 1 ) );
                    break;
                case  0:
                    fcc = - 1 / N * M * ( N +2 ) / ( N + 1 );
                    break;
                case -1:
                    fcc = - 1 / N * ( M * ( N +2 ) / ( N + 1 ) - 1 );
                    break;
                default:
                    fcc = 0;
                    break;
            };
            break;
        default:
            fcc = 0;
            break;
    };
    return KH*fcc;
}


/*!
 *  Return the frequency change of the split Zeeman line parts as found from
 *  g_s x M - g_s' x M'
 * 
 *  \param  n       In:     Main quantum number.
 *  \param  m       In:     Secondary rotational quantum number.
 *  \param  j       In:     Coupled rotational quantum number.
 *  \param  s       In:     Electron total spin quantum number.
 *  \param  DJ      In:     Change in the coupled rotational quantum number.
 *  \param  DM      In:     Change in the secondary rotational quantum number.
 *  \param  DN      In:     Change in the main rotational quantum number.
 *  \param  H_mag   In:     Magnitude of the magnetic field in Tesla.
 * 
 *  \author Richard Larsson
 *  \date   2012-09-25
 */
Numeric FrequencyChangeDirect(Rational n, Rational m, Rational j, Numeric S, Index DJ, Index DM, Index DN, Numeric H_mag, Numeric)
{
    const Numeric       N_up = (Numeric)n;
    const Numeric       N_lo = N_up - (Numeric)DN;
    const Numeric       M_lo = (Numeric)m;
    const Numeric       M_up = M_lo + (Numeric)DM; // This will probably confuse people even though it is correct. Should fix so outer loop is over M_up.
    const Numeric       J_up = (Numeric)j;
    const Numeric       J_lo = J_up - (Numeric)DJ;
    
    Numeric fcc;
    
    //NEEDED: Special rule for which g-value to use.
    
    // special case is for fitting with Lenoir. 
    fcc = (!(j == 0 && DJ == -1)) ? 
        - H_mag * (M_lo) * EFFECTIVE_LANDE_G(N_lo,J_lo,S,LANDE_GS) / PLANCK_CONST * BOHR_MAGNETON + 
          H_mag * (M_up) * EFFECTIVE_LANDE_G(N_up,J_up,S,LANDE_GS) / PLANCK_CONST * BOHR_MAGNETON : 
          H_mag * (M_lo) * EFFECTIVE_LANDE_G(N_lo,J_lo,S,LANDE_GS) / PLANCK_CONST * BOHR_MAGNETON ;
          
    return fcc;
}

/*!
    Return the frequency change of the split Zeeman line parts as found from
    g_s x M - g_s' x M'
    
    This version takes into account non-free electron GS constant
   
    \param  n       In:     Main quantum number.
    \param  m       In:     Secondary rotational quantum number.
    \param  j       In:     Coupled rotational quantum number.
    \param  s       In:     Electron total spin quantum number.
    \param  DJ      In:     Change in the coupled rotational quantum number.
    \param  DM      In:     Change in the secondary rotational quantum number.
    \param  DN      In:     Change in the main rotational quantum number.
    \param  H_mag   In:     Magnitude of the magnetic field in Tesla.
    \param  GS      In:     G-consant from data for the moleculer.
   
    \author Richard Larsson
    \date   2012-11-13
 */
Numeric FrequencyChangeDirectWithGS(Rational n, Rational m, Rational j, Numeric S, Index DJ, Index DM, Index DN, Numeric H_mag, Numeric GS)
{
    const Numeric       N_up = (Numeric)n;
    const Numeric       N_lo = N_up - (Numeric)DN;
    const Numeric       M_lo = (Numeric)m;
    const Numeric       M_up = M_lo + (Numeric)DM; // This will probably confuse people even though it is correct. Should fix so outer loop is over M_up.
    const Numeric       J_up = (Numeric)j;
    const Numeric       J_lo = J_up - (Numeric)DJ;
    
    Numeric fcc;
    
    //NEEDED: Special rule for which g-value to use.
    
    // special case is for fitting with Lenoir. 
    fcc = (!(j == 0 && DJ == -1)) ? 
    - H_mag * (M_lo) * EFFECTIVE_LANDE_G(N_lo,J_lo,S,GS) / PLANCK_CONST * BOHR_MAGNETON + 
    H_mag * (M_up)   * EFFECTIVE_LANDE_G(N_up,J_up,S,GS) / PLANCK_CONST * BOHR_MAGNETON : 
    H_mag * (M_lo)   * EFFECTIVE_LANDE_G(N_lo,J_lo,S,GS) / PLANCK_CONST * BOHR_MAGNETON ;
    
    return fcc;
}

/*!
Helper function. This is the only place where m_zeeman interacts with other absorption protocols.
*/
void Part_Return_Zeeman(  Tensor3View part_abs_mat, const ArrayOfArrayOfSpeciesTag& abs_species, const ArrayOfLineshapeSpec& abs_lineshape, 
                          const ArrayOfLineRecord& lr, const SpeciesAuxData& isotopologue_ratios, const Matrix& abs_vmrs, const Vector& abs_p,
                          const Vector& abs_t, const Vector& f_grid, const Numeric& theta, const Numeric& eta, const Index& DM, const Index& this_species,
                          const Verbosity& verbosity )
{
    assert( part_abs_mat.npages() == f_grid.nelem() && part_abs_mat.ncols() == 4 && part_abs_mat.nrows() == 4 );
    
    Matrix temp1(f_grid.nelem(), 1, 0.);
    
    for ( Index i=0; i<abs_species[this_species].nelem(); ++i )
    {
        Matrix temp2(f_grid.nelem(), 1, 0.);
        
        xsec_species( temp2, f_grid, abs_p, abs_t, abs_vmrs, abs_species, this_species, lr,
                    abs_lineshape[i].Ind_ls(), abs_lineshape[i].Ind_lsn(), abs_lineshape[i].Cutoff(),
                    isotopologue_ratios, verbosity ); // Now in cross section
        
        temp2 *= abs_vmrs(this_species, 0) * number_density( abs_p[0],abs_t[0]); // Now in absorption coef.
        
        temp1 += temp2;
    }
    Matrix  K(4,4); K_MATRIX(K, theta*DEG2RAD, eta*DEG2RAD, DM);
    mult(part_abs_mat,temp1(joker,0), K);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_mat_per_speciesAddZeeman(Tensor4& abs_mat_per_species,
				      const Vector& f_grid,
				      const ArrayOfArrayOfSpeciesTag& abs_species,
				      const ArrayOfArrayOfLineRecord& abs_lines_per_species,
				      const ArrayOfLineshapeSpec& abs_lineshape,
                                     const SpeciesAuxData& isotopologue_ratios,
                                     const SpeciesAuxData& isotopologue_quantum,
				      const Numeric& rte_pressure,
				      const Numeric& rte_temperature,
				      const Vector& rte_vmr_list,
				      const Numeric& rte_doppler,
				      const Vector& rte_mag,
				      const Vector& ppath_los,
				      const Index& atmosphere_dim,
				      const Verbosity& verbosity)
{
    CREATE_OUT3;
    Vector R_path_los;
    mirror_los(R_path_los, ppath_los, atmosphere_dim);
    
    bool do_zeeman = false;

    /*
        This function will, for each Zeeman species, make a local
        ArrayOfLineRecord(s) for the various transition types with Zeeman
        altered LineRecord(s).  These are then composed into a single
        ArrayOfArrayOfLineRecord which is processed as per the scalar case.
        
        The line broadened absorption coefficients are finally multiplied with
        the transition type rotation matrix and the new variable is returned.
        
        Note that between 55 GHz and 65 GHz there is usually ~700 O_2 lines,
        however, when this Zeeman splitting method is used, the number of
        lines is increased to about 45,000. This is a time consuming method...
    */

    // Check that correct isotopologue ratios are defined for the species
    // we want to calculate
    checkIsotopologueRatios(abs_species, isotopologue_ratios);

    // Begin TEST(s)
    for(Index II = 0; II<abs_species.nelem(); II++)
        if(is_zeeman(abs_species[II])) { do_zeeman = true; break; } // If any species is Zeeman, do it.
    if( abs_mat_per_species.ncols()  != 4 )
        throw runtime_error("Zeeman Effect is only implemented for Stokes dimension 4.");
    if( abs_mat_per_species.nrows()  != 4 )
        throw runtime_error("Zeeman Effect is only implemented for Stokes dimension 4.");
    if( abs_mat_per_species.npages() != f_grid.nelem() )
        throw runtime_error("Frequency dimension of abs_mat_per_species not equal to length of f_grid.");
    if( abs_mat_per_species.nbooks() != abs_species.nelem() )
        throw runtime_error("Species dimension of abs_mat_per_species not equal to length of abs_species.");
    if( rte_mag.nelem() != 3 )
      throw runtime_error("*rte_mag* must have length 3.");
    // End   TEST(s)
    Vector local_f_grid;
    // Make pointer point to original.
    const Vector* f_grid_pointer = &f_grid;

    /*
        Doppler treatment, do this only if there is a non-zero Doppler
        shift. We do this after the frequency selection, so in the case
        that we have only a single frequency, we have to shift only that!

        Unfortunately, we need yet another local copy of f_grid. In
        contrast to the frequency selection, we here want to modify the
        actual frequency values inside!
    */
    Vector local_doppler_f_grid;
    if (rte_doppler==0)
    {
        out3 << "  Doppler shift: None\n";
    }
    else
    {
        ostringstream os;
        os << "  Doppler shift: " << rte_doppler << " Hz\n";
        out3 << os.str();

        Numeric local_doppler;
        NumericScale( local_doppler, rte_doppler, -1, verbosity );
        // I could just have multiplied by -1 directly, but I like using
        // the WSM here.

        VectorAddScalar( local_doppler_f_grid,  *f_grid_pointer, local_doppler, verbosity );

        // Make pointer point to the doppler shifted frequency grid.
        f_grid_pointer = &local_doppler_f_grid;
    }

    // Using the standard scalar absorption functions to get physics parameters,
    Vector abs_p, abs_t; Matrix abs_vmrs;
    AbsInputFromRteScalars( abs_p, abs_t, abs_vmrs,                        // Output
                            rte_pressure, rte_temperature, rte_vmr_list,  //Input
                            verbosity);                                  // Verbose!
    if( do_zeeman  && ( rte_mag[0]!=0 || rte_mag[1]!=0 || rte_mag[2]!=0 ) )
    {
        //Get the magnitude of the magnetic field and store a local unit Vector for simplified angle calculations.
        const Numeric H_mag = sqrt( rte_mag * rte_mag );
        
        Numeric theta, eta;
        { // Getting angles from coordinate system
        Vector H(3);
            H  = rte_mag;
            H /= H_mag;
        // Direction vector of radiation
        Numeric dx, dy, dz;
            zaaa2cart(dx,dy,dz,R_path_los[0],R_path_los[1]);
        // Radiation path direction as per Mishchenko.
        Vector R_path(3);
            R_path[0] = dx;
            R_path[1] = dy;
            R_path[2] = dz;
        // Vertical polarization direction as per Mishchenko.
        Vector e_v(3);
            e_v[0] =  cos(R_path_los[0] * DEG2RAD) * cos(R_path_los[1] * DEG2RAD);
            e_v[1] =  cos(R_path_los[0] * DEG2RAD) * sin(R_path_los[1] * DEG2RAD);
            e_v[2] = -sin(R_path_los[0] * DEG2RAD);
        // Horizontal polarization direction as per Mishchenko.
        Vector e_h(3);
            e_h[0] = -sin(R_path_los[1] * DEG2RAD);
            e_h[1] =  cos(R_path_los[1] * DEG2RAD);
            e_h[2] =  0;
        // Get the coordinate system used by Lenoir.
        Vector tmp(3);
            proj(tmp, R_path, H);
        Vector R11(3);
            R11  = H;
            R11 -= tmp;
            R11 *= sqrt(R11*R11);
        Vector R22(3);
            cross3(R22, R11, R_path);
        // Test if the rotation is clockwise or counterclockwise.
        const Numeric eta_test = vector_angle(R22, e_h);
        // Find the angle between Mishchenko vertical/horizontal and Lenoir vertical/horizontal
        eta      = (eta_test > 90.0)?-vector_angle(R22, e_v):vector_angle(R22, e_v);

        // Angle between radiation propagation and magnetic field for determining how the radiation is polarized..
        theta = vector_angle(H, R_path);
        }
        
        Numeric DF, RS; // Delta Frequency and Relative Strength

        // For all species
        for(Index II = 0; II<abs_species.nelem(); II++)
        {
            // Reinitialize every loop to empty the set.
            ArrayOfLineRecord temp_abs_lines_sm, temp_abs_lines_sp, //sigma minus, sigma plus
                              temp_abs_lines_pi, temp_abs_lines_dt; // pi, default

            // If the species isn't Zeeman, look at the next species
            if(!is_zeeman(abs_species[II])) continue;
            // Else loop over all the lines in the species.
            for (Index ii = 0; ii< abs_lines_per_species[II].nelem(); ii++)
            {
                // local LineRecord
                LineRecord temp_LR = abs_lines_per_species[II][ii];
                const Rational J  = temp_LR.Lower_J();
                const Rational N  = temp_LR.Lower_N();
                const Numeric S   = isotopologue_quantum.getParam(temp_LR.Species(), temp_LR.Isotopologue(), 1);
                const Numeric GS  = isotopologue_quantum.getParam(temp_LR.Species(), temp_LR.Isotopologue(), 0);
                const Index DJ    = J - temp_LR.Upper_J();
                const Index DN    = N - temp_LR.Upper_N();
                Numeric RS_sum    = 0; //Sum relative strength (which ought be close to one by the end)
                // Only look at lines which have no change in the main rotational number
                
                if (J != -1 && N != -1 && S>0 ) // This means the lines are considered erroneous.
                {
                    if ( true )
                    {

                        for ( Rational M = -J+DJ; M<=J-DJ; M++ )
                        {
                            /*
                                Note that:
                                    sp := sigma plus,  which means DM =  1
                                    sm := sigma minus, which means DM = -1
                                    pi := planar,      which means DM =  0
                            */
                            if ( DJ ==  1 )
                            { // Then all DM transitions possible for all M
                                DF =  FREQUENCY_CHANGE(N, M, J, S, DJ, -1, DN, H_mag, GS);
                                RS = RELATIVE_STRENGTH(N, M, J, DJ, -1);
                                RS_sum += RS;
                                temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                temp_abs_lines_sm.push_back(temp_LR);
                                
                                DF =  FREQUENCY_CHANGE(N, M, J, S, DJ,  0, DN, H_mag, GS);
                                RS = RELATIVE_STRENGTH(N, M, J, DJ,  0);
                                RS_sum += RS;
                                
                                temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                temp_abs_lines_pi.push_back(temp_LR);
                                
                                DF =  FREQUENCY_CHANGE(N, M, J, S, DJ, +1, DN, H_mag, GS);
                                RS = RELATIVE_STRENGTH(N, M, J, DJ, +1);
                                RS_sum += RS;
                                temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                temp_abs_lines_sp.push_back(temp_LR);
                            }
                            else if ( DJ ==  0 )
                            { // Then all DM transitions possible for all M
                            DF =  FREQUENCY_CHANGE(N, M, J, S, DJ, -1, DN, H_mag, GS);
                            RS = RELATIVE_STRENGTH(N, M, J, DJ, -1);
                            RS_sum += RS;
                            temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                            temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                            temp_abs_lines_sm.push_back(temp_LR);
                            
                            if( ! (M == 0) )
                            {
                                DF =  FREQUENCY_CHANGE(N, M, J, S, DJ,  0, DN, H_mag, GS);
                                RS = RELATIVE_STRENGTH(N, M, J, DJ,  0);
                                RS_sum += RS;
                                temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                temp_abs_lines_pi.push_back(temp_LR);
                            }
                            
                            DF =  FREQUENCY_CHANGE(N, M, J, S, DJ, +1, DN, H_mag, GS);
                            RS = RELATIVE_STRENGTH(N, M, J, DJ, +1);
                            RS_sum += RS;
                            temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                            temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                            temp_abs_lines_sp.push_back(temp_LR);
                            }
                            else if ( DJ == -1 )
                            { // Then certain M results in blocked DM transitions
                                if ( M == -J + DJ && M!=0 )
                                { // Lower limit M only allows DM = 1
                                    DF =  FREQUENCY_CHANGE(N, M, J, S, DJ, +1, DN, H_mag, GS);
                                    RS = RELATIVE_STRENGTH(N, M, J, DJ, +1);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sp.push_back(temp_LR);
                                    
                                }
                                else if ( M == -J + DJ + 1 && M!=0 )
                                { // Next to lower limit M can only allow DM = 1, 0
                                    DF =  FREQUENCY_CHANGE(N, M, J, S, DJ, +1, DN, H_mag, GS);
                                    RS = RELATIVE_STRENGTH(N, M, J, DJ, +1);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sp.push_back(temp_LR);

                                    DF =  FREQUENCY_CHANGE(N, M, J, S, DJ,  0, DN, H_mag, GS);
                                    RS = RELATIVE_STRENGTH(N, M, J, DJ,  0);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_pi.push_back(temp_LR);
                                }
                                else if ( M ==  J - DJ - 1 && M!=0 )
                                { // Next to upper limit M can only allow DM = 0, -1
                                    DF =  FREQUENCY_CHANGE(N, M, J, S, DJ,  0, DN, H_mag, GS);
                                    RS = RELATIVE_STRENGTH(N, M, J, DJ,  0);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_pi.push_back(temp_LR);

                                    DF =  FREQUENCY_CHANGE(N, M, J, S, DJ, -1, DN, H_mag, GS);
                                    RS = RELATIVE_STRENGTH(N, M, J, DJ, -1);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sm.push_back(temp_LR);
                                }
                                else if ( M == J - DJ && M!=0 )
                                { // Upper limit M only allow DM = -1
                                    DF =  FREQUENCY_CHANGE(N, M, J, S, DJ, -1, DN, H_mag, GS);
                                    RS = RELATIVE_STRENGTH(N, M, J, DJ, -1);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sm.push_back(temp_LR);
                                }
                                else if( (-J + DJ + 1) ==  (J - DJ - 1) && M == 0)
                                { // Special case for N=1, J=0, M=0. Only allows DM = 0
                                    DF =  FREQUENCY_CHANGE(N, M, J, S, DJ,  0, DN, H_mag, GS);
                                    RS = RELATIVE_STRENGTH(N, M, J, DJ,  0);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_pi.push_back(temp_LR);
                                }
                                else
                                { // All DM transitions are possible for these M(s)
                                    DF =  FREQUENCY_CHANGE(N, M, J, S, DJ, +1, DN, H_mag, GS);
                                    RS = RELATIVE_STRENGTH(N, M, J, DJ, +1);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sp.push_back(temp_LR);

                                    DF =  FREQUENCY_CHANGE(N, M, J, S, DJ,  0, DN, H_mag, GS);
                                    RS = RELATIVE_STRENGTH(N, M, J, DJ,  0);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_pi.push_back(temp_LR);

                                    DF =  FREQUENCY_CHANGE(N, M, J, S, DJ, -1, DN, H_mag, GS);
                                    RS = RELATIVE_STRENGTH(N, M, J, DJ, -1);
                                    RS_sum += RS;
                                    temp_LR.setF(  abs_lines_per_species[II][ii].F()  + DF );
                                    temp_LR.setI0( abs_lines_per_species[II][ii].I0() * RS );
                                    temp_abs_lines_sm.push_back(temp_LR);
                                }
                            }
                            else
                            { // The tests above failed and catastrophe follows
                                throw runtime_error("If this happens, something is horribly wrong... did some of the tests above fail?");
                            }
                        }
                        if (abs(RS_sum-1.)>1e-4) //Reasonable confidence?
                        {
                            ostringstream os;
                            os << "The sum of relative strengths is not close to one. This is severly problematic and"
                            " you should look into why this happens.\nIt is currently " << RS_sum << " with DJ: "<<DJ<<", DN: "<<DN<<"\n";
                            throw runtime_error(os.str());
                        }
                    }
                    else
                        temp_abs_lines_dt.push_back(temp_LR);
                }
                else
                    temp_abs_lines_dt.push_back(temp_LR);
               
            }

            { // Sort ArrayOfLineRecord(s) by frequency [low → high].
                sort(temp_abs_lines_pi.begin(), temp_abs_lines_pi.end(), sortF);
                sort(temp_abs_lines_sp.begin(), temp_abs_lines_sp.end(), sortF);
                sort(temp_abs_lines_sm.begin(), temp_abs_lines_sm.end(), sortF);
                sort(temp_abs_lines_dt.begin(), temp_abs_lines_dt.end(), sortF);
            }
            
            Tensor3 part_abs_mat((*f_grid_pointer).nelem(), 4, 4);
            
            // Add Pi contribution to final abs_mat_per_species
            Part_Return_Zeeman( part_abs_mat, abs_species, abs_lineshape,
                              temp_abs_lines_pi, isotopologue_ratios,
                              abs_vmrs, abs_p, abs_t, *f_grid_pointer,
                              theta, eta, 0, II, verbosity );
            abs_mat_per_species(II, joker, joker, joker) += part_abs_mat;
        
            // Add Sigma minus contribution to final abs_mat_per_species
            Part_Return_Zeeman( part_abs_mat, abs_species, abs_lineshape,
                              temp_abs_lines_sm, isotopologue_ratios,
                              abs_vmrs, abs_p, abs_t, *f_grid_pointer,
                              theta, eta, -1, II, verbosity );
            abs_mat_per_species(II, joker, joker, joker) += part_abs_mat;
            
            // Add Sigma plus contribution to final abs_mat_per_species
            Part_Return_Zeeman( part_abs_mat, abs_species, abs_lineshape,
                                temp_abs_lines_sp, isotopologue_ratios,
                                abs_vmrs, abs_p, abs_t, *f_grid_pointer,
                                theta, eta, 1, II, verbosity );
            abs_mat_per_species(II, joker, joker, joker) += part_abs_mat;
            
            // Add Default contribution to final abs_mat_per_species
            Part_Return_Zeeman( part_abs_mat, abs_species, abs_lineshape,
                              temp_abs_lines_dt, isotopologue_ratios,
                              abs_vmrs, abs_p, abs_t, *f_grid_pointer,
                              theta, eta, 1023, II, verbosity );
            abs_mat_per_species(II, joker, joker, joker) += part_abs_mat;
        }
    }
    else // if the magnetic field is ignored
    {
        for(Index II = 0; II<abs_species.nelem(); II++)
        {
            // If the species isn't Zeeman, look at the next species.
            if(!is_zeeman(abs_species[II])) continue;
            
            Tensor3 part_abs_mat((*f_grid_pointer).nelem(), 4, 4);
            Part_Return_Zeeman(   part_abs_mat, abs_species, abs_lineshape,
                                abs_lines_per_species[II], isotopologue_ratios,
                                abs_vmrs, abs_p, abs_t, *f_grid_pointer,
                                0,0,1023, II, verbosity );
            abs_mat_per_species(II, joker, joker, joker) += part_abs_mat;
        }
    }
}
