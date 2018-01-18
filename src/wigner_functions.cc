/* Copyright (C) 2012
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

// #include "wigxjpf/inc/wigxjpf.h"
#include "wigner_functions.h"
#include <stdexcept>
#include <sstream>
#include <algorithm>

/*
extern "C"
{
  void wig_table_init(int, int);
  void wig_temp_init(int);
  void wig_temp_free();
  void wig_table_free();
  double wig3jj(int, int, int, int, int, int);
}*/

/*!
  Equation 34.2.4 in http://dlmf.nist.gov/34.2
  
  Some test cases to simplify the calculations are made for m_i=0.
*/
Numeric wigner3j(const Rational j1,const Rational j2,const Rational j3,
                 const Rational m1,const Rational m2,const Rational m3)
{
/*
  const int a = int((j1.Denom() == 2)? j1.Nom() : 2 * j1.Nom()),
            b = int((j2.Denom() == 2)? j2.Nom() : 2 * j2.Nom()), 
            c = int((j3.Denom() == 2)? j3.Nom() : 2 * j3.Nom()),
            d = int((m1.Denom() == 2)? m1.Nom() : 2 * m1.Nom()),
            e = int((m2.Denom() == 2)? m2.Nom() : 2 * m2.Nom()),
            f = int((m3.Denom() == 2)? m3.Nom() : 2 * m3.Nom());
  const int g = std::max(std::max(std::max(std::max(std::max(a, b), c), d), e), f);
  double h;
  
  wig_table_init(g, 3);
  wig_temp_init(g);
            
  h = wig3jj(a, b, c, d, e, f);
  
  wig_temp_free();
  wig_table_free();
  
  return Numeric(h);*/

    Rational J = j1 + j2 + j3;
    J.Simplify();
    
    {// TEST Area
        if( !( ( m1 + m2 + m3 ).toIndex() == 0 ) )
            return 0;
        
        if( !( triangular_inequality( j1, j2, j3 ) ) )
            return 0;
        
        if( !( ( J ).Denom() == 1 ) )
            return 0;
        
        if( !( abs( m1 ) <= j1 && abs( m2 ) <= j2 && abs( m3 ) <= j3 ) )
            return 0;
        
        if(m1.toIndex() ==0 && m2.toIndex() == 0 && m3.toIndex() ==0)
        {
            if( !( (( J ).toIndex() % 2) == 0 ) )
                return 0;
            else
            {
                return ((((J/2)%2)==0)?(1.):(-1.)) * sqrt( factorials(
                        {(J-2*j1).toIndex(),
                         (J-2*j2).toIndex(),
                         (J-2*j3).toIndex(),
                         (J/2).toIndex(),
                         (J/2).toIndex()},
                        {(J+1).toIndex(),
                         (J/2-j1).toIndex(),
                         (J/2-j1).toIndex(),
                         (J/2-j2).toIndex(),
                         (J/2-j2).toIndex(),
                         (J/2-j3).toIndex(),
                         (J/2-j3).toIndex()}));
            }
        }
    }// End of TEST Area
    
    
    
    // This is the first row of equation 32.2.4 on the webpage specified.
    const Numeric constant = ((((j1-j2-m3)%2)==0)?(1.):(-1.)) * sqrt(factorials(
                              {(+j1+j2-j3).toIndex(),
                               (+j1-j2+j3).toIndex(),
                               (-j1+j2+j3).toIndex(),
                               (j1+m1).toIndex(),
                               (j1-m1).toIndex(),
                               (j2+m2).toIndex(),
                               (j2-m2).toIndex(),
                               (j3+m3).toIndex(),
                               (j3-m3).toIndex()},
                              {(J+1).toIndex()}));
    
    Numeric sum=0; Rational s=0;
        
    // This is the sum of the same function.
    while(      (j1+j2-j3-s).toIndex() >= 0 && 
                (j1-m1-s).toIndex()    >= 0 && 
                (j2+m2-s).toIndex()    >= 0 )
    {
        if(     (j3-j2+m1+s).toIndex() >= 0 && 
                (j3-j1-m2+s).toIndex() >= 0 )
        {
            sum += ((((s)%2)==0)?(1.):(-1.)) * factorials(
                    {1},
                    {s.toIndex(),
                     (j1+j2-j3-s).toIndex(),
                     (j1-m1-s).toIndex(),
                     (j2+m2-s).toIndex(),
                     (j3-j2+m1+s).toIndex(),
                     (j3-j1-m2+s).toIndex()});
                }
        s++;
    }
    
    return constant*sum;
}


/*!
  Equation 34.4.2 in http://dlmf.nist.gov/34.4 may be invalid.
  
  Using Equation 4 in http://mathworld.wolfram.com/Wigner6j-Symbol.html instead
  
  Simplifications are made for J1 = 1.
*/
Numeric wigner6j(const Rational j1,const Rational j2,const Rational j3,
                 const Rational J1,const Rational J2,const Rational J3)
{
    
    if( !( triangular_inequality( j1, j2, j3 ) ) )
        return 0;
    
    if( !( triangular_inequality( j1, J2, J3 ) ) )
        return 0;
    
    if( !( triangular_inequality( J1, j2, J3 ) ) )
        return 0;
    
    if( !( triangular_inequality( J1, J2, j3 ) ) )
        return 0;
        
    //Fancy simplification
    if(J1.toIndex()==1)
        if(j2==J3)
            if(j3 == J2)
                return 2 * ((((j1+j2+j3+1).toIndex()%2)==0)?(1.):(-1.)) * 
                ( j2.toNumeric()*(j2.toNumeric()+1) + j3.toNumeric()*(j3.toNumeric()+1) - j1.toNumeric()*(j1.toNumeric()+1) ) / 
                sqrt( 2*j2.toNumeric()*(2*j2.toNumeric()+1)*(2*j2.toNumeric()+2)*2*j3.toNumeric()*(2*j3.toNumeric()+1)*(2*j3.toNumeric()+2) );
    
    Numeric sum=0; Rational s=0;
    //std::cout<<std::endl<<"Wigner6j (" << j1 <<" "<<j2 <<" "<<j3<<"; " << J1 <<" " << J2 << " "<<J3<<")"<<std::endl;
    // Complicated sum
    while(      (j1+j2+J1+J2-s).toIndex() >= 0 &&
                (j2+j3+J2+J3-s).toIndex() >= 0 &&
                (j3+j1+J3+J1-s).toIndex() >= 0 )
    {
        if( (s-j1-j2-j3).toIndex() >= 0 &&
            (s-j1-J2-J3).toIndex() >= 0 &&
            (s-J1-j2-J3).toIndex() >= 0 &&
            (s-J1-J2-j3).toIndex() >= 0 )
        {
            sum += ((((s)%2)==0)?(1.):(-1.))*factorials(
                    {s.toIndex()+1},
                    {(s-j1-j2-j3).toIndex(),
                     (s-j1-J2-J3).toIndex(),
                     (s-J1-j2-J3).toIndex(),
                     (s-J1-J2-j3).toIndex(),
                     (j1+j2+J1+J2-s).toIndex(),
                     (j2+j3+J2+J3-s).toIndex(),
                     (j3+j1+J3+J1-s).toIndex()});
        }
        s++;
    }
    //std::cout<<"s became " << s-1 << std::endl;
    
    const Numeric tc = sqrt(
    triangle_coefficient(j1,j2,j3) *
    triangle_coefficient(j1,J2,J3) *
    triangle_coefficient(J1,j2,J3) *
    triangle_coefficient(J1,J2,j3)
    );
    
    return tc * sum;
}

// /*!
//  * Equation 34.2.5 in http://dlmf.nist.gov/34.2
//  */
Numeric wigner9j(){return NAN;}//Placeholder


Numeric ECS_wigner(Rational L, Rational Nl, Rational Nk, 
                   Rational Jk_lower, Rational Jl_lower, 
                   Rational Jk_upper, Rational Jl_upper)
{
    const Numeric A1 = wigner3j(Nl,Nk,L,0,0,0);
    if( A1 == 0)
        return 0;
    //std::cout << "input to wigner3j is: (" << Nl <<", " << Nk << ", " << L << ")\n";
    
    const Numeric A2 = wigner6j(L,Jk_upper,Jl_upper,1,Nl,Nk);
    if( A2 == 0)
        return 0;
    
    const Numeric A3 = wigner6j(L,Jk_lower,Jl_lower,1,Nl,Nk);
    if( A3 == 0)
        return 0;
    
    const Numeric A4 = wigner6j(L,Jk_upper,Jl_upper,1,Jl_lower,Jk_lower);
    if( A4 == 0)
        return 0;
    
    //std::cout << "A1="<<A1<<", A2="<<A2 << ", A3="<<A3<<", A4="<<A4<<std::endl;
    
    return A1*A2*A3*A4;
}


Numeric triangle_coefficient(const Rational a, const Rational b, const Rational c)
{    
    Rational nom1=(a+b-c), nom2=(a-b+c), nom3=(-a+b+c), denom1=(a+b+c+1);
    return factorials({nom1.toIndex(),nom2.toIndex(),nom3.toIndex()},
                      {denom1.toIndex()});
}

bool triangular_inequality(const Rational x, const Rational y, const Rational z)
{
    if( abs(x-y) > z )
        return false;
    if( (x+y) < z )
        return false;
    
    return true;
}


// Helper function that returns all primes lower than the input number. 
void primes(ArrayOfIndex& output, const Index input)
{
    output.resize(0);
    
    if(input<2)
        return;
    
    output.push_back(2);
    
    if(input==2) return;
    
    Index ii = 3,JJ=1;
    
    while( ii<=input )
    {
        Index test = 1;
        for(Index jj=0;jj<JJ;jj++)
        {
            // If it is divisible by a previous prime number then it is not a prime number.
            if((ii%output[jj])==0)
            {
                test = 0;
                break;
            }
        }
        // If it is divisible by a previous prime number then it is not a prime number otherwise it is.
        if( test ) 
        {
            output.push_back(ii);
            JJ++;
        }
        
        ii++;
    }
    
}


// Helper function that returns the powers for describing factorial of input
void powers(ArrayOfIndex& output, const ArrayOfIndex primes, const Index input)
{
    Index numbers = primes.nelem();
    
    output.resize(numbers);
    
    for(Index ii=0;ii<numbers;ii++)
    {
        output[ii] = 0;
        Index num = primes[ii];
        while( input/num != 0 )
        {
            output[ii] += input/num;
            num *= primes[ii];
        }
    }
}


// Function that calculates prod(NomFac!) / prod(DenomFac!)
Numeric factorials(const ArrayOfIndex& NomFac, const ArrayOfIndex& DenomFac)
{
    //std::cout<<NomFac<<", "<<DenomFac<<std::endl;
    Index max_nom_ind, max_denom_ind;
    
    {//Test for Errors
        std::ostringstream os;
        Index test_if_input_is_crazy = 0;
        
    if(NomFac.nelem()>0)
    {
        max_nom_ind = max(NomFac);
        if( min(NomFac)<0 )
        {
            test_if_input_is_crazy++;
            os << "Negative values in the factorial nominator.\n This is: " << NomFac << std::endl;
        }
    }
    else
        max_nom_ind = 0;
    if(DenomFac.nelem()>0)
    {
        max_denom_ind = max(DenomFac);
        if( min(DenomFac)<0 )
        {
            test_if_input_is_crazy++;
            os << "Negative values in the factorial denominator.\n This is: " << DenomFac << std::endl;
        }
    }
    else
        max_denom_ind = 0;
    if(test_if_input_is_crazy>0)
        throw std::runtime_error(os.str());
    }

    const Index max_ind = (max_denom_ind<max_nom_ind) ? max_nom_ind : max_denom_ind;
    
    // All primes are based on the largest input
    ArrayOfIndex all_primes, all_powers, temp;
    primes(all_primes, max_ind);
    all_powers.resize(all_primes.nelem());
    
    // Get powers and add to all_powers
    for( Index jj = 0; jj < NomFac.nelem(); jj++ )
    {
        powers(temp,all_primes,NomFac[jj]);
        
        for( Index ii = 0; ii < temp.nelem(); ii++ )
            all_powers[ii] += temp[ii];
    }
    
    // Get powers and subtract from all_powers
    for( Index jj = 0; jj < DenomFac.nelem(); jj++ )
    {
        powers(temp,all_primes,DenomFac[jj]);
        
        for( Index ii = 0; ii < temp.nelem(); ii++ )
            all_powers[ii] -= temp[ii];
    }
    
    Numeric total=1;
    bool still_possible=true; Index num=all_primes.nelem();
    
    //std::cout<<all_primes<<"; "<<all_powers<<std::endl;
    
    // This loop might be a time consumer.
    while( still_possible )
    {
        Index test=1;
        for( Index jj = 0 ; jj < num; jj++)
        {
            if(all_powers[jj]>0)
            {
                if(total > 1e20 && min(all_powers)<0)
                    continue;
                total *= (Numeric) all_primes[jj];
                all_powers[jj]--;
            }
            else if(all_powers[jj]<0)
            {
                if(total < 1e-20 && max(all_powers)>0)
                    continue;
                total /= (Numeric) all_primes[jj];
                all_powers[jj]++;
            }
            else
            {
                test++;
            }
        }
        if(test>=num)
            still_possible=false;
    }
    
    return total;
}

