// Template Numerical Toolkit (TNT) for Linear Algebra
//
// BETA VERSION INCOMPLETE AND SUBJECT TO CHANGE
// Please see http://math.nist.gov/tnt for updates
//
// R. Pozo
// Mathematical and Computational Sciences Division
// National Institute of Standards and Technology


#ifndef VECADAPTOR_H
#define VECADAPTOR_H

#include <cstdlib>
#include <iostream>
#include <strstream>
#include <cassert>

#include "subscrpt.h"

#ifdef TNT_USE_REGIONS
#include "region1d.h"
#endif

namespace TNT
{

//  see "tntreq.h" for TNT requirements for underlying vector
//  class.  This need NOT be the STL vector<> class, but a subset
//  that provides minimal services.
//
//  This is a container adaptor that provides the following services.
//
//      o)  adds 1-offset operator() access ([] is always 0 offset)
//      o)  adds TNT_BOUNDS_CHECK to () and []
//      o)  adds initialization from strings, e.g.  "1.0 2.0 3.0";
//      o)  adds newsize(N) function (does not preserve previous values)
//      o)  adds dim() and dim(1)
//      o)  adds free() function to release memory used by vector
//      o)  adds regions, e.g. A(Index(1,10)) = ... 
//      o)  add getVector() method to return adapted container
//      o)  adds simple I/O for ostreams

template <class BBVec>
class Vector_Adaptor
{

  public:
    typedef   typename BBVec::value_type T;
    typedef         T   value_type;
    typedef         T   element_type;
    typedef         T*  pointer;
    typedef         T*  iterator;
    typedef         T&  reference;
    typedef const   T*  const_iterator;
    typedef const   T&  const_reference;
    
    Subscript lbound() const { return 1; }

  protected:
    BBVec v_;
    T* vm1_;

  public:

    Subscript size() const { return v_.size(); }

    // These were removed so that the ANSI C++ valarray class
    // would work as a possible storage container.
    //
    // SAB 30.01.2000 Commented them in again, I don´t use std::valarray,  
    //                but std::vector     
    //
    iterator begin() { return v_.begin();}
    //iterator begin() { return &v_[0];}
    
    iterator end()   { return v_.end(); }
    //iterator end()   { return &v_[0] + v_.size(); }
    
    const_iterator begin() const { return v_.begin();}
    //const_iterator begin() const { return &v_[0];}
    
    const_iterator end()  const { return v_.end(); }
    //const_iterator end()  const { return &v_[0] + v_.size(); }

    // SAB 30.01.2000 Added push_back method, for compatibility with std::vector.
    void push_back(const T& x) 
    {
        v_.push_back(x);
	// vm1_ needs to be set after push_back!! 
	// (push_back changes memory allocation)
	vm1_ = &(v_[0]) - 1;
    }

    // SAB 11.02.2000 Added clear method, for compatibility with std::vector.
    void clear() 
    {
        v_.clear();
	// vm1_ is now no longer valid.
	vm1_ = NULL;
    }

    // SAB 22.03.2000 Added erase method, for compatibility with std::vector.
    void erase(iterator q)
    {
        v_.erase(q);
	// To be on the safe side, re-set vm1_
	vm1_ = &(v_[0]) - 1;
    }

    BBVec& getVector() { return v_; }
    Subscript dim() const { return v_.size(); }
    Subscript dim(Subscript i)
    {
#ifdef TNT_BOUNDS_CHECK
        assert(i==TNT_BASE_OFFSET);
#endif
        return (i==TNT_BASE_OFFSET ? v_.size() : 0 );
    }
    Vector_Adaptor() : v_() {};
    Vector_Adaptor(const Vector_Adaptor<BBVec> &A) : v_(A.v_) 
    { 
        vm1_ = ( v_.size() > 0 ? &(v_[0]) -1 : NULL); 

    } 

    Vector_Adaptor(Subscript N, /*const*/ char *s) : v_(N) 
    {
        istrstream ins(s);
        for (Subscript i=0; i<N; i++)
            ins >> v_[i] ;

        vm1_ = ( v_.size() > 0 ? &(v_[0]) -1 : NULL); 
    }; 

    Vector_Adaptor(Subscript N, const T& value = T()) : v_(N)
    {
        for (Subscript i=0; i<N; i++)
             v_[i]  = value;
        
        vm1_ = ( v_.size() > 0 ? &(v_[0]) -1 : NULL); 
    }

    Vector_Adaptor(Subscript N, const T* values) : v_(N)
    {
        for (Subscript i=0; i<N; i++)
             v_[i]  = values[i];
        vm1_ = ( v_.size() > 0 ? &(v_[0]) -1 : NULL); 
    } 
    Vector_Adaptor(const BBVec & A) : v_(A) 
    {
        vm1_ = ( v_.size() > 0 ? &(v_[0]) -1 : NULL); 
    }

    // NOTE: this assumes that BBVec(0) constructor creates an 
    //  null vector that does not take up space...  It would be
    //  great to require that BBVec have a corresponding free()
    //  function, but in particular STL vectors do not.
    //
    Vector_Adaptor<BBVec>& free()
    {
        return *this = Vector_Adaptor<BBVec>(0);
    }

    Vector_Adaptor<BBVec>& operator=(const Vector_Adaptor<BBVec> &A) 
    { 
        v_ = A.v_ ; 
        vm1_ = ( v_.size() > 0 ? &(v_[0]) -1 : NULL); 
        return *this;
    }

    Vector_Adaptor<BBVec>& newsize(Subscript N)
    {
        // NOTE: this is not as efficient as it could be
        // but to retain compatiblity with STL interface
        // we cannot assume underlying implementation
        // has a newsize() function.

        return *this = Vector_Adaptor<BBVec>(N);

    }

    Vector_Adaptor<BBVec>& operator=(const T &a) 
    {
        Subscript i;
        Subscript N = v_.size();    
        for (i=0; i<N; i++)
            v_[i] = a;

        return *this;
    }

  // SAB 22.05.2000 Added second argument for resizing with
  // initialization to value c.
    Vector_Adaptor<BBVec>& resize(Subscript N, T c=T()) 
    { 
        if (N == size()) return *this;

        Vector_Adaptor<BBVec> tmp(N);
        Subscript n =  (N < size() ? N : size());  // min(N, size());
        Subscript i;

        for (i=0; i<n; i++)
            tmp[i] = v_[i];

	// SAB: Initialize new elements with c:
        for (i=n; i<N; i++)
            tmp[i] = c;
            

        return (*this = tmp);

    }


    reference operator()(Subscript i)
    { 
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i<=dim());
#endif
        return vm1_[i]; 
    }

    const_reference operator()(Subscript i) const
    { 
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i<=dim());
#endif
        return vm1_[i]; 
    }

    reference operator[](Subscript i)
    { 
#ifdef TNT_BOUNDS_CHECK
        assert(0<=i);
        assert(i<dim());
#endif
        return v_[i]; 
    }

    const_reference operator[](Subscript i) const
    { 
#ifdef TNT_BOUNDS_CHECK
        assert(0<=i);
        assert(i<dim());
#endif
        return v_[i]; 
    }


#ifdef TNT_USE_REGIONS
    // "index-aware" features, all of these are 1-based offsets

    typedef Region1D<Vector_Adaptor<BBVec> > Region;

    typedef const_Region1D< Vector_Adaptor<BBVec> > const_Region;

    Region operator()(const Index1D &I)
    {   return Region(*this, I); }

    Region operator()(const Subscript i1, Subscript i2)
    {   return Region(*this, i1, i2); }

    const_Region operator()(const Index1D &I) const
    {   return const_Region(*this, I); }

    const_Region operator()(const Subscript i1, Subscript i2) const
    {   return const_Region(*this, i1, i2); }
#endif
// TNT_USE_REGIONS


};

#include <iostream>

template <class BBVec>
std::ostream& operator<<(std::ostream &s, const Vector_Adaptor<BBVec> &A)
{
    Subscript M=A.size();

    s << M << endl;
    for (Subscript i=1; i<=M; i++)
            s << A(i) << endl;
    return s;
}

template <class BBVec>
std::istream& operator>>(std::istream &s, Vector_Adaptor<BBVec> &A)
{
    Subscript N;
    
    s >> N;

    A.resize(N);

    for (Subscript i=1; i<=N; i++)
        s >> A(i);

    return s;
}

} // namespace TNT

#endif
