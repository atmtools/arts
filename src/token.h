#ifndef token_h
#define token_h

#include "vecmat.h"

/** The different token value types. These are the types that keyword
    parameters in the controlfile can have. */
enum TokValType { str_,    int_,    num_,
	          strvec_, intvec_, numvec_,
                  undefined_ };

/** This stores arbitrary token values and remembers the type. Only
    the correct type can be extracted again. */
class TokVal {
public:

  /** Default Constructor. (Sets type to undefined_) */
  TokVal() {
    mtype = undefined_;
  }

  /** To set TokVal from string (C - style). */
  TokVal(const char s[]) {
    mtype = str_;
    ms = s;
  }

  /** To set TokVal from string (C++ - style). */
  TokVal(const string& s) {
    mtype = str_;
    ms = s;
  }

  /** To set TokVal from an integer. */
  TokVal(int n) {
    mtype = int_;
    mn = n;
  }

  /** To set TokVal from a Numeric. */
  TokVal(Numeric x) {
    mtype = num_;
    mx = x;
  }

  /** To set TokVal from an array of strings. */
  TokVal(ARRAY<string> sv) {
    mtype = strvec_;
    msv = sv;
  }

  /** To set TokVal from an array of integers. */
  TokVal(ARRAY<int> nv) {
    mtype = intvec_;
    mnv = nv;
  }

  /** To set TokVal from an array of Numerics. */
  TokVal(ARRAY<Numeric> xv) {
    mtype = numvec_;
    mxv = xv;
  }

  // Conversion functions to return TokVal for the 6 different types: 
  
  /** Return string. */
  operator string() const;
  /** Return int. */
  operator int() const;
  /** Return Numeric. */
  operator Numeric() const;

  /** Return array of strings. */
  operator ARRAY<string>() const;
  /** Return array of integers. */
  operator ARRAY<int>() const;
  /** Return array of Numerics. */
  operator ARRAY<Numeric>() const;

  /** Output operator. */
  friend ostream& operator<<(ostream& os, const TokVal& a);

private:
  TokValType mtype;
  string       ms;
  int          mn;
  Numeric      mx;   
  ARRAY<string>  msv;
  ARRAY<int>     mnv;
  ARRAY<Numeric> mxv;
};


// typedef ARRAY<TokValType> TokValTypeVector;
// typedef ARRAY<TokVal>     TokValVector;



#endif  // token_h
