#ifndef tokval_io
#define tokval_io

#include <ostream>

class TokVal;

class TokValPrinter {
  const TokVal& ref;
public:
  TokValPrinter(const TokVal& v) : ref(v) {}
  friend std::ostream& operator<<(std::ostream&, const TokValPrinter&);
};

#endif  // tokval_io
