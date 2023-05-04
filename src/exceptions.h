/*!
  \file   exceptions.h
  \brief  The declarations of all the exception classes.

  Dedicated exception classes are used only by the parser. Otherwise,
  only runtime_error is used. Furthermore, errors that corresponds to
  bugs in the program rather than incorrect input should be handled by
  assertions.

  \author Stefan Buehler
  \date   1999-09-24
*/

#ifndef exceptions_h
#define exceptions_h

#include <stdexcept>

#include "matpack_concepts.h"
#include "mystring.h"

// Special stuff for the parser:

class ParseError : public std::runtime_error {
 public:
  ParseError(const String& s = "", String f = "", Index l = 0, Index c = 0)
      : runtime_error(s),
        mFile(std::move(f)),
        mLine(l),
        mColumn(c) { /* Nothing to do here. */
  }

  [[nodiscard]] String file() const { return mFile; }
  [[nodiscard]] Index line() const { return mLine; }
  [[nodiscard]] Index column() const { return mColumn; }

 private:
  /** Filename associated with this part of the text */
  String mFile;
  /** Line where the error occured. */
  Index mLine;
  /** Column where the error occured. */
  Index mColumn;
};

class Eot : public ParseError {
 public:
  Eot(const String& s = "", const String& f = "", Index l = 0, Index c = 0)
      : ParseError(s, f, l, c) { /* Nothing to do here. */
  }
};

class UnexpectedChar : public ParseError {
 public:
  UnexpectedChar(const String& s = "",
                 const String& f = "",
                 Index l = 0,
                 Index c = 0)
      : ParseError(s, f, l, c) { /* Nothing to do here. */
  }
};

class IllegalLinebreak : public ParseError {
 public:
  IllegalLinebreak(const String& s = "",
                   const String& f = "",
                   Index l = 0,
                   Index c = 0)
      : ParseError(s, f, l, c) { /* Nothing to do here. */
  }
};

class UnknownMethod : public ParseError {
 public:
  UnknownMethod(const String& s = "",
                const String& f = "",
                Index l = 0,
                Index c = 0)
      : ParseError(s, f, l, c) { /* Nothing to do here. */
  }
};

class UnknownWsv : public ParseError {
 public:
  UnknownWsv(const String& s = "",
             const String& f = "",
             Index l = 0,
             Index c = 0)
      : ParseError(s, f, l, c) { /* Nothing to do here. */
  }
};

class WsvAlreadyExists : public ParseError {
 public:
  WsvAlreadyExists(const String& s = "",
                   const String& f = "",
                   Index l = 0,
                   Index c = 0)
      : ParseError(s, f, l, c) { /* Nothing to do here. */
  }
};

class WrongWsvGroup : public ParseError {
 public:
  WrongWsvGroup(const String& s = "",
                const String& f = "",
                Index l = 0,
                Index c = 0)
      : ParseError(s, f, l, c) { /* Nothing to do here. */
  }
};

#endif  // exceptions_h
