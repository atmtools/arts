/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>

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

// Special stuff for the parser:

class ParseError : public runtime_error {
public:
  ParseError( const string& s="",
	      const string& f="",
	      int l = 0,
	      int c = 0 ) :
    runtime_error(s),
    mFile(f),
    mLine(l),
    mColumn(c) { /* Nothing to do here. */ }

  virtual string file()		   const { return mFile; }
  virtual int line()   		   const { return mLine; }
  virtual int column() 		   const { return mColumn; }

private:
  /** Filename associated with this part of the text */
  string mFile;
  /** Line where the error occured. */
  int mLine;
  /** Column where the error occured. */
  int mColumn;
};


class Eot : public ParseError {
public:
  Eot( const string& s="",
       const string& f="",
       int l = 0,
       int c = 0 ) :
    ParseError(s,f,l,c) { /* Nothing to do here. */ }
};

class UnexpectedChar : public ParseError {
public:
  UnexpectedChar( const string& s="",
		  const string& f="",
		  int l = 0,
		  int c = 0 ) :
    ParseError(s,f,l,c) { /* Nothing to do here. */ }
};

class IllegalLinebreak : public ParseError {
public:
  IllegalLinebreak( const string& s="",
		    const string& f="",
		    int l = 0,
		    int c = 0 ) :
    ParseError(s,f,l,c) { /* Nothing to do here. */ }
};

class UnknownMethod : public ParseError {
public:
  UnknownMethod( const string& s="",
		 const string& f="",
		 int l = 0,
		 int c = 0 ) :
    ParseError(s,f,l,c) { /* Nothing to do here. */ }
};

class UnknownWsv : public ParseError {
public:
  UnknownWsv( const string& s="",
	      const string& f="",
	      int l = 0,
	      int c = 0 ) :
    ParseError(s,f,l,c) { /* Nothing to do here. */ }
};

class WrongWsvGroup : public ParseError {
public:
  WrongWsvGroup( const string& s="",
		 const string& f="",
		 int l = 0,
		 int c = 0 ) :
    ParseError(s,f,l,c) { /* Nothing to do here. */ }
};

class UnexpectedKeyword : public ParseError {
public:
  UnexpectedKeyword( const string& s="",
		     const string& f="",
		     int l = 0,
		     int c = 0 ) :
    ParseError(s,f,l,c) { /* Nothing to do here. */ }
};




#endif // exceptions_h
