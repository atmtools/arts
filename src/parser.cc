/* Copyright (C) 2000, 2001 Stefan Buehler <sbuehler@uni-bremen.de>

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

#include "arts.h"
#include <map>
#include "messages.h"
#include "exceptions.h"
#include "file.h"
#include "auto_wsv.h"
#include "methods.h"
#include "parser.h"
#include "wsv_aux.h"

void SourceText::AppendFile(const String& name) 
{
  mSfLine.push_back(mText.nelem());
  mSfName.push_back(name);

  read_text_from_file(mText, name);    
}

void SourceText::AdvanceChar() 
{
  if ( mColumn < mText[mLine].nelem()-1 )
    {
      ++mColumn;
    }
  else
    {
      mLineBreak = true;
      do
	{
	  if (mLine>=mText.nelem()-1)
	    {
	      throw Eot( "",
			 this->File(),
			 this->Line(),
			 this->Column() ); 
	    }
	  else
	    {
	      ++mLine;
	      mColumn = 0;
	    }
	}
      while ( 1 > mText[mLine].nelem() ); // Skip empty lines.
    }
}


void SourceText::AdvanceLine() 
{
  mLineBreak = true;
  mColumn = 0;
  do
    {
      if (mLine>=mText.nelem()-1)
	{
	  throw Eot( "",
		     this->File(),
		     this->Line(),
		     this->Column() ); 
	}
      else
	{
	  ++mLine;
	}
    }
  while ( 1 > mText[mLine].nelem() ); // Skip empty lines.
}


const String& SourceText::File()
{
  Index i    = 0;
  bool   stop = false;

  while ( i<mSfLine.nelem()-1 && !stop )
    {
      if (mLine>=mSfLine[i+1]) ++i;
      else                     stop = true;
    }

  return mSfName[i];
}


Index SourceText::Line()
{
  Index i    = 0;
  bool   stop = false;

  while ( i<mSfLine.nelem()-1 && !stop )
    {
      if (mLine>=mSfLine[i+1]) ++i;
      else                     stop = true;
    }

  return mLine - mSfLine[i] + 1; 
}


void SourceText::Init()
{
  mLine   = 0;
  mColumn = 0;
    
  if ( 1 > mText.nelem() )
    {
      throw Eot( "Empty text!",
		 this->File(),
		 this->Line(),
		 this->Column() ); 
    }
  else
    {
      // Skip empty lines:
      while ( 1 > mText[mLine].nelem() )
	{
	  if (mLine>=mText.nelem()-1)
	    {
	      throw Eot( "",
			 this->File(),
			 this->Line(),
			 this->Column() ); 
	    }
	  else
	    {
	      mLineBreak = true;
	      ++mLine;
	    }
	}
    }
}


ostream& operator << (ostream& os, const SourceText& text)
{
  for (Index i=0; i<text.mText.nelem();++i)
    cout << i
	 << "(" << text.mText[i].nelem() << ")"
	 << ": " << text.mText[i] << '\n';
  return(os);
}


/** Returns true if this character is considered whitespace. This
    includes the comment sign `\#'. This function is used by other
    functions to test for delimiting whitespace. 

    The whitespace cases implemented here must be consistent with
    eat_whitespace! 
    @see eat_whitespace  */
bool is_whitespace(const char c)
{
  switch (c) 
    {
    case ' ':
    case '\t':
    case '#':
      return true;
      break;
    default:
      return false;
      break;
    }
}

/** Eats whitespace. Comments are a special case of
    whitespace. Everything from the `\#' to the end of the line is
    eaten. 
  
    The whitespace cases implemented here must be consistent with
    is_whitespace! 
    @see is_whitespace    
  
    @exceptions UnexpectedChar Non-whitespace character encountered. */
void eat_whitespace(SourceText& text)
{
  char dummy;

  while (is_whitespace(dummy=text.Current())) 
    {
      switch (dummy)
	{
	case ' ':
	case '\t':
	  text.AdvanceChar();
	  break;
	case '#':
	  text.AdvanceLine();
	  break;
	default:
	  {
	    ostringstream os;
	    os << "Expected whitespace, but got `" << dummy << "'.";
	    throw UnexpectedChar( os.str(),
				  text.File(),
				  text.Line(),
				  text.Column() ); 
	    break;
	  }
	}
    }
}


/** Reads name of method, keyword, or workspace variable.
  
    These names may consist only of letters (case matters!), numbers, and
    underscores.
    Line break or any other character ends the name.
  
    Whitespace has to have been eaten before. Scanns text for the
    name, starting at position specified by line and column. */
void read_name(String& name, SourceText& text)
{
  bool stop = false;
  name = "";

  while (!stop) 
    {
      char dummy = text.Current();

      if ( isalnum(dummy) || '_'==dummy )
	{
	  name += dummy;
	  // AdvanceChar sets LineBreak if a line break occured.
	  text.LineBreak() = false;
	  text.AdvanceChar();
	  if ( text.LineBreak() ) stop = true;
	}
      else
	{
	  stop = true;
	}
    }

  //  cout << "Name: " << name << '\n';
}

/** Make sure that the current character is equal to c and go to the
    next character.
  
    @exception UnexpectedChar The character is not right. */
void assertain_character(char c, SourceText& text)
{
  if ( c != text.Current() )
    {
      ostringstream os;
      os << "Expected `" << c << "', but got `" << text.Current() << "'.";
      throw UnexpectedChar( os.str(),
			    text.File(),
			    text.Line(),
			    text.Column() ); 
    }
  
  text.AdvanceChar();
}

/** Reads a String, complete with quotation marks. Whitespace has to
    have been eaten before, that is, the current character must be
    the quotation mark (").
    Quotation marks inside Strings are currently not possible.
  
    Line breaks inside Strings are not allowed. 

    @exception IllegalLinebreak An illegal linebreak has occured. */
void parse_String(String& res, SourceText& text)
{
  bool stop = false;
  res = "";

  text.LineBreak() = false;
  assertain_character('"',text);
  if ( text.LineBreak() )
    throw IllegalLinebreak( "Line break before end of String.",
			    text.File(),
			    text.Line(),
			    text.Column() ); 

  while (!stop) 
    {
      char dummy = text.Current();
      if ( dummy != '"' )
	{
	  res += dummy;
	  text.AdvanceChar();

	  if ( text.LineBreak() )
	    throw IllegalLinebreak( "Line break before end of String.",
				    text.File(),
				    text.Line(),
				    text.Column() ); 
	}
      else
	{
	  stop = true;
	  text.AdvanceChar();
	}
    }
}

/** Reads an integer. Whitespace has to
    have been eaten before, that is, the current character must be
    a number or `+' or `-'.
  
    Whitespace or line breaks terminate the scanning!
    There are no whitespaces allowed anywhere, consisten with ANSI C
    scanf. 

    @exception IllegalLinebreak An illegal linebreak has occured. 
    @exception UnexpectedChar Unexpected character encountered. */
void read_integer(String& res, SourceText& text)
{
  bool stop = false;
  res = "";
  char dummy;
  text.LineBreak() = false;

  dummy = text.Current();
  if ( '+' == dummy || '-' == dummy )
    {
      res += dummy;
      text.AdvanceChar();
      if ( text.LineBreak() )
	throw IllegalLinebreak( "Line break after sign.",
				text.File(),
				text.Line(),
				text.Column() ); 
    }

  if (!isdigit(text.Current()))
    {
      ostringstream os;
      os << "Expected digit, but got `" << text.Current() << "'.";
      throw UnexpectedChar(os.str(),
			   text.File(),
			   text.Line(),
			   text.Column());
    }

  while (!stop) 
    {
      char dummy = text.Current();
      if ( isdigit(dummy) )
	{
	  res += dummy;
	  text.AdvanceChar();
	  if ( text.LineBreak() ) stop = true;
	}
      else
	{
	  stop = true;
	}
    }
}

/** Reads a floating point number. Whitespace has to
    have been eaten before, that is, the current character must be
    a number or `+' or `-'.
  
    Example numbers: 23.,  1.0,  -.3,  -3.3e5,  +3e8,  1.0E-9
  
    Illegal numbers: ., 3e, e3, 2e-
  
    Whitespace is not allowed inside the number.
    Line breaks or whitespace terminates the scanning. 

    @exception IllegalLinebreak Illegal line break.
    @exception ParseError Cannot parse this as a number. */
void read_numeric(String& res, SourceText& text)
{
  bool stop;
  res = "";
  char dummy;
  text.LineBreak() = false;

  // To make sure that there is at least one digit:
  bool found_digit = false;

  // Check if there is a sign:
  dummy = text.Current();
  if ( '+' == dummy || '-' == dummy )
    {
      res += dummy;
      text.AdvanceChar();
      if ( text.LineBreak() )
	throw IllegalLinebreak( "Linebreak after sign.",
				text.File(),
				text.Line(),
				text.Column() ); 
    }

  // There could be some digits here:
  stop = false;
  while (!stop) 
    {
      char dummy = text.Current();
      if ( isdigit(dummy) )
	{
	  found_digit = true;
	  res += dummy;
	  text.AdvanceChar();
	  if ( text.LineBreak() ) return; // Line break ends scanning immediately.
	}
      else
	{
	  stop = true;
	}
    }

  // Next there can be a decimal point
  if ( '.' == text.Current() )
    {
      res += ".";
      text.AdvanceChar();
      if ( text.LineBreak() )
	if (found_digit)
	  {
	    // Line break ends scanning immediately, if we have
	    // already found at least one digit.
	    return;
	  }
	else
	  {
	    throw IllegalLinebreak("Expected at least one digit.",
				   text.File(),
				   text.Line(),
				   text.Column());
	  }

      // ... followed by optional more digits
      stop = false;
      while (!stop) 
	{
	  char dummy = text.Current();
	  if ( isdigit(dummy) )
	    {
	      found_digit = true;
	      res += dummy;
	      text.AdvanceChar();
	      if ( text.LineBreak() ) return; // Line break ends scanning immediately.
	    }
	  else
	    {
	      stop = true;
	    }
	}    
    }

  // At this point, we must have found at least one digit.
  if (!found_digit)
    throw ParseError("Expected at least one digit.",
		     text.File(),
		     text.Line(),
		     text.Column());

  // Now there could be a `e' or `E':
  dummy = text.Current();
  if ( 'e' == dummy || 'E' == dummy )
    {
      res += dummy;
      text.AdvanceChar();
      if ( text.LineBreak() )
	throw IllegalLinebreak( "Linebreak after e/E.",
				text.File(),
				text.Line(),
				text.Column() );
      
      // Now there must be an integer (with optional sign)
      {
	String s;
	read_integer(s,text);
	res += s;
      }
    }
}

/** Use a String stream to parse an integer number. */
void parse_integer(Index& n, SourceText& text)
{
  String res;
  read_integer(res, text);
  istringstream is(res);
  is >> n;
}

/** Use a String stream to parse a floating point number. */
void parse_numeric(Numeric& n, SourceText& text)
{
  String res;
  read_numeric(res, text);
  istringstream is(res);
  is >> n;
}

/** Read a vector of Strings. This looks as follows in the control
    file: ["String1","String2"]

    Whitespace has to have been eaten before, that is, the current
    character must be `['.
  
    The empty vector is allowed.
  
    Quotation marks inside Strings are currently not possible.
  
    Line breaks are allowed before and after each String. Line breaks
    inside Strings are not allowed. 
   
    @see parse_String */
void parse_Stringvector(ArrayOfString& res, SourceText& text)
{
  bool first = true;		// To skip the first comma.
  res.resize(0);			// Clear the result vector (just in case).

  // Make sure that the current character really is `[' and proceed.
  assertain_character('[',text);
  // There might have occured a linebreak, which is fine.
  
  eat_whitespace(text);

  // Read the elements of the vector (`]' means that we have
  // reached the end):
  while ( ']' != text.Current() )
    {
      String dummy;

      if (first)
	first = false;
      else
	{
	  assertain_character(',',text);
	  eat_whitespace(text);
	}

      parse_String(dummy, text);
      res.push_back(dummy);
      eat_whitespace(text);
    }

  text.AdvanceChar();
}

/** Read a vector of integers. This looks as follows in the control
    file: [123,5,334]
    Whitespace has to have been eaten before, that is, the current
    character must be `['.
  
    The empty vector is allowed.
    
    Line breaks are allowed before and after each number. Line breaks
    inside numbers are not allowed. 
   
    @see parse_integer */
void parse_intvector(ArrayOfIndex& res, SourceText& text)
{
  bool first = true;		// To skip the first comma.
  res.resize(0);			// Clear the result vector (just in case).

  // Make sure that the current character really is `[' and proceed.
  assertain_character('[',text);
  // There might have occured a linebreak, which is fine.
  
  eat_whitespace(text);

  // Read the elements of the vector (`]' means that we have
  // reached the end):
  while ( ']' != text.Current() )
    {
      Index dummy;

      if (first)
	first = false;
      else
	{
	  assertain_character(',',text);
	  eat_whitespace(text);
	}

      parse_integer(dummy, text);
      res.push_back(dummy);
      eat_whitespace(text);
    }

  text.AdvanceChar();
}

/** Read a vector of Numerics. This looks as follows in the control
    file: [1.3, 5, 3.4]
    Whitespace has to have been eaten before, that is, the current
    character must be `['.
  
    The empty vector is allowed.
  
    Line breaks are allowed before and after each number. Line breaks
    inside numbers are not allowed. 
   
    @see parse_numeric */
void parse_numvector(Vector& res, SourceText& text)
{
  bool first = true;		// To skip the first comma.

  // We need a temporary Array<Numeric>, so that we can use push_back
  // to store the values. FIXME: Need also constructor for Vector from
  // Array<Numeric>.
  Array<Numeric> tres;

  // Make sure that the current character really is `[' and proceed.
  assertain_character('[',text);
  // There might have occured a linebreak, which is fine.
  
  eat_whitespace(text);

  // Read the elements of the vector (`]' means that we have
  // reached the end):
  while ( ']' != text.Current() )
    {
      Numeric dummy;

      if (first)
	first = false;
      else
	{
	  assertain_character(',',text);
	  eat_whitespace(text);
	}

      parse_numeric(dummy, text);
      tres.push_back(dummy);
      eat_whitespace(text);
    }

  // Copy tres to res:
  res.resize(tres.nelem());
  res = tres;

  text.AdvanceChar();
}

/** Parse the Contents of text as ARTS control input. 

    Either values or tasks will be empty.

    @param id     Output. Method id.
    @param values Output. Keyword parameter values for this method.
    @param output Output. Output workspace variables (for generic methods).
    @param input  Output. Input workspace variables (for generic methods).
    @param tasks  Output. A list of other methods.
    @param text The input to parse.

    @param no_eot Suppress throwing an error on EOT after the closing
    curly brace.
   
    @see read_name
    @see eat_whitespace
    @see assertain_character
    @see parse_String
    @see parse_integer
    @see parse_numeric
    @see parse_Stringvector
    @see parse_intvector
    @see parse_numvector
   
   @exception UnknownMethod
   @exception UnknownWsv
   @exception WrongWsvGroup
   @exception UnexpectedKeyword

   @author Stefan Buehler  */
void parse_method(Index& id, 
		  Array<TokVal>& values,
		  ArrayOfIndex& output,
		  ArrayOfIndex& input,
		  Agenda&       tasks,
		  SourceText& text,
		  bool no_eot=false)
{
  extern const Array<WsvRecord> wsv_data;
  extern const Array<MdRecord> md_data;
  extern const Array<MdRecord> md_data_raw;
  extern const ArrayOfString wsv_group_names;

  extern const std::map<String, Index> MdMap;
  extern const std::map<String, Index> MdRawMap;
  extern const std::map<String, Index> WsvMap;

  Index wsvid;			// Workspace variable id, is used to
				// access data in wsv_data.

  String methodname;		// We need this out here, since it is
				// set once and later modified.

  const MdRecord* mdd;		// Handle on the method record. Needed here,
				// because it is modified.

  bool still_supergeneric=true; // Flag that our MdRecord still is
				// from md_data_raw, not from
				// md_data. 

  // Clear all output variables:
  id = 0;
  values.resize( 0 );
  output.resize( 0 );
  input.resize(  0 );
  tasks.resize(  0 );
  
  {
    read_name(methodname, text);

    {
      // Find method raw id in raw map:
      const map<String, Index>::const_iterator i = MdRawMap.find(methodname);
      if ( i == MdRawMap.end() )
	throw UnknownMethod(methodname,
			    text.File(),
			    text.Line(),
			    text.Column());

      id = i->second;

      // Get a convenient handle on the data record for this method. We
      // have to use a pointer here, not a reference, because we later
      // want to change where mdd is pointing!
      mdd = &md_data_raw[id];

//       cout << "id=" << id << '\n';   
//       cout << "Method: " << mdd->Name() << '\n';

      // Is this a supergeneric method? If not, take the record in
      // md_data, rather than in md_data_raw:
      if ( !mdd->Supergeneric() )
	{
	  // Find explicit method id in MdMap:
	  const map<String, Index>::const_iterator i = MdMap.find(methodname);
	  assert ( i != MdMap.end() );
	  id = i->second;	      
	  
	  mdd = &md_data[id];

	  still_supergeneric = false;

// 	  cout << "Adjusted id=" << id << '\n';   
// 	  cout << "Adjusted Method: " << mdd->Name() << '\n';
	}
    }
  }

  eat_whitespace(text);

  // For generic methods the output and input workspace variables have 
  // to be parsed (given in round brackets).
  if ( 0 < mdd->GOutput().nelem() + mdd->GInput().nelem() )
    {
      //      cout << "Generic!" << id << mdd->Name() << '\n';
      String wsvname;
      bool first = true;	// To skip the first comma.

      assertain_character('(',text);
      eat_whitespace(text);

      // First read all output Wsvs:
      for ( Index j=0 ; j<mdd->GOutput().nelem() ; ++j )
	{
	  if (first)
	    first = false;
	  else
	    {
	      assertain_character(',',text);
	      eat_whitespace(text);
	    }

	  read_name(wsvname, text);

	  {
	    // Find Wsv id:
	    const map<String, Index>::const_iterator i = WsvMap.find(wsvname);
	    if ( i == WsvMap.end() )
	      throw UnknownWsv( wsvname,
				text.File(),
				text.Line(),
				text.Column() );

	    wsvid = i->second;
	  }

	  // If this is a supergeneric method, now is the time to find
	  // out the actual group of the argument(s)!
	  if ( still_supergeneric )
	    {
// 	      cout << "wsvid = " << wsvid << "\n";
// 	      cout << "wsv_group_names[wsv_data[wsvid].Group()] = "
// 		   << wsv_group_names[wsv_data[wsvid].Group()] << "\n";
	      ostringstream os;
	      os << mdd->Name() << "_sg_"
		 << wsv_group_names[wsv_data[wsvid].Group()];
	      methodname = os.str();

	      // Find explicit method id in MdMap:
	      const map<String, Index>::const_iterator i = MdMap.find(methodname);
	      assert ( i != MdMap.end() );
	      id = i->second;	      
	      
	      mdd = &md_data[id];

	      still_supergeneric = false;

// 	      cout << "Adjusted id=" << id << '\n';   
// 	      cout << "Adjusted Method: " << mdd->Name() << '\n';
	    }

	  // Now we have explicitly the method record for the right
	  // group. From now on no special treatment of supergeneric
	  // methods should be necessary.

	  // Check that this Wsv belongs to the correct group:
	  if ( wsv_data[wsvid].Group() != mdd->GOutput()[j] )
	    {
       	    throw WrongWsvGroup( wsvname+" is not "+
	                wsv_group_names[mdd->GOutput()[j]]+", it is "+ 
                        wsv_group_names[wsv_data[wsvid].Group()],
				 text.File(),
				 text.Line(),
				 text.Column() );
	    }

	  // Add this one to the list of output workspace variables:
	  output.push_back(wsvid);
	  
	  eat_whitespace(text);
	}

      // Then read all input Wsvs:
      for ( Index j=0 ; j<mdd->GInput().nelem() ; ++j )
	{
	  if (first)
	    first = false;
	  else
	    {
	      assertain_character(',',text);
	      eat_whitespace(text);
	    }

	  read_name(wsvname, text);

	  {
	    // Find Wsv id:
	    const map<String, Index>::const_iterator i = WsvMap.find(wsvname);
 	    if ( i == WsvMap.end() )
 	      throw UnknownWsv( wsvname,
 				text.File(),
 				text.Line(),
 				text.Column() );

	    wsvid = i->second;
	  }

	  // Is the method data record still supergeneric? This could
	  // be the case if there are no output arguments, only input
	  // arguments. In that case, let's find out the actual group!
	  if ( still_supergeneric )
	    {
	      ostringstream os;
	      os << mdd->Name() << "_sg_"
		 << wsv_group_names[wsv_data[wsvid].Group()];
	      methodname = os.str();

	      // Find explicit method id in MdMap:
	      const map<String, Index>::const_iterator i = MdMap.find(methodname);
	      assert ( i != MdMap.end() );
	      id = i->second;	      

	      mdd = &md_data[id];

	      still_supergeneric = false;

// 	      cout << "Adjusted id=" << id << '\n';   
// 	      cout << "Adjusted Method: " << mdd->Name() << '\n';
	    }

	  // Now we have explicitly the method record for the right
	  // group. From now on no special treatment of supergeneric
	  // methods should be necessary.

	  // Check that this Wsv belongs to the correct group:
	  if ( wsv_data[wsvid].Group() != mdd->GInput()[j] )
       	    throw WrongWsvGroup( wsvname+" is not "+
	                wsv_group_names[mdd->GInput()[j]]+", it is "+ 
                        wsv_group_names[wsv_data[wsvid].Group()],
				 text.File(),
				 text.Line(),
				 text.Column() );

	  // Add this one to the list of input workspace variables:
	  input.push_back(wsvid);
	  
	  eat_whitespace(text);
	}

      assertain_character(')',text);
      eat_whitespace(text);
    }

  // Now look for the curly braces:
  assertain_character('{',text);
  eat_whitespace(text);

  // There are two kind of methods, agenda methods, which have other
  // methods in the body, and normal methods, expecting keywords and
  // values. Let's take the agenda case first...
  if ( mdd->AgendaMethod() )
    {
      out3 << "- " << mdd->Name() << "\n";

      out3 << "{\n";

      parse_agenda(tasks,text);

      out3 << "}\n";
    }
  else
    {
      // Now we have to deal with two different cases: Keywords with
      // parameters, or (optionally) only a parameter without a keyword
      // for methods that have only a single argument.
      //
      // We can distinguish the two cases if we check whether the current
      // character is a letter. (If the parameter is specified directly it
      // must be either a number, a +- sign or a quotation mark)
      //
      // KEYWORDS THAT START WITH A NUMBER WILL BREAK THIS CODE!!
      //
      for ( Index i=0 ; i<mdd->Keywords().nelem() ; ++i )
	{
	  if (!isalpha(text.Current()) && 1==mdd->Keywords().nelem())
	    {
	      // Parameter specified directly, without a keyword. This is only
	      // allowed for single parameter methods!

	      // We don't have to do anything here.
	    }
	  else  
	    {      // Look for the keywords and read the parameters:
	  
	      String keyname;
	      read_name(keyname,text);

	      // Is the keyname the expected keyname?
	      if ( keyname != mdd->Keywords()[i] )
		{
		  throw UnexpectedKeyword( keyname,
					   text.File(),
					   text.Line(),
					   text.Column());
		}

	      eat_whitespace(text);

	      // Look for '='
	      assertain_character('=',text);
	      eat_whitespace(text);
	    }

	  // Now parse the key value. This can be:
	  // String_t,    Index_t,    Numeric_t,
	  // Array_String_t, Array_Index_t, Vector_t,
	  switch (mdd->Types()[i]) 
	    {
	    case String_t:
	      {
		String dummy;
		parse_String(dummy, text);
		values.push_back(dummy);
		break;
	      }
	    case Index_t:
	      {
		Index n;
		parse_integer(n, text);
		values.push_back(n);
		break;
	      }
	    case Numeric_t:
	      {
		Numeric n;
		parse_numeric(n, text);
		values.push_back(n);
		break;
	      }
	    case Array_String_t:
	      {
		ArrayOfString dummy;
		parse_Stringvector(dummy, text);
		values.push_back(dummy);
		break;
	      }
	    case Array_Index_t:
	      {
		ArrayOfIndex dummy;
		parse_intvector(dummy, text);
		values.push_back(dummy);
		break;
	      }
	    case Vector_t:
	      {
		Vector dummy;
		parse_numvector(dummy, text);
		values.push_back(dummy);
		break;
	      }
	    default:
	      throw logic_error("Impossible parameter type.");
	      break;
	    }

	  eat_whitespace(text);

	  // Check:
	  //      cout << "Value: " << mdd->Values()[i] << '\n';
	}
    }

  // Now look for the closing curly braces.  We have to catch Eot,
  // because after a method description may be a good place to end
  // the control file.

  try
    {
      assertain_character('}',text);
    }
  catch (const Eot x)
    {
      //      cout << "EOT!!!!" << endl;
      // Re-trow the error if the no_eot flag is not set:
      if (!no_eot) throw Eot(x);
    }
}

/** Parse the Contents of text as ARTS control input. 
  
    This method is used to parse the list of methods given in the
    curly braces of an agenda method. So the end is marked by a
    closing curly brace ahead.

    @param tasklist Output. The ids and keyword parameter values for the
                    methods to run.
  
    @param text   The input to parse.
   
    @see eat_whitespace
    @see parse_method
          
    @author Stefan Buehler */
void parse_agenda( Agenda& tasklist,
		   SourceText& text )
{
  extern const std::map<String, Index> MdMap;
  extern const std::map<String, Index> WsvMap;
  extern const Array<MdRecord> md_data;

  // For method ids:
  Index id;		
  // For keyword parameter values:
  Array<TokVal> values;
  // Output workspace variables (for generic methods):
  ArrayOfIndex output;		
  // Input workspace variables (for generic methods):
  ArrayOfIndex input;
  // For Agenda, if ther is any:
  Agenda tasks;

  eat_whitespace(text);

  while ( '}' != text.Current() )
    {
      parse_method(id,values,output,input,tasks,text);

      // Append taks to task list:      
      tasklist.push_back(MRecord(id,values,output,input,tasks));

      {
	// Everything in this block is just to generate some
	// informative output.  
	extern const Array<WsvRecord> wsv_data;

	out3 << "- " << md_data[id].Name() << "\n";

	for ( Index j=0 ; j<values.nelem() ; ++j )
	  {
	    out3 << "   " 
		 << md_data[id].Keywords()[j] << ": "
		 << values[j] << '\n';
	  }
	  
	// Output workspace variables for generic methods:
	if ( 0 < md_data[id].GOutput().nelem() + md_data[id].GInput().nelem() )
	  {
	    out3 << "   Output: ";
	    for ( Index j=0 ; j<output.nelem() ; ++j )
	      {
		out3 << wsv_data[output[j]].Name() << " ";
	      }
	    out3 << "\n";

	    out3 << "   Input: ";
	    for ( Index j=0 ; j<input.nelem() ; ++j )
	      {
		out3 << wsv_data[input[j]].Name() << " ";
	      }
	    out3 << "\n";
	  }
      }
      
      eat_whitespace(text);
    }
}

/** The main function of the parser. This will parse the entire
    text. FIXME: Add more documentation here.

    @param tasklist Output. The method runtime data that can be used to
                    execute the methods.
    @param text The control text

    @author Stefan Buehler */
void parse_main(Agenda& tasklist, SourceText& text)
{
  extern const std::map<String, Index> MdMap;
  extern const std::map<String, Index> WsvMap;

  try 
    {
      extern const Array<MdRecord> md_data;
      extern const Array<WsvRecord> wsv_data;

      // For method ids:
      Index id;		
      // For keyword parameter values:
      Array<TokVal> values;
      // Output workspace variables (for generic methods):
      ArrayOfIndex output;		
      // Input workspace variables (for generic methods):
      ArrayOfIndex input;
      // For Agenda, if ther is any:
      Agenda tasks;

      out3 << "\nParsing control text:\n";

      text.Init();
      eat_whitespace(text);

      parse_method(id,values,output,input,tasklist,text,true);
	  
      if ( "Main" != md_data[id].Name() )
	{
	  out0 << "The outermost method must be Main!\n"
	       << "(But it seems to be " << md_data[id].Name() << ".)\n";
	  exit(1);
	}
    }
  catch (const Eot x)
    {
      // Unexpected end of the source text:
      out0 << "Unexpected end of control script.\n";
      out0 << "File: " << x.file() << '\n';
      out0 << "Line: " << x.line() << '\n';
      exit(true);
    }
  catch (const UnexpectedChar x)
    {
      // Unexpected Character:
      out0 << "Unexpected character:\n";
      out0 << x.what()   << '\n';
      out0 << "File: "   << x.file() << '\n';
      out0 << "Line: "   << x.line() << '\n';
      out0 << "Column: " << x.column() << '\n';
      exit(true);
    }
  catch (const IllegalLinebreak x)
    {
      // A line break in an illegal position:
      out0 << "Illegal Line break:\n";
      out0 << x.what()   << '\n';
      out0 << "File: "   << x.file() << '\n';
      out0 << "Line: "   << x.line() << '\n';
      exit(true);
    }
  catch (const UnknownMethod x)
    {
      // Method unknown:
      // [**This should give a hint on how to obtain a list of allowed 
      // methods.]
      out0 << "Unknown Method:\n";
      out0 << x.what()   << '\n';
      out0 << "File: "   << x.file() << '\n';
      out0 << "Line: "   << x.line() << '\n';
      out3 << "Column: " << x.column() << '\n';
      exit(true);
    }
  catch (const UnknownWsv x)
    {
      // Workspace variable unknown:
      // [**This should give a hint on how to obtain a list of allowed 
      // Wsvs.]
      out0 << "Unknown workspace variable:\n";
      out0 << x.what()   << '\n';
      out0 << "File: "   << x.file() << '\n';
      out0 << "Line: "   << x.line() << '\n';
      out3 << "Column: " << x.column() << '\n';
      exit(true);
    }
  catch (const WrongWsvGroup x)
    {
      // Workspace variable unknown:
      // [**This should give a hint on how to obtain a list of Wsvs in 
      // this group.
      out0 << "Workspace variable belongs to the wrong group:\n";
      out0 << x.what()   << '\n';
      out0 << "File: "   << x.file() << '\n';
      out0 << "Line: "   << x.line() << '\n';
      out3 << "Column: " << x.column() << '\n';
      exit(true);
    }
  catch (const UnexpectedKeyword x)
    {
      // Keyword unknown:
      // [**This should give a hint on how to obtain a list of allowed 
      // keywords.]
      out0 << "Unknown keyword:\n";
      out0 << x.what()   << '\n';
      out0 << "File: "   << x.file() << '\n';
      out0 << "Line: "   << x.line() << '\n';
      out3 << "Column: " << x.column() << '\n';
      exit(true);
    }
  catch (const ParseError x)
    {
      // General Parse Error (parent of all the above):
      out0 << "Parse error:\n";
      out0 << x.what()   << '\n';
      out0 << "File: "   << x.file() << '\n';
      out0 << "Line: "   << x.line() << '\n';
      out0 << "Column: " << x.column() << '\n';
      exit(true);
    }
  catch (const runtime_error x)
    {
      cout << "Runtime error: ";
      cout << x.what() << '\n';
    }
  catch (const logic_error x)
    {
      cout << "Logic error: ";
      cout << x.what() << '\n';
    }

}
