/* Copyright (C) 2012 Oliver Lemke <olemke@core-dump.info>

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

#ifndef sourcetext_h
#define sourcetext_h

#include "array.h"
#include "mystring.h"
#include "exceptions.h"


/** A smart class to hold the text for parsing. A variable of this
    class can hold not only the text of the ARTS Control file, but
    also a position in the text. This is handy for parsing. There is
    also a function to return the current character and functions to
    advance the position. (AdvanceChar advances the position to the
    next character, doing a line break if necessary; AdvanceLine goes
    to the next line.)
    
    mLine and mColumn are 0 based, but Line and Column are 1 based.
    
    @author Stefan Buehler */
class SourceText {
public:

  /** Default constructor. */
  SourceText() :  mLine(0), mColumn(0), mLineBreak(false)
  { /* Nothing to be done here. */ }

  /** Appends contents of file to the source text.
      @see read_text_from_file */
  void AppendFile(const String& name);

  /** Return the current character. */
  char Current() {
    if (reachedEot())
      throw Eot( "", this->File(), this->Line(), this->Column() ); 

    return mText[mLine][mColumn];
  }

  /** Check if the current position reached the end. */
  bool reachedEot() {
    return (mLine >= mText.nelem()
            || (mLine == mText.nelem()-1 && mColumn >= mText[mLine].nelem()));
  }
  
  /** Advance position pointer by one character. Sets mLineBreak if a
      line break occured.  
   
      @exception Eot The end of text is reached. */
  void AdvanceChar();

  /** Advances position pointer by one line.

      @exception Eot The end of the text is reached. */
  void AdvanceLine();

  /** Return the filename associated with the current position. */
  const String& File();

  /** Return the line number, but for the file that is associated
      with the current position. */
  Index Line() { return GetSourceLine(mLine); };

  /** Return the marked line number, but for the file that is associated
      with the current position. */
  Index MarkedLine() { return GetSourceLine(mMarkedLine); };

  /** Return the line index. */
  Index LineRaw() { return mLine; }

  /** Return the column index. */
  Index ColumnRaw() { return mColumn; }

  /** Return the current column. */
  Index Column() { return mColumn+1; }

  /** Return the current marked column. */
  Index MarkedColumn() { return mMarkedColumn+1; }

  /** Set current position. */
  void SetPosition(Index line, Index column)
    {
        mLine = line;
        mColumn = column;
    }

  /** Mark current position. */
  void SetMark()
    {
        mMarkedLine = mLine;
        mMarkedColumn = mColumn;
    }
    
  /** This sets the pointer to the first existing character in the
      text. (First few lines could be empty). */
  void Init();

  /** Read the line break flag. Set this to false before an operation
      that you want to monitor and check it afterwards. */
  bool& LineBreak() { return mLineBreak; }

  /** Const version of LineBreak
      @see LineBreak */
  bool  LineBreak() const { return mLineBreak; }

  /** Output operator for SourceText. (Only used for debugging) */
  friend std::ostream& operator << (std::ostream& os, const SourceText& text);

private:

  /** Return the line number, but for the file that is associated
      with the given position. */
  Index GetSourceLine(const Index line);

  /** The text. */
  ArrayOfString mText;

  /** Line position in the text. (0 based!) */
  Index mLine;

  /** Column position in the text. (0 based!) */
  Index mColumn;

  /** Marked line position in the text. (0 based!) */
  Index mMarkedLine;

  /** Marked column position in the text. (0 based!) */
  Index mMarkedColumn;

  /** Remember where which source file starts. */
  ArrayOfIndex mSfLine;

  /** Names associated with @see mSfLine. */
  ArrayOfString mSfName;

  /** Is set to true if the last operation caused a line
      break. Has to be cleared explicitly! */
  bool mLineBreak;
};

#endif /* sourcetext_h */

