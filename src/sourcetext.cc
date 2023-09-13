#include "sourcetext.h"
#include <iostream>
#include "file.h"

void SourceText::AppendFile(const String& name) {
  mSfLine.push_back(mText.size());
  mSfName.push_back(name);

  mText = read_text_from_file(name);
}

void SourceText::AdvanceChar() {
  if (mColumn < mText[mLine].size() - 1) {
    ++mColumn;
  } else {
    mLineBreak = true;
    do {
      if (mLine >= mText.size()) {
        throw Eot("", this->File(), this->Line(), this->Column());
      } else if (mLine == mText.size() - 1) {
        mColumn++;
        break;
      } else {
        ++mLine;
        mColumn = 0;
      }
    } while (1 > mText[mLine].size());  // Skip empty lines.
  }
}

void SourceText::AdvanceLine() {
  mLineBreak = true;
  mColumn = 0;
  do {
    if (mLine >= mText.size() - 1) {
      throw Eot("", this->File(), this->Line(), this->Column());
    } else {
      ++mLine;
    }
  } while (1 > mText[mLine].size());  // Skip empty lines.
}

const String& SourceText::File() {
  Index i = 0;
  bool stop = false;

  while (i < mSfLine.size() - 1 && !stop) {
    if (mLine >= mSfLine[i + 1])
      ++i;
    else
      stop = true;
  }

  return mSfName[i];
}

void SourceText::Init() {
  mLine = 0;
  mColumn = 0;

  if (1 > mText.size()) {
    throw Eot("Empty text!", this->File(), this->Line(), this->Column());
  } else {
    // Skip empty lines:
    while (1 > mText[mLine].size()) {
      if (mLine >= mText.size() - 1) {
        throw Eot("", this->File(), this->Line(), this->Column());
      } else {
        mLineBreak = true;
        ++mLine;
      }
    }
  }
}

Index SourceText::GetSourceLine(const Index line) {
  Index i = 0;
  bool stop = false;

  while (i < mSfLine.size() - 1 && !stop) {
    if (line >= mSfLine[i + 1])
      ++i;
    else
      stop = true;
  }

  return line - mSfLine[i] + 1;
}

std::ostream& operator<<(std::ostream& os, const SourceText& text) {
  for (Index i = 0; i < text.mText.size(); ++i)
    os << i << "(" << text.mText[i].size() << ")"
       << ": " << text.mText[i] << '\n';
  return (os);
}
