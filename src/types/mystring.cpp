#include "mystring.h"

namespace std {
std::ostream& operator<<(std::ostream& os, const ArrayOfString& x) {
  for (auto& a : x) os << a << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfString& x) {
  for (auto& a : x) os << a << '\n';
  return os;
}
}  // namespace std

void tolower(String& x) {
  std::transform(x.begin(), x.end(), x.begin(), [](unsigned char c) {
    return ::tolower(c);
  });
}

String tolower(const String& x) {
  String out = x;
  tolower(out);
  return out;
}

void toupper(String& x) {
  std::transform(x.begin(), x.end(), x.begin(), [](unsigned char c) {
    return ::toupper(c);
  });
}

String toupper(const String& x) {
  String out = x;
  toupper(out);
  return out;
}

void split(ArrayOfString& aos, const String& x, const String& delim) {
  Size pos, oldpos;
  pos = oldpos = 0;
  aos.resize(0);

  while (oldpos < x.size() && (pos = x.find(delim, oldpos)) != x.npos) {
    if (pos && pos - oldpos) aos.push_back(x.substr(oldpos, pos - oldpos));
    oldpos = pos + delim.size();
  }

  if (oldpos < x.size()) aos.push_back(x.substr(oldpos));
}

ArrayOfString split(const String& x, const String& delim) {
  ArrayOfString out;
  split(out, x, delim);
  return out;
}

void trim(String& x) {
  while (nonstd::isspace(x.front())) x.erase(x.begin());
  while (nonstd::isspace(x.back())) x.pop_back();
}

String trim(const String& x) {
  String out = x;
  trim(out);
  return out;
}
