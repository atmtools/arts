#include "mystring.h"

#include <nonstd.h>

#include <algorithm>
#include <format>
#include <ranges>
#include <string>

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

String comma(bool& first, const String& spaces) {
  if (first) {
    first = false;
    return "";
  }
  return std::format("{}{}{}", ',', (spaces.size() ? '\n' : ' '), spaces);
}

void join(String& res, const ArrayOfString& list, const String& with) {
  res.clear();

  if (list.empty()) return;

  for (auto& v : list | std::views::drop(1)) {
    std::format_to(std::back_inserter(res), "{}{}", v, with);
  }

  std::format_to(std::back_inserter(res), "{}", list.back());
}

String join(const ArrayOfString& list, const String& with) {
  String res;
  join(res, list, with);
  return res;
}

void replace(String& x, const String& from, const String& to) {
  const auto arr = split(x, from);
  join(x, arr, to);
}
