#include "mystring.h"

#include <nonstd.h>

#include <algorithm>
#include <format>
#include <string>

namespace stdr = std::ranges;

void tolower(String& x) {
  stdr::transform(x, x.begin(), [](unsigned char c) { return ::tolower(c); });
}

void toupper(String& x) {
  stdr::transform(x, x.begin(), [](unsigned char c) { return ::toupper(c); });
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

String comma(bool& first, const String& spaces) {
  if (first) {
    first = false;
    return "";
  }
  return std::format("{}{}{}", ',', (spaces.size() ? '\n' : ' '), spaces);
}

void join(String& res,
          const std::span<const String>& list,
          const String& with) {
  res.clear();

  if (list.empty()) return;

  res = list.front();
  for (auto& v : list.subspan(1)) res += std::format("{}{}", with, v);
}

String join(const std::span<const String>& list, const String& with) {
  String res;
  join(res, list, with);
  return res;
}

void replace(String& x, const String& from, const String& to) {
  const auto arr = split(x, from);
  join(x, arr, to);
}
