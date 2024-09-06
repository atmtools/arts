#include "debug.h"

std::string artsformat(std::format_string<> fmt) {
  return std::format(fmt);
}

std::string artsformat() { return ""; }
