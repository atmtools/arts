#pragma once

#include <cstddef>
#include <cstdlib>

#include "debug.h"
#include "mystring.h"

//----------< First of all, include the configuration header >----------
#include "config.h"

#ifdef HAVE_NAMESPACES
// We need those to support ansi-compliant compilers (gcc-3x)
namespace std {}
using namespace std;
#endif

