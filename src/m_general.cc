/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   m_general.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2002-05-08 

  \brief  Workspace functions of a general and overall character.

  This file is for general functions that do not fit in any other "m_"-file.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <chrono>
#include <cstdlib>
#include <ratio>
#include <stdexcept>

#include "array.h"
#include "arts_constants.h"
#include "check_input.h"
#include "m_general.h"
#include "mystring.h"

#include "math_funcs.h"
#include <workspace.h>

#include "sensor.h"

#include "fastem.h"
#include "tessem.h"
#include "arts_omp.h"

inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void INCLUDE() {}

/* Workspace method: Doxygen documentation will be auto-generated */
void PrintWorkspace(const Workspace& ws) {
  std::cout << ws << '\n';
}

/* Workspace method: Doxygen documentation will be auto-generated */
void StringJoin(String& out,
                const String& in1,
                const String& in2,
                const String& in3,
                const String& in4,
                const String& in5,
                const String& in6,
                const String& in7,
                const String& in8,
                const String& in9,
                const String& in10) {
  out = in1 + in2 + in3 + in4 + in5 + in6 + in7 + in8 + in9 + in10;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Error(const String& msg) {
  ARTS_USER_ERROR(msg);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Exit() {
  std::exit(EXIT_SUCCESS);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GetEnvironmentVariable(  // WS Generic Output:
    String& str,
    // WS Generic Input:
    const String& envvar) {
  char* cstr;
  cstr = std::getenv(envvar.c_str());
  ARTS_USER_ERROR_IF (cstr == NULL,
    "Environment variable ", envvar, " does not exist.")
  str = cstr != NULL ? String(cstr) : "";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GetNumberOfThreads(Index& nthreads) {
  nthreads = arts_omp_get_max_threads();
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GetEnvironmentVariable(  // WS Generic Output:
    Index& i,
    // WS Generic Input:
    const String& envvar) {
  char* cstr;
  cstr = std::getenv(envvar.c_str());
  ARTS_USER_ERROR_IF (cstr == NULL || std::strlen(cstr) == 0,
    "Environment variable ", envvar, " "
    "is empty or does not exist.")
  std::istringstream is(cstr);
  is >> i;
  ARTS_USER_ERROR_IF (!is.eof(),
      "Cannot convert environment variable ", envvar, " "
      "to Index: ", cstr)
}

/* Workspace method: Doxygen documentation will be auto-generated */
#ifdef _OPENMP
void SetNumberOfThreads(const Index& nthreads) {
  omp_set_num_threads((int)nthreads);
}
#else
void SetNumberOfThreads(const Index&) {
}
#endif
