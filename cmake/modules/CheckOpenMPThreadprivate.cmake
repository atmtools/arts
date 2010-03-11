# - Check for OpenMP with threadprivate support
#
# Copyright (c) 2010, Oliver Lemke, <olemke@core-dump.info>

include (CheckCSourceCompiles)
include (FindOpenMP)

if (OPENMP_FOUND)

# sample openmp source code to test
set (THREADPRIVATE_C_TEST_SOURCE 
"
#include <omp.h>
int tpvar = 0;
#pragma omp threadprivate(tpvar)

int main() { 
  return 0; 
}
")

set (CMAKE_REQUIRED_FLAGS "${OpenMP_C_FLAGS}")
check_c_source_compiles ("${THREADPRIVATE_C_TEST_SOURCE}" THREADPRIVATE_DETECTED)

if (NOT THREADPRIVATE_DETECTED)
  message (STATUS "Threadprivate unsupported: Disabling OpenMP")
  set (OPENMP_FOUND false)
endif (NOT THREADPRIVATE_DETECTED)

mark_as_advanced(
  OPENMP_FOUND
)

endif (OPENMP_FOUND)

