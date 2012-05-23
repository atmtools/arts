include (CheckCCompilerFlag)
include (CheckCXXCompilerFlag)

macro (ARTS_ADD_COMPILER_FLAG FLAG)
  Check_C_Compiler_Flag(-${FLAG} CCFLAG_${FLAG})
  if (CCFLAG_${FLAG})
    set (ARTS_C_FLAGS "${ARTS_C_FLAGS} -${FLAG} ${ARGN}")
  endif (CCFLAG_${FLAG})
  Check_CXX_Compiler_Flag(-${FLAG} CXXFLAG_${FLAG})
  if (CXXFLAG_${FLAG})
    set (ARTS_CXX_FLAGS "${ARTS_CXX_FLAGS} -${FLAG} ${ARGN}")
  endif (CXXFLAG_${FLAG})
endmacro (ARTS_ADD_COMPILER_FLAG FLAG)

