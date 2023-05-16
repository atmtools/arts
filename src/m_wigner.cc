/**
 * @file m_wigner.cc
 * @author Richard Larsson
 * @date 2018-04-03
 * 
 * @brief Wigner symbol interactions
 */

#include "messages.h"
#include "wigner_functions.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void Wigner6Init(Index& wigner_initialized,
                 const Index& fast_wigner_stored_symbols,
                 const Index& largest_wigner_symbol_parameter,
                 const Verbosity&) {
  wigner_initialized = make_wigner_ready(int(largest_wigner_symbol_parameter), int(fast_wigner_stored_symbols), 6);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Wigner3Init(Index& wigner_initialized,
                 const Index& fast_wigner_stored_symbols,
                 const Index& largest_wigner_symbol_parameter,
                 const Verbosity&) {
  wigner_initialized = make_wigner_ready(int(largest_wigner_symbol_parameter), int(fast_wigner_stored_symbols), 3);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void WignerFastInfoPrint(const Index& wigner_initialized, const Verbosity&) {
  ARTS_USER_ERROR_IF (not wigner_initialized,
                      "Must first initialize wigner...");

#if DO_FAST_WIGNER
  fastwigxj_print_stats();
#else
  ARTS_USER_ERROR (
      "You cannot do this without having compiled with fast wigner.");
#endif
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Wigner6Unload(Index& wigner_initialized, const Verbosity&) {
  ARTS_USER_ERROR_IF (not wigner_initialized,
                      "Must first initialize wigner...");
  wigner_initialized = 0;

#if DO_FAST_WIGNER
  fastwigxj_unload(3);
  fastwigxj_unload(6);
#endif
  wig_table_free();
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Wigner3Unload(Index& wigner_initialized, const Verbosity&) {
  ARTS_USER_ERROR_IF (not wigner_initialized,
                      "Must first initialize wigner...");
  wigner_initialized = 0;

#if DO_FAST_WIGNER
  fastwigxj_unload(3);
#endif
  wig_table_free();
}
