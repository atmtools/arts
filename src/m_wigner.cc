/* Copyright (C) 2018
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

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
  ARTS_USER_ERROR_IF (true,
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
