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

#include "wigner_functions.h"
#include "messages.h"

void Wigner6Init(Index& wigner_initialized,
                 #if DO_FAST_WIGNER
                 const Index& fast_wigner_stored_symbols,
                 #else
                 const Index&,
                 #endif
                 const Index& largest_wigner_symbol_parameter,
                 const Verbosity&)
{
  wigner_initialized = largest_wigner_symbol_parameter;
  
  #if DO_FAST_WIGNER
    fastwigxj_load(FAST_WIGNER_PATH_3J, 3, NULL);
    fastwigxj_load(FAST_WIGNER_PATH_6J, 6, NULL);
    #ifdef _OPENMP
      fastwigxj_thread_dyn_init(3, fast_wigner_stored_symbols);
      fastwigxj_thread_dyn_init(6, fast_wigner_stored_symbols);
    #else
      fastwigxj_dyn_init(3, fast_wigner_stored_symbols);
      fastwigxj_dyn_init(6, fast_wigner_stored_symbols);
    #endif
  #endif
  wig_table_init(int(largest_wigner_symbol_parameter*2), 6);
}

void Wigner3Init(Index& wigner_initialized,
                 #if DO_FAST_WIGNER
                 const Index& fast_wigner_stored_symbols,
                 #else
                 const Index&,
                 #endif
                 const Index& largest_wigner_symbol_parameter,
                 const Verbosity&)
{
  wigner_initialized = largest_wigner_symbol_parameter;
  
  #if DO_FAST_WIGNER
  fastwigxj_load(FAST_WIGNER_PATH_3J, 3, NULL);
  #ifdef _OPENMP
  fastwigxj_thread_dyn_init(3, fast_wigner_stored_symbols);
  #else
  fastwigxj_dyn_init(3, fast_wigner_stored_symbols);
  #endif
  #endif
  wig_table_init(int(largest_wigner_symbol_parameter*2), 3);
}

void WignerFastInfoPrint(const Index& wigner_initialized,
                         const Verbosity&)
{
  if(not wigner_initialized)
    throw std::runtime_error("Must first initialize wigner...");
  
  #if DO_FAST_WIGNER
    fastwigxj_print_stats();
  #else
    throw std::runtime_error("You cannot do this without having compiled with fast wigner.");
  #endif
}

void Wigner6Unload(Index& wigner_initialized,
                   const Verbosity&)
{
  if(not wigner_initialized)
    throw std::runtime_error("Must first initialize wigner...");
  wigner_initialized = 0;
  
  #if DO_FAST_WIGNER
    fastwigxj_unload(3);
    fastwigxj_unload(6);
  #endif
  wig_table_free();
}

void Wigner3Unload(Index& wigner_initialized,
                   const Verbosity&)
{
  if(not wigner_initialized)
    throw std::runtime_error("Must first initialize wigner...");
  wigner_initialized = 0;
  
  #if DO_FAST_WIGNER
  fastwigxj_unload(3);
  #endif
  wig_table_free();
}
