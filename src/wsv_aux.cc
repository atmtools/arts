/* Copyright (C) 2022 Oliver Lemke <oliver.lemke@uni-hamburg.de>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*!
  \file   wsv_aux.cc
  \author Oliver Lemke <oliver.lemke@uni-hamburg.de>
  \date   Thu Jul 14 14:43:37 CEST 2022
  
  \brief  Implementation of WSV aux functions.
*/

#include "wsv_aux.h"

bool WsvRecord::has_defaults() const {
  return not std::holds_alternative<std::unique_ptr<Any>>(defval.value);
}

std::shared_ptr<void> WsvRecord::get_copy() const {
  if (has_defaults())
    return std::visit(
        [](auto&& val) -> std::shared_ptr<void> {
          using value_type =
              std::remove_cv_t<std::remove_pointer_t<decltype(val.get())>>;
          return std::make_shared<value_type>(*val);
        },
        defval.value);
  return nullptr;
}
