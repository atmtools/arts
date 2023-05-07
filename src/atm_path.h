#pragma once

#include "atm.h"
#include "jacobian.h"
#include "matpack_arrays.h"
#include "ppath_struct.h"
#include "propagationmatrix.h"
#include "transmissionmatrix.h"
#include <tuple>

ArrayOfAtmPoint& atm_path_resize(ArrayOfAtmPoint&, const Ppath&);

void forward_atm_path(ArrayOfAtmPoint&, const Ppath&, const AtmField&);

ArrayOfAtmPoint forward_atm_path(const Ppath&, const AtmField&);

ArrayOfVector& path_freq_resize(ArrayOfVector&, const Vector&, const ArrayOfAtmPoint&);

void forward_path_freq(ArrayOfVector&, const Vector&, const Ppath&, const ArrayOfAtmPoint&, const Numeric);

ArrayOfVector forward_path_freq(const Vector&, const Ppath&, const ArrayOfAtmPoint&, const Numeric);

AtmField forward_1d_atm_field(const ArrayOfAtmPoint&, const Ppath&);
