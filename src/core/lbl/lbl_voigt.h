#pragma once

#include "lbl_data.h"
#include "matpack_mdspan_helpers_grid_t.h"

bool is_voigt(LineByLineLineshape lsm);

namespace lbl {
/** Get the Voigt profile for a line assuming no Zeeman effect
 * 
 * @param y Output profile, overwritten
 * @param l The line
 * @param f The frequency grid
 * @param atm The atmospheric point
 * @param mass The mass of the molecule (in u)
 */
void compute_voigt(VectorView y,
                   const line& l,
                   const AscendingGrid& f,
                   const AtmPoint& atm,
                   const Numeric mass);
}  // namespace lbl
