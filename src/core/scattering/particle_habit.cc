#include "particle_habit.h"
#include "matpack.h"
#include "util/sorting.h"


namespace scattering {

std::pair<Numeric, Numeric> derive_scat_species_a_and_b(const Vector& sizes,
                                                        const Vector& masses,
                                                        const Numeric& fit_start,
                                                        const Numeric& fit_end) {
  const Index nse = sizes.size();
  ARTS_ASSERT(nse > 1);

  ArrayOfIndex intarr_sort, intarr_unsort(0);
  Vector x_unsorted(nse), m_unsorted(nse);
  Vector q;
  Index nsev = 0;

  for (Index i = 0; i < nse; i++) {
    if (std::isnan(sizes[i]))
      ARTS_USER_ERROR("NaN found in selected size grid data.");
    if (std::isnan(masses[i]))
      ARTS_USER_ERROR("NaN found among particle mass data.");

    if (sizes[i] >= fit_start && sizes[i] <= fit_end) {
      x_unsorted[nsev] = sizes[i];
      m_unsorted[nsev] = masses[i];
      nsev += 1;
    }
  }

  if (nsev < 2)
    ARTS_USER_ERROR("Less than two size points found in the range "
                    "[fit_start, fit_end]. It is then not possible "
                    "to determine the a and b parameters.");

  get_sorted_indexes(intarr_sort, x_unsorted[Range(0, nsev)]);
  Vector log_x(nsev), log_m(nsev);

  for (Index i = 0; i < nsev; i++) {
    log_x[i] = log(x_unsorted[intarr_sort[i]]);
    log_m[i] = log(m_unsorted[intarr_sort[i]]);
  }

  linreg(q, log_x, log_m);
  return std::pair<Numeric, Numeric>(exp(q[0]), q[1]);
}


}
