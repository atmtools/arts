#include <mutex>

#include <energylevelmap.h>
#include <jacobian.h>
#include <propagationmatrix.h>

#include "atm.h"
#include "enums.h"
#include "gui.h"
#include "transmissionmatrix.h"

namespace ARTSGUI {
namespace PropmatClearsky {
struct ComputeValues {
  PropagationMatrix pm;
  ArrayOfPropagationMatrix aopm;
  StokesVector sv;
  ArrayOfStokesVector aosv;

  ArrayOfRetrievalQuantity jacobian_quantities;
  ArrayOfSpeciesTag select_abs_species;
  Vector f_grid;
  Vector rtp_los;
  AtmPoint atm_point;

  Numeric transmission_distance;
  TransmissionMatrix tm;
  ArrayOfTransmissionMatrix aotm;
};

struct Control {
  std::mutex copy;
  std::string errmsg{};
  std::atomic_int pos{0};
  std::atomic_bool run{false}, exit{false}, error{false};

  static_assert(
      std::atomic_int::is_always_lock_free,
      "Can only compile with GUI if lock-free integers are supported");
  static_assert(std::atomic_bool::is_always_lock_free,
                "Can only compile with GUI if lock-free bools are supported");
};

constexpr Index n = 3;

struct Results {
  bool auto_update{false}, auto_f_grid{false};
  std::atomic_bool ok{false};
  ComputeValues value{};
};

using ResultsArray = std::array<Results, n>;

ENUMCLASS(XScaling, char,
  Hz,
  GHz,
  THz,
  Angcm,
  Kaycm,
  m,
  nm,
  Angfreq)

ENUMCLASS(PropmatScaling, char, None, Normalize, CrossSection)

ENUMCLASS(TramatScaling, char, None, dB)

struct DisplayOptions {
  XScaling xscale{XScaling::Hz};
  
  Index jacobian_target{-1};

  PropmatScaling propmat_scale{PropmatScaling::None};
  Numeric propmat_scale_const{1.0};
  bool inverse_propmat_scale{false};

  TramatScaling tramat_scale{TramatScaling::None};
  bool inverse_tramat_scale{false};

  int smooth_counter{1};

  bool transmission{false};
};
}  // namespace PropmatClearsky

void propmat(PropmatClearsky::ResultsArray& res,
             PropmatClearsky::Control& ctrl,
             ArrayOfRetrievalQuantity& jacobian_quantities,
             ArrayOfSpeciesTag& select_abs_species,
             Vector& f_grid,
             Vector& rtp_los,
             AtmPoint& atm_point,
             Numeric& transmission_distance,
             const ArrayOfArrayOfSpeciesTag&& abs_species);
}  // namespace ARTSGUI