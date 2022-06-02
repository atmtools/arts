#include <mutex>

#include <energylevelmap.h>
#include <jacobian.h>
#include <propagationmatrix.h>

#include "enums.h"
#include "gui.h"

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
  Vector rtp_mag;
  Vector rtp_los;
  Numeric rtp_pressure;
  Numeric rtp_temperature;
  EnergyLevelMap rtp_nlte;
  Vector rtp_vmr;
};

struct Control {
  std::mutex copy;
  std::string error{};
  std::atomic_int pos{0};
  std::atomic_bool run{false}, exit{false};

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

ENUMCLASS(YScaling, char, None, Normalize, CrossSection)
struct DisplayOptions {
  XScaling xscale{XScaling::Hz};
  
  Index jacobian_target{-1};

  YScaling yscale{YScaling::None};
  Numeric yscale_const{1.0};
  bool inverse_yscale{false};

  int smooth_counter{1};
};
}  // namespace PropmatClearsky

void propmat(PropmatClearsky::ResultsArray& res,
             PropmatClearsky::Control& ctrl,
             ArrayOfRetrievalQuantity& jacobian_quantities,
             ArrayOfSpeciesTag& select_abs_species,
             Vector& f_grid,
             Vector& rtp_mag,
             Vector& rtp_los,
             Numeric& rtp_pressure,
             Numeric& rtp_temperature,
             EnergyLevelMap& rtp_nlte,
             Vector& rtp_vmr,
             const ArrayOfArrayOfSpeciesTag&& abs_species);
}  // namespace ARTSGUI