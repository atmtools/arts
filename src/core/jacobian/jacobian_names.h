#include <atm.h>
#include <lbl.h>
#include <obsel.h>
#include <subsurf.h>
#include <surf.h>

#include <functional>

struct ErrorKey;

namespace Jacobian {
struct AtmTarget;
struct SurfaceTarget;
struct SubsurfaceTarget;
struct LineTarget;
struct SensorTarget;
struct ErrorTarget;

using cv          = ConstVectorView;
using cm          = ConstMatrixView;
using atm_vec     = std::function<Vector(cv, const AtmField&)>;
using atm_mat     = std::function<Matrix(cm, cv, const AtmField&)>;
using surf_vec    = std::function<Vector(cv, const SurfaceField&)>;
using surf_mat    = std::function<Matrix(cm, cv, const SurfaceField&)>;
using subsurf_vec = std::function<Vector(cv, const SubsurfaceField&)>;
using subsurf_mat = std::function<Matrix(cm, cv, const SubsurfaceField&)>;
using line_vec    = std::function<Vector(cv, const AbsorptionBands&)>;
using line_mat    = std::function<Matrix(cm, cv, const AbsorptionBands&)>;
using sensor_vec  = std::function<Vector(cv, const ArrayOfSensorObsel&)>;
using sensor_mat  = std::function<Matrix(cm, cv, const ArrayOfSensorObsel&)>;
using error_vec   = std::function<Vector(cv, cv)>;
using error_mat   = std::function<Matrix(cm, cv, cv)>;
}  // namespace Jacobian
