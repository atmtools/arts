#pragma once

#include <configtypes.h>

namespace Body {
namespace Earth {
inline constexpr Numeric a  = 6378137.0;
inline constexpr Numeric b  = 6356752.314245;
inline constexpr Numeric GM = 398600.435507;
}  // namespace Earth

namespace Jupiter {
inline constexpr Numeric a  = 71492e3;
inline constexpr Numeric b  = 66854e3;
inline constexpr Numeric GM = 126712764.1;
}  // namespace Jupiter

namespace Mars {
inline constexpr Numeric a  = 3396.19e3;
inline constexpr Numeric b  = 3376.20e3;
inline constexpr Numeric GM = 42828.375816;
}  // namespace Mars

namespace Moon {
inline constexpr Numeric a  = 1738.1e3;
inline constexpr Numeric b  = 1736.0e3;
inline constexpr Numeric GM = 4902.800118;
}  // namespace Moon

namespace Mercury {
inline constexpr Numeric a  = 2439.7e3;
inline constexpr Numeric b  = 2439.7e3;
inline constexpr Numeric GM = 22031.868551;
}  // namespace Mercury

namespace Venus {
inline constexpr Numeric a  = 6051.8e3;
inline constexpr Numeric b  = 6051.8e3;
inline constexpr Numeric GM = 324858.592;
}  // namespace Venus

namespace Saturn {
inline constexpr Numeric a  = 0.5 * 120536e3;
inline constexpr Numeric b  = 0.5 * 108728e3;
inline constexpr Numeric GM = 37940584.841800;
}  // namespace Saturn
}  // namespace Body
