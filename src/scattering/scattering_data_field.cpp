#include "scattering/scattering_data_field.h"

namespace scattering {

std::ostream &operator<<(std::ostream &output,
                         const scattering::DataFormat format) {
    switch (format) {
    case scattering::DataFormat::Gridded:
        output << "gridded";
        break;
    case scattering::DataFormat::Spectral:
        output << "spectral";
        break;
    case scattering::DataFormat::FullySpectral:
        output << "fully spectral";
    }
    return output;
}

std::ostream &operator<<(std::ostream &output,
                         const scattering::ParticleType format) {
    switch (format) {
    case scattering::ParticleType::Random:
        output << "random";
        break;
    case scattering::ParticleType::AzimuthallyRandom:
        output << "azimuthally random";
        break;
    case scattering::ParticleType::General:
        output << "general";
        break;
    }
    return output;
}

}
