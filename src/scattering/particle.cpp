#include <scattering/particle.h>
#include <scattering/arts_ssdb.h>

namespace scattering {
    Particle Particle::from_ssdb(std::string path) {
        return arts_ssdb::ParticleFile(path).to_particle();
    }
}
