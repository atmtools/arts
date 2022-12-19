#include <scattering/particle.h>
#include <scattering/arts_ssdb.h>

namespace scattering {
    Particle Particle::from_ssdb(std::string path) {
        return arts_ssdb::ParticleFile(path).to_particle();
    }

  std::ostream &operator<<(std::ostream &output,
                          const scattering::Particle &particle) {
    auto name = particle.get_name();
    auto source = particle.get_source();
    auto refractive_index = particle.get_refractive_index();

    output << "Scattering particle: " << name;
    if (source != "") {
      output << " (" << source << ")";
    }
    output << std::endl;
    if (refractive_index != "") {
      output << std::setw(30)
            << "\tRefractive index: " << particle.get_refractive_index()
            << std::endl;
    }
    output << std::setw(30) << std::left << "\tParticle type: " << particle.get_particle_type() << std::endl;
    output << std::setw(30) << std::left << "\tScattering data format: " << particle.get_data_format()
          << std::endl;
    output << std::setw(30) << std::left << "\tVolume-equivalent diameter: " << particle.get_d_eq() << " m" << std::endl;
    output << std::setw(30) << std::left << "\tMaximum diameter: " << particle.get_d_max() << " m" << std::endl;
    output << std::setw(30) << std::left << "\tMass: " << particle.get_mass() << " kg" << std::endl;
    output << std::endl;
    return output;
  }
}
