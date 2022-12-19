#include <scattering/particle_habit.h>
#include <scattering/arts_ssdb.h>

namespace scattering {
    ParticleHabit ParticleHabit::from_ssdb(std::string path) {
        return arts_ssdb::HabitFolder(path).to_particle_habit();
    }

    std::ostream& operator<<(std::ostream& out, const ParticleHabit& habit) {
        out << "Particle habit containing " << habit.particles_.size() << " particles: " << std::endl;
        for (const auto &particle : habit.particles_) {
            out << particle << std::endl;
        }
        out << std::endl;
        return out;
    }
}
