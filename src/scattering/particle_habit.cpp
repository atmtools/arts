#include <scattering/particle_habit.h>
#include <scattering/arts_ssdb.h>

namespace scattering {
    ParticleHabit ParticleHabit::from_ssdb(std::string path) {
        return arts_ssdb::HabitFolder(path).to_particle_habit();
    }
}
