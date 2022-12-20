#include <numbers>
#include <filesystem>
#include "scattering/particle_habit.h"

using std::numbers::pi_v;

// Test serialization of particles habits.
bool test_io() {
    scattering::math::Vector<double> frequencies{2};
    frequencies << 89e9, 183e9;
    scattering::math::Vector<double> temperatures{2};
    temperatures << 240, 260;
    scattering::math::Vector<double> radii{3};
    radii << 0.1e-3, 0.5e-3, 1e-2;
    scattering::math::Vector<double> scattering_angles = scattering::math::Vector<double>::LinSpaced(901, 0, pi_v<double>);
    scattering::math::Vector<double> pnd{3};
    pnd << 1e2, 1e3, 1e1;

    scattering::ParticleHabit habit = scattering::ParticleHabit::liquid_spheres(
        frequencies,
        temperatures,
        scattering_angles,
        radii
        );


    std::ofstream output("test_particle_habit.bin", ios::binary);
    habit.serialize(output);
    output.close();
    std::ifstream input("test_particle_habit.bin", ios::binary);
    auto habit_other = scattering::ParticleHabit::deserialize(input);
    input.close();
    std::filesystem::remove("test_particle_habit.bin");

    auto pm = habit.get_phase_matrix(183e9, 240.0, 0.0, 0.0, 0.0, pi_v<double> / 2.0, pnd, 4);
    auto pm_other = habit_other.get_phase_matrix(183e9, 240.0, 0.0, 0.0, 0.0, pi_v<double> / 2.0, pnd, 4);
    auto delta = (pm - pm_other).array().abs().maxCoeff();
    if (delta > 1e-6) return false;

    auto habit_spectral = habit.to_spectral();
    output.open("test_particle_habit_spectral.bin", ios::binary);
    habit_spectral.serialize(output);
    output.close();
    input.open("test_particle_habit_spectral.bin", ios::binary);
    auto habit_spectral_other = scattering::ParticleHabit::deserialize(input);
    input.close();
    std::filesystem::remove("test_particle_habit_spectral.bin");

    auto pm_spectral = habit_spectral.get_phase_matrix(183e9, 240.0, 0.0, 0.0, 0.0, pi_v<double> / 2.0, pnd, 4);
    auto pm_spectral_other = habit_spectral_other.get_phase_matrix(183e9, 240.0, 0.0, 0.0, 0.0, pi_v<double> / 2.0, pnd, 4);
    delta = (pm - pm_other).array().abs().maxCoeff();
    if (delta > 1e-6) return false;

    pm_other = habit_spectral_other.get_phase_matrix(183e9, 240.0, 0.0, 0.0, 0.0, pi_v<double> / 2.0, pnd, 4);
    delta = (pm - pm_other).array().abs().maxCoeff();
    if (delta > 1e-4) return false;

    return true;
}


int main(int /*nargs*/, char **/*argv*/) {
    bool passed;

    passed = test_io();
    std::cout << "test_io: ";
    if (passed) {
        std::cout << "PASSED" << std::endl;
    } else {
        std::cout << "FAILED" << std::endl;
        return 1;
    }
}
