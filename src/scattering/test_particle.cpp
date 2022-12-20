#include <numbers>
#include <filesystem>
#include "scattering/particle.h"

using std::numbers::pi_v;

// Test serialization of particles.
bool test_io() {
    scattering::math::Vector<double> frequencies{2};
    frequencies << 89e9, 183e9;
    scattering::math::Vector<double> temperatures{2};
    temperatures << 240, 260;
    scattering::math::Vector<double> scattering_angles = scattering::math::Vector<double>::LinSpaced(901, 0, pi_v<double>);

    scattering::Particle part = scattering::Particle::liquid_sphere(
        frequencies,
        temperatures,
        scattering_angles,
        0.5e-3
        );

    std::cout << part << std::endl;

    std::ofstream output("test_particle.bin", ios::binary);
    part.serialize(output);
    output.close();
    std::ifstream input("test_particle.bin", ios::binary);
    auto part_other = scattering::Particle::deserialize(input);
    input.close();
    std::filesystem::remove("test_particle.bin");
    assert(part_other.get_data_format() == scattering::DataFormat::Gridded);

    auto pm = part.get_phase_matrix_data();
    auto pm_other = part.get_phase_matrix_data();
    scattering::math::Tensor<double, 0> delta = (pm - pm_other).abs().maximum();
    if (delta.coeff() > 1e-6) return false;

    auto part_spectral = part.to_spectral();
    output.open("test_particle_spectral.bin", ios::binary);
    part_spectral.serialize(output);
    output.close();
    input.open("test_particle_spectral.bin", ios::binary);
    auto part_other_spectral = scattering::Particle::deserialize(input);
    input.close();
    std::filesystem::remove("test_particle_spectral.bin");
    assert(part_other_spectral.get_data_format() == scattering::DataFormat::Spectral);

    auto pm_spectral = part_spectral.get_phase_matrix_data_spectral();
    auto pm_spectral_other = part_spectral.get_phase_matrix_data_spectral();
    delta = (pm - pm_other).abs().maximum();
    if (delta.coeff() > 1e-6) return false;

    pm_other = part_spectral.get_phase_matrix_data();
    delta = (pm - pm_other).abs().maximum();
    if (delta.coeff() > 1e-6) return false;

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
