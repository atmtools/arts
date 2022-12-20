#include <numbers>
#include <filesystem>
#include "scattering/mie.h"
#include "scattering/maths.h"
#include "scattering/scattering_data_field.h"
#include "scattering/single_scattering_data.h"

using std::numbers::pi_v;

// Test serialization of gridded scattering data field.
bool test_io() {
    scattering::math::Vector<double> frequencies{2};
    frequencies << 89e9, 183e9;
    scattering::math::Vector<double> temperatures{2};
    temperatures << 240, 260;
    scattering::math::Vector<double> scattering_angles = scattering::math::Vector<double>::LinSpaced(901, 0, pi_v<double>);

    scattering::SingleScatteringData sd = scattering::SingleScatteringData::liquid_sphere(
        frequencies,
        temperatures,
        scattering_angles,
        0.5e-3
        );

    std::ofstream output("test_single_scattering_data.bin", ios::binary);
    sd.serialize(output);
    output.close();
    std::ifstream input("test_single_scattering_data.bin", ios::binary);
    auto sd_other = scattering::SingleScatteringData::deserialize(input);
    input.close();
    std::filesystem::remove("test_single_scattering_data.bin");
    assert(sd_other.get_data_format() == scattering::DataFormat::Gridded);

    auto pm = sd.get_phase_matrix_data();
    auto pm_other = sd.get_phase_matrix_data();
    scattering::math::Tensor<double, 0> delta = (pm - pm_other).abs().maximum();
    if (delta.coeff() > 1e-6) return false;

    auto sd_spectral = sd.to_spectral();
    output.open("test_single_scattering_data_spectral.bin", ios::binary);
    sd_spectral.serialize(output);
    output.close();
    input.open("test_single_scattering_data_spectral.bin", ios::binary);
    auto sd_other_spectral = scattering::SingleScatteringData::deserialize(input);
    input.close();
    std::filesystem::remove("test_single_scattering_data_spectral.bin");
    assert(sd_other_spectral.get_data_format() == scattering::DataFormat::Spectral);

    auto pm_spectral = sd_spectral.get_phase_matrix_data_spectral();
    auto pm_spectral_other = sd_spectral.get_phase_matrix_data_spectral();
    delta = (pm - pm_other).abs().maximum();
    if (delta.coeff() > 1e-6) return false;

    pm_other = sd_spectral.get_phase_matrix_data();
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
