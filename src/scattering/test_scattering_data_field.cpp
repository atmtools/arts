#include <numbers>
#include <filesystem>
#include "scattering/mie.h"
#include "scattering/maths.h"
#include "scattering/scattering_data_field.h"

using std::numbers::pi_v;

namespace scattering {

    /// Create a scattering data field for test purposes.
    ScatteringDataFieldGridded<double> create_scattering_data_field(
        math::Vector<double> frequencies,
        math::Vector<double> temperatures,
        math::Vector<double> scattering_angles,
        double radius
        ) {

        std::array<math::Index, 7> dims = {
            frequencies.size(),
            temperatures.size(),
            1,
            1,
            1,
            scattering_angles.size(),
            6
        };
        math::Tensor<double, 7> data{dims};

        for (math::Index i_f = 0; i_f < frequencies.size(); ++i_f) {
            for (math::Index i_t = 0; i_t < temperatures.size(); ++i_t) {
                MieSphere<double> sphere = MieSphere<double>::Ice(
                    frequencies[i_f],
                    temperatures[i_t],
                    radius,
                    scattering_angles
                    );
                auto sm = sphere.get_scattering_matrix_compact();
                for (math::Index i_a = 0; i_a < scattering_angles.size(); ++i_a) {
                    for (math::Index i_s = 0; i_s < 6; ++i_s) {
                        std::array<math::Index, 7> index{i_f, i_t, 0, 0, 0, i_a, i_s};
                        data.coeffRef(index) = sm(i_a, i_s);
                    }
                }
            }
        }

        return ScatteringDataFieldGridded<double>(
            frequencies,
            temperatures,
            math::Vector<double>::Zero(1),
            math::Vector<double>::Zero(1),
            math::Vector<double>::Zero(1),
            scattering_angles,
            data
            );
    }
}

// Test serialization of gridded scattering data field.
bool test_io() {
    scattering::math::Vector<double> frequencies{2};
    frequencies << 89e9, 183e9;
    scattering::math::Vector<double> temperatures{2};
    temperatures << 240, 260;
    scattering::math::Vector<double> scattering_angles = scattering::math::Vector<double>::LinSpaced(901, 0, pi_v<double>);

    auto sdf = scattering::create_scattering_data_field(
        frequencies,
        temperatures,
        scattering_angles,
        1e-3
        );
    auto sdf_spectral = sdf.to_spectral();

    // Write and read gridded data.
    std::ofstream output("test_scattering_data_field.bin", ios::binary);
    sdf.serialize(output);
    output.close();
    std::ifstream input("test_scattering_data_field.bin", ios::binary);
    auto sdf_other = scattering::ScatteringDataFieldGridded<double>::deserialize(input);
    input.close();
    auto data = sdf.get_data();
    auto data_other = sdf_other.get_data();
    std::filesystem::remove("test_scattering_data_field.bin");
    scattering::math::Tensor<double, 0> delta = (data - data_other).abs().maximum();
    if (delta.coeff() > 1e-6) return false;

    // Write and read spectral data.
    output.open("test_scattering_data_field_spectral.bin", ios::binary);
    sdf_spectral.serialize(output);
    output.close();
    input.open("test_scattering_data_field_spectral.bin", ios::binary);
    auto sdf_other_spectral = scattering::ScatteringDataFieldSpectral<double>::deserialize(input);
    input.close();
    auto data_spectral = sdf_spectral.get_data();
    auto data_other_spectral = sdf_other_spectral.get_data();
    std::filesystem::remove("test_scattering_data_field_spectral.bin");
    delta = (data_spectral - data_other_spectral).abs().maximum();
    if (delta.coeff() > 1e-6) return false;

    data_other = sdf_other_spectral.to_gridded().get_data();
    delta = (data - data_other).abs().maximum();
    std::cout << "DELTA :: " << delta << std::endl;
    if (delta.coeff() > 1e-6) return false;

    // Write and read fully spectral data.
    // NOTE: This is trivial because there is no dependency on the incoming
    // angle.
    output.open("test_scattering_data_field_fully_spectral.bin", ios::binary);
    auto sdf_fully_spectral = sdf_spectral.to_fully_spectral();
    sdf_fully_spectral.serialize(output);
    output.close();
    input.open("test_scattering_data_field_fully_spectral.bin", ios::binary);
    auto sdf_other_fully_spectral = scattering::ScatteringDataFieldFullySpectral<double>::deserialize(input);
    input.close();
    auto data_fully_spectral = sdf_fully_spectral.get_data();
    auto data_other_fully_spectral = sdf_other_fully_spectral.get_data();
    std::filesystem::remove("test_scattering_data_field_fully_spectral.bin");
    delta = (data_fully_spectral - data_other_fully_spectral).abs().maximum();
    if (delta.coeff() > 1e-6) return false;

    data_other = sdf_other_fully_spectral.to_spectral().to_gridded().get_data();
    delta = (data - data_other).abs().maximum();
    std::cout << "DELTA :: " << delta << std::endl;
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
