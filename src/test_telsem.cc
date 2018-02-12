/* Copyright (C) 2018
   Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*!
  \file   test_telsem.cc
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   Mon February 12 12:34:56 2018

  \brief  Tests for TELSEM2 interface.
*/

#include <iostream>
#include <string>
#include <fstream>

#include "arts.h"
#include "telsem.h"

/** Test reading of TELSEM emissivity interpolation
 *
 * This function tests the interpolation of TELSEM emissivities
 * by comparison with the results from the test that is included
 * in the fortran module.
 *
 * @param atlas_file The path to the atlas file which was used in the FORTRAN
 *        module tests
 * @param result_path The folder containing the results of the FORTRAN module
 *        tests
 * @param resolution The resolution to use for the lat/lon map.
 * @param theta The incidence angle to which to interpolate the frequencies.
 * @param frequencies The frequencies [GHz] (!!!) for which to interpolate the emissivities
 */
Numeric test_telsem_interpolate(std::string atlas_file,
                                std::string result_path,
                                Numeric resolution,
                                Numeric theta,
                                Vector frequencies)
{
    TelsemAtlas atlas(atlas_file);

    Index n_freqs = frequencies.nelem();
    std::vector<std::ifstream> results_h(n_freqs), results_v(n_freqs);

    for (Index i = 0; i < n_freqs; ++i) {
        std::string filename_h = result_path;
        std::string filename_v = result_path;
        filename_h += "/emisH_IND_MULT" + std::to_string(i + 1) + ".txt";
        filename_v += "/emisV_IND_MULT" + std::to_string(i + 1) + ".txt";

        results_h[i] = std::ifstream(filename_h, std::ifstream::in);
        results_v[i] = std::ifstream(filename_v, std::ifstream::in);
    }

    Index n_lat = static_cast<Index>(180.0 / resolution);
    Index n_lon = static_cast<Index>(360.0 / resolution);

    Numeric error = 0.0;

    for (Index i = n_lat - 1; i >= 0; --i) {
        for (Index j = 0; j < n_lon; ++j) {

            // An offset of half of the latitude resolution is added to avoid hitting
            // the boundary between to cells.

            Numeric lat = 0.125 + resolution / 2.0 -90.0 + static_cast<Numeric>(i) * resolution;
            Numeric lon = 0.125 + resolution / 2.0       + static_cast<Numeric>(j) * resolution;

            Index cellnumber = atlas.calc_cellnum(lat, lon);

            Vector emis_h(3), emis_h_interp(n_freqs), emis_h_interp_ref(n_freqs),
                   emis_v(3), emis_v_interp(n_freqs), emis_v_interp_ref(n_freqs);

            // Read reference emissivities.
            for (Index k = 0; k < n_freqs; ++k) {
                results_h[k] >> emis_h_interp_ref[k];
                results_v[k] >> emis_v_interp_ref[k];
            }

            emis_h_interp = 0.0;
            emis_v_interp = 0.0;

            // Read emissivities from atlas.
            if (atlas.contains(cellnumber)) {

                emis_h = atlas.get_emis_h(cellnumber);
                emis_v = atlas.get_emis_v(cellnumber);

                Index class1 = atlas.get_class1(cellnumber);
                Index class2 = atlas.get_class2(cellnumber);

                for (Index k = 0; k < n_freqs; ++k) {
                    std::tie(emis_v_interp[k], emis_h_interp[k]) = atlas.emis_interp(theta,
                                                                                     frequencies[k],
                                                                                     class1,
                                                                                     class2,
                                                                                     emis_v,
                                                                                     emis_h);
                }
            }

            for (Index k = 0; k < n_freqs; ++k) {
                error = std::max(error, std::fabs(emis_h_interp[k] - emis_h_interp_ref[k]));
                error = std::max(error, std::fabs(emis_v_interp[k] - emis_v_interp_ref[k]));
            }
        }
    }
    return error;
}

/** Test reading of TELSEM emissivities
 *
 * This function tests the reading of the telsem atlas by creating a lat/lon map
 * of emissivities in the atlas and comparing with the results from the test
 * that is included in the fortran module.
 *
 * @param atlas_file The path to the atlas file which was used in the FORTRAN
 *        module tests
 * @param result_path The folder containing the results of the FORTRAN module
 *        tests
 * @param resolution The resolution used to generate the lat/lon map.
 */
Numeric test_telsem_read(String atlas_file,
                         String result_path,
                         Numeric resolution)
{
    TelsemAtlas atlas(atlas_file);

    std::vector<std::ifstream> results_h(3), results_v(3);

    for (Index i = 0; i < 3; ++i) {
        std::string filename_h = result_path;
        std::string filename_v = result_path;
        filename_h += "/emisH" + std::to_string(i + 1) + ".txt";
        filename_v += "/emisV" + std::to_string(i + 1) + ".txt";

        results_h[i] = std::ifstream(filename_h, std::ifstream::in);
        results_v[i] = std::ifstream(filename_v, std::ifstream::in);
    }

    Index n_lat = static_cast<Index>(180.0 / resolution);
    Index n_lon = static_cast<Index>(360.0 / resolution);

    Numeric error = 0.0;

    for (Index i = n_lat - 1; i >= 0; --i) {
        for (Index j = 0; j < n_lon; ++j) {

            // An offset of half of the latitude resolution is added to avoid hitting
            // the boundary between to cells.

            Numeric lat = 0.125 + resolution / 2.0 -90.0 + static_cast<Numeric>(i) * resolution;
            Numeric lon = 0.125 + resolution / 2.0       + static_cast<Numeric>(j) * resolution;

            Index cellnumber = atlas.calc_cellnum(lat, lon);

            Vector emis_h(3), emis_h_ref(3), emis_v(3), emis_v_ref(3);

            // Read reference emissivities.
            for (Index k = 0; k < 3; ++k) {
                results_h[k] >> emis_h_ref[k];
                results_v[k] >> emis_v_ref[k];
            }

            emis_h = 0.0;
            emis_v = 0.0;

            // Read emissivities from atlas.
            if (atlas.contains(cellnumber)) {
                emis_h = atlas.get_emis_h(cellnumber);
                emis_v = atlas.get_emis_v(cellnumber);
            }

            for (Index k = 0; k < 3; ++k) {
                error = std::max(error, std::fabs(emis_h[k] - emis_h_ref[k]));
                error = std::max(error, std::fabs(emis_v[k] - emis_v_ref[k]));
            }
        }
    }
    return error;
}

int main(int argc, const char ** argv) {

    if (argc != 4) {
        std::cout <<
            "\nThis test uses the test results of the TELSEM2 fortran\n"
            "module to test the ARTS TELSEM interface. \n\n Usage:\n"
            "./test_telsem <atlas_path> <results_folder> <resolution> \n"
            "where:\n - atlas_path: Path to the TELSEM2 atlas to load \n"
            " - results_folder: The folder containing the result files \n"
            "\t of the TELSEM2 tests.\n - resolution: The resolution used"
            "to generate the lat/lon map.\n\n Note that for the interpolation"
            "test the frequencies must match.\n\n";
        return 0;
    }

    String atlas_file  = argv[1];
    String result_path = argv[2];
    Numeric resolution = std::stoi(argv[3]);

    std::cout << "Atlas file:  " << atlas_file << std::endl;
    std::cout << "Result path: " << result_path << std::endl;

    // Reading of emissivities.

    Numeric error = test_telsem_read(atlas_file, result_path, 2.0);
    std::cout << "Maximum error reading emissivities:       " << error << std::endl;

    // Interpolation of emissivities.

    Vector frequencies = {6.0, 25.0, 31.4, 60.0, 190.0};
    Numeric theta      = 15.0; // Incidence angle
    error = test_telsem_interpolate(atlas_file, result_path, 2.0, theta, frequencies);
    std::cout << "Maximum error interpolating emissivities: " << error << std::endl;
    return 0;
}
