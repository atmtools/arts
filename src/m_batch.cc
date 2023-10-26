/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   m_batch.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2004-09-15 

  \brief  Workspace functions for doing batch calculations.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include "gridded_fields.h"
using namespace std;

#include "arts_omp.h"
#include <workspace.h>
#include "math_funcs.h"
#include "physics_funcs.h"
#include "rte.h"
#include "xml_io.h"

inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);
using GriddedFieldGrids::GFIELD3_P_GRID;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated 

   2008-07-21 Stefan Buehler */
void ForLoop(const Workspace& ws,
             // WS Input:
             const Agenda& forloop_agenda,
             // Control Parameters:
             const Index& start,
             const Index& stop,
             const Index& step) {
  for (Index i = start; i <= stop; i += step) {
    forloop_agendaExecute(ws, i, forloop_agenda);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ybatchCalc(const Workspace& ws,
                // WS Output:
                ArrayOfVector& ybatch,
                ArrayOfArrayOfVector& ybatch_aux,
                ArrayOfMatrix& ybatch_jacobians,
                // WS Input:
                const Index& ybatch_start,
                const Index& ybatch_n,
                const Agenda& ybatch_calc_agenda,
                // Control Parameters:
                const Index& robust) {
  Index first_ybatch_index = 0;

  ArrayOfString fail_msg;
  bool do_abort = false;

  // We allow a start index ybatch_start that is different from 0. We
  // will calculate ybatch_n jobs starting at the start
  // index. Internally, we count from zero, which is the right
  // counting for the output array ybatch. When we call
  // ybatch_calc_agenda, we add ybatch_start to the internal index
  // count.

  // We create a counter, so that we can generate nice output about
  // how many jobs are already done. (All parallel threads later will
  // increment this, so that we really get an accurate total count!)
  Index job_counter = 0;

  // Resize the output arrays:
  ybatch.resize(ybatch_n);
  ybatch_aux.resize(ybatch_n);
  ybatch_jacobians.resize(ybatch_n);
  for (Index i = 0; i < ybatch_n; i++) {
    ybatch[i].resize(0);
    ybatch_aux[i].resize(0);
    ybatch_jacobians[i].resize(0, 0);
  }

  // Go through the batch:

  if (ybatch_n) {
#pragma omp parallel for schedule(dynamic) if (!arts_omp_in_parallel() && ybatch_n > 1)
    for (Index ybatch_index = first_ybatch_index; ybatch_index < ybatch_n;
         ybatch_index++) {
      Index l_job_counter;  // Thread-local copy of job counter.

      if (do_abort) continue;
#pragma omp critical(ybatchCalc_job_counter)
      { l_job_counter = ++job_counter; }

      {
        ostringstream os;
        os << "  Job " << l_job_counter << " of " << ybatch_n << ", Index "
           << ybatch_start + ybatch_index << ", Thread-Id "
           << arts_omp_get_thread_num() << "\n";
      }

      try {
        Vector y;
        ArrayOfVector y_aux;
        Matrix jacobian;

        ybatch_calc_agendaExecute(ws,
                                  y,
                                  y_aux,
                                  jacobian,
                                  ybatch_start + ybatch_index,
                                  ybatch_calc_agenda);

        if (y.nelem()) {
#pragma omp critical(ybatchCalc_assign_y)
          ybatch[ybatch_index] = y;
#pragma omp critical(ybatchCalc_assign_y_aux)
          ybatch_aux[ybatch_index] = y_aux;

          // Dimensions of Jacobian:
          const Index Knr = jacobian.nrows();
          const Index Knc = jacobian.ncols();

          if (Knr != 0 || Knc != 0) {
            if (Knr != y.nelem()) {
              ostringstream os;
              os << "First dimension of Jacobian must have same length as the measurement *y*.\n"
                 << "Length of *y*: " << y.nelem() << "\n"
                 << "Dimensions of *jacobian*: (" << Knr << ", " << Knc
                 << ")\n";
              // A mismatch of the Jacobian dimension is a fatal error
              // and should result in program termination. By setting abort
              // to true, this will result in a runtime error in the catch
              // block even if robust == 1
#pragma omp critical(ybatchCalc_setabort)
              do_abort = true;

              throw runtime_error(os.str());
            }

            ybatch_jacobians[ybatch_index] = jacobian;

            // After creation, all individual Jacobi matrices in the array will be
            // empty (size zero). No need for explicit initialization.
          }
        }
      } catch (const std::exception& e) {
        if (robust && !do_abort) {
          ostringstream os;
          os << "WARNING! Job at ybatch_index " << ybatch_start + ybatch_index
             << " failed.\n"
             << "y Vector in output variable ybatch will be empty for this job.\n"
             << "The runtime error produced was:\n"
             << e.what() << "\n";
        } else {
          // The user wants the batch job to fail if one of the
          // jobs goes wrong.
#pragma omp critical(ybatchCalc_setabort)
          do_abort = true;

          ostringstream os;
          os << "  Job at ybatch_index " << ybatch_start + ybatch_index
             << " failed. Aborting...\n";
        }
        ostringstream os;
        os << "Run-time error at ybatch_index " << ybatch_start + ybatch_index
           << ": \n"
           << e.what();
#pragma omp critical(ybatchCalc_push_fail_msg)
        fail_msg.push_back(os.str());
      }
    }
  }

  if (fail_msg.size()) {
    ostringstream os;

    if (!do_abort) os << "\nError messages from failed batch cases:\n";
    for (ArrayOfString::const_iterator it = fail_msg.begin();
         it != fail_msg.end();
         it++)
      os << *it << '\n';

    if (do_abort)
      throw runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ybatchMetProfilesClear(const Workspace& ws,
                            //Output
                            ArrayOfVector& ybatch,
                            //Input
                            const ArrayOfArrayOfSpeciesTag& abs_species,
                            const Agenda& met_profile_calc_agenda,
                            const Vector& f_grid,
                            const Matrix& met_amsu_data,
                            const Matrix& sensor_pos,
                            const SurfaceField& surface_field,
                            //Keyword
                            const Index& nelem_p_grid,
                            const String& met_profile_path) {
  AtmField atm_field;
  ArrayOfGriddedField3 pnd_field_raw;
  Vector z_grid;
  Matrix sensor_los;
  Index cloudbox_on = 0;
  ArrayOfIndex cloudbox_limits;
  Vector y;
  Index no_profiles = met_amsu_data.nrows();
  //Index no_profiles = met_profile_basenames.nelem();
  // The humidity data is stored as  an ArrayOfTensor3 whereas
  // vmr_field_raw is an ArrayOfArrayOfTensor3
  GriddedField3 vmr_field_raw_h2o;

  y.resize(f_grid.nelem());
  ybatch.resize(no_profiles);

  Vector sat_za_from_profile;
  sat_za_from_profile = met_amsu_data(Range(joker), 3);
  Numeric sat_za;

  sensor_los.resize(1, 1);

  Vector lat, lon;
  lat = met_amsu_data(Range(joker), 0);
  lon = met_amsu_data(Range(joker), 1);

//   Vector oro_height;
//   oro_height = met_amsu_data(Range(joker), 5);

  for (Index i = 0; i < no_profiles; ++i) {
    ostringstream lat_os, lon_os;

    Index lat_prec = 3;
    if (lat[i] < 0) lat_prec--;
    if (abs(lat[i]) >= 10) {
      lat_prec--;
      if (abs(lat[i]) >= 100) lat_prec--;
    }

    lat_os.setf(ios::showpoint | ios::fixed);
    lat_os << setprecision((int)lat_prec) << lat[i];

    Index lon_prec = 4;
    if (lon[i] < 0) lon_prec--;
    if (abs(lon[i]) >= 10) {
      lon_prec--;
      if (abs(lon[i]) >= 100) lon_prec--;
    }
    lon_os.setf(ios::showpoint | ios::fixed);
    lon_os << setprecision((int)lon_prec) << lon[i];
    cout << lat_os.str() << endl;
    cout << lon_os.str() << endl;

    sat_za = sat_za_from_profile[i];

    sensor_los(Range(joker), 0) =
        180.0 -
        (asin(surface_field.ellipsoid[0] * sin(sat_za * PI / 180.) / sensor_pos(0, 0))) *
            180. / PI;
    cout << "sensor_los" << sat_za_from_profile[i] << endl;
    cout << "sensor_los" << sat_za << endl;
    cout << "sensor_los" << sensor_los << endl;
    //Reads the t_field_raw from file

    GriddedField3 gf3;
    xml_read_from_file(met_profile_path + "profile.lat_" + lat_os.str() +
                           ".lon_" + lon_os.str() + ".t.xml",
                       gf3);
    atm_field[Atm::Key::t] = gf3;

    //Reads the z_field_raw from file
    xml_read_from_file(met_profile_path + "profile.lat_" + lat_os.str() +
                           ".lon_" + lon_os.str() + ".p.xml",
                       gf3);
    atm_field[Atm::Key::p] = gf3;

    //Reads the humidity from file - it is only an ArrayofTensor3
    // The vmr_field_raw is an ArrayofArrayofTensor3 where the outer
    // array is for species
    xml_read_from_file(met_profile_path + "profile.lat_" + lat_os.str() +
                           ".lon_" + lon_os.str() + ".H2O.xml",
                       vmr_field_raw_h2o);
    //xml_read_from_file("/home/home01/rekha/uk/profiles/sat_vmr/profile.lat_"+lat_os.str()//+".lon_"+lon_os.str() + ".H2O_es.xml",
    //                   vmr_field_raw_h2o);

    cout
        << "--------------------------------------------------------------------------"
        << endl;
    cout << "The file"
         << met_profile_path + "profile.lat_" + lat_os.str() + ".lon_" +
                lon_os.str()
         << "is executed now" << endl;
    cout
        << "--------------------------------------------------------------------------"
        << endl;
    xml_write_to_file("profile_number.xml", i, FILE_TYPE_ASCII, 0);
    // the first element of the species is water vapour.

    // N_p is the number of elements in the pressure grid
    //z_surface(0,0) = oro_height[i]+ 0.01;
    SurfaceField sf2=surface_field;
    sf2[Surf::Key::h] = atm_field[Atm::Key::p].get<const GriddedField3&>().get_numeric_grid(0)[0];
    cout << "z_surface" << sf2 << endl;
    const Vector& tfr_z_grid =
        atm_field[Atm::Key::p].get<const GriddedField3&>().get_numeric_grid(GFIELD3_P_GRID);
    Index N_p = tfr_z_grid.nelem();

    atm_field[abs_species[0].Species()] = vmr_field_raw_h2o;

    // the second element of the species.  the first 3 Tensors in the
    //array are the same .  They are pressure grid, latitude grid and
    // longitude grid.  The third tensor which is the vmr is set to a
    // constant value of 0.782.
    atm_field[abs_species[1].Species()] = 0.782;

    // the second element of the species.  the first 3 Tensors in the
    //array are the same .  They are pressure grid, latitude grid and
    // longitude grid.  The third tensor which is the vmr is set to a
    // constant value of 0.209.
    atm_field[abs_species[2].Species()] = 0.209;

    //xml_write_to_file(met_profile_basenames[i]+ ".N2.xml", vmr_field_raw[1]);
    //xml_write_to_file(met_profile_basenames[i]+ ".O2.xml", vmr_field_raw[2]);

    //Making a p_grid with the first and last element taken from the profile.
    // this is because of the extrapolation problem.

    VectorNLogSpace(
        z_grid, nelem_p_grid, tfr_z_grid[0], tfr_z_grid[N_p - 1]);
    cout << "t_field_raw[0](0,0,0)" << tfr_z_grid[0] << endl;
    cout << "t_field_raw[0](N_p -1,0,0)" << tfr_z_grid[N_p - 1] << endl;
    xml_write_to_file("z_grid.xml", z_grid, FILE_TYPE_ASCII, 0);

    // executing the met_profile_calc_agenda
    met_profile_calc_agendaExecute(ws,
                                   y,
                                   atm_field,
                                   pnd_field_raw,
                                   sensor_los,
                                   cloudbox_on,
                                   cloudbox_limits,
                                   sf2,
                                   met_profile_calc_agenda);

    //putting in the spectra *y* for each profile
    ybatch[i] = y;

  }  // closing the loop over profile basenames
}

/* Workspace method: Doxygen documentation will be auto-generated */
void DOBatchCalc(const Workspace& ws,
                 ArrayOfTensor7& dobatch_cloudbox_field,
                 ArrayOfTensor5& dobatch_radiance_field,
                 ArrayOfTensor4& dobatch_irradiance_field,
                 ArrayOfTensor5& dobatch_spectral_irradiance_field,
                 const Index& ybatch_start,
                 const Index& ybatch_n,
                 const Agenda& dobatch_calc_agenda,
                 const Index& robust) {
  Index first_ybatch_index = 0;

  ArrayOfString fail_msg;
  bool do_abort = false;

  // We allow a start index ybatch_start that is different from 0. We
  // will calculate ybatch_n jobs starting at the start
  // index. Internally, we count from zero, which is the right
  // counting for the output array ybatch. When we call
  // ybatch_calc_agenda, we add ybatch_start to the internal index
  // count.

  // We create a counter, so that we can generate nice output about
  // how many jobs are already done. (All parallel threads later will
  // increment this, so that we really get an accurate total count!)
  Index job_counter = 0;

  // Resize the output arrays:
  dobatch_cloudbox_field.resize(ybatch_n);
  dobatch_radiance_field.resize(ybatch_n);
  dobatch_irradiance_field.resize(ybatch_n);
  dobatch_spectral_irradiance_field.resize(ybatch_n);

  // Go through the batch:

  if (ybatch_n) {  
#pragma omp parallel for schedule(dynamic) if (!arts_omp_in_parallel() && \
                                               ybatch_n > 1)
    for (Index ybatch_index = first_ybatch_index; ybatch_index < ybatch_n;
         ybatch_index++) {
      Index l_job_counter;  // Thread-local copy of job counter.

      if (do_abort) continue;
#pragma omp critical(dobatchCalc_job_counter)
      { l_job_counter = ++job_counter; }

      {
        ostringstream os;
        os << "  Job " << l_job_counter << " of " << ybatch_n << ", Index "
           << ybatch_start + ybatch_index << ", Thread-Id "
           << arts_omp_get_thread_num() << "\n";
      }

      try {
        Tensor7 cloudbox_field;
        Tensor5 radiance_field;
        Tensor4 irradiance_field;
        Tensor5 spectral_irradiance_field;

        dobatch_calc_agendaExecute(ws,
                                   cloudbox_field,
                                   radiance_field,
                                   irradiance_field,
                                   spectral_irradiance_field,
                                   ybatch_start + ybatch_index,
                                   dobatch_calc_agenda);

#pragma omp critical(dobatchCalc_assign_cloudbox_field)
        dobatch_cloudbox_field[ybatch_index] = cloudbox_field;
#pragma omp critical(dobatchCalc_assign_radiance_field)
        dobatch_radiance_field[ybatch_index] = radiance_field;
#pragma omp critical(dobatchCalc_assign_irradiance_field)
        dobatch_irradiance_field[ybatch_index] = irradiance_field;
#pragma omp critical(dobatchCalc_assign_spectral_irradiance_field)
        dobatch_spectral_irradiance_field[ybatch_index] =
            spectral_irradiance_field;

      } catch (const std::exception& e) {
        if (robust && !do_abort) {
          ostringstream os;
          os << "WARNING! Job at ybatch_index " << ybatch_start + ybatch_index
             << " failed.\n"
             << "element in output variables will be empty for this job.\n"
             << "The runtime error produced was:\n"
             << e.what() << "\n";
        } else {
          // The user wants the batch job to fail if one of the
          // jobs goes wrong.
#pragma omp critical(dobatchCalc_setabort)
          do_abort = true;

          ostringstream os;
          os << "  Job at ybatch_index " << ybatch_start + ybatch_index
             << " failed. Aborting...\n";
        }
        ostringstream os;
        os << "Run-time error at ybatch_index " << ybatch_start + ybatch_index
           << ": \n"
           << e.what();
#pragma omp critical(dobatchCalc_push_fail_msg)
        fail_msg.push_back(os.str());
      }
    }
  }

  if (fail_msg.size()) {
    ostringstream os;

    if (!do_abort) os << "\nError messages from failed batch cases:\n";
    for (ArrayOfString::const_iterator it = fail_msg.begin();
         it != fail_msg.end();
         it++)
      os << *it << '\n';

    if (do_abort)
      throw runtime_error(os.str());
  }
}
