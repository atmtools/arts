#ifndef montecarlo_h
#define montecarlo_h

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include "messages.h"
#include "arts.h"
#include "ppath.h"
#include "matpackI.h"
#include "special_interp.h"
#include "check_input.h"
#include <stdexcept>
#include <cmath>
#include "rte.h"
#include "lin_alg.h"
#include "auto_md.h"
#include "logic.h"
#include "physics_funcs.h"
#include "xml_io.h"
#include "rng.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;

void Cloudbox_ppath_rteCalc(
			     Ppath&                ppathcloud,
			     Ppath&                ppath,
			     Ppath&                ppath_step,
			     Vector&               rte_pos,
			     Vector&               rte_los,
			     Vector&               cum_l_step,
			     ArrayOfMatrix&        TArray,
			     ArrayOfMatrix&        ext_matArray,
			     ArrayOfVector&        abs_vecArray,
			     Vector&               t_ppath,
			     Vector&               scat_za_grid,
			     Vector&               scat_aa_grid,
			     Tensor3&              ext_mat,
			     Matrix&               abs_vec,
			     Numeric&              rte_pressure,
			     Numeric&              rte_temperature,
			     Vector&               rte_vmr_list,
			     Matrix&               i_rte,
			     GridPos&              rte_gp_p,
			     GridPos&              rte_gp_lat,
			     GridPos&              rte_gp_lon,
			     Matrix&               i_space,
			     Matrix&               ground_emission,
			     Matrix&               ground_los, 
			     Tensor4&              ground_refl_coeffs,
			     Index&                f_index,
			     Index&                scat_za_index,
			     Index&                scat_aa_index,
			     const Agenda&         ppath_step_agenda,
			     const Index&          atmosphere_dim,
			     const Vector&         p_grid,
			     const Vector&         lat_grid,
			     const Vector&         lon_grid,
			     const Tensor3&        z_field,
			     const Matrix&         r_geoid,
			     const Matrix&         z_ground,
			     const ArrayOfIndex&   cloudbox_limits,
			     const Index&          record_ppathcloud,
			     const Index&          record_ppath,
			     const Agenda&         opt_prop_gas_agenda,
			     const Agenda&         opt_prop_part_agenda,
			     const Agenda&         scalar_gas_absorption_agenda,
			     const Index&          stokes_dim,
			     const Tensor3&        t_field,
			     const Tensor4&        vmr_field,
			     const Agenda&         rte_agenda,
			     const Agenda&         i_space_agenda,
			     const Agenda&         ground_refl_agenda,
			     const Vector&         f_grid,
			     const Index&          photon_number,
			     const Index&          scattering_order);


void cloudbox_ppath_start_stepping(
				   Ppath&          ppath,
				   const Index&          atmosphere_dim,
				   ConstVectorView       p_grid,
				   ConstVectorView       lat_grid,
				   ConstVectorView       lon_grid,
				   ConstTensor3View      z_field,
				   ConstMatrixView       r_geoid,
				   ConstMatrixView       z_ground,
				   ConstVectorView       rte_pos,
				   ConstVectorView       rte_los );
	  

void cum_l_stepCalc(
		      Vector& cum_l_step,
		      const Ppath& ppath
		      );


Matrix interp( ConstVectorView itw,
	       ArrayOfMatrix a,    
	       const GridPos&  tc );

Vector interp( ConstVectorView itw,
                ArrayOfVector a,    
	       const GridPos&  tc );

void interpTArray(Matrix& T,
		  Vector& Kabs,
		  Numeric& temperature,
		  Vector& rte_pos,
		  Vector& rte_los,
		  ArrayOfGridPos& gp,
		  const ArrayOfMatrix& TArray,
		  const ArrayOfMatrix& ext_matArray,
		  const ArrayOfVector& abs_vecArray,
		  const Vector& t_ppath,
		  const Vector& cum_l_step,
		  const Numeric& pathlength,
		  const Index& stokes_dim,
		  const Ppath& ppath
		 );

void Sample_los (
		 Vector& rte_los,
		 Rng& rng
		 );

void Sample_ppathlength (
			 Numeric& pathlength, 
			 Numeric& g,
			 Rng& rng,
			 const ArrayOfMatrix& ext_matArray
			 );

void TArrayCalc(
		//output
		ArrayOfMatrix& TArray,
		ArrayOfMatrix& ext_matArray,
		ArrayOfVector& abs_vecArray,
		Vector& t_ppath,
		Vector& scat_za_grid,
		Vector& scat_aa_grid,
		Tensor3& ext_mat,
		Matrix& abs_vec,
		Numeric&   rte_pressure,
		Numeric&   rte_temperature,
		Vector&    rte_vmr_list,
		Index&    scat_za_index,
		Index&    scat_aa_index,
		//input
		const Ppath& ppath,
		const Agenda& opt_prop_gas_agenda,
		const Agenda& opt_prop_part_agenda,
		const Agenda& scalar_gas_absorption_agenda,
		const Index& stokes_dim,
		const Vector&    p_grid,
		const Vector&    lat_grid,
		const Vector&    lon_grid,
		const Tensor3&   t_field,
		const Tensor4&   vmr_field,
		const Index&     atmosphere_dim
 		);

#endif  // montecarlo_h
