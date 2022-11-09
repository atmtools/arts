# Trailing spaces required to be parsed as list.
set(
  generated_files
  SH_to_spat_fly.c
  spat_to_SH_fly.c
  SH_to_spat_omp.c
  spat_to_SH_omp.c
  SH_to_spat_kernel.c
  spat_to_SH_kernel.c
  )

add_custom_command(
  OUTPUT spat_to_SH_fly.c
  COMMAND make -f shtns_make spat_to_SH_fly.c SFX=fly SED=sed
)
add_custom_command(
  OUTPUT SH_to_spat_fly.c
  COMMAND make -f shtns_make SH_to_spat_fly.c SFX=fly SED=sed
)
add_custom_command(
  OUTPUT spat_to_SH_omp.c
  COMMAND make -f shtns_make spat_to_SH_omp.c SFX=omp SED=sed
)
add_custom_command(
  OUTPUT SH_to_spat_omp.c
  COMMAND make -f shtns_make SH_to_spat_omp.c SFX=omp SED=sed
)
add_custom_command(
  OUTPUT spat_to_SH_kernel.c
  COMMAND make -f shtns_make spat_to_SH_kernel.c SFX=kernel SED=sed
)
add_custom_command(
  OUTPUT SH_to_spat_kernel.c
  COMMAND make -f shtns_make SH_to_spat_kernel.c SFX=kernel SED=sed
)

add_custom_target(shtns_generated_files DEPENDS ${generated_files})
