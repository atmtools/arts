#include <auto_md.h>

int main(int argc, char ** argv) {
  define_wsv_groups();
  define_wsv_data();
  define_wsv_map();

  const bool do_openmp = argc > 1 and argv[1][0] == 'y';

  Workspace ws{};
  WorkspaceOmpParallelCopyGuard wss{ws, do_openmp};

  std::cerr << "arts_omp_in_parallel(): " << arts_omp_in_parallel() <<'\n';
  std::cerr << "arts_omp_get_max_threads() not_eq 1: " << (arts_omp_get_max_threads() not_eq 1) << '\n';

  #pragma omp parallel for firstprivate(wss) if (do_openmp)
  for (Index i=0; i<100; i++) {}
}
