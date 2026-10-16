# HOWTO: Using ThreadSanitizer with OpenMP (Archer)

ThreadSanitizer (TSan) does not understand the synchronization inside the
OpenMP runtime (`libomp`).  Because the runtime is not TSan-instrumented, the
implicit barriers around `#pragma omp parallel` regions are invisible to TSan,
and it then reports data races between accesses that are in fact ordered by
those barriers.  This is what makes an otherwise clean TSan run of
`tests/core/vmr/abs.py` produce dozens of `data race` warnings, e.g. a write
inside `rtepack::SourceVector::init()` "racing" with a read inside a later
`rtepack::constant()` loop.

The fix is [Archer](https://openmp.llvm.org/design/Runtimes.html), an OMPT tool
from the LLVM OpenMP project.  It is loaded into the process next to `libomp`,
observes OpenMP synchronization events through the OMPT interface, and forwards
them to TSan via TSan's dynamic annotation interface (`AnnotateHappensBefore`,
`AnnotateHappensAfter`, ...).  With the tool active, the false positives
disappear while genuine races are still reported.

Archer is not shipped by conda, Homebrew, or Xcode, so it has to be built once.
This document describes the build that was verified on macOS/arm64 with the
conda toolchain used by this repository (clang 23.1.1, CMake 4.4.3, Ninja
1.13.2, LLVM `release/23.x`).

See also the "Thread sanitizer" section in `README.md` for how to build ARTS
itself with `ENABLE_TSAN` and how to preload the TSan runtime for Python.

## Let an LLM do it

The setup below is mechanical and can be delegated to a coding agent
(Claude Code, OpenCode, Codex, ...).  Point it at this file and at the scratch
directory to use:

```text
Set up TSan archer by following the instructions in HOWTO-ARCHER.md.
Use $HOME/myprojects/tsan-archer as TSAN_ROOT.
```

An agent that can read files and run shell commands will do the checkout and
build, run the small verification programs from Step 3, and report the result.
The path in the prompt is arbitrary; it simply becomes `TSAN_ROOT` below.

## Background

- Up to LLVM 12, libomp could be built with TSan annotations
  (`-DLIBOMP_TSAN_SUPPORT=ON`).  That support was removed in LLVM 13
  ([commit 681055ea694b](https://github.com/llvm/llvm-project/commit/681055ea694b),
  "[OpenMP] Remove TSAN annotations from libomp") because Archer supersedes it.
  Do not try the old option on a current LLVM; it no longer exists.
- When no `OMP_TOOL_LIBRARIES` is set, libomp falls back to loading a tool named
  `libarcher.so` (`openmp/runtime/src/ompt-general.cpp`).  The name is hardcoded,
  so on macOS, where the library is called `libarcher.dylib`, the automatic
  fallback never finds it.  **Always set `OMP_TOOL_LIBRARIES` on macOS.**
- OMPT support is only compiled into libomp on supported platforms.  On
  aarch64/Apple it is on by default, and the conda `libomp.dylib` used by ARTS
  does have it (check with `nm -gU .../lib/libomp.dylib | grep -c ompt`).

## Prerequisites

- A Clang that supports `-fsanitize=thread` and an OpenMP source tree of the
  **same major version** (clang 23 -> `release/23.x`).
- `cmake` >= 3.20, `ninja`, `git` >= 2.25 (for sparse checkouts).
- About 100 MB of scratch disk space and a few minutes of wall time (almost all
  of it is the checkout; the build itself takes seconds).

```sh
# Scratch directory for source, build tree and install prefix
export TSAN_ROOT=$HOME/tsan-archer
mkdir -p "$TSAN_ROOT"
```

## Step 1: Get the OpenMP sources

The legacy standalone build (`cmake <llvm-project>/openmp`) was removed in
LLVM 17; current runtimes are built from `<llvm-project>/runtimes`.

```sh
cd "$TSAN_ROOT"
git clone --depth 1 --branch release/23.x --filter=blob:none --sparse \
    https://github.com/llvm/llvm-project.git llvm-project-23
cd llvm-project-23
git sparse-checkout set openmp runtimes cmake llvm/cmake third-party llvm/utils
```

`third-party` (gtest) and `llvm/utils` (llvm-lit) are needed by the CMake
configure step; without them configuration fails with an error about
`third-party/unittest` not being an existing directory.  A full
`git clone --depth 1` also works if disk space is not a concern.

## Step 2: Build libomp and Archer

```sh
cd "$TSAN_ROOT"
cmake -G Ninja -S llvm-project-23/runtimes -B build \
  -DLLVM_ENABLE_RUNTIMES=openmp \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=$(which clang) \
  -DCMAKE_CXX_COMPILER=$(which clang++) \
  -DLIBOMP_OMPT_SUPPORT=ON \
  -DOPENMP_ENABLE_OMPT_TOOLS=ON \
  -DLIBOMP_ARCHER_SUPPORT=ON \
  -DCMAKE_INSTALL_PREFIX="$TSAN_ROOT/install"
cmake --build build -j
cmake --build build --target install
```

Notes:

- `-DOPENMP_ENABLE_OMPT_TOOLS=ON` is required: the OpenMP build system defaults
  it to `OFF` on Apple hosts ("tools are not tested well on Windows or MacOS X").
  `LIBOMP_ARCHER_SUPPORT` already defaults to `ON`; passing it is just explicit.
- Warnings during configure about missing `FileCheck` / disabled `check-*`
  targets are harmless (they only mean the OpenMP test suite is unavailable).
- The results are:
  - `build/openmp/tools/archer/libarcher.dylib`
  - `build/openmp/runtime/src/libomp.dylib` (OMPT-enabled)
  - `install/lib/libarcher.dylib`, `install/lib/libomp.dylib`
- If your existing `libomp.dylib` already has OMPT support (the conda one does),
  you only need the tool, and can build it alone with
  `cmake --build build --target archer`.

## Step 3: Verify the tool

Two tiny programs distinguish the barrier false positive from a real race.
Compile each with `clang -fsanitize=thread -fopenmp -O1 FILE -o OUT`; the TSan
runtime is linked into the executable, so no preloading is needed here.

A sequentially correct program that TSan mis-reports without Archer
(worker writes in region A, worker reads in region B):

```c
/* test_cross.c */
#include <omp.h>
int x = 0, y = 0;
int main(void) {
#pragma omp parallel
  {
    if (omp_get_thread_num() == 1) x = 1; /* worker 1 writes - region A */
  }
#pragma omp parallel
  {
    if (omp_get_thread_num() == 2) y = x; /* worker 2 reads  - region B */
  }
  return y;
}
```

A genuinely racy program, which must still be reported:

```c
/* test_real.c */
#include <omp.h>
int y[4];
int main(void) {
#pragma omp parallel for
  for (int i = 0; i < 100000; i++) y[i % 4]++; /* real race, no reduction */
  return y[0];
}
```

Run them with `OMP_NUM_THREADS=4`:

| Program    | without Archer               | with Archer              |
|------------|------------------------------|--------------------------|
| test_cross | 1 false-positive race report | clean                    |
| test_real  | race report                  | race report (not hidden) |

```sh
export OMP_NUM_THREADS=4
export OMP_TOOL_LIBRARIES=$TSAN_ROOT/install/lib/libarcher.dylib

./test_cross # no output: correct
./test_real  # still reports: SUMMARY: ThreadSanitizer: data race ...
```

A TSan abort exits with code 134; a clean run exits with the program's own
return code.  With `ARCHER_OPTIONS=verbose=1` the tool prints
`Archer detected OpenMP application with TSan, supplying OpenMP
synchronization semantics` on stderr when it attaches.

## Step 4: Run the ARTS test suite with the tool

Environment must be set in the shell that launches `ctest`; ctest forwards it to
the test process (which already gets `DYLD_INSERT_LIBRARIES` from the test
properties created by the `ENABLE_TSAN` build).

```sh
export OMP_TOOL_LIBRARIES=$TSAN_ROOT/install/lib/libarcher.dylib
# optional: export ARCHER_OPTIONS=verbose=1

# from the ARTS source directory (with the TSan build in build/)
ctest --test-dir build -R 'pyarts.arts.tests.core.vmr.abs$' -V
```

Expected: the test passes and TSan prints no `WARNING: ThreadSanitizer` lines.
For reference, without Archer the same command fails with
`ThreadSanitizer: reported 60 warnings` and `Subprocess aborted`.

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `TSAN functionality requested but not available` | You used the removed `-DLIBOMP_TSAN_SUPPORT` option.  Drop it and use Archer instead. |
| Configure fails at `third-party/unittest` or `llvm/utils/llvm-lit` | Widen the sparse checkout (see Step 1). |
| Archer does not attach (no message with `ARCHER_OPTIONS=verbose=1`) | `OMP_TOOL_LIBRARIES` not visible to the test process, wrong path, or `libomp` without OMPT.  Check `nm -gU $(clang -print-file-name=libomp.dylib) \| grep -c ompt`; if that is 0, use the `libomp.dylib` built in Step 2 via `DYLD_LIBRARY_PATH` (consulted before `@rpath`; only possible for non-SIP-restricted binaries such as the conda Python). |
| Still many warnings, some involving uninstrumented modules | Archer suggests `TSAN_OPTIONS="ignore_noninstrumented_modules=1"`.  Note that this suppresses *all* races where one side is in an uninstrumented module - including Python itself in pyarts runs - so treat it as a last resort. |
| Warnings remain that do not look like barrier artifacts | Take them seriously: they are candidates for real races, unless an uninstrumented runtime is involved. |
| Archer segfaults or behaves oddly | Archer on macOS is not tested upstream.  Fall back to a `-DNO_OPENMP=1` TSan build (see `README.md`), which removes OpenMP altogether, or run with `OMP_NUM_THREADS=1` to silence OpenMP while keeping the rest of the code under TSan. |

Useful Archer options (`ARCHER_OPTIONS="opt=value,..."`): `verbose=1`,
`disabled=1`, `ignore_serial=1`, `tasking=1`, `print_max_rss=1`.

## Reference

- Verified with: macOS arm64, clang 23.1.1 (conda), CMake 4.4.3, Ninja 1.13.2,
  LLVM `release/23.x` (commit `0cd2817`), conda `libomp.dylib` used unchanged.
- Before: `ctest ... -R pyarts.arts.tests.core.vmr.abs$` ->
  `ThreadSanitizer: reported 60 warnings`, `Subprocess aborted`.
- After (with Archer): 0 warnings, `100% tests passed`.
- LLVM OpenMP build documentation:
  <https://openmp.llvm.org/Building.html>
