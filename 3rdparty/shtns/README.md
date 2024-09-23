**SHTns is a high performance library for Spherical Harmonic Transform written in C,
aimed at numerical simulation (fluid flows, mhd, ...) in spherical geometries.**

Copyright (c) 2010-2021 Centre National de la Recherche Scientifique.
written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
SHTns is distributed under the open source [CeCILL License](http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html)
(GPL compatible) located in the LICENSE file.

FEATURES:
---------

- **blazingly fast**.
- both **scalar and vector transforms**.
- backward and forward (synthesis and analysis) functions.
- flexible truncation (degree, order, azimuthal periodicity).
- spatial data can be stored in latitude-major or longitude-major arrays.
- various conventions (normalization and Condon-Shortley phase).
- can be used from **Fortran, c/c++, Python, and Java** programs.
- highly efficient Gauss algorithm working with Gauss nodes (based on
  Gauss-Legendre quadrature).
- support for **regular grids** (but they require twice the number of nodes than Gauss grid)
- support for SSE2, SSE3, **AVX, AVX2, AVX-512** vectorization, as well as 
  Xeon Phi (KNL), AltiVec VSX, and Neon.
- **parallel transforms with OpenMP**.
- *Beta*: automatic **GPU off-loading** with nvidia cuda (Kepler, Pascal & Volta).
- synthesis (inverse transform) at any coordinate (not constrained to a grid).
- **on-the-fly transforms** : saving memory and bandwidth, they are even faster
  on modern architectures.
- accurate up to spherical harmonic degree l=16383 (at least).
- **rotation** functions to rotate spherical harmonics.
- special spectral operator functions that do not require a transform
  (multiply by cos(theta)...).
- transforms and rotations for complex-valued spatial fields.
- SHT at fixed m (without fft, aka Legendre transform).


INSTALL:
--------

Briefly, the shell commands `./configure; make; make install` should
configure, build, and install this package. `./configure --help` will
list available options (among which `--enable-openmp` and `--enable-python`).
However, in order to get the best performance, it is highly recommended to
compile and install the FFTW library yourself, because many distributions
include a non-optimized FFTW library.

DOCUMENTATION:
--------------

- On-line doc is available: <https://nschaeff.bitbucket.io/shtns/>
- You can build it locally: Run `make docs` to generate documentation
  (requires doxygen). 
  Then browse the html documentation starting with `doc/html/index.html`
- A related research paper has been published:
  [Efficient Spherical Harmonic Transforms aimed at pseudo-spectral numerical simulations](http://dx.doi.org/10.1002/ggge.20071),
  also [available from arXiv](http://arxiv.org/abs/1202.6522).
- If you use SHTns for research work, please **cite this paper**:

        @article {shtns,
          author = {Schaeffer, Nathanael},
          title = {Efficient spherical harmonic transforms aimed at pseudospectral numerical simulations},
          journal = {Geochemistry, Geophysics, Geosystems},
          doi = {10.1002/ggge.20071}, volume = {14}, number = {3}, pages = {751--758},
          year = {2013},
        }

  If you use Ishioka's recurrence (the default since SHTns v3.4), you may also want to cite his paper:

        @article {ishioka2018,
          author={Ishioka, Keiichi},
          title={A New Recurrence Formula for Efficient Computation of Spherical Harmonic Transform},
          journal={Journal of the Meteorological Society of Japan},
          doi = {10.2151/jmsj.2018-019}, volume={96}, number={2}, pages={241--249},
          year={2018},
        }
