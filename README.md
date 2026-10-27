> [!NOTE]
> **For post-Newtonian waveforms, use
> [PostNewtonian.jl](https://github.com/moble/PostNewtonian.jl).** It
> is far more complete, accurate, and flexible than the post-Newtonian
> code in this package.  Most of the rest of the functionality here
> has moved to [`scri`](https://github.com/moble/scri) and
> [`sxs`](https://github.com/sxs-collaboration/sxs).  GWFrames is kept
> working mainly as a reference implementation to compare other codes
> against.

[![Tests](https://github.com/moble/GWFrames/actions/workflows/tests.yml/badge.svg)](https://github.com/moble/GWFrames/actions/workflows/tests.yml)
[![Wheels](https://github.com/moble/GWFrames/actions/workflows/wheels.yml/badge.svg)](https://github.com/moble/GWFrames/actions/workflows/wheels.yml)

GWFrames
========

Manipulate gravitational waveforms — changing frames, and so on.

The code in this project extends code written for the paper "[Angular
velocity of gravitational radiation from precessing binaries and the
corotating frame](https://arxiv.org/abs/1302.2919)", giving an
explicit implementation of the methods discussed in it.  The first
commit of this project is the last commit of the paper's project.

The license for using this software is basically open (see the LICENSE
file in this directory), though citations to the original paper are
appreciated, where relevant.  Also, if your work depends on features
found in this module that have not been described in a separate
publication of mine, I would appreciate the opportunity to be a
coauthor.

Note that this code is not optimized.  In many cases, obvious
optimizations were rejected in favor of clearer or simpler code, or
code that simply reflected the discussion in the paper more directly.
The code is generally more than fast enough for interactive use on
individual waveforms, but there is plenty of room for improvement.  In
particular, rotations may be painfully slow for large data sets.


Installation
============

GWFrames requires Python 3.11 or newer, on Linux or macOS.  Binary
wheels for Linux and macOS (x86_64 and arm64) are published on PyPI.

The HDF5 functions (`ReadFromH5`, `OutputToH5`, `ReadFromNRAR`, and
`OutputToNRAR`) and the `GWFrames.Extrapolation` module need
[h5py](https://www.h5py.org/).  PyPI has no h5py wheels for Python
3.15, and there may never be any, so on Python 3.15 and newer GWFrames
is installed without h5py; everything else works without it.  If you
need those functions, either use Python 3.14 or older with `uv` or
`pip`, or get h5py from conda-forge with `pixi` or `mamba`, which do
build it for Python 3.15.

With `uv` or `pip`
------------------

```bash
uv add GWFrames          # in a uv project
uv pip install GWFrames  # or in any virtual environment
pip install GWFrames     # or with pip
```

With `uv`, adding `--python 3.14` to `uv init` or `uv venv` selects a
Python for which h5py is available.

To install the latest code from GitHub, or on any platform without a
wheel, GWFrames is compiled from source.  That needs a C++ compiler,
[GSL](https://www.gnu.org/software/gsl/), and
[FFTW](https://www.fftw.org/), which are available from most package
managers:

```bash
sudo apt-get install libgsl-dev libfftw3-dev   # Debian or Ubuntu
brew install gsl fftw                          # macOS with Homebrew
```

after which

```bash
uv pip install "git+https://github.com/moble/GWFrames"
```

builds and installs it.  (The submodules are fetched automatically.)
If GSL or FFTW is installed somewhere unusual, set `CMAKE_PREFIX_PATH`
to its prefix.

With `pixi`
-----------

[pixi](https://pixi.sh/) can take NumPy, SciPy, and h5py from
conda-forge and GWFrames itself from PyPI:

```bash
pixi init my-project
cd my-project
pixi add python numpy scipy h5py
pixi add --pypi GWFrames
pixi run python -c 'import GWFrames; print(GWFrames.__version__)'
```

With `mamba` or `conda`
-----------------------

Similarly, the dependencies can come from conda-forge and GWFrames
from PyPI:

```bash
mamba create -n gwframes -c conda-forge python numpy scipy h5py
mamba activate gwframes
python -m pip install GWFrames
```

To build the latest code from GitHub instead, all of the build
requirements are available from conda-forge:

```bash
mamba create -n gwframes -c conda-forge python numpy scipy h5py gsl fftw \
    cxx-compiler cmake ninja pkg-config pybind11 scikit-build-core
mamba activate gwframes
python -m pip install --no-deps --no-build-isolation "git+https://github.com/moble/GWFrames"
```

On Apple Silicon, make sure that the environment is native: `mamba
info` should show `platform : osx-arm64`.  Older instructions for
GWFrames recommended Intel environments (with `CONDA_SUBDIR=osx-64` or
a `subdir: osx-64` setting), which run under Rosetta.  The first build
in such an environment can fail with "Accessing CMake timed out" while
Rosetta translates CMake; running the same command again works.

Checking the installation
-------------------------

```bash
python -c 'import GWFrames; W = GWFrames.PNWaveform("TaylorT4", 0.0, [0., 0., 0.9], [0.0, 0.0, 0.9], 0.01); print("That worked!")'
```


Usage
=====

The primary objects are `Waveform` and `PNWaveform`, with various
methods defined for each; quaternions and functions of arrays of
quaternions are in `GWFrames.Quaternions`.  Detailed documentation of
most functions may be found through Python's `help` function, as in
`help(GWFrames.Waveform)`.  A Jupyter notebook in the `Docs` directory
contains extensive examples, following the outline of the paper.

Changes from the SWIG-based versions
------------------------------------

Version 2 replaces the SWIG bindings with
[pybind11](https://pybind11.readthedocs.io/), and the old `setup.py`
with CMake and
[scikit-build-core](https://scikit-build-core.readthedocs.io/).  The
Python interface is almost entirely unchanged, apart from these
differences:

* The separate `Quaternions` module is now `GWFrames.Quaternions`, so
  `import Quaternions` becomes `from GWFrames import Quaternions`, and
  `GWFrames.Quaternion` is also available.
* Methods that modify a `Waveform` in place and return it (such as
  `TransformToCorotatingFrame`) return the very same Python object,
  including for `PNWaveform`.

Several bugs were fixed along the way, notably that `Quaternion`
objects compared equal exactly when they were *not* equal, and that
`Hybridize` read past the end of an array and could not handle
waveforms in an inertial frame.


Development
===========

Clone the repository with its submodules, then build and test it with
`uv`:

```bash
git clone --recursive https://github.com/moble/GWFrames.git
cd GWFrames
uv sync --no-editable
uv run --no-sync pytest
```

The C++ code is in `cpp/`, the pybind11 bindings in `cpp/bindings/`,
and the Python package in `src/GWFrames/`.  The submodules
`Quaternions`, `SphericalFunctions`, `SpacetimeAlgebra`, and
`PostNewtonian`, along with a vendored copy of
[spinsfast](https://github.com/moble/spinsfast), are in `extern/`.
The tests include a comparison with results recorded from the last
SWIG-based build; see `tests/reference/generate_swig_reference.py`.

Every push to `main` that passes the wheel builds and tests is
released automatically: the version is bumped, tagged, released on
GitHub, and published to PyPI.  The bump is a patch release unless a
commit message in the push contains `#minor` or `#major`; adding
`#dev` makes a development release (such as `2.1.0.dev0`), and a later
plain push completes it.  A head commit message containing `[no
release]` builds and tests without releasing, and `[skip ci]` skips CI
entirely.

Note on licensing: spinsfast, GSL, and FFTW are licensed under the GNU
GPL, and the binary wheels link all three statically, so the wheels as
a whole are covered by the GPL.


Contributions
=============

The code is primarily written by me (Mike Boyle).  Dan Hemberger
helped by porting some of my older code from
[Triton](https://github.com/moble/Triton) to perform noise-weighted
overlap calculations, and with numerous bug reports and helpful
suggestions.  Serguei Ossokine also helped substantially by
cross-checking the post-Newtonian formulas and results.  Dante Iozzo
helped to keep things working as this code entered its dotage.  Claude
(Opus 5.5) resuscitated the package for Python 3.
