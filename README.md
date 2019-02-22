| NOTE: This module is no longer *actively* maintained.  Much of the functionality has been moved to [the `scri` package](https://github.com/moble/scri).  Pull requests are welcome, and I'll generally reply to issues.  But because of the age of this code and its dependencies (especially SWIG), I may not be able to help very much if you run into problems.  —Mike |
| --- |


GWFrames
========
Manipulate gravitational waveforms—changing frames, and so on.

The code in this project extends code written for the paper "[Angular
velocity of gravitational radiation from precessing binaries and the
corotating frame](http://arxiv.org/abs/1302.2919)", giving an explicit
implementation of the methods discussed in it.  The first commit of
this project is the last commit of the paper's project.

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


Quick Start
===========
The easiest way to get GWFrames running is to
[install anaconda](https://www.anaconda.com/distribution/)
and run the following commands
```bash
# Update conda itself first
conda update -y -n base -c defaults conda

# Create a new conda env just for GWFrames
# You may want to add "matplotlib ipython jupyter" to the end of the next line
conda create -y -n GWFrames python=2 swig==3.0.10 numpy gsl hdf5 fftw h5py
conda activate GWFrames

# Get the code
git clone https://github.com/moble/GWFrames.git
cd GWFrames
git submodule init
git submodule update

# Compile and install the code into this conda env
cd Code
python setup.py install

# Test that it works
cd ~
python -c 'import GWFrames; print("That worked!")'
```
This gives you a separate `GWFrames` environment with working code, and
without polluting other parts of your more modern python installation.


Build requirements
==================
To build just the C++ code:
* C++ compiler
* [GNU Scientific Library](http://www.gnu.org/software/gsl/), built as a shared library

To use the optional—but highly recommended—Python interface:
* [SWIG](http://www.swig.org/) v3.0.10
* [HDF5](http://www.hdfgroup.org/HDF5/) v1.8.10 or greater, built as a shared library with development headers (libhdf5-dev or similar)
* [FFTW](http://www.fftw.org/) v3.2 or greater, built as a shared library

And, of course, python and a few of its goodies, for which I cannot
recommend [anaconda](http://continuum.io/downloads) highly enough --
it makes installation and maintenance of the python software stack
vastly easier.
* [Python](http://www.python.org/getit/) v2.7.4 or greater (but less than v3), with development headers (python-dev or similar)
* [NumPy](http://www.numpy.org/) v1.7 or greater (`conda install numpy`)
* [Matplotlib](http://matplotlib.org/) v1.2 or greater (`conda install matplotlib`)
* [h5py](http://code.google.com/p/h5py/) v2.1 or greater (`conda install h5py`)
* [IPython](http://ipython.org/) v2.0 or greater, with notebook (`conda install ipython`)

All of the above are reasonably standard, and can be installed easily
through package managers such as apt and homebrew.  The notebook
interface for IPython provides an environment much like the
Mathematica notebook interface.  The examples are provided in a
notebook (`.ipynb`) file.


Getting started
===============
The first step is to clone the git repo and its submodules:
```
git clone https://github.com/moble/GWFrames.git
cd GWFrames
git submodule init
git submodule update
```

The code is mostly in the form of C++, for speed.  All functions are
provided in the `GWFrames` namespace.  The interface is reasonably
straightforward, and can easily be included into and compiled with
other code.  Look for examples of how to do this in the `C++Example`
directory.  (See below for a few more details.)

However, it is mostly intended to be used through the interface to
python.  With the above-listed requirements, build the code and
install to the user directory.  To do this, run `make` or

    python setup.py install --user

After this, you should be able to start a new python session (in any
directory) and `import GWFrames`, which is the basic module provided
by this code.  The primary objects are `Quaternion` and `Waveform`
objects, with various methods defined for each.  An IPython notebook
in the `Docs` directory contains extensive examples, following the
outline of the paper.  To use it, run

    ipython notebook Docs/AngVelOfGravRadFromPrecessingBinaries_CorotatingFrame.ipynb --pylab

Alternatively, this code may be used directly through its C++
interface.  A simple example is provided in the `C++Example`
directory, along with the Makefiles needed to build all the necessary
code.  If it does not compile easily, make sure the various paths in
_both_ Makefiles are set properly.

Detailed documentation of most functions may be found through python's
`help` function, or by running `make` in the `Docs` subdirectory, and
reading `Docs/html/index.html`.


Contributions
=============
The code is primarily written by me (Mike Boyle).  Dan Hemberger
helped by porting some of my older code from
[Triton](https://github.com/moble/Triton) to perform noise-weighted
overlap calculations, and with numerous bug reports and helpful
suggestions.  Serguei Ossokine also helped substantially by
cross-checking the post-Newtonian formulas and results.  Dante Iozzo
has done his best to keep things working as this code enters its dotage.

Other contributions are entirely welcome.  The preferred method is via
github's excellent interface.  If you have a bug report, just go to
the [issues](https://github.com/moble/GWFrames/issues) page for this
repo, check for related issues, and if this is new click "New Issue".
If you want to contribute code, it is easiest for everyone if you use
the [fork & pull method described here]
(https://help.github.com/articles/using-pull-requests).
