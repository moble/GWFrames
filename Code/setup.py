#! /usr/bin/env python

"""
Installation file for python code associated with the paper "Angular
velocity of gravitational radiation from precessing binaries and the
corotating frame".

To build this code and run it in place, run
    python setup.py build_ext --inplace
then open python and type 'import GWFrames' in
the current directory.

To install in the user's directory, run
    python setup.py install --user
Now, 'import GWFrames' may be run from a python
instance started in any directory on the system.
"""

from os.path import isdir, isfile, exists, abspath, join
from subprocess import check_call
from sys import argv, executable
from distutils.sysconfig import get_python_lib
from os import makedirs, pardir, devnull, environ, getcwd
from shutil import copyfile
from distutils.command.build import build
from distutils.command.build_ext import build_ext as _build_ext
from distutils.command.install_lib import install_lib as _install_lib
from distutils.command.install_scripts import install_scripts as _install_scripts
from distutils.file_util import copy_file
from distutils.core import setup, Extension
from subprocess import check_output, CalledProcessError
from numpy import get_include
import glob

# Make sure submoduls are checked out
if not (isdir('Quaternions') and isfile('Quaternions/setup.py')):
    raise EnvironmentError("Can't find `Quaternions` module.  Did you forget to `git submodule init` and `git submodule update`?")
if not (isdir('SphericalFunctions') and isfile('SphericalFunctions/setup.py')):
    raise EnvironmentError("Can't find `SphericalFunctions` module.  Did you forget to `git submodule init` and `git submodule update`?")
if not (isdir('SpacetimeAlgebra') and isfile('SpacetimeAlgebra/setup.py')):
    raise EnvironmentError("Can't find `SpacetimeAlgebra` module.  Did you forget to `git submodule init` and `git submodule update`?")
if not (isdir('PostNewtonian/C++') and isfile('PostNewtonian/C++/PNEvolution_Q.cpp')):
    raise EnvironmentError("Can't find `PostNewtonian/C++` module.  Did you forget to `git submodule init` and `git submodule update`?")

## Need to build the spinsfast.so first, in its own setup, because
## it's all C code, whereas GWFrames is C++ code.  But we will just
## include the object files below, when we link the _GWFrames.so
print("\nBuilding spinsfast first")
cmd = 'cd spinsfast && {0} python/setup.py build --build-temp build/temp install --install-lib lib/'.format(executable)
print(cmd)
check_call(cmd, shell=True)
print("Finished building spinsfast")

## Build Quaternions and SphericalFunctions
print("\nInstalling SphericalFunctions")
cmd = ' '.join(['cd SphericalFunctions && ', executable]+argv)
print(cmd)
check_call(cmd, shell=True)
print("Finished installing SphericalFunctions\n")

## If PRD won't let me keep a subdirectory, make one
if not exists('GWFrames') :
    makedirs('GWFrames')
if not exists('GWFrames/plot.py') :
    copyfile('plot.py', 'GWFrames/plot.py')

## distutils doesn't build swig modules in the correct order by
## default -- the python module is installed first.  This will pop
## 'build_ext' to the beginning of the command list.
build.sub_commands = sorted(build.sub_commands, key=lambda sub_command: int(sub_command[0]!='build_ext'))

## We also need to copy the SWIG-generated python script GWFrames.py
## to GWFrames/__init__.py so that it gets installed correctly.
class build_ext(_build_ext):
    """Specialized Python source builder for moving SWIG module."""
    def run(self):
        _build_ext.run(self)
        copy_file('SWIG/GWFrames.py', 'GWFrames/__init__.py')
class install_lib(_install_lib):
    """Specialized Python lib installer for moving spinsfast.so."""
    def run(self):
        from distutils.sysconfig import get_python_lib
        import os
        ## Also include a hack to install in the correct directory if necessary
        if not os.access(self.install_dir, os.W_OK):
            if('--user' in argv):
                import site
                self.install_dir = site.USER_SITE
            else:
                self.install_dir = get_python_lib()
        _install_lib.run(self)
class install_scripts(_install_scripts):
    """Hack to install scripts in the correct directory if necessary"""
    def run(self):
        import os
        if not os.access(self.install_dir, os.W_OK):
            if('--user' in argv):
                import site
                self.install_dir = site.USER_BASE+'/bin'
            else:
                from distutils.sysconfig import get_python_lib;
                self.install_dir = abspath(join(get_python_lib(), pardir, pardir, pardir, 'bin'))
        _install_scripts.run(self)

# Add directories for numpy and other inclusions
IncDirs = ['spinsfast/include',
           'SphericalFunctions',
           'Quaternions',
           'SpacetimeAlgebra',
           get_include()]
LibDirs = []

## See if GSL_HOME is set; if so, use it
if "GSL_HOME" in environ :
    IncDirs += [environ["GSL_HOME"]+'/include']
    LibDirs += [environ["GSL_HOME"]+'/lib']

## See if FFTW3_HOME is set; if so, use it
if "FFTW3_HOME" in environ :
    IncDirs += [environ["FFTW3_HOME"]+'/include']
    LibDirs += [environ["FFTW3_HOME"]+'/lib']

# If /opt/local directories exist, use them
if isdir('/opt/local/include'):
    IncDirs += ['/opt/local/include']
if isdir('/opt/local/lib'):
    LibDirs += ['/opt/local/lib']

## Remove a compiler flag that doesn't belong there for C++
import distutils.sysconfig as ds
cfs=ds.get_config_vars()
for key, value in cfs.items() :
    if(type(cfs[key])==str) :
        cfs[key] = value.replace('-Wstrict-prototypes', '').replace('-Wshorten-64-to-32', '')

## Try to determine an automatic version number for this
try :
    with open(devnull, "w") as fnull :
        GitRev = check_output('git rev-parse HEAD', shell=True, stderr=fnull)[:-1]
        CodeRevision = '"{0}"'.format(GitRev)
        PackageVersion = GitRev[:9]
except (NameError, CalledProcessError) :
    CodeRevision = '"PaperVersion3"'
    PackageVersion = '3'

## Read in the license
try :
    with open('LICENSE', 'r') as myfile :
        License=myfile.read()
except IOError :
    License = 'See LICENSE file in the source code for details.'

swig_opts=['-globals', 'constants', '-c++', '-builtin', '-outdir', 'SWIG/', '-DUSE_GSL']
## Add -py3 to swig_opts if this is python3
try:
    import sys
    python_major = sys.version_info.major
    if(python_major==3) :
        swig_opts += ['-py3']
except AttributeError:
    print("Your version of python is SO old.  'How old is it?'  So old I can't even tell how old it is.")
    print("No, seriously.  You should think about upgrading your python because I don't support this version.")
    print("You can try to make this run by removing the assertion error you're about to get, but don't")
    print("come crying to me when print statements fail or when division gives the wrong answer.")
    raise AssertionError("Wake up grandpa!  You were dreaming of ancient pythons again.")

## This does the actual work
setup(name="GWFrames",
      # version=PackageVersion,
      description='Angular velocity of gravitational radiation from precessing binaries and the corotating frame',
      long_description="""
      This package implements various methods described in the paper
      "Angular velocity of gravitational radiation from precessing
      binaries and the corotating frame".  In particular, it
      implements three useful classes -- Quaternion, Waveform, and
      PNWaveform -- along with various methods for each.""",
      author='Michael Boyle',
      author_email='boyle@astro.cornell.edu',
      url='https://github.com/MOBle',
      license=License,
      packages = ['GWFrames'],
      # py_modules = ['GWFrames'],
      scripts = ['Scripts/RunExtrapolations.py', 'Scripts/ExtrapolateAnnex.py'],
      ext_modules = [
        Extension('_GWFrames',
                  sources = ['Quaternions/Quaternions.cpp',
                             'Quaternions/IntegrateAngularVelocity.cpp',
                             'Quaternions/QuaternionUtilities.cpp',
                             'PostNewtonian/C++/PNEvolution.cpp',
                             'PostNewtonian/C++/PNEvolution_Q.cpp',
                             'PostNewtonian/C++/PNWaveformModes.cpp',
                             'SphericalFunctions/Combinatorics.cpp',
                             'SphericalFunctions/WignerDMatrices.cpp',
                             'SphericalFunctions/SWSHs.cpp',
                             'SpacetimeAlgebra/SpacetimeAlgebra.cpp',
                             'Utilities.cpp',
                             'Waveforms.cpp',
                             'PNWaveforms.cpp',
                             'WaveformsAtAPointFT.cpp',
                             'fft.cpp',
                             'NoiseCurves.cpp',
                             'Interpolate.cpp',
                             'Scri.cpp',
                             'SWIG/GWFrames.i'],
                  depends = ['Quaternions/Quaternions.hpp',
                             'Quaternions/IntegrateAngularVelocity.hpp',
                             'Quaternions/QuaternionUtilities.hpp',
                             'PostNewtonian/C++/PNEvolution.hpp',
                             'PostNewtonian/C++/PNWaveformModes.hpp',
                             'SphericalFunctions/Combinatorics.hpp',
                             'SphericalFunctions/WignerDMatrices.hpp',
                             'SphericalFunctions/SWSHs.hpp',
                             'SpacetimeAlgebra/SpacetimeAlgebra.hpp',
                             'Utilities.hpp',
                             'Waveforms.hpp',
                             'PNWaveforms.hpp',
                             'WaveformsAtAPointFT.hpp',
                             'fft.hpp',
                             'NoiseCurves.hpp',
                             'Interpolate.hpp',
                             'Scri.hpp',
                             'Errors.hpp',
                             'GWFrames_Doc.i'],
                  include_dirs=IncDirs,
                  library_dirs=LibDirs,
                  libraries=['gsl', 'gslcblas', 'fftw3'],
                  define_macros = [('CodeRevision', CodeRevision)],
                  language='c++',
                  swig_opts=swig_opts,
                  extra_objects = glob.glob('spinsfast/build/temp/*/*.o'),
                  extra_link_args = ['-fPIC',],
                  # extra_link_args=['-lgomp', '-fPIC', '-Wl,-undefined,error'], # `-undefined,error` tells the linker to fail on undefined symbols
                  extra_compile_args=['-Wno-deprecated', '-Wno-unused-variable', '-DUSE_GSL'] #'-fopenmp',
                  # extra_compile_args=['-ffast-math'] # DON'T USE fast-math!!!  It makes it impossible to detect NANs
                  ),
        ],
      # classifiers = ,
      # distclass = ,
      # script_name = ,
      # script_args = ,
      # options = ,
      # license = ,
      # keywords = ,
      # platforms = ,
      # cmdclass = ,
      cmdclass={'build_ext': build_ext, 'install_lib': install_lib, 'install_scripts': install_scripts},
      # data_files = ,
      # package_dir =
      )
