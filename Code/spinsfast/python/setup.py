from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib
import os

numpy_inc = os.path.join(get_python_lib(plat_specific=1), 'numpy/core/include')

## The following block is added for nicer behavior with `module`s on clusters
from os import environ
IncDirs = [numpy_inc, "include",]
LibDirs = ["lib",]
## See if GSL_HOME is set; if so, use it
if "GSL_HOME" in environ :
    IncDirs += [environ["GSL_HOME"]+'/include']
    LibDirs += [environ["GSL_HOME"]+'/lib']
## See if FFTW3_HOME is set; if so, use it
if "FFTW3_HOME" in environ :
    IncDirs += [environ["FFTW3_HOME"]+'/include']
    LibDirs += [environ["FFTW3_HOME"]+'/lib']

module1 = Extension('spinsfast',
                    sources = ['python/spinsfast_module.c'],
                    # include_dirs = [numpy_inc,'include'],
                    include_dirs=IncDirs,
                    libraries=['spinsfast','fftw3'],
                    # library_dirs = ["lib"],
                    library_dirs=LibDirs,
                    extra_compile_args=['-std=c99','-fPIC'],
                    depends=['lib/libspinsfast.a'],
                    )

setup (name = 'spinsfast',
       version = '104',
       description = 'Fast, exact spin-s spherical harmonic transforms',
       ext_modules = [module1],
       )
