from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib
import os

numpy_inc = os.path.join(get_python_lib(plat_specific=1), 'numpy/core/include')

module1 = Extension('spinsfast',
                    sources = ['python/spinsfast_module.c'],
                    include_dirs = [numpy_inc,'include'],
                    libraries=['spinsfast','fftw3'],
                    library_dirs = ["lib"],
                    extra_compile_args=['-std=c99','-fPIC'],
                    depends=['lib/libspinsfast.a']
                    )

setup (name = 'spinsfast',
       version = '104',
       description = 'Fast, exact spin-s spherical harmonic transforms',
       ext_modules = [module1],
       )
