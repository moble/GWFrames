##
## NOTE: When compiling python, this Makefile is only used to call the
## python build process.  No variables set in this file affect python.
## If the python build is having trouble, look in setup.py to fix it.
##





# By default, just call the python build process
all : spinsfast
	python setup.py install --user





















############################################
## Flags for building only the c++ files: ##
############################################
# The compiler needs to be able to find the GSL (GNU Scientific
# Library) headers and libraries.  The following paths are the most
# common places for these to be installed.  If compilation doesn't
# work, correct these paths.
INCFLAGS = -I/opt/local/include -I/usr/local/include -Ispinsfast/include -IQuaternions
LIBFLAGS = -L/opt/local/lib -L/usr/local/lib -Lspinsfast/lib
ifdef GSL_HOME
	INCFLAGS = -I${GSL_HOME}/include ${INCFLAGS}
	LIBFLAGS = -L${GSL_HOME}/lib ${LIBFLAGS}
endif
# Set compiler name and optimization flags here, if desired
C++ = g++
OPT = -O3 -fopenmp -Wall -Wno-deprecated
## DON'T USE -ffast-math in OPT


#############################################################################
## The following are pretty standard and probably won't need to be changed ##
#############################################################################

# Tell 'make' not to look for files with the following names
.PHONY : all cpp clean allclean realclean swig spinsfast

# If needed, we can also make object files to use in other C++ programs
cpp : Utilities.o Quaternions/Quaternions.o Waveforms.o PNWaveforms.o Scri.o SphericalHarmonics.o

# This is how to build those object files
%.o : %.cpp %.hpp Errors.hpp
	$(C++) $(OPT) -DCodeRevision=3 -c $(INCFLAGS) $< -o $@

# The following are just handy targets for removing compiled stuff
clean :
	-/bin/rm -f *.o
allclean : clean
	-/bin/rm -rf build
	-build=build/config.mk make -C spinsfast clean
realclean : allclean

# This just runs SWIG, which can be handy for committing pre-swigged files
swig :
	swig -python -globals constants -c++ -o GWFrames_wrap.cpp GWFrames.i

# This attempts to build spinsfast
spinsfast :
	build=build/config.mk make -C spinsfast
	cp spinsfast/lib/spinsfast.so `python -m site --user-site`/
