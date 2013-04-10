# NOTE: This Makefile calls ../Makefile.  If compilation is not
# working, check the options in ../Makefile also.

# The compiler needs to be able to find the GSL (GNU Scientific
# Library) headers and libraries.  The following paths are the most
# common places for these to be installed.  If compilation doesn't
# work, correct these paths.
INCFLAGS = -I/opt/local/include -I/usr/local/include
LIBFLAGS = -L/opt/local/lib -L/usr/local/lib

# Set compiler name and optimization flags here, if desired
C++ = g++
OPT = -O3 -Wall -fopenmp # Don't use -ffast-math
LINK_OPENMP = -lgomp

#############################################################################
## The following are pretty standard and probably won't need to be changed ##
#############################################################################

# Tell 'make' not to look for files with the following names
.PHONY : all clean allclean realclean sources

# Default target calls the targets listed here
all : sources main

# Call 'make cpp' in the directory above this one
sources :
	$(MAKE) -C .. cpp

# Compile main.cpp into an executable
main : sources main.cpp
	$(C++) $(OPT) -I.. $(INCFLAGS) $(LIBFLAGS) $(LINK_OPENMP) -lgsl -lgslcblas ../Utilities.o ../Quaternions.o ../Waveforms.o main.cpp -o main

# The following are just handy targets for removing compiled stuff
clean :
	-/bin/rm -f main
allclean : clean
	@$(MAKE) -C .. clean
realclean : allclean
