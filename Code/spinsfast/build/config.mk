INC+=

CFLAGS:=-Wall  -g -std=c99 -O3 -fpic -fstrict-aliasing -Wstrict-aliasing -D_GNU_SOURCE

# Optional portions of code
#CFLAGS+= -DUSE_FITSIO
#CFLAGS+= -DUSE_HEALPIX

CC:=gcc $(CFLAGS)

LIBFFTW=-lfftw3

# Libraries for optional portions
#LIBGSL=-lgsl -lgslcblas
#LIBFITSIO=-lcfitsio
#LIBCHEALPIX=-lchealpix
# In e.g. RHEL/Fedora this is provided by the gsl-devel, cfitsio-devel, and chealpix-devel package