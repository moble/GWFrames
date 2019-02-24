FROM continuumio/miniconda2

# Install gcc and g++
RUN apt-get update --fix-missing && \
    apt-get install -y --reinstall build-essential && \
    apt-get clean

# Get conda itself set up with our requirements
RUN /opt/conda/bin/conda update -y -n base -c defaults conda && \
    /opt/conda/bin/conda install -y -q -n base swig==3.0.10 numpy gsl hdf5 fftw h5py matplotlib ipython jupyter && \
    /opt/conda/bin/conda clean -y -a

# Get the code, compile, and install
RUN git clone https://github.com/moble/GWFrames.git && \
    cd GWFrames && \
    git submodule init && \
    git submodule update && \
    cd Code && \
    CPATH=/opt/conda/include/ LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:/opt/conda/lib/ /opt/conda/bin/python setup.py install
