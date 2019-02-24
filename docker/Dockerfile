FROM continuumio/miniconda2

# Get conda itself set up with our requirements
RUN /opt/conda/bin/conda update -y -n base -c defaults conda && \
    /opt/conda/bin/conda install -y -q -n base swig==3.0.10 numpy gsl hdf5 fftw h5py matplotlib ipython jupyter && \
    /opt/conda/bin/conda clean -y -a

# # Get the code
# git clone https://github.com/moble/GWFrames.git
# cd GWFrames
# git submodule init
# git submodule update

# Compile and install the code
RUN cd Code && \
    /opt/conda/bin/python setup.py install
