Docker container with everything necessary to run the GWFrames code.

Usage
-----

You can download and run this image using the following commands:

    docker pull moble/gwframes
    docker run -i -t moble/gwframes bash

You'll probably want to attach some local directory to the container as a "volume", by adding such a flag to the command line:

    docker run -i -t moble/gwframes --volume=/path/to/data:/data bash

This uses the directory found on your local machine in `/path/to/data`, and makes it available in the docker container as `/data`.  Note that this allows the directory on your local machine to be both read from and written to, so that your output will survive even after you exit the container.

A convenient way to use this code is to start a Jupyter Notebook server and interact with it via your browser:

    docker run -i -t -p 8888:8888 moble/gwframes --volume=/path/to/data:/data bash -c "mkdir /opt/notebooks && /opt/conda/bin/jupyter notebook --notebook-dir=/opt/notebooks --ip='*' --port=8888 --no-browser"

You can then view the Jupyter Notebook by opening `http://localhost:8888` in your browser.
