# Dockerfile for the TerraRef hyperspectral image conversion extractor
FROM ubuntu:16.04
MAINTAINER Max Burnette <mburnet2@illinois.edu>

# Install any programs needed
RUN useradd -u 49044 extractor \
    && mkdir /home/extractor \
    && chown extractor /home/extractor

RUN apt-get update && apt-get upgrade -y -q \
    && apt-get install -y  --no-install-recommends \
        curl \
        bzip2 \
        python3 \
        python3-pip \
        git \
        libpng-dev \
        libfreetype6-dev \
        python-matplotlib \
        libudunits2-dev

# install conda
USER extractor
RUN cd ~ \
    && curl https://repo.anaconda.com/archive/Anaconda2-5.3.1-Linux-x86_64.sh > Anaconda2-5.3.1-Linux-x86_64.sh \
    && bash Anaconda2-5.3.1-Linux-x86_64.sh -b

#install conda-forge packages
RUN ~/anaconda2/bin/conda config --add channels conda-forge \
    && ~/anaconda2/bin/conda install -y \
         libnetcdf \
         hdf5 \
         netcdf4 \
         nco \
         "gdal>2.2.4" \
         libiconv \
         xerces-c \
         geos \
         udunits2

ENV PATH="/home/extractor/anaconda2/bin:${PATH}"

RUN pip install --user pyclowder terrautils laspy spectral pika pandas scikit-learn h5py

USER root
RUN apt-get update && apt-get upgrade -y -q \
    && apt-get install -y  --no-install-recommends vim libnetcdf-dev
USER extractor

COPY logging_config.json /var/log/
COPY *.sh extractor_info.json *.py *.nc *.nco *.sh /home/extractor/
COPY calibration_new /home/extractor/calibration_new
COPY s9_boundary.csv /home/extractor/s9_boundary.csv

RUN pip install --user cryptography

ENTRYPOINT ["/home/extractor/entrypoint.sh"]
CMD ["extractor"]

# Setup environment variables. These are passed into the container. You can change
# these to your setup. If RABBITMQ_URI is not set, it will try and use the rabbitmq
# server that is linked into the container. MAIN_SCRIPT is set to the script to be
# executed by entrypoint.sh
ENV RABBITMQ_EXCHANGE="terra" \
    RABBITMQ_VHOST="%2F" \
    RABBITMQ_QUEUE="terra.hyperspectral" \
    MAIN_SCRIPT="terra_hyperspectral.py" \
    CLOWDER_SPACE="5bdc8f174f0cb2fdaaf3148e" \
    UDUNITS2_XML_PATH="/usr/share/xml/udunits/udunits2.xml" \
    HDF5_USE_FILE_LOCKING="FALSE"