# Dockerfile for the TerraRef hyperspectral image conversion extractor
FROM terraref/terrautils:1.4
MAINTAINER Max Burnette <mburnet2@illinois.edu>

# Install any programs needed
RUN useradd -u 49044 extractor \
    && mkdir /home/extractor \
    && chown extractor /home/extractor

# install common libraries and python modules
RUN apt-get update && apt-get upgrade -y -q && \
    apt-get install -y \
        build-essential \
        libhdf5-dev \
        libnetcdf-dev \
        libgdal-dev \
        python-gdal \
        gdal-bin \
        libudunits2-dev \
        proj-bin \
        # nco \
    && rm -rf /var/lib/apt/lists/*

# install NCO (apt-get version is only 4.6, need 4.7+)
RUN apt-get install -y \
        wget \
        dpkg \
    && wget https://launchpad.net/ubuntu/+archive/primary/+files/libgslcblas0_2.5+dfsg-5_amd64.deb \
    && dpkg -i libgslcblas0_2.5+dfsg-5_amd64.deb \
    && wget https://launchpad.net/ubuntu/+archive/primary/+files/libgsl23_2.5+dfsg-5_amd64.deb \
    && dpkg -i libgsl23_2.5+dfsg-5_amd64.deb \
    && wget http://launchpadlibrarian.net/365856922/libc6_2.27-3ubuntu1_amd64.deb \
    && dpkg -i libc6_2.27-3ubuntu1_amd64.deb \
    && wget http://launchpadlibrarian.net/365926929/gcc-8-base_8-20180414-1ubuntu2_amd64.deb \
    && dpkg -i gcc-8-base_8-20180414-1ubuntu2_amd64.deb \
    && wget http://launchpadlibrarian.net/365927026/libgfortran5_8-20180414-1ubuntu2_amd64.deb \
    && dpkg -i libgfortran5_8-20180414-1ubuntu2_amd64.deb \
    && wget https://launchpad.net/ubuntu/+archive/primary/+files/libhdf5-100_1.10.0-patch1+docs-4build2_amd64.deb \
    && dpkg -i libhdf5-100_1.10.0-patch1+docs-4build2_amd64.deb \
    && wget https://launchpad.net/ubuntu/+archive/primary/+files/libnetcdf13_4.6.1-3_amd64.deb \
    && dpkg -i libnetcdf13_4.6.1-3_amd64.deb \
    && wget https://launchpad.net/ubuntu/+archive/primary/+files/nco_4.7.6-1_amd64.deb \
    && dpkg -i nco_4.7.6-1_amd64.deb

COPY *.sh extractor_info.json *.py *.nc *.nco *.sh /home/extractor/

USER extractor
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
    CLOWDER_SPACE="5bdc8f174f0cb2fdaaf3148e"