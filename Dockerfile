# 18.04 or 20.04+ should work
FROM ubuntu:latest
#FROM freesurfer/freesurfer:7.1.1
ENV DEBIAN_FRONTEND noninteractive

ARG FSL_VERSION=6.0.4

# SHELL ["/bin/bash", "-c"]
LABEL authors="Christopher Coogan <c.coogan2201@gmail.com>, Adam Li <ali39@jhu.edu>"

# Install dependencies
RUN apt-get update && apt-get -y install bc binutils libgomp1 perl psmisc sudo tar tcsh unzip \
    uuid-dev vim-common curl octave gawk jq locales apt-utils
#    libjpeg62-turbo-dev
RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    locale-gen
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en

# Download Blender and set path
RUN mkdir /usr/local/blender \
    && curl -SL https://download.blender.org/release/Blender2.82/blender-2.82a-linux64.tar.xz -o blender.tar.xz \
    && tar -xf blender.tar.xz -C /usr/local/blender --strip-components=1 \
    && rm blender.tar.xz \
    && curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
    && /usr/local/blender/2.82/python/bin/python3.7m get-pip.py \
    && /usr/local/blender/2.82/python/bin/python3.7m -m pip install pandas

############# SYSTEM LEVEL INSTALLS #############
# install basic ubuntu utilities
RUN apt-get update && apt-get install -y \
    build-essential wget \
    python apt-utils \
    libglu1-mesa \
    git g++ gcc \
    libxmu-dev libxi6 libgconf-2-4 \
    libfontconfig1 libxrender1 \
    octave gawk \
    unzip curl \
    libxss1 libjpeg62 \
    tcsh bc dialog

# acpc detect
RUN mkdir /usr/local/art
ENV ARTHOME /usr/local/art
COPY .data/acpcdetect_v2.0_LinuxCentOS6.7.tar.gz /usr/local/art
RUN tar -xvzf /usr/local/art/acpcdetect_v2.0_LinuxCentOS6.7.tar.gz --no-same-owner -C $ARTHOME
RUN rm /usr/local/art/acpcdetect_v2.0_LinuxCentOS6.7.tar.gz
# doesn't work yet cuz we need to wget from a login page... :(
#RUN wget -O- https://www.nitrc.org/frs/download.php/10595/acpcdetect_v2.0_LinuxCentOS6.7.tar.gz |
#    tar -xvzf --no-same-owner -C $ARTHOME
ENV PATH $ARTHOME/bin:$PATH

# Freesurfer built locally
RUN wget -O- https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/dev/freesurfer-linux-centos7_x86_64-dev.tar.gz | \
    tar zx --no-same-owner -C /usr/local/ \
    --exclude='freesurfer/trctrain' \
    --exclude='freesurfer/subjects/fsaverage_sym' \
    --exclude='freesurfer/subjects/fsaverage3' \
    --exclude='freesurfer/subjects/fsaverage4' \
    --exclude='freesurfer/subjects/fsaverage5' \
    --exclude='freesurfer/subjects/fsaverage6' \
    --exclude='freesurfer/subjects/cvs_avg35' \
    --exclude='freesurfer/subjects/cvs_avg35_inMNI152' \
    --exclude='freesurfer/subjects/bert' \
    --exclude='freesurfer/subjects/V1_average' \
    --exclude='freesurfer/average/mult-comp-cor' \
    --exclude='freesurfer/lib/cuda' 
ENV PATH=/usr/local/freesurfer/bin:/usr/local/freesurfer/mni/bin:$PATH
ENV FREESURFER_HOME /usr/local/freesurfer
RUN $FREESURFER_HOME/SetUpFreeSurfer.sh
COPY ./freesurferlicense.txt /usr/local/freesurfer/.license
ENV SUBJECTS_DIR=/data/derivatives/freesurfer

# FSL builder

# Node & bids-validator
RUN curl -sL https://deb.nodesource.com/setup_13.x | bash - && \
    apt-get install -y nodejs && \
    npm i -g bids-validator

# Seek dependencies as a python virtual environment
FROM python:3.8
RUN pip install --upgrade pip
RUN pip install pipenv

# copy over Pipfiles
COPY Pipfile* /tmp/

# Copy locking -> requirements.txt and install using pip
#RUN cd /tmp && pipenv lock --requirements > requirements.txt
#RUN pip3 install -r /tmp/requirements.txt

# or install by copying entire directory
#COPY . /tmp/myapp
#RUN pip3 install /tmp/myapp

# install via pipenv
RUN cd /tmp && pipenv install --skip-lock --system

############# KEEP BELOW SYSTEM LEVEL INSTALLS #############
# setup working directories
WORKDIR /seek

# set environment variable for where analysis takes place
ENV SEEKHOME /seek

# copy over data files
COPY data /data

# copy over code and workflows
COPY ./seek/ /seek/
COPY ./workflow/ /workflow/
COPY ./config/ /config/
