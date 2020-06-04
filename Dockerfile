FROM ubuntu:18.04
FROM python:latest
# SHELL ["/bin/bash", "-c"]
LABEL authors="Christopher Coogan <c.coogan2201@gmail.com>, Adam Li <ali39@jhu.edu>"

############# SYSTEM LEVEL INSTALLS #############
# install basic ubuntu utilities
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    python \
    apt-utils \
    libglu1-mesa \
    git \
    g++ \
    gcc \
    libxmu-dev \
    libxi6 \       
    libgconf-2-4 \
    libfontconfig1 \
    libxrender1 \
    octave \ 
    gawk \
    unzip \
    curl \
    libxss1 \
    libjpeg62 \
    tcsh \
    bc \
    dialog

# Freesurfer
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
COPY data_examples/scripts/freesurferlicense.txt /usr/local/freesurfer/.license
ENV SUBJECTS_DIR=/data/derivatives/freesurfer

# acpc detect
RUN mkdir /usr/local/art
ENV ARTHOME /usr/local/art
COPY data_examples/scripts/acpcdetect_v2.0_LinuxCentOS6.7.tar.gz /usr/local/art
RUN tar -xvzf /usr/local/art/acpcdetect_v2.0_LinuxCentOS6.7.tar.gz --no-same-owner -C $ARTHOME
RUN rm /usr/local/art/acpcdetect_v2.0_LinuxCentOS6.7.tar.gz
# doesn't work yet cuz we need to wget from a login page... :(
#RUN wget -O- https://www.nitrc.org/frs/download.php/10595/acpcdetect_v2.0_LinuxCentOS6.7.tar.gz |
#    tar -xvzf --no-same-owner -C $ARTHOME
ENV PATH $ARTHOME/bin:$PATH

# Neuroimgpipe dependencies
RUN pip3 install --upgrade pip
RUN pip3 install snakemake mne-bids numpy scipy mne dicom2nifti

# Node & bids-validator
RUN curl -sL https://deb.nodesource.com/setup_13.x | bash - && \
    apt-get install -y nodejs && \
    npm i -g bids-validator

############# KEEP BELOW SYSTEM LEVEL INSTALLS #############
# setup working directories
WORKDIR /seek

# set environment variable for where analysis takes place
ENV SEEKHOME /seek

# copy over data files
COPY data_examples /data
COPY ./seek/pipeline/01-prep /seek/pipeline/01-prep
COPY ./seek/pipeline/02-reconstruction /seek/pipeline/02-reconstruction
COPY ./seek/pipeline/03-coregistration /seek/pipeline/03-coregistration

# copy over files and functions
COPY ./seek/pipeline/fileutils.py /seek/pipeline/fileutils.py
COPY ./seek/format /seek/format
COPY ./seek/pipeline/config/localconfig.yaml /seek/pipeline/config/localconfig.yaml
