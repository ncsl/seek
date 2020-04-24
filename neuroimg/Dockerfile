FROM ubuntu:18.04
FROM python:latest
# SHELL ["/bin/bash", "-c"]
LABEL authors="Christopher Coogan <c.coogan2201@gmail.com>, Adam Li <ali39@jhu.edu>"

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
    bc

WORKDIR /neuroimg
COPY ./Data /data
COPY ./neuroimg/pipeline/01-prep /neuroimg/pipeline/01-prep
COPY ./neuroimg/pipeline/02-reconstruction /neuroimg/pipeline/02-reconstruction
COPY ./neuroimg/pipeline/03-coregistration /neuroimg/pipeline/03-coregistration

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
COPY ./neuroimg/freesurferlicense.txt /usr/local/freesurfer/.license
ENV SUBJECTS_DIR=/data/derivatives/freesurfer


# Neuroimgpipe dependencies
RUN pip3 install snakemake mne-bids numpy scipy mne dicom2nifti 

# Node & bids-validator
RUN curl -sL https://deb.nodesource.com/setup_13.x | bash - && apt-get install -y nodejs && npm i -g bids-validator


COPY ./neuroimg/pipeline/fileutils.py /neuroimg/pipeline/fileutils.py
COPY ./neuroimg/format /neuroimg/format
COPY ./neuroimg/pipeline/config/localconfig.yaml /neuroimg/pipeline/config/localconfig.yaml
