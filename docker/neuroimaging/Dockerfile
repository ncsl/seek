FROM ubuntu:18.04
# FROM freesurfer/freesurfer:6.0
MAINTAINER Adam Li <ali39@jhu.edu>

#############################################################################
############# FOR CREATING A FREESURFER, FSL, MRTRIX3 Pipeline ##############
#############################################################################

# Core system capabilities required
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    nano \
    nodejs \
    npm \
    perl-modules \
    python \
    tar \
    unzip \
    wget \
    libeigen3-dev \
    zlib1g-dev \
    libqt4-opengl-dev \
    libgl1-mesa-dev \
    libfftw3-dev \
    libtiff5-dev \
    g++ \
    python-numpy 

RUN apt-get install git g++ python python-numpy libeigen3-dev zlib1g-dev libqt4-opengl-dev libgl1-mesa-dev libfftw3-dev libtiff5-dev
RUN apt-get upgrade -y

# ANTs installs tzdata as a dependency; however its installer is interactive
# Therefore we need to do some shenanigans here to force it though
RUN DEBIAN_FRONTEND=noninteractive \
    apt-get install -y tzdata

# NeuroDebian setup
RUN wget -qO- http://neuro.debian.net/lists/bionic.us-nh.libre | \
    tee /etc/apt/sources.list.d/neurodebian.sources.list

# Silences unnecessary warning regarding stdout not being a terminal
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1
RUN apt-key adv --recv-keys --keyserver pgp.mit.edu 2649A5A9
RUN apt-get update

###################### For FreeSurfer Setup: ######################

RUN wget -qO- https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/6.0.0/freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.0.tar.gz | \
    tar zx -C /opt \
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
    --exclude='freesurfer/lib/cuda' \
    --exclude='freesurfer/lib/qt'
RUN apt-get install -y ants

RUN apt-get update && ls

# COPY fslinstaller.py /
RUN wget -q https://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py -O /fslinstaller.py
RUN python ./fslinstaller.py -d /opt/fsl
# RUN apt-get install -y fsl-5.0-complete
# RUN apt-get install -y fsl-first-data
# RUN apt-get install -y fsl-mni152-templates

COPY license.txt /opt/freesurfer/.license

# Make ANTS happy
ENV ANTSPATH=/usr/lib/ants
ENV PATH=/usr/lib/ants:$PATH

# Make Freesurfer happy
ENV PATH=/opt/freesurfer/bin:/opt/freesurfer/mni/bin:$PATH
ENV OS Linux
ENV SUBJECTS_DIR /opt/freesurfer/subjects
ENV FSF_OUTPUT_FORMAT nii.gz
ENV MNI_DIR /opt/freesurfer/mni
ENV LOCAL_DIR /opt/freesurfer/local
ENV FREESURFER_HOME /opt/freesurfer
ENV FSFAST_HOME /opt/freesurfer/fsfast
ENV MINC_BIN_DIR /opt/freesurfer/mni/bin
ENV MINC_LIB_DIR /opt/freesurfer/mni/lib
ENV MNI_DATAPATH /opt/freesurfer/mni/data
ENV PERL5LIB /opt/freesurfer/mni/lib/perl5/5.8.5
ENV MNI_PERL5LIB /opt/freesurfer/mni/lib/perl5/5.8.5

###################### For FSL Setup: ######################

#install fsl
RUN apt-get update 

# Make FSL happy
ENV PATH=/usr/lib/fsl/6.0:/opt/eddy:$PATH
ENV FSLDIR=/usr/share/fsl/6.0
ENV FSLMULTIFILEQUIT=TRUE
ENV FSLOUTPUTTYPE=NIFTI
ENV LD_LIBRARY_PATH=/usr/lib/fsl/6.0
# RUN chmod 775 /opt/eddy/eddy_openmp

###################### For MRTRIX3 Setup: ######################

RUN git clone https://github.com/MRtrix3/mrtrix3.git mrtrix3 
RUN cd mrtrix3 && \
	./build clean && \
	git pull && \
	./configure -static -nogui
RUN cd mrtrix3 && \
	./build

# git checkout 3.0_RC3 && \
# ./build clean \
# ./configure -nogui && NUMBER_OF_PROCESSORS=1 \
# git pull \
# ./configure -static [-nogui] \
# ./build 
# git describe --tags > /mrtrix3_version

# Setup environment variables for MRtrix3
ENV PATH=/mrtrix3/bin:$PATH
ENV PYTHONPATH=/mrtrix3/lib:$PYTHONPATH

# remove unnecessary files copied over
RUN rm fslinstaller.py
COPY version /version
