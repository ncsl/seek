FROM ubuntu:18.04
SHELL ["/bin/bash", "-c"]
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
    curl

# Anaconda and Neuroimgpipe's dependencies
RUN wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh -O /home/anaconda.sh
RUN bash /home/anaconda.sh -b -p /usr/local/anaconda
ENV PATH=/usr/local/anaconcda/bin:$PATH
RUN source /usr/local/anaconda/etc/profile.d/conda.sh && \
    conda create -n neuroimgpipe && \
    conda init bash && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install numpy scipy matplotlib scikit-learn scikit-image pandas seaborn nibabel mne snakemake mne-bids flask && \
    conda install pytest black check-manifest pytest-cov pydocstyle
RUN source ~/.bashrc


# Freesurfer
# COPY docker/temp/freesurfer-linux-centos7_x86_64-dev.tar.gz /usr/local/freesurfer.tar.gz
RUN wget -O- https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/dev/freesurfer-linux-centos7_x86_64-dev.tar.gz | \
    tar zx --no-same-owner -C /usr/local/ \
    # RUN tar -C /usr/local -xzvf /usr/local/freesurfer.tar.gz \
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
ENV SUBJECTS_DIR /usr/local/freesurfer/subjects
ENV FREESURFER_HOME /usr/local/freesurfer
RUN $FREESURFER_HOME/SetUpFreeSurfer.sh
COPY docker/license.txt /usr/local/freesurfer/.license

# Blender
RUN mkdir /usr/local/blender
# COPY docker/temp/blender-2.81a-linux-glibc217-x86_64.tar.bz2 /usr/local/blender.tar.bz2
RUN curl -SL https://www.blender.org/download/Blender2.81/blender-2.81a-linux-glibc217-x86_64.tar.bz2/ -o \
    RUN tar -jxvf /usr/local/blender.tar.bz2 -C /usr/local/blender --strip-components=1 && rm /usr/local/blender.tar.bz2
ENV PATH /usr/local/blender:$PATH

# Display variable (freeview, blender, matlab)
ENV DISPLAY=192.168.1.124:0.0

# FSL
# COPY docker/fslinstaller.py /home/fslinstaller.py
RUN wget -q https://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py -O /fslinstaller.py
RUN python /home/fslinstaller.py -d /opt/fsl

ENV FSLDIR=/opt/fsl
ENV USER=me
RUN source $FSLDIR/etc/fslconf/fsl.sh

# Uncomment these when we switch back to wget as opposed to just copying the files over

COPY neuroimg /home/neuroimg

COPY matlab.zip /home/matlab.zip
RUN cd /home && gzip matlab.zip
RUN rm /home/anaconda.sh
RUN rm /home/fslinstaller.py
RUN rm /usr/local/freesurfer.tar.gz
