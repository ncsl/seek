FROM ubuntu:18.04
SHELL ["/bin/bash", "-c"]
LABEL authors="Christopher Coogan <c.coogan2201@gmail.com>, Adam Li <ali39@jhu.edu>"
WORKDIR /home

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
    curl \ 
    gawk \
    unzip

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
COPY ./docker/license.txt /usr/local/freesurfer/.license
ENV SUBJECTS_DIR=/home/data

# Blender
RUN mkdir /usr/local/blender
RUN wget https://download.blender.org/release/Blender2.82/blender-2.82a-linux64.tar.xz -O /home/blender.tar.xz
RUN tar -xf /home/blender.tar.xz -C /usr/local/blender
ENV PATH /usr/local/blender:$PATH

# FSL
RUN wget -q https://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py -O /home/fslinstaller.py
RUN python /home/fslinstaller.py -d /opt/fsl
ENV FSLDIR=/opt/fsl
ENV USER=me
RUN source $FSLDIR/etc/fslconf/fsl.sh

# Neuroimg_pipeline
COPY ./neuroimg /home/neuroimg

# SPM
RUN wget https://www.fil.ion.ucl.ac.uk/spm/download/restricted/eldorado/spm12.zip -O /home/spm12.zip
RUN unzip -q /home/spm12.zip -d /home/spm12

# Fieldtrip
RUN wget ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/fieldtrip-20200302.zip -O /home/fieldtrip

# Matlab
COPY ./docker/matlab.zip /home/matlab.zip
RUN unzip -q /home/matlab.zip -d /home/matlab
RUN /home/matlab/install

# Remove temporary files
RUN rm blender.tar.xz fslinstaller.py matlab.zip spm12.zip anaconda.sh

# Display variable for X11
ENV DISPLAY=192.168.1.156:0.0

# Webserver
    # Expose viz/webserver to the host
EXPOSE 5000/tcp
# Copy completed files from recon folder to static webserver folder (move this to a snakemake rule)

# Set matlab paths