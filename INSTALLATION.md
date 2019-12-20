# INSTALLATION GUIDE

There are a couple of tools that you need to install in your system before everything is working. You ar recommended to use a Linux based OS. 
Follow links and tutorials on each respective tool to install. Preferably this is done via Docker, or Singularity, but if not, then:

Anaconda and Python3.6+ :
   * Conda (https://docs.anaconda.com/anaconda/install/)
   * This is mainly necessary to run img_pipe (ECoG localization with Chang Lab repo), snakemake, and any Python wrapper code
    
    i. Conda env
    
    
        # create environment
        conda create -n neuroimgpipe
        conda activate neuroimgpipe
        
        # optionally separate install
        conda config --add channels bioconda
        conda config --add channels conda-forge
        conda install numpy scipy matplotlib scikit-learn scikit-image pandas seaborn nibabel mne snakemake mne-bids flask
        conda install pytest black check-manifest pytest-cov pydocstyle
        conda env export > environment.yml
        
        # check if installation worked
        cd neuroimg/pipeline/reconstruction/
        snakemake -n    
   

## 3rd Party Modules to Install 
0. Octave
    Runs open-source. This runs various scripts for converting output files to object files for rendering visualizations.
    Follow: https://www.gnu.org/software/octave/#install
    
   
    brew install octave
       
0. Gawk
    Runs command line tools.
    https://brewinstall.org/Install-gawk-on-Mac-with-Brew/

0. Blender
    https://www.blender.org/download/Blender2.81/blender-2.81-linux-glibc217-x86_64.tar.bz2/

1. Reconstruction
    * Freesurfer (https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall)
    * This step is necessary to generate a parcellation and surface reconstruction of the patient's brain. The general requirements is just a 
    Linux, or OSX computer with enough RAM. Currently, this repo is designed to work with FreeSurfer.
    
2. Coregistration
    * FSL Flirt (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/)
    * This step is necessary to map different imaging sessions together. Specifically, for this pipeline, we need it to map CT images to T1 MRI
    * Note that as of 2019, installation still requires Python2, which should come in any Linux distribution.
        
            python2 <run_installer>
    
3. Utility
    * MRTrix3 (https://mrtrix.readthedocs.io/en/latest/installation/linux_install.html)

4. SPM 
    * SPM install (preferably 12): https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

5. Contact-Localization Software (FieldTripToolbox, Img_Pipe, MATLAB)
    * FieldTripToolbox (http://www.fieldtriptoolbox.org/download/)
    * Img_Pipe from the Chang Lab at UCSF will come as a submodule in this git repo. This heavily handles ECoG data only.
 
6. (Optional) Cloud Reconstruction (MRICLOUD):
    * MRICloud (cloud based soln; just send images here) (https://mricloud.org/)
    * the nice thing is that this usually works even when FS fails (e.g. the T1 MRI image isn't good enough quality, or there is a major lesion, etc.).

7. (Optional) ACPC Auto Detection:
    * https://www.nitrc.org/projects/art/
    
8. (Optional) Nonlinear Registration NDREG:
    * NDReg (https://github.com/neurodata/ndreg)
