% REFERENCE: http://www.fieldtriptoolbox.org/tutorial/human_ecog
%% Setup Paths and Directories for Registration Results and Raw EEG Data
addpath("/home/adam2392/Documents/Dropbox/fieldtrip-20181108");
addpath("/Users/adam2392/Dropbox/fieldtrip-20181108");
addpath("/Users/adam2392/Documents/MATLAB/spm12");
% addpath("../../.lib/fieldtrip-20181108");

% run setup of global variables to make tool GUI work
ft_defaults

% neuroimaging output data dir
RESULTS_DIR = '/home/adam2392/hdd/data/neuroimaging/freesurfer_output/outputfiles/';
% RESULTS_DIR = '/Users/adam2392/Dropbox/phd_research/data/neuroimaging_results/freesurfer_output/outputfiles';
% RESULTS_DIR = '/Users/adam2392/Downloads/neuroimaging_results/';
% RESULTS_DIR = '/home/adam2392/hdd/data/neuroimaging/freesurfer_output/";

% raw data EEG dir
RECORD_DATADIR = '/home/adam2392/hdd/data/rawdata';
% RECORD_DATADIR = '/Users/adam2392/Downloads/tngpipeline/';

% filepath to your raw data (the header of the edf file has the channel
% labels and stuff)
rawdatafilepath = fullfile(RECORD_DATADIR, 'cleveland/', subjID, '/seeg/edf/', strcat(subjID, '_ictal.edf'))

% subj id to analyze
subjID = 'la03';

% results directory
SUBJDIR = fullfile(RESULTS_DIR, subjID);

% output filepath for your electrode localization
elec_coords_filepath = [subjID '_elec_f.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE EITHER COREGISTERED CT FILE, OR ORIGINAL CT IMAGE FILE.
% NAME OUTPUT ACCORDINGLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ctimgfile = fullfile(strcat('CT_IN_T1_mgz.nii.gz'))
ctimgfile = fullfile(strcat('CT.nii.gz'))
t1imgfile = fullfile(strcat('T1.nii.gz'))

%% Read in Original MRI Scan
t1_img = ft_read_mri(fullfile(SUBJDIR, t1imgfile));

%% Read in Original CT Scan
% read in CT in original format
ct_img = ft_read_mri(fullfile(SUBJDIR, ctimgfile));

%% Electrode Placement - Runs GUI
hdr = ft_read_header(rawdatafilepath);

% configuration for your GUI to read
cfg         = [];
cfg.channel = hdr.label;
elecf = ft_electrodeplacement(cfg, ct_img);

% print out the structure that stores your electrode contact points - xyz
elecf

% save it into .mat file - name it accordingly
save(elec_coords_filepath, 'elecf');