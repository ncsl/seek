% REFERENCE: http://www.fieldtriptoolbox.org/tutorial/human_ecog
%% Setup Paths and Directories for Registration Results and Raw EEG Data
% ------- CHANGE THIS ------
addpath("/home/ksr/Downloads/fieldtrip-20190615"); %fieldtrip package
%20190129/");
% addpath("/Users/adam2392/Dropbox/fieldtrip-20181108");
%addpath("/Users/adam2392/Dropbox/fieldtrip-20190129");
addpath("/home/ksr/Downloads/spm12") % SPM12 pacakge 
%addpath("/Users/adam2392/Documents/MATLAB/spm12");
% addpath("../../.lib/fieldtrip-20181108");

% ------- CHANGE THIS ------
% neuroimaging output data dir
RESULTS_DIR = '/home/ksr/Desktop/contact_localization/freesurfer_output';
% RESULTS_DIR = '/Users/adam2392/Dropbox/phd_research/data/neuroimaging_results/freesurfer_output/outputfiles';
% RESULTS_DIR = '/Users/adam2392/Downloads/neuroimaging_results/';
% RESULTS_DIR = '/home/adam2392/hdd/data/neuroimaging/freesurfer_output/";

% ------- CHANGE THIS ------
% raw data EEG dir
RECORD_DATADIR = '/home/ksr/Desktop/contact_localization/raw';%'/home/adam2392/hdd/data/rawdata';
% RECORD_DATADIR = '/Users/adam2392/Downloads/tngpipeline/';

% ------- CHANGE THIS ------
% subj id to analyze
subjID = 'efri06';      % e.g. efri01, pt01, umf001, tvb1
center = 'efri';         % e.g. cc, efri, nih, etc.
%modality = 'ieeg';
%Reading electrode names
%Here reading from eeg files
%dataset_id = strcat('SEIZURE.edf')
%dataset_id = strcat('Seizure.edf')
data_id= [upper(subjID) '_WAR_SES1_Setup.mat'];
%Here reading from saved setup files
%% Initialize Variablees and Parameters
% run setup of global variables to make tool GUI work
ft_defaults

% filepath to your raw data (the header of the edf file has the channel
% labels and stuff)
rawdatafilepath = fullfile(RECORD_DATADIR, subjID, data_id);

% results directory
SUBJDIR = fullfile(RESULTS_DIR, subjID);
CTDIR = fullfile(SUBJDIR, 'CT');
MRIDIR = fullfile(SUBJDIR, 'mri');

% output filepath for your electrode localization
elec_coords_filepath = fullfile(CTDIR, [subjID '_elec_initialize.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE EITHER COREGISTERED CT FILE, OR ORIGINAL CT IMAGE FILE.
% NAME OUTPUT ACCORDINGLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ctimgfile = fullfile(strcat('CT_IN_T1_mgz.nii.gz'))
ctimgfile = fullfile(strcat('CT.nii'));
t1imgfile = fullfile(strcat('T1.nii'));

%% Read in Original MRI Scan
t1_img = ft_read_mri(fullfile(MRIDIR, t1imgfile));

%% Read in Original CT Scan
ct_img = ft_read_mri(fullfile(CTDIR, ctimgfile));

%% Electrode Placement - Runs GUI
%EEG labels
%hdr = ft_read_header(rawdatafilepath);

% configuration for your GUI to read
cfg         = [];
cfg.channel =elec_name;%Labels from analyzed data not eeg. If from eeg file hdr.label;
elecf = ft_electrodeplacement(cfg, ct_img);

% print out the structure that stores your electrode contact points - xyz
elecf

% save it into .mat file - name it accordingly
save(elec_coords_filepath, 'elecf');