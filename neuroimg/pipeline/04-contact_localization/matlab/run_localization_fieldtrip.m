clear all
close all
clc

% REFERENCE: http://www.fieldtriptoolbox.org/tutorial/human_ecog
%% Setup Paths and Directories for Registration Results and Raw EEG Data
% ------- CHANGE THIS ------
addpath("/home/ksr/Downloads/fieldtrip-20190615"); %fieldtrip package
addpath("/home/ksr/Downloads/spm12") % SPM12 pacakge 
% addpath("/Users/adam2392/Dropbox/fieldtrip-20181108");
addpath("/Users/adam2392/Dropbox/fieldtrip-20190129");

% ------- CHANGE THIS ------
% neuroimaging output data dir
studyname = 'epilepsy'; % 'efri'
RESULTS_DIR = fullfile('/home/ksr/Desktop/contact_localization/freesurfer_output/', studyname);
% RESULTS_DIR = '/home/adam2392/hdd/data/neuroimaging/freesurfer_output/";
RESULTS_DIR = fullfile('/Users/adam2392/Dropbox/epilepsy_bids/derivatives/freesurfer/');

% ------- CHANGE THIS ------
% raw data EEG dir
RECORD_DATADIR = '/home/ksr/Desktop/contact_localization/freesurfer_output/epilepsy/'; % '/home/ksr/Desktop/contact_localization/freesurfer_output/epilepsy/';

% ------- CHANGE THIS ------
% subj id to analyze
subjID = 'nl12';      % e.g. efri01, pt01, umf001, tvb1
your_initials = 'CH';

% Reading electrode names, or reading from eeg files -> get channel names 
% data_id= [upper(subjID) '_WAR_SES1_Setup.mat']; 
data_id = [subjID '_sz_6p_chanlabels.mat'];
chandata = load(fullfile(RECORD_DATADIR, data_id)); % load(fullfile('/home/ksr/Desktop/contact_localization/raw/', data_id));
elec_names = chandata.chanlabels;
elec_labels = cellstr(elec_names);
% for i=1:size(elec_names,2)
%     elec_labels{i} = convertCharsToStrings(elec_names(:, i)); %, '\W', '');
% end
% elec_names = num2cell(elec_names, 1);
% size(elec_names, 1), size(elec_names, 2));

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
elec_coords_filepath = fullfile(SUBJDIR, 'elecs', [subjID your_initials '_elec_initialize.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE EITHER COREGISTERED CT FILE, OR ORIGINAL CT IMAGE FILE.
% NAME OUTPUT ACCORDINGLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ctimgfile = fullfile(strcat('CT_IN_T1_mgz.nii.gz'))
ctimgfile = fullfile(strcat('CT.nii'));
% t1imgfile = fullfile(strcat('T1.nii'));

%% Read in Original MRI Scan
% t1_img = ft_read_mri(fullfile(MRIDIR, t1imgfile));

%% Read in Original CT Scan
ct_img = ft_read_mri(fullfile(CTDIR, ctimgfile));

%% Electrode Placement - Runs GUI
% EEG labels read from an eeg file
% rawdatafilepath = './la04_ictal.edf';
% hdr = ft_read_header(rawdatafilepath);

% configuration for your GUI to read
cfg         = [];
cfg.channel = elec_labels; % Labels from analyzed data not eeg. If from eeg file hdr.label;
elecf = ft_electrodeplacement(cfg, ct_img);

% print out the structure that stores your electrode contact points - xyz
elecf

% save it into .mat file - name it accordingly
save(elec_coords_filepath, 'elecf');