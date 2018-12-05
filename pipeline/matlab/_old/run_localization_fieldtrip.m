% REFERENCE: http://www.fieldtriptoolbox.org/tutorial/human_ecog

addpath("/home/adam2392/Documents/Dropbox/fieldtrip-20181108");
addpath("/Users/adam2392/Dropbox/fieldtrip-20181108");
addpath("/Users/adam2392/Documents/MATLAB/spm12");
ft_defaults

RESULTS_DIR = '/home/adam2392/hdd/data/neuroimaging/freesurfer_output/outputfiles';
% RESULTS_DIR = '/Users/adam2392/Downloads/neuroimaging_results/';
% RESULTS_DIR = '/home/adam2392/hdd/data/neuroimaging/freesurfer_output/";

RECORD_DATADIR = '/home/adam2392/hdd/data/rawdata';
% RECORD_DATADIR = '/Users/adam2392/Downloads/tngpipeline/';

% subj id
subjID = 'la03';

% results directory
SUBJDIR = fullfile(RESULTS_DIR, subjID);


%% Option 1 - MRI
% % read in mri img
mri_img = fullfile(SUBJDIR, 'testT1.mgz')
mri = ft_read_mri(mri_img);

% reorient mri img
cfg           = [];
cfg.method    = 'interactive';
cfg.coordsys  = 'acpc';
mri_acpc = ft_volumerealign(cfg, mri);

%% Option 2 - MRI
% read in the freesurfer generated mri img
% mri_img = fullfile(SUBJDIR, 'T1.mgz')
% fsmri_acpc = ft_read_mri(mri_img); 
% fsmri_acpc.coordsys = 'ctf';

% read in CT
ct = ft_read_mri(fullfile(SUBJDIR, strcat(subjID, '_CT.nii')));

% realign ct
cfg           = [];
cfg.method    = 'interactive';
cfg.coordsys  = 'ctf';
ct_ctf = ft_volumerealign(cfg, ct);

% convert coord system
ct_acpc = ft_convert_coordsys(ct_ctf, 'ctf');

%% Fuse CT With MRI
cfg             = [];
cfg.method      = 'spm';
cfg.spmversion  = 'spm12';
cfg.coordsys    = 'ctf';
cfg.viewresult  = 'yes';
ct_acpc_f = ft_volumerealign(cfg, ct_acpc, fsmri_acpc);

cfg           = [];
cfg.filename  = [subjID '_CT_acpc_f'];
cfg.filetype  = 'nifti';
cfg.parameter = 'anatomy'; 
ft_volumewrite(cfg, ct_acpc_f);

%% Read in Converted CT In T1 Image 
% % TODO: possibly need to convert this first outside to TAL
% % this means then also the MRI image needs to be converted to TAL
% ct_acpc_f = ft_read_mri(fullfile(SUBJDIR, strcat(subjID,'_CT_IN_T1.nii'))); 
% ct_acpc_f.coordsys = 'acpc';

%% Electrode Placement
hdr = ft_read_header(fullfile(RECORD_DATADIR, 'cleveland/', subjID, '/seeg/edf/', strcat(subjID, '_ictal.edf')));

cfg         = [];
cfg.channel = hdr.label;
elec_acpc_f = ft_electrodeplacement(cfg, ct_acpc_f, fsmri_acpc);

elec_acpc_f


ft_plot_ortho(fsmri_acpc.anatomy, 'transform', fsmri_acpc.transform, 'style', 'intersect');
ft_plot_sens(elec_acpc_f, 'label', 'on', 'fontcolor', 'w');


save([subjID '_elec_acpc_f.mat'], 'elec_acpc_f');

