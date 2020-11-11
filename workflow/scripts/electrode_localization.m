% electrode_localization.m
%
% Script that performs electrode localization using FieldTrip (ft) toolbox as described in:
%
% Stolk 2018 - Integrated analysis of anatomical and electrophysiological human intracranial data.pdf
% http://www.fieldtriptoolbox.org/tutorial/human_ecog/
%
%
% Macauley Breault, Adam Li, Kristin Gunnarsdottir
% Created:  10-09-2019
clear; clc;
% BIDS root of the data
bids_root = 'C:\Users\vanik\Johns Hopkins\Adam Li - epilepsy_bids\';  % change this 
deriv_path = fullfile(bids_root, 'derivatives');
source_path = fullfile(bids_root, 'sourcedata');

% add fieldtrip toolbox
addpath(fullfile(bids_root, 'sample_scripts', 'contact_localization', 'fieldtrip-20190129'));
run('initialize_fieldtrip.m')

save_it  = 1; % Boolean to save results
plot_it  = 1; % Boolean to plot
place_it = 1; % Boolean to implement placement GUI
align_mri = 0; % Boolean to manually align the MRI or not?

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Specification of subject ID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Step 1
% dataset BIDS entities
subjID = 'la04';  % change this for all patients
sessionID = 'presurgery';
space = 'fs';
extension = '.nii';
datatype = 'anat';

% setup the BIDS paths according to entities
bids_sub = ['sub-', subjID];
bids_ses = ['ses-', sessionID];

% setup session path
ses_path = fullfile(bids_root, bids_sub, bids_ses);

% setup file names for the different images
mri_bids_basename = [bids_sub, '_', ...
                bids_ses, '_', ...
                'space-' , space , '_' , 'T1w' , extension];
ct_bids_basename = [bids_sub , '_'  , ...
                bids_ses , '_' , ...
                'space-' , 'orig' , '_' , 'CT' , extension];
stolk_ct_bids_basename = [bids_sub , '_'  , ...
                bids_ses , '_' , ...
                'space-' , space , '_' , ... 
                'proc-' , 'spm12' , '_' , 'CT' , extension];

% use this if you want the raw dicoms
% mri_fname = fullfile(namePath(mfilename('fullpath'),fullfile(save_path,'premri')),'0000.dcm')
mri_fname = fullfile(ses_path, 'anat', mri_bids_basename);
ct_fname = fullfile(ses_path, 'ct', ct_bids_basename);

% the output filename for the CT file after transforming to ACPC
% out_ct_fname = fullfile(save_path,'stolk',[subjID,'_CT_acpc_f']);
out_ct_fname = fullfile(ses_path, 'ct', stolk_ct_bids_basename);

% raw EEG filepath
% setup_elecs_fname = fullfile(save_path,listing(cellfun(@(name) contains(name,'_Setup.mat'),{listing.name})).name);
setup_elecs_fname = fullfile(source_path, 'electrodes localized', 'setup', ['sub-', subjID, '_setup.mat']);

% output electrodes mat file
elecs_mat_fpath = fullfile(source_path, 'electrodes localized', 'stolk', [subjID '_elec_acpc_f.mat']);
%should I add this? save_path = namePath(mfilename('fullpath'),fullfile('electrode_localization','raw',subjID));

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Preprocessing of the anatomical MRI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
close all force

% Step 2
mri_fs = ft_read_mri(mri_fname);


if align_mri
    % Step 3
    ft_determine_coordsys(mri_fs)
    %{
    Box 3. Step 2. X-axis
    Box 3. Step 3. right-to-left orientation
    Determine the orientation of the left?right axis.
    If the values on the left?right axis increase to the right (indicated by a ?+? sign), then the scan has a left-to-right orientation.
    If the values on the left?right axis increase to the left, then the scan has a right-to-left orientation.
    %}
    
    % Step 4

    % Do you want to manually align the T1 MRI to ACPC? If already done,
    % do not need to run this step
    cfg = [];
    cfg.coordsys = 'acpc';
    cfg.method = 'fiducial';
    switch subjID
    %     case 'efri02'
    %         % EFRI02
    %         cfg.fiducial.ac = [260 256 138];
    %         cfg.fiducial.pc = [260 278 113];
    %         cfg.fiducial.xzpoint = [260  83 100];
    %         cfg.fiducial.right = [160 253 100];
        otherwise
            cfg.method = 'interactive';
    end
        
    mri_acpc = ft_volumerealign(cfg, mri_fs);
    
    %{
    1. To change the slice viewed in one plane, either:
       a. click (left mouse) in the image on a different plane. Eg, to view a more
          superior slice in the horizontal plane, click on a superior position in the
          coronal plane, or
       b. use the arrow keys to increase or decrease the slice number by one
    2. To mark a fiducial position or anatomical landmark, do BOTH:
       a. select the position by clicking on it in any slice with the left mouse button
       b. identify it by pressing the letter corresponding to the fiducial/landmark:
          press a for ac, p for pc, z for xzpoint
          press r for an extra control point that should be on the right side
       You can mark the fiducials multiple times, until you are satisfied with the positions.
    3. To change the display:
       a. press c on keyboard to toggle crosshair visibility
       b. press f on keyboard to toggle fiducial visibility
       c. press + or - on (numeric) keyboard to change the color range's upper limit
    4. To finalize markers and quit interactive mode, press q on keyboard
    contrast limits updated to [0.100 0.500]

    % For ac and pc reference, use:
     http://www.fieldtriptoolbox.org/assets/img/faq/anterior_commissure/acpcline.png
    %}

    % Step 5
    if save_it
        cfg = [];
        cfg.filename = fullfile(save_path,'stolk',[subjID,'_MR_acpc']);
        cfg.filetype = 'nifti';
        cfg.parameter = 'anatomy';
        ft_volumewrite(cfg, mri_acpc);
    end
else
    % assume the read in T1 MRI is already acpc and FreeSurfer aligned
    mri_acpc = mri_fs;
end

fsmri_acpc = mri_acpc;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Preprocessing of the anatomical CT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
close all force

% Step 9
disp('ON STEP 9 - Reading in CT image!');
ct = ft_read_mri(ct_fname);

% Step 10
disp('ON STEP 10 - Determining Coordinate system for CT image!');
ft_determine_coordsys(ct)
%{
Box 3. Step 2. X-axis
Box 3. Step 3. right-to-left orientation
%}

% Step 11
disp('ON STEP 11 - Re-aligning CT image as best as possible to ACPC!');
cfg = [];
cfg.coordsys = 'ctf';
cfg.method = 'fiducial';
switch subjID
%     case 'efri02'
%         % % EFRI02
%         cfg.fiducial.nas = [271 101   3];
%         cfg.fiducial.lpa = [401 290   3];
%         cfg.fiducial.rpa = [116 267   7];
%         cfg.fiducial.zpoint = [252 290 133];
    otherwise
        cfg.method = 'interactive';
end

ct_ctf = ft_volumerealign(cfg, ct);
%{
1. To change the slice viewed in one plane, either:
   a. click (left mouse) in the image on a different plane. Eg, to view a more
      superior slice in the horizontal plane, click on a superior position in the
      coronal plane, or
   b. use the arrow keys to increase or decrease the slice number by one
2. To mark a fiducial position or anatomical landmark, do BOTH:
   a. select the position by clicking on it in any slice with the left mouse button
   b. identify it by pressing the letter corresponding to the fiducial/landmark:
      press n for nas, l for lpa, r for rpa
      press z for an extra control point that should have a positive z-value
   You can mark the fiducials multiple times, until you are satisfied with the positions.
3. To change the display:
   a. press c on keyboard to toggle crosshair visibility
   b. press f on keyboard to toggle fiducial visibility
   c. press + or - on (numeric) keyboard to change the color range's upper limit
4. To finalize markers and quit interactive mode, press q on keyboard

 crosshair: voxel  35755737, index = [217 204 137], head = [-17.8 -40.8 64.0] mm
       nas: voxel   7632644, index = [260  60  30], head = [2.3 -80.1 -56.3] mm
       lpa: voxel   5647253, index = [405 278  22], head = [70.3 21.0 -39.1] mm
       rpa: voxel   2749522, index = [ 82 251  11], head = [-81.1 11.4 -52.9] mm
    zpoint: voxel  35755737, index = [217 204 137], head = [-17.8 -40.8 64.0] mm
%}

% Step 12
disp('ON STEP 12 - Converting coordsystem of CT image to approx ACPC!');
ct_acpc_f = ft_convert_coordsys(ct_ctf, 'acpc', 0);

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fusion of the CT with the MRI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
close all force

% Step 13
cfg = [];
cfg.method = 'spm';
cfg.spmversion = 'spm12';
cfg.coordsys = 'acpc';
if plot_it
    cfg.viewresult = 'yes';
else
    cfg.viewresult = 'no';
end
ct_acpc_f = ft_volumerealign(cfg, ct_acpc_f, fsmri_acpc);

% Step 15
if save_it
    cfg           = [];
    cfg.filename  = out_ct_fname;
    cfg.filetype  = 'nifti';
    cfg.parameter = 'anatomy';
    ft_volumewrite(cfg, ct_acpc_f);
end

if plot_it
    pause
end

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Electrode placement ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%

close all force

% Step 16
warning('Warning in %s: Replace my file with electrode names with your file of electrode names here',mfilename)
% listing = dir(save_path);
stp = load(setup_elecs_fname);

% EEG labels read from an eeg file
% stp = ft_read_header(rawdatafilepath);

cfg         = [];

% Step 17
if place_it
    cfg = [];
    cfg.channel = stp.elec_name;
%     cfg.channel = stp.label;   %  If from eeg file hdr.label;

    if exist(elecs_mat_fpath,'file')
        load(elecs_mat_fpath, 'elec_acpc_f')
        cfg.elec = elec_acpc_f;
    end
    
    % run electrode placement
    elec_acpc_f = ft_electrodeplacement(cfg, ct_acpc_f, fsmri_acpc);
    
    if save_it
        if exist(elecs_mat_fpath,'file')
            error('Error in %s: File already exists. Please save by hand if you are sure.',mfilename)
        else
            save(elecs_mat_fpath, 'elec_acpc_f')
        end
    end
end


% Step 18
load(elecs_mat_fpath, 'elec_acpc_f')

% Step 19
if plot_it
    ft_plot_ortho(fsmri_acpc.anatomy, 'transform', fsmri_acpc.transform, 'style', 'intersect');
    ft_plot_sens(elec_acpc_f, 'label', 'on', 'fontcolor', 'w','facecolor',[1 1 1]);
    pause
end


