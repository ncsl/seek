% [electrode_localization.m]
%
% Script that performs electrode localization using FieldTrip (ft) toolbox 
% as described in [1].
%
% All output should be quality-checked.
%
% Required Files
% --------------
% FieldTrip Toolbox     : 
%   The toolbox from Fieldtrip.
% T1w NIFTI file from FreeSurfer    :
%   The T1w image corresponding to ``T1.mgz`` from FreeSurfer (6.0 or 7.x)
%   in NIFTI format. 
% CT NIFTI file 
%   The CT image with iEEG implantations.
% electrodes file
%   Can be the ``*channels.tsv`` file, or a matlab struct with 
%   ``elec_name`` field. 
%
% User Inputs
% -----------
% fieldtrip_toolbox_path    : path
%   The path to the fieldtrip toolbox. Downloadable from FieldTrip website.
% bids_root                 : path
%   The path to the dataset organized according to BIDS. For information 
%   on how to organize your data, see [2].
% subjID                    : str
%   The subject ID.
% sessionID                 : str
%   The session ID for the input/output files. Corresponds to ``session`` 
%   entity defined in BIDS.
% space                     : str
%   The ``space`` entity corresponding to BIDS. Default is 'fs' for 
%   subject-specific FreeSurfer space.
% extension                 : str
%   The file extension of the image files. Default is '.nii' for
%   uncompressed Nifti files.
% setup_elecs_fname         : path
%   The path to the file that will contain the electrode names to be
%   localized on the CT image. See ``electrodes file`` in ``Required
%   Inputs``.
% 
% Outputs
% -------
% elec_acpc_f               : struct
%   The corresponding matlab structure for the electrodes containing
%   ``label``, ``chanpos`` fields for the channel names and xyz channel 
%   positions (in mm) respectively.
% out_ct_fname              : Nifti file
%   The output CT image that is coregistered to the input T1 image. 
%   Is named
%   ``sub-<subjID>_ses-<sessionID>_space-<space>_proc-spm12_CT.nii`` where 
%   a processing entity for ``spm12`` is attached to the filename to depict
%   that coregistration used SPM12 algorithm [3].
% 
% Notes
% -----
% 
% 
% References
% ----------
% [1] Stolk 2018 - Integrated analysis of anatomical and 
% electrophysiological human intracranial data.pdf
% http://www.fieldtriptoolbox.org/tutorial/human_ecog/
% [2] https://bids-specification.readthedocs.io/en/stable/
% [3] SPM12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
% 
% Macauley Breault, Adam Li, Kristin Gunnarsdottir
% Christopher Coogan, Chester Huynh, Raphael Bechtold
% Created:  10-09-2019
clear; clc;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ User specified inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% path to the FieldTrip toolbox
fieldtrip_toolbox_path = fullfile(bids_root, 'sample_scripts', 'contact_localization', 'fieldtrip-20190129')

% BIDS root of the data
bids_root = 'C:\Users\vanik\Johns Hopkins\Adam Li - epilepsy_bids\';  % change this 

% dataset BIDS entities
subjID = 'la04';  % change this for all patients
sessionID = 'presurgery';
space = 'fs';
extension = '.nii';
datatype = 'anat';

% raw EEG filepath
% setup_elecs_fname = fullfile(save_path,listing(cellfun(@(name) contains(name,'_Setup.mat'),{listing.name})).name);
setup_elecs_fname = fullfile(source_path, 'electrodes localized', 'setup', ['sub-', subjID, '_setup.mat']);

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Specification of BIDS file paths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% path to `derivatives` and `sourcedata` folder
deriv_path = fullfile(bids_root, 'derivatives');
source_path = fullfile(bids_root, 'sourcedata');

% add fieldtrip toolbox
addpath(fieldtrip_toolbox_path);

% initializes fieldtrip
if ~contains(path,'fieldtrip-20190129')
    warning('Warning in %s: Replace my path to fieldtrip with your path here. I used fieldtrip-20191008',mfilename)
    fthome = 'C:\Users\vanik\Johns Hopkins\Adam Li - epilepsy_bids\sample_scripts\contact_localization\fieldtrip-20190129';
    path(path,fthome);
    clear fthome;
end
ft_defaults

save_it  = 1; % Boolean to save results
plot_it  = 1; % Boolean to plot
place_it = 1; % Boolean to implement placement GUI
align_mri = 0; % Boolean to manually align the MRI or not?

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Specification of BIDS file paths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Step 1
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

function exportTSV(elec, file)
    % Creates a `.tsv` file from Fieldtrip `cfg.elec` data structure.
    %
    % Meant for saving data in `*electrodes.tsv` BIDS format.
    %
    % Inputs:
    %   elec    : a struct with `label`, and `chanpos` fields.
    %   file    : The filepath to save.

    if ~contains(file, 'electrodes.tsv')
        error(strcat('Filename should be of the form "*electrodes.tsv". ',
        'The filename does not have electrodes.tsv substring inside.'));

    % create separator
    sep = regexp(elec.label,'\d');
    
    % get the vector of x, y, and z coordinates
    xcoord = elec.chanpos(:,1);
    ycoord = elec.chanpos(:,2);
    zcoord = elec.chanpos(:,3);

    % write file to 'name'
    fid = fopen(file,'wt');
    fprintf(fid,'name\tx\ty\tz\n');
    for i = 1:length(elec.label)
        elecName = strcat(elec.label{i}(1:sep{i}-1),"'",elec.label{i}((sep{i}:length(elec.label{i}))));
        fprintf(fid, '%s\t%-.3f\t%-.3f\t%-.3f\n', elecName, xcoord(i), ycoord(i), zcoord(i));
    end
    fclose(fid)
end

