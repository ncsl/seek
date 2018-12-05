% REFERENCE: http://www.fieldtriptoolbox.org/tutorial/human_ecog

%% Setup Paths and Directories for Registration Results and Raw EEG Data
addpath("/home/adam2392/Documents/Dropbox/fieldtrip-20181108");
addpath("/Users/adam2392/Dropbox/fieldtrip-20181108");
addpath("/Users/adam2392/Documents/MATLAB/spm12");
ft_defaults

% neuroimaging output data dir
RESULTS_DIR = '/home/adam2392/hdd/data/neuroimaging/freesurfer_output/outputfiles/';
RESULTS_DIR = '/Users/adam2392/Dropbox/phd_research/data/neuroimaging_results/freesurfer_output/outputfiles';
% RESULTS_DIR = '/Users/adam2392/Downloads/neuroimaging_results/';
% RESULTS_DIR = '/home/adam2392/hdd/data/neuroimaging/freesurfer_output/";

% raw data EEG dir
RECORD_DATADIR = '/home/adam2392/hdd/data/rawdata';
% RECORD_DATADIR = '/Users/adam2392/Downloads/tngpipeline/';

% subj id to analyze
subjID = 'la03';

% results directory
SUBJDIR = fullfile(RESULTS_DIR, subjID);

FSDIR = '/Users/adam2392/Downloads/';

% filepaths
elecfile = fullfile(strcat(subjID, '_elecint1_f.mat'));
ctimgfile = fullfile(strcat('CT.nii.gz'))
t1imgfile = fullfile(strcat('T1.nii.gz'))

figdir = fullfile('./', subjID);

%% Load Data
% load elecs
elec_data = load(elecfile);
elec_f = elec_data.elecf;

%% Read in Original MRI Scan
t1_img = ft_read_mri(fullfile(SUBJDIR, t1imgfile));

%% Read in Original CT Scan
% read in CT in original format
ct_img = ft_read_mri(fullfile(SUBJDIR, ctimgfile));

fig = figure('units','normalized','outerposition',[0 0 1 1]);
ft_plot_ortho(t1_img.anatomy, 'transform', t1_img.transform, 'style', 'intersect');
ft_plot_sens(elec_f, 'label', 'on', 'fontcolor', 'w');

outputfilename = fullfile(figdir, strcat(subjID, '_sensors'));
% first save 3d screenshot
toSave3DFigFile = strcat(outputfilename, '_3d.png');
saveas(fig, toSave3DFigFile);
   
pause(0.05);
% save a screenshot of the data in multiple views
for i=1:4
    el = 0;
    az = 90*(i-1);
    view([az el])
    toSaveFigFile = strcat(outputfilename, '_', int2str(i), '.fig');
    savefig(fig, toSaveFigFile);
    pause(0.05);
%     print(toSaveFigFile, '-dpng', '-r0')
end
close all

% make the entire figure one image
entirefig = figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:4
   entireax(i) = subplot(2,2,i); 
end

for i=1:4
    currsubplot = subplot(entireax(i));
    toSaveFigFile = strcat(outputfilename, '_', int2str(i), '.fig');
    h = openfig(toSaveFigFile, 'reuse'); % open figure
    ax = gca(h); % get handle of axes of figure
    el = 0;
    az = 90*(i-1);
    view([az el])
    fig = get(ax, 'children');
    copyobj(fig, currsubplot);
    pause(0.05);
    close(h);
end
toSaveFinalFig = strcat(outputfilename, '.png');
saveas(entirefig, toSaveFinalFig);
