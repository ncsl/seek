% REFERENCE: http://www.fieldtriptoolbox.org/tutorial/human_ecog

%% Setup Paths and Directories for Registration Results and Raw EEG data_examples
addpath("/home/adam2392/Documents/Dropbox/fieldtrip-20181108");
addpath("/Users/adam2392/Dropbox/fieldtrip-20181108");
addpath("/Users/adam2392/Documents/MATLAB/spm12");

% add MNE-Matlab
addpath("$MNE_ROOT/share/matlab");
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
% FSDIR = '/Users/adam2392/Downloads/';

eegfilename = strcat(subjID, '_sz_1p_fragmodel.mat');
eegfile = fullfile('./', eegfilename);
atlasfile = fullfile(SUBJDIR, 'aparc+aseg.mgz')

testdata = load('./SubjectUCI29_data.mat');
testdata = testdata.data;

%% Load Brain data_examples
% load in atlas-labeled brain volume
atlas = ft_read_atlas(atlasfile); 
% atlas.coordsys = 'acpc';
atlas.coordsys = 'mni';
cfg            = [];
% cfg.inputcoord = 'acpc';
cfg.inputcoord = 'mni';
cfg.atlas      = atlas;
cfg.roi        = {
%                 'ctx-rh-superiorparietal',
                    'ctx-rh-parahippocampal',
                    'ctx-rh-paracentral',
                    'ctx-rh-frontalpole',
%                 'Right-Hippocampus', 
%                 'Right-Amygdala',
%                 'Left-Hippocampus', 
%                 'Left-Amygdala',
                };
% cfg.roi = atlas.aparclabel;
mask_rha = ft_volumelookup(cfg, atlas);

% create mesh
seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'}); 
seg.brain = mask_rha;
cfg             = [];
cfg.method      = 'iso2mesh';
cfg.radbound    = 2;
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 1000;
cfg.smooth      = 3;
cfg.spmversion = 'spm12';
mesh_rha = ft_prepare_mesh(cfg, seg);

%% Example Freq data_examples on testdata
% compute FFT
% cfg            = [];
% cfg.method     = 'mtmconvol';
% cfg.toi        = -.3:0.01:.8;
% cfg.foi        = 5:5:200;
% cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.2;
% cfg.taper      = 'hanning'; 
% cfg.output     = 'pow';
% cfg.keeptrials = 'no';
% freq = ft_freqanalysis(cfg, testdata);
% 
% % relative change to baseline
% cfg              = [];
% cfg.baseline     = [-.3 -.1]; 
% cfg.baselinetype = 'relchange'; 
% freq_blc = ft_freqbaseline(cfg, freq);
% 
% % average over time
% cfg             = []; 
% cfg.frequency   = [70 150]; 
% cfg.avgoverfreq = 'yes'; 
% cfg.latency     = [0 0.8]; 
% cfg.avgovertime = 'yes';
% freq_sel = ft_selectdata(cfg, freq_blc);

%% Example Fragility data_examples
data = load(eegfile);

% filepaths
elecfile = fullfile(strcat(subjID, '_elecint1_f.mat'));
% load elecs
elec_data = load(elecfile);
elec_f = elec_data.elecf;

data.elec = elec_f
data.label = cellstr(data.chanlabels);

chaninds = [];
% loop through chanlabels and get indices of xyz labels to keep
for i=1:length(data.elec.label),
    eleclabel = lower(data.elec.label{i})
    if ~any(strcmp(data.label, eleclabel))
        chaninds = [chaninds; i];
    end
    data.elec.label{i} = eleclabel;
end

for i=1:length(data.elec.cfg.channel)
    eleclabel = lower(data.elec.cfg.channel{i});
    data.elec.cfg.channel{i} = eleclabel;
end

chaninds
data.elec.elecpos(chaninds,:) = [];
data.elec.chanpos(chaninds,:) = [];
data.elec.label(chaninds) = []
data.elec.tra(chaninds,:) = [];
data.elec.tra(:,chaninds) = [];

for i=1:length(data.label)
    data.elec.chanunit{i} = 'V';
    data.elec.chantype{i} = 'eeg';
end
data.elec.chantype = transpose(data.elec.chantype)
data.elec.chanunit = transpose(data.elec.chanunit)
data.elec.coordsys = 'mni';

% construct result data structure
result = struct();
result.label = data.label;
result.elec = data.elec
result.dimord = 'chan_freq_time';
result.freq = 1;
result.powspctrm = data.fragmat(:,:);
result.time = linspace(1,size(data.fragmat,2),size(data.fragmat,2));

% % average over time
cfg             = []; 
cfg.avgovertime = 'yes';
result = ft_selectdata(cfg, result);
frag_sel = result;
% get subcortical electrodes of interest
% cfg         = [];
% cfg.channel = {'l*', 'w*'}; 
% frag_sel = ft_selectdata(cfg, result);
% frag_sel

% freq_sel.powspctrm = result.powspctrm;
% freq_sel.label = result.label;
% freq_sel.elec. = result.elec

%% Plotting Setup of SEEG
cfg              = []; 
cfg.funparameter = 'powspctrm';
cfg.funcolorlim  = [0 1];
cfg.method       = 'cloud';
cfg.funcolormap = 'jet';
cfg.slice        = '3d';
cfg.nslices      = 6;
cfg.facealpha    = .25;
% fig = figure;
ft_sourceplot(cfg, frag_sel, mesh_rha);
view([120 40]); 
lighting gouraud; 
camlight;

%% Plot SEEG
% fig = figure;
cfg.slice        = '2d';
ft_sourceplot(cfg, frag_sel, mesh_rha);


%% PLOT ALL IN A VIDEO
% create the video writer with 1 fps
videoname = './myVideo.avi'
writerObj = VideoWriter(videoname);
writerObj.FrameRate = 10;
% open the video writer
open(writerObj);

for i=1:size(data.fragmat,2)
    % construct result data structure
    result = struct();
    result.label = data.label;
    result.elec = data.elec
    result.dimord = 'chan_freq_time';
    result.freq = 1;
    result.powspctrm = data.fragmat(:,i);
    result.time = 1
%     linspace(1,size(data.fragmat,2),size(data.fragmat,2));

    % % average over time
    % cfg             = []; 
    % cfg.avgovertime = 'yes';
    % % cfg.frequency   = 1; 
    % % cfg.avgoverfreq = 'yes'; 
    % result = ft_selectdata(cfg, result);

    % get subcortical electrodes of interest
    cfg         = [];
    cfg.channel = {'l*', 'w*'}; 
    frag_sel = ft_selectdata(cfg, result);
    % frag_sel

    % freq_sel.powspctrm = result.powspctrm;
    % freq_sel.label = result.label;
    % freq_sel.elec. = result.elec

    %% Plotting Setup of SEEG
    cfg              = []; 
    cfg.funparameter = 'powspctrm';
    cfg.funcolorlim  = [0 1];
    cfg.method       = 'cloud';
    cfg.funcolormap = 'jet';
    cfg.slice        = '3d';
    cfg.nslices      = 2;
    cfg.facealpha    = .25;
    
    ft_sourceplot(cfg, frag_sel, mesh_rha);
    fig = gcf();
    ax = gca();
    view([120 40]); 
    lighting gouraud; 
    camlight;

    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    close all
    %% Plot SEEG
%     cfg.slice        = '2d';
%     ft_sourceplot(cfg, frag_sel, mesh_rha);
end
% close the writer object
close(writerObj);