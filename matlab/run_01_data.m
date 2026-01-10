%% run_01_data.m
% Build dataset Zt = (Y, P, U, R, Slope, M2, Pcom, SP500)' from FRED.
% Matches variable definitions in the paper.

clear; clc;

% --- root/matlab detection (robuste) ---
if isfolder(fullfile(pwd,'matlab')) && isfolder(fullfile(pwd,'data'))
    rootDir = pwd;                       % launched from root
elseif isfolder(fullfile(pwd,'+cf'))
    rootDir = fileparts(pwd);            % launched from root/matlab folder
else
    error('Launch from root or root/matlab.');
end

% --- Load configs (package calls) ---
paths = cf.config.paths_config(rootDir);
dconf = cf.config.data_config();
mconf = cf.config.model_config();
pconf = cf.config.priors_config();

fprintf('[run_01_data] Download/build data...\n');

% --- Build dataset ---
ds = cf.data.build_dataset_canova_forero(dconf, paths);
ds_full = ds; % for plots 1960-2023, without restrictions

% Keep sample window
mask = (ds.dates >= dconf.sample_start) & (ds.dates <= dconf.sample_end);
ds.dates = ds.dates(mask);
ds.Z     = ds.Z(mask,:);
ds.names = ds.names;

% Drop initial NaNs due to YoY transforms (12 lags) etc.
good = all(isfinite(ds.Z),2);
ds.dates = ds.dates(good);
ds.Z     = ds.Z(good,:);

% Save (to data/processed by your paths_config)
outFile = fullfile(paths.data_processed, 'dataset_canova_forero.mat');
if ~exist(paths.data_processed,'dir'); mkdir(paths.data_processed); end
save(outFile,'ds','dconf','mconf','pconf','paths','-v7.3');

fprintf('[run_01_data] Saved: %s\n', outFile);

% --- Charts of the data ---

figDir = fullfile(rootDir, 'output', 'figures');
if ~exist(figDir,'dir'); mkdir(figDir); end

% 1) Overview (1960-2023)
cf.plot.dataset_overview_C15(ds_full, figDir, ...
    'XLim', [datetime(1960,1,1) datetime(2023,12,31)], ...
    'FileStem', 'dataset_overview_C15_style');

% 2) Inflation histogram (ici sur 1960-2023 via ds_full)
cf.plot.inflation_hist_normalfit(ds_full, figDir, ...
    'Title', 'PCE Inflation', ...
    'FileStem', 'inflation_hist_normalfit');

% Si tu veux l'histogramme SUR TON SAMPLE D'ESTIMATION, remplace ds_full par ds :
% cf.plot.inflation_hist_normalfit(ds, figDir, 'Title','PCE Inflation (estimation sample)');

