%% run_01_data.m
% Build dataset Zt = (Y, P, U, R, Slope, M2, Pcom, SP500)' from FRED.
% Matches variable definitions in the paper.

clear; clc;

rootDir = fileparts(pwd); % if pwd = repo/matlab
% you have to position yourself in the root folder of the repo ! 
% TODO implement something more robust

% --- Load configs (package calls) ---
paths = cf.config.paths_config(rootDir);
dconf = cf.config.data_config();
mconf = cf.config.model_config();
pconf = cf.config.priors_config();

fprintf('[run_01_data] Download/build data...\n');

% --- Build dataset ---
ds = cf.data.build_dataset_canova_forero(dconf, paths);

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