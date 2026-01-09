function paths = paths_config(rootDir)
%CF.CONFIG.PATHS_CONFIG  Centralize project folders based on repo structure.
%
% Expected repo structure:
%   root/
%     data/raw/
%     data/processed/
%     output/mcmc/
%     output/figures/
%     output/irf/
%     matlab/+cf/...

paths = struct();
paths.root = rootDir;

% Data
paths.data_root      = fullfile(rootDir, "data");
paths.data_raw       = fullfile(rootDir, "data", "raw");
paths.data_processed = fullfile(rootDir, "data", "processed");

% (Optional) FRED-MD dedicated subfolders (match your tree)
% paths.fredmd_root    = fullfile(paths.data_raw, "fredmd");
% paths.fredmd_vintage = fullfile(paths.fredmd_root, "vintage");
% paths.fredmd_meta    = fullfile(paths.fredmd_root, "meta");

% Output
paths.output_root    = fullfile(rootDir, "output");
paths.output_mcmc    = fullfile(paths.output_root, "mcmc");
paths.output_figures = fullfile(paths.output_root, "figures");
paths.output_irf     = fullfile(paths.output_root, "irf");

% Create folders if missing
mk = {paths.data_raw, paths.data_processed, ...
      paths.output_root, paths.output_mcmc, paths.output_figures, paths.output_irf};

% paths.fredmd_root, paths.fredmd_vintage, paths.fredmd_meta, ...

for i = 1:numel(mk)
    if ~exist(mk{i}, "dir")
        mkdir(mk{i});
    end
end

end
