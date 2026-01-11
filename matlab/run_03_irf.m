%% run_03_irf.m
% Compute regime-dependent IRFs (Appendix B) and reproduce Figures (regime comparison).

clear; clc;

% --- root/matlab detection (robuste) ---
if isfolder(fullfile(pwd,'matlab')) && isfolder(fullfile(pwd,'data'))
    rootDir = pwd;                       % launched from root
elseif isfolder(fullfile(pwd,'+cf'))
    rootDir = fileparts(pwd);            % launched from root/matlab
else
    error('Launch from root or root/matlab.');
end

% --- Configs ---
paths = cf.config.paths_config(rootDir);
mconf = cf.config.model_config();

% Load dataset
dataFile = fullfile(paths.data_processed, "dataset_canova_forero.mat");
S = load(dataFile);  % expects ds, mconf, paths in there
ds = S.ds;

% Load posterior draws from Gibbs
postFile = fullfile(paths.output_mcmc, "posterior_draws.mat");
P = load(postFile);
post = P.res.draws;

% IRF config
iconf = cf.config.irf_config(mconf);

fprintf("[run_03_irf] Computing IRFs (this can be long: S=%d, L=%d)...\n", iconf.S, iconf.L);

out = cf.irf.compute_irfs(ds, post, mconf, iconf);

% Save
if ~exist(paths.output_irf,'dir'); mkdir(paths.output_irf); end
save(fullfile(paths.output_irf, "irf_results.mat"), "out", "iconf", "mconf", "paths", "-v7.3");

% Plot: conventional vs liquidity (regime comparison)
cf.plot.make_fig_irf_regimes(out.conventional, ds, paths, "fig_irf_conventional_regimes", ...
    "Contractionary conventional monetary policy shock (Low vs High)");
cf.plot.make_fig_irf_regimes(out.liquidity, ds, paths, "fig_irf_liquidity_regimes", ...
    "Expansionary liquidity shock (Low vs High) - FFR locked");

fprintf("[run_03_irf] Done. Saved IRFs to output/irf and figures to output/figures.\n");
