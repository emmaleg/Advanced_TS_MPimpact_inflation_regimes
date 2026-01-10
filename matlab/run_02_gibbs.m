%% run_02_gibbs.m
% Run Gibbs sampler for threshold-BVAR with SV (single scalar lambda_t).
% Steps follow section 4 and Appendix A.

clear; clc;

% --- root/matlab detection (robuste) ---
if isfolder(fullfile(pwd,'matlab')) && isfolder(fullfile(pwd,'data'))
    rootDir = pwd;                       % if launched from the root
elseif isfolder(fullfile(pwd,'+cf'))
    rootDir = fileparts(pwd);            % if launched from root/matlab
else
    error('Launch from root or from root/matlab');
end

addpath(genpath(fullfile(rootDir,'matlab')));

% --- Configs ---
paths = cf.config.paths_config(rootDir);
mcmc  = cf.config.mcmc_config();
mconf = cf.config.model_config();
pconf = cf.config.priors_config();

% --- Load dataset ---
inFile = fullfile(paths.data_processed, 'dataset_canova_forero.mat');
load(inFile,'ds','dconf');

fprintf('[run_02_gibbs] Dataset loaded: T=%d, n=%d\n', size(ds.Z,1), size(ds.Z,2));

% (petites sécurités)
assert(size(ds.Z,2) == mconf.n, 'Incohérence: size(ds.Z,2) ~= mconf.n');
assert(mconf.inflation_index_in_Z == 2, 'Inflation has to be in column 2.');

% --- Gibbs ---
res = cf.mcmc.gibbs_threshold_bvar_sv(ds, mconf, pconf, mcmc);

% --- Save ---
if ~exist(paths.output_mcmc,'dir'); mkdir(paths.output_mcmc); end
outFile = fullfile(paths.output_mcmc, 'posterior_draws.mat');
save(outFile,'res','mconf','pconf','mcmc','dconf','-v7.3');
fprintf('[run_02_gibbs] Saved: %s\n', outFile);

% --- Quick plots ---
if isfield(mcmc,'make_quick_plots') && mcmc.make_quick_plots
    figDir = fullfile(paths.output_figures, 'quickplots');
    if ~exist(figDir,'dir'); mkdir(figDir); end

    cf.plot.plot_inflation_high_regime(ds, mconf, res, ...
        'savePath', fullfile(figDir,'inflation_high_regime.png'));

    cf.plot.plot_posterior_Pstar(res, 'savePath', fullfile(figDir,'posterior_Pstar.png'));
    cf.plot.plot_posterior_dstar(res, 'savePath', fullfile(figDir,'posterior_dstar.png'));
    cf.plot.plot_posterior_F(res,     'savePath', fullfile(figDir,'posterior_F.png'));
    cf.plot.plot_posterior_mu(res,    'savePath', fullfile(figDir,'posterior_mu.png'));
    cf.plot.plot_posterior_Q(res,     'savePath', fullfile(figDir,'posterior_Q.png'));
end


% % --- Quick plots ---
% if isfield(mcmc,'make_quick_plots') && mcmc.make_quick_plots
%     Sm = res.S_mean;
%     high = (Sm <= 0.5);
% 
%     figure('Name','Fig3-like: Inflation and High-regime indicator');
%     yyaxis left;
%     plot(ds.dates, ds.Z(:, mconf.inflation_index_in_Z), 'LineWidth', 1.2);
%     ylabel('Inflation');
% 
%     yyaxis right;
%     stem(ds.dates(high), ones(nnz(high),1), 'filled');
%     ylim([0 1.2]);
%     ylabel('High regime indicator');
% 
%     figure('Name','Fig4-like: Posterior of P*');
%     histogram(res.draws.Pstar, 50);
%     xlabel('P*'); ylabel('count');
% end
