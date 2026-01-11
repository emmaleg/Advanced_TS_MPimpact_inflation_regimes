function mcmc = mcmc_config()
% MCMC settings (paper: K=100,000; burn 50,000; thin 10).
% Here: "debug-stable" preset: fast but enough kept draws for stable IRFs.

% --- Main iterations ---
% keep (K-burn)/thin = 700 draws
mcmc.K    = 5000;   % 100000
mcmc.burn = 1500;   % 50000
mcmc.thin = 5;      % 10

% --- Storage ---
mcmc.store_full_draws   = true;
mcmc.store_lambda_path  = true; % IMPORTANT: keep false in debug (huge memory)

% --- MH tuning ---
mcmc.Pstar_init = 0.05;
mcmc.Pstar_prop_sd0 = 0.0005;
mcmc.Pstar_adapt_start = 50;           % 200 % avoid early instability)
mcmc.Pstar_target_accept = [0.2 0.4];

% alpha MH: proposal covariance is V_alpha; scale controls acceptance rate
mcmc.alpha_prop_scale = 0.8;            % target ~0.2–0.4

% lambda/h MH (used in the improved single-move below)
mcmc.lambda_prop_scale = 1.0;           % multiplies the Laplace variance

% before : 
% mcmc.alpha_prop_scale = 1.0;     % scales proposal covariance in alpha MH
% mcmc.lambda_prop_sd = 0.15;      % for fallback RW; main proposal uses conditional prior

% Quick plots
mcmc.make_quick_plots = true;

% Seed
mcmc.seed = 123;

end
