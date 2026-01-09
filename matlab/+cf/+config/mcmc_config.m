function mcmc = mcmc_config()
% MCMC settings: paper uses K=100,000; burn 50,000; keep 1 every 10. 
% For debugging, reduce K drastically.

mcmc.K    = 1000; % for tests now 
mcmc.burn = 500;
mcmc.thin = 10;

mcmc.store_full_draws = true; % store Phi, alpha, sigma, (optionally lambda)
mcmc.store_lambda_path = false; % huge memory if true

% MH tuning (P* and alpha, lambda)
mcmc.Pstar_init = 0.05;
mcmc.Pstar_prop_sd0 = 0.0005;    % initial RW sd
mcmc.Pstar_adapt_start = 2000;   % start adapting after some iterations
mcmc.Pstar_target_accept = [0.2 0.4]; % per appendix 
mcmc.alpha_prop_scale = 1.0;     % scales proposal covariance in alpha MH
mcmc.lambda_prop_sd = 0.15;      % for fallback RW; main proposal uses conditional prior

% Quick plots at end of run_02_gibbs
mcmc.make_quick_plots = true;

% Random seed
mcmc.seed = 123;

end

