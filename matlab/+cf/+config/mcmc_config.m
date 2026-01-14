function mcmc = mcmc_config(preset)
%CF.CONFIG.MCMC_CONFIG  MCMC settings.
% Paper: K=100000, burn=50000, thin=10.
% Default preset here is "debug" (fast) but still posterior-stable enough
% for IRFs once the code is correct.

if nargin < 1 || isempty(preset)
    preset = "debug";  % "debug" | "paper"
end
preset = lower(string(preset));

mcmc = struct();

switch preset
    case "paper"
        mcmc.K    = 100000;
        mcmc.burn = 50000;
        mcmc.thin = 10;

    otherwise % "debug"
        % keep (K-burn)/thin = 1000 draws
        mcmc.K    = 2000;
        mcmc.burn = 500;
        mcmc.thin = 1;
end

% --------------------
% Storage
% --------------------
mcmc.store_full_draws  = true;
mcmc.store_lambda_path = true;   % needed for Appendix-B IRF algorithm (lambda_t path)

% --------------------
% Metropolis / tuning
% --------------------
% P* RW-MH
mcmc.Pstar_init        = 0.05;
mcmc.Pstar_prop_sd0    = 5e-4;
mcmc.Pstar_adapt_start = 200;           % start adapting after burn-in-ish
mcmc.Pstar_target_accept = [0.2 0.4];

% alpha RW-MH: proposal covariance computed in draw_alpha_mh; this scales it.
mcmc.alpha_prop_scale = 0.2;            % typical target acceptance ~ 0.2–0.4

% h_t single-move: multiplier for conditional-prior proposal variance (if used)
mcmc.lambda_prop_scale = 1.0;
% h_t single-move (Jacquier-style accept/reject): cap the number of proposals per time index.
mcmc.lambda_max_tries = 50;

% --------------------
% Utilities
% --------------------
mcmc.make_quick_plots = true;
mcmc.seed = 123;

end
