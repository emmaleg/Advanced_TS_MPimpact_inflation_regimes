function iconf = irf_config(mconf)
%CF.CONFIG.IRF_CONFIG  Settings for IRF computation (Appendix B).
% Paper: H=36, S=1000, L=200, delta=1.

iconf = struct();

% Sample with replacement from posterior draws (recommended)
iconf.remise = true;

% Horizons / simulation counts (paper)
iconf.H = 36;      % months
iconf.S = 1000;    % number of posterior draws used for IRFs
iconf.L = 200;     % Monte Carlo reps per draw

% Shock magnitude parameter (Appendix B): delta=1
iconf.delta = 1;

% -------------------------------------------------------------------------
% Shock normalization
% Paper algorithm sets structural innovation e_1^δ = δ and e_1^0 = 0.
% Start with "structural" to match Appendix B exactly.
%
% Options in your codebase (if implemented):
%   "structural"  : set e_j = +/- delta in the structural shock itself
%   "target_abs"  : rescale so impact(target_var) = +/- delta
%   "target_1sd"  : rescale so impact(target_var) = +/- delta * sd(innovation(target))
% -------------------------------------------------------------------------
iconf.shock_norm_mode = "structural";
iconf.min_norm_impact = 1e-6;     % safeguard if you switch to target_* modes
iconf.target_var_conventional = 4;
iconf.target_var_liquidity    = 6;

% Liquidity shock: keep short-term rate fixed for Hlock months (paper: 24)
iconf.ffr_lock_h = 24;  % Appendix H sometimes uses 12; main benchmark is 24

% Variable indices (paper ordering)
iconf.idx_ip    = 1;
iconf.idx_inf   = 2;
iconf.idx_u     = 3;
iconf.idx_ffr   = 4;
iconf.idx_slope = 5;
iconf.idx_m2    = 6;
iconf.idx_pcom  = 7;
iconf.idx_sp500 = 8;

% Targets (only used if you select target_* modes)
iconf.target_var_conventional = iconf.idx_ffr;
iconf.target_var_liquidity    = iconf.idx_m2;

% Shock indices (as in your identification: columns 4 and 6)
iconf.shock_conventional = 4;
iconf.shock_liquidity    = 6;

% Control shock used to enforce FFR lock (the code uses this to offset FFR)
iconf.shock_ffr_control  = 4;

% Model lags
iconf.p = mconf.p;
iconf.J = mconf.J;

% Needed for safe sampling of t* (since regime uses inflation_{t-d})
iconf.dmax = mconf.dmax;

% RNG
iconf.seed = 1;

end



