function iconf = irf_config(mconf)
%CF.CONFIG.IRF_CONFIG  Settings for IRF computation (Appendix B).

iconf = struct();

iconf.remise = false; % for the draws 

% Horizons / simulation counts (paper: H=36, S=1000, L=200)
iconf.H     = 36;
iconf.S     = 1000;   % number of IRF draws (outer loop)
iconf.L     = 200;    % Monte Carlo reps per draw to approximate expectations

% Shock magnitude parameter:
% - if shock_norm_mode == "structural"   : e_j = +/- delta
% - if shock_norm_mode == "target_abs"   : impact(target) = +/- delta
% - if shock_norm_mode == "target_1sd"   : impact(target) = +/- delta * sd(innovation(target))
% Paper conventional shock is "one standard deviation surprise increase in the short-term rate".
iconf.delta = 1;

% --- shock normalization ---
iconf.shock_norm_mode = "structural";      % "structural" | "target_abs" | "target_1sd"
iconf.min_norm_impact = 1e-6;              % safeguard when B(target,shock) is tiny
iconf.target_var_conventional = 4;         % set below after indices are defined
iconf.target_var_liquidity    = 6;         % default: normalize liquidity by M2 impact

% Liquidity shock: keep short-term rate fixed for Hlock months (paper: 24; appendix H uses 12)
iconf.ffr_lock_h = 24;

% Variable indices (paper order)
iconf.idx_ip    = 1;
iconf.idx_inf   = 2;
iconf.idx_u     = 3;
iconf.idx_ffr   = 4;
iconf.idx_slope = 5;
iconf.idx_m2    = 6;
iconf.idx_pcom  = 7;
iconf.idx_sp500 = 8;

% Now that indices exist, set targets cleanly:
iconf.target_var_conventional = iconf.idx_ffr;  % paper normalization
iconf.target_var_liquidity    = iconf.idx_m2;   % reasonable default; change if desired

% Shock indices (Appendix A matrix A, columns 4 and 6 are the policy shocks)
iconf.shock_conventional = 4;
iconf.shock_liquidity    = 6;

% Control shock used to enforce FFR lock (use conventional shock #4)
iconf.shock_ffr_control  = 4;

% Model lags
iconf.p = mconf.p;      % paper: p=3
iconf.J = mconf.J;      % paper: J=2

% RNG
iconf.seed = 1;

% For safe sampling of t* (need p lags, J lags of lambda, and d up to 6)
iconf.dmax = mconf.dmax;

end



