function mconf = model_config()
%CF.CONFIG.MODEL_CONFIG  Model config for the Threshold-BVAR-SV.
% Paper (v7): p = 3, J = 2, dmax = 6, dim(Z)=8.

mconf = struct();

% Dimensionality
mconf.n = 8;

% VAR lags (paper: p=3)
mconf.p = 3;

% Lags of log(lambda_t) in the mean (paper: J=2, includes j=0,...,J)
mconf.J = 2;

% Threshold regime indicator:
%   S_t = 1  <=>  Pi_{t-d} <= P*   (low inflation regime)
mconf.inflation_index_in_Z = 2;  % inflation is the 2nd variable in Z
mconf.dmax = 6;                  % paper: dmax = 6

% ------------------------------------------------------------
% Identification / impact restrictions (Table 1)
% ------------------------------------------------------------
mconf.id = struct();
mconf.id.n_alpha = 22;                 % from Appendix A.8
mconf.id.enforce_sign_zero = true;     % enforce Table-1 constraints DURING sampling
mconf.id.zero_tol = 1e-3;              % relative tolerance for "zero" (robust)

% ------------------------------------------------------------
% Stationarity truncation for Phi draws (Appendix A)
% ------------------------------------------------------------
mconf.enforce_stationarity = true;
mconf.stationarity_max_tries = 200;

end


