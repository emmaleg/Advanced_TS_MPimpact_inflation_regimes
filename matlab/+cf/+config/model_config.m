function mconf = model_config()
% Model config for eq. (1)-(6): threshold VAR with volatility-in-mean and SV. 

mconf.n = 8;

% --- Not explicitly stated in the BCRP PDF: choose defaults ---
% p: monthly VAR lag length. (Common choice: 12)
mconf.p = 3;

% J: lags of ln(lambda) in mean (includes j=0). Keep small.
mconf.J = 2;

% Threshold regime indicator: S_t = 1 <=> Pi_{t-d} <= P*  (low inflation regime) 
mconf.inflation_index_in_Z = 2; % P is 2nd variable in Z
mconf.dmax = 6;                % paper sets dmax=6 

% Identification / restrictions
mconf.id.n_alpha = 22; % implied by A.8 
mconf.id.enforce_sign_zero = true; % you can switch on later
mconf.id.zero_tol = 1e-4;          % final : 1e-10 // for debug 1e-4 ou 1e-6 // numerical tolerance for "zero" checks

% Stationarity truncation for Phi draws (Appendix A mentions truncation I(phi)) 
mconf.enforce_stationarity = true;
mconf.stationarity_max_tries = 200;

end
