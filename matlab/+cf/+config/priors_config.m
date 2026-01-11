function pconf = priors_config()
% Prior hyperparameters, aligned with Table 2

pconf = struct();

%% --- P* : Truncated Normal around cP * Pbar ---
% Pbar is the sample mean of inflation. For now, fixed to 0.0341 (given in
% the article)
pi_bar = 0.0341;
cP     = 1.5;

muP  = cP * pi_bar;
Pmin = 0.95 * muP;
Pmax = 1.05 * muP;
sdP  = (Pmax - Pmin) / 6;

pconf.Pstar.cP   = cP;
pconf.Pstar.Pbar = pi_bar;

pconf.Pstar.mu   = muP;
pconf.Pstar.sd   = sdP;
pconf.Pstar.Pmin = Pmin;
pconf.Pstar.Pmax = Pmax;

% NOTE: if your draw_Pstar_mh assumes UNIFORM prior, you must add the log prior
% term there to be strictly paper-consistent.

%% --- d : uniform on {1,...,dmax} ---
pconf.d.dmax = 6;  % paper: dmax=6

%% --- Phi (VAR coefficients): Minnesota diagonal, mean 0 ---
% Paper: phi_i = 0, V_i diagonal Minnesota.
pconf.Phi.mean      = 0.0;
pconf.Phi.tight     = 0.2;   % reasonable Minnesota overall tightness (paper says "Minnesota", not the exact number)
pconf.Phi.lag_decay = 1.0;

%% --- alpha (A.8 free parameters): N(0, I) ---
na = 22; % paper: dim(alpha)=22 from A.8
pconf.alpha.mu    = zeros(na,1);
pconf.alpha.Omega = eye(na);

%% --- sigma^2_j : IG with d_sigma=10, sigma=0.01 ---
% Use IG(a,b) with density proportional x^(-a-1) exp(-b/x).
% Match scaled-inv-chi-square: a = d/2, b = d*s^2/2.
d_sig = 10;
s2_sig = 0.01;

pconf.sigma2.a0 = d_sig/2;
pconf.sigma2.b0 = d_sig*s2_sig/2;

%% --- h_t = log(lambda_t) AR(1): mu, F, Q ---
% mu ~ N(0,1)
pconf.lambda.mu0 = 0.0;
pconf.lambda.V0  = 1.0;

% F ~ TruncN(0.8, 0.01) on (0,1)  -> sd = 0.1
pconf.lambda.F0  = 0.8;
pconf.lambda.VF0 = 0.01;

% Q ~ IG(dQ/2, dQ*Qbar/2) with dQ=10, Qbar=0.01 (common choice consistent with d_sigma, sigma)
dQ   = 10;
Qbar = 0.01;

pconf.lambda.aQ0 = dQ/2;
pconf.lambda.bQ0 = dQ*Qbar/2;

end