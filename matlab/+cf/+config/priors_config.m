function pconf = priors_config(pi_bar)
%CF.CONFIG.PRIORS_CONFIG  Prior hyperparameters (Table 2 style).
% pi_bar is the sample mean inflation in YOUR dataset (decimal units).

pconf = struct();

% -------------------------
% P* prior: Truncated Normal around cP * mean(inflation)
% -------------------------
if nargin < 1 || isempty(pi_bar) || ~isfinite(pi_bar)
    % fallback (should not happen if you pass pi_bar from ds)
    pi_bar = 0.035;
end

cP     = 1.5;
muP    = cP * pi_bar;
Pmin   = 0.95 * muP;
Pmax   = 1.05 * muP;
sdP    = (Pmax - Pmin) / 6;

pconf.Pstar = struct();
pconf.Pstar.cP   = cP;
pconf.Pstar.Pbar = pi_bar;
pconf.Pstar.mu   = muP;
pconf.Pstar.sd   = sdP;
pconf.Pstar.Pmin = Pmin;
pconf.Pstar.Pmax = Pmax;

% -------------------------
% d: uniform on {1,...,dmax}
% -------------------------
pconf.d = struct();
pconf.d.dmax = 6;

% -------------------------
% Phi: Minnesota diagonal prior (simple diagonal implementation)
% -------------------------
pconf.Phi = struct();
pconf.Phi.mean      = 0.0;
pconf.Phi.tight     = 0.2;
pconf.Phi.lag_decay = 1.0;

% -------------------------
% alpha: N(0,I)
% -------------------------
na = 22;
pconf.alpha = struct();
pconf.alpha.mu    = zeros(na,1);
pconf.alpha.Omega = eye(na);

% -------------------------
% sigma^2: IG(d/2, d*sigma0^2/2)
% -------------------------
d_sig   = 10;
sig_sd  = 0.01;
sig_var = sig_sd^2;

pconf.sigma2 = struct();
pconf.sigma2.a0 = d_sig/2;
pconf.sigma2.b0 = d_sig*sig_var/2;

% -------------------------
% SV: h_t = log(lambda_t) AR(1)
% mu ~ N(0,1), F ~ N(0.8,0.01) truncated (0,1), Q ~ IG
% IMPORTANT: set Qbar around 1e-3 to match paper scale
% -------------------------
pconf.lambda = struct();

% mu
pconf.lambda.mu0 = 0.0;
pconf.lambda.V0  = 1.0;

% F
pconf.lambda.F0  = 0.8;
pconf.lambda.VF0 = 0.01;  % var (sd=0.1)

% Q (innovation variance of h)
dQ   = 10;
Qbar = 1.0e-3;   % << key: paper posterior around 1.2e-3
pconf.lambda.aQ0 = dQ/2;
pconf.lambda.bQ0 = dQ*Qbar/2;

end