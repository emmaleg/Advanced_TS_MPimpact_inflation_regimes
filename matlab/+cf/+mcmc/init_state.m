function st = init_state(Z, dates, mconf, pconf, mcmc)
% Initialize parameters.

[T,n] = size(Z);
st = struct();

st.Pstar = mcmc.Pstar_init;
st.d = 1;

% SV: h_t = ln lambda_t
st.h = zeros(T,1); % lambda=1 initially
st.mu = pconf.lambda.mu0;
st.F  = min(max(pconf.lambda.F0, 0.01), 0.99);
st.Q  = pconf.lambda.bQ0;

% sigma^2 diagonal
st.sigma2 = 0.01*ones(n,1);

% Phi1, Phi2
kPhi = 1 + n*mconf.p + (mconf.J+1);
st.Phi1 = zeros(n, kPhi);
st.Phi2 = zeros(n, kPhi);

% alpha -> A
st.alpha1 = zeros(mconf.id.n_alpha,1);
st.alpha2 = zeros(mconf.id.n_alpha,1);
st.A1 = cf.id.alpha_to_A(st.alpha1);
st.A2 = cf.id.alpha_to_A(st.alpha2);

% initial regimes
Pi = Z(:, mconf.inflation_index_in_Z);
st.S = cf.model.regime_indicator(Pi, st.Pstar, st.d);

% MH bookkeeping
st.accept.Pstar = 0;
st.accept.Pstar_trials = 0;
st.accept.alpha1 = 0; st.accept.alpha1_trials=0;
st.accept.alpha2 = 0; st.accept.alpha2_trials=0;
st.accept.h = 0; st.accept.h_trials = 0;

% For adaptive P* proposal
st.Pstar_prop_sd = mcmc.Pstar_prop_sd0;
st.Pstar_hist = NaN(mcmc.K,1);

end
