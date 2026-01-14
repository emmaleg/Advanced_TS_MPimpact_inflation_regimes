function st = init_state(Z, dates, mconf, pconf, mcmc)
% Initialize parameters.

[T,n] = size(Z);
st = struct();

st.Pstar = mcmc.Pstar_init;
st.d     = 1;

% SV: h_t = ln lambda_t
st.h  = zeros(T,1); % lambda=1 initially
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

% initial regimes (STORE Svalid!)
Pi = Z(:, mconf.inflation_index_in_Z);
[st.S, st.Svalid] = cf.model.regime_indicator(Pi, st.Pstar, st.d);

% If we enforce identification during estimation, we MUST start inside support
if isfield(mconf,'id') && isfield(mconf.id,'enforce_sign_zero') && mconf.id.enforce_sign_zero
    tol     = mconf.id.zero_tol;
    lambda0 = exp(mean(st.h)); % >0, scaling doesn't affect signs

    [st.alpha1, st.A1] = draw_feasible_alpha(pconf, st.sigma2, lambda0, tol, mconf.id.n_alpha, 5000);
    [st.alpha2, st.A2] = draw_feasible_alpha(pconf, st.sigma2, lambda0, tol, mconf.id.n_alpha, 5000);
end

% MH bookkeeping
st.accept.Pstar = 0;
st.accept.Pstar_trials = 0;
st.accept.alpha1 = 0; st.accept.alpha1_trials=0;
st.accept.alpha2 = 0; st.accept.alpha2_trials=0;
st.accept.h = 0; st.accept.h_trials = 0;
st.accept.h_updates = 0;

% For adaptive P* proposal
st.Pstar_prop_sd = mcmc.Pstar_prop_sd0;
st.Pstar_hist = NaN(mcmc.K,1);

end

% -------- local helper --------
function [alpha, A] = draw_feasible_alpha(pconf, sigma2, lambda0, tol, na, maxTries)
mu = pconf.alpha.mu(:);
Om = (pconf.alpha.Omega + pconf.alpha.Omega')/2;

% robust chol
jitter = 1e-12;
for it=1:8
    [L,p] = chol(Om + jitter*eye(size(Om)), 'lower');
    if p==0, break; end
    jitter = jitter*10;
end
if p~=0
    L = chol(Om + 1e-6*eye(size(Om)), 'lower');
end

for k=1:maxTries
    alpha = mu + L*randn(na,1);
    A     = cf.id.alpha_to_A(alpha);
    if cf.id.check_impact_restrictions(A, sigma2, lambda0, tol)
        return;
    end
end

error("init_state: cannot find feasible alpha after %d tries (tol=%g).", maxTries, tol);
end

