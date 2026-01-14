function st = draw_lambda_single_move(st, Z, mconf, pconf, mcmc)
%CF.MCMC.DRAW_LAMBDA_SINGLE_MOVE  Jacquier et al. (1994) single-move draw of h_t = ln(lambda_t).
%
% Implements Canova & Perez-Forero Online Appendix A, Block 6, eqs. (A.15)-(A.24). :contentReference[oaicite:5]{index=5}
%
% KEY DETAIL for the multivariate system:
%   Appendix formulas are written in univariate form with y_t^2.
%   In the VAR with common scalar volatility, use:
%       y_t^2 = (1/n) * u_t' * Sigma^{-1} * u_t
%   where u_t = A_i * eps_t are the structural innovations and Sigma = diag(sigma2).
%   This keeps the scale consistent with (A.20)-(A.24). :contentReference[oaicite:6]{index=6}
%
% Note: We follow the Appendix "single-move" accept/reject step with convexity bound.

if nargin < 5 || isempty(mcmc)
    mcmc = struct();
end
if ~isfield(mcmc, 'lambda_max_tries') || isempty(mcmc.lambda_max_tries)
    mcmc.lambda_max_tries = 50;  % safety cap (rarely hit)
end

[T, n] = size(Z);

% --- initialize acceptance counters if missing ---
if ~isfield(st,'accept') || isempty(st.accept)
    st.accept = struct();
end
if ~isfield(st.accept,'h');         st.accept.h = 0; end
if ~isfield(st.accept,'h_trials');  st.accept.h_trials = 0; end
if ~isfield(st.accept,'h_updates'); st.accept.h_updates = 0; end

% -------------------------------------------------------------------------
% Precompute y2_t = (1/n) * u_t' Sigma^{-1} u_t, holding u_t FIXED during the sweep.
% This matches the Appendix derivation where y_t is treated as "data" in the bound. :contentReference[oaicite:7]{index=7}
% -------------------------------------------------------------------------
y2 = compute_y2_path(st, Z, mconf); % T x 1, NaN when undefined

mu = st.mu;
F  = st.F;
Q  = st.Q;

% --- Single-move accept/reject per time index t ---
for t = 1:T

    % Conditional prior (A.16)-(A.18) for interior points; edges handled with stationary prior
    [h_star, v2] = cond_prior_h(t, st.h, mu, F, Q);

    % If y2(t) missing (early sample / missing regressors), draw from conditional prior
    if ~isfinite(y2(t))
        st.h(t) = h_star + sqrt(max(v2, 1e-12))*randn();
        continue
    end

    % (A.23): mu_t = h*_t + (v2/2) * ( y_t^2 * exp(-h*_t) - 1 ) :contentReference[oaicite:8]{index=8}
    mu_t = h_star + 0.5*v2*( y2(t)*exp(-h_star) - 1 );

    % Draw candidate and accept with prob alpha = min(1, f*_t / g*_t) (A.24)
    accepted = false;
    tries = 0;

    while ~accepted
        tries = tries + 1;
        st.accept.h_trials = st.accept.h_trials + 1;

        % Candidate h_t^c ~ N(mu_t, v2)  (A.22)-(A.23)
        h_can = mu_t + sqrt(max(v2, 1e-12))*randn();

        % log f*(.) from (A.20) left-hand side:  -1/2 h - y^2/2 exp(-h)
        ln_f = -0.5*h_can - 0.5*y2(t)*exp(-h_can);

        % log g*(.) from (A.20) right-hand side using convexity bound:
        % exp(-h) <= exp(-h*) (1 + h* - h)
        ln_g = -0.5*h_can - 0.5*y2(t)*exp(-h_star)*(1 + h_star - h_can);

        log_alpha = ln_f - ln_g; % should be <= 0

        if log(rand()) < min(0, log_alpha)
            st.h(t) = h_can;
            st.accept.h = st.accept.h + 1;
            accepted = true;
        else
            if tries >= mcmc.lambda_max_tries
                % Safety fallback: keep current value (should almost never happen)
                accepted = true;
            end
        end
    end

end

st.accept.h_updates = st.accept.h_updates + 1;

end

% =========================================================================
% Helpers
% =========================================================================

function y2 = compute_y2_path(st, Z, mconf)
% Compute y2_t = (1/n) * u_t' Sigma^{-1} u_t, where u_t = A_i eps_t.
% eps_t = Z_t - Phi_i x_t. x_t includes lags of Z and lags of h (since ln lambda enters mean). :contentReference[oaicite:9]{index=9}

[T, n] = size(Z);
y2 = NaN(T,1);

% Build regressors x_t = [1, lagged Z, h_t, h_{t-1}, ..., h_{t-J}]
Xlag = cf.model.make_lag_matrix(Z, mconf.p);

Hlags = NaN(T, mconf.J+1);
for j = 0:mconf.J
    Hlags(:, j+1) = lagmatrix(st.h(:), j);
end
Xfull = [ones(T,1), Xlag, Hlags];

% For skipping invalid times
Svalid = true(T,1);
if isfield(st,'Svalid') && ~isempty(st.Svalid)
    Svalid = st.Svalid;
end

for t = 1:T
    if ~Svalid(t), continue; end
    if ~all(isfinite(Xfull(t,:))) || ~all(isfinite(Z(t,:)))
        continue
    end

    if st.S(t)==1
        Phi = st.Phi1;
        A   = st.A1;
    else
        Phi = st.Phi2;
        A   = st.A2;
    end

    eps_t = Z(t,:)' - Phi*Xfull(t,:)';
    u_t   = A*eps_t;

    ss_t = sum((u_t.^2) ./ st.sigma2(:));  % u' Sigma^{-1} u
    y2(t) = ss_t / n;                      % normalization to match Appendix scalar form
end

end


function [m, v] = cond_prior_h(t, h, mu, F, Q)
% Conditional prior moments for AR(1):
% h_t - mu = F (h_{t-1} - mu) + eta_t, eta_t ~ N(0,Q)   :contentReference[oaicite:10]{index=10}

h  = h(:);
T  = numel(h);
F  = min(max(F, 0.001), 0.9999);
Q  = max(Q, 1e-12);

V0 = Q / max(1e-12, (1 - F^2)); % stationary var for h1

if t == 1
    if T == 1
        m = mu; v = V0; return
    end
    % Combine: h1 ~ N(mu,V0) and h2|h1 ~ N(mu + F(h1-mu), Q)
    v = 1 / (1/V0 + (F^2)/Q);
    m = v * (mu/V0 + (F*(h(2) - (1-F)*mu))/Q);
    return
end

if t == T
    m = mu + F*(h(T-1) - mu);
    v = Q;
    return
end

% Interior: (A.17)-(A.18) form
v = Q / (1 + F^2);
m = mu + (F*((h(t-1)-mu) + (h(t+1)-mu))) / (1 + F^2);

end

