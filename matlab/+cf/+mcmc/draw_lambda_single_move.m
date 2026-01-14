function st = draw_lambda_single_move(st, Z, mconf, pconf, mcmc)
%CF.MCMC.DRAW_LAMBDA_SINGLE_MOVE  Draw h_t = log(lambda_t) via Jacquier et al. (1994) single-move.
%
% This implements Appendix A, Block 6 of Canova & Perez Forero (2024)
% (eqs. (A.15)–(A.24)) using the convexity bound for exp(-h_t).  
%
% Notation in the code:
%   - h_t is st.h(t) = log(lambda_t)
%   - AR(1) prior: h_t - mu = F (h_{t-1}-mu) + eta_t, eta_t ~ N(0,Q)   
%   - For each t, compute conditional prior moments (h*_t, v2) as in (A.17)-(A.18)
%   - Let ss_t = sum_j u_{j,t}^2 / sigma_j^2, where u_t are the structural residuals.
%     In a multivariate system with common scale lambda_t, the log-likelihood contribution is:
%         ln f(u_t|h_t) = const - (n/2) h_t - (1/2) ss_t * exp(-h_t)
%     The bound and the accept/reject step generalize directly.
%
% IMPORTANT (consistency with the Appendix derivation):
%   The bound (A.20) assumes ss_t is treated as fixed data when drawing h_t.  
%   We therefore compute ss_t once at the beginning of the block from the current state and
%   keep it fixed while updating the whole path h_1..h_T.

if ~isfield(mcmc, 'lambda_max_tries'); mcmc.lambda_max_tries = 50; end

[T, n] = size(Z);

% -------------------------------------------------------------------------
% Precompute ss_t = sum_j u_{j,t}^2 / sigma_j^2 using the current state.
% Residual definition follows Appendix A.14 (u_{j,t} = A1 e~_{1,t} S_t + A2 e~_{2,t} (1-S_t)).
% -------------------------------------------------------------------------
ss = compute_ss_path(st, Z, mconf); % T x 1, NaN where undefined

mu = st.mu;
F  = st.F;
Q  = st.Q;

% Update each h_t sequentially (single-move).
for t = 1:T
    % Conditional prior for h_t given neighbors under AR(1)
    [h_star, v2] = cond_prior_h(t, st.h, mu, F, Q);

    % If ss_t is missing (early periods / missing regressors), draw from the conditional prior.
    if ~isfinite(ss(t))
        st.h(t) = h_star + sqrt(max(v2, 1e-12))*randn();
        continue
    end

    % Proposal mean mu_t in (A.23) generalized to multivariate common-scale case:
    %   mu_t = h*_t + (v2/2) * ( ss_t * exp(-h*_t) - n )
    % (equivalently, n * (ybar_t^2 exp(-h*_t) - 1) if ybar^2 = ss/n).   
    mu_t = h_star + 0.5*v2*(ss(t)*exp(-h_star) - n);

    % Accept/Reject using envelope implied by (A.20)-(A.24).            
    accepted = false;
    tries = 0;
    while ~accepted
        tries = tries + 1;
        st.accept.h_trials = st.accept.h_trials + 1;

        h_can = mu_t + sqrt(max(v2, 1e-12))*randn();

        % log f*(.) and log g*(.) (multivariate generalization)
        %   ln f* = -(n/2) h - 0.5 ss exp(-h)
        %   ln g* = -(n/2) h - 0.5 ss exp(-h*) (1 + h* - h)
        ln_f = -0.5*n*h_can - 0.5*ss(t)*exp(-h_can);
        ln_g = -0.5*n*h_can - 0.5*ss(t)*exp(-h_star)*(1 + h_star - h_can);

        log_alpha = ln_f - ln_g; % <= 0 by convexity
        if log(rand()) < min(0, log_alpha)
            st.h(t) = h_can;
            st.accept.h = st.accept.h + 1;
            accepted = true;
        else
            % True rejection sampling would keep drawing until accepted.
            % To avoid pathological infinite loops in rare numerical cases, cap the tries.
            if tries >= mcmc.lambda_max_tries
                % Fallback: keep the previous value (MH-style). This should almost never trigger.
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

function ss = compute_ss_path(st, Z, mconf)
% Compute ss_t = sum_j u_{j,t}^2 / sigma_j^2 for t=1..T using the current state.

[T, n] = size(Z);
ss = NaN(T,1);

% Build regressors X_t = [1, Z_{t-1}..Z_{t-p}, h_t..h_{t-J}] (rows aligned with Z).
Xlag = cf.model.make_lag_matrix(Z, mconf.p);
Hlags = NaN(T, mconf.J+1);
for j = 0:mconf.J
    Hlags(:, j+1) = lagmatrix(st.h, j);
end
Xfull = [ones(T,1), Xlag, Hlags];

% Structural matrices for both regimes
A1 = st.A1;
A2 = st.A2;

% Compute u_t for all t where regressors are finite and regime indicator is valid.
for t = 1:T
    if ~st.Svalid(t), continue; end
    if ~all(isfinite(Xfull(t,:))) || ~all(isfinite(Z(t,:)))
        continue
    end

    if st.S(t)==1
        Phi = st.Phi1;
        A   = A1;
    else
        Phi = st.Phi2;
        A   = A2;
    end

    eps_t = Z(t,:)' - Phi*Xfull(t,:)';
    u_t   = A*eps_t;

    % ss_t = sum_j u_{j,t}^2 / sigma_j^2
    ss(t) = sum((u_t.^2) ./ st.sigma2(:));
end

end

function [m, v] = cond_prior_h(t, h, mu, F, Q)
% Conditional prior moments for h_t given neighbors under AR(1) with drift.
%
% h_t = (1-F)*mu + F*h_{t-1} + eta_t, eta_t ~ N(0,Q)

T = numel(h);

% Stationary variance for the initial state
V0 = Q / max(1e-12, (1 - F^2));

if t==1
    % h1 | h2
    if T==1
        m = mu;
        v = V0;
        return
    end
    % Prior: h1 ~ N(mu, V0); likelihood from transition to h2.
    v = 1 / (1/V0 + (F^2)/Q);
    m = v * (mu/V0 + (F*(h(2) - (1-F)*mu))/Q);
    return
end

if t==T
    % hT | h_{T-1}
    m = (1-F)*mu + F*h(T-1);
    v = Q;
    return
end

% Interior: h_t | h_{t-1}, h_{t+1}
v = Q / (1 + F^2);
m = mu + (F*((h(t-1)-mu) + (h(t+1)-mu))) / (1 + F^2);

end

