function st = draw_alpha_mh(st, Z, mconf, pconf, mcmc)
% Draw alpha_i (i=1,2) via MH, following Appendix A (eq. A.9-A.12).
% Proposal: Random-walk N(alpha_cur, c * V_alpha), with V_alpha from A.11.
% Acceptance uses full posterior (likelihood includes log|A|).

[T,n] = size(Z);
na    = numel(pconf.alpha.mu);

% --- Build X_t = [1, lags(Z), lags(h)] for all t ---
Xlag = cf.model.make_lag_matrix(Z, mconf.p);

Hlags = NaN(T, mconf.J+1);
for j=0:mconf.J
    Hlags(:,j+1) = lagmatrix(st.h, j);
end

Xfull = [ones(T,1), Xlag, Hlags];

start = mconf.p + 1;
good  = all(isfinite(Xfull(start:end,:)),2) & all(isfinite(Z(start:end,:)),2);
idx_all = (start:T)'; idx_all = idx_all(good);

% --- SA,sA such that vec(A)=SA*alpha + sA (build from alpha_to_A) ---
[SA, sA] = local_build_SA_sA(na);

OmegaInv = inv(pconf.alpha.Omega);
muA      = pconf.alpha.mu;

% Update both regimes
st = one_regime_alpha_update(st, 1);
st = one_regime_alpha_update(st, 2);

    function st2 = one_regime_alpha_update(st2, reg)
        if reg==1
            alpha_cur = st2.alpha1;
            Phi       = st2.Phi1;
        else
            alpha_cur = st2.alpha2;
            Phi       = st2.Phi2;
        end

        % indices in regime
        if reg==1
            use = (st2.S(idx_all)==1);
        else
            use = (st2.S(idx_all)==0);
        end
        id = idx_all(use);

        if numel(id) < 10
            % too few obs -> skip (debug safeguard)
            return;
        end

        % Build V_alpha (A.11): V = (Omega^{-1} + sum X' H^{-1} X)^{-1}
        Vinv = OmegaInv;
        bvec = OmegaInv * muA;

        Sigma = diag(st2.sigma2(:)); % n x n

        for ii=1:numel(id)
            t = id(ii);

            eps_t = Z(t,:)' - (Phi * Xfull(t,:)');
            E     = kron(eps_t', eye(n));         % n x (n^2)

            e_tilde = E * sA;                     % n x 1
            x_tilde = -E * SA;                    % n x na

            lam  = exp(st2.h(t));
            Hinv = diag(1./(lam * st2.sigma2(:))); % n x n (since H_t = lam*diag(sigma2))

            Vinv = Vinv + (x_tilde' * Hinv * x_tilde);
            bvec = bvec + (x_tilde' * Hinv * e_tilde);
        end

        V = inv(symm(Vinv));

        % RW proposal: alpha_can ~ N(alpha_cur, c*V)
        c = mcmc.alpha_prop_scale;
        L = chol_psd(symm(V));
        alpha_can = alpha_cur + sqrt(c) * (L * randn(na,1));

        Acur = cf.id.alpha_to_A(alpha_cur);
        Acan = cf.id.alpha_to_A(alpha_can);

        % Optional sign/zero restrictions
        if isfield(mconf,'id') && isfield(mconf.id,'enforce_sign_zero') && mconf.id.enforce_sign_zero
            tol = mconf.id.zero_tol;
            ok1 = cf.id.check_impact_restrictions(Acan, st2.sigma2, exp(mean(st2.h)), tol);
            if ~ok1
                % reject immediately
                if reg==1
                    st2.accept.alpha1_trials = st2.accept.alpha1_trials + 1;
                else
                    st2.accept.alpha2_trials = st2.accept.alpha2_trials + 1;
                end
                return;
            end
        end

        % MH accept with full posterior
        lp_cur = logpost_alpha(alpha_cur, Acur, Phi, id, Z, Xfull, st2, pconf);
        lp_can = logpost_alpha(alpha_can, Acan, Phi, id, Z, Xfull, st2, pconf);

        % RW with constant covariance -> symmetric -> q cancels
        acc = min(1, exp(lp_can - lp_cur));

        if reg==1
            st2.accept.alpha1_trials = st2.accept.alpha1_trials + 1;
        else
            st2.accept.alpha2_trials = st2.accept.alpha2_trials + 1;
        end

        if rand() < acc
            if reg==1
                st2.alpha1 = alpha_can;
                st2.A1     = Acan;
                st2.accept.alpha1 = st2.accept.alpha1 + 1;
            else
                st2.alpha2 = alpha_can;
                st2.A2     = Acan;
                st2.accept.alpha2 = st2.accept.alpha2 + 1;
            end
        else
            if reg==1
                st2.A1 = Acur;
            else
                st2.A2 = Acur;
            end
        end
    end
end

function lp = logpost_alpha(alpha, A, Phi, id, Z, Xfull, st, pconf)
% Log posterior up to constant: log prior + sum_t loglik(z_t | ...)
da = alpha - pconf.alpha.mu;
lp = -0.5 * (da' * (pconf.alpha.Omega \ da));

for ii=1:numel(id)
    t = id(ii);
    z = Z(t,:)';
    x = Xfull(t,:)';
    lp = lp + cf.model.loglik_obs(z, x, Phi, A, st.sigma2, st.h(t));
end
end

function [SA, sA] = local_build_SA_sA(na)
% Build SA,sA numerically from cf.id.alpha_to_A, so you don't depend on a separate file.
n = 8;
z = zeros(na,1);
A0 = cf.id.alpha_to_A(z);
sA = A0(:);

SA = zeros(n*n, na);
for k=1:na
    e = zeros(na,1);
    e(k) = 1;
    Ak = cf.id.alpha_to_A(e);
    SA(:,k) = Ak(:) - sA;
end
end

function M = symm(M)
M = (M + M')/2;
end

function L = chol_psd(V)
% robust chol with jitter
V = (V + V')/2;
jitter = 1e-10;
for it=1:8
    [L,p] = chol(V + jitter*eye(size(V)), 'lower');
    if p==0, return; end
    jitter = jitter * 10;
end
% last resort
L = chol(V + 1e-4*eye(size(V)), 'lower');
end

