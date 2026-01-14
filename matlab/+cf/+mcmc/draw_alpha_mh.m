function st = draw_alpha_mh(st, Z, mconf, pconf, mcmc)
% Draw alpha_i (i=1,2) via MH.

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

% IMPORTANT: exclude undefined-regime periods (t<=d)
if isfield(st,'Svalid') && ~isempty(st.Svalid)
    idx_all = idx_all(st.Svalid(idx_all));
end

% --- SA,sA such that vec(A)=SA*alpha + sA ---
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

        if reg==1
            use = (st2.S(idx_all)==1);
        else
            use = (st2.S(idx_all)==0);
        end
        id = idx_all(use);

        if numel(id) < 10
            % too few obs -> skip
            return;
        end

        % Build V_alpha
        Vinv = OmegaInv;

        Sigma = diag(st2.sigma2(:));

        for ii=1:numel(id)
            t = id(ii);

            eps_t = Z(t,:)' - (Phi * Xfull(t,:)');
            E     = kron(eps_t', eye(n));         % n x (n^2)

            e_tilde = E * sA;                     % n x 1
            x_tilde = -E * SA;                    % n x na

            lam  = exp(st2.h(t));
            Hinv = diag(1./(lam * st2.sigma2(:)));

            Vinv = Vinv + (x_tilde' * Hinv * x_tilde);
        end

        V = inv(symm(Vinv));

        % Current A
        Acur = cf.id.alpha_to_A(alpha_cur);

        % If identification is enforced, ensure current point is in support
        if isfield(mconf,'id') && isfield(mconf.id,'enforce_sign_zero') && mconf.id.enforce_sign_zero
            tol = mconf.id.zero_tol;
            ok_cur = cf.id.check_impact_restrictions(Acur, st2.sigma2, exp(mean(st2.h)), tol);
            if ~ok_cur
                % repair by drawing a feasible alpha from the PRIOR
                [alpha_cur, Acur] = draw_feasible_alpha(pconf, st2.sigma2, exp(mean(st2.h)), tol, na, 5000);
                if reg==1
                    st2.alpha1 = alpha_cur; st2.A1 = Acur;
                else
                    st2.alpha2 = alpha_cur; st2.A2 = Acur;
                end
            end
        end

        % RW proposal
        c = mcmc.alpha_prop_scale;
        L = chol_psd(symm(V));
        alpha_can = alpha_cur + sqrt(c) * (L * randn(na,1));
        Acan = cf.id.alpha_to_A(alpha_can);

        % Optional sign/zero restrictions on the candidate
        if isfield(mconf,'id') && isfield(mconf.id,'enforce_sign_zero') && mconf.id.enforce_sign_zero
            tol = mconf.id.zero_tol;
            ok1 = cf.id.check_impact_restrictions(Acan, st2.sigma2, exp(mean(st2.h)), tol);
            if ~ok1
                if reg==1
                    st2.accept.alpha1_trials = st2.accept.alpha1_trials + 1;
                else
                    st2.accept.alpha2_trials = st2.accept.alpha2_trials + 1;
                end
                return;
            end
        end

        lp_cur = logpost_alpha(alpha_cur, Acur, Phi, id, Z, Xfull, st2, pconf);
        lp_can = logpost_alpha(alpha_can, Acan, Phi, id, Z, Xfull, st2, pconf);

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
V = (V + V')/2;
jitter = 1e-10;
for it=1:8
    [L,p] = chol(V + jitter*eye(size(V)), 'lower');
    if p==0, return; end
    jitter = jitter * 10;
end
L = chol(V + 1e-4*eye(size(V)), 'lower');
end

function [alpha, A] = draw_feasible_alpha(pconf, sigma2, lambda0, tol, na, maxTries)
mu = pconf.alpha.mu(:);
Om = (pconf.alpha.Omega + pconf.alpha.Omega')/2;

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

error("draw_alpha_mh: cannot find feasible alpha after %d tries.", maxTries);
end
