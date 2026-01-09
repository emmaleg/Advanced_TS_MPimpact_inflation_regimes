function st = draw_alpha_mh(st, Z, mconf, pconf, mcmc)
% MH draw for alpha_i using A.8 mapping and SA,sA regression trick. 

[SA, sA] = cf.id.build_SA_sA();
OmegaA_inv = inv(pconf.alpha.Omega);

% usable indices for regression part
[T0,~] = size(Z);
Xlag = cf.model.make_lag_matrix(Z,mconf.p);
Hlags = NaN(T0,mconf.J+1);
for j=0:mconf.J, Hlags(:,j+1)=lagmatrix(st.h,j); end
Xfull = [ones(T0,1), Xlag, Hlags];
start = mconf.p+1;
good = all(isfinite(Xfull(start:end,:)),2) & all(isfinite(Z(start:end,:)),2);
idx = (start:T0)'; idx = idx(good);

% Draw alpha1 and alpha2 separately
[st.alpha1, st.A1, st.accept.alpha1, st.accept.alpha1_trials] = ...
    one_alpha_update(1, st.alpha1, st.Phi1, st.A1, idx, Z, Xfull, st.S, st.h, st.sigma2, SA, sA, OmegaA_inv, pconf, mconf, mcmc, st.accept.alpha1, st.accept.alpha1_trials);

[st.alpha2, st.A2, st.accept.alpha2, st.accept.alpha2_trials] = ...
    one_alpha_update(2, st.alpha2, st.Phi2, st.A2, idx, Z, Xfull, st.S, st.h, st.sigma2, SA, sA, OmegaA_inv, pconf, mconf, mcmc, st.accept.alpha2, st.accept.alpha2_trials);

end

function [alpha_new, A_new, acc_cnt, trial_cnt] = one_alpha_update(reg, alpha_cur, Phi, Acur, idx, Z, Xfull, S, h, sigma2, SA, sA, OmegaA_inv, pconf, mconf, mcmc, acc_cnt, trial_cnt)
% Build residuals for times in this regime
if reg==1
    use = (S(idx)==1);
else
    use = (S(idx)==0);
end
id = idx(use);
if numel(id)<20
    alpha_new = alpha_cur;
    A_new = Acur;
    return;
end

n = size(Z,2);
Sigma = diag(sigma2);

% Compute posterior approx (A.11-A.12) for proposal covariance
Vinv = OmegaA_inv;
b = OmegaA_inv * pconf.alpha.mu;

for ii=1:numel(id)
    t = id(ii);
    zt = Z(t,:)';
    xt = Xfull(t,:)';
    eps = zt - (Phi*xt);

    % Build x_tilde and e_tilde
    K = kron(eps', eye(n)); % (eps' ⊗ I)
    etilde = K * sA;        % 64x1
    xtilde = -K * SA;       % 64x22

    % H_t^{-1} with standardization by lambda: u/sqrt(lambda) ~ N(0,Sigma)
    lam = exp(h(t));
    Hinv = (1/lam) * inv(Sigma); % for u ~ N(0,lam*Sigma) (approx)
    % But regression error is in R^n "stacked": use kron(Hinv, I) approx.
    % Practical simplification: weight each row-block by Hinv via kron.
    W = kron(Hinv, eye(8)); 
    % Simpler: since etilde corresponds to vec(A*eps), we weight by inv(lam*Sigma) on each equation component.
    % Use kron(inv(lam*Sigma), I_n?) as coarse. Here we approximate with scalar weight:
    w = 1/lam;

    Vinv = Vinv + w*(xtilde' * xtilde);
    b    = b    + w*(xtilde' * etilde);
end

V = inv((Vinv+Vinv')/2);
mu = V * b;

% Proposal
trial_cnt = trial_cnt + 1;
L = chol((V+V')/2,'lower');
alpha_can = alpha_cur + mcmc.alpha_prop_scale*(L*randn(size(alpha_cur)));

Acan = cf.id.alpha_to_A(alpha_can);

% Optional sign/zero restrictions
if mconf.id.enforce_sign_zero
    ok = check_impact_restrictions(Acan, sigma2, exp(mean(h)), mconf.id.zero_tol);
    if ~ok
        alpha_new = alpha_cur; A_new = Acur;
        return;
    end
end

% MH acceptance ratio on posterior (likelihood via u=A*eps + Jacobian)
lp_cur = logpost_alpha(alpha_cur, Acur, Phi, reg, idx, Z, Xfull, S, h, Sigma, pconf);
lp_can = logpost_alpha(alpha_can, Acan, Phi, reg, idx, Z, Xfull, S, h, Sigma, pconf);

acc = min(1, exp(lp_can - lp_cur));
if rand()<acc
    alpha_new = alpha_can; A_new = Acan;
    acc_cnt = acc_cnt + 1;
else
    alpha_new = alpha_cur; A_new = Acur;
end

end

function lp = logpost_alpha(alpha, A, Phi, reg, idx, Z, Xfull, S, h, Sigma, pconf)
% Log posterior up to constant: loglik + logprior
lp = 0;

% prior
da = alpha - pconf.alpha.mu;
lp = lp - 0.5 * (da' * (pconf.alpha.Omega \ da));

% likelihood
if reg==1
    use = (S(idx)==1);
else
    use = (S(idx)==0);
end
id = idx(use);

logdetA = log(abs(det(A)));
logdetSigma = sum(log(diag(Sigma)));

for ii=1:numel(id)
    t = id(ii);
    zt = Z(t,:)';
    xt = Xfull(t,:)';
    eps = zt - (Phi*xt);
    u = A * eps;

    lam = exp(h(t));
    logdetH = size(Z,2)*log(lam) + logdetSigma;
    quad = u' * ((1/lam) * (Sigma \ u));

    lp = lp + logdetA - 0.5*logdetH - 0.5*quad;
end
end
