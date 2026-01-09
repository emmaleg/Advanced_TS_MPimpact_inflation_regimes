function st = draw_lambda_single_move(st, Z, mconf, pconf, mcmc)
% Single-move MH updates for h_t = ln lambda_t.
% Appendix A discusses single-move due to nonlinearity (lambda enters mean and variance).
%
% Here: independence proposal from conditional AR(1) prior (A.16-A.18), accept by local likelihood ratio.
% We recompute a LOCAL likelihood window [t, t+J] because h_t enters mean with lags up to J.

[T0,~] = size(Z);

% Precompute regressors container
Xlag = cf.model.make_lag_matrix(Z,mconf.p);
Hlags = NaN(T0,mconf.J+1);
for j=0:mconf.J, Hlags(:,j+1)=lagmatrix(st.h,j); end
Xfull = [ones(T0,1), Xlag, Hlags];

start = mconf.p+1;

for t=start:T0
    % Proposal from conditional prior
    [mprior, vprior] = cond_ar1_prior(t, st.h, st.mu, st.F, st.Q);
    h_can = mprior + sqrt(vprior)*randn();

    % local indices affected by h_t through mean: tau in [t, t+J]
    tau1 = t;
    tau2 = min(T0, t+mconf.J);

    ll_cur = local_ll(Z, Xfull, st, tau1, tau2, st.h(t));
    h_old = st.h(t);
    st.h(t) = h_can;
    % update Hlags for impacted rows (cheap update)
    for j=0:mconf.J
        if (t+j) <= T0
            col = 2 + size(Xlag,2) + j; 
            Xfull(t+j, col) = st.h(t);
        end
    end
    % update Hlags for impacted rows (cheap update)
    %for j=0:mconf.J
    %    if (t+j) <= T0
    %        Xfull(t+j, 2 + size(Xlag,2) + (j+1) ) = st.h(t); % column for h_{tau-j} at tau=t+j, j=(tau-t)
    %    end
    %end
    ll_can = local_ll(Z, Xfull, st, tau1, tau2, h_can);

    % acceptance (proposal symmetric in conditional-prior space? It's independence but same both sides if using same prior)
    st.accept.h_trials = st.accept.h_trials + 1;
    acc = min(1, exp(ll_can - ll_cur));
    if rand() < acc
        st.accept.h = st.accept.h + 1;
    else
        st.h(t) = h_old;
        % restore Xfull impacted rows
        for j=0:mconf.J
            if (t+j) <= T0
                col = 2 + size(Xlag,2) + j;   % <-- CORRECT
                Xfull(t+j, col) = h_old;
            end
        end
    end
end
end

function [m,v] = cond_ar1_prior(t, h, mu, F, Q)
% Conditional prior for interior points (A.17-A.18). Endpoints handled simply.
T = numel(h);
if t==1
    m = mu + F*(h(2)-mu); v = Q; return;
elseif t==T
    m = mu + F*(h(T-1)-mu); v = Q; return;
else
    m = mu + (F*((h(t-1)-mu) + (h(t+1)-mu))) / (1+F^2);
    v = Q / (1+F^2);
end
end


function ll = local_ll(Z, Xfull, st, t1, t2, h_t)
% Compute log-likelihood contributions for tau=t1..t2 (only if regressors finite and tau>=p+1)
ll = 0;
for tau=t1:t2
    if any(~isfinite(Xfull(tau,:))) || any(~isfinite(Z(tau,:)))
        continue;
    end

    z = Z(tau,:)';
    x = Xfull(tau,:)';

    if st.S(tau)==1
        ll = ll + cf.model.loglik_obs(z, x, st.Phi1, st.A1, st.sigma2, st.h(tau));
    else
        ll = ll + cf.model.loglik_obs(z, x, st.Phi2, st.A2, st.sigma2, st.h(tau));
    end
end
end

