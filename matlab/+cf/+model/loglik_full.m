function ll = loglik_full(Z, S, Svalid, st, mconf)
% Full log-likelihood sum_t log p(Z_t | X_t, regime, params)
% skipping invalid periods (t<=d) and rows with NaNs.

[T,~] = size(Z);

% build Xfull = [1, lags(Z), lags(h)]
Xlag = cf.model.make_lag_matrix(Z, mconf.p);

Hlags = NaN(T, mconf.J+1);
for j=0:mconf.J
    Hlags(:,j+1) = lagmatrix(st.h, j);
end
Xfull = [ones(T,1), Xlag, Hlags];

start = mconf.p + 1;
ll = 0;

for t = start:T
    if ~Svalid(t), continue; end
    if any(~isfinite(Z(t,:))) || any(~isfinite(Xfull(t,:))), continue; end

    z = Z(t,:)';
    x = Xfull(t,:)';

    if S(t)==1
        ll = ll + cf.model.loglik_obs(z, x, st.Phi1, st.A1, st.sigma2, st.h(t));
    else
        ll = ll + cf.model.loglik_obs(z, x, st.Phi2, st.A2, st.sigma2, st.h(t));
    end
end

end
