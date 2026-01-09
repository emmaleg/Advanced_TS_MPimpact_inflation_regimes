function ll = loglik_full(Z, h, S, Phi1, Phi2, A1, A2, sigma2, p, J)
% Full log-likelihood over usable t after building regressors.

[YY, XX] = cf.model.build_regressors(Z, h, p, J);

% Map these rows back to original time indices:
% build_regressors drops first p rows and then possibly drops NaNs from h lags.
% For simplicity: rebuild an index mask similarly.
[T,~] = size(Z);
Xlag = cf.model.make_lag_matrix(Z,p);
Hlags = NaN(T,J+1);
for j=0:J, Hlags(:,j+1)= lagmatrix(h,j); end
Xfull = [ones(T,1), Xlag, Hlags];
start = p+1;
Y0 = Z(start:end,:);
X0 = Xfull(start:end,:);
good = all(isfinite(X0),2) & all(isfinite(Y0),2);
idx = (start:T)';
idx = idx(good);

ll = 0;
for ii=1:numel(idx)
    t = idx(ii);
    z_t = Z(t,:)';
    x_t = Xfull(t,:)';
    if S(t)==1
        ll = ll + cf.model.loglik_obs(z_t, x_t, Phi1, A1, sigma2, h(t));
    else
        ll = ll + cf.model.loglik_obs(z_t, x_t, Phi2, A2, sigma2, h(t));
    end
end
end
