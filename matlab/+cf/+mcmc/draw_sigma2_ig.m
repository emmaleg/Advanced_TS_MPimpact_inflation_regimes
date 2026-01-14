function st = draw_sigma2_ig(st, Z, mconf, pconf)
% Draw sigma2_j from inverse-gamma.

[T0,n] = size(Z);

Xlag = cf.model.make_lag_matrix(Z,mconf.p);
Hlags = NaN(T0,mconf.J+1);
for j=0:mconf.J, Hlags(:,j+1)=lagmatrix(st.h,j); end
Xfull = [ones(T0,1), Xlag, Hlags];

start = mconf.p+1;
good = all(isfinite(Xfull(start:end,:)),2) & all(isfinite(Z(start:end,:)),2);
idx = (start:T0)'; idx = idx(good);

% IMPORTANT: exclude undefined-regime periods (t<=d)
if isfield(st,'Svalid') && ~isempty(st.Svalid)
    idx = idx(st.Svalid(idx));
end

ss = zeros(n,1);
cnt = 0;

for ii=1:numel(idx)
    t = idx(ii);
    zt = Z(t,:)';
    xt = Xfull(t,:)';

    if st.S(t)==1
        eps = zt - (st.Phi1*xt);
        u = st.A1 * eps;
    else
        eps = zt - (st.Phi2*xt);
        u = st.A2 * eps;
    end

    lam = exp(st.h(t));
    w = u / sqrt(lam);
    ss = ss + (w.^2);
    cnt = cnt + 1;
end

for j=1:n
    a = pconf.sigma2.a0 + 0.5*cnt;
    b = pconf.sigma2.b0 + 0.5*ss(j);
    st.sigma2(j) = cf.mcmc.rand_ig(a,b);
end
end
