function st = draw_Pstar_mh(st, Z, mconf, pconf, mcmc, iter)
% RW Metropolis for P* with simple adaptation (Haario-like in 1D).

Pcur = st.Pstar;
sd = st.Pstar_prop_sd;

Pcan = Pcur + sd*randn();

% prior truncation
if Pcan <= pconf.Pstar.Pmin || Pcan >= pconf.Pstar.Pmax
    st.Pstar_hist(iter) = Pcur;
    return;
end

% ensure both regimes have at least some obs (practical safeguard)
Pi = Z(:, mconf.inflation_index_in_Z);
Scan = cf.model.regime_indicator(Pi, Pcan, st.d);
if sum(Scan==1) < 30 || sum(Scan==0) < 30
    st.Pstar_hist(iter) = Pcur;
    return;
end

st.accept.Pstar_trials = st.accept.Pstar_trials + 1;

% log posterior ratio: uniform prior cancels inside bounds
ll_cur = cf.model.loglik_full(Z, st.h, st.S, st.Phi1, st.Phi2, st.A1, st.A2, st.sigma2, mconf.p, mconf.J);
ll_can = cf.model.loglik_full(Z, st.h, Scan,   st.Phi1, st.Phi2, st.A1, st.A2, st.sigma2, mconf.p, mconf.J);

acc = min(1, exp(ll_can - ll_cur));
if rand() < acc
    st.Pstar = Pcan;
    st.S = Scan;
    st.accept.Pstar = st.accept.Pstar + 1;
end

% store for adaptation
st.Pstar_hist(iter) = st.Pstar;

% Simple adaptation after some iterations
if iter > mcmc.Pstar_adapt_start
    ph = st.Pstar_hist(1:iter);
    ph = ph(isfinite(ph));
    if numel(ph) > 50
        v = var(ph);
        st.Pstar_prop_sd = max(1e-6, sqrt((2.38^2)*v + 1e-12));
    end
end
end
