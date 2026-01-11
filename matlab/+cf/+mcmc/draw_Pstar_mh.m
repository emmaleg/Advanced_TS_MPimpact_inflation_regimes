function st = draw_Pstar_mh(st, Z, mconf, pconf, mcmc, iter)
% RW Metropolis for P* with simple adaptation (Haario-like in 1D).
% Prior: Truncated Normal N(mu,sd^2) truncated to [Pmin,Pmax] (Table 2).

Pcur = st.Pstar;
sd_prop = st.Pstar_prop_sd;

Pcan = Pcur + sd_prop*randn();

% --- prior truncation (hard bounds) ---
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

% --- likelihood under current/candidate regimes ---
ll_cur = cf.model.loglik_full(Z, st.h, st.S,  st.Phi1, st.Phi2, st.A1, st.A2, st.sigma2, mconf.p, mconf.J);
ll_can = cf.model.loglik_full(Z, st.h, Scan,  st.Phi1, st.Phi2, st.A1, st.A2, st.sigma2, mconf.p, mconf.J);

% --- Truncated Normal prior contribution (normalization cancels in ratio) ---
muP = pconf.Pstar.mu;
sdP = pconf.Pstar.sd;

lp_cur = -0.5 * ((Pcur - muP)/sdP)^2;
lp_can = -0.5 * ((Pcan - muP)/sdP)^2;

log_acc = (ll_can + lp_can) - (ll_cur + lp_cur);
acc = exp(min(0, log_acc));

if rand() < acc
    st.Pstar = Pcan;
    st.S     = Scan;
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