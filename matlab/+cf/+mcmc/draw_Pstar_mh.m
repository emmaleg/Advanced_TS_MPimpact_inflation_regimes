function st = draw_Pstar_mh(st, Z, mconf, pconf, mcmc, iter)
%CF.MCMC.DRAW_PSTAR_MH  Random-walk Metropolis update for P*.
% Prior: Truncated Normal N(mu,sd^2) truncated to [Pmin,Pmax].
% Likelihood uses loglik_full(Z,S,Svalid,st,mconf) and excludes t<=d via Svalid.

Pcur    = st.Pstar;
sd_prop = st.Pstar_prop_sd;

% Propose
Pcan = Pcur + sd_prop*randn();

% Hard truncation prior support
if Pcan <= pconf.Pstar.Pmin || Pcan >= pconf.Pstar.Pmax
    st.Pstar_hist(iter) = Pcur;
    return;
end

% Regime indicator uses inflation_{t-d}
Pi = Z(:, mconf.inflation_index_in_Z);

% Current regimes (should already be consistent in st, but recompute Svalid safely)
[Scur, Svalid] = cf.model.regime_indicator(Pi, Pcur, st.d);
% Candidate regimes use SAME d => same Svalid by construction
[Scan, ~]      = cf.model.regime_indicator(Pi, Pcan, st.d);

% Require enough observations in each regime ON VALID PERIODS ONLY
n_low  = sum(Scan(Svalid) == 1);
n_high = sum(Scan(Svalid) == 0);
if n_low < 30 || n_high < 30
    st.Pstar_hist(iter) = Pcur;
    return;
end

% Likelihoods (new signature)
ll_cur = cf.model.loglik_full(Z, Scur, Svalid, st, mconf);
ll_can = cf.model.loglik_full(Z, Scan, Svalid, st, mconf);

% Truncated Normal prior (normalizing constant cancels in MH ratio)
muP = pconf.Pstar.mu;
sdP = pconf.Pstar.sd;

lp_cur = -0.5 * ((Pcur - muP)/sdP)^2;
lp_can = -0.5 * ((Pcan - muP)/sdP)^2;

% MH acceptance
st.accept.Pstar_trials = st.accept.Pstar_trials + 1;

log_acc = (ll_can + lp_can) - (ll_cur + lp_cur);
acc = exp(min(0, log_acc));

if rand() < acc
    st.Pstar  = Pcan;
    st.S      = Scan;
    st.Svalid = Svalid;
    st.accept.Pstar = st.accept.Pstar + 1;
else
    % keep current
    st.Pstar  = Pcur;
    st.S      = Scur;
    st.Svalid = Svalid;
end

% Store for adaptation
st.Pstar_hist(iter) = st.Pstar;

% Simple adaptive scaling (optional)
if isfield(mcmc,'Pstar_adapt_start') && iter > mcmc.Pstar_adapt_start
    ph = st.Pstar_hist(1:iter);
    ph = ph(isfinite(ph));
    if numel(ph) > 50
        v = var(ph);
        st.Pstar_prop_sd = max(1e-6, sqrt((2.38^2) * v + 1e-12));
    end
end

end
