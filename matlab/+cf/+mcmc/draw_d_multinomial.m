function st = draw_d_multinomial(st, Z, mconf, pconf)
%CF.MCMC.DRAW_D_MULTINOMIAL  Draw delay d from discrete posterior over {1,...,dmax}.
% Uses likelihood loglik_full(Z,S,Svalid,st,mconf) and uniform prior on d.

Pi   = Z(:, mconf.inflation_index_in_Z);
dmax = pconf.d.dmax;

logw = -Inf(dmax,1);

for d = 1:dmax
    [Sd, Sd_valid] = cf.model.regime_indicator(Pi, st.Pstar, d);

    % Optional: ensure enough obs per regime on valid periods
    n_low  = sum(Sd(Sd_valid) == 1);
    n_high = sum(Sd(Sd_valid) == 0);
    if n_low < 30 || n_high < 30
        logw(d) = -Inf;
        continue;
    end

    logw(d) = cf.model.loglik_full(Z, Sd, Sd_valid, st, mconf);
end

mx = max(logw);
if ~isfinite(mx)
    error('draw_d_multinomial: all candidate d had -Inf log weight. Check Svalid and likelihood.');
end

w = exp(logw - mx);
w = w / sum(w);

u = rand();
cdf = cumsum(w);
dnew = find(u <= cdf, 1, 'first');

st.d = dnew;
[st.S, st.Svalid] = cf.model.regime_indicator(Pi, st.Pstar, st.d);

end
