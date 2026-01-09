function st = draw_d_multinomial(st, Z, mconf, pconf)
% Draw delay d from multinomial over {1,...,dmax}. 

Pi = Z(:, mconf.inflation_index_in_Z);

logw = -Inf(mconf.dmax,1);
for d=1:mconf.dmax
    Sd = cf.model.regime_indicator(Pi, st.Pstar, d);
    % (uniform prior on d) => weights proportional to likelihood 
    logw(d) = cf.model.loglik_full(Z, st.h, Sd, st.Phi1, st.Phi2, st.A1, st.A2, st.sigma2, mconf.p, mconf.J);
end

% normalize
mx = max(logw);
w = exp(logw - mx);
w = w / sum(w);

u = rand();
cdf = cumsum(w);
dnew = find(u <= cdf, 1, 'first');

st.d = dnew;
st.S = cf.model.regime_indicator(Pi, st.Pstar, st.d);

end
