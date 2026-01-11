function st = draw_d_multinomial(st, Z, mconf, pconf)
% Draw delay d from multinomial over {1,...,dmax}. 

Pi = Z(:, mconf.inflation_index_in_Z);
dmax = pconf.d.dmax;

logw = -Inf(dmax,1);
for d = 1:dmax
    [Sd, Sd_valid] = cf.model.regime_indicator(Pi, st.Pstar, d);

    % weights ∝ likelihood (uniform prior on d)
    logw(d) = cf.model.loglik_full(Z, st.h, Sd, st.Phi1, st.Phi2, ...
                                   st.A1, st.A2, st.sigma2, mconf.p, mconf.J, Sd_valid);
end

% normalize
mx = max(logw);
if ~isfinite(mx)
    error('draw_d_multinomial: all log-weights are -Inf. Check likelihood inputs.');
end

w = exp(logw - mx);
w = w / sum(w);

u = rand();
cdf = cumsum(w);
dnew = find(u <= cdf, 1, 'first');

st.d = dnew;
[st.S, st.Svalid] = cf.model.regime_indicator(Pi, st.Pstar, st.d);

end
