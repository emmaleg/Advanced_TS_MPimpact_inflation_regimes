function st = draw_Q_ig(st, pconf)
% Q | h, mu, F is Inverse-Gamma with stationary initial term.

h  = st.h(:);
mu = st.mu;
F  = min(max(st.F, 0.001), 0.9999);

% innovations:
e1 = sqrt(1 - F^2) * (h(1) - mu);
e  = h(2:end) - mu - F*(h(1:end-1) - mu);

ss = e1^2 + sum(e.^2);

a0 = pconf.lambda.aQ0;
b0 = pconf.lambda.bQ0;

T = numel(h);
a = a0 + 0.5*T;
b = b0 + 0.5*ss;

st.Q = cf.mcmc.rand_ig(a,b);

end

