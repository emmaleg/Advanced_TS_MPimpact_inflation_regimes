function st = draw_mu_normal(st, pconf)
% mu | h, F, Q is Normal (with stationary initial density)

h  = st.h(:);
T  = numel(h);
F  = min(max(st.F, 0.001), 0.9999);
Q  = max(st.Q, 1e-10);

mu0 = pconf.lambda.mu0;
V0  = pconf.lambda.V0;

% Using:
% h1 ~ N(mu, Q/(1-F^2))
% ht = mu + F(h_{t-1}-mu) + eta  => (ht - F h_{t-1}) = (1-F) mu + eta

y  = h(2:end) - F*h(1:end-1);
s1 = (1 - F^2) / Q;
s2 = (1 - F)^2 * (T-1) / Q;

prec = (1/V0) + s1 + s2;

rhs = (mu0/V0) + s1*h(1) + ((1 - F)/Q)*sum(y);

Vpost = 1/prec;
mpost = Vpost * rhs;

st.mu = mpost + sqrt(max(Vpost,1e-12))*randn();

end

