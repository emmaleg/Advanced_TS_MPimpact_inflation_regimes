function st = draw_mu_normal(st, pconf)
% Draw mu in h_t - mu = F(h_{t-1}-mu)+eta_t, eta~N(0,Q).

h = st.h;
F = st.F;
Q = st.Q;

% Regression: (h_t - F h_{t-1}) = (1-F)*mu + eta_t
y = h(2:end) - F*h(1:end-1);
X = (1-F)*ones(numel(y),1);

% Prior mu ~ N(mu0, V0)
V0 = pconf.lambda.V0;
mu0 = pconf.lambda.mu0;

Vpost = 1 / ( (X'*X)/Q + 1/V0 );
mpost = Vpost * ( (X'*y)/Q + mu0/V0 );

st.mu = mpost + sqrt(Vpost)*randn();
end
