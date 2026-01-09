function st = draw_Q_ig(st, pconf)
% Draw Q from inverse-gamma, based on eta_t = (h_t-mu)-F(h_{t-1}-mu). 

h = st.h; mu = st.mu; F = st.F;

eta = (h(2:end)-mu) - F*(h(1:end-1)-mu);
ss = sum(eta.^2);

a = pconf.lambda.aQ0 + 0.5*numel(eta);
b = pconf.lambda.bQ0 + 0.5*ss;

st.Q = cf.mcmc.rand_ig(a,b);
end
