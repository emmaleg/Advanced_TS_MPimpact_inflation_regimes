function st = draw_F_truncnorm(st, pconf)
% Draw F from truncated normal on (0,1). 

h = st.h;
mu = st.mu;
Q  = st.Q;

x = (h(1:end-1) - mu);
y = (h(2:end)   - mu);

% Likelihood: y = F x + eta
V = 1 / ( (x'*x)/Q + 1/pconf.lambda.VF0 );
m = V * ( (x'*y)/Q + pconf.lambda.F0/pconf.lambda.VF0 );

Fdraw = cf.mcmc.rand_truncnorm(m, sqrt(V), 0.0001, 0.9999);
st.F = Fdraw;
end
