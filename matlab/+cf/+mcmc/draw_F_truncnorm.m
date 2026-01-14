function st = draw_F_truncnorm(st, pconf)
% F | h, mu, Q is Normal truncated to (0,1)

h  = st.h(:);
mu = st.mu;
Q  = max(st.Q, 1e-10);

x = h(1:end-1) - mu;
y = h(2:end)   - mu;

F0  = pconf.lambda.F0;
VF0 = pconf.lambda.VF0;

prec = (1/VF0) + (sum(x.^2)/Q);
Vpost = 1/prec;
mpost = Vpost * (F0/VF0 + sum(x.*y)/Q);

% Truncated normal by rejection (works because posterior tight)
for k=1:5000
    Fcan = mpost + sqrt(max(Vpost,1e-12))*randn();
    if Fcan > 0 && Fcan < 1
        st.F = Fcan;
        return;
    end
end

% fallback: clamp (rare)
st.F = min(max(mpost, 1e-3), 0.999);

end

