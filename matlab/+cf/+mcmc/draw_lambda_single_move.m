function st = draw_lambda_single_move(st, Z, mconf, pconf, mcmc)
% Block update for h_t = log(lambda_t) using a Gaussian approximation
% to log-chi2 and FFBS (Kalman simulation smoother).
%
% Key idea:
%   s_t = sum_j (u_{j,t}^2 / sigma2_j)  ~  lambda_t * chi2(n)
% => log(s_t) = h_t + log(chi2(n))  approximately Gaussian for n>=5.

[T,n] = size(Z);

% Build Xfull from current h (we treat it as fixed in this block; standard Gibbs approx)
Xlag = cf.model.make_lag_matrix(Z, mconf.p);
Hlags = NaN(T, mconf.J+1);
for j=0:mconf.J
    Hlags(:,j+1) = lagmatrix(st.h, j);
end
Xfull = [ones(T,1), Xlag, Hlags];

start = mconf.p + 1;

% Moments of log(chi2_n)
m_v = psi(n/2) + log(2);     % E[log chi2_n]
v_v = psi(1, n/2);           % Var[log chi2_n] (trigamma)
R   = max(v_v, 1e-8);        % measurement variance

% Construct pseudo-observation y_t ≈ h_t + e_t
y = NaN(T,1);

for t = start:T
    if isfield(st,'Svalid') && ~isempty(st.Svalid) && ~st.Svalid(t), continue; end
    if any(~isfinite(Z(t,:))) || any(~isfinite(Xfull(t,:))), continue; end

    x = Xfull(t,:)';
    z = Z(t,:)';

    if st.S(t)==1
        eps = z - (st.Phi1 * x);
        u   = st.A1 * eps;
    else
        eps = z - (st.Phi2 * x);
        u   = st.A2 * eps;
    end

    s = sum( (u.^2) ./ st.sigma2(:) );
    s = max(s, 1e-12);

    y(t) = log(s) - m_v;
end

% State equation: h_t = c + F*h_{t-1} + eta_t
F  = min(max(st.F, 0.001), 0.9999);
Q  = max(st.Q, 1e-10);
mu = st.mu;
c  = (1 - F) * mu;

% Stationary initial distribution
a1 = mu;
P1 = Q / max(1 - F^2, 1e-8);

% Kalman filter
a_pred = zeros(T,1); P_pred = zeros(T,1);
a_filt = zeros(T,1); P_filt = zeros(T,1);

a_pred(1) = a1; P_pred(1) = P1;

if isfinite(y(1))
    v  = y(1) - a_pred(1);
    S  = P_pred(1) + R;
    K  = P_pred(1) / S;
    a_filt(1) = a_pred(1) + K*v;
    P_filt(1) = (1 - K)*P_pred(1);
else
    a_filt(1) = a_pred(1);
    P_filt(1) = P_pred(1);
end

for t = 2:T
    a_pred(t) = c + F*a_filt(t-1);
    P_pred(t) = F^2 * P_filt(t-1) + Q;

    if isfinite(y(t))
        v  = y(t) - a_pred(t);
        S  = P_pred(t) + R;
        K  = P_pred(t) / S;
        a_filt(t) = a_pred(t) + K*v;
        P_filt(t) = (1 - K)*P_pred(t);
    else
        a_filt(t) = a_pred(t);
        P_filt(t) = P_pred(t);
    end
end

% Backward simulation (Carter-Kohn)
h_draw = zeros(T,1);
h_draw(T) = a_filt(T) + sqrt(max(P_filt(T),1e-12))*randn();

for t = T-1:-1:1
    denom = max(P_pred(t+1), 1e-12);
    Jt    = P_filt(t) * F / denom;

    m = a_filt(t) + Jt*(h_draw(t+1) - a_pred(t+1));
    V = P_filt(t) - (Jt^2)*P_pred(t+1);
    V = max(V, 1e-12);

    h_draw(t) = m + sqrt(V)*randn();
end

st.h = h_draw;

% For book-keeping : counting the nbr of updates
if ~isfield(st,'accept') || isempty(st.accept)
    st.accept = struct();
end
if ~isfield(st.accept,'h_updates')
    st.accept.h_updates = 0;
end
st.accept.h_updates = st.accept.h_updates + 1;

end

