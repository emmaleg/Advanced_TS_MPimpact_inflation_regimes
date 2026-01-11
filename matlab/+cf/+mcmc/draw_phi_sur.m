function st = draw_phi_sur(st, Z, mconf, pconf)
% Draw Phi1 and Phi2 using multivariate regression with known Omega_bar
% (Omega_t = lambda_t * Omega_bar). Appendix A.6-A.7 structure.

n = size(Z,2);

% Build full X with h-lags
[T0,~] = size(Z);
Xlag = cf.model.make_lag_matrix(Z, mconf.p);
Hlags = NaN(T0, mconf.J+1);
for j=0:mconf.J
    Hlags(:,j+1) = lagmatrix(st.h, j);
end
Xfull = [ones(T0,1), Xlag, Hlags];

start = mconf.p + 1;
Y0 = Z(start:end,:);
X0 = Xfull(start:end,:);

good = all(isfinite(X0),2) & all(isfinite(Y0),2);
idx  = (start:T0)'; idx = idx(good);

% Determine regime for these indices
Sidx = st.S(idx);

% Omega_bar per regime: inv(A)*Sigma*inv(A)'
Sigma = diag(st.sigma2(:));  % n x n
A1 = st.A1; A2 = st.A2;

OmegaBar1 = A1 \ (Sigma / (A1'));
OmegaBar2 = A2 \ (Sigma / (A2'));

OmegaBarInv1 = inv((OmegaBar1+OmegaBar1')/2);
OmegaBarInv2 = inv((OmegaBar2+OmegaBar2')/2);

% Draw each regime separately
st.Phi1 = draw_one_regime_phi(idx, Sidx==1, Z, Xfull, st.h, OmegaBarInv1, mconf, pconf);
st.Phi2 = draw_one_regime_phi(idx, Sidx==0, Z, Xfull, st.h, OmegaBarInv2, mconf, pconf);

end

function Phi = draw_one_regime_phi(idx, mask, Z, Xfull, h, OmegaInv, mconf, pconf)
% Create scaled regression: divide y and x by sqrt(lambda_t)
% so error cov is OmegaBar constant.

n  = size(Z,2);
id = idx(mask);

k = size(Xfull,2);

if numel(id) < 20
    % too few obs -> keep a safe fallback (zeros is dangerous)
    Phi = zeros(n, k);
    return;
end

Y = Z(id,:);
X = Xfull(id,:);

lambda = exp(h(id));
w      = 1 ./ sqrt(lambda);

Yt = Y .* w;   % scale rows
Xt = X .* w;

% Prior on vec(B), B is k x n in Y = X B + E
b0 = pconf.Phi.mean * ones(k*n,1);
V0 = diag(build_minnesota_diag(n, mconf.p, mconf.J, pconf.Phi.tight, pconf.Phi.lag_decay));
V0inv = diag(1./diag(V0));

% Posterior precision: V0^{-1} + (OmegaInv ⊗ X'X)
XX   = Xt' * Xt;
Prec = V0inv + kron(OmegaInv, XX);

% Posterior covariance
Vpost = inv((Prec + Prec')/2);

% RHS: V0^{-1} b0 + vec(X'Y * OmegaInv)
XY  = Xt' * Yt;                     % k x n
rhs = V0inv*b0 + reshape(XY * OmegaInv, [], 1);

m = Vpost * rhs;

% Draw + stationarity truncation
tries = 0;
while true
    tries = tries + 1;

    % robust chol
    C = (Vpost + Vpost')/2;
    jitter = 1e-10;
    for it=1:6
        [L,p] = chol(C + jitter*eye(size(C)), 'lower');
        if p==0, break; end
        jitter = jitter * 10;
    end
    if p~=0
        L = chol(C + 1e-6*eye(size(C)), 'lower');
    end

    bdraw = m + L*randn(k*n,1);

    B   = reshape(bdraw, [k n]);  % k x n
    Phi = B';                     % n x k

    if ~mconf.enforce_stationarity
        break;
    end
    if cf.model.is_stable_var(Phi, n, mconf.p)
        break;
    end
    if tries >= mconf.stationarity_max_tries
        break;
    end
end
end

function d = build_minnesota_diag(n, p, J, tight, lag_decay)
% Diagonal prior variances for vec(B) (k*n elements), simplistic Minnesota-like.
% k = 1 + n*p + (J+1)
k = 1 + n*p + (J+1);
d = zeros(k*n,1);

for eq=1:n
    for col=1:k
        ii = (eq-1)*k + col;
        if col==1
            d(ii) = 10*tight^2; % loose on intercept
        elseif col <= 1+n*p
            lagnum = ceil((col-1)/n);
            d(ii)  = (tight^2) / (lagnum^lag_decay);
        else
            d(ii)  = (tight^2) / 4;
        end
    end
end
d = max(d, 1e-8);
end
