function st = draw_phi_sur(st, Z, mconf, pconf)
% Draw Phi1 and Phi2 using multivariate regression with known Omega_bar (time-varying scalar lambda).
% Appendix A.6-A.7 structure. 

n = size(Z,2);

% Build (Y,X) for all t consistent with h lags etc.
[YY, XX] = cf.model.build_regressors(Z, st.h, mconf.p, mconf.J);

% Recover time indices used
[T0,~] = size(Z);
Xlag = cf.model.make_lag_matrix(Z,mconf.p);
Hlags = NaN(T0,mconf.J+1);
for j=0:mconf.J, Hlags(:,j+1)=lagmatrix(st.h,j); end
Xfull = [ones(T0,1), Xlag, Hlags];
start = mconf.p+1;
Y0 = Z(start:end,:);
X0 = Xfull(start:end,:);
good = all(isfinite(X0),2) & all(isfinite(Y0),2);
idx = (start:T0)'; idx = idx(good);

% Determine regime for these indices
Sidx = st.S(idx);

% Precompute Omega_bar^{-1} per regime (Omega_t = lambda_t * Omega_bar)
Sigma = diag(st.sigma2);
A1 = st.A1; A2 = st.A2;
OmegaBar1 = (A1 \ (Sigma / (A1')));  % inv(A1)*Sigma*inv(A1)'
OmegaBar2 = (A2 \ (Sigma / (A2')));

OmegaBarInv1 = inv(OmegaBar1);
OmegaBarInv2 = inv(OmegaBar2);

% Draw each regime separately
st.Phi1 = draw_one_regime_phi(idx, Sidx==1, Z, Xfull, st.h, OmegaBarInv1, mconf, pconf, 1);
st.Phi2 = draw_one_regime_phi(idx, Sidx==0, Z, Xfull, st.h, OmegaBarInv2, mconf, pconf, 2);

end

function Phi = draw_one_regime_phi(idx, mask, Z, Xfull, h, OmegaInv, mconf, pconf, rid)
% Create scaled regression: divide y and x by sqrt(lambda_t) so error cov is OmegaBar constant.
n = size(Z,2);
id = idx(mask);
if numel(id) < 20
    % too few obs: keep previous (caller already has old, but we return zeros if needed)
    Phi = zeros(n, size(Xfull,2));
    return;
end

Y = Z(id,:);
X = Xfull(id,:);
lambda = exp(h(id));
w = 1 ./ sqrt(lambda);

Yt = Y .* w;       % elementwise scale rows
Xt = X .* w;

k = size(Xt,2);

% Prior on vec(B), where B is k x n in Y = X B + E
b0 = pconf.Phi.mean * ones(k*n,1);
V0 = diag( build_minnesota_diag(n, mconf.p, mconf.J, pconf.Phi.tight, pconf.Phi.lag_decay) );

% Posterior precision: V0^{-1} + kron(X'X, OmegaInv)
XX = Xt' * Xt;
Prec = inv(V0) + kron(XX, OmegaInv);

Vpost = inv(Prec);

% Right-hand term: vec(X'Y*OmegaInv)
XY = Xt' * Yt;         % k x n
rhs = inv(V0)*b0 + reshape(XY * OmegaInv, [], 1);

m = Vpost * rhs;

% Draw + stationarity truncation
tries = 0;
while true
    tries = tries + 1;
    L = chol((Vpost+Vpost')/2, 'lower');
    bdraw = m + L*randn(k*n,1);

    B = reshape(bdraw, [k n]);   % k x n
    Phi = B';                    % n x k

    if ~mconf.enforce_stationarity
        break;
    end
    if cf.model.is_stable_var(Phi, n, mconf.p)
        break;
    end
    if tries >= mconf.stationarity_max_tries
        % fallback: accept last draw even if unstable
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
        idx = (eq-1)*k + col;
        if col==1
            d(idx) = 10*tight^2; % loose on intercept
        elseif col <= 1+n*p
            % lagged Z
            lagnum = ceil((col-1)/n);
            d(idx) = (tight^2) / (lagnum^lag_decay);
        else
            % h-lags
            d(idx) = (tight^2) / 4;
        end
    end
end
d = max(d, 1e-8);
end
