function st = draw_phi_sur(st, Z, mconf, pconf)
% Draw Phi1 and Phi2 using multivariate regression with known Omega_bar.

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

% IMPORTANT: exclude undefined-regime periods (t<=d)
if isfield(st,'Svalid') && ~isempty(st.Svalid)
    idx = idx(st.Svalid(idx));
end

Sidx = st.S(idx);

Sigma = diag(st.sigma2(:));
A1 = st.A1; A2 = st.A2;

OmegaBar1 = A1 \ (Sigma / (A1'));
OmegaBar2 = A2 \ (Sigma / (A2'));

OmegaBarInv1 = inv((OmegaBar1+OmegaBar1')/2);
OmegaBarInv2 = inv((OmegaBar2+OmegaBar2')/2);

st.Phi1 = draw_one_regime_phi(idx, Sidx==1, Z, Xfull, st.h, OmegaBarInv1, mconf, pconf, st.Phi1);
st.Phi2 = draw_one_regime_phi(idx, Sidx==0, Z, Xfull, st.h, OmegaBarInv2, mconf, pconf, st.Phi2);

end

function Phi = draw_one_regime_phi(idx, mask, Z, Xfull, h, OmegaInv, mconf, pconf, Phi_prev)
n  = size(Z,2);
id = idx(mask);
k  = size(Xfull,2);

if numel(id) < 20
    % keep previous draw (better than zeros)
    Phi = Phi_prev;
    return;
end

Y = Z(id,:);
X = Xfull(id,:);

lambda = exp(h(id));
w      = 1 ./ sqrt(lambda);

Yt = Y .* w;
Xt = X .* w;

b0 = pconf.Phi.mean * ones(k*n,1);
V0 = diag(build_minnesota_diag(n, mconf.p, mconf.J, pconf.Phi.tight, pconf.Phi.lag_decay));
V0inv = diag(1./diag(V0));

XX   = Xt' * Xt;
Prec = V0inv + kron(OmegaInv, XX);

Vpost = inv((Prec + Prec')/2);

XY  = Xt' * Yt;
rhs = V0inv*b0 + reshape(XY * OmegaInv, [], 1);

m = Vpost * rhs;

tries = 0;
while true
    tries = tries + 1;

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

    B   = reshape(bdraw, [k n]);
    Phi = B';

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
k = 1 + n*p + (J+1);
d = zeros(k*n,1);

for eq=1:n
    for col=1:k
        ii = (eq-1)*k + col;
        if col==1
            d(ii) = 10*tight^2;
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
