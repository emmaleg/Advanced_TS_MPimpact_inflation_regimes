function [Z0_out, Zs_out] = simulate_pair(Zhist, hHist, D, mconf, iconf, shockType, e_all, eta_all)
% Simulate baseline and shock scenario for H horizons, allowing regimes to switch endogenously.

H = iconf.H;
p = iconf.p;                 % VAR lag order
p0 = size(Zhist,1);          % stored history length (>= max(p,d))
J = iconf.J;
n = size(Zhist,2);

% allocate full histories (p0 + H)
Z0 = zeros(p0+H, n);
Zs = zeros(p0+H, n);
Z0(1:p0,:) = Zhist;
Zs(1:p0,:) = Zhist;

% lambda (log) history for regressor: [h_t, h_{t-1}, h_{t-2}]
h0 = hHist(:);  % length J+1, ordered [h_t, h_{t-1}, h_{t-2}]
h_path = zeros(H,1);
h_prev = h0(1);

for t=1:H
    % forecast log-lambda (eq. (5))
    h_prev = D.mu + D.F*(h_prev - D.mu) + sqrt(D.Q)*eta_all(t);
    h_path(t) = h_prev;
end

% --- baseline first ---
for t=1:H
    idx = p0 + t;

    % determine regime using inflation at (t - d) in this scenario (eq. (2))
    S_low = cf.irf.is_low_regime(Z0, idx, D.Pstar, D.d, iconf.idx_inf);

    % select parameters by regime
    if S_low
        Phi = D.Phi1; A = D.A1;
    else
        Phi = D.Phi2; A = D.A2;
    end

    % build regressor x_t = [1; Z_{t-1}; ...; Z_{t-p}; h_t; h_{t-1}; h_{t-2}]
    x = cf.irf.make_xt(Z0, idx, p, h_path, h0, t, J);

    muZ = Phi * x; % n x 1

    % impact matrix B = A^{-1} * sqrt(lambda_t) * diag(sigma)
    lam_t = exp(h_path(t));
    B = (A \ eye(n)) * (sqrt(lam_t) * diag(D.sigma));

    e = e_all(t,:)';
    % no-shock: zero out target component at t=1 only
    if t==1
        if shockType=="conventional"
            e(iconf.shock_conventional) = 0;
        else
            e(iconf.shock_liquidity) = 0;
        end
    end

    Z0(idx,:) = (muZ + B*e).';

    % guard against numerical explosion (explosive draws can blow up late horizons)
    if any(~isfinite(Z0(idx,:))) || max(abs(Z0(idx,:))) > 1e6
        Z0_out = nan(H,n); Zs_out = nan(H,n);
        return
    end
end

% --- shock scenario ---
for t=1:H
    idx = p0 + t;

    S_low = cf.irf.is_low_regime(Zs, idx, D.Pstar, D.d, iconf.idx_inf);

    if S_low
        Phi = D.Phi1; A = D.A1;
    else
        Phi = D.Phi2; A = D.A2;
    end

    x = cf.irf.make_xt(Zs, idx, p, h_path, h0, t, J);
    muZ = Phi * x;

    lam_t = exp(h_path(t));
    B = (A \ eye(n)) * (sqrt(lam_t) * diag(D.sigma));

    e = e_all(t,:)';

    % sign normalization + delta at t=1
    if t==1
        if shockType=="conventional"
            sgn = cf.irf.normalize_sign_conventional(B, iconf);
            e(iconf.shock_conventional) = sgn * iconf.delta;
        else
            sgn = cf.irf.normalize_sign_liquidity(B, iconf);
            e(iconf.shock_liquidity) = sgn * iconf.delta;
        end
    end

    % liquidity: enforce FFR locked for first ffr_lock_h months (Appendix B step 5c)
    if shockType=="liquidity" && t <= iconf.ffr_lock_h
        r_desired = Z0(idx, iconf.idx_ffr); % match baseline path => IRF(R)=0
        e = cf.irf.enforce_ffr(B, muZ, e, r_desired, iconf);
    end

    Zs(idx,:) = (muZ + B*e).';

    if any(~isfinite(Zs(idx,:))) || max(abs(Zs(idx,:))) > 1e6
        Z0_out = nan(H,n); Zs_out = nan(H,n);
        return
    end
end

Z0_out = Z0(p0+1:p0+H, :);
Zs_out = Zs(p0+1:p0+H, :);

end
