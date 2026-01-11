function [Z0_out, Zs_out] = simulate_pair(Zhist, hHist, D, mconf, iconf, shockType, e_all, eta_all)
% Simulate baseline and shock scenario for H horizons, allowing regimes to switch endogenously.
%
% IMPORTANT (Appendix B): at t=1 we must have
%   e_1^0 = 0 (vector)
%   e_1^δ = δ in the targeted shock component, and 0 elsewhere.
% For t>=2, use common random numbers (same e_t in baseline and shock).

H = iconf.H;
p = iconf.p;
p0 = size(Zhist,1);
J  = iconf.J;
n  = size(Zhist,2);

% allocate full histories (p0 + H)
Z0 = zeros(p0+H, n);
Zs = zeros(p0+H, n);
Z0(1:p0,:) = Zhist;
Zs(1:p0,:) = Zhist;

% lambda (log) history for regressor: h0 ordered [h_t, h_{t-1}, ..., h_{t-J}]
h0     = hHist(:);
h_path = zeros(H,1);

% forecast log-lambda (eq. (5)), common to both scenarios
h_prev = h0(1);
Qpos   = max(D.Q, 1e-12);
for t=1:H
    h_prev    = D.mu + D.F*(h_prev - D.mu) + sqrt(Qpos)*eta_all(t);
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

    % build regressor x_t = [1; Z_{t-1}; ...; Z_{t-p}; h_t; ...; h_{t-J}]
    x   = cf.irf.make_xt(Z0, idx, p, h_path, h0, t, J);
    muZ = Phi * x; % n x 1

    % impact matrix B = A^{-1} * sqrt(lambda_t) * diag(sigma)
    lam_t = exp(h_path(t));
    B     = A \ (sqrt(lam_t) * diag(D.sigma));

    % shocks:
    % t=1: e1^0 = 0 vector
    % t>=2: common random shocks
    if t==1
        e = zeros(n,1);
    else
        e = e_all(t,:)';
    end

    Z0(idx,:) = (muZ + B*e).';

    % guard against numerical explosion
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

    x   = cf.irf.make_xt(Zs, idx, p, h_path, h0, t, J);
    muZ = Phi * x;

    lam_t = exp(h_path(t));
    B     = A \ (sqrt(lam_t) * diag(D.sigma));

    % common random numbers:
    if t==1
        e = zeros(n,1);  % start from 0 vector at impact
    else
        e = e_all(t,:)';
    end

    % shock injection only at t=1
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
    % We implement "IRF(FFR)=0" by matching the BASELINE FFR path.
    if shockType=="liquidity" && t <= iconf.ffr_lock_h
        % If the impact of the control shock on FFR is too small, enforcing the lock
        % requires a gigantic control shock -> explosive scales. Skip such draws.
        jctrl = iconf.shock_ffr_control;
        denom = B(iconf.idx_ffr, jctrl);

        minImpact = 1e-6;
        if isfield(iconf,'min_ffr_control_impact') && ~isempty(iconf.min_ffr_control_impact)
            minImpact = iconf.min_ffr_control_impact;
        end

        if ~isfinite(denom) || abs(denom) < minImpact
            Z0_out = nan(H,n); Zs_out = nan(H,n);
            return
        end

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